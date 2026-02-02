function enhanced_data = enhance_cryosat2_data(txt_file, config)
% ENHANCE_CRYOSAT2_DATA 增强CryoSat-2数据，添加地形信息
%
% 关键差异：CryoSat-2使用11×11邻域平均NASADEM（vs ICESat-2的3×3）
% 参考：HMA_footprint_track_20230213.m 第32-56行
%
% 功能：
%   1. 添加NASADEM高程信息（11×11邻域平均）
%   2. 计算坡度信息
%   3. 添加区域标识（HMA22, RGI, 网格）
%   4. 计算最终高程变化
%
% 第1列 lon
% 第2列 lat
% 第3 ele
% 第4 days, since 2000-01-01-00:00:00
% 第5 year
% 第6 month
% 第7 geoid %EGM96
% 第8 uncert
% 第9 HMA22区域编号
% 第10 GTN编号  
% 第11 0.5°网格编号
% 第12 1°网格编号
% 第13 青藏高原标识
% 输入:
%   txt_file - TXT数据文件路径
%   config - 配置结构体
%
% 输出:
%   enhanced_data - 增强后的数据矩阵
%                  列14: NASADEM参考高程（11×11邻域平均）
%                  列15: NASADEM参考slope（11×11邻域平均）
%                  列16: 高程变化

fprintf('开始CryoSat-2数据增强处理...\n');
warning('off', 'map:removing:latlon2pix');
%% 加载数据
fprintf('加载数据文件...\n');
try
    data = load(txt_file);
    fprintf('✓ 成功加载 %d 个数据点\n', size(data, 1));
catch ME
    error('数据文件加载失败: %s', ME.message);
end

%% 加载地形数据
fprintf('加载地形数据...\n');
try
    [srtm, ref_srtm] = geotiffread(config.files.nasadem);
    [slope_data, ref_slope] = geotiffread(config.files.slope);
    
    % 加载区域数据
    %HMA22 = shaperead(config.files.hma22_regions);
    %rgiregion = shaperead(config.files.gtn_regions);
    glacier_grid = shaperead(config.files.grid10);
    
    fprintf('✓ 地形数据加载完成\n');
catch ME
    error('地形数据加载失败: %s', ME.message);
end

%% 添加NASADEM高程信息（11×11邻域平均）
fprintf('添加NASADEM高程信息（11×11邻域平均）...\n');
enhanced_data = add_nasadem_elevation_11x11(data, srtm, ref_srtm, config);

%% 添加坡度信息
fprintf('添加坡度信息...\n');
enhanced_data = add_slope_information(enhanced_data, slope_data, ref_slope);

%% 计算高程变化
fprintf('计算高程变化...\n');
enhanced_data = calculate_height_change(enhanced_data, glacier_grid);

fprintf('✓ 数据增强完成，输出 %d 个数据点\n', size(enhanced_data, 1));

end

function enhanced_data = add_nasadem_elevation_11x11(data, srtm, ref_srtm, config)
% ADD_NASADEM_ELEVATION_11X11 添加NASADEM高程（11×11邻域平均）
%
% 这是CryoSat-2与ICESat-2的关键差异！
% CryoSat-2: 11×11邻域（121个像元）
% ICESat-2:  3×3邻域（9个像元）
%
% 原因：CryoSat-2足迹尺寸~300m，远大于ICESat-2的~17m

neighbor_size = config.spatial.nasadem_neighbor_size;  % 11
half_size = floor(neighbor_size / 2);  % 5

% 扩展数据矩阵
enhanced_data = data;
if size(enhanced_data, 2) < 14
    enhanced_data(:, 14) = nan;
end

fprintf('  使用%d×%d邻域平均（%d个像元）...\n', neighbor_size, neighbor_size, neighbor_size^2);
% 0代表处理异常
for i = 1:size(data, 1)
    % 获取像元坐标
    [row, col] = latlon2pix(ref_srtm, data(i, 2), data(i, 1));
    row = ceil(row);
    col = ceil(col);
    
    % 检查边界
    if row <= half_size || row > size(srtm, 1) - half_size || ...
       col <= half_size || col > size(srtm, 2) - half_size
        continue;
    end
    
    % 11×11邻域平均,可能会有缺值
    ele_list = [];
    for j = -half_size:half_size  % -5:5
        for k = -half_size:half_size
            ele_list = [ele_list, (srtm(row+j, col+k))];
        end
    end
    
    enhanced_data(i, 14) = nanmean(ele_list);  % ele/121
    
    % 显示进度
    if mod(i, 10000) == 0
        fprintf('  处理进度: %d/%d (%.1f%%)\n', i, size(data, 1), 100*i/size(data,1));
    end
end

fprintf('✓ NASADEM高程（11×11邻域）添加完成\n');

end

function enhanced_data = add_slope_information(data, slope_data, ref_slope)
% ADD_SLOPE_INFORMATION 添加坡度信息
%列15

enhanced_data = data;
if size(enhanced_data, 2) < 15
    enhanced_data(:, 15) = nan;
end

half_size = 3;
for i = 1:size(data, 1)
    [row, col] = latlon2pix(ref_slope, data(i, 2), data(i, 1));
    row = ceil(row);
    col = ceil(col);

    % 检查边界
    if row <= half_size || row > size(slope_data, 1) - half_size || ...
       col <= half_size || col > size(slope_data, 2) - half_size
        continue;
    end

    slope_list = [];
    for j = -half_size:half_size
        for k = -half_size:half_size
            slope_list = [slope_list, slope_data(row+j, col+k)];
        end
    end
    enhanced_data(i, 15) = nanmean(slope_list);
end

fprintf('✓ 坡度信息添加完成\n');

end

function enhanced_data = calculate_height_change(data, glacier_grid)
% CALCULATE_HEIGHT_CHANGE 计算高程变化
%列16
% 公式：地表高程 - NASADEM高程（11×11邻域） - 大地水准面

enhanced_data = data;

% 计算高程变化
for i=1:length(enhanced_data)
    lon_indicator = floor(enhanced_data(i, 1));
    lat_indicator = floor(enhanced_data(i, 2));
    idx = ([glacier_grid.x] >= lon_indicator) & ([glacier_grid.x] < lon_indicator+1) & ([glacier_grid.y] >= lat_indicator) & ([glacier_grid.y] < lat_indicator+1);
    pdd = glacier_grid(idx).pdd;

    enhanced_data(i, 16) = data(i, 3) - data(i, 7) - (data(i, 14) + pdd);%H-NASADEM
end

fprintf('✓ 高程变化计算完成\n');

end

