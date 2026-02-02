function enhanced_data = enhance_icesat_data(txt_file, config)
% ENHANCE_ICESAT_DATA 增强ICESat数据，添加地形信息
%
% 关键差异：ICESat使用3×3邻域平均NASADEM（vs CryoSat-2的11×11）
% 参考：HMA_footprint_track_20230213.m 第3-28行
%
% 功能：
%   1. 提取地形参数（3×3邻域）：NASADEM、slope、aspect
%   2. 计算h_li2000（WGS84高程）
%   3. 添加pdd校正
%   4. 计算高程变化（dh和dh_with_pdd）
%   5. 添加HMA4、HMA6、TP区域标识
%
% ICESat数据列定义（输入14列，注意lon在前lat在后）:
% 列1: lon (经度)
% 列2: lat (纬度)
% 列3: elevation (观测高程)
% 列4: beam_azimuth (波束方位角)
% 列5: days (天数，since 2000-01-01)
% 列6: geoid (EGM96大地水准面)
% 列7: satcorr (卫星高程校正)
% 列8: satflag (卫星校正标识)
% 列9: glacier_flag (冰川标识, 0=冰川内)
% 列10: year (年份)
% 列11: hma_region (HMA区域编号)
% 列12: hma22_region (HMA22区域编号)
% 列13: gtn_region (GTN15区域编号)
% 列14: grid_region (grid496区域编号)
%
% 输出（扩展到24列）:
% 列15: nasadem (3×3邻域平均) ⭐
% 列16: slope (3×3邻域平均)
% 列17: aspect (3×3邻域平均)
% 列18: pdd_corr (PDD校正值)
% 列19: h_li2000 (WGS84高程，校正后)
% 列20: dh = h_li2000 - nasadem ⭐用于统计
% 列21: dh_with_pdd = h_li2000 - nasadem - pdd
% 列22: HMA4区域ID
% 列23: HMA6区域ID（季节模式）
% 列24: TP标识

fprintf('开始ICESat数据增强处理...\n');
warning('off', 'MATLAB:polyshape:repairedBySimplify');
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
    [nasadem, ref_nasadem] = geotiffread(config.files.nasadem);
    [slope_data, ref_slope] = geotiffread(config.files.slope);
    [aspect_data, ref_aspect] = geotiffread(config.files.aspect);
    
    % 加载区域数据
    glacier_grid = shaperead(config.files.grid10);
    HMA4 = shaperead(config.files.hma4_regions);
    HMA6 = shaperead(config.files.hma6_regions);
    TP = shaperead(config.files.tp_boundary);
    
    fprintf('✓ 地形数据加载完成\n');
catch ME
    error('地形数据加载失败: %s', ME.message);
end

%% 添加地形参数（3×3邻域平均）
fprintf('添加地形参数（3×3邻域）：NASADEM、slope、aspect...\n');
enhanced_data = add_terrain_parameters_3x3(data, nasadem, ref_nasadem, slope_data, ref_slope, aspect_data, ref_aspect, config);

%% 计算h_li2000和高程变化
fprintf('计算h_li2000和高程变化...\n');
enhanced_data = calculate_h_li2000_and_dh(enhanced_data);

%% 添加pdd校正并计算dh_with_pdd
fprintf('添加pdd校正并计算dh_with_pdd...\n');
enhanced_data = add_pdd_and_calculate_dh_with_pdd(enhanced_data, glacier_grid);

%% 添加区域标识
fprintf('添加区域标识（HMA4, HMA6, TP）...\n');
enhanced_data = add_region_identifiers(enhanced_data, HMA4, HMA6, TP);

fprintf('✓ 数据增强完成，输出 %d 个数据点\n', size(enhanced_data, 1));

end

function enhanced_data = add_terrain_parameters_3x3(data, nasadem, ref_nasadem, slope_data, ref_slope, aspect_data, ref_aspect, config)
% ADD_TERRAIN_PARAMETERS_3X3 添加地形参数（3×3邻域平均）
%
% 这是ICESat与CryoSat-2的关键差异！
% ICESat:    3×3邻域（9个像元）
% CryoSat-2: 11×11邻域（121个像元）
%
% 原因：ICESat足迹尺寸~70m，小于CryoSat-2的~300m
%
% 输出：
% 列15: nasadem (3×3邻域)
% 列16: slope (3×3邻域)
% 列17: aspect (3×3邻域)

neighbor_size = config.spatial.nasadem_neighbor_size;  % 3
half_size = floor(neighbor_size / 2);  % 1

% 扩展数据矩阵
enhanced_data = data;
if size(enhanced_data, 2) < 17
    enhanced_data(:, 15) = nan;  % nasadem
    enhanced_data(:, 16) = nan;  % slope
    enhanced_data(:, 17) = nan;  % aspect
end

fprintf('  使用%d×%d邻域平均（%d个像元）...\n', ...
    neighbor_size, neighbor_size, neighbor_size^2);

for i = 1:size(data, 1)
    % 获取像元坐标（注意：lon在列1，lat在列2）
    [row, col] = latlon2pix(ref_nasadem, data(i, 2), data(i, 1));
    row = ceil(row);
    col = ceil(col);
    
    % 检查边界
    if row <= half_size || row > size(nasadem, 1) - half_size || ...
       col <= half_size || col > size(nasadem, 2) - half_size
        continue;
    end
    
    % 3×3邻域平均 - NASADEM
    ele_sum = [];
    slope_sum = [];
    aspect_sum = [];
    
    for j = -half_size:half_size  % -1:1
        for k = -half_size:half_size
            ele_sum = [ele_sum, single(nasadem(row+j, col+k))];
            slope_sum = [slope_sum, single(slope_data(row+j, col+k))];
            aspect_sum = [aspect_sum, single(aspect_data(row+j, col+k))];
        end
    end
    
    enhanced_data(i, 15) = nanmean(ele_sum);      % nasadem
    enhanced_data(i, 16) = nanmean(slope_sum);    % slope
    enhanced_data(i, 17) = nanmean(aspect_sum);   % aspect
    
    % 显示进度
    if mod(i, 100000) == 0
        fprintf('  处理进度: %d/%d (%.1f%%)\n', i, size(data, 1), 100*i/size(data,1));
    end
end

fprintf('✓ 地形参数（3×3邻域）添加完成\n');

end

function enhanced_data = calculate_h_li2000_and_dh(data)
% CALCULATE_H_LI2000_AND_DH 计算h_li2000和高程变化
%
% 输出：
% 列19: h_li2000 (WGS84高程)
% 列20: dh = h_li2000 - nasadem

enhanced_data = data;
if size(enhanced_data, 2) < 20
    enhanced_data(:, 19) = nan;  % h_li2000
    enhanced_data(:, 20) = nan;  % dh
end

for i = 1:size(data, 1)
    % 计算h_li2000（参考HMA_Icesattrack_full.m）
    latitude = enhanced_data(i, 2);
    offset = 0.8337-0.3848*(sind(latitude))^2;
    if data(i, 8) == 2  % satflag
        enhanced_data(i, 19) = data(i, 3) - data(i, 6) + data(i, 7) - offset;
    else
        enhanced_data(i, 19) = data(i, 3) - data(i, 6) - offset;
    end
    
    % 计算dh = h_li2000 - nasadem,不考虑pdd
    enhanced_data(i, 20) = enhanced_data(i, 19) - data(i, 15);
end

fprintf('✓ h_li2000和dh计算完成\n');

end

function enhanced_data = add_pdd_and_calculate_dh_with_pdd(data, glacier_grid)
% ADD_PDD_AND_CALCULATE_DH_WITH_PDD 添加pdd校正并计算dh_with_pdd
%
% 输出：
% 列18: pdd_corr (从grid_region获取)
% 列21: dh_with_pdd = h_li2000 - nasadem - pdd

enhanced_data = data;
if size(enhanced_data, 2) < 21
    enhanced_data(:, 18) = nan;  % pdd_corr
    enhanced_data(:, 21) = nan;  % dh_with_pdd
end

% 根据列14 (grid_region)获取pdd校正值
for i = 1:size(data, 1)
    grid_id = data(i, 14);
    if ~isnan(grid_id) && grid_id > 0 && grid_id <= length(glacier_grid)
        enhanced_data(i, 18) = glacier_grid(grid_id).pdd;
    end
end

% 计算dh_with_pdd = h_li2000 - nasadem - pdd
enhanced_data(:, 21) = data(:, 19) - data(:, 15) - enhanced_data(:, 18);

fprintf('✓ pdd校正和dh_with_pdd计算完成\n');

end

function enhanced_data = add_region_identifiers(data, HMA4, HMA6, TP)
% ADD_REGION_IDENTIFIERS 添加区域标识
% 列22: HMA4区域ID
% 列23: HMA6区域ID（季节模式）
% 列24: TP标识

enhanced_data = data;
if size(enhanced_data, 2) < 24
    enhanced_data(:, 22) = 0;
    enhanced_data(:, 23) = 0;
    enhanced_data(:, 24) = 0;
end

% 添加HMA4区域标识（注意：lon在列1，lat在列2）
for k = 1:length(HMA4)
    [in, ~] = inpolygon(data(:,1), data(:,2), HMA4(k).X, HMA4(k).Y);
    enhanced_data(in==1, 22) = k;
end

% 添加HMA6季节模式区域标识
for k = 1:length(HMA6)
    [in, ~] = inpolygon(data(:,1), data(:,2), HMA6(k).X, HMA6(k).Y);
    enhanced_data(in==1, 23) = k;
end

% 添加TP标识
[in, ~] = inpolygon(data(:,1), data(:,2), TP.X, TP.Y);
enhanced_data(in==1, 24) = 1;

fprintf('✓ 区域标识添加完成\n');

end

