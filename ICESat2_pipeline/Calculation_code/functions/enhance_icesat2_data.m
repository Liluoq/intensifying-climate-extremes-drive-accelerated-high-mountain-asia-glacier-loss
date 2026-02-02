function enhanced_data = enhance_icesat2_data(csv_file, config)
% ENHANCE_ICESAT2_DATA 增强ICESat-2数据，添加地形信息和穿透深度校正
%
% 功能：
%   1. 添加NASADEM高程信息
%   2. 计算坡度和坡向
%   3. 应用穿透深度校正
%   4. 计算最终高程变化
%
% 输入:
%   csv_file - CSV数据文件路径
%   config - 配置结构体
%
% 输出:
%   track_data - 提取的轨道数据矩阵
%               列1-2: 经度、纬度
%               列3: 时间（天数，相对于参考日期）
%               列4: 观测高程 (h_li)
%               列5: EGM96大地水准面高度
%               列6: EGM2008大地水准面高度
%               列7: 质量标识
%               列8: HMA22区域编号
%               列9: GTN编号  
%               列10: 0.5°网格编号
%               列11: 1°网格编号
%               列12: 青藏高原标识
%               列13: 冰川标识
%   enhanced_data - 增强后的数据矩阵
%                  在原有13列基础上新增：
%                  列14: NASADEM参考高程
%                  列15: 坡度
%                  列16: 坡向  
%                  列17: 观测高程与参考DEM差值
%                  列18: 穿透深度校正值
%                  列19: 校正后的最终高程变化
%                  列20: coregisteration v shift
% 关闭LATLON2PIX的警告
warning('off', 'map:removing:latlon2pix');

fprintf('开始数据增强处理...\n');

%% 加载数据
fprintf('加载CSV数据文件...\n');
try
    data = readmatrix(csv_file);
    fprintf('✓ 成功加载 %d 个数据点\n', size(data, 1));
catch ME
    error('CSV文件加载失败: %s', ME.message);
end

%% 加载地形数据
fprintf('加载地形数据...\n');
try
    [srtm, ref_srtm] = geotiffread(config.files.nasadem);
    [slope_data, ref_slope] = geotiffread(config.files.slope);
    [aspect_data, ref_aspect] = geotiffread(config.files.aspect);
    
    % 加载网格数据（用于穿透深度校正）
    grid10 = shaperead(config.files.grid10);
    
    fprintf('✓ 地形数据加载完成\n');
catch ME
    error('地形数据加载失败: %s', ME.message);
end

%% 添加地形信息
fprintf('添加地形信息...\n');
enhanced_data = add_terrain_information(data, srtm, ref_srtm, slope_data, ...
    ref_slope, aspect_data, ref_aspect);

%% 应用穿透深度校正
fprintf('应用穿透深度校正...\n');
enhanced_data = apply_penetration_depth_correction(enhanced_data, grid10);

%% 计算最终高程变化
fprintf('计算最终高程变化...\n');
enhanced_data = calculate_final_height_change(enhanced_data);

fprintf('✓ 数据增强完成，输出 %d 个数据点\n', size(enhanced_data, 1));

end

function enhanced_data = add_terrain_information(data, srtm, ref_srtm, slope_data, ...
    ref_slope, aspect_data, ref_aspect)
% ADD_TERRAIN_INFORMATION 添加地形信息
%
% 输入:
%   data - 原始数据
%   srtm, ref_srtm - NASADEM数据和参考信息
%   slope_data, ref_slope - 坡度数据和参考信息
%   aspect_data, ref_aspect - 坡向数据和参考信息
%
% 输出:
%   enhanced_data - 添加地形信息后的数据

% 初始化输出数据矩阵
enhanced_data = zeros(size(data, 1), size(data, 2) + 3);
enhanced_data(:, 1:size(data, 2)) = data;

% 添加地形信息
for i = 1:size(data, 1)
    % 获取NASADEM高程
    [row, col] = latlon2pix(ref_srtm, data(i, 2), data(i, 1));
    row = ceil(row);
    col = ceil(col);
    
    if row > 0 && row <= size(srtm, 1) && col > 0 && col <= size(srtm, 2)
        enhanced_data(i, 14) = double(srtm(row, col));
    end
    
    % 获取坡度信息
    [row, col] = latlon2pix(ref_slope, data(i, 2), data(i, 1));
    row = ceil(row);
    col = ceil(col);
    
    if row > 0 && row <= size(slope_data, 1) && col > 0 && col <= size(slope_data, 2)
        enhanced_data(i, 15) = double(slope_data(row, col));
    end
    
    % 获取坡向信息
    [row, col] = latlon2pix(ref_aspect, data(i, 2), data(i, 1));
    row = ceil(row);
    col = ceil(col);
    
    if row > 0 && row <= size(aspect_data, 1) && col > 0 && col <= size(aspect_data, 2)
        enhanced_data(i, 16) = double(aspect_data(row, col));
    end

    if mod(i, 100000) == 0
        fprintf('处理了 %d/%d (%.1f%%)\n', i, size(enhanced_data, 1), i/size(enhanced_data, 1)*100);
    end
end

% 计算观测高程与参考DEM的差值
enhanced_data(:, 17) = enhanced_data(:, 4) - enhanced_data(:, 5) - enhanced_data(:, 14);

end

function enhanced_data = apply_penetration_depth_correction(enhanced_data, grid10)
% APPLY_PENETRATION_DEPTH_CORRECTION 应用穿透深度校正
%
% 输入:
%   enhanced_data - 增强数据
%   grid10 - 1°网格数据（包含穿透深度校正值）
%
% 输出:
%   enhanced_data - 应用校正后的数据

% 初始化穿透深度校正列
enhanced_data(:, 18) = 0;
% 初始化水平配准列
enhanced_data(:, 20) = 0;

% 为每个网格分配穿透深度校正值
for i = 1:length(grid10)
    grid_indices = (enhanced_data(:, 11) == i);
    if any(grid_indices)
        enhanced_data(grid_indices, 18) = grid10(i).pdd;
        enhanced_data(grid_indices, 20) = grid10(i).coreg_vshi;
    end
end

end

function enhanced_data = calculate_final_height_change(enhanced_data)
% CALCULATE_FINAL_HEIGHT_CHANGE 计算最终的高程变化
%
% 输入:
%   enhanced_data - 增强数据
%
% 输出:
%   enhanced_data - 包含最终高程变化的数据

% 计算校正后的最终高程变化
% 公式：观测高程 - EGM96大地水准面 - NASADEM高程 - 穿透深度校正
enhanced_data(:, 19) = enhanced_data(:, 4) - enhanced_data(:, 5) - ...
                      (enhanced_data(:, 14) + enhanced_data(:, 18) + enhanced_data(:, 20));

% 应用质量控制过滤
to_delete = abs(enhanced_data(:, 19)) > 200;
enhanced_data(to_delete, :) = [];

fprintf('✓ 质量控制完成，排除 %d 个异常值\n', sum(to_delete));

end
