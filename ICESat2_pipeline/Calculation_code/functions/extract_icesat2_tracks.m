function track_data = extract_icesat2_tracks(config)
% EXTRACT_ICESAT2_TRACKS 从ICESat-2 ATL06 HDF5文件中提取轨道数据
%
% 功能：
%   1. 读取HDF5格式的ICESat-2 ATL06 v6数据
%   2. 进行空间过滤和区域标识
%   3. 分配网格编号和区域编号
%
% 输入:
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
%               列9: RGI区域编号  
%               列10: 0.5°网格编号
%               列11: 1°网格编号
%               列12: 青藏高原标识
%               列13: 冰川标识

% 关闭HDF5警告
warning('off','MATLAB:imagesci:hdf5dataset:libraryError');
warning('off','MATLAB:imagesci:h5info:libraryError');
warning('off','MATLAB:imagesci:geotiffread:libraryError');

fprintf('开始提取ICESat-2轨道数据...\n');

%% 加载空间参考数据
fprintf('加载空间参考数据...\n');

% 加载边界和网格数据
HMA_boundary = shaperead(config.files.hma_boundary);
HMA22_regions = shaperead(config.files.hma22_regions);
GTN_regions = shaperead(config.files.gtn_regions);
TP_boundary = shaperead(config.files.tp_boundary);
grid05 = shaperead(config.files.grid05);
grid10 = shaperead(config.files.grid10);

% 加载冰川掩膜
[glacier_mask, ref_glacier] = geotiffread(config.files.glacier_mask);
info_glacier = geotiffinfo(config.files.glacier_mask);
[p, q] = size(glacier_mask);

fprintf('空间参考数据加载完成。\n');

%% 获取HDF5文件列表
h5_files = dir(fullfile(config.paths.input_h5_dir, '*.h5'));
fprintf('找到 %d 个HDF5文件需要处理。\n', length(h5_files));

%% 定义激光轨道
ground_tracks = ['gt1l'; 'gt1r'; 'gt2l'; 'gt2r'; 'gt3l'; 'gt3r'];

%% 初始化数据存储
track_data = [];
processed_files = 0;

%% 逐文件处理
for i = 1:length(h5_files)
    try
        fprintf('处理文件 %d/%d: %s\n', i, length(h5_files), h5_files(i).name);
        
        % 处理单个HDF5文件
        file_data = process_single_h5_file(h5_files(i), config.paths.input_h5_dir, ...
            ground_tracks, HMA_boundary, HMA22_regions, GTN_regions, TP_boundary, ...
            grid05, grid10, glacier_mask, ref_glacier, config);
        
        if ~isempty(file_data)
            track_data = [track_data; file_data];
            processed_files = processed_files + 1;
        end
        
    catch ME
        fprintf('❌ 文件处理失败 %s: %s\n', h5_files(i).name, ME.message);
        continue;
    end
end

fprintf('✓ 成功处理 %d/%d 个文件，提取 %d 个数据点。\n', ...
    processed_files, length(h5_files), size(track_data, 1));

end

function file_data = process_single_h5_file(h5_file, input_dir, ground_tracks, ...
    HMA_boundary, HMA22_regions, GTN_regions, TP_boundary, grid05, grid10, ...
    glacier_mask, ref_glacier, config)
% PROCESS_SINGLE_H5_FILE 处理单个HDF5文件
%
% 输入:
%   h5_file - 文件信息结构
%   input_dir - 输入目录
%   ground_tracks - 激光轨道名称
%   其他 - 空间参考数据
%   config - 配置参数，读取参数
%
% 输出:
%   file_data - 文件数据矩阵

file_path = fullfile(input_dir, h5_file.name);
file_data = [];

try
    % 获取HDF5文件信息
    info = h5info(file_path);
    num_groups = length(info.Groups);
    
    % 处理6个激光轨道（跳过前后两个元数据组）
    for j = 3:(num_groups-2)
        try
            track_data = extract_single_track(file_path, info.Groups(j).Name, config);
            
            if isempty(track_data)
                continue;
            end
            
            % 空间过滤和区域标识
            track_data = apply_spatial_filters(track_data, HMA_boundary, config);
            
            if isempty(track_data)
                continue;
            end
            
            % 冰川掩膜检查
            track_data = apply_glacier_mask(track_data, glacier_mask, ref_glacier);
            
            if isempty(track_data)
                continue;
            end
            
            % 区域标识
            track_data = assign_spatial_regions(track_data, HMA22_regions, GTN_regions, ...
                TP_boundary, grid05, grid10);
            
            % 合并数据
            file_data = [file_data; track_data];
            
        catch ME
            fprintf('   轨道 %s 处理失败: %s\n', info.Groups(j).Name, ME.message);
            continue;
        end
    end
    
catch ME
    fprintf('   HDF5文件读取失败: %s\n', ME.message);
end

end

function track_data = extract_single_track(file_path, group_name, config)
% EXTRACT_SINGLE_TRACK 提取单个激光轨道数据
%
% 输入:
%   file_path - HDF5文件路径
%   group_name - 轨道组名
%   config - 配置参数
%
% 输出:
%   track_data - 轨道数据矩阵

try
    % 读取基本测量数据
    flag = h5read(file_path, [group_name '/land_ice_segments/atl06_quality_summary']);
    lon = h5read(file_path, [group_name '/land_ice_segments/longitude']);
    lat = h5read(file_path, [group_name '/land_ice_segments/latitude']);
    time = h5read(file_path, [group_name '/land_ice_segments/delta_time']);
    h_li = h5read(file_path, [group_name '/land_ice_segments/h_li']);
    geoid = h5read(file_path, [group_name '/land_ice_segments/dem/geoid_h']);
    
    % 构建数据矩阵
    track_data = zeros(length(lon), 7);
    track_data(:,1) = lon;
    track_data(:,2) = lat;
    track_data(:,3) = floor(time/86400);  % 转换为天数
    track_data(:,4) = h_li;
    track_data(:,5) = geoidheight(lat, lon, 'EGM96');  % EGM96大地水准面高度
    track_data(:,6) = geoid;  % EGM2008大地水准面高度
    track_data(:,7) = flag;
    
catch ME
    fprintf('      数据读取失败: %s\n', ME.message);
    track_data = [];
end

end

function track_data = apply_spatial_filters(track_data, HMA_boundary, config)
% APPLY_SPATIAL_FILTERS 应用空间过滤器
%
% 输入:
%   track_data - 轨道数据
%   HMA_boundary - HMA边界
%   config - 配置参数
%
% 输出:
%   track_data - 过滤后的数据

if isempty(track_data)
    return;
end

% HMA边界内的点
[in, ~] = inpolygon(track_data(:,1), track_data(:,2), HMA_boundary.X, HMA_boundary.Y);
track_data(:,8) = in;

% 排除HMA边界外的点
to_delete = (track_data(:,8) == 0);
track_data(to_delete, :) = [];

% 排除高程异常值
to_delete = (track_data(:,4) > config.quality.max_elevation);
track_data(to_delete, :) = [];

end

function track_data = apply_glacier_mask(track_data, glacier_mask, ref_glacier)
% APPLY_GLACIER_MASK 应用冰川掩膜
%
% 输入:
%   track_data - 轨道数据
%   glacier_mask - 冰川掩膜栅格
%   ref_glacier - 栅格参考信息
%
% 输出:
%   track_data - 标识冰川的数据

if isempty(track_data)
    return;
end

[p, q] = size(glacier_mask);
track_data(:,13) = 0;  % 初始化冰川标识

for ii = 1:size(track_data, 1)
    try
        [row, col] = latlon2pix(ref_glacier, track_data(ii, 2), track_data(ii, 1));
        row = ceil(row);
        col = ceil(col);
        
        if row > p || row < 1 || col > q || col < 1
            track_data(ii, 13) = 1;
            continue;
        end
        
        if glacier_mask(row, col) == -2147483647  % 非冰川区域标识
            track_data(ii, 13) = 1;  % 标记为非冰川
        end
        
    catch
        continue;
    end
end

to_delete = (track_data(:,13) == 1);
track_data(to_delete, :) = [];

end

function track_data = assign_spatial_regions(track_data, HMA22_regions, GTN_regions, ...
    TP_boundary, grid05, grid10)
% ASSIGN_SPATIAL_REGIONS 分配空间区域编号
%
% 输入:
%   track_data - 轨道数据
%   各种区域边界数据
%
% 输出:
%   track_data - 带有区域编号的数据

if isempty(track_data)
    return;
end

% 初始化区域编号列
track_data(:,[8:12]) = 0;

% HMA22区域编号
for k = 1:length(HMA22_regions)
    [in, ~] = inpolygon(track_data(:,1), track_data(:,2), ...
        HMA22_regions(k).X, HMA22_regions(k).Y);
    track_data(in==1, 8) = k;
end

% GTN区域编号
for k = 1:length(GTN_regions)
    [in, ~] = inpolygon(track_data(:,1), track_data(:,2), ...
        GTN_regions(k).X, GTN_regions(k).Y);
    track_data(in==1, 9) = k;
end

% 0.5°网格编号
for k = 1:length(grid05)
    [in, ~] = inpolygon(track_data(:,1), track_data(:,2), ...
        grid05(k).X, grid05(k).Y);
    track_data(in==1, 10) = k;
end

% 1°网格编号
for k = 1:length(grid10)
    [in, ~] = inpolygon(track_data(:,1), track_data(:,2), ...
        grid10(k).X, grid10(k).Y);
    track_data(in==1, 11) = k;
end

% 青藏高原标识
[in, ~] = inpolygon(track_data(:,1), track_data(:,2), ...
    TP_boundary.X, TP_boundary.Y);
track_data(:,12) = in;

end
