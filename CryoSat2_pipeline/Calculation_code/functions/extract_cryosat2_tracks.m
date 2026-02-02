function track_data = extract_cryosat2_tracks(config)
% EXTRACT_CRYOSAT2_TRACKS 从CryoSat-2 NetCDF文件中提取轨道数据（RGI冰川区域）
%
% 参考：revised_ICESat2_code/functions/extract_icesat2_tracks.m
%
% 功能：
%   1. 读取NetCDF格式的CryoSat-2 L2I/TEMPO数据
%   2. 进行空间过滤和区域标识
%   3. 仅保留RGI冰川区域内的数据
%
% 输入:
%   config - 配置结构体
%
% 输出:
%   track_data - 提取的轨道数据矩阵
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
% 第14 是否在HMA
% 第15 是否在rgi

warning('off','MATLAB:imagesci:geotiffread:libraryError');

fprintf('开始提取CryoSat-2轨道数据（RGI冰川区域）...\n');

%% 加载空间参考数据
fprintf('加载空间参考数据...\n');
HMA_boundary = shaperead(config.files.hma_boundary);
HMA22_regions = shaperead(config.files.hma22_regions);
GTN_regions = shaperead(config.files.gtn_regions);
TP_boundary = shaperead(config.files.tp_boundary);
grid05 = shaperead(config.files.grid05);
grid10 = shaperead(config.files.grid10);

[glacier_mask, ref_glacier] = geotiffread(config.files.glacier_mask);
[p, q] = size(glacier_mask);

fprintf('空间参考数据加载完成。\n');

%% 初始化数据存储
track_data = [];
processed_files = 0;

% %% 处理L2I数据
% if config.data_source.use_l2i
%     fprintf('\n=== 处理L2I数据 ===\n');
%     l2i_data = process_l2i_files(config, HMA_boundary, HMA22_regions, GTN_regions, ...
%         TP_boundary, grid05, grid10, glacier_mask, ref_glacier, p, q);
%     track_data = [track_data; l2i_data];
% end

%% 处理CryoTEMPO数据
if config.data_source.use_tempo
    fprintf('\n=== 处理CryoTEMPO数据 ===\n');
    tempo_data = process_tempo_files(config, HMA_boundary, HMA22_regions, GTN_regions, ...
        TP_boundary, grid05, grid10, glacier_mask, ref_glacier, p, q);
    track_data = [track_data; tempo_data];
end

fprintf('\n✓ 数据提取完成，共 %d 个数据点\n', size(track_data, 1));
track_data = track_data(:, 1:13);
end


function tempo_data = process_tempo_files(config, HMA_boundary, HMA22_regions, GTN_regions, ...
    TP_boundary, grid05, grid10, glacier_mask, ref_glacier, p, q)
% PROCESS_TEMPO_FILES 处理CryoTEMPO格式数据
path = config.paths.input_tempo_dir;
file_list = cell(length(path), 1);  % 预分配cell数组
for i = 1:length(path)
    file_list{i} = dir(fullfile(path{i}, '*.nc'));
end

tempo_data = [];

nc_files = vertcat(file_list{:});% struct形式
fprintf('找到 %d 个TEMPO NetCDF文件\n', length(nc_files));

for i = 1:length(nc_files)
    try
        fprintf('处理文件 %d/%d: %s\n', i, length(nc_files), nc_files(i).name);
        
        file_path = fullfile(nc_files(i).folder, nc_files(i).name);
        
        % 读取NetCDF数据
        lon = ncread(file_path, 'lon');
        lat = ncread(file_path, 'lat');
        elevation = ncread(file_path, 'elevation');
        time_total = ncread(file_path, 'time');
        time_total = floor(time_total / 86400);% seconds->days, since 2000-01-01-00:00:00
        uncertainty = ncread(file_path, 'uncertainty');
        
        % 从文件名提取年份和月份
        yr = nc_files(i).name(30:33);
        mon = nc_files(i).name(35:36);
        
        % 计算大地水准面高度
        geoid = geoidheight(lat, lon);%EGM96
        
        % 组织数据
        point = zeros(length(lon), 15);
        point(:,1) = lon;
        point(:,2) = lat;
        point(:,3) = elevation;
        point(:,4) = time_total; %days, since 2000-01-01-00:00:00
        point(:,5) = ones(length(lat),1)*str2num(yr);
        point(:,6) = ones(length(lat),1)*str2num(mon);
        point(:,7) = geoid;
        point(:,8) = uncertainty;
        
        % 空间过滤
        [in, ~] = inpolygon(lon, lat, HMA_boundary.X, HMA_boundary.Y);
        point(:,14) = in;
        
        toDelete = (point(:,14) == 0);
        point(toDelete, :) = [];
        
        if isempty(point)
            continue;
        end
        
        % 应用冰川掩膜
        point = apply_glacier_mask_rgi(point, glacier_mask, ref_glacier, p, q);
        
        if isempty(point)
            continue;
        end
        
        % 分配区域编号
        point = assign_spatial_regions(point, HMA22_regions, GTN_regions, ...
            TP_boundary, grid05, grid10);
        
        tempo_data = [tempo_data; point];
        
    catch ME
        fprintf('文件处理失败: %s\n', ME.message);
        continue;
    end
end

end

function point = apply_glacier_mask_rgi(point, glacier_mask, ref_glacier, p, q)
% APPLY_GLACIER_MASK_RGI 应用冰川掩膜，仅保留RGI冰川区域内的点

point(:,15) = 1;  % 初始化冰川标识，1代表在冰川内

for ii = 1:size(point, 1)
    [row, col] = latlon2pix(ref_glacier, point(ii, 2), point(ii, 1));
    row = ceil(row);
    col = ceil(col);
    
    if row > p || row < 1 || col > q || col < 1
        point(ii, 15) = 0;  % 标记为非冰川
        continue;
    end
    
    % 注意：不同数据集的非冰川标识值可能不同
    if glacier_mask(row, col) == -2147483647
        point(ii, 15) = 0;  % 标记为非冰川
    end

end

% 排除非冰川区域
to_delete = (point(:,15) == 0);
point(to_delete, :) = [];

end

function point = assign_spatial_regions(point, HMA22_regions, GTN_regions, ...
    TP_boundary, grid05, grid10)
% ASSIGN_SPATIAL_REGIONS 分配空间区域编号

if isempty(point)
    return;
end

% HMA22区域编号
for k = 1:length(HMA22_regions)
    [in, ~] = inpolygon(point(:,1), point(:,2), HMA22_regions(k).X, HMA22_regions(k).Y);
    point(in==1, 9) = k;
end

% GTN区域编号
for k = 1:length(GTN_regions)
    [in, ~] = inpolygon(point(:,1), point(:,2), GTN_regions(k).X, GTN_regions(k).Y);
    point(in==1, 10) = k;
end

% 0.5°网格编号
for k = 1:length(grid05)
    [in, ~] = inpolygon(point(:,1), point(:,2), grid05(k).X, grid05(k).Y);
    point(in==1, 11) = k;
end

% 1°网格编号
for k = 1:length(grid10)
    [in, ~] = inpolygon(point(:,1), point(:,2), grid10(k).X, grid10(k).Y);
    point(in==1, 12) = k;
end

% 青藏高原标识
[in, ~] = inpolygon(point(:,1), point(:,2), TP_boundary.X, TP_boundary.Y);
point(:,13) = in;

end

