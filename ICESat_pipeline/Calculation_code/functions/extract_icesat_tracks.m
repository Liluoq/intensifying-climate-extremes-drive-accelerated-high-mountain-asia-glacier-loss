function [track_data, track_data_outrgi] = extract_icesat_tracks(config)
% EXTRACT_ICESAT_TRACKS 从ICESat GLAH14 HDF文件中提取轨道数据
%
% 功能：从ICESat GLAH14产品（HDF格式）中提取RGI冰川区域内的轨道数据
% 参考：HMA_Icesattrack_full.m
%
% 输入:
%   config - 配置结构体
%
% 输出:
%   track_data - 提取的轨道数据矩阵（14列）
%
% 输出数据列定义（注意：第1列lon，第2列lat）:
%   列1: lon (经度)
%   列2: lat (纬度)
%   列3: elevation (观测高程)
%   列4: beam_azimuth (波束方位角)
%   列5: days (天数，since 2000-01-01)
%   列6: geoid (EGM96大地水准面)
%   列7: satcorr (卫星高程校正)
%   列8: satflag (卫星校正标识)
%   列9: glacier_flag (冰川标识, 0=冰川内, 1=冰川外)
%   列10: year (年份)
%   列11: hma_region (HMA区域编号)
%   列12: hma22_region (HMA22区域编号)
%   列13: gtn_region (GTN15区域编号)
%   列14: grid_region (grid496区域编号)
warning('off', 'map:removing:latlon2pix');
fprintf('=== ICESat GLAH14数据提取 ===\n');

%% 初始化
track_data = [];
track_data_outrgi = [];

% 检查输入目录
if ~exist(config.paths.input_dir, 'dir')
    error('输入目录不存在: %s', config.paths.input_dir);
end

%% 加载地形数据
fprintf('加载地形数据...\n');
try
    [srtm, ref_srtm] = geotiffread(config.files.nasadem);
    [slope_data, ref_slope] = geotiffread(config.files.slope);
    [aspect_data, ref_aspect] = geotiffread(config.files.aspect);
    [glacier_mask, ref_glacier] = geotiffread(config.files.glacier_mask);
    
    fprintf('✓ 地形数据加载完成\n');
catch ME
    error('地形数据加载失败: %s', ME.message);
end

%% 加载区域shapefile
fprintf('加载区域shapefile...\n');
try
    HMA = shaperead(config.files.hma_boundary);
    HMA22 = shaperead(config.files.hma22_regions);
    GTN = shaperead(config.files.gtn_regions);
    grid10 = shaperead(config.files.grid10);
    
    fprintf('✓ 区域数据加载完成\n');
catch ME
    error('区域数据加载失败: %s', ME.message);
end

%% 查找HDF文件
fprintf('搜索HDF文件...\n');
hdf_files = dir(fullfile(config.paths.input_dir, '**', 'GLAH14*.H5'));

if isempty(hdf_files)
    warning('未找到GLAH14 HDF文件');
    fprintf('请检查输入目录: %s\n', config.paths.input_dir);
    fprintf('预期文件格式: GLAH14*.H5\n');
    return;
end

fprintf('✓ 找到 %d 个HDF文件\n', length(hdf_files));

%% 提取数据
fprintf('开始提取数据...\n');

for file_idx = 1:length(hdf_files)
    hdf_file = fullfile(hdf_files(file_idx).folder, hdf_files(file_idx).name);
    
    fprintf('  处理文件 %d/%d: %s\n', file_idx, length(hdf_files), hdf_files(file_idx).name);

    % 打开HDF5文件并读取数据
    latitude = h5read(hdf_file, '/Data_40HZ/Geolocation/d_lat');
    longtitude = h5read(hdf_file, '/Data_40HZ/Geolocation/d_lon');
    geoid = h5read(hdf_file, '/Data_40HZ/Geophysical/d_gdHt');
    time_tai = h5read(hdf_file, '/Data_40HZ/DS_UTCTime_40'); %seconds since 2000-1-1
    time_tai = round(time_tai/86400, 0);  % 转换为天数
    elevation = h5read(hdf_file, '/Data_40HZ/Elevation_Surfaces/d_elev');
    satcorr = h5read(hdf_file, '/Data_40HZ/Elevation_Corrections/d_satElevCorr');
    satflag = h5read(hdf_file, '/Data_40HZ/Quality/sat_corr_flg');
    angles = h5read(hdf_file, '/Data_40HZ/Elevation_Angles/d_beamAzimuth');
    
    % 初始化数据矩阵（19列）
    point = zeros(numel(latitude), 14);
    point(:,1) = longtitude;  % 列1: lon
    point(:,2) = latitude;    % 列2: lat
    point(:,3) = elevation;
    point(:,4) = angles;
    point(:,5) = time_tai; %days
    point(:,6) = geoid;
    point(:,7) = satcorr;
    point(:,8) = satflag;
    
    % 边界框筛选（扩展0.1度）
    lat_min = 25 - 0.1;
    lat_max = 47 + 0.1;
    lon_min = 66 - 0.1;
    lon_max = 106 + 0.1;
    
    toDelete = (point(:,1) < lon_min) | (point(:,1) > lon_max) | ...
               (point(:,2) < lat_min) | (point(:,2) > lat_max);
    point(toDelete,:) = [];
    
    if isempty(point)
        continue;
    end
    
    % 计算年份
    point(:,10) = year(datetime(2000,1,1) + days(point(:,5)));
    
    % 对每个点添加地形信息
    for ii = 1:size(point, 1)
        % 检查是否在冰川内
        [m1, n1] = size(glacier_mask);
        [row1, col1] = latlon2pix(ref_glacier, point(ii,2), point(ii,1));
        row1 = ceil(row1);
        col1 = ceil(col1);
        
        if row1 > m1 || row1 < 1 || col1 > n1 || col1 < 1
            point(ii,9) = 1;  % 超出范围，标记为非冰川
            continue;
        end
        
        % 检查冰川掩膜值
        if glacier_mask(row1, col1) == -2147483647 || glacier_mask(row1, col1) == 65535
            point(ii,9) = 1;  % 非冰川
        else
            point(ii,9) = 0;  % 冰川内
        end
    end
    
    % 筛选出冰川与非冰川点
    mask = (point(:,9) == 0);
    glacier_point = add_spatial_flag(point, mask, HMA, HMA22, GTN, grid10);

    mask = (point(:,9) == 1);
    out_glacier_point = add_spatial_flag(point, mask, HMA, HMA22, GTN, grid10);
    
    if ~isempty(glacier_point)
        track_data = [track_data; glacier_point];
    end
    
    if ~isempty(out_glacier_point)
        track_data_outrgi = [track_data_outrgi; out_glacier_point];
    end
end

fprintf('✓ 数据提取完成，共 %d 个冰川数据点\n', size(track_data, 1));
fprintf('✓ 数据提取完成，共 %d 个非冰川数据点\n', size(track_data_outrgi, 1));

end

function spatial_point = add_spatial_flag(point, mask, HMA, HMA22, gtn, grid10)
    spatial_point = point(mask, :);

    if isempty(spatial_point)
        return;
    end

    % 添加区域标识
    for k = 1:length(HMA)
        [in, ~] = inpolygon(spatial_point(:,1), spatial_point(:,2), HMA(k).X, HMA(k).Y);
        spatial_point(in==1, 11) = k;
    end
    
    for k = 1:length(HMA22)
        [in, ~] = inpolygon(spatial_point(:,1), spatial_point(:,2), HMA22(k).X, HMA22(k).Y);
        spatial_point(in==1, 12) = k;
    end

    for k = 1:length(gtn)
        [in, ~] = inpolygon(spatial_point(:,1), spatial_point(:,2), gtn(k).X, gtn(k).Y);
        spatial_point(in==1, 13) = k;
    end

    for k = 1:length(grid10)
        [in, ~] = inpolygon(spatial_point(:,1), spatial_point(:,2), grid10(k).X, grid10(k).Y);
        spatial_point(in==1, 14) = k;
    end
end

