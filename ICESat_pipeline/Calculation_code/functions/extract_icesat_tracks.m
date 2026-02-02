function [track_data, track_data_outrgi] = extract_icesat_tracks(config)
% EXTRACT_ICESAT_TRACKS Extract track data from ICESat GLAH14 HDF files
%
% Function: extract track data within RGI glacier regions from ICESat GLAH14 products (HDF format)
% Reference: HMA_Icesattrack_full.m
%
% Input:
%   config - configuration structure
%
% Output:
%   track_data - extracted track data matrix (14 columns)
%
% Output column definitions (note: column 1 is lon, column 2 is lat):
%   Col 1: lon (longitude)
%   Col 2: lat (latitude)
%   Col 3: elevation (observed elevation)
%   Col 4: beam_azimuth (beam azimuth)
%   Col 5: days (days since 2000-01-01)
%   Col 6: geoid (EGM96 geoid height)
%   Col 7: satcorr (satellite elevation correction)
%   Col 8: satflag (satellite correction flag)
%   Col 9: glacier_flag (glacier flag, 0 = inside glacier, 1 = outside glacier)
%   Col10: year
%   Col11: hma_region (HMA region ID)
%   Col12: hma22_region (HMA22 region ID)
%   Col13: gtn_region (GTN15 region ID)
%   Col14: grid_region (grid496 region ID)
warning('off', 'map:removing:latlon2pix');
fprintf('=== ICESat GLAH14数据提取 ===\n');

%% Initialization
track_data = [];
track_data_outrgi = [];

% Check input directory
if ~exist(config.paths.input_dir, 'dir')
    error('输入目录不存在: %s', config.paths.input_dir);
end

%% Load terrain data
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

%% Load region shapefiles
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

%% Search for HDF files
fprintf('搜索HDF文件...\n');
hdf_files = dir(fullfile(config.paths.input_dir, '**', 'GLAH14*.H5'));

if isempty(hdf_files)
    warning('未找到GLAH14 HDF文件');
    fprintf('请检查输入目录: %s\n', config.paths.input_dir);
    fprintf('预期文件格式: GLAH14*.H5\n');
    return;
end

fprintf('✓ 找到 %d 个HDF文件\n', length(hdf_files));

%% Extract data
fprintf('开始提取数据...\n');

for file_idx = 1:length(hdf_files)
    hdf_file = fullfile(hdf_files(file_idx).folder, hdf_files(file_idx).name);
    
    fprintf('  处理文件 %d/%d: %s\n', file_idx, length(hdf_files), hdf_files(file_idx).name);

    % Open HDF5 file and read data
    latitude = h5read(hdf_file, '/Data_40HZ/Geolocation/d_lat');
    longtitude = h5read(hdf_file, '/Data_40HZ/Geolocation/d_lon');
    geoid = h5read(hdf_file, '/Data_40HZ/Geophysical/d_gdHt');
    time_tai = h5read(hdf_file, '/Data_40HZ/DS_UTCTime_40'); %seconds since 2000-1-1
    time_tai = round(time_tai/86400, 0);  % 转换为天数
    elevation = h5read(hdf_file, '/Data_40HZ/Elevation_Surfaces/d_elev');
    satcorr = h5read(hdf_file, '/Data_40HZ/Elevation_Corrections/d_satElevCorr');
    satflag = h5read(hdf_file, '/Data_40HZ/Quality/sat_corr_flg');
    angles = h5read(hdf_file, '/Data_40HZ/Elevation_Angles/d_beamAzimuth');
    
    % Initialize data matrix (19 columns)
    point = zeros(numel(latitude), 14);
    point(:,1) = longtitude;  % Col 1: lon
    point(:,2) = latitude;    % Col 2: lat
    point(:,3) = elevation;
    point(:,4) = angles;
    point(:,5) = time_tai; % days
    point(:,6) = geoid;
    point(:,7) = satcorr;
    point(:,8) = satflag;
    
    % Bounding-box filter (extended by 0.1°)
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
    
    % Compute year
    point(:,10) = year(datetime(2000,1,1) + days(point(:,5)));
    
    % Add terrain information for each point
    for ii = 1:size(point, 1)
        % Check whether the point lies within a glacier
        [m1, n1] = size(glacier_mask);
        [row1, col1] = latlon2pix(ref_glacier, point(ii,2), point(ii,1));
        row1 = ceil(row1);
        col1 = ceil(col1);
        
        if row1 > m1 || row1 < 1 || col1 > n1 || col1 < 1
            point(ii,9) = 1;  % Out of bounds, mark as non-glacier
            continue;
        end
        
        % Check glacier mask value
        if glacier_mask(row1, col1) == -2147483647 || glacier_mask(row1, col1) == 65535
            point(ii,9) = 1;  % Non-glacier
        else
            point(ii,9) = 0;  % Inside glacier
        end
    end
    
    % Separate glacier and non-glacier points
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

    % Add region identifiers
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

