function extract_nonrgi_cryosat2_tracks(config)
% EXTRACT_NONRGI_CRYOSAT2_TRACKS 提取CryoSat-2非RGI冰川区域数据
%
% 参考：revised_ICESat2_code/functions/extract_nonrgi_icesat2_tracks.m
%
% 功能：提取RGI冰川区域外的数据，用于对比研究
%
% 输入:
%   config - 配置结构体
%
% 输出:
% 第1列 lon
% 第2列 lat
% 第3 ele
% 第4 days, since 2000-01-01-00:00:00
% 第5 year
% 第6 month
% 第7 geoid %EGM96
% 第8 uncert
%   按1°网格分块保存为CSV文件

warning('off','MATLAB:imagesci:geotiffread:libraryError');

fprintf('开始提取CryoSat-2非冰川区域数据...\n');

%% 加载空间参考数据
fprintf('加载空间参考数据...\n');
grid10 = shaperead(config.files.grid10);
fprintf('空间参考数据加载完成。\n');

%% 获取NetCDF文件列表
path = config.paths.input_tempo_dir;
file_list = cell(length(path), 1);  % 预分配cell数组
for i = 1:length(path)
    file_list{i} = dir(fullfile(path{i}, '*.nc'));
end

nc_files = vertcat(file_list{:});% struct形式
fprintf('找到 %d 个NetCDF文件需要处理。\n', length(nc_files));

%% 并行处理文件
for i = 1:length(nc_files)
    process_single_nc_file(nc_files(i), grid10, config, i);
end

fprintf('✓ 非冰川区域数据提取完成\n');

end

function process_single_nc_file(nc_file, grid10, config, i)
% PROCESS_SINGLE_NC_FILE 处理单个NetCDF文件

file_path = fullfile(nc_file.folder, nc_file.name);
file_data = [];

try
    % 读取NetCDF数据
    lon = ncread(file_path, 'lon');
    lat = ncread(file_path, 'lat');
    elevation = ncread(file_path, 'elevation');
    time_total = ncread(file_path, 'time');
    time_total = floor(time_total / 86400);% seconds->days, since 2000-01-01-00:00:00
    uncertainty = ncread(file_path, 'uncertainty');
    
    % 从文件名提取年份
    yr = nc_file.name(30:33);
    mon = nc_file.name(35:36);
    geoid = geoidheight(lat, lon);%EGM96
    
    % 组织数据
    file_data = zeros(length(lon), 8);
    file_data(:,1) = lon;
    file_data(:,2) = lat;
    file_data(:,3) = elevation;
    file_data(:,4) = time_total; %days, since 2000-01-01-00:00:00
    file_data(:,5) = ones(length(lat),1)*str2num(yr);
    file_data(:,6) = ones(length(lat),1)*str2num(mon);
    file_data(:,7) = geoid;
    file_data(:,8) = uncertainty;
    
catch ME
    log_fid = fopen(config.files.process_log, 'at');
    fprintf(log_fid, '【%s】file 【%05d】 %s read failed: %s\n', ...
        datetime('now'), i, nc_file.name, ME.message);
    fclose(log_fid);
    return;
end

% 如果数据点太少，跳过
if size(file_data, 1) < 1
    log_fid = fopen(config.files.process_log, 'at');
    fprintf(log_fid, '【%s】file 【%05d】 %s too few points\n', ...
        datetime('now'), i, nc_file.name);
    fclose(log_fid);
    return;
end

%% 按网格处理
for grid_i = 1:length(grid10)
    lon_indicator = floor(grid10(grid_i).x);
    lat_indicator = floor(grid10(grid_i).y);
    
    rgi_filename = sprintf('N%02dE%03d_rgi.shp', lat_indicator, lon_indicator);
    grid_rgi_path = fullfile(config.paths.glacier_mask_dir, rgi_filename);
    
    % 筛选该网格的数据
    file_data_temp = file_data((file_data(:,1) >= lon_indicator) & ...
                               (file_data(:,1) < lon_indicator+1) & ...
                               (file_data(:,2) >= lat_indicator) & ...
                               (file_data(:,2) < lat_indicator+1), :);
    
    if isempty(file_data_temp)
        log_fid = fopen(config.files.process_log, 'at');
        fprintf(log_fid, '【%s】file 【%05d】 %s no save\n', datetime('now'), i, nc_file.name);
        fclose(log_fid);
        continue
    end
    
    % 如果存在冰川掩膜，排除冰川区域
    if exist(grid_rgi_path, "file")
        grid_rgi = shaperead(grid_rgi_path);
        
        lon_points = file_data_temp(:,1);
        lat_points = file_data_temp(:,2);
        
        % 合并所有冰川多边形
        all_glacier_poly_x_cells = cell(1, 2*length(grid_rgi) - 1);
        all_glacier_poly_y_cells = cell(1, 2*length(grid_rgi) - 1);
        
        for k_poly = 1:length(grid_rgi)
            all_glacier_poly_x_cells{2*k_poly-1} = grid_rgi(k_poly).X;
            all_glacier_poly_y_cells{2*k_poly-1} = grid_rgi(k_poly).Y;
            if k_poly < length(grid_rgi)
                all_glacier_poly_x_cells{2*k_poly} = NaN;
                all_glacier_poly_y_cells{2*k_poly} = NaN;
            end
        end
        
        all_glacier_poly_x = cell2mat(all_glacier_poly_x_cells);
        all_glacier_poly_y = cell2mat(all_glacier_poly_y_cells);
        
        % 排除冰川区域内的点
        is_inside_any_glacier = inpolygon(lon_points, lat_points, ...
            all_glacier_poly_x, all_glacier_poly_y);
        
        file_data_temp = file_data_temp(~is_inside_any_glacier, :);
    end
    
    % 保存数据
    if ~isempty(file_data_temp)
        save_path = fullfile(config.paths.foot_out_glacier_dir, ...
            sprintf('N%02dE%03d_footprint_outrgi.csv', lat_indicator, lon_indicator));
        writematrix(file_data_temp, save_path, 'WriteMode', 'append');
        
        log_fid = fopen(config.files.process_log, 'at');
        fprintf(log_fid, '【%s】file 【%05d】 %s; %d footprints save to %s\n', ...
            datetime('now'), i, nc_file.name, size(file_data_temp, 1), ...
            sprintf('N%02dE%03d_footprint_outrgi.csv', lat_indicator, lon_indicator));
        fclose(log_fid);
    else
        log_fid = fopen(config.files.process_log, 'at');
        fprintf(log_fid, '【%s】file 【%05d】 %s no save\n', datetime('now'), i, nc_file.name);
        fclose(log_fid);
    end
end

end

