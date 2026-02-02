function extract_nonrgi_icesat2_tracks(config)
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

% 关闭HDF5警告
warning('off','MATLAB:imagesci:hdf5dataset:libraryError');
warning('off','MATLAB:imagesci:h5info:libraryError');
warning('off','MATLAB:imagesci:geotiffread:libraryError');

fprintf('开始提取ICESat-2轨道数据...\n');

%% 加载空间参考数据
fprintf('加载空间参考数据...\n');

% 加载边界和网格数据
grid10 = shaperead(config.files.grid10);
fprintf('空间参考数据加载完成。\n');

%% 获取HDF5文件列表
h5_files = dir(fullfile(config.paths.input_h5_dir, '*.h5'));
fprintf('找到 %d 个HDF5文件需要处理。\n', length(h5_files));

%% 定义激光轨道
ground_tracks = ['gt1l'; 'gt1r'; 'gt2l'; 'gt2r'; 'gt3l'; 'gt3r'];

%% 逐文件处理
parfor i = 1:length(h5_files)
    process_single_h5_file(h5_files(i), config.paths.input_h5_dir, ...
        ground_tracks, grid10, config, i);
    
end

end

function process_single_h5_file(h5_file, input_dir, ground_tracks, ...
    grid10, config, i)
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
% 获取HDF5文件信息
info = h5info(file_path);
num_groups = length(info.Groups);

% 处理6个激光轨道（跳过前后两个元数据组）
for j = 3:(num_groups-2)
    track_data = extract_single_track(file_path, info.Groups(j).Name, config);
    
    if isempty(track_data)
        continue;
    end
    
    file_data = [file_data; track_data];
end

if size(file_data, 1)>100
    for grid_i = 1:length(grid10)
        lon_indicator = floor(grid10(grid_i).x);
        lat_indicator = floor(grid10(grid_i).y);
    
        rgi_filename = sprintf('N%02dE%03d_rgi.shp', lat_indicator, lon_indicator);
        grid_rgi_path = fullfile(config.paths.glacier_mask_dir, rgi_filename);
        
        file_data_temp = file_data((file_data(:,1) >= lon_indicator) & (file_data(:,1) < lon_indicator+1) ...
            & (file_data(:,2) >= lat_indicator) & (file_data(:,2) < lat_indicator+1), :);
        if isempty(file_data_temp)
            continue
        end

        if exist(grid_rgi_path, "file")
            grid_rgi = shaperead(grid_rgi_path);

            lon_points = file_data_temp(:,1);
            lat_points = file_data_temp(:,2);

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

            is_inside_any_glacier = inpolygon(lon_points, lat_points, all_glacier_poly_x, all_glacier_poly_y);

            file_data_temp = file_data_temp(~is_inside_any_glacier, :);
        end
        
        if ~isempty(file_data_temp)
            save_path = fullfile(config.paths.foot_out_glacier_dir, sprintf('N%02dE%03d_footprint_outrgi.csv', lat_indicator, lon_indicator));
            writematrix(file_data_temp, save_path, 'WriteMode', 'append');
            log_fid = fopen(config.files.process_log, 'at');
            fprintf(log_fid, '【%s】file 【%04d】 %s; %d footprints save to %s\n', datetime('now'), i, ...
                h5_file.name, size(file_data_temp, 1), ...
                sprintf('N%02dE%03d_footprint_outrgi.csv', lat_indicator, lon_indicator));
            fclose(log_fid);
        else
            log_fid = fopen(config.files.process_log, 'at');
            fprintf(log_fid, '【%s】file 【%04d】 %s no save', datetime('now'), i, h5_file.name);
            fclose(log_fid);
        end
    end
else
    log_fid = fopen(config.files.process_log, 'at');
    fprintf(log_fid, '【%s】file 【%04d】 %s no save', datetime('now'), i, h5_file.name);
    fclose(log_fid);
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

end