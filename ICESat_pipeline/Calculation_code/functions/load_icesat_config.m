function config = load_icesat_config()
% LOAD_ICESAT_CONFIG Load configuration parameters for ICESat data processing
%
% Output:
%   config - structure containing all processing parameters
%
% Example:
%   config = load_icesat_config();

%% Path configuration
config.paths.input_dir = 'D:\HMA_ICESat\file';                  % Directory of original ICESat data
config.paths.output_dir = 'E:\revised_NCC_data\ICESat_result';          % Output results directory
config.paths.shapefile_dir = 'E:\HMA_subregion';           % Shapefile directory
config.paths.dem_dir = 'E:\HMA';                           % DEM data directory
config.paths.temp_dir = 'temp';                            % Temporary file directory
config.paths.ele_band_dir = 'E:\revised_NCC_data\ICESat_result\ele_band_results';  % Elevation-band results directory

%% Input data files
config.files.hma_boundary = 'E:\HMA_subregion\HMA_boundary\HMA.shp';
config.files.hma22_regions = 'E:\HMA_subregion\regions_hma_v03_zheng\boundary_mountain_regions_hma_v3_zheng_20200601.shp';
config.files.gtn_regions = 'E:\HMA_subregion\GlacReg_2017\GTN_HMA.shp';
config.files.tp_boundary = 'E:\HMA_subregion\TP_Boundary\TP_boundary.shp';
config.files.hma4_regions = 'E:\HMA_subregion\HMA_four_parts.shp';
config.files.hma6_regions = 'E:\HMA_subregion\HMA_seasonal_pattern.shp';
config.files.grid10 = 'E:\revised_NCC_data\grid10_pdd\HMA_fishnet_10_clip_coreg_filled.shp';  % 1° grid
config.files.glacier_mask = 'E:\HMA\glacier_mask\rgi60_hma.tif';
config.files.nasadem = 'E:\HMA\NASADEM_HMA.tif';
config.files.slope = 'E:\HMA\NASADEM_HMA_slope.tif';
config.files.aspect = 'E:\HMA\NASADEM_HMA_aspect.tif';
config.files.elevation_area = 'E:\HMA\HMA高度带.xlsx';
config.files.process_log = 'process_log.txt';

%% Output file names
config.files.track_txt = 'HMA_ICESat_rgiregion.txt';           % Extracted track data
config.files.enhanced_txt = 'HMA_ICESAT_rgiregion_update.txt'; % Enhanced data
config.files.results_mat = 'ICESat_analysis_results.mat';      % Analysis results

%% Processing parameters
% Time parameters
config.time.reference_date = datetime(2000, 1, 1);     % Reference time (ICESat uses 2000-01-01)
config.time.start_date = datetime(2003, 1, 1);         % Start time
config.time.end_date = datetime(2009, 12, 1);         % End time

% Quality control parameters
config.quality.max_elevation = 8900;                   % Maximum elevation threshold (m)
config.quality.max_height_change = 200;                % Maximum elevation change threshold (m) - first screening
config.quality.max_slope = 40;                         % Maximum slope threshold (degrees)
config.quality.height_change_filter1 = 100;            % First-round elevation-change filtering (m)
config.quality.height_change_filter2 = 50;             % Second-round elevation-change filtering (m) - relative to the median
config.quality.outlier_percentile_low = 1;             % Lower percentile for outliers
config.quality.outlier_percentile_high = 99;           % Upper percentile for outliers
config.quality.max_swb = 255;                          % Maximum swb threshold (reserved, not actually used)

% Spatial parameters
config.spatial.grid_resolution = 1.0;                  % Grid resolution (degrees) - 1° grid
config.spatial.elevation_band_width = 100;              % Elevation-band width (m)
config.spatial.elevation_range = [2250, 8750];         % Elevation range (m)
config.spatial.nasadem_neighbor_size = 3;              % NASADEM neighborhood size (ICESat uses 3x3)

% Temporal analysis parameters
config.temporal.monthly_analysis = true;               % Monthly analysis switch
config.temporal.annual_analysis = true;                % Annual analysis switch
config.temporal.start_year = 2003;                     % Start year
config.temporal.end_year = 2009;                       % End year

% ICESat-specific time windows (campaigns)
% Defined according to HMA_ICESat_COMPUTE_20230603.m
config.temporal.campaigns = [
    2003.191781, datenum(2003,02,20)-datenum(2000,1,1), datenum(2003,03,29)-datenum(2000,1,1)+1;
    2003.810959, datenum(2003,09,25)-datenum(2000,1,1), datenum(2003,11,19)-datenum(2000,1,1)+1;
    2004.178082, datenum(2004,02,17)-datenum(2000,1,1), datenum(2004,03,21)-datenum(2000,1,1)+1;
    2004.428767, datenum(2004,05,18)-datenum(2000,1,1), datenum(2004,06,21)-datenum(2000,1,1)+1;
    2004.809589, datenum(2004,10,03)-datenum(2000,1,1), datenum(2004,11,08)-datenum(2000,1,1)+1;
    2005.183562, datenum(2005,02,17)-datenum(2000,1,1), datenum(2005,03,24)-datenum(2000,1,1)+1;
    2005.434247, datenum(2005,05,20)-datenum(2000,1,1), datenum(2005,06,23)-datenum(2000,1,1)+1;
    2005.856164, datenum(2005,10,21)-datenum(2000,1,1), datenum(2005,11,24)-datenum(2000,1,1)+1;
    2006.195890, datenum(2006,02,22)-datenum(2000,1,1), datenum(2006,03,28)-datenum(2000,1,1)+1;
    2006.443836, datenum(2006,05,24)-datenum(2000,1,1), datenum(2006,06,26)-datenum(2000,1,1)+1;
    2006.865753, datenum(2006,10,25)-datenum(2000,1,1), datenum(2006,11,27)-datenum(2000,1,1)+1;
    2007.243836, datenum(2007,03,12)-datenum(2000,1,1), datenum(2007,04,14)-datenum(2000,1,1)+1;
    2007.536986, datenum(2007,05,15)-datenum(2000,1,1), datenum(2007,09,11)-datenum(2000,1,1)+1;
    2007.804110, datenum(2007,10,02)-datenum(2000,1,1), datenum(2007,11,05)-datenum(2000,1,1)+1;
    2008.180822, datenum(2008,02,17)-datenum(2000,1,1), datenum(2008,03,21)-datenum(2000,1,1)+1;
    2008.786301, datenum(2008,10,04)-datenum(2000,1,1), datenum(2008,10,19)-datenum(2000,1,1)+1;
    2008.938356, datenum(2008,11,25)-datenum(2000,1,1), datenum(2008,12,17)-datenum(2000,1,1)+1;
    2009.238356, datenum(2009,03,09)-datenum(2000,1,1), datenum(2009,04,11)-datenum(2000,1,1)+1;
    2009.769863, datenum(2009,09,30)-datenum(2000,1,1), datenum(2009,10,11)-datenum(2000,1,1)+1;
];

%% Processing options
config.options.parallel_processing = false;            % Parallel processing switch
config.options.save_intermediate = true;               % Save intermediate results
config.options.generate_plots = true;                  % Generate plots
config.options.generate_report = true;                 % Generate report
config.options.verbose = true;                         % Verbose output

%% Validate configuration
config = validate_config(config);

fprintf('ICESat配置加载完成。\n');

end

function config = validate_config(config)
% VALIDATE_CONFIG Validate configuration parameter correctness

% Create output directory
if ~exist(config.paths.output_dir, 'dir')
    mkdir(config.paths.output_dir);
    fprintf('创建输出目录: %s\n', config.paths.output_dir);
end

% Create temporary directory
temp_full_path = fullfile(config.paths.output_dir, config.paths.temp_dir);
if ~exist(temp_full_path, 'dir')
    mkdir(temp_full_path);
end
config.paths.temp_dir = temp_full_path;

% Create elevation-band results directory
if ~exist(config.paths.ele_band_dir, 'dir')
    mkdir(config.paths.ele_band_dir);
end

% Validate time parameters
if config.time.start_date >= config.time.end_date
    error('开始时间必须早于结束时间');
end

fprintf('配置验证通过。\n');

end

