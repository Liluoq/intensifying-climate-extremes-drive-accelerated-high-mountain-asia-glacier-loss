function config = load_processing_config()
% LOAD_PROCESSING_CONFIG 加载ICESat-2数据处理的配置参数
%
% 输出:
%   config - 包含所有处理参数的结构体
%
% 示例:
%   config = load_processing_config();

%% 路径配置
config.paths.input_h5_dir = 'D:\HMA_ICESat2\h5file';  % ICESat-2 HDF5文件目录
config.paths.output_dir = 'E:\revised_NCC_data\ICESat2_result';    % 输出结果目录
config.paths.shapefile_dir = 'E:\HMA_subregion';      % Shapefile目录
config.paths.dem_dir = 'E:\HMA';                      % DEM数据目录
config.paths.temp_dir = 'temp';                       % 临时文件目录
config.paths.ele_band_dir = 'E:\revised_NCC_data\ICESat2_result\ele_band_results';   %存储各个区域的高度带结果

%% 输入数据文件
config.files.hma_boundary = 'E:\HMA_subregion\HMA_boundary\HMA.shp';
config.files.hma22_regions = 'E:\HMA_subregion\regions_hma_v03_zheng\boundary_mountain_regions_hma_v3_zheng_20200601.shp';
config.files.gtn_regions = 'E:\HMA_subregion\GlacReg_2017\GTN_HMA.shp';
config.files.tp_boundary = 'E:\HMA_subregion\TP_Boundary\TP_boundary.shp';
config.files.hma_four_parts = 'E:\HMA_subregion\HMA_four_parts.shp';
config.files.hma_seasonal_pattern = 'E:\HMA_subregion\HMA_seasonal_pattern.shp';
config.files.grid05 = 'E:\HMA\HMA_boundary\HMA_fishnet_05_clip.shp';
config.files.grid10 = 'E:\revised_NCC_data\grid10_pdd\HMA_fishnet_10_clip_coreg_filled.shp';
config.files.glacier_mask = 'E:\HMA\glacier_mask\rgi60_hma.tif';
config.files.nasadem = 'E:\HMA\NASADEM_HMA.tif';
config.files.slope = 'E:\HMA\NASADEM_HMA_slope.tif';
config.files.aspect = 'E:\HMA\NASADEM_HMA_aspect.tif';
config.files.elevation_area = 'E:\HMA\HMA高度带.xlsx';

%% 输出文件名
config.files.track_csv = 'HMA_ICESat2_glacier_V6.csv';
config.files.enhanced_txt = 'HMA_ICESat2_glacier_V6.txt';
config.files.results_mat = 'ICESat2_analysis_results.mat';

%% 处理参数
% 时间参数
config.time.reference_date = datetime(2018, 1, 1);     % 参考时间
config.time.start_date = datetime(2018, 10, 1);        % 开始时间
config.time.end_date = datetime(2025, 10, 1);          % 结束时间

% 质量控制参数
config.quality.max_elevation = 8900;                   % 最大高程阈值 (m)
config.quality.max_height_change = 200;                % 最大高程变化阈值 (m)
config.quality.max_slope = 40;                         % 最大坡度阈值 (度)
config.quality.max_curvature = 4;                      % 最大曲率阈值
config.quality.height_change_filter1 = 150;            % 第一轮高程变化过滤 (m)
config.quality.height_change_filter2 = 75;             % 第二轮高程变化过滤 (m)
config.quality.outlier_percentile = 1;                 % 异常值百分位数

% 空间参数
config.spatial.grid_resolution = [0.5, 1.0];           % 网格分辨率 (度)
config.spatial.elevation_band_width = [50, 100];       % 高程带宽度 (m)
config.spatial.elevation_range = [2250, 8750];         % 高程范围 (m)

% 时间分析参数
config.temporal.monthly_analysis = true;               % 月度分析开关
config.temporal.seasonal_analysis = true;              % 季节分析开关
config.temporal.annual_analysis = true;                % 年度分析开关
config.temporal.cycle_months = 3;                      % 循环分析月数

%% 处理选项
config.options.parallel_processing = false;            % 并行处理开关
config.options.save_intermediate = true;               % 保存中间结果
config.options.generate_plots = true;                  % 生成图表
config.options.generate_report = true;                 % 生成报告
config.options.verbose = true;                         % 详细输出

%% 验证配置
config = validate_config(config);

fprintf('配置验证完成。\n');

end

function config = validate_config(config)
% VALIDATE_CONFIG 验证配置参数的有效性
%
% 输入:
%   config - 配置结构体
% 输出:
%   config - 验证后的配置结构体

% 创建输出目录
if ~exist(config.paths.output_dir, 'dir')
    mkdir(config.paths.output_dir);
    fprintf('创建输出目录: %s\n', config.paths.output_dir);
end

% 创建临时目录
temp_full_path = fullfile(config.paths.output_dir, config.paths.temp_dir);
if ~exist(temp_full_path, 'dir')
    mkdir(temp_full_path);
end
config.paths.temp_dir = temp_full_path;

% 验证关键文件存在
critical_files = {
    config.files.hma_boundary,
    config.files.hma22_regions,
    config.files.gtn_regions,
    config.files.glacier_mask
};

for i = 1:length(critical_files)
    if ~exist(critical_files{i}, 'file')
        warning('关键文件不存在: %s', critical_files{i});
    end
end

% 验证时间参数
if config.time.start_date >= config.time.end_date
    error('开始时间必须早于结束时间');
end

fprintf('配置验证通过。\n');

end
