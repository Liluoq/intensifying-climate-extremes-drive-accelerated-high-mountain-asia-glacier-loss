function config = load_icesat_config()
% LOAD_ICESAT_CONFIG 加载ICESat数据处理的配置参数
%
% 输出:
%   config - 包含所有处理参数的结构体
%
% 示例:
%   config = load_icesat_config();

%% 路径配置
config.paths.input_dir = 'D:\HMA_ICESat\file';                  % ICESat原始数据目录
config.paths.output_dir = 'E:\revised_NCC_data\ICESat_result';          % 输出结果目录
config.paths.shapefile_dir = 'E:\HMA_subregion';           % Shapefile目录
config.paths.dem_dir = 'E:\HMA';                           % DEM数据目录
config.paths.temp_dir = 'temp';                            % 临时文件目录
config.paths.ele_band_dir = 'E:\revised_NCC_data\ICESat_result\ele_band_results';  % 高程带结果目录

%% 输入数据文件
config.files.hma_boundary = 'E:\HMA_subregion\HMA_boundary\HMA.shp';
config.files.hma22_regions = 'E:\HMA_subregion\regions_hma_v03_zheng\boundary_mountain_regions_hma_v3_zheng_20200601.shp';
config.files.gtn_regions = 'E:\HMA_subregion\GlacReg_2017\GTN_HMA.shp';
config.files.tp_boundary = 'E:\HMA_subregion\TP_Boundary\TP_boundary.shp';
config.files.hma4_regions = 'E:\HMA_subregion\HMA_four_parts.shp';
config.files.hma6_regions = 'E:\HMA_subregion\HMA_seasonal_pattern.shp';
config.files.grid10 = 'E:\revised_NCC_data\grid10_pdd\HMA_fishnet_10_clip_coreg_filled.shp';  % 1°网格
config.files.glacier_mask = 'E:\HMA\glacier_mask\rgi60_hma.tif';
config.files.nasadem = 'E:\HMA\NASADEM_HMA.tif';
config.files.slope = 'E:\HMA\NASADEM_HMA_slope.tif';
config.files.aspect = 'E:\HMA\NASADEM_HMA_aspect.tif';
config.files.elevation_area = 'E:\HMA\HMA高度带.xlsx';
config.files.process_log = 'process_log.txt';

%% 输出文件名
config.files.track_txt = 'HMA_ICESat_rgiregion.txt';           % 提取的轨道数据
config.files.enhanced_txt = 'HMA_ICESAT_rgiregion_update.txt'; % 增强后的数据
config.files.results_mat = 'ICESat_analysis_results.mat';      % 分析结果

%% 处理参数
% 时间参数
config.time.reference_date = datetime(2000, 1, 1);     % 参考时间（ICESat使用2000-01-01）
config.time.start_date = datetime(2003, 1, 1);         % 开始时间
config.time.end_date = datetime(2009, 12, 1);         % 结束时间

% 质量控制参数
config.quality.max_elevation = 8900;                   % 最大高程阈值(m)
config.quality.max_height_change = 200;                % 最大高程变化阈值(m) - 第一轮筛选
config.quality.max_slope = 40;                         % 最大坡度阈值(度)
config.quality.height_change_filter1 = 100;            % 第一轮高程变化过滤(m)
config.quality.height_change_filter2 = 50;             % 第二轮高程变化过滤(m) - 相对中值
config.quality.outlier_percentile_low = 1;             % 异常值下百分位数
config.quality.outlier_percentile_high = 99;           % 异常值上百分位数
config.quality.max_swb = 255;                          % 最大swb阈值（预留，实际未使用）

% 空间参数
config.spatial.grid_resolution = 1.0;                  % 网格分辨率(度) - 1°网格
config.spatial.elevation_band_width = 100;              % 高程带宽度(m)
config.spatial.elevation_range = [2250, 8750];         % 高程范围(m)
config.spatial.nasadem_neighbor_size = 3;              % NASADEM邻域大小 (ICESat使用3x3)

% 时间分析参数
config.temporal.monthly_analysis = true;               % 月度分析开关
config.temporal.annual_analysis = true;                % 年度分析开关
config.temporal.start_year = 2003;                     % 起始年份
config.temporal.end_year = 2009;                       % 结束年份

% ICESat特定的时间窗口（campaigns）
% 根据HMA_ICESat_COMPUTE_20230603.m中的定义
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

%% 处理选项
config.options.parallel_processing = false;            % 并行处理开关
config.options.save_intermediate = true;               % 保存中间结果
config.options.generate_plots = true;                  % 生成图表
config.options.generate_report = true;                 % 生成报告
config.options.verbose = true;                         % 详细输出

%% 验证配置
config = validate_config(config);

fprintf('ICESat配置加载完成。\n');

end

function config = validate_config(config)
% VALIDATE_CONFIG 验证配置参数的有效性

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

% 创建高程带结果目录
if ~exist(config.paths.ele_band_dir, 'dir')
    mkdir(config.paths.ele_band_dir);
end

% 验证时间参数
if config.time.start_date >= config.time.end_date
    error('开始时间必须早于结束时间');
end

fprintf('配置验证通过。\n');

end

