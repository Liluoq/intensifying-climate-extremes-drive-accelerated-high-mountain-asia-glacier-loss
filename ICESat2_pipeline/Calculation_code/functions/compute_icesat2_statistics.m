function results = compute_icesat2_statistics(txt_file, config)
% COMPUTE_ICESAT2_STATISTICS ICESat-2数据统计分析主函数
%
% 功能：
%   1. 月度、季节、年度时间序列分析
%   2. 空间网格统计
%   3. 高程带分析
%   4. 面积加权计算
%
% 输入:
%   txt_file - 增强数据文件路径
%   config - 配置结构体
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
%                  在原有13列基础上新增：
%                  列14: NASADEM参考高程
%                  列15: 坡度
%                  列16: 坡向  
%                  列17: 观测高程与参考DEM差值
%                  列18: 穿透深度校正值
%                  列19: 校正后的最终高程变化
% 输出:
%   results - 包含所有分析结果的结构体

fprintf('开始统计分析计算...\n');

%% 加载数据和空间参考
fprintf('加载数据和空间参考信息...\n');
try
    data = load(txt_file);
    
    % 加载空间边界数据
    glacier_grid = shaperead(config.files.grid10);
    rgi_regions = shaperead(config.files.gtn_regions);
    hma22_regions = shaperead(config.files.hma22_regions);
    hma4_regions = shaperead(config.files.hma_four_parts);
    hma6_regions = shaperead(config.files.hma_seasonal_pattern);
    tp_boundary = shaperead(config.files.tp_boundary);
    
    fprintf('✓ 数据加载完成，共 %d 个观测点\n', size(data, 1));
    
catch ME
    error('数据加载失败: %s', ME.message);
end
config.spatial.num_regions = 39;  % HMA + 15个GTN + 22个HMA22 + TP 

%% 数据预处理和质量控制,剔除dh>200m的点并计算HMA的分位数
fprintf('执行数据预处理和质量控制...\n');
[data, quality_stats] = apply_quality_control(data, config);

%% 分配区域编号（如果需要），暂时不需要
%fprintf('更新区域编号...\n');
%data = update_region_assignments(data, hma4_regions, hma6_regions);

%% 时间序列分析
fprintf('执行时间序列分析...\n');

% 月度分析，没有面积加权
if config.temporal.monthly_analysis
    results.monthly = compute_monthly_statistics(data, config);
    fprintf('✓ 月度分析完成\n');
end

% 年度分析
if config.temporal.annual_analysis
    results.annual = compute_annual_statistics(data, config);
    fprintf('✓ 年度分析完成\n');
end

% 季节分析,3个月滑动窗口
if config.temporal.seasonal_analysis
    results.seasonal = compute_seasonal_statistics(data, config);
    fprintf('✓ 季节分析完成\n');
end

%% 空间统计分析，不执行面积加权
fprintf('执行空间统计分析...\n');

% 网格级统计
results.grid = compute_grid_statistics(data, glacier_grid, config);
fprintf('✓ 网格统计完成\n');

% 区域统计，只统计各个区域的分位数
results.regions = compute_regional_statistics(data, rgi_regions, hma22_regions, tp_boundary, config);
fprintf('✓ 区域统计完成\n');

%% 高程带分析
fprintf('执行高程带分析...\n');
results.elevation_bands = compute_elevation_band_analysis(data, glacier_grid, config);
fprintf('✓ 高程带分析完成\n');

%% 面积加权计算(3个月季节平均)
fprintf('执行面积加权计算，输出季节结果...\n');
results.area_weighted_season = compute_area_weighted_analysis(results, glacier_grid, config);
fprintf('✓ 面积加权计算完成\n');
%% 面积加权计算(多年)
fprintf('执行面积加权计算，输出年尺度结果...\n');
results.area_weighted_annual = compute_area_weighted_analysis_year(results, glacier_grid, config);
fprintf('✓ 面积加权计算完成\n');

%% 添加处理元数据
results.metadata.processing_date = datetime('now');
results.metadata.data_points = size(data, 1);
results.metadata.quality_stats = quality_stats;
results.metadata.config = config;

fprintf('统计分析计算完成。\n');

end

function [data, quality_stats] = apply_quality_control(data, config)
% APPLY_QUALITY_CONTROL 应用质量控制过滤器
%
% 输入:
%   data - 原始数据
%   config - 配置参数
%
% 输出:
%   data - 质量控制后的数据
%   quality_stats - 质量统计信息

original_count = size(data, 1);

% 记录过滤步骤
quality_stats.original_points = original_count;

% 1. 排除绝对高程变化过大的点
to_delete = abs(data(:, 19)) > config.quality.max_height_change;%200 m
data(to_delete, :) = [];
quality_stats.after_height_filter = size(data, 1);
quality_stats.height_filtered = sum(to_delete);

% 2. 计算分位数界限
dif = sort(data(:, 19));
quality_stats.percentile_bounds = [prctile(dif, config.quality.outlier_percentile), ...
                                  prctile(dif, 100 - config.quality.outlier_percentile)];

fprintf('✓ 质量控制完成：%d → %d 个点 (移除 %d 个异常值)\n', ...
    original_count, size(data, 1), original_count - size(data, 1));

end

function data = update_region_assignments(data, hma4_regions, hma6_regions)
% UPDATE_REGION_ASSIGNMENTS 更新区域编号分配
%
% 输入:
%   data - 数据矩阵
%   hma4_regions - HMA四分区
%   hma6_regions - HMA六分区季节模式
%
% 输出:
%   data - 更新区域编号后的数据

% 为HMA四分区分配编号
for k = 1:length(hma4_regions)
    [in, ~] = inpolygon(data(:,1), data(:,2), hma4_regions(k).X, hma4_regions(k).Y);
    data(in==1, 20) = k;
end

% 为HMA六分区分配编号
for k = 1:length(hma6_regions)
    [in, ~] = inpolygon(data(:,1), data(:,2), hma6_regions(k).X, hma6_regions(k).Y);
    data(in==1, 21) = k;
end

end

function monthly_results = compute_monthly_statistics(data, config)
% COMPUTE_MONTHLY_STATISTICS 计算月度统计
%
% 输入:
%   data - 数据矩阵
%   config - 配置参数
%
% 输出:
%   monthly_results - 月度统计结果

% 生成月度时间窗口
start_date = config.time.start_date;
end_date = config.time.end_date;
ref_date = config.time.reference_date;

t = start_date:calmonths(1):end_date;
t2 = t + calmonths(1);

% 转换为相对天数
time_windows = [daysact(ref_date, t)', daysact(ref_date, t2)'];

% 初始化结果矩阵
num_periods = length(t);
num_regions = config.spatial.num_regions;  % HMA + 15个GTN + 22个HMA22 + TP 
monthly_results.values = zeros(num_periods, num_regions);
monthly_results.counts = zeros(num_periods, num_regions);
monthly_results.errors = zeros(num_periods, num_regions);
monthly_results.time_windows = time_windows;
monthly_results.mid_dates = ref_date + caldays(fix((time_windows(:, 1)+time_windows(:, 2))/2));

% 逐月计算统计
for i = 1:num_periods
    monthly_results = compute_period_statistics(data, time_windows(i,:), ...
        monthly_results, i, config);
end

monthly_results.description = '月度高程变化统计结果';

end

function annual_results = compute_annual_statistics(data, config)
% COMPUTE_ANNUAL_STATISTICS 计算年度统计
%
% 输入:
%   data - 数据矩阵
%   config - 配置参数
%
% 输出:
%   annual_results - 年度统计结果

ref_date = config.time.reference_date;

% 定义年度时间窗口
years = [2019:year(config.time.end_date);
        2020:year(config.time.end_date)+1];

% 两年期窗口
two_year_windows = [2019:year(config.time.end_date)-1;
                    2021:year(config.time.end_date)+1];
% 三年期窗口
three_year_windows = [2019:year(config.time.end_date)-2;
                    2022:year(config.time.end_date)+1];
whole_window = [2018;year(config.time.end_date)+1];

% 组合所有年度窗口
all_windows = [years, two_year_windows, three_year_windows, whole_window];

% 转换为相对天数
time_windows = zeros(size(all_windows, 2), 2);
for i = 1:size(all_windows, 2)
    time_windows(i, 1) = daysact(ref_date, datetime(all_windows(1, i), 1, 1));
    time_windows(i, 2) = daysact(ref_date, datetime(all_windows(2, i), 1, 1));
end

% 初始化结果
num_regions = config.spatial.num_regions;  % HMA + 15个GTN + 22个HMA22 + TP 
annual_results.values = zeros(size(time_windows, 1), num_regions);
annual_results.counts = zeros(size(time_windows, 1), num_regions);
annual_results.errors = zeros(size(time_windows, 1), num_regions);
annual_results.time_windows = time_windows;
annual_results.year_labels = all_windows';
annual_results.mid_dates = ref_date + caldays(fix((time_windows(:, 1)+time_windows(:, 2))/2));

% 计算年度统计
for i = 1:size(time_windows, 1)
    annual_results = compute_period_statistics(data, time_windows(i,:), ...
        annual_results, i, config);
end

annual_results.description = '年度和多年期高程变化统计结果';

end

function period_results = compute_period_statistics(data, time_window, period_results, period_idx, config)
% COMPUTE_PERIOD_STATISTICS 计算特定时间段的统计
%
% 输入:
%   data - 数据矩阵
%   time_window - 时间窗口 [开始天数, 结束天数]
%   period_results - 结果结构
%   period_idx - 时间段索引
%   config - 配置参数
%
% 输出:
%   period_results - 更新后的结果结构

% 提取时间段内的数据
time_mask = data(:,3) >= time_window(1) & data(:,3) < time_window(2);
period_data = data(time_mask, :);

if isempty(period_data)
    return;
end

col_idx = 1;

% HMA整体统计
[value, count, error] = compute_robust_statistics(period_data, [], config);
period_results.values(period_idx, col_idx) = value;
period_results.counts(period_idx, col_idx) = count;
period_results.errors(period_idx, col_idx) = error;
col_idx = col_idx + 1;

% RGI区域统计（15个区域）
for j = 1:15
    region_mask = period_data(:, 9) == j;
    [value, count, error] = compute_robust_statistics(period_data, region_mask, config);
    period_results.values(period_idx, col_idx) = value;
    period_results.counts(period_idx, col_idx) = count;
    period_results.errors(period_idx, col_idx) = error;
    col_idx = col_idx + 1;
end

% HMA22区域统计（22个区域）
for j = 1:22
    region_mask = period_data(:, 8) == j;
    [value, count, error] = compute_robust_statistics(period_data, region_mask, config);
    period_results.values(period_idx, col_idx) = value;
    period_results.counts(period_idx, col_idx) = count;
    period_results.errors(period_idx, col_idx) = error;
    col_idx = col_idx + 1;
end

% 青藏高原统计
tp_mask = period_data(:, 12) == 1;
[value, count, error] = compute_robust_statistics(period_data, tp_mask, config);
period_results.values(period_idx, col_idx) = value;
period_results.counts(period_idx, col_idx) = count;
period_results.errors(period_idx, col_idx) = error;
col_idx = col_idx + 1;

end

function [mean_value, count, std_error] = compute_robust_statistics(data, region_mask, config)
% COMPUTE_ROBUST_STATISTICS 计算鲁棒统计量，没有面积加权
%
% 输入:
%   data - 数据矩阵
%   region_mask - 区域掩膜（空则使用全部数据）
%   config - 配置参数
%
% 输出:
%   mean_value - 平均值
%   count - 有效数据点数
%   std_error - 标准误差

if isempty(region_mask)
    subset_data = data;
else
    subset_data = data(region_mask, :);
end

if isempty(subset_data)
    mean_value = NaN;
    count = 0;
    std_error = NaN;
    return;
end

% 使用第19列作为高程变化值
dh = subset_data(:, 19);

% 质量过滤
quality_mask = subset_data(:, 15) <= config.quality.max_slope & ...
               subset_data(:, 7) == 0 & ...
               abs(dh) < config.quality.height_change_filter1;

if ~any(quality_mask)
    mean_value = NaN;
    count = 0;
    std_error = NaN;
    return;
end

dh_filtered = dh(quality_mask);

% 计算中位数作为参考
median_dh = median(dh_filtered);

% 进一步过滤异常值
final_mask = abs(dh_filtered - median_dh) < config.quality.height_change_filter2;
dh_final = dh_filtered(final_mask);

% 计算统计量
if length(dh_final) > 100
    mean_value = mean(dh_final);
    count = length(dh_final);
    std_error = std(dh_final) / sqrt(count);
else
    mean_value = NaN;
    count = 0;
    std_error = NaN;
end

end

function seasonal_results = compute_seasonal_statistics(data, config)
% COMPUTE_SEASONAL_STATISTICS 计算季节统计
%
% 输入:
%   data - 数据矩阵
%   config - 配置参数
%
% 输出:
%   seasonal_results - 季节统计结果

ref_date = config.time.reference_date;
start_date = config.time.start_date;
end_date = config.time.end_date;
cycle_months = config.temporal.cycle_months;

% 生成3个月滑动窗口
t = start_date:calmonths(1):(end_date - calmonths(3));
t2 = t + calmonths(cycle_months);

time_windows = [daysact(ref_date, t)', daysact(ref_date, t2)'];

% 初始化结果
num_regions = config.spatial.num_regions;
seasonal_results.values = zeros(length(t), num_regions);
seasonal_results.counts = zeros(length(t), num_regions);
seasonal_results.errors = zeros(length(t), num_regions);
seasonal_results.time_windows = time_windows;
seasonal_results.mid_dates = ref_date + caldays(fix((time_windows(:, 1)+time_windows(:, 2))/2));

% 计算季节统计
for i = 1:length(t)
    seasonal_results = compute_period_statistics(data, time_windows(i,:), ...
        seasonal_results, i, config);
end

seasonal_results.description = '3个月滑动窗口季节统计结果';

end

function grid_results = compute_grid_statistics(data, glacier_grid, config)
% COMPUTE_GRID_STATISTICS 计算网格级统计
%
% 输入:
%   data - 数据矩阵
%   glacier_grid - 网格边界数据
%   config - 配置参数
%
% 输出:
%   grid_results - 网格统计结果

% 定义时间窗口（年度）
ref_date = config.time.reference_date;
years = [2019:year(config.time.end_date);
        2020:year(config.time.end_date)+1];
% 两年期窗口
two_year_windows = [2019:year(config.time.end_date)-1;
                    2021:year(config.time.end_date)+1];
% 三年期窗口
three_year_windows = [2019:year(config.time.end_date)-2;
                    2022:year(config.time.end_date)+1];
whole_window = [2018;year(config.time.end_date)+1];

% 组合所有年度窗口
all_windows = [years, two_year_windows, three_year_windows, whole_window];

% 转换为相对天数
year_windows = zeros(size(all_windows, 2), 2);
for i = 1:size(all_windows, 2)
    year_windows(i, 1) = daysact(ref_date, datetime(all_windows(1, i), 1, 1));
    year_windows(i, 2) = daysact(ref_date, datetime(all_windows(2, i), 1, 1));
end

num_grids = length(glacier_grid);
grid_results.values = zeros(num_grids, size(year_windows, 1));  % 多个时间窗口的结果
grid_results.metadata = zeros(num_grids, 4); % 网格元数据
grid_results.time_windows = year_windows;
grid_results.year_labels = all_windows';
grid_results.mid_dates = ref_date + caldays(fix((year_windows(:, 1)+year_windows(:, 2))/2));

for i = 1:num_grids
    % 网格元数据
    grid_results.metadata(i, 1) = glacier_grid(i).id;
    grid_results.metadata(i, 2) = glacier_grid(i).x;
    grid_results.metadata(i, 3) = glacier_grid(i).y;
    grid_results.metadata(i, 4) = glacier_grid(i).Area;%冰川面积
    
    % 网格内数据
    grid_mask = data(:, 11) == i;
    grid_data = data(grid_mask, :);
    
    if isempty(grid_data)
        continue;
    end
    
    % 计算各年度统计
    for k = 1:size(year_windows, 1)
        time_mask = grid_data(:,3) >= year_windows(k,1) & grid_data(:,3) < year_windows(k,2);
        
        if any(time_mask)
            subset_data = grid_data(time_mask, :);
            dh = subset_data(:, 19);
            
            % 质量过滤
            quality_mask = subset_data(:, 15) <= config.quality.max_slope & ...
                          subset_data(:, 7) == 0 & ...
                          abs(dh) < config.quality.height_change_filter1;
            
            if any(quality_mask)
                dh_filtered = dh(quality_mask);
                median_dh = median(dh_filtered);
                final_mask = abs(dh_filtered - median_dh) < config.quality.height_change_filter2;
                dh_final = dh_filtered(final_mask);
                
                if ~isempty(dh_final)
                    grid_results.values(i, k) = mean(dh_final);
                end
            end
        end
    end
end

grid_results.description = '1度网格级年度统计结果';

end

function regional_results = compute_regional_statistics(data, rgi_regions, hma22_regions, tp_boundary, config)
% COMPUTE_REGIONAL_STATISTICS 计算区域统计
%
% 输入:
%   data - 数据矩阵
%   rgi_regions, hma22_regions, tp_boundary - 区域边界
%   config - 配置参数
%
% 输出:
%   regional_results - 区域统计结果

% 计算各区域的统计界限
regional_results.bounds = calculate_regional_bounds(data, config);
regional_results.bounds_lable = ['1 percentile', '99 percentile', 'median'];

regional_results.description = '区域长期统计结果';

end

function bounds = calculate_regional_bounds(data, config)
% CALCULATE_REGIONAL_BOUNDS 计算各区域的统计界限
%
% 输入:
%   data - 数据矩阵
%   config - 配置参数
%
% 输出:
%   bounds - 区域界限矩阵

var = config.quality.outlier_percentile;% 1
bounds = zeros(config.spatial.num_regions, 3);  % HMA + 15个RGI + 22个HMA22 + TP

% 整体界限
dif = sort(data(:, 19));
bounds(1, 1) = prctile(dif, var);
bounds(1, 2) = prctile(dif, 100-var);
bounds(1, 3) = prctile(dif, 50);

% GTN区域界限
for i = 1:15
    region_data = data(data(:, 9) == i, :);
    if ~isempty(region_data)
        dif0 = sort(region_data(:, 19));
        bounds(i+1, 1) = prctile(dif0, var);
        bounds(i+1, 2) = prctile(dif0, 100-var);
        bounds(i+1, 3) = prctile(dif0, 50);
    end
end

% HMA22区域界限
for i = 1:22
    region_data = data(data(:, 8) == i, :);
    if ~isempty(region_data)
        dif0 = sort(region_data(:, 19));
        bounds(i+16, 1) = prctile(dif0, var);
        bounds(i+16, 2) = prctile(dif0, 100-var);
        bounds(i+16, 3) = prctile(dif0, 50);
    end
end

end

function [slope, intercept, p_value, uncertainity] = calculate_regional_trends(x, y)
% CALCULATE_REGIONAL_TRENDS 计算各区域趋势
%
% 输入:
%   data - 子区域，数据矩阵
%   config - 配置参数
%
% 输出:
%   [slope, intercept, p_value, uncertainity] 鲁棒线性拟合结果

% 为简化起见，这里返回基本的区域平均值
% 在实际应用中可以扩展为详细的趋势分析
[slope, intercept, p_value, uncertainity] = deal(nan);

x_final = x;
dh_final = y;

if size(dh_final, 1) < 10
    return
end

%计算趋势，x是dh(:, 3), y是dh(:, 19)，鲁棒最小二乘拟合，输出p值和slope的95%不确定区间
try
    % 使用 fitlm 进行鲁棒拟合
    mdl = fitlm(x_final, dh_final, 'RobustOpts', 'on');
    
    % 斜率和截距
    intercept = mdl.Coefficients.Estimate(1);
    slope     = mdl.Coefficients.Estimate(2);
    
    % 斜率的显著性 p 值
    p_value   = mdl.Coefficients.pValue(2);

    % 斜率的 95% 置信区间
    CI        = coefCI(mdl, 0.05);
    uncertainity = (CI(2,2) - CI(2,1)) / 2; % 半宽
catch
    return
end
end

function elevation_results = compute_elevation_band_analysis(data, glacier_grid, config)
% COMPUTE_ELEVATION_BAND_ANALYSIS 计算高程带分析
%
% 输入:
%   data - 数据矩阵
%   config - 配置参数
%
% 输出:
%   elevation_results - 高程带分析结果

fprintf('  开始高程带分析...\n');

% 定义高程带
elevation_range = config.spatial.elevation_range; %[2250, 8750]
band_width = config.spatial.elevation_band_width(2);  % 使用100m带宽
elevation_centers = elevation_range(1):band_width:elevation_range(2);
num_bands = length(elevation_centers);

% 时间窗口设置
ref_date = config.time.reference_date;
start_date = config.time.start_date;
end_date = config.time.end_date;
% 是0的都不是好东西
% 3个月滑动窗口
t = start_date:calmonths(1):(end_date - calmonths(3));
t2 = t + calmonths(3);
time_windows = [daysact(ref_date, t)', daysact(ref_date, t2)'];

% 初始化结果矩阵
num_periods = length(t);
num_regions = config.spatial.num_regions + 496;  % HMA + RGI + HMA22 + TP + grid496
elevation_results.values = zeros(num_bands, num_regions, num_periods);
elevation_results.counts = zeros(num_bands, num_regions, num_periods);
elevation_results.errors = zeros(num_bands, num_regions, num_periods);

% 逐时间段分析
for m = 1:num_periods
    time_mask = data(:,3) >= time_windows(m,1) & data(:,3) < time_windows(m,2);
    period_data = data(time_mask, :);
    
    if isempty(period_data)
        continue;
    end
    
    % 质量过滤
    quality_mask = period_data(:, 15) <= config.quality.max_slope & ...
                   period_data(:, 7) == 0 & ...
                   abs(period_data(:, 19)) < 120;
    
    period_data = period_data(quality_mask, :);
    
    if isempty(period_data)
        continue;
    end
    
    % 逐高程带分析
    for k = 1:num_bands
        elevation_mask = period_data(:, 14) >= (elevation_centers(k) - band_width/2) & ...
                        period_data(:, 14) < (elevation_centers(k) + band_width/2);
        
        if any(elevation_mask)
            band_data = period_data(elevation_mask, :);
            dh = band_data(:, 19);
            
            % HMA整体
            elevation_results.values(k, 1, m) = mean(dh, 'omitnan');
            elevation_results.counts(k, 1, m) = length(dh);
            elevation_results.errors(k, 1, m) = std(dh, 'omitnan') / sqrt(length(dh));
            
            % RGI区域
            for j = 1:15
                region_mask = band_data(:, 9) == j;
                if any(region_mask)
                    region_dh = dh(region_mask);
                    elevation_results.values(k, j+1, m) = mean(region_dh, 'omitnan');
                    elevation_results.counts(k, j+1, m) = length(region_dh);
                    elevation_results.errors(k, j+1, m) = std(region_dh, 'omitnan') / sqrt(length(region_dh));
                end
            end
            
            % HMA22区域
            for j = 1:22
                region_mask = band_data(:, 8) == j;
                if any(region_mask)
                    region_dh = dh(region_mask);
                    elevation_results.values(k, j+16, m) = mean(region_dh, 'omitnan');
                    elevation_results.counts(k, j+16, m) = length(region_dh);
                    elevation_results.errors(k, j+16, m) = std(region_dh, 'omitnan') / sqrt(length(region_dh));
                end
            end
            
            % 青藏高原
            tp_mask = band_data(:, 12) == 1;
            if any(tp_mask)
                tp_dh = dh(tp_mask);
                elevation_results.values(k, 39, m) = mean(tp_dh, 'omitnan');
                elevation_results.counts(k, 39, m) = length(tp_dh);
                elevation_results.errors(k, 39, m) = std(tp_dh, 'omitnan') / sqrt(length(tp_dh));
            end

            % grid
            for j = 1:length(glacier_grid)
                region_mask = band_data(:, 11)==j;
                if any(region_mask)
                    region_dh = dh(region_mask);
                    elevation_results.values(k, j+39, m) = mean(region_dh, 'omitnan');
                    elevation_results.counts(k, j+39, m) = length(region_dh);
                    elevation_results.errors(k, j+39, m) = std(region_dh, 'omitnan') / sqrt(length(region_dh));
                end
            end
        end
    end
end

elevation_results.elevation_centers = elevation_centers';
elevation_results.time_windows = time_windows;
elevation_results.mid_dates = ref_date + caldays(fix((time_windows(:, 1)+time_windows(:, 2))/2));
elevation_results.description = '高程带时间序列分析结果';
% 单年
years = [2019:year(config.time.end_date);
        2020:year(config.time.end_date)+1];
% 两年期窗口
two_year_windows = [2019:year(config.time.end_date)-1;
                    2021:year(config.time.end_date)+1];
% 三年期窗口
three_year_windows = [2019:year(config.time.end_date)-2;
                    2022:year(config.time.end_date)+1];
whole_window = [2018;year(config.time.end_date)+1];

% 组合所有年度窗口
all_windows = [years, two_year_windows, three_year_windows, whole_window];

% 转换为相对天数
year_windows = zeros(size(all_windows, 2), 2);
for i = 1:size(all_windows, 2)
    year_windows(i, 1) = daysact(ref_date, datetime(all_windows(1, i), 1, 1));
    year_windows(i, 2) = daysact(ref_date, datetime(all_windows(2, i), 1, 1));
end
num_periods = size(year_windows, 1);
num_regions = config.spatial.num_regions + 496;  % HMA + RGI + HMA22 + TP + grid496
elevation_results.values_year = zeros(num_bands, num_regions, num_periods);
elevation_results.counts_year = zeros(num_bands, num_regions, num_periods);
elevation_results.errors_year = zeros(num_bands, num_regions, num_periods);
elevation_results.year_label = all_windows';
elevation_results.year_windows = year_windows;

% 逐年段分析
for m = 1:num_periods
    time_mask = data(:,3) >= year_windows(m,1) & data(:,3) < year_windows(m,2);
    period_data = data(time_mask, :);
    
    if isempty(period_data)
        continue;
    end
    
    % 质量过滤
    quality_mask = period_data(:, 15) <= config.quality.max_slope & ...
                   period_data(:, 7) == 0 & ...
                   abs(period_data(:, 19)) < 120;
    
    period_data = period_data(quality_mask, :);
    
    if isempty(period_data)
        continue;
    end
    
    % 逐高程带分析
    for k = 1:num_bands
        elevation_mask = period_data(:, 14) >= (elevation_centers(k) - band_width/2) & ...
                        period_data(:, 14) < (elevation_centers(k) + band_width/2);
        
        if any(elevation_mask)
            band_data = period_data(elevation_mask, :);
            dh = band_data(:, 19);
            
            % HMA整体
            elevation_results.values_year(k, 1, m) = mean(dh, 'omitnan');
            elevation_results.counts_year(k, 1, m) = length(dh);
            elevation_results.errors_year(k, 1, m) = std(dh, 'omitnan') / sqrt(length(dh));
            
            % RGI区域
            for j = 1:15
                region_mask = band_data(:, 9) == j;
                if any(region_mask)
                    region_dh = dh(region_mask);
                    elevation_results.values_year(k, j+1, m) = mean(region_dh, 'omitnan');
                    elevation_results.counts_year(k, j+1, m) = length(region_dh);
                    elevation_results.errors_year(k, j+1, m) = std(region_dh, 'omitnan') / sqrt(length(region_dh));
                end
            end
            
            % HMA22区域
            for j = 1:22
                region_mask = band_data(:, 8) == j;
                if any(region_mask)
                    region_dh = dh(region_mask);
                    elevation_results.values_year(k, j+16, m) = mean(region_dh, 'omitnan');
                    elevation_results.counts_year(k, j+16, m) = length(region_dh);
                    elevation_results.errors_year(k, j+16, m) = std(region_dh, 'omitnan') / sqrt(length(region_dh));
                end
            end
            
            % 青藏高原
            tp_mask = band_data(:, 12) == 1;
            if any(tp_mask)
                tp_dh = dh(tp_mask);
                elevation_results.values_year(k, 39, m) = mean(tp_dh, 'omitnan');
                elevation_results.counts_year(k, 39, m) = length(tp_dh);
                elevation_results.errors_year(k, 39, m) = std(tp_dh, 'omitnan') / sqrt(length(tp_dh));
            end

            % grid
            for j = 1:length(glacier_grid)
                region_mask = band_data(:, 11)==j;
                if any(region_mask)
                    region_dh = dh(region_mask);
                    elevation_results.values_year(k, j+39, m) = mean(region_dh, 'omitnan');
                    elevation_results.counts_year(k, j+39, m) = length(region_dh);
                    elevation_results.errors_year(k, j+39, m) = std(region_dh, 'omitnan') / sqrt(length(region_dh));
                end
            end
        end
    end
end

end

function area_weighted_results = compute_area_weighted_analysis(results, glacier_grid, config)
% COMPUTE_AREA_WEIGHTED_ANALYSIS 计算面积加权分析
%
% 输入:
%   results - 之前的分析结果
%   config - 配置参数
%
% 输出:
%   area_weighted_results - 面积加权结果

fprintf('  计算面积加权平均值...\n');

try
    % 加载高程带面积数据
    area_data_rgi = readmatrix(config.files.elevation_area, 'Sheet', 'HMA_RGI_100', 'Range', 'A1:Q66');%ele+hma+15 gtn
    area_data_hma22 = readmatrix(config.files.elevation_area, 'Sheet', 'HMA22_100', 'Range', 'A1:Y66');%ele+hma+tp + 22
    area_data_grid = readmatrix(config.files.elevation_area, 'Sheet', 'grid10_100_496', 'Range', 'A2:SC67');%ele+496 grid
    
    % 合并面积数据
    combined_area = [area_data_rgi(:, 2:end), area_data_hma22(:, 4:25), area_data_hma22(:, 3), area_data_grid(:, 2:end)];%hma+15 gtn+22 hma+tp+496 grid
    
    % 检查高程带分析结果是否存在
    if ~isfield(results, 'elevation_bands')
        warning('高程带分析结果不存在，跳过面积加权计算');
        area_weighted_results = [];
        return;
    end
    
    elevation_values = results.elevation_bands.values;% num_bands*num_regions*n_time
    num_periods = size(elevation_values, 3);
    num_regions = config.spatial.num_regions+496;  % HMA + RGI + HMA22 + TP + grid496
    
    % 初始化面积加权结果
    area_weighted_results.values = zeros(num_regions, num_periods);
    area_weighted_results.total_counts = zeros(num_regions, num_periods);
    area_weighted_results.weighted_errors = zeros(num_regions, num_periods);
    
    % 逐时间段计算面积加权平均
    for j = 1:num_periods
        ice_values = elevation_values(:, :, j);
        count_values = results.elevation_bands.counts(:, :, j);
        error_values = results.elevation_bands.errors(:, :, j);
        
        for i = 1:num_regions
            area_column = combined_area(:, i);
            ice_column = ice_values(:, i);
            count_column = count_values(:, i);
            error_column = error_values(:, i);
            
            % 构建面积加权数据
            area_weight_data = [results.elevation_bands.elevation_centers, area_column, ...
                               area_column .* ice_column, count_column, ...
                               ice_column, error_column];
            
            % 移除NaN值
            valid_mask = ~isnan(area_weight_data(:, 3));
            area_weight_data = area_weight_data(valid_mask, :);
            
            % 质量过滤
            quality_mask = area_weight_data(:, 4) >= 10 & abs(area_weight_data(:, 5)) <= 50 & abs(area_weight_data(:, 5)) > 0;
            area_weight_data = area_weight_data(quality_mask, :);
            
            % 异常值过滤
            % outlier_mask = abs(area_weight_data(:, 5)) <= 80;
            % area_weight_data = area_weight_data(outlier_mask, :);
            
            if ~isempty(area_weight_data) && sum(area_weight_data(:, 2)) > 0
                % 面积加权平均
                area_weighted_results.values(i, j) = sum(area_weight_data(:, 3)) / sum(area_weight_data(:, 2));
                area_weighted_results.total_counts(i, j) = sum(area_weight_data(:, 4));
                area_weighted_results.weighted_errors(i, j) = mean(area_weight_data(:, 6)) / sqrt(length(area_weight_data(:, 6)));
            end
        end
    end
    area_weighted_results.values = area_weighted_results.values';
    area_weighted_results.total_counts = area_weighted_results.total_counts';
    area_weighted_results.weighted_errors = area_weighted_results.weighted_errors';%time*num_regions

    area_weighted_results.mid_dates = results.elevation_bands.mid_dates;
    area_weighted_results.description = '面积加权区域平均结果';
    
    ref_date = config.time.reference_date;
    days = daysact(ref_date, area_weighted_results.mid_dates);%列向量
    area_weighted_results.mid_days = days;

    %计算趋势
    area_weighted_results.trends = nan(num_regions,1);
    area_weighted_results.intecepts = nan(num_regions,1);
    area_weighted_results.pvalues = nan(num_regions,1);
    area_weighted_results.trend_uncerts = nan(num_regions,1);
    for i = 1:num_regions
        area_weight_dh = area_weighted_results.values(:, i);%列向量
        area_weight_count = area_weighted_results.total_counts(:, i);

        valid_day_mask = (area_weight_dh ~= 0) & (~isnan(area_weight_dh)) & (area_weight_count > 0);

        valid_days = days(valid_day_mask);
        valid_dh = area_weight_dh(valid_day_mask);

        if ~isempty(valid_dh)
            [slope, intercept, p_value, uncertainity] = calculate_regional_trends(valid_days, valid_dh);
            area_weighted_results.trends(i, 1) = slope; %m/day
            area_weighted_results.intecepts(i, 1) = intercept;
            area_weighted_results.pvalues(i, 1) = p_value;
            area_weighted_results.trend_uncerts(i, 1) = uncertainity;%m/day
        end
    end
    
catch ME
    fprintf('❌ 面积加权计算失败: %s\n', ME.message);
    area_weighted_results = [];
end

end

function area_weighted_results = compute_area_weighted_analysis_year(results, glacier_grid, config)
% COMPUTE_AREA_WEIGHTED_ANALYSIS 计算面积加权分析
%
% 输入:
%   results - 之前的分析结果
%   config - 配置参数
%
% 输出:
%   area_weighted_results - 面积加权结果

fprintf('  计算面积加权平均值...\n');

% 加载高程带面积数据
area_data_rgi = readmatrix(config.files.elevation_area, 'Sheet', 'HMA_RGI_100', 'Range', 'A1:Q66');%ele+hma+15 gtn
area_data_hma22 = readmatrix(config.files.elevation_area, 'Sheet', 'HMA22_100', 'Range', 'A1:Y66');%ele+hma+tp + 22
area_data_grid = readmatrix(config.files.elevation_area, 'Sheet', 'grid10_100_496', 'Range', 'A2:SC67');%ele+496 grid

% 合并面积数据
combined_area = [area_data_rgi(:, 2:end), area_data_hma22(:, 4:25), area_data_hma22(:, 3), area_data_grid(:, 2:end)];%hma+15 gtn+22 hma+tp+496 grid

% 检查高程带分析结果是否存在
if ~isfield(results, 'elevation_bands')
    warning('高程带分析结果不存在，跳过面积加权计算');
    area_weighted_results = [];
    return;
end

elevation_values = results.elevation_bands.values_year;% num_bands*num_regions*n_time
num_periods = size(elevation_values, 3);
num_regions = config.spatial.num_regions+496;  % HMA + RGI + HMA22 + TP + grid496

% 初始化面积加权结果
area_weighted_results.values = zeros(num_regions, num_periods);
area_weighted_results.total_counts = zeros(num_regions, num_periods);
area_weighted_results.weighted_errors = zeros(num_regions, num_periods);

% 逐时间段计算面积加权平均
for j = 1:num_periods
    ice_values = elevation_values(:, :, j);
    count_values = results.elevation_bands.counts_year(:, :, j);
    error_values = results.elevation_bands.errors_year(:, :, j);
    
    for i = 1:num_regions
        area_column = combined_area(:, i);
        ice_column = ice_values(:, i);
        count_column = count_values(:, i);
        error_column = error_values(:, i);
        
        % 构建面积加权数据
        area_weight_data = [results.elevation_bands.elevation_centers, area_column, ...
                           area_column .* ice_column, count_column, ...
                           ice_column, error_column];
        
        % 移除NaN值
        valid_mask = ~isnan(area_weight_data(:, 3));
        area_weight_data = area_weight_data(valid_mask, :);
        
        % 质量过滤
        quality_mask = area_weight_data(:, 4) >= 10 & abs(area_weight_data(:, 5)) <= 50 & abs(area_weight_data(:, 5)) > 0;
        area_weight_data = area_weight_data(quality_mask, :);
        
        % 异常值过滤
        % outlier_mask = abs(area_weight_data(:, 5)) <= 80;
        % area_weight_data = area_weight_data(outlier_mask, :);
        
        if ~isempty(area_weight_data) && sum(area_weight_data(:, 2)) > 0
            % 面积加权平均
            area_weighted_results.values(i, j) = sum(area_weight_data(:, 3)) / sum(area_weight_data(:, 2));
            area_weighted_results.total_counts(i, j) = sum(area_weight_data(:, 4));
            area_weighted_results.weighted_errors(i, j) = mean(area_weight_data(:, 6)) / sqrt(length(area_weight_data(:, 6)));
        end
    end
end
area_weighted_results.values = area_weighted_results.values';
area_weighted_results.total_counts = area_weighted_results.total_counts';
area_weighted_results.weighted_errors = area_weighted_results.weighted_errors';%time*num_regions

area_weighted_results.year_label = results.elevation_bands.year_label;
area_weighted_results.description = '年尺度面积加权区域平均结果';

ref_date = config.time.reference_date;
year_windows = results.elevation_bands.year_windows;
dates = ref_date + caldays(fix((year_windows(:, 1)+year_windows(:, 2))/2));%列向量
area_weighted_results.mid_dates = dates;

end