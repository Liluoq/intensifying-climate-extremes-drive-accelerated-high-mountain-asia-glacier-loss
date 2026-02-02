function results = compute_cryosat2_statistics(txt_file, config)
% COMPUTE_CRYOSAT2_STATISTICS CryoSat-2数据统计分析主函数
%
% 参考：
%   - revised_ICESat2_code/functions/compute_icesat2_statistics.m
%   - HMA_CryoSat2_COMPUTE_20230423m.m
%
% 数据列定义（16列）：
%   列1: lon, 列2: lat, 列3: elevation, 列4: days(1970-1-1), 列5: year, 列6: month
%   列7: geoid, 列8: uncertainty, 列9: HMA22, 列10: GTN/RGI, 列11: grid05
%   列12: grid10, 列13: TP, 列14: NASADEM(11×11), 列15: slope(7×7), 列16: dh
%
% 输入:
%   txt_file - 增强数据文件路径（16列）
%   config - 配置结构体
%
% 输出:
%   results - 包含所有分析结果的结构体

fprintf('开始CryoSat-2统计分析计算...\n');

%% 加载数据
fprintf('加载数据...\n');
try
    data = load(txt_file);
    fprintf('✓ 数据加载完成，共 %d 个观测点\n', size(data, 1));
    fprintf('  数据维度: %d × %d\n', size(data, 1), size(data, 2));
catch ME
    error('数据加载失败: %s', ME.message);
end

%% 数据预处理和质量控制
fprintf('执行质量控制...\n');
[data, quality_stats] = apply_quality_control(data, config);

%% 加载空间参考数据
fprintf('加载空间参考数据...\n');
try
    glacier_grid = shaperead(config.files.grid10);
    HMA22 = shaperead(config.files.hma22_regions);
    rgiregion = shaperead(config.files.gtn_regions);
    fprintf('✓ 空间参考加载完成\n');
catch ME
    warning('空间参考加载失败: %s', ME.message);
    glacier_grid = [];
    HMA22 = [];
    rgiregion = [];
end

%% 年度统计分析（单年，2010-2025, 整个HMA）
fprintf('执行年度统计分析...\n');
results.annual = compute_annual_statistics(data, config);
fprintf('✓ 年度分析完成\n');

%% 多年期统计（3年窗口，2010-2025滑动，HMA子区域，不面积加权）
fprintf('执行多年期统计（3年窗口）...\n');
results.multiyear = compute_multiyear_statistics(data, config);
fprintf('✓ 多年期分析完成\n');

%% 月度统计（2010年7月开始，不滑动）
fprintf('执行月度统计分析...\n');
results.monthly = compute_monthly_statistics(data, config);
fprintf('✓ 月度分析完成\n');

%% 网格统计
fprintf('执行网格统计...\n');
if ~isempty(glacier_grid)
    results.grid = compute_grid_statistics(data, glacier_grid, config);
    fprintf('✓ 网格统计完成\n');
end

%% 高程带分析（ele×区域*num_periods）
fprintf('执行高程带分析（季节+年尺度）...\n');
results.elevation_bands = compute_elevation_band_seasonal(data, HMA22, rgiregion, glacier_grid, config);
fprintf('✓ 季节/年高程带分析完成\n');

%% 面积加权计算
fprintf('执行面积加权计算...\n');
results.area_weighted = compute_area_weighted_analysis(results, config);
fprintf('✓ 面积加权完成\n');

%% 添加元数据
results.metadata.processing_date = datetime('now');
results.metadata.data_points = size(data, 1);
results.metadata.quality_stats = quality_stats;
results.metadata.config = config;

fprintf('✓ 统计分析计算完成\n');

end

function [data, quality_stats] = apply_quality_control(data, config)
% APPLY_QUALITY_CONTROL 应用质量控制

original_count = size(data, 1);
quality_stats.original_points = original_count;

% 排除异常值（列16是高程变化）
to_delete = abs(data(:, 16)) > config.quality.max_height_change;
data(to_delete, :) = [];
quality_stats.after_height_filter = size(data, 1);

% 百分位数界限
dif = sort(data(:, 16));
quality_stats.percentile_bounds = [prctile(dif, config.quality.outlier_percentile_low), ...
                                   prctile(dif, config.quality.outlier_percentile_high)];
quality_stats.percentile_descrp = '5和95分位数';
fprintf('✓ 质量控制：%d → %d 个点\n', original_count, size(data, 1));

end

function annual_results = compute_annual_statistics(data, config)
% COMPUTE_ANNUAL_STATISTICS 计算年度统计
%
% 使用列5(year)和列16(dh)

n_years = config.temporal.end_year - config.temporal.start_year + 1;
annual_results.values = zeros(n_years, 10);
annual_results.values(:,1) = linspace(config.temporal.start_year, config.temporal.end_year, n_years);

for i = 1:n_years
    year = annual_results.values(i, 1);
    
    % 筛选该年数据（使用列5: year，列15: slope，列16: dh）
    year_mask = (data(:,5) == year) & (data(:,15) <= config.quality.max_slope);
    
    if sum(year_mask) == 0
        continue;
    end
    
    year_dh = data(year_mask, 16);  % 列16是高程变化
    
    % 鲁棒统计：第一轮过滤
    valid1 = abs(year_dh) < 150;
    if sum(valid1) == 0
        continue;
    end
    
    % 第二轮：基于中位数过滤
    med = median(year_dh(valid1));
    valid = abs(year_dh - med) < 75;
    
    if sum(valid) > 0
        annual_results.values(i, 2) = mean(year_dh(valid), 'omitnan');
        annual_results.values(i, 3) = median(year_dh(valid));
        annual_results.values(i, 4) = std(year_dh(valid), 'omitnan');
        annual_results.values(i, 5) = std(year_dh(valid), 'omitnan') / sqrt(sum(valid));
        annual_results.values(i, 6) = sum(valid);
        annual_results.values(i, 7) = prctile(year_dh(valid), 25);
        annual_results.values(i, 8) = prctile(year_dh(valid), 75);
        annual_results.values(i, 9) = min(year_dh(valid));
        annual_results.values(i, 10) = max(year_dh(valid));
    end
    annual_results.descrp = 'HMA annual results';
    annual_results.col_name = {'year', 'mean', 'median', 'std', 'std error',...
                                'count', '25 percentile', '75 percentile', 'min', 'max'};
end

end

function spatial_results = compute_spatial_statistics(data, config)
% COMPUTE_SPATIAL_STATISTICS 计算空间统计
%
% 使用列9(HMA22)、列10(RGI)、列15(slope)、列16(dh)

% HMA整体统计
spatial_results.hma_overall = compute_region_stats(data, [], config);

% HMA22子区域统计（列9）
spatial_results.hma22 = zeros(22, 10);
for i = 1:22
    region_mask = (data(:,9) == i);  % 列9是HMA22
    spatial_results.hma22(i, :) = compute_region_stats(data, region_mask, config);
end

% RGI区域统计（列10）
spatial_results.rgi = zeros(15, 10);
for i = 1:15
    region_mask = (data(:,10) == i);  % 列10是RGI
    spatial_results.rgi(i, :) = compute_region_stats(data, region_mask, config);
end

end

function region_stats = compute_region_stats(data, mask, config)
% COMPUTE_REGION_STATS 计算单个区域统计
%
% 关键：列15是坡度，列16是高程变化

region_stats = zeros(1, 10);

if isempty(mask)
    mask = true(size(data, 1), 1);
end

% 质量控制（列15: slope, 列16: dh）
valid_mask = mask & (data(:,15) <= config.quality.max_slope) & ...
             (abs(data(:,16)) < 150);

if sum(valid_mask) == 0
    return;
end

dh = data(valid_mask, 16);  % 列16是高程变化

% 鲁棒统计：基于中位数过滤
med = median(dh);
valid = abs(dh - med) < 75;

if sum(valid) > 0
    region_stats(1) = mean(dh(valid), 'omitnan');
    region_stats(2) = median(dh(valid));
    region_stats(3) = std(dh(valid), 'omitnan');
    region_stats(4) = std(dh(valid), 'omitnan') / sqrt(sum(valid));
    region_stats(5) = sum(valid);
    region_stats(6) = prctile(dh(valid), 25);
    region_stats(7) = prctile(dh(valid), 75);
    region_stats(8) = min(dh(valid));
    region_stats(9) = max(dh(valid));
    region_stats(10) = sum(valid_mask);
end

end

function multiyear_results = compute_multiyear_statistics(data, config)
% COMPUTE_MULTIYEAR_STATISTICS 计算多年期统计（3年窗口）
%
% 参考：HMA_CryoSat2_COMPUTE_20230423m.m 第134-210行
% 使用列4(days)、列5(year)、列15(slope)、列16(dh)、列9-10(区域)

n_years = config.temporal.end_year - config.temporal.start_year + 1;
multiyear_results = struct();
multiyear_results.values = zeros(n_years, 39);  % HMA+15RGI+22HMA22+TP
multiyear_results.counts = zeros(n_years, 39);
multiyear_results.errors = zeros(n_years, 39);
multiyear_results.years = (config.temporal.start_year: config.temporal.end_year)';
multiyear_results.time_windows = zeros(n_years, 3);

for i = 1:n_years
    year_start = 2010 + i - 1;
    
    % 时间窗口（相对1970-01-01，因为原始代码用这个）
    start_day = datenum(year_start - 1, 1, 1) - datenum(1970, 1, 1);
    end_day = datenum(year_start + 2, 1, 1) - datenum(1970, 1, 1);
    
    multiyear_results.time_windows(i,:) = [year_start, start_day, end_day];
    
    % 时间窗口筛选
    time_mask = (data(:, 4) >= start_day) & (data(:, 4) < end_day) & ...
                (data(:,15) <= config.quality.max_slope);  % 坡度过滤
    
    col_idx = 1;
    
    % HMA整体
    [val, cnt, err] = compute_robust_stats(data, time_mask, 16); %自带dh阈值过滤
    multiyear_results.values(i, col_idx) = val;
    multiyear_results.counts(i, col_idx) = cnt;
    multiyear_results.errors(i, col_idx) = err;
    col_idx = col_idx + 1;
    
    % 15个RGI区域（列10）
    for j = 1:15
        region_mask = time_mask & (data(:,10) == j);
        [val, cnt, err] = compute_robust_stats(data, region_mask, 16);
        multiyear_results.values(i, col_idx) = val;
        multiyear_results.counts(i, col_idx) = cnt;
        multiyear_results.errors(i, col_idx) = err;
        col_idx = col_idx + 1;
    end
    
    % 22个HMA22区域（列9）
    for j = 1:22
        region_mask = time_mask & (data(:,9) == j);
        [val, cnt, err] = compute_robust_stats(data, region_mask, 16);
        multiyear_results.values(i, col_idx) = val;
        multiyear_results.counts(i, col_idx) = cnt;
        multiyear_results.errors(i, col_idx) = err;
        col_idx = col_idx + 1;
    end
    
    % TP（列13）
    tp_mask = time_mask & (data(:,13) == 1);
    [val, cnt, err] = compute_robust_stats(data, tp_mask, 16);
    multiyear_results.values(i, col_idx) = val;
    multiyear_results.counts(i, col_idx) = cnt;
    multiyear_results.errors(i, col_idx) = err;
end
multiyear_results.descrp = '三年滑动平均结果';
end

function monthly_results = compute_monthly_statistics(data, config)
% COMPUTE_MONTHLY_STATISTICS 计算月度统计
%
% 参考：HMA_CryoSat2_COMPUTE_20230423m.m 第214-285行
% 183个月，从2010年7月开始

t = [datetime(2010, 7, 1):calmonths(1):(config.time.end_date-calmonths(1))]';
t2 = t + calmonths(1);
time_windows = [daysact(datetime(1970,1,1), t), daysact(datetime(1970,1,1), t2)];

n_periods = size(time_windows, 1);
monthly_results = struct();
monthly_results.values = zeros(n_periods, 39);
monthly_results.counts = zeros(n_periods, 39);
monthly_results.errors = zeros(n_periods, 39);
monthly_results.time_windows = time_windows;
monthly_results.descrp = '月尺度结果(不平滑)';

for i = 1:n_periods
    start_day = time_windows(i, 1);
    end_day = time_windows(i, 2);
    
    time_mask = (data(:,4) >= start_day) & (data(:,4) < end_day) & ...
                (data(:,15) <= 40);
    
    col_idx = 1;
    
    % HMA整体
    [val, cnt, err] = compute_robust_stats(data, time_mask, 16);
    monthly_results.values(i, col_idx) = val;
    monthly_results.counts(i, col_idx) = cnt;
    monthly_results.errors(i, col_idx) = err;
    col_idx = col_idx + 1;
    
    % 15个RGI区域
    for j = 1:15
        region_mask = time_mask & (data(:,10) == j);
        [val, cnt, err] = compute_robust_stats(data, region_mask, 16);
        monthly_results.values(i, col_idx) = val;
        monthly_results.counts(i, col_idx) = cnt;
        monthly_results.errors(i, col_idx) = err;
        col_idx = col_idx + 1;
    end
    
    % 22个HMA22区域
    for j = 1:22
        region_mask = time_mask & (data(:,9) == j);
        [val, cnt, err] = compute_robust_stats(data, region_mask, 16);
        monthly_results.values(i, col_idx) = val;
        monthly_results.counts(i, col_idx) = cnt;
        monthly_results.errors(i, col_idx) = err;
        col_idx = col_idx + 1;
    end
    
    % TP（列13）
    tp_mask = time_mask & (data(:,13) == 1);
    [val, cnt, err] = compute_robust_stats(data, tp_mask, 16);
    monthly_results.values(i, col_idx) = val;
    monthly_results.counts(i, col_idx) = cnt;
    monthly_results.errors(i, col_idx) = err;
end

monthly_results.mid_date = t + caldays(floor(daysact(t, t2)/2));
end

function grid_results = compute_grid_statistics(data, glacier_grid, config)
% COMPUTE_GRID_STATISTICS 网格统计
%
% 参考：HMA_CryoSat2_COMPUTE_20230423m.m 第287-320行
% 使用列12(grid10)、列15(slope)、列16(dh)
num_periods = config.temporal.end_year - config.temporal.start_year + 1;

num_grids = length(glacier_grid);
grid_results = struct();
grid_results.metadata = zeros(num_grids, 4);
grid_results.values_3yr = zeros(num_periods, num_grids);
grid_results.counts_3yr = zeros(num_periods, num_grids);
grid_results.errors_3yr = zeros(num_periods, num_grids);
grid_results.values_1yr = zeros(num_periods, num_grids);
grid_results.counts_1yr = zeros(num_periods, num_grids);
grid_results.errors_1yr = zeros(num_periods, num_grids);

% 网格元数据
for i = 1:num_grids
    grid_results.metadata(i,1) = glacier_grid(i).id;
    grid_results.metadata(i,2) = glacier_grid(i).x;
    grid_results.metadata(i,3) = glacier_grid(i).y;
    grid_results.metadata(i,4) = glacier_grid(i).Area;
end
grid_results.descrp = '网格尺度结果';
grid_results.year_labe = [config.temporal.start_year:config.temporal.end_year]';

% 3年窗口（2010-2025）
for i = 1:num_grids
    grid_mask = (data(:,12) == i) & (data(:,15) <= 40);  % 列12是grid10
    
    for j = config.temporal.start_year:config.temporal.end_year
        start_day = datenum(j-1, 1, 1) - datenum(1970, 1, 1);
        end_day = datenum(j+2, 1, 1) - datenum(1970, 1, 1);
        
        time_mask = grid_mask & (data(:,4) >= start_day) & (data(:,4) < end_day);
        
        [val, cnt, err] = compute_robust_stats(data, time_mask, 16);
        grid_results.values_3yr(j-config.temporal.start_year+1, i) = val;
        grid_results.counts_3yr(j-config.temporal.start_year+1, i) = cnt;
        grid_results.errors_3yr(j-config.temporal.start_year+1, i) = err;
    end
    
    % 1年窗口（2010-2025）
    for j = config.temporal.start_year:config.temporal.end_year
        start_day = datenum(j, 1, 1) - datenum(1970, 1, 1);
        end_day = datenum(j + 1, 1, 1) - datenum(1970, 1, 1);
        
        time_mask = grid_mask & (data(:, 4) >= start_day) & (data(:, 4) < end_day);
        
        [val, cnt, err] = compute_robust_stats(data, time_mask, 16);
        grid_results.values_1yr(j-config.temporal.start_year+1, i) = val;
        grid_results.counts_1yr(j-config.temporal.start_year+1, i) = cnt;
        grid_results.errors_1yr(j-config.temporal.start_year+1, i) = err;
    end
end

end

function elevation_seasonal = compute_elevation_band_seasonal(data, HMA22, rgiregion, glacier_grid, config)
% COMPUTE_ELEVATION_BAND_SEASONAL 高程带×季节×区域分析
%
% 参考：HMA_CryoSat2_COMPUTE_20230423m.m 第324-402行
% 使用列4(days)、列9(HMA22)、列10(RGI)、列13(TP)、列14(NASADEM)、列15(slope)、列16(dh)

fprintf('  计算季节高程带...\n');

% 季节窗口（3个月）
t = [datetime(2010, 7, 1):calmonths(1):(config.time.end_date-calmonths(3))]';%181*1
t2 = t + calmonths(3);
time_windows = [daysact(datetime(1970,1,1), t), daysact(datetime(1970,1,1), t2)];

n_periods = length(t);
n_bands = (config.spatial.elevation_range(2) - config.spatial.elevation_range(1))/config.spatial.elevation_band_width(2)+1;  % 2250-8750, 100m间隔
n_regions = 39+496;  % HMA+15RGI+22HMA22+TP

elevation_seasonal = struct();
elevation_seasonal.values = nan(n_bands, n_regions, n_periods);
elevation_seasonal.counts = nan(n_bands, n_regions, n_periods);
elevation_seasonal.errors = nan(n_bands, n_regions, n_periods);
elevation_seasonal.elevation_centers = linspace(config.spatial.elevation_range(1), config.spatial.elevation_range(2), n_bands)';
elevation_seasonal.time_windows = time_windows;

elevation_seasonal.mid_date = t + caldays(floor(daysact(t, t2)/2)); 
elevation_seasonal.descrp = '高度带结果';

half_band_width = config.spatial.elevation_band_width(2)/2;
% HMA整体
for i = 1:n_periods
    time_mask = (data(:,4) >= time_windows(i,1)) & ...
                (data(:,4) < time_windows(i,2)) & ...
                (data(:,15) < 40) & (abs(data(:,16)) < 120);
    
    for k = 1:n_bands
        elev_center = elevation_seasonal.elevation_centers(k);
        band_mask = time_mask & (data(:,14) >= elev_center-half_band_width) & (data(:,14) < elev_center+half_band_width);
        
        if any(band_mask)
            dh = data(band_mask, 16);
            elevation_seasonal.values(k, 1, i) = mean(dh, 'omitnan');
            elevation_seasonal.counts(k, 1, i) = length(dh);
            elevation_seasonal.errors(k, 1, i) = std(dh, 'omitnan') / sqrt(length(dh));
        end
    end
end

% RGI区域（列10）
if ~isempty(rgiregion)
    for j = 1:length(rgiregion)
        region_mask = (data(:,10) == j) & (data(:,15) < 40) & (abs(data(:,16)) < 120);
        data1 = data(region_mask, :);
        
        for i = 1:n_periods
            time_mask = (data1(:,4) >= time_windows(i,1)) & (data1(:,4) < time_windows(i,2));
            
            for k = 1:n_bands
                elev_center = elevation_seasonal.elevation_centers(k);
                band_mask = time_mask & (data1(:,14) >= elev_center-half_band_width) & (data1(:,14) < elev_center+half_band_width);
                
                if any(band_mask)
                    dh = data1(band_mask, 16);
                    elevation_seasonal.values(k, j+1, i) = mean(dh, 'omitnan');
                    elevation_seasonal.counts(k, j+1, i) = length(dh);
                    elevation_seasonal.errors(k, j+1, i) = std(dh, 'omitnan') / sqrt(length(dh));
                end
            end
        end
    end
end

% HMA22区域（列9）
if ~isempty(HMA22)
    for j = 1:length(HMA22)
        region_mask = (data(:,9) == j) & (data(:,15) < 40) & (abs(data(:,16)) < 120);
        data1 = data(region_mask, :);
        
        for i = 1:n_periods
            time_mask = (data1(:,4) >= time_windows(i,1)) & (data1(:,4) < time_windows(i,2));
            
            for k = 1:n_bands
                elev_center = elevation_seasonal.elevation_centers(k);
                band_mask = time_mask & (data1(:,14) >= elev_center-half_band_width) & (data1(:,14) < elev_center+half_band_width);
                
                if any(band_mask)
                    dh = data1(band_mask, 16);
                    elevation_seasonal.values(k, j+16, i) = mean(dh, 'omitnan');
                    elevation_seasonal.counts(k, j+16, i) = length(dh);
                    elevation_seasonal.errors(k, j+16, i) = std(dh, 'omitnan') / sqrt(length(dh));
                end
            end
        end
    end
end

% TP（列13）
tp_mask = (data(:,13) == 1) & (data(:,15) < 40) & (abs(data(:,16)) < 120);
data1 = data(tp_mask, :);
if ~isempty(data1)
    for i = 1:n_periods
        time_mask = (data1(:,4) >= time_windows(i,1)) & (data1(:,4) < time_windows(i,2));
        
        for k = 1:n_bands
            elev_center = elevation_seasonal.elevation_centers(k);
            band_mask = time_mask & (data1(:,14) >= elev_center-half_band_width) & (data1(:,14) < elev_center+half_band_width);
            
            if any(band_mask)
                dh = data1(band_mask, 16);
                elevation_seasonal.values(k, 39, i) = mean(dh, 'omitnan');
                elevation_seasonal.counts(k, 39, i) = length(dh);
                elevation_seasonal.errors(k, 39, i) = std(dh, 'omitnan') / sqrt(length(dh));
            end
        end
    end
end

for j = 1:length(glacier_grid)
    region_mask = (data(:,12) == j) & (data(:,15) < 40) & (abs(data(:,16)) < 120);
    data1 = data(region_mask, :);

    if isempty(data1)
        continue
    end
    
    for i = 1:n_periods
        time_mask = (data1(:,4) >= time_windows(i,1)) & (data1(:,4) < time_windows(i,2));
        
        for k = 1:n_bands
            elev_center = elevation_seasonal.elevation_centers(k);
            band_mask = time_mask & (data1(:,14) >= elev_center-half_band_width) & (data1(:,14) < elev_center+half_band_width);
            
            if any(band_mask)
                dh = data1(band_mask, 16);
                elevation_seasonal.values(k, j+39, i) = mean(dh, 'omitnan');
                elevation_seasonal.counts(k, j+39, i) = length(dh);
                elevation_seasonal.errors(k, j+39, i) = std(dh, 'omitnan') / sqrt(length(dh));
            end
        end
    end
end

% 年窗口（3年）
years = [2010:year(config.time.end_date);
        2011:year(config.time.end_date)+1];
% 两年期窗口
two_year_windows = [2010:year(config.time.end_date)-1;
                    2012:year(config.time.end_date)+1];
% 三年期窗口
three_year_windows = [2010:year(config.time.end_date)-2;
                    2013:year(config.time.end_date)+1];
whole_window = [2010;year(config.time.end_date)+1];
ref_date = config.time.reference_date;
% 组合所有年度窗口
all_windows = [years, two_year_windows, three_year_windows, whole_window];

% 转换为相对天数
year_windows = zeros(size(all_windows, 2), 2);
for i = 1:size(all_windows, 2)
    year_windows(i, 1) = daysact(ref_date, datetime(all_windows(1, i), 1, 1));
    year_windows(i, 2) = daysact(ref_date, datetime(all_windows(2, i), 1, 1));
end
n_periods = size(year_windows, 1);

elevation_seasonal.values_year = nan(n_bands, n_regions, n_periods);
elevation_seasonal.counts_year = nan(n_bands, n_regions, n_periods);
elevation_seasonal.errors_year = nan(n_bands, n_regions, n_periods);

% HMA整体
for i = 1:n_periods
    time_mask = (data(:,4) >= year_windows(i,1)) & ...
                (data(:,4) < year_windows(i,2)) & ...
                (data(:,15) < 40) & (abs(data(:,16)) < 120);
    
    for k = 1:n_bands
        elev_center = elevation_seasonal.elevation_centers(k);
        band_mask = time_mask & (data(:,14) >= elev_center-half_band_width) & (data(:,14) < elev_center+half_band_width);
        
        if any(band_mask)
            dh = data(band_mask, 16);
            elevation_seasonal.values_year(k, 1, i) = mean(dh, 'omitnan');
            elevation_seasonal.counts_year(k, 1, i) = length(dh);
            elevation_seasonal.errors_year(k, 1, i) = std(dh, 'omitnan') / sqrt(length(dh));
        end
    end
end

% RGI区域（列10）
if ~isempty(rgiregion)
    for j = 1:length(rgiregion)
        region_mask = (data(:,10) == j) & (data(:,15) < 40) & (abs(data(:,16)) < 120);
        data1 = data(region_mask, :);
        
        for i = 1:n_periods
            time_mask = (data1(:,4) >= year_windows(i,1)) & (data1(:,4) < year_windows(i,2));
            
            for k = 1:n_bands
                elev_center = elevation_seasonal.elevation_centers(k);
                band_mask = time_mask & (data1(:,14) >= elev_center-half_band_width) & (data1(:,14) < elev_center+half_band_width);
                
                if any(band_mask)
                    dh = data1(band_mask, 16);
                    elevation_seasonal.values_year(k, j+1, i) = mean(dh, 'omitnan');
                    elevation_seasonal.counts_year(k, j+1, i) = length(dh);
                    elevation_seasonal.errors_year(k, j+1, i) = std(dh, 'omitnan') / sqrt(length(dh));
                end
            end
        end
    end
end

% HMA22区域（列9）
if ~isempty(HMA22)
    for j = 1:length(HMA22)
        region_mask = (data(:,9) == j) & (data(:,15) < 40) & (abs(data(:,16)) < 120);
        data1 = data(region_mask, :);
        
        for i = 1:n_periods
            time_mask = (data1(:,4) >= year_windows(i,1)) & (data1(:,4) < year_windows(i,2));
            
            for k = 1:n_bands
                elev_center = elevation_seasonal.elevation_centers(k);
                band_mask = time_mask & (data1(:,14) >= elev_center-half_band_width) & (data1(:,14) < elev_center+half_band_width);
                
                if any(band_mask)
                    dh = data1(band_mask, 16);
                    elevation_seasonal.values_year(k, j+16, i) = mean(dh, 'omitnan');
                    elevation_seasonal.counts_year(k, j+16, i) = length(dh);
                    elevation_seasonal.errors_year(k, j+16, i) = std(dh, 'omitnan') / sqrt(length(dh));
                end
            end
        end
    end
end

% TP（列13）
tp_mask = (data(:,13) == 1) & (data(:,15) < 40) & (abs(data(:,16)) < 120);
data1 = data(tp_mask, :);
if ~isempty(data1)
    for i = 1:n_periods
        time_mask = (data1(:,4) >= year_windows(i,1)) & (data1(:,4) < year_windows(i,2));
        
        for k = 1:n_bands
            elev_center = elevation_seasonal.elevation_centers(k);
            band_mask = time_mask & (data1(:,14) >= elev_center-half_band_width) & (data1(:,14) < elev_center+half_band_width);
            
            if any(band_mask)
                dh = data1(band_mask, 16);
                elevation_seasonal.values_year(k, 39, i) = mean(dh, 'omitnan');
                elevation_seasonal.counts_year(k, 39, i) = length(dh);
                elevation_seasonal.errors_year(k, 39, i) = std(dh, 'omitnan') / sqrt(length(dh));
            end
        end
    end
end
% grid
for j = 1:length(glacier_grid)
    region_mask = (data(:,12) == j) & (data(:,15) < 40) & (abs(data(:,16)) < 120);
    data1 = data(region_mask, :);

    if isempty(data1)
        continue
    end
    
    for i = 1:n_periods
        time_mask = (data1(:,4) >= year_windows(i,1)) & (data1(:,4) < year_windows(i,2));
        
        for k = 1:n_bands
            elev_center = elevation_seasonal.elevation_centers(k);
            band_mask = time_mask & (data1(:,14) >= elev_center-half_band_width) & (data1(:,14) < elev_center+half_band_width);
            
            if any(band_mask)
                dh = data1(band_mask, 16);
                elevation_seasonal.values_year(k, j+39, i) = mean(dh, 'omitnan');
                elevation_seasonal.counts_year(k, j+39, i) = length(dh);
                elevation_seasonal.errors_year(k, j+39, i) = std(dh, 'omitnan') / sqrt(length(dh));
            end
        end
    end
end

elevation_seasonal.year_label = all_windows';
elevation_seasonal.year_windows = year_windows;
end

function area_weighted = compute_area_weighted_analysis(results, config)
% COMPUTE_AREA_WEIGHTED_ANALYSIS 面积加权计算
%
% 参考：HMA_CryoSat2_COMPUTE_20230423m.m 第405-495行

fprintf('  面积加权计算...\n');

% 加载面积数据
area1 = xlsread(config.files.elevation_area, 'HMA_RGI_100', 'A1:Q66');
area2 = xlsread(config.files.elevation_area, 'HMA22_100', 'A1:Y66');
area_grid = xlsread(config.files.elevation_area, 'grid10_100_496', 'A2:SC67');
ele = [area1, area2(:,4:25), area2(:,3), area_grid(:,2:end)];% ele_center, HMA, GTN15, HMA22, TP, 496grid

% 季节高程带面积加权
if isfield(results, 'elevation_bands')
    ice_values = results.elevation_bands.values;
    len_values = results.elevation_bands.counts;
    err_values = results.elevation_bands.errors;
    n_periods = size(ice_values, 3);
    n_regions = size(ice_values, 2);
    n_bands = size(ice_values, 1);
    
    area_weighted.seasonal_values = nan(n_regions, n_periods);
    area_weighted.seasonal_counts = nan(n_regions, n_periods);
    area_weighted.seasonal_errors = nan(n_regions, n_periods);
    
    for j = 1:n_periods
        ice = ice_values(:,:,j);
        len = len_values(:,:,j);
        err = err_values(:,:,j);
        
        for i = 1:n_regions
            area_data = [ele(:,1), ele(:,i+1), ...
                        ele(:,i+1).*ice(:,i), len(:,i), ice(:,i), err(:,i)];
            
            % 移除NaN
            area_data(isnan(area_data(:,3)), :) = [];
            
            % 质量过滤
            area_data(area_data(:,4)<10 & abs(area_data(:,5))>60, :) = [];
            %area_data(abs(area_data(:,5))>80, :) = [];
            
            if ~isempty(area_data) && sum(area_data(:,2)) > 0
                area_weighted.seasonal_values(i,j) = sum(area_data(:,3)) / sum(area_data(:,2));
                area_weighted.seasonal_counts(i,j) = sum(area_data(:,4));
                area_weighted.seasonal_errors(i,j) = mean(area_data(:,6))/sqrt(length(area_data(:,6)));
            end
        end
    end
end

area_weighted.seasonal_values = area_weighted.seasonal_values';
area_weighted.seasonal_counts = area_weighted.seasonal_counts';
area_weighted.seasonal_errors = area_weighted.seasonal_errors';
area_weighted.mid_dates = results.elevation_bands.mid_date;
ref_date = config.time.reference_date;
area_weighted.mid_days = daysact(ref_date, area_weighted.mid_dates);

% 年高程带面积加权
if isfield(results, 'elevation_bands')
    ice_values = results.elevation_bands.values_year;
    len_values = results.elevation_bands.counts_year;
    err_values = results.elevation_bands.errors_year;
    n_periods = size(ice_values, 3);
    n_regions = size(ice_values, 2);
    n_bands = size(ice_values, 1);
    
    area_weighted.yearly_values = nan(n_regions, n_periods);
    area_weighted.yearly_counts = nan(n_regions, n_periods);
    area_weighted.yearly_errors = nan(n_regions, n_periods);
    
    for j = 1:n_periods
        ice = ice_values(:,:,j);
        len = len_values(:,:,j);
        err = err_values(:,:,j);
        
        for i = 1:n_regions
            area_data = [ele(:,1), ele(:,i+1), ...
                        ele(:,i+1).*ice(:,i), len(:,i), ice(:,i), err(:,i)];
            
            % 移除NaN
            area_data(isnan(area_data(:,3)), :) = [];
            
            % 质量过滤
            area_data(area_data(:,4)<10 & abs(area_data(:,5))>60, :) = [];
            %area_data(abs(area_data(:,5))>80, :) = [];
            
            if ~isempty(area_data) && sum(area_data(:,2)) > 0
                area_weighted.yearly_values(i,j) = sum(area_data(:,3)) / sum(area_data(:,2));
                area_weighted.yearly_counts(i,j) = sum(area_data(:,4));
                area_weighted.yearly_errors(i,j) = mean(area_data(:,6))/sqrt(length(area_data(:,6)));
            end
        end
    end
end

area_weighted.yearly_values = area_weighted.yearly_values';
area_weighted.yearly_counts = area_weighted.yearly_counts';
area_weighted.yearly_errors = area_weighted.yearly_errors';
area_weighted.year_label = results.elevation_bands.year_label;

ref_date = config.time.reference_date;
year_windows = results.elevation_bands.year_windows;
dates = ref_date + caldays(fix((year_windows(:, 1)+year_windows(:, 2))/2));%列向量
area_weighted.mid_dates_year = dates;
end

function [mean_val, count, std_err] = compute_robust_stats(data, mask, dh_col)
% COMPUTE_ROBUST_STATS 计算鲁棒统计量
%
% 输入:
%   data - 数据矩阵
%   mask - 筛选掩膜
%   dh_col - 高程变化所在列（通常是16）
%
% 输出:
%   mean_val - 平均值
%   count - 数据点数
%   std_err - 标准误

if ~any(mask)
    mean_val = NaN;
    count = 0;
    std_err = NaN;
    return;
end

% 第一轮过滤
row = find(mask & abs(data(:,dh_col)) < 150);

if isempty(row)
    mean_val = NaN;
    count = 0;
    std_err = NaN;
    return;
end

dh = data(row, dh_col);
med = median(dh);

% 第二轮过滤
row0 = find(mask & abs(data(:,dh_col) - med) < 75);

if length(row0)<30
    mean_val = NaN;
    count = 0;
    std_err = NaN;
    return;
end

dh_final = data(row0, dh_col);
mean_val = mean(dh_final, 'omitnan');
count = length(row0);
std_err = std(dh_final, 'omitnan') / sqrt(count);

end

