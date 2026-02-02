function results = compute_icesat_statistics(txt_file, config)
% COMPUTE_ICESAT_STATISTICS ICESat data statistical analysis
%
% Reference: HMA_ICESat_COMPUTE_20230603.m
%
% Functions:
%   1. Annual statistics (2003–2009, ending in October each year)
%   2. 3-year moving-window statistics
%   3. Campaign statistics (19 observation windows)
%   4. Multi-region statistics (HMA, 15 RGI, 22 HMA22, TP, HMA4, HMA6)
%
% Quality-control conditions (following the original code):
%   - Column 16 (slope) <= 40 degrees
%   - Column 9 (glacier_flag) == 0 (inside glacier)
%   - abs(column 21 (dh)) < 200 (overall filter)
%   - abs(column 21 (dh)) < 100 (initial filter)
%   - abs(column 21 - median) < 50 (refined filter)
%
% Column references (note: lon is in column 1, lat in column 2):
%   - Col 5: days (since 2000-01-01)
   %   - 列9: glacier_flag (0=冰川内)
%   - Col12: hma22_region (HMA22 region ID)
%   - Col13: gtn_region (GTN15 region ID)
%   - Col16: slope (3×3 neighborhood)
%   - Col20: dh (elevation change, without PDD) ⭐ used for statistics
%   - Col21: dh_with_pdd (elevation change, with PDD)
%   - Col22: HMA4 region ID
%   - Col23: HMA6 region ID
%   - Col24: TP flag

fprintf('===ICESat数据统计分析开始===\n');

%% Load data
fprintf('加载增强数据...\n');
try
    data = load(txt_file);
    fprintf('✓ 成功加载 %d 个数据点\n', size(data, 1));
catch ME
    error('数据文件加载失败: %s', ME.message);
end

%% Load region data
fprintf('加载区域Shape文件...\n');
HMA22 = shaperead(config.files.hma22_regions);
rgiregion = shaperead(config.files.gtn_regions);
glacier_grid = shaperead(config.files.grid10);

%% Data quality control
fprintf('执行质量控制...\n');
data = apply_quality_control(data, config);
fprintf('✓ 质量控制后保留 %d 个数据点\n', size(data, 1));

%% Annual statistics (2003–2009)
fprintf('计算年度统计...\n');
results.annual = compute_annual_statistics(data, rgiregion, HMA22, config);
fprintf('✓ 年度统计完成\n');

%% 3-year moving-window statistics
fprintf('计算3年滑动窗口...\n');
results.multiyear = compute_multiyear_statistics(data, rgiregion, HMA22, config);
fprintf('✓ 3年滑动窗口统计完成\n');

%% Campaign statistics (19 observation windows)
fprintf('计算Campaign统计...\n');
results.campaigns = compute_campaign_statistics(data, rgiregion, HMA22, config);
fprintf('✓ Campaign统计完成\n');

%% Grid statistics
fprintf('计算网格统计...\n');
results.grid = compute_grid_statistics(data, glacier_grid, config);
fprintf('✓ 网格统计完成\n');

%% Elevation-band analysis (seasonal + annual scales)
fprintf('执行高程带分析（季节+年尺度）...\n');
results.elevation_bands = compute_elevation_band_seasonal(data, HMA22, rgiregion, glacier_grid, config);
fprintf('✓ 季节/年高程带分析完成\n');

%% Area-weighted calculations
fprintf('执行面积加权计算...\n');
results.area_weighted = compute_area_weighted_analysis(results, config);
fprintf('✓ 面积加权完成\n');

%% 保存元数据
results.metadata.processing_date = datetime('now');
results.metadata.data_points_original = size(load(txt_file), 1);
results.metadata.data_points_qc = size(data, 1);
results.metadata.config = config;

fprintf('===ICESat数据统计分析完成===\n');

end

function data_qc = apply_quality_control(data, config)
% APPLY_QUALITY_CONTROL Apply quality-control filters
% Reference: HMA_ICESat_COMPUTE_20230603.m, lines 48–52

fprintf('  应用质量控制条件...\n');

% Remove data with excessively large elevation changes (>200 m), using column 21 (dh)
todelete = abs(data(:,21)) > config.quality.max_height_change;
data(todelete,:) = [];
fprintf('  - 删除极端高程变化数据: %d点\n', sum(todelete));

data_qc = data;

end

function annual = compute_annual_statistics(data, rgiregion, HMA22, config)
% COMPUTE_ANNUAL_STATISTICS Annual statistics
% Reference: HMA_ICESat_COMPUTE_20230603.m, lines 85–172

years = config.temporal.start_year:config.temporal.end_year;
n_years = length(years);
n_rgi = length(rgiregion);
n_hma22 = length(HMA22);

% Initialize result matrices
% Column structure: year | start day | end day | HMA | 15 RGI | 22 HMA22 | TP
total_cols = 1 + n_rgi + n_hma22 + 1;  % 39 columns
results_values = nan(n_years, total_cols);
results_counts = nan(n_years, total_cols);
results_errors = nan(n_years, total_cols);

annual.year_label = years';
annual.time_windows = zeros(n_years, 2);

for i = 1:n_years
    year = years(i);
    
    % Time window: from year-01-01 to (year+1)-10-01 (same as original code)
    start_day = datenum(year, 1, 1) - datenum(2000, 1, 1);
    end_day = datenum(year+1, 10, 1) - datenum(2000, 1, 1);
    
    annual.time_windows(i,1) = start_day;
    annual.time_windows(i,2) = end_day;
    
    col_idx = 1;
    
    % HMA overall statistics
    [val, cnt, err] = compute_region_stats(data, start_day, end_day, [], [], config);
    results_values(i, col_idx) = val;
    results_counts(i, col_idx) = cnt;
    results_errors(i, col_idx) = err;
    col_idx = col_idx + 1;
    
    % Statistics for 15 RGI regions
    for j = 1:n_rgi
        [val, cnt, err] = compute_region_stats(data, start_day, end_day, 'rgi', j, config);
        results_values(i, col_idx) = val;
        results_counts(i, col_idx) = cnt;
        results_errors(i, col_idx) = err;
        col_idx = col_idx + 1;
    end
    
    % Statistics for 22 HMA22 regions
    for k = 1:n_hma22
        [val, cnt, err] = compute_region_stats(data, start_day, end_day, 'hma22', k, config);
        results_values(i, col_idx) = val;
        results_counts(i, col_idx) = cnt;
        results_errors(i, col_idx) = err;
        col_idx = col_idx + 1;
    end
    
    % TP region statistics
    [val, cnt, err] = compute_region_stats(data, start_day, end_day, 'tp', 1, config);
    results_values(i, col_idx) = val;
    results_counts(i, col_idx) = cnt;
    results_errors(i, col_idx) = err;
    %col_idx = col_idx + 1;
end

annual.values = results_values;
annual.counts = results_counts;
annual.errors = results_errors;

end

function multiyear = compute_multiyear_statistics(data, rgiregion, HMA22, config)
% COMPUTE_MULTIYEAR_STATISTICS 3-year moving-window statistics
% Reference: HMA_ICESat_COMPUTE_20230603.m, lines 174–260

years = config.temporal.start_year:config.temporal.end_year;
n_years = length(years);
n_rgi = length(rgiregion);
n_hma22 = length(HMA22);

% Initialize result matrices
total_cols = 1 + n_rgi + n_hma22 + 1;  % 39 columns
results_values = nan(n_years, total_cols);
results_counts = nan(n_years, total_cols);
results_errors = nan(n_years, total_cols);

multiyear.year_label = zeros(n_years, 2);
multiyear.time_windows = zeros(n_years, 2);

for i = 1:n_years
    year = years(i);

    multiyear.year_label(i,1) = year-1;
    multiyear.year_label(i,2) = year+2;
    
    % 3-year window: from (year-1)-01-01 to (year+2)-01-01
    start_day = datenum(year-1, 1, 1) - datenum(2000, 1, 1);
    end_day = datenum(year+2, 1, 1) - datenum(2000, 1, 1);
    
    multiyear.time_windows(i,1) = start_day;
    multiyear.time_windows(i,2) = end_day;
    
    col_idx = 1;
    
    % HMA overall
    [val, cnt, err] = compute_region_stats(data, start_day, end_day, [], [], config);
    results_values(i, col_idx) = val;
    results_counts(i, col_idx) = cnt;
    results_errors(i, col_idx) = err;
    col_idx = col_idx + 1;
    
    % RGI regions
    for k = 1:n_rgi
        [val, cnt, err] = compute_region_stats(data, start_day, end_day, 'rgi', k, config);
        results_values(i, col_idx) = val;
        results_counts(i, col_idx) = cnt;
        results_errors(i, col_idx) = err;
        col_idx = col_idx + 1;
    end
    
    % HMA22 regions
    for j = 1:n_hma22
        [val, cnt, err] = compute_region_stats(data, start_day, end_day, 'hma22', j, config);
        results_values(i, col_idx) = val;
        results_counts(i, col_idx) = cnt;
        results_errors(i, col_idx) = err;
        col_idx = col_idx + 1;
    end
    
    % TP region
    [val, cnt, err] = compute_region_stats(data, start_day, end_day, 'tp', 1, config);
    results_values(i, col_idx) = val;
    results_counts(i, col_idx) = cnt;
    results_errors(i, col_idx) = err;
    %col_idx = col_idx + 1;
end

multiyear.values = results_values;
multiyear.counts = results_counts;
multiyear.errors = results_errors;
multiyear.descrp = '3-year moving window';
end

function campaigns = compute_campaign_statistics(data, rgiregion, HMA22, config)
% COMPUTE_CAMPAIGN_STATISTICS Campaign statistics (19 observation windows)
% Reference: HMA_ICESat_COMPUTE_20230603.m, lines 263–359

campaign_times = config.temporal.campaigns;
n_campaigns = size(campaign_times, 1);
n_rgi = length(rgiregion);
n_hma22 = length(HMA22);

% Initialize result matrices
total_cols = 1 + n_rgi + n_hma22 + 1;
results_values = nan(n_campaigns, total_cols);
results_counts = nan(n_campaigns, total_cols);
results_errors = nan(n_campaigns, total_cols);

campaigns.date_windows = cell(n_campaigns, 2);
campaigns.mid_date = cell(n_campaigns,1);
campaigns.time_windows = zeros(n_campaigns, 2);

for i = 1:n_campaigns
    campaigns.date_windows{i,1} = config.time.reference_date+caldays(campaign_times(i,2));
    campaigns.date_windows{i,2} = config.time.reference_date+caldays(campaign_times(i,3));
    campaigns.mid_date{i,1} = config.time.reference_date+caldays(floor((campaign_times(i,2)+campaign_times(i,3))/2));

    campaigns.time_windows(i,:) = campaign_times(i,2:3);
    
    start_day = campaign_times(i, 2);
    end_day = campaign_times(i, 3);
    
    col_idx = 1;
    
    % HMA overall
    [val, cnt, err] = compute_region_stats(data, start_day, end_day, [], [], config);
    results_values(i, col_idx) = val;
    results_counts(i, col_idx) = cnt;
    results_errors(i, col_idx) = err;
    col_idx = col_idx + 1;
    
    % RGI regions
    for k = 1:n_rgi
        [val, cnt, err] = compute_region_stats(data, start_day, end_day, 'rgi', k, config);
        results_values(i, col_idx) = val;
        results_counts(i, col_idx) = cnt;
        results_errors(i, col_idx) = err;
        col_idx = col_idx + 1;
    end
    
    % HMA22 regions
    for j = 1:n_hma22
        [val, cnt, err] = compute_region_stats(data, start_day, end_day, 'hma22', j, config);
        results_values(i, col_idx) = val;
        results_counts(i, col_idx) = cnt;
        results_errors(i, col_idx) = err;
        col_idx = col_idx + 1;
    end
    
    % TP region
    [val, cnt, err] = compute_region_stats(data, start_day, end_day, 'tp', 1, config);
    results_values(i, col_idx) = val;
    results_counts(i, col_idx) = cnt;
    results_errors(i, col_idx) = err;
    %col_idx = col_idx + 1;
end

campaigns.values = results_values;
campaigns.counts = results_counts;
campaigns.errors = results_errors;
campaigns.descrp = 'Campaign results';

end

function grid_results = compute_grid_statistics(data, glacier_grid, config)
% COMPUTE_GRID_STATISTICS Grid-based statistics
%
% Reference: compute_cryosat2_statistics.m, function compute_grid_statistics
% Uses columns 5 (days), 14 (grid_region), 16 (slope), 20 (dh)

fprintf('  开始网格统计（3年和1年窗口）...\n');

num_periods = config.temporal.end_year - config.temporal.start_year + 1;
num_grids = length(glacier_grid);

grid_results = struct();
grid_results.metadata = zeros(num_grids, 4);
grid_results.values_3yr = nan(num_periods, num_grids);
grid_results.counts_3yr = nan(num_periods, num_grids);
grid_results.errors_3yr = nan(num_periods, num_grids);
grid_results.values_1yr = nan(num_periods, num_grids);
grid_results.counts_1yr = nan(num_periods, num_grids);
grid_results.errors_1yr = nan(num_periods, num_grids);

% Grid metadata
for i = 1:num_grids
    grid_results.metadata(i,1) = glacier_grid(i).id;
    grid_results.metadata(i,2) = glacier_grid(i).x;
    grid_results.metadata(i,3) = glacier_grid(i).y;
    grid_results.metadata(i,4) = glacier_grid(i).Area;
end

grid_results.year_label = [config.temporal.start_year:config.temporal.end_year]';

% 3-year window statistics
for i = 1:num_grids
    grid_mask = (data(:,14) == i) & (data(:,16) <= 40);  % Col14 is grid_region, Col16 is slope
    
    for j = config.temporal.start_year:config.temporal.end_year
        start_day = datenum(j-1, 1, 1) - datenum(2000, 1, 1);
        end_day = datenum(j+2, 1, 1) - datenum(2000, 1, 1);
        
        time_mask = grid_mask & (data(:,5) >= start_day) & (data(:,5) < end_day);
        
        [val, cnt, err] = compute_robust_stats(data, time_mask, 21);
        grid_results.values_3yr(j-config.temporal.start_year+1, i) = val;
        grid_results.counts_3yr(j-config.temporal.start_year+1, i) = cnt;
        grid_results.errors_3yr(j-config.temporal.start_year+1, i) = err;
    end
    
    % 1-year window statistics
    for j = config.temporal.start_year:config.temporal.end_year
        start_day = datenum(j, 1, 1) - datenum(2000, 1, 1);
        end_day = datenum(j + 1, 1, 1) - datenum(2000, 1, 1);
        
        time_mask = grid_mask & (data(:, 5) >= start_day) & (data(:, 5) < end_day);
        
        [val, cnt, err] = compute_robust_stats(data, time_mask, 21);
        grid_results.values_1yr(j-config.temporal.start_year+1, i) = val;
        grid_results.counts_1yr(j-config.temporal.start_year+1, i) = cnt;
        grid_results.errors_1yr(j-config.temporal.start_year+1, i) = err;
    end
end

grid_results.descrp = 'Grid-scale results';

end

function elevation_seasonal = compute_elevation_band_seasonal(data, HMA22, rgiregion, glacier_grid, config)
% COMPUTE_ELEVATION_BAND_SEASONAL Elevation-band × season × region analysis
%
% Reference: compute_cryosat2_statistics.m, function compute_elevation_band_seasonal
% Uses columns 5 (days), 12 (hma22), 13 (gtn), 14 (grid), 15 (nasadem),
% 16 (slope), 21 (dh), 24 (TP)

fprintf('  计算季节高程带（基于campaigns）...\n');

% ICESat uses campaigns as time windows (19 observation windows)
campaign_times = config.temporal.campaigns;
n_periods = size(campaign_times, 1);

n_bands = (config.spatial.elevation_range(2) - config.spatial.elevation_range(1))/config.spatial.elevation_band_width + 1;
n_regions = 39 + 496;  % HMA + 15 RGI + 22 HMA22 + TP + 496 grid cells

elevation_seasonal = struct();
elevation_seasonal.values = nan(n_bands, n_regions, n_periods);
elevation_seasonal.counts = nan(n_bands, n_regions, n_periods);
elevation_seasonal.errors = nan(n_bands, n_regions, n_periods);
elevation_seasonal.elevation_centers = linspace(config.spatial.elevation_range(1), config.spatial.elevation_range(2), n_bands)';
elevation_seasonal.time_windows = campaign_times;
elevation_seasonal.mid_date = datetime(2000,1,1) + days(floor((campaign_times(:,2) + campaign_times(:,3))/2));
elevation_seasonal.date_windows(:,1) = datetime(2000,1,1) + days(campaign_times(:,2));
elevation_seasonal.date_windows(:,2) = datetime(2000,1,1) + days(campaign_times(:,3));
elevation_seasonal.descrp = 'Elevation-band results (based on campaigns)';

half_band_width = config.spatial.elevation_band_width / 2;

% HMA overall
for i = 1:n_periods
    time_mask = (data(:,5) >= campaign_times(i,2)) & ...
                (data(:,5) < campaign_times(i,3)) & ...
                (data(:,16) < 40) & (abs(data(:,21)) < 100);  % Col16: slope, Col21: dh
    
    for k = 1:n_bands
        elev_center = elevation_seasonal.elevation_centers(k);
        band_mask = time_mask & (data(:,15) >= elev_center-half_band_width) & (data(:,15) < elev_center+half_band_width);
        
        if any(band_mask)
            dh = data(band_mask, 21);
            elevation_seasonal.values(k, 1, i) = mean(dh, 'omitnan');
            elevation_seasonal.counts(k, 1, i) = length(dh);
            elevation_seasonal.errors(k, 1, i) = std(dh, 'omitnan') / sqrt(length(dh));
        end
    end
end

% RGI regions (column 13)
if ~isempty(rgiregion)
    for j = 1:length(rgiregion)
        region_mask = (data(:,13) == j) & (data(:,16) < 40) & (abs(data(:,21)) < 100);
        data1 = data(region_mask, :);
        
        for i = 1:n_periods
            time_mask = (data1(:,5) >= campaign_times(i,2)) & (data1(:,5) < campaign_times(i,3));
            
            for k = 1:n_bands
                elev_center = elevation_seasonal.elevation_centers(k);
                band_mask = time_mask & (data1(:,15) >= elev_center-half_band_width) & (data1(:,15) < elev_center+half_band_width);
                
                if any(band_mask)
                    dh = data1(band_mask, 21);
                    elevation_seasonal.values(k, j+1, i) = mean(dh, 'omitnan');
                    elevation_seasonal.counts(k, j+1, i) = length(dh);
                    elevation_seasonal.errors(k, j+1, i) = std(dh, 'omitnan') / sqrt(length(dh));
                end
            end
        end
    end
end

% HMA22 regions (column 12)
if ~isempty(HMA22)
    for j = 1:length(HMA22)
        region_mask = (data(:,12) == j) & (data(:,16) < 40) & (abs(data(:,21)) < 100);
        data1 = data(region_mask, :);
        
        for i = 1:n_periods
            time_mask = (data1(:,5) >= campaign_times(i,2)) & (data1(:,5) < campaign_times(i,3));
            
            for k = 1:n_bands
                elev_center = elevation_seasonal.elevation_centers(k);
                band_mask = time_mask & (data1(:,15) >= elev_center-half_band_width) & (data1(:,15) < elev_center+half_band_width);
                
                if any(band_mask)
                    dh = data1(band_mask, 21);
                    elevation_seasonal.values(k, j+16, i) = mean(dh, 'omitnan');
                    elevation_seasonal.counts(k, j+16, i) = length(dh);
                    elevation_seasonal.errors(k, j+16, i) = std(dh, 'omitnan') / sqrt(length(dh));
                end
            end
        end
    end
end

% TP (column 24)
tp_mask = (data(:,24) == 1) & (data(:,16) < 40) & (abs(data(:,21)) < 100);
data1 = data(tp_mask, :);
if ~isempty(data1)
    for i = 1:n_periods
        time_mask = (data1(:,5) >= campaign_times(i,2)) & (data1(:,5) < campaign_times(i,3));
        
        for k = 1:n_bands
            elev_center = elevation_seasonal.elevation_centers(k);
            band_mask = time_mask & (data1(:,15) >= elev_center-half_band_width) & (data1(:,15) < elev_center+half_band_width);
            
            if any(band_mask)
                dh = data1(band_mask, 21);
                elevation_seasonal.values(k, 39, i) = mean(dh, 'omitnan');
                elevation_seasonal.counts(k, 39, i) = length(dh);
                elevation_seasonal.errors(k, 39, i) = std(dh, 'omitnan') / sqrt(length(dh));
            end
        end
    end
end

% Grid496 (column 14)
for j = 1:length(glacier_grid)
    region_mask = (data(:,14) == j) & (data(:,16) < 40) & (abs(data(:,21)) < 100);
    data1 = data(region_mask, :);
    
    if isempty(data1)
        continue
    end
    
    for i = 1:n_periods
        time_mask = (data1(:,5) >= campaign_times(i,2)) & (data1(:,5) < campaign_times(i,3));
        
        for k = 1:n_bands
            elev_center = elevation_seasonal.elevation_centers(k);
            band_mask = time_mask & (data1(:,15) >= elev_center-half_band_width) & (data1(:,15) < elev_center+half_band_width);
            
            if any(band_mask)
                dh = data1(band_mask, 21);
                elevation_seasonal.values(k, j+39, i) = mean(dh, 'omitnan');
                elevation_seasonal.counts(k, j+39, i) = length(dh);
                elevation_seasonal.errors(k, j+39, i) = std(dh, 'omitnan') / sqrt(length(dh));
            end
        end
    end
end

% Elevation bands on annual scales (using annual windows)
fprintf('  计算年尺度高程带...\n');

% Annual windows (3-year)
years = [2003:year(config.time.end_date);
        2004:year(config.time.end_date)+1];
% Two-year windows
two_year_windows = [2003:year(config.time.end_date)-1;
                    2005:year(config.time.end_date)+1];
% Three-year windows
three_year_windows = [2003:year(config.time.end_date)-2;
                    2006:year(config.time.end_date)+1];
whole_window = [2003;year(config.time.end_date)+1];
ref_date = config.time.reference_date;
% Combine all annual windows
all_windows = [years, two_year_windows, three_year_windows, whole_window];

% Convert to days relative to reference date
year_windows = zeros(size(all_windows, 2), 2);
for i = 1:size(all_windows, 2)
    year_windows(i, 1) = daysact(ref_date, datetime(all_windows(1, i), 1, 1));
    year_windows(i, 2) = daysact(ref_date, datetime(all_windows(2, i), 1, 1));
end
n_year_periods = size(year_windows, 1);

elevation_seasonal.values_year = nan(n_bands, n_regions, n_year_periods);
elevation_seasonal.counts_year = nan(n_bands, n_regions, n_year_periods);
elevation_seasonal.errors_year = nan(n_bands, n_regions, n_year_periods);

% HMA overall
for i = 1:n_year_periods
    time_mask = (data(:,5) >= year_windows(i,1)) & ...
                (data(:,5) < year_windows(i,2)) & ...
                (data(:,16) < 40) & (abs(data(:,21)) < 100);
    
    for k = 1:n_bands
        elev_center = elevation_seasonal.elevation_centers(k);
        band_mask = time_mask & (data(:,15) >= elev_center-half_band_width) & (data(:,15) < elev_center+half_band_width);
        
        if any(band_mask)
            dh = data(band_mask, 21);
            elevation_seasonal.values_year(k, 1, i) = mean(dh, 'omitnan');
            elevation_seasonal.counts_year(k, 1, i) = length(dh);
            elevation_seasonal.errors_year(k, 1, i) = std(dh, 'omitnan') / sqrt(length(dh));
        end
    end
end

% RGI regions
for j = 1:length(rgiregion)
    region_mask = (data(:,13) == j) & (data(:,16) < 40) & (abs(data(:,21)) < 100);
    data1 = data(region_mask, :);
    
    for i = 1:n_year_periods
        time_mask = (data1(:,5) >= year_windows(i,1)) & (data1(:,5) < year_windows(i,2));
        
        for k = 1:n_bands
            elev_center = elevation_seasonal.elevation_centers(k);
            band_mask = time_mask & (data1(:,15) >= elev_center-half_band_width) & (data1(:,15) < elev_center+half_band_width);
            
            if any(band_mask)
                dh = data1(band_mask, 21);
                elevation_seasonal.values_year(k, j+1, i) = mean(dh, 'omitnan');
                elevation_seasonal.counts_year(k, j+1, i) = length(dh);
                elevation_seasonal.errors_year(k, j+1, i) = std(dh, 'omitnan') / sqrt(length(dh));
            end
        end
    end
end

% HMA22 regions
for j = 1:length(HMA22)
    region_mask = (data(:,12) == j) & (data(:,16) < 40) & (abs(data(:,21)) < 100);
    data1 = data(region_mask, :);
    
    for i = 1:n_year_periods
        time_mask = (data1(:,5) >= year_windows(i,1)) & (data1(:,5) < year_windows(i,2));
        
        for k = 1:n_bands
            elev_center = elevation_seasonal.elevation_centers(k);
            band_mask = time_mask & (data1(:,15) >= elev_center-half_band_width) & (data1(:,15) < elev_center+half_band_width);
            
            if any(band_mask)
                dh = data1(band_mask, 21);
                elevation_seasonal.values_year(k, j+16, i) = mean(dh, 'omitnan');
                elevation_seasonal.counts_year(k, j+16, i) = length(dh);
                elevation_seasonal.errors_year(k, j+16, i) = std(dh, 'omitnan') / sqrt(length(dh));
            end
        end
    end
end

% TP
tp_mask = (data(:,24) == 1) & (data(:,16) < 40) & (abs(data(:,21)) < 100);
data1 = data(tp_mask, :);
if ~isempty(data1)
    for i = 1:n_year_periods
        time_mask = (data1(:,5) >= year_windows(i,1)) & (data1(:,5) < year_windows(i,2));
        
        for k = 1:n_bands
            elev_center = elevation_seasonal.elevation_centers(k);
            band_mask = time_mask & (data1(:,15) >= elev_center-half_band_width) & (data1(:,15) < elev_center+half_band_width);
            
            if any(band_mask)
                dh = data1(band_mask, 21);
                elevation_seasonal.values_year(k, 39, i) = mean(dh, 'omitnan');
                elevation_seasonal.counts_year(k, 39, i) = length(dh);
                elevation_seasonal.errors_year(k, 39, i) = std(dh, 'omitnan') / sqrt(length(dh));
            end
        end
    end
end

% Grid496
for j = 1:length(glacier_grid)
    region_mask = (data(:,14) == j) & (data(:,16) < 40) & (abs(data(:,21)) < 100);
    data1 = data(region_mask, :);
    
    if isempty(data1)
        continue
    end
    
    for i = 1:n_year_periods
        time_mask = (data1(:,5) >= year_windows(i,1)) & (data1(:,5) < year_windows(i,2));
        
        for k = 1:n_bands
            elev_center = elevation_seasonal.elevation_centers(k);
            band_mask = time_mask & (data1(:,15) >= elev_center-half_band_width) & (data1(:,15) < elev_center+half_band_width);
            
            if any(band_mask)
                dh = data1(band_mask, 21);
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
% 参考：compute_cryosat2_statistics.m 的 compute_area_weighted_analysis

fprintf('  面积加权计算...\n');

% 加载面积数据
area1 = xlsread(config.files.elevation_area, 'HMA_RGI_100', 'A1:Q66');
area2 = xlsread(config.files.elevation_area, 'HMA22_100', 'A1:Y66');
area_grid = xlsread(config.files.elevation_area, 'grid10_100_496', 'A2:SC67');
ele = [area1, area2(:,4:25), area2(:,3), area_grid(:,2:end)];  % ele_center, HMA, GTN15, HMA22, TP, 496grid

% 季节（campaigns）高程带面积加权
if isfield(results, 'elevation_bands')
    ice_values = results.elevation_bands.values;
    len_values = results.elevation_bands.counts;
    err_values = results.elevation_bands.errors;
    n_periods = size(ice_values, 3);
    n_regions = size(ice_values, 2);
    
    area_weighted.seasonal_values = nan(n_periods, n_regions);
    area_weighted.seasonal_counts = nan(n_periods, n_regions);
    area_weighted.seasonal_errors = nan(n_periods, n_regions);
    
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
            area_data(area_data(:,4)<20 & abs(area_data(:,5))>60, :) = [];
            
            if ~isempty(area_data) && sum(area_data(:,2)) > 0
                area_weighted.seasonal_values(j,i) = sum(area_data(:,3)) / sum(area_data(:,2));
                area_weighted.seasonal_counts(j,i) = sum(area_data(:,4));
                area_weighted.seasonal_errors(j,i) = mean(area_data(:,6))/sqrt(length(area_data(:,6)));
            end
        end
    end
end

area_weighted.mid_dates = results.elevation_bands.mid_date;
area_weighted.mid_days = days(area_weighted.mid_dates - datetime(2000,1,1));

% 年高程带面积加权
if isfield(results, 'elevation_bands') && isfield(results.elevation_bands, 'values_year')
    ice_values = results.elevation_bands.values_year;
    len_values = results.elevation_bands.counts_year;
    err_values = results.elevation_bands.errors_year;
    n_periods = size(ice_values, 3);
    n_regions = size(ice_values, 2);
    
    area_weighted.yearly_values = nan(n_periods, n_regions);
    area_weighted.yearly_counts = nan(n_periods, n_regions);
    area_weighted.yearly_errors = nan(n_periods, n_regions);
    
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
            area_data(area_data(:,4)<20 & abs(area_data(:,5))>60, :) = [];
            
            if ~isempty(area_data) && sum(area_data(:,2)) > 0
                area_weighted.yearly_values(j,i) = sum(area_data(:,3)) / sum(area_data(:,2));
                area_weighted.yearly_counts(j,i) = sum(area_data(:,4));
                area_weighted.yearly_errors(j,i) = mean(area_data(:,6))/sqrt(length(area_data(:,6)));
            end
        end
    end
end

area_weighted.year_label = results.elevation_bands.year_label;

end

function [mean_val, count, std_err] = compute_robust_stats(data, mask, dh_col)
% COMPUTE_ROBUST_STATS 计算鲁棒统计量
%
% 输入:
%   data - 数据矩阵
%   mask - 筛选掩膜
%   dh_col - 高程变化所在列（ICESat为21）
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

% 第一轮过滤（abs(dh) < 100）
row = find(mask & abs(data(:,dh_col)) < 100);

if isempty(row)
    mean_val = NaN;
    count = 0;
    std_err = NaN;
    return;
end

dh = data(row, dh_col);
med = median(dh);

% 第二轮过滤（abs(dh - median) < 50）
row0 = find(mask & abs(data(:,dh_col) - med) < 50);

if length(row0) < 10  % 最少10个点
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

function [mean_val, count_val, error_val] = compute_region_stats(data, start_day, end_day, region_type, region_id, config)
% COMPUTE_REGION_STATS 计算单个区域的统计值
% 两步筛选：先初步筛选abs(dh)<100，再相对中值筛选abs(dh-median)<50

% 基础质量控制条件
base_mask = (data(:,5) >= start_day) & (data(:,5) < end_day) & ...
            (data(:,16) <= config.quality.max_slope) & ...  % 列16: slope (3×3)
            (data(:,9) == 0);  % 列9: glacier_flag==0表示在冰川内

% 区域筛选
if isempty(region_type)
    % HMA整体，无需额外筛选
    region_mask = true(size(data, 1), 1);
elseif strcmp(region_type, 'rgi')
    % GTN/RGI区域，使用列13
    region_mask = (data(:,13) == region_id);
elseif strcmp(region_type, 'hma22')
    region_mask = (data(:,12) == region_id);  % 列12: hma22_region
elseif strcmp(region_type, 'tp')
    region_mask = (data(:,24) == 1);  % 列24: TP标识
elseif strcmp(region_type, 'hma4')
    region_mask = (data(:,22) == region_id);  % 列22: HMA4
elseif strcmp(region_type, 'hma6')
    region_mask = (data(:,23) == region_id);  % 列23: HMA6
else
    region_mask = false(size(data, 1), 1);
end

% 第一步：初步筛选（abs(dh) < 100），使用列21 (dh)
row = find(base_mask & region_mask & abs(data(:,21)) < config.quality.height_change_filter1);

if isempty(row)
    mean_val = NaN;
    count_val = 0;
    error_val = NaN;
    return;
end

dh = data(row, 21);
med = median(dh);

% 第二步：相对中值筛选（abs(dh - median) < 50）
row0 = find(base_mask & region_mask & abs(data(:,21) - med) < config.quality.height_change_filter2);

if isempty(row0)
    mean_val = NaN;
    count_val = 0;
    error_val = NaN;
else
    dh0 = data(row0, 21);
    mean_val = mean(dh0);
    count_val = length(row0);
    error_val = std(dh0) / sqrt(length(row0));
end

end

