%% ICESat数据处理主流程脚本
%
% 参考：HMA_ICESat_COMPUTE_20230603.m + revised_CryoSat2/CryoSat2_Main_Processing_Pipeline.m
%
% 功能：高亚洲山地（HMA）ICESat数据完整处理流程
% 重构版本：v1.0
% 数据时间：2003-2009年
%
% 关键差异：ICESat使用3×3邻域平均NASADEM（vs CryoSat-2的11×11）
%
% 处理流程：
% 1. 原始数据提取（ICESat二进制 → TXT） - 需要手动运行或已有数据
% 2. 数据增强（3×3邻域NASADEM + 区域标识）
% 3. 统计分析（年度/3年窗口/Campaign）
% 4. 结果输出

%% 初始化
clear; clc; close all;

% 添加函数路径
addpath('functions');

% 加载配置
config = load_icesat_config();

fprintf('=== ICESat 数据处理流程开始 ===\n');
fprintf('重要提示：ICESat使用3×3邻域平均NASADEM\n');
fprintf('数据时间范围：2003-2009年\n');
fprintf('配置加载完成，开始数据处理...\n\n');

%% 第一阶段：原始数据提取
%
% 注意：ICESat原始数据提取通常需要从GLA14二进制文件中提取
% 如果已有提取好的TXT文件，可以跳过此步骤
%
% 此处假设已有文件：HMA_ICESat_rgiregion.txt

fprintf('--- 第一阶段：原始数据提取 ---\n');

track_txt_path = fullfile(config.paths.output_dir, config.files.track_txt);

if ~exist(track_txt_path, 'file')
    warning('未找到原始数据文件: %s\n', track_txt_path);
    [track_data, track_data_outrgi] = extract_icesat_tracks(config);
    dlmwrite(track_txt_path, track_data, 'precision', 10);

    dlmwrite(fullfile(config.paths.output_dir, 'HMA_ICESat_norgiregion.txt'), track_data_outrgi, 'precision', 10);

    fprintf('完成ICESat冰川数据提取至 %s\n\n', track_txt_path);
else
    fprintf('✓ 找到原始数据文件\n');
    fprintf('文件路径: %s\n\n', track_txt_path);
end

%% 第二阶段：数据增强
%
% 关键步骤：使用3×3邻域平均NASADEM（ICESat特有）
% 添加网格、坡度、区域标识，计算高程变化

fprintf('--- 第二阶段：数据增强（3×3邻域NASADEM） ---\n');

enhanced_txt_path = fullfile(config.paths.output_dir, config.files.enhanced_txt);

if exist(enhanced_txt_path, 'file') && ~config.options.save_intermediate
    fprintf('检测到已存在的增强数据文件，跳过增强步骤\n');
    fprintf('文件路径: %s\n\n', enhanced_txt_path);
else
    enhanced_data = enhance_icesat_data(track_txt_path, config);
    
    fprintf('✓ 地形参数添加完成（3×3邻域）：NASADEM、slope、aspect\n');
    fprintf('✓ h_li2000计算完成\n');
    fprintf('✓ pdd校正添加完成\n');
    fprintf('✓ 高程变化计算完成（dh和dh_with_pdd）\n');
    fprintf('✓ 区域标识完成（HMA4, HMA6, TP）\n');
    
    % 保存
    if config.options.save_intermediate
        dlmwrite(enhanced_txt_path, enhanced_data, 'precision', 10);
        fprintf('✓ 增强数据已保存至: %s\n', enhanced_txt_path);
    end
end

fprintf('\n');

%% 第三阶段：统计分析
%
% 多时间尺度统计分析：
% - 年度统计（2003-2009，每年到次年10月）
% - 3年滑动窗口
% - Campaign统计（19个观测窗口）

fprintf('--- 第三阶段：统计分析 ---\n');

results = compute_icesat_statistics(enhanced_txt_path, config);

fprintf('✓ 年度统计分析完成（2003-2009）\n');
fprintf('✓ 3年滑动窗口统计完成\n');
fprintf('✓ Campaign统计完成（19个观测窗口）\n');

% 保存结果
results_file = fullfile(config.paths.output_dir, config.files.results_mat);
save(results_file, 'results', '-v7.3');
fprintf('✓ 分析结果已保存至: %s\n', results_file);

fprintf('\n');

%% 第四阶段：结果导出为Excel
%
% 将统计结果导出为Excel文件，便于后续分析

fprintf('--- 第四阶段：结果导出 ---\n');

%% 季节结果保存为excel
gtn = shaperead(config.files.gtn_regions);
hma22 = shaperead(config.files.hma22_regions);
header = {'start_day', 'end_day', 'mid_day', 'HMA', gtn.fullname, hma22.Name, 'TP'};
for i = 1:496
    header{length(header)+1} = i;
end

season_file = fullfile(config.paths.output_dir, 'ICESat_seasonal_elevation_results.xlsx');

numeric_data = [results.elevation_bands.time_windows(:,2:3), results.area_weighted.mid_days, results.area_weighted.seasonal_values];
data_as_cell = num2cell(numeric_data);
output_data = [header; data_as_cell];
writecell(output_data, season_file, 'Sheet', 'campaign_value');

numeric_data = [results.elevation_bands.time_windows(:,2:3), results.area_weighted.mid_days, results.area_weighted.seasonal_counts];
data_as_cell = num2cell(numeric_data);
output_data = [header; data_as_cell];
writecell(output_data, season_file, 'Sheet', 'campaign_count');

numeric_data = [results.elevation_bands.time_windows(:,2:3), results.area_weighted.mid_days, results.area_weighted.seasonal_errors];
data_as_cell = num2cell(numeric_data);
output_data = [header; data_as_cell];
writecell(output_data, season_file, 'Sheet', 'campaign_uncert');
%% 年结果保存为excel
header = {'start_year', 'end_year(not included)', 'HMA', gtn.fullname, hma22.Name, 'TP'};
for i = 1:496
    header{length(header)+1} = i;
end

year_file = fullfile(config.paths.output_dir, 'ICESat_yearly_elevation_results.xlsx');

numeric_data = [results.area_weighted.year_label, results.area_weighted.yearly_values];
data_as_cell = num2cell(numeric_data);
output_data = [header; data_as_cell];
writecell(output_data, year_file, 'Sheet', 'year_value');

numeric_data = [results.area_weighted.year_label, results.area_weighted.yearly_counts];
data_as_cell = num2cell(numeric_data);
output_data = [header; data_as_cell];
writecell(output_data, year_file, 'Sheet', 'year_count');

numeric_data = [results.area_weighted.year_label, results.area_weighted.yearly_errors];
data_as_cell = num2cell(numeric_data);
output_data = [header; data_as_cell];
writecell(output_data, year_file, 'Sheet', 'year_uncert');

%% 输出各个区域的高度带结果：ele_bands * year
all_region_names = {'HMA', gtn.fullname, hma22.Name, 'TP'};
for i = 1:length(all_region_names)
    region_name = all_region_names{i};
    region_name = strrep(region_name, '/', '-');
    data_3d = results.elevation_bands.values_year(:, i, :);
    data = squeeze(data_3d);
    
    data = [results.elevation_bands.elevation_centers, data];
    data_as_cell = num2cell(data);

    header = cell(3, size(data_as_cell, 2));
    header{1,1} = 'start_year';
    header{2,1} = 'end_year(not included)';
    header{3,1} = 'mid_year';
    for j = 2:size(data_as_cell, 2)
        header{1,j} = results.elevation_bands.year_label(j-1, 1);
        header{2,j} = results.elevation_bands.year_label(j-1, 2);
        header{3,j} = (header{1,j}+header{2,j})/2;
    end
    
    output_data = [header; data_as_cell];
    save_path = fullfile(config.paths.ele_band_dir, ['【', num2str(i),'】',region_name, '.xlsx']);
    writecell(output_data, save_path, 'Sheet', 'ele_band_value');

    data_3d = results.elevation_bands.counts_year(:, i, :);
    data = squeeze(data_3d);
    data = [results.elevation_bands.elevation_centers, data];
    data_as_cell = num2cell(data);
    output_data = [header; data_as_cell];
    writecell(output_data, save_path, 'Sheet', 'ele_band_num');

    data_3d = results.elevation_bands.errors_year(:, i, :);
    data = squeeze(data_3d);
    data = [results.elevation_bands.elevation_centers, data];
    data_as_cell = num2cell(data);
    output_data = [header; data_as_cell];
    writecell(output_data, save_path, 'Sheet', 'ele_band_uncert');
end

%% 处理完成总结

fprintf('=== ICESat 数据处理流程完成 ===\n\n');

% 显示处理摘要
display_processing_summary(results, config);

fprintf('\n所有输出文件已保存在: %s\n', config.paths.output_dir);
fprintf('处理流程执行完毕。\n');

%% 清理
fprintf('\n清理临时变量...\n');
clear enhanced_data;
fprintf('流程结束。\n');

%% 辅助函数

function display_processing_summary(results, config)
% DISPLAY_PROCESSING_SUMMARY 显示处理摘要

fprintf('--- 处理摘要 ---\n');

% 数据信息
if isfield(results, 'metadata')
    fprintf('数据信息:\n');
    fprintf('  原始数据点数: %d\n', results.metadata.data_points_original);
    fprintf('  质量控制后: %d\n', results.metadata.data_points_qc);
    fprintf('  数据保留率: %.1f%%\n', ...
        100 * results.metadata.data_points_qc / results.metadata.data_points_original);
end

% 年度统计
if isfield(results, 'annual')
    fprintf('\n年度统计:\n');
    fprintf('  处理年份: %d - %d\n', config.temporal.start_year, config.temporal.end_year);
    
    % HMA整体趋势
    valid = ~isnan(results.annual.values(:,1));
    if sum(valid) > 1
        years = results.annual.year_label;
        dh = results.annual.values(valid, 1);
        p = polyfit(years, dh, 1);
        fprintf('  HMA整体趋势: %.3f m/yr\n', p(1));
        fprintf('  平均数据点数: %.0f 点/年\n', mean(results.annual.counts(valid, 1)));
    end
end

% 3年窗口统计
if isfield(results, 'multiyear')
    fprintf('\n3年窗口统计:\n');
    fprintf('  窗口数: %d\n', size(results.multiyear.values, 1));
end

% Campaign统计
if isfield(results, 'campaigns')
    fprintf('\nCampaign统计:\n');
    fprintf('  观测窗口数: %d\n', size(results.campaigns.values, 1));
end

% 处理时间
if isfield(results, 'metadata')
    fprintf('\n处理信息:\n');
    fprintf('  处理时间: %s\n', results.metadata.processing_date);
end

fprintf('\n');

end

