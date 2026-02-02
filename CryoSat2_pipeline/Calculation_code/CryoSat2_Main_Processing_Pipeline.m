%% CryoSat-2数据处理主流程脚本
% 
% 参考：revised_ICESat2_code/ICESat2_Main_Processing_Pipeline.m
%
% 功能：高亚洲山地（HMA）CryoSat-2数据完整处理流程
% 重构版本：v2.0（正确版本）
% 参考代码：revised_ICESat2_code + HMA_footprint_track_20230213.m
% 
% 关键差异：CryoSat-2使用11×11邻域平均NASADEM（vs ICESat-2的3×3）
%
% 处理流程：
% 1. 原始数据提取（L2I/TEMPO NetCDF → TXT）
% 2. 数据增强（11×11邻域NASADEM + 区域标识）
% 3. 统计分析（年度/空间/高程带）
% 4. 可视化输出
if isempty(gcp('nocreate')) % 如果没有并行池在运行
    parpool(12); % 创建一个包含5个 worker 的并行池
end
%% 初始化
clear; clc; close all;

% 添加函数路径
addpath('functions');

% 加载配置
config = load_cryosat2_config();

fprintf('=== CryoSat-2 数据处理流程开始 ===\n');
fprintf('重要提示：CryoSat-2使用11×11邻域平均NASADEM\n');
fprintf('配置加载完成，开始数据处理...\n\n');

%% 第一阶段：原始数据提取
%
% 从CryoSat-2 L2I/TEMPO NetCDF文件中提取轨道数据
% 分为RGI冰川区域内和区域外两部分

fprintf('--- 第一阶段：原始数据提取 ---\n');

% 提取RGI冰川区域数据
track_txt_path = fullfile(config.paths.output_dir, config.files.track_txt_tempo);

if exist(track_txt_path, 'file') && ~config.options.save_intermediate
    fprintf('检测到已存在的轨道数据文件，跳过提取步骤\n');
    fprintf('文件路径: %s\n\n', track_txt_path);
else
    try
        fprintf('提取RGI冰川区域内的足迹数据...\n');
        track_data = extract_cryosat2_tracks(config);
        fprintf('✓ 轨道数据提取完成，共 %d 个数据点\n', size(track_data, 1));
        
        % 保存
        if config.options.save_intermediate
            dlmwrite(track_txt_path, track_data, 'precision', 10);%每个数值保留 10 位有效数字
            fprintf('✓ 数据已保存至: %s\n', track_txt_path);
        end
        
    catch ME
        fprintf('❌ 轨道数据提取失败: %s\n', ME.message);
        fprintf('错误堆栈:\n');
        disp(ME.stack);
        return;
    end
end

%% 可选：提取非RGI区域数据（用于对比研究）
log_fid = fopen(config.files.process_log, 'w');
fprintf(log_fid, '【%s】 Start to process\n', datetime('now'));

if config.options.save_intermediate
    fprintf('\n提取非RGI区域数据（可选）...\n');
    try
        extract_nonrgi_cryosat2_tracks(config);
        fprintf('✓ 非RGI区域数据提取完成\n');
    catch ME
        fprintf('⚠ 非RGI数据提取失败（可忽略）: %s\n', ME.message);
    end
end

fprintf('\n');

%% 第二阶段：数据增强
%
% 关键步骤：使用11×11邻域平均NASADEM（CryoSat-2特有）
% 添加坡度、区域标识，计算高程变化

fprintf('--- 第二阶段：数据增强（11×11邻域NASADEM） ---\n');

enhanced_txt_path = fullfile(config.paths.output_dir, config.files.enhanced_txt);

if exist(enhanced_txt_path, 'file') && ~config.options.save_intermediate
    fprintf('检测到已存在的增强数据文件，跳过增强步骤\n');
    fprintf('文件路径: %s\n\n', enhanced_txt_path);
else
    try
        enhanced_data = enhance_cryosat2_data(track_txt_path, config);
        
        fprintf('✓ NASADEM高程添加完成（11×11邻域）\n');
        fprintf('✓ 坡度信息添加完成\n');
        fprintf('✓ 区域标识完成\n');
        fprintf('✓ 高程变化计算完成\n');
        
        % 保存
        if config.options.save_intermediate
            dlmwrite(enhanced_txt_path, enhanced_data, 'precision', 10);
            fprintf('✓ 增强数据已保存至: %s\n', enhanced_txt_path);
        end
        
    catch ME
        fprintf('❌ 数据增强失败: %s\n', ME.message);
        fprintf('错误堆栈:\n');
        disp(ME.stack);
        return;
    end
end

fprintf('\n');

%% 第三阶段：统计分析
%
% 多时间尺度、多空间尺度统计分析

fprintf('--- 第三阶段：统计分析 ---\n');

try
    results = compute_cryosat2_statistics(enhanced_txt_path, config);
    
    fprintf('✓ 年度统计分析完成\n');
    fprintf('✓ 多年期统计完成（3年窗口）\n');
    fprintf('✓ 月度统计完成（152个月）\n');
    fprintf('✓ 网格统计完成（3年和1年窗口）\n');
    fprintf('✓ 高程带分析完成（基础+季节）\n');
    fprintf('✓ 网格高程带分析完成\n');
    fprintf('✓ 面积加权计算完成\n');
    
    % 保存结果
    results_file = fullfile(config.paths.output_dir, config.files.results_mat);
    save(results_file, 'results', '-v7.3');
    fprintf('✓ 分析结果已保存至: %s\n', results_file);
    
catch ME
    fprintf('❌ 统计分析失败: %s\n', ME.message);
    fprintf('错误堆栈:\n');
    disp(ME.stack);
    return;
end

fprintf('\n');
%% 季节结果保存为excel
gtn = shaperead(config.files.gtn_regions);
hma22 = shaperead(config.files.hma22_regions);
header = {'days','HMA', gtn.fullname, hma22.Name, 'TP'};
for i = 1:496
    header{length(header)+1} = i;
end

season_file = fullfile(config.paths.output_dir, 'CryoSat2_seasonal_elevation_results.xlsx');

numeric_data = [results.area_weighted.mid_days, results.area_weighted.seasonal_values];
data_as_cell = num2cell(numeric_data);
output_data = [header; data_as_cell];
writecell(output_data, season_file, 'Sheet', 'month_value');

numeric_data = [results.area_weighted.mid_days, results.area_weighted.seasonal_counts];
data_as_cell = num2cell(numeric_data);
output_data = [header; data_as_cell];
writecell(output_data, season_file, 'Sheet', 'month_count');

numeric_data = [results.area_weighted.mid_days, results.area_weighted.seasonal_errors];
data_as_cell = num2cell(numeric_data);
output_data = [header; data_as_cell];
writecell(output_data, season_file, 'Sheet', 'month_uncert');
%% 年结果保存为excel
header = {'start_year', 'end_year(not included)', 'HMA', gtn.fullname, hma22.Name, 'TP'};
for i = 1:496
    header{length(header)+1} = i;
end

year_file = fullfile(config.paths.output_dir, 'CryoSat2_yearly_elevation_results.xlsx');

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

fprintf('=== CryoSat-2 数据处理流程完成 ===\n\n');

% 显示处理摘要
display_processing_summary(results, config);

fprintf('\n所有输出文件已保存在: %s\n', config.paths.output_dir);
fprintf('处理流程执行完毕。\n');

%% 清理
fprintf('\n清理临时变量...\n');
clear track_data enhanced_data;
fprintf('流程结束。\n');

%% 辅助函数

function display_processing_summary(results, config)
% DISPLAY_PROCESSING_SUMMARY 显示处理摘要

fprintf('--- 处理摘要 ---\n');

% 年度统计
if isfield(results, 'annual') && ~isempty(results.annual)
    fprintf('年度统计:\n');
    fprintf('  处理年份: %d - %d\n', ...
        config.temporal.start_year, config.temporal.end_year);
    
    % 计算趋势
    valid = ~isnan(results.annual.values(:,2));
    if sum(valid) > 1
        years = results.annual.values(valid, 1);
        dh = results.annual.values(valid, 2);
        p = polyfit(years, dh, 1);
        fprintf('  整体趋势: %.3f m/yr\n', p(1));
        fprintf('  平均数据点数: %.0f 点/年\n', mean(results.annual.values(valid, 6)));
    end
end

% 多年期统计
if isfield(results, 'multiyear')
    fprintf('\n多年期统计:\n');
    fprintf('  3年窗口数: %d (2010-2025)\n', length(results.multiyear.years));
end

% 月度统计
if isfield(results, 'monthly')
    fprintf('\n月度统计:\n');
    fprintf('  月度数据点数: %d\n', size(results.monthly.values, 1));
end

% 季节高程带
if isfield(results, 'elevation_bands')
    fprintf('\n季节高程带统计:\n');
    fprintf('  高程带数: %d\n', length(results.elevation_bands.elevation_centers));
    fprintf('  时间段数: %d (3个月窗口)\n', size(results.elevation_bands.values, 3));
end

% 面积加权
if isfield(results, 'area_weighted') && ~isempty(results.area_weighted)
    fprintf('\n面积加权:\n');
    if isfield(results.area_weighted, 'seasonal_values')
        fprintf('  季节面积加权: %d时段 × %d区域\n', ...
            size(results.area_weighted.seasonal_values, 1), ...
            size(results.area_weighted.seasonal_values, 2));
    end
    if isfield(results.area_weighted, 'yearly_values')
        fprintf('  年面积加权: %d时段 × %d区域\n', ...
            size(results.area_weighted.yearly_values, 1), ...
            size(results.area_weighted.yearly_values, 2));
    end
end

% 元数据
if isfield(results, 'metadata')
    fprintf('\n处理信息:\n');
    fprintf('  处理时间: %s\n', results.metadata.processing_date);
    fprintf('  数据点数: %d\n', results.metadata.data_points);
end

fprintf('\n');

end
