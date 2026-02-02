%% ICESat-2数据处理主流程脚本
% 
% 功能：高亚洲山地（HMA）ICESat-2数据完整处理流程
% 作者：重构版本
% 日期：2025年
% 
% 处理流程：
% 1. 原始数据提取与轨道处理
% 2. 数据后处理与质量增强  
% 3. 统计分析与计算
% 4. 可视化与结果输出

%% 初始化和配置设置
clear; clc; close all;

% 添加函数路径
addpath('functions');

% 加载配置参数
config = load_processing_config();

fprintf('=== ICESat-2 数据处理流程开始 ===\n');
fprintf('配置加载完成，开始数据处理...\n\n');

%% 第一阶段：原始数据提取与轨道处理
% 
% 从ICESat-2 ATL06 v6 HDF5文件中提取轨道数据，进行空间过滤和区域标识

fprintf('--- 第一阶段：原始数据提取与轨道处理 ---\n');

% 检查输入数据路径
if ~exist(config.paths.input_h5_dir, 'dir')
    error('输入HDF5数据目录不存在: %s', config.paths.input_h5_dir);
end

% 执行轨道数据提取
try
    track_data = extract_icesat2_tracks(config);
    fprintf('✓ 轨道数据提取完成，共处理 %d 个数据点\n', size(track_data, 1));
    
    % 保存中间结果
    csv_output = fullfile(config.paths.output_dir, config.files.track_csv);
    writematrix(track_data, csv_output);
    fprintf('✓ 数据已保存至: %s\n\n', csv_output);
    
catch ME
    fprintf('❌ 轨道数据提取失败: %s\n', ME.message);
    return;
end

%% 第二阶段：数据后处理与质量增强
% 
% 添加地形信息、穿透深度校正，计算最终高程变化

fprintf('--- 第二阶段：数据后处理与质量增强 ---\n');

try
    % 加载第一阶段数据
    csv_file = fullfile(config.paths.output_dir, config.files.track_csv);
    enhanced_data = enhance_icesat2_data(csv_file, config);
    
    fprintf('✓ 地形信息添加完成\n');
    fprintf('✓ 穿透深度校正完成\n');
    fprintf('✓ 高程变化计算完成\n');
    
    % 保存增强后的数据
    txt_output = fullfile(config.paths.output_dir, config.files.enhanced_txt);
    writematrix(enhanced_data, txt_output);
    fprintf('✓ 增强数据已保存至: %s\n\n', txt_output);
    
catch ME
    fprintf('❌ 数据后处理失败: %s\n', ME.message);
    return;
end

%% 第三阶段：统计分析与计算
% 
% 进行多时间尺度、多空间尺度的统计分析

fprintf('--- 第三阶段：统计分析与计算 ---\n');

try
    % 加载增强后的数据
    txt_file = fullfile(config.paths.output_dir, config.files.enhanced_txt);
    
    % 执行统计分析
    results = compute_icesat2_statistics(txt_file, config);
    
    fprintf('✓ 时间序列分析完成\n');
    fprintf('✓ 空间统计分析完成\n');
    fprintf('✓ 高程带分析完成\n');
    fprintf('✓ 面积加权计算完成\n');
    
    % 保存分析结果
    results_file = fullfile(config.paths.output_dir, 'ICESat2_analysis_results.mat');
    save(results_file, 'results', '-v7.3');
    fprintf('✓ 分析结果已保存至: %s\n\n', results_file);
    
catch ME
    fprintf('❌ 统计分析失败: %s\n', ME.message);
    return;
end
%% 季节结果保存为excel
gtn = shaperead(config.files.gtn_regions);
hma22 = shaperead(config.files.hma22_regions);
header = {'days','HMA', gtn.fullname, hma22.Name, 'TP'};
for i = 1:496
    header{length(header)+1} = i;
end

season_file = fullfile(config.paths.output_dir, 'ICESat2_seasonal_elevation_results.xlsx');

numeric_data = [results.area_weighted_season.mid_days, results.area_weighted_season.values];
data_as_cell = num2cell(numeric_data);
output_data = [header; data_as_cell];
writecell(output_data, season_file, 'Sheet', 'month_value');

numeric_data = [results.area_weighted_season.mid_days, results.area_weighted_season.total_counts];
data_as_cell = num2cell(numeric_data);
output_data = [header; data_as_cell];
writecell(output_data, season_file, 'Sheet', 'month_count');

numeric_data = [results.area_weighted_season.mid_days, results.area_weighted_season.weighted_errors];
data_as_cell = num2cell(numeric_data);
output_data = [header; data_as_cell];
writecell(output_data, season_file, 'Sheet', 'month_uncert');
%% 年结果保存为excel
header = {'start_year', 'end_year(not included)', 'HMA', gtn.fullname, hma22.Name, 'TP'};
for i = 1:496
    header{length(header)+1} = i;
end

year_file = fullfile(config.paths.output_dir, 'ICESat2_yearly_elevation_results.xlsx');

numeric_data = [results.area_weighted_annual.year_label, results.area_weighted_annual.values];
data_as_cell = num2cell(numeric_data);
output_data = [header; data_as_cell];
writecell(output_data, year_file, 'Sheet', 'year_value');

numeric_data = [results.area_weighted_annual.year_label, results.area_weighted_annual.total_counts];
data_as_cell = num2cell(numeric_data);
output_data = [header; data_as_cell];
writecell(output_data, year_file, 'Sheet', 'year_count');

numeric_data = [results.area_weighted_annual.year_label, results.area_weighted_annual.weighted_errors];
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

%% 可选：生成处理报告
% 
% 生成详细的处理报告文档

if config.options.generate_report
    fprintf('\n--- 生成处理报告 ---\n');
    
    try
        report_file = generate_processing_report(results, config);
        fprintf('✓ 处理报告已生成: %s\n', report_file);
    catch ME
        fprintf('❌ 报告生成失败: %s\n', ME.message);
    end
end

%% 清理和资源释放

fprintf('\n清理临时变量...\n');
clear track_data enhanced_data;
fprintf('流程结束。\n');
