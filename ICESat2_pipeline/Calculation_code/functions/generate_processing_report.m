function report_file = generate_processing_report(results, config)
% GENERATE_PROCESSING_REPORT 生成处理报告
%
% 功能：
%   生成详细的ICESat-2数据处理报告，包括：
%   1. 处理参数总结
%   2. 数据质量统计
%   3. 主要分析结果
%   4. 图表引用
%
% 输入:
%   results - 分析结果结构体
%   config - 配置参数
%
% 输出:
%   report_file - 报告文件路径

fprintf('生成处理报告...\n');

% 报告文件路径
report_file = fullfile(config.paths.output_dir, 'ICESat2_Processing_Report.txt');

% 打开文件写入
fid = fopen(report_file, 'w');
if fid == -1
    error('无法创建报告文件: %s', report_file);
end

try
    % 写入报告头部
    write_report_header(fid, config);
    
    % 写入处理配置
    write_processing_config(fid, config);
    
    % 写入数据质量统计
    write_quality_statistics(fid, results);
    
    % 写入主要分析结果
    write_analysis_results(fid, results);
    
    % 写入结论和建议
    write_conclusions(fid, results);
    
    fclose(fid);
    
    fprintf('✓ 处理报告已生成: %s\n', report_file);
    
catch ME
    fclose(fid);
    error('报告生成失败: %s', ME.message);
end

end

function write_report_header(fid, config)
% WRITE_REPORT_HEADER 写入报告头部
%
% 输入:
%   fid - 文件句柄
%   config - 配置参数

fprintf(fid, '=====================================\n');
fprintf(fid, '     ICESat-2数据处理报告\n');
fprintf(fid, '=====================================\n\n');

fprintf(fid, '生成时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '处理版本: 重构版本 v1.0\n');
fprintf(fid, '处理区域: 高亚洲山地 (HMA)\n\n');

end

function write_processing_config(fid, config)
% WRITE_PROCESSING_CONFIG 写入处理配置信息
%
% 输入:
%   fid - 文件句柄
%   config - 配置参数

fprintf(fid, '1. 处理配置参数\n');
fprintf(fid, '================\n\n');

fprintf(fid, '时间范围:\n');
fprintf(fid, '  参考日期: %s\n', datestr(config.time.reference_date, 'yyyy-mm-dd'));
fprintf(fid, '  开始日期: %s\n', datestr(config.time.start_date, 'yyyy-mm-dd'));
fprintf(fid, '  结束日期: %s\n', datestr(config.time.end_date, 'yyyy-mm-dd'));

fprintf(fid, '\n质量控制参数:\n');
fprintf(fid, '  最大高程阈值: %d m\n', config.quality.max_elevation);
fprintf(fid, '  最大高程变化阈值: %d m\n', config.quality.max_height_change);
fprintf(fid, '  最大坡度阈值: %d 度\n', config.quality.max_slope);
fprintf(fid, '  高程变化过滤阈值1: %d m\n', config.quality.height_change_filter1);
fprintf(fid, '  高程变化过滤阈值2: %d m\n', config.quality.height_change_filter2);

fprintf(fid, '\n空间参数:\n');
fprintf(fid, '  网格分辨率: %.1f°, %.1f°\n', config.spatial.grid_resolution);
fprintf(fid, '  高程带宽度: %d m, %d m\n', config.spatial.elevation_band_width);
fprintf(fid, '  高程范围: %d - %d m\n', config.spatial.elevation_range);

fprintf(fid, '\n');

end

function write_quality_statistics(fid, results)
% WRITE_QUALITY_STATISTICS 写入数据质量统计
%
% 输入:
%   fid - 文件句柄
%   results - 分析结果

fprintf(fid, '2. 数据质量统计\n');
fprintf(fid, '================\n\n');

if isfield(results, 'metadata') && isfield(results.metadata, 'quality_stats')
    quality = results.metadata.quality_stats;
    
    fprintf(fid, '数据点统计:\n');
    fprintf(fid, '  原始数据点数: %d\n', quality.original_points);
    
    if isfield(quality, 'after_height_filter')
        fprintf(fid, '  高程过滤后: %d\n', quality.after_height_filter);
        fprintf(fid, '  过滤移除点数: %d\n', quality.height_filtered);
        fprintf(fid, '  数据保留率: %.1f%%\n', ...
            100 * quality.after_height_filter / quality.original_points);
    end
    
    if isfield(quality, 'percentile_bounds')
        fprintf(fid, '  1-99百分位数范围: [%.2f, %.2f] m\n', quality.percentile_bounds);
    end
else
    fprintf(fid, '质量统计信息不可用。\n');
end

fprintf(fid, '\n');

end

function write_analysis_results(fid, results)
% WRITE_ANALYSIS_RESULTS 写入主要分析结果
%
% 输入:
%   fid - 文件句柄
%   results - 分析结果

fprintf(fid, '3. 主要分析结果\n');
fprintf(fid, '================\n\n');

% 月度分析结果
if isfield(results, 'monthly')
    fprintf(fid, '月度分析:\n');
    monthly = results.monthly;
    if ~isempty(monthly.values)
        hma_values = monthly.values(:, 1);
        valid_hma = (~isnan(hma_values)) & (hma_values ~= 0);
        
        if any(valid_hma)
            fprintf(fid, '  HMA整体月度变化范围: [%.3f, %.3f] m\n', ...
                min(hma_values(valid_hma)), max(hma_values(valid_hma)));
            fprintf(fid, '  HMA整体平均变化: %.3f m\n', mean(hma_values(valid_hma)));
            fprintf(fid, '  有效月度观测: %d/%d 个月\n', sum(valid_hma), length(hma_values));
        end
    end
    fprintf(fid, '\n');
end

% 年度分析结果
if isfield(results, 'annual')
    fprintf(fid, '年度分析:\n');
    annual = results.annual;
    if ~isempty(annual.values) && size(annual.values, 1) >= 6
        % 显示单年度结果（前6行）
        for i = 1:min(6, size(annual.year_labels, 1))
            year = annual.year_labels(i, 1);
            year2 = annual.year_labels(i, 2);
            value = annual.values(i, 1);  % HMA整体
            if (~isnan(value)) & (value ~= 0)
                fprintf(fid, '  %d-%d年HMA平均变化: %.3f m\n', year, year2, value);
            end
        end
    end
    fprintf(fid, '\n');
end

% 高程带分析结果
if isfield(results, 'elevation_bands')
    fprintf(fid, '高程带分析:\n');
    elevation = results.elevation_bands;
    
    if ~isempty(elevation.values)
        % 计算高程带平均
        dh_values = elevation.values(:, 1, :);
        dh_values(dh_values==0) = nan;
        hma_elevation_profile = mean(dh_values, 3, 'omitnan');
        valid_elevations = ~isnan(hma_elevation_profile) & hma_elevation_profile ~= 0;
        
        if any(valid_elevations)
            fprintf(fid, '  有效高程带数量: %d/%d\n', sum(valid_elevations), length(valid_elevations));
            fprintf(fid, '  高程变化范围: [%.3f, %.3f] m\n', ...
                min(hma_elevation_profile(valid_elevations)), ...
                max(hma_elevation_profile(valid_elevations)));
        end
    end
    fprintf(fid, '\n');
end

end

function write_conclusions(fid, results)
% WRITE_CONCLUSIONS 写入结论和建议
%
% 输入:
%   fid - 文件句柄
%   results - 分析结果

fprintf(fid, '4. 处理结论和建议\n');
fprintf(fid, '==================\n\n');

fprintf(fid, '处理质量评估:\n');

% 基于数据点数量评估
if isfield(results, 'metadata') && isfield(results.metadata, 'data_points')
    data_points = results.metadata.data_points;
    if data_points > 1000000
        fprintf(fid, '  ✓ 数据量充足 (%d 个观测点)\n', data_points);
    elseif data_points > 100000
        fprintf(fid, '  △ 数据量中等 (%d 个观测点)\n', data_points);
    else
        fprintf(fid, '  ⚠ 数据量较少 (%d 个观测点)\n', data_points);
    end
end

fprintf(fid, '\n建议和注意事项:\n');
fprintf(fid, '  1. 结果解释时需考虑ICESat-2激光穿透深度的影响\n');
fprintf(fid, '  2. 高坡度区域的结果不确定性较大\n');
fprintf(fid, '  3. 季节性变化可能受积雪影响\n');
fprintf(fid, '  4. 建议结合其他卫星数据进行交叉验证\n');

fprintf(fid, '\n处理完成。\n');

end
