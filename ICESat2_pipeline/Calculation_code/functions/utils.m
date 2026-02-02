function varargout = utils(varargin)
% UTILS ICESat-2处理实用工具函数集合
%
% 使用方法:
%   colormap_data = utils('customcolormap', [0 0.4 1], [0 0 1; 1 1 1; 1 0 0]);
%   valid_data = utils('remove_outliers', data, threshold);

if nargin < 1
    error('需要指定工具函数名称');
end

func_name = varargin{1};

switch lower(func_name)
    case 'customcolormap'
        % 自定义颜色映射
        if nargin < 3
            error('customcolormap需要positions和colors参数');
        end
        positions = varargin{2};
        colors = varargin{3};
        varargout{1} = create_custom_colormap(positions, colors);
        
    case 'remove_outliers'
        % 移除异常值
        if nargin < 3
            error('remove_outliers需要data和threshold参数');
        end
        data = varargin{2};
        threshold = varargin{3};
        varargout{1} = remove_outliers(data, threshold);
        
    case 'robust_mean'
        % 鲁棒平均值计算
        if nargin < 2
            error('robust_mean需要data参数');
        end
        data = varargin{2};
        if nargin >= 3
            percentile = varargin{3};
        else
            percentile = 5;  % 默认移除5%异常值
        end
        [varargout{1}, varargout{2}] = compute_robust_mean(data, percentile);
        
    case 'date_to_days'
        % 日期转换为相对天数
        if nargin < 3
            error('date_to_days需要date和reference_date参数');
        end
        date = varargin{2};
        ref_date = varargin{3};
        varargout{1} = daysact(ref_date, date);
        
    case 'days_to_date'
        % 相对天数转换为日期
        if nargin < 3
            error('days_to_date需要days和reference_date参数');
        end
        days = varargin{2};
        ref_date = varargin{3};
        varargout{1} = ref_date + days(days);
        
    otherwise
        error('未知的工具函数: %s', func_name);
end

end

function cmap = create_custom_colormap(positions, colors)
% CREATE_CUSTOM_COLORMAP 创建自定义颜色映射
%
% 输入:
%   positions - 颜色位置向量
%   colors - 颜色矩阵
%
% 输出:
%   cmap - 颜色映射矩阵

if length(positions) ~= size(colors, 1)
    error('位置向量长度必须与颜色矩阵行数相等');
end

% 创建256级颜色映射
map_size = 256;
cmap = zeros(map_size, 3);

for i = 1:3  % RGB三个通道
    cmap(:, i) = interp1(positions, colors(:, i), linspace(0, 1, map_size));
end

% 确保颜色值在[0,1]范围内
cmap = max(0, min(1, cmap));

end

function clean_data = remove_outliers(data, threshold)
% REMOVE_OUTLIERS 移除异常值
%
% 输入:
%   data - 数据向量
%   threshold - 阈值（标准差倍数）
%
% 输出:
%   clean_data - 清理后的数据

if isempty(data)
    clean_data = data;
    return;
end

% 计算统计量
data_mean = mean(data, 'omitnan');
data_std = std(data, 'omitnan');

% 标识异常值
outlier_mask = abs(data - data_mean) > threshold * data_std;

% 移除异常值
clean_data = data(~outlier_mask);

end

function [robust_mean, robust_std] = compute_robust_mean(data, percentile)
% COMPUTE_ROBUST_MEAN 计算鲁棒平均值
%
% 输入:
%   data - 数据向量
%   percentile - 要移除的百分位数（两端各移除该百分比）
%
% 输出:
%   robust_mean - 鲁棒平均值
%   robust_std - 鲁棒标准差

if isempty(data) || all(isnan(data))
    robust_mean = NaN;
    robust_std = NaN;
    return;
end

% 移除NaN值
valid_data = data(~isnan(data));

if isempty(valid_data)
    robust_mean = NaN;
    robust_std = NaN;
    return;
end

% 计算百分位数界限
lower_bound = prctile(valid_data, percentile);
upper_bound = prctile(valid_data, 100 - percentile);

% 过滤数据
filtered_data = valid_data(valid_data >= lower_bound & valid_data <= upper_bound);

if isempty(filtered_data)
    robust_mean = NaN;
    robust_std = NaN;
else
    robust_mean = mean(filtered_data);
    robust_std = std(filtered_data);
end

end
