function varargout = utils(varargin)
% UTILS Utility function collection for ICESat processing
%
% Reference: revised_CryoSat2/functions/utils.m

if nargin < 1
    error('需要指定工具函数名称');
end

func_name = varargin{1};

switch lower(func_name)
    case 'customcolormap'
        if nargin < 3
            error('customcolormap需要positions和colors参数');
        end
        positions = varargin{2};
        colors = varargin{3};
        varargout{1} = create_custom_colormap(positions, colors);
        
    case 'robust_mean'
        if nargin < 2
            error('robust_mean需要data参数');
        end
        data = varargin{2};
        percentile = 1;  % ICESat uses 1% percentile (vs CryoSat2 using 5%)
        if nargin >= 3
            percentile = varargin{3};
        end
        [varargout{1}, varargout{2}] = compute_robust_mean(data, percentile);
        
    case 'median_filter'
        if nargin < 2
            error('median_filter需要data参数');
        end
        data = varargin{2};
        threshold = 50;  % ICESat uses a 50 m threshold (vs 75 m for CryoSat2)
        if nargin >= 3
            threshold = varargin{3};
        end
        varargout{1} = median_filter(data, threshold);
        
    otherwise
        error('未知的工具函数: %s', func_name);
end

end

function cmap = create_custom_colormap(positions, colors)
% CREATE_CUSTOM_COLORMAP Create a custom color map

if length(positions) ~= size(colors, 1)
    error('位置向量长度必须与颜色矩阵行数相等');
end

map_size = 256;
cmap = zeros(map_size, 3);

for i = 1:3
    cmap(:, i) = interp1(positions, colors(:, i), linspace(0, 1, map_size));
end

cmap = max(0, min(1, cmap));

end

function [robust_mean, robust_std] = compute_robust_mean(data, percentile)
% COMPUTE_ROBUST_MEAN Compute a robust mean

if isempty(data) || all(isnan(data))
    robust_mean = NaN;
    robust_std = NaN;
    return;
end

valid_data = data(~isnan(data));

if isempty(valid_data)
    robust_mean = NaN;
    robust_std = NaN;
    return;
end

lower_bound = prctile(valid_data, percentile);
upper_bound = prctile(valid_data, 100 - percentile);

filtered_data = valid_data(valid_data >= lower_bound & valid_data <= upper_bound);

if isempty(filtered_data)
    robust_mean = NaN;
    robust_std = NaN;
else
    robust_mean = mean(filtered_data);
    robust_std = std(filtered_data);
end

end

function filtered_data = median_filter(data, threshold)
% MEDIAN_FILTER Outlier filtering based on the median

if isempty(data)
    filtered_data = data;
    return;
end

med = median(data, 'omitnan');
valid_mask = abs(data - med) < threshold;
filtered_data = data(valid_mask);

end

