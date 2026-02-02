function enhanced_data = enhance_icesat_data(txt_file, config)
% ENHANCE_ICESAT_DATA Enhance ICESat data by adding terrain information
%
% Key difference: ICESat uses a 3×3 neighborhood mean of NASADEM (vs CryoSat-2's 11×11)
% Reference: HMA_footprint_track_20230213.m, lines 3–28
%
% Functions:
%   1. Extract terrain parameters (3×3 neighborhood): NASADEM, slope, aspect
%   2. Compute h_li2000 (WGS84 elevation)
%   3. Add PDD correction
%   4. Compute elevation change (dh and dh_with_pdd)
%   5. Add HMA4, HMA6, and TP region identifiers
%
% ICESat data column definitions (input has 14 columns, note lon is before lat):
% Col 1: lon (longitude)
% Col 2: lat (latitude)
% Col 3: elevation (observed elevation)
% Col 4: beam_azimuth (beam azimuth)
% Col 5: days (days since 2000-01-01)
% Col 6: geoid (EGM96 geoid height)
% Col 7: satcorr (satellite elevation correction)
% Col 8: satflag (satellite correction flag)
% Col 9: glacier_flag (glacier flag, 0 = inside glacier)
% Col10: year
% Col11: hma_region (HMA region ID)
% Col12: hma22_region (HMA22 region ID)
% Col13: gtn_region (GTN15 region ID)
% Col14: grid_region (grid496 region ID)
%
% Output (expanded to 24 columns):
% Col15: nasadem (3×3 neighborhood mean) ⭐
% Col16: slope (3×3 neighborhood mean)
% Col17: aspect (3×3 neighborhood mean)
% Col18: pdd_corr (PDD correction value)
% Col19: h_li2000 (WGS84 elevation, corrected)
% Col20: dh = h_li2000 - nasadem ⭐ used for statistics
% Col21: dh_with_pdd = h_li2000 - nasadem - pdd
% Col22: HMA4 region ID
% Col23: HMA6 region ID (seasonal pattern)
% Col24: TP flag

fprintf('开始ICESat数据增强处理...\n');
warning('off', 'MATLAB:polyshape:repairedBySimplify');
warning('off', 'map:removing:latlon2pix');

%% Load data
fprintf('加载数据文件...\n');
try
    data = load(txt_file);
    fprintf('✓ 成功加载 %d 个数据点\n', size(data, 1));
catch ME
    error('数据文件加载失败: %s', ME.message);
end

%% Load terrain data
fprintf('加载地形数据...\n');
try
    [nasadem, ref_nasadem] = geotiffread(config.files.nasadem);
    [slope_data, ref_slope] = geotiffread(config.files.slope);
    [aspect_data, ref_aspect] = geotiffread(config.files.aspect);
    
    % Load region data
    glacier_grid = shaperead(config.files.grid10);
    HMA4 = shaperead(config.files.hma4_regions);
    HMA6 = shaperead(config.files.hma6_regions);
    TP = shaperead(config.files.tp_boundary);
    
    fprintf('✓ 地形数据加载完成\n');
catch ME
    error('地形数据加载失败: %s', ME.message);
end

%% Add terrain parameters (3×3 neighborhood mean)
fprintf('添加地形参数（3×3邻域）：NASADEM、slope、aspect...\n');
enhanced_data = add_terrain_parameters_3x3(data, nasadem, ref_nasadem, slope_data, ref_slope, aspect_data, ref_aspect, config);

%% Compute h_li2000 and elevation change
fprintf('计算h_li2000和高程变化...\n');
enhanced_data = calculate_h_li2000_and_dh(enhanced_data);

%% Add PDD correction and compute dh_with_pdd
fprintf('添加pdd校正并计算dh_with_pdd...\n');
enhanced_data = add_pdd_and_calculate_dh_with_pdd(enhanced_data, glacier_grid);

%% Add region identifiers
fprintf('添加区域标识（HMA4, HMA6, TP）...\n');
enhanced_data = add_region_identifiers(enhanced_data, HMA4, HMA6, TP);

fprintf('✓ 数据增强完成，输出 %d 个数据点\n', size(enhanced_data, 1));

end

function enhanced_data = add_terrain_parameters_3x3(data, nasadem, ref_nasadem, slope_data, ref_slope, aspect_data, ref_aspect, config)
% ADD_TERRAIN_PARAMETERS_3X3 Add terrain parameters (3×3 neighborhood mean)
%
% This is the key difference between ICESat and CryoSat-2!
% ICESat:    3×3 neighborhood (9 pixels)
% CryoSat-2: 11×11 neighborhood (121 pixels)
%
% Reason: ICESat footprint size is ~70 m, smaller than CryoSat-2's ~300 m
%
% Output:
% Col15: nasadem (3×3 neighborhood)
% Col16: slope (3×3 neighborhood)
% Col17: aspect (3×3 neighborhood)

neighbor_size = config.spatial.nasadem_neighbor_size;  % 3
half_size = floor(neighbor_size / 2);  % 1

% Extend data matrix
enhanced_data = data;
if size(enhanced_data, 2) < 17
    enhanced_data(:, 15) = nan;  % nasadem
    enhanced_data(:, 16) = nan;  % slope
    enhanced_data(:, 17) = nan;  % aspect
end

fprintf('  使用%d×%d邻域平均（%d个像元）...\n', ...
    neighbor_size, neighbor_size, neighbor_size^2);

for i = 1:size(data, 1)
    % Get pixel coordinates (note: lon is in column 1, lat is in column 2)
    [row, col] = latlon2pix(ref_nasadem, data(i, 2), data(i, 1));
    row = ceil(row);
    col = ceil(col);
    
    % Check boundaries
    if row <= half_size || row > size(nasadem, 1) - half_size || ...
       col <= half_size || col > size(nasadem, 2) - half_size
        continue;
    end
    
    % 3×3 neighborhood mean - NASADEM
    ele_sum = [];
    slope_sum = [];
    aspect_sum = [];
    
    for j = -half_size:half_size  % -1:1
        for k = -half_size:half_size
            ele_sum = [ele_sum, single(nasadem(row+j, col+k))];
            slope_sum = [slope_sum, single(slope_data(row+j, col+k))];
            aspect_sum = [aspect_sum, single(aspect_data(row+j, col+k))];
        end
    end
    
    enhanced_data(i, 15) = nanmean(ele_sum);      % nasadem
    enhanced_data(i, 16) = nanmean(slope_sum);    % slope
    enhanced_data(i, 17) = nanmean(aspect_sum);   % aspect
    
    % Show progress
    if mod(i, 100000) == 0
        fprintf('  处理进度: %d/%d (%.1f%%)\n', i, size(data, 1), 100*i/size(data,1));
    end
end

fprintf('✓ 地形参数（3×3邻域）添加完成\n');

end

function enhanced_data = calculate_h_li2000_and_dh(data)
% CALCULATE_H_LI2000_AND_DH Compute h_li2000 and elevation change
%
% Output:
% Col19: h_li2000 (WGS84 elevation)
% Col20: dh = h_li2000 - nasadem

enhanced_data = data;
if size(enhanced_data, 2) < 20
    enhanced_data(:, 19) = nan;  % h_li2000
    enhanced_data(:, 20) = nan;  % dh
end

for i = 1:size(data, 1)
    % Compute h_li2000 (see HMA_Icesattrack_full.m)
    latitude = enhanced_data(i, 2);
    offset = 0.8337-0.3848*(sind(latitude))^2;
    if data(i, 8) == 2  % satflag
        enhanced_data(i, 19) = data(i, 3) - data(i, 6) + data(i, 7) - offset;
    else
        enhanced_data(i, 19) = data(i, 3) - data(i, 6) - offset;
    end
    
    % Compute dh = h_li2000 - nasadem, without considering PDD
    enhanced_data(i, 20) = enhanced_data(i, 19) - data(i, 15);
end

fprintf('✓ h_li2000和dh计算完成\n');

end

function enhanced_data = add_pdd_and_calculate_dh_with_pdd(data, glacier_grid)
% ADD_PDD_AND_CALCULATE_DH_WITH_PDD Add PDD correction and compute dh_with_pdd
%
% Output:
% Col18: pdd_corr (obtained from grid_region)
% Col21: dh_with_pdd = h_li2000 - nasadem - pdd

enhanced_data = data;
if size(enhanced_data, 2) < 21
    enhanced_data(:, 18) = nan;  % pdd_corr
    enhanced_data(:, 21) = nan;  % dh_with_pdd
end

% Obtain PDD correction value from column 14 (grid_region)
for i = 1:size(data, 1)
    grid_id = data(i, 14);
    if ~isnan(grid_id) && grid_id > 0 && grid_id <= length(glacier_grid)
        enhanced_data(i, 18) = glacier_grid(grid_id).pdd;
    end
end

% Compute dh_with_pdd = h_li2000 - nasadem - pdd
enhanced_data(:, 21) = data(:, 19) - data(:, 15) - enhanced_data(:, 18);

fprintf('✓ pdd校正和dh_with_pdd计算完成\n');

end

function enhanced_data = add_region_identifiers(data, HMA4, HMA6, TP)
% ADD_REGION_IDENTIFIERS Add region identifiers
% Col22: HMA4 region ID
% Col23: HMA6 region ID (seasonal pattern)
% Col24: TP flag

enhanced_data = data;
if size(enhanced_data, 2) < 24
    enhanced_data(:, 22) = 0;
    enhanced_data(:, 23) = 0;
    enhanced_data(:, 24) = 0;
end

% Add HMA4 region identifiers (note: lon is in column 1, lat is in column 2)
for k = 1:length(HMA4)
    [in, ~] = inpolygon(data(:,1), data(:,2), HMA4(k).X, HMA4(k).Y);
    enhanced_data(in==1, 22) = k;
end

% Add HMA6 seasonal-pattern region identifiers
for k = 1:length(HMA6)
    [in, ~] = inpolygon(data(:,1), data(:,2), HMA6(k).X, HMA6(k).Y);
    enhanced_data(in==1, 23) = k;
end

% Add TP flag
[in, ~] = inpolygon(data(:,1), data(:,2), TP.X, TP.Y);
enhanced_data(in==1, 24) = 1;

fprintf('✓ 区域标识添加完成\n');

end

