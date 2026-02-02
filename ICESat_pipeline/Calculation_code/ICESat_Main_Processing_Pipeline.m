--- START OF FILE ICESat_Main_Processing_Pipeline.m ---

%% ICESat Data Processing Main Pipeline Script
%
% Reference: HMA_ICESat_COMPUTE_20230603.m + revised_CryoSat2/CryoSat2_Main_Processing_Pipeline.m
%
% Function: Complete processing pipeline for High Mountain Asia (HMA) ICESat data
% Refactored Version: v1.0
% Data Period: 2003-2009
%
% Key Difference: ICESat uses 3x3 neighborhood average NASADEM (vs CryoSat-2's 11x11)
%
% Processing Steps:
% 1. Raw Data Extraction (ICESat binary -> TXT) - Manual run or existing data required
% 2. Data Enhancement (3x3 neighborhood NASADEM + Region ID)
% 3. Statistical Analysis (Annual / 3-year window / Campaign)
% 4. Result Output

%% Initialization
clear; clc; close all;

% Add function paths
addpath('functions');

% Load configuration
config = load_icesat_config();

fprintf('=== ICESat Data Processing Pipeline Started ===\n');
fprintf('Important Note: ICESat uses 3x3 neighborhood average NASADEM\n');
fprintf('Data Time Range: 2003-2009\n');
fprintf('Configuration loaded, starting data processing...\n\n');

%% Phase 1: Raw Data Extraction
%
% Note: ICESat raw data extraction usually requires extraction from GLA14 binary files.
% If extracted TXT files already exist, this step can be skipped.
%
% Assuming existing file: HMA_ICESat_rgiregion.txt

fprintf('--- Phase 1: Raw Data Extraction ---\n');

track_txt_path = fullfile(config.paths.output_dir, config.files.track_txt);

if ~exist(track_txt_path, 'file')
    warning('Raw data file not found: %s\n', track_txt_path);
    [track_data, track_data_outrgi] = extract_icesat_tracks(config);
    dlmwrite(track_txt_path, track_data, 'precision', 10);

    dlmwrite(fullfile(config.paths.output_dir, 'HMA_ICESat_norgiregion.txt'), track_data_outrgi, 'precision', 10);

    fprintf('ICESat glacier data extraction completed to %s\n\n', track_txt_path);
else
    fprintf('✓ Raw data file found\n');
    fprintf('File path: %s\n\n', track_txt_path);
end

%% Phase 2: Data Enhancement
%
% Key Step: Use 3x3 neighborhood average NASADEM (Specific to ICESat)
% Add grid, slope, region IDs, calculate elevation change

fprintf('--- Phase 2: Data Enhancement (3x3 NASADEM) ---\n');

enhanced_txt_path = fullfile(config.paths.output_dir, config.files.enhanced_txt);

if exist(enhanced_txt_path, 'file') && ~config.options.save_intermediate
    fprintf('Existing enhanced data file detected, skipping enhancement step\n');
    fprintf('File path: %s\n\n', enhanced_txt_path);
else
    enhanced_data = enhance_icesat_data(track_txt_path, config);
    
    fprintf('✓ Terrain parameters added (3x3 neighborhood): NASADEM, slope, aspect\n');
    fprintf('✓ h_li2000 calculation completed\n');
    fprintf('✓ PDD correction added\n');
    fprintf('✓ Elevation change calculation completed (dh and dh_with_pdd)\n');
    fprintf('✓ Region identification completed (HMA4, HMA6, TP)\n');
    
    % Save
    if config.options.save_intermediate
        dlmwrite(enhanced_txt_path, enhanced_data, 'precision', 10);
        fprintf('✓ Enhanced data saved to: %s\n', enhanced_txt_path);
    end
end

fprintf('\n');

%% Phase 3: Statistical Analysis
%
% Multi-temporal scale statistical analysis:
% - Annual statistics (2003-2009, each year to next October)
% - 3-year sliding window
% - Campaign statistics (19 observation windows)

fprintf('--- Phase 3: Statistical Analysis ---\n');

results = compute_icesat_statistics(enhanced_txt_path, config);

fprintf('✓ Annual statistical analysis completed (2003-2009)\n');
fprintf('✓ 3-year sliding window statistics completed\n');
fprintf('✓ Campaign statistics completed (19 observation windows)\n');

% Save results
results_file = fullfile(config.paths.output_dir, config.files.results_mat);
save(results_file, 'results', '-v7.3');
fprintf('✓ Analysis results saved to: %s\n', results_file);

fprintf('\n');

%% Phase 4: Result Export to Excel
%
% Export statistical results to Excel files for further analysis

fprintf('--- Phase 4: Result Export ---\n');

%% Save seasonal results to Excel
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
%% Save yearly results to Excel
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

%% Output elevation band results for each region: ele_bands * year
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

%% Processing Completion Summary

fprintf('=== ICESat Data Processing Pipeline Completed ===\n\n');

% Display processing summary
display_processing_summary(results, config);

fprintf('\nAll output files saved in: %s\n', config.paths.output_dir);
fprintf('Processing pipeline finished.\n');

%% Cleanup
fprintf('\nCleaning up temporary variables...\n');
clear enhanced_data;
fprintf('Pipeline ended.\n');

%% Helper Functions

function display_processing_summary(results, config)
% DISPLAY_PROCESSING_SUMMARY Displays processing summary

fprintf('--- Processing Summary ---\n');

% Data Information
if isfield(results, 'metadata')
    fprintf('Data Information:\n');
    fprintf('  Original data points: %d\n', results.metadata.data_points_original);
    fprintf('  After Quality Control: %d\n', results.metadata.data_points_qc);
    fprintf('  Data retention rate: %.1f%%\n', ...
        100 * results.metadata.data_points_qc / results.metadata.data_points_original);
end

% Annual Statistics
if isfield(results, 'annual')
    fprintf('\nAnnual Statistics:\n');
    fprintf('  Processing years: %d - %d\n', config.temporal.start_year, config.temporal.end_year);
    
    % HMA Overall Trend
    valid = ~isnan(results.annual.values(:,1));
    if sum(valid) > 1
        years = results.annual.year_label;
        dh = results.annual.values(valid, 1);
        p = polyfit(years, dh, 1);
        fprintf('  HMA Overall Trend: %.3f m/yr\n', p(1));
        fprintf('  Average data points: %.0f points/year\n', mean(results.annual.counts(valid, 1)));
    end
end

% 3-Year Window Statistics
if isfield(results, 'multiyear')
    fprintf('\n3-Year Window Statistics:\n');
    fprintf('  Number of windows: %d\n', size(results.multiyear.values, 1));
end

% Campaign Statistics
if isfield(results, 'campaigns')
    fprintf('\nCampaign Statistics:\n');
    fprintf('  Number of observation windows: %d\n', size(results.campaigns.values, 1));
end

% Processing Time
if isfield(results, 'metadata')
    fprintf('\nProcessing Information:\n');
    fprintf('  Processing Date: %s\n', results.metadata.processing_date);
end

fprintf('\n');

end