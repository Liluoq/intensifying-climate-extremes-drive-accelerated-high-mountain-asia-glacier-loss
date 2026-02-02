--- START OF FILE ICESat2_Main_Processing_Pipeline.m ---

%% ICESat-2 Data Processing Main Pipeline Script
% 
% Function: Complete processing pipeline for High Mountain Asia (HMA) ICESat-2 data
% Author: Refactored Version
% Date: 2025
% 
% Processing Steps:
% 1. Raw data extraction and track processing
% 2. Data post-processing and quality enhancement  
% 3. Statistical analysis and calculation
% 4. Visualization and result output

%% Initialization and Configuration Settings
clear; clc; close all;

% Add function paths
addpath('functions');

% Load configuration parameters
config = load_processing_config();

fprintf('=== ICESat-2 Data Processing Pipeline Started ===\n');
fprintf('Configuration loaded, starting data processing...\n\n');

%% Phase 1: Raw Data Extraction and Track Processing
% 
% Extract track data from ICESat-2 ATL06 v6 HDF5 files, perform spatial filtering and region identification

fprintf('--- Phase 1: Raw Data Extraction and Track Processing ---\n');

% Check input data path
if ~exist(config.paths.input_h5_dir, 'dir')
    error('Input HDF5 data directory does not exist: %s', config.paths.input_h5_dir);
end

% Execute track data extraction
try
    track_data = extract_icesat2_tracks(config);
    fprintf('✓ Track data extraction completed, processed %d data points\n', size(track_data, 1));
    
    % Save intermediate results
    csv_output = fullfile(config.paths.output_dir, config.files.track_csv);
    writematrix(track_data, csv_output);
    fprintf('✓ Data saved to: %s\n\n', csv_output);
    
catch ME
    fprintf('❌ Track data extraction failed: %s\n', ME.message);
    return;
end

%% Phase 2: Data Post-processing and Quality Enhancement
% 
% Add terrain information, penetration depth correction, and calculate final elevation change

fprintf('--- Phase 2: Data Post-processing and Quality Enhancement ---\n');

try
    % Load Phase 1 data
    csv_file = fullfile(config.paths.output_dir, config.files.track_csv);
    enhanced_data = enhance_icesat2_data(csv_file, config);
    
    fprintf('✓ Terrain information added\n');
    fprintf('✓ Penetration depth correction completed\n');
    fprintf('✓ Elevation change calculation completed\n');
    
    % Save enhanced data
    txt_output = fullfile(config.paths.output_dir, config.files.enhanced_txt);
    writematrix(enhanced_data, txt_output);
    fprintf('✓ Enhanced data saved to: %s\n\n', txt_output);
    
catch ME
    fprintf('❌ Data post-processing failed: %s\n', ME.message);
    return;
end

%% Phase 3: Statistical Analysis and Calculation
% 
% Perform multi-temporal and multi-spatial scale statistical analysis

fprintf('--- Phase 3: Statistical Analysis and Calculation ---\n');

try
    % Load enhanced data
    txt_file = fullfile(config.paths.output_dir, config.files.enhanced_txt);
    
    % Execute statistical analysis
    results = compute_icesat2_statistics(txt_file, config);
    
    fprintf('✓ Time series analysis completed\n');
    fprintf('✓ Spatial statistical analysis completed\n');
    fprintf('✓ Elevation band analysis completed\n');
    fprintf('✓ Area-weighted calculation completed\n');
    
    % Save analysis results
    results_file = fullfile(config.paths.output_dir, 'ICESat2_analysis_results.mat');
    save(results_file, 'results', '-v7.3');
    fprintf('✓ Analysis results saved to: %s\n\n', results_file);
    
catch ME
    fprintf('❌ Statistical analysis failed: %s\n', ME.message);
    return;
end
%% Save seasonal results to Excel
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
%% Save yearly results to Excel
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
    save_path = fullfile(config.paths.ele_band_dir, ['[', num2str(i),']',region_name, '.xlsx']);
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

%% Optional: Generate Processing Report
% 
% Generate detailed processing report document

if config.options.generate_report
    fprintf('\n--- Generating Processing Report ---\n');
    
    try
        report_file = generate_processing_report(results, config);
        fprintf('✓ Processing report generated: %s\n', report_file);
    catch ME
        fprintf('❌ Report generation failed: %s\n', ME.message);
    end
end

%% Cleanup and Resource Release

fprintf('\nCleaning up temporary variables...\n');
clear track_data enhanced_data;
fprintf('Pipeline ended.\n');