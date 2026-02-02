--- START OF FILE text/plain ---

%% CryoSat-2 Data Processing Main Pipeline Script
% 
% Reference: revised_ICESat2_code/ICESat2_Main_Processing_Pipeline.m
%
% Function: Complete processing pipeline for High Mountain Asia (HMA) CryoSat-2 data
% Refactored Version: v2.0 (Correct Version)
% Reference Code: revised_ICESat2_code + HMA_footprint_track_20230213.m
% 
% Key Difference: CryoSat-2 uses 11x11 neighborhood average NASADEM (vs ICESat-2's 3x3)
%
% Processing Steps:
% 1. Raw Data Extraction (L2I/TEMPO NetCDF -> TXT)
% 2. Data Enhancement (11x11 neighborhood NASADEM + Region ID)
% 3. Statistical Analysis (Annual / Spatial / Elevation Bands)
% 4. Visualization Output

if isempty(gcp('nocreate')) % If no parallel pool is running
    parpool(12); % Create a parallel pool with 12 workers
end

%% Initialization
clear; clc; close all;

% Add function paths
addpath('functions');

% Load configuration
config = load_cryosat2_config();

fprintf('=== CryoSat-2 Data Processing Pipeline Started ===\n');
fprintf('Important Note: CryoSat-2 uses 11x11 neighborhood average NASADEM\n');
fprintf('Configuration loaded, starting data processing...\n\n');

%% Phase 1: Raw Data Extraction
%
% Extract track data from CryoSat-2 L2I/TEMPO NetCDF files
% Divided into RGI glacier regions (inside) and outside regions

fprintf('--- Phase 1: Raw Data Extraction ---\n');

% Extract RGI glacier region data
track_txt_path = fullfile(config.paths.output_dir, config.files.track_txt_tempo);

if exist(track_txt_path, 'file') && ~config.options.save_intermediate
    fprintf('Existing track data file detected, skipping extraction step\n');
    fprintf('File path: %s\n\n', track_txt_path);
else
    try
        fprintf('Extracting footprint data within RGI glacier regions...\n');
        track_data = extract_cryosat2_tracks(config);
        fprintf('✓ Track data extraction completed, total %d data points\n', size(track_data, 1));
        
        % Save
        if config.options.save_intermediate
            dlmwrite(track_txt_path, track_data, 'precision', 10); % Keep 10 significant digits
            fprintf('✓ Data saved to: %s\n', track_txt_path);
        end
        
    catch ME
        fprintf('❌ Track data extraction failed: %s\n', ME.message);
        fprintf('Error Stack:\n');
        disp(ME.stack);
        return;
    end
end

%% Optional: Extract non-RGI region data (for comparative study)
log_fid = fopen(config.files.process_log, 'w');
fprintf(log_fid, '[%s] Start to process\n', datetime('now'));

if config.options.save_intermediate
    fprintf('\nExtracting non-RGI region data (optional)...\n');
    try
        extract_nonrgi_cryosat2_tracks(config);
        fprintf('✓ Non-RGI region data extraction completed\n');
    catch ME
        fprintf('⚠ Non-RGI data extraction failed (can be ignored): %s\n', ME.message);
    end
end

fprintf('\n');

%% Phase 2: Data Enhancement
%
% Key Step: Use 11x11 neighborhood average NASADEM (Specific to CryoSat-2)
% Add slope, region IDs, calculate elevation change

fprintf('--- Phase 2: Data Enhancement (11x11 NASADEM) ---\n');

enhanced_txt_path = fullfile(config.paths.output_dir, config.files.enhanced_txt);

if exist(enhanced_txt_path, 'file') && ~config.options.save_intermediate
    fprintf('Existing enhanced data file detected, skipping enhancement step\n');
    fprintf('File path: %s\n\n', enhanced_txt_path);
else
    try
        enhanced_data = enhance_cryosat2_data(track_txt_path, config);
        
        fprintf('✓ NASADEM elevation added (11x11 neighborhood)\n');
        fprintf('✓ Slope information added\n');
        fprintf('✓ Region identification completed\n');
        fprintf('✓ Elevation change calculation completed\n');
        
        % Save
        if config.options.save_intermediate
            dlmwrite(enhanced_txt_path, enhanced_data, 'precision', 10);
            fprintf('✓ Enhanced data saved to: %s\n', enhanced_txt_path);
        end
        
    catch ME
        fprintf('❌ Data enhancement failed: %s\n', ME.message);
        fprintf('Error Stack:\n');
        disp(ME.stack);
        return;
    end
end

fprintf('\n');

%% Phase 3: Statistical Analysis
%
% Multi-temporal and multi-spatial scale statistical analysis

fprintf('--- Phase 3: Statistical Analysis ---\n');

try
    results = compute_cryosat2_statistics(enhanced_txt_path, config);
    
    fprintf('✓ Annual statistical analysis completed\n');
    fprintf('✓ Multi-year statistics completed (3-year window)\n');
    fprintf('✓ Monthly statistics completed (152 months)\n');
    fprintf('✓ Grid statistics completed (3-year and 1-year windows)\n');
    fprintf('✓ Elevation band analysis completed (Basic + Seasonal)\n');
    fprintf('✓ Grid elevation band analysis completed\n');
    fprintf('✓ Area-weighted calculation completed\n');
    
    % Save results
    results_file = fullfile(config.paths.output_dir, config.files.results_mat);
    save(results_file, 'results', '-v7.3');
    fprintf('✓ Analysis results saved to: %s\n', results_file);
    
catch ME
    fprintf('❌ Statistical analysis failed: %s\n', ME.message);
    fprintf('Error Stack:\n');
    disp(ME.stack);
    return;
end

fprintf('\n');
%% Save seasonal results to Excel
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
%% Save yearly results to Excel
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

%% Processing Completion Summary

fprintf('=== CryoSat-2 Data Processing Pipeline Completed ===\n\n');

% Display processing summary
display_processing_summary(results, config);

fprintf('\nAll output files saved in: %s\n', config.paths.output_dir);
fprintf('Processing pipeline finished.\n');

%% Cleanup
fprintf('\nCleaning up temporary variables...\n');
clear track_data enhanced_data;
fprintf('Pipeline ended.\n');

%% Helper Functions

function display_processing_summary(results, config)
% DISPLAY_PROCESSING_SUMMARY Displays processing summary

fprintf('--- Processing Summary ---\n');

% Annual Statistics
if isfield(results, 'annual') && ~isempty(results.annual)
    fprintf('Annual Statistics:\n');
    fprintf('  Processing Years: %d - %d\n', ...
        config.temporal.start_year, config.temporal.end_year);
    
    % Calculate Trend
    valid = ~isnan(results.annual.values(:,2));
    if sum(valid) > 1
        years = results.annual.values(valid, 1);
        dh = results.annual.values(valid, 2);
        p = polyfit(years, dh, 1);
        fprintf('  Overall Trend: %.3f m/yr\n', p(1));
        fprintf('  Average data points: %.0f points/year\n', mean(results.annual.values(valid, 6)));
    end
end

% Multi-year Statistics
if isfield(results, 'multiyear')
    fprintf('\nMulti-year Statistics:\n');
    fprintf('  3-Year Windows: %d (2010-2025)\n', length(results.multiyear.years));
end

% Monthly Statistics
if isfield(results, 'monthly')
    fprintf('\nMonthly Statistics:\n');
    fprintf('  Monthly data points: %d\n', size(results.monthly.values, 1));
end

% Seasonal Elevation Bands
if isfield(results, 'elevation_bands')
    fprintf('\nSeasonal Elevation Band Statistics:\n');
    fprintf('  Number of Elevation Bands: %d\n', length(results.elevation_bands.elevation_centers));
    fprintf('  Number of Time Periods: %d (3-month window)\n', size(results.elevation_bands.values, 3));
end

% Area Weighting
if isfield(results, 'area_weighted') && ~isempty(results.area_weighted)
    fprintf('\nArea Weighting:\n');
    if isfield(results.area_weighted, 'seasonal_values')
        fprintf('  Seasonal Area Weighted: %d periods x %d regions\n', ...
            size(results.area_weighted.seasonal_values, 1), ...
            size(results.area_weighted.seasonal_values, 2));
    end
    if isfield(results.area_weighted, 'yearly_values')
        fprintf('  Yearly Area Weighted: %d periods x %d regions\n', ...
            size(results.area_weighted.yearly_values, 1), ...
            size(results.area_weighted.yearly_values, 2));
    end
end

% Metadata
if isfield(results, 'metadata')
    fprintf('\nProcessing Information:\n');
    fprintf('  Processing Time: %s\n', results.metadata.processing_date);
    fprintf('  Data Points: %d\n', results.metadata.data_points);
end

fprintf('\n');

end