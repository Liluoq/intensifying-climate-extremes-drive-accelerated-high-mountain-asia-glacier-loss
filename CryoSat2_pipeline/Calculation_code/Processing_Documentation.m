--- START OF FILE Processing_Documentation.m ---

%% CryoSat-2 Data Processing Detailed Documentation
%
% Reference: revised_ICESat2_code/Processing_Documentation.m
% This script provides technical details for CryoSat-2 data processing.

%% 1. Core Technical Difference: 11x11 Neighborhood Average

fprintf('=== Key Technical Features of CryoSat-2 ===\n\n');

% Neighborhood Size Comparison
neighbor_comparison = [
    "Data Source", "Footprint Size", "Neighborhood Size", "Pixel Count", "Code";
    "ICESat-2", "~17m", "3x3", "9", "for j=-1:1; for k=-1:1";
    "CryoSat-2", "~300m", "11x11", "121", "for j=-5:5; for k=-5:5"
];

disp('Neighborhood Size Comparison:');
disp(neighbor_comparison);

fprintf('\nCalculation Process:\n');
fprintf('CryoSat-2 Footprint: ~300m\n');
fprintf('NASADEM Resolution: 30m\n');
fprintf('Required Pixels: 300m/30m ≈ 10\n');
fprintf('Neighborhood Size: 11x11 (to ensure coverage)\n');
fprintf('Total Pixels: 11x11 = 121\n\n');

%% 2. Code Implementation Example

fprintf('=== 11x11 Neighborhood Averaging Code Implementation ===\n\n');

% Simulation Example
fprintf('Example Code:\n');
fprintf('-----------------------------------------------------\n');
fprintf('neighbor_size = 11;\n');
fprintf('half_size = 5;\n');
fprintf('\n');
fprintf('for i = 1:size(data, 1)\n');
fprintf('    [row, col] = latlon2pix(ref_srtm, data(i,2), data(i,3));\n');
fprintf('    row = ceil(row);\n');
fprintf('    col = ceil(col);\n');
fprintf('    \n');
fprintf('    ele = 0.0;\n');
fprintf('    for j = -5:5  %% 11 pixels\n');
fprintf('        for k = -5:5  %% 11 pixels\n');
fprintf('            ele = ele + double(srtm(row+j, col+k));\n');
fprintf('        end\n');
fprintf('    end\n');
fprintf('    \n');
fprintf('    data(i, 23) = ele / 121;  %% Average\n');
fprintf('end\n');
fprintf('-----------------------------------------------------\n\n');

%% 3. Data Format Description

fprintf('=== Data Format Description ===\n\n');

% Phase 1 Output Format
stage1_format = [
    "Col No", "Field", "Description";
    "1", "Year", "Extracted from filename";
    "2", "Latitude", "WGS84";
    "3", "Longitude", "WGS84";
    "4", "Obs Elevation", "CryoSat-2 altimetry value";
    "5", "Time", "Days relative to 2000-01-01";
    "6", "Geoid", "Geoid Height";
    "7", "Surface Elev", "elevation - geoid";
    "8", "HMA22 Region", "1-22";
    "9", "GTN Region", "RGI Region";
    "10", "0.5° Grid", "Grid ID";
    "11", "1° Grid", "Grid ID";
    "12", "TP Flag", "Tibetan Plateau";
    "13", "Glacier Flag", "0=Glacier, 1=Non-glacier"
];

disp('Phase 1 Output Format (HMA_CryoSat2_glacier.txt):');
disp(stage1_format);

% Phase 2 Output Format
stage2_format = [
    "Col No", "Field", "Description";
    "23", "NASADEM Elev", "11x11 Neighborhood Average (Key!)";
    "24", "Elev Change", "Surface Elev - NASADEM(11x11)"
];

disp('Phase 2 New Columns (HMA_CryoSat2_update.txt):');
disp(stage2_format);

%% 4. Quality Control Strategy

fprintf('\n=== Quality Control Strategy ===\n\n');

quality_control_steps = [
    "Step", "Filter Condition", "Threshold";
    "1. Elev Range", "Obs Elevation", "<8900m";
    "2. Boundary Filter", "HMA Boundary", "inpolygon";
    "3. Glacier Mask", "RGI Mask", "glacier_flag=0";
    "4. Elev Change", "abs(dh)", "<400m";
    "5. Percentile", "dh Range", "5%-85%";
    "6. Median Filter", "abs(dh-median)", "<75m"
];

disp('Quality Control Steps:');
disp(quality_control_steps);

%% 5. Time Processing

fprintf('\n=== Time Processing Description ===\n\n');

fprintf('Time Reference Base: 2000-01-01\n');
fprintf('Data Period: 2010-2021\n');
fprintf('Time Unit: Days\n\n');

fprintf('Time Conversion Formula:\n');
fprintf('  NetCDF Timestamp (seconds) -> Days:\n');
fprintf('  days = floor(time_20_ku / 86400)\n\n');
fprintf('  Days -> Date:\n');
fprintf('  date = datetime(2000,1,1) + days(days_value)\n\n');

fprintf('Annual Time Windows:\n');
for year = 2010:2:2020
    year_start = datenum(year, 1, 1) - datenum(2000, 1, 1);
    year_end = datenum(year+1, 1, 1) - datenum(2000, 1, 1);
    fprintf('  Year %d: Days %d - %d\n', year, year_start, year_end);
end

%% 6. Statistical Methods

fprintf('\n=== Robust Statistical Methods ===\n\n');

fprintf('Median Filter Method:\n');
fprintf('  med = median(dh);\n');
fprintf('  valid = abs(dh - med) < 75;\n');
fprintf('  mean_robust = mean(dh(valid));\n\n');

fprintf('Percentile Method:\n');
fprintf('  lower = prctile(dh, 5);\n');
fprintf('  upper = prctile(dh, 85);\n');
fprintf('  valid = (dh >= lower) & (dh <= upper);\n');
fprintf('  Retain approx 88%% of data\n\n');

%% 7. Complete Comparison with ICESat-2

fprintf('=== Complete Comparison: CryoSat-2 vs ICESat-2 ===\n\n');

full_comparison = [
    "Feature", "CryoSat-2", "ICESat-2";
    "Altimetry Principle", "Radar Altimeter", "Laser Altimeter";
    "Footprint Size", "~300m", "~17m";
    "NASADEM Neighborhood", "11x11 (121)", "3x3 (9)";
    "Data Period", "2010-2021", "2018-2024+";
    "Time Reference", "2000-01-01", "2018-01-01";
    "Repeat Cycle", "369 days", "91 days";
    "Penetration Correction", "Not required", "Required";
    "Data Density", "Lower", "Higher";
    "Quality Control", "Relatively loose", "Relatively strict"
];

disp(full_comparison);

%% 8. Processing Flowchart

fprintf('\n=== Data Processing Workflow ===\n\n');

fprintf('L2I/TEMPO NetCDF Files\n');
fprintf('        |\n');
fprintf('[extract_cryosat2_tracks] - Extract RGI Glacier Regions\n');
fprintf('        |\n');
fprintf('HMA_CryoSat2_glacier.txt\n');
fprintf('        |\n');
fprintf('[enhance_cryosat2_data] - 11x11 Neighborhood NASADEM\n');
fprintf('        |\n');
fprintf('HMA_CryoSat2_update.txt\n');
fprintf('        |\n');
fprintf('[compute_cryosat2_statistics]\n');
fprintf('        |\n');
fprintf('Analysis Results (MAT)\n');
fprintf('        |\n');
fprintf('[generate_cryosat2_plots]\n');
fprintf('        |\n');
fprintf('Charts (PNG/FIG)\n\n');

%% 9. Frequently Asked Questions

fprintf('=== Frequently Asked Questions ===\n\n');

fprintf('Q1: Why does CryoSat-2 use an 11x11 neighborhood?\n');
fprintf('A1: Because the CryoSat-2 footprint is ~300m, much larger than ICESat-2 (~17m),\n');
fprintf('    requiring a larger neighborhood to match its spatial resolution.\n\n');

fprintf('Q2: Can I use a 3x3 neighborhood for CryoSat-2?\n');
fprintf('A2: No! This will cause spatial mismatch and unreliable results.\n\n');

fprintf('Q3: What if NetCDF reading fails?\n');
fprintf('A3: Check if variable names are correct (L2I uses _poca_20_ku suffix).\n\n');

fprintf('Q4: How to verify if neighborhood averaging is correct?\n');
fprintf('A4: Check that NASADEM values in column 23 are smooth,\n');
fprintf('    and close to single-point NASADEM values but smoother.\n\n');

%% 10. Completion Note

fprintf('=== Documentation Explanation Completed ===\n');
fprintf('Please run the main script: CryoSat2_Main_Processing_Pipeline.m\n');
fprintf('Remember key point: 11x11 neighborhood, not 3x3!\n');