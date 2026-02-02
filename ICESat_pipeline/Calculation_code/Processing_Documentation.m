%% ICESat Data Processing Workflow Documentation
%
% This document describes the complete ICESat data processing workflow in detail.
% Reference: revised_CryoSat2/Processing_Documentation.m
%
% Creation date: 2025-10-18
% Version: v1.0

%% ========================================================================
%  I. Data Overview
%  ========================================================================

% ICESat (Ice, Cloud, and land Elevation Satellite)
% - Launch date: January 2003
% - End of mission: October 2009
% - Laser altimetry system: GLAS (Geoscience Laser Altimeter System)
% - Data product: GLA14 (L2 Global Land Surface Altimetry Data)
% - Footprint size: ~70 m
% - Along-track spacing: ~170 m
% - Repeat cycle: 91 days (theoretical)

%% ========================================================================
%  II. Data Characteristics
%  ========================================================================

% 1. Temporal characteristics
%    - Discontinuous observations: 19 campaigns
%    - Each campaign lasts about 1–4 months
%    - Several months gap between campaigns
%
% 2. Spatial characteristics
%    - Small footprint (~70 m) → high spatial resolution
%    - Large track spacing → limited spatial coverage
%    - Suitable for monitoring glacier elevation changes
%
% 3. Quality control
%    - Cloud flag: must be 0 (cloud-free)
%    - Saturation: <= 255
%    - Slope: <= 40 degrees
%    - Curvature: abs < 4

%% ========================================================================
%  III. Key Differences Compared to CryoSat-2
%  ========================================================================

% +------------------+------------------+------------------+
% | Item             | ICESat           | CryoSat-2        |
% +------------------+------------------+------------------+
% | Time span        | 2003-2009        | 2010-2025        |
% | Footprint size   | ~70m             | ~300m            |
% | NASADEM window   | 3×3 (9 pixels)   | 11×11 (121 pixels)|
% | Observation mode | 19 campaigns     | Continuous        |
% | Data format      | HDF5             | NetCDF            |
% | Lon/Lat order    | lat(1), lon(2)   | lon(1), lat(2)    |
% +------------------+------------------+------------------+

%% ========================================================================
%  IV. Processing Workflow
%  ========================================================================

% Stage 1: Data extraction (optional)
% ┌─────────────────────┐
% │ GLA14 HDF files     │
% └──────────┬──────────┘
%            │ extract_icesat_tracks.m
%            ↓
% ┌─────────────────────┐
% │ Raw TXT (~20 columns)   │ HMA_ICESat_rgiregion.txt
% └──────────┬──────────┘
%
% Stage 2: Data enhancement
%            │ enhance_icesat_data.m
%            ├─→ Add NASADEM (3×3 neighborhood)
%            ├─→ Add grid and PDD
%            ├─→ Add region identifiers
%            ├─→ Compute elevation change
%            ↓
% ┌─────────────────────┐
% │ Enhanced TXT (28 columns)     │ HMA_ICESAT_rgiregion_update.txt
% └──────────┬──────────┘
%
% Stage 3: Statistical analysis
%            │ compute_icesat_statistics.m
%            ├─→ Annual statistics (2003-2009)
%            ├─→ 3-year moving windows
%            ├─→ Campaign statistics (19 windows)
%            ↓
% ┌─────────────────────┐
% │ Statistical results (MAT+Excel)│
% └─────────────────────┘

%% ========================================================================
%  V. Core Technique: 3×3 Neighborhood Mean
%  ========================================================================

% Why does ICESat use a 3×3 neighborhood?
%
% 1. Footprint matching
%    - ICESat footprint diameter ~70 m
%    - NASADEM resolution 30 m
%    - 3×3 neighborhood = 90 m × 90 m ≈ footprint area
%
% 2. Spatial representativeness
%    - Terrain variability within the footprint
%    - Spatial integration effect of the laser return
%
% 3. Comparison: CryoSat-2 uses an 11×11 neighborhood
%    - CryoSat-2 footprint ~300 m
%    - Requires a larger neighborhood for matching

% 实现代码：
neighbor_size = 3;
half_size = 1;

ele = 0.0;
for j = -1:1
    for k = -1:1
        ele = ele + single(srtm(row+j, col+k));
    end
end
nasadem_avg = ele / 9;

%% ========================================================================
%  VI. Quality-Control Strategy
%  ========================================================================

% Two-step filtering strategy:

% Step 1: Overall quality control
% - Remove saturation_index > 255
% - Remove abs(dh) > 200 m

% Step 2: Conditional filters in statistics
% Basic conditions:
%   - Within the time window
%   - Slope <= 40 degrees
%   - abs(curvature) < 4
%   - cloud_flag == 0

% Two-step filtering:
%   Step 1: abs(dh) < 100 m  (initial filter)
%   Step 2: abs(dh - median) < 50 m  (median-relative filter)

% Example code:
% row1 = find(base_conditions & abs(dh) < 100);
% med = median(dh(row1));
% row2 = find(base_conditions & abs(dh - med) < 50);
% final_mean = mean(dh(row2));

%% ========================================================================
%  VII. Definition of Time Windows
%  ========================================================================

% 1. Annual statistics
%    Time window: from year-01-01 to year+1-10-01
%    Example: 2003 data = 2003-01-01 to 2004-10-01
%
% 2. 3-year moving window
%    Time window: from year-01-01 to year+3-01-01
%    Example: 2003-2005 = 2003-01-01 to 2006-01-01
%
% 3. Campaign statistics
%    19 predefined observation windows
%    Example: Campaign 1 = 2003-02-20 to 2003-03-29

campaigns = [
    2003.191781, datenum(2003,02,20)-datenum(2000,1,1), datenum(2003,03,29)-datenum(2000,1,1)+1;
    2003.810959, datenum(2003,09,25)-datenum(2000,1,1), datenum(2003,11,19)-datenum(2000,1,1)+1;
    % ... total 19 windows
];

%% ========================================================================
%  VIII. Regional Classification System
%  ========================================================================

% Regional classifications used in ICESat processing:
%
% 1. HMA22 - 22 mountain subregions
%    - Based on mountain ranges
%    - Used for comparing elevation changes across regions
%
% 2. GTN/RGI - 15 glacier regions
%    - Based on the RGI glacier inventory
%    - Used for glacier elevation-change statistics
%
% 3. HMA4 - 4 large regions
%    - Macro-scale regional division
%    - Used for large-scale analysis
%
% 4. HMA6 - 6 seasonal-pattern regions
%    - Based on climate and seasonal characteristics
%    - Used for seasonal analysis
%
% 5. TP - Tibetan Plateau flag
%    - 0/1 indicator
%    - Used to compare inside vs outside the Plateau

%% ========================================================================
%  IX. Output Products
%  ========================================================================

% 1. Enhanced data file
%    HMA_ICESAT_rgiregion_update.txt
%    - 28 columns
%    - Includes NASADEM, region identifiers, elevation change, etc.
%
% 2. MAT-format results
%    ICESat_analysis_results.mat
%    - results.annual: annual statistics
%    - results.multiyear: 3-year window statistics
%    - results.campaigns: campaign statistics
%    - results.metadata: metadata
%
% 3. Excel-format results
%    - ICESat_yearly_results.xlsx
%      * Sheet 1: year_value (elevation-change values)
%      * Sheet 2: year_count (number of data points)
%      * Sheet 3: year_uncert (uncertainty)
%    
%    - ICESat_3year_window_results.xlsx
%      * 3-year window statistics
%    
%    - ICESat_campaign_results.xlsx
%      * Statistics for the 19 campaigns

%% ========================================================================
%  X. Usage Examples
%  ========================================================================

% Basic usage:
cd('revised_ICESat');
run ICESat_Main_Processing_Pipeline.m

% Custom configuration:
config = load_icesat_config();
config.quality.max_slope = 35;  % Modify slope threshold
config.temporal.start_year = 2004;  % Modify start year

% Run data enhancement only:
config = load_icesat_config();
input_file = 'G:\HMA\ICESat_Result\HMA_ICESat_rgiregion.txt';
enhanced_data = enhance_icesat_data(input_file, config);

% Run statistical analysis only:
config = load_icesat_config();
enhanced_file = 'G:\HMA\ICESat_Result\HMA_ICESAT_rgiregion_update.txt';
results = compute_icesat_statistics(enhanced_file, config);

%% ========================================================================
%  XI. Frequently Asked Questions (FAQ)
%  ========================================================================

% Q1: Why do ICESat and CryoSat-2 use different longitude/latitude column orders?
% A1: This is due to historical reasons: the original extraction scripts defined the columns differently. The refactored code preserves the original order.
%     ICESat: lat (column 1), lon (column 2)
%     CryoSat-2: lon (column 1), lat (column 2)

% Q2: How should missing raw data be handled?
% A2: If GLA14 raw data are not available:
%     1. Use an already extracted TXT file
%     2. Skip the extract_icesat_tracks step
%     3. Start directly from enhance_icesat_data

% Q3: Why use column 22 instead of column 25 for statistics?
% A3: Column 22 is the initial elevation change (without PDD), which more directly reflects surface change.
%     Column 25 includes the PDD correction and is mainly used for penetration-depth analysis.

% Q4: Can the campaign time windows be modified?
% A4: Yes, they can be changed in config.temporal.campaigns.
%     However, it is recommended to use the predefined 19 windows, which correspond to the standard ICESat observation periods.

%% ========================================================================
%  XII. References
%  ========================================================================

% 1. ICESat data products
%    Zwally, H. J., et al. (2014). GLAS/ICESat L2 Global Land Surface
%    Altimetry Data (GLA14), Version 34. NASA NSIDC DAAC.
%
% 2. Refactored reference code
%    - revised_CryoSat2/: refactored CryoSat-2 code
%    - HMA_ICESat_COMPUTE_20230603.m: original statistics code
%    - HMA_footprint_track_20230213.m: data enhancement reference
%
% 3. DEM data
%    NASA JPL (2020). NASADEM Merged DEM Global 1 arc second V001.
%
% 4. Glacier inventory
%    RGI Consortium (2017). Randolph Glacier Inventory - A Dataset
%    of Global Glacier Outlines, Version 6.

%% ========================================================================
%  XIII. Version History
%  ========================================================================

% v1.0 (2025-10-18)
% - Initial version
% - Core processing workflow completed
% - Refactored structure based on CryoSat-2
% - Implemented 3×3 NASADEM neighborhood processing
% - Implemented annual / 3-year window / campaign statistics

%% ========================================================================
%  XIV. Planned Extensions
%  ========================================================================

% 1. Elevation-band analysis
%    - 50 m or 100 m elevation bands
%    - Elevation change in different bands
%    - Seasonal elevation-band analysis
%
% 2. Spatial grid statistics
%    - 0.5° or 1° grids
%    - Gridded elevation-change trends
%
% 3. Visualization capabilities
%    - Spatial distribution maps
%    - Time-series plots
%    - Regional comparison plots
%
% 4. Comparison with other datasets
%    - Linkage with CryoSat-2 data
%    - DEM differencing comparisons
%    - Validation with in-situ measurements

%% ========================================================================
%  XV. Acknowledgements
%  ========================================================================

% Thanks to the original code developers
% Thanks to the refactored CryoSat-2 code for providing a reference framework
% Thanks to NSIDC for providing the ICESat data products

fprintf('ICESat data processing workflow documentation loaded.\n');
fprintf('For detailed usage instructions please refer to README.md and Data_Column_Description.md\n');

