# Libraries
import os
from time import gmtime, strftime
import geopandas as gpd
import shapely.geometry as shpg
import xarray as xr
import numpy as np
import argparse
# import sys

# OGGM imports
from oggm import utils, workflow, tasks, graphics, global_tasks, shop
import oggm.cfg as cfg
from oggm.shop import gcm_climate, ecmwf
from oggm.core import massbalance, climate
from oggm.core.massbalance import MultipleFlowlineMassBalance
from oggm.shop import create_scieno

# Plotting
# import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
# from statsmodels.tsa.seasonal import STL
# import statsmodels.api as sm
# from statsmodels.robust.robust_linear_model import RLM
# from statsmodels.robust.norms import HuberT
import pickle
import multiprocessing
import logging

# Set plotting style
sns.set_style('ticks')
sns.set_context('notebook')

parser = argparse.ArgumentParser(description='OGGM glacier modeling for HMA region')
parser.add_argument('process_rgi_index_start', type=int, 
                    help='Start index for processing RGI IDs (1-based)')
parser.add_argument('process_rgi_index_end', type=int,
                    help='End index for processing RGI IDs (inclusive)')
parser.add_argument('use_core_num', type=int,
                    help='Number of cores to use for multiprocessing')
args = parser.parse_args()

process_rgi_index_start = args.process_rgi_index_start # strats from 1
process_rgi_index_end = args.process_rgi_index_end # 1-10, 11-20, 21-30, etc.
use_core_num = args.use_core_num

max_cores = multiprocessing.cpu_count()
if use_core_num > max_cores:
    # logger.warning(f"【{process_rgi_index_start}-{process_rgi_index_end}】Requested {use_core_num} cores exceeds available {max_cores} cores. Using {max_cores} cores instead.")
    use_core_num = max_cores
elif use_core_num < 1:
    # logger.error(f"【{process_rgi_index_start}-{process_rgi_index_end}】use_core_num must be >= 1, got {use_core_num}")
    raise ValueError(f"【{process_rgi_index_start}-{process_rgi_index_end}】use_core_num must be >= 1, got {use_core_num}")

# process_rgi_idxs, input from cmd
### Configuration for HMA region
subregion = 'HMA'
boundary_shapefile = './boundary/HMA_one.shp'
working_dir = './climate_extreme_HMA_glacier_modeling_CMIP6_ts1000'
altimetry_filepath = './Altimetry_model_input/Altimetrymb_rgiregion_result_for_OGGM0024.csv'
ref_period = '2000-01-01_2024-01-01'

era5_local_dir = {
        'tmp':'./ERA5_Land_monthly/era5_land_monthly_t2m_1950-2025_flat_HMA.nc',
        'pre':'./ERA5_Land_monthly/era5_land_monthly_prcp_1950-2025_flat_HMA.nc',
        'inv':'./ERA5_Land_monthly/era5_land_invariant_flat_HMA.nc'
}
# CMIP6 scenarios to analyze
cmip6_scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp585']
cmip6_models = ['BCC-CSM2-MR', 'CAMS-CSM1-0', 'CESM2', 'CESM2-WACCM', 
                'EC-Earth3', 'EC-Earth3-Veg', 'FGOALS-f3-L', 'GFDL-ESM4',
                'INM-CM4-8', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NorESM2-MM']
# Base URLs for CMIP6 data (update as needed)
base_url_temp = './CMIP6_monthly/GCM/{}/{}_{}_r1i1p1f1_tas.nc'
base_url_precip = './CMIP6_monthly/GCM/{}/{}_{}_r1i1p1f1_pr.nc'

extreme_period_1_start = 2005
extreme_period_1_end = 2014
extreme_period_2_start = 2015
extreme_period_2_end = 2024

rgi_regions = [13, 14, 15]
rgi_version = '62'
projection_start_year = 2025
projection_end_year = 2100
SPINUP_START_YR = 1979

### Initialize OGGM
cfg.initialize(logging_level='WARNING')

# Enable multiprocessing
cfg.PARAMS['use_multiprocessing'] = True
cfg.PARAMS['mp_processes'] = use_core_num
# Enable geometry storage for analysis
cfg.PARAMS['store_model_geometry'] = True
cfg.PARAMS['store_fl_diagnostics'] = True

# Set working directory
cfg.PATHS['working_dir'] = utils.mkdir(working_dir, reset=False) #reset to True for first use
cfg.PARAMS['cfl_min_dt'] = 30
cfg.PARAMS['continue_on_error'] = True
cfg.PATHS['rgi_dir'] = './rgi'
# cfg.PARAMS['use_intersects'] = False

# Configure our own logger
log_filename = './OGGM_ts1000-v2.log'
logger = logging.getLogger('OGGM_Custom')
logger.setLevel(logging.INFO)
for handler in logger.handlers[:]:
    logger.removeHandler(handler)

file_handler = logging.FileHandler(log_filename, mode='a')
file_handler.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
console_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.addHandler(console_handler)

logger.propagate = False # Important! If not set, the logger will propagate to the root logger

gdir_file = cfg.PATHS['working_dir'] + '/gdir_list.pkl'
if not os.path.exists(gdir_file):
    logger.error(f"【{process_rgi_index_start}-{process_rgi_index_end}】GDIR LIST file not found: {gdir_file}")
    raise FileNotFoundError(f"【{process_rgi_index_start}-{process_rgi_index_end}】GDIR LIST file not found: {gdir_file}")

with open(gdir_file, 'rb') as f:
    gdirs = pickle.load(f) # store all glacier dictionaries

warming_trend_file = cfg.PATHS['working_dir'] + '/warming_trend_0025.xlsx'
if not os.path.exists(warming_trend_file):
    logger.error(f"【{process_rgi_index_start}-{process_rgi_index_end}】Warming trend file not found: {warming_trend_file}")
    raise FileNotFoundError(f"【{process_rgi_index_start}-{process_rgi_index_end}】Warming trend file not found: {warming_trend_file}")
    
if process_rgi_index_start < 1:
    logger.error(f"【{process_rgi_index_start}-{process_rgi_index_end}】process_rgi_index_start must be >= 1, got {process_rgi_index_start}")
    raise ValueError(f"【{process_rgi_index_start}-{process_rgi_index_end}】process_rgi_index_start must be >= 1, got {process_rgi_index_start}")
if process_rgi_index_end > len(gdirs):
    logger.error(f"【{process_rgi_index_start}-{process_rgi_index_end}】process_rgi_index_end ({process_rgi_index_end}) exceeds total RGI IDs ({len(gdirs)})")
    raise ValueError(f"【{process_rgi_index_start}-{process_rgi_index_end}】process_rgi_index_end ({process_rgi_index_end}) exceeds total RGI IDs ({len(gdirs)})")
if process_rgi_index_start > process_rgi_index_end:
    logger.error(f"【{process_rgi_index_start}-{process_rgi_index_end}】process_rgi_index_start ({process_rgi_index_start}) must be <= process_rgi_index_end ({process_rgi_index_end})")
    raise ValueError(f"【{process_rgi_index_start}-{process_rgi_index_end}】process_rgi_index_start ({process_rgi_index_start}) must be <= process_rgi_index_end ({process_rgi_index_end})")


# Read HMA boundary shapefile
# basin = gpd.read_file(boundary_shapefile)
# ### Get RGI glaciers within the boundary for all specified regions
# gdf_sel = gpd.GeoDataFrame()
# for region in rgi_regions:
#     fr = utils.get_rgi_region_file(region, version=rgi_version) #local
#     gdf = gpd.read_file(fr)
#     in_bas = [basin.geometry.contains(shpg.Point(x, y))[0] for (x, y) in zip(gdf.CenLon, gdf.CenLat)]
#     gdf_region = gdf.loc[in_bas]
#     gdf_sel = gdf_sel._append(gdf_region)

# if len(gdf_sel) == 0:
#     raise ValueError(f"No glaciers found within the boundary for the specified regions")
# Initialize glacier directories from pre-processed data, withour reset
# We directly import the directories from pkl files, make sure that the rudimentary files are present in each directory dir
# gdirs = workflow.init_glacier_directories(gdf_sel) # if first use, set from_prepro_level=2 and download url

gdirs = gdirs[process_rgi_index_start-1: process_rgi_index_end] # selected the glaciers to process
logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Processing glaciers: {len(gdirs)}")

### Climate Data Processing
# Process ERA5-Land data for HMA
if not os.path.exists(os.path.join(gdirs[-1].dir, 'climate_historical.nc')):
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Processing ERA5-Land data for HMA (this may take a while)...")
    workflow.execute_entity_task(ecmwf.process_ecmwf_data, gdirs, dataset='ERA5L-LATEST',
                                download=False, local_path_dict=era5_local_dir)
    # Remove glaciers without valid climate data
    valid_gdirs = []
    for gdir in gdirs:
        fpath = gdir.dir + '/climate_historical.nc'
        try:
            tmp = xr.open_dataset(fpath)
            valid_gdirs.append(gdir)
            tmp.close()
        except:
            pass
    gdirs = valid_gdirs
# Process CMIP6 data for HMA
if not os.path.exists(os.path.join(gdirs[-1].dir, f'gcm_data_{cmip6_models[-1]}_{cmip6_scenarios[-1]}.nc')):
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Processing CMIP6 data for HMA (this may take a while)...")
    for gcm in cmip6_models:
        for ssp in cmip6_scenarios:
            # File paths for temperature and precipitation
            ft = base_url_temp.format(gcm, gcm, ssp)
            fp = base_url_precip.format(gcm, gcm, ssp)
            # Process CMIP6 data for glaciers
            workflow.execute_entity_task(
                gcm_climate.process_cmip_data, gdirs,
                filesuffix='_{}_{}'.format(gcm, ssp),
                year_range=('2000', '2019'),  # for bias correction
                fpath_temp=ft,
                fpath_precip=fp
            )
### Calibration
# Get geodetic mass balance reference data
df_ref_dmdtda = utils.get_geodetic_mb_dataframe(file_path=altimetry_filepath)

## Calibrate each glacier
logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Calibrating glaciers (this will take a while for large region)...")
for gdir in gdirs:
    df_ref_dmdtda0 = df_ref_dmdtda.loc[gdir.rgi_id] #series
    #df_ref_dmdtda0 = df_ref_dmdtda0.loc[df_ref_dmdtda0['period'] == ref_period]
    dmdtda_reference = df_ref_dmdtda0['dmdtda'] * 1000 #m/yr -> mm/yr
    dmdtda_reference_error = df_ref_dmdtda0['err_dmdtda'] * 1000
    
    # gdf_temp = gdf_sel[gdf_sel['RGIId'] == gdir.rgi_id] #gdf
    # diff = (gdf_temp['Zmax'] - gdf_temp['Zmin']).values[0]
    
    try:
        climate.mu_star_calibration_from_geodetic_mb(gdir, 
                                                    ignore_hydro_months=True,
                                                    ref_mb=dmdtda_reference,
                                                    ref_period=ref_period,
                                                    step_height_for_corr=10,
                                                    max_height_change_for_corr=5000,
                                                    min_mu_star=20,
                                                    max_mu_star=500)
    except Exception:
        pass  # Skip glaciers that fail calibration
# Remove glaciers without valid calibration
valid_gdirs = []
for gdir in gdirs:
    try:
        df = gdir.read_json('local_mustar')
        mu_star = df['mu_star_glacierwide']
        if not pd.isna(mu_star):
            valid_gdirs.append(gdir)
    except:
        pass
gdirs = valid_gdirs
logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully calibrated glaciers: {len(gdirs)}")

# Create mass balance model
# mbmod = workflow.execute_entity_task(MultipleFlowlineMassBalance, gdirs, use_inversion_flowlines=True)
# Compile mass balance data
# df = utils.compile_fixed_geometry_mass_balance(gdirs, path=False)
# Compute apparent mass balance, the MB that assumes the glacier was steady (average MB = 0)
for gdir in gdirs:
    climate.apparent_mb_from_any_mb(gdir)

# Run inversion tasks (Important!!!)
global_tasks.inversion_tasks(gdirs, glen_a=None, fs=None, filter_inversion_output=True)
# invension thickness for the initial time
# Initialize present-time glacier
workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs, filesuffix='')

gdf_sel2 = []
for gdir in gdirs:
    inversion_file = os.path.join(gdir.dir, 'inversion_output.pkl')
    if not os.path.exists(inversion_file):
        gdf_sel2.append(gdir.rgi_id)
# Remove incomplete glaciers
if len(gdf_sel2) > 0:
    gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully computed ice thickness for glaciers: {len(gdirs)}")

##### Start Dimulation
### Historical Simulation (no spinup)
# Run historical simulation
file_id = f'_hist_{process_rgi_index_start}-{process_rgi_index_end}'

if not os.path.exists(os.path.join(cfg.PATHS['working_dir'], f'run_output{file_id}.nc')):
    for gdir in gdirs:
        if os.path.exists(os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc'))
        if os.path.exists(os.path.join(gdir.dir, f'model_geometry{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'model_geometry{file_id}.nc'))
        if os.path.exists(os.path.join(gdir.dir, f'fl_diagnostics{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'fl_diagnostics{file_id}.nc'))

    workflow.execute_entity_task(tasks.run_with_hydro, gdirs,
                                run_task=tasks.run_from_climate_data,
                                climate_filename='climate_historical',
                                fixed_geometry_spinup_yr=2000,
                                ref_area_from_y0=True,
                                output_filesuffix=file_id,
                                store_monthly_hydro=False)
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    ds_hist = utils.compile_run_output(gdirs, input_filesuffix=file_id)
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully run historical simulation (no spinup) for glaciers: {len(gdirs)}")
else:
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Historical simulation (no spinup) for glaciers: {len(gdirs)} has been run before, skip it")

# Simu 1, remove warming trend for the second decade
climate_filesuffix = '_warm_rm'
file_id = f'_hist_warm_rm_{process_rgi_index_start}-{process_rgi_index_end}'

if not os.path.exists(os.path.join(cfg.PATHS['working_dir'], f'run_output{file_id}.nc')):
    for gdir in gdirs:
        if os.path.exists(os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc'))
        if os.path.exists(os.path.join(gdir.dir, f'model_geometry{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'model_geometry{file_id}.nc'))
        if os.path.exists(os.path.join(gdir.dir, f'fl_diagnostics{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'fl_diagnostics{file_id}.nc'))
        if os.path.exists(os.path.join(gdir.dir, f'climate_historical{climate_filesuffix}.nc')):
            os.remove(os.path.join(gdir.dir, f'climate_historical{climate_filesuffix}.nc'))

    workflow.execute_entity_task(
        create_scieno.remove_warming_trend, gdirs,
        ys_start=extreme_period_2_start,
        ys_end=extreme_period_2_end,
        output_filesuffix=climate_filesuffix      
    )

    workflow.execute_entity_task(
                                tasks.run_with_hydro, gdirs,
                                run_task=tasks.run_from_climate_data,
                                climate_filename='climate_historical',
                                climate_input_filesuffix=climate_filesuffix,
                                fixed_geometry_spinup_yr=2000,
                                ref_area_from_y0=True,
                                output_filesuffix=file_id,
                                store_monthly_hydro=False
                                )
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    # Compile output
    ds_hist = utils.compile_run_output(gdirs, input_filesuffix=file_id)
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully conduct hist simuation (remove warming trend for the second decade, no spinup) for glaciers: {len(gdirs)}")
else:
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Historical simulation (remove warming trend for the second decade, no spinup) for glaciers: {len(gdirs)} has been run before, skip it")

# Simu 2, add warming trend for the first decade
climate_filesuffix = '_warm_add'
file_id = f'_hist_warm_add_{process_rgi_index_start}-{process_rgi_index_end}'

if not os.path.exists(os.path.join(cfg.PATHS['working_dir'], f'run_output{file_id}.nc')):
    for gdir in gdirs:
        if os.path.exists(os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc'))
        if os.path.exists(os.path.join(gdir.dir, f'model_geometry{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'model_geometry{file_id}.nc'))
        if os.path.exists(os.path.join(gdir.dir, f'fl_diagnostics{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'fl_diagnostics{file_id}.nc'))
        if os.path.exists(os.path.join(gdir.dir, f'climate_historical{climate_filesuffix}.nc')):
            os.remove(os.path.join(gdir.dir, f'climate_historical{climate_filesuffix}.nc'))

    workflow.execute_entity_task(
        create_scieno.add_warming_trend, gdirs,
        ys_start=extreme_period_1_start,
        ys_end=extreme_period_1_end,
        output_filesuffix=climate_filesuffix      
    )

    workflow.execute_entity_task(
                                tasks.run_with_hydro, gdirs,
                                run_task=tasks.run_from_climate_data,
                                climate_filename='climate_historical',
                                climate_input_filesuffix=climate_filesuffix,
                                fixed_geometry_spinup_yr=2000,
                                ref_area_from_y0=True,
                                output_filesuffix=file_id,
                                store_monthly_hydro=False
                                )
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    # Compile output
    ds_hist = utils.compile_run_output(gdirs, input_filesuffix=file_id)
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully conduct hist simuation (add warming trend for the first decade, no spinup) for glaciers: {len(gdirs)}")
else:
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Historical simulation (add warming trend for the first decade, no spinup) for glaciers: {len(gdirs)} has been run before, skip it")

### Projection
# Normal (CMIP6)
if not os.path.exists(os.path.join(cfg.PATHS['working_dir'], f'run_output_normal_{process_rgi_index_start}-{process_rgi_index_end}.nc')):
    for gdir in gdirs:
        for ssp in cmip6_scenarios:
            for gcm in cmip6_models:
                projection_id = f'_proj_{gcm}_{ssp}_normal_{process_rgi_index_start}-{process_rgi_index_end}'
                if os.path.exists(os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc'))
                if os.path.exists(os.path.join(gdir.dir, f'model_geometry{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'model_geometry{projection_id}.nc'))
                if os.path.exists(os.path.join(gdir.dir, f'fl_diagnostics{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'fl_diagnostics{projection_id}.nc'))
    
    workflow.execute_entity_task(
                    create_scieno.repulicate_cmip6,
                    gdirs,
                    start_year=1950,
                    end_year=2100,
                    output_filesuffix='_normal',
                )

    projection_results = []
    for ssp in cmip6_scenarios:
        for gcm in cmip6_models:
            # Define file suffix for this model-scenario combination
            extreme_filesuffix = f'_{gcm}_{ssp}_normal'
            projection_id = f'_proj_{gcm}_{ssp}_normal_{process_rgi_index_start}-{process_rgi_index_end}'
                
            workflow.execute_entity_task(
                    tasks.run_with_hydro, gdirs,
                    run_task=tasks.run_from_climate_data,
                    climate_filename='gcm_data',
                    climate_input_filesuffix=extreme_filesuffix,
                    fixed_geometry_spinup_yr=2000,
                    ref_area_from_y0=True,
                    output_filesuffix=projection_id,
                    store_monthly_hydro=False
                )
            gdf_sel2 = []
            for gdir in gdirs:
                result_file = os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')
                if not os.path.exists(result_file):
                    gdf_sel2.append(gdir.rgi_id)
            if len(gdf_sel2) > 0:
                gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]

            ds_projection = utils.compile_run_output(gdirs, input_filesuffix=projection_id, path=False)
            ds_projection = ds_projection.assign_coords(GCM=gcm, SSP=ssp)
            ds_projection = ds_projection.expand_dims(['GCM', 'SSP'])
            projection_results.append(ds_projection)

    ds_all_projections = xr.combine_by_coords(projection_results, fill_value=np.nan, combine_attrs='override')
    ds_all_projections = ds_all_projections.sortby('rgi_id')
    output_path = cfg.PATHS['working_dir'] + f'/run_output_cmip6_normal_{process_rgi_index_start}-{process_rgi_index_end}.nc'
    ds_all_projections.to_netcdf(output_path)
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully conduct normal (CMIP6) projection (no spinup) for glaciers: {len(gdirs)}")
else:
    projection_id = f'_proj_{cmip6_models[-1]}_{cmip6_scenarios[-1]}_normal_{process_rgi_index_start}-{process_rgi_index_end}'
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Normal (CMIP6) projection (no spinup) for glaciers: {len(gdirs)} has been run before, skip it")

# Extreme using QDM
if not os.path.exists(os.path.join(cfg.PATHS['working_dir'], f'run_output_cmip6_repu_his_extremes_QDM_{process_rgi_index_start}-{process_rgi_index_end}.nc')):
    for gdir in gdirs:
        for ssp in cmip6_scenarios:
            for gcm in cmip6_models:
                extreme_filesuffix = f'_{gcm}_{ssp}_repu_his_extremes_QDM'
                projection_id = f'_proj_{gcm}_{ssp}_repu_his_extremes_QDM_{process_rgi_index_start}-{process_rgi_index_end}'
                if os.path.exists(os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc'))
                if os.path.exists(os.path.join(gdir.dir, f'model_geometry{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'model_geometry{projection_id}.nc'))
                if os.path.exists(os.path.join(gdir.dir, f'fl_diagnostics{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'fl_diagnostics{projection_id}.nc'))
                if os.path.exists(os.path.join(gdir.dir, f'gcm_data{extreme_filesuffix}.nc')):
                    os.remove(os.path.join(gdir.dir, f'gcm_data{extreme_filesuffix}.nc'))
    
    workflow.execute_entity_task(
                    create_scieno.simulate_future_extremes_QDM,
                    gdirs,
                    start_year=1950,
                    end_year=2100,
                    future_start_year=projection_start_year,
                    future_cooling_factor=None,
                    output_filesuffix='_repu_his_extremes_QDM',
                )

    projection_results = []
    for ssp in cmip6_scenarios:
        for gcm in cmip6_models:
            # Define file suffix for this model-scenario combination
            extreme_filesuffix = f'_{gcm}_{ssp}_repu_his_extremes_QDM'
            projection_id = f'_proj_{gcm}_{ssp}_repu_his_extremes_QDM_{process_rgi_index_start}-{process_rgi_index_end}'
                
            workflow.execute_entity_task(
                    tasks.run_with_hydro, gdirs,
                    run_task=tasks.run_from_climate_data,
                    climate_filename='gcm_data',
                    climate_input_filesuffix=extreme_filesuffix,
                    fixed_geometry_spinup_yr=2000,
                    ref_area_from_y0=True,
                    output_filesuffix=projection_id,
                    store_monthly_hydro=False
                )
            gdf_sel2 = []
            for gdir in gdirs:
                result_file = os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')
                if not os.path.exists(result_file):
                    gdf_sel2.append(gdir.rgi_id)
            if len(gdf_sel2) > 0:
                gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]

            ds_projection = utils.compile_run_output(gdirs, input_filesuffix=projection_id, path=False)
            ds_projection = ds_projection.assign_coords(GCM=gcm, SSP=ssp)
            ds_projection = ds_projection.expand_dims(['GCM', 'SSP'])
            projection_results.append(ds_projection)

    ds_all_projections = xr.combine_by_coords(projection_results, fill_value=np.nan, combine_attrs='override')
    ds_all_projections = ds_all_projections.sortby('rgi_id')
    output_path = cfg.PATHS['working_dir'] + f'/run_output_cmip6_repu_his_extremes_QDM_{process_rgi_index_start}-{process_rgi_index_end}.nc'
    ds_all_projections.to_netcdf(output_path)
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully conduct extreme projection (using QDM, no spinup) for glaciers: {len(gdirs)}")
else:
    projection_id = f'_proj_{cmip6_models[-1]}_{cmip6_scenarios[-1]}_repu_his_extremes_QDM_{process_rgi_index_start}-{process_rgi_index_end}'
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Extreme projection (using QDM, no spinup) for glaciers: {len(gdirs)} has been run before, skip it")

# Extreme using QM
if not os.path.exists(os.path.join(cfg.PATHS['working_dir'], f'run_output_cmip6_repu_his_extremes_QM_{process_rgi_index_start}-{process_rgi_index_end}.nc')):
    for gdir in gdirs:
        for ssp in cmip6_scenarios:
            for gcm in cmip6_models:
                extreme_filesuffix = f'_{gcm}_{ssp}_repu_his_extremes_QM'
                projection_id = f'_proj_{gcm}_{ssp}_repu_his_extremes_QM_{process_rgi_index_start}-{process_rgi_index_end}'
                if os.path.exists(os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc'))
                if os.path.exists(os.path.join(gdir.dir, f'model_geometry{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'model_geometry{projection_id}.nc'))
                if os.path.exists(os.path.join(gdir.dir, f'fl_diagnostics{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'fl_diagnostics{projection_id}.nc'))
                if os.path.exists(os.path.join(gdir.dir, f'gcm_data{extreme_filesuffix}.nc')):
                    os.remove(os.path.join(gdir.dir, f'gcm_data{extreme_filesuffix}.nc'))
    
    workflow.execute_entity_task(
                    create_scieno.simulate_future_extremes_Detrend_QM,
                    gdirs,
                    start_year=1950,
                    end_year=2100,
                    future_start_year=projection_start_year,
                    future_cooling_factor=None,
                    output_filesuffix='_repu_his_extremes_QM',
                )

    projection_results = []
    for ssp in cmip6_scenarios:
        for gcm in cmip6_models:
            # Define file suffix for this model-scenario combination
            extreme_filesuffix = f'_{gcm}_{ssp}_repu_his_extremes_QM'
            projection_id = f'_proj_{gcm}_{ssp}_repu_his_extremes_QM_{process_rgi_index_start}-{process_rgi_index_end}'
                
            workflow.execute_entity_task(
                    tasks.run_with_hydro, gdirs,
                    run_task=tasks.run_from_climate_data,
                    climate_filename='gcm_data',
                    climate_input_filesuffix=extreme_filesuffix,
                    fixed_geometry_spinup_yr=2000,
                    ref_area_from_y0=True,
                    output_filesuffix=projection_id,
                    store_monthly_hydro=False
                )
            gdf_sel2 = []
            for gdir in gdirs:
                result_file = os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')
                if not os.path.exists(result_file):
                    gdf_sel2.append(gdir.rgi_id)
            if len(gdf_sel2) > 0:
                gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]

            ds_projection = utils.compile_run_output(gdirs, input_filesuffix=projection_id, path=False)
            ds_projection = ds_projection.assign_coords(GCM=gcm, SSP=ssp)
            ds_projection = ds_projection.expand_dims(['GCM', 'SSP'])
            projection_results.append(ds_projection)

    ds_all_projections = xr.combine_by_coords(projection_results, fill_value=np.nan, combine_attrs='override')
    ds_all_projections = ds_all_projections.sortby('rgi_id')
    output_path = cfg.PATHS['working_dir'] + f'/run_output_cmip6_repu_his_extremes_QM_{process_rgi_index_start}-{process_rgi_index_end}.nc'
    ds_all_projections.to_netcdf(output_path)
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully conduct extreme projection (using QM, no spinup) for glaciers: {len(gdirs)}")
else:
    projection_id = f'_proj_{cmip6_models[-1]}_{cmip6_scenarios[-1]}_repu_his_extremes_QM_{process_rgi_index_start}-{process_rgi_index_end}'
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Extreme projection (using QM, no spinup) for glaciers: {len(gdirs)} has been run before, skip it")

### Spin up
## Historical Simulation (spinup)
# Dynamic spining up
spinup_filesuffix = f'_spinup'
for gdir in gdirs:
    if os.path.exists(os.path.join(gdir.dir, f'model_diagnostics{spinup_filesuffix}.nc')):
        os.remove(os.path.join(gdir.dir, f'model_diagnostics{spinup_filesuffix}.nc'))
    if os.path.exists(os.path.join(gdir.dir, f'model_geometry{spinup_filesuffix}.nc')):
        os.remove(os.path.join(gdir.dir, f'model_geometry{spinup_filesuffix}.nc'))
    if os.path.exists(os.path.join(gdir.dir, f'fl_diagnostics{spinup_filesuffix}.nc')):
        os.remove(os.path.join(gdir.dir, f'fl_diagnostics{spinup_filesuffix}.nc'))

workflow.execute_entity_task(
    tasks.run_dynamic_spinup, gdirs,
    climate_input_filesuffix=None, 
    spinup_start_yr=SPINUP_START_YR,         
    minimise_for='area',              
    precision_percent=3,
    precision_absolute=3,       
    add_fixed_geometry_spinup=True,              
    output_filesuffix=spinup_filesuffix          
)

path = cfg.PATHS['working_dir'] + '/per_glacier'
gdf_sel2 = []
for gdir in gdirs:
    result_file = os.path.join(gdir.dir, 'model_diagnostics'+spinup_filesuffix+'.nc')
    if not os.path.exists(result_file):
        gdf_sel2.append(gdir.rgi_id)
if len(gdf_sel2) > 0:
    gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully spin up for glaciers: {len(gdirs)}")

# Normal
file_id = f'_hist_spinup_{process_rgi_index_start}-{process_rgi_index_end}'
if not os.path.exists(os.path.join(cfg.PATHS['working_dir'], f'run_output{file_id}.nc')):
    for gdir in gdirs:
        if os.path.exists(os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc'))
        if os.path.exists(os.path.join(gdir.dir, f'model_geometry{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'model_geometry{file_id}.nc'))
        if os.path.exists(os.path.join(gdir.dir, f'fl_diagnostics{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'fl_diagnostics{file_id}.nc'))
    
    workflow.execute_entity_task(tasks.run_with_hydro, gdirs,
                                run_task=tasks.run_from_climate_data,
                                climate_filename='climate_historical',
                                init_model_filesuffix=spinup_filesuffix,
                                init_model_yr=1999,
                                ref_area_from_y0=True,
                                output_filesuffix=file_id,
                                store_monthly_hydro=False)
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    ds_hist = utils.compile_run_output(gdirs, input_filesuffix=file_id)
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully run historical simulation (spinup) for glaciers: {len(gdirs)}")
else:
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Historical simulation (spinup) for glaciers: {len(gdirs)} has been run before, skip it")

# Simu 1, remove warming trend for the second decade
file_id = f'_hist_warm_rm_spinup_{process_rgi_index_start}-{process_rgi_index_end}'
climate_filesuffix = '_warm_rm'
if not os.path.exists(os.path.join(cfg.PATHS['working_dir'], f'run_output{file_id}.nc')):
    for gdir in gdirs:
        if os.path.exists(os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc'))
        if os.path.exists(os.path.join(gdir.dir, f'model_geometry{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'model_geometry{file_id}.nc'))
        if os.path.exists(os.path.join(gdir.dir, f'fl_diagnostics{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'fl_diagnostics{file_id}.nc'))
    
    workflow.execute_entity_task(
                                tasks.run_with_hydro, gdirs,
                                run_task=tasks.run_from_climate_data,
                                climate_filename='climate_historical',
                                climate_input_filesuffix=climate_filesuffix,
                                init_model_filesuffix=spinup_filesuffix,
                                init_model_yr=1999,
                                ref_area_from_y0=True,
                                output_filesuffix=file_id,
                                store_monthly_hydro=False
                                )
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    # Compile output
    ds_hist = utils.compile_run_output(gdirs, input_filesuffix=file_id)
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully conduct hist simuation (remove warming trend for the second decade, spinup) for glaciers: {len(gdirs)}")
else:
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Historical simulation (remove warming trend for the second decade, spinup) for glaciers: {len(gdirs)} has been run before, skip it")

# Simu 2, add warming trend for the first decade
climate_filesuffix = '_warm_add'
file_id = f'_hist_warm_add_spinup_{process_rgi_index_start}-{process_rgi_index_end}'
if not os.path.exists(os.path.join(cfg.PATHS['working_dir'], f'run_output{file_id}.nc')):
    for gdir in gdirs:
        if os.path.exists(os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc'))
        if os.path.exists(os.path.join(gdir.dir, f'model_geometry{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'model_geometry{file_id}.nc'))
        if os.path.exists(os.path.join(gdir.dir, f'fl_diagnostics{file_id}.nc')):
            os.remove(os.path.join(gdir.dir, f'fl_diagnostics{file_id}.nc'))

    workflow.execute_entity_task(
                                tasks.run_with_hydro, gdirs,
                                run_task=tasks.run_from_climate_data,
                                climate_filename='climate_historical',
                                climate_input_filesuffix=climate_filesuffix,
                                init_model_filesuffix=spinup_filesuffix,
                                init_model_yr=1999,
                                ref_area_from_y0=True,
                                output_filesuffix=file_id,
                                store_monthly_hydro=False
                                )
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    # Compile output
    ds_hist = utils.compile_run_output(gdirs, input_filesuffix=file_id)
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully conduct hist simuation (add warming trend for the first decade, spinup) for glaciers: {len(gdirs)}")
else:
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{file_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Historical simulation (add warming trend for the first decade, spinup) for glaciers: {len(gdirs)} has been run before, skip it")

### Projection
# Normal (CMIP6)
if not os.path.exists(os.path.join(cfg.PATHS['working_dir'], f'run_output_cmip6_normal_spinup_{process_rgi_index_start}-{process_rgi_index_end}.nc')):
    for gdir in gdirs:
        for ssp in cmip6_scenarios:
            for gcm in cmip6_models:
                projection_id = f'_proj_{gcm}_{ssp}_normal_spinup_{process_rgi_index_start}-{process_rgi_index_end}'
                if os.path.exists(os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc'))
                if os.path.exists(os.path.join(gdir.dir, f'model_geometry{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'model_geometry{projection_id}.nc'))
                if os.path.exists(os.path.join(gdir.dir, f'fl_diagnostics{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'fl_diagnostics{projection_id}.nc'))
    
    projection_results = []
    for ssp in cmip6_scenarios:
        for gcm in cmip6_models:
            # Define file suffix for this model-scenario combination
            extreme_filesuffix = f'_{gcm}_{ssp}_normal'
            projection_id = f'_proj_{gcm}_{ssp}_normal_spinup_{process_rgi_index_start}-{process_rgi_index_end}'
                
            workflow.execute_entity_task(
                    tasks.run_with_hydro, gdirs,
                    run_task=tasks.run_from_climate_data,
                    climate_filename='gcm_data',
                    climate_input_filesuffix=extreme_filesuffix,
                    init_model_filesuffix=spinup_filesuffix,
                    init_model_yr=1999,
                    ref_area_from_y0=True,
                    output_filesuffix=projection_id,
                    store_monthly_hydro=False
                )
            gdf_sel2 = []
            for gdir in gdirs:
                result_file = os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')
                if not os.path.exists(result_file):
                    gdf_sel2.append(gdir.rgi_id)
            if len(gdf_sel2) > 0:
                gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]

            ds_projection = utils.compile_run_output(gdirs, input_filesuffix=projection_id, path=False)
            ds_projection = ds_projection.assign_coords(GCM=gcm, SSP=ssp)
            ds_projection = ds_projection.expand_dims(['GCM', 'SSP'])
            projection_results.append(ds_projection)

    ds_all_projections = xr.combine_by_coords(projection_results, fill_value=np.nan, combine_attrs='override')
    ds_all_projections = ds_all_projections.sortby('rgi_id')
    output_path = cfg.PATHS['working_dir'] + f'/run_output_cmip6_normal_spinup_{process_rgi_index_start}-{process_rgi_index_end}.nc'
    ds_all_projections.to_netcdf(output_path)
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully conduct normal (CMIP6) projection (spinup) for glaciers: {len(gdirs)}")
else:
    projection_id = f'_proj_{cmip6_models[-1]}_{cmip6_scenarios[-1]}_normal_spinup_{process_rgi_index_start}-{process_rgi_index_end}'
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Normal (CMIP6) projection (spinup) for glaciers: {len(gdirs)} has been run before, skip it")

# Extreme using QDM
if not os.path.exists(os.path.join(cfg.PATHS['working_dir'], f'run_output_cmip6_repu_his_extremes_QDM_spinup_{process_rgi_index_start}-{process_rgi_index_end}.nc')):
    for gdir in gdirs:
        for ssp in cmip6_scenarios:
            for gcm in cmip6_models:
                projection_id = f'_proj_{gcm}_{ssp}_repu_his_extremes_QDM_spinup_{process_rgi_index_start}-{process_rgi_index_end}'
                if os.path.exists(os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc'))
                if os.path.exists(os.path.join(gdir.dir, f'model_geometry{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'model_geometry{projection_id}.nc'))
                if os.path.exists(os.path.join(gdir.dir, f'fl_diagnostics{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'fl_diagnostics{projection_id}.nc'))
    
    projection_results = []
    for ssp in cmip6_scenarios:
        for gcm in cmip6_models:
            # Define file suffix for this model-scenario combination
            extreme_filesuffix = f'_{gcm}_{ssp}_repu_his_extremes_QDM'
            projection_id = f'_proj_{gcm}_{ssp}_repu_his_extremes_QDM_spinup_{process_rgi_index_start}-{process_rgi_index_end}'
                
            workflow.execute_entity_task(
                    tasks.run_with_hydro, gdirs,
                    run_task=tasks.run_from_climate_data,
                    climate_filename='gcm_data',
                    climate_input_filesuffix=extreme_filesuffix,
                    init_model_filesuffix=spinup_filesuffix,
                    init_model_yr=1999,
                    ref_area_from_y0=True,
                    output_filesuffix=projection_id,
                    store_monthly_hydro=False
                )
            gdf_sel2 = []
            for gdir in gdirs:
                result_file = os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')
                if not os.path.exists(result_file):
                    gdf_sel2.append(gdir.rgi_id)
            if len(gdf_sel2) > 0:
                gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]

            ds_projection = utils.compile_run_output(gdirs, input_filesuffix=projection_id, path=False)
            ds_projection = ds_projection.assign_coords(GCM=gcm, SSP=ssp)
            ds_projection = ds_projection.expand_dims(['GCM', 'SSP'])
            projection_results.append(ds_projection)

    ds_all_projections = xr.combine_by_coords(projection_results, fill_value=np.nan, combine_attrs='override')
    ds_all_projections = ds_all_projections.sortby('rgi_id')
    output_path = cfg.PATHS['working_dir'] + f'/run_output_cmip6_repu_his_extremes_QDM_spinup_{process_rgi_index_start}-{process_rgi_index_end}.nc'
    ds_all_projections.to_netcdf(output_path)
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully conduct extreme projection (using QDM, spinup) for glaciers: {len(gdirs)}")
else:
    projection_id = f'_proj_{cmip6_models[-1]}_{cmip6_scenarios[-1]}_repu_his_extremes_QDM_spinup_{process_rgi_index_start}-{process_rgi_index_end}'
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Extreme projection (using QDM, spinup) for glaciers: {len(gdirs)} has been run before, skip it")

# Extreme using QM
if not os.path.exists(os.path.join(cfg.PATHS['working_dir'], f'run_output_cmip6_repu_his_extremes_QM_spinup_{process_rgi_index_start}-{process_rgi_index_end}.nc')):
    for gdir in gdirs:
        for ssp in cmip6_scenarios:
            for gcm in cmip6_models:
                projection_id = f'_proj_{gcm}_{ssp}_repu_his_extremes_QM_spinup_{process_rgi_index_start}-{process_rgi_index_end}'
                if os.path.exists(os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc'))
                if os.path.exists(os.path.join(gdir.dir, f'model_geometry{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'model_geometry{projection_id}.nc'))
                if os.path.exists(os.path.join(gdir.dir, f'fl_diagnostics{projection_id}.nc')):
                    os.remove(os.path.join(gdir.dir, f'fl_diagnostics{projection_id}.nc'))
    
    projection_results = []
    for ssp in cmip6_scenarios:
        for gcm in cmip6_models:
            # Define file suffix for this model-scenario combination
            extreme_filesuffix = f'_{gcm}_{ssp}_repu_his_extremes_QM'
            projection_id = f'_proj_{gcm}_{ssp}_repu_his_extremes_QM_spinup_{process_rgi_index_start}-{process_rgi_index_end}'
                
            workflow.execute_entity_task(
                    tasks.run_with_hydro, gdirs,
                    run_task=tasks.run_from_climate_data,
                    climate_filename='gcm_data',
                    climate_input_filesuffix=extreme_filesuffix,
                    init_model_filesuffix=spinup_filesuffix,
                    init_model_yr=1999,
                    ref_area_from_y0=True,
                    output_filesuffix=projection_id,
                    store_monthly_hydro=False
                )
            gdf_sel2 = []
            for gdir in gdirs:
                result_file = os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')
                if not os.path.exists(result_file):
                    gdf_sel2.append(gdir.rgi_id)
            if len(gdf_sel2) > 0:
                gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]

            ds_projection = utils.compile_run_output(gdirs, input_filesuffix=projection_id, path=False)
            ds_projection = ds_projection.assign_coords(GCM=gcm, SSP=ssp)
            ds_projection = ds_projection.expand_dims(['GCM', 'SSP'])
            projection_results.append(ds_projection)

    ds_all_projections = xr.combine_by_coords(projection_results, fill_value=np.nan, combine_attrs='override')
    ds_all_projections = ds_all_projections.sortby('rgi_id')
    output_path = cfg.PATHS['working_dir'] + f'/run_output_cmip6_repu_his_extremes_QM_spinup_{process_rgi_index_start}-{process_rgi_index_end}.nc'
    ds_all_projections.to_netcdf(output_path)
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Successfully conduct extreme projection (using QM, spinup) for glaciers: {len(gdirs)}")
else:
    projection_id = f'_proj_{cmip6_models[-1]}_{cmip6_scenarios[-1]}_repu_his_extremes_QM_spinup_{process_rgi_index_start}-{process_rgi_index_end}'
    gdf_sel2 = []
    for gdir in gdirs:
        result_file = os.path.join(gdir.dir, f'model_diagnostics{projection_id}.nc')
        if not os.path.exists(result_file):
            gdf_sel2.append(gdir.rgi_id)
    if len(gdf_sel2) > 0:
        gdirs[:] = [gdir for gdir in gdirs if gdir.rgi_id not in gdf_sel2]
    logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】Extreme projection (using QM, spinup) for glaciers: {len(gdirs)} has been run before, skip it")


# Export summary statistics
logger.info(f"【{process_rgi_index_start}-{process_rgi_index_end}】HMA glacier modeling complete!")