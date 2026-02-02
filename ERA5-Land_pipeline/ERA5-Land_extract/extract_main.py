import extract_util
import geopandas as gpd
import rasterio as rio
import zarr
import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
import json
import logging
from typing import List, Dict, Tuple
import os

# Get logger for this module (will be configured by setup_logging)
logger = logging.getLogger(__name__)


def load_config(config_path: str = "config.json") -> dict:
    """Load configuration from JSON file."""
    config_file = Path(__file__).parent.parent / config_path
    with open(config_file, 'r') as f:
        config = json.load(f)
    logger.info(f"Loaded configuration from {config_file}")
    return config


def scan_raster_files(raster_dir: str, file_pattern: str) -> List[Tuple[str, datetime]]:
    """
    Scan directory for raster files and extract dates from filenames.
    
    Parameters:
    -----------
    raster_dir : str
        Directory containing raster files
    file_pattern : str
        Glob pattern for raster files (e.g., 'ERA5-Land_*.tif')
        
    Returns:
    --------
    List[Tuple[str, datetime]]
        List of (file_path, date) tuples sorted by date
    """
    logger.info(f"Scanning for raster files in {raster_dir} with pattern {file_pattern}")
    
    raster_dir_path = Path(raster_dir)
    raster_files = []
    
    for file_path in raster_dir_path.glob(file_pattern):
        # Extract date from filename (e.g., ERA5-Land_1990-01-01.tif -> 1990-01-01)
        filename = file_path.stem
        try:
            date_str = filename.split('_')[1]  # Assumes format: ERA5-Land_YYYY-MM-DD
            date = datetime.strptime(date_str, '%Y-%m-%d')
            raster_files.append((str(file_path), date))
        except Exception as e:
            logger.error(f"Error parsing date from filename {filename}: {e}")
            raise e
    
    # Sort by date
    raster_files.sort(key=lambda x: x[1])
    
    logger.info(f"Found {len(raster_files)} raster files")
    if raster_files:
        logger.info(f"Date range: {raster_files[0][1]} to {raster_files[-1][1]}")
    
    return raster_files


def initialize_zarr_store(
    zarr_path: str,
    glacier_ids: List[str],
    variable_names: List[str],
    timestamps: List[datetime],
    chunk_size_time: int = 100,
    chunk_size_glacier: int = 1,
) -> zarr.Group:
    """
    Initialize zarr store with pre-allocated arrays and coordinates.
    
    Parameters:
    -----------
    zarr_path : str
        Path to zarr store
    glacier_ids : List[str]
        List of glacier IDs
    variable_names : List[str]
        List of variable names
    timestamps : List[datetime]
        List of timestamps
    chunk_size_time : int
        Chunk size along time dimension
        
    Returns:
    --------
    zarr.Group
        Initialized zarr group
    """
    logger.info(f"Initializing zarr store at {zarr_path}")
    
    # Create zarr store (overwrite if exists)
    store = zarr.DirectoryStore(zarr_path)
    root = zarr.group(store, overwrite=True)
    
    n_glaciers = len(glacier_ids)
    n_variables = len(variable_names)
    n_times = len(timestamps)
    
    logger.info(f"Dimensions: glaciers={n_glaciers}, variables={n_variables}, time={n_times}")
    
    # Define chunks: all glaciers and variables, but chunked in time
    chunks = (chunk_size_glacier, n_variables, chunk_size_time)
    
    # Create main data array
    data_array = root.create_dataset(
        'data',
        shape=(n_glaciers, n_variables, n_times),
        chunks=chunks,
        dtype='float32',
        fill_value=np.nan,
        compressor=zarr.Blosc(cname='zstd', clevel=3, shuffle=2)
    )
    
    logger.info(f"Created data array with shape {data_array.shape} and chunks {data_array.chunks}")
    
    # Create coordinate arrays
    max_glacier_id_length = max(len(glacier_id) for glacier_id in glacier_ids)
    root.create_dataset('glacier_id', data=np.array(glacier_ids, dtype=f'U{max_glacier_id_length}'), overwrite=True)
    max_variable_name_length = max(len(variable_name) for variable_name in variable_names)
    root.create_dataset('variable', data=np.array(variable_names, dtype=f'U{max_variable_name_length}'), overwrite=True)
    
    # Convert timestamps to numpy datetime64
    timestamps_np = np.array([pd.Timestamp(ts) for ts in timestamps], dtype='datetime64[ns]')
    root.create_dataset('time', data=timestamps_np, overwrite=True)
    
    # Store metadata
    root.attrs['dimensions'] = ['glacier_id', 'variable', 'time']
    root.attrs['created'] = datetime.now().isoformat()
    root.attrs['description'] = 'Extracted ERA5-Land statistics for HMA glaciers'
    
    logger.info("Zarr store initialization complete")
    
    return root


def run_parallel_extraction(config_path: str = "config.json", test_mode: bool = False):
    """
    Main pipeline for parallel extraction and zarr storage.
    
    Parameters:
    -----------
    config_path : str
        Path to configuration file
    test_mode : bool
        If True, only process first 10 files for testing
    """
    # Load configuration
    config = load_config(config_path)
    
    # Set up logging (to both console and file) with multiprocessing support
    log_file = config.get('log_file')
    if log_file:
        # Make log file path relative to the config file location
        log_file_path = Path(__file__).parent.parent / log_file
    else:
        log_file_path = None
    
    # Setup returns log_queue for worker processes
    log_queue = extract_util.setup_logging(log_file=str(log_file_path) if log_file_path else None)
    
    logger.info("Logging initialized")
    if log_file_path:
        logger.info(f"Log file: {log_file_path}")
    
    # Extract configuration parameters
    raster_dir = config['input_raster_dir']
    glacier_shapefile = config['glacier_shapefile']
    zarr_path = config['output_zarr_path']
    variable_names = config['variable_names']
    stats_op = config['stats_operation']
    batch_size = config['batch_size']
    max_workers = config['max_workers']
    file_pattern = config['raster_file_pattern']
    glacier_id_col = config['glacier_id_column']
    chunk_size_glacier = config['chunk_size_glacier']
    
    logger.info("="*80)
    logger.info("STARTING PARALLEL ERA5-LAND EXTRACTION PIPELINE")
    logger.info("="*80)
    
    # PHASE 1: INITIALIZATION
    logger.info("\n[PHASE 1] Initialization")
    logger.info("-"*80)
    
    # Scan for raster files
    raster_files = scan_raster_files(raster_dir, file_pattern)
    
    if not raster_files:
        logger.error("No raster files found! Exiting.")
        return
    
    # Test mode: limit to first 10 files
    if test_mode:
        logger.warning("TEST MODE: Processing only first 10 files")
        raster_files = raster_files[:10]
    
    file_paths = [f[0] for f in raster_files]
    timestamps = [f[1] for f in raster_files]
    
    # Load glacier shapefile
    logger.info(f"Loading glacier shapefile: {glacier_shapefile}")
    glaciers = gpd.read_file(glacier_shapefile)
    
    # Ensure CRS is set
    if glaciers.crs is None:
        glaciers = glaciers.set_crs("EPSG:4326")
        logger.info("Set CRS to EPSG:4326")
    
    glacier_ids = glaciers[glacier_id_col].tolist()
    n_glaciers = len(glacier_ids)
    logger.info(f"Loaded {n_glaciers} glaciers")
    
    # Validate number of variables with first raster
    logger.info("Validating raster structure...")
    with rio.open(file_paths[0]) as src:
        n_bands = src.count
    
    if n_bands != len(variable_names):
        logger.error(f"Mismatch: Raster has {n_bands} bands but config specifies {len(variable_names)} variables")
        logger.error(f"Please update variable_names in config.json")
        return
    
    logger.info(f"Validated: {n_bands} bands match {len(variable_names)} variable names")
    
    # Initialize zarr store
    zarr_root = initialize_zarr_store(
        zarr_path=zarr_path,
        glacier_ids=glacier_ids,
        variable_names=variable_names,
        timestamps=timestamps,
        chunk_size_time=batch_size,
        chunk_size_glacier=chunk_size_glacier
    )
    
    # PHASE 2: BATCH PROCESSING
    logger.info("\n[PHASE 2] Batch Processing")
    logger.info("-"*80)
    
    n_files = len(file_paths)
    n_batches = (n_files + batch_size - 1) // batch_size
    
    logger.info(f"Processing {n_files} files in {n_batches} batches (batch_size={batch_size})")
    
    for batch_idx in range(n_batches):
        batch_start = batch_idx * batch_size
        batch_end = min(batch_start + batch_size, n_files)
        batch_files = file_paths[batch_start:batch_end]
        batch_timestamps = timestamps[batch_start:batch_end]
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Batch {batch_idx+1}/{n_batches} - Processing files {batch_start} to {batch_end-1}")
        logger.info(f"{'='*60}")
        
        # Extract batch in parallel (pass log_queue for worker logging)
        batch_results = extract_util.extract_batch_parallel(
            raster_paths=batch_files,
            vectors=glaciers,
            operation=[stats_op],
            extract_strategy='feature-sequential',
            include_cols=[glacier_id_col],
            max_workers=max_workers,
            log_queue=log_queue
        )
        
        if not batch_results:
            logger.warning(f"Batch {batch_idx+1} produced no results!")
            continue
        
        # Sort results by filename to match timestamp order
        results_dict = {path: df for path, df in batch_results}
        sorted_results = [
            results_dict.get(path) 
            for path in batch_files 
            if path in results_dict
        ]
        assert len(sorted_results) == len(batch_files), f"Number of results does not match number of files, got {len(sorted_results)} results for {len(batch_files)} files, check the extraction results"
        
        # Convert results to array and write to zarr
        logger.info(f"Converting {len(sorted_results)} results to array format...")
        
        for i, (file_path, result_df) in enumerate(zip(batch_files, sorted_results)):
            if result_df is None:
                logger.warning(f"Skipping {file_path} - no data")
                continue
            
            time_idx = batch_start + i
            
            # Ensure glacier order matches zarr glacier_id coordinate
            result_df = result_df.set_index(glacier_id_col).reindex(glacier_ids)
            
            # Extract mean columns for each band
            stat_cols = [col for col in result_df.columns if col not in [glacier_id_col]]
            
            if len(stat_cols) != len(variable_names):
                logger.error(f"Mismatch in {file_path}: {len(stat_cols)} stat columns vs {len(variable_names)} variables")
                continue
            
            # Convert to array: shape (n_glaciers, n_variables)
            data_slice = result_df[stat_cols].values
            
            # Write to zarr: [:, :, time_idx]
            zarr_root['data'][:, :, time_idx] = data_slice
            
        logger.info(f"Batch {batch_idx+1} written to zarr (time indices {batch_start} to {batch_end-1})")
        logger.info(f"Progress: {batch_end}/{n_files} files ({100*batch_end/n_files:.1f}%)")
    
    # PHASE 3: COMPLETION
    logger.info("\n[PHASE 3] Completion")
    logger.info("-"*80)
    
    logger.info("Verifying zarr store...")
    data_array = zarr_root['data']
    logger.info(f"Final data shape: {data_array.shape}")
    logger.info(f"Data type: {data_array.dtype}")
    logger.info(f"Fill value: {data_array.fill_value}")
    logger.info(f"Compression: {data_array.compressor}")
    
    # Calculate statistics
    sample_data = data_array[:, :, :min(10, data_array.shape[2])]
    n_valid = np.sum(~np.isnan(sample_data))
    n_total = sample_data.size
    logger.info(f"Sample validation (first 10 timesteps): {n_valid}/{n_total} valid values ({100*n_valid/n_total:.1f}%)")
    
    logger.info("\n" + "="*80)
    logger.info("EXTRACTION PIPELINE COMPLETE!")
    logger.info("="*80)
    logger.info(f"Output saved to: {zarr_path}")
    logger.info(f"Access with: zarr.open('{zarr_path}', mode='r')")
    logger.info(f"Or with xarray: xr.open_zarr('{zarr_path}')")


def test_extract():
    """Legacy test function - kept for compatibility."""
    rgi_hma_path = "/WORK/Data/ncc_glacier/glacier_boundary/rgi_hma/rgi_hma.shp"
    rgi_hma = gpd.read_file(rgi_hma_path).set_crs("EPSG:4326")
    
    test_era5_land_path = "/WORK/Data/ncc_glacier/era5_land/raw/ERA5-Land_1990-01-01.tif"
    
    test_result = extract_util.extract_single(
        raster_path=test_era5_land_path,
        vectors=rgi_hma,
        extract_strategy='feature-sequential',
        progress=True,
        include_cols=['RGIId'],
    )
    
    print(f"Test extraction complete. Result shape: {test_result.shape}")
    return test_result


if __name__ == "__main__":
    # Run in test mode by default (processes only 10 files)
    # Set test_mode=False to process all files
    run_parallel_extraction(config_path="config.json", test_mode=False)