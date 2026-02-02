"""
Aggregate daily ERA5-Land glacier statistics to monthly means.

This script reads the zarr output from the extraction pipeline, computes
calendar month means, and saves results to an Excel file with one sheet
per variable.
"""

import zarr
import numpy as np
import pandas as pd
from pathlib import Path
import json
import logging
from datetime import datetime
from typing import Tuple, Dict, List
import openpyxl
from openpyxl.utils.dataframe import dataframe_to_rows

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)
logger.addHandler(logging.FileHandler('aggregate_to_monthly.log'))
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)
logger.info("Starting monthly aggregation...")
logger.info("="*80)


def load_config(config_path: str = "config.json") -> dict:
    """Load configuration from JSON file."""
    config_file = Path(__file__).parent.parent / config_path
    with open(config_file, 'r') as f:
        config = json.load(f)
    logger.info(f"Loaded configuration from {config_file}")
    return config


def load_zarr_data(zarr_path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, pd.DatetimeIndex]:
    """
    Load data from zarr store.
    
    Parameters:
    -----------
    zarr_path : str
        Path to zarr store
        
    Returns:
    --------
    tuple
        (data, glacier_ids, variable_names, timestamps)
    """
    logger.info(f"Loading zarr data from {zarr_path}")
    
    # Open zarr store
    root = zarr.open(zarr_path, mode='r')
    
    # Load data and coordinates
    data = root['data'][:]
    glacier_ids = root['glacier_id'][:]
    variable_names = root['variable'][:]
    timestamps = pd.to_datetime(root['time'][:])
    
    logger.info(f"Loaded data shape: {data.shape}")
    logger.info(f"  Glaciers: {len(glacier_ids)}")
    logger.info(f"  Variables: {len(variable_names)}")
    logger.info(f"  Time steps: {len(timestamps)}")
    logger.info(f"  Date range: {timestamps[0]} to {timestamps[-1]}")
    
    return data, glacier_ids, variable_names, timestamps


def aggregate_to_monthly(
    data: np.ndarray,
    timestamps: pd.DatetimeIndex,
    log_file: str = None
) -> Tuple[np.ndarray, pd.PeriodIndex, List[Tuple]]:
    """
    Aggregate daily data to monthly means.
    
    Parameters:
    -----------
    data : np.ndarray
        Daily data with shape (n_glaciers, n_variables, n_times)
    timestamps : pd.DatetimeIndex
        Daily timestamps
    log_file : str, optional
        Path to log file for NaN warnings
        
    Returns:
    --------
    tuple
        (monthly_data, month_periods, nan_warnings)
        - monthly_data: shape (n_glaciers, n_variables, n_months)
        - month_periods: PeriodIndex with month labels
        - nan_warnings: list of (glacier_idx, var_idx, month) tuples with all NaN
    """
    logger.info("Aggregating to monthly means...")
    
    n_glaciers, n_variables, n_times = data.shape
    
    # Create month periods from timestamps
    month_periods = timestamps.to_period('M').unique().sort_values()
    n_months = len(month_periods)
    
    logger.info(f"Number of unique months: {n_months}")
    logger.info(f"Month range: {month_periods[0]} to {month_periods[-1]}")
    
    # Initialize output array
    monthly_data = np.full((n_glaciers, n_variables, n_months), np.nan, dtype=np.float32)
    
    # Track NaN warnings
    nan_warnings = []
    
    # Convert timestamps to month periods for grouping
    time_months = timestamps.to_period('M')
    
    # Aggregate for each month
    for month_idx, month in enumerate(month_periods):
        # Find indices for this month
        month_mask = time_months == month
        month_data = data[:, :, month_mask]
        
        # Compute mean for each glacier-variable pair
        for glacier_idx in range(n_glaciers):
            for var_idx in range(n_variables):
                values = month_data[glacier_idx, var_idx, :]
                
                # Check if all values are NaN
                if np.all(np.isnan(values)):
                    nan_warnings.append((glacier_idx, var_idx, month))
                    # Keep as NaN in output
                else:
                    # Compute mean, ignoring NaN
                    monthly_data[glacier_idx, var_idx, month_idx] = np.nanmean(values)
        
        if (month_idx + 1) % 12 == 0:
            logger.info(f"  Processed {month_idx + 1}/{n_months} months ({100*(month_idx+1)/n_months:.1f}%)")
    
    logger.info(f"Aggregation complete")
    logger.info(f"Found {len(nan_warnings)} glacier-variable-month combinations with all NaN values")
    
    # Log NaN warnings to file
    if log_file and nan_warnings:
        log_nan_warnings(nan_warnings, log_file)
    
    return monthly_data, month_periods, nan_warnings


def log_nan_warnings(
    nan_warnings: List[Tuple],
    log_file: str
):
    """
    Log NaN warnings to file.
    
    Parameters:
    -----------
    nan_warnings : list
        List of (glacier_idx, var_idx, month) tuples
    log_file : str
        Path to log file
    """
    logger.info(f"Writing NaN warnings to {log_file}")
    
    with open(log_file, 'w') as f:
        f.write("NaN Warnings - Months with All NaN Values\n")
        f.write("=" * 80 + "\n")
        f.write(f"Total warnings: {len(nan_warnings)}\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("Format: glacier_index, variable_index, month\n\n")
        
        for glacier_idx, var_idx, month in nan_warnings:
            f.write(f"Glacier {glacier_idx}, Variable {var_idx}, Month {month}\n")
        
        f.write("\n" + "=" * 80 + "\n")
        f.write(f"End of warnings. Total: {len(nan_warnings)}\n")


def write_to_excel(
    monthly_data: np.ndarray,
    glacier_ids: np.ndarray,
    variable_names: np.ndarray,
    month_periods: pd.PeriodIndex,
    output_path: str
):
    """
    Write monthly data to Excel file with one sheet per variable.
    
    Parameters:
    -----------
    monthly_data : np.ndarray
        Monthly data with shape (n_glaciers, n_variables, n_months)
    glacier_ids : np.ndarray
        Glacier identifiers
    variable_names : np.ndarray
        Variable names
    month_periods : pd.PeriodIndex
        Month labels
    output_path : str
        Path to output Excel file
    """
    logger.info(f"Writing data to Excel file: {output_path}")
    
    # Convert month periods to strings (YYYY-MM format)
    month_labels = [str(m) for m in month_periods]
    
    # Create Excel writer
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        
        # Create one sheet per variable
        for var_idx, var_name in enumerate(variable_names):
            logger.info(f"  Writing sheet: {var_name}")
            
            # Extract data for this variable: (n_glaciers, n_months)
            var_data = monthly_data[:, var_idx, :]
            
            # Create DataFrame
            df = pd.DataFrame(
                var_data,
                index=glacier_ids,
                columns=month_labels
            )
            
            # Set index name
            df.index.name = 'glacier_id'
            
            # Write to Excel
            # Clean variable name for sheet name (Excel has limits)
            sheet_name = var_name[:31]  # Excel sheet name max length
            df.to_excel(writer, sheet_name=sheet_name)
    
    logger.info(f"Excel file written successfully")


def main(config_path: str = "config.json"):
    """
    Main function to orchestrate monthly aggregation.
    
    Parameters:
    -----------
    config_path : str
        Path to configuration file
    """
    logger.info("="*80)
    logger.info("STARTING MONTHLY AGGREGATION")
    logger.info("="*80)
    
    # Load configuration
    config = load_config(config_path)
    
    zarr_path = config['output_zarr_path']
    monthly_output = config.get('monthly_output_xlsx', 
                                 '/WORK/Data/ncc_glacier/era5_land/extracted/glacier_stats_monthly.xlsx')
    nan_log_file = config.get('monthly_aggregation_log',
                              '/WORK/Data/ncc_glacier/era5_land/extracted/monthly_nan_warnings.log')
    variables_to_aggregate = config.get('variables_to_aggregate', None)
    
    logger.info(f"Input zarr: {zarr_path}")
    logger.info(f"Output Excel: {monthly_output}")
    logger.info(f"NaN log: {nan_log_file}")
    
    # Create output directory if needed
    output_dir = Path(monthly_output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    data, glacier_ids, variable_names, timestamps = load_zarr_data(zarr_path)
    
    # Filter variables if specified
    if variables_to_aggregate is not None:
        logger.info(f"Filtering variables to aggregate: {variables_to_aggregate}")
        
        # Find indices of selected variables
        var_indices = []
        selected_var_names = []
        for var_name in variables_to_aggregate:
            if var_name in variable_names:
                idx = list(variable_names).index(var_name)
                var_indices.append(idx)
                selected_var_names.append(var_name)
            else:
                logger.warning(f"Variable '{var_name}' not found in zarr data, skipping")
        
        if not var_indices:
            logger.error("No valid variables selected! Using all variables.")
        else:
            # Filter data to selected variables
            data = data[:, var_indices, :]
            variable_names = np.array(selected_var_names)
            logger.info(f"Selected {len(var_indices)} variables: {list(variable_names)}")
    else:
        logger.info("No variable filter specified, aggregating all variables")
    
    # Aggregate to monthly
    monthly_data, month_periods, nan_warnings = aggregate_to_monthly(
        data, timestamps, nan_log_file
    )
    
    # delete the zarr in RAM
    del data
    
    # Write to Excel
    write_to_excel(
        monthly_data,
        glacier_ids,
        variable_names,
        month_periods,
        monthly_output
    )
    
    # Summary
    logger.info("\n" + "="*80)
    logger.info("AGGREGATION COMPLETE")
    logger.info("="*80)
    logger.info(f"Input: {len(timestamps)} daily timesteps")
    logger.info(f"Output: {len(month_periods)} monthly timesteps")
    logger.info(f"Glaciers: {len(glacier_ids)}")
    logger.info(f"Variables: {len(variable_names)}")
    logger.info(f"NaN warnings: {len(nan_warnings)}")
    logger.info(f"\nOutput saved to: {monthly_output}")
    if nan_warnings:
        logger.info(f"NaN warnings logged to: {nan_log_file}")


if __name__ == "__main__":
    main()

