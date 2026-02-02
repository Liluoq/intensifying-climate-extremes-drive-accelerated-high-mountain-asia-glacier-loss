"""
Example script demonstrating how to read and analyze the zarr output.
"""

import zarr
import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import logging

# Simple logging setup for standalone script
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def read_zarr_basic(zarr_path: str):
    """
    Read zarr store using basic zarr interface.
    
    Parameters:
    -----------
    zarr_path : str
        Path to zarr store
    """
    print("=" * 80)
    print("READING ZARR STORE - BASIC INTERFACE")
    print("=" * 80)
    
    # Open zarr store
    root = zarr.open(zarr_path, mode='r')
    
    print("\nZarr structure:")
    print(root.tree())
    
    print("\nMetadata:")
    for key, value in root.attrs.items():
        print(f"  {key}: {value}")
    
    # Access arrays
    data = root['data']
    glacier_ids = root['glacier_id'][:]
    variables = root['variable'][:]
    times = root['time'][:]
    
    print(f"\nData shape: {data.shape}")
    print(f"  Glaciers: {len(glacier_ids)}")
    print(f"  Variables: {len(variables)}")
    print(f"  Timesteps: {len(times)}")
    
    print(f"\nVariable names: {list(variables)}")
    
    # Example: Access data for first glacier, all variables, first timestep
    first_glacier_data = data[0, :, 0]
    print(f"\nFirst glacier, first timestep:")
    for var_idx, var_name in enumerate(variables):
        print(f"  {var_name}: {first_glacier_data[var_idx]:.4f}")
    
    return root


def read_zarr_xarray(zarr_path: str):
    """
    Read zarr store using xarray for more convenient data access.
    
    Parameters:
    -----------
    zarr_path : str
        Path to zarr store
    """
    print("\n" + "=" * 80)
    print("READING ZARR STORE - XARRAY INTERFACE")
    print("=" * 80)
    
    # Open zarr store
    root = zarr.open(zarr_path, mode='r')
    
    # Create xarray Dataset
    data = root['data'][:]
    glacier_ids = root['glacier_id'][:]
    variables = root['variable'][:]
    times = pd.to_datetime(root['time'][:])
    
    # Build xarray Dataset
    ds = xr.Dataset(
        data_vars={
            'statistics': (['glacier_id', 'variable', 'time'], data)
        },
        coords={
            'glacier_id': glacier_ids,
            'variable': variables,
            'time': times
        },
        attrs=dict(root.attrs)
    )
    
    print("\nDataset structure:")
    print(ds)
    
    # Example queries
    print("\n" + "-" * 80)
    print("Example 1: Select specific variable (temperature_2m)")
    print("-" * 80)
    temp_data = ds['statistics'].sel(variable='temperature_2m')
    print(f"Shape: {temp_data.shape}")
    print(f"Mean temperature across all glaciers and times: {temp_data.mean().values:.2f}")
    
    print("\n" + "-" * 80)
    print("Example 2: Select specific glacier")
    print("-" * 80)
    first_glacier_id = glacier_ids[0]
    glacier_data = ds['statistics'].sel(glacier_id=first_glacier_id)
    print(f"Glacier ID: {first_glacier_id}")
    print(f"Shape: {glacier_data.shape}")
    print(f"Variables: {list(glacier_data.coords['variable'].values)}")
    
    print("\n" + "-" * 80)
    print("Example 3: Time series for one glacier and one variable")
    print("-" * 80)
    timeseries = ds['statistics'].sel(
        glacier_id=first_glacier_id,
        variable='temperature_2m'
    )
    print(f"Time series length: {len(timeseries)}")
    print(f"First 5 values:")
    for t, val in zip(times[:5], timeseries.values[:5]):
        print(f"  {t.strftime('%Y-%m-%d')}: {val:.4f}")
    
    print("\n" + "-" * 80)
    print("Example 4: Compute statistics across glaciers")
    print("-" * 80)
    mean_by_var = ds['statistics'].mean(dim=['glacier_id', 'time'])
    print("Mean values by variable (across all glaciers and times):")
    for var, val in zip(variables, mean_by_var.values):
        print(f"  {var}: {val:.4f}")
    
    print("\n" + "-" * 80)
    print("Example 5: Export subset to pandas DataFrame")
    print("-" * 80)
    # Convert a slice to DataFrame
    df = ds['statistics'].isel(time=slice(0, 5)).to_dataframe().reset_index()
    print(f"DataFrame shape: {df.shape}")
    print(f"\nFirst few rows:")
    print(df.head(10))
    
    return ds


def analyze_zarr(zarr_path: str):
    """
    Perform basic analysis on the zarr data.
    
    Parameters:
    -----------
    zarr_path : str
        Path to zarr store
    """
    print("\n" + "=" * 80)
    print("DATA QUALITY ANALYSIS")
    print("=" * 80)
    
    root = zarr.open(zarr_path, mode='r')
    data = root['data'][:]
    
    # Check for missing values
    total_values = data.size
    missing_values = np.sum(np.isnan(data))
    valid_values = total_values - missing_values
    
    print(f"\nData completeness:")
    print(f"  Total values: {total_values:,}")
    print(f"  Valid values: {valid_values:,} ({100 * valid_values / total_values:.2f}%)")
    print(f"  Missing values: {missing_values:,} ({100 * missing_values / total_values:.2f}%)")
    
    # Statistics by variable
    variables = root['variable'][:]
    print(f"\nStatistics by variable (excluding NaN):")
    print(f"{'Variable':<30} {'Mean':>12} {'Std':>12} {'Min':>12} {'Max':>12}")
    print("-" * 80)
    
    for var_idx, var_name in enumerate(variables):
        var_data = data[:, var_idx, :]
        var_data_valid = var_data[~np.isnan(var_data)]
        
        if len(var_data_valid) > 0:
            mean_val = np.mean(var_data_valid)
            std_val = np.std(var_data_valid)
            min_val = np.min(var_data_valid)
            max_val = np.max(var_data_valid)
            
            print(f"{var_name:<30} {mean_val:>12.4f} {std_val:>12.4f} {min_val:>12.4f} {max_val:>12.4f}")
    
    # Check data integrity
    print(f"\nData integrity checks:")
    print(f"  Contains inf: {np.any(np.isinf(data))}")
    print(f"  All finite (excluding NaN): {np.all(np.isfinite(data[~np.isnan(data)]))}")
    
    return data


def main():
    """Main function to demonstrate zarr reading."""
    # Update this path to match your config
    zarr_path = "/WORK/Data/ncc_glacier/era5_land/extracted/glacier_stats.zarr"
    
    # Check if zarr store exists
    if not Path(zarr_path).exists():
        print(f"Error: Zarr store not found at {zarr_path}")
        print("Please run extract_main.py first to generate the data.")
        return
    
    # Demonstrate different reading methods
    root = read_zarr_basic(zarr_path)
    ds = read_zarr_xarray(zarr_path)
    data = analyze_zarr(zarr_path)
    
    print("\n" + "=" * 80)
    print("READING COMPLETE")
    print("=" * 80)
    print("\nTo use this data in your analysis:")
    print("  1. Use zarr.open() for low-level access")
    print("  2. Use xarray for high-level, labeled data operations")
    print("  3. Convert to pandas DataFrame for tabular analysis")
    print("  4. Use .sel() method for intuitive data selection")


if __name__ == "__main__":
    main()

