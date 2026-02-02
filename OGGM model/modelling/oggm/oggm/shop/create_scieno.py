from math import floor, ceil
import xarray as xr
import pandas as pd
import numpy as np
# import copy
from oggm.utils import entity_task
from oggm.workflow import execute_entity_task
import logging
import oggm.cfg as cfg
import os
from netCDF4 import Dataset

log = logging.getLogger(__name__)

@entity_task(log)
def sensitivity_hist_modelling(gdir, 
                                year=2023,
                                var='temp',
                                months=np.arange(1, 13),
                                change=-0.5,
                                output_filesuffix=''):
    """
    Modify climate variable (temperature or precipitation) for specified year and months.
    
    Parameters
    ----------
    gdir : GlacierDirectory
        Glacier directory object
    year : int
        Year to modify
    var : str
        Variable to modify ('temp' or 'prcp')
    months : array
        Months to modify (1-12)
    change : float
        Amount to add to the variable
    output_filesuffix : str
        Suffix for output filename
    
    Returns
    -------
    bool
        True if successful
    """
    if os.path.exists(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix)):
        os.remove(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))

    climate_file = gdir.get_filepath('climate_historical')
    tmp = xr.open_dataset(climate_file)

    for month in months:
        tmp[var].loc[dict(time=slice(f'{year}-{month:02d}-01', f'{year}-{month:02d}-28'))] = \
                    tmp[var].sel(time=slice(f'{year}-{month:02d}-01', f'{year}-{month:02d}-28')) + change
    
    tmp.to_netcdf(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))
    tmp.close()
    return True


def _get_warming_trends():
    """Get warming trends, loading from cache if needed."""

    cache_file = os.path.join(cfg.PATHS['working_dir'], 'warming_trend_0025.xlsx')
    
    if os.path.exists(cache_file):
        warming_trends = pd.read_excel(cache_file, index_col=0)
        return warming_trends
    else:
        raise FileNotFoundError(f"Warming trends file not found: {cache_file}")

def _ERA5_hist_scores(sheet_name):
    """Get ERA5 historical scores, loading from cache if needed."""

    cache_file = os.path.join(cfg.PATHS['working_dir'], 'glacier_extreme_scores_ERA5.xlsx')
    
    if os.path.exists(cache_file):
        scores = pd.read_excel(cache_file, index_col=0, header=0, sheet_name=sheet_name) # The first row is the header, the first column is the index
        return scores
    else:
        raise FileNotFoundError(f"Score file not found: {cache_file}")

@entity_task(log)
def first_dec_hist_modelling(gdir, 
                            first_dec_start=2005,
                            first_dec_end=2014,
                            second_dec_start=2015,
                            second_dec_end=2024,
                            output_filesuffix=''):
    """
    Create extreme scenario: replace second decade with first decade data, then apply warming trend.
    
    Parameters
    ----------
    gdir : GlacierDirectory
        Glacier directory object
    first_dec_start : int
        First decade start year
    first_dec_end : int
        First decade end year
    second_dec_start : int
        Second decade start year
    second_dec_end : int
        Second decade end year
    output_filesuffix : str
        Suffix for output filename
    
    Returns
    -------
    bool
        True if successful
    """
    warming_trends = _get_warming_trends()

    if os.path.exists(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix)):
        os.remove(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))
                            
    climate_file = gdir.get_filepath('climate_historical')
    tmp = xr.open_dataset(climate_file)
    trend = warming_trends.loc[gdir.rgi_id] # from 1 to 12

    tmp['prcp'].loc[dict(time=slice(f'{second_dec_start}-01', f'{second_dec_end}-12'))] = \
                tmp['prcp'].sel(time=slice(f'{first_dec_start}-01', f'{first_dec_end}-12')).values
    tmp['temp'].loc[dict(time=slice(f'{second_dec_start}-01', f'{second_dec_end}-12'))] = \
                tmp['temp'].sel(time=slice(f'{first_dec_start}-01', f'{first_dec_end}-12')).values
    # consider warming trend for the second ten years
    for year in range(second_dec_start, second_dec_end+1):
        years_offset = year - second_dec_start 
        
        for month in range(1, 13):
            month_data = tmp['temp'].sel(time=slice(f'{year}-{month:02d}-01', f'{year}-{month:02d}-28')).values
            warming_per_month = trend[month] * years_offset
            tmp['temp'].loc[dict(time=slice(f'{year}-{month:02d}-01', f'{year}-{month:02d}-28'))] = \
                month_data + warming_per_month

    tmp.to_netcdf(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))
    tmp.close()
    return True

@entity_task(log)
def second_dec_hist_modelling(gdir, 
                            first_dec_start=2005,
                            first_dec_end=2014,
                            second_dec_start=2015,
                            second_dec_end=2024,
                            output_filesuffix=''):
    """
    Create extreme scenario: replace first decade with second decade data, then remove warming trend.
    
    Parameters
    ----------
    gdir : GlacierDirectory
        Glacier directory object
    first_dec_start : int
        First decade start year
    first_dec_end : int
        First decade end year
    second_dec_start : int
        Second decade start year
    second_dec_end : int
        Second decade end year
    output_filesuffix : str
        Suffix for output filename
    
    Returns
    -------
    bool
        True if successful
    """
    warming_trends = _get_warming_trends()

    if os.path.exists(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix)):
        os.remove(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))
    
    climate_file = gdir.get_filepath('climate_historical')
    tmp = xr.open_dataset(climate_file)
    trend = warming_trends.loc[gdir.rgi_id] # from 1 to 12

    tmp['prcp'].loc[dict(time=slice(f'{first_dec_start}-01', f'{first_dec_end}-12'))] = \
                tmp['prcp'].sel(time=slice(f'{second_dec_start}-01', f'{second_dec_end}-12')).values
    tmp['temp'].loc[dict(time=slice(f'{first_dec_start}-01', f'{first_dec_end}-12'))] = \
                tmp['temp'].sel(time=slice(f'{second_dec_start}-01', f'{second_dec_end}-12')).values
    # elinimate warming trend for the first ten years
    for year in range(first_dec_start, first_dec_end+1):
        years_offset = year - first_dec_start
        
        for month in range(1, 13):
            month_data = tmp['temp'].sel(time=slice(f'{year}-{month:02d}-01', f'{year}-{month:02d}-28')).values
            warming_per_month = trend[month] * years_offset
            tmp['temp'].loc[dict(time=slice(f'{year}-{month:02d}-01', f'{year}-{month:02d}-28'))] = \
                month_data - warming_per_month

    tmp.to_netcdf(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))
    tmp.close()
    return True

@entity_task(log)
def clip_spinup_file_to_giventime(gdir, 
                                given_time=2000,
                                input_filesuffix='_spinup',
                                output_filesuffix=''):
    """Clip spinup files to a given time.
    
    This function truncates the model_diagnostics and model_geometry files
    from a spinup run to end at a specified year.
    
    Parameters
    ----------
    gdir : :py:class:`oggm.GlacierDirectory`
        the glacier directory to process
    given_time : int
        the year to clip the spinup files to (default: 2000)
    input_filesuffix : str
        filesuffix of the input spinup files (default: '_spinup')
    output_filesuffix : str
        filesuffix for the output clipped files (default: '')
    
    Returns
    -------
    bool
        True if successful
    """
    if os.path.exists(gdir.get_filepath('model_diagnostics', filesuffix=output_filesuffix)):
        os.remove(gdir.get_filepath('model_diagnostics', filesuffix=output_filesuffix))
    
    spinup_result_path = gdir.get_filepath('model_diagnostics', filesuffix=input_filesuffix)
    tmp = xr.open_dataset(spinup_result_path)
    tmp = tmp.sel(time=slice(1950, given_time)) # force the end time of spinup file is 2000
    tmp.to_netcdf(gdir.get_filepath('model_diagnostics', filesuffix=output_filesuffix)) # 1979 to 2000
    tmp.close()

    # 2. Clip model_geometry (complex, with groups) - 需要特殊处理
    input_geom_path = gdir.get_filepath('model_geometry', filesuffix=input_filesuffix)
    output_geom_path = gdir.get_filepath('model_geometry', filesuffix=output_filesuffix)
    
    if os.path.exists(output_geom_path):
        os.remove(output_geom_path)
    
    # 打开输入文件并复制所有 groups
    with Dataset(input_geom_path, 'r') as src:
        with Dataset(output_geom_path, 'w') as dst:
            # 复制全局属性
            dst.setncatts(src.__dict__)
            
            # 复制根 group 的维度
            for name, dimension in src.dimensions.items():
                if name == 'time':
                    # 找到 given_time 对应的索引
                    time_var = src.variables['time'][:]
                    time_idx = np.where(time_var <= given_time)[0]
                    if len(time_idx) > 0:
                        new_len = time_idx[-1] + 1
                    else:
                        new_len = 0
                    dst.createDimension(name, new_len if not dimension.isunlimited() else None)
                else:
                    dst.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
            
            # 复制根 group 的变量
            for name, variable in src.variables.items():
                out_var = dst.createVariable(name, variable.datatype, variable.dimensions)
                out_var.setncatts(variable.__dict__)
                if 'time' in variable.dimensions:
                    out_var[:] = variable[:new_len]
                else:
                    out_var[:] = variable[:]
            
            # 复制所有 flowline groups (fl_0, fl_1, ...)
            for grp_name in src.groups:
                src_grp = src.groups[grp_name]
                dst_grp = dst.createGroup(grp_name)
                
                # 复制 group 属性
                dst_grp.setncatts(src_grp.__dict__)
                
                # 复制 group 的维度
                for dim_name, dimension in src_grp.dimensions.items():
                    if dim_name == 'time':
                        dst_grp.createDimension(dim_name, new_len if not dimension.isunlimited() else None)
                    else:
                        dst_grp.createDimension(dim_name, len(dimension) if not dimension.isunlimited() else None)
                
                # 复制 group 的变量
                for var_name, variable in src_grp.variables.items():
                    out_var = dst_grp.createVariable(var_name, variable.datatype, variable.dimensions)
                    out_var.setncatts(variable.__dict__)
                    if 'time' in variable.dimensions:
                        out_var[:] = variable[:new_len]
                    else:
                        out_var[:] = variable[:]

    return True

@entity_task(log)
def remove_warming_trend(gdir, 
                        ys_start=2015,
                        ys_end=2024,
                        output_filesuffix=''):
    """
    Remove warming trend for certain year range.
    
    Returns
    -------
    bool
        True if successful
    """
    warming_trends = _get_warming_trends()

    if os.path.exists(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix)):
        os.remove(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))
    
    climate_file = gdir.get_filepath('climate_historical')
    tmp = xr.open_dataset(climate_file)
    trend = warming_trends.loc[gdir.rgi_id] # from 1 to 12

    # elinimate warming trend for the given year range
    for year in range(ys_start, ys_end+1):
        years_offset = year - ys_start
        for month in range(1, 13):
            month_data = tmp['temp'].sel(time=slice(f'{year}-{month:02d}-01', f'{year}-{month:02d}-28')).values
            warming_per_month = trend[month] * years_offset
            tmp['temp'].loc[dict(time=slice(f'{year}-{month:02d}-01', f'{year}-{month:02d}-28'))] = \
                month_data - warming_per_month

    tmp.to_netcdf(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))
    tmp.close()
    return True

@entity_task(log)
def add_warming_trend(gdir, 
                    ys_start=2005,
                    ys_end=2014,
                    output_filesuffix=''):
    """
    Add warming trend for certain year range.
    
    Returns
    -------
    bool
        True if successful
    """
    warming_trends = _get_warming_trends()

    if os.path.exists(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix)):
        os.remove(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))
                            
    climate_file = gdir.get_filepath('climate_historical')
    tmp = xr.open_dataset(climate_file)
    trend = warming_trends.loc[gdir.rgi_id] # from 1 to 12

    # add warming trend for the given year range
    for year in range(ys_start, ys_end+1):
        years_offset = year - ys_start
        
        for month in range(1, 13):
            month_data = tmp['temp'].sel(time=slice(f'{year}-{month:02d}-01', f'{year}-{month:02d}-28')).values
            warming_per_month = trend[month] * years_offset
            tmp['temp'].loc[dict(time=slice(f'{year}-{month:02d}-01', f'{year}-{month:02d}-28'))] = \
                month_data + warming_per_month

    tmp.to_netcdf(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))
    tmp.close()
    return True

@entity_task(log)
def extend_climate_with_baseline_period(gdir, 
                                      start_year=1950,
                                      end_year=2100,
                                      baseline_start_year=1991,
                                      baseline_end_year=2020,
                                      repeat_start_year=2025,
                                      output_filesuffix=''):
    """
    Extend climate data and fill with baseline period averages.
    
    Parameters
    ----------
    gdir : GlacierDirectory
        Glacier directory object
    start_year : int
        Starting year of the climate time series,
    end_year : int
        Ending year of the extended climate time series, 2100
    baseline_start_year : int
        Start year of baseline period for calculating monthly averages
    baseline_end_year : int
        End year of baseline period
    repeat_start_year : int
        Year when repetition of baseline averages begins, 2025
    output_filesuffix : str
        Suffix for output filename
    
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    if os.path.exists(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix)):
        os.remove(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))

    climate_file = gdir.get_filepath('climate_historical')
    tmp = xr.open_dataset(climate_file)
    
    # Calculate number of months for original and extended periods
    original_n_months = len(tmp['time'])
    new_n_months = (end_year - start_year + 1) * 12  # Total months from 1950 to 2100
    
    # Create new time coordinate
    new_time = pd.date_range(start=f'{start_year}-01-01', end=f'{end_year}-12-01', freq='MS')
    
    # Get baseline period data for calculating monthly averages
    reference_climate = tmp.sel(time=slice(f'{baseline_start_year}-01-01', f'{baseline_end_year}-12-01'))
    reference_climate_monthly_mean = reference_climate.groupby('time.month').mean(dim='time')
    
    # Get the dimensions of temp and prcp (excluding time dimension)
    temp_dims = list(tmp['temp'].dims)
    prcp_dims = list(tmp['prcp'].dims)
    
    # Get shape of data arrays (excluding time dimension)
    temp_other_dims = tmp['temp'].shape[1:]  # All dimensions except time
    prcp_other_dims = tmp['prcp'].shape[1:]  # All dimensions except time
    
    # Create new arrays with extended time dimension
    temp_new = np.zeros((new_n_months,) + temp_other_dims, dtype=tmp['temp'].dtype)
    prcp_new = np.zeros((new_n_months,) + prcp_other_dims, dtype=tmp['prcp'].dtype)
    
    # Copy original data to the beginning
    temp_new[:original_n_months] = tmp['temp'].values
    prcp_new[:original_n_months] = tmp['prcp'].values
    
    # Fill extended period (from repeat_start_year to end_year) with baseline monthly averages
    for year in range(repeat_start_year, end_year + 1):
        # Calculate index for this year's January
        year_start_idx = (year - start_year) * 12
        
        # Fill 12 months with baseline monthly averages
        for month in range(1, 13):
            month_idx = year_start_idx + (month - 1)
            if month_idx < new_n_months:
                # Get baseline average for this month
                temp_new[month_idx] = reference_climate_monthly_mean['temp'].sel(month=month).values
                prcp_new[month_idx] = reference_climate_monthly_mean['prcp'].sel(month=month).values
    
    # Create new dataset
    new_tmp = xr.Dataset()
    
    # Copy all dimensions except time
    for dim_name, dim_size in tmp.dims.items():
        if dim_name != 'time':
            new_tmp[dim_name] = tmp[dim_name]
    
    # Set new time coordinate
    new_tmp['time'] = xr.DataArray(new_time, dims='time')
    
    # Set new temp and prcp data arrays
    new_tmp['temp'] = xr.DataArray(temp_new, dims=temp_dims)
    new_tmp['prcp'] = xr.DataArray(prcp_new, dims=prcp_dims)
    
    # Copy attributes
    new_tmp['temp'].attrs = tmp['temp'].attrs.copy()
    new_tmp['prcp'].attrs = tmp['prcp'].attrs.copy()
    new_tmp.attrs = tmp.attrs.copy()
    new_tmp.attrs['hydro_yr_1'] = end_year
    
    # Save output
    new_tmp.to_netcdf(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))
    
    tmp.close()
    new_tmp.close()
    
    return True

@entity_task(log)
def extend_climate_loop(gdir, 
                        start_year=1950,
                        end_year=2100,
                        loop_start_year=2011,
                        loop_end_year=2020,
                        repeat_start_year=2025,
                        output_filesuffix=''):
    """
    Extend climate data and fill with looping period data.
    
    Parameters
    ----------
    gdir : GlacierDirectory
        Glacier directory object
    start_year : int
        Starting year of the climate time series
    end_year : int
        Ending year of the extended climate time series
    loop_start_year : int
        Start year of loop period (data to be repeated)
    loop_end_year : int
        End year of loop period (inclusive)
    repeat_start_year : int
        Year when looping begins
    output_filesuffix : str
        Suffix for output filename
    
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    if os.path.exists(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix)):
        os.remove(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))

    climate_file = gdir.get_filepath('climate_historical')
    tmp = xr.open_dataset(climate_file)
    
    # Calculate number of months for original and extended periods
    original_n_months = len(tmp['time'])
    new_n_months = (end_year - start_year + 1) * 12  # Total months from 1950 to 2100
    
    # Create new time coordinate
    new_time = pd.date_range(start=f'{start_year}-01-01', end=f'{end_year}-12-01', freq='MS')
    
    # Get baseline period data for calculating monthly averages
    loop_climate = tmp.sel(time=slice(f'{loop_start_year}-01-01', f'{loop_end_year}-12-01'))
    loop_n_months = len(loop_climate['time'])

    # Get the dimensions of temp and prcp (excluding time dimension)
    temp_dims = list(tmp['temp'].dims)
    prcp_dims = list(tmp['prcp'].dims)
    
    # Get shape of data arrays (excluding time dimension)
    temp_other_dims = tmp['temp'].shape[1:]  # All dimensions except time
    prcp_other_dims = tmp['prcp'].shape[1:]  # All dimensions except time
    
    # Create new arrays with extended time dimension
    temp_new = np.zeros((new_n_months,) + temp_other_dims, dtype=tmp['temp'].dtype)
    prcp_new = np.zeros((new_n_months,) + prcp_other_dims, dtype=tmp['prcp'].dtype)
    
    # Copy original data to the beginning
    temp_new[:original_n_months] = tmp['temp'].values
    prcp_new[:original_n_months] = tmp['prcp'].values

    repeat_start_idx = (repeat_start_year - start_year) * 12

    # Fill extended period (from repeat_start_year to end_year) with loop monthly averages
    for month_idx in range(repeat_start_idx, new_n_months):
        # Calculate which year and month this index represents
        year = start_year + month_idx // 12
        month = (month_idx % 12) + 1  # 1-12
        
        # Calculate corresponding month in loop period
        # Find the position of this month within the loop period
        loop_month_idx = (month_idx - repeat_start_idx) % loop_n_months
        
        # Get loop period data for this month
        temp_new[month_idx] = loop_climate['temp'].values[loop_month_idx]
        prcp_new[month_idx] = loop_climate['prcp'].values[loop_month_idx]
    
    # Create new dataset
    new_tmp = xr.Dataset()
    
    # Copy all dimensions except time
    for dim_name, dim_size in tmp.dims.items():
        if dim_name != 'time':
            new_tmp[dim_name] = tmp[dim_name]
    
    # Set new time coordinate
    new_tmp['time'] = xr.DataArray(new_time, dims='time')
    
    # Set new temp and prcp data arrays
    new_tmp['temp'] = xr.DataArray(temp_new, dims=temp_dims)
    new_tmp['prcp'] = xr.DataArray(prcp_new, dims=prcp_dims)
    
    # Copy attributes
    new_tmp['temp'].attrs = tmp['temp'].attrs.copy()
    new_tmp['prcp'].attrs = tmp['prcp'].attrs.copy()
    new_tmp.attrs = tmp.attrs.copy()
    new_tmp.attrs['hydro_yr_1'] = end_year
    
    # Save output
    new_tmp.to_netcdf(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))
    
    tmp.close()
    new_tmp.close()
    
    return True

@entity_task(log)
def extend_climate_loop_with_warming(gdir, 
                                start_year=1950,
                                end_year=2100,
                                loop_start_year=2011,
                                loop_end_year=2020,
                                repeat_start_year=2025,
                                output_filesuffix=''):
    """
    Extend climate data, fill with looping period data, and apply warming trend.
    
    Parameters
    ----------
    gdir : GlacierDirectory
        Glacier directory object
    start_year : int
        Starting year
    end_year : int
        Ending year
    loop_start_year : int
        Start year of loop period
    loop_end_year : int
        End year of loop period (inclusive)
    repeat_start_year : int
        Year when looping begins
    output_filesuffix : str
        Suffix for output filename
    
    Returns
    -------
    bool
        True if successful
    """
    warming_trends = _get_warming_trends()

    if os.path.exists(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix)):
        os.remove(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))

    climate_file = gdir.get_filepath('climate_historical')
    tmp = xr.open_dataset(climate_file)
    trend = warming_trends.loc[gdir.rgi_id]  # from 1 to 12

    # Calculate number of months for original and extended periods
    original_n_months = len(tmp['time'])
    new_n_months = (end_year - start_year + 1) * 12  # Total months from 1950 to 2100
    
    # Create new time coordinate
    new_time = pd.date_range(start=f'{start_year}-01-01', end=f'{end_year}-12-01', freq='MS')
    
    # Get baseline period data for calculating monthly averages
    loop_climate = tmp.sel(time=slice(f'{loop_start_year}-01-01', f'{loop_end_year}-12-01'))
    loop_n_months = len(loop_climate['time'])

    # Get the dimensions of temp and prcp (excluding time dimension)
    temp_dims = list(tmp['temp'].dims)
    prcp_dims = list(tmp['prcp'].dims)
    
    # Get shape of data arrays (excluding time dimension)
    temp_other_dims = tmp['temp'].shape[1:]  # All dimensions except time
    prcp_other_dims = tmp['prcp'].shape[1:]  # All dimensions except time
    
    # Create new arrays with extended time dimension
    temp_new = np.zeros((new_n_months,) + temp_other_dims, dtype=tmp['temp'].dtype)
    prcp_new = np.zeros((new_n_months,) + prcp_other_dims, dtype=tmp['prcp'].dtype)
    
    # Copy original data to the beginning
    temp_new[:original_n_months] = tmp['temp'].values
    prcp_new[:original_n_months] = tmp['prcp'].values

    repeat_start_idx = (repeat_start_year - start_year) * 12

    # Fill extended period (from repeat_start_year to end_year) with loop monthly averages
    for month_idx in range(repeat_start_idx, new_n_months):
        # Calculate which year and month this index represents
        year = start_year + month_idx // 12
        month = (month_idx % 12) + 1  # 1-12
        offset_years = year - repeat_start_year
        
        # Calculate corresponding month in loop period
        # Find the position of this month within the loop period
        loop_month_idx = (month_idx - repeat_start_idx) % loop_n_months
        
        # Get loop period data for this month
        temp_new[month_idx] = loop_climate['temp'].values[loop_month_idx] + trend[month] * offset_years
        prcp_new[month_idx] = loop_climate['prcp'].values[loop_month_idx]
    
    # Create new dataset
    new_tmp = xr.Dataset()
    
    # Copy all dimensions except time
    for dim_name, dim_size in tmp.dims.items():
        if dim_name != 'time':
            new_tmp[dim_name] = tmp[dim_name]
    
    # Set new time coordinate
    new_tmp['time'] = xr.DataArray(new_time, dims='time')
    
    # Set new temp and prcp data arrays
    new_tmp['temp'] = xr.DataArray(temp_new, dims=temp_dims)
    new_tmp['prcp'] = xr.DataArray(prcp_new, dims=prcp_dims)
    
    # Copy attributes
    new_tmp['temp'].attrs = tmp['temp'].attrs.copy()
    new_tmp['prcp'].attrs = tmp['prcp'].attrs.copy()
    new_tmp.attrs = tmp.attrs.copy()
    new_tmp.attrs['hydro_yr_1'] = end_year
    
    # Save output
    new_tmp.to_netcdf(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))
    
    tmp.close()
    new_tmp.close()
    
    return True

@entity_task(log)
def extend_climate_with_baseline_warming(gdir, 
                                      start_year=1950,
                                      end_year=2100,
                                      baseline_start_year=1991,
                                      baseline_end_year=2020,
                                      repeat_start_year=2025,
                                      output_filesuffix=''):
    """
    Extend climate data, fill with baseline period averages, and apply warming trend.
    
    Parameters
    ----------
    gdir : GlacierDirectory
        Glacier directory object
    start_year : int
        Starting year
    end_year : int
        Ending year
    baseline_start_year : int
        Start year of baseline period
    baseline_end_year : int
        End year of baseline period (inclusive)
    repeat_start_year : int
        Year when repetition begins
    output_filesuffix : str
        Suffix for output filename
    
    Returns
    -------
    bool
        True if successful
    """
    warming_trends = _get_warming_trends()

    if os.path.exists(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix)):
        os.remove(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))
    
    climate_file = gdir.get_filepath('climate_historical')
    tmp = xr.open_dataset(climate_file)
    trend = warming_trends.loc[gdir.rgi_id] # from 1 to 12
    
    # Calculate number of months for original and extended periods
    original_n_months = len(tmp['time'])
    new_n_months = (end_year - start_year + 1) * 12  # Total months from 1950 to 2100
    
    # Create new time coordinate
    new_time = pd.date_range(start=f'{start_year}-01-01', end=f'{end_year}-12-01', freq='MS')
    
    # Get baseline period data for calculating monthly averages
    reference_climate = tmp.sel(time=slice(f'{baseline_start_year}-01-01', f'{baseline_end_year}-12-01'))
    reference_climate_monthly_mean = reference_climate.groupby('time.month').mean(dim='time')
    
    # Get the dimensions of temp and prcp (excluding time dimension)
    temp_dims = list(tmp['temp'].dims)
    prcp_dims = list(tmp['prcp'].dims)
    
    # Get shape of data arrays (excluding time dimension)
    temp_other_dims = tmp['temp'].shape[1:]  # All dimensions except time
    prcp_other_dims = tmp['prcp'].shape[1:]  # All dimensions except time
    
    # Create new arrays with extended time dimension
    temp_new = np.zeros((new_n_months,) + temp_other_dims, dtype=tmp['temp'].dtype)
    prcp_new = np.zeros((new_n_months,) + prcp_other_dims, dtype=tmp['prcp'].dtype)
    
    # Copy original data to the beginning
    temp_new[:original_n_months] = tmp['temp'].values
    prcp_new[:original_n_months] = tmp['prcp'].values
    
    # Fill extended period (from repeat_start_year to end_year) with baseline monthly averages
    for year in range(repeat_start_year, end_year + 1):
        # Calculate index for this year's January
        year_start_idx = (year - start_year) * 12
        years_offset = year - repeat_start_year
        
        # Fill 12 months with baseline monthly averages
        for month in range(1, 13):
            month_idx = year_start_idx + (month - 1)
            if month_idx < new_n_months:
                # Get baseline average for this month
                temp_new[month_idx] = reference_climate_monthly_mean['temp'].sel(month=month).values + trend[month] * years_offset
                prcp_new[month_idx] = reference_climate_monthly_mean['prcp'].sel(month=month).values
    
    # Create new dataset
    new_tmp = xr.Dataset()
    
    # Copy all dimensions except time
    for dim_name, dim_size in tmp.dims.items():
        if dim_name != 'time':
            new_tmp[dim_name] = tmp[dim_name]
    
    # Set new time coordinate
    new_tmp['time'] = xr.DataArray(new_time, dims='time')
    
    # Set new temp and prcp data arrays
    new_tmp['temp'] = xr.DataArray(temp_new, dims=temp_dims)
    new_tmp['prcp'] = xr.DataArray(prcp_new, dims=prcp_dims)
    
    # Copy attributes
    new_tmp['temp'].attrs = tmp['temp'].attrs.copy()
    new_tmp['prcp'].attrs = tmp['prcp'].attrs.copy()
    new_tmp.attrs = tmp.attrs.copy()
    new_tmp.attrs['hydro_yr_1'] = end_year
    
    # Save output
    new_tmp.to_netcdf(gdir.get_filepath('climate_historical', filesuffix=output_filesuffix))
    
    tmp.close()
    new_tmp.close()
    
    return True

@entity_task(log)
def simulate_future_extremes_window_scaling(gdir, 
                                        start_year=1950,
                                        end_year=2100,
                                        ref_start_year=2010,
                                        ref_end_year=2024, # inclusive
                                        future_start_year=2025,
                                        scale_factor=None, #0-0.5
                                        output_filesuffix=''):
    """
    Repulicate extremes during 2010-2024 to future period.
    
    Parameters
    ----------
    gdir : GlacierDirectory
        Glacier directory object
    start_year : int
        Starting year for simulation
    end_year : int
        Ending year for simulation
    ref_start_year : int
        Start year of reference period (for repulicating extremes)
    ref_end_year : int
        End year of reference period (inclusive, for repulicating extremes)
    future_start_year : int
        Start year of future period
    output_filesuffix : str
        Suffix for output filename
    
    Returns
    -------
    bool
        True if successful
    """
    from statsmodels.tsa.seasonal import STL
    # CMIP6 scenarios to analyze
    cmip6_scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp585']
    cmip6_models = ['BCC-CSM2-MR', 'CAMS-CSM1-0', 'CESM2', 'CESM2-WACCM', 
                    'EC-Earth3', 'EC-Earth3-Veg', 'FGOALS-f3-L', 'GFDL-ESM4',
                    'INM-CM4-8', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NorESM2-MM']
    
    if os.path.exists(os.path.join(gdir.dir, 'climate_historical.nc')):
        ref_era5_path = os.path.join(gdir.dir, 'climate_historical.nc')
        ref_era5 = xr.open_dataset(ref_era5_path).sel(time=slice(f"{start_year}-01-01", f"{future_start_year-1}-12-01"))
    else:
        raise FileNotFoundError(f"ERA5 climate file not found for {gdir.rgi_id}")

    hist_date_range = pd.date_range(start=f"01/01/{start_year}", end=f"12/01/{future_start_year-1}", freq="MS")
    ref_era5['time'] = hist_date_range

    # adjust_cmip = []
    # rolling window in the future
    selected_date_range = pd.date_range(start=f"01/01/{ref_start_year}", end=f"12/01/{ref_end_year}", freq="MS")
    window_length = len(selected_date_range)
    
    # STL for long-term warming trend from ERA5 Data
    ref_era5_df = ref_era5.to_dataframe()
    prcp_trend = STL(ref_era5_df['prcp'], period=12, trend=241, seasonal=13).fit().trend # long-term warming trend
    temp_trend = STL(ref_era5_df['temp'], period=12, trend=241, seasonal=13).fit().trend # long-term warming trend

    selected_prcp_fluc = (ref_era5_df['prcp'] - prcp_trend).loc[selected_date_range]
    selected_temp_fluc = (ref_era5_df['temp'] - temp_trend).loc[selected_date_range]

    for GCM in cmip6_models:  # loop through all GCMs
        for SSP in cmip6_scenarios:  # loop through all SSPs
            cmip_path = os.path.join(gdir.dir, 'gcm_data_{}_{}.nc'.format(GCM, SSP))
            target_cmip = xr.open_dataset(cmip_path).sel(time=slice(f"{start_year}-01-01", f"{end_year}-12-31"))

            target_cmip['time'] = pd.date_range(start=f"01/01/{start_year}", end=target_cmip.time.max().values, freq="MS")
            target_cmip_df = target_cmip.to_dataframe() # structure same to selected_era5，2001-01-01
            target_cmip_df['year'] = target_cmip_df.index.year
            # STL
            cmip_prcp_trend = STL(target_cmip_df['prcp'], period=12, trend=241, seasonal=13).fit().trend
            cmip_temp_trend = STL(target_cmip_df['temp'], period=12, trend=241, seasonal=13).fit().trend
            # fluc
            selected_cmip_prcp_fluc = (target_cmip_df['prcp'] - cmip_prcp_trend).loc[selected_date_range]
            selected_cmip_temp_fluc = (target_cmip_df['temp'] - cmip_temp_trend).loc[selected_date_range]

            ratio_prcp = selected_prcp_fluc/selected_cmip_prcp_fluc
            ratio_temp = selected_temp_fluc/selected_cmip_temp_fluc

            ratio_prcp = np.where(ratio_prcp<0, 1, ratio_prcp)
            ratio_prcp = np.where(ratio_prcp>2, 2, ratio_prcp)
            ratio_temp = np.where(ratio_temp<0, 1, ratio_temp)
            ratio_temp = np.where(ratio_temp>2, 2, ratio_temp) # [0, 2]

            # apply ratio to CMIP6 data during selected range and roll future
            adjusted_prcp = target_cmip_df['prcp'].copy()
            adjusted_temp = target_cmip_df['temp'].copy()
            # selected range
            adjusted_prcp.loc[selected_date_range] = cmip_prcp_trend.loc[selected_date_range] + \
                        selected_cmip_prcp_fluc * ratio_prcp
            adjusted_temp.loc[selected_date_range] = cmip_temp_trend.loc[selected_date_range] + \
                        selected_cmip_temp_fluc * ratio_temp

            # future
            future_start = pd.Timestamp(f"{future_start_year}-01-01")
            future_end = target_cmip['time'].max().values #2099-12-1
            current_start = future_start #2025/1/1
            while current_start <= future_end:
                current_end = current_start + pd.DateOffset(months=window_length-1)
                if current_end > future_end:
                    current_end = future_end

                current_window = pd.date_range(start=current_start, end=current_end, freq='MS')
                current_cmip_prcp_trend = cmip_prcp_trend.loc[current_window]
                current_cmip_temp_trend = cmip_temp_trend.loc[current_window]

                current_cmip_prcp_fluc = (target_cmip_df['prcp'] - cmip_prcp_trend).loc[current_window] # within current window
                current_cmip_temp_fluc = (target_cmip_df['temp'] - cmip_temp_trend).loc[current_window]

                # trend + fluc * ratio
                if scale_factor is None:
                    adjusted_prcp.loc[current_window] = current_cmip_prcp_trend + current_cmip_prcp_fluc \
                                                        * ratio_prcp[:len(current_window)]
                    adjusted_temp.loc[current_window] = current_cmip_temp_trend + current_cmip_temp_fluc \
                                                        * ratio_temp[:len(current_window)]
                else:
                    # scale ratio_prcp and ratio_temp
                    # current_cmip_prcp_fluc>0, need to enlarge; current_cmip_prcp_fluc<0, need to shrink, for temp
                    ratio_temp_window = ratio_temp[:len(current_window)] # match current window length
                    ratio_prcp_window = ratio_prcp[:len(current_window)]
                    dynamic_scale_factor = (target_cmip_df.loc[current_window, 'year'].values - future_start_year)/(end_year - future_start_year)*scale_factor # array
                    scaled_ratio_temp = np.where(current_cmip_temp_fluc>0, ratio_temp_window * (1+dynamic_scale_factor), ratio_temp_window * (1-dynamic_scale_factor)) # array
                    scaled_ratio_prcp = np.where(current_cmip_prcp_fluc>0, ratio_prcp_window * (1-dynamic_scale_factor), ratio_prcp_window * (1+dynamic_scale_factor)) # array

                    adjusted_prcp.loc[current_window] = current_cmip_prcp_trend + current_cmip_prcp_fluc \
                                                        * scaled_ratio_prcp
                    adjusted_temp.loc[current_window] = current_cmip_temp_trend + current_cmip_temp_fluc \
                                                        * scaled_ratio_temp
                # move to next month
                current_start = current_end + pd.DateOffset(months=1)

            adjusted_prcp = adjusted_prcp.mask(adjusted_prcp<0, 0)
            adjusted_df = pd.DataFrame({
                'prcp': adjusted_prcp,
                'temp': adjusted_temp
            })
            
            # turn to nc
            adjusted_ds = target_cmip.copy(deep=False)
            
            # 更新 prcp 和 temp 的值
            adjusted_ds['prcp'] = xr.DataArray(
                adjusted_df['prcp'].values,
                dims=target_cmip.prcp.dims,
                coords={dim: target_cmip.coords[dim] for dim in target_cmip.prcp.dims},
                attrs=target_cmip.prcp.attrs if hasattr(target_cmip.prcp, 'attrs') else {}
            )
            adjusted_ds['temp'] = xr.DataArray(
                adjusted_df['temp'].values,
                dims=target_cmip.temp.dims,
                coords={dim: target_cmip.coords[dim] for dim in target_cmip.temp.dims},
                attrs=target_cmip.temp.attrs if hasattr(target_cmip.temp, 'attrs') else {}
            )
            adjusted_ds['temp'].loc[{'time': hist_date_range}] = ref_era5['temp'].values
            adjusted_ds['prcp'].loc[{'time': hist_date_range}] = ref_era5['prcp'].values
            adjusted_ds.attrs['ref_hgt'] = ref_era5.attrs['ref_hgt']
            adjusted_ds.attrs['ref_pix_lon'] = ref_era5.attrs['ref_pix_lon']
            adjusted_ds.attrs['ref_pix_lat'] = ref_era5.attrs['ref_pix_lat']
            adjusted_ds.attrs['ref_pix_dis'] = ref_era5.attrs['ref_pix_dis']
            adjusted_ds.attrs['valid_future_years'] = len(adjusted_ds['time'])/12 - 75 # 75 years from 2025 to 2099

            if np.any(np.isnan(adjusted_ds.temp.values)) or np.any(np.isnan(adjusted_ds.prcp.values)):
                raise ValueError(f"NaN values found in {GCM}_{SSP} for {gdir.rgi_id}")
            # dims_to_expand = []
            # for coord_name in ['rgi_id', 'GCM', 'SSP']:
            #     if coord_name in adjusted_ds.coords and coord_name not in adjusted_ds.dims:
            #         dims_to_expand.append(coord_name)
            
            # if dims_to_expand:
            #     adjusted_ds = adjusted_ds.expand_dims(dims_to_expand)
            out_path = os.path.join(gdir.dir, 'gcm_data_{}_{}{}.nc'.format(GCM, SSP, output_filesuffix))
            adjusted_ds.to_netcdf(out_path)
            adjusted_ds.close()

    return True

@entity_task(log)
def repulicate_cmip6(gdir, 
                    start_year=1950,
                    end_year=2100,
                    future_start_year=2025,
                    output_filesuffix=''):
    """
    Redesign the time index of original CMIP6 and update the historical data using ERA5.
    
    Parameters
    ----------
    gdir : GlacierDirectory
        Glacier directory object
    start_year : int
        Starting year for simulation
    end_year : int
        Ending year for simulation
    future_start_year : int
        Start year of future period
    output_filesuffix : str
        Suffix for output filename
    
    Returns
    -------
    bool
        True if successful
    """
    # from statsmodels.tsa.seasonal import STL
    # CMIP6 scenarios to analyze
    cmip6_scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp585']
    cmip6_models = ['BCC-CSM2-MR', 'CAMS-CSM1-0', 'CESM2', 'CESM2-WACCM', 
                    'EC-Earth3', 'EC-Earth3-Veg', 'FGOALS-f3-L', 'GFDL-ESM4',
                    'INM-CM4-8', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NorESM2-MM']
    
    # if os.path.exists(os.path.join(gdir.dir, 'climate_historical.nc')):
    #     ref_era5_path = os.path.join(gdir.dir, 'climate_historical.nc')
    #     ref_era5 = xr.open_dataset(ref_era5_path).sel(time=slice("1950-01-01", "2024-12-01"))
    # else:
    #     raise FileNotFoundError(f"ERA5 climate file not found for {gdir.rgi_id}")

    # ref_era5['time'] = pd.date_range(start="01/01/1950", end="12/01/2024", freq="MS")
    # ref_era5_df = ref_era5.to_dataframe()
    # hist data of CMIP6 should align with ERA5
    ds_ref_path = os.path.join(gdir.dir, 'climate_historical.nc')
    ds_ref = xr.open_dataset(ds_ref_path).sel(time=slice(f"{start_year}-01-01", f"{future_start_year-1}-12-31"))
    hist_time_range = pd.date_range(start=f"01/01/{start_year}", end=f"12/01/{future_start_year-1}", freq="MS")
    ds_ref['time'] = hist_time_range

    for GCM in cmip6_models:  # loop through all GCMs
        for SSP in cmip6_scenarios:  # loop through all SSPs
            cmip_path = os.path.join(gdir.dir, 'gcm_data_{}_{}.nc'.format(GCM, SSP))
            target_cmip = xr.open_dataset(cmip_path).sel(time=slice(f"{start_year}-01-01", f"{end_year}-12-31"))

            target_cmip['time'] = pd.date_range(start=f"01/01/{start_year}", end=target_cmip.time.max().values, freq="MS")
            target_cmip['temp'].loc[{'time': hist_time_range}] = ds_ref['temp'].values
            target_cmip['prcp'].loc[{'time': hist_time_range}] = ds_ref['prcp'].values

            target_cmip.attrs['ref_hgt'] = ds_ref.attrs['ref_hgt']
            target_cmip.attrs['ref_pix_lon'] = ds_ref.attrs['ref_pix_lon']
            target_cmip.attrs['ref_pix_lat'] = ds_ref.attrs['ref_pix_lat']
            target_cmip.attrs['ref_pix_dis'] = ds_ref.attrs['ref_pix_dis']
            target_cmip.attrs['valid_future_years'] = len(target_cmip['time'])/12 - 75 # 75 years from 2025 to 2099

            out_path = os.path.join(gdir.dir, 'gcm_data_{}_{}{}.nc'.format(GCM, SSP, output_filesuffix))
            target_cmip.to_netcdf(out_path)
            target_cmip.close()
    return True


@entity_task(log)
def simulate_future_extremes_rank(gdir, 
                                    start_year=1950,
                                    end_year=2100,
                                    future_start_year=2025,
                                    future_cooling_factor=None, #>1, for sensitive test
                                    output_filesuffix=''):
    """
    Rank future years based on scores derived from CMIP6 climate projections.
    The future rank corresponds to historical ranks derived from ERA5 data.
    The same rank indicates the same scaling upon average-aboved temperatures and average-below precipitation.
    Parameters
    ----------
    gdir : GlacierDirectory
        Glacier directory object
    start_year : int
        Starting year for simulation
    end_year : int
        Ending year for simulation
    future_start_year : int
        Start year of future period
    output_filesuffix : str
        Suffix for output filename
    
    Returns
    -------
    bool
        True if successful
    """
    from statsmodels.tsa.seasonal import STL
    # CMIP6 scenarios to analyze
    cmip6_scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp585']
    cmip6_models = ['BCC-CSM2-MR', 'CAMS-CSM1-0', 'CESM2', 'CESM2-WACCM', 
                    'EC-Earth3', 'EC-Earth3-Veg', 'FGOALS-f3-L', 'GFDL-ESM4',
                    'INM-CM4-8', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NorESM2-MM']
    
    if os.path.exists(os.path.join(gdir.dir, 'climate_historical.nc')):
        ref_era5_path = os.path.join(gdir.dir, 'climate_historical.nc')
        ref_era5 = xr.open_dataset(ref_era5_path).sel(time=slice(f"{start_year}-01-01", f"{future_start_year-1}-12-01"))
    else:
        raise FileNotFoundError(f"ERA5 climate file not found for {gdir.rgi_id}")

    hist_date_range = pd.date_range(start=f"01/01/{start_year}", end=f"12/01/{future_start_year-1}", freq="MS")
    ref_era5['time'] = hist_date_range

    # adjust_cmip = []
    era5_temp_scores = _ERA5_hist_scores(sheet_name='Temperature_Extreme_Scores')
    era5_prcp_scores = _ERA5_hist_scores(sheet_name='Precipitation_Extreme_Scores')
    era5_temp_scores = era5_temp_scores.loc[[gdir.rgi_id], :].T
    era5_temp_scores.columns = ['temp_scores'] # from 1990 to 2024
    era5_temp_scores.index = era5_temp_scores.index.astype(int)
    era5_prcp_scores = era5_prcp_scores.loc[[gdir.rgi_id], :].T
    era5_prcp_scores.columns = ['prcp_scores']
    era5_prcp_scores.index = era5_prcp_scores.index.astype(int)
    # cal the rank
    era5_temp_scores_reversed = era5_temp_scores[::-1].copy() # from 2024 to 1990
    era5_temp_scores['temp_rank'] = era5_temp_scores_reversed['temp_scores'].rank(
        ascending=False, method='first'
    )
    era5_prcp_scores_reversed = era5_prcp_scores[::-1].copy() # from 2024 to 1990
    era5_prcp_scores['prcp_rank'] = era5_prcp_scores_reversed['prcp_scores'].rank(
        ascending=False, method='first'
    ) # unique rank, 1-35
    
    # STL for long-term warming trend from ERA5 Data
    ref_era5_df = ref_era5.to_dataframe()
    # prcp_trend = STL(ref_era5_df['prcp'], period=12, trend=241, seasonal=13).fit().trend # long-term warming trend
    # temp_trend = STL(ref_era5_df['temp'], period=12, trend=241, seasonal=13).fit().trend # long-term warming trend

    annual_mean_era5_temp = ref_era5_df['temp'].resample('AS').mean()
    annual_mean_era5_prcp = ref_era5_df['prcp'].resample('AS').mean()

    for GCM in cmip6_models:  # loop through all GCMs
        for SSP in cmip6_scenarios:  # loop through all SSPs
            cmip_path = os.path.join(gdir.dir, 'gcm_data_{}_{}.nc'.format(GCM, SSP))
            target_cmip = xr.open_dataset(cmip_path).sel(time=slice(f"{start_year}-01-01", f"{end_year}-12-31"))

            target_cmip['time'] = pd.date_range(start=f"01/01/{start_year}", end=target_cmip.time.max().values, freq="MS")
            target_cmip_df = target_cmip.to_dataframe() # structure same to selected_era5，2001-01-01
            target_cmip_df['year'] = target_cmip_df.index.year
            # score for future years
            future_cmip_df = target_cmip_df.loc[target_cmip_df['year']>=future_start_year] #monthly data
            future_annual_df = future_cmip_df.groupby('year')
            future_scores = pd.DataFrame(index=future_annual_df.groups.keys(), columns=['temp', 'prcp'])

            # calculate critical thresholds for future years
            temp_99 = np.nanpercentile(future_cmip_df['temp'], 99)
            temp_95 = np.nanpercentile(future_cmip_df['temp'], 95)
            temp_90 = np.nanpercentile(future_cmip_df['temp'], 90)
            prcp_01 = np.nanpercentile(future_cmip_df['prcp'], 1)
            prcp_05 = np.nanpercentile(future_cmip_df['prcp'], 5)
            prcp_10 = np.nanpercentile(future_cmip_df['prcp'], 10)

            for year, year_data in future_annual_df:
                year_score_temp = 4 * (year_data['temp'] >= temp_99).sum() + \
                                    2 * ((year_data['temp'] >= temp_95) & (year_data['temp'] < temp_99)).sum() + \
                                    ((year_data['temp'] >= temp_90) & (year_data['temp'] < temp_95)).sum()
                year_score_prcp = 4 * (year_data['prcp'] <= prcp_01).sum() + \
                                    2 * ((year_data['prcp'] <= prcp_05) & (year_data['prcp'] > prcp_01)).sum() + \
                                    ((year_data['prcp'] <= prcp_10) & (year_data['prcp'] > prcp_05)).sum()
                future_scores.loc[year, 'temp'] = year_score_temp
                future_scores.loc[year, 'prcp'] = year_score_prcp
            # # STL
            # cmip_prcp_trend = STL(target_cmip_df['prcp'], period=12, trend=241, seasonal=13).fit().trend
            # cmip_temp_trend = STL(target_cmip_df['temp'], period=12, trend=241, seasonal=13).fit().trend

            # cal future rank
            future_scores_reversed = future_scores[::-1].copy() # from 2100 to 2025
            future_scores.loc[2025:2094, 'temp_rank'] = future_scores_reversed.loc[2094:2025, 'temp'].rank(
                ascending=False, method='first'
            )
            future_scores.loc[2095:, 'temp_rank'] = future_scores_reversed['temp'].rank(
                ascending=False, method='first'
            ).loc[:2095]
            future_scores.loc[2025:2094, 'prcp_rank'] = future_scores_reversed.loc[2094:2025, 'prcp'].rank(
                ascending=False, method='first'
            )
            future_scores.loc[2095:, 'prcp_rank'] = future_scores_reversed['prcp'].rank(
                ascending=False, method='first'
            ).loc[:2095] # reserved

            future_scores['temp_rank'] = np.where(future_scores['temp_rank']>70, 70, future_scores['temp_rank'])
            future_scores['prcp_rank'] = np.where(future_scores['prcp_rank']>70, 70, future_scores['prcp_rank'])
            # cal scaling factor
            future_scores['temp_scaling'] = 1
            future_scores['prcp_scaling'] = 1

            for year, year_info in future_scores.iterrows():
                year_temp_rank = ceil(year_info['temp_rank']/2)
                hist_year_mapped = era5_temp_scores[era5_temp_scores['temp_rank']==year_temp_rank].index[0] # year
                hist_year_temp = ref_era5_df.loc[ref_era5_df.index.year==hist_year_mapped, 'temp']
                hist_year_temp_anomaly = hist_year_temp - hist_year_temp.mean()
                hist_year_temp_culmutive = hist_year_temp_anomaly[hist_year_temp_anomaly>=0].sum() #>0

                hist_cmip_temp = target_cmip_df.loc[target_cmip_df.index.year==hist_year_mapped, 'temp']
                hist_cmip_temp_anomaly = hist_cmip_temp - hist_cmip_temp.mean() # monthly data
                hist_cmip_temp_culmutive = hist_cmip_temp_anomaly[hist_cmip_temp_anomaly>=0].sum()
                temp_scaling = hist_year_temp_culmutive / hist_cmip_temp_culmutive # ERA5/CMIP6, only anomalies for >0

                if temp_scaling > 2:
                    temp_scaling = 2
                future_scores.loc[year, 'temp_scaling'] = temp_scaling

                year_prcp_rank = ceil(year_info['prcp_rank']/2)
                hist_year_mapped = era5_prcp_scores[era5_prcp_scores['prcp_rank']==year_prcp_rank].index[0] # year
                hist_year_prcp = ref_era5_df.loc[ref_era5_df.index.year==hist_year_mapped, 'prcp']
                hist_year_prcp_anomaly = hist_year_prcp - hist_year_prcp.mean()
                hist_year_prcp_culmutive = hist_year_prcp_anomaly[hist_year_prcp_anomaly<=0].sum() #<0

                hist_cmip_prcp = target_cmip_df.loc[target_cmip_df.index.year==hist_year_mapped, 'prcp']
                hist_cmip_prcp_anomaly = hist_cmip_prcp - hist_cmip_prcp.mean() # monthly data
                hist_cmip_prcp_culmutive = hist_cmip_prcp_anomaly[hist_cmip_prcp_anomaly<=0].sum()
                prcp_scaling = hist_year_prcp_culmutive / hist_cmip_prcp_culmutive # ERA5/CMIP6, only anomalies for <0

                if prcp_scaling > 2:
                    prcp_scaling = 2
                future_scores.loc[year, 'prcp_scaling'] = prcp_scaling

            # apply ratio to CMIP6 data during selected range and roll future
            adjusted_prcp = target_cmip_df['prcp'].copy()
            adjusted_temp = target_cmip_df['temp'].copy()

            # future adjustment
            for year, year_info in future_scores.iterrows():
                cmip_year_temp = adjusted_temp.loc[adjusted_temp.index.year==year]
                cmip_year_temp_anomaly = cmip_year_temp - cmip_year_temp.mean()
                cmip_year_temp_adjusted = np.where(cmip_year_temp_anomaly>=0, 
                                                    cmip_year_temp.mean()+cmip_year_temp_anomaly*year_info['temp_scaling'],
                                                     cmip_year_temp)
                
                if future_cooling_factor is not None:
                    cmip_year_temp_adjusted = np.where(cmip_year_temp_anomaly<0, 
                                                    cmip_year_temp.mean()+cmip_year_temp_anomaly*future_cooling_factor,
                                                    cmip_year_temp_adjusted)

                adjusted_temp.loc[adjusted_temp.index.year==year] = cmip_year_temp_adjusted

                cmip_year_prcp = adjusted_prcp.loc[adjusted_prcp.index.year==year]
                cmip_year_prcp_anomaly = cmip_year_prcp - cmip_year_prcp.mean()
                cmip_year_prcp_adjusted = np.where(cmip_year_prcp_anomaly<=0,
                                                    cmip_year_prcp.mean()+cmip_year_prcp_anomaly*year_info['prcp_scaling'],
                                                    cmip_year_prcp)
                
                adjusted_prcp.loc[adjusted_prcp.index.year==year] = cmip_year_prcp_adjusted

            adjusted_prcp = adjusted_prcp.mask(adjusted_prcp<0, 0)
            adjusted_df = pd.DataFrame({
                'prcp': adjusted_prcp,
                'temp': adjusted_temp
            })
            
            # turn to nc
            adjusted_ds = target_cmip.copy(deep=False)
            
            # 更新 prcp 和 temp 的值
            adjusted_ds['prcp'] = xr.DataArray(
                adjusted_df['prcp'].values,
                dims=target_cmip.prcp.dims,
                coords={dim: target_cmip.coords[dim] for dim in target_cmip.prcp.dims},
                attrs=target_cmip.prcp.attrs if hasattr(target_cmip.prcp, 'attrs') else {}
            )
            adjusted_ds['temp'] = xr.DataArray(
                adjusted_df['temp'].values,
                dims=target_cmip.temp.dims,
                coords={dim: target_cmip.coords[dim] for dim in target_cmip.temp.dims},
                attrs=target_cmip.temp.attrs if hasattr(target_cmip.temp, 'attrs') else {}
            )
            adjusted_ds['temp'].loc[{'time': hist_date_range}] = ref_era5['temp'].values # for hist periods, using ERA5 data
            adjusted_ds['prcp'].loc[{'time': hist_date_range}] = ref_era5['prcp'].values
            adjusted_ds.attrs['ref_hgt'] = ref_era5.attrs['ref_hgt']
            adjusted_ds.attrs['ref_pix_lon'] = ref_era5.attrs['ref_pix_lon']
            adjusted_ds.attrs['ref_pix_lat'] = ref_era5.attrs['ref_pix_lat']
            adjusted_ds.attrs['ref_pix_dis'] = ref_era5.attrs['ref_pix_dis']
            adjusted_ds.attrs['valid_future_years'] = len(adjusted_ds['time'])/12 - 75 # 75 years from 2025 to 2099

            if np.any(np.isnan(adjusted_ds.temp.values)) or np.any(np.isnan(adjusted_ds.prcp.values)):
                raise ValueError(f"NaN values found in {GCM}_{SSP} for {gdir.rgi_id}")
            # dims_to_expand = []
            # for coord_name in ['rgi_id', 'GCM', 'SSP']:
            #     if coord_name in adjusted_ds.coords and coord_name not in adjusted_ds.dims:
            #         dims_to_expand.append(coord_name)
            
            # if dims_to_expand:
            #     adjusted_ds = adjusted_ds.expand_dims(dims_to_expand)
            out_path = os.path.join(gdir.dir, 'gcm_data_{}_{}{}.nc'.format(GCM, SSP, output_filesuffix))
            adjusted_ds.to_netcdf(out_path)
            adjusted_ds.close()

    return True

@entity_task(log)
def simulate_future_extremes_QDM(gdir, 
                                start_year=1950,
                                end_year=2100,
                                future_start_year=2025,
                                future_cooling_factor=None, #>1, for sensitive test
                                output_filesuffix=''):
    """
    Utilize Quantile Delta Mapping (QDM) to adjust future extremes.
    x ̂_(m,p) (t)=F_(o,h)^(-1) (τ_t )+(x_(m,p) (t)-F_(m,h)^(-1) (τ_t))
    That is, the adjusted value is the quantile of the historical distribution (ERA5-Land)
    at the same quantile of the future distribution, with warming added provided by CMIP6 data.
    Must groupby months!
    Parameters
    ----------
    gdir : GlacierDirectory
        Glacier directory object
    start_year : int
        Starting year for simulation
    end_year : int
        Ending year for simulation
    future_start_year : int
        Start year of future period
    output_filesuffix : str
        Suffix for output filename
    
    Returns
    -------
    bool
        True if successful
    """
    from statsmodels.tsa.seasonal import STL
    # CMIP6 scenarios to analyze
    cmip6_scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp585']
    cmip6_models = ['BCC-CSM2-MR', 'CAMS-CSM1-0', 'CESM2', 'CESM2-WACCM', 
                    'EC-Earth3', 'EC-Earth3-Veg', 'FGOALS-f3-L', 'GFDL-ESM4',
                    'INM-CM4-8', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NorESM2-MM']
    
    if os.path.exists(os.path.join(gdir.dir, 'climate_historical.nc')):
        ref_era5_path = os.path.join(gdir.dir, 'climate_historical.nc')
        ref_era5 = xr.open_dataset(ref_era5_path).sel(time=slice(f"{start_year}-01-01", f"{future_start_year-1}-12-01"))
    else:
        raise FileNotFoundError(f"ERA5 climate file not found for {gdir.rgi_id}")

    hist_date_range = pd.date_range(start=f"01/01/{start_year}", end=f"12/01/{future_start_year-1}", freq="MS")
    ref_era5['time'] = hist_date_range

    ref_date_range = pd.date_range(start=f"01/01/1990", end=f"12/01/{future_start_year-1}", freq="MS")
    ref_era5_df = ref_era5.to_dataframe().loc[ref_date_range] # 1990-2024

    for GCM in cmip6_models:  # loop through all GCMs
        for SSP in cmip6_scenarios:  # loop through all SSPs
            cmip_path = os.path.join(gdir.dir, 'gcm_data_{}_{}.nc'.format(GCM, SSP))
            target_cmip = xr.open_dataset(cmip_path).sel(time=slice(f"{start_year}-01-01", f"{end_year}-12-31"))

            target_cmip['time'] = pd.date_range(start=f"01/01/{start_year}", end=target_cmip.time.max().values, freq="MS")
            target_cmip_df = target_cmip.to_dataframe() # structure same to selected_era5，2001-01-01
            target_cmip_df['year'] = target_cmip_df.index.year
            target_cmip_df['month'] = target_cmip_df.index.month

            future_date_range = pd.date_range(start=f"01/01/{future_start_year}", end=target_cmip.time.max().values, freq="MS")
            future_cmip_df = target_cmip_df.loc[future_date_range]
            hist_cmip_df = target_cmip_df.loc[ref_date_range]

            adjusted_temp = target_cmip_df['temp'].copy()
            adjusted_prcp = target_cmip_df['prcp'].copy()
            adjusted_temp.loc[hist_date_range] = ref_era5.to_dataframe().loc[hist_date_range, 'temp']
            adjusted_prcp.loc[hist_date_range] = ref_era5.to_dataframe().loc[hist_date_range, 'prcp']

            # future adjustment using QDM
            for m in range(1, 13): # Loop 1 to 12
                m_hist_era5_df = ref_era5_df[ref_era5_df.index.month==m]
                m_hist_cmip_df = hist_cmip_df[hist_cmip_df.index.month==m]
                m_future_cmip_df = future_cmip_df[future_cmip_df.index.month==m].copy()

                m_future_cmip_df['τ_temp'] = m_future_cmip_df['temp'].rank(method='average', pct=True)
                m_future_cmip_df['τ_prcp'] = m_future_cmip_df['prcp'].rank(method='average', pct=True)

                for date, date_info in m_future_cmip_df.iterrows():
                    cmip_temp = date_info['temp']
                    τ_temp = date_info['τ_temp'] # quantile of future temp
                    
                    hist_era5_τ_temp = np.nanpercentile(m_hist_era5_df['temp'], τ_temp*100)
                    hist_cmip6_τ_temp = np.nanpercentile(m_hist_cmip_df['temp'], τ_temp*100)
                    corrected_cmip_temp = hist_era5_τ_temp + (cmip_temp - hist_cmip6_τ_temp)
                    # for prcp, times
                    cmip_prcp = date_info['prcp']
                    τ_prcp = date_info['τ_prcp'] # quantile of future prcp

                    hist_era5_τ_prcp = np.nanpercentile(m_hist_era5_df['prcp'], τ_prcp*100)
                    hist_cmip6_τ_prcp = np.nanpercentile(m_hist_cmip_df['prcp'], τ_prcp*100)
                    scalling_factor = cmip_prcp / hist_cmip6_τ_prcp
                    if (scalling_factor > 4) or (np.isinf(scalling_factor)) or (np.isnan(scalling_factor)):
                        scalling_factor = np.nanpercentile(m_future_cmip_df['prcp'], 50) / np.nanpercentile(m_hist_cmip_df['prcp'], 50) # use general trend
                        if (scalling_factor > 4) or (np.isinf(scalling_factor)) or (np.isnan(scalling_factor)):
                            scalling_factor = np.nanpercentile(future_cmip_df['prcp'], 50) / np.nanpercentile(hist_cmip_df['prcp'], 50)

                    corrected_cmip_prcp = hist_era5_τ_prcp * scalling_factor

                    adjusted_temp.loc[date] = corrected_cmip_temp
                    adjusted_prcp.loc[date] = corrected_cmip_prcp
        
            adjusted_prcp = adjusted_prcp.mask(adjusted_prcp<0, 0)
            adjusted_df = pd.DataFrame({
                'prcp': adjusted_prcp,
                'temp': adjusted_temp
            })

            if future_cooling_factor is not None:
                future_temp_series = adjusted_temp.loc[future_date_range] # future series
                future_temp_annual = future_temp_series.groupby(future_temp_series.index.year)
                for year, year_info in future_temp_annual:
                    year_temp_anomaly = year_info - year_info.mean()
                    year_temp_adjusted = np.where(year_temp_anomaly<0, 
                                                year_info.mean()+year_temp_anomaly*future_cooling_factor,
                                                year_info)
                    
                    adjusted_df.loc[year_info.index, 'temp'] = year_temp_adjusted
            
            # turn to nc
            adjusted_ds = target_cmip.copy(deep=False)
            
            # 更新 prcp 和 temp 的值
            adjusted_ds['prcp'] = xr.DataArray(
                adjusted_df['prcp'].values,
                dims=target_cmip.prcp.dims,
                coords={dim: target_cmip.coords[dim] for dim in target_cmip.prcp.dims},
                attrs=target_cmip.prcp.attrs if hasattr(target_cmip.prcp, 'attrs') else {}
            )
            adjusted_ds['temp'] = xr.DataArray(
                adjusted_df['temp'].values,
                dims=target_cmip.temp.dims,
                coords={dim: target_cmip.coords[dim] for dim in target_cmip.temp.dims},
                attrs=target_cmip.temp.attrs if hasattr(target_cmip.temp, 'attrs') else {}
            )
            adjusted_ds.attrs['ref_hgt'] = ref_era5.attrs['ref_hgt']
            adjusted_ds.attrs['ref_pix_lon'] = ref_era5.attrs['ref_pix_lon']
            adjusted_ds.attrs['ref_pix_lat'] = ref_era5.attrs['ref_pix_lat']
            adjusted_ds.attrs['ref_pix_dis'] = ref_era5.attrs['ref_pix_dis']
            adjusted_ds.attrs['valid_future_years'] = len(adjusted_ds['time'])/12 - 75 # 75 years from 2025 to 2099

            if np.any(np.isnan(adjusted_ds.temp.values)) or np.any(np.isnan(adjusted_ds.prcp.values)):
                raise ValueError(f"NaN values found in {GCM}_{SSP} for {gdir.rgi_id}")

            out_path = os.path.join(gdir.dir, 'gcm_data_{}_{}{}.nc'.format(GCM, SSP, output_filesuffix))
            adjusted_ds.to_netcdf(out_path)
            adjusted_ds.close()

    return True

@entity_task(log)
def simulate_future_extremes_Detrend_QM(gdir, 
                                        start_year=1950,
                                        end_year=2100,
                                        future_start_year=2025,
                                        future_cooling_factor=None, #>1, for sensitive test
                                        output_filesuffix=''):
    """
    Utilize Quantile Mapping (QM) to adjust future extremes.
    x_(m,p) (t)=〖trend〗_(m,p) (t)+〖flu〗_(m,p) (t)
    τ_t=F_(m,p) (〖flu〗_(m,p) (t))
    x ̂_(m,p) (t)=〖trend〗_(m,p) (t)+F_(o,h)^(-1) (τ_t )
    Must groupby months!
    Parameters
    ----------
    gdir : GlacierDirectory
        Glacier directory object
    start_year : int
        Starting year for simulation
    end_year : int
        Ending year for simulation
    future_start_year : int
        Start year of future period
    output_filesuffix : str
        Suffix for output filename
    
    Returns
    -------
    bool
        True if successful
    """
    from statsmodels.tsa.seasonal import STL
    # CMIP6 scenarios to analyze
    cmip6_scenarios = ['ssp126', 'ssp245', 'ssp370', 'ssp585']
    cmip6_models = ['BCC-CSM2-MR', 'CAMS-CSM1-0', 'CESM2', 'CESM2-WACCM', 
                    'EC-Earth3', 'EC-Earth3-Veg', 'FGOALS-f3-L', 'GFDL-ESM4',
                    'INM-CM4-8', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'NorESM2-MM']
    
    if os.path.exists(os.path.join(gdir.dir, 'climate_historical.nc')):
        era5_path = os.path.join(gdir.dir, 'climate_historical.nc')
        era5_ds = xr.open_dataset(era5_path).sel(time=slice(f"{start_year}-01-01", f"{future_start_year-1}-12-01"))
    else:
        raise FileNotFoundError(f"ERA5 climate file not found for {gdir.rgi_id}")

    hist_date_range = pd.date_range(start=f"01/01/{start_year}", end=f"12/01/{future_start_year-1}", freq="MS") # 1950-2024
    era5_ds['time'] = hist_date_range
    era5_df = era5_ds.to_dataframe() # df
    era5_temp_trend = STL(era5_df['temp'], period=12, trend=241, seasonal=13).fit().trend
    era5_prcp_trend = STL(era5_df['prcp'], period=12, trend=241, seasonal=13).fit().trend
    era5_df['temp_trend'] = era5_temp_trend
    era5_df['temp_leftover'] = era5_df['temp'] - era5_df['temp_trend']
    era5_df['prcp_trend'] = era5_prcp_trend
    era5_df['prcp_leftover'] = era5_df['prcp'] - era5_df['prcp_trend'] # fluc

    ref_date_range = pd.date_range(start=f"01/01/1990", end=f"12/01/{future_start_year-1}", freq="MS")
    ref_era5_df = era5_df.loc[ref_date_range] # 1990-2024

    for GCM in cmip6_models:  # loop through all GCMs
        for SSP in cmip6_scenarios:  # loop through all SSPs
            cmip_path = os.path.join(gdir.dir, 'gcm_data_{}_{}.nc'.format(GCM, SSP))
            target_cmip = xr.open_dataset(cmip_path).sel(time=slice(f"{start_year}-01-01", f"{end_year}-12-31"))

            target_cmip['time'] = pd.date_range(start=f"01/01/{start_year}", end=target_cmip.time.max().values, freq="MS")
            target_cmip_df = target_cmip.to_dataframe() # structure same to selected_era5，2001-01-01
            target_cmip_df['year'] = target_cmip_df.index.year
            target_cmip_df['month'] = target_cmip_df.index.month
            # STL trend
            cmip_temp_trend = STL(target_cmip_df['temp'], period=12, trend=241, seasonal=13).fit().trend
            cmip_prcp_trend = STL(target_cmip_df['prcp'], period=12, trend=241, seasonal=13).fit().trend
            target_cmip_df['temp_trend'] = cmip_temp_trend
            target_cmip_df['temp_leftover'] = target_cmip_df['temp'] - target_cmip_df['temp_trend']
            target_cmip_df['prcp_trend'] = cmip_prcp_trend
            target_cmip_df['prcp_leftover'] = target_cmip_df['prcp'] - target_cmip_df['prcp_trend'] # fluc

            future_date_range = pd.date_range(start=f"01/01/{future_start_year}", end=target_cmip.time.max().values, freq="MS")
            future_cmip_df = target_cmip_df.loc[future_date_range]
            hist_cmip_df = target_cmip_df.loc[ref_date_range]

            adjusted_temp = target_cmip_df['temp'].copy()
            adjusted_prcp = target_cmip_df['prcp'].copy()
            adjusted_temp.loc[hist_date_range] = era5_df.loc[hist_date_range, 'temp']
            adjusted_prcp.loc[hist_date_range] = era5_df.loc[hist_date_range, 'prcp']

            for m in range(1, 13): # Loop 1 to 12, by months
                # specifically for month m
                m_hist_era5_df = ref_era5_df[ref_era5_df.index.month==m]
                m_hist_cmip_df = hist_cmip_df[hist_cmip_df.index.month==m]
                m_future_cmip_df = future_cmip_df[future_cmip_df.index.month==m].copy()

                m_future_cmip_df['τ_temp'] = m_future_cmip_df['temp_leftover'].rank(method='average', pct=True)
                m_future_cmip_df['τ_prcp'] = m_future_cmip_df['prcp_leftover'].rank(method='average', pct=True)

                for date, date_info in m_future_cmip_df.iterrows():
                    cmip_temp_trend = date_info['temp_trend']
                    τ_temp = date_info['τ_temp']
                    hist_era5_τ_temp_flu = np.nanpercentile(m_hist_era5_df['temp_leftover'], τ_temp*100)
                    corrected_cmip_temp = cmip_temp_trend + hist_era5_τ_temp_flu
                    # for prcp
                    cmip_prcp_trend = date_info['prcp_trend']
                    τ_prcp = date_info['τ_prcp']
                    hist_era5_τ_prcp_flu = np.nanpercentile(m_hist_era5_df['prcp_leftover'], τ_prcp*100)
                    corrected_cmip_prcp = cmip_prcp_trend + hist_era5_τ_prcp_flu
                    if corrected_cmip_prcp<0:
                        corrected_cmip_prcp = 0
                    
                    adjusted_temp.loc[date] = corrected_cmip_temp
                    adjusted_prcp.loc[date] = corrected_cmip_prcp
            
            adjusted_prcp = adjusted_prcp.mask(adjusted_prcp<0, 0)
            adjusted_df = pd.DataFrame({
                'prcp': adjusted_prcp,
                'temp': adjusted_temp
            })

            if future_cooling_factor is not None:
                future_temp_series = adjusted_temp.loc[future_date_range] # future series
                future_temp_annual = future_temp_series.groupby(future_temp_series.index.year)
                for year, year_info in future_temp_annual:
                    year_temp_anomaly = year_info - year_info.mean()
                    year_temp_adjusted = np.where(year_temp_anomaly<0, 
                                                year_info.mean()+year_temp_anomaly*future_cooling_factor,
                                                year_info)
                    
                    adjusted_df.loc[year_info.index, 'temp'] = year_temp_adjusted
            
            # turn to nc
            adjusted_ds = target_cmip.copy(deep=False)
            
            # 更新 prcp 和 temp 的值
            adjusted_ds['prcp'] = xr.DataArray(
                adjusted_df['prcp'].values,
                dims=target_cmip.prcp.dims,
                coords={dim: target_cmip.coords[dim] for dim in target_cmip.prcp.dims},
                attrs=target_cmip.prcp.attrs if hasattr(target_cmip.prcp, 'attrs') else {}
            )
            adjusted_ds['temp'] = xr.DataArray(
                adjusted_df['temp'].values,
                dims=target_cmip.temp.dims,
                coords={dim: target_cmip.coords[dim] for dim in target_cmip.temp.dims},
                attrs=target_cmip.temp.attrs if hasattr(target_cmip.temp, 'attrs') else {}
            )
            adjusted_ds.attrs['ref_hgt'] = era5_ds.attrs['ref_hgt']
            adjusted_ds.attrs['ref_pix_lon'] = era5_ds.attrs['ref_pix_lon']
            adjusted_ds.attrs['ref_pix_lat'] = era5_ds.attrs['ref_pix_lat']
            adjusted_ds.attrs['ref_pix_dis'] = era5_ds.attrs['ref_pix_dis']
            adjusted_ds.attrs['valid_future_years'] = len(adjusted_ds['time'])/12 - 75 # 75 years from 2025 to 2099

            if np.any(np.isnan(adjusted_ds.temp.values)) or np.any(np.isnan(adjusted_ds.prcp.values)):
                raise ValueError(f"NaN values found in {GCM}_{SSP} for {gdir.rgi_id}")

            out_path = os.path.join(gdir.dir, 'gcm_data_{}_{}{}.nc'.format(GCM, SSP, output_filesuffix))
            adjusted_ds.to_netcdf(out_path)
            adjusted_ds.close()
    return True