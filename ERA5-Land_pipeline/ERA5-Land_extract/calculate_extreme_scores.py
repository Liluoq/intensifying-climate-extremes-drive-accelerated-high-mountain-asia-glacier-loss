"""
Calculate extreme climate scores for each glacier from monthly aggregated data.

This script analyzes temperature (hot extremes) and precipitation (dry extremes)
to compute annual extreme scores for each glacier.

Methodology:
- Temperature: 99th, 95th, 90th percentiles (hot extremes)
- Precipitation: 10th, 5th, 1st percentiles (dry extremes)
- Annual compound scores combining both metrics
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
import logging
from typing import Tuple, Dict

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def load_config(config_path: str = "config.json") -> dict:
    """Load configuration from JSON file."""
    config_file = Path(__file__).parent.parent / config_path
    with open(config_file, 'r') as f:
        config = json.load(f)
    logger.info(f"Loaded configuration from {config_file}")
    return config


def load_monthly_data(
    xlsx_path: str, 
    start_date: str = None, 
    end_date: str = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load monthly data from Excel file with optional time filtering.
    
    Parameters:
    -----------
    xlsx_path : str
        Path to Excel file with monthly data
    start_date : str, optional
        Start date in format 'YYYY-MM' (inclusive)
    end_date : str, optional
        End date in format 'YYYY-MM' (inclusive)
        
    Returns:
    --------
    tuple
        (temperature_df, precipitation_df)
        Each DataFrame has glacier_ids as index, months as columns
    """
    logger.info(f"Loading monthly data from {xlsx_path}")
    
    # Read temperature data
    temp_df = pd.read_excel(xlsx_path, sheet_name='temperature_2m', index_col=0)
    logger.info(f"  Temperature data: {temp_df.shape} (glaciers × months)")
    
    # Read precipitation data
    precip_df = pd.read_excel(xlsx_path, sheet_name='total_precipitation_sum', index_col=0)
    logger.info(f"  Precipitation data: {precip_df.shape} (glaciers × months)")
    
    # Filter by time period if specified
    if start_date or end_date:
        logger.info(f"Filtering time period: {start_date or 'start'} to {end_date or 'end'}")
        
        # Convert column names to Period objects for robust comparison
        # Columns should be in format "YYYY-MM" or similar
        temp_periods = pd.to_datetime(temp_df.columns, format='%Y-%m', errors='coerce').to_period('M')
        precip_periods = pd.to_datetime(precip_df.columns, format='%Y-%m', errors='coerce').to_period('M')
        
        # Convert start_date and end_date to Period objects
        start_period = pd.Period(start_date, freq='M') if start_date else None
        end_period = pd.Period(end_date, freq='M') if end_date else None
        
        # Create mask for columns to keep
        temp_mask = pd.Series(True, index=range(len(temp_periods)))
        precip_mask = pd.Series(True, index=range(len(precip_periods)))
        
        if start_period:
            temp_mask = temp_mask & (temp_periods >= start_period)
            precip_mask = precip_mask & (precip_periods >= start_period)
        
        if end_period:
            temp_mask = temp_mask & (temp_periods < end_period)
            precip_mask = precip_mask & (precip_periods < end_period)
        
        # Apply filter
        temp_df = temp_df.iloc[:, temp_mask.values]
        precip_df = precip_df.iloc[:, precip_mask.values]
        
        logger.info(f"  After filtering - Temperature: {temp_df.shape}, Precipitation: {precip_df.shape}")
        if temp_df.shape[1] > 0:
            logger.info(f"  Date range: {temp_df.columns[0]} to {temp_df.columns[-1]}")
        else:
            logger.warning("  No data remaining after filtering!")
    
    return temp_df, precip_df


def calculate_anomalies(data: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate anomalies by subtracting mean across all months for each glacier.
    
    Parameters:
    -----------
    data : pd.DataFrame
        Monthly data (glaciers × months)
        
    Returns:
    --------
    pd.DataFrame
        Anomaly values (glaciers × months)
    """
    logger.info("Calculating anomalies...")
    
    # For each glacier (row), subtract its mean across all months
    anomalies = data.sub(data.mean(axis=1), axis=0)
    
    logger.info(f"  Anomalies shape: {anomalies.shape}")
    logger.info(f"  Mean of anomalies: {anomalies.values.mean():.6f} (should be ~0)")
    
    return anomalies


def calculate_quantiles(data: pd.DataFrame, percentiles: list) -> pd.DataFrame:
    """
    Calculate quantiles for each glacier.
    
    Parameters:
    -----------
    data : pd.DataFrame
        Anomaly data (glaciers × months)
    percentiles : list
        List of percentiles to calculate (e.g., [90, 95, 99])
        
    Returns:
    --------
    pd.DataFrame
        Quantile values (glaciers × percentiles)
    """
    logger.info(f"Calculating quantiles: {percentiles}")
    
    # Calculate quantiles for each glacier (row-wise)
    quantiles = data.quantile(q=[p/100 for p in percentiles], axis=1).T
    quantiles.columns = [f'p{p}' for p in percentiles]
    
    logger.info(f"  Quantiles shape: {quantiles.shape}")
    
    return quantiles


def assign_scores_temperature(anomalies: pd.DataFrame, quantiles: pd.DataFrame) -> pd.DataFrame:
    """
    Assign temperature extreme scores based on quantile thresholds.
    
    Hot extremes: 99th=4, 95th=2, 90th=1, else=0
    
    Parameters:
    -----------
    anomalies : pd.DataFrame
        Temperature anomalies (glaciers × months)
    quantiles : pd.DataFrame
        Quantile thresholds (glaciers × [p90, p95, p99])
        
    Returns:
    --------
    pd.DataFrame
        Monthly scores (glaciers × months)
    """
    logger.info("Assigning temperature extreme scores...")
    
    scores = pd.DataFrame(0.0, index=anomalies.index, columns=anomalies.columns)
    
    # For each glacier
    for glacier_id in anomalies.index:
        glacier_anomalies = anomalies.loc[glacier_id]
        p90 = quantiles.loc[glacier_id, 'p90']
        p95 = quantiles.loc[glacier_id, 'p95']
        p99 = quantiles.loc[glacier_id, 'p99']
        
        # Assign scores (check from tightest to loosest)
        for month in anomalies.columns:
            value = glacier_anomalies[month]
            
            if pd.isna(value):
                scores.loc[glacier_id, month] = 0
            elif value >= p99:
                scores.loc[glacier_id, month] = 4
            elif value >= p95:
                scores.loc[glacier_id, month] = 2
            elif value >= p90:
                scores.loc[glacier_id, month] = 1
            else:
                scores.loc[glacier_id, month] = 0
    
    total_score = scores.sum().sum()
    logger.info(f"  Total temperature scores: {total_score:.1f}")
    
    return scores


def assign_scores_precipitation(anomalies: pd.DataFrame, quantiles: pd.DataFrame) -> pd.DataFrame:
    """
    Assign precipitation extreme scores based on quantile thresholds.
    
    Dry extremes: 1st=2, 5th=1, 10th=0.5, else=0
    
    Parameters:
    -----------
    anomalies : pd.DataFrame
        Precipitation anomalies (glaciers × months)
    quantiles : pd.DataFrame
        Quantile thresholds (glaciers × [p1, p5, p10])
        
    Returns:
    --------
    pd.DataFrame
        Monthly scores (glaciers × months)
    """
    logger.info("Assigning precipitation extreme scores...")
    
    scores = pd.DataFrame(0.0, index=anomalies.index, columns=anomalies.columns)
    
    # For each glacier
    for glacier_id in anomalies.index:
        glacier_anomalies = anomalies.loc[glacier_id]
        p1 = quantiles.loc[glacier_id, 'p1']
        p5 = quantiles.loc[glacier_id, 'p5']
        p10 = quantiles.loc[glacier_id, 'p10']
        
        # Assign scores (check from tightest to loosest)
        for month in anomalies.columns:
            value = glacier_anomalies[month]
            
            if pd.isna(value):
                scores.loc[glacier_id, month] = 0
            elif value <= p1:
                scores.loc[glacier_id, month] = 2
            elif value <= p5:
                scores.loc[glacier_id, month] = 1
            elif value <= p10:
                scores.loc[glacier_id, month] = 0.5
            else:
                scores.loc[glacier_id, month] = 0
    
    total_score = scores.sum().sum()
    logger.info(f"  Total precipitation scores: {total_score:.1f}")
    
    return scores


def aggregate_to_annual(monthly_scores: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate monthly scores to annual scores by summing within each year.
    
    Parameters:
    -----------
    monthly_scores : pd.DataFrame
        Monthly scores (glaciers × months with YYYY-MM column names)
        
    Returns:
    --------
    pd.DataFrame
        Annual scores (glaciers × years)
    """
    logger.info("Aggregating monthly scores to annual...")
    
    # Ensure column names are strings and extract year
    # Column names should be in format "YYYY-MM" (e.g., "1990-01")
    years = []
    for col in monthly_scores.columns:
        col_str = str(col)
        # Extract year (first 4 characters if format is YYYY-MM)
        if '-' in col_str:
            year = col_str.split('-')[0]
        else:
            # If no hyphen, try to get first 4 characters
            year = col_str[:4]
        years.append(year)
    
    # Group columns by year and sum
    annual_scores = monthly_scores.T.groupby(years).sum().T
    
    logger.info(f"  Annual scores shape: {annual_scores.shape}")
    logger.info(f"  Years: {list(annual_scores.columns)[:5]}...{list(annual_scores.columns)[-5:]}")
    
    return annual_scores


def calculate_compound_scores(temp_annual: pd.DataFrame, precip_annual: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate compound scores by adding temperature and precipitation annual scores.
    
    Parameters:
    -----------
    temp_annual : pd.DataFrame
        Annual temperature scores (glaciers × years)
    precip_annual : pd.DataFrame
        Annual precipitation scores (glaciers × years)
        
    Returns:
    --------
    pd.DataFrame
        Annual compound scores (glaciers × years)
    """
    logger.info("Calculating compound scores...")
    
    # Add the two DataFrames
    compound = temp_annual + precip_annual
    
    logger.info(f"  Compound scores shape: {compound.shape}")
    logger.info(f"  Mean compound score: {compound.values.mean():.2f}")
    logger.info(f"  Max compound score: {compound.values.max():.2f}")
    
    return compound


def write_scores_to_excel(
    temp_scores: pd.DataFrame,
    precip_scores: pd.DataFrame,
    compound_scores: pd.DataFrame,
    output_path: str
):
    """
    Write scores to Excel file with 3 sheets.
    
    Parameters:
    -----------
    temp_scores : pd.DataFrame
        Annual temperature scores (glaciers × years)
    precip_scores : pd.DataFrame
        Annual precipitation scores (glaciers × years)
    compound_scores : pd.DataFrame
        Annual compound scores (glaciers × years)
    output_path : str
        Path to output Excel file
    """
    logger.info(f"Writing scores to Excel: {output_path}")
    
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        # Sheet 1: Temperature scores
        temp_scores.to_excel(writer, sheet_name='Temperature_Extreme_Scores')
        logger.info("  Written: Temperature_Extreme_Scores")
        
        # Sheet 2: Precipitation scores
        precip_scores.to_excel(writer, sheet_name='Precipitation_Extreme_Scores')
        logger.info("  Written: Precipitation_Extreme_Scores")
        
        # Sheet 3: Compound scores
        compound_scores.to_excel(writer, sheet_name='Compound_Extreme_Scores')
        logger.info("  Written: Compound_Extreme_Scores")
    
    logger.info("Excel file written successfully")


def main(config_path: str = "config.json"):
    """
    Main function to calculate and save extreme scores.
    
    Parameters:
    -----------
    config_path : str
        Path to configuration file
    """
    logger.info("="*80)
    logger.info("STARTING EXTREME SCORE CALCULATION")
    logger.info("="*80)
    
    # Load configuration
    config = load_config(config_path)
    
    monthly_xlsx = config.get('monthly_output_xlsx')
    output_xlsx = config.get('extreme_scores_output',
                            '/WORK/Data/ncc_glacier/era5_land/extracted/glacier_extreme_scores.xlsx')
    start_date = config.get('extreme_scores_start_date', None)
    end_date = config.get('extreme_scores_end_date', None)
    
    logger.info(f"Input: {monthly_xlsx}")
    logger.info(f"Output: {output_xlsx}")
    if start_date or end_date:
        logger.info(f"Time period filter: {start_date or 'start'} to {end_date or 'end'}")
    
    # Create output directory if needed
    output_dir = Path(output_xlsx).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # STEP 1: Load monthly data with time filtering
    temp_df, precip_df = load_monthly_data(monthly_xlsx, start_date, end_date)
    
    # STEP 2: Calculate anomalies
    logger.info("\n" + "-"*80)
    logger.info("TEMPERATURE PROCESSING")
    logger.info("-"*80)
    temp_anomalies = calculate_anomalies(temp_df)
    
    logger.info("\n" + "-"*80)
    logger.info("PRECIPITATION PROCESSING")
    logger.info("-"*80)
    precip_anomalies = calculate_anomalies(precip_df)
    
    # STEP 3: Calculate quantiles
    temp_quantiles = calculate_quantiles(temp_anomalies, [90, 95, 99])
    precip_quantiles = calculate_quantiles(precip_anomalies, [1, 5, 10])
    
    # STEP 4: Assign monthly scores
    temp_monthly_scores = assign_scores_temperature(temp_anomalies, temp_quantiles)
    precip_monthly_scores = assign_scores_precipitation(precip_anomalies, precip_quantiles)
    
    # STEP 5: Aggregate to annual
    logger.info("\n" + "-"*80)
    logger.info("ANNUAL AGGREGATION")
    logger.info("-"*80)
    temp_annual = aggregate_to_annual(temp_monthly_scores)
    precip_annual = aggregate_to_annual(precip_monthly_scores)
    
    # STEP 6: Calculate compound scores
    compound_annual = calculate_compound_scores(temp_annual, precip_annual)
    
    # STEP 7: Write to Excel
    logger.info("\n" + "-"*80)
    logger.info("WRITING OUTPUT")
    logger.info("-"*80)
    write_scores_to_excel(temp_annual, precip_annual, compound_annual, output_xlsx)
    
    # Summary statistics
    logger.info("\n" + "="*80)
    logger.info("EXTREME SCORE CALCULATION COMPLETE")
    logger.info("="*80)
    logger.info(f"Glaciers processed: {len(temp_df)}")
    logger.info(f"Months processed: {len(temp_df.columns)}")
    logger.info(f"Years in output: {len(temp_annual.columns)}")
    logger.info(f"\nTemperature scores:")
    logger.info(f"  Mean annual: {temp_annual.values.mean():.2f}")
    logger.info(f"  Max annual: {temp_annual.values.max():.2f}")
    logger.info(f"\nPrecipitation scores:")
    logger.info(f"  Mean annual: {precip_annual.values.mean():.2f}")
    logger.info(f"  Max annual: {precip_annual.values.max():.2f}")
    logger.info(f"\nCompound scores:")
    logger.info(f"  Mean annual: {compound_annual.values.mean():.2f}")
    logger.info(f"  Max annual: {compound_annual.values.max():.2f}")
    logger.info(f"\nOutput saved to: {output_xlsx}")


if __name__ == "__main__":
    main()

