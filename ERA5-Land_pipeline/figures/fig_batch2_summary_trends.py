"""
Figure Batch 3: Summary Trend Plots (ED_fig3d-e)
================================================

Creates two figures for Extended Data Figure 3:

ED_fig3d: Time Series of HMA-wide Mean Compound Score
    - Shows annual mean compound score across all glaciers
    - Includes linear trend line with slope and p-value

ED_fig3e: Bar Chart of Top-5 Extreme Year Distribution by Period
    - For each glacier, identifies its top 5 extreme years
    - Shows what % of those top-5 years fall in each 5-year period
    - Increasing bars toward recent periods = increasing frequency

These complement the spatial maps from Batch 1 (ED_fig3a-c).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Font settings for Adobe Illustrator compatibility
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


# =============================================================================
# CONFIGURATION
# =============================================================================

class Config:
    """Configuration parameters for Figure Batch 3"""

    # --- Data period ---
    START_YEAR = 1990
    END_YEAR = 2024

    # --- Data paths ---
    EXTREME_SCORES_XLSX = "/WORK/Data/ncc_glacier/era5_land/extracted/glacier_extreme_scores.xlsx"

    # --- Output ---
    OUTPUT_DIR = Path(__file__).parent / "output"
    DPI = 600

    # --- Score types to process (compound only for ED_fig3d-e) ---
    SCORE_TYPES = {
        'compound': 'Compound_Extreme_Scores'
    }

    # --- Top N extreme years (for frequency analysis) ---
    TOP_N_EXTREME_YEARS = 5

    # --- Period definitions for bar chart ---
    PERIODS = [
        (1990, 1994),
        (1995, 1999),
        (2000, 2004),
        (2005, 2009),
        (2010, 2014),
        (2015, 2019),
        (2020, 2024),
    ]

    # --- Styling ---
    FIGURE_SIZE = (8, 5)

    # Colors
    LINE_COLOR = '#2166AC'  # Blue for time series
    TREND_COLOR = '#B2182B'  # Red for trend line
    BAR_COLOR_EARLY = '#67A9CF'  # Light blue for early period bars
    BAR_COLOR_RECENT = '#B2182B'  # Red for recent period bars
    RECENT_PERIOD_START = 2015  # Bars from this year onward are colored red


# =============================================================================
# DATA LOADING
# =============================================================================

def load_extreme_scores(score_type):
    """Load extreme scores from Excel file for a given score type."""
    sheet_name = Config.SCORE_TYPES[score_type]
    print(f"Loading {score_type} extreme scores...")
    df = pd.read_excel(
        Config.EXTREME_SCORES_XLSX,
        sheet_name=sheet_name
    )
    print(f"  Shape: {df.shape}")

    # Extract glacier_id and year columns
    glacier_id_col = df.columns[0]
    year_cols = [c for c in df.columns if isinstance(c, int) or (isinstance(c, str) and c.isdigit())]
    year_cols = [int(c) for c in year_cols]

    # Filter to configured period
    year_cols = [y for y in year_cols if Config.START_YEAR <= y <= Config.END_YEAR]
    year_cols = sorted(year_cols)

    print(f"  Years: {year_cols[0]} to {year_cols[-1]} ({len(year_cols)} years)")
    print(f"  Glaciers: {len(df)}")

    return df, glacier_id_col, year_cols


# =============================================================================
# FIGURE 3A: TIME SERIES WITH TREND (INTENSITY)
# =============================================================================

def create_timeseries_figure(df, glacier_id_col, year_cols, score_type):
    """
    Create time series of HMA-wide mean score with trend line.
    """
    score_label = score_type.capitalize()
    print(f"\nCreating {score_label} Intensity Time Series...")

    # Calculate annual mean across all glaciers
    annual_means = []
    for year in year_cols:
        year_col = year if year in df.columns else str(year)
        if year_col in df.columns:
            mean_val = df[year_col].mean()
            annual_means.append({'year': year, 'mean_score': mean_val})

    ts_df = pd.DataFrame(annual_means)

    # Calculate linear trend
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        ts_df['year'], ts_df['mean_score']
    )
    trend_line = slope * ts_df['year'] + intercept

    # Create figure
    fig, ax = plt.subplots(figsize=Config.FIGURE_SIZE)

    # Plot annual means
    ax.plot(ts_df['year'], ts_df['mean_score'],
            color=Config.LINE_COLOR, linewidth=2, marker='o',
            markersize=5, label='Annual Mean')

    # Plot trend line
    trend_label = f'Trend: {slope:+.3f}/yr (p={p_value:.3f})'
    ax.plot(ts_df['year'], trend_line,
            color=Config.TREND_COLOR, linewidth=2, linestyle='--',
            label=trend_label)

    # Formatting
    ax.set_xlabel('Year', fontsize=11)
    ax.set_ylabel(f'Mean {score_label} Extreme Score', fontsize=11)
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Set x-axis ticks
    ax.set_xlim(Config.START_YEAR - 1, Config.END_YEAR + 1)
    ax.set_xticks(range(1990, 2025, 5))

    plt.tight_layout()

    # Save as PDF only
    output_path = Config.OUTPUT_DIR / 'ED_fig3d.pdf'
    fig.savefig(output_path, dpi=Config.DPI, bbox_inches='tight')
    print(f"  Saved: {output_path}")

    # Print summary statistics
    print(f"\n  Summary Statistics ({score_label}):")
    print(f"    Trend slope: {slope:+.4f} per year")
    print(f"    R-squared: {r_value**2:.4f}")
    print(f"    P-value: {p_value:.6f}")

    # Period comparisons
    early_mean = ts_df[ts_df['year'] < 2010]['mean_score'].mean()
    recent_mean = ts_df[ts_df['year'] >= 2010]['mean_score'].mean()
    print(f"    Mean (1990-2009): {early_mean:.3f}")
    print(f"    Mean (2010-2024): {recent_mean:.3f}")
    print(f"    Difference: {recent_mean - early_mean:+.3f}")

    plt.close(fig)
    return ts_df, slope, p_value


# =============================================================================
# FIGURE 3B: BAR CHART OF TOP-5 DISTRIBUTION (FREQUENCY)
# =============================================================================

def create_frequency_barchart(df, glacier_id_col, year_cols, score_type):
    """
    Create bar chart showing distribution of top-N extreme years by period.

    For each glacier, find its top N extreme years, then count how many
    of those years fall in each 5-year period.
    """
    score_label = score_type.capitalize()
    print(f"\nCreating {score_label} Frequency Bar Chart...")

    n_glaciers = len(df)
    top_n = Config.TOP_N_EXTREME_YEARS

    # For each glacier, find its top N extreme years
    period_counts = {f'{s}-{e}': 0 for s, e in Config.PERIODS}
    total_top_n_selections = 0

    for idx, row in df.iterrows():
        # Get scores for this glacier
        scores = []
        for year in year_cols:
            year_col = year if year in df.columns else str(year)
            if year_col in df.columns:
                val = row[year_col]
                if pd.notna(val):
                    scores.append((year, val))

        if len(scores) < top_n:
            continue

        # Sort by score (descending) and get top N years
        scores_sorted = sorted(scores, key=lambda x: x[1], reverse=True)
        top_years = [s[0] for s in scores_sorted[:top_n]]

        # Count which period each top year falls into
        for year in top_years:
            for start, end in Config.PERIODS:
                if start <= year <= end:
                    period_counts[f'{start}-{end}'] += 1
                    total_top_n_selections += 1
                    break

    # Calculate percentages
    period_stats = []
    for start, end in Config.PERIODS:
        key = f'{start}-{end}'
        count = period_counts[key]
        pct = (count / total_top_n_selections * 100) if total_top_n_selections > 0 else 0
        period_stats.append({
            'period': key,
            'start': start,
            'end': end,
            'count': count,
            'percentage': pct
        })

    period_df = pd.DataFrame(period_stats)

    # Expected percentage if random
    n_periods = len(Config.PERIODS)
    expected_pct = 100 / n_periods

    # Create figure
    fig, ax = plt.subplots(figsize=Config.FIGURE_SIZE)

    # Color bars: recent periods in red
    colors = []
    for _, row in period_df.iterrows():
        if row['start'] >= Config.RECENT_PERIOD_START:
            colors.append(Config.BAR_COLOR_RECENT)
        else:
            colors.append(Config.BAR_COLOR_EARLY)

    # Plot bars
    bars = ax.bar(period_df['period'], period_df['percentage'],
                  color=colors, edgecolor='black', linewidth=0.5)

    # Add value labels on bars
    for bar, pct in zip(bars, period_df['percentage']):
        height = bar.get_height()
        ax.annotate(f'{pct:.1f}%',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=9)

    # Add horizontal line at expected percentage
    ax.axhline(y=expected_pct, color='gray', linestyle='--', linewidth=1.5,
               label=f'Expected if random ({expected_pct:.1f}%)')

    # Formatting
    ax.set_xlabel('Time Period', fontsize=11)
    ax.set_ylabel(f'% of Top-{top_n} {score_label} Extreme Years', fontsize=11)
    ax.legend(loc='upper left', fontsize=10)

    # Rotate x-axis labels for readability
    plt.xticks(rotation=45, ha='right')

    # Set y-axis to start at 0
    ax.set_ylim(0, max(period_df['percentage']) * 1.25)

    plt.tight_layout()

    # Save as PDF only
    output_path = Config.OUTPUT_DIR / 'ED_fig3e.pdf'
    fig.savefig(output_path, dpi=Config.DPI, bbox_inches='tight')
    print(f"  Saved: {output_path}")

    # Print summary
    print(f"\n  Period Statistics ({score_label}, Top-{top_n} extreme years):")
    print(f"    Total selections: {total_top_n_selections:,}")
    print(f"    Expected per period (random): {expected_pct:.1f}%")
    print()
    for _, row in period_df.iterrows():
        marker = "**" if row['percentage'] > expected_pct * 1.5 else ""
        print(f"    {row['period']}: {row['percentage']:5.1f}% ({row['count']:,} selections) {marker}")

    plt.close(fig)
    return period_df


# =============================================================================
# MAIN
# =============================================================================

def process_score_type(score_type):
    """Process all figures for one score type."""
    print(f"\n{'='*60}")
    print(f"Processing: {score_type.upper()}")
    print('='*60)

    # Load data
    df, glacier_id_col, year_cols = load_extreme_scores(score_type)

    # Generate figures
    ts_df, slope, p_value = create_timeseries_figure(df, glacier_id_col, year_cols, score_type)
    period_df = create_frequency_barchart(df, glacier_id_col, year_cols, score_type)

    return {
        'score_type': score_type,
        'slope': slope,
        'p_value': p_value,
        'period_df': period_df
    }


def main():
    """Generate all Figure Batch 3 outputs."""
    print("=" * 60)
    print("FIGURE BATCH 3: SUMMARY TREND PLOTS (ED_fig3d-e)")
    print("=" * 60)

    # Create output directory
    Config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Process each score type
    results = {}
    for score_type in Config.SCORE_TYPES.keys():
        results[score_type] = process_score_type(score_type)

    # Summary for paper
    print("\n" + "=" * 60)
    print("KEY FINDINGS FOR PAPER:")
    print("=" * 60)

    for score_type, result in results.items():
        score_label = score_type.capitalize()
        slope = result['slope']
        p_value = result['p_value']
        period_df = result['period_df']

        sig = "significant" if p_value < 0.05 else "not significant"
        early_pct = period_df[period_df['start'] < 2010]['percentage'].sum()
        recent_pct = period_df[period_df['start'] >= 2010]['percentage'].sum()

        print(f"\n{score_label.upper()}:")
        print(f"  Intensity: slope={slope:+.4f}/yr, p={p_value:.4f} ({sig})")
        print(f"  Frequency: Early={early_pct:.1f}%, Recent={recent_pct:.1f}%")

    print("\n" + "=" * 60)
    print("Generated 2 figures: ED_fig3d.pdf, ED_fig3e.pdf")
    print("COMPLETED")
    print("=" * 60)


if __name__ == '__main__':
    main()
