"""
Figure Batch 6: Distribution of Extreme Scores by Year (ED_fig3f)
==================================================================

Shows how the distribution of extreme scores across glaciers changes over time.

ED_fig3f: Box plot (side-by-side boxes, face colored by median score)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

# Font settings for Adobe Illustrator compatibility
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


# =============================================================================
# CONFIGURATION
# =============================================================================

class Config:
    """Configuration parameters for Figure Batch 6"""

    # --- Data period ---
    START_YEAR = 1990
    END_YEAR = 2024

    # --- Data paths ---
    GLACIER_SHAPEFILE = "/WORK/Data/ncc_glacier/glacier_boundary/rgi_hma/rgi_hma.shp"
    EXTREME_SCORES_XLSX = "/WORK/Data/ncc_glacier/era5_land/extracted/glacier_extreme_scores.xlsx"

    # --- Output ---
    OUTPUT_DIR = Path(__file__).parent / "output"
    DPI = 600

    # --- Styling ---
    BOXPLOT_FIGURE_SIZE = (18, 6)
    LABEL_ROTATION = 45
    COLORMAP = 'RdYlBu_r'  # Blue=low scores, Red=high scores


# =============================================================================
# DATA LOADING
# =============================================================================

def load_compound_scores():
    """Load compound extreme scores."""
    print("Loading compound extreme scores...")
    df = pd.read_excel(Config.EXTREME_SCORES_XLSX, sheet_name='Compound_Extreme_Scores')

    glacier_id_col = df.columns[0]
    year_cols = [c for c in df.columns if isinstance(c, int) or (isinstance(c, str) and c.isdigit())]
    year_cols = [int(c) for c in year_cols]
    year_cols = [y for y in year_cols if Config.START_YEAR <= y <= Config.END_YEAR]
    year_cols = sorted(year_cols)

    if not year_cols:
        raise ValueError(f"No valid year columns found between {Config.START_YEAR} and {Config.END_YEAR}")

    print(f"  Years: {year_cols[0]} to {year_cols[-1]} ({len(year_cols)} years)")
    print(f"  Glaciers: {len(df)}")
    return df, glacier_id_col, year_cols


def prepare_long_format(df, glacier_id_col, year_cols):
    """Convert wide format to long format for plotting."""
    print("Preparing data in long format...")

    records = []
    for year in year_cols:
        year_col = year if year in df.columns else str(year)
        if year_col in df.columns:
            scores = df[year_col].dropna()
            for score in scores:
                records.append({'year': year, 'score': score})

    long_df = pd.DataFrame(records)
    print(f"  Total records: {len(long_df):,}")
    return long_df


# =============================================================================
# VISUALIZATION: BOX PLOT
# =============================================================================

def create_box_plot(long_df, year_cols):
    """Create box plot of extreme scores by year (ED_fig3f)."""
    print("\nCreating box plot (ED_fig3f)...")

    years = sorted(year_cols)

    # Calculate median score per year for coloring
    medians = long_df.groupby('year')['score'].median().to_dict()

    # Determine score range for colormap normalization
    all_medians = list(medians.values())
    score_min, score_max = min(all_medians), max(all_medians)

    # Set up colormap (use plt.colormaps[] instead of deprecated get_cmap)
    cmap = plt.colormaps[Config.COLORMAP]
    norm = mcolors.Normalize(vmin=score_min, vmax=score_max)

    # Create figure
    fig, ax = plt.subplots(figsize=Config.BOXPLOT_FIGURE_SIZE)

    # Prepare box plot data
    box_data = [long_df[long_df['year'] == y]['score'].values for y in years]

    # Create box plot
    bp = ax.boxplot(box_data, positions=range(len(years)), patch_artist=True,
                    widths=0.7, showfliers=False)

    # Color each box by its median
    for i, (patch, year) in enumerate(zip(bp['boxes'], years)):
        median_val = medians[year]
        color = cmap(norm(median_val))
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
        patch.set_alpha(0.8)

    # Style median lines
    for median in bp['medians']:
        median.set_color('black')
        median.set_linewidth(1.5)

    # X-axis labels
    ax.set_xticks(range(len(years)))
    ax.set_xticklabels(years, rotation=Config.LABEL_ROTATION, ha='right', fontsize=8)

    # Y-axis
    ax.set_xlabel('Year', fontsize=11)
    ax.set_ylabel('Compound Extreme Score', fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label('Median Score', fontsize=10)

    plt.tight_layout()

    # Save as PDF only
    output_path = Config.OUTPUT_DIR / 'ED_fig3f.pdf'
    fig.savefig(output_path, dpi=Config.DPI, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Generate Figure Batch 6 output (ED_fig3f)."""
    print("=" * 60)
    print("FIGURE BATCH 6: EXTREME SCORES BY YEAR (ED_fig3f)")
    print("=" * 60)

    # Create output directory
    Config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    scores_df, glacier_id_col, year_cols = load_compound_scores()

    # Prepare long format data
    long_df = prepare_long_format(scores_df, glacier_id_col, year_cols)

    # Create box plot
    create_box_plot(long_df, year_cols)

    # Summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    # Print median scores by year (first 5 and last 5)
    medians = long_df.groupby('year')['score'].median().sort_index()
    print("\nMedian compound score by year:")
    print("  Early years:")
    for year in list(medians.index)[:5]:
        print(f"    {year}: {medians[year]:.3f}")
    print("  ...")
    print("  Recent years:")
    for year in list(medians.index)[-5:]:
        print(f"    {year}: {medians[year]:.3f}")

    # Trend (use already imported stats module)
    years_arr = np.array(list(medians.index))
    medians_arr = np.array(list(medians.values))
    slope, intercept, r_value, p_value, std_err = stats.linregress(years_arr, medians_arr)
    print(f"\n  Trend in median scores:")
    print(f"    Slope: {slope:+.4f} per year")
    print(f"    P-value: {p_value:.4f}")

    print("\n" + "=" * 60)
    print("COMPLETED")
    print("=" * 60)


if __name__ == '__main__':
    main()
