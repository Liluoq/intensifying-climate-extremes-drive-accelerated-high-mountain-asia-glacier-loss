"""
Figure Batch 1: Gridded Difference Maps (ED_fig3a-c)
====================================================

Creates three figures for Extended Data Figure 3:
- ED_fig3a: Temperature Change (Post-2021 vs Pre-2021)
- ED_fig3b: Precipitation Change (Post-2021 vs Pre-2021)
- ED_fig3c: Compound Extreme Score Difference

Styling matches reference notebook exactly: plot_grid_bubble_map()
from map_figure_style_ref/ICESat2_Visualization_Notebook_grid.ipynb
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# CONFIGURATION
# =============================================================================

class Config:
    """Configuration parameters for Figure Batch 1"""

    # --- Time period division ---
    THRESHOLD_YEAR = 2021
    EARLY_PERIOD = (1990, 2020)  # inclusive
    RECENT_PERIOD = (2021, 2024)  # inclusive

    # --- Grid resolution ---
    GRID_SIZE_DEG = 1.0  # 1° × 1° grid cells

    # --- Data paths ---
    GLACIER_SHAPEFILE = "/WORK/Data/ncc_glacier/glacier_boundary/rgi_hma/rgi_hma.shp"
    EXTREME_SCORES_XLSX = "/WORK/Data/ncc_glacier/era5_land/extracted/glacier_extreme_scores.xlsx"
    MONTHLY_STATS_XLSX = "/WORK/Data/ncc_glacier/era5_land/extracted/glacier_stats_monthly.xlsx"

    # --- Boundary shapefiles (USER: Update these paths) ---
    HMA_BOUNDARY_SHP = "/WORK/Data/ncc_glacier/regional_boundary/hma/HMA_one.shp"  # e.g., "/WORK/Data/HMA/HMA_boundary/HMA_one.shp"
    SUBREGION_BOUNDARY_SHP = "/WORK/Data/ncc_glacier/regional_boundary/subregions/boundary_mountain_regions_hma_v3_zheng_20200601.shp"  # e.g., subregion boundary shapefile

    # --- Output ---
    OUTPUT_DIR = Path(__file__).parent / "output"
    DPI = 600  # Matches reference notebook

    # --- Colorbar settings ---
    # Set custom discrete bins for each figure (None = continuous colormap)
    # Bins define the boundaries, e.g., [-1.5, -1, -0.5, 0, 0.5, 1, 1.5] creates 6 color levels
    TEMP_DIFF_BINS = [-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5]
    PRECIP_DIFF_BINS = [-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50]  # mm difference
    TEMP_SCORE_BINS = [-3, -2, -1, 0, 1, 2, 3]
    PRECIP_SCORE_BINS = [-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5]
    COMPOUND_SCORE_BINS = [-5, -3, -1, 0, 1, 3, 5]

    # Colorbar extend: 'neither', 'both', 'min', 'max'
    CBAR_EXTEND = 'both'

    # --- Unit conversion ---
    # ERA5-Land precipitation is in m/day (monthly mean of daily values)
    # Convert to mm/year: m/day × 1000 (m→mm) × 365 (days/year)
    PRECIP_M_PER_DAY_TO_MM_PER_YEAR = 1000 * 365

    # --- Map projection (Albers Equal Area for HMA) - matches reference exactly ---
    CENTRAL_LON = 87.5
    CENTRAL_LAT = 34.7
    STANDARD_PARALLELS = (24.5, 47.5)
    EXTENT = [63.6, 108.3, 24.5, 47.5]  # [lon_min, lon_max, lat_min, lat_max]


# =============================================================================
# DATA LOADING AND PROCESSING
# =============================================================================

def load_glacier_data():
    """Load glacier shapefile with coordinates and area."""
    print("Loading glacier shapefile...")
    gdf = gpd.read_file(Config.GLACIER_SHAPEFILE)
    glaciers = gdf[['RGIId', 'CenLon', 'CenLat', 'Area']].copy()
    glaciers = glaciers.rename(columns={'RGIId': 'glacier_id'})
    print(f"  Loaded {len(glaciers)} glaciers")
    return glaciers


def load_extreme_scores():
    """Load annual extreme scores from Excel file."""
    print("Loading extreme scores...")
    scores = {}

    for sheet, key in [
        ('Temperature_Extreme_Scores', 'temp_score'),
        ('Precipitation_Extreme_Scores', 'precip_score'),
        ('Compound_Extreme_Scores', 'compound_score')
    ]:
        df = pd.read_excel(Config.EXTREME_SCORES_XLSX, sheet_name=sheet)
        scores[key] = df
        print(f"  {sheet}: {df.shape[0]} glaciers, years {df.columns[1]}-{df.columns[-1]}")

    return scores


def load_monthly_temperature():
    """Load monthly temperature data."""
    print("Loading monthly temperature data...")
    df = pd.read_excel(Config.MONTHLY_STATS_XLSX, sheet_name='temperature_2m')
    print(f"  Shape: {df.shape}")
    return df


def load_monthly_precipitation():
    """Load monthly precipitation data."""
    print("Loading monthly precipitation data...")
    df = pd.read_excel(Config.MONTHLY_STATS_XLSX, sheet_name='total_precipitation_sum')
    print(f"  Shape: {df.shape}")
    return df


def load_boundaries():
    """Load boundary shapefiles if available."""
    hma_boundary = None
    subregion_boundary = None

    if Config.HMA_BOUNDARY_SHP and Path(Config.HMA_BOUNDARY_SHP).exists():
        hma_boundary = gpd.read_file(Config.HMA_BOUNDARY_SHP)
        print(f"  Loaded HMA boundary")

    if Config.SUBREGION_BOUNDARY_SHP and Path(Config.SUBREGION_BOUNDARY_SHP).exists():
        subregion_boundary = gpd.read_file(Config.SUBREGION_BOUNDARY_SHP)
        print(f"  Loaded subregion boundaries")

    return hma_boundary, subregion_boundary


def calculate_period_means(df, value_columns, period):
    """Calculate mean values for a time period."""
    start_year, end_year = period

    period_cols = []
    for col in value_columns:
        if isinstance(col, (int, float)):
            year = int(col)
        else:
            year = int(str(col).split('-')[0])

        if start_year <= year <= end_year:
            period_cols.append(col)

    if not period_cols:
        raise ValueError(f"No columns found for period {period}")

    return df[period_cols].mean(axis=1)


def create_grid_cells(glaciers, grid_size=1.0):
    """Assign glaciers to grid cells based on their coordinates."""
    glaciers = glaciers.copy()
    glaciers['grid_lon'] = np.floor(glaciers['CenLon'] / grid_size) * grid_size
    glaciers['grid_lat'] = np.floor(glaciers['CenLat'] / grid_size) * grid_size

    # Cell center for plotting (x, y coordinates like reference notebook)
    glaciers['x'] = glaciers['grid_lon'] + grid_size / 2
    glaciers['y'] = glaciers['grid_lat'] + grid_size / 2

    return glaciers


def aggregate_to_grid(data, glaciers, var_col, agg_func='mean'):
    """Aggregate glacier-level data to grid cells."""
    merged = glaciers.merge(data[['glacier_id', var_col]], on='glacier_id', how='left')

    grid_stats = merged.groupby(['grid_lon', 'grid_lat', 'x', 'y']).agg({
        'Area': 'sum',
        var_col: agg_func,
        'glacier_id': 'count'
    }).reset_index()

    grid_stats = grid_stats.rename(columns={'glacier_id': 'n_glaciers'})

    return grid_stats


def prepare_figure_data():
    """Prepare all data needed for the figures."""
    glaciers = load_glacier_data()
    scores = load_extreme_scores()
    temp_monthly = load_monthly_temperature()
    precip_monthly = load_monthly_precipitation()

    glaciers = create_grid_cells(glaciers, Config.GRID_SIZE_DEG)

    figure_data = {}

    # --- Temperature Difference ---
    print("\nProcessing temperature difference...")
    value_cols = [c for c in temp_monthly.columns if c != 'glacier_id']

    early_mean = calculate_period_means(temp_monthly, value_cols, Config.EARLY_PERIOD)
    recent_mean = calculate_period_means(temp_monthly, value_cols, Config.RECENT_PERIOD)

    temp_diff = pd.DataFrame({
        'glacier_id': temp_monthly['glacier_id'],
        'value': recent_mean - early_mean
    })

    grid_temp = aggregate_to_grid(temp_diff, glaciers, 'value', 'mean')
    figure_data['temp_diff'] = grid_temp
    print(f"  Grid cells: {len(grid_temp)}, Range: {grid_temp['value'].min():.2f} to {grid_temp['value'].max():.2f} K")

    # --- Precipitation Difference ---
    # Note: Monthly precipitation data is in m/day (mean of daily values)
    # Convert to mm/year for meaningful interpretation
    print("\nProcessing precipitation difference...")
    value_cols = [c for c in precip_monthly.columns if c != 'glacier_id']

    early_mean = calculate_period_means(precip_monthly, value_cols, Config.EARLY_PERIOD)
    recent_mean = calculate_period_means(precip_monthly, value_cols, Config.RECENT_PERIOD)

    # Calculate difference and convert from m/day to mm/year
    precip_diff_raw = recent_mean - early_mean  # in m/day
    precip_diff_mm_year = precip_diff_raw * Config.PRECIP_M_PER_DAY_TO_MM_PER_YEAR

    precip_diff = pd.DataFrame({
        'glacier_id': precip_monthly['glacier_id'],
        'value': precip_diff_mm_year
    })

    grid_precip = aggregate_to_grid(precip_diff, glaciers, 'value', 'mean')
    figure_data['precip_diff'] = grid_precip
    print(f"  Grid cells: {len(grid_precip)}, Range: {grid_precip['value'].min():.1f} to {grid_precip['value'].max():.1f} mm/year")

    # --- Compound Extreme Score Difference ---
    print("\nProcessing compound extreme score difference...")
    compound_score_df = scores['compound_score']
    year_cols = [c for c in compound_score_df.columns if c != 'glacier_id']

    early_mean = calculate_period_means(compound_score_df, year_cols, Config.EARLY_PERIOD)
    recent_mean = calculate_period_means(compound_score_df, year_cols, Config.RECENT_PERIOD)

    compound_score_diff = pd.DataFrame({
        'glacier_id': compound_score_df['glacier_id'],
        'value': recent_mean - early_mean
    })

    grid_compound_score = aggregate_to_grid(compound_score_diff, glaciers, 'value', 'mean')
    figure_data['compound_score'] = grid_compound_score
    print(f"  Grid cells: {len(grid_compound_score)}, Range: {grid_compound_score['value'].min():.2f} to {grid_compound_score['value'].max():.2f}")

    return figure_data


# =============================================================================
# VISUALIZATION - Matches reference notebook EXACTLY
# =============================================================================

def plot_grid_bubble_map(grids, var, title_text, vmin, vmax, cbar_label,
                         studyregion=None, subregion=None,
                         save_path=None, colormap=None,
                         bins=None, extend='both'):
    """
    Plot grid-scale bubble map - EXACT COPY of reference notebook style.

    Parameters:
    -----------
    grids : DataFrame
        Grid data with 'x', 'y', 'Area', and var columns
    var : str
        Column name for the variable to plot
    title_text : str
        Figure title
    vmin, vmax : float
        Colormap value range (used for continuous colormap)
    cbar_label : str
        Colorbar label text
    studyregion : GeoDataFrame, optional
        HMA boundary
    subregion : GeoDataFrame, optional
        Subregion boundaries
    save_path : str, optional
        Output file path
    colormap : Colormap, optional
        Custom colormap (if None, uses RdBu_r)
    bins : list, optional
        Custom bin boundaries for discrete colormap, e.g., [-3, -2, -1, 0, 1, 2, 3]
        If None, uses continuous colormap with vmin/vmax
    extend : str, optional
        Colorbar extension: 'neither', 'both', 'min', 'max' (default: 'both')
    """
    # Set up fonts - EXACT match to reference (with fallback for Linux)
    try:
        import matplotlib.font_manager as fm
        if any('Arial' in f.name for f in fm.fontManager.ttflist):
            plt.rcParams['font.family'] = 'Arial'
        else:
            plt.rcParams['font.family'] = 'DejaVu Sans'  # Similar to Arial, available on Linux
    except:
        plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    # Basemap tiles - EXACT match to reference
    arcgis_url = ('https://server.arcgisonline.com/arcgis/rest/services/'
                  'World_Terrain_Base/MapServer/tile/{z}/{y}/{x}.jpg')
    tiles = cimgt.GoogleTiles(url=arcgis_url)

    # Figure size - EXACT match to reference: (18, 12)
    fig = plt.figure(figsize=(18, 12))

    # Projection - EXACT match to reference
    map_proj = ccrs.AlbersEqualArea(
        central_longitude=87.5,
        central_latitude=34.7,
        standard_parallels=(24.5, 47.5)
    )

    ax = fig.add_subplot(1, 1, 1, projection=map_proj)
    extent = [63.6, 108.3, 24.5, 47.5]
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    # Add basemap and features - EXACT match to reference
    ax.add_image(tiles, 9)  # zoom level 9
    ax.add_feature(cfeat.RIVERS.with_scale('10m'), linewidth=2, zorder=1)
    ax.add_feature(cfeat.LAKES.with_scale('10m'), zorder=1)

    # Plot boundaries - EXACT match to reference
    if studyregion is not None:
        studyregion.plot(ax=ax, color="none", edgecolor="maroon",
                         linewidth=3, transform=ccrs.PlateCarree(), zorder=2)

    if subregion is not None:
        subregion.plot(ax=ax, color="none", edgecolor="maroon",
                       linewidth=2, transform=ccrs.PlateCarree(), zorder=5)

    # Default colormap - diverging for difference maps
    if colormap is None:
        colormap = 'RdBu_r'

    # Get the colormap object
    if isinstance(colormap, str):
        cmap = plt.cm.get_cmap(colormap)
    else:
        cmap = colormap

    # Process data - EXACT match to reference approach
    grids_copy = grids.copy()
    mask = (grids_copy[var] == 0) | (pd.isna(grids_copy[var]))
    grids_copy.loc[mask, "Area"] = 0

    # Set up colormap normalization
    if bins is not None:
        # Discrete colormap with custom bins
        norm = mcolors.BoundaryNorm(bins, cmap.N, extend=extend)
        scatter_vmin, scatter_vmax = None, None  # norm handles the range
    else:
        # Continuous colormap with vmin/vmax
        norm = None
        scatter_vmin, scatter_vmax = vmin, vmax

    # Plot scatter - EXACT match to reference (s=Area directly, no scaling!)
    sc = plt.scatter(grids_copy['x'], grids_copy['y'], s=grids_copy["Area"], c=grids_copy[var],
                     transform=ccrs.PlateCarree(), cmap=cmap, norm=norm,
                     vmin=scatter_vmin, vmax=scatter_vmax,
                     edgecolor='k', zorder=2)

    ax.grid('both', linestyle='--', linewidth=0.5)

    # Text labels - EXACT positions and style from reference
    plt.text(x=63.5, y=31.3, s='Glacier area', rotation=0, ha='left',
             va='baseline', transform=ccrs.PlateCarree(),
             fontdict=dict(fontsize=25, color='k', family=plt.rcParams['font.family'], weight='normal'))

    plt.text(x=65.9, y=30.6, s='(km$^2$)', rotation=0, ha='left', va='baseline',
             transform=ccrs.PlateCarree(),
             fontdict=dict(fontsize=25, color='k', family=plt.rcParams['font.family'], weight='normal'))

    # Colorbar label (rotated, on right side)
    plt.text(x=107.6, y=38.5, s=cbar_label, rotation=90,
             ha='left', va='baseline', transform=ccrs.PlateCarree(),
             fontdict=dict(fontsize=25, color='k', family=plt.rcParams['font.family'], weight='normal'))

    # Title with white background - EXACT match to reference
    t = plt.text(x=56.7, y=45, s=title_text, rotation=0, ha='left', va='baseline',
                 transform=ccrs.PlateCarree(),
                 fontdict=dict(fontsize=27, color='k', family=plt.rcParams['font.family'], weight='normal'))
    t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='white'))

    # Gridlines - EXACT match to reference
    gl = ax.gridlines(draw_labels=True, linestyle=":", linewidth=0.3,
                      x_inline=False, y_inline=False, color='k')
    gl.top_labels = True
    gl.right_labels = True
    gl.bottom_labels = False
    gl.left_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(65, 108, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(25, 47, 5))
    gl.xlabel_style = {'size': 23}
    gl.ylabel_style = {'size': 23}

    # Colorbar - position from reference, with extend triangles
    position = fig.add_axes([0.75, 0.58, 0.015, 0.28])  # Slightly larger to fit triangles
    b = plt.colorbar(sc, cax=position, extend=extend,
                     orientation='vertical')
    ax2 = b.ax
    ax2.yaxis.set_ticks_position('left')
    ax2.tick_params(labelsize=23, left=True, right=True)

    # Set colorbar ticks at bin boundaries if discrete
    if bins is not None:
        b.set_ticks(bins)

    # Size legend - EXACT match to reference
    marker1 = ax.scatter([], [], s=50, c='k', alpha=0.3)
    marker2 = ax.scatter([], [], s=500, c='k', alpha=0.3)
    marker3 = ax.scatter([], [], s=1000, c='k', alpha=0.3)
    marker4 = ax.scatter([], [], s=2000, c='k', alpha=0.3)
    ax.legend(handles=[marker1, marker2, marker3, marker4],
              labels=['50', '500', '1000', '2000'],
              scatterpoints=1, frameon=False, fontsize=23,
              labelspacing=1.35, handletextpad=1,
              bbox_to_anchor=(0.198, 0.325), ncol=1,
              facecolor="white", edgecolor="None")

    # Save figure as PDF only
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
        print(f"  Saved: {save_path}")

    plt.show()

    return fig


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Main execution function."""
    print("=" * 60)
    print("Figure Batch 1: Gridded Difference Maps (ED_fig3a-c)")
    print("=" * 60)

    # Create output directory
    Config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Prepare data
    figure_data = prepare_figure_data()

    # Load boundaries
    print("\nLoading boundaries...")
    studyregion, subregion = load_boundaries()

    # Figure configurations
    # bins=None uses continuous colormap; bins=[...] uses discrete colormap
    figures = [
        {
            'data': figure_data['temp_diff'],
            'var': 'value',
            'title': f'a Temperature Change ({Config.RECENT_PERIOD[0]}-{Config.RECENT_PERIOD[1]} vs {Config.EARLY_PERIOD[0]}-{Config.EARLY_PERIOD[1]})',
            'vmin': -1.5, 'vmax': 1.5,
            'bins': None,  # Set to None for continuous
            'cbar_label': 'Temperature Diff (K)',
            'colormap': 'RdBu_r',
            'filename': 'ED_fig3a.pdf'
        },
        {
            'data': figure_data['precip_diff'],
            'var': 'value',
            'title': f'b Precipitation Change ({Config.RECENT_PERIOD[0]}-{Config.RECENT_PERIOD[1]} vs {Config.EARLY_PERIOD[0]}-{Config.EARLY_PERIOD[1]})',
            'vmin': -100, 'vmax': 100,
            'bins': None,  # Set to None for continuous
            'cbar_label': 'Precip Diff (mm)',
            'colormap': 'BrBG',  # Brown=drier (negative), Green=wetter (positive)
            'filename': 'ED_fig3b.pdf'
        },
        {
            'data': figure_data['compound_score'],
            'var': 'value',
            'title': f'c Compound Extreme Score Change',
            'vmin': -5, 'vmax': 5,
            'bins': Config.COMPOUND_SCORE_BINS,  # Set to None for continuous
            'cbar_label': 'Score Diff (+=more)',
            'colormap': 'RdYlBu_r',
            'filename': 'ED_fig3c.pdf'
        }
    ]

    # Generate each figure
    for i, fig_config in enumerate(figures):
        print(f"\n{'='*60}")
        print(f"Generating Figure {i+1}: {fig_config['title']}")
        print("=" * 60)

        save_path = Config.OUTPUT_DIR / fig_config['filename']

        plot_grid_bubble_map(
            grids=fig_config['data'],
            var=fig_config['var'],
            title_text=fig_config['title'],
            vmin=fig_config['vmin'],
            vmax=fig_config['vmax'],
            cbar_label=fig_config['cbar_label'],
            studyregion=studyregion,
            subregion=subregion,
            save_path=save_path,
            colormap=fig_config['colormap'],
            bins=fig_config['bins'],
            extend=Config.CBAR_EXTEND
        )

    print("\n" + "=" * 60)
    print("Complete! All figures saved to:", Config.OUTPUT_DIR)
    print("=" * 60)


if __name__ == "__main__":
    main()
