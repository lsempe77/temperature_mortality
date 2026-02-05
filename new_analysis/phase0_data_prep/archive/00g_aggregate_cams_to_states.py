"""
00g: AGGREGATE CAMS POLLUTION DATA TO STATE-DAY LEVEL
======================================================
Spatially aggregate gridded CAMS pollution data to Brazilian states.

Input: 
- results/cams_pollution_brazil_daily.parquet (0.75° grid)
- Brazil state boundaries (for spatial join)

Output:
- results/cams_pollution_state_daily.parquet

Method:
- Load state boundaries from geopackage
- Assign each CAMS grid cell to its state (by centroid)
- Compute population-weighted or simple mean by state-day

For use in Phase 3 (pollution confounding) analysis.
"""

import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np

print("="*70)
print("00g: AGGREGATE CAMS POLLUTION DATA TO STATE-DAY LEVEL")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# PATHS
# =============================================================================

BASE_DIR = Path(__file__).parent
INPUT_DIR = Path(__file__).parent.parent.parent / 'Input_data'
OUTPUT_DIR = BASE_DIR / 'results'

CAMS_FILE = OUTPUT_DIR / 'cams_pollution_brazil_daily.parquet'
STATES_FILE = INPUT_DIR / 'brazil_municipalities_2022.gpkg'
OUTPUT_FILE = OUTPUT_DIR / 'cams_pollution_state_daily.parquet'

# =============================================================================
# LOAD DATA
# =============================================================================

print("\n" + "-"*70)
print("Loading Data")
print("-"*70)

# Load CAMS data
df_cams = pd.read_parquet(CAMS_FILE)
print(f"✓ CAMS data: {df_cams.shape[0]:,} observations")
print(f"  Grid points: {len(df_cams.drop_duplicates(['latitude', 'longitude'])):,}")
print(f"  Date range: {df_cams['date'].min()} to {df_cams['date'].max()}")

# Get unique grid points
grid_points = df_cams[['latitude', 'longitude']].drop_duplicates().reset_index(drop=True)
print(f"  Unique grid points: {len(grid_points)}")

# =============================================================================
# ASSIGN GRID POINTS TO STATES
# =============================================================================

print("\n" + "-"*70)
print("Assigning Grid Points to States")
print("-"*70)

try:
    import geopandas as gpd
    from shapely.geometry import Point
    
    # Load state boundaries
    print("Loading state boundaries...")
    gdf_munic = gpd.read_file(STATES_FILE)
    print(f"  Loaded {len(gdf_munic)} municipalities")
    
    # Get state-level boundaries by dissolving
    # Column is 'abbrev_state' not 'sigla_uf'
    gdf_states = gdf_munic.dissolve(by='abbrev_state').reset_index()
    print(f"  Dissolved to {len(gdf_states)} states")
    
    # Create GeoDataFrame from grid points
    geometry = [Point(lon, lat) for lat, lon in zip(grid_points['latitude'], grid_points['longitude'])]
    gdf_grid = gpd.GeoDataFrame(grid_points, geometry=geometry, crs="EPSG:4326")
    
    # Spatial join to find which state each point is in
    print("Performing spatial join...")
    gdf_grid_states = gpd.sjoin(gdf_grid, gdf_states[['abbrev_state', 'geometry']], 
                                 how='left', predicate='within')
    
    # Create mapping dictionary
    grid_state_map = dict(zip(
        zip(gdf_grid_states['latitude'], gdf_grid_states['longitude']),
        gdf_grid_states['abbrev_state']
    ))
    
    print(f"  Assigned {gdf_grid_states['abbrev_state'].notna().sum()} of {len(gdf_grid)} points to states")
    
    # Points outside Brazil (ocean, neighboring countries)
    outside = gdf_grid_states['abbrev_state'].isna().sum()
    if outside > 0:
        print(f"  {outside} points outside Brazil boundaries (will be excluded)")
    
    USE_SPATIAL = True
    
except ImportError as e:
    print(f"⚠ GeoPandas not available: {e}")
    print("  Using approximate state assignment based on lat/lon...")
    USE_SPATIAL = False

# Fallback: simple lat/lon-based assignment
if not USE_SPATIAL:
    # Very rough state boundaries for Brazil
    # This is a simplified approach - use only if geopandas not available
    
    def assign_state_approx(lat, lon):
        """Approximate state assignment based on lat/lon."""
        # This is a VERY rough approximation - for demonstration only
        # Should use proper spatial join with state boundaries
        
        # Northern states
        if lat > 0:
            if lon < -60:
                return 'AM'  # Amazonas
            elif lon < -50:
                return 'PA'  # Pará
            else:
                return 'AP'  # Amapá
        
        # Northeast
        elif lat > -10:
            if lon > -40:
                return 'RN'  # Rio Grande do Norte
            elif lon > -45:
                return 'BA'  # Bahia
            else:
                return 'PI'  # Piauí
        
        # Southeast
        elif lat > -25:
            if lon > -45:
                return 'RJ'  # Rio de Janeiro
            elif lon > -50:
                return 'SP'  # São Paulo
            else:
                return 'MG'  # Minas Gerais
        
        # South
        else:
            if lon > -50:
                return 'SC'  # Santa Catarina
            else:
                return 'RS'  # Rio Grande do Sul
        
        return None
    
    grid_state_map = {
        (row['latitude'], row['longitude']): assign_state_approx(row['latitude'], row['longitude'])
        for _, row in grid_points.iterrows()
    }
    
    print("  ⚠ Using approximate state boundaries - install geopandas for accuracy")

# =============================================================================
# AGGREGATE TO STATE-DAY
# =============================================================================

print("\n" + "-"*70)
print("Aggregating to State-Day Level")
print("-"*70)

# Add state to CAMS data
df_cams['state'] = df_cams.apply(
    lambda row: grid_state_map.get((row['latitude'], row['longitude'])), 
    axis=1
)

# Remove observations outside Brazil
n_before = len(df_cams)
df_cams = df_cams[df_cams['state'].notna()]
n_after = len(df_cams)
print(f"  Removed {n_before - n_after:,} observations outside Brazil")

# Ensure date is proper
df_cams['date'] = pd.to_datetime(df_cams['date'])

# Aggregate to state-day level (simple mean across grid cells)
print("Computing state-day means...")

agg_cols = {
    'pm25_ugm3': 'mean',
    'pm10_ugm3': 'mean',
    'o3_column_kgm2': 'mean',
    'latitude': 'count',  # Count of grid cells
}

df_state_day = df_cams.groupby(['state', 'date']).agg(agg_cols).reset_index()
df_state_day = df_state_day.rename(columns={'latitude': 'n_grid_cells'})

# Rename for clarity
df_state_day = df_state_day.rename(columns={
    'pm25_ugm3': 'pm25_mean_ugm3',
    'pm10_ugm3': 'pm10_mean_ugm3',
    'o3_column_kgm2': 'o3_column_mean_kgm2',
    'state': 'state_abbrev'
})

print(f"\n✓ Aggregated data shape: {df_state_day.shape}")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

print("\n" + "-"*70)
print("Summary Statistics")
print("-"*70)

print(f"States: {df_state_day['state_abbrev'].nunique()}")
print(f"Date range: {df_state_day['date'].min()} to {df_state_day['date'].max()}")
print(f"Total state-days: {len(df_state_day):,}")

print("\nPollution statistics:")
print(f"  PM2.5: {df_state_day['pm25_mean_ugm3'].min():.2f} - {df_state_day['pm25_mean_ugm3'].max():.2f} µg/m³")
print(f"         Mean: {df_state_day['pm25_mean_ugm3'].mean():.2f} µg/m³")
print(f"  PM10:  {df_state_day['pm10_mean_ugm3'].min():.2f} - {df_state_day['pm10_mean_ugm3'].max():.2f} µg/m³")
print(f"         Mean: {df_state_day['pm10_mean_ugm3'].mean():.2f} µg/m³")

print("\nGrid cells per state:")
cells_per_state = df_cams.groupby('state')[['latitude', 'longitude']].apply(
    lambda x: len(x.drop_duplicates())
).reset_index(name='n_cells')
print(cells_per_state.sort_values('n_cells', ascending=False).head(10).to_string(index=False))

# =============================================================================
# SAVE OUTPUT
# =============================================================================

print("\n" + "-"*70)
print("Saving Output")
print("-"*70)

df_state_day.to_parquet(OUTPUT_FILE, index=False)
print(f"✓ Saved: {OUTPUT_FILE}")

# Also save as CSV for inspection
csv_file = OUTPUT_DIR / 'cams_pollution_state_daily.csv'
df_state_day.to_csv(csv_file, index=False)
print(f"✓ Saved: {csv_file}")

print(f"\n{'='*70}")
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*70)
