"""
00h2: AGGREGATE ERA5 TEMPERATURE TO REGIONS (INTERMEDIATE AND IMMEDIATE)
=========================================================================
Spatially aggregate gridded ERA5 temperature data to Brazilian geographic regions.

Input: 
- Input_data/era5_brazil_hourly_*.nc (hourly, 0.25° grid)
- brazil_municipalities_2022.gpkg (municipality boundaries)
- results/municipality_region_mapping.csv (municipality to region mapping)

Output:
- results/era5_intermediate_daily.parquet (133 intermediate regions × ~5,500 days)
- results/era5_immediate_daily.parquet (510 immediate regions × ~5,500 days)

Method:
1. Load each yearly ERA5 file
2. Assign each grid cell to nearest municipality centroid (using KD-tree)
3. Aggregate municipality-level to both intermediate (133) and immediate (510) regions
4. Compute daily statistics: mean, min, max temperature

This creates temperature data at the same levels as our mortality analysis.
"""

import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np
import xarray as xr

print("="*70)
print("00h2: AGGREGATE ERA5 TO INTERMEDIATE AND IMMEDIATE REGIONS")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# PATHS
# =============================================================================

BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR.parent.parent.parent / 'Input_data'  # sim_data/Input_data
OUTPUT_DIR = BASE_DIR.parent / 'results'  # Main results folder
MAIN_RESULTS = BASE_DIR.parent / 'results'

STATES_FILE = INPUT_DIR / 'brazil_municipalities_2022.gpkg'
MAPPING_FILE = MAIN_RESULTS / 'municipality_to_all_regions_map.csv'  # Both intermediate and immediate codes

# =============================================================================
# CHECK DEPENDENCIES
# =============================================================================

try:
    import geopandas as gpd
    from scipy.spatial import cKDTree
    print("✓ geopandas, scipy installed")
except ImportError as e:
    print(f"✗ Missing dependency: {e}")
    exit(1)

# =============================================================================
# LOAD MUNICIPALITY BOUNDARIES AND REGION MAPPING
# =============================================================================

print("\n" + "-"*70)
print("Loading Geographic Data")
print("-"*70)

# Load municipality boundaries
print("Loading municipality boundaries...")
gdf_munic = gpd.read_file(STATES_FILE)
print(f"  Loaded {len(gdf_munic)} municipalities")

# Ensure CRS is WGS84 (EPSG:4326) for consistency with ERA5 lon/lat grid
if gdf_munic.crs is None:
    print("  ⚠ No CRS defined, assuming EPSG:4326")
    gdf_munic = gdf_munic.set_crs('EPSG:4326')
elif gdf_munic.crs.to_epsg() != 4326:
    print(f"  Converting CRS from {gdf_munic.crs} to EPSG:4326")
    gdf_munic = gdf_munic.to_crs('EPSG:4326')
else:
    print(f"  ✓ CRS verified: {gdf_munic.crs}")

# Get municipality centroids using representative_point (safer for multi-part geometries)
gdf_munic['centroid'] = gdf_munic.geometry.representative_point()
centroids = np.array([[p.x, p.y] for p in gdf_munic['centroid']])
muni_codes = gdf_munic['code_muni'].values

# Build KD-tree for fast nearest-neighbor lookup
tree = cKDTree(centroids)
print(f"  Built KD-tree with {len(centroids)} municipality centroids")

# Load region mapping
df_mapping = pd.read_csv(MAPPING_FILE)

# Create mappings for both levels
# Use 7-digit codes (matching geopackage)
muni_to_intermediate = dict(zip(df_mapping['code_muni'], df_mapping['intermediate_code']))
muni_to_immediate = dict(zip(df_mapping['code_muni'], df_mapping['immediate_code']))

print(f"  Loaded mapping: {len(df_mapping)} municipalities")
print(f"    → {df_mapping['intermediate_code'].nunique()} intermediate regions")
print(f"    → {df_mapping['immediate_code'].nunique()} immediate regions")

# =============================================================================
# GET LIST OF ERA5 FILES
# =============================================================================

print("\n" + "-"*70)
print("Finding ERA5 Files")
print("-"*70)

era5_files = sorted(INPUT_DIR.glob('era5_brazil_hourly_*.nc'))
print(f"Found {len(era5_files)} ERA5 annual files")
for f in era5_files[:3]:
    print(f"  - {f.name}")
if len(era5_files) > 3:
    print(f"  ... and {len(era5_files)-3} more")

# =============================================================================
# PROCESS ERA5 FILES
# =============================================================================

print("\n" + "-"*70)
print("Processing ERA5 Files")
print("-"*70)

all_daily = []
grid_region_map = None  # Will be computed once and reused

for i, era5_file in enumerate(era5_files):
    year = int(era5_file.stem.split('_')[-1])
    print(f"\n[{i+1}/{len(era5_files)}] Processing {year}...")
    
    try:
        ds = xr.open_dataset(era5_file)
        
        # Get coordinates
        lats = ds['latitude'].values
        lons = ds['longitude'].values
        
        # Handle different time variable names (time vs valid_time)
        if 'time' in ds.coords or 'time' in ds.dims:
            times = pd.to_datetime(ds['time'].values)
            time_dim = 'time'
        elif 'valid_time' in ds.coords or 'valid_time' in ds.dims:
            times = pd.to_datetime(ds['valid_time'].values)
            time_dim = 'valid_time'
        elif 'times' in ds.coords or 'times' in ds.dims:  # Additional fallback
            times = pd.to_datetime(ds['times'].values)
            time_dim = 'times'
        else:
            print(f"  ✗ Could not find time coordinate. Available: {list(ds.coords)}")
            print(f"     Data variables: {list(ds.data_vars)}")
            continue
        
        print(f"  Grid: {len(lats)}×{len(lons)} = {len(lats)*len(lons)} points")
        print(f"  Time steps: {len(times)}")
        
        # Build grid-to-region mappings (only first time)
        if grid_region_map is None:
            print("  Building grid-to-region mappings...")
            lon_grid, lat_grid = np.meshgrid(lons, lats)
            grid_points = np.column_stack([lon_grid.ravel(), lat_grid.ravel()])
            
            # Assign each grid point to nearest municipality
            distances, indices = tree.query(grid_points)
            grid_muni_codes = muni_codes[indices]
            
            # Distance diagnostics (flag grid points far from any municipality)
            dist_quantiles = np.percentile(distances, [50, 90, 95, 99])
            print(f"  Grid-to-muni distances (degrees): p50={dist_quantiles[0]:.3f}, p90={dist_quantiles[1]:.3f}, p95={dist_quantiles[2]:.3f}, p99={dist_quantiles[3]:.3f}")
            far_points = distances > 0.5  # >0.5 degrees (~55km) is suspicious
            if far_points.sum() > 0:
                print(f"  ⚠ Warning: {far_points.sum()} grid points are >0.5° from nearest municipality")
            
            # Create mappings for both levels
            grid_intermediate_map = np.array([muni_to_intermediate.get(int(m), -1) for m in grid_muni_codes])
            grid_immediate_map = np.array([muni_to_immediate.get(int(m), -1) for m in grid_muni_codes])
            grid_region_map = {'intermediate': grid_intermediate_map, 'immediate': grid_immediate_map}
            
            valid_intermediate = grid_intermediate_map >= 0
            valid_immediate = grid_immediate_map >= 0
            print(f"  Valid grid points: intermediate={valid_intermediate.sum()}/{len(grid_intermediate_map)}")
            print(f"                     immediate={valid_immediate.sum()}/{len(grid_immediate_map)}")
        
        # Get temperature variable
        if 't2m' in ds.data_vars:
            temp_var = 't2m'
        elif '2m_temperature' in ds.data_vars:
            temp_var = '2m_temperature'
        else:
            temp_var = list(ds.data_vars)[0]
            print(f"  ⚠ Using variable: {temp_var}")
        
        # Check if dewpoint is available
        dewpoint_var = None
        if 'd2m' in ds.data_vars:
            dewpoint_var = 'd2m'
        elif '2m_dewpoint_temperature' in ds.data_vars:
            dewpoint_var = '2m_dewpoint_temperature'
        
        # Process in daily chunks to manage memory
        dates = pd.Series(times).dt.date.unique()
        print(f"  Processing {len(dates)} days...")
        
        for date in dates:
            day_mask = pd.Series(times).dt.date == date
            day_indices = np.where(day_mask)[0]
            
            # Get hourly temperatures for this day (use dynamic time dimension)
            temp_day = ds[temp_var].isel({time_dim: day_indices}).values  # shape: (hours, lat, lon)
            temp_flat = temp_day.reshape(len(day_indices), -1)  # shape: (hours, grid_points)
            
            # Convert from Kelvin to Celsius if needed
            # First check units attribute, then use heuristic fallback
            temp_units = ds[temp_var].attrs.get('units', '').lower()
            if 'k' in temp_units or 'kelvin' in temp_units:
                temp_flat = temp_flat - 273.15
            elif temp_flat.mean() > 200:  # Heuristic fallback for missing/unclear units
                print(f"    ⚠ Day {date}: units unclear ('{temp_units}'), converting from K to C based on values")
                temp_flat = temp_flat - 273.15
            
            # Get dewpoint if available
            if dewpoint_var:
                dew_day = ds[dewpoint_var].isel({time_dim: day_indices}).values
                dew_flat = dew_day.reshape(len(day_indices), -1)
                if dew_flat.mean() > 200:
                    dew_flat = dew_flat - 273.15
            
            # Compute daily stats per grid point
            temp_mean = np.nanmean(temp_flat, axis=0)
            temp_min = np.nanmin(temp_flat, axis=0)
            temp_max = np.nanmax(temp_flat, axis=0)
            
            # Create dataframe for aggregation with both region levels
            df_grid = pd.DataFrame({
                'intermediate_code': grid_region_map['intermediate'],
                'immediate_code': grid_region_map['immediate'],
                'temp_mean': temp_mean,
                'temp_min': temp_min,
                'temp_max': temp_max
            })
            
            if dewpoint_var:
                df_grid['dewpoint_mean'] = np.nanmean(dew_flat, axis=0)
            
            # Aggregate to INTERMEDIATE region level
            df_intermediate = df_grid[df_grid['intermediate_code'] >= 0].groupby('intermediate_code').agg({
                'temp_mean': 'mean',
                'temp_min': 'mean',
                'temp_max': 'mean',
                **({'dewpoint_mean': 'mean'} if dewpoint_var else {})
            }).reset_index()
            df_intermediate['date'] = date
            df_intermediate['region_type'] = 'intermediate'
            
            # Aggregate to IMMEDIATE region level
            df_immediate = df_grid[df_grid['immediate_code'] >= 0].groupby('immediate_code').agg({
                'temp_mean': 'mean',
                'temp_min': 'mean',
                'temp_max': 'mean',
                **({'dewpoint_mean': 'mean'} if dewpoint_var else {})
            }).reset_index()
            df_immediate['date'] = date
            df_immediate['region_type'] = 'immediate'
            
            # Rename columns for consistency
            df_intermediate = df_intermediate.rename(columns={'intermediate_code': 'region_code'})
            df_immediate = df_immediate.rename(columns={'immediate_code': 'region_code'})
            
            all_daily.append(('intermediate', df_intermediate))
            all_daily.append(('immediate', df_immediate))
        
        ds.close()
        print(f"  ✓ Completed {year}")
        
    except Exception as e:
        print(f"  ✗ Error processing {era5_file.name}: {e}")
        import traceback
        traceback.print_exc()
        continue

# =============================================================================
# COMBINE AND SAVE
# =============================================================================

print("\n" + "-"*70)
print("Combining Results")
print("-"*70)

if all_daily:
    # Separate by region type
    intermediate_dfs = [df for level, df in all_daily if level == 'intermediate']
    immediate_dfs = [df for level, df in all_daily if level == 'immediate']
    
    # Process intermediate regions
    print("\n--- INTERMEDIATE REGIONS (133) ---")
    df_intermediate = pd.concat(intermediate_dfs, ignore_index=True)
    cols = ['date', 'region_code', 'temp_mean', 'temp_min', 'temp_max']
    if 'dewpoint_mean' in df_intermediate.columns:
        cols.append('dewpoint_mean')
    df_intermediate = df_intermediate[cols]
    df_intermediate = df_intermediate.sort_values(['region_code', 'date']).reset_index(drop=True)
    
    print(f"Total observations: {len(df_intermediate):,}")
    print(f"Regions: {df_intermediate['region_code'].nunique()}")
    print(f"Date range: {df_intermediate['date'].min()} to {df_intermediate['date'].max()}")
    print(f"Days per region: {len(df_intermediate) / df_intermediate['region_code'].nunique():.0f}")
    print(f"Temperature: mean={df_intermediate['temp_mean'].mean():.1f}°C, range={df_intermediate['temp_min'].min():.1f} to {df_intermediate['temp_max'].max():.1f}°C")
    
    output_intermediate = OUTPUT_DIR / 'era5_intermediate_daily.parquet'
    df_intermediate.to_parquet(output_intermediate, index=False)
    print(f"✓ Saved: {output_intermediate}")
    
    # Process immediate regions
    print("\n--- IMMEDIATE REGIONS (510) ---")
    df_immediate = pd.concat(immediate_dfs, ignore_index=True)
    cols = ['date', 'region_code', 'temp_mean', 'temp_min', 'temp_max']
    if 'dewpoint_mean' in df_immediate.columns:
        cols.append('dewpoint_mean')
    df_immediate = df_immediate[cols]
    df_immediate = df_immediate.sort_values(['region_code', 'date']).reset_index(drop=True)
    
    print(f"Total observations: {len(df_immediate):,}")
    print(f"Regions: {df_immediate['region_code'].nunique()}")
    print(f"Date range: {df_immediate['date'].min()} to {df_immediate['date'].max()}")
    print(f"Days per region: {len(df_immediate) / df_immediate['region_code'].nunique():.0f}")
    print(f"Temperature: mean={df_immediate['temp_mean'].mean():.1f}°C, range={df_immediate['temp_min'].min():.1f} to {df_immediate['temp_max'].max():.1f}°C")
    
    output_immediate = OUTPUT_DIR / 'era5_immediate_daily.parquet'
    df_immediate.to_parquet(output_immediate, index=False)
    print(f"✓ Saved: {output_immediate}")
    
    # Save CSV samples for inspection
    df_intermediate.head(1000).to_csv(OUTPUT_DIR / 'era5_intermediate_daily_sample.csv', index=False)
    df_immediate.head(1000).to_csv(OUTPUT_DIR / 'era5_immediate_daily_sample.csv', index=False)
    print(f"\n✓ Saved CSV samples for both levels")
else:
    print("✗ No data processed!")

print(f"\n{'='*70}")
print("DONE!")
print("="*70)
