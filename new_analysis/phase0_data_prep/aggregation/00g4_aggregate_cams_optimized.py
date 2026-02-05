"""
00g4: AGGREGATE CAMS POLLUTION DATA - OPTIMIZED VERSION
=========================================================
Fast vectorized aggregation of CAMS pollution data to regions.

Key optimizations:
1. Process entire files at once (not timestep-by-timestep)
2. Use vectorized groupby operations (not apply)
3. Pre-compute aggregation weights once
4. Simple mean instead of weighted mean (area differences are small within regions)

Input: 
- temp_cams/*.zip (quarterly CAMS files, 0.75° grid)
- brazil_municipalities_2022.gpkg (municipality boundaries)
- results/municipality_to_all_regions_map.csv

Output:
- results/cams_intermediate_daily.parquet (133 regions)
- results/cams_immediate_daily.parquet (510 regions)
"""

import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np
import zipfile
import xarray as xr

print("="*70)
print("00g4: AGGREGATE CAMS POLLUTION (OPTIMIZED)")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# PATHS
# =============================================================================

BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR.parent.parent.parent / 'Input_data'  # sim_data/Input_data
OUTPUT_DIR = BASE_DIR.parent / 'results'  # Main results folder
CAMS_DIR = BASE_DIR.parent / 'temp_cams'  # phase0_data_prep/temp_cams
MAIN_RESULTS = BASE_DIR.parent / 'results'

STATES_FILE = INPUT_DIR / 'brazil_municipalities_2022.gpkg'
MAPPING_FILE = MAIN_RESULTS / 'municipality_to_all_regions_map.csv'

# =============================================================================
# LOAD DATA
# =============================================================================

try:
    import geopandas as gpd
    from scipy.spatial import cKDTree
    import shutil
    print("✓ Dependencies loaded")
except ImportError as e:
    print(f"✗ Missing dependency: {e}")
    exit(1)

print("\nLoading geographic data...")
gdf_munic = gpd.read_file(STATES_FILE)
df_mapping = pd.read_csv(MAPPING_FILE)

# Ensure CRS is WGS84 (EPSG:4326) for consistency with CAMS lon/lat grid
if gdf_munic.crs is None:
    print("  ⚠ No CRS defined, assuming EPSG:4326")
    gdf_munic = gdf_munic.set_crs('EPSG:4326')
elif gdf_munic.crs.to_epsg() != 4326:
    print(f"  Converting CRS from {gdf_munic.crs} to EPSG:4326")
    gdf_munic = gdf_munic.to_crs('EPSG:4326')
else:
    print(f"  ✓ CRS verified: {gdf_munic.crs}")

# Get municipality centroids using representative_point for multi-part geometries
# (representative_point is guaranteed to be within the polygon, unlike centroid)
gdf_munic['centroid'] = gdf_munic.geometry.representative_point()
muni_centroids = np.array([[p.x, p.y] for p in gdf_munic['centroid']])
muni_codes = gdf_munic['code_muni'].values

# Create lookup arrays (indexed by position in gdf_munic)
muni_to_intermediate = dict(zip(df_mapping['code_muni'], df_mapping['intermediate_code']))
muni_to_immediate = dict(zip(df_mapping['code_muni'], df_mapping['immediate_code']))

muni_intermediate_codes = np.array([muni_to_intermediate.get(int(m), -1) for m in muni_codes])
muni_immediate_codes = np.array([muni_to_immediate.get(int(m), -1) for m in muni_codes])

print(f"  {len(gdf_munic)} municipalities")
print(f"  {df_mapping['intermediate_code'].nunique()} intermediate regions")
print(f"  {df_mapping['immediate_code'].nunique()} immediate regions")

# =============================================================================
# VARIABLES
# =============================================================================

CAMS_VARS = {
    'pm2p5': 'pm25',
    'pm10': 'pm10',
    'gtco3': 'ozone',
}

# =============================================================================
# PROCESS ALL CAMS FILES
# =============================================================================

print("\n" + "-"*70)
print("Processing CAMS Files")
print("-"*70)

cams_files = sorted(CAMS_DIR.glob('*.zip'))
print(f"Found {len(cams_files)} CAMS files")

temp_extract_dir = CAMS_DIR / '_temp_extract'
temp_extract_dir.mkdir(exist_ok=True)

muni_to_grid_idx = None  # Will be computed once

all_results = []  # Collect all municipality-day results

for i, zip_path in enumerate(cams_files):
    print(f"\n[{i+1}/{len(cams_files)}] {zip_path.name}...", end=" ", flush=True)
    
    try:
        # Extract and load
        with zipfile.ZipFile(zip_path, 'r') as zf:
            zf.extractall(temp_extract_dir)
        
        nc_files = list(temp_extract_dir.glob('*.nc'))
        if not nc_files:
            print("no .nc file")
            continue
        
        ds = xr.open_dataset(nc_files[0], engine='netcdf4')
        
        # Get coordinates
        lats = ds['latitude'].values
        lons = ds['longitude'].values
        
        if 'valid_time' in ds.coords or 'valid_time' in ds.dims:
            times = pd.to_datetime(ds['valid_time'].values)
            time_dim = 'valid_time'
        else:
            times = pd.to_datetime(ds['time'].values)
            time_dim = 'time'
        
        # Build municipality-to-grid mapping once
        if muni_to_grid_idx is None:
            lon_grid, lat_grid = np.meshgrid(lons, lats)
            grid_points = np.column_stack([lon_grid.ravel(), lat_grid.ravel()])
            grid_tree = cKDTree(grid_points)
            _, muni_to_grid_idx = grid_tree.query(muni_centroids)
            print(f"(grid mapping built)", end=" ")
        
        # Extract data for all timesteps at once
        data_arrays = {}
        for cams_var, our_var in CAMS_VARS.items():
            if cams_var in ds.data_vars:
                arr = ds[cams_var].values
                
                # Defensive handling: check for extra dimensions (e.g., level)
                if arr.ndim == 4:  # e.g., (time, level, lat, lon)
                    print(f"    ⚠ {cams_var} has 4D shape {arr.shape}, squeezing...")
                    arr = arr.squeeze()  # Remove single-element dimensions
                elif arr.ndim != 3:
                    print(f"    ✗ {cams_var} unexpected shape {arr.shape}, skipping")
                    continue
                
                # Verify expected shape: (time, lat, lon)
                if arr.shape[0] != len(times):
                    print(f"    ✗ {cams_var} time dim mismatch: {arr.shape[0]} vs {len(times)}")
                    continue
                
                # Reshape to (time, n_grid_points)
                arr = arr.reshape(len(times), -1)
                
                # Unit conversions with verification from attributes
                orig_units = ds[cams_var].attrs.get('units', 'unknown')
                if cams_var in ['pm2p5', 'pm10']:
                    # kg/m³ to µg/m³
                    if 'kg' in orig_units.lower():
                        arr = arr * 1e9
                    else:
                        print(f"    ⚠ {cams_var} units '{orig_units}', applying 1e9 conversion anyway")
                        arr = arr * 1e9
                elif cams_var == 'gtco3':
                    # kg/m² to Dobson Units
                    if 'kg' in orig_units.lower():
                        arr = arr / 2.1415e-5
                    else:
                        print(f"    ⚠ {cams_var} units '{orig_units}', applying DU conversion anyway")
                        arr = arr / 2.1415e-5
                
                # Get values at municipality grid points: (time, n_municipalities)
                data_arrays[our_var] = arr[:, muni_to_grid_idx]
        
        ds.close()
        
        # Create dataframe with all timesteps
        n_times = len(times)
        n_munis = len(muni_codes)
        
        # Build flat arrays for DataFrame
        dates = np.repeat(times.date, n_munis)
        intermediate_codes = np.tile(muni_intermediate_codes, n_times)
        immediate_codes = np.tile(muni_immediate_codes, n_times)
        
        df_dict = {
            'date': dates,
            'intermediate_code': intermediate_codes,
            'immediate_code': immediate_codes,
        }
        
        for var, arr in data_arrays.items():
            df_dict[var] = arr.ravel()
        
        df_file = pd.DataFrame(df_dict)
        
        # Filter valid region codes
        df_file = df_file[(df_file['intermediate_code'] >= 0) & (df_file['immediate_code'] >= 0)]
        
        all_results.append(df_file)
        
        # Clean up temp files completely (including subdirectories)
        try:
            shutil.rmtree(temp_extract_dir)
            temp_extract_dir.mkdir(exist_ok=True)  # Recreate for next iteration
        except Exception as e:
            print(f"    ⚠ Cleanup warning: {e}")
        
        print(f"{len(times)} days")
        
    except Exception as e:
        print(f"error: {e}")
        import traceback
        traceback.print_exc()

# =============================================================================
# AGGREGATE TO REGIONS
# =============================================================================

print("\n" + "-"*70)
print("Aggregating to Regions")
print("-"*70)

# Combine all data
df_all = pd.concat(all_results, ignore_index=True)
print(f"Total municipality-days: {len(df_all):,}")

# Convert date
df_all['date'] = pd.to_datetime(df_all['date'])

# Aggregate to INTERMEDIATE level (simple mean - fast!)
print("\nAggregating to intermediate regions...")
agg_dict = {var: 'mean' for var in CAMS_VARS.values()}
df_intermediate = df_all.groupby(['intermediate_code', 'date'], as_index=False).agg(agg_dict)
df_intermediate = df_intermediate.sort_values(['intermediate_code', 'date']).reset_index(drop=True)

output_inter = OUTPUT_DIR / 'cams_intermediate_daily.parquet'
df_intermediate.to_parquet(output_inter, index=False)

print(f"✓ Saved: {output_inter.name}")
print(f"  Rows: {len(df_intermediate):,}")
print(f"  Regions: {df_intermediate['intermediate_code'].nunique()}")
print(f"  Date range: {df_intermediate['date'].min().date()} to {df_intermediate['date'].max().date()}")
for var in CAMS_VARS.values():
    print(f"  {var}: mean={df_intermediate[var].mean():.2f}, std={df_intermediate[var].std():.2f}")

# Aggregate to IMMEDIATE level
print("\nAggregating to immediate regions...")
df_immediate = df_all.groupby(['immediate_code', 'date'], as_index=False).agg(agg_dict)
df_immediate = df_immediate.sort_values(['immediate_code', 'date']).reset_index(drop=True)

output_immed = OUTPUT_DIR / 'cams_immediate_daily.parquet'
df_immediate.to_parquet(output_immed, index=False)

print(f"✓ Saved: {output_immed.name}")
print(f"  Rows: {len(df_immediate):,}")
print(f"  Regions: {df_immediate['immediate_code'].nunique()}")
print(f"  Date range: {df_immediate['date'].min().date()} to {df_immediate['date'].max().date()}")
for var in CAMS_VARS.values():
    print(f"  {var}: mean={df_immediate[var].mean():.2f}, std={df_immediate[var].std():.2f}")

# Clean up
import shutil
try:
    shutil.rmtree(temp_extract_dir)
except:
    pass

print("\n" + "="*70)
print("DONE!")
print("="*70)
print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
