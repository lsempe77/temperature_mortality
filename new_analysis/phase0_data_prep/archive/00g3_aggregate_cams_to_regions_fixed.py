"""
00g3: AGGREGATE CAMS POLLUTION DATA TO REGION-DAY LEVEL (FIXED ALGORITHM)
==========================================================================
Spatially aggregate gridded CAMS pollution data to Brazilian geographic regions.

FIX: Uses municipality centroid → nearest grid point assignment
(instead of grid point → nearest municipality centroid)
This guarantees 100% coverage of all regions.

Input: 
- temp_cams/*.zip (quarterly CAMS files, 0.75° grid)
- brazil_municipalities_2022.gpkg (municipality boundaries)
- results/municipality_to_all_regions_map.csv (municipality to region mapping)

Output:
- results/cams_intermediate_daily.parquet (133 intermediate regions × ~5,500 days)
- results/cams_immediate_daily.parquet (510 immediate regions × ~5,500 days)

Method:
1. Extract and process each quarterly CAMS file
2. Assign each MUNICIPALITY to its nearest grid cell (reversed!)
3. Aggregate municipality-level to both intermediate (133) and immediate (510) regions
4. Area-weighted mean across municipalities within each region

Units:
- PM2.5: converted from kg/m³ to µg/m³ (×1e9)
- PM10: converted from kg/m³ to µg/m³ (×1e9)
- Ozone: converted from kg/m² to Dobson Units (÷2.1415e-5)
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
print("00g3: AGGREGATE CAMS POLLUTION (FIXED ALGORITHM)")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# PATHS
# =============================================================================

BASE_DIR = Path(__file__).parent
INPUT_DIR = Path(__file__).parent.parent.parent / 'Input_data'
OUTPUT_DIR = BASE_DIR / 'results'
CAMS_DIR = BASE_DIR / 'temp_cams'
MAIN_RESULTS = BASE_DIR.parent / 'results'

STATES_FILE = INPUT_DIR / 'brazil_municipalities_2022.gpkg'
MAPPING_FILE = MAIN_RESULTS / 'municipality_to_all_regions_map.csv'

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

# Get municipality centroids
gdf_munic['centroid'] = gdf_munic.geometry.centroid
muni_centroids = np.array([[p.x, p.y] for p in gdf_munic['centroid']])
muni_codes = gdf_munic['code_muni'].values

# Get municipality areas for weighted aggregation
gdf_proj = gdf_munic.to_crs('EPSG:5880')  # Brazil Polyconic
gdf_munic['area_km2'] = gdf_proj.geometry.area / 1e6
muni_areas = gdf_munic['area_km2'].values

# Load region mapping
df_mapping = pd.read_csv(MAPPING_FILE)

# Create mappings for both levels (code_muni to region codes)
muni_to_intermediate = dict(zip(df_mapping['code_muni'], df_mapping['intermediate_code']))
muni_to_immediate = dict(zip(df_mapping['code_muni'], df_mapping['immediate_code']))

print(f"  Loaded mapping: {len(df_mapping)} municipalities")
print(f"    → {df_mapping['intermediate_code'].nunique()} intermediate regions")
print(f"    → {df_mapping['immediate_code'].nunique()} immediate regions")

# =============================================================================
# VARIABLES OF INTEREST
# =============================================================================

# Actual CAMS variable names
CAMS_VARS = {
    'pm2p5': 'pm25',       # PM2.5
    'pm10': 'pm10',        # PM10
    'gtco3': 'ozone',      # Total column ozone
}

# =============================================================================
# PROCESS CAMS FILES
# =============================================================================

print("\n" + "-"*70)
print("Processing CAMS Quarterly Files")
print("-"*70)

# Get list of CAMS zip files
cams_files = sorted(CAMS_DIR.glob('*.zip'))
print(f"Found {len(cams_files)} CAMS files to process")

# Create a persistent temp directory for extraction
temp_extract_dir = CAMS_DIR / '_temp_extract'
temp_extract_dir.mkdir(exist_ok=True)

# Build municipality-to-grid mapping on first file
muni_to_grid_idx = None
grid_points = None

all_intermediate_data = []
all_immediate_data = []

for i, zip_path in enumerate(cams_files):
    print(f"\n[{i+1}/{len(cams_files)}] Processing {zip_path.name}...")
    
    try:
        # Extract zip to temp directory  
        with zipfile.ZipFile(zip_path, 'r') as zf:
            zf.extractall(temp_extract_dir)
        
        # Find the .nc file
        nc_files = list(temp_extract_dir.glob('*.nc'))
        if not nc_files:
            print(f"  ⚠ No .nc file found in {zip_path.name}")
            continue
        
        nc_file = nc_files[0]
        
        # Load with xarray
        ds = xr.open_dataset(nc_file, engine='netcdf4')
        
        # Get available variables
        available_vars = [v for v in CAMS_VARS.keys() if v in ds.data_vars]
        if not available_vars:
            print(f"  ⚠ No pollution variables found in {zip_path.name}")
            print(f"    Available: {list(ds.data_vars)}")
            ds.close()
            for f in temp_extract_dir.glob('*'):
                f.unlink()
            continue
        
        # Get coordinates
        lats = ds['latitude'].values
        lons = ds['longitude'].values
        
        # Handle valid_time vs time
        if 'valid_time' in ds.coords or 'valid_time' in ds.dims:
            times = pd.to_datetime(ds['valid_time'].values)
            time_dim = 'valid_time'
        elif 'time' in ds.coords or 'time' in ds.dims:
            times = pd.to_datetime(ds['time'].values)
            time_dim = 'time'
        else:
            print(f"  ⚠ No time coordinate found")
            ds.close()
            for f in temp_extract_dir.glob('*'):
                f.unlink()
            continue
        
        print(f"  Grid: {len(lats)}×{len(lons)} = {len(lats)*len(lons)} points")
        print(f"  Time steps: {len(times)}")
        print(f"  Variables: {available_vars}")
        
        # Build municipality-to-grid mapping on first file
        # FIX: Assign each MUNICIPALITY to its nearest GRID POINT
        if muni_to_grid_idx is None:
            print("  Building municipality-to-grid mapping (FIXED algorithm)...")
            lon_grid, lat_grid = np.meshgrid(lons, lats)
            grid_points = np.column_stack([lon_grid.ravel(), lat_grid.ravel()])
            
            # Build KD-tree of GRID POINTS (not municipality centroids!)
            grid_tree = cKDTree(grid_points)
            
            # For each municipality, find its nearest grid point
            _, muni_to_grid_idx = grid_tree.query(muni_centroids)
            
            # Create lookup: municipality index -> region codes
            muni_intermediate_codes = np.array([muni_to_intermediate.get(int(m), -1) for m in muni_codes])
            muni_immediate_codes = np.array([muni_to_immediate.get(int(m), -1) for m in muni_codes])
            
            # Verify coverage
            unique_intermediate = set(muni_intermediate_codes[muni_intermediate_codes >= 0])
            unique_immediate = set(muni_immediate_codes[muni_immediate_codes >= 0])
            print(f"  Coverage: {len(unique_intermediate)}/133 intermediate, {len(unique_immediate)}/510 immediate")
        
        # Process each timestep
        for t_idx, t in enumerate(times):
            # Get grid values for this timestep
            grid_values = {}
            for cams_var, our_var in CAMS_VARS.items():
                if cams_var in ds.data_vars:
                    vals = ds[cams_var].isel({time_dim: t_idx}).values.ravel()
                    
                    # Unit conversions
                    if cams_var in ['pm2p5', 'pm10']:
                        # kg/m³ to µg/m³
                        vals = vals * 1e9
                    elif cams_var == 'gtco3':
                        # kg/m² to Dobson Units
                        vals = vals / 2.1415e-5
                    
                    grid_values[our_var] = vals
                else:
                    grid_values[our_var] = np.full(len(grid_points), np.nan)
            
            # For each municipality, get its grid cell value
            muni_values = {var: grid_values[var][muni_to_grid_idx] for var in CAMS_VARS.values()}
            
            # Create municipality-level dataframe
            df_muni = pd.DataFrame({
                'date': t.date(),
                'intermediate_code': muni_intermediate_codes,
                'immediate_code': muni_immediate_codes,
                'area_km2': muni_areas,
                **muni_values
            })
            
            # Filter out invalid region codes
            df_muni = df_muni[(df_muni['intermediate_code'] >= 0) & (df_muni['immediate_code'] >= 0)]
            
            # Aggregate to INTERMEDIATE region level (area-weighted mean)
            def weighted_mean(group):
                result = {'date': group['date'].iloc[0]}
                total_area = group['area_km2'].sum()
                for var in CAMS_VARS.values():
                    valid = ~np.isnan(group[var])
                    if valid.any():
                        weights = group.loc[valid, 'area_km2'] / group.loc[valid, 'area_km2'].sum()
                        result[var] = (group.loc[valid, var] * weights).sum()
                    else:
                        result[var] = np.nan
                return pd.Series(result)
            
            df_intermediate = df_muni.groupby('intermediate_code').apply(weighted_mean, include_groups=False).reset_index()
            all_intermediate_data.append(df_intermediate)
            
            # Aggregate to IMMEDIATE region level (area-weighted mean)
            df_immediate = df_muni.groupby('immediate_code').apply(weighted_mean, include_groups=False).reset_index()
            all_immediate_data.append(df_immediate)
        
        ds.close()
        
        # Clean up temp files
        for f in temp_extract_dir.glob('*'):
            try:
                f.unlink()
            except:
                pass
        
        print(f"  ✓ Processed {len(times)} timesteps")
        
    except Exception as e:
        print(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()

# =============================================================================
# COMBINE AND SAVE
# =============================================================================

print("\n" + "-"*70)
print("Saving Results")
print("-"*70)

# Combine intermediate data
if all_intermediate_data:
    df_inter_all = pd.concat(all_intermediate_data, ignore_index=True)
    df_inter_all['date'] = pd.to_datetime(df_inter_all['date'])
    df_inter_all = df_inter_all.sort_values(['intermediate_code', 'date']).reset_index(drop=True)
    
    output_inter = OUTPUT_DIR / 'cams_intermediate_daily.parquet'
    df_inter_all.to_parquet(output_inter, index=False)
    print(f"✓ Saved intermediate: {output_inter.name}")
    print(f"  Rows: {len(df_inter_all):,}")
    print(f"  Regions: {df_inter_all['intermediate_code'].nunique()}")
    print(f"  Date range: {df_inter_all['date'].min()} to {df_inter_all['date'].max()}")
    
    # Summary stats
    print("\n  Variable statistics (intermediate):")
    for var in CAMS_VARS.values():
        mean_val = df_inter_all[var].mean()
        std_val = df_inter_all[var].std()
        print(f"    {var}: mean={mean_val:.2f}, std={std_val:.2f}")

# Combine immediate data
if all_immediate_data:
    df_immed_all = pd.concat(all_immediate_data, ignore_index=True)
    df_immed_all['date'] = pd.to_datetime(df_immed_all['date'])
    df_immed_all = df_immed_all.sort_values(['immediate_code', 'date']).reset_index(drop=True)
    
    output_immed = OUTPUT_DIR / 'cams_immediate_daily.parquet'
    df_immed_all.to_parquet(output_immed, index=False)
    print(f"\n✓ Saved immediate: {output_immed.name}")
    print(f"  Rows: {len(df_immed_all):,}")
    print(f"  Regions: {df_immed_all['immediate_code'].nunique()}")
    print(f"  Date range: {df_immed_all['date'].min()} to {df_immed_all['date'].max()}")
    
    # Summary stats
    print("\n  Variable statistics (immediate):")
    for var in CAMS_VARS.values():
        mean_val = df_immed_all[var].mean()
        std_val = df_immed_all[var].std()
        print(f"    {var}: mean={mean_val:.2f}, std={std_val:.2f}")

# Clean up temp directory
import shutil
try:
    shutil.rmtree(temp_extract_dir)
except:
    pass

print("\n" + "="*70)
print("CAMS AGGREGATION COMPLETE (FIXED)")
print("="*70)
print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
