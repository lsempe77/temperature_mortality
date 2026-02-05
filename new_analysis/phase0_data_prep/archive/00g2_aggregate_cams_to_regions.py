"""
00g2: AGGREGATE CAMS POLLUTION DATA TO REGION-DAY LEVEL
=========================================================
Spatially aggregate gridded CAMS pollution data to Brazilian geographic regions.

Input: 
- temp_cams/*.zip (quarterly CAMS files, 0.75° grid)
- brazil_municipalities_2022.gpkg (municipality boundaries)
- results/municipality_region_mapping.csv (municipality to region mapping)

Output:
- results/cams_intermediate_daily.parquet (133 intermediate regions × ~5,500 days)
- results/cams_immediate_daily.parquet (510 immediate regions × ~5,500 days)

Method:
1. Extract and process each quarterly CAMS file
2. Assign each grid cell to nearest municipality centroid
3. Aggregate municipality-level to both intermediate (133) and immediate (510) regions
4. Simple mean across grid cells within each region

This creates the pollution data at the same levels as our mortality analysis.
"""

import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np
import zipfile
import tempfile
import xarray as xr

print("="*70)
print("00g2: AGGREGATE CAMS POLLUTION TO INTERMEDIATE AND IMMEDIATE REGIONS")
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
    from shapely.geometry import Point
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
centroids = np.array([[p.x, p.y] for p in gdf_munic['centroid']])
muni_codes = gdf_munic['code_muni'].values

# Build KD-tree for fast nearest-neighbor lookup
tree = cKDTree(centroids)
print(f"  Built KD-tree with {len(centroids)} municipality centroids")

# Load region mapping
df_mapping = pd.read_csv(MAPPING_FILE)

# Create mappings for both levels
muni_to_intermediate = dict(zip(df_mapping['code_muni'], df_mapping['intermediate_code']))
muni_to_immediate = dict(zip(df_mapping['code_muni'], df_mapping['immediate_code']))

print(f"  Loaded mapping: {len(df_mapping)} municipalities")
print(f"    → {df_mapping['intermediate_code'].nunique()} intermediate regions")
print(f"    → {df_mapping['immediate_code'].nunique()} immediate regions")

# =============================================================================
# PROCESS CAMS FILES
# =============================================================================

print("\n" + "-"*70)
print("Processing CAMS Quarterly Files")
print("-"*70)

# Get list of CAMS zip files
cams_files = sorted(CAMS_DIR.glob('*.zip'))
print(f"Found {len(cams_files)} CAMS files to process")

# Variables of interest - actual CAMS variable names
CAMS_VARS = {
    'pm2p5': 'pm25',       # PM2.5 (µg/m³)
    'pm10': 'pm10',        # PM10 (µg/m³)
    'gtco3': 'ozone',      # Total column ozone
}

all_data = []

# Create a persistent temp directory for extraction
temp_extract_dir = CAMS_DIR / '_temp_extract'
temp_extract_dir.mkdir(exist_ok=True)

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
            # Clean up
            for f in temp_extract_dir.glob('*'):
                f.unlink()
            continue
        
        # Get coordinates - handle valid_time vs time
        lats = ds['latitude'].values
        lons = ds['longitude'].values
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
        
        # Create grid points
        lon_grid, lat_grid = np.meshgrid(lons, lats)
        grid_points = np.column_stack([lon_grid.ravel(), lat_grid.ravel()])
        
        # Assign each grid point to nearest municipality
        _, indices = tree.query(grid_points)
        grid_muni_codes = muni_codes[indices]
        grid_intermediate = np.array([muni_to_intermediate.get(int(m), np.nan) for m in grid_muni_codes])
        grid_immediate = np.array([muni_to_immediate.get(int(m), np.nan) for m in grid_muni_codes])
        
        # Process each day
        for t_idx, t in enumerate(times):
            day_data = {'date': t.date(), 'time': t}
            
            # Extract values for each variable (use dynamic time dimension)
            for cams_var, our_var in CAMS_VARS.items():
                if cams_var in ds.data_vars:
                    vals = ds[cams_var].isel({time_dim: t_idx}).values.ravel()
                    day_data[our_var] = vals
                else:
                    day_data[our_var] = np.full(len(grid_intermediate), np.nan)
            
            # Create dataframe for this timestep with both region levels
            df_day = pd.DataFrame({
                'date': day_data['date'],
                'intermediate_code': grid_intermediate,
                'immediate_code': grid_immediate,
                **{k: day_data[k] for k in CAMS_VARS.values()}
            })
            
            # Aggregate to INTERMEDIATE region level
            df_inter = df_day.dropna(subset=['intermediate_code']).copy()
            df_inter['intermediate_code'] = df_inter['intermediate_code'].astype(int)
            df_inter_agg = df_inter.groupby(['date', 'intermediate_code']).agg({
                'pm25': 'mean', 'pm10': 'mean', 'ozone': 'mean'
            }).reset_index()
            df_inter_agg = df_inter_agg.rename(columns={'intermediate_code': 'region_code'})
            
            # Aggregate to IMMEDIATE region level
            df_imm = df_day.dropna(subset=['immediate_code']).copy()
            df_imm['immediate_code'] = df_imm['immediate_code'].astype(int)
            df_imm_agg = df_imm.groupby(['date', 'immediate_code']).agg({
                'pm25': 'mean', 'pm10': 'mean', 'ozone': 'mean'
            }).reset_index()
            df_imm_agg = df_imm_agg.rename(columns={'immediate_code': 'region_code'})
            
            all_data.append(('intermediate', df_inter_agg))
            all_data.append(('immediate', df_imm_agg))
        
        ds.close()
        del ds  # Ensure file handle is released on Windows
        
        # Clean up extracted files
        for f in temp_extract_dir.glob('*'):
            f.unlink()
        
        print(f"  ✓ Processed {len(times)} timesteps")
        
    except Exception as e:
        print(f"  ✗ Error processing {zip_path.name}: {e}")
        # Clean up on error
        for f in temp_extract_dir.glob('*'):
            try:
                f.unlink()
            except:
                pass
        continue

# Clean up temp directory
try:
    temp_extract_dir.rmdir()
except:
    pass

# =============================================================================
# COMBINE AND SAVE
# =============================================================================

print("\n" + "-"*70)
print("Combining Results")
print("-"*70)

if all_data:
    # Separate by region type
    intermediate_dfs = [df for level, df in all_data if level == 'intermediate']
    immediate_dfs = [df for level, df in all_data if level == 'immediate']
    
    # Process INTERMEDIATE regions
    print("\n--- INTERMEDIATE REGIONS (133) ---")
    df_inter_all = pd.concat(intermediate_dfs, ignore_index=True)
    
    # Aggregate to daily (if hourly data)
    df_inter_daily = df_inter_all.groupby(['date', 'region_code']).agg({
        'pm25': 'mean', 'pm10': 'mean', 'ozone': 'mean'
    }).reset_index()
    
    print(f"Daily observations: {len(df_inter_daily):,}")
    print(f"Regions: {df_inter_daily['region_code'].nunique()}")
    print(f"Date range: {df_inter_daily['date'].min()} to {df_inter_daily['date'].max()}")
    
    output_inter = OUTPUT_DIR / 'cams_intermediate_daily.parquet'
    df_inter_daily.to_parquet(output_inter, index=False)
    print(f"✓ Saved: {output_inter}")
    
    # Process IMMEDIATE regions
    print("\n--- IMMEDIATE REGIONS (510) ---")
    df_imm_all = pd.concat(immediate_dfs, ignore_index=True)
    
    # Aggregate to daily (if hourly data)
    df_imm_daily = df_imm_all.groupby(['date', 'region_code']).agg({
        'pm25': 'mean', 'pm10': 'mean', 'ozone': 'mean'
    }).reset_index()
    
    print(f"Daily observations: {len(df_imm_daily):,}")
    print(f"Regions: {df_imm_daily['region_code'].nunique()}")
    print(f"Date range: {df_imm_daily['date'].min()} to {df_imm_daily['date'].max()}")
    
    output_imm = OUTPUT_DIR / 'cams_immediate_daily.parquet'
    df_imm_daily.to_parquet(output_imm, index=False)
    print(f"✓ Saved: {output_imm}")
    
    # Summary stats
    print("\nPollution Summary - Intermediate (mean ± std):")
    for var in ['pm25', 'pm10', 'ozone']:
        if var in df_inter_daily.columns:
            mean_val = df_inter_daily[var].mean()
            std_val = df_inter_daily[var].std()
            print(f"  {var}: {mean_val:.2f} ± {std_val:.2f} µg/m³")
else:
    print("✗ No data processed!")

print(f"\n{'='*70}")
print("DONE!")
print("="*70)
