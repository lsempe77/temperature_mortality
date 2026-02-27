"""
00h2_fast: FAST ERA5 aggregation to intermediate regions
=========================================================
Vectorized version: loads entire year at once, uses numpy bincount
for region aggregation instead of pandas day-by-day loop.

Only produces intermediate-level output (133 regions) for meta-regression.
"""
import sys, warnings
warnings.filterwarnings('ignore')
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd

print("="*70)
print("FAST ERA5 AGGREGATION -> INTERMEDIATE REGIONS")
print("="*70)
print(f"Started: {datetime.now():%Y-%m-%d %H:%M:%S}")
sys.stdout.flush()

BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR.parent.parent.parent / 'Input_data'
OUTPUT_DIR = BASE_DIR.parent / 'results'

# --- Load geographic mapping ---
print("\nLoading geographic data...")
sys.stdout.flush()
import geopandas as gpd
from scipy.spatial import cKDTree

gdf = gpd.read_file(INPUT_DIR / 'brazil_municipalities_2022.gpkg')
if gdf.crs and gdf.crs.to_epsg() != 4326:
    gdf = gdf.to_crs('EPSG:4326')

centroids = np.array([[p.x, p.y] for p in gdf.geometry.representative_point()])
muni_codes = gdf['code_muni'].values
tree = cKDTree(centroids)

mapping = pd.read_csv(OUTPUT_DIR / 'municipality_to_all_regions_map.csv')
muni_to_inter = dict(zip(mapping['code_muni'], mapping['intermediate_code']))
n_inter = mapping['intermediate_code'].nunique()
print(f"  {len(gdf)} municipalities -> {n_inter} intermediate regions")
sys.stdout.flush()

# --- Get ERA5 files ---
era5_files = sorted(INPUT_DIR.glob('era5_brazil_hourly_*.nc'))
print(f"  {len(era5_files)} ERA5 files found")
sys.stdout.flush()

import xarray as xr

grid_inter_codes = None   # computed once
all_results = []

for fi, era5_file in enumerate(era5_files):
    year = int(era5_file.stem.split('_')[-1])
    print(f"\n[{fi+1}/{len(era5_files)}] {year}...", end=" ")
    sys.stdout.flush()
    
    try:
        ds = xr.open_dataset(era5_file)
        
        lats = ds['latitude'].values
        lons = ds['longitude'].values
        
        # Find time dim
        for td in ('valid_time', 'time', 'times'):
            if td in ds.coords or td in ds.dims:
                time_dim = td
                break
        else:
            print("SKIP (no time dim)")
            ds.close()
            continue
        
        times = pd.DatetimeIndex(ds[time_dim].values)
        
        # Grid-to-region mapping (built once)
        if grid_inter_codes is None:
            lon_grid, lat_grid = np.meshgrid(lons, lats)
            grid_pts = np.column_stack([lon_grid.ravel(), lat_grid.ravel()])
            _, indices = tree.query(grid_pts)
            grid_muni = muni_codes[indices]
            grid_inter_codes = np.array([muni_to_inter.get(int(m), -1) for m in grid_muni])
            
            # Pre-compute region indices for bincount
            unique_codes = np.sort(np.unique(grid_inter_codes[grid_inter_codes >= 0]))
            code_to_idx = {c: i for i, c in enumerate(unique_codes)}
            grid_idx = np.array([code_to_idx.get(c, -1) for c in grid_inter_codes])
            valid_mask = grid_idx >= 0
            n_regions = len(unique_codes)
            print(f"[grid: {len(lats)}x{len(lons)}, {n_regions} regions]", end=" ")
            sys.stdout.flush()
        
        # Determine temperature variable
        for tvar in ('t2m', '2m_temperature'):
            if tvar in ds.data_vars:
                temp_var = tvar
                break
        else:
            temp_var = list(ds.data_vars)[0]
        
        # Check for dewpoint
        dew_var = None
        for dv in ('d2m', '2m_dewpoint_temperature'):
            if dv in ds.data_vars:
                dew_var = dv
                break
        
        # === VECTORIZED: Load entire year, compute daily stats ===
        # Load all temperature data at once: shape (n_times, n_lat, n_lon)
        temp_all = ds[temp_var].values  # load into memory
        n_times = temp_all.shape[0]
        n_grid = temp_all.shape[1] * temp_all.shape[2]
        temp_flat = temp_all.reshape(n_times, n_grid)  # (n_times, n_grid)
        
        # Convert K -> C if needed
        if temp_flat.mean() > 200:
            temp_flat = temp_flat - 273.15
        
        # Load dewpoint if available
        dew_flat = None
        if dew_var:
            dew_all = ds[dew_var].values
            dew_flat = dew_all.reshape(n_times, n_grid)
            if dew_flat.mean() > 200:
                dew_flat = dew_flat - 273.15
            del dew_all
        
        del temp_all
        ds.close()
        
        # Get dates for each timestep
        dates = times.date
        unique_dates = np.unique(dates)
        n_days = len(unique_dates)
        
        # Vectorized: compute daily stats per grid point, then aggregate to regions
        valid_idx = grid_idx[valid_mask]  # region indices for valid grid points
        counts = np.bincount(valid_idx, minlength=n_regions).astype(float)
        counts[counts == 0] = np.nan
        
        # Pre-allocate arrays for this year
        year_rmean = np.empty((n_days, n_regions))
        year_rmin  = np.empty((n_days, n_regions))
        year_rmax  = np.empty((n_days, n_regions))
        year_dates = []
        
        for di, day in enumerate(unique_dates):
            day_mask = dates == day
            day_temps = temp_flat[day_mask][:, valid_mask]  # (hours, valid_grid)
            
            gmean = np.nanmean(day_temps, axis=0)
            gmin  = np.nanmin(day_temps, axis=0)
            gmax  = np.nanmax(day_temps, axis=0)
            
            year_rmean[di] = np.bincount(valid_idx, weights=gmean, minlength=n_regions) / counts
            year_rmin[di]  = np.bincount(valid_idx, weights=gmin,  minlength=n_regions) / counts
            year_rmax[di]  = np.bincount(valid_idx, weights=gmax,  minlength=n_regions) / counts
            year_dates.append(str(day))
        
        # Build DataFrame for this year in one shot
        date_arr = np.repeat(year_dates, n_regions)
        code_arr = np.tile(unique_codes, n_days)
        
        year_df = pd.DataFrame({
            'date': date_arr,
            'region_code': code_arr.astype(int),
            'temp_mean': year_rmean.ravel(),
            'temp_min': year_rmin.ravel(),
            'temp_max': year_rmax.ravel()
        })
        all_results.append(year_df)
        
        del temp_flat, dew_flat, year_rmean, year_rmin, year_rmax
        print(f"OK ({n_days} days)")
        sys.stdout.flush()
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback; traceback.print_exc()
        sys.stdout.flush()
        continue

# --- Save ---
print(f"\nCombining {len(all_results)} yearly DataFrames...")
sys.stdout.flush()

df = pd.concat(all_results, ignore_index=True)
df = df.sort_values(['region_code', 'date']).reset_index(drop=True)

print(f"  Observations: {len(df):,}")
print(f"  Regions: {df['region_code'].nunique()}")
print(f"  Date range: {df['date'].min()} to {df['date'].max()}")
print(f"  Temperature: mean={df['temp_mean'].mean():.1f} C")
sys.stdout.flush()

outfile = OUTPUT_DIR / 'era5_intermediate_daily.parquet'
df.to_parquet(outfile, index=False)
print(f"  Saved: {outfile}")

# Also save a CSV sample
df.head(1000).to_csv(OUTPUT_DIR / 'era5_intermediate_daily_sample.csv', index=False)
print("  Saved CSV sample")

print(f"\nDONE! {datetime.now():%H:%M:%S}")
sys.stdout.flush()
