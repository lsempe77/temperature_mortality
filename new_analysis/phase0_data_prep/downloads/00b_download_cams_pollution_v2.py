"""
00b: DOWNLOAD CAMS AIR POLLUTION DATA (CORRECTED)
==================================================
Download PM2.5 and O3 data from Copernicus Atmosphere Data Store (ADS).

IMPORTANT: This uses the ADS API, NOT the CDS API!
The .cdsapirc file must contain:
    url: https://ads.atmosphere.copernicus.eu/api
    key: YOUR_PERSONAL_ACCESS_TOKEN

Data source: CAMS Global Reanalysis (EAC4)
- Dataset ID: cams-global-reanalysis-eac4
- Resolution: 0.75° x 0.75° (~83 km)
- Variables: PM2.5 (kg/m³), PM10 (kg/m³), Total Column Ozone (kg/m²)
- Period: 2003-2024

Why control for air pollution:
- PM2.5 and O3 are correlated with temperature (photochemistry, stagnation)
- Both independently affect respiratory and cardiovascular mortality
- Not controlling for pollution may bias heat effect estimates upward

Output: results/cams_pollution_brazil_daily.parquet
Documentation: See CAMS_DATA_DOWNLOAD_GUIDE.md
"""

import os
import sys
from datetime import datetime
from pathlib import Path
import warnings

print("="*70)
print("00b: DOWNLOAD CAMS AIR POLLUTION DATA (CORRECTED)")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# CHECK DEPENDENCIES
# =============================================================================

try:
    import cdsapi
    print("✓ cdsapi installed")
except ImportError:
    print("✗ cdsapi not installed")
    print("  Install with: pip install 'cdsapi>=0.7.7'")
    sys.exit(1)

try:
    import xarray as xr
    print("✓ xarray installed")
except ImportError:
    print("✗ xarray not installed")
    print("  Install with: pip install xarray netcdf4")
    sys.exit(1)

try:
    import pandas as pd
    import numpy as np
    print("✓ pandas/numpy installed")
except ImportError:
    print("✗ pandas/numpy not installed")
    sys.exit(1)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Brazil bounding box with buffer
# Format: [North, West, South, East]
BRAZIL_AREA = [5, -75, -35, -30]

# Years to download (matching our mortality data)
# Extended to 2010-2024 for long-term analysis
YEARS = list(range(2010, 2025))

# Variables to download (single-level, fast-access)
# See: https://confluence.ecmwf.int/x/OIX4B
VARIABLES = [
    'particulate_matter_2.5um',    # PM2.5 in kg/m³
    'particulate_matter_10um',     # PM10 in kg/m³  
    'total_column_ozone',          # Total O3 in kg/m²
]

# Analysis times (3-hourly)
TIMES = ['00:00', '06:00', '12:00', '18:00']  # 4 times per day

# Output directories
BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'results'
TEMP_DIR = BASE_DIR / 'temp_cams'
OUTPUT_DIR.mkdir(exist_ok=True)
TEMP_DIR.mkdir(exist_ok=True)

# Dataset name
DATASET_ID = 'cams-global-reanalysis-eac4'

# =============================================================================
# API CONFIGURATION CHECK
# =============================================================================

def check_api_config():
    """
    Verify the .cdsapirc file exists and get API key.
    
    Note: Copernicus has unified CDS and ADS authentication systems.
    A CDS key now works for ADS too!
    
    Returns tuple (success: bool, key: str or None)
    """
    
    cdsapirc = Path.home() / '.cdsapirc'
    
    print("\n" + "-"*70)
    print("Checking API Configuration")
    print("-"*70)
    
    if not cdsapirc.exists():
        print(f"✗ Config file not found: {cdsapirc}")
        print_setup_instructions()
        return False, None
    
    # Read and parse the config
    with open(cdsapirc, 'r') as f:
        content = f.read()
    
    print(f"✓ Found config at: {cdsapirc}")
    
    # Extract the key from the config
    api_key = None
    for line in content.split('\n'):
        if line.strip().startswith('key:'):
            api_key = line.split(':', 1)[1].strip()
            break
    
    if not api_key:
        print("✗ Could not find API key in config")
        return False, None
    
    # Note about CDS vs ADS
    if 'cds.climate.copernicus.eu' in content:
        print("ℹ Config uses CDS URL, but Copernicus has unified auth systems")
        print("  Your CDS key will work for ADS too!")
    
    # Try to create client with ADS URL
    try:
        # Force ADS URL regardless of what's in config
        client = cdsapi.Client(
            url='https://ads.atmosphere.copernicus.eu/api',
            key=api_key
        )
        print(f"✓ API client created successfully for ADS")
        print(f"  URL: {client.url}")
        return True, api_key
    except Exception as e:
        print(f"✗ Error creating API client: {e}")
        return False, None


def print_setup_instructions():
    """Print detailed setup instructions."""
    
    print("""
╔══════════════════════════════════════════════════════════════════════╗
║  ADS API CONFIGURATION REQUIRED                                      ║
╚══════════════════════════════════════════════════════════════════════╝

CAMS data is on the Atmosphere Data Store (ADS), NOT the Climate Data Store (CDS).

STEP 1: Register at ADS
   https://ads.atmosphere.copernicus.eu/

STEP 2: Get your Personal Access Token
   - Log in to ADS
   - Go to your profile
   - Copy the API key

STEP 3: Create or update ~/.cdsapirc
   
   On Windows, create: C:\\Users\\<YourName>\\.cdsapirc
   
   Contents (two lines):
   ┌────────────────────────────────────────────────────────────────────┐
   │ url: https://ads.atmosphere.copernicus.eu/api                     │
   │ key: YOUR_PERSONAL_ACCESS_TOKEN_HERE                              │
   └────────────────────────────────────────────────────────────────────┘

STEP 4: Accept license terms
   - Go to: https://ads.atmosphere.copernicus.eu/datasets/cams-global-reanalysis-eac4
   - Click "Download" tab
   - Accept the CC-BY licence at the bottom

STEP 5: Re-run this script
   python 00b_download_cams_pollution.py

For detailed documentation, see:
   docs/CAMS_DATA_DOWNLOAD_GUIDE.md
""")


# =============================================================================
# DOWNLOAD FUNCTIONS
# =============================================================================

def download_cams_quarter(year: int, quarter: int, output_file: Path) -> bool:
    """
    Download CAMS EAC4 data for one quarter.
    
    Splitting by quarter avoids timeout on large requests.
    
    Parameters:
    -----------
    year : int
        Year to download
    quarter : int
        Quarter (1-4)
    output_file : Path
        Output file path (should end in .zip for netcdf_zip format)
    
    Returns:
    --------
    bool : True if successful
    """
    
    if output_file.exists():
        print(f"  ✓ Already exists: {output_file.name}")
        return True
    
    # Define quarter date ranges
    quarter_dates = {
        1: (f'{year}-01-01', f'{year}-03-31'),
        2: (f'{year}-04-01', f'{year}-06-30'),
        3: (f'{year}-07-01', f'{year}-09-30'),
        4: (f'{year}-10-01', f'{year}-12-31'),
    }
    
    start_date, end_date = quarter_dates[quarter]
    
    print(f"  Downloading Q{quarter} ({start_date} to {end_date})...")
    
    try:
        # Get API key from config file
        cdsapirc = Path.home() / '.cdsapirc'
        api_key = None
        if cdsapirc.exists():
            with open(cdsapirc, 'r') as f:
                for line in f:
                    if line.strip().startswith('key:'):
                        api_key = line.split(':', 1)[1].strip()
                        break
        
        # Force ADS URL regardless of config file setting
        client = cdsapi.Client(
            url='https://ads.atmosphere.copernicus.eu/api',
            key=api_key
        )
        
        request = {
            'variable': VARIABLES,
            'date': f'{start_date}/{end_date}',
            'time': TIMES,
            'area': BRAZIL_AREA,  # [N, W, S, E]
            'data_format': 'netcdf_zip',  # Zipped NetCDF
        }
        
        client.retrieve(DATASET_ID, request, str(output_file))
        
        print(f"  ✓ Downloaded: {output_file.name}")
        return True
        
    except Exception as e:
        error_msg = str(e)
        
        if '404' in error_msg:
            print(f"  ✗ 404 Error - Dataset not found")
            print(f"    Check that you accepted the license at:")
            print(f"    https://ads.atmosphere.copernicus.eu/datasets/{DATASET_ID}")
        elif 'terms' in error_msg.lower() or 'licence' in error_msg.lower():
            print(f"  ✗ License not accepted")
            print(f"    Accept at: https://ads.atmosphere.copernicus.eu/datasets/{DATASET_ID}")
        else:
            print(f"  ✗ Error: {error_msg[:200]}")
        
        return False


def download_year(year: int) -> list:
    """
    Download all quarters for a year.
    
    Returns list of downloaded files.
    """
    
    print(f"\nYear {year}:")
    
    downloaded = []
    
    for quarter in [1, 2, 3, 4]:
        output_file = TEMP_DIR / f'cams_brazil_{year}_Q{quarter}.zip'
        
        if download_cams_quarter(year, quarter, output_file):
            downloaded.append(output_file)
    
    return downloaded


# =============================================================================
# DATA PROCESSING
# =============================================================================

def extract_and_process_files(zip_files: list) -> pd.DataFrame:
    """
    Extract zip files and process to daily means.
    
    Parameters:
    -----------
    zip_files : list
        List of downloaded zip file paths
    
    Returns:
    --------
    pd.DataFrame with daily pollution data
    """
    
    import zipfile
    
    print("\n" + "-"*70)
    print("Processing Downloaded Data")
    print("-"*70)
    
    all_data = []
    
    for zf in zip_files:
        if not zf.exists():
            continue
            
        print(f"\nProcessing: {zf.name}")
        
        # Extract zip
        extract_dir = TEMP_DIR / zf.stem
        extract_dir.mkdir(exist_ok=True)
        
        with zipfile.ZipFile(zf, 'r') as z:
            z.extractall(extract_dir)
        
        # Find NetCDF files
        nc_files = list(extract_dir.glob('*.nc'))
        
        if not nc_files:
            print(f"  ⚠ No NetCDF files found in {zf.name}")
            continue
        
        for nc_file in nc_files:
            print(f"  Loading: {nc_file.name}")
            
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    ds = xr.open_dataset(nc_file)
                
                # Print variable info
                print(f"    Variables: {list(ds.data_vars)}")
                
                # Check which time dimension exists
                if 'valid_time' in ds.dims:
                    time_dim = 'valid_time'
                elif 'time' in ds.dims:
                    time_dim = 'time'
                else:
                    print(f"    ⚠ No time dimension found! Dims: {ds.dims}")
                    continue
                
                print(f"    Time dimension: {time_dim} ({ds.dims[time_dim]} timesteps)")
                
                # Compute daily means using the correct dimension
                ds_daily = ds.resample({time_dim: '1D'}).mean()
                
                # Convert to dataframe
                df = ds_daily.to_dataframe().reset_index()
                
                all_data.append(df)
                
                ds.close()
                
            except Exception as e:
                print(f"    Error: {e}")
                import traceback
                traceback.print_exc()
    
    if not all_data:
        print("✗ No data processed!")
        return None
    
    # Combine all data
    df = pd.concat(all_data, ignore_index=True)
    
    # Standardize column names
    # Note: CAMS EAC4 uses 'gtco3' for total column ozone, not 'tco3'
    rename_map = {
        'pm2p5': 'pm25_kgm3',
        'pm10': 'pm10_kgm3',
        'gtco3': 'o3_column_kgm2',  # Fixed: was 'tco3'
        'tco3': 'o3_column_kgm2',   # Keep for backwards compatibility
        'valid_time': 'date',
        'time': 'date',
    }
    
    # Apply renames for columns that exist
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})
    
    # Ensure we have a date column
    if 'date' not in df.columns and 'time' in df.columns:
        df['date'] = pd.to_datetime(df['time']).dt.date
    elif 'date' in df.columns:
        df['date'] = pd.to_datetime(df['date']).dt.date
    
    # Convert units
    # PM2.5 and PM10: kg/m³ → µg/m³ (multiply by 1e9)
    for col in ['pm25_kgm3', 'pm10_kgm3']:
        if col in df.columns:
            new_col = col.replace('_kgm3', '_ugm3')
            df[new_col] = df[col] * 1e9
    
    # Sort and remove duplicates
    df = df.sort_values(['date', 'latitude', 'longitude'])
    df = df.drop_duplicates(subset=['date', 'latitude', 'longitude'])
    
    print(f"\nProcessed data shape: {df.shape}")
    
    return df


# =============================================================================
# CREATE PLACEHOLDER IF DOWNLOAD FAILS
# =============================================================================

def create_placeholder():
    """
    Create a placeholder file with instructions if download fails.
    This allows the pipeline to continue while the user fixes API config.
    """
    
    placeholder = pd.DataFrame({
        'message': [
            'CAMS pollution data could not be downloaded.',
            'API configuration needed - see CAMS_DATA_DOWNLOAD_GUIDE.md',
            '',
            'As an alternative for pilot analysis:',
            '1. Use fixed effects to absorb regional pollution patterns',
            '2. Download data manually from ADS web interface',
            '3. Use published PM2.5 estimates (e.g., van Donkelaar et al.)',
        ]
    })
    
    output = OUTPUT_DIR / 'cams_pollution_PLACEHOLDER.csv'
    placeholder.to_csv(output, index=False)
    print(f"\n✓ Created placeholder: {output}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Main execution function."""
    
    # Step 1: Check API configuration
    success, api_key = check_api_config()
    if not success:
        create_placeholder()
        return
    
    # Step 2: Download data by year/quarter
    print("\n" + "-"*70)
    print("Downloading CAMS EAC4 Data")
    print("-"*70)
    
    all_files = []
    
    for year in YEARS:
        files = download_year(year)
        all_files.extend(files)
    
    if not all_files:
        print("\n✗ No files downloaded successfully")
        create_placeholder()
        return
    
    print(f"\n✓ Downloaded {len(all_files)} files")
    
    # Step 3: Process data
    df = extract_and_process_files(all_files)
    
    if df is None:
        create_placeholder()
        return
    
    # Step 4: Save results
    print("\n" + "-"*70)
    print("Saving Results")
    print("-"*70)
    
    output_file = OUTPUT_DIR / 'cams_pollution_brazil_daily.parquet'
    df.to_parquet(output_file, index=False)
    print(f"✓ Saved: {output_file}")
    
    # Also save as NetCDF for spatial analysis
    output_nc = OUTPUT_DIR / 'cams_pollution_brazil_daily.nc'
    
    try:
        df_for_nc = df.copy()
        # Convert date back to datetime for xarray
        df_for_nc['date'] = pd.to_datetime(df_for_nc['date'])
        
        # Set multi-index and convert to xarray
        df_indexed = df_for_nc.set_index(['date', 'latitude', 'longitude'])
        ds = df_indexed.to_xarray()
        ds.to_netcdf(output_nc)
        print(f"✓ Saved: {output_nc}")
    except Exception as e:
        print(f"⚠ Could not save NetCDF: {e}")
    
    # Step 5: Print summary
    print("\n" + "-"*70)
    print("Summary")
    print("-"*70)
    
    date_col = 'date' if 'date' in df.columns else 'time'
    print(f"Date range: {df[date_col].min()} to {df[date_col].max()}")
    print(f"Grid points: {df[['latitude', 'longitude']].drop_duplicates().shape[0]}")
    print(f"Days: {df[date_col].nunique()}")
    
    for col in df.columns:
        if 'pm25' in col or 'pm10' in col or 'o3' in col:
            if df[col].notna().any():
                print(f"{col}: {df[col].min():.3f} to {df[col].max():.3f}")
    
    # Cleanup temp files (optional)
    # import shutil
    # shutil.rmtree(TEMP_DIR)
    
    print(f"\n{'='*70}")
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*70)


if __name__ == '__main__':
    main()
