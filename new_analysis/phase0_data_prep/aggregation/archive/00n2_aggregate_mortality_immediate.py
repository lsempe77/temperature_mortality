"""
00n2: AGGREGATE MORTALITY TO IMMEDIATE REGIONS (510)
=====================================================
Re-aggregate mortality data to immediate regions (510 units) for finer
spatial resolution analysis.

This complements the intermediate region (133) aggregation for sensitivity analysis.

Input: 
- Input_data/DO*.csv (yearly mortality files, 2010-2024)
- results/municipality_to_all_regions_map.csv (municipality to region mapping)

Output:
- results/mortality_immediate_daily.parquet (510 regions × ~5,500 days)
- results/mortality_immediate_daily_elderly.parquet (elderly 60+ only)
"""

import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np

print("="*70)
print("00n2: AGGREGATE MORTALITY TO IMMEDIATE REGIONS (510)")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# PATHS
# =============================================================================

BASE_DIR = Path(__file__).parent
INPUT_DIR = Path(__file__).parent.parent.parent / 'Input_data'
OUTPUT_DIR = BASE_DIR / 'results'
MAIN_RESULTS = BASE_DIR.parent / 'results'

MAPPING_FILE = MAIN_RESULTS / 'municipality_to_all_regions_map.csv'
OUTPUT_FILE = OUTPUT_DIR / 'mortality_immediate_daily.parquet'
OUTPUT_ELDERLY = OUTPUT_DIR / 'mortality_immediate_daily_elderly.parquet'

# =============================================================================
# LOAD REGION MAPPING
# =============================================================================

print("\n" + "-"*70)
print("Loading Region Mapping")
print("-"*70)

df_mapping = pd.read_csv(MAPPING_FILE)
print(f"Municipalities: {len(df_mapping)}")
print(f"Immediate regions: {df_mapping['immediate_code'].nunique()}")
print(f"Intermediate regions: {df_mapping['intermediate_code'].nunique()}")

# Create 6-digit mapping (SIM uses 6-digit codes)
df_mapping['code_muni_6'] = df_mapping['code_muni'].astype(str).str[:6].astype(int)
muni_to_immediate = dict(zip(df_mapping['code_muni_6'], df_mapping['immediate_code']))
muni_to_immediate_name = dict(zip(df_mapping['code_muni_6'], df_mapping['immediate_name']))

# Also keep intermediate for reference
muni_to_intermediate = dict(zip(df_mapping['code_muni_6'], df_mapping['intermediate_code']))

print(f"Mapping ready: {len(muni_to_immediate)} municipalities → 510 immediate regions")

# =============================================================================
# FIND MORTALITY FILES
# =============================================================================

print("\n" + "-"*70)
print("Finding Mortality Files")
print("-"*70)

mortality_files = []
for year in range(2010, 2025):
    yy = str(year)[-2:]
    file_path = INPUT_DIR / f'DO{yy}OPEN.csv'
    if file_path.exists():
        mortality_files.append((year, file_path))

print(f"Found {len(mortality_files)} mortality files")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def parse_sim_age(idade):
    """Parse SIM IDADE field to age in years.
    
    Format: XYY where X is unit and YY is value
    - 4YY: years (e.g., 465 = 65 years)
    - 5YY: 100+ years (e.g., 505 = 105 years)
    
    Note: Ages <1 year (codes 0-3xx) are collapsed to 0 for adult mortality analyses.
    """
    try:
        idade = int(idade)
        if idade >= 400 and idade < 500:
            return idade - 400
        elif idade >= 500:
            return idade - 400
        else:
            return 0  # Under 1 year (infant)
    except:
        return np.nan

def parse_sim_date(dtobito):
    """Parse SIM DTOBITO field to date."""
    try:
        s = str(dtobito).strip()
        if '-' in s:
            return pd.Timestamp(s).date()
        s = str(int(float(dtobito))).zfill(8)
        day = int(s[0:2])
        month = int(s[2:4])
        year = int(s[4:8])
        return pd.Timestamp(year=year, month=month, day=day).date()
    except:
        return None

# =============================================================================
# PROCESS MORTALITY FILES
# =============================================================================

print("\n" + "-"*70)
print("Processing Mortality Files")
print("-"*70)

all_daily = []
all_daily_elderly = []

for year, file_path in mortality_files:
    print(f"\n[{year}] Processing {file_path.name}...")
    
    try:
        # Detect separator
        with open(file_path, 'r', encoding='latin1') as f:
            first_line = f.readline()
        sep = ',' if '","' in first_line or first_line.count(',') > first_line.count(';') else ';'
        
        # Read CSV
        df = pd.read_csv(
            file_path, sep=sep, encoding='latin1',
            usecols=['DTOBITO', 'CODMUNRES', 'IDADE', 'CAUSABAS', 'SEXO'],
            dtype={'DTOBITO': str, 'CODMUNRES': str, 'IDADE': str, 'CAUSABAS': str}
        )
        
        print(f"  Loaded: {len(df):,} deaths")
        
        # Parse date
        df['date'] = df['DTOBITO'].apply(parse_sim_date)
        date_failures = df['date'].isna().sum()
        if date_failures > 0:
            print(f"  ⚠ Date parse failures: {date_failures:,} ({date_failures/len(df)*100:.2f}%)")
        df = df.dropna(subset=['date'])
        
        # Parse age
        df['age'] = df['IDADE'].apply(parse_sim_age)
        
        # Parse municipality code
        df['muni_code'] = pd.to_numeric(df['CODMUNRES'], errors='coerce')
        df = df.dropna(subset=['muni_code'])
        df['muni_code'] = df['muni_code'].astype(int)
        
        # Filter sentinel values (e.g., 999999 for unknown municipality)
        sentinel_count = (df['muni_code'] >= 999990).sum()
        if sentinel_count > 0:
            print(f"  ⚠ Filtering {sentinel_count:,} records with sentinel muni codes (>=999990)")
            df = df[df['muni_code'] < 999990]
        
        # Map to immediate region
        df['immediate_code'] = df['muni_code'].map(muni_to_immediate)
        df['immediate_name'] = df['muni_code'].map(muni_to_immediate_name)
        df['intermediate_code'] = df['muni_code'].map(muni_to_intermediate)
        
        matched = df['immediate_code'].notna().sum()
        total_before_map = len(df)
        match_rate = matched/total_before_map*100
        print(f"  Matched to regions: {matched:,} ({match_rate:.1f}%)")
        
        if match_rate < 90:
            print(f"  ⚠ Warning: Low mapping rate may indicate code format issues")
        
        df = df.dropna(subset=['immediate_code'])
        df['immediate_code'] = df['immediate_code'].astype(int)
        
        # Classify causes
        df['respiratory'] = df['CAUSABAS'].str.startswith('J', na=False).astype(int)
        df['cardiovascular'] = df['CAUSABAS'].str.startswith('I', na=False).astype(int)
        df['heat_direct'] = (
            df['CAUSABAS'].str.startswith('T67', na=False) |
            df['CAUSABAS'].str.startswith('X30', na=False)
        ).astype(int)
        
        # All-cause by immediate region-day
        df_all = df.groupby(['date', 'immediate_code', 'immediate_name']).agg(
            deaths_all=('muni_code', 'count'),
            deaths_respiratory=('respiratory', 'sum'),
            deaths_cardiovascular=('cardiovascular', 'sum'),
            deaths_heat_direct=('heat_direct', 'sum')
        ).reset_index()
        
        all_daily.append(df_all)
        
        # Elderly (60+)
        df_elderly = df[df['age'] >= 60].copy()
        df_elderly_agg = df_elderly.groupby(['date', 'immediate_code', 'immediate_name']).agg(
            deaths_elderly=('muni_code', 'count'),
            deaths_elderly_resp=('respiratory', 'sum'),
            deaths_elderly_cvd=('cardiovascular', 'sum'),
            deaths_elderly_heat=('heat_direct', 'sum')
        ).reset_index()
        
        all_daily_elderly.append(df_elderly_agg)
        
        print(f"  Region-days: {len(df_all):,} (all), {len(df_elderly_agg):,} (elderly)")
        
    except Exception as e:
        print(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()

# =============================================================================
# COMBINE AND SAVE
# =============================================================================

print("\n" + "-"*70)
print("Combining Results")
print("-"*70)

if all_daily:
    df_all = pd.concat(all_daily, ignore_index=True)
    df_all = df_all.groupby(['date', 'immediate_code', 'immediate_name']).sum().reset_index()
    df_all = df_all.sort_values(['immediate_code', 'date']).reset_index(drop=True)
    
    print(f"\nAll-ages mortality:")
    print(f"  Total observations: {len(df_all):,}")
    print(f"  Immediate regions: {df_all['immediate_code'].nunique()}")
    print(f"  Date range: {df_all['date'].min()} to {df_all['date'].max()}")
    print(f"  Total deaths: {df_all['deaths_all'].sum():,}")
    print(f"  Mean deaths/region-day: {df_all['deaths_all'].mean():.1f}")
    
    df_all.to_parquet(OUTPUT_FILE, index=False)
    print(f"  ✓ Saved: {OUTPUT_FILE}")

if all_daily_elderly:
    df_elderly = pd.concat(all_daily_elderly, ignore_index=True)
    df_elderly = df_elderly.groupby(['date', 'immediate_code', 'immediate_name']).sum().reset_index()
    df_elderly = df_elderly.sort_values(['immediate_code', 'date']).reset_index(drop=True)
    
    print(f"\nElderly (60+) mortality:")
    print(f"  Total observations: {len(df_elderly):,}")
    print(f"  Immediate regions: {df_elderly['immediate_code'].nunique()}")
    print(f"  Total deaths: {df_elderly['deaths_elderly'].sum():,}")
    print(f"  Mean deaths/region-day: {df_elderly['deaths_elderly'].mean():.1f}")
    
    df_elderly.to_parquet(OUTPUT_ELDERLY, index=False)
    print(f"  ✓ Saved: {OUTPUT_ELDERLY}")

# Summary by year
print("\n" + "-"*70)
print("Summary by Year")
print("-"*70)

if all_daily:
    df_all['year'] = pd.to_datetime(df_all['date']).dt.year
    yearly = df_all.groupby('year').agg({
        'deaths_all': 'sum',
        'date': 'nunique',
        'immediate_code': 'nunique'
    }).rename(columns={'date': 'days', 'immediate_code': 'regions'})
    print(yearly)

print(f"\n{'='*70}")
print("DONE!")
print("="*70)
