"""
00n: AGGREGATE MORTALITY DATA TO INTERMEDIATE AND IMMEDIATE REGIONS
===================================================================
Process SIM (Sistema de Informação sobre Mortalidade) data and aggregate
to both intermediate (133) and immediate (510) region-day levels for DLNM analysis.

Input: 
- Input_data/DO*.csv (yearly mortality files, 2010-2024)
- results/municipality_to_all_regions_map.csv (municipality to both region levels mapping)

Output:
- results/mortality_intermediate_daily.parquet (133 intermediate regions × ~5,500 days)
- results/mortality_intermediate_daily_elderly.parquet (elderly 60+ only, intermediate)
- results/mortality_immediate_daily.parquet (510 immediate regions × ~5,500 days)
- results/mortality_immediate_daily_elderly.parquet (elderly 60+ only, immediate)

Key variables:
- DTOBITO: Death date (DDMMYYYY format)
- CODMUNRES: Municipality of residence (7-digit IBGE code)
- IDADE: Age (4xx = years, e.g., 465 = 65 years old)
- CAUSABAS: Underlying cause (ICD-10)

Focus: All-cause mortality (J00-J99, I00-I99 for sensitivity)
Target: Elderly population (60+)
"""

import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np

print("="*70)
print("00n: AGGREGATE MORTALITY TO INTERMEDIATE AND IMMEDIATE REGIONS")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# PATHS
# =============================================================================

BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR.parent.parent.parent / 'Input_data'
OUTPUT_DIR = BASE_DIR.parent / 'results'  # Main results folder
MAIN_RESULTS = BASE_DIR.parent / 'results'

MAPPING_FILE = MAIN_RESULTS / 'municipality_to_all_regions_map.csv'

# Output files - intermediate (133 regions)
OUTPUT_INTERMEDIATE = OUTPUT_DIR / 'mortality_intermediate_daily.parquet'
OUTPUT_INTERMEDIATE_ELDERLY = OUTPUT_DIR / 'mortality_intermediate_daily_elderly.parquet'
# Output files - immediate (510 regions)
OUTPUT_IMMEDIATE = OUTPUT_DIR / 'mortality_immediate_daily.parquet'
OUTPUT_IMMEDIATE_ELDERLY = OUTPUT_DIR / 'mortality_immediate_daily_elderly.parquet'

# =============================================================================
# LOAD REGION MAPPING
# =============================================================================

print("\n" + "-"*70)
print("Loading Region Mapping")
print("-"*70)

df_mapping = pd.read_csv(MAPPING_FILE)

# Create mapping dicts for both region levels
# SIM uses 6-digit codes, mapping has 7-digit codes

# Check code format
sample_codes = df_mapping['code_muni'].head()
print(f"Mapping code format sample: {sample_codes.tolist()}")

# Create 6-digit version (drop last digit - verification digit)
df_mapping['code_muni_6'] = df_mapping['code_muni'].astype(str).str[:6].astype(int)

# Create mappings for both region levels using 6-digit codes
muni_to_intermediate_6 = dict(zip(df_mapping['code_muni_6'], df_mapping['intermediate_code']))
muni_to_immediate_6 = dict(zip(df_mapping['code_muni_6'], df_mapping['immediate_code']))

# Also create 7-digit mappings for fallback
muni_to_intermediate_7 = dict(zip(df_mapping['code_muni'], df_mapping['intermediate_code']))
muni_to_immediate_7 = dict(zip(df_mapping['code_muni'], df_mapping['immediate_code']))

print(f"Mappings created:")
print(f"  Intermediate: {len(muni_to_intermediate_6)} (6-digit), {len(muni_to_intermediate_7)} (7-digit)")
print(f"  Immediate: {len(muni_to_immediate_6)} (6-digit), {len(muni_to_immediate_7)} (7-digit)")
print(f"  Unique intermediate regions: {df_mapping['intermediate_code'].nunique()}")
print(f"  Unique immediate regions: {df_mapping['immediate_code'].nunique()}")

# =============================================================================
# FIND MORTALITY FILES
# =============================================================================

print("\n" + "-"*70)
print("Finding Mortality Files")
print("-"*70)

# Pattern: DO{YY}OPEN.csv where YY = 10, 11, ..., 24
mortality_files = []
for year in range(2010, 2025):
    yy = str(year)[-2:]  # Last 2 digits
    file_path = INPUT_DIR / f'DO{yy}OPEN.csv'
    if file_path.exists():
        mortality_files.append((year, file_path))

print(f"Found {len(mortality_files)} mortality files:")
for year, fp in mortality_files[:3]:
    print(f"  - {year}: {fp.name}")
if len(mortality_files) > 3:
    print(f"  ... and {len(mortality_files)-3} more")

# =============================================================================
# PROCESS MORTALITY FILES
# =============================================================================

print("\n" + "-"*70)
print("Processing Mortality Files")
print("-"*70)

def parse_sim_age(idade):
    """Parse SIM IDADE field to age in years.
    
    Format: XYY where X is unit and YY is value
    - 0YY: minutes
    - 1YY: hours  
    - 2YY: days
    - 3YY: months
    - 4YY: years (e.g., 465 = 65 years)
    - 5YY: 100+ years (e.g., 505 = 105 years)
    
    2021 format: Plain age numbers (e.g., 65 = 65 years)
    
    Note: Ages <1 year (codes 0-3xx) are collapsed to 0 for adult mortality analyses.
    """
    try:
        idade = int(idade)
        if idade >= 400 and idade < 500:
            return idade - 400  # Years (4xx format)
        elif idade >= 500:
            return idade - 400  # 100+ years (5xx format: 500 = 100, etc)
        elif idade < 200:
            # 2021 format: plain numbers 0-150 representing years
            # Distinguish from 0YY/1YY (minutes/hours) by range
            # If < 200, assume it's years (not 1YY hours format)
            return idade
        else:
            return 0  # Under 1 year (infant) or unknown format
    except:
        return np.nan

def parse_sim_date(dtobito):
    """Parse SIM DTOBITO field to date.
    
    Handles multiple formats:
    - DDMMYYYY (numeric, e.g., 21042022)
    - YYYY-MM-DD (ISO format, e.g., 2021-03-23)
    """
    try:
        s = str(dtobito).strip()
        
        # Check if ISO format (contains hyphen)
        if '-' in s:
            return pd.Timestamp(s).date()
        
        # Otherwise parse as DDMMYYYY
        s = str(int(float(dtobito))).zfill(8)
        day = int(s[0:2])
        month = int(s[2:4])
        year = int(s[4:8])
        return pd.Timestamp(year=year, month=month, day=day).date()
    except:
        return None

all_daily = []
all_daily_elderly = []
all_immediate_daily = []
all_immediate_elderly = []

for year, file_path in mortality_files:
    print(f"\n[{year}] Processing {file_path.name}...")
    
    try:
        # Detect separator by reading first line
        with open(file_path, 'r', encoding='latin1') as f:
            first_line = f.readline()
        
        # Choose separator based on content
        if '","' in first_line or first_line.count(',') > first_line.count(';'):
            sep = ','
        else:
            sep = ';'
        
        # Read CSV (SIM uses latin1 encoding)
        # Read all columns first to handle all file formats
        df = pd.read_csv(
            file_path, 
            sep=sep, 
            encoding='latin1',
            dtype=str,
            low_memory=False
        )
        
        # Clean column names
        df.columns = df.columns.str.strip().str.strip('"')
        
        # Select only needed columns
        required_cols = ['DTOBITO', 'CODMUNRES', 'IDADE', 'CAUSABAS', 'SEXO']
        if not all(col in df.columns for col in required_cols):
            print(f"  ⚠ ERROR: Missing required columns")
            continue
        
        df = df[required_cols].copy()
        
        print(f"  Loaded: {len(df):,} deaths")
        
        # Parse date
        df['date'] = df['DTOBITO'].apply(parse_sim_date)
        date_failures = df['date'].isna().sum()
        if date_failures > 0:
            print(f"  ⚠ Date parse failures: {date_failures:,} ({date_failures/len(df)*100:.2f}%)")
        df = df.dropna(subset=['date'])
        print(f"  Valid dates: {len(df):,}")
        
        # Parse age
        df['age'] = df['IDADE'].apply(parse_sim_age)
        
        # Parse municipality code (try as integer)
        df['muni_code'] = pd.to_numeric(df['CODMUNRES'], errors='coerce')
        df = df.dropna(subset=['muni_code'])
        df['muni_code'] = df['muni_code'].astype(int)
        
        # Filter sentinel values (e.g., 999999 for unknown municipality)
        sentinel_count = (df['muni_code'] >= 999990).sum()
        if sentinel_count > 0:
            print(f"  ⚠ Filtering {sentinel_count:,} records with sentinel muni codes (>=999990)")
            df = df[df['muni_code'] < 999990]
        
        # Map to INTERMEDIATE region (try 6-digit first, then 7-digit)
        df['intermediate_code'] = df['muni_code'].map(muni_to_intermediate_6)
        missing = df['intermediate_code'].isna()
        if missing.sum() > 0:
            df.loc[missing, 'intermediate_code'] = df.loc[missing, 'muni_code'].map(muni_to_intermediate_7)
        
        # Map to IMMEDIATE region (try 6-digit first, then 7-digit)
        df['immediate_code'] = df['muni_code'].map(muni_to_immediate_6)
        missing_imm = df['immediate_code'].isna()
        if missing_imm.sum() > 0:
            df.loc[missing_imm, 'immediate_code'] = df.loc[missing_imm, 'muni_code'].map(muni_to_immediate_7)
        
        matched_int = df['intermediate_code'].notna().sum()
        matched_imm = df['immediate_code'].notna().sum()
        total_before_map = len(df)
        print(f"  Matched to intermediate: {matched_int:,} ({matched_int/total_before_map*100:.1f}%)")
        print(f"  Matched to immediate: {matched_imm:,} ({matched_imm/total_before_map*100:.1f}%)")
        
        if matched_int/total_before_map*100 < 90:
            print(f"  ⚠ Warning: Low mapping rate may indicate code format issues")
        
        # Classify causes (for future sensitivity analysis)
        df['respiratory'] = df['CAUSABAS'].str.startswith('J', na=False).astype(int)
        df['cardiovascular'] = df['CAUSABAS'].str.startswith('I', na=False).astype(int)
        df['heat_direct'] = (
            df['CAUSABAS'].str.startswith('T67', na=False) |  # Heat effects
            df['CAUSABAS'].str.startswith('X30', na=False)    # Exposure to excessive heat
        ).astype(int)
        
        # --- INTERMEDIATE REGION AGGREGATION ---
        df_int = df.dropna(subset=['intermediate_code']).copy()
        df_int['intermediate_code'] = df_int['intermediate_code'].astype(int)
        
        # All-cause mortality by intermediate region-day
        df_int_all = df_int.groupby(['date', 'intermediate_code']).agg(
            deaths_all=('muni_code', 'count'),
            deaths_respiratory=('respiratory', 'sum'),
            deaths_cardiovascular=('cardiovascular', 'sum'),
            deaths_heat_direct=('heat_direct', 'sum')
        ).reset_index()
        
        all_daily.append(df_int_all)
        
        # Elderly (60+) mortality - intermediate
        df_int_elderly = df_int[df_int['age'] >= 60].copy()
        df_int_elderly_agg = df_int_elderly.groupby(['date', 'intermediate_code']).agg(
            deaths_elderly=('muni_code', 'count'),
            deaths_elderly_resp=('respiratory', 'sum'),
            deaths_elderly_cvd=('cardiovascular', 'sum'),
            deaths_elderly_heat=('heat_direct', 'sum')
        ).reset_index()
        
        all_daily_elderly.append(df_int_elderly_agg)
        
        # --- IMMEDIATE REGION AGGREGATION ---
        df_imm = df.dropna(subset=['immediate_code']).copy()
        df_imm['immediate_code'] = df_imm['immediate_code'].astype(int)
        
        # All-cause mortality by immediate region-day
        df_imm_all = df_imm.groupby(['date', 'immediate_code']).agg(
            deaths_all=('muni_code', 'count'),
            deaths_respiratory=('respiratory', 'sum'),
            deaths_cardiovascular=('cardiovascular', 'sum'),
            deaths_heat_direct=('heat_direct', 'sum')
        ).reset_index()
        
        all_immediate_daily.append(df_imm_all)
        
        # Elderly (60+) mortality - immediate
        df_imm_elderly = df_imm[df_imm['age'] >= 60].copy()
        df_imm_elderly_agg = df_imm_elderly.groupby(['date', 'immediate_code']).agg(
            deaths_elderly=('muni_code', 'count'),
            deaths_elderly_resp=('respiratory', 'sum'),
            deaths_elderly_cvd=('cardiovascular', 'sum'),
            deaths_elderly_heat=('heat_direct', 'sum')
        ).reset_index()
        
        all_immediate_elderly.append(df_imm_elderly_agg)
        
        print(f"  Intermediate region-days: {len(df_int_all):,} (all), {len(df_int_elderly_agg):,} (elderly)")
        print(f"  Immediate region-days: {len(df_imm_all):,} (all), {len(df_imm_elderly_agg):,} (elderly)")
        print(f"  Deaths: {df_int_all['deaths_all'].sum():,} (all), {df_int_elderly_agg['deaths_elderly'].sum():,} (60+)")
        
    except Exception as e:
        print(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        continue

# =============================================================================
# COMBINE AND SAVE
# =============================================================================

print("\n" + "-"*70)
print("Combining Results")
print("-"*70)

# --- INTERMEDIATE REGION OUTPUTS ---
print("\n=== INTERMEDIATE REGIONS (133) ===")

if all_daily:
    # All ages - intermediate
    df_int_all = pd.concat(all_daily, ignore_index=True)
    
    # Aggregate duplicates (same date-region across files)
    df_int_all = df_int_all.groupby(['date', 'intermediate_code']).sum().reset_index()
    df_int_all = df_int_all.sort_values(['intermediate_code', 'date']).reset_index(drop=True)
    
    print(f"\nIntermediate - All-ages mortality:")
    print(f"  Total observations: {len(df_int_all):,}")
    print(f"  Regions: {df_int_all['intermediate_code'].nunique()}")
    print(f"  Date range: {df_int_all['date'].min()} to {df_int_all['date'].max()}")
    print(f"  Total deaths: {df_int_all['deaths_all'].sum():,}")
    print(f"  Mean deaths/region-day: {df_int_all['deaths_all'].mean():.1f}")
    
    df_int_all.to_parquet(OUTPUT_INTERMEDIATE, index=False)
    print(f"  ✓ Saved: {OUTPUT_INTERMEDIATE.name}")

if all_daily_elderly:
    # Elderly only - intermediate
    df_int_elderly = pd.concat(all_daily_elderly, ignore_index=True)
    df_int_elderly = df_int_elderly.groupby(['date', 'intermediate_code']).sum().reset_index()
    df_int_elderly = df_int_elderly.sort_values(['intermediate_code', 'date']).reset_index(drop=True)
    
    print(f"\nIntermediate - Elderly (60+) mortality:")
    print(f"  Total observations: {len(df_int_elderly):,}")
    print(f"  Regions: {df_int_elderly['intermediate_code'].nunique()}")
    print(f"  Date range: {df_int_elderly['date'].min()} to {df_int_elderly['date'].max()}")
    print(f"  Total deaths: {df_int_elderly['deaths_elderly'].sum():,}")
    print(f"  Mean deaths/region-day: {df_int_elderly['deaths_elderly'].mean():.1f}")
    
    df_int_elderly.to_parquet(OUTPUT_INTERMEDIATE_ELDERLY, index=False)
    print(f"  ✓ Saved: {OUTPUT_INTERMEDIATE_ELDERLY.name}")

# --- IMMEDIATE REGION OUTPUTS ---
print("\n=== IMMEDIATE REGIONS (510) ===")

if all_immediate_daily:
    # All ages - immediate
    df_imm_all = pd.concat(all_immediate_daily, ignore_index=True)
    
    df_imm_all = df_imm_all.groupby(['date', 'immediate_code']).sum().reset_index()
    df_imm_all = df_imm_all.sort_values(['immediate_code', 'date']).reset_index(drop=True)
    
    print(f"\nImmediate - All-ages mortality:")
    print(f"  Total observations: {len(df_imm_all):,}")
    print(f"  Regions: {df_imm_all['immediate_code'].nunique()}")
    print(f"  Date range: {df_imm_all['date'].min()} to {df_imm_all['date'].max()}")
    print(f"  Total deaths: {df_imm_all['deaths_all'].sum():,}")
    print(f"  Mean deaths/region-day: {df_imm_all['deaths_all'].mean():.1f}")
    
    df_imm_all.to_parquet(OUTPUT_IMMEDIATE, index=False)
    print(f"  ✓ Saved: {OUTPUT_IMMEDIATE.name}")

if all_immediate_elderly:
    # Elderly only - immediate
    df_imm_elderly = pd.concat(all_immediate_elderly, ignore_index=True)
    df_imm_elderly = df_imm_elderly.groupby(['date', 'immediate_code']).sum().reset_index()
    df_imm_elderly = df_imm_elderly.sort_values(['immediate_code', 'date']).reset_index(drop=True)
    
    print(f"\nImmediate - Elderly (60+) mortality:")
    print(f"  Total observations: {len(df_imm_elderly):,}")
    print(f"  Regions: {df_imm_elderly['immediate_code'].nunique()}")
    print(f"  Date range: {df_imm_elderly['date'].min()} to {df_imm_elderly['date'].max()}")
    print(f"  Total deaths: {df_imm_elderly['deaths_elderly'].sum():,}")
    print(f"  Mean deaths/region-day: {df_imm_elderly['deaths_elderly'].mean():.1f}")
    
    df_imm_elderly.to_parquet(OUTPUT_IMMEDIATE_ELDERLY, index=False)
    print(f"  ✓ Saved: {OUTPUT_IMMEDIATE_ELDERLY.name}")

# =============================================================================
# SUMMARY BY YEAR
# =============================================================================

print("\n" + "-"*70)
print("Summary by Year (Intermediate Regions)")
print("-"*70)

if all_daily:
    df_int_all['year'] = pd.to_datetime(df_int_all['date']).dt.year
    yearly = df_int_all.groupby('year').agg({
        'deaths_all': 'sum',
        'date': 'nunique',
        'intermediate_code': 'nunique'
    }).rename(columns={'date': 'days', 'intermediate_code': 'regions'})
    print(yearly)

print(f"\n{'='*70}")
print("DONE!")
print(f"  Intermediate outputs: {OUTPUT_INTERMEDIATE.name}, {OUTPUT_INTERMEDIATE_ELDERLY.name}")
print(f"  Immediate outputs: {OUTPUT_IMMEDIATE.name}, {OUTPUT_IMMEDIATE_ELDERLY.name}")
print("="*70)
