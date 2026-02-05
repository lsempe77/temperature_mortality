"""
00m3: PROCESS ADDITIONAL INFLUENZA FILES (2019-2024)
=====================================================
Process the INFLUENZA*.csv files that were added separately.
These cover 2019-2024 and have a slightly different structure.

Input:
- temp_srag/INFLUENZA19.csv through INFLUENZA24.csv

Output:
- Append to influenza_daily_by_intermediate_region.parquet
"""

import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np

print("="*70)
print("00m3: PROCESS ADDITIONAL INFLUENZA FILES (2019-2024)")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# PATHS
# =============================================================================

BASE_DIR = Path(__file__).resolve().parent
SRAG_DIR = BASE_DIR.parent / 'temp_srag'  # phase0_data_prep/temp_srag
RESULTS_DIR = BASE_DIR.parent / 'results'  # phase0_data_prep/results
MAIN_RESULTS = BASE_DIR.parent / 'results'

MAPPING_FILE = MAIN_RESULTS / 'municipality_to_all_regions_map.csv'

# =============================================================================
# LOAD MAPPING
# =============================================================================

print("\n" + "-"*70)
print("Loading Region Mapping")
print("-"*70)

df_mapping = pd.read_csv(MAPPING_FILE)

# Create 6-digit mapping (INFLUENZA files use 6-digit codes)
df_mapping['code_muni_6'] = df_mapping['code_muni'].astype(str).str[:6].astype(int)

# Intermediate region mapping
muni_to_inter = dict(zip(df_mapping['code_muni_6'], df_mapping['intermediate_code']))
muni_to_inter_name = dict(zip(df_mapping['code_muni_6'], df_mapping['intermediate_name']))

# Immediate region mapping
muni_to_imm = dict(zip(df_mapping['code_muni_6'], df_mapping['immediate_code']))
muni_to_imm_name = dict(zip(df_mapping['code_muni_6'], df_mapping['immediate_name']))

print(f"Mapping: {len(muni_to_inter)} municipalities")
print(f"  → {df_mapping['intermediate_code'].nunique()} intermediate regions")
print(f"  → {df_mapping['immediate_code'].nunique()} immediate regions")

# =============================================================================
# PROCESS INFLUENZA FILES
# =============================================================================

print("\n" + "-"*70)
print("Processing INFLUENZA Files")
print("-"*70)

influenza_files = sorted(SRAG_DIR.glob('INFLUENZA*.csv'))
print(f"Found {len(influenza_files)} files")

import re

all_records = []

for f in influenza_files:
    # Robust year parsing with regex
    match = re.search(r'INFLUENZA(\d{2})', f.stem, re.IGNORECASE)
    if not match:
        print(f"\n⚠ Skipping {f.name}: cannot parse year from filename")
        continue
    year = int(match.group(1)) + 2000
    print(f"\n[{year}] Processing {f.name}...")
    
    try:
        # Read CSV
        df = pd.read_csv(f, encoding='latin1', sep=';', low_memory=False)
        print(f"  Loaded: {len(df):,} records")
        
        # Verify required columns exist
        required_cols = ['NU_IDADE_N', 'CO_MUN_RES', 'DT_SIN_PRI']
        missing_cols = [c for c in required_cols if c not in df.columns]
        if missing_cols:
            print(f"  ✗ Missing columns: {missing_cols}")
            print(f"    Available: {list(df.columns)[:20]}...")
            continue
        
        # Key columns
        # NU_IDADE_N: age in years (directly)
        # TP_IDADE: type of age (3 = years usually)
        # CO_MUN_RES: municipality of residence (6-digit)
        # DT_SIN_PRI: date of first symptoms
        
        # Filter to elderly (60+)
        df['age'] = pd.to_numeric(df['NU_IDADE_N'], errors='coerce')
        df_elderly = df[df['age'] >= 60].copy()
        print(f"  Elderly (60+): {len(df_elderly):,}")
        
        # Parse date
        df_elderly['date'] = pd.to_datetime(df_elderly['DT_SIN_PRI'], errors='coerce')
        df_elderly = df_elderly.dropna(subset=['date'])
        print(f"  Valid dates: {len(df_elderly):,}")
        
        # Map to region
        df_elderly['muni_code'] = pd.to_numeric(df_elderly['CO_MUN_RES'], errors='coerce')
        df_elderly = df_elderly.dropna(subset=['muni_code'])
        df_elderly['muni_code'] = df_elderly['muni_code'].astype(int)
        
        # Map to both region levels
        df_elderly['intermediate_code'] = df_elderly['muni_code'].map(muni_to_inter)
        df_elderly['intermediate_name'] = df_elderly['muni_code'].map(muni_to_inter_name)
        df_elderly['immediate_code'] = df_elderly['muni_code'].map(muni_to_imm)
        df_elderly['immediate_name'] = df_elderly['muni_code'].map(muni_to_imm_name)
        
        matched = df_elderly['intermediate_code'].notna().sum()
        print(f"  Matched to regions: {matched:,} ({matched/len(df_elderly)*100:.1f}%)")
        
        df_elderly = df_elderly.dropna(subset=['intermediate_code'])
        
        # Keep relevant columns
        df_out = df_elderly[['date', 'muni_code', 'intermediate_code', 'intermediate_name', 
                              'immediate_code', 'immediate_name', 'age']].copy()
        df_out['year'] = year
        
        all_records.append(df_out)
        print(f"  ✓ Added {len(df_out):,} elderly records")
        
    except Exception as e:
        print(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()

# =============================================================================
# AGGREGATE
# =============================================================================

print("\n" + "-"*70)
print("Aggregating Results")
print("-"*70)

if all_records:
    df_all = pd.concat(all_records, ignore_index=True)
    print(f"Total elderly records: {len(df_all):,}")
    print(f"Date range: {df_all['date'].min()} to {df_all['date'].max()}")
    
    # Ensure types
    df_all['intermediate_code'] = df_all['intermediate_code'].astype(int)
    df_all['immediate_code'] = df_all['immediate_code'].astype(int)
    
    # =============================================================================
    # INTERMEDIATE REGION AGGREGATION
    # =============================================================================
    print("\n--- INTERMEDIATE REGIONS (133) ---")
    
    df_daily_inter = df_all.groupby(['date', 'intermediate_code', 'intermediate_name']).agg(
        srag_cases=('muni_code', 'count')
    ).reset_index()
    # Keep intermediate_code column name for consistency with 00m2 output
    
    # Add placeholder columns for compatibility with original format
    df_daily_inter['srag_influenza'] = df_daily_inter['srag_cases']  # INFLUENZA files are confirmed flu
    df_daily_inter['srag_covid'] = 0
    df_daily_inter['srag_deaths'] = 0
    
    print(f"Daily intermediate: {len(df_daily_inter):,} region-day observations")
    print(f"Regions: {df_daily_inter['intermediate_code'].nunique()}")
    
    # Weekly aggregation
    df_all['year'] = df_all['date'].dt.isocalendar().year
    df_all['week'] = df_all['date'].dt.isocalendar().week
    
    df_weekly_inter = df_all.groupby(['year', 'week', 'intermediate_code', 'intermediate_name']).agg(
        srag_cases=('muni_code', 'count'),
        start_date=('date', 'min'),
        end_date=('date', 'max')
    ).reset_index()
    # Keep intermediate_code column name for consistency
    df_weekly_inter['srag_influenza'] = df_weekly_inter['srag_cases']
    df_weekly_inter['srag_covid'] = 0
    df_weekly_inter['srag_deaths'] = 0
    
    print(f"Weekly intermediate: {len(df_weekly_inter):,} region-week observations")
    
    # =============================================================================
    # IMMEDIATE REGION AGGREGATION
    # =============================================================================
    print("\n--- IMMEDIATE REGIONS (510) ---")
    
    df_daily_imm = df_all.groupby(['date', 'immediate_code', 'immediate_name']).agg(
        srag_cases=('muni_code', 'count')
    ).reset_index()
    
    df_daily_imm['srag_influenza'] = df_daily_imm['srag_cases']
    df_daily_imm['srag_covid'] = 0
    df_daily_imm['srag_deaths'] = 0
    
    print(f"Daily immediate: {len(df_daily_imm):,} region-day observations")
    print(f"Regions: {df_daily_imm['immediate_code'].nunique()}")
    
    # =============================================================================
    # MERGE WITH EXISTING DATA - INTERMEDIATE
    # =============================================================================
    
    print("\n" + "-"*70)
    print("Merging with Existing Data - INTERMEDIATE")
    print("-"*70)
    
    # Load existing intermediate
    existing_daily_inter = RESULTS_DIR / 'influenza_daily_by_intermediate_region.parquet'
    existing_weekly_inter = RESULTS_DIR / 'influenza_weekly_by_intermediate_region.parquet'
    
    if existing_daily_inter.exists():
        df_existing_daily = pd.read_parquet(existing_daily_inter)
        print(f"Existing daily: {len(df_existing_daily):,} rows ({df_existing_daily['date'].min()} to {df_existing_daily['date'].max()})")
        
        # Ensure consistent column types
        df_existing_daily['intermediate_code'] = df_existing_daily['intermediate_code'].astype(int)
        df_daily_inter['intermediate_code'] = df_daily_inter['intermediate_code'].astype(int)
        
        # Combine
        df_daily_inter_combined = pd.concat([df_existing_daily, df_daily_inter], ignore_index=True)
        
        # Aggregate any duplicates - sum all numeric columns
        agg_dict = {
            'srag_cases': 'sum',
            'srag_influenza': 'sum',
            'srag_covid': 'sum',
            'srag_deaths': 'sum'
        }
        df_daily_inter_combined = df_daily_inter_combined.groupby(['date', 'intermediate_code', 'intermediate_name']).agg(agg_dict).reset_index()
        
        df_daily_inter_combined = df_daily_inter_combined.sort_values(['intermediate_code', 'date']).reset_index(drop=True)
        print(f"Combined daily: {len(df_daily_inter_combined):,} rows ({df_daily_inter_combined['date'].min()} to {df_daily_inter_combined['date'].max()})")
    else:
        df_daily_inter_combined = df_daily_inter
        print("No existing daily data, using new data only")
    
    if existing_weekly_inter.exists():
        df_existing_weekly = pd.read_parquet(existing_weekly_inter)
        print(f"Existing weekly: {len(df_existing_weekly):,} rows")
        
        # Rename columns if needed for compatibility
        if 'epi_year' in df_existing_weekly.columns:
            df_existing_weekly = df_existing_weekly.rename(columns={'epi_year': 'year', 'epi_week': 'week'})
        if 'week_start' in df_existing_weekly.columns:
            df_existing_weekly = df_existing_weekly.rename(columns={'week_start': 'start_date'})
        # Ensure intermediate_code column exists (handle legacy region_code)
        if 'region_code' in df_existing_weekly.columns and 'intermediate_code' not in df_existing_weekly.columns:
            df_existing_weekly = df_existing_weekly.rename(columns={'region_code': 'intermediate_code', 'region_name': 'intermediate_name'})
        
        # Add end_date if missing
        if 'end_date' not in df_existing_weekly.columns and 'start_date' in df_existing_weekly.columns:
            df_existing_weekly['end_date'] = df_existing_weekly['start_date'] + pd.Timedelta(days=6)
        
        # Ensure consistent types
        df_existing_weekly['intermediate_code'] = df_existing_weekly['intermediate_code'].astype(int)
        df_weekly_inter['intermediate_code'] = df_weekly_inter['intermediate_code'].astype(int)
        
        # Combine
        df_weekly_inter_combined = pd.concat([df_existing_weekly, df_weekly_inter], ignore_index=True)
        
        # Aggregate any duplicates
        agg_cols = {
            'srag_cases': 'sum',
            'srag_influenza': 'sum',
            'srag_covid': 'sum',
            'srag_deaths': 'sum'
        }
        if 'start_date' in df_weekly_inter_combined.columns:
            agg_cols['start_date'] = 'min'
            agg_cols['end_date'] = 'max'
            
        df_weekly_inter_combined = df_weekly_inter_combined.groupby(['year', 'week', 'intermediate_code', 'intermediate_name']).agg(agg_cols).reset_index()
        df_weekly_inter_combined = df_weekly_inter_combined.sort_values(['intermediate_code', 'year', 'week']).reset_index(drop=True)
        print(f"Combined weekly: {len(df_weekly_inter_combined):,} rows")
    else:
        df_weekly_inter_combined = df_weekly_inter
        print("No existing weekly data, using new data only")
    
    # =============================================================================
    # MERGE WITH EXISTING DATA - IMMEDIATE
    # =============================================================================
    
    print("\n" + "-"*70)
    print("Merging with Existing Data - IMMEDIATE")
    print("-"*70)
    
    existing_daily_imm = RESULTS_DIR / 'influenza_daily_by_immediate_region.parquet'
    
    if existing_daily_imm.exists():
        df_existing_imm = pd.read_parquet(existing_daily_imm)
        print(f"Existing daily: {len(df_existing_imm):,} rows ({df_existing_imm['date'].min()} to {df_existing_imm['date'].max()})")
        
        # Ensure consistent column types
        df_existing_imm['immediate_code'] = df_existing_imm['immediate_code'].astype(int)
        df_daily_imm['immediate_code'] = df_daily_imm['immediate_code'].astype(int)
        
        # Combine
        df_daily_imm_combined = pd.concat([df_existing_imm, df_daily_imm], ignore_index=True)
        
        # Aggregate any duplicates
        agg_dict = {
            'srag_cases': 'sum',
            'srag_influenza': 'sum',
            'srag_covid': 'sum',
            'srag_deaths': 'sum'
        }
        df_daily_imm_combined = df_daily_imm_combined.groupby(['date', 'immediate_code', 'immediate_name']).agg(agg_dict).reset_index()
        df_daily_imm_combined = df_daily_imm_combined.sort_values(['immediate_code', 'date']).reset_index(drop=True)
        print(f"Combined daily: {len(df_daily_imm_combined):,} rows ({df_daily_imm_combined['date'].min()} to {df_daily_imm_combined['date'].max()})")
    else:
        df_daily_imm_combined = df_daily_imm
        print("No existing daily data, using new data only")
    
    # =============================================================================
    # SAVE
    # =============================================================================
    
    print("\n" + "-"*70)
    print("Saving Results")
    print("-"*70)
    
    # Save intermediate
    existing_daily_inter.parent.mkdir(parents=True, exist_ok=True)
    df_daily_inter_combined.to_parquet(existing_daily_inter, index=False)
    print(f"✓ Saved: {existing_daily_inter}")
    
    df_weekly_inter_combined.to_parquet(existing_weekly_inter, index=False)
    print(f"✓ Saved: {existing_weekly_inter}")
    
    # Save immediate
    df_daily_imm_combined.to_parquet(existing_daily_imm, index=False)
    print(f"✓ Saved: {existing_daily_imm}")
    
    # Summary by year
    print("\nIntermediate Region Summary by Year:")
    df_daily_inter_combined['year'] = pd.to_datetime(df_daily_inter_combined['date']).dt.year
    yearly = df_daily_inter_combined.groupby('year').agg({
        'srag_cases': 'sum',
        'srag_influenza': 'sum',
        'date': 'nunique',
        'intermediate_code': 'nunique'
    }).rename(columns={'date': 'days', 'intermediate_code': 'regions'})
    print(yearly)
    
    print("\nImmediate Region Summary:")
    print(f"  Total rows: {len(df_daily_imm_combined):,}")
    print(f"  Date range: {df_daily_imm_combined['date'].min()} to {df_daily_imm_combined['date'].max()}")
    print(f"  Regions: {df_daily_imm_combined['immediate_code'].nunique()}")
    print(yearly)

print(f"\n{'='*70}")
print("DONE!")
print("="*70)
