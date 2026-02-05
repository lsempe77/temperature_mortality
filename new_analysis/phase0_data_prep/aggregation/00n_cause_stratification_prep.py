"""
00n_cause_stratification: CREATE CAUSE-STRATIFIED MORTALITY FILES
================================================================
Process SIM mortality data and create cause-stratified files for Phase 4 analysis.

Cause Groups (ICD-10):
- CVD (Cardiovascular): I00-I99
- Respiratory: J00-J99
- External: V01-Y98
- Other: All other causes

Output:
- mortality_intermediate_daily_cvd.parquet
- mortality_intermediate_daily_respiratory.parquet
- mortality_intermediate_daily_external.parquet
- mortality_intermediate_daily_other.parquet
- mortality_immediate_daily_cvd.parquet
- mortality_immediate_daily_respiratory.parquet
- mortality_immediate_daily_external.parquet
- mortality_immediate_daily_other.parquet
"""

import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np

print("="*70)
print("00n_cause_stratification: CREATE CAUSE-STRATIFIED MORTALITY FILES")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# Paths
BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR.parent.parent.parent / 'Input_data'
OUTPUT_DIR = BASE_DIR.parent / 'results'
MAPPING_FILE = OUTPUT_DIR / 'municipality_to_all_regions_map.csv'

# Load mapping
print("\n[1] Loading region mapping...")
df_mapping = pd.read_csv(MAPPING_FILE)
df_mapping['code_muni_6'] = df_mapping['code_muni'].astype(str).str[:6].astype(int)

muni_to_intermediate_6 = dict(zip(df_mapping['code_muni_6'], df_mapping['intermediate_code']))
muni_to_immediate_6 = dict(zip(df_mapping['code_muni_6'], df_mapping['immediate_code']))

print(f"  Municipalities: {len(df_mapping)}")
print(f"  Intermediate regions: {df_mapping['intermediate_code'].nunique()}")
print(f"  Immediate regions: {df_mapping['immediate_code'].nunique()}")

# Find mortality files
print("\n[2] Processing mortality files...")
mort_files = sorted(INPUT_DIR.glob('DO*.csv'))
print(f"  Found {len(mort_files)} files: {min([f.stem for f in mort_files])} to {max([f.stem for f in mort_files])}")

# Initialize storage
def classify_cause(icd_code):
    """Classify ICD-10 code into cause groups"""
    if pd.isna(icd_code):
        return 'other'
    
    code = str(icd_code).strip().upper()
    if not code:
        return 'other'
    
    # Extract letter prefix
    if len(code) >= 1:
        letter = code[0]
        
        # CVD: I00-I99
        if letter == 'I':
            return 'cvd'
        
        # Respiratory: J00-J99
        elif letter == 'J':
            return 'respiratory'
        
        # External: V01-Y98
        elif letter in ['V', 'W', 'X', 'Y']:
            return 'external'
    
    return 'other'

cause_groups = ['cvd', 'respiratory', 'external', 'other']

# Storage for both region levels
intermediate_data = {cause: [] for cause in cause_groups}
immediate_data = {cause: [] for cause in cause_groups}

# Process each file
for i, file in enumerate(mort_files, 1):
    year = file.stem[2:4]
    full_year = f"20{year}"
    
    print(f"\n  [{i}/{len(mort_files)}] Processing {file.name} (year {full_year})...")
    
    try:
        # Auto-detect separator (2021/2024 use comma, others use semicolon)
        with open(file, 'r', encoding='latin1') as f:
            first_line = f.readline().strip()
        
        # Detect separator
        if first_line.count(',') > first_line.count(';'):
            sep = ','
        else:
            sep = ';'
        
        # Read CSV with detected separator
        df = pd.read_csv(
            file,
            encoding='latin1',
            sep=sep,
            dtype=str,
            low_memory=False
        )
        
        # Clean column names (remove quotes, spaces, unnamed columns)
        df.columns = df.columns.str.strip().str.strip('"')
        
        # Select only needed columns
        required_cols = ['DTOBITO', 'CODMUNRES', 'IDADE', 'CAUSABAS']
        if not all(col in df.columns for col in required_cols):
            print(f"      ERROR: Missing required columns. Available: {df.columns.tolist()[:10]}")
            continue
        
        df = df[required_cols].copy()
        
        print(f"      Loaded {len(df):,} records")
        
        # Parse date (handle both DDMMYYYY and YYYY-MM-DD formats)
        df['date'] = pd.to_datetime(df['DTOBITO'], format='%d%m%Y', errors='coerce')
        if df['date'].isna().all():
            # Try YYYY-MM-DD format (2021 format)
            df['date'] = pd.to_datetime(df['DTOBITO'], errors='coerce')
        df = df.dropna(subset=['date'])
        
        # Parse age (handle both 4xx format and plain numbers)
        df['age_str'] = df['IDADE'].astype(str).str.strip()
        # If starts with 4, extract last 2 digits (4xx format)
        # Otherwise use the number as-is (2021 format)
        df['age'] = df['age_str'].apply(lambda x: 
            int(x[1:]) if (len(x) >= 3 and x[0] == '4') else 
            (int(x) if x.isdigit() else None)
        )
        
        # Keep only elderly (60+)
        df = df[df['age'] >= 60].copy()
        
        # Classify cause
        df['cause'] = df['CAUSABAS'].apply(classify_cause)
        
        # Parse municipality code (6 digits)
        df['muni_code_6'] = pd.to_numeric(
            df['CODMUNRES'].str[:6],
            errors='coerce'
        )
        df = df.dropna(subset=['muni_code_6'])
        df['muni_code_6'] = df['muni_code_6'].astype(int)
        
        # Map to regions
        df['intermediate_code'] = df['muni_code_6'].map(muni_to_intermediate_6)
        df['immediate_code'] = df['muni_code_6'].map(muni_to_immediate_6)
        
        # Drop unmapped
        df_inter = df.dropna(subset=['intermediate_code']).copy()
        df_immed = df.dropna(subset=['immediate_code']).copy()
        
        print(f"      Elderly deaths (60+): {len(df):,}")
        print(f"      Mapped to intermediate: {len(df_inter):,}")
        print(f"      Mapped to immediate: {len(df_immed):,}")
        
        # Show cause distribution
        cause_counts = df['cause'].value_counts()
        print(f"      Cause distribution: CVD={cause_counts.get('cvd', 0):,}, "
              f"Respiratory={cause_counts.get('respiratory', 0):,}, "
              f"External={cause_counts.get('external', 0):,}, "
              f"Other={cause_counts.get('other', 0):,}")
        
        # Aggregate by cause for each region level
        for cause_code in cause_groups:
            # Intermediate
            df_cause_inter = df_inter[df_inter['cause'] == cause_code]
            if len(df_cause_inter) > 0:
                agg_inter = df_cause_inter.groupby(['date', 'intermediate_code']).size().reset_index(name='deaths')
                agg_inter.columns = ['date', 'region_code', 'deaths']
                agg_inter['region_code'] = agg_inter['region_code'].astype(int)
                intermediate_data[cause_code].append(agg_inter)
            
            # Immediate
            df_cause_immed = df_immed[df_immed['cause'] == cause_code]
            if len(df_cause_immed) > 0:
                agg_immed = df_cause_immed.groupby(['date', 'immediate_code']).size().reset_index(name='deaths')
                agg_immed.columns = ['date', 'region_code', 'deaths']
                agg_immed['region_code'] = agg_immed['region_code'].astype(int)
                immediate_data[cause_code].append(agg_immed)
        
    except Exception as e:
        print(f"      ERROR: {e}")
        continue

# Combine and save
print("\n[3] Combining and saving cause-stratified files...")

cause_labels = {
    'cvd': 'Cardiovascular (I00-I99)',
    'respiratory': 'Respiratory (J00-J99)',
    'external': 'External (V01-Y98)',
    'other': 'Other causes'
}

for cause_code in cause_groups:
    cause_label = cause_labels[cause_code]
    
    # Intermediate
    if intermediate_data[cause_code]:
        df_combined = pd.concat(intermediate_data[cause_code], ignore_index=True)
        df_combined = df_combined.groupby(['date', 'region_code'])['deaths'].sum().reset_index()
        
        # Fill date range
        date_range = pd.date_range('2010-01-01', '2024-12-31', freq='D')
        regions = sorted(df_combined['region_code'].unique())
        
        full_grid = pd.DataFrame([
            (date, region)
            for date in date_range
            for region in regions
        ], columns=['date', 'region_code'])
        
        df_full = full_grid.merge(df_combined, on=['date', 'region_code'], how='left')
        df_full['deaths'] = df_full['deaths'].fillna(0).astype(int)
        df_full['region_code'] = df_full['region_code'].astype(int)
        df_full['date'] = df_full['date'].dt.strftime('%Y-%m-%d')  # Convert to string to match ERA5
        
        output_file = OUTPUT_DIR / f'mortality_intermediate_daily_{cause_code}.parquet'
        df_full.to_parquet(output_file, index=False)
        
        print(f"  {cause_label}: {output_file.name}")
        print(f"    Regions: {len(regions)}, Records: {len(df_full):,}, Deaths: {df_full['deaths'].sum():,}")
    
    # Immediate
    if immediate_data[cause_code]:
        df_combined = pd.concat(immediate_data[cause_code], ignore_index=True)
        df_combined = df_combined.groupby(['date', 'region_code'])['deaths'].sum().reset_index()
        
        # Fill date range
        date_range = pd.date_range('2010-01-01', '2024-12-31', freq='D')
        regions = sorted(df_combined['region_code'].unique())
        
        full_grid = pd.DataFrame([
            (date, region)
            for date in date_range
            for region in regions
        ], columns=['date', 'region_code'])
        
        df_full = full_grid.merge(df_combined, on=['date', 'region_code'], how='left')
        df_full['deaths'] = df_full['deaths'].fillna(0).astype(int)
        df_full['region_code'] = df_full['region_code'].astype(int)
        df_full['date'] = df_full['date'].dt.strftime('%Y-%m-%d')  # Convert to string to match ERA5
        
        output_file = OUTPUT_DIR / f'mortality_immediate_daily_{cause_code}.parquet'
        df_full.to_parquet(output_file, index=False)
        
        print(f"  {cause_label}: {output_file.name}")
        print(f"    Regions: {len(regions)}, Records: {len(df_full):,}, Deaths: {df_full['deaths'].sum():,}")

print("\n" + "="*70)
print(f"Done! {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*70)
