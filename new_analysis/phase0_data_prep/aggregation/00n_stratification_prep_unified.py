"""
Unified Stratification Data Preparation
Generates mortality data stratified by age, sex, and cause for Phase 4 heterogeneity analyses
Handles 2021/2024 format variations (comma separator, different date/age formats)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import pyarrow.parquet as pq
import pyarrow as pa

# Paths
BASE_DIR = Path(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data")
INPUT_DIR = BASE_DIR / "Input_data"
OUTPUT_DIR = BASE_DIR / "new_analysis" / "phase0_data_prep" / "results"
MAPPING_FILE = OUTPUT_DIR / "municipality_to_all_regions_map.csv"

# Load mappings
print("\n" + "="*80)
print("UNIFIED STRATIFICATION DATA PREPARATION")
print("="*80)
print("\nLoading municipality mappings...")
mappings = pd.read_csv(MAPPING_FILE)
print(f"  Loaded {len(mappings):,} municipality mappings")

def parse_sim_age(idade_str):
    """Parse SIM IDADE field to age in years.
    
    Handles multiple formats:
    - 4xx format: 400-499 = years, 500+ = 100+ years
    - Plain numbers: 0-150 = years directly
    - NaN/invalid: returns None
    """
    try:
        # Remove any whitespace and convert to string
        idade_str = str(idade_str).strip()
        
        # Handle NaN
        if idade_str in ['nan', 'None', '']:
            return None
            
        # Remove comma if present (e.g., "0," becomes "0")
        idade_str = idade_str.replace(',', '')
        
        idade = int(idade_str)
        
        if idade >= 400 and idade < 500:
            return idade - 400  # Years (4xx format)
        elif idade >= 500:
            return idade - 400  # 100+ years (5xx format: 500 = 100, etc)
        elif idade < 200:
            # Plain numbers 0-150 representing years
            return idade
        else:
            return None  # Unknown format
    except:
        return None

def classify_cause(causabas):
    """Classify cause of death by ICD-10 code.
    
    Categories:
    - CVD: I00-I99 (Cardiovascular diseases)
    - Respiratory: J00-J99 (Respiratory diseases)
    - External: V01-Y89 (External causes: accidents, violence)
    - Other: All other causes
    """
    if pd.isna(causabas):
        return 'other'
    
    causabas = str(causabas).strip().upper()
    
    if causabas.startswith('I') and len(causabas) >= 3:
        # Cardiovascular I00-I99
        try:
            num = int(causabas[1:3])
            if 0 <= num <= 99:
                return 'cvd'
        except:
            pass
    elif causabas.startswith('J') and len(causabas) >= 3:
        # Respiratory J00-J99
        try:
            num = int(causabas[1:3])
            if 0 <= num <= 99:
                return 'respiratory'
        except:
            pass
    elif len(causabas) >= 3:
        # External causes V01-Y89
        first_char = causabas[0]
        if first_char in ['V', 'W', 'X', 'Y']:
            return 'external'
    
    return 'other'

# Initialize storage for all stratifications
age_data = {level: [] for level in ['intermediate', 'immediate']}
sex_data = {level: [] for level in ['intermediate', 'immediate']}
cause_data = {level: [] for level in ['intermediate', 'immediate']}

years = range(2010, 2025)  # 2010-2024

print(f"\nProcessing {len(years)} years of mortality data...")
print("="*80)

for year_idx, year in enumerate(years, 1):
    csv_file = INPUT_DIR / f"DO{year % 100:02d}OPEN.csv"
    
    if not csv_file.exists():
        print(f"  [{year_idx}/{len(years)}] Skipping {csv_file.name} (not found)")
        continue
    
    print(f"\n  [{year_idx}/{len(years)}] Processing {csv_file.name} (year {year})...")
    
    try:
        # Auto-detect separator (comma vs semicolon)
        with open(csv_file, 'r', encoding='latin1') as f:
            first_line = f.readline()
            comma_count = first_line.count(',')
            semicolon_count = first_line.count(';')
            sep = ',' if comma_count > semicolon_count else ';'
        
        # Read CSV
        df = pd.read_csv(
            csv_file,
            sep=sep,
            encoding='latin1',
            usecols=['DTOBITO', 'CODMUNRES', 'IDADE', 'SEXO', 'CAUSABAS'],
            low_memory=False
        )
        
        df.columns = df.columns.str.strip().str.strip('"')
        print(f"      Loaded {len(df):,} records")
        
        # Parse date (handle both DDMMYYYY and YYYY-MM-DD formats)
        df['date'] = pd.to_datetime(df['DTOBITO'], format='%d%m%Y', errors='coerce')
        if df['date'].isna().all():
            df['date'] = pd.to_datetime(df['DTOBITO'], errors='coerce')
        df = df.dropna(subset=['date'])
        
        # Parse age
        df['age'] = df['IDADE'].apply(parse_sim_age)
        
        # Keep only elderly (60+)
        df = df[df['age'] >= 60].copy()
        print(f"      Elderly deaths (60+): {len(df):,}")
        
        if len(df) == 0:
            print(f"      WARNING: No elderly deaths for {year}")
            continue
        
        # Parse municipality code (6 digits)
        df['muni_code_6'] = pd.to_numeric(
            df['CODMUNRES'].astype(str).str[:6],
            errors='coerce'
        )
        
        # Merge with mappings
        df = df.merge(
            mappings[['muni_code_6', 'region_code', 'immediate_code']],
            on='muni_code_6',
            how='left'
        )
        
        # Filter to mapped municipalities
        df_intermediate = df[df['region_code'].notna()].copy()
        df_immediate = df[df['immediate_code'].notna()].copy()
        
        print(f"      Mapped to intermediate: {len(df_intermediate):,}")
        print(f"      Mapped to immediate: {len(df_immediate):,}")
        
        # ========== AGE STRATIFICATION ==========
        # Categorize ages: 60-69, 70-79, 80+
        df_intermediate['age_group'] = pd.cut(
            df_intermediate['age'],
            bins=[60, 70, 80, 150],
            labels=['60_69', '70_79', '80plus'],
            right=False
        )
        
        df_immediate['age_group'] = pd.cut(
            df_immediate['age'],
            bins=[60, 70, 80, 150],
            labels=['60_69', '70_79', '80plus'],
            right=False
        )
        
        # Aggregate by region, date, age group
        for level, df_level in [('intermediate', df_intermediate), ('immediate', df_immediate)]:
            region_col = 'region_code' if level == 'intermediate' else 'immediate_code'
            
            age_agg = df_level.groupby([region_col, 'date', 'age_group']).size().reset_index(name='deaths')
            age_agg = age_agg.pivot_table(
                index=[region_col, 'date'],
                columns='age_group',
                values='deaths',
                fill_value=0
            ).reset_index()
            
            age_agg.columns.name = None
            age_agg.rename(columns={
                region_col: 'region_code',
                '60_69': 'deaths_60_69',
                '70_79': 'deaths_70_79',
                '80plus': 'deaths_80plus'
            }, inplace=True)
            
            age_data[level].append(age_agg)
        
        # ========== SEX STRATIFICATION ==========
        # SEXO: 1=Male, 2=Female
        df_intermediate['sex'] = df_intermediate['SEXO'].apply(
            lambda x: 'male' if str(x).strip() == '1' else ('female' if str(x).strip() == '2' else None)
        )
        df_immediate['sex'] = df_immediate['SEXO'].apply(
            lambda x: 'male' if str(x).strip() == '1' else ('female' if str(x).strip() == '2' else None)
        )
        
        # Aggregate by region, date, sex
        for level, df_level in [('intermediate', df_intermediate), ('immediate', df_immediate)]:
            region_col = 'region_code' if level == 'intermediate' else 'immediate_code'
            
            df_sex = df_level[df_level['sex'].notna()].copy()
            sex_agg = df_sex.groupby([region_col, 'date', 'sex']).size().reset_index(name='deaths')
            sex_agg = sex_agg.pivot_table(
                index=[region_col, 'date'],
                columns='sex',
                values='deaths',
                fill_value=0
            ).reset_index()
            
            sex_agg.columns.name = None
            sex_agg.rename(columns={
                region_col: 'region_code',
                'male': 'deaths_male',
                'female': 'deaths_female'
            }, inplace=True)
            
            sex_data[level].append(sex_agg)
        
        # ========== CAUSE STRATIFICATION ==========
        df_intermediate['cause'] = df_intermediate['CAUSABAS'].apply(classify_cause)
        df_immediate['cause'] = df_immediate['CAUSABAS'].apply(classify_cause)
        
        # Aggregate by region, date, cause
        for level, df_level in [('intermediate', df_intermediate), ('immediate', df_immediate)]:
            region_col = 'region_code' if level == 'intermediate' else 'immediate_code'
            
            cause_agg = df_level.groupby([region_col, 'date', 'cause']).size().reset_index(name='deaths')
            cause_agg = cause_agg.pivot_table(
                index=[region_col, 'date'],
                columns='cause',
                values='deaths',
                fill_value=0
            ).reset_index()
            
            cause_agg.columns.name = None
            cause_agg.rename(columns={
                region_col: 'region_code',
                'cvd': 'deaths_cvd',
                'respiratory': 'deaths_respiratory',
                'external': 'deaths_external',
                'other': 'deaths_other'
            }, inplace=True)
            
            cause_data[level].append(cause_agg)
        
    except Exception as e:
        print(f"      ERROR: {e}")
        continue

print("\n" + "="*80)
print("COMBINING AND SAVING RESULTS")
print("="*80)

# ========== SAVE AGE STRATIFICATION ==========
print("\nAge Stratification:")
for level in ['intermediate', 'immediate']:
    if age_data[level]:
        df_combined = pd.concat(age_data[level], ignore_index=True)
        
        # Aggregate duplicates (if any)
        df_combined = df_combined.groupby(['region_code', 'date']).agg({
            'deaths_60_69': 'sum',
            'deaths_70_79': 'sum',
            'deaths_80plus': 'sum'
        }).reset_index()
        
        output_file = OUTPUT_DIR / f"mortality_by_age_{level}.parquet"
        df_combined.to_parquet(output_file, index=False)
        
        print(f"  {level.capitalize()}: {len(df_combined):,} rows saved to {output_file.name}")
        print(f"    Total deaths 60-69: {df_combined['deaths_60_69'].sum():,}")
        print(f"    Total deaths 70-79: {df_combined['deaths_70_79'].sum():,}")
        print(f"    Total deaths 80+: {df_combined['deaths_80plus'].sum():,}")

# ========== SAVE SEX STRATIFICATION ==========
print("\nSex Stratification:")
for level in ['intermediate', 'immediate']:
    if sex_data[level]:
        df_combined = pd.concat(sex_data[level], ignore_index=True)
        
        # Aggregate duplicates
        df_combined = df_combined.groupby(['region_code', 'date']).agg({
            'deaths_male': 'sum',
            'deaths_female': 'sum'
        }).reset_index()
        
        output_file = OUTPUT_DIR / f"mortality_by_sex_{level}.parquet"
        df_combined.to_parquet(output_file, index=False)
        
        print(f"  {level.capitalize()}: {len(df_combined):,} rows saved to {output_file.name}")
        print(f"    Total deaths male: {df_combined['deaths_male'].sum():,}")
        print(f"    Total deaths female: {df_combined['deaths_female'].sum():,}")

# ========== SAVE CAUSE STRATIFICATION ==========
print("\nCause Stratification:")
for level in ['intermediate', 'immediate']:
    if cause_data[level]:
        df_combined = pd.concat(cause_data[level], ignore_index=True)
        
        # Aggregate duplicates
        df_combined = df_combined.groupby(['region_code', 'date']).agg({
            'deaths_cvd': 'sum',
            'deaths_respiratory': 'sum',
            'deaths_external': 'sum',
            'deaths_other': 'sum'
        }).reset_index()
        
        output_file = OUTPUT_DIR / f"mortality_by_cause_{level}.parquet"
        df_combined.to_parquet(output_file, index=False)
        
        print(f"  {level.capitalize()}: {len(df_combined):,} rows saved to {output_file.name}")
        print(f"    Total deaths CVD: {df_combined['deaths_cvd'].sum():,}")
        print(f"    Total deaths respiratory: {df_combined['deaths_respiratory'].sum():,}")
        print(f"    Total deaths external: {df_combined['deaths_external'].sum():,}")
        print(f"    Total deaths other: {df_combined['deaths_other'].sum():,}")

print("\n" + "="*80)
print("UNIFIED STRATIFICATION PREP COMPLETE")
print("="*80)
print("\nGenerated 6 parquet files:")
print("  - mortality_by_age_{intermediate,immediate}.parquet")
print("  - mortality_by_sex_{intermediate,immediate}.parquet")
print("  - mortality_by_cause_{intermediate,immediate}.parquet")
