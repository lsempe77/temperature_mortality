"""
================================================================================
UNIFIED MORTALITY DATA PROCESSOR
================================================================================
Single source of truth for all mortality parquet files.

This script handles ALL format variations in the raw SIM (Sistema de Informação
sobre Mortalidade) CSV files and generates ALL required parquet outputs.

FORMAT VARIATIONS HANDLED:
--------------------------
| Year      | Separator | Date Format  | Age Format | SEXO Format           |
|-----------|-----------|--------------|------------|-----------------------|
| 2010-2020 | ;         | DDMMYYYY     | 4xx coded  | Numeric (1=M, 2=F)    |
| 2021      | ,         | YYYY-MM-DD   | 4xx coded  | Text (Masculino/...)  |
| 2022-2023 | ;         | DDMMYYYY     | 4xx coded  | Numeric (1=M, 2=F)    |
| 2024      | ,         | DDMMYYYY     | 4xx coded  | Numeric (1=M, 2=F)    |

OUTPUTS GENERATED:
------------------
Main files:
  - mortality_intermediate_daily.parquet
  - mortality_intermediate_daily_elderly.parquet
  - mortality_immediate_daily.parquet
  - mortality_immediate_daily_elderly.parquet

Age stratified (elderly only):
  - mortality_{level}_daily_age_60_69.parquet
  - mortality_{level}_daily_age_70_79.parquet
  - mortality_{level}_daily_age_80plus.parquet

Sex stratified (elderly only):
  - mortality_{level}_daily_male.parquet
  - mortality_{level}_daily_female.parquet

Cause stratified (elderly only):
  - mortality_{level}_daily_cvd.parquet
  - mortality_{level}_daily_respiratory.parquet
  - mortality_{level}_daily_external.parquet
  - mortality_{level}_daily_other.parquet

Where {level} = intermediate (133 regions) or immediate (510 regions)

Author: Analysis Pipeline
Created: December 17, 2025
================================================================================
"""

import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np
from typing import Optional, Dict, List, Tuple

print("="*80)
print("UNIFIED MORTALITY DATA PROCESSOR")
print("="*80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR.parent.parent.parent / 'Input_data'
OUTPUT_DIR = BASE_DIR.parent / 'results'
MAPPING_FILE = OUTPUT_DIR / 'municipality_to_all_regions_map.csv'

# Years to process
YEARS = range(2010, 2025)

# Age groups for stratification
AGE_GROUPS = {
    '60_69': (60, 70),
    '70_79': (70, 80),
    '80plus': (80, 200),
}

# Cause of death categories (ICD-10)
CAUSE_CATEGORIES = {
    'cvd': lambda x: str(x).upper().startswith('I'),  # I00-I99
    'respiratory': lambda x: str(x).upper().startswith('J'),  # J00-J99
    'external': lambda x: str(x).upper()[0] in ['V', 'W', 'X', 'Y'],  # V01-Y89
}

# =============================================================================
# PARSING FUNCTIONS
# =============================================================================

def detect_csv_format(file_path: Path) -> Dict:
    """Detect the format of a mortality CSV file."""
    with open(file_path, 'r', encoding='latin1') as f:
        header_line = f.readline()
        first_data_line = f.readline()
    
    # Detect separator
    comma_count = header_line.count(',')
    semicolon_count = header_line.count(';')
    sep = ',' if comma_count > semicolon_count else ';'
    
    # Detect date format from first data line
    parts = first_data_line.split(sep)
    # Find DTOBITO position in header
    header_parts = [h.strip().strip('"') for h in header_line.split(sep)]
    try:
        dtobito_idx = header_parts.index('DTOBITO')
        dtobito_val = parts[dtobito_idx].strip().strip('"')
        date_format = 'iso' if '-' in dtobito_val else 'dmy'
    except:
        date_format = 'dmy'  # Default
    
    # Detect SEXO format
    try:
        sexo_idx = header_parts.index('SEXO')
        sexo_val = parts[sexo_idx].strip().strip('"').lower()
        sexo_format = 'text' if sexo_val in ['masculino', 'feminino', 'ignorado'] else 'numeric'
    except:
        sexo_format = 'numeric'  # Default
    
    return {
        'separator': sep,
        'date_format': date_format,
        'sexo_format': sexo_format,
    }


def parse_age(idade_str) -> Optional[int]:
    """
    Parse SIM IDADE field to age in years.
    
    Format (all years):
    - 0YY: minutes (infant)
    - 1YY: hours (infant)
    - 2YY: days (infant)
    - 3YY: months (infant)
    - 4YY: years (e.g., 465 = 65 years)
    - 5YY: 100+ years (e.g., 505 = 105 years)
    """
    try:
        # Clean and convert
        s = str(idade_str).strip().replace(',', '').replace('"', '')
        if s in ['', 'nan', 'None', 'NA']:
            return None
        
        idade = int(float(s))
        
        if idade >= 500:
            return idade - 400  # 100+ years
        elif idade >= 400:
            return idade - 400  # 0-99 years
        elif idade >= 300:
            return 0  # Months (infant)
        elif idade >= 200:
            return 0  # Days (infant)
        elif idade >= 100:
            return 0  # Hours (infant)
        else:
            return 0  # Minutes (infant) or unknown
            
    except (ValueError, TypeError):
        return None


def parse_date(dtobito_str, date_format: str) -> Optional[pd.Timestamp]:
    """
    Parse SIM DTOBITO field to date.
    
    Handles:
    - DDMMYYYY (e.g., 21042022)
    - YYYY-MM-DD (e.g., 2021-03-23)
    """
    try:
        s = str(dtobito_str).strip().replace('"', '')
        
        if date_format == 'iso' or '-' in s:
            return pd.Timestamp(s)
        else:
            # DDMMYYYY format
            s = s.zfill(8)
            day = int(s[0:2])
            month = int(s[2:4])
            year = int(s[4:8])
            return pd.Timestamp(year=year, month=month, day=day)
    except:
        return None


def parse_sex(sexo_str, sexo_format: str) -> Optional[str]:
    """
    Parse SIM SEXO field to standardized sex value.
    
    Returns: 'male', 'female', or None
    """
    try:
        s = str(sexo_str).strip().lower().replace('"', '')
        
        if sexo_format == 'text':
            if s.startswith('masc'):
                return 'male'
            elif s.startswith('fem'):
                return 'female'
            else:
                return None
        else:
            # Numeric format
            if s == '1':
                return 'male'
            elif s == '2':
                return 'female'
            else:
                return None
    except:
        return None


def classify_cause(causabas_str) -> str:
    """
    Classify cause of death by ICD-10 code.
    
    Returns: 'cvd', 'respiratory', 'external', or 'other'
    """
    if pd.isna(causabas_str):
        return 'other'
    
    s = str(causabas_str).strip().upper()
    
    if s.startswith('I'):
        return 'cvd'
    elif s.startswith('J'):
        return 'respiratory'
    elif s[0:1] in ['V', 'W', 'X', 'Y']:
        return 'external'
    else:
        return 'other'


# =============================================================================
# LOAD REGION MAPPINGS
# =============================================================================

print("\n" + "-"*80)
print("Loading Region Mappings")
print("-"*80)

df_mapping = pd.read_csv(MAPPING_FILE)

# Create 6-digit version (drop last digit - verification digit)
df_mapping['muni_code_6'] = df_mapping['code_muni'].astype(str).str[:6].astype(int)

# Create mapping dictionaries (6-digit)
muni_to_intermediate = dict(zip(df_mapping['muni_code_6'], df_mapping['intermediate_code']))
muni_to_immediate = dict(zip(df_mapping['muni_code_6'], df_mapping['immediate_code']))

# Also create 7-digit mappings for fallback
muni_to_intermediate_7 = dict(zip(df_mapping['code_muni'], df_mapping['intermediate_code']))
muni_to_immediate_7 = dict(zip(df_mapping['code_muni'], df_mapping['immediate_code']))

print(f"  Municipalities mapped: {len(muni_to_intermediate):,}")
print(f"  Intermediate regions: {df_mapping['intermediate_code'].nunique()}")
print(f"  Immediate regions: {df_mapping['immediate_code'].nunique()}")


# =============================================================================
# PROCESS ALL MORTALITY FILES
# =============================================================================

print("\n" + "-"*80)
print("Processing Mortality Files")
print("-"*80)

# Storage for all data
all_records = []

for year in YEARS:
    yy = str(year)[-2:]
    file_path = INPUT_DIR / f'DO{yy}OPEN.csv'
    
    if not file_path.exists():
        print(f"  [{year}] ⚠ File not found: {file_path.name}")
        continue
    
    print(f"\n  [{year}] Processing {file_path.name}...")
    
    # Detect format
    fmt = detect_csv_format(file_path)
    print(f"         Format: sep='{fmt['separator']}', date={fmt['date_format']}, sexo={fmt['sexo_format']}")
    
    try:
        # Read CSV
        df = pd.read_csv(
            file_path,
            sep=fmt['separator'],
            encoding='latin1',
            usecols=['DTOBITO', 'CODMUNRES', 'IDADE', 'SEXO', 'CAUSABAS'],
            dtype=str,
            low_memory=False
        )
        
        # Clean column names
        df.columns = df.columns.str.strip().str.strip('"')
        
        print(f"         Raw records: {len(df):,}")
        
        # Parse date
        df['date'] = df['DTOBITO'].apply(lambda x: parse_date(x, fmt['date_format']))
        date_fails = df['date'].isna().sum()
        if date_fails > 0:
            print(f"         ⚠ Date parse failures: {date_fails:,} ({date_fails/len(df)*100:.2f}%)")
        
        # Parse age
        df['age'] = df['IDADE'].apply(parse_age)
        
        # Parse sex
        df['sex'] = df['SEXO'].apply(lambda x: parse_sex(x, fmt['sexo_format']))
        
        # Parse municipality code
        df['muni_code'] = pd.to_numeric(
            df['CODMUNRES'].astype(str).str.replace('"', '').str[:6],
            errors='coerce'
        )
        
        # Classify cause
        df['cause'] = df['CAUSABAS'].apply(classify_cause)
        
        # Drop records with invalid date or municipality
        df = df.dropna(subset=['date', 'muni_code'])
        df['muni_code'] = df['muni_code'].astype(int)
        
        # Filter sentinel municipality codes (999999 = unknown)
        df = df[df['muni_code'] < 999990]
        
        # Map to regions
        df['intermediate_code'] = df['muni_code'].map(muni_to_intermediate)
        df['immediate_code'] = df['muni_code'].map(muni_to_immediate)
        
        # Fallback to 7-digit mapping if needed
        missing_int = df['intermediate_code'].isna()
        if missing_int.sum() > 0:
            df.loc[missing_int, 'intermediate_code'] = df.loc[missing_int, 'muni_code'].map(
                lambda x: muni_to_intermediate_7.get(int(str(x) + '0'), np.nan)
            )
        
        missing_imm = df['immediate_code'].isna()
        if missing_imm.sum() > 0:
            df.loc[missing_imm, 'immediate_code'] = df.loc[missing_imm, 'muni_code'].map(
                lambda x: muni_to_immediate_7.get(int(str(x) + '0'), np.nan)
            )
        
        # Keep only mapped records
        mapped_int = df['intermediate_code'].notna().sum()
        mapped_imm = df['immediate_code'].notna().sum()
        print(f"         Mapped: {mapped_int:,} intermediate, {mapped_imm:,} immediate")
        
        # Select final columns
        df_final = df[['date', 'age', 'sex', 'cause', 'intermediate_code', 'immediate_code']].copy()
        df_final['year'] = year
        
        all_records.append(df_final)
        
        # Summary stats
        elderly = (df_final['age'] >= 60).sum()
        print(f"         Total valid: {len(df_final):,}, Elderly (60+): {elderly:,}")
        
    except Exception as e:
        print(f"         ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        continue

# =============================================================================
# COMBINE ALL DATA
# =============================================================================

print("\n" + "-"*80)
print("Combining All Records")
print("-"*80)

df_all = pd.concat(all_records, ignore_index=True)
print(f"  Total records: {len(df_all):,}")

# Ensure date is datetime
df_all['date'] = pd.to_datetime(df_all['date'])

# Summary by year
yearly_summary = df_all.groupby('year').agg(
    total=('date', 'count'),
    elderly=('age', lambda x: (x >= 60).sum()),
    male=('sex', lambda x: (x == 'male').sum()),
    female=('sex', lambda x: (x == 'female').sum()),
).reset_index()

print("\n  Deaths by Year:")
print(yearly_summary.to_string(index=False))


# =============================================================================
# GENERATE ALL OUTPUTS
# =============================================================================

def save_aggregated(df: pd.DataFrame, region_col: str, output_name: str, 
                    filter_func=None, deaths_col_name: str = 'deaths') -> int:
    """
    Aggregate and save a parquet file.
    
    Returns: total deaths in output
    """
    # Apply filter if provided
    if filter_func is not None:
        df_filtered = df[filter_func(df)].copy()
    else:
        df_filtered = df.copy()
    
    # Drop records without region mapping
    df_filtered = df_filtered.dropna(subset=[region_col])
    df_filtered[region_col] = df_filtered[region_col].astype(int)
    
    # Aggregate by date and region
    df_agg = df_filtered.groupby(['date', region_col]).size().reset_index(name=deaths_col_name)
    df_agg = df_agg.sort_values([region_col, 'date']).reset_index(drop=True)
    
    # Rename region column for consistency
    if region_col == 'intermediate_code':
        df_agg = df_agg.rename(columns={region_col: 'region_code'})
    
    # Save
    output_path = OUTPUT_DIR / output_name
    df_agg.to_parquet(output_path, index=False)
    
    total = df_agg[deaths_col_name].sum()
    return total


print("\n" + "-"*80)
print("Generating Output Files")
print("-"*80)

output_stats = []

# --- MAIN FILES ---
print("\n  === MAIN FILES ===")

for level, region_col in [('intermediate', 'intermediate_code'), ('immediate', 'immediate_code')]:
    # All ages
    total = save_aggregated(
        df_all, region_col, 
        f'mortality_{level}_daily.parquet',
        deaths_col_name='deaths_all'
    )
    output_stats.append((f'mortality_{level}_daily.parquet', total))
    print(f"    mortality_{level}_daily.parquet: {total:,} deaths")
    
    # Elderly (60+)
    total = save_aggregated(
        df_all, region_col,
        f'mortality_{level}_daily_elderly.parquet',
        filter_func=lambda df: df['age'] >= 60,
        deaths_col_name='deaths_elderly'
    )
    output_stats.append((f'mortality_{level}_daily_elderly.parquet', total))
    print(f"    mortality_{level}_daily_elderly.parquet: {total:,} deaths")


# --- AGE STRATIFIED (elderly only) ---
print("\n  === AGE STRATIFIED ===")

for level, region_col in [('intermediate', 'intermediate_code'), ('immediate', 'immediate_code')]:
    for age_group, (age_min, age_max) in AGE_GROUPS.items():
        total = save_aggregated(
            df_all, region_col,
            f'mortality_{level}_daily_age_{age_group}.parquet',
            filter_func=lambda df, amin=age_min, amax=age_max: (df['age'] >= amin) & (df['age'] < amax),
            deaths_col_name='deaths'
        )
        output_stats.append((f'mortality_{level}_daily_age_{age_group}.parquet', total))
        print(f"    mortality_{level}_daily_age_{age_group}.parquet: {total:,} deaths")


# --- SEX STRATIFIED (elderly only) ---
print("\n  === SEX STRATIFIED ===")

for level, region_col in [('intermediate', 'intermediate_code'), ('immediate', 'immediate_code')]:
    for sex in ['male', 'female']:
        total = save_aggregated(
            df_all, region_col,
            f'mortality_{level}_daily_{sex}.parquet',
            filter_func=lambda df, s=sex: (df['age'] >= 60) & (df['sex'] == s),
            deaths_col_name='deaths'
        )
        output_stats.append((f'mortality_{level}_daily_{sex}.parquet', total))
        print(f"    mortality_{level}_daily_{sex}.parquet: {total:,} deaths")


# --- CAUSE STRATIFIED (elderly only) ---
print("\n  === CAUSE STRATIFIED ===")

for level, region_col in [('intermediate', 'intermediate_code'), ('immediate', 'immediate_code')]:
    for cause in ['cvd', 'respiratory', 'external', 'other']:
        total = save_aggregated(
            df_all, region_col,
            f'mortality_{level}_daily_{cause}.parquet',
            filter_func=lambda df, c=cause: (df['age'] >= 60) & (df['cause'] == c),
            deaths_col_name='deaths'
        )
        output_stats.append((f'mortality_{level}_daily_{cause}.parquet', total))
        print(f"    mortality_{level}_daily_{cause}.parquet: {total:,} deaths")


# =============================================================================
# VALIDATION
# =============================================================================

print("\n" + "-"*80)
print("Validation Summary")
print("-"*80)

# Check all files have 2021 and 2024
print("\n  Checking 2021 and 2024 data presence:")

validation_errors = []

for fname, _ in output_stats:
    fpath = OUTPUT_DIR / fname
    df_check = pd.read_parquet(fpath)
    df_check['date'] = pd.to_datetime(df_check['date'])
    df_check['year'] = df_check['date'].dt.year
    
    deaths_col = [c for c in df_check.columns if 'deaths' in c.lower()][0]
    yearly = df_check.groupby('year')[deaths_col].sum()
    
    has_2021 = yearly.get(2021, 0) > 0
    has_2024 = yearly.get(2024, 0) > 0
    
    if not has_2021 or not has_2024:
        validation_errors.append(fname)
        status = "❌ MISSING: "
        if not has_2021:
            status += "2021 "
        if not has_2024:
            status += "2024"
    else:
        status = f"✅ OK (2021: {yearly.get(2021, 0):,}, 2024: {yearly.get(2024, 0):,})"
    
    print(f"    {fname}: {status}")

if validation_errors:
    print(f"\n  ⚠ {len(validation_errors)} files have missing years!")
else:
    print(f"\n  ✅ All {len(output_stats)} files validated successfully!")


# =============================================================================
# FINAL SUMMARY
# =============================================================================

print("\n" + "="*80)
print("PROCESSING COMPLETE")
print("="*80)
print(f"  Total input files processed: {len(YEARS)}")
print(f"  Total records: {len(df_all):,}")
print(f"  Total elderly (60+): {(df_all['age'] >= 60).sum():,}")
print(f"  Output files generated: {len(output_stats)}")
print(f"  Output directory: {OUTPUT_DIR}")
print(f"  Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)
