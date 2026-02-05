"""
00m4: AGGREGATE INFLUENZA DATA TO IMMEDIATE REGIONS (510)
==========================================================
Re-aggregate the raw influenza/SRAG data at immediate region level.

Input: results/influenza_raw_elderly.parquet (from 00m_download_influenza_data.py)
       results/municipality_region_mapping.csv (from 00j_download_ibge_mapping.py)

Output:
- results/influenza_daily_by_immediate_region.parquet
"""

import sys
from datetime import datetime
from pathlib import Path

print("=" * 70)
print("00m4: AGGREGATE INFLUENZA DATA TO IMMEDIATE REGIONS")
print("=" * 70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# DEPENDENCIES
# =============================================================================

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

BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'results'
RESULTS_DIR = BASE_DIR.parent / 'results'  # For mapping file

# Input files
RAW_FILE = OUTPUT_DIR / 'influenza_raw_elderly.parquet'
MAPPING_FILE = RESULTS_DIR / 'municipality_to_all_regions_map.csv'

# =============================================================================
# LOAD DATA
# =============================================================================

print("\n" + "-" * 70)
print("STEP 1: Load raw influenza data and municipality mapping")
print("-" * 70)

# Load raw data
if not RAW_FILE.exists():
    print(f"✗ Raw file not found: {RAW_FILE}")
    print("  Run 00m_download_influenza_data.py first")
    sys.exit(1)

raw = pd.read_parquet(RAW_FILE)
print(f"✓ Loaded raw data: {len(raw):,} elderly SRAG records")

# Load mapping
if not MAPPING_FILE.exists():
    print(f"✗ Mapping file not found: {MAPPING_FILE}")
    print("  Run 00j_download_ibge_mapping.py first")
    sys.exit(1)

mapping = pd.read_csv(MAPPING_FILE)
print(f"✓ Loaded mapping: {len(mapping):,} municipalities")
print(f"  → {mapping['immediate_code'].nunique()} immediate regions")
print(f"  → {mapping['intermediate_code'].nunique()} intermediate regions")

# =============================================================================
# STEP 2: MAP TO IMMEDIATE REGIONS
# =============================================================================

print("\n" + "-" * 70)
print("STEP 2: Map municipality codes to immediate regions")
print("-" * 70)

# The raw data has ID_MN_RESI (municipality of residence)
# This is a 6-digit IBGE code (first 2 = state, next 4 = municipality within state)
# The mapping has 7-digit codes (with check digit)

# Check the format
print(f"Raw ID_MN_RESI sample: {raw['ID_MN_RESI'].dropna().head(5).tolist()}")
print(f"Mapping code_muni sample: {mapping['code_muni'].head(5).tolist()}")

# IBGE codes in SRAG might be 6-digit (without check digit)
# Need to match by first 6 digits
raw['code_muni_6'] = pd.to_numeric(raw['ID_MN_RESI'], errors='coerce').astype('Int64')
mapping['code_muni_6'] = (mapping['code_muni'] // 10).astype('Int64')  # Remove check digit

# Verify required columns exist in mapping before merge
required_cols = ['code_muni_6', 'code_muni', 'abbrev_state', 'immediate_code', 'immediate_name']
missing_cols = [c for c in required_cols if c not in mapping.columns]
if missing_cols:
    print(f"✗ ERROR: Mapping file missing columns: {missing_cols}")
    print(f"  Available columns: {list(mapping.columns)}")
    sys.exit(1)

# Merge with mapping
raw_mapped = raw.merge(
    mapping[required_cols],
    on='code_muni_6',
    how='left'
)

# Check merge success
matched = raw_mapped['immediate_code'].notna().sum()
match_rate = 100*matched/len(raw_mapped)
print(f"\nMatched to immediate regions: {matched:,} / {len(raw_mapped):,} ({match_rate:.1f}%)")

MIN_MATCH_RATE = 85.0
if match_rate < MIN_MATCH_RATE:
    print(f"✗ ERROR: Match rate {match_rate:.1f}% is below threshold {MIN_MATCH_RATE}%")
    unmatched_codes = raw_mapped[raw_mapped['immediate_code'].isna()]['code_muni_6'].dropna().unique()
    # Log unmatched codes to file for investigation
    unmatched_file = OUTPUT_DIR / 'unmatched_influenza_imm_muni_codes.txt'
    with open(unmatched_file, 'w') as f:
        for code in unmatched_codes:
            f.write(f"{code}\n")
    print(f"  Logged {len(unmatched_codes)} unmatched codes to: {unmatched_file.name}")
    sys.exit(1)
elif match_rate < 95:
    print("⚠ Warning: Some unmatched records")
    unmatched_codes = raw_mapped[raw_mapped['immediate_code'].isna()]['code_muni_6'].dropna().unique()[:10]
    print(f"  Sample unmatched codes: {unmatched_codes}")

# =============================================================================
# STEP 3: AGGREGATE TO IMMEDIATE REGION LEVEL
# =============================================================================

print("\n" + "-" * 70)
print("STEP 3: Aggregate to immediate region level")
print("-" * 70)

# Ensure date is datetime
raw_mapped['date'] = pd.to_datetime(raw_mapped['date'])

# Filter to matched records
df = raw_mapped[raw_mapped['immediate_code'].notna()].copy()
df['immediate_code'] = df['immediate_code'].astype(int)

# Immediate region-day aggregation
print("Creating daily aggregation...")

region_daily = df.groupby(['immediate_code', 'immediate_name', 'date']).agg(
    srag_cases=('is_elderly', 'sum'),
    srag_influenza=('is_influenza', 'sum'),
    srag_covid=('is_covid', 'sum'),
    srag_deaths=('is_death', 'sum'),
).reset_index()

print(f"  Created {len(region_daily):,} region-day observations")
print(f"  Immediate regions with data: {region_daily['immediate_code'].nunique()}")

# Stats
print(f"\nDate range: {region_daily['date'].min()} to {region_daily['date'].max()}")
print(f"Total SRAG cases: {region_daily['srag_cases'].sum():,}")
print(f"Total SRAG influenza: {region_daily['srag_influenza'].sum():,}")
print(f"Total SRAG covid: {region_daily['srag_covid'].sum():,}")

# =============================================================================
# STEP 4: SAVE
# =============================================================================

print("\n" + "-" * 70)
print("STEP 4: Save output")
print("-" * 70)

# Save daily data
output_file = OUTPUT_DIR / 'influenza_daily_by_immediate_region.parquet'
region_daily.to_parquet(output_file, index=False)
print(f"✓ Saved: {output_file}")

# Save sample CSV
sample_file = OUTPUT_DIR / 'influenza_daily_by_immediate_region_sample.csv'
region_daily.head(1000).to_csv(sample_file, index=False)
print(f"✓ Saved sample: {sample_file}")

print(f"\n{'='*70}")
print("DONE!")
print("="*70)
