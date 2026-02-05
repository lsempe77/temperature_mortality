"""
00m2: AGGREGATE INFLUENZA DATA AT MUNICIPAL/INTERMEDIATE REGION LEVEL
======================================================================
Re-aggregate the raw influenza/SRAG data at municipality and intermediate 
region levels for micro-level analysis.

Input: results/influenza_raw_elderly.parquet (from 00m_download_influenza_data.py)
       results/municipality_to_region_map.csv (from 00j_download_ibge_mapping.py)

Output:
- results/influenza_daily_by_municipality.parquet
- results/influenza_daily_by_intermediate_region.parquet
- results/influenza_weekly_by_intermediate_region.parquet
"""

import sys
from datetime import datetime
from pathlib import Path

print("=" * 70)
print("00m2: AGGREGATE INFLUENZA DATA AT MUNICIPAL LEVEL")
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

BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR.parent / 'results'  # Main results folder
RESULTS_DIR = BASE_DIR.parent / 'results'  # For mapping file

# Input files
RAW_FILE = RESULTS_DIR / 'influenza_raw_elderly.parquet'
# Use the unified mapping file with both intermediate AND immediate codes
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
print(f"  → {mapping['intermediate_code'].nunique()} intermediate regions")
print(f"  → {mapping['immediate_code'].nunique()} immediate regions")

# =============================================================================
# STEP 2: ADD MUNICIPALITY CODE
# =============================================================================

print("\n" + "-" * 70)
print("STEP 2: Map municipality codes to intermediate regions")
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
required_cols = ['code_muni_6', 'code_muni', 'abbrev_state', 'intermediate_code', 'intermediate_name', 'immediate_code', 'immediate_name']
missing_cols = [c for c in required_cols if c not in mapping.columns]
if missing_cols:
    print(f"✗ ERROR: Mapping file missing columns: {missing_cols}")
    print(f"  Available columns: {list(mapping.columns)}")
    import sys
    sys.exit(1)

# Merge with mapping
raw_mapped = raw.merge(
    mapping[required_cols],
    on='code_muni_6',
    how='left'
)

# Check merge success
matched = raw_mapped['intermediate_code'].notna().sum()
match_rate = 100*matched/len(raw_mapped)
print(f"\nMatched to intermediate regions: {matched:,} / {len(raw_mapped):,} ({match_rate:.1f}%)")

MIN_MATCH_RATE = 85.0
if match_rate < MIN_MATCH_RATE:
    print(f"✗ ERROR: Match rate {match_rate:.1f}% is below threshold {MIN_MATCH_RATE}%")
    unmatched_codes = raw_mapped[raw_mapped['intermediate_code'].isna()]['code_muni_6'].dropna().unique()
    # Log unmatched codes to file for investigation
    unmatched_file = OUTPUT_DIR / 'unmatched_influenza_muni_codes.txt'
    with open(unmatched_file, 'w') as f:
        for code in unmatched_codes:
            f.write(f"{code}\n")
    print(f"  Logged {len(unmatched_codes)} unmatched codes to: {unmatched_file.name}")
    import sys
    sys.exit(1)
elif matched < len(raw_mapped) * 0.95:
    print("⚠ Warning: Some unmatched records")
    unmatched_codes = raw_mapped[raw_mapped['intermediate_code'].isna()]['code_muni_6'].dropna().unique()[:10]
    print(f"  Sample unmatched codes: {unmatched_codes}")

# =============================================================================
# STEP 3: AGGREGATE AT DIFFERENT LEVELS
# =============================================================================

print("\n" + "-" * 70)
print("STEP 3: Aggregate to different geographic levels")
print("-" * 70)

# Ensure date is datetime
raw_mapped['date'] = pd.to_datetime(raw_mapped['date'])

# --- 3a: Municipality-day aggregation ---
print("\n3a. Municipality-day aggregation...")

muni_daily = raw_mapped.groupby(['code_muni', 'date']).agg(
    srag_cases=('is_elderly', 'sum'),
    srag_influenza=('is_influenza', 'sum'),
    srag_covid=('is_covid', 'sum'),
    srag_deaths=('is_death', 'sum'),
).reset_index()

print(f"  Created {len(muni_daily):,} municipality-day observations")
print(f"  Municipalities with data: {muni_daily['code_muni'].nunique():,}")

# --- 3b: Intermediate region-day aggregation ---
print("\n3b. Intermediate region-day aggregation...")

intermediate_daily = raw_mapped[raw_mapped['intermediate_code'].notna()].groupby(['intermediate_code', 'intermediate_name', 'date']).agg(
    srag_cases=('is_elderly', 'sum'),
    srag_influenza=('is_influenza', 'sum'),
    srag_covid=('is_covid', 'sum'),
    srag_deaths=('is_death', 'sum'),
).reset_index()

print(f"  Created {len(intermediate_daily):,} intermediate region-day observations")
print(f"  Intermediate regions with data: {intermediate_daily['intermediate_code'].nunique()}")

# --- 3c: Immediate region-day aggregation ---
print("\n3c. Immediate region-day aggregation...")

immediate_daily = raw_mapped[raw_mapped['immediate_code'].notna()].groupby(['immediate_code', 'immediate_name', 'date']).agg(
    srag_cases=('is_elderly', 'sum'),
    srag_influenza=('is_influenza', 'sum'),
    srag_covid=('is_covid', 'sum'),
    srag_deaths=('is_death', 'sum'),
).reset_index()

print(f"  Created {len(immediate_daily):,} immediate region-day observations")
print(f"  Immediate regions with data: {immediate_daily['immediate_code'].nunique()}")

# --- 3d: Intermediate region-week aggregation ---
print("\n3d. Intermediate region-week aggregation...")

raw_mapped['epi_week'] = raw_mapped['date'].dt.isocalendar().week
raw_mapped['epi_year'] = raw_mapped['date'].dt.isocalendar().year
raw_mapped['week_start'] = raw_mapped['date'] - pd.to_timedelta(raw_mapped['date'].dt.dayofweek, unit='D')

intermediate_weekly = raw_mapped[raw_mapped['intermediate_code'].notna()].groupby(['intermediate_code', 'intermediate_name', 'epi_year', 'epi_week', 'week_start']).agg(
    srag_cases=('is_elderly', 'sum'),
    srag_influenza=('is_influenza', 'sum'),
    srag_covid=('is_covid', 'sum'),
    srag_deaths=('is_death', 'sum'),
).reset_index()

print(f"  Created {len(intermediate_weekly):,} intermediate region-week observations")

# --- 3e: Immediate region-week aggregation ---
print("\n3e. Immediate region-week aggregation...")

immediate_weekly = raw_mapped[raw_mapped['immediate_code'].notna()].groupby(['immediate_code', 'immediate_name', 'epi_year', 'epi_week', 'week_start']).agg(
    srag_cases=('is_elderly', 'sum'),
    srag_influenza=('is_influenza', 'sum'),
    srag_covid=('is_covid', 'sum'),
    srag_deaths=('is_death', 'sum'),
).reset_index()

print(f"  Created {len(immediate_weekly):,} immediate region-week observations")

# =============================================================================
# STEP 4: SAVE OUTPUTS
# =============================================================================

print("\n" + "-" * 70)
print("STEP 4: Save output files")
print("-" * 70)

# Municipality-day
muni_file = OUTPUT_DIR / 'influenza_daily_by_municipality.parquet'
muni_daily.to_parquet(muni_file, index=False)
print(f"✓ Saved: {muni_file.name} ({len(muni_daily):,} rows)")

# Intermediate region-day
intermediate_daily_file = OUTPUT_DIR / 'influenza_daily_by_intermediate_region.parquet'
intermediate_daily.to_parquet(intermediate_daily_file, index=False)
print(f"✓ Saved: {intermediate_daily_file.name} ({len(intermediate_daily):,} rows)")

# Immediate region-day
immediate_daily_file = OUTPUT_DIR / 'influenza_daily_by_immediate_region.parquet'
immediate_daily.to_parquet(immediate_daily_file, index=False)
print(f"✓ Saved: {immediate_daily_file.name} ({len(immediate_daily):,} rows)")

# Intermediate region-week
intermediate_weekly_file = OUTPUT_DIR / 'influenza_weekly_by_intermediate_region.parquet'
intermediate_weekly.to_parquet(intermediate_weekly_file, index=False)
print(f"✓ Saved: {intermediate_weekly_file.name} ({len(intermediate_weekly):,} rows)")

# Immediate region-week
immediate_weekly_file = OUTPUT_DIR / 'influenza_weekly_by_immediate_region.parquet'
immediate_weekly.to_parquet(immediate_weekly_file, index=False)
print(f"✓ Saved: {immediate_weekly_file.name} ({len(immediate_weekly):,} rows)")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"\nGeographic coverage:")
print(f"  Municipalities: {muni_daily['code_muni'].nunique():,}")
print(f"  Intermediate regions: {intermediate_daily['intermediate_code'].nunique()}")
print(f"  Immediate regions: {immediate_daily['immediate_code'].nunique()}")

print(f"\nTemporal coverage:")
print(f"  Date range: {intermediate_daily['date'].min().date()} to {intermediate_daily['date'].max().date()}")

print(f"\nTotal influenza cases (elderly):")
print(f"  Daily aggregates: {muni_daily['srag_influenza'].sum():,}")

print(f"\n✓ Done! {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
