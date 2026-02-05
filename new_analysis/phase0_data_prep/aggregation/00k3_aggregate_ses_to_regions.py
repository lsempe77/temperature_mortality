"""
00k3: Aggregate SES Data to Intermediate and Immediate Regions
===============================================================
Aggregates municipality-level SES covariates to both:
- 133 Intermediate Regions
- 510 Immediate Regions

Input:
- results/municipality_ses_covariates.csv (5,570 municipalities)
- results/municipality_to_region_map.csv (municipality → region mapping)

Output:
- results/ses_intermediate_covariates.csv (133 regions)
- results/ses_immediate_covariates.csv (510 regions)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

print("="*70)
print("00k3: AGGREGATE SES TO REGIONS")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# PATHS
# =============================================================================

BASE_DIR = Path(__file__).resolve().parent
# Data files and output both in main results folder
DATA_DIR = BASE_DIR.parent / 'results'
OUTPUT_DIR = BASE_DIR.parent / 'results'

SES_FILE = DATA_DIR / 'municipality_ses_covariates.csv'
MAPPING_FILE = DATA_DIR / 'municipality_to_all_regions_map.csv'  # Has both intermediate and immediate

OUTPUT_INTERMEDIATE = OUTPUT_DIR / 'ses_intermediate_covariates.csv'
OUTPUT_IMMEDIATE = OUTPUT_DIR / 'ses_immediate_covariates.csv'

# =============================================================================
# LOAD DATA
# =============================================================================

print("\n" + "-"*70)
print("Loading Data")
print("-"*70)

df_ses = pd.read_csv(SES_FILE)
print(f"SES data: {len(df_ses)} municipalities")
print(f"Columns: {list(df_ses.columns)}")

df_mapping = pd.read_csv(MAPPING_FILE)
print(f"Mapping: {len(df_mapping)} municipalities")

# Check columns in mapping
print(f"Mapping columns: {list(df_mapping.columns)}")

# =============================================================================
# MERGE SES WITH MAPPING
# =============================================================================

print("\n" + "-"*70)
print("Merging SES with Region Mapping")
print("-"*70)

# Merge on municipality code
df = pd.merge(
    df_ses,
    df_mapping[['code_muni', 'intermediate_code', 'intermediate_name', 'immediate_code', 'immediate_name']],
    on='code_muni',
    how='left'
)

matched = df['intermediate_code'].notna().sum()
match_rate = matched/len(df)*100
print(f"Matched to intermediate regions: {matched} / {len(df)} ({match_rate:.1f}%)")

matched_imm = df['immediate_code'].notna().sum()
match_rate_imm = matched_imm/len(df)*100
print(f"Matched to immediate regions: {matched_imm} / {len(df)} ({match_rate_imm:.1f}%)")

# Early fail if match rate is too low
MIN_MATCH_RATE = 90.0
if match_rate < MIN_MATCH_RATE:
    print(f"\n✗ ERROR: Match rate {match_rate:.1f}% is below threshold {MIN_MATCH_RATE}%")
    print("  Check that code_muni formats match between SES and mapping files")
    import sys
    sys.exit(1)

# Drop unmatched
df = df.dropna(subset=['intermediate_code'])
df['intermediate_code'] = df['intermediate_code'].astype(int)
df['immediate_code'] = df['immediate_code'].astype(int)

# Guard against zero population (would cause division errors)
zero_pop = (df['pop_total'] == 0).sum()
if zero_pop > 0:
    print(f"\n⚠ Warning: {zero_pop} municipalities have zero population, replacing with NaN")
    df.loc[df['pop_total'] == 0, 'pop_total'] = np.nan

# =============================================================================
# AGGREGATE TO INTERMEDIATE REGIONS (133)
# =============================================================================

print("\n" + "-"*70)
print("Aggregating to Intermediate Regions (133)")
print("-"*70)

# For aggregation:
# - Population metrics: SUM
# - Rates/percentages: weighted average by population
# - GDP per capita: weighted average by population

df_intermediate = df.groupby(['intermediate_code', 'intermediate_name']).agg(
    pop_total=('pop_total', 'sum'),
    pop_elderly=('pop_elderly', 'sum'),
    gdp_total=('gdp_total', 'sum'),
    pop_urban=('pop_urban', 'sum'),
    n_municipalities=('code_muni', 'nunique')  # Use nunique to avoid counting duplicates
).reset_index()

# Calculate rates at region level
df_intermediate['pct_elderly'] = (df_intermediate['pop_elderly'] / df_intermediate['pop_total']) * 100
df_intermediate['gdp_per_capita'] = (df_intermediate['gdp_total'] * 1000) / df_intermediate['pop_total']
df_intermediate['urbanization_rate'] = (df_intermediate['pop_urban'] / df_intermediate['pop_total']) * 100

# Cap urbanization at 100% (some 2010→2022 population changes)
df_intermediate['urbanization_rate'] = df_intermediate['urbanization_rate'].clip(upper=100)

print(f"Intermediate regions: {len(df_intermediate)}")
print(f"\nVariable statistics:")
for col in ['pop_total', 'pop_elderly', 'pct_elderly', 'gdp_per_capita', 'urbanization_rate']:
    print(f"  {col}: mean={df_intermediate[col].mean():.2f}, median={df_intermediate[col].median():.2f}")

# =============================================================================
# AGGREGATE TO IMMEDIATE REGIONS (510)
# =============================================================================

print("\n" + "-"*70)
print("Aggregating to Immediate Regions (510)")
print("-"*70)

df_immediate = df.groupby(['immediate_code', 'immediate_name']).agg(
    pop_total=('pop_total', 'sum'),
    pop_elderly=('pop_elderly', 'sum'),
    gdp_total=('gdp_total', 'sum'),
    pop_urban=('pop_urban', 'sum'),
    n_municipalities=('code_muni', 'nunique')  # Use nunique to avoid counting duplicates
).reset_index()

# Calculate rates at region level
df_immediate['pct_elderly'] = (df_immediate['pop_elderly'] / df_immediate['pop_total']) * 100
df_immediate['gdp_per_capita'] = (df_immediate['gdp_total'] * 1000) / df_immediate['pop_total']
df_immediate['urbanization_rate'] = (df_immediate['pop_urban'] / df_immediate['pop_total']) * 100

# Cap urbanization at 100%
df_immediate['urbanization_rate'] = df_immediate['urbanization_rate'].clip(upper=100)

print(f"Immediate regions: {len(df_immediate)}")
print(f"\nVariable statistics:")
for col in ['pop_total', 'pop_elderly', 'pct_elderly', 'gdp_per_capita', 'urbanization_rate']:
    print(f"  {col}: mean={df_immediate[col].mean():.2f}, median={df_immediate[col].median():.2f}")

# =============================================================================
# SAVE
# =============================================================================

print("\n" + "-"*70)
print("Saving Results")
print("-"*70)

# Ensure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df_intermediate.to_csv(OUTPUT_INTERMEDIATE, index=False)
print(f"Saved: {OUTPUT_INTERMEDIATE}")
print(f"  Shape: {df_intermediate.shape}")

df_immediate.to_csv(OUTPUT_IMMEDIATE, index=False)
print(f"Saved: {OUTPUT_IMMEDIATE}")
print(f"  Shape: {df_immediate.shape}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print(f"\nIntermediate Regions (133):")
print(df_intermediate.head(5).to_string())

print(f"\nImmediate Regions (510):")
print(df_immediate.head(5).to_string())

print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
