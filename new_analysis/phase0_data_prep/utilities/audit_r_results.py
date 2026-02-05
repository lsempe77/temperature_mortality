"""
Audit R Script Results and Input Data
Check if Phase 1, 2, 3 R scripts used complete data (including 2021 and 2024)
"""

import json
import pandas as pd
from pathlib import Path

BASE_DIR = Path(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis")
RESULTS_DIR = BASE_DIR / "phase0_data_prep" / "results"

print("=" * 80)
print("AUDIT: R SCRIPT RESULTS AND INPUT DATA COMPLETENESS")
print("=" * 80)

# =============================================================================
# 1. CHECK INPUT PARQUET FILES (used by ALL R scripts)
# =============================================================================
print("\n" + "=" * 80)
print("SECTION 1: INPUT PARQUET FILES (Source of truth for all R scripts)")
print("=" * 80)

for level in ['intermediate', 'immediate']:
    print(f"\n--- {level.upper()} ---")
    
    # Mortality
    mort_file = RESULTS_DIR / f'mortality_{level}_daily_elderly.parquet'
    if mort_file.exists():
        mort = pd.read_parquet(mort_file)
        mort['date'] = pd.to_datetime(mort['date'])
        mort['year'] = mort['date'].dt.year
        yearly = mort.groupby('year')['deaths_elderly'].sum()
        
        print(f"Mortality file: {mort_file.name}")
        print(f"  Date range: {mort['date'].min().date()} to {mort['date'].max().date()}")
        print(f"  2021 deaths: {yearly.get(2021, 0):,}")
        print(f"  2024 deaths: {yearly.get(2024, 0):,}")
        print(f"  Total deaths: {yearly.sum():,}")
        
        if yearly.get(2021, 0) == 0:
            print("  ⚠️  WARNING: 2021 MISSING!")
        if yearly.get(2024, 0) == 0:
            print("  ⚠️  WARNING: 2024 MISSING!")
    else:
        print(f"  ❌ File not found: {mort_file.name}")
    
    # ERA5
    era5_file = RESULTS_DIR / f'era5_{level}_daily.parquet'
    if era5_file.exists():
        era5 = pd.read_parquet(era5_file)
        era5['date'] = pd.to_datetime(era5['date'])
        print(f"ERA5 file: {era5_file.name}")
        print(f"  Date range: {era5['date'].min().date()} to {era5['date'].max().date()}")
    else:
        print(f"  ❌ ERA5 file not found: {era5_file.name}")

# =============================================================================
# 2. CHECK PHASE 1 R RESULTS
# =============================================================================
print("\n" + "=" * 80)
print("SECTION 2: PHASE 1 R SCRIPT RESULTS")
print("=" * 80)

phase1_dir = BASE_DIR / "phase1_r" / "results"

# DLNM results
for level in ['intermediate', 'immediate']:
    fname = f'dlnm_r_{level}_results_v2.json'
    fpath = phase1_dir / fname
    print(f"\n--- {fname} ---")
    if fpath.exists():
        with open(fpath, 'r') as f:
            data = json.load(f)
        meta = data.get('metadata', {})
        print(f"  Date range: {meta.get('date_range', 'N/A')}")
        print(f"  Total deaths: {meta.get('total_deaths', 'N/A'):,}" if isinstance(meta.get('total_deaths'), (int, float)) else f"  Total deaths: {meta.get('total_deaths', 'N/A')}")
        print(f"  Observations: {meta.get('n_observations', 'N/A'):,}" if isinstance(meta.get('n_observations'), (int, float)) else f"  Observations: {meta.get('n_observations', 'N/A')}")
        print(f"  Regions: {meta.get('n_regions', 'N/A')}")
    else:
        print(f"  ❌ File not found")

# Attributable burden
for level in ['intermediate', 'immediate']:
    fname = f'attributable_burden_r_{level}.json'
    fpath = phase1_dir / fname
    print(f"\n--- {fname} ---")
    if fpath.exists():
        with open(fpath, 'r') as f:
            data = json.load(f)
        summary = data.get('summary', {})
        print(f"  Total deaths: {summary.get('total_deaths', 'N/A'):,}" if isinstance(summary.get('total_deaths'), (int, float)) else f"  Total deaths: {summary.get('total_deaths', 'N/A')}")
        print(f"  Attributable fraction: {summary.get('total_af', 'N/A')}")
    else:
        print(f"  ❌ File not found")

# YLL
for level in ['intermediate', 'immediate']:
    fname = f'yll_r_{level}.json'
    fpath = phase1_dir / fname
    print(f"\n--- {fname} ---")
    if fpath.exists():
        with open(fpath, 'r') as f:
            data = json.load(f)
        summary = data.get('summary', {})
        print(f"  Total YLL: {summary.get('total_yll', 'N/A'):,.0f}" if isinstance(summary.get('total_yll'), (int, float)) else f"  Total YLL: {summary.get('total_yll', 'N/A')}")
    else:
        print(f"  ❌ File not found")

# Excess mortality
for level in ['intermediate', 'immediate']:
    fname = f'excess_mortality_r_{level}.json'
    fpath = phase1_dir / fname
    print(f"\n--- {fname} ---")
    if fpath.exists():
        with open(fpath, 'r') as f:
            data = json.load(f)
        meta = data.get('metadata', {})
        print(f"  Date range: {meta.get('date_range', 'N/A')}")
        print(f"  Total deaths: {meta.get('total_deaths', 'N/A'):,}" if isinstance(meta.get('total_deaths'), (int, float)) else f"  Total deaths: {meta.get('total_deaths', 'N/A')}")
    else:
        print(f"  ❌ File not found")

# Case-crossover
for level in ['intermediate', 'immediate']:
    fname = f'case_crossover_r_{level}.json'
    fpath = phase1_dir / fname
    print(f"\n--- {fname} ---")
    if fpath.exists():
        with open(fpath, 'r') as f:
            data = json.load(f)
        meta = data.get('metadata', {})
        print(f"  Sample size: {meta.get('n_cases', 'N/A'):,}" if isinstance(meta.get('n_cases'), (int, float)) else f"  Sample size: {meta.get('n_cases', 'N/A')}")
    else:
        print(f"  ❌ File not found")

# =============================================================================
# 3. CHECK PHASE 2 R RESULTS
# =============================================================================
print("\n" + "=" * 80)
print("SECTION 3: PHASE 2 R SCRIPT RESULTS")
print("=" * 80)

phase2_dir = BASE_DIR / "phase2_r" / "results"

for script in ['sensitivity', 'harvesting', 'heatwave']:
    for level in ['intermediate', 'immediate']:
        fname = f'{script}_r_{level}.json'
        fpath = phase2_dir / fname
        print(f"\n--- {fname} ---")
        if fpath.exists():
            with open(fpath, 'r') as f:
                data = json.load(f)
            meta = data.get('metadata', {})
            print(f"  Date range: {meta.get('date_range', 'N/A')}")
            print(f"  Total deaths: {meta.get('total_deaths', 'N/A'):,}" if isinstance(meta.get('total_deaths'), (int, float)) else f"  Total deaths: {meta.get('total_deaths', 'N/A')}")
        else:
            print(f"  ❌ File not found (may still be running)")

# =============================================================================
# 4. CHECK PHASE 3 R RESULTS
# =============================================================================
print("\n" + "=" * 80)
print("SECTION 4: PHASE 3 R SCRIPT RESULTS")
print("=" * 80)

phase3_dir = BASE_DIR / "phase3_r" / "results"

for level in ['intermediate', 'immediate']:
    fname = f'supplementary_r_{level}.json'
    fpath = phase3_dir / fname
    print(f"\n--- {fname} ---")
    if fpath.exists():
        with open(fpath, 'r') as f:
            data = json.load(f)
        meta = data.get('metadata', {})
        print(f"  Date range: {meta.get('date_range', 'N/A')}")
        print(f"  Total deaths: {meta.get('total_deaths', 'N/A'):,}" if isinstance(meta.get('total_deaths'), (int, float)) else f"  Total deaths: {meta.get('total_deaths', 'N/A')}")
    else:
        print(f"  ❌ File not found")

# =============================================================================
# 5. SUMMARY
# =============================================================================
print("\n" + "=" * 80)
print("SUMMARY: DATA COMPLETENESS CHECK")
print("=" * 80)

# Expected deaths for 2021 and 2024 (from the audit we ran earlier)
expected_2021 = 1_249_162
expected_2024 = 1_093_604
expected_total = 13_762_741

print(f"""
Expected values (from current parquet files):
  - 2021 elderly deaths: {expected_2021:,}
  - 2024 elderly deaths: {expected_2024:,}
  - Total elderly deaths: {expected_total:,}

If R script results show significantly different totals, they need to be re-run
with the updated parquet files.

Check the date ranges above:
  - Should be: 2010-01-01 to 2024-12-31
  - If ending before 2024, scripts need re-run
""")

print("=" * 80)
