"""
Audit R Script Results - Version 2
Check actual JSON structure and extract relevant data
"""

import json
import pandas as pd
from pathlib import Path

BASE_DIR = Path(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis")
RESULTS_DIR = BASE_DIR / "phase0_data_prep" / "results"

print("=" * 80)
print("AUDIT: R SCRIPT RESULTS - DETAILED INSPECTION")
print("=" * 80)

# =============================================================================
# 1. INPUT PARQUET FILES
# =============================================================================
print("\n" + "=" * 80)
print("SECTION 1: INPUT PARQUET FILES")
print("=" * 80)

for level in ['intermediate', 'immediate']:
    print(f"\n--- {level.upper()} ---")
    mort_file = RESULTS_DIR / f'mortality_{level}_daily_elderly.parquet'
    if mort_file.exists():
        mort = pd.read_parquet(mort_file)
        mort['date'] = pd.to_datetime(mort['date'])
        mort['year'] = mort['date'].dt.year
        yearly = mort.groupby('year')['deaths_elderly'].sum()
        print(f"  Date range: {mort['date'].min().date()} to {mort['date'].max().date()}")
        print(f"  2021: {yearly.get(2021, 0):,}")
        print(f"  2024: {yearly.get(2024, 0):,}")
        print(f"  Total: {yearly.sum():,}")

# =============================================================================
# 2. PHASE 1 RESULTS - Inspect JSON structure
# =============================================================================
print("\n" + "=" * 80)
print("SECTION 2: PHASE 1 R RESULTS - JSON STRUCTURE INSPECTION")
print("=" * 80)

phase1_dir = BASE_DIR / "phase1_r" / "results"

def inspect_json(fpath, name):
    """Inspect JSON file and extract key info."""
    print(f"\n--- {name} ---")
    if not fpath.exists():
        print("  ❌ File not found")
        return
    
    with open(fpath, 'r') as f:
        data = json.load(f)
    
    # Show top-level keys
    print(f"  Top-level keys: {list(data.keys())}")
    
    # Try to find death counts or date info
    def search_for_key(d, target_keys, prefix=""):
        """Recursively search for keys containing certain strings."""
        results = []
        if isinstance(d, dict):
            for k, v in d.items():
                for target in target_keys:
                    if target.lower() in k.lower():
                        results.append((f"{prefix}{k}", v))
                if isinstance(v, (dict, list)):
                    results.extend(search_for_key(v, target_keys, f"{prefix}{k}."))
        elif isinstance(d, list) and len(d) > 0:
            results.extend(search_for_key(d[0], target_keys, f"{prefix}[0]."))
        return results
    
    # Search for relevant fields
    found = search_for_key(data, ['death', 'total', 'date', 'range', 'observation', 'region', 'n_'])
    
    if found:
        print("  Found relevant fields:")
        for key, val in found[:10]:  # Limit to first 10
            if isinstance(val, (int, float)):
                print(f"    {key}: {val:,.0f}" if val > 100 else f"    {key}: {val}")
            elif isinstance(val, str) and len(val) < 100:
                print(f"    {key}: {val}")
    else:
        # Just show first few items
        print("  Sample content:")
        for k, v in list(data.items())[:3]:
            if isinstance(v, (int, float, str)):
                print(f"    {k}: {v}")
            elif isinstance(v, dict):
                print(f"    {k}: {{...}} ({len(v)} keys)")
            elif isinstance(v, list):
                print(f"    {k}: [...] ({len(v)} items)")

# Check all Phase 1 files
for fname in ['dlnm_r_intermediate_results_v2.json', 'dlnm_r_immediate_results_v2.json',
              'attributable_burden_r_intermediate.json', 'attributable_burden_r_immediate.json',
              'yll_r_intermediate.json', 'yll_r_immediate.json',
              'excess_mortality_r_intermediate.json', 'excess_mortality_r_immediate.json']:
    inspect_json(phase1_dir / fname, fname)

# =============================================================================
# 3. CHECK CSV FILES FOR ACTUAL DATA
# =============================================================================
print("\n" + "=" * 80)
print("SECTION 3: PHASE 1 CSV FILES (may have more detail)")
print("=" * 80)

csv_files = list(phase1_dir.glob("*.csv"))
for csv_file in csv_files:
    print(f"\n--- {csv_file.name} ---")
    try:
        df = pd.read_csv(csv_file)
        print(f"  Shape: {df.shape}")
        print(f"  Columns: {list(df.columns)[:5]}...")
        
        # Check for date columns
        for col in df.columns:
            if 'date' in col.lower() or 'year' in col.lower():
                if df[col].dtype == 'object':
                    vals = df[col].dropna().unique()
                    print(f"  {col} range: {min(vals)} to {max(vals)}")
                else:
                    print(f"  {col} range: {df[col].min()} to {df[col].max()}")
        
        # Check for death columns
        for col in df.columns:
            if 'death' in col.lower() or 'total' in col.lower():
                if pd.api.types.is_numeric_dtype(df[col]):
                    print(f"  {col} sum: {df[col].sum():,.0f}")
    except Exception as e:
        print(f"  Error: {e}")

# =============================================================================
# 4. PHASE 2 & 3 QUICK CHECK
# =============================================================================
print("\n" + "=" * 80)
print("SECTION 4: PHASE 2 & 3 RESULTS")
print("=" * 80)

for phase, phase_dir in [('Phase 2', BASE_DIR / "phase2_r" / "results"),
                          ('Phase 3', BASE_DIR / "phase3_r" / "results")]:
    print(f"\n=== {phase} ===")
    if phase_dir.exists():
        json_files = list(phase_dir.glob("*.json"))
        for jf in json_files:
            inspect_json(jf, jf.name)
    else:
        print("  Directory not found")

print("\n" + "=" * 80)
print("DONE")
print("=" * 80)
