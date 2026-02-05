"""
Check Phase 4 stratified data completeness
"""
import pandas as pd
from pathlib import Path

RESULTS_DIR = Path(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis\phase0_data_prep\results")

print("="*80)
print("PHASE 4 STRATIFIED DATA CHECK")
print("="*80)

# Age stratified
print("\n=== AGE STRATIFIED ===")
for level in ['intermediate', 'immediate']:
    print(f"\n--- {level.upper()} ---")
    for age in ['age_60_69', 'age_70_79', 'age_80plus']:
        f = RESULTS_DIR / f'mortality_{level}_daily_{age}.parquet'
        if f.exists():
            df = pd.read_parquet(f)
            df['date'] = pd.to_datetime(df['date'])
            df['year'] = df['date'].dt.year
            deaths_col = [c for c in df.columns if 'deaths' in c.lower()][0]
            yearly = df.groupby('year')[deaths_col].sum()
            status = "✅" if yearly.get(2021,0) > 0 and yearly.get(2024,0) > 0 else "❌ MISSING DATA"
            print(f"  {age}: Total={yearly.sum():,}, 2021={yearly.get(2021,0):,}, 2024={yearly.get(2024,0):,} {status}")
        else:
            print(f"  {age}: ❌ FILE NOT FOUND")

# Sex stratified
print("\n=== SEX STRATIFIED ===")
for level in ['intermediate', 'immediate']:
    print(f"\n--- {level.upper()} ---")
    for sex in ['male', 'female']:
        f = RESULTS_DIR / f'mortality_{level}_daily_{sex}.parquet'
        if f.exists():
            df = pd.read_parquet(f)
            df['date'] = pd.to_datetime(df['date'])
            df['year'] = df['date'].dt.year
            deaths_col = [c for c in df.columns if 'deaths' in c.lower()][0]
            yearly = df.groupby('year')[deaths_col].sum()
            status = "✅" if yearly.get(2021,0) > 0 and yearly.get(2024,0) > 0 else "❌ MISSING DATA"
            print(f"  {sex}: Total={yearly.sum():,}, 2021={yearly.get(2021,0):,}, 2024={yearly.get(2024,0):,} {status}")
        else:
            print(f"  {sex}: ❌ FILE NOT FOUND")

# Cause stratified
print("\n=== CAUSE STRATIFIED ===")
for level in ['intermediate', 'immediate']:
    print(f"\n--- {level.upper()} ---")
    for cause in ['cvd', 'respiratory', 'external', 'other']:
        f = RESULTS_DIR / f'mortality_{level}_daily_{cause}.parquet'
        if f.exists():
            df = pd.read_parquet(f)
            df['date'] = pd.to_datetime(df['date'])
            df['year'] = df['date'].dt.year
            deaths_col = [c for c in df.columns if 'deaths' in c.lower()][0]
            yearly = df.groupby('year')[deaths_col].sum()
            status = "✅" if yearly.get(2021,0) > 0 and yearly.get(2024,0) > 0 else "❌ MISSING DATA"
            print(f"  {cause}: Total={yearly.sum():,}, 2021={yearly.get(2021,0):,}, 2024={yearly.get(2024,0):,} {status}")
        else:
            print(f"  {cause}: ❌ FILE NOT FOUND")

print("\n" + "="*80)
