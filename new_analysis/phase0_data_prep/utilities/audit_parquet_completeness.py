"""
Audit all mortality parquet files to check data completeness.
"""
import pandas as pd
from pathlib import Path

RESULTS_DIR = Path(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis\phase0_data_prep\results")

parquet_files = [
    # Main files
    "mortality_intermediate_daily.parquet",
    "mortality_intermediate_daily_elderly.parquet",
    "mortality_immediate_daily.parquet",
    "mortality_immediate_daily_elderly.parquet",
    # Age stratified
    "mortality_intermediate_daily_age_60_69.parquet",
    "mortality_intermediate_daily_age_70_79.parquet",
    "mortality_intermediate_daily_age_80plus.parquet",
    "mortality_immediate_daily_age_60_69.parquet",
    "mortality_immediate_daily_age_70_79.parquet",
    "mortality_immediate_daily_age_80plus.parquet",
    # Sex stratified
    "mortality_intermediate_daily_male.parquet",
    "mortality_intermediate_daily_female.parquet",
    "mortality_immediate_daily_male.parquet",
    "mortality_immediate_daily_female.parquet",
    # Cause stratified
    "mortality_intermediate_daily_cvd.parquet",
    "mortality_intermediate_daily_respiratory.parquet",
    "mortality_intermediate_daily_external.parquet",
    "mortality_immediate_daily_cvd.parquet",
    "mortality_immediate_daily_respiratory.parquet",
    "mortality_immediate_daily_external.parquet",
]

print("="*100)
print("MORTALITY PARQUET COMPLETENESS AUDIT")
print("="*100)
print(f"{'File':<50} {'Total':>12} {'2021':>12} {'2024':>12} {'Missing?'}")
print("-"*100)

for fname in parquet_files:
    fpath = RESULTS_DIR / fname
    if not fpath.exists():
        print(f"{fname:<50} {'NOT FOUND':>12}")
        continue
    
    try:
        df = pd.read_parquet(fpath)
        
        # Find date column
        date_col = [c for c in df.columns if 'date' in c.lower()][0]
        df[date_col] = pd.to_datetime(df[date_col])
        df['year'] = df[date_col].dt.year
        
        # Find deaths column
        deaths_col = [c for c in df.columns if 'deaths' in c.lower()][0]
        
        # Get totals
        yearly = df.groupby('year')[deaths_col].sum()
        total = yearly.sum()
        d2021 = yearly.get(2021, 0)
        d2024 = yearly.get(2024, 0)
        
        missing = []
        if d2021 == 0:
            missing.append("2021")
        if d2024 == 0:
            missing.append("2024")
        
        missing_str = ", ".join(missing) if missing else "OK"
        
        print(f"{fname:<50} {total:>12,} {d2021:>12,} {d2024:>12,} {missing_str}")
        
    except Exception as e:
        print(f"{fname:<50} ERROR: {e}")

print("="*100)
