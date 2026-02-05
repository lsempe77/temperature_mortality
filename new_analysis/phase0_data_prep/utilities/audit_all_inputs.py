"""
COMPREHENSIVE DATA AUDIT
Audit ALL input data sources used by R scripts
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

BASE_DIR = Path(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis")
RESULTS_DIR = BASE_DIR / "phase0_data_prep" / "results"

print("=" * 90)
print("COMPREHENSIVE DATA AUDIT - ALL INPUT SOURCES")
print(f"Run at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 90)

# =============================================================================
# 1. MORTALITY DATA
# =============================================================================
print("\n" + "=" * 90)
print("1. MORTALITY DATA")
print("=" * 90)

for level in ['intermediate', 'immediate']:
    print(f"\n--- {level.upper()} ---")
    
    # Main elderly mortality
    mort_file = RESULTS_DIR / f'mortality_{level}_daily_elderly.parquet'
    if mort_file.exists():
        mort = pd.read_parquet(mort_file)
        mort['date'] = pd.to_datetime(mort['date'])
        mort['year'] = mort['date'].dt.year
        
        region_col = [c for c in mort.columns if 'code' in c.lower()][0]
        deaths_col = [c for c in mort.columns if 'deaths' in c.lower()][0]
        
        yearly = mort.groupby('year')[deaths_col].sum()
        
        print(f"  File: {mort_file.name}")
        print(f"  Date range: {mort['date'].min().date()} to {mort['date'].max().date()}")
        print(f"  Regions: {mort[region_col].nunique()}")
        print(f"  Total observations: {len(mort):,}")
        print(f"  Total deaths: {yearly.sum():,}")
        print(f"  Deaths by year:")
        for y in [2010, 2011, 2020, 2021, 2022, 2023, 2024]:
            d = yearly.get(y, 0)
            flag = "⚠️ ZERO!" if d == 0 else ""
            print(f"    {y}: {d:>10,} {flag}")
    else:
        print(f"  ❌ File not found: {mort_file.name}")

# =============================================================================
# 2. ERA5 TEMPERATURE DATA
# =============================================================================
print("\n" + "=" * 90)
print("2. ERA5 TEMPERATURE DATA")
print("=" * 90)

for level in ['intermediate', 'immediate']:
    print(f"\n--- {level.upper()} ---")
    
    era5_file = RESULTS_DIR / f'era5_{level}_daily.parquet'
    if era5_file.exists():
        era5 = pd.read_parquet(era5_file)
        era5['date'] = pd.to_datetime(era5['date'])
        era5['year'] = era5['date'].dt.year
        
        print(f"  File: {era5_file.name}")
        print(f"  Date range: {era5['date'].min().date()} to {era5['date'].max().date()}")
        print(f"  Regions: {era5['region_code'].nunique()}")
        print(f"  Total observations: {len(era5):,}")
        print(f"  Columns: {list(era5.columns)}")
        
        # Check for temp_mean column
        if 'temp_mean' in era5.columns:
            print(f"  temp_mean range: {era5['temp_mean'].min():.1f}°C to {era5['temp_mean'].max():.1f}°C")
            print(f"  temp_mean NaN: {era5['temp_mean'].isna().sum():,}")
        
        # Check coverage by year
        yearly_obs = era5.groupby('year').size()
        print(f"  Observations by year:")
        for y in [2010, 2021, 2024]:
            obs = yearly_obs.get(y, 0)
            flag = "⚠️ ZERO!" if obs == 0 else ""
            print(f"    {y}: {obs:>10,} {flag}")
    else:
        print(f"  ❌ File not found: {era5_file.name}")

# =============================================================================
# 3. CAMS POLLUTION DATA
# =============================================================================
print("\n" + "=" * 90)
print("3. CAMS POLLUTION DATA")
print("=" * 90)

for level in ['intermediate', 'immediate']:
    print(f"\n--- {level.upper()} ---")
    
    cams_file = RESULTS_DIR / f'cams_{level}_daily.parquet'
    if cams_file.exists():
        cams = pd.read_parquet(cams_file)
        cams['date'] = pd.to_datetime(cams['date'])
        cams['year'] = cams['date'].dt.year
        
        print(f"  File: {cams_file.name}")
        print(f"  Date range: {cams['date'].min().date()} to {cams['date'].max().date()}")
        
        # Find region column (could be region_code, intermediate_code, or immediate_code)
        region_cols = [c for c in cams.columns if 'code' in c.lower()]
        if region_cols:
            print(f"  Regions: {cams[region_cols[0]].nunique()}")
        
        print(f"  Total observations: {len(cams):,}")
        print(f"  Columns: {list(cams.columns)}")
        
        # Check for PM2.5 and O3
        for col in ['pm25_mean', 'o3_mean', 'pm2p5', 'o3', 'pm25', 'ozone', 'pm10']:
            if col in cams.columns:
                non_null = cams[col].notna().sum()
                print(f"  {col}: {non_null:,} non-null values ({non_null/len(cams)*100:.1f}%)")
        
        yearly_obs = cams.groupby('year').size()
        print(f"  Observations by year:")
        for y in [2010, 2021, 2024]:
            obs = yearly_obs.get(y, 0)
            flag = "⚠️ ZERO!" if obs == 0 else ""
            print(f"    {y}: {obs:>10,} {flag}")
    else:
        print(f"  ❌ File not found: {cams_file.name}")

# =============================================================================
# 4. INFLUENZA DATA
# =============================================================================
print("\n" + "=" * 90)
print("4. INFLUENZA DATA")
print("=" * 90)

for level in ['intermediate', 'immediate']:
    print(f"\n--- {level.upper()} ---")
    
    flu_file = RESULTS_DIR / f'influenza_daily_by_{level}_region.parquet'
    if flu_file.exists():
        flu = pd.read_parquet(flu_file)
        
        # Find date column
        date_cols = [c for c in flu.columns if 'date' in c.lower()]
        if date_cols:
            flu['date'] = pd.to_datetime(flu[date_cols[0]])
            flu['year'] = flu['date'].dt.year
            
            print(f"  File: {flu_file.name}")
            print(f"  Date range: {flu['date'].min().date()} to {flu['date'].max().date()}")
            print(f"  Total observations: {len(flu):,}")
            print(f"  Columns: {list(flu.columns)}")
            
            yearly_obs = flu.groupby('year').size()
            print(f"  Observations by year:")
            for y in [2010, 2019, 2020, 2021, 2024]:
                obs = yearly_obs.get(y, 0)
                flag = "⚠️ ZERO!" if obs == 0 else ""
                print(f"    {y}: {obs:>10,} {flag}")
        else:
            print(f"  Columns: {list(flu.columns)}")
            print(f"  ⚠️ No date column found")
    else:
        print(f"  ❌ File not found: {flu_file.name}")

# =============================================================================
# 5. REGION MAPPING
# =============================================================================
print("\n" + "=" * 90)
print("5. REGION MAPPING")
print("=" * 90)

mapping_file = RESULTS_DIR / 'municipality_to_all_regions_map.csv'
if mapping_file.exists():
    mapping = pd.read_csv(mapping_file)
    print(f"  File: {mapping_file.name}")
    print(f"  Total municipalities: {len(mapping):,}")
    print(f"  Columns: {list(mapping.columns)}")
    print(f"  Intermediate regions: {mapping['intermediate_code'].nunique()}")
    print(f"  Immediate regions: {mapping['immediate_code'].nunique()}")
    
    # Check for missing mappings
    int_missing = mapping['intermediate_code'].isna().sum()
    imm_missing = mapping['immediate_code'].isna().sum()
    print(f"  Municipalities without intermediate mapping: {int_missing}")
    print(f"  Municipalities without immediate mapping: {imm_missing}")
else:
    print(f"  ❌ File not found: {mapping_file.name}")

# =============================================================================
# 6. IMMEDIATE vs INTERMEDIATE DEATH DISCREPANCY
# =============================================================================
print("\n" + "=" * 90)
print("6. INVESTIGATING IMMEDIATE vs INTERMEDIATE DEATH DISCREPANCY")
print("=" * 90)

mort_int = pd.read_parquet(RESULTS_DIR / 'mortality_intermediate_daily_elderly.parquet')
mort_imm = pd.read_parquet(RESULTS_DIR / 'mortality_immediate_daily_elderly.parquet')

total_int = mort_int['deaths_elderly'].sum()
total_imm = mort_imm['deaths_elderly'].sum()

print(f"  Intermediate total: {total_int:,}")
print(f"  Immediate total: {total_imm:,}")
print(f"  Difference: {total_int - total_imm:,}")

if total_int != total_imm:
    print(f"\n  ⚠️ DISCREPANCY DETECTED!")
    print(f"  The immediate level is missing {total_int - total_imm:,} deaths")
    print(f"  This is {(total_int - total_imm) / total_int * 100:.2f}% of total")
    
    # Compare by year
    mort_int['date'] = pd.to_datetime(mort_int['date'])
    mort_imm['date'] = pd.to_datetime(mort_imm['date'])
    mort_int['year'] = mort_int['date'].dt.year
    mort_imm['year'] = mort_imm['date'].dt.year
    
    yearly_int = mort_int.groupby('year')['deaths_elderly'].sum()
    yearly_imm = mort_imm.groupby('year')['deaths_elderly'].sum()
    
    print(f"\n  Year-by-year comparison:")
    print(f"  {'Year':<6} {'Intermediate':>12} {'Immediate':>12} {'Diff':>10}")
    print(f"  {'-'*42}")
    for y in range(2010, 2025):
        i = yearly_int.get(y, 0)
        m = yearly_imm.get(y, 0)
        d = i - m
        flag = "⚠️" if d != 0 else ""
        print(f"  {y:<6} {i:>12,} {m:>12,} {d:>10,} {flag}")
else:
    print(f"\n  ✅ Both levels have identical death counts")

# =============================================================================
# 7. SES COVARIATES
# =============================================================================
print("\n" + "=" * 90)
print("7. SES COVARIATES")
print("=" * 90)

for level in ['intermediate', 'immediate']:
    print(f"\n--- {level.upper()} ---")
    
    ses_file = RESULTS_DIR / f'ses_{level}_covariates.csv'
    if ses_file.exists():
        ses = pd.read_csv(ses_file)
        print(f"  File: {ses_file.name}")
        print(f"  Regions: {len(ses)}")
        print(f"  Columns: {list(ses.columns)}")
        
        # Check for NaN values
        nan_counts = ses.isna().sum()
        if nan_counts.sum() > 0:
            print(f"  Columns with NaN:")
            for col in nan_counts[nan_counts > 0].index:
                print(f"    {col}: {nan_counts[col]} NaN values")
    else:
        print(f"  ❌ File not found: {ses_file.name}")

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 90)
print("AUDIT SUMMARY")
print("=" * 90)

print("""
KEY FINDINGS:
-------------
1. Check date ranges - all should cover 2010-2024
2. Check 2021 and 2024 specifically (known format issues)
3. Check intermediate vs immediate discrepancy
4. Check NaN values in each dataset

NEXT STEPS:
-----------
- If mortality discrepancy exists: check mapping coverage
- If ERA5/CAMS missing years: check download scripts
- If influenza incomplete: may need to re-download
- After unified processor finishes: re-run audit
""")

print("=" * 90)
