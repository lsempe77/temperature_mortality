"""
Diagnose why certain regions failed to converge in DLNM models.
"""

import pandas as pd
from pathlib import Path

BASE = Path(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis")

# Load data
print("Loading data...")
mort_int = pd.read_parquet(BASE / "phase0_data_prep/results/mortality_intermediate_daily.parquet")
era5_int = pd.read_parquet(BASE / "phase0_data_prep/results/era5_intermediate_daily.parquet")

mort_imm = pd.read_parquet(BASE / "phase0_data_prep/results/mortality_immediate_daily.parquet")
mort_imm = mort_imm.rename(columns={'immediate_code': 'region_code', 'deaths_all': 'deaths'})
mort_int = mort_int.rename(columns={'deaths_all': 'deaths'})

era5_imm = pd.read_parquet(BASE / "phase0_data_prep/results/era5_immediate_daily.parquet")
if 'immediate_code' in era5_imm.columns:
    era5_imm = era5_imm.rename(columns={'immediate_code': 'region_code'})

# Failed regions (as integers to match data)
failed_int = [1201, 4202, 4203, 4205]
failed_imm = [110004, 120002, 310053, 410003, 410004, 410010, 410029, 
              420003, 420005, 420006, 420008, 420014, 420015, 420017, 430039]

print("\n" + "="*70)
print("INTERMEDIATE FAILED REGIONS DIAGNOSTICS")
print("="*70)

for r in failed_int:
    print(f"\n--- Region {r} ---")
    
    # Mortality stats
    mort_sub = mort_int[mort_int['region_code'] == r]
    print(f"Mortality: {len(mort_sub)} days, deaths: {mort_sub['deaths'].min()}-{mort_sub['deaths'].max()} (mean: {mort_sub['deaths'].mean():.1f})")
    
    # Check for zeros
    zero_days = (mort_sub['deaths'] == 0).sum()
    low_days = (mort_sub['deaths'] <= 1).sum()
    print(f"  Zero-death days: {zero_days} ({100*zero_days/len(mort_sub):.1f}%)")
    print(f"  Low-death days (<=1): {low_days} ({100*low_days/len(mort_sub):.1f}%)")
    
    # Temperature stats
    era5_sub = era5_int[era5_int['region_code'] == r]
    print(f"Temperature: range {era5_sub['temp_mean'].min():.1f}-{era5_sub['temp_mean'].max():.1f}°C, std: {era5_sub['temp_mean'].std():.2f}")
    
    # Check temperature variance
    temp_p1 = era5_sub['temp_mean'].quantile(0.01)
    temp_p99 = era5_sub['temp_mean'].quantile(0.99)
    print(f"  P1-P99 range: {temp_p1:.1f} - {temp_p99:.1f}°C (spread: {temp_p99-temp_p1:.1f})")
    
    # Check for NAs
    mort_na = mort_sub['deaths'].isna().sum()
    era5_na = era5_sub['temp_mean'].isna().sum()
    print(f"  Missing: mortality {mort_na}, temperature {era5_na}")

print("\n" + "="*70)
print("IMMEDIATE FAILED REGIONS DIAGNOSTICS")
print("="*70)

for r in failed_imm:
    print(f"\n--- Region {r} ---")
    
    # Mortality stats
    mort_sub = mort_imm[mort_imm['region_code'] == r]
    print(f"Mortality: {len(mort_sub)} days, deaths: {mort_sub['deaths'].min()}-{mort_sub['deaths'].max()} (mean: {mort_sub['deaths'].mean():.1f})")
    
    # Check for zeros
    zero_days = (mort_sub['deaths'] == 0).sum()
    low_days = (mort_sub['deaths'] <= 1).sum()
    print(f"  Zero-death days: {zero_days} ({100*zero_days/len(mort_sub):.1f}%)")
    print(f"  Low-death days (<=1): {low_days} ({100*low_days/len(mort_sub):.1f}%)")
    
    # Temperature stats
    era5_sub = era5_imm[era5_imm['region_code'] == r]
    print(f"Temperature: range {era5_sub['temp_mean'].min():.1f}-{era5_sub['temp_mean'].max():.1f}°C, std: {era5_sub['temp_mean'].std():.2f}")
    
    # Check temperature variance
    temp_p1 = era5_sub['temp_mean'].quantile(0.01)
    temp_p99 = era5_sub['temp_mean'].quantile(0.99)
    print(f"  P1-P99 range: {temp_p1:.1f} - {temp_p99:.1f}°C (spread: {temp_p99-temp_p1:.1f})")
    
    # Check for NAs
    mort_na = mort_sub['deaths'].isna().sum()
    era5_na = era5_sub['temp_mean'].isna().sum()
    print(f"  Missing: mortality {mort_na}, temperature {era5_na}")

# Compare to successfully converged regions
print("\n" + "="*70)
print("COMPARISON: Converged vs Non-Converged (Intermediate)")
print("="*70)

# Get all region stats
int_stats = []
for r in mort_int['region_code'].unique():
    mort_sub = mort_int[mort_int['region_code'] == r]
    era5_sub = era5_int[era5_int['region_code'] == r]
    
    int_stats.append({
        'region_code': r,
        'n_days': len(mort_sub),
        'mean_deaths': mort_sub['deaths'].mean(),
        'zero_pct': 100 * (mort_sub['deaths'] == 0).sum() / len(mort_sub) if len(mort_sub) > 0 else 0,
        'temp_std': era5_sub['temp_mean'].std() if len(era5_sub) > 0 else 0,
        'temp_range': era5_sub['temp_mean'].max() - era5_sub['temp_mean'].min() if len(era5_sub) > 0 else 0,
        'failed': r in failed_int
    })

int_stats_df = pd.DataFrame(int_stats)

print("\nMean stats for converged regions:")
conv = int_stats_df[~int_stats_df['failed']]
print(f"  Mean deaths/day: {conv['mean_deaths'].mean():.1f}")
print(f"  Zero-death days: {conv['zero_pct'].mean():.1f}%")
print(f"  Temp std: {conv['temp_std'].mean():.2f}")
print(f"  Temp range: {conv['temp_range'].mean():.1f}°C")

print("\nMean stats for NON-converged regions:")
fail = int_stats_df[int_stats_df['failed']]
print(f"  Mean deaths/day: {fail['mean_deaths'].mean():.1f}")
print(f"  Zero-death days: {fail['zero_pct'].mean():.1f}%")
print(f"  Temp std: {fail['temp_std'].mean():.2f}")
print(f"  Temp range: {fail['temp_range'].mean():.1f}°C")

print("\n" + "="*70)
print("Diagnosis complete!")
print("="*70)
