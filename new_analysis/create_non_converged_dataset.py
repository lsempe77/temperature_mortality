"""
Create comprehensive dataset of non-converged regions for re-analysis.
Includes diagnostic information to help choose alternative model specifications.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
from datetime import datetime

BASE = Path(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis")

print("="*70)
print("Creating Comprehensive Non-Converged Regions Dataset")
print("="*70)

# Load source data
print("\nLoading source data...")
mort_int = pd.read_parquet(BASE / "phase0_data_prep/results/mortality_intermediate_daily.parquet")
era5_int = pd.read_parquet(BASE / "phase0_data_prep/results/era5_intermediate_daily.parquet")
mort_imm = pd.read_parquet(BASE / "phase0_data_prep/results/mortality_immediate_daily.parquet")
era5_imm = pd.read_parquet(BASE / "phase0_data_prep/results/era5_immediate_daily.parquet")

# Rename columns for consistency
mort_imm = mort_imm.rename(columns={'immediate_code': 'region_code', 'deaths_all': 'deaths'})
mort_int = mort_int.rename(columns={'deaths_all': 'deaths'})
if 'immediate_code' in era5_imm.columns:
    era5_imm = era5_imm.rename(columns={'immediate_code': 'region_code'})

# Define failed regions from DLNM first stage
failed_int = [1201, 4202, 4203, 4205]
failed_imm = [110004, 120002, 310053, 410003, 410004, 410010, 410029, 
              420003, 420005, 420006, 420008, 420014, 420015, 420017, 430039]

# State lookup
state_codes = {
    '11': 'Rondonia', '12': 'Acre', '13': 'Amazonas', '14': 'Roraima', '15': 'Para',
    '16': 'Amapa', '17': 'Tocantins', '21': 'Maranhao', '22': 'Piaui', '23': 'Ceara',
    '24': 'Rio Grande do Norte', '25': 'Paraiba', '26': 'Pernambuco', '27': 'Alagoas',
    '28': 'Sergipe', '29': 'Bahia', '31': 'Minas Gerais', '32': 'Espirito Santo',
    '33': 'Rio de Janeiro', '35': 'Sao Paulo', '41': 'Parana', '42': 'Santa Catarina',
    '43': 'Rio Grande do Sul', '50': 'Mato Grosso do Sul', '51': 'Mato Grosso',
    '52': 'Goias', '53': 'Distrito Federal'
}

# Build comprehensive dataset
print("\nBuilding diagnostic dataset for non-converged regions...")

results = []

# Process intermediate regions
for r in failed_int:
    mort_sub = mort_int[mort_int['region_code'] == r]
    era5_sub = era5_int[era5_int['region_code'] == r]
    
    state = str(r)[:2]
    
    # Calculate diagnostic metrics
    n_days = len(mort_sub)
    total_deaths = mort_sub['deaths'].sum()
    mean_deaths = mort_sub['deaths'].mean()
    max_deaths = mort_sub['deaths'].max()
    min_deaths = mort_sub['deaths'].min()
    zero_days = (mort_sub['deaths'] == 0).sum()
    low_death_days = (mort_sub['deaths'] <= 1).sum()
    low_death_pct = 100 * low_death_days / n_days if n_days > 0 else 0
    
    temp_mean = era5_sub['temp_mean'].mean()
    temp_std = era5_sub['temp_mean'].std()
    temp_min = era5_sub['temp_mean'].min()
    temp_max = era5_sub['temp_mean'].max()
    temp_p1 = era5_sub['temp_mean'].quantile(0.01)
    temp_p99 = era5_sub['temp_mean'].quantile(0.99)
    temp_spread = temp_p99 - temp_p1
    
    # Determine likely cause of non-convergence
    cause = []
    if mean_deaths < 5:
        cause.append("low_deaths")
    if low_death_pct > 10:
        cause.append("sparse_data")
    if temp_spread < 10:
        cause.append("low_temp_variance")
    if not cause:
        cause.append("unknown")
    
    # Suggest alternative specifications
    suggestions = []
    if "low_deaths" in cause or "sparse_data" in cause:
        suggestions.append("aggregate_to_larger_region")
        suggestions.append("use_weekly_instead_of_daily")
        suggestions.append("reduce_spline_df")
    if "low_temp_variance" in cause:
        suggestions.append("use_linear_model")
        suggestions.append("reduce_temp_df")
    
    results.append({
        'region_code': r,
        'exposure_type': 'intermediate',
        'state_code': state,
        'state_name': state_codes.get(state, 'Unknown'),
        'n_days': n_days,
        'total_deaths': total_deaths,
        'mean_daily_deaths': mean_deaths,
        'max_daily_deaths': max_deaths,
        'min_daily_deaths': min_deaths,
        'zero_death_days': zero_days,
        'low_death_days_pct': low_death_pct,
        'temp_mean': temp_mean,
        'temp_std': temp_std,
        'temp_min': temp_min,
        'temp_max': temp_max,
        'temp_p1': temp_p1,
        'temp_p99': temp_p99,
        'temp_p1_p99_spread': temp_spread,
        'likely_cause': '; '.join(cause),
        'suggested_alternatives': '; '.join(suggestions)
    })

# Process immediate regions
for r in failed_imm:
    mort_sub = mort_imm[mort_imm['region_code'] == r]
    era5_sub = era5_imm[era5_imm['region_code'] == r]
    
    state = str(r)[:2]
    
    n_days = len(mort_sub)
    total_deaths = mort_sub['deaths'].sum()
    mean_deaths = mort_sub['deaths'].mean()
    max_deaths = mort_sub['deaths'].max()
    min_deaths = mort_sub['deaths'].min()
    zero_days = (mort_sub['deaths'] == 0).sum()
    low_death_days = (mort_sub['deaths'] <= 1).sum()
    low_death_pct = 100 * low_death_days / n_days if n_days > 0 else 0
    
    temp_mean = era5_sub['temp_mean'].mean()
    temp_std = era5_sub['temp_mean'].std()
    temp_min = era5_sub['temp_mean'].min()
    temp_max = era5_sub['temp_mean'].max()
    temp_p1 = era5_sub['temp_mean'].quantile(0.01)
    temp_p99 = era5_sub['temp_mean'].quantile(0.99)
    temp_spread = temp_p99 - temp_p1
    
    # Determine likely cause
    cause = []
    if mean_deaths < 5:
        cause.append("low_deaths")
    if low_death_pct > 10:
        cause.append("sparse_data")
    if temp_spread < 10:
        cause.append("low_temp_variance")
    if not cause:
        cause.append("unknown")
    
    # Suggest alternatives
    suggestions = []
    if "low_deaths" in cause or "sparse_data" in cause:
        suggestions.append("aggregate_to_larger_region")
        suggestions.append("use_weekly_instead_of_daily")
        suggestions.append("reduce_spline_df")
    if "low_temp_variance" in cause:
        suggestions.append("use_linear_model")
        suggestions.append("reduce_temp_df")
    
    results.append({
        'region_code': r,
        'exposure_type': 'immediate',
        'state_code': state,
        'state_name': state_codes.get(state, 'Unknown'),
        'n_days': n_days,
        'total_deaths': total_deaths,
        'mean_daily_deaths': mean_deaths,
        'max_daily_deaths': max_deaths,
        'min_daily_deaths': min_deaths,
        'zero_death_days': zero_days,
        'low_death_days_pct': low_death_pct,
        'temp_mean': temp_mean,
        'temp_std': temp_std,
        'temp_min': temp_min,
        'temp_max': temp_max,
        'temp_p1': temp_p1,
        'temp_p99': temp_p99,
        'temp_p1_p99_spread': temp_spread,
        'likely_cause': '; '.join(cause),
        'suggested_alternatives': '; '.join(suggestions)
    })

# Create DataFrame
df = pd.DataFrame(results)
df = df.sort_values(['exposure_type', 'region_code'])

# Save to CSV
output_path = BASE / "results" / "non_converged_regions_for_reanalysis.csv"
df.to_csv(output_path, index=False)
print(f"\nSaved comprehensive dataset to: {output_path}")

# Print summary
print("\n" + "="*70)
print("SUMMARY OF NON-CONVERGED REGIONS")
print("="*70)

print(f"\nTotal non-converged regions: {len(df)}")
print(f"  Intermediate: {len(df[df['exposure_type']=='intermediate'])}")
print(f"  Immediate: {len(df[df['exposure_type']=='immediate'])}")

print("\nBy State:")
for state, count in df['state_name'].value_counts().items():
    print(f"  {state}: {count}")

print("\nBy Likely Cause:")
cause_counts = {}
for causes in df['likely_cause']:
    for c in causes.split('; '):
        cause_counts[c] = cause_counts.get(c, 0) + 1
for cause, count in sorted(cause_counts.items(), key=lambda x: -x[1]):
    print(f"  {cause}: {count}")

print("\n" + "="*70)
print("DETAILED NON-CONVERGED REGIONS")
print("="*70)

print("\n--- INTERMEDIATE REGIONS (4 regions) ---")
int_df = df[df['exposure_type']=='intermediate']
for _, row in int_df.iterrows():
    print(f"\nRegion {row['region_code']} ({row['state_name']})")
    print(f"  Deaths: {row['mean_daily_deaths']:.1f}/day (range: {row['min_daily_deaths']}-{row['max_daily_deaths']})")
    print(f"  Low-death days: {row['low_death_days_pct']:.1f}%")
    print(f"  Temp spread (P1-P99): {row['temp_p1_p99_spread']:.1f}°C")
    print(f"  Likely cause: {row['likely_cause']}")
    print(f"  Suggested: {row['suggested_alternatives']}")

print("\n--- IMMEDIATE REGIONS (15 regions) ---")
imm_df = df[df['exposure_type']=='immediate']
for _, row in imm_df.iterrows():
    print(f"\nRegion {row['region_code']} ({row['state_name']})")
    print(f"  Deaths: {row['mean_daily_deaths']:.1f}/day (range: {row['min_daily_deaths']}-{row['max_daily_deaths']})")
    print(f"  Low-death days: {row['low_death_days_pct']:.1f}%")
    print(f"  Temp spread (P1-P99): {row['temp_p1_p99_spread']:.1f}°C")
    print(f"  Likely cause: {row['likely_cause']}")
    print(f"  Suggested: {row['suggested_alternatives']}")

# Also create a summary JSON for programmatic use
summary = {
    "generated": datetime.now().isoformat(),
    "total_non_converged": len(df),
    "intermediate": {
        "count": len(df[df['exposure_type']=='intermediate']),
        "region_codes": df[df['exposure_type']=='intermediate']['region_code'].tolist()
    },
    "immediate": {
        "count": len(df[df['exposure_type']=='immediate']),
        "region_codes": df[df['exposure_type']=='immediate']['region_code'].tolist()
    },
    "causes": cause_counts,
    "states_affected": df['state_name'].value_counts().to_dict(),
    "recommendations": {
        "low_deaths": "Consider aggregating to larger regions or using weekly data",
        "sparse_data": "Reduce spline degrees of freedom or use simpler model",
        "low_temp_variance": "Use linear temperature term instead of non-linear spline"
    }
}

summary_path = BASE / "results" / "non_converged_regions_summary.json"
with open(summary_path, 'w') as f:
    json.dump(summary, f, indent=2)
print(f"\nSaved summary JSON to: {summary_path}")

print("\n" + "="*70)
print("Script complete!")
print("="*70)
