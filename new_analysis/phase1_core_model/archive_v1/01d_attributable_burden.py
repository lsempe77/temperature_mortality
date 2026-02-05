"""
01d: ATTRIBUTABLE BURDEN CALCULATION
=====================================
Calculates temperature-attributable mortality using DLNM results from 01a.

Method (Gasparrini & Leone, 2014):
1. Load fitted DLNM models from 01a (intermediate and immediate levels)
2. For each region, reconstruct the E-R curve using stored coefficients
3. Calculate attributable fraction: AF = (RR - 1) / RR for each day
4. Sum attributable deaths:
   - Heat burden: AF × deaths for days > P75 (or where temp > MMT)
   - Cold burden: AF × deaths for days < P25 (or where temp < MMT)
5. Aggregate to national level with uncertainty

Key Outputs:
- Attributable fraction (AF) for heat and cold
- Attributable number (AN) of deaths
- Attributable deaths per 100,000 population
- Annual attributable mortality

References:
- Gasparrini & Leone (2014) - Attributable risk from DLNM
- Gasparrini et al. (2015) Lancet - Temperature-mortality burden
"""

import warnings
warnings.filterwarnings('ignore')

import os
import json
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
from patsy import dmatrix

# =============================================================================
# PATHS
# =============================================================================

SCRIPT_DIR = Path(__file__).parent
PHASE0_RESULTS = SCRIPT_DIR.parent / 'phase0_data_prep' / 'results'
PHASE1_RESULTS = SCRIPT_DIR / 'results'
OUTPUT_DIR = PHASE1_RESULTS
OUTPUT_DIR.mkdir(exist_ok=True)

print("=" * 70)
print("01d: ATTRIBUTABLE BURDEN CALCULATION")
print("=" * 70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print()

# =============================================================================
# LOAD DATA - INTERMEDIATE LEVEL (133 regions)
# =============================================================================
print("[1] Loading data for INTERMEDIATE level (133 regions)...")

# Load temperature data
df_temp_int = pd.read_parquet(PHASE0_RESULTS / 'era5_intermediate_daily.parquet')
df_temp_int['date'] = pd.to_datetime(df_temp_int['date'])
print(f"  ERA5 Intermediate: {len(df_temp_int):,} rows")

# Load mortality data  
df_mort_int = pd.read_parquet(PHASE0_RESULTS / 'mortality_regional_daily_elderly.parquet')
df_mort_int['date'] = pd.to_datetime(df_mort_int['date'])
print(f"  Mortality Intermediate: {len(df_mort_int):,} rows")

# Load SES for population
df_ses_int = pd.read_csv(PHASE0_RESULTS / 'ses_intermediate_covariates.csv')
pop_map_int = dict(zip(df_ses_int['intermediate_code'].astype(str), df_ses_int['pop_elderly']))
print(f"  SES Intermediate: {len(df_ses_int)} regions")

# Merge intermediate data
df_int = pd.merge(
    df_mort_int,
    df_temp_int[['date', 'region_code', 'temp_mean']],
    on=['date', 'region_code'],
    how='inner'
)
df_int['region_code'] = df_int['region_code'].astype(str)
df_int['pop_elderly'] = df_int['region_code'].map(pop_map_int)
print(f"  Merged Intermediate: {len(df_int):,} rows, {df_int['region_code'].nunique()} regions")

# =============================================================================
# LOAD DATA - IMMEDIATE LEVEL (510 regions)
# =============================================================================
print("\n[1b] Loading data for IMMEDIATE level (510 regions)...")

# Load temperature data
df_temp_imm = pd.read_parquet(PHASE0_RESULTS / 'era5_immediate_daily.parquet')
df_temp_imm['date'] = pd.to_datetime(df_temp_imm['date'])
# Rename region_code to immediate_code for consistency
df_temp_imm = df_temp_imm.rename(columns={'region_code': 'immediate_code'})
print(f"  ERA5 Immediate: {len(df_temp_imm):,} rows")

# Load mortality data
df_mort_imm = pd.read_parquet(PHASE0_RESULTS / 'mortality_immediate_daily_elderly.parquet')
df_mort_imm['date'] = pd.to_datetime(df_mort_imm['date'])
print(f"  Mortality Immediate: {len(df_mort_imm):,} rows")

# Load SES for population
df_ses_imm = pd.read_csv(PHASE0_RESULTS / 'ses_immediate_covariates.csv')
pop_map_imm = dict(zip(df_ses_imm['immediate_code'].astype(str), df_ses_imm['pop_elderly']))
print(f"  SES Immediate: {len(df_ses_imm)} regions")

# Merge immediate data
df_imm = pd.merge(
    df_mort_imm,
    df_temp_imm[['date', 'immediate_code', 'temp_mean']],
    on=['date', 'immediate_code'],
    how='inner'
)
df_imm['immediate_code'] = df_imm['immediate_code'].astype(str)
df_imm['pop_elderly'] = df_imm['immediate_code'].map(pop_map_imm)
print(f"  Merged Immediate: {len(df_imm):,} rows, {df_imm['immediate_code'].nunique()} regions")

print(f"\n  Date range: {df_int['date'].min().date()} to {df_int['date'].max().date()}")

# Load DLNM results
print("\n[2] Loading DLNM results from 01a...")

with open(PHASE1_RESULTS / 'dlnm_intermediate_results.json', 'r') as f:
    intermediate_results = json.load(f)
    
with open(PHASE1_RESULTS / 'dlnm_immediate_results.json', 'r') as f:
    immediate_results = json.load(f)

print(f"  Intermediate regions: {intermediate_results['n_successful']} / {intermediate_results['n_regions']}")
print(f"  Immediate regions: {immediate_results['n_successful']} / {immediate_results['n_regions']}")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def compute_cumulative_rr(target_temp, ref_temp, cb_coefs, col_names, temp_center, temp_scale):
    """
    Compute cumulative RR for a temperature relative to reference.
    Uses the polynomial cross-basis structure from 01a.
    
    The cross-basis has columns: lag{0-21}_poly{1-3}
    - For each lag, we have 3 polynomial terms (temp^1, temp^2, temp^3)
    - Cumulative effect sums the contrast across all lags
    """
    temp_std_target = (target_temp - temp_center) / temp_scale
    temp_std_ref = (ref_temp - temp_center) / temp_scale
    
    # Create contrast vector
    contrast = np.zeros(len(col_names))
    
    for i, name in enumerate(col_names):
        # Parse lag and polynomial from column name: lag{L}_poly{P}
        parts = name.split('_')
        poly = int(parts[1].replace('poly', ''))
        
        # Difference in polynomial terms
        contrast[i] = temp_std_target**poly - temp_std_ref**poly
    
    # Cumulative effect (sum over all lags)
    log_rr = np.dot(contrast, cb_coefs)
    
    return np.exp(log_rr)


def compute_attributable_fraction(rr):
    """Calculate attributable fraction from relative risk."""
    if rr <= 0:
        return 0
    af = (rr - 1) / rr
    return np.clip(af, -1, 1)


def calculate_burden_for_region(df_region, region_results, region_code):
    """
    Calculate attributable burden for a single region using multiple thresholds.
    
    Uses P2.5/P97.5 as primary thresholds and P1/P99 as secondary (comparison).
    
    Parameters:
    -----------
    df_region : pd.DataFrame
        Daily mortality/temperature data for this region
    region_results : dict
        DLNM results for this region from 01a
    region_code : str
        Region identifier
        
    Returns:
    --------
    dict with burden estimates at both threshold levels
    """
    if region_code not in region_results:
        return None
    
    model = region_results[region_code]
    
    # Get model parameters
    cb_coefs = np.array(model['cb_coefs'])
    col_names = model['col_names']
    temp_center = model['temp_center']
    temp_scale = model['temp_scale']
    temp_percentiles = model['temp_percentiles']
    
    # Reference temperature (P50 - median)
    ref_temp = temp_percentiles['p50']
    
    # Get temperatures and deaths
    temps = df_region['temp_mean'].values
    deaths = df_region['deaths_elderly'].values
    
    # Calculate temperature percentiles from daily data for this region
    # Use P2.5/P97.5 as primary thresholds (extreme but not rare)
    # Use P1/P99 as secondary thresholds (very extreme)
    p1 = np.percentile(temps, 1)
    p2_5 = np.percentile(temps, 2.5)
    p97_5 = np.percentile(temps, 97.5)
    p99 = np.percentile(temps, 99)
    
    # Calculate RR and AF for each day
    rr_values = np.array([
        compute_cumulative_rr(t, ref_temp, cb_coefs, col_names, temp_center, temp_scale) 
        for t in temps
    ])
    af_values = np.array([compute_attributable_fraction(rr) for rr in rr_values])
    
    # Totals
    total_deaths = deaths.sum()
    n_days = len(df_region)
    n_years = n_days / 365.25
    
    # Population (use elderly population)
    if 'pop_elderly' in df_region.columns:
        mean_pop = df_region['pop_elderly'].mean()
    else:
        mean_pop = None
    
    # === P2.5/P97.5 THRESHOLDS (PRIMARY) ===
    heat_mask_97_5 = temps > p97_5
    cold_mask_2_5 = temps < p2_5
    
    # Heat burden (only positive AF - excess risk)
    heat_af_97_5 = np.where((af_values > 0) & heat_mask_97_5, af_values, 0)
    heat_an_97_5 = np.sum(heat_af_97_5 * deaths)
    
    # Cold burden (only positive AF - excess risk)
    cold_af_2_5 = np.where((af_values > 0) & cold_mask_2_5, af_values, 0)
    cold_an_2_5 = np.sum(cold_af_2_5 * deaths)
    
    # === P1/P99 THRESHOLDS (COMPARISON) ===
    heat_mask_99 = temps > p99
    cold_mask_1 = temps < p1
    
    heat_af_99 = np.where((af_values > 0) & heat_mask_99, af_values, 0)
    heat_an_99 = np.sum(heat_af_99 * deaths)
    
    cold_af_1 = np.where((af_values > 0) & cold_mask_1, af_values, 0)
    cold_an_1 = np.sum(cold_af_1 * deaths)
    
    result = {
        'region_code': region_code,
        'n_days': n_days,
        'n_years': n_years,
        'total_deaths': float(total_deaths),
        'ref_temp': float(ref_temp),
        'mean_pop': float(mean_pop) if mean_pop else None,
        
        # P2.5/P97.5 thresholds (PRIMARY)
        'p2_5': float(p2_5),
        'p97_5': float(p97_5),
        'n_heat_days_97_5': int(heat_mask_97_5.sum()),
        'n_cold_days_2_5': int(cold_mask_2_5.sum()),
        'heat_an_97_5': float(heat_an_97_5),
        'cold_an_2_5': float(cold_an_2_5),
        'heat_af_97_5': float(heat_an_97_5 / total_deaths * 100) if total_deaths > 0 else 0,
        'cold_af_2_5': float(cold_an_2_5 / total_deaths * 100) if total_deaths > 0 else 0,
        'heat_annual_97_5': float(heat_an_97_5 / n_years),
        'cold_annual_2_5': float(cold_an_2_5 / n_years),
        'heat_rate_per_100k_97_5': float(heat_an_97_5 / n_years / mean_pop * 100000) if mean_pop else None,
        'cold_rate_per_100k_2_5': float(cold_an_2_5 / n_years / mean_pop * 100000) if mean_pop else None,
        
        # P1/P99 thresholds (COMPARISON)
        'p1': float(p1),
        'p99': float(p99),
        'n_heat_days_99': int(heat_mask_99.sum()),
        'n_cold_days_1': int(cold_mask_1.sum()),
        'heat_an_99': float(heat_an_99),
        'cold_an_1': float(cold_an_1),
        'heat_af_99': float(heat_an_99 / total_deaths * 100) if total_deaths > 0 else 0,
        'cold_af_1': float(cold_an_1 / total_deaths * 100) if total_deaths > 0 else 0,
        'heat_annual_99': float(heat_an_99 / n_years),
        'cold_annual_1': float(cold_an_1 / n_years),
        'heat_rate_per_100k_99': float(heat_an_99 / n_years / mean_pop * 100000) if mean_pop else None,
        'cold_rate_per_100k_1': float(cold_an_1 / n_years / mean_pop * 100000) if mean_pop else None,
    }
    
    return result


def calculate_burden_at_level(df, dlnm_results, level_name, region_col):
    """
    Calculate attributable burden at a given aggregation level.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Full mortality/temperature dataset
    dlnm_results : dict
        DLNM results from 01a
    level_name : str
        'intermediate' or 'immediate'
    region_col : str
        Column name for region code ('region_code' or 'immediate_code')
    """
    print(f"\n[{level_name.upper()}] Calculating attributable burden...")
    
    region_results = dlnm_results['region_results']
    region_codes = list(region_results.keys())
    
    # Make sure region column is string
    df[region_col] = df[region_col].astype(str)
    
    results = []
    successful = 0
    
    for region_code in region_codes:
        df_region = df[df[region_col] == region_code].copy()
        
        if len(df_region) < 365:
            continue
        
        burden = calculate_burden_for_region(df_region, region_results, region_code)
        
        if burden:
            results.append(burden)
            successful += 1
    
    print(f"  Calculated burden for {successful} / {len(region_codes)} regions")
    
    # Aggregate to national level
    results_df = pd.DataFrame(results)
    
    # Primary thresholds: P2.5 / P97.5
    national = {
        'level': level_name,
        'n_regions': len(results),
        'total_deaths': results_df['total_deaths'].sum(),
        
        # P2.5/P97.5 (PRIMARY)
        'heat_an_97_5': results_df['heat_an_97_5'].sum(),
        'cold_an_2_5': results_df['cold_an_2_5'].sum(),
        'heat_af_97_5': results_df['heat_an_97_5'].sum() / results_df['total_deaths'].sum() * 100,
        'cold_af_2_5': results_df['cold_an_2_5'].sum() / results_df['total_deaths'].sum() * 100,
        'heat_annual_97_5': results_df['heat_annual_97_5'].sum(),
        'cold_annual_2_5': results_df['cold_annual_2_5'].sum(),
        
        # P1/P99 (COMPARISON)
        'heat_an_99': results_df['heat_an_99'].sum(),
        'cold_an_1': results_df['cold_an_1'].sum(),
        'heat_af_99': results_df['heat_an_99'].sum() / results_df['total_deaths'].sum() * 100,
        'cold_af_1': results_df['cold_an_1'].sum() / results_df['total_deaths'].sum() * 100,
        'heat_annual_99': results_df['heat_annual_99'].sum(),
        'cold_annual_1': results_df['cold_annual_1'].sum(),
    }
    
    # Calculate rates if population available
    if results_df['mean_pop'].notna().all():
        total_pop = results_df['mean_pop'].sum()
        n_years = results_df['n_years'].mean()
        national['total_pop_elderly'] = total_pop
        # P2.5/P97.5 rates
        national['heat_rate_per_100k_97_5'] = national['heat_annual_97_5'] / total_pop * 100000
        national['cold_rate_per_100k_2_5'] = national['cold_annual_2_5'] / total_pop * 100000
        # P1/P99 rates
        national['heat_rate_per_100k_99'] = national['heat_annual_99'] / total_pop * 100000
        national['cold_rate_per_100k_1'] = national['cold_annual_1'] / total_pop * 100000
    
    return results_df, national


# =============================================================================
# CALCULATE BURDEN AT BOTH LEVELS
# =============================================================================

print("\n[3] Calculating attributable burden...")

# Intermediate level (133 regions)
intermediate_burden_df, intermediate_national = calculate_burden_at_level(
    df=df_int,
    dlnm_results=intermediate_results,
    level_name='intermediate',
    region_col='region_code'
)

# Immediate level (510 regions)
immediate_burden_df, immediate_national = calculate_burden_at_level(
    df=df_imm,
    dlnm_results=immediate_results,
    level_name='immediate',
    region_col='immediate_code'
)

# =============================================================================
# DISPLAY RESULTS
# =============================================================================

print("\n" + "=" * 70)
print("ATTRIBUTABLE BURDEN RESULTS")
print("=" * 70)

print("\n[INTERMEDIATE LEVEL - 133 Regions]")
print(f"  Total deaths: {intermediate_national['total_deaths']:,.0f}")
print("  --- P2.5/P97.5 Thresholds (PRIMARY) ---")
print(f"  Heat-attributable deaths (>P97.5): {intermediate_national['heat_an_97_5']:,.0f} ({intermediate_national['heat_af_97_5']:.2f}%)")
print(f"  Cold-attributable deaths (<P2.5): {intermediate_national['cold_an_2_5']:,.0f} ({intermediate_national['cold_af_2_5']:.2f}%)")
print(f"  Annual heat deaths: {intermediate_national['heat_annual_97_5']:,.0f}")
print(f"  Annual cold deaths: {intermediate_national['cold_annual_2_5']:,.0f}")
print("  --- P1/P99 Thresholds (COMPARISON) ---")
print(f"  Heat-attributable deaths (>P99): {intermediate_national['heat_an_99']:,.0f} ({intermediate_national['heat_af_99']:.2f}%)")
print(f"  Cold-attributable deaths (<P1): {intermediate_national['cold_an_1']:,.0f} ({intermediate_national['cold_af_1']:.2f}%)")
print(f"  Annual heat deaths: {intermediate_national['heat_annual_99']:,.0f}")
print(f"  Annual cold deaths: {intermediate_national['cold_annual_1']:,.0f}")

print("\n[IMMEDIATE LEVEL - 510 Regions]")
print(f"  Total deaths: {immediate_national['total_deaths']:,.0f}")
print("  --- P2.5/P97.5 Thresholds (PRIMARY) ---")
print(f"  Heat-attributable deaths (>P97.5): {immediate_national['heat_an_97_5']:,.0f} ({immediate_national['heat_af_97_5']:.2f}%)")
print(f"  Cold-attributable deaths (<P2.5): {immediate_national['cold_an_2_5']:,.0f} ({immediate_national['cold_af_2_5']:.2f}%)")
print(f"  Annual heat deaths: {immediate_national['heat_annual_97_5']:,.0f}")
print(f"  Annual cold deaths: {immediate_national['cold_annual_2_5']:,.0f}")
print("  --- P1/P99 Thresholds (COMPARISON) ---")
print(f"  Heat-attributable deaths (>P99): {immediate_national['heat_an_99']:,.0f} ({immediate_national['heat_af_99']:.2f}%)")
print(f"  Cold-attributable deaths (<P1): {immediate_national['cold_an_1']:,.0f} ({immediate_national['cold_af_1']:.2f}%)")
print(f"  Annual heat deaths: {immediate_national['heat_annual_99']:,.0f}")
print(f"  Annual cold deaths: {immediate_national['cold_annual_1']:,.0f}")

# =============================================================================
# REGIONAL SUMMARY
# =============================================================================

print("\n[4] Regional burden distribution...")

for level_name, burden_df in [('Intermediate', intermediate_burden_df), ('Immediate', immediate_burden_df)]:
    print(f"\n  {level_name} Level (P2.5/P97.5 thresholds):")
    print(f"    Heat AF range: {burden_df['heat_af_97_5'].min():.2f}% - {burden_df['heat_af_97_5'].max():.2f}%")
    print(f"    Cold AF range: {burden_df['cold_af_2_5'].min():.2f}% - {burden_df['cold_af_2_5'].max():.2f}%")
    print(f"    Mean Heat AF: {burden_df['heat_af_97_5'].mean():.2f}%")
    print(f"    Mean Cold AF: {burden_df['cold_af_2_5'].mean():.2f}%")

# =============================================================================
# SAVE RESULTS
# =============================================================================

print("\n[5] Saving results...")

# Convert numpy types for JSON
def convert_for_json(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: convert_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_for_json(i) for i in obj]
    return obj

# Save regional burden CSVs
intermediate_burden_df.to_csv(OUTPUT_DIR / 'burden_intermediate_regions.csv', index=False)
immediate_burden_df.to_csv(OUTPUT_DIR / 'burden_immediate_regions.csv', index=False)

# Save national summary
national_summary = {
    'intermediate': convert_for_json(intermediate_national),
    'immediate': convert_for_json(immediate_national),
    'timestamp': datetime.now().isoformat(),
    'method': 'Attributable fraction (RR-1)/RR × deaths, thresholds P2.5/P97.5 (primary) and P1/P99 (comparison), reference P50',
    'source': 'DLNM results from 01a_intermediate_dlnm.py and 01a_immediate_dlnm.py'
}

with open(OUTPUT_DIR / 'attributable_burden_national.json', 'w') as f:
    json.dump(national_summary, f, indent=2)

print(f"  Saved: burden_intermediate_regions.csv")
print(f"  Saved: burden_immediate_regions.csv")
print(f"  Saved: attributable_burden_national.json")

# =============================================================================
# SUMMARY TABLE
# =============================================================================

print("\n" + "=" * 70)
print("SUMMARY TABLE - PRIMARY (P2.5/P97.5 thresholds)")
print("=" * 70)

summary_data = [
    {
        'Level': 'Intermediate (133)',
        'Total Deaths': f"{intermediate_national['total_deaths']:,.0f}",
        'Heat AN (>P97.5)': f"{intermediate_national['heat_an_97_5']:,.0f}",
        'Cold AN (<P2.5)': f"{intermediate_national['cold_an_2_5']:,.0f}",
        'Heat AF (%)': f"{intermediate_national['heat_af_97_5']:.2f}",
        'Cold AF (%)': f"{intermediate_national['cold_af_2_5']:.2f}",
    },
    {
        'Level': 'Immediate (510)',
        'Total Deaths': f"{immediate_national['total_deaths']:,.0f}",
        'Heat AN (>P97.5)': f"{immediate_national['heat_an_97_5']:,.0f}",
        'Cold AN (<P2.5)': f"{immediate_national['cold_an_2_5']:,.0f}",
        'Heat AF (%)': f"{immediate_national['heat_af_97_5']:.2f}",
        'Cold AF (%)': f"{immediate_national['cold_af_2_5']:.2f}",
    }
]

summary_df = pd.DataFrame(summary_data)
print(summary_df.to_string(index=False))

print("\n" + "-" * 70)
print("COMPARISON TABLE - (P1/P99 thresholds)")
print("-" * 70)

comparison_data = [
    {
        'Level': 'Intermediate (133)',
        'Heat AN (>P99)': f"{intermediate_national['heat_an_99']:,.0f}",
        'Cold AN (<P1)': f"{intermediate_national['cold_an_1']:,.0f}",
        'Heat AF (%)': f"{intermediate_national['heat_af_99']:.2f}",
        'Cold AF (%)': f"{intermediate_national['cold_af_1']:.2f}",
    },
    {
        'Level': 'Immediate (510)',
        'Heat AN (>P99)': f"{immediate_national['heat_an_99']:,.0f}",
        'Cold AN (<P1)': f"{immediate_national['cold_an_1']:,.0f}",
        'Heat AF (%)': f"{immediate_national['heat_af_99']:.2f}",
        'Cold AF (%)': f"{immediate_national['cold_af_1']:.2f}",
    }
]

comparison_df = pd.DataFrame(comparison_data)
print(comparison_df.to_string(index=False))

# Save summary
all_summary = pd.concat([summary_df, comparison_df], axis=1)
all_summary.to_csv(OUTPUT_DIR / 'attributable_burden_summary.csv', index=False)

print("\n" + "=" * 70)
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)
