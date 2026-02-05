"""
01d_attributable_burden_v2.py
==============================
Attributable Burden Calculation for DLNM v2 (Natural Spline Cross-Basis)

CRITICAL: This script must match the basis construction in 01a_*_dlnm_v2.py.
If the DLNM basis type changes, this reconstruction will produce GARBAGE.

Key Changes from v1:
1. Auto-detects basis type (polynomial vs natural spline) from DLNM results
2. Uses region-specific MMT as reference (not P50)
3. Clearly documents AF sign handling
4. Adds edge case guards and validation
5. Validates threshold separation

Method (Gasparrini & Leone, 2014):
1. Load fitted DLNM models from 01a_v2 
2. For each region, reconstruct E-R curve using stored coefficients
3. Calculate AF = (RR - 1) / RR for each day relative to MMT
4. Sum attributable deaths for heat (temp > MMT) and cold (temp < MMT)
5. Aggregate to national level with uncertainty

Author: Climate-Health Analysis Pipeline
Date: December 2025
"""

import warnings
warnings.filterwarnings('ignore')

import os
import sys
import json
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path

# Import shared DLNM utilities
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.dlnm_module import (
    # Natural spline basis
    ns_basis,
    # Burden calculation functions
    compute_cumulative_rr_for_burden,
    compute_attributable_fraction,
    detect_basis_type,
    convert_to_json_serializable,
)

# =============================================================================
# PATHS
# =============================================================================

SCRIPT_DIR = Path(__file__).parent
PHASE0_RESULTS = SCRIPT_DIR.parent / 'phase0_data_prep' / 'results'
PHASE1_RESULTS = SCRIPT_DIR / 'results'
OUTPUT_DIR = PHASE1_RESULTS
OUTPUT_DIR.mkdir(exist_ok=True)

print("=" * 70)
print("01d: ATTRIBUTABLE BURDEN CALCULATION v2")
print("       (Compatible with Natural Spline DLNM)")
print("=" * 70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# NATURAL SPLINE BASIS FUNCTIONS (imported from utils.dlnm_module)
# =============================================================================
# The following functions are imported from utils/dlnm_module.py:
#   - ns_basis
#   - compute_cumulative_rr_for_burden
#   - compute_attributable_fraction  
#   - detect_basis_type
#   - convert_to_json_serializable
# This ensures consistency across all analysis scripts.


def compute_cumulative_rr_poly(target_temp, ref_temp, cb_coefs, col_names, temp_center, temp_scale):
    """
    Compute cumulative RR for polynomial cross-basis (v1 compatibility).
    This is kept locally as it's only needed for backward compatibility with v1 results.
    """
    temp_std_target = (target_temp - temp_center) / temp_scale
    temp_std_ref = (ref_temp - temp_center) / temp_scale
    
    contrast = np.zeros(len(col_names))
    
    for i, name in enumerate(col_names):
        parts = name.split('_')
        poly = int(parts[1].replace('poly', ''))
        contrast[i] = temp_std_target**poly - temp_std_ref**poly
    
    log_rr = np.dot(contrast, cb_coefs)
    return float(np.exp(log_rr))


def calculate_burden_for_region_v2(df_region, model_result, region_code, basis_type):
    """
    Calculate attributable burden for a single region using MMT as reference.
    
    Parameters:
    -----------
    df_region : pd.DataFrame
        Daily mortality/temperature data
    model_result : dict
        DLNM model results for this region
    region_code : str
        Region identifier
    basis_type : str
        'ns' for natural spline, 'poly' for polynomial
        
    Returns:
    --------
    dict with burden estimates
    """
    # Get temperature and deaths
    temps = df_region['temp_mean'].values
    deaths = df_region['deaths_elderly'].values
    
    # === EDGE CASE GUARDS ===
    n_days = len(df_region)
    if n_days == 0:
        return None
    
    n_years = n_days / 365.25
    if n_years <= 0:
        return None
    
    total_deaths = deaths.sum()
    if total_deaths <= 0:
        return None
    
    # Population (with guard)
    if 'pop_elderly' in df_region.columns:
        mean_pop = df_region['pop_elderly'].mean()
        if mean_pop <= 0 or np.isnan(mean_pop):
            mean_pop = None
    else:
        mean_pop = None
    
    # === GET REFERENCE TEMPERATURE (MMT, not P50!) ===
    if 'mmt' in model_result:
        ref_temp = model_result['mmt']
        mmt_percentile = model_result.get('mmt_percentile', 50)
    else:
        # Fallback to P50 if MMT not available (v1 compatibility)
        ref_temp = model_result['temp_percentiles']['p50']
        mmt_percentile = 50
        print(f"  WARNING: Region {region_code} using P50 fallback (no MMT)")
    
    # === COMPUTE RR FOR EACH DAY ===
    if basis_type == 'ns':
        cb_info = model_result['crossbasis_info']
        cb_coefs = np.array(model_result['cb_coefs'])
        
        rr_values = np.array([
            compute_cumulative_rr_for_burden(t, ref_temp, cb_coefs, cb_info)
            for t in temps
        ])
    else:  # poly
        cb_coefs = np.array(model_result['cb_coefs'])
        col_names = model_result['col_names']
        temp_center = model_result['temp_center']
        temp_scale = model_result['temp_scale']
        
        rr_values = np.array([
            compute_cumulative_rr_poly(t, ref_temp, cb_coefs, col_names, temp_center, temp_scale)
            for t in temps
        ])
    
    # === COMPUTE ATTRIBUTABLE FRACTIONS ===
    af_values = np.array([compute_attributable_fraction(rr) for rr in rr_values])
    
    # === TEMPERATURE THRESHOLDS ===
    # Region-specific percentiles
    p1 = np.percentile(temps, 1)
    p2_5 = np.percentile(temps, 2.5)
    p97_5 = np.percentile(temps, 97.5)
    p99 = np.percentile(temps, 99)
    
    # Check threshold separation
    temp_range = p99 - p1
    if temp_range < 5:  # Very narrow range
        print(f"  WARNING: Region {region_code} has narrow temp range ({temp_range:.1f}°C)")
    
    # === HEAT/COLD CLASSIFICATION BASED ON MMT ===
    # Heat: days where temperature > MMT (not just > P97.5)
    # Cold: days where temperature < MMT (not just < P2.5)
    heat_mask_mmt = temps > ref_temp
    cold_mask_mmt = temps < ref_temp
    
    # Threshold masks for reporting
    heat_mask_97_5 = temps > p97_5
    cold_mask_2_5 = temps < p2_5
    heat_mask_99 = temps > p99
    cold_mask_1 = temps < p1
    
    # === BURDEN CALCULATION ===
    # Option 1: Only count positive AF (excess mortality)
    # This is the standard approach - we attribute deaths to temperature
    # only when RR > 1 (temperature increases risk)
    
    # Heat burden: positive AF on days above threshold
    heat_af_97_5 = np.where((af_values > 0) & heat_mask_97_5, af_values, 0)
    heat_an_97_5 = np.sum(heat_af_97_5 * deaths)
    
    cold_af_2_5 = np.where((af_values > 0) & cold_mask_2_5, af_values, 0)
    cold_an_2_5 = np.sum(cold_af_2_5 * deaths)
    
    heat_af_99 = np.where((af_values > 0) & heat_mask_99, af_values, 0)
    heat_an_99 = np.sum(heat_af_99 * deaths)
    
    cold_af_1 = np.where((af_values > 0) & cold_mask_1, af_values, 0)
    cold_an_1 = np.sum(cold_af_1 * deaths)
    
    # Total attributable to non-optimal temperature (any temp != MMT)
    # Only positive AF
    total_heat_af = np.where((af_values > 0) & heat_mask_mmt, af_values, 0)
    total_heat_an = np.sum(total_heat_af * deaths)
    
    total_cold_af = np.where((af_values > 0) & cold_mask_mmt, af_values, 0)
    total_cold_an = np.sum(total_cold_af * deaths)
    
    # Option 2: Also track protective effects (for sensitivity analysis)
    # This shows where RR < 1 (temperature is protective relative to MMT)
    protective_heat = np.where((af_values < 0) & heat_mask_mmt, -af_values * deaths, 0).sum()
    protective_cold = np.where((af_values < 0) & cold_mask_mmt, -af_values * deaths, 0).sum()
    
    # === BUILD RESULT ===
    result = {
        'region_code': region_code,
        'n_days': int(n_days),
        'n_years': float(n_years),
        'total_deaths': float(total_deaths),
        'mmt': float(ref_temp),
        'mmt_percentile': float(mmt_percentile),
        'mean_pop': float(mean_pop) if mean_pop else None,
        'temp_range': float(temp_range),
        
        # Total attributable to non-optimal temperature
        'total_heat_an': float(total_heat_an),
        'total_cold_an': float(total_cold_an),
        'total_heat_af_pct': float(total_heat_an / total_deaths * 100),
        'total_cold_af_pct': float(total_cold_an / total_deaths * 100),
        
        # P2.5/P97.5 thresholds (PRIMARY)
        'p2_5': float(p2_5),
        'p97_5': float(p97_5),
        'n_heat_days_97_5': int(heat_mask_97_5.sum()),
        'n_cold_days_2_5': int(cold_mask_2_5.sum()),
        'heat_an_97_5': float(heat_an_97_5),
        'cold_an_2_5': float(cold_an_2_5),
        'heat_af_pct_97_5': float(heat_an_97_5 / total_deaths * 100),
        'cold_af_pct_2_5': float(cold_an_2_5 / total_deaths * 100),
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
        'heat_af_pct_99': float(heat_an_99 / total_deaths * 100),
        'cold_af_pct_1': float(cold_an_1 / total_deaths * 100),
        'heat_annual_99': float(heat_an_99 / n_years),
        'cold_annual_1': float(cold_an_1 / n_years),
        'heat_rate_per_100k_99': float(heat_an_99 / n_years / mean_pop * 100000) if mean_pop else None,
        'cold_rate_per_100k_1': float(cold_an_1 / n_years / mean_pop * 100000) if mean_pop else None,
        
        # Protective effects (for sensitivity analysis)
        'protective_heat_deaths': float(protective_heat),
        'protective_cold_deaths': float(protective_cold),
    }
    
    return result


# =============================================================================
# LOAD DATA
# =============================================================================

print("\n[1] Loading data...")

# Load intermediate data
df_temp_int = pd.read_parquet(PHASE0_RESULTS / 'era5_intermediate_daily.parquet')
df_temp_int['date'] = pd.to_datetime(df_temp_int['date'])

df_mort_int = pd.read_parquet(PHASE0_RESULTS / 'mortality_regional_daily_elderly.parquet')
df_mort_int['date'] = pd.to_datetime(df_mort_int['date'])

df_ses_int = pd.read_csv(PHASE0_RESULTS / 'ses_intermediate_covariates.csv')
pop_map_int = dict(zip(df_ses_int['intermediate_code'].astype(str), df_ses_int['pop_elderly']))

df_int = pd.merge(
    df_mort_int,
    df_temp_int[['date', 'region_code', 'temp_mean']],
    on=['date', 'region_code'],
    how='inner'
)
df_int['region_code'] = df_int['region_code'].astype(str)
df_int['pop_elderly'] = df_int['region_code'].map(pop_map_int)
print(f"  Intermediate: {len(df_int):,} rows, {df_int['region_code'].nunique()} regions")

# Load immediate data
df_temp_imm = pd.read_parquet(PHASE0_RESULTS / 'era5_immediate_daily.parquet')
df_temp_imm['date'] = pd.to_datetime(df_temp_imm['date'])
df_temp_imm = df_temp_imm.rename(columns={'region_code': 'immediate_code'})

df_mort_imm = pd.read_parquet(PHASE0_RESULTS / 'mortality_immediate_daily_elderly.parquet')
df_mort_imm['date'] = pd.to_datetime(df_mort_imm['date'])

df_ses_imm = pd.read_csv(PHASE0_RESULTS / 'ses_immediate_covariates.csv')
pop_map_imm = dict(zip(df_ses_imm['immediate_code'].astype(str), df_ses_imm['pop_elderly']))

df_imm = pd.merge(
    df_mort_imm,
    df_temp_imm[['date', 'immediate_code', 'temp_mean']],
    on=['date', 'immediate_code'],
    how='inner'
)
df_imm['immediate_code'] = df_imm['immediate_code'].astype(str)
df_imm['pop_elderly'] = df_imm['immediate_code'].map(pop_map_imm)
print(f"  Immediate: {len(df_imm):,} rows, {df_imm['immediate_code'].nunique()} regions")

# =============================================================================
# LOAD DLNM RESULTS - TRY V2 FIRST, FALLBACK TO V1
# =============================================================================

print("\n[2] Loading DLNM results...")

# Try v2 first (natural splines with MMT)
try:
    with open(PHASE1_RESULTS / 'dlnm_v2_intermediate_results.json', 'r') as f:
        intermediate_results = json.load(f)
    print(f"  Intermediate: Loaded v2 results (natural spline)")
except FileNotFoundError:
    with open(PHASE1_RESULTS / 'dlnm_intermediate_results.json', 'r') as f:
        intermediate_results = json.load(f)
    print(f"  Intermediate: Loaded v1 results (polynomial)")

try:
    with open(PHASE1_RESULTS / 'dlnm_v2_immediate_results.json', 'r') as f:
        immediate_results = json.load(f)
    print(f"  Immediate: Loaded v2 results (natural spline)")
except FileNotFoundError:
    with open(PHASE1_RESULTS / 'dlnm_immediate_results.json', 'r') as f:
        immediate_results = json.load(f)
    print(f"  Immediate: Loaded v1 results (polynomial)")

# Detect basis type
int_basis = detect_basis_type(intermediate_results)
imm_basis = detect_basis_type(immediate_results)
print(f"  Intermediate basis type: {int_basis}")
print(f"  Immediate basis type: {imm_basis}")

print(f"  Intermediate regions: {intermediate_results['n_successful']} / {intermediate_results['n_regions']}")
print(f"  Immediate regions: {immediate_results['n_successful']} / {immediate_results['n_regions']}")


# =============================================================================
# CALCULATE BURDEN
# =============================================================================

def calculate_burden_at_level(df, dlnm_results, level_name, region_col, basis_type):
    """Calculate burden for all regions at a given level."""
    print(f"\n[{level_name.upper()}] Calculating attributable burden...")
    
    region_results = dlnm_results['region_results']
    region_codes = list(region_results.keys())
    
    df[region_col] = df[region_col].astype(str)
    
    results = []
    successful = 0
    
    for region_code in region_codes:
        df_region = df[df[region_col] == region_code].copy()
        
        if len(df_region) < 365:
            continue
        
        burden = calculate_burden_for_region_v2(
            df_region, 
            region_results[region_code], 
            region_code,
            basis_type
        )
        
        if burden:
            results.append(burden)
            successful += 1
    
    print(f"  Calculated burden for {successful} / {len(region_codes)} regions")
    
    return results


# Run calculations
int_burdens = calculate_burden_at_level(
    df_int, intermediate_results, 'intermediate', 'region_code', int_basis
)
imm_burdens = calculate_burden_at_level(
    df_imm, immediate_results, 'immediate', 'immediate_code', imm_basis
)

# =============================================================================
# AGGREGATE TO NATIONAL LEVEL
# =============================================================================

def aggregate_national(burdens, level_name):
    """Aggregate region burdens to national level."""
    if not burdens:
        return {}
    
    df = pd.DataFrame(burdens)
    
    total_deaths = df['total_deaths'].sum()
    total_pop = df['mean_pop'].dropna().sum()
    total_years = df['n_years'].mean()  # Average years of data
    
    national = {
        'level': level_name,
        'n_regions': len(df),
        'total_deaths': float(total_deaths),
        'total_pop': float(total_pop) if total_pop > 0 else None,
        'avg_years': float(total_years),
        
        # Mean MMT across regions
        'mean_mmt': float(df['mmt'].mean()),
        'mean_mmt_percentile': float(df['mmt_percentile'].mean()),
        
        # Total attributable (all days above/below MMT)
        'total_heat_an': float(df['total_heat_an'].sum()),
        'total_cold_an': float(df['total_cold_an'].sum()),
        'total_heat_af_pct': float(df['total_heat_an'].sum() / total_deaths * 100),
        'total_cold_af_pct': float(df['total_cold_an'].sum() / total_deaths * 100),
        
        # P2.5/P97.5 (PRIMARY)
        'heat_an_97_5': float(df['heat_an_97_5'].sum()),
        'cold_an_2_5': float(df['cold_an_2_5'].sum()),
        'heat_af_pct_97_5': float(df['heat_an_97_5'].sum() / total_deaths * 100),
        'cold_af_pct_2_5': float(df['cold_an_2_5'].sum() / total_deaths * 100),
        'heat_annual_97_5': float(df['heat_annual_97_5'].sum()),
        'cold_annual_2_5': float(df['cold_annual_2_5'].sum()),
        
        # P1/P99 (COMPARISON)
        'heat_an_99': float(df['heat_an_99'].sum()),
        'cold_an_1': float(df['cold_an_1'].sum()),
        'heat_af_pct_99': float(df['heat_an_99'].sum() / total_deaths * 100),
        'cold_af_pct_1': float(df['cold_an_1'].sum() / total_deaths * 100),
        'heat_annual_99': float(df['heat_annual_99'].sum()),
        'cold_annual_1': float(df['cold_annual_1'].sum()),
        
        # Protective effects
        'protective_heat': float(df['protective_heat_deaths'].sum()),
        'protective_cold': float(df['protective_cold_deaths'].sum()),
    }
    
    # Add rates per 100k if population available
    if total_pop and total_pop > 0:
        national['heat_rate_97_5'] = float(df['heat_annual_97_5'].sum() / total_pop * 100000)
        national['cold_rate_2_5'] = float(df['cold_annual_2_5'].sum() / total_pop * 100000)
    
    return national


national_int = aggregate_national(int_burdens, 'intermediate')
national_imm = aggregate_national(imm_burdens, 'immediate')

# =============================================================================
# PRINT RESULTS
# =============================================================================

print("\n" + "=" * 70)
print("NATIONAL ATTRIBUTABLE BURDEN SUMMARY")
print("=" * 70)

for name, nat in [('INTERMEDIATE (133 regions)', national_int), 
                   ('IMMEDIATE (510 regions) - PRIMARY', national_imm)]:
    if not nat:
        continue
    
    print(f"\n{name}")
    print("-" * 50)
    print(f"  Regions: {nat['n_regions']}")
    print(f"  Total deaths: {nat['total_deaths']:,.0f}")
    print(f"  Mean MMT: P{nat['mean_mmt_percentile']:.1f} ({nat['mean_mmt']:.1f}°C)")
    
    print(f"\n  TOTAL (all days ≠ MMT):")
    print(f"    Heat attributable: {nat['total_heat_an']:,.0f} ({nat['total_heat_af_pct']:.2f}%)")
    print(f"    Cold attributable: {nat['total_cold_an']:,.0f} ({nat['total_cold_af_pct']:.2f}%)")
    
    print(f"\n  EXTREME HEAT (>P97.5) - PRIMARY:")
    print(f"    Attributable deaths: {nat['heat_an_97_5']:,.0f} ({nat['heat_af_pct_97_5']:.2f}%)")
    print(f"    Annual: {nat['heat_annual_97_5']:,.0f}")
    
    print(f"\n  EXTREME COLD (<P2.5) - PRIMARY:")
    print(f"    Attributable deaths: {nat['cold_an_2_5']:,.0f} ({nat['cold_af_pct_2_5']:.2f}%)")
    print(f"    Annual: {nat['cold_annual_2_5']:,.0f}")
    
    print(f"\n  VERY EXTREME (P1/P99) - COMPARISON:")
    print(f"    Heat (>P99): {nat['heat_an_99']:,.0f} ({nat['heat_af_pct_99']:.2f}%)")
    print(f"    Cold (<P1): {nat['cold_an_1']:,.0f} ({nat['cold_af_pct_1']:.2f}%)")
    
    if nat['protective_heat'] > 0 or nat['protective_cold'] > 0:
        print(f"\n  PROTECTIVE EFFECTS (RR < 1):")
        print(f"    Heat days: {nat['protective_heat']:,.0f} deaths avoided")
        print(f"    Cold days: {nat['protective_cold']:,.0f} deaths avoided")

# =============================================================================
# SAVE RESULTS
# =============================================================================

print("\n" + "=" * 70)
print("SAVING RESULTS")
print("=" * 70)

# Region-level results
for burdens, name in [(int_burdens, 'intermediate'), (imm_burdens, 'immediate')]:
    if burdens:
        df_out = pd.DataFrame(burdens)
        out_file = OUTPUT_DIR / f'burden_v2_{name}_regions.csv'
        df_out.to_csv(out_file, index=False)
        print(f"  Saved: {out_file}")

# National summary
summary = {
    'intermediate': national_int,
    'immediate': national_imm,
    'timestamp': datetime.now().isoformat(),
    'basis_types': {
        'intermediate': int_basis,
        'immediate': imm_basis
    }
}

out_file = OUTPUT_DIR / 'burden_v2_national_summary.json'
with open(out_file, 'w') as f:
    json.dump(summary, f, indent=2)
print(f"  Saved: {out_file}")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("\nNotes:")
print("  - Using region-specific MMT as reference (not P50)")
print("  - Only positive AF counted (excess mortality)")
print("  - P2.5/P97.5 are PRIMARY thresholds")
print("  - Protective effects tracked separately for sensitivity")
