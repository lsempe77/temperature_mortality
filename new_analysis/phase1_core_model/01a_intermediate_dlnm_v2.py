"""
01a_intermediate_dlnm_v2.py
============================
Two-Stage DLNM Analysis with TRUE Natural Cubic Spline Cross-Basis

Implements the methodology from Gasparrini et al. (2010, 2015):
- Natural cubic splines for temperature dimension
- Natural cubic splines for lag dimension  
- Tensor product cross-basis (12-16 parameters vs 66 in polynomial)

Key differences from v1 (polynomial):
1. Uses ns() basis for both dimensions (smoother, fewer params)
2. Proper knot placement at percentiles
3. Centering at MMT for RR computation
4. Compatible with MCC methodology

Input Data (from phase0_data_prep/results/):
- era5_intermediate_daily.parquet: Temperature data
- mortality_regional_daily_elderly.parquet: Elderly mortality  
- ses_intermediate_covariates.csv: Population for offset
- brazilian_holidays_daily.parquet: Holiday controls

Output:
- results/dlnm_v2_intermediate_results.json
- results/dlnm_v2_intermediate_summary.csv
- results/dlnm_v2_intermediate_pooled.json

Author: Climate-Health Analysis Pipeline
Date: December 2025
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod.families import Poisson
from patsy import dmatrix
from scipy import stats
from scipy.interpolate import BSpline
import os
import sys
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Import shared DLNM utilities
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.dlnm_module import (
    # Natural spline basis
    ns_basis, 
    create_lag_matrix, 
    create_crossbasis_ns,
    # Phase 1 compatible functions
    compute_cumulative_rr_ns_with_se,
    find_mmt_from_coefficients,
    compute_effects_relative_to_mmt,
    create_time_spline,
    random_effects_meta_analysis,
    pool_region_results,
    convert_to_json_serializable,
    # MVMeta for proper coefficient pooling (Gasparrini methodology)
    mvmeta_pool_regions,
)

print("="*70)
print("01a: INTERMEDIATE REGION DLNM v2 (Natural Spline Cross-Basis)")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PHASE0_RESULTS = os.path.join(BASE_DIR, 'new_analysis', 'phase0_data_prep', 'results')
OUTPUT_DIR = os.path.join(BASE_DIR, 'new_analysis', 'phase1_core_model', 'results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# DLNM parameters (following Gasparrini methodology)
MAX_LAG = 21              # Maximum lag in days
TEMP_DF = 4               # DF for temperature (natural spline)
LAG_DF = 4                # DF for lag structure (natural spline with log scale)
TIME_SPLINE_DF_PER_YEAR = 1  # DF per year for long-term trend (coarse trend; month dummies handle seasonality)
MIN_OBS = 1000            # Minimum observations per region

# Temperature knot percentiles (Gasparrini standard)
TEMP_KNOT_PCTS = [10, 75, 90]  # Interior knots at these percentiles
TEMP_BOUNDARY_PCTS = [1, 99]   # Boundary knots

# Lag knots (log scale for declining effect)
LAG_KNOTS = [1, 3, 7, 14]  # Days - captures fast initial + slow decay

print(f"\nConfiguration (Gasparrini methodology):")
print(f"  Max lag: {MAX_LAG} days")
print(f"  Temperature DF: {TEMP_DF} (ns with knots at P{TEMP_KNOT_PCTS})")
print(f"  Lag DF: {LAG_DF} (ns with knots at days {LAG_KNOTS})")
print(f"  Time trend: {TIME_SPLINE_DF_PER_YEAR} df/year (natural spline)")
print(f"  Minimum obs per region: {MIN_OBS}")

# =============================================================================
# NATURAL CUBIC SPLINE FUNCTIONS (imported from utils.dlnm_module)
# =============================================================================
# All core DLNM functions are imported from utils/dlnm_module.py:
#   - ns_basis, create_lag_matrix, create_crossbasis_ns
#   - compute_cumulative_rr_ns_with_se, find_mmt_from_coefficients
#   - compute_effects_relative_to_mmt, create_time_spline
#   - random_effects_meta_analysis, pool_region_results
#   - convert_to_json_serializable
# This ensures consistency across all analysis scripts.


# =============================================================================
# MODEL FITTING
# =============================================================================

def fit_region_dlnm(df_region, region_code, max_lag=MAX_LAG):
    """
    Fit DLNM with natural spline cross-basis for a single region.
    """
    # Basic size guard
    if len(df_region) < MIN_OBS:
        print(f"    [DEBUG] Region {region_code}: n_rows={len(df_region)} < MIN_OBS={MIN_OBS} (skipping)")
        return None
    
    df_region = df_region.sort_values('date').reset_index(drop=True)
    
    # Compute temperature percentiles for this region
    temp = df_region['temp_mean'].values
    temp_valid = temp[~np.isnan(temp)]
    
    temp_knots = np.percentile(temp_valid, TEMP_KNOT_PCTS)
    temp_boundary = (np.percentile(temp_valid, TEMP_BOUNDARY_PCTS[0]),
                     np.percentile(temp_valid, TEMP_BOUNDARY_PCTS[1]))
    
    # Create cross-basis
    X_cb, cb_info = create_crossbasis_ns(temp, max_lag, temp_knots, temp_boundary, LAG_KNOTS)
    
    # Find valid rows (drop initial rows with insufficient lag history or missing data)
    valid = ~np.isnan(X_cb).any(axis=1)
    n_valid = int(valid.sum())
    if n_valid < MIN_OBS:
        print(f"    [DEBUG] Region {region_code}: valid_rows={n_valid} < MIN_OBS={MIN_OBS} after lag construction (skipping)")
        return None
    
    X_cb_valid = X_cb[valid]
    df_valid = df_region.loc[valid].copy()
    y = df_valid['deaths_elderly'].values
    
    # Population offset
    pop = df_valid['pop_elderly'].values
    offset = np.log(pop / 100000)
    
    # Control variables
    month_dummies = pd.get_dummies(df_valid['month'], prefix='month', drop_first=True)
    dow_dummies = pd.get_dummies(df_valid['day_of_week'], prefix='dow', drop_first=True)
    n_years = (df_valid['date'].max() - df_valid['date'].min()).days / 365.25
    time_spline = create_time_spline(df_valid['time_index'].values, n_years, TIME_SPLINE_DF_PER_YEAR)
    holiday = df_valid['is_holiday'].values.reshape(-1, 1)
    
    X_controls = np.column_stack([
        month_dummies.values,
        dow_dummies.values,
        time_spline,
        holiday
    ])
    
    X_full = np.column_stack([X_cb_valid, X_controls])
    X_full = sm.add_constant(X_full)
    
    # Fit Quasi-Poisson GLM (allow more iterations for challenging regions)
    try:
        model = GLM(y, X_full, family=Poisson(), offset=offset).fit(scale='X2', maxiter=200)
        
        if not model.converged:
            print(f"    [DEBUG] Region {region_code}: GLM did not converge (n_valid={len(y)}, total_deaths={int(y.sum())}, dispersion={getattr(model, 'scale', float('nan')):.3f})")
            return None
        
        dispersion = model.scale
        
    except Exception as e:
        print(f"    [DEBUG] Region {region_code}: ERROR during GLM fit: {type(e).__name__}: {e}")
        return None
    
    # Extract cross-basis coefficients
    n_cb = cb_info['n_params']
    cb_coefs = model.params[1:n_cb+1]
    
    vcov_full = model.cov_params()
    if hasattr(vcov_full, 'iloc'):
        cb_vcov = vcov_full.iloc[1:n_cb+1, 1:n_cb+1].values
    else:
        cb_vcov = vcov_full[1:n_cb+1, 1:n_cb+1]
    
    # Compute temperature percentiles
    temps_valid = df_valid['temp_mean'].dropna()
    percentiles = [1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99]
    temp_pcts = {f'p{p}'.replace('.', '_'): float(temps_valid.quantile(p/100)) for p in percentiles}
    temp_pcts = {k.replace('_', '.'): v for k, v in temp_pcts.items()}
    
    # Step 1: Create fine temperature grid for MMT search (P1 to P99)
    # Following Gasparrini: search the FULL range, no constraint
    temp_min = np.percentile(temps_valid, 1)
    temp_max = np.percentile(temps_valid, 99)
    temp_grid = np.linspace(temp_min, temp_max, 100)
    
    # Step 2: Find MMT using P50 as initial computational centering
    # No constraint - let MMT fall wherever it naturally is (Gasparrini methodology)
    centering_temp = temp_pcts['p50']
    mmt, mmt_percentile_grid, rr_curve = find_mmt_from_coefficients(
        temp_grid, cb_coefs, cb_vcov, cb_info, centering_temp,
        constrain_to_interior=False  # Gasparrini: no constraint on MMT
    )
    
    # Convert MMT temperature into an empirical percentile of the region's
    # actual temperature distribution
    mmt_percentile_emp = float((temps_valid <= mmt).mean() * 100.0)
    
    # Step 3: Compute all effects relative to MMT (standard approach)
    effects = compute_effects_relative_to_mmt(temp_pcts, mmt, cb_coefs, cb_vcov, cb_info)
    
    # Store RR curve for visualization (subsample to reduce size)
    rr_curve_subset = rr_curve[::5]  # Every 5th point
    
    return {
        'region_code': int(region_code),
        'n_obs': int(valid.sum()),
        'n_years': float(n_years),
        'mean_deaths': float(y.mean()),
        'total_deaths': int(y.sum()),
        'dispersion': float(dispersion),
        'aic': float(model.aic),
        'temp_percentiles': temp_pcts,
        'mmt': float(mmt),
        'mmt_percentile': float(mmt_percentile_emp),
        'crossbasis_info': cb_info,
        'effects': effects,
        'rr_curve': rr_curve_subset,
        'cb_coefs': cb_coefs.tolist(),
        'cb_vcov': cb_vcov.tolist()
    }


# =============================================================================
# META-ANALYSIS (imported from utils.dlnm_module)
# =============================================================================
# random_effects_meta_analysis, pool_region_results, and 
# convert_to_json_serializable are imported from utils.dlnm_module


# =============================================================================
# DATA LOADING
# =============================================================================

print("\n" + "-"*70)
print("Loading Data")
print("-"*70)

era5_file = os.path.join(PHASE0_RESULTS, 'era5_intermediate_daily.parquet')
df_temp = pd.read_parquet(era5_file)
df_temp['date'] = pd.to_datetime(df_temp['date'])
print(f"ERA5: {len(df_temp):,} rows, {df_temp['region_code'].nunique()} regions")

mort_file = os.path.join(PHASE0_RESULTS, 'mortality_regional_daily_elderly.parquet')
df_mort = pd.read_parquet(mort_file)
df_mort['date'] = pd.to_datetime(df_mort['date'])
print(f"Mortality: {len(df_mort):,} rows")

ses_file = os.path.join(PHASE0_RESULTS, 'ses_intermediate_covariates.csv')
df_ses = pd.read_csv(ses_file)
pop_map = dict(zip(df_ses['intermediate_code'], df_ses['pop_elderly']))

holiday_file = os.path.join(PHASE0_RESULTS, 'brazilian_holidays_daily.parquet')
df_holidays = pd.read_parquet(holiday_file)
df_holidays['date'] = pd.to_datetime(df_holidays['date'])

# Merge data
print("\nMerging data...")
df = pd.merge(df_mort, df_temp[['date', 'region_code', 'temp_mean']], 
              on=['date', 'region_code'], how='inner')
df = pd.merge(df, df_holidays[['date', 'is_holiday']], on='date', how='left')
df['is_holiday'] = df['is_holiday'].fillna(0).astype(int)
df['pop_elderly'] = df['region_code'].map(pop_map)
df['year'] = df['date'].dt.year
df['month'] = df['date'].dt.month
df['day_of_week'] = df['date'].dt.dayofweek
df['time_index'] = (df['date'] - df['date'].min()).dt.days

df = df.dropna(subset=['temp_mean', 'deaths_elderly', 'pop_elderly'])
print(f"Final dataset: {len(df):,} rows, {df['region_code'].nunique()} regions")

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("FITTING REGION-SPECIFIC DLNM MODELS (Natural Spline Cross-Basis)")
print("="*70)

regions = sorted(df['region_code'].unique())
print(f"\nTotal regions: {len(regions)}")

region_results = {}
successful = 0

mmt_values = []

for i, region_code in enumerate(regions):
    df_region = df[df['region_code'] == region_code].copy()
    result = fit_region_dlnm(df_region, region_code)
    region_results[region_code] = result
    
    if result is not None:
        successful += 1
        mmt_values.append(result['mmt_percentile'])
        if (i + 1) % 20 == 0 or i < 5:
            eff_p99 = result['effects']['p99']
            eff_p1 = result['effects']['p1']
            print(f"  Region {region_code}: MMT=P{result['mmt_percentile']:.0f}, "
                  f"Heat RR(P99)={eff_p99['rr']:.3f}, Cold RR(P1)={eff_p1['rr']:.3f}")

print(f"\nSuccessfully fitted: {successful} / {len(regions)} regions")

# =============================================================================
# EXCLUSION DIAGNOSTICS & SANITY CHECKS
# =============================================================================

print("\n" + "-"*70)
print("EXCLUSION DIAGNOSTICS")
print("-"*70)

excluded_regions = [r for r, res in region_results.items() if res is None]
exclusion_rate = len(excluded_regions) / len(regions) * 100

print(f"\nExclusion Summary:")
print(f"  Total regions: {len(regions)}")
print(f"  Successful: {successful} ({100 - exclusion_rate:.1f}%)")
print(f"  Excluded: {len(excluded_regions)} ({exclusion_rate:.1f}%)")

# SANITY CHECK: Warn if exclusion rate is too high
if exclusion_rate > 20:
    print(f"\n  ⚠️  WARNING: High exclusion rate ({exclusion_rate:.1f}%)!")
    print(f"      Consider lowering MIN_OBS (currently {MIN_OBS})")
    print(f"      Or check data quality in excluded regions")
elif exclusion_rate > 10:
    print(f"\n  ⚠  CAUTION: Moderate exclusion rate ({exclusion_rate:.1f}%)")
else:
    print(f"\n  ✓ Exclusion rate acceptable ({exclusion_rate:.1f}%)")

# Analyze excluded regions
if excluded_regions:
    excluded_stats = []
    for r in excluded_regions:
        df_r = df[df['region_code'] == r]
        excluded_stats.append({
            'region': r,
            'n_obs': len(df_r),
            'total_deaths': df_r['deaths_elderly'].sum() if 'deaths_elderly' in df_r.columns else 0,
            'mean_deaths': df_r['deaths_elderly'].mean() if 'deaths_elderly' in df_r.columns else 0
        })
    
    exc_df = pd.DataFrame(excluded_stats)
    print(f"\nExcluded Region Characteristics:")
    print(f"  Obs range: {exc_df['n_obs'].min():.0f} - {exc_df['n_obs'].max():.0f} (threshold: {MIN_OBS})")
    print(f"  Total deaths range: {exc_df['total_deaths'].min():.0f} - {exc_df['total_deaths'].max():.0f}")
    print(f"  Mean deaths/day: {exc_df['mean_deaths'].mean():.2f}")
    
    # Check if any excluded regions have substantial data
    substantial = exc_df[exc_df['total_deaths'] > 500]
    if len(substantial) > 0:
        print(f"\n  ⚠️  {len(substantial)} excluded regions have >500 total deaths!")
        print(f"      These may have failed due to convergence, not sample size")
        print(f"      Region codes: {substantial['region'].tolist()[:5]}...")

# Coverage check: What % of total deaths are we capturing?
included_deaths = sum(r['total_deaths'] for r in region_results.values() if r is not None)
total_deaths = df['deaths_elderly'].sum()
coverage = included_deaths / total_deaths * 100 if total_deaths > 0 else 0

print(f"\nDeath Coverage:")
print(f"  Included: {included_deaths:,.0f} deaths ({coverage:.1f}%)")
print(f"  Excluded: {total_deaths - included_deaths:,.0f} deaths ({100-coverage:.1f}%)")

if coverage < 80:
    print(f"  ⚠️  WARNING: Capturing <80% of deaths - results may not be representative!")
elif coverage < 90:
    print(f"  ⚠  CAUTION: Capturing <90% of deaths")
else:
    print(f"  ✓ Good coverage (≥90% of deaths)")

# MMT distribution
if mmt_values:
    print(f"\nMMT Distribution across {len(mmt_values)} regions:")
    print(f"  Mean: P{np.mean(mmt_values):.1f}")
    print(f"  Median: P{np.median(mmt_values):.1f}")
    print(f"  Range: P{np.min(mmt_values):.0f} - P{np.max(mmt_values):.0f}")
    print(f"  (Expected for Brazil: P75-P85, warmer climates have higher MMT)")

# =============================================================================
# META-ANALYSIS POOLING (MVMeta - Gasparrini Methodology)
# =============================================================================

print("\n" + "="*70)
print("POOLING RESULTS VIA MULTIVARIATE META-ANALYSIS (MVMeta)")
print("="*70)
print("  Following Gasparrini 2015 Lancet methodology:")
print("  - Pool DLNM coefficients (not RRs at percentiles)")
print("  - Account for full covariance structure")
print("  - Find pooled MMT from pooled exposure-response curve")

# MVMeta pooling of coefficients
mvmeta_results = mvmeta_pool_regions(region_results, min_obs=MIN_OBS)

if 'error' not in mvmeta_results:
    print(f"\n  Regions included: {mvmeta_results['n_regions']}")
    print(f"  Pooled MMT: {mvmeta_results['mmt']:.1f}°C (P{mvmeta_results['mmt_percentile']:.0f})")
    
    print("\nHEAT EFFECTS (MVMeta pooled, vs pooled MMT)")
    print("-"*50)
    for pct in ['p75', 'p90', 'p95', 'p97.5', 'p99']:
        if pct in mvmeta_results['effects']:
            eff = mvmeta_results['effects'][pct]
            print(f"  {pct.upper()}: RR = {eff['rr']:.3f} "
                  f"[{eff['rr_lower']:.3f}-{eff['rr_upper']:.3f}]")
    
    print("\nCOLD EFFECTS (MVMeta pooled, vs pooled MMT)")
    print("-"*50)
    for pct in ['p25', 'p10', 'p5', 'p2.5', 'p1']:
        if pct in mvmeta_results['effects']:
            eff = mvmeta_results['effects'][pct]
            print(f"  {pct.upper()}: RR = {eff['rr']:.3f} "
                  f"[{eff['rr_lower']:.3f}-{eff['rr_upper']:.3f}]")
else:
    print(f"\n  MVMeta error: {mvmeta_results.get('error')}")

# Also compute simple RR pooling for comparison
print("\n" + "-"*50)
print("COMPARISON: Simple RR Pooling (for reference)")
print("-"*50)

pooled_results = {}
for pct in ['p1', 'p2.5', 'p5', 'p10', 'p25', 'p50', 'p75', 'p90', 'p95', 'p97.5', 'p99']:
    pooled_results[pct] = pool_region_results(region_results, pct)

for pct in ['p99', 'p1']:
    p = pooled_results[pct]
    if 'pooled_rr' in p and not np.isnan(p.get('pooled_rr', np.nan)):
        print(f"  {pct.upper()}: RR = {p['pooled_rr']:.3f} "
              f"[{p['pooled_rr_lower']:.3f}-{p['pooled_rr_upper']:.3f}], I² = {p['I2']:.1f}%")

# =============================================================================
# SAVE RESULTS
# =============================================================================

print("\n" + "="*70)
print("SAVING RESULTS")
print("="*70)

# Compute pooled MMT statistics
mmt_stats = {
    'mean_percentile': float(np.mean(mmt_values)) if mmt_values else None,
    'median_percentile': float(np.median(mmt_values)) if mmt_values else None,
    'min_percentile': float(np.min(mmt_values)) if mmt_values else None,
    'max_percentile': float(np.max(mmt_values)) if mmt_values else None,
    'sd_percentile': float(np.std(mmt_values)) if mmt_values else None
}

full_results = {
    'analysis_level': 'intermediate',
    'method': 'DLNM with natural cubic spline cross-basis',
    'pooling_method': 'MVMeta (Multivariate Meta-Analysis of coefficients)',
    'reference': 'Pooled MMT from MVMeta curve',
    'n_regions': len(regions),
    'n_successful': successful,
    'mmt_statistics': mmt_stats,
    'parameters': {
        'max_lag': MAX_LAG,
        'temp_df': TEMP_DF,
        'lag_df': LAG_DF,
        'temp_knot_pcts': TEMP_KNOT_PCTS,
        'lag_knots': LAG_KNOTS,
        'time_spline_df_per_year': TIME_SPLINE_DF_PER_YEAR,
        'min_obs': MIN_OBS
    },
    'region_results': {str(k): v for k, v in region_results.items() if v is not None},
    'mvmeta_results': mvmeta_results if 'error' not in mvmeta_results else None,
    'simple_pooled_results': pooled_results,  # Keep for comparison
    'timestamp': datetime.now().isoformat()
}

full_results = convert_to_json_serializable(full_results)

output_file = os.path.join(OUTPUT_DIR, 'dlnm_v2_intermediate_results.json')
with open(output_file, 'w') as f:
    json.dump(full_results, f, indent=2)
print(f"Saved: {output_file}")

# Summary CSV
summary_rows = []
for region_code, result in region_results.items():
    if result is None:
        continue
    row = {
        'region_code': region_code,
        'n_obs': result['n_obs'],
        'n_years': result['n_years'],
        'mean_deaths': result['mean_deaths'],
        'dispersion': result['dispersion']
    }
    for pct in ['p1', 'p5', 'p50', 'p95', 'p99']:
        if pct in result['effects']:
            eff = result['effects'][pct]
            row[f'rr_{pct}'] = eff['rr']
            row[f'rr_{pct}_lower'] = eff['rr_lower']
            row[f'rr_{pct}_upper'] = eff['rr_upper']

            # Wald-style p-value from log-RR and its standard error
            log_rr = eff.get('log_rr')
            se = eff.get('log_rr_se')
            if se is not None and se > 0:
                z = log_rr / se
                p_value = 2 * (1 - stats.norm.cdf(abs(z)))
            else:
                p_value = np.nan
            row[f'rr_{pct}_p'] = p_value
    summary_rows.append(row)

summary_df = pd.DataFrame(summary_rows)
summary_file = os.path.join(OUTPUT_DIR, 'dlnm_v2_intermediate_summary.csv')
summary_df.to_csv(summary_file, index=False)
print(f"Saved: {summary_file}")

pooled_file = os.path.join(OUTPUT_DIR, 'dlnm_v2_intermediate_pooled.json')
with open(pooled_file, 'w') as f:
    json.dump(convert_to_json_serializable(pooled_results), f, indent=2)
print(f"Saved: {pooled_file}")

# Save MVMeta results separately
if 'error' not in mvmeta_results:
    mvmeta_file = os.path.join(OUTPUT_DIR, 'dlnm_v2_intermediate_mvmeta.json')
    # Don't save the full coefficient arrays to reduce file size
    mvmeta_save = {
        'method': mvmeta_results['method'],
        'n_regions': mvmeta_results['n_regions'],
        'mmt': mvmeta_results['mmt'],
        'mmt_percentile': mvmeta_results['mmt_percentile'],
        'global_percentiles': mvmeta_results['global_percentiles'],
        'effects': mvmeta_results['effects'],
    }
    with open(mvmeta_file, 'w') as f:
        json.dump(convert_to_json_serializable(mvmeta_save), f, indent=2)
    print(f"Saved: {mvmeta_file}")

# =============================================================================
# COMPARISON WITH V1 (if available)
# =============================================================================

v1_file = os.path.join(OUTPUT_DIR, 'dlnm_intermediate_pooled.json')
if os.path.exists(v1_file):
    print("\n" + "="*70)
    print("COMPARISON: V1 (Polynomial) vs V2 (Natural Spline)")
    print("="*70)
    
    with open(v1_file, 'r') as f:
        v1_pooled = json.load(f)
    
    print(f"\n{'Percentile':<12} {'V1 (Poly)':<20} {'V2 (NS)':<20} {'Diff':<10}")
    print("-"*62)
    
    for pct in ['p1', 'p5', 'p95', 'p99']:
        v1_rr = v1_pooled.get(pct, {}).get('pooled_rr', np.nan)
        v2_rr = pooled_results.get(pct, {}).get('pooled_rr', np.nan)
        
        if not np.isnan(v1_rr) and not np.isnan(v2_rr):
            diff = (v2_rr - v1_rr) / v1_rr * 100
            print(f"{pct:<12} {v1_rr:.3f}{'':<15} {v2_rr:.3f}{'':<15} {diff:+.1f}%")

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
