"""
04a: META-REGRESSION (REGIONAL DLNM) - v2
==========================================
Fits region-specific DLNMs and uses meta-regression to explain heterogeneity
in heat/cold effects across Brazil's intermediate regions.

v2 IMPROVEMENTS (Dec 2025):
---------------------------
1. Uses natural spline cross-basis from utils (not polynomial)
2. Population offset in all Poisson models
3. Full covariance delta method for CIs
4. Region-specific percentiles for thresholds
5. Proper DerSimonian-Laird meta-analysis with tau²
6. Named parameter lookup (not positional indexing)

Predictors tested:
- AC ownership (heat adaptation)
- HDI / GDP per capita (socioeconomic)
- Climate zone / mean temperature
- Urbanization
- Healthcare capacity (physicians, hospital beds)

Two-stage approach:
1. Fit DLNM in each of 133 regions → extract heat/cold RRs
2. Meta-regression: test which covariates explain RR variation
"""

import warnings
warnings.filterwarnings('ignore')

import sys
import json
import pandas as pd
import numpy as np
import statsmodels.api as sm
from datetime import datetime
from pathlib import Path
from scipy import stats
from typing import Dict, Any, Optional, List, Tuple
import argparse

# Add parent directory to path for utils import
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.dlnm_module import (
    fit_region_dlnm,
    predict_cumulative_rr,
    meta_random_effects,
    pool_region_effects,
    convert_to_json_serializable,
    mvmeta_pool_coefficients,
    compute_pooled_rr_from_mvmeta,
    find_mmt,
)

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = Path(__file__).parent.parent
PHASE0_RESULTS = BASE_DIR / 'phase0_data_prep' / 'results'
PHASE1_RESULTS = BASE_DIR / 'phase1_core_model' / 'results'
OUTPUT_DIR = Path(__file__).parent / 'results'
OUTPUT_DIR.mkdir(exist_ok=True)

# DLNM parameters (consistent with Phase 1 v2)
MAX_LAG = 21
VAR_DF = 4  # Natural spline df for temperature
LAG_DF = 4  # Natural spline df for lag

# Minimum observations per region
MIN_OBS = 500

parser = argparse.ArgumentParser(description='Meta-regression of regional DLNM heat/cold effects')
parser.add_argument('--level', type=str, default='intermediate',
                    choices=['intermediate', 'immediate'],
                    help='Spatial level to analyze')
args, _ = parser.parse_known_args()

level = args.level
if level == 'intermediate':
    temp_filename = 'era5_intermediate_daily.parquet'
    mort_filename = 'mortality_regional_daily_elderly.parquet'
    covariates_filename = 'regional_covariates.csv'
    cov_region_col = 'region_code'
    level_label = 'INTERMEDIATE (133 regions)'
    output_suffix = ''
else:
    temp_filename = 'era5_immediate_daily.parquet'
    mort_filename = 'mortality_immediate_daily_elderly.parquet'
    covariates_filename = 'ses_immediate_covariates.csv'
    cov_region_col = 'immediate_code'
    level_label = 'IMMEDIATE (510 regions)'
    output_suffix = '_immediate'

print("="*70)
print("04a: META-REGRESSION (REGIONAL DLNM) - v2")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Using: Natural spline cross-basis (var_df={VAR_DF}, lag_df={LAG_DF})")
print(f"Level: {level} [{level_label}]")

# =============================================================================
# 1. LOAD DATA
# =============================================================================
print("\n[1] Loading data...")

# Temperature data
temp_df = pd.read_parquet(PHASE0_RESULTS / temp_filename)
temp_df['date'] = pd.to_datetime(temp_df['date'])
print(f"  Temperature records: {len(temp_df):,}")

# Mortality data (with population)
mort_df = pd.read_parquet(PHASE0_RESULTS / mort_filename)
mort_df['date'] = pd.to_datetime(mort_df['date'])
print(f"  Mortality records: {len(mort_df):,}")

# Regional covariates
covariates = pd.read_csv(PHASE0_RESULTS / covariates_filename)
print(f"  Regions with covariates: {len(covariates)}")
print(f"  Covariates: {list(covariates.columns)[:10]}...")

# Standardize region column names for merging
if level == 'immediate':
    mort_df = mort_df.rename(columns={'immediate_code': 'region_code'})
elif level == 'intermediate':
    if 'intermediate_code' in mort_df.columns:
        mort_df = mort_df.rename(columns={'intermediate_code': 'region_code'})

# Merge data
df = mort_df.merge(temp_df, on=['date', 'region_code'], how='inner')
print(f"  Merged dataset: {len(df):,} observations")

# Ensure population column exists
if 'pop_elderly' not in df.columns:
    # Try to get from covariates
    if 'pop_elderly' in covariates.columns:
        pop_map = covariates.set_index(cov_region_col)['pop_elderly'].to_dict()
        df['pop_elderly'] = df['region_code'].map(pop_map)
        print("  Added pop_elderly from covariates table")
    else:
        # Try to load from SES file
        if level == 'intermediate':
            ses_file = PHASE0_RESULTS / 'ses_intermediate_covariates.csv'
            ses_col = 'intermediate_code'
        else:
            ses_file = PHASE0_RESULTS / 'ses_immediate_covariates.csv'
            ses_col = 'immediate_code'
        
        if ses_file.exists():
            ses_df = pd.read_csv(ses_file)
            if 'pop_elderly' in ses_df.columns:
                pop_map = dict(zip(ses_df[ses_col], ses_df['pop_elderly']))
                df['pop_elderly'] = df['region_code'].map(pop_map)
                print(f"  Added pop_elderly from SES file ({len(pop_map)} regions)")
            else:
                print("  WARNING: No population column found - will use no offset")
                df['pop_elderly'] = 1  # Placeholder
        else:
            print("  WARNING: No population column found - will use no offset")
            df['pop_elderly'] = 1  # Placeholder

# Add time variables needed by fit_region_dlnm
df['dow'] = df['date'].dt.dayofweek
df['month'] = df['date'].dt.month
df['year'] = df['date'].dt.year

# Count regions
regions = df['region_code'].unique()
print(f"  Unique regions: {len(regions)}")

# =============================================================================
# 2. COMPUTE REGION-SPECIFIC PERCENTILES
# =============================================================================
print("\n[2] Computing region-specific temperature percentiles...")

region_percentiles = {}
for region in regions:
    temps = df.loc[df['region_code'] == region, 'temp_mean']
    if len(temps) < 100:
        continue
    region_percentiles[region] = {
        'p1': float(temps.quantile(0.01)),
        'p2_5': float(temps.quantile(0.025)),
        'p50': float(temps.quantile(0.50)),
        'p97_5': float(temps.quantile(0.975)),
        'p99': float(temps.quantile(0.99)),
        'mean': float(temps.mean()),
        'std': float(temps.std()),
    }

print(f"  Computed percentiles for {len(region_percentiles)} regions")

# Also compute global for comparison
global_temps = df['temp_mean']
global_percentiles = {
    'p1': float(global_temps.quantile(0.01)),
    'p2_5': float(global_temps.quantile(0.025)),
    'p50': float(global_temps.quantile(0.50)),
    'p97_5': float(global_temps.quantile(0.975)),
    'p99': float(global_temps.quantile(0.99)),
}
print(f"  Global: P1={global_percentiles['p1']:.1f}°C, P2.5={global_percentiles['p2_5']:.1f}°C, "
      f"P50={global_percentiles['p50']:.1f}°C, P97.5={global_percentiles['p97_5']:.1f}°C, "
      f"P99={global_percentiles['p99']:.1f}°C")

# =============================================================================
# 3. FIT REGION-SPECIFIC DLNMS
# =============================================================================
print("\n[3] Fitting region-specific DLNMs...")

def fit_and_extract_rr(df_region: pd.DataFrame, pctiles: Dict) -> Optional[Dict]:
    """
    Fit DLNM for a region and extract heat/cold RRs.
    Uses utils/dlnm_module functions for proper implementation.
    
    Note: Meta-regression uses P50 as reference to compare effects across regions.
    This is different from MMT-based analysis in Phase 1, but appropriate for
    explaining heterogeneity in effects.
    """
    try:
        # Fit DLNM using shared utilities
        fit_res = fit_region_dlnm(
            df_region=df_region,
            temp_col='temp_mean',
            deaths_col='deaths_elderly',
            pop_col='pop_elderly',  # Population offset!
            max_lag=MAX_LAG,
            var_df=VAR_DF,
            lag_df=LAG_DF,
            family='quasi-poisson',
            min_obs=MIN_OBS,
        )
        
        if fit_res is None:
            return None
        
        # Extract RR at region-specific percentiles vs P50 (reference)
        ref_temp = pctiles['p50']
        
        # Heat effects (P97.5, P99)
        heat_97_5 = predict_cumulative_rr(fit_res, pctiles['p97_5'], ref_temp)
        heat_99 = predict_cumulative_rr(fit_res, pctiles['p99'], ref_temp)
        
        # Cold effects (P2.5, P1)
        cold_2_5 = predict_cumulative_rr(fit_res, pctiles['p2_5'], ref_temp)
        cold_1 = predict_cumulative_rr(fit_res, pctiles['p1'], ref_temp)
        
        # Derive SEs consistently from confidence intervals
        def se_from_ci(rr, lo, hi):
            if rr <= 0 or lo <= 0 or hi <= 0 or np.isnan(rr) or np.isnan(lo) or np.isnan(hi):
                return np.nan
            return (np.log(hi) - np.log(lo)) / (2 * 1.96)

        heat_se_97_5 = se_from_ci(heat_97_5[0], heat_97_5[1], heat_97_5[2])
        cold_se_2_5 = se_from_ci(cold_2_5[0], cold_2_5[1], cold_2_5[2])
        heat_se_99 = se_from_ci(heat_99[0], heat_99[1], heat_99[2])
        cold_se_1 = se_from_ci(cold_1[0], cold_1[1], cold_1[2])
        
        # Extract cross-basis coefficients for potential MVMeta pooling
        cb_colnames = fit_res.get('cb_colnames', [])
        params = fit_res.get('params', pd.Series())
        cov = fit_res.get('cov', pd.DataFrame())
        
        if len(cb_colnames) > 0 and all(c in params.index for c in cb_colnames):
            cb_coefs = params[cb_colnames].values
            cb_vcov = cov.loc[cb_colnames, cb_colnames].values
        else:
            cb_coefs = np.array([])
            cb_vcov = np.array([[]])

        return {
            'n_obs': fit_res['n_obs'],
            'dispersion': fit_res.get('dispersion', np.nan),
            
            # Primary thresholds (P97.5/P2.5 vs P50)
            'heat_rr_97_5': heat_97_5[0],
            'heat_rr_97_5_lo': heat_97_5[1],
            'heat_rr_97_5_hi': heat_97_5[2],
            'heat_log_rr_97_5': heat_97_5[3],
            'heat_se_97_5': heat_se_97_5,
            
            'cold_rr_2_5': cold_2_5[0],
            'cold_rr_2_5_lo': cold_2_5[1],
            'cold_rr_2_5_hi': cold_2_5[2],
            'cold_log_rr_2_5': cold_2_5[3],
            'cold_se_2_5': cold_se_2_5,
            
            # Comparison thresholds (P99/P1 vs P50)
            'heat_rr_99': heat_99[0],
            'heat_rr_99_lo': heat_99[1],
            'heat_rr_99_hi': heat_99[2],
            'heat_log_rr_99': heat_99[3],
            'heat_se_99': heat_se_99,
            
            'cold_rr_1': cold_1[0],
            'cold_rr_1_lo': cold_1[1],
            'cold_rr_1_hi': cold_1[2],
            'cold_log_rr_1': cold_1[3],
            'cold_se_1': cold_se_1,
            
            # Temperature info
            'temp_p50': pctiles['p50'],
            'temp_p97_5': pctiles['p97_5'],
            'temp_p99': pctiles['p99'],
            'temp_p2_5': pctiles['p2_5'],
            'temp_p1': pctiles['p1'],
            
            # For potential MVMeta
            'cb_coefs': cb_coefs.tolist() if len(cb_coefs) > 0 else None,
            'cb_vcov': cb_vcov.tolist() if cb_vcov.size > 0 else None,
            'cb_meta': fit_res.get('cb_meta'),
        }
        
    except Exception as e:
        return None

# Fit each region
region_results = {}
successful = 0
failed_regions = []

for i, region in enumerate(regions):
    if (i + 1) % 20 == 0 or i == 0:
        print(f"  Processing region {i+1}/{len(regions)}...")
    
    if region not in region_percentiles:
        failed_regions.append((region, "no_percentiles"))
        continue
    
    df_region = df[df['region_code'] == region].copy()
    pctiles = region_percentiles[region]
    
    result = fit_and_extract_rr(df_region, pctiles)
    
    if result is not None:
        region_results[int(region)] = result
        successful += 1
    else:
        failed_regions.append((region, "fit_failed"))

print(f"\n  Successfully fit: {successful}/{len(regions)} regions")
if failed_regions:
    print(f"  Failed regions: {len(failed_regions)}")

# =============================================================================
# 4. META-ANALYSIS (POOLED EFFECTS)
# =============================================================================
print("\n[4] Meta-analysis of regional effects...")

def compute_log_rr_se(result: Dict, threshold: str) -> tuple:
    """Extract log-RR and SE for a threshold from result dict."""
    rr = result.get(f'{threshold}_rr_{threshold.split("_")[1]}', np.nan)
    if np.isnan(rr) or rr <= 0:
        return np.nan, np.nan
    log_rr = np.log(rr)
    
    # Compute SE from CI width
    lo = result.get(f'{threshold}_rr_{threshold.split("_")[1]}_lo', np.nan)
    hi = result.get(f'{threshold}_rr_{threshold.split("_")[1]}_hi', np.nan)
    if not np.isnan(lo) and not np.isnan(hi) and lo > 0 and hi > 0:
        se = (np.log(hi) - np.log(lo)) / (2 * 1.96)
    else:
        se = np.nan
    
    return log_rr, se

# Pool heat effects (P97.5)
heat_effects = []
heat_variances = []
for region, res in region_results.items():
    log_rr = np.log(res['heat_rr_97_5']) if res['heat_rr_97_5'] > 0 else np.nan
    lo, hi = res['heat_rr_97_5_lo'], res['heat_rr_97_5_hi']
    if not np.isnan(log_rr) and lo > 0 and hi > 0:
        se = (np.log(hi) - np.log(lo)) / (2 * 1.96)
        heat_effects.append(log_rr)
        heat_variances.append(se ** 2)

heat_meta = meta_random_effects(np.array(heat_effects), np.array(heat_variances))
print(f"\n  Heat (P97.5 vs P50):")
print(f"    Pooled RR: {np.exp(heat_meta['effect']):.3f} "
      f"[{np.exp(heat_meta['effect'] - 1.96*heat_meta['se']):.3f}-"
      f"{np.exp(heat_meta['effect'] + 1.96*heat_meta['se']):.3f}]")
print(f"    I²: {heat_meta['I2']:.1f}%, tau²: {heat_meta['tau2']:.4f}, k={heat_meta['k']}")

# Pool cold effects (P2.5)
cold_effects = []
cold_variances = []
for region, res in region_results.items():
    log_rr = np.log(res['cold_rr_2_5']) if res['cold_rr_2_5'] > 0 else np.nan
    lo, hi = res['cold_rr_2_5_lo'], res['cold_rr_2_5_hi']
    if not np.isnan(log_rr) and lo > 0 and hi > 0:
        se = (np.log(hi) - np.log(lo)) / (2 * 1.96)
        cold_effects.append(log_rr)
        cold_variances.append(se ** 2)

cold_meta = meta_random_effects(np.array(cold_effects), np.array(cold_variances))
print(f"\n  Cold (P2.5 vs P50):")
print(f"    Pooled RR: {np.exp(cold_meta['effect']):.3f} "
      f"[{np.exp(cold_meta['effect'] - 1.96*cold_meta['se']):.3f}-"
      f"{np.exp(cold_meta['effect'] + 1.96*cold_meta['se']):.3f}]")
print(f"    I²: {cold_meta['I2']:.1f}%, tau²: {cold_meta['tau2']:.4f}, k={cold_meta['k']}")

# =============================================================================
# 5. META-REGRESSION
# =============================================================================
print("\n[5] Meta-regression analysis...")

# Prepare data for meta-regression
meta_data = []
for region_code, result in region_results.items():
    cov_row = covariates[covariates[cov_region_col] == region_code] if cov_region_col in covariates.columns else covariates[covariates['region_code'] == region_code]
    if len(cov_row) == 0:
        continue
    
    cov_row = cov_row.iloc[0]
    
    # Extract heat log-RR and SE
    heat_log_rr = np.log(result['heat_rr_97_5']) if result['heat_rr_97_5'] > 0 else np.nan
    heat_lo, heat_hi = result['heat_rr_97_5_lo'], result['heat_rr_97_5_hi']
    if heat_lo > 0 and heat_hi > 0:
        heat_se = (np.log(heat_hi) - np.log(heat_lo)) / (2 * 1.96)
    else:
        heat_se = np.nan
    
    # Extract cold log-RR and SE
    cold_log_rr = np.log(result['cold_rr_2_5']) if result['cold_rr_2_5'] > 0 else np.nan
    cold_lo, cold_hi = result['cold_rr_2_5_lo'], result['cold_rr_2_5_hi']
    if cold_lo > 0 and cold_hi > 0:
        cold_se = (np.log(cold_hi) - np.log(cold_lo)) / (2 * 1.96)
    else:
        cold_se = np.nan
    
    meta_data.append({
        'region_code': region_code,
        'region_name': cov_row.get('region_name', ''),
        'heat_log_rr': heat_log_rr,
        'heat_se': heat_se,
        'cold_log_rr': cold_log_rr,
        'cold_se': cold_se,
        'n_obs': result['n_obs'],
        'temp_p50': result['temp_p50'],
        'ac_pct': cov_row.get('ac_pct', np.nan),
        'hdi': cov_row.get('hdi', np.nan),
        # GDP per capita: regional table uses gdp_per_capita_brl, SES immediate uses gdp_per_capita
        'gdp_per_capita': cov_row.get('gdp_per_capita_brl', cov_row.get('gdp_per_capita', np.nan)),
        # Urbanization: regional table has urban_pct, SES immediate has urbanization_rate
        'urban_pct': cov_row.get('urban_pct', cov_row.get('urbanization_rate', np.nan)),
        # Elderly share: regional has elderly_pct, SES immediate has pct_elderly
        'elderly_pct': cov_row.get('elderly_pct', cov_row.get('pct_elderly', np.nan)),
        'hospital_beds': cov_row.get('hospital_beds_per_1000', np.nan),
        'physicians': cov_row.get('physicians_per_1000', np.nan),
        'mean_temp': cov_row.get('mean_temp_annual', np.nan),
        'macro_region': cov_row.get('macro_region', ''),
    })

meta_df = pd.DataFrame(meta_data)

# Drop rows with missing key variables
meta_df_clean = meta_df.dropna(subset=['heat_log_rr', 'heat_se', 'cold_log_rr', 'cold_se'])
print(f"  Regions in meta-regression: {len(meta_df_clean)}")

# Standardize predictors
predictors = ['ac_pct', 'hdi', 'urban_pct', 'elderly_pct', 'mean_temp', 'hospital_beds']
for pred in predictors:
    if pred in meta_df_clean.columns:
        col = meta_df_clean[pred]
        if col.notna().sum() > 0:
            meta_df_clean[f'{pred}_z'] = (col - col.mean()) / col.std()
        else:
            meta_df_clean[f'{pred}_z'] = np.nan

# Run meta-regressions with proper random-effects variance
meta_results = {'heat': {}, 'cold': {}}

for outcome in ['heat', 'cold']:
    print(f"\n  === {outcome.upper()} META-REGRESSION ===")
    
    y = meta_df_clean[f'{outcome}_log_rr'].values
    within_var = meta_df_clean[f'{outcome}_se'].values ** 2
    
    # Add between-study variance (tau²) to weights
    # First, estimate tau² from the data
    k = len(y)
    w_fixed = 1.0 / within_var
    fixed_mean = np.sum(w_fixed * y) / np.sum(w_fixed)
    Q = np.sum(w_fixed * (y - fixed_mean) ** 2)
    c = np.sum(w_fixed) - (np.sum(w_fixed ** 2) / np.sum(w_fixed))
    tau2 = max(0.0, (Q - (k - 1)) / c) if c > 0 else 0.0
    
    # Random-effects weights
    weights = 1.0 / (within_var + tau2)
    
    # Pooled estimate
    pooled_mean = np.sum(weights * y) / np.sum(weights)
    pooled_se = np.sqrt(1.0 / np.sum(weights))
    pooled_rr = np.exp(pooled_mean)
    pooled_rr_lo = np.exp(pooled_mean - 1.96 * pooled_se)
    pooled_rr_hi = np.exp(pooled_mean + 1.96 * pooled_se)
    
    # I² 
    I2 = max(0, (Q - (k-1)) / Q * 100) if Q > 0 else 0
    
    print(f"    Pooled RR: {pooled_rr:.3f} [{pooled_rr_lo:.3f}-{pooled_rr_hi:.3f}]")
    print(f"    I²: {I2:.1f}%, tau²: {tau2:.4f}")
    
    meta_results[outcome]['pooled'] = {
        'rr': float(pooled_rr),
        'rr_lo': float(pooled_rr_lo),
        'rr_hi': float(pooled_rr_hi),
        'log_rr': float(pooled_mean),
        'se': float(pooled_se),
        'i_squared': float(I2),
        'tau_squared': float(tau2),
        'n_regions': k,
    }
    
    # Meta-regression for each predictor (using RE weights)
    meta_results[outcome]['predictors'] = {}
    
    for pred in predictors:
        pred_z = f'{pred}_z'
        if pred_z not in meta_df_clean.columns:
            continue
        
        x = meta_df_clean[pred_z].values
        valid = ~np.isnan(x)
        if valid.sum() < 10:
            continue
        
        X = sm.add_constant(x[valid])
        
        try:
            model = sm.WLS(y[valid], X, weights=weights[valid])
            result = model.fit()
            
            coef = result.params[1]
            se = result.bse[1]
            pval = result.pvalues[1]
            
            # Effect per 1 SD increase in predictor
            rr_ratio = np.exp(coef)
            rr_ratio_lo = np.exp(coef - 1.96 * se)
            rr_ratio_hi = np.exp(coef + 1.96 * se)
            
            sig_marker = "*" if pval < 0.05 else ""
            print(f"    {pred}: RR ratio = {rr_ratio:.3f} [{rr_ratio_lo:.3f}-{rr_ratio_hi:.3f}], "
                  f"p = {pval:.3f}{sig_marker}")
            
            meta_results[outcome]['predictors'][pred] = {
                'coef': float(coef),
                'se': float(se),
                'rr_ratio_per_sd': float(rr_ratio),
                'rr_ratio_lo': float(rr_ratio_lo),
                'rr_ratio_hi': float(rr_ratio_hi),
                'pvalue': float(pval),
                'significant': bool(pval < 0.05),
                'n_regions': int(valid.sum()),
            }
        except Exception as e:
            print(f"    {pred}: Failed - {str(e)[:50]}")

# =============================================================================
# 6. MULTIVARIATE META-REGRESSION
# =============================================================================
print("\n[6] Multivariate meta-regression...")

for outcome in ['heat', 'cold']:
    print(f"\n  === {outcome.upper()} MULTIVARIATE ===")
    
    y = meta_df_clean[f'{outcome}_log_rr'].values
    within_var = meta_df_clean[f'{outcome}_se'].values ** 2
    
    # Estimate tau² for RE weights
    w_fixed = 1.0 / within_var
    fixed_mean = np.sum(w_fixed * y) / np.sum(w_fixed)
    Q = np.sum(w_fixed * (y - fixed_mean) ** 2)
    c = np.sum(w_fixed) - (np.sum(w_fixed ** 2) / np.sum(w_fixed))
    tau2 = max(0.0, (Q - (len(y) - 1)) / c) if c > 0 else 0.0
    weights = 1.0 / (within_var + tau2)
    
    # Key predictors for multivariate model
    X_vars = ['ac_pct_z', 'mean_temp_z', 'urban_pct_z']
    X_vars_present = [v for v in X_vars if v in meta_df_clean.columns]
    
    if len(X_vars_present) < 2:
        print("    Insufficient predictors for multivariate model")
        continue
    
    X_data = meta_df_clean[X_vars_present].values
    valid = ~np.isnan(X_data).any(axis=1)
    
    if valid.sum() < 20:
        print(f"    Insufficient observations: {valid.sum()}")
        continue
    
    X = sm.add_constant(X_data[valid])
    
    try:
        model = sm.WLS(y[valid], X, weights=weights[valid])
        result = model.fit()
        
        print(f"    R² = {result.rsquared:.3f}")
        print(f"    N = {valid.sum()}")
        
        meta_results[outcome]['multivariate'] = {
            'r_squared': float(result.rsquared),
            'n_regions': int(valid.sum()),
            'predictors': {},
        }
        
        for i, var in enumerate(X_vars_present):
            coef = result.params[i+1]
            se = result.bse[i+1]
            pval = result.pvalues[i+1]
            
            sig_marker = "*" if pval < 0.05 else ""
            print(f"    {var}: β = {coef:.4f} (SE={se:.4f}), p = {pval:.3f}{sig_marker}")
            
            meta_results[outcome]['multivariate']['predictors'][var] = {
                'coef': float(coef),
                'se': float(se),
                'pvalue': float(pval),
                'rr_ratio_per_sd': float(np.exp(coef)),
            }
    except Exception as e:
        print(f"    Failed: {e}")

# =============================================================================
# 7. SAVE RESULTS
# =============================================================================
print("\n[7] Saving results...")

# Summary statistics
summary = {
    'analysis_date': datetime.now().isoformat(),
    'level': level,
    'n_regions_total': len(regions),
    'n_regions_fitted': successful,
    'n_regions_meta': len(meta_df_clean),
    'dlnm_params': {
        'max_lag': MAX_LAG,
        'var_df': VAR_DF,
        'lag_df': LAG_DF,
        'basis_type': 'natural_spline',
    },
    'percentile_method': 'region_specific',
    'offset': 'log(pop_elderly)',
    'meta_results': meta_results,
    'global_percentiles': global_percentiles,
}

# Region-level results
region_output = {
    int(k): convert_to_json_serializable(v) 
    for k, v in region_results.items()
}

# Save
with open(OUTPUT_DIR / f'meta_regression_v2_results{output_suffix}.json', 'w') as f:
    json.dump(convert_to_json_serializable(summary), f, indent=2)
    
with open(OUTPUT_DIR / f'meta_regression_v2_region_effects{output_suffix}.json', 'w') as f:
    json.dump(region_output, f, indent=2)

meta_df_clean.to_csv(OUTPUT_DIR / f'meta_regression_v2_data{output_suffix}.csv', index=False)

print(f"\n  Saved to {OUTPUT_DIR}:")
print(f"    - meta_regression_v2_results{output_suffix}.json")
print(f"    - meta_regression_v2_region_effects{output_suffix}.json")
print(f"    - meta_regression_v2_data{output_suffix}.csv")

print(f"\n{'='*70}")
print("04a: META-REGRESSION COMPLETE")
print(f"{'='*70}")
print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
