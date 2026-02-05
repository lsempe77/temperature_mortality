"""
04d: CAUSE-SPECIFIC STRATIFICATION (REGIONAL DLNM) - v2
========================================================
Tests effect modification by cause of death:
- All-cause
- Cardiovascular (CVD)
- Respiratory
- Other causes

v2 IMPROVEMENTS (Dec 2025):
---------------------------
1. CRITICAL FIX: Keep zero-count days (v1 incorrectly filtered them out!)
2. Per-region DLNM fits then meta-analysis (not pooled model with dummies)
3. Uses natural spline cross-basis from utils
4. Population offset in all models
5. Full covariance delta method for CIs
6. Region-specific percentiles
7. Proper heterogeneity testing between causes

Two-stage approach (per cause):
1. Fit DLNM in each region → extract heat/cold RRs
2. Meta-analyze across regions (DerSimonian-Laird)
3. Compare pooled effects between causes
"""

import warnings
warnings.filterwarnings('ignore')

import sys
import json
import pandas as pd
import numpy as np
from datetime import datetime
from pathlib import Path
from scipy import stats
from typing import Dict, Any, Optional, Tuple
import argparse

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.dlnm_module import (
    fit_region_dlnm,
    predict_cumulative_rr,
    find_mmt,
    meta_random_effects,
    convert_to_json_serializable,
    mvmeta_pool_coefficients,
    compute_pooled_rr_from_mvmeta,
)

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = Path(__file__).parent.parent
PHASE0_RESULTS = BASE_DIR / 'phase0_data_prep' / 'results'
OUTPUT_DIR = Path(__file__).parent / 'results'
OUTPUT_DIR.mkdir(exist_ok=True)

MAX_LAG = 21
VAR_DF = 4
LAG_DF = 4
MIN_OBS = 365  # Per region-cause

parser = argparse.ArgumentParser(description='Cause-specific regional DLNM meta-analysis (v2)')
parser.add_argument('--level', type=str, default='intermediate',
                    choices=['intermediate', 'immediate'],
                    help='Spatial level to analyze')
args, _ = parser.parse_known_args()

level = args.level
if level == 'intermediate':
    era5_filename = 'era5_intermediate_daily.parquet'
    mort_filename = 'mortality_regional_daily_elderly.parquet'
    pop_covariates_filename = 'regional_covariates.csv'
    pop_region_col = 'region_code'
    level_label = 'INTERMEDIATE (133 regions)'
else:
    era5_filename = 'era5_immediate_daily.parquet'
    mort_filename = 'mortality_immediate_daily_elderly.parquet'
    pop_covariates_filename = 'ses_immediate_covariates.csv'
    pop_region_col = 'immediate_code'
    level_label = 'IMMEDIATE (510 regions)'

# Cause categories
CAUSES = {
    'all_cause': {
        'col': 'deaths_elderly',
        'label': 'All-cause',
    },
    'cardiovascular': {
        'col': 'deaths_elderly_cvd',
        'label': 'Cardiovascular (I00-I99)',
    },
    'respiratory': {
        'col': 'deaths_elderly_resp',
        'label': 'Respiratory (J00-J99)',
    },
}

print("="*70)
print("04d: CAUSE-SPECIFIC STRATIFICATION (REGIONAL DLNM) - v2")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Using: Per-region fits + meta-analysis")
print(f"CRITICAL: Zero-count days are KEPT (v1 bug fixed)")
print(f"Level: {level} [{level_label}]")

# =============================================================================
# 1. LOAD DATA
# =============================================================================
print("\n[1] Loading data...")

# Temperature data
temp_df = pd.read_parquet(PHASE0_RESULTS / era5_filename)
temp_df['date'] = pd.to_datetime(temp_df['date'])
print(f"  Temperature records: {len(temp_df):,}")

# Mortality data with cause stratification
mort_df = pd.read_parquet(PHASE0_RESULTS / mort_filename)
mort_df['date'] = pd.to_datetime(mort_df['date'])
print(f"  Mortality records: {len(mort_df):,}")

# Check available columns
available_cols = mort_df.columns.tolist()
print(f"  Available columns: {available_cols}")

# Standardize region column name
if level == 'immediate':
    mort_df = mort_df.rename(columns={'immediate_code': 'region_code'})
elif level == 'intermediate':
    if 'intermediate_code' in mort_df.columns:
        mort_df = mort_df.rename(columns={'intermediate_code': 'region_code'})

# Merge data
df = mort_df.merge(temp_df, on=['date', 'region_code'], how='inner')
print(f"  Merged dataset: {len(df):,} observations")

# Add time variables
df['dow'] = df['date'].dt.dayofweek
df['month'] = df['date'].dt.month
df['year'] = df['date'].dt.year

# Calculate "other" causes (all-cause minus CVD and respiratory)
if 'deaths_elderly_resp' in df.columns and 'deaths_elderly_cvd' in df.columns:
    df['deaths_elderly_other'] = (df['deaths_elderly'] - 
                                   df['deaths_elderly_resp'] - 
                                   df['deaths_elderly_cvd'])
    # NOTE: We allow negative values to become zero, but we DO NOT filter them out!
    df['deaths_elderly_other'] = df['deaths_elderly_other'].clip(lower=0)
    CAUSES['other'] = {
        'col': 'deaths_elderly_other',
        'label': 'Other causes',
    }
    print(f"  Created 'other causes' column")

# Check for heat-related deaths
if 'deaths_elderly_heat' in df.columns:
    CAUSES['heat_related'] = {
        'col': 'deaths_elderly_heat',
        'label': 'Heat-related (X30)',
    }
    print(f"  Heat-related deaths column found")

# Ensure population column exists
if 'pop_elderly' not in df.columns:
    try:
        covariates = pd.read_csv(PHASE0_RESULTS / pop_covariates_filename)
        pop_map = covariates.set_index(pop_region_col)['pop_elderly'].to_dict()
        df['pop_elderly'] = df['region_code'].map(pop_map)
        df['pop_elderly'] = df['pop_elderly'].fillna(10000)
        print(f"  Added pop_elderly from covariates table")
    except Exception:
        df['pop_elderly'] = 10000
        print(f"  WARNING: Using default population")

# Print death counts by cause
print("\n  Death counts by cause:")
for cause_key, cause_info in CAUSES.items():
    col = cause_info['col']
    if col in df.columns:
        total = df[col].sum()
        mean_daily = df[col].mean()
        zero_days = (df[col] == 0).sum()
        pct_zero = 100 * zero_days / len(df)
        print(f"    {cause_info['label']}: {total:,.0f} total, {mean_daily:.2f} daily mean, "
              f"{pct_zero:.1f}% zero days (KEPT)")

# =============================================================================
# 2. COMPUTE REGION-SPECIFIC PERCENTILES
# =============================================================================
print("\n[2] Computing region-specific temperature percentiles...")

regions = df['region_code'].unique()
region_percentiles = {}

for region in regions:
    temps = df.loc[df['region_code'] == region, 'temp_mean']
    if len(temps) < 365:
        continue
    region_percentiles[region] = {
        'p1': float(temps.quantile(0.01)),
        'p2_5': float(temps.quantile(0.025)),
        'p10': float(temps.quantile(0.10)),
        'p50': float(temps.quantile(0.50)),
        'p90': float(temps.quantile(0.90)),
        'p97_5': float(temps.quantile(0.975)),
        'p99': float(temps.quantile(0.99)),
    }

print(f"  Computed percentiles for {len(region_percentiles)} regions")

# =============================================================================
# 3. FIT PER-REGION DLNM FOR EACH CAUSE
# =============================================================================
print("\n[3] Fitting per-region DLNMs for each cause...")

def fit_region_cause(df_region: pd.DataFrame, deaths_col: str, pctiles: Dict) -> Optional[Dict]:
    """
    Fit DLNM for a region and cause, extract effects relative to MMT.
    
    Following Gasparrini methodology:
    - Find MMT (minimum mortality temperature) from the fitted curve
    - Compare extreme percentiles (P99, P1) to MMT
    - Store coefficients for MVMeta pooling
    
    IMPORTANT: We keep all rows including zero deaths!
    Poisson GLM can handle zeros - they are valid observations.
    """
    try:
        # Check if we have enough non-zero days for stable estimation
        # But we DO NOT remove zero days - just check for feasibility
        n_nonzero = (df_region[deaths_col] > 0).sum()
        if n_nonzero < 100:
            return None
        
        fit_res = fit_region_dlnm(
            df_region=df_region,
            temp_col='temp_mean',
            deaths_col=deaths_col,
            pop_col='pop_elderly',
            max_lag=MAX_LAG,
            var_df=VAR_DF,
            lag_df=LAG_DF,
            family='quasi-poisson',
            min_obs=MIN_OBS,
        )
        
        if fit_res is None:
            return None
        
        # Find MMT following Gasparrini methodology
        # Search over temperature range from P1 to P99
        temp_range = np.linspace(pctiles['p1'], pctiles['p99'], 100)
        mmt = find_mmt(fit_res, temp_range, pctiles['p50'])
        
        # Calculate MMT percentile (position in temperature distribution)
        # Approximate based on position in the P1-P99 range
        mmt_percentile = 1 + 98 * (mmt - pctiles['p1']) / (pctiles['p99'] - pctiles['p1'])
        mmt_percentile = max(1, min(99, mmt_percentile))  # Clamp to 1-99
        
        # Use MMT as reference (Gasparrini standard)
        ref_temp = mmt
        
        # Compare P99 (heat) and P1 (cold) to MMT
        heat = predict_cumulative_rr(fit_res, pctiles['p99'], ref_temp)
        cold = predict_cumulative_rr(fit_res, pctiles['p1'], ref_temp)
        
        # Extract cross-basis coefficients for MVMeta pooling
        cb_colnames = fit_res.get('cb_colnames', [])
        params = fit_res.get('params', pd.Series())
        cov = fit_res.get('cov', pd.DataFrame())
        
        # Get CB coefficients and vcov
        if len(cb_colnames) > 0 and all(c in params.index for c in cb_colnames):
            cb_coefs = params[cb_colnames].values
            cb_vcov = cov.loc[cb_colnames, cb_colnames].values
        else:
            cb_coefs = np.array([])
            cb_vcov = np.array([[]])
        
        return {
            'n_obs': fit_res['n_obs'],
            'n_nonzero': int(n_nonzero),
            'mmt': float(mmt),
            'mmt_percentile': float(mmt_percentile),
            'heat_rr': heat[0],
            'heat_rr_lo': heat[1],
            'heat_rr_hi': heat[2],
            'heat_log_rr': heat[3],
            'heat_log_rr_se': heat[4] if len(heat) > 4 else np.nan,
            'cold_rr': cold[0],
            'cold_rr_lo': cold[1],
            'cold_rr_hi': cold[2],
            'cold_log_rr': cold[3],
            'cold_log_rr_se': cold[4] if len(cold) > 4 else np.nan,
            'dispersion': fit_res.get('dispersion', np.nan),
            # For MVMeta pooling
            'cb_coefs': cb_coefs.tolist() if len(cb_coefs) > 0 else None,
            'cb_vcov': cb_vcov.tolist() if cb_vcov.size > 0 else None,
            'cb_meta': fit_res.get('cb_meta'),
            'temp_percentiles': pctiles,
        }
    except Exception as e:
        return None

cause_results = {}

for cause_key, cause_info in CAUSES.items():
    col = cause_info['col']
    label = cause_info['label']
    
    if col not in df.columns:
        print(f"\n  Skipping {label} - column not found")
        continue
    
    print(f"\n  === {label} ===")
    
    region_results = {}
    successful = 0
    skipped_sparse = 0
    
    for i, region in enumerate(regions):
        if (i + 1) % 30 == 0:
            print(f"    Region {i+1}/{len(regions)}...")
        
        if region not in region_percentiles:
            continue
        
        df_region = df[df['region_code'] == region].copy()
        result = fit_region_cause(df_region, col, region_percentiles[region])
        
        if result is not None:
            region_results[int(region)] = result
            successful += 1
        else:
            skipped_sparse += 1
    
    cause_results[cause_key] = region_results
    print(f"    Fitted: {successful}/{len(regions)} regions")
    if skipped_sparse > 0:
        print(f"    Skipped (sparse data): {skipped_sparse}")

# =============================================================================
# 4. META-ANALYZE WITHIN EACH CAUSE USING MVMETA
# =============================================================================
print("\n[4] Meta-analyzing results within each cause using MVMeta...")

pooled_results = {}

for cause_key, cause_info in CAUSES.items():
    label = cause_info['label']
    region_results = cause_results.get(cause_key, {})
    
    if len(region_results) < 5:
        print(f"\n  {label}: Insufficient regions ({len(region_results)})")
        continue
    
    print(f"\n  === {label} ===")
    
    # Collect coefficients for MVMeta pooling
    mvmeta_coefs = []
    mvmeta_vcovs = []
    mvmeta_pctiles = []
    valid_regions = []
    
    for region_code, res in region_results.items():
        if res.get('cb_coefs') is not None and res.get('cb_vcov') is not None:
            coefs = np.array(res['cb_coefs'])
            vcov = np.array(res['cb_vcov'])
            if coefs.shape[0] > 0 and vcov.shape[0] == coefs.shape[0]:
                mvmeta_coefs.append(coefs)
                mvmeta_vcovs.append(vcov)
                mvmeta_pctiles.append(res.get('temp_percentiles', {}))
                valid_regions.append(region_code)
    
    if len(mvmeta_coefs) < 5:
        print(f"    Insufficient valid coefficients for MVMeta: {len(mvmeta_coefs)}")
        continue
    
    print(f"    Pooling {len(mvmeta_coefs)} regions via MVMeta...")
    
    try:
        mvmeta_result = mvmeta_pool_coefficients(mvmeta_coefs, mvmeta_vcovs)
        
        if mvmeta_result is None:
            print(f"    MVMeta pooling failed, falling back to simple meta-analysis")
            # Fallback to simple meta-analysis
            heat_effects = []
            heat_variances = []
            cold_effects = []
            cold_variances = []
            for res in region_results.values():
                if res['heat_rr'] > 0 and res['heat_rr_lo'] > 0 and res['heat_rr_hi'] > 0:
                    log_rr = np.log(res['heat_rr'])
                    se = (np.log(res['heat_rr_hi']) - np.log(res['heat_rr_lo'])) / (2 * 1.96)
                    if se > 0.01 and se < 2.0:
                        heat_effects.append(log_rr)
                        heat_variances.append(se ** 2)
                if res['cold_rr'] > 0 and res['cold_rr_lo'] > 0 and res['cold_rr_hi'] > 0:
                    log_rr = np.log(res['cold_rr'])
                    se = (np.log(res['cold_rr_hi']) - np.log(res['cold_rr_lo'])) / (2 * 1.96)
                    if se > 0.01 and se < 2.0:
                        cold_effects.append(log_rr)
                        cold_variances.append(se ** 2)
            
            heat_meta = meta_random_effects(np.array(heat_effects), np.array(heat_variances))
            cold_meta = meta_random_effects(np.array(cold_effects), np.array(cold_variances))
            
            pooled_results[cause_key] = {
                'label': label,
                'method': 'simple_meta',
                'heat': {
                    'rr': float(np.exp(heat_meta['effect'])),
                    'rr_lo': float(np.exp(heat_meta['effect'] - 1.96 * heat_meta['se'])),
                    'rr_hi': float(np.exp(heat_meta['effect'] + 1.96 * heat_meta['se'])),
                    'I2': float(heat_meta['I2']),
                    'k': int(heat_meta['k']),
                },
                'cold': {
                    'rr': float(np.exp(cold_meta['effect'])),
                    'rr_lo': float(np.exp(cold_meta['effect'] - 1.96 * cold_meta['se'])),
                    'rr_hi': float(np.exp(cold_meta['effect'] + 1.96 * cold_meta['se'])),
                    'I2': float(cold_meta['I2']),
                    'k': int(cold_meta['k']),
                },
                'n_regions': len(region_results),
            }
        else:
            # Compute pooled RRs from MVMeta result
            all_p1 = [p.get('p1', 10) for p in mvmeta_pctiles if p]
            all_p50 = [p.get('p50', 20) for p in mvmeta_pctiles if p]
            all_p99 = [p.get('p99', 35) for p in mvmeta_pctiles if p]
            
            pooled_pctiles = {
                'p1': float(np.median(all_p1)) if all_p1 else 10.0,
                'p50': float(np.median(all_p50)) if all_p50 else 20.0,
                'p99': float(np.median(all_p99)) if all_p99 else 35.0,
            }
            
            # Get cb_meta from first valid result
            cb_meta = None
            for res in region_results.values():
                if res.get('cb_meta'):
                    cb_meta = res['cb_meta']
                    break
            
            if cb_meta is None:
                print(f"    Missing cb_meta, using defaults")
                cb_meta = {'var_df': VAR_DF, 'lag_df': LAG_DF, 'max_lag': MAX_LAG}
            
            rr_results = compute_pooled_rr_from_mvmeta(
                mvmeta_result, pooled_pctiles, cb_meta
            )
            
            pooled_results[cause_key] = {
                'label': label,
                'method': 'mvmeta',
                'heat': {
                    'rr': rr_results['heat_rr'],
                    'rr_lo': rr_results['heat_rr_lo'],
                    'rr_hi': rr_results['heat_rr_hi'],
                    'k': len(mvmeta_coefs),
                },
                'cold': {
                    'rr': rr_results['cold_rr'],
                    'rr_lo': rr_results['cold_rr_lo'],
                    'rr_hi': rr_results['cold_rr_hi'],
                    'k': len(mvmeta_coefs),
                },
                'mmt': rr_results['mmt'],
                'mmt_percentile': rr_results['mmt_percentile'],
                'n_regions': len(region_results),
                'n_pooled': len(mvmeta_coefs),
                'pooled_percentiles': pooled_pctiles,
            }
            
    except Exception as e:
        print(f"    MVMeta error: {e}")
        continue
    
    print(f"    Heat RR: {pooled_results[cause_key]['heat']['rr']:.3f} "
          f"[{pooled_results[cause_key]['heat']['rr_lo']:.3f}-{pooled_results[cause_key]['heat']['rr_hi']:.3f}] "
          f"(k={pooled_results[cause_key]['heat']['k']})")
    print(f"    Cold RR: {pooled_results[cause_key]['cold']['rr']:.3f} "
          f"[{pooled_results[cause_key]['cold']['rr_lo']:.3f}-{pooled_results[cause_key]['cold']['rr_hi']:.3f}] "
          f"(k={pooled_results[cause_key]['cold']['k']})")

# =============================================================================
# 5. HETEROGENEITY ACROSS CAUSES
# =============================================================================
print("\n[5] Testing heterogeneity across causes...")

def test_heterogeneity_causes(results: Dict, effect_type: str) -> Dict:
    """Cochran's Q test for heterogeneity across causes.
    Works with MVMeta results that have RR and CIs."""
    log_rrs = []
    ses = []
    causes = []
    
    for cause_key, res in results.items():
        if res is None or effect_type not in res:
            continue
        eff = res[effect_type]
        rr = eff.get('rr', 1.0)
        rr_lo = eff.get('rr_lo', rr)
        rr_hi = eff.get('rr_hi', rr)
        
        if rr > 0 and rr_lo > 0 and rr_hi > 0:
            log_rr = np.log(rr)
            se = (np.log(rr_hi) - np.log(rr_lo)) / (2 * 1.96)
            if se > 0.001:  # Valid SE
                log_rrs.append(log_rr)
                ses.append(se)
                causes.append(cause_key)
    
    if len(log_rrs) < 2:
        return {'Q': np.nan, 'df': 0, 'p': np.nan, 'I2': np.nan, 'causes': causes}
    
    log_rrs = np.array(log_rrs)
    variances = np.array(ses) ** 2
    weights = 1.0 / variances
    
    pooled = np.sum(weights * log_rrs) / np.sum(weights)
    Q = np.sum(weights * (log_rrs - pooled) ** 2)
    df = len(log_rrs) - 1
    p = 1 - stats.chi2.cdf(Q, df)
    I2 = max(0, (Q - df) / Q * 100) if Q > 0 else 0
    
    return {'Q': float(Q), 'df': int(df), 'p': float(p), 'I2': float(I2), 'causes': causes}

het_heat = test_heterogeneity_causes(pooled_results, 'heat')
het_cold = test_heterogeneity_causes(pooled_results, 'cold')

print(f"\n  Heat effect heterogeneity across causes:")
print(f"    Q = {het_heat['Q']:.2f}, df = {het_heat['df']}, p = {het_heat['p']:.3f}")
print(f"    I² = {het_heat['I2']:.1f}%")

print(f"\n  Cold effect heterogeneity across causes:")
print(f"    Q = {het_cold['Q']:.2f}, df = {het_cold['df']}, p = {het_cold['p']:.3f}")
print(f"    I² = {het_cold['I2']:.1f}%")

# =============================================================================
# 6. PAIRWISE COMPARISONS
# =============================================================================
print("\n[6] Pairwise comparisons between causes...")

comparisons = {}

# CVD vs Respiratory
if 'cardiovascular' in pooled_results and 'respiratory' in pooled_results:
    comparisons['cvd_vs_respiratory'] = {}
    
    for effect_type in ['heat', 'cold']:
        cvd = pooled_results['cardiovascular'][effect_type]
        resp = pooled_results['respiratory'][effect_type]
        
        diff = cvd['log_rr'] - resp['log_rr']
        se_diff = np.sqrt(cvd['se']**2 + resp['se']**2)
        z = diff / se_diff if se_diff > 0 else 0
        p = 2 * (1 - stats.norm.cdf(abs(z)))
        
        rr_ratio = np.exp(diff)
        rr_ratio_lo = np.exp(diff - 1.96 * se_diff)
        rr_ratio_hi = np.exp(diff + 1.96 * se_diff)
        
        comparisons['cvd_vs_respiratory'][effect_type] = {
            'rr_ratio': float(rr_ratio),
            'rr_ratio_lo': float(rr_ratio_lo),
            'rr_ratio_hi': float(rr_ratio_hi),
            'z': float(z),
            'p': float(p),
            'significant': bool(p < 0.05),
        }
        
        sig = "*" if p < 0.05 else ""
        print(f"\n  {effect_type.upper()}: CVD vs Respiratory")
        print(f"    CVD RR: {cvd['rr']:.3f} [{cvd['rr_lo']:.3f}-{cvd['rr_hi']:.3f}]")
        print(f"    Respiratory RR: {resp['rr']:.3f} [{resp['rr_lo']:.3f}-{resp['rr_hi']:.3f}]")
        print(f"    RR ratio (CVD/Resp): {rr_ratio:.3f} [{rr_ratio_lo:.3f}-{rr_ratio_hi:.3f}], p = {p:.3f}{sig}")

# CVD vs All-cause (to see if CVD-specific is stronger)
if 'cardiovascular' in pooled_results and 'all_cause' in pooled_results:
    comparisons['cvd_vs_all_cause'] = {}
    
    for effect_type in ['heat', 'cold']:
        cvd = pooled_results['cardiovascular'][effect_type]
        all_cause = pooled_results['all_cause'][effect_type]
        
        diff = cvd['log_rr'] - all_cause['log_rr']
        se_diff = np.sqrt(cvd['se']**2 + all_cause['se']**2)
        z = diff / se_diff if se_diff > 0 else 0
        p = 2 * (1 - stats.norm.cdf(abs(z)))
        
        rr_ratio = np.exp(diff)
        rr_ratio_lo = np.exp(diff - 1.96 * se_diff)
        rr_ratio_hi = np.exp(diff + 1.96 * se_diff)
        
        comparisons['cvd_vs_all_cause'][effect_type] = {
            'rr_ratio': float(rr_ratio),
            'rr_ratio_lo': float(rr_ratio_lo),
            'rr_ratio_hi': float(rr_ratio_hi),
            'z': float(z),
            'p': float(p),
            'significant': bool(p < 0.05),
        }

# =============================================================================
# 7. SAVE RESULTS
# =============================================================================
print("\n[7] Saving results...")

summary = {
    'analysis_date': datetime.now().isoformat(),
    'level': level,
    'causes': {k: v['label'] for k, v in CAUSES.items()},
    'dlnm_params': {
        'max_lag': MAX_LAG,
        'var_df': VAR_DF,
        'lag_df': LAG_DF,
        'basis_type': 'natural_spline',
    },
    'approach': 'per_region_dlnm_then_mvmeta_pooling',
    'methodology': 'Gasparrini MVMeta - pool DLNM coefficients not RRs',
    'offset': 'log(pop_elderly)',
    'percentile_method': 'region_specific',
    'zero_count_handling': 'KEPT (v1 bug fixed)',
    'pooled_results': pooled_results,
    'heterogeneity_across_causes': {
        'heat': het_heat,
        'cold': het_cold,
    },
    'pairwise_comparisons': comparisons,
}

# Region-level results by cause
region_output = {
    cause: {int(k): convert_to_json_serializable(v) for k, v in results.items()}
    for cause, results in cause_results.items()
}

# Use level suffix for output files
suffix = '' if level == 'intermediate' else '_immediate'

with open(OUTPUT_DIR / f'cause_stratification_v2_results{suffix}.json', 'w') as f:
    json.dump(convert_to_json_serializable(summary), f, indent=2)

with open(OUTPUT_DIR / f'cause_stratification_v2_region_effects{suffix}.json', 'w') as f:
    json.dump(region_output, f, indent=2)

print(f"\n  Saved to {OUTPUT_DIR}:")
print(f"    - cause_stratification_v2_results{suffix}.json")
print(f"    - cause_stratification_v2_region_effects{suffix}.json")

print(f"\n{'='*70}")
print("04d: CAUSE-SPECIFIC STRATIFICATION COMPLETE")
print(f"{'='*70}")
print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
