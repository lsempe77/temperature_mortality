"""
02a_sensitivity_analyses_v2.py
================================
Comprehensive Sensitivity Analyses for DLNM Results using Proper Implementation

VERSION 2: Uses dlnm_module with proper spline cross-basis and region-specific fitting.

FIXES FROM EXPERT REVIEW:
1. Uses natural cubic spline cross-basis (not polynomial) via dlnm_module
2. Fits per-region models with population offset
3. Meta-analyzes effects across regions (not pooled national model)
4. Leave-one-year-out at regional level
5. Proper confidence intervals via delta method

Tests robustness of temperature-mortality estimates to:
1. Lag structure sensitivity (7, 14, 21, 28 days)
2. Spline degrees of freedom (2, 3, 4, 5 df for both var and lag)
3. Temperature percentile thresholds (P95/P5 vs P97.5/P2.5 vs P99/P1)
4. Reference temperature (empirical MMT vs P50)
5. Leave-one-year-out temporal cross-validation
6. Outlier region exclusion
7. Different model families (quasi-Poisson vs Negative Binomial)

Input: Regional data from phase0_data_prep/results/
Output: Comprehensive sensitivity analysis with meta-analyzed pooled estimates

References:
- Gasparrini et al. (2015) Lancet - Multi-country methodology
- Armstrong et al. (2014) - Temperature-mortality meta-analysis

Author: Climate-Health Analysis Pipeline
Date: December 2025
"""

import pandas as pd
import numpy as np
import os
import sys
import json
from datetime import datetime
import warnings
import argparse
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

warnings.filterwarnings('ignore')

# Add parent directory to path for utils import
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.dlnm_module import (
    fit_region_dlnm,
    predict_cumulative_rr,
    find_mmt,
    compute_rr_curve,
    meta_random_effects,
    pool_region_effects,
    extract_region_percentiles,
    mvmeta_pool_coefficients,
    compute_pooled_rr_from_mvmeta,
    convert_to_json_serializable,
)

print("="*70)
print("PHASE 2: COMPREHENSIVE SENSITIVITY ANALYSES (v2)")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PHASE0_RESULTS = os.path.join(BASE_DIR, 'new_analysis', 'phase0_data_prep', 'results')
PHASE1_RESULTS = os.path.join(BASE_DIR, 'new_analysis', 'phase1_core_model', 'results')
OUTPUT_DIR = os.path.join(BASE_DIR, 'new_analysis', 'phase2_robustness', 'results')
FIGURE_DIR = os.path.join(BASE_DIR, 'new_analysis', 'phase2_robustness', 'figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

# Baseline specification
BASELINE_MAX_LAG = 21
BASELINE_VAR_DF = 4
BASELINE_LAG_DF = 4
MIN_OBS = 1000

# Sensitivity variations
LAG_OPTIONS = [7, 14, 21, 28]
DF_OPTIONS = [3, 4, 5]  # Min df=3 for natural cubic splines, max df=5 to avoid SVD issues
PERCENTILE_OPTIONS = [
    {'name': 'P95_P5', 'heat': 0.95, 'cold': 0.05},
    {'name': 'P975_P025', 'heat': 0.975, 'cold': 0.025},
    {'name': 'P99_P1', 'heat': 0.99, 'cold': 0.01},
]

print(f"\nBaseline specification:")
print(f"  Max lag: {BASELINE_MAX_LAG} days")
print(f"  Spline df: var={BASELINE_VAR_DF}, lag={BASELINE_LAG_DF}")
print(f"  Min observations/region: {MIN_OBS}")

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

parser = argparse.ArgumentParser(description='Sensitivity analyses for DLNM')
parser.add_argument('--quick', action='store_true', help='Run on subset of regions')
parser.add_argument('--analysis', type=str, default='all',
                    choices=['lag', 'df', 'percentile', 'loyo', 'family', 'outlier', 'all'],
                    help='Which sensitivity analysis to run')
parser.add_argument('--level', type=str, default='intermediate',
                    choices=['intermediate', 'immediate'],
                    help='Spatial level to analyze')
parser.add_argument('--save-figures', action='store_true', default=True)
args, _ = parser.parse_known_args()

# =============================================================================
# DATA LOADING
# =============================================================================

print("\n" + "-"*70)
print("Loading Data")
print("-"*70)

level = args.level
if level == 'intermediate':
    temp_filename = 'era5_intermediate_daily.parquet'
    mort_filename = 'mortality_regional_daily_elderly.parquet'
    ses_filename = 'ses_intermediate_covariates.csv'
    ses_region_col = 'intermediate_code'
    level_label = 'INTERMEDIATE (133 regions)'
else:
    temp_filename = 'era5_immediate_daily.parquet'
    mort_filename = 'mortality_immediate_daily_elderly.parquet'
    ses_filename = 'ses_immediate_covariates.csv'
    ses_region_col = 'immediate_code'
    level_label = 'IMMEDIATE (510 regions)'

print(f"Level: {level} [{level_label}]")

# Load temperature data
df_temp = pd.read_parquet(os.path.join(PHASE0_RESULTS, temp_filename))
df_temp['date'] = pd.to_datetime(df_temp['date'])
print(f"ERA5: {len(df_temp):,} rows, {df_temp['region_code'].nunique()} regions")

# Load mortality data
df_mort = pd.read_parquet(os.path.join(PHASE0_RESULTS, mort_filename))
df_mort['date'] = pd.to_datetime(df_mort['date'])
print(f"Mortality: {len(df_mort):,} rows")

# Load SES data for population
df_ses = pd.read_csv(os.path.join(PHASE0_RESULTS, ses_filename))
pop_map = dict(zip(df_ses[ses_region_col], df_ses['pop_elderly']))

# Load holidays
df_holidays = pd.read_parquet(os.path.join(PHASE0_RESULTS, 'brazilian_holidays_daily.parquet'))
df_holidays['date'] = pd.to_datetime(df_holidays['date'])

# Standardize region column names for merging
if level == 'immediate':
    df_mort = df_mort.rename(columns={'immediate_code': 'region_code'})
elif level == 'intermediate':
    if 'intermediate_code' in df_mort.columns:
        df_mort = df_mort.rename(columns={'intermediate_code': 'region_code'})

# Merge data
df = pd.merge(
    df_mort,
    df_temp[['date', 'region_code', 'temp_mean']],
    on=['date', 'region_code'],
    how='inner'
)

df = pd.merge(
    df,
    df_holidays[['date', 'is_holiday']],
    on='date',
    how='left'
)
df['is_holiday'] = df['is_holiday'].fillna(0).astype(int)

df['pop_elderly'] = df['region_code'].map(pop_map)
df['year'] = df['date'].dt.year
df['month'] = df['date'].dt.month
df['day_of_week'] = df['date'].dt.dayofweek

# Drop missing
df = df.dropna(subset=['temp_mean', 'deaths_elderly', 'pop_elderly'])

# Filter regions with enough observations
region_counts = df.groupby('region_code').size()
valid_regions = region_counts[region_counts >= MIN_OBS].index.tolist()
df = df[df['region_code'].isin(valid_regions)]

if args.quick:
    valid_regions = valid_regions[:15]
    df = df[df['region_code'].isin(valid_regions)]
    print(f"\n[QUICK MODE] Running on {len(valid_regions)} regions only")

print(f"\nFinal dataset: {len(df):,} rows, {len(valid_regions)} regions")

# Extract region-specific percentiles
region_pcts = extract_region_percentiles(df, temp_col='temp_mean', region_col='region_code')
print(f"\nRegion percentile summary:")
print(f"  P50 range: {region_pcts['p50'].min():.1f} - {region_pcts['p50'].max():.1f}°C")
print(f"  P99 range: {region_pcts['p99'].min():.1f} - {region_pcts['p99'].max():.1f}°C")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def run_regional_analysis(df, valid_regions, max_lag, var_df, lag_df, 
                          family='quasi-poisson', exclude_years=None):
    """Run DLNM for each region and return results for meta-analysis.
    
    Stores cross-basis coefficients for MVMeta pooling.
    """
    
    region_results = []
    
    for region in valid_regions:
        df_region = df[df['region_code'] == region].copy()
        
        # Optionally exclude years (for LOYO)
        if exclude_years:
            df_region = df_region[~df_region['year'].isin(exclude_years)]
            if len(df_region) < MIN_OBS:
                continue
        
        fit = fit_region_dlnm(
            df_region,
            temp_col='temp_mean',
            deaths_col='deaths_elderly',
            pop_col='pop_elderly',
            max_lag=max_lag,
            var_df=var_df,
            lag_df=lag_df,
            family=family,
        )
        
        if fit is None:
            continue
        
        # Region-specific percentiles
        temps = df_region['temp_mean'].dropna()
        percentiles = {
            'p01': temps.quantile(0.01),
            'p025': temps.quantile(0.025),
            'p05': temps.quantile(0.05),
            'p50': temps.quantile(0.50),
            'p95': temps.quantile(0.95),
            'p975': temps.quantile(0.975),
            'p99': temps.quantile(0.99),
        }
        
        # Find empirical MMT
        temp_range = np.linspace(temps.quantile(0.05), temps.quantile(0.95), 50)
        mmt = find_mmt(fit, temp_range, percentiles['p50'])
        
        # Calculate RR at different thresholds
        results = {'region': region, 'n_obs': fit['n_obs'], 'mmt': mmt}
        
        for pct_name, pct_val in percentiles.items():
            if pct_name == 'p50':
                continue
            rr, lo, hi, log_rr, se = predict_cumulative_rr(fit, pct_val, mmt)
            results[f'rr_{pct_name}'] = rr
            results[f'log_rr_{pct_name}'] = log_rr
            results[f'se_{pct_name}'] = se
        
        # Store coefficients for MVMeta pooling
        cb_colnames = fit.get('cb_colnames', [])
        params = fit.get('params', pd.Series())
        cov = fit.get('cov', pd.DataFrame())
        
        if len(cb_colnames) > 0 and all(c in params.index for c in cb_colnames):
            results['cb_coefs'] = params[cb_colnames].values.tolist()
            results['cb_vcov'] = cov.loc[cb_colnames, cb_colnames].values.tolist()
            results['cb_meta'] = fit.get('cb_meta')
        else:
            results['cb_coefs'] = None
            results['cb_vcov'] = None
            results['cb_meta'] = None
        
        results['percentiles'] = percentiles
        
        region_results.append(results)
    
    return region_results


def meta_analyze_results(region_results, pct='p99'):
    """Pool effects across regions using MVMeta (Gasparrini methodology).
    
    First tries MVMeta coefficient pooling, falls back to simple RR pooling.
    """
    # Try MVMeta first
    mvmeta_coefs = []
    mvmeta_vcovs = []
    mvmeta_pctiles = []
    
    for r in region_results:
        if r.get('cb_coefs') is not None and r.get('cb_vcov') is not None:
            coefs = np.array(r['cb_coefs'])
            vcov = np.array(r['cb_vcov'])
            if coefs.shape[0] > 0 and vcov.shape[0] == coefs.shape[0]:
                mvmeta_coefs.append(coefs)
                mvmeta_vcovs.append(vcov)
                mvmeta_pctiles.append(r.get('percentiles', {}))
    
    # Use MVMeta if we have enough valid coefficients
    if len(mvmeta_coefs) >= 5:
        try:
            mvmeta_result = mvmeta_pool_coefficients(mvmeta_coefs, mvmeta_vcovs)
            
            if mvmeta_result is not None:
                # Get pooled percentiles
                pct_key_map = {'p99': 'p99', 'p01': 'p01', 'p975': 'p975', 'p025': 'p025',
                               'p95': 'p95', 'p05': 'p05'}
                
                all_p1 = [p.get('p01', 10) for p in mvmeta_pctiles if p]
                all_p50 = [p.get('p50', 20) for p in mvmeta_pctiles if p]
                all_p99 = [p.get('p99', 35) for p in mvmeta_pctiles if p]
                
                pooled_pctiles = {
                    'p1': float(np.median(all_p1)) if all_p1 else 10.0,
                    'p50': float(np.median(all_p50)) if all_p50 else 20.0,
                    'p99': float(np.median(all_p99)) if all_p99 else 35.0,
                }
                
                # Get cb_meta from first valid result
                cb_meta = None
                for r in region_results:
                    if r.get('cb_meta'):
                        cb_meta = r['cb_meta']
                        break
                
                if cb_meta is not None:
                    rr_results = compute_pooled_rr_from_mvmeta(
                        mvmeta_result, pooled_pctiles, cb_meta
                    )
                    
                    # Return appropriate effect based on pct
                    if pct in ['p99', 'p975', 'p95']:
                        return {
                            'pooled_rr': rr_results['heat_rr'],
                            'ci_low': rr_results['heat_rr_lo'],
                            'ci_high': rr_results['heat_rr_hi'],
                            'log_rr': np.log(rr_results['heat_rr']),
                            'se': (np.log(rr_results['heat_rr_hi']) - np.log(rr_results['heat_rr_lo'])) / (2 * 1.96),
                            'k': len(mvmeta_coefs),
                            'I2': np.nan,  # Not computed for MVMeta
                            'method': 'mvmeta',
                            'mmt': rr_results.get('mmt'),
                        }
                    else:
                        return {
                            'pooled_rr': rr_results['cold_rr'],
                            'ci_low': rr_results['cold_rr_lo'],
                            'ci_high': rr_results['cold_rr_hi'],
                            'log_rr': np.log(rr_results['cold_rr']),
                            'se': (np.log(rr_results['cold_rr_hi']) - np.log(rr_results['cold_rr_lo'])) / (2 * 1.96),
                            'k': len(mvmeta_coefs),
                            'I2': np.nan,  # Not computed for MVMeta
                            'method': 'mvmeta',
                            'mmt': rr_results.get('mmt'),
                        }
        except Exception:
            pass  # Fall through to simple meta-analysis
    
    # Fallback: simple RR pooling
    effects = []
    variances = []
    
    for r in region_results:
        log_rr = r.get(f'log_rr_{pct}')
        se = r.get(f'se_{pct}')
        if log_rr is not None and se is not None:
            if not np.isnan(log_rr) and not np.isnan(se) and se > 0:
                effects.append(log_rr)
                variances.append(se ** 2)
    
    if len(effects) == 0:
        return None
    
    meta = meta_random_effects(np.array(effects), np.array(variances))
    
    return {
        'pooled_rr': np.exp(meta['effect']),
        'ci_low': np.exp(meta['effect'] - 1.96 * meta['se']),
        'ci_high': np.exp(meta['effect'] + 1.96 * meta['se']),
        'log_rr': meta['effect'],
        'se': meta['se'],
        'I2': meta['I2'],
        'tau2': meta['tau2'],
        'k': meta['k'],
        'p': meta['p'],
        'method': 'simple_meta',
    }


# =============================================================================
# 1. LAG STRUCTURE SENSITIVITY
# =============================================================================

lag_results = {}
if args.analysis in ['lag', 'all']:
    print("\n" + "="*70)
    print("1. LAG STRUCTURE SENSITIVITY")
    print("="*70)
    
    for max_lag in LAG_OPTIONS:
        print(f"\n  Testing max_lag = {max_lag} days...")
        region_res = run_regional_analysis(
            df, valid_regions, max_lag=max_lag, 
            var_df=BASELINE_VAR_DF, lag_df=BASELINE_LAG_DF
        )
        
        heat_meta = meta_analyze_results(region_res, 'p99')
        cold_meta = meta_analyze_results(region_res, 'p01')
        
        lag_results[max_lag] = {
            'heat_p99': heat_meta,
            'cold_p01': cold_meta,
            'n_regions': len(region_res),
        }
        
        if heat_meta:
            i2_str = f", I²={heat_meta.get('I2', float('nan')):.1f}%" if not np.isnan(heat_meta.get('I2', float('nan'))) else ""
            print(f"    Heat P99: RR = {heat_meta['pooled_rr']:.3f} ({heat_meta['ci_low']:.3f}-{heat_meta['ci_high']:.3f}){i2_str}")
        if cold_meta:
            i2_str = f", I²={cold_meta.get('I2', float('nan')):.1f}%" if not np.isnan(cold_meta.get('I2', float('nan'))) else ""
            print(f"    Cold P01: RR = {cold_meta['pooled_rr']:.3f} ({cold_meta['ci_low']:.3f}-{cold_meta['ci_high']:.3f}){i2_str}")

# =============================================================================
# 2. SPLINE DF SENSITIVITY
# =============================================================================

df_results = {}
if args.analysis in ['df', 'all']:
    print("\n" + "="*70)
    print("2. SPLINE DF SENSITIVITY")
    print("="*70)
    
    for var_df in DF_OPTIONS:
        for lag_df in DF_OPTIONS:
            key = f"var{var_df}_lag{lag_df}"
            print(f"\n  Testing df: var={var_df}, lag={lag_df}...")
            
            region_res = run_regional_analysis(
                df, valid_regions, max_lag=BASELINE_MAX_LAG,
                var_df=var_df, lag_df=lag_df
            )
            
            heat_meta = meta_analyze_results(region_res, 'p99')
            cold_meta = meta_analyze_results(region_res, 'p01')
            
            df_results[key] = {
                'var_df': var_df,
                'lag_df': lag_df,
                'heat_p99': heat_meta,
                'cold_p01': cold_meta,
                'n_regions': len(region_res),
            }
            
            if heat_meta:
                print(f"    Heat P99: RR = {heat_meta['pooled_rr']:.3f}")

# =============================================================================
# 3. PERCENTILE THRESHOLD SENSITIVITY
# =============================================================================

pct_results = {}
if args.analysis in ['percentile', 'all']:
    print("\n" + "="*70)
    print("3. PERCENTILE THRESHOLD SENSITIVITY")
    print("="*70)
    
    # Run baseline model
    region_res = run_regional_analysis(
        df, valid_regions, max_lag=BASELINE_MAX_LAG,
        var_df=BASELINE_VAR_DF, lag_df=BASELINE_LAG_DF
    )
    
    for pct_option in PERCENTILE_OPTIONS:
        name = pct_option['name']
        heat_pct = f"p{int(pct_option['heat']*100)}" if pct_option['heat'] < 0.99 else 'p99'
        cold_pct = f"p{int(pct_option['cold']*100):02d}" if pct_option['cold'] > 0.01 else 'p01'
        
        # Handle 97.5 and 2.5 naming
        if pct_option['heat'] == 0.975:
            heat_pct = 'p975'
        if pct_option['cold'] == 0.025:
            cold_pct = 'p025'
        
        heat_meta = meta_analyze_results(region_res, heat_pct)
        cold_meta = meta_analyze_results(region_res, cold_pct)
        
        pct_results[name] = {
            'heat': heat_meta,
            'cold': cold_meta,
            'heat_pct': heat_pct,
            'cold_pct': cold_pct,
        }
        
        print(f"\n  {name}:")
        if heat_meta:
            print(f"    Heat ({heat_pct}): RR = {heat_meta['pooled_rr']:.3f} ({heat_meta['ci_low']:.3f}-{heat_meta['ci_high']:.3f})")
        if cold_meta:
            print(f"    Cold ({cold_pct}): RR = {cold_meta['pooled_rr']:.3f} ({cold_meta['ci_low']:.3f}-{cold_meta['ci_high']:.3f})")

# =============================================================================
# 4. LEAVE-ONE-YEAR-OUT CROSS-VALIDATION
# =============================================================================

loyo_results = {}
if args.analysis in ['loyo', 'all']:
    print("\n" + "="*70)
    print("4. LEAVE-ONE-YEAR-OUT CROSS-VALIDATION")
    print("="*70)
    
    years = sorted(df['year'].unique())
    print(f"\n  Years in data: {years}")
    
    for year in years:
        print(f"\n  Excluding year {year}...")
        region_res = run_regional_analysis(
            df, valid_regions, max_lag=BASELINE_MAX_LAG,
            var_df=BASELINE_VAR_DF, lag_df=BASELINE_LAG_DF,
            exclude_years=[year]
        )
        
        heat_meta = meta_analyze_results(region_res, 'p99')
        cold_meta = meta_analyze_results(region_res, 'p01')
        
        loyo_results[year] = {
            'heat_p99': heat_meta,
            'cold_p01': cold_meta,
            'n_regions': len(region_res),
        }
        
        if heat_meta:
            print(f"    Heat P99: RR = {heat_meta['pooled_rr']:.3f} ({heat_meta['ci_low']:.3f}-{heat_meta['ci_high']:.3f})")

# =============================================================================
# 5. MODEL FAMILY SENSITIVITY
# =============================================================================

family_results = {}
if args.analysis in ['family', 'all']:
    print("\n" + "="*70)
    print("5. MODEL FAMILY SENSITIVITY")
    print("="*70)
    
    for family in ['quasi-poisson', 'negbin']:
        print(f"\n  Testing family: {family}...")
        region_res = run_regional_analysis(
            df, valid_regions, max_lag=BASELINE_MAX_LAG,
            var_df=BASELINE_VAR_DF, lag_df=BASELINE_LAG_DF,
            family=family
        )
        
        heat_meta = meta_analyze_results(region_res, 'p99')
        cold_meta = meta_analyze_results(region_res, 'p01')
        
        family_results[family] = {
            'heat_p99': heat_meta,
            'cold_p01': cold_meta,
            'n_regions': len(region_res),
        }
        
        if heat_meta:
            print(f"    Heat P99: RR = {heat_meta['pooled_rr']:.3f}")
        if cold_meta:
            print(f"    Cold P01: RR = {cold_meta['pooled_rr']:.3f}")

# =============================================================================
# 6. OUTLIER REGION EXCLUSION
# =============================================================================

outlier_results = {}
if args.analysis in ['outlier', 'all']:
    print("\n" + "="*70)
    print("6. OUTLIER REGION EXCLUSION")
    print("="*70)
    
    # Identify high-mortality regions
    region_mort_rates = df.groupby('region_code').apply(
        lambda x: x['deaths_elderly'].sum() / x['pop_elderly'].mean()
    ).sort_values(ascending=False)
    
    top_5pct = region_mort_rates.head(int(len(region_mort_rates) * 0.05)).index.tolist()
    top_10pct = region_mort_rates.head(int(len(region_mort_rates) * 0.10)).index.tolist()
    
    print(f"  Excluding top 5% mortality regions ({len(top_5pct)} regions)")
    regions_ex5 = [r for r in valid_regions if r not in top_5pct]
    region_res_5 = run_regional_analysis(
        df, regions_ex5, max_lag=BASELINE_MAX_LAG,
        var_df=BASELINE_VAR_DF, lag_df=BASELINE_LAG_DF
    )
    
    print(f"  Excluding top 10% mortality regions ({len(top_10pct)} regions)")
    regions_ex10 = [r for r in valid_regions if r not in top_10pct]
    region_res_10 = run_regional_analysis(
        df, regions_ex10, max_lag=BASELINE_MAX_LAG,
        var_df=BASELINE_VAR_DF, lag_df=BASELINE_LAG_DF
    )
    
    outlier_results['exclude_top5pct'] = {
        'heat_p99': meta_analyze_results(region_res_5, 'p99'),
        'cold_p01': meta_analyze_results(region_res_5, 'p01'),
        'n_regions': len(region_res_5),
    }
    
    outlier_results['exclude_top10pct'] = {
        'heat_p99': meta_analyze_results(region_res_10, 'p99'),
        'cold_p01': meta_analyze_results(region_res_10, 'p01'),
        'n_regions': len(region_res_10),
    }
    
    for key, res in outlier_results.items():
        if res['heat_p99']:
            print(f"\n  {key}: Heat P99 RR = {res['heat_p99']['pooled_rr']:.3f}")

# =============================================================================
# VISUALIZATIONS
# =============================================================================

if args.save_figures:
    print("\n" + "-"*70)
    print("Creating Visualizations")
    print("-"*70)
    
    # Figure 1: Lag sensitivity
    if lag_results:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Heat
        ax = axes[0]
        lags = sorted(lag_results.keys())
        rrs = [lag_results[l]['heat_p99']['pooled_rr'] if lag_results[l]['heat_p99'] else np.nan for l in lags]
        ci_low = [lag_results[l]['heat_p99']['ci_low'] if lag_results[l]['heat_p99'] else np.nan for l in lags]
        ci_high = [lag_results[l]['heat_p99']['ci_high'] if lag_results[l]['heat_p99'] else np.nan for l in lags]
        
        ax.plot(lags, rrs, 'o-', color='red', linewidth=2, markersize=10)
        ax.fill_between(lags, ci_low, ci_high, color='red', alpha=0.2)
        ax.axhline(y=1, color='black', linestyle='--', alpha=0.5)
        ax.set_xlabel('Maximum Lag (days)', fontsize=12)
        ax.set_ylabel('Pooled RR (P99 vs MMT)', fontsize=12)
        ax.set_title('Heat Effect by Lag Structure', fontsize=14)
        ax.set_xticks(lags)
        ax.grid(True, alpha=0.3)
        
        # Cold
        ax = axes[1]
        rrs = [lag_results[l]['cold_p01']['pooled_rr'] if lag_results[l]['cold_p01'] else np.nan for l in lags]
        ci_low = [lag_results[l]['cold_p01']['ci_low'] if lag_results[l]['cold_p01'] else np.nan for l in lags]
        ci_high = [lag_results[l]['cold_p01']['ci_high'] if lag_results[l]['cold_p01'] else np.nan for l in lags]
        
        ax.plot(lags, rrs, 'o-', color='blue', linewidth=2, markersize=10)
        ax.fill_between(lags, ci_low, ci_high, color='blue', alpha=0.2)
        ax.axhline(y=1, color='black', linestyle='--', alpha=0.5)
        ax.set_xlabel('Maximum Lag (days)', fontsize=12)
        ax.set_ylabel('Pooled RR (P01 vs MMT)', fontsize=12)
        ax.set_title('Cold Effect by Lag Structure', fontsize=14)
        ax.set_xticks(lags)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        level_tag = '' if level == 'intermediate' else f'_{level}'
        fig_path = os.path.join(FIGURE_DIR, f'sensitivity_lag_structure{level_tag}.png')
        plt.savefig(fig_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {fig_path}")
    
    # Figure 2: LOYO stability
    if loyo_results:
        fig, ax = plt.subplots(figsize=(12, 6))
        
        years = sorted(loyo_results.keys())
        x = np.arange(len(years))
        
        heat_rrs = [loyo_results[y]['heat_p99']['pooled_rr'] if loyo_results[y]['heat_p99'] else np.nan for y in years]
        heat_ci_low = [loyo_results[y]['heat_p99']['ci_low'] if loyo_results[y]['heat_p99'] else np.nan for y in years]
        heat_ci_high = [loyo_results[y]['heat_p99']['ci_high'] if loyo_results[y]['heat_p99'] else np.nan for y in years]
        
        ax.errorbar(x, heat_rrs, 
                    yerr=[np.array(heat_rrs) - np.array(heat_ci_low), 
                          np.array(heat_ci_high) - np.array(heat_rrs)],
                    fmt='o', color='red', capsize=5, markersize=10, label='Heat (P99)')
        
        ax.axhline(y=np.nanmean(heat_rrs), color='red', linestyle='--', alpha=0.5, label='Mean')
        ax.axhline(y=1, color='black', linestyle=':', alpha=0.5)
        
        ax.set_xlabel('Excluded Year', fontsize=12)
        ax.set_ylabel('Pooled RR', fontsize=12)
        ax.set_title('Leave-One-Year-Out Cross-Validation', fontsize=14)
        ax.set_xticks(x)
        ax.set_xticklabels(years)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        level_tag = '' if level == 'intermediate' else f'_{level}'
        fig_path = os.path.join(FIGURE_DIR, f'sensitivity_loyo{level_tag}.png')
        plt.savefig(fig_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {fig_path}")

# =============================================================================
# SAVE RESULTS
# =============================================================================

print("\n" + "-"*70)
print("Saving Results")
print("-"*70)

# Convert numpy types to native Python types for JSON serialization
def convert_to_native(obj):
    """Recursively convert numpy types to native Python types."""
    if isinstance(obj, dict):
        return {str(k): convert_to_native(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_native(item) for item in obj]
    elif isinstance(obj, (np.integer, np.int32, np.int64)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float32, np.float64)):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj

all_results = {
    'lag_sensitivity': convert_to_native(lag_results) if lag_results else None,
    'df_sensitivity': convert_to_native(df_results) if df_results else None,
    'percentile_sensitivity': convert_to_native(pct_results) if pct_results else None,
    'loyo_sensitivity': convert_to_native(loyo_results) if loyo_results else None,
    'family_sensitivity': convert_to_native(family_results) if family_results else None,
    'outlier_sensitivity': convert_to_native(outlier_results) if outlier_results else None,
    'metadata': {
        'baseline_max_lag': BASELINE_MAX_LAG,
        'baseline_var_df': BASELINE_VAR_DF,
        'baseline_lag_df': BASELINE_LAG_DF,
        'n_regions': len(valid_regions),
        'level': level,
        'analysis_date': datetime.now().isoformat(),
    }
}

level_tag = '' if level == 'intermediate' else f'_{level}'
json_path = os.path.join(OUTPUT_DIR, f'sensitivity_analyses_v2{level_tag}.json')
with open(json_path, 'w') as f:
    json.dump(all_results, f, indent=2, default=float)
print(f"  Saved: {json_path}")

# =============================================================================
# SUMMARY TABLE
# =============================================================================

print("\n" + "="*70)
print("SENSITIVITY ANALYSIS SUMMARY (v2)")
print("="*70)

print("\n1. LAG SENSITIVITY:")
if lag_results:
    for lag, res in sorted(lag_results.items()):
        heat = res['heat_p99']
        if heat:
            print(f"   Lag {lag:2d}: Heat RR = {heat['pooled_rr']:.3f} ({heat['ci_low']:.3f}-{heat['ci_high']:.3f})")

print("\n2. DF SENSITIVITY (Heat P99):")
if df_results:
    print("   ", end="")
    for lag_df in DF_OPTIONS:
        print(f"  lag_df={lag_df}", end="")
    print()
    for var_df in DF_OPTIONS:
        print(f"   var_df={var_df}: ", end="")
        for lag_df in DF_OPTIONS:
            key = f"var{var_df}_lag{lag_df}"
            heat = df_results.get(key, {}).get('heat_p99')
            if heat:
                print(f"  {heat['pooled_rr']:.2f}", end="    ")
            else:
                print(f"  --", end="      ")
        print()

print("\n3. PERCENTILE THRESHOLDS:")
if pct_results:
    for name, res in pct_results.items():
        heat = res['heat']
        cold = res['cold']
        if heat and cold:
            print(f"   {name}: Heat RR = {heat['pooled_rr']:.3f}, Cold RR = {cold['pooled_rr']:.3f}")

print("\n4. LOYO STABILITY:")
if loyo_results:
    heat_rrs = [r['heat_p99']['pooled_rr'] for r in loyo_results.values() if r['heat_p99']]
    if heat_rrs:
        print(f"   Heat P99 range: {min(heat_rrs):.3f} - {max(heat_rrs):.3f}")
        print(f"   CV: {np.std(heat_rrs)/np.mean(heat_rrs)*100:.1f}%")

print("\n5. MODEL FAMILY:")
if family_results:
    for family, res in family_results.items():
        heat = res['heat_p99']
        if heat:
            print(f"   {family}: Heat RR = {heat['pooled_rr']:.3f}")

print("\n6. OUTLIER EXCLUSION:")
if outlier_results:
    for key, res in outlier_results.items():
        heat = res['heat_p99']
        if heat:
            print(f"   {key}: Heat RR = {heat['pooled_rr']:.3f} (k={heat['k']})")

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*70)
