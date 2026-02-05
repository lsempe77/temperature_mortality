"""
02b_harvesting_analysis_v2.py
==============================
Harvesting (Mortality Displacement) Analysis using Proper DLNM

VERSION 2: Uses dlnm_module with proper spline cross-basis and region-specific fitting.

FIXES FROM EXPERT REVIEW:
1. Uses natural cubic spline cross-basis (not polynomial) via dlnm_module
2. Fits per-region models with population offset
3. Meta-analyzes cumulative RR at multiple horizons
4. Proper delta-method confidence intervals for cumulative RR
5. Harvesting ratio computed from Excess Relative Risk (ERR = RR - 1)

Tests whether temperature-attributable deaths represent:
- TRUE EXCESS: Deaths that would not have occurred otherwise (ERR stable across lags)
- HARVESTING: Short-term displacement (ERR decreases at longer lags)

Method (Armstrong et al., 2014; Gasparrini et al., 2015):
1. Fit DLNM with extended lag (up to 35 days) for each region
2. Calculate CUMULATIVE RR at different lag horizons (7, 14, 21, 28, 35)
3. Pool via DerSimonian-Laird meta-analysis
4. Harvesting ratio = 1 - (ERR_35 / ERR_7)
   - Positive ratio → harvesting (mortality displacement)
   - Zero/negative → true excess mortality

References:
- Armstrong et al. (2014) Am J Epidemiol - Forward displacement
- Gasparrini et al. (2015) Lancet - Multi-country Collaborative

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
    harvesting_for_region,
    compute_harvesting_ratio,
    pool_region_effects,
    meta_random_effects,
    extract_region_percentiles,
    convert_to_json_serializable,
)

print("="*70)
print("02b: HARVESTING / MORTALITY DISPLACEMENT ANALYSIS (v2)")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PHASE0_RESULTS = os.path.join(BASE_DIR, 'new_analysis', 'phase0_data_prep', 'results')
OUTPUT_DIR = os.path.join(BASE_DIR, 'new_analysis', 'phase2_robustness', 'results')
FIGURE_DIR = os.path.join(BASE_DIR, 'new_analysis', 'phase2_robustness', 'figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

# Extended lag to detect harvesting
MAX_LAG = 35  # Beyond typical 21-day window
LAG_HORIZONS = [7, 14, 21, 28, 35]  # Cumulative effect evaluation points
VAR_DF = 4   # Natural spline df for temperature
LAG_DF = 4   # Natural spline df for lag
MIN_OBS = 1000

print(f"\nConfiguration:")
print(f"  Extended max lag: {MAX_LAG} days")
print(f"  Lag horizons: {LAG_HORIZONS}")
print(f"  Spline df: var={VAR_DF}, lag={LAG_DF}")
print(f"  Min observations per region: {MIN_OBS}")

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

parser = argparse.ArgumentParser(description='Harvesting analysis for DLNM')
parser.add_argument('--quick', action='store_true', help='Run on subset of regions for testing')
parser.add_argument('--level', type=str, default='intermediate',
                    choices=['intermediate', 'immediate'],
                    help='Spatial level to analyze')
parser.add_argument('--save-figures', action='store_true', default=True, help='Save figures')
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
    valid_regions = valid_regions[:10]
    df = df[df['region_code'].isin(valid_regions)]
    print(f"\n[QUICK MODE] Running on {len(valid_regions)} regions only")

print(f"\nFinal dataset: {len(df):,} rows, {len(valid_regions)} regions")

# =============================================================================
# RUN HARVESTING ANALYSIS PER REGION
# =============================================================================

print("\n" + "-"*70)
print("Running Harvesting Analysis per Region")
print("-"*70)

region_results = []

for i, region in enumerate(valid_regions):
    if (i + 1) % 10 == 0 or i == 0:
        print(f"  Processing region {i+1}/{len(valid_regions)}: {region}")
    
    df_region = df[df['region_code'] == region].copy()
    
    result = harvesting_for_region(
        df_region,
        max_lag=MAX_LAG,
        horizons=LAG_HORIZONS,
        var_df=VAR_DF,
        lag_df=LAG_DF,
        temp_col='temp_mean',
        deaths_col='deaths_elderly',
        pop_col='pop_elderly',
        family='quasi-poisson',
    )
    
    if result is not None:
        region_results.append(result)

print(f"\nSuccessfully fitted {len(region_results)}/{len(valid_regions)} regions")

# =============================================================================
# META-ANALYSIS AT EACH HORIZON
# =============================================================================

print("\n" + "-"*70)
print("Meta-Analysis at Each Horizon")
print("-"*70)

def meta_analyze_horizon(region_results, horizon, percentile='p99'):
    """Pool log_rr at given horizon across regions."""
    effects = []
    variances = []
    region_names = []
    
    for r in region_results:
        try:
            rr_data = r['rrs_by_horizon'][horizon][percentile]
            log_rr = rr_data['log_rr']
            se = rr_data['se']
            
            if not np.isnan(log_rr) and not np.isnan(se) and se > 0:
                effects.append(log_rr)
                variances.append(se ** 2)
                region_names.append(r['region'])
        except (KeyError, TypeError):
            continue
    
    if len(effects) == 0:
        return None
    
    meta = meta_random_effects(np.array(effects), np.array(variances))
    
    return {
        'horizon': horizon,
        'percentile': percentile,
        'pooled_log_rr': meta['effect'],
        'pooled_se': meta['se'],
        'pooled_rr': np.exp(meta['effect']),
        'pooled_ci_low': np.exp(meta['effect'] - 1.96 * meta['se']),
        'pooled_ci_high': np.exp(meta['effect'] + 1.96 * meta['se']),
        'tau2': meta['tau2'],
        'I2': meta['I2'],
        'k': meta['k'],
        'p': meta['p'],
        'pooled_err': np.exp(meta['effect']) - 1,  # Excess Relative Risk
    }


# Analyze heat (P99 and P97.5)
heat_results = []
for pct in ['p99', 'p975']:
    for h in LAG_HORIZONS:
        result = meta_analyze_horizon(region_results, h, pct)
        if result:
            heat_results.append(result)
            if pct == 'p99':
                rr = result['pooled_rr']
                ci = (result['pooled_ci_low'], result['pooled_ci_high'])
                print(f"  Heat {pct.upper()} @ lag {h:2d}: RR={rr:.3f} ({ci[0]:.3f}-{ci[1]:.3f}), I²={result['I2']:.1f}%")

# Analyze cold (P01 and P025)
cold_results = []
for pct in ['p01', 'p025']:
    for h in LAG_HORIZONS:
        result = meta_analyze_horizon(region_results, h, pct)
        if result:
            cold_results.append(result)
            if pct == 'p01':
                rr = result['pooled_rr']
                ci = (result['pooled_ci_low'], result['pooled_ci_high'])
                print(f"  Cold {pct.upper()} @ lag {h:2d}: RR={rr:.3f} ({ci[0]:.3f}-{ci[1]:.3f}), I²={result['I2']:.1f}%")

# =============================================================================
# COMPUTE HARVESTING RATIOS
# =============================================================================

print("\n" + "-"*70)
print("Computing Harvesting Ratios")
print("-"*70)

def compute_meta_harvesting(meta_results, short_horizon=7, long_horizon=35, pct='p99'):
    """Compute harvesting ratio from pooled ERR at short vs long horizon."""
    short = next((r for r in meta_results if r['horizon'] == short_horizon and r['percentile'] == pct), None)
    long = next((r for r in meta_results if r['horizon'] == long_horizon and r['percentile'] == pct), None)
    
    if short is None or long is None:
        return None
    
    err_short = short['pooled_err']
    err_long = long['pooled_err']
    
    if abs(err_short) < 0.001:
        return {'ratio': np.nan, 'err_short': err_short, 'err_long': err_long}
    
    ratio = 1 - (err_long / err_short)
    
    # Approximate CI via bootstrap in real implementation
    # For now, just point estimate
    return {
        'ratio': ratio,
        'rr_short': short['pooled_rr'],
        'rr_long': long['pooled_rr'],
        'err_short': err_short,
        'err_long': err_long,
        'ci_short': (short['pooled_ci_low'], short['pooled_ci_high']),
        'ci_long': (long['pooled_ci_low'], long['pooled_ci_high']),
    }


heat_harvesting = {}
cold_harvesting = {}

for pct in ['p99', 'p975']:
    hr = compute_meta_harvesting(heat_results, 7, 35, pct)
    if hr:
        heat_harvesting[pct] = hr
        print(f"\nHeat harvesting ({pct}):")
        print(f"  ERR @ lag 7:  {hr['err_short']:.4f}")
        print(f"  ERR @ lag 35: {hr['err_long']:.4f}")
        print(f"  Ratio: {hr['ratio']:.3f}")
        if hr['ratio'] > 0.3:
            print(f"  → Substantial harvesting detected")
        elif hr['ratio'] > 0.1:
            print(f"  → Moderate harvesting detected")
        else:
            print(f"  → Minimal harvesting (mostly true excess)")

for pct in ['p01', 'p025']:
    hr = compute_meta_harvesting(cold_results, 7, 35, pct)
    if hr:
        cold_harvesting[pct] = hr
        print(f"\nCold harvesting ({pct}):")
        print(f"  ERR @ lag 7:  {hr['err_short']:.4f}")
        print(f"  ERR @ lag 35: {hr['err_long']:.4f}")
        print(f"  Ratio: {hr['ratio']:.3f}")
        if hr['ratio'] < -0.3:
            print(f"  → Delayed cold effects (ERR increases over time)")
        elif hr['ratio'] > 0.3:
            print(f"  → Cold harvesting detected")
        else:
            print(f"  → Cold effects mostly persistent")

# =============================================================================
# VISUALIZATIONS
# =============================================================================

if args.save_figures:
    print("\n" + "-"*70)
    print("Creating Visualizations")
    print("-"*70)
    
    # Figure 1: Cumulative RR curves by horizon (heat and cold)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Heat panel
    ax = axes[0]
    for pct, color, label in [('p99', 'darkred', 'P99'), ('p975', 'salmon', 'P97.5')]:
        horizons = []
        rrs = []
        cis = []
        for h in LAG_HORIZONS:
            result = next((r for r in heat_results if r['horizon'] == h and r['percentile'] == pct), None)
            if result:
                horizons.append(h)
                rrs.append(result['pooled_rr'])
                cis.append((result['pooled_ci_low'], result['pooled_ci_high']))
        
        if horizons:
            ax.plot(horizons, rrs, 'o-', color=color, linewidth=2, markersize=8, label=label)
            ci_low = [c[0] for c in cis]
            ci_high = [c[1] for c in cis]
            ax.fill_between(horizons, ci_low, ci_high, color=color, alpha=0.2)
    
    ax.axhline(y=1, color='black', linestyle='--', alpha=0.5)
    ax.set_xlabel('Lag Horizon (days)', fontsize=12)
    ax.set_ylabel('Pooled Cumulative RR', fontsize=12)
    ax.set_title('Heat Effects: Cumulative RR by Horizon', fontsize=14)
    ax.legend(loc='best')
    ax.set_xticks(LAG_HORIZONS)
    ax.grid(True, alpha=0.3)
    
    # Cold panel
    ax = axes[1]
    for pct, color, label in [('p01', 'darkblue', 'P01'), ('p025', 'lightblue', 'P2.5')]:
        horizons = []
        rrs = []
        cis = []
        for h in LAG_HORIZONS:
            result = next((r for r in cold_results if r['horizon'] == h and r['percentile'] == pct), None)
            if result:
                horizons.append(h)
                rrs.append(result['pooled_rr'])
                cis.append((result['pooled_ci_low'], result['pooled_ci_high']))
        
        if horizons:
            ax.plot(horizons, rrs, 'o-', color=color, linewidth=2, markersize=8, label=label)
            ci_low = [c[0] for c in cis]
            ci_high = [c[1] for c in cis]
            ax.fill_between(horizons, ci_low, ci_high, color=color, alpha=0.2)
    
    ax.axhline(y=1, color='black', linestyle='--', alpha=0.5)
    ax.set_xlabel('Lag Horizon (days)', fontsize=12)
    ax.set_ylabel('Pooled Cumulative RR', fontsize=12)
    ax.set_title('Cold Effects: Cumulative RR by Horizon', fontsize=14)
    ax.legend(loc='best')
    ax.set_xticks(LAG_HORIZONS)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    level_tag = '' if level == 'intermediate' else f'_{level}'
    fig_path = os.path.join(FIGURE_DIR, f'harvesting_cumulative_rr_by_horizon{level_tag}.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {fig_path}")
    
    # Figure 2: ERR trajectories 
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for pct, color, label in [('p99', 'red', 'Heat P99'), ('p01', 'blue', 'Cold P01')]:
        results = heat_results if 'p9' in pct else cold_results
        horizons = []
        errs = []
        for h in LAG_HORIZONS:
            result = next((r for r in results if r['horizon'] == h and r['percentile'] == pct), None)
            if result:
                horizons.append(h)
                errs.append(result['pooled_err'] * 100)  # Convert to percentage
        
        if horizons:
            ax.plot(horizons, errs, 'o-', color=color, linewidth=2.5, markersize=10, label=label)
    
    ax.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    ax.set_xlabel('Lag Horizon (days)', fontsize=12)
    ax.set_ylabel('Excess Relative Risk (%)', fontsize=12)
    ax.set_title('Mortality Displacement Assessment:\nERR Trajectory by Lag Horizon', fontsize=14)
    ax.legend(loc='best')
    ax.set_xticks(LAG_HORIZONS)
    ax.grid(True, alpha=0.3)
    
    # Add interpretation annotation
    ax.annotate('↓ Declining heat ERR suggests\nmortality displacement (harvesting)',
                xy=(28, 1), fontsize=10, style='italic', color='red', alpha=0.8)
    ax.annotate('→ Stable/increasing cold ERR suggests\ntrue excess mortality',
                xy=(14, 3), fontsize=10, style='italic', color='blue', alpha=0.8)
    
    plt.tight_layout()
    level_tag = '' if level == 'intermediate' else f'_{level}'
    fig_path = os.path.join(FIGURE_DIR, f'harvesting_err_trajectory{level_tag}.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {fig_path}")

# =============================================================================
# SAVE RESULTS
# =============================================================================

print("\n" + "-"*70)
print("Saving Results")
print("-"*70)

# Save region-level results
region_df_list = []
for r in region_results:
    row = {
        'region_code': r['region'],
        'n_obs': r['n_obs'],
    }
    for h in LAG_HORIZONS:
        for pct in ['p01', 'p025', 'p975', 'p99']:
            try:
                rr_data = r['rrs_by_horizon'][h][pct]
                row[f'rr_{pct}_lag{h}'] = rr_data['rr']
                row[f'rr_{pct}_lag{h}_low'] = rr_data['low']
                row[f'rr_{pct}_lag{h}_high'] = rr_data['high']
            except (KeyError, TypeError):
                pass
    region_df_list.append(row)

region_df = pd.DataFrame(region_df_list)
level_tag = '' if level == 'intermediate' else f'_{level}'
region_path = os.path.join(OUTPUT_DIR, f'harvesting_region_results_v2{level_tag}.csv')
region_df.to_csv(region_path, index=False)
print(f"  Saved: {region_path}")

# Save pooled results
pooled_results = {
    'heat': heat_results,
    'cold': cold_results,
    'heat_harvesting': heat_harvesting,
    'cold_harvesting': cold_harvesting,
    'metadata': {
        'max_lag': MAX_LAG,
        'horizons': LAG_HORIZONS,
        'var_df': VAR_DF,
        'lag_df': LAG_DF,
        'n_regions': len(region_results),
        'level': level,
        'analysis_date': datetime.now().isoformat(),
    }
}

pooled_path = os.path.join(OUTPUT_DIR, f'harvesting_pooled_results_v2{level_tag}.json')
with open(pooled_path, 'w') as f:
    json.dump(pooled_results, f, indent=2, default=float)
print(f"  Saved: {pooled_path}")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("HARVESTING ANALYSIS SUMMARY (v2)")
print("="*70)

print(f"\nRegions analyzed: {len(region_results)}/{len(valid_regions)}")
print(f"Maximum lag: {MAX_LAG} days")

if heat_harvesting.get('p99'):
    hr = heat_harvesting['p99']
    print(f"\nHeat (P99 vs MMT):")
    print(f"  RR @ lag 7:  {hr['rr_short']:.3f} {hr['ci_short']}")
    print(f"  RR @ lag 35: {hr['rr_long']:.3f} {hr['ci_long']}")
    print(f"  Harvesting ratio: {hr['ratio']:.3f}")
    
if cold_harvesting.get('p01'):
    hr = cold_harvesting['p01']
    print(f"\nCold (P01 vs MMT):")
    print(f"  RR @ lag 7:  {hr['rr_short']:.3f} {hr['ci_short']}")
    print(f"  RR @ lag 35: {hr['rr_long']:.3f} {hr['ci_long']}")
    print(f"  Harvesting ratio: {hr['ratio']:.3f}")

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*70)
