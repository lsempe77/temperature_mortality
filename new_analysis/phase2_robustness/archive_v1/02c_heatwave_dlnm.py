"""
02c_heatwave_dlnm.py
=====================
Heatwave Effect Modification Analysis

Tests whether sustained extreme temperatures (heatwaves) have ADDITIONAL
mortality effects beyond the daily temperature-mortality relationship.

Key Question: Does a day of 35°C have a stronger effect when it's part 
              of a multi-day heatwave vs. an isolated hot day?

Method:
1. Define heatwaves: ≥2 consecutive days with temp > region-specific P90
2. Fit DLNM with heatwave interaction term
3. Compare mortality effects on heatwave vs non-heatwave extreme days

This is an EFFECT MODIFICATION analysis - heatwave status modifies 
the temperature-mortality relationship.

References:
- Gasparrini & Armstrong (2011) - Heatwave effect modification
- Anderson & Bell (2009) - Heatwave definition sensitivity

Author: Climate-Health Analysis Pipeline
Date: December 2025
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod.families import Poisson
import os
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("02c: HEATWAVE EFFECT MODIFICATION ANALYSIS")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PHASE0_RESULTS = os.path.join(BASE_DIR, 'new_analysis', 'phase0_data_prep', 'results')
OUTPUT_DIR = os.path.join(BASE_DIR, 'new_analysis', 'phase2_robustness', 'results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# DLNM parameters
MAX_LAG = 21
POLY_DEGREE = 3
MIN_OBS = 1000

# Heatwave definitions to test
HEATWAVE_DEFINITIONS = [
    {'name': 'strict', 'threshold_pct': 95, 'min_days': 3},
    {'name': 'moderate', 'threshold_pct': 90, 'min_days': 2},
    {'name': 'lenient', 'threshold_pct': 85, 'min_days': 2},
]

print(f"\nConfiguration:")
print(f"  Max lag: {MAX_LAG} days")
print(f"  Heatwave definitions: {len(HEATWAVE_DEFINITIONS)}")

# =============================================================================
# DATA LOADING
# =============================================================================

print("\n" + "-"*70)
print("Loading Data")
print("-"*70)

# Load temperature data
df_temp = pd.read_parquet(os.path.join(PHASE0_RESULTS, 'era5_intermediate_daily.parquet'))
df_temp['date'] = pd.to_datetime(df_temp['date'])
print(f"ERA5: {len(df_temp):,} rows, {df_temp['region_code'].nunique()} regions")

# Load mortality data
df_mort = pd.read_parquet(os.path.join(PHASE0_RESULTS, 'mortality_regional_daily_elderly.parquet'))
df_mort['date'] = pd.to_datetime(df_mort['date'])
print(f"Mortality: {len(df_mort):,} rows")

# Load SES data for population
df_ses = pd.read_csv(os.path.join(PHASE0_RESULTS, 'ses_intermediate_covariates.csv'))
pop_map = dict(zip(df_ses['intermediate_code'], df_ses['pop_elderly']))

# Load holidays
df_holidays = pd.read_parquet(os.path.join(PHASE0_RESULTS, 'brazilian_holidays_daily.parquet'))
df_holidays['date'] = pd.to_datetime(df_holidays['date'])

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
valid_regions = region_counts[region_counts >= MIN_OBS].index
df = df[df['region_code'].isin(valid_regions)]

print(f"\nFinal dataset: {len(df):,} rows, {len(valid_regions)} regions")

# =============================================================================
# DEFINE HEATWAVES
# =============================================================================

print("\n" + "-"*70)
print("Identifying Heatwaves")
print("-"*70)

# Sort by region and date
df = df.sort_values(['region_code', 'date']).reset_index(drop=True)

# Calculate region-specific temperature thresholds
region_temp_stats = df.groupby('region_code')['temp_mean'].agg(['mean', 'std'])
region_temp_stats['p85'] = df.groupby('region_code')['temp_mean'].quantile(0.85)
region_temp_stats['p90'] = df.groupby('region_code')['temp_mean'].quantile(0.90)
region_temp_stats['p95'] = df.groupby('region_code')['temp_mean'].quantile(0.95)

print(f"  Region-specific thresholds calculated")
print(f"  National P90: {df['temp_mean'].quantile(0.90):.1f}°C")
print(f"  National P95: {df['temp_mean'].quantile(0.95):.1f}°C")

def identify_heatwaves(df, threshold_col, min_days):
    """
    Identify heatwave days: consecutive days above threshold.
    
    Returns boolean series where True = day is part of a heatwave.
    """
    df = df.copy()
    
    # Mark days above threshold
    df['above_threshold'] = df.apply(
        lambda x: x['temp_mean'] > region_temp_stats.loc[x['region_code'], threshold_col],
        axis=1
    )
    
    # Identify consecutive runs
    df['run_start'] = (df['above_threshold'] != df.groupby('region_code')['above_threshold'].shift()).cumsum()
    
    # Count run lengths
    run_lengths = df[df['above_threshold']].groupby(['region_code', 'run_start']).size()
    
    # Mark days that are part of runs >= min_days
    df['heatwave'] = False
    for (region, run_start), length in run_lengths.items():
        if length >= min_days:
            mask = (df['region_code'] == region) & (df['run_start'] == run_start) & df['above_threshold']
            df.loc[mask, 'heatwave'] = True
    
    return df['heatwave']

# Apply each heatwave definition
for hw_def in HEATWAVE_DEFINITIONS:
    col_name = f"heatwave_{hw_def['name']}"
    threshold_col = f"p{hw_def['threshold_pct']}"
    df[col_name] = identify_heatwaves(df, threshold_col, hw_def['min_days'])
    n_days = df[col_name].sum()
    n_events = df[df[col_name]].groupby(['region_code']).apply(
        lambda x: (x[col_name] != x[col_name].shift()).sum() // 2
    ).sum()
    print(f"  {hw_def['name'].capitalize()}: {n_days:,} days ({n_days/len(df)*100:.2f}%), ~{n_events:,} events")

# =============================================================================
# CREATE LAGGED TEMPERATURES
# =============================================================================

print("\n" + "-"*70)
print("Creating Lag Matrix")
print("-"*70)

for lag in range(1, MAX_LAG + 1):
    df[f'temp_lag{lag}'] = df.groupby('region_code')['temp_mean'].shift(lag)

print(f"  Created lags 1-{MAX_LAG}")

# =============================================================================
# DLNM FUNCTIONS WITH HEATWAVE INTERACTION
# =============================================================================

def fit_dlnm_with_heatwave(df_input, heatwave_col, max_lag=21, poly_degree=3):
    """
    Fit DLNM with heatwave interaction term.
    
    Model: log(deaths) = baseline_temp_effect + heatwave_main + temp×heatwave_interaction
    
    Returns:
    - Base temperature effect (non-heatwave days)
    - Heatwave main effect
    - Heatwave modification (additional effect on extreme hot days during heatwave)
    """
    df = df_input.copy()
    
    # Temperature stats
    temp_center = df['temp_mean'].median()
    temp_scale = df['temp_mean'].std()
    
    # Standardize temperatures
    df['temp_std'] = (df['temp_mean'] - temp_center) / temp_scale
    for lag in range(1, max_lag + 1):
        df[f'temp_lag{lag}_std'] = (df[f'temp_lag{lag}'] - temp_center) / temp_scale
    
    # Build cross-basis for temperature
    cb_list = []
    cb_names = []
    
    for lag in range(max_lag + 1):
        if lag == 0:
            temp_col = 'temp_std'
        else:
            temp_col = f'temp_lag{lag}_std'
        
        for p in range(1, poly_degree + 1):
            col_name = f'temp_pow{p}_lag{lag}'
            df[col_name] = df[temp_col].values ** p
            cb_list.append(col_name)
            cb_names.append(col_name)
    
    # Heatwave indicator
    df['heatwave'] = df[heatwave_col].astype(int)
    
    # Interaction terms: heatwave × temperature cross-basis
    interaction_names = []
    for cb_name in cb_names:
        int_name = f'{cb_name}_hw_int'
        df[int_name] = df[cb_name] * df['heatwave']
        interaction_names.append(int_name)
    
    # Remove rows with NaN
    all_cols = cb_names + interaction_names
    valid_mask = ~df[all_cols].isna().any(axis=1)
    df_valid = df[valid_mask].copy()
    
    if len(df_valid) < 1000:
        return None
    
    # Build design matrix
    X = pd.get_dummies(
        df_valid[['day_of_week', 'month', 'year', 'region_code']],
        columns=['day_of_week', 'month', 'year', 'region_code'],
        drop_first=True
    )
    
    # Add temperature cross-basis
    for name in cb_names:
        X[name] = df_valid[name].astype(float)
    
    # Add heatwave main effect
    X['heatwave'] = df_valid['heatwave']
    
    # Add interactions
    for name in interaction_names:
        X[name] = df_valid[name].astype(float)
    
    X['is_holiday'] = df_valid['is_holiday']
    X = sm.add_constant(X)
    
    y = df_valid['deaths_elderly'].values.astype(float)
    
    # Fit model
    try:
        model = GLM(y, X.astype(float), family=Poisson())
        result = model.fit(scale='X2', cov_type='HC1')
    except Exception as e:
        print(f"    Model fitting failed: {e}")
        return None
    
    # Extract key coefficients
    hw_idx = list(X.columns).index('heatwave')
    hw_coef = result.params.iloc[hw_idx]
    hw_se = result.bse.iloc[hw_idx]
    hw_rr = np.exp(hw_coef)
    hw_rr_lower = np.exp(hw_coef - 1.96 * hw_se)
    hw_rr_upper = np.exp(hw_coef + 1.96 * hw_se)
    
    # Get interaction coefficient sums (total additional effect during heatwave)
    int_coef_sum = 0
    int_var_sum = 0
    for name in interaction_names:
        if name in X.columns:
            idx = list(X.columns).index(name)
            int_coef_sum += result.params.iloc[idx]
            int_var_sum += result.bse.iloc[idx] ** 2
    
    # Calculate RR at P97.5 and P99 on heatwave vs non-heatwave days
    temp_p97_5 = df['temp_mean'].quantile(0.975)
    temp_p99 = df['temp_mean'].quantile(0.99)
    temp_p50 = df['temp_mean'].quantile(0.50)
    
    # Calculate for both thresholds
    def calc_heat_rr(target_temp):
        target_std = (target_temp - temp_center) / temp_scale
        ref_std = (temp_p50 - temp_center) / temp_scale
        
        # Base effect (non-heatwave)
        log_rr_base = 0
        var_base = 0
        for lag in range(max_lag + 1):
            for p in range(1, poly_degree + 1):
                cb_name = f'temp_pow{p}_lag{lag}'
                if cb_name in X.columns:
                    idx = list(X.columns).index(cb_name)
                    coef = result.params.iloc[idx]
                    se = result.bse.iloc[idx]
                    diff = target_std ** p - ref_std ** p
                    log_rr_base += coef * diff
                    var_base += (se * diff) ** 2
        
        rr_base = np.exp(log_rr_base)
        rr_base_lower = np.exp(log_rr_base - 1.96 * np.sqrt(var_base))
        rr_base_upper = np.exp(log_rr_base + 1.96 * np.sqrt(var_base))
        
        # Heatwave additional effect
        log_rr_int = 0
        var_int = 0
        for lag in range(max_lag + 1):
            for p in range(1, poly_degree + 1):
                int_name = f'temp_pow{p}_lag{lag}_hw_int'
                if int_name in X.columns:
                    idx = list(X.columns).index(int_name)
                    coef = result.params.iloc[idx]
                    se = result.bse.iloc[idx]
                    diff = target_std ** p - ref_std ** p
                    log_rr_int += coef * diff
                    var_int += (se * diff) ** 2
        
        # Total effect during heatwave = base + main + interaction
        log_rr_heatwave = log_rr_base + hw_coef + log_rr_int
        rr_heatwave = np.exp(log_rr_heatwave)
        
        # Combined variance (simplified)
        var_combined = var_base + hw_se**2 + var_int
        rr_heatwave_lower = np.exp(log_rr_heatwave - 1.96 * np.sqrt(var_combined))
        rr_heatwave_upper = np.exp(log_rr_heatwave + 1.96 * np.sqrt(var_combined))
        
        return {
            'non_heatwave': {'rr': float(rr_base), 'rr_lower': float(rr_base_lower), 'rr_upper': float(rr_base_upper)},
            'during_heatwave': {'rr': float(rr_heatwave), 'rr_lower': float(rr_heatwave_lower), 'rr_upper': float(rr_heatwave_upper)}
        }
    
    # Calculate for P97.5 (primary) and P99 (comparison)
    rr_p97_5 = calc_heat_rr(temp_p97_5)
    rr_p99 = calc_heat_rr(temp_p99)
    
    return {
        'n_obs': len(df_valid),
        'n_heatwave_days': int(df_valid['heatwave'].sum()),
        'pct_heatwave_days': float(df_valid['heatwave'].mean() * 100),
        'aic': result.aic,
        'heatwave_main_effect': {
            'rr': float(hw_rr),
            'rr_lower': float(hw_rr_lower),
            'rr_upper': float(hw_rr_upper),
            'significant': bool(hw_rr_lower > 1 or hw_rr_upper < 1)
        },
        'rr_p97_5_non_heatwave': rr_p97_5['non_heatwave'],
        'rr_p97_5_during_heatwave': rr_p97_5['during_heatwave'],
        'rr_p99_non_heatwave': rr_p99['non_heatwave'],
        'rr_p99_during_heatwave': rr_p99['during_heatwave'],
        'effect_modification_ratio_p97_5': float(rr_p97_5['during_heatwave']['rr'] / rr_p97_5['non_heatwave']['rr']) if rr_p97_5['non_heatwave']['rr'] > 0 else None,
        'effect_modification_ratio_p99': float(rr_p99['during_heatwave']['rr'] / rr_p99['non_heatwave']['rr']) if rr_p99['non_heatwave']['rr'] > 0 else None,
        'temp_center': float(temp_center),
        'temp_scale': float(temp_scale)
    }

# =============================================================================
# RUN ANALYSES FOR EACH HEATWAVE DEFINITION
# =============================================================================

print("\n" + "="*70)
print("HEATWAVE EFFECT MODIFICATION ANALYSIS")
print("="*70)

results = {}

for hw_def in HEATWAVE_DEFINITIONS:
    name = hw_def['name']
    col_name = f"heatwave_{name}"
    
    print(f"\n--- {name.upper()} DEFINITION ---")
    print(f"  Threshold: P{hw_def['threshold_pct']}, Min days: {hw_def['min_days']}")
    
    result = fit_dlnm_with_heatwave(df, col_name, max_lag=MAX_LAG, poly_degree=POLY_DEGREE)
    
    if result:
        results[name] = result
        print(f"\n  Heatwave days: {result['n_heatwave_days']:,} ({result['pct_heatwave_days']:.2f}%)")
        print(f"  Heatwave main effect: RR = {result['heatwave_main_effect']['rr']:.4f} "
              f"[{result['heatwave_main_effect']['rr_lower']:.4f}-{result['heatwave_main_effect']['rr_upper']:.4f}]")
        print(f"  RR at P97.5 (non-heatwave): {result['rr_p97_5_non_heatwave']['rr']:.4f} <- PRIMARY")
        print(f"  RR at P97.5 (during heatwave): {result['rr_p97_5_during_heatwave']['rr']:.4f}")
        print(f"  RR at P99 (non-heatwave): {result['rr_p99_non_heatwave']['rr']:.4f}")
        print(f"  RR at P99 (during heatwave): {result['rr_p99_during_heatwave']['rr']:.4f}")
        print(f"  Effect modification ratio (P97.5): {result['effect_modification_ratio_p97_5']:.2f}x")
        print(f"  Effect modification ratio (P99): {result['effect_modification_ratio_p99']:.2f}x")
        
        if result['effect_modification_ratio_p97_5'] > 1.2:
            print(f"  → Heatwaves AMPLIFY heat effects by {(result['effect_modification_ratio_p97_5']-1)*100:.0f}%")
        elif result['effect_modification_ratio_p97_5'] < 0.8:
            print(f"  → Heatwaves ATTENUATE heat effects")
        else:
            print(f"  → No significant effect modification")
    else:
        print("  Model fitting failed")

# =============================================================================
# SAVE RESULTS
# =============================================================================

print("\n" + "="*70)
print("SAVING RESULTS")
print("="*70)

output_results = {
    'analysis_date': datetime.now().isoformat(),
    'configuration': {
        'max_lag': MAX_LAG,
        'poly_degree': POLY_DEGREE,
        'n_regions': int(len(valid_regions)),
        'definitions_tested': HEATWAVE_DEFINITIONS
    },
    'temperature_percentiles': {
        'p50': float(df['temp_mean'].quantile(0.50)),
        'p85': float(df['temp_mean'].quantile(0.85)),
        'p90': float(df['temp_mean'].quantile(0.90)),
        'p95': float(df['temp_mean'].quantile(0.95)),
        'p99': float(df['temp_mean'].quantile(0.99))
    },
    'results_by_definition': results
}

# Save results
output_file = os.path.join(OUTPUT_DIR, 'heatwave_effect_modification.json')
with open(output_file, 'w') as f:
    json.dump(output_results, f, indent=2, default=str)
print(f"\nSaved: {output_file}")

# Create summary table
summary_data = []
for name, res in results.items():
    summary_data.append({
        'definition': name,
        'heatwave_days': res['n_heatwave_days'],
        'pct_days': res['pct_heatwave_days'],
        'heatwave_main_rr': res['heatwave_main_effect']['rr'],
        'heatwave_main_significant': res['heatwave_main_effect']['significant'],
        'rr_p97_5_normal': res['rr_p97_5_non_heatwave']['rr'],
        'rr_p97_5_heatwave': res['rr_p97_5_during_heatwave']['rr'],
        'rr_p99_normal': res['rr_p99_non_heatwave']['rr'],
        'rr_p99_heatwave': res['rr_p99_during_heatwave']['rr'],
        'effect_modification_p97_5': res['effect_modification_ratio_p97_5'],
        'effect_modification_p99': res['effect_modification_ratio_p99']
    })

summary_df = pd.DataFrame(summary_data)
summary_file = os.path.join(OUTPUT_DIR, 'heatwave_summary.csv')
summary_df.to_csv(summary_file, index=False)
print(f"Saved: {summary_file}")

# =============================================================================
# PRINT SUMMARY
# =============================================================================

print("\n" + "="*70)
print("HEATWAVE ANALYSIS SUMMARY")
print("="*70)

print("\n" + summary_df.to_string(index=False))

print("\n" + "-"*70)
print("KEY FINDINGS:")
print("-"*70)

# Find the definition with strongest effect modification
if results:
    max_mod = max(results.items(), key=lambda x: x[1]['effect_modification_ratio_p97_5'])
    min_mod = min(results.items(), key=lambda x: x[1]['effect_modification_ratio_p97_5'])
    
    print(f"\n1. STRONGEST MODIFICATION: {max_mod[0].upper()}")
    print(f"   Effect amplification (P97.5): {(max_mod[1]['effect_modification_ratio_p97_5']-1)*100:.0f}%")
    print(f"   Effect amplification (P99): {(max_mod[1]['effect_modification_ratio_p99']-1)*100:.0f}%")
    
    print(f"\n2. HEATWAVE MAIN EFFECTS:")
    for name, res in results.items():
        sig = "✓" if res['heatwave_main_effect']['significant'] else "✗"
        print(f"   {name}: RR = {res['heatwave_main_effect']['rr']:.3f} [{sig}]")
    
    # Check if heatwaves consistently amplify effects
    all_amplify = all(r['effect_modification_ratio_p97_5'] > 1.1 for r in results.values())
    if all_amplify:
        print(f"\n3. CONCLUSION: Heatwaves consistently AMPLIFY temperature effects")
        print(f"   Policy implication: Heat warnings during multi-day events are critical")
    else:
        print(f"\n3. CONCLUSION: Mixed evidence for heatwave effect modification")

print("\n" + "="*70)
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*70)
