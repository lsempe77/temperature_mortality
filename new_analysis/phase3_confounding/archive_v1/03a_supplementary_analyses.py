"""
03a_supplementary_analyses.py
==============================
Supplementary Analyses for DLNM Results

These are NOT primary results but SUPPLEMENTARY material showing:
1. Apparent temperature (humidity-adjusted) vs dry-bulb
2. Pollution-adjusted model (PM2.5, O3 controls)
3. Influenza-adjusted model (flu season control)
4. Holiday-adjusted model (holiday effects)

NOTE: Pollution adjustment is SUPPLEMENTARY because pollution may be a 
MEDIATOR of heat effects (heat → stagnation → pollution), not just a 
confounder. Over-adjustment would bias results.

Input: Regional data from phase0_data_prep/results/
Output: Supplementary analysis results for comparison

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
print("PHASE 3: SUPPLEMENTARY ANALYSES")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PHASE0_RESULTS = os.path.join(BASE_DIR, 'new_analysis', 'phase0_data_prep', 'results')
OUTPUT_DIR = os.path.join(BASE_DIR, 'new_analysis', 'phase3_confounding', 'results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# DLNM parameters (same as Phase 1/2)
MAX_LAG = 21
POLY_DEGREE = 3
MIN_OBS = 1000

print(f"\nConfiguration:")
print(f"  Max lag: {MAX_LAG} days")
print(f"  Polynomial degree: {POLY_DEGREE}")

# =============================================================================
# DATA LOADING
# =============================================================================

print("\n" + "-"*70)
print("Loading Data")
print("-"*70)

# Load temperature data (includes dewpoint for apparent temp)
df_temp = pd.read_parquet(os.path.join(PHASE0_RESULTS, 'era5_intermediate_daily.parquet'))
df_temp['date'] = pd.to_datetime(df_temp['date'])
print(f"ERA5: {len(df_temp):,} rows, {df_temp['region_code'].nunique()} regions")
print(f"  Columns: {list(df_temp.columns)}")

# Load pollution data
df_poll = pd.read_parquet(os.path.join(PHASE0_RESULTS, 'cams_intermediate_daily.parquet'))
df_poll['date'] = pd.to_datetime(df_poll['date'])
print(f"CAMS Pollution: {len(df_poll):,} rows")
print(f"  Columns: {list(df_poll.columns)}")

# Load influenza data
df_flu = pd.read_parquet(os.path.join(PHASE0_RESULTS, 'influenza_daily_by_intermediate_region.parquet'))
df_flu['date'] = pd.to_datetime(df_flu['date'])
print(f"Influenza: {len(df_flu):,} rows")
print(f"  Columns: {list(df_flu.columns)}")

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
print(f"Holidays: {len(df_holidays)} days")

# =============================================================================
# MERGE DATA
# =============================================================================

print("\n" + "-"*70)
print("Merging Data")
print("-"*70)

# Start with mortality
df = df_mort.copy()

# Merge temperature
df = pd.merge(
    df,
    df_temp[['date', 'region_code', 'temp_mean', 'dewpoint_mean']],
    on=['date', 'region_code'],
    how='inner'
)
print(f"After temp merge: {len(df):,} rows")

# Merge pollution
# Handle different column naming in CAMS data
if 'intermediate_code' in df_poll.columns:
    df_poll = df_poll.rename(columns={'intermediate_code': 'region_code'})
if 'ozone' in df_poll.columns:
    df_poll = df_poll.rename(columns={'ozone': 'o3'})

poll_cols = [c for c in df_poll.columns if c not in ['date', 'region_code']]
df = pd.merge(
    df,
    df_poll[['date', 'region_code'] + poll_cols],
    on=['date', 'region_code'],
    how='left'
)
print(f"After pollution merge: {len(df):,} rows")

# Merge influenza
flu_cols = [c for c in df_flu.columns if c not in ['date', 'region_code', 'intermediate_code']]
# Rename region column if needed
if 'intermediate_code' in df_flu.columns:
    df_flu = df_flu.rename(columns={'intermediate_code': 'region_code'})
df = pd.merge(
    df,
    df_flu[['date', 'region_code'] + flu_cols] if 'region_code' in df_flu.columns else df_flu,
    on=['date', 'region_code'],
    how='left'
)
print(f"After influenza merge: {len(df):,} rows")

# Merge holidays
df = pd.merge(
    df,
    df_holidays[['date', 'is_holiday', 'is_holiday_week']],
    on='date',
    how='left'
)
df['is_holiday'] = df['is_holiday'].fillna(0).astype(int)

# Add population and time variables
df['pop_elderly'] = df['region_code'].map(pop_map)
df['year'] = df['date'].dt.year
df['month'] = df['date'].dt.month
df['day_of_week'] = df['date'].dt.dayofweek

# Calculate apparent temperature
# AT = Ta + 0.33 × e - 0.70 × ws - 4.00 (where e = vapor pressure from dewpoint)
# Simplified: AT ≈ Ta + 0.33 × (6.105 × exp((17.27 × Td)/(237.7 + Td))) - 4.0
# Or use heat index approximation
if 'dewpoint_mean' in df.columns:
    # Convert dewpoint to relative humidity approximation
    # Then calculate apparent temperature
    Td = df['dewpoint_mean'].values
    Ta = df['temp_mean'].values
    
    # Vapor pressure from dewpoint (Magnus formula)
    e = 6.112 * np.exp((17.67 * Td) / (Td + 243.5))
    
    # Apparent temperature (simplified)
    df['apparent_temp'] = Ta + 0.33 * e - 4.0
    print(f"Apparent temperature calculated")
else:
    df['apparent_temp'] = df['temp_mean']
    print(f"WARNING: No dewpoint data, using dry-bulb as apparent temp")

# Fill missing pollution/flu with 0 or interpolate
for col in ['pm25', 'pm10', 'o3', 'srag_cases', 'confirmed_flu']:
    if col in df.columns:
        df[col] = df[col].fillna(0)

print(f"\nFinal dataset: {len(df):,} rows")
print(f"Missing values:")
for col in ['temp_mean', 'apparent_temp', 'pm25', 'srag_cases', 'deaths_elderly']:
    if col in df.columns:
        print(f"  {col}: {df[col].isna().sum():,}")

# Filter regions with enough observations
df = df.dropna(subset=['temp_mean', 'deaths_elderly', 'pop_elderly'])
region_counts = df.groupby('region_code').size()
valid_regions = region_counts[region_counts >= MIN_OBS].index
df = df[df['region_code'].isin(valid_regions)]
print(f"\nAfter filtering: {len(df):,} rows, {len(valid_regions)} regions")

# Sort and create lags
df = df.sort_values(['region_code', 'date']).reset_index(drop=True)

# Create lagged temperatures
for lag in range(1, MAX_LAG + 1):
    df[f'temp_lag{lag}'] = df.groupby('region_code')['temp_mean'].shift(lag)
    df[f'apparent_lag{lag}'] = df.groupby('region_code')['apparent_temp'].shift(lag)

# Temperature percentiles (including P2.5/P97.5 for consistency)
temp_p1 = df['temp_mean'].quantile(0.01)
temp_p2_5 = df['temp_mean'].quantile(0.025)
temp_p5 = df['temp_mean'].quantile(0.05)
temp_p50 = df['temp_mean'].quantile(0.50)
temp_p95 = df['temp_mean'].quantile(0.95)
temp_p97_5 = df['temp_mean'].quantile(0.975)
temp_p99 = df['temp_mean'].quantile(0.99)

print(f"\nDry-bulb percentiles:")
print(f"  P1={temp_p1:.1f}°C, P2.5={temp_p2_5:.1f}°C, P50={temp_p50:.1f}°C, P97.5={temp_p97_5:.1f}°C, P99={temp_p99:.1f}°C")

# Apparent temperature percentiles
app_p1 = df['apparent_temp'].quantile(0.01)
app_p50 = df['apparent_temp'].quantile(0.50)
app_p99 = df['apparent_temp'].quantile(0.99)
print(f"Apparent temp percentiles: P1={app_p1:.1f}°C, P50={app_p50:.1f}°C, P99={app_p99:.1f}°C")

# =============================================================================
# DLNM HELPER FUNCTIONS
# =============================================================================

def fit_dlnm_model(df_input, temp_col='temp_mean', extra_controls=None, max_lag=21, poly_degree=3):
    """
    Fit DLNM with specified temperature variable and optional extra controls.
    
    Returns RR at P99 and P1 relative to P50.
    """
    df = df_input.copy()
    
    # Temperature stats
    temp_center = df[temp_col].median()
    temp_scale = df[temp_col].std()
    
    # Temperature percentiles for this variable
    t_p1 = df[temp_col].quantile(0.01)
    t_p2_5 = df[temp_col].quantile(0.025)
    t_p50 = df[temp_col].quantile(0.50)
    t_p97_5 = df[temp_col].quantile(0.975)
    t_p99 = df[temp_col].quantile(0.99)
    
    # Standardize temperatures
    df['temp_std'] = (df[temp_col] - temp_center) / temp_scale
    
    # Create lags if not already present
    for lag in range(1, max_lag + 1):
        lag_col = f'{temp_col}_lag{lag}' if temp_col != 'temp_mean' else f'temp_lag{lag}'
        if temp_col == 'apparent_temp':
            lag_col = f'apparent_lag{lag}'
        if lag_col in df.columns:
            df[f'temp_std_lag{lag}'] = (df[lag_col] - temp_center) / temp_scale
        else:
            df[f'temp_std_lag{lag}'] = df.groupby('region_code')['temp_std'].shift(lag)
    
    # Build cross-basis
    cb_list = []
    cb_names = []
    
    for lag in range(max_lag + 1):
        if lag == 0:
            temp_col_std = 'temp_std'
        else:
            temp_col_std = f'temp_std_lag{lag}'
        
        for p in range(1, poly_degree + 1):
            col_name = f'cb_lag{lag}_pow{p}'
            df[col_name] = df[temp_col_std].values ** p
            cb_names.append(col_name)
    
    # Remove rows with NaN
    valid_mask = ~df[cb_names].isna().any(axis=1)
    df_valid = df[valid_mask].copy()
    
    if len(df_valid) < 1000:
        return None
    
    # Build design matrix
    X = pd.get_dummies(
        df_valid[['day_of_week', 'month', 'year', 'region_code']],
        columns=['day_of_week', 'month', 'year', 'region_code'],
        drop_first=True
    )
    
    for name in cb_names:
        X[name] = df_valid[name].astype(float)
    
    X['is_holiday'] = df_valid['is_holiday']
    
    # Add extra controls if specified
    if extra_controls:
        for ctrl_name in extra_controls:
            if ctrl_name in df_valid.columns:
                X[ctrl_name] = df_valid[ctrl_name].astype(float)
    
    X = sm.add_constant(X)
    y = df_valid['deaths_elderly'].values.astype(float)
    
    # Fit model
    try:
        model = GLM(y, X.astype(float), family=Poisson())
        result = model.fit(scale='X2', cov_type='HC1')
    except Exception as e:
        print(f"    Model fitting failed: {e}")
        return None
    
    # Calculate RR at key percentiles
    def calc_rr(target_temp, ref_temp):
        target_std = (target_temp - temp_center) / temp_scale
        ref_std = (ref_temp - temp_center) / temp_scale
        
        log_rr = 0
        var_sum = 0
        
        for lag in range(max_lag + 1):
            for p in range(1, poly_degree + 1):
                col_name = f'cb_lag{lag}_pow{p}'
                if col_name in X.columns:
                    idx = list(X.columns).index(col_name)
                    coef = result.params.iloc[idx]
                    se = result.bse.iloc[idx]
                    diff = target_std ** p - ref_std ** p
                    log_rr += coef * diff
                    var_sum += (se * diff) ** 2
        
        rr = np.exp(log_rr)
        se_log_rr = np.sqrt(var_sum)
        rr_lower = np.exp(log_rr - 1.96 * se_log_rr)
        rr_upper = np.exp(log_rr + 1.96 * se_log_rr)
        
        return rr, rr_lower, rr_upper
    
    rr_p99, rr_p99_lo, rr_p99_hi = calc_rr(t_p99, t_p50)
    rr_p97_5, rr_p97_5_lo, rr_p97_5_hi = calc_rr(t_p97_5, t_p50)
    rr_p1, rr_p1_lo, rr_p1_hi = calc_rr(t_p1, t_p50)
    rr_p2_5, rr_p2_5_lo, rr_p2_5_hi = calc_rr(t_p2_5, t_p50)
    
    return {
        'n_obs': len(df_valid),
        'aic': result.aic,
        'bic': result.bic,
        'dispersion': result.scale,
        'temp_variable': temp_col,
        'extra_controls': extra_controls,
        'rr_p97_5': {
            'rr': float(rr_p97_5),
            'rr_lower': float(rr_p97_5_lo),
            'rr_upper': float(rr_p97_5_hi)
        },
        'rr_p99': {
            'rr': float(rr_p99),
            'rr_lower': float(rr_p99_lo),
            'rr_upper': float(rr_p99_hi)
        },
        'rr_p2_5': {
            'rr': float(rr_p2_5),
            'rr_lower': float(rr_p2_5_lo),
            'rr_upper': float(rr_p2_5_hi)
        },
        'rr_p1': {
            'rr': float(rr_p1),
            'rr_lower': float(rr_p1_lo),
            'rr_upper': float(rr_p1_hi)
        },
        'temp_percentiles': {
            'p1': float(t_p1),
            'p2_5': float(t_p2_5),
            'p50': float(t_p50),
            'p97_5': float(t_p97_5),
            'p99': float(t_p99)
        }
    }

# =============================================================================
# ANALYSIS 1: DRY-BULB VS APPARENT TEMPERATURE
# =============================================================================

print("\n" + "="*70)
print("1. DRY-BULB VS APPARENT TEMPERATURE")
print("="*70)

# Baseline: dry-bulb temperature
print("\n  Fitting dry-bulb temperature model...")
result_drybulb = fit_dlnm_model(df, temp_col='temp_mean')

# Apparent temperature
print("  Fitting apparent temperature model...")
result_apparent = fit_dlnm_model(df, temp_col='apparent_temp')

if result_drybulb and result_apparent:
    print(f"\n  DRY-BULB TEMPERATURE:")
    print(f"    Heat RR(P97.5): {result_drybulb['rr_p97_5']['rr']:.4f} [{result_drybulb['rr_p97_5']['rr_lower']:.4f}-{result_drybulb['rr_p97_5']['rr_upper']:.4f}] <- PRIMARY")
    print(f"    Heat RR(P99):   {result_drybulb['rr_p99']['rr']:.4f} [{result_drybulb['rr_p99']['rr_lower']:.4f}-{result_drybulb['rr_p99']['rr_upper']:.4f}]")
    print(f"    Cold RR(P2.5):  {result_drybulb['rr_p2_5']['rr']:.4f} [{result_drybulb['rr_p2_5']['rr_lower']:.4f}-{result_drybulb['rr_p2_5']['rr_upper']:.4f}] <- PRIMARY")
    print(f"    Cold RR(P1):    {result_drybulb['rr_p1']['rr']:.4f} [{result_drybulb['rr_p1']['rr_lower']:.4f}-{result_drybulb['rr_p1']['rr_upper']:.4f}]")
    print(f"    AIC: {result_drybulb['aic']:.1f}")
    
    print(f"\n  APPARENT TEMPERATURE:")
    print(f"    Heat RR(P97.5): {result_apparent['rr_p97_5']['rr']:.4f} [{result_apparent['rr_p97_5']['rr_lower']:.4f}-{result_apparent['rr_p97_5']['rr_upper']:.4f}] <- PRIMARY")
    print(f"    Heat RR(P99):   {result_apparent['rr_p99']['rr']:.4f} [{result_apparent['rr_p99']['rr_lower']:.4f}-{result_apparent['rr_p99']['rr_upper']:.4f}]")
    print(f"    Cold RR(P2.5):  {result_apparent['rr_p2_5']['rr']:.4f} [{result_apparent['rr_p2_5']['rr_lower']:.4f}-{result_apparent['rr_p2_5']['rr_upper']:.4f}] <- PRIMARY")
    print(f"    Cold RR(P1):    {result_apparent['rr_p1']['rr']:.4f} [{result_apparent['rr_p1']['rr_lower']:.4f}-{result_apparent['rr_p1']['rr_upper']:.4f}]")
    print(f"    AIC: {result_apparent['aic']:.1f}")
    
    # Compare using primary thresholds P97.5/P2.5
    heat_diff = (result_apparent['rr_p97_5']['rr'] - result_drybulb['rr_p97_5']['rr']) / result_drybulb['rr_p97_5']['rr'] * 100
    cold_diff = (result_apparent['rr_p2_5']['rr'] - result_drybulb['rr_p2_5']['rr']) / result_drybulb['rr_p2_5']['rr'] * 100
    
    print(f"\n  COMPARISON:")
    print(f"    Heat effect change: {heat_diff:+.1f}%")
    print(f"    Cold effect change: {cold_diff:+.1f}%")
    print(f"    Better fit (lower AIC): {'Apparent' if result_apparent['aic'] < result_drybulb['aic'] else 'Dry-bulb'}")

# =============================================================================
# ANALYSIS 2: POLLUTION-ADJUSTED MODEL
# =============================================================================

print("\n" + "="*70)
print("2. POLLUTION-ADJUSTED MODEL (SUPPLEMENTARY)")
print("="*70)
print("  NOTE: Pollution may be a MEDIATOR of heat effects, not just a confounder.")
print("        These results go in supplementary materials only.")

# Check available pollution variables
poll_vars = []
for var in ['pm25', 'pm10', 'o3']:
    if var in df.columns and df[var].notna().sum() > 10000:
        poll_vars.append(var)
        print(f"    Available: {var} ({df[var].notna().sum():,} non-null)")

if poll_vars:
    # Unadjusted
    print("\n  Fitting unadjusted model...")
    result_unadj = fit_dlnm_model(df, temp_col='temp_mean', extra_controls=None)
    
    # Pollution-adjusted
    print(f"  Fitting pollution-adjusted model (controls: {poll_vars})...")
    result_poll = fit_dlnm_model(df, temp_col='temp_mean', extra_controls=poll_vars)
    
    if result_unadj and result_poll:
        print(f"\n  UNADJUSTED:")
        print(f"    Heat RR(P99): {result_unadj['rr_p99']['rr']:.4f}")
        print(f"    Cold RR(P1):  {result_unadj['rr_p1']['rr']:.4f}")
        
        print(f"\n  POLLUTION-ADJUSTED:")
        print(f"    Heat RR(P99): {result_poll['rr_p99']['rr']:.4f}")
        print(f"    Cold RR(P1):  {result_poll['rr_p1']['rr']:.4f}")
        
        heat_change = (result_poll['rr_p99']['rr'] - result_unadj['rr_p99']['rr']) / result_unadj['rr_p99']['rr'] * 100
        cold_change = (result_poll['rr_p1']['rr'] - result_unadj['rr_p1']['rr']) / result_unadj['rr_p1']['rr'] * 100
        
        print(f"\n  CHANGE AFTER ADJUSTMENT:")
        print(f"    Heat effect: {heat_change:+.1f}%")
        print(f"    Cold effect: {cold_change:+.1f}%")
        
        if abs(heat_change) > 10:
            print("    ⚠️ Substantial change - may indicate confounding OR mediation")
        else:
            print("    ✓ Minimal confounding by pollution")
else:
    result_poll = None
    print("  WARNING: No pollution data available")

# =============================================================================
# ANALYSIS 3: INFLUENZA-ADJUSTED MODEL
# =============================================================================

print("\n" + "="*70)
print("3. INFLUENZA-ADJUSTED MODEL")
print("="*70)
print("  Purpose: Check if 'cold' effect is confounded by winter flu season")

# Check available influenza variables
flu_vars = []
for var in ['srag_cases', 'confirmed_flu', 'srag_elderly']:
    if var in df.columns and df[var].sum() > 1000:
        flu_vars.append(var)
        print(f"    Available: {var} (total: {df[var].sum():,.0f})")

if flu_vars:
    flu_control = flu_vars[0]  # Use first available
    
    print(f"\n  Fitting flu-adjusted model (control: {flu_control})...")
    result_flu = fit_dlnm_model(df, temp_col='temp_mean', extra_controls=[flu_control])
    
    if result_unadj and result_flu:
        print(f"\n  UNADJUSTED:")
        print(f"    Cold RR(P1): {result_unadj['rr_p1']['rr']:.4f}")
        
        print(f"\n  FLU-ADJUSTED:")
        print(f"    Cold RR(P1): {result_flu['rr_p1']['rr']:.4f}")
        
        cold_change = (result_flu['rr_p1']['rr'] - result_unadj['rr_p1']['rr']) / result_unadj['rr_p1']['rr'] * 100
        
        print(f"\n  COLD EFFECT CHANGE: {cold_change:+.1f}%")
        
        if cold_change < -10:
            print("    ⚠️ Cold effect attenuated - some confounding by flu season")
        else:
            print("    ✓ Cold effect robust to flu adjustment")
else:
    result_flu = None
    print("  WARNING: No influenza data available")

# =============================================================================
# ANALYSIS 4: FULLY ADJUSTED MODEL (SUPPLEMENTARY)
# =============================================================================

print("\n" + "="*70)
print("4. FULLY ADJUSTED MODEL (ALL CONTROLS)")
print("="*70)

all_controls = poll_vars + flu_vars if flu_vars else poll_vars
if all_controls:
    print(f"  Controls: {all_controls}")
    
    result_full = fit_dlnm_model(df, temp_col='temp_mean', extra_controls=all_controls)
    
    if result_full and result_unadj:
        print(f"\n  FULLY ADJUSTED:")
        print(f"    Heat RR(P99): {result_full['rr_p99']['rr']:.4f} [{result_full['rr_p99']['rr_lower']:.4f}-{result_full['rr_p99']['rr_upper']:.4f}]")
        print(f"    Cold RR(P1):  {result_full['rr_p1']['rr']:.4f} [{result_full['rr_p1']['rr_lower']:.4f}-{result_full['rr_p1']['rr_upper']:.4f}]")
        
        heat_change = (result_full['rr_p99']['rr'] - result_unadj['rr_p99']['rr']) / result_unadj['rr_p99']['rr'] * 100
        cold_change = (result_full['rr_p1']['rr'] - result_unadj['rr_p1']['rr']) / result_unadj['rr_p1']['rr'] * 100
        
        print(f"\n  TOTAL CHANGE VS UNADJUSTED:")
        print(f"    Heat: {heat_change:+.1f}%")
        print(f"    Cold: {cold_change:+.1f}%")
else:
    result_full = None
    print("  No additional controls available")

# =============================================================================
# SAVE RESULTS
# =============================================================================

print("\n" + "="*70)
print("SAVING RESULTS")
print("="*70)

results = {
    'analysis_date': datetime.now().isoformat(),
    'configuration': {
        'max_lag': MAX_LAG,
        'poly_degree': POLY_DEGREE,
        'n_regions': int(len(valid_regions)),
    },
    'temperature_comparison': {
        'drybulb': result_drybulb,
        'apparent': result_apparent,
        'better_fit': 'apparent' if result_apparent and result_drybulb and result_apparent['aic'] < result_drybulb['aic'] else 'drybulb'
    },
    'pollution_adjustment': {
        'unadjusted': result_unadj,
        'adjusted': result_poll,
        'controls_used': poll_vars
    },
    'influenza_adjustment': {
        'unadjusted': result_unadj,
        'adjusted': result_flu,
        'controls_used': flu_vars
    },
    'fully_adjusted': {
        'result': result_full,
        'all_controls': all_controls
    }
}

output_file = os.path.join(OUTPUT_DIR, 'supplementary_analyses_regional.json')
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved: {output_file}")

# Create summary table
summary_rows = []

if result_drybulb:
    summary_rows.append({
        'model': 'Dry-bulb (baseline)',
        'heat_rr_p99': result_drybulb['rr_p99']['rr'],
        'heat_ci': f"[{result_drybulb['rr_p99']['rr_lower']:.3f}-{result_drybulb['rr_p99']['rr_upper']:.3f}]",
        'cold_rr_p1': result_drybulb['rr_p1']['rr'],
        'cold_ci': f"[{result_drybulb['rr_p1']['rr_lower']:.3f}-{result_drybulb['rr_p1']['rr_upper']:.3f}]",
        'aic': result_drybulb['aic']
    })

if result_apparent:
    summary_rows.append({
        'model': 'Apparent temperature',
        'heat_rr_p99': result_apparent['rr_p99']['rr'],
        'heat_ci': f"[{result_apparent['rr_p99']['rr_lower']:.3f}-{result_apparent['rr_p99']['rr_upper']:.3f}]",
        'cold_rr_p1': result_apparent['rr_p1']['rr'],
        'cold_ci': f"[{result_apparent['rr_p1']['rr_lower']:.3f}-{result_apparent['rr_p1']['rr_upper']:.3f}]",
        'aic': result_apparent['aic']
    })

if result_poll:
    summary_rows.append({
        'model': f'Pollution-adjusted ({", ".join(poll_vars)})',
        'heat_rr_p99': result_poll['rr_p99']['rr'],
        'heat_ci': f"[{result_poll['rr_p99']['rr_lower']:.3f}-{result_poll['rr_p99']['rr_upper']:.3f}]",
        'cold_rr_p1': result_poll['rr_p1']['rr'],
        'cold_ci': f"[{result_poll['rr_p1']['rr_lower']:.3f}-{result_poll['rr_p1']['rr_upper']:.3f}]",
        'aic': result_poll['aic']
    })

if result_flu:
    summary_rows.append({
        'model': f'Flu-adjusted ({flu_vars[0]})',
        'heat_rr_p99': result_flu['rr_p99']['rr'],
        'heat_ci': f"[{result_flu['rr_p99']['rr_lower']:.3f}-{result_flu['rr_p99']['rr_upper']:.3f}]",
        'cold_rr_p1': result_flu['rr_p1']['rr'],
        'cold_ci': f"[{result_flu['rr_p1']['rr_lower']:.3f}-{result_flu['rr_p1']['rr_upper']:.3f}]",
        'aic': result_flu['aic']
    })

if result_full:
    summary_rows.append({
        'model': 'Fully adjusted',
        'heat_rr_p99': result_full['rr_p99']['rr'],
        'heat_ci': f"[{result_full['rr_p99']['rr_lower']:.3f}-{result_full['rr_p99']['rr_upper']:.3f}]",
        'cold_rr_p1': result_full['rr_p1']['rr'],
        'cold_ci': f"[{result_full['rr_p1']['rr_lower']:.3f}-{result_full['rr_p1']['rr_upper']:.3f}]",
        'aic': result_full['aic']
    })

summary_df = pd.DataFrame(summary_rows)
summary_file = os.path.join(OUTPUT_DIR, 'supplementary_summary.csv')
summary_df.to_csv(summary_file, index=False)
print(f"Saved: {summary_file}")

# =============================================================================
# PRINT SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUPPLEMENTARY ANALYSES SUMMARY")
print("="*70)

print("\n" + summary_df.to_string(index=False))

print("\n" + "-"*70)
print("KEY FINDINGS:")
print("-"*70)

if result_drybulb and result_apparent:
    print(f"\n1. TEMPERATURE METRIC:")
    if result_apparent['aic'] < result_drybulb['aic']:
        print(f"   Apparent temperature provides better fit (ΔAIC = {result_drybulb['aic'] - result_apparent['aic']:.1f})")
    else:
        print(f"   Dry-bulb temperature provides adequate fit")

if result_unadj and result_poll:
    heat_change = (result_poll['rr_p99']['rr'] - result_unadj['rr_p99']['rr']) / result_unadj['rr_p99']['rr'] * 100
    print(f"\n2. POLLUTION ADJUSTMENT (supplementary):")
    print(f"   Heat effect change: {heat_change:+.1f}%")
    if abs(heat_change) < 5:
        print(f"   → Minimal confounding - unadjusted results valid for main analysis")
    else:
        print(f"   → Report in supplementary; may indicate mediation")

if result_unadj and result_flu:
    cold_change = (result_flu['rr_p1']['rr'] - result_unadj['rr_p1']['rr']) / result_unadj['rr_p1']['rr'] * 100
    print(f"\n3. INFLUENZA ADJUSTMENT:")
    print(f"   Cold effect change: {cold_change:+.1f}%")
    if cold_change < -10:
        print(f"   → Some cold effect may be flu-attributable")
    else:
        print(f"   → Cold effect robust to flu season adjustment")

print("\n" + "="*70)
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*70)
