"""
02b_harvesting_analysis.py
===========================
Harvesting (Mortality Displacement) Analysis for DLNM

Tests whether temperature-attributable deaths represent:
- TRUE EXCESS: Deaths that would not have occurred otherwise
- HARVESTING: Short-term displacement (deaths advanced by days/weeks)

Method (Armstrong et al., 2014; Gasparrini et al., 2015):
1. Fit DLNM with extended lag (up to 35 days)
2. Calculate CUMULATIVE effects at different lag horizons (7, 14, 21, 28, 35)
3. If cumulative effect decreases at longer lags → harvesting present
4. Harvesting ratio = 1 - (cumRR_35 / cumRR_7)

Input: Regional data from phase0_data_prep/results/
Output: Harvesting analysis results with lag-specific cumulative effects

References:
- Armstrong et al. (2014) Am J Epidemiol - Forward displacement
- Gasparrini et al. (2015) Lancet - Short-term displacement

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
print("02b: HARVESTING / MORTALITY DISPLACEMENT ANALYSIS")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PHASE0_RESULTS = os.path.join(BASE_DIR, 'new_analysis', 'phase0_data_prep', 'results')
OUTPUT_DIR = os.path.join(BASE_DIR, 'new_analysis', 'phase2_robustness', 'results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Extended lag to detect harvesting
MAX_LAG = 35  # Beyond typical 21-day window
LAG_HORIZONS = [7, 14, 21, 28, 35]  # Cumulative effect points
POLY_DEGREE = 3
MIN_OBS = 1000

print(f"\nConfiguration:")
print(f"  Extended max lag: {MAX_LAG} days")
print(f"  Lag horizons: {LAG_HORIZONS}")
print(f"  Polynomial degree: {POLY_DEGREE}")

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

# Temperature percentiles (including P2.5/P97.5 for consistency)
temp_p1 = df['temp_mean'].quantile(0.01)
temp_p2_5 = df['temp_mean'].quantile(0.025)
temp_p5 = df['temp_mean'].quantile(0.05)
temp_p50 = df['temp_mean'].quantile(0.50)
temp_p95 = df['temp_mean'].quantile(0.95)
temp_p97_5 = df['temp_mean'].quantile(0.975)
temp_p99 = df['temp_mean'].quantile(0.99)

print(f"\nTemperature: P1={temp_p1:.1f}°C, P2.5={temp_p2_5:.1f}°C, P50={temp_p50:.1f}°C, P97.5={temp_p97_5:.1f}°C, P99={temp_p99:.1f}°C")

# =============================================================================
# CREATE LAG MATRIX
# =============================================================================

print("\n" + "-"*70)
print("Creating Lag Matrix")
print("-"*70)

# Sort by region and date
df = df.sort_values(['region_code', 'date']).reset_index(drop=True)

# Create lagged temperatures
print(f"Creating lags up to {MAX_LAG} days...")
for lag in range(1, MAX_LAG + 1):
    df[f'temp_lag{lag}'] = df.groupby('region_code')['temp_mean'].shift(lag)
    if lag % 10 == 0:
        print(f"  Lag {lag} complete")

print(f"  All {MAX_LAG} lags created")

# =============================================================================
# DLNM FUNCTIONS (with expert-recommended fixes)
# =============================================================================

def fit_dlnm_with_lags(df_input, max_lag, poly_degree=3):
    """
    Fit DLNM with specified max lag and return lag-specific coefficients.
    
    Returns coefficients that allow calculating RR at any temperature
    for any lag horizon.
    """
    df = df_input.copy()
    
    # Get temperature stats
    temp_center = df['temp_mean'].median()
    temp_scale = df['temp_mean'].std()
    
    # Standardize current and lagged temperatures
    df['temp_std'] = (df['temp_mean'] - temp_center) / temp_scale
    for lag in range(1, max_lag + 1):
        df[f'temp_lag{lag}_std'] = (df[f'temp_lag{lag}'] - temp_center) / temp_scale
    
    # Build cross-basis: polynomial in temperature × lag weights
    cb_list = []
    cb_names = []
    
    for p in range(1, poly_degree + 1):
        for lag in range(max_lag + 1):
            if lag == 0:
                col = df['temp_std'].values ** p
            else:
                col = df[f'temp_lag{lag}_std'].values ** p
            cb_list.append(col)
            cb_names.append(f'temp_pow{p}_lag{lag}')
    
    # Add cross-basis to dataframe
    for i, name in enumerate(cb_names):
        df[name] = cb_list[i]
    
    # Remove rows with NaN
    valid_mask = ~df[cb_names].isna().any(axis=1)
    df_valid = df[valid_mask].copy()
    
    if len(df_valid) < 1000:
        return None
    
    # Build design matrix with controls
    X = pd.get_dummies(
        df_valid[['day_of_week', 'month', 'year', 'region_code']],
        columns=['day_of_week', 'month', 'year', 'region_code'],
        drop_first=True
    )
    for name in cb_names:
        X[name] = df_valid[name].astype(float)
    X['is_holiday'] = df_valid['is_holiday']
    X = sm.add_constant(X)
    
    y = df_valid['deaths_elderly'].values.astype(float)
    
    # Fit quasi-Poisson (scale='X2' estimates dispersion via Pearson chi-sq)
    try:
        model = GLM(y, X.astype(float), family=Poisson())
        result = model.fit(scale='X2', cov_type='HC1')
    except Exception as e:
        print(f"    Model fitting failed: {e}")
        return None
    
    # Extract coefficients using names (not position - more robust)
    # Expert fix: use result.params[name] instead of iloc[idx]
    coef_dict = {}
    for name in cb_names:
        if name in result.params.index:
            coef_dict[name] = float(result.params[name])
    
    # Get full covariance matrix for delta method CI calculation
    # Expert fix: store cov matrix to properly account for coefficient correlations
    cov_matrix = result.cov_params()
    
    return {
        'coefs': coef_dict,
        'cov': cov_matrix,  # Full covariance matrix for delta method
        'cb_names': cb_names,  # Cross-basis term names
        'temp_center': temp_center,
        'temp_scale': temp_scale,
        'n_obs': len(df_valid),
        'aic': result.aic,
        'dispersion': result.scale,  # Estimated dispersion parameter
        'poly_degree': poly_degree,
        'max_lag': max_lag
    }


def calc_cumulative_rr(model_results, target_temp, ref_temp, up_to_lag):
    """
    Calculate CUMULATIVE RR at target_temp vs ref_temp, summing over lags 0 to up_to_lag.
    
    Uses the DELTA METHOD with full covariance matrix for proper CI calculation.
    Expert fix: var(log RR) = d' Σ d, not sum of (se * diff)^2
    
    The cumulative effect captures the total mortality impact including delayed effects.
    """
    coefs = model_results['coefs']
    cov = model_results['cov']  # Full covariance matrix
    temp_center = model_results['temp_center']
    temp_scale = model_results['temp_scale']
    poly_degree = model_results['poly_degree']
    max_lag = model_results['max_lag']
    
    # Standardize temperatures
    target_std = (target_temp - temp_center) / temp_scale
    ref_std = (ref_temp - temp_center) / temp_scale
    
    # Build coefficient names and difference vector for delta method
    names = []
    diffs = []
    betas = []
    
    for lag in range(min(up_to_lag + 1, max_lag + 1)):
        for p in range(1, poly_degree + 1):
            coef_name = f'temp_pow{p}_lag{lag}'
            if coef_name in coefs:
                names.append(coef_name)
                diff = (target_std ** p) - (ref_std ** p)
                diffs.append(diff)
                betas.append(coefs[coef_name])
    
    # Convert to arrays
    d = np.array(diffs)  # Gradient vector
    beta_vec = np.array(betas)
    
    # Calculate log(RR) = beta' * d
    log_rr = float(np.dot(beta_vec, d))
    
    # DELTA METHOD: var(log RR) = d' Σ d
    # Extract sub-covariance matrix for the relevant coefficients
    try:
        if isinstance(cov, pd.DataFrame):
            # Reindex to get only our coefficient subset, fill missing with 0
            cov_sub = cov.reindex(index=names, columns=names).fillna(0).to_numpy()
        else:
            # If ndarray, assume order matches
            cov_sub = cov
        
        var_log_rr = float(d @ cov_sub @ d)
        se_log_rr = np.sqrt(max(var_log_rr, 0))  # Guard against numerical issues
    except Exception as e:
        # Fallback to independence assumption if covariance extraction fails
        print(f"    Warning: Covariance extraction failed ({e}), using independence assumption")
        var_log_rr = sum((cov.loc[n, n] if n in cov.index else 0) * (d[i]**2) 
                         for i, n in enumerate(names))
        se_log_rr = np.sqrt(max(var_log_rr, 0))
    
    rr = np.exp(log_rr)
    rr_lower = np.exp(log_rr - 1.96 * se_log_rr)
    rr_upper = np.exp(log_rr + 1.96 * se_log_rr)
    
    return rr, rr_lower, rr_upper

# =============================================================================
# FIT EXTENDED LAG MODEL
# =============================================================================

print("\n" + "="*70)
print(f"FITTING DLNM WITH {MAX_LAG}-DAY LAG")
print("="*70)

model_results = fit_dlnm_with_lags(df, max_lag=MAX_LAG, poly_degree=POLY_DEGREE)

if model_results is None:
    print("ERROR: Model fitting failed!")
    exit(1)

print(f"\n  Model fitted successfully")
print(f"  N observations: {model_results['n_obs']:,}")
print(f"  AIC: {model_results['aic']:.1f}")

# =============================================================================
# CALCULATE CUMULATIVE EFFECTS AT DIFFERENT LAG HORIZONS
# =============================================================================

print("\n" + "="*70)
print("CUMULATIVE EFFECTS AT DIFFERENT LAG HORIZONS")
print("="*70)

# Heat effects (P97.5 as primary, P99 as comparison)
print("\n--- HEAT (P97.5 vs P50) - PRIMARY ---")
heat_cumulative_97_5 = []
for horizon in LAG_HORIZONS:
    rr, rr_lo, rr_hi = calc_cumulative_rr(model_results, temp_p97_5, temp_p50, horizon)
    heat_cumulative_97_5.append({
        'lag_horizon': horizon,
        'cumulative_rr': rr,
        'rr_lower': rr_lo,
        'rr_upper': rr_hi
    })
    print(f"  Lag 0-{horizon:2d}: RR = {rr:.4f} [{rr_lo:.4f}-{rr_hi:.4f}]")

print("\n--- HEAT (P99 vs P50) - COMPARISON ---")
heat_cumulative = []
for horizon in LAG_HORIZONS:
    rr, rr_lo, rr_hi = calc_cumulative_rr(model_results, temp_p99, temp_p50, horizon)
    heat_cumulative.append({
        'lag_horizon': horizon,
        'cumulative_rr': rr,
        'rr_lower': rr_lo,
        'rr_upper': rr_hi
    })
    print(f"  Lag 0-{horizon:2d}: RR = {rr:.4f} [{rr_lo:.4f}-{rr_hi:.4f}]")

# Cold effects (P2.5 as primary, P1 as comparison)
print("\n--- COLD (P2.5 vs P50) - PRIMARY ---")
cold_cumulative_2_5 = []
for horizon in LAG_HORIZONS:
    rr, rr_lo, rr_hi = calc_cumulative_rr(model_results, temp_p2_5, temp_p50, horizon)
    cold_cumulative_2_5.append({
        'lag_horizon': horizon,
        'cumulative_rr': rr,
        'rr_lower': rr_lo,
        'rr_upper': rr_hi
    })
    print(f"  Lag 0-{horizon:2d}: RR = {rr:.4f} [{rr_lo:.4f}-{rr_hi:.4f}]")

print("\n--- COLD (P1 vs P50) - COMPARISON ---")
cold_cumulative = []
for horizon in LAG_HORIZONS:
    rr, rr_lo, rr_hi = calc_cumulative_rr(model_results, temp_p1, temp_p50, horizon)
    cold_cumulative.append({
        'lag_horizon': horizon,
        'cumulative_rr': rr,
        'rr_lower': rr_lo,
        'rr_upper': rr_hi
    })
    print(f"  Lag 0-{horizon:2d}: RR = {rr:.4f} [{rr_lo:.4f}-{rr_hi:.4f}]")

# =============================================================================
# CALCULATE HARVESTING RATIO (with expert-recommended guards)
# =============================================================================

def calc_harvesting_ratio(err_short, err_long, min_err_threshold=0.01):
    """
    Calculate harvesting ratio with proper guards.
    
    Expert fix: Guard against division by zero and handle edge cases.
    
    Parameters:
    -----------
    err_short : float - Excess relative risk at short lag (e.g., 7 days)
    err_long : float - Excess relative risk at long lag (e.g., 35 days)
    min_err_threshold : float - Minimum ERR to consider meaningful (default 0.01 = 1%)
    
    Returns:
    --------
    tuple: (harvesting_ratio, interpretation, is_valid)
    """
    # Guard 1: Check for NaN
    if np.isnan(err_short) or np.isnan(err_long):
        return np.nan, "Invalid (NaN values)", False
    
    # Guard 2: Check if short-term ERR is too small (near zero)
    if abs(err_short) < min_err_threshold:
        return np.nan, "Invalid (short-term ERR near zero)", False
    
    # Guard 3: Check for protective effect (ERR < 0)
    if err_short < 0:
        return np.nan, "Invalid (protective effect at short lag)", False
    
    # Calculate ratio
    ratio = 1 - (err_long / err_short)
    
    # Interpret
    if ratio > 1.0:
        interp = "Overcorrection (ERR_long negative) - model instability"
        valid = False
    elif ratio > 0.5:
        interp = "High harvesting (>50% displaced)"
        valid = True
    elif ratio > 0.2:
        interp = "Moderate harvesting (20-50% displaced)"
        valid = True
    elif ratio >= 0:
        interp = "Low harvesting (<20%) - mostly true excess"
        valid = True
    else:  # ratio < 0
        interp = "Negative harvesting - effects INCREASE over time (persistent burden)"
        valid = True
    
    return ratio, interp, valid

print("\n" + "="*70)
print("HARVESTING ANALYSIS")
print("="*70)

# Extract RRs and calculate ERRs for P99/P1 (comparison thresholds)
heat_rr_7 = heat_cumulative[0]['cumulative_rr']
heat_rr_35 = heat_cumulative[-1]['cumulative_rr']
cold_rr_7 = cold_cumulative[0]['cumulative_rr']
cold_rr_35 = cold_cumulative[-1]['cumulative_rr']

heat_err_7 = heat_rr_7 - 1
heat_err_35 = heat_rr_35 - 1
cold_err_7 = cold_rr_7 - 1
cold_err_35 = cold_rr_35 - 1

# Calculate harvesting ratios with guards
heat_harvesting_ratio, heat_interp, heat_valid = calc_harvesting_ratio(heat_err_7, heat_err_35)
cold_harvesting_ratio, cold_interp, cold_valid = calc_harvesting_ratio(cold_err_7, cold_err_35)

# Also calculate for P97.5/P2.5 (PRIMARY thresholds)
heat_err_7_97_5 = heat_cumulative_97_5[0]['cumulative_rr'] - 1
heat_err_35_97_5 = heat_cumulative_97_5[-1]['cumulative_rr'] - 1
cold_err_7_2_5 = cold_cumulative_2_5[0]['cumulative_rr'] - 1
cold_err_35_2_5 = cold_cumulative_2_5[-1]['cumulative_rr'] - 1

heat_harvesting_ratio_97_5, heat_interp_97_5, heat_valid_97_5 = calc_harvesting_ratio(heat_err_7_97_5, heat_err_35_97_5)
cold_harvesting_ratio_2_5, cold_interp_2_5, cold_valid_2_5 = calc_harvesting_ratio(cold_err_7_2_5, cold_err_35_2_5)

# Print detailed results for PRIMARY thresholds
print(f"\n--- HEAT HARVESTING (PRIMARY P97.5) ---")
print(f"  RR at lag 0-7:   {heat_cumulative_97_5[0]['cumulative_rr']:.4f} [{heat_cumulative_97_5[0]['rr_lower']:.4f}-{heat_cumulative_97_5[0]['rr_upper']:.4f}]")
print(f"  RR at lag 0-35:  {heat_cumulative_97_5[-1]['cumulative_rr']:.4f} [{heat_cumulative_97_5[-1]['rr_lower']:.4f}-{heat_cumulative_97_5[-1]['rr_upper']:.4f}]")
print(f"  ERR at lag 0-7:  {heat_err_7_97_5:.4f}")
print(f"  ERR at lag 0-35: {heat_err_35_97_5:.4f}")
if heat_valid_97_5:
    print(f"  Harvesting ratio: {heat_harvesting_ratio_97_5:.1%}")
else:
    print(f"  Harvesting ratio: N/A")
print(f"  Interpretation: {heat_interp_97_5}")

print(f"\n--- HEAT HARVESTING (COMPARISON P99) ---")
print(f"  RR at lag 0-7:   {heat_rr_7:.4f} [{heat_cumulative[0]['rr_lower']:.4f}-{heat_cumulative[0]['rr_upper']:.4f}]")
print(f"  RR at lag 0-35:  {heat_rr_35:.4f} [{heat_cumulative[-1]['rr_lower']:.4f}-{heat_cumulative[-1]['rr_upper']:.4f}]")
print(f"  ERR at lag 0-7:  {heat_err_7:.4f}")
print(f"  ERR at lag 0-35: {heat_err_35:.4f}")
if heat_valid:
    print(f"  Harvesting ratio: {heat_harvesting_ratio:.1%}")
else:
    print(f"  Harvesting ratio: N/A")
print(f"  Interpretation: {heat_interp}")

print(f"\n--- COLD HARVESTING (PRIMARY P2.5) ---")
print(f"  RR at lag 0-7:   {cold_cumulative_2_5[0]['cumulative_rr']:.4f} [{cold_cumulative_2_5[0]['rr_lower']:.4f}-{cold_cumulative_2_5[0]['rr_upper']:.4f}]")
print(f"  RR at lag 0-35:  {cold_cumulative_2_5[-1]['cumulative_rr']:.4f} [{cold_cumulative_2_5[-1]['rr_lower']:.4f}-{cold_cumulative_2_5[-1]['rr_upper']:.4f}]")
print(f"  ERR at lag 0-7:  {cold_err_7_2_5:.4f}")
print(f"  ERR at lag 0-35: {cold_err_35_2_5:.4f}")
if cold_valid_2_5:
    print(f"  Harvesting ratio: {cold_harvesting_ratio_2_5:.1%}")
else:
    print(f"  Harvesting ratio: N/A")
print(f"  Interpretation: {cold_interp_2_5}")

print(f"\n--- COLD HARVESTING (COMPARISON P1) ---")
print(f"  RR at lag 0-7:   {cold_rr_7:.4f} [{cold_cumulative[0]['rr_lower']:.4f}-{cold_cumulative[0]['rr_upper']:.4f}]")
print(f"  RR at lag 0-35:  {cold_rr_35:.4f} [{cold_cumulative[-1]['rr_lower']:.4f}-{cold_cumulative[-1]['rr_upper']:.4f}]")
print(f"  ERR at lag 0-7:  {cold_err_7:.4f}")
print(f"  ERR at lag 0-35: {cold_err_35:.4f}")
if cold_valid:
    print(f"  Harvesting ratio: {cold_harvesting_ratio:.1%}")
else:
    print(f"  Harvesting ratio: N/A")
print(f"  Interpretation: {cold_interp}")

# =============================================================================
# SAVE RESULTS
# =============================================================================

print("\n" + "="*70)
print("SAVING RESULTS")
print("="*70)

results = {
    'analysis_date': datetime.now().isoformat(),
    'methodology_notes': {
        'ci_method': 'Delta method with full covariance matrix',
        'harvesting_formula': '1 - (ERR_35day / ERR_7day)',
        'expert_review': 'Fixes applied per Dec 2025 review: covariance-aware CIs, robust coefficient extraction, harvesting guards'
    },
    'configuration': {
        'max_lag': MAX_LAG,
        'poly_degree': POLY_DEGREE,
        'lag_horizons': LAG_HORIZONS,
        'n_regions': int(len(valid_regions)),
        'n_observations': int(model_results['n_obs']),
        'dispersion': float(model_results.get('dispersion', np.nan)),
    },
    'temperature_percentiles': {
        'p1': float(temp_p1),
        'p2_5': float(temp_p2_5),
        'p5': float(temp_p5),
        'p50': float(temp_p50),
        'p95': float(temp_p95),
        'p97_5': float(temp_p97_5),
        'p99': float(temp_p99),
    },
    'heat_cumulative_p97_5': heat_cumulative_97_5,
    'heat_cumulative_p99': heat_cumulative,
    'cold_cumulative_p2_5': cold_cumulative_2_5,
    'cold_cumulative_p1': cold_cumulative,
    'harvesting': {
        'heat_p99': {
            'rr_7day': float(heat_rr_7),
            'rr_35day': float(heat_rr_35),
            'err_short_7day': float(heat_err_7),
            'err_long_35day': float(heat_err_35),
            'harvesting_ratio': float(heat_harvesting_ratio) if heat_valid else None,
            'percent_displaced': f"{heat_harvesting_ratio:.1%}" if heat_valid else None,
            'is_valid': heat_valid,
            'interpretation': heat_interp
        },
        'heat_p97_5': {
            'rr_7day': float(heat_cumulative_97_5[0]['cumulative_rr']),
            'rr_35day': float(heat_cumulative_97_5[-1]['cumulative_rr']),
            'err_short_7day': float(heat_err_7_97_5),
            'err_long_35day': float(heat_err_35_97_5),
            'harvesting_ratio': float(heat_harvesting_ratio_97_5) if heat_valid_97_5 else None,
            'percent_displaced': f"{heat_harvesting_ratio_97_5:.1%}" if heat_valid_97_5 else None,
            'is_valid': heat_valid_97_5,
            'interpretation': heat_interp_97_5
        },
        'cold_p1': {
            'rr_7day': float(cold_rr_7),
            'rr_35day': float(cold_rr_35),
            'err_short_7day': float(cold_err_7),
            'err_long_35day': float(cold_err_35),
            'harvesting_ratio': float(cold_harvesting_ratio) if cold_valid else None,
            'percent_displaced': f"{cold_harvesting_ratio:.1%}" if cold_valid else None,
            'is_valid': cold_valid,
            'interpretation': cold_interp
        },
        'cold_p2_5': {
            'rr_7day': float(cold_cumulative_2_5[0]['cumulative_rr']),
            'rr_35day': float(cold_cumulative_2_5[-1]['cumulative_rr']),
            'err_short_7day': float(cold_err_7_2_5),
            'err_long_35day': float(cold_err_35_2_5),
            'harvesting_ratio': float(cold_harvesting_ratio_2_5) if cold_valid_2_5 else None,
            'percent_displaced': f"{cold_harvesting_ratio_2_5:.1%}" if cold_valid_2_5 else None,
            'is_valid': cold_valid_2_5,
            'interpretation': cold_interp_2_5
        }
    },
    'interpretation': {
        'heat_p99': heat_interp,
        'heat_p97_5': heat_interp_97_5,
        'cold_p1': cold_interp,
        'cold_p2_5': cold_interp_2_5
    }
}

# Save results
output_file = os.path.join(OUTPUT_DIR, 'harvesting_analysis_regional.json')
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved: {output_file}")

# Create summary table
summary_data = []
for h in heat_cumulative:
    summary_data.append({
        'type': 'heat',
        'temperature': f'P99 ({temp_p99:.1f}°C) vs P50',
        'lag_horizon': h['lag_horizon'],
        'cumulative_rr': h['cumulative_rr'],
        'ci_lower': h['rr_lower'],
        'ci_upper': h['rr_upper']
    })
for c in cold_cumulative:
    summary_data.append({
        'type': 'cold',
        'temperature': f'P1 ({temp_p1:.1f}°C) vs P50',
        'lag_horizon': c['lag_horizon'],
        'cumulative_rr': c['cumulative_rr'],
        'ci_lower': c['rr_lower'],
        'ci_upper': c['rr_upper']
    })

summary_df = pd.DataFrame(summary_data)
summary_file = os.path.join(OUTPUT_DIR, 'harvesting_summary.csv')
summary_df.to_csv(summary_file, index=False)
print(f"Saved: {summary_file}")

# =============================================================================
# PRINT SUMMARY
# =============================================================================

print("\n" + "="*70)
print("HARVESTING ANALYSIS SUMMARY")
print("="*70)

print("\n" + summary_df.to_string(index=False))

print("\n" + "-"*70)
print("KEY FINDINGS:")
print("-"*70)

print(f"\n1. HEAT (PRIMARY P97.5): {heat_interp_97_5}")
if heat_valid_97_5:
    print(f"   - Harvesting ratio: {heat_harvesting_ratio_97_5:.1%}")
    print(f"   - True excess heat deaths ≈ {(1-heat_harvesting_ratio_97_5)*100:.0f}% of unadjusted estimate")

print(f"\n2. HEAT (COMPARISON P99): {heat_interp}")
if heat_valid:
    print(f"   - Harvesting ratio: {heat_harvesting_ratio:.1%}")

print(f"\n3. COLD (PRIMARY P2.5): {cold_interp_2_5}")
if cold_valid_2_5:
    print(f"   - Harvesting ratio: {cold_harvesting_ratio_2_5:.1%}")

print(f"\n4. COLD (COMPARISON P1): {cold_interp}")
if cold_valid:
    print(f"   - Harvesting ratio: {cold_harvesting_ratio:.1%}")

print("\n5. IMPLICATIONS FOR BURDEN ESTIMATES:")
# Use PRIMARY threshold (P97.5) for heat adjustment
if heat_valid_97_5 and heat_harvesting_ratio_97_5 > 0.5:
    print("   - Heat burden estimates should be adjusted downward for harvesting")
    print(f"   - Using PRIMARY (P97.5) ratio: {heat_harvesting_ratio_97_5:.1%} displaced")
elif heat_valid_97_5:
    print(f"   - Heat harvesting ratio: {heat_harvesting_ratio_97_5:.1%}")
else:
    print("   - Heat harvesting could not be reliably estimated")

# Cold typically shows persistent effects
if cold_valid_2_5 and cold_harvesting_ratio_2_5 < 0:
    print("   - Cold effects PERSISTENT (negative harvesting) - no adjustment needed")
elif cold_valid_2_5:
    print(f"   - Cold harvesting ratio: {cold_harvesting_ratio_2_5:.1%}")

# =============================================================================
# HARVESTING-ADJUSTED BURDEN ESTIMATES
# =============================================================================

print("\n" + "="*70)
print("HARVESTING-ADJUSTED ATTRIBUTABLE BURDEN")
print("="*70)

# Load Phase 1 burden estimates
try:
    phase1_burden_file = os.path.join(BASE_DIR, 'new_analysis', 'phase1_core_model', 'results', 'attributable_burden_national.json')
    with open(phase1_burden_file, 'r') as f:
        phase1_burden = json.load(f)
    
    for level in ['intermediate', 'immediate']:
        level_data = phase1_burden.get(level, {})
        if not level_data:
            continue
            
        print(f"\n{level.upper()} LEVEL ({level_data.get('n_regions', 'N/A')} regions):")
        print("-" * 50)
        
        # P2.5/P97.5 (PRIMARY)
        heat_an_97_5 = level_data.get('heat_an_97_5', 0)
        cold_an_2_5 = level_data.get('cold_an_2_5', 0)
        
        # P1/P99 (COMPARISON)
        heat_an_99 = level_data.get('heat_an_99', 0)
        cold_an_1 = level_data.get('cold_an_1', 0)
        
        # Calculate harvesting-adjusted burden using PRIMARY threshold ratios
        # Heat: reduce by harvesting ratio (displaced deaths)
        # Use P97.5 harvesting ratio for P97.5 burden, P99 for P99 burden
        if heat_valid_97_5 and not np.isnan(heat_harvesting_ratio_97_5):
            heat_adj_factor_97_5 = 1 - max(0, min(heat_harvesting_ratio_97_5, 1))  # Clip to [0,1]
        else:
            heat_adj_factor_97_5 = 1.0  # No adjustment if ratio invalid
            
        if heat_valid and not np.isnan(heat_harvesting_ratio):
            heat_adj_factor_99 = 1 - max(0, min(heat_harvesting_ratio, 1))
        else:
            heat_adj_factor_99 = 1.0
        
        # For cold: if negative harvesting (effects persistent), no reduction needed
        if cold_valid_2_5 and not np.isnan(cold_harvesting_ratio_2_5):
            if cold_harvesting_ratio_2_5 < 0:
                cold_adj_factor_2_5 = 1.0  # Persistent effects - no reduction
            else:
                cold_adj_factor_2_5 = 1 - max(0, min(cold_harvesting_ratio_2_5, 1))
        else:
            cold_adj_factor_2_5 = 1.0
            
        if cold_valid and not np.isnan(cold_harvesting_ratio):
            if cold_harvesting_ratio < 0:
                cold_adj_factor_1 = 1.0
            else:
                cold_adj_factor_1 = 1 - max(0, min(cold_harvesting_ratio, 1))
        else:
            cold_adj_factor_1 = 1.0
        
        heat_adj_97_5 = heat_an_97_5 * heat_adj_factor_97_5
        heat_adj_99 = heat_an_99 * heat_adj_factor_99
        cold_adj_2_5 = cold_an_2_5 * cold_adj_factor_2_5
        cold_adj_1 = cold_an_1 * cold_adj_factor_1
        
        print(f"\nP2.5/P97.5 (PRIMARY):")
        print(f"  Heat (>P97.5):")
        print(f"    Unadjusted: {heat_an_97_5:,.0f} deaths")
        print(f"    Harvesting-adjusted: {heat_adj_97_5:,.0f} deaths ({heat_adj_factor_97_5:.1%} retained)")
        print(f"  Cold (<P2.5):")
        print(f"    Unadjusted: {cold_an_2_5:,.0f} deaths")
        if cold_valid_2_5 and cold_harvesting_ratio_2_5 < 0:
            print(f"    Harvesting-adjusted: {cold_adj_2_5:,.0f} deaths (persistent effects - no reduction)")
        else:
            print(f"    Harvesting-adjusted: {cold_adj_2_5:,.0f} deaths ({cold_adj_factor_2_5:.1%} retained)")
        
        print(f"\nP1/P99 (COMPARISON):")
        print(f"  Heat (>P99):")
        print(f"    Unadjusted: {heat_an_99:,.0f} deaths")
        print(f"    Harvesting-adjusted: {heat_adj_99:,.0f} deaths ({heat_adj_factor_99:.1%} retained)")
        print(f"  Cold (<P1):")
        print(f"    Unadjusted: {cold_an_1:,.0f} deaths")
        if cold_valid and cold_harvesting_ratio < 0:
            print(f"    Harvesting-adjusted: {cold_adj_1:,.0f} deaths (persistent effects - no reduction)")
        else:
            print(f"    Harvesting-adjusted: {cold_adj_1:,.0f} deaths ({cold_adj_factor_1:.1%} retained)")
        
        # Update results with adjusted burden
        results['adjusted_burden'] = results.get('adjusted_burden', {})
        results['adjusted_burden'][level] = {
            'heat_adj_factor_97_5': float(heat_adj_factor_97_5),
            'heat_adj_factor_99': float(heat_adj_factor_99),
            'cold_adj_factor_2_5': float(cold_adj_factor_2_5),
            'cold_adj_factor_1': float(cold_adj_factor_1),
            'heat_an_97_5_unadj': float(heat_an_97_5),
            'heat_an_97_5_adj': float(heat_adj_97_5),
            'cold_an_2_5_unadj': float(cold_an_2_5),
            'cold_an_2_5_adj': float(cold_adj_2_5),
            'heat_an_99_unadj': float(heat_an_99),
            'heat_an_99_adj': float(heat_adj_99),
            'cold_an_1_unadj': float(cold_an_1),
            'cold_an_1_adj': float(cold_adj_1),
            'heat_annual_97_5_adj': float(heat_adj_97_5 / 14),  # 14 years
            'cold_annual_2_5_adj': float(cold_adj_2_5 / 14),
            'heat_annual_99_adj': float(heat_adj_99 / 14),
            'cold_annual_1_adj': float(cold_adj_1 / 14),
        }
    
    # Re-save results with adjusted burden
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nUpdated: {output_file}")
    
except FileNotFoundError:
    print("\n  Phase 1 burden file not found - run 01d_attributable_burden.py first")
except Exception as e:
    print(f"\n  Error loading Phase 1 burden: {e}")

print("\n" + "="*70)
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*70)
