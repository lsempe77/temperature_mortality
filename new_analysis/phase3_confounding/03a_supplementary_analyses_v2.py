"""
03a_supplementary_analyses_v2.py
=================================
Supplementary Analyses for DLNM Results (v2 - Corrected Implementation)

Addresses issues in v1:
1. Uses natural spline cross-basis (not polynomial) from utils/dlnm_module
2. Fits REGION-SPECIFIC models then meta-analyzes (not pooled GLM)
3. Includes population offset in all models
4. Uses delta method with full covariance for CIs
5. Handles missing pollution/flu properly (interpolation, not zeros)
6. Uses region-specific percentiles
7. Validates apparent temperature formula

These are NOT primary results but SUPPLEMENTARY material showing:
1. Apparent temperature (humidity-adjusted) vs dry-bulb
2. Pollution-adjusted model (PM2.5, O3 controls)
3. Influenza-adjusted model (flu season control)
4. Combined adjusted model

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
from patsy import dmatrix
from scipy import stats
import os
import sys
import json
from datetime import datetime
from pathlib import Path
import warnings
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
warnings.filterwarnings('ignore')

# Import shared DLNM utilities
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.dlnm_module import (
    ns_basis,
    create_crossbasis_ns,
    create_lag_matrix,
    compute_cumulative_rr_ns_with_se,
    find_mmt_from_coefficients,
    compute_effects_relative_to_mmt,
    create_time_spline,
    meta_random_effects,
    random_effects_meta_analysis,
    pool_region_results,
    convert_to_json_serializable,
    # MVMeta for proper coefficient pooling (Gasparrini methodology)
    mvmeta_pool_regions,
)

print("="*70)
print("PHASE 3: SUPPLEMENTARY ANALYSES v2")
print("       (Region-Specific Models + Meta-Analysis)")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = Path(__file__).parent.parent.parent
PHASE0_RESULTS = BASE_DIR / 'new_analysis' / 'phase0_data_prep' / 'results'
OUTPUT_DIR = Path(__file__).parent / 'results'
OUTPUT_DIR.mkdir(exist_ok=True)

# DLNM parameters (aligned with Phase 1/2 v2 core specs)
MAX_LAG = 21
TEMP_DF = 4                   # Natural spline df for temperature
LAG_DF = 4                    # Natural spline df for lag
TIME_SPLINE_DF_PER_YEAR = 1   # df per year for long-term trend (coarse; month dummies handle seasonality)
MIN_OBS = 1000                # Minimum observations per region
MIN_POLLUTION_COVERAGE = 0.7  # Require 70% non-missing pollution for adjusted models

# Temperature knots (Gasparrini methodology)
TEMP_KNOT_PCTS = [10, 75, 90]
TEMP_BOUNDARY_PCTS = [1, 99]
LAG_KNOTS = [1, 3, 7, 14]

print(f"\nConfiguration:")
print(f"  Max lag: {MAX_LAG} days")
print(f"  Temperature df: {TEMP_DF} (natural splines)")
print(f"  Lag df: {LAG_DF} (natural splines)")
print(f"  Minimum observations: {MIN_OBS}")
print(f"  Minimum pollution coverage: {MIN_POLLUTION_COVERAGE:.0%}")


# =============================================================================
# DATA LOADING
# =============================================================================

print("\n" + "-"*70)
print("Loading Data")
print("-"*70)

parser = argparse.ArgumentParser(description='Supplementary DLNM analyses (apparent temp, pollution, flu, full-adjusted)')
parser.add_argument('--level', type=str, default='intermediate',
                    choices=['intermediate', 'immediate'],
                    help='Spatial level to analyze')
args, _ = parser.parse_known_args()

level = args.level
if level == 'intermediate':
    temp_filename = 'era5_intermediate_daily.parquet'
    poll_filename = 'cams_intermediate_daily.parquet'
    flu_filename = 'influenza_daily_by_intermediate_region.parquet'
    mort_filename = 'mortality_regional_daily_elderly.parquet'
    ses_filename = 'ses_intermediate_covariates.csv'
    ses_region_col = 'intermediate_code'
    level_label = 'INTERMEDIATE (133 regions)'
else:
    temp_filename = 'era5_immediate_daily.parquet'
    poll_filename = 'cams_immediate_daily.parquet'
    flu_filename = 'influenza_daily_by_immediate_region.parquet'
    mort_filename = 'mortality_immediate_daily_elderly.parquet'
    ses_filename = 'ses_immediate_covariates.csv'
    ses_region_col = 'immediate_code'
    level_label = 'IMMEDIATE (510 regions)'

print(f"Level: {level} [{level_label}]")

# Load temperature data (includes dewpoint for apparent temp)
df_temp = pd.read_parquet(PHASE0_RESULTS / temp_filename)
df_temp['date'] = pd.to_datetime(df_temp['date'])
print(f"ERA5: {len(df_temp):,} rows, {df_temp['region_code'].nunique()} regions")

# Load pollution data
df_poll = pd.read_parquet(PHASE0_RESULTS / poll_filename)
df_poll['date'] = pd.to_datetime(df_poll['date'])
print(f"CAMS Pollution: {len(df_poll):,} rows")

# Load influenza data
df_flu = pd.read_parquet(PHASE0_RESULTS / flu_filename)
df_flu['date'] = pd.to_datetime(df_flu['date'])
print(f"Influenza: {len(df_flu):,} rows")

# Load mortality data
df_mort = pd.read_parquet(PHASE0_RESULTS / mort_filename)
df_mort['date'] = pd.to_datetime(df_mort['date'])
print(f"Mortality: {len(df_mort):,} rows")

# Load SES data for population
df_ses = pd.read_csv(PHASE0_RESULTS / ses_filename)
pop_map = dict(zip(df_ses[ses_region_col], df_ses['pop_elderly']))

# Load holidays
df_holidays = pd.read_parquet(PHASE0_RESULTS / 'brazilian_holidays_daily.parquet')
df_holidays['date'] = pd.to_datetime(df_holidays['date'])
print(f"Holidays: {len(df_holidays)} days")


# =============================================================================
# DATA MERGING AND PREPARATION
# =============================================================================

print("\n" + "-"*70)
print("Merging and Preparing Data")
print("-"*70)

# Standardize region column names for merging
if level == 'immediate':
    df_mort = df_mort.rename(columns={'immediate_code': 'region_code'})
elif level == 'intermediate':
    if 'intermediate_code' in df_mort.columns:
        df_mort = df_mort.rename(columns={'intermediate_code': 'region_code'})

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

# Standardize pollution file region column
if level == 'immediate':
    if 'immediate_code' in df_poll.columns:
        df_poll = df_poll.rename(columns={'immediate_code': 'region_code'})
elif level == 'intermediate':
    if 'intermediate_code' in df_poll.columns:
        df_poll = df_poll.rename(columns={'intermediate_code': 'region_code'})

# Standardize flu file region column
if level == 'immediate':
    if 'immediate_code' in df_flu.columns:
        df_flu = df_flu.rename(columns={'immediate_code': 'region_code'})
elif level == 'intermediate':
    if 'intermediate_code' in df_flu.columns:
        df_flu = df_flu.rename(columns={'intermediate_code': 'region_code'})

# Merge pollution (handle column naming)
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
flu_cols = [c for c in df_flu.columns if c not in ['date', 'region_code']]
df = pd.merge(
    df,
    df_flu[['date', 'region_code'] + flu_cols],
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
df['time_index'] = (df['date'] - df['date'].min()).dt.days


# =============================================================================
# CALCULATE APPARENT TEMPERATURE (Validated Formula)
# =============================================================================

print("\n" + "-"*70)
print("Calculating Apparent Temperature")
print("-"*70)

if 'dewpoint_mean' in df.columns:
    Td = df['dewpoint_mean'].values
    Ta = df['temp_mean'].values
    
    # Calculate relative humidity from dewpoint and air temperature
    # Using Magnus-Tetens approximation
    # es(T) = 6.112 × exp((17.67 × T) / (T + 243.5))
    es_Td = 6.112 * np.exp((17.67 * Td) / (Td + 243.5))  # Saturation vapor pressure at dewpoint
    es_Ta = 6.112 * np.exp((17.67 * Ta) / (Ta + 243.5))  # Saturation vapor pressure at air temp
    
    # Relative humidity (%)
    RH = 100.0 * (es_Td / es_Ta)
    RH = np.clip(RH, 0, 100)  # Ensure physical bounds
    
    # Apparent Temperature (Australian Bureau of Meteorology formula)
    # AT = Ta + 0.33 × e − 0.70 × ws − 4.00
    # Where e = water vapor pressure (hPa), ws = wind speed (m/s)
    # Without wind speed, simplified: AT = Ta + 0.33 × e − 4.00
    # e = (RH/100) × es_Ta
    e = (RH / 100.0) * es_Ta
    df['apparent_temp'] = Ta + 0.33 * e - 4.0
    
    # Validation: check distribution
    at_diff = df['apparent_temp'] - df['temp_mean']
    print(f"  Apparent temp calculated")
    print(f"  AT - Ta difference: mean={at_diff.mean():.1f}°C, "
          f"sd={at_diff.std():.1f}°C, range=[{at_diff.min():.1f}, {at_diff.max():.1f}]")
    
    if at_diff.mean() < -10 or at_diff.mean() > 20:
        print(f"  WARNING: AT-Ta difference may indicate formula issues")
else:
    df['apparent_temp'] = df['temp_mean']
    print(f"  WARNING: No dewpoint data, using dry-bulb as apparent temp")


# =============================================================================
# HANDLE MISSING POLLUTION & INFLUENZA (Interpolation, not zeros)
# =============================================================================

print("\n" + "-"*70)
print("Handling Missing Confounders")
print("-"*70)

# Sort by region and date for proper interpolation
df = df.sort_values(['region_code', 'date']).reset_index(drop=True)

# Pollution: interpolate within region (max 14 day gap)
poll_vars = ['pm25', 'pm10', 'o3']
for col in poll_vars:
    if col in df.columns:
        n_missing_before = df[col].isna().sum()
        # Interpolate within region, limit to 14 days
        df[col] = df.groupby('region_code')[col].transform(
            lambda x: x.interpolate(method='linear', limit=14)
        )
        # Forward/backward fill for remaining gaps at edges
        df[col] = df.groupby('region_code')[col].transform(
            lambda x: x.fillna(method='ffill').fillna(method='bfill')
        )
        n_missing_after = df[col].isna().sum()
        coverage = 1 - n_missing_after / len(df)
        print(f"  {col}: {n_missing_before:,} → {n_missing_after:,} missing "
              f"(coverage: {coverage:.1%})")
        
        # Create missingness indicator for sensitivity
        df[f'{col}_was_missing'] = df[col].isna().astype(int)

# Influenza: interpolate within region
flu_vars = ['srag_cases', 'confirmed_flu']
for col in flu_vars:
    if col in df.columns:
        n_missing_before = df[col].isna().sum()
        df[col] = df.groupby('region_code')[col].transform(
            lambda x: x.interpolate(method='linear', limit=14)
        )
        df[col] = df.groupby('region_code')[col].transform(
            lambda x: x.fillna(method='ffill').fillna(method='bfill')
        )
        # Fill remaining with 0 (no cases is different from no data)
        df[col] = df[col].fillna(0)
        n_missing_after = df[col].isna().sum()
        print(f"  {col}: {n_missing_before:,} → {n_missing_after:,} missing")


# =============================================================================
# FILTER VALID REGIONS
# =============================================================================

print("\n" + "-"*70)
print("Filtering Valid Regions")
print("-"*70)

# Require complete core data
df = df.dropna(subset=['temp_mean', 'deaths_elderly', 'pop_elderly'])
region_counts = df.groupby('region_code').size()
valid_regions = region_counts[region_counts >= MIN_OBS].index
df = df[df['region_code'].isin(valid_regions)]
print(f"Valid regions (≥{MIN_OBS} obs): {len(valid_regions)}")
print(f"Total observations: {len(df):,}")

# Identify regions with sufficient pollution coverage for adjusted models
poll_coverage = df.groupby('region_code').apply(
    lambda g: 1 - g[poll_vars].isna().any(axis=1).mean() if all(c in g.columns for c in poll_vars) else 0
)
pollution_valid_regions = poll_coverage[poll_coverage >= MIN_POLLUTION_COVERAGE].index
print(f"Regions with ≥{MIN_POLLUTION_COVERAGE:.0%} pollution coverage: {len(pollution_valid_regions)}")


# =============================================================================
# REGION-SPECIFIC DLNM FITTING FUNCTION
# =============================================================================

def fit_region_dlnm_supplementary(
    df_region: pd.DataFrame,
    region_code: int,
    temp_col: str = 'temp_mean',
    extra_controls: list = None,
    max_lag: int = MAX_LAG
) -> dict:
    """
    Fit DLNM for a single region with optional extra controls.
    
    Uses natural spline cross-basis and proper delta method for CIs.
    Includes population offset for rate modeling.
    
    Parameters:
    -----------
    df_region : DataFrame
        Daily data for this region
    region_code : int
        Region identifier
    temp_col : str
        Temperature column to use ('temp_mean' or 'apparent_temp')
    extra_controls : list
        Additional control variables to include (e.g., ['pm25', 'o3'])
    max_lag : int
        Maximum lag days
        
    Returns:
    --------
    dict with model results, RRs at key percentiles, and diagnostics
    """
    if len(df_region) < MIN_OBS:
        return None
    
    df_region = df_region.sort_values('date').reset_index(drop=True)
    
    # Compute region-specific temperature percentiles
    temp = df_region[temp_col].values
    temp_valid = temp[~np.isnan(temp)]
    
    if len(temp_valid) < MIN_OBS:
        return None
    
    temp_knots = np.percentile(temp_valid, TEMP_KNOT_PCTS)
    temp_boundary = (np.percentile(temp_valid, TEMP_BOUNDARY_PCTS[0]),
                     np.percentile(temp_valid, TEMP_BOUNDARY_PCTS[1]))
    
    # Create cross-basis (natural spline)
    X_cb, cb_info = create_crossbasis_ns(temp, max_lag, temp_knots, temp_boundary, LAG_KNOTS)
    
    # Find valid rows
    valid = ~np.isnan(X_cb).any(axis=1)
    if valid.sum() < MIN_OBS:
        return None
    
    X_cb_valid = X_cb[valid]
    df_valid = df_region.loc[valid].copy()
    y = df_valid['deaths_elderly'].values
    
    # Population offset (required for rate modeling)
    pop = df_valid['pop_elderly'].values
    if np.any(pop <= 0):
        pop = np.maximum(pop, 1)  # Guard against zero population
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
    
    # Add extra controls if specified
    extra_control_cols = []
    if extra_controls:
        for ctrl in extra_controls:
            if ctrl in df_valid.columns and df_valid[ctrl].notna().all():
                X_controls = np.column_stack([X_controls, df_valid[ctrl].values])
                extra_control_cols.append(ctrl)
    
    X_full = np.column_stack([X_cb_valid, X_controls])
    X_full = sm.add_constant(X_full)
    
    # Fit Quasi-Poisson GLM with population offset
    try:
        model = GLM(y, X_full, family=Poisson(), offset=offset).fit(scale='X2', maxiter=200)
        
        if not model.converged:
            print(f"    [DEBUG] Region {region_code}: GLM did not fully converge in supplementary model "
                  f"(n_valid={len(y)}, total_deaths={int(y.sum())}, dispersion={getattr(model, 'scale', float('nan')):.3f})")
        
        dispersion = model.scale
        
    except Exception as e:
        print(f"    [DEBUG] Region {region_code}: ERROR during supplementary GLM fit: {type(e).__name__}: {e}")
        return None
    
    # Extract cross-basis coefficients
    n_cb = cb_info['n_params']
    cb_coefs = model.params[1:n_cb+1]
    
    vcov_full = model.cov_params()
    if hasattr(vcov_full, 'iloc'):
        cb_vcov = vcov_full.iloc[1:n_cb+1, 1:n_cb+1].values
    else:
        cb_vcov = vcov_full[1:n_cb+1, 1:n_cb+1]
    
    # Compute temperature percentiles for this region
    temps_valid = df_valid[temp_col].dropna()
    temp_pcts = {
        'p1': float(temps_valid.quantile(0.01)),
        'p2.5': float(temps_valid.quantile(0.025)),
        'p5': float(temps_valid.quantile(0.05)),
        'p10': float(temps_valid.quantile(0.10)),
        'p25': float(temps_valid.quantile(0.25)),
        'p50': float(temps_valid.quantile(0.50)),
        'p75': float(temps_valid.quantile(0.75)),
        'p90': float(temps_valid.quantile(0.90)),
        'p95': float(temps_valid.quantile(0.95)),
        'p97.5': float(temps_valid.quantile(0.975)),
        'p99': float(temps_valid.quantile(0.99)),
    }
    
    # Find MMT
    temp_min = np.percentile(temps_valid, 1)
    temp_max = np.percentile(temps_valid, 99)
    temp_grid = np.linspace(temp_min, temp_max, 100)
    centering_temp = temp_pcts['p50']
    
    mmt, mmt_percentile, rr_curve = find_mmt_from_coefficients(
        temp_grid, cb_coefs, cb_vcov, cb_info, centering_temp
    )
    
    # Compute effects at key percentiles (relative to MMT)
    effects = compute_effects_relative_to_mmt(temp_pcts, mmt, cb_coefs, cb_vcov, cb_info)
    
    return {
        'region_code': int(region_code),
        'n_obs': int(valid.sum()),
        'n_years': float(n_years),
        'mean_deaths': float(y.mean()),
        'total_deaths': int(y.sum()),
        'dispersion': float(dispersion),
        'aic': float(model.aic),
        'temp_variable': temp_col,
        'extra_controls': extra_control_cols,
        'temp_percentiles': temp_pcts,
        'mmt': float(mmt),
        'mmt_percentile': float(mmt_percentile),
        'crossbasis_info': cb_info,
        'effects': effects,
        'cb_coefs': cb_coefs.tolist(),
        'cb_vcov': cb_vcov.tolist(),
    }


def run_analysis_and_pool(
    df: pd.DataFrame,
    analysis_name: str,
    temp_col: str = 'temp_mean',
    extra_controls: list = None,
    region_subset: list = None,
    use_parallel: bool = True,
    n_workers: int = None
) -> dict:
    """
    Run region-specific DLNM analysis and pool results via meta-analysis.
    
    Parameters:
    -----------
    df : DataFrame
        Full dataset
    analysis_name : str
        Name of this analysis (for logging)
    temp_col : str
        Temperature column
    extra_controls : list
        Extra control variables
    region_subset : list
        Subset of regions to use (e.g., for pollution-valid regions)
    use_parallel : bool
        Whether to use parallel processing
    n_workers : int
        Number of parallel workers (default: CPU count - 1)
        
    Returns:
    --------
    dict with region-level and pooled results
    """
    print(f"\n  Running: {analysis_name}")
    
    regions = list(region_subset) if region_subset is not None else list(df['region_code'].unique())
    n_regions = len(regions)
    
    if n_workers is None:
        n_workers = max(1, multiprocessing.cpu_count() - 1)
    
    region_results = {}
    n_success = 0
    n_failed = 0
    
    # Sequential processing (more reliable on Windows)
    print(f"    Fitting {n_regions} regions sequentially...")
    for i, region_code in enumerate(regions, 1):
        df_region = df[df['region_code'] == region_code].copy()
        
        result = fit_region_dlnm_supplementary(
            df_region,
            region_code,
            temp_col=temp_col,
            extra_controls=extra_controls
        )
        
        if result is not None:
            region_results[region_code] = result
            n_success += 1
        else:
            n_failed += 1
        
        if i % 20 == 0 or i == n_regions:
            print(f"      [{i}/{n_regions}] Success: {n_success}, Failed: {n_failed}")
    
    print(f"    Final: Success: {n_success}, Failed: {n_failed}")
    
    if n_success == 0:
        return {'error': 'No successful region fits', 'analysis_name': analysis_name}
    
    # Pool results via MVMeta (Gasparrini methodology)
    # This pools DLNM coefficients, not RRs at percentiles
    mvmeta_result = mvmeta_pool_regions(region_results, min_obs=MIN_OBS)
    
    # Also compute simple RR pooling for comparison
    simple_pooled = {}
    for pct in ['p1', 'p2.5', 'p97.5', 'p99']:
        simple_pooled[pct] = pool_region_results(region_results, pct)
    
    # Summary statistics
    mmts = [r['mmt'] for r in region_results.values()]
    dispersions = [r['dispersion'] for r in region_results.values()]
    
    return {
        'analysis_name': analysis_name,
        'temp_variable': temp_col,
        'extra_controls': extra_controls,
        'n_regions': n_success,
        'n_regions_failed': n_failed,
        'region_results': region_results,
        'mvmeta': mvmeta_result,  # Primary pooled results (Gasparrini methodology)
        'simple_pooled': simple_pooled,  # For comparison
        'mmt_median': float(np.median(mmts)),
        'mmt_iqr': [float(np.percentile(mmts, 25)), float(np.percentile(mmts, 75))],
        'dispersion_median': float(np.median(dispersions)),
    }


# =============================================================================
# RUN ANALYSES
# =============================================================================

print("\n" + "="*70)
print("RUNNING SUPPLEMENTARY ANALYSES")
print("="*70)

all_results = {}

# -----------------------------------------------------------------------------
# 1. BASE MODEL (Dry-bulb temperature, no extra controls)
# -----------------------------------------------------------------------------
print("\n[1/5] Base Model (Dry-bulb, no pollution/flu adjustment)")
all_results['base'] = run_analysis_and_pool(
    df, 
    analysis_name='Base (Dry-bulb)',
    temp_col='temp_mean',
    extra_controls=None
)

# -----------------------------------------------------------------------------
# 2. APPARENT TEMPERATURE MODEL
# -----------------------------------------------------------------------------
print("\n[2/5] Apparent Temperature Model")
all_results['apparent_temp'] = run_analysis_and_pool(
    df,
    analysis_name='Apparent Temperature',
    temp_col='apparent_temp',
    extra_controls=None
)

# -----------------------------------------------------------------------------
# 3. POLLUTION-ADJUSTED MODEL
# -----------------------------------------------------------------------------
print("\n[3/5] Pollution-Adjusted Model (PM2.5 + O3)")
# Only use regions with sufficient pollution coverage
pollution_controls = [c for c in ['pm25', 'o3'] if c in df.columns]
if pollution_controls and len(pollution_valid_regions) > 0:
    all_results['pollution_adjusted'] = run_analysis_and_pool(
        df,
        analysis_name='Pollution-Adjusted',
        temp_col='temp_mean',
        extra_controls=pollution_controls,
        region_subset=list(pollution_valid_regions)
    )
else:
    print("  Skipped: insufficient pollution data")
    all_results['pollution_adjusted'] = {'error': 'Insufficient pollution data'}

# -----------------------------------------------------------------------------
# 4. INFLUENZA-ADJUSTED MODEL
# -----------------------------------------------------------------------------
print("\n[4/5] Influenza-Adjusted Model")
flu_controls = [c for c in ['srag_cases', 'confirmed_flu'] if c in df.columns]
if flu_controls:
    all_results['flu_adjusted'] = run_analysis_and_pool(
        df,
        analysis_name='Influenza-Adjusted',
        temp_col='temp_mean',
        extra_controls=flu_controls
    )
else:
    print("  Skipped: no influenza data")
    all_results['flu_adjusted'] = {'error': 'No influenza data'}

# -----------------------------------------------------------------------------
# 5. FULLY ADJUSTED MODEL (Pollution + Influenza)
# -----------------------------------------------------------------------------
print("\n[5/5] Fully Adjusted Model (Pollution + Influenza)")
all_controls = pollution_controls + flu_controls
if all_controls and len(pollution_valid_regions) > 0:
    all_results['fully_adjusted'] = run_analysis_and_pool(
        df,
        analysis_name='Fully Adjusted',
        temp_col='temp_mean',
        extra_controls=all_controls,
        region_subset=list(pollution_valid_regions)
    )
else:
    print("  Skipped: insufficient data")
    all_results['fully_adjusted'] = {'error': 'Insufficient data for full adjustment'}


# =============================================================================
# COMPARISON SUMMARY (Using MVMeta - Gasparrini Methodology)
# =============================================================================

print("\n" + "="*70)
print("COMPARISON SUMMARY (MVMeta - Gasparrini Methodology)")
print("="*70)

def format_mvmeta_rr(mvmeta_result, percentile):
    """Format RR with 95% CI from MVMeta for display."""
    if 'error' in mvmeta_result:
        return "Error"
    if 'effects' not in mvmeta_result or percentile not in mvmeta_result['effects']:
        return "N/A"
    eff = mvmeta_result['effects'][percentile]
    return f"{eff['rr']:.3f} ({eff['rr_lower']:.3f}-{eff['rr_upper']:.3f})"

print("\n" + "-"*70)
print("HEAT EFFECTS (P99 vs Pooled MMT)")
print("-"*70)
print(f"{'Analysis':<30} {'RR (95% CI)':<25} {'N Regions':<10} {'Pooled MMT':<12}")
print("-"*70)
for name, result in all_results.items():
    if 'error' in result:
        print(f"{name:<30} {'Error: ' + result['error'][:40]}")
        continue
    mvmeta = result.get('mvmeta', {})
    rr_str = format_mvmeta_rr(mvmeta, 'p99')
    n_reg = mvmeta.get('n_regions', result.get('n_regions', 0))
    mmt = mvmeta.get('mmt', np.nan)
    mmt_pct = mvmeta.get('mmt_percentile', np.nan)
    mmt_str = f"P{mmt_pct:.0f}" if not np.isnan(mmt_pct) else "N/A"
    print(f"{result['analysis_name']:<30} {rr_str:<25} {n_reg:<10} {mmt_str:<12}")

print("\n" + "-"*70)
print("COLD EFFECTS (P1 vs Pooled MMT)")
print("-"*70)
print(f"{'Analysis':<30} {'RR (95% CI)':<25} {'N Regions':<10} {'Pooled MMT':<12}")
print("-"*70)
for name, result in all_results.items():
    if 'error' in result:
        continue
    mvmeta = result.get('mvmeta', {})
    rr_str = format_mvmeta_rr(mvmeta, 'p1')
    n_reg = mvmeta.get('n_regions', result.get('n_regions', 0))
    mmt = mvmeta.get('mmt', np.nan)
    mmt_pct = mvmeta.get('mmt_percentile', np.nan)
    mmt_str = f"P{mmt_pct:.0f}" if not np.isnan(mmt_pct) else "N/A"
    print(f"{result['analysis_name']:<30} {rr_str:<25} {n_reg:<10} {mmt_str:<12}")


# =============================================================================
# SAVE RESULTS
# =============================================================================

print("\n" + "-"*70)
print("Saving Results")
print("-"*70)

# Use level suffix for output files
suffix = '' if level == 'intermediate' else '_immediate'

# Full results (JSON)
results_file = OUTPUT_DIR / f'supplementary_analyses_v2{suffix}.json'
with open(results_file, 'w') as f:
    json.dump(convert_to_json_serializable(all_results), f, indent=2)
print(f"Full results: {results_file}")

# Summary table (CSV) - using MVMeta results
summary_rows = []
for name, result in all_results.items():
    if 'error' in result:
        continue
    mvmeta = result.get('mvmeta', {})
    if 'error' in mvmeta or 'effects' not in mvmeta:
        continue
    for pct in ['p1', 'p2.5', 'p97.5', 'p99']:
        if pct in mvmeta['effects']:
            eff = mvmeta['effects'][pct]
            summary_rows.append({
                'analysis': result['analysis_name'],
                'temp_variable': result.get('temp_variable', 'temp_mean'),
                'extra_controls': ', '.join(result.get('extra_controls', []) or []),
                'pooling_method': 'MVMeta',
                'percentile': pct,
                'pooled_rr': eff.get('rr'),
                'pooled_rr_lower': eff.get('rr_lower'),
                'pooled_rr_upper': eff.get('rr_upper'),
                'log_rr': eff.get('log_rr'),
                'log_rr_se': eff.get('log_rr_se'),
                'n_regions': mvmeta.get('n_regions'),
                'pooled_mmt': mvmeta.get('mmt'),
                'pooled_mmt_percentile': mvmeta.get('mmt_percentile'),
            })

summary_df = pd.DataFrame(summary_rows)
summary_file = OUTPUT_DIR / f'supplementary_analyses_v2_summary{suffix}.csv'
summary_df.to_csv(summary_file, index=False)
print(f"Summary table: {summary_file}")

# Metadata
metadata = {
    'script': '03a_supplementary_analyses_v2.py',
    'timestamp': datetime.now().isoformat(),
    'level': level,
    'configuration': {
        'max_lag': MAX_LAG,
        'temp_df': TEMP_DF,
        'lag_df': LAG_DF,
        'min_obs': MIN_OBS,
        'min_pollution_coverage': MIN_POLLUTION_COVERAGE,
        'temp_knot_pcts': TEMP_KNOT_PCTS,
        'lag_knots': LAG_KNOTS,
    },
    'data': {
        'total_observations': len(df),
        'n_regions_all': df['region_code'].nunique(),
        'n_regions_pollution_valid': len(pollution_valid_regions),
    },
    'analyses_run': list(all_results.keys()),
}

metadata_file = OUTPUT_DIR / f'supplementary_analyses_v2_metadata{suffix}.json'
with open(metadata_file, 'w') as f:
    json.dump(metadata, f, indent=2)
print(f"Metadata: {metadata_file}")

print("\n" + "="*70)
print("PHASE 3 v2 COMPLETE")
print("="*70)
print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
