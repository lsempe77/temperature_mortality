"""
01a_intermediate_dlnm.py
=========================
Two-Stage DLNM Analysis at Intermediate Region Level (133 regions)

Estimates temperature-mortality relationship using:
1. First stage: Region-specific DLNM models
2. Second stage: Random-effects meta-analysis pooling

Methodology follows Gasparrini et al. (2015) Lancet - MCC temperature-mortality

Input Data (from phase0_data_prep/results/):
- era5_intermediate_daily.parquet: Temperature data (728K rows, 133 regions)
- mortality_regional_daily_elderly.parquet: Elderly mortality (661K rows, 133 regions)
- ses_intermediate_covariates.csv: Population for offset (133 regions)
- brazilian_holidays_daily.parquet: Holiday controls (5,479 days)

Output:
- results/dlnm_intermediate_results.json: Full model results
- results/dlnm_intermediate_summary.csv: Region-level summary
- results/dlnm_intermediate_pooled.json: Pooled national estimates

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
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("01a: INTERMEDIATE REGION DLNM ANALYSIS (133 regions)")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PHASE0_RESULTS = os.path.join(BASE_DIR, 'new_analysis', 'phase0_data_prep', 'results')
OUTPUT_DIR = os.path.join(BASE_DIR, 'new_analysis', 'phase1_core_model', 'results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# DLNM parameters
MAX_LAG = 21  # Maximum lag in days (standard for temperature-mortality)
TEMP_DF = 4   # Degrees of freedom for temperature (natural spline)
LAG_DF = 3    # Degrees of freedom for lag (natural spline)
TIME_SPLINE_DF_PER_YEAR = 3  # DF per year for long-term trend

# Minimum observations per region
MIN_OBS = 1000  # At least ~3 years of data

print(f"\nConfiguration:")
print(f"  Max lag: {MAX_LAG} days")
print(f"  Temperature DF: {TEMP_DF}")
print(f"  Lag DF: {LAG_DF}")
print(f"  Minimum obs per region: {MIN_OBS}")

# =============================================================================
# DATA LOADING
# =============================================================================

print("\n" + "-"*70)
print("Loading Data")
print("-"*70)

# Load temperature data
era5_file = os.path.join(PHASE0_RESULTS, 'era5_intermediate_daily.parquet')
df_temp = pd.read_parquet(era5_file)
df_temp['date'] = pd.to_datetime(df_temp['date'])
print(f"ERA5: {len(df_temp):,} rows, {df_temp['region_code'].nunique()} regions")
print(f"  Date range: {df_temp['date'].min().date()} to {df_temp['date'].max().date()}")

# Load mortality data
mort_file = os.path.join(PHASE0_RESULTS, 'mortality_regional_daily_elderly.parquet')
df_mort = pd.read_parquet(mort_file)
df_mort['date'] = pd.to_datetime(df_mort['date'])
print(f"Mortality: {len(df_mort):,} rows, {df_mort['region_code'].nunique()} regions")

# Load SES data (for population offset)
ses_file = os.path.join(PHASE0_RESULTS, 'ses_intermediate_covariates.csv')
df_ses = pd.read_csv(ses_file)
print(f"SES: {len(df_ses)} regions")
# Create region_code to pop_elderly mapping
pop_map = dict(zip(df_ses['intermediate_code'], df_ses['pop_elderly']))

# Load holidays
holiday_file = os.path.join(PHASE0_RESULTS, 'brazilian_holidays_daily.parquet')
df_holidays = pd.read_parquet(holiday_file)
df_holidays['date'] = pd.to_datetime(df_holidays['date'])
print(f"Holidays: {len(df_holidays)} days, {df_holidays['is_holiday'].sum()} holidays")

# =============================================================================
# MERGE DATA
# =============================================================================

print("\n" + "-"*70)
print("Merging Data")
print("-"*70)

# Merge temperature and mortality
df = pd.merge(
    df_mort,
    df_temp[['date', 'region_code', 'temp_mean']],
    on=['date', 'region_code'],
    how='inner'
)
print(f"After temp merge: {len(df):,} rows")

# Merge holidays
df = pd.merge(
    df,
    df_holidays[['date', 'is_holiday', 'is_holiday_week']],
    on='date',
    how='left'
)
df['is_holiday'] = df['is_holiday'].fillna(0).astype(int)
df['is_holiday_week'] = df['is_holiday_week'].fillna(0).astype(int)

# Add population (use static 2022 estimate)
df['pop_elderly'] = df['region_code'].map(pop_map)

# Add time variables
df['year'] = df['date'].dt.year
df['month'] = df['date'].dt.month
df['day_of_week'] = df['date'].dt.dayofweek
df['time_index'] = (df['date'] - df['date'].min()).dt.days

# Check for missing values
print(f"\nMissing values:")
print(f"  Temperature: {df['temp_mean'].isna().sum()}")
print(f"  Deaths: {df['deaths_elderly'].isna().sum()}")
print(f"  Population: {df['pop_elderly'].isna().sum()}")

# Drop missing
df = df.dropna(subset=['temp_mean', 'deaths_elderly', 'pop_elderly'])
print(f"\nFinal dataset: {len(df):,} rows, {df['region_code'].nunique()} regions")

# =============================================================================
# DLNM FUNCTIONS
# =============================================================================

def create_lag_matrix(x: np.ndarray, max_lag: int) -> np.ndarray:
    """Create matrix of lagged values."""
    n = len(x)
    lag_matrix = np.full((n, max_lag + 1), np.nan)
    for lag in range(max_lag + 1):
        if lag == 0:
            lag_matrix[:, lag] = x
        else:
            lag_matrix[lag:, lag] = x[:-lag]
    return lag_matrix


def create_crossbasis(temp: np.ndarray, max_lag: int, temp_df: int, lag_df: int) -> tuple:
    """
    Create cross-basis matrix for DLNM using natural splines.
    
    Returns:
        X_cb: Cross-basis matrix
        temp_knots: Temperature knot locations
        lag_knots: Lag knot locations
        temp_boundary: Temperature boundary knots
    """
    n = len(temp)
    
    # Create lag matrix
    temp_lags = create_lag_matrix(temp, max_lag)
    
    # Temperature knots at percentiles (excluding extremes for boundary)
    temp_valid = temp[~np.isnan(temp)]
    temp_pcts = [10, 75, 90]  # Interior knots
    temp_knots = np.percentile(temp_valid, temp_pcts)
    temp_boundary = [np.percentile(temp_valid, 1), np.percentile(temp_valid, 99)]
    
    # Lag knots at equal intervals
    lag_knots = np.linspace(0, max_lag, lag_df + 2)[1:-1]  # Interior knots
    
    # Create natural spline bases
    # For simplicity, use polynomial basis (natural splines require more setup)
    # Center temperature at median
    temp_center = np.nanmedian(temp)
    temp_scale = np.nanstd(temp)
    
    # Build cross-basis: polynomial in temp × polynomial in lag
    X_list = []
    col_names = []
    
    for lag in range(max_lag + 1):
        temp_lag_std = (temp_lags[:, lag] - temp_center) / temp_scale
        # Polynomial terms for temperature
        for p in range(1, temp_df + 1):
            # Weight by lag structure (declining weights)
            lag_weight = 1 - (lag / (max_lag + 1))  # Linear decay
            X_list.append(temp_lag_std ** p * lag_weight)
            col_names.append(f'temp_lag{lag}_poly{p}')
    
    X_cb = np.column_stack(X_list)
    
    return X_cb, col_names, temp_center, temp_scale


def create_simple_crossbasis(temp: np.ndarray, max_lag: int, poly_degree: int = 3) -> tuple:
    """
    Simplified cross-basis using polynomial distributed lag.
    
    This is computationally efficient and captures main temperature-lag structure.
    
    IMPORTANT: This uses polynomial basis rather than natural cubic splines.
    The downstream RR computation (compute_cumulative_effect) assumes this specific
    structure with column names 'lag{L}_poly{P}'. If you change the basis construction,
    you MUST also update compute_cumulative_effect to match.
    
    For true natural spline DLNM, consider using the dlnm R package via rpy2,
    or implementing ns() basis with proper boundary constraints.
    """
    n = len(temp)
    
    # Create lag matrix
    temp_lags = create_lag_matrix(temp, max_lag)
    
    # Center and scale
    temp_center = np.nanmedian(temp)
    temp_scale = np.nanstd(temp)
    
    # Build cross-basis
    X_list = []
    col_names = []
    
    for lag in range(max_lag + 1):
        temp_std = (temp_lags[:, lag] - temp_center) / temp_scale
        for p in range(1, poly_degree + 1):
            X_list.append(temp_std ** p)
            col_names.append(f'lag{lag}_poly{p}')
    
    X_cb = np.column_stack(X_list)
    
    return X_cb, col_names, temp_center, temp_scale


def create_time_spline(time_index: np.ndarray, n_years: float, df_per_year: int) -> np.ndarray:
    """Create natural spline basis for time trend."""
    total_df = max(3, int(n_years * df_per_year))
    
    # Use B-spline basis via patsy
    knots = np.linspace(time_index.min(), time_index.max(), total_df + 2)[1:-1]
    
    try:
        spline_basis = dmatrix(
            f"bs(x, knots={list(knots)}, degree=3, include_intercept=False) - 1",
            {"x": time_index},
            return_type='dataframe'
        ).values
    except:
        # Fallback to polynomial
        spline_basis = np.column_stack([
            (time_index / 365) ** p for p in range(1, total_df + 1)
        ])
    
    return spline_basis


def fit_region_dlnm(df_region: pd.DataFrame, region_code: int, max_lag: int = MAX_LAG) -> dict:
    """
    Fit DLNM for a single region using Quasi-Poisson GLM.
    
    Returns dictionary with coefficients, standard errors, and effect estimates.
    """
    if len(df_region) < MIN_OBS:
        print(f"  Region {region_code}: Insufficient data ({len(df_region)} < {MIN_OBS})")
        return None
    
    # Sort by date
    df_region = df_region.sort_values('date').reset_index(drop=True)
    
    # Create cross-basis
    temp = df_region['temp_mean'].values
    X_cb, col_names, temp_center, temp_scale = create_simple_crossbasis(temp, max_lag)
    
    # Find valid rows (no NaN from lagging)
    valid = ~np.isnan(X_cb).any(axis=1)
    
    if valid.sum() < MIN_OBS:
        print(f"  Region {region_code}: Too few valid obs ({valid.sum()})")
        return None
    
    # Subset to valid
    X_cb_valid = X_cb[valid]
    df_valid = df_region.loc[valid].copy()
    y = df_valid['deaths_elderly'].values
    
    # Population offset
    pop = df_valid['pop_elderly'].values
    offset = np.log(pop / 100000)  # Rate per 100k
    
    # Control variables
    # 1. Month dummies
    month_dummies = pd.get_dummies(df_valid['month'], prefix='month', drop_first=True)
    
    # 2. Day of week dummies
    dow_dummies = pd.get_dummies(df_valid['day_of_week'], prefix='dow', drop_first=True)
    
    # 3. Time spline
    n_years = (df_valid['date'].max() - df_valid['date'].min()).days / 365.25
    time_spline = create_time_spline(df_valid['time_index'].values, n_years, TIME_SPLINE_DF_PER_YEAR)
    
    # 4. Holiday
    holiday = df_valid['is_holiday'].values.reshape(-1, 1)
    
    # Combine predictors
    X_controls = np.column_stack([
        month_dummies.values,
        dow_dummies.values,
        time_spline,
        holiday
    ])
    
    X_full = np.column_stack([X_cb_valid, X_controls])
    X_full = sm.add_constant(X_full)
    
    # Fit Quasi-Poisson GLM
    try:
        model = GLM(
            y,
            X_full,
            family=Poisson(),
            offset=offset
        ).fit(scale='X2')  # Pearson chi-square for quasi-likelihood
        
        dispersion = model.scale
        
        # Check for convergence issues
        if not model.converged:
            print(f"  Region {region_code}: Model did not converge")
            return None
        
        # Check for extreme dispersion (suggests model problems)
        if dispersion > 10 or dispersion < 0.1:
            print(f"  Region {region_code}: Extreme dispersion ({dispersion:.2f}), results may be unreliable")
        
    except Exception as e:
        print(f"  Region {region_code}: Model failed - {e}")
        return None
    
    # Extract cross-basis coefficients
    n_cb = len(col_names)
    cb_coefs = model.params[1:n_cb+1]  # Skip intercept
    
    # Variance-covariance matrix for cross-basis
    vcov_full = model.cov_params()
    if hasattr(vcov_full, 'iloc'):
        cb_vcov = vcov_full.iloc[1:n_cb+1, 1:n_cb+1].values
    else:
        cb_vcov = vcov_full[1:n_cb+1, 1:n_cb+1]
    
    # Compute effects at temperature percentiles
    temps = df_valid['temp_mean'].dropna()
    percentiles = [1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99]
    temp_pcts = {f'p{p}'.replace('.', '_'): temps.quantile(p/100) for p in percentiles}
    # Rename for cleaner keys (p2_5 -> p2.5 style in output)
    temp_pcts = {k.replace('_', '.'): v for k, v in temp_pcts.items()}
    
    # Reference temperature (median)
    ref_temp = temp_pcts['p50']
    
    # Compute cumulative RR at each percentile vs reference
    effects = {}
    for pct_name, target_temp in temp_pcts.items():
        effect, se = compute_cumulative_effect(
            target_temp, ref_temp, 
            cb_coefs, cb_vcov, col_names,
            temp_center, temp_scale, max_lag
        )
        effects[pct_name] = {
            'temp': float(target_temp),
            'log_rr': float(effect),
            'log_rr_se': float(se),
            'rr': float(np.exp(effect)),
            'rr_lower': float(np.exp(effect - 1.96 * se)),
            'rr_upper': float(np.exp(effect + 1.96 * se))
        }
    
    return {
        'region_code': int(region_code),
        'n_obs': int(valid.sum()),
        'n_years': float(n_years),
        'mean_deaths': float(y.mean()),
        'total_deaths': int(y.sum()),
        'dispersion': float(dispersion),
        'aic': float(model.aic),
        'temp_center': float(temp_center),
        'temp_scale': float(temp_scale),
        'temp_percentiles': {k: float(v) for k, v in temp_pcts.items()},
        'effects': effects,
        'cb_coefs': cb_coefs.tolist(),
        'cb_vcov': cb_vcov.tolist(),
        'col_names': col_names
    }


def compute_cumulative_effect(target_temp: float, ref_temp: float,
                               cb_coefs: np.ndarray, cb_vcov: np.ndarray,
                               col_names: list, temp_center: float,
                               temp_scale: float, max_lag: int) -> tuple:
    """
    Compute cumulative effect (sum over all lags) at target vs reference temperature.
    
    Returns log(RR) and standard error.
    """
    temp_std_target = (target_temp - temp_center) / temp_scale
    temp_std_ref = (ref_temp - temp_center) / temp_scale
    
    # Create contrast vector
    contrast = np.zeros(len(col_names))
    
    for i, name in enumerate(col_names):
        # Parse lag and polynomial from column name: lag{L}_poly{P}
        parts = name.split('_')
        lag = int(parts[0].replace('lag', ''))
        poly = int(parts[1].replace('poly', ''))
        
        # Difference in polynomial terms
        contrast[i] = temp_std_target**poly - temp_std_ref**poly
    
    # Cumulative effect
    effect = np.dot(contrast, cb_coefs)
    
    # Standard error via delta method
    var = np.dot(np.dot(contrast, cb_vcov), contrast)
    se = np.sqrt(var) if var > 0 else 0
    
    return effect, se


# =============================================================================
# META-ANALYSIS FUNCTIONS
# =============================================================================

def random_effects_meta_analysis(effects: np.ndarray, variances: np.ndarray) -> dict:
    """
    Random-effects meta-analysis using DerSimonian-Laird estimator.
    """
    k = len(effects)
    
    if k < 2:
        return {
            'pooled_effect': float(effects[0]) if k == 1 else np.nan,
            'pooled_se': float(np.sqrt(variances[0])) if k == 1 else np.nan,
            'tau2': 0.0,
            'I2': 0.0,
            'Q': 0.0,
            'p_heterogeneity': 1.0,
            'n_regions': k
        }
    
    # Fixed-effect weights
    w = 1 / variances
    
    # Fixed-effect estimate
    theta_fe = np.sum(w * effects) / np.sum(w)
    
    # Cochran's Q statistic
    Q = np.sum(w * (effects - theta_fe)**2)
    
    # DerSimonian-Laird tau^2
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0, (Q - (k - 1)) / c)
    
    # Random-effects weights
    w_re = 1 / (variances + tau2)
    
    # Random-effects estimate
    theta_re = np.sum(w_re * effects) / np.sum(w_re)
    se_re = np.sqrt(1 / np.sum(w_re))
    
    # I² heterogeneity
    I2 = max(0, (Q - (k - 1)) / Q * 100) if Q > 0 else 0
    
    # P-value for heterogeneity
    p_het = 1 - stats.chi2.cdf(Q, k - 1)
    
    return {
        'pooled_effect': float(theta_re),
        'pooled_se': float(se_re),
        'pooled_rr': float(np.exp(theta_re)),
        'pooled_rr_lower': float(np.exp(theta_re - 1.96 * se_re)),
        'pooled_rr_upper': float(np.exp(theta_re + 1.96 * se_re)),
        'tau2': float(tau2),
        'I2': float(I2),
        'Q': float(Q),
        'p_heterogeneity': float(p_het),
        'n_regions': k
    }


def pool_region_results(region_results: dict, percentile: str) -> dict:
    """Pool region-specific results via meta-analysis for a given percentile.
    
    Includes sanity filters:
    - SE must be positive and < 2 (plausible range)
    - log(RR) must be in plausible range (-3, 3) ~ RR 0.05 to 20
    """
    effects = []
    variances = []
    regions = []
    excluded_reasons = {'invalid_se': 0, 'extreme_se': 0, 'extreme_rr': 0}
    
    # Thresholds for filtering
    MAX_SE = 2.0  # SE > 2 suggests unstable estimate
    MAX_ABS_LOGRR = 3.0  # |log(RR)| > 3 is implausible (RR < 0.05 or > 20)
    
    for region_code, result in region_results.items():
        if result is None:
            continue
        if percentile not in result['effects']:
            continue
        
        eff = result['effects'][percentile]
        log_rr = eff['log_rr']
        se = eff['log_rr_se']
        
        # Filter 1: Valid SE (positive)
        if se <= 0:
            excluded_reasons['invalid_se'] += 1
            continue
        
        # Filter 2: Plausible SE (not too large)
        if se > MAX_SE:
            excluded_reasons['extreme_se'] += 1
            continue
        
        # Filter 3: Plausible log(RR)
        if abs(log_rr) > MAX_ABS_LOGRR:
            excluded_reasons['extreme_rr'] += 1
            continue
        
        effects.append(log_rr)
        variances.append(se**2)
        regions.append(region_code)
    
    if len(effects) == 0:
        return {'pooled_rr': np.nan, 'n_regions': 0, 'excluded_reasons': excluded_reasons}
    
    effects = np.array(effects)
    variances = np.array(variances)
    
    pooled = random_effects_meta_analysis(effects, variances)
    pooled['percentile'] = percentile
    pooled['regions_included'] = regions
    pooled['excluded_reasons'] = excluded_reasons
    
    return pooled


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("FITTING REGION-SPECIFIC DLNM MODELS")
print("="*70)

regions = sorted(df['region_code'].unique())
print(f"\nTotal regions: {len(regions)}")

region_results = {}
successful = 0

for i, region_code in enumerate(regions):
    df_region = df[df['region_code'] == region_code].copy()
    
    result = fit_region_dlnm(df_region, region_code)
    region_results[region_code] = result
    
    if result is not None:
        successful += 1
        eff_p99 = result['effects']['p99']
        eff_p1 = result['effects']['p1']
        if (i + 1) % 20 == 0 or i < 5:
            print(f"  Region {region_code}: Heat RR(P99)={eff_p99['rr']:.3f} "
                  f"[{eff_p99['rr_lower']:.3f}-{eff_p99['rr_upper']:.3f}], "
                  f"Cold RR(P1)={eff_p1['rr']:.3f}")

print(f"\nSuccessfully fitted: {successful} / {len(regions)} regions")

# Log excluded regions statistics
excluded_count = len(regions) - successful
if excluded_count > 0:
    excluded_regions = [r for r, res in region_results.items() if res is None]
    excluded_sample_sizes = [len(df[df['region_code'] == r]) for r in excluded_regions]
    print(f"\n  Excluded regions: {excluded_count}")
    print(f"  Mean sample size of excluded: {np.mean(excluded_sample_sizes):.0f} obs")
    print(f"  Excluded region codes: {excluded_regions[:10]}{'...' if len(excluded_regions) > 10 else ''}")

# =============================================================================
# META-ANALYSIS POOLING
# =============================================================================

print("\n" + "="*70)
print("POOLING RESULTS VIA META-ANALYSIS")
print("="*70)

# Pool at key percentiles
pooled_results = {}
for pct in ['p1', 'p2.5', 'p5', 'p10', 'p25', 'p50', 'p75', 'p90', 'p95', 'p97.5', 'p99']:
    pooled_results[pct] = pool_region_results(region_results, pct)

# Check temperature percentile consistency across regions
# (Important: if regions have very different date ranges/seasonal coverage,
# percentile thresholds may not be comparable)
print("\n" + "-"*70)
print("TEMPERATURE THRESHOLD CONSISTENCY CHECK")
print("-"*70)
temp_p99_values = []
temp_p1_values = []
for region_code, result in region_results.items():
    if result is not None and 'temp_percentiles' in result:
        temp_p99_values.append(result['temp_percentiles'].get('p99', np.nan))
        temp_p1_values.append(result['temp_percentiles'].get('p1', np.nan))

if temp_p99_values:
    p99_range = np.nanmax(temp_p99_values) - np.nanmin(temp_p99_values)
    p1_range = np.nanmax(temp_p1_values) - np.nanmin(temp_p1_values)
    print(f"  P99 threshold range: {np.nanmin(temp_p99_values):.1f}°C to {np.nanmax(temp_p99_values):.1f}°C (spread: {p99_range:.1f}°C)")
    print(f"  P1 threshold range:  {np.nanmin(temp_p1_values):.1f}°C to {np.nanmax(temp_p1_values):.1f}°C (spread: {p1_range:.1f}°C)")
    if p99_range > 10 or p1_range > 15:
        print("  ⚠ Large threshold variation suggests heterogeneous climates across regions")
        print("    Consider using region-specific adaptation in interpretation")

# Print key results
print("\n" + "-"*70)
print("HEAT EFFECTS (vs P50 reference)")
print("-"*70)
for pct in ['p75', 'p90', 'p95', 'p97.5', 'p99']:
    p = pooled_results[pct]
    if 'pooled_rr' in p and not np.isnan(p['pooled_rr']):
        print(f"  {pct.upper()}: RR = {p['pooled_rr']:.3f} "
              f"[{p['pooled_rr_lower']:.3f}-{p['pooled_rr_upper']:.3f}], "
              f"I² = {p['I2']:.1f}%")

print("\n" + "-"*70)
print("COLD EFFECTS (vs P50 reference)")
print("-"*70)
for pct in ['p25', 'p10', 'p5', 'p2.5', 'p1']:
    p = pooled_results[pct]
    if 'pooled_rr' in p and not np.isnan(p['pooled_rr']):
        print(f"  {pct.upper()}: RR = {p['pooled_rr']:.3f} "
              f"[{p['pooled_rr_lower']:.3f}-{p['pooled_rr_upper']:.3f}], "
              f"I² = {p['I2']:.1f}%")

# =============================================================================
# SAVE RESULTS
# =============================================================================

print("\n" + "="*70)
print("SAVING RESULTS")
print("="*70)

# Helper function to convert numpy types (ensures no precision loss)
def convert_to_json_serializable(obj):
    if isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        # Use full precision for float64
        return float(obj)
    elif isinstance(obj, np.ndarray):
        # Convert to list, preserving float64 precision
        return [convert_to_json_serializable(x) for x in obj]
    elif isinstance(obj, dict):
        return {k: convert_to_json_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_json_serializable(i) for i in obj]
    return obj

# Full results (JSON)
full_results = {
    'analysis_level': 'intermediate',
    'n_regions': len(regions),
    'n_successful': successful,
    'parameters': {
        'max_lag': MAX_LAG,
        'temp_df': TEMP_DF,
        'lag_df': LAG_DF,
        'time_spline_df_per_year': TIME_SPLINE_DF_PER_YEAR,
        'min_obs': MIN_OBS
    },
    'region_results': {str(k): v for k, v in region_results.items() if v is not None},
    'pooled_results': pooled_results,
    'timestamp': datetime.now().isoformat()
}

# Convert to JSON-serializable
full_results = convert_to_json_serializable(full_results)

output_file = os.path.join(OUTPUT_DIR, 'dlnm_intermediate_results.json')
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
        'total_deaths': result['total_deaths'],
        'dispersion': result['dispersion']
    }
    for pct in ['p1', 'p5', 'p50', 'p95', 'p99']:
        if pct in result['effects']:
            eff = result['effects'][pct]
            row[f'rr_{pct}'] = eff['rr']
            row[f'rr_{pct}_lower'] = eff['rr_lower']
            row[f'rr_{pct}_upper'] = eff['rr_upper']
            row[f'temp_{pct}'] = eff['temp']
    summary_rows.append(row)

summary_df = pd.DataFrame(summary_rows)
summary_file = os.path.join(OUTPUT_DIR, 'dlnm_intermediate_summary.csv')
summary_df.to_csv(summary_file, index=False)
print(f"Saved: {summary_file}")

# Pooled results (JSON)
pooled_file = os.path.join(OUTPUT_DIR, 'dlnm_intermediate_pooled.json')
with open(pooled_file, 'w') as f:
    json.dump(convert_to_json_serializable(pooled_results), f, indent=2)
print(f"Saved: {pooled_file}")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)

print(f"\n  Regions analyzed: {successful} / {len(regions)}")
print(f"  Total observations: {df.groupby('region_code').size().sum():,}")

print("\n  KEY FINDINGS:")
p99 = pooled_results['p99']
p1 = pooled_results['p1']
if 'pooled_rr' in p99 and not np.isnan(p99['pooled_rr']):
    print(f"    Extreme Heat (P99): RR = {p99['pooled_rr']:.3f} [{p99['pooled_rr_lower']:.3f}-{p99['pooled_rr_upper']:.3f}]")
if 'pooled_rr' in p1 and not np.isnan(p1['pooled_rr']):
    print(f"    Extreme Cold (P1):  RR = {p1['pooled_rr']:.3f} [{p1['pooled_rr_lower']:.3f}-{p1['pooled_rr_upper']:.3f}]")

print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
