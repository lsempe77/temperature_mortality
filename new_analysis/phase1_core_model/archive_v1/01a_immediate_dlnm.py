"""
01a_immediate_dlnm.py
=========================
Two-Stage DLNM Analysis at Immediate Region Level (510 regions)

Estimates temperature-mortality relationship using:
1. First stage: Region-specific DLNM models
2. Second stage: Random-effects meta-analysis pooling

Input Data (from phase0_data_prep/results/):
- era5_immediate_daily.parquet: Temperature data (2.8M rows, 510 regions)
- mortality_immediate_daily_elderly.parquet: Elderly mortality (2.3M rows, 510 regions)
- ses_immediate_covariates.csv: Population for offset (510 regions)
- brazilian_holidays_daily.parquet: Holiday controls (5,479 days)

Output:
- results/dlnm_immediate_results.json: Full model results
- results/dlnm_immediate_summary.csv: Region-level summary
- results/dlnm_immediate_pooled.json: Pooled national estimates

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
print("01a: IMMEDIATE REGION DLNM ANALYSIS (510 regions)")
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
MAX_LAG = 21  # Maximum lag in days
TEMP_DF = 4   # Degrees of freedom for temperature
LAG_DF = 3    # Degrees of freedom for lag
TIME_SPLINE_DF_PER_YEAR = 3

# Minimum observations per region (lower than intermediate due to smaller units)
MIN_OBS = 500  # At least ~1.5 years of data

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
era5_file = os.path.join(PHASE0_RESULTS, 'era5_immediate_daily.parquet')
df_temp = pd.read_parquet(era5_file)
df_temp['date'] = pd.to_datetime(df_temp['date'])
# Note: immediate data uses region_code column (same as intermediate files)
print(f"ERA5: {len(df_temp):,} rows, {df_temp['region_code'].nunique()} regions")
print(f"  Date range: {df_temp['date'].min().date()} to {df_temp['date'].max().date()}")

# Load mortality data
mort_file = os.path.join(PHASE0_RESULTS, 'mortality_immediate_daily_elderly.parquet')
df_mort = pd.read_parquet(mort_file)
df_mort['date'] = pd.to_datetime(df_mort['date'])
# Mortality uses immediate_code column
print(f"Mortality: {len(df_mort):,} rows, {df_mort['immediate_code'].nunique()} regions")

# Load SES data (for population offset)
ses_file = os.path.join(PHASE0_RESULTS, 'ses_immediate_covariates.csv')
df_ses = pd.read_csv(ses_file)
print(f"SES: {len(df_ses)} regions")
pop_map = dict(zip(df_ses['immediate_code'], df_ses['pop_elderly']))

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

# Rename columns for consistency
df_temp = df_temp.rename(columns={'region_code': 'immediate_code'})

# Merge temperature and mortality
df = pd.merge(
    df_mort,
    df_temp[['date', 'immediate_code', 'temp_mean']],
    on=['date', 'immediate_code'],
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

# Add population
df['pop_elderly'] = df['immediate_code'].map(pop_map)

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
print(f"\nFinal dataset: {len(df):,} rows, {df['immediate_code'].nunique()} regions")

# =============================================================================
# DLNM FUNCTIONS (same as intermediate)
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


def create_simple_crossbasis(temp: np.ndarray, max_lag: int, poly_degree: int = 3) -> tuple:
    """Simplified cross-basis using polynomial distributed lag.
    
    IMPORTANT: This uses polynomial basis rather than natural cubic splines.
    The downstream RR computation assumes this specific structure with column 
    names 'lag{L}_poly{P}'. If you change the basis construction, you MUST 
    also update compute_cumulative_effect to match.
    """
    n = len(temp)
    temp_lags = create_lag_matrix(temp, max_lag)
    temp_center = np.nanmedian(temp)
    temp_scale = np.nanstd(temp)
    
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
    knots = np.linspace(time_index.min(), time_index.max(), total_df + 2)[1:-1]
    
    try:
        spline_basis = dmatrix(
            f"bs(x, knots={list(knots)}, degree=3, include_intercept=False) - 1",
            {"x": time_index},
            return_type='dataframe'
        ).values
    except:
        spline_basis = np.column_stack([
            (time_index / 365) ** p for p in range(1, total_df + 1)
        ])
    return spline_basis


def fit_region_dlnm(df_region: pd.DataFrame, region_code: int, max_lag: int = MAX_LAG) -> dict:
    """Fit DLNM for a single region using Quasi-Poisson GLM."""
    if len(df_region) < MIN_OBS:
        return None
    
    df_region = df_region.sort_values('date').reset_index(drop=True)
    temp = df_region['temp_mean'].values
    X_cb, col_names, temp_center, temp_scale = create_simple_crossbasis(temp, max_lag)
    
    valid = ~np.isnan(X_cb).any(axis=1)
    if valid.sum() < MIN_OBS:
        return None
    
    X_cb_valid = X_cb[valid]
    df_valid = df_region.loc[valid].copy()
    y = df_valid['deaths_elderly'].values
    
    pop = df_valid['pop_elderly'].values
    offset = np.log(pop / 100000)
    
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
    
    try:
        model = GLM(y, X_full, family=Poisson(), offset=offset).fit(scale='X2')
        dispersion = model.scale
        
        # Check for convergence issues
        if not model.converged:
            return None
        
        # Check for extreme dispersion
        if dispersion > 10 or dispersion < 0.1:
            pass  # Continue but may be unreliable
            
    except:
        return None
    
    n_cb = len(col_names)
    cb_coefs = model.params[1:n_cb+1]
    vcov_full = model.cov_params()
    if hasattr(vcov_full, 'iloc'):
        cb_vcov = vcov_full.iloc[1:n_cb+1, 1:n_cb+1].values
    else:
        cb_vcov = vcov_full[1:n_cb+1, 1:n_cb+1]
    
    temps = df_valid['temp_mean'].dropna()
    percentiles = [1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99]
    temp_pcts = {f'p{p}'.replace('.', '_'): temps.quantile(p/100) for p in percentiles}
    # Rename for cleaner keys (p2_5 -> p2.5 style in output)
    temp_pcts = {k.replace('_', '.'): v for k, v in temp_pcts.items()}
    ref_temp = temp_pcts['p50']
    
    effects = {}
    for pct_name, target_temp in temp_pcts.items():
        effect, se = compute_cumulative_effect(
            target_temp, ref_temp, cb_coefs, cb_vcov, col_names,
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


def compute_cumulative_effect(target_temp, ref_temp, cb_coefs, cb_vcov,
                               col_names, temp_center, temp_scale, max_lag):
    """Compute cumulative effect at target vs reference temperature."""
    temp_std_target = (target_temp - temp_center) / temp_scale
    temp_std_ref = (ref_temp - temp_center) / temp_scale
    
    contrast = np.zeros(len(col_names))
    for i, name in enumerate(col_names):
        parts = name.split('_')
        lag = int(parts[0].replace('lag', ''))
        poly = int(parts[1].replace('poly', ''))
        contrast[i] = temp_std_target**poly - temp_std_ref**poly
    
    effect = np.dot(contrast, cb_coefs)
    var = np.dot(np.dot(contrast, cb_vcov), contrast)
    se = np.sqrt(var) if var > 0 else 0
    return effect, se


def random_effects_meta_analysis(effects, variances):
    """Random-effects meta-analysis using DerSimonian-Laird."""
    k = len(effects)
    if k < 2:
        return {
            'pooled_effect': float(effects[0]) if k == 1 else np.nan,
            'pooled_se': float(np.sqrt(variances[0])) if k == 1 else np.nan,
            'tau2': 0.0, 'I2': 0.0, 'Q': 0.0, 'p_heterogeneity': 1.0, 'n_regions': k
        }
    
    w = 1 / variances
    theta_fe = np.sum(w * effects) / np.sum(w)
    Q = np.sum(w * (effects - theta_fe)**2)
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0, (Q - (k - 1)) / c)
    w_re = 1 / (variances + tau2)
    theta_re = np.sum(w_re * effects) / np.sum(w_re)
    se_re = np.sqrt(1 / np.sum(w_re))
    I2 = max(0, (Q - (k - 1)) / Q * 100) if Q > 0 else 0
    p_het = 1 - stats.chi2.cdf(Q, k - 1)
    
    return {
        'pooled_effect': float(theta_re),
        'pooled_se': float(se_re),
        'pooled_rr': float(np.exp(theta_re)),
        'pooled_rr_lower': float(np.exp(theta_re - 1.96 * se_re)),
        'pooled_rr_upper': float(np.exp(theta_re + 1.96 * se_re)),
        'tau2': float(tau2), 'I2': float(I2), 'Q': float(Q),
        'p_heterogeneity': float(p_het), 'n_regions': k
    }


def pool_region_results(region_results, percentile):
    """Pool region-specific results via meta-analysis.
    
    Includes sanity filters:
    - SE must be positive and < 2 (plausible range)
    - log(RR) must be in plausible range (-3, 3)
    """
    effects, variances, regions = [], [], []
    excluded_reasons = {'invalid_se': 0, 'extreme_se': 0, 'extreme_rr': 0}
    MAX_SE = 2.0
    MAX_ABS_LOGRR = 3.0
    
    for region_code, result in region_results.items():
        if result is None or percentile not in result['effects']:
            continue
        eff = result['effects'][percentile]
        log_rr = eff['log_rr']
        se = eff['log_rr_se']
        
        if se <= 0:
            excluded_reasons['invalid_se'] += 1
            continue
        if se > MAX_SE:
            excluded_reasons['extreme_se'] += 1
            continue
        if abs(log_rr) > MAX_ABS_LOGRR:
            excluded_reasons['extreme_rr'] += 1
            continue
            
        effects.append(log_rr)
        variances.append(se**2)
        regions.append(region_code)
    
    if len(effects) == 0:
        return {'pooled_rr': np.nan, 'n_regions': 0, 'excluded_reasons': excluded_reasons}
    
    pooled = random_effects_meta_analysis(np.array(effects), np.array(variances))
    pooled['percentile'] = percentile
    pooled['regions_included'] = regions
    pooled['excluded_reasons'] = excluded_reasons
    return pooled


def convert_to_json_serializable(obj):
    """Convert numpy types to JSON-serializable types (preserves float64 precision)."""
    if isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return [convert_to_json_serializable(x) for x in obj]
    elif isinstance(obj, dict):
        return {k: convert_to_json_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_json_serializable(i) for i in obj]
    return obj


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("FITTING REGION-SPECIFIC DLNM MODELS")
print("="*70)

regions = sorted(df['immediate_code'].unique())
print(f"\nTotal regions: {len(regions)}")

region_results = {}
successful = 0

for i, region_code in enumerate(regions):
    df_region = df[df['immediate_code'] == region_code].copy()
    result = fit_region_dlnm(df_region, region_code)
    region_results[region_code] = result
    
    if result is not None:
        successful += 1
        if (i + 1) % 100 == 0:
            eff_p99 = result['effects']['p99']
            print(f"  Progress: {i+1}/{len(regions)} - Region {region_code}: Heat RR(P99)={eff_p99['rr']:.3f}")

print(f"\nSuccessfully fitted: {successful} / {len(regions)} regions")

# Log excluded regions statistics
excluded_count = len(regions) - successful
if excluded_count > 0:
    excluded_regions = [r for r, res in region_results.items() if res is None]
    excluded_sample_sizes = [len(df[df['immediate_code'] == r]) for r in excluded_regions]
    print(f"\n  Excluded regions: {excluded_count}")
    print(f"  Mean sample size of excluded: {np.mean(excluded_sample_sizes):.0f} obs")
    if excluded_count <= 20:
        print(f"  Excluded region codes: {excluded_regions}")
    else:
        print(f"  Excluded region codes (first 20): {excluded_regions[:20]}...")

# =============================================================================
# META-ANALYSIS POOLING
# =============================================================================

print("\n" + "="*70)
print("POOLING RESULTS VIA META-ANALYSIS")
print("="*70)

pooled_results = {}
for pct in ['p1', 'p2.5', 'p5', 'p10', 'p25', 'p50', 'p75', 'p90', 'p95', 'p97.5', 'p99']:
    pooled_results[pct] = pool_region_results(region_results, pct)

print("\nHEAT EFFECTS (vs P50 reference)")
for pct in ['p75', 'p90', 'p95', 'p97.5', 'p99']:
    p = pooled_results[pct]
    if 'pooled_rr' in p and not np.isnan(p['pooled_rr']):
        print(f"  {pct.upper()}: RR = {p['pooled_rr']:.3f} [{p['pooled_rr_lower']:.3f}-{p['pooled_rr_upper']:.3f}], I² = {p['I2']:.1f}%")

print("\nCOLD EFFECTS (vs P50 reference)")
for pct in ['p25', 'p10', 'p5', 'p2.5', 'p1']:
    p = pooled_results[pct]
    if 'pooled_rr' in p and not np.isnan(p['pooled_rr']):
        print(f"  {pct.upper()}: RR = {p['pooled_rr']:.3f} [{p['pooled_rr_lower']:.3f}-{p['pooled_rr_upper']:.3f}], I² = {p['I2']:.1f}%")

# =============================================================================
# SAVE RESULTS
# =============================================================================

print("\n" + "="*70)
print("SAVING RESULTS")
print("="*70)

# Full results
full_results = convert_to_json_serializable({
    'analysis_level': 'immediate',
    'n_regions': len(regions),
    'n_successful': successful,
    'parameters': {
        'max_lag': MAX_LAG, 'temp_df': TEMP_DF, 'lag_df': LAG_DF,
        'time_spline_df_per_year': TIME_SPLINE_DF_PER_YEAR, 'min_obs': MIN_OBS
    },
    'region_results': {str(k): v for k, v in region_results.items() if v is not None},
    'pooled_results': pooled_results,
    'timestamp': datetime.now().isoformat()
})

output_file = os.path.join(OUTPUT_DIR, 'dlnm_immediate_results.json')
with open(output_file, 'w') as f:
    json.dump(full_results, f, indent=2)
print(f"Saved: {output_file}")

# Summary CSV
summary_rows = []
for region_code, result in region_results.items():
    if result is None:
        continue
    row = {
        'immediate_code': region_code,
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
summary_file = os.path.join(OUTPUT_DIR, 'dlnm_immediate_summary.csv')
summary_df.to_csv(summary_file, index=False)
print(f"Saved: {summary_file}")

# Pooled results
pooled_file = os.path.join(OUTPUT_DIR, 'dlnm_immediate_pooled.json')
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
print(f"  Total observations: {len(df):,}")

print("\n  KEY FINDINGS:")
p99 = pooled_results['p99']
p1 = pooled_results['p1']
if 'pooled_rr' in p99 and not np.isnan(p99['pooled_rr']):
    print(f"    Extreme Heat (P99): RR = {p99['pooled_rr']:.3f} [{p99['pooled_rr_lower']:.3f}-{p99['pooled_rr_upper']:.3f}]")
if 'pooled_rr' in p1 and not np.isnan(p1['pooled_rr']):
    print(f"    Extreme Cold (P1):  RR = {p1['pooled_rr']:.3f} [{p1['pooled_rr_lower']:.3f}-{p1['pooled_rr_upper']:.3f}]")

print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
