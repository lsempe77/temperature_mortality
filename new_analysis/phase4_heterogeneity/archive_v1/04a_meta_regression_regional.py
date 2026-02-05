"""
04a: META-REGRESSION (REGIONAL DLNM)
====================================
Fits region-specific DLNMs and uses meta-regression to explain heterogeneity
in heat/cold effects across Brazil's intermediate regions.

Predictors:
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

import json
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.genmod.generalized_linear_model import GLM
from datetime import datetime
from pathlib import Path
from scipy import stats

# Paths
BASE_DIR = Path(__file__).parent.parent
PHASE0_RESULTS = BASE_DIR / 'phase0_data_prep' / 'results'
PHASE1_RESULTS = BASE_DIR / 'phase1_core_model' / 'results'
OUTPUT_DIR = Path(__file__).parent / 'results'
OUTPUT_DIR.mkdir(exist_ok=True)

MAX_LAG = 21
POLY_DEGREE = 3

print("="*70)
print("04a: META-REGRESSION (REGIONAL DLNM)")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# ============================================================================
# 1. LOAD DATA
# ============================================================================
print("\n[1] Loading data...")

# Temperature data
temp_df = pd.read_parquet(PHASE0_RESULTS / 'era5_intermediate_daily.parquet')
temp_df['date'] = pd.to_datetime(temp_df['date'])
print(f"  Temperature records: {len(temp_df):,}")

# Mortality data
mort_df = pd.read_parquet(PHASE0_RESULTS / 'mortality_regional_daily_elderly.parquet')
mort_df['date'] = pd.to_datetime(mort_df['date'])
print(f"  Mortality records: {len(mort_df):,}")

# Regional covariates
covariates = pd.read_csv(PHASE0_RESULTS / 'regional_covariates.csv')
print(f"  Regions with covariates: {len(covariates)}")
print(f"  Covariates: {covariates.columns.tolist()}")

# Merge data
df = mort_df.merge(temp_df, on=['date', 'region_code'], how='inner')
print(f"  Merged dataset: {len(df):,} observations")

# Add time variables
df['year'] = df['date'].dt.year
df['month'] = df['date'].dt.month
df['dow'] = df['date'].dt.dayofweek
df['doy'] = df['date'].dt.dayofyear

# Calculate temperature lags for each region
print("\n[2] Creating lagged temperatures...")
for lag in range(1, MAX_LAG + 1):
    df[f'temp_lag{lag}'] = df.groupby('region_code')['temp_mean'].shift(lag)

# Drop rows with missing lags
df = df.dropna(subset=[f'temp_lag{MAX_LAG}'])
print(f"  After lag creation: {len(df):,} observations")

# Global temperature distribution (for consistent percentiles - including P2.5/P97.5)
temp_p1 = df['temp_mean'].quantile(0.01)
temp_p2_5 = df['temp_mean'].quantile(0.025)
temp_p5 = df['temp_mean'].quantile(0.05)
temp_p50 = df['temp_mean'].quantile(0.50)
temp_p95 = df['temp_mean'].quantile(0.95)
temp_p97_5 = df['temp_mean'].quantile(0.975)
temp_p99 = df['temp_mean'].quantile(0.99)
temp_mean = df['temp_mean'].mean()
temp_std = df['temp_mean'].std()

print(f"  Temperature distribution:")
print(f"    P1: {temp_p1:.1f}°C, P2.5: {temp_p2_5:.1f}°C, P50: {temp_p50:.1f}°C, P97.5: {temp_p97_5:.1f}°C, P99: {temp_p99:.1f}°C")

# ============================================================================
# 2. REGION-SPECIFIC DLNM FITS
# ============================================================================
print("\n[3] Fitting region-specific DLNMs...")

def create_cb(df_sub, t_mean, t_std):
    """Create cross-basis matrix for DLNM."""
    n = len(df_sub)
    temp = df_sub['temp_mean'].values.astype(float)
    temp_z = (temp - t_mean) / t_std
    
    # Lagged temperatures
    temp_lagged = np.zeros((n, MAX_LAG + 1))
    temp_lagged[:, 0] = temp_z
    for lag in range(1, MAX_LAG + 1):
        temp_lagged[:, lag] = (df_sub[f'temp_lag{lag}'].values - t_mean) / t_std
    
    # Polynomial basis for temperature
    temp_basis = []
    for d in range(1, POLY_DEGREE + 1):
        temp_basis.append(temp_lagged ** d)
    
    # Lag basis (natural cubic spline-like with 3 knots)
    lag_knots = np.array([3, 7, 14])
    lag_idx = np.arange(MAX_LAG + 1).reshape(1, -1)
    
    # Build cross-basis
    cb_cols = []
    for d in range(1, POLY_DEGREE + 1):
        cb_cols.append(temp_basis[d-1].sum(axis=1))  # Cumulative effect
        for k in lag_knots:
            lag_weight = np.maximum(0, 1 - np.abs(lag_idx - k) / 4.0)
            weighted = (temp_basis[d-1] * lag_weight).sum(axis=1)
            cb_cols.append(weighted)
    
    return np.column_stack(cb_cols)

def fit_region_dlnm(df_region, t_mean, t_std):
    """Fit DLNM for a single region and extract heat/cold RRs."""
    try:
        n = len(df_region)
        if n < 500:  # Need enough data
            return None
        
        # Create cross-basis
        cb = create_cb(df_region, t_mean, t_std)
        
        # Control variables
        month_dummies = pd.get_dummies(df_region['month'], prefix='month', drop_first=True)
        dow_dummies = pd.get_dummies(df_region['dow'], prefix='dow', drop_first=True)
        year_dummies = pd.get_dummies(df_region['year'], prefix='year', drop_first=True)
        
        # Time trend
        time_idx = np.arange(n) / 365.25
        time_controls = np.column_stack([
            np.sin(2 * np.pi * time_idx),
            np.cos(2 * np.pi * time_idx),
            np.sin(4 * np.pi * time_idx),
            np.cos(4 * np.pi * time_idx)
        ])
        
        # Design matrix
        X = np.column_stack([
            np.ones(n),
            cb,
            month_dummies.values,
            dow_dummies.values,
            year_dummies.values,
            time_controls
        ])
        
        y = df_region['deaths_elderly'].values.astype(float)
        
        # Fit model
        model = GLM(y, X, family=sm.families.Poisson())
        result = model.fit()
        
        # Extract temperature effect coefficients (skip intercept)
        n_cb = cb.shape[1]
        temp_coefs = result.params[1:n_cb+1]
        temp_se = result.bse[1:n_cb+1]
        temp_cov = result.cov_params()[1:n_cb+1, 1:n_cb+1]
        
        # Calculate RR at heat (P99) and cold (P1) vs reference (P50)
        def get_rr_at_temp(temp_val, ref_val, coefs, cov):
            """Calculate RR at temp vs reference."""
            # Simplified: use cumulative effect approximation
            temp_z = (temp_val - t_mean) / t_std
            ref_z = (ref_val - t_mean) / t_std
            
            # Effect basis (just polynomial terms summed over lags)
            effect = 0
            for d in range(1, POLY_DEGREE + 1):
                idx = (d - 1) * 4  # Index into coefs
                effect += coefs[idx] * (temp_z**d - ref_z**d) * (MAX_LAG + 1)
            
            rr = np.exp(effect)
            
            # Variance (simplified)
            var = 0
            for d in range(1, POLY_DEGREE + 1):
                idx = (d - 1) * 4
                var += (temp_z**d - ref_z**d)**2 * (MAX_LAG + 1)**2 * cov[idx, idx]
            
            se_log = np.sqrt(max(var, 1e-10))
            rr_lo = np.exp(np.log(rr) - 1.96 * se_log)
            rr_hi = np.exp(np.log(rr) + 1.96 * se_log)
            
            return float(rr), float(rr_lo), float(rr_hi), float(se_log)
        
        # Calculate RR at heat (P97.5, P99) and cold (P2.5, P1) vs reference (P50)
        heat_rr_99, heat_lo_99, heat_hi_99, heat_se_99 = get_rr_at_temp(temp_p99, temp_p50, temp_coefs, temp_cov)
        heat_rr_97_5, heat_lo_97_5, heat_hi_97_5, heat_se_97_5 = get_rr_at_temp(temp_p97_5, temp_p50, temp_coefs, temp_cov)
        cold_rr_1, cold_lo_1, cold_hi_1, cold_se_1 = get_rr_at_temp(temp_p1, temp_p50, temp_coefs, temp_cov)
        cold_rr_2_5, cold_lo_2_5, cold_hi_2_5, cold_se_2_5 = get_rr_at_temp(temp_p2_5, temp_p50, temp_coefs, temp_cov)
        
        return {
            'n_obs': n,
            'mean_deaths': float(df_region['deaths_elderly'].mean()),
            # Primary thresholds (P97.5/P2.5)
            'heat_rr_97_5': heat_rr_97_5,
            'heat_rr_97_5_lo': heat_lo_97_5,
            'heat_rr_97_5_hi': heat_hi_97_5,
            'heat_log_se_97_5': heat_se_97_5,
            'cold_rr_2_5': cold_rr_2_5,
            'cold_rr_2_5_lo': cold_lo_2_5,
            'cold_rr_2_5_hi': cold_hi_2_5,
            'cold_log_se_2_5': cold_se_2_5,
            # Comparison thresholds (P99/P1)
            'heat_rr': heat_rr_99,
            'heat_rr_lo': heat_lo_99,
            'heat_rr_hi': heat_hi_99,
            'heat_log_se': heat_se_99,
            'cold_rr': cold_rr_1,
            'cold_rr_lo': cold_lo_1,
            'cold_rr_hi': cold_hi_1,
            'cold_log_se': cold_se_1
        }
        
    except Exception as e:
        return None

# Fit each region
region_results = {}
regions = df['region_code'].unique()
successful = 0

for i, region in enumerate(regions):
    if (i + 1) % 20 == 0:
        print(f"  Processing region {i+1}/{len(regions)}...")
    
    df_region = df[df['region_code'] == region].sort_values('date').reset_index(drop=True)
    result = fit_region_dlnm(df_region, temp_mean, temp_std)
    
    if result is not None:
        region_results[int(region)] = result
        successful += 1

print(f"  Successfully fit: {successful}/{len(regions)} regions")

# ============================================================================
# 3. META-REGRESSION
# ============================================================================
print("\n[4] Meta-regression analysis...")

# Prepare data for meta-regression
meta_data = []
for region_code, result in region_results.items():
    cov_row = covariates[covariates['region_code'] == region_code]
    if len(cov_row) == 0:
        continue
    
    cov_row = cov_row.iloc[0]
    meta_data.append({
        'region_code': region_code,
        'region_name': cov_row.get('region_name', ''),
        'heat_log_rr': np.log(result['heat_rr']),
        'heat_se': result['heat_log_se'],
        'cold_log_rr': np.log(result['cold_rr']),
        'cold_se': result['cold_log_se'],
        'n_obs': result['n_obs'],
        'ac_pct': cov_row.get('ac_pct', np.nan),
        'hdi': cov_row.get('hdi', np.nan),
        'gdp_per_capita': cov_row.get('gdp_per_capita_brl', np.nan),
        'urban_pct': cov_row.get('urban_pct', np.nan),
        'elderly_pct': cov_row.get('elderly_pct', np.nan),
        'hospital_beds': cov_row.get('hospital_beds_per_1000', np.nan),
        'physicians': cov_row.get('physicians_per_1000', np.nan),
        'mean_temp': cov_row.get('mean_temp_annual', np.nan),
        'macro_region': cov_row.get('macro_region', '')
    })

meta_df = pd.DataFrame(meta_data)
meta_df = meta_df.dropna()
print(f"  Regions in meta-regression: {len(meta_df)}")

# Standardize predictors
predictors = ['ac_pct', 'hdi', 'urban_pct', 'elderly_pct', 'mean_temp', 'hospital_beds']
for pred in predictors:
    meta_df[f'{pred}_z'] = (meta_df[pred] - meta_df[pred].mean()) / meta_df[pred].std()

# Run meta-regressions (weighted by inverse variance)
meta_results = {'heat': {}, 'cold': {}}

for outcome in ['heat', 'cold']:
    print(f"\n  === {outcome.upper()} META-REGRESSION ===")
    
    y = meta_df[f'{outcome}_log_rr'].values
    weights = 1 / (meta_df[f'{outcome}_se'].values ** 2 + 1e-6)
    
    # Random effects pooled estimate (simple weighted mean)
    pooled_mean = np.average(y, weights=weights)
    pooled_var = 1 / np.sum(weights)
    pooled_rr = np.exp(pooled_mean)
    pooled_rr_lo = np.exp(pooled_mean - 1.96 * np.sqrt(pooled_var))
    pooled_rr_hi = np.exp(pooled_mean + 1.96 * np.sqrt(pooled_var))
    
    print(f"    Pooled RR: {pooled_rr:.3f} [{pooled_rr_lo:.3f}-{pooled_rr_hi:.3f}]")
    
    # I² heterogeneity
    Q = np.sum(weights * (y - pooled_mean)**2)
    df_q = len(y) - 1
    I2 = max(0, (Q - df_q) / Q * 100) if Q > 0 else 0
    print(f"    I² heterogeneity: {I2:.1f}%")
    
    meta_results[outcome]['pooled'] = {
        'rr': float(pooled_rr),
        'rr_lo': float(pooled_rr_lo),
        'rr_hi': float(pooled_rr_hi),
        'i_squared': float(I2),
        'n_regions': len(y)
    }
    
    # Meta-regression for each predictor
    meta_results[outcome]['predictors'] = {}
    
    for pred in predictors:
        X = sm.add_constant(meta_df[f'{pred}_z'].values)
        
        try:
            model = sm.WLS(y, X, weights=weights)
            result = model.fit()
            
            coef = result.params[1]
            se = result.bse[1]
            pval = result.pvalues[1]
            
            # Effect per 1 SD increase
            rr_change = np.exp(coef)
            rr_lo = np.exp(coef - 1.96 * se)
            rr_hi = np.exp(coef + 1.96 * se)
            
            print(f"    {pred}: RR ratio = {rr_change:.3f} [{rr_lo:.3f}-{rr_hi:.3f}], p = {pval:.3f}")
            
            meta_results[outcome]['predictors'][pred] = {
                'coef': float(coef),
                'se': float(se),
                'rr_ratio_per_sd': float(rr_change),
                'rr_lo': float(rr_lo),
                'rr_hi': float(rr_hi),
                'pvalue': float(pval),
                'significant': bool(pval < 0.05)
            }
        except Exception as e:
            print(f"    {pred}: Failed - {str(e)[:50]}")

# ============================================================================
# 4. MULTIVARIATE META-REGRESSION
# ============================================================================
print("\n[5] Multivariate meta-regression...")

for outcome in ['heat', 'cold']:
    print(f"\n  === {outcome.upper()} MULTIVARIATE ===")
    
    y = meta_df[f'{outcome}_log_rr'].values
    weights = 1 / (meta_df[f'{outcome}_se'].values ** 2 + 1e-6)
    
    # Key predictors
    X_vars = ['ac_pct_z', 'mean_temp_z', 'urban_pct_z']
    X = sm.add_constant(meta_df[X_vars].values)
    
    try:
        model = sm.WLS(y, X, weights=weights)
        result = model.fit()
        
        print(f"    R² = {result.rsquared:.3f}")
        
        meta_results[outcome]['multivariate'] = {
            'r_squared': float(result.rsquared),
            'predictors': {}
        }
        
        for i, var in enumerate(X_vars):
            coef = result.params[i+1]
            se = result.bse[i+1]
            pval = result.pvalues[i+1]
            
            print(f"    {var}: β = {coef:.4f}, p = {pval:.3f}")
            
            meta_results[outcome]['multivariate']['predictors'][var] = {
                'coef': float(coef),
                'se': float(se),
                'pvalue': float(pval)
            }
    except Exception as e:
        print(f"    Failed: {e}")

# ============================================================================
# 5. REGIONAL EFFECT SUMMARY BY MACRO-REGION
# ============================================================================
print("\n[6] Effects by macro-region...")

macro_summary = meta_df.groupby('macro_region').agg({
    'heat_log_rr': 'mean',
    'cold_log_rr': 'mean',
    'ac_pct': 'mean',
    'mean_temp': 'mean',
    'region_code': 'count'
}).rename(columns={'region_code': 'n_regions'})

macro_summary['heat_rr'] = np.exp(macro_summary['heat_log_rr'])
macro_summary['cold_rr'] = np.exp(macro_summary['cold_log_rr'])

print("\n  Macro-region summary:")
print(f"  {'Region':<15} {'N':<5} {'Heat RR':<10} {'Cold RR':<10} {'AC %':<8} {'Mean T':<8}")
for idx, row in macro_summary.iterrows():
    print(f"  {idx:<15} {row['n_regions']:<5.0f} {row['heat_rr']:<10.3f} {row['cold_rr']:<10.3f} {row['ac_pct']:<8.1f} {row['mean_temp']:<8.1f}")

meta_results['macro_region_summary'] = macro_summary.reset_index().to_dict('records')

# ============================================================================
# 6. SAVE RESULTS
# ============================================================================
print("\n[7] Saving results...")

# Full results
results = {
    'analysis': 'Meta-regression of regional heat/cold effects',
    'timestamp': datetime.now().isoformat(),
    'n_regions_fitted': len(region_results),
    'n_regions_meta': len(meta_df),
    'temperature_reference': {
        'p1': float(temp_p1),
        'p50': float(temp_p50),
        'p99': float(temp_p99)
    },
    'region_specific_results': {str(k): v for k, v in region_results.items()},
    'meta_regression': meta_results
}

with open(OUTPUT_DIR / 'meta_regression_regional.json', 'w') as f:
    json.dump(results, f, indent=2)

# Summary table
summary_rows = []
for outcome in ['heat', 'cold']:
    pooled = meta_results[outcome]['pooled']
    summary_rows.append({
        'outcome': outcome,
        'metric': 'pooled_rr',
        'value': pooled['rr'],
        'ci_lower': pooled['rr_lo'],
        'ci_upper': pooled['rr_hi'],
        'n_regions': pooled['n_regions'],
        'i_squared': pooled['i_squared']
    })
    
    for pred, res in meta_results[outcome]['predictors'].items():
        summary_rows.append({
            'outcome': outcome,
            'metric': f'{pred}_effect',
            'value': res['rr_ratio_per_sd'],
            'ci_lower': res['rr_lo'],
            'ci_upper': res['rr_hi'],
            'pvalue': res['pvalue'],
            'significant': res['significant']
        })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(OUTPUT_DIR / 'meta_regression_summary.csv', index=False)

print(f"\n  Saved: meta_regression_regional.json")
print(f"  Saved: meta_regression_summary.csv")

# ============================================================================
# 7. KEY FINDINGS
# ============================================================================
print("\n" + "="*70)
print("KEY FINDINGS: META-REGRESSION")
print("="*70)

print(f"\nPooled Effects (133 regions):")
print(f"  Heat (P99 vs P50): RR = {meta_results['heat']['pooled']['rr']:.3f} "
      f"[{meta_results['heat']['pooled']['rr_lo']:.3f}-{meta_results['heat']['pooled']['rr_hi']:.3f}]")
print(f"  Cold (P1 vs P50):  RR = {meta_results['cold']['pooled']['rr']:.3f} "
      f"[{meta_results['cold']['pooled']['rr_lo']:.3f}-{meta_results['cold']['pooled']['rr_hi']:.3f}]")

print(f"\nHeterogeneity:")
print(f"  Heat I²: {meta_results['heat']['pooled']['i_squared']:.1f}%")
print(f"  Cold I²: {meta_results['cold']['pooled']['i_squared']:.1f}%")

print(f"\nSignificant predictors of heat effect variation:")
for pred, res in meta_results['heat']['predictors'].items():
    if res['significant']:
        direction = 'higher' if res['rr_ratio_per_sd'] > 1 else 'lower'
        print(f"  - {pred}: {direction} heat RR in regions with higher {pred}")

print(f"\nSignificant predictors of cold effect variation:")
for pred, res in meta_results['cold']['predictors'].items():
    if res['significant']:
        direction = 'higher' if res['rr_ratio_per_sd'] > 1 else 'lower'
        print(f"  - {pred}: {direction} cold RR in regions with higher {pred}")

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
