"""
04d: CAUSE-SPECIFIC STRATIFICATION (REGIONAL DLNM)
==================================================
Tests effect modification by cause of death:
- All-cause
- Cardiovascular (CVD)
- Respiratory
- Heat-related (external causes)

Uses pre-aggregated regional mortality data with cause stratification.
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
OUTPUT_DIR = Path(__file__).parent / 'results'
OUTPUT_DIR.mkdir(exist_ok=True)

MAX_LAG = 21
POLY_DEGREE = 3

print("="*70)
print("04d: CAUSE-SPECIFIC STRATIFICATION (REGIONAL DLNM)")
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

# Mortality data with cause stratification
mort_df = pd.read_parquet(PHASE0_RESULTS / 'mortality_regional_daily_elderly.parquet')
mort_df['date'] = pd.to_datetime(mort_df['date'])
print(f"  Mortality records: {len(mort_df):,}")
print(f"  Mortality columns: {mort_df.columns.tolist()}")

# Merge data
df = mort_df.merge(temp_df, on=['date', 'region_code'], how='inner')
print(f"  Merged dataset: {len(df):,} observations")

# Add time variables
df['year'] = df['date'].dt.year
df['month'] = df['date'].dt.month
df['dow'] = df['date'].dt.dayofweek

# Create lags
print("\n[2] Creating lagged temperatures...")
for lag in range(1, MAX_LAG + 1):
    df[f'temp_lag{lag}'] = df.groupby('region_code')['temp_mean'].shift(lag)

df = df.dropna(subset=[f'temp_lag{MAX_LAG}'])
print(f"  After lag creation: {len(df):,} observations")

# Calculate "other" causes (all-cause minus CVD and respiratory)
if 'deaths_elderly_resp' in df.columns and 'deaths_elderly_cvd' in df.columns:
    df['deaths_elderly_other'] = (df['deaths_elderly'] - 
                                   df['deaths_elderly_resp'] - 
                                   df['deaths_elderly_cvd'])
    df['deaths_elderly_other'] = df['deaths_elderly_other'].clip(lower=0)

# Temperature distribution (including P2.5/P97.5)
temp_p1 = df['temp_mean'].quantile(0.01)
temp_p2_5 = df['temp_mean'].quantile(0.025)
temp_p50 = df['temp_mean'].quantile(0.50)
temp_p97_5 = df['temp_mean'].quantile(0.975)
temp_p99 = df['temp_mean'].quantile(0.99)
temp_mean = df['temp_mean'].mean()
temp_std = df['temp_mean'].std()

print(f"\n  Temperature: P1={temp_p1:.1f}, P2.5={temp_p2_5:.1f}, P50={temp_p50:.1f}, P97.5={temp_p97_5:.1f}, P99={temp_p99:.1f}")

# ============================================================================
# 2. DEFINE CAUSES
# ============================================================================
print("\n[3] Defining cause categories...")

causes = {
    'all_cause': {
        'col': 'deaths_elderly',
        'label': 'All-cause'
    },
    'cardiovascular': {
        'col': 'deaths_elderly_cvd',
        'label': 'Cardiovascular (I00-I99)'
    },
    'respiratory': {
        'col': 'deaths_elderly_resp',
        'label': 'Respiratory (J00-J99)'
    }
}

# Check if heat-related causes exist
if 'deaths_elderly_heat' in df.columns:
    causes['heat_related'] = {
        'col': 'deaths_elderly_heat',
        'label': 'Heat-related (X30)'
    }

if 'deaths_elderly_other' in df.columns:
    causes['other'] = {
        'col': 'deaths_elderly_other',
        'label': 'Other causes'
    }

# Print death counts
for cause_key, cause_info in causes.items():
    col = cause_info['col']
    if col in df.columns:
        total = df[col].sum()
        mean_daily = df[col].mean()
        print(f"  {cause_info['label']}: {total:,.0f} total, {mean_daily:.1f} daily mean")

# ============================================================================
# 3. FIT DLNM FOR EACH CAUSE
# ============================================================================
print("\n[4] Fitting DLNM for each cause...")

def create_crossbasis(df, t_mean, t_std):
    """Create cross-basis for DLNM."""
    n = len(df)
    temp = df['temp_mean'].values.astype(float)
    temp_z = (temp - t_mean) / t_std
    
    temp_lagged = np.zeros((n, MAX_LAG + 1))
    temp_lagged[:, 0] = temp_z
    for lag in range(1, MAX_LAG + 1):
        temp_lagged[:, lag] = (df[f'temp_lag{lag}'].values - t_mean) / t_std
    
    temp_basis = []
    for d in range(1, POLY_DEGREE + 1):
        temp_basis.append(temp_lagged ** d)
    
    lag_knots = np.array([3, 7, 14])
    lag_idx = np.arange(MAX_LAG + 1).reshape(1, -1)
    
    cb_cols = []
    for d in range(1, POLY_DEGREE + 1):
        cb_cols.append(temp_basis[d-1].sum(axis=1))
        for k in lag_knots:
            lag_weight = np.maximum(0, 1 - np.abs(lag_idx - k) / 4.0)
            weighted = (temp_basis[d-1] * lag_weight).sum(axis=1)
            cb_cols.append(weighted)
    
    return np.column_stack(cb_cols)

def fit_dlnm(df, deaths_col, t_mean, t_std):
    """Fit DLNM and extract effects."""
    n = len(df)
    
    # Filter out rows with zero deaths for this cause (avoid log(0))
    df_valid = df[df[deaths_col] > 0].copy().reset_index(drop=True)
    n_valid = len(df_valid)
    
    if n_valid < 1000:
        return None
    
    cb = create_crossbasis(df_valid, t_mean, t_std)
    n_cb = cb.shape[1]
    
    month_dummies = pd.get_dummies(df_valid['month'], prefix='month', drop_first=True)
    dow_dummies = pd.get_dummies(df_valid['dow'], prefix='dow', drop_first=True)
    year_dummies = pd.get_dummies(df_valid['year'], prefix='year', drop_first=True)
    region_dummies = pd.get_dummies(df_valid['region_code'], prefix='region', drop_first=True)
    
    time_idx = np.arange(n_valid) / 365.25
    time_controls = np.column_stack([
        np.sin(2 * np.pi * time_idx),
        np.cos(2 * np.pi * time_idx),
        np.sin(4 * np.pi * time_idx),
        np.cos(4 * np.pi * time_idx)
    ])
    
    X = np.column_stack([
        np.ones(n_valid),
        cb,
        month_dummies.values,
        dow_dummies.values,
        year_dummies.values,
        region_dummies.values,
        time_controls
    ])
    
    y = df_valid[deaths_col].values.astype(float)
    
    model = GLM(y, X, family=sm.families.Poisson())
    result = model.fit(scale='X2')
    
    temp_coefs = result.params[1:n_cb+1]
    temp_cov = result.cov_params()[1:n_cb+1, 1:n_cb+1]
    
    def calc_rr(temp_val, ref_val):
        temp_z = (temp_val - t_mean) / t_std
        ref_z = (ref_val - t_mean) / t_std
        
        effect = 0
        for d in range(1, POLY_DEGREE + 1):
            idx = (d - 1) * 4
            effect += temp_coefs[idx] * (temp_z**d - ref_z**d) * (MAX_LAG + 1)
        
        rr = np.exp(effect)
        
        var = 0
        for d in range(1, POLY_DEGREE + 1):
            idx = (d - 1) * 4
            var += (temp_z**d - ref_z**d)**2 * (MAX_LAG + 1)**2 * temp_cov[idx, idx]
        
        se_log = np.sqrt(max(var, 1e-10))
        rr_lo = np.exp(np.log(rr) - 1.96 * se_log)
        rr_hi = np.exp(np.log(rr) + 1.96 * se_log)
        
        return float(rr), float(rr_lo), float(rr_hi), float(se_log)
    
    # Calculate RRs at all thresholds
    heat_rr_99 = calc_rr(temp_p99, temp_p50)
    heat_rr_97_5 = calc_rr(temp_p97_5, temp_p50)
    cold_rr_1 = calc_rr(temp_p1, temp_p50)
    cold_rr_2_5 = calc_rr(temp_p2_5, temp_p50)
    
    return {
        # Primary thresholds (P97.5/P2.5)
        'heat_rr_97_5': heat_rr_97_5[0],
        'heat_rr_97_5_lo': heat_rr_97_5[1],
        'heat_rr_97_5_hi': heat_rr_97_5[2],
        'heat_log_se_97_5': heat_rr_97_5[3],
        'cold_rr_2_5': cold_rr_2_5[0],
        'cold_rr_2_5_lo': cold_rr_2_5[1],
        'cold_rr_2_5_hi': cold_rr_2_5[2],
        'cold_log_se_2_5': cold_rr_2_5[3],
        # Comparison thresholds (P99/P1)
        'heat_rr': heat_rr_99[0],
        'heat_rr_lo': heat_rr_99[1],
        'heat_rr_hi': heat_rr_99[2],
        'heat_log_se': heat_rr_99[3],
        'cold_rr': cold_rr_1[0],
        'cold_rr_lo': cold_rr_1[1],
        'cold_rr_hi': cold_rr_1[2],
        'cold_log_se': cold_rr_1[3],
        'n_obs': n_valid,
        'n_original': n,
        'overdispersion': float(result.scale)
    }

results = {}
for cause_key, cause_info in causes.items():
    col = cause_info['col']
    label = cause_info['label']
    
    print(f"\n  Fitting {label}...")
    
    if col not in df.columns:
        print(f"    Column {col} not found, skipping")
        continue
    
    try:
        res = fit_dlnm(df, col, temp_mean, temp_std)
        if res is not None:
            results[cause_key] = res
            results[cause_key]['label'] = label
            print(f"    Heat RR: {res['heat_rr']:.3f} [{res['heat_rr_lo']:.3f}-{res['heat_rr_hi']:.3f}]")
            print(f"    Cold RR: {res['cold_rr']:.3f} [{res['cold_rr_lo']:.3f}-{res['cold_rr_hi']:.3f}]")
        else:
            print(f"    Insufficient data")
    except Exception as e:
        print(f"    Error: {e}")

# ============================================================================
# 4. TEST FOR HETEROGENEITY
# ============================================================================
print("\n[5] Testing heterogeneity across causes...")

# Compare CVD vs Respiratory
if 'cardiovascular' in results and 'respiratory' in results:
    cvd = results['cardiovascular']
    resp = results['respiratory']
    
    # Heat
    diff_heat = np.log(cvd['heat_rr']) - np.log(resp['heat_rr'])
    se_diff_heat = np.sqrt(cvd['heat_log_se']**2 + resp['heat_log_se']**2)
    z_heat = diff_heat / se_diff_heat
    p_heat = 2 * (1 - stats.norm.cdf(abs(z_heat)))
    
    print(f"\n  CVD vs Respiratory (Heat):")
    print(f"    Log-RR difference: {diff_heat:.3f}")
    print(f"    Z = {z_heat:.2f}, p = {p_heat:.3f}")
    
    # Cold
    diff_cold = np.log(cvd['cold_rr']) - np.log(resp['cold_rr'])
    se_diff_cold = np.sqrt(cvd['cold_log_se']**2 + resp['cold_log_se']**2)
    z_cold = diff_cold / se_diff_cold
    p_cold = 2 * (1 - stats.norm.cdf(abs(z_cold)))
    
    print(f"\n  CVD vs Respiratory (Cold):")
    print(f"    Log-RR difference: {diff_cold:.3f}")
    print(f"    Z = {z_cold:.2f}, p = {p_cold:.3f}")
    
    heterogeneity = {
        'cvd_vs_resp': {
            'heat': {
                'diff_log_rr': float(diff_heat),
                'z': float(z_heat),
                'pvalue': float(p_heat)
            },
            'cold': {
                'diff_log_rr': float(diff_cold),
                'z': float(z_cold),
                'pvalue': float(p_cold)
            }
        }
    }
else:
    heterogeneity = None

# ============================================================================
# 5. ATTRIBUTABLE BURDEN BY CAUSE
# ============================================================================
print("\n[6] Calculating attributable burden by cause...")

burden_results = {}
for cause_key, cause_info in causes.items():
    if cause_key not in results:
        continue
    
    res = results[cause_key]
    col = cause_info['col']
    label = cause_info['label']
    
    print(f"\n  {label}:")
    
    # Total deaths
    total_deaths = df[col].sum()
    
    # PRIMARY THRESHOLDS (P2.5/P97.5)
    heat_rr_97_5 = res['heat_rr_97_5']
    cold_rr_2_5 = res['cold_rr_2_5']
    
    # Days above P97.5 (heat days)
    heat_days = df[df['temp_mean'] >= temp_p97_5]
    heat_deaths_97_5 = heat_days[col].sum()
    
    # Days below P2.5 (cold days)
    cold_days = df[df['temp_mean'] <= temp_p2_5]
    cold_deaths_2_5 = cold_days[col].sum()
    
    # AF = (RR - 1) / RR
    heat_af_97_5 = (heat_rr_97_5 - 1) / heat_rr_97_5 if heat_rr_97_5 > 1 else 0
    cold_af_2_5 = (cold_rr_2_5 - 1) / cold_rr_2_5 if cold_rr_2_5 > 1 else 0
    
    # AN = AF * deaths on extreme days
    heat_an_97_5 = heat_af_97_5 * heat_deaths_97_5
    cold_an_2_5 = cold_af_2_5 * cold_deaths_2_5
    
    # COMPARISON THRESHOLDS (P1/P99)
    heat_rr_99 = res['heat_rr']
    cold_rr_1 = res['cold_rr']
    
    # Days above P99
    heat_days_99 = df[df['temp_mean'] >= temp_p99]
    heat_deaths_99 = heat_days_99[col].sum()
    
    # Days below P1
    cold_days_1 = df[df['temp_mean'] <= temp_p1]
    cold_deaths_1 = cold_days_1[col].sum()
    
    heat_af_99 = (heat_rr_99 - 1) / heat_rr_99 if heat_rr_99 > 1 else 0
    cold_af_1 = (cold_rr_1 - 1) / cold_rr_1 if cold_rr_1 > 1 else 0
    
    heat_an_99 = heat_af_99 * heat_deaths_99
    cold_an_1 = cold_af_1 * cold_deaths_1
    
    burden_results[cause_key] = {
        'label': label,
        'total_deaths': int(total_deaths),
        # PRIMARY (P97.5/P2.5)
        'heat_days_97_5': len(heat_days),
        'heat_deaths_97_5': int(heat_deaths_97_5),
        'heat_af_97_5': float(heat_af_97_5),
        'heat_an_97_5': float(heat_an_97_5),
        'cold_days_2_5': len(cold_days),
        'cold_deaths_2_5': int(cold_deaths_2_5),
        'cold_af_2_5': float(cold_af_2_5),
        'cold_an_2_5': float(cold_an_2_5),
        # COMPARISON (P99/P1)
        'heat_days_99': len(heat_days_99),
        'heat_deaths_99': int(heat_deaths_99),
        'heat_af_99': float(heat_af_99),
        'heat_an_99': float(heat_an_99),
        'cold_days_1': len(cold_days_1),
        'cold_deaths_1': int(cold_deaths_1),
        'cold_af_1': float(cold_af_1),
        'cold_an_1': float(cold_an_1)
    }
    
    print(f"    Total deaths: {total_deaths:,.0f}")
    print(f"    Heat burden (P97.5): AF={heat_af_97_5:.3f}, AN={heat_an_97_5:,.0f}")
    print(f"    Cold burden (P2.5):  AF={cold_af_2_5:.3f}, AN={cold_an_2_5:,.0f}")
    print(f"    Heat burden (P99):   AF={heat_af_99:.3f}, AN={heat_an_99:,.0f}")
    print(f"    Cold burden (P1):    AF={cold_af_1:.3f}, AN={cold_an_1:,.0f}")

# ============================================================================
# 6. SAVE RESULTS
# ============================================================================
print("\n[7] Saving results...")

output = {
    'analysis': 'Cause-specific stratification DLNM',
    'timestamp': datetime.now().isoformat(),
    'causes': list(causes.keys()),
    'temperature_reference': {
        'p1': float(temp_p1),
        'p2_5': float(temp_p2_5),
        'p50': float(temp_p50),
        'p97_5': float(temp_p97_5),
        'p99': float(temp_p99)
    },
    'results': results,
    'heterogeneity': heterogeneity,
    'burden': burden_results
}

with open(OUTPUT_DIR / 'cause_stratification_regional.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)

# Summary CSV
summary_rows = []
for cause_key, res in results.items():
    if res is not None:
        b = burden_results.get(cause_key, {})
        summary_rows.append({
            'cause': cause_key,
            'label': res.get('label', ''),
            'n_obs': res['n_obs'],
            'heat_rr_97_5': res['heat_rr_97_5'],
            'heat_rr_97_5_lo': res['heat_rr_97_5_lo'],
            'heat_rr_97_5_hi': res['heat_rr_97_5_hi'],
            'cold_rr_2_5': res['cold_rr_2_5'],
            'cold_rr_2_5_lo': res['cold_rr_2_5_lo'],
            'cold_rr_2_5_hi': res['cold_rr_2_5_hi'],
            'heat_rr': res['heat_rr'],
            'heat_rr_lo': res['heat_rr_lo'],
            'heat_rr_hi': res['heat_rr_hi'],
            'cold_rr': res['cold_rr'],
            'cold_rr_lo': res['cold_rr_lo'],
            'cold_rr_hi': res['cold_rr_hi'],
            'heat_af_97_5_pct': b.get('heat_af_97_5', 0) * 100,
            'cold_af_2_5_pct': b.get('cold_af_2_5', 0) * 100,
            'heat_an_97_5': b.get('heat_an_97_5', 0),
            'cold_an_2_5': b.get('cold_an_2_5', 0)
        })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(OUTPUT_DIR / 'cause_stratification_summary.csv', index=False)

print(f"\n  Saved: cause_stratification_regional.json")
print(f"  Saved: cause_stratification_summary.csv")

# ============================================================================
# 7. KEY FINDINGS
# ============================================================================
print("\n" + "="*70)
print("KEY FINDINGS: CAUSE-SPECIFIC STRATIFICATION")
print("="*70)

print("\nRelative Risks by Cause (PRIMARY P2.5/P97.5):")
print(f"\n{'Cause':<25} {'Heat RR P97.5 [95% CI]':<28} {'Cold RR P2.5 [95% CI]':<28}")
print("-"*81)
for cause_key, res in results.items():
    if res:
        label = res.get('label', cause_key)[:24]
        heat_str = f"{res['heat_rr_97_5']:.3f} [{res['heat_rr_97_5_lo']:.3f}-{res['heat_rr_97_5_hi']:.3f}]"
        cold_str = f"{res['cold_rr_2_5']:.3f} [{res['cold_rr_2_5_lo']:.3f}-{res['cold_rr_2_5_hi']:.3f}]"
        print(f"{label:<25} {heat_str:<28} {cold_str:<28}")

print("\nAttributable Burden by Cause (PRIMARY P2.5/P97.5):")
print(f"\n{'Cause':<25} {'Total Deaths':<15} {'Heat AN':<12} {'Cold AN':<12}")
print("-"*70)
for cause_key, b in burden_results.items():
    label = b.get('label', cause_key)[:24]
    print(f"{label:<25} {b['total_deaths']:>14,} {b['heat_an_97_5']:>11,.0f} {b['cold_an_2_5']:>11,.0f}")

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
