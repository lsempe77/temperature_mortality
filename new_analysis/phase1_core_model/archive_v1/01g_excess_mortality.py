"""
01g: EXCESS MORTALITY VALIDATION
=================================
Validates DLNM attributable burden using counterfactual excess mortality approach.

Method:
1. Fit baseline model predicting expected deaths (without temperature)
2. Calculate excess = observed - expected on extreme temperature days
3. Compare with DLNM attributable estimates

This provides a simple, intuitive validation of the DLNM burden estimates.

References:
- Fouillet et al. (2006) - Excess mortality approach for 2003 heatwave
- Gasparrini et al. (2022) - Comparison of burden estimation methods
"""

import warnings
warnings.filterwarnings('ignore')

import json
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
import statsmodels.api as sm
from statsmodels.genmod.generalized_linear_model import GLM

# =============================================================================
# PATHS
# =============================================================================

SCRIPT_DIR = Path(__file__).parent
PHASE0_RESULTS = SCRIPT_DIR.parent / 'phase0_data_prep' / 'results'
PHASE1_RESULTS = SCRIPT_DIR / 'results'
OUTPUT_DIR = PHASE1_RESULTS
OUTPUT_DIR.mkdir(exist_ok=True)

print("=" * 70)
print("01g: EXCESS MORTALITY VALIDATION")
print("=" * 70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print()

# =============================================================================
# LOAD DATA
# =============================================================================
print("[1] Loading data...")

# Load temperature data
df_temp = pd.read_parquet(PHASE0_RESULTS / 'era5_intermediate_daily.parquet')
df_temp['date'] = pd.to_datetime(df_temp['date'])

# Load mortality data  
df_mort = pd.read_parquet(PHASE0_RESULTS / 'mortality_regional_daily_elderly.parquet')
df_mort['date'] = pd.to_datetime(df_mort['date'])

# Merge
df = pd.merge(
    df_mort,
    df_temp[['date', 'region_code', 'temp_mean']],
    on=['date', 'region_code'],
    how='inner'
)

# Aggregate to national daily
df_national = df.groupby('date').agg({
    'deaths_elderly': 'sum',
    'temp_mean': 'mean'
}).reset_index()

df_national['year'] = df_national['date'].dt.year
df_national['month'] = df_national['date'].dt.month
df_national['dow'] = df_national['date'].dt.dayofweek
df_national['doy'] = df_national['date'].dt.dayofyear
df_national['time_idx'] = (df_national['date'] - df_national['date'].min()).dt.days

print(f"  National daily observations: {len(df_national):,}")
print(f"  Total deaths: {df_national['deaths_elderly'].sum():,}")
print(f"  Date range: {df_national['date'].min().date()} to {df_national['date'].max().date()}")

# Temperature percentiles
temp_p25 = df_national['temp_mean'].quantile(0.25)
temp_p50 = df_national['temp_mean'].quantile(0.50)
temp_p75 = df_national['temp_mean'].quantile(0.75)
temp_p95 = df_national['temp_mean'].quantile(0.95)
temp_p05 = df_national['temp_mean'].quantile(0.05)

print(f"\n  Temperature percentiles:")
print(f"    P5:  {temp_p05:.1f}°C")
print(f"    P25: {temp_p25:.1f}°C")
print(f"    P50: {temp_p50:.1f}°C")
print(f"    P75: {temp_p75:.1f}°C")
print(f"    P95: {temp_p95:.1f}°C")

# =============================================================================
# FIT BASELINE MODEL (WITHOUT TEMPERATURE)
# =============================================================================
print("\n[2] Fitting baseline model (expected deaths without temperature)...")

# Create design matrix - seasonal and trend controls only
dow_dummies = pd.get_dummies(df_national['dow'], prefix='dow', drop_first=True)
month_dummies = pd.get_dummies(df_national['month'], prefix='month', drop_first=True)

X_baseline = pd.concat([dow_dummies, month_dummies], axis=1)
X_baseline['time_trend'] = df_national['time_idx'] / 365.25
X_baseline['time_trend_sq'] = X_baseline['time_trend'] ** 2
X_baseline = sm.add_constant(X_baseline)
X_baseline = X_baseline.astype(float)

y = df_national['deaths_elderly'].values.astype(float)

# Fit Quasi-Poisson model
model_baseline = GLM(y, X_baseline, family=sm.families.Poisson()).fit(scale='X2')

# Predicted (expected) deaths
df_national['expected'] = model_baseline.predict(X_baseline)
df_national['excess'] = df_national['deaths_elderly'] - df_national['expected']

print(f"  Model dispersion: {model_baseline.scale:.2f}")
print(f"  Mean observed: {df_national['deaths_elderly'].mean():.0f}")
print(f"  Mean expected: {df_national['expected'].mean():.0f}")

# =============================================================================
# CALCULATE EXCESS ON EXTREME TEMPERATURE DAYS
# =============================================================================
print("\n[3] Calculating excess mortality on extreme temperature days...")

# Define extreme temperature periods
heat_mask = df_national['temp_mean'] > temp_p75
cold_mask = df_national['temp_mean'] < temp_p25
normal_mask = (df_national['temp_mean'] >= temp_p25) & (df_national['temp_mean'] <= temp_p75)

extreme_heat_mask = df_national['temp_mean'] > temp_p95
extreme_cold_mask = df_national['temp_mean'] < temp_p05

# Calculate excess for each period
heat_excess = df_national.loc[heat_mask, 'excess'].sum()
cold_excess = df_national.loc[cold_mask, 'excess'].sum()
normal_excess = df_national.loc[normal_mask, 'excess'].sum()

extreme_heat_excess = df_national.loc[extreme_heat_mask, 'excess'].sum()
extreme_cold_excess = df_national.loc[extreme_cold_mask, 'excess'].sum()

n_heat_days = heat_mask.sum()
n_cold_days = cold_mask.sum()
n_extreme_heat = extreme_heat_mask.sum()
n_extreme_cold = extreme_cold_mask.sum()

print(f"\n  Hot days (>P75): {n_heat_days} days")
print(f"    Excess deaths: {heat_excess:,.0f}")
print(f"    Mean excess/day: {heat_excess/n_heat_days:.1f}")

print(f"\n  Cold days (<P25): {n_cold_days} days")
print(f"    Excess deaths: {cold_excess:,.0f}")
print(f"    Mean excess/day: {cold_excess/n_cold_days:.1f}")

print(f"\n  Extreme heat (>P95): {n_extreme_heat} days")
print(f"    Excess deaths: {extreme_heat_excess:,.0f}")
print(f"    Mean excess/day: {extreme_heat_excess/n_extreme_heat:.1f}")

print(f"\n  Extreme cold (<P5): {n_extreme_cold} days")
print(f"    Excess deaths: {extreme_cold_excess:,.0f}")
print(f"    Mean excess/day: {extreme_cold_excess/n_extreme_cold:.1f}")

# =============================================================================
# COMPARE WITH DLNM BURDEN
# =============================================================================
print("\n[4] Comparing with DLNM attributable burden...")

# Load DLNM burden results
with open(PHASE1_RESULTS / 'attributable_burden_national.json', 'r') as f:
    dlnm_burden = json.load(f)

dlnm_heat_an = dlnm_burden['intermediate']['heat_an']
dlnm_cold_an = dlnm_burden['intermediate']['cold_an']

print(f"\n  DLNM Attributable (P75/P25 thresholds):")
print(f"    Heat AN: {dlnm_heat_an:,.0f}")
print(f"    Cold AN: {dlnm_cold_an:,.0f}")

print(f"\n  Excess Mortality (P75/P25 thresholds):")
print(f"    Heat excess: {heat_excess:,.0f}")
print(f"    Cold excess: {cold_excess:,.0f}")

# Calculate ratio
heat_ratio = heat_excess / dlnm_heat_an if dlnm_heat_an > 0 else np.nan
cold_ratio = cold_excess / dlnm_cold_an if dlnm_cold_an > 0 else np.nan

print(f"\n  Ratio (Excess / DLNM):")
print(f"    Heat: {heat_ratio:.2f}")
print(f"    Cold: {cold_ratio:.2f}")

# =============================================================================
# ANNUAL EXCESS MORTALITY
# =============================================================================
print("\n[5] Annual excess mortality by temperature category...")

annual_excess = df_national.groupby('year').apply(
    lambda x: pd.Series({
        'total_deaths': x['deaths_elderly'].sum(),
        'expected': x['expected'].sum(),
        'total_excess': x['excess'].sum(),
        'heat_excess': x.loc[x['temp_mean'] > temp_p75, 'excess'].sum(),
        'cold_excess': x.loc[x['temp_mean'] < temp_p25, 'excess'].sum(),
        'n_heat_days': (x['temp_mean'] > temp_p75).sum(),
        'n_cold_days': (x['temp_mean'] < temp_p25).sum()
    })
).reset_index()

print("\n  Annual Summary:")
print(annual_excess[['year', 'total_deaths', 'heat_excess', 'cold_excess']].to_string(index=False))

# =============================================================================
# SAVE RESULTS
# =============================================================================
print("\n[6] Saving results...")

def convert_for_json(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: convert_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_for_json(i) for i in obj]
    return obj

excess_results = {
    'method': 'Counterfactual excess mortality (observed - expected from baseline model)',
    'baseline_model': 'Quasi-Poisson with month, dow, time trend controls (no temperature)',
    'temperature_thresholds': {
        'p05': float(temp_p05),
        'p25': float(temp_p25),
        'p50': float(temp_p50),
        'p75': float(temp_p75),
        'p95': float(temp_p95)
    },
    'excess_mortality': {
        'heat_p75': {
            'n_days': int(n_heat_days),
            'excess_deaths': float(heat_excess),
            'mean_excess_per_day': float(heat_excess / n_heat_days)
        },
        'cold_p25': {
            'n_days': int(n_cold_days),
            'excess_deaths': float(cold_excess),
            'mean_excess_per_day': float(cold_excess / n_cold_days)
        },
        'extreme_heat_p95': {
            'n_days': int(n_extreme_heat),
            'excess_deaths': float(extreme_heat_excess),
            'mean_excess_per_day': float(extreme_heat_excess / n_extreme_heat)
        },
        'extreme_cold_p05': {
            'n_days': int(n_extreme_cold),
            'excess_deaths': float(extreme_cold_excess),
            'mean_excess_per_day': float(extreme_cold_excess / n_extreme_cold)
        }
    },
    'comparison_with_dlnm': {
        'dlnm_heat_an': dlnm_heat_an,
        'dlnm_cold_an': dlnm_cold_an,
        'excess_heat': float(heat_excess),
        'excess_cold': float(cold_excess),
        'ratio_heat': float(heat_ratio),
        'ratio_cold': float(cold_ratio)
    },
    'annual_summary': annual_excess.to_dict('records'),
    'timestamp': datetime.now().isoformat()
}

with open(OUTPUT_DIR / 'excess_mortality_validation.json', 'w') as f:
    json.dump(convert_for_json(excess_results), f, indent=2)

annual_excess.to_csv(OUTPUT_DIR / 'excess_mortality_annual.csv', index=False)

print(f"  Saved: excess_mortality_validation.json")
print(f"  Saved: excess_mortality_annual.csv")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("VALIDATION SUMMARY")
print("=" * 70)

print(f"""
Excess Mortality Approach vs DLNM Attributable Burden:

1. HEAT (days > P75):
   - Excess mortality: {heat_excess:,.0f} deaths
   - DLNM attributable: {dlnm_heat_an:,.0f} deaths
   - Ratio: {heat_ratio:.2f}

2. COLD (days < P25):
   - Excess mortality: {cold_excess:,.0f} deaths
   - DLNM attributable: {dlnm_cold_an:,.0f} deaths
   - Ratio: {cold_ratio:.2f}

3. INTERPRETATION:
   - Ratio > 1: Excess approach finds MORE deaths than DLNM
   - Ratio < 1: Excess approach finds FEWER deaths than DLNM
   - Ratio ≈ 1: Methods agree closely

   The excess mortality approach is simpler but doesn't account for:
   - Distributed lag effects (mortality spread over days)
   - Non-linear temperature-mortality relationship
   - Harvesting/displacement effects
   
   DLNM is more sophisticated and accounts for these factors.
   
4. COLD > HEAT PATTERN: {'✓ CONFIRMED' if cold_excess > heat_excess else '✗ NOT CONFIRMED'}
""")

print("=" * 70)
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)
