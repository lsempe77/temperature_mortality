"""
04c: SEX STRATIFICATION (REGIONAL DLNM)
=======================================
Tests effect modification by sex: Male vs Female elderly

This script:
1. Processes raw SIM mortality by sex
2. Aggregates to intermediate regions
3. Fits DLNM for each sex
4. Tests for heterogeneity between sexes
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
INPUT_DATA_DIR = BASE_DIR.parent / 'Input_data'
PHASE0_RESULTS = BASE_DIR / 'phase0_data_prep' / 'results'
OUTPUT_DIR = Path(__file__).parent / 'results'
OUTPUT_DIR.mkdir(exist_ok=True)

MAX_LAG = 21
POLY_DEGREE = 3

print("="*70)
print("04c: SEX STRATIFICATION (REGIONAL DLNM)")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# ============================================================================
# 1. LOAD MUNICIPALITY TO REGION MAPPING
# ============================================================================
print("\n[1] Loading municipality-region mapping...")

try:
    # Use municipal covariates file which has the mapping
    muni_df = pd.read_parquet(PHASE0_RESULTS / 'municipal_covariates.parquet')
    muni_to_region = dict(zip(
        muni_df['code_muni'].astype(str).str[:6],  # 6-digit code
        muni_df['region_code'].astype(int)  # Intermediate region code
    ))
    print(f"  Loaded mapping for {len(muni_to_region)} municipalities")
except Exception as e:
    print(f"  Could not load mapping: {e}")
    muni_to_region = None

# ============================================================================
# 2. PROCESS MORTALITY DATA BY SEX
# ============================================================================
print("\n[2] Processing mortality data by sex...")

def parse_age(idade):
    """Parse SIM age format."""
    if pd.isna(idade):
        return None
    try:
        idade = str(int(float(idade)))
        if len(idade) < 2:
            return None
        unit = int(idade[0])
        value = int(idade[1:])
        if unit == 4:
            return value
        elif unit == 5:
            return 100 + value
        return None
    except:
        return None

def parse_sex(sexo):
    """Parse SIM sex: 1=Male, 2=Female."""
    try:
        s = int(float(sexo))
        if s == 1:
            return 'Male'
        elif s == 2:
            return 'Female'
        return None
    except:
        return None

def parse_date(dtobito):
    """Parse SIM date format DDMMYYYY."""
    try:
        if pd.isna(dtobito):
            return None
        s = str(int(float(dtobito))).zfill(8)
        day, month, year = int(s[:2]), int(s[2:4]), int(s[4:])
        return pd.Timestamp(year=year, month=month, day=day)
    except:
        return None

# Process each year
sex_mortality = []
years = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24]

for year_suffix in years:
    fname = f'DO{year_suffix}OPEN.csv'
    fpath = INPUT_DATA_DIR / fname
    
    if not fpath.exists():
        continue
    
    full_year = 2000 + year_suffix if year_suffix < 50 else 1900 + year_suffix
    print(f"  Processing {fname} ({full_year})...")
    
    try:
        df = pd.read_csv(fpath, sep=';', low_memory=False, encoding='latin-1',
                        usecols=['DTOBITO', 'IDADE', 'SEXO', 'CODMUNRES'])
        
        df['date'] = df['DTOBITO'].apply(parse_date)
        df['age'] = df['IDADE'].apply(parse_age)
        df['sex'] = df['SEXO'].apply(parse_sex)
        
        # Only elderly (60+)
        df = df[(df['age'].notna()) & (df['age'] >= 60)].copy()
        df = df[df['sex'].notna()].copy()
        df = df[df['date'].notna()].copy()
        
        # Map to region
        df['muni_code'] = df['CODMUNRES'].astype(str).str[:6]
        
        if muni_to_region is not None:
            df['region_code'] = df['muni_code'].map(muni_to_region)
        else:
            df['region_code'] = df['muni_code'].str[:2].astype(int) * 100
        
        df = df.dropna(subset=['region_code'])
        df['region_code'] = df['region_code'].astype(int)
        
        sex_mortality.append(df[['date', 'region_code', 'sex']])
        
    except Exception as e:
        print(f"    Error: {e}")

if sex_mortality:
    mort_all = pd.concat(sex_mortality, ignore_index=True)
    print(f"\n  Total elderly deaths: {len(mort_all):,}")
    print(f"  Male: {(mort_all['sex']=='Male').sum():,}")
    print(f"  Female: {(mort_all['sex']=='Female').sum():,}")
else:
    raise ValueError("No mortality data loaded!")

# Aggregate to region-day-sex
print("\n[3] Aggregating to region-day-sex...")

mort_agg = mort_all.groupby(['date', 'region_code', 'sex']).size().reset_index(name='deaths')
print(f"  Aggregated records: {len(mort_agg):,}")

# ============================================================================
# 3. LOAD TEMPERATURE DATA
# ============================================================================
print("\n[4] Loading temperature data...")

temp_df = pd.read_parquet(PHASE0_RESULTS / 'era5_intermediate_daily.parquet')
temp_df['date'] = pd.to_datetime(temp_df['date'])
print(f"  Temperature records: {len(temp_df):,}")

# ============================================================================
# 4. PREPARE DATA FOR EACH SEX
# ============================================================================
print("\n[5] Preparing analysis datasets...")

sexes = ['Male', 'Female']
datasets = {}

for sex in sexes:
    mort_sex = mort_agg[mort_agg['sex'] == sex].copy()
    mort_sex = mort_sex.rename(columns={'deaths': f'deaths_{sex}'})
    
    # Merge with temperature
    df_sex = mort_sex.merge(temp_df, on=['date', 'region_code'], how='inner')
    
    # Add time variables
    df_sex['year'] = df_sex['date'].dt.year
    df_sex['month'] = df_sex['date'].dt.month
    df_sex['dow'] = df_sex['date'].dt.dayofweek
    
    # Create lags
    for lag in range(1, MAX_LAG + 1):
        df_sex[f'temp_lag{lag}'] = df_sex.groupby('region_code')['temp_mean'].shift(lag)
    
    df_sex = df_sex.dropna(subset=[f'temp_lag{MAX_LAG}'])
    
    datasets[sex] = df_sex
    print(f"  {sex}: {len(df_sex):,} obs, mean deaths = {df_sex[f'deaths_{sex}'].mean():.1f}")

# Global temperature percentiles (including P2.5/P97.5)
all_temps = temp_df['temp_mean']
temp_p1 = all_temps.quantile(0.01)
temp_p2_5 = all_temps.quantile(0.025)
temp_p50 = all_temps.quantile(0.50)
temp_p97_5 = all_temps.quantile(0.975)
temp_p99 = all_temps.quantile(0.99)
temp_mean = all_temps.mean()
temp_std = all_temps.std()

print(f"\n  Temperature: P1={temp_p1:.1f}, P2.5={temp_p2_5:.1f}, P50={temp_p50:.1f}, P97.5={temp_p97_5:.1f}, P99={temp_p99:.1f}")

# ============================================================================
# 5. FIT DLNM FOR EACH SEX
# ============================================================================
print("\n[6] Fitting DLNM for each sex...")

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
    
    cb = create_crossbasis(df, t_mean, t_std)
    n_cb = cb.shape[1]
    
    month_dummies = pd.get_dummies(df['month'], prefix='month', drop_first=True)
    dow_dummies = pd.get_dummies(df['dow'], prefix='dow', drop_first=True)
    year_dummies = pd.get_dummies(df['year'], prefix='year', drop_first=True)
    region_dummies = pd.get_dummies(df['region_code'], prefix='region', drop_first=True)
    
    time_idx = np.arange(n) / 365.25
    time_controls = np.column_stack([
        np.sin(2 * np.pi * time_idx),
        np.cos(2 * np.pi * time_idx),
        np.sin(4 * np.pi * time_idx),
        np.cos(4 * np.pi * time_idx)
    ])
    
    X = np.column_stack([
        np.ones(n),
        cb,
        month_dummies.values,
        dow_dummies.values,
        year_dummies.values,
        region_dummies.values,
        time_controls
    ])
    
    y = df[deaths_col].values.astype(float)
    
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
        'n_obs': n,
        'overdispersion': float(result.scale)
    }

results = {}
for sex in sexes:
    print(f"\n  Fitting {sex}...")
    df_sex = datasets[sex]
    deaths_col = f'deaths_{sex}'
    
    try:
        res = fit_dlnm(df_sex, deaths_col, temp_mean, temp_std)
        results[sex] = res
        print(f"    Heat RR: {res['heat_rr']:.3f} [{res['heat_rr_lo']:.3f}-{res['heat_rr_hi']:.3f}]")
        print(f"    Cold RR: {res['cold_rr']:.3f} [{res['cold_rr_lo']:.3f}-{res['cold_rr_hi']:.3f}]")
    except Exception as e:
        print(f"    Error: {e}")
        results[sex] = None

# ============================================================================
# 6. TEST FOR HETEROGENEITY
# ============================================================================
print("\n[7] Testing heterogeneity between sexes...")

if all(results.get(s) for s in sexes):
    # Wald test for difference
    male_log = np.log(results['Male']['heat_rr'])
    female_log = np.log(results['Female']['heat_rr'])
    male_se = results['Male']['heat_log_se']
    female_se = results['Female']['heat_log_se']
    
    diff = male_log - female_log
    se_diff = np.sqrt(male_se**2 + female_se**2)
    z = diff / se_diff
    p_diff = 2 * (1 - stats.norm.cdf(abs(z)))
    
    ratio_rr = np.exp(diff)
    ratio_lo = np.exp(diff - 1.96 * se_diff)
    ratio_hi = np.exp(diff + 1.96 * se_diff)
    
    print(f"\n  Heat effect: Male vs Female")
    print(f"    RR ratio (Male/Female): {ratio_rr:.3f} [{ratio_lo:.3f}-{ratio_hi:.3f}]")
    print(f"    Z = {z:.2f}, p = {p_diff:.3f}")
    
    heterogeneity = {
        'heat': {
            'male_female_ratio': float(ratio_rr),
            'ratio_lo': float(ratio_lo),
            'ratio_hi': float(ratio_hi),
            'z_statistic': float(z),
            'pvalue': float(p_diff)
        },
        'cold': {}
    }
    
    # Same for cold
    male_log_c = np.log(results['Male']['cold_rr'])
    female_log_c = np.log(results['Female']['cold_rr'])
    male_se_c = results['Male']['cold_log_se']
    female_se_c = results['Female']['cold_log_se']
    
    diff_c = male_log_c - female_log_c
    se_diff_c = np.sqrt(male_se_c**2 + female_se_c**2)
    z_c = diff_c / se_diff_c
    p_diff_c = 2 * (1 - stats.norm.cdf(abs(z_c)))
    
    ratio_rr_c = np.exp(diff_c)
    ratio_lo_c = np.exp(diff_c - 1.96 * se_diff_c)
    ratio_hi_c = np.exp(diff_c + 1.96 * se_diff_c)
    
    print(f"\n  Cold effect: Male vs Female")
    print(f"    RR ratio (Male/Female): {ratio_rr_c:.3f} [{ratio_lo_c:.3f}-{ratio_hi_c:.3f}]")
    print(f"    Z = {z_c:.2f}, p = {p_diff_c:.3f}")
    
    heterogeneity['cold'] = {
        'male_female_ratio': float(ratio_rr_c),
        'ratio_lo': float(ratio_lo_c),
        'ratio_hi': float(ratio_hi_c),
        'z_statistic': float(z_c),
        'pvalue': float(p_diff_c)
    }
else:
    heterogeneity = None

# ============================================================================
# 7. ATTRIBUTABLE BURDEN BY SEX
# ============================================================================
print("\n[8] Calculating attributable burden by sex...")

burden_results = {}

for sex in sexes:
    if results.get(sex) is None:
        burden_results[sex] = None
        continue
    
    print(f"\n  {sex}:")
    df_sex = datasets[sex]
    deaths_col = f'deaths_{sex}'
    
    # Total deaths for this sex
    total_deaths = df_sex[deaths_col].sum()
    
    # PRIMARY THRESHOLDS (P2.5/P97.5)
    heat_rr_97_5 = results[sex]['heat_rr_97_5']
    cold_rr_2_5 = results[sex]['cold_rr_2_5']
    
    # Days above P97.5 (heat days)
    heat_days = df_sex[df_sex['temp_mean'] >= temp_p97_5]
    heat_deaths_97_5 = heat_days[deaths_col].sum()
    
    # Days below P2.5 (cold days)
    cold_days = df_sex[df_sex['temp_mean'] <= temp_p2_5]
    cold_deaths_2_5 = cold_days[deaths_col].sum()
    
    # AF = (RR - 1) / RR
    heat_af_97_5 = (heat_rr_97_5 - 1) / heat_rr_97_5 if heat_rr_97_5 > 1 else 0
    cold_af_2_5 = (cold_rr_2_5 - 1) / cold_rr_2_5 if cold_rr_2_5 > 1 else 0
    
    # AN = AF * deaths on extreme days
    heat_an_97_5 = heat_af_97_5 * heat_deaths_97_5
    cold_an_2_5 = cold_af_2_5 * cold_deaths_2_5
    
    # COMPARISON THRESHOLDS (P1/P99)
    heat_rr_99 = results[sex]['heat_rr']
    cold_rr_1 = results[sex]['cold_rr']
    
    # Days above P99 (extreme heat days)
    heat_days_99 = df_sex[df_sex['temp_mean'] >= temp_p99]
    heat_deaths_99 = heat_days_99[deaths_col].sum()
    
    # Days below P1 (extreme cold days)
    cold_days_1 = df_sex[df_sex['temp_mean'] <= temp_p1]
    cold_deaths_1 = cold_days_1[deaths_col].sum()
    
    heat_af_99 = (heat_rr_99 - 1) / heat_rr_99 if heat_rr_99 > 1 else 0
    cold_af_1 = (cold_rr_1 - 1) / cold_rr_1 if cold_rr_1 > 1 else 0
    
    heat_an_99 = heat_af_99 * heat_deaths_99
    cold_an_1 = cold_af_1 * cold_deaths_1
    
    burden_results[sex] = {
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
# 8. SAVE RESULTS
# ============================================================================
print("\n[9] Saving results...")

output = {
    'analysis': 'Sex stratification DLNM',
    'timestamp': datetime.now().isoformat(),
    'sexes': sexes,
    'temperature_reference': {
        'p1': float(temp_p1),
        'p2_5': float(temp_p2_5),
        'p50': float(temp_p50),
        'p97_5': float(temp_p97_5),
        'p99': float(temp_p99)
    },
    'results': results,
    'burden': burden_results,
    'heterogeneity': heterogeneity
}

with open(OUTPUT_DIR / 'sex_stratification_regional.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)

# Summary CSV
summary_rows = []
for sex, res in results.items():
    if res is not None:
        summary_rows.append({
            'sex': sex,
            'n_obs': res['n_obs'],
            'heat_rr': res['heat_rr'],
            'heat_rr_lo': res['heat_rr_lo'],
            'heat_rr_hi': res['heat_rr_hi'],
            'cold_rr': res['cold_rr'],
            'cold_rr_lo': res['cold_rr_lo'],
            'cold_rr_hi': res['cold_rr_hi']
        })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(OUTPUT_DIR / 'sex_stratification_summary.csv', index=False)

print(f"\n  Saved: sex_stratification_regional.json")
print(f"  Saved: sex_stratification_summary.csv")

# ============================================================================
# 9. KEY FINDINGS
# ============================================================================
print("\n" + "="*70)
print("KEY FINDINGS: SEX STRATIFICATION")
print("="*70)

print("\nRelative Risks by Sex:")
print(f"\n{'Sex':<12} {'Heat RR P97.5 [95% CI]':<28} {'Cold RR P2.5 [95% CI]':<28}")
print("-"*68)
for sex in sexes:
    if results.get(sex):
        r = results[sex]
        heat_str = f"{r['heat_rr_97_5']:.3f} [{r['heat_rr_97_5_lo']:.3f}-{r['heat_rr_97_5_hi']:.3f}]"
        cold_str = f"{r['cold_rr_2_5']:.3f} [{r['cold_rr_2_5_lo']:.3f}-{r['cold_rr_2_5_hi']:.3f}]"
        print(f"{sex:<12} {heat_str:<28} {cold_str:<28}")

print("\nAttributable Burden by Sex (PRIMARY P2.5/P97.5):")
print(f"\n{'Sex':<12} {'Total Deaths':<15} {'Heat AN':<12} {'Cold AN':<12} {'Total AN':<12}")
print("-"*65)
for sex in sexes:
    if burden_results.get(sex):
        b = burden_results[sex]
        total_an = b['heat_an_97_5'] + b['cold_an_2_5']
        print(f"{sex:<12} {b['total_deaths']:>14,} {b['heat_an_97_5']:>11,.0f} {b['cold_an_2_5']:>11,.0f} {total_an:>11,.0f}")

# Sum across sexes
total_heat_an = sum(b['heat_an_97_5'] for b in burden_results.values() if b)
total_cold_an = sum(b['cold_an_2_5'] for b in burden_results.values() if b)
print("-"*65)
print(f"{'TOTAL':<12} {'':<15} {total_heat_an:>11,.0f} {total_cold_an:>11,.0f} {total_heat_an + total_cold_an:>11,.0f}")

if heterogeneity:
    print(f"\nMale vs Female comparison:")
    print(f"  Heat RR ratio: {heterogeneity['heat']['male_female_ratio']:.3f}, "
          f"p = {heterogeneity['heat']['pvalue']:.3f}")
    print(f"  Cold RR ratio: {heterogeneity['cold']['male_female_ratio']:.3f}, "
          f"p = {heterogeneity['cold']['pvalue']:.3f}")

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
