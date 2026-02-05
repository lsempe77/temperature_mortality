"""
04c: SEX STRATIFICATION (REGIONAL DLNM) - v2
=============================================
Tests effect modification by sex: Male vs Female elderly

v2 IMPROVEMENTS (Dec 2025):
---------------------------
1. Per-region DLNM fits then meta-analysis (not pooled model with dummies)
2. Uses natural spline cross-basis from utils
3. Population offset in all models
4. Full covariance delta method for CIs
5. Region-specific percentiles
6. Proper heterogeneity testing between sexes

Two-stage approach (per sex):
1. Fit DLNM in each region → extract heat/cold RRs
2. Meta-analyze across regions (DerSimonian-Laird)
3. Compare pooled effects between sexes
"""

import warnings
warnings.filterwarnings('ignore')

import sys
import json
import pandas as pd
import numpy as np
from datetime import datetime
from pathlib import Path
from scipy import stats
from typing import Dict, Any, Optional, Tuple
import argparse

sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.dlnm_module import (
    fit_region_dlnm,
    predict_cumulative_rr,
    find_mmt,
    meta_random_effects,
    convert_to_json_serializable,
    mvmeta_pool_coefficients,
    compute_pooled_rr_from_mvmeta,
)

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = Path(__file__).parent.parent
INPUT_DATA_DIR = BASE_DIR.parent / 'Input_data'
PHASE0_RESULTS = BASE_DIR / 'phase0_data_prep' / 'results'
OUTPUT_DIR = Path(__file__).parent / 'results'
OUTPUT_DIR.mkdir(exist_ok=True)

MAX_LAG = 21
VAR_DF = 4
LAG_DF = 4
MIN_OBS = 365

SEXES = ['Male', 'Female']

parser = argparse.ArgumentParser(description='Sex-stratified regional DLNM meta-analysis (v2)')
parser.add_argument('--level', type=str, default='intermediate',
                    choices=['intermediate', 'immediate'],
                    help='Spatial level to analyze')
args, _ = parser.parse_known_args()

level = args.level
if level == 'intermediate':
    era5_filename = 'era5_intermediate_daily.parquet'
    pop_covariates_filename = 'regional_covariates.csv'
    pop_region_col = 'region_code'
    level_label = 'INTERMEDIATE (133 regions)'
    output_suffix = ''
else:
    era5_filename = 'era5_immediate_daily.parquet'
    pop_covariates_filename = 'ses_immediate_covariates.csv'
    pop_region_col = 'immediate_code'
    level_label = 'IMMEDIATE (510 regions)'
    output_suffix = '_immediate'

print("="*70)
print("04c: SEX STRATIFICATION (REGIONAL DLNM) - v2")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Using: Per-region fits + meta-analysis (not pooled model)")
print(f"Level: {level} [{level_label}]")

# =============================================================================
# 1. LOAD MUNICIPALITY TO REGION MAPPING
# =============================================================================
print("\n[1] Loading municipality-region mapping...")

muni_to_region = None
try:
    if level == 'intermediate':
        muni_df = pd.read_parquet(PHASE0_RESULTS / 'municipal_covariates.parquet')
        muni_to_region = dict(zip(
            muni_df['code_muni'].astype(str).str[:6],
            muni_df['region_code'].astype(int)
        ))
    else:
        muni_df = pd.read_csv(PHASE0_RESULTS / 'municipality_to_all_regions_map.csv')
        muni_to_region = dict(zip(
            muni_df['code_muni'].astype(str).str[:6],
            muni_df['immediate_code'].astype(int)
        ))
    print(f"  Loaded mapping for {len(muni_to_region)} municipalities")
except Exception as e:
    print(f"  Could not load mapping: {e}")
    muni_to_region = None

# =============================================================================
# 2. PROCESS MORTALITY DATA BY SEX
# =============================================================================
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
    print(f"  Processing {fname} ({full_year})...", end='')
    
    try:
        # Auto-detect separator (DO24 uses comma, earlier years use semicolon)
        with open(fpath, 'r', encoding='latin-1') as fh:
            first_line = fh.readline()
        sep = ',' if ('",' in first_line or first_line.count(',') > first_line.count(';')) else ';'
        
        df = pd.read_csv(fpath, sep=sep, low_memory=False, encoding='latin-1',
                        usecols=['DTOBITO', 'IDADE', 'SEXO', 'CODMUNRES'])
        
        df['date'] = df['DTOBITO'].apply(parse_date)
        df['age'] = df['IDADE'].apply(parse_age)
        df['sex'] = df['SEXO'].apply(parse_sex)
        
        # Only elderly (60+)
        df = df[(df['age'].notna()) & (df['age'] >= 60)].copy()
        df = df[df['sex'].notna() & df['date'].notna()].copy()
        
        df['muni_code'] = df['CODMUNRES'].astype(str).str[:6]
        
        if muni_to_region is not None:
            df['region_code'] = df['muni_code'].map(muni_to_region)
        else:
            df['region_code'] = df['muni_code'].str[:2].astype(int) * 100
        
        df = df.dropna(subset=['region_code'])
        df['region_code'] = df['region_code'].astype(int)
        
        sex_mortality.append(df[['date', 'region_code', 'sex']])
        print(f" {len(df):,} elderly deaths")
        
    except Exception as e:
        print(f" Error: {e}")

if sex_mortality:
    mort_all = pd.concat(sex_mortality, ignore_index=True)
    print(f"\n  Total elderly deaths: {len(mort_all):,}")
    for sex in SEXES:
        n = (mort_all['sex'] == sex).sum()
        print(f"    {sex}: {n:,}")
else:
    raise ValueError("No mortality data loaded!")

# Aggregate to region-day-sex
print("\n[3] Aggregating to region-day-sex...")

mort_agg = mort_all.groupby(['date', 'region_code', 'sex']).size().reset_index(name='deaths')
print(f"  Aggregated records: {len(mort_agg):,}")

# =============================================================================
# 3. LOAD TEMPERATURE DATA
# =============================================================================
print("\n[4] Loading temperature data...")

temp_df = pd.read_parquet(PHASE0_RESULTS / era5_filename)
temp_df['date'] = pd.to_datetime(temp_df['date'])
print(f"  Temperature records: {len(temp_df):,}")

# =============================================================================
# 4. LOAD POPULATION BY SEX (approximate)
# =============================================================================
print("\n[5] Loading population data...")

try:
    covariates = pd.read_csv(PHASE0_RESULTS / pop_covariates_filename)
    pop_elderly_total = covariates.set_index(pop_region_col)['pop_elderly'].to_dict()

    # Approximate sex proportions (based on IBGE data - more females in elderly)
    sex_proportions = {'Male': 0.44, 'Female': 0.56}
    print("  Using approximate sex proportions for population offset")
except Exception:
    pop_elderly_total = {}
    sex_proportions = {'Male': 0.44, 'Female': 0.56}
    print("  WARNING: No population data found; using default offsets")

# =============================================================================
# 5. PREPARE DATA FOR EACH SEX
# =============================================================================
print("\n[6] Preparing analysis datasets...")

datasets = {}

for sex in SEXES:
    mort_sex = mort_agg[mort_agg['sex'] == sex].copy()
    mort_sex = mort_sex.rename(columns={'deaths': 'deaths_elderly'})
    
    # Merge with temperature
    df_sex = mort_sex.merge(temp_df, on=['date', 'region_code'], how='inner')
    
    # Add time variables
    df_sex['dow'] = df_sex['date'].dt.dayofweek
    df_sex['month'] = df_sex['date'].dt.month
    df_sex['year'] = df_sex['date'].dt.year
    
    # Add population (proportional)
    if pop_elderly_total:
        df_sex['pop_elderly'] = df_sex['region_code'].map(pop_elderly_total) * sex_proportions[sex]
        df_sex['pop_elderly'] = df_sex['pop_elderly'].fillna(10000)
    else:
        df_sex['pop_elderly'] = 10000
    
    datasets[sex] = df_sex
    print(f"  {sex}: {len(df_sex):,} obs, mean deaths = {df_sex['deaths_elderly'].mean():.2f}")

# =============================================================================
# 6. COMPUTE REGION-SPECIFIC PERCENTILES
# =============================================================================
print("\n[7] Computing region-specific temperature percentiles...")

all_regions = set()
for sex in SEXES:
    all_regions.update(datasets[sex]['region_code'].unique())

region_percentiles = {}
for region in all_regions:
    temps = temp_df.loc[temp_df['region_code'] == region, 'temp_mean']
    if len(temps) < 365:
        continue
    region_percentiles[region] = {
        'p1': float(temps.quantile(0.01)),
        'p2_5': float(temps.quantile(0.025)),
        'p50': float(temps.quantile(0.50)),
        'p97_5': float(temps.quantile(0.975)),
        'p99': float(temps.quantile(0.99)),
    }

print(f"  Computed percentiles for {len(region_percentiles)} regions")

# =============================================================================
# 7. FIT PER-REGION DLNM FOR EACH SEX
# =============================================================================
print("\n[8] Fitting per-region DLNMs for each sex...")

def fit_region_sex(df_region: pd.DataFrame, pctiles: Dict) -> Optional[Dict]:
    """
    Fit DLNM for a region and extract effects relative to MMT.
    
    Following Gasparrini methodology:
    - Find MMT from fitted curve
    - Compare P99 (heat) and P1 (cold) to MMT
    - Store coefficients for MVMeta pooling
    """
    try:
        fit_res = fit_region_dlnm(
            df_region=df_region,
            temp_col='temp_mean',
            deaths_col='deaths_elderly',
            pop_col='pop_elderly',
            max_lag=MAX_LAG,
            var_df=VAR_DF,
            lag_df=LAG_DF,
            family='quasi-poisson',
            min_obs=MIN_OBS,
        )
        
        if fit_res is None:
            return None
        
        # Find MMT (Gasparrini methodology)
        temp_range = np.linspace(pctiles['p1'], pctiles['p99'], 100)
        mmt = find_mmt(fit_res, temp_range, pctiles['p50'])
        mmt_percentile = 1 + 98 * (mmt - pctiles['p1']) / (pctiles['p99'] - pctiles['p1'])
        mmt_percentile = max(1, min(99, mmt_percentile))
        
        ref_temp = mmt  # Use MMT as reference
        
        # Heat (P99 vs MMT)
        heat = predict_cumulative_rr(fit_res, pctiles['p99'], ref_temp)
        # Cold (P1 vs MMT)
        cold = predict_cumulative_rr(fit_res, pctiles['p1'], ref_temp)
        
        # Extract cross-basis coefficients for MVMeta pooling
        cb_colnames = fit_res.get('cb_colnames', [])
        params = fit_res.get('params', pd.Series())
        cov = fit_res.get('cov', pd.DataFrame())
        
        # Get CB coefficients and vcov
        if len(cb_colnames) > 0 and all(c in params.index for c in cb_colnames):
            cb_coefs = params[cb_colnames].values
            cb_vcov = cov.loc[cb_colnames, cb_colnames].values
        else:
            cb_coefs = np.array([])
            cb_vcov = np.array([[]])
        
        return {
            'n_obs': fit_res['n_obs'],
            'mmt': float(mmt),
            'mmt_percentile': float(mmt_percentile),
            'heat_rr': heat[0],
            'heat_rr_lo': heat[1],
            'heat_rr_hi': heat[2],
            'heat_log_rr': heat[3],
            'heat_log_rr_se': heat[4] if len(heat) > 4 else np.nan,
            'cold_rr': cold[0],
            'cold_rr_lo': cold[1],
            'cold_rr_hi': cold[2],
            'cold_log_rr': cold[3],
            'cold_log_rr_se': cold[4] if len(cold) > 4 else np.nan,
            # For MVMeta pooling
            'cb_coefs': cb_coefs.tolist() if len(cb_coefs) > 0 else None,
            'cb_vcov': cb_vcov.tolist() if cb_vcov.size > 0 else None,
            'cb_meta': fit_res.get('cb_meta'),
            'temp_percentiles': pctiles,
        }
    except Exception:
        return None

sex_results = {}

for sex in SEXES:
    print(f"\n  === {sex} ===")
    df_sex = datasets[sex]
    regions = df_sex['region_code'].unique()
    
    region_results = {}
    successful = 0
    
    for i, region in enumerate(regions):
        if (i + 1) % 30 == 0:
            print(f"    Region {i+1}/{len(regions)}...")
        
        if region not in region_percentiles:
            continue
        
        df_region = df_sex[df_sex['region_code'] == region].copy()
        result = fit_region_sex(df_region, region_percentiles[region])
        
        if result is not None:
            region_results[int(region)] = result
            successful += 1
    
    sex_results[sex] = region_results
    print(f"    Fitted: {successful}/{len(regions)} regions")

# =============================================================================
# 8. META-ANALYZE WITHIN EACH SEX USING MVMETA
# =============================================================================
print("\n[9] Meta-analyzing results within each sex using MVMeta...")

pooled_results = {}

for sex in SEXES:
    print(f"\n  === {sex} ===")
    region_results = sex_results[sex]
    
    if len(region_results) < 5:
        print(f"    Insufficient regions: {len(region_results)}")
        continue
    
    # Collect coefficients for MVMeta pooling
    mvmeta_coefs = []
    mvmeta_vcovs = []
    mvmeta_pctiles = []
    valid_regions = []
    
    for region_code, res in region_results.items():
        if res.get('cb_coefs') is not None and res.get('cb_vcov') is not None:
            coefs = np.array(res['cb_coefs'])
            vcov = np.array(res['cb_vcov'])
            if coefs.shape[0] > 0 and vcov.shape[0] == coefs.shape[0]:
                mvmeta_coefs.append(coefs)
                mvmeta_vcovs.append(vcov)
                mvmeta_pctiles.append(res.get('temp_percentiles', {}))
                valid_regions.append(region_code)
    
    if len(mvmeta_coefs) < 5:
        print(f"    Insufficient valid coefficients for MVMeta: {len(mvmeta_coefs)}")
        continue
    
    print(f"    Pooling {len(mvmeta_coefs)} regions via MVMeta...")
    
    # Pool coefficients using MVMeta
    try:
        mvmeta_result = mvmeta_pool_coefficients(mvmeta_coefs, mvmeta_vcovs)
        
        if mvmeta_result is None:
            print(f"    MVMeta pooling failed, falling back to simple meta-analysis")
            # Fallback to simple meta-analysis
            heat_effects = []
            heat_variances = []
            cold_effects = []
            cold_variances = []
            for res in region_results.values():
                if res['heat_rr'] > 0 and res['heat_rr_lo'] > 0 and res['heat_rr_hi'] > 0:
                    log_rr = np.log(res['heat_rr'])
                    se = (np.log(res['heat_rr_hi']) - np.log(res['heat_rr_lo'])) / (2 * 1.96)
                    if se > 0.01 and se < 2.0:
                        heat_effects.append(log_rr)
                        heat_variances.append(se ** 2)
                if res['cold_rr'] > 0 and res['cold_rr_lo'] > 0 and res['cold_rr_hi'] > 0:
                    log_rr = np.log(res['cold_rr'])
                    se = (np.log(res['cold_rr_hi']) - np.log(res['cold_rr_lo'])) / (2 * 1.96)
                    if se > 0.01 and se < 2.0:
                        cold_effects.append(log_rr)
                        cold_variances.append(se ** 2)
            
            heat_meta = meta_random_effects(np.array(heat_effects), np.array(heat_variances))
            cold_meta = meta_random_effects(np.array(cold_effects), np.array(cold_variances))
            
            pooled_results[sex] = {
                'method': 'simple_meta',
                'heat': {
                    'rr': float(np.exp(heat_meta['effect'])),
                    'rr_lo': float(np.exp(heat_meta['effect'] - 1.96 * heat_meta['se'])),
                    'rr_hi': float(np.exp(heat_meta['effect'] + 1.96 * heat_meta['se'])),
                    'I2': float(heat_meta['I2']),
                    'k': int(heat_meta['k']),
                },
                'cold': {
                    'rr': float(np.exp(cold_meta['effect'])),
                    'rr_lo': float(np.exp(cold_meta['effect'] - 1.96 * cold_meta['se'])),
                    'rr_hi': float(np.exp(cold_meta['effect'] + 1.96 * cold_meta['se'])),
                    'I2': float(cold_meta['I2']),
                    'k': int(cold_meta['k']),
                },
                'n_regions': len(region_results),
            }
        else:
            # Compute pooled RRs from MVMeta result
            # Get representative percentiles
            all_p1 = [p.get('p1', 10) for p in mvmeta_pctiles if p]
            all_p50 = [p.get('p50', 20) for p in mvmeta_pctiles if p]
            all_p99 = [p.get('p99', 35) for p in mvmeta_pctiles if p]
            
            pooled_pctiles = {
                'p1': float(np.median(all_p1)) if all_p1 else 10.0,
                'p50': float(np.median(all_p50)) if all_p50 else 20.0,
                'p99': float(np.median(all_p99)) if all_p99 else 35.0,
            }
            
            # Get cb_meta from first valid result
            cb_meta = None
            for res in region_results.values():
                if res.get('cb_meta'):
                    cb_meta = res['cb_meta']
                    break
            
            if cb_meta is None:
                print(f"    Missing cb_meta, using defaults")
                cb_meta = {'var_df': VAR_DF, 'lag_df': LAG_DF, 'max_lag': MAX_LAG}
            
            rr_results = compute_pooled_rr_from_mvmeta(
                mvmeta_result, pooled_pctiles, cb_meta
            )
            
            pooled_results[sex] = {
                'method': 'mvmeta',
                'heat': {
                    'rr': rr_results['heat_rr'],
                    'rr_lo': rr_results['heat_rr_lo'],
                    'rr_hi': rr_results['heat_rr_hi'],
                    'k': len(mvmeta_coefs),
                },
                'cold': {
                    'rr': rr_results['cold_rr'],
                    'rr_lo': rr_results['cold_rr_lo'],
                    'rr_hi': rr_results['cold_rr_hi'],
                    'k': len(mvmeta_coefs),
                },
                'mmt': rr_results['mmt'],
                'mmt_percentile': rr_results['mmt_percentile'],
                'n_regions': len(region_results),
                'n_pooled': len(mvmeta_coefs),
                'pooled_percentiles': pooled_pctiles,
            }
            
    except Exception as e:
        print(f"    MVMeta error: {e}")
        continue
    
    print(f"    Heat RR: {pooled_results[sex]['heat']['rr']:.3f} "
          f"[{pooled_results[sex]['heat']['rr_lo']:.3f}-{pooled_results[sex]['heat']['rr_hi']:.3f}] "
          f"(k={pooled_results[sex]['heat']['k']})")
    print(f"    Cold RR: {pooled_results[sex]['cold']['rr']:.3f} "
          f"[{pooled_results[sex]['cold']['rr_lo']:.3f}-{pooled_results[sex]['cold']['rr_hi']:.3f}] "
          f"(k={pooled_results[sex]['cold']['k']})")

# =============================================================================
# 9. COMPARE MALE VS FEMALE
# =============================================================================
print("\n[10] Comparing Male vs Female...")

def get_log_rr_se(eff: Dict) -> Tuple[float, float]:
    """Extract log_rr and se from MVMeta or simple meta results."""
    rr = eff.get('rr', 1.0)
    rr_lo = eff.get('rr_lo', rr)
    rr_hi = eff.get('rr_hi', rr)
    if rr > 0 and rr_lo > 0 and rr_hi > 0:
        log_rr = np.log(rr)
        se = (np.log(rr_hi) - np.log(rr_lo)) / (2 * 1.96)
        return log_rr, se
    return 0.0, 1.0

comparisons = {}

if 'Male' in pooled_results and 'Female' in pooled_results:
    for effect_type in ['heat', 'cold']:
        male = pooled_results['Male'][effect_type]
        female = pooled_results['Female'][effect_type]
        
        # Get log RR and SE for comparison
        male_log_rr, male_se = get_log_rr_se(male)
        female_log_rr, female_se = get_log_rr_se(female)
        
        # Difference in log-RR (Male - Female)
        diff = male_log_rr - female_log_rr
        se_diff = np.sqrt(male_se**2 + female_se**2)
        z = diff / se_diff if se_diff > 0 else 0
        p = 2 * (1 - stats.norm.cdf(abs(z)))
        
        rr_ratio = np.exp(diff)
        rr_ratio_lo = np.exp(diff - 1.96 * se_diff)
        rr_ratio_hi = np.exp(diff + 1.96 * se_diff)
        
        comparisons[effect_type] = {
            'male_vs_female': {
                'rr_ratio': float(rr_ratio),
                'rr_ratio_lo': float(rr_ratio_lo),
                'rr_ratio_hi': float(rr_ratio_hi),
                'diff_log_rr': float(diff),
                'se_diff': float(se_diff),
                'z': float(z),
                'p': float(p),
                'significant': bool(p < 0.05),
            }
        }
        
        sig = "*" if p < 0.05 else ""
        print(f"\n  {effect_type.upper()}: Male vs Female")
        print(f"    Male RR: {male['rr']:.3f} [{male['rr_lo']:.3f}-{male['rr_hi']:.3f}]")
        print(f"    Female RR: {female['rr']:.3f} [{female['rr_lo']:.3f}-{female['rr_hi']:.3f}]")
        print(f"    RR ratio (M/F): {rr_ratio:.3f} [{rr_ratio_lo:.3f}-{rr_ratio_hi:.3f}], p = {p:.3f}{sig}")

# =============================================================================
# 10. PAIRED WITHIN-REGION TEST (more powerful)
# =============================================================================
print("\n[11] Paired within-region comparison...")

# For regions with both sexes, compute difference and meta-analyze
paired_diffs = {'heat': [], 'cold': []}
paired_vars = {'heat': [], 'cold': []}

common_regions = set(sex_results.get('Male', {}).keys()) & set(sex_results.get('Female', {}).keys())
print(f"  Regions with both sexes: {len(common_regions)}")

for region in common_regions:
    male = sex_results['Male'].get(region)
    female = sex_results['Female'].get(region)
    
    if male is None or female is None:
        continue
    
    for effect_type in ['heat', 'cold']:
        rr_key = f'{effect_type}_rr'
        lo_key = f'{effect_type}_rr_lo'
        hi_key = f'{effect_type}_rr_hi'
        
        if male[rr_key] > 0 and female[rr_key] > 0:
            if male[lo_key] > 0 and male[hi_key] > 0 and female[lo_key] > 0 and female[hi_key] > 0:
                diff = np.log(male[rr_key]) - np.log(female[rr_key])
                se_m = (np.log(male[hi_key]) - np.log(male[lo_key])) / (2 * 1.96)
                se_f = (np.log(female[hi_key]) - np.log(female[lo_key])) / (2 * 1.96)
                var = se_m**2 + se_f**2  # Conservative (assumes independence)
                
                if var > 0:
                    paired_diffs[effect_type].append(diff)
                    paired_vars[effect_type].append(var)

for effect_type in ['heat', 'cold']:
    if len(paired_diffs[effect_type]) < 5:
        continue
    
    meta = meta_random_effects(
        np.array(paired_diffs[effect_type]), 
        np.array(paired_vars[effect_type])
    )
    
    rr_ratio = np.exp(meta['effect'])
    rr_lo = np.exp(meta['effect'] - 1.96 * meta['se'])
    rr_hi = np.exp(meta['effect'] + 1.96 * meta['se'])
    
    sig = "*" if meta['p'] < 0.05 else ""
    print(f"\n  {effect_type.upper()} Paired comparison (k={meta['k']}):")
    print(f"    RR ratio (M/F): {rr_ratio:.3f} [{rr_lo:.3f}-{rr_hi:.3f}], p = {meta['p']:.3f}{sig}")
    
    comparisons[effect_type]['paired'] = {
        'rr_ratio': float(rr_ratio),
        'rr_ratio_lo': float(rr_lo),
        'rr_ratio_hi': float(rr_hi),
        'p': float(meta['p']),
        'k': int(meta['k']),
        'I2': float(meta['I2']),
    }

# =============================================================================
# 11. SAVE RESULTS
# =============================================================================
print("\n[12] Saving results...")

summary = {
    'analysis_date': datetime.now().isoformat(),
    'level': level,
    'sexes': SEXES,
    'dlnm_params': {
        'max_lag': MAX_LAG,
        'var_df': VAR_DF,
        'lag_df': LAG_DF,
        'basis_type': 'natural_spline',
    },
    'approach': 'per_region_dlnm_then_mvmeta_pooling',
    'methodology': 'Gasparrini MVMeta - pool DLNM coefficients not RRs',
    'offset': 'log(pop_elderly * sex_proportion)',
    'percentile_method': 'region_specific',
    'pooled_results': pooled_results,
    'comparisons': comparisons,
}

region_output = {
    sex: {int(k): convert_to_json_serializable(v) for k, v in results.items()}
    for sex, results in sex_results.items()
}

with open(OUTPUT_DIR / f'sex_stratification_v2_results{output_suffix}.json', 'w') as f:
    json.dump(convert_to_json_serializable(summary), f, indent=2)

with open(OUTPUT_DIR / f'sex_stratification_v2_region_effects{output_suffix}.json', 'w') as f:
    json.dump(region_output, f, indent=2)

print(f"\n  Saved to {OUTPUT_DIR}:")
print(f"    - sex_stratification_v2_results{output_suffix}.json")
print(f"    - sex_stratification_v2_region_effects{output_suffix}.json")

print(f"\n{'='*70}")
print("04c: SEX STRATIFICATION COMPLETE")
print(f"{'='*70}")
print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
