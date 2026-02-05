"""
04b: AGE STRATIFICATION (REGIONAL DLNM) - v2
=============================================
Tests effect modification by elderly age groups: 60-69, 70-79, 80+

v2 IMPROVEMENTS (Dec 2025):
---------------------------
1. Per-region DLNM fits then meta-analysis (not pooled model with dummies)
2. Uses natural spline cross-basis from utils
3. Population offset in all models
4. Full covariance delta method for CIs
5. Region-specific percentiles
6. Proper heterogeneity testing (Q-test between age groups)

Two-stage approach (per age group):
1. Fit DLNM in each region → extract heat/cold RRs
2. Meta-analyze across regions (DerSimonian-Laird)
3. Compare pooled effects between age groups
"""

import warnings
warnings.filterwarnings('ignore')

import sys
import json
import pandas as pd
import numpy as np
import statsmodels.api as sm
from datetime import datetime
from pathlib import Path
from scipy import stats
from typing import Dict, Any, Optional, List, Tuple
import argparse

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.dlnm_module import (
    fit_region_dlnm,
    predict_cumulative_rr,
    find_mmt,
    meta_random_effects,
    convert_to_json_serializable,
    # MVMeta for proper coefficient pooling (Gasparrini methodology)
    mvmeta_pool_regions,
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
MIN_OBS = 365  # Minimum obs per region-age group

AGE_GROUPS = ['60-69', '70-79', '80+']

parser = argparse.ArgumentParser(description='Age-stratified regional DLNM meta-analysis (v2)')
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
print("04b: AGE STRATIFICATION (REGIONAL DLNM) - v2")
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
        # Use explicit municipality → immediate mapping from phase0 results
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
# 2. PROCESS MORTALITY DATA BY AGE GROUP
# =============================================================================
print("\n[2] Processing mortality data by age group...")

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

def assign_age_group(age):
    """Assign elderly age groups."""
    if age is None or age < 60:
        return None
    elif age < 70:
        return '60-69'
    elif age < 80:
        return '70-79'
    else:
        return '80+'

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
age_mortality = []
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
                        usecols=['DTOBITO', 'IDADE', 'CODMUNRES'])
        
        df['date'] = df['DTOBITO'].apply(parse_date)
        df['age'] = df['IDADE'].apply(parse_age)
        df['age_group'] = df['age'].apply(assign_age_group)
        
        df = df[df['age_group'].notna() & df['date'].notna()].copy()
        df['muni_code'] = df['CODMUNRES'].astype(str).str[:6]
        
        if muni_to_region is not None:
            df['region_code'] = df['muni_code'].map(muni_to_region)
        else:
            df['region_code'] = df['muni_code'].str[:2].astype(int) * 100
        
        df = df.dropna(subset=['region_code'])
        df['region_code'] = df['region_code'].astype(int)
        
        age_mortality.append(df[['date', 'region_code', 'age_group']])
        print(f" {len(df):,} elderly deaths")
        
    except Exception as e:
        print(f" Error: {e}")

if age_mortality:
    mort_all = pd.concat(age_mortality, ignore_index=True)
    print(f"\n  Total elderly deaths: {len(mort_all):,}")
    for ag in AGE_GROUPS:
        n = (mort_all['age_group'] == ag).sum()
        print(f"    {ag}: {n:,}")
else:
    raise ValueError("No mortality data loaded!")

# Aggregate to region-day-age_group
print("\n[3] Aggregating to region-day-age_group...")

mort_agg = mort_all.groupby(['date', 'region_code', 'age_group']).size().reset_index(name='deaths')
print(f"  Aggregated records: {len(mort_agg):,}")

# =============================================================================
# 3. LOAD TEMPERATURE DATA
# =============================================================================
print("\n[4] Loading temperature data...")

temp_df = pd.read_parquet(PHASE0_RESULTS / era5_filename)
temp_df['date'] = pd.to_datetime(temp_df['date'])
print(f"  Temperature records: {len(temp_df):,}")

# =============================================================================
# 4. LOAD POPULATION BY AGE GROUP (approximate)
# =============================================================================
print("\n[5] Loading population data...")

# Try to load elderly population at chosen spatial level, fallback to proportional allocation
try:
    covariates = pd.read_csv(PHASE0_RESULTS / pop_covariates_filename)
    pop_elderly_total = covariates.set_index(pop_region_col)['pop_elderly'].to_dict()

    # Approximate age group proportions (based on IBGE data)
    # 60-69: ~50%, 70-79: ~30%, 80+: ~20%
    age_proportions = {'60-69': 0.50, '70-79': 0.30, '80+': 0.20}
    print("  Using approximate age group proportions for population offset")
except Exception:
    pop_elderly_total = {}
    age_proportions = {'60-69': 0.50, '70-79': 0.30, '80+': 0.20}
    print("  WARNING: No population data found; using default offsets")

# =============================================================================
# 5. PREPARE DATA FOR EACH AGE GROUP
# =============================================================================
print("\n[6] Preparing analysis datasets...")

datasets = {}

for ag in AGE_GROUPS:
    mort_ag = mort_agg[mort_agg['age_group'] == ag].copy()
    mort_ag = mort_ag.rename(columns={'deaths': 'deaths_elderly'})
    
    # Merge with temperature
    df_ag = mort_ag.merge(temp_df, on=['date', 'region_code'], how='inner')
    
    # Add time variables
    df_ag['dow'] = df_ag['date'].dt.dayofweek
    df_ag['month'] = df_ag['date'].dt.month
    df_ag['year'] = df_ag['date'].dt.year
    
    # Add population (proportional)
    if pop_elderly_total:
        df_ag['pop_elderly'] = df_ag['region_code'].map(pop_elderly_total) * age_proportions[ag]
        df_ag['pop_elderly'] = df_ag['pop_elderly'].fillna(10000)  # Default
    else:
        df_ag['pop_elderly'] = 10000  # Placeholder
    
    datasets[ag] = df_ag
    print(f"  {ag}: {len(df_ag):,} obs, mean deaths = {df_ag['deaths_elderly'].mean():.2f}")

# =============================================================================
# 6. COMPUTE REGION-SPECIFIC PERCENTILES
# =============================================================================
print("\n[7] Computing region-specific temperature percentiles...")

# Use combined temperature data
all_regions = set()
for ag in AGE_GROUPS:
    all_regions.update(datasets[ag]['region_code'].unique())

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
# 7. FIT PER-REGION DLNM FOR EACH AGE GROUP
# =============================================================================
print("\n[8] Fitting per-region DLNMs for each age group...")

def fit_region_age_group(df_region: pd.DataFrame, pctiles: Dict) -> Optional[Dict]:
    """
    Fit DLNM for a region and extract effects relative to MMT.
    
    Following Gasparrini methodology:
    - Find MMT from fitted curve
    - Compare P99 (heat) and P1 (cold) to MMT
    - Store cross-basis coefficients for MVMeta pooling
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

age_group_results = {}

for ag in AGE_GROUPS:
    print(f"\n  === Age Group {ag} ===")
    df_ag = datasets[ag]
    regions = df_ag['region_code'].unique()
    
    region_results = {}
    successful = 0
    
    for i, region in enumerate(regions):
        if (i + 1) % 30 == 0:
            print(f"    Region {i+1}/{len(regions)}...")
        
        if region not in region_percentiles:
            continue
        
        df_region = df_ag[df_ag['region_code'] == region].copy()
        result = fit_region_age_group(df_region, region_percentiles[region])
        
        if result is not None:
            region_results[int(region)] = result
            successful += 1
    
    age_group_results[ag] = region_results
    print(f"    Fitted: {successful}/{len(regions)} regions")

# =============================================================================
# 8. META-ANALYZE WITHIN EACH AGE GROUP USING MVMETA
# =============================================================================
print("\n[9] Meta-analyzing results within each age group using MVMeta...")

pooled_results = {}

for ag in AGE_GROUPS:
    print(f"\n  === {ag} ===")
    region_results = age_group_results[ag]
    
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
            
            pooled_results[ag] = {
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
            
            pooled_results[ag] = {
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
    
    print(f"    Heat RR: {pooled_results[ag]['heat']['rr']:.3f} "
          f"[{pooled_results[ag]['heat']['rr_lo']:.3f}-{pooled_results[ag]['heat']['rr_hi']:.3f}] "
          f"(k={pooled_results[ag]['heat']['k']})")
    print(f"    Cold RR: {pooled_results[ag]['cold']['rr']:.3f} "
          f"[{pooled_results[ag]['cold']['rr_lo']:.3f}-{pooled_results[ag]['cold']['rr_hi']:.3f}] "
          f"(k={pooled_results[ag]['cold']['k']})")

# =============================================================================
# 9. TEST FOR HETEROGENEITY ACROSS AGE GROUPS
# =============================================================================
print("\n[10] Testing heterogeneity across age groups...")

def test_heterogeneity(group_results: Dict, effect_type: str) -> Dict:
    """Cochran's Q test for heterogeneity across groups.
    Works with MVMeta results that have RR and CIs."""
    log_rrs = []
    ses = []
    groups = []
    
    for ag, res in group_results.items():
        if res is None or effect_type not in res:
            continue
        eff = res[effect_type]
        rr = eff.get('rr', 1.0)
        rr_lo = eff.get('rr_lo', rr)
        rr_hi = eff.get('rr_hi', rr)
        
        if rr > 0 and rr_lo > 0 and rr_hi > 0:
            log_rr = np.log(rr)
            se = (np.log(rr_hi) - np.log(rr_lo)) / (2 * 1.96)
            if se > 0.001:  # Valid SE
                log_rrs.append(log_rr)
                ses.append(se)
                groups.append(ag)
    
    if len(log_rrs) < 2:
        return {'Q': np.nan, 'df': 0, 'p': np.nan, 'I2': np.nan, 'groups': groups}
    
    log_rrs = np.array(log_rrs)
    variances = np.array(ses) ** 2
    weights = 1.0 / variances
    
    # Fixed-effect pooled estimate
    pooled = np.sum(weights * log_rrs) / np.sum(weights)
    
    # Q statistic
    Q = np.sum(weights * (log_rrs - pooled) ** 2)
    df = len(log_rrs) - 1
    p = 1 - stats.chi2.cdf(Q, df)
    I2 = max(0, (Q - df) / Q * 100) if Q > 0 else 0
    
    return {'Q': float(Q), 'df': int(df), 'p': float(p), 'I2': float(I2), 'groups': groups}

het_heat = test_heterogeneity(pooled_results, 'heat')
het_cold = test_heterogeneity(pooled_results, 'cold')

print(f"\n  Heat effect heterogeneity across age groups:")
print(f"    Q = {het_heat['Q']:.2f}, df = {het_heat['df']}, p = {het_heat['p']:.3f}")
print(f"    I² = {het_heat['I2']:.1f}%")

print(f"\n  Cold effect heterogeneity across age groups:")
print(f"    Q = {het_cold['Q']:.2f}, df = {het_cold['df']}, p = {het_cold['p']:.3f}")
print(f"    I² = {het_cold['I2']:.1f}%")

# =============================================================================
# 10. PAIRWISE COMPARISONS
# =============================================================================
print("\n[11] Pairwise comparisons between age groups...")

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

for effect_type in ['heat', 'cold']:
    comparisons[effect_type] = {}
    
    # 80+ vs 60-69
    if '80+' in pooled_results and '60-69' in pooled_results:
        r1, r2 = pooled_results['80+'][effect_type], pooled_results['60-69'][effect_type]
        log_rr1, se1 = get_log_rr_se(r1)
        log_rr2, se2 = get_log_rr_se(r2)
        
        diff = log_rr1 - log_rr2
        se_diff = np.sqrt(se1**2 + se2**2)
        z = diff / se_diff if se_diff > 0 else 0
        p = 2 * (1 - stats.norm.cdf(abs(z)))
        
        rr_ratio = np.exp(diff)
        rr_ratio_lo = np.exp(diff - 1.96 * se_diff)
        rr_ratio_hi = np.exp(diff + 1.96 * se_diff)
        
        comparisons[effect_type]['80+_vs_60-69'] = {
            'rr_ratio': float(rr_ratio),
            'rr_ratio_lo': float(rr_ratio_lo),
            'rr_ratio_hi': float(rr_ratio_hi),
            'z': float(z),
            'p': float(p),
            'significant': bool(p < 0.05),
        }
        
        sig = "*" if p < 0.05 else ""
        print(f"\n  {effect_type.upper()}: 80+ vs 60-69")
        print(f"    RR ratio: {rr_ratio:.3f} [{rr_ratio_lo:.3f}-{rr_ratio_hi:.3f}], p = {p:.3f}{sig}")

# =============================================================================
# 11. SAVE RESULTS
# =============================================================================
print("\n[12] Saving results...")

summary = {
    'analysis_date': datetime.now().isoformat(),
    'level': level,
    'age_groups': AGE_GROUPS,
    'dlnm_params': {
        'max_lag': MAX_LAG,
        'var_df': VAR_DF,
        'lag_df': LAG_DF,
        'basis_type': 'natural_spline',
    },
    'approach': 'per_region_dlnm_then_mvmeta_pooling',
    'methodology': 'Gasparrini MVMeta - pool DLNM coefficients not RRs',
    'offset': 'log(pop_elderly * age_proportion)',
    'percentile_method': 'region_specific',
    'pooled_results': pooled_results,
    'heterogeneity': {
        'heat': het_heat,
        'cold': het_cold,
    },
    'pairwise_comparisons': comparisons,
}

# Region-level results by age group
region_output = {
    ag: {int(k): convert_to_json_serializable(v) for k, v in results.items()}
    for ag, results in age_group_results.items()
}

with open(OUTPUT_DIR / f'age_stratification_v2_results{output_suffix}.json', 'w') as f:
    json.dump(convert_to_json_serializable(summary), f, indent=2)

with open(OUTPUT_DIR / f'age_stratification_v2_region_effects{output_suffix}.json', 'w') as f:
    json.dump(region_output, f, indent=2)

print(f"\n  Saved to {OUTPUT_DIR}:")
print(f"    - age_stratification_v2_results{output_suffix}.json")
print(f"    - age_stratification_v2_region_effects{output_suffix}.json")

print(f"\n{'='*70}")
print("04b: AGE STRATIFICATION COMPLETE")
print(f"{'='*70}")
print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
