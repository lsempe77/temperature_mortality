"""
01f_v2: PROPER CASE-CROSSOVER VALIDATION
=========================================
Validates DLNM temperature effects using TRUE case-crossover design.

v2 CRITICAL FIX (Dec 2025 expert review):
The v1 script was NOT a proper case-crossover - it aggregated deaths to daily
counts and compared high vs low mortality days. That's an ecological analysis,
not a self-matched case-crossover.

TRUE CASE-CROSSOVER DESIGN:
- Unit of analysis: INDIVIDUAL deaths (not daily counts)
- For each death (case) at date t, select control days t' with same:
  - Year-month
  - Day-of-week
- Fit conditional logistic regression: case vs control on temperature
- This removes ALL time-invariant confounding via self-matching

IMPLEMENTATION:
Option A: Build matched dataset, export to R for clogit (gold standard)
Option B: Conditional Poisson approximation (statsmodels, faster)
Option C: Fixed-effects logistic (approximates conditional logit for rare events)

This script implements Option B (Conditional Poisson) with Option C fallback,
plus exports matched data for R validation.

References:
- Maclure (1991) - Case-crossover design
- Levy et al. (2001) - Time-stratified case-crossover
- Armstrong et al. (2014) - Case-crossover for temperature-mortality
- Lu & Zeger (2007) - Conditional Poisson equivalence
"""

import warnings
warnings.filterwarnings('ignore')

import json
import re
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
from typing import Dict, Tuple, Optional
from scipy import stats
import statsmodels.api as sm
from collections import defaultdict
import argparse

# =============================================================================
# CONFIGURATION
# =============================================================================

class Config:
    """Configuration for case-crossover analysis."""
    
    # Control selection
    CONTROL_STRATEGY = 'time_stratified'  # same month, same DOW
    
    # Sample size for individual-level analysis
    # With millions of deaths, we may need to sample for computational feasibility
    MAX_CASES_FULL = 100000      # Use all if fewer than this
    SAMPLE_SIZE = 50000         # Sample this many if more cases
    
    # Temperature exposure
    LAG_DAYS = [0, 1, 2, 3]     # Lag days to consider (0 = same day)
    
    # Bootstrap
    N_BOOTSTRAP = 500
    CI_LEVEL = 0.95
    
    # Minimum controls per case
    MIN_CONTROLS = 2

# =============================================================================
# PATHS
# =============================================================================

SCRIPT_DIR = Path(__file__).parent
PHASE0_RESULTS = SCRIPT_DIR.parent / 'phase0_data_prep' / 'results'
INPUT_DATA = SCRIPT_DIR.parent.parent / 'Input_data'
PHASE1_RESULTS = SCRIPT_DIR / 'results'
OUTPUT_DIR = PHASE1_RESULTS
OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def _resolve_sim_columns(df: pd.DataFrame) -> Optional[Dict[str, str]]:
    """Map SIM raw column names to canonical ones.

    Returns dict with keys: date, age, muni, cause.
    If any cannot be resolved, returns None.
    """
    def normalize(name: str) -> str:
        return re.sub(r"[^0-9A-Za-z]", "", str(name)).lower()

    norm_map = {normalize(c): c for c in df.columns}

    def find(target_keys):
        for key in target_keys:
            if key in norm_map:
                return norm_map[key]
        return None

    # Common SIM patterns
    date_col = find(['dtobito', 'datobito', 'dataobito', 'dto'])
    age_col = find(['idade', 'idadeobito'])
    muni_col = find(['codmunocor', 'codmunres', 'codmunocorr'])
    cause_col = find(['causabas', 'causabasic'])

    if not all([date_col, age_col, muni_col, cause_col]):
        missing = []
        if not date_col:
            missing.append('date')
        if not age_col:
            missing.append('age')
        if not muni_col:
            missing.append('municipality')
        if not cause_col:
            missing.append('cause')
        print(f"      Skipping file: could not find columns for {', '.join(missing)}")
        return None

    return {
        'date': date_col,
        'age': age_col,
        'muni': muni_col,
        'cause': cause_col,
    }


def load_individual_deaths(years: list = None) -> pd.DataFrame:
    """
    Load individual-level death records from SIM data.
    
    Returns DataFrame with columns: date, age, region_code, cause
    """
    print("  Loading individual death records from SIM...")
    
    all_deaths = []
    
    # Find DO files
    do_files = sorted(INPUT_DATA.glob('DO*OPEN.csv'))
    
    if not do_files:
        print("  WARNING: No SIM files found, falling back to aggregated data")
        return None
    
    for f in do_files:
        # Extract year from filename (e.g., DO10OPEN.csv -> 2010)
        year_str = f.stem.replace('DO', '').replace('OPEN', '')
        try:
            year = int('20' + year_str) if len(year_str) == 2 else int(year_str)
        except ValueError:
            continue
        
        if years and year not in years:
            continue
        
        print(f"    Loading {f.name}...", end=' ')

        try:
            # Load full header and infer relevant columns to be robust to naming.
            # Delimiter can be either ';' or ',' depending on year; mimic the
            # robust detection used in Phase 0 aggregation (00n_aggregate_...).

            with open(f, 'r', encoding='latin-1') as fh:
                first_line = fh.readline()

            if '","' in first_line or first_line.count(',') > first_line.count(';'):
                sep = ','
            else:
                sep = ';'

            df = pd.read_csv(
                f,
                sep=sep,
                encoding='latin-1',
                dtype=str,
                low_memory=False
            )

            col_map = _resolve_sim_columns(df)
            if col_map is None:
                continue

            date_raw = col_map['date']
            age_raw = col_map['age']
            muni_raw = col_map['muni']
            cause_raw = col_map['cause']

            # Parse date - handle both DDMMYYYY and ISO (YYYY-MM-DD) formats
            # Mirrors Phase 0 logic from 00n_aggregate_mortality_to_regions.py
            def parse_sim_date(dtobito):
                try:
                    s = str(dtobito).strip()
                    if '-' in s:  # ISO format
                        return pd.Timestamp(s)
                    # DDMMYYYY format
                    s = str(int(float(dtobito))).zfill(8)
                    day, month, year = int(s[0:2]), int(s[2:4]), int(s[4:8])
                    return pd.Timestamp(year=year, month=month, day=day)
                except:
                    return pd.NaT

            df['date'] = df[date_raw].apply(parse_sim_date)

            # Parse age - filter elderly (60+)
            # IDADE format: 4XX = years (e.g., 465 = 65 years), 5XX = 100+ years
            # Matches Phase 0 logic: if IDADE >= 400, age = IDADE - 400
            def parse_age(x):
                try:
                    code = int(float(str(x).split(',')[0].strip()))
                    if code >= 400:
                        return code - 400
                    return np.nan  # Under 1 year or invalid
                except:
                    return np.nan

            df['age'] = df[age_raw].apply(parse_age)
            df = df[df['age'] >= 60]

            # Get municipality code (first 6 digits)
            df['muni_code'] = df[muni_raw].str[:6]

            # Keep needed columns
            df = df[['date', 'age', 'muni_code', cause_raw]].dropna(subset=['date'])
            df.rename(columns={cause_raw: 'cause'}, inplace=True)

            all_deaths.append(df)
            print(f"{len(df):,} elderly deaths")

        except Exception as e:
            print(f"Error: {e}")
            continue
    
    if not all_deaths:
        return None
    
    df_all = pd.concat(all_deaths, ignore_index=True)
    print(f"  Total elderly deaths loaded: {len(df_all):,}")
    
    return df_all


def load_municipality_to_region_mapping() -> Dict[str, str]:
    """Load mapping from municipality to intermediate region."""
    mapping_file = PHASE0_RESULTS / 'municipality_to_intermediate_region.csv'
    
    if mapping_file.exists():
        df = pd.read_csv(mapping_file, dtype=str)
        return dict(zip(df['muni_code'], df['region_code']))
    
    # Try alternative file
    mapping_file = SCRIPT_DIR.parent / 'results' / 'municipality_to_intermediate.csv'
    if mapping_file.exists():
        df = pd.read_csv(mapping_file, dtype=str)
        return dict(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1].astype(str)))
    
    print("  WARNING: Municipality mapping not found")
    return {}


def build_matched_dataset(df_deaths: pd.DataFrame, 
                          df_temp: pd.DataFrame,
                          sample_size: int = None) -> pd.DataFrame:
    """
    Build case-crossover matched dataset.
    
    For each death (case), find control days with same year-month-DOW.
    
    Parameters:
    -----------
    df_deaths : Individual death records with 'date' column
    df_temp : Daily temperature with 'date', 'temp_mean' columns
    sample_size : If specified, sample this many cases
    
    Returns:
    --------
    Matched dataset with columns: case_id, date, is_case, temp_mean
    """
    print("\n[2] Building matched case-control dataset...")
    
    # Sample if needed
    if sample_size and len(df_deaths) > sample_size:
        print(f"  Sampling {sample_size:,} cases from {len(df_deaths):,}")
        df_deaths = df_deaths.sample(n=sample_size, random_state=42)
    
    # Pre-compute temperature lookup
    temp_lookup = df_temp.set_index('date')['temp_mean'].to_dict()
    
    # Pre-compute control days by stratum (year-month-dow)
    # This is much faster than computing per case
    print("  Pre-computing control days by stratum...")
    stratum_controls = defaultdict(list)
    
    for date in df_temp['date'].unique():
        dt = pd.Timestamp(date)
        stratum = (dt.year, dt.month, dt.dayofweek)
        stratum_controls[stratum].append(dt)
    
    # Build matched pairs
    print("  Building matched pairs...")
    matched_records = []
    cases_with_no_controls = 0
    
    # Process in batches for memory efficiency
    batch_size = 10000
    n_cases = len(df_deaths)
    
    for batch_start in range(0, n_cases, batch_size):
        batch_end = min(batch_start + batch_size, n_cases)
        batch = df_deaths.iloc[batch_start:batch_end]
        
        for idx, (_, row) in enumerate(batch.iterrows()):
            case_id = batch_start + idx
            case_date = row['date']
            
            if pd.isna(case_date):
                continue
            
            case_date = pd.Timestamp(case_date)
            stratum = (case_date.year, case_date.month, case_date.dayofweek)
            
            # Get control days (same stratum, excluding case date)
            control_dates = [d for d in stratum_controls.get(stratum, []) 
                            if d != case_date]
            
            if len(control_dates) < Config.MIN_CONTROLS:
                cases_with_no_controls += 1
                continue
            
            # Get case temperature
            case_temp = temp_lookup.get(case_date, np.nan)
            if pd.isna(case_temp):
                continue
            
            # Add case record
            matched_records.append({
                'case_id': case_id,
                'date': case_date,
                'is_case': 1,
                'temp_mean': case_temp,
                'stratum': f"{stratum[0]}_{stratum[1]:02d}_{stratum[2]}"
            })
            
            # Add control records
            for control_date in control_dates:
                control_temp = temp_lookup.get(control_date, np.nan)
                if pd.isna(control_temp):
                    continue
                    
                matched_records.append({
                    'case_id': case_id,
                    'date': control_date,
                    'is_case': 0,
                    'temp_mean': control_temp,
                    'stratum': f"{stratum[0]}_{stratum[1]:02d}_{stratum[2]}"
                })
        
        if (batch_end) % 20000 == 0:
            print(f"    Processed {batch_end:,}/{n_cases:,} cases...")
    
    df_matched = pd.DataFrame(matched_records)
    
    n_cases_included = df_matched['case_id'].nunique()
    n_controls_total = len(df_matched) - n_cases_included
    
    print(f"  Cases included: {n_cases_included:,}")
    print(f"  Cases dropped (insufficient controls): {cases_with_no_controls:,}")
    print(f"  Total control records: {n_controls_total:,}")
    print(f"  Mean controls per case: {n_controls_total / n_cases_included:.1f}")
    
    return df_matched


def fit_conditional_poisson(df_matched: pd.DataFrame) -> Dict:
    """
    Fit conditional Poisson model (approximates conditional logistic).
    
    Model: is_case ~ temp_mean + temp_mean^2 | case_id (fixed effects)
    
    Uses the Chamberlain (1980) result that conditional Poisson = conditional logistic
    for binary outcomes.
    """
    print("\n[3] Fitting conditional Poisson model...")
    
    # Create temperature variables
    df = df_matched.copy()
    
    # Center temperature for numerical stability
    temp_mean = df['temp_mean'].mean()
    temp_std = df['temp_mean'].std()
    df['temp_c'] = (df['temp_mean'] - temp_mean) / temp_std
    df['temp_c2'] = df['temp_c'] ** 2
    
    # Compute temperature percentiles for effect estimation
    temp_percentiles = {
        'p01': df_matched['temp_mean'].quantile(0.01),
        'p025': df_matched['temp_mean'].quantile(0.025),
        'p25': df_matched['temp_mean'].quantile(0.25),
        'p50': df_matched['temp_mean'].quantile(0.50),
        'p75': df_matched['temp_mean'].quantile(0.75),
        'p975': df_matched['temp_mean'].quantile(0.975),
        'p99': df_matched['temp_mean'].quantile(0.99)
    }
    
    # Approach 1: Conditional Poisson via fixed effects
    # We use case_id as fixed effects (strata)
    
    try:
        # Try statsmodels conditional logit if available
        from statsmodels.discrete.conditional_models import ConditionalLogit
        
        # Prepare data
        y = df['is_case'].values
        X = df[['temp_c', 'temp_c2']].values
        groups = df['case_id'].values
        
        model = ConditionalLogit(y, X, groups=groups)
        result = model.fit(disp=False)
        
        print("  Using statsmodels ConditionalLogit")
        
        coef = result.params
        se = result.bse
        
        # Extract coefficients
        beta_temp = coef[0]
        beta_temp2 = coef[1] if len(coef) > 1 else 0
        se_temp = se[0]
        se_temp2 = se[1] if len(se) > 1 else 0
        
        converged = True
        method = 'ConditionalLogit'
        
    except (ImportError, Exception) as e:
        print(f"  ConditionalLogit not available ({e}), using fixed-effects approximation")
        
        # Fallback: Fixed-effects logistic approximation
        # Include case_id dummies (feasible for modest sample sizes)
        
        # For large datasets, use within-transformation
        # Subtract case-specific means (within-case centering)
        
        case_means = df.groupby('case_id')[['temp_c', 'temp_c2']].transform('mean')
        df['temp_c_within'] = df['temp_c'] - case_means['temp_c']
        df['temp_c2_within'] = df['temp_c2'] - case_means['temp_c2']
        
        # Fit logistic on within-transformed data
        X = sm.add_constant(df[['temp_c_within', 'temp_c2_within']])
        y = df['is_case']
        
        try:
            model = sm.Logit(y, X)
            result = model.fit(disp=False, maxiter=100)
            
            beta_temp = result.params.get('temp_c_within', 0)
            beta_temp2 = result.params.get('temp_c2_within', 0)
            se_temp = result.bse.get('temp_c_within', np.nan)
            se_temp2 = result.bse.get('temp_c2_within', np.nan)
            
            converged = result.converged
            method = 'FixedEffectsLogit'
            
        except Exception as e2:
            print(f"  Fixed-effects also failed ({e2}), using simple within-stratum regression")
            
            # Ultimate fallback: within-stratum regression
            from scipy.stats import linregress
            
            slope, intercept, r, p, se = linregress(df['temp_c_within'], df['is_case'])
            beta_temp = slope
            beta_temp2 = 0
            se_temp = se
            se_temp2 = np.nan
            converged = True
            method = 'WithinStratumOLS'
    
    # Calculate odds ratios at key percentiles (relative to P50)
    def calc_or(temp, temp_ref, beta1, beta2, temp_mean, temp_std):
        """Calculate OR for temp vs temp_ref."""
        t = (temp - temp_mean) / temp_std
        t_ref = (temp_ref - temp_mean) / temp_std
        log_or = beta1 * (t - t_ref) + beta2 * (t**2 - t_ref**2)
        return np.exp(log_or)
    
    ref_temp = temp_percentiles['p50']
    
    or_results = {}
    for pname, ptemp in temp_percentiles.items():
        if pname == 'p50':
            continue
        or_val = calc_or(ptemp, ref_temp, beta_temp, beta_temp2, temp_mean, temp_std)
        or_results[pname] = {
            'temperature': float(ptemp),
            'or': float(or_val)
        }
    
    print(f"  Model: {method}")
    print(f"  Converged: {converged}")
    print(f"  β(temp): {beta_temp:.4f} (SE: {se_temp:.4f})")
    if beta_temp2 != 0:
        print(f"  β(temp²): {beta_temp2:.4f} (SE: {se_temp2:.4f})")
    
    return {
        'method': method,
        'converged': converged,
        'coefficients': {
            'temp_linear': float(beta_temp),
            'temp_quadratic': float(beta_temp2),
            'se_linear': float(se_temp),
            'se_quadratic': float(se_temp2) if not np.isnan(se_temp2) else None
        },
        'centering': {
            'temp_mean': float(temp_mean),
            'temp_std': float(temp_std)
        },
        'temperature_percentiles': {k: float(v) for k, v in temp_percentiles.items()},
        'odds_ratios_vs_p50': or_results
    }


def compute_stratified_or(df_matched: pd.DataFrame) -> Dict:
    """
    Compute odds ratios using simple stratified 2x2 tables.
    
    This is an alternative estimation approach that's more transparent.
    """
    print("\n[4] Computing stratified odds ratios (Mantel-Haenszel)...")
    
    # Temperature percentiles
    temp_p50 = df_matched['temp_mean'].quantile(0.50)
    temp_p25 = df_matched['temp_mean'].quantile(0.25)
    temp_p75 = df_matched['temp_mean'].quantile(0.75)
    temp_p975 = df_matched['temp_mean'].quantile(0.975)
    temp_p025 = df_matched['temp_mean'].quantile(0.025)
    
    results = {}
    
    # Heat effect: P97.5 vs P50
    def compute_mh_or(high_thresh, low_thresh, comparison_name):
        """Compute Mantel-Haenszel OR for high vs low temp."""
        
        # For each case, compare case temp to control temps
        case_records = df_matched[df_matched['is_case'] == 1]
        control_records = df_matched[df_matched['is_case'] == 0]
        
        # Aggregate by case_id
        numerator = 0
        denominator = 0
        
        for case_id in case_records['case_id'].unique():
            case_temp = case_records[case_records['case_id'] == case_id]['temp_mean'].values[0]
            control_temps = control_records[control_records['case_id'] == case_id]['temp_mean'].values
            
            if len(control_temps) == 0:
                continue
            
            # Is case exposed (above high threshold)?
            case_exposed = 1 if case_temp > high_thresh else 0
            
            # Proportion of controls exposed
            n_controls = len(control_temps)
            controls_exposed = np.sum(control_temps > high_thresh)
            controls_unexposed = n_controls - controls_exposed
            
            if case_exposed:
                numerator += controls_unexposed / n_controls
            else:
                denominator += controls_exposed / n_controls
        
        if denominator > 0:
            mh_or = numerator / denominator
        else:
            mh_or = np.nan
        
        return mh_or
    
    # Heat: >P97.5 vs <P97.5
    heat_or = compute_mh_or(temp_p975, temp_p50, 'heat')
    print(f"  Heat OR (>P97.5): {heat_or:.3f}")
    
    # Cold: <P2.5 vs >P2.5 (invert interpretation)
    # For cold, we look at case being below threshold
    cold_or_inv = compute_mh_or(temp_p025, temp_p50, 'cold')
    cold_or = 1 / cold_or_inv if cold_or_inv > 0 else np.nan
    print(f"  Cold OR (<P2.5): {cold_or:.3f}")
    
    return {
        'heat': {
            'threshold': float(temp_p975),
            'or': float(heat_or) if not np.isnan(heat_or) else None,
            'comparison': '>P97.5 vs ≤P97.5'
        },
        'cold': {
            'threshold': float(temp_p025),
            'or': float(cold_or) if not np.isnan(cold_or) else None,
            'comparison': '<P2.5 vs ≥P2.5'
        }
    }


def run_diagnostics(df_matched: pd.DataFrame) -> Dict:
    """Run case-crossover diagnostics."""
    print("\n[5] Running diagnostics...")
    
    # Controls per case distribution
    controls_per_case = df_matched.groupby('case_id').apply(
        lambda x: (x['is_case'] == 0).sum()
    )
    
    # Temperature distribution: cases vs controls
    case_temps = df_matched[df_matched['is_case'] == 1]['temp_mean']
    control_temps = df_matched[df_matched['is_case'] == 0]['temp_mean']
    
    # KS test: are case temperatures different from control temperatures?
    ks_stat, ks_pval = stats.ks_2samp(case_temps, control_temps)
    
    # Strata with zero controls
    zero_control_cases = (controls_per_case == 0).sum()
    
    diagnostics = {
        'n_cases': int(df_matched['case_id'].nunique()),
        'n_total_records': len(df_matched),
        'controls_per_case': {
            'mean': float(controls_per_case.mean()),
            'median': float(controls_per_case.median()),
            'min': int(controls_per_case.min()),
            'max': int(controls_per_case.max()),
            'std': float(controls_per_case.std())
        },
        'zero_control_cases': int(zero_control_cases),
        'temperature_distribution': {
            'case_mean': float(case_temps.mean()),
            'control_mean': float(control_temps.mean()),
            'case_std': float(case_temps.std()),
            'control_std': float(control_temps.std()),
            'ks_statistic': float(ks_stat),
            'ks_pvalue': float(ks_pval)
        }
    }
    
    print(f"  Cases analyzed: {diagnostics['n_cases']:,}")
    print(f"  Mean controls per case: {diagnostics['controls_per_case']['mean']:.1f}")
    print(f"  Case mean temp: {case_temps.mean():.1f}°C")
    print(f"  Control mean temp: {control_temps.mean():.1f}°C")
    print(f"  KS test p-value: {ks_pval:.4f} {'✓ significant' if ks_pval < 0.05 else ''}")
    
    return diagnostics


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def main(use_individual: bool = True, sample_size: int = None):
    """
    Run proper case-crossover analysis.
    
    Parameters:
    -----------
    use_individual : Use individual death records (True) or fall back to aggregated
    sample_size : Sample size for individual analysis (None = use config)
    """
    
    print("=" * 70)
    print("01f_v2: PROPER CASE-CROSSOVER VALIDATION")
    print("=" * 70)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Using individual-level data: {use_individual}")
    print()
    
    # =========================================================================
    # LOAD DATA
    # =========================================================================
    print("[1] Loading data...")
    
    # Load temperature data (national average)
    df_temp = pd.read_parquet(PHASE0_RESULTS / 'era5_intermediate_daily.parquet')
    df_temp['date'] = pd.to_datetime(df_temp['date'])
    
    # Aggregate to national daily temperature
    df_temp_national = df_temp.groupby('date').agg({
        'temp_mean': 'mean'
    }).reset_index()
    
    print(f"  Temperature data: {len(df_temp_national):,} days")
    print(f"  Date range: {df_temp_national['date'].min().date()} to {df_temp_national['date'].max().date()}")
    
    # Load death records
    if use_individual:
        df_deaths = load_individual_deaths()
        
        if df_deaths is None or len(df_deaths) == 0:
            print("  Falling back to aggregated approach")
            use_individual = False
    
    if not use_individual:
        # Fall back to aggregated (but with proper warning)
        print("\n  ⚠️ WARNING: Using aggregated data - this is NOT a true case-crossover!")
        print("     Results will be labeled as 'ecological within-stratum correlation'")
        
        df_mort = pd.read_parquet(PHASE0_RESULTS / 'mortality_regional_daily_elderly.parquet')
        df_mort['date'] = pd.to_datetime(df_mort['date'])
        
        # Create pseudo-individual records (one per death on each day)
        # This is still ecological but better than pure aggregation
        df_national = df_mort.groupby('date')['deaths_elderly'].sum().reset_index()
        
        pseudo_deaths = []
        for _, row in df_national.iterrows():
            for _ in range(int(row['deaths_elderly'])):
                pseudo_deaths.append({'date': row['date']})
        
        df_deaths = pd.DataFrame(pseudo_deaths)
        print(f"  Created {len(df_deaths):,} pseudo-individual records from aggregated data")
    
    print(f"  Total deaths for analysis: {len(df_deaths):,}")
    
    # =========================================================================
    # BUILD MATCHED DATASET
    # =========================================================================
    
    if sample_size is None:
        sample_size = Config.SAMPLE_SIZE if len(df_deaths) > Config.MAX_CASES_FULL else None
    
    df_matched = build_matched_dataset(df_deaths, df_temp_national, sample_size)
    
    # =========================================================================
    # FIT CONDITIONAL MODEL
    # =========================================================================
    
    model_results = fit_conditional_poisson(df_matched)
    
    # =========================================================================
    # COMPUTE STRATIFIED OR
    # =========================================================================
    
    stratified_or = compute_stratified_or(df_matched)
    
    # =========================================================================
    # DIAGNOSTICS
    # =========================================================================
    
    diagnostics = run_diagnostics(df_matched)
    
    # =========================================================================
    # COMPARE WITH DLNM
    # =========================================================================
    print("\n[6] Comparing with DLNM results...")
    
    comparison = {}
    try:
        # Prefer compact pooled files where available, with v2 first
        pooled_v2_file = PHASE1_RESULTS / 'dlnm_v2_intermediate_pooled.json'
        results_v2_file = PHASE1_RESULTS / 'dlnm_v2_intermediate_results.json'
        legacy_pooled_file = PHASE1_RESULTS / 'dlnm_intermediate_pooled.json'

        if pooled_v2_file.exists():
            dlnm_file = pooled_v2_file
        elif results_v2_file.exists():
            dlnm_file = results_v2_file
        else:
            dlnm_file = legacy_pooled_file

        with open(dlnm_file, 'r') as f:
            dlnm_results = json.load(f)

        dlnm_heat = None
        dlnm_cold = None

        # v2 structure: { ..., 'pooled_results': { 'p99': {'pooled_rr': ...}, ... } }
        if isinstance(dlnm_results, dict) and 'pooled_results' in dlnm_results:
            pooled_block = dlnm_results.get('pooled_results', {})
            dlnm_heat = pooled_block.get('p99', {}).get('pooled_rr') or \
                        pooled_block.get('p99', {}).get('rr')
            dlnm_cold = pooled_block.get('p1', {}).get('pooled_rr') or \
                        pooled_block.get('p1', {}).get('rr')
        else:
            # Legacy pooled structure: top-level percentiles
            dlnm_heat = dlnm_results.get('p99', {}).get('pooled_rr') or \
                        dlnm_results.get('p99', {}).get('rr')
            dlnm_cold = dlnm_results.get('p1', {}).get('pooled_rr') or \
                        dlnm_results.get('p1', {}).get('rr')

            # Older flat naming convention
            if dlnm_heat is None and 'pooled_rr_p99' in dlnm_results:
                dlnm_heat = dlnm_results['pooled_rr_p99']
            if dlnm_cold is None and 'pooled_rr_p1' in dlnm_results:
                dlnm_cold = dlnm_results['pooled_rr_p1']
        
        # Get case-crossover ORs
        cc_heat = model_results['odds_ratios_vs_p50'].get('p99', {}).get('or', None)
        cc_cold = model_results['odds_ratios_vs_p50'].get('p01', {}).get('or', None)
        
        comparison = {
            'dlnm_heat_rr_p99': dlnm_heat,
            'dlnm_cold_rr_p1': dlnm_cold,
            'case_crossover_heat_or_p99': cc_heat,
            'case_crossover_cold_or_p01': cc_cold,
            'heat_ratio': cc_heat / dlnm_heat if cc_heat and dlnm_heat else None,
            'cold_ratio': cc_cold / dlnm_cold if cc_cold and dlnm_cold else None
        }
        
        print(f"\n  DLNM Heat RR (P99): {dlnm_heat:.3f}" if dlnm_heat else "\n  DLNM Heat: N/A")
        print(f"  Case-Crossover Heat OR (P99): {cc_heat:.3f}" if cc_heat else "  CC Heat: N/A")
        print(f"\n  DLNM Cold RR (P1): {dlnm_cold:.3f}" if dlnm_cold else "\n  DLNM Cold: N/A")
        print(f"  Case-Crossover Cold OR (P1): {cc_cold:.3f}" if cc_cold else "  CC Cold: N/A")
        
    except FileNotFoundError:
        print("  DLNM results file not found")
    
    # =========================================================================
    # SAVE RESULTS
    # =========================================================================
    print("\n[7] Saving results...")
    
    def convert_for_json(obj):
        if isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, pd.Timestamp):
            return obj.isoformat()
        elif isinstance(obj, dict):
            return {k: convert_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_for_json(i) for i in obj]
        return obj
    
    results = {
        'method': 'Time-stratified case-crossover (individual-level)' if use_individual 
                  else 'Ecological within-stratum analysis (NOT true case-crossover)',
        'version': 'v2',
        'is_true_case_crossover': use_individual,
        'design': {
            'control_strategy': Config.CONTROL_STRATEGY,
            'description': 'Same year-month, same day-of-week',
            'sample_size_used': sample_size or len(df_deaths)
        },
        'model': model_results,
        'stratified_or': stratified_or,
        'diagnostics': diagnostics,
        'comparison_with_dlnm': comparison,
        'timestamp': datetime.now().isoformat()
    }
    
    with open(OUTPUT_DIR / 'case_crossover_v2.json', 'w') as f:
        json.dump(convert_for_json(results), f, indent=2)
    
    # Export matched dataset for R validation
    export_file = OUTPUT_DIR / 'case_crossover_matched_data.csv'
    df_matched_export = df_matched[['case_id', 'date', 'is_case', 'temp_mean', 'stratum']].copy()
    df_matched_export['date'] = df_matched_export['date'].dt.strftime('%Y-%m-%d')
    df_matched_export.to_csv(export_file, index=False)
    
    print(f"  Saved: case_crossover_v2.json")
    print(f"  Saved: case_crossover_matched_data.csv (for R validation)")
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    
    print("\n" + "=" * 70)
    print("CASE-CROSSOVER VALIDATION SUMMARY (v2)")
    print("=" * 70)
    
    if not use_individual:
        print("\n⚠️ WARNING: This analysis used aggregated data pseudo-expanded to")
        print("   individual records. This is NOT a true case-crossover design.")
        print("   Interpret as ecological within-stratum association only.\n")
    
    cc_heat = model_results['odds_ratios_vs_p50'].get('p99', {}).get('or', None)
    cc_cold = model_results['odds_ratios_vs_p50'].get('p01', {}).get('or', None)
    
    # Format values with None handling
    cc_heat_str = f"{cc_heat:.3f}" if cc_heat is not None else "N/A"
    cc_cold_str = f"{cc_cold:.3f}" if cc_cold is not None else "N/A"
    heat_or = stratified_or['heat']['or']
    cold_or = stratified_or['cold']['or']
    heat_or_str = f"{heat_or:.3f}" if heat_or is not None else "N/A"
    cold_or_str = f"{cold_or:.3f}" if cold_or is not None else "N/A"
    
    print(f"""
Case-Crossover Design:
  - Cases: {diagnostics['n_cases']:,} deaths
  - Controls: Same year-month, same day-of-week
  - Mean controls per case: {diagnostics['controls_per_case']['mean']:.1f}
  - Model: {model_results['method']}

Results (Odds Ratios vs P50 reference):
  - Heat (P99 vs P50): OR = {cc_heat_str}
  - Cold (P1 vs P50): OR = {cc_cold_str}

Stratified OR (Mantel-Haenszel approximation):
  - Heat (>P97.5): OR = {heat_or_str}
  - Cold (<P2.5): OR = {cold_or_str}

Temperature Distribution:
  - Case mean: {diagnostics['temperature_distribution']['case_mean']:.1f}°C
  - Control mean: {diagnostics['temperature_distribution']['control_mean']:.1f}°C
  - KS test p-value: {diagnostics['temperature_distribution']['ks_pvalue']:.4f}

Validation:
  - Cold > Heat pattern: {'✓ CONFIRMED' if (cc_cold and cc_heat and cc_cold > cc_heat) else '?'}
  - Consistent with DLNM: {comparison.get('heat_ratio', 'N/A')}
""")
    
    print("=" * 70)
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)
    
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Case-Crossover Validation v2')
    parser.add_argument('--aggregated', action='store_true',
                        help='Use aggregated data (NOT recommended, for testing only)')
    parser.add_argument('--sample', type=int, default=None,
                        help='Sample size for analysis (default: auto)')
    
    args = parser.parse_args()
    
    main(use_individual=not args.aggregated, sample_size=args.sample)
