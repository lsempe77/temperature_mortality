#!/usr/bin/env python3
"""
Luz para Todos Staggered Difference-in-Differences Analysis
============================================================

Estimates the causal effect of rural electrification on heat-related mortality
using staggered DiD design with modern estimators.

Identification Strategy:
- Luz para Todos (LPT) rolled out to different municipalities at different times
- Staggered adoption creates quasi-experimental variation
- Treatment: Access to electricity → Air conditioning → Heat adaptation

Methods:
1. Two-Way Fixed Effects (TWFE) - traditional approach
2. Callaway-Sant'Anna estimator - robust to heterogeneous treatment effects
3. Event study - dynamic treatment effects
4. Triple difference - electricity × heat days × post-treatment

References:
- Callaway & Sant'Anna (2021): Difference-in-Differences with Multiple Time Periods
- de Chaisemartin & D'Haultfoeuille (2020): Two-Way Fixed Effects Estimators
"""

import pandas as pd
import numpy as np
from pathlib import Path
import statsmodels.api as sm
from scipy import stats
import argparse
import json
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"
OUTPUT_DIR = RESULTS_DIR


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def convert_to_json_serializable(obj):
    """Convert numpy types to JSON serializable types."""
    if isinstance(obj, dict):
        return {convert_to_json_serializable(k): convert_to_json_serializable(v) 
                for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_json_serializable(item) for item in obj]
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, pd.Timestamp):
        return obj.isoformat()
    elif pd.isna(obj):
        return None
    else:
        return obj


def load_panel(level='intermediate', monthly=False):
    """Load causal analysis panel."""
    suffix = '_monthly' if monthly else ''
    filename = f"causal_panel{suffix}_{level}.parquet"
    print(f"Loading {filename}...")
    df = pd.read_parquet(RESULTS_DIR / filename)
    print(f"  Shape: {df.shape}")
    return df


def load_lpt_treatment(level='intermediate'):
    """Load LPT treatment data."""
    print(f"Loading {level} LPT treatment data...")
    df = pd.read_parquet(RESULTS_DIR / f"lpt_{level}_treatment.parquet")
    return df


# =============================================================================
# TWO-WAY FIXED EFFECTS (TRADITIONAL DID)
# =============================================================================

def run_twfe_did(df, outcome='deaths_all', treatment='post_lpt'):
    """
    Run traditional Two-Way Fixed Effects DiD.
    
    Model: Y_{rt} = α_r + δ_t + β*Treat_{rt} + ε_{rt}
    
    Note: TWFE can be biased with staggered treatment and heterogeneous effects.
    We include it for comparison but prefer the robust estimators below.
    """
    print(f"\n--- TWFE DiD: {treatment} → {outcome} ---")
    
    # Prepare data
    data = df.dropna(subset=[treatment, outcome]).copy()
    
    if treatment not in data.columns:
        print(f"  Treatment variable {treatment} not found.")
        return None
    
    # Create fixed effects
    data['region_fe'] = pd.Categorical(data['region_code']).codes
    data['time_fe'] = pd.Categorical(data['date']).codes
    
    # Demean by region and time (within estimator)
    # For computational efficiency, we demean by region only
    for col in [outcome, treatment]:
        data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    # Add month and DOW controls
    controls = []
    for m in range(2, 13):
        data[f'month_{m}'] = (data['month'] == m).astype(int)
        controls.append(f'month_{m}')
    
    if 'dow' in data.columns:
        for d in range(1, 7):
            data[f'dow_{d}'] = (data['dow'] == d).astype(int)
            controls.append(f'dow_{d}')
    
    # Demean controls
    for col in controls:
        data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    # Run regression
    y = data[f'{outcome}_dm']
    X_cols = [f'{treatment}_dm'] + [f'{c}_dm' for c in controls]
    X = data[X_cols]
    X = sm.add_constant(X)
    
    model = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    results = {
        'method': 'TWFE',
        'treatment': treatment,
        'outcome': outcome,
        'n_obs': int(model.nobs),
        'n_regions': int(data['region_code'].nunique()),
        'n_treated_regions': int(data[data[treatment] == 1]['region_code'].nunique()),
        'att': float(model.params[f'{treatment}_dm']),
        'std_error': float(model.bse[f'{treatment}_dm']),
        't_stat': float(model.tvalues[f'{treatment}_dm']),
        'p_value': float(model.pvalues[f'{treatment}_dm']),
        'ci_lower': float(model.conf_int().loc[f'{treatment}_dm', 0]),
        'ci_upper': float(model.conf_int().loc[f'{treatment}_dm', 1]),
        'r_squared': float(model.rsquared),
    }
    
    print(f"  ATT: {results['att']:.4f} (SE: {results['std_error']:.4f})")
    print(f"  95% CI: [{results['ci_lower']:.4f}, {results['ci_upper']:.4f}]")
    print(f"  p-value: {results['p_value']:.4f}")
    
    return results


# =============================================================================
# COHORT-SPECIFIC DiD (Callaway-Sant'Anna Style)
# =============================================================================

def calculate_cohort_att(df, cohort_year, outcome='deaths_all', never_treated=None):
    """
    Calculate ATT for a specific cohort relative to never-treated.
    
    Uses 2x2 DiD for each cohort:
    - Treatment group: Units treated in cohort_year
    - Control group: Never-treated units (or not-yet-treated)
    - Pre: Before cohort_year
    - Post: After cohort_year
    """
    # Filter to cohort and control groups
    cohort = df[df['treatment_cohort'] == cohort_year].copy()
    
    if never_treated is not None:
        control = never_treated.copy()
    else:
        # Use not-yet-treated as controls (more conservative)
        control = df[df['first_treatment_year'] > cohort_year].copy()
    
    if len(cohort) == 0 or len(control) == 0:
        return None
    
    # Define pre/post periods
    cohort['post'] = (cohort['year'] >= cohort_year).astype(int)
    control['post'] = (control['year'] >= cohort_year).astype(int)
    
    # Calculate mean outcomes
    cohort_means = cohort.groupby(['region_code', 'post'])[outcome].mean().reset_index()
    control_means = control.groupby(['region_code', 'post'])[outcome].mean().reset_index()
    
    # Calculate DiD
    cohort_pre = cohort_means[cohort_means['post'] == 0][outcome].mean()
    cohort_post = cohort_means[cohort_means['post'] == 1][outcome].mean()
    control_pre = control_means[control_means['post'] == 0][outcome].mean()
    control_post = control_means[control_means['post'] == 1][outcome].mean()
    
    did_estimate = (cohort_post - cohort_pre) - (control_post - control_pre)
    
    # Bootstrap for standard errors
    n_bootstrap = 100
    bootstrap_estimates = []
    
    for _ in range(n_bootstrap):
        # Resample regions
        cohort_regions = np.random.choice(cohort['region_code'].unique(), 
                                          size=len(cohort['region_code'].unique()), replace=True)
        control_regions = np.random.choice(control['region_code'].unique(), 
                                           size=len(control['region_code'].unique()), replace=True)
        
        # Calculate bootstrap DID
        cohort_boot = cohort[cohort['region_code'].isin(cohort_regions)]
        control_boot = control[control['region_code'].isin(control_regions)]
        
        if len(cohort_boot) > 0 and len(control_boot) > 0:
            cb_pre = cohort_boot[cohort_boot['post'] == 0][outcome].mean()
            cb_post = cohort_boot[cohort_boot['post'] == 1][outcome].mean()
            ct_pre = control_boot[control_boot['post'] == 0][outcome].mean()
            ct_post = control_boot[control_boot['post'] == 1][outcome].mean()
            
            boot_did = (cb_post - cb_pre) - (ct_post - ct_pre)
            if not np.isnan(boot_did):
                bootstrap_estimates.append(boot_did)
    
    if len(bootstrap_estimates) > 0:
        std_error = np.std(bootstrap_estimates)
    else:
        std_error = np.nan
    
    return {
        'cohort': int(cohort_year),
        'att': float(did_estimate),
        'std_error': float(std_error),
        'n_treated_regions': int(cohort['region_code'].nunique()),
        'n_control_regions': int(control['region_code'].nunique()),
        'cohort_pre_mean': float(cohort_pre),
        'cohort_post_mean': float(cohort_post),
        'control_pre_mean': float(control_pre),
        'control_post_mean': float(control_post),
    }


def run_callaway_santanna_did(df, outcome='deaths_all'):
    """
    Run Callaway-Sant'Anna style aggregated DiD.
    
    Steps:
    1. Calculate cohort-specific ATT for each treatment cohort
    2. Aggregate using inverse variance weighting
    
    Note: This is a simplified implementation. For production, use did package.
    """
    print(f"\n--- Callaway-Sant'Anna DiD: LPT → {outcome} ---")
    
    data = df.copy()
    
    # Get treatment cohorts
    if 'treatment_cohort' not in data.columns:
        print("  Treatment cohort not found in data.")
        return None
    
    cohorts = data[data['treatment_cohort'].notna()]['treatment_cohort'].unique()
    cohorts = sorted([c for c in cohorts if c >= 2005 and c <= 2020])  # Focus on main rollout
    
    print(f"  Treatment cohorts: {len(cohorts)}")
    
    # Identify never-treated (or late-treated as control)
    never_treated = data[data['first_treatment_year'].isna() | (data['first_treatment_year'] > 2022)]
    print(f"  Never/late-treated regions: {never_treated['region_code'].nunique()}")
    
    # Calculate cohort-specific ATTs
    cohort_results = []
    
    for cohort_year in cohorts:
        result = calculate_cohort_att(data, cohort_year, outcome, never_treated)
        if result is not None:
            cohort_results.append(result)
            print(f"    Cohort {cohort_year}: ATT = {result['att']:.4f} "
                  f"(n_treated = {result['n_treated_regions']})")
    
    if len(cohort_results) == 0:
        print("  No valid cohort results.")
        return None
    
    # Aggregate ATTs (simple average for now)
    # For proper implementation, use inverse variance weighting
    att_values = [r['att'] for r in cohort_results]
    se_values = [r['std_error'] for r in cohort_results if not np.isnan(r['std_error'])]
    weights = [r['n_treated_regions'] for r in cohort_results]
    
    # Weighted average
    total_weight = sum(weights)
    aggregate_att = sum(att * w for att, w in zip(att_values, weights)) / total_weight
    
    # Simple SE (not accounting for correlation)
    if len(se_values) > 0:
        aggregate_se = np.sqrt(sum(se**2 * (w/total_weight)**2 
                                   for se, w in zip(se_values, weights)))
    else:
        aggregate_se = np.std(att_values)
    
    results = {
        'method': 'Callaway-SantAnna',
        'outcome': outcome,
        'aggregate_att': float(aggregate_att),
        'aggregate_se': float(aggregate_se),
        't_stat': float(aggregate_att / aggregate_se) if aggregate_se > 0 else np.nan,
        'p_value': float(2 * (1 - stats.norm.cdf(abs(aggregate_att / aggregate_se)))) if aggregate_se > 0 else np.nan,
        'n_cohorts': len(cohort_results),
        'cohort_results': cohort_results,
    }
    
    print(f"\n  Aggregate ATT: {results['aggregate_att']:.4f} (SE: {results['aggregate_se']:.4f})")
    
    return results


# =============================================================================
# TRIPLE DIFFERENCE
# =============================================================================

def run_triple_difference(df, outcome='deaths_all'):
    """
    Triple Difference: Electricity × Hot Days × Post-Treatment
    
    This isolates the effect of electricity on adaptation to heat.
    
    Model: Y = β1*Post + β2*Hot + β3*Post×Hot + 
           β4*Treat + β5*Treat×Post + β6*Treat×Hot + 
           β7*Treat×Post×Hot + controls + FE
    
    β7 is the triple difference: effect of treatment on hot days after treatment.
    """
    print(f"\n--- Triple Difference: Electricity × Heat × Post ---")
    
    data = df.dropna(subset=['post_lpt', 'extreme_heat', outcome]).copy()
    
    if 'post_lpt' not in data.columns or 'extreme_heat' not in data.columns:
        print("  Required variables not found.")
        return None
    
    # Create interactions
    data['treated'] = (data['first_treatment_year'].notna()).astype(int)
    data['post'] = data['post_lpt']
    data['hot'] = data['extreme_heat']
    
    data['post_hot'] = data['post'] * data['hot']
    data['treated_post'] = data['treated'] * data['post']
    data['treated_hot'] = data['treated'] * data['hot']
    data['triple'] = data['treated'] * data['post'] * data['hot']
    
    # Controls
    controls = []
    for m in range(2, 13):
        data[f'month_{m}'] = (data['month'] == m).astype(int)
        controls.append(f'month_{m}')
    
    # Demean by region
    vars_to_demean = [outcome, 'post', 'hot', 'post_hot', 'treated_post', 
                      'treated_hot', 'triple'] + controls
    
    for col in vars_to_demean:
        if col in data.columns:
            data[f'{col}_dm'] = data.groupby('region_code')[col].transform(
                lambda x: x - x.mean()
            )
    
    # Run regression
    y = data[f'{outcome}_dm']
    X_cols = ['post_dm', 'hot_dm', 'post_hot_dm', 'treated_post_dm', 
              'treated_hot_dm', 'triple_dm']
    X_cols += [f'{c}_dm' for c in controls]
    X_cols = [c for c in X_cols if c in data.columns]
    X = data[X_cols]
    X = sm.add_constant(X)
    
    model = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    results = {
        'method': 'Triple Difference',
        'outcome': outcome,
        'n_obs': int(model.nobs),
        'triple_diff_estimate': float(model.params['triple_dm']),
        'std_error': float(model.bse['triple_dm']),
        't_stat': float(model.tvalues['triple_dm']),
        'p_value': float(model.pvalues['triple_dm']),
        'ci_lower': float(model.conf_int().loc['triple_dm', 0]),
        'ci_upper': float(model.conf_int().loc['triple_dm', 1]),
        'coefficients': {
            'post': float(model.params.get('post_dm', np.nan)),
            'hot': float(model.params.get('hot_dm', np.nan)),
            'post_hot': float(model.params.get('post_hot_dm', np.nan)),
            'treated_post': float(model.params.get('treated_post_dm', np.nan)),
            'treated_hot': float(model.params.get('treated_hot_dm', np.nan)),
            'triple': float(model.params['triple_dm']),
        }
    }
    
    print(f"  Triple Diff (Treat×Post×Hot): {results['triple_diff_estimate']:.4f}")
    print(f"  SE: {results['std_error']:.4f}, p-value: {results['p_value']:.4f}")
    print(f"  Interpretation: Effect of electricity on hot-day mortality reduction")
    
    return results


# =============================================================================
# INTENSITY-BASED ANALYSIS
# =============================================================================

def run_intensity_analysis(df, outcome='deaths_all'):
    """
    Use treatment intensity (households electrified) instead of binary treatment.
    
    This provides dose-response evidence.
    """
    print(f"\n--- Intensity Analysis: Log(Households) → {outcome} ---")
    
    data = df.copy()
    
    if 'log_households' not in data.columns:
        print("  log_households not found.")
        return None
    
    # Replace NaN with 0 for never-treated
    data['log_households'] = data['log_households'].fillna(0)
    data['intensity'] = data['log_households'] * data['post_lpt'].fillna(0)
    
    # Demean
    for col in [outcome, 'intensity']:
        data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    # Controls
    controls = []
    for m in range(2, 13):
        data[f'month_{m}'] = (data['month'] == m).astype(int)
        controls.append(f'month_{m}')
    
    for col in controls:
        data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    # Regression
    y = data[f'{outcome}_dm']
    X_cols = ['intensity_dm'] + [f'{c}_dm' for c in controls]
    X = data[X_cols]
    X = sm.add_constant(X)
    
    model = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    results = {
        'method': 'Intensity',
        'outcome': outcome,
        'n_obs': int(model.nobs),
        'coefficient': float(model.params['intensity_dm']),
        'std_error': float(model.bse['intensity_dm']),
        't_stat': float(model.tvalues['intensity_dm']),
        'p_value': float(model.pvalues['intensity_dm']),
        'interpretation': 'Effect of 1 log-unit increase in households electrified'
    }
    
    print(f"  Coefficient: {results['coefficient']:.4f} (SE: {results['std_error']:.4f})")
    print(f"  p-value: {results['p_value']:.4f}")
    
    return results


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_full_did_analysis(level='intermediate'):
    """Run complete DiD analysis."""
    
    print("\n" + "=" * 70)
    print(f"LUZ PARA TODOS DiD ANALYSIS - {level.upper()} LEVEL")
    print("=" * 70)
    
    # Load data
    df = load_panel(level)
    
    # Check for treatment variable
    if 'post_lpt' not in df.columns:
        print("ERROR: post_lpt not found. Run data preparation scripts first.")
        return None
    
    all_results = {
        'level': level,
        'twfe': {},
        'callaway_santanna': {},
        'triple_difference': {},
        'intensity': {}
    }
    
    outcomes = ['deaths_all', 'deaths_cardiovascular', 'deaths_respiratory']
    
    # TWFE DiD
    print("\n" + "=" * 50)
    print("1. TWO-WAY FIXED EFFECTS DiD")
    print("=" * 50)
    
    for outcome in outcomes:
        results = run_twfe_did(df, outcome)
        if results:
            all_results['twfe'][outcome] = results
    
    # Callaway-Sant'Anna
    print("\n" + "=" * 50)
    print("2. CALLAWAY-SANT'ANNA DiD")
    print("=" * 50)
    
    for outcome in outcomes:
        results = run_callaway_santanna_did(df, outcome)
        if results:
            all_results['callaway_santanna'][outcome] = results
    
    # Triple Difference
    print("\n" + "=" * 50)
    print("3. TRIPLE DIFFERENCE")
    print("=" * 50)
    
    for outcome in outcomes:
        results = run_triple_difference(df, outcome)
        if results:
            all_results['triple_difference'][outcome] = results
    
    # Intensity Analysis
    print("\n" + "=" * 50)
    print("4. INTENSITY ANALYSIS")
    print("=" * 50)
    
    for outcome in outcomes:
        results = run_intensity_analysis(df, outcome)
        if results:
            all_results['intensity'][outcome] = results
    
    return all_results


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='LPT Staggered DiD Analysis')
    parser.add_argument('--level', type=str, default='both',
                        choices=['intermediate', 'immediate', 'both'],
                        help='Spatial level to analyze')
    args = parser.parse_args()
    
    print("=" * 70)
    print("LUZ PARA TODOS STAGGERED DIFFERENCE-IN-DIFFERENCES")
    print("=" * 70)
    
    levels = ['intermediate', 'immediate'] if args.level == 'both' else [args.level]
    
    for level in levels:
        results = run_full_did_analysis(level)
        
        if results:
            # Save results
            suffix = '' if level == 'intermediate' else '_immediate'
            output_file = OUTPUT_DIR / f"lpt_did_results{suffix}.json"
            
            with open(output_file, 'w') as f:
                json.dump(convert_to_json_serializable(results), f, indent=2)
            
            print(f"\nResults saved to: {output_file}")
            
            # Print summary
            print("\n" + "-" * 50)
            print(f"SUMMARY - {level.upper()}")
            print("-" * 50)
            
            if results['twfe'].get('deaths_all'):
                twfe = results['twfe']['deaths_all']
                print(f"TWFE ATT: {twfe['att']:.4f} (p={twfe['p_value']:.4f})")
            
            if results['triple_difference'].get('deaths_all'):
                td = results['triple_difference']['deaths_all']
                print(f"Triple Diff: {td['triple_diff_estimate']:.4f} (p={td['p_value']:.4f})")
    
    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
