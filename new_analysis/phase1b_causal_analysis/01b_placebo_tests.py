#!/usr/bin/env python3
"""
Placebo and Falsification Tests
================================

Comprehensive battery of tests to validate causal identification:

1. Placebo outcomes - Effects on outcomes that shouldn't be affected
2. Placebo timing - Fake treatment dates before actual treatment
3. Permutation tests - Random treatment assignment
4. Exclusion restriction tests - ENSO → direct mortality pathways
5. Balance tests - Pre-treatment characteristics

These tests help establish internal validity of the quasi-experimental designs.
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
    elif isinstance(obj, (np.bool_, bool)):
        return bool(obj)
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


def load_panel(level='intermediate'):
    """Load causal analysis panel."""
    print(f"Loading {level} panel...")
    df = pd.read_parquet(RESULTS_DIR / f"causal_panel_{level}.parquet")
    print(f"  Shape: {df.shape}")
    return df


# =============================================================================
# PLACEBO OUTCOME TESTS
# =============================================================================

def run_placebo_outcome_test(df, placebo_outcomes=['deaths_cardiovascular']):
    """
    Test for effects on placebo outcomes.
    
    For ENSO: Temperature should affect heat-related deaths more than
    non-temperature-related causes.
    
    For LPT: Electrification should affect heat deaths but not unrelated causes.
    """
    print("\n" + "=" * 50)
    print("PLACEBO OUTCOME TESTS")
    print("=" * 50)
    
    results = {}
    
    # Test 1: ENSO effect on different cause-specific mortality
    print("\n--- ENSO effects by cause of death ---")
    
    if 'oni' in df.columns:
        data = df.dropna(subset=['oni']).copy()
        
        outcomes = ['deaths_all', 'deaths_cardiovascular', 'deaths_respiratory']
        
        for outcome in outcomes:
            if outcome not in data.columns:
                continue
            
            for col in ['oni', outcome]:
                data[f'{col}_dm'] = data.groupby('region_code')[col].transform(
                    lambda x: x - x.mean()
                )
            
            y = data[f'{outcome}_dm']
            X = sm.add_constant(data['oni_dm'])
            
            try:
                model = sm.OLS(y, X).fit(cov_type='cluster', 
                                          cov_kwds={'groups': data['region_code']})
                
                results[f'enso_{outcome}'] = {
                    'coefficient': float(model.params['oni_dm']),
                    'std_error': float(model.bse['oni_dm']),
                    'p_value': float(model.pvalues['oni_dm']),
                    'n_obs': int(model.nobs)
                }
                
                print(f"  {outcome}: β = {results[f'enso_{outcome}']['coefficient']:.4f} "
                      f"(p = {results[f'enso_{outcome}']['p_value']:.4f})")
            except Exception as e:
                print(f"  Error for {outcome}: {str(e)}")
    
    # Test 2: LPT effect on different causes
    print("\n--- LPT effects by cause of death ---")
    
    if 'post_lpt' in df.columns:
        data = df.dropna(subset=['post_lpt']).copy()
        
        for outcome in outcomes:
            if outcome not in data.columns:
                continue
            
            for col in ['post_lpt', outcome]:
                data[f'{col}_dm'] = data.groupby('region_code')[col].transform(
                    lambda x: x - x.mean()
                )
            
            y = data[f'{outcome}_dm']
            X = sm.add_constant(data['post_lpt_dm'])
            
            try:
                model = sm.OLS(y, X).fit(cov_type='cluster',
                                          cov_kwds={'groups': data['region_code']})
                
                results[f'lpt_{outcome}'] = {
                    'coefficient': float(model.params['post_lpt_dm']),
                    'std_error': float(model.bse['post_lpt_dm']),
                    'p_value': float(model.pvalues['post_lpt_dm']),
                    'n_obs': int(model.nobs)
                }
                
                print(f"  {outcome}: β = {results[f'lpt_{outcome}']['coefficient']:.4f} "
                      f"(p = {results[f'lpt_{outcome}']['p_value']:.4f})")
            except Exception as e:
                print(f"  Error for {outcome}: {str(e)}")
    
    return results


# =============================================================================
# PLACEBO TIMING TESTS
# =============================================================================

def run_placebo_timing_test(df, outcome='deaths_all', fake_leads=[2, 3, 4]):
    """
    Test for effects at fake treatment times before actual treatment.
    
    If we see effects before treatment, something is wrong with the design.
    """
    print("\n" + "=" * 50)
    print("PLACEBO TIMING TESTS")
    print("=" * 50)
    
    results = {}
    
    data = df.copy()
    
    if 'first_treatment_year' not in data.columns:
        print("  first_treatment_year not found")
        return results
    
    for lead in fake_leads:
        print(f"\n--- Fake treatment {lead} years before actual ---")
        
        # Create fake treatment
        data['fake_treatment_year'] = data['first_treatment_year'] - lead
        data['fake_post'] = (data['year'] >= data['fake_treatment_year']).astype(int)
        data.loc[data['first_treatment_year'].isna(), 'fake_post'] = 0
        
        # Only use pre-treatment period for treated units
        pre_treatment = data[
            (data['year'] < data['first_treatment_year']) | 
            (data['first_treatment_year'].isna())
        ].copy()
        
        if len(pre_treatment) == 0:
            continue
        
        # Demean
        for col in [outcome, 'fake_post']:
            pre_treatment[f'{col}_dm'] = pre_treatment.groupby('region_code')[col].transform(
                lambda x: x - x.mean()
            )
        
        y = pre_treatment[f'{outcome}_dm']
        X = sm.add_constant(pre_treatment['fake_post_dm'])
        
        try:
            model = sm.OLS(y, X).fit(cov_type='cluster',
                                      cov_kwds={'groups': pre_treatment['region_code']})
            
            results[f'lead_{lead}'] = {
                'coefficient': float(model.params['fake_post_dm']),
                'std_error': float(model.bse['fake_post_dm']),
                'p_value': float(model.pvalues['fake_post_dm']),
                'n_obs': int(model.nobs),
                'significant': model.pvalues['fake_post_dm'] < 0.05
            }
            
            print(f"  Lead {lead}: β = {results[f'lead_{lead}']['coefficient']:.4f} "
                  f"(p = {results[f'lead_{lead}']['p_value']:.4f})")
            
            if results[f'lead_{lead}']['significant']:
                print(f"  ⚠️ WARNING: Significant effect at fake timing!")
        except Exception as e:
            print(f"  Error: {str(e)}")
    
    # Summary
    any_significant = any(r.get('significant', False) for r in results.values())
    results['summary'] = {
        'any_significant_placebo': any_significant,
        'interpretation': 'Placebo test FAILED - pre-trends detected' if any_significant 
                         else 'Placebo test PASSED - no pre-trends'
    }
    
    return results


# =============================================================================
# PERMUTATION TESTS
# =============================================================================

def run_permutation_test(df, outcome='deaths_all', n_permutations=100):
    """
    Permutation test: randomly reassign treatment across regions.
    
    If the effect is truly causal, randomly assigned treatment should
    show no effect on average.
    """
    print("\n" + "=" * 50)
    print(f"PERMUTATION TEST (n={n_permutations})")
    print("=" * 50)
    
    data = df.copy()
    
    if 'post_lpt' not in data.columns:
        print("  post_lpt not found")
        return None
    
    # Get actual effect first
    print("\n--- Calculating actual effect ---")
    
    data_clean = data.dropna(subset=['post_lpt', outcome]).copy()
    
    for col in [outcome, 'post_lpt']:
        data_clean[f'{col}_dm'] = data_clean.groupby('region_code')[col].transform(
            lambda x: x - x.mean()
        )
    
    y = data_clean[f'{outcome}_dm']
    X = sm.add_constant(data_clean['post_lpt_dm'])
    
    actual_model = sm.OLS(y, X).fit()
    actual_effect = actual_model.params['post_lpt_dm']
    
    print(f"  Actual effect: {actual_effect:.4f}")
    
    # Run permutations
    print(f"\n--- Running {n_permutations} permutations ---")
    
    regions = data_clean['region_code'].unique()
    
    # Get treatment status by region
    region_treatment = data_clean.groupby('region_code')['post_lpt'].max().to_dict()
    treated_regions = [r for r, t in region_treatment.items() if t > 0]
    n_treated = len(treated_regions)
    
    permutation_effects = []
    
    for i in range(n_permutations):
        # Randomly assign same number of treated regions
        fake_treated = np.random.choice(regions, size=n_treated, replace=False)
        
        # Create fake treatment
        data_clean['fake_treated'] = data_clean['region_code'].isin(fake_treated).astype(int)
        
        # Demean
        data_clean['fake_dm'] = data_clean.groupby('region_code')['fake_treated'].transform(
            lambda x: x - x.mean()
        )
        
        # Estimate
        try:
            y = data_clean[f'{outcome}_dm']
            X = sm.add_constant(data_clean['fake_dm'])
            model = sm.OLS(y, X).fit()
            permutation_effects.append(model.params['fake_dm'])
        except:
            pass
        
        if (i + 1) % 20 == 0:
            print(f"    Completed {i + 1}/{n_permutations}")
    
    # Calculate p-value
    permutation_effects = np.array(permutation_effects)
    
    # Two-sided test
    p_value = np.mean(np.abs(permutation_effects) >= np.abs(actual_effect))
    
    results = {
        'actual_effect': float(actual_effect),
        'n_permutations': n_permutations,
        'mean_permutation_effect': float(np.mean(permutation_effects)),
        'std_permutation_effect': float(np.std(permutation_effects)),
        'p_value': float(p_value),
        'percentile_5': float(np.percentile(permutation_effects, 5)),
        'percentile_95': float(np.percentile(permutation_effects, 95)),
        'interpretation': 'Effect is significant beyond random chance' if p_value < 0.05
                         else 'Effect could be due to chance'
    }
    
    print(f"\n  Permutation test p-value: {p_value:.4f}")
    print(f"  {results['interpretation']}")
    
    return results


# =============================================================================
# EXCLUSION RESTRICTION TESTS (FOR IV)
# =============================================================================

def run_exclusion_restriction_test(df, outcome='deaths_all'):
    """
    Test exclusion restriction for ENSO as instrument.
    
    ENSO should affect mortality ONLY through temperature.
    We test by controlling for temperature - residual ENSO effect should be zero.
    """
    print("\n" + "=" * 50)
    print("EXCLUSION RESTRICTION TEST (ENSO IV)")
    print("=" * 50)
    
    data = df.copy()
    
    required = ['oni', 'temp_mean', outcome]
    if not all(v in data.columns for v in required):
        print("  Required variables not found")
        return None
    
    data = data.dropna(subset=required).copy()
    
    results = {}
    
    # Model 1: ENSO → Mortality (reduced form)
    print("\n--- Model 1: ENSO → Mortality (reduced form) ---")
    
    for col in ['oni', outcome]:
        data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    y = data[f'{outcome}_dm']
    X = sm.add_constant(data['oni_dm'])
    
    model1 = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    results['reduced_form'] = {
        'enso_effect': float(model1.params['oni_dm']),
        'std_error': float(model1.bse['oni_dm']),
        'p_value': float(model1.pvalues['oni_dm'])
    }
    
    print(f"  ENSO → Mortality: {results['reduced_form']['enso_effect']:.4f} "
          f"(p = {results['reduced_form']['p_value']:.4f})")
    
    # Model 2: ENSO → Mortality | Temperature (residual effect)
    print("\n--- Model 2: ENSO → Mortality | Temperature ---")
    
    for col in ['oni', 'temp_mean', outcome]:
        data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    y = data[f'{outcome}_dm']
    X = sm.add_constant(data[['oni_dm', 'temp_mean_dm']])
    
    model2 = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    results['conditional'] = {
        'enso_effect_conditional': float(model2.params['oni_dm']),
        'std_error': float(model2.bse['oni_dm']),
        'p_value': float(model2.pvalues['oni_dm']),
        'temp_effect': float(model2.params['temp_mean_dm']),
    }
    
    print(f"  ENSO → Mortality | Temp: {results['conditional']['enso_effect_conditional']:.4f} "
          f"(p = {results['conditional']['p_value']:.4f})")
    print(f"  Temperature → Mortality: {results['conditional']['temp_effect']:.4f}")
    
    # Interpretation
    reduced_effect = abs(results['reduced_form']['enso_effect'])
    conditional_effect = abs(results['conditional']['enso_effect_conditional'])
    
    if conditional_effect < reduced_effect * 0.5:
        interpretation = ("Exclusion restriction SUPPORTED: ENSO effect reduced by >50% "
                         "when controlling for temperature")
    elif results['conditional']['p_value'] > 0.1:
        interpretation = ("Exclusion restriction SUPPORTED: ENSO effect is insignificant "
                         "conditional on temperature")
    else:
        interpretation = ("Exclusion restriction QUESTIONABLE: ENSO has residual effect "
                         "not explained by temperature")
    
    results['interpretation'] = interpretation
    print(f"\n  {interpretation}")
    
    return results


# =============================================================================
# BALANCE TESTS
# =============================================================================

def run_balance_test(df):
    """
    Test for balance in pre-treatment characteristics between treated and control.
    
    Treated and never-treated regions should be similar before LPT.
    """
    print("\n" + "=" * 50)
    print("BALANCE TESTS")
    print("=" * 50)
    
    data = df.copy()
    
    if 'first_treatment_year' not in data.columns:
        print("  first_treatment_year not found")
        return None
    
    # Define treatment groups
    data['ever_treated'] = data['first_treatment_year'].notna().astype(int)
    
    # Use pre-treatment period (before 2005, main LPT rollout)
    pre_treatment = data[data['year'] < 2005].copy()
    
    if len(pre_treatment) == 0:
        print("  No pre-treatment data available")
        return None
    
    # Variables to test
    balance_vars = ['temp_mean', 'deaths_all', 'deaths_cardiovascular', 'deaths_respiratory']
    balance_vars = [v for v in balance_vars if v in pre_treatment.columns]
    
    results = {}
    
    for var in balance_vars:
        treated = pre_treatment[pre_treatment['ever_treated'] == 1].groupby('region_code')[var].mean()
        control = pre_treatment[pre_treatment['ever_treated'] == 0].groupby('region_code')[var].mean()
        
        # T-test
        t_stat, p_value = stats.ttest_ind(treated, control, nan_policy='omit')
        
        results[var] = {
            'treated_mean': float(treated.mean()),
            'control_mean': float(control.mean()),
            'difference': float(treated.mean() - control.mean()),
            't_stat': float(t_stat),
            'p_value': float(p_value),
            'balanced': p_value > 0.1
        }
        
        status = '✓ Balanced' if results[var]['balanced'] else '⚠️ Imbalanced'
        print(f"  {var}: Treated={results[var]['treated_mean']:.2f}, "
              f"Control={results[var]['control_mean']:.2f}, p={p_value:.4f} {status}")
    
    # Summary
    all_balanced = all(r['balanced'] for r in results.values())
    results['summary'] = {
        'all_balanced': all_balanced,
        'interpretation': 'Pre-treatment balance SATISFIED' if all_balanced
                         else 'Pre-treatment balance VIOLATED - consider matching'
    }
    
    print(f"\n  {results['summary']['interpretation']}")
    
    return results


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_full_placebo_analysis(level='intermediate'):
    """Run complete placebo and falsification tests."""
    
    print("\n" + "=" * 70)
    print(f"PLACEBO AND FALSIFICATION TESTS - {level.upper()} LEVEL")
    print("=" * 70)
    
    # Load data
    df = load_panel(level)
    
    all_results = {
        'level': level,
        'placebo_outcomes': {},
        'placebo_timing': {},
        'permutation': {},
        'exclusion_restriction': {},
        'balance': {}
    }
    
    # Placebo outcome tests
    try:
        results = run_placebo_outcome_test(df)
        all_results['placebo_outcomes'] = results
    except Exception as e:
        print(f"  Error in placebo outcome test: {str(e)}")
    
    # Placebo timing tests
    try:
        results = run_placebo_timing_test(df)
        all_results['placebo_timing'] = results
    except Exception as e:
        print(f"  Error in placebo timing test: {str(e)}")
    
    # Permutation test (reduced iterations for speed)
    try:
        results = run_permutation_test(df, n_permutations=50)
        all_results['permutation'] = results
    except Exception as e:
        print(f"  Error in permutation test: {str(e)}")
    
    # Exclusion restriction test
    try:
        results = run_exclusion_restriction_test(df)
        all_results['exclusion_restriction'] = results
    except Exception as e:
        print(f"  Error in exclusion restriction test: {str(e)}")
    
    # Balance tests
    try:
        results = run_balance_test(df)
        all_results['balance'] = results
    except Exception as e:
        print(f"  Error in balance test: {str(e)}")
    
    return all_results


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Placebo and Falsification Tests')
    parser.add_argument('--level', type=str, default='both',
                        choices=['intermediate', 'immediate', 'both'],
                        help='Spatial level to analyze')
    args = parser.parse_args()
    
    print("=" * 70)
    print("PLACEBO AND FALSIFICATION TESTS")
    print("=" * 70)
    
    levels = ['intermediate', 'immediate'] if args.level == 'both' else [args.level]
    
    for level in levels:
        results = run_full_placebo_analysis(level)
        
        # Save results
        suffix = '' if level == 'intermediate' else '_immediate'
        output_file = OUTPUT_DIR / f"placebo_tests_results{suffix}.json"
        
        with open(output_file, 'w') as f:
            json.dump(convert_to_json_serializable(results), f, indent=2)
        
        print(f"\nResults saved to: {output_file}")
        
        # Summary
        print("\n" + "=" * 50)
        print(f"VALIDITY SUMMARY - {level.upper()}")
        print("=" * 50)
        
        if results['placebo_timing'] and results['placebo_timing'].get('summary'):
            print(f"  Placebo timing: {results['placebo_timing']['summary']['interpretation']}")
        
        if results['exclusion_restriction'] and results['exclusion_restriction'].get('interpretation'):
            print(f"  Exclusion restriction: {results['exclusion_restriction']['interpretation'][:60]}...")
        
        if results['balance'] and results['balance'].get('summary'):
            print(f"  Balance: {results['balance']['summary']['interpretation']}")
    
    print("\n" + "=" * 70)
    print("All tests complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
