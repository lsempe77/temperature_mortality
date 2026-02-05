#!/usr/bin/env python3
"""
Luz para Todos Event Study Analysis
====================================

Dynamic treatment effects analysis for the Luz para Todos program.
Estimates pre-trends and post-treatment effects over time.

Methods:
1. Standard event study with time relative to treatment
2. Heterogeneity by treatment intensity
3. Pre-trend tests for parallel trends assumption
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


def load_panel(level='intermediate'):
    """Load causal analysis panel."""
    print(f"Loading {level} panel...")
    df = pd.read_parquet(RESULTS_DIR / f"causal_panel_{level}.parquet")
    print(f"  Shape: {df.shape}")
    return df


# =============================================================================
# EVENT STUDY SETUP
# =============================================================================

def create_event_time(df):
    """
    Create event time variable (years relative to treatment).
    
    Event time k:
    - k < 0: k years before treatment
    - k = 0: year of treatment
    - k > 0: k years after treatment
    """
    print("\nCreating event time indicators...")
    
    data = df.copy()
    
    # Calculate years since treatment
    if 'first_treatment_year' not in data.columns:
        print("  first_treatment_year not found")
        return data
    
    data['event_year'] = data['year'] - data['first_treatment_year']
    
    # For never-treated, set to NaN (they will be controls)
    data.loc[data['first_treatment_year'].isna(), 'event_year'] = np.nan
    
    # Cap event time to avoid sparse cells
    data['event_year_capped'] = data['event_year'].clip(-5, 10)
    
    print(f"  Event time range: {data['event_year'].min():.0f} to {data['event_year'].max():.0f}")
    print(f"  Capped range: -5 to 10")
    
    return data


# =============================================================================
# EVENT STUDY ESTIMATION
# =============================================================================

def run_event_study(df, outcome='deaths_all', reference_period=-1, 
                   pre_periods=5, post_periods=10):
    """
    Run event study regression.
    
    Model: Y_{rt} = Σ_{k≠-1} β_k * 1{EventTime=k} + X_{rt}γ + α_r + δ_t + ε_{rt}
    
    Args:
        df: Panel data with event_year
        outcome: Dependent variable
        reference_period: Period to normalize to (usually -1)
        pre_periods: Number of pre-treatment periods
        post_periods: Number of post-treatment periods
    """
    print(f"\n--- Event Study: LPT → {outcome} ---")
    
    # Prepare data
    data = df.dropna(subset=[outcome]).copy()
    
    # Need event time for treated units
    treated = data['first_treatment_year'].notna()
    never_treated = ~treated
    
    print(f"  Treated observations: {treated.sum():,}")
    print(f"  Never-treated observations: {never_treated.sum():,}")
    
    # Create event time dummies
    event_times = list(range(-pre_periods, post_periods + 1))
    event_times.remove(reference_period)
    
    for t in event_times:
        data[f'event_{t}'] = ((data['event_year_capped'] == t) & treated).astype(int)
    
    # Add controls
    controls = []
    for m in range(2, 13):
        data[f'month_{m}'] = (data['month'] == m).astype(int)
        controls.append(f'month_{m}')
    
    if 'dow' in data.columns:
        for d in range(1, 7):
            data[f'dow_{d}'] = (data['dow'] == d).astype(int)
            controls.append(f'dow_{d}')
    
    # Demean by region (within estimator)
    vars_to_demean = [outcome] + [f'event_{t}' for t in event_times] + controls
    
    for col in vars_to_demean:
        if col in data.columns:
            data[f'{col}_dm'] = data.groupby('region_code')[col].transform(
                lambda x: x - x.mean()
            )
    
    # Run regression
    y = data[f'{outcome}_dm']
    X_cols = [f'event_{t}_dm' for t in event_times]
    X_cols += [f'{c}_dm' for c in controls]
    X = data[X_cols]
    X = sm.add_constant(X)
    
    model = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    # Extract coefficients
    event_coefficients = {}
    for t in event_times:
        col = f'event_{t}_dm'
        event_coefficients[t] = {
            'coefficient': float(model.params[col]),
            'std_error': float(model.bse[col]),
            'ci_lower': float(model.conf_int().loc[col, 0]),
            'ci_upper': float(model.conf_int().loc[col, 1]),
            'p_value': float(model.pvalues[col])
        }
    
    # Add reference period
    event_coefficients[reference_period] = {
        'coefficient': 0.0,
        'std_error': 0.0,
        'ci_lower': 0.0,
        'ci_upper': 0.0,
        'p_value': 1.0
    }
    
    # Pre-trend tests
    pre_coefs = [event_coefficients[t]['coefficient'] for t in event_times if t < 0]
    pre_pvals = [event_coefficients[t]['p_value'] for t in event_times if t < 0]
    
    # Joint F-test for pre-trends
    pre_vars = [f'event_{t}_dm' for t in event_times if t < 0]
    r_matrix = np.zeros((len(pre_vars), len(model.params)))
    for i, var in enumerate(pre_vars):
        r_matrix[i, list(model.params.index).index(var)] = 1
    
    try:
        f_test = model.f_test(r_matrix)
        pre_trend_f = float(f_test.fvalue)
        pre_trend_p = float(f_test.pvalue)
    except:
        pre_trend_f = np.nan
        pre_trend_p = np.nan
    
    # Post-treatment average
    post_coefs = [event_coefficients[t]['coefficient'] for t in event_times if t >= 0]
    avg_post_effect = np.mean(post_coefs)
    
    results = {
        'outcome': outcome,
        'reference_period': reference_period,
        'pre_periods': pre_periods,
        'post_periods': post_periods,
        'n_obs': int(model.nobs),
        'n_regions': int(data['region_code'].nunique()),
        'r_squared': float(model.rsquared),
        'event_coefficients': event_coefficients,
        'pre_trend_test': {
            'f_stat': pre_trend_f,
            'p_value': pre_trend_p,
            'any_significant': any(p < 0.05 for p in pre_pvals),
            'interpretation': 'Parallel trends supported' if pre_trend_p > 0.05 
                             else 'Pre-trends detected - caution'
        },
        'avg_post_effect': float(avg_post_effect),
    }
    
    print(f"  Observations: {results['n_obs']:,}")
    print(f"  Pre-trend F-test p-value: {pre_trend_p:.4f}")
    print(f"  Pre-trends: {results['pre_trend_test']['interpretation']}")
    print(f"  Average post-treatment effect: {avg_post_effect:.4f}")
    
    return results


def run_event_study_by_intensity(df, outcome='deaths_all'):
    """
    Event study stratified by treatment intensity.
    
    Compares high vs low intensity treatment.
    """
    print(f"\n--- Event Study by Intensity: LPT → {outcome} ---")
    
    data = df.copy()
    
    if 'log_households' not in data.columns:
        print("  log_households not found")
        return None
    
    # Calculate intensity terciles among treated
    treated_intensity = data[data['first_treatment_year'].notna()].groupby('region_code')['log_households'].first()
    
    q33 = treated_intensity.quantile(0.33)
    q67 = treated_intensity.quantile(0.67)
    
    intensity_map = {}
    for region, intensity in treated_intensity.items():
        if intensity <= q33:
            intensity_map[region] = 'low'
        elif intensity <= q67:
            intensity_map[region] = 'medium'
        else:
            intensity_map[region] = 'high'
    
    data['intensity_group'] = data['region_code'].map(intensity_map)
    
    results = {}
    
    for intensity in ['low', 'medium', 'high']:
        print(f"\n  Intensity group: {intensity}")
        
        # Filter to this intensity group and never-treated
        subset = data[
            (data['intensity_group'] == intensity) | 
            (data['first_treatment_year'].isna())
        ].copy()
        
        if len(subset) > 0:
            try:
                es_results = run_event_study(subset, outcome)
                results[intensity] = es_results
            except Exception as e:
                print(f"    Error: {str(e)}")
    
    return results


def run_event_study_by_heat_exposure(df, outcome='deaths_all'):
    """
    Event study stratified by baseline heat exposure.
    
    Tests if electrification has larger effects in hotter regions.
    """
    print(f"\n--- Event Study by Heat Exposure: LPT → {outcome} ---")
    
    data = df.copy()
    
    # Calculate baseline heat exposure (pre-treatment average extreme heat days)
    if 'extreme_heat' not in data.columns:
        print("  extreme_heat not found")
        return None
    
    # Use pre-2005 (before main LPT rollout) as baseline
    baseline = data[data['year'] < 2005].groupby('region_code')['extreme_heat'].mean()
    
    median_heat = baseline.median()
    heat_map = {r: 'high' if h >= median_heat else 'low' for r, h in baseline.items()}
    
    data['heat_group'] = data['region_code'].map(heat_map)
    
    results = {}
    
    for heat in ['low', 'high']:
        print(f"\n  Heat exposure group: {heat}")
        
        subset = data[data['heat_group'] == heat].copy()
        
        if len(subset) > 0:
            try:
                es_results = run_event_study(subset, outcome)
                results[heat] = es_results
            except Exception as e:
                print(f"    Error: {str(e)}")
    
    return results


# =============================================================================
# DYNAMIC EFFECTS ON HOT DAYS
# =============================================================================

def run_hot_day_event_study(df, outcome='deaths_all'):
    """
    Event study separately for hot days vs normal days.
    
    If electrification enables AC access, effects should be larger on hot days.
    """
    print(f"\n--- Hot Day vs Normal Day Event Study ---")
    
    data = df.copy()
    
    if 'extreme_heat' not in data.columns:
        print("  extreme_heat not found")
        return None
    
    results = {}
    
    for day_type in ['hot', 'normal']:
        print(f"\n  Day type: {day_type}")
        
        if day_type == 'hot':
            subset = data[data['extreme_heat'] == 1].copy()
        else:
            subset = data[data['extreme_heat'] == 0].copy()
        
        if len(subset) > 0:
            try:
                es_results = run_event_study(subset, outcome)
                results[day_type] = es_results
            except Exception as e:
                print(f"    Error: {str(e)}")
    
    # Compare effects
    if 'hot' in results and 'normal' in results:
        hot_avg = results['hot']['avg_post_effect']
        normal_avg = results['normal']['avg_post_effect']
        
        results['differential_effect'] = {
            'hot_days_effect': hot_avg,
            'normal_days_effect': normal_avg,
            'difference': hot_avg - normal_avg,
            'interpretation': 'Larger effect on hot days suggests AC mechanism' 
                             if hot_avg < normal_avg else 'No differential effect by temperature'
        }
        
        print(f"\n  Differential effect (hot - normal): {hot_avg - normal_avg:.4f}")
    
    return results


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_full_event_study_analysis(level='intermediate'):
    """Run complete event study analysis."""
    
    print("\n" + "=" * 70)
    print(f"LUZ PARA TODOS EVENT STUDY - {level.upper()} LEVEL")
    print("=" * 70)
    
    # Load and prepare data
    df = load_panel(level)
    df = create_event_time(df)
    
    if 'event_year_capped' not in df.columns:
        print("ERROR: Could not create event time. Check treatment data.")
        return None
    
    all_results = {
        'level': level,
        'main_event_study': {},
        'by_intensity': {},
        'by_heat_exposure': {},
        'hot_vs_normal': {}
    }
    
    outcomes = ['deaths_all', 'deaths_cardiovascular', 'deaths_respiratory']
    
    # Main event study
    print("\n" + "=" * 50)
    print("1. MAIN EVENT STUDY")
    print("=" * 50)
    
    for outcome in outcomes:
        try:
            results = run_event_study(df, outcome)
            all_results['main_event_study'][outcome] = results
        except Exception as e:
            print(f"  Error in {outcome}: {str(e)}")
    
    # By intensity
    print("\n" + "=" * 50)
    print("2. EVENT STUDY BY TREATMENT INTENSITY")
    print("=" * 50)
    
    for outcome in ['deaths_all']:  # Focus on main outcome
        try:
            results = run_event_study_by_intensity(df, outcome)
            if results:
                all_results['by_intensity'][outcome] = results
        except Exception as e:
            print(f"  Error: {str(e)}")
    
    # By heat exposure
    print("\n" + "=" * 50)
    print("3. EVENT STUDY BY BASELINE HEAT EXPOSURE")
    print("=" * 50)
    
    for outcome in ['deaths_all']:
        try:
            results = run_event_study_by_heat_exposure(df, outcome)
            if results:
                all_results['by_heat_exposure'][outcome] = results
        except Exception as e:
            print(f"  Error: {str(e)}")
    
    # Hot vs normal days
    print("\n" + "=" * 50)
    print("4. HOT DAYS VS NORMAL DAYS")
    print("=" * 50)
    
    for outcome in ['deaths_all']:
        try:
            results = run_hot_day_event_study(df, outcome)
            if results:
                all_results['hot_vs_normal'][outcome] = results
        except Exception as e:
            print(f"  Error: {str(e)}")
    
    return all_results


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='LPT Event Study Analysis')
    parser.add_argument('--level', type=str, default='both',
                        choices=['intermediate', 'immediate', 'both'],
                        help='Spatial level to analyze')
    args = parser.parse_args()
    
    print("=" * 70)
    print("LUZ PARA TODOS EVENT STUDY ANALYSIS")
    print("=" * 70)
    
    levels = ['intermediate', 'immediate'] if args.level == 'both' else [args.level]
    
    for level in levels:
        results = run_full_event_study_analysis(level)
        
        if results:
            # Save results
            suffix = '' if level == 'intermediate' else '_immediate'
            output_file = OUTPUT_DIR / f"lpt_event_study_results{suffix}.json"
            
            with open(output_file, 'w') as f:
                json.dump(convert_to_json_serializable(results), f, indent=2)
            
            print(f"\nResults saved to: {output_file}")
            
            # Print summary
            print("\n" + "-" * 50)
            print(f"SUMMARY - {level.upper()}")
            print("-" * 50)
            
            if results['main_event_study'].get('deaths_all'):
                es = results['main_event_study']['deaths_all']
                print(f"Pre-trends: {es['pre_trend_test']['interpretation']}")
                print(f"Avg post-treatment effect: {es['avg_post_effect']:.4f}")
    
    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
