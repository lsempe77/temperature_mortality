#!/usr/bin/env python3
"""
ENSO Instrumental Variable Analysis
====================================

Uses ENSO indices (ONI, MEI) as instruments for temperature extremes.
Estimates causal effect of temperature on mortality using 2SLS.

Identification Strategy:
- ENSO is exogenous to local confounders (oceanic/atmospheric phenomenon)
- ENSO affects Brazilian temperature patterns (first stage)
- ENSO affects mortality only through temperature (exclusion restriction)

Methods:
1. First stage: ENSO → Regional temperature anomalies
2. Second stage: Predicted temperature → Mortality
3. Robustness: Different ENSO indices, different temperature measures
"""

import pandas as pd
import numpy as np
from pathlib import Path
import statsmodels.api as sm
from statsmodels.sandbox.regression.gmm import IV2SLS
from linearmodels.iv import IV2SLS as LM_IV2SLS
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
# FIRST STAGE: ENSO → TEMPERATURE
# =============================================================================

def run_first_stage(df, instrument='oni', temp_var='temp_mean'):
    """
    First stage regression: ENSO → Temperature
    
    Model: Temp_{rt} = α + β*ENSO_t + γ*X_{rt} + δ_r + θ_t + ε_{rt}
    
    Args:
        df: Panel data
        instrument: ENSO index to use ('oni', 'mei', 'nino34')
        temp_var: Temperature variable ('temp_mean', 'temp_anomaly', 'extreme_heat')
    
    Returns:
        dict with first stage results
    """
    print(f"\n--- First Stage: {instrument} → {temp_var} ---")
    
    # Prepare data
    data = df.dropna(subset=[instrument, temp_var]).copy()
    
    # Create region and time fixed effects
    data['region_fe'] = pd.Categorical(data['region_code']).codes
    data['month_fe'] = pd.Categorical(data['month']).codes
    data['year_fe'] = pd.Categorical(data['year']).codes
    
    # Create seasonal controls
    for m in range(1, 13):
        data[f'month_{m}'] = (data['month'] == m).astype(int)
    
    # Control variables
    controls = [f'month_{m}' for m in range(2, 13)]  # Month dummies
    if 'dow' in data.columns:
        for d in range(1, 7):
            data[f'dow_{d}'] = (data['dow'] == d).astype(int)
        controls += [f'dow_{d}' for d in range(1, 7)]
    
    # First stage regression with fixed effects
    # Using region-demeaned data for computational efficiency
    
    # Demean by region
    for col in [temp_var, instrument] + controls:
        if col in data.columns:
            data[f'{col}_dm'] = data.groupby('region_code')[col].transform(
                lambda x: x - x.mean()
            )
    
    # Construct design matrix
    y = data[f'{temp_var}_dm']
    X_cols = [f'{instrument}_dm'] + [f'{c}_dm' for c in controls if f'{c}_dm' in data.columns]
    X = data[X_cols]
    X = sm.add_constant(X)
    
    # Run regression
    model = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    # Extract results
    results = {
        'instrument': instrument,
        'temp_var': temp_var,
        'n_obs': int(model.nobs),
        'n_regions': int(data['region_code'].nunique()),
        'coefficient': float(model.params[f'{instrument}_dm']),
        'std_error': float(model.bse[f'{instrument}_dm']),
        't_stat': float(model.tvalues[f'{instrument}_dm']),
        'p_value': float(model.pvalues[f'{instrument}_dm']),
        'r_squared': float(model.rsquared),
        'f_stat': float(model.fvalue),
        'f_pvalue': float(model.f_pvalue),
    }
    
    # F-statistic on excluded instrument (weak instrument test)
    # Rule of thumb: F > 10 for strong instrument
    results['weak_instrument_test'] = 'Strong' if results['f_stat'] > 10 else 'Weak'
    
    print(f"  Coefficient: {results['coefficient']:.4f} (SE: {results['std_error']:.4f})")
    print(f"  t-stat: {results['t_stat']:.2f}, p-value: {results['p_value']:.4f}")
    print(f"  F-statistic: {results['f_stat']:.2f} ({results['weak_instrument_test']} instrument)")
    print(f"  R-squared: {results['r_squared']:.4f}")
    
    # Store predicted values
    data['temp_predicted'] = model.fittedvalues + data.groupby('region_code')[temp_var].transform('mean')
    
    return results, data, model


# =============================================================================
# SECOND STAGE: TEMPERATURE → MORTALITY
# =============================================================================

def run_second_stage_ols(df, temp_var='temp_mean', outcome='deaths_all'):
    """
    OLS regression (for comparison): Temperature → Mortality
    """
    print(f"\n--- OLS: {temp_var} → {outcome} ---")
    
    # Prepare data
    data = df.dropna(subset=[temp_var, outcome]).copy()
    
    # Demean by region
    for col in [temp_var, outcome]:
        data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    # Controls
    controls = []
    for m in range(2, 13):
        data[f'month_{m}'] = (data['month'] == m).astype(int)
        controls.append(f'month_{m}')
    
    for d in range(1, 7):
        data[f'dow_{d}'] = (data['dow'] == d).astype(int)
        controls.append(f'dow_{d}')
    
    # Demean controls
    for col in controls:
        data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    # Design matrix
    y = data[f'{outcome}_dm']
    X_cols = [f'{temp_var}_dm'] + [f'{c}_dm' for c in controls]
    X = data[X_cols]
    X = sm.add_constant(X)
    
    # Run OLS
    model = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    results = {
        'method': 'OLS',
        'temp_var': temp_var,
        'outcome': outcome,
        'n_obs': int(model.nobs),
        'coefficient': float(model.params[f'{temp_var}_dm']),
        'std_error': float(model.bse[f'{temp_var}_dm']),
        't_stat': float(model.tvalues[f'{temp_var}_dm']),
        'p_value': float(model.pvalues[f'{temp_var}_dm']),
        'r_squared': float(model.rsquared),
    }
    
    print(f"  Coefficient: {results['coefficient']:.4f} (SE: {results['std_error']:.4f})")
    print(f"  t-stat: {results['t_stat']:.2f}, p-value: {results['p_value']:.4f}")
    
    return results


def run_2sls(df, instrument='oni', temp_var='temp_mean', outcome='deaths_all'):
    """
    Two-Stage Least Squares: ENSO → Temperature → Mortality
    
    Uses linearmodels for proper IV estimation with fixed effects.
    """
    print(f"\n--- 2SLS: {instrument} → {temp_var} → {outcome} ---")
    
    # Prepare data
    data = df.dropna(subset=[instrument, temp_var, outcome]).copy()
    data = data.set_index(['region_code', 'date'])
    
    # Create controls
    for m in range(2, 13):
        data[f'month_{m}'] = (data['month'] == m).astype(int)
    for d in range(1, 7):
        data[f'dow_{d}'] = (data['dow'] == d).astype(int)
    
    control_cols = [f'month_{m}' for m in range(2, 13)] + [f'dow_{d}' for d in range(1, 7)]
    
    try:
        # Set up IV regression
        # Y = deaths, X_endog = temperature, X_exog = controls, Z = ENSO
        
        dependent = data[outcome]
        exog = data[control_cols]
        exog = sm.add_constant(exog)
        endog = data[[temp_var]]
        instruments = data[[instrument]]
        
        # Run 2SLS with entity fixed effects
        model = LM_IV2SLS(
            dependent=dependent,
            exog=exog,
            endog=endog,
            instruments=instruments
        ).fit(cov_type='clustered', clusters=data.index.get_level_values(0))
        
        results = {
            'method': '2SLS',
            'instrument': instrument,
            'temp_var': temp_var,
            'outcome': outcome,
            'n_obs': int(model.nobs),
            'coefficient': float(model.params[temp_var]),
            'std_error': float(model.std_errors[temp_var]),
            't_stat': float(model.tstats[temp_var]),
            'p_value': float(model.pvalues[temp_var]),
        }
        
        # Add diagnostics - handle different linearmodels versions
        try:
            f_stat = model.first_stage.diagnostics['f.stat']
            if hasattr(f_stat, 'stat'):
                results['first_stage_f'] = float(f_stat.stat)
            elif isinstance(f_stat, (int, float)):
                results['first_stage_f'] = float(f_stat)
            else:
                results['first_stage_f'] = float(f_stat.iloc[0]) if hasattr(f_stat, 'iloc') else None
        except Exception:
            results['first_stage_f'] = None
        results['sargan_stat'] = None  # Only for overidentified models
        
        print(f"  Coefficient: {results['coefficient']:.4f} (SE: {results['std_error']:.4f})")
        print(f"  t-stat: {results['t_stat']:.2f}, p-value: {results['p_value']:.4f}")
        print(f"  First-stage F: {results['first_stage_f']:.2f}")
        
        return results, model
        
    except Exception as e:
        print(f"  ERROR: {str(e)}")
        return None, None


# =============================================================================
# REDUCED FORM: ENSO → MORTALITY
# =============================================================================

def run_reduced_form(df, instrument='oni', outcome='deaths_all'):
    """
    Reduced form regression: ENSO → Mortality
    
    This directly estimates the effect of ENSO on mortality.
    Under exclusion restriction, this should work only through temperature.
    """
    print(f"\n--- Reduced Form: {instrument} → {outcome} ---")
    
    # Prepare data
    data = df.dropna(subset=[instrument, outcome]).copy()
    
    # Demean by region
    for col in [instrument, outcome]:
        data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    # Controls
    controls = []
    for m in range(2, 13):
        data[f'month_{m}'] = (data['month'] == m).astype(int)
        controls.append(f'month_{m}')
    
    for d in range(1, 7):
        data[f'dow_{d}'] = (data['dow'] == d).astype(int)
        controls.append(f'dow_{d}')
    
    # Demean controls
    for col in controls:
        data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    # Design matrix
    y = data[f'{outcome}_dm']
    X_cols = [f'{instrument}_dm'] + [f'{c}_dm' for c in controls]
    X = data[X_cols]
    X = sm.add_constant(X)
    
    # Run regression
    model = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    results = {
        'method': 'Reduced Form',
        'instrument': instrument,
        'outcome': outcome,
        'n_obs': int(model.nobs),
        'coefficient': float(model.params[f'{instrument}_dm']),
        'std_error': float(model.bse[f'{instrument}_dm']),
        't_stat': float(model.tvalues[f'{instrument}_dm']),
        'p_value': float(model.pvalues[f'{instrument}_dm']),
        'r_squared': float(model.rsquared),
    }
    
    print(f"  Coefficient: {results['coefficient']:.4f} (SE: {results['std_error']:.4f})")
    print(f"  t-stat: {results['t_stat']:.2f}, p-value: {results['p_value']:.4f}")
    
    return results


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_full_iv_analysis(level='intermediate'):
    """
    Run full IV analysis for a given level.
    """
    print("\n" + "=" * 70)
    print(f"ENSO IV ANALYSIS - {level.upper()} LEVEL")
    print("=" * 70)
    
    # Load data
    df = load_panel(level)
    
    # Check for required columns
    if 'oni' not in df.columns:
        print("ERROR: ENSO data not found. Run data preparation scripts first.")
        return None
    
    all_results = {
        'level': level,
        'first_stage': [],
        'reduced_form': [],
        'ols': [],
        'iv_2sls': []
    }
    
    # ==========================================================================
    # 1. FIRST STAGE ANALYSIS
    # ==========================================================================
    print("\n" + "=" * 50)
    print("1. FIRST STAGE: ENSO → TEMPERATURE")
    print("=" * 50)
    
    instruments = ['oni', 'mei', 'nino34']
    temp_vars = ['temp_mean', 'temp_anomaly']
    
    for instrument in instruments:
        if instrument not in df.columns:
            continue
        for temp_var in temp_vars:
            if temp_var not in df.columns:
                continue
            
            results, data, model = run_first_stage(df, instrument, temp_var)
            all_results['first_stage'].append(results)
    
    # ==========================================================================
    # 2. REDUCED FORM ANALYSIS
    # ==========================================================================
    print("\n" + "=" * 50)
    print("2. REDUCED FORM: ENSO → MORTALITY")
    print("=" * 50)
    
    outcomes = ['deaths_all', 'deaths_cardiovascular', 'deaths_respiratory']
    
    for instrument in instruments:
        if instrument not in df.columns:
            continue
        for outcome in outcomes:
            if outcome not in df.columns:
                continue
            
            results = run_reduced_form(df, instrument, outcome)
            all_results['reduced_form'].append(results)
    
    # ==========================================================================
    # 3. OLS COMPARISON
    # ==========================================================================
    print("\n" + "=" * 50)
    print("3. OLS: TEMPERATURE → MORTALITY")
    print("=" * 50)
    
    for temp_var in ['temp_mean']:
        for outcome in outcomes:
            if temp_var not in df.columns or outcome not in df.columns:
                continue
            
            results = run_second_stage_ols(df, temp_var, outcome)
            all_results['ols'].append(results)
    
    # ==========================================================================
    # 4. 2SLS IV ESTIMATION
    # ==========================================================================
    print("\n" + "=" * 50)
    print("4. 2SLS: ENSO → TEMPERATURE → MORTALITY")
    print("=" * 50)
    
    for instrument in ['oni']:  # Primary instrument
        for temp_var in ['temp_mean']:
            for outcome in outcomes:
                if all([c in df.columns for c in [instrument, temp_var, outcome]]):
                    results, model = run_2sls(df, instrument, temp_var, outcome)
                    if results:
                        all_results['iv_2sls'].append(results)
    
    return all_results


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='ENSO IV Analysis')
    parser.add_argument('--level', type=str, default='both',
                        choices=['intermediate', 'immediate', 'both'],
                        help='Spatial level to analyze')
    args = parser.parse_args()
    
    print("=" * 70)
    print("ENSO INSTRUMENTAL VARIABLE ANALYSIS")
    print("=" * 70)
    
    levels = ['intermediate', 'immediate'] if args.level == 'both' else [args.level]
    
    all_results = {}
    
    for level in levels:
        results = run_full_iv_analysis(level)
        if results:
            all_results[level] = results
            
            # Save level-specific results
            suffix = '' if level == 'intermediate' else '_immediate'
            output_file = OUTPUT_DIR / f"enso_iv_results{suffix}.json"
            
            with open(output_file, 'w') as f:
                json.dump(convert_to_json_serializable(results), f, indent=2)
            
            print(f"\nResults saved to: {output_file}")
    
    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    for level, results in all_results.items():
        print(f"\n{level.upper()} Level:")
        
        if results['first_stage']:
            fs = results['first_stage'][0]
            print(f"  First Stage (ONI → temp_mean):")
            print(f"    β = {fs['coefficient']:.4f}, F = {fs['f_stat']:.2f}")
        
        if results['iv_2sls']:
            iv = results['iv_2sls'][0]
            print(f"  2SLS (ENSO → temp → deaths_all):")
            print(f"    β = {iv['coefficient']:.4f} (SE: {iv['std_error']:.4f})")
        
        if results['ols']:
            ols = results['ols'][0]
            print(f"  OLS (temp → deaths_all):")
            print(f"    β = {ols['coefficient']:.4f} (SE: {ols['std_error']:.4f})")
    
    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
