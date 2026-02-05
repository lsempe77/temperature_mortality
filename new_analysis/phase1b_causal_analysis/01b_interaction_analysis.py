#!/usr/bin/env python3
"""
ENSO × Electrification Interaction Analysis
=============================================

Combines both sources of quasi-experimental variation to estimate:
1. Whether electrification attenuates mortality during ENSO events
2. Heterogeneous treatment effects by climate conditions
3. Complementary evidence for causal effects

Identification:
- ENSO provides exogenous temperature shocks
- LPT provides exogenous variation in adaptation capacity (AC access)
- Interaction isolates causal effect of adaptation on climate vulnerability
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
# INTERACTION MODELS
# =============================================================================

def run_enso_electrification_interaction(df, outcome='deaths_all'):
    """
    Estimate interaction between ENSO and electrification.
    
    Model:
    Y_{rt} = β1*ENSO_t + β2*Electrified_{rt} + β3*(ENSO × Electrified)_{rt}
             + X_{rt}γ + α_r + δ_t + ε_{rt}
    
    β3: Differential effect of ENSO on mortality for electrified vs non-electrified regions
    """
    print(f"\n--- ENSO × Electrification Interaction: {outcome} ---")
    
    data = df.copy()
    
    # Check required variables
    required = ['oni', 'post_lpt', outcome]
    for var in required:
        if var not in data.columns:
            print(f"  {var} not found")
            return None
    
    data = data.dropna(subset=['oni', outcome])
    
    # Create variables
    data['enso_positive'] = (data['oni'] > 0).astype(float)  # El Niño tendency
    data['enso_negative'] = (data['oni'] < 0).astype(float)  # La Niña tendency
    data['is_el_nino'] = (data['oni'] >= 0.5).astype(int)
    data['is_la_nina'] = (data['oni'] <= -0.5).astype(int)
    data['electrified'] = data['post_lpt'].fillna(0)
    
    # Interactions
    data['oni_x_electrified'] = data['oni'] * data['electrified']
    data['el_nino_x_electrified'] = data['is_el_nino'] * data['electrified']
    data['la_nina_x_electrified'] = data['is_la_nina'] * data['electrified']
    
    # Controls
    controls = []
    for m in range(2, 13):
        data[f'month_{m}'] = (data['month'] == m).astype(int)
        controls.append(f'month_{m}')
    
    if 'dow' in data.columns:
        for d in range(1, 7):
            data[f'dow_{d}'] = (data['dow'] == d).astype(int)
            controls.append(f'dow_{d}')
    
    results = {}
    
    # Model 1: Continuous ONI interaction
    print("\n  Model 1: Continuous ONI × Electrification")
    
    vars_to_demean = [outcome, 'oni', 'electrified', 'oni_x_electrified'] + controls
    for col in vars_to_demean:
        if col in data.columns:
            data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    y = data[f'{outcome}_dm']
    X_cols = ['oni_dm', 'electrified_dm', 'oni_x_electrified_dm']
    X_cols += [f'{c}_dm' for c in controls]
    X = data[X_cols]
    X = sm.add_constant(X)
    
    model1 = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    results['continuous_model'] = {
        'oni_effect': {
            'coefficient': float(model1.params['oni_dm']),
            'std_error': float(model1.bse['oni_dm']),
            'p_value': float(model1.pvalues['oni_dm'])
        },
        'electrified_effect': {
            'coefficient': float(model1.params['electrified_dm']),
            'std_error': float(model1.bse['electrified_dm']),
            'p_value': float(model1.pvalues['electrified_dm'])
        },
        'interaction': {
            'coefficient': float(model1.params['oni_x_electrified_dm']),
            'std_error': float(model1.bse['oni_x_electrified_dm']),
            'p_value': float(model1.pvalues['oni_x_electrified_dm'])
        },
        'n_obs': int(model1.nobs),
        'r_squared': float(model1.rsquared)
    }
    
    print(f"    ONI effect: {results['continuous_model']['oni_effect']['coefficient']:.4f}")
    print(f"    Electrified effect: {results['continuous_model']['electrified_effect']['coefficient']:.4f}")
    print(f"    ONI × Electrified: {results['continuous_model']['interaction']['coefficient']:.4f} "
          f"(p={results['continuous_model']['interaction']['p_value']:.4f})")
    
    # Model 2: Binary ENSO phases
    print("\n  Model 2: El Niño/La Niña × Electrification")
    
    vars_to_demean = [outcome, 'is_el_nino', 'is_la_nina', 'electrified', 
                      'el_nino_x_electrified', 'la_nina_x_electrified'] + controls
    for col in vars_to_demean:
        if col in data.columns:
            data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    y = data[f'{outcome}_dm']
    X_cols = ['is_el_nino_dm', 'is_la_nina_dm', 'electrified_dm',
              'el_nino_x_electrified_dm', 'la_nina_x_electrified_dm']
    X_cols += [f'{c}_dm' for c in controls]
    X = data[X_cols]
    X = sm.add_constant(X)
    
    model2 = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    results['binary_model'] = {
        'el_nino_effect': {
            'coefficient': float(model2.params['is_el_nino_dm']),
            'std_error': float(model2.bse['is_el_nino_dm']),
            'p_value': float(model2.pvalues['is_el_nino_dm'])
        },
        'la_nina_effect': {
            'coefficient': float(model2.params['is_la_nina_dm']),
            'std_error': float(model2.bse['is_la_nina_dm']),
            'p_value': float(model2.pvalues['is_la_nina_dm'])
        },
        'el_nino_x_electrified': {
            'coefficient': float(model2.params['el_nino_x_electrified_dm']),
            'std_error': float(model2.bse['el_nino_x_electrified_dm']),
            'p_value': float(model2.pvalues['el_nino_x_electrified_dm'])
        },
        'la_nina_x_electrified': {
            'coefficient': float(model2.params['la_nina_x_electrified_dm']),
            'std_error': float(model2.bse['la_nina_x_electrified_dm']),
            'p_value': float(model2.pvalues['la_nina_x_electrified_dm'])
        },
        'n_obs': int(model2.nobs),
        'r_squared': float(model2.rsquared)
    }
    
    print(f"    El Niño × Electrified: {results['binary_model']['el_nino_x_electrified']['coefficient']:.4f} "
          f"(p={results['binary_model']['el_nino_x_electrified']['p_value']:.4f})")
    print(f"    La Niña × Electrified: {results['binary_model']['la_nina_x_electrified']['coefficient']:.4f} "
          f"(p={results['binary_model']['la_nina_x_electrified']['p_value']:.4f})")
    
    return results


def run_temperature_electrification_interaction(df, outcome='deaths_all'):
    """
    Estimate interaction between temperature extremes and electrification.
    
    Model:
    Y_{rt} = β1*ExtremeHeat_{rt} + β2*Electrified_{rt} 
             + β3*(ExtremeHeat × Electrified)_{rt}
             + X_{rt}γ + α_r + δ_t + ε_{rt}
    
    β3: Protective effect of electrification during heat extremes
    """
    print(f"\n--- Temperature × Electrification Interaction: {outcome} ---")
    
    data = df.copy()
    
    required = ['extreme_heat', 'extreme_cold', 'post_lpt', outcome]
    for var in required:
        if var not in data.columns:
            print(f"  {var} not found")
            return None
    
    data = data.dropna(subset=[outcome])
    
    # Create variables
    data['electrified'] = data['post_lpt'].fillna(0)
    
    # Interactions
    data['heat_x_electrified'] = data['extreme_heat'] * data['electrified']
    data['cold_x_electrified'] = data['extreme_cold'] * data['electrified']
    
    # Controls
    controls = []
    for m in range(2, 13):
        data[f'month_{m}'] = (data['month'] == m).astype(int)
        controls.append(f'month_{m}')
    
    if 'dow' in data.columns:
        for d in range(1, 7):
            data[f'dow_{d}'] = (data['dow'] == d).astype(int)
            controls.append(f'dow_{d}')
    
    # Demean
    vars_to_demean = [outcome, 'extreme_heat', 'extreme_cold', 'electrified',
                      'heat_x_electrified', 'cold_x_electrified'] + controls
    for col in vars_to_demean:
        if col in data.columns:
            data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    # Run regression
    y = data[f'{outcome}_dm']
    X_cols = ['extreme_heat_dm', 'extreme_cold_dm', 'electrified_dm',
              'heat_x_electrified_dm', 'cold_x_electrified_dm']
    X_cols += [f'{c}_dm' for c in controls]
    X = data[X_cols]
    X = sm.add_constant(X)
    
    model = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    results = {
        'extreme_heat_effect': {
            'coefficient': float(model.params['extreme_heat_dm']),
            'std_error': float(model.bse['extreme_heat_dm']),
            'p_value': float(model.pvalues['extreme_heat_dm'])
        },
        'extreme_cold_effect': {
            'coefficient': float(model.params['extreme_cold_dm']),
            'std_error': float(model.bse['extreme_cold_dm']),
            'p_value': float(model.pvalues['extreme_cold_dm'])
        },
        'heat_x_electrified': {
            'coefficient': float(model.params['heat_x_electrified_dm']),
            'std_error': float(model.bse['heat_x_electrified_dm']),
            'p_value': float(model.pvalues['heat_x_electrified_dm']),
            'ci_lower': float(model.conf_int().loc['heat_x_electrified_dm', 0]),
            'ci_upper': float(model.conf_int().loc['heat_x_electrified_dm', 1])
        },
        'cold_x_electrified': {
            'coefficient': float(model.params['cold_x_electrified_dm']),
            'std_error': float(model.bse['cold_x_electrified_dm']),
            'p_value': float(model.pvalues['cold_x_electrified_dm']),
            'ci_lower': float(model.conf_int().loc['cold_x_electrified_dm', 0]),
            'ci_upper': float(model.conf_int().loc['cold_x_electrified_dm', 1])
        },
        'n_obs': int(model.nobs),
        'r_squared': float(model.rsquared),
    }
    
    print(f"  Extreme heat effect: {results['extreme_heat_effect']['coefficient']:.4f}")
    print(f"  Heat × Electrified: {results['heat_x_electrified']['coefficient']:.4f} "
          f"(p={results['heat_x_electrified']['p_value']:.4f})")
    print(f"  Extreme cold effect: {results['extreme_cold_effect']['coefficient']:.4f}")
    print(f"  Cold × Electrified: {results['cold_x_electrified']['coefficient']:.4f} "
          f"(p={results['cold_x_electrified']['p_value']:.4f})")
    
    # Interpretation
    if results['heat_x_electrified']['coefficient'] < 0:
        results['interpretation'] = {
            'heat_protection': 'Electrification reduces heat-related mortality',
            'mechanism': 'Consistent with AC access hypothesis',
            'magnitude': f"{abs(results['heat_x_electrified']['coefficient']):.2f} fewer deaths per extreme heat day"
        }
    else:
        results['interpretation'] = {
            'heat_protection': 'No protective effect of electrification on heat days',
            'mechanism': 'AC mechanism not supported'
        }
    
    return results


def run_triple_interaction(df, outcome='deaths_all'):
    """
    Triple interaction: ENSO × Temperature × Electrification
    
    Model with three-way interaction to test:
    - Does electrification protect more during ENSO-induced temperature extremes?
    """
    print(f"\n--- Triple Interaction: ENSO × Heat × Electrification ---")
    
    data = df.copy()
    
    required = ['oni', 'extreme_heat', 'post_lpt', outcome]
    for var in required:
        if var not in data.columns:
            print(f"  {var} not found")
            return None
    
    data = data.dropna(subset=['oni', outcome])
    
    # Create variables
    data['is_el_nino'] = (data['oni'] >= 0.5).astype(int)
    data['heat'] = data['extreme_heat'].fillna(0)
    data['electrified'] = data['post_lpt'].fillna(0)
    
    # All interactions
    data['en_heat'] = data['is_el_nino'] * data['heat']
    data['en_elec'] = data['is_el_nino'] * data['electrified']
    data['heat_elec'] = data['heat'] * data['electrified']
    data['triple'] = data['is_el_nino'] * data['heat'] * data['electrified']
    
    # Controls
    controls = []
    for m in range(2, 13):
        data[f'month_{m}'] = (data['month'] == m).astype(int)
        controls.append(f'month_{m}')
    
    # Demean
    vars_to_demean = [outcome, 'is_el_nino', 'heat', 'electrified',
                      'en_heat', 'en_elec', 'heat_elec', 'triple'] + controls
    for col in vars_to_demean:
        if col in data.columns:
            data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    # Regression
    y = data[f'{outcome}_dm']
    X_cols = ['is_el_nino_dm', 'heat_dm', 'electrified_dm',
              'en_heat_dm', 'en_elec_dm', 'heat_elec_dm', 'triple_dm']
    X_cols += [f'{c}_dm' for c in controls]
    X = data[X_cols]
    X = sm.add_constant(X)
    
    model = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    results = {
        'main_effects': {
            'el_nino': float(model.params.get('is_el_nino_dm', np.nan)),
            'heat': float(model.params.get('heat_dm', np.nan)),
            'electrified': float(model.params.get('electrified_dm', np.nan)),
        },
        'two_way_interactions': {
            'el_nino_x_heat': float(model.params.get('en_heat_dm', np.nan)),
            'el_nino_x_electrified': float(model.params.get('en_elec_dm', np.nan)),
            'heat_x_electrified': float(model.params.get('heat_elec_dm', np.nan)),
        },
        'triple_interaction': {
            'coefficient': float(model.params['triple_dm']),
            'std_error': float(model.bse['triple_dm']),
            'p_value': float(model.pvalues['triple_dm']),
            'ci_lower': float(model.conf_int().loc['triple_dm', 0]),
            'ci_upper': float(model.conf_int().loc['triple_dm', 1])
        },
        'n_obs': int(model.nobs),
        'r_squared': float(model.rsquared)
    }
    
    print(f"  Triple interaction (EN × Heat × Elec): {results['triple_interaction']['coefficient']:.4f}")
    print(f"  p-value: {results['triple_interaction']['p_value']:.4f}")
    
    # Interpretation
    if results['triple_interaction']['coefficient'] < 0 and results['triple_interaction']['p_value'] < 0.1:
        results['interpretation'] = (
            "Electrification provides additional protection during El Niño-induced heat. "
            "This supports the causal mechanism: ENSO → temperature → mortality, "
            "with electrification (AC access) as a moderator."
        )
    else:
        results['interpretation'] = (
            "No significant triple interaction. Electrification effect is not "
            "differentially larger during ENSO-induced heat."
        )
    
    return results


def run_heterogeneity_analysis(df, outcome='deaths_all'):
    """
    Analyze heterogeneity in electrification effects by region characteristics.
    """
    print(f"\n--- Heterogeneity Analysis ---")
    
    data = df.copy()
    
    if 'post_lpt' not in data.columns:
        return None
    
    results = {}
    
    # By baseline temperature (hot vs cool regions)
    if 'temp_mean' in data.columns:
        print("\n  By baseline temperature:")
        
        baseline_temp = data.groupby('region_code')['temp_mean'].mean()
        median_temp = baseline_temp.median()
        
        data['hot_region'] = data['region_code'].map(
            {r: 1 if t >= median_temp else 0 for r, t in baseline_temp.items()}
        )
        
        hot_regions = data[data['hot_region'] == 1]
        cool_regions = data[data['hot_region'] == 0]
        
        for subset, name in [(hot_regions, 'hot_regions'), (cool_regions, 'cool_regions')]:
            if len(subset) > 0 and 'post_lpt' in subset.columns:
                subset = subset.copy()
                for col in [outcome, 'post_lpt']:
                    if col in subset.columns:
                        subset[f'{col}_dm'] = subset.groupby('region_code')[col].transform(
                            lambda x: x - x.mean()
                        )
                
                if f'{outcome}_dm' in subset.columns and 'post_lpt_dm' in subset.columns:
                    y = subset[f'{outcome}_dm'].dropna()
                    X = sm.add_constant(subset.loc[y.index, 'post_lpt_dm'])
                    
                    try:
                        model = sm.OLS(y, X).fit(cov_type='cluster', 
                                                  cov_kwds={'groups': subset.loc[y.index, 'region_code']})
                        results[name] = {
                            'att': float(model.params['post_lpt_dm']),
                            'std_error': float(model.bse['post_lpt_dm']),
                            'p_value': float(model.pvalues['post_lpt_dm']),
                            'n_obs': int(model.nobs)
                        }
                        print(f"    {name}: ATT = {results[name]['att']:.4f} "
                              f"(p = {results[name]['p_value']:.4f})")
                    except:
                        pass
    
    # By urbanization
    if 'urban_rate' in data.columns:
        print("\n  By urbanization:")
        
        baseline_urban = data.groupby('region_code')['urban_rate'].first()
        median_urban = baseline_urban.median()
        
        data['urban_region'] = data['region_code'].map(
            {r: 1 if u >= median_urban else 0 for r, u in baseline_urban.items()}
        )
        
        for urban_val, name in [(1, 'urban_regions'), (0, 'rural_regions')]:
            subset = data[data['urban_region'] == urban_val].copy()
            
            if len(subset) > 0 and 'post_lpt' in subset.columns:
                for col in [outcome, 'post_lpt']:
                    if col in subset.columns:
                        subset[f'{col}_dm'] = subset.groupby('region_code')[col].transform(
                            lambda x: x - x.mean()
                        )
                
                if f'{outcome}_dm' in subset.columns and 'post_lpt_dm' in subset.columns:
                    y = subset[f'{outcome}_dm'].dropna()
                    X = sm.add_constant(subset.loc[y.index, 'post_lpt_dm'])
                    
                    try:
                        model = sm.OLS(y, X).fit(cov_type='cluster',
                                                  cov_kwds={'groups': subset.loc[y.index, 'region_code']})
                        results[name] = {
                            'att': float(model.params['post_lpt_dm']),
                            'std_error': float(model.bse['post_lpt_dm']),
                            'p_value': float(model.pvalues['post_lpt_dm']),
                            'n_obs': int(model.nobs)
                        }
                        print(f"    {name}: ATT = {results[name]['att']:.4f} "
                              f"(p = {results[name]['p_value']:.4f})")
                    except:
                        pass
    
    return results


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_full_interaction_analysis(level='intermediate'):
    """Run complete interaction analysis."""
    
    print("\n" + "=" * 70)
    print(f"INTERACTION ANALYSIS - {level.upper()} LEVEL")
    print("=" * 70)
    
    # Load data
    df = load_panel(level)
    
    all_results = {
        'level': level,
        'enso_electrification': {},
        'temperature_electrification': {},
        'triple_interaction': {},
        'heterogeneity': {}
    }
    
    outcomes = ['deaths_all', 'deaths_cardiovascular', 'deaths_respiratory']
    
    # ENSO × Electrification
    print("\n" + "=" * 50)
    print("1. ENSO × ELECTRIFICATION")
    print("=" * 50)
    
    for outcome in outcomes:
        try:
            results = run_enso_electrification_interaction(df, outcome)
            if results:
                all_results['enso_electrification'][outcome] = results
        except Exception as e:
            print(f"  Error in {outcome}: {str(e)}")
    
    # Temperature × Electrification
    print("\n" + "=" * 50)
    print("2. TEMPERATURE × ELECTRIFICATION")
    print("=" * 50)
    
    for outcome in outcomes:
        try:
            results = run_temperature_electrification_interaction(df, outcome)
            if results:
                all_results['temperature_electrification'][outcome] = results
        except Exception as e:
            print(f"  Error in {outcome}: {str(e)}")
    
    # Triple interaction
    print("\n" + "=" * 50)
    print("3. TRIPLE INTERACTION")
    print("=" * 50)
    
    for outcome in ['deaths_all']:  # Focus on main outcome
        try:
            results = run_triple_interaction(df, outcome)
            if results:
                all_results['triple_interaction'][outcome] = results
        except Exception as e:
            print(f"  Error: {str(e)}")
    
    # Heterogeneity
    print("\n" + "=" * 50)
    print("4. HETEROGENEITY ANALYSIS")
    print("=" * 50)
    
    for outcome in ['deaths_all']:
        try:
            results = run_heterogeneity_analysis(df, outcome)
            if results:
                all_results['heterogeneity'][outcome] = results
        except Exception as e:
            print(f"  Error: {str(e)}")
    
    return all_results


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='ENSO × Electrification Interaction Analysis')
    parser.add_argument('--level', type=str, default='both',
                        choices=['intermediate', 'immediate', 'both'],
                        help='Spatial level to analyze')
    args = parser.parse_args()
    
    print("=" * 70)
    print("ENSO × ELECTRIFICATION INTERACTION ANALYSIS")
    print("=" * 70)
    
    levels = ['intermediate', 'immediate'] if args.level == 'both' else [args.level]
    
    for level in levels:
        results = run_full_interaction_analysis(level)
        
        if results:
            # Save results
            suffix = '' if level == 'intermediate' else '_immediate'
            output_file = OUTPUT_DIR / f"interaction_results{suffix}.json"
            
            with open(output_file, 'w') as f:
                json.dump(convert_to_json_serializable(results), f, indent=2)
            
            print(f"\nResults saved to: {output_file}")
            
            # Summary
            print("\n" + "-" * 50)
            print(f"KEY FINDINGS - {level.upper()}")
            print("-" * 50)
            
            if results['temperature_electrification'].get('deaths_all'):
                te = results['temperature_electrification']['deaths_all']
                print(f"Heat × Electrified: {te['heat_x_electrified']['coefficient']:.4f} "
                      f"(p={te['heat_x_electrified']['p_value']:.4f})")
                if 'interpretation' in te:
                    print(f"  → {te['interpretation'].get('heat_protection', '')}")
    
    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
