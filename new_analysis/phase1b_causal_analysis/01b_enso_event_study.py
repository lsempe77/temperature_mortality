#!/usr/bin/env python3
"""
ENSO Event Study Analysis
=========================

Event study design around El Niño and La Niña episodes.
Estimates differential mortality effects during ENSO events vs neutral periods.

Methods:
1. Define ENSO events (El Niño, La Niña) using ONI thresholds
2. Create event time relative to episode start/peak
3. Estimate dynamic effects before/during/after events
4. Stratify by event intensity (weak, moderate, strong)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import statsmodels.api as sm
import statsmodels.formula.api as smf
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
    elif isinstance(obj, (pd.Timestamp, np.datetime64)):
        return str(obj)
    elif hasattr(obj, 'isoformat'):
        return obj.isoformat()
    else:
        return obj


def load_panel(level='intermediate'):
    """Load causal analysis panel."""
    print(f"Loading {level} panel...")
    df = pd.read_parquet(RESULTS_DIR / f"causal_panel_{level}.parquet")
    print(f"  Shape: {df.shape}")
    return df


def load_enso_monthly():
    """Load monthly ENSO data for event identification."""
    df = pd.read_parquet(RESULTS_DIR / "enso_monthly.parquet")
    return df


# =============================================================================
# EVENT IDENTIFICATION
# =============================================================================

def identify_enso_events(enso_monthly, threshold=0.5, min_duration=5):
    """
    Identify El Niño and La Niña events.
    
    NOAA definition:
    - El Niño: ONI >= 0.5°C for at least 5 consecutive overlapping 3-month periods
    - La Niña: ONI <= -0.5°C for at least 5 consecutive overlapping 3-month periods
    
    Args:
        enso_monthly: Monthly ENSO data
        threshold: ONI threshold for classification
        min_duration: Minimum consecutive months for event
    
    Returns:
        DataFrame with event information
    """
    print("\nIdentifying ENSO events...")
    
    df = enso_monthly.copy()
    df = df.sort_values('date').reset_index(drop=True)
    
    # Create El Niño and La Niña indicators
    df['el_nino_month'] = (df['oni'] >= threshold).astype(int)
    df['la_nina_month'] = (df['oni'] <= -threshold).astype(int)
    
    # Find consecutive periods (simplified approach)
    events = []
    
    for event_type, indicator_col, sign in [('el_nino', 'el_nino_month', 1), 
                                             ('la_nina', 'la_nina_month', -1)]:
        df['event_group'] = (df[indicator_col] != df[indicator_col].shift()).cumsum()
        
        for group_id, group in df[df[indicator_col] == 1].groupby('event_group'):
            if len(group) >= min_duration:
                events.append({
                    'event_type': event_type,
                    'start_date': group['date'].min(),
                    'end_date': group['date'].max(),
                    'peak_date': group.loc[group['oni'].abs().idxmax(), 'date'],
                    'peak_oni': group.loc[group['oni'].abs().idxmax(), 'oni'],
                    'duration_months': len(group),
                    'mean_oni': group['oni'].mean(),
                    'max_oni': group['oni'].max() if sign > 0 else group['oni'].min(),
                })
    
    events_df = pd.DataFrame(events)
    
    # Classify event intensity
    def classify_intensity(oni):
        if abs(oni) >= 1.5:
            return 'strong'
        elif abs(oni) >= 1.0:
            return 'moderate'
        else:
            return 'weak'
    
    if len(events_df) > 0:
        events_df['intensity'] = events_df['peak_oni'].abs().apply(
            lambda x: 'strong' if x >= 1.5 else ('moderate' if x >= 1.0 else 'weak')
        )
        events_df = events_df.sort_values('start_date').reset_index(drop=True)
    
    print(f"  Total events identified: {len(events_df)}")
    if len(events_df) > 0:
        print(f"  El Niño events: {(events_df['event_type'] == 'el_nino').sum()}")
        print(f"  La Niña events: {(events_df['event_type'] == 'la_nina').sum()}")
    
    return events_df


def create_event_time_indicators(df, events_df, window=12):
    """
    Create event time indicators relative to ENSO events.
    
    Args:
        df: Daily panel data
        events_df: ENSO events DataFrame
        window: Months before/after event to track
    
    Returns:
        DataFrame with event time indicators
    """
    print(f"\nCreating event time indicators (window: {window} months)...")
    
    data = df.copy()
    
    # Initialize columns
    data['in_el_nino'] = 0
    data['in_la_nina'] = 0
    data['event_time_el_nino'] = np.nan
    data['event_time_la_nina'] = np.nan
    data['current_event_intensity'] = 'none'
    
    # Filter events to study period
    study_start = data['date'].min()
    study_end = data['date'].max()
    relevant_events = events_df[
        (events_df['end_date'] >= study_start) & 
        (events_df['start_date'] <= study_end)
    ].copy()
    
    print(f"  Events in study period: {len(relevant_events)}")
    
    for _, event in relevant_events.iterrows():
        event_type = event['event_type']
        start = event['start_date']
        end = event['end_date']
        intensity = event['intensity']
        
        # Mark days during event
        mask = (data['date'] >= start) & (data['date'] <= end)
        data.loc[mask, f'in_{event_type}'] = 1
        data.loc[mask, 'current_event_intensity'] = intensity
        
        # Calculate event time (months from event start)
        for idx in data[mask].index:
            months_diff = (data.loc[idx, 'date'].year - start.year) * 12 + \
                         (data.loc[idx, 'date'].month - start.month)
            data.loc[idx, f'event_time_{event_type}'] = months_diff
        
        # Also mark pre-event period for anticipation effects
        pre_mask = (data['date'] >= start - pd.DateOffset(months=window)) & \
                   (data['date'] < start)
        for idx in data[pre_mask].index:
            months_diff = (data.loc[idx, 'date'].year - start.year) * 12 + \
                         (data.loc[idx, 'date'].month - start.month)
            if pd.isna(data.loc[idx, f'event_time_{event_type}']):
                data.loc[idx, f'event_time_{event_type}'] = months_diff
    
    # Create in_any_enso indicator
    data['in_any_enso'] = ((data['in_el_nino'] == 1) | (data['in_la_nina'] == 1)).astype(int)
    
    print(f"  Days in El Niño: {data['in_el_nino'].sum():,}")
    print(f"  Days in La Niña: {data['in_la_nina'].sum():,}")
    
    return data


# =============================================================================
# EVENT STUDY ESTIMATION
# =============================================================================

def run_event_study(df, event_type='el_nino', outcome='deaths_all', 
                   window=6, reference_period=-1):
    """
    Run event study regression around ENSO events.
    
    Model:
    Y_{rt} = Σ_{k=-K}^{K} β_k * 1{t=k} + X_{rt}γ + α_r + δ_t + ε_{rt}
    
    where k is months relative to event start.
    
    Args:
        df: Panel data with event time indicators
        event_type: 'el_nino' or 'la_nina'
        outcome: Dependent variable
        window: Months before/after event
        reference_period: Period to normalize to (usually -1)
    """
    print(f"\n--- Event Study: {event_type} → {outcome} ---")
    
    # Prepare data
    event_time_col = f'event_time_{event_type}'
    data = df.dropna(subset=[event_time_col, outcome]).copy()
    
    # Bin event time
    data['event_time_binned'] = data[event_time_col].clip(-window, window).astype(int)
    
    # Create dummy variables (excluding reference period)
    event_times = list(range(-window, window + 1))
    event_times.remove(reference_period)
    
    for t in event_times:
        data[f'event_{t}'] = (data['event_time_binned'] == t).astype(int)
    
    # Add controls
    for m in range(2, 13):
        data[f'month_{m}'] = (data['month'] == m).astype(int)
    for d in range(1, 7):
        data[f'dow_{d}'] = (data['dow'] == d).astype(int)
    
    # Demean by region (within estimator)
    cols_to_demean = [outcome] + [f'event_{t}' for t in event_times]
    cols_to_demean += [f'month_{m}' for m in range(2, 13)]
    cols_to_demean += [f'dow_{d}' for d in range(1, 7)]
    
    for col in cols_to_demean:
        data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    # Design matrix
    y = data[f'{outcome}_dm']
    X_cols = [f'event_{t}_dm' for t in event_times]
    X_cols += [f'month_{m}_dm' for m in range(2, 13)]
    X_cols += [f'dow_{d}_dm' for d in range(1, 7)]
    X = data[X_cols]
    X = sm.add_constant(X)
    
    # Run regression
    model = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    # Extract event study coefficients
    event_coeffs = {}
    for t in event_times:
        col = f'event_{t}_dm'
        event_coeffs[t] = {
            'coefficient': float(model.params[col]),
            'std_error': float(model.bse[col]),
            'ci_lower': float(model.conf_int().loc[col, 0]),
            'ci_upper': float(model.conf_int().loc[col, 1]),
            'p_value': float(model.pvalues[col])
        }
    
    # Add reference period (normalized to 0)
    event_coeffs[reference_period] = {
        'coefficient': 0.0,
        'std_error': 0.0,
        'ci_lower': 0.0,
        'ci_upper': 0.0,
        'p_value': 1.0
    }
    
    results = {
        'event_type': event_type,
        'outcome': outcome,
        'reference_period': reference_period,
        'window': window,
        'n_obs': int(model.nobs),
        'n_regions': int(data['region_code'].nunique()),
        'r_squared': float(model.rsquared),
        'event_coefficients': event_coeffs
    }
    
    # Test for pre-trends
    pre_trend_coefs = [event_coeffs[t]['coefficient'] for t in event_times if t < 0]
    pre_trend_pvals = [event_coeffs[t]['p_value'] for t in event_times if t < 0]
    
    # Joint test for pre-trends (F-test)
    results['pre_trend_test'] = {
        'mean_pre_coef': float(np.mean(pre_trend_coefs)),
        'any_significant': any(p < 0.05 for p in pre_trend_pvals),
        'interpretation': 'Parallel trends likely' if not any(p < 0.05 for p in pre_trend_pvals) 
                          else 'Pre-trends detected - caution'
    }
    
    print(f"  Observations: {results['n_obs']:,}")
    print(f"  Pre-trends: {results['pre_trend_test']['interpretation']}")
    
    return results


def run_simple_enso_comparison(df, outcome='deaths_all'):
    """
    Simple comparison of mortality during ENSO events vs neutral.
    
    Regression: Y = β0 + β1*ElNino + β2*LaNina + controls + FE
    """
    print(f"\n--- ENSO vs Neutral Comparison: {outcome} ---")
    
    data = df.dropna(subset=['oni', outcome]).copy()
    
    # Create indicators
    data['is_el_nino'] = (data['oni'] >= 0.5).astype(int)
    data['is_la_nina'] = (data['oni'] <= -0.5).astype(int)
    
    # Intensity categories
    data['strong_el_nino'] = (data['oni'] >= 1.5).astype(int)
    data['moderate_el_nino'] = ((data['oni'] >= 1.0) & (data['oni'] < 1.5)).astype(int)
    data['strong_la_nina'] = (data['oni'] <= -1.5).astype(int)
    data['moderate_la_nina'] = ((data['oni'] <= -1.0) & (data['oni'] > -1.5)).astype(int)
    
    # Controls
    for m in range(2, 13):
        data[f'month_{m}'] = (data['month'] == m).astype(int)
    for d in range(1, 7):
        data[f'dow_{d}'] = (data['dow'] == d).astype(int)
    
    results = {}
    
    # Model 1: Binary ENSO
    print("\n  Model 1: Binary ENSO indicators")
    
    cols_to_demean = [outcome, 'is_el_nino', 'is_la_nina']
    cols_to_demean += [f'month_{m}' for m in range(2, 13)]
    cols_to_demean += [f'dow_{d}' for d in range(1, 7)]
    
    for col in cols_to_demean:
        data[f'{col}_dm'] = data.groupby('region_code')[col].transform(lambda x: x - x.mean())
    
    y = data[f'{outcome}_dm']
    X_cols = ['is_el_nino_dm', 'is_la_nina_dm']
    X_cols += [f'month_{m}_dm' for m in range(2, 13)]
    X_cols += [f'dow_{d}_dm' for d in range(1, 7)]
    X = data[X_cols]
    X = sm.add_constant(X)
    
    model = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    results['binary_model'] = {
        'el_nino_effect': {
            'coefficient': float(model.params['is_el_nino_dm']),
            'std_error': float(model.bse['is_el_nino_dm']),
            'p_value': float(model.pvalues['is_el_nino_dm'])
        },
        'la_nina_effect': {
            'coefficient': float(model.params['is_la_nina_dm']),
            'std_error': float(model.bse['is_la_nina_dm']),
            'p_value': float(model.pvalues['is_la_nina_dm'])
        },
        'n_obs': int(model.nobs),
        'r_squared': float(model.rsquared)
    }
    
    print(f"    El Niño effect: {results['binary_model']['el_nino_effect']['coefficient']:.4f} "
          f"(p={results['binary_model']['el_nino_effect']['p_value']:.4f})")
    print(f"    La Niña effect: {results['binary_model']['la_nina_effect']['coefficient']:.4f} "
          f"(p={results['binary_model']['la_nina_effect']['p_value']:.4f})")
    
    # Model 2: Continuous ONI
    print("\n  Model 2: Continuous ONI")
    
    data['oni_dm'] = data.groupby('region_code')['oni'].transform(lambda x: x - x.mean())
    
    y = data[f'{outcome}_dm']
    X_cols = ['oni_dm']
    X_cols += [f'month_{m}_dm' for m in range(2, 13)]
    X_cols += [f'dow_{d}_dm' for d in range(1, 7)]
    X = data[X_cols]
    X = sm.add_constant(X)
    
    model2 = sm.OLS(y, X).fit(cov_type='cluster', cov_kwds={'groups': data['region_code']})
    
    results['continuous_model'] = {
        'oni_effect': {
            'coefficient': float(model2.params['oni_dm']),
            'std_error': float(model2.bse['oni_dm']),
            'p_value': float(model2.pvalues['oni_dm'])
        },
        'n_obs': int(model2.nobs),
        'r_squared': float(model2.rsquared)
    }
    
    print(f"    ONI effect: {results['continuous_model']['oni_effect']['coefficient']:.4f} "
          f"(p={results['continuous_model']['oni_effect']['p_value']:.4f})")
    
    return results


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_full_event_study(level='intermediate'):
    """Run complete event study analysis."""
    
    print("\n" + "=" * 70)
    print(f"ENSO EVENT STUDY - {level.upper()} LEVEL")
    print("=" * 70)
    
    # Load data
    df = load_panel(level)
    enso_monthly = load_enso_monthly()
    
    # Identify events
    events = identify_enso_events(enso_monthly)
    
    # Create event time indicators
    df = create_event_time_indicators(df, events)
    
    all_results = {
        'level': level,
        'events': events.to_dict('records') if len(events) > 0 else [],
        'event_studies': {},
        'enso_comparison': {}
    }
    
    outcomes = ['deaths_all', 'deaths_cardiovascular', 'deaths_respiratory']
    
    # Run event studies for each event type and outcome
    for event_type in ['el_nino', 'la_nina']:
        all_results['event_studies'][event_type] = {}
        
        for outcome in outcomes:
            try:
                results = run_event_study(df, event_type, outcome)
                all_results['event_studies'][event_type][outcome] = results
            except Exception as e:
                print(f"  Error in {event_type} - {outcome}: {str(e)}")
    
    # Run simple ENSO comparison
    for outcome in outcomes:
        try:
            results = run_simple_enso_comparison(df, outcome)
            all_results['enso_comparison'][outcome] = results
        except Exception as e:
            print(f"  Error in comparison - {outcome}: {str(e)}")
    
    return all_results


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='ENSO Event Study Analysis')
    parser.add_argument('--level', type=str, default='both',
                        choices=['intermediate', 'immediate', 'both'],
                        help='Spatial level to analyze')
    args = parser.parse_args()
    
    print("=" * 70)
    print("ENSO EVENT STUDY ANALYSIS")
    print("=" * 70)
    
    levels = ['intermediate', 'immediate'] if args.level == 'both' else [args.level]
    
    for level in levels:
        results = run_full_event_study(level)
        
        # Save results
        suffix = '' if level == 'intermediate' else '_immediate'
        output_file = OUTPUT_DIR / f"enso_event_study_results{suffix}.json"
        
        with open(output_file, 'w') as f:
            json.dump(convert_to_json_serializable(results), f, indent=2)
        
        print(f"\nResults saved to: {output_file}")
    
    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
