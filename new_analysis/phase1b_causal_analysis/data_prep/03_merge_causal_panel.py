#!/usr/bin/env python3
"""
Panel Data Merge Script
=======================

Merges mortality, temperature, ENSO, and LPT data into analysis-ready panels.
Creates panels at both intermediate and immediate region levels.

Outputs:
- causal_panel_intermediate.parquet: Daily panel for intermediate regions
- causal_panel_immediate.parquet: Daily panel for immediate regions
- causal_panel_monthly_intermediate.parquet: Monthly aggregated panel
- causal_panel_monthly_immediate.parquet: Monthly aggregated panel
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = Path(__file__).resolve().parent.parent.parent.parent
PHASE0_RESULTS = BASE_DIR / "new_analysis" / "phase0_data_prep" / "results"
PHASE1B_RESULTS = Path(__file__).resolve().parent.parent / "results"
OUTPUT_DIR = PHASE1B_RESULTS


# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

def load_mortality_data(level='intermediate'):
    """Load daily mortality data."""
    print(f"Loading {level} mortality data...")
    
    df = pd.read_parquet(PHASE0_RESULTS / f"mortality_{level}_daily.parquet")
    
    # Standardize region column
    region_col = f'{level}_code'
    if region_col in df.columns:
        df = df.rename(columns={region_col: 'region_code'})
    
    df['date'] = pd.to_datetime(df['date'])
    
    print(f"  Shape: {df.shape}")
    print(f"  Date range: {df['date'].min().date()} to {df['date'].max().date()}")
    
    return df


def load_temperature_data(level='intermediate'):
    """Load daily temperature data."""
    print(f"Loading {level} temperature data...")
    
    df = pd.read_parquet(PHASE0_RESULTS / f"era5_{level}_daily.parquet")
    df['date'] = pd.to_datetime(df['date'])
    
    print(f"  Shape: {df.shape}")
    
    return df


def load_enso_data():
    """Load daily ENSO data."""
    print("Loading ENSO data...")
    
    try:
        df = pd.read_parquet(PHASE1B_RESULTS / "enso_daily.parquet")
        df['date'] = pd.to_datetime(df['date'])
        print(f"  Shape: {df.shape}")
        return df
    except FileNotFoundError:
        print("  WARNING: ENSO data not found. Run 01_prepare_enso_data.py first.")
        return None


def load_lpt_treatment(level='intermediate'):
    """Load LPT treatment data."""
    print(f"Loading {level} LPT treatment data...")
    
    try:
        df = pd.read_parquet(PHASE1B_RESULTS / f"lpt_{level}_treatment.parquet")
        print(f"  Shape: {df.shape}")
        return df
    except FileNotFoundError:
        print("  WARNING: LPT data not found. Run 02_prepare_lpt_data.py first.")
        return None


def load_lpt_panel(level='intermediate'):
    """Load LPT monthly panel."""
    print(f"Loading {level} LPT panel...")
    
    try:
        df = pd.read_parquet(PHASE1B_RESULTS / f"lpt_{level}_panel.parquet")
        df['date'] = pd.to_datetime(df['date'])
        print(f"  Shape: {df.shape}")
        return df
    except FileNotFoundError:
        print("  WARNING: LPT panel not found. Run 02_prepare_lpt_data.py first.")
        return None


def load_covariates(level='intermediate'):
    """Load regional covariates."""
    print(f"Loading {level} covariates...")
    
    try:
        df = pd.read_csv(PHASE0_RESULTS / f"ses_{level}_covariates.csv")
        
        # Standardize region column
        region_col = f'{level}_code'
        if region_col in df.columns:
            df = df.rename(columns={region_col: 'region_code'})
        
        print(f"  Shape: {df.shape}")
        return df
    except FileNotFoundError:
        print("  WARNING: Covariates not found.")
        return None


def load_pollution_data(level='intermediate'):
    """Load pollution data if available."""
    print(f"Loading {level} pollution data...")
    
    try:
        df = pd.read_parquet(PHASE0_RESULTS / f"cams_{level}_daily.parquet")
        df['date'] = pd.to_datetime(df['date'])
        
        # Standardize region column
        if 'region_code' not in df.columns:
            for col in df.columns:
                if 'code' in col.lower():
                    df = df.rename(columns={col: 'region_code'})
                    break
        
        print(f"  Shape: {df.shape}")
        return df
    except FileNotFoundError:
        print("  WARNING: Pollution data not found.")
        return None


def load_flu_data(level='intermediate'):
    """Load influenza data if available."""
    print(f"Loading {level} influenza data...")
    
    try:
        df = pd.read_parquet(PHASE0_RESULTS / f"influenza_daily_by_{level}_region.parquet")
        df['date'] = pd.to_datetime(df['date'])
        
        # Standardize region column
        region_col = f'{level}_code'
        if region_col in df.columns:
            df = df.rename(columns={region_col: 'region_code'})
        
        print(f"  Shape: {df.shape}")
        return df
    except FileNotFoundError:
        print("  WARNING: Influenza data not found.")
        return None


# =============================================================================
# PANEL CONSTRUCTION
# =============================================================================

def create_daily_panel(level='intermediate'):
    """
    Create daily panel with all variables merged.
    """
    print("\n" + "=" * 50)
    print(f"Creating daily panel for {level} level")
    print("=" * 50)
    
    # Load core data
    mortality = load_mortality_data(level)
    temperature = load_temperature_data(level)
    enso = load_enso_data()
    lpt_treatment = load_lpt_treatment(level)
    covariates = load_covariates(level)
    pollution = load_pollution_data(level)
    flu = load_flu_data(level)
    
    # Start with mortality as base
    panel = mortality.copy()
    
    # Merge temperature
    print("\nMerging temperature...")
    panel = panel.merge(
        temperature,
        on=['date', 'region_code'],
        how='left'
    )
    print(f"  After temperature merge: {len(panel):,} rows")
    
    # Merge ENSO (same for all regions)
    if enso is not None:
        print("Merging ENSO...")
        enso_cols = ['date', 'oni', 'mei', 'nino34', 'tna', 'tsa', 'whwp',
                     'enso_phase', 'enso_binary', 'is_el_nino', 'is_la_nina', 'is_strong_enso']
        enso_subset = enso[[c for c in enso_cols if c in enso.columns]]
        panel = panel.merge(enso_subset, on='date', how='left')
        print(f"  After ENSO merge: {len(panel):,} rows")
    
    # Merge LPT treatment (time-invariant characteristics)
    if lpt_treatment is not None:
        print("Merging LPT treatment...")
        lpt_cols = ['region_code', 'first_treatment_date', 'first_treatment_year',
                    'total_households_electrified', 'log_households', 'treatment_cohort',
                    'n_treated_municipalities']
        lpt_subset = lpt_treatment[[c for c in lpt_cols if c in lpt_treatment.columns]]
        panel = panel.merge(lpt_subset, on='region_code', how='left')
        
        # Create time-varying treatment indicators
        if 'first_treatment_date' in panel.columns:
            panel['post_lpt'] = (panel['date'] >= panel['first_treatment_date']).astype(int)
            panel['post_lpt'] = panel['post_lpt'].fillna(0)
            
            # Years since treatment
            panel['years_since_lpt'] = np.where(
                panel['post_lpt'] == 1,
                (panel['date'] - panel['first_treatment_date']).dt.days / 365.25,
                np.nan
            )
        
        print(f"  After LPT merge: {len(panel):,} rows")
    
    # Merge covariates (time-invariant)
    if covariates is not None:
        print("Merging covariates...")
        panel = panel.merge(covariates, on='region_code', how='left')
        print(f"  After covariates merge: {len(panel):,} rows")
    
    # Merge pollution
    if pollution is not None:
        print("Merging pollution...")
        pollution_cols = ['date', 'region_code'] + [c for c in pollution.columns 
                                                     if c not in ['date', 'region_code']]
        panel = panel.merge(pollution[pollution_cols], on=['date', 'region_code'], how='left')
        print(f"  After pollution merge: {len(panel):,} rows")
    
    # Merge flu
    if flu is not None:
        print("Merging influenza...")
        flu_cols = ['date', 'region_code'] + [c for c in flu.columns 
                                               if c not in ['date', 'region_code']]
        panel = panel.merge(flu[flu_cols], on=['date', 'region_code'], how='left')
        print(f"  After flu merge: {len(panel):,} rows")
    
    # Add time variables
    print("Adding time variables...")
    panel['year'] = panel['date'].dt.year
    panel['month'] = panel['date'].dt.month
    panel['day'] = panel['date'].dt.day
    panel['dow'] = panel['date'].dt.dayofweek
    panel['doy'] = panel['date'].dt.dayofyear
    panel['week'] = panel['date'].dt.isocalendar().week.astype(int)
    
    # Create seasonal indicators
    panel['summer'] = panel['month'].isin([12, 1, 2]).astype(int)  # Southern hemisphere
    panel['winter'] = panel['month'].isin([6, 7, 8]).astype(int)
    
    # Create temperature extremes
    if 'temp_mean' in panel.columns:
        # Calculate regional percentiles
        percentiles = panel.groupby('region_code')['temp_mean'].agg(
            ['mean', 'std', lambda x: x.quantile(0.05), lambda x: x.quantile(0.95)]
        )
        percentiles.columns = ['temp_mean_region', 'temp_std_region', 'temp_p05', 'temp_p95']
        panel = panel.merge(percentiles, on='region_code', how='left')
        
        # Extreme temperature indicators
        panel['extreme_cold'] = (panel['temp_mean'] < panel['temp_p05']).astype(int)
        panel['extreme_heat'] = (panel['temp_mean'] > panel['temp_p95']).astype(int)
        
        # Temperature anomaly
        panel['temp_anomaly'] = panel['temp_mean'] - panel['temp_mean_region']
    
    # Sort and clean
    panel = panel.sort_values(['region_code', 'date']).reset_index(drop=True)
    
    print(f"\nFinal panel shape: {panel.shape}")
    print(f"Columns: {panel.columns.tolist()}")
    
    return panel


def create_monthly_panel(daily_panel, level='intermediate'):
    """
    Create monthly aggregated panel for DiD analysis.
    """
    print(f"\nCreating monthly panel for {level} level...")
    
    # Add year-month
    daily_panel['year_month'] = daily_panel['date'].dt.to_period('M')
    
    # Define aggregation
    agg_dict = {
        'deaths_all': 'sum',
        'deaths_respiratory': 'sum',
        'deaths_cardiovascular': 'sum',
    }
    
    # Add temperature aggregations if available
    if 'temp_mean' in daily_panel.columns:
        agg_dict['temp_mean'] = 'mean'
        agg_dict['extreme_cold'] = 'sum'
        agg_dict['extreme_heat'] = 'sum'
    
    # Add ENSO if available
    if 'oni' in daily_panel.columns:
        agg_dict['oni'] = 'mean'
        agg_dict['is_el_nino'] = 'max'
        agg_dict['is_la_nina'] = 'max'
    
    # Add pollution if available
    for col in ['pm25', 'pm10', 'o3', 'no2']:
        if col in daily_panel.columns:
            agg_dict[col] = 'mean'
    
    # Aggregate
    monthly = daily_panel.groupby(['region_code', 'year_month', 'year', 'month']).agg(agg_dict).reset_index()
    
    # Convert year_month back to datetime
    monthly['date'] = monthly['year_month'].dt.to_timestamp()
    monthly = monthly.drop(columns=['year_month'])
    
    # Add time-invariant variables
    time_invariant = ['first_treatment_date', 'first_treatment_year', 'treatment_cohort',
                      'total_households_electrified', 'log_households']
    
    existing_invariant = [c for c in time_invariant if c in daily_panel.columns]
    if existing_invariant:
        invariant_df = daily_panel.groupby('region_code')[existing_invariant].first().reset_index()
        monthly = monthly.merge(invariant_df, on='region_code', how='left')
    
    # Re-calculate post treatment indicator
    if 'first_treatment_date' in monthly.columns:
        monthly['post_lpt'] = (monthly['date'] >= monthly['first_treatment_date']).astype(int)
        monthly['post_lpt'] = monthly['post_lpt'].fillna(0)
    
    # Add covariates
    covariate_cols = ['population', 'gdp_pc', 'hdi', 'urban_rate', 'elderly_share']
    existing_covars = [c for c in covariate_cols if c in daily_panel.columns]
    if existing_covars:
        covar_df = daily_panel.groupby('region_code')[existing_covars].first().reset_index()
        monthly = monthly.merge(covar_df, on='region_code', how='left')
    
    # Calculate mortality rates
    if 'population' in monthly.columns:
        monthly['death_rate'] = monthly['deaths_all'] / monthly['population'] * 100000
    
    # Sort
    monthly = monthly.sort_values(['region_code', 'date']).reset_index(drop=True)
    
    print(f"  Monthly panel shape: {monthly.shape}")
    
    return monthly


def main():
    """Main function to create all panels."""
    parser = argparse.ArgumentParser(description='Create causal analysis panels')
    parser.add_argument('--level', type=str, default='both',
                        choices=['intermediate', 'immediate', 'both'],
                        help='Spatial level to process')
    args = parser.parse_args()
    
    print("=" * 70)
    print("CAUSAL ANALYSIS PANEL CREATION")
    print("=" * 70)
    
    levels = ['intermediate', 'immediate'] if args.level == 'both' else [args.level]
    
    for level in levels:
        print(f"\n{'='*70}")
        print(f"Processing {level.upper()} level")
        print("=" * 70)
        
        # Create daily panel
        daily_panel = create_daily_panel(level)
        
        # Create monthly panel
        monthly_panel = create_monthly_panel(daily_panel, level)
        
        # Save outputs
        print(f"\nSaving {level} panels...")
        daily_panel.to_parquet(OUTPUT_DIR / f"causal_panel_{level}.parquet", index=False)
        monthly_panel.to_parquet(OUTPUT_DIR / f"causal_panel_monthly_{level}.parquet", index=False)
        
        # Save sample CSV for inspection
        daily_panel.head(1000).to_csv(OUTPUT_DIR / f"causal_panel_{level}_sample.csv", index=False)
        
        print(f"  Saved: causal_panel_{level}.parquet ({daily_panel.shape})")
        print(f"  Saved: causal_panel_monthly_{level}.parquet ({monthly_panel.shape})")
    
    print("\n" + "=" * 70)
    print("Panel creation complete!")
    print(f"Output files saved to: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
