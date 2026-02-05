#!/usr/bin/env python3
"""
Luz para Todos (LPT) Data Preparation Script
=============================================

Prepares the Luz para Todos rural electrification program data for causal analysis.
The program provides staggered treatment timing across municipalities for DiD design.

Key Variables:
- Treatment date: First electrification date per municipality
- Treatment intensity: Cumulative households electrified
- Program phases: Different implementation waves

Outputs:
- lpt_municipality_treatment.parquet: Municipality-level treatment timing
- lpt_municipality_panel.parquet: Monthly panel of cumulative electrification
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = Path(__file__).resolve().parent.parent.parent.parent
INPUT_DIR = BASE_DIR / "Input_data" / "data_causal_analysis"
RESULTS_DIR = BASE_DIR / "new_analysis" / "phase0_data_prep" / "results"
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# MAIN FUNCTIONS
# =============================================================================

def load_lpt_data():
    """Load Luz para Todos program data."""
    print("Loading Luz para Todos data...")
    
    df = pd.read_csv(
        INPUT_DIR / "domicilios_atendidos.csv",
        sep=';',
        encoding='latin1'
    )
    
    print(f"  Raw records: {len(df):,}")
    print(f"  Columns: {df.columns.tolist()}")
    
    # Rename columns for clarity
    df.columns = ['programa', 'qtd_domicilios', 'mes', 'ano', 'municipio', 'estado', 'dt_homologacao']
    
    # Parse homologation date
    df['dt_homologacao'] = pd.to_datetime(df['dt_homologacao'], format='%d/%m/%Y', errors='coerce')
    
    # Create year-month date
    df['year'] = df['ano'].astype(int)
    df['month'] = df['mes'].astype(int)
    df['date'] = pd.to_datetime(df[['year', 'month']].assign(day=1))
    
    # Clean municipality names
    df['municipio'] = df['municipio'].str.strip()
    df['estado'] = df['estado'].str.strip()
    
    # Create municipality identifier
    df['muni_state'] = df['municipio'] + '_' + df['estado']
    
    print(f"  Programs: {df['programa'].unique()}")
    print(f"  Date range: {df['date'].min().date()} to {df['date'].max().date()}")
    print(f"  Unique municipalities: {df['muni_state'].nunique():,}")
    
    return df


def load_municipality_mapping():
    """Load municipality to region mapping."""
    print("\nLoading municipality mapping...")
    
    mapping = pd.read_csv(RESULTS_DIR / "municipality_to_all_regions_map.csv")
    
    # State name to abbreviation mapping
    state_abbrev = {
        'Acre': 'AC', 'Alagoas': 'AL', 'Amapá': 'AP', 'Amazonas': 'AM',
        'Bahia': 'BA', 'Ceará': 'CE', 'Distrito Federal': 'DF',
        'Espírito Santo': 'ES', 'Goiás': 'GO', 'Maranhão': 'MA',
        'Mato Grosso': 'MT', 'Mato Grosso do Sul': 'MS', 'Minas Gerais': 'MG',
        'Pará': 'PA', 'Paraíba': 'PB', 'Paraná': 'PR', 'Pernambuco': 'PE',
        'Piauí': 'PI', 'Rio de Janeiro': 'RJ', 'Rio Grande do Norte': 'RN',
        'Rio Grande do Sul': 'RS', 'Rondônia': 'RO', 'Roraima': 'RR',
        'Santa Catarina': 'SC', 'São Paulo': 'SP', 'Sergipe': 'SE', 'Tocantins': 'TO'
    }
    mapping['state_abbrev_to_name'] = {v: k for k, v in state_abbrev.items()}
    
    # Create name-state identifier using abbreviation
    mapping['muni_state'] = mapping['name_muni'].str.strip() + '_' + mapping['abbrev_state'].str.strip()
    
    # Also create normalized version for fuzzy matching
    mapping['muni_normalized'] = (
        mapping['name_muni'].str.strip()
        .str.upper()
        .str.replace("'", "", regex=False)
        .str.replace("-", " ", regex=False)
    )
    
    print(f"  Municipalities in mapping: {len(mapping):,}")
    
    return mapping, state_abbrev


def create_municipality_treatment(lpt_df, mapping, state_abbrev):
    """
    Create municipality-level treatment variables.
    
    Treatment defined as first electrification under Luz para Todos.
    """
    print("\nCreating municipality treatment data...")
    
    # Focus on main LPT program
    lpt_rural = lpt_df[lpt_df['programa'].str.contains('LPT', na=False)].copy()
    print(f"  LPT records: {len(lpt_rural):,}")
    
    # Convert state names to abbreviations in LPT data
    lpt_rural['state_abbrev'] = lpt_rural['estado'].map(state_abbrev)
    
    # Create matching key
    lpt_rural['muni_normalized'] = (
        lpt_rural['municipio'].str.strip()
        .str.upper()
        .str.replace("'", "", regex=False)
        .str.replace("-", " ", regex=False)
    )
    lpt_rural['match_key'] = lpt_rural['muni_normalized'] + '_' + lpt_rural['state_abbrev'].fillna('')
    
    # Create mapping match key
    mapping['match_key'] = mapping['muni_normalized'] + '_' + mapping['abbrev_state']
    
    # Aggregate by municipality and date
    muni_monthly = lpt_rural.groupby(['match_key', 'municipio', 'estado', 'date']).agg({
        'qtd_domicilios': 'sum'
    }).reset_index()
    
    # Calculate first treatment date per municipality
    first_treatment = muni_monthly.groupby('match_key').agg({
        'date': 'min',
        'municipio': 'first',
        'estado': 'first'
    }).reset_index()
    first_treatment.columns = ['match_key', 'first_treatment_date', 'municipio', 'estado']
    
    # Calculate total households and intensity
    total_hh = lpt_rural.groupby('match_key')['qtd_domicilios'].sum().reset_index()
    total_hh.columns = ['match_key', 'total_households_electrified']
    
    # Merge
    treatment = first_treatment.merge(total_hh, on='match_key', how='left')
    
    # Match to municipality codes
    treatment = treatment.merge(
        mapping[['match_key', 'code_muni', 'intermediate_code', 'immediate_code']],
        on='match_key',
        how='left'
    )
    
    print(f"  Municipalities with treatment data: {len(treatment):,}")
    print(f"  Matched to IBGE codes: {treatment['code_muni'].notna().sum():,}")
    
    # Create treatment timing variables
    treatment['first_treatment_year'] = treatment['first_treatment_date'].dt.year
    treatment['first_treatment_month'] = treatment['first_treatment_date'].dt.month
    treatment['first_treatment_ym'] = (
        treatment['first_treatment_year'] * 12 + treatment['first_treatment_month']
    )
    
    # Create cohort groups (year of first treatment)
    treatment['treatment_cohort'] = treatment['first_treatment_year']
    
    # Calculate treatment intensity (log households)
    treatment['log_households'] = np.log1p(treatment['total_households_electrified'])
    
    # Create early/late treatment indicator
    median_date = treatment['first_treatment_date'].median()
    treatment['early_treatment'] = (treatment['first_treatment_date'] <= median_date).astype(int)
    
    print(f"  Treatment cohorts: {treatment['treatment_cohort'].value_counts().head(10).to_dict()}")
    
    return treatment


def create_municipality_panel(lpt_df, mapping, state_abbrev):
    """
    Create monthly panel of cumulative electrification by municipality.
    For time-varying treatment intensity analysis.
    """
    print("\nCreating municipality monthly panel...")
    
    # Focus on LPT program
    lpt_rural = lpt_df[lpt_df['programa'].str.contains('LPT', na=False)].copy()
    
    # Convert state names to abbreviations
    lpt_rural['state_abbrev'] = lpt_rural['estado'].map(state_abbrev)
    lpt_rural['muni_normalized'] = (
        lpt_rural['municipio'].str.strip()
        .str.upper()
        .str.replace("'", "", regex=False)
        .str.replace("-", " ", regex=False)
    )
    lpt_rural['match_key'] = lpt_rural['muni_normalized'] + '_' + lpt_rural['state_abbrev'].fillna('')
    
    # Create full panel of municipality x months
    municipalities = lpt_rural['match_key'].unique()
    months = pd.date_range(start='2004-01-01', end='2024-12-01', freq='MS')
    
    panel = pd.MultiIndex.from_product(
        [municipalities, months],
        names=['match_key', 'date']
    ).to_frame(index=False)
    
    # Aggregate monthly households
    monthly_hh = lpt_rural.groupby(['match_key', 'date'])['qtd_domicilios'].sum().reset_index()
    
    # Merge to panel
    panel = panel.merge(monthly_hh, on=['match_key', 'date'], how='left')
    panel['qtd_domicilios'] = panel['qtd_domicilios'].fillna(0)
    
    # Calculate cumulative households
    panel = panel.sort_values(['match_key', 'date'])
    panel['cumulative_households'] = panel.groupby('match_key')['qtd_domicilios'].cumsum()
    
    # Create treatment indicator (ever treated)
    panel['ever_treated'] = (panel['cumulative_households'] > 0).astype(int)
    
    # Calculate log cumulative households
    panel['log_cumulative_hh'] = np.log1p(panel['cumulative_households'])
    
    # Add year/month
    panel['year'] = panel['date'].dt.year
    panel['month'] = panel['date'].dt.month
    
    # Match to municipality codes
    mapping['match_key'] = mapping['muni_normalized'] + '_' + mapping['abbrev_state']
    muni_info = mapping[['match_key', 'code_muni', 'intermediate_code', 'immediate_code']].drop_duplicates()
    panel = panel.merge(muni_info, on='match_key', how='left')
    
    print(f"  Panel observations: {len(panel):,}")
    print(f"  Date range: {panel['date'].min().date()} to {panel['date'].max().date()}")
    
    return panel


def aggregate_to_regions(treatment, level='intermediate'):
    """
    Aggregate treatment data to intermediate or immediate region level.
    
    For each region:
    - First treatment date: earliest treatment in any municipality
    - Treatment intensity: sum of households electrified
    - Coverage: proportion of municipalities treated
    """
    print(f"\nAggregating to {level} region level...")
    
    code_col = f'{level}_code'
    
    # Filter to matched municipalities
    matched = treatment[treatment[code_col].notna()].copy()
    
    # Aggregate by region
    region_treatment = matched.groupby(code_col).agg({
        'first_treatment_date': 'min',
        'total_households_electrified': 'sum',
        'code_muni': 'count'  # Number of treated municipalities
    }).reset_index()
    
    region_treatment.columns = [
        code_col,
        'first_treatment_date',
        'total_households_electrified',
        'n_treated_municipalities'
    ]
    
    # Calculate region-level variables
    region_treatment['first_treatment_year'] = region_treatment['first_treatment_date'].dt.year
    region_treatment['log_households'] = np.log1p(region_treatment['total_households_electrified'])
    
    # Create treatment cohort
    region_treatment['treatment_cohort'] = region_treatment['first_treatment_year']
    
    # Rename for consistency
    region_treatment = region_treatment.rename(columns={code_col: 'region_code'})
    
    print(f"  Regions with treatment: {len(region_treatment):,}")
    
    return region_treatment


def create_region_panel(muni_panel, level='intermediate'):
    """
    Create regional panel from municipality panel.
    Aggregates treatment intensity over time.
    """
    print(f"\nCreating {level} region panel...")
    
    code_col = f'{level}_code'
    
    # Filter to matched municipalities
    matched = muni_panel[muni_panel[code_col].notna()].copy()
    
    # Aggregate by region and date
    region_panel = matched.groupby([code_col, 'date', 'year', 'month']).agg({
        'qtd_domicilios': 'sum',
        'cumulative_households': 'sum',
        'ever_treated': 'max',  # 1 if any municipality treated
        'code_muni': 'count'  # Number of municipalities in panel
    }).reset_index()
    
    region_panel.columns = [
        code_col, 'date', 'year', 'month',
        'new_households', 'cumulative_households',
        'region_ever_treated', 'n_municipalities'
    ]
    
    # Calculate log cumulative households
    region_panel['log_cumulative_hh'] = np.log1p(region_panel['cumulative_households'])
    
    # Rename for consistency
    region_panel = region_panel.rename(columns={code_col: 'region_code'})
    
    print(f"  Panel observations: {len(region_panel):,}")
    
    return region_panel


def main():
    """Main function to prepare LPT data."""
    print("=" * 70)
    print("LUZ PARA TODOS DATA PREPARATION")
    print("=" * 70)
    
    # Load data
    lpt_df = load_lpt_data()
    mapping, state_abbrev = load_municipality_mapping()
    
    # Create municipality-level treatment
    muni_treatment = create_municipality_treatment(lpt_df, mapping, state_abbrev)
    
    # Create municipality panel
    muni_panel = create_municipality_panel(lpt_df, mapping, state_abbrev)
    
    # Create regional aggregations
    print("\n" + "-" * 50)
    print("Creating regional aggregations...")
    
    # Intermediate level
    intermediate_treatment = aggregate_to_regions(muni_treatment, level='intermediate')
    intermediate_panel = create_region_panel(muni_panel, level='intermediate')
    
    # Immediate level
    immediate_treatment = aggregate_to_regions(muni_treatment, level='immediate')
    immediate_panel = create_region_panel(muni_panel, level='immediate')
    
    # Save outputs
    print("\n" + "-" * 50)
    print("Saving outputs...")
    
    # Municipality level
    muni_treatment.to_parquet(OUTPUT_DIR / "lpt_municipality_treatment.parquet", index=False)
    muni_panel.to_parquet(OUTPUT_DIR / "lpt_municipality_panel.parquet", index=False)
    
    # Intermediate level
    intermediate_treatment.to_parquet(OUTPUT_DIR / "lpt_intermediate_treatment.parquet", index=False)
    intermediate_treatment.to_csv(OUTPUT_DIR / "lpt_intermediate_treatment.csv", index=False)
    intermediate_panel.to_parquet(OUTPUT_DIR / "lpt_intermediate_panel.parquet", index=False)
    
    # Immediate level
    immediate_treatment.to_parquet(OUTPUT_DIR / "lpt_immediate_treatment.parquet", index=False)
    immediate_treatment.to_csv(OUTPUT_DIR / "lpt_immediate_treatment.csv", index=False)
    immediate_panel.to_parquet(OUTPUT_DIR / "lpt_immediate_panel.parquet", index=False)
    
    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    
    print("\nMunicipality Treatment Data:")
    print(f"  Total municipalities: {len(muni_treatment):,}")
    print(f"  Matched to IBGE codes: {muni_treatment['code_muni'].notna().sum():,}")
    print(f"  Total households electrified: {muni_treatment['total_households_electrified'].sum():,.0f}")
    
    print("\nTreatment Cohort Distribution:")
    print(muni_treatment['treatment_cohort'].value_counts().sort_index().head(15))
    
    print("\nIntermediate Region Treatment:")
    print(f"  Regions treated: {len(intermediate_treatment):,}")
    print(f"  Cohort range: {intermediate_treatment['treatment_cohort'].min()} - {intermediate_treatment['treatment_cohort'].max()}")
    
    print("\nImmediate Region Treatment:")
    print(f"  Regions treated: {len(immediate_treatment):,}")
    print(f"  Cohort range: {immediate_treatment['treatment_cohort'].min()} - {immediate_treatment['treatment_cohort'].max()}")
    
    print("\n" + "=" * 70)
    print("LPT data preparation complete!")
    print(f"Output files saved to: {OUTPUT_DIR}")
    print("=" * 70)
    
    return {
        'municipality_treatment': muni_treatment,
        'municipality_panel': muni_panel,
        'intermediate_treatment': intermediate_treatment,
        'intermediate_panel': intermediate_panel,
        'immediate_treatment': immediate_treatment,
        'immediate_panel': immediate_panel
    }


if __name__ == "__main__":
    results = main()
