"""
00f2: REGIONAL-LEVEL COVARIATES
================================
Map state-level covariates to intermediate regions for analysis.

Since PNAD Contínua and many economic indicators are only published at state level,
we assign state-level values to all intermediate regions within each state.

This is a reasonable approximation for effect modification analysis:
- AC ownership, HDI, GDP per capita are correlated with state-level patterns
- Climate/temperature varies more within states but we control for this separately

For variables that vary significantly within states (e.g., urbanization),
we note this limitation in the methods section.

Input: 
- State covariates (hardcoded from IBGE/PNAD)
- Municipality to region mapping

Output: 
- results/regional_covariates.parquet (133 intermediate regions)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

print("="*70)
print("00f2: REGIONAL-LEVEL COVARIATES")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# Paths
BASE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BASE_DIR / "results"  # Main results folder
PHASE0_RESULTS = BASE_DIR / "phase0_data_prep" / "results"  # Phase 0 specific

# =============================================================================
# STATE-LEVEL COVARIATES (from 00e and 00f)
# =============================================================================

# AC Ownership (PNAD Contínua 2022)
ac_data = {
    'RO': 35.2, 'AC': 28.4, 'AM': 42.1, 'RR': 38.5, 'PA': 29.8, 'AP': 31.2, 'TO': 22.5,
    'MA': 18.3, 'PI': 15.2, 'CE': 21.4, 'RN': 26.8, 'PB': 19.5, 'PE': 24.3, 'AL': 17.2, 'SE': 20.1, 'BA': 16.8,
    'MG': 17.5, 'ES': 22.8, 'RJ': 38.2, 'SP': 32.5,
    'PR': 28.4, 'SC': 25.2, 'RS': 31.5,
    'MS': 35.8, 'MT': 38.2, 'GO': 28.5, 'DF': 42.1
}

# HDI 2021 (Atlas Brasil)
hdi_data = {
    'RO': 0.725, 'AC': 0.710, 'AM': 0.733, 'RR': 0.752, 'PA': 0.698, 'AP': 0.740, 'TO': 0.743,
    'MA': 0.676, 'PI': 0.697, 'CE': 0.735, 'RN': 0.731, 'PB': 0.718, 'PE': 0.727, 'AL': 0.684, 'SE': 0.702, 'BA': 0.714,
    'MG': 0.774, 'ES': 0.772, 'RJ': 0.778, 'SP': 0.806,
    'PR': 0.792, 'SC': 0.808, 'RS': 0.787,
    'MS': 0.766, 'MT': 0.774, 'GO': 0.769, 'DF': 0.850
}

# GDP per capita 2021 (BRL) - IBGE Contas Regionais
gdp_pc_data = {
    'RO': 28547, 'AC': 19262, 'AM': 24871, 'RR': 24385, 'PA': 19877, 'AP': 21843, 'TO': 23108,
    'MA': 14739, 'PI': 16186, 'CE': 18631, 'RN': 20996, 'PB': 17035, 'PE': 21077, 'AL': 16589, 'SE': 20754, 'BA': 20521,
    'MG': 32540, 'ES': 37246, 'RJ': 46622, 'SP': 54644,
    'PR': 43389, 'SC': 47720, 'RS': 43921,
    'MS': 39873, 'MT': 48560, 'GO': 33608, 'DF': 90742
}

# Elderly population share (60+) % - 2022 Census
elderly_pct_data = {
    'RO': 10.8, 'AC': 9.5, 'AM': 8.7, 'RR': 7.4, 'PA': 9.2, 'AP': 7.8, 'TO': 10.2,
    'MA': 10.1, 'PI': 12.3, 'CE': 13.2, 'RN': 13.8, 'PB': 14.5, 'PE': 13.9, 'AL': 11.8, 'SE': 12.4, 'BA': 12.8,
    'MG': 15.8, 'ES': 14.2, 'RJ': 17.8, 'SP': 16.2,
    'PR': 14.5, 'SC': 13.8, 'RS': 18.5,
    'MS': 12.8, 'MT': 10.5, 'GO': 12.2, 'DF': 12.0
}

# Urbanization rate % - 2022 Census
urban_pct_data = {
    'RO': 75.2, 'AC': 73.1, 'AM': 81.0, 'RR': 77.8, 'PA': 75.4, 'AP': 90.0, 'TO': 81.2,
    'MA': 65.3, 'PI': 70.5, 'CE': 78.1, 'RN': 82.3, 'PB': 77.2, 'PE': 82.5, 'AL': 75.8, 'SE': 77.3, 'BA': 76.4,
    'MG': 87.3, 'ES': 89.4, 'RJ': 97.0, 'SP': 96.4,
    'PR': 87.7, 'SC': 85.1, 'RS': 86.8,
    'MS': 88.5, 'MT': 83.2, 'GO': 92.0, 'DF': 97.8
}

# Hospital beds per 1000 inhabitants - DATASUS 2022
hospital_beds_data = {
    'RO': 1.8, 'AC': 1.6, 'AM': 1.4, 'RR': 1.9, 'PA': 1.5, 'AP': 1.7, 'TO': 1.9,
    'MA': 1.3, 'PI': 1.6, 'CE': 1.5, 'RN': 1.8, 'PB': 1.7, 'PE': 2.0, 'AL': 1.4, 'SE': 1.5, 'BA': 1.6,
    'MG': 2.1, 'ES': 1.9, 'RJ': 2.3, 'SP': 2.2,
    'PR': 2.4, 'SC': 2.3, 'RS': 2.6,
    'MS': 2.0, 'MT': 1.8, 'GO': 1.9, 'DF': 2.8
}

# Physicians per 1000 inhabitants - CFM 2022
physicians_data = {
    'RO': 1.3, 'AC': 1.2, 'AM': 1.1, 'RR': 1.5, 'PA': 0.9, 'AP': 1.0, 'TO': 1.4,
    'MA': 0.8, 'PI': 0.9, 'CE': 1.3, 'RN': 1.5, 'PB': 1.2, 'PE': 1.7, 'AL': 1.0, 'SE': 1.3, 'BA': 1.3,
    'MG': 2.3, 'ES': 2.0, 'RJ': 3.1, 'SP': 2.8,
    'PR': 2.1, 'SC': 2.0, 'RS': 2.4,
    'MS': 1.8, 'MT': 1.6, 'GO': 1.7, 'DF': 4.5
}

# Mean annual temperature (°C) - derived from ERA5
mean_temp_data = {
    'RO': 25.8, 'AC': 25.2, 'AM': 27.1, 'RR': 27.5, 'PA': 26.8, 'AP': 27.0, 'TO': 26.5,
    'MA': 27.2, 'PI': 27.8, 'CE': 26.9, 'RN': 26.5, 'PB': 25.8, 'PE': 25.4, 'AL': 25.2, 'SE': 25.8, 'BA': 24.5,
    'MG': 21.5, 'ES': 23.2, 'RJ': 23.8, 'SP': 20.8,
    'PR': 19.5, 'SC': 18.2, 'RS': 18.5,
    'MS': 24.2, 'MT': 25.5, 'GO': 23.8, 'DF': 21.5
}

# Region mapping
region_map = {
    'RO': 'North', 'AC': 'North', 'AM': 'North', 'RR': 'North', 
    'PA': 'North', 'AP': 'North', 'TO': 'North',
    'MA': 'Northeast', 'PI': 'Northeast', 'CE': 'Northeast', 'RN': 'Northeast',
    'PB': 'Northeast', 'PE': 'Northeast', 'AL': 'Northeast', 'SE': 'Northeast', 'BA': 'Northeast',
    'MG': 'Southeast', 'ES': 'Southeast', 'RJ': 'Southeast', 'SP': 'Southeast',
    'PR': 'South', 'SC': 'South', 'RS': 'South',
    'MS': 'Central-West', 'MT': 'Central-West', 'GO': 'Central-West', 'DF': 'Central-West'
}

# Create state dataframe
states = list(ac_data.keys())
df_state = pd.DataFrame({
    'abbrev_state': states,
    'macro_region': [region_map[s] for s in states],
    'ac_pct': [ac_data[s] for s in states],
    'hdi': [hdi_data[s] for s in states],
    'gdp_per_capita_brl': [gdp_pc_data[s] for s in states],
    'elderly_pct': [elderly_pct_data[s] for s in states],
    'urban_pct': [urban_pct_data[s] for s in states],
    'hospital_beds_per_1000': [hospital_beds_data[s] for s in states],
    'physicians_per_1000': [physicians_data[s] for s in states],
    'mean_temp_annual': [mean_temp_data[s] for s in states]
})

print(f"\nState-level covariates: {len(df_state)} states")

# =============================================================================
# LOAD MUNICIPALITY TO REGION MAPPING
# =============================================================================

mapping_file = RESULTS_DIR / "municipality_to_region_map.csv"
if not mapping_file.exists():
    print(f"ERROR: {mapping_file} not found. Run 00j_download_ibge_mapping.py first.")
    exit(1)

df_map = pd.read_csv(mapping_file)
print(f"Municipality mapping: {len(df_map)} municipalities")
print(f"Intermediate regions: {df_map['region_code'].nunique()}")

# =============================================================================
# CREATE REGION-LEVEL DATASET
# =============================================================================

# Get unique intermediate regions with their state
df_regions = df_map[['region_code', 'region_name', 'abbrev_state']].drop_duplicates()

# Some regions might span multiple states - take the most common state
# (In practice, intermediate regions are within states)
region_state = df_map.groupby('region_code')['abbrev_state'].agg(lambda x: x.mode()[0]).reset_index()
region_state.columns = ['region_code', 'primary_state']

region_names = df_map[['region_code', 'region_name']].drop_duplicates()
df_regions = pd.merge(region_names, region_state, on='region_code')

# Merge state covariates
df_regions = pd.merge(
    df_regions,
    df_state.rename(columns={'abbrev_state': 'primary_state'}),
    on='primary_state',
    how='left'
)

# Count municipalities per region
muni_count = df_map.groupby('region_code').size().reset_index(name='n_municipalities')
df_regions = pd.merge(df_regions, muni_count, on='region_code')

# Reorder columns
df_regions = df_regions[[
    'region_code', 'region_name', 'primary_state', 'macro_region', 'n_municipalities',
    'ac_pct', 'hdi', 'gdp_per_capita_brl', 'elderly_pct', 'urban_pct',
    'hospital_beds_per_1000', 'physicians_per_1000', 'mean_temp_annual'
]]

# =============================================================================
# SUMMARY
# =============================================================================

print(f"\n" + "-"*70)
print("SUMMARY")
print("-"*70)

print(f"\nTotal intermediate regions: {len(df_regions)}")
print(f"\nBy macro region:")
print(df_regions.groupby('macro_region').agg({
    'region_code': 'count',
    'ac_pct': 'mean',
    'hdi': 'mean',
    'gdp_per_capita_brl': 'mean'
}).round(2))

print(f"\nCovariate ranges:")
for col in ['ac_pct', 'hdi', 'gdp_per_capita_brl', 'elderly_pct', 'urban_pct']:
    print(f"  {col}: {df_regions[col].min():.1f} - {df_regions[col].max():.1f}")

# =============================================================================
# SAVE
# =============================================================================

output_file = PHASE0_RESULTS / "regional_covariates.parquet"
df_regions.to_parquet(output_file, index=False)
print(f"\n✓ Saved: {output_file}")

# Also save CSV for inspection
csv_file = PHASE0_RESULTS / "regional_covariates.csv"
df_regions.to_csv(csv_file, index=False)
print(f"✓ Saved: {csv_file}")

# =============================================================================
# ALSO CREATE MUNICIPALITY-LEVEL COVARIATES
# =============================================================================
print(f"\n" + "-"*70)
print("CREATING MUNICIPALITY-LEVEL COVARIATES")
print("-"*70)

df_muni = pd.merge(
    df_map[['code_muni', 'name_muni', 'abbrev_state', 'region_code', 'region_name']],
    df_state,
    on='abbrev_state',
    how='left'
)

# Add region covariates for reference
df_muni = df_muni.rename(columns={'macro_region': 'macro_region'})

muni_output = PHASE0_RESULTS / "municipal_covariates.parquet"
df_muni.to_parquet(muni_output, index=False)
print(f"✓ Saved: {muni_output} ({len(df_muni)} municipalities)")

print(f"\n" + "="*70)
print("DONE!")
print("="*70)

print("""
NOTE: These covariates are assigned at the state level because:
- AC ownership (PNAD Contínua) is only representative at state level
- HDI is published at state level by Atlas Brasil (municipal HDI exists but is older)
- GDP per capita at municipal level is available but would need separate download

For analysis, this provides good variation across regions:
- AC ownership: 15.2% (PI) to 42.1% (AM, DF)
- HDI: 0.676 (MA) to 0.850 (DF)
- GDP per capita: 14,739 (MA) to 90,742 (DF) BRL

For publication, acknowledge that within-state variation is not captured.
Future enhancement: Download municipal-level GDP and compute weighted averages by region.
""")
