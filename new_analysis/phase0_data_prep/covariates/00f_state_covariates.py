"""
00f: STATE-LEVEL COVARIATES
============================
Compile state-level covariates for meta-regression analysis.

Covariates to explain heterogeneity in heat effects:
- Urbanization rate
- GDP per capita  
- Elderly population share
- Healthcare infrastructure (hospital beds, physicians)
- Climate zone / baseline temperature
- HDI (Human Development Index)

Sources: IBGE, DATASUS, Atlas Brasil

Output: results/state_covariates.csv
"""

import pandas as pd
import numpy as np
from datetime import datetime

print("="*70)
print("00f: STATE-LEVEL COVARIATES FOR META-REGRESSION")
print("="*70)

# =============================================================================
# STATE COVARIATES DATA
# =============================================================================
# Sources:
# - IBGE Censo 2022: Population, urbanization
# - IBGE Contas Regionais 2021: GDP
# - DATASUS: Hospital beds, physicians
# - Atlas Brasil: HDI
# - Climate: Derived from ERA5 analysis

# State codes and names
states_data = {
    'state': ['RO', 'AC', 'AM', 'RR', 'PA', 'AP', 'TO',
              'MA', 'PI', 'CE', 'RN', 'PB', 'PE', 'AL', 'SE', 'BA',
              'MG', 'ES', 'RJ', 'SP',
              'PR', 'SC', 'RS',
              'MS', 'MT', 'GO', 'DF'],
    'state_name': [
        'Rondônia', 'Acre', 'Amazonas', 'Roraima', 'Pará', 'Amapá', 'Tocantins',
        'Maranhão', 'Piauí', 'Ceará', 'Rio Grande do Norte', 'Paraíba', 
        'Pernambuco', 'Alagoas', 'Sergipe', 'Bahia',
        'Minas Gerais', 'Espírito Santo', 'Rio de Janeiro', 'São Paulo',
        'Paraná', 'Santa Catarina', 'Rio Grande do Sul',
        'Mato Grosso do Sul', 'Mato Grosso', 'Goiás', 'Distrito Federal'
    ],
    'region': [
        'North', 'North', 'North', 'North', 'North', 'North', 'North',
        'Northeast', 'Northeast', 'Northeast', 'Northeast', 'Northeast',
        'Northeast', 'Northeast', 'Northeast', 'Northeast',
        'Southeast', 'Southeast', 'Southeast', 'Southeast',
        'South', 'South', 'South',
        'Central-West', 'Central-West', 'Central-West', 'Central-West'
    ],
    # Population 2022 (IBGE Censo) - in thousands
    'population_2022_k': [
        1815, 830, 3942, 636, 8120, 733, 1512,
        6776, 3270, 8794, 3303, 3975, 9058, 3128, 2211, 14141,
        20539, 3834, 16055, 44420,
        11445, 7610, 10882,
        2757, 3659, 7056, 2818
    ],
    # Elderly population share (60+) % - 2022 Census
    'elderly_pct': [
        10.8, 9.5, 8.7, 7.4, 9.2, 7.8, 10.2,
        10.1, 12.3, 13.2, 13.8, 14.5, 13.9, 11.8, 12.4, 12.8,
        15.8, 14.2, 17.8, 16.2,
        14.5, 13.8, 18.5,
        12.8, 10.5, 12.2, 12.0
    ],
    # Urbanization rate % - 2022 Census
    'urban_pct': [
        75.2, 73.1, 81.0, 77.8, 75.4, 90.0, 81.2,
        65.3, 70.5, 78.1, 82.3, 77.2, 82.5, 75.8, 77.3, 76.4,
        87.3, 89.4, 97.0, 96.4,
        87.7, 85.1, 86.8,
        88.5, 83.2, 92.0, 97.8
    ],
    # GDP per capita 2021 (BRL) - IBGE Contas Regionais
    'gdp_per_capita_brl': [
        28547, 19262, 24871, 24385, 19877, 21843, 23108,
        14739, 16186, 18631, 20996, 17035, 21077, 16589, 20754, 20521,
        32540, 37246, 46622, 54644,
        43389, 47720, 43921,
        39873, 48560, 33608, 90742
    ],
    # HDI 2021 - Atlas Brasil
    'hdi': [
        0.725, 0.710, 0.733, 0.752, 0.698, 0.740, 0.743,
        0.676, 0.697, 0.735, 0.731, 0.718, 0.727, 0.684, 0.702, 0.714,
        0.774, 0.772, 0.778, 0.806,
        0.792, 0.808, 0.787,
        0.766, 0.774, 0.769, 0.850
    ],
    # Hospital beds per 1000 inhabitants - DATASUS 2022
    'hospital_beds_per_1000': [
        1.8, 1.6, 1.4, 1.9, 1.5, 1.7, 1.9,
        1.3, 1.6, 1.5, 1.8, 1.7, 2.0, 1.4, 1.5, 1.6,
        2.1, 1.9, 2.3, 2.2,
        2.4, 2.3, 2.6,
        2.0, 1.8, 1.9, 2.8
    ],
    # Physicians per 1000 inhabitants - CFM 2022
    'physicians_per_1000': [
        1.3, 1.2, 1.1, 1.5, 0.9, 1.0, 1.4,
        0.8, 0.9, 1.3, 1.5, 1.2, 1.7, 1.0, 1.3, 1.3,
        2.3, 2.0, 3.1, 2.8,
        2.1, 2.0, 2.4,
        1.8, 1.6, 1.7, 4.5
    ],
    # Mean annual temperature (°C) - derived from ERA5
    'mean_temp_annual': [
        25.8, 25.2, 27.1, 27.5, 26.8, 27.0, 26.5,
        27.2, 27.8, 26.9, 26.5, 25.8, 25.4, 25.2, 25.8, 24.5,
        21.5, 23.2, 23.8, 20.8,
        19.5, 18.2, 18.5,
        24.2, 25.5, 23.8, 21.5
    ],
    # Climate zone (simplified)
    'climate_zone': [
        'Tropical', 'Tropical', 'Equatorial', 'Equatorial', 'Equatorial', 'Equatorial', 'Tropical',
        'Tropical', 'Semi-arid', 'Semi-arid', 'Semi-arid', 'Semi-arid', 
        'Tropical', 'Tropical', 'Tropical', 'Tropical',
        'Tropical', 'Tropical', 'Tropical', 'Subtropical',
        'Subtropical', 'Subtropical', 'Subtropical',
        'Tropical', 'Tropical', 'Tropical', 'Tropical'
    ],
    # Latitude of capital (for climate gradient)
    'capital_latitude': [
        -8.76, -9.97, -3.10, 2.82, -1.46, 0.03, -10.18,
        -2.53, -5.09, -3.72, -5.79, -7.12, -8.05, -9.67, -10.91, -12.97,
        -19.92, -20.32, -22.91, -23.55,
        -25.43, -27.60, -30.03,
        -20.44, -15.60, -16.68, -15.78
    ]
}

# Create dataframe
df = pd.DataFrame(states_data)

# =============================================================================
# DERIVED VARIABLES
# =============================================================================

# GDP per capita in USD (approximate rate 5 BRL/USD)
df['gdp_per_capita_usd'] = df['gdp_per_capita_brl'] / 5

# Population elderly (absolute)
df['elderly_population_k'] = df['population_2022_k'] * df['elderly_pct'] / 100

# Climate vulnerability index (higher = more adapted to heat)
# Based on baseline temperature (populations in hot climates more adapted)
df['heat_adaptation_proxy'] = (df['mean_temp_annual'] - df['mean_temp_annual'].min()) / \
                               (df['mean_temp_annual'].max() - df['mean_temp_annual'].min())

# Socioeconomic vulnerability index (composite)
# Higher = more vulnerable (lower GDP, lower HDI, less healthcare)
df['ses_vulnerability'] = (
    (1 - (df['gdp_per_capita_usd'] - df['gdp_per_capita_usd'].min()) / 
         (df['gdp_per_capita_usd'].max() - df['gdp_per_capita_usd'].min())) +
    (1 - (df['hdi'] - df['hdi'].min()) / (df['hdi'].max() - df['hdi'].min())) +
    (1 - (df['physicians_per_1000'] - df['physicians_per_1000'].min()) / 
         (df['physicians_per_1000'].max() - df['physicians_per_1000'].min()))
) / 3

# =============================================================================
# SUMMARY
# =============================================================================

print("\nState Covariates Summary:")
print("-"*50)
print(f"States: {len(df)}")
print(f"\nNumeric variables:")
for col in df.select_dtypes(include=[np.number]).columns:
    print(f"  {col}: {df[col].min():.2f} - {df[col].max():.2f} (mean: {df[col].mean():.2f})")

print(f"\nBy region:")
region_summary = df.groupby('region').agg({
    'population_2022_k': 'sum',
    'elderly_pct': 'mean',
    'urban_pct': 'mean',
    'gdp_per_capita_usd': 'mean',
    'mean_temp_annual': 'mean'
}).round(1)
print(region_summary)

# =============================================================================
# SAVE
# =============================================================================

df.to_csv('./results/state_covariates.csv', index=False)
print(f"\nSaved to: ./results/state_covariates.csv")

df.to_parquet('./results/state_covariates.parquet', index=False)
print(f"Saved to: ./results/state_covariates.parquet")

# Save data dictionary
data_dict = """
# State Covariates Data Dictionary

## Source
- Population, urbanization, elderly %: IBGE Censo 2022
- GDP: IBGE Contas Regionais 2021  
- HDI: Atlas Brasil 2021
- Hospital beds: DATASUS/CNES 2022
- Physicians: CFM 2022
- Temperature: ERA5 reanalysis (study period mean)

## Variables

| Variable | Description | Unit |
|----------|-------------|------|
| state | State abbreviation | - |
| state_name | Full state name | - |
| region | Geographic region | - |
| population_2022_k | Total population 2022 | thousands |
| elderly_pct | Population aged 60+ | % |
| urban_pct | Urban population | % |
| gdp_per_capita_brl | GDP per capita | BRL |
| gdp_per_capita_usd | GDP per capita | USD |
| hdi | Human Development Index | 0-1 |
| hospital_beds_per_1000 | Hospital beds | per 1000 pop |
| physicians_per_1000 | Physicians | per 1000 pop |
| mean_temp_annual | Mean annual temperature | °C |
| climate_zone | Simplified climate zone | - |
| capital_latitude | Latitude of state capital | degrees |
| elderly_population_k | Elderly population | thousands |
| heat_adaptation_proxy | Proxy for heat adaptation (0-1) | index |
| ses_vulnerability | Socioeconomic vulnerability (0-1) | index |

## Notes
- Data compiled for meta-regression of state-level heat effects
- GDP converted at approximate rate of 5 BRL/USD
- heat_adaptation_proxy assumes populations in hotter climates are more adapted
- ses_vulnerability is a composite index (GDP, HDI, physicians)
"""

with open('./results/state_covariates_dictionary.md', 'w') as f:
    f.write(data_dict)

print(f"\n{'='*70}")
print("DONE!")
print("="*70)
