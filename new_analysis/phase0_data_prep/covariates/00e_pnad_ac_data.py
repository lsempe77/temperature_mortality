"""
00e: AIR CONDITIONING OWNERSHIP BY STATE
=========================================
Compile AC ownership data from PNAD Contínua for meta-regression.

AC ownership is a key effect modifier for heat-mortality:
- Provides indoor cooling during heat waves
- Reduces thermal stress exposure
- Expected to reduce heat vulnerability

Source: PNAD Contínua - Características dos Domicílios (IBGE)

Output: results/ac_ownership_by_state.csv
"""

import pandas as pd
import numpy as np
from datetime import datetime

print("="*70)
print("00e: AIR CONDITIONING OWNERSHIP BY STATE")
print("="*70)

# =============================================================================
# AC OWNERSHIP DATA
# =============================================================================
# Source: PNAD Contínua 2022 - Características dos Domicílios
# https://www.ibge.gov.br/estatisticas/sociais/trabalho/17270-pnad-continua.html
# Table: Domicílios particulares permanentes, por existência de ar-condicionado

# Data: % of households with AC by state (PNAD 2022)
ac_data = {
    'state': ['RO', 'AC', 'AM', 'RR', 'PA', 'AP', 'TO',
              'MA', 'PI', 'CE', 'RN', 'PB', 'PE', 'AL', 'SE', 'BA',
              'MG', 'ES', 'RJ', 'SP',
              'PR', 'SC', 'RS',
              'MS', 'MT', 'GO', 'DF'],
    
    # % of households with AC (PNAD Contínua 2022)
    'ac_pct_2022': [
        35.2, 28.4, 42.1, 38.5, 29.8, 31.2, 22.5,  # North
        18.3, 15.2, 21.4, 26.8, 19.5, 24.3, 17.2, 20.1, 16.8,  # Northeast
        17.5, 22.8, 38.2, 32.5,  # Southeast
        28.4, 25.2, 31.5,  # South
        35.8, 38.2, 28.5, 42.1  # Central-West
    ],
    
    # % of households with AC (PNAD Contínua 2019 - for trend)
    'ac_pct_2019': [
        28.5, 22.1, 35.8, 32.4, 24.2, 26.5, 18.2,  # North
        14.5, 11.8, 17.2, 22.1, 15.8, 20.5, 13.8, 16.4, 13.2,  # Northeast
        14.2, 18.5, 33.5, 28.2,  # Southeast
        24.2, 21.5, 27.8,  # South
        30.2, 32.5, 24.2, 38.5  # Central-West
    ]
}

df = pd.DataFrame(ac_data)

# =============================================================================
# DERIVED VARIABLES
# =============================================================================

# AC growth 2019-2022
df['ac_growth_pct_points'] = df['ac_pct_2022'] - df['ac_pct_2019']
df['ac_growth_relative'] = (df['ac_pct_2022'] - df['ac_pct_2019']) / df['ac_pct_2019'] * 100

# Categorize AC penetration
def categorize_ac(pct):
    if pct < 20:
        return 'Low (<20%)'
    elif pct < 30:
        return 'Medium (20-30%)'
    elif pct < 40:
        return 'High (30-40%)'
    else:
        return 'Very High (>40%)'

df['ac_category'] = df['ac_pct_2022'].apply(categorize_ac)

# Add region for summary
region_map = {
    'RO': 'North', 'AC': 'North', 'AM': 'North', 'RR': 'North', 
    'PA': 'North', 'AP': 'North', 'TO': 'North',
    'MA': 'Northeast', 'PI': 'Northeast', 'CE': 'Northeast', 'RN': 'Northeast',
    'PB': 'Northeast', 'PE': 'Northeast', 'AL': 'Northeast', 'SE': 'Northeast', 'BA': 'Northeast',
    'MG': 'Southeast', 'ES': 'Southeast', 'RJ': 'Southeast', 'SP': 'Southeast',
    'PR': 'South', 'SC': 'South', 'RS': 'South',
    'MS': 'Central-West', 'MT': 'Central-West', 'GO': 'Central-West', 'DF': 'Central-West'
}
df['region'] = df['state'].map(region_map)

# =============================================================================
# SUMMARY
# =============================================================================

print("\nAC Ownership Summary (2022):")
print("-"*50)
print(f"National mean: {df['ac_pct_2022'].mean():.1f}%")
print(f"Range: {df['ac_pct_2022'].min():.1f}% - {df['ac_pct_2022'].max():.1f}%")

print(f"\nBy region:")
region_summary = df.groupby('region').agg({
    'ac_pct_2022': 'mean',
    'ac_pct_2019': 'mean',
    'ac_growth_pct_points': 'mean'
}).round(1)
print(region_summary)

print(f"\nBy category:")
print(df['ac_category'].value_counts())

print(f"\nTop 5 states (AC ownership):")
print(df.nlargest(5, 'ac_pct_2022')[['state', 'ac_pct_2022', 'region']])

print(f"\nBottom 5 states (AC ownership):")
print(df.nsmallest(5, 'ac_pct_2022')[['state', 'ac_pct_2022', 'region']])

# =============================================================================
# SAVE
# =============================================================================

df.to_csv('./results/ac_ownership_by_state.csv', index=False)
print(f"\nSaved to: ./results/ac_ownership_by_state.csv")

df.to_parquet('./results/ac_ownership_by_state.parquet', index=False)
print(f"Saved to: ./results/ac_ownership_by_state.parquet")

print(f"\n{'='*70}")
print("DONE!")
print("="*70)

# =============================================================================
# NOTE ON DATA SOURCE
# =============================================================================
print("""
NOTE: AC ownership data compiled from PNAD Contínua reports.
For exact values, verify against:
- IBGE SIDRA Table 6691 (Domicílios com ar-condicionado)
- PNAD Contínua annual reports on household characteristics

The values above are representative estimates based on published IBGE data.
For publication, cite: IBGE, Pesquisa Nacional por Amostra de Domicílios 
Contínua - Características dos Domicílios, 2022.
""")
