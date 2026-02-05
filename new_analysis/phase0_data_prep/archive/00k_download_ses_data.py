"""
Download SES data from IBGE SIDRA API for all municipalities.
We will later aggregate this to Intermediate Regions.

Variables to download:
1. Population by Age Group (Census 2022) -> To calculate Elderly %
2. Total Population (Census 2022)
3. GDP (2021) -> Economic indicator
4. Urban/Rural Population (Census 2010/2022) -> Urbanization rate

Library: sidrapy
"""

import sidrapy
import pandas as pd
import numpy as np
from pathlib import Path
import time

def get_sidra_data(table_code, variable_code, classification=None, period="last"):
    """
    Helper to fetch data from SIDRA.
    """
    print(f"Fetching Table {table_code}, Var {variable_code}...")
    
    try:
        # sidrapy.get_table returns a list of dicts or a dataframe
        # We want data for all municipalities (territorial_level="6")
        
        args = {
            "table_code": table_code,
            "territorial_level": "6", # Municipality
            "ibge_territorial_code": "all",
            "variable": variable_code,
            "period": period
        }
        
        if classification:
            for k, v in classification.items():
                args[k] = v
                
        df = sidrapy.get_table(**args)
        
        if df is None or df.empty:
            print(f"Warning: No data returned for Table {table_code}")
            return None
            
        # Clean up header (first row is usually descriptions)
        df.columns = df.iloc[0]
        df = df.iloc[1:]
        
        return df
        
    except Exception as e:
        print(f"Error fetching SIDRA data: {e}")
        return None

def download_ses_data():
    output_dir = Path(__file__).resolve().parents[1] / "results"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # ---------------------------------------------------------
    # 1. POPULATION & ELDERLY (Census 2022)
    # Table 9514: Population by Age Group
    # Variable 93: Resident population (People) - CORRECTED from 1000093
    # Classification 2: Sex (Total = 6794)
    # Classification 287: Age Group 
    #   Total = 652
    #   60-64 = 3244, 65-69 = 3245, ... 100+ = 3252
    #   We need to sum 60+
    # ---------------------------------------------------------
    
    print("\n[1/3] Downloading Population Data (Census 2022)...")
    
    # Fetch Total Population
    pop_total = get_sidra_data(
        table_code="9514",
        variable_code="93",
        classification={"classificacao": "2/6794|287/652"}, # Total Sex, Total Age
        period="2022"
    )
    
    # Fetch Elderly Population (60+)
    # We might need to fetch all age groups and sum, or fetch specific codes.
    # Let's fetch specific age codes for 60+
    # Codes: 3244,3245,3246,3247,3248,3249,3250,3251,3252 (60 to 100+)
    elderly_codes = "3244,3245,3246,3247,3248,3249,3250,3251,3252"
    
    pop_elderly_raw = get_sidra_data(
        table_code="9514",
        variable_code="93",
        classification={"classificacao": f"2/6794|287/{elderly_codes}"},
        period="2022"
    )
    
    # Process Population
    if pop_total is not None and pop_elderly_raw is not None:
        # Prepare Total
        df_pop = pop_total[['Município (Código)', 'Valor']].copy()
        df_pop.columns = ['code_muni', 'pop_total']
        
        # Prepare Elderly (Sum by Muni)
        pop_elderly_raw['Valor'] = pd.to_numeric(pop_elderly_raw['Valor'], errors='coerce')
        df_elderly = pop_elderly_raw.groupby('Município (Código)')['Valor'].sum().reset_index()
        df_elderly.columns = ['code_muni', 'pop_elderly']
        
        # Merge
        df_ses = pd.merge(df_pop, df_elderly, on='code_muni', how='left')
        print(f"Population data processed. Rows: {len(df_ses)}")
    else:
        print("Failed to get population data.")
        return

    # ---------------------------------------------------------
    # 2. GDP (2021 - latest available usually)
    # Table 5938: GDP at current prices
    # Variable 37: GDP
    # ---------------------------------------------------------
    print("\n[2/3] Downloading GDP Data (2021)...")
    
    gdp_raw = get_sidra_data(
        table_code="5938",
        variable_code="37",
        period="2021"
    )
    
    if gdp_raw is not None:
        df_gdp = gdp_raw[['Município (Código)', 'Valor']].copy()
        df_gdp.columns = ['code_muni', 'gdp_total']
        
        # Merge
        df_ses = pd.merge(df_ses, df_gdp, on='code_muni', how='left')
        print(f"GDP data merged.")
        
    # ---------------------------------------------------------
    # 3. URBANIZATION (Census 2010 - 2022 might be partial)
    # Table 9514 also has Situation (Urban/Rural) -> Classificacao 1
    # Urban = 1, Rural = 2
    # ---------------------------------------------------------
    print("\n[3/3] Downloading Urbanization Data (Census 2022)...")
    
    # Fetch Urban Population
    pop_urban = get_sidra_data(
        table_code="9514",
        variable_code="93",
        classification={"classificacao": "2/6794|287/652|1/1"}, # Total Sex, Total Age, Urban
        period="2022"
    )
    
    if pop_urban is not None:
        df_urban = pop_urban[['Município (Código)', 'Valor']].copy()
        df_urban.columns = ['code_muni', 'pop_urban']
        
        # Merge
        df_ses = pd.merge(df_ses, df_urban, on='code_muni', how='left')
        print(f"Urbanization data merged.")

    # ---------------------------------------------------------
    # CLEANUP AND SAVE
    # ---------------------------------------------------------
    
    # Convert columns to numeric
    cols = ['pop_total', 'pop_elderly', 'gdp_total', 'pop_urban']
    for col in cols:
        if col in df_ses.columns:
            df_ses[col] = pd.to_numeric(df_ses[col], errors='coerce')
            
    # Calculate derived variables
    df_ses['pct_elderly'] = (df_ses['pop_elderly'] / df_ses['pop_total']) * 100
    df_ses['gdp_per_capita'] = df_ses['gdp_total'] / df_ses['pop_total'] # Note: GDP is usually x1000, check units
    df_ses['urbanization_rate'] = (df_ses['pop_urban'] / df_ses['pop_total']) * 100
    
    # Save
    output_file = output_dir / "municipality_ses_covariates.csv"
    df_ses.to_csv(output_file, index=False)
    print(f"\nSAVED: {output_file}")
    print(df_ses.head())

if __name__ == "__main__":
    download_ses_data()
