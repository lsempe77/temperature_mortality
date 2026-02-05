"""
00k2: Download SES data from IBGE SIDRA API - FIXED VERSION
============================================================
Uses direct API calls instead of sidrapy for better control over classifications.

Variables:
1. Total Population (Census 2022)
2. Elderly Population 60+ (Census 2022) 
3. GDP (2021)
4. Urban Population (Census 2022)

Output:
- results/municipality_ses_covariates.csv
"""

import requests
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import time

print("="*70)
print("00k2: DOWNLOAD SES DATA (FIXED)")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# PATHS
# =============================================================================

BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR.parent / 'results'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def fetch_sidra_direct(url, max_retries=3):
    """Fetch data from SIDRA API with retries."""
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=120)
            if response.status_code == 200:
                return response.json()
            else:
                print(f"  API returned status {response.status_code}")
                if attempt < max_retries - 1:
                    time.sleep(5)
        except Exception as e:
            print(f"  Error: {e}")
            if attempt < max_retries - 1:
                time.sleep(5)
    return None


def fetch_all_municipalities(base_url, batch_size=500):
    """
    Fetch data for all municipalities in batches.
    SIDRA has limits on how many results can be returned at once.
    """
    # First, get total count with a small request
    test_url = base_url.replace("/n6/all/", "/n6/3550308/")  # Just São Paulo for test
    test_data = fetch_sidra_direct(test_url)
    if not test_data:
        return None
    
    # Now fetch all municipalities
    # The API supports "all" for municipalities
    full_url = base_url
    print(f"  Fetching all municipalities...")
    data = fetch_sidra_direct(full_url)
    
    if data:
        print(f"  Got {len(data)} rows")
        return data
    return None


# =============================================================================
# 1. TOTAL POPULATION (Census 2022)
# =============================================================================

print("\n" + "-"*70)
print("[1/4] Total Population (Census 2022)")
print("-"*70)

# Table 9514, Variable 93, Classification 2 (Sex)=6794 (Total), Classification 287 (Age)=100362 (Total)
# URL format: /t/TABLE/n6/all/v/VAR/p/PERIOD/c2/SEX/c287/AGE/f/c (f/c = codes only)
url_pop_total = "https://apisidra.ibge.gov.br/values/t/9514/n6/all/v/93/p/2022/c2/6794/c287/100362/f/c"

pop_data = fetch_sidra_direct(url_pop_total)

if pop_data:
    # Convert to DataFrame (first row is header with column descriptions)
    df_pop = pd.DataFrame(pop_data[1:])  # Skip header row
    # With /f/c format: D1C = municipality code, V = value
    df_pop = df_pop[['D1C', 'V']].copy()
    df_pop.columns = ['code_muni', 'pop_total']
    df_pop['code_muni'] = pd.to_numeric(df_pop['code_muni'], errors='coerce')
    df_pop['pop_total'] = pd.to_numeric(df_pop['pop_total'], errors='coerce')
    df_pop = df_pop.dropna()
    print(f"  Total population: {len(df_pop)} municipalities")
    print(f"  Sample: {df_pop.head(3).to_string()}")
else:
    print("  ERROR: Failed to fetch total population")
    df_pop = None

# =============================================================================
# 2. ELDERLY POPULATION 60+ (Census 2022)
# =============================================================================

print("\n" + "-"*70)
print("[2/4] Elderly Population 60+ (Census 2022)")
print("-"*70)

# Age group codes for 60+ (5-year intervals) from Census 2022:
# Split into 3 batches to stay under SIDRA's 50,000 value limit
# (5570 municipalities × 3 age groups = 16,710 per batch < 50,000)
age_batches = [
    "93095,93096,93097",  # 60-64, 65-69, 70-74
    "93098,49108,49109",  # 75-79, 80-84, 85-89
    "60040,60041,6653"    # 90-94, 95-99, 100+
]

print(f"  Fetching elderly pop (60+) in {len(age_batches)} batches...")
all_elderly_rows = []
for i, codes in enumerate(age_batches):
    url = f"https://apisidra.ibge.gov.br/values/t/9514/n6/all/v/93/p/2022/c2/6794/c287/{codes}/f/c"
    data = fetch_sidra_direct(url)
    if data:
        print(f"    Batch {i+1}: {len(data)-1} rows")
        all_elderly_rows.extend(data[1:])
    else:
        print(f"    Batch {i+1}: FAILED")
    time.sleep(1)  # Small delay between requests

if all_elderly_rows:
    df_elderly = pd.DataFrame(all_elderly_rows)
    df_elderly['code_muni'] = pd.to_numeric(df_elderly['D1C'], errors='coerce')
    df_elderly['value'] = pd.to_numeric(df_elderly['V'], errors='coerce')
    
    # Sum by municipality (across all 60+ age groups)
    df_elderly_sum = df_elderly.groupby('code_muni')['value'].sum().reset_index()
    df_elderly_sum.columns = ['code_muni', 'pop_elderly']
    print(f"  Elderly population: {len(df_elderly_sum)} municipalities")
else:
    print("  ERROR: Failed to fetch elderly population")
    df_elderly_sum = None

# =============================================================================
# 3. GDP (2021)
# =============================================================================

print("\n" + "-"*70)
print("[3/4] GDP (2021)")
print("-"*70)

# Table 5938, Variable 37 (GDP at current prices)
url_gdp = "https://apisidra.ibge.gov.br/values/t/5938/n6/all/v/37/p/2021/f/c"

gdp_data = fetch_sidra_direct(url_gdp)

if gdp_data:
    df_gdp = pd.DataFrame(gdp_data[1:])
    # With /f/c format: D1C = municipality code, V = value
    df_gdp = df_gdp[['D1C', 'V']].copy()
    df_gdp.columns = ['code_muni', 'gdp_total']
    df_gdp['code_muni'] = pd.to_numeric(df_gdp['code_muni'], errors='coerce')
    df_gdp['gdp_total'] = pd.to_numeric(df_gdp['gdp_total'], errors='coerce')
    df_gdp = df_gdp.dropna()
    print(f"  GDP data: {len(df_gdp)} municipalities")
else:
    print("  ERROR: Failed to fetch GDP")
    df_gdp = None

# =============================================================================
# 4. URBAN POPULATION (Census 2010)
# =============================================================================

print("\n" + "-"*70)
print("[4/4] Urban Population (Census 2010)")
print("-"*70)

# Table 1378: Population by situation (urban/rural) - Census 2010
# Classification 1 = 1 (Urban)
# Note: 2022 Census doesn't have this classification published yet for municipalities
url_urban = "https://apisidra.ibge.gov.br/values/t/1378/n6/all/v/93/p/2010/c1/1/f/c"

print("  Fetching urban population (Census 2010)...")
urban_data = fetch_sidra_direct(url_urban)

if urban_data:
    df_urban = pd.DataFrame(urban_data[1:])
    # With /f/c format: D1C = municipality code, V = value
    df_urban = df_urban[['D1C', 'V']].copy()
    df_urban.columns = ['code_muni', 'pop_urban']
    df_urban['code_muni'] = pd.to_numeric(df_urban['code_muni'], errors='coerce')
    df_urban['pop_urban'] = pd.to_numeric(df_urban['pop_urban'], errors='coerce')
    df_urban = df_urban.dropna()
    print(f"  Urban population: {len(df_urban)} municipalities")
else:
    print("  ERROR: Failed to fetch urban population")
    df_urban = None

# =============================================================================
# MERGE AND SAVE
# =============================================================================

print("\n" + "-"*70)
print("Merging and Saving")
print("-"*70)

if df_pop is not None:
    df_ses = df_pop.copy()
    
    # Merge elderly
    if df_elderly_sum is not None:
        df_ses = pd.merge(df_ses, df_elderly_sum, on='code_muni', how='left')
        print(f"  Merged elderly: {df_ses['pop_elderly'].notna().sum()} matched")
    
    # Merge GDP
    if df_gdp is not None:
        df_ses = pd.merge(df_ses, df_gdp, on='code_muni', how='left')
        print(f"  Merged GDP: {df_ses['gdp_total'].notna().sum()} matched")
    
    # Merge urban
    if df_urban is not None:
        df_ses = pd.merge(df_ses, df_urban, on='code_muni', how='left')
        print(f"  Merged urban: {df_ses['pop_urban'].notna().sum()} matched")
    
    # Calculate derived variables
    if 'pop_elderly' in df_ses.columns:
        df_ses['pct_elderly'] = (df_ses['pop_elderly'] / df_ses['pop_total']) * 100
    
    if 'gdp_total' in df_ses.columns:
        # GDP is in R$ thousands
        df_ses['gdp_per_capita'] = (df_ses['gdp_total'] * 1000) / df_ses['pop_total']
    
    if 'pop_urban' in df_ses.columns:
        df_ses['urbanization_rate'] = (df_ses['pop_urban'] / df_ses['pop_total']) * 100
        # Cap at 100% (some municipalities may have grown)
        df_ses['urbanization_rate'] = df_ses['urbanization_rate'].clip(upper=100)
    
    # Save
    output_file = OUTPUT_DIR / "municipality_ses_covariates.csv"
    df_ses.to_csv(output_file, index=False)
    
    print(f"\n{'='*70}")
    print(f"SAVED: {output_file}")
    print(f"{'='*70}")
    print(f"\nDataset summary:")
    print(f"  Municipalities: {len(df_ses)}")
    print(f"\nVariable statistics:")
    for col in ['pop_total', 'pop_elderly', 'pct_elderly', 'gdp_per_capita', 'urbanization_rate']:
        if col in df_ses.columns:
            print(f"  {col}: mean={df_ses[col].mean():.2f}, median={df_ses[col].median():.2f}")
    
    print(f"\nSample rows:")
    print(df_ses.head(5).to_string())
else:
    print("ERROR: Could not create SES dataset - no population data")

print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
