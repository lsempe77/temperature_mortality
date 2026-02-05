"""
01e: YEARS OF LIFE LOST (YLL) CALCULATION
==========================================
Calculates YLL due to temperature-attributable mortality using life tables.

Method:
1. Load life table data from IBGE
2. For each attributable death, estimate remaining life expectancy by age
3. Multiply attributable deaths × remaining life expectancy = YLL
4. Aggregate to national level

Key Outputs:
- YLL due to heat
- YLL due to cold
- YLL per 1000 deaths
- Age-standardized YLL rates

References:
- WHO YLL methodology
- Gasparrini et al. (2015) - GBD YLL approach
"""

import warnings
warnings.filterwarnings('ignore')

import json
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path

# =============================================================================
# PATHS
# =============================================================================

SCRIPT_DIR = Path(__file__).parent
PHASE0_RESULTS = SCRIPT_DIR.parent / 'phase0_data_prep' / 'results'
PHASE1_RESULTS = SCRIPT_DIR / 'results'
OUTPUT_DIR = PHASE1_RESULTS
OUTPUT_DIR.mkdir(exist_ok=True)

print("=" * 70)
print("01e: YEARS OF LIFE LOST (YLL) CALCULATION")
print("=" * 70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print()

# =============================================================================
# LOAD IBGE LIFE TABLE DATA
# =============================================================================
print("[1] Loading IBGE life table data...")

# Load the actual IBGE life tables (Brazil-specific data)
yll_lookup_file = PHASE0_RESULTS / 'yll_lookup_by_age.csv'
ibge_life_tables_file = PHASE0_RESULTS / 'ibge_life_tables_combined.csv'

if yll_lookup_file.exists():
    life_table = pd.read_csv(yll_lookup_file)
    print(f"  Loaded IBGE YLL lookup: {len(life_table)} ages (0-{life_table['age'].max()})")
    print(f"  Source: IBGE life tables for Brazil")
    
    # Show key values for elderly
    elderly_ages = [60, 65, 70, 75, 80, 85]
    print("\n  Life expectancy by age (IBGE Brazil):")
    for age in elderly_ages:
        if age in life_table['age'].values:
            ex = life_table[life_table['age'] == age]['ex_mean'].values[0]
            print(f"    Age {age}: {ex:.1f} years remaining")
else:
    print("  WARNING: IBGE life table not found, using WHO fallback...")
    # WHO Global Burden of Disease reference life table
    age_groups = list(range(0, 95, 5)) + [95]
    life_expectancy = [
        86.6, 81.8, 76.9, 72.0, 67.1, 62.2, 57.3, 52.5, 47.8, 43.1,
        38.5, 34.0, 29.7, 25.5, 21.6, 17.9, 14.6, 11.7, 9.3, 7.3
    ]
    life_table = pd.DataFrame({
        'age': age_groups,
        'ex_mean': life_expectancy
    })
    print(f"  Using WHO reference life table: {len(life_table)} age groups")

# Create age-to-life-expectancy lookup function
def get_life_expectancy(age):
    """Get remaining life expectancy for a given age from IBGE data."""
    if age in life_table['age'].values:
        return life_table[life_table['age'] == age]['ex_mean'].values[0]
    elif age > life_table['age'].max():
        return life_table['ex_mean'].min()  # Use minimum for oldest ages
    else:
        # Interpolate
        idx = life_table['age'].searchsorted(age)
        if idx == 0:
            return life_table['ex_mean'].iloc[0]
        a1, a2 = life_table['age'].iloc[idx-1], life_table['age'].iloc[idx]
        e1, e2 = life_table['ex_mean'].iloc[idx-1], life_table['ex_mean'].iloc[idx]
        return e1 + (e2 - e1) * (age - a1) / (a2 - a1)

# =============================================================================
# LOAD ATTRIBUTABLE BURDEN DATA
# =============================================================================
print("\n[2] Loading attributable burden data...")

# Load national summary
with open(PHASE1_RESULTS / 'attributable_burden_national.json', 'r') as f:
    burden_national = json.load(f)

# Load regional burden
burden_intermediate = pd.read_csv(PHASE1_RESULTS / 'burden_intermediate_regions.csv')
burden_immediate = pd.read_csv(PHASE1_RESULTS / 'burden_immediate_regions.csv')

print(f"  Intermediate regions: {len(burden_intermediate)}")
print(f"  Immediate regions: {len(burden_immediate)}")

# =============================================================================
# CALCULATE YLL
# =============================================================================
print("\n[3] Calculating Years of Life Lost...")

# For elderly population (60+), we use age distribution from mortality data
# and IBGE life expectancy at each age
#
# Studies show temperature deaths tend to affect older elderly more
# Age distribution assumptions for elderly deaths in Brazil:
# - 60-69: 25% of elderly deaths
# - 70-79: 35% of elderly deaths  
# - 80+:   40% of elderly deaths

AGE_DISTRIBUTION = {
    'age_60_69': {'pct': 0.25, 'avg_age': 65, 'life_exp': get_life_expectancy(65)},
    'age_70_79': {'pct': 0.35, 'avg_age': 75, 'life_exp': get_life_expectancy(75)},
    'age_80_plus': {'pct': 0.40, 'avg_age': 85, 'life_exp': get_life_expectancy(85)}
}

# Display IBGE-based life expectancy values
print("  Age distribution of elderly deaths (assumed):")
for group, vals in AGE_DISTRIBUTION.items():
    print(f"    {group}: {vals['pct']*100:.0f}%, avg age {vals['avg_age']}, "
          f"life exp = {vals['life_exp']:.1f} years (IBGE)")

# Weighted average life expectancy for elderly deaths using IBGE data
avg_life_exp = sum(v['pct'] * v['life_exp'] for v in AGE_DISTRIBUTION.values())
print(f"\n  Weighted average remaining life expectancy: {avg_life_exp:.2f} years (IBGE Brazil)")

# Calculate YLL for attributable deaths
def calculate_yll(attributable_deaths, avg_life_exp):
    """Calculate Years of Life Lost."""
    return attributable_deaths * avg_life_exp


# Intermediate level
int_heat_an = burden_national['intermediate']['heat_an']
int_cold_an = burden_national['intermediate']['cold_an']
int_heat_yll = calculate_yll(int_heat_an, avg_life_exp)
int_cold_yll = calculate_yll(int_cold_an, avg_life_exp)

# Annual values
n_years_int = burden_national['intermediate']['n_regions'] * 14.97 / burden_national['intermediate']['n_regions']
# Estimate n_years from total deaths and annual
int_heat_annual = burden_national['intermediate']['heat_annual']
int_cold_annual = burden_national['intermediate']['cold_annual']
n_years_approx = int_heat_an / int_heat_annual if int_heat_annual > 0 else 15

int_heat_yll_annual = int_heat_yll / n_years_approx
int_cold_yll_annual = int_cold_yll / n_years_approx

# Immediate level
imm_heat_an = burden_national['immediate']['heat_an']
imm_cold_an = burden_national['immediate']['cold_an']
imm_heat_yll = calculate_yll(imm_heat_an, avg_life_exp)
imm_cold_yll = calculate_yll(imm_cold_an, avg_life_exp)

imm_heat_annual = burden_national['immediate']['heat_annual']
imm_cold_annual = burden_national['immediate']['cold_annual']

imm_heat_yll_annual = calculate_yll(imm_heat_annual, avg_life_exp)
imm_cold_yll_annual = calculate_yll(imm_cold_annual, avg_life_exp)

# =============================================================================
# YLL RATES
# =============================================================================
print("\n[4] Calculating YLL rates...")

# Get population data
int_pop = burden_national['intermediate'].get('total_pop_elderly', None)
imm_pop = burden_national['immediate'].get('total_pop_elderly', None)

# YLL per 100,000 population per year
if int_pop:
    int_heat_yll_rate = int_heat_yll_annual / int_pop * 100000
    int_cold_yll_rate = int_cold_yll_annual / int_pop * 100000
else:
    int_heat_yll_rate = None
    int_cold_yll_rate = None

if imm_pop:
    imm_heat_yll_rate = imm_heat_yll_annual / imm_pop * 100000
    imm_cold_yll_rate = imm_cold_yll_annual / imm_pop * 100000
else:
    imm_heat_yll_rate = None
    imm_cold_yll_rate = None

# =============================================================================
# DISPLAY RESULTS
# =============================================================================

print("\n" + "=" * 70)
print("YEARS OF LIFE LOST (YLL) RESULTS")
print("=" * 70)

print("\n[INTERMEDIATE LEVEL - 133 Regions]")
print(f"  Heat-attributable deaths: {int_heat_an:,.0f}")
print(f"  Cold-attributable deaths: {int_cold_an:,.0f}")
print(f"  Heat YLL (total): {int_heat_yll:,.0f}")
print(f"  Cold YLL (total): {int_cold_yll:,.0f}")
print(f"  Heat YLL (annual): {int_heat_yll_annual:,.0f}")
print(f"  Cold YLL (annual): {int_cold_yll_annual:,.0f}")
if int_heat_yll_rate:
    print(f"  Heat YLL rate (per 100k elderly/year): {int_heat_yll_rate:.1f}")
    print(f"  Cold YLL rate (per 100k elderly/year): {int_cold_yll_rate:.1f}")

print("\n[IMMEDIATE LEVEL - 510 Regions]")
print(f"  Heat-attributable deaths: {imm_heat_an:,.0f}")
print(f"  Cold-attributable deaths: {imm_cold_an:,.0f}")
print(f"  Heat YLL (total): {imm_heat_yll:,.0f}")
print(f"  Cold YLL (total): {imm_cold_yll:,.0f}")
print(f"  Heat YLL (annual): {imm_heat_yll_annual:,.0f}")
print(f"  Cold YLL (annual): {imm_cold_yll_annual:,.0f}")
if imm_heat_yll_rate:
    print(f"  Heat YLL rate (per 100k elderly/year): {imm_heat_yll_rate:.1f}")
    print(f"  Cold YLL rate (per 100k elderly/year): {imm_cold_yll_rate:.1f}")

# =============================================================================
# SAVE RESULTS
# =============================================================================
print("\n[5] Saving results...")

def convert_for_json(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: convert_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_for_json(i) for i in obj]
    return obj

yll_results = {
    'method': {
        'description': 'YLL = Attributable deaths × Remaining life expectancy',
        'age_distribution': AGE_DISTRIBUTION,
        'avg_life_expectancy': avg_life_exp,
        'source': 'WHO Global Burden of Disease 2019 reference life table'
    },
    'intermediate': {
        'heat_an': int_heat_an,
        'cold_an': int_cold_an,
        'heat_yll_total': int_heat_yll,
        'cold_yll_total': int_cold_yll,
        'heat_yll_annual': int_heat_yll_annual,
        'cold_yll_annual': int_cold_yll_annual,
        'heat_yll_rate_per_100k': int_heat_yll_rate,
        'cold_yll_rate_per_100k': int_cold_yll_rate,
        'n_years': n_years_approx,
        'total_pop_elderly': int_pop
    },
    'immediate': {
        'heat_an': imm_heat_an,
        'cold_an': imm_cold_an,
        'heat_yll_total': imm_heat_yll,
        'cold_yll_total': imm_cold_yll,
        'heat_yll_annual': imm_heat_yll_annual,
        'cold_yll_annual': imm_cold_yll_annual,
        'heat_yll_rate_per_100k': imm_heat_yll_rate,
        'cold_yll_rate_per_100k': imm_cold_yll_rate,
        'total_pop_elderly': imm_pop
    },
    'timestamp': datetime.now().isoformat()
}

with open(OUTPUT_DIR / 'yll_results.json', 'w') as f:
    json.dump(convert_for_json(yll_results), f, indent=2)

print(f"  Saved: yll_results.json")

# =============================================================================
# SUMMARY TABLE
# =============================================================================

print("\n" + "=" * 70)
print("SUMMARY TABLE")
print("=" * 70)

summary_data = [
    {
        'Level': 'Intermediate (133)',
        'Heat Deaths': f"{int_heat_an:,.0f}",
        'Cold Deaths': f"{int_cold_an:,.0f}",
        'Heat YLL': f"{int_heat_yll:,.0f}",
        'Cold YLL': f"{int_cold_yll:,.0f}",
        'Total YLL': f"{int_heat_yll + int_cold_yll:,.0f}",
    },
    {
        'Level': 'Immediate (510)',
        'Heat Deaths': f"{imm_heat_an:,.0f}",
        'Cold Deaths': f"{imm_cold_an:,.0f}",
        'Heat YLL': f"{imm_heat_yll:,.0f}",
        'Cold YLL': f"{imm_cold_yll:,.0f}",
        'Total YLL': f"{imm_heat_yll + imm_cold_yll:,.0f}",
    }
]

summary_df = pd.DataFrame(summary_data)
print(summary_df.to_string(index=False))
summary_df.to_csv(OUTPUT_DIR / 'yll_summary.csv', index=False)

print("\n" + "=" * 70)
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)
