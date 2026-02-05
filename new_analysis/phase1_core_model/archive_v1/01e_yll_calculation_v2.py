"""
01e_yll_calculation_v2.py
==========================
Years of Life Lost (YLL) Calculation for Attributable Burden v2

Key Changes from v1:
1. Compatible with burden_v2 outputs (using P2.5/P97.5 thresholds)
2. Better edge case handling (division by zero, NaN guards)
3. Clearer n_years calculation
4. Uses IBGE life tables (Brazil-specific)

Method:
1. Load IBGE life table data
2. For each attributable death, estimate remaining life expectancy by age
3. Multiply attributable deaths × remaining life expectancy = YLL
4. Aggregate to national level

References:
- WHO YLL methodology
- Gasparrini et al. (2015) - GBD YLL approach
- IBGE Brazilian life tables
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
print("01e: YEARS OF LIFE LOST (YLL) CALCULATION v2")
print("=" * 70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# LOAD IBGE LIFE TABLE DATA
# =============================================================================
print("\n[1] Loading IBGE life table data...")

yll_lookup_file = PHASE0_RESULTS / 'yll_lookup_by_age.csv'
ibge_life_tables_file = PHASE0_RESULTS / 'ibge_life_tables_combined.csv'

if yll_lookup_file.exists():
    life_table = pd.read_csv(yll_lookup_file)
    print(f"  Loaded IBGE YLL lookup: {len(life_table)} ages")
    life_table_source = 'IBGE Brazil'
    
    # Show key values for elderly
    elderly_ages = [60, 65, 70, 75, 80, 85]
    print("\n  Life expectancy by age (IBGE Brazil):")
    for age in elderly_ages:
        if age in life_table['age'].values:
            ex = life_table[life_table['age'] == age]['ex_mean'].values[0]
            print(f"    Age {age}: {ex:.1f} years remaining")
            
elif ibge_life_tables_file.exists():
    life_table = pd.read_csv(ibge_life_tables_file)
    print(f"  Loaded IBGE combined life tables")
    life_table_source = 'IBGE Brazil'
else:
    print("  WARNING: IBGE life table not found, using WHO GBD 2019 reference...")
    # WHO Global Burden of Disease 2019 reference life table
    age_groups = list(range(0, 100, 5))
    life_expectancy = [
        88.9, 84.0, 79.0, 74.1, 69.1, 64.1, 59.2, 54.2, 49.3, 44.4,
        39.5, 34.7, 30.0, 25.5, 21.2, 17.2, 13.6, 10.5, 7.9, 5.8
    ]
    life_table = pd.DataFrame({
        'age': age_groups,
        'ex_mean': life_expectancy
    })
    life_table_source = 'WHO GBD 2019 reference'
    print(f"  Using WHO reference life table: {len(life_table)} age groups")

# Create age-to-life-expectancy lookup function
def get_life_expectancy(age):
    """Get remaining life expectancy for a given age."""
    if life_table is None or len(life_table) == 0:
        return 10.0  # Fallback for missing data
    
    if 'age' not in life_table.columns:
        return 10.0
    
    if age in life_table['age'].values:
        return float(life_table[life_table['age'] == age]['ex_mean'].values[0])
    elif age > life_table['age'].max():
        return float(life_table['ex_mean'].min())
    else:
        # Interpolate
        idx = life_table['age'].searchsorted(age)
        if idx == 0:
            return float(life_table['ex_mean'].iloc[0])
        if idx >= len(life_table):
            return float(life_table['ex_mean'].iloc[-1])
        a1, a2 = life_table['age'].iloc[idx-1], life_table['age'].iloc[idx]
        e1, e2 = life_table['ex_mean'].iloc[idx-1], life_table['ex_mean'].iloc[idx]
        if a2 == a1:
            return float(e1)
        return float(e1 + (e2 - e1) * (age - a1) / (a2 - a1))


# =============================================================================
# AGE DISTRIBUTION OF ELDERLY DEATHS
# =============================================================================
print("\n[2] Setting up age distribution of elderly deaths...")

# Age distribution assumptions for elderly deaths in Brazil
# Based on typical mortality patterns
AGE_DISTRIBUTION = {
    'age_60_69': {'pct': 0.25, 'avg_age': 65, 'life_exp': get_life_expectancy(65)},
    'age_70_79': {'pct': 0.35, 'avg_age': 75, 'life_exp': get_life_expectancy(75)},
    'age_80_plus': {'pct': 0.40, 'avg_age': 85, 'life_exp': get_life_expectancy(85)}
}

print("  Age distribution of elderly deaths (assumed):")
for group, vals in AGE_DISTRIBUTION.items():
    print(f"    {group}: {vals['pct']*100:.0f}%, avg age {vals['avg_age']}, "
          f"life exp = {vals['life_exp']:.1f} years")

# Weighted average life expectancy for elderly deaths
avg_life_exp = sum(v['pct'] * v['life_exp'] for v in AGE_DISTRIBUTION.values())
print(f"\n  Weighted average remaining life expectancy: {avg_life_exp:.2f} years")
print(f"  Source: {life_table_source}")


# =============================================================================
# LOAD ATTRIBUTABLE BURDEN DATA (v2)
# =============================================================================
print("\n[3] Loading attributable burden data...")

# Try v2 first, fallback to v1
burden_summary_file = PHASE1_RESULTS / 'burden_v2_national_summary.json'
if not burden_summary_file.exists():
    burden_summary_file = PHASE1_RESULTS / 'attributable_burden_national.json'
    print(f"  Using v1 burden results (v2 not found)")
else:
    print(f"  Using v2 burden results")

with open(burden_summary_file, 'r') as f:
    burden_national = json.load(f)

# Load regional burden data
int_burden_file = PHASE1_RESULTS / 'burden_v2_intermediate_regions.csv'
imm_burden_file = PHASE1_RESULTS / 'burden_v2_immediate_regions.csv'

if not int_burden_file.exists():
    int_burden_file = PHASE1_RESULTS / 'burden_intermediate_regions.csv'
if not imm_burden_file.exists():
    imm_burden_file = PHASE1_RESULTS / 'burden_immediate_regions.csv'

burden_int = pd.read_csv(int_burden_file) if int_burden_file.exists() else None
burden_imm = pd.read_csv(imm_burden_file) if imm_burden_file.exists() else None

print(f"  Intermediate regions: {len(burden_int) if burden_int is not None else 0}")
print(f"  Immediate regions: {len(burden_imm) if burden_imm is not None else 0}")


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def safe_get(d, key, default=0):
    """Safely get a value from dict, returning default if missing or None."""
    val = d.get(key, default)
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return default
    return val


def safe_divide(numerator, denominator, default=0):
    """Safely divide, returning default if denominator is 0 or invalid."""
    if denominator is None or denominator == 0:
        return default
    if np.isnan(denominator):
        return default
    return numerator / denominator


def calculate_yll(attributable_deaths, avg_life_exp):
    """Calculate Years of Life Lost with guard."""
    if attributable_deaths is None or np.isnan(attributable_deaths):
        return 0.0
    if avg_life_exp is None or np.isnan(avg_life_exp):
        avg_life_exp = 10.0  # Fallback
    return float(attributable_deaths * avg_life_exp)


# =============================================================================
# CALCULATE YLL
# =============================================================================
print("\n[4] Calculating Years of Life Lost...")

yll_results = {
    'method': {
        'description': 'YLL = Attributable deaths × Remaining life expectancy',
        'age_distribution': {k: {kk: float(vv) if isinstance(vv, (int, float)) else vv 
                                  for kk, vv in v.items()} 
                              for k, v in AGE_DISTRIBUTION.items()},
        'avg_life_expectancy': float(avg_life_exp),
        'source': life_table_source
    },
    'intermediate': {},
    'immediate': {},
    'timestamp': datetime.now().isoformat()
}

# Process each level
for level_name, nat_data in [('intermediate', burden_national.get('intermediate', {})),
                              ('immediate', burden_national.get('immediate', {}))]:
    if not nat_data:
        print(f"  {level_name.upper()}: No data available")
        continue
    
    print(f"\n  {level_name.upper()}:")
    
    # Get attributable deaths (try v2 field names first, then v1)
    # V2 uses P2.5/P97.5 as primary thresholds
    heat_an = safe_get(nat_data, 'heat_an_97_5', safe_get(nat_data, 'heat_an', 0))
    cold_an = safe_get(nat_data, 'cold_an_2_5', safe_get(nat_data, 'cold_an', 0))
    
    # Also get total attributable (all days above/below MMT)
    total_heat_an = safe_get(nat_data, 'total_heat_an', heat_an)
    total_cold_an = safe_get(nat_data, 'total_cold_an', cold_an)
    
    print(f"    Heat deaths (>P97.5): {heat_an:,.0f}")
    print(f"    Cold deaths (<P2.5): {cold_an:,.0f}")
    print(f"    Total heat (>MMT): {total_heat_an:,.0f}")
    print(f"    Total cold (<MMT): {total_cold_an:,.0f}")
    
    # Calculate YLL
    heat_yll = calculate_yll(heat_an, avg_life_exp)
    cold_yll = calculate_yll(cold_an, avg_life_exp)
    total_heat_yll = calculate_yll(total_heat_an, avg_life_exp)
    total_cold_yll = calculate_yll(total_cold_an, avg_life_exp)
    
    print(f"    Heat YLL (>P97.5): {heat_yll:,.0f}")
    print(f"    Cold YLL (<P2.5): {cold_yll:,.0f}")
    
    # Get n_years (with guard)
    avg_years = safe_get(nat_data, 'avg_years', 0)
    if avg_years <= 0:
        # Try to estimate from annual values
        heat_annual = safe_get(nat_data, 'heat_annual_97_5', safe_get(nat_data, 'heat_annual', 0))
        if heat_annual > 0 and heat_an > 0:
            avg_years = heat_an / heat_annual
        else:
            avg_years = 14.0  # Default assumption
    
    print(f"    Years of data: {avg_years:.1f}")
    
    # Annual values
    heat_yll_annual = safe_divide(heat_yll, avg_years, 0)
    cold_yll_annual = safe_divide(cold_yll, avg_years, 0)
    
    # Population and rates
    total_pop = safe_get(nat_data, 'total_pop', None)
    if total_pop and total_pop > 0:
        heat_yll_rate = heat_yll_annual / total_pop * 100000
        cold_yll_rate = cold_yll_annual / total_pop * 100000
        print(f"    Heat YLL rate: {heat_yll_rate:.1f} per 100k/year")
        print(f"    Cold YLL rate: {cold_yll_rate:.1f} per 100k/year")
    else:
        heat_yll_rate = None
        cold_yll_rate = None
    
    # Store results
    yll_results[level_name] = {
        # Primary thresholds (P2.5/P97.5)
        'heat_an': float(heat_an),
        'cold_an': float(cold_an),
        'heat_yll_total': float(heat_yll),
        'cold_yll_total': float(cold_yll),
        'heat_yll_annual': float(heat_yll_annual),
        'cold_yll_annual': float(cold_yll_annual),
        'heat_yll_rate_per_100k': float(heat_yll_rate) if heat_yll_rate else None,
        'cold_yll_rate_per_100k': float(cold_yll_rate) if cold_yll_rate else None,
        
        # Total attributable (all days above/below MMT)
        'total_heat_an': float(total_heat_an),
        'total_cold_an': float(total_cold_an),
        'total_heat_yll': float(total_heat_yll),
        'total_cold_yll': float(total_cold_yll),
        
        # Metadata
        'n_years': float(avg_years),
        'total_pop': float(total_pop) if total_pop else None,
        'n_regions': int(safe_get(nat_data, 'n_regions', 0))
    }


# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("YEARS OF LIFE LOST (YLL) SUMMARY")
print("=" * 70)

for level in ['intermediate', 'immediate']:
    data = yll_results[level]
    if not data:
        continue
    
    suffix = " (PRIMARY)" if level == 'immediate' else ""
    print(f"\n{level.upper()}{suffix}")
    print("-" * 50)
    
    print(f"  Extreme heat (>P97.5):")
    print(f"    Deaths: {data['heat_an']:,.0f}")
    print(f"    YLL total: {data['heat_yll_total']:,.0f}")
    print(f"    YLL annual: {data['heat_yll_annual']:,.0f}")
    if data.get('heat_yll_rate_per_100k'):
        print(f"    YLL rate: {data['heat_yll_rate_per_100k']:.1f} per 100k elderly/year")
    
    print(f"\n  Extreme cold (<P2.5):")
    print(f"    Deaths: {data['cold_an']:,.0f}")
    print(f"    YLL total: {data['cold_yll_total']:,.0f}")
    print(f"    YLL annual: {data['cold_yll_annual']:,.0f}")
    if data.get('cold_yll_rate_per_100k'):
        print(f"    YLL rate: {data['cold_yll_rate_per_100k']:.1f} per 100k elderly/year")
    
    print(f"\n  Total attributable to non-optimal temperature:")
    print(f"    Heat (>MMT): {data['total_heat_an']:,.0f} deaths, {data['total_heat_yll']:,.0f} YLL")
    print(f"    Cold (<MMT): {data['total_cold_an']:,.0f} deaths, {data['total_cold_yll']:,.0f} YLL")
    
    total_yll = data['heat_yll_total'] + data['cold_yll_total']
    total_all_yll = data['total_heat_yll'] + data['total_cold_yll']
    print(f"\n  Total YLL (extreme only): {total_yll:,.0f}")
    print(f"  Total YLL (all non-optimal): {total_all_yll:,.0f}")


# =============================================================================
# SAVE RESULTS
# =============================================================================
print("\n" + "=" * 70)
print("SAVING RESULTS")
print("=" * 70)

out_file = OUTPUT_DIR / 'yll_v2_results.json'
with open(out_file, 'w') as f:
    json.dump(yll_results, f, indent=2)
print(f"  Saved: {out_file}")

# Summary CSV
summary_rows = []
for level in ['intermediate', 'immediate']:
    data = yll_results[level]
    if not data:
        continue
    summary_rows.append({
        'level': level,
        'n_regions': data.get('n_regions', 0),
        'heat_deaths': data['heat_an'],
        'cold_deaths': data['cold_an'],
        'heat_yll': data['heat_yll_total'],
        'cold_yll': data['cold_yll_total'],
        'total_yll': data['heat_yll_total'] + data['cold_yll_total'],
        'heat_yll_annual': data['heat_yll_annual'],
        'cold_yll_annual': data['cold_yll_annual'],
        'heat_yll_rate': data.get('heat_yll_rate_per_100k'),
        'cold_yll_rate': data.get('cold_yll_rate_per_100k'),
    })

summary_df = pd.DataFrame(summary_rows)
out_file = OUTPUT_DIR / 'yll_v2_summary.csv'
summary_df.to_csv(out_file, index=False)
print(f"  Saved: {out_file}")

print("\n" + "=" * 70)
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)
print("\nNotes:")
print("  - YLL based on IBGE Brazilian life tables")
print("  - Using P2.5/P97.5 thresholds (primary)")
print("  - Total also includes all days above/below MMT")
