"""
01e_yll_unified.py
==================
Unified Years of Life Lost (YLL) Calculation

MERGES:
- 01e_yll_calculation.py (assumed age distribution)
- 01e2_yll_with_age_distribution.py (actual ages from SIM)

Key Features:
1. Single source of truth for life table loading
2. Robust SIM age parsing with quality diagnostics
3. CLI flags for mode selection (--mode actual|assumed)
4. Explicit documentation of assumptions vs data
5. Comparison output between modes for sensitivity analysis

Usage:
    python 01e_yll_unified.py --mode actual      # Use actual SIM ages (default, preferred)
    python 01e_yll_unified.py --mode assumed     # Use assumed age distribution
    python 01e_yll_unified.py --mode both        # Run both and compare

References:
- WHO YLL methodology
- Gasparrini et al. (2015) - GBD YLL approach
- IBGE Brazilian life tables

Author: Climate-Health Analysis Pipeline
Date: December 2025
"""

import warnings
warnings.filterwarnings('ignore')

import argparse
import json
import re
import logging
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
from glob import glob

# =============================================================================
# LOGGING SETUP
# =============================================================================

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)

# =============================================================================
# PATHS
# =============================================================================

SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent.parent
INPUT_DATA = BASE_DIR / 'Input_data'
PHASE0_RESULTS = SCRIPT_DIR.parent / 'phase0_data_prep' / 'results'
PHASE1_RESULTS = SCRIPT_DIR / 'results'
OUTPUT_DIR = PHASE1_RESULTS
OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# LIFE TABLE LOADING (SINGLE SOURCE OF TRUTH)
# =============================================================================

def load_life_table(custom_path=None):
    """
    Load life table with consistent behavior across all modes.
    
    Priority:
    1. Custom path if provided
    2. IBGE Brazil-specific life table (yll_lookup_by_age.csv)
    3. IBGE combined tables (ibge_life_tables_combined.csv)
    4. WHO GBD 2019 reference (fallback)
    
    Returns:
    --------
    life_table : DataFrame with columns ['age', 'ex']
    source : str describing the source
    """
    # Try custom path first
    if custom_path and Path(custom_path).exists():
        lt = pd.read_csv(custom_path)
        # Normalize column names
        lt.columns = lt.columns.str.lower()
        if 'ex_mean' in lt.columns:
            lt = lt.rename(columns={'ex_mean': 'ex'})
        lt['age'] = lt['age'].astype(int)
        lt = lt.sort_values('age').reset_index(drop=True)
        logger.info(f"Loaded custom life table from {custom_path}")
        return lt[['age', 'ex']], f"Custom: {custom_path}"
    
    # Try IBGE YLL lookup
    ibge_yll = PHASE0_RESULTS / 'yll_lookup_by_age.csv'
    if ibge_yll.exists():
        lt = pd.read_csv(ibge_yll)
        lt.columns = lt.columns.str.lower()
        if 'ex_mean' in lt.columns:
            lt = lt.rename(columns={'ex_mean': 'ex'})
        lt['age'] = lt['age'].astype(int)
        lt = lt.sort_values('age').reset_index(drop=True)
        logger.info(f"Loaded IBGE YLL lookup: {len(lt)} ages")
        return lt[['age', 'ex']], "IBGE Brazil life table"
    
    # Try IBGE combined
    ibge_combined = PHASE0_RESULTS / 'ibge_life_tables_combined.csv'
    if ibge_combined.exists():
        lt = pd.read_csv(ibge_combined)
        lt.columns = lt.columns.str.lower()
        if 'ex_mean' in lt.columns:
            lt = lt.rename(columns={'ex_mean': 'ex'})
        lt['age'] = lt['age'].astype(int)
        lt = lt.sort_values('age').reset_index(drop=True)
        logger.info(f"Loaded IBGE combined life tables")
        return lt[['age', 'ex']], "IBGE Brazil combined"
    
    # Fallback to WHO GBD 2019 reference
    logger.warning("IBGE life tables not found. Using WHO GBD 2019 reference.")
    age_groups = list(range(0, 100, 5))
    life_expectancy = [
        88.9, 84.0, 79.0, 74.1, 69.1, 64.1, 59.2, 54.2, 49.3, 44.4,
        39.5, 34.7, 30.0, 25.5, 21.2, 17.2, 13.6, 10.5, 7.9, 5.8
    ]
    lt = pd.DataFrame({'age': age_groups, 'ex': life_expectancy})
    return lt, "WHO GBD 2019 reference (fallback)"


def get_life_expectancy_at_age(age, life_table):
    """Get remaining life expectancy for a specific age.

    Uses linear interpolation for ages not exactly in the table.
    """
    if life_table is None or len(life_table) == 0:
        return 10.0  # Ultimate fallback

    age = int(age)
    age_vals = life_table['age'].values
    ex_vals = life_table['ex'].values

    # Exact match
    if age in age_vals:
        return float(life_table[life_table['age'] == age]['ex'].values[0])

    # Beyond range
    if age > age_vals.max():
        return float(ex_vals[-1])  # Use minimum (oldest age)
    if age < age_vals.min():
        return float(ex_vals[0])

    # Interpolate within range
    return float(np.interp(age, age_vals, ex_vals))


# =============================================================================
# ROBUST SIM AGE PARSING (ALIGNED WITH PHASE 0 AGGREGATION)
# =============================================================================

def parse_sim_age(raw):
    """Parse SIM IDADE field to age in years (Phase 0–consistent).

    Mirrors the canonical logic used in
    `phase0_data_prep/aggregation/00n_aggregate_mortality_to_regions.py`.

    IDADE is coded as XYY where X is the unit and YY is the value:
    - 0YY: minutes
    - 1YY: hours
    - 2YY: days
    - 3YY: months
    - 4YY: years (e.g., 465 -> 65 years)
    - 5YY: 100+ years (e.g., 505 -> 105 years)

    Implementation in Phase 0:
    - If IDADE >= 400: age_years = IDADE - 400
    - Else: age_years = 0 (collapsed infants/children for adult analyses)
    """

    if pd.isna(raw):
        return None

    # Normalize to an integer code, tolerating strings and decimals
    try:
        s = str(raw).strip()
        if '.' in s:
            s = s.split('.')[0]
        digits = re.sub(r"\D", "", s)
        if digits == "":
            return None
        code = int(digits)
    except Exception:
        return None

    if code >= 400:
        # Years: 4YY (0–99) and 5YY (100+) and any higher codes treated
        # the same way as in Phase 0 (age = code - 400).
        return code - 400

    # For codes < 400, Phase 0 collapses these to 0 (under 1 year).
    # For elderly (60+), this distinction does not matter; we return 0
    # rather than None so that "parse success" reflects Phase 0 logic.
    if code >= 0:
        return 0

    return None


def extract_ages_from_sim(sim_dir=None, min_elderly_age=60):
    """
    Extract actual age distribution from SIM mortality data.

    Returns:
    --------
    ages_series : pd.Series of ages for each elderly death
    parse_stats : dict with parsing diagnostics
    """
    if sim_dir is None:
        sim_dir = INPUT_DATA

    sim_dir = Path(sim_dir)

    # Find all DO files
    do_files = sorted(glob(str(sim_dir / 'DO*.csv')))

    if not do_files:
        logger.warning(f"No DO*.csv files found in {sim_dir}")
        return None, {'success_rate': 0, 'total': 0}

    logger.info(f"Found {len(do_files)} SIM mortality files")

    all_ages = []
    total_rows = 0
    parsed_rows = 0
    unparsed_examples = []
    years_processed = []

    for file_path in do_files:
        file_name = Path(file_path).name

        # Extract year from filename (DO10OPEN.csv -> 2010)
        try:
            year_code = file_name[2:4]
            year = 2000 + int(year_code)
        except Exception:
            year = 0

        try:
            # Auto-detect separator (DO24 uses comma, earlier years use semicolon)
            with open(file_path, 'r', encoding='latin-1') as fh:
                first_line = fh.readline()
            sep = ',' if ('",' in first_line or first_line.count(',') > first_line.count(';')) else ';'
            
            # Read file
            df = pd.read_csv(file_path, sep=sep, encoding='latin-1', low_memory=False)

            # Find age column
            age_col = None
            for col in df.columns:
                if 'IDADE' in col.upper():
                    age_col = col
                    break

            if age_col is None:
                logger.warning(f"{file_name}: No IDADE column found")
                continue

            # Parse ages using Phase 0–consistent logic
            df['age_parsed'] = df[age_col].apply(parse_sim_age)

            total_rows += len(df)
            n_parsed = df['age_parsed'].notna().sum()
            parsed_rows += n_parsed

            # Keep only elderly
            elderly = df[df['age_parsed'] >= min_elderly_age]['age_parsed'].dropna()
            all_ages.extend(elderly.tolist())

            # Save unparsed examples for debugging
            unparsed = df[df['age_parsed'].isna()][age_col].head(10).tolist()
            unparsed_examples.extend(unparsed)

            years_processed.append(year)
            logger.info(
                f"  {file_name}: {len(elderly):,} elderly deaths "
                f"(parse rate: {100 * n_parsed / len(df):.1f}%)"
            )

        except Exception as e:
            logger.error(f"  {file_name}: Error - {e}")
            continue

    if not all_ages:
        return None, {'success_rate': 0, 'total': 0}

    # Parse statistics
    success_rate = 100 * parsed_rows / total_rows if total_rows > 0 else 0

    parse_stats = {
        'total_rows': total_rows,
        'parsed_rows': parsed_rows,
        'success_rate': success_rate,
        'n_elderly_deaths': len(all_ages),
        'years_processed': years_processed,
        'year_range': f"{min(years_processed)}-{max(years_processed)}" if years_processed else "N/A",
        'unparsed_examples': unparsed_examples[:20],  # Save some examples
    }

    # Log quality check
    if success_rate < 95:
        logger.warning(f"Age parsing success rate is low: {success_rate:.1f}%")
        logger.warning("Consider using assumed distribution or inspecting SIM data")

        # Save debug file
        debug_file = OUTPUT_DIR / 'debug_unparsed_ages_sample.csv'
        pd.DataFrame({'raw_age': unparsed_examples[:50]}).to_csv(debug_file, index=False)
        logger.info(f"Saved unparsed examples to {debug_file}")
    else:
        logger.info(f"Age parsing success rate: {success_rate:.1f}% \u2713")

    # Clean implausible ages for YLL calculation to avoid overflow and
    # exclude obviously invalid SIM codes (e.g., ages >> 120 years).
    clean_ages = []
    dropped_out_of_range = 0
    for a in all_ages:
        if a is None:
            continue
        try:
            val = float(a)
        except Exception:
            dropped_out_of_range += 1
            continue

        # Keep only plausible human ages
        if val < 0 or val > 120:
            dropped_out_of_range += 1
            continue

        clean_ages.append(int(round(val)))

    if dropped_out_of_range > 0:
        logger.warning(
            f"Dropped {dropped_out_of_range} implausible ages outside [0, 120] "
            f"for YLL calculation out of {len(all_ages)} elderly deaths (SIM anomalies)."
        )

    parse_stats['n_ages_used_for_yll'] = len(clean_ages)

    return pd.Series(clean_ages, dtype=int), parse_stats


# =============================================================================
# YLL COMPUTATION ENGINE (SINGLE SOURCE)
# =============================================================================

def compute_yll_from_ages(ages, life_table, group_ages=True):
    """
    Compute YLL from a series of ages at death.
    
    Parameters:
    -----------
    ages : array-like of integer ages
    life_table : DataFrame with columns ['age', 'ex']
    group_ages : bool, if True also return breakdown by age group
    
    Returns:
    --------
    dict with total_yll, n_deaths, avg_life_exp, and optionally age_group_breakdown
    """
    ages = np.array([a for a in ages if a is not None and not np.isnan(a)])
    
    if len(ages) == 0:
        return {'total_yll': 0, 'n_deaths': 0, 'avg_life_exp': 0}
    
    # Get life expectancy for each age
    age_vals = life_table['age'].values
    ex_vals = life_table['ex'].values
    ex_at_age = np.interp(ages, age_vals, ex_vals, left=ex_vals[0], right=ex_vals[-1])
    
    total_yll = float(np.sum(ex_at_age))
    avg_life_exp = float(np.mean(ex_at_age))
    
    result = {
        'total_yll': total_yll,
        'n_deaths': int(len(ages)),
        'avg_life_exp': avg_life_exp,
        'yll_per_death': total_yll / len(ages) if len(ages) > 0 else 0
    }
    
    if group_ages:
        # Breakdown by age group
        age_groups = [
            (60, 64, '60-64'),
            (65, 69, '65-69'),
            (70, 74, '70-74'),
            (75, 79, '75-79'),
            (80, 84, '80-84'),
            (85, 89, '85-89'),
            (90, 120, '90+')
        ]
        
        breakdown = []
        for start, end, label in age_groups:
            mask = (ages >= start) & (ages <= end)
            group_ages_arr = ages[mask]
            group_ex = ex_at_age[mask]
            
            if len(group_ages_arr) > 0:
                breakdown.append({
                    'age_group': label,
                    'n_deaths': int(len(group_ages_arr)),
                    'pct_deaths': float(100 * len(group_ages_arr) / len(ages)),
                    'total_yll': float(np.sum(group_ex)),
                    'avg_life_exp': float(np.mean(group_ex))
                })
        
        result['age_group_breakdown'] = breakdown
    
    return result


def compute_yll_assumed(n_deaths, life_table, assumed_distribution=None):
    """
    Compute YLL using assumed age distribution.
    
    Parameters:
    -----------
    n_deaths : int, total number of deaths
    life_table : DataFrame with ['age', 'ex']
    assumed_distribution : dict mapping age group -> proportion
    
    Returns:
    --------
    dict with total_yll, breakdown, and explicit documentation that this is assumed
    """
    if assumed_distribution is None:
        # Default assumed distribution for elderly deaths
        # Based on typical mortality patterns in Brazil
        assumed_distribution = {
            65: 0.10,   # 60-69: 25% (split into 65)
            70: 0.15,   # 
            75: 0.25,   # 70-79: 35%
            80: 0.25,   # 80-84: 25%
            85: 0.15,   # 85-89: 15%
            90: 0.10    # 90+: 10%
        }
    
    total_yll = 0
    breakdown = []
    
    for age, proportion in assumed_distribution.items():
        group_deaths = n_deaths * proportion
        life_exp = get_life_expectancy_at_age(age, life_table)
        group_yll = group_deaths * life_exp
        total_yll += group_yll
        
        breakdown.append({
            'representative_age': age,
            'proportion': proportion,
            'n_deaths': group_deaths,
            'life_exp': life_exp,
            'yll': group_yll
        })
    
    avg_life_exp = total_yll / n_deaths if n_deaths > 0 else 0
    
    return {
        'total_yll': total_yll,
        'n_deaths': n_deaths,
        'avg_life_exp': avg_life_exp,
        'yll_per_death': avg_life_exp,
        'assumed_distribution': assumed_distribution,
        'breakdown': breakdown,
        'is_assumed': True,
        'warning': 'YLL computed using ASSUMED age distribution, not actual ages'
    }


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_yll_analysis(mode='actual', life_table_path=None, sim_dir=None):
    """
    Run YLL analysis in specified mode.
    
    Parameters:
    -----------
    mode : str, one of 'actual', 'assumed', 'both'
    life_table_path : str, optional custom life table path
    sim_dir : str, optional custom SIM data directory
    """
    print("=" * 70)
    print(f"01e: YEARS OF LIFE LOST (YLL) CALCULATION")
    print(f"     Mode: {mode.upper()}")
    print("=" * 70)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Load life table (same for all modes)
    print("\n[1] Loading life table...")
    life_table, lt_source = load_life_table(life_table_path)
    print(f"    Source: {lt_source}")
    print(f"    Ages: {life_table['age'].min()}-{life_table['age'].max()}")
    print(f"    Life exp at 60: {get_life_expectancy_at_age(60, life_table):.1f} years")
    print(f"    Life exp at 80: {get_life_expectancy_at_age(80, life_table):.1f} years")
    
    # Load burden data
    print("\n[2] Loading attributable burden data...")
    
    burden_file = PHASE1_RESULTS / 'burden_v2_national_summary.json'
    if not burden_file.exists():
        burden_file = PHASE1_RESULTS / 'attributable_burden_national.json'
    
    if not burden_file.exists():
        logger.error("No burden file found. Run 01d_attributable_burden first.")
        return None
    
    with open(burden_file, 'r') as f:
        burden_data = json.load(f)
    
    print(f"    Loaded: {burden_file.name}")
    
    results = {
        'mode': mode,
        'life_table_source': lt_source,
        'timestamp': datetime.now().isoformat(),
        'actual': None,
        'assumed': None,
        'comparison': None
    }
    
    # === ACTUAL AGE MODE ===
    if mode in ['actual', 'both']:
        print("\n[3a] Computing YLL using ACTUAL ages from SIM data...")
        
        ages, parse_stats = extract_ages_from_sim(sim_dir)
        
        if ages is not None and len(ages) > 0:
            # Quality check
            if parse_stats['success_rate'] < 95:
                logger.warning(f"Low parse rate ({parse_stats['success_rate']:.1f}%). Results may be biased.")
            
            # Compute overall YLL per death
            yll_result = compute_yll_from_ages(ages, life_table, group_ages=True)
            
            print(f"\n    Elderly deaths: {yll_result['n_deaths']:,}")
            print(f"    Average life expectancy: {yll_result['avg_life_exp']:.2f} years")
            print(f"    YLL per death: {yll_result['yll_per_death']:.2f}")
            
            # Apply to attributable deaths from burden data
            for level in ['intermediate', 'immediate']:
                level_data = burden_data.get(level, {})
                if not level_data:
                    continue
                
                # Get attributable deaths (try v2 names first)
                heat_an = level_data.get('heat_an_97_5', level_data.get('heat_an', 0))
                cold_an = level_data.get('cold_an_2_5', level_data.get('cold_an', 0))
                
                if heat_an is None: heat_an = 0
                if cold_an is None: cold_an = 0
                
                # Compute YLL using per-death rate from actual ages
                yll_per_death = yll_result['yll_per_death']
                
                level_yll = {
                    'heat_an': float(heat_an),
                    'cold_an': float(cold_an),
                    'heat_yll': float(heat_an * yll_per_death),
                    'cold_yll': float(cold_an * yll_per_death),
                    'total_yll': float((heat_an + cold_an) * yll_per_death),
                    'yll_per_death': float(yll_per_death)
                }
                
                yll_result[f'{level}_yll'] = level_yll
            
            yll_result['parse_stats'] = parse_stats
            results['actual'] = yll_result
            
        else:
            logger.warning("Could not extract ages from SIM data. Falling back to assumed mode.")
            mode = 'assumed' if mode == 'actual' else mode
    
    # === ASSUMED MODE ===
    if mode in ['assumed', 'both']:
        print("\n[3b] Computing YLL using ASSUMED age distribution...")
        
        # Get total deaths from burden data
        for level in ['intermediate', 'immediate']:
            level_data = burden_data.get(level, {})
            if not level_data:
                continue
            
            heat_an = level_data.get('heat_an_97_5', level_data.get('heat_an', 0))
            cold_an = level_data.get('cold_an_2_5', level_data.get('cold_an', 0))
            
            if heat_an is None: heat_an = 0
            if cold_an is None: cold_an = 0
            
            total_an = heat_an + cold_an
            
            yll_assumed = compute_yll_assumed(total_an, life_table)
            
            heat_yll = heat_an * yll_assumed['yll_per_death']
            cold_yll = cold_an * yll_assumed['yll_per_death']
            
            if results['assumed'] is None:
                results['assumed'] = {
                    'avg_life_exp': yll_assumed['avg_life_exp'],
                    'yll_per_death': yll_assumed['yll_per_death'],
                    'assumed_distribution': yll_assumed['assumed_distribution'],
                    'is_assumed': True
                }
            
            results['assumed'][f'{level}_yll'] = {
                'heat_an': float(heat_an),
                'cold_an': float(cold_an),
                'heat_yll': float(heat_yll),
                'cold_yll': float(cold_yll),
                'total_yll': float(heat_yll + cold_yll)
            }
        
        print(f"    Assumed YLL per death: {results['assumed']['yll_per_death']:.2f}")
    
    # === COMPARISON ===
    if mode == 'both' and results['actual'] and results['assumed']:
        print("\n[4] Comparing actual vs assumed modes...")
        
        comparison = {
            'yll_per_death_actual': results['actual']['yll_per_death'],
            'yll_per_death_assumed': results['assumed']['yll_per_death'],
            'difference_pct': 100 * (results['actual']['yll_per_death'] - results['assumed']['yll_per_death']) / results['assumed']['yll_per_death']
        }
        
        for level in ['intermediate', 'immediate']:
            actual_yll = results['actual'].get(f'{level}_yll', {}).get('total_yll', 0)
            assumed_yll = results['assumed'].get(f'{level}_yll', {}).get('total_yll', 0)
            
            if assumed_yll > 0:
                comparison[f'{level}_diff_pct'] = 100 * (actual_yll - assumed_yll) / assumed_yll
            
        results['comparison'] = comparison
        
        print(f"    YLL per death (actual): {comparison['yll_per_death_actual']:.2f}")
        print(f"    YLL per death (assumed): {comparison['yll_per_death_assumed']:.2f}")
        print(f"    Difference: {comparison['difference_pct']:+.1f}%")
    
    # === SUMMARY ===
    print("\n" + "=" * 70)
    print("YLL SUMMARY")
    print("=" * 70)
    
    primary_result = results['actual'] if results['actual'] else results['assumed']
    method_label = "ACTUAL ages" if results['actual'] else "ASSUMED distribution"
    
    print(f"\nMethod: {method_label}")
    print(f"Life table: {lt_source}")
    print(f"YLL per elderly death: {primary_result['yll_per_death']:.2f} years")
    
    for level in ['intermediate', 'immediate']:
        level_yll = primary_result.get(f'{level}_yll', {})
        if level_yll:
            suffix = " (PRIMARY)" if level == 'immediate' else ""
            print(f"\n{level.upper()}{suffix}:")
            print(f"  Heat deaths: {level_yll.get('heat_an', 0):,.0f}")
            print(f"  Cold deaths: {level_yll.get('cold_an', 0):,.0f}")
            print(f"  Heat YLL: {level_yll.get('heat_yll', 0):,.0f}")
            print(f"  Cold YLL: {level_yll.get('cold_yll', 0):,.0f}")
            print(f"  Total YLL: {level_yll.get('total_yll', 0):,.0f}")
    
    # === SAVE RESULTS ===
    print("\n" + "=" * 70)
    print("SAVING RESULTS")
    print("=" * 70)
    
    out_file = OUTPUT_DIR / 'yll_unified_results.json'
    
    # Clean for JSON
    def clean_for_json(obj):
        if isinstance(obj, (np.integer, np.int64, np.int32)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, pd.Series):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: clean_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [clean_for_json(i) for i in obj]
        return obj
    
    with open(out_file, 'w') as f:
        json.dump(clean_for_json(results), f, indent=2)
    print(f"  Saved: {out_file}")
    
    # Summary CSV
    summary_rows = []
    for method in ['actual', 'assumed']:
        method_data = results.get(method)
        if not method_data:
            continue
        for level in ['intermediate', 'immediate']:
            level_yll = method_data.get(f'{level}_yll', {})
            if level_yll:
                summary_rows.append({
                    'method': method,
                    'level': level,
                    'heat_deaths': level_yll.get('heat_an', 0),
                    'cold_deaths': level_yll.get('cold_an', 0),
                    'heat_yll': level_yll.get('heat_yll', 0),
                    'cold_yll': level_yll.get('cold_yll', 0),
                    'total_yll': level_yll.get('total_yll', 0),
                    'yll_per_death': method_data.get('yll_per_death', 0)
                })
    
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_file = OUTPUT_DIR / 'yll_unified_summary.csv'
        summary_df.to_csv(summary_file, index=False)
        print(f"  Saved: {summary_file}")
    
    print("\n" + "=" * 70)
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)
    
    return results


# =============================================================================
# CLI INTERFACE
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Compute Years of Life Lost (YLL) for temperature-attributable deaths',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python 01e_yll_unified.py --mode actual
  python 01e_yll_unified.py --mode assumed
  python 01e_yll_unified.py --mode both
  python 01e_yll_unified.py --mode actual --life_table path/to/custom_lt.csv
        """
    )
    
    parser.add_argument(
        '--mode',
        choices=['actual', 'assumed', 'both'],
        default='both',
        help='YLL computation mode: actual (from SIM ages), assumed (fixed distribution), or both (for comparison)'
    )
    
    parser.add_argument(
        '--life_table',
        type=str,
        default=None,
        help='Path to custom life table CSV (must have columns: age, ex)'
    )
    
    parser.add_argument(
        '--sim_dir',
        type=str,
        default=None,
        help='Path to directory containing SIM DO*.csv files'
    )
    
    args = parser.parse_args()
    
    run_yll_analysis(
        mode=args.mode,
        life_table_path=args.life_table,
        sim_dir=args.sim_dir
    )


if __name__ == '__main__':
    main()
