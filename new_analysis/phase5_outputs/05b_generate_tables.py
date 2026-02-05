"""
05b: GENERATE TABLES FOR ALL PHASES
====================================
Creates publication-ready tables from Phase 1-4 results.

Outputs:
--------
- Table 1: Main effects summary (both levels)
- Table 2: Attributable burden by region
- Table 3: Sensitivity analyses summary
- Table 4: Confounder adjustment comparison
- Table 5: Stratified analyses (age, sex, cause)
- Table 6: Meta-regression coefficients

Author: Heat-Mortality Brazil Analysis Pipeline
Date: December 2025
"""

import warnings
warnings.filterwarnings('ignore')

import json
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = Path(__file__).parent.parent
PHASE1_RESULTS = BASE_DIR / 'phase1_core_model' / 'results'
PHASE2_RESULTS = BASE_DIR / 'phase2_robustness' / 'results'
PHASE3_RESULTS = BASE_DIR / 'phase3_confounding' / 'results'
PHASE4_RESULTS = BASE_DIR / 'phase4_heterogeneity' / 'results'
OUTPUT_DIR = Path(__file__).parent / 'tables'
OUTPUT_DIR.mkdir(exist_ok=True)

print("="*70)
print("05b: GENERATE TABLES FOR ALL PHASES")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_death_column(df):
    """Identify the death column in a dataframe."""
    for col in ['deaths_all', 'deaths_elderly', 'deaths', 'total_deaths']:
        if col in df.columns:
            return col
    raise KeyError(f"Column not found: deaths. Available: {list(df.columns)}")

def get_region_column(df, level='intermediate'):
    """Identify the region column in a dataframe."""
    candidates = ['region_code', f'{level}_code', 'intermediate_code', 'immediate_code']
    for col in candidates:
        if col in df.columns:
            return col
    raise KeyError(f"Region column not found. Available: {list(df.columns)}")

def load_json(filepath):
    """Load JSON file safely."""
    try:
        with open(filepath, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"  Warning: Could not load {filepath}: {e}")
        return None

def format_rr(rr, lo, hi, decimals=3):
    """Format RR with 95% CI."""
    if rr is None or np.isnan(rr):
        return "—"
    return f"{rr:.{decimals}f} ({lo:.{decimals}f}–{hi:.{decimals}f})"

def format_number(n, decimals=0):
    """Format number with thousands separator."""
    if n is None or (isinstance(n, float) and np.isnan(n)):
        return "—"
    if decimals == 0:
        return f"{int(n):,}"
    return f"{n:,.{decimals}f}"

def save_table(df, name, formats=['csv', 'xlsx']):
    """Save table in multiple formats."""
    for fmt in formats:
        filepath = OUTPUT_DIR / f"{name}.{fmt}"
        if fmt == 'csv':
            df.to_csv(filepath, index=False)
        elif fmt == 'xlsx':
            df.to_excel(filepath, index=False)
    print(f"  Saved: {name}")

# =============================================================================
# DESCRIPTIVE TABLES
# =============================================================================

INPUT_DATA_DIR = BASE_DIR.parent / 'Input_data'
PHASE0_RESULTS = BASE_DIR / 'phase0_data_prep' / 'results'

def create_table_descriptive_mortality():
    """Create descriptive statistics for mortality data."""
    print("\n[Descriptive] Mortality statistics...")
    
    rows = []
    
    for level in ['intermediate', 'immediate']:
        # Load mortality summary
        mort_file = PHASE0_RESULTS / f'mortality_{level}_daily.parquet'
        if not mort_file.exists():
            continue
        
        try:
            mort_df = pd.read_parquet(mort_file)
            mort_df['date'] = pd.to_datetime(mort_df['date'])
            
            # Identify death column
            death_col = get_death_column(mort_df)
            region_col = get_region_column(mort_df, level)
            
            level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
            
            # Overall statistics
            total_deaths = mort_df[death_col].sum()
            n_regions = mort_df[region_col].nunique()
            n_days = mort_df['date'].nunique()
            date_range = f"{mort_df['date'].min().strftime('%Y-%m-%d')} to {mort_df['date'].max().strftime('%Y-%m-%d')}"
            
            # Daily deaths per region
            daily_mean = mort_df[death_col].mean()
            daily_sd = mort_df[death_col].std()
            daily_min = mort_df[death_col].min()
            daily_max = mort_df[death_col].max()
            
            # Annual statistics
            mort_df['year'] = mort_df['date'].dt.year
            annual_deaths = mort_df.groupby('year')[death_col].sum()
            
            rows.append({
                'Level': level_label,
                'N Regions': n_regions,
                'N Days': n_days,
                'Date Range': date_range,
                'Total Deaths (Elderly)': format_number(total_deaths),
                'Daily Deaths Mean (SD)': f"{daily_mean:.1f} ({daily_sd:.1f})",
                'Daily Deaths Range': f"{daily_min} – {daily_max}",
                'Annual Deaths Mean': format_number(annual_deaths.mean()),
                'Annual Deaths Range': f"{format_number(annual_deaths.min())} – {format_number(annual_deaths.max())}",
            })
            
        except Exception as e:
            print(f"  Error loading {level}: {e}")
    
    if rows:
        df = pd.DataFrame(rows)
        save_table(df, 'tableD1_mortality_descriptive')
        return df
    return None

def create_table_descriptive_temperature():
    """Create descriptive statistics for temperature data."""
    print("\n[Descriptive] Temperature statistics...")
    
    rows = []
    
    for level in ['intermediate', 'immediate']:
        temp_file = PHASE0_RESULTS / f'era5_{level}_daily.parquet'
        if not temp_file.exists():
            continue
        
        try:
            temp_df = pd.read_parquet(temp_file)
            temp_df['date'] = pd.to_datetime(temp_df['date'])
            
            level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
            
            # Overall temperature statistics
            temp_mean = temp_df['temp_mean'].mean()
            temp_sd = temp_df['temp_mean'].std()
            temp_min = temp_df['temp_mean'].min()
            temp_max = temp_df['temp_mean'].max()
            
            # Percentiles
            p1 = temp_df['temp_mean'].quantile(0.01)
            p5 = temp_df['temp_mean'].quantile(0.05)
            p25 = temp_df['temp_mean'].quantile(0.25)
            p50 = temp_df['temp_mean'].quantile(0.50)
            p75 = temp_df['temp_mean'].quantile(0.75)
            p95 = temp_df['temp_mean'].quantile(0.95)
            p99 = temp_df['temp_mean'].quantile(0.99)
            
            # Seasonal variation
            temp_df['month'] = temp_df['date'].dt.month
            monthly_means = temp_df.groupby('month')['temp_mean'].mean()
            
            rows.append({
                'Level': level_label,
                'Mean (SD)': f"{temp_mean:.1f} ({temp_sd:.1f})",
                'Min': f"{temp_min:.1f}",
                'P1': f"{p1:.1f}",
                'P5': f"{p5:.1f}",
                'P25': f"{p25:.1f}",
                'Median': f"{p50:.1f}",
                'P75': f"{p75:.1f}",
                'P95': f"{p95:.1f}",
                'P99': f"{p99:.1f}",
                'Max': f"{temp_max:.1f}",
                'Coldest Month': f"{monthly_means.idxmin()} ({monthly_means.min():.1f}°C)",
                'Warmest Month': f"{monthly_means.idxmax()} ({monthly_means.max():.1f}°C)",
            })
            
        except Exception as e:
            print(f"  Error loading {level}: {e}")
    
    if rows:
        df = pd.DataFrame(rows)
        save_table(df, 'tableD2_temperature_descriptive')
        return df
    return None

def create_table_descriptive_by_region():
    """Create regional descriptive statistics."""
    print("\n[Descriptive] Regional statistics...")
    
    all_rows = []
    
    for level in ['intermediate', 'immediate']:
        mort_file = PHASE0_RESULTS / f'mortality_{level}_daily.parquet'
        temp_file = PHASE0_RESULTS / f'era5_{level}_daily.parquet'
        
        if not mort_file.exists() or not temp_file.exists():
            continue
        
        try:
            mort_df = pd.read_parquet(mort_file)
            temp_df = pd.read_parquet(temp_file)
            
            death_col = get_death_column(mort_df)
            region_col = get_region_column(mort_df, level)
            temp_region_col = get_region_column(temp_df, level)
            
            level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
            
            # Aggregate by region
            mort_by_region = mort_df.groupby(region_col).agg({
                death_col: ['sum', 'mean', 'std']
            }).reset_index()
            mort_by_region.columns = ['region_code', 'total_deaths', 'daily_mean', 'daily_sd']
            
            temp_by_region = temp_df.groupby(temp_region_col).agg({
                'temp_mean': ['mean', 'std', 'min', 'max']
            }).reset_index()
            temp_by_region.columns = ['region_code', 'temp_mean', 'temp_sd', 'temp_min', 'temp_max']
            
            # Merge
            regional = mort_by_region.merge(temp_by_region, on='region_code', how='inner')
            
            # Top 10 by deaths
            top10 = regional.nlargest(10, 'total_deaths')
            
            for _, row in top10.iterrows():
                all_rows.append({
                    'Level': level_label,
                    'Region Code': int(row['region_code']),
                    'Total Deaths': format_number(row['total_deaths']),
                    'Daily Deaths Mean (SD)': f"{row['daily_mean']:.1f} ({row['daily_sd']:.1f})",
                    'Mean Temp (°C)': f"{row['temp_mean']:.1f}",
                    'Temp Range (°C)': f"{row['temp_min']:.1f} – {row['temp_max']:.1f}",
                })
            
        except Exception as e:
            print(f"  Error: {e}")
    
    if all_rows:
        df = pd.DataFrame(all_rows)
        save_table(df, 'tableD3_regional_descriptive')
        return df
    return None

def create_table_annual_summary():
    """Create annual summary statistics."""
    print("\n[Descriptive] Annual summary...")
    
    rows = []
    
    for level in ['intermediate', 'immediate']:
        mort_file = PHASE0_RESULTS / f'mortality_{level}_daily.parquet'
        temp_file = PHASE0_RESULTS / f'era5_{level}_daily.parquet'
        
        if not mort_file.exists() or not temp_file.exists():
            continue
        
        try:
            mort_df = pd.read_parquet(mort_file)
            temp_df = pd.read_parquet(temp_file)
            mort_df['date'] = pd.to_datetime(mort_df['date'])
            temp_df['date'] = pd.to_datetime(temp_df['date'])
            
            death_col = get_death_column(mort_df)
            
            mort_df['year'] = mort_df['date'].dt.year
            temp_df['year'] = temp_df['date'].dt.year
            
            level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
            
            # Annual aggregation
            annual_deaths = mort_df.groupby('year')[death_col].sum()
            annual_temp = temp_df.groupby('year')['temp_mean'].mean()
            
            # Days above/below thresholds
            p99 = temp_df['temp_mean'].quantile(0.99)
            p1 = temp_df['temp_mean'].quantile(0.01)
            
            for year in sorted(annual_deaths.index):
                year_temp = temp_df[temp_df['year'] == year]['temp_mean']
                hot_days = (year_temp > p99).sum()
                cold_days = (year_temp < p1).sum()
                
                rows.append({
                    'Level': level_label,
                    'Year': year,
                    'Total Deaths': format_number(annual_deaths.get(year, 0)),
                    'Mean Temp (°C)': f"{annual_temp.get(year, np.nan):.1f}",
                    'Extreme Hot Days (>P99)': hot_days,
                    'Extreme Cold Days (<P1)': cold_days,
                })
                
        except Exception as e:
            print(f"  Error: {e}")
    
    if rows:
        df = pd.DataFrame(rows)
        save_table(df, 'tableD4_annual_summary')
        return df
    return None

def create_table_covariates():
    """Create covariate/SES descriptive statistics."""
    print("\n[Descriptive] Covariates/SES...")
    
    rows = []
    
    for level in ['intermediate', 'immediate']:
        if level == 'intermediate':
            cov_file = PHASE0_RESULTS / 'regional_covariates.csv'
        else:
            cov_file = PHASE0_RESULTS / 'ses_immediate_covariates.csv'
        
        if not cov_file.exists():
            continue
        
        try:
            cov_df = pd.read_csv(cov_file)
            level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
            
            # Identify key covariates
            covariate_cols = {
                'pop_elderly': 'Elderly Population',
                'pop_total': 'Total Population',
                'hdi': 'Human Development Index',
                'ac_penetration': 'AC Penetration (%)',
                'urban_pct': 'Urban Population (%)',
                'mean_temp': 'Mean Temperature (°C)',
                'temp_range': 'Temperature Range (°C)',
            }
            
            for col, label in covariate_cols.items():
                if col in cov_df.columns:
                    values = cov_df[col].dropna()
                    if len(values) > 0:
                        rows.append({
                            'Level': level_label,
                            'Variable': label,
                            'N': len(values),
                            'Mean': f"{values.mean():.2f}",
                            'SD': f"{values.std():.2f}",
                            'Min': f"{values.min():.2f}",
                            'Median': f"{values.median():.2f}",
                            'Max': f"{values.max():.2f}",
                        })
                        
        except Exception as e:
            print(f"  Error: {e}")
    
    if rows:
        df = pd.DataFrame(rows)
        save_table(df, 'tableD5_covariates')
        return df
    return None

# =============================================================================
# TABLE 1: MAIN EFFECTS SUMMARY
# =============================================================================

def create_table1_main_effects():
    """Create main effects summary table."""
    print("\n[Table 1] Main effects summary...")
    
    rows = []
    
    for level in ['intermediate', 'immediate']:
        suffix = '' if level == 'intermediate' else '_immediate'
        level_label = 'Intermediate (133)' if level == 'intermediate' else 'Immediate (510)'
        
        # Try cause stratification first (has all-cause)
        cause_results = load_json(PHASE4_RESULTS / f'cause_stratification_v2_results{suffix}.json')
        
        if cause_results and 'pooled_results' in cause_results:
            all_cause = cause_results['pooled_results'].get('all_cause', {})
            
            for effect_type in ['heat', 'cold']:
                effect = all_cause.get(effect_type, {})
                if not effect:
                    continue
                    
                row = {
                    'Spatial Level': level_label,
                    'Effect Type': effect_type.title(),
                    'Percentile': 'P99 vs MMT' if effect_type == 'heat' else 'P1 vs MMT',
                    'RR': effect.get('rr', np.nan),
                    'RR Lower': effect.get('rr_lo', np.nan),
                    'RR Upper': effect.get('rr_hi', np.nan),
                    'I²': effect.get('I2', np.nan),
                    'τ²': effect.get('tau2', np.nan),
                    'N Regions': effect.get('k', np.nan),
                }
                
                row['RR (95% CI)'] = format_rr(row['RR'], row['RR Lower'], row['RR Upper'])
                row['I² (%)'] = f"{row['I²']:.1f}" if not np.isnan(row.get('I²', np.nan)) else "—"
                
                rows.append(row)
        else:
            # Fallback to pooled DLNM results
            pooled = load_json(PHASE1_RESULTS / f'dlnm_v2_{level}_pooled.json')
            
            if pooled is None:
                continue
            
            effects_map = {
                'heat': pooled.get('p99', {}),
                'cold': pooled.get('p1', {}),
            }
            
            for effect_type, effect in effects_map.items():
                if not effect:
                    continue
                    
                row = {
                    'Spatial Level': level_label,
                    'Effect Type': effect_type.title(),
                    'Percentile': 'P99 vs MMT' if effect_type == 'heat' else 'P1 vs MMT',
                    'RR': effect.get('pooled_rr', effect.get('rr', np.nan)),
                    'RR Lower': effect.get('pooled_rr_lower', effect.get('rr_lo', np.nan)),
                    'RR Upper': effect.get('pooled_rr_upper', effect.get('rr_hi', np.nan)),
                    'I²': effect.get('I2', np.nan),
                    'τ²': effect.get('tau2', np.nan),
                    'N Regions': effect.get('n_regions', effect.get('k', np.nan)),
                }
                
                row['RR (95% CI)'] = format_rr(row['RR'], row['RR Lower'], row['RR Upper'])
                row['I² (%)'] = f"{row['I²']:.1f}" if not np.isnan(row.get('I²', np.nan)) else "—"
                
                rows.append(row)
    
    if rows:
        df = pd.DataFrame(rows)
        df = df[['Spatial Level', 'Effect Type', 'Percentile', 'RR (95% CI)', 
                 'I² (%)', 'N Regions']]
        save_table(df, 'table1_main_effects')
        return df
    return None

# =============================================================================
# TABLE 2: ATTRIBUTABLE BURDEN
# =============================================================================

def create_table2_attributable_burden():
    """Create attributable burden table."""
    print("\n[Table 2] Attributable burden...")
    
    # Load national summary
    burden = load_json(PHASE1_RESULTS / 'burden_v2_national_summary.json')
    
    if burden is None:
        print("  Skipping: No burden data")
        return None
    
    rows = []
    
    # National summary rows - one per spatial level
    for level in ['intermediate', 'immediate']:
        level_data = burden.get(level, {})
        if not level_data:
            continue
            
        level_label = 'Intermediate (133 regions)' if level == 'intermediate' else 'Immediate (510 regions)'
        
        # Handle different key naming conventions
        heat_deaths = level_data.get('total_heat_an', level_data.get('total_heat_attributable_deaths', 0))
        cold_deaths = level_data.get('total_cold_an', level_data.get('total_cold_attributable_deaths', 0))
        total_deaths = level_data.get('total_deaths', heat_deaths + cold_deaths)
        heat_af = level_data.get('total_heat_af_pct', level_data.get('heat_af_percent', 0))
        cold_af = level_data.get('total_cold_af_pct', level_data.get('cold_af_percent', 0))
        
        rows.append({
            'Region': f'Brazil ({level_label})',
            'Heat Deaths': format_number(heat_deaths),
            'Cold Deaths': format_number(cold_deaths),
            'Total Deaths': format_number(total_deaths),
            'Heat AF (%)': f"{heat_af:.2f}",
            'Cold AF (%)': f"{cold_af:.2f}",
        })
    
    # Regional burden if available
    for level in ['intermediate', 'immediate']:
        regional_file = PHASE1_RESULTS / f'burden_v2_{level}_regions.csv'
        if regional_file.exists():
            regional_df = pd.read_csv(regional_file)
            
            # Find heat deaths column
            heat_col = None
            for col in ['total_heat_an', 'heat_deaths', 'heat_attributable_deaths']:
                if col in regional_df.columns:
                    heat_col = col
                    break
            
            if heat_col:
                top5 = regional_df.nlargest(5, heat_col)
                
                for _, r in top5.iterrows():
                    region_name = r.get('region_name', r.get('region_code', 'Unknown'))
                    
                    # Handle different column names
                    cold_col = 'total_cold_an' if 'total_cold_an' in regional_df.columns else 'cold_deaths'
                    total_col = 'total_deaths' if 'total_deaths' in regional_df.columns else None
                    heat_af_col = 'total_heat_af_pct' if 'total_heat_af_pct' in regional_df.columns else 'heat_af'
                    cold_af_col = 'total_cold_af_pct' if 'total_cold_af_pct' in regional_df.columns else 'cold_af'
                    
                    heat_val = r.get(heat_col, 0)
                    cold_val = r.get(cold_col, 0)
                    total_val = r.get(total_col, heat_val + cold_val) if total_col else heat_val + cold_val
                    
                    rows.append({
                        'Region': f"{region_name} ({level})",
                        'Heat Deaths': format_number(heat_val),
                        'Cold Deaths': format_number(cold_val),
                        'Total Deaths': format_number(total_val),
                        'Heat AF (%)': f"{r.get(heat_af_col, 0):.2f}",
                        'Cold AF (%)': f"{r.get(cold_af_col, 0):.2f}",
                    })
    
    if rows:
        df = pd.DataFrame(rows)
        save_table(df, 'table2_attributable_burden')
        return df
    return None

# =============================================================================
# TABLE 3: SENSITIVITY ANALYSES
# =============================================================================

def create_table3_sensitivity():
    """Create sensitivity analyses summary table."""
    print("\n[Table 3] Sensitivity analyses...")
    
    rows = []
    
    for level in ['intermediate', 'immediate']:
        suffix = '' if level == 'intermediate' else '_immediate'
        level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
        
        results = load_json(PHASE2_RESULTS / f'sensitivity_analyses_v2{suffix}.json')
        
        if results is None:
            continue
        
        for analysis_name, analysis_data in results.items():
            if not isinstance(analysis_data, dict):
                continue
            
            pooled = analysis_data.get('pooled', {})
            
            # Heat effect
            heat = pooled.get('p99', pooled.get('heat', {}))
            if isinstance(heat, dict) and 'pooled_rr' in heat:
                rows.append({
                    'Spatial Level': level_label,
                    'Analysis': analysis_name.replace('_', ' ').title(),
                    'Effect': 'Heat (P99)',
                    'RR': heat.get('pooled_rr', np.nan),
                    'RR Lower': heat.get('pooled_rr_lower', np.nan),
                    'RR Upper': heat.get('pooled_rr_upper', np.nan),
                    'I²': heat.get('I2', np.nan),
                    'N Regions': heat.get('n_regions', np.nan),
                })
            
            # Cold effect
            cold = pooled.get('p1', pooled.get('cold', {}))
            if isinstance(cold, dict) and 'pooled_rr' in cold:
                rows.append({
                    'Spatial Level': level_label,
                    'Analysis': analysis_name.replace('_', ' ').title(),
                    'Effect': 'Cold (P1)',
                    'RR': cold.get('pooled_rr', np.nan),
                    'RR Lower': cold.get('pooled_rr_lower', np.nan),
                    'RR Upper': cold.get('pooled_rr_upper', np.nan),
                    'I²': cold.get('I2', np.nan),
                    'N Regions': cold.get('n_regions', np.nan),
                })
    
    if rows:
        df = pd.DataFrame(rows)
        df['RR (95% CI)'] = df.apply(lambda r: format_rr(r['RR'], r['RR Lower'], r['RR Upper']), axis=1)
        df['I² (%)'] = df['I²'].apply(lambda x: f"{x:.1f}" if pd.notna(x) else "—")
        
        df = df[['Spatial Level', 'Analysis', 'Effect', 'RR (95% CI)', 'I² (%)', 'N Regions']]
        save_table(df, 'table3_sensitivity_analyses')
        return df
    return None

# =============================================================================
# TABLE 4: CONFOUNDER ADJUSTMENT
# =============================================================================

def create_table4_confounder():
    """Create confounder adjustment comparison table for both levels."""
    print("\n[Table 4] Confounder adjustment...")
    
    rows = []
    
    analysis_labels = {
        'base': 'Base Model (Dry-bulb)',
        'apparent_temp': 'Apparent Temperature',
        'pollution_adjusted': 'Pollution Adjusted (PM2.5 + O3)',
        'influenza_adjusted': 'Influenza Adjusted',
        'fully_adjusted': 'Fully Adjusted (All)'
    }
    
    for level in ['intermediate', 'immediate']:
        suffix = '' if level == 'intermediate' else '_immediate'
        level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
        
        results = load_json(PHASE3_RESULTS / f'supplementary_analyses_v2{suffix}.json')
        
        if results is None:
            continue
        
        for analysis_key, label in analysis_labels.items():
            if analysis_key not in results:
                continue
            
            analysis_data = results[analysis_key]
            if not isinstance(analysis_data, dict):
                continue
            
            pooled = analysis_data.get('pooled', {})
            
            # Use P97.5 for heat and P2.5 for cold (better pooled estimates than P99/P1)
            for pct, pct_label in [('p97.5', 'Heat (P97.5)'), ('p2.5', 'Cold (P2.5)')]:
                if pct in pooled:
                    effect = pooled[pct]
                    rows.append({
                        'Spatial Level': level_label,
                        'Model': label,
                        'Effect': pct_label,
                        'RR': effect.get('pooled_rr', np.nan),
                        'RR Lower': effect.get('pooled_rr_lower', np.nan),
                        'RR Upper': effect.get('pooled_rr_upper', np.nan),
                        'I²': effect.get('I2', np.nan),
                        'N Regions': effect.get('n_regions', np.nan),
                    })
    
    if rows:
        df = pd.DataFrame(rows)
        df['RR (95% CI)'] = df.apply(lambda r: format_rr(r['RR'], r['RR Lower'], r['RR Upper']), axis=1)
        df['I² (%)'] = df['I²'].apply(lambda x: f"{x:.1f}" if pd.notna(x) else "—")
        
        df = df[['Spatial Level', 'Model', 'Effect', 'RR (95% CI)', 'I² (%)', 'N Regions']]
        save_table(df, 'table4_confounder_adjustment')
        return df
    return None

# =============================================================================
# TABLE 5: STRATIFIED ANALYSES
# =============================================================================

def create_table5_stratification():
    """Create stratified analyses table (age, sex, cause)."""
    print("\n[Table 5] Stratified analyses...")
    
    rows = []
    
    # Age stratification
    for level in ['intermediate', 'immediate']:
        suffix = '' if level == 'intermediate' else '_immediate'
        level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
        
        results = load_json(PHASE4_RESULTS / f'age_stratification_v2_results{suffix}.json')
        if results and 'pooled_results' in results:
            for age_group, data in results['pooled_results'].items():
                for effect_type in ['heat', 'cold']:
                    if effect_type in data:
                        effect = data[effect_type]
                        rows.append({
                            'Stratification': 'Age Group',
                            'Category': age_group,
                            'Spatial Level': level_label,
                            'Effect': effect_type.title(),
                            'RR': effect.get('rr', np.nan),
                            'RR Lower': effect.get('rr_lo', np.nan),
                            'RR Upper': effect.get('rr_hi', np.nan),
                            'I²': effect.get('I2', np.nan),
                        })
    
    # Sex stratification
    for level in ['intermediate', 'immediate']:
        suffix = '' if level == 'intermediate' else '_immediate'
        level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
        
        results = load_json(PHASE4_RESULTS / f'sex_stratification_v2_results{suffix}.json')
        if results and 'pooled_results' in results:
            for sex, data in results['pooled_results'].items():
                for effect_type in ['heat', 'cold']:
                    if effect_type in data:
                        effect = data[effect_type]
                        rows.append({
                            'Stratification': 'Sex',
                            'Category': sex,
                            'Spatial Level': level_label,
                            'Effect': effect_type.title(),
                            'RR': effect.get('rr', np.nan),
                            'RR Lower': effect.get('rr_lo', np.nan),
                            'RR Upper': effect.get('rr_hi', np.nan),
                            'I²': effect.get('I2', np.nan),
                        })
    
    # Cause stratification - both levels
    for level in ['intermediate', 'immediate']:
        suffix = '' if level == 'intermediate' else '_immediate'
        level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
        
        results = load_json(PHASE4_RESULTS / f'cause_stratification_v2_results{suffix}.json')
        if results and 'pooled_results' in results:
            for cause, data in results['pooled_results'].items():
                for effect_type in ['heat', 'cold']:
                    if effect_type in data:
                        effect = data[effect_type]
                        rows.append({
                            'Stratification': 'Cause of Death',
                            'Category': cause.replace('_', ' ').title(),
                            'Spatial Level': level_label,
                            'Effect': effect_type.title(),
                            'RR': effect.get('rr', np.nan),
                            'RR Lower': effect.get('rr_lo', np.nan),
                            'RR Upper': effect.get('rr_hi', np.nan),
                            'I²': effect.get('I2', np.nan),
                        })
    
    if rows:
        df = pd.DataFrame(rows)
        df['RR (95% CI)'] = df.apply(lambda r: format_rr(r['RR'], r['RR Lower'], r['RR Upper']), axis=1)
        df['I² (%)'] = df['I²'].apply(lambda x: f"{x:.1f}" if pd.notna(x) else "—")
        
        df = df[['Stratification', 'Category', 'Spatial Level', 'Effect', 'RR (95% CI)', 'I² (%)']]
        save_table(df, 'table5_stratified_analyses')
        return df
    return None

# =============================================================================
# TABLE 6: META-REGRESSION
# =============================================================================

def create_table6_meta_regression():
    """Create meta-regression coefficients table."""
    print("\n[Table 6] Meta-regression...")
    
    rows = []
    
    for level in ['intermediate', 'immediate']:
        suffix = '' if level == 'intermediate' else '_immediate'
        level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
        
        results = load_json(PHASE4_RESULTS / f'meta_regression_v2_results{suffix}.json')
        
        if results is None:
            continue
        
        for effect_type in ['heat', 'cold']:
            if effect_type not in results:
                continue
            
            effect_data = results[effect_type]
            
            if 'moderators' in effect_data:
                for mod_name, mod_data in effect_data['moderators'].items():
                    if isinstance(mod_data, dict):
                        rows.append({
                            'Spatial Level': level_label,
                            'Effect': effect_type.title(),
                            'Moderator': mod_name.replace('_', ' ').title(),
                            'Coefficient': mod_data.get('coefficient', mod_data.get('coef', np.nan)),
                            'SE': mod_data.get('se', np.nan),
                            'p-value': mod_data.get('p', mod_data.get('pvalue', np.nan)),
                            'Significant': 'Yes' if mod_data.get('p', mod_data.get('pvalue', 1)) < 0.05 else 'No',
                        })
    
    if rows:
        df = pd.DataFrame(rows)
        df['Coefficient'] = df['Coefficient'].apply(lambda x: f"{x:.4f}" if pd.notna(x) else "—")
        df['SE'] = df['SE'].apply(lambda x: f"{x:.4f}" if pd.notna(x) else "—")
        df['p-value'] = df['p-value'].apply(lambda x: f"{x:.4f}" if pd.notna(x) else "—")
        
        save_table(df, 'table6_meta_regression')
        return df
    return None

# =============================================================================
# TABLE S1: REGIONAL EFFECTS SUMMARY
# =============================================================================

def create_tableS1_regional_summary():
    """Create regional effects summary table."""
    print("\n[Table S1] Regional effects summary...")
    
    rows = []
    
    for level in ['intermediate', 'immediate']:
        summary_file = PHASE1_RESULTS / f'dlnm_v2_{level}_summary.csv'
        if summary_file.exists():
            df = pd.read_csv(summary_file)
            
            # Use correct column names (rr_p99 for heat, rr_p1 for cold)
            heat_col = 'rr_p99' if 'rr_p99' in df.columns else 'heat_rr'
            cold_col = 'rr_p1' if 'rr_p1' in df.columns else 'cold_rr'
            
            # Summary statistics
            rows.append({
                'Level': level.title(),
                'N Regions': len(df),
                'Heat RR Mean': df[heat_col].mean() if heat_col in df.columns else np.nan,
                'Heat RR Median': df[heat_col].median() if heat_col in df.columns else np.nan,
                'Heat RR Min': df[heat_col].min() if heat_col in df.columns else np.nan,
                'Heat RR Max': df[heat_col].max() if heat_col in df.columns else np.nan,
                'Cold RR Mean': df[cold_col].mean() if cold_col in df.columns else np.nan,
                'Cold RR Median': df[cold_col].median() if cold_col in df.columns else np.nan,
                'Cold RR Min': df[cold_col].min() if cold_col in df.columns else np.nan,
                'Cold RR Max': df[cold_col].max() if cold_col in df.columns else np.nan,
            })
    
    if rows:
        df = pd.DataFrame(rows)
        for col in df.columns:
            if 'RR' in col:
                df[col] = df[col].apply(lambda x: f"{x:.3f}" if pd.notna(x) else "—")
        save_table(df, 'tableS1_regional_summary')
        return df
    return None

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("\nCreating descriptive tables...")
    create_table_descriptive_mortality()
    create_table_descriptive_temperature()
    create_table_descriptive_by_region()
    create_table_annual_summary()
    create_table_covariates()
    
    print("\nCreating analytical tables...")
    create_table1_main_effects()
    create_table2_attributable_burden()
    create_table3_sensitivity()
    create_table4_confounder()
    create_table5_stratification()
    create_table6_meta_regression()
    create_tableS1_regional_summary()
    
    print("\n" + "="*70)
    print(f"TABLE GENERATION COMPLETE")
    print(f"Output directory: {OUTPUT_DIR}")
    print("="*70)
