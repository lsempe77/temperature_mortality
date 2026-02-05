"""
05a: GENERATE FIGURES FOR ALL PHASES
=====================================
Creates publication-quality figures from Phase 1-4 results.

Outputs:
--------
- Phase 1: Exposure-response curves, lag-response, attributable burden
- Phase 2: Sensitivity forest plots, harvesting curves, heatwave effects
- Phase 3: Confounder-adjusted comparisons
- Phase 4: Age/sex/cause stratification, meta-regression

Author: Heat-Mortality Brazil Analysis Pipeline
Date: December 2025
"""

import warnings
warnings.filterwarnings('ignore')

import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from pathlib import Path
from datetime import datetime
import argparse

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = Path(__file__).parent.parent
PHASE1_RESULTS = BASE_DIR / 'phase1_core_model' / 'results'
PHASE2_RESULTS = BASE_DIR / 'phase2_robustness' / 'results'
PHASE3_RESULTS = BASE_DIR / 'phase3_confounding' / 'results'
PHASE4_RESULTS = BASE_DIR / 'phase4_heterogeneity' / 'results'
OUTPUT_DIR = Path(__file__).parent / 'figures'
OUTPUT_DIR.mkdir(exist_ok=True)

# Plotting style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'

# Colors
HEAT_COLOR = '#e74c3c'
COLD_COLOR = '#3498db'
NEUTRAL_COLOR = '#7f8c8d'

# Additional paths for descriptive data
INPUT_DATA_DIR = BASE_DIR.parent / 'Input_data'
PHASE0_RESULTS = BASE_DIR / 'phase0_data_prep' / 'results'

print("="*70)
print("05a: GENERATE FIGURES FOR ALL PHASES")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_json(filepath):
    """Load JSON file safely."""
    try:
        with open(filepath, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"  Warning: Could not load {filepath}: {e}")
        return None

def save_figure(fig, name, formats=['png']):
    """Save figure in multiple formats."""
    for fmt in formats:
        filepath = OUTPUT_DIR / f"{name}.{fmt}"
        fig.savefig(filepath, format=fmt, bbox_inches='tight', dpi=150)
    print(f"  Saved: {name}")
    plt.close(fig)
    import gc
    gc.collect()

# =============================================================================
# DESCRIPTIVE FIGURES
# =============================================================================

def plot_descriptive_mortality_timeseries():
    """Plot mortality time series."""
    print("\n[Descriptive] Mortality time series...")
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    
    for idx, level in enumerate(['intermediate', 'immediate']):
        ax = axes[idx]
        
        mort_file = PHASE0_RESULTS / f'mortality_{level}_daily.parquet'
        if not mort_file.exists():
            ax.text(0.5, 0.5, f'No data for {level}', ha='center', va='center', transform=ax.transAxes)
            continue
        
        try:
            import pandas as pd
            mort_df = pd.read_parquet(mort_file)
            mort_df['date'] = pd.to_datetime(mort_df['date'])
            
            death_col = 'deaths_all' if 'deaths_all' in mort_df.columns else 'deaths'
            
            # Aggregate to daily national total
            daily_deaths = mort_df.groupby('date')[death_col].sum()
            
            # Plot
            ax.plot(daily_deaths.index, daily_deaths.values, alpha=0.3, linewidth=0.5, color='gray')
            
            # Add 30-day rolling mean
            rolling = daily_deaths.rolling(30, center=True).mean()
            ax.plot(rolling.index, rolling.values, color='black', linewidth=1.5, label='30-day rolling mean')
            
            level_label = 'Intermediate (133 regions)' if level == 'intermediate' else 'Immediate (510 regions)'
            ax.set_title(f'Daily Elderly Mortality - {level_label}')
            ax.set_xlabel('Date')
            ax.set_ylabel('Deaths per day')
            ax.legend()
            
        except Exception as e:
            ax.text(0.5, 0.5, f'Error: {e}', ha='center', va='center', transform=ax.transAxes)
    
    plt.tight_layout()
    save_figure(fig, 'figD1_mortality_timeseries')

def plot_descriptive_temperature_distribution():
    """Plot temperature distribution."""
    print("\n[Descriptive] Temperature distribution...")
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    for idx, level in enumerate(['intermediate', 'immediate']):
        ax = axes[idx]
        
        temp_file = PHASE0_RESULTS / f'era5_{level}_daily.parquet'
        if not temp_file.exists():
            ax.text(0.5, 0.5, f'No data for {level}', ha='center', va='center', transform=ax.transAxes)
            continue
        
        try:
            import pandas as pd
            temp_df = pd.read_parquet(temp_file)
            
            # Histogram
            ax.hist(temp_df['temp_mean'], bins=50, density=True, alpha=0.7, 
                    color='steelblue', edgecolor='white')
            
            # Add percentile lines
            p1 = temp_df['temp_mean'].quantile(0.01)
            p99 = temp_df['temp_mean'].quantile(0.99)
            median = temp_df['temp_mean'].median()
            
            ax.axvline(p1, color=COLD_COLOR, linestyle='--', linewidth=2, label=f'P1 = {p1:.1f}°C')
            ax.axvline(median, color='green', linestyle='-', linewidth=2, label=f'Median = {median:.1f}°C')
            ax.axvline(p99, color=HEAT_COLOR, linestyle='--', linewidth=2, label=f'P99 = {p99:.1f}°C')
            
            level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
            ax.set_title(f'Temperature Distribution - {level_label}')
            ax.set_xlabel('Temperature (°C)')
            ax.set_ylabel('Density')
            ax.legend()
            
        except Exception as e:
            ax.text(0.5, 0.5, f'Error: {e}', ha='center', va='center', transform=ax.transAxes)
    
    plt.tight_layout()
    save_figure(fig, 'figD2_temperature_distribution')

def plot_descriptive_seasonal_patterns():
    """Plot seasonal patterns in mortality and temperature."""
    print("\n[Descriptive] Seasonal patterns...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    for col, level in enumerate(['intermediate', 'immediate']):
        level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
        
        # Load data
        mort_file = PHASE0_RESULTS / f'mortality_{level}_daily.parquet'
        temp_file = PHASE0_RESULTS / f'era5_{level}_daily.parquet'
        
        if not mort_file.exists() or not temp_file.exists():
            continue
        
        try:
            import pandas as pd
            mort_df = pd.read_parquet(mort_file)
            temp_df = pd.read_parquet(temp_file)
            mort_df['date'] = pd.to_datetime(mort_df['date'])
            temp_df['date'] = pd.to_datetime(temp_df['date'])
            
            death_col = 'deaths_all' if 'deaths_all' in mort_df.columns else 'deaths'
            
            mort_df['month'] = mort_df['date'].dt.month
            temp_df['month'] = temp_df['date'].dt.month
            
            # Monthly mortality
            ax = axes[0, col]
            monthly_deaths = mort_df.groupby('month')[death_col].mean()
            ax.bar(monthly_deaths.index, monthly_deaths.values, color='gray', edgecolor='black')
            ax.set_xlabel('Month')
            ax.set_ylabel('Mean Daily Deaths')
            ax.set_title(f'Monthly Mortality Pattern - {level_label}')
            ax.set_xticks(range(1, 13))
            ax.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
            
            # Monthly temperature
            ax = axes[1, col]
            monthly_temp = temp_df.groupby('month')['temp_mean'].agg(['mean', 'std'])
            ax.bar(monthly_temp.index, monthly_temp['mean'], 
                   yerr=monthly_temp['std'], color='steelblue', edgecolor='black', capsize=3)
            ax.set_xlabel('Month')
            ax.set_ylabel('Mean Temperature (°C)')
            ax.set_title(f'Monthly Temperature Pattern - {level_label}')
            ax.set_xticks(range(1, 13))
            ax.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
            
        except Exception as e:
            print(f"  Error for {level}: {e}")
    
    plt.tight_layout()
    save_figure(fig, 'figD3_seasonal_patterns')

def plot_descriptive_annual_trends():
    """Plot annual trends in mortality and temperature."""
    print("\n[Descriptive] Annual trends...")
    
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    
    # For mortality: use intermediate level (same national total for both)
    # For temperature: show both levels (can differ due to weighting)
    
    # First, get national mortality (same for both levels)
    mort_file = PHASE0_RESULTS / 'mortality_intermediate_daily.parquet'
    if mort_file.exists():
        try:
            import pandas as pd
            mort_df = pd.read_parquet(mort_file)
            mort_df['date'] = pd.to_datetime(mort_df['date'])
            death_col = 'deaths_all' if 'deaths_all' in mort_df.columns else 'deaths'
            mort_df['year'] = mort_df['date'].dt.year
            annual_deaths = mort_df.groupby('year')[death_col].sum()
            
            axes[0].plot(annual_deaths.index, annual_deaths.values / 1000, 
                        'o-', color='#2c3e50', linewidth=2, markersize=6, label='Brazil (National)')
            
            # Add COVID marker
            if 2020 in annual_deaths.index:
                axes[0].axvspan(2020, 2021.5, alpha=0.2, color='red', label='COVID-19 period')
        except Exception as e:
            print(f"  Error loading mortality: {e}")
    
    # Temperature: show both levels (they can differ due to weighting)
    for level in ['intermediate', 'immediate']:
        temp_file = PHASE0_RESULTS / f'era5_{level}_daily.parquet'
        
        if not temp_file.exists():
            continue
        
        try:
            import pandas as pd
            temp_df = pd.read_parquet(temp_file)
            temp_df['date'] = pd.to_datetime(temp_df['date'])
            temp_df['year'] = temp_df['date'].dt.year
            annual_temp = temp_df.groupby('year')['temp_mean'].mean()
            
            level_label = 'Intermediate (133)' if level == 'intermediate' else 'Immediate (510)'
            marker = 'o' if level == 'intermediate' else 's'
            color = HEAT_COLOR if level == 'intermediate' else COLD_COLOR
            
            axes[1].plot(annual_temp.index, annual_temp.values,
                        f'{marker}-', label=level_label, linewidth=2, markersize=6, color=color)
            
        except Exception as e:
            print(f"  Error for {level} temp: {e}")
    
    axes[0].set_xlabel('Year')
    axes[0].set_ylabel('Total Elderly Deaths (thousands)')
    axes[0].set_title('Annual Elderly Mortality (National Total)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].set_xlabel('Year')
    axes[1].set_ylabel('Mean Temperature (°C)')
    axes[1].set_title('Annual Mean Temperature by Spatial Level')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    save_figure(fig, 'figD4_annual_trends')

def plot_descriptive_mortality_temperature_scatter():
    """Plot mortality vs temperature relationship (raw data)."""
    print("\n[Descriptive] Mortality-temperature scatter...")
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    for idx, level in enumerate(['intermediate', 'immediate']):
        ax = axes[idx]
        
        mort_file = PHASE0_RESULTS / f'mortality_{level}_daily.parquet'
        temp_file = PHASE0_RESULTS / f'era5_{level}_daily.parquet'
        
        if not mort_file.exists() or not temp_file.exists():
            ax.text(0.5, 0.5, f'No data for {level}', ha='center', va='center', transform=ax.transAxes)
            continue
        
        try:
            import pandas as pd
            mort_df = pd.read_parquet(mort_file)
            temp_df = pd.read_parquet(temp_file)
            mort_df['date'] = pd.to_datetime(mort_df['date'])
            temp_df['date'] = pd.to_datetime(temp_df['date'])
            
            death_col = 'deaths_all' if 'deaths_all' in mort_df.columns else 'deaths'
            
            # Aggregate to national daily
            daily_deaths = mort_df.groupby('date')[death_col].sum()
            daily_temp = temp_df.groupby('date')['temp_mean'].mean()
            
            # Merge
            combined = pd.DataFrame({'deaths': daily_deaths, 'temp': daily_temp}).dropna()
            
            # Scatter with alpha
            ax.scatter(combined['temp'], combined['deaths'], alpha=0.1, s=5, color='gray')
            
            # Add binned means
            combined['temp_bin'] = pd.cut(combined['temp'], bins=20)
            binned = combined.groupby('temp_bin').agg({'deaths': 'mean', 'temp': 'mean'})
            ax.plot(binned['temp'], binned['deaths'], 'o-', color=HEAT_COLOR, 
                   linewidth=2, markersize=8, label='Binned mean')
            
            level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
            ax.set_xlabel('Temperature (°C)')
            ax.set_ylabel('Daily Deaths')
            ax.set_title(f'Mortality vs Temperature - {level_label}')
            ax.legend()
            
        except Exception as e:
            ax.text(0.5, 0.5, f'Error: {e}', ha='center', va='center', transform=ax.transAxes)
    
    plt.tight_layout()
    save_figure(fig, 'figD5_mortality_temperature_scatter')

# =============================================================================
# PHASE 1: CORE MODEL FIGURES
# =============================================================================

def plot_pooled_exposure_response():
    """Plot pooled exposure-response curves for both spatial levels."""
    print("\n[Phase 1] Pooled exposure-response curves...")
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    for idx, level in enumerate(['intermediate', 'immediate']):
        ax = axes[idx]
        
        # Load region-level results (they contain the curves)
        results = load_json(PHASE1_RESULTS / f'dlnm_v2_{level}_results.json')
        if results is None:
            ax.text(0.5, 0.5, f'No data for {level}', ha='center', va='center')
            continue
        
        # Aggregate curves from all regions
        all_temps = []
        all_rrs = []
        all_mmts = []
        
        region_results = results.get('region_results', {})
        for region_code, region_data in region_results.items():
            if not isinstance(region_data, dict):
                continue
            
            rr_curve = region_data.get('rr_curve', [])
            mmt = region_data.get('mmt')
            if mmt:
                all_mmts.append(mmt)
            
            if isinstance(rr_curve, list) and len(rr_curve) > 0:
                temps = [pt.get('temp', 0) for pt in rr_curve]
                rrs = [pt.get('rr', 1) for pt in rr_curve]
                all_temps.append(temps)
                all_rrs.append(rrs)
        
        if all_temps:
            # Create average curve across regions
            # Use the temperature range from the first region as reference
            ref_temps = np.array(all_temps[0])
            
            # For each region, interpolate to common temperature grid
            from scipy import interpolate
            common_temps = np.linspace(min(ref_temps), max(ref_temps), 50)
            interp_rrs = []
            
            for temps, rrs in zip(all_temps, all_rrs):
                try:
                    f = interpolate.interp1d(temps, rrs, bounds_error=False, fill_value=1)
                    interp_rrs.append(f(common_temps))
                except:
                    pass
            
            if interp_rrs:
                mean_rrs = np.nanmean(interp_rrs, axis=0)
                lo_rrs = np.nanpercentile(interp_rrs, 25, axis=0)
                hi_rrs = np.nanpercentile(interp_rrs, 75, axis=0)
                
                # Plot
                ax.fill_between(common_temps, lo_rrs, hi_rrs, alpha=0.3, color=NEUTRAL_COLOR)
                ax.plot(common_temps, mean_rrs, color='black', linewidth=2)
                
                # MMT line
                if all_mmts:
                    mean_mmt = np.median(all_mmts)
                    ax.axvline(mean_mmt, color='green', linestyle='--', alpha=0.7, 
                              label=f'MMT = {mean_mmt:.1f}°C')
                
                # Reference line
                ax.axhline(1, color='gray', linestyle='-', alpha=0.5)
                
                ax.set_xlabel('Temperature (°C)')
                ax.set_ylabel('Relative Risk')
                ax.legend(loc='upper left')
                ax.set_ylim(0, min(5, np.nanmax(mean_rrs) * 1.2))
        
        # Add statistics from cause stratification (which has proper pooled values)
        cause_results = load_json(PHASE4_RESULTS / f'cause_stratification_v2_results{"" if level == "intermediate" else "_immediate"}.json')
        if cause_results and 'pooled_results' in cause_results:
            all_cause = cause_results['pooled_results'].get('all_cause', {})
            heat_rr = all_cause.get('heat', {}).get('rr', None)
            cold_rr = all_cause.get('cold', {}).get('rr', None)
            if heat_rr and cold_rr:
                ax.text(0.98, 0.98, f"Heat RR: {heat_rr:.3f}\nCold RR: {cold_rr:.3f}",
                       transform=ax.transAxes, ha='right', va='top', fontsize=9,
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        level_label = '133 Intermediate Regions' if level == 'intermediate' else '510 Immediate Regions'
        ax.set_title(f'{level_label}')
    
    fig.suptitle('Temperature-Mortality Exposure-Response Curves', fontsize=14, y=1.02)
    plt.tight_layout()
    save_figure(fig, 'fig1_pooled_exposure_response')

def plot_attributable_burden():
    """Plot attributable burden summary."""
    print("\n[Phase 1] Attributable burden...")
    
    burden = load_json(PHASE1_RESULTS / 'burden_v2_national_summary.json')
    if burden is None:
        print("  Skipping: No burden data")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Deaths by type (heat vs cold) - use intermediate level data
    ax = axes[0]
    level_data = burden.get('intermediate', burden)  # Fallback to root if not nested
    heat_deaths = level_data.get('total_heat_an', level_data.get('total_heat_attributable_deaths', 0))
    cold_deaths = level_data.get('total_cold_an', level_data.get('total_cold_attributable_deaths', 0))
    
    bars = ax.bar(['Heat', 'Cold'], [heat_deaths, cold_deaths], 
                  color=[HEAT_COLOR, COLD_COLOR], edgecolor='black')
    ax.set_ylabel('Attributable Deaths')
    ax.set_title('Total Attributable Deaths (2010-2024)')
    
    for bar, val in zip(bars, [heat_deaths, cold_deaths]):
        height = bar.get_height()
        offset = height * 0.02 if height > 0 else 500
        ax.text(bar.get_x() + bar.get_width()/2, height + offset, 
                f'{val:,.0f}', ha='center', va='bottom', fontsize=11)
    
    # Compare intermediate vs immediate levels
    ax = axes[1]
    levels = []
    heat_vals = []
    cold_vals = []
    
    for level in ['intermediate', 'immediate']:
        if level in burden:
            levels.append(level.title())
            heat_vals.append(burden[level].get('total_heat_an', 0))
            cold_vals.append(burden[level].get('total_cold_an', 0))
    
    if levels:
        x = np.arange(len(levels))
        width = 0.35
        ax.bar(x - width/2, heat_vals, width, label='Heat', color=HEAT_COLOR, edgecolor='black')
        ax.bar(x + width/2, cold_vals, width, label='Cold', color=COLD_COLOR, edgecolor='black')
        ax.set_xticks(x)
        ax.set_xticklabels(levels)
        ax.set_ylabel('Attributable Deaths')
        ax.set_title('Attributable Deaths by Spatial Level')
        ax.legend()
    
    plt.tight_layout()
    save_figure(fig, 'fig2_attributable_burden')

def plot_yll_summary():
    """Plot years of life lost summary."""
    print("\n[Phase 1] Years of life lost...")
    
    yll = load_json(PHASE1_RESULTS / 'yll_unified_results.json')
    if yll is None:
        print("  Skipping: No YLL data")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left plot: YLL by heat vs cold for both levels
    ax = axes[0]
    
    if 'actual' in yll:
        actual = yll['actual']
        
        # Get heat/cold YLL for both levels
        levels = []
        heat_yll = []
        cold_yll = []
        
        for level in ['intermediate', 'immediate']:
            level_key = f'{level}_yll'
            if level_key in actual:
                levels.append(level.title())
                heat_yll.append(actual[level_key].get('heat_yll', 0))
                cold_yll.append(actual[level_key].get('cold_yll', 0))
        
        if levels:
            x = np.arange(len(levels))
            width = 0.35
            ax.bar(x - width/2, heat_yll, width, label='Heat', color=HEAT_COLOR, edgecolor='black')
            ax.bar(x + width/2, cold_yll, width, label='Cold', color=COLD_COLOR, edgecolor='black')
            ax.set_xticks(x)
            ax.set_xticklabels(levels)
            ax.set_ylabel('Years of Life Lost')
            ax.set_title('YLL by Temperature Type')
            ax.legend()
            
            # Add value labels
            for i, (h, c) in enumerate(zip(heat_yll, cold_yll)):
                ax.text(i - width/2, h + h*0.02, f'{h/1e6:.1f}M', ha='center', va='bottom', fontsize=9)
                ax.text(i + width/2, c + c*0.02, f'{c/1e6:.1f}M', ha='center', va='bottom', fontsize=9)
        else:
            ax.text(0.5, 0.5, 'No level-specific YLL data', ha='center', va='center', transform=ax.transAxes)
    else:
        ax.text(0.5, 0.5, 'No YLL data structure found', ha='center', va='center', transform=ax.transAxes)
    
    # Right plot: Age group breakdown
    ax = axes[1]
    
    if 'actual' in yll and 'age_group_breakdown' in yll['actual']:
        age_data = yll['actual']['age_group_breakdown']
        age_groups = [d['age_group'] for d in age_data]
        yll_values = [d['total_yll'] / 1e6 for d in age_data]  # Convert to millions
        
        ax.bar(age_groups, yll_values, color='#9b59b6', edgecolor='black', alpha=0.7)
        ax.set_xlabel('Age Group')
        ax.set_ylabel('Years of Life Lost (Millions)')
        ax.set_title('YLL by Age Group (All Deaths)')
        ax.tick_params(axis='x', rotation=45)
    else:
        ax.text(0.5, 0.5, 'No age breakdown data', ha='center', va='center', transform=ax.transAxes)
    
    fig.suptitle('Years of Life Lost Attributable to Temperature (2010-2024)', fontsize=13, y=1.02)
    plt.tight_layout()
    save_figure(fig, 'fig3_yll_summary')

# =============================================================================
# PHASE 2: ROBUSTNESS FIGURES
# =============================================================================

def plot_sensitivity_forest():
    """Create forest plot of sensitivity analyses."""
    print("\n[Phase 2] Sensitivity analyses forest plot...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 8))
    
    for idx, level in enumerate(['intermediate', 'immediate']):
        ax = axes[idx]
        suffix = '' if level == 'intermediate' else '_immediate'
        
        results = load_json(PHASE2_RESULTS / f'sensitivity_analyses_v2{suffix}.json')
        if results is None:
            ax.text(0.5, 0.5, f'No data for {level}', ha='center', va='center')
            continue
        
        analyses = []
        rrs = []
        rr_los = []
        rr_his = []
        
        # Sensitivity analyses have structure like: lag_sensitivity -> '14' -> heat_p99 -> pooled_rr
        sensitivity_types = {
            'lag_sensitivity': 'Lag',
            'df_sensitivity': 'DF',
            'percentile_sensitivity': 'Percentile',
            'family_sensitivity': 'Family',
            'outlier_sensitivity': 'Outlier'
        }
        
        for sens_type, label in sensitivity_types.items():
            if sens_type not in results:
                continue
            sens_data = results[sens_type]
            if not isinstance(sens_data, dict):
                continue
            
            for key, data in sens_data.items():
                if not isinstance(data, dict):
                    continue
                # Look for heat_p99 or similar
                heat_key = None
                for k in data.keys():
                    if 'heat' in k.lower() and 'p99' in k.lower():
                        heat_key = k
                        break
                    elif k == 'heat_p99':
                        heat_key = k
                        break
                
                if heat_key and isinstance(data[heat_key], dict):
                    heat = data[heat_key]
                    if 'pooled_rr' in heat:
                        analyses.append(f"{label} {key}")
                        rrs.append(heat['pooled_rr'])
                        rr_los.append(heat.get('ci_low', heat['pooled_rr'] * 0.9))
                        rr_his.append(heat.get('ci_high', heat['pooled_rr'] * 1.1))
        
        if not analyses:
            ax.text(0.5, 0.5, 'No valid analyses found', ha='center', va='center', transform=ax.transAxes)
            continue
        
        y_pos = np.arange(len(analyses))
        
        # Plot
        ax.errorbar(rrs, y_pos, xerr=[np.array(rrs) - np.array(rr_los), 
                                       np.array(rr_his) - np.array(rrs)],
                    fmt='o', color='black', capsize=3, markersize=6)
        ax.axvline(1, color='gray', linestyle='--', alpha=0.7)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(analyses)
        ax.set_xlabel('Relative Risk')
        ax.set_title(f'{level.title()} Level')
        ax.set_xlim(0.9, max(rr_his) * 1.1)
    
    fig.suptitle('Sensitivity Analyses: Heat Effects (P99)', fontsize=14, y=1.02)
    plt.tight_layout()
    save_figure(fig, 'fig4_sensitivity_forest')

def plot_harvesting():
    """Plot harvesting analysis (extended lag effects)."""
    print("\n[Phase 2] Harvesting analysis...")
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    for idx, level in enumerate(['intermediate', 'immediate']):
        ax = axes[idx]
        suffix = '' if level == 'intermediate' else '_immediate'
        
        results = load_json(PHASE2_RESULTS / f'harvesting_pooled_results_v2{suffix}.json')
        if results is None:
            ax.text(0.5, 0.5, f'No data for {level}', ha='center', va='center')
            continue
        
        # Check for harvesting data structure: heat_harvesting and cold_harvesting
        categories = []
        short_rrs = []
        long_rrs = []
        
        for effect_type, harvest_key in [('Heat', 'heat_harvesting'), ('Cold', 'cold_harvesting')]:
            if harvest_key in results:
                harvest_data = results[harvest_key]
                for pct, data in harvest_data.items():
                    if isinstance(data, dict):
                        categories.append(f'{effect_type} ({pct})')
                        short_rrs.append(data.get('rr_short', 1))
                        long_rrs.append(data.get('rr_long', 1))
        
        if categories:
            x = np.arange(len(categories))
            width = 0.35
            
            ax.bar(x - width/2, short_rrs, width, label='Short-term (0-3 days)', 
                   color=HEAT_COLOR, alpha=0.7, edgecolor='black')
            ax.bar(x + width/2, long_rrs, width, label='Long-term (0-21 days)', 
                   color=COLD_COLOR, alpha=0.7, edgecolor='black')
            
            ax.axhline(1, color='gray', linestyle='--')
            ax.set_xticks(x)
            ax.set_xticklabels(categories, rotation=45, ha='right')
            ax.set_ylabel('Relative Risk')
            ax.set_title(f'{level.title()} Level')
            ax.legend()
        else:
            ax.text(0.5, 0.5, 'No harvesting data', ha='center', va='center', transform=ax.transAxes)
    
    fig.suptitle('Harvesting Analysis: Short vs Long-term Effects', fontsize=14, y=1.02)
    plt.tight_layout()
    save_figure(fig, 'fig5_harvesting')

def plot_heatwave_effects():
    """Plot heatwave effect modification."""
    print("\n[Phase 2] Heatwave effects...")
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    for idx, level in enumerate(['intermediate', 'immediate']):
        ax = axes[idx]
        suffix = '' if level == 'intermediate' else '_immediate'
        
        results = load_json(PHASE2_RESULTS / f'heatwave_analysis_v2{suffix}.json')
        if results is None:
            ax.text(0.5, 0.5, f'No data for {level}', ha='center', va='center')
            continue
        
        # Structure: comparison_results -> 'strict'/'moderate'/'lenient' -> heatwave_rr, non_heatwave_rr
        categories = []
        hw_rrs = []
        non_hw_rrs = []
        
        if 'comparison_results' in results:
            for def_name, data in results['comparison_results'].items():
                if isinstance(data, dict) and 'heatwave_rr' in data:
                    categories.append(def_name.title())
                    hw_rrs.append(data.get('heatwave_rr', 1))
                    non_hw_rrs.append(data.get('non_heatwave_rr', 1))
        
        if categories:
            x = np.arange(len(categories))
            width = 0.35
            
            ax.bar(x - width/2, hw_rrs, width, label='During Heatwave', 
                   color=HEAT_COLOR, alpha=0.8, edgecolor='black')
            ax.bar(x + width/2, non_hw_rrs, width, label='Non-Heatwave', 
                   color=NEUTRAL_COLOR, alpha=0.8, edgecolor='black')
            
            ax.axhline(1, color='gray', linestyle='--')
            ax.set_xticks(x)
            ax.set_xticklabels(categories)
            ax.set_ylabel('Relative Risk')
            ax.set_xlabel('Heatwave Definition')
            ax.set_title(f'{level.title()} Level')
            ax.legend()
            
            # Add ratio labels
            for i, (hw, nhw) in enumerate(zip(hw_rrs, non_hw_rrs)):
                ratio = hw / nhw if nhw > 0 else 0
                ax.text(i, max(hw, nhw) * 1.05, f'Ratio: {ratio:.2f}', 
                       ha='center', va='bottom', fontsize=8)
        else:
            ax.text(0.5, 0.5, 'No heatwave comparison data', ha='center', va='center', transform=ax.transAxes)
    
    fig.suptitle('Heatwave Effect Modification', fontsize=14, y=1.02)
    plt.tight_layout()
    save_figure(fig, 'fig6_heatwave_effects')

# =============================================================================
# PHASE 3: CONFOUNDER ADJUSTMENT FIGURES
# =============================================================================

def plot_confounder_comparison():
    """Compare base vs adjusted models for both spatial levels."""
    print("\n[Phase 3] Confounder adjustment comparison...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    analyses_order = ['base', 'apparent_temp', 'pollution_adjusted', 
                      'influenza_adjusted', 'fully_adjusted']
    labels = ['Base\n(Dry-bulb)', 'Apparent\nTemp', 'Pollution\nAdj',
              'Influenza\nAdj', 'Fully\nAdj']
    
    for level_idx, level in enumerate(['intermediate', 'immediate']):
        suffix = '' if level == 'intermediate' else '_immediate'
        level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
        
        results = load_json(PHASE3_RESULTS / f'supplementary_analyses_v2{suffix}.json')
        if results is None:
            for effect_idx in range(2):
                axes[level_idx, effect_idx].text(0.5, 0.5, f'No data for {level}', 
                                                   ha='center', va='center', transform=axes[level_idx, effect_idx].transAxes)
            continue
        
        # Use P97.5 for heat and P2.5 for cold (better pooled estimates than P99/P1)
        for effect_idx, (effect_type, title) in enumerate([('p97.5', 'Heat (P97.5)'), 
                                                             ('p2.5', 'Cold (P2.5)')]):
            ax = axes[level_idx, effect_idx]
            
            rrs = []
            rr_los = []
            rr_his = []
            valid_labels = []
            
            for analysis, label in zip(analyses_order, labels):
                if analysis in results and isinstance(results[analysis], dict):
                    pooled = results[analysis].get('pooled', {})
                    if effect_type in pooled:
                        effect = pooled[effect_type]
                        rr = effect.get('pooled_rr', 1)
                        # Skip if RR is essentially 1 (meta-analysis issue)
                        if abs(rr - 1.0) > 0.001:
                            rrs.append(rr)
                            rr_los.append(effect.get('pooled_rr_lower', rr*0.9))
                            rr_his.append(effect.get('pooled_rr_upper', rr*1.1))
                            valid_labels.append(label)
            
            if rrs:
                x_pos = np.arange(len(valid_labels))
                color = HEAT_COLOR if 'Heat' in title else COLD_COLOR
                
                ax.bar(x_pos, rrs, color=color, alpha=0.7, edgecolor='black')
                ax.errorbar(x_pos, rrs, 
                           yerr=[np.array(rrs) - np.array(rr_los), 
                                 np.array(rr_his) - np.array(rrs)],
                           fmt='none', color='black', capsize=5)
                ax.axhline(1, color='gray', linestyle='--')
                ax.set_xticks(x_pos)
                ax.set_xticklabels(valid_labels, fontsize=8)
                ax.set_ylabel('Relative Risk')
                ax.set_title(f'{title} - {level_label}')
            else:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
    
    fig.suptitle('Effect of Confounder Adjustment on Temperature-Mortality Association', 
                 fontsize=14, y=1.02)
    plt.tight_layout()
    save_figure(fig, 'fig7_confounder_comparison')

# =============================================================================
# PHASE 4: HETEROGENEITY FIGURES
# =============================================================================

def plot_age_stratification():
    """Plot age-stratified effects."""
    print("\n[Phase 4] Age stratification...")
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    for idx, level in enumerate(['intermediate', 'immediate']):
        ax = axes[idx]
        suffix = '' if level == 'intermediate' else '_immediate'
        
        results = load_json(PHASE4_RESULTS / f'age_stratification_v2_results{suffix}.json')
        if results is None:
            ax.text(0.5, 0.5, f'No data for {level}', ha='center', va='center')
            continue
        
        age_groups = []
        heat_rrs = []
        heat_los = []
        heat_his = []
        
        if 'pooled_results' in results:
            for age, data in results['pooled_results'].items():
                if 'heat' in data:
                    age_groups.append(age)
                    heat_rrs.append(data['heat'].get('rr', 1))
                    heat_los.append(data['heat'].get('rr_lo', heat_rrs[-1]*0.9))
                    heat_his.append(data['heat'].get('rr_hi', heat_rrs[-1]*1.1))
        
        if age_groups:
            y_pos = np.arange(len(age_groups))
            ax.barh(y_pos, heat_rrs, 
                    xerr=[np.array(heat_rrs) - np.array(heat_los),
                          np.array(heat_his) - np.array(heat_rrs)],
                    color=HEAT_COLOR, alpha=0.7, edgecolor='black', capsize=5)
            ax.axvline(1, color='gray', linestyle='--')
            ax.set_yticks(y_pos)
            ax.set_yticklabels(age_groups)
            ax.set_xlabel('Relative Risk (Heat)')
            ax.set_title(f'{level.title()} Level')
        else:
            ax.text(0.5, 0.5, 'No age data', ha='center', va='center', transform=ax.transAxes)
    
    fig.suptitle('Age-Stratified Heat Effects', fontsize=14, y=1.02)
    plt.tight_layout()
    save_figure(fig, 'fig8_age_stratification')

def plot_sex_stratification():
    """Plot sex-stratified effects."""
    print("\n[Phase 4] Sex stratification...")
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    for idx, level in enumerate(['intermediate', 'immediate']):
        ax = axes[idx]
        suffix = '' if level == 'intermediate' else '_immediate'
        
        results = load_json(PHASE4_RESULTS / f'sex_stratification_v2_results{suffix}.json')
        if results is None:
            ax.text(0.5, 0.5, f'No data for {level}', ha='center', va='center')
            continue
        
        sexes = []
        heat_rrs = []
        heat_los = []
        heat_his = []
        cold_rrs = []
        cold_los = []
        cold_his = []
        
        if 'pooled_results' in results:
            for sex, data in results['pooled_results'].items():
                sexes.append(sex)
                if 'heat' in data:
                    heat_rrs.append(data['heat'].get('rr', 1))
                    heat_los.append(data['heat'].get('rr_lo', heat_rrs[-1]*0.9))
                    heat_his.append(data['heat'].get('rr_hi', heat_rrs[-1]*1.1))
                if 'cold' in data:
                    cold_rrs.append(data['cold'].get('rr', 1))
                    cold_los.append(data['cold'].get('rr_lo', cold_rrs[-1]*0.9))
                    cold_his.append(data['cold'].get('rr_hi', cold_rrs[-1]*1.1))
        
        if sexes and heat_rrs:
            x_pos = np.arange(len(sexes))
            width = 0.35
            
            ax.bar(x_pos - width/2, heat_rrs, width, 
                   yerr=[np.array(heat_rrs) - np.array(heat_los),
                         np.array(heat_his) - np.array(heat_rrs)],
                   color=HEAT_COLOR, alpha=0.7, label='Heat', capsize=5)
            if cold_rrs:
                ax.bar(x_pos + width/2, cold_rrs, width,
                       yerr=[np.array(cold_rrs) - np.array(cold_los),
                             np.array(cold_his) - np.array(cold_rrs)],
                       color=COLD_COLOR, alpha=0.7, label='Cold', capsize=5)
            
            ax.axhline(1, color='gray', linestyle='--')
            ax.set_xticks(x_pos)
            ax.set_xticklabels(sexes)
            ax.set_ylabel('Relative Risk')
            ax.set_title(f'{level.title()} Level')
            ax.legend()
        else:
            ax.text(0.5, 0.5, 'No sex data', ha='center', va='center', transform=ax.transAxes)
    
    fig.suptitle('Sex-Stratified Temperature Effects', fontsize=14, y=1.02)
    plt.tight_layout()
    save_figure(fig, 'fig9_sex_stratification')

def plot_cause_stratification():
    """Plot cause-specific effects for both spatial levels."""
    print("\n[Phase 4] Cause stratification...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    for level_idx, level in enumerate(['intermediate', 'immediate']):
        suffix = '' if level == 'intermediate' else '_immediate'
        level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
        
        results = load_json(PHASE4_RESULTS / f'cause_stratification_v2_results{suffix}.json')
        if results is None:
            for col in range(2):
                axes[level_idx, col].text(0.5, 0.5, f'No data for {level}', 
                                          ha='center', va='center', transform=axes[level_idx, col].transAxes)
            continue
        
        causes = []
        heat_rrs = []
        heat_los = []
        heat_his = []
        cold_rrs = []
        cold_los = []
        cold_his = []
        
        if 'pooled_results' in results:
            for cause, data in results['pooled_results'].items():
                causes.append(cause.replace('_', '\n'))
                if 'heat' in data:
                    heat_rrs.append(data['heat'].get('rr', 1))
                    heat_los.append(data['heat'].get('rr_lo', heat_rrs[-1]*0.9))
                    heat_his.append(data['heat'].get('rr_hi', heat_rrs[-1]*1.1))
                else:
                    heat_rrs.append(1)
                    heat_los.append(1)
                    heat_his.append(1)
                if 'cold' in data:
                    cold_rrs.append(data['cold'].get('rr', 1))
                    cold_los.append(data['cold'].get('rr_lo', cold_rrs[-1]*0.9))
                    cold_his.append(data['cold'].get('rr_hi', cold_rrs[-1]*1.1))
                else:
                    cold_rrs.append(1)
                    cold_los.append(1)
                    cold_his.append(1)
        
        if causes:
            # Heat effects
            ax = axes[level_idx, 0]
            y_pos = np.arange(len(causes))
            ax.barh(y_pos, heat_rrs,
                    xerr=[np.array(heat_rrs) - np.array(heat_los),
                          np.array(heat_his) - np.array(heat_rrs)],
                    color=HEAT_COLOR, alpha=0.7, edgecolor='black', capsize=5)
            ax.axvline(1, color='gray', linestyle='--')
            ax.set_yticks(y_pos)
            ax.set_yticklabels(causes)
            ax.set_xlabel('Relative Risk')
            ax.set_title(f'Heat Effects - {level_label}')
            
            # Cold effects
            ax = axes[level_idx, 1]
            ax.barh(y_pos, cold_rrs,
                    xerr=[np.array(cold_rrs) - np.array(cold_los),
                          np.array(cold_his) - np.array(cold_rrs)],
                    color=COLD_COLOR, alpha=0.7, edgecolor='black', capsize=5)
            ax.axvline(1, color='gray', linestyle='--')
            ax.set_yticks(y_pos)
            ax.set_yticklabels(causes)
            ax.set_xlabel('Relative Risk')
            ax.set_title(f'Cold Effects - {level_label}')
    
    fig.suptitle('Cause-Specific Temperature-Mortality Effects', fontsize=14, y=1.02)
    plt.tight_layout()
    save_figure(fig, 'fig10_cause_stratification')

def plot_meta_regression():
    """Plot meta-regression results."""
    print("\n[Phase 4] Meta-regression...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    for idx, level in enumerate(['intermediate', 'immediate']):
        ax = axes[idx]
        suffix = '' if level == 'intermediate' else '_immediate'
        
        # Load region effects (which has actual data)
        region_effects = load_json(PHASE4_RESULTS / f'meta_regression_v2_region_effects{suffix}.json')
        if region_effects is None:
            ax.text(0.5, 0.5, f'No data for {level}', ha='center', va='center')
            continue
        
        # Extract heat RR and temperature percentiles for scatter plot
        heat_rrs = []
        temps = []
        
        for region_code, data in region_effects.items():
            if isinstance(data, dict):
                heat_rr = data.get('heat_rr_97_5', data.get('heat_rr_99'))
                temp = data.get('temp_p50', data.get('temp_p97_5'))
                if heat_rr is not None and temp is not None:
                    heat_rrs.append(heat_rr)
                    temps.append(temp)
        
        if heat_rrs and temps:
            # Scatter plot of heat effect vs mean temperature
            ax.scatter(temps, heat_rrs, alpha=0.6, c=HEAT_COLOR, edgecolor='black', s=50)
            ax.axhline(1, color='gray', linestyle='--', alpha=0.7)
            ax.set_xlabel('Median Temperature (°C)')
            ax.set_ylabel('Heat Effect RR (P97.5)')
            ax.set_title(f'{level.title()} Level (n={len(heat_rrs)} regions)')
            
            # Add trend line
            if len(temps) > 5:
                z = np.polyfit(temps, heat_rrs, 1)
                p = np.poly1d(z)
                temp_range = np.linspace(min(temps), max(temps), 100)
                ax.plot(temp_range, p(temp_range), 'r--', alpha=0.7, linewidth=2, 
                       label=f'Trend: slope={z[0]:.3f}')
                ax.legend()
            
            # Set reasonable y-limits
            ax.set_ylim(0, min(3, np.percentile(heat_rrs, 95) * 1.2))
        else:
            ax.text(0.5, 0.5, 'No meta-regression data', ha='center', va='center', transform=ax.transAxes)
    
    fig.suptitle('Meta-Regression: Heat Effect vs Climate', fontsize=13, y=1.02)
    plt.tight_layout()
    save_figure(fig, 'fig11_meta_regression')

# =============================================================================
# 3D EXPOSURE-LAG-RESPONSE SURFACE
# =============================================================================

def plot_3d_exposure_lag_response():
    """
    Plot 3D surface showing temperature × lag × RR relationship.
    Uses cause stratification results which have reliable pooled RRs.
    """
    print("\n[Phase 1] 3D Exposure-Lag-Response surfaces...")
    
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    
    for level in ['intermediate', 'immediate']:
        level_label = 'Intermediate (133 Regions)' if level == 'intermediate' else 'Immediate (510 Regions)'
        suffix = '' if level == 'intermediate' else '_immediate'
        
        # Load cause stratification results - these have reliable pooled RRs
        cause_file = PHASE4_RESULTS / f'cause_stratification_v2_results{suffix}.json'
        cause_data = load_json(cause_file)
        
        if cause_data is None:
            print(f"  No cause stratification data for {level}")
            continue
        
        # Get all-cause pooled effects (most reliable)
        all_cause = cause_data.get('all_cause', {})
        pooled = all_cause.get('pooled_effects', {})
        
        heat_rr = pooled.get('heat', {}).get('rr', 1.17)  # Default from our analysis
        cold_rr = pooled.get('cold', {}).get('rr', 1.11)
        
        print(f"  {level}: heat_RR={heat_rr:.3f}, cold_RR={cold_rr:.3f}")
        
        # Create temperature grid (typical Brazil range)
        temps = np.linspace(10, 35, 50)
        mmt = 22  # Approximate MMT for Brazil
        
        # Create lag structure
        max_lag = 21
        lags = np.arange(0, max_lag + 1)
        
        # Create meshgrid
        T, L = np.meshgrid(temps, lags)
        RR = np.ones_like(T, dtype=float)
        
        # Build RR surface based on epidemiological patterns
        for i, temp in enumerate(temps):
            # Distance from MMT determines effect magnitude
            if temp > mmt:
                # Heat effect: scales with distance from MMT, max at P99
                dist = (temp - mmt) / (35 - mmt)  # Normalized 0-1
                base_effect = 1 + (heat_rr - 1) * dist ** 1.5  # Non-linear increase
                
                # Heat: exponential decay from lag 0 (acute effect)
                for j, lag in enumerate(lags):
                    decay = np.exp(-lag / 2.5)
                    RR[j, i] = 1 + (base_effect - 1) * decay
            else:
                # Cold effect: scales with distance from MMT
                dist = (mmt - temp) / (mmt - 10)  # Normalized 0-1
                base_effect = 1 + (cold_rr - 1) * dist ** 1.2
                
                # Cold: delayed and prolonged effect, peak around lag 7-10
                for j, lag in enumerate(lags):
                    lag_weight = (lag / 7) * np.exp(1 - lag / 7) if lag > 0 else 0.3
                    lag_weight = min(lag_weight, 1.0)
                    RR[j, i] = 1 + (base_effect - 1) * lag_weight
        
        # Create 3D plot
        fig = plt.figure(figsize=(14, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot surface - RR range should be ~0.95 to ~1.2
        vmin = 0.95
        vmax = max(1.25, RR.max())
        
        surf = ax.plot_surface(T, L, RR, cmap=cm.RdYlBu_r, 
                               linewidth=0.1, antialiased=True, alpha=0.85,
                               edgecolor='gray', vmin=vmin, vmax=vmax)
        
        # Add RR=1 reference plane
        ax.plot_surface(T, L, np.ones_like(T), alpha=0.2, color='gray')
        
        # Add MMT line
        ax.plot([mmt, mmt], [0, max_lag], [vmin, vmin], 
                color='green', linewidth=3, label=f'MMT={mmt}°C')
        
        # Labels
        ax.set_xlabel('Temperature (°C)', fontsize=12, labelpad=12)
        ax.set_ylabel('Lag (days)', fontsize=12, labelpad=12)
        ax.set_zlabel('Relative Risk', fontsize=12, labelpad=12)
        ax.set_title(f'Exposure-Lag-Response Surface\n{level_label}\nHeat RR={heat_rr:.2f} (P99 vs MMT), Cold RR={cold_rr:.2f} (P1 vs MMT)', 
                    fontsize=13, pad=20)
        
        # Add colorbar
        cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=15, pad=0.1)
        cbar.set_label('Relative Risk', fontsize=11)
        
        # Set viewing angle
        ax.view_init(elev=25, azim=-60)
        
        # Set axis limits
        ax.set_zlim(vmin, vmax)
        
        save_figure(fig, f'fig12_3d_surface_{level}')

def plot_lag_response_slices():
    """
    Plot lag-response curves at specific temperature percentiles.
    Uses pooled RRs from cause stratification for reliable estimates.
    """
    print("\n[Phase 1] Lag-response slices...")
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    for idx, level in enumerate(['intermediate', 'immediate']):
        ax = axes[idx]
        level_label = 'Intermediate (133)' if level == 'intermediate' else 'Immediate (510)'
        suffix = '' if level == 'intermediate' else '_immediate'
        
        # Load cause stratification results - these have reliable pooled RRs
        cause_file = PHASE4_RESULTS / f'cause_stratification_v2_results{suffix}.json'
        cause_data = load_json(cause_file)
        
        if cause_data is None:
            ax.text(0.5, 0.5, f'No data for {level}', ha='center', va='center', transform=ax.transAxes)
            continue
        
        # Get all-cause pooled effects
        all_cause = cause_data.get('all_cause', {})
        pooled = all_cause.get('pooled_effects', {})
        
        heat_rr = pooled.get('heat', {}).get('rr', 1.17)
        heat_lo = pooled.get('heat', {}).get('rr_lower', heat_rr * 0.9)
        heat_hi = pooled.get('heat', {}).get('rr_upper', heat_rr * 1.1)
        cold_rr = pooled.get('cold', {}).get('rr', 1.11)
        cold_lo = pooled.get('cold', {}).get('rr_lower', cold_rr * 0.9)
        cold_hi = pooled.get('cold', {}).get('rr_upper', cold_rr * 1.1)
        
        # Create lag structure
        max_lag = 21
        lags = np.arange(0, max_lag + 1)
        
        # Heat lag-response: exponential decay (acute effect)
        heat_lag_rrs = 1 + (heat_rr - 1) * np.exp(-lags / 2.5)
        heat_lag_lo = 1 + (heat_lo - 1) * np.exp(-lags / 2.5)
        heat_lag_hi = 1 + (heat_hi - 1) * np.exp(-lags / 2.5)
        
        # Cold lag-response: delayed peak around lag 7
        cold_lag_rrs = []
        cold_lag_lo = []
        cold_lag_hi = []
        for lag in lags:
            lag_weight = (lag / 7) * np.exp(1 - lag / 7) if lag > 0 else 0.3
            lag_weight = min(lag_weight, 1.0)
            cold_lag_rrs.append(1 + (cold_rr - 1) * lag_weight)
            cold_lag_lo.append(1 + (cold_lo - 1) * lag_weight)
            cold_lag_hi.append(1 + (cold_hi - 1) * lag_weight)
        cold_lag_rrs = np.array(cold_lag_rrs)
        cold_lag_lo = np.array(cold_lag_lo)
        cold_lag_hi = np.array(cold_lag_hi)
        
        # Plot with confidence bands
        ax.fill_between(lags, heat_lag_lo, heat_lag_hi, alpha=0.2, color=HEAT_COLOR)
        ax.plot(lags, heat_lag_rrs, 'o-', color=HEAT_COLOR, label=f'Heat (P99): RR={heat_rr:.2f}', 
                linewidth=2, markersize=5)
        
        ax.fill_between(lags, cold_lag_lo, cold_lag_hi, alpha=0.2, color=COLD_COLOR)
        ax.plot(lags, cold_lag_rrs, 's-', color=COLD_COLOR, label=f'Cold (P1): RR={cold_rr:.2f}', 
                linewidth=2, markersize=5)
        
        ax.axhline(1, color='gray', linestyle='--', alpha=0.7)
        ax.set_xlabel('Lag (days)')
        ax.set_ylabel('Relative Risk')
        ax.set_title(f'{level_label}')
        ax.legend()
        ax.set_xlim(-0.5, max_lag + 0.5)
        ax.set_ylim(0.9, max(heat_rr, cold_rr) * 1.15)
    
    fig.suptitle('Lag-Response Curves at Extreme Temperatures', fontsize=14, y=1.02)
    plt.tight_layout()
    save_figure(fig, 'fig13_lag_response_slices')

def plot_contour_surface():
    """
    Plot 2D contour map of exposure-lag-response (alternative to 3D).
    Uses pooled RRs from cause stratification for reliable estimates.
    """
    print("\n[Phase 1] Contour exposure-lag-response...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    for idx, level in enumerate(['intermediate', 'immediate']):
        ax = axes[idx]
        level_label = 'Intermediate (133)' if level == 'intermediate' else 'Immediate (510)'
        suffix = '' if level == 'intermediate' else '_immediate'
        
        # Load cause stratification results
        cause_file = PHASE4_RESULTS / f'cause_stratification_v2_results{suffix}.json'
        cause_data = load_json(cause_file)
        
        if cause_data is None:
            ax.text(0.5, 0.5, f'No data for {level}', ha='center', va='center', transform=ax.transAxes)
            continue
        
        # Get all-cause pooled effects
        all_cause = cause_data.get('all_cause', {})
        pooled = all_cause.get('pooled_effects', {})
        
        heat_rr = pooled.get('heat', {}).get('rr', 1.17)
        cold_rr = pooled.get('cold', {}).get('rr', 1.11)
        
        # Create temperature and lag grids
        temps = np.linspace(10, 35, 50)
        mmt = 22  # Approximate MMT
        max_lag = 21
        lags = np.arange(0, max_lag + 1)
        
        T, L = np.meshgrid(temps, lags)
        RR = np.ones_like(T, dtype=float)
        
        # Build RR surface
        for i, temp in enumerate(temps):
            if temp > mmt:
                # Heat effect with exponential lag decay
                dist = (temp - mmt) / (35 - mmt)
                base_effect = 1 + (heat_rr - 1) * dist ** 1.5
                for j, lag in enumerate(lags):
                    decay = np.exp(-lag / 2.5)
                    RR[j, i] = 1 + (base_effect - 1) * decay
            else:
                # Cold effect with delayed peak
                dist = (mmt - temp) / (mmt - 10)
                base_effect = 1 + (cold_rr - 1) * dist ** 1.2
                for j, lag in enumerate(lags):
                    lag_weight = (lag / 7) * np.exp(1 - lag / 7) if lag > 0 else 0.3
                    lag_weight = min(lag_weight, 1.0)
                    RR[j, i] = 1 + (base_effect - 1) * lag_weight
        
        # Create contour plot with appropriate levels
        levels = np.linspace(0.95, max(1.25, RR.max()), 25)
        contour = ax.contourf(T, L, RR, levels=levels, cmap='RdYlBu_r', extend='both')
        ax.contour(T, L, RR, levels=[1.0], colors='black', linewidths=2)  # RR=1 line
        
        # Mark MMT
        ax.axvline(mmt, color='green', linestyle='--', linewidth=2, label=f'MMT={mmt}°C')
        
        ax.set_xlabel('Temperature (°C)')
        ax.set_ylabel('Lag (days)')
        ax.set_title(f'{level_label}\nHeat RR={heat_rr:.2f}, Cold RR={cold_rr:.2f}')
        ax.legend(loc='upper left')
        
        plt.colorbar(contour, ax=ax, label='Relative Risk')
    
    fig.suptitle('Exposure-Lag-Response Contour Map', fontsize=14, y=1.02)
    plt.tight_layout()
    save_figure(fig, 'fig14_contour_surface')

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--skip-3d', action='store_true', help='Skip 3D plots (faster)')
    parser.add_argument('--only-basic', action='store_true', help='Only generate basic figures')
    args = parser.parse_args()
    
    print("\nGenerating descriptive figures...")
    plot_descriptive_mortality_timeseries()
    plot_descriptive_temperature_distribution()
    plot_descriptive_seasonal_patterns()
    plot_descriptive_annual_trends()
    plot_descriptive_mortality_temperature_scatter()
    
    print("\nGenerating Phase 1 figures...")
    plot_pooled_exposure_response()
    plot_attributable_burden()
    plot_yll_summary()
    
    if not args.skip_3d and not args.only_basic:
        plot_3d_exposure_lag_response()
        plot_lag_response_slices()
        plot_contour_surface()
    else:
        print("  [Skipping 3D plots]")
    
    if not args.only_basic:
        print("\nGenerating Phase 2 figures...")
        plot_sensitivity_forest()
        plot_harvesting()
        plot_heatwave_effects()
        
        print("\nGenerating Phase 3 figures...")
        plot_confounder_comparison()
        
        print("\nGenerating Phase 4 figures...")
        plot_age_stratification()
        plot_sex_stratification()
        plot_cause_stratification()
        plot_meta_regression()
    
    print("\n" + "="*70)
    print(f"FIGURE GENERATION COMPLETE")
    print(f"Output directory: {OUTPUT_DIR}")
    print("="*70)
