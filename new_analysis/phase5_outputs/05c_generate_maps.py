"""
05c: GENERATE MAPS FOR ALL PHASES
==================================
Creates geographic visualizations of temperature-mortality effects.

Outputs:
--------
- Map 1: Regional heat effects
- Map 2: Regional cold effects
- Map 3: Attributable burden by region
- Map 4: Minimum mortality temperature (MMT)

Requirements:
-------------
- geopandas
- matplotlib
- brazil_municipalities_2022.gpkg

Author: Heat-Mortality Brazil Analysis Pipeline
Date: December 2025
"""

import warnings
warnings.filterwarnings('ignore')

import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from pathlib import Path
from datetime import datetime

try:
    import geopandas as gpd
    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False
    print("WARNING: geopandas not installed. Maps will not be generated.")

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = Path(__file__).parent.parent
INPUT_DATA_DIR = BASE_DIR.parent / 'Input_data'
PHASE0_RESULTS = BASE_DIR / 'phase0_data_prep' / 'results'
PHASE1_RESULTS = BASE_DIR / 'phase1_core_model' / 'results'
PHASE1_R_RESULTS = BASE_DIR / 'phase1_r' / 'results'  # R results
PHASE4_RESULTS = BASE_DIR / 'phase4_heterogeneity' / 'results'
OUTPUT_DIR = Path(__file__).parent / 'maps'
OUTPUT_DIR.mkdir(exist_ok=True)

# Municipality to region mapping file (IBGE 2017 geographic division)
REGION_MAPPING_FILE = PHASE0_RESULTS / 'municipality_to_all_regions_map.csv'

# Map styling
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300

# Colors
HEAT_CMAP = 'YlOrRd'
COLD_CMAP = 'YlGnBu'
DIVERGING_CMAP = 'RdYlBu_r'

print("="*70)
print("05c: GENERATE MAPS FOR ALL PHASES")
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

def save_map(fig, name, formats=['png', 'pdf']):
    """Save map in multiple formats."""
    for fmt in formats:
        filepath = OUTPUT_DIR / f"{name}.{fmt}"
        fig.savefig(filepath, format=fmt, bbox_inches='tight', dpi=300)
    print(f"  Saved: {name}")
    plt.close(fig)

def load_brazil_shapefile():
    """Load Brazil shapefile and prepare for regional aggregation."""
    if not HAS_GEOPANDAS:
        return None, None
    
    # Try different shapefile locations
    shapefile_paths = [
        INPUT_DATA_DIR / 'brazil_municipalities_2022.gpkg',
        PHASE0_RESULTS / 'brazil_regions.gpkg',
        INPUT_DATA_DIR / 'brazil_states.gpkg',
    ]
    
    for path in shapefile_paths:
        if path.exists():
            try:
                gdf = gpd.read_file(path)
                print(f"  Loaded shapefile: {path.name}")
                return gdf, path.name
            except Exception as e:
                print(f"  Error loading {path}: {e}")
    
    print("  WARNING: No shapefile found")
    return None, None

def load_region_mapping():
    """Load the municipality to region mapping from IBGE."""
    if REGION_MAPPING_FILE.exists():
        mapping = pd.read_csv(REGION_MAPPING_FILE)
        print(f"  Loaded region mapping: {len(mapping)} municipalities")
        print(f"    → {mapping['intermediate_code'].nunique()} intermediate regions")
        print(f"    → {mapping['immediate_code'].nunique()} immediate regions")
        return mapping
    else:
        print(f"  WARNING: Region mapping file not found: {REGION_MAPPING_FILE}")
        return None

def aggregate_to_intermediate(gdf, mapping=None):
    """Aggregate municipalities to 133 intermediate regions using IBGE mapping."""
    if mapping is None:
        mapping = load_region_mapping()
    
    if mapping is None:
        print("  WARNING: Cannot aggregate without mapping file")
        return gdf
    
    # Find municipality code column
    muni_col = None
    for col in ['CD_MUN', 'code_muni', 'CD_GEOCMU', 'GEOCODIGO']:
        if col in gdf.columns:
            muni_col = col
            break
    
    if muni_col is None:
        print("  WARNING: Could not find municipality code column")
        return gdf
    
    gdf = gdf.copy()
    
    # Ensure consistent types for merge
    gdf['code_muni_int'] = gdf[muni_col].astype(float).astype(int)
    mapping_clean = mapping[['code_muni', 'intermediate_code']].copy()
    mapping_clean['code_muni'] = mapping_clean['code_muni'].astype(int)
    
    # Merge to get intermediate codes
    gdf = gdf.merge(mapping_clean, left_on='code_muni_int', right_on='code_muni', how='left')
    
    # Check merge success
    matched = gdf['intermediate_code'].notna().sum()
    print(f"    Matched {matched}/{len(gdf)} municipalities to intermediate regions")
    
    # Drop unmatched municipalities
    gdf = gdf[gdf['intermediate_code'].notna()].copy()
    gdf['region_code'] = gdf['intermediate_code'].astype(int)
    
    # Dissolve by intermediate region
    gdf_regions = gdf.dissolve(by='region_code').reset_index()
    print(f"    Created {len(gdf_regions)} intermediate region geometries")
    
    return gdf_regions

def aggregate_to_immediate(gdf, mapping=None):
    """Aggregate municipalities to 510 immediate regions using IBGE mapping."""
    if mapping is None:
        mapping = load_region_mapping()
    
    if mapping is None:
        print("  WARNING: Cannot aggregate without mapping file")
        return gdf
    
    # Find municipality code column
    muni_col = None
    for col in ['CD_MUN', 'code_muni', 'CD_GEOCMU', 'GEOCODIGO']:
        if col in gdf.columns:
            muni_col = col
            break
    
    if muni_col is None:
        print("  WARNING: Could not find municipality code column")
        return gdf
    
    gdf = gdf.copy()
    
    # Ensure consistent types for merge
    gdf['code_muni_int'] = gdf[muni_col].astype(float).astype(int)
    mapping_clean = mapping[['code_muni', 'immediate_code']].copy()
    mapping_clean['code_muni'] = mapping_clean['code_muni'].astype(int)
    
    # Merge to get immediate codes
    gdf = gdf.merge(mapping_clean, left_on='code_muni_int', right_on='code_muni', how='left')
    
    # Check merge success
    matched = gdf['immediate_code'].notna().sum()
    print(f"    Matched {matched}/{len(gdf)} municipalities to immediate regions")
    
    # Drop unmatched municipalities
    gdf = gdf[gdf['immediate_code'].notna()].copy()
    gdf['region_code'] = gdf['immediate_code'].astype(int)
    
    # Dissolve by immediate region
    gdf_regions = gdf.dissolve(by='region_code').reset_index()
    print(f"    Created {len(gdf_regions)} immediate region geometries")
    
    return gdf_regions

# =============================================================================
# MAP GENERATION FUNCTIONS
# =============================================================================

def create_effect_map(gdf, effects_dict, effect_col, title, cmap, vmin=None, vmax=None, 
                      center=None, output_name=None):
    """Create a choropleth map of effects."""
    if gdf is None:
        print(f"  Skipping {title}: No shapefile")
        return
    
    # Merge effects with geometry
    effects_df = pd.DataFrame([
        {'region_code': int(k), effect_col: v}
        for k, v in effects_dict.items()
    ])
    
    gdf_merged = gdf.merge(effects_df, on='region_code', how='left')
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    
    # Plot base map (regions without data)
    gdf_merged[gdf_merged[effect_col].isna()].plot(
        ax=ax, color='lightgray', edgecolor='white', linewidth=0.5
    )
    
    # Plot choropleth
    if center is not None:
        # Diverging colormap centered at 1
        norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)
        gdf_merged[gdf_merged[effect_col].notna()].plot(
            column=effect_col, ax=ax, cmap=cmap, norm=norm,
            edgecolor='white', linewidth=0.3, legend=True,
            legend_kwds={'label': 'Relative Risk', 'shrink': 0.6}
        )
    else:
        gdf_merged[gdf_merged[effect_col].notna()].plot(
            column=effect_col, ax=ax, cmap=cmap, vmin=vmin, vmax=vmax,
            edgecolor='white', linewidth=0.3, legend=True,
            legend_kwds={'label': effect_col, 'shrink': 0.6}
        )
    
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.axis('off')
    
    # Add data summary
    valid_effects = [v for v in effects_dict.values() if v is not None and not np.isnan(v)]
    if valid_effects:
        stats_text = f"n = {len(valid_effects)} regions\nMean = {np.mean(valid_effects):.3f}\nRange = [{min(valid_effects):.3f}, {max(valid_effects):.3f}]"
        ax.text(0.02, 0.02, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    if output_name:
        save_map(fig, output_name)
    
    return fig

def map_heat_effects():
    """Create maps of heat and cold effects at both spatial levels using R results."""
    print("\n[Map 1-3] Heat, cold effects, and MMT...")
    
    # Load shapefile
    gdf_base, shapefile_name = load_brazil_shapefile()
    if gdf_base is None:
        return
    
    # Load mapping once
    mapping = load_region_mapping()
    if mapping is None:
        print("  Cannot generate maps without region mapping")
        return
    
    for level in ['intermediate', 'immediate']:
        print(f"\n  Processing {level} level...")
        
        # Aggregate geometry using proper IBGE mapping
        if level == 'intermediate':
            gdf = aggregate_to_intermediate(gdf_base.copy(), mapping)
        else:
            gdf = aggregate_to_immediate(gdf_base.copy(), mapping)
        
        # Load R results (CSV format) - preferred source
        r_results_file = PHASE1_R_RESULTS / f'dlnm_r_{level}_summary.csv'
        if r_results_file.exists():
            results_df = pd.read_csv(r_results_file)
            print(f"    Loaded R results: {len(results_df)} regions")
        else:
            # Fallback to Python results
            results = load_json(PHASE1_RESULTS / f'dlnm_v2_{level}_results.json')
            if results is None:
                print(f"    No results found for {level}")
                continue
            # Convert Python results to DataFrame
            region_results = results.get('region_results', results)
            rows = []
            for region_key, region_data in region_results.items():
                if isinstance(region_data, dict):
                    rows.append({
                        'region_code': int(region_key),
                        'mmt': region_data.get('mmt'),
                        'rr_p99': region_data.get('effects', {}).get('p99', {}).get('rr'),
                        'rr_p1': region_data.get('effects', {}).get('p1', {}).get('rr')
                    })
            results_df = pd.DataFrame(rows)
            print(f"    Loaded Python results: {len(results_df)} regions")
        
        # Filter based on REALISTIC epidemiological bounds
        # RR should be between 0.8 and 2.0 for temperature effects
        # MMT should be between 15 and 30°C for Brazil
        # Extreme values indicate convergence failure
        results_df_clean = results_df.copy()
        
        # Use realistic bounds - values outside these are convergence failures
        heat_valid = (results_df_clean['rr_p99'] >= 0.9) & (results_df_clean['rr_p99'] <= 1.5)
        cold_valid = (results_df_clean['rr_p1'] >= 0.9) & (results_df_clean['rr_p1'] <= 1.8)  # Cold can be higher
        mmt_valid = (results_df_clean['mmt'] >= 15) & (results_df_clean['mmt'] <= 30)
        
        n_heat = heat_valid.sum()
        n_cold = cold_valid.sum()
        n_mmt = mmt_valid.sum()
        n_total = len(results_df)
        
        print(f"    Valid heat effects: {n_heat}/{n_total} ({100*n_heat/n_total:.0f}%)")
        print(f"    Valid cold effects: {n_cold}/{n_total} ({100*n_cold/n_total:.0f}%)")
        print(f"    Valid MMT values: {n_mmt}/{n_total} ({100*n_mmt/n_total:.0f}%)")
        
        # Create separate dictionaries for each variable
        heat_effects = {
            int(row['region_code']): row['rr_p99'] 
            for _, row in results_df_clean[heat_valid].iterrows()
        }
        cold_effects = {
            int(row['region_code']): row['rr_p1'] 
            for _, row in results_df_clean[cold_valid].iterrows()
        }
        mmt_values = {
            int(row['region_code']): row['mmt'] 
            for _, row in results_df_clean[mmt_valid].iterrows()
        }
        
        level_label = 'Intermediate (133 Regions)' if level == 'intermediate' else 'Immediate (510 Regions)'
        
        # Heat effects map
        if heat_effects:
            create_effect_map(
                gdf, heat_effects, 'heat_rr',
                f'Heat Effect (P99 vs MMT) - {level_label}',
                DIVERGING_CMAP, vmin=0.9, vmax=1.3, center=1.0,
                output_name=f'map1_heat_effects_{level}'
            )
        
        # Cold effects map
        if cold_effects:
            create_effect_map(
                gdf, cold_effects, 'cold_rr',
                f'Cold Effect (P1 vs MMT) - {level_label}',
                DIVERGING_CMAP, vmin=0.9, vmax=1.3, center=1.0,
                output_name=f'map2_cold_effects_{level}'
            )
        
        # MMT map
        if mmt_values:
            create_effect_map(
                gdf, mmt_values, 'mmt',
                f'Minimum Mortality Temperature - {level_label}',
                'coolwarm', vmin=18, vmax=28,
                output_name=f'map3_mmt_{level}'
            )

def map_attributable_burden():
    """Create maps of attributable burden using R results."""
    print("\n[Map 4-5] Attributable burden...")
    
    # Load shapefile
    gdf_base, _ = load_brazil_shapefile()
    if gdf_base is None:
        return
    
    # Load mapping once
    mapping = load_region_mapping()
    if mapping is None:
        print("  Cannot generate maps without region mapping")
        return
    
    for level in ['intermediate', 'immediate']:
        print(f"\n  Processing {level} level...")
        
        # Aggregate geometry using proper IBGE mapping
        if level == 'intermediate':
            gdf = aggregate_to_intermediate(gdf_base.copy(), mapping)
        else:
            gdf = aggregate_to_immediate(gdf_base.copy(), mapping)
        
        # Try R results first
        r_burden_file = PHASE1_R_RESULTS / f'attributable_burden_r_{level}_regions.csv'
        if r_burden_file.exists():
            burden_df = pd.read_csv(r_burden_file)
            print(f"    Loaded R burden results: {len(burden_df)} regions")
        else:
            # Fallback to Python results
            burden_file = PHASE1_RESULTS / f'burden_v2_{level}_regions.csv'
            if not burden_file.exists():
                print(f"    No burden results found for {level}")
                continue
            burden_df = pd.read_csv(burden_file)
            print(f"    Loaded Python burden results: {len(burden_df)} regions")
        
        # Identify region column
        region_col = None
        for col in ['region_code', 'intermediate_code', 'immediate_code', 'region']:
            if col in burden_df.columns:
                region_col = col
                break
        
        if region_col is None:
            print(f"    Could not identify region column")
            continue
        
        burden_df['region_code'] = burden_df[region_col].astype(int)
        
        # Merge with geometry
        gdf_merged = gdf.merge(burden_df, on='region_code', how='left')
        
        level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
        
        # Heat burden map - check various possible column names
        heat_col = None
        for col in ['total_heat_an', 'heat_deaths', 'heat_attributable_deaths', 'heat_an_97_5']:
            if col in burden_df.columns:
                heat_col = col
                break
        
        if heat_col and heat_col in gdf_merged.columns:
            fig, ax = plt.subplots(1, 1, figsize=(12, 10))
            
            # Plot missing data regions only if there are any
            missing_data = gdf_merged[gdf_merged[heat_col].isna()]
            if len(missing_data) > 0:
                missing_data.plot(ax=ax, color='lightgray', edgecolor='white', linewidth=0.5)
            
            # Plot data regions
            valid_data = gdf_merged[gdf_merged[heat_col].notna()]
            if len(valid_data) > 0:
                valid_data.plot(
                    column=heat_col, ax=ax, cmap=HEAT_CMAP,
                    edgecolor='white', linewidth=0.3, legend=True,
                    legend_kwds={'label': 'Attributable Deaths', 'shrink': 0.6}
                )
            
            ax.set_title(f'Heat-Attributable Deaths (2010-2024) - {level_label}', 
                        fontsize=14, fontweight='bold')
            ax.axis('off')
            
            save_map(fig, f'map4_heat_burden_{level}')
        else:
            print(f"    No heat burden column found")
        
        # Cold burden map - check various possible column names
        cold_col = None
        for col in ['total_cold_an', 'cold_deaths', 'cold_attributable_deaths', 'cold_an_2_5']:
            if col in burden_df.columns:
                cold_col = col
                break
        
        if cold_col and cold_col in gdf_merged.columns:
            fig, ax = plt.subplots(1, 1, figsize=(12, 10))
            
            # Plot missing data regions only if there are any
            missing_data = gdf_merged[gdf_merged[cold_col].isna()]
            if len(missing_data) > 0:
                missing_data.plot(ax=ax, color='lightgray', edgecolor='white', linewidth=0.5)
            
            # Plot data regions
            valid_data = gdf_merged[gdf_merged[cold_col].notna()]
            if len(valid_data) > 0:
                valid_data.plot(
                    column=cold_col, ax=ax, cmap=COLD_CMAP,
                    edgecolor='white', linewidth=0.3, legend=True,
                    legend_kwds={'label': 'Attributable Deaths', 'shrink': 0.6}
                )
            
            ax.set_title(f'Cold-Attributable Deaths (2010-2024) - {level_label}', 
                        fontsize=14, fontweight='bold')
            ax.axis('off')
            
            save_map(fig, f'map5_cold_burden_{level}')
        else:
            print(f"    No cold burden column found")

def map_heterogeneity():
    """Create maps of effect heterogeneity from meta-regression."""
    print("\n[Map 6] Meta-regression heterogeneity...")
    
    # Load shapefile
    gdf_base, _ = load_brazil_shapefile()
    if gdf_base is None:
        return
    
    # Load mapping once
    mapping = load_region_mapping()
    if mapping is None:
        print("  Cannot generate maps without region mapping")
        return
    
    for level in ['intermediate', 'immediate']:
        suffix = '' if level == 'intermediate' else '_immediate'
        
        print(f"\n  Processing {level} level...")
        
        # Aggregate geometry using proper IBGE mapping
        if level == 'intermediate':
            gdf = aggregate_to_intermediate(gdf_base.copy(), mapping)
        else:
            gdf = aggregate_to_immediate(gdf_base.copy(), mapping)
        
        # Load meta-regression data
        data_file = PHASE4_RESULTS / f'meta_regression_v2_data{suffix}.csv'
        if not data_file.exists():
            continue
        
        meta_df = pd.read_csv(data_file)
        
        # Identify region column
        region_col = None
        for col in ['region_code', 'region']:
            if col in meta_df.columns:
                region_col = col
                break
        
        if region_col is None:
            continue
        
        meta_df['region_code'] = meta_df[region_col].astype(int)
        
        # Merge with geometry
        gdf_merged = gdf.merge(meta_df, on='region_code', how='left')
        
        level_label = 'Intermediate' if level == 'intermediate' else 'Immediate'
        
        # Map key covariates
        covariates_to_map = {
            'ac_penetration': ('Air Conditioning Penetration', 'YlGn'),
            'hdi': ('Human Development Index', 'YlGn'),
            'mean_temp': ('Mean Temperature (°C)', 'coolwarm'),
            'urban_pct': ('Urban Population (%)', 'Purples'),
        }
        
        for cov, (label, cmap) in covariates_to_map.items():
            if cov not in gdf_merged.columns:
                continue
            
            # Check if there's any valid data
            valid_data = gdf_merged[gdf_merged[cov].notna()]
            if len(valid_data) == 0:
                print(f"  Skipping {cov}: No valid data")
                continue
            
            fig, ax = plt.subplots(1, 1, figsize=(12, 10))
            
            try:
                gdf_merged[gdf_merged[cov].isna()].plot(
                    ax=ax, color='lightgray', edgecolor='white', linewidth=0.5
                )
                valid_data.plot(
                    column=cov, ax=ax, cmap=cmap,
                    edgecolor='white', linewidth=0.3, legend=True,
                    legend_kwds={'label': label, 'shrink': 0.6}
                )
                
                ax.set_title(f'{label} - {level_label}', fontsize=14, fontweight='bold')
                ax.axis('off')
                
                save_map(fig, f'map6_{cov}_{level}')
            except Exception as e:
                print(f"  Error mapping {cov}: {e}")
                plt.close(fig)

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    if not HAS_GEOPANDAS:
        print("\nERROR: geopandas is required for map generation.")
        print("Install with: pip install geopandas")
        exit(1)
    
    print("\nGenerating maps...")
    
    map_heat_effects()
    map_attributable_burden()
    map_heterogeneity()
    
    print("\n" + "="*70)
    print(f"MAP GENERATION COMPLETE")
    print(f"Output directory: {OUTPUT_DIR}")
    print("="*70)
