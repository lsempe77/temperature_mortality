"""
Check spatial data for aggregation levels.
Inspects the municipality shapefile to see if we can aggregate to:
1. Intermediate Regions (Regiões Intermediárias) - Preferred
2. Micro-regions (Microrregiões) - Alternative
3. Immediate Regions (Regiões Imediatas)

And creates a mapping file for the analysis.
"""

import geopandas as gpd
import pandas as pd
import os
from pathlib import Path

def check_spatial_data():
    # Paths
    base_dir = Path(__file__).resolve().parents[2]
    input_path = base_dir / "Input_data" / "brazil_municipalities_2022.gpkg"
    output_path = base_dir / "new_analysis" / "phase0_data_prep" / "results" / "municipality_to_region_map.csv"
    
    print(f"Inspecting: {input_path}")
    
    try:
        # Load data (just the attribute table, no geometry needed for mapping)
        gdf = gpd.read_file(input_path, ignore_geometry=True)
        
        print("\nColumns found:")
        print(gdf.columns.tolist())
        
        # Check for aggregation columns (Standard IBGE codes)
        # CD_RGINT / NM_RGINT = Intermediate Region
        # CD_MICRO / NM_MICRO = Microregion
        # CD_RGI / NM_RGI = Immediate Region
        
        has_intermediate = 'CD_RGINT' in gdf.columns
        has_micro = 'CD_MICRO' in gdf.columns
        has_immediate = 'CD_RGI' in gdf.columns
        
        print("\nAggregation Levels Available:")
        print(f"- Intermediate Regions: {has_intermediate}")
        print(f"- Micro-regions: {has_micro}")
        print(f"- Immediate Regions: {has_immediate}")
        
        # Determine best aggregation level
        target_col = None
        target_name_col = None
        level_name = None
        
        if has_intermediate:
            target_col = 'CD_RGINT'
            target_name_col = 'NM_RGINT'
            level_name = 'intermediate_region'
        elif has_micro:
            target_col = 'CD_MICRO'
            target_name_col = 'NM_MICRO'
            level_name = 'micro_region'
        elif has_immediate:
            target_col = 'CD_RGI'
            target_name_col = 'NM_RGI'
            level_name = 'immediate_region'
            
        if target_col:
            n_units = gdf[target_col].nunique()
            print(f"\nRecommended Aggregation: {level_name.upper()}")
            print(f"Number of units: {n_units} (vs 27 States)")
            print("This will significantly reduce Berkson error while maintaining statistical power.")
            
            # Create mapping file
            # We need: code_muni -> code_region
            # Assuming CD_MUN is the municipality code
            
            if 'CD_MUN' in gdf.columns:
                mapping = gdf[['CD_MUN', 'NM_MUN', 'SIGLA_UF', target_col, target_name_col]].copy()
                mapping.rename(columns={
                    target_col: 'region_code', 
                    target_name_col: 'region_name'
                }, inplace=True)
                
                mapping['aggregation_level'] = level_name
                
                # Ensure output directory exists
                output_path.parent.mkdir(parents=True, exist_ok=True)
                
                mapping.to_csv(output_path, index=False)
                print(f"\nMapping saved to: {output_path}")
                print(mapping.head())
            else:
                print("Error: CD_MUN column not found. Cannot create mapping.")
                
        else:
            print("\nWARNING: No suitable aggregation columns found (RGINT, MICRO, or RGI).")
            print("We may need to download a richer shapefile from IBGE.")

    except Exception as e:
        print(f"\nError reading file: {e}")
        print("Make sure 'geopandas' is installed: pip install geopandas")

if __name__ == "__main__":
    check_spatial_data()
