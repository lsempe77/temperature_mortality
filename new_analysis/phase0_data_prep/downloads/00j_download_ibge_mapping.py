"""
00j: Download Municipality to Region Mapping from IBGE API

Downloads municipality-to-region mapping that includes BOTH:
- Intermediate regions (133 regions - Regiões Geográficas Intermediárias)
- Immediate regions (510 regions - Regiões Geográficas Imediatas)

These are the 2017 IBGE geographic division, replacing the older
meso/micro region structure.

Source: IBGE Service Data API (https://servicodados.ibge.gov.br/api/docs/localidades)

Output:
    results/municipality_to_all_regions_map.csv
        - code_muni: Municipality IBGE code (7 digits)
        - name_muni: Municipality name
        - abbrev_state: State abbreviation (e.g., SP, RJ)
        - intermediate_code: Intermediate region code (4 digits)
        - intermediate_name: Intermediate region name
        - immediate_code: Immediate region code (6 digits)
        - immediate_name: Immediate region name
"""

import requests
import pandas as pd
from pathlib import Path


def download_ibge_mapping():
    """Download municipality to region mapping from IBGE API."""
    
    print("=" * 70)
    print("00j: DOWNLOAD IBGE MUNICIPALITY-REGION MAPPING")
    print("=" * 70)
    
    # Use the 'nivelado' view which flattens the hierarchy
    url = "https://servicodados.ibge.gov.br/api/v1/localidades/municipios?view=nivelado"
    
    print(f"\nFetching data from IBGE API...")
    print(f"URL: {url}")
    
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        data = response.json()
        
        print(f"✓ Retrieved {len(data)} municipalities")
        
        # Convert to DataFrame
        df = pd.DataFrame(data)
        
        print(f"\nColumns received from API:")
        for col in sorted(df.columns):
            print(f"  - {col}")
        
        # Check for 2017 geographic division columns
        has_intermediate = 'regiao-intermediaria-id' in df.columns
        has_immediate = 'regiao-imediata-id' in df.columns
        
        if not has_intermediate or not has_immediate:
            print("\n⚠️ Warning: 2017 geographic division not found in API response.")
            print("  Expected columns: 'regiao-intermediaria-id', 'regiao-imediata-id'")
            
            # Try alternative approach - fetch from regioes-imediatas endpoint
            print("\nTrying alternative endpoint...")
            return download_via_immediate_regions()
        
        # Extract required columns
        print("\n✓ Found 2017 geographic division (Intermediate + Immediate regions)")
        
        required_cols = [
            'municipio-id', 
            'municipio-nome', 
            'UF-sigla',
            'regiao-intermediaria-id', 
            'regiao-intermediaria-nome',
            'regiao-imediata-id',
            'regiao-imediata-nome'
        ]
        
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            print(f"⚠️ Missing columns: {missing}")
            return None
        
        # Create mapping DataFrame
        mapping = df[required_cols].copy()
        
        mapping.rename(columns={
            'municipio-id': 'code_muni',
            'municipio-nome': 'name_muni',
            'UF-sigla': 'abbrev_state',
            'regiao-intermediaria-id': 'intermediate_code',
            'regiao-intermediaria-nome': 'intermediate_name',
            'regiao-imediata-id': 'immediate_code',
            'regiao-imediata-nome': 'immediate_name'
        }, inplace=True)
        
        # Ensure codes are integers
        mapping['code_muni'] = mapping['code_muni'].astype(int)
        mapping['intermediate_code'] = mapping['intermediate_code'].astype(int)
        mapping['immediate_code'] = mapping['immediate_code'].astype(int)
        
        # Save
        output_dir = Path(__file__).resolve().parents[1] / "results"
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "municipality_to_all_regions_map.csv"
        
        mapping.to_csv(output_path, index=False)
        
        # Summary statistics
        print(f"\n" + "=" * 70)
        print("MAPPING SUMMARY")
        print("=" * 70)
        print(f"Total municipalities: {len(mapping):,}")
        print(f"Unique intermediate regions: {mapping['intermediate_code'].nunique()}")
        print(f"Unique immediate regions: {mapping['immediate_code'].nunique()}")
        print(f"States covered: {mapping['abbrev_state'].nunique()}")
        
        print(f"\nSample data:")
        print(mapping.head(10).to_string(index=False))
        
        print(f"\n✓ Saved: {output_path}")
        
        return mapping
        
    except requests.exceptions.RequestException as e:
        print(f"❌ Error fetching data: {e}")
        return None
    except Exception as e:
        print(f"❌ Error processing data: {e}")
        import traceback
        traceback.print_exc()
        return None


def download_via_immediate_regions():
    """
    Alternative approach: fetch immediate regions first, then get municipalities.
    
    Uses endpoint: /localidades/regioes-imediatas/{id}/municipios
    """
    print("\n" + "-" * 70)
    print("Using alternative approach via immediate regions endpoint")
    print("-" * 70)
    
    # First, get all immediate regions
    url_imediatas = "https://servicodados.ibge.gov.br/api/v1/localidades/regioes-imediatas"
    
    try:
        print("Fetching immediate regions...")
        response = requests.get(url_imediatas, timeout=60)
        response.raise_for_status()
        imediatas = response.json()
        
        print(f"✓ Found {len(imediatas)} immediate regions")
        
        # Build mapping by iterating through immediate regions
        all_mappings = []
        
        for i, imediata in enumerate(imediatas):
            immediate_code = imediata['id']
            immediate_name = imediata['nome']
            
            # Get parent intermediate region
            intermediaria = imediata.get('regiao-intermediaria', {})
            intermediate_code = intermediaria.get('id')
            intermediate_name = intermediaria.get('nome')
            
            # Get state
            uf = intermediaria.get('UF', {})
            abbrev_state = uf.get('sigla', '')
            
            # Get municipalities for this immediate region
            url_munis = f"https://servicodados.ibge.gov.br/api/v1/localidades/regioes-imediatas/{immediate_code}/municipios"
            
            try:
                resp_munis = requests.get(url_munis, timeout=30)
                resp_munis.raise_for_status()
                municipios = resp_munis.json()
                
                for muni in municipios:
                    all_mappings.append({
                        'code_muni': muni['id'],
                        'name_muni': muni['nome'],
                        'abbrev_state': abbrev_state,
                        'intermediate_code': intermediate_code,
                        'intermediate_name': intermediate_name,
                        'immediate_code': immediate_code,
                        'immediate_name': immediate_name
                    })
                    
            except Exception as e:
                print(f"  ⚠️ Error fetching municipalities for {immediate_name}: {e}")
                continue
            
            # Progress update every 50 regions
            if (i + 1) % 50 == 0:
                print(f"  Processed {i + 1}/{len(imediatas)} immediate regions...")
        
        if not all_mappings:
            print("❌ No mappings collected")
            return None
        
        # Create DataFrame
        mapping = pd.DataFrame(all_mappings)
        
        # Ensure integer types
        mapping['code_muni'] = mapping['code_muni'].astype(int)
        mapping['intermediate_code'] = mapping['intermediate_code'].astype(int)
        mapping['immediate_code'] = mapping['immediate_code'].astype(int)
        
        # Save
        output_dir = Path(__file__).resolve().parents[1] / "results"
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "municipality_to_all_regions_map.csv"
        
        mapping.to_csv(output_path, index=False)
        
        # Summary
        print(f"\n" + "=" * 70)
        print("MAPPING SUMMARY")
        print("=" * 70)
        print(f"Total municipalities: {len(mapping):,}")
        print(f"Unique intermediate regions: {mapping['intermediate_code'].nunique()}")
        print(f"Unique immediate regions: {mapping['immediate_code'].nunique()}")
        print(f"States covered: {mapping['abbrev_state'].nunique()}")
        
        print(f"\nSample data:")
        print(mapping.head(10).to_string(index=False))
        
        print(f"\n✓ Saved: {output_path}")
        
        return mapping
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return None


if __name__ == "__main__":
    download_ibge_mapping()
