"""
00l: PROCESS IBGE LIFE TABLES FOR YEARS OF LIFE LOST (YLL) CALCULATION
=========================================================================
Process complete life tables from IBGE (Brazilian Institute of Geography 
and Statistics) for calculating Years of Life Lost attributable to heat.

Data Source: IBGE - Tábuas Completas de Mortalidade
URL: https://www.ibge.gov.br/estatisticas/sociais/populacao/9126-tabuas-completas-de-mortalidade.html
FTP: https://ftp.ibge.gov.br/Tabuas_Completas_de_Mortalidade/

Files: Manually downloaded to Input_data/tabua{YYYY}.xls(x)

Key Variables:
- Column 0: Age (X) - Idades Exatas
- Column 1: qx - Probability of dying between ages x and x+1 (per 1000)
- Column 2: dx - Deaths between ages x and x+1
- Column 3: lx - Number of survivors to age x (from 100,000 births)
- Column 4: Lx - Person-years lived between x and x+1
- Column 5: Tx - Total person-years above age x
- Column 6: ex - Life expectancy at age x (E(X))

Output:
- results/ibge_life_tables_combined.parquet: All years combined
- results/ibge_life_expectancy_by_age.csv: Simplified ex by age for YLL
- results/yll_lookup_by_age.csv: Age-specific life expectancy lookup
- results/yll_lookup_by_age_group.csv: Age-group specific lookup

YLL Calculation Method:
- For each death at age x, YLL = ex (remaining life expectancy at that age)
- Total YLL = Σ (deaths_age_x × ex)
- Heat-attributable YLL = Total YLL × Attributable Fraction
"""

import sys
from datetime import datetime
from pathlib import Path
import warnings

print("=" * 70)
print("00l: PROCESS IBGE LIFE TABLES")
print("=" * 70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# DEPENDENCIES
# =============================================================================

try:
    import pandas as pd
    import numpy as np
    print("✓ pandas/numpy installed")
except ImportError:
    print("✗ pandas/numpy not installed")
    sys.exit(1)

try:
    import openpyxl
    print("✓ openpyxl installed (for .xlsx files)")
except ImportError:
    print("⚠️ openpyxl not installed - install with: pip install openpyxl")

try:
    import xlrd
    print("✓ xlrd installed (for .xls files)")
except ImportError:
    print("⚠️ xlrd not installed - install with: pip install xlrd")
    print("  Needed for reading older .xls files (2010-2020)")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Directories
BASE_DIR = Path(__file__).parent
PROJECT_DIR = BASE_DIR.parent.parent  # sim_data folder
INPUT_DIR = PROJECT_DIR / 'Input_data'
OUTPUT_DIR = BASE_DIR / 'results'
OUTPUT_DIR.mkdir(exist_ok=True)

# Years to process (2010-2024 to match mortality data)
YEARS = list(range(2010, 2025))

# Column indices in IBGE life tables (consistent across years)
# Based on inspection of actual files
COLUMN_INDICES = {
    'age': 0,      # Idades Exatas (X)
    'qx': 1,       # Probability of dying (per 1000)
    'dx': 2,       # Deaths
    'lx': 3,       # Survivors
    'Lx': 4,       # Person-years lived
    'Tx': 5,       # Total person-years above age
    'ex': 6        # Life expectancy at age x - E(X)
}

# Number of header rows to skip (data starts at row 6, 0-indexed)
HEADER_ROWS = 6

# =============================================================================
# PARSING FUNCTIONS
# =============================================================================

def find_life_table_file(year: int) -> Path:
    """Find the life table file for a given year."""
    
    # Try both extensions
    for ext in ['.xlsx', '.xls']:
        filepath = INPUT_DIR / f"tabua{year}{ext}"
        if filepath.exists():
            return filepath
    
    return None


def parse_life_table(filepath: Path, year: int) -> pd.DataFrame:
    """
    Parse IBGE life table Excel file.
    
    Structure (consistent across 2010-2024):
    - Rows 0-5: Headers/metadata
    - Row 6+: Data (ages 0 to 90)
    - Column 0: Age
    - Column 6: Life expectancy E(X)
    """
    
    try:
        # Read raw data, skipping header rows
        df = pd.read_excel(
            filepath,
            header=None,
            skiprows=HEADER_ROWS,
            usecols=list(COLUMN_INDICES.values())
        )
        
        # Rename columns
        col_names = list(COLUMN_INDICES.keys())
        df.columns = col_names
        
        # Clean data
        # Age should be numeric (0-90)
        df['age'] = pd.to_numeric(df['age'], errors='coerce')
        df = df.dropna(subset=['age'])
        df['age'] = df['age'].astype(int)
        
        # Keep only valid ages (0-90)
        df = df[(df['age'] >= 0) & (df['age'] <= 90)].copy()
        
        # Convert numeric columns
        for col in ['qx', 'dx', 'lx', 'Lx', 'Tx', 'ex']:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Add year
        df['year'] = year
        
        # Sort by age
        df = df.sort_values('age').reset_index(drop=True)
        
        return df
        
    except Exception as e:
        print(f"    ✗ Parse error: {str(e)[:80]}")
        return None


def create_yll_lookup(combined_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create simplified lookup table for YLL calculation.
    
    Returns: age, ex (life expectancy at that age), averaged across years
    """
    
    yll_lookup = combined_df.groupby('age').agg(
        ex_mean=('ex', 'mean'),
        ex_min=('ex', 'min'),
        ex_max=('ex', 'max'),
        n_years=('year', 'nunique')
    ).reset_index()
    
    return yll_lookup


def create_age_group_yll(combined_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create YLL lookup for common mortality age groups.
    
    Age groups typically used in mortality data:
    - 0, 1-4, 5-9, ..., 75-79, 80+
    """
    
    # Define age groups
    age_bins = [0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 200]
    age_labels = ['0', '1-4', '5-9', '10-14', '15-19', '20-24', '25-29', '30-34',
                  '35-39', '40-44', '45-49', '50-54', '55-59', '60-64', '65-69',
                  '70-74', '75-79', '80+']
    
    df = combined_df.copy()
    df['age_group'] = pd.cut(df['age'], bins=age_bins, labels=age_labels, right=False)
    
    # Average life expectancy within each age group
    grouped = df.groupby(['age_group', 'year']).agg(
        ex_mean=('ex', 'mean'),
        age_midpoint=('age', 'median')
    ).reset_index()
    
    # Then average across years
    final = grouped.groupby('age_group').agg(
        ex=('ex_mean', 'mean'),
        age_midpoint=('age_midpoint', 'first')
    ).reset_index()
    
    return final


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Process all IBGE life tables from local files."""
    
    print(f"\nInput directory: {INPUT_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    print("\n" + "-" * 70)
    print("STEP 1: Find and parse IBGE life tables")
    print("-" * 70)
    
    all_data = []
    years_found = []
    years_missing = []
    
    for year in YEARS:
        filepath = find_life_table_file(year)
        
        if filepath is None:
            print(f"  ⚠️ {year}: File not found (tabua{year}.xls or .xlsx)")
            years_missing.append(year)
            continue
        
        print(f"  → {year}: Found {filepath.name}")
        
        df = parse_life_table(filepath, year)
        
        if df is not None and len(df) > 0:
            all_data.append(df)
            years_found.append(year)
            
            # Quick summary - life expectancy at birth
            e0 = df[df['age'] == 0]['ex'].values
            if len(e0) > 0:
                print(f"    ✓ {len(df)} ages, E(0) = {e0[0]:.2f} years")
        else:
            print(f"    ✗ Failed to parse")
            years_missing.append(year)
    
    if not all_data:
        print("\n✗ No data parsed successfully!")
        return
    
    # Summary
    print(f"\n  Found: {len(years_found)} years")
    if years_missing:
        print(f"  Missing: {years_missing}")
    
    # Combine all years
    print("\n" + "-" * 70)
    print("STEP 2: Combine and create lookup tables")
    print("-" * 70)
    
    combined = pd.concat(all_data, ignore_index=True)
    print(f"\nTotal records: {len(combined):,}")
    print(f"Years: {combined['year'].min()} to {combined['year'].max()}")
    print(f"Ages: {combined['age'].min()} to {combined['age'].max()}")
    
    # Life expectancy at birth trend
    print("\nLife expectancy at birth by year:")
    e0_by_year = combined[combined['age'] == 0][['year', 'ex']].sort_values('year')
    for _, row in e0_by_year.iterrows():
        print(f"  {int(row['year'])}: {row['ex']:.2f} years")
    
    # Create outputs
    print("\n" + "-" * 70)
    print("STEP 3: Save output files")
    print("-" * 70)
    
    # Full combined data (Parquet)
    combined_parquet = OUTPUT_DIR / 'ibge_life_tables_combined.parquet'
    combined.to_parquet(combined_parquet, index=False)
    print(f"✓ Saved: {combined_parquet.name}")
    
    # Full combined data (CSV for easy inspection)
    combined_csv = OUTPUT_DIR / 'ibge_life_tables_combined.csv'
    combined.to_csv(combined_csv, index=False)
    print(f"✓ Saved: {combined_csv.name}")
    
    # Create YLL lookup table (by single year of age)
    yll_lookup = create_yll_lookup(combined)
    yll_file = OUTPUT_DIR / 'yll_lookup_by_age.csv'
    yll_lookup.to_csv(yll_file, index=False)
    print(f"✓ Saved: {yll_file.name}")
    
    # Create age group lookup
    age_group_yll = create_age_group_yll(combined)
    age_group_file = OUTPUT_DIR / 'yll_lookup_by_age_group.csv'
    age_group_yll.to_csv(age_group_file, index=False)
    print(f"✓ Saved: {age_group_file.name}")
    
    # Summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    print(f"\nYears covered: {sorted(combined['year'].unique())}")
    print(f"Age range: {int(combined['age'].min())} to {int(combined['age'].max())}")
    
    # Life expectancy at key ages (for elderly analysis)
    print("\nLife expectancy at key ages (average across all years):")
    for age in [60, 65, 70, 75, 80, 85, 90]:
        ex_vals = combined[combined['age'] == age]['ex']
        if len(ex_vals) > 0:
            print(f"  Age {age}: {ex_vals.mean():.2f} remaining years")
    
    # Change over time
    print("\nLife expectancy trend at age 65:")
    e65_by_year = combined[combined['age'] == 65][['year', 'ex']].sort_values('year')
    e65_start = e65_by_year.iloc[0]['ex']
    e65_end = e65_by_year.iloc[-1]['ex']
    print(f"  {int(e65_by_year.iloc[0]['year'])}: {e65_start:.2f} years")
    print(f"  {int(e65_by_year.iloc[-1]['year'])}: {e65_end:.2f} years")
    print(f"  Change: {e65_end - e65_start:+.2f} years")
    
    # Example YLL calculation
    print("\n" + "-" * 70)
    print("YLL CALCULATION EXAMPLE")
    print("-" * 70)
    
    # Using 2024 data (most recent)
    latest = combined[combined['year'] == combined['year'].max()]
    
    print(f"\nUsing {int(combined['year'].max())} life table:")
    print("  If 100 people die at age 70:")
    e70 = latest[latest['age'] == 70]['ex'].values[0]
    yll = 100 * e70
    print(f"    Life expectancy at 70: {e70:.2f} years")
    print(f"    YLL = 100 × {e70:.2f} = {yll:.0f} years lost")
    
    print(f"\n✓ Done! {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    return combined


if __name__ == "__main__":
    main()
