"""
00a: DOCUMENT EXISTING DATA
============================
Inspect all available data files and create a comprehensive data dictionary.

This script:
1. Lists all parquet files in processed_data/
2. Documents columns, types, ranges, missing values
3. Creates a data dictionary for the methods section
4. Identifies what additional data we need

Output: results/data_dictionary.md
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime
from pathlib import Path

# =============================================================================
# CONFIGURATION
# =============================================================================

PROCESSED_DATA_DIR = '../../processed_data'
INPUT_DATA_DIR = '../../Input_data'
RESULTS_DIR = './results'

print("="*70)
print("00a: DOCUMENT EXISTING DATA")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def describe_dataframe(df, name):
    """Generate detailed description of a dataframe."""
    
    info = {
        'name': name,
        'rows': len(df),
        'columns': len(df.columns),
        'memory_mb': df.memory_usage(deep=True).sum() / 1e6,
        'columns_info': []
    }
    
    for col in df.columns:
        col_info = {
            'name': col,
            'dtype': str(df[col].dtype),
            'non_null': df[col].notna().sum(),
            'null_pct': df[col].isna().mean() * 100,
            'unique': df[col].nunique()
        }
        
        # Add range info for numeric columns
        if pd.api.types.is_numeric_dtype(df[col]):
            col_info['min'] = df[col].min()
            col_info['max'] = df[col].max()
            col_info['mean'] = df[col].mean()
            col_info['median'] = df[col].median()
        
        # Add sample values for categorical/object columns
        elif df[col].dtype == 'object' or df[col].dtype.name == 'category':
            sample_vals = df[col].dropna().unique()[:5]
            col_info['sample_values'] = list(sample_vals)
        
        # Add date range for datetime columns
        elif pd.api.types.is_datetime64_any_dtype(df[col]):
            col_info['min_date'] = str(df[col].min())
            col_info['max_date'] = str(df[col].max())
        
        info['columns_info'].append(col_info)
    
    return info

def format_markdown_table(info):
    """Format column info as markdown table."""
    
    lines = []
    lines.append(f"### {info['name']}")
    lines.append(f"- **Rows**: {info['rows']:,}")
    lines.append(f"- **Columns**: {info['columns']}")
    lines.append(f"- **Memory**: {info['memory_mb']:.1f} MB")
    lines.append("")
    lines.append("| Column | Type | Non-Null | Null% | Details |")
    lines.append("|--------|------|----------|-------|---------|")
    
    for col in info['columns_info']:
        details = ""
        if 'min' in col:
            details = f"Range: [{col['min']:.2f}, {col['max']:.2f}], Mean: {col['mean']:.2f}"
        elif 'sample_values' in col:
            vals = [str(v)[:20] for v in col['sample_values'][:3]]
            details = f"e.g., {', '.join(vals)}"
        elif 'min_date' in col:
            details = f"{col['min_date']} to {col['max_date']}"
        
        lines.append(f"| {col['name']} | {col['dtype']} | {col['non_null']:,} | {col['null_pct']:.1f}% | {details} |")
    
    lines.append("")
    return "\n".join(lines)

# =============================================================================
# INSPECT PARQUET FILES
# =============================================================================

print("\n" + "-"*70)
print("INSPECTING PARQUET FILES")
print("-"*70)

parquet_files = list(Path(PROCESSED_DATA_DIR).glob("*.parquet"))
print(f"Found {len(parquet_files)} parquet files")

all_info = []
markdown_output = []
markdown_output.append("# Data Dictionary")
markdown_output.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
markdown_output.append("---\n")
markdown_output.append("## Processed Data Files\n")

for pq_file in parquet_files:
    print(f"\n  Loading {pq_file.name}...")
    df = pd.read_parquet(pq_file)
    info = describe_dataframe(df, pq_file.name)
    all_info.append(info)
    markdown_output.append(format_markdown_table(info))
    print(f"    → {info['rows']:,} rows × {info['columns']} columns")

# =============================================================================
# INSPECT ERA5 NETCDF FILES
# =============================================================================

print("\n" + "-"*70)
print("INSPECTING ERA5 NETCDF FILES")
print("-"*70)

markdown_output.append("## ERA5 Climate Data Files\n")

nc_files = list(Path(INPUT_DATA_DIR).glob("era5_brazil_*.nc"))
print(f"Found {len(nc_files)} NetCDF files")

try:
    import xarray as xr
    
    for nc_file in sorted(nc_files):
        print(f"\n  {nc_file.name}:")
        ds = xr.open_dataset(nc_file)
        
        lines = []
        lines.append(f"### {nc_file.name}")
        lines.append(f"- **Variables**: {list(ds.data_vars)}")
        lines.append(f"- **Dimensions**: {dict(ds.dims)}")
        
        # Time range
        if 'time' in ds.dims:
            time_min = pd.Timestamp(ds.time.values.min())
            time_max = pd.Timestamp(ds.time.values.max())
            lines.append(f"- **Time range**: {time_min.strftime('%Y-%m-%d')} to {time_max.strftime('%Y-%m-%d')}")
        
        # Spatial extent
        if 'latitude' in ds.dims:
            lines.append(f"- **Latitude**: {float(ds.latitude.min()):.2f}° to {float(ds.latitude.max()):.2f}°")
            lines.append(f"- **Longitude**: {float(ds.longitude.min()):.2f}° to {float(ds.longitude.max()):.2f}°")
        
        lines.append("")
        markdown_output.append("\n".join(lines))
        
        ds.close()
        print(f"    Variables: {list(ds.data_vars)}")

except ImportError:
    print("  xarray not available - skipping NetCDF inspection")
    markdown_output.append("*xarray not available for NetCDF inspection*\n")

# =============================================================================
# INSPECT RAW MORTALITY FILES
# =============================================================================

print("\n" + "-"*70)
print("INSPECTING RAW MORTALITY FILES (first 1000 rows)")
print("-"*70)

markdown_output.append("## Raw Mortality Data (DO files)\n")

csv_files = list(Path(INPUT_DATA_DIR).glob("DO*.csv"))
print(f"Found {len(csv_files)} mortality CSV files")

for csv_file in sorted(csv_files):
    print(f"\n  {csv_file.name}:")
    
    # Try different encodings and delimiters
    try:
        # Try semicolon first (common in Brazilian data)
        df = pd.read_csv(csv_file, nrows=1000, sep=';', encoding='latin-1', low_memory=False)
        if len(df.columns) < 3:
            df = pd.read_csv(csv_file, nrows=1000, sep=',', encoding='latin-1', low_memory=False)
    except:
        try:
            df = pd.read_csv(csv_file, nrows=1000, encoding='utf-8', low_memory=False)
        except Exception as e:
            print(f"    Error reading file: {e}")
            continue
    
    # Count total rows
    with open(csv_file, 'r', encoding='latin-1', errors='ignore') as f:
        total_rows = sum(1 for _ in f) - 1  # Subtract header
    
    print(f"    → {total_rows:,} rows, {len(df.columns)} columns")
    
    lines = []
    lines.append(f"### {csv_file.name}")
    lines.append(f"- **Total rows**: {total_rows:,}")
    lines.append(f"- **Columns**: {len(df.columns)}")
    lines.append("")
    lines.append("**Key columns** (sample from first 1000 rows):")
    lines.append("")
    
    # Identify key columns
    key_cols = ['DTOBITO', 'CODMUNRES', 'IDADE', 'SEXO', 'RACACOR', 'ESC', 'CAUSABAS', 'LOCOCOR']
    for col in key_cols:
        if col in df.columns:
            unique = df[col].nunique()
            sample = df[col].dropna().iloc[:3].tolist() if len(df[col].dropna()) > 0 else []
            lines.append(f"- `{col}`: {unique} unique values (e.g., {sample[:2]})")
    
    lines.append("")
    markdown_output.append("\n".join(lines))

# =============================================================================
# IDENTIFY MISSING DATA NEEDS
# =============================================================================

print("\n" + "-"*70)
print("DATA GAPS ANALYSIS")
print("-"*70)

markdown_output.append("## Data Gaps & Needs\n")
markdown_output.append("### Currently Available\n")
markdown_output.append("- ✅ Daily mortality by municipality (2022-2024)\n")
markdown_output.append("- ✅ ERA5 temperature data (hourly → aggregated to daily)\n")
markdown_output.append("- ✅ Demographic stratification (age, sex, race, education)\n")
markdown_output.append("- ✅ Cause of death (ICD-10 major categories)\n")
markdown_output.append("- ✅ Place of death (hospital vs home)\n")
markdown_output.append("\n")

markdown_output.append("### Need to Create/Download\n")
markdown_output.append("| Data | Source | Priority | Notes |\n")
markdown_output.append("|------|--------|----------|-------|\n")
markdown_output.append("| Humidity/Apparent Temp | Calculate from ERA5 dewpoint | HIGH | Already have dewpoint data |\n")
markdown_output.append("| PM2.5 daily | CAMS reanalysis | MEDIUM | Copernicus ADS |\n")
markdown_output.append("| O3 daily | CAMS reanalysis | MEDIUM | Copernicus ADS |\n")
markdown_output.append("| INMET station data | INMET API | LOW | For ERA5 validation |\n")
markdown_output.append("| Brazilian holidays | Manual | LOW | ~15 national holidays/year |\n")
markdown_output.append("| AC ownership by state | PNAD 2022-2023 | MEDIUM | For meta-regression |\n")
markdown_output.append("| State covariates | IBGE SIDRA | MEDIUM | Urbanization, GDP, etc. |\n")

# =============================================================================
# CHECK FOR HUMIDITY DATA IN ERA5
# =============================================================================

print("\nChecking for dewpoint temperature in ERA5 (for humidity calculation)...")

has_dewpoint = False
for nc_file in nc_files:
    try:
        ds = xr.open_dataset(nc_file)
        if 'd2m' in ds.data_vars or 'dewpoint' in str(ds.data_vars).lower():
            has_dewpoint = True
            print(f"  ✓ Found dewpoint in {nc_file.name}")
        ds.close()
    except:
        pass

if has_dewpoint:
    print("  → Can calculate relative humidity and apparent temperature from existing data!")
    markdown_output.append("\n### Humidity Status\n")
    markdown_output.append("✅ Dewpoint temperature available in ERA5 files - can calculate humidity without new downloads.\n")
else:
    print("  → May need to download dewpoint data separately")
    markdown_output.append("\n### Humidity Status\n")
    markdown_output.append("⚠️ Dewpoint not found in current ERA5 files - may need additional download.\n")

# =============================================================================
# SAVE DATA DICTIONARY
# =============================================================================

output_file = f'{RESULTS_DIR}/data_dictionary.md'
with open(output_file, 'w', encoding='utf-8') as f:
    f.write("\n".join(markdown_output))

print(f"\n{'='*70}")
print(f"Data dictionary saved to: {output_file}")
print(f"{'='*70}")

# Also save a summary JSON
import json

summary = {
    'generated': datetime.now().isoformat(),
    'parquet_files': [info['name'] for info in all_info],
    'total_parquet_rows': sum(info['rows'] for info in all_info),
    'era5_files': [f.name for f in nc_files],
    'mortality_files': [f.name for f in csv_files],
    'has_dewpoint': has_dewpoint
}

with open(f'{RESULTS_DIR}/data_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
