"""
00a: COMPREHENSIVE DATA DOCUMENTATION
======================================
Generate complete documentation for all Phase 0 aggregated datasets.

Outputs:
1. data_dictionary.md - Column-by-column documentation for all datasets
2. schema_specification.md - Unified schema for appendix publication  
3. validation_report.md - Completeness, missingness, alignment checks
4. metadata.json - Structured metadata (rows, dates, regions, types)

Run after any Phase 0 data updates to regenerate documentation.
"""

import pandas as pd
import numpy as np
import json
from datetime import datetime
from pathlib import Path
from collections import OrderedDict

print("=" * 80)
print("00a: COMPREHENSIVE DATA DOCUMENTATION")
print("=" * 80)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = Path(__file__).parent
RESULTS_DIR = BASE_DIR / 'results'
PHASE0_RESULTS = BASE_DIR.parent / 'results'
OUTPUT_DIR = RESULTS_DIR

# Ensure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Define all Phase 0 datasets with metadata
DATASETS = {
    # === EXPOSURE DATA ===
    'era5_intermediate_daily': {
        'file': PHASE0_RESULTS / 'era5_intermediate_daily.parquet',
        'category': 'Exposure',
        'description': 'ERA5 temperature data aggregated to 133 intermediate regions',
        'spatial_unit': 'intermediate',
        'temporal_unit': 'daily',
        'key_columns': {
            'region_code': 'Intermediate region identifier (IBGE code)',
            'date': 'Date of observation',
            'temp_mean': 'Daily mean temperature (°C)',
            'temp_min': 'Daily minimum temperature (°C)',
            'temp_max': 'Daily maximum temperature (°C)',
            'dewpoint_mean': 'Daily mean dewpoint temperature (°C)',
        }
    },
    'era5_immediate_daily': {
        'file': PHASE0_RESULTS / 'era5_immediate_daily.parquet',
        'category': 'Exposure',
        'description': 'ERA5 temperature data aggregated to 510 immediate regions',
        'spatial_unit': 'immediate',
        'temporal_unit': 'daily',
        'key_columns': {
            'region_code': 'Immediate region identifier (IBGE code)',  # Note: uses region_code not immediate_code
            'date': 'Date of observation',
            'temp_mean': 'Daily mean temperature (°C)',
            'temp_min': 'Daily minimum temperature (°C)',
            'temp_max': 'Daily maximum temperature (°C)',
            'dewpoint_mean': 'Daily mean dewpoint temperature (°C)',
        }
    },
    
    # === POLLUTION DATA ===
    'cams_intermediate_daily': {
        'file': PHASE0_RESULTS / 'cams_intermediate_daily.parquet',
        'category': 'Exposure',
        'description': 'CAMS pollution data aggregated to 133 intermediate regions',
        'spatial_unit': 'intermediate',
        'temporal_unit': 'daily',
        'key_columns': {
            'intermediate_code': 'Intermediate region identifier',
            'date': 'Date of observation',
            'pm25': 'Daily mean PM2.5 concentration (μg/m³)',
            'ozone': 'Daily mean O3 concentration (μg/m³)',
        }
    },
    'cams_immediate_daily': {
        'file': PHASE0_RESULTS / 'cams_immediate_daily.parquet',
        'category': 'Exposure',
        'description': 'CAMS pollution data aggregated to 510 immediate regions',
        'spatial_unit': 'immediate',
        'temporal_unit': 'daily',
        'key_columns': {
            'immediate_code': 'Immediate region identifier',
            'date': 'Date of observation',
            'pm25': 'Daily mean PM2.5 concentration (μg/m³)',
            'ozone': 'Daily mean O3 concentration (μg/m³)',
        }
    },
    
    # === MORTALITY DATA ===
    'mortality_intermediate_daily': {
        'file': PHASE0_RESULTS / 'mortality_intermediate_daily.parquet',
        'category': 'Outcome',
        'description': 'All-ages mortality aggregated to 133 intermediate regions',
        'spatial_unit': 'intermediate',
        'temporal_unit': 'daily',
        'key_columns': {
            'intermediate_code': 'Intermediate region identifier',
            'date': 'Date of death',
            'deaths_all': 'Total deaths (all ages)',
            'deaths_respiratory': 'Respiratory deaths (ICD J00-J99)',
            'deaths_cardiovascular': 'Cardiovascular deaths (ICD I00-I99)',
            'deaths_heat_direct': 'Direct heat deaths (ICD T67, X30)',
        }
    },
    'mortality_intermediate_daily_elderly': {
        'file': PHASE0_RESULTS / 'mortality_intermediate_daily_elderly.parquet',
        'category': 'Outcome',
        'description': 'Elderly (60+) mortality aggregated to 133 intermediate regions',
        'spatial_unit': 'intermediate',
        'temporal_unit': 'daily',
        'key_columns': {
            'intermediate_code': 'Intermediate region identifier',
            'date': 'Date of death',
            'deaths_elderly': 'Total elderly deaths (60+)',
            'deaths_elderly_resp': 'Elderly respiratory deaths',
            'deaths_elderly_cvd': 'Elderly cardiovascular deaths',
            'deaths_elderly_heat': 'Elderly direct heat deaths',
        }
    },
    'mortality_immediate_daily': {
        'file': PHASE0_RESULTS / 'mortality_immediate_daily.parquet',
        'category': 'Outcome',
        'description': 'All-ages mortality aggregated to 510 immediate regions',
        'spatial_unit': 'immediate',
        'temporal_unit': 'daily',
        'key_columns': {
            'immediate_code': 'Immediate region identifier',
            'date': 'Date of death',
            'deaths_all': 'Total deaths (all ages)',
            'deaths_respiratory': 'Respiratory deaths (ICD J00-J99)',
            'deaths_cardiovascular': 'Cardiovascular deaths (ICD I00-I99)',
            'deaths_heat_direct': 'Direct heat deaths (ICD T67, X30)',
        }
    },
    'mortality_immediate_daily_elderly': {
        'file': PHASE0_RESULTS / 'mortality_immediate_daily_elderly.parquet',
        'category': 'Outcome',
        'description': 'Elderly (60+) mortality aggregated to 510 immediate regions',
        'spatial_unit': 'immediate',
        'temporal_unit': 'daily',
        'key_columns': {
            'immediate_code': 'Immediate region identifier',
            'date': 'Date of death',
            'deaths_elderly': 'Total elderly deaths (60+)',
            'deaths_elderly_resp': 'Elderly respiratory deaths',
            'deaths_elderly_cvd': 'Elderly cardiovascular deaths',
            'deaths_elderly_heat': 'Elderly direct heat deaths',
        }
    },
    # Legacy file names (for backwards compatibility check)
    'mortality_regional_daily': {
        'file': PHASE0_RESULTS / 'mortality_regional_daily.parquet',
        'category': 'Outcome',
        'description': '[LEGACY] All-ages mortality - intermediate regions (old naming)',
        'spatial_unit': 'intermediate',
        'temporal_unit': 'daily',
        'key_columns': {
            'region_code': 'Intermediate region identifier',
            'date': 'Date of death',
            'deaths_all': 'Total deaths (all ages)',
        }
    },
    'mortality_regional_daily_elderly': {
        'file': PHASE0_RESULTS / 'mortality_regional_daily_elderly.parquet',
        'category': 'Outcome',
        'description': '[LEGACY] Elderly mortality - intermediate regions (old naming)',
        'spatial_unit': 'intermediate',
        'temporal_unit': 'daily',
        'key_columns': {
            'region_code': 'Intermediate region identifier',
            'date': 'Date of death',
            'deaths_elderly': 'Total elderly deaths (60+)',
        }
    },
    
    # === CONFOUNDER DATA ===
    'influenza_daily_by_intermediate_region': {
        'file': PHASE0_RESULTS / 'influenza_daily_by_intermediate_region.parquet',
        'category': 'Confounder',
        'description': 'SRAG/Influenza cases aggregated to 133 intermediate regions',
        'spatial_unit': 'intermediate',
        'temporal_unit': 'daily',
        'key_columns': {
            'intermediate_code': 'Intermediate region identifier',
            'date': 'Date of notification',
            'srag_cases': 'Total SRAG cases (elderly)',
            'srag_influenza': 'Confirmed influenza cases',
            'srag_covid': 'Confirmed COVID-19 cases',
            'srag_deaths': 'SRAG-related deaths',
        }
    },
    'influenza_daily_by_immediate_region': {
        'file': PHASE0_RESULTS / 'influenza_daily_by_immediate_region.parquet',
        'category': 'Confounder',
        'description': 'SRAG/Influenza cases aggregated to 510 immediate regions',
        'spatial_unit': 'immediate',
        'temporal_unit': 'daily',
        'key_columns': {
            'immediate_code': 'Immediate region identifier',
            'date': 'Date of notification',
            'srag_cases': 'Total SRAG cases (elderly)',
            'srag_influenza': 'Confirmed influenza cases',
            'srag_covid': 'Confirmed COVID-19 cases',
            'srag_deaths': 'SRAG-related deaths',
        }
    },
    'brazilian_holidays_daily': {
        'file': PHASE0_RESULTS / 'brazilian_holidays_daily.parquet',
        'category': 'Confounder',
        'description': 'Brazilian national holidays 2010-2024',
        'spatial_unit': 'national',
        'temporal_unit': 'daily',
        'key_columns': {
            'date': 'Date',
            'is_holiday': 'Binary indicator for national holiday',
            'holiday_name_pt': 'Name of holiday in Portuguese (if applicable)',
            'holiday_type': 'Type of holiday',
        },
        # Columns that are expected to be sparse (only filled on holiday days)
        'sparse_columns': ['holiday_name_pt', 'holiday_type']
    },
    
    # === COVARIATE DATA ===
    'ses_intermediate_covariates': {
        'file': PHASE0_RESULTS / 'ses_intermediate_covariates.csv',
        'category': 'Covariate',
        'description': 'Socioeconomic covariates for 133 intermediate regions',
        'spatial_unit': 'intermediate',
        'temporal_unit': 'cross-sectional',
        'key_columns': {
            'intermediate_code': 'Intermediate region identifier',
            'intermediate_name': 'Region name',
            'pop_total': 'Total population',
            'pop_elderly': 'Elderly population (60+)',
            'pct_elderly': 'Percentage elderly',
            'gdp_total': 'Total GDP (R$)',
            'pop_urban': 'Urban population',
        }
    },
    'ses_immediate_covariates': {
        'file': PHASE0_RESULTS / 'ses_immediate_covariates.csv',
        'category': 'Covariate',
        'description': 'Socioeconomic covariates for 510 immediate regions',
        'spatial_unit': 'immediate',
        'temporal_unit': 'cross-sectional',
        'key_columns': {
            'immediate_code': 'Immediate region identifier',
            'immediate_name': 'Region name',
            'pop_total': 'Total population',
            'pop_elderly': 'Elderly population (60+)',
            'pct_elderly': 'Percentage elderly',
            'gdp_total': 'Total GDP (R$)',
            'pop_urban': 'Urban population',
        }
    },
    'regional_covariates': {
        'file': PHASE0_RESULTS / 'regional_covariates.csv',
        'category': 'Covariate',
        'description': 'Additional regional covariates for meta-regression',
        'spatial_unit': 'intermediate',
        'temporal_unit': 'cross-sectional',
        'key_columns': {
            'region_code': 'Intermediate region identifier',
            'region_name': 'Region name',
            'ac_pct': 'Air conditioning ownership rate (%)',
            'hdi': 'Human Development Index',
            'gdp_per_capita_brl': 'GDP per capita (R$)',
        }
    },
    
    # === MAPPING FILES ===
    'municipality_to_all_regions_map': {
        'file': PHASE0_RESULTS / 'municipality_to_all_regions_map.csv',
        'category': 'Reference',
        'description': 'Municipality to intermediate/immediate region mapping',
        'spatial_unit': 'municipality',
        'temporal_unit': 'cross-sectional',
        'key_columns': {
            'code_muni': 'Municipality code (7-digit IBGE)',
            'name_muni': 'Municipality name',
            'abbrev_state': 'State abbreviation',
            'intermediate_code': 'Intermediate region code',
            'intermediate_name': 'Intermediate region name',
            'immediate_code': 'Immediate region code',
            'immediate_name': 'Immediate region name',
        }
    },
    
    # === LIFE TABLES ===
    'ibge_life_tables_combined': {
        'file': PHASE0_RESULTS / 'ibge_life_tables_combined.parquet',
        'category': 'Reference',
        'description': 'IBGE life tables for YLL calculation',
        'spatial_unit': 'national',
        'temporal_unit': 'annual',
        'key_columns': {
            'year': 'Reference year',
            'age': 'Age in years',
            'ex': 'Remaining life expectancy',
            'lx': 'Number surviving to age x',
            'qx': 'Probability of dying between age x and x+1',
        }
    },
    'yll_lookup_by_age': {
        'file': PHASE0_RESULTS / 'yll_lookup_by_age.csv',
        'category': 'Reference',
        'description': 'YLL lookup table by single year of age',
        'spatial_unit': 'national',
        'temporal_unit': 'cross-sectional',
        'key_columns': {
            'age': 'Age in years (0-89)',
            'ex_mean': 'Mean remaining life expectancy across years',
            'ex_min': 'Minimum life expectancy',
            'ex_max': 'Maximum life expectancy',
        }
    },
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_dataset(filepath):
    """Load dataset from parquet or CSV."""
    if not filepath.exists():
        return None
    
    try:
        if filepath.suffix == '.parquet':
            return pd.read_parquet(filepath)
        elif filepath.suffix == '.csv':
            return pd.read_csv(filepath)
    except Exception as e:
        print(f"    Error loading {filepath.name}: {e}")
        return None
    return None


def get_column_stats(df, col):
    """Get comprehensive statistics for a column."""
    stats = {
        'name': col,
        'dtype': str(df[col].dtype),
        'non_null': int(df[col].notna().sum()),
        'null_count': int(df[col].isna().sum()),
        'null_pct': round(df[col].isna().mean() * 100, 2),
        'unique': int(df[col].nunique()),
    }
    
    # Numeric columns
    if pd.api.types.is_numeric_dtype(df[col]):
        valid = df[col].dropna()
        if len(valid) > 0:
            stats['min'] = float(valid.min())
            stats['max'] = float(valid.max())
            stats['mean'] = round(float(valid.mean()), 4)
            stats['median'] = float(valid.median())
            stats['std'] = round(float(valid.std()), 4)
            stats['p25'] = float(valid.quantile(0.25))
            stats['p75'] = float(valid.quantile(0.75))
    
    # Datetime columns
    elif pd.api.types.is_datetime64_any_dtype(df[col]):
        valid = df[col].dropna()
        if len(valid) > 0:
            stats['min_date'] = str(valid.min().date())
            stats['max_date'] = str(valid.max().date())
            stats['date_range_days'] = (valid.max() - valid.min()).days
    
    # Categorical/object columns
    elif df[col].dtype == 'object' or df[col].dtype.name == 'category':
        sample = df[col].dropna().unique()[:5]
        stats['sample_values'] = [str(v)[:50] for v in sample]
    
    return stats


def validate_dataset(df, name, meta):
    """Validate a dataset for completeness and quality."""
    issues = []
    warnings = []
    
    # Check for empty dataset
    if len(df) == 0:
        issues.append("Dataset is empty")
        return {'issues': issues, 'warnings': warnings, 'passed': False}
    
    # Check for expected columns
    expected_cols = list(meta.get('key_columns', {}).keys())
    missing_cols = [c for c in expected_cols if c not in df.columns]
    if missing_cols:
        issues.append(f"Missing expected columns: {missing_cols}")
    
    # Get columns expected to be sparse (e.g., holiday names only on ~3.6% of days)
    sparse_columns = meta.get('sparse_columns', [])
    
    # Check for high missingness (excluding expected sparse columns)
    for col in df.columns:
        miss_pct = df[col].isna().mean() * 100
        if col in sparse_columns:
            # Sparse columns are expected to have high missingness - just info
            continue
        elif miss_pct > 50:
            issues.append(f"Column '{col}' has {miss_pct:.1f}% missing values")
        elif miss_pct > 10:
            warnings.append(f"Column '{col}' has {miss_pct:.1f}% missing values")
    
    # Check date column if temporal
    if meta.get('temporal_unit') == 'daily':
        date_cols = [c for c in df.columns if 'date' in c.lower()]
        for date_col in date_cols:
            if date_col in df.columns:
                try:
                    dates = pd.to_datetime(df[date_col])
                    if dates.min().year < 2000:
                        warnings.append(f"Date column '{date_col}' has dates before 2000")
                    if dates.max().year > 2025:
                        issues.append(f"Date column '{date_col}' has future dates")
                except:
                    warnings.append(f"Could not parse dates in '{date_col}'")
    
    # Check region codes
    region_cols = [c for c in df.columns if 'code' in c.lower() and 'muni' not in c.lower()]
    for col in region_cols:
        if df[col].isna().any():
            n_missing = df[col].isna().sum()
            pct_missing = n_missing / len(df) * 100
            if pct_missing > 5:
                issues.append(f"Region column '{col}' has {pct_missing:.1f}% missing values")
            else:
                warnings.append(f"Region column '{col}' has {n_missing} missing values")
    
    # Check for duplicate rows (key columns)
    if meta.get('temporal_unit') == 'daily':
        key_cols = [c for c in df.columns if ('code' in c.lower() and 'muni' not in c.lower()) or c == 'date']
        if len(key_cols) >= 2:
            dupes = df.duplicated(subset=key_cols).sum()
            if dupes > 0:
                warnings.append(f"Found {dupes} duplicate rows on key columns {key_cols}")
    
    passed = len(issues) == 0
    return {'issues': issues, 'warnings': warnings, 'passed': passed}


# =============================================================================
# GENERATE DATA DICTIONARY
# =============================================================================

def generate_data_dictionary():
    """Generate comprehensive data dictionary."""
    
    print("\n" + "-" * 80)
    print("GENERATING DATA DICTIONARY")
    print("-" * 80)
    
    lines = []
    lines.append("# Phase 0 Data Dictionary")
    lines.append("")
    lines.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")
    lines.append("This document provides column-by-column documentation for all Phase 0 aggregated datasets.")
    lines.append("")
    lines.append("---")
    lines.append("")
    
    # Group by category
    categories = ['Exposure', 'Outcome', 'Confounder', 'Covariate', 'Reference']
    
    for category in categories:
        cat_datasets = {k: v for k, v in DATASETS.items() if v['category'] == category}
        if not cat_datasets:
            continue
        
        lines.append(f"## {category} Data")
        lines.append("")
        
        for name, meta in cat_datasets.items():
            # Skip legacy files in main dictionary
            if '[LEGACY]' in meta['description']:
                continue
                
            print(f"  Processing {name}...")
            
            df = load_dataset(meta['file'])
            if df is None:
                lines.append(f"### {name}")
                lines.append(f"**Status:** ❌ File not found: `{meta['file'].name}`")
                lines.append("")
                continue
            
            lines.append(f"### {name}")
            lines.append("")
            lines.append(f"**Description:** {meta['description']}")
            lines.append("")
            lines.append(f"**File:** `{meta['file'].name}`")
            lines.append("")
            lines.append(f"- **Rows:** {len(df):,}")
            lines.append(f"- **Columns:** {len(df.columns)}")
            lines.append(f"- **Spatial Unit:** {meta['spatial_unit']}")
            lines.append(f"- **Temporal Unit:** {meta['temporal_unit']}")
            lines.append(f"- **Memory:** {df.memory_usage(deep=True).sum() / 1e6:.1f} MB")
            lines.append("")
            
            # Column table
            lines.append("#### Columns")
            lines.append("")
            lines.append("| Column | Type | Non-Null | Null% | Description | Statistics |")
            lines.append("|--------|------|----------|-------|-------------|------------|")
            
            for col in df.columns:
                stats = get_column_stats(df, col)
                desc = meta['key_columns'].get(col, '')
                
                # Format statistics
                stat_str = ""
                if 'min' in stats:
                    stat_str = f"Range: [{stats['min']:.2f}, {stats['max']:.2f}], Mean: {stats['mean']:.2f}"
                elif 'min_date' in stats:
                    stat_str = f"{stats['min_date']} to {stats['max_date']}"
                elif 'sample_values' in stats:
                    stat_str = f"e.g., {', '.join(stats['sample_values'][:2])}"
                
                lines.append(f"| `{col}` | {stats['dtype']} | {stats['non_null']:,} | {stats['null_pct']:.1f}% | {desc} | {stat_str} |")
            
            lines.append("")
    
    # Write to file
    output_path = OUTPUT_DIR / 'data_dictionary.md'
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    
    print(f"  ✓ Saved: {output_path}")
    return output_path


# =============================================================================
# GENERATE SCHEMA SPECIFICATION
# =============================================================================

def generate_schema_specification():
    """Generate unified schema specification for appendix."""
    
    print("\n" + "-" * 80)
    print("GENERATING SCHEMA SPECIFICATION")
    print("-" * 80)
    
    lines = []
    lines.append("# Unified Data Schema Specification")
    lines.append("")
    lines.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")
    lines.append("This document provides a formal schema specification for all datasets used in the")
    lines.append("temperature-mortality analysis. Suitable for publication in supplementary materials.")
    lines.append("")
    lines.append("---")
    lines.append("")
    
    # Overview table
    lines.append("## Dataset Overview")
    lines.append("")
    lines.append("| Dataset | Category | Spatial Unit | Temporal Unit | Rows | Columns |")
    lines.append("|---------|----------|--------------|---------------|------|---------|")
    
    for name, meta in DATASETS.items():
        # Skip legacy files
        if '[LEGACY]' in meta['description']:
            continue
            
        df = load_dataset(meta['file'])
        if df is not None:
            lines.append(f"| {name} | {meta['category']} | {meta['spatial_unit']} | {meta['temporal_unit']} | {len(df):,} | {len(df.columns)} |")
        else:
            lines.append(f"| {name} | {meta['category']} | {meta['spatial_unit']} | {meta['temporal_unit']} | — | — |")
    
    lines.append("")
    
    # Spatial units
    lines.append("## Spatial Units")
    lines.append("")
    lines.append("| Level | IBGE Name | Count | Description |")
    lines.append("|-------|-----------|-------|-------------|")
    lines.append("| Intermediate | Região Geográfica Intermediária | 133 | Larger regions for stable estimates |")
    lines.append("| Immediate | Região Geográfica Imediata | 510 | Finer spatial resolution |")
    lines.append("| National | Brasil | 1 | Country-level data |")
    lines.append("")
    
    # Variable definitions
    lines.append("## Variable Definitions")
    lines.append("")
    
    lines.append("### Identifier Variables")
    lines.append("")
    lines.append("| Variable | Description | Format |")
    lines.append("|----------|-------------|--------|")
    lines.append("| `intermediate_code` | Intermediate region identifier | Integer (IBGE code) |")
    lines.append("| `immediate_code` | Immediate region identifier | Integer (IBGE code) |")
    lines.append("| `region_code` | [Legacy] Intermediate region identifier | Integer (IBGE code) |")
    lines.append("| `code_muni` | Municipality code | 7-digit integer |")
    lines.append("| `date` | Date of observation | YYYY-MM-DD |")
    lines.append("")
    
    lines.append("### Exposure Variables (Temperature)")
    lines.append("")
    lines.append("| Variable | Unit | Source | Definition |")
    lines.append("|----------|------|--------|------------|")
    lines.append("| `temp_mean` | °C | ERA5 | Daily mean 2m temperature, population-weighted |")
    lines.append("| `temp_min` | °C | ERA5 | Daily minimum 2m temperature |")
    lines.append("| `temp_max` | °C | ERA5 | Daily maximum 2m temperature |")
    lines.append("| `dewpoint_mean` | °C | ERA5 | Daily mean 2m dewpoint for humidity calculation |")
    lines.append("")
    
    lines.append("### Exposure Variables (Pollution)")
    lines.append("")
    lines.append("| Variable | Unit | Source | Definition |")
    lines.append("|----------|------|--------|------------|")
    lines.append("| `pm25_mean` | μg/m³ | CAMS | Daily mean PM2.5 concentration |")
    lines.append("| `o3_mean` | μg/m³ | CAMS | Daily mean O3 (ozone) concentration |")
    lines.append("")
    
    lines.append("### Outcome Variables (Mortality)")
    lines.append("")
    lines.append("| Variable | Definition | ICD-10 Codes |")
    lines.append("|----------|------------|--------------|")
    lines.append("| `deaths_all` | All-cause mortality, all ages | All |")
    lines.append("| `deaths_elderly` | All-cause mortality, age ≥60 | All |")
    lines.append("| `deaths_respiratory` / `deaths_elderly_resp` | Respiratory mortality | J00-J99 |")
    lines.append("| `deaths_cardiovascular` / `deaths_elderly_cvd` | Cardiovascular mortality | I00-I99 |")
    lines.append("| `deaths_heat_direct` / `deaths_elderly_heat` | Direct heat-related mortality | T67, X30 |")
    lines.append("")
    
    lines.append("### Confounder Variables")
    lines.append("")
    lines.append("| Variable | Unit | Source | Definition |")
    lines.append("|----------|------|--------|------------|")
    lines.append("| `srag_cases` | Count | SIVEP-Gripe | Severe acute respiratory infection cases |")
    lines.append("| `srag_influenza` | Count | SIVEP-Gripe | Laboratory-confirmed influenza |")
    lines.append("| `srag_covid` | Count | SIVEP-Gripe | Laboratory-confirmed COVID-19 |")
    lines.append("| `is_holiday` | Binary | holidays-br | National holiday indicator |")
    lines.append("| `is_holiday_week` | Binary | Derived | Week contains a national holiday |")
    lines.append("")
    
    lines.append("### Covariate Variables")
    lines.append("")
    lines.append("| Variable | Unit | Source | Definition |")
    lines.append("|----------|------|--------|------------|")
    lines.append("| `pop_total` | Count | IBGE Census | Total population |")
    lines.append("| `pop_elderly` | Count | IBGE Census | Population aged ≥60 |")
    lines.append("| `pct_elderly` | % | Derived | Percentage of population ≥60 |")
    lines.append("| `gdp_per_capita` | R$ | IBGE | Regional GDP per capita |")
    lines.append("| `urbanization_rate` | % | IBGE | Percentage urban population |")
    lines.append("| `ac_ownership` | % | PNAD | Air conditioning ownership rate |")
    lines.append("")
    
    lines.append("### Reference Variables (Life Tables)")
    lines.append("")
    lines.append("| Variable | Unit | Source | Definition |")
    lines.append("|----------|------|--------|------------|")
    lines.append("| `life_expectancy` | Years | IBGE | Remaining life expectancy at given age |")
    lines.append("| `yll` | Years | Derived | Years of life lost if death at given age |")
    lines.append("")
    
    # Data sources
    lines.append("## Data Sources")
    lines.append("")
    lines.append("| Source | Dataset | Time Period | URL |")
    lines.append("|--------|---------|-------------|-----|")
    lines.append("| ECMWF | ERA5 Reanalysis | 2010-2024 | https://cds.climate.copernicus.eu |")
    lines.append("| ECMWF | CAMS Reanalysis | 2010-2024 | https://ads.atmosphere.copernicus.eu |")
    lines.append("| DATASUS | SIM (Mortality) | 2010-2024 | https://datasus.saude.gov.br |")
    lines.append("| DATASUS | SIVEP-Gripe | 2010-2024 | https://opendatasus.saude.gov.br |")
    lines.append("| IBGE | Census/SIDRA | 2010-2022 | https://sidra.ibge.gov.br |")
    lines.append("| IBGE | Life Tables | 2010-2024 | https://ibge.gov.br |")
    lines.append("| IBGE | Geographic Regions | 2017 | https://ibge.gov.br/geociencias |")
    lines.append("")
    
    # Write to file
    output_path = OUTPUT_DIR / 'schema_specification.md'
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    
    print(f"  ✓ Saved: {output_path}")
    return output_path


# =============================================================================
# GENERATE VALIDATION REPORT
# =============================================================================

def generate_validation_report():
    """Generate validation report checking completeness and alignment."""
    
    print("\n" + "-" * 80)
    print("GENERATING VALIDATION REPORT")
    print("-" * 80)
    
    lines = []
    lines.append("# Data Validation Report")
    lines.append("")
    lines.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")
    lines.append("This report validates completeness, missingness, and alignment across all datasets.")
    lines.append("")
    lines.append("---")
    lines.append("")
    
    # Summary
    total_datasets = len([d for d in DATASETS.values() if '[LEGACY]' not in d['description']])
    loaded_datasets = 0
    passed_datasets = 0
    all_validations = {}
    
    lines.append("## Validation Summary")
    lines.append("")
    
    for name, meta in DATASETS.items():
        # Skip legacy files in validation
        if '[LEGACY]' in meta['description']:
            continue
            
        print(f"  Validating {name}...")
        df = load_dataset(meta['file'])
        
        if df is None:
            all_validations[name] = {'status': 'NOT_FOUND', 'issues': ['File not found'], 'warnings': []}
            continue
        
        loaded_datasets += 1
        validation = validate_dataset(df, name, meta)
        validation['status'] = 'PASS' if validation['passed'] else 'FAIL'
        all_validations[name] = validation
        
        if validation['passed']:
            passed_datasets += 1
    
    lines.append(f"- **Datasets defined:** {total_datasets}")
    lines.append(f"- **Datasets found:** {loaded_datasets}")
    lines.append(f"- **Datasets passed validation:** {passed_datasets}")
    lines.append(f"- **Datasets with issues:** {loaded_datasets - passed_datasets}")
    lines.append("")
    
    # Status table
    lines.append("## Dataset Status")
    lines.append("")
    lines.append("| Dataset | Status | Issues | Warnings |")
    lines.append("|---------|--------|--------|----------|")
    
    for name, val in all_validations.items():
        status_icon = "✅" if val['status'] == 'PASS' else ("❌" if val['status'] == 'FAIL' else "⚠️")
        n_issues = len(val.get('issues', []))
        n_warnings = len(val.get('warnings', []))
        lines.append(f"| {name} | {status_icon} {val['status']} | {n_issues} | {n_warnings} |")
    
    lines.append("")
    
    # Detailed issues
    lines.append("## Detailed Issues")
    lines.append("")
    
    has_issues = False
    for name, val in all_validations.items():
        if val.get('issues') or val.get('warnings'):
            has_issues = True
            lines.append(f"### {name}")
            lines.append("")
            
            if val.get('issues'):
                lines.append("**Issues (must fix):**")
                for issue in val['issues']:
                    lines.append(f"- ❌ {issue}")
                lines.append("")
            
            if val.get('warnings'):
                lines.append("**Warnings (review):**")
                for warning in val['warnings']:
                    lines.append(f"- ⚠️ {warning}")
                lines.append("")
    
    if not has_issues:
        lines.append("*No issues or warnings found.*")
        lines.append("")
    
    # Cross-dataset alignment
    lines.append("## Cross-Dataset Alignment")
    lines.append("")
    
    # Check date ranges align
    lines.append("### Temporal Coverage")
    lines.append("")
    lines.append("| Dataset | Start Date | End Date | Days |")
    lines.append("|---------|------------|----------|------|")
    
    for name, meta in DATASETS.items():
        if meta.get('temporal_unit') != 'daily' or '[LEGACY]' in meta['description']:
            continue
        df = load_dataset(meta['file'])
        if df is None:
            continue
        
        date_col = None
        for col in df.columns:
            if 'date' in col.lower():
                date_col = col
                break
        
        if date_col:
            try:
                dates = pd.to_datetime(df[date_col])
                lines.append(f"| {name} | {dates.min().date()} | {dates.max().date()} | {(dates.max() - dates.min()).days:,} |")
            except:
                lines.append(f"| {name} | — | — | — |")
    
    lines.append("")
    
    # Check region coverage
    lines.append("### Spatial Coverage")
    lines.append("")
    lines.append("| Dataset | Spatial Unit | Unique Regions | Expected | Status |")
    lines.append("|---------|--------------|----------------|----------|--------|")
    
    for name, meta in DATASETS.items():
        if '[LEGACY]' in meta['description']:
            continue
            
        df = load_dataset(meta['file'])
        if df is None:
            continue
        
        spatial = meta['spatial_unit']
        expected = {'intermediate': 133, 'immediate': 510, 'national': 1, 'municipality': 5570}.get(spatial, None)
        
        # Find region column
        region_col = None
        for col in df.columns:
            if 'code' in col.lower() and ('intermediate' in col.lower() or 'immediate' in col.lower() or col == 'region_code'):
                region_col = col
                break
        
        if region_col and expected:
            unique = df[region_col].nunique()
            match = "✅" if unique >= expected * 0.95 else "⚠️"
            lines.append(f"| {name} | {spatial} | {unique} | {expected} | {match} |")
        elif expected:
            lines.append(f"| {name} | {spatial} | — | {expected} | ⚠️ |")
    
    lines.append("")
    
    # Missingness summary
    lines.append("### Missingness Summary")
    lines.append("")
    lines.append("*Only showing columns with >5% missing values*")
    lines.append("")
    lines.append("| Dataset | Column | Missing % | Assessment |")
    lines.append("|---------|--------|-----------|------------|")
    
    has_missing = False
    for name, meta in DATASETS.items():
        if '[LEGACY]' in meta['description']:
            continue
            
        df = load_dataset(meta['file'])
        if df is None:
            continue
        
        for col in df.columns:
            miss_pct = df[col].isna().mean() * 100
            if miss_pct > 5:  # Only report columns with >5% missing
                has_missing = True
                assessment = "❌ High" if miss_pct > 20 else "⚠️ Moderate"
                lines.append(f"| {name} | {col} | {miss_pct:.1f}% | {assessment} |")
    
    if not has_missing:
        lines.append("| *None* | — | — | All columns <5% missing |")
    
    lines.append("")
    
    # Write to file
    output_path = OUTPUT_DIR / 'validation_report.md'
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    
    print(f"  ✓ Saved: {output_path}")
    return output_path, all_validations


# =============================================================================
# GENERATE METADATA JSON
# =============================================================================

def generate_metadata():
    """Generate structured metadata JSON."""
    
    print("\n" + "-" * 80)
    print("GENERATING METADATA")
    print("-" * 80)
    
    metadata = {
        'generated': datetime.now().isoformat(),
        'generator': '00a_document_data.py',
        'datasets': {}
    }
    
    for name, meta in DATASETS.items():
        # Include legacy files in metadata for completeness
        print(f"  Processing {name}...")
        
        df = load_dataset(meta['file'])
        
        ds_meta = {
            'file': meta['file'].name,
            'file_path': str(meta['file']),
            'file_exists': df is not None,
            'category': meta['category'],
            'description': meta['description'],
            'spatial_unit': meta['spatial_unit'],
            'temporal_unit': meta['temporal_unit'],
            'is_legacy': '[LEGACY]' in meta['description'],
        }
        
        if df is not None:
            ds_meta['rows'] = len(df)
            ds_meta['columns'] = len(df.columns)
            ds_meta['column_names'] = list(df.columns)
            ds_meta['memory_mb'] = round(df.memory_usage(deep=True).sum() / 1e6, 2)
            
            # Column types
            ds_meta['column_types'] = {col: str(df[col].dtype) for col in df.columns}
            
            # Date range
            for col in df.columns:
                if 'date' in col.lower():
                    try:
                        dates = pd.to_datetime(df[col])
                        ds_meta['date_range'] = {
                            'start': str(dates.min().date()),
                            'end': str(dates.max().date()),
                            'days': int((dates.max() - dates.min()).days)
                        }
                    except:
                        pass
                    break
            
            # Region coverage
            for col in df.columns:
                if 'code' in col.lower() and 'muni' not in col.lower():
                    expected = {
                        'intermediate': 133, 
                        'immediate': 510, 
                        'national': 1
                    }.get(meta['spatial_unit'], None)
                    
                    ds_meta['region_coverage'] = {
                        'column': col,
                        'unique_regions': int(df[col].nunique()),
                        'expected': expected
                    }
                    break
            
            # Summary statistics for key numeric columns
            ds_meta['summary_stats'] = {}
            for col in df.columns:
                if pd.api.types.is_numeric_dtype(df[col]):
                    valid = df[col].dropna()
                    if len(valid) > 0:
                        ds_meta['summary_stats'][col] = {
                            'min': float(valid.min()),
                            'max': float(valid.max()),
                            'mean': round(float(valid.mean()), 4),
                            'std': round(float(valid.std()), 4),
                            'missing_pct': round(df[col].isna().mean() * 100, 2)
                        }
        
        metadata['datasets'][name] = ds_meta
    
    # Summary (exclude legacy files)
    active_datasets = {k: v for k, v in metadata['datasets'].items() if not v.get('is_legacy', False)}
    
    metadata['summary'] = {
        'total_datasets': len(active_datasets),
        'datasets_found': sum(1 for d in active_datasets.values() if d['file_exists']),
        'total_rows': sum(d.get('rows', 0) for d in active_datasets.values()),
        'categories': list(set(d['category'] for d in active_datasets.values())),
        'legacy_files_count': sum(1 for d in metadata['datasets'].values() if d.get('is_legacy', False)),
    }
    
    # Write to file
    output_path = OUTPUT_DIR / 'metadata.json'
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)
    
    print(f"  ✓ Saved: {output_path}")
    return output_path, metadata


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == '__main__':
    
    print("\n" + "=" * 80)
    print("STARTING COMPREHENSIVE DATA DOCUMENTATION")
    print("=" * 80)
    
    # Generate all outputs
    dict_path = generate_data_dictionary()
    schema_path = generate_schema_specification()
    validation_path, validations = generate_validation_report()
    metadata_path, metadata = generate_metadata()
    
    # Final summary
    print("\n" + "=" * 80)
    print("DOCUMENTATION COMPLETE")
    print("=" * 80)
    
    print("\nGenerated files:")
    print(f"  1. {dict_path.name} - Column-by-column documentation")
    print(f"  2. {schema_path.name} - Unified schema for appendix")
    print(f"  3. {validation_path.name} - Completeness and alignment checks")
    print(f"  4. {metadata_path.name} - Structured metadata")
    
    print(f"\nDataset summary:")
    print(f"  - Datasets defined: {metadata['summary']['total_datasets']}")
    print(f"  - Datasets found: {metadata['summary']['datasets_found']}")
    print(f"  - Total rows: {metadata['summary']['total_rows']:,}")
    print(f"  - Legacy files: {metadata['summary']['legacy_files_count']}")
    
    # Validation summary
    passed = sum(1 for v in validations.values() if v.get('status') == 'PASS')
    failed = sum(1 for v in validations.values() if v.get('status') == 'FAIL')
    missing = sum(1 for v in validations.values() if v.get('status') == 'NOT_FOUND')
    
    print(f"\nValidation summary:")
    print(f"  - ✅ Passed: {passed}")
    print(f"  - ❌ Failed: {failed}")
    print(f"  - ⚠️ Not found: {missing}")
    
    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
