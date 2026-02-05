"""
00m: DOWNLOAD INFLUENZA/SRAG DATA FROM OPENDATASUS
===================================================
Download Severe Acute Respiratory Syndrome (SRAG) data from Brazil's OpenDataSUS.
This data includes confirmed influenza cases used to control for flu confounding
in temperature-mortality analysis.

Data Source: SIVEP-Gripe (Sistema de Informação de Vigilância Epidemiológica da Gripe)
Historical: Sinan Web Influenza (2009-2018)

Datasets:
- SRAG 2009-2012: https://opendatasus.saude.gov.br/dataset/srag-2009-2012 (CSV)
- SRAG 2013-2018: https://opendatasus.saude.gov.br/dataset/srag-2013-2018 (CSV)
- SRAG 2019-2024: https://opendatasus.saude.gov.br/dataset/srag-2021-a-2024 (Parquet)

Why control for influenza:
- Cold-mortality effects may be confounded by flu season (peaks in winter)
- Flu causes ~20,000-30,000 deaths/year in Brazil, concentrated in cold months
- Without control, cold effects may be overestimated

Variables of interest:
- CLASSI_FIN: Final classification (1=SRAG by influenza, 2=other virus, etc.)
- DT_NOTIFIC: Notification date
- DT_SIN_PRI: Symptom onset date  
- SG_UF_NOT: State of notification
- EVOLUCAO: Outcome (1=cure, 2=death, 3=death other cause)
- NU_IDADE_N: Age in years
- CS_SEXO: Sex (M/F)

Output: 
- results/influenza_weekly_by_state.parquet (aggregated for modeling)
- results/influenza_raw_elderly.parquet (filtered raw data for elderly 60+)

Documentation: See data dictionary at OpenDataSUS
"""

import os
import sys
from datetime import datetime
from pathlib import Path
import warnings

print("=" * 70)
print("00m: DOWNLOAD INFLUENZA/SRAG DATA (2010-2024)")
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
    import requests
    print("✓ requests installed")
except ImportError:
    print("✗ requests not installed - install with: pip install requests")
    sys.exit(1)

try:
    import pyarrow.parquet as pq
    print("✓ pyarrow installed")
except ImportError:
    print("✗ pyarrow not installed - install with: pip install pyarrow")
    sys.exit(1)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Output directory
BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'results'
TEMP_DIR = BASE_DIR / 'temp_srag'
OUTPUT_DIR.mkdir(exist_ok=True)
TEMP_DIR.mkdir(exist_ok=True)

# =============================================================================
# DATA SOURCE URLS
# =============================================================================

# Historical data (2009-2018) uses CSV format from Sinan Web Influenza
# The CSV files are accessed via CKAN API - we need resource IDs
# Format: https://opendatasus.saude.gov.br/dataset/{dataset}/resource/{resource_id}

# CSV files are served via CloudFront CDN (not direct S3 access)
SRAG_CSV_URLS = {
    # 2009-2012 dataset (from CloudFront)
    2009: "https://d26692udehoye.cloudfront.net/SRAG/2009-2012/INFLUD09.csv",
    2010: "https://d26692udehoye.cloudfront.net/SRAG/2009-2012/INFLUD10.csv",
    2011: "https://d26692udehoye.cloudfront.net/SRAG/2009-2012/INFLUD11.csv",
    2012: "https://d26692udehoye.cloudfront.net/SRAG/2009-2012/INFLUD12.csv",
    # 2013-2018 dataset (from CloudFront)
    2013: "https://d26692udehoye.cloudfront.net/SRAG/2013-2018/INFLUD13.csv",
    2014: "https://d26692udehoye.cloudfront.net/SRAG/2013-2018/INFLUD14.csv",
    2015: "https://d26692udehoye.cloudfront.net/SRAG/2013-2018/INFLUD15.csv",
    2016: "https://d26692udehoye.cloudfront.net/SRAG/2013-2018/INFLUD16.csv",
    2017: "https://d26692udehoye.cloudfront.net/SRAG/2013-2018/INFLUD17.csv",
    2018: "https://d26692udehoye.cloudfront.net/SRAG/2013-2018/INFLUD18.csv",
}

# Modern data (2019+) uses Parquet format from SIVEP-Gripe
SRAG_PARQUET_URLS = {
    2019: "https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2019/INFLUD19-26-06-2025.parquet",
    2020: "https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2020/INFLUD20-26-06-2025.parquet",
    2021: "https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2021/INFLUD21-26-06-2025.parquet",
    2022: "https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2022/INFLUD22-26-06-2025.parquet",
    2023: "https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2023/INFLUD23-26-06-2025.parquet",
    2024: "https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2024/INFLUD24-26-06-2025.parquet",
}

# Years to download (matching our mortality data period: 2010-2024)
YEARS = list(range(2010, 2025))

# Classification codes for influenza
# CLASSI_FIN values:
# 1 = SRAG por influenza
# 2 = SRAG por outro vírus respiratório  
# 3 = SRAG por outro agente etiológico
# 4 = SRAG não especificado
# 5 = SRAG por COVID-19
INFLUENZA_CODE = 1

# State codes mapping
STATE_CODES = {
    'AC': 12, 'AL': 27, 'AM': 13, 'AP': 16, 'BA': 29, 'CE': 23, 'DF': 53,
    'ES': 32, 'GO': 52, 'MA': 21, 'MG': 31, 'MS': 50, 'MT': 51, 'PA': 15,
    'PB': 25, 'PE': 26, 'PI': 22, 'PR': 41, 'RJ': 33, 'RN': 24, 'RO': 11,
    'RR': 14, 'RS': 43, 'SC': 42, 'SE': 28, 'SP': 35, 'TO': 17
}

CODE_TO_STATE = {v: k for k, v in STATE_CODES.items()}

# =============================================================================
# DOWNLOAD FUNCTIONS
# =============================================================================

def download_srag_year(year: int) -> pd.DataFrame:
    """
    Download SRAG data for one year from OpenDataSUS.
    
    - 2010-2018: CSV format (Sinan Web Influenza)
    - 2019-2024: Parquet format (SIVEP-Gripe)
    
    Returns DataFrame with relevant columns only.
    """
    
    # Determine format and URL
    if year in SRAG_CSV_URLS:
        url = SRAG_CSV_URLS[year]
        file_format = 'csv'
        local_file = TEMP_DIR / f"srag_{year}.csv"
    elif year in SRAG_PARQUET_URLS:
        url = SRAG_PARQUET_URLS[year]
        file_format = 'parquet'
        local_file = TEMP_DIR / f"srag_{year}.parquet"
    else:
        print(f"  ⚠️ {year}: No data URL configured")
        return None
    
    # Check if already downloaded
    if local_file.exists():
        print(f"  ✓ {year}: Using cached file ({file_format})")
        try:
            if file_format == 'csv':
                # CSV files may have encoding issues - try multiple encodings
                for encoding in ['latin1', 'utf-8', 'cp1252']:
                    try:
                        return pd.read_csv(local_file, encoding=encoding, 
                                          low_memory=False, sep=';')
                    except UnicodeDecodeError:
                        continue
                # If semicolon doesn't work, try comma
                return pd.read_csv(local_file, encoding='latin1', low_memory=False)
            else:
                return pd.read_parquet(local_file)
        except Exception as e:
            print(f"    ⚠️ Cached file error: {e}, re-downloading...")
            local_file.unlink()
    
    print(f"  ↓ {year}: Downloading ({file_format})...")
    
    try:
        # Download with streaming to handle large files
        response = requests.get(url, stream=True, timeout=600)
        response.raise_for_status()
        
        # Save to temp file
        total_size = int(response.headers.get('content-length', 0))
        downloaded = 0
        
        with open(local_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded += len(chunk)
                if total_size > 0:
                    pct = downloaded / total_size * 100
                    print(f"\r    Progress: {pct:.1f}% ({downloaded // 1024 // 1024} MB)", end='')
        
        print(f"\n  ✓ {year}: Downloaded successfully ({downloaded // 1024 // 1024} MB)")
        
        # Read and return
        if file_format == 'csv':
            # Try multiple encodings and separators for CSV
            for encoding in ['latin1', 'utf-8', 'cp1252']:
                for sep in [';', ',']:
                    try:
                        df = pd.read_csv(local_file, encoding=encoding, 
                                        low_memory=False, sep=sep)
                        # Verify we got columns (not just one column with delimiters)
                        if len(df.columns) > 5:
                            return df
                    except:
                        continue
            # Last resort
            df = pd.read_csv(local_file, encoding='latin1', low_memory=False)
            return df
        else:
            df = pd.read_parquet(local_file)
            return df
        
    except Exception as e:
        print(f"\n  ✗ {year}: Error - {str(e)[:100]}")
        if local_file.exists():
            local_file.unlink()
        return None


def process_srag_data(df: pd.DataFrame, year: int) -> pd.DataFrame:
    """
    Process raw SRAG data to extract relevant information.
    
    Handles both:
    - 2010-2018: Sinan Web format (column names may differ slightly)
    - 2019-2024: SIVEP-Gripe format
    
    Focus on:
    - Confirmed influenza cases (CLASSI_FIN == 1)
    - Elderly population (60+)
    - State and date for aggregation
    """
    
    if df is None or len(df) == 0:
        return None
    
    print(f"  Processing {year}: {len(df):,} total SRAG records")
    
    # Normalize column names (uppercase, strip whitespace)
    df.columns = [str(c).upper().strip() for c in df.columns]
    
    # Column name variations across years
    # Date columns
    date_col_options = ['DT_SIN_PRI', 'DT_NOTIFIC', 'DT_NOTIFICACAO', 'DT_DIGITA']
    # State columns  
    state_col_options = ['SG_UF_NOT', 'SG_UF', 'UF_NOT', 'CO_UF_NOT', 'UF']
    # Classification columns
    class_col_options = ['CLASSI_FIN', 'CLASSIFIN', 'CLASSI_GES', 'CLASSIFICACAO']
    # Age columns
    age_col_options = ['NU_IDADE_N', 'IDADE', 'NU_IDADE']
    age_type_options = ['TP_IDADE', 'CS_IDADE', 'COD_IDADE']
    # Outcome columns
    outcome_col_options = ['EVOLUCAO', 'EVOLUACAO', 'CS_EVOLUCA']
    # Sex columns
    sex_col_options = ['CS_SEXO', 'SEXO']
    
    def find_column(df, options):
        """Find first matching column from options list."""
        for opt in options:
            if opt in df.columns:
                return opt
        return None
    
    # Find available columns
    date_col = find_column(df, date_col_options)
    state_col = find_column(df, state_col_options)
    class_col = find_column(df, class_col_options)
    age_col = find_column(df, age_col_options)
    age_type_col = find_column(df, age_type_options)
    outcome_col = find_column(df, outcome_col_options)
    sex_col = find_column(df, sex_col_options)
    
    # Report what we found
    print(f"    Columns found: date={date_col}, state={state_col}, class={class_col}, age={age_col}, age_type={age_type_col}")
    
    if date_col is None:
        print(f"    ✗ No date column found. Available: {list(df.columns)[:20]}")
        return None
    
    # Parse date
    df['date'] = pd.to_datetime(df[date_col], errors='coerce', dayfirst=True)
    
    # Also try secondary date column if primary has too many NaTs
    if df['date'].isna().sum() > len(df) * 0.5:
        for backup_date in date_col_options:
            if backup_date in df.columns and backup_date != date_col:
                backup = pd.to_datetime(df[backup_date], errors='coerce', dayfirst=True)
                df['date'] = df['date'].fillna(backup)
    
    # Filter to valid dates
    df = df[df['date'].notna()].copy()
    
    if len(df) == 0:
        print(f"    ✗ No valid dates after parsing")
        return None
    
    # Extract state
    if state_col:
        df['state'] = df[state_col].astype(str).str.upper().str.strip()
        # Convert numeric codes to state abbreviations if needed
        if df['state'].str.match(r'^\d+$').any():
            df['state'] = df['state'].apply(
                lambda x: CODE_TO_STATE.get(int(float(x)), x) if x.isdigit() else x
            )
    else:
        df['state'] = 'UNKNOWN'
    
    # Calculate age in years
    # Two different formats exist:
    # 1) CSV 2009-2018: NU_IDADE_N contains encoded value like 4042 (first digit=type, rest=value)
    #    Type: 1=hours, 2=days, 3=months, 4=years
    #    Example: 4042 = 42 years, 3011 = 11 months
    # 2) Parquet 2019+: NU_IDADE_N and TP_IDADE are separate columns
    #    BUT the type codes are DIFFERENT: 1=hours, 2=days, 3=YEARS (not months!)
    #    This is confirmed by median age ~61 in type 3 records
    
    if age_col:
        age_raw = pd.to_numeric(df[age_col], errors='coerce')
        
        if age_type_col:
            # Format 2: Separate columns (Parquet 2019+)
            # IMPORTANT: In SIVEP-Gripe (2019+), type codes are:
            # 1=hours, 2=days, 3=YEARS (different from CSV!)
            age_type = pd.to_numeric(df[age_type_col], errors='coerce')
            age_value = age_raw.astype(float)
            
            # Debug: show distribution
            type_counts = age_type.value_counts().head(5)
            print(f"    Age type distribution: {dict(type_counts)}")
            
            # Convert based on type (SIVEP-Gripe format)
            # Type 3 = years, Type 2 = days, Type 1 = hours
            df['age_years'] = np.nan
            mask_years = age_type == 3  # In SIVEP-Gripe, 3 = years!
            mask_days = age_type == 2
            mask_hours = age_type == 1
            
            df.loc[mask_years, 'age_years'] = age_value[mask_years]
            df.loc[mask_days, 'age_years'] = age_value[mask_days] / 365.0
            df.loc[mask_hours, 'age_years'] = 0.0
            
        else:
            # Format 1: Encoded in single column (CSV 2009-2018)
            # First digit is type, remaining digits are value
            # Type: 1=hours, 2=days, 3=months, 4=years
            # Example: 4042 -> type=4, value=42 -> 42 years
            print(f"    Using encoded age format (first digit=type)")
            
            # Extract type (first digit) and value (remaining digits)
            age_type = (age_raw // 1000).astype(float)  # First digit
            age_value = (age_raw % 1000).astype(float)  # Remaining digits
            
            # Debug
            type_counts = age_type.value_counts().head(5)
            print(f"    Encoded age type distribution: {dict(type_counts)}")
            
            # Convert based on type (Sinan format)
            # Type 4 = years, Type 3 = months, Type 2 = days, Type 1 = hours
            df['age_years'] = np.nan
            mask_years = age_type == 4
            mask_months = age_type == 3
            mask_days = age_type == 2
            mask_hours = age_type == 1
            
            df.loc[mask_years, 'age_years'] = age_value[mask_years]
            df.loc[mask_months, 'age_years'] = age_value[mask_months] / 12.0
            df.loc[mask_days, 'age_years'] = age_value[mask_days] / 365.0
            df.loc[mask_hours, 'age_years'] = 0.0
    else:
        df['age_years'] = np.nan
    
    # Ensure age_years is float for comparison
    df['age_years'] = pd.to_numeric(df['age_years'], errors='coerce').astype(float)
    
    # Classification flags
    if class_col:
        class_val = pd.to_numeric(df[class_col], errors='coerce')
        df['is_influenza'] = class_val == INFLUENZA_CODE
        df['is_covid'] = class_val == 5
        df['is_other_virus'] = class_val == 2
    else:
        # For older years without classification, treat all as potential influenza
        df['is_influenza'] = True  # Conservative: count all SRAG
        df['is_covid'] = False
        df['is_other_virus'] = False
        print(f"    ⚠️ No classification column - counting all SRAG as potential influenza")
    
    # Outcome flag
    if outcome_col:
        outcome_val = pd.to_numeric(df[outcome_col], errors='coerce')
        df['is_death'] = outcome_val.isin([2, 3])  # 2=death by SRAG, 3=death other
    else:
        df['is_death'] = False
    
    # Elderly flag (60+) - explicit comparison with float
    df['is_elderly'] = (df['age_years'].notna()) & (df['age_years'] >= 60.0)
    
    # Add year for reference
    df['year'] = year
    
    # Summary stats
    n_influenza = df['is_influenza'].sum()
    n_elderly = df['is_elderly'].sum()
    n_elderly_flu = ((df['is_influenza']) & (df['is_elderly'])).sum()
    
    print(f"    → Total records (with valid date): {len(df):,}")
    print(f"    → Influenza/SRAG cases: {n_influenza:,}")
    print(f"    → Elderly (60+): {n_elderly:,}")
    print(f"    → Elderly influenza: {n_elderly_flu:,}")
    
    return df


def aggregate_to_state_week(all_data: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate SRAG data to state-week level for modeling.
    
    Creates weekly counts of:
    - Total SRAG cases
    - Confirmed influenza cases
    - COVID-19 cases
    - Deaths
    """
    
    print("\nAggregating to state-week level...")
    
    # Create week column (epidemiological week)
    all_data['epi_week'] = all_data['date'].dt.isocalendar().week
    all_data['epi_year'] = all_data['date'].dt.isocalendar().year
    
    # Also create a date-based week start for easier joining
    all_data['week_start'] = all_data['date'] - pd.to_timedelta(all_data['date'].dt.dayofweek, unit='D')
    
    # Aggregate by state and week
    weekly = all_data.groupby(['state', 'week_start', 'epi_year', 'epi_week']).agg(
        srag_total=('date', 'count'),
        srag_influenza=('is_influenza', 'sum'),
        srag_covid=('is_covid', 'sum'),
        srag_other_virus=('is_other_virus', 'sum'),
        srag_deaths=('is_death', 'sum'),
        srag_elderly=('is_elderly', 'sum'),
        srag_elderly_influenza=pd.NamedAgg(
            column='is_influenza',
            aggfunc=lambda x: ((all_data.loc[x.index, 'is_influenza']) & 
                              (all_data.loc[x.index, 'is_elderly'])).sum()
        ),
    ).reset_index()
    
    print(f"  Created {len(weekly):,} state-week observations")
    
    return weekly


def aggregate_to_state_day(all_data: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate SRAG data to state-day level for direct joining with mortality data.
    """
    
    print("\nAggregating to state-day level...")
    
    # Aggregate by state and date
    daily = all_data.groupby(['state', 'date']).agg(
        srag_total=('date', 'count'),
        srag_influenza=('is_influenza', 'sum'),
        srag_covid=('is_covid', 'sum'),
        srag_deaths=('is_death', 'sum'),
        srag_elderly=('is_elderly', 'sum'),
    ).reset_index()
    
    # Rename date for joining
    daily = daily.rename(columns={'date': 'date'})
    
    print(f"  Created {len(daily):,} state-day observations")
    
    return daily


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Download and process all SRAG data."""
    
    print("\n" + "-" * 70)
    print("STEP 1: Download SRAG data from OpenDataSUS")
    print("-" * 70)
    
    all_data = []
    
    for year in YEARS:
        print(f"\n{year}:")
        
        # Download
        df_raw = download_srag_year(year)
        
        if df_raw is not None:
            # Process
            df_processed = process_srag_data(df_raw, year)
            
            if df_processed is not None:
                all_data.append(df_processed)
    
    if not all_data:
        print("\n✗ No data downloaded successfully!")
        return
    
    # Combine all years
    print("\n" + "-" * 70)
    print("STEP 2: Combine and aggregate data")
    print("-" * 70)
    
    combined = pd.concat(all_data, ignore_index=True)
    print(f"\nTotal records: {len(combined):,}")
    
    # Summary by year
    print("\nRecords by year:")
    print(combined.groupby('year').size())
    
    print("\nInfluenza cases by year:")
    print(combined[combined['is_influenza']].groupby('year').size())
    
    # Create aggregated outputs
    print("\n" + "-" * 70)
    print("STEP 3: Create output files")
    print("-" * 70)
    
    # State-week aggregation
    weekly = aggregate_to_state_week(combined)
    weekly_file = OUTPUT_DIR / 'influenza_weekly_by_state.parquet'
    weekly.to_parquet(weekly_file, index=False)
    print(f"\n✓ Saved: {weekly_file}")
    
    # State-day aggregation (for direct joining with mortality)
    daily = aggregate_to_state_day(combined)
    daily_file = OUTPUT_DIR / 'influenza_daily_by_state.parquet'
    daily.to_parquet(daily_file, index=False)
    print(f"✓ Saved: {daily_file}")
    
    # Also save filtered elderly-only raw data
    elderly_data = combined[combined['is_elderly']].copy()
    elderly_file = OUTPUT_DIR / 'influenza_raw_elderly.parquet'
    elderly_data.to_parquet(elderly_file, index=False)
    print(f"✓ Saved: {elderly_file} ({len(elderly_data):,} elderly records)")
    
    # Summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    print(f"\nTotal SRAG records: {len(combined):,}")
    print(f"  - Influenza: {combined['is_influenza'].sum():,}")
    print(f"  - COVID-19: {combined['is_covid'].sum():,}")
    print(f"  - Other virus: {combined['is_other_virus'].sum():,}")
    print(f"  - Elderly (60+): {combined['is_elderly'].sum():,}")
    print(f"  - Deaths: {combined['is_death'].sum():,}")
    
    print(f"\nDate range: {combined['date'].min()} to {combined['date'].max()}")
    print(f"States covered: {combined['state'].nunique()}")
    
    print(f"\n✓ Done! {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    return combined


if __name__ == "__main__":
    main()
