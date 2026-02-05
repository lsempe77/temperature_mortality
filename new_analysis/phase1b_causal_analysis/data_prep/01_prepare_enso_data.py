#!/usr/bin/env python3
"""
ENSO Data Preparation Script
============================

Prepares ENSO (El Niño Southern Oscillation) indices for causal analysis.
Loads multiple ENSO indices, harmonizes dates, and creates phase classifications.

ENSO Indices Available:
- ONI (Oceanic Niño Index): 3-month running mean of SST anomalies in Niño 3.4 region
- MEI.v2 (Multivariate ENSO Index): Combines multiple atmospheric/oceanic variables
- Niño 3.4 anomalies: Raw SST anomalies
- TNA/TSA: Tropical North/South Atlantic indices (regional confounders)
- WHWP: Western Hemisphere Warm Pool

Outputs:
- enso_monthly.parquet: Monthly ENSO indices with phase classifications
- enso_daily.parquet: Daily interpolated values for merging with mortality data
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = Path(__file__).resolve().parent.parent.parent.parent
INPUT_DIR = BASE_DIR / "Input_data" / "data_causal_analysis" / "ocean_data"
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_oni():
    """Load Oceanic Niño Index data."""
    print("Loading ONI data...")
    df = pd.read_csv(INPUT_DIR / "oni.csv", skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    
    # Parse date column
    date_col = df.columns[0]
    df['date'] = pd.to_datetime(df[date_col])
    
    # Get ONI value column
    value_col = [c for c in df.columns if c != date_col and c != 'date'][0]
    df['oni'] = pd.to_numeric(df[value_col], errors='coerce')
    
    df = df[['date', 'oni']].dropna()
    df['year'] = df['date'].dt.year
    df['month'] = df['date'].dt.month
    
    print(f"  ONI: {len(df)} months, {df['date'].min().date()} to {df['date'].max().date()}")
    return df[['date', 'year', 'month', 'oni']]


def load_mei():
    """Load Multivariate ENSO Index v2 data."""
    print("Loading MEI.v2 data...")
    df = pd.read_csv(INPUT_DIR / "meiv2.csv", skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    
    # Parse date column
    date_col = df.columns[0]
    df['date'] = pd.to_datetime(df[date_col])
    
    # Get MEI value column
    value_col = [c for c in df.columns if c != date_col and c != 'date'][0]
    df['mei'] = pd.to_numeric(df[value_col], errors='coerce')
    df['mei'] = df['mei'].replace(-999, np.nan)
    
    df = df[['date', 'mei']].dropna()
    print(f"  MEI: {len(df)} months, {df['date'].min().date()} to {df['date'].max().date()}")
    return df


def load_nino34():
    """Load Niño 3.4 SST anomaly data."""
    print("Loading Niño 3.4 data...")
    df = pd.read_csv(INPUT_DIR / "nina34.anom.csv", skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    
    date_col = df.columns[0]
    df['date'] = pd.to_datetime(df[date_col])
    
    value_col = [c for c in df.columns if c != date_col and c != 'date'][0]
    df['nino34'] = pd.to_numeric(df[value_col], errors='coerce')
    df['nino34'] = df['nino34'].replace(-99.99, np.nan).replace(-999, np.nan)
    
    df = df[['date', 'nino34']].dropna()
    print(f"  Niño 3.4: {len(df)} months, {df['date'].min().date()} to {df['date'].max().date()}")
    return df


def load_tna():
    """Load Tropical North Atlantic index."""
    print("Loading TNA data...")
    df = pd.read_csv(INPUT_DIR / "tna.csv", skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    
    date_col = df.columns[0]
    df['date'] = pd.to_datetime(df[date_col])
    
    value_col = [c for c in df.columns if c != date_col and c != 'date'][0]
    df['tna'] = pd.to_numeric(df[value_col], errors='coerce')
    df['tna'] = df['tna'].replace(-99.99, np.nan).replace(-999, np.nan)
    
    df = df[['date', 'tna']].dropna()
    print(f"  TNA: {len(df)} months, {df['date'].min().date()} to {df['date'].max().date()}")
    return df


def load_tsa():
    """Load Tropical South Atlantic index."""
    print("Loading TSA data...")
    df = pd.read_csv(INPUT_DIR / "tsa.csv", skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    
    date_col = df.columns[0]
    df['date'] = pd.to_datetime(df[date_col])
    
    value_col = [c for c in df.columns if c != date_col and c != 'date'][0]
    df['tsa'] = pd.to_numeric(df[value_col], errors='coerce')
    df['tsa'] = df['tsa'].replace(-99.99, np.nan).replace(-999, np.nan)
    
    df = df[['date', 'tsa']].dropna()
    print(f"  TSA: {len(df)} months, {df['date'].min().date()} to {df['date'].max().date()}")
    return df


def load_whwp():
    """Load Western Hemisphere Warm Pool index."""
    print("Loading WHWP data...")
    df = pd.read_csv(INPUT_DIR / "whwp.csv", skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    
    date_col = df.columns[0]
    df['date'] = pd.to_datetime(df[date_col])
    
    value_col = [c for c in df.columns if c != date_col and c != 'date'][0]
    df['whwp'] = pd.to_numeric(df[value_col], errors='coerce')
    df['whwp'] = df['whwp'].replace(-99.99, np.nan).replace(-999, np.nan)
    
    df = df[['date', 'whwp']].dropna()
    print(f"  WHWP: {len(df)} months, {df['date'].min().date()} to {df['date'].max().date()}")
    return df


def classify_enso_phase(oni_value):
    """
    Classify ENSO phase based on ONI value.
    
    NOAA definitions:
    - El Niño: ONI >= 0.5 for 5+ consecutive overlapping 3-month periods
    - La Niña: ONI <= -0.5 for 5+ consecutive overlapping 3-month periods
    - Neutral: -0.5 < ONI < 0.5
    
    For monthly classification (simplified):
    - Strong El Niño: ONI >= 1.5
    - Moderate El Niño: 1.0 <= ONI < 1.5
    - Weak El Niño: 0.5 <= ONI < 1.0
    - Neutral: -0.5 < ONI < 0.5
    - Weak La Niña: -1.0 < ONI <= -0.5
    - Moderate La Niña: -1.5 < ONI <= -1.0
    - Strong La Niña: ONI <= -1.5
    """
    if pd.isna(oni_value):
        return 'unknown'
    elif oni_value >= 1.5:
        return 'strong_el_nino'
    elif oni_value >= 1.0:
        return 'moderate_el_nino'
    elif oni_value >= 0.5:
        return 'weak_el_nino'
    elif oni_value > -0.5:
        return 'neutral'
    elif oni_value > -1.0:
        return 'weak_la_nina'
    elif oni_value > -1.5:
        return 'moderate_la_nina'
    else:
        return 'strong_la_nina'


def classify_enso_binary(oni_value):
    """Binary classification for simpler analysis."""
    if pd.isna(oni_value):
        return np.nan
    elif oni_value >= 0.5:
        return 1  # El Niño
    elif oni_value <= -0.5:
        return -1  # La Niña
    else:
        return 0  # Neutral


def interpolate_to_daily(monthly_df, value_cols, start_date='2010-01-01', end_date='2024-12-31'):
    """
    Interpolate monthly values to daily for merging with mortality data.
    Uses linear interpolation between monthly midpoints.
    """
    print(f"Interpolating to daily ({start_date} to {end_date})...")
    
    # Create daily date range
    daily_dates = pd.date_range(start=start_date, end=end_date, freq='D')
    daily_df = pd.DataFrame({'date': daily_dates})
    
    # Set monthly dates to mid-month for interpolation
    monthly = monthly_df.copy()
    monthly['date'] = monthly['date'] + pd.Timedelta(days=14)  # Mid-month
    
    # Merge and interpolate
    daily_df = daily_df.merge(monthly[['date'] + value_cols], on='date', how='left')
    
    # Linear interpolation
    for col in value_cols:
        daily_df[col] = daily_df[col].interpolate(method='linear')
        # Forward/backward fill for edges
        daily_df[col] = daily_df[col].ffill().bfill()
    
    daily_df['year'] = daily_df['date'].dt.year
    daily_df['month'] = daily_df['date'].dt.month
    
    return daily_df


def main():
    """Main function to prepare ENSO data."""
    print("=" * 70)
    print("ENSO DATA PREPARATION")
    print("=" * 70)
    
    # Load all indices
    oni = load_oni()
    mei = load_mei()
    nino34 = load_nino34()
    tna = load_tna()
    tsa = load_tsa()
    whwp = load_whwp()
    
    # Merge all indices on date
    print("\nMerging all ENSO indices...")
    enso_monthly = oni.copy()
    
    for df, name in [(mei, 'mei'), (nino34, 'nino34'), (tna, 'tna'), (tsa, 'tsa'), (whwp, 'whwp')]:
        enso_monthly = enso_monthly.merge(df, on='date', how='left')
    
    # Classify ENSO phases
    print("Classifying ENSO phases...")
    enso_monthly['enso_phase'] = enso_monthly['oni'].apply(classify_enso_phase)
    enso_monthly['enso_binary'] = enso_monthly['oni'].apply(classify_enso_binary)
    
    # Create dummy variables for phases
    enso_monthly['is_el_nino'] = (enso_monthly['oni'] >= 0.5).astype(int)
    enso_monthly['is_la_nina'] = (enso_monthly['oni'] <= -0.5).astype(int)
    enso_monthly['is_strong_enso'] = ((enso_monthly['oni'] >= 1.5) | (enso_monthly['oni'] <= -1.5)).astype(int)
    
    # Save monthly data
    print(f"\nSaving monthly ENSO data...")
    enso_monthly.to_parquet(OUTPUT_DIR / "enso_monthly.parquet", index=False)
    enso_monthly.to_csv(OUTPUT_DIR / "enso_monthly.csv", index=False)
    
    # Filter to study period and create daily interpolation
    study_period = enso_monthly[(enso_monthly['year'] >= 2009) & (enso_monthly['year'] <= 2024)]
    
    value_cols = ['oni', 'mei', 'nino34', 'tna', 'tsa', 'whwp']
    enso_daily = interpolate_to_daily(study_period, value_cols)
    
    # Add phase classifications to daily data
    enso_daily['enso_phase'] = enso_daily['oni'].apply(classify_enso_phase)
    enso_daily['enso_binary'] = enso_daily['oni'].apply(classify_enso_binary)
    enso_daily['is_el_nino'] = (enso_daily['oni'] >= 0.5).astype(int)
    enso_daily['is_la_nina'] = (enso_daily['oni'] <= -0.5).astype(int)
    enso_daily['is_strong_enso'] = ((enso_daily['oni'] >= 1.5) | (enso_daily['oni'] <= -1.5)).astype(int)
    
    # Save daily data
    print(f"Saving daily ENSO data...")
    enso_daily.to_parquet(OUTPUT_DIR / "enso_daily.parquet", index=False)
    
    # Print summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    
    print("\nMonthly ENSO data shape:", enso_monthly.shape)
    print("\nENSO phase distribution (monthly, all years):")
    print(enso_monthly['enso_phase'].value_counts())
    
    print("\nENSO phase distribution (study period 2010-2024):")
    study_phases = enso_monthly[(enso_monthly['year'] >= 2010) & (enso_monthly['year'] <= 2024)]
    print(study_phases['enso_phase'].value_counts())
    
    print("\nDaily ENSO data shape:", enso_daily.shape)
    print(f"Date range: {enso_daily['date'].min().date()} to {enso_daily['date'].max().date()}")
    
    print("\nCorrelation between indices:")
    corr_cols = ['oni', 'mei', 'nino34', 'tna', 'tsa']
    print(study_period[corr_cols].corr().round(3))
    
    print("\n" + "=" * 70)
    print("ENSO data preparation complete!")
    print(f"Output files saved to: {OUTPUT_DIR}")
    print("=" * 70)
    
    return enso_monthly, enso_daily


if __name__ == "__main__":
    enso_monthly, enso_daily = main()
