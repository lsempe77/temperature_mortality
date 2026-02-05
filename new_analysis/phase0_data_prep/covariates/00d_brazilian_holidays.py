"""
00d: BRAZILIAN HOLIDAYS 2010-2024
==================================
Create a dataset of Brazilian national holidays for use as control variables.

Uses the `holidays` Python library to generate accurate holiday dates
for all years programmatically (handles Easter, Carnival, Corpus Christi etc.)

Holidays can affect mortality through:
- Changed healthcare access (clinics closed)
- Behavioral changes (travel, alcohol)
- Reporting delays

Output: 
- results/brazilian_holidays_list.csv (list of holidays)
- results/brazilian_holidays_daily.csv (daily indicator for all days)
- results/brazilian_holidays_daily.parquet
"""

import pandas as pd
from datetime import datetime
from pathlib import Path

print("="*70)
print("00d: BRAZILIAN HOLIDAYS 2010-2024")
print("="*70)

# =============================================================================
# DEPENDENCIES
# =============================================================================

try:
    import holidays as holidays_lib
    print("✓ holidays library installed")
except ImportError:
    print("✗ holidays library not installed")
    print("  Install with: pip install holidays")
    import sys
    sys.exit(1)

# =============================================================================
# CONFIGURATION
# =============================================================================

YEARS = list(range(2010, 2025))  # 2010-2024
# Output to phase0_data_prep/results/ (parent directory of covariates/)
OUTPUT_DIR = Path(__file__).parent.parent / 'results'
OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# GENERATE BRAZILIAN HOLIDAYS
# =============================================================================

print(f"\nGenerating holidays for {YEARS[0]}-{YEARS[-1]}...")

# Create Brazil holidays object for all years
br_holidays = holidays_lib.Brazil(years=YEARS)

# =============================================================================
# ADD CARNIVAL DAYS (Saturday through Tuesday before Ash Wednesday)
# =============================================================================
# Carnival is NOT a national holiday in Brazil (except in some states/cities),
# but it's a de facto holiday that affects mortality patterns significantly.
# We need to calculate it from Easter, as Carnival is 47 days before Easter Sunday.

from dateutil.easter import easter
from datetime import timedelta

print("\nAdding Carnival days (Sábado de Carnaval through Terça-feira de Carnaval)...")

carnival_records = []
for year in YEARS:
    easter_date = easter(year)
    # Ash Wednesday is 46 days before Easter
    ash_wednesday = easter_date - timedelta(days=46)
    
    # Carnival days:
    # - Sábado de Carnaval (Saturday): 4 days before Ash Wednesday
    # - Domingo de Carnaval (Sunday): 3 days before Ash Wednesday  
    # - Segunda-feira de Carnaval (Monday): 2 days before Ash Wednesday
    # - Terça-feira de Carnaval (Tuesday): 1 day before Ash Wednesday
    
    carnival_saturday = ash_wednesday - timedelta(days=4)
    carnival_sunday = ash_wednesday - timedelta(days=3)
    carnival_monday = ash_wednesday - timedelta(days=2)
    carnival_tuesday = ash_wednesday - timedelta(days=1)
    
    carnival_records.extend([
        {'date': carnival_saturday, 'name': 'Sábado de Carnaval'},
        {'date': carnival_sunday, 'name': 'Domingo de Carnaval'},
        {'date': carnival_monday, 'name': 'Segunda-feira de Carnaval'},
        {'date': carnival_tuesday, 'name': 'Terça-feira de Carnaval'},
    ])

print(f"  Added {len(carnival_records)} Carnival days ({len(YEARS)} years × 4 days)")

# Build holiday list from holidays library
holiday_records = []
for date, name in sorted(br_holidays.items()):
    holiday_records.append({
        'date': pd.Timestamp(date),  # Convert to pandas Timestamp
        'holiday_name_pt': name,
        'year': date.year,
        'month': date.month,
        'day': date.day,
        'day_of_week': date.strftime('%A'),
    })

# Add Carnival days (avoiding duplicates if library already includes some)
existing_dates = {pd.Timestamp(r['date']) for r in holiday_records}
for carnival in carnival_records:
    carnival_date = pd.Timestamp(carnival['date'])
    if carnival_date not in existing_dates:
        holiday_records.append({
            'date': carnival_date,
            'holiday_name_pt': carnival['name'],
            'year': carnival_date.year,
            'month': carnival_date.month,
            'day': carnival_date.day,
            'day_of_week': carnival_date.strftime('%A'),
        })
        existing_dates.add(carnival_date)

holidays_df = pd.DataFrame(holiday_records)
holidays_df['date'] = pd.to_datetime(holidays_df['date'])  # Ensure datetime type
holidays_df = holidays_df.sort_values('date').reset_index(drop=True)  # Sort by date
print(f"  Found {len(holidays_df)} holidays (including Carnival days)")

# Print summary by year
print("\nHolidays per year:")
for year in YEARS:
    count = len(holidays_df[holidays_df['year'] == year])
    print(f"  {year}: {count} holidays")

# =============================================================================
# CATEGORIZE HOLIDAYS
# =============================================================================

def categorize_holiday(name):
    """Categorize holiday by type for analysis."""
    name_lower = name.lower()
    if 'carnaval' in name_lower or 'cinzas' in name_lower:
        return 'carnival'
    elif 'sexta-feira santa' in name_lower or 'páscoa' in name_lower or 'paixão' in name_lower:
        return 'easter'
    elif 'natal' in name_lower or 'confraternização' in name_lower or 'ano novo' in name_lower:
        return 'year_end'
    elif 'corpus christi' in name_lower:
        return 'corpus_christi'
    else:
        return 'civic'

holidays_df['holiday_type'] = holidays_df['holiday_name_pt'].apply(categorize_holiday)
holidays_df['is_holiday'] = 1

print(f"\nHolidays by type:")
print(holidays_df['holiday_type'].value_counts())

# =============================================================================
# CREATE DAILY HOLIDAY INDICATOR
# =============================================================================

print("\nCreating daily holiday indicators...")

# Create full date range for 2010-2024
date_range = pd.date_range('2010-01-01', '2024-12-31', freq='D')
daily_df = pd.DataFrame({'date': date_range})

print(f"  Date range: {daily_df['date'].min().date()} to {daily_df['date'].max().date()}")
print(f"  Total days: {len(daily_df)}")

# Merge holidays
daily_df = daily_df.merge(
    holidays_df[['date', 'holiday_name_pt', 'holiday_type', 'is_holiday']], 
    on='date', 
    how='left'
)
daily_df['is_holiday'] = daily_df['is_holiday'].fillna(0).astype(int)

# Add day before/after holiday indicators
daily_df['is_day_before_holiday'] = daily_df['is_holiday'].shift(-1).fillna(0).astype(int)
daily_df['is_day_after_holiday'] = daily_df['is_holiday'].shift(1).fillna(0).astype(int)

# Add holiday week indicator (for carnival week, etc.)
daily_df['is_holiday_week'] = (
    daily_df['is_holiday'] | 
    daily_df['is_day_before_holiday'] | 
    daily_df['is_day_after_holiday']
).astype(int)

# Add year and month for convenience
daily_df['year'] = daily_df['date'].dt.year
daily_df['month'] = daily_df['date'].dt.month
daily_df['day_of_week'] = daily_df['date'].dt.dayofweek

print(f"\nDaily dataset summary:")
print(f"  Total days: {len(daily_df)}")
print(f"  Holiday days: {daily_df['is_holiday'].sum()}")
print(f"  Days before/after holiday: {daily_df['is_day_before_holiday'].sum() + daily_df['is_day_after_holiday'].sum()}")

# =============================================================================
# SAVE
# =============================================================================

print("\n" + "-"*70)
print("Saving output files...")
print("-"*70)

# Save holiday list
holidays_df.to_csv(OUTPUT_DIR / 'brazilian_holidays_list.csv', index=False)
print(f"✓ Saved: brazilian_holidays_list.csv ({len(holidays_df)} holidays)")

# Save daily indicators
daily_df.to_csv(OUTPUT_DIR / 'brazilian_holidays_daily.csv', index=False)
print(f"✓ Saved: brazilian_holidays_daily.csv ({len(daily_df)} days)")

# Also save as parquet for faster loading
daily_df.to_parquet(OUTPUT_DIR / 'brazilian_holidays_daily.parquet', index=False)
print(f"✓ Saved: brazilian_holidays_daily.parquet")

print(f"\n{'='*70}")
print("DONE!")
print(f"{'='*70}")
print("="*70)
