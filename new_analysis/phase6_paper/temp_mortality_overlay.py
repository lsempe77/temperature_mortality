"""
Plot temperature and mortality trends overlaid to visualize associations.
Creates a dual-axis plot showing annual temperature metrics alongside mortality.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Paths
DATA_DIR = Path(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis\phase0_data_prep\results")
MORTALITY_DATA = Path(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis\phase1_r\results\attributable_burden_r_intermediate_annual.csv")
OUTPUT_DIR = Path(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis\phase6_paper")

print("Loading data...")

# Load mortality data (already have annual data)
mort_df = pd.read_csv(MORTALITY_DATA)
print(f"  Mortality data: {len(mort_df)} years ({mort_df['year'].min()}-{mort_df['year'].max()})")

# Load ERA5 temperature data (intermediate daily)
era5_df = pd.read_parquet(DATA_DIR / 'era5_intermediate_daily.parquet')
print(f"  ERA5 data: {len(era5_df)} rows")

# Check columns
print(f"  ERA5 columns: {list(era5_df.columns)}")

# Parse date and extract year
if 'date' in era5_df.columns:
    era5_df['date'] = pd.to_datetime(era5_df['date'])
    era5_df['year'] = era5_df['date'].dt.year
elif 'time' in era5_df.columns:
    era5_df['time'] = pd.to_datetime(era5_df['time'])
    era5_df['year'] = era5_df['time'].dt.year

# Find temperature column
temp_col = None
for col in ['temp_mean', 'tmean', 't2m', 'temperature', 'temp']:
    if col in era5_df.columns:
        temp_col = col
        break

if temp_col is None:
    print("Available columns:", era5_df.columns.tolist())
    raise ValueError("No temperature column found!")

print(f"  Using temperature column: {temp_col}")

# Calculate annual temperature statistics (national average across all regions)
annual_temp = era5_df.groupby('year').agg({
    temp_col: ['mean', 'std', 'min', 'max', 
               lambda x: np.percentile(x, 1),   # P1 (cold threshold)
               lambda x: np.percentile(x, 99),  # P99 (heat threshold)
               lambda x: (x > np.percentile(x, 97.5)).sum(),  # Number of heat days
               lambda x: (x < np.percentile(x, 2.5)).sum()]   # Number of cold days
}).reset_index()

# Flatten column names
annual_temp.columns = ['year', 'mean_temp', 'std_temp', 'min_temp', 'max_temp', 
                       'p1_temp', 'p99_temp', 'n_heat_days', 'n_cold_days']

print(f"\nAnnual temperature summary:")
print(annual_temp.head())

# Merge with mortality data
combined = pd.merge(mort_df, annual_temp, on='year', how='inner')
print(f"\nCombined data: {len(combined)} years")

# Create the visualization
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# ============================================================================
# Panel A: Heat mortality vs P99 temperature
# ============================================================================
ax1 = axes[0, 0]
ax1_twin = ax1.twinx()

# Heat deaths (left axis)
line1 = ax1.plot(combined['year'], combined['heat_an']/1000, 'o-', 
                  color='#E74C3C', linewidth=2.5, markersize=8, label='Heat deaths (thousands)')
ax1.fill_between(combined['year'], combined['heat_an']/1000, alpha=0.2, color='#E74C3C')
ax1.set_ylabel('Heat-Attributable Deaths (thousands)', color='#E74C3C', fontsize=11)
ax1.tick_params(axis='y', labelcolor='#E74C3C')

# P99 temperature (right axis)
line2 = ax1_twin.plot(combined['year'], combined['p99_temp'], 's--', 
                       color='#FF6B35', linewidth=2, markersize=7, label='P99 Temperature')
ax1_twin.set_ylabel('P99 Temperature (°C)', color='#FF6B35', fontsize=11)
ax1_twin.tick_params(axis='y', labelcolor='#FF6B35')

ax1.set_xlabel('Year', fontsize=11)
ax1.set_title('A. Heat Mortality vs P99 Temperature', fontsize=13, fontweight='bold')
ax1.grid(True, alpha=0.3)

# Combined legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax1_twin.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

# ============================================================================
# Panel B: Cold mortality vs P1 temperature
# ============================================================================
ax2 = axes[0, 1]
ax2_twin = ax2.twinx()

# Cold deaths (left axis)
line1 = ax2.plot(combined['year'], combined['cold_an']/1000, 'o-', 
                  color='#2E86AB', linewidth=2.5, markersize=8, label='Cold deaths (thousands)')
ax2.fill_between(combined['year'], combined['cold_an']/1000, alpha=0.2, color='#2E86AB')
ax2.set_ylabel('Cold-Attributable Deaths (thousands)', color='#2E86AB', fontsize=11)
ax2.tick_params(axis='y', labelcolor='#2E86AB')

# P1 temperature (right axis) - inverted because lower P1 = colder year
line2 = ax2_twin.plot(combined['year'], combined['p1_temp'], 's--', 
                       color='#1A5276', linewidth=2, markersize=7, label='P1 Temperature')
ax2_twin.set_ylabel('P1 Temperature (°C)', color='#1A5276', fontsize=11)
ax2_twin.tick_params(axis='y', labelcolor='#1A5276')
ax2_twin.invert_yaxis()  # Invert so lower = more cold exposure

ax2.set_xlabel('Year', fontsize=11)
ax2.set_title('B. Cold Mortality vs P1 Temperature (inverted)', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3)

lines1, labels1 = ax2.get_legend_handles_labels()
lines2, labels2 = ax2_twin.get_legend_handles_labels()
ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

# ============================================================================
# Panel C: Scatter plot - Heat deaths vs Mean Temperature
# ============================================================================
ax3 = axes[1, 0]
ax3.scatter(combined['mean_temp'], combined['heat_an']/1000, 
            s=100, c='#E74C3C', edgecolor='black', alpha=0.7)

# Add trend line
z = np.polyfit(combined['mean_temp'], combined['heat_an']/1000, 1)
p = np.poly1d(z)
x_line = np.linspace(combined['mean_temp'].min(), combined['mean_temp'].max(), 100)
ax3.plot(x_line, p(x_line), '--', color='#C0392B', linewidth=2)

# Calculate correlation
corr_heat = combined['mean_temp'].corr(combined['heat_an'])
ax3.set_xlabel('Mean Annual Temperature (°C)', fontsize=11)
ax3.set_ylabel('Heat-Attributable Deaths (thousands)', fontsize=11)
ax3.set_title(f'C. Heat Deaths vs Mean Temperature (r = {corr_heat:.2f})', fontsize=13, fontweight='bold')
ax3.grid(True, alpha=0.3)

# Add year labels
for _, row in combined.iterrows():
    ax3.annotate(str(int(row['year'])), (row['mean_temp'], row['heat_an']/1000),
                 textcoords="offset points", xytext=(5,5), fontsize=8, alpha=0.7)

# ============================================================================
# Panel D: Scatter plot - Cold deaths vs Mean Temperature  
# ============================================================================
ax4 = axes[1, 1]
ax4.scatter(combined['mean_temp'], combined['cold_an']/1000, 
            s=100, c='#2E86AB', edgecolor='black', alpha=0.7)

# Add trend line
z = np.polyfit(combined['mean_temp'], combined['cold_an']/1000, 1)
p = np.poly1d(z)
ax4.plot(x_line, p(x_line), '--', color='#1A5276', linewidth=2)

# Calculate correlation
corr_cold = combined['mean_temp'].corr(combined['cold_an'])
ax4.set_xlabel('Mean Annual Temperature (°C)', fontsize=11)
ax4.set_ylabel('Cold-Attributable Deaths (thousands)', fontsize=11)
ax4.set_title(f'D. Cold Deaths vs Mean Temperature (r = {corr_cold:.2f})', fontsize=13, fontweight='bold')
ax4.grid(True, alpha=0.3)

# Add year labels
for _, row in combined.iterrows():
    ax4.annotate(str(int(row['year'])), (row['mean_temp'], row['cold_an']/1000),
                 textcoords="offset points", xytext=(5,5), fontsize=8, alpha=0.7)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'temp_mortality_overlay.png', dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nSaved: temp_mortality_overlay.png")

# ============================================================================
# Additional: Simple time series overlay
# ============================================================================
fig2, ax = plt.subplots(figsize=(14, 6))

# Normalize for comparison (z-scores)
heat_z = (combined['heat_an'] - combined['heat_an'].mean()) / combined['heat_an'].std()
cold_z = (combined['cold_an'] - combined['cold_an'].mean()) / combined['cold_an'].std()
temp_z = (combined['mean_temp'] - combined['mean_temp'].mean()) / combined['mean_temp'].std()

ax.plot(combined['year'], heat_z, 'o-', color='#E74C3C', linewidth=2.5, markersize=8, label='Heat Deaths (z-score)')
ax.plot(combined['year'], cold_z, 's-', color='#2E86AB', linewidth=2.5, markersize=8, label='Cold Deaths (z-score)')
ax.plot(combined['year'], temp_z, '^--', color='#27AE60', linewidth=2, markersize=7, label='Mean Temperature (z-score)')

ax.axhline(y=0, color='gray', linestyle=':', linewidth=1)
ax.fill_between(combined['year'], heat_z, alpha=0.15, color='#E74C3C')
ax.fill_between(combined['year'], cold_z, alpha=0.15, color='#2E86AB')

ax.set_xlabel('Year', fontsize=12)
ax.set_ylabel('Standardized Value (z-score)', fontsize=12)
ax.set_title('Temperature and Mortality Trends (Standardized)', fontsize=14, fontweight='bold')
ax.legend(loc='upper left')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'temp_mortality_zscore.png', dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved: temp_mortality_zscore.png")

# Print correlation summary
print("\n" + "="*60)
print("CORRELATION SUMMARY")
print("="*60)
print(f"\nHeat deaths vs Mean temperature:  r = {corr_heat:.3f}")
print(f"Cold deaths vs Mean temperature:  r = {corr_cold:.3f}")
print(f"Heat deaths vs P99 temperature:   r = {combined['p99_temp'].corr(combined['heat_an']):.3f}")
print(f"Cold deaths vs P1 temperature:    r = {combined['p1_temp'].corr(combined['cold_an']):.3f}")
print(f"\nHeat deaths vs Cold deaths:       r = {combined['heat_an'].corr(combined['cold_an']):.3f}")

plt.show()
