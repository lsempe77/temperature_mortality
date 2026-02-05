"""
Plot the trend of cold-to-heat attributable deaths ratio over time (2010-2024)
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the annual attributable burden data
data_path = r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis\phase1_r\results\attributable_burden_r_intermediate_annual.csv"
df = pd.read_csv(data_path)

# Calculate the cold-to-heat ratio
df['cold_heat_ratio'] = df['cold_an'] / df['heat_an']

# Calculate the percentage of cold in total temperature-attributable deaths
df['cold_pct'] = df['cold_an'] / df['total_an'] * 100
df['heat_pct'] = df['heat_an'] / df['total_an'] * 100

# Create figure with multiple panels
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Calculate the overall ratio from totals (this is what the manuscript uses)
overall_ratio = df['cold_an'].sum() / df['heat_an'].sum()

# Panel A: Cold-to-Heat Ratio over time
ax1 = axes[0, 0]
ax1.plot(df['year'], df['cold_heat_ratio'], 'o-', color='#2E86AB', linewidth=2, markersize=8)
ax1.axhline(y=overall_ratio, color='#E74C3C', linestyle='--', linewidth=1.5, 
            label=f'Overall ratio (total cold/total heat): {overall_ratio:.2f}:1')
ax1.fill_between(df['year'], df['cold_heat_ratio'], alpha=0.2, color='#2E86AB')
ax1.set_xlabel('Year', fontsize=12)
ax1.set_ylabel('Cold:Heat Ratio', fontsize=12)
ax1.set_title('A. Cold-to-Heat Attributable Deaths Ratio Over Time', fontsize=14, fontweight='bold')
ax1.legend(loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.set_xticks(df['year'])
ax1.set_xticklabels(df['year'], rotation=45)

# Panel B: Stacked bar chart of heat and cold deaths
ax2 = axes[0, 1]
width = 0.7
ax2.bar(df['year'], df['heat_an']/1000, width, label='Heat-attributable', color='#E74C3C', alpha=0.8)
ax2.bar(df['year'], df['cold_an']/1000, width, bottom=df['heat_an']/1000, label='Cold-attributable', color='#2E86AB', alpha=0.8)
ax2.set_xlabel('Year', fontsize=12)
ax2.set_ylabel('Attributable Deaths (thousands)', fontsize=12)
ax2.set_title('B. Annual Attributable Deaths by Type', fontsize=14, fontweight='bold')
ax2.legend(loc='upper left')
ax2.grid(True, alpha=0.3, axis='y')
ax2.set_xticks(df['year'])
ax2.set_xticklabels(df['year'], rotation=45)

# Panel C: Percentage contribution (stacked area)
ax3 = axes[1, 0]
ax3.stackplot(df['year'], df['heat_pct'], df['cold_pct'], 
              labels=['Heat (%)', 'Cold (%)'],
              colors=['#E74C3C', '#2E86AB'], alpha=0.8)
ax3.axhline(y=50, color='white', linestyle='--', linewidth=1)
ax3.set_xlabel('Year', fontsize=12)
ax3.set_ylabel('Percentage of Total Burden', fontsize=12)
ax3.set_title('C. Relative Contribution of Heat vs Cold', fontsize=14, fontweight='bold')
ax3.legend(loc='center right')
ax3.set_ylim(0, 100)
ax3.set_xticks(df['year'])
ax3.set_xticklabels(df['year'], rotation=45)
ax3.grid(True, alpha=0.3, axis='y')

# Panel D: Summary statistics table
ax4 = axes[1, 1]
ax4.axis('off')

# Calculate summary statistics
summary_stats = {
    'Metric': ['Mean Cold:Heat Ratio', 'Min Ratio (Year)', 'Max Ratio (Year)', 
               'Mean Cold Deaths/Year', 'Mean Heat Deaths/Year', 
               'Mean Cold %', 'Mean Heat %'],
    'Value': [
        f'{df["cold_heat_ratio"].mean():.2f}:1',
        f'{df["cold_heat_ratio"].min():.2f}:1 ({df.loc[df["cold_heat_ratio"].idxmin(), "year"]})',
        f'{df["cold_heat_ratio"].max():.2f}:1 ({df.loc[df["cold_heat_ratio"].idxmax(), "year"]})',
        f'{df["cold_an"].mean():,.0f}',
        f'{df["heat_an"].mean():,.0f}',
        f'{df["cold_pct"].mean():.1f}%',
        f'{df["heat_pct"].mean():.1f}%'
    ]
}
summary_df = pd.DataFrame(summary_stats)

# Create table
table = ax4.table(cellText=summary_df.values,
                  colLabels=summary_df.columns,
                  cellLoc='left',
                  loc='center',
                  colWidths=[0.5, 0.5])
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1.2, 1.8)
# Style header
for i in range(len(summary_df.columns)):
    table[(0, i)].set_facecolor('#2E86AB')
    table[(0, i)].set_text_props(color='white', fontweight='bold')
ax4.set_title('D. Summary Statistics', fontsize=14, fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis\phase6_paper\cold_heat_ratio_trend.png", 
            dpi=300, bbox_inches='tight', facecolor='white')
plt.show()

# Print detailed yearly data
print("\n" + "="*80)
print("COLD-TO-HEAT ATTRIBUTABLE DEATHS RATIO BY YEAR")
print("="*80)
print(f"\n{'Year':<8} {'Heat Deaths':>15} {'Cold Deaths':>15} {'Ratio (Cold:Heat)':>20} {'Cold %':>10}")
print("-"*70)
for _, row in df.iterrows():
    print(f"{int(row['year']):<8} {row['heat_an']:>15,.0f} {row['cold_an']:>15,.0f} {row['cold_heat_ratio']:>20.2f}:1 {row['cold_pct']:>9.1f}%")

print("-"*70)
print(f"{'MEAN':<8} {df['heat_an'].mean():>15,.0f} {df['cold_an'].mean():>15,.0f} {df['cold_heat_ratio'].mean():>20.2f}:1 {df['cold_pct'].mean():>9.1f}%")
print(f"{'TOTAL':<8} {df['heat_an'].sum():>15,.0f} {df['cold_an'].sum():>15,.0f}")

# Key distinction: two ways to calculate the ratio
total_cold = df['cold_an'].sum()
total_heat = df['heat_an'].sum()
ratio_from_totals = total_cold / total_heat

print("\n" + "="*80)
print("KEY INSIGHT: Two ways to calculate the cold:heat ratio")
print("="*80)
print(f"\n1. Mean of yearly ratios:        {df['cold_heat_ratio'].mean():.2f}:1")
print(f"2. Ratio of total sums:          {ratio_from_totals:.2f}:1  <-- This is the 2.7:1 in the manuscript")
print(f"\n   Total cold deaths (2010-2024): {total_cold:,.0f}")
print(f"   Total heat deaths (2010-2024): {total_heat:,.0f}")
print(f"   {total_cold:,.0f} / {total_heat:,.0f} = {ratio_from_totals:.2f}")
print("\nThe manuscript uses method 2 (ratio of totals), not the mean of yearly ratios.")
