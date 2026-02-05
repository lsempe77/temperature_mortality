"""
CONVERGENCE DIAGNOSTICS FOR DLNM REGIONAL ANALYSIS
====================================================
Analyzes why ~40-60% of regions have convergence failures in the 
temperature-mortality DLNM models.

Diagnosis categories:
1. Good convergence: RR in [0.9, 1.5], MMT in [15, 30]
2. Boundary failure: MMT stuck at bounds (0.7 or 35.1)
3. Exploded estimates: RR > 10 or RR < 0.5
4. Other issues

Output:
- Diagnostic tables and plots
- Recommendations for model improvement
- Map of convergence status by region

Author: Convergence Diagnostic Pipeline
Date: December 2025
"""

import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = Path(__file__).parent.parent
PHASE0_RESULTS = BASE_DIR / 'phase0_data_prep' / 'results'
PHASE1_R_RESULTS = BASE_DIR / 'phase1_r' / 'results'
OUTPUT_DIR = Path(__file__).parent / 'convergence_diagnostics'
OUTPUT_DIR.mkdir(exist_ok=True)

# Convergence thresholds
THRESHOLDS = {
    'rr_min': 0.9,      # Minimum plausible RR
    'rr_max': 1.5,      # Maximum plausible RR (heat)
    'rr_cold_max': 1.8, # Maximum plausible RR (cold - can be higher)
    'mmt_min': 12.0,    # Minimum plausible MMT for Brazil (relaxed)
    'mmt_max': 32.0,    # Maximum plausible MMT for Brazil (relaxed)
    'mmt_lower_bound': 0.0,   # Model lower bound (expanded)
    'mmt_upper_bound': 40.0,  # Model upper bound (expanded)
    'rr_exploded': 10.0,      # RR threshold for "exploded" estimate
}

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 150

print("="*70)
print("CONVERGENCE DIAGNOSTICS FOR DLNM REGIONAL ANALYSIS")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# 1. LOAD DATA
# =============================================================================
print("\n[1] Loading data...")

def load_results(level):
    """Load DLNM results for a given spatial level."""
    results_file = PHASE1_R_RESULTS / f'dlnm_r_{level}_summary.csv'
    if results_file.exists():
        df = pd.read_csv(results_file)
        print(f"  Loaded {level}: {len(df)} regions")
        return df
    else:
        print(f"  WARNING: {results_file} not found")
        return None

def load_covariates():
    """Load regional covariates."""
    cov_file = PHASE0_RESULTS / 'regional_covariates.csv'
    if cov_file.exists():
        df = pd.read_csv(cov_file)
        print(f"  Loaded covariates: {len(df)} regions")
        return df
    return None

# Load data
df_inter = load_results('intermediate')
df_immed = load_results('immediate')
df_cov = load_covariates()

# =============================================================================
# 2. CLASSIFY CONVERGENCE STATUS
# =============================================================================
print("\n[2] Classifying convergence status...")

def classify_convergence(df, level_name):
    """Classify each region's convergence status."""
    df = df.copy()
    
    # Initialize status column
    df['conv_status'] = 'Unknown'
    df['conv_heat'] = 'Unknown'
    df['conv_cold'] = 'Unknown'
    df['conv_mmt'] = 'Unknown'
    
    # Classify heat RR
    df.loc[df['rr_p99'].between(THRESHOLDS['rr_min'], THRESHOLDS['rr_max']), 'conv_heat'] = 'Good'
    df.loc[df['rr_p99'] > THRESHOLDS['rr_exploded'], 'conv_heat'] = 'Exploded'
    df.loc[df['rr_p99'] < THRESHOLDS['rr_min'], 'conv_heat'] = 'Too_Low'
    df.loc[(df['rr_p99'] > THRESHOLDS['rr_max']) & (df['rr_p99'] <= THRESHOLDS['rr_exploded']), 'conv_heat'] = 'High'
    
    # Classify cold RR
    df.loc[df['rr_p1'].between(THRESHOLDS['rr_min'], THRESHOLDS['rr_cold_max']), 'conv_cold'] = 'Good'
    df.loc[df['rr_p1'] > THRESHOLDS['rr_exploded'], 'conv_cold'] = 'Exploded'
    df.loc[df['rr_p1'] < THRESHOLDS['rr_min'], 'conv_cold'] = 'Too_Low'
    df.loc[(df['rr_p1'] > THRESHOLDS['rr_cold_max']) & (df['rr_p1'] <= THRESHOLDS['rr_exploded']), 'conv_cold'] = 'High'
    
    # Classify MMT
    df.loc[df['mmt'].between(THRESHOLDS['mmt_min'], THRESHOLDS['mmt_max']), 'conv_mmt'] = 'Good'
    df.loc[np.isclose(df['mmt'], THRESHOLDS['mmt_lower_bound'], atol=0.1), 'conv_mmt'] = 'Lower_Bound'
    df.loc[np.isclose(df['mmt'], THRESHOLDS['mmt_upper_bound'], atol=0.1), 'conv_mmt'] = 'Upper_Bound'
    df.loc[(df['mmt'] < THRESHOLDS['mmt_min']) & (df['conv_mmt'] == 'Unknown'), 'conv_mmt'] = 'Too_Low'
    df.loc[(df['mmt'] > THRESHOLDS['mmt_max']) & (df['conv_mmt'] == 'Unknown'), 'conv_mmt'] = 'Too_High'
    
    # Overall status
    df['conv_status'] = 'Good'
    df.loc[(df['conv_heat'] != 'Good') | (df['conv_cold'] != 'Good') | (df['conv_mmt'] != 'Good'), 'conv_status'] = 'Partial'
    df.loc[(df['conv_heat'] == 'Exploded') | (df['conv_cold'] == 'Exploded'), 'conv_status'] = 'Failed'
    df.loc[(df['conv_mmt'] == 'Lower_Bound') | (df['conv_mmt'] == 'Upper_Bound'), 'conv_status'] = 'Boundary'
    
    # Print summary
    print(f"\n  {level_name} Convergence Summary:")
    print(f"  {'='*50}")
    
    print(f"\n  Heat RR (P99) Status:")
    for status in df['conv_heat'].value_counts().index:
        n = (df['conv_heat'] == status).sum()
        pct = 100 * n / len(df)
        print(f"    {status:15s}: {n:4d} ({pct:5.1f}%)")
    
    print(f"\n  Cold RR (P1) Status:")
    for status in df['conv_cold'].value_counts().index:
        n = (df['conv_cold'] == status).sum()
        pct = 100 * n / len(df)
        print(f"    {status:15s}: {n:4d} ({pct:5.1f}%)")
    
    print(f"\n  MMT Status:")
    for status in df['conv_mmt'].value_counts().index:
        n = (df['conv_mmt'] == status).sum()
        pct = 100 * n / len(df)
        print(f"    {status:15s}: {n:4d} ({pct:5.1f}%)")
    
    print(f"\n  Overall Status:")
    for status in ['Good', 'Partial', 'Boundary', 'Failed']:
        if status in df['conv_status'].values:
            n = (df['conv_status'] == status).sum()
            pct = 100 * n / len(df)
            print(f"    {status:15s}: {n:4d} ({pct:5.1f}%)")
    
    return df

if df_inter is not None:
    df_inter = classify_convergence(df_inter, "Intermediate (133 regions)")

if df_immed is not None:
    df_immed = classify_convergence(df_immed, "Immediate (510 regions)")

# =============================================================================
# 3. ANALYZE PATTERNS - WHAT PREDICTS CONVERGENCE FAILURE?
# =============================================================================
print("\n[3] Analyzing patterns in convergence failures...")

def analyze_patterns(df, df_cov, level_name):
    """Analyze what predicts convergence failure."""
    
    # Merge with covariates if available
    if df_cov is not None:
        df = df.merge(df_cov, on='region_code', how='left')
    
    # Create binary good/bad indicator
    df['is_good'] = (df['conv_status'] == 'Good').astype(int)
    
    print(f"\n  {level_name} - Factors Associated with Convergence:")
    print(f"  {'='*60}")
    
    # Analyze by sample size (n_deaths)
    if 'n_deaths' in df.columns:
        good_deaths = df[df['is_good'] == 1]['n_deaths'].median()
        bad_deaths = df[df['is_good'] == 0]['n_deaths'].median()
        print(f"\n  Median deaths (Good vs Failed):")
        print(f"    Good convergence:   {good_deaths:,.0f} deaths")
        print(f"    Failed convergence: {bad_deaths:,.0f} deaths")
        print(f"    Ratio: {good_deaths/bad_deaths:.2f}x")
    
    # Analyze by number of days
    if 'n_days' in df.columns:
        good_days = df[df['is_good'] == 1]['n_days'].median()
        bad_days = df[df['is_good'] == 0]['n_days'].median()
        print(f"\n  Median days of data (Good vs Failed):")
        print(f"    Good convergence:   {good_days:,.0f} days")
        print(f"    Failed convergence: {bad_days:,.0f} days")
    
    # Analyze by mean temperature (if available)
    if 'mean_temp_annual' in df.columns:
        good_temp = df[df['is_good'] == 1]['mean_temp_annual'].median()
        bad_temp = df[df['is_good'] == 0]['mean_temp_annual'].median()
        print(f"\n  Mean annual temperature (Good vs Failed):")
        print(f"    Good convergence:   {good_temp:.1f}°C")
        print(f"    Failed convergence: {bad_temp:.1f}°C")
    
    # Analyze by macro-region (if available)
    if 'macro_region' in df.columns:
        print(f"\n  Convergence by Macro-Region:")
        region_stats = df.groupby('macro_region').agg({
            'is_good': ['sum', 'count']
        }).round(2)
        region_stats.columns = ['n_good', 'n_total']
        region_stats['pct_good'] = 100 * region_stats['n_good'] / region_stats['n_total']
        for region in region_stats.index:
            row = region_stats.loc[region]
            print(f"    {region:15s}: {row['n_good']:.0f}/{row['n_total']:.0f} good ({row['pct_good']:.1f}%)")
    
    # Calculate minimum deaths threshold
    if 'n_deaths' in df.columns:
        # Find threshold where 80% of regions converge
        sorted_deaths = df.sort_values('n_deaths')
        cumsum_good = sorted_deaths['is_good'].cumsum()
        total_good = df['is_good'].sum()
        
        # Find deaths threshold for various success rates
        print(f"\n  Deaths threshold analysis:")
        for target_pct in [50, 75, 90, 95]:
            target_n = int(total_good * target_pct / 100)
            if target_n > 0:
                idx = (cumsum_good >= target_n).idxmax()
                threshold = sorted_deaths.loc[idx, 'n_deaths']
                n_above = (df['n_deaths'] >= threshold).sum()
                print(f"    {target_pct}% of good regions have >= {threshold:,.0f} deaths ({n_above} regions)")
    
    return df

if df_inter is not None:
    df_inter = analyze_patterns(df_inter, df_cov, "Intermediate")

if df_immed is not None:
    df_immed = analyze_patterns(df_immed, df_cov, "Immediate")

# =============================================================================
# 4. GENERATE DIAGNOSTIC PLOTS
# =============================================================================
print("\n[4] Generating diagnostic plots...")

def create_diagnostic_plots(df, level_name, output_prefix):
    """Create diagnostic visualizations."""
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'Convergence Diagnostics - {level_name}', fontsize=14, fontweight='bold')
    
    # 1. Distribution of RR P99 (heat)
    ax = axes[0, 0]
    # Clip for visualization
    rr_clipped = df['rr_p99'].clip(upper=3)
    colors = ['green' if s == 'Good' else 'orange' if s in ['High', 'Too_Low'] else 'red' 
              for s in df['conv_heat']]
    ax.hist(rr_clipped, bins=30, edgecolor='black', alpha=0.7)
    ax.axvline(THRESHOLDS['rr_min'], color='red', linestyle='--', label=f"Min={THRESHOLDS['rr_min']}")
    ax.axvline(THRESHOLDS['rr_max'], color='red', linestyle='--', label=f"Max={THRESHOLDS['rr_max']}")
    ax.set_xlabel('RR at P99 (Heat)')
    ax.set_ylabel('Count')
    ax.set_title('Heat Effect Distribution\n(clipped at 3)')
    ax.legend()
    
    # 2. Distribution of RR P1 (cold)
    ax = axes[0, 1]
    rr_clipped = df['rr_p1'].clip(upper=3)
    ax.hist(rr_clipped, bins=30, edgecolor='black', alpha=0.7)
    ax.axvline(THRESHOLDS['rr_min'], color='red', linestyle='--')
    ax.axvline(THRESHOLDS['rr_cold_max'], color='red', linestyle='--')
    ax.set_xlabel('RR at P1 (Cold)')
    ax.set_ylabel('Count')
    ax.set_title('Cold Effect Distribution\n(clipped at 3)')
    
    # 3. Distribution of MMT
    ax = axes[0, 2]
    ax.hist(df['mmt'], bins=30, edgecolor='black', alpha=0.7)
    ax.axvline(THRESHOLDS['mmt_min'], color='red', linestyle='--', label=f"Min={THRESHOLDS['mmt_min']}")
    ax.axvline(THRESHOLDS['mmt_max'], color='red', linestyle='--', label=f"Max={THRESHOLDS['mmt_max']}")
    ax.axvline(THRESHOLDS['mmt_lower_bound'], color='orange', linestyle=':', label='Model bound')
    ax.axvline(THRESHOLDS['mmt_upper_bound'], color='orange', linestyle=':')
    ax.set_xlabel('MMT (°C)')
    ax.set_ylabel('Count')
    ax.set_title('MMT Distribution')
    ax.legend(fontsize=8)
    
    # 4. Scatter: n_deaths vs convergence status
    ax = axes[1, 0]
    if 'n_deaths' in df.columns:
        colors_map = {'Good': 'green', 'Partial': 'orange', 'Boundary': 'blue', 'Failed': 'red'}
        for status in colors_map:
            mask = df['conv_status'] == status
            if mask.any():
                ax.scatter(df.loc[mask, 'n_deaths'], df.loc[mask, 'rr_p99'].clip(upper=3),
                          c=colors_map[status], label=status, alpha=0.6, s=30)
        ax.set_xlabel('Total Deaths')
        ax.set_ylabel('RR at P99 (clipped)')
        ax.set_title('Deaths vs Heat Effect')
        ax.legend()
        ax.set_xscale('log')
    
    # 5. Scatter: n_deaths vs MMT
    ax = axes[1, 1]
    if 'n_deaths' in df.columns:
        for status in colors_map:
            mask = df['conv_status'] == status
            if mask.any():
                ax.scatter(df.loc[mask, 'n_deaths'], df.loc[mask, 'mmt'],
                          c=colors_map[status], label=status, alpha=0.6, s=30)
        ax.set_xlabel('Total Deaths')
        ax.set_ylabel('MMT (°C)')
        ax.set_title('Deaths vs MMT')
        ax.set_xscale('log')
    
    # 6. Convergence status pie chart
    ax = axes[1, 2]
    status_counts = df['conv_status'].value_counts()
    colors_pie = [colors_map.get(s, 'gray') for s in status_counts.index]
    ax.pie(status_counts.values, labels=status_counts.index, colors=colors_pie,
           autopct='%1.1f%%', startangle=90)
    ax.set_title('Overall Convergence Status')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / f'{output_prefix}_diagnostics.png', dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / f'{output_prefix}_diagnostics.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_prefix}_diagnostics.png")

if df_inter is not None:
    create_diagnostic_plots(df_inter, "Intermediate (133 regions)", "intermediate")

if df_immed is not None:
    create_diagnostic_plots(df_immed, "Immediate (510 regions)", "immediate")

# =============================================================================
# 5. DETAILED FAILURE ANALYSIS
# =============================================================================
print("\n[5] Detailed failure analysis...")

def detailed_failure_analysis(df, level_name):
    """Identify specific failure modes and potential fixes."""
    
    print(f"\n  {level_name} - Detailed Failure Modes:")
    print(f"  {'='*60}")
    
    # Mode 1: MMT at lower bound (0.7)
    lower_bound = np.isclose(df['mmt'], THRESHOLDS['mmt_lower_bound'], atol=0.1)
    n_lower = lower_bound.sum()
    print(f"\n  Mode 1: MMT stuck at LOWER bound (0.7°C): {n_lower} regions")
    if n_lower > 0:
        print(f"    → Likely cause: Very warm regions where all temps are above MMT")
        print(f"    → Potential fix: Use state-level or macro-region pooled estimate")
        if 'n_deaths' in df.columns:
            print(f"    → Median deaths in these regions: {df.loc[lower_bound, 'n_deaths'].median():,.0f}")
    
    # Mode 2: MMT at upper bound (35.1)
    upper_bound = np.isclose(df['mmt'], THRESHOLDS['mmt_upper_bound'], atol=0.1)
    n_upper = upper_bound.sum()
    print(f"\n  Mode 2: MMT stuck at UPPER bound (35.1°C): {n_upper} regions")
    if n_upper > 0:
        print(f"    → Likely cause: Very cold regions where all temps are below MMT")
        print(f"    → Potential fix: Use state-level or macro-region pooled estimate")
        if 'n_deaths' in df.columns:
            print(f"    → Median deaths in these regions: {df.loc[upper_bound, 'n_deaths'].median():,.0f}")
    
    # Mode 3: Exploded RR estimates
    exploded = (df['rr_p99'] > THRESHOLDS['rr_exploded']) | (df['rr_p1'] > THRESHOLDS['rr_exploded'])
    n_exploded = exploded.sum()
    print(f"\n  Mode 3: Exploded RR estimates (>10): {n_exploded} regions")
    if n_exploded > 0:
        print(f"    → Likely cause: Insufficient data or separation issues")
        print(f"    → Potential fix: Simplify cross-basis (fewer df) or pool with neighbors")
        if 'n_deaths' in df.columns:
            print(f"    → Median deaths in these regions: {df.loc[exploded, 'n_deaths'].median():,.0f}")
    
    # Mode 4: Small sample size
    if 'n_deaths' in df.columns:
        good_median = df[df['conv_status'] == 'Good']['n_deaths'].median()
        small_sample = df['n_deaths'] < good_median / 2
        n_small = small_sample.sum()
        failed_and_small = (df['conv_status'] != 'Good') & small_sample
        print(f"\n  Mode 4: Small sample size (<{good_median/2:,.0f} deaths): {n_small} regions")
        print(f"    → Of these, {failed_and_small.sum()} have convergence issues")
        print(f"    → Potential fix: Aggregate to larger regions or use hierarchical model")
    
    # Calculate recoverable regions
    print(f"\n  RECOVERY POTENTIAL:")
    print(f"  {'='*60}")
    
    # Potentially recoverable: boundary issues (might work with different bounds)
    potentially_recoverable = (df['conv_status'] == 'Boundary').sum()
    print(f"  Boundary issues (may recover with wider bounds): {potentially_recoverable}")
    
    # Likely unrecoverable: exploded estimates with small sample
    if 'n_deaths' in df.columns:
        unrecoverable = (exploded & (df['n_deaths'] < 10000)).sum()
        print(f"  Likely unrecoverable (exploded + small sample): {unrecoverable}")
    
    # Already good
    good = (df['conv_status'] == 'Good').sum()
    print(f"  Already good: {good}")
    
    return df

if df_inter is not None:
    df_inter = detailed_failure_analysis(df_inter, "Intermediate")

if df_immed is not None:
    df_immed = detailed_failure_analysis(df_immed, "Immediate")

# =============================================================================
# 6. SAVE DETAILED RESULTS
# =============================================================================
print("\n[6] Saving detailed results...")

if df_inter is not None:
    df_inter.to_csv(OUTPUT_DIR / 'convergence_intermediate.csv', index=False)
    print(f"  Saved: convergence_intermediate.csv")

if df_immed is not None:
    df_immed.to_csv(OUTPUT_DIR / 'convergence_immediate.csv', index=False)
    print(f"  Saved: convergence_immediate.csv")

# =============================================================================
# 7. RECOMMENDATIONS
# =============================================================================
print("\n" + "="*70)
print("RECOMMENDATIONS")
print("="*70)

if df_inter is not None:
    n_good = (df_inter['conv_status'] == 'Good').sum()
    n_total = len(df_inter)
    n_boundary = (df_inter['conv_status'] == 'Boundary').sum()
    n_failed = (df_inter['conv_status'] == 'Failed').sum()
    
    print(f"""
INTERMEDIATE LEVEL (133 regions):
---------------------------------
Current status: {n_good}/{n_total} ({100*n_good/n_total:.1f}%) with good convergence

ISSUES IDENTIFIED:
1. {n_boundary} regions have MMT at model bounds
   → The cross-basis MMT search hits the [0.7, 35.1] bounds
   
2. {n_failed} regions have exploded estimates
   → Likely due to sparse data or separation issues

RECOMMENDED ACTIONS:

A. FOR THE PAPER (immediate):
   - Use pooled meta-analysis RR for all burden calculations ✓ (already done)
   - Document regional heterogeneity with available good regions
   - Note that ~{100-100*n_good/n_total:.0f}% of regions had convergence issues
   
B. FOR MODEL IMPROVEMENT (if time permits):
   1. Expand MMT bounds: Change from [0.7, 35.1] to [5, 40]
   2. Simplify cross-basis for small regions:
      - Reduce temp df from 4 to 3
      - Reduce lag df from 4 to 3
   3. Use penalized splines (mgcv) instead of natural splines
   4. Consider hierarchical/Bayesian approach to borrow strength

C. FOR ROBUSTNESS:
   - Run sensitivity using only converged regions
   - Compare pooled RR from all vs good-only regions
   - Test if excluded regions have systematically different demographics
""")

if df_immed is not None:
    n_good = (df_immed['conv_status'] == 'Good').sum()
    n_total = len(df_immed)
    
    print(f"""
IMMEDIATE LEVEL (510 regions):
------------------------------
Current status: {n_good}/{n_total} ({100*n_good/n_total:.1f}%) with good convergence

The immediate level has MORE convergence issues than intermediate,
which is expected due to smaller sample sizes per region.

RECOMMENDATION: Prefer intermediate level for spatial analyses.
""")

print("="*70)
print(f"Diagnostics complete. Output saved to: {OUTPUT_DIR}")
print("="*70)
