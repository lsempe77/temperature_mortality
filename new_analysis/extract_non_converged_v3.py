"""
Comprehensive extraction of non-converged regions across ALL analysis results (v3).
Handles various JSON structures from R and Python outputs.
"""

import json
import os
import pandas as pd
from pathlib import Path
from datetime import datetime

BASE_DIR = Path(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis")

# Load reference data for all regions
print("Loading reference region lists...")
mortality_int = pd.read_parquet(BASE_DIR / "phase0_data_prep/results/mortality_intermediate_daily.parquet")
mortality_imm = pd.read_parquet(BASE_DIR / "phase0_data_prep/results/mortality_immediate_daily.parquet")

# Rename columns for consistency
mortality_int = mortality_int.rename(columns={'deaths_all': 'deaths'})
mortality_imm = mortality_imm.rename(columns={'immediate_code': 'region_code', 'deaths_all': 'deaths'})

all_intermediate = set(mortality_int['region_code'].unique().astype(str))
all_immediate = set(mortality_imm['region_code'].unique().astype(str))

print(f"Total intermediate regions in data: {len(all_intermediate)}")
print(f"Total immediate regions in data: {len(all_immediate)}")

# Store all non-converged results
non_converged = []

def find_missing_regions(successful_regions, ref_regions, source_file, exposure_type, analysis_type, stratum=None):
    """Identify missing regions and create issue records."""
    issues = []
    successful_set = set(str(r) for r in successful_regions)
    missing = ref_regions - successful_set
    
    for region in missing:
        issue = {
            "source_file": str(source_file),
            "exposure_type": exposure_type,
            "analysis_type": analysis_type,
            "region_code": region,
            "issue": "did_not_converge",
            "n_success": len(successful_set),
            "n_total": len(ref_regions)
        }
        if stratum:
            issue["stratum"] = stratum
        issues.append(issue)
    
    return issues

def extract_regions_from_results(data):
    """Extract region codes from various result structures."""
    regions = []
    
    if isinstance(data, list):
        for item in data:
            if isinstance(item, dict):
                rc = item.get('region_code', item.get('region_id', item.get('region', '')))
                if rc:
                    regions.append(str(rc))
    elif isinstance(data, dict):
        regions = [str(k) for k in data.keys() if k.isdigit() or (isinstance(k, str) and k.replace('-', '').isdigit())]
    
    return regions

# ============================================================
# Process Phase 1 DLNM Results
# ============================================================
print("\n" + "="*60)
print("Phase 1: DLNM Core Models")
print("="*60)

# R DLNM results
for exposure in ["intermediate", "immediate"]:
    r_file = BASE_DIR / f"phase1_r/results/dlnm_r_{exposure}_results_v2.json"
    ref_regions = all_intermediate if exposure == "intermediate" else all_immediate
    
    if r_file.exists():
        with open(r_file) as f:
            data = json.load(f)
        
        n_total = data.get("n_regions_total", 0)
        n_success = data.get("n_regions_success", 0)
        
        regional = data.get("region_results", [])
        successful = extract_regions_from_results(regional)
        
        n_failed = len(ref_regions) - len(successful)
        print(f"  {r_file.name}: {len(successful)}/{len(ref_regions)} converged ({n_failed} failed)")
        
        issues = find_missing_regions(successful, ref_regions, r_file.relative_to(BASE_DIR), 
                                       exposure, "dlnm_first_stage")
        non_converged.extend(issues)

# ============================================================
# Process Phase 4 Heterogeneity: Age Stratification
# ============================================================
print("\n" + "="*60)
print("Phase 4: Age Stratification")
print("="*60)

for exposure in ["intermediate", "immediate"]:
    r_file = BASE_DIR / f"phase4_r/results/age_stratification_{exposure}.json"
    ref_regions = all_intermediate if exposure == "intermediate" else all_immediate
    
    if r_file.exists():
        with open(r_file) as f:
            data = json.load(f)
        
        age_groups = data.get("age_groups", {})
        for age_group, age_data in age_groups.items():
            n_regions = age_data.get("n_regions", 0)
            converged = age_data.get("converged", True)
            
            n_failed = len(ref_regions) - n_regions
            if n_failed > 0:
                print(f"  {r_file.name} [{age_group}]: {n_regions}/{len(ref_regions)} converged ({n_failed} failed)")
                
                # We don't have individual region info, so mark as unknown which specific regions failed
                # Create placeholder issues
                for i in range(n_failed):
                    non_converged.append({
                        "source_file": str(r_file.relative_to(BASE_DIR)),
                        "exposure_type": exposure,
                        "analysis_type": f"age_{age_group}",
                        "region_code": f"unknown_{i+1}",
                        "issue": "did_not_converge",
                        "stratum": age_group,
                        "n_success": n_regions,
                        "n_total": len(ref_regions),
                        "note": "Individual region IDs not available in output"
                    })

# ============================================================
# Process Phase 4 Heterogeneity: Sex Stratification
# ============================================================
print("\n" + "="*60)
print("Phase 4: Sex Stratification")
print("="*60)

for exposure in ["intermediate", "immediate"]:
    r_file = BASE_DIR / f"phase4_r/results/sex_stratification_{exposure}.json"
    ref_regions = all_intermediate if exposure == "intermediate" else all_immediate
    
    if r_file.exists():
        with open(r_file) as f:
            data = json.load(f)
        
        sex_groups = data.get("sex_groups", {})
        for sex, sex_data in sex_groups.items():
            n_regions = sex_data.get("n_regions", 0)
            
            n_failed = len(ref_regions) - n_regions
            if n_failed > 0:
                print(f"  {r_file.name} [{sex}]: {n_regions}/{len(ref_regions)} converged ({n_failed} failed)")
                
                for i in range(n_failed):
                    non_converged.append({
                        "source_file": str(r_file.relative_to(BASE_DIR)),
                        "exposure_type": exposure,
                        "analysis_type": f"sex_{sex}",
                        "region_code": f"unknown_{i+1}",
                        "issue": "did_not_converge",
                        "stratum": sex,
                        "n_success": n_regions,
                        "n_total": len(ref_regions),
                        "note": "Individual region IDs not available in output"
                    })

# ============================================================
# Process Phase 4 Heterogeneity: Cause Stratification
# ============================================================
print("\n" + "="*60)
print("Phase 4: Cause Stratification")
print("="*60)

for exposure in ["intermediate", "immediate"]:
    r_file = BASE_DIR / f"phase4_r/results/cause_stratification_{exposure}.json"
    ref_regions = all_intermediate if exposure == "intermediate" else all_immediate
    
    if r_file.exists():
        with open(r_file) as f:
            data = json.load(f)
        
        cause_groups = data.get("cause_groups", {})
        for cause, cause_data in cause_groups.items():
            n_regions = cause_data.get("n_regions", 0)
            
            n_failed = len(ref_regions) - n_regions
            if n_failed > 0:
                print(f"  {r_file.name} [{cause}]: {n_regions}/{len(ref_regions)} converged ({n_failed} failed)")
                
                for i in range(n_failed):
                    non_converged.append({
                        "source_file": str(r_file.relative_to(BASE_DIR)),
                        "exposure_type": exposure,
                        "analysis_type": f"cause_{cause}",
                        "region_code": f"unknown_{i+1}",
                        "issue": "did_not_converge",
                        "stratum": cause,
                        "n_success": n_regions,
                        "n_total": len(ref_regions),
                        "note": "Individual region IDs not available in output"
                    })

# ============================================================
# Process Phase 2 Sensitivity Results
# ============================================================
print("\n" + "="*60)
print("Phase 2: Sensitivity Analysis")
print("="*60)

for exposure in ["intermediate", "immediate"]:
    r_file = BASE_DIR / f"phase2_r/results/sensitivity_r_{exposure}.json"
    ref_regions = all_intermediate if exposure == "intermediate" else all_immediate
    
    if r_file.exists():
        with open(r_file) as f:
            data = json.load(f)
        
        # Check various sensitivity tests
        for test_type in ["df_var", "lag_structure", "knot_placement"]:
            if test_type in data:
                test_data = data[test_type]
                if isinstance(test_data, dict):
                    for variant, variant_data in test_data.items():
                        if isinstance(variant_data, dict):
                            n_regions = variant_data.get("n_regions", variant_data.get("n_regions_success", 0))
                            
                            n_failed = len(ref_regions) - n_regions
                            if n_failed > 0:
                                print(f"  {r_file.name} [{test_type}/{variant}]: {n_regions}/{len(ref_regions)} converged ({n_failed} failed)")
                                
                                for i in range(n_failed):
                                    non_converged.append({
                                        "source_file": str(r_file.relative_to(BASE_DIR)),
                                        "exposure_type": exposure,
                                        "analysis_type": f"sensitivity_{test_type}_{variant}",
                                        "region_code": f"unknown_{i+1}",
                                        "issue": "did_not_converge",
                                        "stratum": f"{test_type}_{variant}",
                                        "n_success": n_regions,
                                        "n_total": len(ref_regions)
                                    })

# ============================================================
# Now identify SPECIFIC regions that failed in the DLNM first stage
# by comparing successful regions to the full set
# ============================================================
print("\n" + "="*60)
print("Identifying Specific Failed Regions from DLNM First Stage")
print("="*60)

dlnm_failed = []

# Intermediate
r_file = BASE_DIR / "phase1_r/results/dlnm_r_intermediate_results_v2.json"
with open(r_file) as f:
    data = json.load(f)
regional = data.get("region_results", [])
successful_int = set(str(r.get("region_code", "")) for r in regional)
failed_int = all_intermediate - successful_int
print(f"Failed intermediate regions: {sorted(failed_int)}")

for region in failed_int:
    dlnm_failed.append({
        "region_code": region,
        "exposure_type": "intermediate",
        "analysis_type": "dlnm_first_stage"
    })

# Immediate
r_file = BASE_DIR / "phase1_r/results/dlnm_r_immediate_results_v2.json"
with open(r_file) as f:
    data = json.load(f)
regional = data.get("region_results", [])
successful_imm = set(str(r.get("region_code", "")) for r in regional)
failed_imm = all_immediate - successful_imm
print(f"Failed immediate regions: {sorted(failed_imm)}")

for region in failed_imm:
    dlnm_failed.append({
        "region_code": region,
        "exposure_type": "immediate",
        "analysis_type": "dlnm_first_stage"
    })

# ============================================================
# Create Output Dataset with KNOWN failed regions
# ============================================================
print("\n" + "="*60)
print("Creating Output Dataset")
print("="*60)

# Create the main output with KNOWN failed regions only (from DLNM first stage)
df_known = pd.DataFrame(dlnm_failed)

# Add region characteristics
region_stats_int = mortality_int.groupby('region_code').agg({
    'deaths': ['sum', 'mean', 'std'],
}).reset_index()
region_stats_int.columns = ['region_code', 'total_deaths', 'mean_daily_deaths', 'std_deaths']
region_stats_int['region_code'] = region_stats_int['region_code'].astype(str)

region_stats_imm = mortality_imm.groupby('region_code').agg({
    'deaths': ['sum', 'mean', 'std'],
}).reset_index()
region_stats_imm.columns = ['region_code', 'total_deaths', 'mean_daily_deaths', 'std_deaths']
region_stats_imm['region_code'] = region_stats_imm['region_code'].astype(str)

region_stats = pd.concat([region_stats_int, region_stats_imm], ignore_index=True)

# Merge region stats
df_known = df_known.merge(region_stats, on='region_code', how='left')

# Load temperature data to add temperature stats for these regions
era5_int = pd.read_parquet(BASE_DIR / "phase0_data_prep/results/era5_intermediate_daily.parquet")
era5_imm = pd.read_parquet(BASE_DIR / "phase0_data_prep/results/era5_immediate_daily.parquet")

# Rename for consistency
if 'immediate_code' in era5_imm.columns:
    era5_imm = era5_imm.rename(columns={'immediate_code': 'region_code'})

# Get temp stats
temp_stats_int = era5_int.groupby('region_code').agg({
    'temp_mean': ['mean', 'std', 'min', 'max']
}).reset_index()
temp_stats_int.columns = ['region_code', 'temp_avg', 'temp_std', 'temp_min', 'temp_max']
temp_stats_int['region_code'] = temp_stats_int['region_code'].astype(str)

temp_stats_imm = era5_imm.groupby('region_code').agg({
    'temp_mean': ['mean', 'std', 'min', 'max']
}).reset_index()
temp_stats_imm.columns = ['region_code', 'temp_avg', 'temp_std', 'temp_min', 'temp_max']
temp_stats_imm['region_code'] = temp_stats_imm['region_code'].astype(str)

temp_stats = pd.concat([temp_stats_int, temp_stats_imm], ignore_index=True)

# Merge temp stats
df_known = df_known.merge(temp_stats, on='region_code', how='left')

# Add number of observations
obs_count_int = mortality_int.groupby('region_code').size().reset_index(name='n_days')
obs_count_int['region_code'] = obs_count_int['region_code'].astype(str)

obs_count_imm = mortality_imm.groupby('region_code').size().reset_index(name='n_days')
obs_count_imm['region_code'] = obs_count_imm['region_code'].astype(str)

obs_count = pd.concat([obs_count_int, obs_count_imm], ignore_index=True)
df_known = df_known.merge(obs_count, on='region_code', how='left')

# Sort
df_known = df_known.sort_values(['exposure_type', 'region_code']).reset_index(drop=True)

print(f"\nTotal KNOWN non-converged regions: {len(df_known)}")
print(f"  Intermediate: {len(df_known[df_known['exposure_type']=='intermediate'])}")
print(f"  Immediate: {len(df_known[df_known['exposure_type']=='immediate'])}")

# Save
output_path = BASE_DIR / "results" / "non_converged_regions.csv"
output_path.parent.mkdir(exist_ok=True)
df_known.to_csv(output_path, index=False)
print(f"\nSaved to: {output_path}")

# Display the data
print("\n" + "="*60)
print("Non-Converged Regions Dataset (for re-analysis)")
print("="*60)
print(df_known.to_string())

# Also save a summary of ALL analysis non-convergences (including unknowns)
df_all = pd.DataFrame(non_converged)
if len(df_all) > 0:
    all_summary = df_all.groupby(['exposure_type', 'analysis_type']).agg({
        'n_success': 'first',
        'n_total': 'first'
    }).reset_index()
    all_summary['n_failed'] = all_summary['n_total'] - all_summary['n_success']
    all_summary['pct_converged'] = (all_summary['n_success'] / all_summary['n_total'] * 100).round(1)
    
    print("\n" + "="*60)
    print("Convergence Summary Across All Analyses")
    print("="*60)
    print(all_summary.to_string())
    
    summary_path = BASE_DIR / "results" / "convergence_summary_all_analyses.csv"
    all_summary.to_csv(summary_path, index=False)
    print(f"\nSaved summary to: {summary_path}")

print("\n" + "="*60)
print("Script complete!")
print("="*60)
