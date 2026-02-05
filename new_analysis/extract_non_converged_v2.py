"""
Comprehensive extraction of non-converged regions across ALL analysis results.
Creates a dataset for re-analysis with different model specifications.
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

def add_missing_regions(result_file, exposure_type, analysis_type, successful_regions, all_regions_set):
    """Identify regions that are in the data but not in successful results."""
    missing = all_regions_set - set(successful_regions)
    issues = []
    for region in missing:
        issues.append({
            "source_file": str(result_file),
            "exposure_type": exposure_type,
            "analysis_type": analysis_type,
            "region_code": region,
            "issue": "missing_from_results",
            "reason": "Region did not converge or was excluded"
        })
    return issues

def process_dlnm_results(file_path, exposure_type):
    """Process DLNM result files to find non-converged regions."""
    issues = []
    
    if not file_path.exists():
        return issues
    
    with open(file_path, 'r') as f:
        data = json.load(f)
    
    # Get reference regions
    ref_regions = all_intermediate if exposure_type == "intermediate" else all_immediate
    
    # Check n_regions counts
    n_total = data.get("n_regions_total", 0)
    n_success = data.get("n_regions_success", 0)
    n_failed = n_total - n_success
    
    print(f"  {file_path.name}: {n_success}/{n_total} converged ({n_failed} failed)")
    
    # Get successful regions
    region_results = data.get("region_results", [])
    if isinstance(region_results, list):
        successful = [str(r.get("region_code", "")) for r in region_results]
    elif isinstance(region_results, dict):
        successful = list(region_results.keys())
    else:
        successful = []
    
    # Find missing regions
    missing = ref_regions - set(successful)
    for region in missing:
        issues.append({
            "source_file": str(file_path.relative_to(BASE_DIR)),
            "exposure_type": exposure_type,
            "analysis_type": "dlnm_first_stage",
            "region_code": region,
            "issue": "did_not_converge",
            "n_total": n_total,
            "n_success": n_success
        })
    
    return issues

def process_stratification_results(file_path, exposure_type, stratification_type):
    """Process stratification results (age, sex, cause) to find non-converged regions."""
    issues = []
    
    if not file_path.exists():
        return issues
    
    with open(file_path, 'r') as f:
        data = json.load(f)
    
    ref_regions = all_intermediate if exposure_type == "intermediate" else all_immediate
    
    # Check different possible structures
    strata_results = data.get("strata_results", data.get("results", {}))
    
    if isinstance(strata_results, dict):
        for stratum, stratum_data in strata_results.items():
            if isinstance(stratum_data, dict):
                # Check for regional results within stratum
                regional = stratum_data.get("regional_results", stratum_data.get("region_results", []))
                
                if isinstance(regional, list):
                    successful = [str(r.get("region_code", "")) for r in regional]
                elif isinstance(regional, dict):
                    successful = list(regional.keys())
                else:
                    continue
                
                n_success = len(successful)
                n_total = len(ref_regions)
                
                if n_success < n_total:
                    missing = ref_regions - set(successful)
                    for region in missing:
                        issues.append({
                            "source_file": str(file_path.relative_to(BASE_DIR)),
                            "exposure_type": exposure_type,
                            "analysis_type": f"{stratification_type}_{stratum}",
                            "region_code": region,
                            "issue": "did_not_converge",
                            "stratum": stratum,
                            "n_success": n_success,
                            "n_total": n_total
                        })
    
    elif isinstance(strata_results, list):
        for item in strata_results:
            if isinstance(item, dict):
                stratum = item.get("stratum", item.get("group", "unknown"))
                regional = item.get("regional_results", item.get("region_results", []))
                
                if isinstance(regional, list):
                    successful = [str(r.get("region_code", "")) for r in regional]
                elif isinstance(regional, dict):
                    successful = list(regional.keys())
                else:
                    continue
                
                n_success = len(successful)
                n_total = len(ref_regions)
                
                if n_success < n_total:
                    missing = ref_regions - set(successful)
                    for region in missing:
                        issues.append({
                            "source_file": str(file_path.relative_to(BASE_DIR)),
                            "exposure_type": exposure_type,
                            "analysis_type": f"{stratification_type}_{stratum}",
                            "region_code": region,
                            "issue": "did_not_converge",
                            "stratum": stratum,
                            "n_success": n_success,
                            "n_total": n_total
                        })
    
    return issues

# ============================================================
# Process Phase 1 R DLNM Results
# ============================================================
print("\n" + "="*60)
print("Phase 1: DLNM Core Models")
print("="*60)

# R DLNM results
non_converged.extend(process_dlnm_results(
    BASE_DIR / "phase1_r/results/dlnm_r_intermediate_results_v2.json", "intermediate"))
non_converged.extend(process_dlnm_results(
    BASE_DIR / "phase1_r/results/dlnm_r_immediate_results_v2.json", "immediate"))

# Python DLNM results
non_converged.extend(process_dlnm_results(
    BASE_DIR / "phase1_core_model/results/dlnm_v2_intermediate_results.json", "intermediate"))
non_converged.extend(process_dlnm_results(
    BASE_DIR / "phase1_core_model/results/dlnm_v2_immediate_results.json", "immediate"))

# ============================================================
# Process Phase 2 Robustness Results
# ============================================================
print("\n" + "="*60)
print("Phase 2: Robustness Checks")
print("="*60)

# Sensitivity analyses
for exposure in ["intermediate", "immediate"]:
    suffix = "" if exposure == "intermediate" else "_immediate"
    
    for analysis in ["sensitivity", "harvesting", "heatwave"]:
        r_file = BASE_DIR / f"phase2_r/results/{analysis}_r_{exposure}.json"
        if r_file.exists():
            try:
                with open(r_file) as f:
                    data = json.load(f)
                
                ref_regions = all_intermediate if exposure == "intermediate" else all_immediate
                
                # Check various result structures
                if "regional_results" in data:
                    regional = data["regional_results"]
                    if isinstance(regional, list):
                        successful = [str(r.get("region_code", "")) for r in regional]
                    elif isinstance(regional, dict):
                        successful = list(regional.keys())
                    else:
                        successful = []
                    
                    n_success = len(successful)
                    n_total = len(ref_regions)
                    
                    if n_success < n_total:
                        print(f"  {r_file.name}: {n_success}/{n_total} converged")
                        missing = ref_regions - set(successful)
                        for region in missing:
                            non_converged.append({
                                "source_file": str(r_file.relative_to(BASE_DIR)),
                                "exposure_type": exposure,
                                "analysis_type": analysis,
                                "region_code": region,
                                "issue": "did_not_converge",
                                "n_success": n_success,
                                "n_total": n_total
                            })
            except Exception as e:
                print(f"  Error processing {r_file}: {e}")

# ============================================================
# Process Phase 4 Heterogeneity Results
# ============================================================
print("\n" + "="*60)
print("Phase 4: Heterogeneity Analysis")
print("="*60)

# Process stratification results
for stratification in ["age", "sex", "cause"]:
    for exposure in ["intermediate", "immediate"]:
        suffix = "" if exposure == "intermediate" else "_immediate"
        
        # R results
        r_file = BASE_DIR / f"phase4_r/results/{stratification}_stratification_{exposure}.json"
        if r_file.exists():
            try:
                with open(r_file) as f:
                    data = json.load(f)
                
                ref_regions = all_intermediate if exposure == "intermediate" else all_immediate
                
                # Check strata
                if "strata" in data:
                    for stratum_data in data["strata"]:
                        stratum_name = stratum_data.get("stratum", stratum_data.get("group", "unknown"))
                        regional = stratum_data.get("regional_results", stratum_data.get("region_results", []))
                        
                        if isinstance(regional, list):
                            successful = [str(r.get("region_code", "")) for r in regional]
                        elif isinstance(regional, dict):
                            successful = list(regional.keys())
                        else:
                            successful = []
                        
                        n_success = len(successful)
                        n_total = len(ref_regions)
                        
                        if n_success < n_total:
                            print(f"  {r_file.name} [{stratum_name}]: {n_success}/{n_total} converged")
                            missing = ref_regions - set(successful)
                            for region in missing:
                                non_converged.append({
                                    "source_file": str(r_file.relative_to(BASE_DIR)),
                                    "exposure_type": exposure,
                                    "analysis_type": f"{stratification}_{stratum_name}",
                                    "region_code": region,
                                    "issue": "did_not_converge",
                                    "stratum": stratum_name,
                                    "n_success": n_success,
                                    "n_total": n_total
                                })
            except Exception as e:
                print(f"  Error processing {r_file}: {e}")
        
        # Python results
        py_file = BASE_DIR / f"phase4_heterogeneity/results/{stratification}_stratification_v2_results{suffix}.json"
        if py_file.exists():
            try:
                with open(py_file) as f:
                    data = json.load(f)
                
                ref_regions = all_intermediate if exposure == "intermediate" else all_immediate
                
                if "strata_results" in data:
                    strata = data["strata_results"]
                    for stratum_name, stratum_data in (strata.items() if isinstance(strata, dict) else []):
                        if isinstance(stratum_data, dict):
                            n_success = stratum_data.get("n_regions_success", 0)
                            n_total = stratum_data.get("n_regions_total", len(ref_regions))
                            
                            if n_success < n_total:
                                print(f"  {py_file.name} [{stratum_name}]: {n_success}/{n_total} converged")
                                
                                # Get regional results to find which ones are missing
                                regional = stratum_data.get("regional_results", stratum_data.get("region_results", []))
                                if isinstance(regional, list):
                                    successful = [str(r.get("region_code", "")) for r in regional]
                                elif isinstance(regional, dict):
                                    successful = list(regional.keys())
                                else:
                                    successful = []
                                
                                missing = ref_regions - set(successful)
                                for region in missing:
                                    non_converged.append({
                                        "source_file": str(py_file.relative_to(BASE_DIR)),
                                        "exposure_type": exposure,
                                        "analysis_type": f"{stratification}_{stratum_name}",
                                        "region_code": region,
                                        "issue": "did_not_converge",
                                        "stratum": stratum_name,
                                        "n_success": n_success,
                                        "n_total": n_total
                                    })
            except Exception as e:
                print(f"  Error processing {py_file}: {e}")

# ============================================================
# Create Output Dataset
# ============================================================
print("\n" + "="*60)
print("Creating Output Dataset")
print("="*60)

if non_converged:
    df = pd.DataFrame(non_converged)
    
    # Add region characteristics
    # Get mean deaths and temperature variance per region for context
    region_stats_int = mortality_int.groupby('region_code').agg({
        'deaths': ['sum', 'mean', 'std'],
    }).reset_index()
    region_stats_int.columns = ['region_code', 'total_deaths', 'mean_daily_deaths', 'std_deaths']
    region_stats_int['region_code'] = region_stats_int['region_code'].astype(str)
    region_stats_int['region_level'] = 'intermediate'
    
    region_stats_imm = mortality_imm.groupby('region_code').agg({
        'deaths': ['sum', 'mean', 'std'],
    }).reset_index()
    region_stats_imm.columns = ['region_code', 'total_deaths', 'mean_daily_deaths', 'std_deaths']
    region_stats_imm['region_code'] = region_stats_imm['region_code'].astype(str)
    region_stats_imm['region_level'] = 'immediate'
    
    region_stats = pd.concat([region_stats_int, region_stats_imm], ignore_index=True)
    
    # Merge region stats
    df = df.merge(region_stats, on='region_code', how='left')
    
    # Sort and clean
    df = df.sort_values(['exposure_type', 'analysis_type', 'region_code']).reset_index(drop=True)
    
    # Summary statistics
    print(f"\nTotal non-converged region-analysis pairs: {len(df)}")
    print(f"\nUnique regions that failed in at least one analysis:")
    print(f"  Intermediate: {df[df['exposure_type']=='intermediate']['region_code'].nunique()}")
    print(f"  Immediate: {df[df['exposure_type']=='immediate']['region_code'].nunique()}")
    
    print(f"\nBy Analysis Type:")
    for analysis, count in df['analysis_type'].value_counts().items():
        print(f"  {analysis}: {count}")
    
    print(f"\nMost problematic regions (failing in multiple analyses):")
    problem_regions = df.groupby(['region_code', 'exposure_type']).size().reset_index(name='n_failed_analyses')
    problem_regions = problem_regions.sort_values('n_failed_analyses', ascending=False)
    print(problem_regions.head(20).to_string())
    
    # Save full dataset
    output_path = BASE_DIR / "results" / "non_converged_regions.csv"
    output_path.parent.mkdir(exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"\nSaved full dataset to: {output_path}")
    
    # Save summary by unique region
    summary_df = df.groupby(['region_code', 'exposure_type']).agg({
        'analysis_type': lambda x: '; '.join(sorted(set(x))),
        'issue': 'count',
        'total_deaths': 'first',
        'mean_daily_deaths': 'first',
        'std_deaths': 'first'
    }).reset_index()
    summary_df.columns = ['region_code', 'exposure_type', 'failed_analyses', 'n_failures', 
                          'total_deaths', 'mean_daily_deaths', 'std_deaths']
    summary_df = summary_df.sort_values(['exposure_type', 'n_failures'], ascending=[True, False])
    
    summary_path = BASE_DIR / "results" / "non_converged_regions_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"Saved summary to: {summary_path}")
    
    # Display summary
    print("\n" + "="*60)
    print("Summary of Non-Converged Regions")
    print("="*60)
    print(summary_df.to_string())
    
else:
    print("\nNo non-converged regions found!")
    
    # Still create empty files for consistency
    output_path = BASE_DIR / "results" / "non_converged_regions.csv"
    output_path.parent.mkdir(exist_ok=True)
    pd.DataFrame(columns=['region_code', 'exposure_type', 'analysis_type', 'issue']).to_csv(output_path, index=False)
    
    summary_path = BASE_DIR / "results" / "non_converged_regions_summary.csv"
    pd.DataFrame(columns=['region_code', 'exposure_type', 'failed_analyses', 'n_failures']).to_csv(summary_path, index=False)

print("\n" + "="*60)
print("Script complete!")
print("="*60)
