"""
Extract non-converged regions from all JSON result files.
Creates a dataset of regions/models that failed to converge for re-analysis.
"""

import json
import os
import pandas as pd
from pathlib import Path
from datetime import datetime

# Base directory
BASE_DIR = Path(r"c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis")

# Collect all JSON files
json_files = list(BASE_DIR.rglob("*.json"))
print(f"Found {len(json_files)} JSON files")

# Store non-converged results
non_converged = []

def check_convergence(data, source_file, prefix=""):
    """Recursively check for convergence issues in nested data structures."""
    issues = []
    
    if isinstance(data, dict):
        # Check for explicit convergence flags
        if "converged" in data:
            if data["converged"] in [False, "false", "FALSE", 0, "0"]:
                issues.append({
                    "source_file": str(source_file),
                    "path": prefix,
                    "issue_type": "converged=False",
                    "details": str(data.get("convergence_message", data.get("message", "")))[:500]
                })
        
        # Check for convergence field
        if "convergence" in data:
            conv = data["convergence"]
            if isinstance(conv, dict):
                if conv.get("converged") in [False, "false", "FALSE", 0, "0"]:
                    issues.append({
                        "source_file": str(source_file),
                        "path": prefix,
                        "issue_type": "convergence.converged=False",
                        "details": str(conv.get("message", conv.get("reason", "")))[:500]
                    })
                if conv.get("status") in ["failed", "error", "not_converged"]:
                    issues.append({
                        "source_file": str(source_file),
                        "path": prefix,
                        "issue_type": f"convergence.status={conv.get('status')}",
                        "details": str(conv)[:500]
                    })
        
        # Check for error/warning flags
        if "error" in data and data["error"]:
            issues.append({
                "source_file": str(source_file),
                "path": prefix,
                "issue_type": "error=True",
                "details": str(data.get("error_message", data.get("error", "")))[:500]
            })
        
        # Check for failed status
        if "status" in data and data["status"] in ["failed", "error", "not_converged", "skipped"]:
            issues.append({
                "source_file": str(source_file),
                "path": prefix,
                "issue_type": f"status={data['status']}",
                "details": str(data.get("message", data.get("reason", "")))[:500]
            })
        
        # Check for fit issues
        if "fit_status" in data and data["fit_status"] not in ["success", "converged", True, 1]:
            issues.append({
                "source_file": str(source_file),
                "path": prefix,
                "issue_type": f"fit_status={data['fit_status']}",
                "details": str(data)[:500]
            })
        
        # Check for region-level results with convergence issues
        if "region_code" in data or "region_id" in data or "region" in data:
            region_id = data.get("region_code", data.get("region_id", data.get("region", "")))
            if "converged" in data and not data["converged"]:
                issues.append({
                    "source_file": str(source_file),
                    "path": prefix,
                    "issue_type": "region_not_converged",
                    "region_code": region_id,
                    "details": str(data)[:500]
                })
        
        # Check for model diagnostics showing issues
        if "diagnostics" in data:
            diag = data["diagnostics"]
            if isinstance(diag, dict):
                if diag.get("convergence_warning"):
                    issues.append({
                        "source_file": str(source_file),
                        "path": prefix + ".diagnostics",
                        "issue_type": "convergence_warning",
                        "details": str(diag.get("convergence_warning"))[:500]
                    })
                if diag.get("singular_hessian") or diag.get("hessian_singular"):
                    issues.append({
                        "source_file": str(source_file),
                        "path": prefix + ".diagnostics",
                        "issue_type": "singular_hessian",
                        "details": str(diag)[:500]
                    })
        
        # Check for warnings containing convergence-related text
        if "warnings" in data:
            warnings = data["warnings"]
            if isinstance(warnings, list):
                for i, w in enumerate(warnings):
                    w_str = str(w).lower()
                    if any(kw in w_str for kw in ["converg", "singular", "hessian", "not converged", "failed"]):
                        issues.append({
                            "source_file": str(source_file),
                            "path": f"{prefix}.warnings[{i}]",
                            "issue_type": "warning",
                            "details": str(w)[:500]
                        })
            elif isinstance(warnings, str):
                w_str = warnings.lower()
                if any(kw in w_str for kw in ["converg", "singular", "hessian", "not converged", "failed"]):
                    issues.append({
                        "source_file": str(source_file),
                        "path": f"{prefix}.warnings",
                        "issue_type": "warning",
                        "details": warnings[:500]
                    })
        
        # Recurse into nested structures
        for key, value in data.items():
            if key not in ["warnings", "diagnostics", "convergence"]:  # Already handled
                nested_issues = check_convergence(value, source_file, f"{prefix}.{key}" if prefix else key)
                issues.extend(nested_issues)
    
    elif isinstance(data, list):
        for i, item in enumerate(data):
            nested_issues = check_convergence(item, source_file, f"{prefix}[{i}]")
            issues.extend(nested_issues)
    
    return issues

# Process each JSON file
for json_file in json_files:
    try:
        with open(json_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        # Get relative path for cleaner output
        rel_path = json_file.relative_to(BASE_DIR)
        
        # Check for convergence issues
        issues = check_convergence(data, rel_path)
        non_converged.extend(issues)
        
        if issues:
            print(f"  {rel_path}: Found {len(issues)} issues")
        
    except json.JSONDecodeError as e:
        print(f"  Error parsing {json_file}: {e}")
    except Exception as e:
        print(f"  Error processing {json_file}: {e}")

print(f"\n{'='*60}")
print(f"Total convergence issues found: {len(non_converged)}")

# Create DataFrame
if non_converged:
    df = pd.DataFrame(non_converged)
    
    # Extract additional info from path
    df['phase'] = df['source_file'].apply(lambda x: x.split('\\')[0] if '\\' in x else x.split('/')[0])
    df['filename'] = df['source_file'].apply(lambda x: os.path.basename(x))
    
    # Clean up details
    df['details'] = df['details'].str.replace('\n', ' ').str.strip()
    
    # Sort by source file and path
    df = df.sort_values(['source_file', 'path']).reset_index(drop=True)
    
    # Save to CSV
    output_path = BASE_DIR / "results" / "non_converged_regions.csv"
    output_path.parent.mkdir(exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"\nSaved to: {output_path}")
    
    # Also save summary
    summary_path = BASE_DIR / "results" / "non_converged_summary.txt"
    with open(summary_path, 'w') as f:
        f.write(f"Non-Convergence Summary Report\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"{'='*60}\n\n")
        
        f.write(f"Total issues found: {len(non_converged)}\n\n")
        
        f.write("Issues by Phase:\n")
        for phase, count in df['phase'].value_counts().items():
            f.write(f"  {phase}: {count}\n")
        
        f.write(f"\nIssues by Type:\n")
        for issue_type, count in df['issue_type'].value_counts().items():
            f.write(f"  {issue_type}: {count}\n")
        
        f.write(f"\nIssues by File:\n")
        for file, count in df['filename'].value_counts().items():
            f.write(f"  {file}: {count}\n")
        
        f.write(f"\n{'='*60}\n")
        f.write("Detailed Issues:\n\n")
        for idx, row in df.iterrows():
            f.write(f"Issue #{idx+1}\n")
            f.write(f"  File: {row['source_file']}\n")
            f.write(f"  Path: {row['path']}\n")
            f.write(f"  Type: {row['issue_type']}\n")
            f.write(f"  Details: {row['details'][:200]}...\n\n")
    
    print(f"Saved summary to: {summary_path}")
    
    # Display summary
    print(f"\n{'='*60}")
    print("Issues by Phase:")
    print(df['phase'].value_counts().to_string())
    print(f"\nIssues by Type:")
    print(df['issue_type'].value_counts().to_string())
    
else:
    print("\nNo convergence issues found!")

# Also look for specific patterns in region results
print(f"\n{'='*60}")
print("Scanning for region-specific non-convergence patterns...")

# Look at DLNM results specifically for first-stage model issues
dlnm_files = [
    BASE_DIR / "phase1_r/results/dlnm_r_intermediate_results_v2.json",
    BASE_DIR / "phase1_r/results/dlnm_r_immediate_results_v2.json",
    BASE_DIR / "phase1_core_model/results/dlnm_v2_intermediate_results.json",
    BASE_DIR / "phase1_core_model/results/dlnm_v2_immediate_results.json"
]

region_issues = []

for dlnm_file in dlnm_files:
    if dlnm_file.exists():
        try:
            with open(dlnm_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            # Check for first_stage results
            if "first_stage" in data:
                first_stage = data["first_stage"]
                if isinstance(first_stage, dict):
                    # Check regional_results
                    if "regional_results" in first_stage:
                        for region_code, region_data in first_stage["regional_results"].items():
                            if isinstance(region_data, dict):
                                converged = region_data.get("converged", True)
                                if not converged:
                                    region_issues.append({
                                        "source_file": str(dlnm_file.relative_to(BASE_DIR)),
                                        "region_code": region_code,
                                        "model_type": "first_stage",
                                        "converged": converged,
                                        "n_obs": region_data.get("n_obs", region_data.get("observations")),
                                        "message": region_data.get("message", region_data.get("error", ""))
                                    })
                    
                    # Check results if it's a list
                    if "results" in first_stage and isinstance(first_stage["results"], list):
                        for r in first_stage["results"]:
                            if isinstance(r, dict) and not r.get("converged", True):
                                region_issues.append({
                                    "source_file": str(dlnm_file.relative_to(BASE_DIR)),
                                    "region_code": r.get("region_code", r.get("region_id", "unknown")),
                                    "model_type": "first_stage",
                                    "converged": False,
                                    "n_obs": r.get("n_obs"),
                                    "message": r.get("message", "")
                                })
            
            # Check region_results at top level
            if "region_results" in data:
                for region_code, region_data in data["region_results"].items():
                    if isinstance(region_data, dict):
                        if not region_data.get("converged", True):
                            region_issues.append({
                                "source_file": str(dlnm_file.relative_to(BASE_DIR)),
                                "region_code": region_code,
                                "model_type": "dlnm",
                                "converged": False,
                                "n_obs": region_data.get("n_obs"),
                                "message": region_data.get("message", "")
                            })
            
            print(f"  Checked {dlnm_file.name}")
            
        except Exception as e:
            print(f"  Error reading {dlnm_file}: {e}")

if region_issues:
    region_df = pd.DataFrame(region_issues)
    region_output = BASE_DIR / "results" / "non_converged_regions_dlnm.csv"
    region_df.to_csv(region_output, index=False)
    print(f"\nSaved {len(region_issues)} region-specific issues to: {region_output}")
    print(region_df.head(20).to_string())
else:
    print("No region-specific DLNM convergence issues found in first-stage models.")

print("\n" + "="*60)
print("Script complete!")
