#!/usr/bin/env python3
"""
Build script for Temperature-Mortality Paper
=============================================

Orchestrates the complete build process:
1. Validates upstream data availability
2. Runs Quarto render for the master document
3. Generates provenance report

Usage:
    python build.py                    # Full build
    python build.py --validate-only    # Just check data availability
    python build.py --paper-only       # Just render academic_paper.qmd
"""

import subprocess
import sys
import json
from pathlib import Path
from datetime import datetime
import argparse


def check_quarto():
    """Check if Quarto is available."""
    try:
        result = subprocess.run(['quarto', '--version'], 
                               capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            print(f"✓ Quarto version: {result.stdout.strip()}")
            return True
    except:
        pass
    print("✗ Quarto not found. Install from https://quarto.org/")
    return False


def validate_data():
    """Validate upstream data availability."""
    print("\n" + "="*60)
    print("VALIDATING UPSTREAM DATA")
    print("="*60 + "\n")
    
    # Navigate to project root
    script_dir = Path(__file__).parent
    project_root = script_dir.parent.parent
    
    checks = {
        'Phase 0 Results': project_root / 'phase0_data_prep' / 'results',
        'Phase 1 Results': project_root / 'phase1_core_model' / 'results',
        'Phase 2 Results': project_root / 'phase2_robustness' / 'results',
        'Phase 3 Results': project_root / 'phase3_confounding' / 'results',
        'Phase 4 Results': project_root / 'phase4_heterogeneity' / 'results',
        'Phase 5 v2 Figures': project_root / 'phase5_outputs' / 'v2' / 'figures',
        'Tables': project_root / 'tables',
        'Bibliography': script_dir / 'references.bib',
        'CSL Style': script_dir / 'lancet-planetary-health.csl',
    }
    
    all_ok = True
    for name, path in checks.items():
        exists = path.exists()
        status = "✓" if exists else "✗"
        print(f"  {status} {name}: {path}")
        if not exists:
            all_ok = False
    
    # Check specific critical files
    critical_files = [
        project_root / 'phase1_core_model' / 'results' / 'dlnm_immediate_pooled.json',
        project_root / 'phase1_core_model' / 'results' / 'attributable_burden_national.json',
    ]
    
    print("\nCritical files:")
    for f in critical_files:
        exists = f.exists()
        status = "✓" if exists else "✗"
        print(f"  {status} {f.name}")
        if not exists:
            all_ok = False
    
    return all_ok


def render_master():
    """Render the master document."""
    print("\n" + "="*60)
    print("RENDERING MASTER DOCUMENT")
    print("="*60 + "\n")
    
    script_dir = Path(__file__).parent
    
    result = subprocess.run(
        ['quarto', 'render', '_master_document.qmd', '--to', 'html'],
        cwd=script_dir,
        capture_output=False  # Show output in real-time
    )
    
    return result.returncode == 0


def render_paper():
    """Render just the academic paper."""
    print("\n" + "="*60)
    print("RENDERING ACADEMIC PAPER")
    print("="*60 + "\n")
    
    script_dir = Path(__file__).parent
    
    result = subprocess.run(
        ['quarto', 'render', 'academic_paper.qmd'],
        cwd=script_dir,
        capture_output=False
    )
    
    return result.returncode == 0


def save_build_report(success: bool):
    """Save build provenance report."""
    script_dir = Path(__file__).parent
    
    # Get git info
    try:
        commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'],
                               capture_output=True, text=True, timeout=5)
        git_commit = commit.stdout.strip() if commit.returncode == 0 else "no-git"
        
        branch = subprocess.run(['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
                               capture_output=True, text=True, timeout=5)
        git_branch = branch.stdout.strip() if branch.returncode == 0 else "unknown"
    except:
        git_commit = "no-git"
        git_branch = "unknown"
    
    report = {
        'build_date': datetime.now().isoformat(),
        'success': success,
        'git_commit': git_commit,
        'git_branch': git_branch,
        'quarto_version': subprocess.run(['quarto', '--version'], 
                                         capture_output=True, text=True).stdout.strip(),
    }
    
    report_path = script_dir / 'build_report.json'
    report_path.write_text(json.dumps(report, indent=2))
    print(f"\nBuild report saved: {report_path}")


def main():
    parser = argparse.ArgumentParser(description="Build Temperature-Mortality Paper")
    parser.add_argument('--validate-only', action='store_true',
                       help="Only validate data availability")
    parser.add_argument('--paper-only', action='store_true',
                       help="Only render academic_paper.qmd")
    parser.add_argument('--skip-validation', action='store_true',
                       help="Skip data validation")
    args = parser.parse_args()
    
    print("="*60)
    print("TEMPERATURE-MORTALITY PAPER BUILD")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*60)
    
    # Check Quarto
    if not check_quarto():
        sys.exit(1)
    
    # Validate data
    if not args.skip_validation:
        data_ok = validate_data()
        if not data_ok:
            print("\n⚠️  Some data sources missing. Build may fail or produce incomplete output.")
    
    if args.validate_only:
        print("\n✓ Validation complete")
        sys.exit(0)
    
    # Render
    if args.paper_only:
        success = render_paper()
    else:
        success = render_master()
    
    # Save report
    save_build_report(success)
    
    if success:
        print("\n" + "="*60)
        print("✓ BUILD COMPLETE")
        print("="*60)
    else:
        print("\n" + "="*60)
        print("✗ BUILD FAILED")
        print("="*60)
        sys.exit(1)


if __name__ == "__main__":
    main()
