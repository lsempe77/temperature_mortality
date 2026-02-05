"""
Python wrapper for R-based MVMeta pooling.

Uses Gasparrini's dlnm and mvmeta R packages which are the gold standard
for distributed lag non-linear models and multivariate meta-analysis.
"""

import subprocess
import json
import tempfile
import os
from pathlib import Path
from typing import Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)

# Path to the R script
R_SCRIPT_PATH = Path(__file__).parent / "mvmeta_pooling.R"


def run_mvmeta_r(
    input_json_path: str,
    output_json_path: Optional[str] = None,
    r_executable: str = "Rscript"
) -> Dict[str, Any]:
    """
    Run MVMeta pooling using R's dlnm and mvmeta packages.
    
    Parameters:
    -----------
    input_json_path : str
        Path to JSON file with DLNM results (from Phase 1)
    output_json_path : str, optional
        Path for output JSON. If None, uses a temp file.
    r_executable : str
        Path to Rscript executable
        
    Returns:
    --------
    dict with pooled results
    """
    # Use temp file if no output path specified
    if output_json_path is None:
        fd, output_json_path = tempfile.mkstemp(suffix=".json")
        os.close(fd)
        cleanup_output = True
    else:
        cleanup_output = False
    
    try:
        # Run R script
        cmd = [
            r_executable,
            str(R_SCRIPT_PATH),
            input_json_path,
            output_json_path
        ]
        
        logger.info(f"Running: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        if result.returncode != 0:
            logger.error(f"R script failed: {result.stderr}")
            raise RuntimeError(f"R script failed: {result.stderr}")
        
        logger.info(f"R output: {result.stdout}")
        
        # Read results
        with open(output_json_path, 'r') as f:
            output = json.load(f)
        
        return output
        
    finally:
        if cleanup_output and os.path.exists(output_json_path):
            os.remove(output_json_path)


def check_r_packages() -> bool:
    """Check if required R packages are installed."""
    check_script = '''
    packages <- c("dlnm", "mvmeta", "jsonlite", "splines")
    missing <- packages[!packages %in% installed.packages()[,"Package"]]
    if (length(missing) > 0) {
        cat("Missing packages:", paste(missing, collapse=", "), "\\n")
        quit(status=1)
    }
    cat("All packages installed\\n")
    '''
    
    result = subprocess.run(
        ["Rscript", "-e", check_script],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        logger.warning(f"R package check: {result.stdout}{result.stderr}")
        return False
    return True


def install_r_packages() -> bool:
    """Install required R packages."""
    install_script = '''
    packages <- c("dlnm", "mvmeta", "jsonlite")
    for (pkg in packages) {
        if (!require(pkg, character.only=TRUE, quietly=TRUE)) {
            install.packages(pkg, repos="https://cloud.r-project.org")
        }
    }
    '''
    
    result = subprocess.run(
        ["Rscript", "-e", install_script],
        capture_output=True,
        text=True,
        timeout=600
    )
    
    return result.returncode == 0


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python mvmeta_r_wrapper.py <input_json> [output_json]")
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else None
    
    # Check R packages
    print("Checking R packages...")
    if not check_r_packages():
        print("Installing R packages...")
        install_r_packages()
    
    # Run pooling
    print(f"Running MVMeta on {input_path}...")
    result = run_mvmeta_r(input_path, output_path)
    
    # Print summary
    print("\n=== Results ===")
    print(f"Converged: {result['mvmeta']['converged']}")
    print(f"I2: {result['mvmeta']['I2']:.1f}%")
    
    rr = result['pooled_rr']
    print(f"\nHeat P99 vs P50: RR={rr['heat_p99_vs_p50']['rr']:.4f} "
          f"({rr['heat_p99_vs_p50']['rr_lo']:.4f}-{rr['heat_p99_vs_p50']['rr_hi']:.4f})")
    print(f"Cold P1 vs P50:  RR={rr['cold_p1_vs_p50']['rr']:.4f} "
          f"({rr['cold_p1_vs_p50']['rr_lo']:.4f}-{rr['cold_p1_vs_p50']['rr_hi']:.4f})")
