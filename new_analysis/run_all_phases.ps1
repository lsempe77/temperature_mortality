# =============================================================================
# run_all_phases.ps1
# Complete Analysis Pipeline - Phase 1 to 4
# =============================================================================
# Usage:
#   .\run_all_phases.ps1                    # Run both intermediate and immediate
#   .\run_all_phases.ps1 -Level intermediate # Run only intermediate
#   .\run_all_phases.ps1 -Level immediate    # Run only immediate
# =============================================================================

param(
    [ValidateSet("intermediate", "immediate", "both")]
    [string]$Level = "both"
)

$ErrorActionPreference = "Continue"
$RSCRIPT = "C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
$BASE_DIR = Split-Path -Parent $MyInvocation.MyCommand.Path

# Determine which levels to run
$levels = @()
if ($Level -eq "both") {
    $levels = @("intermediate", "immediate")
} else {
    $levels = @($Level)
}

Write-Host "=============================================================" -ForegroundColor Cyan
Write-Host "COMPLETE ANALYSIS PIPELINE - PHASES 1-4" -ForegroundColor Cyan
Write-Host "=============================================================" -ForegroundColor Cyan
Write-Host "Started: $(Get-Date)" -ForegroundColor Gray
Write-Host "Levels to run: $($levels -join ', ')" -ForegroundColor Gray
Write-Host ""

# Track timing
$startTime = Get-Date
$results = @{}

foreach ($LEVEL in $levels) {
    Write-Host ""
    Write-Host "=============================================================" -ForegroundColor Yellow
    Write-Host "RUNNING ALL PHASES FOR: $LEVEL" -ForegroundColor Yellow
    Write-Host "=============================================================" -ForegroundColor Yellow
    
    $levelStart = Get-Date
    $results[$LEVEL] = @{}

    # =========================================================================
    # PHASE 1: Core DLNM Analysis
    # =========================================================================
    Write-Host ""
    Write-Host "--- PHASE 1: Core DLNM Analysis ---" -ForegroundColor Magenta
    
    # 1a. Main DLNM with meta-analysis
    Write-Host "  [1/5] 01_dlnm_analysis_v2.R..." -ForegroundColor White
    $t1 = Get-Date
    Set-Location "$BASE_DIR\phase1_r"
    & $RSCRIPT "01_dlnm_analysis_v2.R" $LEVEL 2>&1 | Tee-Object -Variable output
    $results[$LEVEL]["01_dlnm"] = @{
        "duration" = ((Get-Date) - $t1).TotalMinutes
        "success" = ($LASTEXITCODE -eq 0)
    }
    Write-Host "    Done in $([math]::Round($results[$LEVEL]['01_dlnm']['duration'], 1)) min" -ForegroundColor $(if($results[$LEVEL]['01_dlnm']['success']){"Green"}else{"Red"})

    # 1b. Attributable Burden
    Write-Host "  [2/5] 01b_attributable_burden.R..." -ForegroundColor White
    $t1 = Get-Date
    & $RSCRIPT "01b_attributable_burden.R" $LEVEL 2>&1 | Tee-Object -Variable output
    $results[$LEVEL]["01b_burden"] = @{
        "duration" = ((Get-Date) - $t1).TotalMinutes
        "success" = ($LASTEXITCODE -eq 0)
    }
    Write-Host "    Done in $([math]::Round($results[$LEVEL]['01b_burden']['duration'], 1)) min" -ForegroundColor $(if($results[$LEVEL]['01b_burden']['success']){"Green"}else{"Red"})

    # 1c. YLL Calculation
    Write-Host "  [3/5] 01c_yll_calculation.R..." -ForegroundColor White
    $t1 = Get-Date
    & $RSCRIPT "01c_yll_calculation.R" $LEVEL 2>&1 | Tee-Object -Variable output
    $results[$LEVEL]["01c_yll"] = @{
        "duration" = ((Get-Date) - $t1).TotalMinutes
        "success" = ($LASTEXITCODE -eq 0)
    }
    Write-Host "    Done in $([math]::Round($results[$LEVEL]['01c_yll']['duration'], 1)) min" -ForegroundColor $(if($results[$LEVEL]['01c_yll']['success']){"Green"}else{"Red"})

    # 1d. Case-Crossover
    Write-Host "  [4/5] 01d_case_crossover.R..." -ForegroundColor White
    $t1 = Get-Date
    & $RSCRIPT "01d_case_crossover.R" $LEVEL 2>&1 | Tee-Object -Variable output
    $results[$LEVEL]["01d_crossover"] = @{
        "duration" = ((Get-Date) - $t1).TotalMinutes
        "success" = ($LASTEXITCODE -eq 0)
    }
    Write-Host "    Done in $([math]::Round($results[$LEVEL]['01d_crossover']['duration'], 1)) min" -ForegroundColor $(if($results[$LEVEL]['01d_crossover']['success']){"Green"}else{"Red"})

    # 1e. Excess Mortality
    Write-Host "  [5/5] 01e_excess_mortality.R..." -ForegroundColor White
    $t1 = Get-Date
    & $RSCRIPT "01e_excess_mortality.R" $LEVEL 2>&1 | Tee-Object -Variable output
    $results[$LEVEL]["01e_excess"] = @{
        "duration" = ((Get-Date) - $t1).TotalMinutes
        "success" = ($LASTEXITCODE -eq 0)
    }
    Write-Host "    Done in $([math]::Round($results[$LEVEL]['01e_excess']['duration'], 1)) min" -ForegroundColor $(if($results[$LEVEL]['01e_excess']['success']){"Green"}else{"Red"})

    # =========================================================================
    # PHASE 2: Robustness Checks
    # =========================================================================
    Write-Host ""
    Write-Host "--- PHASE 2: Robustness Checks ---" -ForegroundColor Magenta
    Set-Location "$BASE_DIR\phase2_r"

    # 2a. Sensitivity Analyses
    Write-Host "  [1/3] 02a_sensitivity.R..." -ForegroundColor White
    $t1 = Get-Date
    & $RSCRIPT "02a_sensitivity.R" $LEVEL 2>&1 | Tee-Object -Variable output
    $results[$LEVEL]["02a_sensitivity"] = @{
        "duration" = ((Get-Date) - $t1).TotalMinutes
        "success" = ($LASTEXITCODE -eq 0)
    }
    Write-Host "    Done in $([math]::Round($results[$LEVEL]['02a_sensitivity']['duration'], 1)) min" -ForegroundColor $(if($results[$LEVEL]['02a_sensitivity']['success']){"Green"}else{"Red"})

    # 2b. Harvesting (Mortality Displacement)
    Write-Host "  [2/3] 02b_harvesting.R..." -ForegroundColor White
    $t1 = Get-Date
    & $RSCRIPT "02b_harvesting.R" $LEVEL 2>&1 | Tee-Object -Variable output
    $results[$LEVEL]["02b_harvesting"] = @{
        "duration" = ((Get-Date) - $t1).TotalMinutes
        "success" = ($LASTEXITCODE -eq 0)
    }
    Write-Host "    Done in $([math]::Round($results[$LEVEL]['02b_harvesting']['duration'], 1)) min" -ForegroundColor $(if($results[$LEVEL]['02b_harvesting']['success']){"Green"}else{"Red"})

    # 2c. Heatwave Analysis
    Write-Host "  [3/3] 02c_heatwave.R..." -ForegroundColor White
    $t1 = Get-Date
    & $RSCRIPT "02c_heatwave.R" $LEVEL 2>&1 | Tee-Object -Variable output
    $results[$LEVEL]["02c_heatwave"] = @{
        "duration" = ((Get-Date) - $t1).TotalMinutes
        "success" = ($LASTEXITCODE -eq 0)
    }
    Write-Host "    Done in $([math]::Round($results[$LEVEL]['02c_heatwave']['duration'], 1)) min" -ForegroundColor $(if($results[$LEVEL]['02c_heatwave']['success']){"Green"}else{"Red"})

    # =========================================================================
    # PHASE 3: Confounding Controls
    # =========================================================================
    Write-Host ""
    Write-Host "--- PHASE 3: Confounding Controls ---" -ForegroundColor Magenta
    Set-Location "$BASE_DIR\phase3_r"

    # 3a. Supplementary Analyses
    Write-Host "  [1/1] 03a_supplementary.R..." -ForegroundColor White
    $t1 = Get-Date
    & $RSCRIPT "03a_supplementary.R" $LEVEL 2>&1 | Tee-Object -Variable output
    $results[$LEVEL]["03a_supplementary"] = @{
        "duration" = ((Get-Date) - $t1).TotalMinutes
        "success" = ($LASTEXITCODE -eq 0)
    }
    Write-Host "    Done in $([math]::Round($results[$LEVEL]['03a_supplementary']['duration'], 1)) min" -ForegroundColor $(if($results[$LEVEL]['03a_supplementary']['success']){"Green"}else{"Red"})

    # =========================================================================
    # PHASE 4: Heterogeneity & Stratification
    # =========================================================================
    Write-Host ""
    Write-Host "--- PHASE 4: Heterogeneity & Stratification ---" -ForegroundColor Magenta
    Set-Location "$BASE_DIR\phase4_r"

    # 4a. Meta-Regression
    Write-Host "  [1/4] 04a_meta_regression.R..." -ForegroundColor White
    $t1 = Get-Date
    & $RSCRIPT "04a_meta_regression.R" $LEVEL 2>&1 | Tee-Object -Variable output
    $results[$LEVEL]["04a_metareg"] = @{
        "duration" = ((Get-Date) - $t1).TotalMinutes
        "success" = ($LASTEXITCODE -eq 0)
    }
    Write-Host "    Done in $([math]::Round($results[$LEVEL]['04a_metareg']['duration'], 1)) min" -ForegroundColor $(if($results[$LEVEL]['04a_metareg']['success']){"Green"}else{"Red"})

    # 4b. Age Stratification
    Write-Host "  [2/4] 04b_age_stratification.R..." -ForegroundColor White
    $t1 = Get-Date
    & $RSCRIPT "04b_age_stratification.R" $LEVEL 2>&1 | Tee-Object -Variable output
    $results[$LEVEL]["04b_age"] = @{
        "duration" = ((Get-Date) - $t1).TotalMinutes
        "success" = ($LASTEXITCODE -eq 0)
    }
    Write-Host "    Done in $([math]::Round($results[$LEVEL]['04b_age']['duration'], 1)) min" -ForegroundColor $(if($results[$LEVEL]['04b_age']['success']){"Green"}else{"Red"})

    # 4c. Sex Stratification
    Write-Host "  [3/4] 04c_sex_stratification.R..." -ForegroundColor White
    $t1 = Get-Date
    & $RSCRIPT "04c_sex_stratification.R" $LEVEL 2>&1 | Tee-Object -Variable output
    $results[$LEVEL]["04c_sex"] = @{
        "duration" = ((Get-Date) - $t1).TotalMinutes
        "success" = ($LASTEXITCODE -eq 0)
    }
    Write-Host "    Done in $([math]::Round($results[$LEVEL]['04c_sex']['duration'], 1)) min" -ForegroundColor $(if($results[$LEVEL]['04c_sex']['success']){"Green"}else{"Red"})

    # 4d. Cause Stratification
    Write-Host "  [4/4] 04d_cause_stratification.R..." -ForegroundColor White
    $t1 = Get-Date
    & $RSCRIPT "04d_cause_stratification.R" $LEVEL 2>&1 | Tee-Object -Variable output
    $results[$LEVEL]["04d_cause"] = @{
        "duration" = ((Get-Date) - $t1).TotalMinutes
        "success" = ($LASTEXITCODE -eq 0)
    }
    Write-Host "    Done in $([math]::Round($results[$LEVEL]['04d_cause']['duration'], 1)) min" -ForegroundColor $(if($results[$LEVEL]['04d_cause']['success']){"Green"}else{"Red"})

    # Level summary
    $levelEnd = Get-Date
    $levelDuration = ($levelEnd - $levelStart).TotalMinutes
    $successCount = ($results[$LEVEL].Values | Where-Object { $_.success }).Count
    $totalCount = $results[$LEVEL].Count
    
    Write-Host ""
    Write-Host "--- $LEVEL COMPLETE ---" -ForegroundColor Yellow
    Write-Host "  Total time: $([math]::Round($levelDuration, 1)) minutes" -ForegroundColor Gray
    Write-Host "  Success: $successCount / $totalCount scripts" -ForegroundColor $(if($successCount -eq $totalCount){"Green"}else{"Yellow"})
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================
$endTime = Get-Date
$totalDuration = ($endTime - $startTime).TotalMinutes

Write-Host ""
Write-Host "=============================================================" -ForegroundColor Cyan
Write-Host "PIPELINE COMPLETE" -ForegroundColor Cyan
Write-Host "=============================================================" -ForegroundColor Cyan
Write-Host "Total runtime: $([math]::Round($totalDuration, 1)) minutes" -ForegroundColor White
Write-Host ""

foreach ($LEVEL in $levels) {
    Write-Host "--- $LEVEL Results ---" -ForegroundColor Yellow
    foreach ($script in $results[$LEVEL].Keys | Sort-Object) {
        $status = if ($results[$LEVEL][$script].success) { "[OK]" } else { "[FAIL]" }
        $color = if ($results[$LEVEL][$script].success) { "Green" } else { "Red" }
        $duration = [math]::Round($results[$LEVEL][$script].duration, 1)
        Write-Host "  $status $script ($duration min)" -ForegroundColor $color
    }
    Write-Host ""
}

Write-Host "Finished: $(Get-Date)" -ForegroundColor Gray
Write-Host "Results saved to phase*/results/ directories" -ForegroundColor Gray

# Return to base directory
Set-Location $BASE_DIR
