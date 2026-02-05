# =============================================================================
# Batch script to run all intermediate-level analyses
# =============================================================================
# Usage: .\run_intermediate_batch.ps1
# =============================================================================

$ErrorActionPreference = "Continue"
$startTime = Get-Date

# Full path to Rscript
$Rscript = "C:\Program Files\R\R-4.4.1\bin\Rscript.exe"

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "RUNNING ALL INTERMEDIATE ANALYSES" -ForegroundColor Cyan
Write-Host "Started: $startTime" -ForegroundColor Cyan
Write-Host "Using R: $Rscript" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

$baseDir = Split-Path -Parent $MyInvocation.MyCommand.Path

# Track results
$results = @()

# -----------------------------------------------------------------------------
# Phase 1: Core Models
# -----------------------------------------------------------------------------
Write-Host "`n[PHASE 1a] DLNM Analysis..." -ForegroundColor Yellow
$t1 = Get-Date
Set-Location "$baseDir\phase1_r"
& $Rscript 01_dlnm_analysis_v2.R intermediate 2>&1 | Tee-Object -Variable dlnmOutput
$dlnmExit = $LASTEXITCODE
$results += [PSCustomObject]@{Script="01_dlnm_analysis_v2.R"; ExitCode=$dlnmExit; Duration=((Get-Date)-$t1).TotalMinutes}
Write-Host "DLNM completed in $([math]::Round(((Get-Date)-$t1).TotalMinutes, 1)) minutes" -ForegroundColor Green

Write-Host "`n[PHASE 1b] Attributable Burden..." -ForegroundColor Yellow
$t2 = Get-Date
& $Rscript 01b_attributable_burden.R intermediate 2>&1 | Tee-Object -Variable burdenOutput
$burdenExit = $LASTEXITCODE
$results += [PSCustomObject]@{Script="01b_attributable_burden.R"; ExitCode=$burdenExit; Duration=((Get-Date)-$t2).TotalMinutes}
Write-Host "Burden completed in $([math]::Round(((Get-Date)-$t2).TotalMinutes, 1)) minutes" -ForegroundColor Green

Write-Host "`n[PHASE 1c] YLL Calculation..." -ForegroundColor Yellow
$t3 = Get-Date
& $Rscript 01c_yll_calculation.R intermediate 2>&1 | Tee-Object -Variable yllOutput
$yllExit = $LASTEXITCODE
$results += [PSCustomObject]@{Script="01c_yll_calculation.R"; ExitCode=$yllExit; Duration=((Get-Date)-$t3).TotalMinutes}
Write-Host "YLL completed in $([math]::Round(((Get-Date)-$t3).TotalMinutes, 1)) minutes" -ForegroundColor Green

# -----------------------------------------------------------------------------
# Phase 2: Robustness
# -----------------------------------------------------------------------------
Write-Host "`n[PHASE 2a] Sensitivity Analyses..." -ForegroundColor Yellow
$t4 = Get-Date
Set-Location "$baseDir\phase2_r"
& $Rscript 02a_sensitivity.R intermediate 2>&1 | Tee-Object -Variable sensOutput
$sensExit = $LASTEXITCODE
$results += [PSCustomObject]@{Script="02a_sensitivity.R"; ExitCode=$sensExit; Duration=((Get-Date)-$t4).TotalMinutes}
Write-Host "Sensitivity completed in $([math]::Round(((Get-Date)-$t4).TotalMinutes, 1)) minutes" -ForegroundColor Green

Write-Host "`n[PHASE 2b] Harvesting Analysis..." -ForegroundColor Yellow
$t5 = Get-Date
& $Rscript 02b_harvesting.R intermediate 2>&1 | Tee-Object -Variable harvOutput
$harvExit = $LASTEXITCODE
$results += [PSCustomObject]@{Script="02b_harvesting.R"; ExitCode=$harvExit; Duration=((Get-Date)-$t5).TotalMinutes}
Write-Host "Harvesting completed in $([math]::Round(((Get-Date)-$t5).TotalMinutes, 1)) minutes" -ForegroundColor Green

Write-Host "`n[PHASE 2c] Heatwave Analysis..." -ForegroundColor Yellow
$t6 = Get-Date
& $Rscript 02c_heatwave.R intermediate 2>&1 | Tee-Object -Variable heatOutput
$heatExit = $LASTEXITCODE
$results += [PSCustomObject]@{Script="02c_heatwave.R"; ExitCode=$heatExit; Duration=((Get-Date)-$t6).TotalMinutes}
Write-Host "Heatwave completed in $([math]::Round(((Get-Date)-$t6).TotalMinutes, 1)) minutes" -ForegroundColor Green

# -----------------------------------------------------------------------------
# Phase 3: Confounding Controls
# -----------------------------------------------------------------------------
Write-Host "`n[PHASE 3a] Supplementary Confounding..." -ForegroundColor Yellow
$t7 = Get-Date
Set-Location "$baseDir\phase3_r"
& $Rscript 03a_supplementary.R intermediate 2>&1 | Tee-Object -Variable suppOutput
$suppExit = $LASTEXITCODE
$results += [PSCustomObject]@{Script="03a_supplementary.R"; ExitCode=$suppExit; Duration=((Get-Date)-$t7).TotalMinutes}
Write-Host "Supplementary completed in $([math]::Round(((Get-Date)-$t7).TotalMinutes, 1)) minutes" -ForegroundColor Green

# -----------------------------------------------------------------------------
# Phase 4: Heterogeneity
# -----------------------------------------------------------------------------
Write-Host "`n[PHASE 4a] Meta-Regression..." -ForegroundColor Yellow
$t8 = Get-Date
Set-Location "$baseDir\phase4_r"
& $Rscript 04a_meta_regression.R intermediate 2>&1 | Tee-Object -Variable metaOutput
$metaExit = $LASTEXITCODE
$results += [PSCustomObject]@{Script="04a_meta_regression.R"; ExitCode=$metaExit; Duration=((Get-Date)-$t8).TotalMinutes}
Write-Host "Meta-regression completed in $([math]::Round(((Get-Date)-$t8).TotalMinutes, 1)) minutes" -ForegroundColor Green

Write-Host "`n[PHASE 4b] Age Stratification..." -ForegroundColor Yellow
$t9 = Get-Date
& $Rscript 04b_age_stratification.R intermediate 2>&1 | Tee-Object -Variable ageOutput
$ageExit = $LASTEXITCODE
$results += [PSCustomObject]@{Script="04b_age_stratification.R"; ExitCode=$ageExit; Duration=((Get-Date)-$t9).TotalMinutes}
Write-Host "Age stratification completed in $([math]::Round(((Get-Date)-$t9).TotalMinutes, 1)) minutes" -ForegroundColor Green

Write-Host "`n[PHASE 4c] Sex Stratification..." -ForegroundColor Yellow
$t10 = Get-Date
& $Rscript 04c_sex_stratification.R intermediate 2>&1 | Tee-Object -Variable sexOutput
$sexExit = $LASTEXITCODE
$results += [PSCustomObject]@{Script="04c_sex_stratification.R"; ExitCode=$sexExit; Duration=((Get-Date)-$t10).TotalMinutes}
Write-Host "Sex stratification completed in $([math]::Round(((Get-Date)-$t10).TotalMinutes, 1)) minutes" -ForegroundColor Green

Write-Host "`n[PHASE 4d] Cause Stratification..." -ForegroundColor Yellow
$t11 = Get-Date
& $Rscript 04d_cause_stratification.R intermediate 2>&1 | Tee-Object -Variable causeOutput
$causeExit = $LASTEXITCODE
$results += [PSCustomObject]@{Script="04d_cause_stratification.R"; ExitCode=$causeExit; Duration=((Get-Date)-$t11).TotalMinutes}
Write-Host "Cause stratification completed in $([math]::Round(((Get-Date)-$t11).TotalMinutes, 1)) minutes" -ForegroundColor Green

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
$endTime = Get-Date
$totalDuration = ($endTime - $startTime).TotalMinutes

Write-Host "`n============================================================" -ForegroundColor Cyan
Write-Host "BATCH COMPLETE" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "Total duration: $([math]::Round($totalDuration, 1)) minutes"
Write-Host "`nResults summary:"
$results | Format-Table -AutoSize

$failed = $results | Where-Object { $_.ExitCode -ne 0 }
if ($failed.Count -gt 0) {
    Write-Host "`nFAILED SCRIPTS:" -ForegroundColor Red
    $failed | Format-Table -AutoSize
} else {
    Write-Host "`nAll scripts completed successfully!" -ForegroundColor Green
}

Set-Location $baseDir
