# Run All R Analysis Scripts - INTERMEDIATE Level
# Temperature-Mortality Analysis Pipeline
# Created: December 18, 2025

$ErrorActionPreference = "Continue"
$RScript = "C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
$BaseDir = Split-Path -Parent $MyInvocation.MyCommand.Path

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "  RUNNING ALL R SCRIPTS - INTERMEDIATE LEVEL (133 regions)" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host ""

$startTime = Get-Date
$results = @()

# Define all scripts to run
$scripts = @(
    @{Phase="Phase 1"; Script="phase1_r/01_dlnm_analysis_v2.R"; Args="intermediate"},
    @{Phase="Phase 1"; Script="phase1_r/01b_attributable_burden.R"; Args="intermediate"},
    @{Phase="Phase 1"; Script="phase1_r/01c_yll_calculation.R"; Args="intermediate"},
    @{Phase="Phase 1"; Script="phase1_r/01d_case_crossover.R"; Args="intermediate"},
    @{Phase="Phase 1"; Script="phase1_r/01e_excess_mortality.R"; Args="intermediate"},
    @{Phase="Phase 2"; Script="phase2_r/02a_sensitivity.R"; Args="intermediate"},
    @{Phase="Phase 2"; Script="phase2_r/02b_harvesting.R"; Args="intermediate"},
    @{Phase="Phase 2"; Script="phase2_r/02c_heatwave.R"; Args="intermediate"},
    @{Phase="Phase 3"; Script="phase3_r/03a_supplementary.R"; Args="intermediate"},
    @{Phase="Phase 4"; Script="phase4_r/04a_meta_regression.R"; Args="intermediate"},
    @{Phase="Phase 4"; Script="phase4_r/04b_age_stratification.R"; Args="intermediate"},
    @{Phase="Phase 4"; Script="phase4_r/04c_sex_stratification.R"; Args="intermediate"},
    @{Phase="Phase 4"; Script="phase4_r/04d_cause_stratification.R"; Args="intermediate"}
)

$totalScripts = $scripts.Count
$currentScript = 0

foreach ($s in $scripts) {
    $currentScript++
    $scriptPath = Join-Path $BaseDir $s.Script
    
    Write-Host ""
    Write-Host "------------------------------------------------------------" -ForegroundColor Yellow
    Write-Host "[$currentScript/$totalScripts] $($s.Phase): $($s.Script)" -ForegroundColor Yellow
    Write-Host "------------------------------------------------------------" -ForegroundColor Yellow
    
    $scriptStart = Get-Date
    
    try {
        & $RScript $scriptPath $s.Args 2>&1 | ForEach-Object { Write-Host $_ }
        $exitCode = $LASTEXITCODE
        $status = if ($exitCode -eq 0) { "SUCCESS" } else { "FAILED (exit $exitCode)" }
    } catch {
        $status = "ERROR: $_"
        $exitCode = 1
    }
    
    $scriptEnd = Get-Date
    $duration = $scriptEnd - $scriptStart
    
    $results += [PSCustomObject]@{
        Phase = $s.Phase
        Script = $s.Script
        Status = $status
        Duration = $duration.ToString("mm\:ss")
    }
    
    $color = if ($exitCode -eq 0) { "Green" } else { "Red" }
    Write-Host ""
    Write-Host "  Status: $status | Duration: $($duration.ToString('mm\:ss'))" -ForegroundColor $color
}

$endTime = Get-Date
$totalDuration = $endTime - $startTime

Write-Host ""
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "  PIPELINE COMPLETE - INTERMEDIATE LEVEL" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Summary:" -ForegroundColor White
Write-Host ""

$results | Format-Table -AutoSize

$successCount = ($results | Where-Object { $_.Status -eq "SUCCESS" }).Count
$failCount = $totalScripts - $successCount

Write-Host ""
Write-Host "Total: $successCount/$totalScripts succeeded" -ForegroundColor $(if ($failCount -eq 0) { "Green" } else { "Yellow" })
Write-Host "Total Duration: $($totalDuration.ToString('hh\:mm\:ss'))" -ForegroundColor Cyan
Write-Host ""
