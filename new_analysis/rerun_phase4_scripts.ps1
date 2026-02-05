# Phase 4 R Scripts Batch Runner
# Run all Phase 4 heterogeneity/stratification R scripts sequentially

$ErrorActionPreference = "Continue"
$baseDir = "c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis"
$Rscript = "C:\Program Files\R\R-4.4.1\bin\Rscript.exe"

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  Phase 4 R Scripts Batch Runner" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

# Define Phase 4 scripts
$phase4Scripts = @(
    "phase4_r\04a_meta_regression.R",
    "phase4_r\04b_age_stratification.R",
    "phase4_r\04c_sex_stratification.R",
    "phase4_r\04d_cause_stratification.R"
)

$totalScripts = $phase4Scripts.Count
$currentScript = 0
$successCount = 0
$failCount = 0

Write-Host "Total scripts to run: $totalScripts" -ForegroundColor Yellow
Write-Host ""

foreach ($script in $phase4Scripts) {
    $currentScript++
    $scriptPath = Join-Path $baseDir $script
    $scriptName = Split-Path $script -Leaf
    
    Write-Host "----------------------------------------" -ForegroundColor Gray
    Write-Host "[$currentScript/$totalScripts] Running: $scriptName" -ForegroundColor Green
    Write-Host "Path: $scriptPath" -ForegroundColor Gray
    Write-Host "Started: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" -ForegroundColor Gray
    Write-Host ""
    
    if (Test-Path $scriptPath) {
        $startTime = Get-Date
        
        # Run the R script
        & $Rscript $scriptPath
        $exitCode = $LASTEXITCODE
        
        $endTime = Get-Date
        $duration = $endTime - $startTime
        
        if ($exitCode -eq 0) {
            Write-Host ""
            Write-Host "SUCCESS: $scriptName completed in $($duration.ToString('hh\:mm\:ss'))" -ForegroundColor Green
            $successCount++
        } else {
            Write-Host ""
            Write-Host "FAILED: $scriptName (exit code: $exitCode)" -ForegroundColor Red
            $failCount++
        }
    } else {
        Write-Host "ERROR: Script not found: $scriptPath" -ForegroundColor Red
        $failCount++
    }
    
    Write-Host ""
}

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  BATCH COMPLETE" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Successful: $successCount / $totalScripts" -ForegroundColor $(if ($successCount -eq $totalScripts) { "Green" } else { "Yellow" })
Write-Host "Failed:     $failCount / $totalScripts" -ForegroundColor $(if ($failCount -eq 0) { "Green" } else { "Red" })
Write-Host ""
Write-Host "Completed at: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" -ForegroundColor Gray
