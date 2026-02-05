# =============================================================================
# BATCH RE-RUN: All Immediate-Level R Scripts
# =============================================================================
# Run this script to re-run all immediate-level R scripts with corrected data
# The intermediate-level results from today (22:00+) are already valid
# =============================================================================

Write-Host "=" * 80 -ForegroundColor Cyan
Write-Host "BATCH RE-RUN: Immediate-Level R Scripts" -ForegroundColor Cyan
Write-Host "Started: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" -ForegroundColor Cyan
Write-Host "=" * 80 -ForegroundColor Cyan

$R_EXE = "C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
$BASE_DIR = "c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis"

# Track timing
$start_time = Get-Date

# =============================================================================
# PHASE 1: Core DLNM
# =============================================================================
Write-Host "`n" + "=" * 80 -ForegroundColor Green
Write-Host "PHASE 1: Core DLNM (Immediate)" -ForegroundColor Green
Write-Host "=" * 80 -ForegroundColor Green

Set-Location "$BASE_DIR\phase1_r"

Write-Host "`n[1/5] Running 01_dlnm_analysis_v2.R immediate..." -ForegroundColor Yellow
& $R_EXE 01_dlnm_analysis_v2.R immediate
if ($LASTEXITCODE -ne 0) { Write-Host "ERROR in DLNM!" -ForegroundColor Red }

Write-Host "`n[2/5] Running 01b_attributable_burden.R immediate..." -ForegroundColor Yellow
& $R_EXE 01b_attributable_burden.R immediate
if ($LASTEXITCODE -ne 0) { Write-Host "ERROR in Attributable Burden!" -ForegroundColor Red }

Write-Host "`n[3/5] Running 01c_yll_calculation.R immediate..." -ForegroundColor Yellow
& $R_EXE 01c_yll_calculation.R immediate
if ($LASTEXITCODE -ne 0) { Write-Host "ERROR in YLL!" -ForegroundColor Red }

Write-Host "`n[4/5] Running 01d_case_crossover.R immediate..." -ForegroundColor Yellow
& $R_EXE 01d_case_crossover.R immediate
if ($LASTEXITCODE -ne 0) { Write-Host "ERROR in Case-Crossover!" -ForegroundColor Red }

Write-Host "`n[5/5] Running 01e_excess_mortality.R immediate..." -ForegroundColor Yellow
& $R_EXE 01e_excess_mortality.R immediate
if ($LASTEXITCODE -ne 0) { Write-Host "ERROR in Excess Mortality!" -ForegroundColor Red }

# =============================================================================
# PHASE 2: Robustness
# =============================================================================
Write-Host "`n" + "=" * 80 -ForegroundColor Green
Write-Host "PHASE 2: Robustness (Immediate)" -ForegroundColor Green
Write-Host "=" * 80 -ForegroundColor Green

Set-Location "$BASE_DIR\phase2_r"

Write-Host "`n[1/3] Running 02a_sensitivity.R immediate..." -ForegroundColor Yellow
& $R_EXE 02a_sensitivity.R immediate
if ($LASTEXITCODE -ne 0) { Write-Host "ERROR in Sensitivity!" -ForegroundColor Red }

Write-Host "`n[2/3] Running 02b_harvesting.R immediate..." -ForegroundColor Yellow
& $R_EXE 02b_harvesting.R immediate
if ($LASTEXITCODE -ne 0) { Write-Host "ERROR in Harvesting!" -ForegroundColor Red }

Write-Host "`n[3/3] Running 02c_heatwave.R immediate..." -ForegroundColor Yellow
& $R_EXE 02c_heatwave.R immediate
if ($LASTEXITCODE -ne 0) { Write-Host "ERROR in Heatwave!" -ForegroundColor Red }

# =============================================================================
# PHASE 3: Confounding
# =============================================================================
Write-Host "`n" + "=" * 80 -ForegroundColor Green
Write-Host "PHASE 3: Confounding (Immediate)" -ForegroundColor Green
Write-Host "=" * 80 -ForegroundColor Green

Set-Location "$BASE_DIR\phase3_r"

Write-Host "`n[1/1] Running 03a_supplementary.R immediate..." -ForegroundColor Yellow
& $R_EXE 03a_supplementary.R immediate
if ($LASTEXITCODE -ne 0) { Write-Host "ERROR in Supplementary!" -ForegroundColor Red }

# =============================================================================
# SUMMARY
# =============================================================================
$end_time = Get-Date
$duration = $end_time - $start_time

Write-Host "`n" + "=" * 80 -ForegroundColor Cyan
Write-Host "BATCH RE-RUN COMPLETE" -ForegroundColor Cyan
Write-Host "=" * 80 -ForegroundColor Cyan
Write-Host "Started:  $($start_time.ToString('yyyy-MM-dd HH:mm:ss'))"
Write-Host "Finished: $($end_time.ToString('yyyy-MM-dd HH:mm:ss'))"
Write-Host "Duration: $($duration.ToString('hh\:mm\:ss'))"
Write-Host "=" * 80 -ForegroundColor Cyan

# Return to base directory
Set-Location $BASE_DIR
