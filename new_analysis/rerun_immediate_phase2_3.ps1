# Rerun Phase 2 & 3 R scripts for IMMEDIATE level
# Sensitivity runs last (longest running)

$RScript = "C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
$BaseDir = Split-Path -Parent $MyInvocation.MyCommand.Path

Write-Host "`n" + ("=" * 80) -ForegroundColor Green
Write-Host "PHASE 2 & 3: Immediate Level Rerun" -ForegroundColor Green
Write-Host ("=" * 80) + "`n" -ForegroundColor Green

# Phase 2: Harvesting
Write-Host "[1/4] Running 02b_harvesting.R immediate..." -ForegroundColor Yellow
& $RScript "$BaseDir\phase2_r\02b_harvesting.R" immediate
if ($LASTEXITCODE -ne 0) { Write-Host "  FAILED!" -ForegroundColor Red } else { Write-Host "  Done!" -ForegroundColor Green }

# Phase 2: Heatwave
Write-Host "`n[2/4] Running 02c_heatwave.R immediate..." -ForegroundColor Yellow
& $RScript "$BaseDir\phase2_r\02c_heatwave.R" immediate
if ($LASTEXITCODE -ne 0) { Write-Host "  FAILED!" -ForegroundColor Red } else { Write-Host "  Done!" -ForegroundColor Green }

# Phase 3: Supplementary
Write-Host "`n[3/4] Running 03a_supplementary.R immediate..." -ForegroundColor Yellow
& $RScript "$BaseDir\phase3_r\03a_supplementary.R" immediate
if ($LASTEXITCODE -ne 0) { Write-Host "  FAILED!" -ForegroundColor Red } else { Write-Host "  Done!" -ForegroundColor Green }

# Phase 2: Sensitivity (LAST - longest running)
Write-Host "`n[4/4] Running 02a_sensitivity.R immediate..." -ForegroundColor Yellow
& $RScript "$BaseDir\phase2_r\02a_sensitivity.R" immediate
if ($LASTEXITCODE -ne 0) { Write-Host "  FAILED!" -ForegroundColor Red } else { Write-Host "  Done!" -ForegroundColor Green }

Write-Host "`n" + ("=" * 80) -ForegroundColor Green
Write-Host "ALL DONE!" -ForegroundColor Green
Write-Host ("=" * 80) + "`n" -ForegroundColor Green
