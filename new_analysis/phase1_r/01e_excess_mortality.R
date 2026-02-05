# =============================================================================
# 01e_excess_mortality.R
# Excess Mortality Validation using Counterfactual Approach
# =============================================================================
#
# Method:
# 1. Fit baseline mortality model (GAM with seasonality + trend)
# 2. Calculate expected deaths under no temperature effect
# 3. Compare observed vs expected during heat/cold periods
# 4. Validate DLNM-derived attributable burden
#
# References:
# - Fouillet et al. (2006) - Excess mortality approach
# - Gasparrini et al. (2022) - Comparison of burden estimation methods
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(jsonlite)
  library(mgcv)         # For GAM
  library(splines)
})

# Get exposure type from command line args
args <- commandArgs(trailingOnly = TRUE)
EXPOSURE_TYPE <- if (length(args) > 0) args[1] else "intermediate"

cat("=======================================================\n")
cat("EXCESS MORTALITY VALIDATION:", EXPOSURE_TYPE, "\n")
cat("=======================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------
# Get script directory for proper path resolution
SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("--file=", "", file_arg)))
  } else {
    getwd()
  }
})

DATA_DIR <- file.path(SCRIPT_DIR, "..", "phase0_data_prep", "results")
DLNM_DIR <- file.path(SCRIPT_DIR, "results")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Script directory:", SCRIPT_DIR, "\n")
cat("Data directory:", normalizePath(DATA_DIR, mustWork = FALSE), "\n")

cat("[1] Loading data...\n")

# Load mortality data
mort_file <- file.path(DATA_DIR, paste0("mortality_", EXPOSURE_TYPE, "_daily_elderly.parquet"))
mort <- as.data.table(read_parquet(mort_file))

# Load ERA5 temperature data
era5_file <- file.path(DATA_DIR, paste0("era5_", EXPOSURE_TYPE, "_daily.parquet"))
era5 <- as.data.table(read_parquet(era5_file))

# Standardize column names - handle both old (intermediate_code) and new (region_code) formats
region_col <- paste0(EXPOSURE_TYPE, "_code")
if (region_col %in% names(mort)) {
  setnames(mort, region_col, "region_code")
} else if (!"region_code" %in% names(mort)) {
  stop("Neither '", region_col, "' nor 'region_code' found in mortality data")
}

# Convert dates for merge compatibility (mort has POSIXct, era5 has Date)
mort[, date := as.character(as.Date(date))]
era5[, date := as.character(date)]

# Merge
setkey(mort, region_code, date)
setkey(era5, region_code, date)
data <- merge(mort, era5, by = c("region_code", "date"))

# Convert date back to Date for analysis
data[, date := as.Date(date)]

# Working columns
data[, deaths := deaths_elderly]
data[, tmean := temp_mean]
data <- data[!is.na(deaths) & !is.na(tmean)]

# Time variables
data[, `:=`(
  year = year(date),
  month = month(date),
  doy = yday(date),
  dow = wday(date),
  time = as.numeric(date - min(date))
)]

cat("  Rows:", format(nrow(data), big.mark = ","), "\n")
cat("  Regions:", uniqueN(data$region_code), "\n")
cat("  Years:", min(data$year), "to", max(data$year), "\n")

# Load DLNM results for MMT
dlnm_file <- file.path(DLNM_DIR, paste0("dlnm_r_", EXPOSURE_TYPE, "_results_v2.json"))
dlnm_results <- fromJSON(dlnm_file)
pooled_mmt <- dlnm_results$pooled$mmt
cat("  Pooled MMT:", round(pooled_mmt, 2), "°C\n")

# -----------------------------------------------------------------------------
# 2. Aggregate to National Level
# -----------------------------------------------------------------------------
cat("\n[2] Aggregating to national level...\n")

national <- data[, .(
  deaths = sum(deaths),
  tmean = mean(tmean),
  tmax = max(tmean),
  tmin = min(tmean)
), by = .(date, year, month, doy, dow)]

national[, time := as.numeric(date - min(date))]
setorder(national, date)

cat("  National days:", nrow(national), "\n")
cat("  Total deaths:", format(sum(national$deaths), big.mark = ","), "\n")
cat("  Mean daily deaths:", round(mean(national$deaths), 1), "\n")

# Define temperature thresholds
p1 <- quantile(national$tmean, 0.01)
p2_5 <- quantile(national$tmean, 0.025)
p97_5 <- quantile(national$tmean, 0.975)
p99 <- quantile(national$tmean, 0.99)

cat("\n  Temperature percentiles:\n")
cat("    P1:", round(p1, 2), "°C\n")
cat("    P2.5:", round(p2_5, 2), "°C\n")
cat("    P97.5:", round(p97_5, 2), "°C\n")
cat("    P99:", round(p99, 2), "°C\n")
cat("    MMT:", round(pooled_mmt, 2), "°C\n")

# -----------------------------------------------------------------------------
# 3. Fit Baseline Mortality Model (GAM)
# -----------------------------------------------------------------------------
cat("\n[3] Fitting baseline mortality model...\n")

# GAM with:
# - Cyclic cubic spline for seasonality (day of year)
# - Thin plate spline for long-term trend
# - Day-of-week factor
# - No temperature term (counterfactual)

n_years <- uniqueN(national$year)
trend_k <- max(10, n_years * 2)  # Flexible trend

cat("  Fitting GAM (seasonality + trend)...\n")

# Fit quasi-Poisson GAM for overdispersion
gam_fit <- tryCatch({
  gam(
    deaths ~ s(doy, bs = "cc", k = 10) +  # Cyclic seasonality
             s(time, bs = "tp", k = trend_k) +  # Long-term trend
             factor(dow),  # Day of week
    data = national,
    family = quasipoisson(link = "log"),
    method = "REML"
  )
}, error = function(e) {
  cat("  GAM failed:", conditionMessage(e), "\n")
  cat("  Trying simpler model...\n")
  
  # Simpler model with fewer knots
  gam(
    deaths ~ s(doy, bs = "cc", k = 6) +
             s(time, bs = "tp", k = n_years + 2) +
             factor(dow),
    data = national,
    family = quasipoisson(link = "log"),
    method = "REML"
  )
})

cat("  GAM deviance explained:", round(summary(gam_fit)$dev.expl * 100, 1), "%\n")
cat("  Scale parameter:", round(summary(gam_fit)$scale, 2), "\n")

# Get expected deaths (baseline without temperature)
national[, expected := predict(gam_fit, type = "response")]

# Residuals
national[, excess := deaths - expected]
national[, excess_pct := (deaths - expected) / expected * 100]

# -----------------------------------------------------------------------------
# 4. Classify Heat/Cold Periods
# -----------------------------------------------------------------------------
cat("\n[4] Classifying temperature periods...\n")

national[, `:=`(
  is_heat = tmean > pooled_mmt,
  is_cold = tmean < pooled_mmt,
  is_extreme_heat = tmean > p97_5,
  is_extreme_cold = tmean < p2_5,
  is_heat_99 = tmean > p99,
  is_cold_1 = tmean < p1,
  is_optimal = abs(tmean - pooled_mmt) < 2  # Within 2°C of MMT
)]

cat("  Heat days (> MMT):", sum(national$is_heat), "\n")
cat("  Cold days (< MMT):", sum(national$is_cold), "\n")
cat("  Extreme heat (> P97.5):", sum(national$is_extreme_heat), "\n")
cat("  Extreme cold (< P2.5):", sum(national$is_extreme_cold), "\n")

# -----------------------------------------------------------------------------
# 5. Calculate Excess Mortality by Period
# -----------------------------------------------------------------------------
cat("\n[5] Calculating excess mortality...\n")

# Function to summarize excess
calc_excess_summary <- function(dt, label) {
  list(
    label = label,
    n_days = nrow(dt),
    total_deaths = sum(dt$deaths),
    total_expected = sum(dt$expected),
    total_excess = sum(dt$excess),
    excess_pct = sum(dt$excess) / sum(dt$expected) * 100,
    mean_daily_excess = mean(dt$excess),
    se_daily_excess = sd(dt$excess) / sqrt(nrow(dt))
  )
}

# Calculate for each period
excess_all <- calc_excess_summary(national, "All days")
excess_heat <- calc_excess_summary(national[is_heat == TRUE], "Heat (> MMT)")
excess_cold <- calc_excess_summary(national[is_cold == TRUE], "Cold (< MMT)")
excess_extreme_heat <- calc_excess_summary(national[is_extreme_heat == TRUE], "Extreme heat (> P97.5)")
excess_extreme_cold <- calc_excess_summary(national[is_extreme_cold == TRUE], "Extreme cold (< P2.5)")
excess_heat_99 <- calc_excess_summary(national[is_heat_99 == TRUE], "Extreme heat (> P99)")
excess_cold_1 <- calc_excess_summary(national[is_cold_1 == TRUE], "Extreme cold (< P1)")
excess_optimal <- calc_excess_summary(national[is_optimal == TRUE], "Near MMT (±2°C)")

# -----------------------------------------------------------------------------
# 6. Annual Excess Analysis
# -----------------------------------------------------------------------------
cat("\n[6] Annual excess analysis...\n")

annual_excess <- national[, .(
  total_deaths = sum(deaths),
  total_expected = sum(expected),
  total_excess = sum(excess),
  heat_excess = sum(excess[is_heat == TRUE]),
  cold_excess = sum(excess[is_cold == TRUE]),
  extreme_heat_excess = sum(excess[is_extreme_heat == TRUE]),
  extreme_cold_excess = sum(excess[is_extreme_cold == TRUE])
), by = year]

annual_excess[, `:=`(
  excess_pct = total_excess / total_expected * 100,
  heat_excess_pct = heat_excess / total_expected * 100,
  cold_excess_pct = cold_excess / total_expected * 100
)]

# -----------------------------------------------------------------------------
# 7. Compare with DLNM Attributable Burden
# -----------------------------------------------------------------------------
cat("\n[7] Comparing with DLNM attributable burden...\n")

# Load attributable burden
burden_file <- file.path(DLNM_DIR, paste0("attributable_burden_r_", EXPOSURE_TYPE, ".json"))
if (file.exists(burden_file)) {
  burden <- fromJSON(burden_file)
  dlnm_heat_an <- burden$national$total_heat_an
  dlnm_cold_an <- burden$national$total_cold_an
  dlnm_total_an <- burden$national$total_an
  
  # Scale DLNM estimates to match national aggregation
  # (DLNM is region-level, excess is national)
  excess_heat_total <- excess_heat$total_excess
  excess_cold_total <- excess_cold$total_excess
  
  comparison <- list(
    method = "Comparison: Excess Mortality vs DLNM AF",
    note = "Excess mortality is from GAM baseline (no temp), DLNM uses pooled E-R curve",
    
    heat = list(
      dlnm_attributable = dlnm_heat_an,
      excess_mortality = excess_heat_total,
      ratio = excess_heat_total / dlnm_heat_an
    ),
    
    cold = list(
      dlnm_attributable = dlnm_cold_an,
      excess_mortality = excess_cold_total,
      ratio = excess_cold_total / dlnm_cold_an
    ),
    
    total = list(
      dlnm_attributable = dlnm_total_an,
      excess_mortality = excess_heat_total + excess_cold_total,
      ratio = (excess_heat_total + excess_cold_total) / dlnm_total_an
    )
  )
} else {
  comparison <- NULL
  cat("  DLNM burden file not found for comparison\n")
}

# -----------------------------------------------------------------------------
# 8. Print Results
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("EXCESS MORTALITY RESULTS (", EXPOSURE_TYPE, ")\n", sep = "")
cat("=======================================================\n")

cat("\n--- BASELINE MODEL ---\n")
cat("Model: GAM quasi-Poisson (seasonality + trend + DOW)\n")
cat("Deviance explained:", round(summary(gam_fit)$dev.expl * 100, 1), "%\n")

cat("\n--- EXCESS BY TEMPERATURE PERIOD ---\n")
print_excess <- function(x) {
  cat(sprintf("%-25s: %8s observed, %8s expected, %+8s excess (%+.1f%%)\n",
              x$label,
              format(round(x$total_deaths), big.mark = ","),
              format(round(x$total_expected), big.mark = ","),
              format(round(x$total_excess), big.mark = ","),
              x$excess_pct))
}

print_excess(excess_all)
print_excess(excess_heat)
print_excess(excess_cold)
print_excess(excess_extreme_heat)
print_excess(excess_extreme_cold)
print_excess(excess_optimal)

if (!is.null(comparison)) {
  cat("\n--- COMPARISON WITH DLNM ---\n")
  cat(sprintf("Heat: DLNM = %s, Excess = %s (ratio: %.2f)\n",
              format(round(comparison$heat$dlnm_attributable), big.mark = ","),
              format(round(comparison$heat$excess_mortality), big.mark = ","),
              comparison$heat$ratio))
  cat(sprintf("Cold: DLNM = %s, Excess = %s (ratio: %.2f)\n",
              format(round(comparison$cold$dlnm_attributable), big.mark = ","),
              format(round(comparison$cold$excess_mortality), big.mark = ","),
              comparison$cold$ratio))
  cat(sprintf("Total: DLNM = %s, Excess = %s (ratio: %.2f)\n",
              format(round(comparison$total$dlnm_attributable), big.mark = ","),
              format(round(comparison$total$excess_mortality), big.mark = ","),
              comparison$total$ratio))
}

cat("\n--- ANNUAL TREND ---\n")
print(annual_excess[order(year), .(year, total_deaths, total_excess = round(total_excess),
                                   excess_pct = round(excess_pct, 2))])

# -----------------------------------------------------------------------------
# 9. Save Results
# -----------------------------------------------------------------------------
cat("\n[9] Saving results...\n")

output <- list(
  exposure_type = EXPOSURE_TYPE,
  analysis_date = as.character(Sys.time()),
  
  model = list(
    type = "GAM quasi-Poisson",
    components = c("cyclic seasonality", "long-term trend", "day of week"),
    deviance_explained = summary(gam_fit)$dev.expl,
    scale_parameter = summary(gam_fit)$scale
  ),
  
  thresholds = list(
    mmt = pooled_mmt,
    p1 = p1,
    p2_5 = p2_5,
    p97_5 = p97_5,
    p99 = p99
  ),
  
  excess_summary = list(
    all = excess_all,
    heat = excess_heat,
    cold = excess_cold,
    extreme_heat = excess_extreme_heat,
    extreme_cold = excess_extreme_cold,
    near_optimal = excess_optimal
  ),
  
  comparison_with_dlnm = comparison,
  
  annual = as.data.frame(annual_excess)
)

# Save JSON
output_file <- file.path(OUTPUT_DIR, paste0("excess_mortality_r_", EXPOSURE_TYPE, ".json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 4)
cat("  JSON saved to:", output_file, "\n")

# Save daily data for plotting
daily_file <- file.path(OUTPUT_DIR, paste0("excess_mortality_r_", EXPOSURE_TYPE, "_daily.csv"))
fwrite(national[, .(date, deaths, expected, excess, excess_pct, tmean, 
                    is_heat, is_cold, is_extreme_heat, is_extreme_cold)],
       daily_file)
cat("  Daily CSV saved to:", daily_file, "\n")

# Save annual
annual_file <- file.path(OUTPUT_DIR, paste0("excess_mortality_r_", EXPOSURE_TYPE, "_annual.csv"))
fwrite(annual_excess, annual_file)
cat("  Annual CSV saved to:", annual_file, "\n")

cat("\n=======================================================\n")
cat("Done!", as.character(Sys.time()), "\n")
cat("=======================================================\n")
