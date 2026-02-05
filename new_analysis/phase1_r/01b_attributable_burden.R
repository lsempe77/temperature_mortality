# =============================================================================
# 01b_attributable_burden.R
# Attributable Burden Calculation using R DLNM Results
# =============================================================================
# 
# Method (Gasparrini & Leone, 2014):
# 1. Load pooled DLNM coefficients from 01_dlnm_analysis_v2.R
# 2. For each region/day, compute RR relative to MMT
# 3. Calculate AF = (RR - 1) / RR
# 4. Sum attributable deaths for heat (temp > MMT) and cold (temp < MMT)
# 5. Aggregate to national level
#
# =============================================================================

suppressPackageStartupMessages({
  library(dlnm)
  library(data.table)
  library(arrow)
  library(jsonlite)
  library(splines)
})

# Get exposure type from command line args
args <- commandArgs(trailingOnly = TRUE)
EXPOSURE_TYPE <- if (length(args) > 0) args[1] else "intermediate"

cat("=======================================================\n")
cat("ATTRIBUTABLE BURDEN CALCULATION:", EXPOSURE_TYPE, "\n")
cat("=======================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# 1. Load DLNM Results
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

cat("[1] Loading DLNM results...\n")

dlnm_file <- file.path(DLNM_DIR, paste0("dlnm_r_", EXPOSURE_TYPE, "_results_v2.json"))
if (!file.exists(dlnm_file)) {
  stop("DLNM results not found: ", dlnm_file, "\nRun 01_dlnm_analysis_v2.R first!")
}

dlnm_results <- fromJSON(dlnm_file)
cat("  Loaded:", dlnm_file, "\n")

# Extract parameters
TEMP_BOUNDARY <- dlnm_results$dlnm_params$temp_boundary
TEMP_KNOTS <- dlnm_results$dlnm_params$temp_knots
MAX_LAG <- dlnm_results$dlnm_params$max_lag
LAG_DF <- dlnm_results$dlnm_params$lag_df

# Get pooled coefficients and vcov
if (is.null(dlnm_results$pooled)) {
  stop("No pooled results available!")
}

pooled_coef <- dlnm_results$pooled$pooled_coef
n_params <- dlnm_results$pooled$n_params
pooled_vcov <- matrix(dlnm_results$pooled$pooled_vcov, nrow = n_params, ncol = n_params)
pooled_mmt <- dlnm_results$pooled$mmt

cat("  Pooled MMT:", round(pooled_mmt, 2), "°C\n")
cat("  N params:", n_params, "\n")

# -----------------------------------------------------------------------------
# 2. Load Daily Data
# -----------------------------------------------------------------------------
cat("\n[2] Loading daily data...\n")

# Mortality data
mort_file <- file.path(DATA_DIR, paste0("mortality_", EXPOSURE_TYPE, "_daily_elderly.parquet"))
mort <- as.data.table(read_parquet(mort_file))

# ERA5 temperature data
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
data[, year := year(date)]

cat("  Rows:", format(nrow(data), big.mark = ","), "\n")
cat("  Regions:", uniqueN(data$region_code), "\n")
cat("  Years:", min(data$year), "to", max(data$year), "\n")

# Load population data
ses_file <- file.path(DATA_DIR, paste0("ses_", EXPOSURE_TYPE, "_covariates.csv"))
ses <- fread(ses_file)
if (EXPOSURE_TYPE == "immediate") {
  setnames(ses, "immediate_code", "region_code", skip_absent = TRUE)
} else {
  setnames(ses, "intermediate_code", "region_code", skip_absent = TRUE)
}
pop_map <- setNames(ses$pop_elderly, as.character(ses$region_code))

# -----------------------------------------------------------------------------
# 3. Create Prediction Cross-Basis
# -----------------------------------------------------------------------------
cat("\n[3] Setting up cross-basis for prediction...\n")

# Get unique temperatures in data
temp_unique <- sort(unique(round(data$tmean, 2)))
temp_for_pred <- seq(min(temp_unique), max(temp_unique), by = 0.1)

# Create cross-basis matching the fitted model
cb_pred <- crossbasis(
  temp_for_pred,
  lag = MAX_LAG,
  argvar = list(fun = "ns", knots = TEMP_KNOTS, Boundary.knots = TEMP_BOUNDARY),
  arglag = list(fun = "ns", df = LAG_DF)
)

# Get pooled predictions centered at MMT
cp_pooled <- crosspred(
  cb_pred,
  coef = pooled_coef,
  vcov = pooled_vcov,
  model.link = "log",
  at = temp_for_pred,
  cen = pooled_mmt
)

# Create lookup table: temperature -> RR
rr_lookup <- data.table(
  temp = temp_for_pred,
  rr = cp_pooled$allRRfit,
  rr_lo = cp_pooled$allRRlow,
  rr_hi = cp_pooled$allRRhigh
)
setkey(rr_lookup, temp)

cat("  Created RR lookup table:", nrow(rr_lookup), "temperature points\n")
cat("  RR range:", round(min(rr_lookup$rr), 4), "to", round(max(rr_lookup$rr), 4), "\n")

# -----------------------------------------------------------------------------
# 4. Calculate Attributable Burden by Region
# -----------------------------------------------------------------------------
cat("\n[4] Calculating attributable burden by region...\n")

# Function to get RR for a temperature (interpolate from lookup)
get_rr <- function(temp_val) {
  # Find nearest temperature in lookup
  idx <- which.min(abs(rr_lookup$temp - temp_val))
  return(rr_lookup$rr[idx])
}

# Vectorized RR lookup using approx
rr_func <- approxfun(rr_lookup$temp, rr_lookup$rr, rule = 2)

# Calculate RR and AF for all rows
data[, rr := rr_func(tmean)]
data[, af := fifelse(rr > 1, (rr - 1) / rr, 0)]  # Only count when RR > 1
data[, an := af * deaths]  # Attributable number

# Classify heat vs cold
data[, is_heat := tmean > pooled_mmt]
data[, is_cold := tmean < pooled_mmt]

# Percentile thresholds
data[, `:=`(
  p1 = quantile(tmean, 0.01, na.rm = TRUE),
  p2_5 = quantile(tmean, 0.025, na.rm = TRUE),
  p97_5 = quantile(tmean, 0.975, na.rm = TRUE),
  p99 = quantile(tmean, 0.99, na.rm = TRUE)
), by = region_code]

data[, is_extreme_heat := tmean > p97_5]
data[, is_extreme_cold := tmean < p2_5]
data[, is_heat_99 := tmean > p99]
data[, is_cold_1 := tmean < p1]

# Calculate burden by region
region_burden <- data[, .(
  n_days = .N,
  n_years = .N / 365.25,
  total_deaths = sum(deaths),
  
  # Total heat (> MMT)
  total_heat_an = sum(an[is_heat == TRUE], na.rm = TRUE),
  
  # Total cold (< MMT)
  total_cold_an = sum(an[is_cold == TRUE], na.rm = TRUE),
  
  # Extreme thresholds
  heat_an_97_5 = sum(an[is_extreme_heat == TRUE], na.rm = TRUE),
  cold_an_2_5 = sum(an[is_extreme_cold == TRUE], na.rm = TRUE),
  heat_an_99 = sum(an[is_heat_99 == TRUE], na.rm = TRUE),
  cold_an_1 = sum(an[is_cold_1 == TRUE], na.rm = TRUE),
  
  # Day counts
  n_heat_days = sum(is_heat),
  n_cold_days = sum(is_cold),
  n_extreme_heat_days = sum(is_extreme_heat),
  n_extreme_cold_days = sum(is_extreme_cold),
  
  # Temperature stats
  mean_temp = mean(tmean),
  p1_temp = unique(p1)[1],
  p99_temp = unique(p99)[1]
  
), by = region_code]

# Add percentages
region_burden[, `:=`(
  total_heat_af_pct = total_heat_an / total_deaths * 100,
  total_cold_af_pct = total_cold_an / total_deaths * 100,
  heat_af_pct_97_5 = heat_an_97_5 / total_deaths * 100,
  cold_af_pct_2_5 = cold_an_2_5 / total_deaths * 100,
  heat_af_pct_99 = heat_an_99 / total_deaths * 100,
  cold_af_pct_1 = cold_an_1 / total_deaths * 100
)]

# Add annual rates
region_burden[, `:=`(
  heat_annual_97_5 = heat_an_97_5 / n_years,
  cold_annual_2_5 = cold_an_2_5 / n_years,
  heat_annual_99 = heat_an_99 / n_years,
  cold_annual_1 = cold_an_1 / n_years
)]

# Add population-based rates
region_burden[, pop_elderly := pop_map[as.character(region_code)]]
region_burden[, `:=`(
  heat_rate_per_100k = fifelse(pop_elderly > 0, heat_annual_97_5 / pop_elderly * 100000, NA_real_),
  cold_rate_per_100k = fifelse(pop_elderly > 0, cold_annual_2_5 / pop_elderly * 100000, NA_real_)
)]

cat("  Processed", nrow(region_burden), "regions\n")

# -----------------------------------------------------------------------------
# 5. National Aggregation
# -----------------------------------------------------------------------------
cat("\n[5] National aggregation...\n")

national <- list(
  # Totals
  n_regions = nrow(region_burden),
  total_days = sum(region_burden$n_days),
  n_years = sum(region_burden$n_days) / uniqueN(data$region_code) / 365.25,
  total_deaths = sum(region_burden$total_deaths),
  total_pop_elderly = sum(region_burden$pop_elderly, na.rm = TRUE),
  
  # Heat burden (all heat days > MMT)
  total_heat_an = sum(region_burden$total_heat_an),
  total_heat_af_pct = sum(region_burden$total_heat_an) / sum(region_burden$total_deaths) * 100,
  
  # Cold burden (all cold days < MMT)
  total_cold_an = sum(region_burden$total_cold_an),
  total_cold_af_pct = sum(region_burden$total_cold_an) / sum(region_burden$total_deaths) * 100,
  
  # Extreme heat (P97.5)
  heat_an_97_5 = sum(region_burden$heat_an_97_5),
  heat_af_pct_97_5 = sum(region_burden$heat_an_97_5) / sum(region_burden$total_deaths) * 100,
  
  # Extreme cold (P2.5)
  cold_an_2_5 = sum(region_burden$cold_an_2_5),
  cold_af_pct_2_5 = sum(region_burden$cold_an_2_5) / sum(region_burden$total_deaths) * 100,
  
  # P99/P1 thresholds
  heat_an_99 = sum(region_burden$heat_an_99),
  heat_af_pct_99 = sum(region_burden$heat_an_99) / sum(region_burden$total_deaths) * 100,
  cold_an_1 = sum(region_burden$cold_an_1),
  cold_af_pct_1 = sum(region_burden$cold_an_1) / sum(region_burden$total_deaths) * 100,
  
  # Total non-optimal temperature burden
  total_an = sum(region_burden$total_heat_an) + sum(region_burden$total_cold_an),
  total_af_pct = (sum(region_burden$total_heat_an) + sum(region_burden$total_cold_an)) / sum(region_burden$total_deaths) * 100,
  
  # Annual averages
  annual_heat_deaths = sum(region_burden$total_heat_an) / (sum(region_burden$n_days) / uniqueN(data$region_code) / 365.25),
  annual_cold_deaths = sum(region_burden$total_cold_an) / (sum(region_burden$n_days) / uniqueN(data$region_code) / 365.25),
  
  # Rates per 100k
  heat_rate_per_100k = sum(region_burden$total_heat_an) / sum(region_burden$n_years) / sum(region_burden$pop_elderly, na.rm = TRUE) * 100000,
  cold_rate_per_100k = sum(region_burden$total_cold_an) / sum(region_burden$n_years) / sum(region_burden$pop_elderly, na.rm = TRUE) * 100000
)

# Calculate annual burden by year
annual_burden <- data[, .(
  total_deaths = sum(deaths),
  heat_an = sum(an[is_heat == TRUE], na.rm = TRUE),
  cold_an = sum(an[is_cold == TRUE], na.rm = TRUE)
), by = year]

annual_burden[, `:=`(
  heat_af_pct = heat_an / total_deaths * 100,
  cold_af_pct = cold_an / total_deaths * 100,
  total_an = heat_an + cold_an,
  total_af_pct = (heat_an + cold_an) / total_deaths * 100
)]

# -----------------------------------------------------------------------------
# 6. Print Results
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("ATTRIBUTABLE BURDEN RESULTS (", EXPOSURE_TYPE, ")\n", sep = "")
cat("=======================================================\n")

cat("\n--- NATIONAL SUMMARY ---\n")
cat("Total deaths (elderly):", format(national$total_deaths, big.mark = ","), "\n")
cat("Years of data:", round(national$n_years, 1), "\n")
cat("Pooled MMT:", round(pooled_mmt, 1), "°C\n")

cat("\n--- HEAT BURDEN (temp > MMT) ---\n")
cat("Total heat deaths:", format(round(national$total_heat_an), big.mark = ","), "\n")
cat("Heat AF:", round(national$total_heat_af_pct, 2), "%\n")
cat("Annual heat deaths:", format(round(national$annual_heat_deaths), big.mark = ","), "\n")

cat("\n--- COLD BURDEN (temp < MMT) ---\n")
cat("Total cold deaths:", format(round(national$total_cold_an), big.mark = ","), "\n")
cat("Cold AF:", round(national$total_cold_af_pct, 2), "%\n")
cat("Annual cold deaths:", format(round(national$annual_cold_deaths), big.mark = ","), "\n")

cat("\n--- TOTAL NON-OPTIMAL TEMPERATURE ---\n")
cat("Total attributable deaths:", format(round(national$total_an), big.mark = ","), "\n")
cat("Total AF:", round(national$total_af_pct, 2), "%\n")

cat("\n--- EXTREME THRESHOLDS ---\n")
cat("Heat (P97.5):", format(round(national$heat_an_97_5), big.mark = ","), 
    "deaths (", round(national$heat_af_pct_97_5, 2), "%)\n")
cat("Cold (P2.5):", format(round(national$cold_an_2_5), big.mark = ","),
    "deaths (", round(national$cold_af_pct_2_5, 2), "%)\n")
cat("Heat (P99):", format(round(national$heat_an_99), big.mark = ","),
    "deaths (", round(national$heat_af_pct_99, 2), "%)\n")
cat("Cold (P1):", format(round(national$cold_an_1), big.mark = ","),
    "deaths (", round(national$cold_af_pct_1, 2), "%)\n")

cat("\n--- ANNUAL TREND ---\n")
print(annual_burden[order(year)])

# -----------------------------------------------------------------------------
# 7. Save Results
# -----------------------------------------------------------------------------
cat("\n[7] Saving results...\n")

# Prepare output
output <- list(
  exposure_type = EXPOSURE_TYPE,
  analysis_date = as.character(Sys.time()),
  pooled_mmt = pooled_mmt,
  national = national,
  annual_burden = as.list(annual_burden),
  region_burden = as.data.frame(region_burden)
)

# Save JSON
output_file <- file.path(OUTPUT_DIR, paste0("attributable_burden_r_", EXPOSURE_TYPE, ".json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("  JSON saved to:", output_file, "\n")

# Save region-level CSV
csv_file <- file.path(OUTPUT_DIR, paste0("attributable_burden_r_", EXPOSURE_TYPE, "_regions.csv"))
fwrite(region_burden, csv_file)
cat("  CSV saved to:", csv_file, "\n")

# Save annual CSV
annual_csv <- file.path(OUTPUT_DIR, paste0("attributable_burden_r_", EXPOSURE_TYPE, "_annual.csv"))
fwrite(annual_burden, annual_csv)
cat("  Annual CSV saved to:", annual_csv, "\n")

cat("\n=======================================================\n")
cat("Done!", as.character(Sys.time()), "\n")
cat("=======================================================\n")
