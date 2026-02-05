# =============================================================================
# 01b_attributable_burden_mc.R
# Attributable Burden Calculation with Monte Carlo Uncertainty Propagation
# =============================================================================
# 
# This script implements proper uncertainty quantification for attributable
# burden estimates using Monte Carlo simulation from the pooled coefficient
# distribution.
#
# Method:
# 1. Load pooled DLNM coefficients and variance-covariance matrix
# 2. Draw N_SIM samples from multivariate normal distribution
# 3. For each simulation:
#    a. Generate exposure-response curve from sampled coefficients
#    b. Calculate attributable fractions and deaths
# 4. Compute percentile-based confidence intervals (2.5%, 97.5%)
#
# =============================================================================

suppressPackageStartupMessages({
  library(dlnm)
  library(data.table)
  library(arrow)
  library(jsonlite)
  library(splines)
  library(MASS)  # For mvrnorm
})

# Configuration
N_SIM <- 1000  # Number of Monte Carlo simulations
set.seed(12345)  # For reproducibility

# Get exposure type from command line args
args <- commandArgs(trailingOnly = TRUE)
EXPOSURE_TYPE <- if (length(args) > 0) args[1] else "intermediate"

cat("=======================================================\n")
cat("ATTRIBUTABLE BURDEN WITH MONTE CARLO UNCERTAINTY\n")
cat("Exposure level:", EXPOSURE_TYPE, "\n")
cat("N simulations:", N_SIM, "\n")
cat("=======================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# 1. Load DLNM Results
# -----------------------------------------------------------------------------
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

# Standardize column names
region_col <- paste0(EXPOSURE_TYPE, "_code")
if (region_col %in% names(mort)) {
  setnames(mort, region_col, "region_code")
}
if (region_col %in% names(era5)) {
  setnames(era5, region_col, "region_code")
}

# Standardize temperature column name
if ("temp_mean" %in% names(era5) && !"tmean" %in% names(era5)) {
  setnames(era5, "temp_mean", "tmean")
}

# Standardize deaths column name
if ("deaths_elderly" %in% names(mort) && !"deaths" %in% names(mort)) {
  setnames(mort, "deaths_elderly", "deaths")
}

# Standardize date columns to Date type
mort[, date := as.Date(date)]
era5[, date := as.Date(date)]

# Merge mortality and temperature
data <- merge(mort, era5[, .(region_code, date, tmean)], 
              by = c("region_code", "date"), all.x = TRUE)

# Filter to valid data
data <- data[!is.na(tmean) & !is.na(deaths)]
cat("  Merged data:", format(nrow(data), big.mark = ","), "region-days\n")

# Load population data
pop_file <- file.path(DATA_DIR, paste0("population_", EXPOSURE_TYPE, "_elderly.parquet"))
if (file.exists(pop_file)) {
  pop <- as.data.table(read_parquet(pop_file))
  if (region_col %in% names(pop)) {
    setnames(pop, region_col, "region_code")
  }
  pop_map <- setNames(pop$pop_elderly, as.character(pop$region_code))
} else {
  pop_map <- setNames(rep(NA_real_, uniqueN(data$region_code)), 
                      as.character(unique(data$region_code)))
}

# Add year
data[, year := year(date)]

# Temperature stats
cat("  Temperature range:", round(min(data$tmean), 1), "to", 
    round(max(data$tmean), 1), "°C\n")
cat("  Total deaths:", format(sum(data$deaths), big.mark = ","), "\n")

# -----------------------------------------------------------------------------
# 3. Set up Cross-basis for Predictions
# -----------------------------------------------------------------------------
cat("\n[3] Setting up prediction cross-basis...\n")

# Create fine temperature grid for predictions
temp_range <- range(data$tmean, na.rm = TRUE)
temp_for_pred <- seq(floor(temp_range[1]), ceiling(temp_range[2]), by = 0.1)

# Create cross-basis for prediction
cb_pred <- crossbasis(
  temp_for_pred,
  lag = MAX_LAG,
  argvar = list(fun = "ns", knots = TEMP_KNOTS, Boundary.knots = TEMP_BOUNDARY),
  arglag = list(fun = "ns", df = LAG_DF)
)

# -----------------------------------------------------------------------------
# 4. Monte Carlo Simulation
# -----------------------------------------------------------------------------
cat("\n[4] Running Monte Carlo simulations (N =", N_SIM, ")...\n")

# Draw coefficient samples from multivariate normal
coef_samples <- mvrnorm(N_SIM, mu = pooled_coef, Sigma = pooled_vcov)
cat("  Generated", N_SIM, "coefficient samples\n")

# Function to calculate burden for one coefficient set
calc_burden_for_coef <- function(coef_vec, cb_pred, temp_for_pred, mmt, data, n_params) {
  
  # Create a dummy vcov (zeros) since we only need point estimates
  dummy_vcov <- matrix(0, nrow = n_params, ncol = n_params)
  
  # Get predictions centered at MMT
  cp <- crosspred(
    cb_pred,
    coef = coef_vec,
    vcov = dummy_vcov,
    model.link = "log",
    at = temp_for_pred,
    cen = mmt
  )
  
  # Create RR lookup
  rr_lookup <- data.table(
    temp = temp_for_pred,
    rr = cp$allRRfit
  )
  
  # Interpolation function
  rr_func <- approxfun(rr_lookup$temp, rr_lookup$rr, rule = 2)
  
  # Calculate RR for all data points
  rr_vals <- rr_func(data$tmean)
  
  # Calculate AF (only when RR > 1)
  af_vals <- ifelse(rr_vals > 1, (rr_vals - 1) / rr_vals, 0)
  
  # Attributable numbers
  an_vals <- af_vals * data$deaths
  
  # Heat vs cold
  is_heat <- data$tmean > mmt
  is_cold <- data$tmean < mmt
  
  # Percentile thresholds (pre-calculated)
  p1 <- quantile(data$tmean, 0.01)
  p2_5 <- quantile(data$tmean, 0.025)
  p97_5 <- quantile(data$tmean, 0.975)
  p99 <- quantile(data$tmean, 0.99)
  
  # Return burden metrics
  list(
    total_deaths = sum(data$deaths),
    
    # Total heat/cold (relative to MMT)
    heat_an_total = sum(an_vals[is_heat], na.rm = TRUE),
    cold_an_total = sum(an_vals[is_cold], na.rm = TRUE),
    
    # Extreme thresholds
    heat_an_97_5 = sum(an_vals[data$tmean > p97_5], na.rm = TRUE),
    cold_an_2_5 = sum(an_vals[data$tmean < p2_5], na.rm = TRUE),
    heat_an_99 = sum(an_vals[data$tmean > p99], na.rm = TRUE),
    cold_an_1 = sum(an_vals[data$tmean < p1], na.rm = TRUE)
  )
}

# Run simulations with progress
sim_results <- vector("list", N_SIM)
pb_interval <- max(1, N_SIM / 20)

cat("  Progress: ")
for (i in 1:N_SIM) {
  sim_results[[i]] <- calc_burden_for_coef(
    coef_vec = coef_samples[i, ],
    cb_pred = cb_pred,
    temp_for_pred = temp_for_pred,
    mmt = pooled_mmt,
    data = data,
    n_params = n_params
  )
  
  if (i %% pb_interval == 0) {
    cat(round(i / N_SIM * 100), "% ", sep = "")
  }
}
cat("Done!\n")

# Convert to data.table
sim_dt <- rbindlist(lapply(sim_results, as.data.table))

# -----------------------------------------------------------------------------
# 5. Calculate Point Estimates and Confidence Intervals
# -----------------------------------------------------------------------------
cat("\n[5] Calculating confidence intervals...\n")

# Function to get mean and percentile CIs
get_ci <- function(x) {
  c(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    lo = quantile(x, 0.025, na.rm = TRUE),
    hi = quantile(x, 0.975, na.rm = TRUE)
  )
}

# Total deaths (constant)
total_deaths <- sim_dt$total_deaths[1]
n_years <- nrow(data) / uniqueN(data$region_code) / 365.25

# Calculate CIs for each burden metric
metrics <- list(
  heat_total = get_ci(sim_dt$heat_an_total),
  cold_total = get_ci(sim_dt$cold_an_total),
  heat_97_5 = get_ci(sim_dt$heat_an_97_5),
  cold_2_5 = get_ci(sim_dt$cold_an_2_5),
  heat_99 = get_ci(sim_dt$heat_an_99),
  cold_1 = get_ci(sim_dt$cold_an_1),
  total = get_ci(sim_dt$heat_an_total + sim_dt$cold_an_total)
)

# Annual versions
annual_metrics <- list(
  heat_total = get_ci(sim_dt$heat_an_total / n_years),
  cold_total = get_ci(sim_dt$cold_an_total / n_years),
  heat_97_5 = get_ci(sim_dt$heat_an_97_5 / n_years),
  cold_2_5 = get_ci(sim_dt$cold_an_2_5 / n_years),
  heat_99 = get_ci(sim_dt$heat_an_99 / n_years),
  cold_1 = get_ci(sim_dt$cold_an_1 / n_years),
  total = get_ci((sim_dt$heat_an_total + sim_dt$cold_an_total) / n_years)
)

# AF percentages
af_metrics <- list(
  heat_total = get_ci(sim_dt$heat_an_total / total_deaths * 100),
  cold_total = get_ci(sim_dt$cold_an_total / total_deaths * 100),
  heat_97_5 = get_ci(sim_dt$heat_an_97_5 / total_deaths * 100),
  cold_2_5 = get_ci(sim_dt$cold_an_2_5 / total_deaths * 100),
  heat_99 = get_ci(sim_dt$heat_an_99 / total_deaths * 100),
  cold_1 = get_ci(sim_dt$cold_an_1 / total_deaths * 100),
  total = get_ci((sim_dt$heat_an_total + sim_dt$cold_an_total) / total_deaths * 100)
)

# -----------------------------------------------------------------------------
# 6. Print Results
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("MONTE CARLO ATTRIBUTABLE BURDEN RESULTS (", EXPOSURE_TYPE, ")\n", sep = "")
cat("=======================================================\n")
cat("N simulations:", N_SIM, "\n")
cat("Total deaths:", format(total_deaths, big.mark = ","), "\n")
cat("Years:", round(n_years, 1), "\n")
cat("Pooled MMT:", round(pooled_mmt, 1), "°C\n")

cat("\n--- HEAT BURDEN (temp > MMT) ---\n")
cat("  Total:", format(round(metrics$heat_total["mean"]), big.mark = ","),
    "(95% CI:", format(round(metrics$heat_total["lo.2.5%"]), big.mark = ","), "-",
    format(round(metrics$heat_total["hi.97.5%"]), big.mark = ","), ")\n")
cat("  Annual:", format(round(annual_metrics$heat_total["mean"]), big.mark = ","),
    "(95% CI:", format(round(annual_metrics$heat_total["lo.2.5%"]), big.mark = ","), "-",
    format(round(annual_metrics$heat_total["hi.97.5%"]), big.mark = ","), ")\n")
cat("  AF%:", round(af_metrics$heat_total["mean"], 2),
    "(95% CI:", round(af_metrics$heat_total["lo.2.5%"], 2), "-",
    round(af_metrics$heat_total["hi.97.5%"], 2), ")\n")

cat("\n--- COLD BURDEN (temp < MMT) ---\n")
cat("  Total:", format(round(metrics$cold_total["mean"]), big.mark = ","),
    "(95% CI:", format(round(metrics$cold_total["lo.2.5%"]), big.mark = ","), "-",
    format(round(metrics$cold_total["hi.97.5%"]), big.mark = ","), ")\n")
cat("  Annual:", format(round(annual_metrics$cold_total["mean"]), big.mark = ","),
    "(95% CI:", format(round(annual_metrics$cold_total["lo.2.5%"]), big.mark = ","), "-",
    format(round(annual_metrics$cold_total["hi.97.5%"]), big.mark = ","), ")\n")
cat("  AF%:", round(af_metrics$cold_total["mean"], 2),
    "(95% CI:", round(af_metrics$cold_total["lo.2.5%"], 2), "-",
    round(af_metrics$cold_total["hi.97.5%"], 2), ")\n")

cat("\n--- TOTAL NON-OPTIMAL TEMPERATURE ---\n")
cat("  Total:", format(round(metrics$total["mean"]), big.mark = ","),
    "(95% CI:", format(round(metrics$total["lo.2.5%"]), big.mark = ","), "-",
    format(round(metrics$total["hi.97.5%"]), big.mark = ","), ")\n")
cat("  Annual:", format(round(annual_metrics$total["mean"]), big.mark = ","),
    "(95% CI:", format(round(annual_metrics$total["lo.2.5%"]), big.mark = ","), "-",
    format(round(annual_metrics$total["hi.97.5%"]), big.mark = ","), ")\n")
cat("  AF%:", round(af_metrics$total["mean"], 2),
    "(95% CI:", round(af_metrics$total["lo.2.5%"], 2), "-",
    round(af_metrics$total["hi.97.5%"], 2), ")\n")

cat("\n--- EXTREME THRESHOLDS (P99/P1) ---\n")
cat("  Heat P99:", format(round(metrics$heat_99["mean"]), big.mark = ","),
    "(95% CI:", format(round(metrics$heat_99["lo.2.5%"]), big.mark = ","), "-",
    format(round(metrics$heat_99["hi.97.5%"]), big.mark = ","), ")\n")
cat("  Cold P1:", format(round(metrics$cold_1["mean"]), big.mark = ","),
    "(95% CI:", format(round(metrics$cold_1["lo.2.5%"]), big.mark = ","), "-",
    format(round(metrics$cold_1["hi.97.5%"]), big.mark = ","), ")\n")

cat("\n--- EXTREME THRESHOLDS (P97.5/P2.5) ---\n")
cat("  Heat P97.5:", format(round(metrics$heat_97_5["mean"]), big.mark = ","),
    "(95% CI:", format(round(metrics$heat_97_5["lo.2.5%"]), big.mark = ","), "-",
    format(round(metrics$heat_97_5["hi.97.5%"]), big.mark = ","), ")\n")
cat("  Cold P2.5:", format(round(metrics$cold_2_5["mean"]), big.mark = ","),
    "(95% CI:", format(round(metrics$cold_2_5["lo.2.5%"]), big.mark = ","), "-",
    format(round(metrics$cold_2_5["hi.97.5%"]), big.mark = ","), ")\n")

# -----------------------------------------------------------------------------
# 7. Save Results
# -----------------------------------------------------------------------------
cat("\n[7] Saving results...\n")

output <- list(
  exposure_type = EXPOSURE_TYPE,
  analysis_date = as.character(Sys.time()),
  n_simulations = N_SIM,
  seed = 12345,
  pooled_mmt = pooled_mmt,
  total_deaths = total_deaths,
  n_years = n_years,
  
  # Total burden with CIs
  burden_total = list(
    heat = as.list(metrics$heat_total),
    cold = as.list(metrics$cold_total),
    total = as.list(metrics$total)
  ),
  
  # Annual burden with CIs
  burden_annual = list(
    heat = as.list(annual_metrics$heat_total),
    cold = as.list(annual_metrics$cold_total),
    total = as.list(annual_metrics$total)
  ),
  
  # AF percentages with CIs
  af_percent = list(
    heat = as.list(af_metrics$heat_total),
    cold = as.list(af_metrics$cold_total),
    total = as.list(af_metrics$total)
  ),
  
  # Extreme thresholds
  extreme_burden = list(
    heat_p99 = as.list(metrics$heat_99),
    cold_p1 = as.list(metrics$cold_1),
    heat_p97_5 = as.list(metrics$heat_97_5),
    cold_p2_5 = as.list(metrics$cold_2_5)
  ),
  
  extreme_annual = list(
    heat_p99 = as.list(annual_metrics$heat_99),
    cold_p1 = as.list(annual_metrics$cold_1),
    heat_p97_5 = as.list(annual_metrics$heat_97_5),
    cold_p2_5 = as.list(annual_metrics$cold_2_5)
  ),
  
  extreme_af = list(
    heat_p99 = as.list(af_metrics$heat_99),
    cold_p1 = as.list(af_metrics$cold_1),
    heat_p97_5 = as.list(af_metrics$heat_97_5),
    cold_p2_5 = as.list(af_metrics$cold_2_5)
  ),
  
  # Raw simulation results (for further analysis)
  simulations = list(
    heat_total = sim_dt$heat_an_total,
    cold_total = sim_dt$cold_an_total,
    heat_99 = sim_dt$heat_an_99,
    cold_1 = sim_dt$cold_an_1
  )
)

# Save JSON
output_file <- file.path(OUTPUT_DIR, paste0("attributable_burden_mc_", EXPOSURE_TYPE, ".json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("  JSON saved to:", output_file, "\n")

# Save summary CSV
summary_df <- data.frame(
  metric = c("heat_total", "cold_total", "total", 
             "heat_p99", "cold_p1", "heat_p97_5", "cold_p2_5"),
  mean = c(metrics$heat_total["mean"], metrics$cold_total["mean"], metrics$total["mean"],
           metrics$heat_99["mean"], metrics$cold_1["mean"], 
           metrics$heat_97_5["mean"], metrics$cold_2_5["mean"]),
  ci_lo = c(metrics$heat_total["lo.2.5%"], metrics$cold_total["lo.2.5%"], metrics$total["lo.2.5%"],
            metrics$heat_99["lo.2.5%"], metrics$cold_1["lo.2.5%"],
            metrics$heat_97_5["lo.2.5%"], metrics$cold_2_5["lo.2.5%"]),
  ci_hi = c(metrics$heat_total["hi.97.5%"], metrics$cold_total["hi.97.5%"], metrics$total["hi.97.5%"],
            metrics$heat_99["hi.97.5%"], metrics$cold_1["hi.97.5%"],
            metrics$heat_97_5["hi.97.5%"], metrics$cold_2_5["hi.97.5%"]),
  annual_mean = c(annual_metrics$heat_total["mean"], annual_metrics$cold_total["mean"], 
                  annual_metrics$total["mean"], annual_metrics$heat_99["mean"], 
                  annual_metrics$cold_1["mean"], annual_metrics$heat_97_5["mean"], 
                  annual_metrics$cold_2_5["mean"]),
  annual_ci_lo = c(annual_metrics$heat_total["lo.2.5%"], annual_metrics$cold_total["lo.2.5%"],
                   annual_metrics$total["lo.2.5%"], annual_metrics$heat_99["lo.2.5%"],
                   annual_metrics$cold_1["lo.2.5%"], annual_metrics$heat_97_5["lo.2.5%"],
                   annual_metrics$cold_2_5["lo.2.5%"]),
  annual_ci_hi = c(annual_metrics$heat_total["hi.97.5%"], annual_metrics$cold_total["hi.97.5%"],
                   annual_metrics$total["hi.97.5%"], annual_metrics$heat_99["hi.97.5%"],
                   annual_metrics$cold_1["hi.97.5%"], annual_metrics$heat_97_5["hi.97.5%"],
                   annual_metrics$cold_2_5["hi.97.5%"]),
  af_mean = c(af_metrics$heat_total["mean"], af_metrics$cold_total["mean"], 
              af_metrics$total["mean"], af_metrics$heat_99["mean"],
              af_metrics$cold_1["mean"], af_metrics$heat_97_5["mean"],
              af_metrics$cold_2_5["mean"]),
  af_ci_lo = c(af_metrics$heat_total["lo.2.5%"], af_metrics$cold_total["lo.2.5%"],
               af_metrics$total["lo.2.5%"], af_metrics$heat_99["lo.2.5%"],
               af_metrics$cold_1["lo.2.5%"], af_metrics$heat_97_5["lo.2.5%"],
               af_metrics$cold_2_5["lo.2.5%"]),
  af_ci_hi = c(af_metrics$heat_total["hi.97.5%"], af_metrics$cold_total["hi.97.5%"],
               af_metrics$total["hi.97.5%"], af_metrics$heat_99["hi.97.5%"],
               af_metrics$cold_1["hi.97.5%"], af_metrics$heat_97_5["hi.97.5%"],
               af_metrics$cold_2_5["hi.97.5%"])
)

csv_file <- file.path(OUTPUT_DIR, paste0("attributable_burden_mc_", EXPOSURE_TYPE, "_summary.csv"))
fwrite(summary_df, csv_file)
cat("  Summary CSV saved to:", csv_file, "\n")

cat("\n=======================================================\n")
cat("MONTE CARLO UNCERTAINTY PROPAGATION COMPLETE\n")
cat("Finished:", as.character(Sys.time()), "\n")
cat("=======================================================\n")
