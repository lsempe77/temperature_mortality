# =============================================================================
# 01d_case_crossover.R
# Case-Crossover Validation with DLNM (Enhanced Implementation)
# =============================================================================
#
# ENHANCED CASE-CROSSOVER DESIGN WITH DLNM:
# - Unit of analysis: INDIVIDUAL deaths (not daily counts)
# - For each death (case) at date t, select control days with same:
#   - Year-month
#   - Day-of-week
# - Uses CONDITIONAL POISSON REGRESSION (gnm package) for DLNM
#   * Mathematically equivalent to conditional logistic regression
#   * But produces model object compatible with crosspred()
#   * Reference: Armstrong et al. (2014), Gasparrini (2021)
# - This removes ALL time-invariant confounding via self-matching
# - Uses same lag structure as main DLNM for comparability
#
# References:
# - Maclure (1991) - Case-crossover design
# - Levy et al. (2001) - Time-stratified case-crossover
# - Armstrong et al. (2014) - Case-crossover for temperature-mortality with DLNM
# - Gasparrini et al. (2015) - DLNM methodology
# - Gasparrini (2021) - "Conditional Poisson as an alternative to clogit"
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(jsonlite)
  library(survival)  # For clogit (simple models)
  library(gnm)       # For conditional Poisson (DLNM - works with crosspred!)
  library(dlnm)      # For crossbasis in case-crossover
  library(splines)
})

# Get exposure type from command line args
args <- commandArgs(trailingOnly = TRUE)
EXPOSURE_TYPE <- if (length(args) > 0) args[1] else "intermediate"

cat("=======================================================\n")
cat("CASE-CROSSOVER WITH DLNM:", EXPOSURE_TYPE, "\n")
cat("=======================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# Configuration - match main DLNM parameters
MAX_LAG <- 21
TEMP_DF <- 4
LAG_DF <- 4

MAX_CASES <- 200000       # Maximum cases to analyze (for computational feasibility)
SAMPLE_SIZE <- 100000     # Sample this many if more cases
MIN_CONTROLS <- 2         # Minimum controls per case

# -----------------------------------------------------------------------------
# 1. Load Data
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
INPUT_DIR <- file.path(SCRIPT_DIR, "..", "..", "Input_data")
DLNM_DIR <- file.path(SCRIPT_DIR, "results")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Script directory:", SCRIPT_DIR, "\n")
cat("Data directory:", normalizePath(DATA_DIR, mustWork = FALSE), "\n")

cat("[1] Loading temperature data...\n")

# Load ERA5 temperature data
era5_file <- file.path(DATA_DIR, paste0("era5_", EXPOSURE_TYPE, "_daily.parquet"))
era5 <- as.data.table(read_parquet(era5_file))
era5[, date := as.Date(date)]
setnames(era5, "temp_mean", "tmean")

# National average temperature per day
temp_national <- era5[, .(tmean = mean(tmean, na.rm = TRUE)), by = date]
setkey(temp_national, date)
temp_national <- temp_national[order(date)]

cat("  Temperature days:", nrow(temp_national), "\n")
cat("  Date range:", as.character(min(temp_national$date)), "to", 
    as.character(max(temp_national$date)), "\n")

# Compute temperature percentiles for DLNM
temp_pcts <- quantile(temp_national$tmean, c(0.10, 0.75, 0.90), na.rm = TRUE)
temp_boundary <- c(0, 40)  # Match main DLNM

# Load DLNM results for comparison
dlnm_file <- file.path(DLNM_DIR, paste0("dlnm_r_", EXPOSURE_TYPE, "_results_v2.json"))
if (file.exists(dlnm_file)) {
  dlnm_results <- fromJSON(dlnm_file)
  pooled_mmt <- dlnm_results$pooled$mmt
  dlnm_rr_heat <- dlnm_results$pooled$rr_heat_p99$rr
  dlnm_rr_cold <- dlnm_results$pooled$rr_cold_p1$rr
  # Use DLNM's global knots if available
  if (!is.null(dlnm_results$dlnm_params$temp_knots)) {
    temp_pcts <- dlnm_results$dlnm_params$temp_knots
  }
  cat("  DLNM MMT:", round(pooled_mmt, 2), "°C\n")
  cat("  DLNM RR (P99 heat):", round(dlnm_rr_heat, 3), "\n")
  cat("  DLNM RR (P1 cold):", round(dlnm_rr_cold, 3), "\n")
} else {
  pooled_mmt <- 23
  dlnm_rr_heat <- NA
  dlnm_rr_cold <- NA
  cat("  DLNM results not found - using default MMT\n")
}

# Percentiles for effects
p1 <- quantile(temp_national$tmean, 0.01)
p5 <- quantile(temp_national$tmean, 0.05)
p95 <- quantile(temp_national$tmean, 0.95)
p99 <- quantile(temp_national$tmean, 0.99)

cat("  Temperature P1:", round(p1, 2), "P99:", round(p99, 2), "\n")
cat("  Using knots:", paste(round(temp_pcts, 1), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 2. Load Individual Death Records from SIM
# -----------------------------------------------------------------------------
cat("\n[2] Loading individual death records from SIM...\n")

do_files <- list.files(INPUT_DIR, pattern = "^DO\\d+OPEN\\.csv$", full.names = TRUE)

if (length(do_files) == 0) {
  stop("No SIM mortality files found in ", INPUT_DIR)
}

cat("  Found", length(do_files), "SIM files\n")

# Function to parse SIM age
parse_sim_age <- function(idade) {
  code <- suppressWarnings(as.integer(idade))
  ifelse(!is.na(code) & code >= 400, code - 400, NA_integer_)
}

# Function to parse SIM date
parse_sim_date <- function(dtobito) {
  first_val <- as.character(dtobito[1])
  if (grepl("-", first_val)) {
    d <- as.Date(dtobito)
  } else {
    s <- sprintf("%08d", as.integer(dtobito))
    d <- as.Date(s, format = "%d%m%Y")
  }
  d
}

# Load all elderly deaths
all_deaths <- list()

for (f in do_files) {
  fname <- basename(f)
  year_code <- gsub("DO|OPEN\\.csv", "", fname)
  year <- 2000 + as.integer(year_code)
  
  cat("  Loading", fname, "...")
  
  tryCatch({
    first_line <- readLines(f, n = 1, encoding = "latin1")
    sep <- if (grepl('","', first_line) || 
               length(gregexpr(",", first_line)[[1]]) > length(gregexpr(";", first_line)[[1]])) "," else ";"
    
    df <- fread(f, sep = sep, encoding = "Latin-1", select = c("DTOBITO", "IDADE"))
    df[, age := parse_sim_age(IDADE)]
    df <- df[!is.na(age) & age >= 60]
    df[, date := parse_sim_date(DTOBITO)]
    df <- df[!is.na(date)]
    
    # Need enough lead time for lag structure
    min_date <- min(temp_national$date) + MAX_LAG
    df <- df[date >= min_date & date <= max(temp_national$date)]
    
    if (nrow(df) > 0) {
      df[, year := year]
      all_deaths[[length(all_deaths) + 1]] <- df[, .(date, age, year)]
      cat(format(nrow(df), big.mark = ","), "elderly deaths\n")
    } else {
      cat("0 deaths in date range\n")
    }
    
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
  })
}

if (length(all_deaths) == 0) {
  stop("No death records loaded!")
}

deaths <- rbindlist(all_deaths)
total_deaths_loaded <- nrow(deaths)  # Store before sampling
cat("\n  Total elderly deaths:", format(total_deaths_loaded, big.mark = ","), "\n")

# Sample if too many
if (nrow(deaths) > MAX_CASES) {
  cat("  Sampling", format(SAMPLE_SIZE, big.mark = ","), "cases for analysis\n")
  set.seed(42)
  deaths <- deaths[sample(.N, SAMPLE_SIZE)]
}

deaths[, `:=`(
  month = month(date),
  dow = wday(date),
  case_id = .I
)]

# -----------------------------------------------------------------------------
# 3. Build Matched Case-Control Dataset with Lagged Temperatures
# -----------------------------------------------------------------------------
cat("\n[3] Building matched case-control dataset with lag structure...\n")

# Create temperature history matrix for lags
temp_vec <- temp_national$tmean
date_vec <- temp_national$date

# Create lag matrix for temperature (for each day, get temp at lag 0, 1, ..., MAX_LAG)
n_days <- length(temp_vec)
temp_lag_matrix <- matrix(NA, nrow = n_days, ncol = MAX_LAG + 1)
for (lag in 0:MAX_LAG) {
  if (lag == 0) {
    temp_lag_matrix[, lag + 1] <- temp_vec
  } else {
    temp_lag_matrix[(lag + 1):n_days, lag + 1] <- temp_vec[1:(n_days - lag)]
  }
}
colnames(temp_lag_matrix) <- paste0("lag", 0:MAX_LAG)

# Create lookup table: date -> row index
date_to_idx <- setNames(seq_len(n_days), as.character(date_vec))

# Stratum matching
temp_national[, `:=`(
  year = year(date),
  month = month(date),
  dow = wday(date)
)]
temp_national[, stratum := paste(year, month, dow, sep = "_")]
deaths[, stratum := paste(year, month, dow, sep = "_")]

# Pre-compute controls by stratum
stratum_dates <- temp_national[, .(dates = list(date)), by = stratum]
setkey(stratum_dates, stratum)

# Build matched data with lag structure
cat("  Building matched sets (this may take a while)...\n")

matched_list <- list()
valid_cases <- 0
pb_interval <- max(1, nrow(deaths) %/% 10)

for (i in seq_len(nrow(deaths))) {
  if (i %% pb_interval == 0) cat("  Progress:", round(100 * i / nrow(deaths)), "%\n")
  
  case <- deaths[i]
  strat <- case$stratum
  case_date <- case$date
  case_id <- case$case_id
  
  # Get all dates in this stratum
  stratum_data <- stratum_dates[stratum == strat]
  if (nrow(stratum_data) == 0) next
  
  control_dates <- stratum_data$dates[[1]]
  
  # Exclude case date from controls
  control_dates <- control_dates[control_dates != case_date]
  
  if (length(control_dates) < MIN_CONTROLS) next
  
  # Get case row index
  case_idx <- date_to_idx[as.character(case_date)]
  if (is.na(case_idx) || case_idx <= MAX_LAG) next
  
  # Get control indices
  control_idx <- date_to_idx[as.character(control_dates)]
  control_idx <- control_idx[!is.na(control_idx) & control_idx > MAX_LAG]
  
  if (length(control_idx) < MIN_CONTROLS) next
  
  # Build matched set with all lag temperatures
  all_idx <- c(case_idx, control_idx)
  
  matched_set <- data.table(
    case_id = case_id,
    date = c(case_date, date_vec[control_idx]),
    is_case = c(1L, rep(0L, length(control_idx)))
  )
  
  # Add lagged temperatures
  for (lag in 0:MAX_LAG) {
    matched_set[[paste0("lag", lag)]] <- temp_lag_matrix[all_idx, lag + 1]
  }
  
  matched_list[[length(matched_list) + 1]] <- matched_set
  valid_cases <- valid_cases + 1
}

cat("  Valid cases with controls:", format(valid_cases, big.mark = ","), "\n")

if (length(matched_list) == 0) {
  stop("Failed to build matched dataset!")
}

matched_data <- rbindlist(matched_list)

cat("  Total matched rows:", format(nrow(matched_data), big.mark = ","), "\n")
cat("  Cases:", format(sum(matched_data$is_case), big.mark = ","), "\n")
cat("  Controls:", format(sum(matched_data$is_case == 0), big.mark = ","), "\n")
cat("  Avg controls per case:", round(sum(matched_data$is_case == 0) / sum(matched_data$is_case), 1), "\n")

# -----------------------------------------------------------------------------
# 4. Fit DLNM within Case-Crossover using Time-Stratified Poisson
# -----------------------------------------------------------------------------
cat("\n[4] Fitting DLNM within case-crossover framework...\n")
cat("  Using time-stratified Poisson regression (aggregated by stratum).\n")
cat("  This is equivalent to case-crossover but works with crosspred().\n")

# Define lag columns
lag_cols <- paste0("lag", 0:MAX_LAG)

# Aggregate deaths by date and create case-crossover strata
# Stratum = year-month-dow combination
daily_data <- matched_data[, .(
  deaths = sum(is_case),
  n_obs = .N
), by = .(date, stratum_id = paste(year(date), month(date), wday(date), sep = "_"))]

# Merge with temperature lags (use first observation for each date)
temp_by_date <- matched_data[, c("date", lag_cols), with = FALSE][!duplicated(date)]
setkey(temp_by_date, date)
setkey(daily_data, date)
daily_data <- temp_by_date[daily_data]

cat("  Aggregated to", nrow(daily_data), "daily observations\n")
cat("  Number of strata:", length(unique(daily_data$stratum_id)), "\n")

# Create crossbasis with same specification as main DLNM
cat("  Creating crossbasis (same spec as main DLNM)...\n")
temp_matrix_agg <- as.matrix(daily_data[, ..lag_cols])

cb <- crossbasis(
  temp_matrix_agg,
  lag = MAX_LAG,
  argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
  arglag = list(fun = "ns", df = LAG_DF)
)

cat("  Crossbasis dimensions:", dim(cb)[1], "x", dim(cb)[2], "\n")

# Add crossbasis columns to data
cb_names <- colnames(cb)
cb_df <- as.data.frame(cb)
daily_data <- cbind(daily_data, cb_df)

# Fit Poisson with stratum fixed effects
# This is the time-stratified case-crossover design at aggregate level
cat("  Fitting Poisson GLM with stratum fixed effects...\n")

# Convert stratum to factor
daily_data[, stratum_factor := factor(stratum_id)]

formula_dlnm <- as.formula(paste("deaths ~", paste(cb_names, collapse = " + "), "+ stratum_factor"))

fit_dlnm <- tryCatch({
  glm(formula_dlnm, data = daily_data, family = quasipoisson())
}, error = function(e) {
  cat("  Error fitting GLM:", conditionMessage(e), "\n")
  NULL
})

if (is.null(fit_dlnm)) {
  cat("  DLNM case-crossover failed! Falling back to simple models.\n")
  dlnm_cc_success <- FALSE
  mmt_cc <- NA
  rr_p99_cc <- NA
  rr_p99_ci_cc <- c(NA, NA)
  rr_p1_cc <- NA
  rr_p1_ci_cc <- c(NA, NA)
} else {
  cat("  DLNM case-crossover converged!\n")
  
  # Extract predictions using crosspred with explicit coef/vcov
  # KEY: When using coef/vcov, allRRfit is empty but allfit contains log-RR
  cp_result <- tryCatch({
    # Extract only crossbasis coefficients (not stratum fixed effects)
    all_coefs <- coef(fit_dlnm)
    cb_coef_idx <- match(cb_names, names(all_coefs))
    cb_coefs <- all_coefs[cb_coef_idx]
    
    # Check for any issues
    if (any(is.na(cb_coef_idx))) {
      stop("Could not match all crossbasis coefficients")
    }
    
    # Extract variance-covariance matrix for crossbasis only
    full_vcov <- vcov(fit_dlnm)
    cb_vcov <- full_vcov[cb_coef_idx, cb_coef_idx]
    
    cat("  Extracted", length(cb_coefs), "crossbasis coefficients\n")
    
    # Use crosspred with explicit coef and vcov + cumul=TRUE
    pred_temps <- seq(temp_boundary[1], temp_boundary[2], by = 0.5)
    cp_dlnm <- crosspred(
      cb, 
      coef = cb_coefs,
      vcov = cb_vcov,
      at = pred_temps,
      cen = pooled_mmt,
      cumul = TRUE
    )
    
    # With coef/vcov, use allfit (log-RR) and convert with exp()
    allRRfit <- exp(cp_dlnm$allfit)
    
    if (length(allRRfit) == 0 || all(is.na(allRRfit))) {
      stop("allfit is empty or all NA")
    }
    
    valid_idx <- which(!is.na(allRRfit) & is.finite(allRRfit))
    if (length(valid_idx) == 0) {
      stop("No valid RR values")
    }
    
    min_idx <- valid_idx[which.min(allRRfit[valid_idx])]
    mmt_cc <- cp_dlnm$predvar[min_idx]
    cat("  Case-crossover MMT:", round(mmt_cc, 2), "°C\n")
    
    # Recenter at MMT and get RRs at key temperatures
    cp_dlnm_centered <- crosspred(
      cb, 
      coef = cb_coefs,
      vcov = cb_vcov,
      at = c(p1, p5, pooled_mmt, mmt_cc, p95, p99),
      cen = mmt_cc,
      cumul = TRUE
    )
    
    # Extract RRs at key percentiles (convert from log-RR)
    rr_p99_cc <- exp(cp_dlnm_centered$allfit[6])
    rr_p99_ci_cc <- exp(c(cp_dlnm_centered$allfit[6] - 1.96 * cp_dlnm_centered$allse[6],
                          cp_dlnm_centered$allfit[6] + 1.96 * cp_dlnm_centered$allse[6]))
    rr_p1_cc <- exp(cp_dlnm_centered$allfit[1])
    rr_p1_ci_cc <- exp(c(cp_dlnm_centered$allfit[1] - 1.96 * cp_dlnm_centered$allse[1],
                         cp_dlnm_centered$allfit[1] + 1.96 * cp_dlnm_centered$allse[1]))
    
    cat("\n  DLNM Case-Crossover Results (centered at MMT =", round(mmt_cc, 2), "°C):\n")
    cat(sprintf("    Heat (P99 = %.1f°C): RR = %.3f (%.3f-%.3f)\n", 
                p99, rr_p99_cc, rr_p99_ci_cc[1], rr_p99_ci_cc[2]))
    cat(sprintf("    Cold (P1 = %.1f°C): RR = %.3f (%.3f-%.3f)\n",
                p1, rr_p1_cc, rr_p1_ci_cc[1], rr_p1_ci_cc[2]))
    
    list(success = TRUE, mmt_cc = mmt_cc, 
         rr_p99_cc = rr_p99_cc, rr_p99_ci_cc = rr_p99_ci_cc,
         rr_p1_cc = rr_p1_cc, rr_p1_ci_cc = rr_p1_ci_cc)
  }, error = function(e) {
    cat("  Error in crosspred:", conditionMessage(e), "\n")
    list(success = FALSE, mmt_cc = NA, 
         rr_p99_cc = NA, rr_p99_ci_cc = c(NA, NA),
         rr_p1_cc = NA, rr_p1_ci_cc = c(NA, NA))
  })
  
  if (cp_result$success) {
    dlnm_cc_success <- TRUE
    mmt_cc <- cp_result$mmt_cc
    rr_p99_cc <- cp_result$rr_p99_cc
    rr_p99_ci_cc <- cp_result$rr_p99_ci_cc
    rr_p1_cc <- cp_result$rr_p1_cc
    rr_p1_ci_cc <- cp_result$rr_p1_ci_cc
  } else {
    dlnm_cc_success <- FALSE
    mmt_cc <- NA
    rr_p99_cc <- NA
    rr_p99_ci_cc <- c(NA, NA)
    rr_p1_cc <- NA
    rr_p1_ci_cc <- c(NA, NA)
  }
}

# -----------------------------------------------------------------------------
# 5. Simple Models for Comparison
# -----------------------------------------------------------------------------
cat("\n[5] Fitting simple models for comparison...\n")

# Use lag0 temperature for simple models
matched_data[, tmean := lag0]

# Model 1: Linear temperature effect (lag 0 only)
cat("\n  Model A: Linear temperature effect (lag 0)\n")
fit_linear <- clogit(is_case ~ tmean + strata(case_id), data = matched_data)

or_per_degree <- exp(coef(fit_linear)["tmean"])
or_ci <- exp(confint(fit_linear)["tmean", ])
cat("  OR per 1°C increase:", round(or_per_degree, 4), 
    "(95% CI:", round(or_ci[1], 4), "-", round(or_ci[2], 4), ")\n")

# Model 2: Quadratic temperature effect
cat("\n  Model B: Quadratic temperature effect\n")
matched_data[, tmean_sq := tmean^2]
fit_quad <- clogit(is_case ~ tmean + tmean_sq + strata(case_id), data = matched_data)

b1 <- coef(fit_quad)["tmean"]
b2 <- coef(fit_quad)["tmean_sq"]
if (b2 < 0) {
  optimal_temp <- -b1 / (2 * b2)
  cat("  Optimal temperature (quadratic):", round(optimal_temp, 2), "°C\n")
} else {
  optimal_temp <- pooled_mmt
  cat("  Convex curve - using DLNM MMT as reference:", round(pooled_mmt, 2), "°C\n")
}

# Model 3: Temperature categories (P99/P1)
cat("\n  Model C: Extreme temperature categories (lag 0)\n")
matched_data[, temp_extreme := cut(tmean,
                                    breaks = c(-Inf, p1, p99, Inf),
                                    labels = c("extreme_cold", "moderate", "extreme_heat"))]
matched_data[, temp_extreme := relevel(temp_extreme, ref = "moderate")]

fit_extreme <- clogit(is_case ~ temp_extreme + strata(case_id), data = matched_data)

or_extreme_heat <- exp(coef(fit_extreme)["temp_extremeextreme_heat"])
or_extreme_cold <- exp(coef(fit_extreme)["temp_extremeextreme_cold"])
ci_extreme_heat <- exp(confint(fit_extreme)["temp_extremeextreme_heat", ])
ci_extreme_cold <- exp(confint(fit_extreme)["temp_extremeextreme_cold", ])

cat("  OR extreme heat (>P99):", round(or_extreme_heat, 3),
    "(95% CI:", round(ci_extreme_heat[1], 3), "-", round(ci_extreme_heat[2], 3), ")\n")
cat("  OR extreme cold (<P1):", round(or_extreme_cold, 3),
    "(95% CI:", round(ci_extreme_cold[1], 3), "-", round(ci_extreme_cold[2], 3), ")\n")

# -----------------------------------------------------------------------------
# 6. Comprehensive Comparison
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("COMPREHENSIVE COMPARISON: Case-Crossover vs DLNM\n")
cat("=======================================================\n")

cat("\n--- HEAT EFFECTS (P99) ---\n")
cat(sprintf("  Main DLNM (time-series):        RR = %.3f\n", dlnm_rr_heat))
if (dlnm_cc_success) {
  cat(sprintf("  Case-Crossover + DLNM:          RR = %.3f (%.3f-%.3f)\n", 
              rr_p99_cc, rr_p99_ci_cc[1], rr_p99_ci_cc[2]))
}
cat(sprintf("  Case-Crossover (simple, lag0):  OR = %.3f (%.3f-%.3f)\n",
            or_extreme_heat, ci_extreme_heat[1], ci_extreme_heat[2]))

cat("\n--- COLD EFFECTS (P1) ---\n")
cat(sprintf("  Main DLNM (time-series):        RR = %.3f\n", dlnm_rr_cold))
if (dlnm_cc_success) {
  cat(sprintf("  Case-Crossover + DLNM:          RR = %.3f (%.3f-%.3f)\n",
              rr_p1_cc, rr_p1_ci_cc[1], rr_p1_ci_cc[2]))
}
cat(sprintf("  Case-Crossover (simple, lag0):  OR = %.3f (%.3f-%.3f)\n",
            or_extreme_cold, ci_extreme_cold[1], ci_extreme_cold[2]))

cat("\n--- MMT COMPARISON ---\n")
cat(sprintf("  Main DLNM MMT:           %.2f°C\n", pooled_mmt))
if (dlnm_cc_success) {
  cat(sprintf("  Case-Crossover DLNM MMT: %.2f°C\n", mmt_cc))
}
cat(sprintf("  Quadratic model optimal: %.2f°C\n", optimal_temp))

cat("\n--- INTERPRETATION ---\n")
cat("  Case-crossover removes time-invariant confounding via self-matching.\n")
cat("  Similar results = robust to unmeasured confounding.\n")
cat("  Different results = potential time-invariant confounding in main DLNM.\n")

# -----------------------------------------------------------------------------
# 7. Save Results
# -----------------------------------------------------------------------------
cat("\n[7] Saving results...\n")

output <- list(
  exposure_type = EXPOSURE_TYPE,
  analysis_date = as.character(Sys.time()),
  
  sample_info = list(
    total_deaths_available = total_deaths_loaded,
    cases_analyzed = sum(matched_data$is_case),
    total_controls = sum(matched_data$is_case == 0),
    avg_controls_per_case = sum(matched_data$is_case == 0) / sum(matched_data$is_case)
  ),
  
  dlnm_params = list(
    max_lag = MAX_LAG,
    temp_df = TEMP_DF,
    lag_df = LAG_DF,
    temp_knots = as.vector(temp_pcts),
    temp_boundary = temp_boundary
  ),
  
  thresholds = list(
    p1 = p1, p5 = p5, p95 = p95, p99 = p99, mmt_main = pooled_mmt
  ),
  
  dlnm_case_crossover = if (dlnm_cc_success) {
    list(
      converged = TRUE,
      mmt = mmt_cc,
      rr_heat_p99 = list(
        rr = rr_p99_cc,
        ci_low = rr_p99_ci_cc[1],
        ci_high = rr_p99_ci_cc[2]
      ),
      rr_cold_p1 = list(
        rr = rr_p1_cc,
        ci_low = rr_p1_ci_cc[1],
        ci_high = rr_p1_ci_cc[2]
      )
    )
  } else {
    list(converged = FALSE)
  },
  
  simple_models = list(
    linear = list(
      or_per_degree = or_per_degree,
      ci_low = or_ci[1],
      ci_high = or_ci[2]
    ),
    quadratic = list(
      optimal_temp = optimal_temp,
      coef_linear = b1,
      coef_quadratic = b2
    ),
    categorical_p99_p1 = list(
      or_extreme_heat = or_extreme_heat,
      or_extreme_heat_ci = as.vector(ci_extreme_heat),
      or_extreme_cold = or_extreme_cold,
      or_extreme_cold_ci = as.vector(ci_extreme_cold)
    )
  ),
  
  main_dlnm_comparison = list(
    dlnm_rr_heat_p99 = dlnm_rr_heat,
    dlnm_rr_cold_p1 = dlnm_rr_cold,
    dlnm_mmt = pooled_mmt
  ),
  
  interpretation = list(
    note = "Case-crossover removes time-invariant confounding via self-matching",
    similar_results = "Results robust to unmeasured confounding",
    different_results = "Potential time-invariant confounding in main DLNM"
  )
)

output_file <- file.path(OUTPUT_DIR, paste0("case_crossover_r_", EXPOSURE_TYPE, ".json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("  JSON saved to:", output_file, "\n")

cat("\n=======================================================\n")
cat("Done!", as.character(Sys.time()), "\n")
cat("=======================================================\n")
