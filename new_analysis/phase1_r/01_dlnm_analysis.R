# =============================================================================
# Phase 1: DLNM Analysis using R's dlnm package
# Following Gasparrini et al. 2015 Lancet methodology
# =============================================================================

library(dlnm)
library(mvmeta)
library(splines)
library(jsonlite)
library(data.table)

# Configuration
args <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- if (length(args) >= 1) args[1] else "../phase0_data_prep/results"
OUTPUT_DIR <- if (length(args) >= 2) args[2] else "results"
EXPOSURE_TYPE <- if (length(args) >= 3) args[3] else "immediate"  # or "intermediate"

cat("=======================================================\n")
cat("DLNM Analysis - R Implementation\n")
cat("=======================================================\n")
cat("Data directory:", DATA_DIR, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Exposure type:", EXPOSURE_TYPE, "\n\n")

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------
cat("Loading data...\n")

# We need arrow package to read parquet
if (!require("arrow", quietly = TRUE)) {
  install.packages("arrow", repos = "https://cloud.r-project.org")
  library(arrow)
}

# Load mortality data (parquet format)
mortality_file <- file.path(DATA_DIR, paste0("mortality_", EXPOSURE_TYPE, "_daily_elderly.parquet"))
if (!file.exists(mortality_file)) {
  stop("Mortality file not found: ", mortality_file)
}
mortality <- as.data.table(read_parquet(mortality_file))
cat("  Mortality records:", nrow(mortality), "\n")
cat("  Mortality columns:", paste(names(mortality), collapse = ", "), "\n")

# Load weather data (parquet format)
weather_file <- file.path(DATA_DIR, paste0("era5_", EXPOSURE_TYPE, "_daily.parquet"))
if (!file.exists(weather_file)) {
  stop("Weather file not found: ", weather_file)
}
weather <- as.data.table(read_parquet(weather_file))
cat("  Weather records:", nrow(weather), "\n")
cat("  Weather columns:", paste(names(weather), collapse = ", "), "\n")

# Standardize column names
# Mortality: date, immediate_code/intermediate_code, deaths_elderly
# Weather: date, region_code, temp_mean

# Rename region code columns
region_col <- if (EXPOSURE_TYPE == "immediate") "immediate_code" else "intermediate_code"
if (region_col %in% names(mortality)) {
  setnames(mortality, region_col, "region_code")
}

# Rename temperature column
if ("temp_mean" %in% names(weather)) {
  setnames(weather, "temp_mean", "tmean")
}

# Ensure date columns are Date type
mortality[, date := as.Date(date)]
weather[, date := as.Date(date)]

# Merge by region and date
data <- merge(mortality, weather[, .(region_code, date, tmean)],
              by = c("region_code", "date"),
              all.x = TRUE)

cat("  Merged records:", nrow(data), "\n")

# Rename deaths column
if ("deaths_elderly" %in% names(data)) {
  setnames(data, "deaths_elderly", "deaths")
}

# Remove missing temperature
data <- data[!is.na(tmean)]
cat("  After removing NA temps:", nrow(data), "\n")

# -----------------------------------------------------------------------------
# 2. DLNM Specification
# -----------------------------------------------------------------------------
cat("\nSetting up DLNM specification...\n")

# Cross-basis parameters following Gasparrini 2015
MAX_LAG <- 21  # 21 days for lag
TEMP_DF <- 4   # df for temperature (natural cubic spline)
LAG_DF <- 4    # df for lag (natural cubic spline of log(lag+1))

# Temperature percentiles for knots
temp_all <- data$tmean
temp_pcts <- quantile(temp_all, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)
temp_boundary <- range(temp_all, na.rm = TRUE)

cat("  Temperature range:", round(temp_boundary[1], 1), "to", round(temp_boundary[2], 1), "\n")
cat("  Temperature knots:", paste(round(temp_pcts, 1), collapse = ", "), "\n")
cat("  Max lag:", MAX_LAG, "days\n")

# -----------------------------------------------------------------------------
# 3. Fit DLNM by Region
# -----------------------------------------------------------------------------
cat("\nFitting DLNM by region...\n")

regions <- unique(data$region_code)
cat("  Number of regions:", length(regions), "\n\n")

# Storage for results
region_results <- list()
coef_list <- list()
vcov_list <- list()
valid_regions <- c()

for (reg in regions) {
  cat("Processing region:", reg, "... ")
  
  # Subset data for this region
  reg_data <- data[region_code == reg]
  
  if (nrow(reg_data) < 365 * 2) {
    cat("SKIP (insufficient data)\n")
    next
  }
  
  # Data is already at region-day level, just ensure sorted
  daily <- reg_data[order(date)]
  setnames(daily, "deaths", "deaths", skip_absent = TRUE)
  
  if (nrow(daily) < 365 * 2) {
    cat("SKIP (insufficient daily records)\n")
    next
  }
  
  # Check temperature variation
  temp_range <- range(daily$tmean, na.rm = TRUE)
  if (diff(temp_range) < 5) {
    cat("SKIP (insufficient temp variation)\n")
    next
  }
  
  tryCatch({
    # Create cross-basis for temperature
    # argvar: natural cubic spline for temperature
    # arglag: natural cubic spline for log(lag+1)
    cb <- crossbasis(
      daily$tmean,
      lag = MAX_LAG,
      argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
      arglag = list(fun = "ns", df = LAG_DF)
    )
    
    # Add time variables for seasonality control
    daily[, `:=`(
      year = year(date),
      month = month(date),
      dow = wday(date),
      doy = yday(date),
      time = as.numeric(date - min(date))
    )]
    
    # Fit quasi-Poisson model with DLNM
    # Control for: long-term trend, seasonality, day of week
    model <- glm(
      deaths ~ cb + 
        ns(time, df = 8 * length(unique(year))) +  # Long-term trend (8 df/year)
        factor(dow),                                 # Day of week
      family = quasipoisson(link = "log"),
      data = daily,
      na.action = na.exclude
    )
    
    # Extract crossbasis coefficients and vcov
    # Get indices for crossbasis terms
    cb_idx <- grep("cb", names(coef(model)))
    
    cb_coef <- coef(model)[cb_idx]
    cb_vcov <- vcov(model)[cb_idx, cb_idx]
    
    # Find MMT (minimum mortality temperature)
    # Predict across temperature range
    temp_seq <- seq(temp_boundary[1], temp_boundary[2], length.out = 100)
    
    # Create prediction crossbasis
    cb_pred <- crossbasis(
      temp_seq,
      lag = MAX_LAG,
      argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
      arglag = list(fun = "ns", df = LAG_DF)
    )
    
    # Cumulative effect (summing over all lags)
    cp <- crosspred(cb, model, at = temp_seq, cumul = TRUE, cen = median(daily$tmean))
    
    # Find MMT (temperature with minimum cumulative RR)
    mmt_idx <- which.min(cp$allRRfit)
    mmt <- temp_seq[mmt_idx]
    
    # Re-center predictions at MMT
    cp_mmt <- crosspred(cb, model, at = temp_seq, cumul = TRUE, cen = mmt)
    
    # Get RRs at key percentiles
    reg_pcts <- quantile(daily$tmean, probs = c(0.01, 0.05, 0.50, 0.95, 0.99), na.rm = TRUE)
    
    cp_pcts <- crosspred(cb, model, at = reg_pcts, cumul = TRUE, cen = mmt)
    
    # Store results
    region_results[[reg]] <- list(
      region_code = reg,
      n_days = nrow(daily),
      n_deaths = sum(daily$deaths),
      temp_range = as.list(temp_range),
      temp_percentiles = list(
        p1 = reg_pcts[1],
        p5 = reg_pcts[2],
        p50 = reg_pcts[3],
        p95 = reg_pcts[4],
        p99 = reg_pcts[5]
      ),
      mmt = mmt,
      cb_coefs = as.list(cb_coef),
      cb_vcov = as.matrix(cb_vcov),
      rr_p99 = list(
        rr = cp_pcts$allRRfit[5],
        rr_lo = cp_pcts$allRRlow[5],
        rr_hi = cp_pcts$allRRhigh[5]
      ),
      rr_p1 = list(
        rr = cp_pcts$allRRfit[1],
        rr_lo = cp_pcts$allRRlow[1],
        rr_hi = cp_pcts$allRRhigh[1]
      ),
      crossbasis_info = list(
        temp_knots = as.list(temp_pcts),
        temp_boundary = as.list(temp_boundary),
        temp_df = TEMP_DF,
        lag_df = LAG_DF,
        max_lag = MAX_LAG,
        n_params = length(cb_coef)
      )
    )
    
    # Quality check for MVMeta inclusion:
    # 1. MMT should not be at boundary
    # 2. RR should be in reasonable range (0.5-5 for P99)
    # 3. vcov should have positive diagonal
    mmt_at_boundary <- (abs(mmt - temp_boundary[1]) < 1) | (abs(mmt - temp_boundary[2]) < 1)
    rr_reasonable <- is.finite(cp_pcts$allRRfit[5]) && 
                     cp_pcts$allRRfit[5] > 0.5 && 
                     cp_pcts$allRRfit[5] < 5
    vcov_ok <- all(diag(cb_vcov) > 0) && all(is.finite(cb_vcov))
    
    if (!mmt_at_boundary && rr_reasonable && vcov_ok) {
      # Store for MVMeta
      coef_list[[length(coef_list) + 1]] <- cb_coef
      vcov_list[[length(vcov_list) + 1]] <- cb_vcov
      valid_regions <- c(valid_regions, reg)
      cat(sprintf("OK (n=%d, MMT=%.1f, RR_P99=%.3f) [MVMETA]\n", 
                  nrow(daily), mmt, cp_pcts$allRRfit[5]))
    } else {
      cat(sprintf("OK (n=%d, MMT=%.1f, RR_P99=%.3f) [excluded from MVMeta]\n", 
                  nrow(daily), mmt, cp_pcts$allRRfit[5]))
    }
    
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
  })
}

cat("\nSuccessfully processed", length(region_results), "regions\n")
cat("Regions valid for MVMeta:", length(valid_regions), "\n")

# -----------------------------------------------------------------------------
# 4. Multivariate Meta-Analysis
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("Running Multivariate Meta-Analysis...\n")
cat("=======================================================\n")

pooled_results <- NULL

if (length(coef_list) >= 10) {
  
  # Stack coefficients
  coef_matrix <- do.call(rbind, coef_list)
  rownames(coef_matrix) <- valid_regions
  
  # Run MVMeta with error handling
  # Try REML first, then ML, then fixed effects
  mv_fit <- NULL
  mv_method <- "reml"
  
  tryCatch({
    mv_fit <- mvmeta(coef_matrix, vcov_list, method = "reml", control = list(maxiter = 200))
    mv_method <- "reml"
  }, error = function(e) {
    cat("REML failed, trying ML...\n")
    tryCatch({
      mv_fit <<- mvmeta(coef_matrix, vcov_list, method = "ml", control = list(maxiter = 200))
      mv_method <<- "ml"
    }, error = function(e) {
      cat("ML failed, trying fixed effects...\n")
      tryCatch({
        mv_fit <<- mvmeta(coef_matrix, vcov_list, method = "fixed")
        mv_method <<- "fixed"
      }, error = function(e) {
        cat("All MVMeta methods failed:", conditionMessage(e), "\n")
      })
    })
  })
  
  if (!is.null(mv_fit)) {
    cat("MVMeta method:", mv_method, "\n")
    cat("MVMeta converged:", mv_fit$converged, "\n")
    if (!is.null(mv_fit$I2)) cat("I-squared:", round(mv_fit$I2, 1), "%\n")
    
    # Get pooled coefficients
    pooled_coef <- coef(mv_fit)
    pooled_vcov <- vcov(mv_fit)
    
    cat("Pooled coefficient range:", 
        round(min(pooled_coef), 4), "to", round(max(pooled_coef), 4), "\n")
    
    # Create pooled predictions
    # We need to rebuild the crossbasis for prediction
    temp_seq <- seq(temp_boundary[1], temp_boundary[2], length.out = 100)
    
    # Get median MMT from valid regions
    mmts <- sapply(region_results[valid_regions], function(x) x$mmt)
    pooled_mmt <- median(mmts, na.rm = TRUE)
    
    cat("Pooled MMT:", round(pooled_mmt, 1), "\n")
  
    # Get pooled percentiles
    all_p1 <- sapply(region_results[valid_regions], function(x) x$temp_percentiles$p1)
    all_p5 <- sapply(region_results[valid_regions], function(x) x$temp_percentiles$p5)
    all_p50 <- sapply(region_results[valid_regions], function(x) x$temp_percentiles$p50)
    all_p95 <- sapply(region_results[valid_regions], function(x) x$temp_percentiles$p95)
    all_p99 <- sapply(region_results[valid_regions], function(x) x$temp_percentiles$p99)
    
    pooled_pcts <- list(
      p1 = median(all_p1),
      p5 = median(all_p5),
      p50 = median(all_p50),
      p95 = median(all_p95),
      p99 = median(all_p99)
    )
    
    # For pooled predictions, we use the BLUP (best linear unbiased prediction)
    # Create a reference region with average temperature distribution
    ref_tmean <- data[, mean(tmean, na.rm = TRUE), by = date]$V1
    
    cb_ref <- crossbasis(
      ref_tmean,
      lag = MAX_LAG,
      argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
      arglag = list(fun = "ns", df = LAG_DF)
    )
    
    # Predict using pooled coefficients via crosspred with coef/vcov
    cp_pooled <- crosspred(
      cb_ref,
      coef = pooled_coef,
      vcov = pooled_vcov,
      at = c(pooled_pcts$p1, pooled_pcts$p5, pooled_pcts$p50, pooled_pcts$p95, pooled_pcts$p99),
      cen = pooled_mmt
    )
    
    cat("\n=== POOLED RESULTS ===\n")
    cat(sprintf("Heat P99 (%.1f°C) vs MMT (%.1f°C): RR=%.4f (%.4f-%.4f)\n",
                pooled_pcts$p99, pooled_mmt,
                cp_pooled$allRRfit[5], cp_pooled$allRRlow[5], cp_pooled$allRRhigh[5]))
    cat(sprintf("Heat P95 (%.1f°C) vs MMT (%.1f°C): RR=%.4f (%.4f-%.4f)\n",
                pooled_pcts$p95, pooled_mmt,
                cp_pooled$allRRfit[4], cp_pooled$allRRlow[4], cp_pooled$allRRhigh[4]))
    cat(sprintf("Cold P1 (%.1f°C) vs MMT (%.1f°C):  RR=%.4f (%.4f-%.4f)\n",
                pooled_pcts$p1, pooled_mmt,
                cp_pooled$allRRfit[1], cp_pooled$allRRlow[1], cp_pooled$allRRhigh[1]))
    cat(sprintf("Cold P5 (%.1f°C) vs MMT (%.1f°C):  RR=%.4f (%.4f-%.4f)\n",
                pooled_pcts$p5, pooled_mmt,
                cp_pooled$allRRfit[2], cp_pooled$allRRlow[2], cp_pooled$allRRhigh[2]))
    
    # Store pooled results
    pooled_results <- list(
      mvmeta = list(
        converged = mv_fit$converged,
        method = mv_method,
        n_regions = length(valid_regions),
        n_params = length(pooled_coef),
        I2 = if(!is.null(mv_fit$I2)) mv_fit$I2 else NA,
        Qstat = if(!is.null(mv_fit$Qstat)) mv_fit$Qstat else NA,
        Qdf = if(!is.null(mv_fit$Qdf)) mv_fit$Qdf else NA,
        Qpval = if(!is.null(mv_fit$pvalue)) mv_fit$pvalue else NA
      ),
      pooled_coef = as.list(pooled_coef),
      pooled_vcov = as.matrix(pooled_vcov),
      pooled_mmt = pooled_mmt,
      pooled_percentiles = pooled_pcts,
      pooled_rr = list(
        heat_p99 = list(
          rr = cp_pooled$allRRfit[5],
          rr_lo = cp_pooled$allRRlow[5],
          rr_hi = cp_pooled$allRRhigh[5],
          temp = pooled_pcts$p99
        ),
        heat_p95 = list(
          rr = cp_pooled$allRRfit[4],
          rr_lo = cp_pooled$allRRlow[4],
          rr_hi = cp_pooled$allRRhigh[4],
          temp = pooled_pcts$p95
        ),
        cold_p1 = list(
          rr = cp_pooled$allRRfit[1],
          rr_lo = cp_pooled$allRRlow[1],
          rr_hi = cp_pooled$allRRhigh[1],
          temp = pooled_pcts$p1
        ),
        cold_p5 = list(
          rr = cp_pooled$allRRfit[2],
          rr_lo = cp_pooled$allRRlow[2],
          rr_hi = cp_pooled$allRRhigh[2],
          temp = pooled_pcts$p5
        )
      ),
      crossbasis_info = list(
        temp_knots = as.list(temp_pcts),
        temp_boundary = as.list(temp_boundary),
        temp_df = TEMP_DF,
        lag_df = LAG_DF,
        max_lag = MAX_LAG
      )
    )
  }  # End of !is.null(mv_fit)
  
} else {
  cat("Not enough regions for meta-analysis (need >= 10, have", length(coef_list), ")\n")
}

# -----------------------------------------------------------------------------
# 5. Save Results
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("Saving results...\n")
cat("=======================================================\n")

output <- list(
  exposure_type = EXPOSURE_TYPE,
  analysis_date = as.character(Sys.time()),
  n_regions_total = length(region_results),
  n_regions_mvmeta = length(valid_regions),
  region_results = region_results,
  pooled = pooled_results
)

output_file <- file.path(OUTPUT_DIR, paste0("dlnm_r_", EXPOSURE_TYPE, "_results.json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 8)
cat("Results saved to:", output_file, "\n")

# Save summary CSV for all processed regions
if (length(region_results) > 0) {
  all_region_codes <- names(region_results)
  summary_df <- data.frame(
    region = all_region_codes,
    n_days = sapply(region_results, function(x) x$n_days),
    n_deaths = sapply(region_results, function(x) x$n_deaths),
    mmt = sapply(region_results, function(x) x$mmt),
    rr_p99 = sapply(region_results, function(x) x$rr_p99$rr),
    rr_p99_lo = sapply(region_results, function(x) x$rr_p99$rr_lo),
    rr_p99_hi = sapply(region_results, function(x) x$rr_p99$rr_hi),
    rr_p1 = sapply(region_results, function(x) x$rr_p1$rr),
    rr_p1_lo = sapply(region_results, function(x) x$rr_p1$rr_lo),
    rr_p1_hi = sapply(region_results, function(x) x$rr_p1$rr_hi),
    in_mvmeta = all_region_codes %in% valid_regions
  )
  
  summary_file <- file.path(OUTPUT_DIR, paste0("dlnm_r_", EXPOSURE_TYPE, "_summary.csv"))
  fwrite(summary_df, summary_file)
  cat("Summary saved to:", summary_file, "\n")
  cat("Regions in MVMeta:", sum(summary_df$in_mvmeta), "of", nrow(summary_df), "\n")
}

cat("\nDone!\n")
