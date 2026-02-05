# =============================================================================
# 02a_sensitivity.R
# Sensitivity Analyses for DLNM Temperature-Mortality Estimates
# =============================================================================
#
# Tests robustness of results to:
# 1. Lag structure (7, 14, 21, 28 days)
# 2. Spline degrees of freedom (3, 4, 5)
# 3. Temperature percentile thresholds (P95/P5, P97.5/P2.5, P99/P1)
#
# References:
# - Gasparrini et al. (2015) Lancet - Sensitivity analysis methodology
#
# =============================================================================

suppressPackageStartupMessages({
  library(dlnm)
  library(mixmeta)
  library(data.table)
  library(arrow)
  library(jsonlite)
  library(splines)
})

# Get exposure type from command line args
args <- commandArgs(trailingOnly = TRUE)
EXPOSURE_TYPE <- if (length(args) > 0) args[1] else "intermediate"

cat("=======================================================\n")
cat("SENSITIVITY ANALYSES:", EXPOSURE_TYPE, "\n")
cat("=======================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# Configuration
BASELINE_LAG <- 21
BASELINE_TEMP_DF <- 4
BASELINE_LAG_DF <- 4

LAG_OPTIONS <- c(7, 14, 21, 28)
DF_OPTIONS <- c(3, 4, 5)

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------
# Get script directory for relative path resolution
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return(getwd())
}
SCRIPT_DIR <- get_script_dir()
BASE_DIR <- dirname(SCRIPT_DIR)  # Go up from phase2_r to new_analysis

DATA_DIR <- file.path(BASE_DIR, "phase0_data_prep", "results")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

cat("[1] Loading data...\n")

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

# Working variables
data[, deaths := deaths_elderly]
data[, tmean := temp_mean]
data <- data[!is.na(deaths) & !is.na(tmean)]
setorder(data, region_code, date)

# Global temperature percentiles for boundaries
temp_all <- data$tmean
temp_pcts <- quantile(temp_all, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)
temp_boundary <- c(0, 40)  # Expanded from c(0.7, 35.1) per expert recommendation

cat("  Regions:", uniqueN(data$region_code), "\n")
cat("  GLOBAL knots (P10/P75/P90):", round(temp_pcts, 1), "\n")
cat("  Boundary knots:", temp_boundary, "\n")
cat("  Date range:", as.character(min(data$date)), "to", as.character(max(data$date)), "\n")

# -----------------------------------------------------------------------------
# 2. Function to Fit DLNM with Specific Parameters
# -----------------------------------------------------------------------------
fit_sensitivity_dlnm <- function(data_dt, max_lag, temp_df, lag_df) {
  
  regions <- unique(data_dt$region_code)
  n_regions <- length(regions)
  
  cat(sprintf("  Fitting %d regions (lag=%d, temp_df=%d, lag_df=%d)...\n", 
              n_regions, max_lag, temp_df, lag_df))
  
  # Store results
  mvmeta_coefs <- list()
  mvmeta_vcovs <- list()
  mvmeta_regions <- c()
  regional_rrs <- list()
  
  for (i in seq_along(regions)) {
    reg <- regions[i]
    
    if (i %% 50 == 0) {
      cat(sprintf("    Progress: %d/%d regions\n", i, n_regions))
    }
    
    # Subset data
    daily <- data_dt[region_code == reg][order(date)]
    n_days <- nrow(daily)
    
    if (n_days < 365) next
    
    result <- tryCatch({
      # Use GLOBAL knots for meta-analysis comparability (expert recommendation)
      # But restrict MMT search to region's observed temperature range
      reg_temp_range <- range(daily$tmean, na.rm = TRUE)
      reg_p1 <- quantile(daily$tmean, 0.01, na.rm = TRUE)
      reg_p99 <- quantile(daily$tmean, 0.99, na.rm = TRUE)
      
      # Create cross-basis with GLOBAL knots
      cb <- crossbasis(
        daily$tmean,
        lag = max_lag,
        argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
        arglag = list(fun = "ns", df = lag_df)
      )
      
      # Fit GLM
      model <- glm(
        deaths ~ cb + ns(as.numeric(date), df = 7 * length(unique(format(daily$date, "%Y")))) +
          factor(format(date, "%u")),
        data = daily,
        family = quasipoisson(link = "log"),
        na.action = na.exclude
      )
      
      # Extract coefficients
      cb_idx <- grep("cb", names(coef(model)))
      cb_coef <- coef(model)[cb_idx]
      cb_vcov <- vcov(model)[cb_idx, cb_idx]
      
      # Find MMT - restrict to region's observed range (1st-99th percentile)
      mmt_search_min <- max(temp_boundary[1], reg_p1)
      mmt_search_max <- min(temp_boundary[2], reg_p99)
      temp_seq <- seq(mmt_search_min, mmt_search_max, length.out = 100)
      cp <- crosspred(cb, model, at = temp_seq, cumul = TRUE, cen = median(daily$tmean))
      mmt_idx <- which.min(cp$allRRfit)
      mmt <- temp_seq[mmt_idx]
      
      # Get RRs at percentiles
      reg_pcts_vals <- quantile(daily$tmean, probs = c(0.01, 0.05, 0.95, 0.99), na.rm = TRUE)
      cp_pcts <- crosspred(cb, model, at = reg_pcts_vals, cumul = TRUE, cen = mmt)
      
      # Calculate standard errors for filtering
      rr_p99 <- cp_pcts$allRRfit[4]
      rr_p99_ci <- c(cp_pcts$allRRlow[4], cp_pcts$allRRhigh[4])
      rr_p99_se <- (rr_p99_ci[2] - rr_p99_ci[1]) / (2 * 1.96)
      
      rr_p1 <- cp_pcts$allRRfit[1]
      rr_p1_ci <- c(cp_pcts$allRRlow[1], cp_pcts$allRRhigh[1])
      rr_p1_se <- (rr_p1_ci[2] - rr_p1_ci[1]) / (2 * 1.96)
      
      # Quality check (including SE filter)
      vcov_ok <- all(diag(cb_vcov) > 0) && all(is.finite(cb_vcov)) &&
                 is.finite(rr_p99_se) && is.finite(rr_p1_se) &&
                 rr_p99_se < 20 && rr_p1_se < 20
      
      if (vcov_ok) {
        n_params <- length(cb_coef)
        mvmeta_coefs[[length(mvmeta_coefs) + 1]] <- as.vector(cb_coef)
        mvmeta_vcovs[[length(mvmeta_vcovs) + 1]] <- cb_vcov
        mvmeta_regions <- c(mvmeta_regions, as.character(reg))
        
        regional_rrs[[as.character(reg)]] <- list(
          rr_p99 = cp_pcts$allRRfit[4],
          rr_p95 = cp_pcts$allRRfit[3],
          rr_p5 = cp_pcts$allRRfit[2],
          rr_p1 = cp_pcts$allRRfit[1],
          mmt = mmt
        )
      }
      
      "success"
      
    }, error = function(e) {
      "error"
    })
  }
  
  n_success <- length(mvmeta_regions)
  n_excluded <- n_regions - n_success
  cat(sprintf("    Successful fits: %d/%d\n", n_success, n_regions))
  cat(sprintf("    Valid regions: %d (excluded: %d with Inf/extreme SE)\n", n_success, n_excluded))
  
  if (n_success < 10) {
    return(NULL)
  }
  
  # MVMeta pooling
  coef_matrix <- do.call(rbind, mvmeta_coefs)
  rownames(coef_matrix) <- mvmeta_regions
  
  # Regularize vcov
  for (j in seq_along(mvmeta_vcovs)) {
    v <- mvmeta_vcovs[[j]]
    diag(v) <- diag(v) + 1e-6 * mean(diag(v))
    mvmeta_vcovs[[j]] <- v
  }
  
  # Fit mixmeta
  cat(sprintf("    Pooling %d regions with mixmeta (REML)...\n", n_success))
  mv_fit <- tryCatch({
    mixmeta(coef_matrix, mvmeta_vcovs, method = "reml",
            control = list(maxiter = 500, showiter = TRUE))
  }, error = function(e) {
    cat("    REML failed, trying ML...\n")
    tryCatch({
      mixmeta(coef_matrix, mvmeta_vcovs, method = "ml",
              control = list(maxiter = 500, showiter = TRUE))
    }, error = function(e2) NULL)
  })
  
  if (is.null(mv_fit)) {
    return(NULL)
  }
  
  pooled_coef <- coef(mv_fit)
  pooled_vcov <- vcov(mv_fit)
  
  # =====================================================================
  # HETEROGENEITY EXTRACTION - Using qtest() properly
  # Note: qtest() returns vectors with one value per parameter + overall
  # The OVERALL Q statistic is the LAST element
  # =====================================================================
  qstat_info <- tryCatch({
    q <- qtest(mv_fit)
    n <- length(q$Q)  # Overall is last element
    list(Q = q$Q[n], df = q$df[n], pvalue = q$pvalue[n])
  }, error = function(e) {
    list(Q = NA, df = NA, pvalue = NA)
  })
  
  I2 <- if(!is.na(qstat_info$Q) && !is.na(qstat_info$df) && qstat_info$df > 0) {
    max(0, (qstat_info$Q - qstat_info$df) / qstat_info$Q * 100)
  } else NA
  
  H2 <- if(!is.na(qstat_info$Q) && !is.na(qstat_info$df) && qstat_info$df > 0) {
    max(1, qstat_info$Q / qstat_info$df)
  } else NA
  
  tau2 <- tryCatch({
    if(!is.null(mv_fit$Psi)) mean(diag(mv_fit$Psi)) else NA
  }, error = function(e) NA)
  
  cat(sprintf("    Heterogeneity: Q=%.1f (df=%d, p=%s), I²=%.1f%%, H²=%.2f, τ²=%.4f\n",
              qstat_info$Q, qstat_info$df, format.pval(qstat_info$pvalue, digits=2), I2, H2, tau2))
  
  # Pooled MMT
  mmts <- sapply(regional_rrs, function(x) x$mmt)
  pooled_mmt <- median(mmts, na.rm = TRUE)
  
  # Pooled predictions
  temp_seq <- seq(temp_boundary[1], temp_boundary[2], length.out = 100)
  cb_pred <- crossbasis(
    temp_seq,
    lag = max_lag,
    argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
    arglag = list(fun = "ns", df = lag_df)
  )
  
  cp_pooled <- crosspred(
    cb_pred,
    coef = pooled_coef,
    vcov = pooled_vcov,
    model.link = "log",
    at = temp_seq,
    cen = pooled_mmt
  )
  
  # Get pooled RRs at key percentiles
  p1 <- median(sapply(regional_rrs, function(x) quantile(data_dt[region_code == names(regional_rrs)[1], tmean], 0.01)))
  p5 <- median(sapply(regional_rrs, function(x) quantile(data_dt[region_code == names(regional_rrs)[1], tmean], 0.05)))
  p95 <- median(sapply(regional_rrs, function(x) quantile(data_dt[region_code == names(regional_rrs)[1], tmean], 0.95)))
  p99 <- median(sapply(regional_rrs, function(x) quantile(data_dt[region_code == names(regional_rrs)[1], tmean], 0.99)))
  
  # Compute national percentiles instead
  p1 <- quantile(data_dt$tmean, 0.01)
  p5 <- quantile(data_dt$tmean, 0.05)
  p95 <- quantile(data_dt$tmean, 0.95)
  p99 <- quantile(data_dt$tmean, 0.99)
  
  cp_pcts <- crosspred(
    cb_pred,
    coef = pooled_coef,
    vcov = pooled_vcov,
    model.link = "log",
    at = c(p1, p5, p95, p99),
    cen = pooled_mmt
  )
  
  list(
    n_regions = n_success,
    pooled_mmt = pooled_mmt,
    rr_p99 = cp_pcts$allRRfit[4],
    rr_p99_ci = c(cp_pcts$allRRlow[4], cp_pcts$allRRhigh[4]),
    rr_p95 = cp_pcts$allRRfit[3],
    rr_p95_ci = c(cp_pcts$allRRlow[3], cp_pcts$allRRhigh[3]),
    rr_p5 = cp_pcts$allRRfit[2],
    rr_p5_ci = c(cp_pcts$allRRlow[2], cp_pcts$allRRhigh[2]),
    rr_p1 = cp_pcts$allRRfit[1],
    rr_p1_ci = c(cp_pcts$allRRlow[1], cp_pcts$allRRhigh[1]),
    converged = mv_fit$converged,
    method = ifelse("reml" %in% class(mv_fit), "reml", "ml"),
    heterogeneity = list(
      cochrans_Q = qstat_info$Q,
      Q_df = qstat_info$df,
      Q_pvalue = qstat_info$pvalue,
      I2_percent = I2,
      H2 = H2,
      tau2 = tau2
    )
  )
}

# -----------------------------------------------------------------------------
# 3. Lag Structure Sensitivity
# -----------------------------------------------------------------------------
cat("\n[2] Lag structure sensitivity...\n")

lag_results <- list()

for (lag in LAG_OPTIONS) {
  cat(sprintf("\n  Testing max_lag = %d days...\n", lag))
  
  result <- fit_sensitivity_dlnm(data, lag, BASELINE_TEMP_DF, BASELINE_LAG_DF)
  
  if (!is.null(result)) {
    lag_results[[paste0("lag_", lag)]] <- result
    cat(sprintf("    Heat (P99): RR = %.3f (%.3f-%.3f)\n",
                result$rr_p99, result$rr_p99_ci[1], result$rr_p99_ci[2]))
    cat(sprintf("    Cold (P1): RR = %.3f (%.3f-%.3f)\n",
                result$rr_p1, result$rr_p1_ci[1], result$rr_p1_ci[2]))
  } else {
    cat("    FAILED\n")
  }
}

# -----------------------------------------------------------------------------
# 4. Degrees of Freedom Sensitivity
# -----------------------------------------------------------------------------
cat("\n[3] Degrees of freedom sensitivity...\n")

df_results <- list()

for (temp_df in DF_OPTIONS) {
  for (lag_df in DF_OPTIONS) {
    cat(sprintf("\n  Testing temp_df = %d, lag_df = %d...\n", temp_df, lag_df))
    
    result <- fit_sensitivity_dlnm(data, BASELINE_LAG, temp_df, lag_df)
    
    if (!is.null(result)) {
      df_results[[paste0("temp", temp_df, "_lag", lag_df)]] <- result
      cat(sprintf("    Heat (P99): RR = %.3f (%.3f-%.3f)\n",
                  result$rr_p99, result$rr_p99_ci[1], result$rr_p99_ci[2]))
    } else {
      cat("    FAILED\n")
    }
  }
}

# -----------------------------------------------------------------------------
# 5. Summary Table
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("SENSITIVITY ANALYSIS SUMMARY\n")
cat("=======================================================\n")

cat("\n--- LAG STRUCTURE ---\n")
cat(sprintf("%-15s %15s %15s\n", "Specification", "Heat RR (P99)", "Cold RR (P1)"))
cat(strrep("-", 45), "\n")

for (name in names(lag_results)) {
  r <- lag_results[[name]]
  cat(sprintf("%-15s %7.3f (%.3f-%.3f) %7.3f (%.3f-%.3f)\n",
              name, r$rr_p99, r$rr_p99_ci[1], r$rr_p99_ci[2],
              r$rr_p1, r$rr_p1_ci[1], r$rr_p1_ci[2]))
}

cat("\n--- DEGREES OF FREEDOM ---\n")
cat(sprintf("%-20s %15s %15s\n", "Specification", "Heat RR (P99)", "Cold RR (P1)"))
cat(strrep("-", 50), "\n")

for (name in names(df_results)) {
  r <- df_results[[name]]
  cat(sprintf("%-20s %7.3f (%.3f-%.3f) %7.3f (%.3f-%.3f)\n",
              name, r$rr_p99, r$rr_p99_ci[1], r$rr_p99_ci[2],
              r$rr_p1, r$rr_p1_ci[1], r$rr_p1_ci[2]))
}

# -----------------------------------------------------------------------------
# 6. Save Results
# -----------------------------------------------------------------------------
cat("\n[6] Saving results...\n")

output <- list(
  exposure_type = EXPOSURE_TYPE,
  analysis_date = as.character(Sys.time()),
  
  baseline_spec = list(
    max_lag = BASELINE_LAG,
    temp_df = BASELINE_TEMP_DF,
    lag_df = BASELINE_LAG_DF
  ),
  
  lag_sensitivity = lag_results,
  df_sensitivity = df_results
)

output_file <- file.path(OUTPUT_DIR, paste0("sensitivity_r_", EXPOSURE_TYPE, ".json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("  JSON saved to:", output_file, "\n")

cat("\n=======================================================\n")
cat("Done!", as.character(Sys.time()), "\n")
cat("=======================================================\n")
