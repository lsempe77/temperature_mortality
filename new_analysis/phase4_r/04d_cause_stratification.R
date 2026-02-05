# =============================================================================
# 04d_cause_stratification.R
# Cause-Specific DLNM Analysis
# =============================================================================
#
# Tests heterogeneity in temperature-mortality associations by cause:
# 1. Cardiovascular diseases (ICD-10: I00-I99)
# 2. Respiratory diseases (ICD-10: J00-J99)
# 3. External causes (ICD-10: V01-Y98)
# 4. All other causes
#
# **NOTE:** This script requires cause-stratified mortality data.
# Run `00n_cause_stratification_prep.py` first to create:
# - mortality_intermediate_daily_cvd.parquet
# - mortality_intermediate_daily_respiratory.parquet
# - mortality_intermediate_daily_external.parquet
# - mortality_intermediate_daily_other.parquet
#
# References:
# - Basu & Samet (2002) Epidemiology - Cause-specific temperature effects
# - Gasparrini et al. (2012) Epidemiology - Cause-specific burden
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
cat("CAUSE-STRATIFIED DLNM ANALYSIS:", EXPOSURE_TYPE, "\n")
cat("=======================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# Configuration
MAX_LAG <- 21
TEMP_DF <- 4
LAG_DF <- 4

CAUSE_GROUPS <- list(
  "cvd" = "Cardiovascular (I00-I99)",
  "respiratory" = "Respiratory (J00-J99)",
  "external" = "External (V01-Y98)",
  "other" = "All Other Causes"
)

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------
# Get script directory for reliable relative paths
SCRIPT_DIR <- tryCatch({
  dirname(normalizePath(sys.frame(1)$ofile))
}, error = function(e) {
  # Fallback: use command line args if sourced
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("--file=", "", file_arg)))
  } else {
    getwd()
  }
})

DATA_DIR <- file.path(SCRIPT_DIR, "..", "phase0_data_prep", "results")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

cat("Script directory:", SCRIPT_DIR, "\n")
cat("Data directory:", normalizePath(DATA_DIR, mustWork = FALSE), "\n")

cat("[1] Loading data...\n")

# ERA5 temperature data (same for all causes)
era5_file <- file.path(DATA_DIR, paste0("era5_", EXPOSURE_TYPE, "_daily.parquet"))
era5 <- as.data.table(read_parquet(era5_file))
era5[, date := as.character(date)]  # Ensure date is character for merge

# Check if cause-stratified mortality files exist
cause_files <- sapply(names(CAUSE_GROUPS), function(cause) {
  file.path(DATA_DIR, paste0("mortality_", EXPOSURE_TYPE, "_daily_", cause, ".parquet"))
})

missing_files <- !file.exists(cause_files)
if (any(missing_files)) {
  cat("\n*** ERROR: Cause-stratified mortality files not found ***\n")
  cat("Missing files:\n")
  for (f in cause_files[missing_files]) {
    cat("  -", f, "\n")
  }
  cat("\nPlease run cause stratification data prep script first:\n")
  cat("  python phase0_data_prep/aggregation/00n_cause_stratification_prep.py\n\n")
  stop("Missing cause-stratified data files")
}

# Global temperature percentiles for boundaries - use GLOBAL knots for meta-analysis
temp_all <- era5$temp_mean
temp_pcts <- quantile(temp_all, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)
temp_boundary <- c(0, 40)  # Expanded from c(0.7, 35.1) per expert recommendation

cat(sprintf("  GLOBAL knots (P10/P75/P90): %.1f / %.1f / %.1f\n", temp_pcts[1], temp_pcts[2], temp_pcts[3]))
cat(sprintf("  Boundary knots: %.0f to %.0f\n", temp_boundary[1], temp_boundary[2]))

# -----------------------------------------------------------------------------
# 2. Function to Fit Cause-Specific DLNM
# -----------------------------------------------------------------------------
fit_cause_dlnm <- function(mort_data, era5_data, cause_label) {
  
  cat(sprintf("\n[CAUSE GROUP: %s]\n", cause_label))
  
  # Standardize column names
  region_col <- paste0(EXPOSURE_TYPE, "_code")
  if (region_col %in% names(mort_data)) {
    setnames(mort_data, region_col, "region_code")
  }
  
  # Merge
  setkey(mort_data, region_code, date)
  setkey(era5_data, region_code, date)
  data <- merge(mort_data, era5_data, by = c("region_code", "date"))
  
  # Convert date back to Date class for model fitting (needed for format() calls)
  data[, date := as.Date(date)]
  
  # Working variables - try to find deaths column
  if ("deaths" %in% names(data)) {
    # Already named
  } else {
    # Try common patterns
    death_cols <- grep("death", names(data), value = TRUE, ignore.case = TRUE)
    if (length(death_cols) > 0) {
      setnames(data, death_cols[1], "deaths")
    }
  }
  
  data[, tmean := temp_mean]
  data <- data[!is.na(deaths) & !is.na(tmean)]
  setorder(data, region_code, date)
  
  regions <- unique(data$region_code)
  n_regions <- length(regions)
  
  cat(sprintf("  Fitting %d regions...\n", n_regions))
  
  # Storage
  mvmeta_coefs <- list()
  mvmeta_vcovs <- list()
  mvmeta_regions <- c()
  regional_rrs <- list()
  
  for (i in seq_along(regions)) {
    reg <- regions[i]
    
    if (i %% 50 == 0) {
      cat(sprintf("    Progress: %d/%d regions\n", i, n_regions))
    }
    
    daily <- data[region_code == reg][order(date)]
    
    if (nrow(daily) < 365) next
    
    # Convert date from character to Date class for modeling
    daily[, date := as.Date(date)]
    
    result <- tryCatch({
      # Use GLOBAL knots for meta-analysis comparability (expert recommendation)
      reg_p1 <- quantile(daily$tmean, 0.01, na.rm = TRUE)
      reg_p99 <- quantile(daily$tmean, 0.99, na.rm = TRUE)
      
      # Create cross-basis with GLOBAL knots
      cb <- crossbasis(
        daily$tmean,
        lag = MAX_LAG,
        argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
        arglag = list(fun = "ns", df = LAG_DF)
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
      reg_pcts_vals <- quantile(daily$tmean, probs = c(0.01, 0.99), na.rm = TRUE)
      cp_pcts <- crosspred(cb, model, at = reg_pcts_vals, cumul = TRUE, cen = mmt)
      
      # Calculate standard errors for filtering
      rr_p99 <- cp_pcts$allRRfit[2]
      rr_p99_ci <- c(cp_pcts$allRRlow[2], cp_pcts$allRRhigh[2])
      rr_p99_se <- (rr_p99_ci[2] - rr_p99_ci[1]) / (2 * 1.96)
      
      rr_p1 <- cp_pcts$allRRfit[1]
      rr_p1_ci <- c(cp_pcts$allRRlow[1], cp_pcts$allRRhigh[1])
      rr_p1_se <- (rr_p1_ci[2] - rr_p1_ci[1]) / (2 * 1.96)
      
      # Quality check
      vcov_ok <- all(diag(cb_vcov) > 0) && all(is.finite(cb_vcov)) &&
                 is.finite(rr_p99_se) && is.finite(rr_p1_se) &&
                 rr_p99_se < 20 && rr_p1_se < 20
      
      if (vcov_ok) {
        mvmeta_coefs[[length(mvmeta_coefs) + 1]] <- as.vector(cb_coef)
        mvmeta_vcovs[[length(mvmeta_vcovs) + 1]] <- cb_vcov
        mvmeta_regions <- c(mvmeta_regions, as.character(reg))
        
        regional_rrs[[as.character(reg)]] <- list(
          rr_p99 = rr_p99,
          rr_p1 = rr_p1,
          mmt = mmt
        )
      }
      
      "success"
      
    }, error = function(e) {
      "error"
    })
  }
  
  n_success <- length(mvmeta_regions)
  cat(sprintf("    Successful fits: %d/%d (%.1f%%)\n", 
              n_success, n_regions, 100 * n_success / n_regions))
  
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
  cat("    Pooling with mixmeta (REML)...\n")
  mv_fit <- tryCatch({
    mixmeta(coef_matrix, mvmeta_vcovs, method = "reml",
            control = list(maxiter = 500, showiter = FALSE))
  }, error = function(e) {
    cat("    REML failed, trying ML...\n")
    tryCatch({
      mixmeta(coef_matrix, mvmeta_vcovs, method = "ml",
              control = list(maxiter = 500, showiter = FALSE))
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
  
  # Calculate I² and H² from Q statistic
  I2 <- if(!is.na(qstat_info$Q) && !is.na(qstat_info$df) && qstat_info$df > 0) {
    max(0, (qstat_info$Q - qstat_info$df) / qstat_info$Q * 100)
  } else NA
  
  H2 <- if(!is.na(qstat_info$Q) && !is.na(qstat_info$df) && qstat_info$df > 0) {
    max(1, qstat_info$Q / qstat_info$df)
  } else NA
  
  # Extract tau² from model's random effects variance (Psi matrix)
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
    lag = MAX_LAG,
    argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
    arglag = list(fun = "ns", df = LAG_DF)
  )
  
  # National percentiles
  p1 <- quantile(data$tmean, 0.01)
  p99 <- quantile(data$tmean, 0.99)
  
  cp_pcts <- crosspred(
    cb_pred,
    coef = pooled_coef,
    vcov = pooled_vcov,
    model.link = "log",
    at = c(p1, p99),
    cen = pooled_mmt
  )
  
  list(
    n_regions = n_success,
    pooled_mmt = pooled_mmt,
    heterogeneity = list(
      cochrans_Q = qstat_info$Q,
      Q_df = qstat_info$df,
      Q_pvalue = qstat_info$pvalue,
      I2_percent = I2
    ),
    heat = list(
      rr = cp_pcts$allRRfit[2],
      ci_low = cp_pcts$allRRlow[2],
      ci_high = cp_pcts$allRRhigh[2],
      log_rr = log(cp_pcts$allRRfit[2]),
      se = (log(cp_pcts$allRRhigh[2]) - log(cp_pcts$allRRlow[2])) / (2 * 1.96)
    ),
    cold = list(
      rr = cp_pcts$allRRfit[1],
      ci_low = cp_pcts$allRRlow[1],
      ci_high = cp_pcts$allRRhigh[1],
      log_rr = log(cp_pcts$allRRfit[1]),
      se = (log(cp_pcts$allRRhigh[1]) - log(cp_pcts$allRRlow[1])) / (2 * 1.96)
    ),
    converged = mv_fit$converged,
    method = ifelse("reml" %in% class(mv_fit), "reml", "ml")
  )
}

# -----------------------------------------------------------------------------
# 3. Fit All Cause Groups
# -----------------------------------------------------------------------------
cat("\n[2] Fitting cause-stratified DLNMs...\n")

cause_results <- list()

for (cause_code in names(CAUSE_GROUPS)) {
  cause_label <- CAUSE_GROUPS[[cause_code]]
  
  # Load cause-specific mortality
  mort_file <- file.path(DATA_DIR, paste0("mortality_", EXPOSURE_TYPE, "_daily_", cause_code, ".parquet"))
  mort <- as.data.table(read_parquet(mort_file))
  mort[, date := as.character(as.Date(date))]  # Convert POSIXct to Date then to character for merge
  
  # Fit DLNM
  result <- fit_cause_dlnm(mort, era5, cause_label)
  
  if (!is.null(result)) {
    cause_results[[cause_code]] <- result
    cat(sprintf("    %s - Heat: RR=%.3f (%.3f-%.3f), Cold: RR=%.3f (%.3f-%.3f)\n",
                cause_label,
                result$heat$rr, result$heat$ci_low, result$heat$ci_high,
                result$cold$rr, result$cold$ci_low, result$cold$ci_high))
  } else {
    cat(sprintf("    %s - FAILED\n", cause_label))
  }
}

# -----------------------------------------------------------------------------
# 4. Test for Cause Differences
# -----------------------------------------------------------------------------
cat("\n[3] Testing for cause differences...\n")

if (length(cause_results) >= 2) {
  # Compare CVD vs Respiratory (if both available)
  if ("cvd" %in% names(cause_results) && "respiratory" %in% names(cause_results)) {
    heat_diff <- cause_results$cvd$heat$log_rr - cause_results$respiratory$heat$log_rr
    heat_diff_se <- sqrt(cause_results$cvd$heat$se^2 + cause_results$respiratory$heat$se^2)
    heat_z <- heat_diff / heat_diff_se
    heat_p <- 2 * (1 - pnorm(abs(heat_z)))
    
    cold_diff <- cause_results$cvd$cold$log_rr - cause_results$respiratory$cold$log_rr
    cold_diff_se <- sqrt(cause_results$cvd$cold$se^2 + cause_results$respiratory$cold$se^2)
    cold_z <- cold_diff / cold_diff_se
    cold_p <- 2 * (1 - pnorm(abs(cold_z)))
    
    cat(sprintf("\nCardiovascular vs Respiratory:\n"))
    cat(sprintf("  Heat: diff=%.4f, SE=%.4f, Z=%.2f, p=%.4f %s\n", 
                heat_diff, heat_diff_se, heat_z, heat_p,
                ifelse(heat_p < 0.05, "*", "ns")))
    cat(sprintf("  Cold: diff=%.4f, SE=%.4f, Z=%.2f, p=%.4f %s\n",
                cold_diff, cold_diff_se, cold_z, cold_p,
                ifelse(cold_p < 0.05, "*", "ns")))
  }
}

# -----------------------------------------------------------------------------
# 5. Summary
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("CAUSE-STRATIFIED ANALYSIS SUMMARY\n")
cat("=======================================================\n")

cat("\nCause of Death                Heat RR (99th %ile)     Cold RR (1st %ile)\n")
cat("--------------------------------------------------------------------------------\n")

for (cause_code in names(cause_results)) {
  r <- cause_results[[cause_code]]
  cause_label <- CAUSE_GROUPS[[cause_code]]
  cat(sprintf("%-28s  %.3f (%.3f-%.3f)     %.3f (%.3f-%.3f)\n",
              cause_label,
              r$heat$rr, r$heat$ci_low, r$heat$ci_high,
              r$cold$rr, r$cold$ci_low, r$cold$ci_high))
}

cat("\nInterpretation:\n")
cat("- CVD and respiratory expected to show strongest temperature associations\n")
cat("- External causes expected to be less temperature-related\n")

# -----------------------------------------------------------------------------
# 6. Save Results
# -----------------------------------------------------------------------------
cat("\n[6] Saving results...\n")

output <- list(
  exposure_type = EXPOSURE_TYPE,
  analysis_date = as.character(Sys.time()),
  
  cause_groups = cause_results
)

output_file <- file.path(OUTPUT_DIR, paste0("cause_stratification_", EXPOSURE_TYPE, ".json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("  JSON saved to:", output_file, "\n")

cat("\n=======================================================\n")
cat("Done!", as.character(Sys.time()), "\n")
cat("=======================================================\n")
