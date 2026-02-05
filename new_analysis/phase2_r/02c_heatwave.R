# =============================================================================
# 02c_heatwave.R
# Heatwave Effect Modification Analysis
# =============================================================================
#
# Tests whether sustained extreme temperatures (heatwaves) have ADDITIONAL
# mortality effects beyond the daily temperature-mortality relationship.
#
# Key Question: Does a day of 35°C during a multi-day heatwave have stronger
#               effects than an isolated 35°C day?
#
# Method (Gasparrini & Armstrong 2011):
# 1. Define heatwaves (3+ consecutive days above P95)
# 2. Fit DLNM with heatwave indicator: deaths ~ cb(temp) + heatwave + cb×heatwave
# 3. Compare RR during heatwave days vs non-heatwave days
# 4. Pool interaction effects via mixmeta
#
# References:
# - Gasparrini & Armstrong (2011) Stat Med
# - Anderson & Bell (2009) Epidemiology
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
cat("HEATWAVE ANALYSIS:", EXPOSURE_TYPE, "\n")
cat("=======================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# Configuration
MAX_LAG <- 21
TEMP_DF <- 4
LAG_DF <- 4

# Heatwave definition
HEATWAVE_THRESHOLD_PCT <- 0.95  # P95
MIN_HEATWAVE_DAYS <- 3          # 3+ consecutive days

cat("Configuration:\n")
cat("  Max lag:", MAX_LAG, "days\n")
cat("  Heatwave: >=", MIN_HEATWAVE_DAYS, "consecutive days above P", 
    HEATWAVE_THRESHOLD_PCT * 100, "\n\n")

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

# Global temperature percentiles
temp_all <- data$tmean
temp_pcts <- quantile(temp_all, probs = c(0.10, 0.75, 0.90, 0.95, 0.99), na.rm = TRUE)
temp_pcts_knots <- quantile(temp_all, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)  # For cross-basis
temp_boundary <- c(0, 40)  # Expanded from c(0.7, 35.1) per expert recommendation

cat("  Regions:", uniqueN(data$region_code), "\n")
cat("  Date range:", as.character(min(data$date)), "to", as.character(max(data$date)), "\n")
cat("  GLOBAL knots (P10/P75/P90):", round(temp_pcts_knots, 1), "\n")
cat("  Boundary knots:", temp_boundary, "\n")
cat("  Temperature P95 (heatwave threshold):", round(temp_pcts["95%"], 2), "°C\n")

# -----------------------------------------------------------------------------
# 2. Identify Heatwaves
# -----------------------------------------------------------------------------
cat("\n[2] Identifying heatwaves...\n")

identify_heatwaves <- function(temp_vector, threshold, min_days) {
  # Identify consecutive days above threshold
  above_threshold <- temp_vector > threshold
  
  # Run-length encoding to find consecutive sequences
  rle_result <- rle(above_threshold)
  
  # Mark sequences of TRUE that are >= min_days
  heatwave_indicator <- rep(FALSE, length(temp_vector))
  
  end_idx <- cumsum(rle_result$lengths)
  start_idx <- c(1, end_idx[-length(end_idx)] + 1)
  
  for (i in seq_along(rle_result$values)) {
    if (rle_result$values[i] && rle_result$lengths[i] >= min_days) {
      heatwave_indicator[start_idx[i]:end_idx[i]] <- TRUE
    }
  }
  
  return(as.integer(heatwave_indicator))
}

# Add heatwave indicator for each region
data[, heatwave := identify_heatwaves(tmean, temp_pcts["95%"], MIN_HEATWAVE_DAYS), 
     by = region_code]

n_heatwave_days <- sum(data$heatwave)
pct_heatwave <- 100 * n_heatwave_days / nrow(data)

cat(sprintf("  Heatwave days: %d (%.2f%% of total)\n", n_heatwave_days, pct_heatwave))
cat(sprintf("  Regions with heatwaves: %d/%d\n", 
            uniqueN(data[heatwave == 1, region_code]),
            uniqueN(data$region_code)))

# -----------------------------------------------------------------------------
# 3. Fit DLNM with Heatwave Interaction
# -----------------------------------------------------------------------------
cat("\n[3] Fitting DLNM with heatwave interaction...\n")

fit_heatwave_dlnm <- function(data_dt, max_lag, temp_df, lag_df) {
  
  regions <- unique(data_dt$region_code)
  n_regions <- length(regions)
  
  cat(sprintf("  Fitting %d regions with heatwave interaction...\n", n_regions))
  
  # Store results
  main_rr_p99 <- list()
  interaction_coefs <- list()
  interaction_vcovs <- list()
  valid_regions <- c()
  
  n_success <- 0
  n_no_heatwave <- 0
  n_errors <- 0
  
  for (i in seq_along(regions)) {
    reg <- regions[i]
    
    if (i %% 50 == 0) {
      cat(sprintf("    Progress: %d/%d regions\n", i, n_regions))
    }
    
    # Subset data
    daily <- data_dt[region_code == reg][order(date)]
    n_days <- nrow(daily)
    
    if (n_days < 365) next
    
    # Check if region has any heatwave days
    if (sum(daily$heatwave) < 10) {
      n_no_heatwave <- n_no_heatwave + 1
      next
    }
    
    result <- tryCatch({
      # Use GLOBAL knots for meta-analysis comparability (expert recommendation)
      reg_p1 <- quantile(daily$tmean, 0.01, na.rm = TRUE)
      reg_p99 <- quantile(daily$tmean, 0.99, na.rm = TRUE)
      reg_temp_pcts <- quantile(daily$tmean, probs = c(0.10, 0.75, 0.90, 0.99), 
                                na.rm = TRUE)
      
      # Create cross-basis with GLOBAL knots
      cb <- crossbasis(
        daily$tmean,
        lag = max_lag,
        argvar = list(fun = "ns", knots = temp_pcts_knots, 
                     Boundary.knots = temp_boundary),
        arglag = list(fun = "ns", df = lag_df)
      )
      
      # Fit GLM with heatwave main effect (no interaction for now - simpler model)
      model <- glm(
        deaths ~ cb + heatwave + 
          ns(as.numeric(date), df = 7 * length(unique(format(daily$date, "%Y")))) +
          factor(format(date, "%u")),
        data = daily,
        family = quasipoisson(link = "log"),
        na.action = na.exclude
      )
      
      # Extract coefficients
      cb_idx <- grep("^cb", names(coef(model)))
      heatwave_idx <- which(names(coef(model)) == "heatwave")
      
      if (length(heatwave_idx) == 0) {
        return("no_heatwave_coef")
      }
      
      heatwave_coef <- coef(model)[heatwave_idx]
      heatwave_se <- sqrt(diag(vcov(model))[heatwave_idx])
      
      # Get main temperature effect at P99
      temp_p99 <- reg_temp_pcts["99%"]
      
      # Find MMT - restrict to region's observed range (1st-99th percentile)
      mmt_search_min <- max(temp_boundary[1], reg_p1)
      mmt_search_max <- min(temp_boundary[2], reg_p99)
      temp_seq <- seq(mmt_search_min, mmt_search_max, length.out = 100)
      cp_full <- crosspred(cb, model, at = temp_seq, cumul = TRUE, cen = median(daily$tmean))
      mmt_idx <- which.min(cp_full$allRRfit)
      mmt <- temp_seq[mmt_idx]
      
      # Get RR at P99
      cp_p99 <- crosspred(cb, model, at = temp_p99, cumul = TRUE, cen = mmt)
      rr_p99 <- cp_p99$allRRfit[1]
      
      # Store results
      if (is.finite(rr_p99) && is.finite(heatwave_coef) && is.finite(heatwave_se)) {
        main_rr_p99[[length(main_rr_p99) + 1]] <- rr_p99
        interaction_coefs[[length(interaction_coefs) + 1]] <- heatwave_coef
        interaction_vcovs[[length(interaction_vcovs) + 1]] <- heatwave_se^2
        valid_regions <- c(valid_regions, as.character(reg))
        n_success <- n_success + 1
      }
      
      "success"
      
    }, error = function(e) {
      n_errors <- n_errors + 1
      "error"
    })
  }
  
  cat(sprintf("  Successful fits: %d/%d\n", n_success, n_regions))
  cat(sprintf("  Regions without heatwaves: %d\n", n_no_heatwave))
  cat(sprintf("  Errors: %d\n", n_errors))
  
  if (n_success < 10) {
    return(NULL)
  }
  
  # Pool heatwave main effect
  cat("\n  Pooling heatwave effects via mixmeta...\n")
  
  hw_coefs <- unlist(interaction_coefs)
  hw_vcovs <- unlist(interaction_vcovs)
  
  mv_hw <- tryCatch({
    mixmeta(hw_coefs, hw_vcovs, method = "reml",
            control = list(maxiter = 500, showiter = FALSE))
  }, error = function(e) {
    tryCatch({
      mixmeta(hw_coefs, hw_vcovs, method = "ml",
              control = list(maxiter = 500, showiter = FALSE))
    }, error = function(e2) NULL)
  })
  
  if (!is.null(mv_hw)) {
    pooled_log_rr_hw <- coef(mv_hw)[1]
    pooled_se_hw <- sqrt(vcov(mv_hw)[1, 1])
    pooled_rr_hw <- exp(pooled_log_rr_hw)
    pooled_ci_hw <- exp(pooled_log_rr_hw + c(-1.96, 1.96) * pooled_se_hw)
    
    # Heterogeneity extraction using qtest() properly
    qstat_hw <- tryCatch({
      q <- qtest(mv_hw)
      n <- length(q$Q)  # Overall is last element
      list(Q = q$Q[n], df = q$df[n], pvalue = q$pvalue[n])
    }, error = function(e) list(Q = NA, df = NA, pvalue = NA))
    
    I2_hw <- if(!is.na(qstat_hw$Q) && !is.na(qstat_hw$df) && qstat_hw$df > 0) {
      max(0, (qstat_hw$Q - qstat_hw$df) / qstat_hw$Q * 100)
    } else NA
    
    cat(sprintf("    Heatwave RR (additive): %.3f (%.3f-%.3f), I²=%.1f%%\n", 
                pooled_rr_hw, pooled_ci_hw[1], pooled_ci_hw[2], I2_hw))
  } else {
    pooled_rr_hw <- NA
    pooled_ci_hw <- c(NA, NA)
    pooled_log_rr_hw <- NA
    pooled_se_hw <- NA
    qstat_hw <- list(Q = NA, df = NA, pvalue = NA)
    I2_hw <- NA
  }
  
  # Average main temperature effect
  avg_rr_p99 <- median(unlist(main_rr_p99), na.rm = TRUE)
  
  list(
    n_regions = n_success,
    n_heatwave_days = n_heatwave_days,
    pct_heatwave = pct_heatwave,
    avg_main_rr_p99 = avg_rr_p99,
    heatwave_rr = pooled_rr_hw,
    heatwave_ci_low = pooled_ci_hw[1],
    heatwave_ci_high = pooled_ci_hw[2],
    heatwave_log_rr = pooled_log_rr_hw,
    heatwave_se = pooled_se_hw,
    heterogeneity = list(
      cochrans_Q = qstat_hw$Q,
      Q_df = qstat_hw$df,
      Q_pvalue = qstat_hw$pvalue,
      I2_percent = I2_hw
    )
  )
}

heatwave_results <- fit_heatwave_dlnm(data, MAX_LAG, TEMP_DF, LAG_DF)

# -----------------------------------------------------------------------------
# 4. Summary
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("HEATWAVE ANALYSIS SUMMARY\n")
cat("=======================================================\n")

if (!is.null(heatwave_results)) {
  
  cat("\n--- HEATWAVE EFFECT ---\n")
  cat(sprintf("Definition: %d+ consecutive days above P95 (%.1f°C)\n",
              MIN_HEATWAVE_DAYS, temp_pcts["95%"]))
  cat(sprintf("Heatwave days: %d (%.2f%% of total)\n",
              heatwave_results$n_heatwave_days, heatwave_results$pct_heatwave))
  cat(sprintf("Regions analyzed: %d\n", heatwave_results$n_regions))
  
  cat("\n--- MORTALITY EFFECTS ---\n")
  cat(sprintf("Main temperature effect (P99): RR = %.3f\n", 
              heatwave_results$avg_main_rr_p99))
  cat(sprintf("Heatwave additive effect: RR = %.3f (%.3f-%.3f)\n",
              heatwave_results$heatwave_rr,
              heatwave_results$heatwave_ci_low,
              heatwave_results$heatwave_ci_high))
  
  cat("\nInterpretation:\n")
  if (heatwave_results$heatwave_rr > 1.0 && heatwave_results$heatwave_ci_low > 1.0) {
    cat("  ✓ SIGNIFICANT heatwave effect: Multi-day heat has additional mortality risk\n")
    pct_increase <- (heatwave_results$heatwave_rr - 1) * 100
    cat(sprintf("  ✓ Heatwave days have %.1f%% higher mortality than expected from temperature alone\n", 
                pct_increase))
  } else if (heatwave_results$heatwave_rr > 1.0) {
    cat("  ~ POSITIVE but non-significant heatwave effect\n")
  } else {
    cat("  ✗ NO additional heatwave effect detected\n")
    cat("  → Mortality explained by daily temperature alone\n")
  }
}

# -----------------------------------------------------------------------------
# 5. Save Results
# -----------------------------------------------------------------------------
cat("\n[5] Saving results...\n")

output <- list(
  exposure_type = EXPOSURE_TYPE,
  analysis_date = as.character(Sys.time()),
  heatwave_definition = list(
    threshold_pct = HEATWAVE_THRESHOLD_PCT,
    threshold_temp = temp_pcts["95%"],
    min_days = MIN_HEATWAVE_DAYS
  ),
  results = heatwave_results
)

output_file <- file.path(OUTPUT_DIR, paste0("heatwave_r_", EXPOSURE_TYPE, ".json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("  JSON saved to:", output_file, "\n")

cat("\n=======================================================\n")
cat("Done!", as.character(Sys.time()), "\n")
cat("=======================================================\n")
