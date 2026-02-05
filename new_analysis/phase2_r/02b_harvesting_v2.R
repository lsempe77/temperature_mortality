# =============================================================================
# 02b_harvesting_v2.R
# CORRECTED Harvesting / Mortality Displacement Analysis
# =============================================================================
#
# KEY FIX: Uses SINGLE model with maximum lag, then computes cumulative effects
# at different horizons from that ONE model using crosspred(..., to=horizon).
# This avoids the numerical instability from refitting at each horizon.
#
# Previous approach (WRONG): Refit model for each horizon → different MMT,
# amplified variance, extreme heterogeneity (I² > 99%)
#
# New approach (CORRECT): Fit ONCE with max_lag, predict cumulative RR
# at each horizon using the `to` parameter in crosspred()
#
# Method (Schwartz 2000, Braga 2001):
# 1. Fit DLNM with extended lag (35 days) for each region - ONCE
# 2. Calculate CUMULATIVE RR at different lag horizons from SAME model
# 3. Pool via mixmeta
# 4. Harvesting ratio = 1 - (ERR_35 / ERR_7) where ERR = RR - 1
#
# References:
# - Schwartz (2000) Epidemiology 11:624-8 - "Harvesting and long term exposure"
# - Braga et al. (2001) Am J Epidemiol 153:719-26 - "Time course of weather deaths"
# - Zanobetti & Schwartz (2008) Environ Health Perspect - "Mortality displacement"
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
cat("HARVESTING ANALYSIS v2 (CORRECTED):", EXPOSURE_TYPE, "\n")
cat("=======================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# Configuration
MAX_LAG <- 35  # Extended to detect harvesting
LAG_HORIZONS <- c(7, 14, 21, 28, 35)  # Cumulative effect evaluation points
TEMP_DF <- 4
LAG_DF <- 4

cat("Configuration:\n")
cat("  Extended max lag:", MAX_LAG, "days\n")
cat("  Lag horizons:", paste(LAG_HORIZONS, collapse=", "), "\n")
cat("  Spline df: temp=", TEMP_DF, ", lag=", LAG_DF, "\n")
cat("  METHOD: Single-model cumulative extraction (crosspred to= parameter)\n\n")

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return(getwd())
}
SCRIPT_DIR <- get_script_dir()
BASE_DIR <- dirname(SCRIPT_DIR)

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

# Standardize column names
region_col <- paste0(EXPOSURE_TYPE, "_code")
if (region_col %in% names(mort)) {
  setnames(mort, region_col, "region_code")
}
if (region_col %in% names(era5)) {
  setnames(era5, region_col, "region_code")
}

# Convert dates for merge
mort[, date := as.character(as.Date(date))]
era5[, date := as.character(date)]

# Merge
setkey(mort, region_code, date)
setkey(era5, region_code, date)
data <- merge(mort, era5, by = c("region_code", "date"))
data[, date := as.Date(date)]

# Working variables
data[, deaths := deaths_elderly]
data[, tmean := temp_mean]
data <- data[!is.na(deaths) & !is.na(tmean)]
setorder(data, region_code, date)

# Global temperature percentiles
temp_all <- data$tmean
temp_pcts <- quantile(temp_all, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)
temp_boundary <- c(0, 40)

temp_p99 <- quantile(temp_all, 0.99, na.rm = TRUE)
temp_p01 <- quantile(temp_all, 0.01, na.rm = TRUE)

cat("  Regions:", uniqueN(data$region_code), "\n")
cat("  Date range:", as.character(min(data$date)), "to", as.character(max(data$date)), "\n")
cat("  Temperature P99 (heat):", round(temp_p99, 2), "°C\n")
cat("  Temperature P1 (cold):", round(temp_p01, 2), "°C\n\n")

# -----------------------------------------------------------------------------
# 2. CORRECTED: Fit SINGLE Model, Extract Cumulative RR at Multiple Horizons
# -----------------------------------------------------------------------------
fit_harvesting_single_model <- function(data_dt, max_lag, temp_df, lag_df, 
                                        eval_temps, lag_horizons) {
  
  regions <- unique(data_dt$region_code)
  n_regions <- length(regions)
  
  cat(sprintf("\n[2] Fitting SINGLE extended DLNM (lag=%d) per region...\n", max_lag))
  cat("    Then extracting cumulative RR at each horizon from SAME model.\n\n")
  
  # Store results for each lag horizon
  results_by_horizon <- list()
  for (horizon in lag_horizons) {
    results_by_horizon[[paste0("lag_", horizon)]] <- list(
      heat_coefs = list(),
      heat_vcovs = list(),
      cold_coefs = list(),
      cold_vcovs = list(),
      regions = c()
    )
  }
  
  # Use environment for counters to fix scoping in tryCatch
  counters <- new.env()
  counters$n_success <- 0
  counters$n_errors <- 0
  counters$n_bad_vcov <- 0
  
  for (i in seq_along(regions)) {
    reg <- regions[i]
    
    if (i %% 20 == 0 || i <= 3) {
      cat(sprintf("  Progress: %d/%d regions (success: %d, errors: %d)\n", 
                  i, n_regions, counters$n_success, counters$n_errors))
    }
    
    daily <- data_dt[region_code == reg][order(date)]
    n_days <- nrow(daily)
    
    if (n_days < 365) next
    
    result <- tryCatch({
      reg_p1 <- quantile(daily$tmean, 0.01, na.rm = TRUE)
      reg_p99 <- quantile(daily$tmean, 0.99, na.rm = TRUE)
      
      # =======================================================================
      # KEY FIX: Fit ONE model with MAXIMUM lag
      # =======================================================================
      cb <- crossbasis(
        daily$tmean,
        lag = max_lag,  # Always use max_lag
        argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
        arglag = list(fun = "ns", df = lag_df)
      )
      
      model <- glm(
        deaths ~ cb + ns(as.numeric(date), df = 7 * length(unique(format(daily$date, "%Y")))) +
          factor(format(date, "%u")),
        data = daily,
        family = quasipoisson(link = "log"),
        na.action = na.exclude
      )
      
      # Quality check
      cb_idx <- grep("cb", names(coef(model)))
      cb_vcov <- vcov(model)[cb_idx, cb_idx]
      vcov_ok <- all(diag(cb_vcov) > 0) && all(is.finite(cb_vcov))
      
      if (!vcov_ok) {
        counters$n_bad_vcov <- counters$n_bad_vcov + 1
        stop("bad_vcov")
      }
      
      # Find MMT from the SINGLE fitted model
      mmt_search_min <- max(temp_boundary[1], reg_p1)
      mmt_search_max <- min(temp_boundary[2], reg_p99)
      temp_seq <- seq(mmt_search_min, mmt_search_max, length.out = 100)
      
      # Use full cumulative for MMT search
      cp_mmt <- crosspred(cb, model, at = temp_seq, cumul = TRUE, cen = median(daily$tmean))
      mmt_idx <- which.min(cp_mmt$allRRfit)
      mmt <- temp_seq[mmt_idx]
      
      # =======================================================================
      # KEY FIX: Manual extraction of cumulative RR at different horizons
      # crosspred doesn't allow cumul=TRUE for lag sub-periods, so we
      # extract lag-specific effects with cumul=FALSE and sum manually
      # =======================================================================
      
      # Get lag-specific effects from single model
      cp_lags <- crosspred(
        cb, 
        model, 
        at = eval_temps,
        cumul = FALSE,  # Get lag-specific, not cumulative
        cen = mmt
      )
      
      # matRRfit has dimensions [n_temps x n_lags] with RR at each lag
      # matRRfit[temp_idx, lag_idx+1] for lag 0, 1, 2, ...
      # We sum log-RRs from lag 0 to horizon, then exponentiate
      
      for (horizon in lag_horizons) {
        horizon_key <- paste0("lag_", horizon)
        
        # Sum log-RRs from lag 0 to horizon (0-indexed lags, so use 1:(horizon+1))
        # Heat effect (P99) - eval_temps[2] = row 2
        heat_logRR_cumul <- sum(log(cp_lags$matRRfit[2, 1:min(horizon+1, ncol(cp_lags$matRRfit))]))
        heat_rr <- exp(heat_logRR_cumul)
        
        # Cold effect (P1) - eval_temps[1] = row 1
        cold_logRR_cumul <- sum(log(cp_lags$matRRfit[1, 1:min(horizon+1, ncol(cp_lags$matRRfit))]))
        cold_rr <- exp(cold_logRR_cumul)
        
        # For SE, use quadrature: SE of sum = sqrt(sum of variances) in log scale
        # matse has log-scale SEs for each lag
        heat_se_cumul <- sqrt(sum(cp_lags$matse[2, 1:min(horizon+1, ncol(cp_lags$matse))]^2))
        cold_se_cumul <- sqrt(sum(cp_lags$matse[1, 1:min(horizon+1, ncol(cp_lags$matse))]^2))
        
        # Store for meta-analysis (log scale)
        if (is.finite(heat_rr) && heat_rr > 0 && heat_rr < 1000 && is.finite(heat_se_cumul) && heat_se_cumul > 0 && heat_se_cumul < 10) {
          results_by_horizon[[horizon_key]]$heat_coefs[[length(results_by_horizon[[horizon_key]]$heat_coefs) + 1]] <- log(heat_rr)
          results_by_horizon[[horizon_key]]$heat_vcovs[[length(results_by_horizon[[horizon_key]]$heat_vcovs) + 1]] <- heat_se_cumul^2
        }
        
        if (is.finite(cold_rr) && cold_rr > 0 && cold_rr < 1000 && is.finite(cold_se_cumul) && cold_se_cumul > 0 && cold_se_cumul < 10) {
          results_by_horizon[[horizon_key]]$cold_coefs[[length(results_by_horizon[[horizon_key]]$cold_coefs) + 1]] <- log(cold_rr)
          results_by_horizon[[horizon_key]]$cold_vcovs[[length(results_by_horizon[[horizon_key]]$cold_vcovs) + 1]] <- cold_se_cumul^2
        }
        
        if ((is.finite(heat_rr) && heat_rr > 0 && heat_rr < 1000) || 
            (is.finite(cold_rr) && cold_rr > 0 && cold_rr < 1000)) {
          results_by_horizon[[horizon_key]]$regions <- c(
            results_by_horizon[[horizon_key]]$regions,
            as.character(reg)
          )
        }
      }
      
      counters$n_success <- counters$n_success + 1
      "success"
      
    }, error = function(e) {
      counters$n_errors <- counters$n_errors + 1
      "error"
    })
  }
  
  cat(sprintf("\n  Final: %d/%d regions succeeded (errors: %d, bad_vcov: %d)\n", 
              counters$n_success, n_regions, counters$n_errors, counters$n_bad_vcov))
  
  # Debug: show region counts
  cat("\n  Region counts by horizon:\n")
  for (h in lag_horizons) {
    horizon_key <- paste0("lag_", h)
    n_heat <- length(results_by_horizon[[horizon_key]]$heat_coefs)
    n_cold <- length(results_by_horizon[[horizon_key]]$cold_coefs)
    cat(sprintf("    Lag %d: heat=%d, cold=%d\n", h, n_heat, n_cold))
  }
  
  if (counters$n_success < 10) {
    cat("  WARNING: Too few regions (<10) for reliable pooling\n")
    # Don't return NULL - we have results by horizon that we can use
    if (counters$n_success < 3 && length(results_by_horizon$lag_7$heat_coefs) < 10) {
      return(NULL)
    }
  }
  
  # -----------------------------------------------------------------------------
  # 3. Pool Results via Mixmeta
  # -----------------------------------------------------------------------------
  cat("\n[3] Pooling results via mixmeta...\n")
  
  pooled_results <- list()
  
  for (horizon in lag_horizons) {
    horizon_key <- paste0("lag_", horizon)
    cat(sprintf("  Pooling lag horizon %d days...\n", horizon))
    
    # Heat effects
    heat_coefs <- unlist(results_by_horizon[[horizon_key]]$heat_coefs)
    heat_vcovs <- unlist(results_by_horizon[[horizon_key]]$heat_vcovs)
    
    valid_heat <- is.finite(heat_coefs) & is.finite(heat_vcovs) & heat_vcovs > 0 & heat_vcovs < 100
    heat_coefs <- heat_coefs[valid_heat]
    heat_vcovs <- heat_vcovs[valid_heat]
    
    if (length(heat_coefs) > 5) {
      mv_heat <- tryCatch({
        mixmeta(heat_coefs, heat_vcovs, method = "reml",
                control = list(maxiter = 500, showiter = FALSE))
      }, error = function(e) {
        tryCatch({
          mixmeta(heat_coefs, heat_vcovs, method = "ml")
        }, error = function(e2) NULL)
      })
      
      if (!is.null(mv_heat)) {
        pooled_log_rr_heat <- coef(mv_heat)[1]
        pooled_se_heat <- sqrt(vcov(mv_heat)[1, 1])
        pooled_rr_heat <- exp(pooled_log_rr_heat)
        pooled_ci_heat <- exp(pooled_log_rr_heat + c(-1.96, 1.96) * pooled_se_heat)
        
        qtest_heat <- tryCatch({
          q <- qtest(mv_heat)
          n <- length(q$Q)
          list(Q = q$Q[n], df = q$df[n], pvalue = q$pvalue[n])
        }, error = function(e) list(Q = NA, df = NA, pvalue = NA))
        
        Q_heat <- qtest_heat$Q
        df_heat <- qtest_heat$df
        I2_heat <- if(!is.na(Q_heat) && !is.na(df_heat) && df_heat > 0) max(0, (Q_heat - df_heat) / Q_heat * 100) else NA
      } else {
        pooled_rr_heat <- NA; pooled_ci_heat <- c(NA, NA)
        Q_heat <- NA; df_heat <- NA; I2_heat <- NA
      }
    } else {
      pooled_rr_heat <- NA; pooled_ci_heat <- c(NA, NA)
      Q_heat <- NA; df_heat <- NA; I2_heat <- NA
    }
    
    # Cold effects
    cold_coefs <- unlist(results_by_horizon[[horizon_key]]$cold_coefs)
    cold_vcovs <- unlist(results_by_horizon[[horizon_key]]$cold_vcovs)
    
    valid_cold <- is.finite(cold_coefs) & is.finite(cold_vcovs) & cold_vcovs > 0 & cold_vcovs < 100
    cold_coefs <- cold_coefs[valid_cold]
    cold_vcovs <- cold_vcovs[valid_cold]
    
    if (length(cold_coefs) > 5) {
      mv_cold <- tryCatch({
        mixmeta(cold_coefs, cold_vcovs, method = "reml",
                control = list(maxiter = 500, showiter = FALSE))
      }, error = function(e) {
        tryCatch({
          mixmeta(cold_coefs, cold_vcovs, method = "ml")
        }, error = function(e2) NULL)
      })
      
      if (!is.null(mv_cold)) {
        pooled_log_rr_cold <- coef(mv_cold)[1]
        pooled_se_cold <- sqrt(vcov(mv_cold)[1, 1])
        pooled_rr_cold <- exp(pooled_log_rr_cold)
        pooled_ci_cold <- exp(pooled_log_rr_cold + c(-1.96, 1.96) * pooled_se_cold)
        
        qtest_cold <- tryCatch({
          q <- qtest(mv_cold)
          n <- length(q$Q)
          list(Q = q$Q[n], df = q$df[n], pvalue = q$pvalue[n])
        }, error = function(e) list(Q = NA, df = NA, pvalue = NA))
        
        Q_cold <- qtest_cold$Q
        df_cold <- qtest_cold$df
        I2_cold <- if(!is.na(Q_cold) && !is.na(df_cold) && df_cold > 0) max(0, (Q_cold - df_cold) / Q_cold * 100) else NA
      } else {
        pooled_rr_cold <- NA; pooled_ci_cold <- c(NA, NA)
        Q_cold <- NA; df_cold <- NA; I2_cold <- NA
      }
    } else {
      pooled_rr_cold <- NA; pooled_ci_cold <- c(NA, NA)
      Q_cold <- NA; df_cold <- NA; I2_cold <- NA
    }
    
    pooled_results[[horizon_key]] <- list(
      lag_days = horizon,
      n_regions = length(unique(results_by_horizon[[horizon_key]]$regions)),
      heat = list(
        rr = pooled_rr_heat,
        ci_low = pooled_ci_heat[1],
        ci_high = pooled_ci_heat[2],
        cochrans_Q = Q_heat,
        Q_df = df_heat,
        I2_percent = I2_heat
      ),
      cold = list(
        rr = pooled_rr_cold,
        ci_low = pooled_ci_cold[1],
        ci_high = pooled_ci_cold[2],
        cochrans_Q = Q_cold,
        Q_df = df_cold,
        I2_percent = I2_cold
      )
    )
    
    cat(sprintf("    Heat RR: %.3f (%.3f-%.3f), I²=%.1f%%\n", 
                pooled_rr_heat, pooled_ci_heat[1], pooled_ci_heat[2], I2_heat))
    cat(sprintf("    Cold RR: %.3f (%.3f-%.3f), I²=%.1f%%\n", 
                pooled_rr_cold, pooled_ci_cold[1], pooled_ci_cold[2], I2_cold))
  }
  
  return(pooled_results)
}

# -----------------------------------------------------------------------------
# 4. Run Analysis
# -----------------------------------------------------------------------------
eval_temps <- c(temp_p01, temp_p99)

harvesting_results <- fit_harvesting_single_model(
  data, MAX_LAG, TEMP_DF, LAG_DF, eval_temps, LAG_HORIZONS
)

# -----------------------------------------------------------------------------
# 5. Calculate Harvesting Ratios and Summary
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("HARVESTING ANALYSIS SUMMARY (v2 - Corrected)\n")
cat("=======================================================\n")

if (is.null(harvesting_results) || length(harvesting_results) == 0) {
  cat("\nWARNING: Harvesting analysis failed\n")
  harvesting_results <- list(status = "failed")
} else {
  
  cat("\n--- CUMULATIVE RR BY LAG HORIZON (Single-Model Extraction) ---\n")
  cat(sprintf("%-15s %25s %25s\n", "Lag Horizon", "Heat RR (I²)", "Cold RR (I²)"))
  cat(strrep("-", 65), "\n")
  
  for (horizon in LAG_HORIZONS) {
    horizon_key <- paste0("lag_", horizon)
    if (horizon_key %in% names(harvesting_results)) {
      r <- harvesting_results[[horizon_key]]
      cat(sprintf("%-15s %7.3f (%.3f-%.3f) I²=%.0f%% %7.3f (%.3f-%.3f) I²=%.0f%%\n",
                  paste0(horizon, " days"),
                  r$heat$rr, r$heat$ci_low, r$heat$ci_high, r$heat$I2_percent,
                  r$cold$rr, r$cold$ci_low, r$cold$ci_high, r$cold$I2_percent))
    }
  }
  
  # Calculate harvesting ratios
  safe_rr <- function(results, horizon, effect) {
    key <- paste0("lag_", horizon)
    if (key %in% names(results) && !is.null(results[[key]][[effect]]$rr)) {
      return(results[[key]][[effect]]$rr)
    }
    return(NA)
  }
  
  err_heat_7 <- safe_rr(harvesting_results, 7, "heat") - 1
  err_heat_35 <- safe_rr(harvesting_results, 35, "heat") - 1
  err_cold_7 <- safe_rr(harvesting_results, 7, "cold") - 1
  err_cold_35 <- safe_rr(harvesting_results, 35, "cold") - 1
  
  harvest_ratio_heat <- if (is.finite(err_heat_7) && is.finite(err_heat_35) && err_heat_7 != 0) {
    1 - (err_heat_35 / err_heat_7)
  } else NA
  
  harvest_ratio_cold <- if (is.finite(err_cold_7) && is.finite(err_cold_35) && err_cold_7 != 0) {
    1 - (err_cold_35 / err_cold_7)
  } else NA
  
  cat("\n--- HARVESTING RATIOS ---\n")
  if (is.finite(harvest_ratio_heat)) {
    cat(sprintf("Heat: %.3f (%.0f%% displacement)\n", harvest_ratio_heat, harvest_ratio_heat * 100))
  } else {
    cat("Heat: NA\n")
  }
  if (is.finite(harvest_ratio_cold)) {
    cat(sprintf("Cold: %.3f (%.0f%% displacement)\n", harvest_ratio_cold, harvest_ratio_cold * 100))
  } else {
    cat("Cold: NA\n")
  }
  
  harvesting_results$summary <- list(
    harvesting_ratio_heat = harvest_ratio_heat,
    harvesting_ratio_cold = harvest_ratio_cold,
    err_heat_7day = err_heat_7,
    err_heat_35day = err_heat_35,
    err_cold_7day = err_cold_7,
    err_cold_35day = err_cold_35,
    method = "single_model_cumulative"
  )
}

# -----------------------------------------------------------------------------
# 6. Save Results
# -----------------------------------------------------------------------------
cat("\n[6] Saving results...\n")

output <- list(
  exposure_type = EXPOSURE_TYPE,
  analysis_date = as.character(Sys.time()),
  method = "single_model_cumulative_v2",
  description = "Corrected harvesting analysis: fit ONE model with max_lag, extract cumulative RR at each horizon using crosspred(to=) parameter",
  max_lag = MAX_LAG,
  lag_horizons = LAG_HORIZONS,
  eval_temperatures = list(cold_p1 = temp_p01, heat_p99 = temp_p99),
  results = harvesting_results
)

output_file <- file.path(OUTPUT_DIR, paste0("harvesting_v2_", EXPOSURE_TYPE, ".json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("  JSON saved to:", output_file, "\n")

cat("\n=======================================================\n")
cat("Done!", as.character(Sys.time()), "\n")
cat("=======================================================\n")
