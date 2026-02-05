# =============================================================================
# 02b_harvesting.R
# Harvesting / Mortality Displacement Analysis
# =============================================================================
#
# Tests whether temperature-attributable deaths represent:
# - TRUE EXCESS: Deaths that would not have occurred otherwise
# - HARVESTING: Short-term mortality displacement
#
# Method (Armstrong et al. 2014, Gasparrini et al. 2015):
# 1. Fit DLNM with extended lag (35 days) for each region
# 2. Calculate CUMULATIVE RR at different lag horizons (7, 14, 21, 28, 35)
# 3. Pool via mixmeta
# 4. Harvesting ratio = 1 - (ERR_35 / ERR_7) where ERR = RR - 1
#    - Positive ratio → harvesting (displacement)
#    - Zero/negative → true excess mortality
#
# References:
# - Armstrong et al. (2014) Am J Epidemiol
# - Gasparrini et al. (2015) Lancet
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
cat("HARVESTING ANALYSIS:", EXPOSURE_TYPE, "\n")
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
cat("  Spline df: temp=", TEMP_DF, ", lag=", LAG_DF, "\n\n")

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
if (region_col %in% names(era5)) {
  setnames(era5, region_col, "region_code")
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

# Temperature thresholds for evaluation
temp_p99 <- quantile(temp_all, 0.99, na.rm = TRUE)
temp_p01 <- quantile(temp_all, 0.01, na.rm = TRUE)

cat("  Regions:", uniqueN(data$region_code), "\n")
cat("  Date range:", as.character(min(data$date)), "to", as.character(max(data$date)), "\n")
cat("  GLOBAL knots (P10/P75/P90):", round(temp_pcts, 1), "\n")
cat("  Boundary knots:", temp_boundary, "\n")
cat("  Temperature P99 (heat):", round(temp_p99, 2), "°C\n")
cat("  Temperature P1 (cold):", round(temp_p01, 2), "°C\n")

# -----------------------------------------------------------------------------
# 2. Function to Fit Extended DLNM for Harvesting Analysis
# -----------------------------------------------------------------------------
fit_harvesting_dlnm <- function(data_dt, max_lag, temp_df, lag_df, 
                                 eval_temps, lag_horizons) {
  
  regions <- unique(data_dt$region_code)
  n_regions <- length(regions)
  
  cat(sprintf("\n[2] Fitting extended DLNM (lag=%d) for %d regions...\n", 
              max_lag, n_regions))
  
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
  
  n_success <- 0
  n_errors <- 0
  n_bad_vcov <- 0
  
  # Wrap entire loop in tryCatch to catch any unexpected errors
  loop_result <- tryCatch({
    for (i in seq_along(regions)) {
      reg <- regions[i]
      
      if (i %% 20 == 0 || i <= 3) {
        cat(sprintf("  Progress: %d/%d regions (success so far: heat=%d, cold=%d, errors: %d)\n", 
                    i, n_regions, 
                    length(results_by_horizon$lag_7$heat_coefs),
                    length(results_by_horizon$lag_7$cold_coefs),
                    n_errors))
      }
      
      # Subset data
      daily <- data_dt[region_code == reg][order(date)]
      n_days <- nrow(daily)
      
      if (n_days < 365) next
      
      result <- tryCatch({
        # Use GLOBAL knots for meta-analysis comparability (expert recommendation)
        reg_p1 <- quantile(daily$tmean, 0.01, na.rm = TRUE)
      reg_p99 <- quantile(daily$tmean, 0.99, na.rm = TRUE)
      
      # Create cross-basis with extended lag using GLOBAL knots
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
      
      # Quality check
      vcov_ok <- all(diag(cb_vcov) > 0) && all(is.finite(cb_vcov))
      
      if (!vcov_ok) {
        n_bad_vcov <<- n_bad_vcov + 1
        stop("bad_vcov")  # Use stop() instead of return() to stay in tryCatch
      }
      
      # Find MMT - restrict to region's observed range (1st-99th percentile)
      mmt_search_min <- max(temp_boundary[1], reg_p1)
      mmt_search_max <- min(temp_boundary[2], reg_p99)
      temp_seq <- seq(mmt_search_min, mmt_search_max, length.out = 100)
      cp_full <- crosspred(cb, model, at = temp_seq, cumul = TRUE, cen = median(daily$tmean))
      mmt_idx <- which.min(cp_full$allRRfit)
      mmt <- temp_seq[mmt_idx]
      
      # For each lag horizon, compute cumulative RR
      for (horizon in lag_horizons) {
        # Create a cross-basis restricted to this lag horizon (using GLOBAL knots)
        cb_horizon <- crossbasis(
          daily$tmean,
          lag = horizon,  # Shorter lag for this horizon
          argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
          arglag = list(fun = "ns", df = lag_df)
        )
        
        # Refit model with this shorter lag
        model_horizon <- glm(
          deaths ~ cb_horizon + ns(as.numeric(date), df = 7 * length(unique(format(daily$date, "%Y")))) +
            factor(format(date, "%u")),
          data = daily,
          family = quasipoisson(link = "log"),
          na.action = na.exclude
        )
        
        # Get cumulative prediction
        cp_horizon <- crosspred(
          cb_horizon, 
          model_horizon, 
          at = eval_temps,
          cumul = TRUE,
          cen = mmt
        )
        
        # Extract coefficients for this cumulative window
        # Heat effect (P99)
        heat_rr <- cp_horizon$allRRfit[2]  # eval_temps[2] = P99
        heat_se <- (cp_horizon$allRRhigh[2] - cp_horizon$allRRlow[2]) / (2 * 1.96)
        
        # Cold effect (P1)
        cold_rr <- cp_horizon$allRRfit[1]  # eval_temps[1] = P1
        cold_se <- (cp_horizon$allRRhigh[1] - cp_horizon$allRRlow[1]) / (2 * 1.96)
        
        # Store for meta-analysis
        horizon_key <- paste0("lag_", horizon)
        
        # Convert RR to log scale for pooling
        if (is.finite(heat_rr) && heat_rr > 0 && is.finite(heat_se) && heat_se > 0) {
          results_by_horizon[[horizon_key]]$heat_coefs[[length(results_by_horizon[[horizon_key]]$heat_coefs) + 1]] <- log(heat_rr)
          results_by_horizon[[horizon_key]]$heat_vcovs[[length(results_by_horizon[[horizon_key]]$heat_vcovs) + 1]] <- heat_se^2
          if (i == 1 && horizon == 7) {
            cat(sprintf("  DEBUG: First region heat stored - RR=%.3f, SE=%.4f\n", heat_rr, heat_se))
          }
        }
        
        if (is.finite(cold_rr) && cold_rr > 0 && is.finite(cold_se) && cold_se > 0) {
          results_by_horizon[[horizon_key]]$cold_coefs[[length(results_by_horizon[[horizon_key]]$cold_coefs) + 1]] <- log(cold_rr)
          results_by_horizon[[horizon_key]]$cold_vcovs[[length(results_by_horizon[[horizon_key]]$cold_vcovs) + 1]] <- cold_se^2
          if (i == 1 && horizon == 7) {
            cat(sprintf("  DEBUG: First region cold stored - RR=%.3f, SE=%.4f\n", cold_rr, cold_se))
          }
        }
        
        # Track region for this horizon
        if ((is.finite(heat_rr) && heat_rr > 0 && is.finite(heat_se) && heat_se > 0) || 
            (is.finite(cold_rr) && cold_rr > 0 && is.finite(cold_se) && cold_se > 0)) {
          results_by_horizon[[horizon_key]]$regions <- c(
            results_by_horizon[[horizon_key]]$regions,
            as.character(reg)
          )
        }
      }
      
      "success"
      
    }, error = function(e) {
      n_errors <<- n_errors + 1
      if (n_errors == 1) {
        cat("  First error message:", conditionMessage(e), "\n")
      }
      "error"
    })
  }
  "loop_completed"
  }, error = function(e) {
    cat("  FATAL ERROR in main loop:", conditionMessage(e), "\n")
    cat("  Error occurred around region", i, "\n")
    "loop_error"
  })
  
  if (loop_result == "loop_error") {
    cat("  Attempting to continue with partial results...\n")
  }
  
  # Count successful regions by checking if they have stored coefficients
  n_success <- length(unique(unlist(lapply(lag_horizons, function(h) {
    results_by_horizon[[paste0("lag_", h)]]$regions
  }))))
  
  # Debug: show how many regions stored for each horizon
  cat("\n  Region counts by horizon:\n")
  for (h in lag_horizons) {
    horizon_key <- paste0("lag_", h)
    n_heat <- length(results_by_horizon[[horizon_key]]$heat_coefs)
    n_cold <- length(results_by_horizon[[horizon_key]]$cold_coefs)
    cat(sprintf("    Lag %d: heat=%d, cold=%d\n", h, n_heat, n_cold))
  }
  
  cat(sprintf("  Unique successful regions: %d/%d (errors: %d, bad_vcov: %d)\n", 
              n_success, n_regions, n_errors, n_bad_vcov))
  
  if (n_success < 10) {
    cat("  WARNING: Too few regions (<10) for reliable pooling\n")
    # Still try to pool if we have at least some data
    if (n_success < 3) {
      return(NULL)
    }
    cat("  Attempting pooling with", n_success, "regions...\n")
  }
  
  # -----------------------------------------------------------------------------
  # 3. Pool Results via Mixmeta for Each Lag Horizon
  # -----------------------------------------------------------------------------
  cat("\n[3] Pooling results via mixmeta for each lag horizon...\n")
  
  pooled_results <- list()
  
  for (horizon in lag_horizons) {
    horizon_key <- paste0("lag_", horizon)
    cat(sprintf("  Pooling lag horizon %d days...\n", horizon))
    
    # Heat effects - filter out Inf/NA values
    heat_coefs <- unlist(results_by_horizon[[horizon_key]]$heat_coefs)
    heat_vcovs <- unlist(results_by_horizon[[horizon_key]]$heat_vcovs)
    
    # Remove Inf and extreme values
    valid_heat <- is.finite(heat_coefs) & is.finite(heat_vcovs) & heat_vcovs > 0 & heat_vcovs < 100
    heat_coefs <- heat_coefs[valid_heat]
    heat_vcovs <- heat_vcovs[valid_heat]
    
    cat(sprintf("    Heat: n=%d valid coefs (after filtering)\n", length(heat_coefs)))
    
    if (length(heat_coefs) > 5) {
      mv_heat <- tryCatch({
        mixmeta(heat_coefs, heat_vcovs, method = "reml",
                control = list(maxiter = 500, showiter = FALSE))
      }, error = function(e) {
        cat(sprintf("      REML failed: %s\n", conditionMessage(e)))
        tryCatch({
          mixmeta(heat_coefs, heat_vcovs, method = "ml",
                  control = list(maxiter = 500, showiter = FALSE))
        }, error = function(e2) {
          cat(sprintf("      ML also failed: %s\n", conditionMessage(e2)))
          NULL
        })
      })
      
      if (!is.null(mv_heat)) {
        pooled_log_rr_heat <- coef(mv_heat)[1]
        pooled_se_heat <- sqrt(vcov(mv_heat)[1, 1])
        pooled_rr_heat <- exp(pooled_log_rr_heat)
        pooled_ci_heat <- exp(pooled_log_rr_heat + c(-1.96, 1.96) * pooled_se_heat)
        # Heterogeneity for heat using qtest() properly
        qtest_heat <- tryCatch({
          q <- qtest(mv_heat)
          n <- length(q$Q)  # Overall is last element
          list(Q = q$Q[n], df = q$df[n], pvalue = q$pvalue[n])
        }, error = function(e) list(Q = NA, df = NA, pvalue = NA))
        Q_heat <- qtest_heat$Q
        df_heat <- qtest_heat$df
        I2_heat <- if(!is.na(Q_heat) && !is.na(df_heat) && df_heat > 0) max(0, (Q_heat - df_heat) / Q_heat * 100) else NA
      } else {
        pooled_rr_heat <- NA
        pooled_ci_heat <- c(NA, NA)
        Q_heat <- NA; df_heat <- NA; I2_heat <- NA
      }
    } else {
      pooled_rr_heat <- NA
      pooled_ci_heat <- c(NA, NA)
      Q_heat <- NA; df_heat <- NA; I2_heat <- NA
    }
    
    # Cold effects - filter out Inf/NA values
    cold_coefs <- unlist(results_by_horizon[[horizon_key]]$cold_coefs)
    cold_vcovs <- unlist(results_by_horizon[[horizon_key]]$cold_vcovs)
    
    # Remove Inf and extreme values
    valid_cold <- is.finite(cold_coefs) & is.finite(cold_vcovs) & cold_vcovs > 0 & cold_vcovs < 100
    cold_coefs <- cold_coefs[valid_cold]
    cold_vcovs <- cold_vcovs[valid_cold]
    
    cat(sprintf("    Cold: n=%d valid coefs (after filtering)\n", length(cold_coefs)))
    
    if (length(cold_coefs) > 5) {
      mv_cold <- tryCatch({
        mixmeta(cold_coefs, cold_vcovs, method = "reml",
                control = list(maxiter = 500, showiter = FALSE))
      }, error = function(e) {
        cat(sprintf("      REML failed: %s\n", conditionMessage(e)))
        tryCatch({
          mixmeta(cold_coefs, cold_vcovs, method = "ml",
                  control = list(maxiter = 500, showiter = FALSE))
        }, error = function(e2) {
          cat(sprintf("      ML also failed: %s\n", conditionMessage(e2)))
          NULL
        })
      })
      
      if (!is.null(mv_cold)) {
        pooled_log_rr_cold <- coef(mv_cold)[1]
        pooled_se_cold <- sqrt(vcov(mv_cold)[1, 1])
        pooled_rr_cold <- exp(pooled_log_rr_cold)
        pooled_ci_cold <- exp(pooled_log_rr_cold + c(-1.96, 1.96) * pooled_se_cold)
        # Heterogeneity for cold using qtest() properly
        qtest_cold <- tryCatch({
          q <- qtest(mv_cold)
          n <- length(q$Q)  # Overall is last element
          list(Q = q$Q[n], df = q$df[n], pvalue = q$pvalue[n])
        }, error = function(e) list(Q = NA, df = NA, pvalue = NA))
        Q_cold <- qtest_cold$Q
        df_cold <- qtest_cold$df
        I2_cold <- if(!is.na(Q_cold) && !is.na(df_cold) && df_cold > 0) max(0, (Q_cold - df_cold) / Q_cold * 100) else NA
      } else {
        pooled_rr_cold <- NA
        pooled_ci_cold <- c(NA, NA)
        Q_cold <- NA; df_cold <- NA; I2_cold <- NA
      }
    } else {
      pooled_rr_cold <- NA
      pooled_ci_cold <- c(NA, NA)
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
    
    cat(sprintf("    Heat RR: %.3f (%.3f-%.3f)\n", 
                pooled_rr_heat, pooled_ci_heat[1], pooled_ci_heat[2]))
    cat(sprintf("    Cold RR: %.3f (%.3f-%.3f)\n", 
                pooled_rr_cold, pooled_ci_cold[1], pooled_ci_cold[2]))
  }
  
  return(pooled_results)
}

# -----------------------------------------------------------------------------
# 4. Run Harvesting Analysis
# -----------------------------------------------------------------------------
eval_temps <- c(temp_p01, temp_p99)  # Cold, Heat

harvesting_results <- fit_harvesting_dlnm(
  data, 
  MAX_LAG, 
  TEMP_DF, 
  LAG_DF,
  eval_temps,
  LAG_HORIZONS
)

# -----------------------------------------------------------------------------
# 5. Calculate Harvesting Ratios
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("HARVESTING ANALYSIS SUMMARY\n")
cat("=======================================================\n")

if (is.null(harvesting_results) || length(harvesting_results) == 0) {
  cat("\nWARNING: Harvesting analysis failed - insufficient valid regions\n")
  cat("This may occur when extended lag (35 days) causes convergence issues.\n")
  cat("Consider using the main DLNM results (21-day lag) for interpretation.\n")
  
  # Create empty results for JSON output
  harvesting_results <- list(
    status = "failed",
    message = "Insufficient valid regions for extended lag analysis"
  )
  
} else if (!all(paste0("lag_", LAG_HORIZONS) %in% names(harvesting_results))) {
  cat("\nWARNING: Some lag horizons failed to pool\n")
  cat("Available horizons:", paste(names(harvesting_results), collapse = ", "), "\n")
  
} else {
  
  cat("\n--- CUMULATIVE RR BY LAG HORIZON ---\n")
  cat(sprintf("%-15s %20s %20s\n", "Lag Horizon", "Heat RR (P99)", "Cold RR (P1)"))
  cat(strrep("-", 55), "\n")
  
  for (horizon in LAG_HORIZONS) {
    horizon_key <- paste0("lag_", horizon)
    if (horizon_key %in% names(harvesting_results)) {
      r <- harvesting_results[[horizon_key]]
      cat(sprintf("%-15s %7.3f (%.3f-%.3f) %7.3f (%.3f-%.3f)\n",
                  paste0(horizon, " days"),
                  r$heat$rr, r$heat$ci_low, r$heat$ci_high,
                  r$cold$rr, r$cold$ci_low, r$cold$ci_high))
    } else {
      cat(sprintf("%-15s %20s %20s\n", paste0(horizon, " days"), "NA", "NA"))
    }
  }
  
  # Calculate harvesting ratios
  cat("\n--- HARVESTING RATIOS ---\n")
  
  # Safe access helper
  safe_rr <- function(results, horizon, effect) {
    key <- paste0("lag_", horizon)
    if (key %in% names(results) && !is.null(results[[key]][[effect]]$rr)) {
      return(results[[key]][[effect]]$rr)
    }
    return(NA)
  }
  
  # ERR = RR - 1
  err_heat_7 <- safe_rr(harvesting_results, 7, "heat") - 1
  err_heat_35 <- safe_rr(harvesting_results, 35, "heat") - 1
  err_cold_7 <- safe_rr(harvesting_results, 7, "cold") - 1
  err_cold_35 <- safe_rr(harvesting_results, 35, "cold") - 1
  
  # Harvesting ratio = 1 - (ERR_35 / ERR_7)
  harvest_ratio_heat <- if (is.finite(err_heat_7) && is.finite(err_heat_35)) {
    1 - (err_heat_35 / err_heat_7)
  } else {
    NA
  }
  
  harvest_ratio_cold <- if (is.finite(err_cold_7) && is.finite(err_cold_35)) {
    1 - (err_cold_35 / err_cold_7)
  } else {
    NA
  }
  
  if (is.finite(harvest_ratio_heat)) {
    cat(sprintf("Heat harvesting ratio (7→35 days): %.3f\n", harvest_ratio_heat))
  } else {
    cat("Heat harvesting ratio (7→35 days): NA (pooling failed for some horizons)\n")
  }
  
  if (is.finite(harvest_ratio_cold)) {
    cat(sprintf("Cold harvesting ratio (7→35 days): %.3f\n", harvest_ratio_cold))
  } else {
    cat("Cold harvesting ratio (7→35 days): NA (pooling failed for some horizons)\n")
  }
  
  cat("\nInterpretation:\n")
  cat("  > 0: Evidence of harvesting (mortality displacement)\n")
  cat("  ≈ 0: True excess mortality (no displacement)\n")
  cat("  < 0: Increasing cumulative effect (no harvesting)\n")
  
  # Determine interpretation
  if (is.finite(harvest_ratio_heat)) {
    if (harvest_ratio_heat > 0.2) {
      cat("\n  Heat: STRONG harvesting detected (>20% displacement)\n")
    } else if (harvest_ratio_heat > 0.05) {
      cat("\n  Heat: MODERATE harvesting detected (5-20% displacement)\n")
    } else if (harvest_ratio_heat >= -0.05) {
      cat("\n  Heat: TRUE EXCESS mortality (minimal displacement)\n")
    } else {
      cat("\n  Heat: INCREASING cumulative effect (no harvesting)\n")
    }
  } else {
    cat("\n  Heat: Cannot determine (insufficient valid data)\n")
  }
  
  if (is.finite(harvest_ratio_cold)) {
    if (harvest_ratio_cold > 0.2) {
      cat("  Cold: STRONG harvesting detected (>20% displacement)\n")
    } else if (harvest_ratio_cold > 0.05) {
      cat("  Cold: MODERATE harvesting detected (5-20% displacement)\n")
    } else if (harvest_ratio_cold >= -0.05) {
      cat("  Cold: TRUE EXCESS mortality (minimal displacement)\n")
    } else {
      cat("  Cold: INCREASING cumulative effect (no harvesting)\n")
    }
  } else {
    cat("  Cold: Cannot determine (insufficient valid data)\n")
  }
  
  # Add to results
  harvesting_results$summary <- list(
    harvesting_ratio_heat = harvest_ratio_heat,
    harvesting_ratio_cold = harvest_ratio_cold,
    err_heat_7day = err_heat_7,
    err_heat_35day = err_heat_35,
    err_cold_7day = err_cold_7,
    err_cold_35day = err_cold_35
  )
}

# -----------------------------------------------------------------------------
# 6. Save Results
# -----------------------------------------------------------------------------
cat("\n[6] Saving results...\n")

output <- list(
  exposure_type = EXPOSURE_TYPE,
  analysis_date = as.character(Sys.time()),
  max_lag = MAX_LAG,
  lag_horizons = LAG_HORIZONS,
  eval_temperatures = list(
    cold_p1 = temp_p01,
    heat_p99 = temp_p99
  ),
  results = harvesting_results
)

output_file <- file.path(OUTPUT_DIR, paste0("harvesting_r_", EXPOSURE_TYPE, ".json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("  JSON saved to:", output_file, "\n")

cat("\n=======================================================\n")
cat("Done!", as.character(Sys.time()), "\n")
cat("=======================================================\n")
