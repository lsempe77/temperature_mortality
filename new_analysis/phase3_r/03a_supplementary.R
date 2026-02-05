# =============================================================================
# 03a_supplementary.R
# Supplementary Confounding Control Analyses
# =============================================================================
#
# Tests robustness of temperature-mortality associations to potential confounders:
# 1. Apparent temperature (humidity-adjusted) vs dry-bulb temperature
# 2. Air pollution controls (PM2.5, O3)
# 3. Influenza season controls
#
# NOTE: These are SUPPLEMENTARY analyses. Pollution may be a MEDIATOR
# (heat → stagnation → pollution), not just a confounder. Over-adjustment
# could bias results downward.
#
# Method:
# - Fit separate DLNM models for each specification
# - Compare pooled RRs across specifications
# - Assess whether confounding substantially changes estimates
#
# References:
# - Gasparrini et al. (2015) Lancet - Multi-country methodology
# - Armstrong (2006) Epidemiology - Pollution as mediator discussion
#
# =============================================================================

suppressPackageStartupMessages({
  library(dlnm)
  library(mixmeta)
  library(data.table)
  library(arrow)
  library(jsonlite)
  library(splines)
  library(zoo)
})

# Get exposure type from command line args
args <- commandArgs(trailingOnly = TRUE)
EXPOSURE_TYPE <- if (length(args) > 0) args[1] else "intermediate"

cat("=======================================================\n")
cat("SUPPLEMENTARY CONFOUNDING ANALYSES:", EXPOSURE_TYPE, "\n")
cat("=======================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# Configuration
MAX_LAG <- 21
TEMP_DF <- 4
LAG_DF <- 4
MIN_POLLUTION_COVERAGE <- 0.7  # Require 70% non-missing pollution data

cat("Configuration:\n")
cat("  Max lag:", MAX_LAG, "days\n")
cat("  Spline df: temp=", TEMP_DF, ", lag=", LAG_DF, "\n")
cat("  Min pollution coverage:", MIN_POLLUTION_COVERAGE * 100, "%\n\n")

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
BASE_DIR <- dirname(SCRIPT_DIR)  # Go up from phase3_r to new_analysis

DATA_DIR <- file.path(BASE_DIR, "phase0_data_prep", "results")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

cat("[1] Loading data...\n")

# Mortality data
mort_file <- file.path(DATA_DIR, paste0("mortality_", EXPOSURE_TYPE, "_daily_elderly.parquet"))
mort <- as.data.table(read_parquet(mort_file))

# ERA5 temperature data (includes dewpoint for apparent temp)
era5_file <- file.path(DATA_DIR, paste0("era5_", EXPOSURE_TYPE, "_daily.parquet"))
era5 <- as.data.table(read_parquet(era5_file))

# CAMS pollution data
cams_file <- file.path(DATA_DIR, paste0("cams_", EXPOSURE_TYPE, "_daily.parquet"))
if (file.exists(cams_file)) {
  cams <- as.data.table(read_parquet(cams_file))
  # Convert datetime to date
  cams[, date := as.Date(date)]
  has_pollution <- TRUE
  cat("  Pollution data loaded\n")
} else {
  has_pollution <- FALSE
  cat("  WARNING: No pollution data found\n")
}

# Influenza data
flu_file <- file.path(DATA_DIR, paste0("influenza_daily_by_", EXPOSURE_TYPE, "_region.parquet"))
if (file.exists(flu_file)) {
  flu <- as.data.table(read_parquet(flu_file))
  has_flu <- TRUE
  cat("  Influenza data loaded\n")
} else {
  has_flu <- FALSE
  cat("  WARNING: No influenza data found\n")
}

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
if (has_pollution && region_col %in% names(cams)) {
  setnames(cams, region_col, "region_code")
  cat(sprintf("  Renamed %s to region_code in CAMS\n", region_col))
}
if (has_flu && region_col %in% names(flu)) {
  setnames(flu, region_col, "region_code")
}

# Convert dates for merge compatibility (mort has POSIXct, era5 has Date)
mort[, date := as.character(as.Date(date))]
era5[, date := as.character(date)]

# Merge mortality + temperature
setkey(mort, region_code, date)
setkey(era5, region_code, date)
data <- merge(mort, era5, by = c("region_code", "date"))

# Convert date back to Date for analysis
data[, date := as.Date(date)]

# Working variables
data[, deaths := deaths_elderly]
data[, tmean := temp_mean]

# Merge pollution if available
if (has_pollution) {
  setkey(cams, region_code, date)
  # Rename ozone column if needed
  if ("ozone" %in% names(cams)) {
    setnames(cams, "ozone", "o3")
  }
  # Select available pollution columns
  poll_cols <- intersect(c("pm25", "pm10", "o3"), names(cams))
  data <- merge(data, cams[, c("region_code", "date", poll_cols), with = FALSE], 
                by = c("region_code", "date"), all.x = TRUE)
}

# Merge influenza if available
if (has_flu) {
  setkey(flu, region_code, date)
  data <- merge(data, flu[, .(region_code, date, srag_cases)], 
                by = c("region_code", "date"), all.x = TRUE)
}

# Keep only complete cases for core variables
data <- data[!is.na(deaths) & !is.na(tmean)]
setorder(data, region_code, date)

cat("  Regions:", uniqueN(data$region_code), "\n")
cat("  Date range:", as.character(min(data$date)), "to", as.character(max(data$date)), "\n")

# -----------------------------------------------------------------------------
# 2. Calculate Apparent Temperature
# -----------------------------------------------------------------------------
cat("\n[2] Calculating apparent temperature...\n")

if ("dewpoint_mean" %in% names(data)) {
  # Calculate relative humidity from dewpoint and air temperature
  # Magnus-Tetens approximation
  Td <- data$dewpoint_mean
  Ta <- data$tmean
  
  # Saturation vapor pressure
  es_Td <- 6.112 * exp((17.67 * Td) / (Td + 243.5))
  es_Ta <- 6.112 * exp((17.67 * Ta) / (Ta + 243.5))
  
  # Relative humidity (%)
  RH <- 100 * (es_Td / es_Ta)
  RH <- pmin(pmax(RH, 0), 100)  # Clip to [0, 100]
  
  # Apparent temperature (Australian Bureau of Meteorology formula)
  # AT = Ta + 0.33 × e - 4.0
  # where e = water vapor pressure
  e <- (RH / 100) * es_Ta
  data[, apparent_temp := tmean + 0.33 * e - 4.0]
  
  at_diff <- data$apparent_temp - data$tmean
  cat(sprintf("  AT - Ta difference: mean=%.1f°C, sd=%.1f°C, range=[%.1f, %.1f]\n",
              mean(at_diff, na.rm=TRUE), sd(at_diff, na.rm=TRUE),
              min(at_diff, na.rm=TRUE), max(at_diff, na.rm=TRUE)))
} else {
  data[, apparent_temp := tmean]
  cat("  WARNING: No dewpoint data, using dry-bulb as apparent temp\n")
}

# -----------------------------------------------------------------------------
# 3. Handle Missing Data (Interpolation)
# -----------------------------------------------------------------------------
cat("\n[3] Handling missing pollution and influenza data...\n")

# Interpolate pollution within regions (max 14 day gap)
# NOTE: CAMS has complete coverage, but keep interpolation for robustness
if (has_pollution) {
  for (pollutant in c("pm25", "pm10", "o3")) {
    if (pollutant %in% names(data)) {
      n_missing_before <- sum(is.na(data[[pollutant]]))
      
      # Simple check - if no missing data, skip interpolation
      if (n_missing_before == 0) {
        cat(sprintf("  %s: No missing data (coverage: 100.0%%)\n", pollutant))
        next
      }
      
      # Interpolate within region (not needed for CAMS but kept for other data)
      data[, (pollutant) := {
        x <- get(pollutant)
        if (sum(!is.na(x)) > 0) {
          zoo::na.approx(x, maxgap = 14, na.rm = FALSE)
        } else {
          x
        }
      }, by = region_code]
      
      n_missing_after <- sum(is.na(data[[pollutant]]))
      coverage <- 1 - n_missing_after / nrow(data)
      cat(sprintf("  %s: %d -> %d missing (coverage: %.1f%%)\n",
                  pollutant, n_missing_before, n_missing_after, coverage * 100))
    }
  }
}

# Interpolate influenza within regions
if (has_flu) {
  if ("srag_cases" %in% names(data)) {
    n_missing_before <- sum(is.na(data$srag_cases))
    
    if (n_missing_before > 0) {
      data[, srag_cases := {
        x <- srag_cases
        if (sum(!is.na(x)) > 0) {
          zoo::na.approx(x, maxgap = 14, na.rm = FALSE)
        } else {
          x
        }
      }, by = region_code]
      
      # Fill remaining with 0 (no cases)
      data[is.na(srag_cases), srag_cases := 0]
      
      n_missing_after <- sum(is.na(data$srag_cases))
      cat(sprintf("  srag_cases: %d -> %d missing\n", n_missing_before, n_missing_after))
    } else {
      cat("  srag_cases: No missing data\n")
    }
  }
}

# -----------------------------------------------------------------------------
# 4. Filter Regions by Pollution Coverage
# -----------------------------------------------------------------------------
cat("\n[4] Assessing pollution coverage by region...\n")

if (has_pollution) {
  # Calculate pollution coverage per region
  poll_coverage <- data[, .(
    pm25_coverage = sum(!is.na(pm25)) / .N,
    o3_coverage = sum(!is.na(o3)) / .N
  ), by = region_code]
  
  poll_coverage[, min_coverage := pmin(pm25_coverage, o3_coverage)]
  
  pollution_valid_regions <- poll_coverage[min_coverage >= MIN_POLLUTION_COVERAGE, region_code]
  
  cat(sprintf("  Regions with >=%.0f%% pollution coverage: %d/%d\n",
              MIN_POLLUTION_COVERAGE * 100, 
              length(pollution_valid_regions),
              uniqueN(data$region_code)))
} else {
  pollution_valid_regions <- character(0)
}

# Global temperature percentiles - use GLOBAL knots for meta-analysis comparability
temp_all <- data$tmean
temp_pcts <- quantile(temp_all, probs = c(0.10, 0.75, 0.90, 0.99), na.rm = TRUE)
temp_pcts_knots <- temp_pcts[c("10%", "75%", "90%")]  # For cross-basis
temp_boundary <- c(0, 40)  # Expanded from c(0.7, 35.1) per expert recommendation

cat("\nTemperature distribution:\n")
cat(sprintf("  GLOBAL knots (P10/P75/P90): %.1f / %.1f / %.1f\n",
            temp_pcts_knots[1], temp_pcts_knots[2], temp_pcts_knots[3]))
cat(sprintf("  Boundary knots: %.0f to %.0f\n", temp_boundary[1], temp_boundary[2]))
cat(sprintf("  P99: %.2f°C\n", temp_pcts["99%"]))

# -----------------------------------------------------------------------------
# 5. Fit DLNM Models with Different Specifications
# -----------------------------------------------------------------------------
cat("\n[5] Fitting DLNM models with different specifications...\n")

fit_supplementary_dlnm <- function(data_dt, temp_col, control_vars = c(), 
                                    valid_regions_filter = NULL) {
  
  regions <- unique(data_dt$region_code)
  
  # Filter to valid regions if specified
  if (!is.null(valid_regions_filter)) {
    regions <- intersect(regions, valid_regions_filter)
  }
  
  n_regions <- length(regions)
  
  cat(sprintf("  Fitting %d regions with %s and controls: %s\n", 
              n_regions, temp_col, 
              ifelse(length(control_vars) > 0, paste(control_vars, collapse=", "), "none")))
  
  # Store results
  rr_p99_list <- list()
  rr_p99_se_list <- list()
  rr_p1_list <- list()
  rr_p1_se_list <- list()
  valid_region_ids <- c()
  
  n_success <- 0
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
    
    # Check if control variables are available
    if (length(control_vars) > 0) {
      if (!all(control_vars %in% names(daily))) next
      if (any(sapply(control_vars, function(v) sum(!is.na(daily[[v]])) < n_days * 0.7))) next
    }
    
    result <- tryCatch({
      # Use GLOBAL knots for meta-analysis comparability (expert recommendation)
      temp_vector <- daily[[temp_col]]
      reg_temp_pcts <- quantile(temp_vector, probs = c(0.01, 0.10, 0.75, 0.90, 0.99), 
                                na.rm = TRUE)
      reg_p1 <- reg_temp_pcts["1%"]
      reg_p99 <- reg_temp_pcts["99%"]
      
      # Create cross-basis with GLOBAL knots
      cb <- crossbasis(
        temp_vector,
        lag = MAX_LAG,
        argvar = list(fun = "ns", knots = temp_pcts_knots, 
                     Boundary.knots = temp_boundary),
        arglag = list(fun = "ns", df = LAG_DF)
      )
      
      # Build formula
      formula_str <- "deaths ~ cb + ns(as.numeric(date), df = 7 * length(unique(format(date, '%Y')))) + factor(format(date, '%u'))"
      
      # Add control variables
      if (length(control_vars) > 0) {
        for (ctrl_var in control_vars) {
          # Use natural spline for continuous controls
          if (ctrl_var %in% c("pm25", "pm10", "o3")) {
            formula_str <- paste0(formula_str, " + ns(", ctrl_var, ", df=3)")
          } else {
            formula_str <- paste0(formula_str, " + ", ctrl_var)
          }
        }
      }
      
      # Fit GLM
      model <- glm(
        as.formula(formula_str),
        data = daily,
        family = quasipoisson(link = "log"),
        na.action = na.exclude
      )
      
      # Extract coefficients
      cb_idx <- grep("^cb", names(coef(model)))
      cb_coef <- coef(model)[cb_idx]
      cb_vcov <- vcov(model)[cb_idx, cb_idx]
      
      # Quality check
      vcov_ok <- all(diag(cb_vcov) > 0) && all(is.finite(cb_vcov))
      if (!vcov_ok) stop("bad_vcov")  # Use stop() to stay in tryCatch
      
      # Find MMT - restrict to region's observed range (1st-99th percentile)
      mmt_search_min <- max(temp_boundary[1], reg_p1)
      mmt_search_max <- min(temp_boundary[2], reg_p99)
      temp_seq <- seq(mmt_search_min, mmt_search_max, length.out = 100)
      cp_full <- crosspred(cb, model, at = temp_seq, cumul = TRUE, cen = median(temp_vector, na.rm=TRUE))
      mmt_idx <- which.min(cp_full$allRRfit)
      mmt <- temp_seq[mmt_idx]
      
      # Get RR at P99 and P1
      cp_extremes <- crosspred(cb, model, 
                               at = reg_temp_pcts[c("1%", "99%")], 
                               cumul = TRUE, cen = mmt)
      
      rr_p1 <- cp_extremes$allRRfit[1]
      rr_p1_se <- (cp_extremes$allRRhigh[1] - cp_extremes$allRRlow[1]) / (2 * 1.96)
      
      rr_p99 <- cp_extremes$allRRfit[2]
      rr_p99_se <- (cp_extremes$allRRhigh[2] - cp_extremes$allRRlow[2]) / (2 * 1.96)
      
      # Debug first region only
      if (i == 1) {
        cat(sprintf("      [Debug Region 1] RR_p99=%.3f, SE=%.3f, RR_p1=%.3f, SE=%.3f\n",
                    rr_p99, rr_p99_se, rr_p1, rr_p1_se))
      }
      
      # Store if valid - exclude Inf, NA, 0, negative values
      # Also exclude extreme SE values (> 20 on log scale protects against numerical failures)
      if (all(is.finite(c(rr_p99, rr_p99_se, rr_p1, rr_p1_se))) &&
          rr_p99 > 0 && rr_p1 > 0 && 
          rr_p99_se > 0 && rr_p1_se > 0 &&
          rr_p99_se < 20 && rr_p1_se < 20) {
        
        rr_p99_list[[length(rr_p99_list) + 1]] <- log(rr_p99)
        rr_p99_se_list[[length(rr_p99_se_list) + 1]] <- rr_p99_se / rr_p99  # SE of log(RR)
        
        rr_p1_list[[length(rr_p1_list) + 1]] <- log(rr_p1)
        rr_p1_se_list[[length(rr_p1_se_list) + 1]] <- rr_p1_se / rr_p1
        
        valid_region_ids <- c(valid_region_ids, as.character(reg))
        n_success <- n_success + 1
      }
      
      "success"
      
    }, error = function(e) {
      n_errors <- n_errors + 1
      "error"
    })
  }
  
  cat(sprintf("    Successful fits: %d/%d (errors: %d)\n", n_success, n_regions, n_errors))
  cat(sprintf("    Valid regions stored: %d (excluded: %d with Inf/extreme SE)\n", 
              length(valid_region_ids), n_regions - length(valid_region_ids)))
  
  if (length(valid_region_ids) < 10) {
    cat("    WARNING: Too few valid regions for pooling\n")
    return(NULL)
  }
  
  # Pool via mixmeta
  pool_effect <- function(log_rr_list, se_list, effect_name) {
    coefs <- unlist(log_rr_list)
    vcovs <- unlist(se_list)^2
    
    cat(sprintf("      [Pool %s] n=%d, coef range=[%.3f, %.3f], vcov range=[%.6f, %.6f]\n",
                effect_name, length(coefs), min(coefs), max(coefs), min(vcovs), max(vcovs)))
    
    mv_fit <- tryCatch({
      mixmeta(coefs, vcovs, method = "reml",
              control = list(maxiter = 500, showiter = FALSE))
    }, error = function(e) {
      cat(sprintf("      [Pool %s] REML failed: %s\n", effect_name, e$message))
      tryCatch({
        mixmeta(coefs, vcovs, method = "ml",
                control = list(maxiter = 500, showiter = FALSE))
      }, error = function(e2) {
        cat(sprintf("      [Pool %s] ML also failed: %s\n", effect_name, e2$message))
        NULL
      })
    })
    
    if (!is.null(mv_fit)) {
      pooled_log_rr <- coef(mv_fit)[1]
      pooled_se <- sqrt(vcov(mv_fit)[1, 1])
      pooled_rr <- exp(pooled_log_rr)
      pooled_ci <- exp(pooled_log_rr + c(-1.96, 1.96) * pooled_se)
      
      # Heterogeneity extraction using qtest() properly
      qstat_info <- tryCatch({
        q <- qtest(mv_fit)
        n <- length(q$Q)  # Overall is last element
        list(Q = q$Q[n], df = q$df[n], pvalue = q$pvalue[n])
      }, error = function(e) list(Q = NA, df = NA, pvalue = NA))
      
      I2 <- if(!is.na(qstat_info$Q) && !is.na(qstat_info$df) && qstat_info$df > 0) {
        max(0, (qstat_info$Q - qstat_info$df) / qstat_info$Q * 100)
      } else NA
      
      cat(sprintf("      [Pool %s] Q=%.1f, I²=%.1f%%\n", effect_name, qstat_info$Q, I2))
      
      list(
        rr = pooled_rr,
        ci_low = pooled_ci[1],
        ci_high = pooled_ci[2],
        log_rr = pooled_log_rr,
        se = pooled_se,
        heterogeneity = list(
          cochrans_Q = qstat_info$Q,
          Q_df = qstat_info$df,
          Q_pvalue = qstat_info$pvalue,
          I2_percent = I2
        )
      )
    } else {
      list(rr = NA, ci_low = NA, ci_high = NA, log_rr = NA, se = NA,
           heterogeneity = list(cochrans_Q = NA, Q_df = NA, Q_pvalue = NA, I2_percent = NA))
    }
  }
  
  heat_pooled <- pool_effect(rr_p99_list, rr_p99_se_list, "heat")
  cold_pooled <- pool_effect(rr_p1_list, rr_p1_se_list, "cold")
  
  cat(sprintf("    Heat (P99): RR = %.3f (%.3f-%.3f)\n",
              heat_pooled$rr, heat_pooled$ci_low, heat_pooled$ci_high))
  cat(sprintf("    Cold (P1): RR = %.3f (%.3f-%.3f)\n",
              cold_pooled$rr, cold_pooled$ci_low, cold_pooled$ci_high))
  
  list(
    n_regions = n_success,
    heat = heat_pooled,
    cold = cold_pooled
  )
}

# Model 1: Baseline (dry-bulb temperature, no extra controls)
cat("\n  MODEL 1: Baseline (dry-bulb temperature)\n")
model1_results <- fit_supplementary_dlnm(data, "tmean")

# Model 2: Apparent temperature
cat("\n  MODEL 2: Apparent temperature\n")
model2_results <- fit_supplementary_dlnm(data, "apparent_temp")

# Model 3: Dry-bulb + pollution controls (subset of regions with good coverage)
if (has_pollution && length(pollution_valid_regions) >= 10) {
  cat("\n  MODEL 3: Dry-bulb + pollution controls\n")
  model3_results <- fit_supplementary_dlnm(
    data, "tmean", 
    control_vars = c("pm25", "o3"),
    valid_regions_filter = pollution_valid_regions
  )
} else {
  model3_results <- NULL
  cat("\n  MODEL 3: SKIPPED (insufficient pollution data)\n")
}

# Model 4: Dry-bulb + influenza controls
if (has_flu) {
  cat("\n  MODEL 4: Dry-bulb + influenza controls\n")
  model4_results <- fit_supplementary_dlnm(
    data, "tmean",
    control_vars = c("srag_cases")
  )
} else {
  model4_results <- NULL
  cat("\n  MODEL 4: SKIPPED (no influenza data)\n")
}

# -----------------------------------------------------------------------------
# 6. Summary
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("SUPPLEMENTARY ANALYSIS SUMMARY\n")
cat("=======================================================\n")

compare_models <- function(baseline, comparison, model_name) {
  if (is.null(comparison)) {
    cat(sprintf("\n%s: NOT AVAILABLE\n", model_name))
    return()
  }
  
  cat(sprintf("\n%s:\n", model_name))
  cat(sprintf("  Regions: %d (baseline: %d)\n", comparison$n_regions, baseline$n_regions))
  cat(sprintf("  Heat RR: %.3f (%.3f-%.3f) vs baseline %.3f (%.3f-%.3f)\n",
              comparison$heat$rr, comparison$heat$ci_low, comparison$heat$ci_high,
              baseline$heat$rr, baseline$heat$ci_low, baseline$heat$ci_high))
  cat(sprintf("  Cold RR: %.3f (%.3f-%.3f) vs baseline %.3f (%.3f-%.3f)\n",
              comparison$cold$rr, comparison$cold$ci_low, comparison$cold$ci_high,
              baseline$cold$rr, baseline$cold$ci_low, baseline$cold$ci_high))
  
  # Calculate percent change
  heat_pct_change <- ((comparison$heat$rr - baseline$heat$rr) / baseline$heat$rr) * 100
  cold_pct_change <- ((comparison$cold$rr - baseline$cold$rr) / baseline$cold$rr) * 100
  
  cat(sprintf("  Change: Heat %.1f%%, Cold %.1f%%\n", heat_pct_change, cold_pct_change))
  
  if (abs(heat_pct_change) < 5 && abs(cold_pct_change) < 5) {
    cat("  → Results ROBUST to this specification\n")
  } else if (abs(heat_pct_change) < 10 && abs(cold_pct_change) < 10) {
    cat("  → Results MODERATELY sensitive to this specification\n")
  } else {
    cat("  → Results SENSITIVE to this specification\n")
  }
}

# Check if results are valid lists with expected structure
is_valid_result <- function(res) {
  !is.null(res) && is.list(res) && "n_regions" %in% names(res) && 
    "heat" %in% names(res) && "cold" %in% names(res)
}

if (is_valid_result(model1_results)) {
  cat("\nBASELINE (Dry-bulb temperature):\n")
  cat(sprintf("  Regions: %d\n", model1_results$n_regions))
  cat(sprintf("  Heat RR (P99): %.3f (%.3f-%.3f)\n",
              model1_results$heat$rr, model1_results$heat$ci_low, model1_results$heat$ci_high))
  cat(sprintf("  Cold RR (P1): %.3f (%.3f-%.3f)\n",
              model1_results$cold$rr, model1_results$cold$ci_low, model1_results$cold$ci_high))
  
  if (is_valid_result(model2_results)) compare_models(model1_results, model2_results, "Apparent Temperature")
  if (is_valid_result(model3_results)) compare_models(model1_results, model3_results, "Pollution-Adjusted")
  if (is_valid_result(model4_results)) compare_models(model1_results, model4_results, "Influenza-Adjusted")
} else {
  cat("\nWARNING: Baseline model failed to produce valid results\n")
  if (!is.null(model1_results)) {
    cat("  Result type:", class(model1_results), "\n")
  }
}

# -----------------------------------------------------------------------------
# 7. Save Results
# -----------------------------------------------------------------------------
cat("\n[7] Saving results...\n")

output <- list(
  exposure_type = EXPOSURE_TYPE,
  analysis_date = as.character(Sys.time()),
  
  baseline = model1_results,
  apparent_temp = model2_results,
  pollution_adjusted = model3_results,
  flu_adjusted = model4_results,
  
  metadata = list(
    pollution_regions = length(pollution_valid_regions),
    min_pollution_coverage = MIN_POLLUTION_COVERAGE
  )
)

output_file <- file.path(OUTPUT_DIR, paste0("supplementary_r_", EXPOSURE_TYPE, ".json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("  JSON saved to:", output_file, "\n")

cat("\n=======================================================\n")
cat("Done!", as.character(Sys.time()), "\n")
cat("=======================================================\n")
