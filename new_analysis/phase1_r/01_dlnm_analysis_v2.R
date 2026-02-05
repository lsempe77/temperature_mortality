# =============================================================================
# DLNM Analysis for Temperature-Mortality using R's dlnm package
# Version 2: Fixed JSON serialization and error handling
# =============================================================================

suppressPackageStartupMessages({
  library(dlnm)
  library(mixmeta)  # Successor to mvmeta (Gasparrini)
  library(data.table)
  library(arrow)
  library(jsonlite)
  library(splines)  # For ns() function
})

# Get exposure type from command line args
args <- commandArgs(trailingOnly = TRUE)
EXPOSURE_TYPE <- if (length(args) > 0) args[1] else "immediate"

cat("=======================================================\n")
cat("DLNM Analysis:", EXPOSURE_TYPE, "exposure\n")
cat("=======================================================\n")

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------
# Get script directory for proper path resolution
get_script_dir <- function() {
  # Try multiple methods to find script location
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  # Fallback: assume running from parent directory
  return("phase1_r")
}
script_dir <- get_script_dir()
cat("Script directory:", script_dir, "\n")

DATA_DIR <- normalizePath(file.path(script_dir, "../phase0_data_prep/results"), mustWork = FALSE)
OUTPUT_DIR <- normalizePath(file.path(script_dir, "results"), mustWork = FALSE)
cat("Data directory:", DATA_DIR, "\n")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\nLoading data...\n")

# Mortality data
mort_file <- file.path(DATA_DIR, paste0("mortality_", EXPOSURE_TYPE, "_daily_elderly.parquet"))
mort <- as.data.table(read_parquet(mort_file))
cat("  Mortality data:", nrow(mort), "rows\n")

# ERA5 temperature data
era5_file <- file.path(DATA_DIR, paste0("era5_", EXPOSURE_TYPE, "_daily.parquet"))
era5 <- as.data.table(read_parquet(era5_file))
cat("  ERA5 data:", nrow(era5), "rows\n")

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

# Merge datasets
setkey(mort, region_code, date)
setkey(era5, region_code, date)

data <- merge(mort, era5, by = c("region_code", "date"))
cat("  Merged data:", nrow(data), "rows\n")

# Convert date back to Date for analysis
data[, date := as.Date(date)]

# Create working variables
data[, deaths := deaths_elderly]
data[, tmean := temp_mean]
setorder(data, region_code, date)

# Remove missing values
data <- data[!is.na(deaths) & !is.na(tmean)]
cat("  After removing NAs:", nrow(data), "rows\n")

# -----------------------------------------------------------------------------
# 2. DLNM Parameters
# -----------------------------------------------------------------------------
# Following Gasparrini et al. (2015) Lancet methodology
MAX_LAG <- 21      # 21 days for lag (standard for temperature-mortality)
TEMP_DF <- 4       # df for temperature variable spline
LAG_DF <- 4        # df for lag spline
# Total parameters: TEMP_DF * LAG_DF = 16

# Global temperature percentiles for knots (3 internal knots for df=4)
temp_all <- data$tmean
temp_pcts <- quantile(temp_all, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)
# EXPANDED BOUNDS: Changed from [0.7, 35.1] to [0, 40] to allow MMT estimation
# at temperature extremes and prevent boundary convergence issues
temp_boundary <- c(0, 40)  # Expanded boundary for better MMT estimation

cat("  Temperature range:", round(min(temp_all), 1), "to", round(max(temp_all), 1), "\n")
cat("  Cross-basis boundary:", temp_boundary[1], "to", temp_boundary[2], "\n")
cat("  Knots:", paste(round(temp_pcts, 1), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 3. Fit DLNM by Region
# -----------------------------------------------------------------------------
regions <- unique(data$region_code)
n_regions <- length(regions)
cat("\nFitting DLNM for", n_regions, "regions...\n\n")

# Pre-allocate lists
all_results <- vector("list", n_regions)
mvmeta_coefs <- list()
mvmeta_vcovs <- list()
mvmeta_regions <- c()
success_count <- 0
mvmeta_count <- 0

for (i in seq_along(regions)) {
  reg <- regions[i]
  if (i %% 100 == 0) cat("Progress:", i, "/", n_regions, "\n")
  if (i <= 5) cat("Region", i, ":", reg, "\n")  # Debug first 5
  
  # Subset data
  daily <- data[region_code == reg][order(date)]
  n_days <- nrow(daily)
  
  # Skip if insufficient data
  if (n_days < 365) {
    all_results[[i]] <- list(status = "skip", reason = "insufficient_data", region_code = as.character(reg))
    next
  }
  
  result <- tryCatch({
    # =========================================================================
    # CRITICAL FIX: Use GLOBAL knots for all regions
    # =========================================================================
    # In two-stage DLNM meta-analysis, coefficients must correspond to the 
    # SAME basis specification across regions to be poolable.
    # Using region-specific knots means coefficients aren't comparable.
    # See Gasparrini et al. methodology - same basis everywhere.
    # =========================================================================
    
    # Region-specific temperature range (for MMT search bounds)
    reg_temp_range <- quantile(daily$tmean, probs = c(0.01, 0.99), na.rm = TRUE)
    
    # Create cross-basis with GLOBAL knots and boundaries (same for all regions)
    cb <- crossbasis(
      daily$tmean,
      lag = MAX_LAG,
      argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),  # GLOBAL knots
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
    
    # =========================================================================
    # FIX #2: Search for MMT only within region's OBSERVED temperature range
    # =========================================================================
    # Extrapolation + splines at tails = goofy minima and huge RR
    # Use region's 1st-99th percentile to avoid extrapolation
    # =========================================================================
    temp_seq <- seq(reg_temp_range[1], reg_temp_range[2], length.out = 100)
    cp <- crosspred(cb, model, at = temp_seq, cumul = TRUE, cen = median(daily$tmean))
    mmt_idx <- which.min(cp$allRRfit)
    mmt <- temp_seq[mmt_idx]
    
    # Get RRs at percentiles (centered at MMT)
    # Also compute number of days at each percentile for quality gate
    reg_pcts_vals <- quantile(daily$tmean, probs = c(0.01, 0.05, 0.50, 0.95, 0.99), na.rm = TRUE)
    cp_pcts <- crosspred(cb, model, at = reg_pcts_vals, cumul = TRUE, cen = mmt)
    
    # =========================================================================
    # FIX #3: Data support gate for tail summaries
    # =========================================================================
    # Count days near each percentile threshold (within 1°C)
    n_days_p1 <- sum(abs(daily$tmean - reg_pcts_vals[1]) <= 1, na.rm = TRUE)
    n_days_p99 <- sum(abs(daily$tmean - reg_pcts_vals[5]) <= 1, na.rm = TRUE)
    
    # Quality checks for MVMeta inclusion
    rr_p99 <- cp_pcts$allRRfit[5]
    rr_p1 <- cp_pcts$allRRfit[1]
    vcov_ok <- all(diag(cb_vcov) > 0) && all(is.finite(cb_vcov))
    include_mvmeta <- vcov_ok  # Include all regions with valid vcov
    
    list(
      status = "success",
      region_code = as.character(reg),
      n_days = n_days,
      n_deaths = sum(daily$deaths),
      mmt = mmt,
      mmt_search_range = as.list(reg_temp_range),  # Range used for MMT search
      temp_pcts = as.list(reg_pcts_vals),
      temp_knots_used = "global",  # Flag that we used global knots
      rr_p99 = rr_p99,
      rr_p99_lo = cp_pcts$allRRlow[5],
      rr_p99_hi = cp_pcts$allRRhigh[5],
      rr_p1 = cp_pcts$allRRfit[1],
      rr_p1_lo = cp_pcts$allRRlow[1],
      rr_p1_hi = cp_pcts$allRRhigh[1],
      # Data support info for tail estimates
      n_days_near_p1 = n_days_p1,    # Days within 1°C of P1
      n_days_near_p99 = n_days_p99,  # Days within 1°C of P99
      cb_coef = as.vector(cb_coef),
      cb_vcov = as.vector(cb_vcov),  # Flatten to vector for JSON
      n_params = length(cb_coef),
      include_mvmeta = include_mvmeta
    )
    
  }, error = function(e) {
    list(status = "error", reason = conditionMessage(e), region_code = as.character(reg))
  })
  
  all_results[[i]] <- result
  
  if (i <= 5) cat("  Status:", result$status, "\n")  # Debug
  
  if (!is.null(result$status) && result$status == "success") {
    success_count <- success_count + 1
    
    if (result$include_mvmeta) {
      mvmeta_count <- mvmeta_count + 1
      mvmeta_coefs[[mvmeta_count]] <- result$cb_coef
      # Reconstruct vcov matrix
      n_params <- result$n_params
      mvmeta_vcovs[[mvmeta_count]] <- matrix(result$cb_vcov, nrow = n_params, ncol = n_params)
      mvmeta_regions <- c(mvmeta_regions, result$region_code)
    }
  }
}

cat("\n=======================================================\n")
cat("First-stage results:\n")
cat("  Successful fits:", success_count, "\n")
cat("  Valid for MVMeta:", mvmeta_count, "\n")
cat("=======================================================\n")

# -----------------------------------------------------------------------------
# 4. MVMeta Pooling
# -----------------------------------------------------------------------------
pooled_result <- NULL

if (mvmeta_count >= 10) {
  cat("\nRunning MVMeta pooling with", mvmeta_count, "regions...\n")
  
  # Stack coefficients
  coef_matrix <- do.call(rbind, mvmeta_coefs)
  rownames(coef_matrix) <- mvmeta_regions
  
  # Regularize vcov matrices (add small ridge to diagonal)
  for (j in seq_along(mvmeta_vcovs)) {
    v <- mvmeta_vcovs[[j]]
    diag(v) <- diag(v) + 1e-6 * mean(diag(v))
    mvmeta_vcovs[[j]] <- v
  }
  
  mv_fit <- NULL
  mv_method <- NULL
  
  # Use REML (recommended), fallback to ML if needed
  cat("  Trying method: reml ... ")
  mv_fit <- tryCatch({
    mixmeta(coef_matrix, mvmeta_vcovs, method = "reml", 
            control = list(maxiter = 500, showiter = TRUE))
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
    NULL
  })
  
  if (!is.null(mv_fit)) {
    mv_method <- "reml"
    cat("SUCCESS\n")
  } else {
    cat("FAILED\n")
    cat("  Trying fallback method: ml ... ")
    mv_fit <- tryCatch({
      mixmeta(coef_matrix, mvmeta_vcovs, method = "ml", 
              control = list(maxiter = 500, showiter = TRUE))
    }, error = function(e) {
      cat("Error:", conditionMessage(e), "\n")
      NULL
    })
    
    if (!is.null(mv_fit)) {
      mv_method <- "ml"
      cat("SUCCESS\n")
    } else {
      cat("FAILED\n")
    }
  }
  
  if (!is.null(mv_fit)) {
    pooled_coef <- coef(mv_fit)
    pooled_vcov <- vcov(mv_fit)
    
    cat("\n  MVMeta method:", mv_method, "\n")
    cat("  Converged:", mv_fit$converged, "\n")
    cat("  Pooled coef range:", round(min(pooled_coef), 4), "to", round(max(pooled_coef), 4), "\n")
    
    # Get median MMT from valid regions
    mmts <- sapply(all_results[sapply(all_results, function(x) x$status == "success" && x$include_mvmeta)], function(x) x$mmt)
    pooled_mmt <- median(mmts, na.rm = TRUE)
    
    # Create pooled predictions
    temp_seq <- seq(temp_boundary[1], temp_boundary[2], length.out = 100)
    
    # Create prediction cross-basis
    cb_pred <- crossbasis(
      temp_seq,
      lag = MAX_LAG,
      argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
      arglag = list(fun = "ns", df = LAG_DF)
    )
    
    # Pooled crosspred
    cp_pooled <- crosspred(
      cb_pred,
      coef = pooled_coef,
      vcov = pooled_vcov,
      model.link = "log",
      at = temp_seq,
      cen = pooled_mmt
    )
    
    # Get pooled RRs at percentiles
    pooled_pcts <- list(
      p1 = median(sapply(all_results[sapply(all_results, function(x) x$status == "success")], function(x) x$temp_pcts[[1]]), na.rm = TRUE),
      p5 = median(sapply(all_results[sapply(all_results, function(x) x$status == "success")], function(x) x$temp_pcts[[2]]), na.rm = TRUE),
      p50 = median(sapply(all_results[sapply(all_results, function(x) x$status == "success")], function(x) x$temp_pcts[[3]]), na.rm = TRUE),
      p95 = median(sapply(all_results[sapply(all_results, function(x) x$status == "success")], function(x) x$temp_pcts[[4]]), na.rm = TRUE),
      p99 = median(sapply(all_results[sapply(all_results, function(x) x$status == "success")], function(x) x$temp_pcts[[5]]), na.rm = TRUE)
    )
    
    cp_pooled_pcts <- crosspred(
      cb_pred,
      coef = pooled_coef,
      vcov = pooled_vcov,
      model.link = "log",
      at = c(pooled_pcts$p1, pooled_pcts$p5, pooled_pcts$p50, pooled_pcts$p95, pooled_pcts$p99),
      cen = pooled_mmt
    )
    
    # Extract heterogeneity statistics from mixmeta
    # For multivariate meta-analysis, we need to compute Q properly
    cat("\n  Extracting heterogeneity statistics...\n")
    
    # Method 1: Try qtest from model
    qstat_from_qtest <- tryCatch({
      q <- qtest(mv_fit)
      # qtest returns vectors: Q, df, pvalue for each parameter + overall
      # The first element is the overall test statistic
      overall_Q <- if ("Q" %in% names(q) && length(q$Q) >= 1) q$Q[1] else NA
      overall_df <- if ("df" %in% names(q) && length(q$df) >= 1) q$df[1] else NA
      overall_pvalue <- if ("pvalue" %in% names(q) && length(q$pvalue) >= 1) q$pvalue[1] else NA
      
      list(Q = overall_Q, df = overall_df, pvalue = overall_pvalue)
    }, error = function(e) NULL)
    
    # Method 2: Compute Q manually if qtest fails
    if (is.null(qstat_from_qtest) || is.na(qstat_from_qtest$Q)) {
      qstat_from_qtest <- tryCatch({
        # Residuals approach: Q = sum of squared standardized residuals
        fitted_vals <- fitted(mv_fit)
        residuals_mat <- coef_matrix - fitted_vals
        
        # Total Q statistic (sum across all studies and parameters)
        Q_total <- 0
        for (j in seq_len(nrow(residuals_mat))) {
          resid_j <- residuals_mat[j, ]
          vcov_j <- mvmeta_vcovs[[j]]
          # Add small regularization to avoid singularity
          vcov_j_reg <- vcov_j + diag(1e-6 * mean(diag(vcov_j)), nrow(vcov_j))
          Q_j <- tryCatch({
            t(resid_j) %*% solve(vcov_j_reg) %*% resid_j
          }, error = function(e) NA)
          if (!is.na(Q_j)) Q_total <- Q_total + as.numeric(Q_j)
        }
        
        # Degrees of freedom: (n_studies - 1) * n_params for multivariate
        n_studies <- nrow(coef_matrix)
        n_params <- ncol(coef_matrix)
        Q_df <- (n_studies - 1) * n_params
        Q_pvalue <- pchisq(Q_total, df = Q_df, lower.tail = FALSE)
        
        list(Q = Q_total, df = Q_df, pvalue = Q_pvalue)
      }, error = function(e) {
        list(Q = NA, df = NA, pvalue = NA)
      })
    }
    
    qstat_info <- qstat_from_qtest
    
    # Calculate I² (proportion of total variability due to heterogeneity)
    # I² = max(0, (Q - df) / Q) * 100
    I2 <- if(!is.na(qstat_info$Q) && !is.na(qstat_info$df) && qstat_info$df > 0 && qstat_info$Q > 0) {
      max(0, (qstat_info$Q - qstat_info$df) / qstat_info$Q * 100)
    } else {
      NA
    }
    
    # Calculate H² (relative excess in Q over df)
    # H² = Q / df
    H2 <- if(!is.na(qstat_info$Q) && !is.na(qstat_info$df) && qstat_info$df > 0) {
      qstat_info$Q / qstat_info$df
    } else {
      NA
    }
    
    # Extract between-study variance (tau²) from mixmeta
    # For multivariate, this is a matrix; we summarize with trace or diagonal
    tau2_info <- tryCatch({
      psi <- mv_fit$Psi  # Between-study covariance matrix
      if (!is.null(psi) && is.matrix(psi)) {
        list(
          tau2_trace = sum(diag(psi)),  # Total between-study variance
          tau2_diag = diag(psi),         # Diagonal elements
          tau2_mean = mean(diag(psi))    # Average variance per parameter
        )
      } else if (!is.null(psi)) {
        list(tau2_trace = sum(psi), tau2_diag = psi, tau2_mean = mean(psi))
      } else {
        list(tau2_trace = NA, tau2_diag = NA, tau2_mean = NA)
      }
    }, error = function(e) {
      list(tau2_trace = NA, tau2_diag = NA, tau2_mean = NA)
    })
    
    cat("\n  Heterogeneity statistics:\n")
    cat("    Cochran's Q:", round(qstat_info$Q, 2), "on", qstat_info$df, "df (p =", format.pval(qstat_info$pvalue, digits = 3), ")\n")
    cat("    I² (approx):", round(I2, 1), "%\n")
    cat("    H² (Q/df):", round(H2, 2), "\n")
    cat("    τ² (trace):", round(tau2_info$tau2_trace, 4), "\n")
    
    pooled_result <- list(
      method = mv_method,
      converged = mv_fit$converged,
      n_regions = mvmeta_count,
      mmt = pooled_mmt,
      # Heterogeneity statistics (comprehensive)
      heterogeneity = list(
        cochrans_Q = qstat_info$Q,
        Q_df = qstat_info$df,
        Q_pvalue = qstat_info$pvalue,
        I2_percent = I2,
        H2 = H2,
        tau2_trace = tau2_info$tau2_trace,
        tau2_mean = tau2_info$tau2_mean
      ),
      pooled_coef = as.vector(pooled_coef),
      pooled_vcov = as.vector(pooled_vcov),
      n_params = length(pooled_coef),
      rr_heat_p99 = list(
        rr = cp_pooled_pcts$allRRfit[5],
        rr_lo = cp_pooled_pcts$allRRlow[5],
        rr_hi = cp_pooled_pcts$allRRhigh[5],
        temp = pooled_pcts$p99
      ),
      rr_heat_p95 = list(
        rr = cp_pooled_pcts$allRRfit[4],
        rr_lo = cp_pooled_pcts$allRRlow[4],
        rr_hi = cp_pooled_pcts$allRRhigh[4],
        temp = pooled_pcts$p95
      ),
      rr_cold_p1 = list(
        rr = cp_pooled_pcts$allRRfit[1],
        rr_lo = cp_pooled_pcts$allRRlow[1],
        rr_hi = cp_pooled_pcts$allRRhigh[1],
        temp = pooled_pcts$p1
      ),
      rr_cold_p5 = list(
        rr = cp_pooled_pcts$allRRfit[2],
        rr_lo = cp_pooled_pcts$allRRlow[2],
        rr_hi = cp_pooled_pcts$allRRhigh[2],
        temp = pooled_pcts$p5
      ),
      temp_curve = list(
        temp = temp_seq,
        rr = as.vector(cp_pooled$allRRfit),
        rr_lo = as.vector(cp_pooled$allRRlow),
        rr_hi = as.vector(cp_pooled$allRRhigh)
      )
    )
    
    cat("\n=== POOLED RESULTS ===\n")
    cat("Heat (P99):", round(pooled_result$rr_heat_p99$rr, 3), 
        "(", round(pooled_result$rr_heat_p99$rr_lo, 3), "-", 
        round(pooled_result$rr_heat_p99$rr_hi, 3), ")\n")
    cat("Cold (P1):", round(pooled_result$rr_cold_p1$rr, 3),
        "(", round(pooled_result$rr_cold_p1$rr_lo, 3), "-",
        round(pooled_result$rr_cold_p1$rr_hi, 3), ")\n")
    cat("MMT:", round(pooled_mmt, 1), "°C\n")
    cat("\n=== HETEROGENEITY ===\n")
    cat("Cochran's Q:", round(pooled_result$heterogeneity$cochrans_Q, 2), 
        "on", pooled_result$heterogeneity$Q_df, "df\n")
    cat("Q p-value:", format.pval(pooled_result$heterogeneity$Q_pvalue, digits = 3), "\n")
    cat("I² (approx):", round(pooled_result$heterogeneity$I2_percent, 1), "%\n")
    cat("H² (Q/df):", round(pooled_result$heterogeneity$H2, 2), "\n")
    cat("τ² (trace):", round(pooled_result$heterogeneity$tau2_trace, 4), "\n")
    
  } else {
    cat("\n  All MVMeta methods failed.\n")
  }
} else {
  cat("\n  Insufficient regions for MVMeta (need 10, have", mvmeta_count, ")\n")
}

# -----------------------------------------------------------------------------
# 5. Save Results
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("Saving results...\n")
cat("=======================================================\n")

# Filter to only successful results for saving
success_results <- Filter(function(x) x$status == "success", all_results)

# Convert to simpler list for JSON
region_list <- lapply(success_results, function(x) {
  list(
    region_code = x$region_code,
    n_days = x$n_days,
    n_deaths = x$n_deaths,
    mmt = x$mmt,
    rr_p99 = x$rr_p99,
    rr_p99_lo = x$rr_p99_lo,
    rr_p99_hi = x$rr_p99_hi,
    rr_p1 = x$rr_p1,
    rr_p1_lo = x$rr_p1_lo,
    rr_p1_hi = x$rr_p1_hi,
    include_mvmeta = x$include_mvmeta,
    cb_coef = x$cb_coef,
    cb_vcov = x$cb_vcov,
    n_params = x$n_params
  )
})

output <- list(
  exposure_type = EXPOSURE_TYPE,
  analysis_date = as.character(Sys.time()),
  dlnm_params = list(
    max_lag = MAX_LAG,
    temp_df = TEMP_DF,
    lag_df = LAG_DF,
    temp_boundary = as.vector(temp_boundary),
    temp_knots = as.vector(temp_pcts)
  ),
  n_regions_total = n_regions,
  n_regions_success = success_count,
  n_regions_mvmeta = mvmeta_count,
  region_results = region_list,
  pooled = pooled_result
)

output_file <- file.path(OUTPUT_DIR, paste0("dlnm_r_", EXPOSURE_TYPE, "_results_v2.json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 8)
cat("Results saved to:", output_file, "\n")

# Save summary CSV
if (length(success_results) > 0) {
  summary_df <- data.frame(
    region_code = sapply(success_results, function(x) x$region_code),
    n_days = sapply(success_results, function(x) x$n_days),
    n_deaths = sapply(success_results, function(x) x$n_deaths),
    mmt = sapply(success_results, function(x) x$mmt),
    rr_p99 = sapply(success_results, function(x) x$rr_p99),
    rr_p99_lo = sapply(success_results, function(x) x$rr_p99_lo),
    rr_p99_hi = sapply(success_results, function(x) x$rr_p99_hi),
    rr_p1 = sapply(success_results, function(x) x$rr_p1),
    rr_p1_lo = sapply(success_results, function(x) x$rr_p1_lo),
    rr_p1_hi = sapply(success_results, function(x) x$rr_p1_hi),
    in_mvmeta = sapply(success_results, function(x) x$include_mvmeta)
  )
  
  summary_file <- file.path(OUTPUT_DIR, paste0("dlnm_r_", EXPOSURE_TYPE, "_summary.csv"))
  fwrite(summary_df, summary_file)
  cat("Summary saved to:", summary_file, "\n")
}

cat("\nDone!\n")
