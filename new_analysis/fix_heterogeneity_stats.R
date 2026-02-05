#!/usr/bin/env Rscript
# fix_heterogeneity_stats.R
# Adds heterogeneity statistics to existing JSON results without re-running the full analysis
# Uses qtest() from mixmeta package to compute Cochran's Q and I²

library(jsonlite)
library(mixmeta)

cat("\n========================================\n")
cat("Adding Heterogeneity Statistics to Results\n")
cat("========================================\n\n")

# Get script directory
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("--file=", "", file_arg)))
} else {
  script_dir <- getwd()
}

# Function to compute Q-test from stored coefficients and vcov
compute_heterogeneity_from_regions <- function(region_list) {
  # Extract valid region results (those with include_mvmeta == TRUE)
  valid_regions <- Filter(function(r) isTRUE(r$include_mvmeta), region_list)
  n_regions <- length(valid_regions)
  
  if (n_regions < 2) {
    return(list(
      cochrans_Q = NA,
      Q_df = NA,
      Q_pvalue = NA,
      I2_percent = NA,
      n_regions = n_regions
    ))
  }
  
  # Get number of parameters from first region
  n_params <- valid_regions[[1]]$n_params
  
  # Extract coefficients - each cb_coef is a list of 16 values
  coef_matrix <- do.call(rbind, lapply(valid_regions, function(r) {
    unlist(r$cb_coef)
  }))
  
  # Extract vcov matrices - each cb_vcov is a list of 256 values
  S_list <- lapply(valid_regions, function(r) {
    v <- unlist(r$cb_vcov)
    matrix(v, nrow = n_params, ncol = n_params)
  })
  
  # Fit mixmeta to get Q-test
  tryCatch({
    # Fit model directly with matrix (simpler approach)
    mv_fit <- mixmeta(
      coef_matrix,
      S = S_list,
      method = "reml",
      control = list(maxiter = 500)
    )
    
    # Get Q-test
    qt <- qtest(mv_fit)
    
    # Extract overall Q statistic (first element named ".all")
    Q <- as.numeric(qt$Q[".all"])
    df <- as.numeric(qt$df[".all"])
    pvalue <- as.numeric(qt$pvalue[".all"])
    
    # Calculate I²
    I2 <- if (!is.na(Q) && !is.na(df) && df > 0 && Q > df) {
      (Q - df) / Q * 100
    } else {
      0
    }
    
    return(list(
      cochrans_Q = Q,
      Q_df = df,
      Q_pvalue = pvalue,
      I2_percent = I2,
      n_regions = n_regions
    ))
    
  }, error = function(e) {
    cat("    Error computing Q-test:", conditionMessage(e), "\n")
    return(list(
      cochrans_Q = NA,
      Q_df = NA,
      Q_pvalue = NA,
      I2_percent = NA,
      n_regions = n_regions
    ))
  })
}

# Process DLNM results for both levels
for (level in c("intermediate", "immediate")) {
  cat("Processing", toupper(level), "level...\n")
  
  json_path <- file.path(script_dir, "phase1_r", "results", 
                         paste0("dlnm_r_", level, "_results_v2.json"))
  
  if (!file.exists(json_path)) {
    cat("  File not found:", json_path, "\n")
    next
  }
  
  # Read JSON
  d <- fromJSON(json_path, simplifyVector = FALSE)
  
  # Get region results list directly
  region_list <- d$region_results
  
  n_valid <- sum(sapply(region_list, function(r) isTRUE(r$include_mvmeta)))
  cat("  Found", n_valid, "regions with valid coefficients\n")
  
  # Compute heterogeneity
  het_stats <- compute_heterogeneity_from_regions(region_list)
  
  cat("  Cochran's Q:", round(het_stats$cochrans_Q, 2), "\n")
  cat("  Q df:", het_stats$Q_df, "\n")
  cat("  Q p-value:", format.pval(het_stats$Q_pvalue, digits = 3), "\n")
  cat("  I²:", round(het_stats$I2_percent, 1), "%\n")
  
  # Update the JSON
  d$pooled$heterogeneity <- list(
    cochrans_Q = het_stats$cochrans_Q,
    Q_df = het_stats$Q_df,
    Q_pvalue = het_stats$Q_pvalue,
    I2_percent = het_stats$I2_percent
  )
  
  # Write back
  write_json(d, json_path, auto_unbox = TRUE, digits = 8, pretty = TRUE)
  cat("  Updated:", json_path, "\n\n")
}

cat("========================================\n")
cat("Heterogeneity statistics added successfully!\n")
cat("========================================\n")
