# MVMeta pooling using Gasparrini's R packages
# This is the proper implementation following the 2015 Lancet methodology

library(dlnm)
library(mvmeta)
library(jsonlite)
library(splines)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript mvmeta_pooling.R <input_json> <output_json>")
}

input_file <- args[1]
output_file <- args[2]

cat("Loading data from:", input_file, "\n")

# Load the Phase 1 results
data <- fromJSON(input_file)
region_results <- data$region_results

# Extract coefficients and vcov matrices
coef_list <- list()
vcov_list <- list()
region_codes <- c()
temp_pcts <- list()

for (code in names(region_results)) {
  r <- region_results[[code]]
  if (!is.null(r$cb_coefs) && !is.null(r$cb_vcov)) {
    coef_list[[length(coef_list) + 1]] <- as.numeric(r$cb_coefs)
    vcov_list[[length(vcov_list) + 1]] <- matrix(unlist(r$cb_vcov), 
                                                   nrow = sqrt(length(unlist(r$cb_vcov))))
    region_codes <- c(region_codes, code)
    temp_pcts[[length(temp_pcts) + 1]] <- r$temp_percentiles
  }
}

n_regions <- length(coef_list)
n_params <- length(coef_list[[1]])
cat("Loaded", n_regions, "regions with", n_params, "parameters each\n")

# Get crossbasis info from first region
cb_info <- region_results[[region_codes[1]]]$crossbasis_info
temp_knots <- as.numeric(cb_info$temp_knots)
temp_boundary <- as.numeric(cb_info$temp_boundary)
lag_basis <- matrix(unlist(cb_info$lag_basis), nrow = cb_info$max_lag + 1)
max_lag <- cb_info$max_lag
temp_df <- cb_info$temp_df
lag_df <- cb_info$lag_df

cat("Cross-basis: temp_df=", temp_df, ", lag_df=", lag_df, ", max_lag=", max_lag, "\n")

# Stack coefficients and vcov into matrices for mvmeta
coef_matrix <- do.call(rbind, coef_list)
rownames(coef_matrix) <- region_codes

# Run multivariate meta-analysis using mvmeta
cat("\nRunning multivariate meta-analysis...\n")

# mvmeta expects a list of vcov matrices
mv_fit <- mvmeta(coef_matrix, vcov_list, method = "reml")

cat("MVMeta converged:", mv_fit$converged, "\n")
cat("I2 (heterogeneity):", round(mv_fit$I2, 1), "%\n")

# Extract pooled coefficients and vcov
pooled_coef <- coef(mv_fit)
pooled_vcov <- vcov(mv_fit)

cat("Pooled coefficient range:", round(min(pooled_coef), 6), "to", round(max(pooled_coef), 6), "\n")

# Create crossbasis for prediction using the same specification
# We need to recreate the basis at specific temperatures

# Define a function to compute the temperature spline basis
compute_temp_basis <- function(temp, knots, boundary) {
  # Natural cubic spline with specified knots
  ns(temp, knots = knots, Boundary.knots = boundary)
}

# Compute lag basis sum (cumulative effect)
lag_sum <- colSums(lag_basis)

# Function to predict RR at a temperature vs reference
predict_rr <- function(target_temp, ref_temp, pooled_coef, pooled_vcov, 
                       temp_knots, temp_boundary, lag_sum) {
  
  # Compute temperature basis at target and reference
  temp_basis_target <- compute_temp_basis(target_temp, temp_knots, temp_boundary)
  temp_basis_ref <- compute_temp_basis(ref_temp, temp_knots, temp_boundary)
  temp_diff <- as.numeric(temp_basis_target - temp_basis_ref)
  
  # Build contrast vector (tensor product)
  # Column order: for each temp basis column, all lag basis columns
  K_var <- length(temp_diff)
  K_lag <- length(lag_sum)
  
  contrast <- numeric(K_var * K_lag)
  for (t_idx in 1:K_var) {
    for (l_idx in 1:K_lag) {
      col_idx <- (t_idx - 1) * K_lag + l_idx
      contrast[col_idx] <- temp_diff[t_idx] * lag_sum[l_idx]
    }
  }
  
  # Log-RR and SE via delta method
  log_rr <- sum(contrast * pooled_coef)
  var_log_rr <- t(contrast) %*% pooled_vcov %*% contrast
  se_log_rr <- sqrt(max(0, var_log_rr))
  
  # RR and CI
  rr <- exp(log_rr)
  rr_lo <- exp(log_rr - 1.96 * se_log_rr)
  rr_hi <- exp(log_rr + 1.96 * se_log_rr)
  
  return(list(
    rr = rr,
    rr_lo = rr_lo,
    rr_hi = rr_hi,
    log_rr = log_rr,
    se_log_rr = as.numeric(se_log_rr)
  ))
}

# Get median percentiles across regions
p1_vals <- sapply(temp_pcts, function(x) x$p1)
p5_vals <- sapply(temp_pcts, function(x) x$p5)
p50_vals <- sapply(temp_pcts, function(x) x$p50)
p95_vals <- sapply(temp_pcts, function(x) x$p95)
p99_vals <- sapply(temp_pcts, function(x) x$p99)

p1 <- median(p1_vals, na.rm = TRUE)
p5 <- median(p5_vals, na.rm = TRUE)
p50 <- median(p50_vals, na.rm = TRUE)
p95 <- median(p95_vals, na.rm = TRUE)
p99 <- median(p99_vals, na.rm = TRUE)

cat("\nPercentiles: P1=", round(p1, 1), ", P50=", round(p50, 1), ", P99=", round(p99, 1), "\n")

# Compute RRs at key percentiles
heat_extreme <- predict_rr(p99, p50, pooled_coef, pooled_vcov, temp_knots, temp_boundary, lag_sum)
heat_moderate <- predict_rr(p95, p50, pooled_coef, pooled_vcov, temp_knots, temp_boundary, lag_sum)
cold_extreme <- predict_rr(p1, p50, pooled_coef, pooled_vcov, temp_knots, temp_boundary, lag_sum)
cold_moderate <- predict_rr(p5, p50, pooled_coef, pooled_vcov, temp_knots, temp_boundary, lag_sum)

cat("\n=== Pooled Results ===\n")
cat(sprintf("Heat P99 vs P50: RR=%.4f (%.4f-%.4f)\n", 
            heat_extreme$rr, heat_extreme$rr_lo, heat_extreme$rr_hi))
cat(sprintf("Heat P95 vs P50: RR=%.4f (%.4f-%.4f)\n", 
            heat_moderate$rr, heat_moderate$rr_lo, heat_moderate$rr_hi))
cat(sprintf("Cold P1 vs P50:  RR=%.4f (%.4f-%.4f)\n", 
            cold_extreme$rr, cold_extreme$rr_lo, cold_extreme$rr_hi))
cat(sprintf("Cold P5 vs P50:  RR=%.4f (%.4f-%.4f)\n", 
            cold_moderate$rr, cold_moderate$rr_lo, cold_moderate$rr_hi))

# Build output
output <- list(
  mvmeta = list(
    converged = mv_fit$converged,
    method = "reml",
    n_regions = n_regions,
    n_params = n_params,
    I2 = mv_fit$I2,
    pooled_coef = as.list(pooled_coef),
    pooled_vcov = pooled_vcov
  ),
  percentiles = list(
    p1 = p1, p5 = p5, p50 = p50, p95 = p95, p99 = p99
  ),
  pooled_rr = list(
    heat_p99_vs_p50 = heat_extreme,
    heat_p95_vs_p50 = heat_moderate,
    cold_p1_vs_p50 = cold_extreme,
    cold_p5_vs_p50 = cold_moderate
  ),
  crossbasis_info = list(
    temp_knots = temp_knots,
    temp_boundary = temp_boundary,
    temp_df = temp_df,
    lag_df = lag_df,
    max_lag = max_lag,
    lag_sum = as.list(lag_sum)
  )
)

# Write output
cat("\nWriting results to:", output_file, "\n")
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE)

cat("Done!\n")
