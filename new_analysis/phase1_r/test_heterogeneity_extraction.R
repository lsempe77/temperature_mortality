# =============================================================================
# test_heterogeneity_extraction.R
# Quick test of heterogeneity extraction on a small sample
# =============================================================================

suppressPackageStartupMessages({
  library(dlnm)
  library(mixmeta)
  library(data.table)
  library(arrow)
  library(jsonlite)
  library(splines)
})

cat("=======================================================\n")
cat("TESTING HETEROGENEITY EXTRACTION (Small Sample)\n")
cat("=======================================================\n")

# -----------------------------------------------------------------------------
# 1. Load Data (subset to 20 regions for speed)
# -----------------------------------------------------------------------------
SCRIPT_DIR <- getwd()
DATA_DIR <- file.path(SCRIPT_DIR, "..", "phase0_data_prep", "results")

cat("[1] Loading data (small sample)...\n")

mort <- as.data.table(read_parquet(file.path(DATA_DIR, "mortality_intermediate_daily_elderly.parquet")))
era5 <- as.data.table(read_parquet(file.path(DATA_DIR, "era5_intermediate_daily.parquet")))

# Standardize column names
if ("intermediate_code" %in% names(mort)) setnames(mort, "intermediate_code", "region_code")

# Merge
mort[, date := as.character(as.Date(date))]
era5[, date := as.character(date)]
setkey(mort, region_code, date)
setkey(era5, region_code, date)
data <- merge(mort, era5, by = c("region_code", "date"))
data[, date := as.Date(date)]
data[, deaths := deaths_elderly]
data[, tmean := temp_mean]
data <- data[!is.na(deaths) & !is.na(tmean)]

# Take only 20 regions with most data
region_counts <- data[, .N, by = region_code][order(-N)]
test_regions <- region_counts[1:20, region_code]
data <- data[region_code %in% test_regions]

cat("  Testing with", length(test_regions), "regions\n")
cat("  Total rows:", nrow(data), "\n")

# Global percentiles
temp_all <- data$tmean
temp_pcts <- quantile(temp_all, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)
temp_boundary <- c(0, 40)

cat("  Global knots:", round(temp_pcts, 1), "\n")

# -----------------------------------------------------------------------------
# 2. Fit Region-Specific DLNMs
# -----------------------------------------------------------------------------
cat("\n[2] Fitting DLNMs for", length(test_regions), "regions...\n")

MAX_LAG <- 21
LAG_DF <- 4

mvmeta_coefs <- list()
mvmeta_vcovs <- list()
mvmeta_regions <- c()

for (i in seq_along(test_regions)) {
  reg <- test_regions[i]
  daily <- data[region_code == reg][order(date)]
  
  if (nrow(daily) < 365) next
  
  result <- tryCatch({
    reg_p1 <- quantile(daily$tmean, 0.01, na.rm = TRUE)
    reg_p99 <- quantile(daily$tmean, 0.99, na.rm = TRUE)
    
    cb <- crossbasis(
      daily$tmean,
      lag = MAX_LAG,
      argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
      arglag = list(fun = "ns", df = LAG_DF)
    )
    
    model <- glm(
      deaths ~ cb + ns(as.numeric(date), df = 7 * length(unique(format(daily$date, "%Y")))) +
        factor(format(date, "%u")),
      data = daily,
      family = quasipoisson(link = "log"),
      na.action = na.exclude
    )
    
    cb_idx <- grep("cb", names(coef(model)))
    cb_coef <- coef(model)[cb_idx]
    cb_vcov <- vcov(model)[cb_idx, cb_idx]
    
    vcov_ok <- all(diag(cb_vcov) > 0) && all(is.finite(cb_vcov))
    
    if (vcov_ok) {
      mvmeta_coefs[[length(mvmeta_coefs) + 1]] <- as.vector(cb_coef)
      mvmeta_vcovs[[length(mvmeta_vcovs) + 1]] <- cb_vcov
      mvmeta_regions <- c(mvmeta_regions, as.character(reg))
    }
    
    "success"
  }, error = function(e) "error")
  
  cat("  Region", i, ":", reg, "-", result, "\n")
}

cat("\n  Successful fits:", length(mvmeta_regions), "\n")

# -----------------------------------------------------------------------------
# 3. Run MVMeta and Extract Heterogeneity
# -----------------------------------------------------------------------------
cat("\n[3] Running MVMeta pooling...\n")

coef_matrix <- do.call(rbind, mvmeta_coefs)
rownames(coef_matrix) <- mvmeta_regions

# Regularize vcov
for (j in seq_along(mvmeta_vcovs)) {
  v <- mvmeta_vcovs[[j]]
  diag(v) <- diag(v) + 1e-6 * mean(diag(v))
  mvmeta_vcovs[[j]] <- v
}

# Fit
mv_fit <- mixmeta(coef_matrix, mvmeta_vcovs, method = "reml",
                  control = list(maxiter = 500, showiter = TRUE))

cat("\n  Converged:", mv_fit$converged, "\n")
cat("  Method:", mv_fit$method, "\n")

# -----------------------------------------------------------------------------
# 4. Extract Heterogeneity Statistics
# -----------------------------------------------------------------------------
cat("\n[4] Extracting heterogeneity statistics...\n")

# Method 1: Try qtest()
qstat_from_qtest <- tryCatch({
  q <- qtest(mv_fit)
  cat("  qtest() output:\n")
  print(q)
  # qtest returns: Q (vector per param + overall), df, pvalue
  # The OVERALL Q is in q$Q[1] or we sum individual Qs
  # Actually the structure is: Q = vector of length (k+1) where first is overall
  # Let's extract the overall values properly
  if (is.list(q)) {
    # Get overall Q (first element or named "Overall")
    overall_Q <- if ("Q" %in% names(q)) {
      if (length(q$Q) > 1) q$Q[1] else q$Q  # First element is overall
    } else NA
    overall_df <- if ("df" %in% names(q)) {
      if (length(q$df) > 1) q$df[1] else q$df
    } else NA
    overall_pvalue <- if ("pvalue" %in% names(q)) {
      if (length(q$pvalue) > 1) q$pvalue[1] else q$pvalue
    } else NA
    
    cat("  Extracted overall Q:", overall_Q, "\n")
    cat("  Extracted overall df:", overall_df, "\n")
    cat("  Extracted overall p:", overall_pvalue, "\n")
    
    list(Q = overall_Q, df = overall_df, pvalue = overall_pvalue)
  } else {
    NULL
  }
}, error = function(e) {
  cat("  qtest() failed:", e$message, "\n")
  NULL
})

# Method 2: Manual computation
cat("\n  Manual Q computation:\n")
fitted_vals <- fitted(mv_fit)
residuals_mat <- coef_matrix - fitted_vals

Q_total <- 0
for (j in seq_len(nrow(residuals_mat))) {
  resid_j <- residuals_mat[j, ]
  vcov_j <- mvmeta_vcovs[[j]]
  vcov_j_reg <- vcov_j + diag(1e-6 * mean(diag(vcov_j)), nrow(vcov_j))
  Q_j <- tryCatch({
    t(resid_j) %*% solve(vcov_j_reg) %*% resid_j
  }, error = function(e) NA)
  if (!is.na(Q_j)) Q_total <- Q_total + as.numeric(Q_j)
}

n_studies <- nrow(coef_matrix)
n_params <- ncol(coef_matrix)
Q_df <- (n_studies - 1) * n_params
Q_pvalue <- pchisq(Q_total, df = Q_df, lower.tail = FALSE)

cat("    Manual Q:", round(Q_total, 2), "\n")
cat("    Manual df:", Q_df, "\n")
cat("    Manual p-value:", format.pval(Q_pvalue, digits = 3), "\n")

# Use qtest if available, else manual
if (!is.null(qstat_from_qtest)) {
  qstat_info <- qstat_from_qtest
  cat("\n  Using qtest() results\n")
} else {
  qstat_info <- list(Q = Q_total, df = Q_df, pvalue = Q_pvalue)
  cat("\n  Using manual computation\n")
}

# Calculate I²
I2 <- if(!is.na(qstat_info$Q) && !is.na(qstat_info$df) && qstat_info$df > 0 && qstat_info$Q > 0) {
  max(0, (qstat_info$Q - qstat_info$df) / qstat_info$Q * 100)
} else {
  NA
}

# Calculate H²
H2 <- if(!is.na(qstat_info$Q) && !is.na(qstat_info$df) && qstat_info$df > 0) {
  qstat_info$Q / qstat_info$df
} else {
  NA
}

# Extract tau² from Psi matrix
cat("\n  Extracting tau² from Psi matrix...\n")
psi <- mv_fit$Psi
if (!is.null(psi)) {
  cat("    Psi dimensions:", dim(psi), "\n")
  tau2_trace <- sum(diag(psi))
  tau2_mean <- mean(diag(psi))
  cat("    tau² (trace):", round(tau2_trace, 6), "\n")
  cat("    tau² (mean per param):", round(tau2_mean, 6), "\n")
} else {
  tau2_trace <- NA
  tau2_mean <- NA
  cat("    Psi is NULL\n")
}

# -----------------------------------------------------------------------------
# 5. Summary
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("HETEROGENEITY SUMMARY\n")
cat("=======================================================\n")
cat("Cochran's Q:", round(qstat_info$Q, 2), "on", qstat_info$df, "df\n")
cat("Q p-value:", format.pval(qstat_info$pvalue, digits = 3), "\n")
cat("I² (%):", round(I2, 1), "\n")
cat("H² (Q/df):", round(H2, 2), "\n")
cat("τ² (trace):", round(tau2_trace, 6), "\n")
cat("τ² (mean):", round(tau2_mean, 6), "\n")

# Final structure
heterogeneity <- list(
  cochrans_Q = qstat_info$Q,
  Q_df = qstat_info$df,
  Q_pvalue = qstat_info$pvalue,
  I2_percent = I2,
  H2 = H2,
  tau2_trace = tau2_trace,
  tau2_mean = tau2_mean
)

cat("\n  JSON output:\n")
cat(toJSON(heterogeneity, auto_unbox = TRUE, pretty = TRUE), "\n")

cat("\n=======================================================\n")
cat("Test complete!\n")
cat("=======================================================\n")
