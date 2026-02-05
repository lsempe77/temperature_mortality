# =============================================================================
# Test Script: Debug DLNM in Case-Crossover Framework
# =============================================================================
# Minimal test to understand why crosspred fails with stratum fixed effects

suppressPackageStartupMessages({
  library(data.table)
  library(dlnm)
  library(splines)
})

cat("=======================================================\n")
cat("TEST: DLNM Case-Crossover Debugging\n")
cat("=======================================================\n\n")

# -----------------------------------------------------------------------------
# 1. Create synthetic data (5 strata, ~100 days each)
# -----------------------------------------------------------------------------
cat("[1] Creating synthetic test data...\n")

set.seed(42)
n_strata <- 5
days_per_stratum <- 100
MAX_LAG <- 21

# Generate dates and strata
dates <- seq(as.Date("2020-01-01"), by = "day", length.out = 500)
strata <- rep(paste0("S", 1:n_strata), each = days_per_stratum)

# Generate temperature (smooth seasonal pattern + noise)
day_of_year <- as.numeric(format(dates, "%j"))
temp_base <- 20 + 8 * sin(2 * pi * (day_of_year - 80) / 365)
temp <- temp_base + rnorm(length(dates), 0, 2)

# Create lag matrix
temp_lag_matrix <- matrix(NA, nrow = length(temp), ncol = MAX_LAG + 1)
for (lag in 0:MAX_LAG) {
  if (lag == 0) {
    temp_lag_matrix[, lag + 1] <- temp
  } else {
    temp_lag_matrix[(lag + 1):length(temp), lag + 1] <- temp[1:(length(temp) - lag)]
  }
}
colnames(temp_lag_matrix) <- paste0("lag", 0:MAX_LAG)

# Generate deaths (Poisson with temperature effect)
# True effect: U-shaped with MMT around 22°C
mmt_true <- 22
beta_cold <- 0.02  # log-RR per degree below MMT
beta_heat <- 0.03  # log-RR per degree above MMT
temp_effect <- ifelse(temp < mmt_true, 
                      beta_cold * (mmt_true - temp),
                      beta_heat * (temp - mmt_true))
lambda <- exp(3 + temp_effect)  # baseline ~20 deaths/day
deaths <- rpois(length(dates), lambda)

# Build data.table
daily_data <- data.table(
  date = dates,
  stratum = strata,
  deaths = deaths,
  temp = temp
)
daily_data <- cbind(daily_data, as.data.table(temp_lag_matrix))

# Remove rows with incomplete lags FIRST
daily_data <- daily_data[(MAX_LAG + 1):nrow(daily_data)]

cat("  Rows:", nrow(daily_data), "\n")
cat("  Strata:", length(unique(daily_data$stratum)), "\n")
cat("  Deaths range:", min(daily_data$deaths), "-", max(daily_data$deaths), "\n")
cat("  Temp range:", round(min(daily_data$temp), 1), "-", round(max(daily_data$temp), 1), "\n")

# -----------------------------------------------------------------------------
# 2. Set up DLNM crossbasis - TWO APPROACHES
# -----------------------------------------------------------------------------
cat("\n[2] Creating crossbasis...\n")

lag_cols <- paste0("lag", 0:MAX_LAG)
temp_matrix <- as.matrix(daily_data[, ..lag_cols])

# Temperature percentiles for knots
temp_pcts <- quantile(daily_data$temp, c(0.10, 0.75, 0.90))
temp_boundary <- c(min(daily_data$temp) - 1, max(daily_data$temp) + 1)

# APPROACH 1: Create crossbasis from VECTOR (lets dlnm build lag internally)
# This is how the main DLNM analysis works
cat("\n  Approach 1: Crossbasis from vector (dlnm builds lags)...\n")
cb_vector <- crossbasis(
  daily_data$temp,  # Just the temperature vector
  lag = MAX_LAG,
  argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
  arglag = list(fun = "ns", df = 4)
)
cat("    Crossbasis dim:", dim(cb_vector)[1], "x", dim(cb_vector)[2], "\n")

# APPROACH 2: Create crossbasis from pre-built matrix
cat("\n  Approach 2: Crossbasis from matrix (pre-built lags)...\n")
cb_matrix <- crossbasis(
  temp_matrix,
  lag = MAX_LAG,
  argvar = list(fun = "ns", knots = temp_pcts, Boundary.knots = temp_boundary),
  arglag = list(fun = "ns", df = 4)
)
cat("    Crossbasis dim:", dim(cb_matrix)[1], "x", dim(cb_matrix)[2], "\n")

# Check if they're different
cat("\n  Comparing approaches:\n")
cat("    Vector approach names:", paste(head(colnames(cb_vector), 4), collapse=", "), "...\n")
cat("    Matrix approach names:", paste(head(colnames(cb_matrix), 4), collapse=", "), "...\n")
cat("    Max abs diff:", max(abs(cb_vector - cb_matrix), na.rm=TRUE), "\n")

# -----------------------------------------------------------------------------
# 3. Test different model specifications
# -----------------------------------------------------------------------------
cat("\n[3] Testing model specifications...\n")

# Use cb_vector (from vector approach) for testing
cb <- cb_vector
cb_names <- colnames(cb)

# Model A: Simple Poisson without strata (should work)
cat("\n  Model A: Simple Poisson with VECTOR crossbasis...\n")

# Create a clean data frame with crossbasis columns
model_data <- data.frame(
  deaths = daily_data$deaths,
  stratum = daily_data$stratum
)
# Add crossbasis columns with proper names
for (i in seq_along(cb_names)) {
  model_data[[cb_names[i]]] <- cb[, i]
}

formula_A <- as.formula(paste("deaths ~", paste(cb_names, collapse = " + ")))
fit_A <- glm(formula_A, data = model_data, family = quasipoisson())
cat("    Converged:", fit_A$converged, "\n")
cat("    N obs in model:", nrow(model_data), "\n")
cat("    N rows in cb:", nrow(cb), "\n")

cp_A <- tryCatch({
  crosspred(cb, fit_A, at = seq(10, 30, by = 1), cen = 22)
}, error = function(e) {
  cat("    crosspred error:", conditionMessage(e), "\n")
  NULL
})
if (!is.null(cp_A)) {
  cat("    SUCCESS! RR at 30°C:", round(cp_A$allRRfit["30"], 3), "\n")
  cat("    MMT:", cp_A$predvar[which.min(cp_A$allRRfit)], "°C\n")
}

# Model B: Poisson with stratum fixed effects (using vector crossbasis)
cat("\n  Model B: Poisson + strata with VECTOR crossbasis...\n")
model_data$stratum_factor <- factor(model_data$stratum)
formula_B <- as.formula(paste("deaths ~", paste(cb_names, collapse = " + "), "+ stratum_factor"))
fit_B <- glm(formula_B, data = model_data, family = quasipoisson())
cat("    Converged:", fit_B$converged, "\n")
cat("    N coefficients:", length(coef(fit_B)), "\n")

cp_B <- tryCatch({
  crosspred(cb, fit_B, at = seq(10, 30, by = 1), cen = 22)
}, error = function(e) {
  cat("    crosspred error:", conditionMessage(e), "\n")
  NULL
})
if (!is.null(cp_B)) {
  cat("    SUCCESS! RR at 30°C:", round(cp_B$allRRfit["30"], 3), "\n")
}
if (!is.null(cp_B)) {
  cat("    SUCCESS! RR at 30°C:", round(cp_B$allRRfit["30"], 3), "\n")
}

# Model C: Extract coefficients manually from Model B
cat("\n  Model C: Manual coef extraction from Model B...\n")
all_coefs <- coef(fit_B)
cb_coef_idx <- match(cb_names, names(all_coefs))
cat("    CB coef indices:", paste(cb_coef_idx, collapse = ", "), "\n")

if (all(!is.na(cb_coef_idx))) {
  cb_coefs <- all_coefs[cb_coef_idx]
  full_vcov <- vcov(fit_B)
  cb_vcov <- full_vcov[cb_coef_idx, cb_coef_idx]
  
  cat("    CB coefs (first 4):", paste(round(cb_coefs[1:4], 4), collapse = ", "), "\n")
  
  cp_C <- tryCatch({
    crosspred(cb, coef = cb_coefs, vcov = cb_vcov, at = seq(10, 30, by = 1), cen = 22, cumul = TRUE)
  }, error = function(e) {
    cat("    crosspred error:", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(cp_C)) {
    # Explore the crosspred object structure
    cat("    SUCCESS! crosspred worked.\n")
    cat("    Object elements:", paste(names(cp_C), collapse=", "), "\n")
    cat("    predvar length:", length(cp_C$predvar), "\n")
    cat("    allRRfit length:", length(cp_C$allRRfit), "\n")
    cat("    allfit (cumulative log-RR) length:", length(cp_C$allfit), "\n")
    
    if (length(cp_C$allfit) > 0) {
      cat("    allfit values:", paste(round(head(cp_C$allfit, 6), 4), collapse=", "), "\n")
      cat("    Converting to RR...\n")
      allRR <- exp(cp_C$allfit)
      cat("    RRs:", paste(round(head(allRR, 6), 3), collapse=", "), "\n")
      
      # Find MMT
      mmt_idx <- which.min(allRR)
      mmt <- cp_C$predvar[mmt_idx]
      cat("    MMT:", mmt, "°C\n")
      
      # RR at 30°C
      rr_30_idx <- which(cp_C$predvar == 30)
      if (length(rr_30_idx) > 0) {
        cat("    RR at 30°C:", round(allRR[rr_30_idx], 3), "\n")
      }
    }
  }
} else {
  cp_C <- NULL
}

# -----------------------------------------------------------------------------
# 4. Summary
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("SUMMARY\n")
cat("=======================================================\n")
cat("Model A (no strata):     ", ifelse(!is.null(cp_A), "WORKS", "FAILS"), "\n")
cat("Model B (strata direct): ", ifelse(!is.null(cp_B), "WORKS", "FAILS"), "\n")
cat("Model C (manual coefs):  ", ifelse(!is.null(cp_C), "WORKS", "FAILS"), "\n")

cat("\nDone!\n")
