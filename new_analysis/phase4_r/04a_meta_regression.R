# =============================================================================
# 04a_meta_regression.R
# Meta-Regression Analysis: Testing Regional Moderators
# =============================================================================
#
# Tests whether temperature-mortality associations vary by:
# 1. Urbanization rate (% urban population)
# 2. GDP per capita (economic development)
# 3. Elderly proportion (demographic aging)
# 4. Climate zone (temperature variability)
#
# References:
# - Gasparrini et al. (2015) Lancet - Meta-regression for heterogeneity
# - Tobias et al. (2017) EHP - Regional moderators in temperature-mortality
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
cat("META-REGRESSION ANALYSIS:", EXPOSURE_TYPE, "\n")
cat("=======================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# Configuration
MAX_LAG <- 21
TEMP_DF <- 4
LAG_DF <- 4

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

# Mortality data
mort_file <- file.path(DATA_DIR, paste0("mortality_", EXPOSURE_TYPE, "_daily_elderly.parquet"))
mort <- as.data.table(read_parquet(mort_file))

# ERA5 temperature data
era5_file <- file.path(DATA_DIR, paste0("era5_", EXPOSURE_TYPE, "_daily.parquet"))
era5 <- as.data.table(read_parquet(era5_file))

# Regional covariates
cov_file <- file.path(DATA_DIR, paste0("ses_", EXPOSURE_TYPE, "_covariates.csv"))
covariates <- fread(cov_file)

# Standardize column names
region_col <- paste0(EXPOSURE_TYPE, "_code")
if (region_col %in% names(mort)) {
  setnames(mort, region_col, "region_code")
}
if (region_col %in% names(covariates)) {
  setnames(covariates, region_col, "region_code")
}

# Convert dates to character for consistent merge (mort has POSIXct, era5 has Date)
era5[, date := as.character(date)]
mort[, date := as.character(as.Date(date))]

# Merge mortality + temperature
setkey(mort, region_code, date)
setkey(era5, region_code, date)
data <- merge(mort, era5, by = c("region_code", "date"))

# Convert date back to Date class for model fitting (needed for format() calls)
data[, date := as.Date(date)]

# Working variables
data[, deaths := deaths_elderly]
data[, tmean := temp_mean]
data <- data[!is.na(deaths) & !is.na(tmean)]
setorder(data, region_code, date)

# Global temperature percentiles for boundaries - use GLOBAL knots for meta-analysis
temp_all <- data$tmean
temp_pcts <- quantile(temp_all, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)
temp_boundary <- c(0, 40)  # Expanded from c(0.7, 35.1) per expert recommendation

cat("  Regions:", uniqueN(data$region_code), "\n")
cat("  Date range:", as.character(min(data$date)), "to", as.character(max(data$date)), "\n")
cat("  GLOBAL knots (P10/P75/P90):", round(temp_pcts, 1), "\n")
cat("  Boundary knots:", temp_boundary, "\n")
cat("  Covariates loaded: urbanization, GDP, elderly %\n")

# -----------------------------------------------------------------------------
# 2. Regional DLNM Fitting
# -----------------------------------------------------------------------------
cat("\n[2] Fitting region-specific DLNMs...\n")

regions <- unique(data$region_code)
n_regions <- length(regions)

# Storage
mvmeta_coefs <- list()
mvmeta_vcovs <- list()
mvmeta_regions <- c()
regional_rrs <- list()

for (i in seq_along(regions)) {
  reg <- regions[i]
  
  if (i %% 50 == 0) {
    cat(sprintf("  Progress: %d/%d regions\n", i, n_regions))
  }
  
  daily <- data[region_code == reg][order(date)]
  
  if (nrow(daily) < 365) next
  
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
        log_rr_p99 = log(rr_p99),
        log_rr_p1 = log(rr_p1),
        se_p99 = rr_p99_se,
        se_p1 = rr_p1_se,
        mmt = mmt
      )
    }
    
    "success"
    
  }, error = function(e) {
    "error"
  })
}

n_success <- length(mvmeta_regions)
cat(sprintf("  Successful fits: %d/%d (%.1f%%)\n", 
            n_success, n_regions, 100 * n_success / n_regions))

if (n_success < 10) {
  stop("Too few regions for meta-regression")
}

# -----------------------------------------------------------------------------
# 3. Prepare Meta-Regression Data
# -----------------------------------------------------------------------------
cat("\n[3] Preparing meta-regression data...\n")

# Extract regional estimates
rr_data <- data.table(
  region_code = names(regional_rrs),
  log_rr_heat = sapply(regional_rrs, function(x) x$log_rr_p99),
  se_heat = sapply(regional_rrs, function(x) x$se_p99),
  log_rr_cold = sapply(regional_rrs, function(x) x$log_rr_p1),
  se_cold = sapply(regional_rrs, function(x) x$se_p1)
)

# Convert region_code to match covariates type
rr_data[, region_code := as.integer(region_code)]

# Merge with covariates
rr_data <- merge(rr_data, covariates[, .(region_code, urbanization_rate, gdp_per_capita, 
                                          pct_elderly, pop_elderly)], 
                 by = "region_code")

# Calculate temperature variability per region
temp_var <- data[, .(temp_sd = sd(tmean, na.rm = TRUE),
                     temp_range = diff(range(tmean, na.rm = TRUE))),
                 by = region_code]
rr_data <- merge(rr_data, temp_var, by = "region_code")

# Standardize covariates (z-scores)
rr_data[, urban_z := (urbanization_rate - mean(urbanization_rate)) / sd(urbanization_rate)]
rr_data[, gdp_z := (gdp_per_capita - mean(gdp_per_capita)) / sd(gdp_per_capita)]
rr_data[, elderly_z := (pct_elderly - mean(pct_elderly)) / sd(pct_elderly)]
rr_data[, temp_var_z := (temp_sd - mean(temp_sd)) / sd(temp_sd)]

cat(sprintf("  Meta-regression data: %d regions\n", nrow(rr_data)))
cat(sprintf("  Covariates: urban (%.1f-%.1f%%), GDP ($%.0f-$%.0f), elderly (%.1f-%.1f%%)\n",
            min(rr_data$urbanization_rate), max(rr_data$urbanization_rate),
            min(rr_data$gdp_per_capita), max(rr_data$gdp_per_capita),
            min(rr_data$pct_elderly), max(rr_data$pct_elderly)))

# -----------------------------------------------------------------------------
# 4. Meta-Regression Models
# -----------------------------------------------------------------------------
cat("\n[4] Running meta-regression models...\n")

meta_reg_results <- list()

# Helper function to safely extract meta-regression results
safe_metareg_results <- function(fit, name) {
  tryCatch({
    summ <- summary(fit)
    list(
      intercept = coef(fit)[1],
      slope = coef(fit)[2],
      intercept_se = sqrt(vcov(fit)[1,1]),
      slope_se = sqrt(vcov(fit)[2,2]),
      p_value = summ$coefficients[2,4]
    )
  }, error = function(e) {
    cat(sprintf("    Warning: Could not extract %s summary: %s\n", name, e$message))
    # Fallback: use Wald test
    slope <- coef(fit)[2]
    se <- tryCatch(sqrt(vcov(fit)[2,2]), error = function(e) NA)
    z <- if (!is.na(se) && se > 0) slope / se else NA
    p <- if (!is.na(z)) 2 * pnorm(-abs(z)) else NA
    list(
      intercept = coef(fit)[1],
      slope = slope,
      intercept_se = tryCatch(sqrt(vcov(fit)[1,1]), error = function(e) NA),
      slope_se = se,
      p_value = p
    )
  })
}

# Model 1: Urbanization
cat("\n  Model 1: Urbanization moderator\n")
meta_urban_heat <- tryCatch(
  mixmeta(log_rr_heat ~ urban_z, S = se_heat^2, data = rr_data, method = "reml"),
  error = function(e) {
    cat("    REML failed, trying ML...\n")
    mixmeta(log_rr_heat ~ urban_z, S = se_heat^2, data = rr_data, method = "ml")
  }
)
meta_urban_cold <- tryCatch(
  mixmeta(log_rr_cold ~ urban_z, S = se_cold^2, data = rr_data, method = "reml"),
  error = function(e) {
    cat("    REML failed, trying ML...\n")
    mixmeta(log_rr_cold ~ urban_z, S = se_cold^2, data = rr_data, method = "ml")
  }
)

meta_reg_results$urbanization <- list(
  heat = safe_metareg_results(meta_urban_heat, "urban_heat"),
  cold = safe_metareg_results(meta_urban_cold, "urban_cold")
)

cat(sprintf("    Heat: slope=%.4f, SE=%.4f, p=%.4f\n",
            meta_reg_results$urbanization$heat$slope,
            meta_reg_results$urbanization$heat$slope_se,
            meta_reg_results$urbanization$heat$p_value))
cat(sprintf("    Cold: slope=%.4f, SE=%.4f, p=%.4f\n",
            meta_reg_results$urbanization$cold$slope,
            meta_reg_results$urbanization$cold$slope_se,
            meta_reg_results$urbanization$cold$p_value))

# Model 2: GDP per capita
cat("\n  Model 2: GDP moderator\n")
meta_gdp_heat <- tryCatch(
  mixmeta(log_rr_heat ~ gdp_z, S = se_heat^2, data = rr_data, method = "reml"),
  error = function(e) mixmeta(log_rr_heat ~ gdp_z, S = se_heat^2, data = rr_data, method = "ml")
)
meta_gdp_cold <- tryCatch(
  mixmeta(log_rr_cold ~ gdp_z, S = se_cold^2, data = rr_data, method = "reml"),
  error = function(e) mixmeta(log_rr_cold ~ gdp_z, S = se_cold^2, data = rr_data, method = "ml")
)

meta_reg_results$gdp <- list(
  heat = safe_metareg_results(meta_gdp_heat, "gdp_heat"),
  cold = safe_metareg_results(meta_gdp_cold, "gdp_cold")
)

cat(sprintf("    Heat: slope=%.4f, SE=%.4f, p=%.4f\n",
            meta_reg_results$gdp$heat$slope,
            meta_reg_results$gdp$heat$slope_se,
            meta_reg_results$gdp$heat$p_value))
cat(sprintf("    Cold: slope=%.4f, SE=%.4f, p=%.4f\n",
            meta_reg_results$gdp$cold$slope,
            meta_reg_results$gdp$cold$slope_se,
            meta_reg_results$gdp$cold$p_value))

# Model 3: Elderly proportion
cat("\n  Model 3: Elderly % moderator\n")
meta_elderly_heat <- tryCatch(
  mixmeta(log_rr_heat ~ elderly_z, S = se_heat^2, data = rr_data, method = "reml"),
  error = function(e) mixmeta(log_rr_heat ~ elderly_z, S = se_heat^2, data = rr_data, method = "ml")
)
meta_elderly_cold <- tryCatch(
  mixmeta(log_rr_cold ~ elderly_z, S = se_cold^2, data = rr_data, method = "reml"),
  error = function(e) mixmeta(log_rr_cold ~ elderly_z, S = se_cold^2, data = rr_data, method = "ml")
)

meta_reg_results$elderly <- list(
  heat = safe_metareg_results(meta_elderly_heat, "elderly_heat"),
  cold = safe_metareg_results(meta_elderly_cold, "elderly_cold")
)

cat(sprintf("    Heat: slope=%.4f, SE=%.4f, p=%.4f\n",
            meta_reg_results$elderly$heat$slope,
            meta_reg_results$elderly$heat$slope_se,
            meta_reg_results$elderly$heat$p_value))
cat(sprintf("    Cold: slope=%.4f, SE=%.4f, p=%.4f\n",
            meta_reg_results$elderly$cold$slope,
            meta_reg_results$elderly$cold$slope_se,
            meta_reg_results$elderly$cold$p_value))

# Model 4: Temperature variability
cat("\n  Model 4: Climate variability moderator\n")
meta_tempvar_heat <- tryCatch(
  mixmeta(log_rr_heat ~ temp_var_z, S = se_heat^2, data = rr_data, method = "reml"),
  error = function(e) mixmeta(log_rr_heat ~ temp_var_z, S = se_heat^2, data = rr_data, method = "ml")
)
meta_tempvar_cold <- tryCatch(
  mixmeta(log_rr_cold ~ temp_var_z, S = se_cold^2, data = rr_data, method = "reml"),
  error = function(e) mixmeta(log_rr_cold ~ temp_var_z, S = se_cold^2, data = rr_data, method = "ml")
)

meta_reg_results$temp_variability <- list(
  heat = safe_metareg_results(meta_tempvar_heat, "tempvar_heat"),
  cold = safe_metareg_results(meta_tempvar_cold, "tempvar_cold")
)

cat(sprintf("    Heat: slope=%.4f, SE=%.4f, p=%.4f\n",
            meta_reg_results$temp_variability$heat$slope,
            meta_reg_results$temp_variability$heat$slope_se,
            meta_reg_results$temp_variability$heat$p_value))
cat(sprintf("    Cold: slope=%.4f, SE=%.4f, p=%.4f\n",
            meta_reg_results$temp_variability$cold$slope,
            meta_reg_results$temp_variability$cold$slope_se,
            meta_reg_results$temp_variability$cold$p_value))

# -----------------------------------------------------------------------------
# 5. Summary
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("META-REGRESSION SUMMARY\n")
cat("=======================================================\n")

interpret_slope <- function(slope, p_value) {
  if (is.na(p_value)) return("NA")
  if (p_value < 0.001) return("***")
  if (p_value < 0.01) return("**")
  if (p_value < 0.05) return("*")
  return("ns")
}

cat("\nModerat or           Heat Effect            Cold Effect\n")
cat("-----------------------------------------------------------\n")
cat(sprintf("Urban (%%)       %+.4f (%s)          %+.4f (%s)\n",
            meta_reg_results$urbanization$heat$slope,
            interpret_slope(meta_reg_results$urbanization$heat$slope, 
                          meta_reg_results$urbanization$heat$p_value),
            meta_reg_results$urbanization$cold$slope,
            interpret_slope(meta_reg_results$urbanization$cold$slope,
                          meta_reg_results$urbanization$cold$p_value)))
cat(sprintf("GDP              %+.4f (%s)          %+.4f (%s)\n",
            meta_reg_results$gdp$heat$slope,
            interpret_slope(meta_reg_results$gdp$heat$slope,
                          meta_reg_results$gdp$heat$p_value),
            meta_reg_results$gdp$cold$slope,
            interpret_slope(meta_reg_results$gdp$cold$slope,
                          meta_reg_results$gdp$cold$p_value)))
cat(sprintf("Elderly (%%)      %+.4f (%s)          %+.4f (%s)\n",
            meta_reg_results$elderly$heat$slope,
            interpret_slope(meta_reg_results$elderly$heat$slope,
                          meta_reg_results$elderly$heat$p_value),
            meta_reg_results$elderly$cold$slope,
            interpret_slope(meta_reg_results$elderly$cold$slope,
                          meta_reg_results$elderly$cold$p_value)))
cat(sprintf("Temp Var (SD)    %+.4f (%s)          %+.4f (%s)\n",
            meta_reg_results$temp_variability$heat$slope,
            interpret_slope(meta_reg_results$temp_variability$heat$slope,
                          meta_reg_results$temp_variability$heat$p_value),
            meta_reg_results$temp_variability$cold$slope,
            interpret_slope(meta_reg_results$temp_variability$cold$slope,
                          meta_reg_results$temp_variability$cold$p_value)))
cat("\nSignificance: *** p<0.001, ** p<0.01, * p<0.05, ns=not significant\n")
cat("Note: Slopes are for standardized covariates (1 SD change)\n")

# -----------------------------------------------------------------------------
# 6. Save Results
# -----------------------------------------------------------------------------
cat("\n[6] Saving results...\n")

output <- list(
  exposure_type = EXPOSURE_TYPE,
  analysis_date = as.character(Sys.time()),
  n_regions = n_success,
  
  moderators = meta_reg_results,
  
  regional_data = rr_data
)

output_file <- file.path(OUTPUT_DIR, paste0("meta_regression_", EXPOSURE_TYPE, ".json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("  JSON saved to:", output_file, "\n")

cat("\n=======================================================\n")
cat("Done!", as.character(Sys.time()), "\n")
cat("=======================================================\n")
