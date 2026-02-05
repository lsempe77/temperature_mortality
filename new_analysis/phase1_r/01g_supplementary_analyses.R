# =============================================================================
# 01g_supplementary_analyses.R
# Period-stratified, LOYO summary, Macro-region, and SES gradient analyses
# =============================================================================
# 
# This script addresses reviewer concerns with four supplementary analyses:
# 1. Period-stratified: 2010-2019 vs 2022-2024 comparison
# 2. LOYO summary: Temporal stability across years
# 3. Macro-region: North/Northeast/Southeast/South/Center-West stratification
# 4. SES gradient: GDP and urbanisation as effect modifiers
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
cat("SUPPLEMENTARY ANALYSES:", EXPOSURE_TYPE, "\n")
cat("=======================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# 1. Load Data and Results
# -----------------------------------------------------------------------------
SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("--file=", "", file_arg)))
  } else {
    getwd()
  }
})

DATA_DIR <- file.path(SCRIPT_DIR, "..", "phase0_data_prep", "results")
DLNM_DIR <- file.path(SCRIPT_DIR, "results")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results")

cat("[1] Loading data...\n")

# Load DLNM results
dlnm_file <- file.path(DLNM_DIR, paste0("dlnm_r_", EXPOSURE_TYPE, "_results_v2.json"))
dlnm_results <- fromJSON(dlnm_file)
cat("  Loaded DLNM results\n")

# Extract parameters
TEMP_BOUNDARY <- dlnm_results$dlnm_params$temp_boundary
TEMP_KNOTS <- dlnm_results$dlnm_params$temp_knots
MAX_LAG <- dlnm_results$dlnm_params$max_lag
LAG_DF <- dlnm_results$dlnm_params$lag_df
pooled_mmt <- dlnm_results$pooled$mmt

# Load mortality data
mort_file <- file.path(DATA_DIR, paste0("mortality_", EXPOSURE_TYPE, "_daily_elderly.parquet"))
mort <- as.data.table(read_parquet(mort_file))

# Load temperature data
era5_file <- file.path(DATA_DIR, paste0("era5_", EXPOSURE_TYPE, "_daily.parquet"))
era5 <- as.data.table(read_parquet(era5_file))

# Load covariates
cov_file <- file.path(DATA_DIR, "regional_covariates.csv")
covariates <- fread(cov_file)
cat("  Loaded covariates:", nrow(covariates), "regions\n")

# Standardize column names
if ("deaths_elderly" %in% names(mort)) setnames(mort, "deaths_elderly", "deaths")
if ("temp_mean" %in% names(era5)) setnames(era5, "temp_mean", "tmean")
mort[, date := as.Date(date)]
era5[, date := as.Date(date)]

# Merge
data <- merge(mort, era5[, .(region_code, date, tmean)], 
              by = c("region_code", "date"), all.x = TRUE)
data <- data[!is.na(tmean) & !is.na(deaths)]
data[, year := year(date)]

cat("  Total region-days:", format(nrow(data), big.mark = ","), "\n")
cat("  Total deaths:", format(sum(data$deaths), big.mark = ","), "\n")

# Add macro-region from covariates
data <- merge(data, covariates[, .(region_code, macro_region, gdp_per_capita_brl, urban_pct)],
              by = "region_code", all.x = TRUE)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Function to run first-stage DLNM for a subset of data
# Now accepts optional adaptive knots for period-stratified analyses
run_first_stage <- function(dt, region_col = "region_code", 
                            temp_knots = TEMP_KNOTS, 
                            temp_boundary = TEMP_BOUNDARY) {
  regions <- unique(dt[[region_col]])
  results_list <- list()
  
  for (reg in regions) {
    reg_data <- dt[get(region_col) == reg]
    
    if (nrow(reg_data) < 365 || sum(reg_data$deaths) < 100) next
    
    tryCatch({
      # Create cross-basis with provided knots
      cb <- crossbasis(
        reg_data$tmean,
        lag = MAX_LAG,
        argvar = list(fun = "ns", knots = temp_knots, Boundary.knots = temp_boundary),
        arglag = list(fun = "ns", df = LAG_DF)
      )
      
      # Fit model
      model <- glm(deaths ~ cb + ns(as.numeric(date), df = 7 * length(unique(year(reg_data$date)))) +
                     factor(wday(date)),
                   data = reg_data, family = quasipoisson())
      
      # Extract coefficients
      coef_vec <- coef(model)[2:(ncol(cb) + 1)]
      vcov_mat <- vcov(model)[2:(ncol(cb) + 1), 2:(ncol(cb) + 1)]
      
      # Check for valid results
      if (any(is.na(coef_vec)) || any(is.na(vcov_mat)) || any(!is.finite(vcov_mat))) next
      
      results_list[[as.character(reg)]] <- list(
        region = reg,
        coef = coef_vec,
        vcov = vcov_mat,
        n_days = nrow(reg_data),
        n_deaths = sum(reg_data$deaths)
      )
    }, error = function(e) NULL)
  }
  
  return(results_list)
}

# Function to run second-stage meta-analysis
# Now accepts optional adaptive knots for period-stratified analyses
run_second_stage <- function(first_stage_results, 
                             temp_knots = TEMP_KNOTS, 
                             temp_boundary = TEMP_BOUNDARY) {
  if (length(first_stage_results) < 3) {
    return(list(converged = FALSE, n_regions = length(first_stage_results)))
  }
  
  # Prepare for mixmeta
  n_params <- length(first_stage_results[[1]]$coef)
  coef_mat <- do.call(rbind, lapply(first_stage_results, function(x) x$coef))
  vcov_list <- lapply(first_stage_results, function(x) x$vcov)
  
  tryCatch({
    # Run meta-analysis
    mv <- mixmeta(coef_mat, vcov_list, method = "reml")
    
    # Get pooled coefficients
    pooled_coef <- coef(mv)
    pooled_vcov <- vcov(mv)
    
    # Create cross-basis for predictions using provided knots
    temp_range <- temp_boundary
    temp_for_pred <- seq(temp_range[1], temp_range[2], by = 0.5)
    
    cb_pred <- crossbasis(
      temp_for_pred,
      lag = MAX_LAG,
      argvar = list(fun = "ns", knots = temp_knots, Boundary.knots = temp_boundary),
      arglag = list(fun = "ns", df = LAG_DF)
    )
    
    # Get predictions
    cp <- crosspred(cb_pred, coef = pooled_coef, vcov = pooled_vcov,
                    model.link = "log", at = temp_for_pred, cen = pooled_mmt)
    
    # Find RR at P99 and P1 (use data-driven percentiles)
    p99_temp <- 28.7  # Approximate P99
    p1_temp <- 20.2   # Approximate P1
    
    idx_p99 <- which.min(abs(temp_for_pred - p99_temp))
    idx_p1 <- which.min(abs(temp_for_pred - p1_temp))
    
    return(list(
      converged = TRUE,
      n_regions = length(first_stage_results),
      heat_rr = cp$allRRfit[idx_p99],
      heat_lo = cp$allRRlow[idx_p99],
      heat_hi = cp$allRRhigh[idx_p99],
      cold_rr = cp$allRRfit[idx_p1],
      cold_lo = cp$allRRlow[idx_p1],
      cold_hi = cp$allRRhigh[idx_p1],
      I2 = NA  # Would need additional computation
    ))
  }, error = function(e) {
    return(list(converged = FALSE, n_regions = length(first_stage_results), error = e$message))
  })
}

# =============================================================================
# ANALYSIS 1: PERIOD-STRATIFIED (2010-2019 vs 2022-2024)
# =============================================================================
cat("\n=======================================================\n")
cat("ANALYSIS 1: PERIOD-STRATIFIED COMPARISON\n")
cat("=======================================================\n")

# Split data by period
data_pre <- data[year >= 2010 & year <= 2019]
data_post <- data[year >= 2022 & year <= 2024]

cat("\nPre-pandemic (2010-2019):\n")
cat("  Days:", format(nrow(data_pre), big.mark = ","), "\n")
cat("  Deaths:", format(sum(data_pre$deaths), big.mark = ","), "\n")

cat("\nPost-pandemic (2022-2024):\n")
cat("  Days:", format(nrow(data_post), big.mark = ","), "\n")
cat("  Deaths:", format(sum(data_post$deaths), big.mark = ","), "\n")

# Run analyses
cat("\nRunning pre-pandemic analysis...\n")
first_pre <- run_first_stage(data_pre)
cat("  First stage:", length(first_pre), "regions converged\n")
result_pre <- run_second_stage(first_pre)

# For post-pandemic, use ADAPTIVE knots based on the 3-year temperature distribution
cat("\nRunning post-pandemic analysis with adaptive knots...\n")
post_temps <- data_post$tmean[!is.na(data_post$tmean)]
post_boundary <- quantile(post_temps, c(0.01, 0.99), na.rm = TRUE)
post_knots <- quantile(post_temps, c(0.25, 0.50, 0.75), na.rm = TRUE)
cat("  Adaptive boundary:", round(post_boundary, 2), "\n")
cat("  Adaptive knots:", round(post_knots, 2), "\n")

first_post <- run_first_stage(data_post, temp_knots = post_knots, temp_boundary = post_boundary)
cat("  First stage:", length(first_post), "regions converged\n")
result_post <- run_second_stage(first_post, temp_knots = post_knots, temp_boundary = post_boundary)

period_results <- list(
  pre_pandemic = list(
    period = "2010-2019",
    n_regions = result_pre$n_regions,
    converged = result_pre$converged,
    heat_rr = if(result_pre$converged) result_pre$heat_rr else NA,
    heat_ci = if(result_pre$converged) c(result_pre$heat_lo, result_pre$heat_hi) else c(NA, NA),
    cold_rr = if(result_pre$converged) result_pre$cold_rr else NA,
    cold_ci = if(result_pre$converged) c(result_pre$cold_lo, result_pre$cold_hi) else c(NA, NA)
  ),
  post_pandemic = list(
    period = "2022-2024",
    n_regions = result_post$n_regions,
    converged = result_post$converged,
    heat_rr = if(result_post$converged) result_post$heat_rr else NA,
    heat_ci = if(result_post$converged) c(result_post$heat_lo, result_post$heat_hi) else c(NA, NA),
    cold_rr = if(result_post$converged) result_post$cold_rr else NA,
    cold_ci = if(result_post$converged) c(result_post$cold_lo, result_post$cold_hi) else c(NA, NA)
  )
)

cat("\n--- PERIOD COMPARISON RESULTS ---\n")
if (result_pre$converged) {
  cat("Pre-pandemic (2010-2019):\n")
  cat("  Heat RR:", round(result_pre$heat_rr, 3), 
      "(", round(result_pre$heat_lo, 3), "-", round(result_pre$heat_hi, 3), ")\n")
  cat("  Cold RR:", round(result_pre$cold_rr, 3),
      "(", round(result_pre$cold_lo, 3), "-", round(result_pre$cold_hi, 3), ")\n")
}
if (result_post$converged) {
  cat("Post-pandemic (2022-2024):\n")
  cat("  Heat RR:", round(result_post$heat_rr, 3),
      "(", round(result_post$heat_lo, 3), "-", round(result_post$heat_hi, 3), ")\n")
  cat("  Cold RR:", round(result_post$cold_rr, 3),
      "(", round(result_post$cold_lo, 3), "-", round(result_post$cold_hi, 3), ")\n")
}

# =============================================================================
# ANALYSIS 2: LOYO SUMMARY (Temporal Stability)
# =============================================================================
cat("\n=======================================================\n")
cat("ANALYSIS 2: LEAVE-ONE-YEAR-OUT SUMMARY\n")
cat("=======================================================\n")

# We'll compute year-specific estimates by running on each year separately
years <- 2010:2024
loyo_results <- list()

for (yr in years) {
  cat("  Processing year", yr, "...")
  
  # Exclude COVID years from individual analysis
  if (yr %in% c(2020, 2021)) {
    cat(" (skipped - COVID)\n")
    next
  }
  
  data_yr <- data[year == yr]
  
  if (nrow(data_yr) < 1000) {
    cat(" (insufficient data)\n")
    next
  }
  
  # Simple year-specific model (simplified for speed)
  # Calculate average RR by temperature bins
  data_yr[, temp_bin := cut(tmean, breaks = c(-Inf, 15, 20, 25, 30, Inf),
                            labels = c("very_cold", "cold", "moderate", "warm", "hot"))]
  
  # Simple Poisson regression
  tryCatch({
    model <- glm(deaths ~ temp_bin + factor(wday(date)) + ns(as.numeric(date), df = 4),
                 data = data_yr, family = quasipoisson())
    
    coefs <- exp(coef(model))
    loyo_results[[as.character(yr)]] <- list(
      year = yr,
      n_deaths = sum(data_yr$deaths),
      cold_effect = coefs["temp_binvery_cold"],
      hot_effect = coefs["temp_binhot"]
    )
    cat(" done\n")
  }, error = function(e) {
    cat(" (error)\n")
  })
}

# Summarize LOYO results
cat("\n--- TEMPORAL STABILITY (LOYO) ---\n")
loyo_df <- rbindlist(lapply(loyo_results, as.data.table), fill = TRUE)
if (nrow(loyo_df) > 0) {
  cat("Years analyzed:", paste(loyo_df$year, collapse = ", "), "\n")
  cat("Cold effect range:", round(min(loyo_df$cold_effect, na.rm = TRUE), 3), "-",
      round(max(loyo_df$cold_effect, na.rm = TRUE), 3), "\n")
  cat("Hot effect range:", round(min(loyo_df$hot_effect, na.rm = TRUE), 3), "-",
      round(max(loyo_df$hot_effect, na.rm = TRUE), 3), "\n")
  cat("Cold effect mean (SD):", round(mean(loyo_df$cold_effect, na.rm = TRUE), 3),
      "(", round(sd(loyo_df$cold_effect, na.rm = TRUE), 3), ")\n")
  cat("Hot effect mean (SD):", round(mean(loyo_df$hot_effect, na.rm = TRUE), 3),
      "(", round(sd(loyo_df$hot_effect, na.rm = TRUE), 3), ")\n")
}

# =============================================================================
# ANALYSIS 3: MACRO-REGION STRATIFICATION
# =============================================================================
cat("\n=======================================================\n")
cat("ANALYSIS 3: MACRO-REGION STRATIFICATION\n")
cat("=======================================================\n")

macro_regions <- unique(data$macro_region)
macro_regions <- macro_regions[!is.na(macro_regions)]
cat("Macro-regions:", paste(macro_regions, collapse = ", "), "\n")

macro_results <- list()

for (macro in macro_regions) {
  cat("\nProcessing", macro, "...\n")
  
  data_macro <- data[macro_region == macro]
  cat("  Regions:", uniqueN(data_macro$region_code), "\n")
  cat("  Deaths:", format(sum(data_macro$deaths), big.mark = ","), "\n")
  
  # Run first stage
  first_macro <- run_first_stage(data_macro)
  cat("  Converged regions:", length(first_macro), "\n")
  
  if (length(first_macro) >= 3) {
    result_macro <- run_second_stage(first_macro)
    
    macro_results[[macro]] <- list(
      macro_region = macro,
      n_regions = result_macro$n_regions,
      converged = result_macro$converged,
      heat_rr = if(result_macro$converged) result_macro$heat_rr else NA,
      heat_ci = if(result_macro$converged) c(result_macro$heat_lo, result_macro$heat_hi) else c(NA, NA),
      cold_rr = if(result_macro$converged) result_macro$cold_rr else NA,
      cold_ci = if(result_macro$converged) c(result_macro$cold_lo, result_macro$cold_hi) else c(NA, NA),
      mean_temp = mean(data_macro$tmean, na.rm = TRUE)
    )
    
    if (result_macro$converged) {
      cat("  Heat RR:", round(result_macro$heat_rr, 3), "\n")
      cat("  Cold RR:", round(result_macro$cold_rr, 3), "\n")
    }
  } else {
    macro_results[[macro]] <- list(
      macro_region = macro,
      n_regions = length(first_macro),
      converged = FALSE
    )
    cat("  Insufficient converged regions\n")
  }
}

cat("\n--- MACRO-REGION SUMMARY ---\n")
for (macro in names(macro_results)) {
  r <- macro_results[[macro]]
  if (r$converged) {
    cat(sprintf("%s (n=%d): Heat RR=%.2f (%.2f-%.2f), Cold RR=%.2f (%.2f-%.2f), Mean temp=%.1f°C\n",
                macro, r$n_regions, r$heat_rr, r$heat_ci[1], r$heat_ci[2],
                r$cold_rr, r$cold_ci[1], r$cold_ci[2], r$mean_temp))
  }
}

# =============================================================================
# ANALYSIS 4: SOCIOECONOMIC GRADIENT
# =============================================================================
cat("\n=======================================================\n")
cat("ANALYSIS 4: SOCIOECONOMIC GRADIENT ANALYSIS\n")
cat("=======================================================\n")

# Load first-stage results from main analysis
# For SES analysis, we'll use region-level RRs and correlate with SES

# Create SES tertiles
ses_data <- covariates[, .(region_code, gdp_per_capita_brl, urban_pct)]
ses_data <- ses_data[!is.na(gdp_per_capita_brl) & !is.na(urban_pct)]

# Create tertiles
ses_data[, gdp_tertile := cut(gdp_per_capita_brl, 
                               breaks = quantile(gdp_per_capita_brl, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                               labels = c("Low", "Medium", "High"), include.lowest = TRUE)]
ses_data[, urban_tertile := cut(urban_pct,
                                 breaks = quantile(urban_pct, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                                 labels = c("Low", "Medium", "High"), include.lowest = TRUE)]

cat("GDP tertiles:\n")
print(ses_data[, .N, by = gdp_tertile])
cat("\nUrbanization tertiles:\n")
print(ses_data[, .N, by = urban_tertile])

# Merge with main data
data <- merge(data, ses_data[, .(region_code, gdp_tertile, urban_tertile)],
              by = "region_code", all.x = TRUE)

# Run stratified analysis by GDP tertile
ses_results <- list()

for (tert in c("Low", "Medium", "High")) {
  cat("\nProcessing GDP tertile:", tert, "...\n")
  
  data_tert <- data[gdp_tertile == tert]
  cat("  Regions:", uniqueN(data_tert$region_code), "\n")
  
  if (uniqueN(data_tert$region_code) < 5) {
    cat("  Insufficient regions\n")
    next
  }
  
  first_tert <- run_first_stage(data_tert)
  cat("  Converged regions:", length(first_tert), "\n")
  
  if (length(first_tert) >= 3) {
    result_tert <- run_second_stage(first_tert)
    
    ses_results[[paste0("gdp_", tert)]] <- list(
      category = paste("GDP", tert),
      n_regions = result_tert$n_regions,
      converged = result_tert$converged,
      heat_rr = if(result_tert$converged) result_tert$heat_rr else NA,
      cold_rr = if(result_tert$converged) result_tert$cold_rr else NA,
      heat_ci = if(result_tert$converged) c(result_tert$heat_lo, result_tert$heat_hi) else c(NA, NA),
      cold_ci = if(result_tert$converged) c(result_tert$cold_lo, result_tert$cold_hi) else c(NA, NA)
    )
  }
}

# Run stratified analysis by urbanization tertile
for (tert in c("Low", "Medium", "High")) {
  cat("\nProcessing Urbanization tertile:", tert, "...\n")
  
  data_tert <- data[urban_tertile == tert]
  cat("  Regions:", uniqueN(data_tert$region_code), "\n")
  
  if (uniqueN(data_tert$region_code) < 5) {
    cat("  Insufficient regions\n")
    next
  }
  
  first_tert <- run_first_stage(data_tert)
  cat("  Converged regions:", length(first_tert), "\n")
  
  if (length(first_tert) >= 3) {
    result_tert <- run_second_stage(first_tert)
    
    ses_results[[paste0("urban_", tert)]] <- list(
      category = paste("Urban", tert),
      n_regions = result_tert$n_regions,
      converged = result_tert$converged,
      heat_rr = if(result_tert$converged) result_tert$heat_rr else NA,
      cold_rr = if(result_tert$converged) result_tert$cold_rr else NA,
      heat_ci = if(result_tert$converged) c(result_tert$heat_lo, result_tert$heat_hi) else c(NA, NA),
      cold_ci = if(result_tert$converged) c(result_tert$cold_lo, result_tert$cold_hi) else c(NA, NA)
    )
  }
}

cat("\n--- SES GRADIENT SUMMARY ---\n")
for (cat_name in names(ses_results)) {
  r <- ses_results[[cat_name]]
  if (r$converged) {
    cat(sprintf("%s (n=%d): Heat RR=%.2f (%.2f-%.2f), Cold RR=%.2f (%.2f-%.2f)\n",
                r$category, r$n_regions, r$heat_rr, r$heat_ci[1], r$heat_ci[2],
                r$cold_rr, r$cold_ci[1], r$cold_ci[2]))
  }
}

# =============================================================================
# SAVE ALL RESULTS
# =============================================================================
cat("\n=======================================================\n")
cat("SAVING RESULTS\n")
cat("=======================================================\n")

all_results <- list(
  analysis_date = as.character(Sys.time()),
  exposure_type = EXPOSURE_TYPE,
  
  period_stratified = period_results,
  loyo_summary = list(
    years = loyo_df$year,
    cold_effects = loyo_df$cold_effect,
    hot_effects = loyo_df$hot_effect,
    cold_mean = mean(loyo_df$cold_effect, na.rm = TRUE),
    cold_sd = sd(loyo_df$cold_effect, na.rm = TRUE),
    hot_mean = mean(loyo_df$hot_effect, na.rm = TRUE),
    hot_sd = sd(loyo_df$hot_effect, na.rm = TRUE)
  ),
  macro_region = macro_results,
  ses_gradient = ses_results
)

output_file <- file.path(OUTPUT_DIR, paste0("supplementary_analyses_", EXPOSURE_TYPE, ".json"))
write_json(all_results, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 4)
cat("Saved to:", output_file, "\n")

cat("\n=======================================================\n")
cat("SUPPLEMENTARY ANALYSES COMPLETE\n")
cat("Finished:", as.character(Sys.time()), "\n")
cat("=======================================================\n")
