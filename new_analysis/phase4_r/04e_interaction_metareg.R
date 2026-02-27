# =============================================================================
# 04e_interaction_metareg.R
# Meta-Regression: Region × Sex × Deprivation Interactions
# =============================================================================
#
# Purpose:
#   Tests whether temperature-mortality associations vary by the interaction
#   of geographic region (macro-region), sex, and socioeconomic deprivation (HDI).
#
# Approach:
#   1. Fit per-region sex-specific DLNMs (re-uses 04c logic, but saves per-region estimates)
#   2. Stack into ~258 rows (129 regions × 2 sexes)
#   3. Merge HDI as deprivation proxy (inverted: lower HDI = higher deprivation)
#   4. Run mixmeta with sex × deprivation × macro_region interactions
#
# Unit of analysis: region-sex (intermediate level, ~129 regions × 2 sexes)
#
# Moderators:
#   - sex (binary: female=1)
#   - deprivation_z (standardized inverse-HDI: higher = more deprived)
#   - macro_region (North, Northeast, Southeast, South, Central-West)
#   - Two-way and three-way interactions
#
# References:
#   - Gasparrini et al. (2015) Lancet - Meta-regression for heterogeneity
#   - Benmarhnia et al. (2015) EHP - Sex differences in heat vulnerability
#   - Sera et al. (2019) IJE - Urban characteristics and vulnerability
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

# Get exposure type from command line args (only intermediate supported here)
args <- commandArgs(trailingOnly = TRUE)
EXPOSURE_TYPE <- if (length(args) > 0) args[1] else "intermediate"

cat("=======================================================\n")
cat("INTERACTION META-REGRESSION: Region × Sex × HDI\n")
cat("Exposure level:", EXPOSURE_TYPE, "\n")
cat("=======================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# Configuration (must match 04c to ensure comparability)
MAX_LAG <- 21
TEMP_DF <- 4
LAG_DF  <- 4

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------
SCRIPT_DIR <- tryCatch({
  dirname(normalizePath(sys.frame(1)$ofile))
}, error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("--file=", "", file_arg)))
  } else {
    getwd()
  }
})

DATA_DIR   <- file.path(SCRIPT_DIR, "..", "phase0_data_prep", "results")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

cat("[1] Loading data...\n")

# ERA5 temperature (shared across sexes)
era5_file <- file.path(DATA_DIR, paste0("era5_", EXPOSURE_TYPE, "_daily.parquet"))
era5 <- as.data.table(read_parquet(era5_file))
era5[, date := as.character(date)]

# Regional covariates (HDI, macro_region)
reg_cov_file <- file.path(DATA_DIR, "regional_covariates.csv")
reg_covariates <- fread(reg_cov_file)
cat("  Regional covariates columns:", paste(names(reg_covariates), collapse=", "), "\n")

# SES covariates (GDP per capita at intermediate level)
ses_cov_file <- file.path(DATA_DIR, paste0("ses_", EXPOSURE_TYPE, "_covariates.csv"))
ses_covariates <- fread(ses_cov_file)

# Standardize region_code column name in SES covariates
region_col_ses <- paste0(EXPOSURE_TYPE, "_code")
if (region_col_ses %in% names(ses_covariates)) {
  setnames(ses_covariates, region_col_ses, "region_code")
}

# Sex-specific mortality files
sex_mort_files <- list(
  male   = file.path(DATA_DIR, paste0("mortality_", EXPOSURE_TYPE, "_daily_male.parquet")),
  female = file.path(DATA_DIR, paste0("mortality_", EXPOSURE_TYPE, "_daily_female.parquet"))
)

missing <- !sapply(sex_mort_files, file.exists)
if (any(missing)) {
  cat("\n*** ERROR: Missing sex-stratified mortality files ***\n")
  for (nm in names(sex_mort_files)[missing]) {
    cat("  -", sex_mort_files[[nm]], "\n")
  }
  stop("Run sex stratification data prep first: 00n_sex_stratification_prep.py")
}

# Global temperature percentiles (consistent with all other scripts)
temp_all  <- era5$temp_mean
temp_pcts <- quantile(temp_all, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)
temp_boundary <- c(0, 40)

cat(sprintf("  GLOBAL knots (P10/P75/P90): %.1f / %.1f / %.1f\n",
            temp_pcts[1], temp_pcts[2], temp_pcts[3]))
cat(sprintf("  Boundary knots: %.0f to %.0f\n", temp_boundary[1], temp_boundary[2]))

# Temperature variability per region (for supplementary models)
temp_var <- as.data.table(read_parquet(era5_file))
temp_var <- temp_var[, .(temp_sd = sd(temp_mean, na.rm = TRUE)), by = region_code]

# -----------------------------------------------------------------------------
# 2. Fit Per-Region Sex-Specific DLNMs (saving per-region estimates)
# -----------------------------------------------------------------------------
cat("\n[2] Fitting per-region sex-specific DLNMs...\n")

fit_region_sex_dlnm <- function(mort_file, sex_label, era5_dt) {
  
  cat(sprintf("\n  --- %s ---\n", sex_label))
  
  mort <- as.data.table(read_parquet(mort_file))
  mort[, date := as.character(as.Date(date))]
  
  # Standardize region column
  region_col <- paste0(EXPOSURE_TYPE, "_code")
  if (region_col %in% names(mort)) {
    setnames(mort, region_col, "region_code")
  }
  
  # Merge
  setkey(mort, region_code, date)
  era5_copy <- copy(era5_dt)
  setkey(era5_copy, region_code, date)
  data <- merge(mort, era5_copy, by = c("region_code", "date"))
  data[, date := as.Date(date)]
  
  # Working variables
  if (!"deaths" %in% names(data)) {
    death_cols <- grep("death", names(data), value = TRUE, ignore.case = TRUE)
    if (length(death_cols) > 0) setnames(data, death_cols[1], "deaths")
  }
  data[, tmean := temp_mean]
  data <- data[!is.na(deaths) & !is.na(tmean)]
  setorder(data, region_code, date)
  
  regions <- unique(data$region_code)
  n_regions <- length(regions)
  
  # Storage for per-region results
  region_results <- list()
  
  for (i in seq_along(regions)) {
    reg <- regions[i]
    if (i %% 50 == 0) cat(sprintf("    Progress: %d/%d\n", i, n_regions))
    
    daily <- data[region_code == reg][order(date)]
    if (nrow(daily) < 365) next
    
    result <- tryCatch({
      reg_p1  <- quantile(daily$tmean, 0.01, na.rm = TRUE)
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
      
      # Find MMT
      mmt_search_min <- max(temp_boundary[1], reg_p1)
      mmt_search_max <- min(temp_boundary[2], reg_p99)
      temp_seq <- seq(mmt_search_min, mmt_search_max, length.out = 100)
      cp <- crosspred(cb, model, at = temp_seq, cumul = TRUE, cen = median(daily$tmean))
      mmt_idx <- which.min(cp$allRRfit)
      mmt <- temp_seq[mmt_idx]
      
      # RRs at P1 (cold) and P99 (heat)
      reg_pcts_vals <- quantile(daily$tmean, probs = c(0.01, 0.99), na.rm = TRUE)
      cp_pcts <- crosspred(cb, model, at = reg_pcts_vals, cumul = TRUE, cen = mmt)
      
      rr_p99    <- cp_pcts$allRRfit[2]
      rr_p99_ci <- c(cp_pcts$allRRlow[2], cp_pcts$allRRhigh[2])
      se_p99    <- (log(rr_p99_ci[2]) - log(rr_p99_ci[1])) / (2 * 1.96)
      
      rr_p1     <- cp_pcts$allRRfit[1]
      rr_p1_ci  <- c(cp_pcts$allRRlow[1], cp_pcts$allRRhigh[1])
      se_p1     <- (log(rr_p1_ci[2]) - log(rr_p1_ci[1])) / (2 * 1.96)
      
      # Quality check
      ok <- is.finite(se_p99) && is.finite(se_p1) && se_p99 < 5 && se_p1 < 5 &&
            se_p99 > 0 && se_p1 > 0
      
      if (ok) {
        list(
          region_code  = as.integer(reg),
          log_rr_heat  = log(rr_p99),
          se_heat      = se_p99,
          rr_heat      = rr_p99,
          rr_heat_lo   = rr_p99_ci[1],
          rr_heat_hi   = rr_p99_ci[2],
          log_rr_cold  = log(rr_p1),
          se_cold      = se_p1,
          rr_cold      = rr_p1,
          rr_cold_lo   = rr_p1_ci[1],
          rr_cold_hi   = rr_p1_ci[2],
          mmt          = mmt,
          n_days       = nrow(daily)
        )
      } else NULL
      
    }, error = function(e) NULL)
    
    if (!is.null(result)) {
      region_results[[as.character(reg)]] <- result
    }
  }
  
  # Convert to data.table
  dt <- rbindlist(region_results, fill = TRUE)
  dt[, sex := sex_label]
  
  cat(sprintf("    Successful: %d/%d regions (%.1f%%)\n",
              nrow(dt), n_regions, 100 * nrow(dt) / n_regions))
  
  return(dt)
}

# Fit male and female
male_dt   <- fit_region_sex_dlnm(sex_mort_files$male,   "male",   era5)
female_dt <- fit_region_sex_dlnm(sex_mort_files$female, "female", era5)

# Stack
stacked <- rbind(male_dt, female_dt)
cat(sprintf("\n  Stacked dataset: %d rows (%d male + %d female)\n",
            nrow(stacked), nrow(male_dt), nrow(female_dt)))

# -----------------------------------------------------------------------------
# 3. Merge Covariates
# -----------------------------------------------------------------------------
cat("\n[3] Merging covariates...\n")

# Merge HDI + macro_region from regional_covariates
stacked <- merge(stacked, reg_covariates[, .(region_code, hdi, macro_region, urban_pct,
                                              gdp_per_capita_brl, hospital_beds_per_1000,
                                              mean_temp_annual)],
                 by = "region_code", all.x = TRUE)

# Merge GDP per capita from SES covariates (more granular)
stacked <- merge(stacked, ses_covariates[, .(region_code, gdp_per_capita, urbanization_rate)],
                 by = "region_code", all.x = TRUE)

# Merge temperature variability
stacked <- merge(stacked, temp_var, by = "region_code", all.x = TRUE)

# Create deprivation proxy: inverse HDI (higher = more deprived)
# Use SES GDP if regional GDP is missing
stacked[, gdp := fifelse(is.na(gdp_per_capita_brl), gdp_per_capita, gdp_per_capita_brl)]

# Standardised moderators
stacked[, deprivation   := 1 - hdi]  # Higher = more deprived
stacked[, deprivation_z := (deprivation - mean(deprivation, na.rm=TRUE)) / sd(deprivation, na.rm=TRUE)]
stacked[, gdp_z         := (gdp - mean(gdp, na.rm=TRUE)) / sd(gdp, na.rm=TRUE)]
stacked[, urban_z        := (urbanization_rate - mean(urbanization_rate, na.rm=TRUE)) / sd(urbanization_rate, na.rm=TRUE)]
stacked[, temp_sd_z     := (temp_sd - mean(temp_sd, na.rm=TRUE)) / sd(temp_sd, na.rm=TRUE)]
stacked[, female        := as.integer(sex == "female")]

# Make macro_region a factor with Southeast as reference (largest population)
stacked[, macro_region := factor(macro_region,
                                  levels = c("Southeast", "South", "Northeast", "North", "Central-West"))]

# Drop rows with missing covariates
n_before <- nrow(stacked)
stacked <- stacked[!is.na(deprivation_z) & !is.na(macro_region)]
n_after <- nrow(stacked)
if (n_before > n_after) {
  cat(sprintf("  Dropped %d rows with missing covariates\n", n_before - n_after))
}

cat(sprintf("  Final dataset: %d rows\n", nrow(stacked)))
cat(sprintf("  HDI range: %.3f - %.3f\n", min(stacked$hdi, na.rm=TRUE), max(stacked$hdi, na.rm=TRUE)))
cat(sprintf("  Deprivation_z range: %.2f - %.2f\n", min(stacked$deprivation_z), max(stacked$deprivation_z)))
cat(sprintf("  Macro-regions: %s\n", paste(levels(stacked$macro_region), collapse=", ")))
cat(sprintf("  Sex: %d male, %d female\n", sum(stacked$female==0), sum(stacked$female==1)))

# -----------------------------------------------------------------------------
# 4. Meta-Regression Models
# -----------------------------------------------------------------------------
cat("\n[4] Running interaction meta-regression models...\n")

# Helper to extract coefficients safely
extract_results <- function(fit, outcome_label) {
  tryCatch({
    s <- summary(fit)
    coefs <- s$coefficients
    list(
      coefficients = as.data.frame(coefs),
      aic  = AIC(fit),
      bic  = BIC(fit),
      loglik = logLik(fit),
      converged = fit$converged,
      label = outcome_label
    )
  }, error = function(e) {
    cat(sprintf("    Warning extracting %s: %s\n", outcome_label, e$message))
    list(
      coefficients = data.frame(
        Estimate = coef(fit),
        Std.Error = tryCatch(sqrt(diag(vcov(fit))), error=function(e) rep(NA, length(coef(fit))))
      ),
      label = outcome_label
    )
  })
}

# Safe model fitting wrapper
safe_mixmeta <- function(formula, S_var, data, label) {
  cat(sprintf("\n  --- %s ---\n", label))
  fit <- tryCatch({
    mixmeta(formula, S = S_var, data = data, method = "reml")
  }, error = function(e) {
    cat(sprintf("    REML failed (%s), trying ML...\n", e$message))
    tryCatch({
      mixmeta(formula, S = S_var, data = data, method = "ml")
    }, error = function(e2) {
      cat(sprintf("    ML also failed: %s\n", e2$message))
      NULL
    })
  })
  
  if (!is.null(fit)) {
    res <- extract_results(fit, label)
    cat("    Coefficients:\n")
    print(res$coefficients, digits=4)
    return(list(fit = fit, results = res))
  }
  return(NULL)
}

all_results <- list()

# ---- HEAT MODELS ----
cat("\n========== HEAT (P99) MODELS ==========\n")

# Model H0: Null (intercept only)
all_results$heat_null <- safe_mixmeta(
  log_rr_heat ~ 1, S_var = stacked$se_heat^2, data = stacked,
  label = "H0: Heat null"
)

# Model H1: Sex only
all_results$heat_sex <- safe_mixmeta(
  log_rr_heat ~ female, S_var = stacked$se_heat^2, data = stacked,
  label = "H1: Heat ~ sex"
)

# Model H2: Sex + deprivation
all_results$heat_sex_dep <- safe_mixmeta(
  log_rr_heat ~ female + deprivation_z, S_var = stacked$se_heat^2, data = stacked,
  label = "H2: Heat ~ sex + deprivation"
)

# Model H3: Sex × deprivation (interaction)
all_results$heat_sex_x_dep <- safe_mixmeta(
  log_rr_heat ~ female * deprivation_z, S_var = stacked$se_heat^2, data = stacked,
  label = "H3: Heat ~ sex × deprivation"
)

# Model H4: Sex + deprivation + macro_region (additive)
all_results$heat_additive <- safe_mixmeta(
  log_rr_heat ~ female + deprivation_z + macro_region, S_var = stacked$se_heat^2, data = stacked,
  label = "H4: Heat ~ sex + deprivation + region"
)

# Model H5: Sex × deprivation + macro_region
all_results$heat_sexdep_region <- safe_mixmeta(
  log_rr_heat ~ female * deprivation_z + macro_region, S_var = stacked$se_heat^2, data = stacked,
  label = "H5: Heat ~ sex × deprivation + region"
)

# Model H6: Sex × macro_region + deprivation
all_results$heat_sexregion_dep <- safe_mixmeta(
  log_rr_heat ~ female * macro_region + deprivation_z, S_var = stacked$se_heat^2, data = stacked,
  label = "H6: Heat ~ sex × region + deprivation"
)

# Model H7: Full three-way interaction
all_results$heat_full <- safe_mixmeta(
  log_rr_heat ~ female * deprivation_z * macro_region, S_var = stacked$se_heat^2, data = stacked,
  label = "H7: Heat ~ sex × deprivation × region"
)

# ---- COLD MODELS ----
cat("\n\n========== COLD (P1) MODELS ==========\n")

# Model C0: Null
all_results$cold_null <- safe_mixmeta(
  log_rr_cold ~ 1, S_var = stacked$se_cold^2, data = stacked,
  label = "C0: Cold null"
)

# Model C1: Sex only
all_results$cold_sex <- safe_mixmeta(
  log_rr_cold ~ female, S_var = stacked$se_cold^2, data = stacked,
  label = "C1: Cold ~ sex"
)

# Model C2: Sex + deprivation
all_results$cold_sex_dep <- safe_mixmeta(
  log_rr_cold ~ female + deprivation_z, S_var = stacked$se_cold^2, data = stacked,
  label = "C2: Cold ~ sex + deprivation"
)

# Model C3: Sex × deprivation
all_results$cold_sex_x_dep <- safe_mixmeta(
  log_rr_cold ~ female * deprivation_z, S_var = stacked$se_cold^2, data = stacked,
  label = "C3: Cold ~ sex × deprivation"
)

# Model C4: Sex + deprivation + region
all_results$cold_additive <- safe_mixmeta(
  log_rr_cold ~ female + deprivation_z + macro_region, S_var = stacked$se_cold^2, data = stacked,
  label = "C4: Cold ~ sex + deprivation + region"
)

# Model C5: Sex × deprivation + region
all_results$cold_sexdep_region <- safe_mixmeta(
  log_rr_cold ~ female * deprivation_z + macro_region, S_var = stacked$se_cold^2, data = stacked,
  label = "C5: Cold ~ sex × deprivation + region"
)

# Model C6: Sex × region + deprivation
all_results$cold_sexregion_dep <- safe_mixmeta(
  log_rr_cold ~ female * macro_region + deprivation_z, S_var = stacked$se_cold^2, data = stacked,
  label = "C6: Cold ~ sex × region + deprivation"
)

# Model C7: Full three-way
all_results$cold_full <- safe_mixmeta(
  log_rr_cold ~ female * deprivation_z * macro_region, S_var = stacked$se_cold^2, data = stacked,
  label = "C7: Cold ~ sex × deprivation × region"
)

# -----------------------------------------------------------------------------
# 5. Model Comparison
# -----------------------------------------------------------------------------
cat("\n\n=======================================================\n")
cat("MODEL COMPARISON\n")
cat("=======================================================\n")

# Build comparison table
model_names <- c(
  "H0: Null", "H1: Sex", "H2: Sex+Dep", "H3: Sex×Dep",
  "H4: Additive", "H5: Sex×Dep+Region", "H6: Sex×Region+Dep", "H7: Full 3-way",
  "C0: Null", "C1: Sex", "C2: Sex+Dep", "C3: Sex×Dep",
  "C4: Additive", "C5: Sex×Dep+Region", "C6: Sex×Region+Dep", "C7: Full 3-way"
)

model_keys <- c(
  "heat_null", "heat_sex", "heat_sex_dep", "heat_sex_x_dep",
  "heat_additive", "heat_sexdep_region", "heat_sexregion_dep", "heat_full",
  "cold_null", "cold_sex", "cold_sex_dep", "cold_sex_x_dep",
  "cold_additive", "cold_sexdep_region", "cold_sexregion_dep", "cold_full"
)

comparison <- data.table(
  model    = model_names,
  outcome  = rep(c("heat", "cold"), each = 8),
  n_params = NA_real_,
  aic      = NA_real_,
  bic      = NA_real_,
  loglik   = NA_real_
)

for (i in seq_along(model_keys)) {
  key <- model_keys[i]
  if (!is.null(all_results[[key]]) && !is.null(all_results[[key]]$fit)) {
    fit <- all_results[[key]]$fit
    comparison$n_params[i] <- length(coef(fit))
    comparison$aic[i]      <- tryCatch(AIC(fit), error = function(e) NA)
    comparison$bic[i]      <- tryCatch(BIC(fit), error = function(e) NA)
    comparison$loglik[i]   <- tryCatch(as.numeric(logLik(fit)), error = function(e) NA)
  }
}

cat("\nHEAT models:\n")
print(comparison[outcome == "heat"], digits = 4)

cat("\nCOLD models:\n")
print(comparison[outcome == "cold"], digits = 4)

# Best model by AIC
best_heat <- comparison[outcome == "heat"][which.min(aic)]
best_cold <- comparison[outcome == "cold"][which.min(aic)]
cat(sprintf("\nBest heat model (AIC): %s (AIC=%.1f)\n", best_heat$model, best_heat$aic))
cat(sprintf("Best cold model (AIC): %s (AIC=%.1f)\n", best_cold$model, best_cold$aic))

# -----------------------------------------------------------------------------
# 6. Supplementary: Descriptive Table of Region-Sex Estimates
# -----------------------------------------------------------------------------
cat("\n[6] Region × Sex descriptive summary...\n")

# Mean RR by macro_region × sex
desc <- stacked[, .(
  n_regions     = .N,
  mean_rr_heat  = mean(rr_heat, na.rm=TRUE),
  median_rr_heat = median(rr_heat, na.rm=TRUE),
  mean_rr_cold  = mean(rr_cold, na.rm=TRUE),
  median_rr_cold = median(rr_cold, na.rm=TRUE),
  mean_hdi      = mean(hdi, na.rm=TRUE)
), by = .(macro_region, sex)]

setorder(desc, macro_region, sex)

cat("\nDescriptive: Mean RR by macro-region × sex\n")
cat("--------------------------------------------------------------\n")
print(desc, digits = 3)

# Also create deprivation tertiles for interpretability
stacked[, dep_tertile := cut(deprivation_z,
                              breaks = quantile(deprivation_z, probs = c(0, 1/3, 2/3, 1)),
                              labels = c("Low deprivation", "Medium", "High deprivation"),
                              include.lowest = TRUE)]

desc2 <- stacked[, .(
  n    = .N,
  mean_rr_heat = mean(rr_heat, na.rm=TRUE),
  mean_rr_cold = mean(rr_cold, na.rm=TRUE)
), by = .(dep_tertile, sex)]

setorder(desc2, dep_tertile, sex)
cat("\nDescriptive: Mean RR by deprivation tertile × sex\n")
cat("--------------------------------------------------------------\n")
print(desc2, digits = 3)

# Three-way cross-tabulation
desc3 <- stacked[, .(
  n = .N,
  mean_rr_heat = mean(rr_heat, na.rm=TRUE),
  mean_rr_cold = mean(rr_cold, na.rm=TRUE),
  mean_hdi     = mean(hdi, na.rm=TRUE)
), by = .(macro_region, sex, dep_tertile)]

setorder(desc3, macro_region, sex, dep_tertile)
cat("\nDescriptive: 3-way cross-tabulation (region × sex × deprivation)\n")
cat("--------------------------------------------------------------\n")
print(desc3, digits = 3)

# -----------------------------------------------------------------------------
# 7. Save Results
# -----------------------------------------------------------------------------
cat("\n[7] Saving results...\n")

# Prepare JSON-safe coefficient tables
json_coefs <- list()
for (key in model_keys) {
  if (!is.null(all_results[[key]]) && !is.null(all_results[[key]]$results)) {
    coef_df <- all_results[[key]]$results$coefficients
    json_coefs[[key]] <- list(
      label       = all_results[[key]]$results$label,
      coefficients = as.list(as.data.frame(t(coef_df))),
      term_names   = rownames(coef_df)
    )
  }
}

output <- list(
  exposure_type    = EXPOSURE_TYPE,
  analysis_date    = as.character(Sys.time()),
  n_regions        = uniqueN(stacked$region_code),
  n_obs            = nrow(stacked),
  deprivation_proxy = "inverse_HDI",
  
  model_comparison = as.list(comparison),
  best_heat_model  = as.list(best_heat),
  best_cold_model  = as.list(best_cold),
  
  coefficients     = json_coefs,
  
  descriptive_region_sex = as.list(desc),
  descriptive_dep_sex    = as.list(desc2),
  descriptive_3way       = as.list(desc3)
)

# Save JSON
output_file <- file.path(OUTPUT_DIR, paste0("interaction_metareg_", EXPOSURE_TYPE, ".json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 6)
cat("  JSON:", output_file, "\n")

# Save stacked dataset as CSV for further analysis / plotting
csv_file <- file.path(OUTPUT_DIR, paste0("interaction_metareg_data_", EXPOSURE_TYPE, ".csv"))
fwrite(stacked[, .(region_code, sex, female, log_rr_heat, se_heat, rr_heat, rr_heat_lo, rr_heat_hi,
                    log_rr_cold, se_cold, rr_cold, rr_cold_lo, rr_cold_hi,
                    mmt, n_days, hdi, macro_region, deprivation, deprivation_z,
                    gdp, gdp_z, urbanization_rate, urban_z, temp_sd, temp_sd_z, dep_tertile)],
       csv_file)
cat("  CSV:", csv_file, "\n")

cat("\n=======================================================\n")
cat("DONE!", as.character(Sys.time()), "\n")
cat("=======================================================\n")
