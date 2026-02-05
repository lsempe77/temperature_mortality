# Debug post-pandemic mixmeta failure
# =============================================================================

suppressPackageStartupMessages({
  library(dlnm)
  library(mixmeta)
  library(data.table)
  library(arrow)
  library(jsonlite)
  library(splines)
})

# Load data
DATA_DIR <- "phase0_data_prep/results"
DLNM_DIR <- "phase1_r/results"

cat("Loading DLNM results...\n")
dlnm_results <- fromJSON(file.path(DLNM_DIR, "dlnm_r_intermediate_results_v2.json"))

# Use ADAPTIVE knots based on the subset data, not fixed knots from full period
MAX_LAG <- dlnm_results$dlnm_params$max_lag
LAG_DF <- dlnm_results$dlnm_params$lag_df
pooled_mmt <- dlnm_results$pooled$mmt

cat("Loading panel data...\n")
# Load mortality and weather separately and merge
mort_file <- file.path(DATA_DIR, "mortality_intermediate_daily_elderly.parquet")
weather_file <- file.path(DATA_DIR, "era5_intermediate_daily.parquet")

mort <- as.data.table(read_parquet(mort_file))
weather <- as.data.table(read_parquet(weather_file))

# CRITICAL: Convert dates to same format before merging
mort$date <- as.Date(mort$date)
weather$date <- as.Date(weather$date)

# Merge on region_code and date
data <- merge(mort, weather, by = c("region_code", "date"), all.x = TRUE)
data[, year := year(date)]
data[, doy := yday(date)]
data[, dow := wday(date)]
data[, region_id := region_code]
data[, deaths := deaths_elderly]
data[, tmean := temp_mean]

cat("Total rows:", nrow(data), "\n")
cat("NA in tmean:", sum(is.na(data$tmean)), "\n")

# Filter to post-pandemic
data_post <- data[year >= 2022 & year <= 2024]
cat("Post-pandemic data:", nrow(data_post), "rows\n")
cat("Regions:", length(unique(data_post$region_id)), "\n")

# Calculate ADAPTIVE knots from the POST-PANDEMIC data
all_temps <- data_post$tmean[!is.na(data_post$tmean)]
TEMP_BOUNDARY_ADAPTIVE <- quantile(all_temps, c(0.01, 0.99), na.rm = TRUE)
TEMP_KNOTS_ADAPTIVE <- quantile(all_temps, c(0.25, 0.50, 0.75), na.rm = TRUE)

cat("\nOriginal knots from full period:\n")
cat("  Boundary:", dlnm_results$dlnm_params$temp_boundary, "\n")
cat("  Knots:", dlnm_results$dlnm_params$temp_knots, "\n")

cat("\nAdaptive knots for 2022-2024:\n")
cat("  Boundary:", TEMP_BOUNDARY_ADAPTIVE, "\n")
cat("  Knots:", TEMP_KNOTS_ADAPTIVE, "\n")

# Filter to post-pandemic
data_post <- data[year >= 2022 & year <= 2024]
cat("Post-pandemic data:", nrow(data_post), "rows\n")
cat("Regions:", length(unique(data_post$region_id)), "\n")

# Run first stage for ALL regions
regions <- unique(data_post$region_id)
first_stage_results <- list()

cat("\nRunning first stage models...\n")
for (i in seq_along(regions)) {
  reg <- regions[i]
  reg_data <- data_post[region_id == reg]
  
  if (nrow(reg_data) < 100) {
    cat("Region", reg, "skipped (too few rows)\n")
    next
  }
  
  tryCatch({
    cb <- crossbasis(
      reg_data$tmean,
      lag = MAX_LAG,
      argvar = list(fun = "ns", knots = TEMP_KNOTS_ADAPTIVE, Boundary.knots = TEMP_BOUNDARY_ADAPTIVE),
      arglag = list(fun = "ns", df = LAG_DF)
    )
    
    model <- glm(deaths ~ cb + ns(doy, df = 4) + factor(dow),
                 family = quasipoisson(), data = reg_data)
    
    red <- crossreduce(cb, model, cen = pooled_mmt)
    first_stage_results[[as.character(reg)]] <- list(coef = coef(red), vcov = vcov(red))
    
    if (i %% 20 == 0) cat("  Processed", i, "of", length(regions), "\n")
  }, error = function(e) {
    cat("Region", reg, "FAILED:", e$message, "\n")
  })
}

cat("\nFirst stage results:", length(first_stage_results), "converged\n")

# Now try mixmeta
if (length(first_stage_results) >= 3) {
  coef_mat <- do.call(rbind, lapply(first_stage_results, function(x) x$coef))
  vcov_list <- lapply(first_stage_results, function(x) x$vcov)
  
  cat("\nCoef matrix dim:", dim(coef_mat), "\n")
  cat("Number of vcov matrices:", length(vcov_list), "\n")
  
  # Check for issues in vcov matrices
  problematic <- c()
  for (i in seq_along(vcov_list)) {
    v <- vcov_list[[i]]
    eig <- eigen(v)$values
    if (any(eig <= 0)) {
      problematic <- c(problematic, i)
      cat("Region", names(first_stage_results)[i], "has non-positive eigenvalues:", min(eig), "\n")
    }
    if (any(is.na(v))) {
      problematic <- c(problematic, i)
      cat("Region", names(first_stage_results)[i], "has NA in vcov\n")
    }
  }
  
  cat("\nProblematic regions:", length(problematic), "\n")
  
  # Try mixmeta with error handling
  cat("\n=== Attempting mixmeta (REML) ===\n")
  tryCatch({
    mv <- mixmeta(coef_mat, vcov_list, method = "reml")
    cat("SUCCESS! mixmeta converged\n")
    print(summary(mv))
  }, error = function(e) {
    cat("MIXMETA REML FAILED:", e$message, "\n")
  })
  
  # Try fixed effects
  cat("\n=== Trying fixed effects (method=fixed) ===\n")
  tryCatch({
    mv_fixed <- mixmeta(coef_mat, vcov_list, method = "fixed")
    cat("Fixed effects SUCCESS!\n")
    print(summary(mv_fixed))
  }, error = function(e) {
    cat("Fixed effects FAILED:", e$message, "\n")
  })
  
  # Try ML instead of REML
  cat("\n=== Trying ML instead of REML ===\n")
  tryCatch({
    mv_ml <- mixmeta(coef_mat, vcov_list, method = "ml")
    cat("ML SUCCESS!\n")
    print(summary(mv_ml))
  }, error = function(e) {
    cat("ML FAILED:", e$message, "\n")
  })
  
  # Try removing problematic regions
  if (length(problematic) > 0) {
    cat("\n=== Trying after removing problematic regions ===\n")
    good_idx <- setdiff(seq_along(vcov_list), problematic)
    coef_mat_clean <- coef_mat[good_idx, , drop = FALSE]
    vcov_list_clean <- vcov_list[good_idx]
    
    cat("Regions remaining:", length(vcov_list_clean), "\n")
    
    tryCatch({
      mv_clean <- mixmeta(coef_mat_clean, vcov_list_clean, method = "reml")
      cat("SUCCESS after cleaning!\n")
      print(summary(mv_clean))
    }, error = function(e) {
      cat("Still FAILED:", e$message, "\n")
    })
  }
}

cat("\n=== DONE ===\n")
