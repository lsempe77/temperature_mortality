# =============================================================================
# diagnose_extreme_se.R
# Investigate regions with extreme/infinite standard errors
# =============================================================================

suppressPackageStartupMessages({
  library(dlnm)
  library(data.table)
  library(arrow)
  library(splines)
})

EXPOSURE_TYPE <- "intermediate"
MAX_LAG <- 21
TEMP_DF <- 4
LAG_DF <- 4

cat("=======================================================\n")
cat("DIAGNOSING EXTREME SE REGIONS\n")
cat("=======================================================\n\n")

# Load data
DATA_DIR <- "../phase0_data_prep/results"
mort_file <- file.path(DATA_DIR, paste0("mortality_", EXPOSURE_TYPE, "_daily_elderly.parquet"))
era5_file <- file.path(DATA_DIR, paste0("era5_", EXPOSURE_TYPE, "_daily.parquet"))

mort <- as.data.table(read_parquet(mort_file))
era5 <- as.data.table(read_parquet(era5_file))

# Standardize column names
region_col <- paste0(EXPOSURE_TYPE, "_code")
if (region_col %in% names(mort)) setnames(mort, region_col, "region_code")
if (region_col %in% names(era5)) setnames(era5, region_col, "region_code")

# Merge
setkey(mort, region_code, date)
setkey(era5, region_code, date)
data <- merge(mort, era5, by = c("region_code", "date"))
data[, deaths := deaths_elderly]
data[, tmean := temp_mean]
data <- data[!is.na(deaths) & !is.na(tmean)]

regions <- unique(data$region_code)
cat(sprintf("Testing %d regions...\n\n", length(regions)))

# Store diagnostics
diagnostics <- list()

for (i in seq_along(regions)) {
  reg <- regions[i]
  
  if (i %% 20 == 0 || i == 1) {
    cat(sprintf("Processing region %d/%d (%s)...\n", i, length(regions), reg))
  }
  
  daily <- data[region_code == reg][order(date)]
  n_days <- nrow(daily)
  
  if (n_days < 365) next
  
  # Temperature characteristics
  temp_vector <- daily$tmean
  temp_range <- range(temp_vector, na.rm = TRUE)
  temp_mean <- mean(temp_vector, na.rm = TRUE)
  temp_sd <- sd(temp_vector, na.rm = TRUE)
  reg_temp_pcts <- quantile(temp_vector, probs = c(0.01, 0.10, 0.75, 0.90, 0.99), na.rm = TRUE)
  
  # Death characteristics
  deaths_total <- sum(daily$deaths)
  deaths_mean <- mean(daily$deaths)
  deaths_zero_pct <- sum(daily$deaths == 0) / n_days * 100
  
  result <- tryCatch({
    # Fit DLNM
    reg_boundary <- c(floor(temp_range[1]), ceiling(temp_range[2]))
    
    cb <- crossbasis(
      temp_vector,
      lag = MAX_LAG,
      argvar = list(fun = "ns", knots = reg_temp_pcts[c("10%", "75%", "90%")], 
                   Boundary.knots = reg_boundary),
      arglag = list(fun = "ns", df = LAG_DF)
    )
    
    formula_str <- "deaths ~ cb + ns(as.numeric(date), df = 7 * length(unique(format(date, '%Y')))) + factor(format(date, '%u'))"
    
    model <- glm(
      as.formula(formula_str),
      data = daily,
      family = quasipoisson(link = "log"),
      na.action = na.exclude
    )
    
    # Check convergence
    converged <- model$converged
    
    # Extract coefficients
    cb_idx <- grep("^cb", names(coef(model)))
    cb_coef <- coef(model)[cb_idx]
    cb_vcov <- vcov(model)[cb_idx, cb_idx]
    
    # Check vcov
    vcov_finite <- all(is.finite(cb_vcov))
    vcov_positive <- all(diag(cb_vcov) > 0)
    
    # Get predictions
    temp_seq <- seq(reg_boundary[1], reg_boundary[2], length.out = 100)
    cp_full <- crosspred(cb, model, at = temp_seq, cumul = TRUE, cen = median(temp_vector, na.rm=TRUE))
    mmt_idx <- which.min(cp_full$allRRfit)
    mmt <- temp_seq[mmt_idx]
    
    cp_extremes <- crosspred(cb, model, 
                             at = reg_temp_pcts[c("1%", "99%")], 
                             cumul = TRUE, cen = mmt)
    
    rr_p1 <- cp_extremes$allRRfit[1]
    rr_p1_se <- (cp_extremes$allRRhigh[1] - cp_extremes$allRRlow[1]) / (2 * 1.96)
    
    rr_p99 <- cp_extremes$allRRfit[2]
    rr_p99_se <- (cp_extremes$allRRhigh[2] - cp_extremes$allRRlow[2]) / (2 * 1.96)
    
    # Categorize SE
    se_category <- "valid"
    if (!is.finite(rr_p99_se) || !is.finite(rr_p1_se)) {
      se_category <- "infinite"
    } else if (rr_p99_se > 20 || rr_p1_se > 20) {
      se_category <- "extreme"
    } else if (rr_p99_se > 10 || rr_p1_se > 10) {
      se_category <- "high"
    }
    
    list(
      region = as.character(reg),
      n_days = n_days,
      temp_range = diff(temp_range),
      temp_mean = temp_mean,
      temp_sd = temp_sd,
      temp_p1 = reg_temp_pcts["1%"],
      temp_p99 = reg_temp_pcts["99%"],
      deaths_total = deaths_total,
      deaths_mean = deaths_mean,
      deaths_zero_pct = deaths_zero_pct,
      converged = converged,
      vcov_finite = vcov_finite,
      vcov_positive = vcov_positive,
      rr_p99 = rr_p99,
      rr_p99_se = rr_p99_se,
      rr_p1 = rr_p1,
      rr_p1_se = rr_p1_se,
      se_category = se_category,
      status = "success"
    )
  }, error = function(e) {
    list(
      region = as.character(reg),
      n_days = n_days,
      temp_range = diff(temp_range),
      temp_mean = temp_mean,
      temp_sd = temp_sd,
      temp_p1 = reg_temp_pcts["1%"],
      temp_p99 = reg_temp_pcts["99%"],
      deaths_total = deaths_total,
      deaths_mean = deaths_mean,
      deaths_zero_pct = deaths_zero_pct,
      converged = NA,
      vcov_finite = NA,
      vcov_positive = NA,
      rr_p99 = NA,
      rr_p99_se = NA,
      rr_p1 = NA,
      rr_p1_se = NA,
      se_category = "error",
      error_msg = e$message,
      status = "error"
    )
  })
  
  diagnostics[[i]] <- result
}

# Convert to data.table
diag_dt <- rbindlist(diagnostics, fill = TRUE)

cat("\n=======================================================\n")
cat("SUMMARY\n")
cat("=======================================================\n\n")

cat("SE Categories:\n")
print(table(diag_dt$se_category, useNA = "ifany"))

cat("\n\nRegions with extreme/infinite SE:\n")
extreme_regions <- diag_dt[se_category %in% c("extreme", "infinite", "high")]
if (nrow(extreme_regions) > 0) {
  print(extreme_regions[order(-rr_p99_se)][, .(
    region, n_days, temp_range, temp_sd, deaths_mean, deaths_zero_pct,
    rr_p99, rr_p99_se, rr_p1, rr_p1_se, se_category
  )])
  
  cat("\n\nDetailed characteristics of extreme SE regions:\n")
  cat("========================================\n")
  
  for (i in 1:min(5, nrow(extreme_regions))) {
    reg_data <- extreme_regions[i]
    cat(sprintf("\nRegion %s (SE category: %s):\n", reg_data$region, reg_data$se_category))
    cat(sprintf("  Days: %d\n", reg_data$n_days))
    cat(sprintf("  Temperature: mean=%.1f°C, sd=%.1f°C, range=%.1f°C\n", 
                reg_data$temp_mean, reg_data$temp_sd, reg_data$temp_range))
    cat(sprintf("  Extremes: P1=%.1f°C, P99=%.1f°C\n", reg_data$temp_p1, reg_data$temp_p99))
    cat(sprintf("  Deaths: total=%d, mean=%.2f/day, zero-days=%.1f%%\n",
                reg_data$deaths_total, reg_data$deaths_mean, reg_data$deaths_zero_pct))
    cat(sprintf("  Heat RR: %.3f (SE=%.3f)\n", reg_data$rr_p99, reg_data$rr_p99_se))
    cat(sprintf("  Cold RR: %.3f (SE=%.3f)\n", reg_data$rr_p1, reg_data$rr_p1_se))
    cat(sprintf("  Converged: %s, Finite vcov: %s, Positive vcov: %s\n",
                reg_data$converged, reg_data$vcov_finite, reg_data$vcov_positive))
  }
  
  # Compare with valid regions
  cat("\n\n=======================================================\n")
  cat("COMPARISON: Extreme vs Valid Regions\n")
  cat("=======================================================\n\n")
  
  valid_regions <- diag_dt[se_category == "valid"]
  
  cat("Temperature SD:\n")
  cat(sprintf("  Extreme: mean=%.2f, median=%.2f, range=[%.2f, %.2f]\n",
              mean(extreme_regions$temp_sd, na.rm=TRUE),
              median(extreme_regions$temp_sd, na.rm=TRUE),
              min(extreme_regions$temp_sd, na.rm=TRUE),
              max(extreme_regions$temp_sd, na.rm=TRUE)))
  cat(sprintf("  Valid: mean=%.2f, median=%.2f, range=[%.2f, %.2f]\n",
              mean(valid_regions$temp_sd, na.rm=TRUE),
              median(valid_regions$temp_sd, na.rm=TRUE),
              min(valid_regions$temp_sd, na.rm=TRUE),
              max(valid_regions$temp_sd, na.rm=TRUE)))
  
  cat("\nDeaths per day:\n")
  cat(sprintf("  Extreme: mean=%.2f, median=%.2f, range=[%.2f, %.2f]\n",
              mean(extreme_regions$deaths_mean, na.rm=TRUE),
              median(extreme_regions$deaths_mean, na.rm=TRUE),
              min(extreme_regions$deaths_mean, na.rm=TRUE),
              max(extreme_regions$deaths_mean, na.rm=TRUE)))
  cat(sprintf("  Valid: mean=%.2f, median=%.2f, range=[%.2f, %.2f]\n",
              mean(valid_regions$deaths_mean, na.rm=TRUE),
              median(valid_regions$deaths_mean, na.rm=TRUE),
              min(valid_regions$deaths_mean, na.rm=TRUE),
              max(valid_regions$deaths_mean, na.rm=TRUE)))
  
  cat("\nZero-death days (%):\n")
  cat(sprintf("  Extreme: mean=%.1f%%, median=%.1f%%, range=[%.1f%%, %.1f%%]\n",
              mean(extreme_regions$deaths_zero_pct, na.rm=TRUE),
              median(extreme_regions$deaths_zero_pct, na.rm=TRUE),
              min(extreme_regions$deaths_zero_pct, na.rm=TRUE),
              max(extreme_regions$deaths_zero_pct, na.rm=TRUE)))
  cat(sprintf("  Valid: mean=%.1f%%, median=%.1f%%, range=[%.1f%%, %.1f%%]\n",
              mean(valid_regions$deaths_zero_pct, na.rm=TRUE),
              median(valid_regions$deaths_zero_pct, na.rm=TRUE),
              min(valid_regions$deaths_zero_pct, na.rm=TRUE),
              max(valid_regions$deaths_zero_pct, na.rm=TRUE)))
  
} else {
  cat("No regions with extreme/infinite SE found!\n")
}

# Save results
output_file <- "results/se_diagnostics.csv"
fwrite(diag_dt, output_file)
cat(sprintf("\n\nDiagnostics saved to: %s\n", output_file))

cat("\n=======================================================\n")
cat("Done!\n")
cat("=======================================================\n")
