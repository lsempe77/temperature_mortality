# =============================================================================
# 01c_yll_calculation.R
# Years of Life Lost (YLL) Calculation using R DLNM Attributable Burden
# =============================================================================
#
# Method:
# YLL = Attributable Deaths × Average Life Expectancy at Death
#
# For elderly population, we use:
# 1. IBGE Brazilian life tables (age-specific)
# 2. Assumed age distribution for elderly (60-69: 50%, 70-79: 35%, 80+: 15%)
# 3. Calculate weighted average YLL per death
#
# References:
# - WHO YLL methodology
# - GBD YLL approach
# - IBGE Brazilian life tables
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# Get exposure type from command line args
args <- commandArgs(trailingOnly = TRUE)
EXPOSURE_TYPE <- if (length(args) > 0) args[1] else "intermediate"

cat("=======================================================\n")
cat("YEARS OF LIFE LOST (YLL) CALCULATION:", EXPOSURE_TYPE, "\n")
cat("=======================================================\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# -----------------------------------------------------------------------------
# 1. Load Life Table
# -----------------------------------------------------------------------------
# Get script directory for proper path resolution
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
BURDEN_DIR <- file.path(SCRIPT_DIR, "results")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Script directory:", SCRIPT_DIR, "\n")
cat("Data directory:", normalizePath(DATA_DIR, mustWork = FALSE), "\n")

cat("[1] Loading life table...\n")

# Load IBGE life table
lt_file <- file.path(DATA_DIR, "yll_lookup_by_age.csv")
if (file.exists(lt_file)) {
  life_table <- fread(lt_file)
  cat("  Loaded IBGE life table:", nrow(life_table), "ages\n")
} else {
  # Fallback to combined life tables
  lt_file_alt <- file.path(DATA_DIR, "ibge_life_tables_combined.csv")
  if (file.exists(lt_file_alt)) {
    life_table <- fread(lt_file_alt)
    cat("  Loaded IBGE combined life tables\n")
  } else {
    # Create WHO GBD 2019 reference fallback
    cat("  WARNING: Using WHO GBD 2019 reference (fallback)\n")
    life_table <- data.table(
      age = seq(0, 95, by = 5),
      ex_mean = c(88.9, 84.0, 79.0, 74.1, 69.1, 64.1, 59.2, 54.2, 49.3, 44.4,
                  39.5, 34.7, 30.0, 25.5, 21.2, 17.2, 13.6, 10.5, 7.9, 5.8)
    )
  }
}

# Function to get life expectancy at a specific age
get_life_expectancy <- function(age_val) {
  if (age_val %in% life_table$age) {
    return(life_table[age == age_val, ex_mean])
  }
  # Interpolate
  if (age_val > max(life_table$age)) {
    return(life_table[age == max(age), ex_mean])
  }
  if (age_val < min(life_table$age)) {
    return(life_table[age == min(age), ex_mean])
  }
  return(approx(life_table$age, life_table$ex_mean, xout = age_val)$y)
}

# Life expectancy at key elderly ages
le_60 <- get_life_expectancy(60)
le_65 <- get_life_expectancy(65)
le_70 <- get_life_expectancy(70)
le_75 <- get_life_expectancy(75)
le_80 <- get_life_expectancy(80)
le_85 <- get_life_expectancy(85)
le_90 <- get_life_expectancy(90)

cat("\n  Life expectancy at key ages:\n")
cat("    60 years:", round(le_60, 1), "years\n")
cat("    70 years:", round(le_70, 1), "years\n")
cat("    80 years:", round(le_80, 1), "years\n")
cat("    90 years:", round(le_90, 1), "years\n")

# -----------------------------------------------------------------------------
# 2. Define Age Distribution for Elderly
# -----------------------------------------------------------------------------
cat("\n[2] Setting up age distribution...\n")

# Assumed age distribution for elderly deaths in Brazil (based on literature)
# More refined distribution based on typical patterns
age_dist <- data.table(
  age_group = c("60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90+"),
  mid_age = c(62, 67, 72, 77, 82, 87, 92),
  weight = c(0.15, 0.18, 0.20, 0.17, 0.15, 0.10, 0.05)  # Weights sum to 1
)

# Get life expectancy for each age group
age_dist[, life_exp := sapply(mid_age, get_life_expectancy)]

# Calculate weighted average life expectancy
weighted_avg_le <- sum(age_dist$weight * age_dist$life_exp)

cat("  Age distribution weights:\n")
print(age_dist[, .(age_group, weight, life_exp = round(life_exp, 1))])
cat("\n  Weighted average LE for elderly deaths:", round(weighted_avg_le, 2), "years\n")

# Alternative: use median elderly age (75) for simpler calculation
median_elderly_le <- get_life_expectancy(75)
cat("  LE at median elderly age (75):", round(median_elderly_le, 2), "years\n")

# -----------------------------------------------------------------------------
# 3. Load Attributable Burden Results
# -----------------------------------------------------------------------------
cat("\n[3] Loading attributable burden results...\n")

burden_file <- file.path(BURDEN_DIR, paste0("attributable_burden_r_", EXPOSURE_TYPE, ".json"))
if (!file.exists(burden_file)) {
  stop("Attributable burden results not found: ", burden_file, "\nRun 01b_attributable_burden.R first!")
}

burden_results <- fromJSON(burden_file)
cat("  Loaded:", burden_file, "\n")

# Extract key values
national <- burden_results$national
total_heat_an <- national$total_heat_an
total_cold_an <- national$total_cold_an
total_an <- national$total_an
n_years <- national$n_years

cat("\n  Attributable deaths:\n")
cat("    Heat:", format(round(total_heat_an), big.mark = ","), "\n")
cat("    Cold:", format(round(total_cold_an), big.mark = ","), "\n")
cat("    Total:", format(round(total_an), big.mark = ","), "\n")

# -----------------------------------------------------------------------------
# 4. Calculate YLL
# -----------------------------------------------------------------------------
cat("\n[4] Calculating Years of Life Lost...\n")

# Method 1: Using weighted average life expectancy
yll_heat_weighted <- total_heat_an * weighted_avg_le
yll_cold_weighted <- total_cold_an * weighted_avg_le
yll_total_weighted <- total_an * weighted_avg_le

# Method 2: Using median age (75) life expectancy
yll_heat_median <- total_heat_an * median_elderly_le
yll_cold_median <- total_cold_an * median_elderly_le
yll_total_median <- total_an * median_elderly_le

# Annual rates
annual_yll_heat <- yll_heat_weighted / n_years
annual_yll_cold <- yll_cold_weighted / n_years
annual_yll_total <- yll_total_weighted / n_years

# Per death YLL
yll_per_heat_death <- weighted_avg_le
yll_per_cold_death <- weighted_avg_le

# Per 100k population
pop_elderly <- national$total_pop_elderly
if (!is.null(pop_elderly) && pop_elderly > 0) {
  yll_rate_heat <- annual_yll_heat / pop_elderly * 100000
  yll_rate_cold <- annual_yll_cold / pop_elderly * 100000
  yll_rate_total <- annual_yll_total / pop_elderly * 100000
} else {
  yll_rate_heat <- NA
  yll_rate_cold <- NA
  yll_rate_total <- NA
}

# -----------------------------------------------------------------------------
# 5. Load Region-Level Burden for Regional YLL
# -----------------------------------------------------------------------------
cat("\n[5] Calculating regional YLL...\n")

region_csv <- file.path(BURDEN_DIR, paste0("attributable_burden_r_", EXPOSURE_TYPE, "_regions.csv"))
region_burden <- fread(region_csv)

# Calculate YLL by region
region_burden[, `:=`(
  yll_heat = total_heat_an * weighted_avg_le,
  yll_cold = total_cold_an * weighted_avg_le,
  yll_total = (total_heat_an + total_cold_an) * weighted_avg_le
)]

region_burden[, `:=`(
  annual_yll_heat = yll_heat / n_years,
  annual_yll_cold = yll_cold / n_years,
  annual_yll_total = yll_total / n_years
)]

# Per 100k rates
region_burden[, `:=`(
  yll_rate_heat = fifelse(pop_elderly > 0, annual_yll_heat / pop_elderly * 100000, NA_real_),
  yll_rate_cold = fifelse(pop_elderly > 0, annual_yll_cold / pop_elderly * 100000, NA_real_),
  yll_rate_total = fifelse(pop_elderly > 0, annual_yll_total / pop_elderly * 100000, NA_real_)
)]

cat("  Processed", nrow(region_burden), "regions\n")

# -----------------------------------------------------------------------------
# 6. Print Results
# -----------------------------------------------------------------------------
cat("\n=======================================================\n")
cat("YEARS OF LIFE LOST RESULTS (", EXPOSURE_TYPE, ")\n", sep = "")
cat("=======================================================\n")

cat("\n--- METHODOLOGY ---\n")
cat("Life expectancy source: IBGE Brazilian life tables\n")
cat("Weighted average LE for elderly:", round(weighted_avg_le, 2), "years\n")
cat("Study period:", round(n_years, 1), "years\n")

cat("\n--- TOTAL YLL (Weighted Age Distribution) ---\n")
cat("Heat YLL:", format(round(yll_heat_weighted), big.mark = ","), "\n")
cat("Cold YLL:", format(round(yll_cold_weighted), big.mark = ","), "\n")
cat("Total YLL:", format(round(yll_total_weighted), big.mark = ","), "\n")

cat("\n--- ANNUAL YLL ---\n")
cat("Heat YLL/year:", format(round(annual_yll_heat), big.mark = ","), "\n")
cat("Cold YLL/year:", format(round(annual_yll_cold), big.mark = ","), "\n")
cat("Total YLL/year:", format(round(annual_yll_total), big.mark = ","), "\n")

cat("\n--- YLL RATES (per 100,000 elderly) ---\n")
if (!is.na(yll_rate_total)) {
  cat("Heat YLL rate:", round(yll_rate_heat, 1), "per 100k/year\n")
  cat("Cold YLL rate:", round(yll_rate_cold, 1), "per 100k/year\n")
  cat("Total YLL rate:", round(yll_rate_total, 1), "per 100k/year\n")
} else {
  cat("Population data not available for rate calculation\n")
}

cat("\n--- COMPARISON: MEDIAN AGE METHOD ---\n")
cat("Heat YLL (median):", format(round(yll_heat_median), big.mark = ","), "\n")
cat("Cold YLL (median):", format(round(yll_cold_median), big.mark = ","), "\n")
cat("Total YLL (median):", format(round(yll_total_median), big.mark = ","), "\n")

cat("\n--- REGIONAL DISTRIBUTION ---\n")
cat("  Top 5 regions by total YLL:\n")
top5 <- region_burden[order(-yll_total)][1:5]
print(top5[, .(region_code, total_deaths, yll_total = round(yll_total), 
               annual_yll = round(annual_yll_total))])

# -----------------------------------------------------------------------------
# 7. Save Results
# -----------------------------------------------------------------------------
cat("\n[7] Saving results...\n")

# Prepare output
output <- list(
  exposure_type = EXPOSURE_TYPE,
  analysis_date = as.character(Sys.time()),
  
  methodology = list(
    life_table_source = "IBGE Brazilian life tables",
    age_distribution = as.data.frame(age_dist),
    weighted_avg_le = weighted_avg_le,
    median_elderly_le = median_elderly_le
  ),
  
  national = list(
    n_years = n_years,
    total_deaths = national$total_deaths,
    pop_elderly = pop_elderly,
    
    # Attributable deaths
    heat_attributable = total_heat_an,
    cold_attributable = total_cold_an,
    total_attributable = total_an,
    
    # YLL (weighted method)
    yll_heat = yll_heat_weighted,
    yll_cold = yll_cold_weighted,
    yll_total = yll_total_weighted,
    
    # Annual YLL
    annual_yll_heat = annual_yll_heat,
    annual_yll_cold = annual_yll_cold,
    annual_yll_total = annual_yll_total,
    
    # YLL rates per 100k
    yll_rate_heat_per_100k = yll_rate_heat,
    yll_rate_cold_per_100k = yll_rate_cold,
    yll_rate_total_per_100k = yll_rate_total,
    
    # YLL per death
    yll_per_death = weighted_avg_le,
    
    # Comparison method (median age)
    yll_heat_median_method = yll_heat_median,
    yll_cold_median_method = yll_cold_median,
    yll_total_median_method = yll_total_median
  ),
  
  regional_summary = list(
    n_regions = nrow(region_burden),
    total_regional_yll = sum(region_burden$yll_total),
    mean_regional_yll = mean(region_burden$yll_total),
    median_regional_yll = median(region_burden$yll_total)
  )
)

# Save JSON
output_file <- file.path(OUTPUT_DIR, paste0("yll_r_", EXPOSURE_TYPE, ".json"))
write_json(output, output_file, auto_unbox = TRUE, pretty = TRUE, digits = 4)
cat("  JSON saved to:", output_file, "\n")

# Save regional CSV
region_yll <- region_burden[, .(
  region_code, n_days, n_years, total_deaths, pop_elderly,
  heat_an = total_heat_an, cold_an = total_cold_an,
  yll_heat, yll_cold, yll_total,
  annual_yll_heat, annual_yll_cold, annual_yll_total,
  yll_rate_heat, yll_rate_cold, yll_rate_total
)]

csv_file <- file.path(OUTPUT_DIR, paste0("yll_r_", EXPOSURE_TYPE, "_regions.csv"))
fwrite(region_yll, csv_file)
cat("  Regional CSV saved to:", csv_file, "\n")

cat("\n=======================================================\n")
cat("Done!", as.character(Sys.time()), "\n")
cat("=======================================================\n")
