# Extract Results Summary from All Phases
library(jsonlite)

cat("\n")
cat("=============================================================\n")
cat("          TEMPERATURE-MORTALITY ANALYSIS RESULTS             \n")
cat("              Brazil Elderly (60+), 2010-2024                 \n")
cat("=============================================================\n")

# PHASE 1: DLNM
cat("\n### PHASE 1: DLNM (Main Analysis) ###\n\n")

for (level in c("intermediate", "immediate")) {
  f <- paste0("phase1_r/results/dlnm_r_", level, "_results_v2.json")
  if (file.exists(f)) {
    d <- fromJSON(f)
    pr <- d$pooled
    cat(toupper(level), "(", d$n_regions_success, "regions ):\n")
    cat("  Heat P99: RR =", round(pr$rr_heat_p99$rr, 3), 
        "(", round(pr$rr_heat_p99$rr_lo, 3), "-", round(pr$rr_heat_p99$rr_hi, 3), ")\n")
    cat("  Cold P1:  RR =", round(pr$rr_cold_p1$rr, 3),
        "(", round(pr$rr_cold_p1$rr_lo, 3), "-", round(pr$rr_cold_p1$rr_hi, 3), ")\n")
    cat("  MMT:", round(pr$mmt, 1), "C\n\n")
  }
}

# PHASE 1: Attributable Burden
cat("### PHASE 1: Attributable Burden ###\n\n")
for (level in c("intermediate", "immediate")) {
  f <- paste0("phase1_r/results/attributable_burden_r_", level, ".json")
  if (file.exists(f)) {
    d <- fromJSON(f)
    n <- d$national
    cat(toupper(level), ":\n")
    cat("  Total attributable deaths:", format(round(n$total_an), big.mark=","), "\n")
    cat("  Heat deaths (P99):", format(round(n$heat_an_99), big.mark=","), 
        "(", round(n$heat_af_pct_99, 2), "% AF)\n")
    cat("  Cold deaths (P1):", format(round(n$cold_an_1), big.mark=","),
        "(", round(n$cold_af_pct_1, 2), "% AF)\n")
    cat("  Annual heat deaths:", format(round(n$annual_heat_deaths), big.mark=","), "\n")
    cat("  Annual cold deaths:", format(round(n$annual_cold_deaths), big.mark=","), "\n\n")
  }
}

# PHASE 1: YLL
cat("### PHASE 1: Years of Life Lost ###\n\n")
for (level in c("intermediate", "immediate")) {
  f <- paste0("phase1_r/results/yll_r_", level, ".json")
  if (file.exists(f)) {
    d <- fromJSON(f)
    cat(toupper(level), ":\n")
    cat("  Total YLL:", format(round(d$summary$total_yll), big.mark=","), "\n")
    cat("  Heat YLL:", format(round(d$summary$heat_yll), big.mark=","), "\n")
    cat("  Cold YLL:", format(round(d$summary$cold_yll), big.mark=","), "\n\n")
  }
}

# PHASE 2: Sensitivity
cat("### PHASE 2: Sensitivity Analysis ###\n\n")
for (level in c("intermediate", "immediate")) {
  f <- paste0("phase2_r/results/sensitivity_r_", level, ".json")
  if (file.exists(f)) {
    d <- fromJSON(f)
    cat(toupper(level), ": Lag structures tested -", length(d$lag_sensitivity), "\n")
  }
}

# PHASE 2: Harvesting
cat("\n### PHASE 2: Harvesting Analysis ###\n\n")
for (level in c("intermediate", "immediate")) {
  f <- paste0("phase2_r/results/harvesting_r_", level, ".json")
  if (file.exists(f)) {
    d <- fromJSON(f)
    cat(toupper(level), ":\n")
    if (!is.null(d$harvesting_ratios)) {
      cat("  Lag-7 Heat RR:", round(d$pooled_by_horizon$lag_7$heat$rr, 3), "\n")
      cat("  Lag-21 Heat RR:", round(d$pooled_by_horizon$lag_21$heat$rr, 3), "\n")
      cat("  Lag-35 Heat RR:", round(d$pooled_by_horizon$lag_35$heat$rr, 3), "\n")
      cat("  Harvesting ratio (35/21):", round(d$harvesting_ratios$heat_35_21, 3), "\n\n")
    }
  }
}

# PHASE 4: Age Stratification
cat("### PHASE 4: Age Stratification ###\n\n")
for (level in c("intermediate", "immediate")) {
  f <- paste0("phase4_r/results/age_stratification_", level, ".json")
  if (file.exists(f)) {
    d <- fromJSON(f)
    cat(toupper(level), ":\n")
    for (age in names(d$results)) {
      r <- d$results[[age]]
      if (!is.null(r$heat$rr)) {
        cat("  ", age, ": Heat RR =", round(r$heat$rr, 3), 
            "(", round(r$heat$ci_low, 3), "-", round(r$heat$ci_high, 3), ")\n")
      }
    }
    cat("\n")
  }
}

# PHASE 4: Sex Stratification
cat("### PHASE 4: Sex Stratification ###\n\n")
for (level in c("intermediate", "immediate")) {
  f <- paste0("phase4_r/results/sex_stratification_", level, ".json")
  if (file.exists(f)) {
    d <- fromJSON(f)
    cat(toupper(level), ":\n")
    for (sex in names(d$results)) {
      r <- d$results[[sex]]
      if (!is.null(r$heat$rr)) {
        cat("  ", sex, ": Heat RR =", round(r$heat$rr, 3),
            "(", round(r$heat$ci_low, 3), "-", round(r$heat$ci_high, 3), ")\n")
      }
    }
    cat("\n")
  }
}

# PHASE 4: Cause Stratification
cat("### PHASE 4: Cause Stratification ###\n\n")
for (level in c("intermediate", "immediate")) {
  f <- paste0("phase4_r/results/cause_stratification_", level, ".json")
  if (file.exists(f)) {
    d <- fromJSON(f)
    cat(toupper(level), ":\n")
    for (cause in names(d$results)) {
      r <- d$results[[cause]]
      if (!is.null(r$heat$rr)) {
        cat("  ", cause, ": Heat RR =", round(r$heat$rr, 3),
            "(", round(r$heat$ci_low, 3), "-", round(r$heat$ci_high, 3), ")\n")
      }
    }
    cat("\n")
  }
}

cat("=============================================================\n")
cat("                      END OF SUMMARY                         \n")
cat("=============================================================\n")
