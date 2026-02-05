# Extract and summarize all results from Phase 1-4
library(jsonlite)

cat("============================================================\n")
cat("COMPLETE RESULTS EXTRACTION - ALL PHASES\n")
cat("============================================================\n")
cat("Date:", as.character(Sys.time()), "\n\n")

# ============================================
# PHASE 1: DLNM CORE
# ============================================
cat("============================================\n")
cat("PHASE 1: DLNM CORE RESULTS\n")
cat("============================================\n")

# DLNM
for (level in c("intermediate", "immediate")) {
  d <- fromJSON(sprintf("phase1_r/results/dlnm_r_%s_results_v2.json", level))
  cat(sprintf("\n[%s DLNM]\n", toupper(level)))
  cat(sprintf("  Regions: %d total, %d success, %d in meta (%.1f%%)\n", 
              d$n_regions_total, d$n_regions_success, d$n_regions_mvmeta,
              100*d$n_regions_mvmeta/d$n_regions_total))
  cat(sprintf("  Pooled MMT: %.2f C\n", d$pooled$mmt))
  cat(sprintf("  Heat RR (P99): %.3f (%.3f-%.3f)\n", 
              d$pooled$rr_heat_p99$rr, d$pooled$rr_heat_p99$rr_lo, d$pooled$rr_heat_p99$rr_hi))
  cat(sprintf("  Cold RR (P1): %.3f (%.3f-%.3f)\n",
              d$pooled$rr_cold_p1$rr, d$pooled$rr_cold_p1$rr_lo, d$pooled$rr_cold_p1$rr_hi))
  cat(sprintf("  I2: %.1f%%, Converged: %s\n", 
              d$pooled$heterogeneity$I2_percent, d$pooled$converged))
}

# Attributable Burden
cat("\n--- Attributable Burden ---\n")
for (level in c("intermediate", "immediate")) {
  d <- fromJSON(sprintf("phase1_r/results/attributable_burden_r_%s.json", level))
  cat(sprintf("[%s] Total AF: %.2f%% (Heat: %.2f%%, Cold: %.2f%%)\n",
              toupper(level), d$national$total_af_pct, 
              d$national$total_heat_af_pct, d$national$total_cold_af_pct))
  cat(sprintf("         Annual deaths: Heat=%.0f, Cold=%.0f, Total=%.0f\n",
              d$national$annual_heat_deaths, d$national$annual_cold_deaths,
              d$national$annual_heat_deaths + d$national$annual_cold_deaths))
}

# YLL
cat("\n--- Years of Life Lost ---\n")
for (level in c("intermediate", "immediate")) {
  d <- fromJSON(sprintf("phase1_r/results/yll_r_%s.json", level))
  cat(sprintf("[%s] Annual YLL: Heat=%.0f, Cold=%.0f, Total=%.0f\n",
              toupper(level), d$national$annual_yll_heat, 
              d$national$annual_yll_cold, d$national$annual_yll_total))
}

# Case-Crossover
cat("\n--- Case-Crossover Validation ---\n")
for (level in c("intermediate", "immediate")) {
  d <- fromJSON(sprintf("phase1_r/results/case_crossover_r_%s.json", level))
  cat(sprintf("[%s] Linear OR/C: %.3f, Extreme Heat OR: %.3f, Extreme Cold OR: %.3f\n",
              toupper(level), 
              d$simple_models$linear$or_per_degree,
              d$simple_models$categorical_p99_p1$or_extreme_heat,
              d$simple_models$categorical_p99_p1$or_extreme_cold))
  cat(sprintf("         DLNM converged: %s\n", d$dlnm_case_crossover$converged))
}

# Excess Mortality
cat("\n--- Excess Mortality ---\n")
for (level in c("intermediate", "immediate")) {
  d <- fromJSON(sprintf("phase1_r/results/excess_mortality_r_%s.json", level))
  cat(sprintf("[%s] Model deviance explained: %.1f%%\n", toupper(level), 
              d$model$deviance_explained * 100))
  cat(sprintf("         Extreme heat excess: %.0f deaths (+%.1f%%)\n",
              d$excess_summary$extreme_heat$total_excess,
              d$excess_summary$extreme_heat$excess_pct))
}

# ============================================
# PHASE 2: ROBUSTNESS
# ============================================
cat("\n============================================\n")
cat("PHASE 2: ROBUSTNESS ANALYSES\n")
cat("============================================\n")

# Sensitivity
cat("\n--- Lag Sensitivity ---\n")
for (level in c("intermediate", "immediate")) {
  d <- fromJSON(sprintf("phase2_r/results/sensitivity_r_%s.json", level))
  cat(sprintf("[%s] Lag 21 (baseline): Heat RR=%.3f, Cold RR=%.3f\n",
              toupper(level), d$lag_sensitivity$lag_21$rr_p99, d$lag_sensitivity$lag_21$rr_p1))
}

# Harvesting
cat("\n--- Harvesting Analysis ---\n")
for (level in c("intermediate", "immediate")) {
  d <- fromJSON(sprintf("phase2_r/results/harvesting_r_%s.json", level))
  cat(sprintf("[%s] Heat harvesting ratio: %.2f, Cold: %.2f\n",
              toupper(level), 
              d$results$summary$harvesting_ratio_heat,
              d$results$summary$harvesting_ratio_cold))
  if (d$results$summary$harvesting_ratio_heat > 1) {
    cat("         -> Heat shows mortality displacement (harvesting)\n")
  }
  if (abs(d$results$summary$harvesting_ratio_cold) < 0.5) {
    cat("         -> Cold represents true excess mortality\n")
  }
}

# Heatwave
cat("\n--- Heatwave Analysis ---\n")
for (level in c("intermediate", "immediate")) {
  d <- fromJSON(sprintf("phase2_r/results/heatwave_r_%s.json", level))
  cat(sprintf("[%s] Heatwave RR: %.3f (%.3f-%.3f), N regions: %d\n",
              toupper(level), d$results$heatwave_rr, 
              d$results$heatwave_ci_low, d$results$heatwave_ci_high,
              d$results$n_regions))
}

# ============================================
# PHASE 3: CONFOUNDING
# ============================================
cat("\n============================================\n")
cat("PHASE 3: CONFOUNDING CONTROL\n")
cat("============================================\n")

for (level in c("intermediate", "immediate")) {
  d <- fromJSON(sprintf("phase3_r/results/supplementary_r_%s.json", level))
  cat(sprintf("\n[%s SUPPLEMENTARY]\n", toupper(level)))
  cat(sprintf("  Baseline: Heat RR=%.3f, Cold RR=%.3f (n=%d regions)\n",
              d$baseline$heat$rr, d$baseline$cold$rr, d$baseline$n_regions))
  
  if (!is.null(d$apparent_temp) && !is.null(d$apparent_temp$heat$rr)) {
    heat_change <- (d$apparent_temp$heat$rr - d$baseline$heat$rr) / d$baseline$heat$rr * 100
    cold_change <- (d$apparent_temp$cold$rr - d$baseline$cold$rr) / d$baseline$cold$rr * 100
    cat(sprintf("  Apparent Temp: Heat %.1f%% change, Cold %.1f%% change -> %s\n",
                heat_change, cold_change,
                ifelse(abs(heat_change) < 10 & abs(cold_change) < 10, "ROBUST", "SENSITIVE")))
  }
  
  if (!is.null(d$pollution_adjusted) && !is.null(d$pollution_adjusted$heat$rr)) {
    heat_change <- (d$pollution_adjusted$heat$rr - d$baseline$heat$rr) / d$baseline$heat$rr * 100
    cold_change <- (d$pollution_adjusted$cold$rr - d$baseline$cold$rr) / d$baseline$cold$rr * 100
    cat(sprintf("  Pollution Adj: Heat %.1f%% change, Cold %.1f%% change -> %s\n",
                heat_change, cold_change,
                ifelse(abs(heat_change) < 10 & abs(cold_change) < 10, "ROBUST", "SENSITIVE")))
  }
  
  if (!is.null(d$flu_adjusted) && !is.null(d$flu_adjusted$heat$rr)) {
    heat_change <- (d$flu_adjusted$heat$rr - d$baseline$heat$rr) / d$baseline$heat$rr * 100
    cold_change <- (d$flu_adjusted$cold$rr - d$baseline$cold$rr) / d$baseline$cold$rr * 100
    cat(sprintf("  Flu Adjusted: Heat %.1f%% change, Cold %.1f%% change -> %s\n",
                heat_change, cold_change,
                ifelse(abs(heat_change) < 10 & abs(cold_change) < 10, "ROBUST", "SENSITIVE")))
  }
}

# ============================================
# PHASE 4: HETEROGENEITY
# ============================================
cat("\n============================================\n")
cat("PHASE 4: HETEROGENEITY & EFFECT MODIFICATION\n")
cat("============================================\n")

# Meta-regression
cat("\n--- Meta-Regression (Effect Modifiers) ---\n")
for (level in c("intermediate", "immediate")) {
  d <- fromJSON(sprintf("phase4_r/results/meta_regression_%s.json", level))
  cat(sprintf("\n[%s] N regions: %d\n", toupper(level), d$n_regions))
  
  mods <- names(d$moderators)
  for (mod in mods) {
    heat_p <- d$moderators[[mod]]$heat$p_value
    cold_p <- d$moderators[[mod]]$cold$p_value
    heat_sig <- ifelse(heat_p < 0.05, "*", "")
    cold_sig <- ifelse(cold_p < 0.05, "*", "")
    cat(sprintf("  %s: Heat p=%.4f%s, Cold p=%.4f%s\n", 
                mod, heat_p, heat_sig, cold_p, cold_sig))
  }
}

# Age Stratification
cat("\n--- Age Stratification ---\n")
for (level in c("intermediate", "immediate")) {
  d <- fromJSON(sprintf("phase4_r/results/age_stratification_%s.json", level))
  cat(sprintf("\n[%s]\n", toupper(level)))
  for (age_grp in names(d$age_groups)) {
    grp <- d$age_groups[[age_grp]]
    cat(sprintf("  %s: Heat RR=%.3f (%.3f-%.3f), Cold RR=%.3f (%.3f-%.3f), Converged=%s\n",
                age_grp, grp$heat$rr, grp$heat$ci_low, grp$heat$ci_high,
                grp$cold$rr, grp$cold$ci_low, grp$cold$ci_high, grp$converged))
  }
}

# Sex Stratification
cat("\n--- Sex Stratification ---\n")
for (level in c("intermediate", "immediate")) {
  d <- fromJSON(sprintf("phase4_r/results/sex_stratification_%s.json", level))
  cat(sprintf("\n[%s]\n", toupper(level)))
  for (sex in names(d$sex_groups)) {
    grp <- d$sex_groups[[sex]]
    cat(sprintf("  %s: Heat RR=%.3f (%.3f-%.3f), Cold RR=%.3f (%.3f-%.3f), Converged=%s\n",
                sex, grp$heat$rr, grp$heat$ci_low, grp$heat$ci_high,
                grp$cold$rr, grp$cold$ci_low, grp$cold$ci_high, grp$converged))
  }
  cat(sprintf("  Sex difference: Heat p=%.5f, Cold p=%.3f\n",
              d$sex_comparison$heat_p, d$sex_comparison$cold_p))
}

# Cause Stratification
cat("\n--- Cause-of-Death Stratification ---\n")
for (level in c("intermediate", "immediate")) {
  d <- fromJSON(sprintf("phase4_r/results/cause_stratification_%s.json", level))
  cat(sprintf("\n[%s]\n", toupper(level)))
  for (cause in names(d$cause_groups)) {
    grp <- d$cause_groups[[cause]]
    cat(sprintf("  %s: Heat RR=%.3f (%.3f-%.3f), Cold RR=%.3f (%.3f-%.3f), Converged=%s\n",
                cause, grp$heat$rr, grp$heat$ci_low, grp$heat$ci_high,
                grp$cold$rr, grp$cold$ci_low, grp$cold$ci_high, grp$converged))
  }
}

cat("\n============================================\n")
cat("EXTRACTION COMPLETE\n")
cat("============================================\n")
