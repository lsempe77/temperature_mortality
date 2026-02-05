#!/usr/bin/env Rscript
# extract_results_v2.R
# Extract key results from all JSON output files

library(jsonlite)

cat("\n========================================\n")
cat("COMPREHENSIVE RESULTS SUMMARY\n")
cat("Temperature-Mortality Analysis (Brazil)\n")
cat("Elderly 60+, 2010-2024\n")
cat("========================================\n\n")

# Get script directory for relative paths
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("--file=", "", file_arg)))
} else {
  script_dir <- getwd()
}

# --- PHASE 1: DLNM Main Results ---
cat("=== PHASE 1: DLNM META-ANALYSIS ===\n\n")

for (level in c("intermediate", "immediate")) {
  results_path <- file.path(script_dir, "phase1_r", "results", 
                            paste0("dlnm_r_", level, "_results_v2.json"))
  
  if (file.exists(results_path)) {
    d <- fromJSON(results_path)
    pr <- d[["pooled"]]
    
    n_reg <- ifelse(level == "intermediate", 133, 510)
    cat(sprintf("%s LEVEL (%d regions)\n", toupper(level), n_reg))
    cat("----------------------------\n")
    
    # Handle both old (rr_lo/rr_hi) and new (ci_low/ci_high) naming
    heat_lo <- if (!is.null(pr$rr_heat_p99$rr_lo)) pr$rr_heat_p99$rr_lo else pr$rr_heat_p99$ci_low
    heat_hi <- if (!is.null(pr$rr_heat_p99$rr_hi)) pr$rr_heat_p99$rr_hi else pr$rr_heat_p99$ci_high
    cold_lo <- if (!is.null(pr$rr_cold_p1$rr_lo)) pr$rr_cold_p1$rr_lo else pr$rr_cold_p1$ci_low
    cold_hi <- if (!is.null(pr$rr_cold_p1$rr_hi)) pr$rr_cold_p1$rr_hi else pr$rr_cold_p1$ci_high
    mmt_val <- if (!is.null(pr$mmt)) pr$mmt else pr$mmt_blup
    
    cat(sprintf("Heat P99: RR = %.3f (%.3f - %.3f)\n", 
                pr$rr_heat_p99$rr, heat_lo, heat_hi))
    cat(sprintf("Cold P1:  RR = %.3f (%.3f - %.3f)\n",
                pr$rr_cold_p1$rr, cold_lo, cold_hi))
    cat(sprintf("MMT: %.1f C\n", mmt_val))
    
    if (!is.null(pr$heterogeneity)) {
      h <- pr$heterogeneity
      if (!is.null(h$cochrans_Q) && !is.na(h$cochrans_Q)) {
        cat(sprintf("Cochran's Q: %.2f (df=%d, p<0.001)\n", h$cochrans_Q, h$Q_df))
      }
      if (!is.null(h$I2_percent) && !is.na(h$I2_percent)) {
        cat(sprintf("I2: %.1f%%\n", h$I2_percent))
      }
    }
    cat("\n")
  } else {
    cat(sprintf("%s: File not found at %s\n\n", toupper(level), results_path))
  }
}

# --- PHASE 1: Attributable Burden ---
cat("=== ATTRIBUTABLE BURDEN ===\n\n")

for (level in c("intermediate", "immediate")) {
  results_path <- file.path(script_dir, "phase1_r", "results",
                            paste0("attributable_burden_r_", level, ".json"))
  
  if (file.exists(results_path)) {
    d <- fromJSON(results_path)
    n <- d[["national"]]
    
    cat(sprintf("%s:\n", toupper(level)))
    cat(sprintf("  Total deaths analyzed: %s\n", format(round(n$total_deaths), big.mark=",")))
    cat(sprintf("  Total attributable: %s (%.1f%% AF)\n", 
                format(round(n$total_an), big.mark=","), n$total_af_pct))
    cat(sprintf("  Heat attributable (P99): %s (%.2f%% AF)\n", 
                format(round(n$heat_an_99), big.mark=","), n$heat_af_pct_99))
    cat(sprintf("  Cold attributable (P1): %s (%.2f%% AF)\n",
                format(round(n$cold_an_1), big.mark=","), n$cold_af_pct_1))
    cat(sprintf("  Annual heat deaths: %s\n", format(round(n$annual_heat_deaths), big.mark=",")))
    cat(sprintf("  Annual cold deaths: %s\n", format(round(n$annual_cold_deaths), big.mark=",")))
    cat("\n")
  }
}

# --- PHASE 1: YLL ---
cat("=== YEARS OF LIFE LOST ===\n\n")

for (level in c("intermediate", "immediate")) {
  results_path <- file.path(script_dir, "phase1_r", "results",
                            paste0("yll_r_", level, ".json"))
  
  if (file.exists(results_path)) {
    d <- fromJSON(results_path)
    ns <- d[["national"]]
    
    cat(sprintf("%s:\n", toupper(level)))
    cat(sprintf("  Annual YLL heat: %s (%.0f per 100k elderly)\n",
                format(round(ns$annual_yll_heat), big.mark=","), ns$yll_rate_heat_per_100k))
    cat(sprintf("  Annual YLL cold: %s (%.0f per 100k elderly)\n",
                format(round(ns$annual_yll_cold), big.mark=","), ns$yll_rate_cold_per_100k))
    cat(sprintf("  YLL per attributable death: %.1f years\n", ns$yll_per_death))
    cat("\n")
  }
}

# --- PHASE 2: Robustness Checks ---
cat("=== PHASE 2: ROBUSTNESS ===\n\n")

# Sensitivity
for (level in c("intermediate", "immediate")) {
  results_path <- file.path(script_dir, "phase2_r", "results",
                            paste0("sensitivity_r_", level, ".json"))
  
  if (file.exists(results_path)) {
    d <- fromJSON(results_path, simplifyVector = FALSE)
    cat(sprintf("%s Sensitivity:\n", toupper(level)))
    
    if (!is.null(d$lag_sensitivity)) {
      for (lag_name in names(d$lag_sensitivity)) {
        ls <- d$lag_sensitivity[[lag_name]]
        if (!is.null(ls$rr_p99) && !is.null(ls$rr_p1)) {
          cat(sprintf("  %s: Heat P99=%.3f, Cold P1=%.3f, MMT=%.1f\n",
                      lag_name, ls$rr_p99, ls$rr_p1, ls$pooled_mmt))
        }
      }
    }
    cat("\n")
  }
}

# Harvesting
cat("--- Harvesting (Lag Extension) ---\n\n")
for (level in c("intermediate", "immediate")) {
  results_path <- file.path(script_dir, "phase2_r", "results",
                            paste0("harvesting_r_", level, ".json"))
  
  if (file.exists(results_path)) {
    d <- fromJSON(results_path, simplifyVector = FALSE)
    cat(sprintf("%s:\n", toupper(level)))
    
    if (!is.null(d$results)) {
      for (lag_name in c("lag_7", "lag_21", "lag_35")) {
        if (!is.null(d$results[[lag_name]])) {
          r <- d$results[[lag_name]]
          heat_rr <- if (!is.null(r$heat$rr)) r$heat$rr else NA
          cold_rr <- if (!is.null(r$cold$rr)) r$cold$rr else NA
          cat(sprintf("  %s: Heat=%.3f, Cold=%.3f\n", lag_name, heat_rr, cold_rr))
        }
      }
      if (!is.null(d$results$summary)) {
        cat(sprintf("  Harvesting ratio (heat): %.3f\n", d$results$summary$harvesting_ratio_heat))
        cat(sprintf("  Harvesting ratio (cold): %.3f\n", d$results$summary$harvesting_ratio_cold))
      }
    }
    cat("\n")
  }
}

# --- PHASE 4: Stratification ---
cat("=== PHASE 4: STRATIFICATION ===\n\n")

# Age
cat("--- Age Groups ---\n")
for (level in c("intermediate", "immediate")) {
  results_path <- file.path(script_dir, "phase4_r", "results",
                            paste0("age_stratification_", level, ".json"))
  
  if (file.exists(results_path)) {
    d <- fromJSON(results_path, simplifyVector = FALSE)
    cat(sprintf("%s:\n", toupper(level)))
    
    if (!is.null(d$age_groups)) {
      for (age in names(d$age_groups)) {
        r <- d$age_groups[[age]]
        if (!is.null(r$heat) && !is.null(r$heat$rr)) {
          cat(sprintf("  %s: Heat=%.3f (%.3f-%.3f), Cold=%.3f (%.3f-%.3f)\n",
                      gsub("age_", "", gsub("_", "-", age)), 
                      r$heat$rr, r$heat$ci_low, r$heat$ci_high,
                      r$cold$rr, r$cold$ci_low, r$cold$ci_high))
        }
      }
    }
    cat("\n")
  }
}

# Sex
cat("--- Sex ---\n")
for (level in c("intermediate", "immediate")) {
  results_path <- file.path(script_dir, "phase4_r", "results",
                            paste0("sex_stratification_", level, ".json"))
  
  if (file.exists(results_path)) {
    d <- fromJSON(results_path, simplifyVector = FALSE)
    cat(sprintf("%s:\n", toupper(level)))
    
    if (!is.null(d$sex_groups)) {
      for (sex in names(d$sex_groups)) {
        r <- d$sex_groups[[sex]]
        if (!is.null(r$heat) && !is.null(r$heat$rr)) {
          cat(sprintf("  %s: Heat=%.3f (%.3f-%.3f), Cold=%.3f (%.3f-%.3f)\n",
                      sex, r$heat$rr, r$heat$ci_low, r$heat$ci_high,
                      r$cold$rr, r$cold$ci_low, r$cold$ci_high))
        }
      }
    }
    cat("\n")
  }
}

cat("========================================\n")
cat("Results extracted:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
