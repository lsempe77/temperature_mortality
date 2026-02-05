################################################################################
# 05b: GENERATE TABLES FOR ALL PHASES (R VERSION)
################################################################################
# Creates publication-ready tables from Phase 1-4 R results
#
# Outputs:
# - Table 1: Main effects summary (both levels)
# - Table 2: Attributable burden summary
# - Table 3: Sensitivity analyses summary
# - Table 4: Stratified analyses (age, sex, cause)
#
# Author: Temperature-Mortality Brazil Analysis Pipeline
# Date: December 2025
################################################################################

library(jsonlite)
library(dplyr)
library(tidyr)
library(openxlsx)

# Directories
script_dir <- tryCatch({
  dirname(rstudioapi::getActiveDocumentContext()$path)
}, error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args[grep("--file=", args)])
  if (length(script_path) > 0) dirname(script_path)
  else "c:/Users/LucasSempe/OneDrive - International Initiative for Impact Evaluation/Desktop/sim_data/new_analysis/phase5_outputs"
})
if (length(script_dir) == 0 || script_dir == "" || script_dir == ".") {
  script_dir <- "c:/Users/LucasSempe/OneDrive - International Initiative for Impact Evaluation/Desktop/sim_data/new_analysis/phase5_outputs"
}
base_dir <- dirname(script_dir)

PHASE1_R <- file.path(base_dir, "phase1_r", "results")
PHASE2_R <- file.path(base_dir, "phase2_r", "results")
PHASE4_R <- file.path(base_dir, "phase4_r", "results")
OUTPUT_DIR <- file.path(script_dir, "tables")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

cat(paste(rep("=", 70), collapse=""), "\n")
cat("05b: GENERATE TABLES (R VERSION)\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Helper functions
load_json <- function(filepath) {
  tryCatch({
    fromJSON(filepath, simplifyVector = FALSE)
  }, error = function(e) {
    cat("  Warning: Could not load", filepath, "\n")
    NULL
  })
}

format_rr <- function(rr, lo, hi, decimals = 3) {
  if (is.null(rr) || is.na(rr)) return("—")
  sprintf("%.*f (%.*f–%.*f)", decimals, rr, decimals, lo, decimals, hi)
}

format_number <- function(n, decimals = 0) {
  if (is.null(n) || is.na(n)) return("—")
  format(round(n, decimals), big.mark = ",", nsmall = decimals)
}

save_table <- function(df, name) {
  # CSV
  write.csv(df, file.path(OUTPUT_DIR, paste0(name, ".csv")), row.names = FALSE)
  # Excel
  write.xlsx(df, file.path(OUTPUT_DIR, paste0(name, ".xlsx")))
  cat("  Saved:", name, "\n")
}

# =============================================================================
# TABLE 1: MAIN EFFECTS SUMMARY
# =============================================================================

generate_table1 <- function() {
  cat("\n[Table 1] Main effects summary...\n")
  
  rows <- list()
  
  for (level in c("intermediate", "immediate")) {
    # DLNM results
    dlnm_file <- file.path(PHASE1_R, paste0("dlnm_r_", level, "_results_v2.json"))
    d <- load_json(dlnm_file)
    
    if (is.null(d)) next
    
    p <- d$pooled
    het <- p$heterogeneity
    
    rows[[level]] <- data.frame(
      Level = tools::toTitleCase(level),
      `N Regions` = p$n_regions,
      `MMT (°C)` = sprintf("%.1f", p$mmt),
      `Heat P99 RR` = format_rr(p$rr_heat_p99$rr, p$rr_heat_p99$rr_lo, p$rr_heat_p99$rr_hi),
      `Heat P95 RR` = format_rr(p$rr_heat_p95$rr, p$rr_heat_p95$rr_lo, p$rr_heat_p95$rr_hi),
      `Cold P1 RR` = format_rr(p$rr_cold_p1$rr, p$rr_cold_p1$rr_lo, p$rr_cold_p1$rr_hi),
      `Cold P5 RR` = format_rr(p$rr_cold_p5$rr, p$rr_cold_p5$rr_lo, p$rr_cold_p5$rr_hi),
      `Cochran Q` = ifelse(is.null(het$cochrans_Q), "—", format_number(het$cochrans_Q, 1)),
      `I² (%)` = ifelse(is.null(het$I_squared), "—", sprintf("%.1f", het$I_squared * 100)),
      check.names = FALSE
    )
  }
  
  df <- bind_rows(rows)
  save_table(df, "table1_main_effects")
  return(df)
}

# =============================================================================
# TABLE 2: ATTRIBUTABLE BURDEN
# =============================================================================

generate_table2 <- function() {
  cat("\n[Table 2] Attributable burden...\n")
  
  rows <- list()
  
  for (level in c("intermediate", "immediate")) {
    burden_file <- file.path(PHASE1_R, paste0("attributable_burden_r_", level, ".json"))
    d <- load_json(burden_file)
    
    if (is.null(d)) next
    
    n <- d$national
    
    rows[[level]] <- data.frame(
      Level = tools::toTitleCase(level),
      `Total Deaths` = format_number(n$total_deaths),
      `Heat AN (P99)` = format_number(n$heat_an_99),
      `Heat AF (%)` = sprintf("%.2f", n$heat_af_pct_99),
      `Cold AN (P1)` = format_number(n$cold_an_1),
      `Cold AF (%)` = sprintf("%.2f", n$cold_af_pct_1),
      `Total AN` = format_number(n$total_an),
      `Total AF (%)` = sprintf("%.2f", n$total_af_pct),
      `Annual Heat` = format_number(n$annual_heat_deaths),
      `Annual Cold` = format_number(n$annual_cold_deaths),
      check.names = FALSE
    )
  }
  
  df <- bind_rows(rows)
  save_table(df, "table2_attributable_burden")
  return(df)
}

# =============================================================================
# TABLE 3: SENSITIVITY ANALYSES
# =============================================================================

generate_table3 <- function() {
  cat("\n[Table 3] Sensitivity analyses...\n")
  
  rows <- list()
  
  for (level in c("intermediate", "immediate")) {
    sens_file <- file.path(PHASE2_R, paste0("sensitivity_r_", level, ".json"))
    d <- load_json(sens_file)
    
    if (is.null(d)) next
    
    for (lag_name in names(d$lag_sensitivity)) {
      lag_data <- d$lag_sensitivity[[lag_name]]
      rows[[paste(level, lag_name)]] <- data.frame(
        Level = tools::toTitleCase(level),
        Analysis = paste("Max Lag", gsub("lag_", "", lag_name), "days"),
        `Heat P99 RR` = format_rr(lag_data$rr_p99, lag_data$rr_p99_ci[[1]], lag_data$rr_p99_ci[[2]]),
        `Cold P1 RR` = format_rr(lag_data$rr_p1, lag_data$rr_p1_ci[[1]], lag_data$rr_p1_ci[[2]]),
        `MMT (°C)` = sprintf("%.1f", lag_data$pooled_mmt),
        `N Regions` = lag_data$n_regions,
        check.names = FALSE
      )
    }
  }
  
  df <- bind_rows(rows)
  save_table(df, "table3_sensitivity")
  return(df)
}

# =============================================================================
# TABLE 4: YLL SUMMARY
# =============================================================================

generate_table4 <- function() {
  cat("\n[Table 4] Years of Life Lost...\n")
  
  rows <- list()
  
  for (level in c("intermediate", "immediate")) {
    yll_file <- file.path(PHASE1_R, paste0("yll_r_", level, ".json"))
    d <- load_json(yll_file)
    
    if (is.null(d)) next
    
    n <- d$national
    
    rows[[level]] <- data.frame(
      Level = tools::toTitleCase(level),
      `Total YLL Heat` = format_number(n$total_yll_heat),
      `Total YLL Cold` = format_number(n$total_yll_cold),
      `Annual YLL Heat` = format_number(n$annual_yll_heat),
      `Annual YLL Cold` = format_number(n$annual_yll_cold),
      `YLL Rate Heat (per 100k)` = sprintf("%.0f", n$yll_rate_heat),
      `YLL Rate Cold (per 100k)` = sprintf("%.0f", n$yll_rate_cold),
      check.names = FALSE
    )
  }
  
  df <- bind_rows(rows)
  save_table(df, "table4_yll")
  return(df)
}

# =============================================================================
# TABLE 5: STRATIFIED ANALYSES
# =============================================================================

generate_table5 <- function() {
  cat("\n[Table 5] Stratified analyses...\n")
  
  rows <- list()
  
  for (level in c("intermediate", "immediate")) {
    # Age stratification
    age_file <- file.path(PHASE4_R, paste0("age_stratification_", level, ".json"))
    d_age <- load_json(age_file)
    
    if (!is.null(d_age)) {
      for (age_group in names(d_age$age_groups)) {
        info <- d_age$age_groups[[age_group]]
        rows[[paste(level, "age", age_group)]] <- data.frame(
          Level = tools::toTitleCase(level),
          Stratum = "Age",
          Group = gsub("age_", "", age_group) %>% gsub("_", "-", .) %>% gsub("plus", "+", .),
          `Heat RR` = format_rr(info$heat$rr, info$heat$ci_low, info$heat$ci_high),
          `Cold RR` = format_rr(info$cold$rr, info$cold$ci_low, info$cold$ci_high),
          check.names = FALSE
        )
      }
    }
    
    # Sex stratification
    sex_file <- file.path(PHASE4_R, paste0("sex_stratification_", level, ".json"))
    d_sex <- load_json(sex_file)
    
    if (!is.null(d_sex)) {
      for (sex in names(d_sex$sex_groups)) {
        info <- d_sex$sex_groups[[sex]]
        rows[[paste(level, "sex", sex)]] <- data.frame(
          Level = tools::toTitleCase(level),
          Stratum = "Sex",
          Group = tools::toTitleCase(sex),
          `Heat RR` = format_rr(info$heat$rr, info$heat$ci_low, info$heat$ci_high),
          `Cold RR` = format_rr(info$cold$rr, info$cold$ci_low, info$cold$ci_high),
          check.names = FALSE
        )
      }
    }
    
    # Cause stratification
    cause_file <- file.path(PHASE4_R, paste0("cause_stratification_", level, ".json"))
    d_cause <- load_json(cause_file)
    
    if (!is.null(d_cause)) {
      for (cause in names(d_cause$cause_groups)) {
        info <- d_cause$cause_groups[[cause]]
        rows[[paste(level, "cause", cause)]] <- data.frame(
          Level = tools::toTitleCase(level),
          Stratum = "Cause",
          Group = gsub("_", " ", cause) %>% tools::toTitleCase(),
          `Heat RR` = format_rr(info$heat$rr, info$heat$ci_low, info$heat$ci_high),
          `Cold RR` = format_rr(info$cold$rr, info$cold$ci_low, info$cold$ci_high),
          check.names = FALSE
        )
      }
    }
  }
  
  df <- bind_rows(rows)
  save_table(df, "table5_stratified")
  return(df)
}

# =============================================================================
# TABLE 6: HARVESTING ANALYSIS
# =============================================================================

generate_table6 <- function() {
  cat("\n[Table 6] Harvesting analysis...\n")
  
  rows <- list()
  
  for (level in c("intermediate", "immediate")) {
    harv_file <- file.path(PHASE2_R, paste0("harvesting_r_", level, ".json"))
    d <- load_json(harv_file)
    
    if (is.null(d)) next
    
    for (lag_name in c("lag_7", "lag_21", "lag_35")) {
      if (!is.null(d$results[[lag_name]])) {
        lag_data <- d$results[[lag_name]]
        heat <- lag_data$heat
        cold <- lag_data$cold
        
        # Handle different field names (rr_lo/rr_hi vs ci_low/ci_high)
        heat_rr_str <- "—"
        if (!is.null(heat) && !is.null(heat$rr) && length(heat$rr) > 0 && !is.na(heat$rr)) {
          lo <- if (!is.null(heat$rr_lo)) heat$rr_lo else heat$ci_low
          hi <- if (!is.null(heat$rr_hi)) heat$rr_hi else heat$ci_high
          if (!is.null(lo) && !is.null(hi)) {
            heat_rr_str <- format_rr(heat$rr, lo, hi)
          }
        }
        
        cold_rr_str <- "—"
        if (!is.null(cold) && !is.null(cold$rr) && length(cold$rr) > 0 && !is.na(cold$rr)) {
          lo <- if (!is.null(cold$rr_lo)) cold$rr_lo else cold$ci_low
          hi <- if (!is.null(cold$rr_hi)) cold$rr_hi else cold$ci_high
          if (!is.null(lo) && !is.null(hi)) {
            cold_rr_str <- format_rr(cold$rr, lo, hi)
          }
        }
        
        rows[[paste(level, lag_name)]] <- data.frame(
          Level = tools::toTitleCase(level),
          `Max Lag` = gsub("lag_", "", lag_name),
          `Heat RR` = heat_rr_str,
          `Cold RR` = cold_rr_str,
          check.names = FALSE
        )
      }
    }
    
    # Add summary if available
    if (!is.null(d$results$summary)) {
      s <- d$results$summary
      heat_ratio <- if (!is.null(s$harvesting_ratio_heat) && !is.na(s$harvesting_ratio_heat)) {
        sprintf("%.3f", s$harvesting_ratio_heat)
      } else "—"
      cold_ratio <- if (!is.null(s$harvesting_ratio_cold) && !is.na(s$harvesting_ratio_cold)) {
        sprintf("%.3f", s$harvesting_ratio_cold)
      } else "—"
      
      rows[[paste(level, "summary")]] <- data.frame(
        Level = tools::toTitleCase(level),
        `Max Lag` = "Ratio",
        `Heat RR` = heat_ratio,
        `Cold RR` = cold_ratio,
        check.names = FALSE
      )
    }
  }
  
  df <- bind_rows(rows)
  save_table(df, "table6_harvesting")
  return(df)
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

cat("\n--- Generating Tables ---\n")
t1 <- generate_table1()
t2 <- generate_table2()
t3 <- generate_table3()
t4 <- generate_table4()
t5 <- generate_table5()
t6 <- generate_table6()

cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("TABLE GENERATION COMPLETE\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat(paste(rep("=", 70), collapse=""), "\n")
