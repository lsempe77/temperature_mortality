################################################################################
# 05a: GENERATE FIGURES FOR ALL PHASES (R VERSION)
################################################################################
# Creates publication-quality figures from Phase 1-4 R results
#
# Outputs:
# - Phase 1: Exposure-response curves, attributable burden, YLL
# - Phase 2: Sensitivity forest plots, harvesting, heatwave effects  
# - Phase 4: Age/sex/cause stratification
#
# Author: Temperature-Mortality Brazil Analysis Pipeline
# Date: December 2025
################################################################################

# =============================================================================
# SETUP
# =============================================================================

library(ggplot2)
library(jsonlite)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# Directories - handle both RStudio and command-line execution
script_dir <- tryCatch({
  dirname(rstudioapi::getActiveDocumentContext()$path)
}, error = function(e) {
  # Running from command line
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args[grep("--file=", args)])
  if (length(script_path) > 0) {
    dirname(script_path)
  } else {
    "c:/Users/LucasSempe/OneDrive - International Initiative for Impact Evaluation/Desktop/sim_data/new_analysis/phase5_outputs"
  }
})

if (length(script_dir) == 0 || script_dir == "" || script_dir == ".") {
  script_dir <- "c:/Users/LucasSempe/OneDrive - International Initiative for Impact Evaluation/Desktop/sim_data/new_analysis/phase5_outputs"
}
base_dir <- dirname(script_dir)

PHASE1_R <- file.path(base_dir, "phase1_r", "results")
PHASE2_R <- file.path(base_dir, "phase2_r", "results")
PHASE4_R <- file.path(base_dir, "phase4_r", "results")
OUTPUT_DIR <- file.path(script_dir, "figures")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Theme for publication
theme_publication <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "sans", size = 11),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 11, face = "bold")
    )
}

# Colors
HEAT_COLOR <- "#e74c3c"
COLD_COLOR <- "#3498db"
NEUTRAL_COLOR <- "#7f8c8d"

cat(paste(rep("=", 70), collapse=""), "\n")
cat("05a: GENERATE FIGURES (R VERSION)\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

load_json <- function(filepath) {
  tryCatch({
    fromJSON(filepath, simplifyVector = FALSE)
  }, error = function(e) {
    cat("  Warning: Could not load", filepath, "\n")
    NULL
  })
}

save_figure <- function(p, name, width = 10, height = 6) {
  filepath <- file.path(OUTPUT_DIR, paste0(name, ".png"))
  ggsave(filepath, p, width = width, height = height, dpi = 300)
  cat("  Saved:", name, "\n")
}

# =============================================================================
# FIGURE 1: POOLED EXPOSURE-RESPONSE CURVES
# =============================================================================

plot_exposure_response <- function() {
  cat("\n[Phase 1] Pooled exposure-response curves...\n")
  
  plots <- list()
  
  for (level in c("intermediate", "immediate")) {
    results_file <- file.path(PHASE1_R, paste0("dlnm_r_", level, "_results_v2.json"))
    d <- load_json(results_file)
    
    if (is.null(d)) next
    
    # Extract curve data
    curve <- d$pooled$temp_curve
    temps <- unlist(curve$temp)
    rrs <- unlist(curve$rr)
    rr_los <- unlist(curve$rr_lo)
    rr_his <- unlist(curve$rr_hi)
    mmt <- d$pooled$mmt
    
    df <- data.frame(temp = temps, rr = rrs, rr_lo = rr_los, rr_hi = rr_his)
    
    # Get key RRs for annotation
    heat_rr <- d$pooled$rr_heat_p99$rr
    heat_lo <- d$pooled$rr_heat_p99$rr_lo
    heat_hi <- d$pooled$rr_heat_p99$rr_hi
    cold_rr <- d$pooled$rr_cold_p1$rr
    cold_lo <- d$pooled$rr_cold_p1$rr_lo
    cold_hi <- d$pooled$rr_cold_p1$rr_hi
    
    level_label <- ifelse(level == "intermediate", 
                          "Intermediate (133 regions)", 
                          "Immediate (510 regions)")
    
    p <- ggplot(df, aes(x = temp)) +
      geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), fill = "gray70", alpha = 0.4) +
      geom_line(aes(y = rr), color = "black", linewidth = 1) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = mmt, linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
      annotate("text", x = mmt + 1, y = max(rrs) * 0.95, 
               label = sprintf("MMT = %.1f°C", mmt), 
               hjust = 0, size = 3.5, color = "darkgreen") +
      annotate("label", x = max(temps) - 2, y = 1.5, 
               label = sprintf("Heat P99: %.3f (%.3f-%.3f)\nCold P1: %.3f (%.3f-%.3f)", 
                               heat_rr, heat_lo, heat_hi, cold_rr, cold_lo, cold_hi),
               hjust = 1, size = 3, fill = "white", alpha = 0.8) +
      scale_x_continuous(breaks = seq(5, 35, 5)) +
      coord_cartesian(ylim = c(0.8, min(max(rrs) * 1.1, 4))) +
      labs(x = "Temperature (°C)", y = "Relative Risk",
           title = level_label) +
      theme_publication()
    
    plots[[level]] <- p
  }
  
  # Combine plots
  combined <- plots$intermediate + plots$immediate +
    plot_annotation(
      title = "Temperature-Mortality Exposure-Response Curves",
      subtitle = "Pooled across regions using multivariate meta-analysis (21-day cumulative lag)",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 11, hjust = 0.5))
    )
  
  save_figure(combined, "fig1_exposure_response", width = 12, height = 5)
}

# =============================================================================
# FIGURE 2: ATTRIBUTABLE BURDEN
# =============================================================================

plot_attributable_burden <- function() {
  cat("\n[Phase 1] Attributable burden...\n")
  
  burden_data <- list()
  
  for (level in c("intermediate", "immediate")) {
    results_file <- file.path(PHASE1_R, paste0("attributable_burden_r_", level, ".json"))
    d <- load_json(results_file)
    
    if (is.null(d)) next
    
    burden_data[[level]] <- data.frame(
      level = tools::toTitleCase(level),
      heat_an = d$national$heat_an_99,
      cold_an = d$national$cold_an_1,
      total_af = d$national$total_af_pct,
      annual_heat = d$national$annual_heat_deaths,
      annual_cold = d$national$annual_cold_deaths
    )
  }
  
  df <- bind_rows(burden_data)
  
  # Plot 1: Total attributable deaths by type
  df_long <- df %>%
    select(level, heat_an, cold_an) %>%
    pivot_longer(cols = c(heat_an, cold_an), names_to = "type", values_to = "deaths") %>%
    mutate(type = ifelse(type == "heat_an", "Heat (P99)", "Cold (P1)"))
  
  p1 <- ggplot(df_long, aes(x = level, y = deaths / 1000, fill = type)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
    geom_text(aes(label = sprintf("%.0fk", deaths/1000)), 
              position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("Heat (P99)" = HEAT_COLOR, "Cold (P1)" = COLD_COLOR)) +
    labs(x = "", y = "Attributable Deaths (thousands)", 
         title = "Total Attributable Deaths (2010-2024)",
         fill = "Temperature Extreme") +
    theme_publication() +
    theme(legend.position = "bottom")
  
  # Plot 2: Annual burden
  df_annual <- df %>%
    select(level, annual_heat, annual_cold) %>%
    pivot_longer(cols = c(annual_heat, annual_cold), names_to = "type", values_to = "deaths") %>%
    mutate(type = ifelse(type == "annual_heat", "Heat", "Cold"))
  
  p2 <- ggplot(df_annual, aes(x = level, y = deaths, fill = type)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
    geom_text(aes(label = sprintf("%.0f", deaths)), 
              position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("Heat" = HEAT_COLOR, "Cold" = COLD_COLOR)) +
    scale_y_continuous(labels = comma) +
    labs(x = "", y = "Annual Attributable Deaths", 
         title = "Average Annual Burden",
         fill = "Type") +
    theme_publication() +
    theme(legend.position = "bottom")
  
  combined <- p1 + p2 +
    plot_annotation(
      title = "Temperature-Attributable Mortality Burden in Brazilian Elderly",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  save_figure(combined, "fig2_attributable_burden", width = 12, height = 5)
}

# =============================================================================
# FIGURE 3: YLL SUMMARY
# =============================================================================

plot_yll_summary <- function() {
  cat("\n[Phase 1] Years of life lost...\n")
  
  yll_data <- list()
  
  for (level in c("intermediate", "immediate")) {
    results_file <- file.path(PHASE1_R, paste0("yll_r_", level, ".json"))
    d <- load_json(results_file)
    
    if (is.null(d)) next
    
    yll_data[[level]] <- data.frame(
      level = tools::toTitleCase(level),
      heat_yll = d$national$annual_yll_heat,
      cold_yll = d$national$annual_yll_cold,
      heat_rate = d$national$yll_rate_heat,
      cold_rate = d$national$yll_rate_cold
    )
  }
  
  df <- bind_rows(yll_data)
  
  # Plot: Annual YLL by type
  df_long <- df %>%
    select(level, heat_yll, cold_yll) %>%
    pivot_longer(cols = c(heat_yll, cold_yll), names_to = "type", values_to = "yll") %>%
    mutate(type = ifelse(type == "heat_yll", "Heat", "Cold"))
  
  p1 <- ggplot(df_long, aes(x = level, y = yll / 1000, fill = type)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
    geom_text(aes(label = sprintf("%.0fk", yll/1000)), 
              position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("Heat" = HEAT_COLOR, "Cold" = COLD_COLOR)) +
    labs(x = "", y = "Annual YLL (thousands)", 
         title = "Annual Years of Life Lost",
         fill = "Temperature") +
    theme_publication() +
    theme(legend.position = "bottom")
  
  # Plot 2: YLL rate per 100k
  df_rate <- df %>%
    select(level, heat_rate, cold_rate) %>%
    pivot_longer(cols = c(heat_rate, cold_rate), names_to = "type", values_to = "rate") %>%
    mutate(type = ifelse(type == "heat_rate", "Heat", "Cold"))
  
  p2 <- ggplot(df_rate, aes(x = level, y = rate, fill = type)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
    geom_text(aes(label = sprintf("%.0f", rate)), 
              position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("Heat" = HEAT_COLOR, "Cold" = COLD_COLOR)) +
    labs(x = "", y = "YLL Rate per 100,000 elderly", 
         title = "YLL Rate by Temperature Type",
         fill = "Type") +
    theme_publication() +
    theme(legend.position = "bottom")
  
  combined <- p1 + p2 +
    plot_annotation(
      title = "Years of Life Lost Attributable to Temperature",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  save_figure(combined, "fig3_yll_summary", width = 12, height = 5)
}

# =============================================================================
# FIGURE 4: SENSITIVITY ANALYSIS (LAG STRUCTURE)
# =============================================================================

plot_sensitivity_lag <- function() {
  cat("\n[Phase 2] Sensitivity analysis - lag structure...\n")
  
  sens_data <- list()
  
  for (level in c("intermediate", "immediate")) {
    results_file <- file.path(PHASE2_R, paste0("sensitivity_r_", level, ".json"))
    d <- load_json(results_file)
    
    if (is.null(d)) next
    
    lag_sens <- d$lag_sensitivity
    for (lag_name in names(lag_sens)) {
      lag_data <- lag_sens[[lag_name]]
      sens_data[[paste(level, lag_name)]] <- data.frame(
        level = tools::toTitleCase(level),
        lag = gsub("lag_", "", lag_name),
        heat_rr = lag_data$rr_p99,
        heat_lo = lag_data$rr_p99_ci[[1]],
        heat_hi = lag_data$rr_p99_ci[[2]],
        cold_rr = lag_data$rr_p1,
        cold_lo = lag_data$rr_p1_ci[[1]],
        cold_hi = lag_data$rr_p1_ci[[2]]
      )
    }
  }
  
  df <- bind_rows(sens_data) %>%
    mutate(lag = factor(lag, levels = c("7", "14", "21", "28")))
  
  # Heat effects
  p1 <- ggplot(df, aes(x = lag, y = heat_rr, color = level, group = level)) +
    geom_point(size = 3, position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = heat_lo, ymax = heat_hi), width = 0.2, 
                  position = position_dodge(width = 0.3)) +
    geom_line(position = position_dodge(width = 0.3), linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Intermediate" = HEAT_COLOR, "Immediate" = "#c0392b")) +
    labs(x = "Maximum Lag (days)", y = "Relative Risk (P99 Heat)",
         title = "Heat Effects by Lag Duration", color = "Level") +
    theme_publication()
  
  # Cold effects
  p2 <- ggplot(df, aes(x = lag, y = cold_rr, color = level, group = level)) +
    geom_point(size = 3, position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = cold_lo, ymax = cold_hi), width = 0.2, 
                  position = position_dodge(width = 0.3)) +
    geom_line(position = position_dodge(width = 0.3), linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Intermediate" = COLD_COLOR, "Immediate" = "#2980b9")) +
    labs(x = "Maximum Lag (days)", y = "Relative Risk (P1 Cold)",
         title = "Cold Effects by Lag Duration", color = "Level") +
    theme_publication()
  
  combined <- p1 + p2 +
    plot_annotation(
      title = "Sensitivity Analysis: Effect of Maximum Lag Specification",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  save_figure(combined, "fig4_sensitivity_lag", width = 12, height = 5)
}

# =============================================================================
# FIGURE 5: HARVESTING ANALYSIS (Updated for v2 methodology)
# =============================================================================

plot_harvesting <- function() {
  cat("\n[Phase 2] Harvesting analysis (v2)...\n")
  
  harv_data <- list()
  
  for (level in c("intermediate")) {
    # Use v2 results file (corrected single-model methodology)
    results_file <- file.path(PHASE2_R, paste0("harvesting_v2_", level, ".json"))
    d <- load_json(results_file)
    
    if (is.null(d)) {
      # Fallback to old file if v2 not available
      results_file <- file.path(PHASE2_R, paste0("harvesting_r_", level, ".json"))
      d <- load_json(results_file)
    }
    
    if (is.null(d)) next
    
    # v2 structure: results$lag_X$heat$ci_low, ci_high (not rr_lo, rr_hi)
    for (lag_name in c("lag_7", "lag_14", "lag_21", "lag_28", "lag_35")) {
      if (!is.null(d$results[[lag_name]])) {
        lag_data <- d$results[[lag_name]]
        harv_data[[paste(level, lag_name)]] <- data.frame(
          level = tools::toTitleCase(level),
          lag = as.numeric(gsub("lag_", "", lag_name)),
          heat_rr = ifelse(is.null(lag_data$heat$rr), NA, lag_data$heat$rr),
          heat_lo = ifelse(is.null(lag_data$heat$ci_low), 
                          ifelse(is.null(lag_data$heat$rr_lo), NA, lag_data$heat$rr_lo),
                          lag_data$heat$ci_low),
          heat_hi = ifelse(is.null(lag_data$heat$ci_high), 
                          ifelse(is.null(lag_data$heat$rr_hi), NA, lag_data$heat$rr_hi),
                          lag_data$heat$ci_high),
          cold_rr = ifelse(is.null(lag_data$cold$rr), NA, lag_data$cold$rr),
          cold_lo = ifelse(is.null(lag_data$cold$ci_low), 
                          ifelse(is.null(lag_data$cold$rr_lo), NA, lag_data$cold$rr_lo),
                          lag_data$cold$ci_low),
          cold_hi = ifelse(is.null(lag_data$cold$ci_high), 
                          ifelse(is.null(lag_data$cold$rr_hi), NA, lag_data$cold$rr_hi),
                          lag_data$cold$ci_high),
          i2 = ifelse(is.null(lag_data$heat$I2_percent), NA, lag_data$heat$I2_percent)
        )
      }
    }
  }
  
  df <- bind_rows(harv_data) %>%
    filter(!is.na(heat_rr)) %>%
    mutate(lag = factor(lag, levels = c(7, 14, 21, 28, 35)))
  
  if (nrow(df) == 0) {
    cat("  No harvesting data available\n")
    return(invisible(NULL))
  }
  
  # Heat effects across lag horizons (line plot for v2)
  p1 <- ggplot(df, aes(x = lag, y = heat_rr, group = level)) +
    geom_point(size = 3, color = HEAT_COLOR) +
    geom_errorbar(aes(ymin = heat_lo, ymax = heat_hi), 
                  width = 0.2, color = HEAT_COLOR) +
    geom_line(color = HEAT_COLOR, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    labs(x = "Cumulative Lag Horizon (days)", y = "Relative Risk (Heat P99)",
         title = "Heat Effects Show Attenuation",
         subtitle = "Pattern consistent with mortality displacement") +
    theme_publication() +
    theme(plot.subtitle = element_text(size = 10, color = "gray40"))
  
  # Cold effects - only show if reasonable (filter extreme values)
  df_cold <- df %>% filter(cold_rr > 0.01 & cold_rr < 100 & cold_hi < 10)
  
  if (nrow(df_cold) > 0) {
    p2 <- ggplot(df_cold, aes(x = lag, y = cold_rr, group = level)) +
      geom_point(size = 3, color = COLD_COLOR) +
      geom_errorbar(aes(ymin = cold_lo, ymax = cold_hi), 
                    width = 0.2, color = COLD_COLOR) +
      geom_line(color = COLD_COLOR, linetype = "dashed") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      labs(x = "Cumulative Lag Horizon (days)", y = "Relative Risk (Cold P1)",
           title = "Cold Effects Unstable",
           subtitle = "High heterogeneity precludes interpretation") +
      theme_publication() +
      theme(plot.subtitle = element_text(size = 10, color = "gray40"))
  } else {
    # If cold results too unstable, show placeholder
    p2 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Cold harvesting estimates\ntoo unstable to display\n(I² > 96%)",
               size = 4, color = "gray50") +
      labs(title = "Cold Effects Not Interpretable") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  }
  
  combined <- p1 + p2 +
    plot_annotation(
      title = "Harvesting Analysis: Cumulative Effects at Extended Lags",
      subtitle = "Single-model methodology (v2): fit once, extract cumulative effects at each horizon",
      caption = "Note: Extreme heterogeneity (I² = 98-99%) indicates regional variation; pooled estimates are exploratory only",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 11, hjust = 0.5),
                    plot.caption = element_text(size = 9, color = "gray40", hjust = 0.5))
    )
  
  save_figure(combined, "fig5_harvesting", width = 12, height = 5)
}

# =============================================================================
# FIGURE 6: HEATWAVE EFFECTS
# =============================================================================

plot_heatwave <- function() {
  cat("\n[Phase 2] Heatwave analysis...\n")
  
  hw_data <- list()
  
  for (level in c("intermediate", "immediate")) {
    results_file <- file.path(PHASE2_R, paste0("heatwave_r_", level, ".json"))
    d <- load_json(results_file)
    
    if (is.null(d)) next
    
    # Structure is under 'results' not 'pooled_results'
    if (!is.null(d$results)) {
      hw_data[[level]] <- data.frame(
        level = tools::toTitleCase(level),
        hw_rr = d$results$heatwave_rr,
        hw_lo = d$results$heatwave_ci_low,
        hw_hi = d$results$heatwave_ci_high,
        main_rr = d$results$avg_main_rr_p99,  # Average main effect for comparison
        n_hw_days = d$results$n_heatwave_days,
        pct_hw = d$results$pct_heatwave
      )
    }
  }
  
  if (length(hw_data) == 0) {
    cat("  No heatwave data available\n")
    return(invisible(NULL))
  }
  
  df <- bind_rows(hw_data)
  
  # Compare heatwave additive effect vs main effect
  df_long <- df %>%
    select(level, hw_rr, main_rr) %>%
    pivot_longer(cols = c(hw_rr, main_rr), names_to = "type", values_to = "rr") %>%
    mutate(type = ifelse(type == "hw_rr", "Heatwave Effect\n(Additive)", "Main Heat Effect\n(P99)"))
  
  p <- ggplot(df_long, aes(x = level, y = rr, fill = type)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("Heatwave Effect\n(Additive)" = HEAT_COLOR, 
                                 "Main Heat Effect\n(P99)" = NEUTRAL_COLOR)) +
    labs(x = "", y = "Relative Risk",
         title = "Heatwave Effect Modification",
         subtitle = sprintf("Heatwave days: Intermediate=%.1f%%, Immediate=%.1f%%", 
                           df$pct_hw[1], df$pct_hw[2]),
         fill = "") +
    theme_publication() +
    theme(legend.position = "bottom")
  
  save_figure(p, "fig6_heatwave", width = 8, height = 6)
}

# =============================================================================
# FIGURE 7: AGE STRATIFICATION
# =============================================================================

plot_age_stratification <- function() {
  cat("\n[Phase 4] Age stratification...\n")
  
  age_data <- list()
  
  for (level in c("intermediate", "immediate")) {
    results_file <- file.path(PHASE4_R, paste0("age_stratification_", level, ".json"))
    d <- load_json(results_file)
    
    if (is.null(d)) next
    
    for (age_group in names(d$age_groups)) {
      age_info <- d$age_groups[[age_group]]
      age_data[[paste(level, age_group)]] <- data.frame(
        level = tools::toTitleCase(level),
        age_group = gsub("age_", "", age_group) %>% gsub("_", "-", .) %>% gsub("plus", "+", .),
        heat_rr = age_info$heat$rr,
        heat_lo = age_info$heat$ci_low,
        heat_hi = age_info$heat$ci_high,
        cold_rr = age_info$cold$rr,
        cold_lo = age_info$cold$ci_low,
        cold_hi = age_info$cold$ci_high
      )
    }
  }
  
  df <- bind_rows(age_data) %>%
    mutate(age_group = factor(age_group, levels = c("60-69", "70-79", "80+")))
  
  # Heat effects by age
  p1 <- ggplot(df, aes(x = age_group, y = heat_rr, fill = level)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = heat_lo, ymax = heat_hi), 
                  position = position_dodge(width = 0.7), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("Intermediate" = HEAT_COLOR, "Immediate" = "#c0392b")) +
    labs(x = "Age Group", y = "Relative Risk",
         title = "Heat Effects (P99)", fill = "Level") +
    theme_publication()
  
  # Cold effects by age
  p2 <- ggplot(df, aes(x = age_group, y = cold_rr, fill = level)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = cold_lo, ymax = cold_hi), 
                  position = position_dodge(width = 0.7), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("Intermediate" = COLD_COLOR, "Immediate" = "#2980b9")) +
    labs(x = "Age Group", y = "Relative Risk",
         title = "Cold Effects (P1)", fill = "Level") +
    theme_publication()
  
  combined <- p1 + p2 +
    plot_annotation(
      title = "Age-Stratified Temperature-Mortality Effects",
      subtitle = "Vulnerability increases with age for both heat and cold",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 11, hjust = 0.5))
    )
  
  save_figure(combined, "fig7_age_stratification", width = 12, height = 5)
}

# =============================================================================
# FIGURE 8: SEX STRATIFICATION
# =============================================================================

plot_sex_stratification <- function() {
  cat("\n[Phase 4] Sex stratification...\n")
  
  sex_data <- list()
  
  for (level in c("intermediate", "immediate")) {
    results_file <- file.path(PHASE4_R, paste0("sex_stratification_", level, ".json"))
    d <- load_json(results_file)
    
    if (is.null(d)) next
    
    for (sex in names(d$sex_groups)) {
      sex_info <- d$sex_groups[[sex]]
      sex_data[[paste(level, sex)]] <- data.frame(
        level = tools::toTitleCase(level),
        sex = tools::toTitleCase(sex),
        heat_rr = sex_info$heat$rr,
        heat_lo = sex_info$heat$ci_low,
        heat_hi = sex_info$heat$ci_high,
        cold_rr = sex_info$cold$rr,
        cold_lo = sex_info$cold$ci_low,
        cold_hi = sex_info$cold$ci_high
      )
    }
  }
  
  df <- bind_rows(sex_data)
  
  # Combined heat and cold
  df_long <- df %>%
    pivot_longer(cols = c(heat_rr, cold_rr), names_to = "temp_type", values_to = "rr") %>%
    mutate(
      temp_type = ifelse(temp_type == "heat_rr", "Heat (P99)", "Cold (P1)"),
      rr_lo = ifelse(temp_type == "Heat (P99)", heat_lo, cold_lo),
      rr_hi = ifelse(temp_type == "Heat (P99)", heat_hi, cold_hi)
    )
  
  p <- ggplot(df_long, aes(x = sex, y = rr, fill = temp_type)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = rr_lo, ymax = rr_hi), 
                  position = position_dodge(width = 0.7), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("Heat (P99)" = HEAT_COLOR, "Cold (P1)" = COLD_COLOR)) +
    facet_wrap(~level) +
    labs(x = "Sex", y = "Relative Risk",
         title = "Sex-Stratified Temperature-Mortality Effects",
         subtitle = "Females show higher heat vulnerability; males show higher cold vulnerability",
         fill = "Temperature") +
    theme_publication() +
    theme(legend.position = "bottom")
  
  save_figure(p, "fig8_sex_stratification", width = 10, height = 6)
}

# =============================================================================
# FIGURE 9: CAUSE STRATIFICATION
# =============================================================================

plot_cause_stratification <- function() {
  cat("\n[Phase 4] Cause stratification...\n")
  
  cause_data <- list()
  
  for (level in c("intermediate", "immediate")) {
    results_file <- file.path(PHASE4_R, paste0("cause_stratification_", level, ".json"))
    d <- load_json(results_file)
    
    if (is.null(d)) next
    
    for (cause in names(d$cause_groups)) {
      cause_info <- d$cause_groups[[cause]]
      cause_data[[paste(level, cause)]] <- data.frame(
        level = tools::toTitleCase(level),
        cause = gsub("_", " ", cause) %>% tools::toTitleCase(),
        heat_rr = cause_info$heat$rr,
        heat_lo = cause_info$heat$ci_low,
        heat_hi = cause_info$heat$ci_high,
        cold_rr = cause_info$cold$rr,
        cold_lo = cause_info$cold$ci_low,
        cold_hi = cause_info$cold$ci_high
      )
    }
  }
  
  df <- bind_rows(cause_data) %>%
    mutate(cause = factor(cause, levels = c("Cardiovascular", "Respiratory", "External", "Other")))
  
  # Heat effects by cause
  p1 <- ggplot(df, aes(x = cause, y = heat_rr, fill = level)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = heat_lo, ymax = heat_hi), 
                  position = position_dodge(width = 0.7), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("Intermediate" = HEAT_COLOR, "Immediate" = "#c0392b")) +
    labs(x = "Cause of Death", y = "Relative Risk",
         title = "Heat Effects (P99)", fill = "Level") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Cold effects by cause
  p2 <- ggplot(df, aes(x = cause, y = cold_rr, fill = level)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = cold_lo, ymax = cold_hi), 
                  position = position_dodge(width = 0.7), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("Intermediate" = COLD_COLOR, "Immediate" = "#2980b9")) +
    labs(x = "Cause of Death", y = "Relative Risk",
         title = "Cold Effects (P1)", fill = "Level") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  combined <- p1 + p2 +
    plot_annotation(
      title = "Cause-Specific Temperature-Mortality Effects",
      subtitle = "Cardiovascular and respiratory causes show strongest associations",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 11, hjust = 0.5))
    )
  
  save_figure(combined, "fig9_cause_stratification", width = 12, height = 6)
}

# =============================================================================
# FIGURE 10: SUMMARY FOREST PLOT
# =============================================================================

plot_summary_forest <- function() {
  cat("\n[Summary] Main effects forest plot...\n")
  
  # Collect main results
  results <- list()
  
  for (level in c("intermediate", "immediate")) {
    # DLNM main effects
    dlnm_file <- file.path(PHASE1_R, paste0("dlnm_r_", level, "_results_v2.json"))
    d <- load_json(dlnm_file)
    
    if (!is.null(d)) {
      results[[paste(level, "main")]] <- data.frame(
        level = tools::toTitleCase(level),
        analysis = "Main Model (21-day lag)",
        heat_rr = d$pooled$rr_heat_p99$rr,
        heat_lo = d$pooled$rr_heat_p99$rr_lo,
        heat_hi = d$pooled$rr_heat_p99$rr_hi,
        cold_rr = d$pooled$rr_cold_p1$rr,
        cold_lo = d$pooled$rr_cold_p1$rr_lo,
        cold_hi = d$pooled$rr_cold_p1$rr_hi
      )
    }
    
    # Sensitivity - different lags
    sens_file <- file.path(PHASE2_R, paste0("sensitivity_r_", level, ".json"))
    s <- load_json(sens_file)
    
    if (!is.null(s)) {
      for (lag in c("7", "14", "28")) {
        lag_key <- paste0("lag_", lag)
        if (!is.null(s$lag_sensitivity[[lag_key]])) {
          lag_data <- s$lag_sensitivity[[lag_key]]
          results[[paste(level, "lag", lag)]] <- data.frame(
            level = tools::toTitleCase(level),
            analysis = paste0("Lag ", lag, " days"),
            heat_rr = lag_data$rr_p99,
            heat_lo = lag_data$rr_p99_ci[[1]],
            heat_hi = lag_data$rr_p99_ci[[2]],
            cold_rr = lag_data$rr_p1,
            cold_lo = lag_data$rr_p1_ci[[1]],
            cold_hi = lag_data$rr_p1_ci[[2]]
          )
        }
      }
    }
  }
  
  df <- bind_rows(results) %>%
    mutate(analysis = factor(analysis, levels = rev(c("Main Model (21-day lag)", 
                                                       "Lag 7 days", "Lag 14 days", "Lag 28 days"))))
  
  # Forest plot - Heat
  p1 <- ggplot(df, aes(x = heat_rr, y = analysis, color = level)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbarh(aes(xmin = heat_lo, xmax = heat_hi), 
                   position = position_dodge(width = 0.5), height = 0.2) +
    scale_color_manual(values = c("Intermediate" = HEAT_COLOR, "Immediate" = "#c0392b")) +
    labs(x = "Relative Risk", y = "", title = "Heat (P99)", color = "Level") +
    theme_publication() +
    theme(legend.position = "bottom")
  
  # Forest plot - Cold
  p2 <- ggplot(df, aes(x = cold_rr, y = analysis, color = level)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbarh(aes(xmin = cold_lo, xmax = cold_hi), 
                   position = position_dodge(width = 0.5), height = 0.2) +
    scale_color_manual(values = c("Intermediate" = COLD_COLOR, "Immediate" = "#2980b9")) +
    labs(x = "Relative Risk", y = "", title = "Cold (P1)", color = "Level") +
    theme_publication() +
    theme(legend.position = "bottom")
  
  combined <- p1 + p2 +
    plot_annotation(
      title = "Summary: Temperature-Mortality Effects Across Sensitivity Analyses",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  save_figure(combined, "fig10_summary_forest", width = 12, height = 6)
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

cat("\n--- Phase 1 Figures ---\n")
plot_exposure_response()
plot_attributable_burden()
plot_yll_summary()

cat("\n--- Phase 2 Figures ---\n")
plot_sensitivity_lag()
plot_harvesting()
plot_heatwave()

cat("\n--- Phase 4 Figures ---\n")
plot_age_stratification()
plot_sex_stratification()
plot_cause_stratification()

cat("\n--- Summary Figure ---\n")
plot_summary_forest()

cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("FIGURE GENERATION COMPLETE\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat(paste(rep("=", 70), collapse=""), "\n")
