################################################################################
# 05c: GENERATE MAPS FOR ALL PHASES (R VERSION)
################################################################################
# Creates geographic visualizations of temperature-mortality effects
# 
# Maps are aggregated at the proper level:
# - Intermediate: 133 health regions
# - Immediate: 510 health regions
#
# Outputs:
# - Map 1: Regional heat effects (P99)
# - Map 2: Regional cold effects (P1)  
# - Map 3: Minimum mortality temperature (MMT)
# - Map 4: Attributable fraction by region
#
# Author: Temperature-Mortality Brazil Analysis Pipeline
# Date: December 2025
################################################################################

library(jsonlite)
library(dplyr)
library(ggplot2)
library(sf)
library(viridis)

# Try to load ggpattern for hatched fills, install if needed
if (!requireNamespace("ggpattern", quietly = TRUE)) {
  cat("Installing ggpattern for hatched fills...\n")
  install.packages("ggpattern", repos = "https://cloud.r-project.org")
}
library(ggpattern)

# Directories - use absolute paths for reliability
base_dir <- "c:/Users/LucasSempe/OneDrive - International Initiative for Impact Evaluation/Desktop/sim_data/new_analysis"
script_dir <- file.path(base_dir, "phase5_outputs")
input_data_dir <- "c:/Users/LucasSempe/OneDrive - International Initiative for Impact Evaluation/Desktop/sim_data/Input_data"
phase0_results <- file.path(base_dir, "phase0_data_prep", "results")

PHASE1_R <- file.path(base_dir, "phase1_r", "results")
OUTPUT_DIR <- file.path(script_dir, "maps")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Shapefile and mapping paths
SHAPEFILE <- file.path(input_data_dir, "brazil_municipalities_2022.gpkg")
REGION_MAP <- file.path(phase0_results, "municipality_to_all_regions_map.csv")

cat(paste(rep("=", 70), collapse=""), "\n")
cat("05c: GENERATE MAPS (R VERSION)\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Helper functions
save_map <- function(p, name, width = 12, height = 10) {
  ggsave(file.path(OUTPUT_DIR, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, paste0(name, ".pdf")), p, width = width, height = height)
  cat("  Saved:", name, "\n")
}

# Theme for maps
theme_map <- function() {
  theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
}

# =============================================================================
# LOAD DATA
# =============================================================================

load_brazil_shapes <- function() {
  cat("\n[Loading Brazil shapefile...]\n")
  
  if (!file.exists(SHAPEFILE)) {
    cat("  ERROR: Shapefile not found:", SHAPEFILE, "\n")
    return(NULL)
  }
  
  brazil <- st_read(SHAPEFILE, quiet = TRUE)
  cat("  Loaded", nrow(brazil), "municipalities\n")
  
  return(brazil)
}

load_region_mapping <- function() {
  cat("\n[Loading region mapping...]\n")
  
  if (!file.exists(REGION_MAP)) {
    cat("  ERROR: Region mapping not found:", REGION_MAP, "\n")
    return(NULL)
  }
  
  mapping <- read.csv(REGION_MAP, stringsAsFactors = FALSE)
  cat("  Loaded mapping for", nrow(mapping), "municipalities\n")
  cat("  Unique intermediate regions:", length(unique(mapping$intermediate_code)), "\n")
  cat("  Unique immediate regions:", length(unique(mapping$immediate_code)), "\n")
  
  return(mapping)
}

# =============================================================================
# PREPARE REGION DATA FROM DLNM RESULTS
# =============================================================================

prepare_region_data <- function(level) {
  cat(sprintf("\n[Preparing %s region data...]\n", level))
  
  dlnm_file <- file.path(PHASE1_R, paste0("dlnm_r_", level, "_results_v2.json"))
  d <- tryCatch({
    fromJSON(dlnm_file)
  }, error = function(e) {
    cat("  ERROR: Could not load DLNM results\n")
    return(NULL)
  })
  
  if (is.null(d)) return(NULL)
  
  df <- d$region_results
  
  if (!is.data.frame(df)) {
    cat("  ERROR: region_results is not a data.frame\n")
    return(NULL)
  }
  
  # Filter to regions included in mvmeta
  if ("include_mvmeta" %in% names(df)) {
    df <- df[df$include_mvmeta == TRUE, ]
  }
  
  df <- df %>%
    select(region_code, mmt, rr_p99, rr_p1, n_deaths) %>%
    rename(heat_rr = rr_p99, cold_rr = rr_p1) %>%
    mutate(region_code = as.integer(region_code))
  
  cat("  Loaded", nrow(df), "regions with valid results\n")
  
  # Get pooled estimates to use for non-converged regions
  pooled <- d$pooled
  # Handle different naming conventions and nested structures
  pooled_heat_raw <- if (!is.null(pooled$rr_heat_p99)) pooled$rr_heat_p99 else pooled$rr_p99
  pooled_cold_raw <- if (!is.null(pooled$rr_cold_p1)) pooled$rr_cold_p1 else pooled$rr_p1
  pooled_mmt <- if (!is.null(pooled$mmt)) pooled$mmt else NA
  
  # Extract RR value from nested structure (list with $rr)
  if (is.list(pooled_heat_raw) && "rr" %in% names(pooled_heat_raw)) {
    pooled_heat <- pooled_heat_raw$rr
  } else if (is.list(pooled_heat_raw) && "estimate" %in% names(pooled_heat_raw)) {
    pooled_heat <- pooled_heat_raw$estimate
  } else {
    pooled_heat <- pooled_heat_raw
  }
  
  if (is.list(pooled_cold_raw) && "rr" %in% names(pooled_cold_raw)) {
    pooled_cold <- pooled_cold_raw$rr
  } else if (is.list(pooled_cold_raw) && "estimate" %in% names(pooled_cold_raw)) {
    pooled_cold <- pooled_cold_raw$estimate
  } else {
    pooled_cold <- pooled_cold_raw
  }
  
  if (length(pooled_heat) > 1) pooled_heat <- pooled_heat[1]
  if (length(pooled_cold) > 1) pooled_cold <- pooled_cold[1]
  
  attr(df, "pooled_heat") <- as.numeric(pooled_heat)
  attr(df, "pooled_cold") <- as.numeric(pooled_cold)
  attr(df, "pooled_mmt") <- as.numeric(pooled_mmt)
  
  cat("  Pooled estimates: Heat RR =", round(as.numeric(pooled_heat), 3), 
      ", Cold RR =", round(as.numeric(pooled_cold), 3), 
      ", MMT =", round(as.numeric(pooled_mmt), 1), "\n")
  
  return(df)
}

# Function to fill non-converged regions with pooled estimates
fill_missing_with_pooled <- function(map_data, data, var_name) {
  pooled_val <- switch(var_name,
    "heat_rr" = attr(data, "pooled_heat"),
    "cold_rr" = attr(data, "pooled_cold"),
    "mmt" = attr(data, "pooled_mmt"),
    NA
  )
  
  if (!is.na(pooled_val)) {
    n_missing <- sum(is.na(map_data[[var_name]]))
    if (n_missing > 0) {
      map_data[[var_name]][is.na(map_data[[var_name]])] <- pooled_val
      # Mark imputed regions
      map_data$imputed <- is.na(map_data$n_deaths)
      cat("    Filled", n_missing, "missing regions with pooled estimate (", round(pooled_val, 3), ")\n")
    }
  }
  return(map_data)
}

# =============================================================================
# AGGREGATE SHAPEFILE TO HEALTH REGIONS
# =============================================================================

aggregate_to_regions <- function(brazil, mapping, level) {
  cat(sprintf("\n[Aggregating shapefile to %s regions...]\n", level))
  
  # Add region codes to municipalities
  code_col <- ifelse(level == "intermediate", "intermediate_code", "immediate_code")
  name_col <- ifelse(level == "intermediate", "intermediate_name", "immediate_name")
  
  # Ensure code_muni is same type in both
  brazil$code_muni <- as.integer(brazil$code_muni)
  mapping$code_muni <- as.integer(mapping$code_muni)
  
  # Join mapping to shapefile
  brazil_with_region <- brazil %>%
    left_join(mapping %>% select(code_muni, !!sym(code_col), !!sym(name_col)), 
              by = "code_muni")
  
  # Check match rate
  n_matched <- sum(!is.na(brazil_with_region[[code_col]]))
  cat("  Matched", n_matched, "of", nrow(brazil), "municipalities to regions\n")
  
  if (n_matched == 0) {
    cat("  ERROR: No municipalities matched\n")
    return(NULL)
  }
  
  # Aggregate to regions
  regions <- brazil_with_region %>%
    filter(!is.na(!!sym(code_col))) %>%
    group_by(region_code = !!sym(code_col)) %>%
    summarise(geometry = st_union(geom), .groups = "drop")
  
  cat("  Aggregated to", nrow(regions), level, "regions\n")
  
  return(regions)
}

# =============================================================================
# CREATE MAPS
# =============================================================================

create_heat_map <- function(regions, data, level) {
  cat("\n[Creating heat effects map...]\n")
  
  # Join DLNM data to regions
  map_data <- regions %>%
    left_join(data, by = "region_code")
  
  n_matched <- sum(!is.na(map_data$heat_rr))
  cat("  Matched", n_matched, "of", nrow(regions), "regions with DLNM data\n")
  
  # Fill missing regions with pooled estimate
  map_data <- fill_missing_with_pooled(map_data, data, "heat_rr")
  
  if (sum(!is.na(map_data$heat_rr)) == 0) {
    cat("  ERROR: No regions with data\n")
    return(NULL)
  }
  
  level_label <- ifelse(level == "intermediate", 
                        paste0("133 Intermediate Health Regions (", n_matched, " converged, ", nrow(regions) - n_matched, " imputed*)"),
                        paste0("510 Immediate Health Regions (", n_matched, " converged, ", nrow(regions) - n_matched, " imputed*)"))
  
  # Identify imputed regions for hatching
  imputed_regions <- map_data %>% filter(imputed == TRUE)
  
  p <- ggplot() +
    # Base layer with color fill
    geom_sf(data = map_data, aes(fill = heat_rr), color = "gray50", size = 0.1) +
    # Hatched overlay on imputed regions
    geom_sf_pattern(
      data = imputed_regions,
      pattern = "stripe",
      pattern_density = 0.3,
      pattern_spacing = 0.015,
      pattern_angle = 45,
      pattern_color = "gray30",
      pattern_fill = NA,
      fill = NA,
      color = "gray30",
      size = 0.3
    ) +
    scale_fill_viridis(
      name = "Heat RR\n(P99 vs MMT)",
      option = "inferno",
      na.value = "gray90",
      limits = c(0.8, 2.0),
      oob = scales::squish
    ) +
    labs(
      title = "Heat Effects on Elderly Mortality (P99)",
      subtitle = level_label,
      caption = "*Hatched areas: non-converged regions using pooled national estimate"
    ) +
    theme_map()
  
  save_map(p, paste0("map1_heat_effects_", level))
  return(p)
}

create_cold_map <- function(regions, data, level) {
  cat("\n[Creating cold effects map...]\n")
  
  map_data <- regions %>%
    left_join(data, by = "region_code")
  
  n_matched <- sum(!is.na(map_data$cold_rr))
  
  # Fill missing regions with pooled estimate
  map_data <- fill_missing_with_pooled(map_data, data, "cold_rr")
  
  level_label <- ifelse(level == "intermediate", 
                        paste0("133 Intermediate Health Regions (", n_matched, " converged, ", nrow(regions) - n_matched, " imputed*)"),
                        paste0("510 Immediate Health Regions (", n_matched, " converged, ", nrow(regions) - n_matched, " imputed*)"))
  
  # Identify imputed regions for hatching
  imputed_regions <- map_data %>% filter(imputed == TRUE)
  
  p <- ggplot() +
    # Base layer with color fill
    geom_sf(data = map_data, aes(fill = cold_rr), color = "gray50", size = 0.1) +
    # Hatched overlay on imputed regions
    geom_sf_pattern(
      data = imputed_regions,
      pattern = "stripe",
      pattern_density = 0.3,
      pattern_spacing = 0.015,
      pattern_angle = 45,
      pattern_color = "gray30",
      pattern_fill = NA,
      fill = NA,
      color = "gray30",
      size = 0.3
    ) +
    scale_fill_viridis(
      name = "Cold RR\n(P1 vs MMT)",
      option = "mako",
      na.value = "gray90",
      limits = c(0.8, 2.5),
      oob = scales::squish
    ) +
    labs(
      title = "Cold Effects on Elderly Mortality (P1)",
      subtitle = level_label,
      caption = "*Hatched areas: non-converged regions using pooled national estimate"
    ) +
    theme_map()
  
  save_map(p, paste0("map2_cold_effects_", level))
  return(p)
}

create_mmt_map <- function(regions, data, level) {
  cat("\n[Creating MMT map...]\n")
  
  map_data <- regions %>%
    left_join(data, by = "region_code")
  
  n_matched <- sum(!is.na(map_data$mmt))
  
  # Fill missing regions with pooled estimate
  map_data <- fill_missing_with_pooled(map_data, data, "mmt")
  
  level_label <- ifelse(level == "intermediate", 
                        paste0("133 Intermediate Health Regions (", n_matched, " converged, ", nrow(regions) - n_matched, " imputed*)"),
                        paste0("510 Immediate Health Regions (", n_matched, " converged, ", nrow(regions) - n_matched, " imputed*)"))
  
  # Identify imputed regions for hatching
  imputed_regions <- map_data %>% filter(imputed == TRUE)
  
  p <- ggplot() +
    # Base layer with color fill
    geom_sf(data = map_data, aes(fill = mmt), color = "gray50", size = 0.1) +
    # Hatched overlay on imputed regions
    geom_sf_pattern(
      data = imputed_regions,
      pattern = "stripe",
      pattern_density = 0.3,
      pattern_spacing = 0.015,
      pattern_angle = 45,
      pattern_color = "gray30",
      pattern_fill = NA,
      fill = NA,
      color = "gray30",
      size = 0.3
    ) +
    scale_fill_viridis(
      name = "MMT (°C)",
      option = "plasma",
      na.value = "gray90",
      limits = c(18, 30)
    ) +
    labs(
      title = "Minimum Mortality Temperature",
      subtitle = level_label,
      caption = "*Hatched areas: non-converged regions using pooled national estimate"
    ) +
    theme_map()
  
  save_map(p, paste0("map3_mmt_", level))
  return(p)
}

create_attributable_map <- function(regions, level) {
  cat("\n[Creating attributable burden map...]\n")
  
  # Load attributable burden data
  burden_file <- file.path(PHASE1_R, paste0("attributable_burden_r_", level, ".json"))
  d <- tryCatch({
    fromJSON(burden_file)
  }, error = function(e) {
    cat("  Warning: Could not load burden file\n")
    return(NULL)
  })
  
  if (is.null(d)) return(NULL)
  
  # Extract region burden
  burden <- as.data.frame(d$region_burden)
  
  if (nrow(burden) == 0) {
    cat("  No burden data\n")
    return(NULL)
  }
  
  burden <- burden %>%
    select(region_code, total_heat_af_pct, total_cold_af_pct) %>%
    mutate(total_af_pct = total_heat_af_pct + total_cold_af_pct)
  
  map_data <- regions %>%
    left_join(burden, by = "region_code")
  
  n_matched <- sum(!is.na(map_data$total_af_pct))
  
  # For attributable burden, mark regions without data as imputed (will use pooled)
  pooled_af <- mean(burden$total_af_pct, na.rm = TRUE)
  map_data <- map_data %>%
    mutate(
      imputed = is.na(total_af_pct),
      total_af_pct = ifelse(is.na(total_af_pct), pooled_af, total_af_pct)
    )
  
  level_label <- ifelse(level == "intermediate", 
                        paste0("133 Intermediate Health Regions (", n_matched, " converged, ", sum(map_data$imputed), " imputed*)"),
                        paste0("510 Immediate Health Regions (", n_matched, " converged, ", sum(map_data$imputed), " imputed*)"))
  
  # Identify imputed regions for hatching
  imputed_regions <- map_data %>% filter(imputed == TRUE)
  
  p <- ggplot() +
    # Base layer with color fill
    geom_sf(data = map_data, aes(fill = total_af_pct), color = "gray50", size = 0.1) +
    # Hatched overlay on imputed regions
    geom_sf_pattern(
      data = imputed_regions,
      pattern = "stripe",
      pattern_density = 0.3,
      pattern_spacing = 0.015,
      pattern_angle = 45,
      pattern_color = "gray30",
      pattern_fill = NA,
      fill = NA,
      color = "gray30",
      size = 0.3
    ) +
    scale_fill_viridis(
      name = "AF (%)",
      option = "turbo",
      na.value = "gray90",
      limits = c(0, 30),
      oob = scales::squish
    ) +
    labs(
      title = "Temperature-Attributable Mortality Fraction",
      subtitle = level_label,
      caption = "*Hatched areas: non-converged regions using pooled national estimate"
    ) +
    theme_map()
  
  save_map(p, paste0("map4_attributable_fraction_", level))
  return(p)
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

# Load base data
brazil <- load_brazil_shapes()
mapping <- load_region_mapping()

if (!is.null(brazil) && !is.null(mapping)) {
  
  for (level in c("intermediate", "immediate")) {
    cat("\n", paste(rep("-", 50), collapse=""), "\n")
    cat("Processing", toupper(level), "level\n")
    cat(paste(rep("-", 50), collapse=""), "\n")
    
    # Aggregate shapefile to this level's regions
    regions <- aggregate_to_regions(brazil, mapping, level)
    
    if (is.null(regions)) {
      cat("  Skipping - could not create regions\n")
      next
    }
    
    # Load DLNM data
    data <- prepare_region_data(level)
    
    if (is.null(data)) {
      cat("  Skipping - no DLNM data\n")
      next
    }
    
    # Create maps at proper regional level
    create_heat_map(regions, data, level)
    create_cold_map(regions, data, level)
    create_mmt_map(regions, data, level)
    create_attributable_map(regions, level)
  }
  
} else {
  cat("\nERROR: Could not load required data. Maps not generated.\n")
}

cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("MAP GENERATION COMPLETE\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat(paste(rep("=", 70), collapse=""), "\n")
