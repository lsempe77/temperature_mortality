#!/usr/bin/env Rscript
# =============================================================================
# 04f: TEST WHETHER AC PENETRATION MODIFIES HEAT-MORTALITY ASSOCIATIONS
# =============================================================================
# Uses state-level AC ownership (PNAD 2022) merged to intermediate regions
# Tests AC as predictor of heat RR in meta-regression

suppressPackageStartupMessages({
  library(mixmeta)
  library(jsonlite)
})

cat("=",.rep=70, "\n")
cat("04f: AC PENETRATION META-REGRESSION\n")
cat("=",.rep=70, "\n\n")

# ---- Load existing meta-regression data ----
data_path <- file.path("results", "interaction_metareg_data_intermediate.csv")
df <- read.csv(data_path, stringsAsFactors = FALSE)
cat("Loaded", nrow(df), "region-sex observations\n")

# ---- Load AC data ----
ac_path <- file.path("..", "phase0_data_prep", "covariates", "results", 
                       "ac_ownership_by_state.csv")
ac <- read.csv(ac_path, stringsAsFactors = FALSE)
cat("Loaded AC data for", nrow(ac), "states\n\n")

# ---- Map intermediate region codes to states ----
# First 2 digits of region_code = state code
state_code_map <- data.frame(
  state_code = c(11,12,13,14,15,16,17,
                 21,22,23,24,25,26,27,28,29,
                 31,32,33,35,
                 41,42,43,
                 50,51,52,53),
  state = c('RO','AC','AM','RR','PA','AP','TO',
            'MA','PI','CE','RN','PB','PE','AL','SE','BA',
            'MG','ES','RJ','SP',
            'PR','SC','RS',
            'MS','MT','GO','DF'),
  stringsAsFactors = FALSE
)

# Extract state code from region code
df$state_code <- as.integer(substr(as.character(df$region_code), 1, 2))
df <- merge(df, state_code_map, by = "state_code", all.x = TRUE)
df <- merge(df, ac[, c("state", "ac_pct_2022")], by = "state", all.x = TRUE)

cat("After merging: ", sum(!is.na(df$ac_pct_2022)), "obs with AC data,",
    sum(is.na(df$ac_pct_2022)), "missing\n\n")

# Standardize AC
df$ac_z <- scale(df$ac_pct_2022)[,1]

# ---- Descriptive: AC by region ----
cat("AC ownership by macro-region:\n")
cat("-",.rep=50, "\n")
agg <- aggregate(ac_pct_2022 ~ macro_region, data = df, 
                 FUN = function(x) round(mean(x, na.rm=TRUE), 1))
print(agg)

cat("\nCorrelations with AC:\n")
cat("  AC vs HDI:    r =", round(cor(df$ac_pct_2022, df$hdi, use="complete"), 3), "\n")
cat("  AC vs GDP:    r =", round(cor(df$ac_pct_2022, df$gdp, use="complete"), 3), "\n")
cat("  AC vs Urban:  r =", round(cor(df$ac_pct_2022, df$urbanization_rate, use="complete"), 3), "\n")
cat("  AC vs log(RR heat):  r =", round(cor(df$ac_pct_2022, df$log_rr_heat, use="complete"), 3), "\n")
cat("  AC vs log(RR cold):  r =", round(cor(df$ac_pct_2022, df$log_rr_cold, use="complete"), 3), "\n")

# ---- Model 1: AC alone predicts heat RR ----
cat("\n", "=",.rep=70, "\n")
cat("MODEL 1: Heat ~ AC (bivariate)\n")
cat("=",.rep=70, "\n")
m1_heat <- mixmeta(log_rr_heat ~ ac_z, S = se_heat^2, data = df, method = "reml")
cat("\n")
print(summary(m1_heat))

# ---- Model 2: AC + sex ----
cat("\n", "=",.rep=70, "\n")
cat("MODEL 2: Heat ~ AC + sex\n")
cat("=",.rep=70, "\n")
m2_heat <- mixmeta(log_rr_heat ~ ac_z + female, S = se_heat^2, data = df, method = "reml")
cat("\n")
print(summary(m2_heat))

# ---- Model 3: AC + sex + deprivation + region ----
cat("\n", "=",.rep=70, "\n")
cat("MODEL 3: Heat ~ AC + sex + deprivation + region\n")
cat("=",.rep=70, "\n")
m3_heat <- mixmeta(log_rr_heat ~ ac_z + female + deprivation_z + 
                     relevel(factor(macro_region), ref="Southeast"),
                   S = se_heat^2, data = df, method = "reml")
cat("\n")
print(summary(m3_heat))

# ---- Model 4: AC × sex interaction ----
cat("\n", "=",.rep=70, "\n")
cat("MODEL 4: Heat ~ AC × sex + deprivation + region\n")
cat("=",.rep=70, "\n")
m4_heat <- mixmeta(log_rr_heat ~ ac_z * female + deprivation_z + 
                     relevel(factor(macro_region), ref="Southeast"),
                   S = se_heat^2, data = df, method = "reml")
cat("\n")
print(summary(m4_heat))

# ---- Model 5: Cold ~ AC ----
cat("\n", "=",.rep=70, "\n")
cat("MODEL 5: Cold ~ AC (bivariate)\n")
cat("=",.rep=70, "\n")
m5_cold <- mixmeta(log_rr_cold ~ ac_z, S = se_cold^2, data = df, method = "reml")
cat("\n")
print(summary(m5_cold))

# ---- Model 6: Cold ~ AC + sex + deprivation + region ----
cat("\n", "=",.rep=70, "\n")
cat("MODEL 6: Cold ~ AC + sex + deprivation + region\n")
cat("=",.rep=70, "\n")
m6_cold <- mixmeta(log_rr_cold ~ ac_z + female + deprivation_z + 
                     relevel(factor(macro_region), ref="Southeast"),
                   S = se_cold^2, data = df, method = "reml")
cat("\n")
print(summary(m6_cold))

# ---- AIC comparison ----
cat("\n", "=",.rep=70, "\n")
cat("AIC COMPARISON\n")
cat("=",.rep=70, "\n")
cat("\nHeat models:\n")
cat("  M1 (AC only):            AIC =", AIC(m1_heat), "\n")
cat("  M2 (AC + sex):           AIC =", AIC(m2_heat), "\n")
cat("  M3 (AC + sex + dep + reg): AIC =", AIC(m3_heat), "\n")
cat("  M4 (AC × sex + dep + reg): AIC =", AIC(m4_heat), "\n")
cat("\nCold models:\n")
cat("  M5 (AC only):            AIC =", AIC(m5_cold), "\n")
cat("  M6 (AC + sex + dep + reg): AIC =", AIC(m6_cold), "\n")

# ---- Key question: does AC add to deprivation model? ----
cat("\n", "=",.rep=70, "\n")
cat("KEY COMPARISON: Does AC add to deprivation + region model?\n")
cat("=",.rep=70, "\n")

# Without AC
m_no_ac <- mixmeta(log_rr_heat ~ female + deprivation_z + 
                     relevel(factor(macro_region), ref="Southeast"),
                   S = se_heat^2, data = df, method = "reml")
cat("\nHeat without AC: AIC =", AIC(m_no_ac), "\n")
cat("Heat with AC:    AIC =", AIC(m3_heat), "\n")
cat("Difference:      ", AIC(m_no_ac) - AIC(m3_heat), 
    "(positive = AC improves fit)\n")

m_no_ac_cold <- mixmeta(log_rr_cold ~ female + deprivation_z + 
                          relevel(factor(macro_region), ref="Southeast"),
                        S = se_cold^2, data = df, method = "reml")
cat("\nCold without AC: AIC =", AIC(m_no_ac_cold), "\n")
cat("Cold with AC:    AIC =", AIC(m6_cold), "\n")
cat("Difference:      ", AIC(m_no_ac_cold) - AIC(m6_cold),
    "(positive = AC improves fit)\n")

# ---- Collinearity check ----
cat("\n", "=",.rep=70, "\n")
cat("COLLINEARITY CHECK\n")
cat("=",.rep=70, "\n")
cat("VIF proxy (correlations among predictors):\n")
pred_cor <- cor(df[, c("ac_z", "deprivation_z", "gdp_z", "urban_z", "temp_sd_z")],
                use = "complete")
print(round(pred_cor, 3))

cat("\n\nDONE!\n")
