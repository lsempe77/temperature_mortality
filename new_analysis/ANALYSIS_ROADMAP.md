# Heat-Mortality Analysis: Brazil 2010-2024

**Last Updated:** December 10, 2025  
**Status:** ✅ Phase 0 Complete — All 17 datasets validated  
**Validation Command:** `python phase0_data_prep/utilities/00a_document_data.py`

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Data Overview](#2-data-overview)
   - 2.1 [Data Summary](#21-data-summary)
   - 2.2 [Geographic Levels](#22-geographic-levels)
   - 2.3 [Daily Panel Data](#23-daily-panel-data)
   - 2.4 [Cross-Sectional Covariates](#24-cross-sectional-covariates)
   - 2.5 [National-Level Data](#25-national-level-data)
   - 2.6 [Key Variables](#26-key-variables)
3. [Phase 0: Data Preparation](#3-phase-0-data-preparation)
   - 3.1 [Aggregation Scripts](#31-aggregation-scripts)
   - 3.2 [All Phase 0 Scripts](#32-all-phase-0-scripts)
   - 3.3 [Pipeline Execution Order](#33-pipeline-execution-order)
   - 3.4 [Output File Summary](#34-output-file-summary)
   - 3.5 [Folder Structure](#35-folder-structure)
4. [Shared Utilities Module](#4-shared-utilities-module)
   - 4.1 [Available Functions](#41-available-functions)
   - 4.2 [Usage Examples](#42-usage-examples)
5. [Phase 1: Core DLNM Analysis](#5-phase-1-core-dlnm-analysis)
   - 5.1 [Scripts](#51-scripts)
   - 5.2 [Output Files](#52-output-files)
6. [Phase 2: Sensitivity Analyses](#6-phase-2-sensitivity-analyses)
   - 6.1 [Scripts](#61-scripts)
   - 6.2 [v2 Improvements](#62-v2-improvements)
7. [Phase 3: Supplementary Analyses](#7-phase-3-supplementary-analyses)
   - 7.1 [Scripts](#71-scripts)
   - 7.2 [v2 Improvements](#72-v2-improvements)
8. [Phase 4: Heterogeneity Analysis](#8-phase-4-heterogeneity-analysis)
   - 8.1 [Scripts](#81-scripts)
   - 8.2 [v2 Improvements](#82-v2-improvements)
   - 8.3 [Output Files](#83-output-files)
9. [Phase 5: Tables & Figures](#9-phase-5-tables--figures)
   - 9.1 [v2 Modular Architecture](#91-v2-modular-architecture)
   - 9.2 [Tables Generated](#92-tables-generated)
   - 9.3 [Usage](#93-usage)
10. [Phase 6: Paper Writing](#10-phase-6-paper-writing)
    - 10.1 [Scripts](#101-scripts)
    - 10.2 [Paper Outputs](#102-paper-outputs)
    - 10.3 [Comprehensive Document](#103-comprehensive-document)

---

## 1. Executive Summary

**Objective:** Estimate heat- and cold-attributable elderly mortality in Brazil using DLNM methodology for publication in a top-tier journal (Lancet Planetary Health, Nature Climate Change, EHP, PNAS).

### Key Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Spatial unit | 133 Intermediate + 510 Immediate Regions | Avoids Berkson error from state-level averaging |
| Time period | 2010-2024 (15 years) | Avoids El Niño dominance of 3-year period |
| Population | Elderly 60+ | Most vulnerable, highest burden |
| Primary model | DLNM (Quasi-Poisson) | Gold standard, MCC-comparable |
| Reference temp | MMT (empirical) | Region-specific minimum mortality temperature |
| Heat/Cold thresholds | P97.5 / P2.5 (PRIMARY), P99/P1 (comparison) | Following Gasparrini et al. (2015) |
| Max lag | 21 days | Captures delayed mortality |
| Pollution adjustment | Supplementary only | Pollution is mediator, not confounder |

---

## 2. Data Overview

### 2.1 Data Summary

| Metric | Value |
|--------|-------|
| Time period | 2010-2024 (15 years, 5,479 days) |
| Geographic coverage | 133 intermediate + 510 immediate regions |
| Municipalities | 5,571 |
| Total validated records | 13,880,407 |
| Total elderly deaths (60+) | 20,358,595 |
| Validation status | ✅ 17/17 datasets passed (December 2025) |

### 2.2 Geographic Levels

| Level | Count | Use Case |
|-------|-------|----------|
| **Intermediate Regions** | 133 | Primary analysis — larger units, more stable estimates |
| **Immediate Regions** | 510 | Secondary analysis — finer spatial resolution |
| Municipalities | 5,571 | Source data level (aggregated upward) |

**File Locations:**
- Phase 0 outputs: `phase0_data_prep/results/`
- Mapping files: `new_analysis/results/`

**CLI usage for spatial levels (Phases 2–4):**
- Most v2 scripts in `phase2_robustness/`, `phase3_confounding/`, and `phase4_heterogeneity/` accept `--level {intermediate, immediate}`.
- Example: `python phase2_robustness/02a_sensitivity_analyses_v2.py --analysis lag --level immediate`.
- Phase 2 results use a `_immediate` suffix for immediate-level outputs; Phase 3–4 store the selected level in the JSON metadata.
 - For copy-paste commands and output names, see `README_LEVELS.md`.

### 2.3 Daily Panel Data

| Dataset | Intermediate File | Rows | Immediate File | Rows |
|---------|-------------------|------|----------------|------|
| Temperature | `era5_intermediate_daily.parquet` | 728,707 | `era5_immediate_daily.parquet` | 2,794,290 |
| Pollution | `cams_intermediate_daily.parquet` | ~726K | `cams_immediate_daily.parquet` | ~2.8M |
| Mortality (All) | `mortality_intermediate_daily.parquet` | 721,096 | `mortality_immediate_daily.parquet` | 2,523,235 |
| Mortality (Elderly) | `mortality_intermediate_daily_elderly.parquet` | 708,767 | `mortality_immediate_daily_elderly.parquet` | 2,302,539 |
| Influenza | `influenza_daily_by_intermediate_region.parquet` | 187,055 | `influenza_daily_by_immediate_region.parquet` | 378,575 |

### 2.4 Cross-Sectional Covariates

| Dataset | Intermediate File | Rows | Immediate File | Rows |
|---------|-------------------|------|----------------|------|
| SES Covariates | `ses_intermediate_covariates.csv` | 133 | `ses_immediate_covariates.csv` | 510 |
| Regional Covariates | `regional_covariates.csv` | 133 | — | — |
| Geography Mapping | `municipality_to_all_regions_map.csv` | 5,571 | (same file) | — |

### 2.5 National-Level Data

| Dataset | File | Rows | Notes |
|---------|------|------|-------|
| Holidays | `brazilian_holidays_daily.parquet` | 5,479 | 2010-2024, includes `is_holiday_week` |
| Life Tables | `ibge_life_tables_combined.parquet` | 1,230 | IBGE Brazil-specific |
| Age-YLL Lookup | `yll_lookup_by_age.csv` | 90 | Life expectancy by single year (0-89) |
| Age-YLL Groups | `yll_lookup_by_age_group.csv` | — | Life expectancy by 5-year groups |

### 2.6 Key Variables

**Temperature (ERA5):**
- `intermediate_code` / `immediate_code`: Region identifier
- `date`: Date
- `temp_mean`, `temp_min`, `temp_max`: Daily temperatures (°C)
- `dewpoint_mean`: For apparent temperature

**Mortality:**
- `deaths_elderly`: Total elderly deaths
- `deaths_elderly_resp`: Respiratory deaths
- `deaths_elderly_cvd`: Cardiovascular deaths
- `deaths_elderly_heat`: Direct heat deaths (ICD T67, X30)

**SES:**
- `pop_total`, `pop_elderly`: Population counts
- `pct_elderly`: % elderly in population
- `gdp_per_capita`: Economic indicator
- `urbanization_rate`: % urban population

---

## 3. Phase 0: Data Preparation

**Status:** Pending  — 17/17 datasets validated (December 2025)

### 3.1 Aggregation Scripts

Located in `phase0_data_prep/aggregation/`. These generate the analysis-ready regional datasets:

| Script | Purpose | Outputs | Status |
|--------|---------|---------|--------|
| `00g4_aggregate_cams_optimized.py` | CAMS pollution → both region levels | `cams_*_daily.parquet` | ✅ VERIFIED |
| `00h2_aggregate_era5_to_regions.py` | ERA5 temperature → both region levels | `era5_*_daily.parquet` | ✅ RAN |
| `00k3_aggregate_ses_to_regions.py` | SES covariates → both region levels | `ses_*_covariates.csv` | ✅ EXISTS |
| `00m2_aggregate_influenza_municipal.py` | Influenza → municipality + both levels | `influenza_daily_by_*.parquet` | ✅ RAN |
| `00m3_process_additional_influenza.py` | Extend influenza 2018→2024 | Updates influenza files | ✅ RAN |
| `00n_aggregate_mortality_to_regions.py` | SIM mortality → both region levels | `mortality_*_daily*.parquet` | ✅ RAN |

### 3.2 All Phase 0 Scripts

| Script | Purpose | Run Frequency |
|--------|---------|---------------|
| `00a_document_data.py` | Validate all datasets | After data updates |
| `00b_download_cams_pollution_v2.py` | Download CAMS PM2.5/O3 | Once per year |
| `00d_brazilian_holidays.py` | Create holiday dataset | Once |
| `00e_pnad_ac_data.py` | Extract AC ownership | Once |
| `00f_state_covariates.py` | State-level covariates | Once |
| `00f2_regional_covariates.py` | Map covariates to regions | Once |
| `00h_download_era5_brazil.py` | Download ERA5 temperature | Once per year |
| `00j_download_ibge_mapping.py` | Download region mappings | Once |
| `00k2_download_ses_data_fixed.py` | Download SES indicators | Once |
| `00l_download_ibge_life_tables.py` | Process life tables | Once |
| `00m_download_influenza_data.py` | Download SRAG/influenza | Once per year |

### 3.3 Pipeline Execution Order

**Prerequisites:**
```powershell
cd "c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis"
.\.venv\Scripts\Activate.ps1
```

**Step 1: Download Raw Data**
```bash
python phase0_data_prep/00b_download_cams_pollution_v2.py
python phase0_data_prep/00h_download_era5_brazil.py
python phase0_data_prep/00j_download_ibge_mapping.py
python phase0_data_prep/00k2_download_ses_data_fixed.py
python phase0_data_prep/00l_download_ibge_life_tables.py
python phase0_data_prep/00m_download_influenza_data.py
```

**Step 2: Build Covariates**
```bash
python phase0_data_prep/00d_brazilian_holidays.py
python phase0_data_prep/00e_pnad_ac_data.py
python phase0_data_prep/00f_state_covariates.py
python phase0_data_prep/00f2_regional_covariates.py
```

**Step 3: Aggregate to Regions**
```bash
python phase0_data_prep/aggregation/00g4_aggregate_cams_optimized.py
python phase0_data_prep/aggregation/00h2_aggregate_era5_to_regions.py
python phase0_data_prep/aggregation/00k3_aggregate_ses_to_regions.py
python phase0_data_prep/aggregation/00m2_aggregate_influenza_municipal.py
python phase0_data_prep/aggregation/00m3_process_additional_influenza.py
python phase0_data_prep/aggregation/00n_aggregate_mortality_to_regions.py
```

**Step 4: Validate**
```bash
python phase0_data_prep/utilities/00a_document_data.py  # Should show 17/17 passed
```

### 3.4 Output File Summary

| Output File | Script | Records | Date Range |
|-------------|--------|---------|------------|
| `era5_intermediate_daily.parquet` | 00h2 | 728,707 | 2010-2024 |
| `era5_immediate_daily.parquet` | 00h2 | 2,794,290 | 2010-2024 |
| `cams_intermediate_daily.parquet` | 00g4 | ~726K | 2010-2024 |
| `cams_immediate_daily.parquet` | 00g4 | ~2.8M | 2010-2024 |
| `mortality_intermediate_daily.parquet` | 00n | 721,096 | 2010-2024 |
| `mortality_intermediate_daily_elderly.parquet` | 00n | 708,767 | 2010-2024 |
| `mortality_immediate_daily.parquet` | 00n | 2,523,235 | 2010-2024 |
| `mortality_immediate_daily_elderly.parquet` | 00n | 2,302,539 | 2010-2024 |
| `influenza_daily_by_intermediate_region.parquet` | 00m2+00m3 | 187,055 | 2010-2024 |
| `influenza_daily_by_immediate_region.parquet` | 00m2+00m3 | 378,575 | 2010-2024 |
| `ses_intermediate_covariates.csv` | 00k3 | 133 | Static |
| `ses_immediate_covariates.csv` | 00k3 | 510 | Static |
| `municipality_to_all_regions_map.csv` | 00j | 5,571 | Static |
| `brazilian_holidays_daily.parquet` | 00d | 5,479 | 2010-2024 |
| `ibge_life_tables_combined.parquet` | 00l | 1,230 | 2010-2022 |

**Total validated records: 13,880,407**

### 3.5 Folder Structure

```
phase0_data_prep/
├── aggregation/           # Spatial aggregation scripts
│   ├── 00g4_aggregate_cams_optimized.py
│   ├── 00h2_aggregate_era5_to_regions.py
│   ├── 00k3_aggregate_ses_to_regions.py
│   ├── 00m2_aggregate_influenza_municipal.py
│   ├── 00m3_process_additional_influenza.py
│   ├── 00n_aggregate_mortality_to_regions.py
│   └── archive/
├── utilities/             # Validation and documentation
│   └── 00a_document_data.py
├── results/               # Output files
├── docs/                  # Documentation
├── temp_cams/             # Temporary CAMS downloads
├── temp_life_tables/      # Temporary life table files
└── temp_srag/             # Temporary influenza files
```

---

## 4. Shared Utilities Module

**Location:** `utils/dlnm_module.py`  
**Purpose:** Centralize all DLNM functions for consistency across all phases

> ⚠️ **Important:** All Phase 1-4 scripts import from `utils/dlnm_module.py`. Do NOT define duplicate functions in individual scripts. Edit `dlnm_module.py` to change core DLNM logic.

### 4.1 Available Functions

| Category | Function | Description |
|----------|----------|-------------|
| **Spline Basis** | `ns_basis()` | R-style natural cubic spline basis |
| | `create_lag_matrix()` | Create lagged exposure matrix |
| | `create_crossbasis_ns()` | Tensor product cross-basis (Phase 1 style) |
| | `create_crossbasis()` | Patsy-based cross-basis (Phase 2 style) |
| **Prediction** | `compute_cumulative_rr_ns()` | Cumulative RR with CI |
| | `compute_cumulative_rr_ns_with_se()` | Cumulative RR with SE (Phase 1) |
| | `predict_cumulative_rr()` | Predict from fit_res object (Phase 2) |
| | `compute_rr_curve()` | Full RR curve across temperature range |
| **MMT & Effects** | `find_mmt()` | Find MMT from fit_res object |
| | `find_mmt_from_coefficients()` | Find MMT from raw coefficients |
| | `compute_effects_relative_to_mmt()` | Effects at all percentiles vs MMT |
| **Model Fitting** | `fit_region_dlnm()` | Fit region-specific DLNM model |
| | `fit_region_dlnm_with_heatwave()` | DLNM with cb×heatwave interaction |
| **Meta-Analysis** | `meta_random_effects()` | DerSimonian-Laird RE meta-analysis |
| | `pool_region_results()` | Pool region effects with filters |
| **Burden** | `compute_attributable_fraction()` | AF = (RR-1)/RR |
| **Harvesting** | `harvesting_for_region()` | Fit DLNM, compute RR at horizons |
| | `compute_harvesting_ratio()` | Harvesting ratio from ERR |
| **Heatwave** | `identify_heatwaves()` | Mark heatwave days (P90, ≥2 consec) |

### 4.2 Usage Examples

**Phase 1 scripts:**
```python
from utils.dlnm_module import (
    ns_basis, create_crossbasis_ns, create_lag_matrix,
    compute_cumulative_rr_ns_with_se, find_mmt_from_coefficients,
    random_effects_meta_analysis, pool_region_results,
)
```

**Phase 2 scripts:**
```python
from utils.dlnm_module import (
    create_crossbasis, fit_region_dlnm, predict_cumulative_rr,
    meta_random_effects, pool_region_effects,
    harvesting_for_region, identify_heatwaves,
)
```

**Key Implementation Notes:**
1. Two cross-basis styles: `create_crossbasis_ns()` (Phase 1) vs `create_crossbasis()` (Phase 2)
2. MMT is empirical — each region has its own minimum mortality temperature
3. Delta method for standard errors using full covariance matrix
4. DerSimonian-Laird random-effects meta-analysis with I² heterogeneity

---

## 5. Phase 1: Core DLNM Analysis

**Goal:** Estimate temperature-mortality relationships at two spatial scales

> All v1 scripts archived → `phase1_core_model/archive_v1/`  
> Use v2 scripts for all new analysis

### 5.1 Scripts

| Script | Description | Status | Summary |
|--------|-------------|--------|---------|
| `01a_intermediate_dlnm_v2.py` | Natural spline DLNM (133 regions) | ✅ COMPLETE | INTERMEDIATE_DLNM_V2_SUMMARY.md
| `01a_immediate_dlnm_v2.py` | Natural spline DLNM (510 regions) | ✅ v2 (primary spatial scale) |
| `01d_attributable_burden_v2.py` | Burden (P2.5/P97.5, MMT reference; **uses both intermediate & immediate DLNM results**) | ✅ v2 |
| `01e_yll_unified.py` | Unified YLL (actual/assumed modes; **applied to both levels using burden_v2_national_summary.json**) | ✅ v2 |
| `01f_case_crossover_v2.py` | Case-crossover validation (national, compares with pooled intermediate DLNM) | ✅ v2 |
| `01g_excess_mortality_v2.py` | Excess mortality validation (**compares excess deaths with both intermediate & immediate attributable burdens**) | ✅ v2 |

**v2 Key Features:**
- Natural cubic spline cross-basis (16 params vs 66 in polynomial v1)
- Empirical MMT per region (not fixed P50)
- Thresholds: P2.5/P97.5 (PRIMARY), P1/P99 (COMPARISON)
- Population offset: `log(pop_elderly / 100000)`
- Controls: month, day_of_week, time_spline (1 df/year), holidays

**Latest DLNM v2 results (Dec 2025, immediate level):**
- Regions: 510/510 fitted successfully; exclusions: 0; elderly death coverage ≈99.7%.
- MMT distribution across immediate regions: mean percentile 42.8, median 40.3, range P1–P99.
- Pooled RRs vs region-specific MMT (immediate level):
   - Heat: P75 RR ≈ 1.43 [1.34–1.52]; P95 RR ≈ 3.52 [2.98–4.15]; P97.5 RR ≈ 3.03 [2.54–3.63].
   - Cold: P10 RR ≈ 1.65 [1.54–1.78]; P5 RR ≈ 1.73 [1.62–1.86]; P2.5 RR ≈ 1.41 [1.35–1.48].
   - P99 and P1 pooled RRs are treated as diagnostic only (too few regions at those extremes for a stable pooled estimate).

**Pipeline:**
```
01a_*_dlnm_v2.py → 01d_attributable_burden_v2.py → 01e_yll_unified.py
```

All Phase 1 burden/YLL/validation scripts now consume **both** intermediate and immediate DLNM v2 results, with immediate regions treated as the primary spatial scale for headline estimates.

### 5.2 Output Files

Located in `phase1_core_model/results/`

- Summary of immediate-level DLNM v2 run: `DLNM_V2_IMMEDIATE_SUMMARY.md`

---

## 6. Phase 2: Sensitivity Analyses

**Purpose:** Test whether main findings are robust to analytical choices

> All v1 scripts archived → `phase2_robustness/archive_v1/`  
> Use v2 scripts for all new analysis

### 6.1 Scripts

| Script | Analysis | Status |
|--------|----------|--------|
| `02a_sensitivity_analyses_v2.py` | Comprehensive sensitivity tests | ✅ v2; `--level` intermediate/immediate |
| `02b_harvesting_analysis_v2.py` | Mortality displacement | ✅ v2; `--level` intermediate/immediate |
| `02c_heatwave_dlnm_v2.py` | Heatwave effect modification | ✅ v2; `--level` intermediate/immediate |

### 6.2 v2 Improvements

1. Uses `utils/dlnm_module.py` with natural cubic spline cross-basis
2. Fits per-region models with population offset, then meta-analyzes
3. Proper delta-method confidence intervals for cumulative RR
4. Harvesting computed from Excess Relative Risk (ERR = RR - 1) at multiple horizons
5. Heatwave analysis uses proper cb × heatwave interaction terms
6. Unified CLI: `--level {intermediate, immediate}` switches between 133-region and 510-region analyses; immediate-level outputs are saved with a `_immediate` suffix (e.g., `sensitivity_analyses_v2_immediate.json`).

---

## 7. Phase 3: Supplementary Analyses

**Purpose:** Additional analyses for supplement (not primary causal estimates)

> All v1 scripts archived → `phase3_confounding/archive_v1/`  
> Use v2 scripts for all new analysis

### 7.1 Scripts

| Script | Analysis | Status |
|--------|----------|--------|
| `03a_supplementary_analyses_v2.py` | All supplementary analyses | ✅ v2; `--level` intermediate/immediate |

**Analyses included:**
1. **Apparent Temperature** — Steadman AT with humidity adjustment
2. **Pollution-Adjusted** — Controls for PM2.5/O3 (supplement only)
3. **Influenza-Adjusted** — Controls for flu deaths
4. **Holiday-Adjusted** — Controls for reduced exposure on holidays

> ⚠️ **Note on Pollution:** Pollution is a **mediator** of heat effects (heat → atmospheric stagnation → pollution → mortality). Primary results are unadjusted; pollution-adjusted results show direct effect only.

### 7.2 v2 Improvements

| Problem in v1 | Fix in v2 |
|---------------|-----------|
| Polynomial cross-basis (66 params) | Natural spline (16 params) |
| Single pooled model | Per-region models → meta-analysis |
| No population offset | `offset=np.log(pop_elderly)` |
| Missing data treated as zeros | Smart interpolation with coverage validation |
| Global percentiles | Region-specific percentiles |

**Output Files:**
- `supplementary_analyses_v2.json` — full nested results for all supplementary models
- `supplementary_analyses_v2_summary.csv` — pooled RR summary table across models/percentiles
- `supplementary_analyses_v2_metadata.json` — configuration and data summary (includes `level`)

---

## 8. Phase 4: Heterogeneity Analysis

**Purpose:** Identify who is most affected and explain regional variation

> All v1 scripts archived → `phase4_heterogeneity/archive_v1/`  
> Use v2 scripts for all new analysis

### 8.1 Scripts

| Script | Analysis | Status |
|--------|----------|--------|
| `04a_meta_regression_v2.py` | Meta-regression | ✅ v2; `--level` intermediate/immediate |
| `04b_age_stratification_v2.py` | Age groups (60-69, 70-79, 80+) | ✅ v2; `--level` intermediate/immediate |
| `04c_sex_stratification_v2.py` | Male vs Female | ✅ v2; `--level` intermediate/immediate |
| `04d_cause_stratification_v2.py` | CVD vs Respiratory vs Other | ✅ v2; `--level` intermediate/immediate |

### 8.2 v2 Improvements

| Problem in v1 | Fix in v2 |
|---------------|-----------|
| Polynomial cross-basis (unstable) | Natural spline (16 params) |
| Pooled GLM with region dummies | Per-region DLNM → meta-analysis |
| No population offset | `offset=log(pop_elderly)` |
| **04d: Filtered zero deaths (BUG)** | Keep all rows including zeros |

**Two-Stage Approach:**
```
For each stratum (age group / sex / cause):
  1. Fit DLNM per region with offset
  2. Extract log-RR and SE at region-specific percentiles
  3. Pool across regions using DerSimonian-Laird meta-analysis
  4. Report I², tau², and heterogeneity tests
```

### 8.3 Output Files

- `meta_regression_v2_results.json`
- `meta_regression_v2_region_effects.json`
- `meta_regression_v2_data.csv`
- `age_stratification_v2_results.json`
- `age_stratification_v2_region_effects.json`
- `sex_stratification_v2_results.json`
- `sex_stratification_v2_region_effects.json`
- `cause_stratification_v2_results.json`
- `cause_stratification_v2_region_effects.json`

---

## 9. Phase 5: Tables & Figures

**Purpose:** Generate all publication-ready outputs

> All v1 scripts archived → `phase5_outputs/archive_v1/`  
> Use v2 scripts for all new analysis

### 9.1 v2 Modular Architecture

**Location:** `phase5_outputs/v2/`

| Module | Purpose |
|--------|---------|
| `loaders.py` | Unified data loading, validation, caching |
| `fig_core.py` | E-R curves, forest plots, sensitivity plots |
| `fig_3d_surface.py` | 3D surfaces from ACTUAL DLNM coefficients |
| `fig_stratification.py` | Age/sex/cause forest plots |
| `tables.py` | Tables 1-5, S1-S2, LaTeX export |
| `generate_all.py` | Main entry point with CLI |

**Key v2 improvement:** Uses actual DLNM cross-basis coefficients (v1 used synthetic data)

### 9.2 Tables Generated

| Table | Description | File |
|-------|-------------|------|
| Table 1 | Descriptive statistics |  |
| Table 2 | MMT and thresholds | |
| Table 3 | Heat/cold effects |  |
| Table 4 | Attributable burden | |
| Table 5 | Sensitivity analysis | `table5_sensitivity.csv` |
| Table S1 | Stratified effects |  |
| Table S2 | Meta-regression |  |

**Output Location:** `phase5_outputs/figures/` and `phase5_outputs/tables/`

### 9.3 Usage

```python
from phase5_outputs.v2.generate_all import Phase5OutputGenerator
generator = Phase5OutputGenerator()
generator.generate_all()
```

**Command line:**
```bash
python phase5_outputs/v2/generate_all.py
python phase5_outputs/v2/generate_all.py --tables-only
python phase5_outputs/v2/generate_all.py --validate-only
```

---

## 10. Phase 6: Paper Writing

**Purpose:** Generate publication-ready manuscript for Lancet Planetary Health

### 10.1 Scripts

| Script | Output | Status |
|--------|--------|--------|
| `06a_generate_paper.py` | Main paper draft |  |
| `06b_supplementary_materials.py` | Supplementary materials | Pending  |
| `06c_verify_consistency.py` | Data consistency verification | Pending  |

### 10.2 Paper Outputs

| Deliverable | File | Status |
|-------------|------|--------|
| Main paper | `paper.qmd` / `paper.pdf` / `paper.html` | Pending |
| Full draft | `full_paper_draft.md` | Pending  |
| Abstract | `abstract.txt` | Pending  |
| Methods | `methods.txt` | Pending  |
| Results | `results.txt` | Pending  |
| Discussion | `discussion.txt` | Pending  |
| Supplement | `supplementary_materials.md` | Pending  |
| LaTeX equations | `equations_latex.md` | Pending  |
| Bibliography | `references.bib` | Pending  |

**Main Paper Tables:**

| Table | File | Description |
|-------|------|-------------|
| Table 1 | `table1_study_characteristics.csv` | Study characteristics |
| Table 2 | `table2_main_results.csv` | Main temperature-mortality results |
| Table 3 | `table3_sensitivity.csv` | Sensitivity analyses |
| Table 4 | `table4_harvesting.csv` | Harvesting analysis |

**Supplementary Tables:**

| Table | File | Description |
|-------|------|-------------|
| Table S1 | `tableS1_percentile_rrs.csv` | RRs at all percentiles |
| Table S2 | `tableS2_lag_sensitivity.csv` | Lag structure sensitivity |
| Table S3 | `tableS3_harvesting_detail.csv` | Detailed harvesting results |
| Table S4 | `tableS4_age_stratification.csv` | Age-stratified results |
| Table S5 | `tableS5_sex_stratification.csv` | Sex-stratified results |
| Table S6 | `tableS6_confounding_analysis.csv` | Confounding assessment |
| Table S7 | `tableS7_yll_age_distribution.csv` | YLL by age distribution |

### 10.3 Comprehensive Document

**Location:** `phase6_paper/comprehensive/`

A complete ~100-page technical report with all methodology and results:

| Section | File | Pages |
|---------|------|-------|
| Master assembly | `_master_document.qmd` | — |
| Front matter | `01_front_matter.qmd` | 3-4 |
| Phase 0 Data | `02_phase0_data.qmd` | 8-10 |
| Phase 1 Methods | `03_phase1_methods_results.qmd` | 15-20 |
| Phase 2 Sensitivity | `04_phase2_sensitivity.qmd` | 12-15 |
| Phase 3 Confounding | `05_phase3_confounding.qmd` | 8-10 |
| Phase 4 Heterogeneity | `06_phase4_heterogeneity.qmd` | 12-15 |
| Phase 5 Figures | `07_phase5_figures_tables.qmd` | 10-12 |
| Discussion | `08_discussion_conclusions.qmd` | 8-10 |
| Appendices | `09_supplementary_appendix.qmd` | 20-25 |
| **TOTAL** | | **96-121 pages** |

**To Render:**
```bash
cd phase6_paper/comprehensive
quarto render _master_document.qmd
```

---

*Document last updated: December 10, 2025. Run `00a_document_data.py` to validate current data status.*
