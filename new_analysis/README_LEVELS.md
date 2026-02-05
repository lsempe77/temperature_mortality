# Spatial Levels & Quick-Run Guide

This short guide summarizes how to run the main **v2** analysis scripts at different spatial levels (intermediate vs immediate) and where their outputs go.

## 1. Spatial Levels

Most v2 scripts in Phases 2–4 accept a `--level` flag:

- `--level intermediate` → 133 **intermediate** regions (default)
- `--level immediate` → 510 **immediate** regions

Scripts that support `--level`:

| Phase | Script | Notes |
|-------|--------|-------|
| 2 | `phase2_robustness/02a_sensitivity_analyses_v2.py` | Comprehensive sensitivity tests |
| 2 | `phase2_robustness/02b_harvesting_analysis_v2.py` | Harvesting (mortality displacement) |
| 2 | `phase2_robustness/02c_heatwave_dlnm_v2.py` | Heatwave effect modification |
| 3 | `phase3_confounding/03a_supplementary_analyses_v2.py` | Apparent temp, pollution, influenza, holidays |
| 4 | `phase4_heterogeneity/04a_meta_regression_v2.py` | Meta-regression of regional effects |
| 4 | `phase4_heterogeneity/04b_age_stratification_v2.py` | Age groups (60–69, 70–79, 80+) |
| 4 | `phase4_heterogeneity/04c_sex_stratification_v2.py` | Male vs Female |
| 4 | `phase4_heterogeneity/04d_cause_stratification_v2.py` | CVD vs Respiratory vs Other |

**Conventions:**
- Phase 2 robustness outputs append `_immediate` for immediate-level results (e.g., `sensitivity_analyses_v2_immediate.json`).
- Phases 3–4 record the chosen level inside their JSON metadata under the `level` field.

Run all commands from the project root:

```powershell
cd "c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data"
.\.venv\Scripts\Activate.ps1
```

---

## 2. Phase 2 – Robustness (Quick Recipes)

### 2.1 Sensitivity Analyses

Intermediate regions (default):

```bash
python new_analysis/phase2_robustness/02a_sensitivity_analyses_v2.py --analysis lag
```

Immediate regions (quick mode on 15 regions):

```bash
python new_analysis/phase2_robustness/02a_sensitivity_analyses_v2.py \
  --analysis lag --level immediate --quick
```

Main outputs (in `new_analysis/phase2_robustness/results/`):
- `sensitivity_analyses_v2.json`
- `sensitivity_analyses_v2_immediate.json`

### 2.2 Harvesting Analysis

Intermediate regions:

```bash
python new_analysis/phase2_robustness/02b_harvesting_analysis_v2.py
```

Immediate regions (quick mode on 10 regions):

```bash
python new_analysis/phase2_robustness/02b_harvesting_analysis_v2.py \
  --level immediate --quick
```

Outputs (in `phase2_robustness/results/`):
- `harvesting_region_results_v2.csv`
- `harvesting_region_results_v2_immediate.csv`
- `harvesting_pooled_results_v2.json`
- `harvesting_pooled_results_v2_immediate.json`

### 2.3 Heatwave Effect Modification

Intermediate regions, moderate definition (default):

```bash
python new_analysis/phase2_robustness/02c_heatwave_dlnm_v2.py
```

Immediate regions, all definitions, quick mode:

```bash
python new_analysis/phase2_robustness/02c_heatwave_dlnm_v2.py \
  --definition all --level immediate --quick
```

Outputs (in `phase2_robustness/results/`):
- `heatwave_region_results_strict_v2.csv`, `..._moderate_...`, `..._lenient_...`
- `heatwave_region_results_strict_v2_immediate.csv`, etc.
- `heatwave_analysis_v2.json`
- `heatwave_analysis_v2_immediate.json`

---

## 3. Phase 3 – Supplementary Analyses

Run all supplementary models (apparent temp, pollution-/flu-/holiday-adjusted):

Intermediate regions:

```bash
python new_analysis/phase3_confounding/03a_supplementary_analyses_v2.py
```

Immediate regions:

```bash
python new_analysis/phase3_confounding/03a_supplementary_analyses_v2.py --level immediate
```

Outputs (in `phase3_confounding/results/`):
- `supplementary_analyses_v2.json`
- `supplementary_analyses_v2_summary.csv`
- `supplementary_analyses_v2_metadata.json` (includes `level`)

---

## 4. Phase 4 – Heterogeneity (Meta-Regression & Stratification)

### 4.1 Meta-Regression

Intermediate regions:

```bash
python new_analysis/phase4_heterogeneity/04a_meta_regression_v2.py
```

Immediate regions:

```bash
python new_analysis/phase4_heterogeneity/04a_meta_regression_v2.py --level immediate
```

Outputs (in `phase4_heterogeneity/results/`):
- `meta_regression_v2_results.json`
- `meta_regression_v2_region_effects.json`
- `meta_regression_v2_data.csv`

### 4.2 Age Stratification

Intermediate regions:

```bash
python new_analysis/phase4_heterogeneity/04b_age_stratification_v2.py
```

Immediate regions:

```bash
python new_analysis/phase4_heterogeneity/04b_age_stratification_v2.py --level immediate
```

Outputs:
- `age_stratification_v2_results.json`
- `age_stratification_v2_region_effects.json`

### 4.3 Sex Stratification

Intermediate regions:

```bash
python new_analysis/phase4_heterogeneity/04c_sex_stratification_v2.py
```

Immediate regions:

```bash
python new_analysis/phase4_heterogeneity/04c_sex_stratification_v2.py --level immediate
```

Outputs:
- `sex_stratification_v2_results.json`
- `sex_stratification_v2_region_effects.json`

### 4.4 Cause-Specific Stratification

Intermediate regions:

```bash
python new_analysis/phase4_heterogeneity/04d_cause_stratification_v2.py
```

Immediate regions:

```bash
python new_analysis/phase4_heterogeneity/04d_cause_stratification_v2.py --level immediate
```

Outputs:
- `cause_stratification_v2_results.json`
- `cause_stratification_v2_region_effects.json`

---

## 5. R Implementation (Primary Analysis Pipeline)

The primary analysis now uses R scripts in the `phase1_r/`, `phase2_r/`, `phase3_r/`, and `phase4_r/` directories, leveraging the `dlnm` and `mixmeta` packages (Gasparrini).

### 5.1 Running the Full Pipeline

Use the batch script to run all phases:

```powershell
cd "c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis"

# Run BOTH intermediate and immediate:
.\run_all_phases.ps1

# Run only intermediate:
.\run_all_phases.ps1 -Level intermediate

# Run only immediate:
.\run_all_phases.ps1 -Level immediate
```

### 5.2 Phase 1 R Scripts

| Script | Purpose | Output |
|--------|---------|--------|
| `01_dlnm_analysis_v2.R` | Main DLNM + meta-analysis | `dlnm_r_{level}_results_v2.json` |
| `01b_attributable_burden.R` | Attributable deaths using DLNM coefficients | `burden_r_{level}.json` |
| `01c_yll_calculation.R` | Years of Life Lost | `yll_r_{level}.json` |
| `01d_case_crossover.R` | **Case-crossover validation with DLNM** | `case_crossover_r_{level}.json` |
| `01e_excess_mortality.R` | Excess mortality validation | `excess_mortality_r_{level}.json` |

### 5.3 Enhanced Case-Crossover Design (01d)

The case-crossover script (`01d_case_crossover.R`) implements an **enhanced validation** following Armstrong et al. (2014):

**Key Features:**
- Uses individual death records from SIM (not aggregated counts)
- Time-stratified matching: same year-month + day-of-week
- **Full 21-day lag structure** (not just lag 0)
- **DLNM crossbasis within clogit** - same spline specification as main analysis
- Direct comparison of RRs between time-series DLNM and case-crossover

**Interpretation:**
- **Similar RRs** → Main DLNM robust to unmeasured time-invariant confounding
- **Different RRs** → Potential confounding in time-series approach

**Output includes:**
```json
{
  "dlnm_case_crossover": {
    "converged": true,
    "mmt": 23.5,
    "rr_heat_p99": {"rr": 1.18, "ci_low": 1.12, "ci_high": 1.24},
    "rr_cold_p1": {"rr": 1.14, "ci_low": 1.10, "ci_high": 1.18}
  },
  "main_dlnm_comparison": {
    "dlnm_rr_heat_p99": 1.16,
    "dlnm_rr_cold_p1": 1.13
  }
}
```

### 5.4 Heterogeneity Statistics

All R scripts now extract comprehensive heterogeneity statistics using `qtest()` from mixmeta:

- **Cochran's Q** - Overall heterogeneity test statistic
- **I²** - Percentage of variability due to heterogeneity (0-100%)
- **H²** - Heterogeneity ratio (Q/df)
- **τ²** - Between-study variance (from Psi matrix)

---

## 6. Notes & Recommended Workflow

1. **Start with intermediate regions** to confirm everything runs end-to-end.
2. Then re-run key scripts at the **immediate level**:
   - `02a`, `02b`, `02c` (robustness)
   - `03a` (supplementary)
   - `04a`–`04d` (heterogeneity)
3. Use `--quick` on robustness scripts for fast checks; drop `--quick` for final production runs.
4. For detailed methodology and data provenance, see `ANALYSIS_ROADMAP.md`.
