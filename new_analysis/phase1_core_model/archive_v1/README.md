# Archive - Superseded v1 Scripts

## ⚠️ DO NOT USE THESE SCRIPTS FOR NEW ANALYSIS

These scripts are archived for reproducibility and audit trail purposes.
They have been superseded by v2 versions with improved methodology.

**Active documentation:** `ANALYSIS_ROADMAP.md` (in `new_analysis/`)

---

## Archive Date: December 2025

## Reason for Archiving:

### 1. DLNM Scripts (01a_*_dlnm.py)
- **Issue**: Used polynomial distributed lag (66 parameters)
- **Fix**: v2 uses natural cubic spline cross-basis (16 parameters)
- **Why it matters**: Natural splines are standard in MCC/Gasparrini methods

### 2. Attributable Burden Script (01d_attributable_burden.py)
- **Issue**: Hard-coded polynomial basis reconstruction
- **Fix**: v2 auto-detects basis type (polynomial vs natural spline)
- **Why it matters**: Must match the DLNM output format

### 3. YLL Scripts (01e_*.py, 01e2_*.py)
- **Issue**: Multiple fragmented scripts with inconsistent approaches
- **Fix**: Unified script with CLI flags (--mode actual/assumed/both)
- **Why it matters**: Robust age parsing, single life table loader

### 4. Excess Mortality Script (01g_excess_mortality.py)
- **Issue**: Baseline too simple (month dummies + quadratic), no population offset,
  no uncertainty quantification, no diagnostics, no flu/holiday adjustment
- **Fix**: v2 uses flexible splines, Negative Binomial, population offset,
  prediction intervals, diagnostics, sensitivity analyses
- **Why it matters**: Credible validation of DLNM estimates requires proper baseline

### 5. Case-Crossover Script (01f_case_crossover.py)
- **Issue**: NOT a true case-crossover - aggregated daily counts, compared high vs low
  mortality days. This is ecological correlation, not self-matched design.
- **Fix**: v2 uses individual death records, builds proper matched dataset
  (same year-month-DOW controls), fits conditional logistic regression
- **Why it matters**: True case-crossover removes all time-invariant confounding

---

## Archived Files:

| Original | Superseded By | Key Difference |
|----------|---------------|----------------|
| `01a_intermediate_dlnm.py` | `01a_intermediate_dlnm_v2.py` | Natural spline, MMT |
| `01a_immediate_dlnm.py` | `01a_immediate_dlnm_v2.py` | Natural spline, MMT |
| `01d_attributable_burden.py` | `01d_attributable_burden_v2.py` | Basis auto-detect |
| `01e_yll_calculation.py` | `01e_yll_unified.py` | Unified, robust |
| `01e2_yll_with_age_distribution.py` | `01e_yll_unified.py` | Unified, robust |
| `01e_yll_calculation_v2.py` | `01e_yll_unified.py` | Unified, robust |
| `01g_excess_mortality.py` | `01g_excess_mortality_v2.py` | Robust baseline |
| `01f_case_crossover.py` | `01f_case_crossover_v2.py` | TRUE case-crossover |

---

## Current Pipeline (USE THESE):

```
01a_*_dlnm_v2.py → 01d_attributable_burden_v2.py → 01e_yll_unified.py
```

---

## If You Need to Reproduce Old Results:

1. These archived scripts can still run
2. Check Git history for exact version used
3. Document which version was used in analysis

---

Last archived by: Automated cleanup
Archive reason: Expert code review - methodology upgrade
