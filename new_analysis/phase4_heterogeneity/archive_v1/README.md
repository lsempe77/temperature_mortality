# Phase 4 Archive - Deprecated v1 Scripts

This folder contains deprecated Phase 4 scripts that have been replaced by improved v2 versions.

## Archived Files

### 04a_meta_regression_regional.py
**Replaced by:** `../04a_meta_regression_v2.py`

### 04b_age_stratification_regional.py
**Replaced by:** `../04b_age_stratification_v2.py`

### 04c_sex_stratification_regional.py
**Replaced by:** `../04c_sex_stratification_v2.py`

### 04d_cause_stratification_regional.py
**Replaced by:** `../04d_cause_stratification_v2.py`

---

## Critical Issues Fixed in v2

### Common Problems (All Files)

| Problem in v1 | Fix in v2 |
|---------------|-----------|
| Polynomial cross-basis (unstable, 66+ params) | Natural spline cross-basis (16 params) |
| Naive variance (diagonal only) | Full covariance delta method |
| No population offset | `offset=log(pop_elderly)` in all Poisson models |
| Global percentiles for thresholds | Region-specific temperature percentiles |
| Pooled GLM with region dummies (04b/04c/04d) | Per-region DLNM → meta-analysis |
| Meta-analysis without tau² | DerSimonian-Laird random-effects with tau² |
| Fragile position-based indexing | Named parameter lookup via `.reindex()` |

### Script-Specific Fixes

#### 04a_meta_regression
- v1 used polynomial CB + simplified variance sum
- v2 uses natural spline CB + full covariance delta method
- v2 adds proper RE meta-analysis with tau²

#### 04b_age_stratification
- v1 used pooled model with region dummies
- v2 fits per-region DLNMs for each age group, then meta-analyzes
- v2 uses age-proportional population offset

#### 04c_sex_stratification
- v1 used pooled model with region dummies
- v2 fits per-region DLNMs for each sex, then meta-analyzes
- v2 includes paired within-region comparison (more powerful)

#### 04d_cause_stratification
- **CRITICAL BUG FIXED:** v1 filtered out zero-death days (`df[df[deaths_col] > 0]`)
  - This biased Poisson rate estimation upward
  - v2 keeps all rows including zeros (Poisson handles zeros correctly)
- v2 uses per-region DLNM → meta-analysis
- v2 properly handles sparse cause-specific data

---

## v2 Key Improvements Summary

1. **Per-region DLNM → Meta-analysis** (two-stage approach)
   - Fits DLNM in each region separately
   - Extracts log-RR with SE using full covariance delta method
   - Pools using DerSimonian-Laird random-effects meta-analysis
   - Properly accounts for between-region heterogeneity (tau²)

2. **Region-specific percentiles**
   - Heat threshold (P97.5) computed per region
   - Cold threshold (P2.5) computed per region
   - Reference (P50) computed per region

3. **Proper heterogeneity testing**
   - Cochran's Q test for heterogeneity across strata
   - I² statistic reported
   - Pairwise comparisons with Z-tests

---

**Do not use these scripts for production analysis.**

**Archived:** December 9, 2025
