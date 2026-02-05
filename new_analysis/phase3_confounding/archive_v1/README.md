# Phase 3 Archive - Deprecated v1 Scripts

This folder contains deprecated Phase 3 scripts that have been replaced by improved v2 versions.

## Archived Files

### 03a_supplementary_analyses.py
**Archived:** 2025-01-09
**Replaced by:** `../03a_supplementary_analyses_v2.py`

**Critical Issues Fixed in v2:**

1. **Polynomial cross-basis + naive variance** → v2 uses natural spline cross-basis (16 params) + full covariance delta method
2. **Single pooled model with region dummies** → v2 fits region-specific models then pools via DerSimonian-Laird meta-analysis
3. **No population offset** → v2 includes `offset=np.log(pop_elderly)` in all Poisson models
4. **Missing pollution/influenza treated as zeros** → v2 uses smart interpolation with coverage validation (MIN_COVERAGE=0.7)
5. **Global percentiles instead of region-specific** → v2 computes region-specific temperature percentiles
6. **Apparent temperature formula uncertainty** → v2 uses validated Steadman AT formula
7. **Fragile CB column mapping** → v2 uses named parameter lookup via `.reindex()`
8. **Filter/MIN_OBS bias** → v2 has explicit eligibility criteria
9. **CI approximation** → v2 uses full covariance delta method via utils module

**Do not use these scripts for production analysis.**
