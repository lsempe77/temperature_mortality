# Phase 2 Archive (v1 Scripts)

**Archived:** December 9, 2025

These scripts were replaced by v2 versions that use the proper DLNM implementation from `utils/dlnm_module.py`.

## Archived Scripts

| Script | Issue | Replacement |
|--------|-------|-------------|
| `02b_harvesting_analysis.py` | Used polynomial basis, pooled national model, incorrect cumulative RR | `02b_harvesting_analysis_v2.py` |
| `02c_heatwave_dlnm.py` | Used polynomial basis, incorrect heatwave interaction | `02c_heatwave_dlnm_v2.py` |

## Key Issues with v1 Scripts (From Expert Review)

1. **Polynomial basis instability**: Polynomial bases become numerically unstable at long lags (35 days)
2. **No population offset**: Models did not include `offset=log(pop_elderly)`
3. **Pooled national model**: Fit one model to all data instead of per-region with meta-analysis
4. **Incorrect cumulative RR**: Did not use delta method with full covariance matrix
5. **Harvesting ratio formula**: Used oversimplified ratio instead of cumulative RR at multiple horizons
6. **Heatwave interaction**: Did not properly specify cb × heatwave interaction terms

## v2 Improvements

The v2 scripts use `utils/dlnm_module.py` which provides:
- Natural cubic spline cross-basis via Patsy `cr()` 
- Region-specific model fitting with population offset
- Delta-method cumulative RR with full covariance
- DerSimonian-Laird meta-analysis for pooling
- Proper harvesting analysis at multiple horizons
- Correct cb × heatwave interaction model

## Do Not Use

These archived scripts should NOT be used for analysis. Use the v2 versions instead.
