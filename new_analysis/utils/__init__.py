"""
Utilities module for Heat-Mortality Brazil Analysis.

Contains:
- dlnm_module: Proper DLNM implementation with natural cubic splines,
  region-specific fitting, meta-analysis, and prediction utilities.
  
Two cross-basis implementations are available:
1. R-style ns_basis: create_crossbasis_ns() - Used by Phase 1 scripts
2. Patsy cr() based: create_crossbasis() - Used by Phase 2 scripts

Both produce valid natural spline bases but with different parameterizations.

Phase 1 Compatible Functions:
- compute_cumulative_rr_ns_with_se: Returns (log_rr, se)
- find_mmt_from_coefficients: Works with raw coefficients
- compute_effects_relative_to_mmt: Effects at all percentiles vs MMT
- create_time_spline: Long-term trend spline
- random_effects_meta_analysis: DerSimonian-Laird meta-analysis
- pool_region_results: Pool region effects via meta-analysis
- convert_to_json_serializable: JSON serialization helper
- compute_attributable_fraction: AF from RR
- detect_basis_type: Detect NS vs polynomial basis
- compute_cumulative_rr_for_burden: RR for burden calculation
"""

from .dlnm_module import (
    # R-style natural spline (Phase 1 compatible)
    ns_basis,
    create_crossbasis_ns,
    compute_cumulative_rr_ns,
    
    # Phase 1 compatible functions
    compute_cumulative_rr_ns_with_se,
    find_mmt_from_coefficients,
    compute_effects_relative_to_mmt,
    create_time_spline,
    random_effects_meta_analysis,
    pool_region_results,
    convert_to_json_serializable,
    
    # Attributable burden functions
    compute_attributable_fraction,
    detect_basis_type,
    compute_cumulative_rr_for_burden,
    
    # Low-level helpers
    create_lag_matrix,
    natural_spline_basis,
    
    # Patsy-based cross-basis (Phase 2)
    create_crossbasis,
    
    # Model fitting
    fit_region_dlnm,
    fit_region_dlnm_with_heatwave,
    
    # Prediction
    predict_cumulative_rr,
    compute_rr_curve,
    find_mmt,
    
    # Meta-analysis
    meta_random_effects,
    pool_region_effects,
    
    # Harvesting
    harvesting_for_region,
    compute_harvesting_ratio,
    
    # Heatwave
    identify_heatwaves,
    
    # Utilities
    extract_region_percentiles,
)

__all__ = [
    # Phase 1 compatible (R-style)
    'ns_basis',
    'create_crossbasis_ns',
    'compute_cumulative_rr_ns',
    'compute_cumulative_rr_ns_with_se',
    'find_mmt_from_coefficients',
    'compute_effects_relative_to_mmt',
    'create_time_spline',
    'random_effects_meta_analysis',
    'pool_region_results',
    'convert_to_json_serializable',
    # Attributable burden
    'compute_attributable_fraction',
    'detect_basis_type',
    'compute_cumulative_rr_for_burden',
    # Low-level
    'create_lag_matrix', 
    'natural_spline_basis',
    # Phase 2
    'create_crossbasis',
    'fit_region_dlnm',
    'fit_region_dlnm_with_heatwave',
    'predict_cumulative_rr',
    'compute_rr_curve',
    'find_mmt',
    'meta_random_effects',
    'pool_region_effects',
    'harvesting_for_region',
    'compute_harvesting_ratio',
    'identify_heatwaves',
    'extract_region_percentiles',
]
