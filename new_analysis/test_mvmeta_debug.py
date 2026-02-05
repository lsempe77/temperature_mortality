"""
Debug MVMeta variance issue
"""
import json
import numpy as np
import sys
sys.path.insert(0, '.')
from utils.dlnm_module import mvmeta_fit, mvmeta_predict_rr

# Load Phase 1 results
print("Loading Phase 1 results...")
with open('phase1_core_model/results/dlnm_v2_intermediate_results.json') as f:
    results = json.load(f)

# Get 10 regions with valid coefficients
coef_list = []
vcov_list = []
cb_meta = None
pcts = []

for code, r in list(results['region_results'].items())[:50]:
    if r.get('cb_coefs') and r.get('cb_vcov'):
        coef_list.append(np.array(r['cb_coefs']))
        vcov_list.append(np.array(r['cb_vcov']))
        pcts.append(r.get('temp_percentiles', {}))
        if cb_meta is None:
            cb_meta = r.get('crossbasis_info', {})
    if len(coef_list) >= 10:
        break

print(f"Regions loaded: {len(coef_list)}")
print(f"Coef dimensions: {coef_list[0].shape}")
print(f"cb_meta keys: {list(cb_meta.keys())}")

# Pool using MVMeta
print("\nPooling with MVMeta...")
mvmeta_result = mvmeta_fit(coef_list, vcov_list)
print(f"MVMeta converged: {mvmeta_result['converged']}")
print(f"Pooled coef shape: {np.array(mvmeta_result['pooled_coef']).shape}")
print(f"Pooled vcov shape: {np.array(mvmeta_result['pooled_vcov']).shape}")

# Check vcov for issues
pooled_coef = np.array(mvmeta_result['pooled_coef'])
pooled_vcov = np.array(mvmeta_result['pooled_vcov'])
print(f"\nPooled coef range: {pooled_coef.min():.6f} to {pooled_coef.max():.6f}")
print(f"Vcov diagonal range: {pooled_vcov.diagonal().min():.6f} to {pooled_vcov.diagonal().max():.6f}")
print(f"Any negative diagonal: {np.any(pooled_vcov.diagonal() < 0)}")
print(f"Any NaN in vcov: {np.any(np.isnan(pooled_vcov))}")
print(f"Any Inf in vcov: {np.any(np.isinf(pooled_vcov))}")

# Get median percentiles
p1_vals = [p.get('p1', 10) for p in pcts]
p50_vals = [p.get('p50', 25) for p in pcts]
p99_vals = [p.get('p99', 35) for p in pcts]
p1 = np.median(p1_vals)
p50 = np.median(p50_vals)
p99 = np.median(p99_vals)
print(f"\nPercentiles: P1={p1:.1f}, P50={p50:.1f}, P99={p99:.1f}")

# Test prediction
print("\nTesting mvmeta_predict_rr...")
print(f"cb_meta has var_formula: {'var_formula' in cb_meta}")
print(f"cb_meta has lag_basis_vals: {'lag_basis_vals' in cb_meta}")

try:
    rr, lo, hi, log_rr, se = mvmeta_predict_rr(pooled_coef, pooled_vcov, cb_meta, p99, p50)
    print(f"\nHeat (P99 vs P50): RR={rr:.4f} ({lo:.4f}-{hi:.4f})")
    print(f"  log_rr={log_rr:.6f}, se={se:.6f}")
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()

try:
    rr, lo, hi, log_rr, se = mvmeta_predict_rr(pooled_coef, pooled_vcov, cb_meta, p1, p50)
    print(f"\nCold (P1 vs P50): RR={rr:.4f} ({lo:.4f}-{hi:.4f})")
    print(f"  log_rr={log_rr:.6f}, se={se:.6f}")
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()

# Debug the contrast computation
print("\n--- Debugging contrast computation ---")
from utils.dlnm_module import _safe_dmatrix

var_formula = cb_meta.get('var_formula')
lag_formula = cb_meta.get('lag_formula')
var_boundary = cb_meta.get('var_boundary', cb_meta.get('temp_boundary'))
lag_basis_vals = cb_meta.get('lag_basis_vals')

print(f"var_formula: {var_formula}")
print(f"lag_formula: {lag_formula}")
print(f"var_boundary: {var_boundary}")
print(f"lag_basis_vals type: {type(lag_basis_vals)}")

if var_formula:
    temp_range = np.linspace(var_boundary[0], var_boundary[1], 100)
    temp_basis_full = _safe_dmatrix(var_formula, temp_range, "x")
    print(f"temp_basis shape: {temp_basis_full.shape}")
    
    target_idx = np.argmin(np.abs(temp_range - p99))
    ref_idx = np.argmin(np.abs(temp_range - p50))
    
    temp_basis_target = temp_basis_full.iloc[target_idx].to_numpy()
    temp_basis_ref = temp_basis_full.iloc[ref_idx].to_numpy()
    temp_diff = temp_basis_target - temp_basis_ref
    
    print(f"temp_diff: {temp_diff}")
    print(f"temp_diff norm: {np.linalg.norm(temp_diff):.6f}")

if lag_basis_vals is not None:
    lag_basis = np.array(lag_basis_vals)
    lag_sum = lag_basis.sum(axis=0)
    print(f"lag_basis shape: {lag_basis.shape}")
    print(f"lag_sum: {lag_sum}")
    print(f"lag_sum norm: {np.linalg.norm(lag_sum):.6f}")
    
    # Build contrast
    K_var = len(temp_diff)
    K_lag = len(lag_sum)
    contrast = np.zeros(len(pooled_coef))
    
    for t_idx in range(K_var):
        for l_idx in range(K_lag):
            col_idx = t_idx * K_lag + l_idx
            if col_idx < len(contrast):
                contrast[col_idx] = temp_diff[t_idx] * lag_sum[l_idx]
    
    print(f"\nContrast shape: {contrast.shape}")
    print(f"Contrast norm: {np.linalg.norm(contrast):.6f}")
    print(f"Contrast non-zero: {np.count_nonzero(contrast)}")
    
    # Compute variance
    log_rr = np.dot(contrast, pooled_coef)
    var_log_rr = contrast @ pooled_vcov @ contrast
    se_log_rr = np.sqrt(max(0, var_log_rr))
    
    print(f"\nlog_rr = {log_rr:.6f}")
    print(f"var_log_rr = {var_log_rr:.6f}")
    print(f"se_log_rr = {se_log_rr:.6f}")
