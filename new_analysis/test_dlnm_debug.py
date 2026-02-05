"""Debug script for DLNM module."""
import pandas as pd
import numpy as np
import sys
sys.path.insert(0, 'new_analysis')
from utils.dlnm_module import fit_region_dlnm, predict_cumulative_rr, find_mmt

PHASE0_RESULTS = 'new_analysis/phase0_data_prep/results'
df_temp = pd.read_parquet(f'{PHASE0_RESULTS}/era5_intermediate_daily.parquet')
df_mort = pd.read_parquet(f'{PHASE0_RESULTS}/mortality_regional_daily_elderly.parquet')
df_ses = pd.read_csv(f'{PHASE0_RESULTS}/ses_intermediate_covariates.csv')
pop_map = dict(zip(df_ses['intermediate_code'], df_ses['pop_elderly']))
df_holidays = pd.read_parquet(f'{PHASE0_RESULTS}/brazilian_holidays_daily.parquet')

df_temp['date'] = pd.to_datetime(df_temp['date'])
df_mort['date'] = pd.to_datetime(df_mort['date'])
df_holidays['date'] = pd.to_datetime(df_holidays['date'])

df = pd.merge(df_mort, df_temp[['date','region_code','temp_mean']], on=['date','region_code'])
df = pd.merge(df, df_holidays[['date','is_holiday']], on='date', how='left')
df['is_holiday'] = df['is_holiday'].fillna(0).astype(int)
df['pop_elderly'] = df['region_code'].map(pop_map)

region = 1301
df_region = df[df['region_code'] == region].copy()
print(f'Region: {region}, N={len(df_region)}')

fit = fit_region_dlnm(df_region, temp_col='temp_mean', deaths_col='deaths_elderly', 
                      pop_col='pop_elderly', max_lag=21, var_df=4, lag_df=4)
print(f'Fit success: {fit is not None}')

if fit:
    print(f"\ncb_colnames from fit: {fit['cb_colnames']}")
    cb_params = [c for c in fit['params'].index if c.startswith('cb_')]
    print(f"cb params in model: {cb_params}")
    print(f"K_var={fit['cb_meta']['K_var']}, K_lag={fit['cb_meta']['K_lag']}")
    print(f"\nvar_knots: {fit['cb_meta'].get('var_knots')}")
    print(f"var_boundary: {fit['cb_meta'].get('var_boundary')}")
    print(f"lag_knots_log: {fit['cb_meta'].get('lag_knots_log')}")
    print(f"lag_boundary_log: {fit['cb_meta'].get('lag_boundary_log')}")
    
    # Compare basis outputs
    from utils.dlnm_module import ns_basis, _safe_dmatrix
    
    temps = df_region['temp_mean'].dropna()
    p50, p99, p01 = temps.quantile(0.50), temps.quantile(0.99), temps.quantile(0.01)
    print(f"\nP01={p01:.1f}, P50={p50:.1f}, P99={p99:.1f}")
    
    # Test temperature at P99
    target_temp = p99
    var_knots = np.array(fit['cb_meta']['var_knots'])
    var_boundary = tuple(fit['cb_meta']['var_boundary'])
    var_formula = fit['cb_meta']['var_formula']
    
    # ns_basis result
    ns_result = ns_basis(np.array([target_temp]), var_knots, var_boundary)
    print(f"\nns_basis at P99={target_temp:.2f}: {ns_result[0]}")
    
    # patsy result (on full temp array to get proper knots, then extract single value)
    all_temps = temps.values
    patsy_full = _safe_dmatrix(var_formula, all_temps, "x")
    # Find the row closest to target_temp
    idx = np.argmin(np.abs(all_temps - target_temp))
    print(f"patsy at similar temp ({all_temps[idx]:.2f}): {patsy_full.iloc[idx].values}")
    
    # The issue: patsy can't evaluate at a single point, so let's verify the knots
    print(f"\nPatsy formula: {var_formula}")
    print(f"Patsy output shape: {patsy_full.shape}")
    print(f"ns_basis output shape: {ns_result.shape}")
    
    # Check coefficient magnitudes
    cb_coefs = [fit['params'][c] for c in cb_params]
    print(f"\nCB coefficients range: {min(cb_coefs):.4f} to {max(cb_coefs):.4f}")
    print(f"CB coefficient mean abs: {np.mean(np.abs(cb_coefs)):.4f}")
    
    # Now test prediction with the fixed method
    print("\n--- Testing prediction ---")
    mmt = find_mmt(fit, np.linspace(temps.quantile(0.05), temps.quantile(0.95), 50), p50)
    print(f"MMT found: {mmt:.2f}")
    
    rr, lo, hi, log_rr = predict_cumulative_rr(fit, p99, mmt)
    print(f"RR(P99 vs MMT): {rr:.4f} ({lo:.4f}-{hi:.4f})")
    print(f"log_rr: {log_rr:.4f}")
    
    se = (np.log(hi) - np.log(lo)) / (2 * 1.96) if hi > lo > 0 else np.nan
    print(f"SE: {se}")
