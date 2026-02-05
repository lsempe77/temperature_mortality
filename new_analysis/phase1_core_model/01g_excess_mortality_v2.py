"""
01g_v2: ROBUST EXCESS MORTALITY VALIDATION
==========================================
Validates DLNM attributable burden using counterfactual excess mortality approach.

v2 IMPROVEMENTS (Dec 2025 expert review):
1. GAM-like flexible baseline: splines for seasonality and long-term trend
2. Negative Binomial to handle overdispersion
3. Population offset for rate modelling
4. Influenza adjustment option
5. Holiday adjustment
6. Predictive intervals for expected counts (parametric bootstrap)
7. Comprehensive diagnostics (dispersion, residual ACF, QQ)
8. Sensitivity to baseline period
9. Proper comparison with DLNM AF (restricted to heat season option)

References:
- Fouillet et al. (2006) - Excess mortality approach for 2003 heatwave
- Gasparrini et al. (2022) - Comparison of burden estimation methods
- Karlinsky & Kobak (2021) - The World Mortality Dataset
"""

import warnings
warnings.filterwarnings('ignore')

import json
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
from typing import Dict, Tuple, Optional
import statsmodels.api as sm
from statsmodels.genmod.generalized_linear_model import GLM
from patsy import dmatrix
from scipy import stats
import argparse

# =============================================================================
# CONFIGURATION
# =============================================================================

class Config:
    """Configuration for excess mortality analysis."""
    
    # Spline degrees of freedom
    SEASONAL_DF = 6          # DoY cyclic spline df
    TREND_DF_PER_YEAR = 0.5  # Time trend spline df per year (more flexible)
    
    # Bootstrap for prediction intervals
    N_BOOTSTRAP = 500
    CI_LEVEL = 0.95
    
    # Diagnostics
    MAX_ACF_LAG = 30
    DISPERSION_THRESHOLD = 1.5  # Use NB if dispersion > this
    
    # Data quality
    MIN_OBS_PER_REGION = 365 * 3  # At least 3 years
    
    # Sensitivity analysis baseline periods
    SENSITIVITY_PERIODS = [
        ('full', None, None),          # Use all data
        ('pre_covid', None, 2019),     # End at 2019
        ('post_2015', 2015, None),     # Start 2015
        ('pre_2020', None, 2019),      # Exclude 2020+
    ]

# =============================================================================
# PATHS
# =============================================================================

SCRIPT_DIR = Path(__file__).parent
PHASE0_RESULTS = SCRIPT_DIR.parent / 'phase0_data_prep' / 'results'
PHASE1_RESULTS = SCRIPT_DIR / 'results'
OUTPUT_DIR = PHASE1_RESULTS
FIGURES_DIR = SCRIPT_DIR / 'figures'
OUTPUT_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def create_spline_basis(df: pd.DataFrame, seasonal_df: int = 6, 
                        trend_df: int = 6) -> pd.DataFrame:
    """
    Create spline basis for seasonal and long-term trend.
    
    Uses B-splines for flexible modelling:
    - Cyclic-like spline for day-of-year (seasonal)
    - Smooth spline for time trend
    """
    # Time numeric (days since start)
    df = df.copy()
    df['time_numeric'] = (df['date'] - df['date'].min()).dt.days
    df['doy'] = df['date'].dt.dayofyear
    
    # Seasonal spline (B-spline on day of year)
    # We use a regular B-spline with knots at key seasonal points
    season_formula = f"bs(doy, df={seasonal_df}, degree=3, include_intercept=False)"
    try:
        season_spline = dmatrix(season_formula, df, return_type='dataframe')
        season_spline.columns = [f'season_{i}' for i in range(season_spline.shape[1])]
    except Exception as e:
        print(f"  Warning: Seasonal spline failed ({e}), using month dummies")
        season_spline = pd.get_dummies(df['date'].dt.month, prefix='month', drop_first=True)
    
    # Trend spline
    trend_formula = f"bs(time_numeric, df={trend_df}, degree=3, include_intercept=False)"
    try:
        trend_spline = dmatrix(trend_formula, df, return_type='dataframe')
        trend_spline.columns = [f'trend_{i}' for i in range(trend_spline.shape[1])]
    except Exception as e:
        print(f"  Warning: Trend spline failed ({e}), using quadratic trend")
        trend_spline = pd.DataFrame({
            'trend_1': df['time_numeric'] / 365.25,
            'trend_2': (df['time_numeric'] / 365.25) ** 2
        })
    
    return pd.concat([season_spline.reset_index(drop=True), 
                      trend_spline.reset_index(drop=True)], axis=1)


def fit_baseline_model(df: pd.DataFrame, 
                       include_flu: bool = False,
                       include_holidays: bool = True,
                       use_nb: bool = True) -> Tuple[object, pd.DataFrame]:
    """
    Fit baseline expected deaths model with flexible splines.
    
    Parameters:
    -----------
    df : DataFrame with date, deaths, pop_elderly columns
    include_flu : Include influenza as covariate
    include_holidays : Include holiday indicator
    use_nb : Use Negative Binomial (vs Poisson with overdispersion)
    
    Returns:
    --------
    model_result, X_design
    """
    # Calculate trend df based on years
    n_years = (df['date'].max() - df['date'].min()).days / 365.25
    trend_df = max(3, int(n_years * Config.TREND_DF_PER_YEAR))
    
    # Create spline basis
    spline_basis = create_spline_basis(df, seasonal_df=Config.SEASONAL_DF, 
                                        trend_df=trend_df)
    
    # Day of week dummies
    dow_dummies = pd.get_dummies(df['date'].dt.dayofweek, prefix='dow', drop_first=True)
    dow_dummies = dow_dummies.reset_index(drop=True)
    
    # Build design matrix
    X_parts = [spline_basis, dow_dummies]
    
    # Optional: holidays
    if include_holidays and 'is_holiday' in df.columns:
        X_parts.append(df[['is_holiday']].reset_index(drop=True).astype(float))
    
    # Optional: influenza
    if include_flu and 'srag_cases' in df.columns:
        flu_col = df[['srag_cases']].reset_index(drop=True).fillna(0).astype(float)
        # Log transform + 1 to handle zeros
        flu_col['log_srag'] = np.log1p(flu_col['srag_cases'])
        X_parts.append(flu_col[['log_srag']])
    
    X = pd.concat(X_parts, axis=1)
    X = sm.add_constant(X)
    X = X.astype(float)
    
    y = df['deaths'].values.astype(float)
    
    # Population offset (log rate modelling)
    if 'pop_elderly' in df.columns:
        offset = np.log(df['pop_elderly'].values.astype(float))
    else:
        offset = None
        print("  Warning: No population column, not using offset")
    
    # First fit Poisson to check dispersion
    if offset is not None:
        model_pois = GLM(y, X, family=sm.families.Poisson(), offset=offset).fit()
    else:
        model_pois = GLM(y, X, family=sm.families.Poisson()).fit()
    
    dispersion = model_pois.pearson_chi2 / model_pois.df_resid
    
    # Decide on model family
    if use_nb and dispersion > Config.DISPERSION_THRESHOLD:
        print(f"  Dispersion = {dispersion:.2f} > {Config.DISPERSION_THRESHOLD}, using Negative Binomial")
        try:
            if offset is not None:
                model = GLM(y, X, family=sm.families.NegativeBinomial(alpha=1.0), 
                           offset=offset).fit()
            else:
                model = GLM(y, X, family=sm.families.NegativeBinomial(alpha=1.0)).fit()
        except Exception as e:
            print(f"  Warning: NB failed ({e}), using Quasi-Poisson")
            if offset is not None:
                model = GLM(y, X, family=sm.families.Poisson(), offset=offset).fit(scale='X2')
            else:
                model = GLM(y, X, family=sm.families.Poisson()).fit(scale='X2')
    else:
        print(f"  Dispersion = {dispersion:.2f}, using Poisson with scale adjustment")
        if offset is not None:
            model = GLM(y, X, family=sm.families.Poisson(), offset=offset).fit(scale='X2')
        else:
            model = GLM(y, X, family=sm.families.Poisson()).fit(scale='X2')
    
    return model, X


def compute_prediction_intervals(model, X: pd.DataFrame, 
                                  n_bootstrap: int = 500,
                                  ci_level: float = 0.95) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute prediction intervals for expected counts using parametric bootstrap.
    
    Returns lower and upper CI bounds.
    """
    # Get predicted mean
    mu = model.predict(X)
    
    # Get scale/dispersion
    scale = getattr(model, 'scale', 1.0)
    
    # Parametric bootstrap
    alpha = 1 - ci_level
    n_obs = len(mu)
    simulated = np.zeros((n_bootstrap, n_obs))
    
    for b in range(n_bootstrap):
        # Sample from Poisson or NB with estimated mean
        if hasattr(model.family, 'alpha'):
            # Negative Binomial
            r = 1 / model.family.alpha  # shape parameter
            p = r / (r + mu)
            simulated[b, :] = stats.nbinom.rvs(r, p)
        else:
            # Poisson (possibly overdispersed)
            simulated[b, :] = stats.poisson.rvs(mu)
    
    lower = np.percentile(simulated, 100 * alpha / 2, axis=0)
    upper = np.percentile(simulated, 100 * (1 - alpha / 2), axis=0)
    
    return lower, upper


def compute_diagnostics(model, y: np.ndarray, X: pd.DataFrame) -> Dict:
    """
    Compute model diagnostics.
    """
    # Predicted values
    mu = model.predict(X)
    
    # Residuals
    raw_resid = y - mu
    pearson_resid = (y - mu) / np.sqrt(mu)
    deviance_resid = model.resid_deviance
    
    # Dispersion
    dispersion = model.pearson_chi2 / model.df_resid
    
    # ACF of residuals (check for autocorrelation)
    from statsmodels.tsa.stattools import acf
    acf_vals = acf(pearson_resid, nlags=Config.MAX_ACF_LAG, fft=True)
    
    # Significant autocorrelation at any lag?
    n = len(pearson_resid)
    acf_ci = 1.96 / np.sqrt(n)
    sig_acf_lags = np.where(np.abs(acf_vals[1:]) > acf_ci)[0] + 1
    
    # Normality of deviance residuals (QQ)
    _, p_shapiro = stats.shapiro(deviance_resid[:min(5000, len(deviance_resid))])
    
    # Cook's distance approximation
    leverage = model.get_hat_matrix_diag() if hasattr(model, 'get_hat_matrix_diag') else None
    if leverage is not None:
        cooks_d = (pearson_resid ** 2 * leverage) / ((1 - leverage) ** 2 * X.shape[1])
        n_influential = np.sum(cooks_d > 4 / n)
    else:
        cooks_d = None
        n_influential = None
    
    return {
        'dispersion': float(dispersion),
        'acf_values': acf_vals[:11].tolist(),  # First 10 lags
        'significant_acf_lags': sig_acf_lags.tolist() if len(sig_acf_lags) < 20 else 'many',
        'shapiro_p_value': float(p_shapiro),
        'n_influential_obs': int(n_influential) if n_influential else None,
        'rmse': float(np.sqrt(np.mean(raw_resid ** 2))),
        'mae': float(np.mean(np.abs(raw_resid))),
        'mean_observed': float(np.mean(y)),
        'mean_predicted': float(np.mean(mu))
    }


def run_sensitivity_analysis(df: pd.DataFrame, 
                             include_flu: bool = False) -> Dict:
    """
    Run sensitivity analyses with different baseline periods.
    """
    results = {}
    
    for name, start_year, end_year in Config.SENSITIVITY_PERIODS:
        df_subset = df.copy()
        
        if start_year:
            df_subset = df_subset[df_subset['date'].dt.year >= start_year]
        if end_year:
            df_subset = df_subset[df_subset['date'].dt.year <= end_year]
        
        if len(df_subset) < Config.MIN_OBS_PER_REGION:
            print(f"    Sensitivity '{name}': Skipped (insufficient data)")
            continue
        
        try:
            model, X = fit_baseline_model(df_subset, include_flu=include_flu)
            expected = model.predict(X)
            excess = df_subset['deaths'].values - expected
            
            results[name] = {
                'years': f"{df_subset['date'].dt.year.min()}-{df_subset['date'].dt.year.max()}",
                'n_obs': len(df_subset),
                'total_excess': float(excess.sum()),
                'mean_daily_excess': float(excess.mean()),
                'dispersion': float(model.pearson_chi2 / model.df_resid)
            }
            print(f"    Sensitivity '{name}': Total excess = {excess.sum():,.0f}")
        except Exception as e:
            print(f"    Sensitivity '{name}': Failed ({e})")
    
    return results


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def main(include_flu: bool = False, 
         restrict_heat_season: bool = False,
         level: str = 'national'):
    """
    Main excess mortality analysis.
    
    Parameters:
    -----------
    include_flu : Adjust for influenza
    restrict_heat_season : Only analyze warm months (Oct-Mar for Brazil)
    level : 'national' or 'regional'
    """
    
    print("=" * 70)
    print("01g_v2: ROBUST EXCESS MORTALITY VALIDATION")
    print("=" * 70)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Settings: flu_adjusted={include_flu}, heat_season_only={restrict_heat_season}")
    print()
    
    # =========================================================================
    # LOAD DATA
    # =========================================================================
    print("[1] Loading data...")
    
    # Load temperature data
    df_temp = pd.read_parquet(PHASE0_RESULTS / 'era5_intermediate_daily.parquet')
    df_temp['date'] = pd.to_datetime(df_temp['date'])
    
    # Load mortality data
    df_mort = pd.read_parquet(PHASE0_RESULTS / 'mortality_regional_daily_elderly.parquet')
    df_mort['date'] = pd.to_datetime(df_mort['date'])
    
    # Load SES data (for population)
    try:
        df_ses = pd.read_csv(PHASE0_RESULTS / 'ses_intermediate_covariates.csv')
        has_pop = True
    except FileNotFoundError:
        print("  Warning: SES file not found, will not use population offset")
        has_pop = False
    
    # Load holidays
    try:
        df_holidays = pd.read_parquet(PHASE0_RESULTS / 'brazilian_holidays_daily.parquet')
        df_holidays['date'] = pd.to_datetime(df_holidays['date'])
        has_holidays = True
    except FileNotFoundError:
        print("  Warning: Holidays file not found")
        has_holidays = False
    
    # Load influenza (if needed)
    has_flu = False
    if include_flu:
        try:
            df_flu = pd.read_parquet(PHASE0_RESULTS / 'influenza_daily_by_intermediate_region.parquet')
            df_flu['date'] = pd.to_datetime(df_flu['date'])
            has_flu = True
        except FileNotFoundError:
            print("  Warning: Influenza file not found, skipping flu adjustment")
    
    # Merge datasets
    df = pd.merge(
        df_mort,
        df_temp[['date', 'region_code', 'temp_mean']],
        on=['date', 'region_code'],
        how='inner'
    )
    
    # Aggregate to national daily
    df_national = df.groupby('date').agg({
        'deaths_elderly': 'sum',
        'temp_mean': 'mean'
    }).reset_index()
    df_national.rename(columns={'deaths_elderly': 'deaths'}, inplace=True)
    
    # Add population (sum across regions for national, or use total)
    if has_pop:
        total_pop = df_ses['pop_elderly'].sum()
        df_national['pop_elderly'] = total_pop
    
    # Add holidays
    if has_holidays:
        df_national = pd.merge(df_national, df_holidays[['date', 'is_holiday']], 
                               on='date', how='left')
        df_national['is_holiday'] = df_national['is_holiday'].fillna(False).astype(int)
    
    # Add influenza (national sum)
    if has_flu:
        flu_national = df_flu.groupby('date')['srag_cases'].sum().reset_index()
        df_national = pd.merge(df_national, flu_national, on='date', how='left')
        df_national['srag_cases'] = df_national['srag_cases'].fillna(0)
    
    # Optional: restrict to heat season (Oct-Mar for Southern Hemisphere)
    if restrict_heat_season:
        heat_months = [10, 11, 12, 1, 2, 3]
        df_national = df_national[df_national['date'].dt.month.isin(heat_months)]
        print(f"  Restricted to heat season months: {heat_months}")
    
    print(f"  National daily observations: {len(df_national):,}")
    print(f"  Total deaths: {df_national['deaths'].sum():,}")
    print(f"  Date range: {df_national['date'].min().date()} to {df_national['date'].max().date()}")
    if has_pop:
        print(f"  Population (elderly): {total_pop:,}")
    
    # Temperature percentiles
    temp_percentiles = {
        'p01': df_national['temp_mean'].quantile(0.01),
        'p025': df_national['temp_mean'].quantile(0.025),
        'p05': df_national['temp_mean'].quantile(0.05),
        'p25': df_national['temp_mean'].quantile(0.25),
        'p50': df_national['temp_mean'].quantile(0.50),
        'p75': df_national['temp_mean'].quantile(0.75),
        'p95': df_national['temp_mean'].quantile(0.95),
        'p975': df_national['temp_mean'].quantile(0.975),
        'p99': df_national['temp_mean'].quantile(0.99)
    }
    
    print(f"\n  Temperature percentiles:")
    for k, v in temp_percentiles.items():
        print(f"    {k}: {v:.1f}°C")
    
    # =========================================================================
    # FIT BASELINE MODEL
    # =========================================================================
    print("\n[2] Fitting baseline model (flexible splines + controls)...")
    
    model, X = fit_baseline_model(
        df_national, 
        include_flu=has_flu and include_flu,
        include_holidays=has_holidays,
        use_nb=True
    )
    
    # Predicted (expected) deaths
    df_national['expected'] = model.predict(X)
    df_national['excess'] = df_national['deaths'] - df_national['expected']
    
    print(f"  Model family: {model.family.__class__.__name__}")
    print(f"  Model dispersion: {model.scale:.2f}")
    print(f"  Mean observed: {df_national['deaths'].mean():.0f}")
    print(f"  Mean expected: {df_national['expected'].mean():.0f}")
    
    # =========================================================================
    # COMPUTE PREDICTION INTERVALS
    # =========================================================================
    print(f"\n[3] Computing prediction intervals ({Config.N_BOOTSTRAP} bootstrap samples)...")
    
    lower, upper = compute_prediction_intervals(model, X, 
                                                 n_bootstrap=Config.N_BOOTSTRAP,
                                                 ci_level=Config.CI_LEVEL)
    df_national['expected_lower'] = lower
    df_national['expected_upper'] = upper
    df_national['excess_lower'] = df_national['deaths'] - upper
    df_national['excess_upper'] = df_national['deaths'] - lower
    
    print(f"  95% PI coverage computed")
    
    # =========================================================================
    # DIAGNOSTICS
    # =========================================================================
    print("\n[4] Computing model diagnostics...")
    
    diagnostics = compute_diagnostics(model, df_national['deaths'].values, X)
    
    print(f"  Dispersion: {diagnostics['dispersion']:.2f}")
    print(f"  RMSE: {diagnostics['rmse']:.1f}")
    print(f"  MAE: {diagnostics['mae']:.1f}")
    if diagnostics['significant_acf_lags']:
        if diagnostics['significant_acf_lags'] == 'many':
            print(f"  ⚠️ Significant residual autocorrelation at many lags")
        else:
            print(f"  ⚠️ Significant ACF at lags: {diagnostics['significant_acf_lags'][:5]}...")
    else:
        print(f"  ✓ No significant residual autocorrelation")
    
    # =========================================================================
    # CALCULATE EXCESS BY TEMPERATURE CATEGORY
    # =========================================================================
    print("\n[5] Calculating excess mortality by temperature category...")
    
    # Define masks using P2.5/P97.5 (PRIMARY) and P1/P99 (COMPARISON)
    thresholds = {
        'primary': {
            'heat': ('p975', temp_percentiles['p975']),
            'cold': ('p025', temp_percentiles['p025'])
        },
        'comparison': {
            'heat': ('p99', temp_percentiles['p99']),
            'cold': ('p01', temp_percentiles['p01'])
        },
        'traditional': {
            'heat': ('p75', temp_percentiles['p75']),
            'cold': ('p25', temp_percentiles['p25'])
        }
    }
    
    excess_results = {}
    
    for thresh_name, thresh_vals in thresholds.items():
        heat_thresh = thresh_vals['heat'][1]
        cold_thresh = thresh_vals['cold'][1]
        
        heat_mask = df_national['temp_mean'] > heat_thresh
        cold_mask = df_national['temp_mean'] < cold_thresh
        
        heat_excess = df_national.loc[heat_mask, 'excess'].sum()
        cold_excess = df_national.loc[cold_mask, 'excess'].sum()
        
        # With uncertainty
        heat_excess_lower = df_national.loc[heat_mask, 'excess_lower'].sum()
        heat_excess_upper = df_national.loc[heat_mask, 'excess_upper'].sum()
        cold_excess_lower = df_national.loc[cold_mask, 'excess_lower'].sum()
        cold_excess_upper = df_national.loc[cold_mask, 'excess_upper'].sum()
        
        n_heat_days = heat_mask.sum()
        n_cold_days = cold_mask.sum()
        
        excess_results[thresh_name] = {
            'heat': {
                'threshold': float(heat_thresh),
                'threshold_name': thresh_vals['heat'][0],
                'n_days': int(n_heat_days),
                'excess': float(heat_excess),
                'excess_95ci': [float(heat_excess_lower), float(heat_excess_upper)],
                'mean_excess_per_day': float(heat_excess / n_heat_days) if n_heat_days > 0 else 0
            },
            'cold': {
                'threshold': float(cold_thresh),
                'threshold_name': thresh_vals['cold'][0],
                'n_days': int(n_cold_days),
                'excess': float(cold_excess),
                'excess_95ci': [float(cold_excess_lower), float(cold_excess_upper)],
                'mean_excess_per_day': float(cold_excess / n_cold_days) if n_cold_days > 0 else 0
            }
        }
        
        print(f"\n  {thresh_name.upper()} thresholds:")
        print(f"    Heat (>{heat_thresh:.1f}°C): {n_heat_days} days, "
              f"excess = {heat_excess:,.0f} [{heat_excess_lower:,.0f}, {heat_excess_upper:,.0f}]")
        print(f"    Cold (<{cold_thresh:.1f}°C): {n_cold_days} days, "
              f"excess = {cold_excess:,.0f} [{cold_excess_lower:,.0f}, {cold_excess_upper:,.0f}]")
    
    # =========================================================================
    # COMPARE WITH DLNM BURDEN
    # =========================================================================
    print("\n[6] Comparing with DLNM attributable burden...")
    
    comparison = {}
    try:
        # Prefer v2 national burden summary (both intermediate & immediate)
        burden_file = PHASE1_RESULTS / 'burden_v2_national_summary.json'
        if not burden_file.exists():
            # Backwards compatibility: older naming
            burden_file = PHASE1_RESULTS / 'attributable_burden_national.json'

        with open(burden_file, 'r') as f:
            dlnm_burden = json.load(f)

        # Use PRIMARY thresholds for comparison
        em_heat = excess_results['primary']['heat']['excess']
        em_cold = excess_results['primary']['cold']['excess']

        level_comparisons = {}

        # Helper to extract heat/cold AN for a given level block
        def _extract_level_an(level_data):
            if not level_data:
                return None, None
            # v2 naming (preferred)
            heat_an = level_data.get('heat_an_97_5')
            cold_an = level_data.get('cold_an_2_5')
            # Fallback to generic names if present
            if heat_an is None:
                heat_an = level_data.get('heat_an')
            if cold_an is None:
                cold_an = level_data.get('cold_an')
            if heat_an is None or cold_an is None:
                return None, None
            return float(heat_an), float(cold_an)

        # Try both intermediate and immediate levels when available
        for level_name in ['intermediate', 'immediate']:
            level_data = dlnm_burden.get(level_name, {}) if isinstance(dlnm_burden, dict) else {}
            heat_an, cold_an = _extract_level_an(level_data)
            if heat_an is None:
                continue

            heat_ratio = em_heat / heat_an if heat_an > 0 else np.nan
            cold_ratio = em_cold / cold_an if cold_an > 0 else np.nan

            level_comparisons[level_name] = {
                'dlnm_heat_an': heat_an,
                'dlnm_cold_an': cold_an,
                'excess_heat': float(em_heat),
                'excess_cold': float(em_cold),
                'ratio_heat': float(heat_ratio),
                'ratio_cold': float(cold_ratio)
            }

        # Backwards-compatible flat structure (no levels key)
        if not level_comparisons and isinstance(dlnm_burden, dict):
            heat_an, cold_an = _extract_level_an(dlnm_burden)
            if heat_an is not None:
                heat_ratio = em_heat / heat_an if heat_an > 0 else np.nan
                cold_ratio = em_cold / cold_an if cold_an > 0 else np.nan
                level_comparisons['national'] = {
                    'dlnm_heat_an': heat_an,
                    'dlnm_cold_an': cold_an,
                    'excess_heat': float(em_heat),
                    'excess_cold': float(em_cold),
                    'ratio_heat': float(heat_ratio),
                    'ratio_cold': float(cold_ratio)
                }

        if level_comparisons:
            comparison = {
                'levels': level_comparisons,
                'threshold_used': 'primary (P2.5/P97.5)'
            }

            for lvl, vals in level_comparisons.items():
                label = f"{lvl.upper()}" if lvl in ['intermediate', 'immediate'] else lvl
                print(f"\n  DLNM Attributable (P2.5/P97.5) - {label}:")
                print(f"    Heat AN: {vals['dlnm_heat_an']:,.0f}")
                print(f"    Cold AN: {vals['dlnm_cold_an']:,.0f}")
                print(f"\n  Excess Mortality (P2.5/P97.5):")
                print(f"    Heat: {vals['excess_heat']:,.0f}")
                print(f"    Cold: {vals['excess_cold']:,.0f}")
                print(f"\n  Ratio (Excess / DLNM):")
                print(f"    Heat: {vals['ratio_heat']:.2f}")
                print(f"    Cold: {vals['ratio_cold']:.2f}")
        else:
            print("  Could not extract DLNM values from burden file")
            
    except FileNotFoundError:
        print("  DLNM burden file not found, skipping comparison")
    
    # =========================================================================
    # SENSITIVITY ANALYSIS
    # =========================================================================
    print("\n[7] Running sensitivity analyses (baseline period)...")
    
    sensitivity = run_sensitivity_analysis(df_national, include_flu=has_flu and include_flu)
    
    # =========================================================================
    # ANNUAL SUMMARY
    # =========================================================================
    print("\n[8] Computing annual excess...")
    
    heat_thresh = temp_percentiles['p975']
    cold_thresh = temp_percentiles['p025']
    
    annual_summary = df_national.groupby(df_national['date'].dt.year).apply(
        lambda x: pd.Series({
            'total_deaths': x['deaths'].sum(),
            'expected': x['expected'].sum(),
            'total_excess': x['excess'].sum(),
            'heat_excess': x.loc[x['temp_mean'] > heat_thresh, 'excess'].sum(),
            'cold_excess': x.loc[x['temp_mean'] < cold_thresh, 'excess'].sum(),
            'n_heat_days': (x['temp_mean'] > heat_thresh).sum(),
            'n_cold_days': (x['temp_mean'] < cold_thresh).sum()
        })
    ).reset_index()
    annual_summary.rename(columns={'date': 'year'}, inplace=True)
    
    print("\n  Annual Summary (P2.5/P97.5):")
    print(annual_summary[['year', 'total_deaths', 'heat_excess', 'cold_excess']].to_string(index=False))
    
    # =========================================================================
    # SAVE RESULTS
    # =========================================================================
    print("\n[9] Saving results...")
    
    def convert_for_json(obj):
        if isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_for_json(i) for i in obj]
        elif isinstance(obj, pd.DataFrame):
            return obj.to_dict('records')
        return obj
    
    results = {
        'method': 'Counterfactual excess mortality with flexible baseline model',
        'version': 'v2',
        'settings': {
            'include_flu': include_flu,
            'restrict_heat_season': restrict_heat_season,
            'seasonal_df': Config.SEASONAL_DF,
            'n_bootstrap': Config.N_BOOTSTRAP,
            'ci_level': Config.CI_LEVEL
        },
        'model': {
            'family': model.family.__class__.__name__,
            'dispersion': float(model.scale),
            'n_parameters': X.shape[1],
            'n_observations': len(df_national),
            'has_population_offset': has_pop,
            'has_holiday_control': has_holidays,
            'has_flu_adjustment': has_flu and include_flu
        },
        'diagnostics': diagnostics,
        'temperature_percentiles': {k: float(v) for k, v in temp_percentiles.items()},
        'excess_mortality': excess_results,
        'comparison_with_dlnm': comparison,
        'sensitivity_analyses': sensitivity,
        'annual_summary': annual_summary.to_dict('records'),
        'timestamp': datetime.now().isoformat()
    }
    
    output_suffix = '_flu' if include_flu else ''
    output_suffix += '_heat_season' if restrict_heat_season else ''
    
    with open(OUTPUT_DIR / f'excess_mortality_v2{output_suffix}.json', 'w') as f:
        json.dump(convert_for_json(results), f, indent=2)
    
    annual_summary.to_csv(OUTPUT_DIR / f'excess_mortality_v2_annual{output_suffix}.csv', index=False)
    
    print(f"  Saved: excess_mortality_v2{output_suffix}.json")
    print(f"  Saved: excess_mortality_v2_annual{output_suffix}.csv")
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY (v2)")
    print("=" * 70)
    
    em_heat = excess_results['primary']['heat']['excess']
    em_cold = excess_results['primary']['cold']['excess']
    em_heat_ci = excess_results['primary']['heat']['excess_95ci']
    em_cold_ci = excess_results['primary']['cold']['excess_95ci']
    
    print(f"""
Excess Mortality Estimates (P2.5/P97.5 thresholds):

1. HEAT (>{temp_percentiles['p975']:.1f}°C):
   - Excess: {em_heat:,.0f} deaths [95% PI: {em_heat_ci[0]:,.0f} to {em_heat_ci[1]:,.0f}]
   - Days: {excess_results['primary']['heat']['n_days']:,}

2. COLD (<{temp_percentiles['p025']:.1f}°C):
   - Excess: {em_cold:,.0f} deaths [95% PI: {em_cold_ci[0]:,.0f} to {em_cold_ci[1]:,.0f}]
   - Days: {excess_results['primary']['cold']['n_days']:,}

3. MODEL QUALITY:
   - Dispersion: {diagnostics['dispersion']:.2f} {'✓' if diagnostics['dispersion'] < 2 else '⚠️'}
   - ACF issues: {'No ✓' if not diagnostics['significant_acf_lags'] else 'Yes ⚠️'}
   - Population offset: {'Yes ✓' if has_pop else 'No ⚠️'}
   - Holiday control: {'Yes ✓' if has_holidays else 'No'}
   - Flu adjustment: {'Yes' if has_flu and include_flu else 'No'}

4. COLD > HEAT PATTERN: {'✓ CONFIRMED' if em_cold > em_heat else '✗ NOT CONFIRMED'}
""")

    if comparison:
        levels = comparison.get('levels', {})
        print("\n5. COMPARISON WITH DLNM:")
        for lvl, vals in levels.items():
            label = f"{lvl.upper()}" if lvl in ['intermediate', 'immediate'] else lvl
            print(f"   - {label} Heat ratio (Excess/DLNM): {vals['ratio_heat']:.2f}")
            print(f"     {label} Cold ratio (Excess/DLNM): {vals['ratio_cold']:.2f}")
        print("""
   Interpretation:
   - Ratio > 1: Excess finds MORE deaths (may include non-temperature causes)
   - Ratio < 1: Excess finds FEWER deaths (expected - DLNM captures lag effects)
   - Ratio ≈ 1: Methods agree
""")
    
    print("=" * 70)
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)
    
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Excess Mortality Validation v2')
    parser.add_argument('--flu', action='store_true', 
                        help='Include influenza adjustment (single run)')
    parser.add_argument('--heat-season', action='store_true',
                        help='Restrict to heat season months (Oct-Mar, single run)')
    parser.add_argument('--single', action='store_true',
                        help='Run single model (default runs all variants)')
    
    args = parser.parse_args()
    
    if args.single:
        # Run single variant based on flags
        main(include_flu=args.flu, restrict_heat_season=args.heat_season)
    else:
        # Default: Run all variants
        print("Running all model variants (use --single for one model)...\n")
        main(include_flu=False, restrict_heat_season=False)
        print("\n" + "=" * 70 + "\n")
        main(include_flu=True, restrict_heat_season=False)
        print("\n" + "=" * 70 + "\n")
        main(include_flu=False, restrict_heat_season=True)
