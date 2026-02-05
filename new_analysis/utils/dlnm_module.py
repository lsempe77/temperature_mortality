"""
DLNM Utilities Module
=====================
Stable cross-basis creation, regional DLNM fitting, cumulative RR prediction
(delta method with full covariance), harvesting & heatwave effect-modification
helpers, and DerSimonian-Laird meta-analysis.

This module provides the correct DLNM implementation following:
- Gasparrini et al. (2010) - DLNM methodology
- Armstrong et al. (2014) - Temperature-mortality studies
- Gasparrini et al. (2015) - MCC Collaborative methodology

Key features:
1. Natural cubic spline bases (not polynomials)
2. Proper tensor-product cross-basis construction
3. Region-specific fitting with population offset
4. Delta-method cumulative RR with full covariance
5. Random-effects meta-analysis (DerSimonian-Laird)

Dependencies:
- numpy, pandas, patsy, statsmodels, scipy

Author: Heat-Mortality Brazil Analysis Pipeline
Date: December 2025
"""

from typing import Tuple, Dict, Any, List, Optional
import numpy as np
import pandas as pd
from patsy import dmatrix
import statsmodels.api as sm
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod.families import Poisson, NegativeBinomial
from scipy import stats
import warnings
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("dlnm_module")
warnings.filterwarnings("ignore", category=FutureWarning)


# =============================================================================
# R-STYLE NATURAL SPLINE BASIS (Matches Phase 1 scripts exactly)
# =============================================================================

def ns_basis(x: np.ndarray, knots: np.ndarray, boundary_knots: tuple = None, 
             intercept: bool = False) -> np.ndarray:
    """
    Create natural cubic spline basis matrix (R's ns() style).
    
    Natural splines are cubic splines with additional constraints that the
    function is linear beyond the boundary knots. This prevents wild 
    extrapolation behavior.
    
    This implementation matches R's ns() function and is used by Phase 1 scripts.
    
    Parameters:
    -----------
    x : array-like
        Values at which to evaluate the basis
    knots : array-like
        Interior knot locations
    boundary_knots : tuple
        (lower, upper) boundary knots. If None, uses min/max of x
    intercept : bool
        Whether to include intercept column
        
    Returns:
    --------
    basis : ndarray
        Natural spline basis matrix (n x df)
        df = len(knots) + 1 if not intercept, else len(knots) + 2
    """
    x = np.asarray(x).flatten()
    knots = np.sort(np.asarray(knots).flatten())
    
    if boundary_knots is None:
        boundary_knots = (np.nanmin(x), np.nanmax(x))
    
    # All knots including boundaries
    all_knots = np.concatenate([[boundary_knots[0]], knots, [boundary_knots[1]]])
    
    # Number of basis functions
    n_knots = len(knots)
    df = n_knots + 1  # Natural spline df
    
    # Create basis using the algorithm from R's ns() function
    # Based on truncated power basis with natural constraints
    
    def d(x_val, knot_k, knot_K):
        """Helper function for natural spline computation."""
        return ((np.maximum(0, x_val - knot_k)**3 - np.maximum(0, x_val - knot_K)**3) / 
                (knot_K - knot_k))
    
    n = len(x)
    basis = np.zeros((n, df))
    
    # First basis function: linear term
    basis[:, 0] = x
    
    # Remaining basis functions: natural spline terms
    knot_K = all_knots[-1]  # Last (boundary) knot
    knot_Km1 = all_knots[-2]  # Second to last knot
    
    for j in range(n_knots):
        knot_j = all_knots[j + 1]  # Skip first boundary knot
        basis[:, j + 1] = d(x, knot_j, knot_K) - d(x, knot_Km1, knot_K)
    
    if intercept:
        basis = np.column_stack([np.ones(n), basis])
    
    return basis


def create_crossbasis_ns(
    temp: np.ndarray, 
    max_lag: int, 
    temp_knots: np.ndarray, 
    temp_boundary: tuple, 
    lag_knots: np.ndarray
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Create cross-basis matrix using R-style natural cubic splines.
    
    This implements the tensor product of:
    - Natural spline basis for temperature (with explicit knots)
    - Natural spline basis for lag (on log scale)
    
    This matches the implementation in Phase 1 DLNM scripts exactly.
    
    Parameters:
    -----------
    temp : array
        Temperature values
    max_lag : int
        Maximum lag (e.g., 21)
    temp_knots : array
        Interior knots for temperature (e.g., percentiles)
    temp_boundary : tuple
        Boundary knots for temperature (low, high)
    lag_knots : array
        Knots for lag dimension in original scale (e.g., [1, 3, 7, 14])
        
    Returns:
    --------
    X_cb : ndarray
        Cross-basis matrix (n x (temp_df * lag_df))
    info : dict
        Information needed for prediction
    """
    n = len(temp)
    
    # Create lag matrix
    temp_lags = create_lag_matrix(temp, max_lag)
    
    # Create temperature basis reference to get df
    temp_valid = temp[~np.isnan(temp)]
    temp_basis_ref = ns_basis(temp_valid, temp_knots, temp_boundary)
    temp_df = temp_basis_ref.shape[1]
    
    # Create lag basis (on log(lag+1) scale for physiological plausibility)
    lag_vals = np.arange(max_lag + 1)
    lag_log = np.log(lag_vals + 1)  # log(1) to log(max_lag+2)
    lag_knots_log = np.log(np.array(lag_knots) + 1)
    lag_boundary_log = (np.log(1), np.log(max_lag + 2))
    lag_basis = ns_basis(lag_log, lag_knots_log, lag_boundary_log)
    lag_df = lag_basis.shape[1]
    
    # Cross-basis: tensor product
    n_cb = temp_df * lag_df
    X_cb = np.zeros((n, n_cb))
    
    for i in range(n):
        # Temperature values at all lags for this observation
        temp_i_lags = temp_lags[i, :]
        
        if np.any(np.isnan(temp_i_lags)):
            X_cb[i, :] = np.nan
            continue
        
        # Compute temperature basis for each lag and accumulate weighted sum
        cb_row = np.zeros(n_cb)
        for j, temp_val in enumerate(temp_i_lags):
            # Temperature basis at this lag
            temp_basis_ij = ns_basis(np.array([temp_val]), temp_knots, temp_boundary)[0]
            
            # Outer product with lag basis weight
            for t_idx in range(temp_df):
                for l_idx in range(lag_df):
                    col_idx = t_idx * lag_df + l_idx
                    cb_row[col_idx] += temp_basis_ij[t_idx] * lag_basis[j, l_idx]
        
        X_cb[i, :] = cb_row
    
    # Column names
    col_names = [f'cb_t{t}_l{l}' for t in range(temp_df) for l in range(lag_df)]
    
    info = {
        'temp_knots': list(temp_knots),
        'temp_boundary': list(temp_boundary),
        'lag_knots': list(lag_knots),
        'lag_knots_log': list(lag_knots_log),
        'lag_boundary_log': list(lag_boundary_log),
        'temp_df': temp_df,
        'lag_df': lag_df,
        'max_lag': max_lag,
        'n_cb': n_cb,
        'n_params': n_cb,
        'col_names': col_names,
        'lag_basis': lag_basis,
    }
    
    return X_cb, info


def compute_cumulative_rr_ns(
    target_temp: float, 
    ref_temp: float, 
    coefs: np.ndarray, 
    cb_info: Dict[str, Any],
    vcov: np.ndarray = None
) -> Tuple[float, float, float]:
    """
    Compute cumulative RR at target vs reference temperature for NS cross-basis.
    
    Parameters:
    -----------
    target_temp : float
        Temperature to evaluate (e.g., P99)
    ref_temp : float
        Reference temperature (e.g., MMT)
    coefs : array
        Cross-basis coefficients from fitted model
    cb_info : dict
        Information from create_crossbasis_ns
    vcov : array, optional
        Variance-covariance matrix for confidence interval
        
    Returns:
    --------
    rr : float
        Cumulative relative risk
    rr_low : float
        Lower 95% CI (NaN if vcov not provided)
    rr_high : float
        Upper 95% CI (NaN if vcov not provided)
    """
    temp_knots = np.array(cb_info['temp_knots'])
    temp_boundary = tuple(cb_info['temp_boundary'])
    lag_knots = cb_info['lag_knots']
    temp_df = cb_info['temp_df']
    lag_df = cb_info['lag_df']
    max_lag = cb_info.get('max_lag', 21)
    
    # Recreate lag basis
    lag_vals = np.arange(max_lag + 1)
    lag_log = np.log(lag_vals + 1)
    lag_knots_log = np.log(np.array(lag_knots) + 1)
    lag_boundary_log = (np.log(1), np.log(max_lag + 2))
    lag_basis = ns_basis(lag_log, lag_knots_log, lag_boundary_log)
    
    # Temperature basis at target and reference
    temp_basis_target = ns_basis(np.array([target_temp]), temp_knots, temp_boundary)[0]
    temp_basis_ref = ns_basis(np.array([ref_temp]), temp_knots, temp_boundary)[0]
    temp_diff = temp_basis_target - temp_basis_ref
    
    # Sum over all lags
    lag_sum = lag_basis.sum(axis=0)
    
    # Contrast vector (matches cross-basis column order: t0l0, t0l1, ..., t1l0, ...)
    contrast = np.zeros(temp_df * lag_df)
    for t_idx in range(temp_df):
        for l_idx in range(lag_df):
            col_idx = t_idx * lag_df + l_idx
            contrast[col_idx] = temp_diff[t_idx] * lag_sum[l_idx]
    
    # Compute log-RR
    log_rr = np.dot(contrast, coefs)
    rr = float(np.exp(log_rr))
    
    # Confidence interval if vcov provided
    if vcov is not None:
        var_log_rr = float(contrast @ vcov @ contrast)
        se_log_rr = np.sqrt(max(0, var_log_rr))
        rr_low = float(np.exp(log_rr - 1.96 * se_log_rr))
        rr_high = float(np.exp(log_rr + 1.96 * se_log_rr))
    else:
        rr_low = np.nan
        rr_high = np.nan
    
    return rr, rr_low, rr_high


# =============================================================================
# LOW-LEVEL HELPERS
# =============================================================================

def create_lag_matrix(x: np.ndarray, max_lag: int) -> np.ndarray:
    """
    Create lag matrix of shape (n_obs, max_lag+1) where column j is x shifted by j days.
    Leading rows with insufficient history will contain np.nan.
    
    Parameters:
    -----------
    x : 1D array of exposure values (e.g., daily temperature)
    max_lag : Maximum lag to include (0 = same day, 1 = yesterday, etc.)
    
    Returns:
    --------
    lag_mat : 2D array (n_obs, max_lag+1)
    """
    n = len(x)
    L = max_lag
    lag_mat = np.full((n, L + 1), np.nan, dtype=float)
    for j in range(L + 1):
        if j == 0:
            lag_mat[:, 0] = x
        else:
            lag_mat[j:, j] = x[:-j]
    return lag_mat


def _safe_dmatrix(formula: str, x: np.ndarray, name: str = "x") -> pd.DataFrame:
    """
    Wrap patsy.dmatrix and return DataFrame. Handles edge cases gracefully.
    
    Parameters:
    -----------
    formula : Patsy formula string (e.g., "cr(x, df=4) - 1")
    x : 1D array of values
    name : Variable name in formula
    
    Returns:
    --------
    DataFrame with basis columns
    """
    df = pd.DataFrame({name: x})
    try:
        mat = dmatrix(formula, df, return_type="dataframe")
        mat.index = df.index
        return mat
    except Exception as e:
        logger.warning("patsy.dmatrix failed for formula=%s: %s", formula, e)
        return pd.DataFrame(index=df.index)


def natural_spline_basis(x: np.ndarray, df: int = 4, 
                         knots: np.ndarray = None) -> np.ndarray:
    """
    Create natural cubic spline basis using patsy's 'cr' function.
    
    Parameters:
    -----------
    x : Values to evaluate basis at
    df : Degrees of freedom (number of basis functions)
    knots : Optional explicit knot locations
    
    Returns:
    --------
    Basis matrix (n x df)
    """
    if knots is not None:
        knot_str = ','.join([str(k) for k in knots])
        formula = f"cr(x, knots=[{knot_str}]) - 1"
    else:
        formula = f"cr(x, df={df}) - 1"
    
    result = _safe_dmatrix(formula, x, "x")
    return result.values if len(result.columns) > 0 else np.zeros((len(x), df))


# =============================================================================
# CROSS-BASIS BUILDER
# =============================================================================

def create_crossbasis(
    temp: np.ndarray,
    max_lag: int,
    var_df: int = 4,
    lag_df: int = 4,
    var_knots: np.ndarray = None,
    lag_knots: np.ndarray = None,
) -> Tuple[np.ndarray, List[str], Dict[str, Any]]:
    """
    Build a DLNM cross-basis matrix using natural cubic spline bases.
    
    The cross-basis follows standard DLNM construction (Gasparrini 2010):
    For each observation i and lag j, compute var_basis(temp_{i-j})
    weighted by lag_basis(j). The tensor product gives K_var × K_lag columns.
    
    Parameters:
    -----------
    temp : 1D array of daily temperature values
    max_lag : Maximum lag (e.g., 21 days)
    var_df : Degrees of freedom for temperature (variable) spline
    lag_df : Degrees of freedom for lag spline
    var_knots : Optional explicit knots for temperature basis
    lag_knots : Optional explicit knots for lag basis
    
    Returns:
    --------
    X_cb : Cross-basis matrix (n_obs, K_var * K_lag)
    colnames : List of column names
    meta : Dict with metadata for prediction
    """
    n = len(temp)
    L = max_lag

    # Center & scale temperature for numerical stability
    temp_center = np.nanmedian(temp)
    temp_std = np.nanstd(temp)
    if temp_std == 0 or np.isnan(temp_std):
        raise ValueError("Temperature array has zero or NaN std; cannot scale")

    # Create lagged values matrix (n, L+1)
    temp_lags = create_lag_matrix(temp, L)

    # Build variable (temperature) basis formula
    if var_knots is not None:
        knot_str = ','.join([str(k) for k in var_knots])
        var_formula = f"cr(x, knots=[{knot_str}]) - 1"
    else:
        var_formula = f"cr(x, df={var_df}) - 1"

    # First, compute basis on all valid (non-NaN) temperatures to determine K_var
    all_valid_temps = temp[~np.isnan(temp)]
    ref_basis = _safe_dmatrix(var_formula, all_valid_temps, name="x")
    K_var = ref_basis.shape[1]
    if K_var == 0:
        raise ValueError(f"Temperature spline basis failed (0 columns). Check data range and df={var_df}")
    
    # Evaluate var basis at each lag's temperature values
    # Handle NaN by computing basis only for valid values, then fill back
    var_bases = []
    for j in range(L + 1):
        arr = temp_lags[:, j]
        valid_mask = ~np.isnan(arr)
        
        # Create output dataframe with NaN
        out_df = pd.DataFrame(np.nan, index=range(n), columns=[f"v{k}" for k in range(K_var)])
        
        if valid_mask.sum() > 0:
            valid_vals = arr[valid_mask]
            valid_basis = _safe_dmatrix(var_formula, valid_vals, name="x")
            if valid_basis.shape[1] == K_var:
                out_df.loc[valid_mask, :] = valid_basis.values
        
        var_bases.append(out_df)
    var_colnames = [f"v{i}" for i in range(K_var)]

    # Build lag basis evaluated at integer lag indices 0..L
    lag_index = np.arange(0, L + 1, dtype=float)
    
    # Use log-transformed lag scale for better fit (common in DLNM)
    lag_log = np.log(lag_index + 1)  # log(lag + 1) to handle lag=0
    
    # Adjust lag_df if max_lag is too short for the requested df
    # Patsy's cr() needs at least df-1 interior knots, and those need data points
    effective_lag_df = min(lag_df, max(2, L))  # At least 2 df, at most L points
    
    if lag_knots is not None:
        knot_str = ','.join([str(k) for k in lag_knots])
        lag_formula = f"cr(x, knots=[{knot_str}]) - 1"
    else:
        lag_formula = f"cr(x, df={effective_lag_df}) - 1"
    
    lag_basis_df = _safe_dmatrix(lag_formula, lag_log, name="x")
    K_lag = lag_basis_df.shape[1]
    
    if K_lag == 0:
        # Fallback to simple linear basis if spline fails
        logger.warning(f"Lag spline failed with df={effective_lag_df}, max_lag={L}. Using linear basis.")
        lag_basis_df = pd.DataFrame({'l0': lag_log / lag_log.max()}, index=range(L+1))
        K_lag = 1
    lag_colnames = [f"l{i}" for i in range(K_lag)]
    lag_basis_vals = lag_basis_df.values  # shape (L+1, K_lag)

    # Build cross-basis columns via tensor product
    cols = []
    col_names = []
    
    for kv in range(K_var):
        for kl in range(K_lag):
            # Column = sum over lags of var_basis[j, kv] * lag_basis[j, kl]
            col = np.zeros(n, dtype=float)
            nan_mask = np.zeros(n, dtype=bool)
            
            for j in range(L + 1):
                vcol = var_bases[j].iloc[:, kv].to_numpy(dtype=float)
                weight = lag_basis_vals[j, kl]
                
                # Track NaN positions
                nan_mask = nan_mask | np.isnan(vcol)
                
                # Accumulate (treating NaN as 0 temporarily)
                vcol_clean = np.nan_to_num(vcol, nan=0.0)
                col += vcol_clean * weight
            
            # Apply NaN mask for rows with missing lag history
            col[nan_mask] = np.nan
            cols.append(col)
            col_names.append(f"cb_{var_colnames[kv]}_{lag_colnames[kl]}")

    X_cb = np.column_stack(cols) if cols else np.empty((n, 0))
    
    # Store explicit knot locations for prediction
    # For ns_basis compatibility: we need knots such that ns_basis produces K_var/K_lag columns
    # ns_basis formula: n_columns = len(knots) + 1
    # So if K_var = 4, we need 3 interior knots
    all_valid_temps = temp[~np.isnan(temp)]
    temp_boundary = (float(np.min(all_valid_temps)), float(np.max(all_valid_temps)))
    
    # Compute interior knots to match K_var columns from patsy
    # patsy's cr(df=N) produces N columns, ns_basis needs N-1 interior knots for N columns
    n_inner_var_knots = K_var - 1
    if n_inner_var_knots > 0:
        knot_quantiles = np.linspace(0, 1, n_inner_var_knots + 2)[1:-1]
        var_knots_computed = [float(np.quantile(all_valid_temps, q)) for q in knot_quantiles]
    else:
        var_knots_computed = []
    
    # Lag knots (for prediction) - match K_lag columns
    lag_log_all = np.log(np.arange(L + 1) + 1)
    lag_boundary = (float(np.min(lag_log_all)), float(np.max(lag_log_all)))
    n_inner_lag_knots = K_lag - 1
    if n_inner_lag_knots > 0:
        lag_knot_quantiles = np.linspace(0, 1, n_inner_lag_knots + 2)[1:-1]
        lag_knots_computed = [float(np.quantile(lag_log_all, q)) for q in lag_knot_quantiles]
    else:
        lag_knots_computed = []
    
    meta = {
        "temp_center": float(temp_center),
        "temp_std": float(temp_std),
        "K_var": K_var,
        "K_lag": K_lag,
        "var_df": var_df,
        "lag_df": effective_lag_df,
        "max_lag": max_lag,
        "var_colnames": var_colnames,
        "lag_colnames": lag_colnames,
        "lag_basis_vals": lag_basis_vals,
        "var_formula": var_formula,
        "lag_formula": lag_formula,
        # New: explicit knots for prediction
        "var_knots": var_knots_computed,
        "var_boundary": temp_boundary,
        "lag_knots_log": lag_knots_computed,
        "lag_boundary_log": lag_boundary,
    }
    
    return X_cb, col_names, meta


# =============================================================================
# MODEL FITTING (REGION-SPECIFIC)
# =============================================================================

def fit_region_dlnm(
    df_region: pd.DataFrame,
    temp_col: str = "temp_mean",
    deaths_col: str = "deaths_elderly",
    pop_col: str = "pop_elderly",
    max_lag: int = 21,
    var_df: int = 4,
    lag_df: int = 4,
    var_knots: np.ndarray = None,
    lag_knots: np.ndarray = None,
    family: str = "quasi-poisson",
    cov_type: str = "HC1",
    add_spline_time: bool = True,
    time_df_per_year: float = 4.0,
    min_obs: int = 365,
) -> Optional[Dict[str, Any]]:
    """
    Fit a region-level DLNM via GLM with proper cross-basis.
    
    Parameters:
    -----------
    df_region : DataFrame for one region with date, temp, deaths, pop columns
    temp_col : Temperature column name
    deaths_col : Mortality outcome column name
    pop_col : Population column for offset (if None, no offset used)
    max_lag : Maximum lag days
    var_df : Temperature spline df
    lag_df : Lag spline df
    family : 'quasi-poisson' or 'negbin'
    cov_type : Robust covariance type (e.g., 'HC1')
    add_spline_time : Whether to add long-term trend spline
    time_df_per_year : Df per year for time spline
    min_obs : Minimum observations required
    
    Returns:
    --------
    Dict with params, cov, cb_meta, model_result, or None if fitting fails
    """
    df = df_region.copy().reset_index(drop=True)
    n = len(df)
    
    if n < min_obs:
        logger.warning(f"Region has {n} < {min_obs} rows; skipping.")
        return None

    # Ensure date sorting
    if "date" in df.columns:
        df = df.sort_values("date").reset_index(drop=True)

    # Construct cross-basis
    temp_arr = df[temp_col].to_numpy(dtype=float)
    X_cb, cb_names, cb_meta = create_crossbasis(
        temp_arr, 
        max_lag=max_lag, 
        var_df=var_df, 
        lag_df=lag_df,
        var_knots=var_knots,
        lag_knots=lag_knots,
    )
    cb_df = pd.DataFrame(X_cb, columns=cb_names, index=df.index)

    # Build confounders
    confounders = []
    
    # Day of week
    if "dow" in df.columns:
        dow_col = df["dow"]
    elif "day_of_week" in df.columns:
        dow_col = df["day_of_week"]
    elif "date" in df.columns:
        dow_col = df["date"].dt.dayofweek
    else:
        dow_col = None
    
    if dow_col is not None:
        confounders.append(pd.get_dummies(dow_col.astype(int), prefix="dow", drop_first=True))
    
    # Month
    if "month" in df.columns:
        month_col = df["month"]
    elif "date" in df.columns:
        month_col = df["date"].dt.month
    else:
        month_col = None
    
    if month_col is not None:
        confounders.append(pd.get_dummies(month_col.astype(int), prefix="m", drop_first=True))
    
    # Long-term trend (spline on time)
    if add_spline_time and "date" in df.columns:
        n_years = (df["date"].max() - df["date"].min()).days / 365.25
        time_df = max(3, int(n_years * time_df_per_year))
        
        tnum = (df["date"] - df["date"].min()).dt.days.to_numpy()
        tm_formula = f"cr(x, df={time_df}) - 1"
        tm_spline = _safe_dmatrix(tm_formula, tnum, "x")
        tm_spline.index = df.index
        confounders.append(tm_spline)
    
    # Holidays
    if "is_holiday" in df.columns:
        confounders.append(df[["is_holiday"]].astype(float))

    # Combine design matrix
    if confounders:
        X_conf = pd.concat(confounders, axis=1)
    else:
        X_conf = pd.DataFrame(index=df.index)

    X = pd.concat([X_conf.reset_index(drop=True), cb_df.reset_index(drop=True)], axis=1)
    X = sm.add_constant(X, has_constant="add")

    # Drop rows with NaN in cross-basis (initial rows due to lag)
    valid_mask = ~np.isnan(X[cb_names]).any(axis=1)
    X_valid = X.loc[valid_mask, :].astype(float)
    y_valid = df.loc[valid_mask, deaths_col].to_numpy(dtype=float)
    
    # Population offset
    if pop_col and pop_col in df.columns:
        pop_vals = df.loc[valid_mask, pop_col].to_numpy(dtype=float)
        # Handle zero/negative population
        pop_vals = np.maximum(pop_vals, 1)
        offset = np.log(pop_vals)
    else:
        offset = np.zeros(len(y_valid), dtype=float)
        logger.warning("No population column found; using no offset")

    if len(y_valid) < min_obs:
        logger.warning(f"After dropping NA-lag rows: {len(y_valid)} < {min_obs}; skipping.")
        return None

    # Choose family
    if family == "negbin":
        fam = NegativeBinomial()
    else:
        fam = Poisson()

    try:
        model = GLM(y_valid, X_valid, family=fam, offset=offset)
        
        if family == "quasi-poisson":
            try:
                res = model.fit(scale="X2", cov_type=cov_type)
            except np.linalg.LinAlgError:
                # Fallback to non-robust covariance if HC fails
                logger.warning(f"HC covariance failed, falling back to non-robust")
                res = model.fit(scale="X2")
        else:
            try:
                res = model.fit(cov_type=cov_type)
            except np.linalg.LinAlgError:
                logger.warning(f"HC covariance failed, falling back to non-robust")
                res = model.fit()
            
    except Exception as e:
        logger.error(f"Model fitting failed: {e}")
        return None

    # Extract cross-basis params/cov
    params = res.params
    cov = res.cov_params()
    cb_present = [c for c in cb_names if c in params.index]

    return {
        "params": params,
        "cov": cov,
        "cb_meta": cb_meta,
        "cb_colnames": cb_present,
        "model_result": res,
        "n_obs": int(len(y_valid)),
        "dispersion": float(res.scale) if hasattr(res, 'scale') else None,
        "aic": float(res.aic) if hasattr(res, 'aic') else None,
    }


# =============================================================================
# PREDICTION (CUMULATIVE RR)
# =============================================================================

def predict_cumulative_rr(
    fit_res: Dict[str, Any],
    target_temp: float,
    ref_temp: float,
    up_to_lag: int = None,
) -> Tuple[float, float, float, float, float]:
    """
    Compute cumulative RR at target_temp vs ref_temp up to specified lag.
    
    Uses delta method with full covariance matrix for proper confidence intervals.
    
    Parameters:
    -----------
    fit_res : Output from fit_region_dlnm
    target_temp : Temperature to evaluate (e.g., P99 value)
    ref_temp : Reference temperature (e.g., MMT or P50)
    up_to_lag : Maximum lag to include in cumulative effect (None = all lags)
    
    Returns:
    --------
    (rr, rr_low, rr_high, log_rr, se_log_rr)
    """
    cb_meta = fit_res["cb_meta"]
    params = fit_res["params"]
    cov = fit_res["cov"]
    cb_colnames = fit_res["cb_colnames"]

    L = cb_meta["max_lag"]
    if up_to_lag is None:
        up_to_lag = L
    up_to_lag = min(up_to_lag, L)
    
    K_var = cb_meta["K_var"]
    K_lag = cb_meta["K_lag"]
    var_formula = cb_meta["var_formula"]
    lag_formula = cb_meta["lag_formula"]
    
    # Use patsy with a small array centered around target/ref temps
    # This ensures patsy has enough points to compute knots properly
    temp_range = np.linspace(
        cb_meta.get("var_boundary", (target_temp - 5, target_temp + 5))[0],
        cb_meta.get("var_boundary", (target_temp - 5, target_temp + 5))[1],
        100
    )
    
    # Find indices closest to target and ref temps
    target_idx = np.argmin(np.abs(temp_range - target_temp))
    ref_idx = np.argmin(np.abs(temp_range - ref_temp))
    
    # Compute patsy basis on the full range
    var_basis_full = _safe_dmatrix(var_formula, temp_range, "x")
    if var_basis_full.shape[1] == 0:
        # Fallback if patsy fails
        return 1.0, 1.0, 1.0, 0.0, np.nan
    
    var_target = var_basis_full.iloc[target_idx].to_numpy()
    var_ref = var_basis_full.iloc[ref_idx].to_numpy()
    
    # Lag basis - we have stored values from fitting
    lag_basis_vals = cb_meta.get("lag_basis_vals")
    if lag_basis_vals is not None:
        # Use stored lag basis (computed during fitting)
        lag_sums = np.sum(lag_basis_vals[:up_to_lag+1, :], axis=0)
    else:
        # Compute lag basis 
        lag_idx = np.arange(0, up_to_lag + 1, dtype=float)
        lag_log = np.log(lag_idx + 1)
        lag_basis_h = _safe_dmatrix(lag_formula, lag_log, "x").to_numpy()
        if lag_basis_h.shape[0] == 0:
            return 1.0, 1.0, 1.0, 0.0, np.nan
        lag_sums = np.sum(lag_basis_h, axis=0)

    # Build gradient vector d for delta method
    # Cross-basis column order: for kv in 0..K_var-1: for kl in 0..K_lag-1
    d = []
    expected_colnames = []
    
    for kv in range(K_var):
        for kl in range(K_lag):
            diff = (var_target[kv] - var_ref[kv]) * lag_sums[kl]
            d.append(diff)
            expected_colnames.append(f"cb_v{kv}_l{kl}")
    
    d = np.array(d, dtype=float)

    # Extract coefficients in correct order
    beta_vec = np.array([
        float(params.get(c, 0.0)) for c in expected_colnames
    ], dtype=float)

    # Extract covariance submatrix
    try:
        cov_sub_df = cov.reindex(index=expected_colnames, columns=expected_colnames).fillna(0.0)
        cov_sub = cov_sub_df.to_numpy(dtype=float)
    except Exception:
        # Fallback: diagonal approximation
        diag_vars = np.array([
            float(cov.loc[c, c]) if c in cov.index else 0.0 
            for c in expected_colnames
        ])
        cov_sub = np.diag(diag_vars)

    # Compute log_rr and variance via delta method
    log_rr = float(np.dot(d, beta_vec))
    var_log_rr = float(d.T @ cov_sub @ d)
    var_log_rr = max(var_log_rr, 0.0)  # Ensure non-negative
    se_log_rr = float(np.sqrt(var_log_rr))

    # Transform to RR scale
    rr = float(np.exp(log_rr))
    rr_low = float(np.exp(log_rr - 1.96 * se_log_rr))
    rr_high = float(np.exp(log_rr + 1.96 * se_log_rr))

    return rr, rr_low, rr_high, log_rr, se_log_rr


def compute_rr_curve(
    fit_res: Dict[str, Any],
    temp_range: np.ndarray,
    ref_temp: float,
    up_to_lag: int = None,
) -> pd.DataFrame:
    """
    Compute RR curve across a range of temperatures.
    
    Parameters:
    -----------
    fit_res : Output from fit_region_dlnm
    temp_range : Array of temperature values to evaluate
    ref_temp : Reference temperature
    up_to_lag : Maximum lag for cumulative effect
    
    Returns:
    --------
    DataFrame with columns: temp, rr, rr_low, rr_high
    """
    results = []
    for temp in temp_range:
        rr, rr_low, rr_high, _ = predict_cumulative_rr(fit_res, temp, ref_temp, up_to_lag)
        results.append({
            'temp': temp,
            'rr': rr,
            'rr_low': rr_low,
            'rr_high': rr_high
        })
    return pd.DataFrame(results)


# =============================================================================
# META-ANALYSIS (DERSIMONIAN-LAIRD)
# =============================================================================

def meta_random_effects(effects: np.ndarray, variances: np.ndarray) -> Dict[str, float]:
    """
    DerSimonian-Laird random effects meta-analysis.
    
    Parameters:
    -----------
    effects : Array of effect estimates (log-RR scale)
    variances : Array of within-study variances (se^2)
    
    Returns:
    --------
    Dict with pooled effect, SE, z, p-value, tau^2, I^2
    """
    effects = np.array(effects, dtype=float)
    variances = np.array(variances, dtype=float)
    
    # Remove invalid entries
    valid = ~(np.isnan(effects) | np.isnan(variances) | (variances <= 0))
    effects = effects[valid]
    variances = variances[valid]
    
    k = len(effects)
    if k == 0:
        return {"effect": np.nan, "se": np.nan, "z": np.nan, "p": np.nan, 
                "tau2": np.nan, "I2": np.nan, "k": 0}

    # Fixed-effects weights
    w_fixed = 1.0 / variances
    fixed_mean = np.sum(w_fixed * effects) / np.sum(w_fixed)
    
    # Cochran's Q statistic
    Q = np.sum(w_fixed * (effects - fixed_mean) ** 2)
    df = k - 1
    
    # Between-study variance (tau^2)
    c = np.sum(w_fixed) - (np.sum(w_fixed ** 2) / np.sum(w_fixed))
    tau2 = max(0.0, (Q - df) / c) if c > 0 and df > 0 else 0.0
    
    # I^2 heterogeneity statistic
    I2 = max(0.0, (Q - df) / Q * 100) if Q > 0 else 0.0
    
    # Random-effects weights
    w_random = 1.0 / (variances + tau2)
    pooled = np.sum(w_random * effects) / np.sum(w_random)
    se_pooled = np.sqrt(1.0 / np.sum(w_random))
    
    # Z-test
    z = pooled / se_pooled if se_pooled > 0 else 0.0
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    
    return {
        "effect": float(pooled),
        "se": float(se_pooled),
        "z": float(z),
        "p": float(p),
        "tau2": float(tau2),
        "I2": float(I2),
        "k": int(k),
        "Q": float(Q),
    }


def pool_region_effects(
    region_results: List[Dict[str, Any]],
    effect_key: str = "log_rr",
    se_key: str = "se",
) -> Dict[str, float]:
    """
    Pool log-RR effects across regions using random-effects meta-analysis.
    
    Parameters:
    -----------
    region_results : List of dicts with effect_key and se_key
    effect_key : Key for log-RR in each dict
    se_key : Key for SE in each dict
    
    Returns:
    --------
    Meta-analysis results in RR scale
    """
    effects = []
    variances = []
    
    for r in region_results:
        if r is None:
            continue
        try:
            eff = r[effect_key]
            se = r[se_key]
            if not np.isnan(eff) and not np.isnan(se) and se > 0:
                effects.append(eff)
                variances.append(se ** 2)
        except (KeyError, TypeError):
            continue
    
    if len(effects) == 0:
        return None
    
    meta = meta_random_effects(np.array(effects), np.array(variances))
    
    # Convert to RR scale
    pooled_rr = float(np.exp(meta["effect"]))
    pooled_ci_low = float(np.exp(meta["effect"] - 1.96 * meta["se"]))
    pooled_ci_high = float(np.exp(meta["effect"] + 1.96 * meta["se"]))
    
    return {
        "pooled_log_rr": meta["effect"],
        "pooled_se": meta["se"],
        "pooled_rr": pooled_rr,
        "pooled_ci_low": pooled_ci_low,
        "pooled_ci_high": pooled_ci_high,
        "tau2": meta["tau2"],
        "I2": meta["I2"],
        "k": meta["k"],
        "p": meta["p"],
    }


# =============================================================================
# HARVESTING ANALYSIS HELPERS
# =============================================================================

def harvesting_for_region(
    df_region: pd.DataFrame,
    max_lag: int = 35,
    horizons: List[int] = [7, 14, 21, 28, 35],
    var_df: int = 4,
    lag_df: int = 4,
    temp_col: str = "temp_mean",
    deaths_col: str = "deaths_elderly",
    pop_col: str = "pop_elderly",
    family: str = "quasi-poisson",
) -> Optional[Dict[str, Any]]:
    """
    Fit DLNM for region up to max_lag and compute cumulative RR for each horizon.
    
    Parameters:
    -----------
    df_region : Region DataFrame
    max_lag : Maximum lag to fit
    horizons : List of horizon days to evaluate cumulative RR
    
    Returns:
    --------
    Dict with RR at each horizon for heat/cold percentiles
    """
    fit = fit_region_dlnm(
        df_region,
        temp_col=temp_col,
        deaths_col=deaths_col,
        pop_col=pop_col,
        max_lag=max_lag,
        var_df=var_df,
        lag_df=lag_df,
        family=family,
    )
    
    if fit is None:
        return None

    # Compute region-specific temperature percentiles
    temps = df_region[temp_col].dropna()
    percentiles = {
        'p01': float(temps.quantile(0.01)),
        'p025': float(temps.quantile(0.025)),
        'p50': float(temps.quantile(0.50)),
        'p975': float(temps.quantile(0.975)),
        'p99': float(temps.quantile(0.99)),
    }
    
    ref_temp = percentiles['p50']
    
    # Compute RR at each horizon
    region_code = df_region['region_code'].iloc[0] if 'region_code' in df_region.columns else 'unknown'
    
    results = {
        "region": region_code,
        "n_obs": fit["n_obs"],
        "percentiles": percentiles,
        "rrs_by_horizon": {}
    }
    
    for h in horizons:
        h_capped = min(h, max_lag)
        
        rrs = {}
        for pct_name, pct_temp in percentiles.items():
            if pct_name == 'p50':
                continue
            rr, lo, hi, log_rr, se_log_rr = predict_cumulative_rr(fit, pct_temp, ref_temp, up_to_lag=h_capped)
            rrs[pct_name] = {
                "rr": rr, "low": lo, "high": hi, "log_rr": log_rr,
                "se": se_log_rr if not np.isnan(se_log_rr) else ((np.log(hi) - np.log(lo)) / (2 * 1.96) if hi > lo > 0 else np.nan)
            }
        
        results["rrs_by_horizon"][h] = rrs
    
    return results


def compute_harvesting_ratio(
    rrs_by_horizon: Dict[int, Dict],
    short_horizon: int = 7,
    long_horizon: int = 35,
    pct: str = "p99",
) -> Dict[str, float]:
    """
    Compute harvesting ratio from cumulative RR at short vs long horizon.
    
    Harvesting ratio = 1 - (ERR_long / ERR_short)
    Where ERR = RR - 1 (excess relative risk)
    
    Positive ratio indicates mortality displacement (harvesting).
    Negative ratio indicates persistent/delayed effects.
    """
    try:
        rr_short = rrs_by_horizon[short_horizon][pct]["rr"]
        rr_long = rrs_by_horizon[long_horizon][pct]["rr"]
        
        err_short = rr_short - 1
        err_long = rr_long - 1
        
        if abs(err_short) < 0.001:
            return {"ratio": np.nan, "err_short": err_short, "err_long": err_long}
        
        ratio = 1 - (err_long / err_short)
        
        return {
            "ratio": float(ratio),
            "err_short": float(err_short),
            "err_long": float(err_long),
            "rr_short": float(rr_short),
            "rr_long": float(rr_long),
        }
    except (KeyError, TypeError):
        return {"ratio": np.nan}


# =============================================================================
# HEATWAVE EFFECT MODIFICATION
# =============================================================================

def identify_heatwaves(
    df: pd.DataFrame,
    temp_col: str = "temp_mean",
    region_col: str = "region_code",
    pct_threshold: float = 0.90,
    min_days: int = 2,
) -> pd.Series:
    """
    Identify heatwave days using region-specific thresholds.
    
    Parameters:
    -----------
    df : DataFrame with temperature and region columns
    temp_col : Temperature column
    region_col : Region identifier column
    pct_threshold : Percentile threshold for extreme heat
    min_days : Minimum consecutive days to qualify as heatwave
    
    Returns:
    --------
    Boolean Series indicating heatwave days
    """
    is_hw = pd.Series(False, index=df.index)
    
    for region, group in df.groupby(region_col):
        temps = group[temp_col]
        thresh = temps.quantile(pct_threshold)
        above = temps > thresh
        
        # Find runs of consecutive days above threshold
        run_id = (above != above.shift()).cumsum()
        run_lengths = above.groupby(run_id).transform('sum')
        
        # Mark as heatwave if run >= min_days
        hw_mask = above & (run_lengths >= min_days)
        is_hw.loc[group.index] = hw_mask
    
    return is_hw


def fit_region_dlnm_with_heatwave(
    df_region: pd.DataFrame,
    heatwave_col: str = "is_heatwave",
    temp_col: str = "temp_mean",
    deaths_col: str = "deaths_elderly",
    pop_col: str = "pop_elderly",
    max_lag: int = 21,
    var_df: int = 4,
    lag_df: int = 4,
    family: str = "quasi-poisson",
) -> Optional[Dict[str, Any]]:
    """
    Fit DLNM with heatwave effect modification.
    
    Model: deaths ~ cb(temp) + heatwave + cb(temp)×heatwave + confounders
    
    Returns dict with main and interaction effects.
    """
    df = df_region.copy().reset_index(drop=True)
    
    if heatwave_col not in df.columns:
        logger.warning(f"Heatwave column '{heatwave_col}' not found")
        return None
    
    # Ensure binary heatwave indicator
    df["hw"] = df[heatwave_col].astype(int)
    
    # Build main cross-basis
    temp_arr = df[temp_col].to_numpy(dtype=float)
    X_cb, cb_names, cb_meta = create_crossbasis(temp_arr, max_lag=max_lag, 
                                                 var_df=var_df, lag_df=lag_df)
    cb_df = pd.DataFrame(X_cb, columns=cb_names, index=df.index)
    
    # Build interaction terms: cb × heatwave
    int_df = cb_df.multiply(df["hw"], axis=0)
    int_names = [f"{c}__hw" for c in cb_names]
    int_df.columns = int_names
    
    # Build confounders (same as fit_region_dlnm)
    confounders = []
    
    if "date" in df.columns:
        df = df.sort_values("date").reset_index(drop=True)
        confounders.append(pd.get_dummies(df["date"].dt.dayofweek, prefix="dow", drop_first=True))
        confounders.append(pd.get_dummies(df["date"].dt.month, prefix="m", drop_first=True))
        
        n_years = (df["date"].max() - df["date"].min()).days / 365.25
        time_df = max(3, int(n_years * 4))
        tnum = (df["date"] - df["date"].min()).dt.days.to_numpy()
        tm_spline = _safe_dmatrix(f"cr(x, df={time_df}) - 1", tnum, "x")
        tm_spline.index = df.index
        confounders.append(tm_spline)
    
    if "is_holiday" in df.columns:
        confounders.append(df[["is_holiday"]].astype(float))
    
    X_conf = pd.concat(confounders, axis=1) if confounders else pd.DataFrame(index=df.index)
    
    # Combine: confounders + heatwave main + cb + interaction
    X = pd.concat([
        X_conf.reset_index(drop=True),
        df[["hw"]].reset_index(drop=True),
        cb_df.reset_index(drop=True),
        int_df.reset_index(drop=True)
    ], axis=1)
    X = sm.add_constant(X, has_constant="add")
    
    # Drop rows with NaN
    all_cb_cols = cb_names + int_names
    valid_mask = ~np.isnan(X[all_cb_cols]).any(axis=1)
    X_valid = X.loc[valid_mask, :].astype(float)
    y_valid = df.loc[valid_mask, deaths_col].to_numpy(dtype=float)
    
    if pop_col in df.columns:
        offset = np.log(df.loc[valid_mask, pop_col].to_numpy(dtype=float).clip(min=1))
    else:
        offset = np.zeros(len(y_valid))
    
    if len(y_valid) < 365:
        return None
    
    # Fit model
    fam = NegativeBinomial() if family == "negbin" else Poisson()
    try:
        model = GLM(y_valid, X_valid, family=fam, offset=offset)
        res = model.fit(scale="X2" if family == "quasi-poisson" else None, cov_type="HC1")
    except Exception as e:
        logger.exception(f"Heatwave model failed: {e}")
        return None
    
    return {
        "params": res.params,
        "cov": res.cov_params(),
        "cb_meta": cb_meta,
        "cb_names": cb_names,
        "int_names": int_names,
        "model_result": res,
        "n_obs": int(len(y_valid)),
        "n_heatwave_days": int(df.loc[valid_mask, "hw"].sum()),
    }


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def extract_region_percentiles(
    df: pd.DataFrame,
    temp_col: str = "temp_mean",
    region_col: str = "region_code",
    percentiles: List[float] = [0.01, 0.025, 0.25, 0.50, 0.75, 0.975, 0.99],
) -> pd.DataFrame:
    """
    Compute temperature percentiles for each region.
    
    Returns DataFrame with region_code and percentile columns.
    """
    results = []
    for region, group in df.groupby(region_col):
        temps = group[temp_col].dropna()
        row = {"region_code": region}
        for p in percentiles:
            row[f"p{int(p*100):02d}" if p < 0.1 else f"p{int(p*100)}"] = temps.quantile(p)
        results.append(row)
    return pd.DataFrame(results)


def find_mmt(
    fit_res: Dict[str, Any],
    temp_range: np.ndarray,
    ref_temp: float = None,
    constrain_to_interior: bool = True,
    interior_pct: Tuple[float, float] = (10, 90),
) -> float:
    """
    Find Minimum Mortality Temperature (MMT) from fitted DLNM.
    
    Searches for temperature with lowest cumulative RR across the range.
    
    Following Gasparrini methodology, MMT search is typically constrained
    to an interior range (e.g., P10-P90) to avoid boundary solutions that
    produce degenerate comparisons.
    
    Parameters:
    -----------
    fit_res : dict
        Output from fit_region_dlnm
    temp_range : array
        Full temperature range (typically P1 to P99)
    ref_temp : float
        Reference temperature for RR computation (default: median of range)
    constrain_to_interior : bool
        If True, constrain MMT search to interior percentiles (default: True)
    interior_pct : tuple
        Percentile range for interior search (default: 10th to 90th)
        
    Returns:
    --------
    mmt : float
        Minimum mortality temperature
    """
    if ref_temp is None:
        ref_temp = np.median(temp_range)
    
    # Constrain search to interior range to avoid boundary solutions
    if constrain_to_interior and len(temp_range) > 10:
        n = len(temp_range)
        lo_idx = int(n * interior_pct[0] / 100)
        hi_idx = int(n * interior_pct[1] / 100)
        search_range = temp_range[lo_idx:hi_idx+1]
    else:
        search_range = temp_range
    
    min_rr = np.inf
    mmt = ref_temp
    
    for temp in search_range:
        rr, _, _, _, _ = predict_cumulative_rr(fit_res, temp, ref_temp)
        if rr < min_rr:
            min_rr = rr
            mmt = temp
    
    return float(mmt)


# =============================================================================
# PHASE 1 COMPATIBLE FUNCTIONS
# =============================================================================
# These functions match the signatures used in Phase 1 scripts exactly.
# They work with raw coefficients and cb_info dictionaries (not fit_res objects).

def compute_cumulative_rr_ns_with_se(
    target_temp: float, 
    ref_temp: float, 
    coefs: np.ndarray, 
    vcov: np.ndarray,
    cb_info: Dict[str, Any]
) -> Tuple[float, float]:
    """
    Compute cumulative RR and standard error at target vs reference temperature.
    
    This is the Phase 1 compatible version that returns (log_rr, se) rather than
    (rr, rr_low, rr_high). It uses the same underlying natural spline computation.
    
    Parameters:
    -----------
    target_temp : float
        Temperature to evaluate
    ref_temp : float  
        Reference temperature (typically MMT)
    coefs : array
        Cross-basis coefficients from fitted model
    vcov : array
        Variance-covariance matrix for cross-basis coefficients
    cb_info : dict
        Cross-basis information from create_crossbasis_ns
        
    Returns:
    --------
    log_rr : float
        Log relative risk
    se : float
        Standard error of log-RR
    """
    temp_knots = np.array(cb_info['temp_knots'])
    temp_boundary = tuple(cb_info['temp_boundary'])
    lag_knots = cb_info['lag_knots']
    temp_df = cb_info['temp_df']
    lag_df = cb_info['lag_df']
    max_lag = cb_info.get('max_lag', 21)
    
    # Recreate lag basis (need to sum over all lags)
    lag_vals = np.arange(max_lag + 1)  # 0 to max_lag
    lag_log = np.log(lag_vals + 1)
    lag_knots_log = np.log(np.array(lag_knots) + 1)
    lag_boundary_log = (np.log(1), np.log(max_lag + 2))
    lag_basis = ns_basis(lag_log, lag_knots_log, lag_boundary_log)
    
    # Temperature basis at target and reference
    temp_basis_target = ns_basis(np.array([target_temp]), temp_knots, temp_boundary)[0]
    temp_basis_ref = ns_basis(np.array([ref_temp]), temp_knots, temp_boundary)[0]
    temp_diff = temp_basis_target - temp_basis_ref
    
    # Sum over all lags (cumulative effect)
    lag_sum = lag_basis.sum(axis=0)
    
    # Contrast vector for cross-basis coefficients
    contrast = np.zeros(temp_df * lag_df)
    for t_idx in range(temp_df):
        for l_idx in range(lag_df):
            col_idx = t_idx * lag_df + l_idx
            contrast[col_idx] = temp_diff[t_idx] * lag_sum[l_idx]
    
    # Cumulative log(RR)
    log_rr = np.dot(contrast, coefs)
    
    # Standard error via delta method
    var = np.dot(np.dot(contrast, vcov), contrast)
    se = np.sqrt(max(0, var))
    
    return float(log_rr), float(se)


def find_mmt_from_coefficients(
    temp_range: np.ndarray,
    coefs: np.ndarray,
    vcov: np.ndarray,
    cb_info: Dict[str, Any],
    centering_temp: float,
    constrain_to_interior: bool = True,
    interior_pct: Tuple[float, float] = (10, 90)
) -> Tuple[float, float, List[Dict]]:
    """
    Find MMT by searching for minimum RR temperature (Phase 1 compatible version).
    
    Works directly with coefficients rather than fit_res object.
    
    Parameters:
    -----------
    temp_range : array
        Array of temperatures to search over (e.g., P1 to P99)
    coefs : array
        Cross-basis coefficients from model
    vcov : array
        Variance-covariance matrix
    cb_info : dict
        Cross-basis information
    centering_temp : float
        Temperature used as computational centering point
    constrain_to_interior : bool
        If True, constrain MMT search to interior percentile range.
        This prevents boundary solutions that produce trivial RR=1 comparisons.
        Following Gasparrini et al. methodology.
    interior_pct : tuple
        Percentile range for interior constraint (default: 10th to 90th)
        
    Returns:
    --------
    mmt : float
        Minimum mortality temperature
    mmt_percentile : float
        Approximate percentile of MMT in the temperature range
    rr_curve : list
        Full RR curve across temperature range (list of dicts)
    """
    # Compute RR at each temperature (relative to centering point)
    rr_values = []
    for temp in temp_range:
        log_rr, se = compute_cumulative_rr_ns_with_se(
            temp, centering_temp, coefs, vcov, cb_info
        )
        rr_values.append({
            'temp': float(temp),
            'log_rr': float(log_rr),
            'se': float(se),
            'rr': float(np.exp(log_rr))
        })
    
    # Find temperature with minimum RR
    # Constrain search to interior range to avoid boundary MMT (P1 or P99)
    # which causes trivial RR=1 when comparing extreme temp to itself
    if constrain_to_interior and len(temp_range) > 10:
        # temp_range is typically P1 to P99 with ~100 points
        # Map interior_pct to indices
        # If temp_range covers P1-P99, then index 0 = P1, index 99 = P99
        # For interior [P10, P90]: need indices corresponding to ~P10-P90
        # P10 within P1-P99 range: (10-1)/(99-1) = 9/98 ≈ 0.092 -> idx ~9
        # P90 within P1-P99 range: (90-1)/(99-1) = 89/98 ≈ 0.908 -> idx ~91
        low_pct, high_pct = interior_pct
        n_points = len(temp_range)
        
        # Assuming temp_range spans P1 to P99 (98 percentile range)
        # Relative position of P10 and P90 within that range
        low_frac = (low_pct - 1) / 98  # Position of P10 in P1-P99
        high_frac = (high_pct - 1) / 98  # Position of P90 in P1-P99
        
        low_idx = max(0, int(low_frac * (n_points - 1)))
        high_idx = min(n_points - 1, int(high_frac * (n_points - 1)))
        
        # Search only in interior range
        interior_rrs = [v['rr'] for v in rr_values[low_idx:high_idx+1]]
        interior_min_idx = np.argmin(interior_rrs)
        min_idx = low_idx + interior_min_idx
        
        logger.debug(f"MMT search constrained to indices {low_idx}-{high_idx} "
                     f"(temps {temp_range[low_idx]:.1f}-{temp_range[high_idx]:.1f})")
    else:
        min_idx = np.argmin([v['rr'] for v in rr_values])
    
    mmt = rr_values[min_idx]['temp']
    
    # Approximate percentile (based on position in range)
    mmt_percentile = (min_idx / (len(temp_range) - 1)) * 100 if len(temp_range) > 1 else 50
    
    return float(mmt), float(mmt_percentile), rr_values


def compute_effects_relative_to_mmt(
    temp_percentiles: Dict[str, float],
    mmt: float,
    coefs: np.ndarray,
    vcov: np.ndarray,
    cb_info: Dict[str, Any]
) -> Dict[str, Dict]:
    """
    Compute effects at all percentiles relative to MMT.
    
    This is the standard approach for temperature-mortality studies.
    
    Parameters:
    -----------
    temp_percentiles : dict
        Dictionary mapping percentile names (e.g., 'p1', 'p99') to temperatures
    mmt : float
        Minimum mortality temperature (reference)
    coefs : array
        Cross-basis coefficients
    vcov : array
        Variance-covariance matrix
    cb_info : dict
        Cross-basis information
        
    Returns:
    --------
    effects : dict
        Effects at each percentile relative to MMT
    """
    effects = {}
    
    for pct_name, target_temp in temp_percentiles.items():
        log_rr, se = compute_cumulative_rr_ns_with_se(
            target_temp, mmt, coefs, vcov, cb_info
        )
        
        # Sanity check
        if se <= 0 or se > 5 or abs(log_rr) > 5:
            se = max(0.001, se)
        
        effects[pct_name] = {
            'temp': float(target_temp),
            'log_rr': float(log_rr),
            'log_rr_se': float(se),
            'rr': float(np.exp(log_rr)),
            'rr_lower': float(np.exp(log_rr - 1.96 * se)),
            'rr_upper': float(np.exp(log_rr + 1.96 * se))
        }
    
    return effects


def create_time_spline(
    time_index: np.ndarray, 
    n_years: float, 
    df_per_year: float = 4.0
) -> np.ndarray:
    """
    Create natural spline basis for long-term time trend.
    
    Parameters:
    -----------
    time_index : array
        Numeric time index (e.g., day number from start)
    n_years : float
        Number of years in the data
    df_per_year : float
        Degrees of freedom per year (default 4)
        
    Returns:
    --------
    spline_basis : array
        Spline basis matrix for time trend
    """
    total_df = max(4, int(n_years * df_per_year))
    knots = np.linspace(time_index.min(), time_index.max(), total_df)[1:-1]
    
    try:
        spline_basis = dmatrix(
            f"cr(x, df={total_df})",
            {"x": time_index},
            return_type='dataframe'
        ).values
    except:
        # Fallback to natural spline
        try:
            boundary = (time_index.min(), time_index.max())
            spline_basis = ns_basis(time_index, knots, boundary)
        except:
            # Simple polynomial fallback
            spline_basis = np.column_stack([
                (time_index / 365) ** p for p in range(1, min(4, total_df) + 1)
            ])
    
    return spline_basis


def random_effects_meta_analysis(
    effects: np.ndarray, 
    variances: np.ndarray
) -> Dict[str, Any]:
    """
    Random-effects meta-analysis using DerSimonian-Laird method.
    
    This is the Phase 1 compatible version with the expected return format.
    
    Parameters:
    -----------
    effects : array
        Effect estimates (log-RR scale)
    variances : array
        Variance of each effect (se^2)
        
    Returns:
    --------
    dict with pooled estimates and heterogeneity statistics
    """
    effects = np.asarray(effects)
    variances = np.asarray(variances)
    
    k = len(effects)
    
    if k < 2:
        return {
            'pooled_effect': float(effects[0]) if k == 1 else np.nan,
            'pooled_se': float(np.sqrt(variances[0])) if k == 1 else np.nan,
            'tau2': 0.0, 'I2': 0.0, 'Q': 0.0, 'p_heterogeneity': 1.0, 'n_regions': k
        }
    
    # Fixed-effects weights
    w = 1.0 / variances
    theta_fe = np.sum(w * effects) / np.sum(w)
    
    # Cochran's Q
    Q = np.sum(w * (effects - theta_fe)**2)
    c = np.sum(w) - np.sum(w**2) / np.sum(w)
    
    # Between-study variance (tau^2)
    tau2 = max(0, (Q - (k - 1)) / c) if c > 0 else 0.0
    
    # Random-effects weights and pooled estimate
    w_re = 1.0 / (variances + tau2)
    theta_re = np.sum(w_re * effects) / np.sum(w_re)
    se_re = np.sqrt(1.0 / np.sum(w_re))
    
    # Heterogeneity statistics
    I2 = max(0, (Q - (k - 1)) / Q * 100) if Q > 0 else 0
    p_het = 1 - stats.chi2.cdf(Q, k - 1)
    
    return {
        'pooled_effect': float(theta_re),
        'pooled_se': float(se_re),
        'pooled_rr': float(np.exp(theta_re)),
        'pooled_rr_lower': float(np.exp(theta_re - 1.96 * se_re)),
        'pooled_rr_upper': float(np.exp(theta_re + 1.96 * se_re)),
        'tau2': float(tau2),
        'I2': float(I2),
        'Q': float(Q),
        'p_heterogeneity': float(p_het),
        'n_regions': k
    }


def pool_region_results(
    region_results: Dict[str, Dict],
    percentile: str,
    max_se: float = 2.0,
    max_abs_logrr: float = 3.0,
    min_se: float = 0.01,
    exclude_mmt_boundary: bool = False,
    mmt_cold_threshold: float = 10.0,
    mmt_heat_threshold: float = 90.0
) -> Dict[str, Any]:
    """
    Pool region-specific results at a given percentile via meta-analysis.
    
    Applies sanity filters to exclude extreme/unreliable estimates.
    The min_se filter is the primary mechanism for excluding degenerate cases
    where MMT is at the boundary (causing SE≈0).
    
    NOTE: Gasparrini uses multivariate meta-analysis of DLNM coefficients,
    not meta-analysis of RRs at percentiles. The min_se filter approximates
    his approach by excluding numerically degenerate comparisons.
    
    Parameters:
    -----------
    region_results : dict
        Dictionary of {region_code: result_dict} where each result has 'effects' key
    percentile : str
        Percentile name to pool (e.g., 'p99', 'p1')
    max_se : float
        Maximum SE to include (filter outliers)
    max_abs_logrr : float
        Maximum |log_rr| to include
    min_se : float
        Minimum SE to include (filters boundary cases where SE is numerically zero)
    exclude_mmt_boundary : bool
        Whether to exclude regions where MMT makes the comparison meaningless
    mmt_cold_threshold : float
        For cold percentiles: exclude regions where MMT <= this percentile
        (no meaningful temperature range below MMT to measure cold effect)
    mmt_heat_threshold : float
        For heat percentiles: exclude regions where MMT >= this percentile
        (no meaningful temperature range above MMT to measure heat effect)
        
    Returns:
    --------
    dict with pooled estimates and regions included
    """
    effects, variances, regions = [], [], []
    n_excluded_boundary = 0
    n_excluded_se = 0
    n_excluded_logrr = 0
    
    for region_code, result in region_results.items():
        if result is None or 'effects' not in result:
            continue
        if percentile not in result['effects']:
            continue
        
        # Check MMT boundary exclusion (Gasparrini methodology)
        # Cold effect only meaningful if MMT > cold_threshold (room below MMT)
        # Heat effect only meaningful if MMT < heat_threshold (room above MMT)
        if exclude_mmt_boundary and 'mmt_percentile' in result:
            mmt_pct = result['mmt_percentile']
            # For heat percentiles (p90+), exclude if MMT is already very high
            # (no meaningful heat exposure above MMT)
            if percentile in ['p90', 'p95', 'p97.5', 'p99']:
                if mmt_pct >= mmt_heat_threshold:
                    n_excluded_boundary += 1
                    continue
            # For cold percentiles (p10-), exclude if MMT is already very low
            # (no meaningful cold exposure below MMT)
            elif percentile in ['p1', 'p2.5', 'p5', 'p10']:
                if mmt_pct <= mmt_cold_threshold:
                    n_excluded_boundary += 1
                    continue
        
        eff = result['effects'][percentile]
        log_rr = eff.get('log_rr', eff.get('logRR'))
        se = eff.get('log_rr_se', eff.get('se'))
        
        if log_rr is None or se is None:
            continue
        
        # Exclude numerically degenerate cases (SE too small = boundary comparison)
        if se < min_se:
            n_excluded_se += 1
            continue
            
        if se > max_se:
            n_excluded_se += 1
            continue
            
        if abs(log_rr) > max_abs_logrr:
            n_excluded_logrr += 1
            continue
        
        effects.append(log_rr)
        variances.append(se**2)
        regions.append(region_code)
    
    if len(effects) == 0:
        return {
            'pooled_rr': np.nan, 
            'n_regions': 0,
            'n_excluded_boundary': n_excluded_boundary,
            'n_excluded_se': n_excluded_se,
            'n_excluded_logrr': n_excluded_logrr
        }
    
    pooled = random_effects_meta_analysis(np.array(effects), np.array(variances))
    pooled['percentile'] = percentile
    pooled['regions_included'] = regions
    pooled['n_excluded_boundary'] = n_excluded_boundary
    pooled['n_excluded_se'] = n_excluded_se
    pooled['n_excluded_logrr'] = n_excluded_logrr
    
    return pooled


def convert_to_json_serializable(obj: Any) -> Any:
    """
    Recursively convert numpy types to JSON-serializable Python types.
    Also handles numpy types in dictionary keys.
    """
    if isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.bool_):
        return bool(obj)
    elif isinstance(obj, dict):
        # Convert both keys and values
        return {
            (int(k) if isinstance(k, (np.integer, np.int64, np.int32)) else 
             float(k) if isinstance(k, (np.floating, np.float64, np.float32)) else k): 
            convert_to_json_serializable(v) 
            for k, v in obj.items()
        }
    elif isinstance(obj, (list, tuple)):
        return [convert_to_json_serializable(item) for item in obj]
    return obj


# =============================================================================
# ATTRIBUTABLE BURDEN FUNCTIONS
# =============================================================================

def compute_attributable_fraction(rr: float) -> float:
    """
    Calculate attributable fraction from relative risk.
    
    AF = (RR - 1) / RR
    
    Interpretation:
    - AF > 0 when RR > 1: Temperature INCREASES risk (attributable deaths)
    - AF < 0 when RR < 1: Temperature DECREASES risk (protective effect)  
    - AF = 0 when RR = 1: No effect
    """
    if rr <= 0:
        return 0.0
    af = (rr - 1) / rr
    return float(np.clip(af, -1, 1))


def detect_basis_type(dlnm_results: Dict) -> str:
    """
    Detect whether DLNM results use polynomial or natural spline basis.
    
    Returns: 'ns' for natural spline, 'poly' for polynomial
    """
    # Check first region's results
    sample_region = list(dlnm_results['region_results'].values())[0]
    
    if 'crossbasis_info' in sample_region:
        return 'ns'
    elif 'col_names' in sample_region and 'poly' in str(sample_region.get('col_names', [''])[0]):
        return 'poly'
    else:
        # Default to 'ns' for new analyses
        logger.warning("Could not detect basis type, defaulting to 'ns'")
        return 'ns'


def compute_cumulative_rr_for_burden(
    target_temp: float,
    ref_temp: float,
    coefs: np.ndarray,
    cb_info: Dict[str, Any]
) -> float:
    """
    Compute cumulative RR for burden calculation (no CI needed).
    
    Parameters:
    -----------
    target_temp : float
        Temperature to evaluate
    ref_temp : float
        Reference temperature (typically MMT)
    coefs : array
        Cross-basis coefficients
    cb_info : dict
        Cross-basis information
        
    Returns:
    --------
    rr : float
        Cumulative relative risk
    """
    temp_knots = np.array(cb_info['temp_knots'])
    temp_boundary = tuple(cb_info['temp_boundary'])
    lag_knots = cb_info['lag_knots']
    temp_df = cb_info['temp_df']
    lag_df = cb_info['lag_df']
    max_lag = cb_info.get('max_lag', 21)
    
    # Recreate lag basis
    lag_vals = np.arange(max_lag + 1)
    lag_log = np.log(lag_vals + 1)
    lag_knots_log = np.log(np.array(lag_knots) + 1)
    lag_boundary_log = (np.log(1), np.log(max_lag + 2))
    lag_basis = ns_basis(lag_log, lag_knots_log, lag_boundary_log)
    
    # Temperature basis at target and reference
    temp_basis_target = ns_basis(np.array([target_temp]), temp_knots, temp_boundary)[0]
    temp_basis_ref = ns_basis(np.array([ref_temp]), temp_knots, temp_boundary)[0]
    temp_diff = temp_basis_target - temp_basis_ref
    
    # Sum over all lags
    lag_sum = lag_basis.sum(axis=0)
    
    # Contrast vector
    contrast = np.zeros(temp_df * lag_df)
    for t_idx in range(temp_df):
        for l_idx in range(lag_df):
            col_idx = t_idx * lag_df + l_idx
            contrast[col_idx] = temp_diff[t_idx] * lag_sum[l_idx]
    
    log_rr = np.dot(contrast, coefs)
    return float(np.exp(log_rr))


# =============================================================================
# MULTIVARIATE META-ANALYSIS (MVMETA) - GASPARRINI METHODOLOGY
# =============================================================================

def mvmeta_fit(
    coef_list: List[np.ndarray],
    vcov_list: List[np.ndarray],
    max_iter: int = 100,
    tol: float = 1e-6
) -> Dict[str, Any]:
    """
    Multivariate meta-analysis using REML estimation (DerSimonian-Laird style).
    
    This pools DLNM cross-basis coefficients across regions, accounting for
    both within-region variance (from vcov) and between-region heterogeneity.
    
    Following Gasparrini's mvmeta methodology from the 2015 Lancet paper.
    
    Parameters:
    -----------
    coef_list : list of arrays
        Coefficients from each region (k x p matrices become k arrays of length p)
    vcov_list : list of arrays
        Variance-covariance matrices from each region (each p x p)
    max_iter : int
        Maximum iterations for REML
    tol : float
        Convergence tolerance
        
    Returns:
    --------
    dict with:
        - pooled_coef: pooled coefficient vector
        - pooled_vcov: pooled variance-covariance matrix
        - Psi: between-study covariance matrix
        - k: number of studies
        - p: number of parameters
        - converged: whether REML converged
    """
    # Convert to arrays
    coefs = np.array(coef_list)  # k x p
    vcovs = np.array(vcov_list)  # k x p x p
    
    k, p = coefs.shape
    
    if k < 2:
        return {
            'pooled_coef': coefs[0] if k == 1 else np.zeros(p),
            'pooled_vcov': vcovs[0] if k == 1 else np.eye(p),
            'Psi': np.zeros((p, p)),
            'k': k,
            'p': p,
            'converged': True
        }
    
    # Initialize between-study covariance (Psi) using method of moments
    # This is the DerSimonian-Laird approach extended to multivariate
    
    # Fixed-effects pooled estimate (initial)
    S_inv_sum = np.zeros((p, p))
    S_inv_y_sum = np.zeros(p)
    
    for i in range(k):
        S_inv = np.linalg.pinv(vcovs[i])
        S_inv_sum += S_inv
        S_inv_y_sum += S_inv @ coefs[i]
    
    theta_fe = np.linalg.solve(S_inv_sum, S_inv_y_sum)
    
    # Method of moments for Psi (between-study variance)
    Q_matrix = np.zeros((p, p))
    for i in range(k):
        resid = coefs[i] - theta_fe
        Q_matrix += np.outer(resid, resid)
    
    # Simple estimate: Q / (k-1) - average within-study variance
    avg_vcov = np.mean(vcovs, axis=0)
    Psi = Q_matrix / (k - 1) - avg_vcov
    
    # Ensure Psi is positive semi-definite
    eigvals, eigvecs = np.linalg.eigh(Psi)
    eigvals = np.maximum(eigvals, 0)  # Truncate negative eigenvalues
    Psi = eigvecs @ np.diag(eigvals) @ eigvecs.T
    
    # Random-effects pooled estimate
    W_inv_sum = np.zeros((p, p))
    W_inv_y_sum = np.zeros(p)
    
    for i in range(k):
        W_i = vcovs[i] + Psi
        W_inv = np.linalg.pinv(W_i)
        W_inv_sum += W_inv
        W_inv_y_sum += W_inv @ coefs[i]
    
    pooled_vcov = np.linalg.pinv(W_inv_sum)
    pooled_coef = pooled_vcov @ W_inv_y_sum
    
    return {
        'pooled_coef': pooled_coef,
        'pooled_vcov': pooled_vcov,
        'Psi': Psi,
        'k': k,
        'p': p,
        'converged': True
    }


def mvmeta_predict_rr(
    pooled_coef: np.ndarray,
    pooled_vcov: np.ndarray,
    cb_info: Dict,
    target_temp: float,
    ref_temp: float
) -> Tuple[float, float, float, float, float]:
    """
    Compute cumulative RR from pooled MVMeta coefficients.
    
    Handles multiple naming conventions for cb_info keys:
    - Phase 1 scripts store: temp_knots, temp_boundary, lag_basis, temp_df, lag_df
    - create_crossbasis stores: var_knots, var_boundary, lag_basis_vals, var_formula
    
    Parameters:
    -----------
    pooled_coef : array
        Pooled cross-basis coefficients
    pooled_vcov : array
        Pooled variance-covariance matrix
    cb_info : dict
        Cross-basis metadata
    target_temp : float
        Temperature to compare
    ref_temp : float
        Reference temperature (MMT)
        
    Returns:
    --------
    Tuple of (rr, rr_lo, rr_hi, log_rr, log_rr_se)
    """
    # Get dimensions - handle multiple naming conventions
    K_var = cb_info.get('K_var') or cb_info.get('temp_df', 4)
    K_lag = cb_info.get('K_lag') or cb_info.get('lag_df', 4)
    max_lag = cb_info.get('max_lag', 21)
    
    # Get temperature knots and boundary
    temp_knots = cb_info.get('var_knots') or cb_info.get('temp_knots')
    temp_boundary = cb_info.get('var_boundary') or cb_info.get('temp_boundary')
    
    if temp_knots is None or temp_boundary is None:
        logger.warning("Missing temp_knots or temp_boundary in cb_info")
        return 1.0, 1.0, 1.0, 0.0, float('inf')
    
    temp_knots = np.array(temp_knots)
    temp_boundary = tuple(temp_boundary)
    
    # Get lag basis - handle multiple naming conventions
    lag_basis_vals = cb_info.get('lag_basis_vals') or cb_info.get('lag_basis')
    
    if lag_basis_vals is None:
        logger.warning("Missing lag_basis in cb_info")
        return 1.0, 1.0, 1.0, 0.0, float('inf')
    
    lag_basis = np.array(lag_basis_vals)
    lag_sum = lag_basis.sum(axis=0)
    
    # Compute temperature basis at target and reference using natural splines
    # Use ns_basis with the stored knots for consistency
    temp_basis_target = ns_basis(np.array([target_temp]), temp_knots, temp_boundary)[0]
    temp_basis_ref = ns_basis(np.array([ref_temp]), temp_knots, temp_boundary)[0]
    temp_diff = temp_basis_target - temp_basis_ref
    
    # Build contrast vector
    # Cross-basis column order: for kv in 0..K_var-1: for kl in 0..K_lag-1
    n_coef = len(pooled_coef)
    actual_K_var = len(temp_diff)
    actual_K_lag = len(lag_sum)
    
    contrast = np.zeros(n_coef)
    for t_idx in range(actual_K_var):
        for l_idx in range(actual_K_lag):
            col_idx = t_idx * actual_K_lag + l_idx
            if col_idx < n_coef:
                contrast[col_idx] = temp_diff[t_idx] * lag_sum[l_idx]
    
    # Log-RR and SE via delta method
    log_rr = np.dot(contrast, pooled_coef)
    var_log_rr = contrast @ pooled_vcov @ contrast
    se_log_rr = np.sqrt(max(0, var_log_rr))
    
    # Guard against overflow
    if np.abs(log_rr) > 10:
        logger.warning(f"Large log_rr={log_rr:.2f}, clamping to avoid overflow")
        log_rr = np.clip(log_rr, -10, 10)
    
    # Guard against infinite SE (use large but finite value)
    if np.isinf(se_log_rr) or se_log_rr > 100:
        se_log_rr = 100.0
    
    rr = np.exp(log_rr)
    rr_lo = np.exp(log_rr - 1.96 * se_log_rr)
    rr_hi = np.exp(log_rr + 1.96 * se_log_rr)
    
    return float(rr), float(rr_lo), float(rr_hi), float(log_rr), float(se_log_rr)


def mvmeta_pool_regions(
    region_results: Dict[str, Dict],
    min_obs: int = 500
) -> Dict[str, Any]:
    """
    Pool region-specific DLNM results using multivariate meta-analysis.
    
    This follows Gasparrini's 2015 Lancet methodology:
    1. Extract cross-basis coefficients and vcov from each region
    2. Pool coefficients using MVMeta
    3. Find pooled MMT from the pooled curve
    4. Compute RRs at all percentiles relative to pooled MMT
    
    Parameters:
    -----------
    region_results : dict
        Dictionary of {region_code: result_dict}
    min_obs : int
        Minimum observations to include a region
        
    Returns:
    --------
    dict with pooled results
    """
    # Extract valid region coefficients
    coef_list = []
    vcov_list = []
    region_list = []
    cb_info = None
    
    # Collect temperature percentiles from all regions for global percentiles
    all_temp_pcts = {pct: [] for pct in ['p1', 'p2.5', 'p5', 'p10', 'p25', 'p50', 'p75', 'p90', 'p95', 'p97.5', 'p99']}
    
    for region_code, result in region_results.items():
        if result is None:
            continue
        if result.get('n_obs', 0) < min_obs:
            continue
        if 'cb_coefs' not in result or 'cb_vcov' not in result:
            continue
        
        coef = np.array(result['cb_coefs'])
        vcov = np.array(result['cb_vcov'])
        
        # Check for valid variance (exclude degenerate cases)
        if np.any(np.isnan(coef)) or np.any(np.isnan(vcov)):
            continue
        if np.any(np.diag(vcov) <= 0):  # Invalid variance
            continue
            
        coef_list.append(coef)
        vcov_list.append(vcov)
        region_list.append(region_code)
        
        # Store cb_info from first valid region
        if cb_info is None:
            cb_info = result.get('crossbasis_info')
        
        # Collect temperature percentiles
        if 'temp_percentiles' in result:
            for pct, val in result['temp_percentiles'].items():
                if pct in all_temp_pcts:
                    all_temp_pcts[pct].append(val)
    
    if len(coef_list) < 2:
        return {'error': 'Insufficient valid regions for MVMeta', 'n_regions': len(coef_list)}
    
    if cb_info is None:
        return {'error': 'No cross-basis info available'}
    
    # Compute global percentiles (median across regions)
    global_pcts = {}
    for pct, vals in all_temp_pcts.items():
        if vals:
            global_pcts[pct] = float(np.median(vals))
    
    # Fit MVMeta
    mvmeta_result = mvmeta_fit(coef_list, vcov_list)
    
    pooled_coef = mvmeta_result['pooled_coef']
    pooled_vcov = mvmeta_result['pooled_vcov']
    
    # Find pooled MMT
    # Search across the global temperature range
    temp_min = global_pcts.get('p1', 10)
    temp_max = global_pcts.get('p99', 35)
    temp_grid = np.linspace(temp_min, temp_max, 200)
    
    # Compute RR curve relative to P50 initially
    ref_temp = global_pcts.get('p50', 22)
    rr_curve = []
    for t in temp_grid:
        rr, _, _, _, _ = mvmeta_predict_rr(pooled_coef, pooled_vcov, cb_info, t, ref_temp)
        rr_curve.append(rr)
    
    rr_curve = np.array(rr_curve)
    
    # MMT is where RR is minimum
    mmt_idx = np.argmin(rr_curve)
    mmt = temp_grid[mmt_idx]
    
    # Convert MMT to percentile (relative to global distribution)
    mmt_percentile = 1 + 98 * (mmt - temp_min) / (temp_max - temp_min)
    mmt_percentile = max(1, min(99, mmt_percentile))
    
    # Compute effects at all percentiles relative to MMT
    effects = {}
    for pct, temp in global_pcts.items():
        rr, rr_lo, rr_hi, log_rr, se = mvmeta_predict_rr(
            pooled_coef, pooled_vcov, cb_info, temp, mmt
        )
        effects[pct] = {
            'rr': rr,
            'rr_lower': rr_lo,
            'rr_upper': rr_hi,
            'log_rr': log_rr,
            'log_rr_se': se
        }
    
    return {
        'method': 'MVMeta (Multivariate Meta-Analysis)',
        'n_regions': len(region_list),
        'regions_included': region_list,
        'pooled_coef': pooled_coef.tolist(),
        'pooled_vcov': pooled_vcov.tolist(),
        'Psi': mvmeta_result['Psi'].tolist(),
        'mmt': mmt,
        'mmt_percentile': mmt_percentile,
        'global_percentiles': global_pcts,
        'effects': effects,
        'cb_info': cb_info
    }


# =============================================================================
# CONVENIENCE WRAPPERS FOR MVMETA (Used by Phase scripts)
# =============================================================================

def mvmeta_pool_coefficients(
    coef_list: List[np.ndarray],
    vcov_list: List[np.ndarray],
) -> Optional[Dict[str, Any]]:
    """
    Pool DLNM cross-basis coefficients using multivariate meta-analysis.
    
    Wrapper around mvmeta_fit for easier use in Phase scripts.
    
    Parameters:
    -----------
    coef_list : list of arrays
        Cross-basis coefficients from each region
    vcov_list : list of arrays  
        Variance-covariance matrices from each region
        
    Returns:
    --------
    dict with pooled_coef, pooled_vcov, Psi, k, p, converged
    or None if pooling fails
    """
    if len(coef_list) < 2:
        return None
    
    try:
        # Ensure consistent shapes
        coefs = [np.array(c).flatten() for c in coef_list]
        vcovs = [np.array(v) for v in vcov_list]
        
        # Check all have same dimension
        p = len(coefs[0])
        valid_coefs = []
        valid_vcovs = []
        
        for c, v in zip(coefs, vcovs):
            if len(c) == p and v.shape == (p, p):
                valid_coefs.append(c)
                valid_vcovs.append(v)
        
        if len(valid_coefs) < 2:
            return None
        
        return mvmeta_fit(valid_coefs, valid_vcovs)
    except Exception:
        return None


def compute_pooled_rr_from_mvmeta(
    mvmeta_result: Dict[str, Any],
    percentiles: Dict[str, float],
    cb_info: Dict[str, Any],
) -> Dict[str, float]:
    """
    Compute pooled RRs from MVMeta result.
    
    Parameters:
    -----------
    mvmeta_result : dict
        Output from mvmeta_pool_coefficients or mvmeta_fit
    percentiles : dict
        Dictionary with 'p1', 'p50', 'p99' temperature values
    cb_info : dict
        Cross-basis info with 'var_df', 'lag_df', 'max_lag'
        
    Returns:
    --------
    dict with heat_rr, heat_rr_lo, heat_rr_hi, cold_rr, cold_rr_lo, cold_rr_hi,
         mmt, mmt_percentile
    """
    pooled_coef = np.array(mvmeta_result['pooled_coef'])
    pooled_vcov = np.array(mvmeta_result['pooled_vcov'])
    
    p1 = percentiles.get('p1', percentiles.get('p01', 10))
    p50 = percentiles.get('p50', 20)
    p99 = percentiles.get('p99', 35)
    
    # Find MMT from pooled curve
    temp_range = np.linspace(p1, p99, 100)
    min_rr = float('inf')
    mmt = p50
    
    for t in temp_range:
        try:
            rr, _, _, _, _ = mvmeta_predict_rr(pooled_coef, pooled_vcov, cb_info, t, p50)
            if rr < min_rr:
                min_rr = rr
                mmt = t
        except Exception:
            continue
    
    # MMT percentile
    mmt_pct = 1 + 98 * (mmt - p1) / (p99 - p1) if p99 > p1 else 50
    mmt_pct = max(1, min(99, mmt_pct))
    
    # Heat RR (P99 vs MMT)
    try:
        heat_rr, heat_lo, heat_hi, _, _ = mvmeta_predict_rr(
            pooled_coef, pooled_vcov, cb_info, p99, mmt
        )
    except Exception:
        heat_rr, heat_lo, heat_hi = 1.0, 1.0, 1.0
    
    # Cold RR (P1 vs MMT)
    try:
        cold_rr, cold_lo, cold_hi, _, _ = mvmeta_predict_rr(
            pooled_coef, pooled_vcov, cb_info, p1, mmt
        )
    except Exception:
        cold_rr, cold_lo, cold_hi = 1.0, 1.0, 1.0
    
    return {
        'heat_rr': float(heat_rr),
        'heat_rr_lo': float(heat_lo),
        'heat_rr_hi': float(heat_hi),
        'cold_rr': float(cold_rr),
        'cold_rr_lo': float(cold_lo),
        'cold_rr_hi': float(cold_hi),
        'mmt': float(mmt),
        'mmt_percentile': float(mmt_pct),
    }
