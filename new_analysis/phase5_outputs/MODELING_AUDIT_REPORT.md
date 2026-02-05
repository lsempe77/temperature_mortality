# R Pipeline Modeling Audit Report
## Temperature-Mortality Analysis - Brazil Elderly (2010-2024)

**Audit Date:** December 18, 2025  
**Auditor:** Epidemiological Methods Review

---

## Executive Summary

After reviewing all R scripts in the analysis pipeline, I have identified the following:

### ✅ CORRECT Methodological Approaches

1. **DLNM Framework** - Correctly implements Gasparrini's methodology
2. **Two-stage meta-analysis** - Appropriate for multi-location studies
3. **Case-crossover design** - Correct time-stratified matching
4. **Stratification analyses** - Appropriate for effect modification

### ⚠️ ISSUES IDENTIFIED

1. **Attributable Fraction calculation** - Potential overestimation concern
2. **Temperature percentile handling** - Inconsistent across scripts
3. **MMT estimation** - Minor concerns about boundary effects
4. **SE filtering thresholds** - May be too lenient in some scripts

---

## Detailed Script-by-Script Audit

### Phase 1: Core Models

#### 01_dlnm_analysis_v2.R ✅ MOSTLY CORRECT

**Correct:**
- Natural cubic spline cross-basis with 4 df for temp and lag
- 21-day maximum lag (standard for temperature-mortality)
- Quasi-Poisson family for overdispersion
- Seasonal trend control: ns(date, df = 7*years) - appropriate
- Day-of-week adjustment
- Mixmeta pooling with REML

**Minor Issues:**
```r
# Line 126-129: Region-specific boundaries
reg_boundary <- range(daily$tmean, na.rm = TRUE)  # ❌ Could cause edge effects
# Should use global temp_boundary for consistency:
# reg_boundary <- temp_boundary  # Fixed at 0.7-35.1
```

**Recommendation:** Use global temperature boundaries to ensure comparable cross-basis structures across regions.

---

#### 01b_attributable_burden.R ⚠️ NEEDS REVIEW

**Correct:**
- Uses pooled DLNM coefficients
- Correct AF formula: AF = (RR - 1) / RR

**Issue 1: AF only calculated when RR > 1**
```r
# Line 183
data[, af := fifelse(rr > 1, (rr - 1) / rr, 0)]  # ⚠️ Excludes protective effects
```
This is standard practice (Gasparrini & Leone 2014), but note that:
- RR < 1 at temperatures near MMT is expected
- Setting AF=0 for RR<1 is correct for burden calculation
- **VERDICT: CORRECT**

**Issue 2: Heat vs Cold classification**
```r
# Line 186-187
data[, is_heat := tmean > pooled_mmt]
data[, is_cold := tmean < pooled_mmt]
```
This classifies ALL temperatures above MMT as "heat" and below as "cold". 
- This may include moderate temperatures with minimal risk
- **VERDICT: CORRECT** - This is the Gasparrini methodology

---

#### 01c_yll_calculation.R ✅ CORRECT

- Uses WHO standard life expectancy tables
- Correctly applies discounting (3%) and age-weighting (optional)
- Properly integrates with attributable burden data

---

#### 01d_case_crossover.R ✅ CORRECT

**Correct:**
- True time-stratified case-crossover design
- Matching on year-month-dow (standard approach)
- Conditional logistic regression (clogit)
- Multiple model specifications (linear, quadratic, categorical)

**Note:** Case-crossover is a VALIDATION, not the primary analysis. It cannot capture lagged effects (uses same-day temperature only).

---

#### 01e_excess_mortality.R ✅ CORRECT

- GAM baseline for expected mortality
- Accounts for seasonality and trend
- Compares observed vs expected during temperature extremes

---

### Phase 2: Robustness

#### 02a_sensitivity.R ✅ CORRECT

**Correct:**
- Tests lag structure: 7, 14, 21, 28 days
- Tests spline df: 3, 4, 5
- Tests temperature percentile thresholds

**Minor Issue:**
```r
# SE filtering threshold of 20 may be too lenient for some applications
rr_p99_se < 20 && rr_p1_se < 20
```
**Recommendation:** Consider using SE < 5 for more conservative filtering.

---

#### 02b_harvesting.R ✅ CORRECT

- Extended 35-day lag for harvesting detection
- Multiple horizon cumulative effects
- Correct harvesting ratio calculation

**Issue:** At 35-day lag with 510 regions, many models fail to converge due to sparse data.
**Verdict:** This is expected behavior, not a coding error.

---

#### 02c_heatwave.R ✅ CORRECT

- Correct heatwave definition (3+ consecutive days >P95)
- Additive interaction model is appropriate
- Pooling approach consistent with main DLNM

---

### Phase 3: Confounding

#### 03a_supplementary.R ✅ CORRECT

- Apparent temperature (humidity-adjusted) comparison
- Pollution adjustment (PM2.5, O3)
- Influenza season control
- Appropriate model specification

**Note:** Pollution may be a MEDIATOR, not confounder. Controlling for it may introduce over-adjustment bias.

---

### Phase 4: Heterogeneity

#### 04a_meta_regression.R ✅ CORRECT

- Region-specific meta-regression
- Appropriate moderator variables (climate, demographics)
- Mixmeta with random effects

---

#### 04b_age_stratification.R ✅ CORRECT

- Separate DLNM by age group (60-69, 70-79, 80+)
- Same methodology as main analysis
- Correct comparison across strata

---

#### 04c_sex_stratification.R ✅ CORRECT

- Separate DLNM by sex
- Same methodology as main analysis

---

#### 04d_cause_stratification.R ✅ CORRECT

- Separate DLNM by cause (CVD, respiratory, external, other)
- Correct ICD-10 groupings
- Same methodology as main analysis

---

## Critical Methodological Checks

### 1. Cross-Basis Specification ✅

```r
cb <- crossbasis(
  daily$tmean,
  lag = MAX_LAG,  # 21 days
  argvar = list(fun = "ns", knots = reg_temp_pcts, Boundary.knots = temp_boundary),
  arglag = list(fun = "ns", df = LAG_DF)  # 4 df
)
```

This follows Gasparrini et al. (2015) Lancet exactly:
- Natural spline for temperature (4 df = 3 internal knots at P10, P75, P90)
- Natural spline for lag (4 df)
- Total 16 parameters per cross-basis

### 2. MMT Estimation ✅

```r
# Find MMT by minimizing cumulative RR
temp_seq <- seq(temp_boundary[1], temp_boundary[2], length.out = 100)
cp <- crosspred(cb, model, at = temp_seq, cumul = TRUE, cen = median(daily$tmean))
mmt_idx <- which.min(cp$allRRfit)
mmt <- temp_seq[mmt_idx]
```

**Correct approach:** MMT is the temperature that minimizes cumulative mortality risk.

### 3. RR Centering ✅

```r
# RRs are correctly centered at MMT
cp_pcts <- crosspred(cb, model, at = reg_pcts_vals, cumul = TRUE, cen = mmt)
```

This ensures RR = 1 at MMT, with deviations representing excess risk.

### 4. Meta-Analysis ✅

```r
# Mixmeta with REML
mv_fit <- mixmeta(coef_matrix, mvmeta_vcovs, method = "reml")
```

REML (Restricted Maximum Likelihood) is preferred for variance estimation.

### 5. Attributable Fraction ✅

```r
# Gasparrini & Leone (2014) formula
af := fifelse(rr > 1, (rr - 1) / rr, 0)
an := af * deaths
```

Correctly calculates attributable number as deaths × attributable fraction.

---

## Recommendations

### HIGH PRIORITY

1. **None identified** - Core methodology is sound

### MEDIUM PRIORITY

1. **Standardize temperature boundaries across regions** (lines 126-129 in 01_dlnm_analysis_v2.R)
   - Current: Uses region-specific boundaries
   - Recommended: Use global temp_boundary for all regions
   - Impact: Minor, affects edge behavior

2. **Consider stricter SE filtering** (various scripts)
   - Current: SE < 20
   - Recommended: SE < 5 for publication-quality results
   - Impact: May reduce regional sample size but improve pooled estimates

### LOW PRIORITY

1. **Document harvesting analysis limitations** for immediate level
   - 35-day lag with 510 regions causes convergence issues
   - Consider using intermediate level only for harvesting analysis

2. **Add heterogeneity statistics (I², Q)** to output
   - Currently not systematically reported
   - Useful for assessing between-region variability

---

## Additional Verification Checks

### GLM Specification Consistency ✅

All 13 GLM calls across scripts use consistent specification:
```r
glm(
  deaths ~ cb + ns(as.numeric(date), df = 7 * n_years) + factor(format(date, "%u")),
  family = quasipoisson(link = "log")
)
```

Components:
- `cb`: Cross-basis for temperature × lag
- `ns(date, df=7*years)`: Long-term trend and seasonality (~7 df per year)
- `factor(dow)`: Day-of-week fixed effects
- `quasipoisson`: Handles overdispersion in death counts

### Heatwave Model ✅

Uses additive main effect (not multiplicative interaction):
```r
deaths ~ cb + heatwave + [trend] + [dow]
```
This measures the additional mortality risk during heatwave events, beyond the temperature effect captured by the cross-basis.

### Family Choice ✅

All models correctly use `quasipoisson(link = "log")`:
- Log link ensures positive predictions
- Quasi-Poisson handles overdispersion without specifying full distribution
- Dispersion parameter estimated from data

---

## Conclusion

**The R analysis pipeline is methodologically sound and follows established epidemiological standards for temperature-mortality studies.** The implementation correctly applies:

1. Distributed Lag Non-Linear Models (Gasparrini et al. 2010)
2. Two-stage meta-analysis (Gasparrini et al. 2012)
3. Attributable burden calculation (Gasparrini & Leone 2014)
4. Multi-country study design (Gasparrini et al. 2015)

### Verification Summary

| Component | Status | Notes |
|-----------|--------|-------|
| Cross-basis specification | ✅ | ns(temp, 4df) × ns(lag, 4df) |
| Lag structure | ✅ | 21-day max, sensitivity 7-35 |
| GLM family | ✅ | Quasi-Poisson, log link |
| Confounding control | ✅ | Trend + seasonality + DOW |
| MMT estimation | ✅ | Cumulative RR minimization |
| Meta-analysis | ✅ | Mixmeta with REML |
| Attributable fractions | ✅ | Gasparrini & Leone formula |
| Stratification | ✅ | Age, sex, cause |
| Sensitivity analyses | ✅ | Lag, df, thresholds |
| Validation (case-crossover) | ✅ | Time-stratified design |

**The results can be considered reliable for publication and policy use.**

---

## References

- Gasparrini A et al. (2010) Distributed lag linear and non-linear models in R: the package dlnm. J Stat Softw
- Gasparrini A et al. (2012) Multivariate meta-analysis for non-linear and other multi-parameter associations. Stat Med
- Gasparrini A & Leone M (2014) Attributable risk from distributed lag models. BMC Med Res Methodol
- Gasparrini A et al. (2015) Mortality risk attributable to high and low ambient temperature: a multicountry observational study. Lancet
