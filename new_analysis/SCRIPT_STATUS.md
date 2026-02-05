# Script Status & Execution Guide

**Last Updated:** December 17, 2025 09:00

---

## Execution Status (Dec 17, 2025)

### Phase 1: Core DLNM ✅ ALL COMPLETE
| Script | Intermediate | Immediate | Last Run |
|--------|--------------|-----------|----------|
| `01_dlnm_analysis_v2.py` | ✅ Complete | ✅ Complete | Dec 16, 00:35 |
| `01d_attributable_burden_v2.py` | ✅ Complete | - | Dec 16, 01:36 |
| `01e_yll_unified.py` | ✅ Complete | ✅ Complete | Dec 16, 08:25 |
| `01f_case_crossover_v2.py` | ✅ Complete | - | Dec 16, 01:22 |
| `01g_excess_mortality_v2.py` | ✅ Complete | - | Dec 16, 00:59 |

### Phase 2: Robustness
| Script | Intermediate | Immediate | Last Run | Notes |
|--------|--------------|-----------|----------|-------|
| `02a_sensitivity_analyses_v2.py` | ✅ Complete | ⚠️ Old | Dec 17, 03:43 / Dec 12 | Immediate needs re-run |
| `02b_harvesting_analysis_v2.py` | ✅ Complete | ✅ Complete | Dec 16, 19:45 / 20:30 | |
| `02c_heatwave_dlnm_v2.py` | ✅ Complete | ✅ Complete | Dec 16, 01:28 / 02:00 | |

### Phase 3: Confounding ✅ ALL COMPLETE
| Script | Intermediate | Immediate | Last Run |
|--------|--------------|-----------|----------|
| `03a_supplementary_analyses_v2.py` | ✅ Complete | ✅ Complete | Dec 17, 01:35 / 06:28 |

### Phase 4: Heterogeneity ⚠️ ISSUES WITH MVMETA CIs
| Script | Intermediate | Immediate | Last Run | Notes |
|--------|--------------|-----------|----------|-------|
| `04a_meta_regression_v2.py` | ⚠️ NaN RRs | ⚠️ NaN RRs | Dec 16, 21:18 / 22:10 | MVMeta output has NaN |
| `04b_age_stratification_v2.py` | ⚠️ Empty pooled | ⚠️ Empty pooled | Dec 16, 22:16 / 23:38 | Pooled results empty |
| `04c_sex_stratification_v2.py` | ⚠️ Inf CIs | ⚠️ Inf CIs | Dec 16, 21:53 / 23:05 | RRs computed but CIs=Inf |
| `04d_cause_stratification_v2.py` | ⚠️ Old | ⚠️ Old | Dec 15, 14:00 / 12:22 | Pre-MVMeta (uses simple meta) |

**Legend:** ✅ Complete & Valid | ⚠️ Issues/Old | ⏳ Running | ❌ Error

---

## Current Issues (Dec 17, 2025)

### MVMeta Confidence Interval Bug
**Symptom:** MVMeta produces valid RRs but `ci_low=0.0`, `ci_high=Infinity`, `se=Infinity`
**Affected scripts:** 02a, 04a, 04b, 04c, 04d
**Status:** Under investigation - variance computation may have numerical issues

### 04a_meta_regression: NaN in pooled results
- Region-level fits work (132/133 regions)
- MVMeta pooling produces NaN RRs
- May be related to `crossbasis_info` structure

### 04b_age_stratification: Empty pooled results
- All age groups fitted successfully (132-133 regions each)
- MVMeta pooling section produced empty results
- Import was missing `mvmeta_pool_coefficients` - FIXED

---

## Scripts That Are Valid (Use Per-Region Prediction)

These scripts use `predict_cumulative_rr` per-region and are NOT affected by MVMeta bugs:
- ✅ Phase 1: 01a, 01d, 01e, 01f, 01g
- ✅ Phase 2: 02b, 02c  
- ✅ Phase 3: 03a (region-level results valid, pooled may have CI issues)

---

## Execution Commands

### Phase 1: Core DLNM
```powershell
# Intermediate (133 regions) - PRIMARY
cd "c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis\phase1_core_model"
python 01a_intermediate_dlnm_v2.py

# Immediate (510 regions) - SECONDARY  
python 01a_immediate_dlnm_v2.py
```

### Phase 2: Robustness
```powershell
cd "c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis\phase2_robustness"

# Sensitivity analyses (tests lag/df/percentile specs)
python 02a_sensitivity_analyses_v2.py --level intermediate
python 02a_sensitivity_analyses_v2.py --level immediate

# Harvesting/mortality displacement
python 02b_harvesting_analysis_v2.py --level intermediate
python 02b_harvesting_analysis_v2.py --level immediate

# Heatwave interaction effects
python 02c_heatwave_dlnm_v2.py --level intermediate
python 02c_heatwave_dlnm_v2.py --level immediate
```

### Phase 3: Confounding Control
```powershell
cd "c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis\phase3_confounding"

# Supplementary analyses (apparent temp, pollution, flu adjustment)
python 03a_supplementary_analyses_v2.py --level intermediate
python 03a_supplementary_analyses_v2.py --level immediate
```

### Phase 4: Heterogeneity
```powershell
cd "c:\Users\LucasSempe\OneDrive - International Initiative for Impact Evaluation\Desktop\sim_data\new_analysis\phase4_heterogeneity"

# Meta-regression (regional moderators)
python 04a_meta_regression_v2.py --level intermediate
python 04a_meta_regression_v2.py --level immediate

# Age stratification (60-69, 70-79, 80+)
python 04b_age_stratification_v2.py --level intermediate
python 04b_age_stratification_v2.py --level immediate

# Sex stratification (Male, Female)
python 04c_sex_stratification_v2.py --level intermediate
python 04c_sex_stratification_v2.py --level immediate

# Cause stratification (Cardiovascular, Respiratory, External)
python 04d_cause_stratification_v2.py --level intermediate
python 04d_cause_stratification_v2.py --level immediate
```

---

## Data Flow

```
Input_data/
├── DO10OPEN.csv ... DO24OPEN.csv  (raw SIM mortality - semicolon EXCEPT DO24 uses comma)
├── era5_brazil_hourly_*.nc         (raw ERA5 climate)
└── brazil_municipalities_2022.gpkg (boundaries)
        │
        ▼ [Phase 0 aggregation scripts]
        
phase0_data_prep/results/
├── era5_intermediate_daily.parquet   (temp by 133 regions)
├── era5_immediate_daily.parquet      (temp by 510 regions)
├── mortality_regional_daily_elderly.parquet
├── mortality_immediate_daily_elderly.parquet
├── mortality_by_cause_intermediate.parquet
├── mortality_by_cause_immediate.parquet
├── cams_intermediate_daily.parquet   (pollution)
├── cams_immediate_daily.parquet
└── ...
        │
        ▼ [Phase 1-4 analysis scripts]
        
phase*/results/
├── dlnm_results_intermediate.json
├── dlnm_results_immediate.json
└── ...
```

---

## Scripts That Read Raw SIM Files

These scripts read directly from `Input_data/DO*OPEN.csv`:

| Script | Columns Used | Separator |
|--------|--------------|-----------|
| `phase0_data_prep/aggregation/00n_aggregate_mortality_to_regions.py` | DTOBITO, CODMUNRES, IDADE, CAUSABAS, SEXO | Auto-detect ✅ |
| `phase1_core_model/01e_yll_unified.py` | IDADE | Auto-detect ✅ |
| `phase1_core_model/01f_case_crossover_v2.py` | DTOBITO, CODMUNRES, IDADE, CAUSABAS | Auto-detect ✅ |
| `phase4_heterogeneity/04b_age_stratification_v2.py` | DTOBITO, IDADE, CODMUNRES | Auto-detect ✅ |
| `phase4_heterogeneity/04c_sex_stratification_v2.py` | DTOBITO, IDADE, SEXO, CODMUNRES | Auto-detect ✅ |

---

## Output Files by Phase

### Phase 1 Core
- `phase1_core_model/results/dlnm_intermediate_mvmeta_results.json`
- `phase1_core_model/results/dlnm_immediate_mvmeta_results.json`
- `phase1_core_model/results/attributable_burden_intermediate.json`
- `phase1_core_model/results/attributable_burden_immediate.json`

### Phase 2 Robustness
- `phase2_robustness/results/sensitivity_analysis_intermediate.json`
- `phase2_robustness/results/sensitivity_analysis_immediate.json`
- `phase2_robustness/results/harvesting_analysis_intermediate.json`
- `phase2_robustness/results/heatwave_dlnm_intermediate.json`

### Phase 3 Confounding
- `phase3_confounding/results/supplementary_analyses_intermediate.json`
- `phase3_confounding/results/supplementary_analyses_immediate.json`

### Phase 4 Heterogeneity
- `phase4_heterogeneity/results/meta_regression_intermediate.json`
- `phase4_heterogeneity/results/age_stratification_intermediate.json`
- `phase4_heterogeneity/results/sex_stratification_intermediate.json`
- `phase4_heterogeneity/results/cause_stratification_intermediate.json`
- (Same pattern for `_immediate`)

---

## Key Parameters (All v2 Scripts)

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Max lag | 21 days | Capture delayed mortality |
| Temp df | 4 (natural spline) | Gasparrini standard |
| Lag df | 4 (natural spline) | Gasparrini standard |
| Heat threshold | P97.5 (primary), P99 (sensitivity) | Literature standard |
| Cold threshold | P2.5 (primary), P1 (sensitivity) | Literature standard |
| Reference | MMT (region-specific) | Minimum mortality temperature |
| Meta-analysis | MVMeta (coefficient pooling) | Avoids boundary MMT issues |

---

## Troubleshooting

### "Usecols do not match columns"
- **Cause:** DO24OPEN.csv uses comma separator, script expected semicolon
- **Fix:** Added auto-detect separator in 04b, 04c, 01e (Dec 16, 2025)

### "mvmeta_pool_coefficients not found"
- **Cause:** Function was named differently in dlnm_module.py
- **Fix:** Added wrapper functions to dlnm_module.py (Dec 15, 2025)

### RR = 1.0 for P99 heat
- **Cause:** MMT at boundary → SE≈0 → infinite weight in meta-analysis
- **Fix:** Use MVMeta coefficient pooling instead of RR pooling
