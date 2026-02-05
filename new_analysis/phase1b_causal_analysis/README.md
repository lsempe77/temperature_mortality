# Phase 1b: Causal Analysis

## Overview

This phase implements quasi-experimental designs to estimate causal effects of temperature on mortality, leveraging two sources of exogenous variation:

1. **ENSO (El Niño Southern Oscillation)** - Oceanic climate phenomenon providing exogenous temperature shocks
2. **Luz para Todos** - Brazil's rural electrification program with staggered rollout

## Data Preparation Scripts

Run these first to prepare the analysis panels:

```bash
# 1. Prepare ENSO indices (ONI, MEI, Niño 3.4)
python data_prep/01_prepare_enso_data.py

# 2. Prepare Luz para Todos treatment data
python data_prep/02_prepare_lpt_data.py

# 3. Merge all data into analysis panel
python data_prep/03_merge_causal_panel.py --level both
```

## Analysis Scripts

### ENSO-Based Designs

| Script | Method | Description |
|--------|--------|-------------|
| `01b_enso_iv_analysis.py` | Instrumental Variables | Uses ENSO as instrument for temperature |
| `01b_enso_event_study.py` | Event Study | Dynamic effects around El Niño/La Niña events |

### Luz para Todos Designs

| Script | Method | Description |
|--------|--------|-------------|
| `01b_lpt_staggered_did.py` | Staggered DiD | TWFE, Callaway-Sant'Anna, Triple Diff |
| `01b_lpt_event_study.py` | Dynamic Effects | Pre-trends and post-treatment dynamics |

### Combined Designs

| Script | Method | Description |
|--------|--------|-------------|
| `01b_interaction_analysis.py` | Interaction | ENSO × Electrification effects |
| `01b_placebo_tests.py` | Validation | Falsification and balance tests |

## Running All Analyses

```bash
# Run for both intermediate (133 regions) and immediate (510 regions) levels

# ENSO analyses
python 01b_enso_iv_analysis.py --level both
python 01b_enso_event_study.py --level both

# Luz para Todos analyses
python 01b_lpt_staggered_did.py --level both
python 01b_lpt_event_study.py --level both

# Interaction and validation
python 01b_interaction_analysis.py --level both
python 01b_placebo_tests.py --level both
```

## Identification Strategies

### 1. ENSO as Instrumental Variable

- **First Stage**: ENSO → Regional temperature anomalies
- **Second Stage**: Predicted temperature → Mortality
- **Exclusion Restriction**: ENSO affects mortality only through temperature
- **Strength**: F-statistic > 10 for strong instrument

### 2. Luz para Todos Staggered DiD

- **Treatment**: First electrification in region (enables AC access)
- **Identification**: Staggered rollout across 5,000+ municipalities (2004-2024)
- **Estimators**: 
  - TWFE (traditional, may be biased)
  - Callaway-Sant'Anna (robust to heterogeneous effects)
  - Triple Difference (electricity × hot days × post-treatment)

### 3. ENSO × Electrification Interaction

- **Key Question**: Does electrification attenuate ENSO-induced mortality?
- **Interpretation**: If β₃ < 0, electrification protects against climate shocks
- **Mechanism**: AC access as adaptation to heat

## Output Files

Results saved to `results/` folder:

| File | Content |
|------|---------|
| `enso_monthly.parquet` | Monthly ENSO indices with phase classification |
| `enso_daily.parquet` | Daily interpolated ENSO for merging |
| `lpt_*_treatment.parquet` | Treatment timing by region |
| `causal_panel_*.parquet` | Merged analysis panel |
| `*_results.json` | Analysis results |
| `*_results_immediate.json` | Immediate level results |

## Key Econometric Considerations

1. **Staggered DiD**: Use Callaway-Sant'Anna over TWFE to avoid negative weighting
2. **Standard Errors**: Clustered at region level for spatial correlation
3. **Pre-trends**: Critical for DiD validity - test in event study
4. **Multiple Testing**: Interpret results as a bundle, not individually
5. **Effect Heterogeneity**: Stratify by temperature exposure, urbanization

## References

- Callaway, B., & Sant'Anna, P. H. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225(2), 200-230.
- de Chaisemartin, C., & D'Haultfoeuille, X. (2020). Two-way fixed effects estimators with heterogeneous treatment effects. *American Economic Review*, 110(9), 2964-96.
