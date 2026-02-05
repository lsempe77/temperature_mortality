# Data Causal Analysis - README

## ⚠️ IMPORTANT: Data Source Documentation

This folder contains auxiliary data for causal analysis and robustness checks.
**Primary temperature data comes from ERA5 reanalysis, NOT from weather stations.**

---

## Folder Contents

### 1. `weather_stations/` - INMET Ground Station Data
**Purpose**: Validation of ERA5 reanalysis data, NOT primary analysis.

- Source: INMET (Instituto Nacional de Meteorologia)
- Years: 2010-2024
- Format: CSV files per station per year
- Naming: `INMET_{region}_{state}_{code}_{city}_{dates}.CSV`

**Usage**:
- Cross-validation of ERA5 temperature estimates
- Sensitivity analysis comparing ground vs. satellite data
- Quality control checks

**DO NOT USE FOR**:
- Primary DLNM analysis (use ERA5 instead)
- Spatial aggregation (coverage is incomplete)

### 2. `ocean_data/` - Climate Indices
**Purpose**: Control for large-scale climate variability.

| File | Index | Description |
|------|-------|-------------|
| `meiv2.csv` | MEI v2 | Multivariate ENSO Index |
| `oni.csv` | ONI | Oceanic Niño Index |
| `nina34.anom.csv` | Niño 3.4 | SST anomaly in Niño 3.4 region |
| `tna.csv` | TNA | Tropical North Atlantic SST |
| `tsa.csv` | TSA | Tropical South Atlantic SST |
| `whwp.csv` | WHWP | Western Hemisphere Warm Pool |

**Usage**:
- Confounding control in regression models
- Heterogeneity analysis by ENSO phase
- Long-term climate trend adjustment

### 3. `domicilios_atendidos.csv` - Household Coverage
**Purpose**: Socioeconomic covariate.

- Contains household access to services by region
- Can be used as SES control variable

---

## Analysis Pipeline Notes

### Primary Data Sources (in `phase0_data_prep/results/`):
- `era5_intermediate_daily.parquet` - ERA5 temperature at 133 intermediate regions
- `era5_immediate_daily.parquet` - ERA5 temperature at 510 immediate regions
- `mortality_*_daily_elderly.parquet` - SIM mortality data

### Scripts Updated (December 2025):

| Script | Version | Key Change |
|--------|---------|------------|
| `01a_*_dlnm_v2.py` | v2 | Natural spline cross-basis, MMT reference |
| `01d_attributable_burden_v2.py` | v2 | Auto-detects basis type, uses MMT |
| `01e_yll_unified.py` | unified | Merged actual/assumed modes, robust age parsing |

### Critical Dependencies:
1. **DLNM → Burden**: Basis type must match (ns vs poly)
2. **Burden → YLL**: File naming conventions must match (v2)
3. **Life tables**: Use IBGE Brazil-specific, not WHO fallback

---

## Quality Control Checklist

Before running analysis:
- [ ] ERA5 data covers study period (2010-2024)
- [ ] Mortality data parsed successfully (>95% parse rate)
- [ ] Life tables loaded from IBGE source
- [ ] MMT distribution is P75-P85 (expected for Brazil)
- [ ] Exclusion rate <20% for intermediate, <30% for immediate
- [ ] Death coverage >90%

---

## Contact

Analysis pipeline maintained by: Climate-Health Research Team
Last updated: December 2025
