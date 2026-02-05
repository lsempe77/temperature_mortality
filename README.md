# Temperature-Mortality Analysis: Brazil 2010-2024

[![DOI](https://img.shields.io/badge/DOI-pending-blue)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the analysis code for studying the relationship between temperature and elderly mortality in Brazil (2010-2024) using Distributed Lag Non-linear Models (DLNM) and multivariate meta-analysis.

**Key Features:**
- Two-stage DLNM with meta-analytic pooling
- Analysis at two spatial scales: 509 Immediate Geographic Regions and 133 Intermediate Geographic Regions
- Comprehensive robustness checks (sensitivity analyses, harvesting, heatwave effects)
- Heterogeneity analysis (age, sex, cause of death)
- Causal inference extensions (ENSO IV, Luz Para Todos DiD)

## Repository Structure

```
new_analysis/
├── phase0_data_prep/     # Data download and preprocessing scripts
├── phase1_core_model/    # Core DLNM analysis (Python)
├── phase1_r/             # Core DLNM analysis (R implementation)
├── phase1b_causal/       # Causal inference extensions
├── phase2_robustness/    # Sensitivity and harvesting analyses
├── phase2_r/             # R robustness checks
├── phase3_confounding/   # Supplementary confounding analyses
├── phase3_r/             # R supplementary analyses
├── phase4_heterogeneity/ # Effect modification analyses
├── phase4_r/             # R heterogeneity analyses
├── phase5_outputs/       # Figure and table generation
└── utils/                # Shared utility modules
```

## Data Sources

**This repository contains CODE ONLY.** Data must be obtained from the original sources:

### Mortality Data
- **Source:** Brazilian Ministry of Health - DATASUS
- **System:** Sistema de Informação sobre Mortalidade (SIM)
- **URL:** https://datasus.saude.gov.br/transferencia-de-arquivos/
- **Files:** DO10OPEN.csv to DO24OPEN.csv (2010-2024)
- **Variables:** Daily mortality counts by municipality, age, sex, ICD-10 cause

### Climate Data
- **Source:** ECMWF ERA5 Reanalysis
- **URL:** https://cds.climate.copernicus.eu/
- **Variables:** 2m temperature (hourly → daily mean)
- **Resolution:** 0.25° × 0.25°, aggregated to geographic regions
- **Period:** 2010-2024

### Air Quality Data
- **Source:** Copernicus Atmosphere Monitoring Service (CAMS)
- **URL:** https://ads.atmosphere.copernicus.eu/
- **Variables:** PM2.5, PM10, O3, NO2 (daily means)
- **Period:** 2010-2024

### Geographic Boundaries
- **Source:** IBGE (Brazilian Institute of Geography and Statistics)
- **URL:** https://www.ibge.gov.br/geociencias/
- **Files:** brazil_municipalities_2022.gpkg
- **Levels:** Municipalities → Immediate Regions → Intermediate Regions → States

### Covariates
| Variable | Source | URL |
|----------|--------|-----|
| Life tables | IBGE | https://www.ibge.gov.br/estatisticas/sociais/populacao/9126-tabuas-completas-de-mortalidade.html |
| Holidays | Brazilian federal calendar | - |
| ENSO indices | NOAA Climate Prediction Center | https://psl.noaa.gov/data/correlation/oni.data |
| Luz Para Todos | Ministry of Mines and Energy | https://www.gov.br/mme/ |
| AC ownership | PNAD Contínua | https://www.ibge.gov.br/estatisticas/sociais/trabalho/17270-pnad-continua.html |

## Requirements

### Python
```bash
pip install pandas numpy scipy statsmodels patsy xarray geopandas matplotlib seaborn
```

### R
```r
install.packages(c("dlnm", "mvmeta", "splines", "survival", "tidyverse", "arrow"))
```

## Reproducibility

1. **Download data** from sources listed above
2. **Place in `Input_data/` folder** with expected filenames
3. **Run Phase 0 scripts** to aggregate data to regional level
4. **Execute analyses** following the pipeline in `ANALYSIS_ROADMAP.md`

### Sample Data

Small sample datasets are included in `*_sample.csv` files to demonstrate data structure and allow code testing without full data downloads.

## Citation

If you use this code, please cite:

```bibtex
@article{temperature_mortality_brazil_2026,
  title={Temperature-Mortality Relationships in Brazil: A Multi-Region Analysis (2010-2024)},
  author={[Authors]},
  journal={[Journal]},
  year={2026}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions about the code or analysis, please open an issue on this repository.
