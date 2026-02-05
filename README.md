# Temperature-Mortality Analysis: Brazil 2010-2024

## Overview

Analysis of temperature-mortality relationships among elderly population in Brazil using Distributed Lag Non-linear Models (DLNM) and multivariate meta-analysis.

**Geographic scales:** 509 Immediate Regions / 133 Intermediate Regions  
**Study period:** 2010-2024  
**Population:** Elderly (60+ years)

---

## Repository Structure

```
new_analysis/
├── phase0_data_prep/          # Step 1: Data preparation
│   ├── downloads/             # Scripts to download raw data
│   ├── aggregation/           # Aggregate to regional level
│   └── covariates/            # Process covariates
│
├── phase1_r/                  # Step 2: Core analysis
│   ├── 01_dlnm_analysis_v2.R  # Two-stage DLNM
│   ├── 01b_attributable_burden.R
│   ├── 01c_yll_calculation.R
│   ├── 01d_case_crossover.R
│   └── 01e_excess_mortality.R
│
├── phase2_r/                  # Step 3: Robustness checks
│   ├── 02a_sensitivity.R      # Sensitivity analyses
│   ├── 02b_harvesting.R       # Mortality displacement
│   └── 02c_heatwave.R         # Heatwave effects
│
├── phase3_r/                  # Step 4: Confounding
│   └── 03a_supplementary.R    # Time-varying confounders
│
└── phase4_r/                  # Step 5: Heterogeneity
    ├── 04a_meta_regression.R  # Effect modifiers
    ├── 04b_age_stratification.R
    ├── 04c_sex_stratification.R
    └── 04d_cause_stratification.R
```

---

## Reproduction Pipeline

### Requirements

```r
install.packages(c(
  "dlnm",        # Distributed lag models
  "mvmeta",      # Multivariate meta-analysis
  "splines",     # Natural splines
  "survival",    # Case-crossover
  "tidyverse",   # Data manipulation
  "arrow"        # Parquet files
))
```

### Step-by-Step Execution

```r
# Set working directory
setwd("path/to/new_analysis")

# PHASE 0: Data Preparation (Python scripts - run once)
# Download data from sources below, then run aggregation scripts

# PHASE 1: Core DLNM Analysis
source("phase1_r/01_dlnm_analysis_v2.R")      # ~2-4 hours
source("phase1_r/01b_attributable_burden.R")   # ~30 min
source("phase1_r/01c_yll_calculation.R")       # ~15 min
source("phase1_r/01d_case_crossover.R")        # ~1 hour
source("phase1_r/01e_excess_mortality.R")      # ~30 min

# PHASE 2: Robustness
source("phase2_r/02a_sensitivity.R")           # ~1 hour
source("phase2_r/02b_harvesting.R")            # ~2 hours
source("phase2_r/02c_heatwave.R")              # ~1 hour

# PHASE 3: Confounding
source("phase3_r/03a_supplementary.R")         # ~1 hour

# PHASE 4: Heterogeneity
source("phase4_r/04a_meta_regression.R")       # ~30 min
source("phase4_r/04b_age_stratification.R")    # ~2 hours
source("phase4_r/04c_sex_stratification.R")    # ~2 hours
source("phase4_r/04d_cause_stratification.R")  # ~2 hours
```

---

## Data Sources

### Mortality Data
| Source | DATASUS - Sistema de Informação sobre Mortalidade (SIM) |
|--------|----------------------------------------------------------|
| URL | https://datasus.saude.gov.br/transferencia-de-arquivos/ |
| Files | DO10OPEN.csv to DO24OPEN.csv |
| Variables | Daily deaths by municipality, age, sex, ICD-10 cause |

### Climate Data
| Source | ECMWF ERA5 Reanalysis |
|--------|----------------------|
| URL | https://cds.climate.copernicus.eu/ |
| Variables | 2m temperature (hourly → daily mean) |
| Resolution | 0.25° × 0.25° |

### Air Quality Data
| Source | Copernicus CAMS |
|--------|-----------------|
| URL | https://ads.atmosphere.copernicus.eu/ |
| Variables | PM2.5, PM10, O3, NO2 |

### Geographic Boundaries
| Source | IBGE |
|--------|------|
| URL | https://www.ibge.gov.br/geociencias/ |
| File | brazil_municipalities_2022.gpkg |

### Other Covariates
| Variable | Source |
|----------|--------|
| Life tables | IBGE (tabuas de mortalidade) |
| Holidays | Brazilian federal calendar |
| ENSO indices | NOAA Climate Prediction Center |

---

## Key Outputs

Results are saved to `phase*_r/results/` folders:

- `dlnm_r_*_summary.csv` - Temperature-mortality coefficients by region
- `attributable_burden_r_*.csv` - Attributable deaths and fractions
- `yll_r_*_regions.csv` - Years of life lost
- `excess_mortality_r_*.csv` - Annual excess deaths

---

## Citation

```bibtex
@article{temperature_mortality_brazil_2026,
  title={Temperature and Elderly Mortality in Brazil: 
         A Multi-Region Time-Series Analysis (2010-2024)},
  author={[Authors]},
  journal={[Journal]},
  year={2026}
}
```

## License

MIT License
