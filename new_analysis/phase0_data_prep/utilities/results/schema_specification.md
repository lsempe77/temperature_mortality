# Unified Data Schema Specification

**Generated:** 2025-12-10 08:11:46

This document provides a formal schema specification for all datasets used in the
temperature-mortality analysis. Suitable for publication in supplementary materials.

---

## Dataset Overview

| Dataset | Category | Spatial Unit | Temporal Unit | Rows | Columns |
|---------|----------|--------------|---------------|------|---------|
| era5_intermediate_daily | Exposure | intermediate | daily | 728,707 | 6 |
| era5_immediate_daily | Exposure | immediate | daily | 2,794,290 | 6 |
| cams_intermediate_daily | Exposure | intermediate | daily | 728,707 | 5 |
| cams_immediate_daily | Exposure | immediate | daily | 2,794,290 | 5 |
| mortality_intermediate_daily | Outcome | intermediate | daily | 721,096 | 6 |
| mortality_intermediate_daily_elderly | Outcome | intermediate | daily | 708,767 | 6 |
| mortality_immediate_daily | Outcome | immediate | daily | 2,523,235 | 6 |
| mortality_immediate_daily_elderly | Outcome | immediate | daily | 2,302,539 | 6 |
| influenza_daily_by_intermediate_region | Confounder | intermediate | daily | 187,055 | 7 |
| influenza_daily_by_immediate_region | Confounder | immediate | daily | 378,575 | 7 |
| brazilian_holidays_daily | Confounder | national | daily | 5,479 | 10 |
| ses_intermediate_covariates | Covariate | intermediate | cross-sectional | 133 | 10 |
| ses_immediate_covariates | Covariate | immediate | cross-sectional | 510 | 10 |
| regional_covariates | Covariate | intermediate | cross-sectional | 133 | 13 |
| municipality_to_all_regions_map | Reference | municipality | cross-sectional | 5,571 | 7 |
| ibge_life_tables_combined | Reference | national | annual | 1,230 | 8 |
| yll_lookup_by_age | Reference | national | cross-sectional | 90 | 5 |

## Spatial Units

| Level | IBGE Name | Count | Description |
|-------|-----------|-------|-------------|
| Intermediate | Região Geográfica Intermediária | 133 | Larger regions for stable estimates |
| Immediate | Região Geográfica Imediata | 510 | Finer spatial resolution |
| National | Brasil | 1 | Country-level data |

## Variable Definitions

### Identifier Variables

| Variable | Description | Format |
|----------|-------------|--------|
| `intermediate_code` | Intermediate region identifier | Integer (IBGE code) |
| `immediate_code` | Immediate region identifier | Integer (IBGE code) |
| `region_code` | [Legacy] Intermediate region identifier | Integer (IBGE code) |
| `code_muni` | Municipality code | 7-digit integer |
| `date` | Date of observation | YYYY-MM-DD |

### Exposure Variables (Temperature)

| Variable | Unit | Source | Definition |
|----------|------|--------|------------|
| `temp_mean` | °C | ERA5 | Daily mean 2m temperature, population-weighted |
| `temp_min` | °C | ERA5 | Daily minimum 2m temperature |
| `temp_max` | °C | ERA5 | Daily maximum 2m temperature |
| `dewpoint_mean` | °C | ERA5 | Daily mean 2m dewpoint for humidity calculation |

### Exposure Variables (Pollution)

| Variable | Unit | Source | Definition |
|----------|------|--------|------------|
| `pm25_mean` | μg/m³ | CAMS | Daily mean PM2.5 concentration |
| `o3_mean` | μg/m³ | CAMS | Daily mean O3 (ozone) concentration |

### Outcome Variables (Mortality)

| Variable | Definition | ICD-10 Codes |
|----------|------------|--------------|
| `deaths_all` | All-cause mortality, all ages | All |
| `deaths_elderly` | All-cause mortality, age ≥60 | All |
| `deaths_respiratory` / `deaths_elderly_resp` | Respiratory mortality | J00-J99 |
| `deaths_cardiovascular` / `deaths_elderly_cvd` | Cardiovascular mortality | I00-I99 |
| `deaths_heat_direct` / `deaths_elderly_heat` | Direct heat-related mortality | T67, X30 |

### Confounder Variables

| Variable | Unit | Source | Definition |
|----------|------|--------|------------|
| `srag_cases` | Count | SIVEP-Gripe | Severe acute respiratory infection cases |
| `srag_influenza` | Count | SIVEP-Gripe | Laboratory-confirmed influenza |
| `srag_covid` | Count | SIVEP-Gripe | Laboratory-confirmed COVID-19 |
| `is_holiday` | Binary | holidays-br | National holiday indicator |
| `is_holiday_week` | Binary | Derived | Week contains a national holiday |

### Covariate Variables

| Variable | Unit | Source | Definition |
|----------|------|--------|------------|
| `pop_total` | Count | IBGE Census | Total population |
| `pop_elderly` | Count | IBGE Census | Population aged ≥60 |
| `pct_elderly` | % | Derived | Percentage of population ≥60 |
| `gdp_per_capita` | R$ | IBGE | Regional GDP per capita |
| `urbanization_rate` | % | IBGE | Percentage urban population |
| `ac_ownership` | % | PNAD | Air conditioning ownership rate |

### Reference Variables (Life Tables)

| Variable | Unit | Source | Definition |
|----------|------|--------|------------|
| `life_expectancy` | Years | IBGE | Remaining life expectancy at given age |
| `yll` | Years | Derived | Years of life lost if death at given age |

## Data Sources

| Source | Dataset | Time Period | URL |
|--------|---------|-------------|-----|
| ECMWF | ERA5 Reanalysis | 2010-2024 | https://cds.climate.copernicus.eu |
| ECMWF | CAMS Reanalysis | 2010-2024 | https://ads.atmosphere.copernicus.eu |
| DATASUS | SIM (Mortality) | 2010-2024 | https://datasus.saude.gov.br |
| DATASUS | SIVEP-Gripe | 2010-2024 | https://opendatasus.saude.gov.br |
| IBGE | Census/SIDRA | 2010-2022 | https://sidra.ibge.gov.br |
| IBGE | Life Tables | 2010-2024 | https://ibge.gov.br |
| IBGE | Geographic Regions | 2017 | https://ibge.gov.br/geociencias |
