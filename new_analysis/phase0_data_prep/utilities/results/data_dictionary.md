# Phase 0 Data Dictionary

**Generated:** 2025-12-10 08:11:36

This document provides column-by-column documentation for all Phase 0 aggregated datasets.

---

## Exposure Data

### era5_intermediate_daily

**Description:** ERA5 temperature data aggregated to 133 intermediate regions

**File:** `era5_intermediate_daily.parquet`

- **Rows:** 728,707
- **Columns:** 6
- **Spatial Unit:** intermediate
- **Temporal Unit:** daily
- **Memory:** 46.6 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `date` | object | 728,707 | 0.0% | Date of observation | e.g., 2010-01-01, 2010-01-02 |
| `region_code` | int64 | 728,707 | 0.0% | Intermediate region identifier (IBGE code) | Range: [1101.00, 5301.00], Mean: 3086.54 |
| `temp_mean` | float32 | 728,707 | 0.0% | Daily mean temperature (°C) | Range: [0.99, 33.91], Mean: 23.99 |
| `temp_min` | float32 | 728,707 | 0.0% | Daily minimum temperature (°C) | Range: [-4.31, 28.16], Mean: 20.40 |
| `temp_max` | float32 | 728,707 | 0.0% | Daily maximum temperature (°C) | Range: [3.91, 41.69], Mean: 28.36 |
| `dewpoint_mean` | float32 | 728,707 | 0.0% | Daily mean dewpoint temperature (°C) | Range: [-6.72, 25.74], Mean: 18.18 |

### era5_immediate_daily

**Description:** ERA5 temperature data aggregated to 510 immediate regions

**File:** `era5_immediate_daily.parquet`

- **Rows:** 2,794,290
- **Columns:** 6
- **Spatial Unit:** immediate
- **Temporal Unit:** daily
- **Memory:** 178.8 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `date` | object | 2,794,290 | 0.0% | Date of observation | e.g., 2010-01-01, 2010-01-02 |
| `region_code` | int64 | 2,794,290 | 0.0% | Immediate region identifier (IBGE code) | Range: [110001.00, 530001.00], Mean: 318584.93 |
| `temp_mean` | float32 | 2,794,290 | 0.0% | Daily mean temperature (°C) | Range: [0.67, 35.14], Mean: 23.72 |
| `temp_min` | float32 | 2,794,290 | 0.0% | Daily minimum temperature (°C) | Range: [-4.77, 30.17], Mean: 19.86 |
| `temp_max` | float32 | 2,794,290 | 0.0% | Daily maximum temperature (°C) | Range: [3.46, 42.70], Mean: 28.45 |
| `dewpoint_mean` | float32 | 2,794,290 | 0.0% | Daily mean dewpoint temperature (°C) | Range: [-8.74, 26.45], Mean: 17.98 |

### cams_intermediate_daily

**Description:** CAMS pollution data aggregated to 133 intermediate regions

**File:** `cams_intermediate_daily.parquet`

- **Rows:** 728,707
- **Columns:** 5
- **Spatial Unit:** intermediate
- **Temporal Unit:** daily
- **Memory:** 20.4 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `intermediate_code` | int64 | 728,707 | 0.0% | Intermediate region identifier | Range: [1101.00, 5301.00], Mean: 3086.54 |
| `date` | datetime64[ns] | 728,707 | 0.0% | Date of observation | 2010-01-01 to 2024-12-31 |
| `pm25` | float32 | 728,707 | 0.0% | Daily mean PM2.5 concentration (μg/m³) | Range: [0.33, 948.62], Mean: 10.90 |
| `pm10` | float32 | 728,707 | 0.0% |  | Range: [0.47, 1288.90], Mean: 15.47 |
| `ozone` | float32 | 728,707 | 0.0% | Daily mean O3 concentration (μg/m³) | Range: [223.69, 370.01], Mean: 264.57 |

### cams_immediate_daily

**Description:** CAMS pollution data aggregated to 510 immediate regions

**File:** `cams_immediate_daily.parquet`

- **Rows:** 2,794,290
- **Columns:** 5
- **Spatial Unit:** immediate
- **Temporal Unit:** daily
- **Memory:** 78.2 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `immediate_code` | int64 | 2,794,290 | 0.0% | Immediate region identifier | Range: [110001.00, 530001.00], Mean: 318584.93 |
| `date` | datetime64[ns] | 2,794,290 | 0.0% | Date of observation | 2010-01-01 to 2024-12-31 |
| `pm25` | float32 | 2,794,290 | 0.0% | Daily mean PM2.5 concentration (μg/m³) | Range: [0.21, 1398.42], Mean: 10.71 |
| `pm10` | float32 | 2,794,290 | 0.0% |  | Range: [0.29, 1850.12], Mean: 15.23 |
| `ozone` | float32 | 2,794,290 | 0.0% | Daily mean O3 concentration (μg/m³) | Range: [222.77, 371.92], Mean: 264.94 |

## Outcome Data

### mortality_intermediate_daily

**Description:** All-ages mortality aggregated to 133 intermediate regions

**File:** `mortality_intermediate_daily.parquet`

- **Rows:** 721,096
- **Columns:** 6
- **Spatial Unit:** intermediate
- **Temporal Unit:** daily
- **Memory:** 57.7 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `date` | object | 721,096 | 0.0% | Date of death | e.g., 2010-01-01, 2010-01-02 |
| `intermediate_code` | int64 | 721,096 | 0.0% | Intermediate region identifier | Range: [1101.00, 5301.00], Mean: 3100.99 |
| `deaths_all` | int64 | 721,096 | 0.0% | Total deaths (all ages) | Range: [1.00, 1106.00], Mean: 28.23 |
| `deaths_respiratory` | int64 | 721,096 | 0.0% | Respiratory deaths (ICD J00-J99) | Range: [0.00, 137.00], Mean: 3.13 |
| `deaths_cardiovascular` | int64 | 721,096 | 0.0% | Cardiovascular deaths (ICD I00-I99) | Range: [0.00, 206.00], Mean: 7.47 |
| `deaths_heat_direct` | int64 | 721,096 | 0.0% | Direct heat deaths (ICD T67, X30) | Range: [0.00, 1.00], Mean: 0.00 |

### mortality_intermediate_daily_elderly

**Description:** Elderly (60+) mortality aggregated to 133 intermediate regions

**File:** `mortality_intermediate_daily_elderly.parquet`

- **Rows:** 708,767
- **Columns:** 6
- **Spatial Unit:** intermediate
- **Temporal Unit:** daily
- **Memory:** 56.7 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `date` | object | 708,767 | 0.0% | Date of death | e.g., 2010-01-01, 2010-01-02 |
| `intermediate_code` | int64 | 708,767 | 0.0% | Intermediate region identifier | Range: [1101.00, 5301.00], Mean: 3120.50 |
| `deaths_elderly` | int64 | 708,767 | 0.0% | Total elderly deaths (60+) | Range: [1.00, 775.00], Mean: 19.30 |
| `deaths_elderly_resp` | int64 | 708,767 | 0.0% | Elderly respiratory deaths | Range: [0.00, 106.00], Mean: 2.66 |
| `deaths_elderly_cvd` | int64 | 708,767 | 0.0% | Elderly cardiovascular deaths | Range: [0.00, 172.00], Mean: 6.07 |
| `deaths_elderly_heat` | int64 | 708,767 | 0.0% | Elderly direct heat deaths | Range: [0.00, 1.00], Mean: 0.00 |

### mortality_immediate_daily

**Description:** All-ages mortality aggregated to 510 immediate regions

**File:** `mortality_immediate_daily.parquet`

- **Rows:** 2,523,235
- **Columns:** 6
- **Spatial Unit:** immediate
- **Temporal Unit:** daily
- **Memory:** 201.9 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `date` | object | 2,523,235 | 0.0% | Date of death | e.g., 2010-01-01, 2010-01-02 |
| `immediate_code` | int64 | 2,523,235 | 0.0% | Immediate region identifier | Range: [110001.00, 530001.00], Mean: 320402.62 |
| `deaths_all` | int64 | 2,523,235 | 0.0% | Total deaths (all ages) | Range: [1.00, 1007.00], Mean: 8.07 |
| `deaths_respiratory` | int64 | 2,523,235 | 0.0% | Respiratory deaths (ICD J00-J99) | Range: [0.00, 132.00], Mean: 0.90 |
| `deaths_cardiovascular` | int64 | 2,523,235 | 0.0% | Cardiovascular deaths (ICD I00-I99) | Range: [0.00, 190.00], Mean: 2.14 |
| `deaths_heat_direct` | int64 | 2,523,235 | 0.0% | Direct heat deaths (ICD T67, X30) | Range: [0.00, 1.00], Mean: 0.00 |

### mortality_immediate_daily_elderly

**Description:** Elderly (60+) mortality aggregated to 510 immediate regions

**File:** `mortality_immediate_daily_elderly.parquet`

- **Rows:** 2,302,539
- **Columns:** 6
- **Spatial Unit:** immediate
- **Temporal Unit:** daily
- **Memory:** 184.2 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `date` | object | 2,302,539 | 0.0% | Date of death | e.g., 2010-01-01, 2010-01-02 |
| `immediate_code` | int64 | 2,302,539 | 0.0% | Immediate region identifier | Range: [110001.00, 530001.00], Mean: 322359.14 |
| `deaths_elderly` | int64 | 2,302,539 | 0.0% | Total elderly deaths (60+) | Range: [1.00, 692.00], Mean: 5.94 |
| `deaths_elderly_resp` | int64 | 2,302,539 | 0.0% | Elderly respiratory deaths | Range: [0.00, 99.00], Mean: 0.82 |
| `deaths_elderly_cvd` | int64 | 2,302,539 | 0.0% | Elderly cardiovascular deaths | Range: [0.00, 154.00], Mean: 1.87 |
| `deaths_elderly_heat` | int64 | 2,302,539 | 0.0% | Elderly direct heat deaths | Range: [0.00, 1.00], Mean: 0.00 |

## Confounder Data

### influenza_daily_by_intermediate_region

**Description:** SRAG/Influenza cases aggregated to 133 intermediate regions

**File:** `influenza_daily_by_intermediate_region.parquet`

- **Rows:** 187,055
- **Columns:** 7
- **Spatial Unit:** intermediate
- **Temporal Unit:** daily
- **Memory:** 21.3 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `date` | datetime64[ns] | 187,055 | 0.0% | Date of notification | 2010-01-03 to 2024-12-28 |
| `intermediate_code` | int64 | 187,055 | 0.0% | Intermediate region identifier | Range: [1101.00, 5301.00], Mean: 3286.82 |
| `intermediate_name` | object | 187,055 | 0.0% |  | e.g., Porto Velho, Ji-Paraná |
| `srag_cases` | int64 | 187,055 | 0.0% | Total SRAG cases (elderly) | Range: [1.00, 1099.00], Mean: 9.74 |
| `srag_influenza` | int64 | 187,055 | 0.0% | Confirmed influenza cases | Range: [0.00, 1099.00], Mean: 9.56 |
| `srag_covid` | int64 | 187,055 | 0.0% | Confirmed COVID-19 cases | Range: [0.00, 0.00], Mean: 0.00 |
| `srag_deaths` | int64 | 187,055 | 0.0% | SRAG-related deaths | Range: [0.00, 12.00], Mean: 0.06 |

### influenza_daily_by_immediate_region

**Description:** SRAG/Influenza cases aggregated to 510 immediate regions

**File:** `influenza_daily_by_immediate_region.parquet`

- **Rows:** 378,575
- **Columns:** 7
- **Spatial Unit:** immediate
- **Temporal Unit:** daily
- **Memory:** 43.8 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `date` | datetime64[ns] | 378,575 | 0.0% | Date of notification | 2010-01-03 to 2024-12-28 |
| `immediate_code` | int64 | 378,575 | 0.0% | Immediate region identifier | Range: [110001.00, 530001.00], Mean: 334877.51 |
| `immediate_name` | object | 378,575 | 0.0% |  | e.g., Porto Velho, Ariquemes |
| `srag_cases` | int64 | 378,575 | 0.0% | Total SRAG cases (elderly) | Range: [1.00, 1039.00], Mean: 4.81 |
| `srag_influenza` | int64 | 378,575 | 0.0% | Confirmed influenza cases | Range: [0.00, 1039.00], Mean: 4.72 |
| `srag_covid` | int64 | 378,575 | 0.0% | Confirmed COVID-19 cases | Range: [0.00, 0.00], Mean: 0.00 |
| `srag_deaths` | int64 | 378,575 | 0.0% | SRAG-related deaths | Range: [0.00, 12.00], Mean: 0.03 |

### brazilian_holidays_daily

**Description:** Brazilian national holidays 2010-2024

**File:** `brazilian_holidays_daily.parquet`

- **Rows:** 5,479
- **Columns:** 10
- **Spatial Unit:** national
- **Temporal Unit:** daily
- **Memory:** 0.6 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `date` | datetime64[ns] | 5,479 | 0.0% | Date | 2010-01-01 to 2024-12-31 |
| `holiday_name_pt` | object | 196 | 96.4% | Name of holiday in Portuguese (if applicable) | e.g., Universal Fraternization Day, Sábado de Carnaval |
| `holiday_type` | object | 196 | 96.4% | Type of holiday | e.g., civic, carnival |
| `is_holiday` | int64 | 5,479 | 0.0% | Binary indicator for national holiday | Range: [0.00, 1.00], Mean: 0.04 |
| `is_day_before_holiday` | int64 | 5,479 | 0.0% |  | Range: [0.00, 1.00], Mean: 0.04 |
| `is_day_after_holiday` | int64 | 5,479 | 0.0% |  | Range: [0.00, 1.00], Mean: 0.04 |
| `is_holiday_week` | int64 | 5,479 | 0.0% |  | Range: [0.00, 1.00], Mean: 0.09 |
| `year` | int32 | 5,479 | 0.0% |  | Range: [2010.00, 2024.00], Mean: 2017.00 |
| `month` | int32 | 5,479 | 0.0% |  | Range: [1.00, 12.00], Mean: 6.52 |
| `day_of_week` | int32 | 5,479 | 0.0% |  | Range: [0.00, 6.00], Mean: 3.00 |

## Covariate Data

### ses_intermediate_covariates

**Description:** Socioeconomic covariates for 133 intermediate regions

**File:** `ses_intermediate_covariates.csv`

- **Rows:** 133
- **Columns:** 10
- **Spatial Unit:** intermediate
- **Temporal Unit:** cross-sectional
- **Memory:** 0.0 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `intermediate_code` | int64 | 133 | 0.0% | Intermediate region identifier | Range: [1101.00, 5301.00], Mean: 3086.54 |
| `intermediate_name` | object | 133 | 0.0% | Region name | e.g., Porto Velho, Ji-Paraná |
| `pop_total` | int64 | 133 | 0.0% | Total population | Range: [90456.00, 22564260.00], Mean: 1526922.98 |
| `pop_elderly` | float64 | 133 | 0.0% | Elderly population (60+) | Range: [8168.00, 3742255.00], Mean: 241454.81 |
| `gdp_total` | int64 | 133 | 0.0% | Total GDP (R$) | Range: [2137659.00, 1469758150.00], Mean: 67760466.24 |
| `pop_urban` | float64 | 133 | 0.0% | Urban population | Range: [38195.00, 21136475.00], Mean: 1209968.45 |
| `n_municipalities` | int64 | 133 | 0.0% |  | Range: [1.00, 146.00], Mean: 41.88 |
| `pct_elderly` | float64 | 133 | 0.0% | Percentage elderly | Range: [6.65, 22.53], Mean: 15.15 |
| `gdp_per_capita` | float64 | 133 | 0.0% |  | Range: [9939.45, 101847.70], Mean: 36935.44 |
| `urbanization_rate` | float64 | 133 | 0.0% |  | Range: [37.95, 100.00], Mean: 70.76 |

### ses_immediate_covariates

**Description:** Socioeconomic covariates for 510 immediate regions

**File:** `ses_immediate_covariates.csv`

- **Rows:** 510
- **Columns:** 10
- **Spatial Unit:** immediate
- **Temporal Unit:** cross-sectional
- **Memory:** 0.1 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `immediate_code` | int64 | 510 | 0.0% | Immediate region identifier | Range: [110001.00, 530001.00], Mean: 318584.93 |
| `immediate_name` | object | 510 | 0.0% | Region name | e.g., Porto Velho, Ariquemes |
| `pop_total` | int64 | 510 | 0.0% | Total population | Range: [30980.00, 20731920.00], Mean: 398197.56 |
| `pop_elderly` | float64 | 510 | 0.0% | Elderly population (60+) | Range: [3090.00, 3391059.00], Mean: 62967.63 |
| `gdp_total` | int64 | 510 | 0.0% | Total GDP (R$) | Range: [488502.00, 1390102770.00], Mean: 17670866.69 |
| `pop_urban` | float64 | 510 | 0.0% | Urban population | Range: [9182.00, 19458888.00], Mean: 315540.79 |
| `n_municipalities` | int64 | 510 | 0.0% |  | Range: [1.00, 47.00], Mean: 10.92 |
| `pct_elderly` | float64 | 510 | 0.0% | Percentage elderly | Range: [5.66, 25.83], Mean: 15.97 |
| `gdp_per_capita` | float64 | 510 | 0.0% |  | Range: [7979.79, 226369.29], Mean: 36059.08 |
| `urbanization_rate` | float64 | 510 | 0.0% |  | Range: [15.06, 100.00], Mean: 69.21 |

### regional_covariates

**Description:** Additional regional covariates for meta-regression

**File:** `regional_covariates.csv`

- **Rows:** 133
- **Columns:** 13
- **Spatial Unit:** intermediate
- **Temporal Unit:** cross-sectional
- **Memory:** 0.0 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `region_code` | int64 | 133 | 0.0% | Intermediate region identifier | Range: [1101.00, 5301.00], Mean: 3086.54 |
| `region_name` | object | 133 | 0.0% | Region name | e.g., Ji-Paraná, Porto Velho |
| `primary_state` | object | 133 | 0.0% |  | e.g., RO, AC |
| `macro_region` | object | 133 | 0.0% |  | e.g., North, Northeast |
| `n_municipalities` | int64 | 133 | 0.0% |  | Range: [1.00, 146.00], Mean: 41.89 |
| `ac_pct` | float64 | 133 | 0.0% | Air conditioning ownership rate (%) | Range: [15.20, 42.10], Mean: 26.30 |
| `hdi` | float64 | 133 | 0.0% | Human Development Index | Range: [0.68, 0.85], Mean: 0.75 |
| `gdp_per_capita_brl` | int64 | 133 | 0.0% | GDP per capita (R$) | Range: [14739.00, 90742.00], Mean: 32185.71 |
| `elderly_pct` | float64 | 133 | 0.0% |  | Range: [7.40, 18.50], Mean: 13.33 |
| `urban_pct` | float64 | 133 | 0.0% |  | Range: [65.30, 97.80], Mean: 83.36 |
| `hospital_beds_per_1000` | float64 | 133 | 0.0% |  | Range: [1.30, 2.80], Mean: 1.91 |
| `physicians_per_1000` | float64 | 133 | 0.0% |  | Range: [0.80, 4.50], Mean: 1.76 |
| `mean_temp_annual` | float64 | 133 | 0.0% |  | Range: [18.20, 27.80], Mean: 23.75 |

## Reference Data

### municipality_to_all_regions_map

**Description:** Municipality to intermediate/immediate region mapping

**File:** `municipality_to_all_regions_map.csv`

- **Rows:** 5,571
- **Columns:** 7
- **Spatial Unit:** municipality
- **Temporal Unit:** cross-sectional
- **Memory:** 1.5 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `code_muni` | int64 | 5,571 | 0.0% | Municipality code (7-digit IBGE) | Range: [1100015.00, 5300108.00], Mean: 3253922.53 |
| `name_muni` | object | 5,571 | 0.0% | Municipality name | e.g., Alta Floresta D'Oeste, Ariquemes |
| `abbrev_state` | object | 5,571 | 0.0% | State abbreviation | e.g., RO, AC |
| `intermediate_code` | int64 | 5,571 | 0.0% | Intermediate region code | Range: [1101.00, 5301.00], Mean: 3242.19 |
| `intermediate_name` | object | 5,571 | 0.0% | Intermediate region name | e.g., Ji-Paraná, Porto Velho |
| `immediate_code` | int64 | 5,571 | 0.0% | Immediate region code | Range: [110001.00, 530001.00], Mean: 323826.89 |
| `immediate_name` | object | 5,571 | 0.0% | Immediate region name | e.g., Cacoal, Ariquemes |

### ibge_life_tables_combined

**Description:** IBGE life tables for YLL calculation

**File:** `ibge_life_tables_combined.parquet`

- **Rows:** 1,230
- **Columns:** 8
- **Spatial Unit:** national
- **Temporal Unit:** annual
- **Memory:** 0.1 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `age` | int64 | 1,230 | 0.0% | Age in years | Range: [0.00, 89.00], Mean: 40.60 |
| `qx` | float64 | 1,230 | 0.0% | Probability of dying between age x and x+1 | Range: [0.19, 131.17], Mean: 10.62 |
| `dx` | float64 | 1,230 | 0.0% |  | Range: [19.18, 3557.53], Mean: 690.46 |
| `lx` | float64 | 1,230 | 0.0% | Number surviving to age x | Range: [21799.28, 100000.00], Mean: 87822.93 |
| `Lx` | float64 | 1,230 | 0.0% |  | Range: [20369.57, 98968.27], Mean: 87470.92 |
| `Tx` | float64 | 1,230 | 0.0% |  | Range: [116034.10, 7702721.34], Mean: 3729914.07 |
| `ex` | float64 | 1,230 | 0.0% | Remaining life expectancy | Range: [5.29, 77.03], Mean: 39.80 |
| `year` | int64 | 1,230 | 0.0% | Reference year | Range: [2010.00, 2024.00], Mean: 2017.15 |

### yll_lookup_by_age

**Description:** YLL lookup table by single year of age

**File:** `yll_lookup_by_age.csv`

- **Rows:** 90
- **Columns:** 5
- **Spatial Unit:** national
- **Temporal Unit:** cross-sectional
- **Memory:** 0.0 MB

#### Columns

| Column | Type | Non-Null | Null% | Description | Statistics |
|--------|------|----------|-------|-------------|------------|
| `age` | int64 | 90 | 0.0% | Age in years (0-89) | Range: [0.00, 89.00], Mean: 44.50 |
| `ex_mean` | float64 | 90 | 0.0% | Mean remaining life expectancy across years | Range: [5.32, 75.71], Mean: 36.88 |
| `ex_min` | float64 | 90 | 0.0% | Minimum life expectancy | Range: [5.29, 74.11], Mean: 35.79 |
| `ex_max` | float64 | 90 | 0.0% | Maximum life expectancy | Range: [5.34, 77.03], Mean: 37.71 |
| `n_years` | int64 | 90 | 0.0% |  | Range: [3.00, 15.00], Mean: 13.67 |
