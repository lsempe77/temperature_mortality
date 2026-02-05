# Data Validation Report

**Generated:** 2025-12-10 08:11:46

This report validates completeness, missingness, and alignment across all datasets.

---

## Validation Summary

- **Datasets defined:** 17
- **Datasets found:** 17
- **Datasets passed validation:** 17
- **Datasets with issues:** 0

## Dataset Status

| Dataset | Status | Issues | Warnings |
|---------|--------|--------|----------|
| era5_intermediate_daily | ✅ PASS | 0 | 0 |
| era5_immediate_daily | ✅ PASS | 0 | 0 |
| cams_intermediate_daily | ✅ PASS | 0 | 0 |
| cams_immediate_daily | ✅ PASS | 0 | 0 |
| mortality_intermediate_daily | ✅ PASS | 0 | 0 |
| mortality_intermediate_daily_elderly | ✅ PASS | 0 | 0 |
| mortality_immediate_daily | ✅ PASS | 0 | 0 |
| mortality_immediate_daily_elderly | ✅ PASS | 0 | 0 |
| influenza_daily_by_intermediate_region | ✅ PASS | 0 | 0 |
| influenza_daily_by_immediate_region | ✅ PASS | 0 | 0 |
| brazilian_holidays_daily | ✅ PASS | 0 | 0 |
| ses_intermediate_covariates | ✅ PASS | 0 | 0 |
| ses_immediate_covariates | ✅ PASS | 0 | 0 |
| regional_covariates | ✅ PASS | 0 | 0 |
| municipality_to_all_regions_map | ✅ PASS | 0 | 0 |
| ibge_life_tables_combined | ✅ PASS | 0 | 0 |
| yll_lookup_by_age | ✅ PASS | 0 | 0 |

## Detailed Issues

*No issues or warnings found.*

## Cross-Dataset Alignment

### Temporal Coverage

| Dataset | Start Date | End Date | Days |
|---------|------------|----------|------|
| era5_intermediate_daily | 2010-01-01 | 2024-12-31 | 5,478 |
| era5_immediate_daily | 2010-01-01 | 2024-12-31 | 5,478 |
| cams_intermediate_daily | 2010-01-01 | 2024-12-31 | 5,478 |
| cams_immediate_daily | 2010-01-01 | 2024-12-31 | 5,478 |
| mortality_intermediate_daily | 2010-01-01 | 2024-12-31 | 5,478 |
| mortality_intermediate_daily_elderly | 2010-01-01 | 2024-12-31 | 5,478 |
| mortality_immediate_daily | 2010-01-01 | 2024-12-31 | 5,478 |
| mortality_immediate_daily_elderly | 2010-01-01 | 2024-12-31 | 5,478 |
| influenza_daily_by_intermediate_region | 2010-01-03 | 2024-12-28 | 5,473 |
| influenza_daily_by_immediate_region | 2010-01-03 | 2024-12-28 | 5,473 |
| brazilian_holidays_daily | 2010-01-01 | 2024-12-31 | 5,478 |

### Spatial Coverage

| Dataset | Spatial Unit | Unique Regions | Expected | Status |
|---------|--------------|----------------|----------|--------|
| era5_intermediate_daily | intermediate | 133 | 133 | ✅ |
| era5_immediate_daily | immediate | 510 | 510 | ✅ |
| cams_intermediate_daily | intermediate | 133 | 133 | ✅ |
| cams_immediate_daily | immediate | 510 | 510 | ✅ |
| mortality_intermediate_daily | intermediate | 133 | 133 | ✅ |
| mortality_intermediate_daily_elderly | intermediate | 133 | 133 | ✅ |
| mortality_immediate_daily | immediate | 510 | 510 | ✅ |
| mortality_immediate_daily_elderly | immediate | 510 | 510 | ✅ |
| influenza_daily_by_intermediate_region | intermediate | 133 | 133 | ✅ |
| influenza_daily_by_immediate_region | immediate | 510 | 510 | ✅ |
| brazilian_holidays_daily | national | — | 1 | ⚠️ |
| ses_intermediate_covariates | intermediate | 133 | 133 | ✅ |
| ses_immediate_covariates | immediate | 510 | 510 | ✅ |
| regional_covariates | intermediate | 133 | 133 | ✅ |
| municipality_to_all_regions_map | municipality | 133 | 5570 | ⚠️ |
| ibge_life_tables_combined | national | — | 1 | ⚠️ |
| yll_lookup_by_age | national | — | 1 | ⚠️ |

### Missingness Summary

*Only showing columns with >5% missing values*

| Dataset | Column | Missing % | Assessment |
|---------|--------|-----------|------------|
| brazilian_holidays_daily | holiday_name_pt | 96.4% | ❌ High |
| brazilian_holidays_daily | holiday_type | 96.4% | ❌ High |
