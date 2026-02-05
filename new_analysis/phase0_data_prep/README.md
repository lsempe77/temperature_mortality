# Phase 0: Data Preparation

**Last Updated:** December 9, 2025

This folder contains all scripts for downloading, processing, and aggregating raw data into analysis-ready datasets for the Brazil Heat-Mortality DLNM study.

## Folder Structure

```
phase0_data_prep/
├── downloads/          # Data download scripts (run once per year)
├── aggregation/        # Spatial aggregation scripts (run after downloads)
├── covariates/         # Covariate construction scripts (run once)
├── utilities/          # Validation and documentation tools
├── archive/            # Superseded/debug scripts (kept for reference)
├── docs/               # Setup guides (CAMS credentials, etc.)
├── results/            # Output files (parquet/csv datasets)
├── temp_cams/          # Temporary CAMS download files
├── temp_life_tables/   # Temporary life table processing files
└── temp_srag/          # Temporary SRAG/influenza files
```

## Quick Reference

### Downloads (Run Once/Yearly)

| Script | Description |
|--------|-------------|
| `00b_download_cams_pollution_v2.py` | Download CAMS PM2.5/O3 from Copernicus ADS |
| `00h_download_era5_brazil.py` | Download ERA5 temperature from CDS API |
| `00j_download_ibge_mapping.py` | Download municipality-region mappings from IBGE |
| `00k2_download_ses_data_fixed.py` | Download SES indicators from SIDRA API |
| `00l_download_ibge_life_tables.py` | Process IBGE life tables for YLL calculations |
| `00m_download_influenza_data.py` | Download SRAG/influenza from OpenDataSUS |

### Aggregation (Run After Downloads)

| Script | Description |
|--------|-------------|
| `00g4_aggregate_cams_optimized.py` | Aggregate CAMS to 133/510 regions (optimized) |
| `00h2_aggregate_era5_to_regions.py` | Aggregate ERA5 to 133/510 regions |
| `00k3_aggregate_ses_to_regions.py` | Aggregate SES indicators to regions |
| `00m2_aggregate_influenza_municipal.py` | Aggregate influenza to municipalities |
| `00m3_process_additional_influenza.py` | Process 2019-2024 influenza files |
| `00m4_aggregate_influenza_immediate.py` | Aggregate influenza to 510 regions |
| `00n_aggregate_mortality_to_regions.py` | Aggregate SIM mortality to 133 regions |
| `00n2_aggregate_mortality_immediate.py` | Aggregate SIM mortality to 510 regions |

### Covariates (Run Once)

| Script | Description |
|--------|-------------|
| `00d_brazilian_holidays.py` | Create Brazilian holidays dataset 2010-2024 |
| `00e_pnad_ac_data.py` | Extract AC ownership by state from PNAD |
| `00f_state_covariates.py` | Compile state-level covariates |
| `00f2_regional_covariates.py` | Map state covariates to 133 intermediate regions |

### Utilities

| Script | Description |
|--------|-------------|
| `00a_document_data.py` | Generate data dictionary and summary |
| `00i_check_spatial_aggregation.py` | Validate region mapping files |

## Execution Order

1. **Download raw data** (downloads/)
2. **Build covariates** (covariates/) - can run in parallel with step 1
3. **Aggregate to regions** (aggregation/)
4. **Validate and document** (utilities/)

## Key Output Files

| File | Location | Description | Rows |
|------|----------|-------------|------|
| `era5_intermediate_daily.parquet` | results/ | Daily temperature (133 regions) | 728,707 |
| `era5_immediate_daily.parquet` | results/ | Daily temperature (510 regions) | 2,794,290 |
| `mortality_regional_daily.parquet` | results/ | Daily mortality (133 regions) | 721,096 |
| `mortality_immediate_daily.parquet` | results/ | Daily mortality (510 regions) | 2,523,235 |
| `ses_intermediate_covariates.csv` | results/ | SES by region (133) | 133 |
| `ses_immediate_covariates.csv` | results/ | SES by region (510) | 510 |

## Archive Folder

Contains superseded scripts kept for reference:
- `00g_*` → Older CAMS aggregation versions (replaced by `00g4_`)
- `00k_*` → Older SES download (replaced by `00k2_`)
- `debug_sidra*.py` → SIDRA API debugging scripts
- `datasus2021fetch.R` → One-time R script for 2021 data

See `ANALYSIS_ROADMAP.md` in the parent folder for complete documentation.
