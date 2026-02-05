# Archived Utilities Scripts

**Archived:** December 9, 2025

---

## Archived Scripts

### `00a_document_data_v1.py`
**Original purpose:** Basic data documentation

**Why archived:** 
- Replaced by comprehensive v2 script that generates 4 complete outputs:
  1. `data_dictionary.md` - Column-by-column documentation
  2. `schema_specification.md` - Unified schema for appendix
  3. `validation_report.md` - Completeness, missingness, alignment checks
  4. `metadata.json` - Structured metadata (rows, dates, regions, types)

### `00i_check_spatial_aggregation.py`
**Original purpose:** Create municipality-to-region mapping file

**Why archived:** 
- Only created `municipality_to_region_map.csv` with intermediate regions (133)
- Superseded by `00j_download_ibge_mapping.py` which creates `municipality_to_all_regions_map.csv` with BOTH:
  - `intermediate_code` / `intermediate_name` (133 regions)
  - `immediate_code` / `immediate_name` (510 regions)

**Replacement:** Use `phase0_data_prep/downloads/00j_download_ibge_mapping.py` instead

---

## Active Utilities Scripts

| Script | Purpose | Outputs |
|--------|---------|---------|
| `00a_document_data.py` | Comprehensive data documentation | `data_dictionary.md`, `schema_specification.md`, `validation_report.md`, `metadata.json` |
