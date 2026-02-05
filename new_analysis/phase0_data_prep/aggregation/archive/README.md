# Archived Aggregation Scripts

**Archived:** December 9, 2025  
**Reason:** Scripts made redundant by consolidating immediate + intermediate outputs into main scripts

---

## Archived Scripts

### `00m4_aggregate_influenza_immediate.py`
**Original purpose:** Aggregate influenza data to 510 immediate regions only

**Why archived:** The main script `00m2_aggregate_influenza_municipal.py` was updated to:
- Use `municipality_to_all_regions_map.csv` (has both region levels)
- Output BOTH intermediate (133) and immediate (510) region files
- New outputs:
  - `influenza_daily_by_intermediate_region.parquet`
  - `influenza_daily_by_immediate_region.parquet`
  - `influenza_weekly_by_intermediate_region.parquet`
  - `influenza_weekly_by_immediate_region.parquet`

### `00n2_aggregate_mortality_immediate.py`
**Original purpose:** Aggregate mortality data to 510 immediate regions only

**Why archived:** The main script `00n_aggregate_mortality_to_regions.py` was updated to:
- Use `municipality_to_all_regions_map.csv` (has both region levels)
- Output BOTH intermediate (133) and immediate (510) region files
- New outputs:
  - `mortality_intermediate_daily.parquet`
  - `mortality_intermediate_daily_elderly.parquet`
  - `mortality_immediate_daily.parquet`
  - `mortality_immediate_daily_elderly.parquet`

---

## Key Changes Made to Main Scripts

### Column Name Updates
The old mapping file (`municipality_to_region_map.csv`) used:
- `region_code`, `region_name`

The new mapping file (`municipality_to_all_regions_map.csv`) uses:
- `intermediate_code`, `intermediate_name`
- `immediate_code`, `immediate_name`

All scripts were updated to use the new column names.

---

## Current Active Aggregation Scripts

| Script | Intermediate (133) | Immediate (510) |
|--------|:------------------:|:---------------:|
| `00g4_aggregate_cams_optimized.py` | ✅ | ✅ |
| `00h2_aggregate_era5_to_regions.py` | ✅ | ✅ |
| `00k3_aggregate_ses_to_regions.py` | ✅ | ✅ |
| `00m2_aggregate_influenza_municipal.py` | ✅ | ✅ |
| `00m3_process_additional_influenza.py` | ✅ | ✅ |
| `00n_aggregate_mortality_to_regions.py` | ✅ | ✅ |

All active scripts now output both region levels from a single run.
