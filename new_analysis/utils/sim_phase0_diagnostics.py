"""
SIM–Phase 0 Sanity Diagnostics
==============================

This script performs non-intrusive checks comparing elderly deaths
in SIM DO*.csv files with the Phase 0 aggregated mortality file
`mortality_immediate_daily_elderly.parquet`.

It DOES NOT modify any inputs or existing analysis outputs.

Usage (from project root):

    python new_analysis/utils/sim_phase0_diagnostics.py

Outputs:
- Prints a small per-year table: SIM elderly vs Phase 0 elderly
- Saves a CSV summary under `new_analysis/phase1_core_model/results/`

This re-implements the SIM age parsing logic from
`01e_yll_unified.py` to avoid changing that script.
"""

import os
import re
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd


# -----------------------------------------------------------------------------
# PATHS
# -----------------------------------------------------------------------------

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
INPUT_DATA = os.path.join(BASE_DIR, "Input_data")
PHASE0_RESULTS = os.path.join(BASE_DIR, "new_analysis", "phase0_data_prep", "results")
PHASE1_RESULTS = os.path.join(BASE_DIR, "new_analysis", "phase1_core_model", "results")
os.makedirs(PHASE1_RESULTS, exist_ok=True)

MAPPING_FILE = os.path.join(PHASE0_RESULTS, "municipality_to_all_regions_map.csv")


# -----------------------------------------------------------------------------
# ROBUST SIM AGE PARSING (ALIGNED WITH PHASE 0 AGGREGATION)
# -----------------------------------------------------------------------------

def parse_sim_age(raw):
    """Parse SIM IDADE field to age in years (Phase 0–consistent).

    Mirrors `parse_sim_age` in
    `phase0_data_prep/aggregation/00n_aggregate_mortality_to_regions.py`:

    - If IDADE >= 400: age_years = IDADE - 400 (4YY years, 5YY 100+ years, etc.).
    - Else: age_years = 0 (under 1 year collapsed to 0).

    We add only minimal robustness for strings/decimals so that the
    aggregated elderly counts from DO*.csv can be compared one-to-one
    with Phase 0 outputs.
    """

    if pd.isna(raw):
        return None

    try:
        s = str(raw).strip()
        if "." in s:
            s = s.split(".")[0]
        digits = re.sub(r"\D", "", s)
        if digits == "":
            return None
        code = int(digits)
    except Exception:
        return None

    if code >= 400:
        return code - 400

    if code >= 0:
        return 0

    return None


def parse_sim_date(raw):
    """Parse SIM DTOBITO field to date (Phase 0–consistent).

    Mirrors the logic in 00n_aggregate_mortality_to_regions.py:
    - If string contains '-', interpret as ISO date (YYYY-MM-DD).
    - Otherwise treat as DDMMYYYY (numeric), tolerating floats/strings.
    Returns a pandas.Timestamp.date or None on failure.
    """
    try:
        s = str(raw).strip()
        if "-" in s:
            return pd.Timestamp(s).date()
        # Otherwise parse as DDMMYYYY
        s = str(int(float(s))).zfill(8)
        day = int(s[0:2])
        month = int(s[2:4])
        year = int(s[4:8])
        return pd.Timestamp(year=year, month=month, day=day).date()
    except Exception:
        return None


# -----------------------------------------------------------------------------
# MAIN DIAGNOSTIC LOGIC
# -----------------------------------------------------------------------------


def summarize_sim_elderly(min_elderly_age: int = 60) -> pd.DataFrame:
    """Summarise elderly deaths per year from SIM DO*.csv files.

    This function now reproduces the *full* Phase 0 filtering
    logic used in 00n_aggregate_mortality_to_regions.py:

    - Read DO*.csv with latin-1 encoding and inferred separator.
    - Parse DTOBITO to dates, drop failures.
    - Parse CODMUNRES to numeric muni_code, drop failures.
    - Filter sentinel municipality codes (>= 999990).
    - Map muni_code to immediate regions using the
      municipality_to_all_regions_map.csv (6- and 7-digit).
    - Parse IDADE to age in years with Phase 0–consistent rules.
    - Count elderly deaths as rows with age >= min_elderly_age and
      a non-null immediate_code.

    The resulting elderly counts should match the totals from
    mortality_immediate_daily_elderly.parquet exactly in years
    where mapping and inputs are consistent.
    """

    # Load mapping once
    if not os.path.exists(MAPPING_FILE):
        raise FileNotFoundError(MAPPING_FILE)

    df_mapping = pd.read_csv(MAPPING_FILE)
    df_mapping["code_muni"] = df_mapping["code_muni"].astype(str)
    df_mapping["code_muni_6"] = df_mapping["code_muni"].str[:6].astype(int)

    muni_to_immediate_6 = dict(zip(df_mapping["code_muni_6"], df_mapping["immediate_code"]))
    muni_to_immediate_7 = dict(zip(df_mapping["code_muni"], df_mapping["immediate_code"]))

    sim_dir = Path(INPUT_DATA)
    do_files = sorted(glob(str(sim_dir / "DO*.csv")))

    rows = []

    for file_path in do_files:
        file_path = Path(file_path)
        file_name = file_path.name

        # Extract year from filename: DO10OPEN.csv -> 2010
        try:
            year_code = file_name[2:4]
            year = 2000 + int(year_code)
        except Exception:
            year = None

        try:
            # Detect separator similarly to 00n
            with open(file_path, "r", encoding="latin1") as f:
                first_line = f.readline()
            if '","' in first_line or first_line.count(",") > first_line.count(";"):
                sep = ","
            else:
                sep = ";"

            df = pd.read_csv(
                file_path,
                sep=sep,
                encoding="latin1",
                dtype={"DTOBITO": str, "CODMUNRES": str, "IDADE": str, "CAUSABAS": str},
                usecols=["DTOBITO", "CODMUNRES", "IDADE", "CAUSABAS"],
                low_memory=False,
            )
        except Exception as e:
            print(f"[SIM] {file_name}: read error -> {e}")
            continue

        total_rows = len(df)

        # Parse date (Phase 0–style)
        df["date"] = df["DTOBITO"].apply(parse_sim_date)
        df = df.dropna(subset=["date"])

        # Municipality code numeric
        df["muni_code"] = pd.to_numeric(df["CODMUNRES"], errors="coerce")
        df = df.dropna(subset=["muni_code"])
        df["muni_code"] = df["muni_code"].astype(int)

        # Filter sentinel municipality codes
        df = df[df["muni_code"] < 999990]

        # Map to immediate regions (6-digit, then 7-digit fallback)
        df["muni_code_str"] = df["muni_code"].astype(str)
        df["muni_code_6"] = df["muni_code_str"].str[:6].astype(int)
        df["immediate_code"] = df["muni_code_6"].map(muni_to_immediate_6)

        missing_imm = df["immediate_code"].isna()
        if missing_imm.any():
            df.loc[missing_imm, "immediate_code"] = df.loc[missing_imm, "muni_code_str"].map(muni_to_immediate_7)

        # Parse age
        df["age_parsed"] = df["IDADE"].apply(parse_sim_age)

        parsed_rows = int(df["age_parsed"].notna().sum())
        parse_rate = 100.0 * parsed_rows / total_rows if total_rows > 0 else 0.0

        # Elderly deaths that actually enter Phase 0: need age>=min_elderly_age
        # and a valid immediate_code
        elderly_mask = (
            df["age_parsed"].notna()
            & (df["age_parsed"] >= min_elderly_age)
            & df["immediate_code"].notna()
        )
        sim_elderly = int(elderly_mask.sum())

        # Also record how many elderly are lost due to mapping/filters
        # for context, although not printed in the compact table yet.
        elderly_any_age = (
            df["age_parsed"].notna()
            & (df["age_parsed"] >= min_elderly_age)
        )
        total_elderly_raw = int(elderly_any_age.sum())

        rows.append(
            {
                "year": year,
                "do_file": file_name,
                "total_rows": total_rows,
                "parsed_rows": parsed_rows,
                "parse_rate": parse_rate,
                "sim_elderly_deaths": sim_elderly,
                "sim_elderly_raw": total_elderly_raw,
            }
        )

    if not rows:
        return pd.DataFrame(
            columns=[
                "year",
                "do_file",
                "total_rows",
                "parsed_rows",
                "parse_rate",
                "sim_elderly_deaths",
                "sim_elderly_raw",
            ]
        )

    return pd.DataFrame(rows)


def summarize_phase0_elderly() -> pd.DataFrame:
    """Summarise elderly deaths per year from Phase 0 parquet.

    Uses `mortality_immediate_daily_elderly.parquet`.

    Returns
    -------
    DataFrame
        Columns: [year, phase0_elderly_deaths]
    """
    mort_file = os.path.join(PHASE0_RESULTS, "mortality_immediate_daily_elderly.parquet")

    if not os.path.exists(mort_file):
        raise FileNotFoundError(mort_file)

    df_mort = pd.read_parquet(mort_file)
    if "date" not in df_mort.columns:
        raise ValueError("mortality file must contain a 'date' column")

    df_mort["date"] = pd.to_datetime(df_mort["date"])
    df_mort["year"] = df_mort["date"].dt.year

    grouped = (
        df_mort.groupby("year")["deaths_elderly"].sum().reset_index()
    )
    grouped = grouped.rename(columns={"deaths_elderly": "phase0_elderly_deaths"})

    return grouped


def main():
    print("=" * 70)
    print("SIM vs Phase 0 Elderly Mortality Diagnostics")
    print("=" * 70)

    # SIM summary
    sim_df = summarize_sim_elderly(min_elderly_age=60)

    if sim_df.empty:
        print("No SIM DO*.csv files could be summarised.")
        return

    # Phase 0 summary
    phase0_df = summarize_phase0_elderly()

    # Merge by year (outer join to see discrepancies clearly)
    merged = pd.merge(sim_df, phase0_df, on="year", how="outer")

    # Compute ratio where possible
    merged["ratio_sim_to_phase0"] = np.where(
        merged["phase0_elderly_deaths"].notna() & (merged["phase0_elderly_deaths"] > 0),
        merged["sim_elderly_deaths"] / merged["phase0_elderly_deaths"],
        np.nan,
    )

    # Sort for readability
    merged = merged.sort_values(["year", "do_file"], na_position="last")

    # Print compact table
    print("\nPer-year comparison (SIM vs Phase 0, with Phase 0 filters applied):")
    print("year  do_file    sim_elderly  phase0_elderly  parse_%   ratio_sim/phase0")
    for _, row in merged.iterrows():
        year = int(row["year"]) if pd.notna(row["year"]) else -1
        do_file = str(row["do_file"]) if pd.notna(row["do_file"]) else "-"
        sim_elderly = int(row["sim_elderly_deaths"]) if pd.notna(row["sim_elderly_deaths"]) else 0
        phase0 = int(row["phase0_elderly_deaths"]) if pd.notna(row["phase0_elderly_deaths"]) else 0
        parse_rate = float(row["parse_rate"]) if pd.notna(row["parse_rate"]) else float("nan")
        ratio = float(row["ratio_sim_to_phase0"]) if pd.notna(row["ratio_sim_to_phase0"]) else float("nan")

        print(
            f"{year:<5} {do_file:<10} "
            f"{sim_elderly:>10} {phase0:>15} "
            f"{parse_rate:7.1f} {ratio:16.3f}"
        )

    # Save CSV summary alongside other Phase 1 outputs
    out_file = os.path.join(PHASE1_RESULTS, "sim_vs_phase0_elderly_by_year.csv")
    merged.to_csv(out_file, index=False)
    print(f"\nSaved detailed summary to: {out_file}")


if __name__ == "__main__":
    main()
