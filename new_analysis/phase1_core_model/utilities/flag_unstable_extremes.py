"""Flag regions with extreme but statistically weak RRs.

This helper inspects the intermediate-level DLNM summary
(`dlnm_v2_intermediate_summary.csv`) and identifies regions
where the estimated RR at very cold (P1) or very hot (P99)
conditions is extremely large but *not* statistically
significant.

The goal is to quickly see which eye‑popping RRs are driven
by very imprecise estimates at the distribution tails.
"""

from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

# Two-sided alpha for p-values
ALPHA = 0.05

# Threshold above which we consider an RR "extreme"
EXTREME_RR_THRESHOLD = 10.0


def main() -> None:
    base_dir = Path(__file__).resolve().parents[1]
    results_dir = base_dir / "results"
    summary_path = results_dir / "dlnm_v2_intermediate_summary.csv"

    if not summary_path.exists():
        raise FileNotFoundError(f"Summary file not found: {summary_path}")

    df = pd.read_csv(summary_path)

    # Ensure required columns exist
    required_cols = [
        "region_code",
        "rr_p1",
        "rr_p1_p",
        "rr_p99",
        "rr_p99_p",
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in summary CSV: {missing}")

    cols_to_keep = [
        "region_code",
        "n_obs",
        "mean_deaths",
        "dispersion",
        "rr_p1",
        "rr_p1_lower",
        "rr_p1_upper",
        "rr_p1_p",
        "rr_p99",
        "rr_p99_lower",
        "rr_p99_upper",
        "rr_p99_p",
    ]
    cols_to_keep = [c for c in cols_to_keep if c in df.columns]

    # Extreme but statistically non-significant cold effects (P1 vs MMT)
    extreme_cold = df[
        (df["rr_p1"] >= EXTREME_RR_THRESHOLD)
        & (df["rr_p1_p"] >= ALPHA)
    ][cols_to_keep].copy()

    # Extreme but statistically non-significant heat effects (P99 vs MMT)
    extreme_heat = df[
        (df["rr_p99"] >= EXTREME_RR_THRESHOLD)
        & (df["rr_p99_p"] >= ALPHA)
    ][cols_to_keep].copy()

    out_path = results_dir / "dlnm_v2_intermediate_unstable_extremes.csv"

    # Tag rows by tail type and combine
    extreme_cold["tail"] = "cold_p1"
    extreme_heat["tail"] = "heat_p99"
    combined = pd.concat([extreme_cold, extreme_heat], ignore_index=True)

    if len(combined) == 0:
        print("No regions with RR >=" f" {EXTREME_RR_THRESHOLD} and p>={ALPHA} at P1/P99.")
        return

    combined = combined.sort_values(["tail", "region_code"]).reset_index(drop=True)
    combined.to_csv(out_path, index=False)

    print(f"Flagged {len(combined)} region-tail combinations.")
    print(f"Saved details to: {out_path}")
    print("\nExample rows:")
    print(combined.head(10))


if __name__ == "__main__":  # pragma: no cover
    main()
