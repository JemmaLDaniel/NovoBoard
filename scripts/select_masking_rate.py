#!/usr/bin/env python3
"""Select the masking rate whose FDR calibration curve is closest to ideal (y=x).

For each origin and decoy rate, reads the per-PSM FDR q-values CSV produced by
`novoboard fdr`, computes the MSE between estimated FDR and true FDR (i.e. the
squared deviation from the ideal diagonal), and reports which masking rate
minimises that error.

Usage:
    python scripts/select_masking_rate.py 1000
    python scripts/select_masking_rate.py 2000
    python scripts/select_masking_rate.py 1000 --fdr-max 0.1
    python scripts/select_masking_rate.py 2000 --fdr-max 1.0   # full range
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


TUNING_DIR = Path("paper_analysis/tuning")
ORIGINS = ["celegans", "sbrodae", "PXD019483"]
RATES = ["0.30", "0.40", "0.50", "0.60", "0.70"]


def qvalues_path(origin: str, sample_n: int, rate: str) -> Path:
    prefix = f"annotated_train_sample{sample_n}"
    fname = f"{prefix}_denovo_fdr_{prefix}_decoy_{rate}_denovo_fdr_qvalues.csv"
    return TUNING_DIR / origin / fname


def compute_mse(
    df: pd.DataFrame,
    tp_metric: str,
    fdr_max: float | None,
) -> float | None:
    """MSE between true FDR and estimated FDR (deviation from y=x diagonal)."""
    true_col = f"true_fdr_{tp_metric}"
    if true_col not in df.columns:
        return None

    mask = pd.notna(df["estimated_fdr"]) & pd.notna(df[true_col])
    if fdr_max is not None:
        mask &= df["estimated_fdr"] <= fdr_max
    est = df.loc[mask, "estimated_fdr"].values
    true = df.loc[mask, true_col].values

    if len(est) == 0:
        return None
    return float(np.mean((est - true) ** 2))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Select masking rate closest to ideal FDR calibration."
    )
    parser.add_argument(
        "sample_n",
        type=int,
        help="Sample size (1000 or 2000)",
    )
    parser.add_argument(
        "--fdr-max",
        type=float,
        default=None,
        help="Only consider PSMs with estimated FDR <= this value. "
        "Default: use full range.",
    )
    parser.add_argument(
        "--tp-metric",
        type=str,
        choices=["peptide", "ion100", "ion90"],
        default="peptide",
        help="True positive metric column suffix (default: peptide)",
    )
    args = parser.parse_args()

    fdr_label = f"0–{args.fdr_max}" if args.fdr_max else "full range"
    print(f"Sample size: {args.sample_n}")
    print(f"FDR range:   {fdr_label}")
    print(f"TP metric:   {args.tp_metric}")
    print()

    for origin in ORIGINS:
        print(f"=== {origin} ===")
        results: list[tuple[str, float]] = []

        for rate in RATES:
            path = qvalues_path(origin, args.sample_n, rate)
            if not path.exists():
                print(f"  rate {rate}: MISSING ({path.name})")
                continue

            df = pd.read_csv(path)
            mse = compute_mse(df, args.tp_metric, args.fdr_max)
            if mse is None:
                print(f"  rate {rate}: no valid data")
                continue

            kept_pct = int(float(rate) * 100)
            masked_pct = 100 - kept_pct
            results.append((rate, mse))
            print(f"  {masked_pct}% masked (kept={rate}): MSE = {mse:.6e}")

        if results:
            best_rate, best_mse = min(results, key=lambda x: x[1])
            best_kept_pct = int(float(best_rate) * 100)
            best_masked_pct = 100 - best_kept_pct
            print(
                f"  --> Best: {best_masked_pct}% masked "
                f"(kept rate = {best_rate}, MSE = {best_mse:.6e})"
            )
        else:
            print("  --> No data found")
        print()


if __name__ == "__main__":
    main()
