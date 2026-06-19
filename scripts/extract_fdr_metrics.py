r"""Extract FDR calibration metrics for generalizability analysis.

This script extracts quantitative metrics from per-origin and combined
FDR validation results, outputting tables suitable for publication.

Tables generated match the format in analysis/fdr_generalizability_analysis.md:
- Table 1: Sample Characteristics by Origin
- Table 2: FDR Calibration at Reference Rate (50% Masked)
- Table 4: Best Masking Rate Per Origin
- Tables 5a-5e: Detailed per-origin calibration at all masking rates
- Table 6: Combined calibration at all masking rates
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Default paths (updated for new workflow)
DEFAULT_PER_ORIGIN_DIR = (
    "/home/j-daniel/Documents/winnow/novoboard/benchmark_data/per_origin_processed"
)
DEFAULT_COMBINED_DIR = (
    "/home/j-daniel/Documents/winnow/novoboard/benchmark_data/novoboard_inputs"
)
DEFAULT_OUTPUT_DIR = "/home/j-daniel/repos/NovoBoard/analysis"

# Decoy rates: fraction of peaks KEPT (sampling rate)
# E.g., 0.10 = 10% kept = 90% masked, 0.90 = 90% kept = 10% masked
DECOY_RATES = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]

# Reference sampling rate for cross-origin comparison (fraction KEPT)
# 0.50 = 50% kept = 50% masked
REFERENCE_RATE = 0.50

# Top origins to generate detailed tables for (by sample size)
TOP_ORIGINS = ["hepg2", "gluc", "tplantibodies", "helaqc"]


@dataclass
class FDRMetrics:
    """FDR calibration metrics for a single origin/rate combination."""

    origin: str
    decoy_rate: float
    n_samples: int
    mae: float | None  # Mean Absolute Error (estimated vs true FDR)
    bias: float | None  # Average (true - estimated), positive = underestimation
    ids_at_1pct: int  # Identifications at 1% estimated FDR
    ids_at_5pct: int  # Identifications at 5% estimated FDR
    ids_at_10pct: int  # Identifications at 10% estimated FDR
    true_fdr_at_1pct_est: float | None  # True FDR at point where estimated_fdr = 1%
    true_fdr_at_5pct_est: float | None  # True FDR at point where estimated_fdr = 5%
    true_fdr_at_10pct_est: float | None  # True FDR at point where estimated_fdr = 10%
    calibration_slope: float | None  # Linear regression slope (estimated vs true)
    calibration_r2: float | None  # R-squared of linear fit
    is_valid: bool  # Whether the FDR curve varies meaningfully


def calculate_metrics(df: pl.DataFrame, origin: str, decoy_rate: float) -> FDRMetrics:
    """Calculate FDR calibration metrics from a q-values DataFrame."""
    n_samples = len(df)

    # Ensure numeric columns are floats (CSV may read them as strings)
    numeric_cols = ["estimated_fdr", "estimated_q_value", "true_fdr_peptide"]
    for col in numeric_cols:
        if col in df.columns and df[col].dtype == pl.Utf8:
            df = df.with_columns(pl.col(col).cast(pl.Float64))

    # Identifications at various FDR thresholds
    ids_at_1pct = int((df["estimated_q_value"] <= 0.01).sum())
    ids_at_5pct = int((df["estimated_q_value"] <= 0.05).sum())
    ids_at_10pct = int((df["estimated_q_value"] <= 0.10).sum())

    # Filter to 0-10% estimated FDR range for metric calculation
    df_filtered = df.filter(pl.col("estimated_fdr") <= 0.10)

    # Check if FDR curve varies meaningfully in the 0-10% range
    est_fdr = df_filtered["estimated_fdr"].to_numpy()
    true_fdr = df_filtered["true_fdr_peptide"].to_numpy()

    # Filter out NaN/inf values
    valid_mask = np.isfinite(est_fdr) & np.isfinite(true_fdr)
    est_fdr_valid = est_fdr[valid_mask]
    true_fdr_valid = true_fdr[valid_mask]

    # Check if FDR varies meaningfully (spans at least 1%)
    est_range = (
        est_fdr_valid.max() - est_fdr_valid.min() if len(est_fdr_valid) > 0 else 0
    )
    is_valid = len(est_fdr_valid) > 10 and est_range >= 0.01

    # Helper function to interpolate true FDR at a given estimated FDR threshold
    def get_true_fdr_at_threshold_interpolated(threshold: float) -> float | None:
        """Get true FDR at threshold using linear interpolation."""
        sorted_df = df.sort("estimated_fdr")
        est_vals = sorted_df["estimated_fdr"].to_numpy()
        true_vals = sorted_df["true_fdr_peptide"].to_numpy()

        # Handle empty arrays
        if len(est_vals) == 0:
            return None

        # Check if threshold is beyond the range
        if threshold > est_vals.max():
            return None
        if threshold <= est_vals.min():
            return float(true_vals[0])

        # Find the two points that bracket the threshold
        above_mask = est_vals >= threshold
        if not above_mask.any():
            return None

        above_idx = np.argmax(above_mask)  # First True index

        if above_idx == 0:
            return float(true_vals[0])

        # Interpolate between (below_idx, above_idx)
        below_idx = above_idx - 1
        x0, x1 = est_vals[below_idx], est_vals[above_idx]
        y0, y1 = true_vals[below_idx], true_vals[above_idx]

        # Linear interpolation: y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
        if x1 == x0:
            return float(y0)

        interpolated = y0 + (y1 - y0) * (threshold - x0) / (x1 - x0)
        return float(interpolated)

    if not is_valid:
        # FDR curve is too flat or sparse - cannot compute meaningful metrics
        # But we can still get True FDR at thresholds if they exist
        return FDRMetrics(
            origin=origin,
            decoy_rate=decoy_rate,
            n_samples=n_samples,
            mae=None,
            bias=None,
            ids_at_1pct=ids_at_1pct,
            ids_at_5pct=ids_at_5pct,
            ids_at_10pct=ids_at_10pct,
            true_fdr_at_1pct_est=get_true_fdr_at_threshold_interpolated(0.01),
            true_fdr_at_5pct_est=get_true_fdr_at_threshold_interpolated(0.05),
            true_fdr_at_10pct_est=get_true_fdr_at_threshold_interpolated(0.10),
            calibration_slope=None,
            calibration_r2=None,
            is_valid=False,
        )

    # Mean Absolute Error (in 0-10% range)
    mae = float(np.abs(est_fdr_valid - true_fdr_valid).mean())

    # Bias: average (true - estimated), positive means underestimation
    bias = float((true_fdr_valid - est_fdr_valid).mean())

    # True FDR at thresholds using interpolation
    true_fdr_at_1pct = get_true_fdr_at_threshold_interpolated(0.01)
    true_fdr_at_5pct = get_true_fdr_at_threshold_interpolated(0.05)
    true_fdr_at_10pct = get_true_fdr_at_threshold_interpolated(0.10)

    # Linear regression for calibration slope
    slope = None
    r2 = None

    if np.std(est_fdr_valid) > 1e-6:
        try:
            # Simple linear regression
            slope, intercept = np.polyfit(est_fdr_valid, true_fdr_valid, 1)
            slope = float(slope)
            # R-squared
            predicted = slope * est_fdr_valid + intercept
            ss_res = np.sum((true_fdr_valid - predicted) ** 2)
            ss_tot = np.sum((true_fdr_valid - np.mean(true_fdr_valid)) ** 2)
            r2 = float(1 - (ss_res / ss_tot)) if ss_tot > 0 else None
        except (np.linalg.LinAlgError, ValueError):
            # Regression failed - keep slope and r2 as None
            pass

    return FDRMetrics(
        origin=origin,
        decoy_rate=decoy_rate,
        n_samples=n_samples,
        mae=mae,
        bias=bias,
        ids_at_1pct=ids_at_1pct,
        ids_at_5pct=ids_at_5pct,
        ids_at_10pct=ids_at_10pct,
        true_fdr_at_1pct_est=true_fdr_at_1pct,
        true_fdr_at_5pct_est=true_fdr_at_5pct,
        true_fdr_at_10pct_est=true_fdr_at_10pct,
        calibration_slope=slope,
        calibration_r2=r2,
        is_valid=True,
    )


def extract_per_origin_metrics(per_origin_dir: Path) -> list[FDRMetrics]:
    """Extract metrics for all origins and decoy rates."""
    all_metrics = []

    for origin_path in sorted(per_origin_dir.iterdir()):
        if not origin_path.is_dir():
            continue

        origin = origin_path.name
        logger.info(f"Processing origin: {origin}")

        for rate in DECOY_RATES:
            rate_str = f"{rate:.2f}"
            qval_files = list(origin_path.glob(f"*decoy_{rate_str}*_fdr_qvalues.csv"))

            if not qval_files:
                logger.warning(f"  No q-values file for rate {rate_str}")
                continue

            df = pl.read_csv(qval_files[0])
            metrics = calculate_metrics(df, origin, rate)
            all_metrics.append(metrics)

    return all_metrics


def extract_combined_metrics(combined_dir: Path) -> list[FDRMetrics]:
    """Extract metrics for combined (all origins) data."""
    all_metrics = []

    logger.info("Processing combined data")

    for rate in DECOY_RATES:
        rate_str = f"{rate:.2f}"
        qval_files = list(combined_dir.glob(f"*decoy_{rate_str}*_fdr_qvalues.csv"))

        if not qval_files:
            logger.warning(f"  No q-values file for rate {rate_str}")
            continue

        df = pl.read_csv(qval_files[0])
        metrics = calculate_metrics(df, "combined", rate)
        all_metrics.append(metrics)

    return all_metrics


def metrics_to_dataframe(metrics: list[FDRMetrics]) -> pl.DataFrame:
    """Convert list of metrics to a Polars DataFrame."""
    return pl.DataFrame(
        [
            {
                "origin": m.origin,
                "decoy_rate": m.decoy_rate,
                "n_samples": m.n_samples,
                "mae": m.mae,
                "bias": m.bias,
                "ids_at_1pct": m.ids_at_1pct,
                "ids_at_5pct": m.ids_at_5pct,
                "ids_at_10pct": m.ids_at_10pct,
                "true_fdr_at_1pct_est": m.true_fdr_at_1pct_est,
                "true_fdr_at_5pct_est": m.true_fdr_at_5pct_est,
                "true_fdr_at_10pct_est": m.true_fdr_at_10pct_est,
                "calibration_slope": m.calibration_slope,
                "calibration_r2": m.calibration_r2,
                "is_valid": m.is_valid,
            }
            for m in metrics
        ]
    )


def find_best_rate_per_origin(df: pl.DataFrame) -> pl.DataFrame:
    """Find the decoy rate with lowest MAE for each origin."""
    # Filter to valid rows only (where MAE is computable)
    valid_df = df.filter(pl.col("mae").is_not_null())

    return (
        valid_df.group_by("origin")
        .agg(
            [
                pl.col("mae").min().alias("best_mae"),
                pl.col("decoy_rate").sort_by("mae").first().alias("best_rate"),
                pl.col("n_samples").first().alias("n_samples"),
            ]
        )
        .sort("origin")
    )


def format_value(val: float | None, fmt: str = ".4f", na_str: str = "—") -> str:
    """Format a value, returning na_str if None."""
    if val is None:
        return na_str
    return f"{val:{fmt}}"


def format_pct(val: float | None, na_str: str = "—") -> str:
    """Format a value as percentage, returning na_str if None."""
    if val is None:
        return na_str
    return f"{val:.1%}"


def generate_detailed_origin_table(df: pl.DataFrame, origin: str) -> str:
    """Generate Table 5x format for a single origin with all masking rates.

    Format matches the markdown:
    | Masked | Kept | MAE | Bias | Slope | R² | True@1% | True@5% | True@10% |
    """
    origin_df = df.filter(pl.col("origin") == origin).sort(
        "decoy_rate", descending=True
    )

    if len(origin_df) == 0:
        return f"No data for origin: {origin}"

    n_samples = origin_df["n_samples"][0]

    lines = [
        f"### {origin} (N={n_samples:,})",
        "",
        "| Masked | Kept | MAE | Bias | Slope | R² | True@1% | True@5% | True@10% |",
        "|--------|------|-----|------|-------|-----|---------|---------|----------|",
    ]

    # Group consecutive invalid rows for cleaner output
    rows = list(origin_df.iter_rows(named=True))
    i = 0
    while i < len(rows):
        row = rows[i]

        if not row["is_valid"]:
            # Find consecutive invalid rows
            j = i
            while j < len(rows) and not rows[j]["is_valid"]:
                j += 1

            # Format as range if multiple consecutive invalid rows
            if j - i > 1:
                start_rate = int(rows[i]["decoy_rate"] * 100)
                end_rate = int(rows[j - 1]["decoy_rate"] * 100)
                # decoy_rate is SAMPLING rate (kept), so swap the labels
                kept_str = f"{end_rate}-{start_rate}%"
                masked_str = f"{100-start_rate}-{100-end_rate}%"
            else:
                # decoy_rate is SAMPLING rate (kept), not masking rate
                kept_str = f"{int(row['decoy_rate'] * 100)}%"
                masked_str = f"{int((1 - row['decoy_rate']) * 100)}%"

            # Check if any True@X% values exist in this range
            true_1pct_vals = [r["true_fdr_at_1pct_est"] for r in rows[i:j]]
            true_5pct_vals = [r["true_fdr_at_5pct_est"] for r in rows[i:j]]
            true_10pct_vals = [r["true_fdr_at_10pct_est"] for r in rows[i:j]]

            # Use the first non-None value if any exist
            true_1pct = next((v for v in true_1pct_vals if v is not None), None)
            true_5pct = next((v for v in true_5pct_vals if v is not None), None)
            true_10pct = next((v for v in true_10pct_vals if v is not None), None)

            lines.append(
                f"| {masked_str} | {kept_str} | — | — | — | — | "
                f"{format_pct(true_1pct)} | {format_pct(true_5pct)} | {format_pct(true_10pct)} |"
            )
            i = j
        else:
            # decoy_rate is the SAMPLING rate (fraction KEPT), not masking rate
            # decoy_0.20.mgf means 20% kept, 80% removed
            kept_pct = int(row["decoy_rate"] * 100)
            masked_pct = int((1 - row["decoy_rate"]) * 100)

            # Check if this is the best row (lowest MAE) for highlighting
            is_best = False
            valid_rows = [r for r in rows if r["is_valid"] and r["mae"] is not None]
            if valid_rows:
                best_mae = min(r["mae"] for r in valid_rows)
                is_best = row["mae"] == best_mae

            # Format values
            mae_str = format_value(row["mae"], ".3f")
            bias_str = (
                f"+{row['bias']:.3f}"
                if row["bias"] and row["bias"] > 0
                else format_value(row["bias"], ".3f")
            )
            slope_str = format_value(row["calibration_slope"], ".2f")
            r2_str = format_value(row["calibration_r2"], ".2f")

            # True FDR values - use N/A if threshold not reached but curve is valid
            true_1pct_str = format_pct(row["true_fdr_at_1pct_est"], na_str="N/A")
            true_5pct_str = format_pct(row["true_fdr_at_5pct_est"], na_str="N/A")
            true_10pct_str = format_pct(row["true_fdr_at_10pct_est"], na_str="N/A")

            if is_best:
                lines.append(
                    f"| **{masked_pct}%** | **{kept_pct}%** | **{mae_str}** | **{bias_str}** | "
                    f"**{slope_str}** | **{r2_str}** | **{true_1pct_str}** | **{true_5pct_str}** | **{true_10pct_str}** |"
                )
            else:
                lines.append(
                    f"| {masked_pct}% | {kept_pct}% | {mae_str} | {bias_str} | "
                    f"{slope_str} | {r2_str} | {true_1pct_str} | {true_5pct_str} | {true_10pct_str} |"
                )
            i += 1

    return "\n".join(lines)


def generate_combined_table(df: pl.DataFrame) -> str:
    """Generate Table 6 format for combined data."""
    # Handle empty DataFrame (no combined data available)
    if len(df) == 0 or "origin" not in df.columns:
        return "No combined data available"

    combined_df = df.filter(pl.col("origin") == "combined").sort(
        "decoy_rate", descending=True
    )

    if len(combined_df) == 0:
        return "No combined data available"

    n_samples = combined_df["n_samples"][0]

    lines = [
        f"### Table 6: COMBINED (N={n_samples:,})",
        "",
        "| Masked | Kept | MAE | Bias | Slope | R² | True@1% | True@5% | True@10% |",
        "|--------|------|-----|------|-------|-----|---------|---------|----------|",
    ]

    rows = list(combined_df.iter_rows(named=True))
    best_mae = min((r["mae"] for r in rows if r["mae"] is not None), default=None)

    for row in rows:
        # decoy_rate is the SAMPLING rate (fraction KEPT), not masking rate
        kept_pct = int(row["decoy_rate"] * 100)
        masked_pct = int((1 - row["decoy_rate"]) * 100)

        is_best = row["mae"] == best_mae and best_mae is not None

        if not row["is_valid"]:
            true_1pct_str = format_pct(row["true_fdr_at_1pct_est"])
            true_5pct_str = format_pct(row["true_fdr_at_5pct_est"])
            true_10pct_str = format_pct(row["true_fdr_at_10pct_est"])
            lines.append(
                f"| {masked_pct}% | {kept_pct}% | — | — | — | — | {true_1pct_str} | {true_5pct_str} | {true_10pct_str} |"
            )
        else:
            mae_str = format_value(row["mae"], ".3f")
            bias_str = (
                f"+{row['bias']:.3f}"
                if row["bias"] and row["bias"] > 0
                else format_value(row["bias"], ".3f")
            )
            slope_str = format_value(row["calibration_slope"], ".2f")
            r2_str = format_value(row["calibration_r2"], ".2f")
            true_1pct_str = format_pct(row["true_fdr_at_1pct_est"], na_str="N/A")
            true_5pct_str = format_pct(row["true_fdr_at_5pct_est"], na_str="N/A")
            true_10pct_str = format_pct(row["true_fdr_at_10pct_est"], na_str="N/A")

            if is_best:
                lines.append(
                    f"| **{masked_pct}%** | **{kept_pct}%** | **{mae_str}** | **{bias_str}** | "
                    f"**{slope_str}** | **{r2_str}** | **{true_1pct_str}** | **{true_5pct_str}** | **{true_10pct_str}** |"
                )
            else:
                lines.append(
                    f"| {masked_pct}% | {kept_pct}% | {mae_str} | {bias_str} | "
                    f"{slope_str} | {r2_str} | {true_1pct_str} | {true_5pct_str} | {true_10pct_str} |"
                )

    return "\n".join(lines)


def generate_markdown_tables(
    per_origin_df: pl.DataFrame,
    combined_df: pl.DataFrame,
    output_path: Path,
) -> dict:
    """Generate markdown tables for the analysis document."""

    # Use reference rate for sample characteristics (any rate has same n_samples)
    ref_rate = REFERENCE_RATE
    ref_df = per_origin_df.filter(pl.col("decoy_rate") == ref_rate)

    # If reference rate not available, use the first available rate
    if len(ref_df) == 0:
        available_rates = per_origin_df["decoy_rate"].unique().to_list()
        if available_rates:
            ref_rate = available_rates[0]
            ref_df = per_origin_df.filter(pl.col("decoy_rate") == ref_rate)

    total_samples = ref_df["n_samples"].sum() if len(ref_df) > 0 else 0

    # Table 1: Sample Characteristics
    sample_df = (
        ref_df.select(["origin", "n_samples"])
        .sort("n_samples", descending=True)
        .with_columns(
            [(pl.col("n_samples") / total_samples * 100).round(1).alias("training_pct")]
        )
    )

    table1_lines = [
        "| Origin | Number of Spectra | Dataset Fraction |",
        "|--------|-------------------|------------------|",
    ]
    for row in sample_df.iter_rows(named=True):
        table1_lines.append(
            f"| `{row['origin']}` | {row['n_samples']:,} | {row['training_pct']:.1f}% |"
        )
    table1 = "\n".join(table1_lines)

    # Table 2: FDR Calibration at Reference Rate
    ref_calib_df = per_origin_df.filter(pl.col("decoy_rate") == REFERENCE_RATE).sort(
        "mae"
    )

    def calibration_quality(mae: float | None, true_fdr_1pct: float | None) -> str:
        if mae is None or true_fdr_1pct is None:
            return "—"
        if mae < 0.03 and true_fdr_1pct < 0.03:
            return "Good"
        elif mae < 0.05:
            return "Moderate"
        elif true_fdr_1pct > 0.10:
            return "Poor, severely underestimates FDR"
        else:
            return "Poor"

    # REFERENCE_RATE is the SAMPLING rate (fraction KEPT), so masked = 1 - REFERENCE_RATE
    ref_masked_pct = int((1 - REFERENCE_RATE) * 100)
    table2_lines = [
        f"### Table 2: FDR Calibration at {ref_masked_pct}% Masked (Reference Rate)",
        "",
        "| Origin | MAE | Slope | True@Est1% | True@Est5% | True@Est10% | Quality |",
        "|--------|-----|-------|------------|------------|-------------|---------|",
    ]
    for row in ref_calib_df.iter_rows(named=True):
        mae_str = format_value(row["mae"], ".4f")
        slope_str = format_value(row["calibration_slope"], ".2f")
        true_1pct_str = format_pct(row["true_fdr_at_1pct_est"], na_str="N/A")
        true_5pct_str = format_pct(row["true_fdr_at_5pct_est"], na_str="N/A")
        true_10pct_str = format_pct(row["true_fdr_at_10pct_est"], na_str="N/A")
        quality = calibration_quality(row["mae"], row["true_fdr_at_1pct_est"])
        table2_lines.append(
            f"| `{row['origin']}` | {mae_str} | {slope_str} | {true_1pct_str} | {true_5pct_str} | {true_10pct_str} | {quality} |"
        )
    table2 = "\n".join(table2_lines)

    # Table 4: Optimal Decoy Rates by Origin (includes Peaks Kept)
    best_rates = find_best_rate_per_origin(per_origin_df)

    table4_lines = [
        "### Table 4: Best Masking Rate Per Origin",
        "",
        "| Origin | Number of Spectra | Best Masking | Peaks Kept | Best MAE | Interpretation |",
        "|--------|-------------------|--------------|------------|----------|----------------|",
    ]
    for row in best_rates.sort("best_mae").iter_rows(named=True):
        if row["best_mae"] < 0.02:
            interp = "Near-perfect calibration"
        elif row["best_mae"] < 0.03:
            interp = "Well-calibrated"
        elif row["best_mae"] < 0.05:
            interp = "Acceptable"
        elif row["best_mae"] < 0.10:
            interp = "Moderate miscalibration"
        else:
            interp = "Severe miscalibration"

        # best_rate is the SAMPLING rate (fraction KEPT), not masking rate
        kept_pct = int(row["best_rate"] * 100)
        masked_pct = 100 - kept_pct
        table4_lines.append(
            f"| `{row['origin']}` | {row['n_samples']:,} | {masked_pct}% | {kept_pct}% | {row['best_mae']:.4f} | {interp} |"
        )
    table4 = "\n".join(table4_lines)

    # Tables 5a-5e: Detailed per-origin tables
    detailed_tables = {}
    origins = per_origin_df["origin"].unique().to_list()

    # Generate for top origins first, then others
    for origin in TOP_ORIGINS:
        if origin in origins:
            detailed_tables[f"table5_{origin}"] = generate_detailed_origin_table(
                per_origin_df, origin
            )

    for origin in sorted(origins):
        if origin not in TOP_ORIGINS:
            detailed_tables[f"table5_{origin}"] = generate_detailed_origin_table(
                per_origin_df, origin
            )

    # Table 6: Combined calibration
    table6 = generate_combined_table(combined_df)

    # Summary table: Best rates comparison
    summary_lines = [
        "### Summary: Best Masking Rates vs Perfect Diagonal",
        "",
        "| Dataset | N | Best Masking | MAE | Bias | Slope | R² | True@1% | True@5% | True@10% |",
        "|---------|---|--------------|-----|------|-------|-----|---------|---------|----------|",
    ]

    # Add best rows for each origin
    for origin in TOP_ORIGINS + ["combined"]:
        if origin == "combined":
            # Skip if combined_df is empty
            if len(combined_df) == 0 or "origin" not in combined_df.columns:
                continue
            origin_df = combined_df.filter(pl.col("origin") == "combined")
        else:
            origin_df = per_origin_df.filter(pl.col("origin") == origin)

        if len(origin_df) == 0:
            continue

        # Find best row
        valid_rows = origin_df.filter(pl.col("mae").is_not_null())
        if len(valid_rows) == 0:
            continue

        best_row = valid_rows.sort("mae").row(0, named=True)
        n = best_row["n_samples"]
        # decoy_rate is SAMPLING rate (kept), so masking = 100 - kept
        kept_pct = int(best_row["decoy_rate"] * 100)
        masked_pct = 100 - kept_pct

        mae_str = format_value(best_row["mae"], ".3f")
        bias_str = (
            f"+{best_row['bias']:.3f}"
            if best_row["bias"] and best_row["bias"] > 0
            else format_value(best_row["bias"], ".3f")
        )
        slope_str = format_value(best_row["calibration_slope"], ".2f")
        r2_str = format_value(best_row["calibration_r2"], ".2f")
        true_1pct_str = format_pct(best_row["true_fdr_at_1pct_est"], na_str="N/A")
        true_5pct_str = format_pct(best_row["true_fdr_at_5pct_est"], na_str="N/A")
        true_10pct_str = format_pct(best_row["true_fdr_at_10pct_est"], na_str="N/A")

        if origin == "combined":
            summary_lines.append(
                f"| **Combined** | {n:,} | **{masked_pct}%** | {mae_str} | {bias_str} | {slope_str} | {r2_str} | {true_1pct_str} | {true_5pct_str} | {true_10pct_str} |"
            )
        else:
            summary_lines.append(
                f"| `{origin}` | {n:,} | {masked_pct}% | {mae_str} | {bias_str} | {slope_str} | {r2_str} | {true_1pct_str} | {true_5pct_str} | {true_10pct_str} |"
            )

    summary_table = "\n".join(summary_lines)

    result = {
        "table1_sample_characteristics": table1,
        "table2_reference_rate_calibration": table2,
        "table4_best_rates": table4,
        "table6_combined": table6,
        "summary_best_rates": summary_table,
        "total_samples": total_samples,
        "n_origins": len(sample_df),
        "reference_rate": REFERENCE_RATE,
    }
    result.update(detailed_tables)

    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract FDR calibration metrics for generalizability analysis"
    )
    parser.add_argument(
        "--per-origin-dir",
        type=Path,
        default=Path(DEFAULT_PER_ORIGIN_DIR),
        help=f"Directory containing per-origin results (default: {DEFAULT_PER_ORIGIN_DIR})",
    )
    parser.add_argument(
        "--combined-dir",
        type=Path,
        default=Path(DEFAULT_COMBINED_DIR),
        help=f"Directory containing combined results (default: {DEFAULT_COMBINED_DIR})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(DEFAULT_OUTPUT_DIR),
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--print-tables",
        action="store_true",
        help="Print generated markdown tables to stdout",
    )

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Extract metrics
    logger.info("Extracting per-origin metrics...")
    per_origin_metrics = extract_per_origin_metrics(args.per_origin_dir)
    per_origin_df = metrics_to_dataframe(per_origin_metrics)

    logger.info("Extracting combined metrics...")
    combined_metrics = extract_combined_metrics(args.combined_dir)
    combined_df = metrics_to_dataframe(combined_metrics)

    # Save raw data
    per_origin_df.write_csv(args.output_dir / "per_origin_metrics.csv")
    combined_df.write_csv(args.output_dir / "combined_metrics.csv")
    logger.info(f"Saved raw metrics to {args.output_dir}")

    # Generate markdown tables
    tables = generate_markdown_tables(per_origin_df, combined_df, args.output_dir)

    # Save tables as JSON for use in document generation
    with open(args.output_dir / "tables.json", "w") as f:
        json.dump(tables, f, indent=2)

    # Save tables as markdown file
    md_output = args.output_dir / "generated_tables.md"
    # REFERENCE_RATE is the SAMPLING rate (fraction KEPT), so masked = 1 - REFERENCE_RATE
    ref_masked_pct = int((1 - REFERENCE_RATE) * 100)
    with open(md_output, "w") as f:
        f.write("# Auto-Generated FDR Metrics Tables\n\n")
        f.write(f"Reference masking rate: {ref_masked_pct}%\n\n")
        f.write("---\n\n")

        f.write("## Table 1: Sample Characteristics\n\n")
        f.write(tables["table1_sample_characteristics"])
        f.write("\n\n---\n\n")

        f.write(tables["table2_reference_rate_calibration"])
        f.write("\n\n---\n\n")

        f.write(tables["table4_best_rates"])
        f.write("\n\n---\n\n")

        f.write(tables["summary_best_rates"])
        f.write("\n\n---\n\n")

        f.write("## Detailed Per-Origin Tables\n\n")
        for key, value in tables.items():
            if key.startswith("table5_"):
                f.write(value)
                f.write("\n\n")

        f.write("---\n\n")
        f.write(tables["table6_combined"])
        f.write("\n")

    logger.info(f"Saved markdown tables to {md_output}")

    # Print summary
    print("\n" + "=" * 60)
    print("FDR METRICS EXTRACTION COMPLETE")
    print("=" * 60)
    print(f"\nTotal samples analyzed: {tables['total_samples']:,}")
    print(f"Number of origins: {tables['n_origins']}")
    # REFERENCE_RATE is the SAMPLING rate (fraction KEPT), so masked = 1 - REFERENCE_RATE
    print(f"Reference masking rate: {int((1 - REFERENCE_RATE) * 100)}%")
    print("\nOutput files:")
    print(f"  - {args.output_dir / 'per_origin_metrics.csv'}")
    print(f"  - {args.output_dir / 'combined_metrics.csv'}")
    print(f"  - {args.output_dir / 'tables.json'}")
    print(f"  - {md_output}")

    if args.print_tables:
        print("\n" + "=" * 60)
        print("GENERATED TABLES")
        print("=" * 60)
        with open(md_output) as f:
            print(f.read())

    print("\n" + "=" * 60)


if __name__ == "__main__":
    main()
