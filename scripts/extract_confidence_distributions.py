r"""Extract confidence score distributions for TPs, FPs, and decoys per origin.

Computes mean/median confidence scores for:
- True Positives (TPs): correct predictions on target spectra
- False Positives (FPs): incorrect predictions on target spectra
- Decoys: predictions on decoy (corrupted) spectra

This analysis helps explain FDR calibration differences:
- Good calibration: FP and decoy confidence distributions are similar
- Underestimation: FPs have higher confidence than decoys (model is fooled by real spectra)
- Overestimation: Decoys have higher confidence than FPs

Usage:
    python scripts/extract_confidence_distributions.py
    python scripts/extract_confidence_distributions.py --retention-rate 0.50
"""

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_PROCESSED_DIR = Path(
    "/home/j-daniel/Documents/winnow/novoboard/benchmark_data/per_origin_processed"
)


def extract_confidence_stats(
    origin: str, retention_rate: float, processed_dir: Path
) -> dict | None:
    """Extract confidence statistics for a single origin at a given retention rate."""
    origin_dir = processed_dir / origin
    if not origin_dir.exists():
        return None

    rate_str = f"{retention_rate:.2f}"

    # Find the combined target-decoy FDR file
    # Pattern: {origin}_de_novo_results_fdr_{origin}_decoy_{rate}_de_novo_results.csv
    fdr_files = list(
        origin_dir.glob(
            f"{origin}_de_novo_results_fdr_{origin}_decoy_{rate_str}_de_novo_results.csv"
        )
    )

    if not fdr_files:
        # Try alternative pattern
        fdr_files = list(
            origin_dir.glob(f"*_fdr_{origin}_decoy_{rate_str}*_de_novo_results.csv")
        )
        # Filter out qvalues, accuracy, and other derived files
        fdr_files = [
            f
            for f in fdr_files
            if not any(
                x in f.name
                for x in ["qvalues", "accuracy", "multifea", "scan2fea", "denovo_only"]
            )
        ]

    if not fdr_files:
        return None

    fdr_file = fdr_files[0]

    # Find the accuracy file (contains correctness info)
    accuracy_files = list(
        origin_dir.glob(f"*_fdr_{origin}_decoy_{rate_str}*_accuracy.csv")
    )
    if not accuracy_files:
        return None

    accuracy_file = accuracy_files[0]

    # Load data
    try:
        fdr_df = pd.read_csv(fdr_file)
        accuracy_df = pd.read_csv(accuracy_file, delimiter="\t")
    except Exception as e:
        print(f"  Error loading {origin}: {e}")
        return None

    # Get score column (try common names)
    score_col = None
    for col in ["ALC (%)", "score", "Score", "confidence"]:
        if col in fdr_df.columns:
            score_col = col
            break

    if score_col is None:
        print(f"  No score column found in {fdr_file}")
        return None

    # Separate targets and decoys
    if "is_target" not in fdr_df.columns:
        print(f"  No is_target column in {fdr_file}")
        return None

    targets = fdr_df[fdr_df["is_target"]].copy()
    decoys = fdr_df[~fdr_df["is_target"]].copy()

    # Merge accuracy info with targets
    if "feature_id" in targets.columns and "feature_id" in accuracy_df.columns:
        targets = targets.set_index("feature_id")
        accuracy_df = accuracy_df.set_index("feature_id")

        # Add correctness info
        targets["matched_aa"] = accuracy_df["matched_amino_acid_count"]
        targets["predicted_len"] = accuracy_df["predicted_sequence_length"]
        targets["is_correct"] = targets["matched_aa"] == targets["predicted_len"]

        # Separate TPs and FPs
        tps = targets[targets["is_correct"]]
        fps = targets[~targets["is_correct"]]
    else:
        print("  No feature_id column for merging")
        return None

    # Compute statistics
    n_targets = len(targets)
    n_tps = len(tps)
    n_fps = len(fps)
    n_decoys = len(decoys)

    if n_tps == 0 or n_fps == 0 or n_decoys == 0:
        return None

    # Mean confidence
    tp_mean = tps[score_col].mean()
    fp_mean = fps[score_col].mean()
    decoy_mean = decoys[score_col].mean()

    # Median confidence
    tp_median = tps[score_col].median()
    fp_median = fps[score_col].median()
    decoy_median = decoys[score_col].median()

    # Accuracy
    accuracy = n_tps / n_targets * 100

    return {
        "origin": origin,
        "masking_rate": 1 - retention_rate,
        "n_targets": n_targets,
        "n_tps": n_tps,
        "n_fps": n_fps,
        "n_decoys": n_decoys,
        "accuracy_pct": accuracy,
        "tp_mean_conf": tp_mean,
        "fp_mean_conf": fp_mean,
        "decoy_mean_conf": decoy_mean,
        "tp_median_conf": tp_median,
        "fp_median_conf": fp_median,
        "decoy_median_conf": decoy_median,
        "fp_decoy_gap": fp_mean
        - decoy_mean,  # Positive = FPs more confident than decoys
    }


def generate_markdown_table(results: list[dict], sampling_rate: float) -> str:
    """Generate markdown table from results.

    Args:
        results: List of confidence statistics dictionaries
        sampling_rate: Fraction of peaks KEPT (0.1 = 10% kept = 90% masked)
    """
    # sampling_rate is fraction KEPT, so masked = 1 - sampling_rate
    masked_pct = int((1 - sampling_rate) * 100)
    lines = [
        f"## Confidence Score Distributions ({masked_pct}% Masked)",
        "",
        "| Origin | N | Accuracy | TP Mean | FP Mean | Decoy Mean | FP−Decoy Gap |",
        "|--------|---|----------|---------|---------|------------|--------------|",
    ]

    # Sort by FP-Decoy gap (descending - worst calibration first)
    sorted_results = sorted(results, key=lambda x: x["fp_decoy_gap"], reverse=True)

    for r in sorted_results:
        gap = r["fp_decoy_gap"]
        gap_str = f"+{gap:.1f}" if gap > 0 else f"{gap:.1f}"

        # Highlight hepg2
        if r["origin"] == "hepg2":
            lines.append(
                f"| **`{r['origin']}`** | **{r['n_targets']:,}** | **{r['accuracy_pct']:.1f}%** | "
                f"**{r['tp_mean_conf']:.1f}** | **{r['fp_mean_conf']:.1f}** | "
                f"**{r['decoy_mean_conf']:.1f}** | **{gap_str}** |"
            )
        else:
            lines.append(
                f"| `{r['origin']}` | {r['n_targets']:,} | {r['accuracy_pct']:.1f}% | "
                f"{r['tp_mean_conf']:.1f} | {r['fp_mean_conf']:.1f} | "
                f"{r['decoy_mean_conf']:.1f} | {gap_str} |"
            )

    lines.extend(
        [
            "",
            "**Interpretation:**",
            "- **FP−Decoy Gap > 0**: FPs have higher confidence than decoys → FDR underestimated (dangerous)",
            "- **FP−Decoy Gap ≈ 0**: FPs and decoys have similar confidence → Good calibration",
            "- **FP−Decoy Gap < 0**: Decoys have higher confidence than FPs → FDR overestimated (conservative)",
        ]
    )

    return "\n".join(lines)


def generate_detailed_table(results: list[dict]) -> str:
    """Generate detailed table with medians."""
    lines = [
        "## Detailed Confidence Statistics",
        "",
        "| Origin | TP Mean | TP Med | FP Mean | FP Med | Decoy Mean | Decoy Med |",
        "|--------|---------|--------|---------|--------|------------|-----------|",
    ]

    for r in sorted(results, key=lambda x: x["origin"]):
        lines.append(
            f"| `{r['origin']}` | {r['tp_mean_conf']:.1f} | {r['tp_median_conf']:.1f} | "
            f"{r['fp_mean_conf']:.1f} | {r['fp_median_conf']:.1f} | "
            f"{r['decoy_mean_conf']:.1f} | {r['decoy_median_conf']:.1f} |"
        )

    return "\n".join(lines)


def generate_csv(results: list[dict]) -> str:
    """Generate CSV from results."""
    lines = [
        "origin,masking_rate,n_targets,n_tps,n_fps,n_decoys,accuracy_pct,"
        "tp_mean_conf,fp_mean_conf,decoy_mean_conf,"
        "tp_median_conf,fp_median_conf,decoy_median_conf,fp_decoy_gap"
    ]

    for r in results:
        lines.append(
            f"{r['origin']},{r['masking_rate']:.2f},{r['n_targets']},{r['n_tps']},"
            f"{r['n_fps']},{r['n_decoys']},{r['accuracy_pct']:.2f},"
            f"{r['tp_mean_conf']:.2f},{r['fp_mean_conf']:.2f},{r['decoy_mean_conf']:.2f},"
            f"{r['tp_median_conf']:.2f},{r['fp_median_conf']:.2f},{r['decoy_median_conf']:.2f},"
            f"{r['fp_decoy_gap']:.2f}"
        )

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Extract confidence distributions per origin"
    )
    parser.add_argument(
        "--processed-dir",
        type=Path,
        default=DEFAULT_PROCESSED_DIR,
        help="Directory containing per-origin processed results",
    )
    parser.add_argument(
        "--retention-rate",
        type=float,
        default=0.50,
        help="Retention rate (fraction of peaks KEPT) to analyze. "
        "E.g., 0.10 = 10%% kept = 90%% masked (default: 0.50)",
    )
    parser.add_argument(
        "--all-rates",
        action="store_true",
        help="Analyze all retention rates (0.1 to 0.9, i.e., 90%% to 10%% masked)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis"),
        help="Output directory for results",
    )
    parser.add_argument(
        "--print-tables",
        action="store_true",
        help="Print markdown tables to stdout",
    )
    args = parser.parse_args()

    # Get list of origins
    origins = [d.name for d in args.processed_dir.iterdir() if d.is_dir()]
    print(f"Found {len(origins)} origins: {origins}")

    # Determine masking rates to analyze
    if args.all_rates:
        retention_rates = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    else:
        retention_rates = [args.retention_rate]

    all_results = []

    for rate in retention_rates:
        # rate is the SAMPLING rate (fraction KEPT), so masked = 1 - rate
        masked_pct = int((1 - rate) * 100)
        print(f"\nAnalyzing {masked_pct}% masked (sampling rate {rate:.0%} kept)...")
        results = []

        for origin in sorted(origins):
            stats = extract_confidence_stats(origin, rate, args.processed_dir)
            if stats:
                results.append(stats)
                print(
                    f"  {origin}: Acc={stats['accuracy_pct']:.1f}%, "
                    f"FP−Decoy gap={stats['fp_decoy_gap']:+.1f}"
                )
            else:
                print(f"  {origin}: No data")

        all_results.extend(results)

        if results:
            # Generate outputs for this rate
            args.output_dir.mkdir(parents=True, exist_ok=True)

            md_table = generate_markdown_table(results, rate)
            md_path = (
                args.output_dir / f"confidence_distributions_{masked_pct}pct_masked.md"
            )
            md_path.write_text(md_table)
            print(f"\nMarkdown saved to: {md_path}")

            if args.print_tables:
                print("\n" + "=" * 80)
                print(md_table)
                print("=" * 80)

    # Save combined CSV
    if all_results:
        csv_content = generate_csv(all_results)
        csv_path = (
            args.output_dir / f"confidence_distributions_{masked_pct}pct_masked.csv"
        )
        csv_path.write_text(csv_content)
        print(f"\nCSV saved to: {csv_path}")


if __name__ == "__main__":
    main()
