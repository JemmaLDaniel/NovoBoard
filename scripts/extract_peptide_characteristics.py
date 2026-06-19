#!/usr/bin/env python3
"""Extract peptide characteristics per origin from parquet data.

Computes:
- Mean peptide length (amino acids only, modifications stripped)
- Charge state distribution (z=2, z=3, z≥4)
- Sample size

Usage:
    python scripts/extract_peptide_characteristics.py --parquet-files /path/to/*.parquet
    python scripts/extract_peptide_characteristics.py  # uses default path
"""

import argparse
import re
from pathlib import Path

import polars as pl


def count_amino_acids(sequence: str) -> int:
    """Count actual amino acids, stripping modification annotations.

    Handles formats like:
    - PEPTIDEC(+57.02)M(+15.99)K -> 10 aa
    - PEPTIDEC(Carbamidomethylation)K -> 10 aa
    - C[UNIMOD:4]PEPTIDE -> 8 aa
    - PEPTIDE -> 7 aa
    """
    if sequence is None:
        return 0
    # Remove modification annotations:
    # - Parentheses: (+15.99), (Carbamidomethylation)
    # - Brackets: [UNIMOD:4], [+57.02]
    stripped = re.sub(r"\([^)]*\)", "", sequence)
    stripped = re.sub(r"\[[^\]]*\]", "", stripped)
    return len(stripped)


def compute_characteristics(df: pl.DataFrame) -> dict:
    """Compute peptide characteristics for a DataFrame."""
    n = len(df)
    if n == 0:
        return {
            "n": 0,
            "mean_length": None,
            "z2_pct": None,
            "z3_pct": None,
            "z4plus_pct": None,
        }

    # Compute amino acid lengths (stripping modifications)
    lengths = [count_amino_acids(seq) for seq in df["sequence"].to_list()]
    mean_length = sum(lengths) / len(lengths)

    # Charge distribution
    charges = df["precursor_charge"].to_list()
    z2 = sum(1 for c in charges if c == 2) / n * 100
    z3 = sum(1 for c in charges if c == 3) / n * 100
    z4plus = sum(1 for c in charges if c >= 4) / n * 100

    return {
        "n": n,
        "mean_length": mean_length,
        "z2_pct": z2,
        "z3_pct": z3,
        "z4plus_pct": z4plus,
    }


def generate_markdown_table(results: dict[str, dict]) -> str:
    """Generate markdown table from results."""
    lines = [
        "## Peptide Characteristics by Origin",
        "",
        "| Origin | N | Mean Length | Doubly-Charged (z=2) | Triply-Charged (z=3) | Higher (z≥4) |",
        "|--------|---|-------------|----------------------|----------------------|--------------|",
    ]

    # Sort by origin name
    for origin in sorted(results.keys()):
        stats = results[origin]
        n = stats["n"]
        mean_len = stats["mean_length"]
        z2 = stats["z2_pct"]
        z3 = stats["z3_pct"]
        z4plus = stats["z4plus_pct"]

        if mean_len is None:
            lines.append(f"| `{origin}` | {n:,} | — | — | — | — |")
        else:
            lines.append(
                f"| `{origin}` | {n:,} | {mean_len:.1f} aa | {z2:.1f}% | {z3:.1f}% | {z4plus:.1f}% |"
            )

    return "\n".join(lines)


def generate_csv(results: dict[str, dict]) -> str:
    """Generate CSV from results."""
    lines = ["origin,n,mean_length,z2_pct,z3_pct,z4plus_pct"]

    for origin in sorted(results.keys()):
        stats = results[origin]
        n = stats["n"]
        mean_len = stats["mean_length"]
        z2 = stats["z2_pct"]
        z3 = stats["z3_pct"]
        z4plus = stats["z4plus_pct"]

        if mean_len is not None:
            lines.append(f"{origin},{n},{mean_len:.2f},{z2:.2f},{z3:.2f},{z4plus:.2f}")
        else:
            lines.append(f"{origin},{n},,,")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Extract peptide characteristics per origin"
    )
    parser.add_argument(
        "--parquet-files",
        nargs="+",
        default=["/home/j-daniel/Documents/winnow/data/paper_datasets/train.parquet"],
        help="Parquet files to analyze (supports glob patterns)",
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

    # Load all parquet files
    dfs = []
    for pattern in args.parquet_files:
        path = Path(pattern)
        if path.exists():
            print(f"Loading: {path}")
            dfs.append(pl.read_parquet(path))
        else:
            # Try as glob
            for f in Path(pattern).parent.glob(Path(pattern).name):
                print(f"Loading: {f}")
                dfs.append(pl.read_parquet(f))

    if not dfs:
        print("No parquet files found!")
        return

    df = pl.concat(dfs)
    print(f"Total spectra: {len(df):,}")

    # Check required columns
    required = ["origin", "sequence", "precursor_charge"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"Missing required columns: {missing}")
        print(f"Available columns: {df.columns}")
        return

    # Compute characteristics per origin
    results = {}
    for origin in sorted(df["origin"].unique().to_list()):
        origin_df = df.filter(pl.col("origin") == origin)
        results[origin] = compute_characteristics(origin_df)
        print(
            f"  {origin}: N={results[origin]['n']:,}, mean_len={results[origin]['mean_length']:.1f} aa"
        )

    # Generate outputs
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Markdown table
    md_table = generate_markdown_table(results)
    md_path = args.output_dir / "peptide_characteristics.md"
    md_path.write_text(md_table)
    print(f"\nMarkdown table saved to: {md_path}")

    # CSV
    csv_content = generate_csv(results)
    csv_path = args.output_dir / "peptide_characteristics.csv"
    csv_path.write_text(csv_content)
    print(f"CSV saved to: {csv_path}")

    if args.print_tables:
        print("\n" + "=" * 80)
        print(md_table)
        print("=" * 80)


if __name__ == "__main__":
    main()
