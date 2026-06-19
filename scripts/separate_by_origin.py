#!/usr/bin/env python3
"""Separate train_sample10 data by origin for per-origin FDR analysis.

This script:
1. Reads the original parquet to get the origin column
2. Joins with preprocessed NovoBoard CSVs on Scan = row_index
3. Filters all files (de novo, database, accuracy) by origin
4. Saves filtered files to per-origin folders with proper naming for cache hits
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import polars as pl

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Default paths
DEFAULT_PARQUET = (
    "/home/j-daniel/Documents/winnow/novoboard/benchmark_data/train_sample10.parquet"
)
DEFAULT_INPUT_DIR = (
    "/home/j-daniel/Documents/winnow/novoboard/benchmark_data/novoboard_inputs"
)
DEFAULT_OUTPUT_DIR = (
    "/home/j-daniel/Documents/winnow/novoboard/benchmark_data/per_origin"
)

# Decoy rates to process
DECOY_RATES = ["0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70", "0.80", "0.90"]


def load_origin_mapping(parquet_path: Path) -> pl.DataFrame:
    """Load parquet and create row_idx to origin mapping."""
    logger.info(f"Loading parquet: {parquet_path}")
    df = pl.read_parquet(parquet_path)
    df = df.with_row_index("row_idx")
    logger.info(f"  Loaded {len(df)} rows")
    logger.info(f"  Unique origins: {df['origin'].unique().to_list()}")
    return df.select(["row_idx", "origin"])


def filter_denovo_by_origin(
    csv_path: Path,
    origin_mapping: pl.DataFrame,
    origin: str,
) -> pl.DataFrame:
    """Filter de novo results CSV by origin using Scan = row_idx."""
    df = pl.read_csv(csv_path)
    df_with_origin = df.join(
        origin_mapping,
        left_on="Scan",
        right_on="row_idx",
        how="left",
    )
    filtered = df_with_origin.filter(pl.col("origin") == origin).drop("origin")
    return filtered


def filter_database_by_origin(
    csv_path: Path,
    origin_mapping: pl.DataFrame,
    origin: str,
) -> pl.DataFrame:
    """Filter database results CSV by origin using Scan = row_idx."""
    df = pl.read_csv(csv_path)
    df_with_origin = df.join(
        origin_mapping,
        left_on="Scan",
        right_on="row_idx",
        how="left",
    )
    filtered = df_with_origin.filter(pl.col("origin") == origin).drop("origin")
    return filtered


def filter_accuracy_by_origin(
    csv_path: Path,
    origin_mapping: pl.DataFrame,
    origin: str,
) -> pl.DataFrame:
    """Filter accuracy CSV by origin using feature_id -> Scan extraction."""
    # Accuracy files are TSV
    df = pl.read_csv(csv_path, separator="\t")

    # Extract scan from feature_id (format: "train_sample10||{scan}")
    df = df.with_columns(
        pl.col("feature_id")
        .str.extract(r"\|\|(\d+)$")
        .cast(pl.Int64)
        .alias("extracted_scan")
    )

    df_with_origin = df.join(
        origin_mapping,
        left_on="extracted_scan",
        right_on="row_idx",
        how="left",
    )

    filtered = df_with_origin.filter(pl.col("origin") == origin).drop(
        ["origin", "extracted_scan"]
    )
    return filtered


def update_source_file_column(df: pl.DataFrame, new_source: str) -> pl.DataFrame:
    """Update the Source File column to match the new origin-based naming."""
    if "Source File" in df.columns:
        df = df.with_columns(pl.lit(new_source).alias("Source File"))
    return df


def update_feature_id_column(df: pl.DataFrame, new_source: str) -> pl.DataFrame:
    """Update the feature_id column to use the new source name."""
    if "feature_id" in df.columns:
        # Replace the source part of feature_id (before ||)
        df = df.with_columns(
            (
                pl.lit(new_source)
                + "||"
                + pl.col("feature_id").str.extract(r"\|\|(\d+)$")
            ).alias("feature_id")
        )
    return df


def process_origin(
    origin: str,
    origin_mapping: pl.DataFrame,
    input_dir: Path,
    output_dir: Path,
) -> None:
    """Process all files for a single origin."""
    logger.info(f"\n{'='*60}")
    logger.info(f"Processing origin: {origin}")
    logger.info(f"{'='*60}")

    origin_output_dir = output_dir / origin
    origin_output_dir.mkdir(parents=True, exist_ok=True)

    # Filter target de novo results
    target_input = input_dir / "train_sample10_de_novo_results.csv"
    target_output = origin_output_dir / f"{origin}_de_novo_results.csv"

    if target_input.exists():
        df = filter_denovo_by_origin(target_input, origin_mapping, origin)
        df = update_source_file_column(df, origin)
        df.write_csv(target_output)
        logger.info(f"  Target de novo: {len(df)} rows -> {target_output.name}")
    else:
        logger.warning(f"  Target file not found: {target_input}")

    # Filter target database results
    db_input = input_dir / "train_sample10_database_results.csv"
    db_output = origin_output_dir / f"{origin}_database_results.csv"

    if db_input.exists():
        df = filter_database_by_origin(db_input, origin_mapping, origin)
        df = update_source_file_column(df, origin)
        df.write_csv(db_output)
        logger.info(f"  Database results: {len(df)} rows -> {db_output.name}")
    else:
        logger.warning(f"  Database file not found: {db_input}")

    # Filter decoy de novo results and their corresponding files
    for rate in DECOY_RATES:
        decoy_denovo_input = (
            input_dir / f"train_sample10_decoy_{rate}_de_novo_results.csv"
        )
        decoy_denovo_output = (
            origin_output_dir / f"{origin}_decoy_{rate}_de_novo_results.csv"
        )

        if decoy_denovo_input.exists():
            df = filter_denovo_by_origin(decoy_denovo_input, origin_mapping, origin)
            df = update_source_file_column(df, f"{origin}_decoy_{rate}")
            df.write_csv(decoy_denovo_output)
            logger.info(
                f"  Decoy {rate} de novo: {len(df)} rows -> {decoy_denovo_output.name}"
            )

        # Filter corresponding accuracy file
        accuracy_input = (
            input_dir
            / f"train_sample10_de_novo_results_fdr_train_sample10_decoy_{rate}_de_novo_results_accuracy.csv"
        )
        accuracy_output = (
            origin_output_dir
            / f"{origin}_de_novo_results_fdr_{origin}_decoy_{rate}_de_novo_results_accuracy.csv"
        )

        if accuracy_input.exists():
            df = filter_accuracy_by_origin(accuracy_input, origin_mapping, origin)
            df = update_feature_id_column(df, origin)
            # Write as TSV to match original format
            df.write_csv(accuracy_output, separator="\t")
            logger.info(f"  Accuracy {rate}: {len(df)} rows -> {accuracy_output.name}")


def validate_row_counts(
    parquet_path: Path,
    input_dir: Path,
) -> None:
    """Validate that parquet and de novo results have compatible row counts."""
    parquet_df = pl.read_parquet(parquet_path)
    parquet_rows = len(parquet_df)

    denovo_path = input_dir / "train_sample10_de_novo_results.csv"
    if not denovo_path.exists():
        raise FileNotFoundError(f"De novo results not found: {denovo_path}")

    denovo_df = pl.read_csv(denovo_path)
    max_scan = denovo_df["Scan"].max()

    # The Scan values should be valid indices into the parquet
    if max_scan >= parquet_rows:
        raise ValueError(
            f"Scan values in de novo results (max={max_scan}) exceed "
            f"parquet row count ({parquet_rows}). Row alignment may be incorrect."
        )

    logger.info(
        f"Validation passed: Scan range [0, {max_scan}] fits within parquet ({parquet_rows} rows)"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Separate train_sample10 data by origin for per-origin FDR analysis"
    )
    parser.add_argument(
        "--parquet",
        type=Path,
        default=Path(DEFAULT_PARQUET),
        help=f"Path to original parquet file (default: {DEFAULT_PARQUET})",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path(DEFAULT_INPUT_DIR),
        help=f"Directory containing preprocessed NovoBoard files (default: {DEFAULT_INPUT_DIR})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(DEFAULT_OUTPUT_DIR),
        help=f"Output directory for per-origin files (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--origins",
        nargs="+",
        default=None,
        help="Specific origins to process (default: all found in parquet)",
    )

    args = parser.parse_args()

    # Validate inputs exist
    if not args.parquet.exists():
        logger.error(f"Parquet file not found: {args.parquet}")
        sys.exit(1)
    if not args.input_dir.exists():
        logger.error(f"Input directory not found: {args.input_dir}")
        sys.exit(1)

    # Validate row alignment
    validate_row_counts(args.parquet, args.input_dir)

    # Load origin mapping
    origin_mapping = load_origin_mapping(args.parquet)

    # Get list of origins to process
    if args.origins:
        origins = args.origins
    else:
        origins = origin_mapping["origin"].unique().to_list()

    logger.info(f"\nWill process {len(origins)} origins: {origins}")

    # Process each origin
    for origin in origins:
        process_origin(origin, origin_mapping, args.input_dir, args.output_dir)

    logger.info(f"\n{'='*60}")
    logger.info("Done! Per-origin files saved to:")
    logger.info(f"  {args.output_dir}")
    logger.info(f"{'='*60}")


if __name__ == "__main__":
    main()
