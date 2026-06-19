#!/usr/bin/env python3
"""Sample spectra from parquet and export to MGF format.

This script reads spectrum data from parquet files (in InstaNovo/Winnow format),
optionally filters by origin, samples a specified number of spectra, and writes
them to MGF format compatible with NovoBoard and other MS tools.

Supports multiple input files via glob patterns.

Usage:
    # Single file
    python scripts/sample_and_export_mgf.py \
        --parquet /path/to/train.parquet \
        --output /path/to/output.mgf \
        --origin hepg2 \
        --sample-n 2000

    # Multiple files via glob
    python scripts/sample_and_export_mgf.py \
        --parquet "/path/to/*.parquet" \
        --output /path/to/output.mgf \
        --origin hepg2 \
        --sample-n 2000

    # Multiple explicit files
    python scripts/sample_and_export_mgf.py \
        --parquet /path/to/train.parquet /path/to/val.parquet /path/to/test.parquet \
        --output /path/to/output.mgf \
        --origin hepg2 \
        --sample-n 2000
"""

from __future__ import annotations

import argparse
import glob
import logging
from pathlib import Path

import polars as pl

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def write_mgf_spectrum(
    f,
    mz_array: list[float],
    intensity_array: list[float],
    precursor_mz: float,
    precursor_charge: int,
    scan_number: int,
    title: str | None = None,
    retention_time: float | None = None,
    peptide_sequence: str | None = None,
) -> None:
    """Write a single spectrum in MGF format.

    Args:
        f: File handle to write to
        mz_array: List of m/z values
        intensity_array: List of intensity values
        precursor_mz: Precursor m/z
        precursor_charge: Precursor charge state
        scan_number: Scan number for identification
        title: Optional spectrum title
        retention_time: Optional retention time in seconds
        peptide_sequence: Optional ground truth peptide sequence (for labeled data)
    """
    f.write("BEGIN IONS\n")

    # Title line
    if title:
        f.write(f"TITLE={title}\n")
    else:
        f.write(f"TITLE=scan={scan_number}\n")

    # Precursor information
    f.write(f"PEPMASS={precursor_mz:.6f}\n")
    f.write(f"CHARGE={precursor_charge}+\n")

    # Retention time if available
    if retention_time is not None:
        f.write(f"RTINSECONDS={retention_time:.2f}\n")

    # Scan number
    f.write(f"SCANS={scan_number}\n")

    # Peptide sequence (ground truth) - required for FDR validation
    if peptide_sequence:
        f.write(f"SEQ={peptide_sequence}\n")

    # Peak list
    for mz, intensity in zip(mz_array, intensity_array):
        f.write(f"{mz:.6f} {intensity:.6f}\n")

    f.write("END IONS\n\n")


def resolve_parquet_paths(parquet_inputs: list[str]) -> list[Path]:
    """Resolve parquet paths, expanding glob patterns.

    Args:
        parquet_inputs: List of paths or glob patterns

    Returns:
        List of resolved Path objects
    """
    resolved: list[Path] = []
    for input_path in parquet_inputs:
        # Check if it's a glob pattern
        if "*" in input_path or "?" in input_path:
            matched = glob.glob(input_path)
            if not matched:
                logger.warning(f"Glob pattern '{input_path}' matched no files")
            resolved.extend(Path(p) for p in sorted(matched))
        else:
            resolved.append(Path(input_path))
    return resolved


def export_mgf(
    parquet_paths: list[Path],
    output_path: Path,
    origin: str | None = None,
    sample_n: int | None = None,
    seed: int = 42,
    shard_idx: int | None = None,
    n_shards: int | None = None,
) -> int:
    """Export spectra from parquet file(s) to MGF format.

    Args:
        parquet_paths: List of paths to input parquet files
        output_path: Path to output MGF file
        origin: Optional origin to filter by
        sample_n: Number of spectra to sample (None = all)
        seed: Random seed for sampling
        shard_idx: 0-based shard index (used together with n_shards)
        n_shards: Total number of shards; rows where row_index % n_shards == shard_idx
                  are selected. Useful for splitting large datasets to avoid RAM OOMs.

    Returns:
        Number of spectra written
    """
    # Validate sharding arguments
    if (shard_idx is None) != (n_shards is None):
        raise ValueError("--shard-idx and --n-shards must be specified together")
    if shard_idx is not None and not (0 <= shard_idx < n_shards):
        raise ValueError(
            f"--shard-idx must be in [0, n_shards): got {shard_idx} with n_shards={n_shards}"
        )

    # Load and concatenate all parquet files
    dfs = []
    for parquet_path in parquet_paths:
        logger.info(f"Reading parquet file: {parquet_path}")
        df = pl.read_parquet(parquet_path)
        logger.info(f"  Loaded {len(df):,} spectra")
        dfs.append(df)

    df = pl.concat(dfs)
    logger.info(f"Total: {len(df):,} spectra from {len(parquet_paths)} file(s)")

    # Filter by origin if specified
    if origin:
        df = df.filter(pl.col("origin") == origin)
        logger.info(f"Filtered to {len(df):,} spectra from origin '{origin}'")

        if len(df) == 0:
            logger.error(f"No spectra found for origin '{origin}'")
            # Get available origins from first file
            available_origins = (
                pl.read_parquet(parquet_paths[0])["origin"].unique().sort().to_list()
            )
            logger.error(f"Available origins: {available_origins}")
            return 0

    # Sample if requested
    if sample_n is not None and sample_n < len(df):
        df = df.sample(n=sample_n, seed=seed)
        logger.info(f"Sampled {len(df):,} spectra")
    elif sample_n is not None:
        logger.info(
            f"Requested {sample_n:,} but only {len(df):,} available - using all"
        )

    # Apply sharding: select every n_shards-th row starting at shard_idx
    if shard_idx is not None:
        df = (
            df.with_row_index("__shard_key__")
            .filter(pl.col("__shard_key__") % n_shards == shard_idx)
            .drop("__shard_key__")
        )
        logger.info(f"Shard {shard_idx + 1}/{n_shards}: {len(df):,} spectra")

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write MGF file
    logger.info(f"Writing MGF file: {output_path}")
    count = 0

    with open(output_path, "w") as f:
        for idx, row in enumerate(df.iter_rows(named=True)):
            mz_array = row["mz_array"]
            intensity_array = row["intensity_array"]

            # Handle polars Series if needed (shouldn't happen with iter_rows but just in case)
            if hasattr(mz_array, "to_list"):
                mz_array = mz_array.to_list()
            if hasattr(intensity_array, "to_list"):
                intensity_array = intensity_array.to_list()

            # Skip empty spectra
            if not mz_array or not intensity_array:
                continue

            # Build title matching InstaNovo's GNPS export style
            # Note: We use original scan_number in title for reference, but SCANS= uses
            # the 0-based index to match InstaNovo's scan_number output
            title = f"{row['experiment_name']} scan={row['scan_number']}"

            # Get peptide sequence for labeled data (ground truth for FDR validation)
            peptide_sequence = row.get("sequence", "")

            write_mgf_spectrum(
                f,
                mz_array=mz_array,
                intensity_array=intensity_array,
                precursor_mz=row["precursor_mz"],
                precursor_charge=row["precursor_charge"],
                # Use 0-based index as scan_number to match InstaNovo's output format
                # InstaNovo outputs scan_number as the spectrum position in the MGF, not SCANS=
                scan_number=idx,
                title=title,
                retention_time=row.get("retention_time"),
                peptide_sequence=peptide_sequence,
            )
            count += 1

    logger.info(f"Wrote {count:,} spectra to {output_path}")
    return count


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Sample spectra from parquet and export to MGF format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Export from single file
    python scripts/sample_and_export_mgf.py \\
        --parquet /path/to/train.parquet \\
        --output /path/to/hepg2.mgf \\
        --origin hepg2 \\
        --sample-n 2000

    # Export from multiple files (glob pattern)
    python scripts/sample_and_export_mgf.py \\
        --parquet "/path/to/*.parquet" \\
        --output /path/to/hepg2.mgf \\
        --origin hepg2 \\
        --sample-n 2000

    # Export from multiple explicit files
    python scripts/sample_and_export_mgf.py \\
        --parquet /path/to/train.parquet /path/to/val.parquet /path/to/test.parquet \\
        --output /path/to/all.mgf

    # List available origins
    python scripts/sample_and_export_mgf.py \\
        --parquet /path/to/train.parquet \\
        --list-origins
        """,
    )

    parser.add_argument(
        "--parquet",
        type=str,
        nargs="+",
        required=True,
        help="Path(s) to input parquet file(s). Supports glob patterns (e.g., '*.parquet')",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Path to output MGF file",
    )
    parser.add_argument(
        "--origin",
        type=str,
        default=None,
        help="Filter to spectra from this origin",
    )
    parser.add_argument(
        "--sample-n",
        type=int,
        default=None,
        help="Number of spectra to sample (default: all)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for sampling (default: 42)",
    )
    parser.add_argument(
        "--list-origins",
        action="store_true",
        help="List available origins and exit",
    )
    parser.add_argument(
        "--n-shards",
        type=int,
        default=None,
        help="Split the dataset into this many shards (use with --shard-idx). "
        "Rows where row_index %% n_shards == shard_idx are written. "
        "Useful for large datasets to avoid RAM OOMs during decoy generation.",
    )
    parser.add_argument(
        "--shard-idx",
        type=int,
        default=None,
        help="0-based index of the shard to export (use with --n-shards).",
    )

    args = parser.parse_args()

    # Resolve parquet paths (expand globs)
    parquet_paths = resolve_parquet_paths(args.parquet)

    if not parquet_paths:
        parser.error(f"No parquet files found matching: {args.parquet}")

    logger.info(f"Found {len(parquet_paths)} parquet file(s)")

    # List origins mode
    if args.list_origins:
        # Load all files and aggregate
        dfs = [pl.read_parquet(p) for p in parquet_paths]
        df = pl.concat(dfs)
        origins = df.group_by("origin").len().sort("origin")
        print(
            f"\nAvailable origins (from {len(parquet_paths)} file(s), {len(df):,} total spectra):"
        )
        for row in origins.iter_rows(named=True):
            print(f"  {row['origin']}: {row['len']:,} spectra")
        return

    # Normal export mode
    if args.output is None:
        parser.error("--output is required unless using --list-origins")

    count = export_mgf(
        parquet_paths=parquet_paths,
        output_path=args.output,
        origin=args.origin,
        sample_n=args.sample_n,
        seed=args.seed,
        shard_idx=args.shard_idx,
        n_shards=args.n_shards,
    )

    if count == 0:
        exit(1)


if __name__ == "__main__":
    main()
