"""Command-line interface for NovoBoard."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys
from typing import Sequence

from novoboard import config
from novoboard.accuracy import WorkerTest
from novoboard.decoy import generate_decoy_mgf
from novoboard.fdr import validate_FDR
from novoboard.plotting import plot_fdr_validation
from novoboard.preprocessing import run_preprocessing

logger = logging.getLogger(__name__)


def download_data(data_dir: Path) -> None:
    """Download example ABRF data from Google Drive using gdown.

    Args:
        data_dir: Directory to download data into
    """
    import gdown

    folder_url = (
        "https://drive.google.com/drive/folders/1_6azR4-YjTUfRYdsXbFhZL9lFvjdrIDh"
    )

    data_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Downloading data to {data_dir}...")
    gdown.download_folder(folder_url, output=str(data_dir), quiet=False)
    logger.info("Download complete.")


def run_accuracy(
    db_file: Path,
    denovo_file: Path,
    spectrum_file: Path,
    col_score: str,
    col_aa_score: str,
) -> None:
    """Calculate accuracy of de novo predictions against database search results.

    Args:
        db_file: Path to database search results CSV
        denovo_file: Path to de novo sequencing results CSV
        spectrum_file: Path to MGF spectrum file
        col_score: Column name for peptide score
        col_aa_score: Column name for AA-level scores
    """
    logger.info(f"Database file: {db_file}")
    logger.info(f"De novo file: {denovo_file}")
    logger.info(f"Spectrum file: {spectrum_file}")

    worker_test = WorkerTest(
        str(db_file),
        str(denovo_file),
        str(spectrum_file),
        col_score,
        col_aa_score,
    )
    worker_test.test_accuracy()


def run_decoy_generation(
    spectrum_files: Sequence[Path],
    peak_sampling: str = "random",
    sampling_rate: float = config.DEFAULT_SAMPLING_RATE,
    seed: int = 99,
) -> None:
    """Generate decoy MGF files for FDR estimation.

    Args:
        spectrum_files: List of input MGF spectrum files
        peak_sampling: Peak sampling strategy
        sampling_rate: Fraction of peaks to sample
        seed: Random seed for reproducibility
    """
    logger.info(f"Generating decoys for {len(spectrum_files)} spectrum file(s)")
    input_mgf_str_list = [str(p) for p in spectrum_files]
    generate_decoy_mgf(input_mgf_str_list, peak_sampling, sampling_rate, seed)


def extract_decoy_label(filename: str) -> str:
    """Extract decoy percentage label from filename.

    Looks for patterns like 'decoy_0.10' or 'decoy_10' in the filename
    and converts to a human-readable percentage label like '10%'.
    Falls back to the filename stem if no pattern is found.
    """
    import re

    # Match patterns like decoy_0.10, decoy_0.5, decoy_10, etc.
    match = re.search(r"decoy_(\d+\.?\d*)", filename, re.IGNORECASE)
    if match:
        value = float(match.group(1))
        # If value is less than 1, assume it's a fraction (0.10 = 10%)
        if value < 1:
            return f"{int(value * 100)}%"
        else:
            return f"{int(value)}%"
    return Path(filename).stem


def run_fdr_validation(
    target_file: Path,
    decoy_files: Sequence[Path],
    db_file: Path,
    spectrum_file: Path,
    output_file: Path,
    col_score: str,
    col_aa_score: str,
    ion_threshold: float = 0.90,
    labels: Sequence[str] | None = None,
    fdr_max: float = 0.05,
    dpi: int = 150,
    monotonic: bool = True,
) -> None:
    """Validate FDR estimation using target-decoy approach.

    Args:
        target_file: Path to target de novo results CSV
        decoy_files: Paths to decoy de novo results CSVs
        db_file: Path to database search results CSV (ground truth)
        spectrum_file: Path to MGF spectrum file
        output_file: Path for output plot
        col_score: Column name for peptide score
        col_aa_score: Column name for AA-level scores
        ion_threshold: Ion matching threshold percentage (default: 0.90)
        labels: Custom labels for each decoy (if None, extracted from filenames)
        fdr_max: Maximum value for FDR axis in plot (default: 0.05)
        dpi: Resolution in dots per inch (default: 150)
        monotonic: If True, filter to monotonically decreasing FDR (default: True)
    """
    logger.info(f"Target file: {target_file}")
    logger.info(f"Decoy files: {len(decoy_files)}")
    logger.info(f"Database file: {db_file}")

    # Generate FDR thresholds up to fdr_max (with 0.001 step size)
    p_decoy = [x / 1000.0 for x in range(0, int(fdr_max * 1000) + 1, 1)]

    results_list = [
        validate_FDR(
            str(target_file),
            str(decoy_file),
            col_score,  # engine_score
            str(db_file),
            str(spectrum_file),
            p_decoy,
            ion_threshold,
            col_score,
            col_aa_score,
            monotonic=monotonic,
        )
        for decoy_file in decoy_files
    ]

    # Use provided labels or extract from filenames
    if labels is None:
        labels = [extract_decoy_label(str(decoy_file)) for decoy_file in decoy_files]
    samples = range(1, len(decoy_files) + 1)
    plot_fdr_validation(
        results_list, samples, str(output_file), labels=labels, fdr_max=fdr_max, dpi=dpi
    )


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for CLI usage.

    Args:
        verbose: If True, set DEBUG level; otherwise INFO level
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def main() -> None:
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        description="NovoBoard - Framework for evaluating de novo peptide sequencing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Preprocess InstaNovo predictions for NovoBoard
  novoboard preprocess --denovo-file data/instanovo_preds.csv --denovo-output results/denovo.csv

  # Preprocess both de novo and database files
  novoboard preprocess --denovo-file data/preds.csv --denovo-output results/denovo.csv \\
                       --db-mgf-file data/labeled.mgf --db-output results/db.csv

  # Calculate accuracy of de novo predictions
  novoboard accuracy --db-file db_results.csv --denovo-file denovo.csv --spectrum-file spectra.mgf

  # Generate decoy spectra
  novoboard decoy --spectrum-file spectra.mgf --sampling-rate 0.5

  # Validate FDR estimation
  novoboard fdr --target-file target.csv --decoy-files decoy1.csv decoy2.csv \\
                --db-file db_results.csv --spectrum-file spectra.mgf

  # Download example ABRF dataset
  novoboard download --output-dir data
        """,
    )

    # Global arguments
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose (debug) logging"
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # =========================================================================
    # DOWNLOAD command
    # =========================================================================
    download_parser = subparsers.add_parser(
        "download", help="Download example ABRF dataset from Google Drive"
    )
    download_parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data"),
        help="Directory to download data into (default: data/)",
    )

    # =========================================================================
    # ACCURACY command
    # =========================================================================
    accuracy_parser = subparsers.add_parser(
        "accuracy",
        help="Calculate accuracy of de novo predictions against database search",
    )
    accuracy_parser.add_argument(
        "--db-file",
        type=Path,
        required=True,
        help="Path to database search results CSV",
    )
    accuracy_parser.add_argument(
        "--denovo-file",
        type=Path,
        required=True,
        help="Path to de novo sequencing results CSV",
    )
    accuracy_parser.add_argument(
        "--spectrum-file", type=Path, required=True, help="Path to MGF spectrum file"
    )
    accuracy_parser.add_argument(
        "--score-column",
        type=str,
        default="ALC (%)",
        help='Column name for peptide score (default: "ALC (%%)")',
    )
    accuracy_parser.add_argument(
        "--aa-score-column",
        type=str,
        default="local confidence (%)",
        help='Column name for AA-level scores (default: "local confidence (%%)")',
    )

    # =========================================================================
    # DECOY command
    # =========================================================================
    decoy_parser = subparsers.add_parser(
        "decoy", help="Generate decoy MGF files for FDR estimation"
    )
    decoy_parser.add_argument(
        "--spectrum-file",
        type=Path,
        required=True,
        nargs="+",
        help="Path(s) to input MGF spectrum file(s)",
    )
    decoy_parser.add_argument(
        "--sampling-strategy",
        type=str,
        default="random",
        choices=[
            "random",
            "intensity",
            "intensity_mass",
            "permutation",
            "500Da",
            "distance",
        ],
        help="Peak sampling strategy (default: random)",
    )
    decoy_parser.add_argument(
        "--sampling-rate",
        type=float,
        default=config.DEFAULT_SAMPLING_RATE,
        help=f"Fraction of peaks to keep (default: {config.DEFAULT_SAMPLING_RATE})",
    )
    decoy_parser.add_argument(
        "--seed",
        type=int,
        default=99,
        help="Random seed for reproducibility (default: 99)",
    )

    # =========================================================================
    # FDR command
    # =========================================================================
    fdr_parser = subparsers.add_parser(
        "fdr", help="Validate FDR estimation using target-decoy approach"
    )
    fdr_parser.add_argument(
        "--target-file",
        type=Path,
        required=True,
        help="Path to target de novo results CSV",
    )
    fdr_parser.add_argument(
        "--decoy-files",
        type=Path,
        required=True,
        nargs="+",
        help="Path(s) to decoy de novo results CSV(s)",
    )
    fdr_parser.add_argument(
        "--db-file",
        type=Path,
        required=True,
        help="Path to database search results CSV (ground truth)",
    )
    fdr_parser.add_argument(
        "--spectrum-file", type=Path, required=True, help="Path to MGF spectrum file"
    )
    fdr_parser.add_argument(
        "--output-file",
        type=Path,
        default=Path("fdr_validation.png"),
        help="Path for output plot (default: fdr_validation.png)",
    )
    fdr_parser.add_argument(
        "--score-column",
        type=str,
        default="ALC (%)",
        help='Column name for peptide score (default: "ALC (%%)")',
    )
    fdr_parser.add_argument(
        "--aa-score-column",
        type=str,
        default="local confidence (%)",
        help='Column name for AA-level scores (default: "local confidence (%%)")',
    )
    fdr_parser.add_argument(
        "--ion-threshold",
        type=float,
        default=0.90,
        help="Ion matching threshold percentage (default: 0.90)",
    )
    fdr_parser.add_argument(
        "--labels",
        type=str,
        nargs="+",
        help="Labels for each decoy file (in same order as --decoy-files). If not provided, extracts from filenames.",
    )

    def fdr_max_type(value: str) -> float:
        """Validate FDR max is between 0 and 1."""
        fval = float(value)
        if not 0 < fval <= 1.0:
            raise argparse.ArgumentTypeError(f"must be between 0 and 1, got {fval}")
        return fval

    fdr_parser.add_argument(
        "--fdr-max",
        type=fdr_max_type,
        default=0.05,
        metavar="FDR_MAX",
        help="Maximum FDR value for plot axes, between 0 and 1 (default: 0.05)",
    )
    fdr_parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Plot resolution in dots per inch (default: 150, use 300 for print quality)",
    )
    fdr_parser.add_argument(
        "--no-monotonic",
        action="store_true",
        help="Disable monotonic filtering to show all FDR data points (may look bumpy)",
    )

    # =========================================================================
    # PREPROCESS command
    # =========================================================================
    preprocess_parser = subparsers.add_parser(
        "preprocess", help="Convert InstaNovo output to NovoBoard format"
    )
    preprocess_parser.add_argument(
        "--denovo-file", type=Path, help="Path to InstaNovo predictions CSV file"
    )
    preprocess_parser.add_argument(
        "--denovo-output", type=Path, help="Path for de novo results output CSV"
    )
    preprocess_parser.add_argument(
        "--db-mgf-file",
        type=Path,
        help="Path to labeled MGF file (for extracting database annotations)",
    )
    preprocess_parser.add_argument(
        "--db-output", type=Path, help="Path for database results output CSV"
    )

    # Parse arguments
    args = parser.parse_args()

    # Set up logging
    setup_logging(verbose=args.verbose)

    # Handle no command
    if args.command is None:
        parser.print_help()
        sys.exit(0)

    # Execute command
    if args.command == "download":
        download_data(args.output_dir)

    elif args.command == "accuracy":
        # Validate files exist
        for f in [args.db_file, args.denovo_file, args.spectrum_file]:
            if not f.exists():
                logger.error(f"File not found: {f}")
                sys.exit(1)

        run_accuracy(
            args.db_file,
            args.denovo_file,
            args.spectrum_file,
            args.score_column,
            args.aa_score_column,
        )

    elif args.command == "decoy":
        # Validate files exist
        for f in args.spectrum_file:
            if not f.exists():
                logger.error(f"File not found: {f}")
                sys.exit(1)

        run_decoy_generation(
            args.spectrum_file,
            args.sampling_strategy,
            args.sampling_rate,
            args.seed,
        )

    elif args.command == "fdr":
        # Validate files exist
        files_to_check = [args.target_file, args.db_file, args.spectrum_file] + list(
            args.decoy_files
        )
        for f in files_to_check:
            if not f.exists():
                logger.error(f"File not found: {f}")
                sys.exit(1)

        # Create output directory if needed
        args.output_file.parent.mkdir(parents=True, exist_ok=True)

        run_fdr_validation(
            args.target_file,
            args.decoy_files,
            args.db_file,
            args.spectrum_file,
            args.output_file,
            args.score_column,
            args.aa_score_column,
            args.ion_threshold,
            labels=getattr(args, "labels", None),
            fdr_max=args.fdr_max,
            dpi=args.dpi,
            monotonic=not args.no_monotonic,
        )

    elif args.command == "preprocess":
        # Validate at least one input/output pair is provided
        if not (args.denovo_file or args.db_mgf_file):
            logger.error("Must specify at least --denovo-file or --db-mgf-file")
            sys.exit(1)

        if args.denovo_file and not args.denovo_output:
            logger.error("--denovo-output is required when --denovo-file is specified")
            sys.exit(1)

        if args.db_mgf_file and not args.db_output:
            logger.error("--db-output is required when --db-mgf-file is specified")
            sys.exit(1)

        # Validate input files exist
        if args.denovo_file and not args.denovo_file.exists():
            logger.error(f"File not found: {args.denovo_file}")
            sys.exit(1)
        if args.db_mgf_file and not args.db_mgf_file.exists():
            logger.error(f"File not found: {args.db_mgf_file}")
            sys.exit(1)

        run_preprocessing(
            denovo_file=args.denovo_file,
            denovo_output=args.denovo_output,
            db_mgf_file=args.db_mgf_file,
            db_output=args.db_output,
        )

    logger.info("NovoBoard analysis complete!")


if __name__ == "__main__":
    main()
