"""Preprocessing utilities for converting InstaNovo outputs to NovoBoard format."""

from __future__ import annotations

import logging
from pathlib import Path

import polars as pl
from matchms.importing import load_from_mgf

logger = logging.getLogger(__name__)

# Column mapping from InstaNovo format to NovoBoard format
COLUMN_MAPPING = {
    "experiment_name": "Source File",
    "scan_number": "Scan",
    "predictions": "Peptide",
    "token_log_probs": "local confidence (%)",
}

# UNIMOD to NovoBoard modification format mapping
UNIMOD_TO_NOVOBOARD = {
    "C[UNIMOD:4]": "C(+57.02)",
    "M[UNIMOD:35]": "M(+15.99)",
    "N[UNIMOD:7]": "N(+0.98)",
    "Q[UNIMOD:7]": "Q(+0.98)",
    "S[UNIMOD:21]": "S(+79.97)",
    "T[UNIMOD:21]": "T(+79.97)",
    "Y[UNIMOD:21]": "Y(+79.97)",
}


def add_score_column(data: pl.DataFrame) -> pl.DataFrame:
    """Create a score column from InstaNovo log-probability (0-100 scale).

    Args:
        data: DataFrame with 'log_probs' column

    Returns:
        DataFrame with 'ALC (%)' column added
    """
    data = data.with_columns(
        pl.col("log_probs").cast(pl.Float64).exp().mul(100).round(1).alias("ALC (%)")
    )
    return data


def convert_unimod_to_novoboard(peptide: str) -> str:
    """Convert InstaNovo UNIMOD notation to NovoBoard format.

    Args:
        peptide: Peptide sequence with UNIMOD modifications

    Returns:
        Peptide sequence with NovoBoard-style modifications
    """
    result = peptide

    # Replace supported residue modifications (works for multiple mods)
    for unimod, novoboard in UNIMOD_TO_NOVOBOARD.items():
        result = result.replace(unimod, novoboard)

    return result


def map_supported_modifications(data: pl.DataFrame) -> pl.DataFrame:
    """Convert UNIMOD notation to NovoBoard format for all peptides.

    Args:
        data: DataFrame with 'Peptide' column

    Returns:
        DataFrame with converted modification notation
    """
    data = data.with_columns(
        pl.col("Peptide")
        .map_elements(convert_unimod_to_novoboard, return_dtype=pl.Utf8)
        .alias("Peptide")
    )
    return data


def filter_unsupported_modifications(data: pl.DataFrame) -> pl.DataFrame:
    """Filter out peptides with unsupported UNIMOD modifications.

    Args:
        data: DataFrame with 'Peptide' column

    Returns:
        DataFrame with unsupported modifications removed
    """
    initial_count = len(data)
    data = data.filter(~pl.col("Peptide").str.contains(r"\[UNIMOD:\d+\]"))
    filtered_count = initial_count - len(data)
    if filtered_count > 0:
        logger.warning(
            f"Filtered {filtered_count} peptides with unsupported modifications"
        )
    return data


def create_de_novo_results_df(
    data_path: Path,
    filter_prefix: str | None = None,
) -> pl.DataFrame:
    """Create de novo results DataFrame from InstaNovo predictions CSV.

    Args:
        data_path: Path to InstaNovo predictions CSV file
        filter_prefix: Optional prefix to filter spectrum_id column by.
            InstaNovo spectrum_id format is "{filename}:{index}".
            If provided, only rows where spectrum_id starts with "{filter_prefix}:"
            are kept. This allows extracting predictions for a single MGF file
            from a combined predictions CSV.

    Returns:
        DataFrame with NovoBoard-compatible columns:
        - Source File: The source file name
        - Scan: The scan number of the spectrum
        - Peptide: The predicted peptide
        - ALC (%): The score of the predicted peptide
        - local confidence (%): The token-level scores
    """
    logger.info(f"Reading de novo predictions from: {data_path}")
    data = pl.read_csv(data_path)

    # Filter by spectrum_id prefix if specified (empty string disables filtering)
    if filter_prefix is not None and filter_prefix != "":
        prefix_pattern = f"{filter_prefix}:"
        initial_count = len(data)
        data = data.filter(pl.col("spectrum_id").str.starts_with(prefix_pattern))
        logger.info(
            f"Filtered to {len(data)} predictions matching prefix '{filter_prefix}' "
            f"(from {initial_count} total)"
        )

        if len(data) == 0:
            logger.warning(f"No predictions found with prefix '{filter_prefix}'")
            # Show some example spectrum_ids to help debug
            sample_data = pl.read_csv(data_path).head(5)
            if "spectrum_id" in sample_data.columns:
                examples = sample_data["spectrum_id"].to_list()
                logger.warning(f"Example spectrum_ids in file: {examples}")

    data = data.rename(COLUMN_MAPPING)
    preprocessed_data = add_score_column(data)
    preprocessed_data = map_supported_modifications(preprocessed_data)
    preprocessed_data = filter_unsupported_modifications(preprocessed_data)
    logger.info(f"Processed {len(preprocessed_data)} peptides")
    return preprocessed_data


def read_mgf(data_path: Path) -> pl.DataFrame:
    """Read labeled MGF file and extract scan numbers and peptide sequences.

    Args:
        data_path: Path to labeled MGF file

    Returns:
        DataFrame with 'Scan' and 'Peptide' columns
    """
    logger.info(f"Reading labeled MGF from: {data_path}")
    spectra = list(load_from_mgf(str(data_path)))
    data: dict[str, list] = {
        "Scan": [],
        "Peptide": [],
    }

    for i, spectrum in enumerate(spectra):
        data["Scan"].append(i)
        data["Peptide"].append(spectrum.metadata.get("peptide_sequence", ""))

    return pl.DataFrame(data)


def add_source_file_column(data: pl.DataFrame, data_path: Path) -> pl.DataFrame:
    """Add source file column to DataFrame.

    Args:
        data: Input DataFrame
        data_path: Path to use for source file name (stem is extracted)

    Returns:
        DataFrame with 'Source File' column added
    """
    exp_name = data_path.stem
    data = data.with_columns(pl.lit(exp_name).alias("Source File").cast(pl.Utf8))
    return data


def create_database_results_df(data_path: Path) -> pl.DataFrame:
    """Create database results DataFrame from labeled MGF file.

    Args:
        data_path: Path to labeled MGF file with peptide annotations

    Returns:
        DataFrame with NovoBoard-compatible columns:
        - Source File: The source file name
        - Scan: The scan number of the spectrum
        - Peptide: The annotated peptide sequence
    """
    data = read_mgf(data_path)
    data = add_source_file_column(data, data_path)
    preprocessed_data = map_supported_modifications(data)
    preprocessed_data = filter_unsupported_modifications(preprocessed_data)
    logger.info(f"Processed {len(preprocessed_data)} database entries")
    return preprocessed_data


def run_preprocessing(
    denovo_file: Path | None = None,
    denovo_output: Path | None = None,
    db_mgf_file: Path | None = None,
    db_output: Path | None = None,
    filter_prefix: str | None = None,
) -> None:
    """Run preprocessing to convert InstaNovo output to NovoBoard format.

    Args:
        denovo_file: Path to InstaNovo predictions CSV
        denovo_output: Path for de novo results output CSV
        db_mgf_file: Path to labeled MGF file (for database results)
        db_output: Path for database results output CSV
        filter_prefix: Prefix to filter spectrum_id by.
            If None and db_mgf_file is provided, defaults to the MGF file stem
            (e.g., "hepg2" for "hepg2.mgf"). This matches InstaNovo's internal
            experiment naming convention.
            Set to empty string "" to disable filtering entirely.
    """
    # Default filter_prefix to MGF file stem if not specified
    effective_filter = filter_prefix
    if effective_filter is None and db_mgf_file is not None:
        effective_filter = db_mgf_file.stem
        logger.info(
            f"Auto-detected filter prefix from MGF filename: '{effective_filter}'"
        )

    if denovo_file and denovo_output:
        de_novo_results = create_de_novo_results_df(
            denovo_file, filter_prefix=effective_filter
        )
        de_novo_results = de_novo_results.select(
            ["Source File", "Scan", "Peptide", "ALC (%)", "local confidence (%)"]
        )
        denovo_output.parent.mkdir(parents=True, exist_ok=True)
        de_novo_results.write_csv(denovo_output)
        logger.info(f"Wrote de novo results to: {denovo_output}")

    if db_mgf_file and db_output:
        database_results = create_database_results_df(db_mgf_file)
        database_results = database_results.select(["Source File", "Scan", "Peptide"])
        db_output.parent.mkdir(parents=True, exist_ok=True)
        database_results.write_csv(db_output)
        logger.info(f"Wrote database results to: {db_output}")
