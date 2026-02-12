"""Decoy FDR calculation and validation."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from novoboard.accuracy import WorkerTest


@dataclass
class FDRValidationResult:
    """Results from FDR validation.

    Attributes:
        denovo_df: DataFrame with de novo sequencing results and accuracy info
        df: Filtered DataFrame with only annotated target spectra
        estimated_fdr: Estimated FDR values
        cumsum: Cumulative count of PSMs
        true_fdr: True FDR based on peptide-level accuracy
        true_fdr_I: True FDR based on ion-level accuracy (100% match)
        true_fdr_T: True FDR based on ion-level accuracy (threshold match)
        estimated_fdr_full: Estimated FDR on all target spectra
        cumsum_full: Cumulative count on all target spectra
    """

    denovo_df: pd.DataFrame
    df: pd.DataFrame
    estimated_fdr: tuple[float, ...]
    cumsum: tuple[int, ...]
    true_fdr: tuple[float, ...]
    true_fdr_I: tuple[float, ...]
    true_fdr_T: tuple[float, ...]
    estimated_fdr_full: tuple[float, ...]
    cumsum_full: tuple[int, ...]


logger = logging.getLogger(__name__)


def read_denovo(
    denovo_csv: str,
    selected_features: set[str] | None = None,
) -> pd.DataFrame:
    """Read de novo sequencing results from CSV file.

    Args:
        denovo_csv: Path to CSV file from PEAKS de novo
        selected_features: Optional set of feature IDs to keep

    Returns:
        DataFrame with feature_id column added
    """
    denovo_psm = pd.read_csv(denovo_csv, keep_default_na=False)
    # Vectorized feature_id creation
    denovo_psm["feature_id"] = (
        denovo_psm["Source File"].str.split(".mgf").str[0]
        + "||"
        + denovo_psm["Scan"].astype(str)
    )
    if selected_features:
        # Vectorized filtering using isin
        denovo_psm = denovo_psm[denovo_psm["feature_id"].isin(selected_features)]
    return denovo_psm


def calculate_FDR(
    target_csv: str,
    decoy_csv: str,
    engine_score: str,
    fdr_list: list[float],
    selected_features: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, list[float], list[int]]:
    """Calculate FDR using target-decoy approach.

    Args:
        target_csv: Path to target (original) de novo results
        decoy_csv: Path to decoy de novo results
        engine_score: Column name for scoring
        fdr_list: List of FDR thresholds to calculate
        selected_features: Optional set of feature IDs to consider

    Returns:
        Tuple of (combined_df, fdr_df, score_thresholds, counts)
    """
    logger.info(f"target_csv = {target_csv}")
    logger.info(f"decoy_csv = {decoy_csv}")
    target_psm = read_denovo(target_csv, selected_features)
    decoy_psm = read_denovo(decoy_csv, selected_features)
    logger.info(f"len(target_psm) = {len(target_psm)}")
    logger.info(f"len(decoy_psm) = {len(decoy_psm)}")
    dfs = (
        pd.concat([target_psm, decoy_psm], keys=["target", "decoy"])
        .reset_index()
        .rename(columns={"level_0": "spectrum"})
    )
    # Vectorized is_target
    dfs["is_target"] = dfs["spectrum"] == "target"

    # target-decoy competition
    dfs.sort_values(
        by=[engine_score, "is_target"], ascending=[False, False], inplace=True
    )
    # competition on whole dataset
    dfs_fdr = dfs.copy()
    # Vectorized feature_id modification for decoys
    dfs_fdr.loc[~dfs_fdr["is_target"], "feature_id"] = (
        dfs_fdr.loc[~dfs_fdr["is_target"], "feature_id"] + "||decoy"
    )
    logger.info(f"len(dfs) = {len(dfs)}")
    logger.info(f"len(dfs_fdr) = {len(dfs_fdr)}")
    logger.info(f"sum(dfs_fdr['is_target']) = {sum(dfs_fdr['is_target'])}")

    # fdr estimation
    cumsum = range(1, len(dfs_fdr) + 1)
    cumsum_target = np.cumsum(np.array(dfs_fdr["is_target"].astype(int)))
    cumsum_decoy = cumsum - cumsum_target
    estimated_fdr = cumsum_decoy / cumsum_target
    dfs_fdr["estimated_fdr"] = estimated_fdr

    score_list: list[float] = []
    count_list: list[int] = []
    for fdr in fdr_list:
        fdr_index = np.flatnonzero(estimated_fdr <= fdr)
        fdr_index = fdr_index[-1] if len(fdr_index) > 0 else 0
        score_list.append(dfs_fdr.iloc[fdr_index][engine_score])
        count_list.append(fdr_index + 1)

    return dfs, dfs_fdr, score_list, count_list


def validate_FDR(
    target_csv: str,
    decoy_csv: str,
    engine_score: str,
    db_csv: str,
    spectrum_file: str,
    p_decoy: list[float],
    T_pct: float,
    col_score: str,
    col_aa_score: str,
    monotonic: bool = True,
) -> FDRValidationResult:
    """Validate FDR estimation against known database matches.

    Args:
        target_csv: Path to target de novo results
        decoy_csv: Path to decoy de novo results
        engine_score: Column name for scoring
        db_csv: Path to database search results (ground truth)
        spectrum_file: Path to MGF spectrum file
        p_decoy: List of FDR thresholds
        T_pct: Threshold percentage for ion matching
        col_score: Score column name for accuracy calculation
        col_aa_score: AA score column name
        monotonic: If True, filter to monotonically decreasing FDR (default: True)

    Returns:
        Dictionary with validation results including DataFrames and FDR curves
    """
    db_psm = pd.read_csv(db_csv, keep_default_na=False)
    # Vectorized feature_id creation
    db_psm["feature_id"] = (
        db_psm["Source File"].str.split(".mgf").str[0]
        + "||"
        + db_psm["Scan"].astype(str)
    )
    selected_features = None  # set(db_psm['feature_id'])

    dfs, dfs_fdr, score_list, count_list = calculate_FDR(
        target_csv, decoy_csv, engine_score, p_decoy, selected_features
    )

    # Build clean output filenames
    target_path = Path(target_csv)
    decoy_path = Path(decoy_csv)
    target_stem = target_path.stem
    decoy_stem = decoy_path.stem
    output_dir = target_path.parent

    # Save combined target-decoy results with estimated FDR
    target_decoy_csv = str(output_dir / f"{target_stem}_fdr_{decoy_stem}.csv")
    dfs_fdr.to_csv(target_decoy_csv, index=False)
    logger.info(f"Saved FDR results to: {target_decoy_csv}")

    # Accuracy file for validation
    accuracy_file = str(output_dir / f"{target_stem}_fdr_{decoy_stem}_accuracy.csv")
    if not Path(accuracy_file).is_file():
        worker_test = WorkerTest(
            db_csv, target_decoy_csv, spectrum_file, col_score, col_aa_score
        )
        worker_test.test_accuracy()

    denovo_df = dfs_fdr.copy()
    denovo_df = denovo_df.set_index("feature_id")
    accuracy_df = pd.read_csv(accuracy_file, delimiter="\t", index_col="feature_id")

    # Vectorized string operations
    denovo_df["db_peptide"] = (
        accuracy_df["target_sequence"]
        .str.replace("C(Carbamidomethylation)", "C(+57.02)", regex=False)
        .str.replace("M(Oxidation)", "M(+15.99)", regex=False)
        .str.replace(",", "", regex=False)
    )
    denovo_df["recall_AA"] = accuracy_df["recall_AA"]
    denovo_df["predicted_len"] = accuracy_df["predicted_len"]
    # Vectorized comparisons
    denovo_df["recall_peptide"] = (
        accuracy_df["recall_AA"] == accuracy_df["predicted_len"]
    )
    denovo_df["target_ion"] = accuracy_df["target_ion"]
    denovo_df["matched_ion"] = accuracy_df["matched_ion"]
    denovo_df["recall_peptide_I"] = (
        accuracy_df["matched_ion"] == accuracy_df["target_ion"]
    )
    denovo_df["recall_peptide_T"] = (
        accuracy_df["matched_ion"] >= accuracy_df["target_ion"] * T_pct
    )
    logger.info(f"len(denovo_df) = {len(denovo_df)}")
    logger.info(f"  with recall_AA = {len(denovo_df[~denovo_df['recall_AA'].isna()])}")
    # Fix the boolean indexing warning
    mask = ~denovo_df["recall_AA"].isna() & denovo_df["is_target"]
    logger.info(f"    is_target = {mask.sum()}")

    # calculate true FDR on annotated target spectra
    df = denovo_df[~denovo_df["recall_AA"].isna() & denovo_df["is_target"]].copy()
    cumsum = range(1, len(df) + 1)
    cumsum_correct = np.cumsum(np.array(df["recall_peptide"].astype(int)))
    cumsum_false = cumsum - cumsum_correct
    true_fdr = cumsum_false / cumsum
    cumsum_correct = np.cumsum(np.array(df["recall_peptide_I"].astype(int)))
    cumsum_false = cumsum - cumsum_correct
    true_fdr_I = cumsum_false / cumsum
    cumsum_correct = np.cumsum(np.array(df["recall_peptide_T"].astype(int)))
    cumsum_false = cumsum - cumsum_correct
    true_fdr_T = cumsum_false / cumsum
    # Optionally filter to monotonically decreasing FDR to avoid bumps
    if monotonic:
        min_est, min_true = 1.0, 1.0
        reported: list[tuple[float, int, float, float, float]] = []
        for x, y, z, v, w in list(
            zip(df["estimated_fdr"], cumsum, true_fdr, true_fdr_I, true_fdr_T)
        )[::-1]:
            if x <= min_est and w <= min_true:
                min_est = x
                min_true = w
                reported.append((x, y, z, v, w))
        estimated_fdr_out, cumsum_out, true_fdr_out, true_fdr_I_out, true_fdr_T_out = (
            zip(*reported)
        )
    else:
        # Return all data points without filtering
        estimated_fdr_out = tuple(df["estimated_fdr"])
        cumsum_out = tuple(cumsum)
        true_fdr_out = tuple(true_fdr)
        true_fdr_I_out = tuple(true_fdr_I)
        true_fdr_T_out = tuple(true_fdr_T)

    # calculate estimated FDR and #PSMs on all target spectra
    df_target = denovo_df[denovo_df["is_target"]].copy()
    cumsum_full = range(1, len(df_target) + 1)
    if monotonic:
        # Only report entries with decreasing fdr from the bottom to avoid bumps
        min_est = 1.0
        reported_full: list[tuple[float, int]] = []
        for x, y in list(zip(df_target["estimated_fdr"], cumsum_full))[::-1]:
            if x <= min_est:
                min_est = x
                reported_full.append((x, y))
        estimated_fdr_full, cumsum_full_out = zip(*reported_full)
    else:
        estimated_fdr_full = tuple(df_target["estimated_fdr"])
        cumsum_full_out = tuple(cumsum_full)

    return FDRValidationResult(
        denovo_df=denovo_df,
        df=df,
        estimated_fdr=estimated_fdr_out,
        cumsum=cumsum_out,
        true_fdr=true_fdr_out,
        true_fdr_I=true_fdr_I_out,
        true_fdr_T=true_fdr_T_out,
        estimated_fdr_full=estimated_fdr_full,
        cumsum_full=cumsum_full_out,
    )
