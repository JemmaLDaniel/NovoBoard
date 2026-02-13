"""False Discovery Rate (FDR) estimation and validation for de novo sequencing.

Background: The Target-Decoy Approach
-------------------------------------
In proteomics, we need to estimate how many of our peptide identifications are false
positives. Traditional database search uses a target-decoy strategy: search against
both real (target) and shuffled/reversed (decoy) protein sequences. Any match to a
decoy is definitionally a false positive, allowing FDR estimation.

The Challenge for De Novo Sequencing
------------------------------------
De novo sequencing predicts peptide sequences directly from spectra without a database.
This breaks the traditional target-decoy paradigm because there's no "decoy database"
to search against. Instead, we use **spectrum-level decoys**: we corrupt the input
spectra themselves to create "impossible" spectra that should not match any real
peptide. Predictions on these corrupted spectra represent false positives.

How This Module Works
---------------------
1. **calculate_FDR()**: Combines target (real) and decoy (corrupted) spectrum
   predictions, sorts by confidence score, and estimates FDR using the classic
   formula: FDR = (decoy hits) / (target hits) at each score threshold.

2. **validate_FDR()**: When we have ground-truth labels (from database search),
   we can validate whether our estimated FDR is accurate. This is crucial for
   benchmarking - if estimated FDR says 1%, is it actually ~1% incorrect?

   We calculate "true FDR" at multiple levels:
   - Peptide-level: exact sequence match
   - Ion-level (100%): all fragment ions match
   - Ion-level (threshold): some fraction of fragment ions match

   Plotting estimated vs true FDR reveals calibration quality - a diagonal line
   means perfect calibration; curves above the diagonal mean FDR is underestimated.

Monotonic Filtering
-------------------
FDR should theoretically decrease as we raise the confidence threshold (keep only
top-scoring predictions). In practice, sampling noise can cause local increases.
Monotonic filtering removes these "bumps" by only reporting points where both
estimated and true FDR decrease, producing cleaner visualization curves.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from novoboard import config
from novoboard.accuracy import WorkerTest


@dataclass
class FDRValidationResult:
    """Container for FDR validation results comparing estimated vs true FDR.

    This class holds all data needed to assess how well target-decoy FDR estimation
    matches reality (ground truth from database search). The key comparison is:
    - estimated_fdr: What the target-decoy method predicts the FDR to be
    - true_fdr variants: What the actual error rate is based on known correct answers

    Plotting estimated vs true FDR reveals calibration quality. Points on the
    diagonal (y=x) indicate perfect calibration; points above mean FDR is
    underestimated (dangerous); points below mean FDR is overestimated (conservative).

    Attributes:
        denovo_df: Full DataFrame with de novo results, accuracy metrics, and FDR
        df: Filtered DataFrame containing only annotated target spectra (those with
            database matches for ground truth comparison)
        estimated_fdr: Estimated FDR from target-decoy competition (x-axis for plots)
        cumsum: Cumulative PSM count at each FDR threshold (for counting identifications)
        true_fdr: True FDR using exact peptide sequence matching (strictest metric)
        true_fdr_I: True FDR using 100% fragment ion matching (slightly more lenient)
        true_fdr_T: True FDR using threshold-based ion matching (most lenient)
        estimated_fdr_full: Estimated FDR on ALL target spectra (not just annotated)
        cumsum_full: Cumulative PSM count on all target spectra
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

    def get_true_fdr(self, metric: str) -> tuple[float, ...]:
        """Get the true FDR tuple for the specified metric.

        Args:
            metric: One of "peptide", "ion-100", or "ion-threshold"

        Returns:
            The corresponding true FDR tuple
        """
        if metric == "peptide":
            return self.true_fdr
        elif metric == "ion-100":
            return self.true_fdr_I
        else:  # ion-threshold
            return self.true_fdr_T


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
    """Estimate FDR using target-decoy competition.

    The Core Idea
    -------------
    We have two sets of predictions:
    - **Target**: Predictions on real, uncorrupted spectra (mix of true/false positives)
    - **Decoy**: Predictions on corrupted spectra (definitionally all false positives)

    By sorting all predictions by score and tracking the ratio of decoy to target hits,
    we estimate what fraction of target hits are likely false positives.

    The assumption: at any score threshold, the rate of false positives among targets
    equals the rate of decoy hits (since both are drawn from the "null" distribution
    of incorrect matches). This is the target-decoy competition principle.

    Algorithm
    ---------
    1. Combine target and decoy predictions
    2. Sort by confidence score (descending) - best predictions first
    3. Walk down the list, counting cumulative target and decoy hits
    4. At each position: estimated_FDR = cumsum_decoy / cumsum_target

    Args:
        target_csv: Path to target (original) de novo results
        decoy_csv: Path to decoy de novo results
        engine_score: Column name for scoring
        fdr_list: List of FDR thresholds to calculate
        selected_features: Optional set of feature IDs to consider

    Returns:
        Tuple of (combined_df, fdr_df, score_thresholds, counts)
    """
    logger.info("Loading predictions for FDR calculation...")
    logger.info(f"  Target predictions: {target_csv}")
    logger.info(f"  Decoy predictions: {decoy_csv}")
    target_psm = read_denovo(target_csv, selected_features)
    decoy_psm = read_denovo(decoy_csv, selected_features)
    logger.info(f"  Loaded {len(target_psm):,d} target PSMs")
    logger.info(f"  Loaded {len(decoy_psm):,d} decoy PSMs")
    # Step 1: Combine target and decoy predictions into a single DataFrame
    # This enables direct competition - we'll rank them together by score
    dfs = (
        pd.concat([target_psm, decoy_psm], keys=["target", "decoy"])
        .reset_index()
        .rename(columns={"level_0": "spectrum"})
    )
    dfs["is_target"] = dfs["spectrum"] == "target"

    # Step 2: Sort by score (descending) - best predictions first
    # Secondary sort by is_target ensures targets win ties (conservative FDR estimate)
    dfs.sort_values(
        by=[engine_score, "is_target"], ascending=[False, False], inplace=True
    )

    # Create a copy for FDR calculation with distinct decoy feature IDs
    dfs_fdr = dfs.copy()
    dfs_fdr.loc[~dfs_fdr["is_target"], "feature_id"] = (
        dfs_fdr.loc[~dfs_fdr["is_target"], "feature_id"] + "||decoy"
    )
    logger.info("Running target-decoy competition...")
    logger.info(f"  Combined predictions: {len(dfs_fdr):,d}")
    logger.info(f"  Target predictions: {sum(dfs_fdr['is_target']):,d}")
    logger.info(f"  Decoy predictions: {len(dfs_fdr) - sum(dfs_fdr['is_target']):,d}")

    # Step 3: Calculate cumulative FDR as we walk down the ranked list
    # At position i: FDR ≈ (# decoys seen) / (# targets seen)
    # Intuition: decoys represent "impossible" matches, so their rate estimates
    # the rate of false matches among targets
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
    tp_metric: str = "ion-threshold",
) -> FDRValidationResult:
    """Validate estimated FDR against ground truth from database search.

    Why Validation Matters
    ----------------------
    Estimated FDR from target-decoy competition is based on assumptions that may
    not hold perfectly:
    - Decoy spectra generate predictions with the same score distribution as false
      positive target predictions
    - The scoring function ranks true positives above false positives consistently

    When we have ground truth labels (from database search on the same spectra),
    we can measure "true FDR" - the actual fraction of incorrect predictions at
    each score threshold. Comparing estimated vs true FDR reveals whether our
    FDR estimation is:
    - Well-calibrated (estimated ≈ true): the diagonal line
    - Underestimated (estimated < true): dangerous, we're overconfident
    - Overestimated (estimated > true): conservative but wasteful

    True FDR Calculation
    --------------------
    We calculate correctness at multiple stringency levels:
    - **Peptide-level**: Exact amino acid sequence match (strictest)
    - **Ion-level 100%**: All theoretical fragment ions found in spectrum
    - **Ion-level T%**: At least T% of theoretical fragment ions found

    The ion-level metrics are more lenient because:
    - Minor mass errors or missing/extra peaks don't invalidate the ID
    - A partial match may still be biologically meaningful

    Args:
        target_csv: Path to target de novo results
        decoy_csv: Path to decoy de novo results
        engine_score: Column name for scoring
        db_csv: Path to database search results (ground truth)
        spectrum_file: Path to MGF spectrum file
        p_decoy: List of FDR thresholds
        T_pct: Threshold percentage for ion matching (e.g., 0.9 = 90% of ions)
        col_score: Score column name for accuracy calculation
        col_aa_score: AA score column name
        monotonic: If True, filter to monotonically decreasing FDR (default: True)
        tp_metric: True positive metric - "peptide", "ion-100", or "ion-threshold"

    Returns:
        FDRValidationResult with estimated and true FDR curves for plotting
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
    logger.info(f"Saved combined target-decoy results to: {target_decoy_csv}")

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

    # Add accuracy metrics from the accuracy file to the de novo results
    denovo_df["database_peptide"] = (
        accuracy_df["target_sequence"]
        .str.replace("C(Carbamidomethylation)", "C(+57.02)", regex=False)
        .str.replace("M(Oxidation)", "M(+15.99)", regex=False)
        .str.replace(",", "", regex=False)
    )
    denovo_df["matched_amino_acid_count"] = accuracy_df["matched_amino_acid_count"]
    denovo_df["predicted_sequence_length"] = accuracy_df["predicted_sequence_length"]
    # Derive boolean correctness flags for FDR calculation
    denovo_df["is_exact_sequence_match"] = (
        accuracy_df["matched_amino_acid_count"]
        == accuracy_df["predicted_sequence_length"]
    )
    denovo_df["target_ion_count"] = accuracy_df["target_ion_count"]
    denovo_df["matched_ion_count"] = accuracy_df["matched_ion_count"]
    denovo_df["is_all_ions_matched"] = (
        accuracy_df["matched_ion_count"] == accuracy_df["target_ion_count"]
    )
    denovo_df["is_threshold_ions_matched"] = (
        accuracy_df["matched_ion_count"] >= accuracy_df["target_ion_count"] * T_pct
    )

    # Log selected metric
    if tp_metric == "peptide":
        metric_description = "peptide-level (Novor algorithm)"
    elif tp_metric == "ion-100":
        metric_description = "ion-level 100%"
    else:  # ion-threshold
        metric_description = f"ion-level threshold ({T_pct:.0%})"
    denovo_df["tp_metric"] = tp_metric

    logger.info(f"Using true positive metric: {metric_description}")
    logger.info("Validating FDR against database ground truth...")
    has_accuracy = ~denovo_df["matched_amino_acid_count"].isna()
    mask = has_accuracy & denovo_df["is_target"]
    logger.info(f"  Total PSMs in analysis: {len(denovo_df):,d}")
    logger.info(f"  PSMs with database ground truth: {has_accuracy.sum():,d}")
    logger.info(f"  Target PSMs with ground truth (for validation): {mask.sum():,d}")

    # Calculate TRUE FDR on annotated target spectra (those with database matches)
    # Unlike estimated FDR which uses decoy counts, true FDR uses actual correctness:
    # true_FDR = (# incorrect predictions) / (# total predictions) at each threshold
    #
    # We filter to only "annotated" spectra - those that have a database match to
    # compare against. Spectra without database matches cannot contribute to true FDR
    # calculation since we don't know their ground truth.
    df = denovo_df[has_accuracy & denovo_df["is_target"]].copy()

    # Drop the first N predictions for FDR stability
    # The highest-confidence predictions have high variance in FDR estimation due to
    # small sample sizes. Dropping them produces more stable FDR curves.
    drop_count = config.FDR_DROP_COUNT
    if drop_count > 0 and len(df) > drop_count:
        df = df.iloc[drop_count:]
        logger.info(f"  Dropped first {drop_count} predictions for FDR stability")

    cumsum = range(1, len(df) + 1)

    # Peptide-level true FDR: prediction is correct only if ALL amino acids match
    cumsum_correct = np.cumsum(np.array(df["is_exact_sequence_match"].astype(int)))
    cumsum_false = cumsum - cumsum_correct
    true_fdr = cumsum_false / cumsum

    # Ion-level 100% true FDR: correct if ALL theoretical fragment ions are matched
    # This is slightly more lenient than peptide-level (allows for equivalent masses
    # like I/L substitution that produce identical fragmentation)
    cumsum_correct = np.cumsum(np.array(df["is_all_ions_matched"].astype(int)))
    cumsum_false = cumsum - cumsum_correct
    true_fdr_I = cumsum_false / cumsum

    # Ion-level threshold true FDR: correct if ≥T% of fragment ions are matched
    # Most lenient metric - accepts partial matches as "correct enough"
    cumsum_correct = np.cumsum(np.array(df["is_threshold_ions_matched"].astype(int)))
    cumsum_false = cumsum - cumsum_correct
    true_fdr_T = cumsum_false / cumsum

    # Store true FDR values as per-PSM columns
    ion_threshold_pct = int(T_pct * 100)
    df["true_fdr_peptide"] = true_fdr
    df["true_fdr_ion100"] = true_fdr_I
    df[f"true_fdr_ion{ion_threshold_pct}"] = true_fdr_T

    # Calculate q-values (running minimum FDR from lowest to highest score)
    # Q-value represents the minimum FDR at which this PSM would be accepted.
    # Unlike FDR which can fluctuate, q-values are monotonically decreasing
    # as confidence increases, making them more suitable for thresholding.
    def compute_q_values(fdr_array: np.ndarray) -> list[float]:
        """Compute q-values using Winnow's algorithm (running minimum FDR)."""
        q_values: list[float] = []
        fdr_min = float("inf")
        # Walk backwards (lowest confidence to highest)
        for current_fdr in reversed(fdr_array):
            if current_fdr > fdr_min:
                q_values.append(fdr_min)
            else:
                q_values.append(current_fdr)
                fdr_min = current_fdr
        q_values.reverse()  # Restore original order (high to low confidence)
        return q_values

    # Compute q-values for each metric (monotonic minimum FDR)
    df["true_q_value_peptide"] = compute_q_values(true_fdr)
    df["true_q_value_ion100"] = compute_q_values(true_fdr_I)
    df[f"true_q_value_ion{ion_threshold_pct}"] = compute_q_values(true_fdr_T)

    # Compute estimated q-value from estimated FDR (doesn't require ground truth)
    df["estimated_q_value"] = compute_q_values(np.array(df["estimated_fdr"]))

    logger.info(f"  Computed q-values for {len(df):,d} PSMs")

    # Save FDR and q-values to CSV
    qvalue_file = str(output_dir / f"{target_stem}_fdr_{decoy_stem}_fdr_qvalues.csv")
    qvalue_cols = [
        "Peptide",
        "estimated_fdr",
        "estimated_q_value",
        "true_fdr_peptide",
        "true_fdr_ion100",
        f"true_fdr_ion{ion_threshold_pct}",
        "true_q_value_peptide",
        "true_q_value_ion100",
        f"true_q_value_ion{ion_threshold_pct}",
        "tp_metric",
    ]
    # Only include columns that exist in df
    qvalue_cols = [c for c in qvalue_cols if c in df.columns]
    df[qvalue_cols].to_csv(qvalue_file)
    logger.info(f"Saved q-values to: {qvalue_file}")

    # Monotonic filtering: remove "bumps" where FDR increases as threshold rises
    #
    # Theoretically, FDR should monotonically decrease as we raise the score threshold
    # (keeping only higher-confidence predictions). However, finite sample sizes cause
    # local fluctuations - you might have a stretch of correct predictions followed
    # by a few incorrect ones, causing temporary FDR increases.
    #
    # For visualization, these bumps create confusing, non-monotonic curves. Monotonic
    # filtering walks backward through the data (from highest to lowest threshold) and
    # only keeps points where BOTH estimated and true FDR are at their minimum so far.
    # This produces cleaner curves that better represent the underlying trend.
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
        if reported:
            (
                estimated_fdr_out,
                cumsum_out,
                true_fdr_out,
                true_fdr_I_out,
                true_fdr_T_out,
            ) = zip(*reported)
        else:
            # No data points after filtering - return empty tuples
            estimated_fdr_out = ()
            cumsum_out = ()
            true_fdr_out = ()
            true_fdr_I_out = ()
            true_fdr_T_out = ()
    else:
        # Return all data points without filtering
        # Reverse to match monotonic branch order (low score → high score)
        # so plotting can reverse back to high score → low score consistently
        estimated_fdr_out = tuple(df["estimated_fdr"].iloc[::-1])
        cumsum_out = tuple(reversed(cumsum))
        true_fdr_out = tuple(reversed(true_fdr))
        true_fdr_I_out = tuple(reversed(true_fdr_I))
        true_fdr_T_out = tuple(reversed(true_fdr_T))

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
        if reported_full:
            estimated_fdr_full, cumsum_full_out = zip(*reported_full)
        else:
            # No data points after filtering - return empty tuples
            estimated_fdr_full = ()
            cumsum_full_out = ()
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
