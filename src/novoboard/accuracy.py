"""Calculate accuracy metrics for de novo peptide sequencing predictions.

Overview
--------
When evaluating de novo peptide sequencing, we need ground truth - known correct
peptide sequences for each spectrum. This typically comes from database search on
the same spectra, which has high confidence when matched to a protein database.

This module compares de novo predictions against database search results at multiple
granularity levels, from lenient to strict.

Accuracy Metrics: A Hierarchy of Stringency
-------------------------------------------

**Fragment Ion Level** (Most lenient)
    Compare theoretical fragment ions of predicted vs. target peptide against the
    observed spectrum. A prediction that produces similar fragmentation patterns
    to the target may be "good enough" even if not sequence-identical.

    Why this matters: Some substitutions (e.g., I↔L with identical mass) cannot be
    distinguished by MS/MS. Ion-level matching catches these.

**Amino Acid Level** (Intermediate)
    Count how many individual amino acids are correctly identified, allowing for
    position shifts. Uses cumulative mass alignment (Novor-style matching).

    Why this matters: Even partial correctness is valuable. A prediction of
    "PEPTIDE" vs. target "PEPTYDE" has 6/7 AA correct - still useful information.

**Peptide Level** (Strictest)
    Binary: is the entire predicted sequence exactly correct (all AAs match)?

    Why this matters: For some applications (e.g., neoantigen discovery), only
    exact sequences are useful. This metric reflects end-to-end accuracy.

The Novor Matching Algorithm
----------------------------
Standard string comparison doesn't work well for peptide matching because:
1. Mass spectrometry has finite mass accuracy (can't distinguish I/L)
2. Predictions may have insertions/deletions that shift positions
3. We care about MASS correctness, not character-by-character matching

The Novor algorithm (Ma, 2015) compares cumulative masses instead:
1. Convert each sequence to cumulative mass arrays (prefix masses)
2. Walk through both arrays, matching positions where cumulative masses align
3. At each aligned position, check if the individual AA masses also match

This elegantly handles mass-equivalent substitutions and small position shifts.

Output Files
------------
- *_accuracy.csv: Per-spectrum accuracy metrics (recall, ion matching, etc.)
- *_denovo_only.csv: Predictions without database matches (novel sequences)
- *_scan2fea.csv: Scan-to-feature mapping (for LC-MS features spanning scans)
- *_multifea.csv: Features assigned to multiple scans (chimeric spectra)

References
----------
- Ma, B. (2015). "Novor: Real-time peptide de novo sequencing software."
"""

from __future__ import annotations

import logging
import re
import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from novoboard import config

# Precompiled regex patterns for performance
_SCAN_SPLIT_PATTERN = re.compile(r";|\r|\n")
_CSV_SPLIT_PATTERN = re.compile(r",|\r|\n")
_MGF_FIELD_PATTERN = re.compile(r"[=\r\n]")
_PEAK_SPLIT_PATTERN = re.compile(r" |\t|\r|\n")


@dataclass
class PredictedPeptide:
    """A predicted peptide from de novo sequencing.

    Attributes:
        feature_id: Unique identifier for the feature (source_file||scan)
        sequence: List of predicted sequences (multiple predictions possible)
        score: List of scores for each prediction
        aa_score: List of amino acid scores for each prediction
        feature_area: Peak area (intensity) of the feature
        scan_list_middle: Middle scan of the feature
        scan_list_original: Original scan list
    """

    feature_id: str
    sequence: list[list[str]] = field(default_factory=list)
    score: list[float] = field(default_factory=list)
    aa_score: list[str] = field(default_factory=list)
    feature_area: float = 0.0
    scan_list_middle: str = ""
    scan_list_original: str = ""


logger = logging.getLogger(__name__)

col_precursor_mz = "m/z"
col_precursor_charge = "z"
col_rt_mean = "RT"
col_raw_sequence = "Peptide"
col_source_file = "Source File"
col_scan_list = "Scan"


def parse_raw_sequence(raw_sequence: str) -> tuple[bool, list[str]]:
    """Parse a peptide sequence string into tokenized amino acids with modifications.

    Peptide Sequence Representation
    --------------------------------
    Mass spectrometry peptide sequences are typically written with modifications
    in parentheses, e.g., "PEPTC(+57.02)IDE" where:
    - Standard amino acids are single letters (P, E, T, I, D)
    - Modified amino acids have mass shifts in parentheses after the residue

    Common modifications:
    - C(+57.02): Carbamidomethylation of cysteine (alkylation artifact)
    - M(+15.99): Oxidation of methionine (common oxidation)
    - N/Q(+0.98): Deamidation (asparagine/glutamine → aspartic/glutamic acid)
    - S/T/Y(+79.97): Phosphorylation (post-translational modification)

    Why Tokenization?
    -----------------
    Downstream processing (mass calculation, vocabulary lookup) needs discrete
    amino acid units. "PEPTC(+57.02)IDE" becomes ["P", "E", "P", "T",
    "C(Carbamidomethylation)", "I", "D", "E"] - a list of 8 tokens.

    Args:
        raw_sequence: Peptide string with modifications (e.g., "PEPTC(+57.02)DE")

    Returns:
        Tuple of (success, peptide_tokens):
        - success: True if all tokens are in the known vocabulary
        - peptide_tokens: List of amino acid strings including modifications
    """
    raw_sequence_len = len(raw_sequence)
    peptide: list[str] = []
    index = 0
    while index < raw_sequence_len:
        if raw_sequence[index] == "(":
            if peptide[-1] == "C" and raw_sequence[index : index + 8] == "(+57.02)":
                peptide[-1] = "C(Carbamidomethylation)"
                index += 8
            elif peptide[-1] == "M" and raw_sequence[index : index + 8] == "(+15.99)":
                peptide[-1] = "M(Oxidation)"
                index += 8
            elif peptide[-1] == "N" and raw_sequence[index : index + 7] == "(+0.98)":
                peptide[-1] = "N(Deamidation)"
                index += 7
            elif peptide[-1] == "Q" and raw_sequence[index : index + 7] == "(+0.98)":
                peptide[-1] = "Q(Deamidation)"
                index += 7
            elif peptide[-1] == "S" and raw_sequence[index : index + 8] == "(+79.97)":
                peptide[-1] = "S(Phosphorylation)"
                index += 8
            elif peptide[-1] == "T" and raw_sequence[index : index + 8] == "(+79.97)":
                peptide[-1] = "T(Phosphorylation)"
                index += 8
            elif peptide[-1] == "Y" and raw_sequence[index : index + 8] == "(+79.97)":
                peptide[-1] = "Y(Phosphorylation)"
                index += 8
            else:  # unknown modification
                return False, peptide
        else:
            peptide.append(raw_sequence[index])
            index += 1

    for aa in peptide:
        if aa not in config.vocab:
            return False, peptide
    return True, peptide


class WorkerTest:
    """Calculate accuracy metrics comparing de novo predictions to database search.

    The Core Task
    -------------
    Given two sets of peptide identifications for the same spectra:
    - Target (ground truth): from database search, assumed correct
    - Predicted: from de novo sequencing, being evaluated

    Calculate how well predictions match targets at multiple stringency levels.

    Why Multiple Metrics?
    ---------------------
    A single accuracy number doesn't capture the full picture:

    - A model predicting "PEPTIDE" when target is "PEPTYDE" gets 0% peptide-level
      accuracy but 86% AA-level accuracy. For some applications, that's great!

    - A model that consistently predicts correct peptides but with I/L swaps
      (indistinguishable by mass) should get credit via ion-level metrics.

    - Different downstream applications have different requirements:
      * Database validation: exact match required
      * Novel peptide discovery: partial matches useful
      * Immunopeptidomics: specific epitope regions matter most

    Output Files
    ------------
    - accuracy_file: Main output with per-spectrum metrics
    - denovo_only_file: Predictions without corresponding database matches
      (potentially novel peptides not in the search database)
    - scan2fea_file: Mapping of scans to features (for LC-MS feature grouping)
    - multifea_file: Features spanning multiple scans (chimeric spectra flags)
    """

    def __init__(
        self,
        target_file: str,
        predicted_file: str,
        spectrum_file: str,
        col_score: str,
        col_aa_score: str,
    ) -> None:
        """Initialize WorkerTest with file paths and column names.

        Args:
            target_file: Path to CSV file containing target (database) peptides
            predicted_file: Path to CSV file containing predicted peptides
            spectrum_file: Path to MGF file containing spectra
            col_score: Column name for peptide score
            col_aa_score: Column name for amino acid scores
        """
        logger.info("=" * 80)
        logger.info("Initializing accuracy calculation...")

        self.MZ_MAX: float = config.MZ_MAX

        self.target_file = target_file
        self.predicted_file = predicted_file
        self.spectrum_file = spectrum_file

        # Build clean output filenames from predicted_file
        input_path = Path(predicted_file)
        stem = input_path.stem
        parent = input_path.parent
        self.accuracy_file = str(parent / f"{stem}_accuracy.csv")
        self.denovo_only_file = str(parent / f"{stem}_denovo_only.csv")
        self.scan2fea_file = str(parent / f"{stem}_scan2fea.csv")
        self.multifea_file = str(parent / f"{stem}_multifea.csv")
        self.col_score = col_score
        self.col_aa_score = col_aa_score
        logger.info("Input files:")
        logger.info(f"  Database (target) peptides: {self.target_file}")
        logger.info(f"  De novo predictions: {self.predicted_file}")
        logger.info(f"  Spectra (MGF): {self.spectrum_file}")
        logger.info("Output files will be written to:")
        logger.info(f"  Accuracy results: {self.accuracy_file}")
        logger.info(f"  De novo only (no DB match): {self.denovo_only_file}")
        logger.info(f"  Scan-to-feature mapping: {self.scan2fea_file}")
        logger.info(f"  Multi-feature entries: {self.multifea_file}")

        self.target_dict: dict[str, list[str]] = {}
        self.predicted_list: list[dict] = []
        self.spectrum_dict: dict[str, list[tuple[float, float]]] = {}

    def _filter_targets_by_db(
        self,
        db_peptide_list: list[list[str]] | None,
    ) -> dict[str, list[str]]:
        """Filter targets to only include those in db_peptide_list.

        Args:
            db_peptide_list: List of peptides to filter by, or None for all targets

        Returns:
            Filtered dictionary of feature_id -> peptide sequence
        """
        if db_peptide_list is None:
            return self.target_dict

        target_dict_db: dict[str, list[str]] = {}
        for feature_id, target in self.target_dict.items():
            # Remove extension 'mod' from variable modifications for comparison
            target_simplified = target[:]
            target_simplified = [
                "M" if x == "M(Oxidation)" else x for x in target_simplified
            ]
            target_simplified = [
                "N" if x == "N(Deamidation)" else x for x in target_simplified
            ]
            target_simplified = [
                "Q" if x == "Q(Deamidation)" else x for x in target_simplified
            ]
            if target_simplified in db_peptide_list:
                target_dict_db[feature_id] = target
            else:
                logger.warning(
                    f"Target peptide not in filter list (skipping): "
                    f"{''.join(target_simplified)}"
                )
        return target_dict_db

    def _filter_targets_by_mass(
        self,
        target_dict: dict[str, list[str]],
    ) -> dict[str, list[str]]:
        """Filter targets to only include those with precursor mass <= MZ_MAX.

        Args:
            target_dict: Dictionary of feature_id -> peptide sequence

        Returns:
            Filtered dictionary
        """
        return {
            feature_id: peptide
            for feature_id, peptide in target_dict.items()
            if self._compute_peptide_mass(peptide) <= self.MZ_MAX
        }

    def _find_best_prediction(
        self,
        target: list[str],
        predicted: dict,
    ) -> tuple[int, str, list[str], float, str]:
        """Find the best matching prediction from a list of predictions.

        Args:
            target: Target peptide sequence
            predicted: Dictionary with 'sequence', 'score', 'aa_score' lists

        Returns:
            Tuple of (recall_AA, aa_match, best_sequence, best_score, best_aa_score)
        """
        best_recall_AA = -1
        best_aa_match = ""
        best_predicted_sequence = predicted["sequence"][0]
        best_predicted_score = predicted["score"][0]
        best_predicted_aa_score = predicted["aa_score"][0]

        for pred_seq, pred_score, pred_aa_score in zip(
            predicted["sequence"], predicted["score"], predicted["aa_score"]
        ):
            predicted_AA_id = [config.vocab[x] for x in pred_seq]
            target_AA_id = [config.vocab[x] for x in target]
            recall_AA, aa_match = self._match_AA_novor(target_AA_id, predicted_AA_id)

            if recall_AA > best_recall_AA or (
                recall_AA == best_recall_AA and pred_score > best_predicted_score
            ):
                best_recall_AA = recall_AA
                best_aa_match = aa_match
                best_predicted_sequence = pred_seq[:]
                best_predicted_score = pred_score
                best_predicted_aa_score = pred_aa_score

        return (
            best_recall_AA,
            best_aa_match,
            best_predicted_sequence,
            best_predicted_score,
            best_predicted_aa_score,
        )

    def _write_scan_files(
        self,
        scan_dict: dict[str, dict],
    ) -> None:
        """Write scan2fea and multifea output files.

        Args:
            scan_dict: Dictionary mapping scan_id -> {feature_count, feature_list}
        """
        # Build multifea_dict from scan_dict
        multifea_dict: dict[str, list[str]] = {}
        for scan_id, value in scan_dict.items():
            feature_count = value["feature_count"]
            feature_list = value["feature_list"]
            if feature_count > 1:
                for feature_id in feature_list:
                    if feature_id in multifea_dict:
                        multifea_dict[feature_id].append(f"{scan_id}:{feature_count}")
                    else:
                        multifea_dict[feature_id] = [f"{scan_id}:{feature_count}"]

        # Write scan2fea file
        with open(self.scan2fea_file, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                ["scan_id", "associated_feature_count", "associated_feature_ids"]
            )
            for scan_id, value in scan_dict.items():
                writer.writerow(
                    [scan_id, value["feature_count"], ";".join(value["feature_list"])]
                )

        # Write multifea file
        with open(self.multifea_file, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["feature_id", "shared_scans_with_counts"])
            for feature_id, scan_list in multifea_dict.items():
                writer.writerow([feature_id, ";".join(scan_list)])

    def _log_metrics(
        self,
        target_count_total: int,
        target_len_total: int,
        target_count_db: int,
        target_len_db: int,
        target_count_db_mass: int,
        target_len_db_mass: int,
        predicted_count_mass: int,
        predicted_count_mass_db: int,
        predicted_len_mass_db: int,
        predicted_only: int,
        recall_AA_total: float,
        recall_peptide_total: float,
        target_ion_total: float,
        matched_ion_total: float,
        predicted_ion_total: float,
        recall_all_peptide_ions_total: float,
    ) -> None:
        """Log all accuracy metrics."""
        logger.info("=" * 80)
        logger.info("ACCURACY RESULTS SUMMARY")
        logger.info("=" * 80)

        logger.info("Dataset statistics:")
        logger.info(f"  Total database peptides: {target_count_total:,d}")
        logger.info(f"  Total amino acids in database: {target_len_total:,d}")
        logger.info(f"  After filtering by peptide list: {target_count_db:,d} peptides")
        logger.info(
            f"  After filtering by mass (<{config.MZ_MAX} Da): "
            f"{target_count_db_mass:,d} peptides"
        )

        logger.info("Prediction statistics:")
        logger.info(f"  Predictions within mass range: {predicted_count_mass:,d}")
        logger.info(f"  Predictions with database match: {predicted_count_mass_db:,d}")
        logger.info(
            f"  Predictions without database match (novel): {predicted_only:,d}"
        )

        logger.info("-" * 40)
        logger.info("AMINO ACID-LEVEL ACCURACY:")
        logger.info(
            f"  Recall (matched AAs / target AAs): "
            f"{recall_AA_total / target_len_db_mass:.2%}"
        )
        logger.info(
            f"  Precision (matched AAs / predicted AAs): "
            f"{recall_AA_total / predicted_len_mass_db:.2%}"
        )

        logger.info("-" * 40)
        logger.info("PEPTIDE-LEVEL ACCURACY:")
        logger.info(
            f"  Recall (exact matches / target peptides): "
            f"{recall_peptide_total / target_count_db_mass:.2%}"
        )
        logger.info(
            f"  Precision (exact matches / predictions): "
            f"{recall_peptide_total / predicted_count_mass_db:.2%}"
        )

        logger.info("-" * 40)
        logger.info("FRAGMENT ION-LEVEL ACCURACY:")
        logger.info(
            f"  Ion recall (matched ions / target ions): "
            f"{matched_ion_total / target_ion_total:.2%}"
        )
        logger.info(
            f"  Ion precision (matched ions / predicted ions): "
            f"{matched_ion_total / predicted_ion_total:.2%}"
        )
        logger.info(
            f"  Perfect ion matches (100% ions correct): "
            f"{recall_all_peptide_ions_total / target_count_db_mass:.2%}"
        )
        logger.info(
            f"  Total target ions: {target_ion_total:,.0f}, "
            f"matched: {matched_ion_total:,.0f}"
        )
        logger.info("=" * 80)

    def test_accuracy(self, db_peptide_list: list[list[str]] | None = None) -> None:
        """Calculate and report all accuracy metrics, writing results to files.

        Main Workflow
        -------------
        1. Load target (database) peptides → ground truth
        2. Load predicted (de novo) peptides → what we're evaluating
        3. Load spectra → for ion-level matching
        4. For each spectrum with both target and prediction:
           - Calculate AA-level accuracy (Novor matching)
           - Calculate ion-level accuracy (fragment overlap)
           - Record to accuracy file
        5. For predictions without targets:
           - These are "de novo only" - potentially novel peptides
           - Record to denovo_only file for manual inspection
        6. Log aggregate statistics

        Optional Filtering
        ------------------
        If db_peptide_list is provided, only targets matching this list are
        considered. This allows focused evaluation on specific peptide sets
        (e.g., only tryptic peptides, only a specific protein's peptides).

        Args:
            db_peptide_list: Optional whitelist of peptides to include
        """
        logger.info("=" * 80)
        logger.info("Starting accuracy calculation...")

        # write the accuracy of predicted peptides
        accuracy_handle = open(self.accuracy_file, "w")
        header_list = [
            "feature_id",
            "feature_area",
            "target_sequence",
            "predicted_sequence",
            "predicted_score",
            "predicted_position_scores",
            "matched_amino_acid_count",
            "amino_acid_match_pattern",
            "predicted_sequence_length",
            "target_sequence_length",
            "target_ion_count",
            "matched_ion_count",
            "unmatched_target_ions_mz",
            "scan_list_middle",
            "scan_list_original",
        ]
        header_row = "\t".join(header_list)
        print(header_row, file=accuracy_handle, end="\n")

        # write denovo_only peptides
        denovo_only_handle = open(self.denovo_only_file, "w")
        header_list = [
            "feature_id",
            "feature_area",
            "predicted_sequence",
            "predicted_score",
            "max_predicted_score",
            "scan_list_middle",
            "scan_list_original",
            "is_target",
        ]
        header_row = "\t".join(header_list)
        print(header_row, file=denovo_only_handle, end="\n")

        self._get_target()
        target_count_total = len(self.target_dict)
        target_len_total = sum(len(x) for x in self.target_dict.values())

        # Filter targets by db_peptide_list (if provided)
        target_dict_db = self._filter_targets_by_db(db_peptide_list)
        target_count_db = len(target_dict_db)
        target_len_db = sum(len(x) for x in target_dict_db.values())

        # Filter targets by precursor mass
        target_dict_db_mass = self._filter_targets_by_mass(target_dict_db)
        target_count_db_mass = len(target_dict_db_mass)
        target_len_db_mass = sum(len(x) for x in target_dict_db_mass.values())

        # read predicted peptides from deepnovo or peaks
        self._get_predicted_peaks_11()

        # note that the prediction has already skipped precursor_mass > MZ_MAX
        # we also skip predicted peptides whose feature_id's are not in target_dict_db_mass
        predicted_count_mass = len(self.predicted_list)
        predicted_count_mass_db = 0
        predicted_len_mass_db = 0
        predicted_only = 0
        # the recall is calculated on remaining peptides
        recall_AA_total = 0.0
        recall_peptide_total = 0.0
        # fragment ion accuracy
        target_ion_total = 0.0
        predicted_ion_total = 0.0
        matched_ion_total = 0.0
        recall_all_peptide_ions_total = 0.0

        # record scan with multiple features
        scan_dict: dict[str, dict] = {}

        # read spectra to calculate fragment ion accuracy
        self._get_spectra()

        id_set: set[str] = set()
        for index, predicted in enumerate(self.predicted_list):
            feature_id = predicted["feature_id"]
            if feature_id in id_set:
                continue
            else:
                id_set.add(feature_id)

            feature_area = str(predicted["feature_area"])
            feature_scan_list_middle = predicted["scan_list_middle"]
            feature_scan_list_original = predicted["scan_list_original"]
            if feature_scan_list_original:
                for scan in _SCAN_SPLIT_PATTERN.split(feature_scan_list_original):
                    if scan in scan_dict:
                        scan_dict[scan]["feature_count"] += 1
                        scan_dict[scan]["feature_list"].append(feature_id)
                    else:
                        scan_dict[scan] = {}
                        scan_dict[scan]["feature_count"] = 1
                        scan_dict[scan]["feature_list"] = [feature_id]

            if feature_id in target_dict_db_mass:
                predicted_count_mass_db += 1

                target = target_dict_db_mass[feature_id]
                target_len = len(target)

                # Find the best matching prediction
                (
                    recall_AA,
                    aa_match,
                    predicted_sequence,
                    predicted_score,
                    predicted_aa_score,
                ) = self._find_best_prediction(target, predicted)

                recall_AA_total += recall_AA
                if recall_AA == target_len:
                    recall_peptide_total += 1
                predicted_len = len(predicted_sequence)
                predicted_len_mass_db += predicted_len

                if feature_id in self.spectrum_dict:
                    target_ion, predicted_ion, matched_ion, unmatched_ion_list = (
                        self._match_ion(
                            target, predicted_sequence, self.spectrum_dict[feature_id]
                        )
                    )
                else:
                    target_ion, predicted_ion, matched_ion, unmatched_ion_list = (
                        0,
                        0,
                        0,
                        "",
                    )
                target_ion_total += target_ion
                predicted_ion_total += predicted_ion
                matched_ion_total += matched_ion
                recall_all_peptide_ions = matched_ion == target_ion
                recall_all_peptide_ions_total += recall_all_peptide_ions

                # convert to string format to print out
                target_sequence = ",".join(target)
                predicted_sequence_str = ",".join(predicted_sequence)
                predicted_score_str = f"{predicted_score:.2f}"
                recall_AA_str = f"{recall_AA:d}"
                predicted_len_str = f"{predicted_len:d}"
                target_len_str = f"{target_len:d}"
                target_ion_str = f"{target_ion:d}"
                matched_ion_str = f"{matched_ion:d}"
                print_list = [
                    feature_id,
                    feature_area,
                    target_sequence,
                    predicted_sequence_str,
                    predicted_score_str,
                    predicted_aa_score,
                    recall_AA_str,
                    aa_match,
                    predicted_len_str,
                    target_len_str,
                    target_ion_str,
                    matched_ion_str,
                    unmatched_ion_list,
                    feature_scan_list_middle,
                    feature_scan_list_original,
                ]
                print_row = "\t".join(print_list)
                print(print_row, file=accuracy_handle, end="\n")
            else:
                predicted_only += 1
                predicted_sequence_out = ";".join(
                    [",".join(x) for x in predicted["sequence"]]
                )
                predicted_score_out = ";".join([f"{x:.2f}" for x in predicted["score"]])
                if predicted["score"]:
                    predicted_score_max = f'{np.max(predicted["score"]):.2f}'
                else:
                    predicted_score_max = ""
                print_list = [
                    feature_id,
                    feature_area,
                    predicted_sequence_out,
                    predicted_score_out,
                    predicted_score_max,
                    feature_scan_list_middle,
                    feature_scan_list_original,
                    str(predicted.get("is_target", "")),
                ]
                print_row = "\t".join(print_list)
                print(print_row, file=denovo_only_handle, end="\n")

        accuracy_handle.close()
        denovo_only_handle.close()

        # Write scan mapping files
        self._write_scan_files(scan_dict)

        # Log accuracy metrics
        self._log_metrics(
            target_count_total=target_count_total,
            target_len_total=target_len_total,
            target_count_db=target_count_db,
            target_len_db=target_len_db,
            target_count_db_mass=target_count_db_mass,
            target_len_db_mass=target_len_db_mass,
            predicted_count_mass=predicted_count_mass,
            predicted_count_mass_db=predicted_count_mass_db,
            predicted_len_mass_db=predicted_len_mass_db,
            predicted_only=predicted_only,
            recall_AA_total=recall_AA_total,
            recall_peptide_total=recall_peptide_total,
            target_ion_total=target_ion_total,
            matched_ion_total=matched_ion_total,
            predicted_ion_total=predicted_ion_total,
            recall_all_peptide_ions_total=recall_all_peptide_ions_total,
        )

    def _compute_peptide_mass(self, peptide: list[str]) -> float:
        """Compute the monoisotopic mass of a peptide.

        Args:
            peptide: List of amino acid strings (including modifications)

        Returns:
            Peptide mass including N/C-terminal groups
        """
        peptide_mass = (
            config.mass_N_terminus
            + sum(config.mass_AA[aa] for aa in peptide)
            + config.mass_C_terminus
        )

        return peptide_mass

    def _get_predicted_peaks_11(self) -> None:
        """Read predicted peptides from PEAKS output CSV file."""
        logger.info("Loading de novo predictions...")

        predicted_list: list[dict] = []
        with open(self.predicted_file, "r") as handle:
            csv_reader = csv.DictReader(handle)
            for row in csv_reader:
                predicted: dict = {}
                predicted["feature_id"] = (
                    row[col_source_file].split(".mgf")[0] + "||" + row[col_scan_list]
                )
                raw_sequence = row["Peptide"]
                assert raw_sequence, "Error: wrong format."
                okay, predicted["sequence"] = parse_raw_sequence(raw_sequence)
                if not okay:
                    # skip unknown mod
                    continue
                # skip peptides with precursor_mass > MZ_MAX
                if self._compute_peptide_mass(predicted["sequence"]) > self.MZ_MAX:
                    continue
                predicted["feature_area"] = 0
                predicted["scan_list_middle"] = ""
                predicted["scan_list_original"] = ""
                predicted["sequence"] = [predicted["sequence"]]
                predicted["score"] = [float(row[self.col_score])]
                predicted["aa_score"] = [row[self.col_aa_score]]
                # Read is_target if present (from FDR command's combined target-decoy file)
                predicted["is_target"] = row.get("is_target", "")
                predicted_list.append(predicted)

        self.predicted_list = predicted_list

    def _get_target(self) -> None:
        """Read target peptides from database search CSV file."""
        logger.info("Loading database (target) peptides...")

        target_dict: dict[str, list[str]] = {}
        with open(self.target_file, "r") as handle:
            header_line = handle.readline()
            header = [x.strip('"') for x in header_line.strip().split(",")]
            logger.debug(f"Header: {header}")
            raw_sequence_index = header.index(col_raw_sequence)
            source_file_index = header.index(col_source_file)
            scan_index = header.index(col_scan_list)

            for line in handle:
                line_parts = [x.strip('"') for x in _CSV_SPLIT_PATTERN.split(line)]
                # Strip .mgf extension from source file for consistent feature_id format
                source_file = line_parts[source_file_index].split(".mgf")[0]
                feature_id = source_file + "||" + line_parts[scan_index]
                raw_sequence = line_parts[raw_sequence_index]
                assert raw_sequence, "Error: wrong target format."
                okay, peptide = parse_raw_sequence(raw_sequence)
                if not okay:
                    # skip unknown mod
                    continue
                target_dict[feature_id] = peptide
        self.target_dict = target_dict

    def _get_spectra(self) -> None:
        """Read spectra from MGF file.

        Feature ID Construction
        -----------------------
        feature_id = "{source_file}||{scan}" where:
        - source_file: MGF filename stem (e.g., "gluc" from "gluc.mgf")
        - scan: 0-based spectrum index in the file (NOT the SCANS= field value)

        This matches InstaNovo's output format where:
        - experiment_name = MGF filename stem
        - scan_number = 0-based position of spectrum in the MGF file
        """
        logger.info("Loading spectra from MGF file...")

        # Use MGF filename stem as source_file (matches InstaNovo's experiment_name)
        source_file = Path(self.spectrum_file).stem
        spectrum_index = 0

        with open(self.spectrum_file, "r") as f_in:
            while True:
                line = f_in.readline()
                if not line:  # end of file
                    break
                if line == "\n":  # empty line
                    continue
                peak_list: list[tuple[float, float]] = []
                # Use 0-based spectrum index as scan (matches InstaNovo's scan_number)
                # Ignore SCANS= field - it contains original scan numbers that don't match
                scan = str(spectrum_index)
                while "END IONS" not in line:
                    # parse header lines - skip them, we only need peak data
                    if "BEGIN IONS" in line or "=" in line:
                        line = f_in.readline()
                        continue
                    # parse ions
                    mz, intensity = _PEAK_SPLIT_PATTERN.split(line)[:2]
                    peak_list.append((float(mz), float(intensity)))
                    line = f_in.readline()
                feature_id = f"{source_file}||{scan}"
                self.spectrum_dict[feature_id] = peak_list
                spectrum_index += 1
        logger.info(f"  Loaded {len(self.spectrum_dict):,d} spectra")

    def _match_AA_novor(
        self,
        target: list[int],
        predicted: list[int],
    ) -> tuple[int, str]:
        """Match amino acids using cumulative mass alignment (Novor algorithm).

        Why Not Simple String Matching?
        -------------------------------
        Consider target "TIDE" vs predicted "TYDE":
        - String matching: 3/4 match (T, D, E), 75% accuracy
        - But I→Y is a mass shift that propagates: after position 1, cumulative
          masses never align again, making downstream matches meaningless.

        The Novor approach uses mass-based alignment instead:
        - Convert sequences to cumulative mass arrays
        - Walk through both, finding positions where masses align (within tolerance)
        - At aligned positions, check if individual AA masses match

        The Algorithm
        -------------
        1. Compute prefix mass sums for both sequences
        2. Two pointers (i for target, j for predicted) start at position 0
        3. At each step:
           - If cumulative masses align (within MASS_TOLERANCE_CUMULATIVE):
             * Check if individual AA masses also align (within MASS_TOLERANCE_AA)
             * If yes, count as match; advance both pointers
             * If no, still advance both (aligned position, wrong AA)
           - If target mass is smaller: advance target pointer (predicted has extra mass)
           - If predicted mass is smaller: advance predicted pointer (target has extra mass)
        4. Return count and match string (1 for match, 0 for mismatch at each position)

        Args:
            target: Target amino acid IDs (integers from config.vocab)
            predicted: Predicted amino acid IDs

        Returns:
            Tuple of (match_count, match_string) where match_string shows "1"/"0"
            at each aligned position
        """
        num_match = 0
        target_len = len(target)
        predicted_len = len(predicted)
        target_mass = [config.mass_ID[x] for x in target]
        target_mass_cum = np.cumsum(target_mass)
        predicted_mass = [config.mass_ID[x] for x in predicted]
        predicted_mass_cum = np.cumsum(predicted_mass)

        i = 0
        j = 0
        aa_match: list[str] = []
        while i < target_len and j < predicted_len:
            if (
                abs(target_mass_cum[i] - predicted_mass_cum[j])
                < config.MASS_TOLERANCE_CUMULATIVE
            ):
                if abs(target_mass[i] - predicted_mass[j]) < config.MASS_TOLERANCE_AA:
                    num_match += 1
                    aa_match.append("1")
                else:
                    aa_match.append("0")
                i += 1
                j += 1
            elif target_mass_cum[i] < predicted_mass_cum[j]:
                i += 1
            else:
                j += 1
                aa_match.append("0")
        aa_match_str = " ".join(aa_match)

        return num_match, aa_match_str

    def _peptide_to_ions(self, peptide: list[str]) -> np.ndarray:
        """Calculate theoretical fragment ion m/z values for a peptide.

        MS/MS Fragmentation Basics
        --------------------------
        In collision-induced dissociation (CID), peptide bonds break to produce:
        - b-ions: N-terminal fragments (keep the amino terminus)
        - y-ions: C-terminal fragments (keep the carboxyl terminus)

        For a peptide of length n, we get n-1 possible cleavage sites, producing
        n-1 b-ions (b1, b2, ..., b_{n-1}) and n-1 y-ions (y1, y2, ..., y_{n-1}).

        Additional complexity:
        - Neutral losses: ions can lose H2O (-18 Da) or NH3 (-17 Da)
        - Multiple charges: ions can be singly (+1) or doubly (+2) charged

        This function generates all ion types: b, b-H2O, b-NH3, y, y-H2O, y-NH3,
        each at charge states +1 and +2. That's 12 ion types per cleavage site.

        Why Generate All Ion Types?
        ---------------------------
        Different peptides and fragmentation conditions produce different ion
        patterns. By generating all possible ions, we can:
        1. Match against the observed spectrum to see which ions are present
        2. Compare target vs predicted peptide ion coverage
        3. Calculate ion-level accuracy metrics

        Args:
            peptide: List of amino acid strings (including modifications)

        Returns:
            2D array of shape (n-1, 12) with m/z values for each ion type
        """
        peptide_mass = self._compute_peptide_mass(peptide)
        prefix_mass = config.mass_AA["_GO"] + np.cumsum(
            [config.mass_AA[aa] for aa in peptide[:-1]]
        )
        suffix_mass = peptide_mass - prefix_mass
        mass_loss = np.array([0, config.mass_H2O, config.mass_NH3])
        b_neutral = prefix_mass.reshape((prefix_mass.size, 1)) - mass_loss.reshape(
            (1, -1)
        )
        y_neutral = suffix_mass.reshape((suffix_mass.size, 1)) - mass_loss.reshape(
            (1, -1)
        )
        by_neutral = np.concatenate([b_neutral, y_neutral], axis=1)
        by_charge1 = by_neutral + config.mass_H
        by_charge2 = (by_neutral + 2 * config.mass_H) / 2
        by_ions = np.concatenate([by_charge1, by_charge2], axis=1)

        return by_ions

    def _match_ion(
        self,
        target: list[str],
        predicted: list[str],
        spectrum: list[tuple[float, float]],
    ) -> tuple[int, int, int, str]:
        """Calculate fragment ion overlap between target, predicted, and spectrum.

        The Three-Way Comparison
        ------------------------
        This function answers: "Do the target and predicted peptides produce
        similar fragmentation in the observed spectrum?"

        We compare three things:
        1. Theoretical ions from target peptide
        2. Theoretical ions from predicted peptide
        3. Observed peaks in the spectrum

        For each observed peak, we check if it matches (within tolerance) any
        theoretical ion from target and/or predicted.

        Key Metrics
        -----------
        - target_ion_count: How many observed peaks match target's theoretical ions
        - predicted_ion_count: How many observed peaks match predicted's theoretical ions
        - matched_ion_count: How many peaks match BOTH target AND predicted

        These enable:
        - Ion recall: matched / target (what fraction of target's signal is explained)
        - Ion precision: matched / predicted (what fraction of prediction is correct)
        - Unmatched list: target ions not covered by prediction (diagnostic)

        Why Ion-Level Matters
        ---------------------
        Two peptides can have different sequences but produce overlapping fragments:
        - PEPTIDE and PEPTYDE differ at position 4, but positions 1-3 produce
          identical b-ions and positions 5-7 produce identical y-ions
        - Ion-level metrics capture this partial correctness

        Args:
            target: Target peptide sequence (list of AA strings)
            predicted: Predicted peptide sequence
            spectrum: Observed peaks as [(m/z, intensity), ...]

        Returns:
            (target_ions, predicted_ions, matched_ions, unmatched_mz_list)
        """
        target_by = self._peptide_to_ions(target).reshape(1, -1)
        predicted_by = self._peptide_to_ions(predicted).reshape(1, -1)
        mz_nby1 = np.array([x[0] for x in spectrum]).reshape(-1, 1)
        target_ion = np.any(np.abs(mz_nby1 - target_by) <= config.ION_TOLERANCE, axis=1)
        predicted_ion = np.any(
            np.abs(mz_nby1 - predicted_by) <= config.ION_TOLERANCE, axis=1
        )
        matched_ion = target_ion * predicted_ion
        mz_nby1 = mz_nby1.flatten()
        unmatched_ion_list = ";".join(
            [
                f"{x:.5f}"
                for x in mz_nby1[np.flatnonzero(target_ion * (1 - predicted_ion))]
            ]
        )
        target_ion_count = int(target_ion.sum())
        predicted_ion_count = int(predicted_ion.sum())
        matched_ion_count = int(matched_ion.sum())

        return (
            target_ion_count,
            predicted_ion_count,
            matched_ion_count,
            unmatched_ion_list,
        )
