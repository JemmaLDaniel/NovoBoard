"""Calculate fragment ion, amino acid, and peptide accuracies."""

from __future__ import annotations

import logging
import re
import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from novoboard import config

# Precompiled regex patterns for performance
_SCAN_SPLIT_PATTERN = re.compile(r';|\r|\n')
_CSV_SPLIT_PATTERN = re.compile(r',|\r|\n')
_MGF_FIELD_PATTERN = re.compile(r'[=\r\n]')
_PEAK_SPLIT_PATTERN = re.compile(r' |\t|\r|\n')


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
    """Parse a raw peptide sequence with modifications.
    
    Args:
        raw_sequence: Peptide sequence string with optional modifications
                     (e.g., "PEPTC(+57.02)DE")
    
    Returns:
        Tuple of (success, peptide_list) where success is True if all
        amino acids and modifications are recognized.
    """
    raw_sequence_len = len(raw_sequence)
    peptide: list[str] = []
    index = 0
    while index < raw_sequence_len:
        if raw_sequence[index] == "(":
            if peptide[-1] == "C" and raw_sequence[index:index + 8] == "(+57.02)":
                peptide[-1] = "C(Carbamidomethylation)"
                index += 8
            elif peptide[-1] == 'M' and raw_sequence[index:index + 8] == "(+15.99)":
                peptide[-1] = 'M(Oxidation)'
                index += 8
            elif peptide[-1] == 'N' and raw_sequence[index:index + 7] == "(+0.98)":
                peptide[-1] = 'N(Deamidation)'
                index += 7
            elif peptide[-1] == 'Q' and raw_sequence[index:index + 7] == "(+0.98)":
                peptide[-1] = 'Q(Deamidation)'
                index += 7
            elif peptide[-1] == 'S' and raw_sequence[index:index + 8] == "(+79.97)":
                peptide[-1] = "S(Phosphorylation)"
                index += 8
            elif peptide[-1] == 'T' and raw_sequence[index:index + 8] == "(+79.97)":
                peptide[-1] = "T(Phosphorylation)"
                index += 8
            elif peptide[-1] == 'Y' and raw_sequence[index:index + 8] == "(+79.97)":
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
    """Calculate accuracy metrics for de novo peptide sequencing.
    
    Compares predicted peptide sequences against target (database) sequences
    and calculates fragment ion, amino acid, and peptide-level accuracy metrics.
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
        logger.info("WorkerTest.__init__()")

        self.MZ_MAX: float = config.MZ_MAX

        self.target_file = target_file
        self.predicted_file = predicted_file
        self.spectrum_file = spectrum_file
        self.accuracy_file = f"{predicted_file}.accuracy"
        self.denovo_only_file = f"{predicted_file}.denovo_only"
        self.scan2fea_file = f"{predicted_file}.scan2fea"
        self.multifea_file = f"{predicted_file}.multifea"
        self.col_score = col_score
        self.col_aa_score = col_aa_score
        logger.info(f"target_file = {self.target_file}")
        logger.info(f"predicted_file = {self.predicted_file}")
        logger.info(f"spectrum_file = {self.spectrum_file}")
        logger.info(f"accuracy_file = {self.accuracy_file}")
        logger.info(f"denovo_only_file = {self.denovo_only_file}")
        logger.info(f"scan2fea_file = {self.scan2fea_file}")
        logger.info(f"multifea_file = {self.multifea_file}")

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
            target_simplified = ['M' if x == 'M(Oxidation)' else x for x in target_simplified]
            target_simplified = ['N' if x == 'N(Deamidation)' else x for x in target_simplified]
            target_simplified = ['Q' if x == 'Q(Deamidation)' else x for x in target_simplified]
            if target_simplified in db_peptide_list:
                target_dict_db[feature_id] = target
            else:
                logger.warning(f"target not found: {target_simplified}")
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
        best_aa_match = ''
        best_predicted_sequence = predicted["sequence"][0]
        best_predicted_score = predicted["score"][0]
        best_predicted_aa_score = predicted["aa_score"][0]
        
        for pred_seq, pred_score, pred_aa_score in zip(
            predicted["sequence"], predicted["score"], predicted["aa_score"]
        ):
            predicted_AA_id = [config.vocab[x] for x in pred_seq]
            target_AA_id = [config.vocab[x] for x in target]
            recall_AA, aa_match = self._match_AA_novor(target_AA_id, predicted_AA_id)
            
            if (recall_AA > best_recall_AA
                    or (recall_AA == best_recall_AA and pred_score > best_predicted_score)):
                best_recall_AA = recall_AA
                best_aa_match = aa_match
                best_predicted_sequence = pred_seq[:]
                best_predicted_score = pred_score
                best_predicted_aa_score = pred_aa_score
        
        return (best_recall_AA, best_aa_match, best_predicted_sequence, 
                best_predicted_score, best_predicted_aa_score)

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
                        multifea_dict[feature_id].append(f'{scan_id}:{feature_count}')
                    else:
                        multifea_dict[feature_id] = [f'{scan_id}:{feature_count}']

        # Write scan2fea file
        with open(self.scan2fea_file, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(["scan_id", "feature_count", "feature_list"])
            for scan_id, value in scan_dict.items():
                writer.writerow([scan_id, value["feature_count"], 
                                 ";".join(value["feature_list"])])

        # Write multifea file
        with open(self.multifea_file, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(["feature_id", "scan_list"])
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
        logger.info(f"target_count_total = {target_count_total:d}")
        logger.info(f"target_len_total = {target_len_total:d}")
        logger.info(f"target_count_db = {target_count_db:d}")
        logger.info(f"target_len_db = {target_len_db:d}")
        logger.info(f"target_count_db_mass: {target_count_db_mass:d}")
        logger.info(f"target_len_db_mass: {target_len_db_mass:d}")

        logger.info(f"predicted_count_mass: {predicted_count_mass:d}")
        logger.info(f"predicted_count_mass_db: {predicted_count_mass_db:d}")
        logger.info(f"predicted_len_mass_db: {predicted_len_mass_db:d}")
        logger.info(f"predicted_only: {predicted_only:d}")

        logger.info(f"recall_AA_total = {recall_AA_total / target_len_total:.4f}")
        logger.info(f"recall_AA_db = {recall_AA_total / target_len_db:.4f}")
        logger.info(f"recall_AA_db_mass = {recall_AA_total / target_len_db_mass:.4f}")
        logger.info(f"recall_peptide_total = {recall_peptide_total / target_count_total:.4f}")
        logger.info(f"recall_peptide_db = {recall_peptide_total / target_count_db:.4f}")
        logger.info(f"recall_peptide_db_mass = {recall_peptide_total / target_count_db_mass:.4f}")
        logger.info(f"precision_AA_mass_db  = {recall_AA_total / predicted_len_mass_db:.4f}")
        logger.info(f"precision_peptide_mass_db  = {recall_peptide_total / predicted_count_mass_db:.4f}")

        logger.info(f"recall_ion = {matched_ion_total / target_ion_total:.4f}")
        logger.info(f"precision_ion = {matched_ion_total / predicted_ion_total:.4f}")
        logger.info(f"recall_all_peptide_ions = {recall_all_peptide_ions_total / target_count_db_mass:.4f}")
        logger.info(f"target_ion_total = {target_ion_total}")
        logger.info(f"matched_ion_total = {matched_ion_total}")

    def test_accuracy(self, db_peptide_list: list[list[str]] | None = None) -> None:
        """Calculate accuracy metrics between predicted and target peptides.
        
        Args:
            db_peptide_list: Optional list of peptides to filter targets
        """
        logger.info("=" * 80)
        logger.info("WorkerTest.test_accuracy()")

        # write the accuracy of predicted peptides
        accuracy_handle = open(self.accuracy_file, 'w')
        header_list = ["feature_id",
                       "feature_area",
                       "target_sequence",
                       "predicted_sequence",
                       "predicted_score",
                       "predicted_aa_score",
                       "recall_AA",
                       "aa_match",
                       "predicted_len",
                       "target_len",
                       "target_ion",
                       "matched_ion",
                       "unmatched_ion_list",
                       "scan_list_middle",
                       "scan_list_original"]
        header_row = "\t".join(header_list)
        print(header_row, file=accuracy_handle, end="\n")

        # write denovo_only peptides
        denovo_only_handle = open(self.denovo_only_file, 'w')
        header_list = ["feature_id",
                       "feature_area",
                       "predicted_sequence",
                       "predicted_score",
                       "predicted_score_max",
                       "scan_list_middle",
                       "scan_list_original"]
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
                recall_AA, aa_match, predicted_sequence, predicted_score, predicted_aa_score = \
                    self._find_best_prediction(target, predicted)

                recall_AA_total += recall_AA
                if recall_AA == target_len:
                    recall_peptide_total += 1
                predicted_len = len(predicted_sequence)
                predicted_len_mass_db += predicted_len

                if feature_id in self.spectrum_dict:
                    target_ion, predicted_ion, matched_ion, unmatched_ion_list = self._match_ion(target, predicted_sequence, self.spectrum_dict[feature_id])
                else:
                    target_ion, predicted_ion, matched_ion, unmatched_ion_list = 0, 0, 0, ''
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
                print_list = [feature_id,
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
                              feature_scan_list_original]
                print_row = "\t".join(print_list)
                print(print_row, file=accuracy_handle, end="\n")
            else:
                predicted_only += 1
                predicted_sequence_out = ';'.join([','.join(x) for x in predicted["sequence"]])
                predicted_score_out = ';'.join([f'{x:.2f}' for x in predicted["score"]])
                if predicted["score"]:
                    predicted_score_max = f'{np.max(predicted["score"]):.2f}'
                else:
                    predicted_score_max = ''
                print_list = [feature_id,
                              feature_area,
                              predicted_sequence_out,
                              predicted_score_out,
                              predicted_score_max,
                              feature_scan_list_middle,
                              feature_scan_list_original]
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
        peptide_mass = (config.mass_N_terminus
                        + sum(config.mass_AA[aa] for aa in peptide)
                        + config.mass_C_terminus)

        return peptide_mass

    def _get_predicted_peaks_11(self) -> None:
        """Read predicted peptides from PEAKS output CSV file."""
        logger.info("=" * 80)
        logger.info("WorkerTest._get_predicted_peaks_11()")

        predicted_list: list[dict] = []
        with open(self.predicted_file, 'r') as handle:
            csv_reader = csv.DictReader(handle)
            for row in csv_reader:
                predicted: dict = {}
                predicted["feature_id"] = row[col_source_file].split('.mgf')[0] + "||" + row[col_scan_list]
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
                predicted_list.append(predicted)

        self.predicted_list = predicted_list

    def _get_target(self) -> None:
        """Read target peptides from database search CSV file."""
        logger.info("=" * 80)
        logger.info("WorkerTest._get_target()")

        target_dict: dict[str, list[str]] = {}
        with open(self.target_file, 'r') as handle:
            header_line = handle.readline()
            header = [x.strip('"') for x in header_line.strip().split(',')]
            logger.debug(f"Header: {header}")
            raw_sequence_index = header.index(col_raw_sequence)
            source_file_index = header.index(col_source_file)
            scan_index = header.index(col_scan_list)

            for line in handle:
                line_parts = [x.strip('"') for x in _CSV_SPLIT_PATTERN.split(line)]
                feature_id = line_parts[source_file_index] + "||" + line_parts[scan_index]
                raw_sequence = line_parts[raw_sequence_index]
                assert raw_sequence, "Error: wrong target format."
                okay, peptide = parse_raw_sequence(raw_sequence)
                if not okay:
                    # skip unknown mod
                    continue
                target_dict[feature_id] = peptide
        self.target_dict = target_dict

    def _get_spectra(self) -> None:
        """Read spectra from MGF file."""
        logger.info("=" * 80)
        logger.info("WorkerTest._get_spectra()")

        # Extract default source file name from spectrum file path
        default_source_file = Path(self.spectrum_file).stem
        spectrum_index = 0

        with open(self.spectrum_file, 'r') as f_in:
            while True:
                line = f_in.readline()
                if not line:  # end of file
                    break
                if line == '\n':  # empty line
                    continue
                peak_list: list[tuple[float, float]] = []
                # Initialize with defaults for each spectrum
                source_file = default_source_file
                scan = str(spectrum_index)
                while "END IONS" not in line:
                    # parse header lines
                    if 'BEGIN IONS' in line or '=' in line:
                        if "TITLE=" in line:
                            title_value = _MGF_FIELD_PATTERN.split(line)[1]
                            # Try to extract source file from title
                            if '\\' in title_value or '.raw' in title_value:
                                source_file = title_value.split('\\')[-1].split('.raw')[0]
                            else:
                                source_file = default_source_file
                        if line[:6] == "SCANS=":
                            scan = _MGF_FIELD_PATTERN.split(line)[1]
                        line = f_in.readline()
                        continue
                    # parse ions
                    mz, intensity = _PEAK_SPLIT_PATTERN.split(line)[:2]
                    peak_list.append((float(mz), float(intensity)))
                    line = f_in.readline()
                feature_id = f'{source_file}||{scan}'
                self.spectrum_dict[feature_id] = peak_list
                spectrum_index += 1
        logger.info(f"len(self.spectrum_dict) = {len(self.spectrum_dict)}")

    def _match_AA_novor(
        self,
        target: list[int],
        predicted: list[int],
    ) -> tuple[int, str]:
        """Match amino acids between target and predicted using cumulative mass.
        
        Uses the Novor-style matching algorithm based on cumulative mass
        alignment with tolerance thresholds.
        
        Args:
            target: List of target amino acid IDs
            predicted: List of predicted amino acid IDs
        
        Returns:
            Tuple of (number of matches, match string)
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
            if abs(target_mass_cum[i] - predicted_mass_cum[j]) < config.MASS_TOLERANCE_CUMULATIVE:
                if abs(target_mass[i] - predicted_mass[j]) < config.MASS_TOLERANCE_AA:
                    num_match += 1
                    aa_match.append('1')
                else:
                    aa_match.append('0')
                i += 1
                j += 1
            elif target_mass_cum[i] < predicted_mass_cum[j]:
                i += 1
            else:
                j += 1
                aa_match.append('0')
        aa_match_str = ' '.join(aa_match)

        return num_match, aa_match_str

    def _peptide_to_ions(self, peptide: list[str]) -> np.ndarray:
        """Calculate theoretical fragment ion m/z values for a peptide.
        
        Generates b and y ions with neutral losses (H2O, NH3) for charges 1 and 2.
        
        Args:
            peptide: List of amino acid strings
        
        Returns:
            2D numpy array of ion m/z values, shape (len(peptide)-1, num_ion_types)
        """
        peptide_mass = self._compute_peptide_mass(peptide)
        prefix_mass = config.mass_AA['_GO'] + np.cumsum([config.mass_AA[aa] for aa in peptide[:-1]])
        suffix_mass = peptide_mass - prefix_mass
        mass_loss = np.array([0, config.mass_H2O, config.mass_NH3])
        b_neutral = prefix_mass.reshape((prefix_mass.size, 1)) - mass_loss.reshape((1, -1))
        y_neutral = suffix_mass.reshape((suffix_mass.size, 1)) - mass_loss.reshape((1, -1))
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
        """Match fragment ions between target, predicted peptides and spectrum.
        
        Args:
            target: Target peptide sequence
            predicted: Predicted peptide sequence
            spectrum: List of (m/z, intensity) tuples
        
        Returns:
            Tuple of (target_ion_count, predicted_ion_count, matched_count, unmatched_list)
        """
        target_by = self._peptide_to_ions(target).reshape(1, -1)
        predicted_by = self._peptide_to_ions(predicted).reshape(1, -1)
        mz_nby1 = np.array([x[0] for x in spectrum]).reshape(-1, 1)
        target_ion = np.any(np.abs(mz_nby1 - target_by) <= config.ION_TOLERANCE, axis=1)
        predicted_ion = np.any(np.abs(mz_nby1 - predicted_by) <= config.ION_TOLERANCE, axis=1)
        matched_ion = target_ion * predicted_ion
        mz_nby1 = mz_nby1.flatten()
        unmatched_ion_list = ';'.join([f'{x:.5f}' for x in mz_nby1[np.flatnonzero(target_ion * (1 - predicted_ion))]])
        target_ion_count = int(target_ion.sum())
        predicted_ion_count = int(predicted_ion.sum())
        matched_ion_count = int(matched_ion.sum())

        return target_ion_count, predicted_ion_count, matched_ion_count, unmatched_ion_list
