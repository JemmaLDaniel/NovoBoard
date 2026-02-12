"""Generate decoy MGF files for FDR estimation."""

from __future__ import annotations

import logging
from pathlib import Path
import re
import random
import numpy as np
from novoboard import config

logger = logging.getLogger(__name__)

# Precompiled regex patterns for performance
_SPACE_PATTERN = re.compile(r" |\r|\n")
_PEPMASS_PATTERN = re.compile(r"=| |\r|\n")
_CHARGE_PATTERN = re.compile(r"=|\+|\r|\n")


def _collect_peaks_distribution(
    input_mgf: str,
    peak_sampling: str,
    sampling_rate: float,
) -> tuple[list[list[float]], list[list[float]]]:
    """Collect peak distribution from MGF file for noise sampling.

    Args:
        input_mgf: Path to input MGF file
        peak_sampling: Sampling strategy name
        sampling_rate: Fraction of peaks to sample

    Returns:
        Tuple of (peaks_distr, removed_peaks_distr)
    """
    peaks_distr: list[list[float]] = []
    removed_peaks_distr: list[list[float]] = []

    if peak_sampling in ("random", "permutation", "distance"):
        with open(input_mgf, "r") as f:
            peaks_distr = [
                [float(x) for x in _SPACE_PATTERN.split(line)[:2]]
                for line in f
                if line[0].isdigit()
            ]
        logger.info(f"len(peaks_distr) = {len(peaks_distr)}")
    elif peak_sampling == "500Da":
        with open(input_mgf, "r") as f:
            peaks_distr = [
                [float(x) for x in _SPACE_PATTERN.split(line)[:2]]
                for line in f
                if line[0].isdigit()
            ]
        logger.info(f"len(peaks_distr) = {len(peaks_distr)}")
        removed_peaks_distr = [[x, y] for x, y in peaks_distr if x < 500]
        logger.info(f"len(removed_peaks_distr) = {len(removed_peaks_distr)}")
    elif peak_sampling in ("intensity", "intensity_mass"):
        with open(input_mgf, "r") as f_in:
            while True:
                line = f_in.readline()
                if not line:  # end of file
                    break
                if line == "\n":  # empty line
                    continue
                peak_list: list[list[float]] = []
                peptide_mass = 0.0
                while "END IONS" not in line:
                    # parse header lines
                    if "BEGIN IONS" in line or "=" in line:
                        line = f_in.readline()
                        if "PEPMASS" in line:
                            mz = float(_PEPMASS_PATTERN.split(line)[1])
                        if "CHARGE" in line:
                            z = float(_CHARGE_PATTERN.split(line)[1])
                            peptide_mass = mz * z - z * config.mass_H
                        continue
                    # parse ions
                    mz_str, intensity_str = _SPACE_PATTERN.split(line)[:2]
                    peak_list.append([float(mz_str), float(intensity_str)])
                    line = f_in.readline()

                # peak removal and noise sampling by intensity
                num_peaks = len(peak_list)
                if peak_sampling == "intensity":
                    num_sampling = int(num_peaks * sampling_rate)
                elif peak_sampling == "intensity_mass":
                    est_len = int(peptide_mass / config.AVG_AA_MASS)
                    num_noise = int(
                        min((est_len - 1) * 2, num_peaks) * (1 - sampling_rate)
                    )
                    num_sampling = num_peaks - num_noise
                peak_list_sorted = sorted(peak_list, key=lambda x: x[1])
                removed_peaks = peak_list_sorted[num_sampling:]
                removed_peaks_distr += removed_peaks
                peaks_distr += peak_list
        logger.info(f"len(peaks_distr) = {len(peaks_distr)}")
        logger.info(f"len(removed_peaks_distr) = {len(removed_peaks_distr)}")

    return peaks_distr, removed_peaks_distr


def _apply_sampling_strategy(
    peak_list: list[list[float]],
    peak_sampling: str,
    sampling_rate: float,
    peptide_mass: float,
    peaks_distr: list[list[float]],
    removed_peaks_distr: list[list[float]],
) -> tuple[list, list]:
    """Apply the specified sampling strategy to peaks.

    Args:
        peak_list: List of [mz, intensity] peaks
        peak_sampling: Sampling strategy name
        sampling_rate: Fraction of peaks to sample
        peptide_mass: Peptide mass for intensity_mass strategy
        peaks_distr: Distribution of peaks for noise sampling
        removed_peaks_distr: Distribution of removed peaks for noise sampling

    Returns:
        Tuple of (sampling_peaks, noise_peaks)
    """
    num_peaks = len(peak_list)
    num_sampling = int(num_peaks * sampling_rate)
    num_noise = num_peaks - num_sampling

    sampling_peaks: list = []
    noise_peaks: list = []

    if peak_sampling == "random":
        sampling_peaks = random.sample(peak_list, num_sampling)
        noise_peaks = random.sample(peaks_distr, num_noise)
    elif peak_sampling == "intensity":
        peak_list_sorted = sorted(peak_list, key=lambda x: x[1])
        sampling_peaks = peak_list_sorted[:num_sampling]
        noise_peaks = random.sample(removed_peaks_distr, num_noise)
    elif peak_sampling == "intensity_mass":
        est_len = int(peptide_mass / config.AVG_AA_MASS)
        num_noise = int(min((est_len - 1) * 2, num_peaks) * (1 - sampling_rate))
        num_sampling = num_peaks - num_noise
        peak_list_sorted = sorted(peak_list, key=lambda x: x[1])
        sampling_peaks = peak_list_sorted[:num_sampling]
        noise_peaks = random.sample(removed_peaks_distr, num_noise)
    elif peak_sampling == "permutation":
        mz_list = [x[0] for x in peak_list]
        intensity_list = [x[1] for x in peak_list]
        sampling_peaks = []
        random.shuffle(intensity_list)
        noise_peaks = list(zip(mz_list, intensity_list))
    elif peak_sampling == "500Da":
        sampling_peaks = [[x, y] for x, y in peak_list if x >= 500]
        num_sampling = len(sampling_peaks)
        num_noise = num_peaks - num_sampling
        noise_peaks = random.sample(removed_peaks_distr, num_noise)
    elif peak_sampling == "distance":
        mz_array = np.array([peak[0] for peak in peak_list])
        pair_distance = np.absolute(
            np.reshape(mz_array, (num_peaks, 1)) - np.reshape(mz_array, (1, num_peaks))
        )
        aa_masses = config.mass_ID_np[3:].reshape(1, -1)
        pair_aa_match = np.absolute(np.expand_dims(pair_distance, axis=2) - aa_masses)
        pair_aa_match = np.any(pair_aa_match <= config.ION_TOLERANCE, axis=2)
        match_peak_indices = list(np.flatnonzero(np.any(pair_aa_match, axis=1)))
        num_noise = int(len(match_peak_indices) * (1 - sampling_rate))
        removed_indices = random.sample(match_peak_indices, num_noise)
        sampling_peaks = [
            peak for index, peak in enumerate(peak_list) if index not in removed_indices
        ]
        noise_peaks = random.sample(peaks_distr, num_noise)
    else:
        sampling_peaks = peak_list
        noise_peaks = []

    return sampling_peaks, noise_peaks


def generate_decoy_mgf(
    input_mgf_list: list[str],
    peak_sampling: str = "random",
    sampling_rate: float = config.DEFAULT_SAMPLING_RATE,
    seed: int = 99,
) -> None:
    """Generate decoy MGF files with specified peak sampling strategy.

    Creates decoy spectra by sampling/shuffling peaks from target spectra,
    used for estimating false discovery rate in de novo sequencing.

    Args:
        input_mgf_list: List of input MGF file paths
        peak_sampling: Sampling strategy. Options:
            - 'random': Random peak sampling
            - 'intensity': Sample based on peak intensity
            - 'intensity_mass': Sample based on intensity and peptide mass
            - 'permutation': Shuffle peak intensities
            - '500Da': Remove peaks under 500 Da
            - 'distance': Remove peaks matching AA mass differences
        sampling_rate: Fraction of peaks to keep (0.0-1.0)
        seed: Random seed for reproducibility (default: 99)
    """
    # Set random seeds at function entry point for reproducibility
    random.seed(seed)
    np.random.seed(seed)

    logger.info(f"peak_sampling = {peak_sampling}")
    logger.info(f"sampling_rate = {sampling_rate}")
    logger.info(f"seed = {seed}")

    for input_mgf in input_mgf_list:
        # Build clean output filename: strip suffix, append decoy info, add .mgf
        input_path = Path(input_mgf)
        stem = input_path.stem  # filename without extension
        parent = input_path.parent

        if peak_sampling == "permutation":
            output_mgf = str(parent / f"{stem}_permutation.mgf")
        elif peak_sampling == "500Da":
            output_mgf = str(parent / f"{stem}_500Da.mgf")
        else:
            output_mgf = str(parent / f"{stem}_decoy_{sampling_rate:.2f}.mgf")

        logger.info(f"Generating decoy: {input_mgf} -> {output_mgf}")

        # Collect peak distributions for noise sampling
        peaks_distr, removed_peaks_distr = _collect_peaks_distribution(
            input_mgf, peak_sampling, sampling_rate
        )

        sampling_peaks_distr: list[list[float]] = []
        noise_peaks_distr: list[list[float]] = []
        decoy_peaks_distr: list[list[float]] = []

        with open(input_mgf, "r") as f_in:
            with open(output_mgf, "w") as f_out:
                while True:
                    line = f_in.readline()
                    if not line:  # end of file
                        break
                    if line == "\n":  # empty line
                        continue
                    peak_list: list[list[float]] = []
                    peptide_mass = 0.0
                    while "END IONS" not in line:
                        # parse header lines
                        if "BEGIN IONS" in line or "=" in line:
                            f_out.write(line)
                            line = f_in.readline()
                            if "PEPMASS" in line:
                                mz = float(_PEPMASS_PATTERN.split(line)[1])
                            if "CHARGE" in line:
                                z = float(_CHARGE_PATTERN.split(line)[1])
                                peptide_mass = mz * z - z * config.mass_H
                            continue
                        # parse ions
                        mz_str, intensity_str = _SPACE_PATTERN.split(line)[:2]
                        peak_list.append([float(mz_str), float(intensity_str)])
                        line = f_in.readline()

                    # Apply sampling strategy
                    sampling_peaks, noise_peaks = _apply_sampling_strategy(
                        peak_list,
                        peak_sampling,
                        sampling_rate,
                        peptide_mass,
                        peaks_distr,
                        removed_peaks_distr,
                    )

                    sampling_peaks_distr += sampling_peaks
                    noise_peaks_distr += noise_peaks
                    decoy_peaks_distr += sampling_peaks + noise_peaks

                    # write ion lines
                    sorted_peaks = sorted(
                        sampling_peaks + noise_peaks, key=lambda x: x[0]
                    )
                    for x, y in sorted_peaks:
                        f_out.write(f"{x:.5f} {y:.5f}\n")
                    f_out.write(line)  # END IONS line
                    f_out.write(f_in.readline())  # empty line between spectra

        logger.info(f"len(sampling_peaks_distr) = {len(sampling_peaks_distr)}")
        logger.info(f"len(noise_peaks_distr) = {len(noise_peaks_distr)}")
        logger.info(f"len(decoy_peaks_distr) = {len(decoy_peaks_distr)}")
