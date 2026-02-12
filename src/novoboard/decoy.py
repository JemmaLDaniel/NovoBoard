"""Generate decoy spectra for FDR estimation in de novo peptide sequencing.

The Problem: FDR Estimation Without a Database
-----------------------------------------------
Traditional proteomics uses database search: match spectra against a protein database.
FDR is estimated using decoy databases (reversed/shuffled proteins) - any match to a
decoy is definitionally a false positive.

De novo sequencing predicts peptide sequences directly from spectra without a database.
This breaks the traditional target-decoy paradigm: there's no "decoy database" because
there's no database at all!

The Solution: Spectrum-Level Decoys
-----------------------------------
Instead of decoy proteins, we create DECOY SPECTRA - deliberately corrupted versions
of the input spectra. The key insight:

1. Real spectra contain systematic patterns (fragment ion series) that arise from
   true peptide fragmentation. These patterns allow successful de novo prediction.

2. If we corrupt these patterns - by adding noise, removing signal peaks, or
   shuffling intensities - the spectrum no longer encodes a coherent peptide.

3. Any "successful" prediction on a corrupted spectrum is definitionally incorrect,
   since no real peptide could produce that spectrum.

4. The rate of predictions on corrupted spectra estimates the false positive rate
   among predictions on real spectra.

Sampling Strategies
-------------------
Different corruption strategies test different failure modes:

- **random**: Remove random peaks, replace with random peaks from the dataset.
  Tests: robustness to random noise and missing information.

- **intensity**: Remove the HIGHEST intensity peaks, replace with low-intensity peaks.
  Tests: reliance on dominant fragment ions (b/y series peaks).

- **intensity_mass**: Like intensity, but scales removal to expected peptide length.
  Tests: same as intensity, but normalized for peptide size.

- **permutation**: Keep m/z values but shuffle intensities across peaks.
  Tests: reliance on intensity patterns (which peak is strongest at each m/z).

- **500Da**: Remove all peaks below 500 Da, replace with noise.
  Tests: reliance on low-mass fragments (often b1, b2, y1, y2 ions).

- **distance**: Remove peaks that could be AA mass differences apart.
  Tests: reliance on the "mass ladder" pattern of fragment ions.

The idea is that different de novo algorithms may fail in different ways. Testing
against multiple decoy types reveals which failure modes affect a given algorithm.

Usage
-----
1. Generate decoy MGF files from your target (real) MGF files
2. Run de novo sequencing on BOTH target and decoy files
3. Use the fdr module to compare predictions and estimate FDR
"""

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
    """Build empirical peak distribution for realistic noise generation.

    Why This Matters
    ----------------
    When we remove signal peaks from a spectrum, we need to replace them with
    "noise" peaks to maintain realistic spectral properties (total peak count,
    intensity distribution). Using random numbers would create obviously fake
    spectra.

    Instead, we sample from the ACTUAL peak distribution in the dataset. This
    means noise peaks have realistic m/z values (in the typical mass range) and
    realistic intensities (following the dataset's dynamic range).

    For intensity-based strategies, we also track which peaks get removed so we
    can use similar-intensity peaks as noise (maintaining intensity statistics).

    Args:
        input_mgf: Path to input MGF file
        peak_sampling: Sampling strategy name (determines what to collect)
        sampling_rate: Fraction of peaks to sample (for removed_peaks calculation)

    Returns:
        Tuple of:
        - peaks_distr: All peaks from dataset [[mz, intensity], ...]
        - removed_peaks_distr: Peaks that would be removed (for intensity strategies)
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
        logger.info(f"  Collected {len(peaks_distr):,d} peaks for noise sampling")
    elif peak_sampling == "500Da":
        with open(input_mgf, "r") as f:
            peaks_distr = [
                [float(x) for x in _SPACE_PATTERN.split(line)[:2]]
                for line in f
                if line[0].isdigit()
            ]
        logger.info(f"  Collected {len(peaks_distr):,d} peaks total")
        removed_peaks_distr = [[x, y] for x, y in peaks_distr if x < 500]
        logger.info(
            f"  Found {len(removed_peaks_distr):,d} peaks below 500 Da for noise pool"
        )
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
        logger.info(f"  Collected {len(peaks_distr):,d} peaks total")
        logger.info(
            f"  Low-intensity peaks for noise pool: {len(removed_peaks_distr):,d}"
        )

    return peaks_distr, removed_peaks_distr


def _apply_sampling_strategy(
    peak_list: list[list[float]],
    peak_sampling: str,
    sampling_rate: float,
    peptide_mass: float,
    peaks_distr: list[list[float]],
    removed_peaks_distr: list[list[float]],
) -> tuple[list, list]:
    """Corrupt a spectrum by removing signal peaks and adding noise.

    The goal is to destroy the peptide's "fingerprint" while maintaining realistic
    spectral properties. Each strategy targets different aspects of the signal:

    - random: Removes arbitrary peaks. Tests tolerance to missing information.
    - intensity: Removes high-intensity peaks (likely fragment ions). Tests reliance
      on dominant signals.
    - intensity_mass: Like intensity, but removal scales with peptide length since
      longer peptides have more fragment ions.
    - permutation: Shuffles intensity values across m/z positions. The m/z ladder
      remains but intensity patterns are destroyed.
    - 500Da: Removes all peaks < 500 Da. Low-mass peaks are often most diagnostic
      (b1, b2 ions) so this tests reliance on N-terminal sequencing.
    - distance: Specifically targets peaks that form amino-acid mass ladders. This
      directly attacks the core de novo sequencing signal.

    Noise peaks are drawn from the dataset's empirical peak distribution so they
    look statistically realistic - same m/z range and intensity distribution as
    real peaks, just not correlated with the peptide.

    Args:
        peak_list: List of [mz, intensity] peaks from the spectrum
        peak_sampling: Sampling strategy name
        sampling_rate: Fraction of peaks to KEEP (0.5 = remove half)
        peptide_mass: Precursor peptide mass (for intensity_mass scaling)
        peaks_distr: Pool of all peaks in dataset (for noise sampling)
        removed_peaks_distr: Pool of removed peaks (for intensity strategy)

    Returns:
        Tuple of (peaks_to_keep, noise_peaks_to_add)
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
    """Generate decoy MGF files by corrupting input spectra.

    Main entry point for decoy generation. For each input MGF file, creates a
    corresponding decoy file where spectra have been corrupted using the specified
    strategy. These decoy files are then used with the fdr module to estimate
    false discovery rates.

    The Workflow
    ------------
    1. First pass: Collect the empirical peak distribution from the dataset.
       This provides realistic "noise" peaks to substitute for removed signal.

    2. Second pass: For each spectrum:
       a. Parse the spectrum header (precursor m/z, charge, etc.)
       b. Apply the sampling strategy: remove some peaks, add noise peaks
       c. Write the corrupted spectrum to the output file

    3. The output file has identical structure to the input (same number of
       spectra, same headers) but with corrupted peak lists.

    Strategy Selection
    ------------------
    Different strategies test different aspects of de novo algorithm robustness:

    - 'random': General noise tolerance (default, most commonly used)
    - 'intensity': Reliance on high-intensity peaks (b/y ion series)
    - 'intensity_mass': Like intensity, but peptide-length normalized
    - 'permutation': Reliance on intensity patterns vs. just peak positions
    - '500Da': Reliance on low-mass diagnostic ions
    - 'distance': Reliance on amino acid mass ladder patterns

    For benchmarking a new algorithm, start with 'random'. Use other strategies
    to diagnose specific failure modes.

    Args:
        input_mgf_list: List of input MGF file paths
        peak_sampling: Sampling strategy (see Strategy Selection above)
        sampling_rate: Fraction of peaks to KEEP (0.5 = remove half)
        seed: Random seed for reproducibility (default: 99)
    """
    # Set random seeds at function entry point for reproducibility
    random.seed(seed)
    np.random.seed(seed)

    logger.info("=" * 80)
    logger.info("Generating decoy spectra...")
    logger.info(f"  Corruption strategy: {peak_sampling}")
    logger.info(f"  Peak retention rate: {sampling_rate:.0%}")
    logger.info(f"  Random seed: {seed}")

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

        logger.info("-" * 40)
        logger.info(f"Processing: {input_mgf}")
        logger.info(f"  Output: {output_mgf}")

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

        kept_pct = len(sampling_peaks_distr) / len(decoy_peaks_distr) * 100
        noise_pct = len(noise_peaks_distr) / len(decoy_peaks_distr) * 100
        logger.info(
            f"  Peaks retained from original: {len(sampling_peaks_distr):,d} ({kept_pct:.1f}%)"
        )
        logger.info(
            f"  Noise peaks added: {len(noise_peaks_distr):,d} ({noise_pct:.1f}%)"
        )
        logger.info(f"  Total peaks in decoy spectra: {len(decoy_peaks_distr):,d}")
