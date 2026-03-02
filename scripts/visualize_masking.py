#!/usr/bin/env python3
"""Visualize original spectra alongside their decoy counterparts from pre-generated files.

This script samples spectra from an original MGF file and displays them alongside
the corresponding spectra from one or more decoy MGF files.

Usage:
    python scripts/visualize_masking.py \
        --mgf-file /path/to/gluc.mgf \
        --decoy-files /path/to/gluc_decoy_0.10.mgf /path/to/gluc_decoy_0.50.mgf \
        --n-spectra 3 --seed 42
"""

from __future__ import annotations

import argparse
import logging
import random
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Default paths
DEFAULT_OUTPUT_DIR = Path("/home/j-daniel/repos/NovoBoard/fig/masking_examples")

# Regex patterns for parsing MGF
SPACE_PATTERN = re.compile(r" |\r|\n")
PEPMASS_PATTERN = re.compile(r"=| |\r|\n")
CHARGE_PATTERN = re.compile(r"=|\+|\r|\n")


def parse_mgf_spectra(mgf_path: Path, max_spectra: int | None = None) -> list[dict]:
    """Parse spectra from an MGF file.

    Args:
        mgf_path: Path to MGF file
        max_spectra: Maximum number of spectra to parse (None for all)

    Returns:
        List of spectrum dictionaries with keys: title, pepmass, charge, peaks
    """
    spectra = []

    with open(mgf_path, "r") as f:
        current_spectrum: dict | None = None

        for line in f:
            line = line.strip()

            if line == "BEGIN IONS":
                current_spectrum = {"peaks": [], "title": "", "pepmass": 0, "charge": 0}
            elif line == "END IONS":
                if current_spectrum and current_spectrum["peaks"]:
                    spectra.append(current_spectrum)
                    if max_spectra is not None and len(spectra) >= max_spectra:
                        break
                current_spectrum = None
            elif current_spectrum is not None:
                if line.startswith("TITLE="):
                    current_spectrum["title"] = line[6:]
                elif line.startswith("PEPMASS="):
                    parts = PEPMASS_PATTERN.split(line)
                    current_spectrum["pepmass"] = float(parts[1])
                elif line.startswith("CHARGE="):
                    parts = CHARGE_PATTERN.split(line)
                    current_spectrum["charge"] = int(parts[1])
                elif line and line[0].isdigit():
                    parts = SPACE_PATTERN.split(line)
                    if len(parts) >= 2:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        current_spectrum["peaks"].append((mz, intensity))

    return spectra


def extract_masking_rate_from_filename(filename: str) -> str:
    """Extract masking rate label from decoy filename.

    Args:
        filename: Decoy filename like 'gluc_decoy_0.10.mgf'

    Returns:
        Label string like '10% masked' or the filename if pattern not found
    """
    # Look for pattern like _decoy_0.XX
    match = re.search(r"_decoy_(\d+\.?\d*)", filename)
    if match:
        rate = float(match.group(1))
        # The rate in filename is fraction of peaks KEPT
        # So 0.10 means 10% kept = 90% masked
        pct_masked = int((1 - rate) * 100)
        return f"{pct_masked}% masked ({int(rate * 100)}% kept)"
    return filename


def plot_spectrum_comparison(
    original_spectrum: dict,
    decoy_spectra: list[dict],
    decoy_labels: list[str],
    spectrum_idx: int,
    mgf_name: str,
    output_path: Path,
) -> None:
    """Plot original spectrum alongside decoy versions.

    Args:
        original_spectrum: Original spectrum dict with peaks, pepmass, charge
        decoy_spectra: List of decoy spectrum dicts
        decoy_labels: Labels for each decoy (e.g., '90% masked')
        spectrum_idx: Index of spectrum in sample
        mgf_name: Name of original MGF file (for title)
        output_path: Output file path
    """
    n_plots = 1 + len(decoy_spectra)
    fig, axes = plt.subplots(n_plots, 1, figsize=(12, 3 * n_plots), sharex=True)

    if n_plots == 1:
        axes = [axes]

    # Color scheme
    original_color = "#2E86AB"
    noise_color = "#F39C12"

    original_peaks = original_spectrum["peaks"]

    def normalize(peaks):
        """Normalize peaks to their own maximum - this is what the model sees."""
        if not peaks:
            return []
        max_int = max(p[1] for p in peaks)
        if max_int == 0:
            return [(mz, 0) for mz, _ in peaks]
        return [(mz, intensity / max_int * 100) for mz, intensity in peaks]

    # Plot original spectrum
    ax = axes[0]
    orig_normalized = normalize(original_peaks)
    mz_values = [p[0] for p in orig_normalized]
    intensities = [p[1] for p in orig_normalized]

    ax.bar(mz_values, intensities, width=1.5, color=original_color, alpha=0.8)
    ax.set_ylabel("Relative\nIntensity (%)", fontsize=10)
    ax.set_title(
        f"Original Spectrum ({len(original_peaks)} peaks) | {mgf_name} | "
        f"m/z: {original_spectrum['pepmass']:.2f} | Charge: {original_spectrum['charge']}+",
        fontsize=11,
        fontweight="bold",
    )
    ax.set_ylim(0, 110)
    ax.axhline(y=0, color="black", linewidth=0.5)

    # Build set of original m/z values for identifying kept vs noise peaks
    original_mz_set = set(round(p[0], 2) for p in original_peaks)

    # Plot decoy spectra
    for i, (decoy_spectrum, label) in enumerate(zip(decoy_spectra, decoy_labels)):
        ax = axes[i + 1]
        decoy_peaks = decoy_spectrum["peaks"]
        decoy_normalized = normalize(decoy_peaks)

        # Separate kept peaks from noise peaks
        kept_mz, kept_int = [], []
        noise_mz, noise_int = [], []

        for mz, intensity in decoy_normalized:
            if round(mz, 2) in original_mz_set:
                kept_mz.append(mz)
                kept_int.append(intensity)
            else:
                noise_mz.append(mz)
                noise_int.append(intensity)

        # Plot kept peaks
        ax.bar(
            kept_mz, kept_int, width=1.5, color=original_color, alpha=0.8, label="Kept"
        )
        # Plot noise peaks
        ax.bar(
            noise_mz, noise_int, width=1.5, color=noise_color, alpha=0.8, label="Noise"
        )

        ax.set_ylabel("Relative\nIntensity (%)", fontsize=10)
        n_kept = len(kept_mz)
        n_noise = len(noise_mz)
        ax.set_title(
            f"Decoy Spectrum | {label} | {n_kept} kept + {n_noise} noise = {n_kept + n_noise} total",
            fontsize=11,
            fontweight="bold",
        )
        ax.set_ylim(0, 110)
        ax.axhline(y=0, color="black", linewidth=0.5)
        ax.legend(loc="upper right", fontsize=9)

    axes[-1].set_xlabel("m/z", fontsize=11)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved: {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Visualize original spectra alongside pre-generated decoy spectra",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Compare original with multiple decoy rates
    python scripts/visualize_masking.py \\
        --mgf-file /path/to/gluc.mgf \\
        --decoy-files /path/to/gluc_decoy_0.10.mgf /path/to/gluc_decoy_0.50.mgf \\
        --n-spectra 3 --seed 42

    # Single decoy comparison
    python scripts/visualize_masking.py \\
        --mgf-file /path/to/hepg2.mgf \\
        --decoy-files /path/to/hepg2_decoy_0.30.mgf \\
        --n-spectra 5
        """,
    )
    parser.add_argument(
        "--mgf-file",
        type=Path,
        required=True,
        help="Path to original MGF file",
    )
    parser.add_argument(
        "--decoy-files",
        type=Path,
        nargs="+",
        required=True,
        help="Path(s) to decoy MGF file(s)",
    )
    parser.add_argument(
        "--n-spectra",
        type=int,
        default=3,
        help="Number of spectra to sample. Default: 3",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory for plots. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility. Default: 42",
    )

    args = parser.parse_args()

    # Validate inputs
    if not args.mgf_file.exists():
        raise FileNotFoundError(f"MGF file not found: {args.mgf_file}")

    for decoy_file in args.decoy_files:
        if not decoy_file.exists():
            raise FileNotFoundError(f"Decoy file not found: {decoy_file}")

    # Set random seed
    random.seed(args.seed)
    np.random.seed(args.seed)

    # Parse original MGF
    logger.info(f"Parsing original MGF: {args.mgf_file}")
    original_spectra = parse_mgf_spectra(args.mgf_file)
    logger.info(f"Parsed {len(original_spectra)} spectra from original file")

    # Parse decoy MGFs
    decoy_spectra_lists = []
    decoy_labels = []
    for decoy_file in args.decoy_files:
        logger.info(f"Parsing decoy MGF: {decoy_file}")
        decoy_spectra = parse_mgf_spectra(decoy_file)
        logger.info(f"Parsed {len(decoy_spectra)} spectra from {decoy_file.name}")

        if len(decoy_spectra) != len(original_spectra):
            logger.warning(
                f"Decoy file {decoy_file.name} has {len(decoy_spectra)} spectra, "
                f"but original has {len(original_spectra)}. Using minimum."
            )

        decoy_spectra_lists.append(decoy_spectra)
        decoy_labels.append(extract_masking_rate_from_filename(decoy_file.name))

    # Determine valid spectrum indices
    min_spectra = min(len(original_spectra), *[len(ds) for ds in decoy_spectra_lists])

    # Sample spectrum indices
    if min_spectra < args.n_spectra:
        logger.warning(f"Only {min_spectra} spectra available, using all")
        sampled_indices = list(range(min_spectra))
    else:
        sampled_indices = sorted(random.sample(range(min_spectra), args.n_spectra))

    logger.info(f"Sampled spectrum indices: {sampled_indices}")

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Get MGF name for plot titles and filenames
    mgf_name = args.mgf_file.stem

    # Build decoy rate string for filename
    decoy_rates_str = "_".join(
        [df.stem.split("_decoy_")[-1] for df in args.decoy_files]
    )

    # Process each sampled spectrum
    for i, idx in enumerate(sampled_indices):
        logger.info(f"Processing spectrum {i+1}/{len(sampled_indices)} (index {idx})")

        original_spectrum = original_spectra[idx]
        decoy_spectrum_list = [ds[idx] for ds in decoy_spectra_lists]

        # Generate output filename using original MGF name
        output_path = (
            args.output_dir / f"{mgf_name}_spectrum{idx}_decoys_{decoy_rates_str}.png"
        )

        # Plot
        plot_spectrum_comparison(
            original_spectrum=original_spectrum,
            decoy_spectra=decoy_spectrum_list,
            decoy_labels=decoy_labels,
            spectrum_idx=idx,
            mgf_name=mgf_name,
            output_path=output_path,
        )

    logger.info(f"\nDone! Plots saved to: {args.output_dir}")
    logger.info(f"Decoy files compared: {[f.name for f in args.decoy_files]}")


if __name__ == "__main__":
    main()
