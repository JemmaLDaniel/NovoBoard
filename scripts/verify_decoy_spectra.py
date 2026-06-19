r"""Verify decoy spectrum generation quality and intensity preservation.

Compares original spectra against decoy spectra at different masking rates to verify:
- Intensity scales are preserved (no unintended normalization)
- Expected fraction of peaks are retained at each masking rate
- Intensity distributions match between original and decoy spectra

Generates:
- Spectrum comparison plots (original vs decoy at various masking rates)
- Intensity distribution histograms
- Summary statistics table

Usage:
    python scripts/verify_decoy_spectra.py
    python scripts/verify_decoy_spectra.py --origin hepg2
    python scripts/verify_decoy_spectra.py --all-origins --output-dir fig/spectrum_distributions
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matchms.importing import load_from_mgf


# Default paths
DEFAULT_MGF_DIR = Path(
    "/home/j-daniel/Documents/winnow/novoboard/benchmark_data/instanovo_inputs/per_origin_mgf"
)
DEFAULT_FIG_OUTPUT_DIR = Path("fig/spectrum_distributions")
DEFAULT_CSV_OUTPUT_DIR = Path("analysis")

# Retention rates to check (fraction of peaks KEPT)
# E.g., 0.20 = 20% kept = 80% masked
DEFAULT_KEPT_RATES = [0.20, 0.50, 0.80]

# Origins available for analysis
ORIGINS = [
    "gluc",
    "hepg2",
    "helaqc",
    "tplantibodies",
    "sbrodae",
    "snakevenoms",
    "woundfluids",
    "immuno",
    "herceptin",
]

# Visualization settings
N_SPECTRA_DEFAULT = 3
SEED = 42


def get_spectrum_stats(spectra: list) -> dict:
    """Get intensity and m/z statistics from all spectra.

    Args:
        spectra: List of matchms Spectrum objects

    Returns:
        Dictionary with intensity and m/z statistics
    """
    all_intensities = []
    all_mz = []
    max_intensities = []
    peak_counts = []

    for s in spectra:
        intensities = s.peaks.intensities
        mz = s.peaks.mz
        all_intensities.extend(intensities)
        all_mz.extend(mz)
        max_intensities.append(max(intensities) if len(intensities) > 0 else 0)
        peak_counts.append(len(intensities))

    return {
        "max_intensity": max(max_intensities) if max_intensities else 0,
        "mean_max_intensity": np.mean(max_intensities) if max_intensities else 0,
        "mean_peak_count": np.mean(peak_counts) if peak_counts else 0,
        "intensity_range": (min(all_intensities), max(all_intensities))
        if all_intensities
        else (0, 0),
        "n_peaks_total": len(all_intensities),
        "mz_range": (min(all_mz), max(all_mz)) if all_mz else (0, 0),
        "mz_mean": np.mean(all_mz) if all_mz else 0,
        "mz_median": np.median(all_mz) if all_mz else 0,
    }


def check_peak_preservation(
    orig_spectra: list, decoy_spectra: list, n: int = 50
) -> tuple[float, float]:
    """Check what fraction of original peaks are preserved in decoys.

    Args:
        orig_spectra: Original spectra list
        decoy_spectra: Decoy spectra list
        n: Number of spectra to check

    Returns:
        Tuple of (mean_overlap_pct, std_overlap_pct)
    """
    overlaps = []
    for orig, dec in zip(orig_spectra[:n], decoy_spectra[:n]):
        orig_mz = set(np.round(orig.peaks.mz, 3))
        dec_mz = set(np.round(dec.peaks.mz, 3))
        overlap = len(orig_mz & dec_mz) / len(orig_mz) * 100 if len(orig_mz) > 0 else 0
        overlaps.append(overlap)
    return float(np.mean(overlaps)), float(np.std(overlaps))


def load_spectra_for_origin(
    origin: str, mgf_dir: Path, kept_rates: list[float]
) -> tuple[list, dict[float, list]]:
    """Load original and decoy spectra for an origin.

    Args:
        origin: Origin name (e.g., 'hepg2')
        mgf_dir: Directory containing MGF files
        kept_rates: List of retention rates to load

    Returns:
        Tuple of (original_spectra, {rate: decoy_spectra})
    """
    # Load original spectra
    original_file = mgf_dir / f"{origin}.mgf"
    if not original_file.exists():
        raise FileNotFoundError(f"Original MGF not found: {original_file}")

    original_spectra = list(load_from_mgf(str(original_file)))
    print(f"Loaded {len(original_spectra)} original spectra from {original_file.name}")

    # Load decoy spectra at each rate
    decoy_spectra = {}
    for rate in kept_rates:
        decoy_file = mgf_dir / f"{origin}_decoy_{rate:.2f}.mgf"
        if decoy_file.exists():
            decoy_spectra[rate] = list(load_from_mgf(str(decoy_file)))
            print(f"Loaded {len(decoy_spectra[rate])} decoys from {decoy_file.name}")
        else:
            print(f"WARNING: {decoy_file.name} not found!")

    return original_spectra, decoy_spectra


def print_summary_statistics(
    original_spectra: list, decoy_spectra: dict[float, list]
) -> dict:
    """Print summary statistics comparing original and decoy spectra.

    Args:
        original_spectra: List of original spectra
        decoy_spectra: Dict mapping retention rate to decoy spectra list

    Returns:
        Dictionary with all statistics for CSV output
    """
    print("\n" + "=" * 80)
    print(f"SUMMARY STATISTICS (all {len(original_spectra):,} spectra)")
    print("=" * 80)

    # Original stats
    orig_stats = get_spectrum_stats(original_spectra)
    print("\nOriginal spectra:")
    print(f"  Max intensity:      {orig_stats['max_intensity']:.4f}")
    print(f"  Mean max intensity: {orig_stats['mean_max_intensity']:.4f}")
    print(
        f"  Intensity range:    [{orig_stats['intensity_range'][0]:.4f}, {orig_stats['intensity_range'][1]:.4f}]"
    )
    print(f"  Mean peak count:    {orig_stats['mean_peak_count']:.1f}")
    print(
        f"  m/z range:          [{orig_stats['mz_range'][0]:.2f}, {orig_stats['mz_range'][1]:.2f}]"
    )
    print(f"  m/z mean:           {orig_stats['mz_mean']:.2f}")
    print(f"  m/z median:         {orig_stats['mz_median']:.2f}")

    results: dict[str | float, dict] = {"original": orig_stats}

    # Decoy stats
    for rate in sorted(decoy_spectra.keys()):
        masked_pct = int((1 - rate) * 100)
        kept_pct = int(rate * 100)
        stats = get_spectrum_stats(decoy_spectra[rate])
        print(f"\nDecoy {rate:.2f} ({masked_pct}% masked, {kept_pct}% kept):")
        print(f"  Max intensity:      {stats['max_intensity']:.4f}")
        print(f"  Mean max intensity: {stats['mean_max_intensity']:.4f}")
        print(
            f"  Intensity range:    [{stats['intensity_range'][0]:.4f}, {stats['intensity_range'][1]:.4f}]"
        )
        print(f"  Mean peak count:    {stats['mean_peak_count']:.1f}")
        print(
            f"  m/z range:          [{stats['mz_range'][0]:.2f}, {stats['mz_range'][1]:.2f}]"
        )
        print(f"  m/z mean:           {stats['mz_mean']:.2f}")
        print(f"  m/z median:         {stats['mz_median']:.2f}")
        results[rate] = stats

    return results


def print_peak_preservation(
    original_spectra: list, decoy_spectra: dict[float, list], n: int = 50
) -> dict:
    """Print peak preservation analysis.

    Args:
        original_spectra: List of original spectra
        decoy_spectra: Dict mapping retention rate to decoy spectra list
        n: Number of spectra to check

    Returns:
        Dictionary with preservation statistics
    """
    print("\n" + "=" * 80)
    print(f"PEAK PRESERVATION CHECK (first {n} spectra)")
    print("=" * 80)

    results = {}
    for rate in sorted(decoy_spectra.keys()):
        mean_overlap, std_overlap = check_peak_preservation(
            original_spectra, decoy_spectra[rate], n
        )
        expected = rate * 100
        masked_pct = int((1 - rate) * 100)
        print(
            f"decoy_{rate:.2f} ({masked_pct}% masked): "
            f"{mean_overlap:.1f}% ± {std_overlap:.1f}% peaks preserved "
            f"(expected ~{expected:.0f}%)"
        )
        results[rate] = {"mean": mean_overlap, "std": std_overlap, "expected": expected}

    return results


def plot_spectrum_comparison(
    original_spectra: list,
    decoy_spectra: dict[float, list],
    sample_indices: np.ndarray,
    output_path: Path,
    kept_rates: list[float],
) -> None:
    """Plot original vs decoy spectra comparison.

    Args:
        original_spectra: List of original spectra
        decoy_spectra: Dict mapping retention rate to decoy spectra list
        sample_indices: Indices of spectra to plot
        output_path: Path to save the figure
        kept_rates: List of retention rates to plot
    """
    n_spectra = len(sample_indices)
    n_cols = len(kept_rates) + 1

    fig, axes = plt.subplots(n_spectra, n_cols, figsize=(4 * n_cols, 3 * n_spectra))

    # Handle single spectrum case
    if n_spectra == 1:
        axes = axes.reshape(1, -1)

    for row, idx in enumerate(sample_indices):
        # Original spectrum
        orig = original_spectra[idx]
        ax = axes[row, 0]
        ax.stem(
            orig.peaks.mz,
            orig.peaks.intensities,
            markerfmt=" ",
            linefmt="b-",
            basefmt=" ",
        )
        ax.set_title(f"Original (idx={idx})\nmax={max(orig.peaks.intensities):.4f}")
        ax.set_xlabel("m/z")
        ax.set_ylabel("Intensity (raw)")

        # Decoy spectra at each rate
        for col, rate in enumerate(kept_rates, start=1):
            if rate in decoy_spectra:
                dec = decoy_spectra[rate][idx]
                ax = axes[row, col]
                ax.stem(
                    dec.peaks.mz,
                    dec.peaks.intensities,
                    markerfmt=" ",
                    linefmt="r-",
                    basefmt=" ",
                )
                masked_pct = int((1 - rate) * 100)
                kept_pct = int(rate * 100)
                max_int = (
                    max(dec.peaks.intensities) if len(dec.peaks.intensities) > 0 else 0
                )
                ax.set_title(
                    f"{masked_pct}% masked ({kept_pct}% kept)\nmax={max_int:.4f}"
                )
                ax.set_xlabel("m/z")

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved spectrum comparison to {output_path}")


def plot_intensity_distributions(
    original_spectra: list,
    decoy_spectra: dict[float, list],
    output_path: Path,
    kept_rates: list[float],
) -> None:
    """Plot intensity distribution histograms using all available spectra.

    Args:
        original_spectra: List of original spectra
        decoy_spectra: Dict mapping retention rate to decoy spectra list
        output_path: Path to save the figure
        kept_rates: List of retention rates to plot
    """
    n_cols = len(kept_rates) + 1
    fig, axes = plt.subplots(1, n_cols, figsize=(4 * n_cols, 4))

    # Collect all original intensities
    orig_intensities = []
    for s in original_spectra:
        orig_intensities.extend(s.peaks.intensities)

    axes[0].hist(orig_intensities, bins=50, alpha=0.7, color="blue", edgecolor="black")
    axes[0].set_title(
        f"Original\nn={len(orig_intensities):,} peaks\n({len(original_spectra):,} spectra)"
    )
    axes[0].set_xlabel("Intensity")
    axes[0].set_ylabel("Count")
    axes[0].set_yscale("log")

    # Decoy intensities (all spectra)
    for col, rate in enumerate(kept_rates, start=1):
        if rate in decoy_spectra:
            dec_intensities = []
            for s in decoy_spectra[rate]:
                dec_intensities.extend(s.peaks.intensities)

            masked_pct = int((1 - rate) * 100)
            axes[col].hist(
                dec_intensities, bins=50, alpha=0.7, color="red", edgecolor="black"
            )
            axes[col].set_title(
                f"{masked_pct}% masked\nn={len(dec_intensities):,} peaks\n({len(decoy_spectra[rate]):,} spectra)"
            )
            axes[col].set_xlabel("Intensity")
            axes[col].set_yscale("log")

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved intensity distribution to {output_path}")


def plot_mz_distributions(
    original_spectra: list,
    decoy_spectra: dict[float, list],
    output_path: Path,
    kept_rates: list[float],
) -> None:
    """Plot m/z distribution histograms using all available spectra.

    Args:
        original_spectra: List of original spectra
        decoy_spectra: Dict mapping retention rate to decoy spectra list
        output_path: Path to save the figure
        kept_rates: List of retention rates to plot
    """
    n_cols = len(kept_rates) + 1
    fig, axes = plt.subplots(1, n_cols, figsize=(4 * n_cols, 4))

    # Collect all original m/z values
    orig_mz = []
    for s in original_spectra:
        orig_mz.extend(s.peaks.mz)

    axes[0].hist(orig_mz, bins=50, alpha=0.7, color="blue", edgecolor="black")
    axes[0].set_title(
        f"Original\nn={len(orig_mz):,} peaks\n({len(original_spectra):,} spectra)"
    )
    axes[0].set_xlabel("m/z")
    axes[0].set_ylabel("Count")
    axes[0].set_yscale("log")

    # Decoy m/z values (all spectra)
    for col, rate in enumerate(kept_rates, start=1):
        if rate in decoy_spectra:
            dec_mz = []
            for s in decoy_spectra[rate]:
                dec_mz.extend(s.peaks.mz)

            masked_pct = int((1 - rate) * 100)
            axes[col].hist(dec_mz, bins=50, alpha=0.7, color="red", edgecolor="black")
            axes[col].set_title(
                f"{masked_pct}% masked\nn={len(dec_mz):,} peaks\n({len(decoy_spectra[rate]):,} spectra)"
            )
            axes[col].set_xlabel("m/z")
            axes[col].set_yscale("log")

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved m/z distribution to {output_path}")


def generate_summary_csv(
    origin: str,
    summary_stats: dict,
    preservation_stats: dict,
    output_path: Path,
) -> None:
    """Generate CSV summary of verification results.

    Args:
        origin: Origin name
        summary_stats: Dictionary with intensity and m/z statistics
        preservation_stats: Dictionary with peak preservation statistics
        output_path: Path to save the CSV
    """
    lines = [
        "origin,type,max_intensity,mean_max_intensity,intensity_min,intensity_max,"
        "mean_peak_count,n_peaks_total,mz_min,mz_max,mz_mean,mz_median,"
        "preservation_mean,preservation_std,preservation_expected"
    ]

    # Original row
    orig = summary_stats["original"]
    lines.append(
        f"{origin},original,{orig['max_intensity']:.6f},{orig['mean_max_intensity']:.6f},"
        f"{orig['intensity_range'][0]:.6f},{orig['intensity_range'][1]:.6f},"
        f"{orig['mean_peak_count']:.1f},{orig['n_peaks_total']},"
        f"{orig['mz_range'][0]:.2f},{orig['mz_range'][1]:.2f},"
        f"{orig['mz_mean']:.2f},{orig['mz_median']:.2f},,,"
    )

    # Decoy rows
    for rate, stats in summary_stats.items():
        if rate == "original":
            continue
        pres = preservation_stats.get(rate, {})
        lines.append(
            f"{origin},decoy_{rate:.2f},{stats['max_intensity']:.6f},{stats['mean_max_intensity']:.6f},"
            f"{stats['intensity_range'][0]:.6f},{stats['intensity_range'][1]:.6f},"
            f"{stats['mean_peak_count']:.1f},{stats['n_peaks_total']},"
            f"{stats['mz_range'][0]:.2f},{stats['mz_range'][1]:.2f},"
            f"{stats['mz_mean']:.2f},{stats['mz_median']:.2f},"
            f"{pres.get('mean', ''):.1f},{pres.get('std', ''):.1f},{pres.get('expected', ''):.0f}"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines))
    print(f"Saved summary CSV to {output_path}")


def verify_origin(
    origin: str,
    mgf_dir: Path,
    fig_output_dir: Path,
    csv_output_dir: Path,
    kept_rates: list[float],
    n_spectra: int = 3,
    save_csv: bool = True,
) -> None:
    """Run full verification for a single origin.

    Args:
        origin: Origin name
        mgf_dir: Directory containing MGF files
        fig_output_dir: Directory to save figures
        csv_output_dir: Directory to save CSV files
        kept_rates: List of retention rates to analyze
        n_spectra: Number of spectra to plot
        save_csv: Whether to save CSV summary
    """
    print(f"\n{'='*80}")
    print(f"VERIFYING DECOY GENERATION FOR: {origin.upper()}")
    print(f"{'='*80}")

    # Load spectra
    original_spectra, decoy_spectra = load_spectra_for_origin(
        origin, mgf_dir, kept_rates
    )

    if not decoy_spectra:
        print(f"No decoy files found for {origin}, skipping...")
        return

    # Sample random spectra for plotting
    np.random.seed(SEED)
    sample_indices = np.random.choice(
        len(original_spectra), size=min(n_spectra, len(original_spectra)), replace=False
    )
    print(f"\nSampling spectra at indices: {sample_indices}")

    # Print statistics
    summary_stats = print_summary_statistics(original_spectra, decoy_spectra)
    preservation_stats = print_peak_preservation(original_spectra, decoy_spectra)

    # Generate plots
    print("\n" + "=" * 80)
    print("GENERATING PLOTS")
    print("=" * 80)

    plot_spectrum_comparison(
        original_spectra,
        decoy_spectra,
        sample_indices,
        fig_output_dir / f"decoy_verification_{origin}.png",
        kept_rates,
    )

    plot_intensity_distributions(
        original_spectra,
        decoy_spectra,
        fig_output_dir / f"intensity_distribution_{origin}.png",
        kept_rates,
    )

    plot_mz_distributions(
        original_spectra,
        decoy_spectra,
        fig_output_dir / f"mz_distribution_{origin}.png",
        kept_rates,
    )

    # Save CSV summary
    if save_csv:
        generate_summary_csv(
            origin,
            summary_stats,
            preservation_stats,
            csv_output_dir / f"decoy_verification_{origin}.csv",
        )


def main():
    parser = argparse.ArgumentParser(
        description="Verify decoy spectrum generation quality"
    )
    parser.add_argument(
        "--origin",
        type=str,
        default="gluc",
        help=f"Origin to analyze (default: gluc). Available: {', '.join(ORIGINS)}",
    )
    parser.add_argument(
        "--all-origins",
        action="store_true",
        help="Analyze all available origins",
    )
    parser.add_argument(
        "--mgf-dir",
        type=Path,
        default=DEFAULT_MGF_DIR,
        help="Directory containing MGF files",
    )
    parser.add_argument(
        "--fig-output-dir",
        type=Path,
        default=DEFAULT_FIG_OUTPUT_DIR,
        help="Directory to save output files",
    )
    parser.add_argument(
        "--csv-output-dir",
        type=Path,
        default=DEFAULT_CSV_OUTPUT_DIR,
        help="Directory to save CSV files",
    )
    parser.add_argument(
        "--kept-rates",
        type=float,
        nargs="+",
        default=DEFAULT_KEPT_RATES,
        help="Retention rates to analyze (fraction of peaks KEPT). "
        "E.g., 0.20 = 20%% kept = 80%% masked (default: 0.20 0.50 0.80)",
    )
    parser.add_argument(
        "--n-spectra",
        type=int,
        default=N_SPECTRA_DEFAULT,
        help="Number of spectra to plot (default: 3)",
    )
    parser.add_argument(
        "--no-csv",
        action="store_true",
        help="Skip CSV summary generation",
    )
    args = parser.parse_args()

    # Determine which origins to process
    if args.all_origins:
        origins_to_process = ORIGINS
    else:
        origins_to_process = [args.origin]

    print(f"MGF directory: {args.mgf_dir}")
    print(f"Fig output directory: {args.fig_output_dir}")
    print(f"CSV output directory: {args.csv_output_dir}")
    print(f"Retention rates: {args.kept_rates}")
    print(f"Origins to process: {origins_to_process}")

    # Process each origin
    for origin in origins_to_process:
        verify_origin(
            origin=origin,
            mgf_dir=args.mgf_dir,
            fig_output_dir=args.fig_output_dir,
            csv_output_dir=args.csv_output_dir,
            kept_rates=args.kept_rates,
            n_spectra=args.n_spectra,
            save_csv=not args.no_csv,
        )

    print("\n" + "=" * 80)
    print("VERIFICATION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
