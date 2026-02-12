"""Plotting functions for FDR validation results."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Sequence

from matplotlib import pyplot

if TYPE_CHECKING:
    from novoboard.fdr import FDRValidationResult

logger = logging.getLogger(__name__)

# Default color palette for plots
DEFAULT_COLORS = ["b", "c", "g", "k", "m", "r", "y", "orange", "pink", "brown"]


def plot_fdr_validation(
    results_list: list[FDRValidationResult],
    samples: Sequence[int] | None = None,
    output_path: str = "fdr_validation.png",
    labels: Sequence[str] | None = None,
    fdr_max: float = 0.05,
    dpi: int = 150,
) -> None:
    """Plot FDR validation results.

    Args:
        results_list: List of FDRValidationResult from validate_FDR
        samples: Sample indices for coloring (optional, defaults to 1..n)
        output_path: Path to save the output figure
        labels: Custom labels for each result (optional)
        fdr_max: Maximum value for FDR axis (default: 0.05)
        dpi: Resolution in dots per inch (default: 150, use 300 for print quality)
    """
    if not results_list:
        logger.warning("No results to plot")
        return

    # Create output directory if needed
    output_dir = Path(output_path).parent
    if output_dir and str(output_dir) != ".":
        output_dir.mkdir(parents=True, exist_ok=True)

    # Default samples to 1..n
    if samples is None:
        samples = range(1, len(results_list) + 1)

    # Generate labels if not provided
    if labels is None:
        labels = [f"Decoy {i}" for i in samples]

    # Assign colors (cycle through palette if needed)
    colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)] for i in range(len(results_list))]

    fig, ax = pyplot.subplots(1, 2, figsize=(9, 4))

    for results, color, label in zip(results_list, colors, labels):
        ax[0].plot(results.estimated_fdr, results.true_fdr_T, color=color, label=label)
        ax[1].plot(results.cumsum, results.estimated_fdr, color=color, label=label)

    # Diagonal reference line
    ax[0].plot(
        [0, fdr_max], [0, fdr_max], color="black", linestyle="--", label="True FDR"
    )
    ax[0].set_xlim(0, fdr_max)
    ax[0].set_ylim(0, fdr_max)
    ax[0].set_xlabel("Estimated FDR")
    ax[0].set_ylabel("True FDR")
    ax[0].legend()

    # Add true FDR line to second plot (use middle result if available)
    mid_idx = len(results_list) // 2
    if results_list:
        ax[1].plot(
            results_list[mid_idx].cumsum,
            results_list[mid_idx].true_fdr_T,
            color="black",
            linestyle="--",
            label="True FDR",
        )
    ax[1].set_ylim(0, fdr_max)
    ax[1].set_xlabel("Number of PSMs")
    ax[1].set_ylabel("FDR")
    ax[1].legend()

    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi)
    logger.info(f"Figure saved to {output_path} (dpi={dpi})")
    pyplot.close(fig)
