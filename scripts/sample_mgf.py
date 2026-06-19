#!/usr/bin/env python3
"""Sample N spectra from an MGF file and write to a new MGF file.

Reads spectra as raw text blocks (preserving all headers and formatting),
randomly samples N of them, and writes the selected spectra to the output.

Usage:
    python scripts/sample_mgf.py \
        --input /path/to/annotated_train.mgf \
        --output /path/to/annotated_train_sample1000.mgf \
        --n 1000 \
        --seed 42
"""

from __future__ import annotations

import argparse
import logging
import random
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def read_mgf_blocks(filepath: Path) -> list[str]:
    """Read an MGF file and return each spectrum as a raw text block."""
    blocks: list[str] = []
    current_lines: list[str] = []
    in_spectrum = False

    with open(filepath) as f:
        for line in f:
            if line.strip() == "BEGIN IONS":
                in_spectrum = True
                current_lines = [line]
            elif line.strip() == "END IONS":
                current_lines.append(line)
                blocks.append("".join(current_lines))
                in_spectrum = False
                current_lines = []
            elif in_spectrum:
                current_lines.append(line)

    return blocks


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Sample N spectra from an MGF file",
    )
    parser.add_argument("--input", type=Path, required=True, help="Input MGF file")
    parser.add_argument("--output", type=Path, required=True, help="Output MGF file")
    parser.add_argument(
        "--n", type=int, required=True, help="Number of spectra to sample"
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="Random seed (default: 42)"
    )

    args = parser.parse_args()

    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        raise SystemExit(1)

    logger.info(f"Reading spectra from {args.input}")
    blocks = read_mgf_blocks(args.input)
    logger.info(f"Found {len(blocks):,} spectra")

    if args.n >= len(blocks):
        logger.info(f"Requested {args.n} but only {len(blocks)} available — using all")
        sampled = blocks
    else:
        rng = random.Random(args.seed)
        sampled = rng.sample(blocks, args.n)
        logger.info(f"Sampled {len(sampled):,} spectra (seed={args.seed})")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with open(args.output, "w") as f:
        for block in sampled:
            f.write(block)
            f.write("\n")

    logger.info(f"Wrote {len(sampled):,} spectra to {args.output}")


if __name__ == "__main__":
    main()
