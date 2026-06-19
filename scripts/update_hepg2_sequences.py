#!/usr/bin/env python3
"""Update SEQ entries in hepg2 MGF files with corrected sequences from parquet."""

import re
from pathlib import Path

import polars as pl


def parse_title(title_line: str) -> tuple[str, int] | None:
    """Extract experiment_name and scan_number from TITLE line.

    Example: TITLE=20190601_QX6_JoMu_SA_uPac200cm_HepG2_f3 scan=93738
    Returns: ("20190601_QX6_JoMu_SA_uPac200cm_HepG2_f3", 93738)
    """
    match = re.match(r"TITLE=(.+)\s+scan=(\d+)", title_line)
    if match:
        return match.group(1), int(match.group(2))
    return None


def update_mgf_sequences(
    mgf_path: Path, lookup: dict[tuple[str, int], str]
) -> tuple[int, int]:
    """Update SEQ entries in an MGF file.

    Returns (updated_count, total_count)
    """
    lines = mgf_path.read_text().splitlines()
    updated_count = 0
    total_count = 0
    current_experiment = None
    current_scan = None

    new_lines = []
    for line in lines:
        if line.startswith("TITLE="):
            parsed = parse_title(line)
            if parsed:
                current_experiment, current_scan = parsed
            new_lines.append(line)
        elif line.startswith("SEQ="):
            total_count += 1
            key = (current_experiment, current_scan)
            if key in lookup:
                new_seq = lookup[key]
                new_lines.append(f"SEQ={new_seq}")
                updated_count += 1
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)

    mgf_path.write_text("\n".join(new_lines) + "\n")
    return updated_count, total_count


def main():
    parquet_path = Path(
        "/home/j-daniel/Documents/winnow/data/external_datasets/hepg2_labelled.parquet"
    )
    mgf_dir = Path(
        "/home/j-daniel/Documents/winnow/novoboard/benchmark_data/instanovo_inputs/per_origin_mgf"
    )

    # Load corrected sequences
    print(f"Loading corrected sequences from {parquet_path}")
    df = pl.read_parquet(parquet_path)

    # Build lookup dict: (experiment_name, scan_number) -> sequence
    lookup = {
        (row["experiment_name"], row["scan_number"]): row["sequence"]
        for row in df.select(["experiment_name", "scan_number", "sequence"]).iter_rows(
            named=True
        )
    }
    print(f"Loaded {len(lookup)} sequences")

    # Find all hepg2 MGF files (including decoys)
    mgf_files = sorted(mgf_dir.glob("hepg2*.mgf"))
    print(f"Found {len(mgf_files)} hepg2 MGF files to process")

    total_updated = 0
    total_spectra = 0

    for mgf_path in mgf_files:
        updated, total = update_mgf_sequences(mgf_path, lookup)
        print(f"  {mgf_path.name}: updated {updated}/{total} sequences")
        total_updated += updated
        total_spectra += total

    print(
        f"\nDone! Updated {total_updated}/{total_spectra} sequences across {len(mgf_files)} files"
    )


if __name__ == "__main__":
    main()
