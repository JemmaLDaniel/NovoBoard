"""MGF file parser for mass spectrometry data."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

from novoboard import config

# Precompiled regex patterns
_SPACE_PATTERN = re.compile(r" |\t|\r|\n")
_FIELD_PATTERN = re.compile(r"[=\r\n]")
_CHARGE_PATTERN = re.compile(r"=|\+|\r|\n")


@dataclass
class Spectrum:
    """A single MS/MS spectrum from an MGF file.

    Attributes:
        title: Spectrum title from the file
        scan: Scan number
        precursor_mz: Precursor m/z value
        charge: Precursor charge state
        peaks: List of (m/z, intensity) tuples
        peptide_mass: Calculated peptide mass
    """

    title: str = ""
    scan: str = ""
    precursor_mz: float = 0.0
    charge: int = 0
    peaks: list[tuple[float, float]] = field(default_factory=list)

    @property
    def peptide_mass(self) -> float:
        """Calculate peptide mass from precursor m/z and charge."""
        return self.precursor_mz * self.charge - self.charge * config.mass_H

    @property
    def feature_id(self) -> str:
        """Generate feature ID from title and scan."""
        return f"{self.title}||{self.scan}"


class MGFParser:
    """Parser for MGF (Mascot Generic Format) files.

    Example usage:
        parser = MGFParser(Path("spectra.mgf"))
        for spectrum in parser:
            print(spectrum.scan, len(spectrum.peaks))
    """

    def __init__(self, filepath: Path | str) -> None:
        """Initialize the parser.

        Args:
            filepath: Path to the MGF file
        """
        self.filepath = Path(filepath)

    def __iter__(self) -> Iterator[Spectrum]:
        """Iterate over spectra in the file.

        Yields:
            Spectrum objects parsed from the file
        """
        with open(self.filepath, "r") as f:
            current_spectrum: Spectrum | None = None

            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line == "BEGIN IONS":
                    current_spectrum = Spectrum()
                elif line == "END IONS":
                    if current_spectrum is not None:
                        yield current_spectrum
                        current_spectrum = None
                elif current_spectrum is not None:
                    if "=" in line:
                        # Parse header field
                        key, value = line.split("=", 1)
                        if key == "TITLE":
                            # Extract source file from title
                            current_spectrum.title = self._extract_source_file(value)
                        elif key == "SCANS":
                            current_spectrum.scan = value
                        elif key == "PEPMASS":
                            # PEPMASS may have intensity after space
                            parts = _SPACE_PATTERN.split(value)
                            current_spectrum.precursor_mz = float(parts[0])
                        elif key == "CHARGE":
                            # Remove trailing + or -
                            charge_str = value.rstrip("+-")
                            current_spectrum.charge = int(charge_str)
                    elif line[0].isdigit():
                        # Parse peak line
                        parts = _SPACE_PATTERN.split(line)
                        if len(parts) >= 2:
                            mz = float(parts[0])
                            intensity = float(parts[1])
                            current_spectrum.peaks.append((mz, intensity))

    def _extract_source_file(self, title: str) -> str:
        """Extract source file name from title.

        Args:
            title: Full title string from MGF

        Returns:
            Cleaned source file name ending with .mgf
        """
        # Handle various title formats
        # e.g., "File: path\\file.raw Scan: 123" -> "file.mgf"
        if "File:" in title:
            parts = title.split("File:")[1].strip()
            filename = parts.split()[0].split("\\")[-1]
            if ".raw" in filename:
                return filename.split(".raw")[0] + ".mgf"
            return filename
        return title.split("||")[0] if "||" in title else title

    def to_list(self) -> list[Spectrum]:
        """Load all spectra into memory.

        Returns:
            List of all Spectrum objects in the file
        """
        return list(self)

    def to_dict(self) -> dict[str, list[tuple[float, float]]]:
        """Load spectra as a dictionary keyed by feature_id.

        Returns:
            Dictionary mapping feature_id -> peaks list
        """
        return {spectrum.feature_id: spectrum.peaks for spectrum in self}

    def count_spectra(self) -> int:
        """Count spectra in the file without loading all into memory.

        Returns:
            Number of spectra in the file
        """
        count = 0
        with open(self.filepath, "r") as f:
            for line in f:
                if line.strip() == "END IONS":
                    count += 1
        return count
