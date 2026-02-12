"""Tests for decoy module."""

import os
import tempfile
from pathlib import Path

import pytest

from novoboard.decoy import generate_decoy_mgf


def get_decoy_path(mgf_path: str, rate: float = 0.5) -> str:
    """Get the expected decoy output path for a given MGF file."""
    p = Path(mgf_path)
    return str(p.parent / f"{p.stem}_decoy_{rate:.2f}.mgf")


class TestGenerateDecoyMgf:
    """Test generate_decoy_mgf function."""

    @pytest.fixture
    def sample_mgf_file(self):
        """Create a minimal sample MGF file for testing."""
        content = """BEGIN IONS
TITLE=test_spectrum.raw-1
PEPMASS=500.5
CHARGE=2+
SCANS=1
100.0 1000.0
200.0 2000.0
300.0 3000.0
400.0 4000.0
END IONS

BEGIN IONS
TITLE=test_spectrum.raw-2
PEPMASS=600.5
CHARGE=3+
SCANS=2
150.0 1500.0
250.0 2500.0
350.0 3500.0
450.0 4500.0
END IONS

"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".mgf", delete=False) as f:
            f.write(content)
            temp_path = f.name
        yield temp_path
        # Cleanup
        os.unlink(temp_path)
        decoy_path = get_decoy_path(temp_path)
        if os.path.exists(decoy_path):
            os.unlink(decoy_path)

    def test_generates_decoy_file(self, sample_mgf_file):
        """Test that generate_decoy_mgf creates an output file."""
        generate_decoy_mgf([sample_mgf_file], peak_sampling="random", sampling_rate=0.5)

        decoy_path = get_decoy_path(sample_mgf_file)
        assert os.path.exists(decoy_path)

    def test_decoy_file_has_content(self, sample_mgf_file):
        """Test that the generated decoy file has content."""
        generate_decoy_mgf([sample_mgf_file], peak_sampling="random", sampling_rate=0.5)

        decoy_path = get_decoy_path(sample_mgf_file)
        with open(decoy_path, "r") as f:
            content = f.read()

        assert "BEGIN IONS" in content
        assert "END IONS" in content

    def test_decoy_preserves_spectrum_count(self, sample_mgf_file):
        """Test that decoy file has same number of spectra."""
        generate_decoy_mgf([sample_mgf_file], peak_sampling="random", sampling_rate=0.5)

        decoy_path = get_decoy_path(sample_mgf_file)
        with open(decoy_path, "r") as f:
            content = f.read()

        # Should have 2 spectra (2 BEGIN IONS markers)
        assert content.count("BEGIN IONS") == 2
        assert content.count("END IONS") == 2
