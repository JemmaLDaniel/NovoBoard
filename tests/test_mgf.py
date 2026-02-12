"""Tests for the MGF parser module."""

import pytest

from novoboard.mgf import MGFParser, Spectrum


class TestSpectrum:
    """Tests for the Spectrum dataclass."""

    def test_peptide_mass_calculation(self):
        """Test peptide mass calculation from precursor m/z and charge."""
        spectrum = Spectrum(
            precursor_mz=500.0,
            charge=2,
        )
        # Expected: 500.0 * 2 - 2 * 1.0078 = 997.9844
        assert abs(spectrum.peptide_mass - 997.9844) < 0.01

    def test_feature_id_generation(self):
        """Test feature ID generation from title and scan."""
        spectrum = Spectrum(
            title="test.mgf",
            scan="123",
        )
        assert spectrum.feature_id == "test.mgf||123"


class TestMGFParser:
    """Tests for the MGFParser class."""

    @pytest.fixture
    def sample_mgf_file(self, tmp_path):
        """Create a sample MGF file for testing."""
        mgf_content = """BEGIN IONS
TITLE=test.mgf
PEPMASS=500.25 1000.0
CHARGE=2+
SCANS=123
100.1 500.0
200.2 300.0
300.3 100.0
END IONS

BEGIN IONS
TITLE=test.mgf
PEPMASS=600.35
CHARGE=3+
SCANS=456
150.1 600.0
250.2 400.0
END IONS
"""
        mgf_file = tmp_path / "test.mgf"
        mgf_file.write_text(mgf_content)
        return mgf_file

    def test_iteration(self, sample_mgf_file):
        """Test that parser can iterate over spectra."""
        parser = MGFParser(sample_mgf_file)
        spectra = list(parser)
        assert len(spectra) == 2

    def test_spectrum_attributes(self, sample_mgf_file):
        """Test that spectrum attributes are parsed correctly."""
        parser = MGFParser(sample_mgf_file)
        spectra = list(parser)

        spectrum1 = spectra[0]
        assert spectrum1.scan == "123"
        assert spectrum1.precursor_mz == 500.25
        assert spectrum1.charge == 2
        assert len(spectrum1.peaks) == 3

    def test_peaks_parsing(self, sample_mgf_file):
        """Test that peaks are parsed correctly."""
        parser = MGFParser(sample_mgf_file)
        spectra = list(parser)

        peaks = spectra[0].peaks
        assert peaks[0] == (100.1, 500.0)
        assert peaks[1] == (200.2, 300.0)
        assert peaks[2] == (300.3, 100.0)

    def test_to_list(self, sample_mgf_file):
        """Test to_list method."""
        parser = MGFParser(sample_mgf_file)
        spectra = parser.to_list()
        assert isinstance(spectra, list)
        assert len(spectra) == 2

    def test_to_dict(self, sample_mgf_file):
        """Test to_dict method."""
        parser = MGFParser(sample_mgf_file)
        spectra_dict = parser.to_dict()
        assert isinstance(spectra_dict, dict)
        assert "test.mgf||123" in spectra_dict
        assert "test.mgf||456" in spectra_dict

    def test_count_spectra(self, sample_mgf_file):
        """Test count_spectra method."""
        parser = MGFParser(sample_mgf_file)
        count = parser.count_spectra()
        assert count == 2

    def test_empty_file(self, tmp_path):
        """Test parsing an empty file."""
        empty_file = tmp_path / "empty.mgf"
        empty_file.write_text("")
        parser = MGFParser(empty_file)
        spectra = list(parser)
        assert len(spectra) == 0
