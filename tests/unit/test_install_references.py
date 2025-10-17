#!/usr/bin/env python3
"""
Unit tests for install_references.py

Tests cover:
- Configuration loading and validation
- File downloading
- MD5 checksum calculation
- Archive extraction
- Aligner detection and indexing
- Error handling
"""

import json
from unittest.mock import Mock, patch

import pytest

from vntyper.scripts.install_references import (
    calculate_md5,
    check_executable_available,
    check_index_exists,
    detect_index_conflicts,
    download_file,
    get_enabled_aligners,
    load_install_config,
)

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def mock_config():
    """Provide a mock configuration dictionary."""
    return {
        "aligners": {
            "bwa": {
                "enabled": True,
                "executable": "bwa",
                "index_files": [".amb", ".ann", ".bwt", ".pac", ".sa"],
                "description": "BWA aligner",
            },
            "minimap2": {
                "enabled": False,
                "executable": "minimap2",
                "index_files": [".mmi"],
                "description": "Minimap2 aligner",
            },
        },
        "ucsc_references": {
            "hg19": {
                "url": "https://example.com/hg19.fa.gz",
                "target_path": "alignment/chr1.hg19.fa.gz",
            }
        },
        "ensembl_references": {
            "hg19_ensembl": {
                "url": "https://ftp.ensembl.org/pub/grch37/chr1.fa.gz",
                "target_path": "alignment/chr1.hg19_ensembl.fa.gz",
                "description": "ENSEMBL GRCh37 chromosome 1",
            }
        },
    }


@pytest.fixture
def temp_config_file(tmp_path, mock_config):
    """Create a temporary config file."""
    config_file = tmp_path / "test_config.json"
    with open(config_file, "w") as f:
        json.dump(mock_config, f)
    return config_file


@pytest.fixture
def temp_dir(tmp_path):
    """Provide a temporary directory for testing."""
    return tmp_path


# =============================================================================
# Test load_install_config()
# =============================================================================


@pytest.mark.unit
class TestLoadInstallConfig:
    """Tests for configuration loading."""

    def test_load_valid_config(self, temp_config_file, mock_config):
        """Test loading a valid configuration file."""
        config = load_install_config(temp_config_file)
        assert config == mock_config
        assert "aligners" in config
        assert "ucsc_references" in config

    def test_missing_config_file(self, tmp_path):
        """Test error when config file doesn't exist."""
        nonexistent = tmp_path / "nonexistent.json"
        with pytest.raises(SystemExit):
            load_install_config(nonexistent)

    def test_invalid_json(self, tmp_path):
        """Test error when config file has invalid JSON."""
        invalid_config = tmp_path / "invalid.json"
        invalid_config.write_text("{invalid json content")

        with pytest.raises(SystemExit):
            load_install_config(invalid_config)

    def test_empty_config(self, tmp_path):
        """Test loading an empty but valid JSON config."""
        empty_config = tmp_path / "empty.json"
        empty_config.write_text("{}")

        config = load_install_config(empty_config)
        assert config == {}


# =============================================================================
# Test calculate_md5()
# =============================================================================


@pytest.mark.unit
class TestCalculateMD5:
    """Tests for MD5 checksum calculation."""

    def test_calculate_md5_simple_file(self, tmp_path):
        """Test MD5 calculation for a simple file."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("Hello, World!")

        md5_hash = calculate_md5(test_file)

        # "Hello, World!" has known MD5
        assert md5_hash == "65a8e27d8879283831b664bd8b7f0ad4"

    def test_calculate_md5_empty_file(self, tmp_path):
        """Test MD5 calculation for an empty file."""
        test_file = tmp_path / "empty.txt"
        test_file.write_text("")

        md5_hash = calculate_md5(test_file)

        # Empty file has known MD5
        assert md5_hash == "d41d8cd98f00b204e9800998ecf8427e"

    def test_calculate_md5_binary_file(self, tmp_path):
        """Test MD5 calculation for a binary file."""
        test_file = tmp_path / "binary.dat"
        test_file.write_bytes(b"\x00\x01\x02\x03\x04")

        md5_hash = calculate_md5(test_file)

        # Binary sequence has deterministic MD5
        assert md5_hash == "d05374dc381d9b52806446a71c8e79b1"

    def test_calculate_md5_large_file(self, tmp_path):
        """Test MD5 calculation for a file larger than chunk size."""
        test_file = tmp_path / "large.txt"
        # Write 10KB of data (larger than 4KB chunk size)
        test_file.write_bytes(b"A" * 10240)

        md5_hash = calculate_md5(test_file)

        # Verify it returns a valid MD5 hash
        assert len(md5_hash) == 32
        assert all(c in "0123456789abcdef" for c in md5_hash)

    def test_calculate_md5_nonexistent_file(self, tmp_path):
        """Test error when file doesn't exist."""
        nonexistent = tmp_path / "nonexistent.txt"

        with pytest.raises(SystemExit):
            calculate_md5(nonexistent)


# =============================================================================
# Test check_executable_available()
# =============================================================================


@pytest.mark.unit
class TestCheckExecutableAvailable:
    """Tests for executable availability checking."""

    @patch("subprocess.run")
    def test_executable_found(self, mock_run):
        """Test when executable is available."""
        mock_run.return_value = Mock(returncode=0, stdout="/usr/bin/bwa\n")

        result = check_executable_available("bwa")

        assert result is True
        mock_run.assert_called_once()

    @patch("subprocess.run")
    def test_executable_not_found(self, mock_run):
        """Test when executable is not available."""
        mock_run.return_value = Mock(returncode=1, stdout="")

        result = check_executable_available("nonexistent_tool")

        assert result is False

    @patch("subprocess.run")
    def test_executable_check_exception(self, mock_run):
        """Test error handling when check fails."""
        mock_run.side_effect = Exception("Command failed")

        result = check_executable_available("bwa")

        assert result is False


# =============================================================================
# Test get_enabled_aligners()
# =============================================================================


@pytest.mark.unit
class TestGetEnabledAligners:
    """Tests for enabled aligner detection."""

    @patch("vntyper.scripts.install_references.check_executable_available")
    def test_get_enabled_aligners_all_available(self, mock_check, mock_config):
        """Test when all enabled aligners are available."""
        mock_check.return_value = True

        enabled = get_enabled_aligners(mock_config["aligners"])

        assert "bwa" in enabled
        assert "minimap2" not in enabled  # Not enabled in config

    @patch("vntyper.scripts.install_references.check_executable_available")
    def test_get_enabled_aligners_missing_executable(self, mock_check, mock_config):
        """Test when enabled aligner executable is missing."""
        mock_check.return_value = False

        enabled = get_enabled_aligners(mock_config["aligners"])

        assert len(enabled) == 0

    @patch("vntyper.scripts.install_references.check_executable_available")
    def test_get_enabled_aligners_mixed_availability(self, mock_check):
        """Test mixed availability of aligners."""
        config = {
            "bwa": {"enabled": True, "executable": "bwa"},
            "minimap2": {"enabled": True, "executable": "minimap2"},
        }

        # BWA available, minimap2 not
        mock_check.side_effect = lambda exe: exe == "bwa"

        enabled = get_enabled_aligners(config)

        assert "bwa" in enabled
        assert "minimap2" not in enabled

    def test_get_enabled_aligners_empty_config(self):
        """Test with empty aligner config."""
        enabled = get_enabled_aligners({})

        assert len(enabled) == 0


# =============================================================================
# Test detect_index_conflicts()
# =============================================================================


@pytest.mark.unit
class TestDetectIndexConflicts:
    """Tests for index file conflict detection."""

    def test_no_conflicts(self):
        """Test when there are no conflicts."""
        aligners = {
            "bwa": {"index_files": [".amb", ".ann", ".bwt"]},
            "minimap2": {"index_files": [".mmi"]},
        }

        conflicts = detect_index_conflicts(aligners)

        assert len(conflicts) == 0

    def test_detect_conflicts(self):
        """Test detection of conflicting index files."""
        aligners = {
            "bwa": {"index_files": [".amb", ".idx"]},
            "bwa-mem2": {"index_files": [".idx", ".pac"]},
        }

        conflicts = detect_index_conflicts(aligners)

        assert len(conflicts) == 1
        assert ".idx" in conflicts[0]

    def test_multiple_conflicts(self):
        """Test detection of multiple conflicts."""
        aligners = {
            "aligner1": {"index_files": [".a", ".b"]},
            "aligner2": {"index_files": [".b", ".c"]},
            "aligner3": {"index_files": [".c", ".d"]},
        }

        conflicts = detect_index_conflicts(aligners)

        # Should detect conflicts for .b and .c
        assert len(conflicts) == 2


# =============================================================================
# Test check_index_exists()
# =============================================================================


@pytest.mark.unit
class TestCheckIndexExists:
    """Tests for index file existence checking."""

    def test_standard_index_exists(self, tmp_path):
        """Test checking standard index files."""
        ref_path = tmp_path / "ref.fa"
        ref_path.write_text(">chr1\nACGT")

        # Create index files
        (tmp_path / "ref.fa.amb").write_text("index")
        (tmp_path / "ref.fa.bwt").write_text("index")

        aligner_info = {"index_files": [".amb", ".bwt"]}

        result = check_index_exists(ref_path, "bwa", aligner_info)

        assert result is True

    def test_standard_index_missing(self, tmp_path):
        """Test when index files are missing."""
        ref_path = tmp_path / "ref.fa"
        ref_path.write_text(">chr1\nACGT")

        aligner_info = {"index_files": [".amb", ".bwt"]}

        result = check_index_exists(ref_path, "bwa", aligner_info)

        assert result is False

    def test_index_base_required(self, tmp_path):
        """Test checking index files with base name."""
        ref_path = tmp_path / "ref.fa"
        ref_path.write_text(">chr1\nACGT")

        # Create index files with base name
        (tmp_path / "ref_bwa.1.bt2").write_text("index")
        (tmp_path / "ref_bwa.2.bt2").write_text("index")

        aligner_info = {
            "requires_index_base": True,
            "index_files": [".1.bt2", ".2.bt2"],
        }

        result = check_index_exists(ref_path, "bwa", aligner_info)

        assert result is True


# =============================================================================
# Test download_file()
# =============================================================================


@pytest.mark.unit
class TestDownloadFile:
    """Tests for file downloading."""

    @patch("vntyper.scripts.install_references.urlretrieve")
    def test_download_success(self, mock_urlretrieve, tmp_path):
        """Test successful file download."""
        dest_path = tmp_path / "download" / "file.txt"
        url = "https://example.com/file.txt"

        download_file(url, dest_path)

        mock_urlretrieve.assert_called_once_with(url, dest_path)
        assert dest_path.parent.exists()

    @patch("vntyper.scripts.install_references.urlretrieve")
    def test_download_file_exists(self, mock_urlretrieve, tmp_path):
        """Test that existing files are not re-downloaded."""
        dest_path = tmp_path / "existing.txt"
        dest_path.write_text("existing content")
        url = "https://example.com/file.txt"

        download_file(url, dest_path)

        # Should not call urlretrieve for existing file
        mock_urlretrieve.assert_not_called()

    @patch("vntyper.scripts.install_references.urlretrieve")
    def test_download_failure(self, mock_urlretrieve, tmp_path):
        """Test download failure handling."""
        mock_urlretrieve.side_effect = Exception("Network error")
        dest_path = tmp_path / "file.txt"
        url = "https://example.com/file.txt"

        with pytest.raises(SystemExit):
            download_file(url, dest_path)


# =============================================================================
# Integration-style Tests
# =============================================================================


@pytest.mark.unit
class TestInstallReferencesIntegration:
    """Integration-style tests for overall functionality."""

    def test_ensembl_references_in_config(self, mock_config):
        """Verify ENSEMBL references are present in config."""
        assert "ensembl_references" in mock_config
        assert "hg19_ensembl" in mock_config["ensembl_references"]

        ensembl_ref = mock_config["ensembl_references"]["hg19_ensembl"]
        assert "url" in ensembl_ref
        assert "target_path" in ensembl_ref
        assert "ftp.ensembl.org" in ensembl_ref["url"]

    def test_config_has_all_reference_types(self, mock_config):
        """Verify config contains all expected reference types."""
        expected_sections = ["ucsc_references", "ensembl_references"]

        for section in expected_sections:
            assert section in mock_config

    @patch("vntyper.scripts.install_references.check_executable_available")
    def test_aligner_configuration(self, mock_check, mock_config):
        """Test aligner configuration and selection."""
        mock_check.return_value = True

        enabled = get_enabled_aligners(mock_config["aligners"])

        # BWA is enabled and available
        assert "bwa" in enabled
        assert enabled["bwa"]["executable"] == "bwa"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
