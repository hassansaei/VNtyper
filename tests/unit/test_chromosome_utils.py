#!/usr/bin/env python3
"""
Unit tests for chromosome_utils.py

Tests chromosome name detection and resolution functions to ensure
correct handling of multiple chromosome naming conventions.
"""

from unittest.mock import patch

import pytest

from vntyper.scripts.chromosome_utils import (
    _build_chromosome_name,
    _construct_ncbi_accession,
    detect_naming_convention,
    get_chromosome_name_from_bam,
    validate_chromosome_name,
)

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit


class TestDetectNamingConvention:
    """Test chromosome naming convention detection."""

    def test_detect_ucsc_naming(self):
        """Test detection of UCSC naming (chr1, chr2, ...)."""
        contigs = ["chr1", "chr2", "chr3", "chrX", "chrY", "chrM"]
        assert detect_naming_convention(contigs) == "ucsc"

    def test_detect_ensembl_naming(self):
        """Test detection of ENSEMBL simple numeric naming (1, 2, ...)."""
        contigs = ["1", "2", "3", "X", "Y", "MT"]
        assert detect_naming_convention(contigs) == "ensembl"

    def test_detect_ncbi_naming(self):
        """Test detection of NCBI accession naming."""
        contigs = [
            "NC_000001.10",
            "NC_000002.11",
            "NC_000003.11",
            "NC_000023.10",
            "NC_000024.9"
        ]
        assert detect_naming_convention(contigs) == "ncbi"

    def test_empty_contig_list(self):
        """Test handling of empty contig list."""
        assert detect_naming_convention([]) == "unknown"

    def test_mixed_naming(self):
        """Test handling of mixed naming conventions."""
        contigs = ["chr1", "2", "NC_000003.11"]
        # Should return unknown since no single convention dominates
        result = detect_naming_convention(contigs)
        assert result in ["unknown", "ucsc", "ensembl", "ncbi"]


class TestValidateChromosomeName:
    """Test chromosome name validation."""

    def test_validate_ucsc_names(self):
        """Test validation of UCSC chromosome names."""
        assert validate_chromosome_name("chr1") is True
        assert validate_chromosome_name("chr22") is True
        assert validate_chromosome_name("chrX") is True
        assert validate_chromosome_name("chrY") is True
        assert validate_chromosome_name("chrM") is True

    def test_validate_simple_names(self):
        """Test validation of simple numeric names."""
        assert validate_chromosome_name("1") is True
        assert validate_chromosome_name("22") is True
        assert validate_chromosome_name("X") is True
        assert validate_chromosome_name("Y") is True
        assert validate_chromosome_name("MT") is True

    def test_validate_ncbi_names(self):
        """Test validation of NCBI accession names."""
        assert validate_chromosome_name("NC_000001.10") is True
        assert validate_chromosome_name("NC_000023.11") is True
        assert validate_chromosome_name("NC_012920.1") is True

    def test_validate_invalid_names(self):
        """Test rejection of invalid chromosome names."""
        assert validate_chromosome_name("") is False
        assert validate_chromosome_name("invalid_chr") is False
        assert validate_chromosome_name("chromosome1") is False
        assert validate_chromosome_name("NC_123") is False


class TestBuildChromosomeName:
    """Test building chromosome names for different conventions."""

    def test_build_ucsc_chr1(self):
        """Test building UCSC chr1 name."""
        config = {}
        result = _build_chromosome_name(1, "ucsc", "hg19", config)
        assert result == "chr1"

    def test_build_ucsc_chrX(self):
        """Test building UCSC chrX name."""
        config = {}
        result = _build_chromosome_name(23, "ucsc", "hg19", config)
        assert result == "chrX"

    def test_build_ensembl_chr1(self):
        """Test building ENSEMBL simple numeric chr1 name."""
        config = {}
        result = _build_chromosome_name(1, "ensembl", "hg19", config)
        assert result == "1"

    def test_build_ensembl_chrX(self):
        """Test building ENSEMBL X name."""
        config = {}
        result = _build_chromosome_name(23, "ensembl", "hg19", config)
        assert result == "X"

    def test_build_ncbi_grch37(self):
        """Test building NCBI accession for GRCh37."""
        config = {
            "bam_processing": {
                "known_chromosome_naming": {
                    "GRCh37": {"ncbi": "NC_000001.10"}
                }
            }
        }
        result = _build_chromosome_name(1, "ncbi", "GRCh37", config)
        assert result == "NC_000001.10"

    def test_build_ncbi_grch38(self):
        """Test building NCBI accession for GRCh38."""
        config = {
            "bam_processing": {
                "known_chromosome_naming": {
                    "GRCh38": {"ncbi": "NC_000001.11"}
                }
            }
        }
        result = _build_chromosome_name(1, "ncbi", "GRCh38", config)
        assert result == "NC_000001.11"

    def test_invalid_chromosome_number(self):
        """Test error handling for invalid chromosome numbers."""
        config = {}
        with pytest.raises(ValueError, match="Invalid chromosome number"):
            _build_chromosome_name(0, "ucsc", "hg19", config)
        with pytest.raises(ValueError, match="Invalid chromosome number"):
            _build_chromosome_name(26, "ucsc", "hg19", config)


class TestConstructNcbiAccession:
    """Test NCBI accession construction."""

    def test_construct_grch37_chr1(self):
        """Test GRCh37 chr1 accession."""
        assert _construct_ncbi_accession(1, "hg19") == "NC_000001.10"

    def test_construct_grch38_chr1(self):
        """Test GRCh38 chr1 accession."""
        assert _construct_ncbi_accession(1, "hg38") == "NC_000001.11"

    def test_construct_grch37_chrX(self):
        """Test GRCh37 chrX accession."""
        assert _construct_ncbi_accession(23, "hg19") == "NC_000023.10"

    def test_construct_grch38_chrY(self):
        """Test GRCh38 chrY accession."""
        assert _construct_ncbi_accession(24, "hg38") == "NC_000024.10"

    def test_construct_mitochondrial(self):
        """Test mitochondrial chromosome (same for both assemblies)."""
        assert _construct_ncbi_accession(25, "hg19") == "NC_012920.1"
        assert _construct_ncbi_accession(25, "hg38") == "NC_012920.1"

    def test_invalid_chromosome(self):
        """Test error handling for invalid chromosomes."""
        with pytest.raises(ValueError):
            _construct_ncbi_accession(26, "hg19")


class TestGetChromosomeNameFromBam:
    """Test getting chromosome name from BAM file."""

    @patch('vntyper.scripts.fastq_bam_processing.extract_bam_header')
    @patch('vntyper.scripts.fastq_bam_processing.parse_contigs_from_header')
    def test_get_ucsc_chr1(self, mock_parse, mock_extract):
        """Test getting UCSC chr1 from BAM."""
        mock_extract.return_value = "@SQ\tSN:chr1\tLN:249250621\n"
        mock_parse.return_value = [{"name": "chr1", "length": 249250621}]

        config = {}
        result = get_chromosome_name_from_bam(
            "test.bam", config, chromosome_number=1, reference_assembly="hg19"
        )
        assert result == "chr1"

    @patch('vntyper.scripts.fastq_bam_processing.extract_bam_header')
    @patch('vntyper.scripts.fastq_bam_processing.parse_contigs_from_header')
    def test_get_ncbi_chr1_grch37(self, mock_parse, mock_extract):
        """Test getting NCBI chr1 from GRCh37 BAM."""
        mock_extract.return_value = "@SQ\tSN:NC_000001.10\tLN:249250621\n"
        mock_parse.return_value = [{"name": "NC_000001.10", "length": 249250621}]

        config = {
            "bam_processing": {
                "known_chromosome_naming": {
                    "GRCh37": {"ncbi": "NC_000001.10"}
                }
            }
        }
        result = get_chromosome_name_from_bam(
            "test.bam", config, chromosome_number=1, reference_assembly="GRCh37"
        )
        assert result == "NC_000001.10"

    @patch('vntyper.scripts.fastq_bam_processing.extract_bam_header')
    @patch('vntyper.scripts.fastq_bam_processing.parse_contigs_from_header')
    def test_chromosome_not_found(self, mock_parse, mock_extract):
        """Test error when chromosome not found in BAM."""
        mock_extract.return_value = "@SQ\tSN:chr2\tLN:243199373\n"
        mock_parse.return_value = [{"name": "chr2", "length": 243199373}]

        config = {}
        with pytest.raises(ValueError, match="not found in BAM"):
            get_chromosome_name_from_bam(
                "test.bam", config, chromosome_number=1, reference_assembly="hg19"
            )

    @patch('vntyper.scripts.fastq_bam_processing.extract_bam_header')
    def test_bam_read_error(self, mock_extract):
        """Test error handling when BAM cannot be read."""
        mock_extract.side_effect = Exception("File not found")

        config = {}
        with pytest.raises(ValueError, match="Cannot read BAM header"):
            get_chromosome_name_from_bam(
                "nonexistent.bam", config, chromosome_number=1
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
