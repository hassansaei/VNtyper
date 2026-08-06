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
        contigs = ["NC_000001.10", "NC_000002.11", "NC_000003.11", "NC_000023.10", "NC_000024.9"]
        assert detect_naming_convention(contigs) == "ncbi"

    def test_empty_contig_list(self):
        """Test handling of empty contig list."""
        assert detect_naming_convention([]) == "unknown"

    def test_a_bare_majority_of_ucsc_names_wins(self):
        """3 of 5 clears the 50% threshold, so the convention is decided."""
        contigs = ["chr1", "chr2", "chr3", "chrUn_KI270302v1", "HLA-A*01:01:01:01"]
        assert detect_naming_convention(contigs) == "ucsc"

    def test_a_minority_of_ucsc_names_decides_nothing(self):
        """2 of 5 misses the threshold and no other convention reaches it either."""
        contigs = ["chr1", "chr2", "chrUn_KI270302v1", "HLA-A*01:01:01:01", "phiX174"]
        assert detect_naming_convention(contigs) == "unknown"

    def test_mixed_naming_with_no_majority_is_unknown(self):
        """Three conventions, one contig each: nothing reaches 50%."""
        contigs = ["chr1", "2", "NC_000003.11"]
        assert detect_naming_convention(contigs) == "unknown"


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
        config = {"bam_processing": {"known_chromosome_naming": {"GRCh37": {"ncbi": "NC_000001.10"}}}}
        result = _build_chromosome_name(1, "ncbi", "GRCh37", config)
        assert result == "NC_000001.10"

    def test_build_ncbi_grch38(self):
        """Test building NCBI accession for GRCh38."""
        config = {"bam_processing": {"known_chromosome_naming": {"GRCh38": {"ncbi": "NC_000001.11"}}}}
        result = _build_chromosome_name(1, "ncbi", "GRCh38", config)
        assert result == "NC_000001.11"

    def test_invalid_chromosome_number(self):
        """Test error handling for invalid chromosome numbers."""
        config = {}
        with pytest.raises(ValueError, match="Invalid chromosome number"):
            _build_chromosome_name(0, "ucsc", "hg19", config)
        with pytest.raises(ValueError, match="Invalid chromosome number"):
            _build_chromosome_name(26, "ucsc", "hg19", config)

    @pytest.mark.parametrize(
        "reference_assembly,expected",
        [
            # The real caller normalises to the build name before calling
            # `_construct_ncbi_accession`, so this is the contract that matters.
            ("GRCh37", "NC_000001.10"),
            ("hg19", "NC_000001.10"),
            ("GRCh38", "NC_000001.11"),
            ("hg38", "NC_000001.11"),
        ],
    )
    def test_ncbi_chr1_without_the_config_shortcut_uses_the_declared_build(self, reference_assembly, expected):
        """No `known_chromosome_naming` in config: the accession must still be the declared build's."""
        assert _build_chromosome_name(1, "ncbi", reference_assembly, {}) == expected

    @pytest.mark.parametrize(
        "reference_assembly,chromosome_number,expected",
        [
            ("GRCh37", 5, "NC_000005.9"),
            ("hg19", 5, "NC_000005.9"),
            ("GRCh38", 5, "NC_000005.10"),
            ("hg38", 5, "NC_000005.10"),
        ],
    )
    def test_ncbi_non_chr1_uses_the_declared_build(self, reference_assembly, chromosome_number, expected):
        """The config shortcut only covers chr1; every other chromosome is constructed."""
        config = {"bam_processing": {"known_chromosome_naming": {"GRCh37": {"ncbi": "NC_000001.10"}}}}
        assert _build_chromosome_name(chromosome_number, "ncbi", reference_assembly, config) == expected


class TestConstructNcbiAccession:
    """Test NCBI accession construction."""

    # Real RefSeq accessions, written out independently of the version tables in
    # `chromosome_utils`. Both spellings of each build must select the same table:
    # `_build_chromosome_name` normalises to "GRCh37"/"GRCh38" before calling,
    # while the docstring and the test builders use "hg19"/"hg38".
    GRCH37_ACCESSIONS = {
        1: "NC_000001.10",
        2: "NC_000002.11",
        5: "NC_000005.9",
        7: "NC_000007.13",
        22: "NC_000022.10",
        23: "NC_000023.10",
        24: "NC_000024.9",
        25: "NC_012920.1",
    }
    GRCH38_ACCESSIONS = {
        1: "NC_000001.11",
        2: "NC_000002.12",
        5: "NC_000005.10",
        7: "NC_000007.14",
        22: "NC_000022.11",
        23: "NC_000023.11",
        24: "NC_000024.10",
        25: "NC_012920.1",
    }

    @pytest.mark.parametrize("assembly", ["hg19", "GRCh37"])
    def test_both_spellings_of_grch37_select_the_grch37_table(self, assembly):
        for chromosome_number, expected in self.GRCH37_ACCESSIONS.items():
            assert _construct_ncbi_accession(chromosome_number, assembly) == expected, (
                f"chromosome {chromosome_number} under assembly spelling {assembly!r}"
            )

    @pytest.mark.parametrize("assembly", ["hg38", "GRCh38"])
    def test_both_spellings_of_grch38_select_the_grch38_table(self, assembly):
        for chromosome_number, expected in self.GRCH38_ACCESSIONS.items():
            assert _construct_ncbi_accession(chromosome_number, assembly) == expected, (
                f"chromosome {chromosome_number} under assembly spelling {assembly!r}"
            )

    def test_the_two_builds_do_not_share_an_accession_except_the_mitochondrion(self):
        """Guards the test above from passing because both tables happen to agree."""
        differing = [n for n in self.GRCH37_ACCESSIONS if self.GRCH37_ACCESSIONS[n] != self.GRCH38_ACCESSIONS[n]]
        assert len(differing) == len(self.GRCH37_ACCESSIONS) - 1
        assert self.GRCH37_ACCESSIONS[25] == self.GRCH38_ACCESSIONS[25] == "NC_012920.1"

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

    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    @patch("vntyper.scripts.fastq_bam_processing.parse_contigs_from_header")
    def test_get_ucsc_chr1(self, mock_parse, mock_extract):
        """Test getting UCSC chr1 from BAM."""
        mock_extract.return_value = "@SQ\tSN:chr1\tLN:249250621\n"
        mock_parse.return_value = [{"name": "chr1", "length": 249250621}]

        config = {}
        result = get_chromosome_name_from_bam("test.bam", config, chromosome_number=1, reference_assembly="hg19")
        assert result == "chr1"

    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    @patch("vntyper.scripts.fastq_bam_processing.parse_contigs_from_header")
    def test_get_ncbi_chr1_grch37(self, mock_parse, mock_extract):
        """Test getting NCBI chr1 from GRCh37 BAM."""
        mock_extract.return_value = "@SQ\tSN:NC_000001.10\tLN:249250621\n"
        mock_parse.return_value = [{"name": "NC_000001.10", "length": 249250621}]

        config = {"bam_processing": {"known_chromosome_naming": {"GRCh37": {"ncbi": "NC_000001.10"}}}}
        result = get_chromosome_name_from_bam("test.bam", config, chromosome_number=1, reference_assembly="GRCh37")
        assert result == "NC_000001.10"

    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    @patch("vntyper.scripts.fastq_bam_processing.parse_contigs_from_header")
    def test_chromosome_not_found(self, mock_parse, mock_extract):
        """Test error when chromosome not found in BAM."""
        mock_extract.return_value = "@SQ\tSN:chr2\tLN:243199373\n"
        mock_parse.return_value = [{"name": "chr2", "length": 243199373}]

        config = {}
        with pytest.raises(ValueError, match="not found in BAM"):
            get_chromosome_name_from_bam("test.bam", config, chromosome_number=1, reference_assembly="hg19")

    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    @patch("vntyper.scripts.fastq_bam_processing.parse_contigs_from_header")
    def test_uppercase_contig_resolves_to_the_headers_own_spelling(self, mock_parse, mock_extract):
        """The region string is handed to samtools, which matches contig names exactly.

        A header spelling `CHR1` must come back as `CHR1`. Returning the
        canonical `chr1` builds a region samtools cannot find, and an empty
        region yields no reads and a confident negative.
        """
        mock_extract.return_value = "@SQ\tSN:CHR1\tLN:249250621\n"
        mock_parse.return_value = [{"name": "CHR1", "length": 249250621}]

        result = get_chromosome_name_from_bam("test.bam", {}, chromosome_number=1, reference_assembly="hg19")
        assert result == "CHR1"

    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    @patch("vntyper.scripts.fastq_bam_processing.parse_contigs_from_header")
    def test_exact_spelling_is_preferred_over_a_case_insensitive_match(self, mock_parse, mock_extract):
        """With both spellings present the exact one wins, whatever the list order."""
        mock_extract.return_value = "@SQ\tSN:CHR1\tLN:249250621\n@SQ\tSN:chr1\tLN:249250621\n"
        mock_parse.return_value = [{"name": "CHR1", "length": 249250621}, {"name": "chr1", "length": 249250621}]

        result = get_chromosome_name_from_bam("test.bam", {}, chromosome_number=1, reference_assembly="hg19")
        assert result == "chr1"

    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    def test_bam_read_error(self, mock_extract):
        """Test error handling when BAM cannot be read."""
        mock_extract.side_effect = Exception("File not found")

        config = {}
        with pytest.raises(ValueError, match="Cannot read BAM header"):
            get_chromosome_name_from_bam("nonexistent.bam", config, chromosome_number=1)


class TestDetectAssemblyChr1Length:
    """Test assembly detection via chr1 length (Phase 1 - Issue #139 fix)."""

    def test_detect_hg38_from_chr1_length_ucsc(self):
        """Chr1 length 248,956,422 should detect hg38/GRCh38 (UCSC naming)."""
        from vntyper.scripts.chromosome_utils import detect_assembly_from_chr1_length

        contigs = [{"name": "chr1", "length": 248956422}]
        assembly = detect_assembly_from_chr1_length(contigs)
        assert assembly in ["hg38", "GRCh38"], f"Expected hg38/GRCh38, got {assembly}"

    def test_detect_hg19_from_chr1_length_ucsc(self):
        """Chr1 length 249,250,621 should detect hg19/GRCh37 (UCSC naming)."""
        from vntyper.scripts.chromosome_utils import detect_assembly_from_chr1_length

        contigs = [{"name": "chr1", "length": 249250621}]
        assembly = detect_assembly_from_chr1_length(contigs)
        assert assembly in ["hg19", "GRCh37"], f"Expected hg19/GRCh37, got {assembly}"

    def test_detect_hg38_from_chr1_length_ensembl(self):
        """Should work with ENSEMBL '1' naming."""
        from vntyper.scripts.chromosome_utils import detect_assembly_from_chr1_length

        contigs = [{"name": "1", "length": 248956422}]
        assembly = detect_assembly_from_chr1_length(contigs)
        assert assembly in ["hg38", "GRCh38"], f"Expected hg38/GRCh38, got {assembly}"

    def test_chr1_not_found_returns_none(self):
        """Should return None if chr1 missing."""
        from vntyper.scripts.chromosome_utils import detect_assembly_from_chr1_length

        contigs = [{"name": "chr2", "length": 242193529}]
        assembly = detect_assembly_from_chr1_length(contigs)
        assert assembly is None, f"Expected None when chr1 missing, got {assembly}"

    def test_unknown_chr1_length_returns_none(self):
        """Should return None for unknown chr1 length."""
        from vntyper.scripts.chromosome_utils import detect_assembly_from_chr1_length

        contigs = [{"name": "chr1", "length": 999999999}]
        assembly = detect_assembly_from_chr1_length(contigs)
        assert assembly is None, f"Expected None for unknown length, got {assembly}"

    def test_hg38_with_195_contigs_uses_chr1_length(self):
        """
        CRITICAL: Real-world hg38 BAM with 195 contigs (Issue #139).

        BAM has 25 canonical + 170 alternates. Chr1 length check should
        reliably detect hg38 regardless of alternate contigs.
        """
        from vntyper.scripts.chromosome_utils import detect_assembly_from_chr1_length

        # Simulate hg38 structure: canonical chr1 + many alternates
        contigs = [
            {"name": "chr1", "length": 248956422},  # hg38 chr1 length
            {"name": "chr1_KI270706v1_random", "length": 175055},
            {"name": "chr1_KI270707v1_random", "length": 32032},
            # ... would be 195 total in real BAM, but chr1 is what matters
        ]
        assembly = detect_assembly_from_chr1_length(contigs)
        assert assembly in ["hg38", "GRCh38"], (
            f"Chr1 length check should detect hg38 despite alternates, got {assembly}"
        )

    def test_hg19_with_ensembl_naming_and_alternates(self):
        """
        Test assembly detection with ENSEMBL naming and alternate contigs.

        Ensures chr1 length detection works robustly with ENSEMBL '1' naming
        even when alternate contigs are present (e.g., from GENCODE references).
        Addresses Sourcery AI suggestion for comprehensive naming convention coverage.
        """
        from vntyper.scripts.chromosome_utils import detect_assembly_from_chr1_length

        # Simulate ENSEMBL structure: canonical '1' + alternates
        contigs = [
            {"name": "1", "length": 249250621},  # hg19/GRCh37 chr1 length
            {"name": "CHR_HSCHR1_1_CTG3", "length": 165050},  # ENSEMBL alternate
            {"name": "CHR_HSCHR1_2_CTG3", "length": 176043},  # ENSEMBL alternate
            # ENSEMBL can have many alternates like UCSC
        ]
        assembly = detect_assembly_from_chr1_length(contigs)
        assert assembly in ["hg19", "GRCh37"], (
            f"Chr1 length check should detect hg19/GRCh37 with ENSEMBL naming despite alternates, got {assembly}"
        )

    def test_case_insensitive_chr1_naming(self):
        """
        Test case-insensitive chr1 detection.

        Ensures detection works with uppercase/mixed case chr1 names (CHR1, Chr1, etc.)
        for robustness against non-standard naming conventions.
        Addresses Sourcery AI suggestion for defensive programming.
        """
        from vntyper.scripts.chromosome_utils import detect_assembly_from_chr1_length

        # Test uppercase CHR1
        contigs_upper = [{"name": "CHR1", "length": 248956422}]
        assembly = detect_assembly_from_chr1_length(contigs_upper)
        assert assembly in ["hg38", "GRCh38"], f"Expected hg38/GRCh38 for 'CHR1', got {assembly}"

        # Test mixed case Chr1
        contigs_mixed = [{"name": "Chr1", "length": 249250621}]
        assembly = detect_assembly_from_chr1_length(contigs_mixed)
        assert assembly in ["hg19", "GRCh37"], f"Expected hg19/GRCh37 for 'Chr1', got {assembly}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
