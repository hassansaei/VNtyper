"""
Unit tests for region_utils.py

Tests region string construction and chromosome-name freshness.
"""

import logging
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts.region_utils import (
    build_region_string,
    get_region_string,
    get_region_string_with_fallback,
    resolve_assembly_alias,
)

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit


class TestBuildRegionString:
    """Test basic region string construction."""

    def test_build_ucsc_region(self):
        """Test building region string with UCSC chromosome name."""
        result = build_region_string("chr1", "155158000-155163000")
        assert result == "chr1:155158000-155163000"

    def test_build_ncbi_region(self):
        """Test building region string with NCBI accession."""
        result = build_region_string("NC_000001.10", "155158000-155163000")
        assert result == "NC_000001.10:155158000-155163000"

    def test_build_simple_region(self):
        """Test building region string with simple numeric chromosome."""
        result = build_region_string("1", "155158000-155163000")
        assert result == "1:155158000-155163000"

    def test_empty_chromosome_name(self):
        """Test error handling for empty chromosome name."""
        with pytest.raises(ValueError, match="Invalid inputs"):
            build_region_string("", "155158000-155163000")

    def test_empty_coordinates(self):
        """Test error handling for empty coordinates."""
        with pytest.raises(ValueError, match="Invalid inputs"):
            build_region_string("chr1", "")

    def test_invalid_coordinate_format(self):
        """Test error handling for invalid coordinate format."""
        with pytest.raises(ValueError, match="Invalid coordinate format"):
            build_region_string("chr1", "155158000")  # Missing hyphen

    def test_invalid_coordinate_values(self):
        """Test error handling for non-numeric coordinates."""
        with pytest.raises(ValueError, match="Invalid coordinate values"):
            build_region_string("chr1", "abc-def")


class TestResolveAssemblyAlias:
    """Test assembly alias resolution to coordinate systems."""

    def test_resolve_hg19_to_grch37(self):
        """Test resolving hg19 to GRCh37 coordinate system."""
        assert resolve_assembly_alias("hg19") == "GRCh37"

    def test_resolve_grch37(self):
        """Test resolving GRCh37 (no change)."""
        assert resolve_assembly_alias("GRCh37") == "GRCh37"

    def test_resolve_hg38_to_grch38(self):
        """Test resolving hg38 to GRCh38 coordinate system."""
        assert resolve_assembly_alias("hg38") == "GRCh38"

    def test_resolve_grch38(self):
        """Test resolving GRCh38 (no change)."""
        assert resolve_assembly_alias("GRCh38") == "GRCh38"

    def test_unknown_assembly(self):
        """Test handling of unknown assembly (defaults to GRCh37)."""
        assert resolve_assembly_alias("unknown") == "GRCh37"


class TestGetRegionString:
    """Test dynamic region string generation."""

    @patch("vntyper.scripts.chromosome_utils.get_chromosome_name_from_bam")
    def test_get_region_hg19_ucsc(self, mock_get_chr):
        """Test getting region string for hg19 with UCSC naming."""
        mock_get_chr.return_value = "chr1"

        config = {
            "bam_processing": {
                "assemblies": {
                    "GRCh37": {
                        "bam_region_coords": "155158000-155163000",
                        "vntr_region_coords": "155160500-155162000",
                        "chromosome": 1,
                    }
                }
            }
        }

        result = get_region_string("test.bam", "hg19", "bam_region_coords", config)
        assert result == "chr1:155158000-155163000"

    @patch("vntyper.scripts.chromosome_utils.get_chromosome_name_from_bam")
    def test_get_region_grch37_ncbi(self, mock_get_chr):
        """Test getting region string for GRCh37 with NCBI naming."""
        mock_get_chr.return_value = "NC_000001.10"

        config = {
            "bam_processing": {
                "assemblies": {"GRCh37": {"bam_region_coords": "155158000-155163000", "chromosome": 1}},
                "known_chromosome_naming": {"GRCh37": {"ncbi": "NC_000001.10"}},
            }
        }

        result = get_region_string("test.bam", "GRCh37", "bam_region_coords", config)
        assert result == "NC_000001.10:155158000-155163000"

    @patch("vntyper.scripts.chromosome_utils.get_chromosome_name_from_bam")
    def test_get_region_hg38(self, mock_get_chr):
        """Test getting region string for hg38."""
        mock_get_chr.return_value = "chr1"

        config = {
            "bam_processing": {"assemblies": {"GRCh38": {"bam_region_coords": "155184000-155194000", "chromosome": 1}}}
        }

        result = get_region_string("test.bam", "hg38", "bam_region_coords", config)
        assert result == "chr1:155184000-155194000"

    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    def test_mutating_naming_policy_for_the_same_path_is_observed(self, mock_extract):
        """Config-driven naming cannot reuse a decision made under old policy."""
        mock_extract.return_value = "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:1\tLN:100\n@SQ\tSN:2\tLN:100\n"
        config = {
            "bam_processing": {
                "assemblies": {
                    "GRCh37": {
                        "bam_region_coords": "155158000-155163000",
                        "vntr_region_coords": "155160500-155162000",
                        "chromosome": 1,
                    }
                }
            }
        }

        assert get_region_string("same.bam", "hg19", "bam_region_coords", config) == "1:155158000-155163000"

        config["assembly_detection"] = {
            "primary_contig_patterns": [r"^chr[0-9XYM]+$", r"^NC_\d{6}\.\d+$", r"^never$"],
        }
        assert get_region_string("same.bam", "hg19", "bam_region_coords", config) == "chr1:155158000-155163000"
        assert mock_extract.call_count == 2

    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    def test_replacing_a_file_at_the_same_path_is_observed(self, mock_extract, tmp_path: Path):
        """A same-path replacement cannot inherit the previous file's naming."""
        alignment = tmp_path / "sample.bam"
        alignment.write_text("@SQ\tSN:chr1\tLN:100\n")
        mock_extract.side_effect = lambda bam_file, _config: Path(bam_file).read_text()
        config = {
            "bam_processing": {"assemblies": {"GRCh37": {"bam_region_coords": "155158000-155163000", "chromosome": 1}}}
        }

        assert get_region_string(str(alignment), "hg19", "bam_region_coords", config) == ("chr1:155158000-155163000")

        alignment.write_text("@SQ\tSN:1\tLN:100\n")
        assert get_region_string(str(alignment), "hg19", "bam_region_coords", config) == "1:155158000-155163000"
        assert mock_extract.call_count == 2

    def test_missing_assembly_config(self):
        """Test error when assembly not found in config."""
        config = {"bam_processing": {}}

        with pytest.raises(KeyError, match="Configuration missing"):
            get_region_string("test.bam", "hg19", "bam_region_coords", config)

    def test_missing_region_type(self):
        """Test error when region type not found."""
        config = {"bam_processing": {"assemblies": {"GRCh37": {"chromosome": 1}}}}

        with pytest.raises(KeyError, match="Region type"):
            get_region_string("test.bam", "hg19", "nonexistent_region", config)


class TestGetRegionStringWithFallback:
    """Test region string resolution with fallback to legacy format."""

    @patch("vntyper.scripts.region_utils.get_region_string")
    def test_new_format_success(self, mock_get_region):
        """Test successful resolution with new format."""
        mock_get_region.return_value = "chr1:155158000-155163000"

        config = {
            "bam_processing": {"assemblies": {"GRCh37": {"bam_region_coords": "155158000-155163000", "chromosome": 1}}}
        }

        result = get_region_string_with_fallback("test.bam", "hg19", "bam_region", config)
        assert result == "chr1:155158000-155163000"

    @patch("vntyper.scripts.region_utils.get_region_string")
    def test_fallback_to_legacy(self, mock_get_region):
        """Test fallback to legacy config format."""
        mock_get_region.side_effect = KeyError("Missing config")

        config = {"bam_processing": {"bam_region_hg19": "chr1:155158000-155163000"}}

        result = get_region_string_with_fallback("test.bam", "hg19", "bam_region", config)
        assert result == "chr1:155158000-155163000"

    @patch("vntyper.scripts.region_utils.get_region_string")
    def test_both_formats_fail(self, mock_get_region):
        """Test error when both new and legacy formats fail."""
        mock_get_region.side_effect = KeyError("Missing config")

        config = {"bam_processing": {}}

        with pytest.raises(ValueError, match="Neither new nor legacy format"):
            get_region_string_with_fallback("test.bam", "hg19", "bam_region", config)

    @pytest.mark.parametrize(
        "header",
        [
            "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:1\tLN:100\n",
            "@SQ\tSN:scaffold_a\tLN:100\n@SQ\tSN:scaffold_b\tLN:100\n",
        ],
        ids=["fifty-fifty", "zero-classifiable"],
    )
    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    def test_naming_ambiguity_is_not_swallowed_by_the_legacy_fallback(self, mock_extract, header):
        """Target resolution must fail rather than select any configured fallback."""
        mock_extract.return_value = header
        config = {
            "bam_processing": {
                "assemblies": {"GRCh37": {"bam_region_coords": "155158000-155163000", "chromosome": 1}},
                "bam_region_hg19": "chr1:155158000-155163000",
            }
        }

        with pytest.raises(ValueError, match="Ambiguous or unclassifiable chromosome naming convention"):
            get_region_string_with_fallback("ambiguous.bam", "hg19", "bam_region", config)


class TestTheLegacyFallbackCannotReturnARegionTheBamDoesNotContain:
    """
    The legacy fallback used to return a region string matching nothing at all.

    ``get_region_string`` raises ``ValueError`` whenever it cannot find the
    chromosome in the BAM header, and the fallback swallowed that and returned
    ``config["bam_processing"]["bam_region_<assembly>"]`` instead. Those legacy
    keys are UCSC-named for the ``hg19``/``hg38`` spellings, so on an NCBI-named
    BAM the pipeline sliced ``chr1:155158000-155163000`` from a file whose only
    chromosome 1 is ``NC_000001.11``. samtools returns zero records for a contig
    it does not have and exits 0. The slice is empty, Kestrel sees no reads, and
    the report states a confident negative.

    Whether it happens turns on the *spelling* of ``--reference-assembly``, which
    is why it went unnoticed: ``GRCh38`` picks the key holding an NCBI accession
    and works, while ``hg38`` picks the UCSC one and silently does not.

    The fallback itself is load-bearing and stays: on an hg38 BAM carrying ~170
    alternate contigs, fewer than half the names match the UCSC pattern, so
    ``detect_naming_convention`` returns ``"unknown"``, the dynamic path raises,
    and the legacy ``chr1:...`` string is correct. The fix is only to refuse a
    fallback region whose contig is demonstrably absent.
    """

    NCBI_HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:NC_000001.11\tLN:248956422\n@SQ\tSN:NC_000002.12\tLN:242193529\n"
    UCSC_HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n@SQ\tSN:chr1_KI270706v1_random\tLN:175055\n"

    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    @patch("vntyper.scripts.region_utils.get_region_string")
    def test_a_ucsc_fallback_region_on_an_ncbi_bam_is_refused(self, mock_get_region, mock_header):
        """
        The silent false negative, made loud.

        Before the guard this returned ``"chr1:155184000-155194000"`` and the
        pipeline went on to slice nothing out of an NCBI-named BAM.
        """
        mock_get_region.side_effect = ValueError("Chromosome 1 not found in BAM file")
        mock_header.return_value = self.NCBI_HEADER

        config = {"bam_processing": {"bam_region_hg38": "chr1:155184000-155194000"}}

        with pytest.raises(ValueError, match="chr1"):
            get_region_string_with_fallback("sample.bam", "hg38", "bam_region", config)

    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    @patch("vntyper.scripts.region_utils.get_region_string")
    def test_a_fallback_region_the_bam_does_contain_is_returned_unchanged(self, mock_get_region, mock_header):
        """
        The hg38-with-alternate-contigs case, which the fallback exists to rescue.

        ``chr1`` is genuinely in this header; the dynamic path failed only because
        the alternate contigs outvoted the primary ones in convention detection.
        Behaviour here must be exactly what it was.
        """
        mock_get_region.side_effect = ValueError("Chromosome 1 not found in BAM file")
        mock_header.return_value = self.UCSC_HEADER

        config = {"bam_processing": {"bam_region_hg38": "chr1:155184000-155194000"}}

        result = get_region_string_with_fallback("sample.bam", "hg38", "bam_region", config)

        assert result == "chr1:155184000-155194000"

    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    @patch("vntyper.scripts.region_utils.get_region_string")
    def test_the_match_is_case_insensitive(self, mock_get_region, mock_header):
        """``get_chromosome_name_from_bam`` matches case-insensitively; so must this."""
        mock_get_region.side_effect = ValueError("Chromosome 1 not found in BAM file")
        mock_header.return_value = "@SQ\tSN:CHR1\tLN:248956422\n"

        config = {"bam_processing": {"bam_region_hg38": "chr1:155184000-155194000"}}

        assert get_region_string_with_fallback("sample.bam", "hg38", "bam_region", config) == (
            "chr1:155184000-155194000"
        )

    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    @patch("vntyper.scripts.region_utils.get_region_string")
    def test_an_unreadable_header_keeps_the_previous_behaviour(self, mock_get_region, mock_header):
        """
        When the header cannot be read there is nothing to check against.

        Refusing here would turn a working run into a failure on the strength of
        no evidence at all, so the legacy region is returned exactly as before and
        the inability to verify is logged.
        """
        mock_get_region.side_effect = ValueError("Chromosome 1 not found in BAM file")
        mock_header.side_effect = OSError("samtools not found")

        config = {"bam_processing": {"bam_region_hg38": "chr1:155184000-155194000"}}

        result = get_region_string_with_fallback("sample.bam", "hg38", "bam_region", config)

        assert result == "chr1:155184000-155194000"

    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    @patch("vntyper.scripts.region_utils.get_region_string")
    def test_a_header_with_no_contigs_keeps_the_previous_behaviour(self, mock_get_region, mock_header):
        """An empty contig list is absence of evidence, not evidence of absence."""
        mock_get_region.side_effect = ValueError("Chromosome 1 not found in BAM file")
        mock_header.return_value = "@HD\tVN:1.6\n"

        config = {"bam_processing": {"bam_region_hg38": "chr1:155184000-155194000"}}

        assert get_region_string_with_fallback("sample.bam", "hg38", "bam_region", config) == (
            "chr1:155184000-155194000"
        )

    @patch("vntyper.scripts.fastq_bam_processing.extract_bam_header")
    @patch("vntyper.scripts.region_utils.get_region_string")
    def test_the_refusal_names_the_contig_and_what_the_bam_actually_has(self, mock_get_region, mock_header, caplog):
        """
        The error must be actionable, because the fix is a CLI flag away.

        Re-running with ``--reference-assembly GRCh38`` selects the legacy key
        holding the NCBI accession and the same BAM works.
        """
        mock_get_region.side_effect = ValueError("Chromosome 1 not found in BAM file")
        mock_header.return_value = self.NCBI_HEADER

        config = {"bam_processing": {"bam_region_hg38": "chr1:155184000-155194000"}}

        with (
            caplog.at_level(logging.INFO, logger="vntyper.scripts.region_utils"),
            pytest.raises(ValueError) as excinfo,
        ):
            get_region_string_with_fallback("sample.bam", "hg38", "bam_region", config)

        message = str(excinfo.value)
        assert "chr1" in message
        assert "NC_000001.11" in message, "the message must show what the BAM does contain"
        assert [r for r in caplog.records if r.levelno == logging.ERROR], (
            "refusing to slice must be logged at ERROR, not swallowed"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
