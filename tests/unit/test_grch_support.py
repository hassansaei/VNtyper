"""Unit tests for GRCh37/GRCh38 support."""

import json
from pathlib import Path

import pytest

from vntyper.scripts.fastq_bam_processing import detect_assembly_from_contigs

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit


class TestGRChSupport:
    """Test NCBI reference assembly support."""

    def setup_method(self):
        """Load config for each test."""
        config_path = Path("vntyper/config.json")
        with open(config_path) as f:
            self.config = json.load(f)

    def test_detect_grch37_from_contigs(self):
        """GRCh37 contigs (no chr prefix) should be detected."""
        header = "@SQ\tSN:1\tLN:249250621\n@SQ\tSN:2\tLN:243199373\n@SQ\tSN:3\tLN:198022430\n"
        result = detect_assembly_from_contigs(header, self.config, threshold=0.1)
        assert result == "GRCh37"

    def test_detect_grch38_from_contigs(self):
        """GRCh38 contigs (no chr prefix) should be detected."""
        header = "@SQ\tSN:1\tLN:248956422\n@SQ\tSN:2\tLN:242193529\n@SQ\tSN:3\tLN:198295559\n"
        result = detect_assembly_from_contigs(header, self.config, threshold=0.1)
        assert result == "GRCh38"

    def test_detect_hg19_still_works(self):
        """Ensure hg19 detection unchanged."""
        header = "@SQ\tSN:chr1\tLN:249250621\n@SQ\tSN:chr2\tLN:243199373\n@SQ\tSN:chr3\tLN:198022430\n"
        result = detect_assembly_from_contigs(header, self.config, threshold=0.1)
        assert result == "hg19"

    def test_detect_hg38_still_works(self):
        """Ensure hg38 detection unchanged."""
        header = "@SQ\tSN:chr1\tLN:248956422\n@SQ\tSN:chr2\tLN:242193529\n@SQ\tSN:chr3\tLN:198295559\n"
        result = detect_assembly_from_contigs(header, self.config, threshold=0.1)
        assert result == "hg38"

    def test_detect_subset_bam_chr1_only(self):
        """Subset BAMs with only chr1 should be detected via chr1 fallback."""
        # UCSC naming
        header_hg38 = "@SQ\tSN:chr1\tLN:248956422\n"
        assert detect_assembly_from_contigs(header_hg38, self.config) == "hg38"

        header_hg19 = "@SQ\tSN:chr1\tLN:249250621\n"
        assert detect_assembly_from_contigs(header_hg19, self.config) == "hg19"

        # NCBI/Ensembl naming
        header_grch38 = "@SQ\tSN:1\tLN:248956422\n"
        assert detect_assembly_from_contigs(header_grch38, self.config) == "GRCh38"

        header_grch37 = "@SQ\tSN:1\tLN:249250621\n"
        assert detect_assembly_from_contigs(header_grch37, self.config) == "GRCh37"

        # Unknown length
        header_unknown = "@SQ\tSN:chr1\tLN:999999999\n"
        assert detect_assembly_from_contigs(header_unknown, self.config) == "Not detected"

    def test_parse_header_pipeline_info_disambiguates_remapped_bams(self, tmp_path):
        """When a BAM header mentions hg19 in an old fastq path but is aligned to hg38 chr1,
        contig evidence should disambiguate assembly_text to hg38 instead of falsely flagging hg19."""
        from vntyper.scripts.fastq_bam_processing import parse_header_pipeline_info

        header = (
            "@HD\tVN:1.5\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:248956422\n"
            "@PG\tID:bwa\tPN:bwa\tVN:0.7.18\tCL:bwa mem reference/alignment/chr1.hg38.fa fastqs/example_hg19_R1.fastq.gz\n"
        )
        parse_header_pipeline_info(header, tmp_path, self.config)
        with open(tmp_path / "pipeline_info.json") as f:
            data = json.load(f)

        assert data["assembly_text"] == "hg38"
        assert data["assembly_contig"] == "hg38"
        assert data["alignment_pipeline"] == "BWA"

    @pytest.mark.parametrize(
        ("header", "expected_text", "expected_contig"),
        [
            (
                "@SQ\tSN:1\tLN:248956422\n@PG\tID:bwa\tCL:bwa mem GRCh38.fa\n",
                "hg38",
                "GRCh38",
            ),
            (
                "@SQ\tSN:1\tLN:249250621\n@PG\tID:bwa\tCL:bwa mem GRCh37.fa\n",
                "hg19",
                "GRCh37",
            ),
            (
                "@SQ\tSN:chr1\tLN:249250621\n@PG\tID:bwa\tCL:bwa mem hg19.fa\n",
                "hg19",
                "hg19",
            ),
            (
                "@HD\tVN:1.5\n@PG\tID:bwa\tCL:bwa mem hg38.fa\n",
                "hg38",
                "Not detected",
            ),
            (
                "@HD\tVN:1.5\n@PG\tID:bwa\tCL:bwa mem hg19.fa\n",
                "hg19",
                "Not detected",
            ),
            (
                "@HD\tVN:1.5\n@PG\tID:bwa\tCL:bwa mem hg38.fa sample_hg19.fq\n",
                "Not detected",
                "Not detected",
            ),
            (
                "@HD\tVN:1.5\n@PG\tID:bwa\tCL:bwa mem unannotated.fa sample.fq\n",
                "Not detected",
                "Not detected",
            ),
        ],
    )
    def test_parse_header_pipeline_info_branches(self, tmp_path, header, expected_text, expected_contig):
        from vntyper.scripts.fastq_bam_processing import parse_header_pipeline_info

        out_file = "info.json"
        parse_header_pipeline_info(header, tmp_path, self.config, output_name=out_file)
        with open(tmp_path / out_file) as f:
            data = json.load(f)

        assert data["assembly_text"] == expected_text
        assert data["assembly_contig"] == expected_contig

    def test_subset_bam_without_text_keeps_text_undetected(self, tmp_path):
        """Header text and contig evidence stay separate, so the report can show both."""
        from vntyper.scripts.fastq_bam_processing import parse_header_pipeline_info

        parse_header_pipeline_info("@SQ\tSN:chr1\tLN:248956422\n", tmp_path, self.config)
        data = json.loads((tmp_path / "pipeline_info.json").read_text())
        assert data["assembly_text"] == "Not detected"
        assert data["assembly_contig"] == "hg38"
