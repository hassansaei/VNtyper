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
        with open(config_path, "r") as f:
            self.config = json.load(f)

    def test_detect_grch37_from_contigs(self):
        """GRCh37 contigs (no chr prefix) should be detected."""
        header = (
            "@SQ\tSN:1\tLN:249250621\n"
            "@SQ\tSN:2\tLN:243199373\n"
            "@SQ\tSN:3\tLN:198022430\n"
        )
        result = detect_assembly_from_contigs(header, self.config, threshold=0.1)
        assert result == "GRCh37"

    def test_detect_grch38_from_contigs(self):
        """GRCh38 contigs (no chr prefix) should be detected."""
        header = (
            "@SQ\tSN:1\tLN:248956422\n"
            "@SQ\tSN:2\tLN:242193529\n"
            "@SQ\tSN:3\tLN:198295559\n"
        )
        result = detect_assembly_from_contigs(header, self.config, threshold=0.1)
        assert result == "GRCh38"

    def test_detect_hg19_still_works(self):
        """Ensure hg19 detection unchanged."""
        header = (
            "@SQ\tSN:chr1\tLN:249250621\n"
            "@SQ\tSN:chr2\tLN:243199373\n"
            "@SQ\tSN:chr3\tLN:198022430\n"
        )
        result = detect_assembly_from_contigs(header, self.config, threshold=0.1)
        assert result == "hg19"

    def test_detect_hg38_still_works(self):
        """Ensure hg38 detection unchanged."""
        header = (
            "@SQ\tSN:chr1\tLN:248956422\n"
            "@SQ\tSN:chr2\tLN:242193529\n"
            "@SQ\tSN:chr3\tLN:198295559\n"
        )
        result = detect_assembly_from_contigs(header, self.config, threshold=0.1)
        assert result == "hg38"
