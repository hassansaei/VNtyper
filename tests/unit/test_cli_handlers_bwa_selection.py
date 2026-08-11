"""BWA reference selection, and the safety guard that shares it.

Two readers used to collapse every assembly onto bwa_reference_hg19/hg38 independently:
cli_handlers for the run, and cli_logging_safety._selected_bwa_reference for the check
that stops --log-file naming an operator input. Because the guard runs before
setup_logging opens the log in append mode, a guard looking at the wrong file means
logging can append into a reference FASTA.
"""

import pytest

from vntyper.scripts.cli_handlers import select_bwa_reference

pytestmark = pytest.mark.unit


def config_with(**reference_data):
    return {"reference_data": dict(reference_data)}


class TestSelection:
    def test_the_physical_key_is_preferred(self):
        cfg = config_with(bwa_reference_hg38_ensembl="/refs/ens.fa", bwa_reference_hg38="/refs/ucsc.fa")
        assert select_bwa_reference(cfg, "hg38_ensembl") == "/refs/ens.fa"

    @pytest.mark.parametrize("label", ["GRCh38", "hg38_ncbi"])
    def test_both_ncbi_labels_select_the_ncbi_reference(self, label):
        cfg = config_with(bwa_reference_GRCh38="/refs/ncbi.fna", bwa_reference_hg38="/refs/ucsc.fa")
        assert select_bwa_reference(cfg, label) == "/refs/ncbi.fna"

    def test_the_ucsc_family_fallback_warns_and_names_both_keys(self, caplog):
        cfg = config_with(bwa_reference_hg38="/refs/ucsc.fa")
        with caplog.at_level("WARNING"):
            assert select_bwa_reference(cfg, "hg38_ensembl") == "/refs/ucsc.fa"
        message = " ".join(record.getMessage() for record in caplog.records)
        assert "hg38_ensembl" in message and "bwa_reference_hg38" in message and "ucsc" in message.lower()

    def test_an_exact_null_is_authoritative_and_fails_closed(self):
        cfg = config_with(bwa_reference_hg38_ensembl=None, bwa_reference_hg38="/refs/ucsc.fa")
        with pytest.raises(ValueError, match="hg38_ensembl"):
            select_bwa_reference(cfg, "hg38_ensembl")

    def test_nothing_configured_fails_closed_naming_every_key_tried(self):
        with pytest.raises(ValueError) as excinfo:
            select_bwa_reference(config_with(), "hg38_ensembl")
        assert "bwa_reference_hg38_ensembl" in str(excinfo.value)
        assert "bwa_reference_hg38" in str(excinfo.value)
