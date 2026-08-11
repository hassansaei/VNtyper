"""BWA reference selection, and the safety guard that shares it.

Two readers used to collapse every assembly onto bwa_reference_hg19/hg38 independently:
cli_handlers for the run, and cli_logging_safety._selected_bwa_reference for the check
that stops --log-file naming an operator input. Because the guard runs before
setup_logging opens the log in append mode, a guard looking at the wrong file means
logging can append into a reference FASTA.
"""

import logging
from unittest import mock

import pytest

from vntyper.scripts import cli_handlers
from vntyper.scripts.artifact_names import PIPELINE_BASENAME
from vntyper.scripts.cli_handlers import select_bwa_reference
from vntyper.scripts.cli_parser import build_parser

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

    def test_a_malformed_present_value_fails_closed_with_value_error_not_type_error(self):
        """A present-but-non-string config value (an int here) must not reach `Path()`
        further down the call chain and crash with a bare `TypeError`; the required run
        path must fail closed the same way a missing key does."""
        cfg = config_with(bwa_reference_hg19=7)
        with pytest.raises(ValueError):
            select_bwa_reference(cfg, "hg19")


NO_BWA_CONFIG: dict = {
    "default_values": {
        "output_dir": "out",
        "threads": 4,
        "reference_assembly": "hg19",
        "output_name": PIPELINE_BASENAME,
        "archive_format": "zip",
    },
    "reference_data": {},
}


def _run_pipeline_handler(argv: list[str]) -> mock.MagicMock:
    """Parse ``argv`` and run ``handle_pipeline`` against :data:`NO_BWA_CONFIG`.

    Args:
        argv: Arguments after ``vntyper``, starting with ``pipeline``.

    Returns:
        unittest.mock.MagicMock: The ``run_pipeline`` stub.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    with mock.patch.object(cli_handlers, "run_pipeline", autospec=True) as stub:
        cli_handlers.handle_pipeline(
            args,
            config=NO_BWA_CONFIG,
            parser=parser,
            log_level_value=logging.INFO,
            log_file_str=None,
        )
    return stub


class TestTheRunPathOnlyRequiresBwaWhenItWillAlign:
    """`pipeline_inputs.py` raises for a missing BWA reference only when
    ``input_type == "FASTQ"``; the BAM and CRAM branches never read it
    (`if input_type == "FASTQ" and not bwa_reference: raise ValueError(...)`,
    `pipeline_inputs.py:153`). Resolving with `required=True` unconditionally would make
    `handle_pipeline` stricter than the consumer that actually needs the value, aborting
    a BAM/CRAM run a config with no BWA keys used to complete.
    """

    def test_a_bam_run_with_no_bwa_reference_configured_does_not_raise(self, tmp_path) -> None:
        stub = _run_pipeline_handler(["pipeline", "-o", str(tmp_path), "--bam", "in.bam"])
        assert stub.call_count == 1
        assert stub.call_args.kwargs["bwa_reference"] is None

    def test_a_cram_run_with_no_bwa_reference_configured_does_not_raise(self, tmp_path) -> None:
        stub = _run_pipeline_handler(["pipeline", "-o", str(tmp_path), "--cram", "in.cram"])
        assert stub.call_count == 1
        assert stub.call_args.kwargs["bwa_reference"] is None

    def test_a_fastq_run_with_no_bwa_reference_configured_still_fails_closed(self, tmp_path) -> None:
        with pytest.raises(ValueError, match="bwa_reference_hg19"):
            _run_pipeline_handler(
                ["pipeline", "-o", str(tmp_path), "--fastq1", "r1.fq", "--fastq2", "r2.fq"],
            )
