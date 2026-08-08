"""Cross-task contracts for the merged CRAM and single-end pipeline changes."""

from __future__ import annotations

import logging
from pathlib import Path
from unittest import mock

import pytest

from tests.support.pipeline_harness import MINIMAL_CONFIG, run_pipeline_under_harness
from vntyper.scripts import cli_handlers
from vntyper.scripts.cli_parser import build_parser

pytestmark = pytest.mark.unit


def test_cli_reference_and_single_fastq_values_survive_one_handler_contract(tmp_path: Path) -> None:
    """Task7's reference option must coexist with Task9's optional second mate."""
    parser = build_parser()
    reference = tmp_path / "reference.fa"
    args = parser.parse_args(
        [
            "pipeline",
            "--fastq1",
            "single.fastq.gz",
            "--reference-fasta",
            str(reference),
            "--output-dir",
            str(tmp_path / "out"),
        ]
    )

    with mock.patch.object(cli_handlers, "run_pipeline", autospec=True) as runner:
        cli_handlers.handle_pipeline(
            args,
            config=MINIMAL_CONFIG,
            parser=parser,
            log_level_value=logging.INFO,
            log_file_str=None,
        )

    assert runner.call_args.kwargs["fastq1"] == "single.fastq.gz"
    assert runner.call_args.kwargs["fastq2"] is None
    assert runner.call_args.kwargs["reference_fasta"] == reference


def test_cram_preflight_transport_and_four_fastq_routing_coexist_in_one_run(tmp_path: Path) -> None:
    """Task7 preflight metadata and Task9 routing must reach their consumers together."""
    cram = tmp_path / "patient.cram"
    cram.touch()
    reference = tmp_path / "reference.fa"
    produced = tuple(f"/reads/output_{name}.fastq.gz" for name in ("R1", "R2", "other", "single"))
    routed_single = produced[2]

    harness = run_pipeline_under_harness(
        tmp_path / "out",
        bam=None,
        cram=str(cram),
        reference_fasta=reference,
        stage_side_effects={
            "process_bam_to_fastq": lambda *args, **kwargs: produced,
            "route_converted_fastqs": lambda paths, config: (routed_single, None),
        },
    )

    preflight = harness.kwargs("run_preflight")
    assert preflight["reference_fasta"] == str(reference)
    assert preflight["error_output_dir"] == str(tmp_path / "out")
    assert harness.positional("route_converted_fastqs")[0] == produced
    kestrel = harness.kwargs("run_kestrel")
    assert kestrel["fastq_1"] == routed_single
    assert kestrel["fastq_2"] is None


def test_merged_pipeline_guards_still_stop_conflicts_and_single_end_shark_before_stages(tmp_path: Path) -> None:
    """Neither merge order may restore input precedence or a literal SHARK ``None`` mate."""
    bam = tmp_path / "patient.bam"
    bam.touch()
    conflict = run_pipeline_under_harness(
        tmp_path / "conflict",
        bam=str(bam),
        fastq1="/reads/single.fastq.gz",
        expect_failure=True,
    )
    shark = run_pipeline_under_harness(
        tmp_path / "shark",
        bam=None,
        cram=None,
        fastq1="/reads/single.fastq.gz",
        fastq2=None,
        extra_modules=["shark"],
        expect_failure=True,
    )

    assert isinstance(conflict.error, ValueError)
    assert isinstance(shark.error, ValueError)
    assert "SHARK requires paired-end FASTQ input" in str(shark.error)
    assert all(not stage.called for stage in conflict.stages.values())
    assert all(not stage.called for stage in shark.stages.values())
