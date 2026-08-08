"""The pipeline routes every conversion output before a FASTQ is consumed."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest import mock

import pytest

from tests.support.pipeline_harness import run_pipeline_under_harness

pytestmark = pytest.mark.unit


def test_alignment_conversion_routes_all_four_outputs_before_kestrel(tmp_path: Path) -> None:
    events: list[str] = []
    routed_single = "/reads/output_other.fastq.gz"

    def record(event: str) -> object:
        events.append(event)
        return mock.DEFAULT

    def route(*args: object, **kwargs: object) -> tuple[str, None]:
        events.append("route")
        return routed_single, None

    harness = run_pipeline_under_harness(
        tmp_path / "out",
        stage_side_effects={
            "process_bam_to_fastq": lambda *args, **kwargs: record("convert"),
            "route_converted_fastqs": route,
            "run_kestrel": lambda *args, **kwargs: record("kestrel"),
        },
    )

    assert events.index("convert") < events.index("route") < events.index("kestrel")
    assert harness.positional("route_converted_fastqs")[0] == harness.stages["process_bam_to_fastq"].return_value
    kestrel_kwargs = harness.kwargs("run_kestrel")
    assert kestrel_kwargs["fastq_1"] == routed_single
    assert kestrel_kwargs["fastq_2"] is None


def test_single_fastq_input_stays_single_through_fastp_bwa_conversion_and_kestrel(tmp_path: Path) -> None:
    routed_single = "/reads/post_alignment_other.fastq.gz"
    harness = run_pipeline_under_harness(
        tmp_path / "out",
        bam=None,
        fastq1="/input/single.fastq.gz",
        fastq2=None,
        stage_side_effects={
            "route_converted_fastqs": lambda *args, **kwargs: (routed_single, None),
        },
    )

    assert harness.positional("validate_fastq_file") == ("/input/single.fastq.gz",)
    assert harness.positional("process_fastq")[1] is None
    assert harness.positional("align_and_sort_fastq")[1] is None
    assert harness.stages["route_converted_fastqs"].call_count == 1
    kestrel_kwargs = harness.kwargs("run_kestrel")
    assert kestrel_kwargs["fastq_1"] == routed_single
    assert kestrel_kwargs["fastq_2"] is None


@pytest.mark.parametrize("alignment_kind", ["bam", "cram"])
def test_fastq_and_alignment_groups_are_rejected_before_validation(
    alignment_kind: str,
    tmp_path: Path,
) -> None:
    alignment = tmp_path / f"input.{alignment_kind}"
    alignment.touch()
    inputs: dict[str, Any] = {
        "bam": None,
        "cram": None,
        "fastq1": "/input/reads.fastq.gz",
        "fastq2": None,
        alignment_kind: str(alignment),
    }

    harness = run_pipeline_under_harness(
        tmp_path / "out",
        expect_failure=True,
        **inputs,
    )

    assert isinstance(harness.error, ValueError)
    assert "not multiples" in str(harness.error)
    assert all(not stage.called for stage in harness.stages.values())


def test_fastq2_without_fastq1_is_rejected_before_validation(tmp_path: Path) -> None:
    harness = run_pipeline_under_harness(
        tmp_path / "out",
        bam=None,
        cram=None,
        fastq1=None,
        fastq2="/input/mate2.fastq.gz",
        expect_failure=True,
    )

    assert isinstance(harness.error, ValueError)
    assert "--fastq1 must be specified" in str(harness.error)
    assert all(not stage.called for stage in harness.stages.values())


def test_single_fastq_with_shark_fails_before_a_none_operand_can_be_emitted(tmp_path: Path) -> None:
    from vntyper.modules.shark import shark_filtering

    commands: list[str] = []

    def record_command(command: str, *args: object, **kwargs: object) -> bool:
        commands.append(command)
        return True

    with mock.patch.object(shark_filtering, "run_command", side_effect=record_command):
        harness = run_pipeline_under_harness(
            tmp_path / "out",
            bam=None,
            cram=None,
            fastq1="/input/single.fastq.gz",
            fastq2=None,
            extra_modules=["shark"],
            expect_failure=True,
        )

    assert commands == [], f"SHARK command escaped the compatibility guard: {commands}"
    assert isinstance(harness.error, ValueError)
    assert "SHARK requires paired-end FASTQ input" in str(harness.error)
    assert all(not stage.called for stage in harness.stages.values())
