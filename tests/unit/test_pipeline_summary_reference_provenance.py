"""``pipeline_summary.json`` must record the reference a run actually used.

MAJOR 5 (milestone-5 PR-2 review): `handle_pipeline` resolves the BWA reference and
`run_pipeline` used to record it into `pipeline_summary.json` (and, from there, into the
HTML report) for every input type - including a BAM run, which never reads a reference
file at all, and a CRAM run, which decodes against whatever `alignment_preflight`
actually resolved and can differ entirely from the configured BWA path. That is a false
provenance claim in a machine-readable artefact, in a milestone whose whole thesis is
that the recorded reference is the one actually used.

These drive the real ``run_pipeline`` through :func:`tests.support.pipeline_harness.
run_pipeline_under_harness` - every stage stubbed except ``summary.record_step`` /
``write_summary``, which are the pipeline's own record of what it did - and assert on
the ``pipeline_summary.json`` that lands on disk, not only on an in-memory dict.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from tests.support.pipeline_harness import run_pipeline_under_harness
from vntyper.scripts.alignment_contract import AlignmentPlan

pytestmark = pytest.mark.unit


def _plan_with(*, reference_path: str | None, reference_source: str):
    """Return a ``run_preflight`` side effect that pins the plan's reference fields.

    Args:
        reference_path: The ``AlignmentPlan.reference_path`` the stubbed preflight call
            should report.
        reference_source: The ``AlignmentPlan.reference_source`` it should report.

    Returns:
        A callable usable as ``stage_side_effects["run_preflight"]``.
    """

    def _run_preflight(**kwargs: Any) -> AlignmentPlan:
        file_format = kwargs["file_format"]
        output_dir = Path(kwargs["output_dir"])
        output_name = kwargs["output_name"]
        view_path = output_dir / f"{output_name}.{file_format}"
        index_suffix = "bai" if file_format == "bam" else "crai"
        return AlignmentPlan(
            input_path=str(kwargs["in_path"]),
            view_path=str(view_path),
            file_format=file_format,
            index_path=f"{view_path}.{index_suffix}",
            reference_path=reference_path,
            reference_source=reference_source,
            uncovered_contigs=(),
            unmapped_scan="indexed",
            binding=kwargs.get("binding"),
        )

    return _run_preflight


def _read_summary(output_dir: Path) -> dict[str, Any]:
    return json.loads((output_dir / "pipeline_summary.json").read_text(encoding="utf-8"))


def test_a_fastq_run_still_records_the_bwa_reference_it_aligned_with(tmp_path: Path) -> None:
    """The one case where the BWA resolution *is* correct, and the reason the fields
    exist at all - a FASTQ run must not lose this when BAM/CRAM stop claiming it."""
    harness = run_pipeline_under_harness(
        tmp_path,
        fastq1="/reads/single.fastq.gz",
        fastq2=None,
        bwa_reference="/refs/hg19.fa",
        reference_key_used="bwa_reference_hg19",
        reference_source_effective="ucsc",
    )

    assert harness.error is None
    summary = _read_summary(harness.output_dir)
    assert summary["reference_assembly_requested"] == "hg19"
    assert summary["reference_key_used"] == "bwa_reference_hg19"
    assert summary["reference_path"] == "/refs/hg19.fa"
    assert summary["reference_source_effective"] == "ucsc"


def test_a_bam_run_does_not_claim_a_reference_it_never_read(tmp_path: Path) -> None:
    """A BAM run's alignment plan already reports `reference_path=None`,
    `reference_source="not-required"` in production - the summary must read that
    through, not fall back to naming the configured (unopened) BWA path."""
    bam = tmp_path.parent / f"{tmp_path.name}_input" / "sample.bam"
    bam.parent.mkdir(parents=True, exist_ok=True)
    bam.touch()

    harness = run_pipeline_under_harness(
        tmp_path,
        bam=str(bam),
        bwa_reference="/refs/hg19.fa",
        reference_key_used="bwa_reference_hg19",
        reference_source_effective="ucsc",
        stage_side_effects={"run_preflight": _plan_with(reference_path=None, reference_source="not-required")},
    )

    assert harness.error is None
    summary = _read_summary(harness.output_dir)
    # The requested assembly is still legitimately recorded - it drives region
    # selection even for a BAM run - but nothing claims a reference file was opened.
    assert summary["reference_assembly_requested"] == "hg19"
    assert summary["reference_key_used"] is None
    assert summary["reference_path"] is None
    assert summary["reference_source_effective"] == "not-required"


def test_a_cram_run_records_the_reference_it_actually_decoded_with(tmp_path: Path) -> None:
    """The literal MAJOR 5 scenario: a CRAM decoded against a candidate distinct from
    the configured BWA path must record *that* reference, not the BWA selection."""
    cram = tmp_path.parent / f"{tmp_path.name}_input" / "sample.cram"
    cram.parent.mkdir(parents=True, exist_ok=True)
    cram.touch()

    harness = run_pipeline_under_harness(
        tmp_path,
        bam=None,
        cram=str(cram),
        bwa_reference="/refs/bwa/hg19.fa",
        reference_key_used="bwa_reference_hg19",
        reference_source_effective="ucsc",
        stage_side_effects={
            "run_preflight": _plan_with(reference_path="/refs/cram/explicit-reference-fasta.fa", reference_source="cli")
        },
    )

    assert harness.error is None
    summary = _read_summary(harness.output_dir)
    assert summary["reference_path"] == "/refs/cram/explicit-reference-fasta.fa"
    assert summary["reference_path"] != "/refs/bwa/hg19.fa", (
        "the CRAM run decoded against a reference distinct from the configured BWA path"
    )
    assert summary["reference_source_effective"] == "cli"
    assert summary["reference_key_used"] is None
