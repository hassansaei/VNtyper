"""Tests for pipeline resume orchestration and output directory warnings (#20)."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import pytest

from tests.support.pipeline_harness import run_pipeline_under_harness
from vntyper.scripts import summary_steps
from vntyper.scripts.run_configuration import resolve_run_configuration
from vntyper.version import __version__ as VERSION

pytestmark = pytest.mark.unit


def test_warning_emitted_when_output_dir_non_empty_without_resume(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    output_dir = tmp_path / "run_dir"
    output_dir.mkdir()
    (output_dir / "pre_existing_file.txt").touch()

    with caplog.at_level(logging.WARNING):
        run_pipeline_under_harness(output_dir=output_dir, run_pipeline_kwargs={"resume": False})

    assert any("is non-empty" in record.message for record in caplog.records)


def test_no_warning_emitted_when_output_dir_empty_without_resume(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    output_dir = tmp_path / "empty_run_dir"
    output_dir.mkdir()

    with caplog.at_level(logging.WARNING):
        run_pipeline_under_harness(output_dir=output_dir, run_pipeline_kwargs={"resume": False})

    assert not any("is non-empty" in record.message for record in caplog.records)


def test_no_warning_emitted_when_resume_is_active(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    output_dir = tmp_path / "resumed_dir"
    output_dir.mkdir()
    (output_dir / "pre_existing_file.txt").touch()

    with caplog.at_level(logging.WARNING):
        # Even if non-empty, resume expects pre-existing files and should not emit this warning
        run_pipeline_under_harness(
            output_dir=output_dir,
            run_pipeline_kwargs={"resume": True},
            expect_failure=True,  # Fails due to missing prior summary, not warning
        )

    assert not any("is non-empty; prior results may be overwritten" in record.message for record in caplog.records)


def test_resume_fatal_refusal_when_prior_summary_missing(tmp_path: Path) -> None:
    output_dir = tmp_path / "empty_dir"
    output_dir.mkdir()

    with pytest.raises(AssertionError, match="run_pipeline exited"):
        run_pipeline_under_harness(
            output_dir=output_dir,
            resume=True,
            expect_failure=False,
        )


def test_resume_fatal_refusal_on_input_files_mismatch(tmp_path: Path) -> None:
    output_dir = tmp_path / "mismatch_dir"
    output_dir.mkdir()

    prior_summary = {
        "schema_version": 3,
        "version": VERSION,
        "input_files": {"bam": "/old/input.bam"},
        "sample_name": "sample",
        "reference_key_used": None,
        "decision_profile_sha256": resolve_run_configuration().decision_profile.digest,
        "pipeline_start": "2026-09-05T09:00:00.000000",
        "steps": [],
    }
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    with pytest.raises(AssertionError, match="run_pipeline exited"):
        run_pipeline_under_harness(
            output_dir=output_dir,
            bam="/new/input.bam",
            resume=True,
            expect_failure=False,
        )


def test_resume_skips_kestrel_and_advntr_when_reusable(tmp_path: Path) -> None:
    output_dir = tmp_path / "success_resume_dir"
    output_dir.mkdir()
    input_dir = tmp_path / "input_dir"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.touch()

    kestrel_dir = output_dir / "kestrel"
    kestrel_dir.mkdir()
    kestrel_tsv = kestrel_dir / "kestrel_result.tsv"
    kestrel_tsv.write_text("kestrel tsv data", encoding="utf-8")
    (kestrel_dir / "output.vcf").write_text("vcf", encoding="utf-8")
    (kestrel_dir / "output.bam").write_text("bam", encoding="utf-8")
    (kestrel_dir / "kestrel_pre_result.tsv").write_text("pre", encoding="utf-8")

    kestrel_md5 = hashlib.md5(b"kestrel tsv data").hexdigest()

    advntr_dir = output_dir / "advntr"
    advntr_dir.mkdir()
    advntr_tsv = advntr_dir / "output_adVNTR_result.tsv"
    advntr_tsv.write_text("advntr data", encoding="utf-8")
    advntr_md5 = hashlib.md5(b"advntr data").hexdigest()

    prior_summary = {
        "schema_version": 3,
        "version": VERSION,
        "input_files": {"bam": input_bam.name},
        "sample_name": "sample",
        "reference_key_used": None,
        "decision_profile_sha256": resolve_run_configuration().decision_profile.digest,
        "pipeline_start": "2026-09-05T08:00:00.000000",
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "start": "2026-09-05T08:05:00",
                "end": "2026-09-05T08:10:00",
                "command": "run_kestrel(...)",
                "result_file": str(kestrel_tsv),
                "file_type": "tsv",
                "md5sum": kestrel_md5,
                "parsed_result": {"comments": [], "data": [{"Motif": "X"}]},
            },
            {
                "step": summary_steps.STEP_ADVNTR,
                "start": "2026-09-05T08:10:00",
                "end": "2026-09-05T08:15:00",
                "command": "run_advntr(...)",
                "result_file": str(advntr_tsv),
                "file_type": "tsv",
                "md5sum": advntr_md5,
                "parsed_result": {"comments": [], "data": []},
            },
        ],
    }
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    harness = run_pipeline_under_harness(
        output_dir=output_dir,
        extra_modules=["advntr"],
        bam=str(input_bam),
        resume=True,
    )

    # Stages should have been skipped!
    harness.stages["run_advntr"].assert_not_called()

    # The resulting summary should have carried forward the steps with reused_from
    final_summary = json.loads((output_dir / "pipeline_summary.json").read_text(encoding="utf-8"))
    step_map = {s["step"]: s for s in final_summary["steps"]}

    assert summary_steps.STEP_KESTREL in step_map
    assert step_map[summary_steps.STEP_KESTREL]["reused_from"] == "2026-09-05T08:00:00.000000"
    assert summary_steps.STEP_ADVNTR in step_map
    assert step_map[summary_steps.STEP_ADVNTR]["reused_from"] == "2026-09-05T08:00:00.000000"


def test_resume_forces_kestrel_rerun_when_corrupted(tmp_path: Path) -> None:
    output_dir = tmp_path / "rerun_dir"
    output_dir.mkdir()
    input_dir = tmp_path / "input_dir"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.touch()

    kestrel_dir = output_dir / "kestrel"
    kestrel_dir.mkdir()
    kestrel_tsv = kestrel_dir / "kestrel_result.tsv"
    kestrel_tsv.write_text("corrupted content", encoding="utf-8")
    (kestrel_dir / "output.vcf").write_text("vcf", encoding="utf-8")
    (kestrel_dir / "output.bam").write_text("bam", encoding="utf-8")
    (kestrel_dir / "kestrel_pre_result.tsv").write_text("pre", encoding="utf-8")

    prior_summary = {
        "schema_version": 3,
        "version": VERSION,
        "input_files": {"bam": input_bam.name},
        "sample_name": "sample",
        "reference_key_used": None,
        "decision_profile_sha256": resolve_run_configuration().decision_profile.digest,
        "pipeline_start": "2026-09-05T08:00:00.000000",
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "start": "2026-09-05T08:05:00",
                "end": "2026-09-05T08:10:00",
                "command": "run_kestrel(...)",
                "result_file": str(kestrel_tsv),
                "file_type": "tsv",
                "md5sum": "original_md5_that_does_not_match",
                "parsed_result": {"comments": [], "data": []},
            }
        ],
    }
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    harness = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        resume=True,
    )
    harness.stages["run_kestrel"].assert_called_once()
