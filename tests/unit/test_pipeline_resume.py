"""Tests for pipeline resume orchestration and output directory warnings (#20)."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Any

import pytest

from tests.support.pipeline_harness import MINIMAL_CONFIG, run_pipeline_under_harness
from vntyper.scripts import summary_steps
from vntyper.scripts.run_configuration import resolve_run_configuration
from vntyper.version import __version__ as VERSION

pytestmark = pytest.mark.unit


def _make_prior_summary(
    *,
    input_files: dict[str, str],
    canonical_input_files: dict[str, str] | None = None,
    sample_name: str = "sample",
    reference_assembly: str = "hg19",
    reference_key_used: str | None = None,
    decision_profile_sha256: str | None = None,
    extra_modules: list[str] | None = None,
    analysis_settings: dict[str, Any] | None = None,
    steps: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    resolved_modules = sorted(extra_modules or [])
    advntr_model = None
    if "advntr" in resolved_modules:
        advntr_ref = MINIMAL_CONFIG.get("reference_data", {}).get(f"advntr_reference_vntr_{reference_assembly}")
        advntr_model = str(Path(advntr_ref).resolve()) if advntr_ref else None

    default_settings = {
        "reference_assembly": reference_assembly,
        "fast_mode": False,
        "custom_regions": None,
        "bed_file": None,
        "advntr_model": advntr_model,
        "advntr_max_coverage": None,
        "bam_processing": dict(MINIMAL_CONFIG.get("bam_processing", {}))
        if "bam_processing" in MINIMAL_CONFIG
        else None,
        "extra_modules": resolved_modules,
    }
    raw_muc1 = MINIMAL_CONFIG.get("reference_data", {}).get("muc1_reference_vntr")
    kestrel_ref_path = str(Path(raw_muc1).resolve()) if raw_muc1 else None
    res = {
        "schema_version": 3,
        "version": VERSION,
        "input_files": input_files,
        "canonical_input_files": canonical_input_files,
        "sample_name": sample_name,
        "reference_assembly_requested": reference_assembly,
        "reference_key_used": reference_key_used,
        "kestrel_reference_path": kestrel_ref_path,
        "decision_profile_sha256": decision_profile_sha256 or resolve_run_configuration().decision_profile.digest,
        "pipeline_start": "2026-09-05T08:00:00.000000",
        "analysis_settings": analysis_settings if analysis_settings is not None else default_settings,
        "steps": steps or [],
    }
    if "advntr" in resolved_modules and advntr_model:
        model_path = Path(advntr_model)
        sha = hashlib.sha256(model_path.read_bytes()).hexdigest() if model_path.is_file() else "0" * 64
        res["advntr_model"] = {
            "sha256": sha,
            "schema_version": "1.0",
            "window_bp": 1000,
            "genomic_interval": "chr1:1-1000",
        }
    return res


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


def test_resume_fatal_refusal_on_input_files_mismatch(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    output_dir = tmp_path / "mismatch_dir"
    output_dir.mkdir()
    input_dir = tmp_path / "input_dir"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.touch()

    prior_summary = _make_prior_summary(
        input_files={"bam": "other_sample.bam"},
        canonical_input_files={"bam": str((input_dir / "other_sample.bam").resolve())},
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    with caplog.at_level(logging.ERROR), pytest.raises(AssertionError, match="run_pipeline exited"):
        run_pipeline_under_harness(
            output_dir=output_dir,
            bam=str(input_bam),
            resume=True,
            expect_failure=False,
        )

    assert any("Resume refused: input files differ" in record.message for record in caplog.records)


def test_resume_refusal_when_basename_matches_but_canonical_path_differs(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Refuses resume when two different inputs share the same display basename (#20 / P1)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    patient_a = tmp_path / "patient_a"
    patient_a.mkdir()
    bam_a = patient_a / "sample.bam"
    bam_a.touch()

    patient_b = tmp_path / "patient_b"
    patient_b.mkdir()
    bam_b = patient_b / "sample.bam"
    bam_b.touch()

    prior_summary = _make_prior_summary(
        input_files={"bam": "sample.bam"},
        canonical_input_files={"bam": str(bam_a.resolve())},
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    with caplog.at_level(logging.ERROR), pytest.raises(AssertionError, match="run_pipeline exited"):
        run_pipeline_under_harness(
            output_dir=output_dir,
            bam=str(bam_b),
            resume=True,
            expect_failure=False,
        )

    assert any("canonical input files differ" in record.message for record in caplog.records)


def test_resume_refusal_when_checkpoint_lacks_canonical_identities(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Refuses legacy checkpoints that cannot establish canonical identity (#20 / P1)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    input_bam = tmp_path / "sample.bam"
    input_bam.touch()

    prior_summary = _make_prior_summary(
        input_files={"bam": "sample.bam"},
        canonical_input_files=None,  # Legacy checkpoint
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    with caplog.at_level(logging.ERROR), pytest.raises(AssertionError, match="run_pipeline exited"):
        run_pipeline_under_harness(
            output_dir=output_dir,
            bam=str(input_bam),
            resume=True,
            expect_failure=False,
        )

    assert any("checkpoint lacks canonical input identities" in record.message for record in caplog.records)


def test_resume_refusal_on_analysis_settings_mismatch(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Refuses resume when result-affecting settings differ from prior run (#20 / P1)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    input_bam = tmp_path / "sample.bam"
    input_bam.touch()

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        analysis_settings={
            "reference_assembly": "hg19",
            "fast_mode": True,  # Prior ran fast_mode=True
            "custom_regions": None,
            "bed_file": None,
            "advntr_model": None,
            "extra_modules": [],
        },
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    with caplog.at_level(logging.ERROR), pytest.raises(AssertionError, match="run_pipeline exited"):
        run_pipeline_under_harness(
            output_dir=output_dir,
            bam=str(input_bam),
            resume=True,
            fast_mode=False,  # Current runs fast_mode=False
            expect_failure=False,
        )

    assert any("analysis setting 'fast_mode' differs" in record.message for record in caplog.records)


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

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        extra_modules=["advntr"],
        steps=[
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
    )
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

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        steps=[
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
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    harness = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        resume=True,
    )
    harness.stages["run_kestrel"].assert_called_once()


def test_resume_skips_bam_conversion_when_reusable(tmp_path: Path) -> None:
    output_dir = tmp_path / "resume_conv_dir"
    output_dir.mkdir()
    input_dir = tmp_path / "input_dir"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.touch()

    fastq_dir = output_dir / "fastq_bam_processing"
    fastq_dir.mkdir()
    r1 = fastq_dir / "output_R1.fastq.gz"
    r2 = fastq_dir / "output_R2.fastq.gz"
    other = fastq_dir / "output_other.fastq.gz"
    single = fastq_dir / "output_single.fastq.gz"
    sliced = fastq_dir / "output_sliced.bam"
    r1.write_bytes(b"r1")
    r2.write_bytes(b"r2")
    other.write_bytes(b"")
    single.write_bytes(b"")
    sliced.write_bytes(b"sliced")

    r1_md5 = hashlib.md5(b"r1").hexdigest()

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_BAM_TO_FASTQ,
                "start": "2026-09-05T08:01:00",
                "end": "2026-09-05T08:05:00",
                "command": "process_bam_to_fastq(...)",
                "result_file": str(r1),
                "file_type": "fastq",
                "md5sum": r1_md5,
            }
        ],
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    harness = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        resume=True,
    )

    harness.stages["process_bam_to_fastq"].assert_not_called()


def test_resume_skips_fastq_alignment_when_reusable(tmp_path: Path) -> None:
    output_dir = tmp_path / "resume_fastq_dir"
    output_dir.mkdir()
    input_dir = tmp_path / "input_dir"
    input_dir.mkdir()
    f1 = input_dir / "in_R1.fastq.gz"
    f2 = input_dir / "in_R2.fastq.gz"
    f1.touch()
    f2.touch()

    fastq_dir = output_dir / "fastq_bam_processing"
    fastq_dir.mkdir()
    r1 = fastq_dir / "output_R1.fastq.gz"
    r2 = fastq_dir / "output_R2.fastq.gz"
    other = fastq_dir / "output_other.fastq.gz"
    single = fastq_dir / "output_single.fastq.gz"
    sliced = fastq_dir / "output_sliced.bam"
    r1.write_bytes(b"r1")
    r2.write_bytes(b"r2")
    other.write_bytes(b"")
    single.write_bytes(b"")
    sliced.write_bytes(b"sliced")

    align_dir = output_dir / "alignment_processing"
    align_dir.mkdir()
    sorted_bam = align_dir / "output_sorted.bam"
    sorted_bam.write_bytes(b"sorted_bam")
    bam_md5 = hashlib.md5(b"sorted_bam").hexdigest()

    r1_md5 = hashlib.md5(b"r1").hexdigest()

    prior_summary = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_FASTQ_ALIGNMENT,
                "start": "2026-09-05T08:01:00",
                "end": "2026-09-05T08:03:00",
                "command": "align_and_sort_fastq(...)",
                "result_file": str(sorted_bam),
                "file_type": "bam",
                "md5sum": bam_md5,
            },
            {
                "step": summary_steps.STEP_BAM_TO_FASTQ_POST_ALIGNMENT,
                "start": "2026-09-05T08:03:00",
                "end": "2026-09-05T08:05:00",
                "command": "process_bam_to_fastq(...)",
                "result_file": str(r1),
                "file_type": "fastq",
                "md5sum": r1_md5,
            },
        ],
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    harness = run_pipeline_under_harness(
        output_dir=output_dir,
        fastq1=str(f1),
        fastq2=str(f2),
        resume=True,
    )

    harness.stages["align_and_sort_fastq"].assert_not_called()


def test_resume_restores_complete_routed_fastq_set_including_other_and_single(
    tmp_path: Path,
) -> None:
    """Restores complete routed FASTQ tuple including other and single outputs (#20 / P1)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.touch()

    fastq_dir = output_dir / "fastq_bam_processing"
    fastq_dir.mkdir()
    r1 = fastq_dir / "output_R1.fastq.gz"
    r2 = fastq_dir / "output_R2.fastq.gz"
    other = fastq_dir / "output_other.fastq.gz"
    single = fastq_dir / "output_single.fastq.gz"
    sliced = fastq_dir / "output_sliced.bam"
    r1.write_bytes(b"r1")
    r2.write_bytes(b"r2")
    other.write_bytes(b"other")
    single.write_bytes(b"single")
    sliced.write_bytes(b"sliced")

    r1_md5 = hashlib.md5(b"r1").hexdigest()

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_BAM_TO_FASTQ,
                "result_file": str(r1),
                "file_type": "fastq",
                "md5sum": r1_md5,
            }
        ],
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    harness = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        resume=True,
        stage_side_effects={
            "route_converted_fastqs": lambda paths, config: (str(r1), str(r2), str(other), str(single))
        },
    )
    harness.stages["process_bam_to_fastq"].assert_not_called()
    assert harness.stages["run_kestrel"].call_count == 1
    call_kwargs = harness.stages["run_kestrel"].call_args.kwargs
    assert len(call_kwargs["fastq_files"]) == 4


def test_resume_reconciles_callers_when_advntr_is_reused(tmp_path: Path) -> None:
    """Ensures cross-caller reconciliation runs when adVNTR results are reused (#20 / P1)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.touch()

    advntr_dir = output_dir / "advntr"
    advntr_dir.mkdir()
    advntr_tsv = advntr_dir / "output_adVNTR_result.tsv"
    advntr_tsv.write_text("advntr data", encoding="utf-8")
    advntr_md5 = hashlib.md5(b"advntr data").hexdigest()

    kestrel_dir = output_dir / "kestrel"
    kestrel_dir.mkdir()
    kestrel_tsv = kestrel_dir / "kestrel_result.tsv"
    kestrel_tsv.write_text("kestrel fresh data", encoding="utf-8")

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        extra_modules=["advntr"],
        steps=[
            {
                "step": summary_steps.STEP_ADVNTR,
                "result_file": str(advntr_tsv),
                "file_type": "tsv",
                "md5sum": advntr_md5,
            }
        ],
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    h = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        extra_modules=["advntr"],
        resume=True,
    )
    h.stages["run_advntr"].assert_not_called()
    h.stages["reconcile_caller_outputs"].assert_called_once()


def test_resume_preserves_conversion_checkpoint_without_overwriting_by_fastq_qc(
    tmp_path: Path,
) -> None:
    """Preserves conversion checkpoints before rerunning FASTQ QC (#20 / P2)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    f1 = tmp_path / "in_R1.fastq.gz"
    f2 = tmp_path / "in_R2.fastq.gz"
    f1.touch()
    f2.touch()

    fastq_dir = output_dir / "fastq_bam_processing"
    fastq_dir.mkdir()
    r1 = fastq_dir / "output_R1.fastq.gz"
    r2 = fastq_dir / "output_R2.fastq.gz"
    other = fastq_dir / "output_other.fastq.gz"
    single = fastq_dir / "output_single.fastq.gz"
    sliced = fastq_dir / "output_sliced.bam"
    r1.write_bytes(b"extracted reads")
    r2.write_bytes(b"extracted reads 2")
    other.write_bytes(b"")
    single.write_bytes(b"")
    sliced.write_bytes(b"sliced")

    align_dir = output_dir / "alignment_processing"
    align_dir.mkdir()
    sorted_bam = align_dir / "output_sorted.bam"
    sorted_bam.write_bytes(b"sorted_bam")
    bam_md5 = hashlib.md5(b"sorted_bam").hexdigest()
    r1_md5 = hashlib.md5(b"extracted reads").hexdigest()

    prior_summary = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_FASTQ_ALIGNMENT,
                "result_file": str(sorted_bam),
                "file_type": "bam",
                "md5sum": bam_md5,
            },
            {
                "step": summary_steps.STEP_BAM_TO_FASTQ_POST_ALIGNMENT,
                "result_file": str(r1),
                "file_type": "fastq",
                "md5sum": r1_md5,
            },
        ],
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    harness = run_pipeline_under_harness(
        output_dir=output_dir,
        fastq1=str(f1),
        fastq2=str(f2),
        resume=True,
    )
    harness.stages["process_fastq"].assert_not_called()
    harness.stages["align_and_sort_fastq"].assert_not_called()
    assert r1.read_bytes() == b"extracted reads"


def test_resume_validates_alignment_artifact_before_skipping_alignment(tmp_path: Path) -> None:
    """Reruns alignment when alignment artifact is missing or corrupted (#20 / P2)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    f1 = tmp_path / "in_R1.fastq.gz"
    f2 = tmp_path / "in_R2.fastq.gz"
    f1.touch()
    f2.touch()

    fastq_dir = output_dir / "fastq_bam_processing"
    fastq_dir.mkdir()
    r1 = fastq_dir / "output_R1.fastq.gz"
    r2 = fastq_dir / "output_R2.fastq.gz"
    sliced = fastq_dir / "output_sliced.bam"
    r1.write_bytes(b"extracted reads")
    r2.write_bytes(b"extracted reads 2")
    sliced.write_bytes(b"sliced")

    align_dir = output_dir / "alignment_processing"
    align_dir.mkdir()
    sorted_bam = align_dir / "output_sorted.bam"
    sorted_bam.write_bytes(b"corrupted_bam")

    r1_md5 = hashlib.md5(b"extracted reads").hexdigest()

    prior_summary = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_FASTQ_ALIGNMENT,
                "result_file": str(sorted_bam),
                "file_type": "bam",
                "md5sum": "recorded_md5_that_differs",
            },
            {
                "step": summary_steps.STEP_BAM_TO_FASTQ_POST_ALIGNMENT,
                "result_file": str(r1),
                "file_type": "fastq",
                "md5sum": r1_md5,
            },
        ],
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    harness = run_pipeline_under_harness(
        output_dir=output_dir,
        fastq1=str(f1),
        fastq2=str(f2),
        resume=True,
    )
    harness.stages["align_and_sort_fastq"].assert_called_once()


def test_incompatible_resume_refuses_before_replacing_provenance_snapshots(
    tmp_path: Path,
) -> None:
    """Incompatible resume rejects before touching provenance snapshots (#20 / P2)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    input_bam = tmp_path / "sample.bam"
    input_bam.touch()

    profile_snapshot = output_dir / "decision_profile.json"
    profile_snapshot.write_text("pre-existing profile snapshot", encoding="utf-8")

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
    )
    prior_summary["version"] = "1.0.0"
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    with pytest.raises(AssertionError, match="run_pipeline exited"):
        run_pipeline_under_harness(
            output_dir=output_dir,
            bam=str(input_bam),
            resume=True,
            expect_failure=False,
        )

    assert profile_snapshot.read_text(encoding="utf-8") == "pre-existing profile snapshot"


def test_resume_prepares_alignment_plan_when_reusing_sorted_bam(tmp_path: Path) -> None:
    """Prepares an alignment plan via preflight when reusing only sorted BAM (#20 / P1)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    f1 = input_dir / "r1.fastq.gz"
    f2 = input_dir / "r2.fastq.gz"
    f1.touch()
    f2.touch()

    align_dir = output_dir / "alignment_processing"
    align_dir.mkdir()
    sorted_bam = align_dir / "output_sorted.bam"
    sorted_bam.write_bytes(b"bam content")
    bam_md5 = hashlib.md5(b"bam content").hexdigest()

    prior_summary = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_FASTQ_ALIGNMENT,
                "result_file": str(sorted_bam),
                "file_type": "bam",
                "md5sum": bam_md5,
            }
        ],
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    harness = run_pipeline_under_harness(
        output_dir=output_dir,
        fastq1=str(f1),
        fastq2=str(f2),
        resume=True,
    )

    harness.stages["align_and_sort_fastq"].assert_not_called()
    assert harness.stages["run_preflight"].call_count >= 1
    plan = harness.stages["process_bam_to_fastq"].call_args.kwargs["plan"]
    assert plan is not None
    assert plan.input_path == str(sorted_bam)


def test_resume_forces_qc_rerun_when_alignment_must_rerun(tmp_path: Path) -> None:
    """Forces fastp QC to rerun when alignment must rerun (#20 / P1)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    f1 = input_dir / "r1.fastq.gz"
    f2 = input_dir / "r2.fastq.gz"
    f1.touch()
    f2.touch()

    fastq_dir = output_dir / "fastq_bam_processing"
    fastq_dir.mkdir()
    qc_json = fastq_dir / "output.json"
    qc_json.write_text("{}", encoding="utf-8")
    qc_md5 = hashlib.md5(b"{}").hexdigest()

    # QC is recorded, but alignment is absent
    prior_summary = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_FASTQ_QC,
                "result_file": str(qc_json),
                "file_type": "json",
                "md5sum": qc_md5,
            }
        ],
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    harness = run_pipeline_under_harness(
        output_dir=output_dir,
        fastq1=str(f1),
        fastq2=str(f2),
        resume=True,
    )

    harness.stages["process_fastq"].assert_called_once()
    harness.stages["align_and_sort_fastq"].assert_called_once()


def test_resume_refusal_on_advntr_max_coverage_mismatch(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Refuses resume when advntr_max_coverage differs from prior run (#20 / P1)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.touch()

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        extra_modules=["advntr"],
        analysis_settings={
            "reference_assembly": "hg19",
            "fast_mode": False,
            "custom_regions": None,
            "bed_file": None,
            "advntr_model": None,
            "advntr_max_coverage": 500,  # Prior ran with max_coverage=500
            "bam_processing": None,
            "extra_modules": ["advntr"],
        },
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    with caplog.at_level(logging.ERROR), pytest.raises(AssertionError, match="run_pipeline exited"):
        run_pipeline_under_harness(
            output_dir=output_dir,
            bam=str(input_bam),
            extra_modules=["advntr"],
            run_pipeline_kwargs={"module_args": {"advntr": {"max_coverage": 1000}}},
            resume=True,
            expect_failure=False,
        )

    assert any("analysis setting 'advntr_max_coverage' differs" in record.message for record in caplog.records)


def test_resume_refuses_reused_stage_when_sibling_artifact_md5_differs(tmp_path: Path) -> None:
    """Rejects reuse when recorded sibling artifact checksum differs (#20 / P2)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.touch()

    kestrel_dir = output_dir / "kestrel"
    kestrel_dir.mkdir()
    kestrel_tsv = kestrel_dir / "kestrel_result.tsv"
    kestrel_tsv.write_text("kestrel tsv data", encoding="utf-8")
    (kestrel_dir / "output.vcf").write_text("vcf", encoding="utf-8")
    output_bam = kestrel_dir / "output.bam"
    output_bam.write_text("corrupted bam content", encoding="utf-8")
    (kestrel_dir / "kestrel_pre_result.tsv").write_text("pre", encoding="utf-8")

    kestrel_md5 = hashlib.md5(b"kestrel tsv data").hexdigest()

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(kestrel_tsv),
                "file_type": "tsv",
                "md5sum": kestrel_md5,
            }
        ],
    )
    prior_summary["stage_artifact_md5s"] = {
        summary_steps.STEP_KESTREL: {
            "output.bam": "original_bam_md5_that_differs",
            "output.vcf": hashlib.md5(b"vcf").hexdigest(),
            "kestrel_pre_result.tsv": hashlib.md5(b"pre").hexdigest(),
        }
    }
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    harness = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        resume=True,
    )

    # Kestrel stage must have rerun because output.bam checksum differed
    harness.stages["run_kestrel"].assert_called_once()


def test_resume_refuses_when_bwa_reference_changes(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Changing bwa_reference to a different file path refuses resume and requires rerun."""
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    f1 = in_dir / "in_R1.fastq.gz"
    f2 = in_dir / "in_R2.fastq.gz"
    f1.write_bytes(b"fq1")
    f2.write_bytes(b"fq2")

    prior_summary = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
    )
    prior_summary["reference_path"] = "/original/ref.fa"
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    with caplog.at_level(logging.ERROR), pytest.raises(AssertionError, match="run_pipeline exited"):
        run_pipeline_under_harness(
            output_dir=output_dir,
            fastq1=str(f1),
            fastq2=str(f2),
            bwa_reference="/different/ref.fa",
            resume=True,
            expect_failure=False,
        )

    assert any("Resume refused: reference path differs" in record.message for record in caplog.records)


def test_resume_preserves_sibling_checksums_across_consecutive_resumes(tmp_path: Path) -> None:
    """A resumed run preserves stage_artifact_md5s so subsequent resumes continue to verify sibling integrity."""
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    input_bam = in_dir / "sample.bam"
    input_bam.write_bytes(b"bam content")

    kestrel_dir = output_dir / "kestrel"
    kestrel_dir.mkdir()
    kestrel_tsv = kestrel_dir / "kestrel_result.tsv"
    kestrel_tsv.write_text("kestrel tsv data", encoding="utf-8")
    output_bam = kestrel_dir / "output.bam"
    output_bam.write_text("bam data", encoding="utf-8")
    output_vcf = kestrel_dir / "output.vcf"
    output_vcf.write_text("vcf data", encoding="utf-8")
    pre_tsv = kestrel_dir / "kestrel_pre_result.tsv"
    pre_tsv.write_text("pre data", encoding="utf-8")

    tsv_md5 = hashlib.md5(b"kestrel tsv data").hexdigest()
    bam_md5 = hashlib.md5(b"bam data").hexdigest()
    vcf_md5 = hashlib.md5(b"vcf data").hexdigest()
    pre_md5 = hashlib.md5(b"pre data").hexdigest()

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(kestrel_tsv),
                "file_type": "tsv",
                "md5sum": tsv_md5,
            }
        ],
    )
    prior_summary["stage_artifact_md5s"] = {
        summary_steps.STEP_KESTREL: {
            "output.bam": bam_md5,
            "output.vcf": vcf_md5,
            "kestrel_pre_result.tsv": pre_md5,
        }
    }
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    # First resume: reuses Kestrel
    harness1 = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        resume=True,
    )
    harness1.stages["run_kestrel"].assert_not_called()

    # Verify that the new summary has preserved stage_artifact_md5s
    resumed_summary = json.loads((output_dir / "pipeline_summary.json").read_text(encoding="utf-8"))
    assert summary_steps.STEP_KESTREL in resumed_summary.get("stage_artifact_md5s", {})
    assert resumed_summary["stage_artifact_md5s"][summary_steps.STEP_KESTREL]["output.bam"] == bam_md5

    # Tamper with output.bam and perform a second resume
    output_bam.write_text("corrupted after resume", encoding="utf-8")
    harness2 = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        resume=True,
    )
    # The second resume must detect tampered sibling and rerun Kestrel!
    harness2.stages["run_kestrel"].assert_called_once()


def test_resume_reruns_advntr_when_model_digest_changes(tmp_path: Path) -> None:
    """adVNTR step is re-run if prior checkpoint model digest does not match active model (#20)."""
    output_dir = tmp_path / "resume_advntr_dir"
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

    advntr_dir = output_dir / "advntr"
    advntr_dir.mkdir()
    advntr_tsv = advntr_dir / "output_adVNTR_result.tsv"
    advntr_tsv.write_text("advntr data", encoding="utf-8")

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        extra_modules=["advntr"],
        steps=[
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(kestrel_tsv),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"kestrel tsv data").hexdigest(),
            },
            {
                "step": summary_steps.STEP_ADVNTR,
                "result_file": str(advntr_tsv),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"advntr data").hexdigest(),
            },
        ],
    )
    # Tamper prior summary model sha256
    prior_summary["advntr_model"] = {"sha256": "0" * 64}
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    harness = run_pipeline_under_harness(
        output_dir=output_dir,
        extra_modules=["advntr"],
        bam=str(input_bam),
        resume=True,
    )
    # adVNTR must be rerun because model digest differed!
    harness.stages["run_advntr"].assert_called_once()
    # Kestrel should still be skipped
    harness.stages["run_kestrel"].assert_not_called()


def test_resume_reruns_kestrel_when_reference_path_changes(tmp_path: Path) -> None:
    """Kestrel is re-run when muc1_reference_vntr path changes across runs (#20)."""
    output_dir = tmp_path / "resume_kestrel_ref_dir"
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

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(kestrel_tsv),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"kestrel tsv data").hexdigest(),
            }
        ],
    )
    # Set a different reference path in prior summary
    prior_summary["kestrel_reference_path"] = "/different/path/muc1.fa"
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    harness = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        resume=True,
    )
    # Kestrel must be rerun because reference path differed!
    harness.stages["run_kestrel"].assert_called_once()


def test_resume_preserves_fastq_conversion_sibling_hashes(tmp_path: Path) -> None:
    """Reused FASTQ conversion stage preserves sibling checksums in new checkpoint (#20)."""
    output_dir = tmp_path / "resume_fq_hashes_dir"
    output_dir.mkdir()
    fq_dir = output_dir / "fastq_bam_processing"
    fq_dir.mkdir()

    r1 = fq_dir / "output_R1.fastq.gz"
    r2 = fq_dir / "output_R2.fastq.gz"
    single = fq_dir / "output_single.fastq.gz"
    other = fq_dir / "output_other.fastq.gz"
    sliced = fq_dir / "output_sliced.bam"
    r1.write_bytes(b"r1_data")
    r2.write_bytes(b"r2_data")
    single.write_bytes(b"single_data")
    other.write_bytes(b"other_data")
    sliced.write_bytes(b"sliced_data")

    input_dir = tmp_path / "input_dir"
    input_dir.mkdir()
    in_f1 = input_dir / "in_1.fastq.gz"
    in_f2 = input_dir / "in_2.fastq.gz"
    in_f1.touch()
    in_f2.touch()

    prior_summary = _make_prior_summary(
        input_files={"fastq1": in_f1.name, "fastq2": in_f2.name},
        canonical_input_files={"fastq1": str(in_f1.resolve()), "fastq2": str(in_f2.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_FASTQ_QC,
                "result_file": str(output_dir / "fastp.json"),
                "file_type": "json",
            },
            {
                "step": summary_steps.STEP_FASTQ_ALIGNMENT,
                "result_file": str(output_dir / "alignment_processing" / "output_sorted.bam"),
                "file_type": "bam",
            },
            {
                "step": summary_steps.STEP_BAM_TO_FASTQ_POST_ALIGNMENT,
                "result_file": str(r1),
                "file_type": "fastq",
                "md5sum": hashlib.md5(b"r1_data").hexdigest(),
            },
        ],
    )
    prior_summary["stage_artifact_md5s"] = {
        summary_steps.STEP_BAM_TO_FASTQ_POST_ALIGNMENT: {
            "output_R2.fastq.gz": hashlib.md5(b"r2_data").hexdigest(),
            "output_single.fastq.gz": hashlib.md5(b"single_data").hexdigest(),
            "output_other.fastq.gz": hashlib.md5(b"other_data").hexdigest(),
            "output_sliced.bam": hashlib.md5(b"sliced_data").hexdigest(),
        }
    }
    align_dir = output_dir / "alignment_processing"
    align_dir.mkdir()
    (align_dir / "output_sorted.bam").touch()
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    run_pipeline_under_harness(
        output_dir=output_dir,
        fastq1=str(in_f1),
        fastq2=str(in_f2),
        reference_assembly="hg19",
        resume=True,
    )

    new_summary = json.loads((output_dir / "pipeline_summary.json").read_text(encoding="utf-8"))
    assert summary_steps.STEP_BAM_TO_FASTQ_POST_ALIGNMENT in new_summary.get("stage_artifact_md5s", {})
    assert (
        new_summary["stage_artifact_md5s"][summary_steps.STEP_BAM_TO_FASTQ_POST_ALIGNMENT]["output_R2.fastq.gz"]
        == hashlib.md5(b"r2_data").hexdigest()
    )
