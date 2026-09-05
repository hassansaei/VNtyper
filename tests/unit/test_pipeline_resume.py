"""Tests for pipeline resume orchestration and output directory warnings (#20)."""

from __future__ import annotations

import copy
import hashlib
import json
import logging
from collections.abc import Mapping
from pathlib import Path
from typing import Any, cast
from unittest import mock

import pytest

from tests.support.pipeline_harness import MINIMAL_CONFIG, run_pipeline_under_harness
from vntyper.scripts import summary_steps
from vntyper.scripts.kestrel_counting import DEFAULT_KANALYZE_PATH
from vntyper.scripts.resume import fingerprint_file, fingerprint_runtime
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
    reference_fingerprint: str | None = None,
    shark_reference_path: str | None = None,
    shark_reference_fingerprint: str | None = None,
    shark_runtime_fingerprint: str | None = None,
    kestrel_runtime_fingerprint: str | None = None,
    kestrel_counting_mode: str | None = None,
) -> dict[str, Any]:
    resolved_modules = sorted(extra_modules or [])
    advntr_model = None
    if "advntr" in resolved_modules:
        advntr_ref = MINIMAL_CONFIG.get("reference_data", {}).get(f"advntr_reference_vntr_{reference_assembly}")
        advntr_model = str(Path(advntr_ref).resolve()) if advntr_ref else None

    project_root = str(Path(__file__).resolve().parent.parent.parent)
    from vntyper.scripts.pipeline_resume_planning import (
        resolve_effective_advntr_runtime,
        resolve_effective_kestrel_runtime,
        resolve_effective_shark_runtime,
    )

    default_mode, _, resolved_kestrel_fp = resolve_effective_kestrel_runtime(
        resolve_run_configuration(),
        MINIMAL_CONFIG,
        project_root,
    )
    default_kestrel_rt_fp: str | None = resolved_kestrel_fp
    default_counting_mode = kestrel_counting_mode or default_mode
    if kestrel_counting_mode and kestrel_counting_mode != default_mode:
        import os

        raw_kestrel = MINIMAL_CONFIG.get("tools", {}).get("kestrel")
        kestrel_path = (
            str(Path(os.path.join(project_root, raw_kestrel)).resolve())
            if raw_kestrel and not os.path.isabs(raw_kestrel)
            else (str(Path(raw_kestrel).resolve()) if raw_kestrel else None)
        )
        kestrel_fp = fingerprint_file(Path(kestrel_path)) if kestrel_path and Path(kestrel_path).is_file() else None
        default_kestrel_rt_fp = fingerprint_runtime(
            {
                **dict(resolve_run_configuration().kestrel_runtime),
                "kestrel_counting_mode": kestrel_counting_mode,
                "kanalyze": MINIMAL_CONFIG.get("tools", {}).get("kanalyze", DEFAULT_KANALYZE_PATH),
                "kestrel_executable": kestrel_path,
                "kestrel_executable_fingerprint": kestrel_fp,
            }
        )

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
    input_fingerprints = None
    if canonical_input_files:
        from vntyper.scripts.resume import fingerprint_file

        input_fingerprints = {
            k: fingerprint_file(Path(v)) for k, v in canonical_input_files.items() if Path(v).is_file()
        }
    bed_candidate = (analysis_settings and analysis_settings.get("bed_file")) or default_settings.get("bed_file")
    if bed_candidate and Path(str(bed_candidate)).is_file():
        if input_fingerprints is None:
            input_fingerprints = {}
        from vntyper.scripts.resume import fingerprint_file

        input_fingerprints["bed_file"] = fingerprint_file(Path(str(bed_candidate)))
    res = {
        "schema_version": 3,
        "version": VERSION,
        "input_files": input_files,
        "canonical_input_files": canonical_input_files,
        "input_fingerprints": input_fingerprints,
        "sample_name": sample_name,
        "reference_assembly_requested": reference_assembly,
        "reference_key_used": reference_key_used,
        "kestrel_reference_path": kestrel_ref_path,
        "decision_profile_sha256": decision_profile_sha256 or resolve_run_configuration().decision_profile.digest,
        "kestrel_runtime_fingerprint": (
            kestrel_runtime_fingerprint if kestrel_runtime_fingerprint is not None else default_kestrel_rt_fp
        ),
        "kestrel_counting_mode": default_counting_mode,
        "advntr_runtime_fingerprint": (
            resolve_effective_advntr_runtime(
                resolve_run_configuration(),
                MINIMAL_CONFIG,
                advntr_version="2.0.4",
            )[1]
            if "advntr" in resolved_modules
            else None
        ),
        "reference_fingerprint": reference_fingerprint,
        "shark_reference_path": shark_reference_path,
        "shark_runtime_fingerprint": shark_runtime_fingerprint
        if shark_runtime_fingerprint is not None
        else (
            resolve_effective_shark_runtime(resolve_run_configuration(), MINIMAL_CONFIG)[1]
            if "shark" in resolved_modules
            else None
        ),
        "pipeline_start": "2026-09-05T08:00:00.000000",
        "analysis_settings": analysis_settings if analysis_settings is not None else default_settings,
        "steps": steps or [],
    }
    if "advntr" in resolved_modules:
        res["tool_versions"] = {"advntr": "2.0.4"}
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
    (kestrel_dir / "output_indel.vcf").write_text("indel vcf", encoding="utf-8")
    (kestrel_dir / "output.bam").write_text("bam", encoding="utf-8")
    (kestrel_dir / "output.bam.bai").write_text("bai", encoding="utf-8")
    (kestrel_dir / "output.bed").write_text("bed", encoding="utf-8")
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
    (kestrel_dir / "output_indel.vcf").write_text("indel vcf", encoding="utf-8")
    (kestrel_dir / "output.bam").write_text("bam", encoding="utf-8")
    (kestrel_dir / "output.bam.bai").write_text("bai", encoding="utf-8")
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

    def overwrite_qc(f1, f2, threads, output, name, config):
        for mate in ("R1", "R2"):
            (Path(output) / f"{name}_{mate}.fastq.gz").write_bytes(b"whole input QC reads")

    harness = run_pipeline_under_harness(
        output_dir=output_dir,
        fastq1=str(f1),
        fastq2=str(f2),
        resume=True,
        stage_side_effects={"process_fastq": overwrite_qc},
    )
    harness.stages["process_fastq"].assert_called_once()
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
    (kestrel_dir / "output_indel.vcf").write_text("indel vcf", encoding="utf-8")
    output_bam = kestrel_dir / "output.bam"
    output_bam.write_text("corrupted bam content", encoding="utf-8")
    (kestrel_dir / "output.bam.bai").write_text("bai", encoding="utf-8")
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
            "output.bam.bai": hashlib.md5(b"bai").hexdigest(),
            "output.vcf": hashlib.md5(b"vcf").hexdigest(),
            "output_indel.vcf": hashlib.md5(b"indel vcf").hexdigest(),
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
    output_bai = kestrel_dir / "output.bam.bai"
    output_bai.write_text("bai data", encoding="utf-8")
    output_vcf = kestrel_dir / "output.vcf"
    output_vcf.write_text("vcf data", encoding="utf-8")
    indel_vcf = kestrel_dir / "output_indel.vcf"
    indel_vcf.write_text("indel data", encoding="utf-8")
    pre_tsv = kestrel_dir / "kestrel_pre_result.tsv"
    pre_tsv.write_text("pre data", encoding="utf-8")

    tsv_md5 = hashlib.md5(b"kestrel tsv data").hexdigest()
    bam_md5 = hashlib.md5(b"bam data").hexdigest()
    bai_md5 = hashlib.md5(b"bai data").hexdigest()
    vcf_md5 = hashlib.md5(b"vcf data").hexdigest()
    indel_md5 = hashlib.md5(b"indel data").hexdigest()
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
            "output.bam.bai": bai_md5,
            "output.vcf": vcf_md5,
            "output_indel.vcf": indel_md5,
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
    (kestrel_dir / "output_indel.vcf").write_text("indel vcf", encoding="utf-8")
    (kestrel_dir / "output.bam").write_text("bam", encoding="utf-8")
    (kestrel_dir / "output.bam.bai").write_text("bai", encoding="utf-8")
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
    (kestrel_dir / "output_indel.vcf").write_text("indel vcf", encoding="utf-8")
    (kestrel_dir / "output.bam").write_text("bam", encoding="utf-8")
    (kestrel_dir / "output.bam.bai").write_text("bai", encoding="utf-8")
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


def test_resume_refuses_when_input_file_content_is_modified_in_place(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Refuses resume when input file path is unchanged but content was modified (#20)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    input_dir = tmp_path / "input_dir"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.write_bytes(b"original alignment contents")

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    # Modify the input file in place
    input_bam.write_bytes(b"replaced different alignment contents")

    with caplog.at_level(logging.ERROR), pytest.raises(AssertionError, match="run_pipeline exited"):
        run_pipeline_under_harness(
            output_dir=output_dir,
            bam=str(input_bam),
            resume=True,
            expect_failure=False,
        )

    assert any("input file contents differ from checkpoint" in record.message for record in caplog.records)


def test_resume_reruns_kestrel_when_identity_aware_and_replay_missing_or_corrupted(tmp_path: Path) -> None:
    """Identity-aware positive Kestrel result cannot be reused without valid replay artifact (#20)."""
    output_dir = tmp_path / "kestrel_replay_dir"
    output_dir.mkdir()
    input_dir = tmp_path / "input_dir"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.touch()

    kestrel_dir = output_dir / "kestrel"
    kestrel_dir.mkdir()
    kestrel_tsv = kestrel_dir / "kestrel_result.tsv"
    tsv_content = "Motif\tVariant\t__Identity_Molecular_Identity\n1\tins\tmol1\n"
    kestrel_tsv.write_text(tsv_content, encoding="utf-8")
    (kestrel_dir / "output.vcf").write_text("vcf", encoding="utf-8")
    (kestrel_dir / "output_indel.vcf").write_text("indel vcf", encoding="utf-8")
    (kestrel_dir / "output.bam").write_text("bam", encoding="utf-8")
    (kestrel_dir / "output.bam.bai").write_text("bai", encoding="utf-8")
    (kestrel_dir / "kestrel_pre_result.tsv").write_text("pre", encoding="utf-8")

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(kestrel_tsv),
                "file_type": "tsv",
                "md5sum": hashlib.md5(tsv_content.encode("utf-8")).hexdigest(),
                "parsed_result": {
                    "comments": [],
                    "data": [{"Motif": "1", "__Identity_Molecular_Identity": "mol1", "Variant": "ins"}],
                },
            }
        ],
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    # Missing replay file: Kestrel must be rerun!
    harness1 = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        resume=True,
    )
    harness1.stages["run_kestrel"].assert_called_once()

    # Now write a corrupted replay file
    (kestrel_dir / "bam_identity_replay.v1.json").write_text("invalid json content", encoding="utf-8")
    harness2 = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        resume=True,
    )
    harness2.stages["run_kestrel"].assert_called_once()


def test_resume_revokes_stale_published_reports_and_exports(tmp_path: Path) -> None:
    """Resume execution immediately revokes stale report and export tables from donor run (#20)."""
    output_dir = tmp_path / "stale_reports_dir"
    output_dir.mkdir()
    input_dir = tmp_path / "input_dir"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.touch()

    # Stale files from previous run
    stale_report = output_dir / "summary_report.html"
    stale_report.write_text("old report", encoding="utf-8")
    stale_csv = output_dir / "pipeline_summary.csv"
    stale_csv.write_text("old csv", encoding="utf-8")

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    # Simulate an error right after preflight/resume start
    from unittest.mock import patch

    from vntyper.scripts import pipeline as pipeline_module

    with (
        patch.object(pipeline_module, "create_output_directories", side_effect=RuntimeError("simulated abort")),
        pytest.raises(AssertionError, match="run_pipeline exited"),
    ):
        run_pipeline_under_harness(
            output_dir=output_dir,
            bam=str(input_bam),
            resume=True,
        )

    # Stale report and CSV must have been revoked even though the pipeline aborted early!
    assert not stale_report.exists()
    assert not stale_csv.exists()


def test_resume_incompatible_advntr_model_preserves_existing_model_snapshot(tmp_path: Path) -> None:
    """Incompatible adVNTR model refusal on resume preserves existing model snapshot (#20)."""
    output_dir = tmp_path / "advntr_preserve_dir"
    output_dir.mkdir()
    input_dir = tmp_path / "input_dir"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.touch()

    advntr_dir = output_dir / "advntr"
    advntr_dir.mkdir()
    snapshot = advntr_dir / "advntr_model.db"
    snapshot.write_bytes(b"original preserved model snapshot")

    # Create an invalid/incompatible model source
    bad_model_dir = tmp_path / "bad_model"
    bad_model_dir.mkdir()
    bad_model = bad_model_dir / "bad_model.db"
    bad_model.write_bytes(b"incompatible model")

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        extra_modules=["advntr"],
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    custom_config = {
        **MINIMAL_CONFIG,
        "reference_data": {
            **MINIMAL_CONFIG["reference_data"],
            "advntr_reference_vntr_hg19": str(bad_model),
        },
    }

    with pytest.raises(AssertionError, match="run_pipeline exited"):
        run_pipeline_under_harness(
            output_dir=output_dir,
            extra_modules=["advntr"],
            bam=str(input_bam),
            config=custom_config,
            resume=True,
        )

    # Destination snapshot must NOT have been replaced by the bad candidate!
    assert snapshot.read_bytes() == b"original preserved model snapshot"


def test_resume_refuses_when_bed_file_content_is_modified_in_place(tmp_path: Path) -> None:
    """Modifying BED file content in-place causes resume refusal (#20 / adversarial review)."""
    output_dir = tmp_path / "bed_resume_dir"
    output_dir.mkdir()
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.write_bytes(b"bam content")
    bed = input_dir / "target.bed"
    bed.write_text("chr1\t100\t200\n", encoding="utf-8")

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        analysis_settings={"bed_file": str(bed.resolve())},
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    # In-place modification of BED contents
    bed.write_text("chr1\t300\t400\n", encoding="utf-8")

    with pytest.raises(AssertionError, match="run_pipeline exited"):
        run_pipeline_under_harness(
            output_dir=output_dir,
            bam=str(input_bam),
            bed_file=str(bed),
            resume=True,
        )


def test_resume_rerun_advntr_invalidates_prior_producer_outputs(tmp_path: Path) -> None:
    """Rerunning adVNTR invalidates prior producer outputs before execution (#20 / adversarial review)."""
    output_dir = tmp_path / "advntr_rerun_dir"
    output_dir.mkdir()
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.touch()

    fastq_dir = output_dir / "fastq_bam_processing"
    fastq_dir.mkdir()
    sliced = fastq_dir / "output_sliced.bam"
    sliced.write_bytes(b"sliced bam")

    advntr_dir = output_dir / "advntr"
    advntr_dir.mkdir()
    stale_raw_tsv = advntr_dir / "output_adVNTR.tsv"
    stale_raw_tsv.write_text("stale raw advntr tsv", encoding="utf-8")
    stale_result = advntr_dir / "output_adVNTR_result.tsv"
    stale_result.write_text("stale result", encoding="utf-8")

    prior_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        extra_modules=["advntr"],
        steps=[],
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    h = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        extra_modules=["advntr"],
        resume=True,
    )
    assert h.stages["run_advntr"].called
    assert not stale_raw_tsv.exists()


def test_resume_reruns_kestrel_when_indel_vcf_missing_or_corrupted(tmp_path: Path) -> None:
    """Kestrel step is re-run if output_indel.vcf is missing or corrupted (#20 / adversarial review)."""
    output_dir = tmp_path / "kestrel_indel_dir"
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
    (kestrel_dir / "output.bam").write_text("bam", encoding="utf-8")
    (kestrel_dir / "kestrel_pre_result.tsv").write_text("pre", encoding="utf-8")
    # Note: output_indel.vcf is intentionally omitted

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
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    h = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        resume=True,
    )
    assert h.stages["run_kestrel"].called


def test_resume_reruns_kestrel_when_reference_fingerprint_changes(tmp_path: Path) -> None:
    """Kestrel step is re-run if kestrel reference content fingerprint changes (#20)."""
    output_dir = tmp_path / "kestrel_fp_dir"
    output_dir.mkdir()
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    input_bam = input_dir / "sample.bam"
    input_bam.touch()

    ref_file = tmp_path / "muc1.fa"
    ref_file.write_text(">muc1\nACGTACGT\n", encoding="utf-8")
    config = dict(MINIMAL_CONFIG)
    config["reference_data"] = dict(MINIMAL_CONFIG.get("reference_data", {}))
    config["reference_data"]["muc1_reference_vntr"] = str(ref_file)

    kestrel_dir = output_dir / "kestrel"
    kestrel_dir.mkdir()
    kestrel_tsv = kestrel_dir / "kestrel_result.tsv"
    kestrel_tsv.write_text("kestrel tsv data", encoding="utf-8")
    (kestrel_dir / "output.vcf").write_text("vcf", encoding="utf-8")
    (kestrel_dir / "output_indel.vcf").write_text("indel vcf", encoding="utf-8")
    (kestrel_dir / "output.bam").write_text("bam", encoding="utf-8")
    (kestrel_dir / "output.bam.bai").write_text("bai", encoding="utf-8")
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
    prior_summary["kestrel_reference_path"] = str(ref_file.resolve())
    prior_summary["kestrel_reference_fingerprint"] = "0" * 64
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    h = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        config=config,
        resume=True,
    )
    assert h.stages["run_kestrel"].called


def test_resume_reruns_process_fastq_when_qc_corrupted_while_reusing_alignment(tmp_path: Path) -> None:
    """Corrupted output.json reruns process_fastq while reusing alignment (#20)."""
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
    qc_json = fastq_dir / "output.json"
    r1.write_bytes(b"extracted reads")
    r2.write_bytes(b"extracted reads 2")
    other.write_bytes(b"")
    single.write_bytes(b"")
    sliced.write_bytes(b"sliced")
    qc_json.write_text("invalid json {", encoding="utf-8")

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
                "step": summary_steps.STEP_FASTQ_QC,
                "result_file": str(qc_json),
                "file_type": "json",
                "md5sum": hashlib.md5(b"invalid json {").hexdigest(),
            },
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
    harness.stages["process_fastq"].assert_called_once()
    harness.stages["align_and_sort_fastq"].assert_not_called()
    assert r1.read_bytes() == b"extracted reads"


def test_resume_skips_process_fastq_when_qc_is_reusable(tmp_path: Path) -> None:
    """Valid output.json skips process_fastq when resuming (#20)."""
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
    qc_json = fastq_dir / "output.json"
    r1.write_bytes(b"extracted reads")
    r2.write_bytes(b"extracted reads 2")
    other.write_bytes(b"")
    single.write_bytes(b"")
    sliced.write_bytes(b"sliced")
    qc_json.write_text("{}", encoding="utf-8")

    align_dir = output_dir / "alignment_processing"
    align_dir.mkdir()
    sorted_bam = align_dir / "output_sorted.bam"
    sorted_bam.write_bytes(b"sorted_bam")
    bam_md5 = hashlib.md5(b"sorted_bam").hexdigest()
    r1_md5 = hashlib.md5(b"extracted reads").hexdigest()
    qc_md5 = hashlib.md5(b"{}").hexdigest()

    prior_summary = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_FASTQ_QC,
                "result_file": str(qc_json),
                "file_type": "json",
                "md5sum": qc_md5,
            },
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


def test_resume_preserves_caller_records_when_interrupted_via_donor(tmp_path: Path) -> None:
    """Interrupted resume preserves prior valid caller records via donor summary (#20)."""
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
    (kestrel_dir / "output_indel.vcf").write_text("indel vcf", encoding="utf-8")
    (kestrel_dir / "output.bam").write_text("bam", encoding="utf-8")
    (kestrel_dir / "output.bam.bai").write_text("bai", encoding="utf-8")
    (kestrel_dir / "kestrel_pre_result.tsv").write_text("pre", encoding="utf-8")

    # A donor summary with a valid Kestrel step exists from an interrupted run
    donor_summary = _make_prior_summary(
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
    (output_dir / "pipeline_summary.donor.json").write_text(json.dumps(donor_summary), encoding="utf-8")

    # The current pipeline_summary.json is truncated or missing Kestrel
    current_summary = _make_prior_summary(
        input_files={"bam": input_bam.name},
        canonical_input_files={"bam": str(input_bam.resolve())},
        steps=[],
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(current_summary), encoding="utf-8")

    h = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        resume=True,
    )
    # Kestrel step was restored from donor and reused without rerunning!
    h.stages["run_kestrel"].assert_not_called()
    assert not (output_dir / "pipeline_summary.donor.json").exists()


def test_interrupted_resume_with_reference_change_does_not_reuse_kestrel(tmp_path: Path) -> None:
    """Interruption before Kestrel does not carry forward invalid old-reference result (#20)."""
    output_dir = tmp_path / "resume_interrupt_dir"
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
    (kestrel_dir / "output_indel.vcf").write_text("indel vcf", encoding="utf-8")
    (kestrel_dir / "output.bam").write_text("bam", encoding="utf-8")
    (kestrel_dir / "output.bam.bai").write_text("bai", encoding="utf-8")
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
    prior_summary["kestrel_reference_path"] = "/different/path/muc1.fa"
    (output_dir / "pipeline_summary.json").write_text(json.dumps(prior_summary), encoding="utf-8")

    # First resume attempt is interrupted during coverage calculation (before Kestrel)
    h1 = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        resume=True,
        expect_failure=True,
        stage_side_effects={"calculate_vntr_coverage": RuntimeError("interrupted before Kestrel")},
    )
    assert h1.error is not None

    # Second resume attempt runs to completion
    h2 = run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(input_bam),
        resume=True,
    )
    # Kestrel MUST be executed because the reference path differed!
    h2.stages["run_kestrel"].assert_called_once()


def test_load_prior_summary_rejects_incompatible_donor_and_fresh_run_clears_it(tmp_path: Path) -> None:
    """Incompatible donor checkpoint is ignored and fresh runs clear stale donors (#20)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    bam_a = input_dir / "sample_a.bam"
    bam_b = input_dir / "sample_b.bam"
    bam_a.touch()
    bam_b.touch()

    # Create primary summary for sample_b
    primary = _make_prior_summary(
        sample_name="sample_b",
        input_files={"bam": bam_b.name},
        canonical_input_files={"bam": str(bam_b.resolve())},
        steps=[],
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(primary), encoding="utf-8")

    # Create donor summary from another sample (sample_a) with Kestrel step
    kestrel_dir = output_dir / "kestrel"
    kestrel_dir.mkdir()
    kestrel_tsv = kestrel_dir / "kestrel_result.tsv"
    kestrel_tsv.write_text("sample_a kestrel", encoding="utf-8")
    (kestrel_dir / "output.vcf").write_text("vcf", encoding="utf-8")
    (kestrel_dir / "output_indel.vcf").write_text("indel", encoding="utf-8")
    (kestrel_dir / "output.bam").write_text("bam", encoding="utf-8")
    (kestrel_dir / "output.bam.bai").write_text("bai", encoding="utf-8")
    (kestrel_dir / "kestrel_pre_result.tsv").write_text("pre", encoding="utf-8")

    donor = _make_prior_summary(
        sample_name="sample_a",
        input_files={"bam": bam_a.name},
        canonical_input_files={"bam": str(bam_a.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(kestrel_tsv),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"sample_a kestrel").hexdigest(),
            }
        ],
    )
    (output_dir / "pipeline_summary.donor.json").write_text(json.dumps(donor), encoding="utf-8")

    # Resuming sample_b must NOT merge Kestrel from sample_a's donor
    from vntyper.scripts.resume import load_prior_summary

    loaded = load_prior_summary(output_dir)
    assert loaded is not None
    assert loaded["sample_name"] == "sample_b"
    assert not any(s.get("step") == summary_steps.STEP_KESTREL for s in loaded.get("steps", []))

    # A fresh run on sample_b (resume=False) must clear the stale donor file
    run_pipeline_under_harness(
        output_dir=output_dir,
        bam=str(bam_b),
        resume=False,
    )
    assert not (output_dir / "pipeline_summary.donor.json").exists()


def test_resume_protects_patient_input_directory_before_revoking_published_reports(tmp_path: Path) -> None:
    """Ownership verification occurs before any file in output_dir is deleted on resume (#20)."""
    p = tmp_path
    bam = p / "input" / "sample.bam"
    bam.parent.mkdir(parents=True)
    bam.touch()

    prior = _make_prior_summary(
        input_files={"bam": bam.name},
        canonical_input_files={"bam": str(bam.resolve())},
        steps=[],
    )
    patient_dir = bam.parent
    (patient_dir / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")
    report = patient_dir / "summary_report.html"
    report.write_text("operator report", encoding="utf-8")

    alias = p / "alias"
    alias.symlink_to(patient_dir, target_is_directory=True)

    harness = run_pipeline_under_harness(
        output_dir=alias,
        bam=str(bam),
        resume=True,
        expect_failure=True,
    )
    assert harness.error is not None
    assert report.is_file(), "Published report was deleted before input ownership verification"
    assert report.read_text(encoding="utf-8") == "operator report"


def test_resume_cram_with_reference_fasta_succeeds_and_preserves_persistent_reference(tmp_path: Path) -> None:
    """CRAM runs record persistent reference paths so resume succeeds after binding cleanup (#20)."""
    out = tmp_path / "out"
    out.mkdir()
    cram = tmp_path / "input" / "sample.cram"
    cram.parent.mkdir(parents=True)
    cram.touch()
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")

    # First run writes persistent reference path into summary
    h1 = run_pipeline_under_harness(
        output_dir=out,
        cram=str(cram),
        reference_fasta=str(ref),
    )
    assert h1.error is None
    summary1 = json.loads((out / "pipeline_summary.json").read_text(encoding="utf-8"))
    assert summary1.get("reference_path") == str(ref.resolve())
    assert summary1.get("persistent_reference_path") == str(ref.resolve())

    # Resume run validates persistent reference and completes cleanly
    h2 = run_pipeline_under_harness(
        output_dir=out,
        cram=str(cram),
        reference_fasta=str(ref),
        resume=True,
    )
    assert h2.error is None


def test_resume_reruns_kestrel_when_additional_motifs_reference_changes(tmp_path: Path) -> None:
    """Changing configured additional motifs invalidates Kestrel checkpoint (#20)."""
    out = tmp_path / "out"
    out.mkdir()
    bam = tmp_path / "input" / "sample.bam"
    bam.parent.mkdir(parents=True)
    bam.touch()

    k_dir = out / "kestrel"
    k_dir.mkdir()
    for name in [
        "kestrel_result.tsv",
        "output.vcf",
        "output_indel.vcf",
        "output.bam",
        "output.bam.bai",
        "kestrel_pre_result.tsv",
    ]:
        (k_dir / name).write_text("x", encoding="utf-8")

    prior = _make_prior_summary(
        input_files={"bam": bam.name},
        canonical_input_files={"bam": str(bam.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(k_dir / "kestrel_result.tsv"),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"x").hexdigest(),
                "parsed_result": {"data": []},
            }
        ],
    )
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    config = copy.deepcopy(MINIMAL_CONFIG)
    new_motifs = tmp_path / "new_motifs.fa"
    new_motifs.write_text(">new\nACGT\n", encoding="utf-8")
    config["reference_data"]["muc1_motifs_rev_com"] = str(new_motifs)

    harness = run_pipeline_under_harness(
        output_dir=out,
        bam=str(bam),
        config=config,
        resume=True,
    )
    assert harness.error is None
    harness.stages["run_kestrel"].assert_called_once()


def test_resume_reruns_advntr_when_code_advntr_rus_reference_changes(tmp_path: Path) -> None:
    """Changing configured code_adVNTR_RUs reference invalidates adVNTR checkpoint (#20)."""
    out = tmp_path / "out"
    out.mkdir()
    bam = tmp_path / "input" / "sample.bam"
    bam.parent.mkdir(parents=True)
    bam.touch()

    adv_dir = out / "advntr"
    adv_dir.mkdir()
    adv_tsv = adv_dir / "advntr_genotype.tsv"
    adv_tsv.write_text("advntr result", encoding="utf-8")

    prior = _make_prior_summary(
        input_files={"bam": bam.name},
        canonical_input_files={"bam": str(bam.resolve())},
        extra_modules=["advntr"],
        steps=[
            {
                "step": summary_steps.STEP_ADVNTR,
                "result_file": str(adv_tsv),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"advntr result").hexdigest(),
                "parsed_result": {"data": []},
            }
        ],
    )
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    config = copy.deepcopy(MINIMAL_CONFIG)
    new_rus = tmp_path / "new_rus.fa"
    new_rus.write_text(">RU1\nACGT\n", encoding="utf-8")
    config["reference_data"]["code_adVNTR_RUs"] = str(new_rus)

    harness = run_pipeline_under_harness(
        output_dir=out,
        bam=str(bam),
        config=config,
        extra_modules=["advntr"],
        resume=True,
    )
    assert harness.error is None
    harness.stages["run_advntr"].assert_called_once()


def test_resume_reruns_kestrel_when_kestrel_runtime_changes(tmp_path: Path) -> None:
    """Changing Kestrel runtime settings reruns Kestrel while preserving conversion (#20)."""
    from dataclasses import replace

    out = tmp_path / "out"
    out.mkdir()
    bam = tmp_path / "input" / "sample.bam"
    bam.parent.mkdir(parents=True)
    bam.touch()

    # Pre-populate converted FASTQs
    fq_dir = out / "fastq_bam_processing"
    fq_dir.mkdir()
    for name in (
        "output_sliced.bam",
        "output_R1.fastq.gz",
        "output_R2.fastq.gz",
        "output_single.fastq.gz",
        "output_other.fastq.gz",
    ):
        (fq_dir / name).write_text("dummy", encoding="utf-8")

    kd = out / "kestrel"
    kd.mkdir()
    for name in (
        "kestrel_result.tsv",
        "output.vcf",
        "output_indel.vcf",
        "output.bam",
        "output.bam.bai",
        "kestrel_pre_result.tsv",
    ):
        (kd / name).write_text("data", encoding="utf-8")

    prior = _make_prior_summary(
        input_files={"bam": bam.name},
        canonical_input_files={"bam": str(bam.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_BAM_TO_FASTQ,
                "result_file": str(fq_dir / "output_R1.fastq.gz"),
                "file_type": "fastq",
                "md5sum": hashlib.md5(b"dummy").hexdigest(),
            },
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(kd / "kestrel_result.tsv"),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"data").hexdigest(),
            },
        ],
    )
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    rc = resolve_run_configuration()
    new_settings = dict(cast(Mapping[str, Any], rc.kestrel_runtime["kestrel_settings"]))
    new_settings["kmer_sizes"] = [31]
    rc_modified = replace(rc, kestrel_runtime={**rc.kestrel_runtime, "kestrel_settings": new_settings})

    h = run_pipeline_under_harness(output_dir=out, bam=str(bam), resume=True, run_configuration=rc_modified)
    assert h.error is None
    # Kestrel must rerun because kmer_sizes changed
    assert h.stages["run_kestrel"].call_count == 1
    # BAM to FASTQ conversion was reused
    assert h.stages["process_bam_to_fastq"].call_count == 0


def test_resume_reruns_advntr_when_advntr_runtime_changes(tmp_path: Path) -> None:
    """Changing adVNTR runtime settings reruns adVNTR (#20)."""
    from dataclasses import replace

    out = tmp_path / "out"
    out.mkdir()
    bam = tmp_path / "input" / "sample.bam"
    bam.parent.mkdir(parents=True)
    bam.touch()

    adv_dir = out / "advntr"
    adv_dir.mkdir()
    adv_tsv = adv_dir / "advntr_genotype.tsv"
    adv_tsv.write_text("advntr result", encoding="utf-8")

    prior = _make_prior_summary(
        input_files={"bam": bam.name},
        canonical_input_files={"bam": str(bam.resolve())},
        extra_modules=["advntr"],
        steps=[
            {
                "step": summary_steps.STEP_ADVNTR,
                "result_file": str(adv_tsv),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"advntr result").hexdigest(),
                "parsed_result": {"data": []},
            }
        ],
    )
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    rc = resolve_run_configuration()
    new_settings = dict(cast(Mapping[str, Any], rc.advntr_runtime["settings"]))
    new_settings["additional_commands"] = "--modified-command"
    rc_modified = replace(rc, advntr_runtime={**rc.advntr_runtime, "settings": new_settings})

    h = run_pipeline_under_harness(
        output_dir=out,
        bam=str(bam),
        extra_modules=["advntr"],
        resume=True,
        run_configuration=rc_modified,
    )
    assert h.error is None
    assert h.stages["run_advntr"].call_count == 1


def test_resume_rejects_advntr_reuse_when_prior_lacks_model_provenance(tmp_path: Path) -> None:
    """If prior checkpoint has STEP_ADVNTR but lacks advntr_model, adVNTR is rerun (#20)."""
    out = tmp_path / "out"
    out.mkdir()
    bam = tmp_path / "input" / "sample.bam"
    bam.parent.mkdir(parents=True)
    bam.touch()

    adv_dir = out / "advntr"
    adv_dir.mkdir()
    adv_tsv = adv_dir / "advntr_genotype.tsv"
    adv_tsv.write_text("advntr result", encoding="utf-8")

    prior = _make_prior_summary(
        input_files={"bam": bam.name},
        canonical_input_files={"bam": str(bam.resolve())},
        extra_modules=["advntr"],
        steps=[
            {
                "step": summary_steps.STEP_ADVNTR,
                "result_file": str(adv_tsv),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"advntr result").hexdigest(),
                "parsed_result": {"data": []},
            }
        ],
    )
    # Strip advntr_model to simulate missing model provenance
    prior.pop("advntr_model", None)
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    h = run_pipeline_under_harness(
        output_dir=out,
        bam=str(bam),
        extra_modules=["advntr"],
        resume=True,
    )
    assert h.error is None
    assert h.stages["run_advntr"].call_count == 1


def test_resume_refuses_when_bwa_reference_content_changes(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """When BWA reference content changes in-place, resume is refused (#20)."""
    out = tmp_path / "out"
    out.mkdir()
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    f1 = input_dir / "r1.fq.gz"
    f2 = input_dir / "r2.fq.gz"
    f1.touch()
    f2.touch()

    ref_file = tmp_path / "ref.fa"
    ref_file.write_text(">chr1\nACGT\n", encoding="utf-8")

    prior = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
        reference_fingerprint="fp_old_reference",
    )
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    with caplog.at_level(logging.ERROR):
        h = run_pipeline_under_harness(
            output_dir=out,
            fastq1=str(f1),
            fastq2=str(f2),
            bwa_reference=str(ref_file),
            resume=True,
            expect_failure=True,
        )
    assert h.error is not None
    assert any("Resume refused: reference content differs" in record.message for record in caplog.records)


def test_resume_refuses_when_shark_reference_changes(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """When SHARK reference changes, resume is refused (#20)."""
    out = tmp_path / "out"
    out.mkdir()
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    f1 = input_dir / "r1.fq.gz"
    f2 = input_dir / "r2.fq.gz"
    f1.touch()
    f2.touch()

    prior = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
        extra_modules=["shark"],
        shark_reference_path="/path/to/old/muc1.fa",
        shark_reference_fingerprint="fp_old_shark",
    )
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    with caplog.at_level(logging.ERROR):
        h = run_pipeline_under_harness(
            output_dir=out,
            fastq1=str(f1),
            fastq2=str(f2),
            extra_modules=["shark"],
            resume=True,
            expect_failure=True,
        )
    assert h.error is not None
    assert any("Resume refused: SHARK reference" in record.message for record in caplog.records)


def test_resume_invalidates_alignment_when_donor_has_bwa_fingerprint_mismatch(tmp_path: Path) -> None:
    """When donor checkpoint has an outdated BWA reference fingerprint, alignment is rerun (#20)."""
    out = tmp_path / "out"
    out.mkdir()
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    f1 = input_dir / "r1.fq.gz"
    f2 = input_dir / "r2.fq.gz"
    f1.touch()
    f2.touch()

    ref_file = tmp_path / "ref.fa"
    ref_file.write_text(">chr1\nACGT\n", encoding="utf-8")
    ref_fp = fingerprint_file(ref_file)

    align_dir = out / "alignment_processing"
    align_dir.mkdir()
    sorted_bam = align_dir / "output_sorted.bam"
    sorted_bam.write_bytes(b"sorted_bam")
    bam_md5 = hashlib.md5(b"sorted_bam").hexdigest()

    # Main summary has current reference_fingerprint but no steps completed
    current_summary = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
        reference_fingerprint=ref_fp,
    )
    # Donor summary has outdated reference_fingerprint and an alignment step
    donor_summary = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
        reference_fingerprint="fp_outdated_ref",
        steps=[
            {
                "step": summary_steps.STEP_FASTQ_ALIGNMENT,
                "result_file": str(sorted_bam),
                "file_type": "bam",
                "md5sum": bam_md5,
            }
        ],
    )
    (out / "pipeline_summary.json").write_text(json.dumps(current_summary), encoding="utf-8")
    (out / "pipeline_summary.donor.json").write_text(json.dumps(donor_summary), encoding="utf-8")

    h = run_pipeline_under_harness(
        output_dir=out,
        fastq1=str(f1),
        fastq2=str(f2),
        bwa_reference=str(ref_file),
        resume=True,
    )
    assert h.error is None
    # Donor's alignment was not grafted; alignment was run
    assert h.stages["align_and_sort_fastq"].call_count == 1


def test_resume_reruns_alignment_and_callers_when_shark_runtime_changes(tmp_path: Path) -> None:
    """When SHARK runtime settings change, previous alignment and caller checkpoints are invalidated (#20 / P1)."""
    from dataclasses import replace

    out = tmp_path / "out"
    out.mkdir()
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    f1 = input_dir / "r1.fq.gz"
    f2 = input_dir / "r2.fq.gz"
    f1.touch()
    f2.touch()

    ref_file = tmp_path / "ref.fa"
    ref_file.write_text(">chr1\nACGT\n", encoding="utf-8")

    align_dir = out / "alignment_processing"
    align_dir.mkdir()
    sorted_bam = align_dir / "output_sorted.bam"
    sorted_bam.write_bytes(b"sorted_bam")
    bam_md5 = hashlib.md5(b"sorted_bam").hexdigest()

    prior = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
        extra_modules=["shark"],
        steps=[
            {
                "step": summary_steps.STEP_FASTQ_ALIGNMENT,
                "result_file": str(sorted_bam),
                "file_type": "bam",
                "md5sum": bam_md5,
            }
        ],
    )
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    rc = resolve_run_configuration()
    new_shark_settings = dict(cast(Mapping[str, Any], rc.shark_runtime["shark_settings"]))
    new_shark_settings["confidence"] = 0.95
    rc_modified = replace(rc, shark_runtime={**rc.shark_runtime, "shark_settings": new_shark_settings})

    with (
        mock.patch(
            "vntyper.modules.shark.shark_filtering.run_shark_filter",
            return_value=(str(f1), str(f2)),
        ),
        mock.patch("vntyper.modules.shark.shark_filtering.write_shark_step_summary"),
    ):
        h = run_pipeline_under_harness(
            output_dir=out,
            fastq1=str(f1),
            fastq2=str(f2),
            bwa_reference=str(ref_file),
            extra_modules=["shark"],
            run_configuration=rc_modified,
            resume=True,
        )
    assert h.error is None
    # Because SHARK runtime changed, alignment is rerun rather than skipped
    assert h.stages["align_and_sort_fastq"].call_count == 1


def test_resume_reruns_kestrel_when_cram_reference_changes_even_without_explicit_cli_flag(
    tmp_path: Path,
) -> None:
    """Removing explicit CRAM reference invalidates prior conversion and callers (#20 / P2)."""
    out = tmp_path / "out"
    out.mkdir()
    inp = tmp_path / "in"
    inp.mkdir()
    cram = inp / "sample.cram"
    cram.touch()
    kd = out / "kestrel"
    kd.mkdir()
    ref = tmp_path / "old.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")

    for name in (
        "kestrel_result.tsv",
        "output.vcf",
        "output_indel.vcf",
        "output.bam",
        "output.bam.bai",
        "kestrel_pre_result.tsv",
    ):
        (kd / name).write_text("data", encoding="utf-8")

    prior = _make_prior_summary(
        input_files={"cram": cram.name},
        canonical_input_files={"cram": str(cram.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(kd / "kestrel_result.tsv"),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"data").hexdigest(),
            }
        ],
    )
    prior.update(
        reference_path=str(ref),
        persistent_reference_path=str(ref),
        reference_fingerprint=hashlib.sha256(ref.read_bytes()).hexdigest(),
    )
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    h = run_pipeline_under_harness(output_dir=out, cram=str(cram), resume=True)
    assert h.error is None
    # Because resolved CRAM reference differs from prior summary's explicit reference,
    # Kestrel is rerun rather than skipped
    assert h.stages["run_kestrel"].call_count == 1
    result = json.loads((out / "pipeline_summary.json").read_text(encoding="utf-8"))
    assert result["reference_path"] != str(ref)


def test_resume_refuses_conversion_reuse_when_advntr_active_and_sliced_bai_missing(
    tmp_path: Path,
) -> None:
    """When adVNTR is enabled, conversion reuse requires output_sliced.bam.bai (#20 / P2)."""
    from vntyper.scripts.resume import step_is_reusable

    out = tmp_path / "out"
    out.mkdir()
    fq = out / "fastq_bam_processing"
    fq.mkdir()
    for name in (
        "output_R1.fastq.gz",
        "output_R2.fastq.gz",
        "output_single.fastq.gz",
        "output_other.fastq.gz",
        "output_sliced.bam",
    ):
        (fq / name).write_bytes(b"data")

    prior = _make_prior_summary(
        input_files={"bam": "sample.bam"},
        extra_modules=["advntr"],
        steps=[
            {
                "step": summary_steps.STEP_BAM_TO_FASTQ,
                "result_file": str(fq / "output_R1.fastq.gz"),
                "file_type": "fastq",
                "md5sum": hashlib.md5(b"data").hexdigest(),
            }
        ],
    )
    # With missing BAI, conversion step is NOT reusable when adVNTR is active
    assert step_is_reusable(prior, summary_steps.STEP_BAM_TO_FASTQ, out, needs_advntr=True) is False

    # Once BAI is present, conversion step is reusable
    (fq / "output_sliced.bam.bai").write_bytes(b"bai_data")
    assert step_is_reusable(prior, summary_steps.STEP_BAM_TO_FASTQ, out, needs_advntr=True) is True


def test_resume_reruns_kestrel_when_kanalyze_path_changes_counting_mode(tmp_path: Path) -> None:
    """Disabling kanalyze switches counting mode and invalidates prior Kestrel checkpoint (#20 / P2)."""
    import copy

    out = tmp_path / "out"
    out.mkdir()
    inp = tmp_path / "in"
    inp.mkdir()
    bam = inp / "sample.bam"
    bam.touch()
    kd = out / "kestrel"
    kd.mkdir()

    for name in (
        "kestrel_result.tsv",
        "output.vcf",
        "output_indel.vcf",
        "output.bam",
        "output.bam.bai",
        "kestrel_pre_result.tsv",
    ):
        (kd / name).write_text("data", encoding="utf-8")

    prior = _make_prior_summary(
        input_files={"bam": bam.name},
        canonical_input_files={"bam": str(bam.resolve())},
        kestrel_counting_mode="split",
        steps=[
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(kd / "kestrel_result.tsv"),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"data").hexdigest(),
            }
        ],
    )
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    config = copy.deepcopy(MINIMAL_CONFIG)
    config["tools"]["kanalyze"] = ""

    h = run_pipeline_under_harness(output_dir=out, bam=str(bam), resume=True, config=config)
    assert h.error is None
    # Because counting mode changed from split to internal, Kestrel is rerun
    assert h.stages["run_kestrel"].call_count == 1
    result = json.loads((out / "pipeline_summary.json").read_text(encoding="utf-8"))
    assert result["kestrel_counting_mode"] == "internal"


def test_load_prior_summary_rejects_donor_with_mismatched_decision_profile(tmp_path: Path) -> None:
    """When a donor checkpoint has a different decision profile digest, its steps are not merged (#20)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    bam = input_dir / "sample.bam"
    bam.touch()

    # Primary summary has profile digest A
    primary = _make_prior_summary(
        sample_name="sample",
        input_files={"bam": bam.name},
        canonical_input_files={"bam": str(bam.resolve())},
        decision_profile_sha256="a" * 64,
        steps=[],
    )
    (output_dir / "pipeline_summary.json").write_text(json.dumps(primary), encoding="utf-8")

    # Donor checkpoint has profile digest B and a Kestrel step
    kestrel_dir = output_dir / "kestrel"
    kestrel_dir.mkdir()
    kestrel_tsv = kestrel_dir / "kestrel_result.tsv"
    kestrel_tsv.write_text("donor kestrel", encoding="utf-8")
    for name in ("output.vcf", "output_indel.vcf", "output.bam", "output.bam.bai", "kestrel_pre_result.tsv"):
        (kestrel_dir / name).write_text("data", encoding="utf-8")

    donor = _make_prior_summary(
        sample_name="sample",
        input_files={"bam": bam.name},
        canonical_input_files={"bam": str(bam.resolve())},
        decision_profile_sha256="b" * 64,
        steps=[
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(kestrel_tsv),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"donor kestrel").hexdigest(),
            }
        ],
    )
    (output_dir / "pipeline_summary.donor.json").write_text(json.dumps(donor), encoding="utf-8")

    from vntyper.scripts.resume import load_prior_summary

    loaded = load_prior_summary(output_dir)
    assert loaded is not None
    assert loaded["decision_profile_sha256"] == "a" * 64
    # Donor Kestrel step must NOT have been merged into the primary summary
    assert not any(s.get("step") == summary_steps.STEP_KESTREL for s in loaded.get("steps", []))


def test_resume_reruns_kestrel_when_configured_kestrel_executable_changes(tmp_path: Path) -> None:
    """When config['tools']['kestrel'] changes, Kestrel is rerun on resume (#20)."""
    import copy

    out = tmp_path / "out"
    out.mkdir()
    inp = tmp_path / "in"
    inp.mkdir()
    bam = inp / "sample.bam"
    bam.touch()
    kd = out / "kestrel"
    kd.mkdir()

    for name in (
        "kestrel_result.tsv",
        "output.vcf",
        "output_indel.vcf",
        "output.bam",
        "output.bam.bai",
        "kestrel_pre_result.tsv",
    ):
        (kd / name).write_text("data", encoding="utf-8")

    jar_v1 = tmp_path / "kestrel_v1.jar"
    jar_v1.write_text("jar v1", encoding="utf-8")
    jar_v2 = tmp_path / "kestrel_v2.jar"
    jar_v2.write_text("jar v2", encoding="utf-8")

    project_root = str(Path(__file__).resolve().parent.parent.parent)
    from vntyper.scripts.pipeline_resume_planning import resolve_effective_kestrel_runtime

    cfg_v1 = copy.deepcopy(MINIMAL_CONFIG)
    cfg_v1["tools"]["kestrel"] = str(jar_v1)
    _, _, rt_fp_v1 = resolve_effective_kestrel_runtime(resolve_run_configuration(), cfg_v1, project_root)

    prior = _make_prior_summary(
        input_files={"bam": bam.name},
        canonical_input_files={"bam": str(bam.resolve())},
        kestrel_runtime_fingerprint=rt_fp_v1,
        steps=[
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(kd / "kestrel_result.tsv"),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"data").hexdigest(),
            }
        ],
    )
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    cfg_v2 = copy.deepcopy(MINIMAL_CONFIG)
    cfg_v2["tools"]["kestrel"] = str(jar_v2)

    h = run_pipeline_under_harness(output_dir=out, bam=str(bam), resume=True, config=cfg_v2)
    assert h.error is None
    # Changing kestrel jar executable causes Kestrel to be rerun instead of reused
    assert h.stages["run_kestrel"].call_count == 1


def test_resume_reruns_kestrel_when_kanalyze_content_changes(tmp_path: Path) -> None:
    """Changing kanalyze.jar content invalidates Kestrel checkpoint on resume (#20)."""
    out = tmp_path / "out"
    out.mkdir()
    inp = tmp_path / "in"
    inp.mkdir()
    bam = inp / "sample.bam"
    bam.touch()
    kd = out / "kestrel"
    kd.mkdir()

    for name in (
        "kestrel_result.tsv",
        "output.vcf",
        "output_indel.vcf",
        "output.bam",
        "output.bam.bai",
        "kestrel_pre_result.tsv",
    ):
        (kd / name).write_text("data", encoding="utf-8")

    kestrel_jar = tmp_path / "kestrel.jar"
    kestrel_jar.write_text("kestrel content", encoding="utf-8")
    kanalyze_jar_v1 = tmp_path / "kanalyze_v1.jar"
    kanalyze_jar_v1.write_text("kanalyze v1", encoding="utf-8")
    kanalyze_jar_v2 = tmp_path / "kanalyze_v2.jar"
    kanalyze_jar_v2.write_text("kanalyze v2", encoding="utf-8")

    project_root = str(Path(__file__).resolve().parent.parent.parent)
    from vntyper.scripts.pipeline_resume_planning import resolve_effective_kestrel_runtime

    cfg_v1 = copy.deepcopy(MINIMAL_CONFIG)
    cfg_v1["tools"]["kestrel"] = str(kestrel_jar)
    cfg_v1["tools"]["kanalyze"] = str(kanalyze_jar_v1)
    _, _, rt_fp_v1 = resolve_effective_kestrel_runtime(resolve_run_configuration(), cfg_v1, project_root)

    prior = _make_prior_summary(
        input_files={"bam": bam.name},
        canonical_input_files={"bam": str(bam.resolve())},
        kestrel_runtime_fingerprint=rt_fp_v1,
        steps=[
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(kd / "kestrel_result.tsv"),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"data").hexdigest(),
            }
        ],
    )
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    cfg_v2 = copy.deepcopy(MINIMAL_CONFIG)
    cfg_v2["tools"]["kestrel"] = str(kestrel_jar)
    cfg_v2["tools"]["kanalyze"] = str(kanalyze_jar_v2)

    h = run_pipeline_under_harness(output_dir=out, bam=str(bam), resume=True, config=cfg_v2)
    assert h.error is None
    # Replacing kanalyze jar causes Kestrel to be rerun instead of reused
    assert h.stages["run_kestrel"].call_count == 1


def test_resume_reruns_qc_and_alignment_when_fastp_changes(tmp_path: Path) -> None:
    """Changing fastp tool invalidates FASTQ QC and alignment on resume (#20)."""
    out = tmp_path / "out"
    out.mkdir()
    inp = tmp_path / "in"
    inp.mkdir()
    f1 = inp / "r1.fq.gz"
    f2 = inp / "r2.fq.gz"
    f1.touch()
    f2.touch()

    fq_dir = out / "fastq_bam_processing"
    fq_dir.mkdir()
    qc_json = fq_dir / "output.json"
    qc_json.write_text("{}", encoding="utf-8")
    r1 = fq_dir / "output_R1.fastq.gz"
    r2 = fq_dir / "output_R2.fastq.gz"
    r1.touch()
    r2.touch()

    align_dir = out / "alignment_processing"
    align_dir.mkdir()
    sorted_bam = align_dir / "output_sorted.bam"
    sorted_bam.touch()

    fastp_v1 = tmp_path / "fastp_v1"
    fastp_v1.write_bytes(b"fastp_v1")
    fastp_v2 = tmp_path / "fastp_v2"
    fastp_v2.write_bytes(b"fastp_v2")
    bwa_bin = tmp_path / "bwa"
    bwa_bin.write_bytes(b"bwa")

    prior = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
        analysis_settings={
            "reference_assembly": "hg19",
            "fast_mode": False,
            "custom_regions": None,
            "bed_file": None,
            "advntr_model": None,
            "advntr_max_coverage": None,
            "bam_processing": copy.deepcopy(MINIMAL_CONFIG["bam_processing"]),
            "extra_modules": [],
            "preprocessing_tools": {
                "fastp": {"command": str(fastp_v1), "executable": str(fastp_v1.resolve()), "fingerprint": "fp_v1"},
                "bwa": {"command": str(bwa_bin), "executable": str(bwa_bin.resolve()), "fingerprint": "fp_bwa"},
            },
        },
        steps=[
            {
                "step": summary_steps.STEP_FASTQ_QC,
                "result_file": str(qc_json),
                "file_type": "json",
                "md5sum": hashlib.md5(b"{}").hexdigest(),
            },
            {
                "step": summary_steps.STEP_FASTQ_ALIGNMENT,
                "result_file": str(sorted_bam),
                "file_type": "bam",
                "md5sum": hashlib.md5(b"").hexdigest(),
            },
        ],
    )
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    cfg_v2 = copy.deepcopy(MINIMAL_CONFIG)
    cfg_v2["tools"]["fastp"] = str(fastp_v2)
    cfg_v2["tools"]["bwa"] = str(bwa_bin)

    h = run_pipeline_under_harness(
        output_dir=out,
        fastq1=str(f1),
        fastq2=str(f2),
        resume=True,
        config=cfg_v2,
    )
    assert h.error is None
    # Changing fastp reruns QC and alignment
    h.stages["process_fastq"].assert_called_once()
    h.stages["align_and_sort_fastq"].assert_called_once()


def test_resume_reruns_alignment_when_bwa_tool_changes(tmp_path: Path) -> None:
    """Changing bwa tool invalidates alignment and downstream callers on resume (#20)."""
    out = tmp_path / "out"
    out.mkdir()
    inp = tmp_path / "in"
    inp.mkdir()
    f1 = inp / "r1.fq.gz"
    f2 = inp / "r2.fq.gz"
    f1.touch()
    f2.touch()

    fq_dir = out / "fastq_bam_processing"
    fq_dir.mkdir()
    qc_json = fq_dir / "output.json"
    qc_json.write_text("{}", encoding="utf-8")
    r1 = fq_dir / "output_R1.fastq.gz"
    r2 = fq_dir / "output_R2.fastq.gz"
    r1.touch()
    r2.touch()

    align_dir = out / "alignment_processing"
    align_dir.mkdir()
    sorted_bam = align_dir / "output_sorted.bam"
    sorted_bam.touch()

    fastp_bin = tmp_path / "fastp"
    fastp_bin.write_bytes(b"fastp")
    bwa_v1 = tmp_path / "bwa_v1"
    bwa_v1.write_bytes(b"bwa_v1")
    bwa_v2 = tmp_path / "bwa_v2"
    bwa_v2.write_bytes(b"bwa_v2")

    prior = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
        analysis_settings={
            "reference_assembly": "hg19",
            "fast_mode": False,
            "custom_regions": None,
            "bed_file": None,
            "advntr_model": None,
            "advntr_max_coverage": None,
            "bam_processing": copy.deepcopy(MINIMAL_CONFIG["bam_processing"]),
            "extra_modules": [],
            "preprocessing_tools": {
                "fastp": {"command": str(fastp_bin), "executable": str(fastp_bin.resolve()), "fingerprint": "fp_fastp"},
                "bwa": {"command": str(bwa_v1), "executable": str(bwa_v1.resolve()), "fingerprint": "fp_v1"},
            },
        },
        steps=[
            {
                "step": summary_steps.STEP_FASTQ_QC,
                "result_file": str(qc_json),
                "file_type": "json",
                "md5sum": hashlib.md5(b"{}").hexdigest(),
            },
            {
                "step": summary_steps.STEP_FASTQ_ALIGNMENT,
                "result_file": str(sorted_bam),
                "file_type": "bam",
                "md5sum": hashlib.md5(b"").hexdigest(),
            },
        ],
    )
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    cfg_v2 = copy.deepcopy(MINIMAL_CONFIG)
    cfg_v2["tools"]["fastp"] = str(fastp_bin)
    cfg_v2["tools"]["bwa"] = str(bwa_v2)

    h = run_pipeline_under_harness(
        output_dir=out,
        fastq1=str(f1),
        fastq2=str(f2),
        resume=True,
        config=cfg_v2,
    )
    assert h.error is None
    # Changing bwa reruns alignment, while QC can be reused
    h.stages["align_and_sort_fastq"].assert_called_once()


def test_resume_reruns_kestrel_when_copied_results_dir_lacks_bam_while_original_exists(tmp_path: Path) -> None:
    """When a results directory is copied and output.bam is missing in the copy, Kestrel reruns (#20 / P2)."""
    orig_dir = tmp_path / "original_out"
    orig_dir.mkdir()
    inp = tmp_path / "in"
    inp.mkdir()
    bam = inp / "sample.bam"
    bam.touch()

    orig_kestrel = orig_dir / "kestrel"
    orig_kestrel.mkdir()
    for name in (
        "kestrel_result.tsv",
        "output.vcf",
        "output_indel.vcf",
        "output.bam",
        "output.bam.bai",
        "output.bed",
        "kestrel_pre_result.tsv",
    ):
        (orig_kestrel / name).write_text("data", encoding="utf-8")

    prior = _make_prior_summary(
        input_files={"bam": bam.name},
        canonical_input_files={"bam": str(bam.resolve())},
        steps=[
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(orig_kestrel / "kestrel_result.tsv"),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"data").hexdigest(),
            }
        ],
    )
    prior["stage_artifact_md5s"] = {
        summary_steps.STEP_KESTREL: {
            "output.bam": hashlib.md5(b"data").hexdigest(),
            "output.bam.bai": hashlib.md5(b"data").hexdigest(),
            "output.vcf": hashlib.md5(b"data").hexdigest(),
            "output_indel.vcf": hashlib.md5(b"data").hexdigest(),
            "output.bed": hashlib.md5(b"data").hexdigest(),
            "kestrel_pre_result.tsv": hashlib.md5(b"data").hexdigest(),
        }
    }
    (orig_dir / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    # Copy to new output directory
    copied_dir = tmp_path / "copied_out"
    copied_kestrel = copied_dir / "kestrel"
    copied_kestrel.mkdir(parents=True)
    (copied_dir / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")
    for name in (
        "kestrel_result.tsv",
        "output.vcf",
        "output_indel.vcf",
        # output.bam intentionally omitted in copied dir!
        "output.bam.bai",
        "output.bed",
        "kestrel_pre_result.tsv",
    ):
        (copied_kestrel / name).write_text("data", encoding="utf-8")

    # Original still has output.bam
    assert (orig_kestrel / "output.bam").is_file()
    assert not (copied_kestrel / "output.bam").is_file()

    h = run_pipeline_under_harness(output_dir=copied_dir, bam=str(bam), resume=True)
    assert h.error is None
    # Because output.bam is missing in copied_dir, Kestrel is rerun rather than reused
    assert h.stages["run_kestrel"].call_count == 1


def test_resume_reruns_alignment_and_callers_when_shark_tool_changes(tmp_path: Path) -> None:
    """Changing shark tool executable invalidates alignment and downstream callers on resume (#20 / P2)."""
    out = tmp_path / "out"
    out.mkdir()
    inp = tmp_path / "in"
    inp.mkdir()
    f1 = inp / "r1.fq.gz"
    f2 = inp / "r2.fq.gz"
    f1.touch()
    f2.touch()

    ref_file = tmp_path / "ref.fa"
    ref_file.write_text(">chr1\nACGT\n", encoding="utf-8")

    align_dir = out / "alignment_processing"
    align_dir.mkdir()
    sorted_bam = align_dir / "output_sorted.bam"
    sorted_bam.touch()

    shark_v1 = tmp_path / "shark_v1"
    shark_v1.write_bytes(b"shark_v1")
    shark_v2 = tmp_path / "shark_v2"
    shark_v2.write_bytes(b"shark_v2")

    cfg_v1 = copy.deepcopy(MINIMAL_CONFIG)
    cfg_v1["tools"]["shark"] = str(shark_v1)
    from vntyper.scripts.pipeline_resume_planning import resolve_effective_shark_runtime

    _, shark_fp_v1 = resolve_effective_shark_runtime(resolve_run_configuration(), cfg_v1)

    prior = _make_prior_summary(
        input_files={"fastq1": f1.name, "fastq2": f2.name},
        canonical_input_files={"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())},
        extra_modules=["shark"],
        shark_runtime_fingerprint=shark_fp_v1,
        analysis_settings={
            "reference_assembly": "hg19",
            "fast_mode": False,
            "custom_regions": None,
            "bed_file": None,
            "advntr_model": None,
            "advntr_max_coverage": None,
            "bam_processing": copy.deepcopy(MINIMAL_CONFIG["bam_processing"]),
            "extra_modules": ["shark"],
            "preprocessing_tools": {
                "fastp": {"command": "fastp", "executable": "/bin/fastp", "fingerprint": "fp_fastp"},
                "bwa": {"command": "bwa", "executable": "/bin/bwa", "fingerprint": "fp_bwa"},
                "shark": {"command": str(shark_v1), "executable": str(shark_v1.resolve()), "fingerprint": "fp_v1"},
            },
        },
        steps=[
            {
                "step": summary_steps.STEP_FASTQ_ALIGNMENT,
                "result_file": str(sorted_bam),
                "file_type": "bam",
                "md5sum": hashlib.md5(b"").hexdigest(),
            },
        ],
    )
    (out / "pipeline_summary.json").write_text(json.dumps(prior), encoding="utf-8")

    cfg_v2 = copy.deepcopy(MINIMAL_CONFIG)
    cfg_v2["tools"]["shark"] = str(shark_v2)

    with (
        mock.patch(
            "vntyper.modules.shark.shark_filtering.run_shark_filter",
            return_value=(str(f1), str(f2)),
        ),
        mock.patch("vntyper.modules.shark.shark_filtering.write_shark_step_summary"),
    ):
        h = run_pipeline_under_harness(
            output_dir=out,
            fastq1=str(f1),
            fastq2=str(f2),
            bwa_reference=str(ref_file),
            extra_modules=["shark"],
            config=cfg_v2,
            resume=True,
        )
    assert h.error is None
    # Changing shark executable causes alignment to rerun
    h.stages["align_and_sort_fastq"].assert_called_once()
