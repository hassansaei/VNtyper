"""Tests for pure resume planning and compatibility evaluation (#20)."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import pytest

from vntyper.scripts import summary_steps
from vntyper.scripts.pipeline_resume_planning import (
    ResumeCompatibility,
    build_analysis_settings,
    build_canonical_inputs_and_fingerprints,
    evaluate_resume_compatibility,
    initial_stage_carry_forward,
    record_reused_stage,
    resolve_effective_advntr_runtime,
    resolve_effective_kestrel_runtime,
    resolve_effective_shark_runtime,
)
from vntyper.scripts.resume import fingerprint_file, fingerprint_runtime
from vntyper.scripts.run_configuration import resolve_run_configuration

pytestmark = pytest.mark.unit


def test_build_analysis_settings(tmp_path: Path) -> None:
    """Build canonical analysis settings dictionary."""
    bed = tmp_path / "regions.bed"
    bed.touch()
    adv_model = tmp_path / "model.db"
    adv_model.touch()

    settings = build_analysis_settings(
        reference_assembly="hg19",
        fast_mode=True,
        custom_regions="chr1:100-200",
        bed_file=str(bed),
        advntr_reference=str(adv_model),
        module_args={"advntr": {"max_coverage": 100}},
        config={"bam_processing": {"threads": 4}},
        extra_modules=("advntr", "shark"),
    )

    assert settings["reference_assembly"] == "hg19"
    assert settings["fast_mode"] is True
    assert settings["custom_regions"] == "chr1:100-200"
    assert settings["bed_file"] == str(bed.resolve())
    assert settings["advntr_model"] == str(adv_model.resolve())
    assert settings["advntr_max_coverage"] == 100
    assert settings["bam_processing"] == {"threads": 4}
    assert settings["extra_modules"] == ["advntr", "shark"]
    assert settings["preprocessing_tools"] is None


def test_build_analysis_settings_fastq_preprocessing_tools(tmp_path: Path) -> None:
    """Build canonical analysis settings captures FASTQ preprocessing tools and fingerprints."""
    fastp_bin = tmp_path / "fastp"
    fastp_bin.write_bytes(b"fastp executable")
    bwa_bin = tmp_path / "bwa"
    bwa_bin.write_bytes(b"bwa executable")

    config = {
        "tools": {
            "fastp": str(fastp_bin),
            "bwa": str(bwa_bin),
        }
    }
    settings = build_analysis_settings(
        reference_assembly="hg19",
        fast_mode=False,
        custom_regions=None,
        bed_file=None,
        advntr_reference=None,
        module_args=None,
        config=config,
        extra_modules=(),
        input_type="FASTQ",
    )
    tools = settings["preprocessing_tools"]
    assert tools is not None
    assert tools["fastp"]["executable"] == str(fastp_bin.resolve())
    assert tools["fastp"]["fingerprint"] == fingerprint_file(fastp_bin)
    assert tools["bwa"]["executable"] == str(bwa_bin.resolve())
    assert tools["bwa"]["fingerprint"] == fingerprint_file(bwa_bin)


def test_build_canonical_inputs_and_fingerprints_bam(tmp_path: Path) -> None:
    """Resolve BAM inputs and compute fingerprint."""
    bam = tmp_path / "test.bam"
    bam.write_text("bam content", encoding="utf-8")
    bed = tmp_path / "test.bed"
    bed.write_text("bed content", encoding="utf-8")

    can, fps = build_canonical_inputs_and_fingerprints("BAM", None, None, str(bam), None, str(bed))
    assert can == {"bam": str(bam.resolve())}
    assert fps["bam"] == fingerprint_file(bam)
    assert fps["bed_file"] == fingerprint_file(bed)


def test_build_canonical_inputs_and_fingerprints_cram(tmp_path: Path) -> None:
    """Resolve CRAM inputs and compute fingerprint."""
    cram = tmp_path / "test.cram"
    cram.write_text("cram content", encoding="utf-8")

    can, fps = build_canonical_inputs_and_fingerprints("CRAM", None, None, None, str(cram), None)
    assert can == {"cram": str(cram.resolve())}
    assert fps["cram"] == fingerprint_file(cram)


def test_build_canonical_inputs_and_fingerprints_fastq(tmp_path: Path) -> None:
    """Resolve paired FASTQ inputs and compute fingerprints."""
    f1 = tmp_path / "r1.fq.gz"
    f1.write_text("r1", encoding="utf-8")
    f2 = tmp_path / "r2.fq.gz"
    f2.write_text("r2", encoding="utf-8")

    can, fps = build_canonical_inputs_and_fingerprints("FASTQ", str(f1), str(f2), None, None, None)
    assert can == {"fastq1": str(f1.resolve()), "fastq2": str(f2.resolve())}
    assert fps["fastq1"] == fingerprint_file(f1)
    assert fps["fastq2"] == fingerprint_file(f2)


def test_resolve_effective_kestrel_runtime(tmp_path: Path) -> None:
    """Resolve effective Kestrel runtime including counting mode and executable fingerprint."""
    jar = tmp_path / "kestrel.jar"
    jar.write_text("jar content", encoding="utf-8")

    config = {
        "tools": {
            "kestrel": str(jar),
            "kanalyze": "/path/to/kanalyze",
        }
    }
    mode, runtime, fp = resolve_effective_kestrel_runtime(
        resolve_run_configuration(),
        config,
        str(tmp_path),
    )
    assert mode == "split"
    assert runtime["kestrel_executable"] == str(jar.resolve())
    assert runtime["kestrel_executable_fingerprint"] == fingerprint_file(jar)
    assert fp == fingerprint_runtime(runtime)


def test_resolve_effective_kestrel_runtime_kanalyze_content_change(tmp_path: Path) -> None:
    """Replacing kanalyze.jar content in place changes Kestrel runtime fingerprint (#20)."""
    kestrel_jar = tmp_path / "kestrel.jar"
    kestrel_jar.write_text("kestrel v1", encoding="utf-8")
    kanalyze_jar = tmp_path / "kanalyze.jar"
    kanalyze_jar.write_text("kanalyze v1", encoding="utf-8")

    config = {
        "tools": {
            "kestrel": str(kestrel_jar),
            "kanalyze": str(kanalyze_jar),
        }
    }
    mode1, rt1, fp1 = resolve_effective_kestrel_runtime(
        resolve_run_configuration(),
        config,
        str(tmp_path),
    )
    assert rt1["kanalyze_fingerprint"] == fingerprint_file(kanalyze_jar)

    # Modify kanalyze.jar in place
    kanalyze_jar.write_text("kanalyze v2 (modified)", encoding="utf-8")
    mode2, rt2, fp2 = resolve_effective_kestrel_runtime(
        resolve_run_configuration(),
        config,
        str(tmp_path),
    )
    assert rt2["kanalyze_fingerprint"] == fingerprint_file(kanalyze_jar)
    assert fp2 != fp1


def test_evaluate_resume_compatibility_matches() -> None:
    """Evaluate compatibility when all prior values match."""
    prior = {
        "kestrel_reference_path": "/ref/kestrel.fa",
        "kestrel_reference_fingerprint": "k_ref_fp",
        "kestrel_motifs_path": "/ref/motifs.fa",
        "kestrel_motifs_fingerprint": "k_mot_fp",
        "kestrel_runtime_fingerprint": "k_rt_fp",
        "kestrel_counting_mode": "split",
        "advntr_model": {"sha256": "adv_sha"},
        "advntr_rus_path": "/ref/rus.txt",
        "advntr_rus_fingerprint": "adv_rus_fp",
        "advntr_runtime_fingerprint": "adv_rt_fp",
        "shark_reference_path": "/ref/shark.fa",
        "shark_reference_fingerprint": "shark_fp",
        "shark_runtime_fingerprint": "shark_rt_fp",
        "reference_fingerprint": "ref_fp",
        "reference_path": "/ref/hg19.fa",
    }

    compat = evaluate_resume_compatibility(
        prior,
        input_type="FASTQ",
        kestrel_reference_path="/ref/kestrel.fa",
        kestrel_reference_fingerprint="k_ref_fp",
        kestrel_motifs_path="/ref/motifs.fa",
        kestrel_motifs_fingerprint="k_mot_fp",
        kestrel_runtime_fingerprint="k_rt_fp",
        kestrel_counting_mode="split",
        advntr_model_sha="adv_sha",
        advntr_rus_path="/ref/rus.txt",
        advntr_rus_fingerprint="adv_rus_fp",
        advntr_runtime_fingerprint="adv_rt_fp",
        shark_reference_path="/ref/shark.fa",
        shark_reference_fingerprint="shark_fp",
        shark_runtime_fingerprint="shark_rt_fp",
        effective_reference_path="/ref/hg19.fa",
        effective_reference_fingerprint="ref_fp",
    )

    assert compat.kestrel_ref_matches is True
    assert compat.advntr_model_matches is True
    assert compat.shark_ref_matches is True
    assert compat.bwa_ref_matches is True
    assert compat.inval_align is False
    assert compat.inval_cram is False


def test_evaluate_resume_compatibility_cram_reference_mismatch() -> None:
    """Evaluate compatibility when CRAM reference mismatches, invalidating callers."""
    prior = {
        "kestrel_runtime_fingerprint": "k_rt_fp",
        "reference_fingerprint": "ref_fp_old",
        "reference_path": "/ref/old.fa",
    }

    compat = evaluate_resume_compatibility(
        prior,
        input_type="CRAM",
        kestrel_reference_path=None,
        kestrel_reference_fingerprint=None,
        kestrel_motifs_path=None,
        kestrel_motifs_fingerprint=None,
        kestrel_runtime_fingerprint="k_rt_fp",
        kestrel_counting_mode=None,
        advntr_model_sha=None,
        advntr_rus_path=None,
        advntr_rus_fingerprint=None,
        advntr_runtime_fingerprint=None,
        shark_reference_path=None,
        shark_reference_fingerprint=None,
        shark_runtime_fingerprint=None,
        effective_reference_path="/ref/new.fa",
        effective_reference_fingerprint="ref_fp_new",
    )

    assert compat.cram_ref_matches is False
    assert compat.inval_cram is True
    assert compat.kestrel_ref_matches is False
    assert compat.advntr_model_matches is False


def test_record_reused_stage() -> None:
    """Record reused stage into summary and carry forward artifact md5s."""
    summary: dict[str, Any] = {"steps": []}
    prior_summary: dict[str, Any] = {
        "pipeline_start": "2026-09-05T08:00:00",
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": "/path/to/kestrel.tsv",
                "file_type": "tsv",
                "md5sum": "1234",
            }
        ],
        "stage_artifact_md5s": {
            summary_steps.STEP_KESTREL: {"kestrel.tsv": "1234"},
        },
    }

    record_reused_stage(summary, prior_summary, summary_steps.STEP_KESTREL)
    assert len(summary["steps"]) == 1
    assert summary["steps"][0]["step"] == summary_steps.STEP_KESTREL
    assert summary["steps"][0]["reused_from"] == "2026-09-05T08:00:00"
    assert summary["stage_artifact_md5s"][summary_steps.STEP_KESTREL] == {"kestrel.tsv": "1234"}


def test_initial_stage_carry_forward(tmp_path: Path) -> None:
    """Carry forward reusable stages and skip invalidated stages."""
    out = tmp_path / "out"
    kd = out / "kestrel"
    kd.mkdir(parents=True)
    k_tsv = kd / "kestrel_result.tsv"
    k_tsv.write_text("data", encoding="utf-8")
    for name in ("output.vcf", "output_indel.vcf", "output.bam", "output.bam.bai", "kestrel_pre_result.tsv"):
        (kd / name).write_text("data", encoding="utf-8")

    summary: dict[str, Any] = {"steps": []}
    prior_summary: dict[str, Any] = {
        "pipeline_start": "2026-09-05T08:00:00",
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(k_tsv),
                "file_type": "tsv",
                "md5sum": hashlib.md5(b"data").hexdigest(),
            },
            {
                "step": summary_steps.STEP_FASTQ_ALIGNMENT,
                "result_file": str(out / "missing.bam"),
                "file_type": "bam",
            },
        ],
        "stage_artifact_md5s": {
            summary_steps.STEP_KESTREL: {"kestrel_result.tsv": hashlib.md5(b"data").hexdigest()},
        },
    }

    compat = ResumeCompatibility(
        kestrel_ref_matches=True,
        advntr_model_matches=False,
        shark_ref_matches=True,
        bwa_ref_matches=True,
        cram_ref_matches=True,
        inval_align=False,
        inval_cram=False,
    )

    initial_stage_carry_forward(
        summary,
        prior_summary,
        out,
        compatibility=compat,
        needs_advntr=False,
        extra_modules=(),
    )

    # Kestrel is reusable and carried forward; alignment result file is missing so it's not carried forward
    reused_step_names = [s["step"] for s in summary["steps"]]
    assert summary_steps.STEP_KESTREL in reused_step_names
    assert summary_steps.STEP_FASTQ_ALIGNMENT not in reused_step_names
    assert summary["stage_artifact_md5s"][summary_steps.STEP_KESTREL] == {
        "kestrel_result.tsv": hashlib.md5(b"data").hexdigest()
    }


def test_resolve_effective_advntr_runtime(tmp_path: Path) -> None:
    """Resolve effective adVNTR runtime including command, file fingerprint, and version."""
    bin_file = tmp_path / "advntr"
    bin_file.write_text("binary content", encoding="utf-8")

    run_config = resolve_run_configuration()
    config1 = {"tools": {"advntr": str(bin_file)}}
    rt1, fp1 = resolve_effective_advntr_runtime(run_config, config1, advntr_version="2.0.4")
    assert rt1["advntr_command"] == str(bin_file)
    assert rt1["advntr_command_fingerprint"] == fingerprint_file(bin_file)
    assert rt1["advntr_version"] == "2.0.4"
    assert fp1 == fingerprint_runtime(rt1)

    # Changing tool command changes fingerprint
    config2 = {"tools": {"advntr": "/different/advntr"}}
    rt2, fp2 = resolve_effective_advntr_runtime(run_config, config2, advntr_version="2.0.4")
    assert fp2 != fp1

    # Changing version changes fingerprint
    rt3, fp3 = resolve_effective_advntr_runtime(run_config, config1, advntr_version=(2, 1, 0))
    assert rt3["advntr_version"] == "2.1.0"
    assert fp3 != fp1


def test_evaluate_resume_compatibility_advntr_version_mismatch() -> None:
    """adVNTR step is not reusable when tool version differs from prior summary."""
    prior = {
        "advntr_runtime_fingerprint": "adv_rt_fp",
        "tool_versions": {"advntr": "2.0.4"},
    }
    compat = evaluate_resume_compatibility(
        prior,
        input_type="BAM",
        kestrel_reference_path=None,
        kestrel_reference_fingerprint=None,
        kestrel_motifs_path=None,
        kestrel_motifs_fingerprint=None,
        kestrel_runtime_fingerprint=None,
        kestrel_counting_mode=None,
        advntr_model_sha=None,
        advntr_rus_path=None,
        advntr_rus_fingerprint=None,
        advntr_runtime_fingerprint="adv_rt_fp",
        shark_reference_path=None,
        shark_reference_fingerprint=None,
        shark_runtime_fingerprint=None,
        effective_reference_path=None,
        effective_reference_fingerprint=None,
        advntr_version="2.1.0",
    )
    assert compat.advntr_model_matches is False


def test_evaluate_resume_compatibility_preprocessing_tools() -> None:
    """Changing fastp invalidates QC and alignment; changing bwa invalidates alignment (#20)."""
    prior_tools = {
        "fastp": {"command": "fastp", "executable": "/bin/fastp", "fingerprint": "fp_fastp_1"},
        "bwa": {"command": "bwa", "executable": "/bin/bwa", "fingerprint": "fp_bwa_1"},
    }
    prior = {
        "analysis_settings": {"preprocessing_tools": prior_tools},
    }
    base_kwargs: dict[str, Any] = {
        "input_type": "FASTQ",
        "kestrel_reference_path": None,
        "kestrel_reference_fingerprint": None,
        "kestrel_motifs_path": None,
        "kestrel_motifs_fingerprint": None,
        "kestrel_runtime_fingerprint": None,
        "kestrel_counting_mode": None,
        "advntr_model_sha": None,
        "advntr_rus_path": None,
        "advntr_rus_fingerprint": None,
        "advntr_runtime_fingerprint": None,
        "shark_reference_path": None,
        "shark_reference_fingerprint": None,
        "shark_runtime_fingerprint": None,
        "effective_reference_path": None,
        "effective_reference_fingerprint": None,
    }

    # All tools match
    compat_match = evaluate_resume_compatibility(
        prior,
        current_preprocessing_tools=prior_tools,
        **base_kwargs,
    )
    assert compat_match.inval_qc is False
    assert compat_match.inval_align is False

    # fastp changes -> inval_qc and inval_align
    fastp_changed = {
        "fastp": {"command": "fastp", "executable": "/bin/fastp", "fingerprint": "fp_fastp_2"},
        "bwa": prior_tools["bwa"],
    }
    compat_fastp = evaluate_resume_compatibility(
        prior,
        current_preprocessing_tools=fastp_changed,
        **base_kwargs,
    )
    assert compat_fastp.inval_qc is True
    assert compat_fastp.inval_align is True
    assert compat_fastp.kestrel_ref_matches is False

    # bwa changes -> inval_align, but inval_qc is False
    bwa_changed = {
        "fastp": prior_tools["fastp"],
        "bwa": {"command": "bwa", "executable": "/bin/bwa", "fingerprint": "fp_bwa_2"},
    }
    compat_bwa = evaluate_resume_compatibility(
        prior,
        current_preprocessing_tools=bwa_changed,
        **base_kwargs,
    )
    assert compat_bwa.inval_qc is False
    assert compat_bwa.inval_align is True
    assert compat_bwa.kestrel_ref_matches is False


def test_resolve_effective_shark_runtime() -> None:
    """Resolve effective SHARK runtime captures command, executable, and fingerprint."""
    import copy

    from tests.unit.test_pipeline_resume import MINIMAL_CONFIG

    rc = resolve_run_configuration()
    cfg = copy.deepcopy(MINIMAL_CONFIG)
    cfg["tools"]["shark"] = "different-shark-v2"

    rt, fp = resolve_effective_shark_runtime(rc, cfg)
    assert rt["shark_command"] == "different-shark-v2"
    assert "shark_executable" in rt
    assert "shark_executable_fingerprint" in rt
    assert fp is not None

    # Changing shark command alters the fingerprint
    cfg2 = copy.deepcopy(MINIMAL_CONFIG)
    cfg2["tools"]["shark"] = "shark-v3"
    _, fp2 = resolve_effective_shark_runtime(rc, cfg2)
    assert fp2 != fp


def test_evaluate_resume_compatibility_shark_tool_mismatch() -> None:
    """Changing shark tool invalidates FASTQ QC, alignment, and downstream callers."""
    prior_tools = {
        "fastp": {"command": "fastp", "executable": "/bin/fastp", "fingerprint": "fp_fastp_1"},
        "bwa": {"command": "bwa", "executable": "/bin/bwa", "fingerprint": "fp_bwa_1"},
        "shark": {"command": "shark_v1", "executable": "/bin/shark_v1", "fingerprint": "fp_shark_1"},
    }
    prior = {
        "analysis_settings": {"preprocessing_tools": prior_tools},
    }
    base_kwargs: dict[str, Any] = {
        "input_type": "FASTQ",
        "kestrel_reference_path": None,
        "kestrel_reference_fingerprint": None,
        "kestrel_motifs_path": None,
        "kestrel_motifs_fingerprint": None,
        "kestrel_runtime_fingerprint": None,
        "kestrel_counting_mode": None,
        "advntr_model_sha": None,
        "advntr_rus_path": None,
        "advntr_rus_fingerprint": None,
        "advntr_runtime_fingerprint": None,
        "shark_reference_path": None,
        "shark_reference_fingerprint": None,
        "shark_runtime_fingerprint": None,
        "effective_reference_path": None,
        "effective_reference_fingerprint": None,
    }

    shark_changed = {
        "fastp": prior_tools["fastp"],
        "bwa": prior_tools["bwa"],
        "shark": {"command": "shark_v2", "executable": "/bin/shark_v2", "fingerprint": "fp_shark_2"},
    }
    compat_shark = evaluate_resume_compatibility(
        prior,
        current_preprocessing_tools=shark_changed,
        **base_kwargs,
    )
    assert compat_shark.inval_qc is True
    assert compat_shark.inval_align is True
    assert compat_shark.kestrel_ref_matches is False
