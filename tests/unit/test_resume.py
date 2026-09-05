"""Tests for resume checkpoint loading and run-identity refusal checks (#20)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

from vntyper.scripts import summary_steps
from vntyper.scripts.resume import (
    STEP_OUTPUT_SIBLINGS,
    _compute_md5,
    caller_advntr_matches,
    caller_kestrel_matches,
    caller_shark_matches,
    fingerprint_runtime,
    load_prior_summary,
    make_reused_step_record,
    reference_content_matches,
    resume_refusals,
    step_is_reusable,
)

pytestmark = pytest.mark.unit

SAMPLE_SUMMARY: dict[str, Any] = {
    "schema_version": 3,
    "version": "2.0.27",
    "input_files": {"bam": "/data/sample.bam"},
    "sample_name": "NA12878",
    "reference_key_used": "hg19",
    "decision_profile_sha256": "a" * 64,
    "pipeline_start": "2026-09-05T10:00:00.000000",
    "steps": [],
}


def test_load_prior_summary_reads_existing_json(tmp_path: Path) -> None:
    summary_file = tmp_path / "pipeline_summary.json"
    summary_file.write_text(json.dumps(SAMPLE_SUMMARY), encoding="utf-8")

    loaded = load_prior_summary(summary_file)
    assert loaded == SAMPLE_SUMMARY


def test_load_prior_summary_accepts_directory_path(tmp_path: Path) -> None:
    summary_file = tmp_path / "pipeline_summary.json"
    summary_file.write_text(json.dumps(SAMPLE_SUMMARY), encoding="utf-8")

    loaded = load_prior_summary(tmp_path)
    assert loaded == SAMPLE_SUMMARY


def test_load_prior_summary_returns_none_when_missing(tmp_path: Path) -> None:
    assert load_prior_summary(tmp_path / "nonexistent.json") is None
    assert load_prior_summary(tmp_path / "empty_dir") is None


def test_load_prior_summary_returns_none_on_corrupt_json(tmp_path: Path) -> None:
    bad_file = tmp_path / "pipeline_summary.json"
    bad_file.write_text("{not valid json", encoding="utf-8")

    assert load_prior_summary(bad_file) is None


def test_resume_refusals_empty_when_identity_matches() -> None:
    refusals = resume_refusals(
        SAMPLE_SUMMARY,
        version="2.0.27",
        input_files={"bam": "/data/sample.bam"},
        sample_name="NA12878",
        reference_key_used="hg19",
        decision_profile_sha256="a" * 64,
    )
    assert refusals == []


def test_resume_refusals_flags_version_mismatch() -> None:
    refusals = resume_refusals(
        SAMPLE_SUMMARY,
        version="2.1.0",
        input_files={"bam": "/data/sample.bam"},
        sample_name="NA12878",
        reference_key_used="hg19",
        decision_profile_sha256="a" * 64,
    )
    assert len(refusals) == 1
    assert "version differs" in refusals[0]
    assert "2.0.27" in refusals[0]
    assert "2.1.0" in refusals[0]


def test_resume_refusals_flags_input_files_mismatch() -> None:
    refusals = resume_refusals(
        SAMPLE_SUMMARY,
        version="2.0.27",
        input_files={"bam": "/data/other_sample.bam"},
        sample_name="NA12878",
        reference_key_used="hg19",
        decision_profile_sha256="a" * 64,
    )
    assert len(refusals) == 1
    assert "input files differ" in refusals[0]


def test_resume_refusals_flags_sample_name_mismatch() -> None:
    refusals = resume_refusals(
        SAMPLE_SUMMARY,
        version="2.0.27",
        input_files={"bam": "/data/sample.bam"},
        sample_name="NA12877",
        reference_key_used="hg19",
        decision_profile_sha256="a" * 64,
    )
    assert len(refusals) == 1
    assert "sample name differs" in refusals[0]


def test_resume_refusals_flags_reference_key_mismatch() -> None:
    refusals = resume_refusals(
        SAMPLE_SUMMARY,
        version="2.0.27",
        input_files={"bam": "/data/sample.bam"},
        sample_name="NA12878",
        reference_key_used="hg38",
        decision_profile_sha256="a" * 64,
    )
    assert len(refusals) == 1
    assert "reference key differs" in refusals[0]


def test_resume_refusals_flags_decision_profile_digest_mismatch() -> None:
    refusals = resume_refusals(
        SAMPLE_SUMMARY,
        version="2.0.27",
        input_files={"bam": "/data/sample.bam"},
        sample_name="NA12878",
        reference_key_used="hg19",
        decision_profile_sha256="b" * 64,
    )
    assert len(refusals) == 1
    assert "decision profile digest differs" in refusals[0]


def test_resume_refusals_collects_multiple_mismatches() -> None:
    refusals = resume_refusals(
        SAMPLE_SUMMARY,
        version="2.1.0",
        input_files={"cram": "/data/sample.cram"},
        sample_name="SAMPLE2",
        reference_key_used=None,
        decision_profile_sha256="b" * 64,
    )
    assert len(refusals) == 5


def _md5(content: bytes) -> str:
    return hashlib.md5(content).hexdigest()


def test_step_is_reusable_kestrel_success(tmp_path: Path) -> None:
    kestrel_dir = tmp_path / "kestrel"
    kestrel_dir.mkdir()
    result_tsv = kestrel_dir / "kestrel_result.tsv"
    result_tsv.write_text("dummy tsv content\n", encoding="utf-8")
    expected_md5 = _md5(b"dummy tsv content\n")

    for sibling in STEP_OUTPUT_SIBLINGS[summary_steps.STEP_KESTREL]:
        (kestrel_dir / sibling).write_text("sibling data", encoding="utf-8")

    prior = {
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(result_tsv),
                "md5sum": expected_md5,
            }
        ]
    }

    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, tmp_path) is True


def test_step_is_reusable_returns_false_when_step_absent() -> None:
    prior: dict[str, Any] = {"steps": []}
    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, Path("/tmp")) is False


def test_step_is_reusable_returns_false_when_result_file_missing(tmp_path: Path) -> None:
    prior = {
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(tmp_path / "kestrel" / "missing.tsv"),
                "md5sum": "d41d8cd98f00b204e9800998ecf8427e",
            }
        ]
    }
    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, tmp_path) is False


def test_step_is_reusable_returns_false_on_md5_mismatch(tmp_path: Path) -> None:
    kestrel_dir = tmp_path / "kestrel"
    kestrel_dir.mkdir()
    result_tsv = kestrel_dir / "kestrel_result.tsv"
    result_tsv.write_text("corrupted content", encoding="utf-8")

    for sibling in STEP_OUTPUT_SIBLINGS[summary_steps.STEP_KESTREL]:
        (kestrel_dir / sibling).write_text("sibling data", encoding="utf-8")

    prior = {
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(result_tsv),
                "md5sum": "different_md5_hash",
            }
        ]
    }
    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, tmp_path) is False


def test_step_is_reusable_returns_false_when_sibling_missing(tmp_path: Path) -> None:
    kestrel_dir = tmp_path / "kestrel"
    kestrel_dir.mkdir()
    result_tsv = kestrel_dir / "kestrel_result.tsv"
    result_tsv.write_text("dummy tsv content\n", encoding="utf-8")

    # Write only two of the three required siblings
    (kestrel_dir / "output.vcf").write_text("vcf data", encoding="utf-8")
    (kestrel_dir / "output.bam").write_text("bam data", encoding="utf-8")
    # kestrel_pre_result.tsv is missing!

    prior = {
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(result_tsv),
                "md5sum": _md5(b"dummy tsv content\n"),
            }
        ]
    }
    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, tmp_path) is False


def test_step_is_reusable_kestrel_returns_false_when_bam_index_missing(tmp_path: Path) -> None:
    kestrel_dir = tmp_path / "kestrel"
    kestrel_dir.mkdir()
    result_tsv = kestrel_dir / "kestrel_result.tsv"
    result_tsv.write_text("dummy tsv content\n", encoding="utf-8")

    for sibling in STEP_OUTPUT_SIBLINGS[summary_steps.STEP_KESTREL]:
        if sibling != "output.bam.bai":
            (kestrel_dir / sibling).write_text("sibling data", encoding="utf-8")

    prior = {
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(result_tsv),
                "md5sum": _md5(b"dummy tsv content\n"),
            }
        ]
    }
    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, tmp_path) is False


def test_step_is_reusable_kestrel_returns_false_when_bed_missing_with_variants(tmp_path: Path) -> None:
    kestrel_dir = tmp_path / "kestrel"
    kestrel_dir.mkdir()
    result_tsv = kestrel_dir / "kestrel_result.tsv"
    result_tsv.write_text("Motif\tPOS\tREF\tALT\n1\t50\tC\tT\n", encoding="utf-8")
    for sibling in STEP_OUTPUT_SIBLINGS[summary_steps.STEP_KESTREL]:
        (kestrel_dir / sibling).write_text("sibling data", encoding="utf-8")

    prior = {
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(result_tsv),
                "md5sum": _md5(b"Motif\tPOS\tREF\tALT\n1\t50\tC\tT\n"),
            }
        ]
    }
    # output.bed is not present, but variants exist in kestrel_result.tsv -> reject reuse
    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, tmp_path) is False

    # Once output.bed is present, reuse succeeds
    (kestrel_dir / "output.bed").write_text("1\t49\t50\n", encoding="utf-8")
    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, tmp_path) is True


def test_step_is_reusable_kestrel_returns_false_when_recorded_bed_deleted(tmp_path: Path) -> None:
    kestrel_dir = tmp_path / "kestrel"
    kestrel_dir.mkdir()
    result_tsv = kestrel_dir / "kestrel_result.tsv"
    result_tsv.write_text("header\n", encoding="utf-8")
    for sibling in STEP_OUTPUT_SIBLINGS[summary_steps.STEP_KESTREL]:
        (kestrel_dir / sibling).write_text("sibling data", encoding="utf-8")
    (kestrel_dir / "output.bed").write_text("bed", encoding="utf-8")

    prior = {
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(result_tsv),
                "md5sum": _md5(b"header\n"),
            }
        ],
        "stage_artifact_md5s": {
            summary_steps.STEP_KESTREL: {
                "output.bed": _md5(b"bed"),
            }
        },
    }
    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, tmp_path) is True
    (kestrel_dir / "output.bed").unlink()
    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, tmp_path) is False


def test_step_is_reusable_advntr_requires_only_result_file(tmp_path: Path) -> None:
    advntr_dir = tmp_path / "advntr"
    advntr_dir.mkdir()
    result_tsv = advntr_dir / "output_adVNTR_result.tsv"
    result_tsv.write_text("advntr data", encoding="utf-8")

    prior = {
        "steps": [
            {
                "step": summary_steps.STEP_ADVNTR,
                "result_file": str(result_tsv),
                "md5sum": _md5(b"advntr data"),
            }
        ]
    }
    assert step_is_reusable(prior, summary_steps.STEP_ADVNTR, tmp_path) is True


def test_step_is_reusable_bam_conversion_requires_sliced_bam_and_mate(tmp_path: Path) -> None:
    conv_dir = tmp_path / "fastq_bam_processing"
    conv_dir.mkdir()
    r1 = conv_dir / "output_R1.fastq.gz"
    r1.write_text("r1", encoding="utf-8")
    r2 = conv_dir / "output_R2.fastq.gz"
    r2.write_text("r2", encoding="utf-8")
    single = conv_dir / "output_single.fastq.gz"
    single.write_text("single", encoding="utf-8")
    other = conv_dir / "output_other.fastq.gz"
    other.write_text("other", encoding="utf-8")
    sliced_bam = conv_dir / "output_sliced.bam"
    sliced_bam.write_text("sliced bam", encoding="utf-8")

    prior = {
        "input_files": {"bam": "/data/in.bam"},
        "steps": [
            {
                "step": summary_steps.STEP_BAM_TO_FASTQ,
                "result_file": str(r1),
                "md5sum": _md5(b"r1"),
            }
        ],
    }
    assert step_is_reusable(prior, summary_steps.STEP_BAM_TO_FASTQ, tmp_path) is True

    # Missing sliced bam forces rerun
    sliced_bam.unlink()
    assert step_is_reusable(prior, summary_steps.STEP_BAM_TO_FASTQ, tmp_path) is False

    # Restore sliced bam, unlink mate -> forces rerun
    sliced_bam.write_text("sliced bam", encoding="utf-8")
    r2.unlink()
    assert step_is_reusable(prior, summary_steps.STEP_BAM_TO_FASTQ, tmp_path) is False


def test_make_reused_step_record_preserves_fields_and_adds_reused_from() -> None:
    original = {
        "step": summary_steps.STEP_KESTREL,
        "start": "2026-09-05T10:00:00",
        "end": "2026-09-05T10:05:00",
        "command": "run_kestrel(...)",
        "result_file": "/results/kestrel/kestrel_result.tsv",
        "file_type": "tsv",
        "md5sum": "abc",
        "parsed_result": {"data": []},
    }
    prior_start = "2026-09-05T09:59:00.000000"

    reused = make_reused_step_record(original, prior_start)

    assert reused["reused_from"] == prior_start
    for k, v in original.items():
        assert reused[k] == v


def test_compute_md5_returns_none_for_non_file(tmp_path: Path) -> None:
    assert _compute_md5(tmp_path / "non_existent") is None


def test_compute_md5_returns_none_on_oserror(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    f = tmp_path / "file.txt"
    f.touch()

    def mock_open(*args: Any, **kwargs: Any) -> Any:
        raise OSError("Permission denied")

    monkeypatch.setattr(Path, "open", mock_open)
    assert _compute_md5(f) is None


def test_load_prior_summary_returns_none_when_data_is_not_dict(tmp_path: Path) -> None:
    p = tmp_path / "pipeline_summary.json"
    p.write_text("[1, 2, 3]", encoding="utf-8")
    assert load_prior_summary(p) is None


def test_step_is_reusable_returns_false_for_unregistered_step(tmp_path: Path) -> None:
    prior: dict[str, Any] = {"steps": []}
    assert step_is_reusable(prior, "Unregistered Step", tmp_path) is False


def test_step_is_reusable_returns_false_when_result_file_missing_flag_set(tmp_path: Path) -> None:
    prior: dict[str, Any] = {
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file_missing": True,
                "result_file": str(tmp_path / "out.tsv"),
                "md5sum": "abc",
            }
        ]
    }
    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, tmp_path) is False


def test_step_is_reusable_returns_false_when_no_result_file(tmp_path: Path) -> None:
    prior: dict[str, Any] = {
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": None,
                "md5sum": "abc",
            }
        ]
    }
    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, tmp_path) is False


def test_step_is_reusable_returns_false_when_no_md5sum(tmp_path: Path) -> None:
    prior: dict[str, Any] = {
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(tmp_path / "out.tsv"),
                "md5sum": None,
            }
        ]
    }
    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, tmp_path) is False


def test_step_is_reusable_relocated_relative_to_output_root(tmp_path: Path) -> None:
    relocated_tsv = tmp_path / "kestrel_result.tsv"
    relocated_tsv.write_text("dummy", encoding="utf-8")
    expected_md5 = _md5(b"dummy")

    for sibling in STEP_OUTPUT_SIBLINGS[summary_steps.STEP_KESTREL]:
        (tmp_path / sibling).write_text("sibling data", encoding="utf-8")

    prior: dict[str, Any] = {
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": "/old/path/kestrel/kestrel_result.tsv",
                "md5sum": expected_md5,
            }
        ]
    }

    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, tmp_path) is True


def test_resume_refusals_detects_reference_path_difference() -> None:
    base = {
        "version": "2.0.28",
        "input_files": {"fastq1": "f1.fq.gz"},
        "sample_name": "s1",
        "reference_path": "/path/to/ref1.fa",
    }
    # Matching reference path
    refusals = resume_refusals(
        base,
        version="2.0.28",
        input_files={"fastq1": "f1.fq.gz"},
        sample_name="s1",
        reference_key_used=None,
        decision_profile_sha256="abc",
        reference_path="/path/to/ref1.fa",
    )
    assert not any("reference path" in r for r in refusals)

    # Differing reference path
    refusals = resume_refusals(
        base,
        version="2.0.28",
        input_files={"fastq1": "f1.fq.gz"},
        sample_name="s1",
        reference_key_used=None,
        decision_profile_sha256="abc",
        reference_path="/different/ref2.fa",
    )
    assert any("reference path differs" in r for r in refusals)


def test_step_is_reusable_verifies_all_four_conversion_fastqs(tmp_path: Path) -> None:
    conv_dir = tmp_path / "fastq_bam_processing"
    conv_dir.mkdir()
    r1 = conv_dir / "output_R1.fastq.gz"
    r1.write_text("r1", encoding="utf-8")
    r2 = conv_dir / "output_R2.fastq.gz"
    r2.write_text("r2", encoding="utf-8")
    single = conv_dir / "output_single.fastq.gz"
    single.write_text("single", encoding="utf-8")
    other = conv_dir / "output_other.fastq.gz"
    other.write_text("other", encoding="utf-8")
    sliced_bam = conv_dir / "output_sliced.bam"
    sliced_bam.write_text("sliced bam", encoding="utf-8")

    prior = {
        "input_files": {"bam": "/data/in.bam"},
        "stage_artifact_md5s": {
            summary_steps.STEP_BAM_TO_FASTQ: {
                "output_sliced.bam": _md5(b"sliced bam"),
                "output_R1.fastq.gz": _md5(b"r1"),
                "output_R2.fastq.gz": _md5(b"r2"),
                "output_single.fastq.gz": _md5(b"single"),
                "output_other.fastq.gz": _md5(b"other"),
            }
        },
        "steps": [
            {
                "step": summary_steps.STEP_BAM_TO_FASTQ,
                "result_file": str(r1),
                "md5sum": _md5(b"r1"),
            }
        ],
    }
    assert step_is_reusable(prior, summary_steps.STEP_BAM_TO_FASTQ, tmp_path) is True

    # Tamper with output_other.fastq.gz -> must refuse
    other.write_text("tampered other", encoding="utf-8")
    assert step_is_reusable(prior, summary_steps.STEP_BAM_TO_FASTQ, tmp_path) is False

    # Restore other, tamper with output_single.fastq.gz -> must refuse
    other.write_text("other", encoding="utf-8")
    single.write_text("tampered single", encoding="utf-8")
    assert step_is_reusable(prior, summary_steps.STEP_BAM_TO_FASTQ, tmp_path) is False


def test_fingerprint_runtime() -> None:
    assert fingerprint_runtime(None) is None
    fp1 = fingerprint_runtime({"kmer_sizes": [20]})
    assert fp1 is not None and len(fp1) == 64
    fp2 = fingerprint_runtime({"kmer_sizes": [31]})
    assert fp2 != fp1


def test_caller_kestrel_matches() -> None:
    assert caller_kestrel_matches(None) is True

    prior: dict[str, Any] = {
        "kestrel_reference_path": "/path/ref.fa",
        "kestrel_reference_fingerprint": "100:1",
        "kestrel_motifs_path": "/path/motifs.fa",
        "kestrel_motifs_fingerprint": "200:2",
        "kestrel_runtime_fingerprint": "a" * 64,
    }
    # Exact match
    assert (
        caller_kestrel_matches(
            prior,
            kestrel_reference_path="/path/ref.fa",
            kestrel_reference_fingerprint="100:1",
            kestrel_motifs_path="/path/motifs.fa",
            kestrel_motifs_fingerprint="200:2",
            kestrel_runtime_fingerprint="a" * 64,
        )
        is True
    )
    # Runtime fingerprint mismatch
    assert (
        caller_kestrel_matches(
            prior,
            kestrel_reference_path="/path/ref.fa",
            kestrel_reference_fingerprint="100:1",
            kestrel_motifs_path="/path/motifs.fa",
            kestrel_motifs_fingerprint="200:2",
            kestrel_runtime_fingerprint="b" * 64,
        )
        is False
    )
    # Reference path mismatch
    assert (
        caller_kestrel_matches(
            prior,
            kestrel_reference_path="/other/ref.fa",
            kestrel_reference_fingerprint="100:1",
            kestrel_motifs_path="/path/motifs.fa",
            kestrel_motifs_fingerprint="200:2",
            kestrel_runtime_fingerprint="a" * 64,
        )
        is False
    )
    # Reference fingerprint mismatch
    assert (
        caller_kestrel_matches(
            prior,
            kestrel_reference_path="/path/ref.fa",
            kestrel_reference_fingerprint="100:999",
            kestrel_motifs_path="/path/motifs.fa",
            kestrel_motifs_fingerprint="200:2",
            kestrel_runtime_fingerprint="a" * 64,
        )
        is False
    )
    # Motifs path mismatch
    assert (
        caller_kestrel_matches(
            prior,
            kestrel_reference_path="/path/ref.fa",
            kestrel_reference_fingerprint="100:1",
            kestrel_motifs_path="/other/motifs.fa",
            kestrel_motifs_fingerprint="200:2",
            kestrel_runtime_fingerprint="a" * 64,
        )
        is False
    )
    # Motifs fingerprint mismatch
    assert (
        caller_kestrel_matches(
            prior,
            kestrel_reference_path="/path/ref.fa",
            kestrel_reference_fingerprint="100:1",
            kestrel_motifs_path="/path/motifs.fa",
            kestrel_motifs_fingerprint="200:999",
            kestrel_runtime_fingerprint="a" * 64,
        )
        is False
    )


def test_caller_advntr_matches() -> None:
    assert caller_advntr_matches(None) is True

    prior: dict[str, Any] = {
        "steps": [{"step": summary_steps.STEP_ADVNTR}],
        "advntr_model": {"sha256": "c" * 64},
        "advntr_rus_path": "/path/rus.fa",
        "advntr_rus_fingerprint": "300:3",
        "advntr_runtime_fingerprint": "d" * 64,
    }
    # Exact match
    assert (
        caller_advntr_matches(
            prior,
            curr_model_sha="c" * 64,
            advntr_rus_path="/path/rus.fa",
            advntr_rus_fingerprint="300:3",
            advntr_runtime_fingerprint="d" * 64,
        )
        is True
    )
    # Model SHA mismatch
    assert (
        caller_advntr_matches(
            prior,
            curr_model_sha="e" * 64,
            advntr_rus_path="/path/rus.fa",
            advntr_rus_fingerprint="300:3",
            advntr_runtime_fingerprint="d" * 64,
        )
        is False
    )
    # Prior lacks model sha
    prior_no_sha = dict(prior, advntr_model={})
    assert (
        caller_advntr_matches(
            prior_no_sha,
            curr_model_sha="c" * 64,
            advntr_rus_path="/path/rus.fa",
            advntr_rus_fingerprint="300:3",
            advntr_runtime_fingerprint="d" * 64,
        )
        is False
    )
    # RUS path mismatch
    assert (
        caller_advntr_matches(
            prior,
            curr_model_sha="c" * 64,
            advntr_rus_path="/other/rus.fa",
            advntr_rus_fingerprint="300:3",
            advntr_runtime_fingerprint="d" * 64,
        )
        is False
    )
    # RUS fingerprint mismatch
    assert (
        caller_advntr_matches(
            prior,
            curr_model_sha="c" * 64,
            advntr_rus_path="/path/rus.fa",
            advntr_rus_fingerprint="300:999",
            advntr_runtime_fingerprint="d" * 64,
        )
        is False
    )
    # Runtime fingerprint mismatch
    assert (
        caller_advntr_matches(
            prior,
            curr_model_sha="c" * 64,
            advntr_rus_path="/path/rus.fa",
            advntr_rus_fingerprint="300:3",
            advntr_runtime_fingerprint="e" * 64,
        )
        is False
    )
    # Version mismatch
    assert (
        caller_advntr_matches(
            {**prior, "tool_versions": {"advntr": "2.0.4"}},
            curr_model_sha="c" * 64,
            advntr_rus_path="/path/rus.fa",
            advntr_rus_fingerprint="300:3",
            advntr_runtime_fingerprint="d" * 64,
            advntr_version="2.1.0",
        )
        is False
    )


def test_caller_shark_matches() -> None:
    assert caller_shark_matches(None) is True

    prior = {
        "shark_reference_path": "/path/shark.fa",
        "shark_reference_fingerprint": "100:abc",
    }
    # Matching
    assert (
        caller_shark_matches(
            prior,
            shark_reference_path="/path/shark.fa",
            shark_reference_fingerprint="100:abc",
        )
        is True
    )
    # Neither has shark
    assert caller_shark_matches({}, shark_reference_path=None, shark_reference_fingerprint=None) is True
    # Prior has shark, current does not
    assert caller_shark_matches(prior, shark_reference_path=None, shark_reference_fingerprint=None) is False
    # Prior has no shark, current does
    assert (
        caller_shark_matches(
            {},
            shark_reference_path="/path/shark.fa",
            shark_reference_fingerprint="100:abc",
        )
        is False
    )
    # Path differs
    assert (
        caller_shark_matches(
            prior,
            shark_reference_path="/other/shark.fa",
            shark_reference_fingerprint="100:abc",
        )
        is False
    )
    # Fingerprint differs
    assert (
        caller_shark_matches(
            prior,
            shark_reference_path="/path/shark.fa",
            shark_reference_fingerprint="100:xyz",
        )
        is False
    )
    # Runtime fingerprint matches
    assert (
        caller_shark_matches(
            dict(prior, shark_runtime_fingerprint="rt_fp"),
            shark_reference_path="/path/shark.fa",
            shark_reference_fingerprint="100:abc",
            shark_runtime_fingerprint="rt_fp",
        )
        is True
    )
    # Runtime fingerprint differs
    assert (
        caller_shark_matches(
            dict(prior, shark_runtime_fingerprint="rt_fp_old"),
            shark_reference_path="/path/shark.fa",
            shark_reference_fingerprint="100:abc",
            shark_runtime_fingerprint="rt_fp_new",
        )
        is False
    )


def test_fingerprint_file_detects_middle_byte_changes_in_large_files(tmp_path: Path) -> None:
    """Fingerprint detects mutations outside the initial and trailing 64 KiB windows (#20 / P1)."""
    from vntyper.scripts.resume import fingerprint_file

    large_file = tmp_path / "large.fq"
    # Create 250 KiB file
    data = bytearray(b"A" * 256000)
    large_file.write_bytes(data)
    fp_before = fingerprint_file(large_file)

    # Mutate a single byte in the middle (byte 128,000)
    data[128000] = ord(b"G")
    large_file.write_bytes(data)
    fp_after = fingerprint_file(large_file)

    assert fp_before != fp_after


def test_reference_content_matches() -> None:
    assert reference_content_matches(None) is True
    assert reference_content_matches({}, reference_fingerprint=None) is True

    prior = {"reference_fingerprint": "500:fp1"}
    assert reference_content_matches(prior, reference_fingerprint="500:fp1") is True
    assert reference_content_matches(prior, reference_fingerprint="500:fp2") is False
    assert reference_content_matches(prior, reference_fingerprint=None) is False
    assert reference_content_matches({}, reference_fingerprint="500:fp1") is False


def test_resume_refusals_reference_fingerprint_and_shark() -> None:
    prior = {
        "version": "2.0.27",
        "input_files": {"fastq_1": "r1.fq"},
        "sample_name": "s1",
        "reference_key_used": "hg19",
        "decision_profile_sha256": "p" * 64,
        "reference_fingerprint": "fp_ref_old",
        "shark_reference_path": "/path/shark.fa",
        "shark_reference_fingerprint": "fp_shark_old",
    }
    base_kwargs: dict[str, Any] = {
        "version": "2.0.27",
        "input_files": {"fastq_1": "r1.fq"},
        "sample_name": "s1",
        "reference_key_used": "hg19",
        "decision_profile_sha256": "p" * 64,
        "reference_fingerprint": "fp_ref_old",
        "shark_reference_path": "/path/shark.fa",
        "shark_reference_fingerprint": "fp_shark_old",
    }
    # All match
    assert resume_refusals(prior, **base_kwargs) == []

    # Reference content mismatch
    mismatch_ref = dict(base_kwargs, reference_fingerprint="fp_ref_new")
    refusals = resume_refusals(prior, **mismatch_ref)
    assert any("reference content differs" in r for r in refusals)

    # SHARK path mismatch
    mismatch_shark_path = dict(base_kwargs, shark_reference_path="/other/shark.fa")
    refusals = resume_refusals(prior, **mismatch_shark_path)
    assert any("SHARK reference path differs" in r for r in refusals)

    # SHARK content mismatch
    mismatch_shark_fp = dict(base_kwargs, shark_reference_fingerprint="fp_shark_new")
    refusals = resume_refusals(prior, **mismatch_shark_fp)
    assert any("SHARK reference content differs" in r for r in refusals)


def test_donor_summary_filters_steps_on_shark_or_bwa_change(tmp_path: Path) -> None:
    main_summary = tmp_path / "pipeline_summary.json"
    donor_summary = tmp_path / "pipeline_summary.donor.json"

    primary_data: dict[str, Any] = {
        "version": "2.0.27",
        "sample_name": "s1",
        "input_files": {"fastq_1": "r1.fq"},
        "canonical_input_files": None,
        "reference_assembly_requested": "hg19",
        "decision_profile": {"digest": "d" * 64},
        "reference_fingerprint": "ref_fp_current",
        "shark_reference_path": "/path/shark.fa",
        "shark_reference_fingerprint": "shark_fp_current",
        "steps": [],
    }
    donor_data: dict[str, Any] = {
        "version": "2.0.27",
        "sample_name": "s1",
        "input_files": {"fastq_1": "r1.fq"},
        "canonical_input_files": None,
        "reference_assembly_requested": "hg19",
        "decision_profile": {"digest": "d" * 64},
        "reference_fingerprint": "ref_fp_old",  # Content mismatch
        "shark_reference_path": "/path/shark.fa",
        "shark_reference_fingerprint": "shark_fp_current",
        "steps": [
            {"step": summary_steps.STEP_FASTQ_ALIGNMENT, "result_file": "/path/bam"},
            {"step": summary_steps.STEP_SHARK, "result_file": "/path/step.json"},
            {"step": summary_steps.STEP_KESTREL, "result_file": "/path/res.tsv"},
            {"step": summary_steps.STEP_FASTQ_QC, "result_file": "/path/qc.json"},
        ],
        "stage_artifact_md5s": {
            summary_steps.STEP_FASTQ_ALIGNMENT: {"out.bam": "md5"},
            summary_steps.STEP_KESTREL: {"output.vcf": "md5"},
            summary_steps.STEP_FASTQ_QC: {"qc.json": "md5"},
        },
    }

    main_summary.write_text(json.dumps(primary_data), encoding="utf-8")
    donor_summary.write_text(json.dumps(donor_data), encoding="utf-8")

    loaded = load_prior_summary(main_summary)
    assert loaded is not None
    loaded_step_names = [s.get("step") for s in loaded.get("steps", [])]
    # Invalidation should filter out FASTQ_ALIGNMENT, SHARK, KESTREL, but allow FASTQ_QC
    assert summary_steps.STEP_FASTQ_ALIGNMENT not in loaded_step_names
    assert summary_steps.STEP_SHARK not in loaded_step_names
    assert summary_steps.STEP_KESTREL not in loaded_step_names
    assert summary_steps.STEP_FASTQ_QC in loaded_step_names
    assert summary_steps.STEP_FASTQ_QC in loaded.get("stage_artifact_md5s", {})
    assert summary_steps.STEP_FASTQ_ALIGNMENT not in loaded.get("stage_artifact_md5s", {})


def test_step_is_reusable_kestrel_allows_negative_result_without_bed(tmp_path: Path) -> None:
    """A negative Kestrel result placeholder does not require output.bed for reuse (#20)."""
    from vntyper.scripts.kestrel_result_artifacts import write_empty_kestrel_artifacts

    kestrel_dir = tmp_path / "kestrel"
    kestrel_dir.mkdir()
    res_path = write_empty_kestrel_artifacts(kestrel_dir, ["## Header"])
    (kestrel_dir / "output.bam").touch()
    (kestrel_dir / "output.bam.bai").touch()
    (kestrel_dir / "output.vcf").touch()
    (kestrel_dir / "output_indel.vcf").touch()

    prior: dict[str, Any] = {
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(res_path),
                "md5sum": _md5(res_path.read_bytes()),
                "parsed_result": {
                    "data": [
                        {
                            "Motif": "None",
                            "Variant": "None",
                            "POS": "None",
                            "REF": "None",
                            "ALT": "None",
                            "Motif_sequence": "None",
                            "Estimated_Depth_AlternateVariant": "None",
                            "Estimated_Depth_Variant_ActiveRegion": "None",
                            "Depth_Score": "None",
                            "Confidence": "Negative",
                        }
                    ]
                },
            }
        ],
        "stage_artifact_md5s": {
            summary_steps.STEP_KESTREL: {
                "output.bam": _md5(b""),
                "output.bam.bai": _md5(b""),
                "output.vcf": _md5(b""),
                "output_indel.vcf": _md5(b""),
                "kestrel_pre_result.tsv": _md5((kestrel_dir / "kestrel_pre_result.tsv").read_bytes()),
            }
        },
    }
    # Reusable even though output.bed does not exist
    assert step_is_reusable(prior, summary_steps.STEP_KESTREL, tmp_path) is True

    # If output.bed is recorded in stage_artifact_md5s, it is required
    prior_with_bed_recorded = dict(prior)
    prior_with_bed_recorded["stage_artifact_md5s"] = {
        summary_steps.STEP_KESTREL: {
            **prior["stage_artifact_md5s"][summary_steps.STEP_KESTREL],
            "output.bed": "some_md5",
        }
    }
    assert step_is_reusable(prior_with_bed_recorded, summary_steps.STEP_KESTREL, tmp_path) is False

    # If result has real variants (not Negative placeholder), output.bed is required
    res_with_variant = kestrel_dir / "kestrel_result.tsv"
    res_with_variant.write_text(
        "Motif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\t"
        "Estimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\n"
        "1\tinsC\t155161000\tA\tAC\tseq\t10\t20\t0.5\tHigh\n",
        encoding="utf-8",
    )
    prior_variant: dict[str, Any] = {
        "steps": [
            {
                "step": summary_steps.STEP_KESTREL,
                "result_file": str(res_with_variant),
                "md5sum": _md5(res_with_variant.read_bytes()),
            }
        ],
        "stage_artifact_md5s": prior["stage_artifact_md5s"],
    }
    assert step_is_reusable(prior_variant, summary_steps.STEP_KESTREL, tmp_path) is False


def test_donor_summary_filters_steps_on_preprocessing_tools_change(tmp_path: Path) -> None:
    main_summary = tmp_path / "pipeline_summary.json"
    donor_summary = tmp_path / "pipeline_summary.donor.json"

    primary_data: dict[str, Any] = {
        "version": "2.0.27",
        "sample_name": "s1",
        "input_files": {"fastq_1": "r1.fq"},
        "analysis_settings": {
            "preprocessing_tools": {
                "fastp": {"command": "fastp", "fingerprint": "fp_fastp_new"},
                "bwa": {"command": "bwa", "fingerprint": "fp_bwa_current"},
            }
        },
        "reference_fingerprint": "ref_fp_current",
        "shark_reference_fingerprint": "shark_fp_current",
        "steps": [],
    }
    donor_data: dict[str, Any] = {
        "version": "2.0.27",
        "sample_name": "s1",
        "input_files": {"fastq_1": "r1.fq"},
        "analysis_settings": {
            "preprocessing_tools": {
                "fastp": {"command": "fastp", "fingerprint": "fp_fastp_old"},
                "bwa": {"command": "bwa", "fingerprint": "fp_bwa_current"},
            }
        },
        "reference_fingerprint": "ref_fp_current",
        "shark_reference_fingerprint": "shark_fp_current",
        "steps": [
            {"step": summary_steps.STEP_FASTQ_QC, "result_file": "/path/qc.json"},
            {"step": summary_steps.STEP_FASTQ_ALIGNMENT, "result_file": "/path/bam"},
            {"step": summary_steps.STEP_KESTREL, "result_file": "/path/res.tsv"},
        ],
        "stage_artifact_md5s": {
            summary_steps.STEP_FASTQ_QC: {"qc.json": "md5"},
            summary_steps.STEP_FASTQ_ALIGNMENT: {"out.bam": "md5"},
            summary_steps.STEP_KESTREL: {"output.vcf": "md5"},
        },
    }

    main_summary.write_text(json.dumps(primary_data), encoding="utf-8")
    donor_summary.write_text(json.dumps(donor_data), encoding="utf-8")

    loaded = load_prior_summary(main_summary)
    assert loaded is not None
    loaded_steps = [s.get("step") for s in loaded.get("steps", [])]
    # Differing fastp should invalidate FASTQ_QC and alignment and Kestrel
    assert summary_steps.STEP_FASTQ_QC not in loaded_steps
    assert summary_steps.STEP_FASTQ_ALIGNMENT not in loaded_steps
    assert summary_steps.STEP_KESTREL not in loaded_steps
