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
    load_prior_summary,
    make_reused_step_record,
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
