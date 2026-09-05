"""Tests for resume checkpoint loading and run-identity refusal checks (#20)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from vntyper.scripts.resume import load_prior_summary, resume_refusals

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
