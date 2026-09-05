"""Tests for pipeline resume orchestration and output directory warnings (#20)."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from tests.support.pipeline_harness import run_pipeline_under_harness

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
