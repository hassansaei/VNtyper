"""Direct pipeline callers must keep their active log outside adVNTR cleanup."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from unittest.mock import patch

import pytest

from tests.support.pipeline_harness import run_pipeline_under_harness
from vntyper.modules.advntr.model_provenance import AdvntrProbeStatus, AdvntrVersionOutcome

pytestmark = pytest.mark.unit


def _unparseable_outcome() -> AdvntrVersionOutcome:
    return AdvntrVersionOutcome(
        AdvntrProbeStatus.UNPARSEABLE_SUCCESS,
        message="adVNTR version command succeeded but its response was unparseable or ambiguous.",
    )


def _verified_outcome() -> AdvntrVersionOutcome:
    return AdvntrVersionOutcome(AdvntrProbeStatus.VERIFIED, version=(2, 0, 4))


@pytest.mark.parametrize("path_case", ["exact", "output-alias"])
def test_direct_pipeline_cannot_replace_an_active_model_snapshot_log(
    tmp_path: Path,
    path_case: str,
) -> None:
    """Snapshot installation must refuse before replacing an open log with SQLite bytes."""
    output = tmp_path / "out"
    output.mkdir()
    pipeline_output: Path = output
    if path_case == "output-alias":
        pipeline_output = tmp_path / "output-alias"
        pipeline_output.symlink_to(output, target_is_directory=True)
    snapshot_log = output / "advntr" / "advntr_model.db"
    snapshot_log.parent.mkdir()
    stale_result = output / "advntr" / "cross_match_results.tsv"
    stale_result.write_text("stale result\n", encoding="utf-8")

    file_handler = logging.FileHandler(snapshot_log, encoding="utf-8")
    root_logger = logging.getLogger()
    root_logger.addHandler(file_handler)
    try:
        with (
            patch(
                "vntyper.scripts.pipeline_advntr_run_context.shutil.copyfileobj",
                wraps=shutil.copyfileobj,
            ) as model_copy,
            patch(
                "vntyper.modules.advntr.model_provenance.detect_advntr_version",
                return_value=_verified_outcome(),
            ) as detector,
        ):
            harness = run_pipeline_under_harness(
                output,
                pipeline_output_dir=pipeline_output,
                extra_modules=["advntr"],
                log_file=str(snapshot_log),
                expect_failure=True,
            )
    finally:
        root_logger.removeHandler(file_handler)
        file_handler.close()

    assert not snapshot_log.read_bytes().startswith(b"SQLite format 3"), (
        "snapshot installation replaced the live log and discarded its diagnostics after close"
    )
    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert "aliases an adVNTR destructive preflight destination" in snapshot_log.read_text(encoding="utf-8")
    assert stale_result.read_text(encoding="utf-8") == "stale result\n"
    model_copy.assert_not_called()
    detector.assert_not_called()
    harness.stages["run_kestrel"].assert_not_called()


@pytest.mark.parametrize("log_spelling", ["exact", "redundant"])
def test_direct_pipeline_preserves_colliding_active_log_and_refuses_before_cleanup_or_probe(
    tmp_path: Path,
    log_spelling: str,
) -> None:
    """Closing the active handler must not make the refusal diagnostics disappear."""
    output = tmp_path / "out"
    output.mkdir()
    stale_result = output / "advntr" / "cross_match_results.tsv"
    stale_result.parent.mkdir()
    stale_result.write_text("stale result\n", encoding="utf-8")
    if log_spelling == "redundant":
        (output / "existing").mkdir()
        log_file = output / "existing" / ".." / "pipeline_summary.json"
    else:
        log_file = output / "pipeline_summary.json"

    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    root_logger = logging.getLogger()
    root_logger.addHandler(file_handler)
    try:
        with patch(
            "vntyper.modules.advntr.model_provenance.detect_advntr_version",
            return_value=_unparseable_outcome(),
        ) as detector:
            harness = run_pipeline_under_harness(
                output,
                extra_modules=["advntr"],
                log_file=str(log_file),
                expect_failure=True,
            )
    finally:
        root_logger.removeHandler(file_handler)
        file_handler.close()

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert log_file.is_file(), "cleanup unlinked the active log, so its diagnostics vanished when the handler closed"
    assert "aliases an adVNTR destructive preflight destination" in log_file.read_text(encoding="utf-8")
    assert stale_result.read_text(encoding="utf-8") == "stale result\n"
    detector.assert_not_called()
    harness.stages["run_kestrel"].assert_not_called()


def test_direct_pipeline_preserves_a_selected_archive_log_before_cleanup_or_probe(tmp_path: Path) -> None:
    """The selected sibling archive is as cleanup-owned as an in-tree summary."""
    output = tmp_path / "out"
    output.mkdir()
    archive_log = Path(f"{output}.zip")
    file_handler = logging.FileHandler(archive_log, encoding="utf-8")
    root_logger = logging.getLogger()
    root_logger.addHandler(file_handler)
    try:
        with patch(
            "vntyper.modules.advntr.model_provenance.detect_advntr_version",
            return_value=_unparseable_outcome(),
        ) as detector:
            harness = run_pipeline_under_harness(
                output,
                extra_modules=["advntr"],
                archive_results=True,
                archive_format="zip",
                log_file=str(archive_log),
                expect_failure=True,
            )
    finally:
        root_logger.removeHandler(file_handler)
        file_handler.close()

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert archive_log.is_file()
    assert "aliases an adVNTR destructive preflight destination" in archive_log.read_text(encoding="utf-8")
    detector.assert_not_called()
    harness.stages["run_kestrel"].assert_not_called()


def test_direct_advntr_pipeline_accepts_an_ordinary_pipeline_log(tmp_path: Path) -> None:
    """A direct caller's ordinary application log remains outside cleanup ownership."""
    output = tmp_path / "out"
    pipeline_log = output / "pipeline.log"
    output.mkdir()
    pipeline_log.write_text("active pipeline log\n", encoding="utf-8")

    harness = run_pipeline_under_harness(output, extra_modules=["advntr"], log_file=str(pipeline_log))

    assert harness.error is None
    assert pipeline_log.read_text(encoding="utf-8") == "active pipeline log\n"
