"""A failed run ends with its cause and where to get help, not with a stack (#338)."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from tests.support.pipeline_harness import run_pipeline_under_harness
from vntyper.scripts.failure_help import (
    KESTREL_ERROR_HINTS,
    TROUBLESHOOTING_URL,
    kestrel_error_hint,
    summarize_failure,
)

pytestmark = pytest.mark.unit


def test_the_help_link_is_the_published_troubleshooting_page():
    assert TROUBLESHOOTING_URL == "https://hassansaei.github.io/VNtyper/user-guide/troubleshooting/"


def test_the_troubleshooting_page_the_link_points_at_exists():
    root = Path(__file__).resolve().parents[2]
    assert (root / "docs" / "user-guide" / "troubleshooting.md").is_file()
    assert "user-guide/troubleshooting.md" in (root / "mkdocs.yml").read_text(encoding="utf-8")


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        ("... ERROR ... - Error reading reference sequence(s): IO error ...", "install-references -d reference"),
        ("... ERROR ... - Cannot open indexed k-mer count (IKC) file: File does not exist", "free disk space"),
        ("... ERROR ... - null: Error setting k-mer counts (in map post-run)", "free disk space"),
        ("... ERROR ... - File not found while variant writer: /out/kestrel/output.vcf", "output directory"),
    ],
)
def test_each_known_kestrel_error_gets_its_fix(line, expected):
    hint = kestrel_error_hint([line])
    assert hint is not None and expected in hint


def test_an_unknown_error_gets_no_hint():
    assert kestrel_error_hint(["... ERROR ... - something nobody has seen"]) is None
    assert kestrel_error_hint([]) is None


def test_the_earliest_known_line_decides():
    """Kestrel repeats one cause through layers; the first matching line is the cause."""
    lines = ["... - unrelated", "... - Error setting k-mer counts", "... - Error reading reference sequence"]
    assert kestrel_error_hint(lines) == kestrel_error_hint([lines[1]])


def test_every_hint_is_actionable_text():
    for marker, hint in KESTREL_ERROR_HINTS:
        assert marker and hint.endswith(".")


def test_the_summary_uses_the_first_line_and_links_help():
    exc = RuntimeError("Kestrel cannot run because files are missing:\n  - a\n  - b")
    summary = summarize_failure(exc)
    assert summary.startswith("VNtyper failed: Kestrel cannot run because files are missing. ")
    assert "- a" not in summary
    assert summary.endswith(f"Help: {TROUBLESHOOTING_URL}")


def test_an_empty_message_falls_back_to_the_exception_type():
    assert summarize_failure(ValueError()).startswith("VNtyper failed: ValueError. ")


def test_a_failed_pipeline_ends_with_the_summary_line(tmp_path, caplog):
    """The operator's last line is the cause and the help link, after the traceback."""

    def kestrel_fails(*_args, **_kwargs):
        raise RuntimeError("Kestrel reported 1 error(s) for k-mer size 20 but exited 0\n  detail")

    with caplog.at_level(logging.INFO):
        harness = run_pipeline_under_harness(
            tmp_path / "out", expect_failure=True, stage_side_effects={"run_kestrel": kestrel_fails}
        )

    assert isinstance(harness.error, SystemExit) and harness.error.code == 1
    critical = [record for record in caplog.records if record.levelno == logging.CRITICAL]
    assert len(critical) == 1
    assert (
        critical[0]
        .getMessage()
        .startswith("VNtyper failed: Kestrel reported 1 error(s) for k-mer size 20 but exited 0. ")
    )
    assert TROUBLESHOOTING_URL in critical[0].getMessage()
    traceback_record = next(record for record in caplog.records if record.getMessage() == "An error occurred")
    assert traceback_record.exc_info is not None, "the traceback must still be logged for bug reports"
    assert caplog.records.index(traceback_record) < caplog.records.index(critical[0])
