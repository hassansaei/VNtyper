"""Tests for the pure alignment preflight contract."""

import logging
from dataclasses import FrozenInstanceError, fields, replace
from typing import cast
from unittest.mock import Mock

import pytest

from vntyper.scripts.alignment_binding import AlignmentBinding
from vntyper.scripts.alignment_contract import (
    FORMAT_BAM,
    FORMAT_CRAM,
    INDEX_SUFFIXES,
    AlignmentPlan,
    TimedOutProbeReason,
    index_candidate_names,
    missing_index_message,
    missing_reference_contig,
    preflight_error_payload,
    unresolvable_reference_message,
)
from vntyper.scripts.reference_binding import ReferenceBinding

pytestmark = pytest.mark.unit


def test_timed_out_probe_reason_preserves_runner_metadata_when_trimmed() -> None:
    """Candidate normalization cannot discard authenticated timeout status."""
    reason = TimedOutProbeReason("  private output\nCommand timed out after 45 seconds.\n", 45)

    stripped = reason.strip()

    assert isinstance(stripped, TimedOutProbeReason)
    assert stripped.timeout_seconds == 45
    assert stripped == "private output\nCommand timed out after 45 seconds."


def _plan(*, reference_path: str | None = "/r/genome.fa") -> AlignmentPlan:
    return AlignmentPlan(
        input_path="/input/sample.cram",
        view_path="/work/sample.bam",
        file_format=FORMAT_CRAM,
        index_path="/work/sample.bam.bai",
        reference_path=reference_path,
        reference_source="cli",
        uncovered_contigs=("chrM",),
        unmapped_scan="indexed",
    )


@pytest.mark.parametrize(
    ("file_format", "expected"),
    [
        (FORMAT_BAM, ("/data/sample.bam.bai", "/data/sample.bai", "/data/sample.bam.csi", "/data/sample.csi")),
        (
            FORMAT_CRAM,
            (
                "/data/sample.cram.crai",
                "/data/sample.crai",
                "/data/sample.cram.csi",
                "/data/sample.csi",
            ),
        ),
    ],
)
def test_index_candidates_try_appended_then_stem_for_each_supported_suffix(
    file_format: str, expected: tuple[str, ...]
) -> None:
    """A resolver can find either spelling for every index suffix."""
    assert index_candidate_names("/data/sample." + file_format, file_format) == expected


def test_bam_bai_only_excludes_csi_candidates() -> None:
    """The legacy BAI-only preflight/protection contract excludes CSI candidates."""
    assert index_candidate_names("/data/sample.bam", FORMAT_BAM, bai_only=True) == (
        "/data/sample.bam.bai",
        "/data/sample.bai",
    )


def test_cram_bai_only_still_offers_all_cram_indexes_and_never_bai() -> None:
    """The BAM-only restriction does not remove htslib-supported CRAM CSI."""
    candidates = index_candidate_names("/data/sample.cram", FORMAT_CRAM, bai_only=True)
    assert candidates == (
        "/data/sample.cram.crai",
        "/data/sample.crai",
        "/data/sample.cram.csi",
        "/data/sample.csi",
    )
    assert all(not candidate.endswith(".bai") for candidate in candidates)


def test_unknown_alignment_format_is_rejected() -> None:
    """A typo must not silently select BAM index rules."""
    with pytest.raises(ValueError, match="unknown alignment format"):
        index_candidate_names("/data/sample.sam", "sam")


def test_index_suffix_contract_has_the_supported_formats() -> None:
    """The public suffix table describes the format-specific index choices."""
    assert INDEX_SUFFIXES == {FORMAT_BAM: ("bai", "csi"), FORMAT_CRAM: ("crai", "csi")}


def test_missing_index_message_names_path_candidates_and_repair_command() -> None:
    """A fresh-index failure names the actual operation, candidates, and diagnostic command."""
    tried = ("/data/sample.bam.bai", "/data/sample.bai", "/data/sample.bam.csi", "/data/sample.csi")
    message = missing_index_message("/data/sample.bam", FORMAT_BAM, tried)
    assert "/data/sample.bam" in message
    assert all(candidate in message for candidate in tried)
    assert "samtools index" in message
    assert "Unable to build a fresh run-local BAM index" in message
    assert "Create one with" not in message


def test_unresolvable_reference_message_reports_m5_and_all_attempts() -> None:
    """Reference diagnostics preserve every source, path, and reason."""
    attempts = (
        ("cli", "/refs/cli.fa", "not found"),
        ("registry", "/refs/registry.fa", "M5 mismatch"),
    )
    message = unresolvable_reference_message("/data/sample.cram", "chr1", "abc123", attempts)
    assert "chr1" in message
    assert "abc123" in message
    assert all(part in message for attempt in attempts for part in attempt)


def test_unresolvable_reference_message_marks_missing_m5() -> None:
    """Missing sequence checksums are explicit instead of being rendered as ``None``."""
    attempts = (("cli", None, "not supplied"),)
    message = unresolvable_reference_message("/data/sample.cram", "chrM", None, attempts)
    assert "chrM" in message
    assert "no M5" in message
    assert "cli" in message
    assert "not supplied" in message


@pytest.mark.parametrize(
    ("diagnostic", "expected"),
    [
        ("[E::cram_decode_slice] Unable to fetch reference chr2:1-12070", "chr2"),
        ("Reference file given, but ref 'chr2' not present", "chr2"),
        (
            "MD5 checksum reference mismatch\n@SQ\tSN:chr2\tLN:20000\tM5:wrong",
            "chr2",
        ),
        ("Unable to fetch reference unknown:1-2", None),
        ("probe exited non-zero", None),
    ],
)
def test_missing_reference_contig_recognizes_only_known_samtools_diagnostics(
    diagnostic: str, expected: str | None
) -> None:
    """Pinned samtools failures identify the actual header contig without trusting arbitrary output."""
    assert missing_reference_contig(diagnostic, ("chr1", "chr2")) == expected


def test_alignment_plan_is_frozen_and_has_no_layout_field() -> None:
    """Layout belongs to conversion, so it cannot leak into this preflight contract."""
    plan = _plan()
    assert {field.name for field in fields(plan)} == {
        "input_path",
        "view_path",
        "file_format",
        "index_path",
        "reference_path",
        "reference_source",
        "uncovered_contigs",
        "unmapped_scan",
        "binding",
        "reference_binding",
    }
    field_name = "index_path"
    with pytest.raises(FrozenInstanceError):
        setattr(plan, field_name, "/other/index.bai")


def test_plan_close_attempts_both_bindings_and_preserves_the_first_failure(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A reference cleanup failure cannot skip alignment cleanup or be replaced by it."""
    reference_binding = Mock(spec=ReferenceBinding)
    reference_failure = RuntimeError("reference cleanup failed")
    reference_binding.close.side_effect = reference_failure
    alignment_binding = Mock(spec=AlignmentBinding)
    alignment_binding.close.side_effect = RuntimeError("alignment cleanup failed")
    plan = replace(
        _plan(),
        binding=cast(AlignmentBinding, alignment_binding),
        reference_binding=cast(ReferenceBinding, reference_binding),
    )

    with caplog.at_level(logging.ERROR), pytest.raises(RuntimeError) as raised:
        plan.close()

    assert raised.value is reference_failure
    reference_binding.close.assert_called_once_with()
    alignment_binding.close.assert_called_once_with()
    assert "Additional alignment-plan cleanup failure: alignment cleanup failed" in caplog.text


def test_the_plan_exposes_a_reference_path_not_a_preformatted_shell_fragment():
    plan = _plan(reference_path="/r/genome ref.fa")

    assert plan.reference_path == "/r/genome ref.fa"
    assert not hasattr(plan, "cram_ref_option")


def test_the_error_payload_carries_no_absolute_worker_paths_beyond_the_input():
    payload = preflight_error_payload("reference_unresolved", "msg", [("cli", None, "not supplied")])
    assert set(payload) == {"code", "message", "candidates"}
    assert payload["code"] == "reference_unresolved"
    assert payload["message"] == "msg"
    assert payload["candidates"] == [("cli", None, "not supplied")]
