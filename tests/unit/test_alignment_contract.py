"""Tests for the pure alignment preflight contract."""

from dataclasses import FrozenInstanceError, fields

import pytest

from vntyper.scripts.alignment_contract import (
    FORMAT_BAM,
    FORMAT_CRAM,
    INDEX_SUFFIXES,
    AlignmentPlan,
    index_candidate_names,
    missing_index_message,
    preflight_error_payload,
    unresolvable_reference_message,
)

pytestmark = pytest.mark.unit


def _plan(*, reference_path: str | None = "/r/genome.fa") -> AlignmentPlan:
    return AlignmentPlan(
        input_path="/input/sample.cram",
        view_path="/work/sample.bam",
        file_format=FORMAT_CRAM,
        index_path="/work/sample.bam.bai",
        reference_path=reference_path,
        reference_source="cli",
        uncovered_contigs=("chrM",),
        unmapped_scan="offset",
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
    """The offset reader accepts BAI only, while general BAM preflight accepts CSI."""
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
    }
    field_name = "index_path"
    with pytest.raises(FrozenInstanceError):
        setattr(plan, field_name, "/other/index.bai")


def test_the_reference_fragment_is_quoted_because_builders_interpolate_it_raw():
    plan = _plan(reference_path="/r/genome ref.fa")
    assert plan.cram_ref_option == "-T '/r/genome ref.fa'"


def test_no_reference_yields_an_empty_fragment_so_no_ref_crams_are_byte_identical():
    assert _plan(reference_path=None).cram_ref_option == ""


def test_the_error_payload_carries_no_absolute_worker_paths_beyond_the_input():
    payload = preflight_error_payload("reference_unresolved", "msg", [("cli", None, "not supplied")])
    assert set(payload) == {"code", "message", "candidates"}
    assert payload["code"] == "reference_unresolved"
    assert payload["message"] == "msg"
    assert payload["candidates"] == [("cli", None, "not supplied")]
