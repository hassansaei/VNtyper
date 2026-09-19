"""Unit tests for cohort summary parsing."""

from __future__ import annotations

import pytest

from vntyper.scripts.cohort_summary_parsing import parse_pipeline_summary
from vntyper.scripts.length_standard_model import LENGTH_WARNING_SENSITIVITY_CUTOFF

pytestmark = pytest.mark.unit


def test_parse_pipeline_summary_empty() -> None:
    kestrel, advntr, stats = parse_pipeline_summary({})
    assert kestrel == []
    assert advntr == []
    assert stats == {
        "runtime": "N/A",
        "version": "N/A",
        "assembly": "N/A",
        "pipeline": "N/A",
        "coverage": {},
    }


def test_parse_pipeline_summary_extracts_length_stats_and_warnings() -> None:
    summary = {
        "estimated_total_repeat_count": 115.5,
        "length_estimation_warnings": [LENGTH_WARNING_SENSITIVITY_CUTOFF],
    }
    _, _, stats = parse_pipeline_summary(summary)
    assert stats["estimated_total_repeat_count"] == 115.5
    assert stats["length_warning"] == LENGTH_WARNING_SENSITIVITY_CUTOFF


def test_parse_pipeline_summary_extracts_length_stats_with_no_warnings() -> None:
    summary = {
        "estimated_total_repeat_count": 95.0,
        "length_estimation_warnings": [],
    }
    _, _, stats = parse_pipeline_summary(summary)
    assert stats["estimated_total_repeat_count"] == 95.0
    assert stats["length_warning"] == "none"


def test_parse_pipeline_summary_extracts_multiple_warnings() -> None:
    summary = {
        "estimated_total_repeat_count": 125.0,
        "length_estimation_warnings": ["feature_A_outside_training_range", LENGTH_WARNING_SENSITIVITY_CUTOFF],
    }
    _, _, stats = parse_pipeline_summary(summary)
    assert stats["estimated_total_repeat_count"] == 125.0
    assert stats["length_warning"] == f"feature_A_outside_training_range;{LENGTH_WARNING_SENSITIVITY_CUTOFF}"


def test_parse_pipeline_summary_extracts_non_list_warning() -> None:
    summary = {
        "estimated_total_repeat_count": 115.0,
        "length_estimation_warnings": "unexpected_string_warning",  # type: ignore[dict-item]
    }
    _, _, stats = parse_pipeline_summary(summary)
    assert stats["estimated_total_repeat_count"] == 115.0
    assert stats["length_warning"] == "unexpected_string_warning"


def test_parse_pipeline_summary_runtime_computation() -> None:
    summary = {
        "pipeline_start": "2026-01-01T00:00:00",
        "pipeline_end": "2026-01-01T00:01:30",
    }
    _, _, stats = parse_pipeline_summary(summary)
    assert stats["runtime"] == "90.00 seconds"


def test_parse_pipeline_summary_invalid_runtime_raises() -> None:
    summary = {
        "pipeline_start": "not-a-timestamp",
        "pipeline_end": "2026-01-01T00:01:30",
    }
    with pytest.raises(ValueError):
        parse_pipeline_summary(summary)


def test_cohort_stats_table_and_frame_integration() -> None:
    from vntyper.scripts.cohort_tables import additional_stats_frame, stats_table_html

    summary = {
        "estimated_total_repeat_count": 115.5,
        "length_estimation_warnings": [LENGTH_WARNING_SENSITIVITY_CUTOFF],
    }
    _, _, stats = parse_pipeline_summary(summary)
    stats["Sample"] = "sample_test"
    frame = additional_stats_frame([stats])
    assert "estimated_total_repeat_count" in frame.columns
    assert "length_warning" in frame.columns
    assert frame["estimated_total_repeat_count"].iloc[0] == 115.5
    assert frame["length_warning"].iloc[0] == LENGTH_WARNING_SENSITIVITY_CUTOFF

    html = stats_table_html(frame)
    assert "115.5" in html
    assert LENGTH_WARNING_SENSITIVITY_CUTOFF in html
