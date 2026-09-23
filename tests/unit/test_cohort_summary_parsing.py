"""Unit tests for cohort summary parsing."""

from __future__ import annotations

import pytest

from vntyper.scripts.cohort_summary_parsing import parse_pipeline_summary

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


@pytest.mark.parametrize(
    ("summary", "expected"),
    [
        ({}, {}),
        (
            {
                "estimated_total_repeat_count": None,
                "length_estimation_status": "unavailable",
                "length_sensitivity_tier": "not-assessed",
            },
            {"estimated_total_repeat_count": None, "length_sensitivity_tier": "not-assessed"},
        ),
        (
            {"estimated_total_repeat_count": 123.456789, "length_sensitivity_tier": "caution"},
            {"estimated_total_repeat_count": 123.5, "length_sensitivity_tier": "caution"},
        ),
        ({"estimated_total_repeat_count": 95.0}, {"estimated_total_repeat_count": 95.0}),
        ({"estimated_total_repeat_count": True}, {"estimated_total_repeat_count": None}),
    ],
)
def test_length_columns(summary: dict[str, object], expected: dict[str, object]) -> None:
    stats = parse_pipeline_summary(summary)[2]
    got = {key: stats[key] for key in ("estimated_total_repeat_count", "length_sensitivity_tier") if key in stats}
    assert got == expected
    assert "length_warning" not in stats


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


def test_mixed_cohort_renders_blanks_not_placeholder_strings() -> None:
    from vntyper.scripts.cohort_tables import additional_stats_frame, stats_table_html

    with_length = parse_pipeline_summary({"estimated_total_repeat_count": 160.25, "length_sensitivity_tier": "high"})[2]
    without_length = parse_pipeline_summary({})[2]
    with_length["Sample"], without_length["Sample"] = "long", "short"
    frame = additional_stats_frame([with_length, without_length])
    assert frame["estimated_total_repeat_count"].iloc[0] == 160.2
    assert frame["length_sensitivity_tier"].iloc[0] == "high"

    html = stats_table_html(frame)
    assert "160.2" in html
    assert ">high<" in html
    for placeholder in (">None<", ">nan<", ">NaN<", ">none<"):
        assert placeholder not in html
