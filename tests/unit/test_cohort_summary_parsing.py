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


def test_length_tier_counts() -> None:
    from vntyper.scripts.cohort_summary_parsing import length_tier_counts

    stats = [
        {"length_sensitivity_tier": "high"},
        {"length_sensitivity_tier": "caution"},
        {"length_sensitivity_tier": "below"},
        {"length_sensitivity_tier": "not-assessed"},
        {},
    ]
    assert length_tier_counts(stats) == {"high": 1, "caution": 1, "assessed": 3}
    assert length_tier_counts([{}, {"length_sensitivity_tier": "not-assessed"}]) is None


def test_cohort_report_shows_length_tier_kpi_only_when_assessed(tmp_path) -> None:
    import pandas as pd

    from vntyper.cli import load_config
    from vntyper.scripts import cohort_summary

    kwargs = {
        "kestrel_df": pd.DataFrame([{"Sample": "s1", "Confidence": "Negative", "Flag": "Not flagged"}]),
        "advntr_df": pd.DataFrame(),
        "summary_file": "cohort_summary.html",
        "config": load_config(None),
    }
    cohort_summary.generate_cohort_summary_report(
        output_dir=str(tmp_path / "with"), length_tier_counts={"high": 2, "caution": 5, "assessed": 9}, **kwargs
    )
    html = (tmp_path / "with" / "cohort_summary.html").read_text()
    assert 'class="kpi-card kpi-length"' in html
    assert ">2<" in html
    assert "Length Tier High" in html
    assert "5 caution · 9 with a length estimate" in html
    cohort_summary.generate_cohort_summary_report(output_dir=str(tmp_path / "without"), **kwargs)
    assert 'class="kpi-card kpi-length"' not in (tmp_path / "without" / "cohort_summary.html").read_text()


def test_unknown_cohort_tier_is_rejected() -> None:
    with pytest.raises(ValueError, match="length sensitivity tier"):
        parse_pipeline_summary({"estimated_total_repeat_count": 120.0, "length_sensitivity_tier": "severe"})
