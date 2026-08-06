"""Unit tests for cross_match.py.

This module had no test file. It evaluates a rule string from config against
DataFrame-derived names (AGENTS.md trap 3) and emits a summary step matched by
exact string literal downstream (trap 5). Both are silent-failure surfaces.

Two behaviours are easy to assume wrong from the name alone, so they get their
own tests rather than being taken on faith: ``determine_variant_type`` returns
``"Other"`` -- not ``"Substitution"`` -- when REF and ALT are the same length
(cross_match.py:48-49), and ``cross_match_variants`` reports ``overall_match``
as the *string* ``"Yes"``/``"No"``, never a boolean (cross_match.py:146).
"""

import csv
import logging

import pytest

from vntyper.scripts.cross_match import (
    compute_allele_change,
    cross_match_variants,
    determine_variant_type,
    extract_results_from_pipeline_summary,
    write_results_tsv,
)
from vntyper.scripts.summary_steps import STEP_ADVNTR, STEP_COVERAGE, STEP_KESTREL

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# determine_variant_type
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("ref", "alt", "expected"),
    [
        ("C", "CC", "Insertion"),
        ("C", "CGGCA", "Insertion"),
        ("CC", "C", "Deletion"),
        # Equal length returns "Other" -- NOT "Substitution". cross_match.py:48-49.
        ("C", "G", "Other"),
        ("C", "C", "Other"),
        # Non-string inputs are coerced at :42-43, so a numeric REF must not crash.
        (1, 22, "Insertion"),
    ],
)
def test_determine_variant_type(ref, alt, expected):
    assert determine_variant_type(ref, alt) == expected


# ---------------------------------------------------------------------------
# compute_allele_change
# ---------------------------------------------------------------------------


def test_compute_allele_change_strips_the_common_prefix_for_an_insertion():
    """cross_match.py:70-72 -- the return is the suffix, with no sign character."""
    assert compute_allele_change("C", "CGG", "Insertion") == "GG"


def test_compute_allele_change_strips_the_common_prefix_for_a_deletion():
    """cross_match.py:74-76."""
    assert compute_allele_change("CGG", "C", "Deletion") == "GG"


def test_compute_allele_change_returns_the_whole_allele_when_there_is_no_common_prefix():
    """cross_match.py:72 and :76 -- the fallback branch, otherwise uncovered."""
    assert compute_allele_change("C", "GGA", "Insertion") == "GGA"
    assert compute_allele_change("GGA", "C", "Deletion") == "GGA"


def test_duplication_is_treated_as_an_insertion():
    """cross_match.py:69 lists "duplication" beside "insertion"."""
    assert compute_allele_change("C", "CC", "Duplication") == "C"


def test_compute_allele_change_returns_empty_for_an_unknown_type():
    """cross_match.py:77 -- the bare `return ""`."""
    assert compute_allele_change("C", "G", "Nonsense") == ""


# ---------------------------------------------------------------------------
# cross_match_variants
# ---------------------------------------------------------------------------


def test_a_matching_pair_is_reported_as_a_match():
    """overall_match is the STRING "Yes", not a boolean. cross_match.py:146."""
    result = cross_match_variants(
        kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        advntr_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
    )
    assert result["overall_match"] == "Yes"
    assert result["matches"][0]["Match"] == "Yes"


def test_a_non_matching_pair_is_not_reported_as_a_match():
    result = cross_match_variants(
        kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        advntr_records=[{"REF": "C", "ALT": "CGGCA", "POS": 67}],
    )
    assert result["overall_match"] == "No"
    assert result["matches"][0]["Match"] == "No"


def test_an_empty_advntr_record_set_produces_no_comparisons():
    """The nested loop at :121-122 never runs, so `matches` is empty."""
    result = cross_match_variants(kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}], advntr_records=[])
    assert result["overall_match"] == "No"
    assert result["matches"] == []


def test_an_explicit_kestrel_variant_field_wins_over_the_inferred_type():
    """cross_match.py:110 -- `Variant` is used when present and non-blank."""
    result = cross_match_variants(
        kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67, "Variant": "Duplication"}],
        advntr_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
    )
    assert result["matches"][0]["Kestrel_Variant_Type"] == "Duplication"


def test_a_rule_naming_a_column_that_does_not_exist_is_silently_no_match_today(caplog):
    """CHARACTERISATION of a live trap-3 fail-open. Do not "fix" this here.

    The eval at cross_match.py:137 is wrapped in `except Exception` at :138-140,
    which sets `match = False` and logs at ERROR. A rule naming a column the
    records do not carry therefore reports "no match" rather than failing --
    the same shape as the `Poylmorhic_Call` typo that disabled a flag for
    months, and the `RU == 7` rule that never fired in its life.

    Turning it fail-closed changes what the cross-match step reports, so it
    needs a domain decision. Filed; pinned here so it cannot drift silently.
    """
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cross_match")

    result = cross_match_variants(
        kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        advntr_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        config={"cross_match": {"match_logic": "Nonexistent_Column == 1"}},
    )

    assert result["overall_match"] == "No"
    assert any("Error evaluating match logic" in r.message for r in caplog.records)


def test_the_default_match_logic_is_used_when_no_config_is_given():
    """cross_match.py:100-103 -- both arms of the config branch."""
    with_none = cross_match_variants(
        kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        advntr_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        config=None,
    )
    with_empty = cross_match_variants(
        kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        advntr_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        config={},
    )
    assert with_none["overall_match"] == with_empty["overall_match"] == "Yes"


# ---------------------------------------------------------------------------
# write_results_tsv
# ---------------------------------------------------------------------------


def test_write_results_tsv_writes_a_tab_delimited_header_and_rows(tmp_path):
    """cross_match.py:161-166 -- fieldnames come from the first record's keys."""
    output_path = tmp_path / "results.tsv"
    results = [
        {"Kestrel_POS": "67", "Kestrel_REF": "C", "Match": "Yes"},
        {"Kestrel_POS": "70", "Kestrel_REF": "G", "Match": "No"},
    ]

    write_results_tsv(results, output_path)

    with open(output_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames

    assert fieldnames == ["Kestrel_POS", "Kestrel_REF", "Match"]
    assert rows == results


def test_write_results_tsv_with_no_results_writes_no_file(tmp_path, caplog):
    """cross_match.py:158-160 -- the empty-input early return, before `open()`."""
    output_path = tmp_path / "results.tsv"
    caplog.set_level(logging.INFO, logger="vntyper.scripts.cross_match")

    write_results_tsv([], output_path)

    assert not output_path.exists()
    assert any("No results to write." in r.message for r in caplog.records)


# ---------------------------------------------------------------------------
# extract_results_from_pipeline_summary
# ---------------------------------------------------------------------------


def test_extract_results_from_pipeline_summary_finds_both_steps():
    """cross_match.py:185-189 matches step names by the shared STEP_* constants."""
    summary = {
        "steps": [
            {"step": STEP_KESTREL, "parsed_result": {"data": [{"REF": "C", "ALT": "CC"}]}},
            {"step": STEP_ADVNTR, "parsed_result": {"data": [{"REF": "C", "ALT": "CGGCA"}]}},
        ]
    }

    kestrel_records, advntr_records = extract_results_from_pipeline_summary(summary)

    assert kestrel_records == [{"REF": "C", "ALT": "CC"}]
    assert advntr_records == [{"REF": "C", "ALT": "CGGCA"}]


def test_extract_results_from_pipeline_summary_with_no_steps_key_returns_none_none():
    """The empty-input path: cross_match.py:183-184 initialise both to None and
    the loop at :185 never runs because `summary.get("steps", [])` is `[]`.
    """
    kestrel_records, advntr_records = extract_results_from_pipeline_summary({})

    assert kestrel_records is None
    assert advntr_records is None


def test_extract_results_from_pipeline_summary_leaves_the_missing_step_as_none():
    """Only Kestrel present -- adVNTR must stay `None`, not `[]`, so a caller can
    tell "step absent" apart from "step present with zero records".
    """
    summary = {"steps": [{"step": STEP_KESTREL, "parsed_result": {"data": [{"REF": "C"}]}}]}

    kestrel_records, advntr_records = extract_results_from_pipeline_summary(summary)

    assert kestrel_records == [{"REF": "C"}]
    assert advntr_records is None


def test_extract_results_from_pipeline_summary_skips_unrelated_steps():
    """A real pipeline_summary.json carries many steps (BAM Header Parsing,
    Coverage Calculation, ...) before Kestrel/adVNTR. Neither branch at
    cross_match.py:186/188 should match them, and the loop must keep going.
    """
    summary = {
        "steps": [
            {"step": STEP_COVERAGE, "parsed_result": {"data": [{"unrelated": "row"}]}},
            {"step": STEP_KESTREL, "parsed_result": {"data": [{"REF": "C"}]}},
            {"step": STEP_ADVNTR, "parsed_result": {"data": [{"REF": "G"}]}},
        ]
    }

    kestrel_records, advntr_records = extract_results_from_pipeline_summary(summary)

    assert kestrel_records == [{"REF": "C"}]
    assert advntr_records == [{"REF": "G"}]
