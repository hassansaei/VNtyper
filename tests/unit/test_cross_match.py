"""Unit tests for cross_match.py.

This module pins the structured comparator contract used for cross-match and emits
a summary step matched by exact string literal downstream (AGENTS.md trap 5).

Two behaviours are easy to assume wrong from the name alone, so they get their
own tests rather than being taken on faith: ``determine_variant_type`` returns
``"Other"`` -- not ``"Substitution"`` -- when REF and ALT are the same length
(cross_match.py:48-49), and ``cross_match_variants`` reports ``overall_match``
as the *string* ``"Yes"``/``"No"``, never a boolean.
"""

import copy
import csv
import json
import logging
from pathlib import Path
from unittest.mock import Mock

import pytest

import vntyper
import vntyper.scripts.cross_match as cross_match_module
from vntyper.scripts.cross_match import (
    CROSS_MATCH_COLUMNS,
    DEFAULT_MATCH_RULE,
    compute_allele_change,
    cross_match_variants,
    determine_variant_type,
    extract_results_from_pipeline_summary,
    write_results_tsv,
)
from vntyper.scripts.summary_steps import STEP_ADVNTR, STEP_COVERAGE, STEP_KESTREL

pytestmark = pytest.mark.unit

LEGACY_MATCH_LOGIC = (
    "Kestrel_Allele_Change == Advntr_Allele_Change and Kestrel_Variant_Type.lower() == Advntr_Variant_Type.lower()"
)


def _shipped_cross_match_config():
    """The ``cross_match`` block of the shipped ``vntyper/config.json``, read from disk.

    Read from the file rather than from a literal so that adding a second rule to the
    shipped config fails the enumeration below instead of silently going untested.
    """
    path = Path(vntyper.__file__).parent / "config.json"
    assert path.exists(), f"config.json not found at {path}"
    return json.loads(path.read_text(encoding="utf-8"))["cross_match"]


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
    """overall_match is the STRING "Yes", not a boolean. cross_match.py:177."""
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
    """The nested loop at :129-130 never runs, so `matches` is empty."""
    result = cross_match_variants(kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}], advntr_records=[])
    assert result["overall_match"] == "No"
    assert result["matches"] == []


def test_an_explicit_kestrel_variant_field_wins_over_the_inferred_type():
    """cross_match.py:118 -- `Variant` is used when present and non-blank."""
    result = cross_match_variants(
        kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67, "Variant": "Duplication"}],
        advntr_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
    )
    assert result["matches"][0]["Kestrel_Variant_Type"] == "Duplication"


@pytest.mark.parametrize(
    ("kestrel_records", "advntr_records"),
    [
        pytest.param([], [], id="both-empty"),
        pytest.param([{"REF": "C", "ALT": "CC", "POS": 67}], [], id="adVNTR-empty"),
        pytest.param(
            [{"REF": "C", "ALT": "CC", "POS": 67}],
            [{"REF": "C", "ALT": "CC", "POS": 67}],
            id="both-nonempty",
        ),
    ],
)
@pytest.mark.parametrize(
    "match_rule",
    [
        pytest.param({"all": "not-a-list"}, id="malformed-nested-rule"),
        pytest.param(
            {
                "all": [
                    {
                        "left": {"column": "Nonexistent_Column"},
                        "operator": "eq",
                        "right": {"literal": 1},
                    }
                ]
            },
            id="unknown-column",
        ),
    ],
)
def test_invalid_rules_raise_before_record_processing_or_mutation(
    kestrel_records, advntr_records, match_rule, monkeypatch
):
    before_kestrel = copy.deepcopy(kestrel_records)
    before_advntr = copy.deepcopy(advntr_records)
    determine = Mock(side_effect=AssertionError("record preprocessing must not run"))
    monkeypatch.setattr(cross_match_module, "determine_variant_type", determine)

    with pytest.raises(ValueError, match=r"cross_match\.match_rule"):
        cross_match_variants(
            kestrel_records=kestrel_records,
            advntr_records=advntr_records,
            config={"cross_match": {"match_rule": match_rule}},
        )

    assert kestrel_records == before_kestrel
    assert advntr_records == before_advntr
    determine.assert_not_called()


@pytest.mark.parametrize(
    "legacy_rule",
    [
        f" {LEGACY_MATCH_LOGIC}",
        f"{LEGACY_MATCH_LOGIC} ",
        LEGACY_MATCH_LOGIC.replace("Advntr_Allele_Change", "Advntr_ALT"),
        f"{LEGACY_MATCH_LOGIC} and Kestrel_POS == Advntr_POS",
        "Kestrel_Variant_Type.upper() == Advntr_Variant_Type.upper()",
        "__import__('os').system('id')",
        "Kestrel_Variant_Type.__class__ == str",
        "Kestrel_Allele_Change[0] == Advntr_Allele_Change[0]",
        "all(x for x in [True])",
        "(lambda: True)()",
    ],
)
def test_only_the_exact_historical_match_logic_is_accepted(legacy_rule):
    with pytest.raises(ValueError, match=r"cross_match\.match_logic uses an unsupported legacy expression"):
        cross_match_variants([], [], config={"cross_match": {"match_logic": legacy_rule}})


def test_the_legacy_match_logic_key_rejects_structured_rules():
    with pytest.raises(ValueError, match=r"cross_match\.match_logic must contain the exact historical string"):
        cross_match_variants([], [], config={"cross_match": {"match_logic": DEFAULT_MATCH_RULE}})


def test_both_cross_match_rule_keys_are_rejected_loudly():
    with pytest.raises(ValueError) as excinfo:
        cross_match_variants(
            [],
            [],
            config={"cross_match": {"match_rule": DEFAULT_MATCH_RULE, "match_logic": LEGACY_MATCH_LOGIC}},
        )

    message = str(excinfo.value)
    assert "match_rule" in message
    assert "match_logic" in message


@pytest.mark.parametrize(
    ("config", "expected"),
    [
        pytest.param([], "cross_match configuration must be a mapping", id="root-not-mapping"),
        pytest.param(
            {"cross_match": []},
            "cross_match configuration block must be a mapping",
            id="block-not-mapping",
        ),
    ],
)
def test_malformed_cross_match_configuration_containers_fail_loudly(config, expected, caplog):
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cross_match")

    with pytest.raises(ValueError, match=expected) as excinfo:
        cross_match_variants([], [], config=config)

    assert str(excinfo.value) in [record.message for record in caplog.records]


@pytest.mark.parametrize(
    ("cross_match_config", "unsupported_keys"),
    [
        pytest.param({"match_rul": DEFAULT_MATCH_RULE}, ["match_rul"], id="misspelled-rule-only"),
        pytest.param(
            {"match_rule": DEFAULT_MATCH_RULE, "match_logc": "typo"},
            ["match_logc"],
            id="extra-key-with-structured-rule",
        ),
        pytest.param(
            {"match_logic": LEGACY_MATCH_LOGIC, "zzz": "secret", "aaa": "secret"},
            ["aaa", "zzz"],
            id="sorted-extra-keys-with-legacy-rule",
        ),
    ],
)
@pytest.mark.parametrize(
    ("kestrel_records", "advntr_records"),
    [
        pytest.param([], [], id="empty-records"),
        pytest.param(
            [{"REF": "C", "ALT": "CC", "POS": 67}],
            [{"REF": "C", "ALT": "CC", "POS": 67}],
            id="nonempty-records",
        ),
    ],
)
def test_unsupported_cross_match_block_keys_raise_before_record_mutation(
    cross_match_config, unsupported_keys, kestrel_records, advntr_records, monkeypatch, caplog
):
    """A non-empty wrapper is valid only when every key belongs to its one-rule schema."""
    before_kestrel = copy.deepcopy(kestrel_records)
    before_advntr = copy.deepcopy(advntr_records)
    determine = Mock(side_effect=AssertionError("record preprocessing must not run"))
    monkeypatch.setattr(cross_match_module, "determine_variant_type", determine)
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cross_match")

    expected = "cross_match configuration block contains unsupported keys: " + ", ".join(
        repr(key) for key in unsupported_keys
    )
    with pytest.raises(ValueError, match="unsupported keys") as excinfo:
        cross_match_variants(
            kestrel_records=kestrel_records,
            advntr_records=advntr_records,
            config={"cross_match": cross_match_config},
        )

    assert str(excinfo.value) == expected
    assert caplog.messages[-1] == expected
    assert "secret" not in str(excinfo.value)
    assert kestrel_records == before_kestrel
    assert advntr_records == before_advntr
    determine.assert_not_called()


def test_null_comparison_is_false():
    rule = {
        "all": [
            {
                "left": {"column": "Kestrel_POS"},
                "operator": "eq",
                "right": {"column": "Advntr_POS"},
            }
        ]
    }
    result = cross_match_variants(
        [{"REF": "C", "ALT": "CC", "POS": None}],
        [{"REF": "C", "ALT": "CC", "POS": None}],
        config={"cross_match": {"match_rule": rule}},
    )

    assert result["matches"][0]["Match"] == "No"


def test_incompatible_record_value_types_fail_loudly_without_record_disclosure(caplog):
    caplog.set_level(logging.ERROR)
    rule = {
        "all": [
            {
                "left": {"column": "Kestrel_POS"},
                "operator": "eq",
                "right": {"column": "Advntr_POS"},
            }
        ]
    }

    with pytest.raises(ValueError, match=r"cross_match\.match_rule\.all\[0\] requires compatible families") as excinfo:
        cross_match_variants(
            [{"REF": "C", "ALT": "CC", "POS": 67}],
            [{"REF": "C", "ALT": "CC", "POS": "67"}],
            config={"cross_match": {"match_rule": rule}},
        )

    assert "Record:" not in str(excinfo.value)
    assert str(excinfo.value) in [record.message for record in caplog.records]


def test_the_default_match_rule_is_used_when_no_config_is_given():
    """cross_match.py:108-111 -- both arms of the config branch."""
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
# cross_match_variants: the correctly-configured path, enumerated
#
# Making the eval fail closed is only safe if it changes nothing for a rule that
# evaluates. These two tables are the proof: every configuration the shipped tree
# can produce, crossed with every shape of input, against expected values captured
# from the code BEFORE the fail-closed change.
# ---------------------------------------------------------------------------


def _row(k_pos, k_ref, k_alt, k_change, k_type, a_pos, a_ref, a_alt, a_change, a_type, match):
    """One `matches` entry, field for field in the order `cross_match_variants` builds it.

    Key order is load-bearing, not cosmetic: `write_results_tsv` takes the TSV
    header from `results[0].keys()`, so a reordering changes the output file.
    """
    return {
        "Kestrel_POS": k_pos,
        "Kestrel_REF": k_ref,
        "Kestrel_ALT": k_alt,
        "Kestrel_Allele_Change": k_change,
        "Kestrel_Variant_Type": k_type,
        "Advntr_POS": a_pos,
        "Advntr_REF": a_ref,
        "Advntr_ALT": a_alt,
        "Advntr_Allele_Change": a_change,
        "Advntr_Variant_Type": a_type,
        "Match": match,
    }


# Every supported way the shipped tree can arrive at a rule. The exact historical
# string is a one-release compatibility adapter, not a general expression surface.
_SHIPPED_CONFIGURATIONS = [
    pytest.param(None, id="config-is-None"),
    pytest.param({}, id="config-without-a-cross_match-key"),
    pytest.param({"cross_match": {}}, id="cross_match-without-a-match_rule-key"),
    pytest.param({"cross_match": {"match_rule": DEFAULT_MATCH_RULE}}, id="explicit-DEFAULT_MATCH_RULE"),
    pytest.param({"cross_match": {"match_logic": LEGACY_MATCH_LOGIC}}, id="exact-historical-match_logic"),
    pytest.param({"cross_match": _shipped_cross_match_config()}, id="shipped-vntyper-config.json"),
]

_MATCH_SCENARIOS = [
    pytest.param(
        [{"REF": "C", "ALT": "CC", "POS": "67"}],
        [{"REF": "C", "ALT": "CC", "POS": "67"}],
        {
            "matches": [_row("67", "C", "CC", "C", "Insertion", "67", "C", "CC", "C", "Insertion", "Yes")],
            "overall_match": "Yes",
        },
        id="identical-duplication-call",
    ),
    pytest.param(
        [{"REF": "C", "ALT": "CC", "POS": "67"}],
        [{"REF": "C", "ALT": "CGGCA", "POS": "67"}],
        {
            "matches": [_row("67", "C", "CC", "C", "Insertion", "67", "C", "CGGCA", "GGCA", "Insertion", "No")],
            "overall_match": "No",
        },
        id="different-inserted-sequence",
    ),
    pytest.param(
        [{"REF": "CC", "ALT": "C", "POS": "67"}],
        [{"REF": "C", "ALT": "CC", "POS": "67"}],
        {
            "matches": [_row("67", "CC", "C", "C", "Deletion", "67", "C", "CC", "C", "Insertion", "No")],
            "overall_match": "No",
        },
        id="same-allele-change-opposite-variant-type",
    ),
    pytest.param(
        [{"REF": "C", "ALT": "CC", "POS": "67", "Variant": "Insertion"}],
        [{"REF": "C", "ALT": "CC", "POS": "67"}],
        {
            "matches": [_row("67", "C", "CC", "C", "Insertion", "67", "C", "CC", "C", "Insertion", "Yes")],
            "overall_match": "Yes",
        },
        id="explicit-Kestrel-Variant-field-agrees",
    ),
    pytest.param(
        [{"REF": "C", "ALT": "CC", "POS": "67", "Variant": "insertion"}],
        [{"REF": "C", "ALT": "CC", "POS": "67"}],
        {
            "matches": [_row("67", "C", "CC", "C", "insertion", "67", "C", "CC", "C", "Insertion", "Yes")],
            "overall_match": "Yes",
        },
        id="variant-type-comparison-is-case-insensitive",
    ),
    pytest.param(
        [{"REF": "C", "ALT": "CC", "POS": "67", "Variant": "Duplication"}],
        [{"REF": "C", "ALT": "CC", "POS": "67"}],
        {
            "matches": [_row("67", "C", "CC", "C", "Duplication", "67", "C", "CC", "C", "Insertion", "No")],
            "overall_match": "No",
        },
        id="explicit-Kestrel-Variant-field-disagrees",
    ),
    pytest.param(
        [{"REF": "C", "ALT": "CC", "POS": "67"}, {"REF": "G", "ALT": "GAT", "POS": "70"}],
        [{"REF": "C", "ALT": "CGGCA", "POS": "67"}, {"REF": "G", "ALT": "GAT", "POS": "70"}],
        {
            "matches": [
                _row("67", "C", "CC", "C", "Insertion", "67", "C", "CGGCA", "GGCA", "Insertion", "No"),
                _row("67", "C", "CC", "C", "Insertion", "70", "G", "GAT", "AT", "Insertion", "No"),
                _row("70", "G", "GAT", "AT", "Insertion", "67", "C", "CGGCA", "GGCA", "Insertion", "No"),
                _row("70", "G", "GAT", "AT", "Insertion", "70", "G", "GAT", "AT", "Insertion", "Yes"),
            ],
            "overall_match": "Yes",
        },
        id="two-by-two-one-matching-pair",
    ),
    pytest.param(
        [{"REF": "C", "ALT": "CC", "POS": "67"}],
        [],
        {"matches": [], "overall_match": "No"},
        id="no-adVNTR-records",
    ),
    pytest.param(
        [],
        [{"REF": "C", "ALT": "CC", "POS": "67"}],
        {"matches": [], "overall_match": "No"},
        id="no-Kestrel-records",
    ),
]


def test_the_shipped_cross_match_configuration_is_the_default_rule():
    """The enumeration above is complete only while this holds.

    `vntyper/config.json` must carry exactly one `cross_match` key, and its
    `match_rule` must equal DEFAULT_MATCH_RULE. Adding a second rule, or diverging
    the shipped one from the default, fails rather than going unexercised.
    """
    shipped = _shipped_cross_match_config()

    assert sorted(shipped) == ["match_rule"]
    assert shipped["match_rule"] == DEFAULT_MATCH_RULE
    assert (
        frozenset(
            {
                "Kestrel_POS",
                "Kestrel_REF",
                "Kestrel_ALT",
                "Kestrel_Allele_Change",
                "Kestrel_Variant_Type",
                "Advntr_POS",
                "Advntr_REF",
                "Advntr_ALT",
                "Advntr_Allele_Change",
                "Advntr_Variant_Type",
            }
        )
        == CROSS_MATCH_COLUMNS
    )


@pytest.mark.parametrize("config", _SHIPPED_CONFIGURATIONS)
@pytest.mark.parametrize(("kestrel_records", "advntr_records", "expected"), _MATCH_SCENARIOS)
def test_a_supported_rule_preserves_the_captured_cross_match_output(config, kestrel_records, advntr_records, expected):
    """The expected values were captured before the comparator migration.

    `cross_match_variants` mutates the records it is given (it writes back
    `Variant_Type` and `Allele_Change`), so each case gets its own deep copy.
    """
    result = cross_match_variants(
        kestrel_records=copy.deepcopy(kestrel_records),
        advntr_records=copy.deepcopy(advntr_records),
        config=copy.deepcopy(config),
    )

    assert result == expected
    # Not just equal contents -- the same column order, which becomes the TSV header.
    assert [list(row) for row in result["matches"]] == [list(row) for row in expected["matches"]]
    # overall_match is the STRING "Yes"/"No", never a boolean; downstream matches on it.
    assert isinstance(result["overall_match"], str)


# ---------------------------------------------------------------------------
# write_results_tsv
# ---------------------------------------------------------------------------


def test_write_results_tsv_writes_a_tab_delimited_header_and_rows(tmp_path):
    """cross_match.py:192-197 -- fieldnames come from the first record's keys."""
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
    """cross_match.py:189-191 -- the empty-input early return, before `open()`."""
    output_path = tmp_path / "results.tsv"
    caplog.set_level(logging.INFO, logger="vntyper.scripts.cross_match")

    write_results_tsv([], output_path)

    assert not output_path.exists()
    assert any("No results to write." in r.message for r in caplog.records)


# ---------------------------------------------------------------------------
# extract_results_from_pipeline_summary
# ---------------------------------------------------------------------------


def test_extract_results_from_pipeline_summary_finds_both_steps():
    """cross_match.py:216-220 matches step names by the shared STEP_* constants."""
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
    """The empty-input path: cross_match.py:214-215 initialise both to None and
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
    cross_match.py:217/219 should match them, and the loop must keep going.
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
