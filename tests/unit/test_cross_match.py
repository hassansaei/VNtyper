"""Unit tests for cross_match.py.

This module had no test file. It evaluates a rule string from config against
DataFrame-derived names (AGENTS.md trap 3) and emits a summary step matched by
exact string literal downstream (trap 5). The eval no longer fails open -- a rule
that cannot be evaluated now raises rather than reporting "no match", see
``test_a_rule_naming_a_column_that_does_not_exist_raises`` below. The step-name
literal is still a silent-failure surface.

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

import pytest

import vntyper
from vntyper.scripts.cross_match import (
    DEFAULT_MATCH_LOGIC,
    compute_allele_change,
    cross_match_variants,
    determine_variant_type,
    extract_results_from_pipeline_summary,
    write_results_tsv,
)
from vntyper.scripts.summary_steps import STEP_ADVNTR, STEP_COVERAGE, STEP_KESTREL

pytestmark = pytest.mark.unit


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


def test_a_rule_naming_a_column_that_does_not_exist_raises(caplog):
    """SPECIFICATION. Was `..._today` characterisation of a live trap-3 fail-open.

    Decided by the repository owner on branch `fix/issue-181-197-followups`: the
    cross-match step must fail loudly on a rule it cannot evaluate rather than
    report "no match". The issue draft is
    `.superpowers/sdd/2026-08-06-issue-181-197-followups-plan/issue-cross-match-fail-open.md`
    (never filed on GitHub, so there is no issue number to cite -- see
    `task-8-report.md`).

    A `NameError` is a *configuration* defect: the ten names the rule can use are
    the fixed dict literal in `cross_match_variants`, identical for every record
    pair, so a rule naming anything else fails for all of them. Reporting "no
    match" made that indistinguishable from genuine discordance -- the same shape
    as the `Poylmorhic_Call` typo that disabled a flag for months, and the
    `RU == 7` rule that never fired in its life.

    The message must name the rule and the available columns, because the whole
    point is that the next person sees which column name is wrong.
    """
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cross_match")

    with pytest.raises(ValueError) as excinfo:
        cross_match_variants(
            kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
            advntr_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
            config={"cross_match": {"match_logic": "Nonexistent_Column == 1"}},
        )

    message = str(excinfo.value)
    # The offending rule string, verbatim.
    assert "Nonexistent_Column == 1" in message
    assert "NameError" in message
    # The columns that WERE available, so the typo is visible without reading the source.
    for column in ("Kestrel_Allele_Change", "Kestrel_Variant_Type", "Advntr_Allele_Change", "Advntr_Variant_Type"):
        assert column in message
    # Classified as a configuration defect, not as a bad record pair: no record dump.
    assert message.startswith("Invalid cross_match match_logic")
    assert "Record:" not in message
    # Repo convention: logger.error(msg) then raise ValueError(msg), same message.
    assert message in [r.message for r in caplog.records]


def test_a_rule_that_is_not_a_valid_expression_raises(caplog):
    """SPECIFICATION, same decision as above.

    A `SyntaxError` is the other purely *structural* failure: like `NameError` it
    is a property of the rule string alone, so it fails identically for every
    record pair and can only ever be a configuration defect. `eval` raises it at
    compile time, before any record value is read.
    """
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cross_match")

    with pytest.raises(ValueError) as excinfo:
        cross_match_variants(
            kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
            advntr_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
            config={"cross_match": {"match_logic": "Kestrel_Allele_Change =="}},
        )

    message = str(excinfo.value)
    assert "Kestrel_Allele_Change ==" in message
    assert "SyntaxError" in message
    assert "Advntr_Allele_Change" in message
    # Same classification as the NameError above, and for the same reason.
    assert message.startswith("Invalid cross_match match_logic")
    assert "Record:" not in message
    assert message in [r.message for r in caplog.records]


@pytest.mark.parametrize(
    ("match_logic", "pos"),
    [
        pytest.param("Kestrel_POS.lower() == Advntr_POS.lower()", 67, id="AttributeError-on-a-numeric-POS"),
        pytest.param("Kestrel_POS > 1", "67", id="TypeError-comparing-a-string-POS-to-an-int"),
    ],
)
def test_a_rule_that_cannot_be_evaluated_against_the_record_values_raises(match_logic, pos, caplog):
    """SPECIFICATION, same decision as above.

    These are the *value*-dependent failures: the rule parses and every name it
    uses exists, but the values in this record pair do not support the operation.
    They take a separate handler because the useful diagnostic is different -- the
    record pair, not just the column list -- but they raise for the same reason:
    a pair that could not be compared is not evidence of discordance.
    """
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cross_match")

    with pytest.raises(ValueError) as excinfo:
        cross_match_variants(
            kestrel_records=[{"REF": "C", "ALT": "CC", "POS": pos}],
            advntr_records=[{"REF": "C", "ALT": "CC", "POS": pos}],
            config={"cross_match": {"match_logic": match_logic}},
        )

    message = str(excinfo.value)
    assert match_logic in message
    assert "Kestrel_POS" in message
    # Classified as a bad record pair, not as a bad rule, and carries the pair itself --
    # this is what stops the two handlers being quietly collapsed back into one.
    assert message.startswith("Could not evaluate cross_match match_logic")
    assert "Record: {" in message
    assert "'Kestrel_Variant_Type': 'Insertion'" in message
    assert message in [r.message for r in caplog.records]


def test_the_default_match_logic_is_used_when_no_config_is_given():
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


# Every way the shipped tree can arrive at a rule. `vntyper/config.json` carries
# exactly one key under `cross_match` and its value is byte-identical to
# DEFAULT_MATCH_LOGIC -- pinned by the test immediately below -- so this is the
# complete set of correct configurations, not a sample of it.
_SHIPPED_CONFIGURATIONS = [
    pytest.param(None, id="config-is-None"),
    pytest.param({}, id="config-without-a-cross_match-key"),
    pytest.param({"cross_match": {}}, id="cross_match-without-a-match_logic-key"),
    pytest.param({"cross_match": {"match_logic": DEFAULT_MATCH_LOGIC}}, id="explicit-DEFAULT_MATCH_LOGIC"),
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


def test_every_shipped_match_logic_configuration_is_the_default_rule():
    """The enumeration above is complete only while this holds.

    `vntyper/config.json` must carry exactly one `cross_match` key, and its
    `match_logic` must be byte-identical to DEFAULT_MATCH_LOGIC. Adding a second
    rule, or diverging the shipped one from the default, fails here rather than
    leaving an unexercised configuration behind.
    """
    shipped = _shipped_cross_match_config()

    assert sorted(shipped) == ["match_logic"]
    assert shipped["match_logic"] == DEFAULT_MATCH_LOGIC


@pytest.mark.parametrize("config", _SHIPPED_CONFIGURATIONS)
@pytest.mark.parametrize(("kestrel_records", "advntr_records", "expected"), _MATCH_SCENARIOS)
def test_a_rule_that_evaluates_is_unaffected_by_the_fail_closed_change(
    config, kestrel_records, advntr_records, expected
):
    """SPECIFICATION. The expected values were captured from the code as it stood
    BEFORE the eval was made fail-closed, so this table passing means the fix
    changed nothing for a configuration that works.

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
