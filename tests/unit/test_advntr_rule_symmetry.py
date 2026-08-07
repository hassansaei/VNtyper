# tests/unit/test_advntr_rule_symmetry.py

"""Symmetry between adVNTR's two positive rules in `report_config.json`'s `algorithm_logic`.

Rule 1 (result `"positive"`) guards on both `VID != "Negative"` and `Flag == "Not
flagged"`. Rule 2 (result `"positive flagged"`) is meant to be rule 1's flagged
counterpart, but as originally shipped it guarded on `Flag` alone and never looked at
`VID` -- so, on paper, a `VID == "Negative"` row carrying any flag outside the three
excused values (`"Not flagged"`, `"Not applicable"`, `"None"`) came out `"positive
flagged"`.

That input was not producible by the pipeline: both sentinel-row constructions in
`advntr_genotyping.py` set `Flag = "Not applicable"`, which rule 2 excludes; the sentinel
never reaches `add_flags`; and even if it did, the shipped flagging rules cannot fire on
sentinel values. See
`.superpowers/sdd/2026-08-06-issue-181-197-followups-plan/issue-cohort-advntr-vid-column-overloaded.md`
for the full reachability argument. The asymmetry was real at the rule-table level even
though it was unreachable at the pipeline level -- this module closes it at the
rule-table level and proves no verdict the pipeline can actually produce moves.

This is a different subject from `tests/unit/test_cohort_rules.py`, which tests the rule
*interpreter* (`compute_algorithm_result`) in the general case. This module tests one
property of the shipped adVNTR *rule table*: its symmetry between rule 1 and rule 2.
"""

import json
from pathlib import Path

import pandas as pd
import pytest

from vntyper.modules.advntr import advntr_genotyping as advntr
from vntyper.scripts.cohort_rules import compute_algorithm_result
from vntyper.scripts.cohort_summary import load_report_config

pytestmark = pytest.mark.unit


def _row(**values: object) -> pd.DataFrame:
    """Build the single-row frame the interpreter is always called with.

    Args:
        **values: The row's columns.

    Returns:
        pd.DataFrame: A one-row frame.
    """
    return pd.DataFrame([values])


def _advntr_flagging_rule_names() -> list[str]:
    """The real flag names adVNTR's shipped flagging rules can write into `Flag`.

    `add_flags` (`vntyper/scripts/flagging.py`) uses each `flagging_rules` key verbatim as
    the flag name it appends when that rule's condition holds -- there is no `"_flagged"`
    suffix. Read from `advntr_config.json` on disk rather than hardcoded, per the brief:
    do not invent the flag strings.

    Returns:
        list[str]: The configured flag names, e.g. `["Low_Coverage", "Polymorphic_Call",
        "Repeat_Unit_7"]`.
    """
    config_path = Path(advntr.__file__).parent / "advntr_config.json"
    with config_path.open() as handle:
        config = json.load(handle)
    return sorted(config["flagging_rules"])


#: The real flag names adVNTR's shipped flagging rules can emit, read from
#: `advntr_config.json` rather than hardcoded.
ADVNTR_FLAG_NAMES = _advntr_flagging_rule_names()

#: `algorithm_logic.advntr` exactly as it shipped before this fix: rule 2 guards on `Flag`
#: alone, with no `VID` condition. Reconstructed here as a literal -- rather than read from
#: disk -- because `report_config.json` now carries the fix; this is what "before" means
#: for the neutrality proof below.
PRE_FIX_ADVNTR_RULES = {
    "rules": [
        {
            "conditions": {
                "VID": {"operator": "!=", "value": "Negative"},
                "Flag": {"operator": "==", "value": "Not flagged"},
            },
            "result": "positive",
        },
        {
            "conditions": {
                "Flag": {"operator": "not in", "value": ["Not flagged", "Not applicable", "None"]},
            },
            "result": "positive flagged",
        },
    ],
    "default": "negative",
}


# ---------------------------------------------------------------------------
# The defect: rule 2 was missing rule 1's VID guard
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("flag", ADVNTR_FLAG_NAMES)
def test_a_negative_vid_with_a_real_flag_is_not_positive_flagged(flag: str) -> None:
    """Specification: rule 2 must guard on `VID` the same way rule 1 does.

    Before the fix, a hand-built `VID == "Negative"` row carrying one of adVNTR's real
    flag names matched rule 2 on the flag alone and came out `"positive flagged"`. A
    sample adVNTR called negative must never be reported positive under any name, so this
    must not happen regardless of whether the pipeline can currently construct such a row.
    """
    logic = load_report_config()["algorithm_logic"]["advntr"]

    verdict = compute_algorithm_result(_row(VID="Negative", Flag=flag), logic)

    assert verdict != "positive flagged"
    assert verdict == "negative"


def test_the_advntr_rule_table_shape_is_unchanged() -> None:
    """Guardrail for "keep the rule order and the default unchanged": two rules, same
    results, same order, same default -- only rule 2's conditions may have grown a guard.
    """
    advntr_logic = load_report_config()["algorithm_logic"]["advntr"]

    assert [rule["result"] for rule in advntr_logic["rules"]] == ["positive", "positive flagged"]
    assert advntr_logic["default"] == "negative"


# ---------------------------------------------------------------------------
# The fix is genotype-neutral: every row the pipeline can actually produce
# ---------------------------------------------------------------------------

#: Every (VID, Flag) combination the real pipeline can hand to `compute_algorithm_result`
#: today, paired with the verdict it must produce both before and after the fix:
#:
#: * the sentinel row both empty-VCF and empty-after-filter paths build
#:   (`advntr_genotyping.py:390-405`, `:425-440`), `VID == "Negative"` with
#:   `Flag == "Not applicable"`;
#: * a real call (`VID == "25561"`, the hardcoded locus ID from `advntr_config.json`'s
#:   `advntr_settings.vid`) with `Flag == "Not flagged"`, i.e. no flagging rule fired;
#: * a real call with each flag `add_flags` can actually write, one per rule in
#:   `advntr_config.json`'s `flagging_rules`.
REACHABLE_ROWS = [
    pytest.param("Negative", "Not applicable", "negative", id="sentinel-row"),
    pytest.param("25561", "Not flagged", "positive", id="real-row-not-flagged"),
    *[pytest.param("25561", flag, "positive flagged", id=f"real-row-{flag}") for flag in ADVNTR_FLAG_NAMES],
]


@pytest.mark.parametrize("vid,flag,expected_verdict", REACHABLE_ROWS)
def test_every_pipeline_reachable_row_yields_the_same_verdict_before_and_after(
    vid: str, flag: str, expected_verdict: str
) -> None:
    """Proof that adding rule 2's `VID` guard changes nothing the pipeline can produce.

    Runs each reachable (VID, Flag) combination through `PRE_FIX_ADVNTR_RULES` (the table
    as it shipped, reconstructed above) and through the live, on-disk `report_config.json`
    (post-fix), and asserts both equal the recorded verdict. The added guard can only
    change rule 2's outcome for a `VID == "Negative"` row, and no reachable row is one --
    this table is that claim, enumerated rather than spot-checked.
    """
    row = _row(VID=vid, Flag=flag)
    live_logic = load_report_config()["algorithm_logic"]["advntr"]

    before = compute_algorithm_result(row, PRE_FIX_ADVNTR_RULES)
    after = compute_algorithm_result(row, live_logic)

    assert before == expected_verdict
    assert after == expected_verdict
