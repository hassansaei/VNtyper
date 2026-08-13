"""The screening interpretation: state in, configured sentence out.

This is the most consequential text VNtyper produces -- one sentence telling a
reader whether a sample screened positive. It is assembled from three axes
(Kestrel result, adVNTR result, quality gate) and a rule table in
``report_config.json``, and until now nothing tested either the axes or the
table.

Two failure modes are covered here that a per-rule test would miss:

* **Emphasis derived from the sentence.** The template decided "is this a
  positive finding?" by testing ``'negative' not in summary_text``. Not one of
  the configured messages contains that word, so every configured message
  rendered as a positive finding -- including "No variant detected by either
  genotyping method". :func:`is_finding` derives it from the state instead, and
  the tests below tie that derivation to the shipped configuration so a renamed
  ``result`` cannot silently invert it.
* **States with no rule.** A state the table does not cover falls through to
  ``screening_summary_default``, which reads "The screening was negative". The
  matrix test below walks the full cross product of reachable states and fails
  on any positive state that falls through.

The interpretive text itself is config-driven and stays that way (AGENTS.md): no
assertion here pins the wording of a message, only which state gets one.
"""

import itertools
import logging
from unittest import mock

import pandas as pd
import pytest

from vntyper.scripts import screening_summary as ss
from vntyper.scripts import summary_steps
from vntyper.scripts.coverage_qc import evaluate_coverage_qc

pytestmark = pytest.mark.unit


@pytest.fixture(scope="module")
def report_config() -> dict:
    """The shipped ``report_config.json``."""
    return ss.load_report_config()


def kestrel_logic(report_config) -> dict:
    return report_config["algorithm_logic"]["kestrel"]


def advntr_logic(report_config) -> dict:
    return report_config["algorithm_logic"]["advntr"]


def declared_results(logic: dict) -> list[str]:
    """Every ``result`` a rule block can produce, plus its ``default``."""
    return [rule["result"] for rule in logic["rules"]] + [logic["default"]]


# ---------------------------------------------------------------------------
# The shipped configuration is well-formed enough for the tests to mean anything
# ---------------------------------------------------------------------------


def test_the_report_config_loads(report_config) -> None:
    """Guard the guard: `load_report_config` swallows every exception and returns
    {}, which would make every assertion in this file vacuous."""
    assert report_config, "report_config.json did not load"
    assert report_config["screening_summary_rules"], "no screening rules loaded"
    assert report_config["screening_summary_default"]


@pytest.mark.parametrize("failure", ["open", "json_load"])
def test_report_config_failure_returns_empty_mapping(monkeypatch, caplog, failure: str) -> None:
    """Unreadable and invalid report configuration both remain visible empty fallbacks."""
    if failure == "open":
        monkeypatch.setattr("builtins.open", mock.Mock(side_effect=OSError("unreadable")))
    else:
        monkeypatch.setattr(ss.json, "load", mock.Mock(side_effect=ValueError("invalid JSON")))

    with caplog.at_level(logging.ERROR, logger="vntyper.scripts.screening_summary"):
        assert ss.load_report_config() == {}

    records = [record for record in caplog.records if record.name == "vntyper.scripts.screening_summary"]
    assert [record.levelno for record in records] == [logging.ERROR]
    assert "Failed to load report config" in records[0].message


def test_the_kestrel_block_declares_the_expected_states(report_config) -> None:
    assert set(declared_results(kestrel_logic(report_config))) == {
        "High_Precision",
        "High_Precision_flagged",
        "Low_Precision",
        "Low_Precision_flagged",
        "negative",
    }


def test_the_advntr_block_declares_the_expected_states(report_config) -> None:
    assert set(declared_results(advntr_logic(report_config))) == {
        "positive",
        "positive flagged",
        "negative",
    }


# ---------------------------------------------------------------------------
# compute_algorithm_result
# ---------------------------------------------------------------------------


def test_an_empty_frame_yields_the_default(report_config) -> None:
    assert compute_kestrel(pd.DataFrame(), report_config) == "negative"


def compute_kestrel(frame, report_config) -> str:
    return ss.compute_algorithm_result(frame, kestrel_logic(report_config))


def compute_advntr(frame, report_config) -> str:
    return ss.compute_algorithm_result(frame, advntr_logic(report_config))


@pytest.mark.parametrize(
    "confidence,flag,expected",
    [
        ("High_Precision", "Not flagged", "High_Precision"),
        ("High_Precision*", "Not flagged", "High_Precision"),
        ("Low_Precision", "Not flagged", "Low_Precision"),
        ("High_Precision", "Low_Depth_Conserved_Motifs", "High_Precision_flagged"),
        ("High_Precision*", "Duplicate", "High_Precision_flagged"),
        ("Low_Precision", "Duplicate", "Low_Precision_flagged"),
        ("Negative", "Not flagged", "negative"),
    ],
)
def test_the_kestrel_state_comes_from_confidence_and_flag(confidence, flag, expected, report_config) -> None:
    """`High_Precision*` is a distinct label the config lists explicitly; a rule
    that only tested `High_Precision` would silently downgrade the starred calls."""
    frame = pd.DataFrame([{"Confidence": confidence, "Flag": flag}])
    assert compute_kestrel(frame, report_config) == expected


def test_only_the_first_row_decides(report_config) -> None:
    """The frame is sorted by depth score, so the leading row is the strongest
    call. A second, weaker row must not change the sample-level state."""
    frame = pd.DataFrame(
        [
            {"Confidence": "High_Precision", "Flag": "Not flagged"},
            {"Confidence": "Low_Precision", "Flag": "Duplicate"},
        ]
    )
    assert compute_kestrel(frame, report_config) == "High_Precision"


@pytest.mark.parametrize(
    "vid,flag,expected",
    [
        ("25561", "Not flagged", "positive"),
        ("25561", "Repeat_Unit_7", "positive flagged"),
        ("Negative", "Not flagged", "negative"),
        ("25561", "Not applicable", "negative"),
    ],
)
def test_the_advntr_state_comes_from_vid_and_flag(vid, flag, expected, report_config) -> None:
    frame = pd.DataFrame([{"VID": vid, "Flag": flag}])
    assert compute_advntr(frame, report_config) == expected


def test_a_missing_column_fails_the_rule_rather_than_raising(report_config) -> None:
    """A frame with no `Flag` column reaches every Kestrel rule and matches none."""
    frame = pd.DataFrame([{"Confidence": "High_Precision"}])
    assert compute_kestrel(frame, report_config) == "negative"


def test_an_unsupported_operator_fails_the_rule() -> None:
    logic = {"rules": [{"conditions": {"C": {"operator": "~=", "value": "x"}}, "result": "hit"}], "default": "miss"}
    assert ss.compute_algorithm_result(pd.DataFrame([{"C": "x"}]), logic) == "miss"


def test_an_implicit_equality_condition_works() -> None:
    """The config supports a bare value as well as an explicit operator mapping."""
    logic = {"rules": [{"conditions": {"C": "x"}, "result": "hit"}], "default": "miss"}
    assert ss.compute_algorithm_result(pd.DataFrame([{"C": "x"}]), logic) == "hit"


def test_an_implicit_membership_condition_works() -> None:
    logic = {"rules": [{"conditions": {"C": ["x", "y"]}, "result": "hit"}], "default": "miss"}
    assert ss.compute_algorithm_result(pd.DataFrame([{"C": "y"}]), logic) == "hit"


# ---------------------------------------------------------------------------
# is_finding -- defect W3's replacement for the substring test
# ---------------------------------------------------------------------------


def test_the_word_negative_appears_in_no_configured_message(report_config) -> None:
    """The evidence for defect W3, pinned. The template tested
    `'negative' not in summary_text` and styled the box as a positive finding
    when that held. Zero of the configured messages contain the word, so the
    test was always true -- only the fallback default contains it."""
    messages = [rule["message"] for rule in report_config["screening_summary_rules"]]
    assert messages, "no messages loaded; this assertion would be vacuous"
    assert not [m for m in messages if "negative" in m.lower()]
    assert "negative" in report_config["screening_summary_default"].lower()


@pytest.mark.parametrize("algorithm", ["kestrel", "advntr"])
def test_every_declared_rule_result_counts_as_a_finding(algorithm, report_config) -> None:
    """Ties `is_finding` to the shipped configuration. Renaming a `result` in
    report_config.json without thinking about emphasis fails here rather than
    silently styling a positive report as negative."""
    logic = report_config["algorithm_logic"][algorithm]
    assert logic["rules"], "no rules loaded; this assertion would be vacuous"
    for rule in logic["rules"]:
        assert ss.is_finding(rule["result"], logic["default"]), f"{rule['result']!r} does not read as a finding"


@pytest.mark.parametrize("algorithm", ["kestrel", "advntr"])
def test_the_declared_default_is_not_a_finding(algorithm, report_config) -> None:
    logic = report_config["algorithm_logic"][algorithm]
    assert not ss.is_finding(logic["default"], logic["default"])


def test_an_algorithm_that_did_not_run_is_not_a_finding() -> None:
    """`none` means adVNTR was never asked, which is not a negative result and
    certainly not a positive one."""
    assert not ss.is_finding(ss.NOT_PERFORMED, "negative")


# ---------------------------------------------------------------------------
# The state matrix -- defect W4
# ---------------------------------------------------------------------------


def reachable_states(report_config) -> list[dict]:
    """Every screening state the pipeline can actually produce.

    adVNTR carries one state the rule blocks do not declare: ``none``, recorded
    when the stage did not run at all.
    """
    kestrel_states = declared_results(kestrel_logic(report_config))
    advntr_states = declared_results(advntr_logic(report_config)) + [ss.NOT_PERFORMED]
    return [
        {"kestrel_result": k, "advntr_result": a, "quality_metrics_pass": q}
        for k, a, q in itertools.product(kestrel_states, advntr_states, (True, False))
    ]


def test_the_state_matrix_is_not_empty(report_config) -> None:
    """Guard the guard: an empty matrix passes the coverage test vacuously."""
    states = reachable_states(report_config)
    assert len(states) == 5 * 4 * 2 == 40


def state_is_positive(state, report_config) -> bool:
    return ss.is_finding(state["kestrel_result"], kestrel_logic(report_config)["default"]) or ss.is_finding(
        state["advntr_result"], advntr_logic(report_config)["default"]
    )


def test_every_positive_state_has_its_own_message(report_config) -> None:
    """Defect W4. `positive flagged` is produced by the adVNTR rule block and
    tested for by none of the screening rules, so a Kestrel-positive,
    adVNTR-positive-but-flagged sample fell through to the default and rendered
    "The screening was negative". Agent E has just made a previously-dead adVNTR
    flagging rule able to fire, which moves more samples into that state."""
    uncovered = [
        state
        for state in reachable_states(report_config)
        if state_is_positive(state, report_config) and ss.find_screening_rule(report_config, state) is None
    ]
    assert not uncovered, "positive states with no configured message:\n" + "\n".join(
        f"  kestrel={s['kestrel_result']!r} advntr={s['advntr_result']!r} qc_pass={s['quality_metrics_pass']}"
        for s in uncovered
    )


def test_every_reachable_state_has_its_own_message(report_config) -> None:
    """Stronger than the contract requires, and cheaper to reason about: with a
    complete table, `screening_summary_default` is unreachable and its appearance
    in a report means the configuration is broken, not that the state is rare."""
    uncovered = [
        state for state in reachable_states(report_config) if ss.find_screening_rule(report_config, state) is None
    ]
    assert not uncovered, "states with no configured message:\n" + "\n".join(
        f"  kestrel={s['kestrel_result']!r} advntr={s['advntr_result']!r} qc_pass={s['quality_metrics_pass']}"
        for s in uncovered
    )


def test_no_two_rules_cover_the_same_state(report_config) -> None:
    """The first match wins, so a duplicated state makes a later rule dead."""
    seen: dict[tuple, int] = {}
    duplicates = []
    for index, rule in enumerate(report_config["screening_summary_rules"]):
        key = tuple(sorted((k, str(v)) for k, v in rule["conditions"].items()))
        if key in seen:
            duplicates.append((seen[key], index, key))
        seen[key] = index
    assert not duplicates, f"duplicate screening rule conditions: {duplicates}"


def test_every_rule_has_a_non_empty_message(report_config) -> None:
    """An empty message counts as "no rule matched" and silently falls back."""
    for index, rule in enumerate(report_config["screening_summary_rules"]):
        assert rule.get("message"), f"screening rule {index} has no message"


# ---------------------------------------------------------------------------
# ScreeningSummary.emphasis -- how the report should style the state
# ---------------------------------------------------------------------------


def _summary(
    *,
    is_positive: bool,
    matched_rule: bool,
    quality_metrics_pass: bool = True,
    text: str = "message",
    kestrel_result: str = "High_Precision",
    advntr_result: str = ss.NOT_PERFORMED,
    kestrel_execution: str = ss.EXECUTION_PERFORMED,
    advntr_execution: str = ss.EXECUTION_PERFORMED,
) -> "ss.ScreeningSummary":
    """Build a `ScreeningSummary` with only the axes each test cares about named."""
    return ss.ScreeningSummary(
        text=text,
        is_positive=is_positive,
        kestrel_result=kestrel_result,
        advntr_result=advntr_result,
        quality_metrics_pass=quality_metrics_pass,
        matched_rule=matched_rule,
        kestrel_execution=kestrel_execution,
        advntr_execution=advntr_execution,
    )


def test_an_unmatched_rule_is_indeterminate_not_negative() -> None:
    """`matched_rule` is the reason this property exists: a state with no configured
    rule is unknown, not a negative (see the module docstring)."""
    summary = _summary(is_positive=False, matched_rule=False)
    assert summary.emphasis == "indeterminate"


def test_an_unmatched_rule_stays_indeterminate_even_when_positive() -> None:
    """The other ordering: `matched_rule is False` wins over `is_positive` too.

    Both algorithms called positive but no rule covers the combination is exactly the
    case `matched_rule` exists to catch -- it must not be reported as a finding just
    because `is_positive` happens to be True.
    """
    summary = _summary(is_positive=True, matched_rule=False)
    assert summary.emphasis == "indeterminate"


def test_a_matched_positive_rule_is_a_finding() -> None:
    summary = _summary(is_positive=True, matched_rule=True)
    assert summary.emphasis == "finding"


def test_a_matched_negative_rule_is_a_no_finding() -> None:
    summary = _summary(is_positive=False, matched_rule=True)
    assert summary.emphasis == "no-finding"


def test_a_finding_with_failing_quality_metrics_is_still_a_finding() -> None:
    """QC is orthogonal to emphasis. It never suppresses a call.

    The rule table describes exactly this combination as a pathogenic finding with
    low-quality metrics (report_config.json), and ``is_positive`` is derived from the
    algorithm calls independently of QC by design (screening_summary.py). An earlier
    draft of this plan let failed QC force "indeterminate", which would have silently
    reclassified a confirmed pathogenic call with poor coverage as "state unknown" --
    contradicting the rule table, which describes exactly that combination as a
    finding with low-quality metrics, and contradicting ``is_positive`` itself.
    """
    summary = _summary(is_positive=True, matched_rule=True, quality_metrics_pass=False)
    assert summary.emphasis == "finding"


def test_a_no_finding_with_failing_quality_metrics_is_still_a_no_finding() -> None:
    """The mirror case: failing QC does not manufacture a finding either."""
    summary = _summary(is_positive=False, matched_rule=True, quality_metrics_pass=False)
    assert summary.emphasis == "no-finding"


# ---------------------------------------------------------------------------
# build_screening_summary, driven the way the report drives it
# ---------------------------------------------------------------------------


def _passing():
    """A coverage verdict that clears both shipped thresholds.

    The quality axis takes a :class:`~vntyper.scripts.coverage_qc.CoverageQC` since
    #172, not a mean and a threshold, so the tests below that are about something
    other than coverage say so by handing it a verdict that passes.
    """
    return evaluate_coverage_qc(250.0, 0.0, 100, 50.0)


KESTREL_FRAMES = {
    "High_Precision": [{"Confidence": "High_Precision", "Flag": "Not flagged"}],
    "High_Precision_flagged": [{"Confidence": "High_Precision", "Flag": "Duplicate"}],
    "Low_Precision": [{"Confidence": "Low_Precision", "Flag": "Not flagged"}],
    "Low_Precision_flagged": [{"Confidence": "Low_Precision", "Flag": "Duplicate"}],
    "negative": [],
}

ADVNTR_FRAMES = {
    "positive": [{"VID": "25561", "Flag": "Not flagged"}],
    "positive flagged": [{"VID": "25561", "Flag": "Repeat_Unit_7"}],
    "negative": [],
}


def test_the_frames_produce_the_states_they_claim(report_config) -> None:
    """Guard the guard: if a frame did not produce its state, the end-to-end
    matrix below would be testing the wrong thing."""
    for state, rows in KESTREL_FRAMES.items():
        assert compute_kestrel(pd.DataFrame(rows), report_config) == state
    for state, rows in ADVNTR_FRAMES.items():
        assert compute_advntr(pd.DataFrame(rows), report_config) == state


@pytest.mark.parametrize(
    "kestrel_state,advntr_state,coverage",
    list(itertools.product(KESTREL_FRAMES, list(ADVNTR_FRAMES) + [ss.NOT_PERFORMED], (250.0, 3.0))),
)
def test_every_state_reached_through_the_real_frames_gets_a_rule(
    kestrel_state, advntr_state, coverage, report_config
) -> None:
    """The same matrix, driven end to end rather than against the config, so a
    state the config covers but the code cannot reach is caught too."""
    advntr_available = advntr_state != ss.NOT_PERFORMED
    summary = ss.build_screening_summary(
        pd.DataFrame(KESTREL_FRAMES[kestrel_state]),
        pd.DataFrame(ADVNTR_FRAMES[advntr_state]) if advntr_available else pd.DataFrame(),
        advntr_available,
        evaluate_coverage_qc(coverage, 0.0, 100, 50.0),
        report_config,
    )

    assert summary.kestrel_result == kestrel_state
    assert summary.advntr_result == advntr_state
    assert summary.quality_metrics_pass is (coverage >= 100)
    assert summary.matched_rule, (
        f"kestrel={kestrel_state!r} advntr={advntr_state!r} qc_pass={coverage >= 100} "
        "fell through to screening_summary_default"
    )
    assert summary.text != report_config["screening_summary_default"]


def test_a_kestrel_negative_advntr_negative_sample_is_not_a_finding(report_config) -> None:
    """The case defect W3 got backwards: "No variant detected by either
    genotyping method" was styled as a positive finding."""
    summary = ss.build_screening_summary(pd.DataFrame(), pd.DataFrame(), True, _passing(), report_config)
    assert summary.is_positive is False


@pytest.mark.parametrize("kestrel_state", ["High_Precision", "High_Precision_flagged", "Low_Precision_flagged"])
def test_a_kestrel_finding_is_positive(kestrel_state, report_config) -> None:
    summary = ss.build_screening_summary(
        pd.DataFrame(KESTREL_FRAMES[kestrel_state]), pd.DataFrame(), False, _passing(), report_config
    )
    assert summary.is_positive is True


def test_an_advntr_only_finding_is_positive(report_config) -> None:
    summary = ss.build_screening_summary(
        pd.DataFrame(), pd.DataFrame(ADVNTR_FRAMES["positive"]), True, _passing(), report_config
    )
    assert summary.is_positive is True


def test_a_flagged_advntr_only_finding_is_positive(report_config) -> None:
    summary = ss.build_screening_summary(
        pd.DataFrame(), pd.DataFrame(ADVNTR_FRAMES["positive flagged"]), True, _passing(), report_config
    )
    assert summary.is_positive is True


def test_an_advntr_that_never_ran_is_not_a_finding(report_config) -> None:
    summary = ss.build_screening_summary(pd.DataFrame(), pd.DataFrame(), False, _passing(), report_config)
    assert summary.advntr_result == ss.NOT_PERFORMED
    assert summary.is_positive is False


def test_coverage_below_the_threshold_fails_the_quality_gate(report_config) -> None:
    summary = ss.build_screening_summary(
        pd.DataFrame(), pd.DataFrame(), False, evaluate_coverage_qc(99.0, 0.0, 100, 50.0), report_config
    )
    assert summary.quality_metrics_pass is False


def test_coverage_that_was_never_measured_passes_the_quality_gate(report_config) -> None:
    """Pinned rather than endorsed: an unmeasured sample is reported as passing
    the coverage gate. Changing that changes a displayed interpretation for
    every run with no Coverage Calculation step, so it is left alone here."""
    summary = ss.build_screening_summary(
        pd.DataFrame(), pd.DataFrame(), False, evaluate_coverage_qc(None, None, 100, 50.0), report_config
    )
    assert summary.quality_metrics_pass is True


def test_a_patchy_sample_fails_the_screening_quality_axis(report_config) -> None:
    """#172: acceptable mean, half the VNTR uncovered. Before this change it passed."""
    qc = evaluate_coverage_qc(250.0, 80.0, 100, 50.0)
    summary = ss.build_screening_summary(pd.DataFrame(), pd.DataFrame(), False, qc, report_config)

    assert summary.quality_metrics_pass is False


def test_widening_the_quality_axis_cannot_change_positivity(report_config) -> None:
    """`is_positive` derives from the two algorithm axes only, never from quality.

    Stated as a test because it is the one way #172 could have become
    genotype-affecting, and it is not.
    """
    failing = evaluate_coverage_qc(1.0, 99.0, 100, 50.0)
    passing = evaluate_coverage_qc(250.0, 1.0, 100, 50.0)

    assert (
        ss.build_screening_summary(pd.DataFrame(), pd.DataFrame(), False, failing, report_config).is_positive
        is ss.build_screening_summary(pd.DataFrame(), pd.DataFrame(), False, passing, report_config).is_positive
    )


def test_a_broken_config_yields_the_unavailable_message(caplog) -> None:
    """An internal screening dependency failure yields the explicit unavailable state.

    This is also the only path through `build_screening_summary` that leaves
    `matched_rule` False in practice: all 40 reachable (kestrel_result, advntr_result,
    quality_metrics_pass) combinations resolve to a configured rule (see
    `test_every_reachable_state_has_its_own_message`), so `emphasis == "indeterminate"`
    is genuinely exceptional and cannot mislabel an ordinary all-negative report.
    """
    with (
        mock.patch.object(ss, "compute_algorithm_result", side_effect=RuntimeError("boom")),
        caplog.at_level(logging.ERROR, logger="vntyper.scripts.screening_summary"),
    ):
        summary = ss.build_screening_summary(pd.DataFrame(), pd.DataFrame(), False, _passing(), {})

    assert summary.text == ss.UNAVAILABLE_SUMMARY_MESSAGE
    assert summary.is_positive is False
    assert summary.kestrel_result == ""
    assert summary.advntr_result == ""
    assert summary.quality_metrics_pass is False
    assert summary.matched_rule is False
    assert summary.emphasis == "indeterminate"
    records = [record for record in caplog.records if record.name == "vntyper.scripts.screening_summary"]
    assert [record.levelno for record in records] == [logging.ERROR]


def test_a_state_with_no_rule_warns_rather_than_going_quiet(caplog, report_config) -> None:
    """The fallback is indistinguishable from a real negative in the HTML, so
    the only trace it leaves is this log line."""
    stripped = {**report_config, "screening_summary_rules": []}
    with caplog.at_level(logging.WARNING, logger="vntyper.scripts.screening_summary"):
        summary = ss.build_screening_summary(
            pd.DataFrame(KESTREL_FRAMES["High_Precision"]), pd.DataFrame(), False, _passing(), stripped
        )

    assert summary.matched_rule is False
    assert summary.text == report_config["screening_summary_default"]
    assert summary.is_positive is True, "the state is still positive even with no message for it"
    assert any(record.levelno >= logging.WARNING for record in caplog.records)


# ---------------------------------------------------------------------------
# Execution state -- a stage that ran and produced nothing readable is not a
# negative, and the message chosen for the state it did not establish is not
# authoritative (#242, Gate B)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("step_state", "expected"),
    [
        (summary_steps.STEP_READ, ss.EXECUTION_PERFORMED),
        (summary_steps.STEP_ABSENT, ss.EXECUTION_NOT_PERFORMED),
        (summary_steps.STEP_UNREADABLE, ss.EXECUTION_FAILED),
    ],
)
def test_each_step_state_maps_to_one_execution_state(step_state, expected) -> None:
    """``get_step_state`` already distinguishes three states; this is the only
    translation of them into the vocabulary the report presents."""
    assert ss.execution_state(step_state) == expected


def test_an_unrecognised_step_state_is_treated_as_a_failure() -> None:
    """The conservative direction. A state this module does not recognise has
    established nothing about the sample, so it must not be presented as a result."""
    assert ss.execution_state("something-new") == ss.EXECUTION_FAILED


#: One failed stage each, with the other reporting normally.
ONE_STAGE_FAILED = [
    (ss.EXECUTION_FAILED, ss.EXECUTION_PERFORMED),
    (ss.EXECUTION_PERFORMED, ss.EXECUTION_FAILED),
]


@pytest.mark.parametrize(("kestrel", "advntr"), ONE_STAGE_FAILED, ids=["kestrel", "advntr"])
def test_a_failed_stage_makes_the_screening_state_indeterminate(kestrel: str, advntr: str) -> None:
    """Either algorithm is enough. A rule was selected by matching a state one of
    the two stages never established, so the state as a whole is not established."""
    summary = _summary(is_positive=False, matched_rule=True, kestrel_execution=kestrel, advntr_execution=advntr)
    assert summary.emphasis == "indeterminate"


@pytest.mark.parametrize(("kestrel", "advntr"), ONE_STAGE_FAILED, ids=["kestrel", "advntr"])
def test_a_failed_stage_outranks_a_positive_call(kestrel: str, advntr: str) -> None:
    """The other ordering, and the one that matters: a Kestrel call cannot make a
    report authoritative about a state adVNTR never established, and vice versa."""
    summary = _summary(is_positive=True, matched_rule=True, kestrel_execution=kestrel, advntr_execution=advntr)
    assert summary.emphasis == "indeterminate"


def test_an_advntr_stage_that_was_never_asked_to_run_is_not_a_failure() -> None:
    """An adVNTR-less run is the commonest run there is, and the rule table covers it.

    Eight of the forty rules are keyed on ``advntr_result == "none"`` and say in words
    that the stage was not performed, so the state *is* established. Treating this like a
    failure would make every adVNTR-less report indeterminate.
    """
    summary = _summary(is_positive=True, matched_rule=True, advntr_execution=ss.EXECUTION_NOT_PERFORMED)
    assert summary.state_is_established is True
    assert summary.emphasis == "finding"


def test_an_absent_kestrel_stage_leaves_the_state_unestablished() -> None:
    """The asymmetry, and where it comes from.

    ``kestrel_result`` has no value for "did not run": an absent Kestrel step hands
    ``compute_algorithm_result`` an empty frame and gets the block's ``default`` back -
    the same ``negative`` a stage that genotyped and called nothing produces. So the rule
    that matches is one keyed on a Kestrel negative, and the report would state a negative
    for a stage that never ran. ``vntyper report`` can legitimately be handed such a
    summary (#207), so this is reachable rather than theoretical.
    """
    summary = _summary(is_positive=False, matched_rule=True, kestrel_execution=ss.EXECUTION_NOT_PERFORMED)
    assert summary.state_is_established is False
    assert summary.emphasis == "indeterminate"


def test_a_failed_kestrel_stage_withholds_the_configured_message(report_config, caplog) -> None:
    """The Gate B defect, at its source.

    ``record_step`` writes an empty ``data`` list for a Kestrel stage whose result
    file is missing (#212) -- the same shape a run that genotyped and called nothing
    produces -- so ``compute_algorithm_result`` returns the block's ``default`` and the
    rule table hands back "No variant detected." A report that prints that sentence is
    asserting a negative the run never reached.
    """
    with caplog.at_level(logging.WARNING, logger="vntyper.scripts.screening_summary"):
        summary = ss.build_screening_summary(
            pd.DataFrame(),
            pd.DataFrame(),
            False,
            _passing(),
            report_config,
            kestrel_execution=ss.EXECUTION_FAILED,
        )

    assert summary.text == ss.UNAVAILABLE_SUMMARY_MESSAGE
    assert summary.emphasis == "indeterminate"
    assert summary.kestrel_execution == ss.EXECUTION_FAILED
    assert any(record.levelno >= logging.WARNING for record in caplog.records)


def test_a_failed_advntr_stage_withholds_the_configured_message(report_config) -> None:
    """The half the first fix missed: only Kestrel's step state reached the report,
    so an unreadable adVNTR stage was reported as "No pathogenic variants identified"."""
    summary = ss.build_screening_summary(
        pd.DataFrame(KESTREL_FRAMES["High_Precision"]),
        pd.DataFrame(),
        False,
        _passing(),
        report_config,
        advntr_execution=ss.EXECUTION_FAILED,
    )

    assert summary.text == ss.UNAVAILABLE_SUMMARY_MESSAGE
    assert summary.emphasis == "indeterminate"


def test_a_run_that_states_its_stages_ran_keeps_its_configured_message(report_config) -> None:
    """Guard the guard: withholding must be reachable *and* avoidable, or the two
    tests above would pass against a module that never renders a message at all."""
    summary = ss.build_screening_summary(
        pd.DataFrame(KESTREL_FRAMES["High_Precision"]),
        pd.DataFrame(),
        False,
        _passing(),
        report_config,
        kestrel_execution=ss.EXECUTION_PERFORMED,
        advntr_execution=ss.EXECUTION_NOT_PERFORMED,
    )

    assert summary.text != ss.UNAVAILABLE_SUMMARY_MESSAGE
    assert summary.emphasis == "finding"


def test_the_execution_states_default_to_the_two_state_model_they_replace(report_config) -> None:
    """A caller that says nothing gets exactly the model that existed before: Kestrel
    ran, and adVNTR ran if and only if ``advntr_available`` says so. No default asserts
    anything the old signature could not."""
    performed = ss.build_screening_summary(pd.DataFrame(), pd.DataFrame(), True, _passing(), report_config)
    absent = ss.build_screening_summary(pd.DataFrame(), pd.DataFrame(), False, _passing(), report_config)

    assert performed.kestrel_execution == absent.kestrel_execution == ss.EXECUTION_PERFORMED
    assert performed.advntr_execution == ss.EXECUTION_PERFORMED
    assert absent.advntr_execution == ss.EXECUTION_NOT_PERFORMED


@pytest.mark.parametrize(
    ("execution", "result", "expected"),
    [
        (ss.EXECUTION_PERFORMED, "High_Precision", "High_Precision"),
        (ss.EXECUTION_PERFORMED, "negative", "negative"),
        (ss.EXECUTION_NOT_PERFORMED, "none", "not performed"),
        (ss.EXECUTION_FAILED, "negative", "failed"),
    ],
)
def test_the_provenance_word_says_the_execution_state_before_the_result(execution, result, expected) -> None:
    """The raw state line prints the computed result only when there is one. The last
    row is the defect: an unreadable stage computes the block's ``default``, so a line
    that printed the result would read "Kestrel: negative"."""
    assert ss.algorithm_state_text(execution, result) == expected
