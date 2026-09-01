"""Safety-first lexicographic calibration objective."""

from dataclasses import replace
from fractions import Fraction
from typing import cast

import pytest

from tests.unit.test_calibration_contract import synthetic_protocol
from vntyper.scripts.calibration_contract import CandidateMetrics, decode_protocol
from vntyper.scripts.calibration_objective import (
    CandidateEvaluation,
    OutcomeObservation,
    calculate_metrics,
    count_free_parameters,
    lexicographic_safety_key,
    select_candidate,
)

pytestmark = pytest.mark.unit


def _metrics(**changes: object) -> CandidateMetrics:
    values: dict[str, object] = {
        "candidate_profile_sha256": "a" * 64,
        "wrong_tier_a_displayed_names": 0,
        "control_findings": 0,
        "wrong_displayed_names_all_tiers": 0,
        "macro_exact_recovery": Fraction(3, 4),
        "binary_detection_sensitivity": Fraction(4, 5),
        "free_parameter_count": 1,
        "abstention_fraction": Fraction(1, 10),
        "tier_a_reachable": True,
        "applicability_matches": True,
        "sha256": "b" * 64,
    }
    values.update(changes)
    values["wrong_displayed_names_all_tiers"] = max(
        cast(int, values["wrong_displayed_names_all_tiers"]),
        cast(int, values["wrong_tier_a_displayed_names"]),
    )
    return CandidateMetrics(**values)  # type: ignore[arg-type]


def _evaluation(metrics: CandidateMetrics, **changes: object) -> CandidateEvaluation:
    values: dict[str, object] = {
        "metrics": metrics,
        "detection_lower_bound": Fraction(0),
        "macro_exact_lower_bound": Fraction(0),
        "stratum_counts": (2, 2),
        "holm_adjusted_p_value": Fraction(1, 100),
    }
    values.update(changes)
    return CandidateEvaluation(**values)  # type: ignore[arg-type]


def _outcome(**changes: object) -> OutcomeObservation:
    values: dict[str, object] = {
        "key": "member-a",
        "assay_class": "assay-a",
        "mutation_class": "duplication",
        "expected_identity": "identity-a",
        "expected_display_name": "59dupC",
        "selected_identity": "identity-a",
        "displayed_name": "59dupC",
        "tier": "A",
        "abstained": False,
        "applicable": True,
        "baseline_applicable": True,
        "baseline_tier": "A",
    }
    values.update(changes)
    return OutcomeObservation(**values)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("earlier", "better_later"),
    [
        ("wrong_tier_a_displayed_names", {"control_findings": 0, "wrong_displayed_names_all_tiers": 0}),
        ("control_findings", {"wrong_displayed_names_all_tiers": 0}),
        ("wrong_displayed_names_all_tiers", {"macro_exact_recovery": Fraction(1)}),
        ("macro_exact_recovery", {"binary_detection_sensitivity": Fraction(1)}),
        ("binary_detection_sensitivity", {"free_parameter_count": 0}),
        ("free_parameter_count", {"candidate_profile_sha256": "0" * 64}),
    ],
)
def test_every_adjacent_objective_boundary_is_decided_by_the_earlier_component(
    earlier: str,
    better_later: dict[str, object],
) -> None:
    preferred = _metrics(candidate_profile_sha256="f" * 64)
    if earlier in {"macro_exact_recovery", "binary_detection_sensitivity"}:
        degraded_value: object = Fraction(1, 2)
    else:
        degraded_value = getattr(preferred, earlier) + 1
    degraded = replace(preferred, **{earlier: degraded_value}, **better_later)  # type: ignore[arg-type]

    assert lexicographic_safety_key(preferred) < lexicographic_safety_key(degraded)


def test_profile_hash_is_the_deterministic_final_tie_break() -> None:
    first = _metrics(candidate_profile_sha256="0" * 64)
    second = _metrics(candidate_profile_sha256="f" * 64)

    assert lexicographic_safety_key(first) < lexicographic_safety_key(second)


@pytest.mark.parametrize(
    "changes",
    [
        {"wrong_tier_a_displayed_names": 1},
        {"control_findings": 1},
        {"tier_a_reachable": False},
        {"applicability_matches": False},
        {"abstention_fraction": Fraction(1, 2)},
        {"free_parameter_count": 5},
    ],
)
def test_hard_safety_and_protocol_constraints_make_candidates_inadmissible(changes: dict[str, object]) -> None:
    protocol = decode_protocol(synthetic_protocol())

    assert select_candidate((_evaluation(_metrics(**changes)),), protocol) is None


def test_noninferiority_bounds_and_minimum_stratum_count_are_hard_constraints() -> None:
    protocol = decode_protocol(synthetic_protocol())
    metrics = _metrics()

    for evaluation in (
        _evaluation(metrics, detection_lower_bound=Fraction(-1, 100)),
        _evaluation(metrics, macro_exact_lower_bound=Fraction(-1, 100)),
        _evaluation(metrics, stratum_counts=(2, 1)),
        _evaluation(metrics, holm_adjusted_p_value=Fraction(6, 100)),
    ):
        assert select_candidate((evaluation,), protocol) is None


def test_one_wrong_tier_a_or_one_control_finding_never_wins_on_later_metrics() -> None:
    protocol = decode_protocol(synthetic_protocol())
    safe = _evaluation(_metrics(macro_exact_recovery=Fraction(1, 10), binary_detection_sensitivity=Fraction(1, 10)))
    wrong_tier_a = _evaluation(
        _metrics(
            candidate_profile_sha256="c" * 64,
            wrong_tier_a_displayed_names=1,
            macro_exact_recovery=Fraction(1),
            binary_detection_sensitivity=Fraction(1),
        )
    )
    control_finding = _evaluation(
        _metrics(
            candidate_profile_sha256="d" * 64,
            control_findings=1,
            macro_exact_recovery=Fraction(1),
            binary_detection_sensitivity=Fraction(1),
        )
    )

    assert select_candidate((wrong_tier_a, safe, control_finding), protocol) == safe


def test_all_abstain_is_inadmissible_and_abstention_stays_in_rate_denominators() -> None:
    protocol = decode_protocol(synthetic_protocol())
    all_abstain = _evaluation(
        _metrics(
            macro_exact_recovery=Fraction(0),
            binary_detection_sensitivity=Fraction(0),
            abstention_fraction=Fraction(1),
            tier_a_reachable=False,
        )
    )

    assert select_candidate((all_abstain,), protocol) is None


def test_free_parameter_count_uses_only_non_neutral_generated_rule_leaves() -> None:
    neutral = {
        "enabled": False,
        "minimum_record_count_margin": 0,
        "minimum_record_share": 0.0,
        "minimum_record_share_margin": 0.0,
        "xd_veto": "disabled",
        "abstain_on_inadmissible_advntr": False,
    }
    active = {
        **neutral,
        "enabled": True,
        "minimum_record_count_margin": 1,
        "minimum_record_share": 0.5,
        "minimum_record_share_margin": 0.25,
        "xd_veto": "missingness",
        "abstain_on_inadmissible_advntr": True,
    }

    assert count_free_parameters(neutral) == 0
    assert count_free_parameters(active) == 6


def test_no_candidate_is_emitted_when_all_violate_a_hard_constraint() -> None:
    protocol = decode_protocol(synthetic_protocol())
    candidates = (
        _evaluation(_metrics(wrong_tier_a_displayed_names=1)),
        _evaluation(_metrics(candidate_profile_sha256="c" * 64, control_findings=1)),
    )

    assert select_candidate(candidates, protocol) is None


def test_selective_abstention_remains_in_exact_and_detection_denominators() -> None:
    rows = (
        OutcomeObservation(
            "a",
            "assay-a",
            "duplication",
            "identity-a",
            "59dupC",
            "identity-a",
            "59dupC",
            "A",
            False,
            True,
            True,
            "A",
        ),
        OutcomeObservation(
            "b",
            "assay-a",
            "duplication",
            "identity-b",
            "60dupA",
            None,
            None,
            None,
            True,
            True,
            True,
            "B",
        ),
    )

    summary = calculate_metrics(
        rows, profile_sha256="a" * 64, free_parameter_count=1, required_strata=("assay-a:duplication",)
    )

    assert summary.metrics.macro_exact_recovery == Fraction(1, 2)
    assert summary.metrics.binary_detection_sensitivity == Fraction(1, 2)
    assert summary.metrics.abstention_fraction == Fraction(1, 2)


def test_wrong_identity_can_count_as_detection_while_display_error_is_independent() -> None:
    rows = (
        OutcomeObservation(
            "a",
            "assay-a",
            "duplication",
            "identity-a",
            "expected-name",
            "wrong-identity",
            "expected-name",
            "B",
            False,
            True,
            True,
            "B",
        ),
        OutcomeObservation(
            "b",
            "assay-a",
            "duplication",
            "identity-b",
            "expected-b",
            "identity-b",
            "wrong-name",
            "A",
            False,
            True,
            True,
            "A",
        ),
    )

    summary = calculate_metrics(
        rows, profile_sha256="a" * 64, free_parameter_count=1, required_strata=("assay-a:duplication",)
    )

    assert summary.metrics.binary_detection_sensitivity == Fraction(1)
    assert summary.metrics.macro_exact_recovery == Fraction(1, 2)
    assert summary.metrics.wrong_displayed_names_all_tiers == 1
    assert summary.metrics.wrong_tier_a_displayed_names == 1


def test_macro_exact_recovery_is_unweighted_across_imbalanced_assay_mutation_strata() -> None:
    rows = [
        OutcomeObservation(
            f"large-{index}",
            "assay-a",
            "duplication",
            f"identity-{index}",
            f"name-{index}",
            f"identity-{index}",
            f"name-{index}",
            "B",
            False,
            True,
            True,
            "B",
        )
        for index in range(9)
    ]
    rows.append(
        OutcomeObservation(
            "small",
            "assay-b",
            "deletion",
            "identity-small",
            "small-name",
            None,
            None,
            None,
            True,
            True,
            True,
            "B",
        )
    )

    summary = calculate_metrics(
        tuple(rows),
        profile_sha256="a" * 64,
        free_parameter_count=1,
        required_strata=("assay-a:duplication", "assay-b:deletion"),
    )

    assert summary.metrics.macro_exact_recovery == Fraction(1, 2)
    assert summary.metrics.binary_detection_sensitivity == Fraction(9, 10)
    assert summary.stratum_counts == (9, 1)


def test_controls_are_counted_separately_and_applicability_mismatch_is_retained() -> None:
    rows = (
        OutcomeObservation(
            "control", "assay-a", "control", None, None, "finding", "59dupC", "B", False, True, True, None
        ),
        OutcomeObservation(
            "mutated",
            "assay-a",
            "duplication",
            "identity",
            "59dupC",
            "identity",
            "59dupC",
            "A",
            False,
            False,
            True,
            "A",
        ),
    )

    summary = calculate_metrics(
        rows, profile_sha256="a" * 64, free_parameter_count=1, required_strata=("assay-a:duplication",)
    )

    assert summary.metrics.control_findings == 1
    assert not summary.metrics.applicability_matches
    assert summary.metrics.tier_a_reachable


def test_objective_value_objects_are_frozen() -> None:
    evaluation = _evaluation(_metrics())
    outcome = _outcome()
    summary = calculate_metrics(
        (outcome,), profile_sha256="a" * 64, free_parameter_count=0, required_strata=("assay-a:duplication",)
    )

    for value in (evaluation, outcome, summary):
        with pytest.raises(AttributeError):
            value.stratum_counts = ()  # type: ignore[attr-defined]


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"metrics": object()}, "CandidateMetrics"),
        ({"detection_lower_bound": 0}, "detection lower bound"),
        ({"macro_exact_lower_bound": 0}, "macro exact lower bound"),
        ({"holm_adjusted_p_value": 0}, "Holm-adjusted"),
        ({"stratum_counts": []}, "non-empty tuple"),
        ({"stratum_counts": ()}, "non-empty tuple"),
        ({"stratum_counts": (True,)}, "non-negative integers"),
        ({"stratum_counts": (1.0,)}, "non-negative integers"),
        ({"stratum_counts": (-1,)}, "non-negative integers"),
        ({"holm_adjusted_p_value": Fraction(-1, 100)}, "between zero and one"),
        ({"holm_adjusted_p_value": Fraction(101, 100)}, "between zero and one"),
    ],
)
def test_candidate_evaluation_rejects_each_invalid_contract_branch(changes: dict[str, object], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        values: dict[str, object] = {
            "metrics": _metrics(),
            "detection_lower_bound": Fraction(0),
            "macro_exact_lower_bound": Fraction(0),
            "stratum_counts": (2, 2),
            "holm_adjusted_p_value": Fraction(1, 100),
        }
        values.update(changes)
        CandidateEvaluation(**values)  # type: ignore[arg-type]


def test_candidate_evaluation_rejects_tier_a_errors_missing_from_all_tier_count() -> None:
    inconsistent = replace(_metrics(), wrong_tier_a_displayed_names=1, wrong_displayed_names_all_tiers=0)

    with pytest.raises(ValueError, match="all-tier wrong"):
        _evaluation(inconsistent)


def test_candidate_evaluation_accepts_zero_members_for_predeclared_stratum() -> None:
    evaluation = _evaluation(_metrics(), stratum_counts=(0,))

    assert evaluation.stratum_counts == (0,)


@pytest.mark.parametrize("field", ["key", "assay_class", "mutation_class"])
@pytest.mark.parametrize("value", ["", 1])
def test_outcome_rejects_invalid_required_text(field: str, value: object) -> None:
    with pytest.raises(ValueError, match="non-empty string"):
        _outcome(**{field: value})


@pytest.mark.parametrize("field", ["abstained", "applicable", "baseline_applicable"])
def test_outcome_rejects_non_boolean_state(field: str) -> None:
    with pytest.raises(ValueError, match="Boolean"):
        _outcome(**{field: 1})


@pytest.mark.parametrize(
    "changes",
    [
        {"expected_identity": None},
        {"expected_display_name": None},
        {"selected_identity": None},
        {"displayed_name": None},
    ],
)
def test_outcome_requires_joint_identity_and_display_name(changes: dict[str, object]) -> None:
    with pytest.raises(ValueError, match="jointly present"):
        _outcome(**changes)


def test_outcome_rejects_selection_on_abstention_and_selection_without_tier() -> None:
    with pytest.raises(ValueError, match="abstention cannot carry"):
        _outcome(abstained=True)
    with pytest.raises(ValueError, match="requires a fixed reconciliation tier"):
        _outcome(tier=None)


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"observations": ()}, "requires outcome observations"),
        ({"observations": (object(),)}, "OutcomeObservation"),
        ({"profile_sha256": "a" * 63}, "64 characters"),
        ({"profile_sha256": "z" * 64}, "lowercase hexadecimal"),
        ({"profile_sha256": 64}, "64 characters"),
        ({"free_parameter_count": True}, "non-negative integer"),
        ({"free_parameter_count": 1.0}, "non-negative integer"),
        ({"free_parameter_count": -1}, "non-negative integer"),
        ({"required_strata": ()}, "predeclared strata"),
        ({"required_strata": ("",)}, "unique non-empty"),
        ({"required_strata": (1,)}, "unique non-empty"),
        ({"required_strata": ("assay-a:duplication", "assay-a:duplication")}, "unique non-empty"),
    ],
)
def test_calculate_metrics_rejects_each_invalid_input_branch(changes: dict[str, object], message: str) -> None:
    arguments: dict[str, object] = {
        "observations": (_outcome(),),
        "profile_sha256": "a" * 64,
        "free_parameter_count": 0,
        "required_strata": ("assay-a:duplication",),
    }
    arguments.update(changes)
    with pytest.raises(ValueError, match=message):
        calculate_metrics(**arguments)  # type: ignore[arg-type]


def test_calculate_metrics_requires_mutated_truth_and_handles_no_applicable_denominator() -> None:
    control = _outcome(
        expected_identity=None,
        expected_display_name=None,
        selected_identity=None,
        displayed_name=None,
        tier=None,
        baseline_tier=None,
    )
    with pytest.raises(ValueError, match="mutated truth"):
        calculate_metrics(
            (control,), profile_sha256="a" * 64, free_parameter_count=0, required_strata=("assay-a:control",)
        )

    summary = calculate_metrics(
        (_outcome(baseline_applicable=False),),
        profile_sha256="a" * 64,
        free_parameter_count=0,
        required_strata=("assay-a:duplication",),
    )
    assert summary.metrics.abstention_fraction == 0


def test_tier_a_reachability_is_required_only_when_the_baseline_requires_it() -> None:
    without_baseline_a = calculate_metrics(
        (_outcome(tier="B", baseline_tier="B"),),
        profile_sha256="a" * 64,
        free_parameter_count=0,
        required_strata=("assay-a:duplication",),
    )
    lost_baseline_a = calculate_metrics(
        (_outcome(selected_identity=None, displayed_name=None, tier=None, abstained=True),),
        profile_sha256="a" * 64,
        free_parameter_count=0,
        required_strata=("assay-a:duplication",),
    )
    demoted_baseline_a = calculate_metrics(
        (_outcome(tier="B"),),
        profile_sha256="a" * 64,
        free_parameter_count=0,
        required_strata=("assay-a:duplication",),
    )

    assert without_baseline_a.metrics.tier_a_reachable
    assert not lost_baseline_a.metrics.tier_a_reachable
    assert not demoted_baseline_a.metrics.tier_a_reachable


def test_exact_recovery_counts_equality_and_empty_predeclared_stratum_as_zero() -> None:
    summary = calculate_metrics(
        (_outcome(),),
        profile_sha256="a" * 64,
        free_parameter_count=0,
        required_strata=("assay-a:duplication", "assay-b:deletion"),
    )

    assert summary.metrics.macro_exact_recovery == Fraction(1, 2)
    assert summary.stratum_counts == (1, 0)


def test_lexicographic_key_and_selection_reject_invalid_types_and_empty_family() -> None:
    protocol = decode_protocol(synthetic_protocol())
    with pytest.raises(ValueError, match="CandidateMetrics"):
        lexicographic_safety_key(object())  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="CalibrationProtocol"):
        select_candidate((), object())  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="CandidateEvaluation"):
        select_candidate((object(),), protocol)  # type: ignore[arg-type]
    assert select_candidate((), protocol) is None


def test_admissibility_accepts_every_constraint_at_its_inclusive_boundary() -> None:
    protocol = decode_protocol(synthetic_protocol())
    boundary = _evaluation(
        _metrics(
            abstention_fraction=protocol.maximum_abstention_fraction,
            free_parameter_count=protocol.maximum_free_parameters,
        ),
        stratum_counts=(protocol.minimum_stratum_count,),
        detection_lower_bound=Fraction(0),
        macro_exact_lower_bound=Fraction(0),
        holm_adjusted_p_value=Fraction(1, 20),
    )

    assert select_candidate((boundary,), protocol) == boundary


def test_free_parameter_contract_rejects_each_invalid_leaf() -> None:
    neutral: dict[str, object] = {
        "enabled": False,
        "minimum_record_count_margin": 0,
        "minimum_record_share": 0.0,
        "minimum_record_share_margin": 0.0,
        "xd_veto": "disabled",
        "abstain_on_inadmissible_advntr": False,
    }
    invalid = (
        {key: value for key, value in neutral.items() if key != "enabled"},
        {**neutral, "extra": 1},
        {**neutral, "enabled": 0},
        {**neutral, "abstain_on_inadmissible_advntr": 0},
        {**neutral, "minimum_record_count_margin": True},
        {**neutral, "minimum_record_count_margin": 1.0},
        {**neutral, "minimum_record_count_margin": -1},
        {**neutral, "minimum_record_share": True},
        {**neutral, "minimum_record_share": "0"},
        {**neutral, "minimum_record_share": -0.1},
        {**neutral, "minimum_record_share": 1.1},
        {**neutral, "minimum_record_share_margin": True},
        {**neutral, "minimum_record_share_margin": "0"},
        {**neutral, "minimum_record_share_margin": -0.1},
        {**neutral, "minimum_record_share_margin": 1.1},
        {**neutral, "xd_veto": "other"},
    )

    for component in invalid:
        with pytest.raises(ValueError):
            count_free_parameters(component)
