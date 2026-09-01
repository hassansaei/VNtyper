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
