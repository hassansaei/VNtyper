"""Deterministic group-level calibration statistics."""

from fractions import Fraction

import pytest

from vntyper.scripts.calibration_statistics import (
    PairedObservation,
    clopper_pearson_interval,
    deterministic_curves,
    holm_adjust,
    joint_surface,
    paired_group_bootstrap,
)

pytestmark = pytest.mark.unit


def test_group_bootstrap_resamples_complete_groups_not_rows() -> None:
    observations = (
        PairedObservation("family-a", "assay-a:dup", Fraction(0), Fraction(1)),
        PairedObservation("family-a", "assay-a:dup", Fraction(0), Fraction(1)),
        PairedObservation("family-b", "assay-a:dup", Fraction(1), Fraction(0)),
    )

    interval = paired_group_bootstrap(observations, iterations=10_000, seed=295)

    assert interval.resampling_unit_count == 2
    assert interval.iterations == 10_000
    assert interval != paired_group_bootstrap(
        tuple(
            PairedObservation(str(index), row.stratum, row.baseline, row.candidate)
            for index, row in enumerate(observations)
        ),
        iterations=10_000,
        seed=295,
    )


def test_bootstrap_is_deterministic_and_identical_pairs_have_zero_bounds() -> None:
    rows = tuple(
        PairedObservation(f"group-{index}", "assay-a:dup", Fraction(index % 2), Fraction(index % 2))
        for index in range(5)
    )

    first = paired_group_bootstrap(rows, iterations=10_000, seed=295)
    second = paired_group_bootstrap(rows, iterations=10_000, seed=295)

    assert first == second
    assert (first.one_sided_lower, first.two_sided_lower, first.two_sided_upper) == (
        Fraction(0),
        Fraction(0),
        Fraction(0),
    )


def test_group_spanning_strata_and_empty_input_fail_closed() -> None:
    crossing = (
        PairedObservation("shared", "assay-a:dup", Fraction(0), Fraction(1)),
        PairedObservation("shared", "assay-b:dup", Fraction(0), Fraction(1)),
    )

    with pytest.raises(ValueError, match="strata"):
        paired_group_bootstrap(crossing, iterations=10_000, seed=295)
    with pytest.raises(ValueError):
        paired_group_bootstrap((), iterations=10_000, seed=295)


def test_clopper_pearson_zero_event_interval_is_exact_boundary_case() -> None:
    interval = clopper_pearson_interval(0, 50)

    assert interval.lower == Fraction(0)
    assert 0.05 < float(interval.upper) < 0.08
    assert clopper_pearson_interval(50, 50).upper == Fraction(1)


def test_roc_and_precision_recall_tables_are_deterministic() -> None:
    scores = (
        (Fraction(9, 10), True),
        (Fraction(4, 5), False),
        (Fraction(1, 2), True),
        (Fraction(1, 10), False),
    )

    points = deterministic_curves(scores)

    assert tuple(point.threshold for point in points) == tuple(sorted({score for score, _ in scores}, reverse=True))
    assert points[0].true_positives == 1 and points[0].false_positives == 0
    assert points[-1].recall == Fraction(1)
    assert points == deterministic_curves(tuple(reversed(scores)))


def test_joint_surface_enumerates_the_complete_declared_cross_product() -> None:
    surface = joint_surface({"count_margin": (1, 2), "share": (Fraction(1, 2), Fraction(3, 4))})

    assert surface == (
        (("count_margin", 1), ("share", Fraction(1, 2))),
        (("count_margin", 1), ("share", Fraction(3, 4))),
        (("count_margin", 2), ("share", Fraction(1, 2))),
        (("count_margin", 2), ("share", Fraction(3, 4))),
    )


def test_holm_adjustment_is_monotone_in_sorted_p_values() -> None:
    adjusted = holm_adjust({"a": Fraction(1, 100), "b": Fraction(3, 100), "c": Fraction(1, 5)})

    assert adjusted == {"a": Fraction(3, 100), "b": Fraction(3, 50), "c": Fraction(1, 5)}
