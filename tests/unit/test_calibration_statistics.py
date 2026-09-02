"""Deterministic group-level calibration statistics."""

from fractions import Fraction
from unittest.mock import patch

import pytest

from vntyper.scripts import calibration_statistics
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

    interval = paired_group_bootstrap(
        observations,
        required_strata=("assay-a:dup",),
        iterations=10_000,
        seed=295,
    )

    assert interval.resampling_unit_count == 2
    assert interval.iterations == 10_000
    assert interval != paired_group_bootstrap(
        tuple(
            PairedObservation(str(index), row.stratum, row.baseline, row.candidate)
            for index, row in enumerate(observations)
        ),
        required_strata=("assay-a:dup",),
        iterations=10_000,
        seed=295,
    )


def test_bootstrap_is_deterministic_and_identical_pairs_have_zero_bounds() -> None:
    rows = tuple(
        PairedObservation(f"group-{index}", "assay-a:dup", Fraction(index % 2), Fraction(index % 2))
        for index in range(5)
    )

    first = paired_group_bootstrap(rows, required_strata=("assay-a:dup",), iterations=10_000, seed=295)
    second = paired_group_bootstrap(rows, required_strata=("assay-a:dup",), iterations=10_000, seed=295)

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
        paired_group_bootstrap(
            crossing,
            required_strata=("assay-a:dup", "assay-b:dup"),
            iterations=10_000,
            seed=295,
        )
    with pytest.raises(ValueError):
        paired_group_bootstrap((), required_strata=("assay-a:dup",), iterations=10_000, seed=295)


def test_bootstrap_rejects_an_empty_cell_in_the_frozen_stratum_vector() -> None:
    rows = (PairedObservation("group", "assay-a:dup", Fraction(0), Fraction(1)),)

    with pytest.raises(ValueError, match="empty declared bootstrap strata.*assay-b:del"):
        paired_group_bootstrap(
            rows,
            required_strata=("assay-a:dup", "assay-b:del"),
            iterations=10_000,
            seed=295,
        )


def test_bootstrap_rejects_observations_outside_the_frozen_stratum_vector() -> None:
    rows = (PairedObservation("group", "assay-b:del", Fraction(0), Fraction(1)),)

    with pytest.raises(ValueError, match="undeclared bootstrap stratum.*assay-b:del"):
        paired_group_bootstrap(
            rows,
            required_strata=("assay-a:dup",),
            iterations=10_000,
            seed=295,
        )


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


def test_paired_observation_is_frozen() -> None:
    row = PairedObservation("group", "assay:dup", Fraction(0), Fraction(1))
    attribute = "group_id"

    with pytest.raises(AttributeError):
        setattr(row, attribute, "changed")


def test_statistical_result_value_objects_are_frozen() -> None:
    bootstrap = paired_group_bootstrap(
        (PairedObservation("group", "assay:dup", Fraction(0), Fraction(1)),),
        required_strata=("assay:dup",),
        iterations=10_000,
        seed=295,
    )
    binomial = clopper_pearson_interval(1, 2)
    curve = deterministic_curves(((Fraction(1), True), (Fraction(0), False)))[0]

    for value, attribute in ((bootstrap, "estimate"), (binomial, "lower"), (curve, "threshold")):
        with pytest.raises(AttributeError):
            setattr(value, attribute, Fraction(0))


@pytest.mark.parametrize(
    "arguments",
    [
        ("", "assay:dup", Fraction(0), Fraction(1)),
        (1, "assay:dup", Fraction(0), Fraction(1)),
        ("group", "", Fraction(0), Fraction(1)),
        ("group", 1, Fraction(0), Fraction(1)),
        ("group", "assay:dup", 0, Fraction(1)),
        ("group", "assay:dup", Fraction(-1, 10), Fraction(1)),
        ("group", "assay:dup", Fraction(0), 1),
        ("group", "assay:dup", Fraction(0), Fraction(11, 10)),
    ],
)
def test_paired_observation_rejects_every_invalid_identifier_and_rate(arguments: tuple[object, ...]) -> None:
    with pytest.raises(ValueError):
        PairedObservation(*arguments)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "changes",
    [
        {"iterations": 9_999},
        {"seed": True},
        {"seed": 1.0},
        {"seed": -1},
        {"observations": ()},
        {"observations": (object(),)},
    ],
)
def test_bootstrap_rejects_each_invalid_control(changes: dict[str, object]) -> None:
    arguments: dict[str, object] = {
        "observations": (PairedObservation("group", "assay:dup", Fraction(0), Fraction(1)),),
        "required_strata": ("assay:dup",),
        "iterations": 10_000,
        "seed": 295,
    }
    arguments.update(changes)

    with pytest.raises(ValueError):
        paired_group_bootstrap(**arguments)  # type: ignore[arg-type]


def test_bootstrap_reports_literal_estimate_and_percentile_bounds() -> None:
    rows = (
        PairedObservation("a", "assay-a:dup", Fraction(0), Fraction(1)),
        PairedObservation("b", "assay-a:dup", Fraction(1), Fraction(0)),
        PairedObservation("c", "assay-b:del", Fraction(0), Fraction(1)),
        PairedObservation("d", "assay-b:del", Fraction(0), Fraction(0)),
    )

    interval = paired_group_bootstrap(
        rows,
        required_strata=("assay-a:dup", "assay-b:del"),
        iterations=10_000,
        seed=295,
    )

    assert interval.estimate == Fraction(1, 4)
    assert interval.one_sided_lower == Fraction(-1, 2)
    assert interval.two_sided_lower == Fraction(-1, 2)
    assert interval.two_sided_upper == Fraction(1)
    assert interval.resampling_unit_count == 4


def test_bootstrap_reports_finite_corrected_one_sided_noninferiority_p_value() -> None:
    rows = tuple(
        PairedObservation(
            f"group-{index:02d}",
            "assay-a:dup",
            Fraction(index >= 13),
            Fraction(index < 13),
        )
        for index in range(20)
    )

    interval = paired_group_bootstrap(rows, required_strata=("assay-a:dup",), iterations=10_000, seed=295)

    assert interval.one_sided_lower == Fraction(0)
    assert interval.one_sided_noninferiority_p_value == Fraction(497, 10_001)
    assert interval.one_sided_noninferiority_p_value != Fraction(0)


def test_bootstrap_accepts_seed_zero_and_calls_the_three_literal_percentiles() -> None:
    rows = (PairedObservation("group", "assay:dup", Fraction(0), Fraction(1)),)

    with patch(
        "vntyper.scripts.calibration_statistics._percentile", wraps=lambda values, probability: values[0]
    ) as percentile:
        interval = paired_group_bootstrap(rows, required_strata=("assay:dup",), iterations=10_000, seed=0)

    assert interval.estimate == Fraction(1)
    assert [call.args[1] for call in percentile.call_args_list] == [
        Fraction(5, 100),
        Fraction(25, 1000),
        Fraction(975, 1000),
    ]


@pytest.mark.parametrize(
    ("events", "total", "confidence"),
    [
        (True, 1, Fraction(95, 100)),
        (1.0, 1, Fraction(95, 100)),
        (0, True, Fraction(95, 100)),
        (0, 1.0, Fraction(95, 100)),
        (0, 0, Fraction(95, 100)),
        (-1, 1, Fraction(95, 100)),
        (2, 1, Fraction(95, 100)),
        (0, 1, 0.95),
        (0, 1, Fraction(0)),
        (0, 1, Fraction(1)),
    ],
)
def test_clopper_pearson_rejects_each_invalid_contract_value(events: object, total: object, confidence: object) -> None:
    with pytest.raises(ValueError):
        clopper_pearson_interval(events, total, confidence=confidence)  # type: ignore[arg-type]


def test_clopper_pearson_pins_non_boundary_estimate_and_interval() -> None:
    interval = clopper_pearson_interval(5, 10)

    assert interval.estimate == Fraction(1, 2)
    assert interval.lower == Fraction(139778878661, 747136917818)
    assert interval.upper == Fraction(607358039157, 747136917818)
    assert (interval.events, interval.total) == (5, 10)

    one_event = clopper_pearson_interval(1, 10)
    assert one_event.lower == Fraction(1087267115, 429991434271)
    assert one_event.upper == Fraction(336279543117, 755656998139)


def test_clopper_pearson_confidence_is_keyword_only() -> None:
    with pytest.raises(TypeError):
        clopper_pearson_interval(1, 2, Fraction(9, 10))  # type: ignore[misc]


@pytest.mark.parametrize(
    "scores",
    [
        (),
        ((0, True),),
        ((Fraction(0), 1),),
        ((Fraction(1), True), (Fraction(0), True)),
        ((Fraction(1), False), (Fraction(0), False)),
    ],
)
def test_curves_reject_invalid_rows_and_single_class_inputs(scores: object) -> None:
    with pytest.raises(ValueError):
        deterministic_curves(scores)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "scores",
    [
        object(),
        ((0, True), (Fraction(0), False)),
        ((Fraction(1), 1), (Fraction(0), False)),
    ],
)
def test_curves_reject_non_sequences_and_each_invalid_row_field(scores: object) -> None:
    with pytest.raises(ValueError, match="calibration curves|calibration curve rows"):
        deterministic_curves(scores)  # type: ignore[arg-type]


def test_curve_rows_pin_every_confusion_count_and_rate() -> None:
    scores = (
        (Fraction(3, 4), True),
        (Fraction(1, 2), False),
        (Fraction(1, 4), True),
    )

    points = deterministic_curves(scores)

    assert points == (
        type(points[0])(Fraction(3, 4), 1, 0, 1, 1, Fraction(1, 2), Fraction(0), Fraction(1), Fraction(1, 2)),
        type(points[0])(Fraction(1, 2), 1, 1, 0, 1, Fraction(1, 2), Fraction(1), Fraction(1, 2), Fraction(1, 2)),
        type(points[0])(Fraction(1, 4), 2, 1, 0, 0, Fraction(1), Fraction(1), Fraction(2, 3), Fraction(1)),
    )


@pytest.mark.parametrize(
    "grid",
    [
        {},
        {"axis": ()},
        {"axis": "not-an-axis"},
        {"axis": 1},
    ],
)
def test_joint_surface_rejects_empty_or_non_sequence_axes(grid: object) -> None:
    with pytest.raises(ValueError):
        joint_surface(grid)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "p_values",
    [
        {},
        {"": Fraction(1, 2)},
        {1: Fraction(1, 2)},
        {"a": 0.5},
        {"a": Fraction(-1, 10)},
        {"a": Fraction(11, 10)},
    ],
)
def test_holm_rejects_invalid_names_and_values(p_values: object) -> None:
    with pytest.raises(ValueError):
        holm_adjust(p_values)  # type: ignore[arg-type]


def test_holm_caps_at_one_preserves_input_order_and_breaks_ties_by_name() -> None:
    adjusted = holm_adjust({"z": Fraction(1, 2), "b": Fraction(1, 100), "a": Fraction(1, 100)})

    assert tuple(adjusted) == ("z", "b", "a")
    assert adjusted == {"z": Fraction(1, 2), "b": Fraction(3, 100), "a": Fraction(3, 100)}
    assert holm_adjust({"a": Fraction(1, 2), "b": Fraction(3, 4)}) == {
        "a": Fraction(1),
        "b": Fraction(1),
    }
    assert holm_adjust({"a": Fraction(3, 4), "b": Fraction(3, 4)}) == {
        "a": Fraction(1),
        "b": Fraction(1),
    }


def test_percentile_uses_the_declared_zero_based_sample_index() -> None:
    assert calibration_statistics._percentile((Fraction(0), Fraction(1), Fraction(2)), Fraction(1, 2)) == Fraction(1)
