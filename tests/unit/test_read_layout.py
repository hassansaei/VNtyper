"""Read-layout decisions must account for every FASTQ conversion output."""

from __future__ import annotations

import pytest

from vntyper.scripts.read_layout import classify_layout, route_fastqs

pytestmark = pytest.mark.unit

PATHS = {
    "r1": "r1.gz",
    "r2": "r2.gz",
    "other": "other.gz",
    "single": "single.gz",
}


def test_r1_and_r2_populated_with_equal_counts_is_paired() -> None:
    assert classify_layout(r1=900, r2=900, other=0, single=0) == "paired"


def test_everything_in_other_is_single_because_samtools_puts_unpaired_reads_there() -> None:
    """``samtools fastq -0`` receives reads whose READ1 and READ2 flags are both unset."""
    assert classify_layout(r1=0, r2=0, other=1800, single=0) == "single"


def test_everything_in_the_singleton_output_is_also_a_single_read_set() -> None:
    assert classify_layout(r1=0, r2=0, other=0, single=17) == "single"


def test_paired_and_unpaired_outputs_together_are_mixed() -> None:
    assert classify_layout(r1=900, r2=900, other=4200, single=0) == "mixed"


def test_two_distinct_unpaired_outputs_are_mixed_not_implicitly_combined() -> None:
    assert classify_layout(r1=0, r2=0, other=10, single=5) == "mixed"


@pytest.mark.parametrize(
    "counts",
    [
        {"r1": 900, "r2": 880, "other": 0, "single": 0},
        {"r1": 10, "r2": 0, "other": 0, "single": 0},
        {"r1": 0, "r2": 10, "other": 0, "single": 0},
    ],
)
def test_invalid_mate_parity_routes_nothing_and_strands_every_nonempty_path(counts: dict[str, int]) -> None:
    assert classify_layout(**counts) == "invalid"
    consumed, stranded = route_fastqs("invalid", PATHS, counts)
    assert consumed == ()
    assert stranded == tuple(PATHS[key] for key in ("r1", "r2") if counts[key] > 0)


def test_nothing_at_all_is_empty_not_single() -> None:
    assert classify_layout(r1=0, r2=0, other=0, single=0) == "empty"


def test_negative_counts_are_rejected_instead_of_becoming_an_empty_layout() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        classify_layout(r1=-1, r2=0, other=0, single=0)


def test_paired_layout_routes_both_mates_in_order() -> None:
    consumed, stranded = route_fastqs("paired", PATHS, {"r1": 900, "r2": 900, "other": 0, "single": 0})

    assert consumed == ("r1.gz", "r2.gz")
    assert stranded == ()


@pytest.mark.parametrize(
    ("counts", "expected"),
    [
        ({"r1": 0, "r2": 0, "other": 12, "single": 0}, ("other.gz",)),
        ({"r1": 0, "r2": 0, "other": 0, "single": 12}, ("single.gz",)),
    ],
)
def test_single_layout_routes_the_one_populated_unpaired_file(
    counts: dict[str, int], expected: tuple[str, ...]
) -> None:
    consumed, stranded = route_fastqs("single", PATHS, counts)

    assert consumed == expected
    assert stranded == ()


def test_a_stranded_non_empty_file_is_reported_never_dropped() -> None:
    _, stranded = route_fastqs("paired", PATHS, {"r1": 900, "r2": 900, "other": 4200, "single": 0})

    assert stranded == ("other.gz",)


@pytest.mark.parametrize(
    ("counts", "verdict", "selected"),
    [
        ({"r1": 14_690, "r2": 14_690, "other": 0, "single": 1}, "mixed", ("r1.gz", "r2.gz", "single.gz")),
        ({"r1": 3_474, "r2": 3_474, "other": 0, "single": 93}, "mixed", ("r1.gz", "r2.gz", "single.gz")),
        (
            {"r1": 9, "r2": 9, "other": 4, "single": 2},
            "mixed",
            ("r1.gz", "r2.gz", "other.gz", "single.gz"),
        ),
        ({"r1": 0, "r2": 0, "other": 4, "single": 2}, "mixed", ("other.gz", "single.gz")),
    ],
)
def test_lossless_layout_routes_every_nonempty_file_once(
    counts: dict[str, int], verdict: str, selected: tuple[str, ...]
) -> None:
    assert classify_layout(**counts) == verdict
    consumed, stranded = route_fastqs(verdict, PATHS, counts)
    assert consumed == selected
    assert stranded == ()


def test_duplicate_selected_paths_are_rejected_instead_of_double_counted() -> None:
    aliased_paths = {**PATHS, "other": "unpaired.gz", "single": "unpaired.gz"}
    counts = {"r1": 0, "r2": 0, "other": 4, "single": 2}

    with pytest.raises(ValueError, match="duplicate"):
        route_fastqs("mixed", aliased_paths, counts)


def test_empty_layout_has_no_consumed_or_stranded_files() -> None:
    assert route_fastqs("empty", PATHS, {"r1": 0, "r2": 0, "other": 0, "single": 0}) == ((), ())


def test_an_unknown_layout_is_rejected_instead_of_dropping_every_file() -> None:
    with pytest.raises(ValueError, match="Unknown FASTQ layout"):
        route_fastqs("typo", PATHS, {"r1": 1, "r2": 1, "other": 0, "single": 0})
