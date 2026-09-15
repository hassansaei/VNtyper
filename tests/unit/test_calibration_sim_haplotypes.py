"""Independent sequence construction for known-truth calibration simulations."""

from dataclasses import replace
from importlib import import_module

import pytest

pytestmark = pytest.mark.unit


def _module():
    return import_module("scripts.calibration_sim.haplotypes")


def _build(units=("ACGT", "TGCA"), **changes):
    arguments = {
        "repeat_units": units,
        "repeat_unit_bp": 4,
        "left_flank": "GG",
        "right_flank": "TT",
        "maximum_haplotype_bp": 2000,
    }
    arguments.update(changes)
    return _module().build_haplotype(**arguments)


def test_haplotype_sequence_and_repeat_truth_are_constructed_independently():
    result = _build()
    assert result.sequence == "GGACGTTGCATT"
    assert result.repeat_count == 2
    assert result.repeat_start == 2
    assert result.repeat_end == 10
    assert result.reference_repeat_bp == 8
    assert result.edits == ()
    assert len(result.sha256) == 64
    assert _build().sha256 == result.sha256
    with pytest.raises(AttributeError):
        result.sequence = "tampered"


def test_insertion_deletion_and_substitution_apply_at_original_repeat_coordinates():
    m = _module()
    result = _build(edits=(m.RepeatEdit(0, 1, "", "C"), m.RepeatEdit(1, 0, "T", ""), m.RepeatEdit(1, 2, "C", "A")))
    assert result.sequence == "GGACCGTGAATT"
    assert result.repeat_count == 2
    assert result.reference_repeat_bp == 8
    assert result.repeat_end == 10
    assert result.sha256 != _build().sha256


def test_frameshift_does_not_redefine_repeat_count_by_dividing_mutated_base_length():
    m = _module()
    inserted = _build(edits=(m.RepeatEdit(0, 1, "", "C"),))
    deleted = _build(edits=(m.RepeatEdit(1, 1, "G", ""),))
    assert inserted.sequence == "GGACCGTTGCATT"
    assert inserted.repeat_end == 11
    assert deleted.sequence == "GGACGTTCATT"
    assert deleted.repeat_end == 9
    assert inserted.repeat_count == deleted.repeat_count == 2


def test_edit_at_repeat_end_precedes_next_unit_and_never_changes_flanks():
    m = _module()
    result = _build(edits=(m.RepeatEdit(0, 4, "", "G"),))
    assert result.sequence == "GGACGTGTGCATT"
    assert result.sequence[: result.repeat_start] == "GG"
    assert result.sequence[result.repeat_end :] == "TT"


def test_diploid_total_preserves_balanced_asymmetric_and_long_alleles():
    m = _module()
    for first, second in ((40, 40), (40, 90), (40, 180)):
        left = _build(units=("ACGT",) * first)
        right = _build(units=("TGCA",) * second)
        truth = m.diploid_truth(left, right)
        assert truth.allele_repeat_counts == (first, second)
        assert truth.total_repeat_count == first + second
        assert truth.haplotype_sha256 == (left.sha256, right.sha256)
        swapped = m.diploid_truth(right, left)
        assert swapped.total_repeat_count == truth.total_repeat_count
        assert swapped.allele_repeat_counts == (second, first)


@pytest.mark.parametrize(
    "changes",
    [
        {"repeat_units": ()},
        {"repeat_units": ["ACGT"]},
        {"repeat_units": ("ACGT", "AAA")},
        {"repeat_units": ("acgt",)},
        {"repeat_units": ("ACGN",)},
        {"repeat_unit_bp": True},
        {"repeat_unit_bp": 0},
        {"maximum_haplotype_bp": True},
        {"maximum_haplotype_bp": 11},
        {"left_flank": "A N"},
        {"right_flank": ""},
        {"edits": []},
    ],
)
def test_invalid_or_over_budget_sequence_contracts_fail(changes):
    with pytest.raises(ValueError):
        _build(**changes)


@pytest.mark.parametrize(
    "edit",
    [
        (True, 0, "", "A"),
        (2, 0, "", "A"),
        (0, -1, "", "A"),
        (0, 5, "", "A"),
        (0, 3, "TT", ""),
        (0, 0, "T", ""),
        (0, 0, "A", "A"),
        (0, 0, "", ""),
        (0, 0, "", "N"),
    ],
)
def test_invalid_or_mismatched_edits_fail(edit):
    with pytest.raises(ValueError):
        _build(edits=(_module().RepeatEdit(*edit),))


def test_overlapping_unsorted_and_same_boundary_edits_are_rejected():
    m = _module()
    cases = (
        (m.RepeatEdit(0, 0, "AC", ""), m.RepeatEdit(0, 1, "", "T")),
        (m.RepeatEdit(1, 0, "", "A"), m.RepeatEdit(0, 0, "", "T")),
        (m.RepeatEdit(0, 4, "", "A"), m.RepeatEdit(1, 0, "", "T")),
        (object(),),
    )
    for edits in cases:
        with pytest.raises(ValueError):
            _build(edits=edits)


def test_insertions_count_against_the_explicit_output_budget():
    m = _module()
    with pytest.raises(ValueError, match="budget"):
        _build(edits=(m.RepeatEdit(0, 0, "", "AAA"),), maximum_haplotype_bp=14)
    assert len(_build(edits=(m.RepeatEdit(0, 0, "", "AAA"),), maximum_haplotype_bp=15).sequence) == 15


def test_diploid_truth_rejects_tampered_sequence_count_and_digest():
    m = _module()
    first = _build()
    for changed in (
        replace(first, sequence="ACGT"),
        replace(first, repeat_count=3),
        replace(first, sha256="0" * 64),
        object(),
    ):
        with pytest.raises(ValueError):
            m.diploid_truth(changed, first)


def test_whole_unit_gain_or_loss_requires_a_different_source_unit_construction():
    m = _module()
    for edits in (
        (m.RepeatEdit(0, 0, "ACGT", ""),),
        (m.RepeatEdit(0, 0, "", "ACGT"),),
        (m.RepeatEdit(0, 0, "AC", ""), m.RepeatEdit(0, 2, "GT", "")),
    ):
        with pytest.raises(ValueError, match="unit count"):
            _build(edits=edits)


def test_diploid_rejects_incompatible_unit_geometry_and_boolean_count():
    m = _module()
    first = _build(units=("ACGT",))
    second = _build(units=("ACGTACGT",), repeat_unit_bp=8)
    with pytest.raises(ValueError, match="unit length"):
        m.diploid_truth(first, second)
    with pytest.raises(ValueError, match="integer"):
        m.diploid_truth(replace(first, repeat_count=True), first)


def test_sixty_base_units_with_long_asymmetric_truth_have_no_implicit_upper_cap():
    m = _module()
    first = _build(units=("ACGT" * 15,) * 40, repeat_unit_bp=60, maximum_haplotype_bp=20000)
    second = _build(units=("TGCA" * 15,) * 180, repeat_unit_bp=60, maximum_haplotype_bp=20000)
    assert len(second.sequence) == 10804
    assert second.reference_repeat_bp == 10800
    assert m.diploid_truth(first, second).total_repeat_count == 220


def test_distributed_small_edits_cannot_silently_create_or_remove_a_whole_unit():
    m = _module()
    for edits in (
        (m.RepeatEdit(0, 1, "", "AA"), m.RepeatEdit(1, 1, "", "AA")),
        (m.RepeatEdit(0, 1, "AA", ""), m.RepeatEdit(1, 1, "AA", "")),
    ):
        with pytest.raises(ValueError, match="unit count"):
            _build(units=("AAAA", "AAAA"), edits=edits)
