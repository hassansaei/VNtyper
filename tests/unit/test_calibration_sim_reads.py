"""Known-origin paired reads generated independently of production alignment."""

import random
from importlib import import_module

import pytest
from calibration_sim.haplotypes import build_haplotype

pytestmark = pytest.mark.unit


def _haplotypes():
    return tuple(
        build_haplotype(
            repeat_units=(unit,) * count,
            repeat_unit_bp=4,
            left_flank="GGGCCC",
            right_flank="TTTAAA",
            maximum_haplotype_bp=200,
        )
        for unit, count in (("ACGT", 3), ("TGCA", 8))
    )


def _generate(**changes):
    m = import_module("calibration_sim.reads")
    arguments = {
        "haplotypes": _haplotypes(),
        "pair_count": 20,
        "read_length": 5,
        "fragment_length": 12,
        "seed": 7,
        "substitution_rate": 0.0,
        "quality_score": 35,
        "allele_copy_weights": (1.0, 1.0),
        "maximum_generated_bases": 1000,
    }
    arguments.update(changes)
    return m.generate_read_pairs(**arguments)


def test_pairs_reconstruct_both_ends_of_the_original_fragment_in_forward_fastq_orientation():
    haplotypes = _haplotypes()
    pairs = list(_generate())
    assert len(pairs) == 20
    assert len({pair.name for pair in pairs}) == 20
    for index, pair in enumerate(pairs):
        source = haplotypes[pair.allele_index].sequence
        assert pair.name == f"synthetic-read-{index:012d}"
        assert pair.fragment_end - pair.fragment_start == 12
        assert 0 <= pair.fragment_start < pair.fragment_end <= len(source)
        assert pair.read1 == source[pair.fragment_start : pair.fragment_start + 5]
        forward_read2 = source[pair.fragment_end - 5 : pair.fragment_end]
        # Literal oracle, independent of the generator's complement function.
        complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
        assert pair.read2 == "".join(complement[base] for base in reversed(forward_read2))
        assert pair.qualities == "D" * 5
        assert pair.read1_substitutions == pair.read2_substitutions == 0


def test_seed_reproducibility_does_not_mutate_process_global_random_state():
    before = random.getstate()
    first = list(_generate())
    assert list(_generate()) == first
    assert list(_generate(seed=8)) != first
    assert random.getstate() == before


def test_copy_weight_dropout_and_zero_coverage_are_explicit():
    assert list(_generate(pair_count=0, maximum_generated_bases=0)) == []
    assert {pair.allele_index for pair in _generate(allele_copy_weights=(0.0, 1.0))} == {1}
    assert {pair.allele_index for pair in _generate(allele_copy_weights=(1.0, 0.0))} == {0}


def test_equal_copy_weights_sample_uniform_fragment_starts_not_equal_read_counts_per_allele():
    # 13 versus33 possible starts; equal copy number must preserve that2.538... ratio.
    pairs = list(_generate(pair_count=10000, maximum_generated_bases=100000))
    fraction_longer = sum(pair.allele_index == 1 for pair in pairs) / len(pairs)
    assert abs(fraction_longer - 33 / 46) < 0.015


def test_copy_weight_scaling_leaves_seeded_reads_identical_and_avoids_weight_overflow():
    assert list(_generate(allele_copy_weights=(1e308, 1e308))) == list(_generate())
    assert list(_generate(allele_copy_weights=(1.0, 2.0))) == list(_generate(allele_copy_weights=(5.0, 10.0)))


def test_substitution_rate_one_changes_every_base_but_preserves_fragment_truth_and_quality():
    haplotypes = _haplotypes()
    pairs = list(_generate(substitution_rate=1.0))
    complement = str.maketrans("ACGT", "TGCA")
    for pair in pairs:
        source = haplotypes[pair.allele_index].sequence
        first = source[pair.fragment_start : pair.fragment_start + 5]
        second = source[pair.fragment_end - 5 : pair.fragment_end].translate(complement)[::-1]
        assert all(observed != original for observed, original in zip(pair.read1, first, strict=True))
        assert all(observed != original for observed, original in zip(pair.read2, second, strict=True))
        assert pair.read1_substitutions == pair.read2_substitutions == 5
        assert pair.qualities == "D" * 5


def test_overlapping_mates_and_single_possible_fragment_are_valid():
    haplotype = _haplotypes()[0]
    pairs = list(_generate(haplotypes=(haplotype, haplotype), read_length=20, fragment_length=24))
    assert all(pair.fragment_start == 0 and pair.fragment_end == 24 for pair in pairs)
    assert all(len(pair.read1) == len(pair.read2) == 20 for pair in pairs)


@pytest.mark.parametrize(
    "changes",
    [
        {"pair_count": True},
        {"pair_count": -1},
        {"read_length": 0},
        {"read_length": 13},
        {"fragment_length": 25},
        {"fragment_length": False},
        {"seed": -1},
        {"seed": True},
        {"substitution_rate": float("nan")},
        {"substitution_rate": -0.01},
        {"substitution_rate": 1.01},
        {"substitution_rate": True},
        {"substitution_rate": 10**400},
        {"quality_score": 94},
        {"quality_score": -1},
        {"quality_score": True},
        {"allele_copy_weights": (0.0, 0.0)},
        {"allele_copy_weights": (float("inf"), 1.0)},
        {"allele_copy_weights": (-1.0, 1.0)},
        {"allele_copy_weights": (True, 1.0)},
        {"allele_copy_weights": [1.0, 1.0]},
        {"maximum_generated_bases": 199},
        {"maximum_generated_bases": True},
        {"haplotypes": ()},
        {"haplotypes": (object(), object())},
    ],
)
def test_invalid_simulation_contract_is_rejected_before_iteration(changes):
    with pytest.raises(ValueError):
        _generate(**changes)


def test_explicit_base_budget_is_inclusive_and_output_is_streamed():
    pairs = _generate(maximum_generated_bases=200)
    assert iter(pairs) is pairs
    assert next(pairs).name == "synthetic-read-000000000000"
    assert len(list(pairs)) == 19


def test_partial_error_rate_tracks_actual_substitutions_against_independent_fragment_bases():
    haplotypes = _haplotypes()
    pairs = list(_generate(pair_count=2000, substitution_rate=0.1, maximum_generated_bases=20000))
    changes = 0
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    for pair in pairs:
        source = haplotypes[pair.allele_index].sequence
        first = source[pair.fragment_start : pair.fragment_start + 5]
        second = "".join(complement[base] for base in reversed(source[pair.fragment_end - 5 : pair.fragment_end]))
        first_changes = sum(a != b for a, b in zip(pair.read1, first, strict=True))
        second_changes = sum(a != b for a, b in zip(pair.read2, second, strict=True))
        assert pair.read1_substitutions == first_changes
        assert pair.read2_substitutions == second_changes
        changes += first_changes + second_changes
    assert abs(changes / 20000 - 0.1) < 0.01


def test_error_rate_changes_preserve_paired_fragment_origins():
    clean = list(_generate(substitution_rate=0.0))
    noisy = list(_generate(substitution_rate=0.2))
    assert [(row.allele_index, row.fragment_start, row.fragment_end) for row in clean] == [
        (row.allele_index, row.fragment_start, row.fragment_end) for row in noisy
    ]


def test_increasing_error_rate_adds_errors_without_moving_existing_substitutions():
    clean = list(_generate(substitution_rate=0.0))
    low = list(_generate(substitution_rate=0.1))
    high = list(_generate(substitution_rate=0.4))
    for original, fewer, more in zip(clean, low, high, strict=True):
        for original_read, fewer_read, more_read in (
            (original.read1, fewer.read1, more.read1),
            (original.read2, fewer.read2, more.read2),
        ):
            for reference, first, second in zip(original_read, fewer_read, more_read, strict=True):
                if first != reference:
                    assert second == first
