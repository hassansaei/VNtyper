"""Public reference-derived MUC1 length measurement definition."""

from __future__ import annotations

import json
from importlib.resources import files

import pytest

from vntyper.scripts.length_annotation import Interval, decode_length_annotation
from vntyper.scripts.length_reference_annotation import (
    grch38_length_annotation,
    verify_grch38_length_reference,
)

pytestmark = pytest.mark.unit


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _target_reference_sequence() -> str:
    pre = (
        "AAGGAGACTTCGGCTACCCAGAGAAGTTCAGTGCCCAGCTCTACTGAGAAGAATGCTGTG"
        "AGTATGACCAGCAGCGTACTCTCCAGCCACAGCCCCGGTTCAGGCTCCTCCACCACTCAG"
        "GGACAGGATGTCACTCTGGCCCCGGCCACGGAACCAGCTTCAGGTTCAGCTGCCACCTGG"
        "GGACAGGATGTCACCTCGGTCCCAGTCACCAGGCCAGCCCTGGGCTCCACCACCCCGCCA"
        "GCCCACGATGTCACCTCAGCCCCGGACAACAAGCCAGCCCCGGGCTCCACCGCCCCCCCA"
    )
    after = (
        "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCGGGCCCCGGGCTCCACCCCGGCCCCG"
        "GGCTCCACCGCCCCCCCAGCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCG"
        "GGCTCCACCGCCCCCCCAGCCCATGGTGTCACCTCGGCCCCGGACAACAGGCCCGCCTTG"
        "GGCTCCACCGCCCCTCCAGTCCACAATGTCACCTCGGCCTCAGGCTCTGCATCAGGCTCA"
    )
    core = "C" * 3213
    return _reverse_complement(after) + core + _reverse_complement(pre)


def test_grch38_definition_separates_target_array_and_external_unique_flanks() -> None:
    annotation = grch38_length_annotation()

    assert annotation.core == (Interval(155188726, 155191939),)
    assert annotation.invariant == (Interval(155188486, 155188726), Interval(155191939, 155192239))
    assert annotation.array == Interval(155188529, 155192010)
    assert annotation.left_flank == Interval(155188296, 155188486)
    assert annotation.right_flank == Interval(155192239, 155192429)
    assert annotation.reference_core_repeat_count is None
    assert annotation.reference_invariant_repeat_count == 9
    assert (annotation.array_only_bp, annotation.target_only_bp) == (0, 272)
    assert not annotation.physical_a_compatible
    assert not annotation.physical_f_compatible


def test_reference_verifier_checks_both_exact_terminal_sequence_witnesses() -> None:
    annotation = grch38_length_annotation()
    sequence = _target_reference_sequence()

    verify_grch38_length_reference(annotation, sequence)

    for position in (0, len(sequence) - 1):
        changed = sequence[:position] + ("A" if sequence[position] != "A" else "C") + sequence[position + 1 :]
        with pytest.raises(ValueError, match="terminal sequence witness"):
            verify_grch38_length_reference(annotation, changed)


def test_reference_verifier_rejects_wrong_geometry_length_type_and_annotation() -> None:
    annotation = grch38_length_annotation()
    sequence = _target_reference_sequence()

    with pytest.raises(ValueError, match="target sequence length"):
        verify_grch38_length_reference(annotation, sequence[:-1])
    with pytest.raises(ValueError, match="target sequence must be DNA text"):
        verify_grch38_length_reference(annotation, 42)  # type: ignore[arg-type]
    synthetic = decode_length_annotation(
        {
            "schema_version": "length-annotation-v1",
            "assembly": "synthetic",
            "contig": "chr1",
            "accepted_contig_aliases": [],
            "reference_fasta_sha256": "a" * 64,
            "coordinate_system": "zero-based-half-open",
            "boundary_definition": "complete-core-plus-invariant-units-v1",
            "repeat_unit_bp": 1,
            "regions": {"CORE": None, "INVARIANT": None, "ARRAY": None, "LEFT_FLANK": None, "RIGHT_FLANK": None},
            "array_boundary_geometry": {"array_only_bp": None, "target_only_bp": None},
            "target_boundary_conversion_sha256": None,
            "physical_hypothesis_compatibility": {"physical_A": False, "physical_F": False},
            "annotation_provenance": "synthetic",
            "annotation_version": "synthetic",
        }
    )
    with pytest.raises(ValueError, match="approved GRCh38 annotation"):
        verify_grch38_length_reference(synthetic, sequence)


def test_packaged_annotation_matches_the_reviewed_builder() -> None:
    path = files("vntyper").joinpath("data/length/grch38-length-annotation-v1.json")
    packaged = json.loads(path.read_text(encoding="utf-8"))

    assert decode_length_annotation(packaged) == grch38_length_annotation()
