"""Logical read identity uses content, orientation, mate and multiplicity."""

from dataclasses import replace
from importlib import import_module

import pytest

pytestmark = pytest.mark.unit


def record(**changes):
    m = import_module("vntyper.scripts.calibration_read_identity")
    values = {
        "name": "synthetic-read",
        "sequence": "ACGN",
        "qualities": (10, 20, 30, 40),
        "mate": 1,
        "flags": 65,
        "mapping_quality": 60,
        "contig": "synthetic-contig",
        "position_zero_based": 10,
        "cigar": ((0, 4),),
        "mate_contig": "synthetic-contig",
        "mate_position_zero_based": 30,
        "template_length": 24,
    }
    values.update(changes)
    return m.PrimaryReadRecord(**values)


def test_forward_oriented_read_content_matches_reverse_alignment_representation():
    m = import_module("vntyper.scripts.calibration_read_identity")
    forward = m.read_identity_tokens(record())
    reverse = m.read_identity_tokens(record(sequence="NCGT", qualities=(40, 30, 20, 10), flags=81))
    assert forward.named_sequence == reverse.named_sequence
    assert forward.unnamed_sequence == reverse.unnamed_sequence
    assert forward.alignment != reverse.alignment
    assert forward.reliable and reverse.reliable


def test_read_renaming_only_changes_named_identity_and_alignment_identity():
    m = import_module("vntyper.scripts.calibration_read_identity")
    first = m.read_identity_tokens(record())
    second = m.read_identity_tokens(record(name="renamed"))
    assert first.named_sequence != second.named_sequence
    assert first.alignment != second.alignment
    assert first.unnamed_sequence == second.unnamed_sequence


@pytest.mark.parametrize(
    "changes",
    [
        {"mate": 2, "flags": 129},
        {"qualities": (10, 20, 30, 41)},
        {"sequence": "ACGT"},
    ],
)
def test_mate_quality_and_sequence_remain_part_of_unnamed_identity(changes):
    m = import_module("vntyper.scripts.calibration_read_identity")
    assert (
        m.read_identity_tokens(record()).unnamed_sequence != m.read_identity_tokens(record(**changes)).unnamed_sequence
    )


@pytest.mark.parametrize(
    "changes",
    [
        {"position_zero_based": 11},
        {"cigar": ((0, 3), (1, 1))},
        {"flags": 1089},
        {"mate_position_zero_based": 31},
        {"template_length": 25},
    ],
)
def test_alignment_fingerprint_keeps_geometry_and_flags(changes):
    m = import_module("vntyper.scripts.calibration_read_identity")
    first = m.read_identity_tokens(record())
    second = m.read_identity_tokens(record(**changes))
    assert first.named_sequence == second.named_sequence
    assert first.alignment != second.alignment


@pytest.mark.parametrize("flags", [321, 2113])
def test_nonprimary_records_are_excluded_from_primary_read_identity(flags):
    m = import_module("vntyper.scripts.calibration_read_identity")
    assert m.read_identity_tokens(record(flags=flags)) is None


def test_hard_clips_and_missing_sequence_do_not_claim_reconstructable_source_reads():
    m = import_module("vntyper.scripts.calibration_read_identity")
    clipped = m.read_identity_tokens(record(cigar=((5, 3), (0, 4))))
    assert not clipped.reliable
    assert clipped.reasons == ("hard_clipped_sequence",)
    missing = m.read_identity_tokens(record(sequence=None, qualities=None))
    assert not missing.reliable
    assert "missing_sequence" in missing.reasons


@pytest.mark.parametrize(
    "changes",
    [
        {"mate": True},
        {"mate": 3},
        {"flags": True},
        {"flags": -1},
        {"name": ""},
        {"sequence": "ACZQ"},
        {"sequence": ""},
        {"qualities": (10, 20)},
        {"qualities": (10, 20, 30, -1)},
        {"qualities": (10, 20, 30, True)},
        {"mate": 2},
        {"position_zero_based": True},
        {"position_zero_based": -2},
        {"cigar": ((0, 0),)},
        {"cigar": ((10, 4),)},
        {"cigar": ((0, 3),)},
    ],
)
def test_malformed_records_are_refused(changes):
    m = import_module("vntyper.scripts.calibration_read_identity")
    with pytest.raises(ValueError):
        m.read_identity_tokens(record(**changes))


def test_sorted_multiset_digest_preserves_multiplicity_and_refuses_unsorted_tokens():
    m = import_module("vntyper.scripts.calibration_read_identity")
    a = m.read_identity_tokens(record()).unnamed_sequence
    b = m.read_identity_tokens(record(sequence="ACGT")).unnamed_sequence
    once = m.digest_sorted_tokens(sorted([a, b]))
    twice = m.digest_sorted_tokens(sorted([a, a, b]))
    assert once.record_count == 2
    assert twice.record_count == 3
    assert once.sha256 != twice.sha256
    assert m.digest_sorted_tokens(iter(sorted([b, a]))) == once
    with pytest.raises(ValueError, match="sort"):
        m.digest_sorted_tokens(sorted([a, b], reverse=True))


@pytest.mark.parametrize("tokens", [[], ["not-sha256"], ["A" * 64], [True]])
def test_incomplete_or_malformed_multisets_have_no_usable_fingerprint(tokens):
    m = import_module("vntyper.scripts.calibration_read_identity")
    with pytest.raises(ValueError):
        m.digest_sorted_tokens(tokens)


def test_mutable_or_forged_record_fields_cannot_bypass_boundary_validation():
    m = import_module("vntyper.scripts.calibration_read_identity")
    original = record()
    with pytest.raises(ValueError):
        m.read_identity_tokens(replace(original, qualities=[10, 20, 30, 40]))
    with pytest.raises(ValueError):
        m.read_identity_tokens(replace(original, flags=129))


def test_unmapped_fastq_pair_content_matches_aligned_read_content():
    m = import_module("vntyper.scripts.calibration_read_identity")
    aligned = m.read_identity_tokens(record())
    unmapped = m.read_identity_tokens(
        record(
            flags=69,
            contig=None,
            position_zero_based=-1,
            cigar=(),
            mate_contig=None,
            mate_position_zero_based=-1,
            template_length=0,
        )
    )
    assert unmapped.named_sequence == aligned.named_sequence
    assert unmapped.unnamed_sequence == aligned.unnamed_sequence
    assert unmapped.alignment != aligned.alignment
    assert unmapped.reliable


@pytest.mark.parametrize(
    "changes",
    [
        {"contig": None},
        {"contig": ""},
        {"contig": None, "position_zero_based": -1},
        {"cigar": []},
        {"cigar": ((0,),)},
        {"cigar": ()},
    ],
)
def test_inconsistent_alignment_geometry_is_rejected(changes):
    m = import_module("vntyper.scripts.calibration_read_identity")
    with pytest.raises(ValueError):
        m.read_identity_tokens(record(**changes))


@pytest.mark.parametrize("tokens", [None, {"a" * 64}, {"a" * 64: 2}])
def test_nonstreams_and_deduplicating_containers_cannot_claim_multiset_identity(tokens):
    m = import_module("vntyper.scripts.calibration_read_identity")
    with pytest.raises(ValueError):
        m.digest_sorted_tokens(tokens)


def test_mapping_quality_is_alignment_evidence_but_not_source_read_content():
    m = import_module("vntyper.scripts.calibration_read_identity")
    high = m.read_identity_tokens(record(mapping_quality=60))
    low = m.read_identity_tokens(record(mapping_quality=0))
    assert high.alignment != low.alignment
    assert high.named_sequence == low.named_sequence


@pytest.mark.parametrize("mapq", [True, -1, 256, 1.5])
def test_invalid_mapping_quality_is_rejected(mapq):
    m = import_module("vntyper.scripts.calibration_read_identity")
    with pytest.raises(ValueError, match="mapping quality"):
        m.read_identity_tokens(record(mapping_quality=mapq))
