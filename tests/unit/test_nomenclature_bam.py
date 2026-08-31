"""Kestrel BAM rescue for alleles its VCF cannot represent.

Kestrel's pinned jar emits only 1-vs-1, 1-vs-N and N-vs-1 records, so a delins is
unrepresentable in its VCF. ``output.bam`` -- already produced, it is the report's
IGV track -- keeps resolved haplotype records aligned to the 120 bp pair reference
with full ``=``/``X``/``I``/``D`` CIGARs, and a ``1X1I`` block *is* a delins.

These tests build small BAMs with pysam rather than reading benchmark files, so the
tier stays in ``make test-unit`` and needs no downloaded data.

Research use only.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import pysam
import pytest

from vntyper.scripts.nomenclature import FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT, Nomenclature
from vntyper.scripts.nomenclature_bam import (
    BamConsensus,
    BamRescuer,
    from_bam,
    is_candidate,
    merge_edits,
    minimum_kmer_depth,
    refine,
    walk_cigar,
)

pytestmark = pytest.mark.unit

#: A 120 bp pair reference, long enough for the spans these tests use.
_REFERENCE = "ACGT" * 30
_TagType = Literal["c", "C", "s", "S", "i", "I", "f", "Z"]


@dataclass(frozen=True)
class _HaplotypeRecordSpec:
    """One synthetic resolved Kestrel haplotype record."""

    name: str
    start: int
    cigar: str
    minimum_kmer_depth: int | float | str | None = 181
    xd_type: _TagType | None = None


class _UnreadableXdRecord(pysam.AlignedSegment):
    """Real pysam record with injected auxiliary-tag decoding failure."""

    xd_error: Exception

    def get_tag(self, *_: object, **__: object) -> object:  # type: ignore[override]
        raise self.xd_error


def _haplotypes(
    prefix: str,
    count: int,
    cigar: str,
    depths: tuple[int | None, ...] | None = None,
) -> list[_HaplotypeRecordSpec]:
    """Build resolved haplotype-record specifications with realistic XD values."""
    observed_depths = depths or tuple((5, 181, 7_416, 8_704)[index % 4] for index in range(count))
    return [_HaplotypeRecordSpec(f"{prefix}{index}", 10, cigar, observed_depths[index]) for index in range(count)]


def _write_bam(path: Path, records: list[_HaplotypeRecordSpec]) -> Path:
    """Write a tiny indexed BAM.

    Args:
        path: Destination ``.bam``.
        records: Resolved haplotype-record specifications.

    Returns:
        Path: The written BAM.
    """
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "K-J", "LN": len(_REFERENCE)}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for specification in records:
            record = pysam.AlignedSegment(handle.header)
            record.query_name = specification.name
            record.reference_id = 0
            record.reference_start = specification.start
            record.mapping_quality = 255
            record.cigarstring = specification.cigar
            tuples = record.cigartuples or []
            length = sum(count for op, count in tuples if op in (0, 1, 4, 7, 8))
            record.query_sequence = "A" * length
            if specification.minimum_kmer_depth is not None:
                record.set_tag("XD", specification.minimum_kmer_depth, value_type=specification.xd_type)
            handle.write(record)
    pysam.index(str(path))  # type: ignore[attr-defined]
    return path


def _tagged_record(value: int | float | str | None, value_type: _TagType | None = None) -> pysam.AlignedSegment:
    """Build a real pysam record carrying one optional XD tag."""
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "K-J", "LN": 120}]})
    record = pysam.AlignedSegment(header)
    if value is not None:
        record.set_tag("XD", value, value_type=value_type)
    return record


def _unreadable_xd_record(error_type: type[Exception]) -> pysam.AlignedSegment:
    """Build a real pysam record whose XD decoding raises at the parser boundary."""
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "K-J", "LN": 120}]})
    record = _UnreadableXdRecord(header)
    record.xd_error = error_type("invalid auxiliary tag")
    return record


def _tier(tier: str, flags: tuple[str, ...] = ()) -> Nomenclature:
    return Nomenclature(
        name="59dupC" if tier != "C" else None,
        event="duplication",
        unit="X",
        tier=tier,
        flags=flags,
        ambiguity=None,
        repeat_form=None,
        net_length=1,
        source="kestrel_vcf",
    )


# ---------------------------------------------------------------------------
# Candidate selection -- the common path must pay nothing
# ---------------------------------------------------------------------------


def test_a_confident_call_is_not_a_candidate() -> None:
    assert is_candidate(_tier("A")) is False


def test_an_undetermined_call_is_a_candidate() -> None:
    assert is_candidate(_tier("C")) is True


def test_a_caller_disagreement_is_a_candidate() -> None:
    assert is_candidate(_tier("B", ("caller-disagreement",))) is True


def test_a_vcf_unrepresentable_allele_is_a_candidate() -> None:
    assert is_candidate(_tier("B", ("allele-unrepresentable-in-vcf",))) is True


def test_an_ordinary_tier_b_call_is_not_a_candidate() -> None:
    """Tier B is the common outcome; consulting the BAM for all of it is the cost
    the design exists to avoid."""
    assert is_candidate(_tier("B", ("position-ambiguous",))) is False


def test_the_rescuer_never_opens_the_bam_without_a_candidate(tmp_path: Path) -> None:
    """A BAM path that does not exist proves it was never opened."""
    rescuer = BamRescuer(tmp_path / "absent.bam")
    assert rescuer.opened is False
    assert rescuer.fetches == 0


# ---------------------------------------------------------------------------
# CIGAR walking and edit merging
# ---------------------------------------------------------------------------


def test_a_substitution_beside_an_insertion_merges_into_one_delins() -> None:
    """This is what recovers the allele the VCF cannot hold.

    ``1X1I`` replaces one reference base with two haplotype-record bases -- the shape of
    ``delinsAT`` -- which is precisely what Kestrel's 1-vs-1 / 1-vs-N / N-vs-1
    record set has no way to write down.
    """
    edits = walk_cigar(10, [(7, 5), (8, 1), (1, 1), (7, 5)])
    merged = merge_edits(edits)
    assert len(merged) == 1
    assert merged[0].kind == "delins"
    assert merged[0].ref_span == 1
    assert merged[0].inserted == 2


def test_separated_edits_are_not_merged() -> None:
    edits = walk_cigar(10, [(7, 5), (8, 1), (7, 5), (1, 1), (7, 5)])
    merged = merge_edits(edits)
    assert len(merged) == 2


def test_a_plain_insertion_stays_an_insertion() -> None:
    merged = merge_edits(walk_cigar(10, [(7, 5), (1, 3), (7, 5)]))
    assert len(merged) == 1
    assert merged[0].kind == "insertion"
    assert merged[0].inserted == 3


def test_a_plain_deletion_stays_a_deletion() -> None:
    merged = merge_edits(walk_cigar(10, [(7, 5), (2, 2), (7, 5)]))
    assert len(merged) == 1
    assert merged[0].kind == "deletion"
    assert merged[0].ref_span == 2


# ---------------------------------------------------------------------------
# XD parsing and record voting
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("value", "value_type"),
    [(5, "c"), (181, "C"), (7_416, "s"), (7_416, "S"), (181, "i"), (2_147_483_647, "I")],
)
def test_minimum_kmer_depth_accepts_pysam_integer_subtypes(value: int, value_type: _TagType) -> None:
    """pysam compacts SAM integer tags, so all BAM integer subtypes are valid."""
    assert minimum_kmer_depth(_tagged_record(value, value_type)) == value


def test_minimum_kmer_depth_returns_none_silently_when_xd_is_absent(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level("DEBUG", logger="vntyper.scripts.nomenclature_bam"):
        assert minimum_kmer_depth(_tagged_record(None)) is None
    assert "XD minimum k-mer depth" not in caplog.text


@pytest.mark.parametrize("error_type", [TypeError, ValueError, OverflowError])
def test_minimum_kmer_depth_contains_present_tag_decoding_errors(
    error_type: type[Exception],
    caplog: pytest.LogCaptureFixture,
) -> None:
    with caplog.at_level("DEBUG", logger="vntyper.scripts.nomenclature_bam"):
        assert minimum_kmer_depth(_unreadable_xd_record(error_type)) is None
    assert "Invalid XD minimum k-mer depth" in caplog.text


@pytest.mark.parametrize(("value", "value_type"), [(1.5, "f"), ("malformed", "Z")])
def test_minimum_kmer_depth_rejects_nonintegral_xd(
    value: float | str,
    value_type: _TagType,
    caplog: pytest.LogCaptureFixture,
) -> None:
    with caplog.at_level("DEBUG", logger="vntyper.scripts.nomenclature_bam"):
        assert minimum_kmer_depth(_tagged_record(value, value_type)) is None
    assert "Invalid XD minimum k-mer depth" in caplog.text


def test_minimum_kmer_depth_preserves_zero() -> None:
    assert minimum_kmer_depth(_tagged_record(0)) == 0


def test_minimum_kmer_depth_rejects_negative_xd(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level("DEBUG", logger="vntyper.scripts.nomenclature_bam"):
        assert minimum_kmer_depth(_tagged_record(-1)) is None
    assert "Invalid XD minimum k-mer depth" in caplog.text


def test_minimum_kmer_depth_preserves_the_largest_producer_integer() -> None:
    assert minimum_kmer_depth(_tagged_record(2_147_483_647)) == 2_147_483_647


def test_minimum_kmer_depth_rejects_unsigned_values_above_the_producer_range(
    caplog: pytest.LogCaptureFixture,
) -> None:
    with caplog.at_level("DEBUG", logger="vntyper.scripts.nomenclature_bam"):
        assert minimum_kmer_depth(_tagged_record(2_147_483_648, "I")) is None
    assert "Invalid XD minimum k-mer depth" in caplog.text


def test_only_length_changing_edits_are_voted_on(tmp_path: Path) -> None:
    """Substitutions in this reference are overwhelmingly motif polymorphism.

    A first attempt that voted on every edit was swamped by them, so the winner
    was whichever motif difference happened to be commonest.
    """
    bam = _write_bam(
        tmp_path / "subs.bam",
        _haplotypes("sub", 12, "20=1X20=") + _haplotypes("ins", 3, "20=1I20="),
    )
    consensus = BamRescuer(bam).rescue("K-J", 30)
    assert consensus is not None
    assert consensus.kind == "insertion"
    assert consensus.supporting_haplotype_records == 3


def test_support_is_reported_not_hidden(tmp_path: Path) -> None:
    """A thin consensus must be visible as thin, never silently accepted."""
    bam = _write_bam(
        tmp_path / "thin.bam",
        _haplotypes("variant", 1, "20=1X1I20=", (5,)) + _haplotypes("reference", 9, "41="),
    )
    consensus = BamRescuer(bam).rescue("K-J", 30)
    assert consensus is not None
    assert consensus.supporting_haplotype_records == 1
    assert consensus.fetched_haplotype_records == 10
    assert consensus.is_thin is True


def test_a_locus_with_no_length_changing_edit_yields_nothing(tmp_path: Path) -> None:
    bam = _write_bam(tmp_path / "flat.bam", _haplotypes("reference", 8, "41="))
    assert BamRescuer(bam).rescue("K-J", 30) is None


def test_distinct_alleles_are_counted(tmp_path: Path) -> None:
    bam = _write_bam(
        tmp_path / "mixed.bam",
        _haplotypes("one-base", 4, "20=1I20=") + _haplotypes("two-base", 2, "20=2I20="),
    )
    consensus = BamRescuer(bam).rescue("K-J", 30)
    assert consensus is not None
    assert consensus.distinct_edit_count == 2
    assert consensus.supporting_haplotype_records == 4


def test_changing_only_xd_does_not_change_a_three_to_two_record_winner(tmp_path: Path) -> None:
    low_winner = _write_bam(
        tmp_path / "low-winner.bam",
        _haplotypes("winner", 3, "20=1I20=", (5, 5, 5)) + _haplotypes("runner-up", 2, "20=2I20=", (8_704, 8_704)),
    )
    high_winner = _write_bam(
        tmp_path / "high-winner.bam",
        _haplotypes("winner", 3, "20=1I20=", (8_704, None, 0))
        + _haplotypes("runner-up", 2, "20=2I20=", (-1, 2_147_483_647)),
    )

    low_consensus = BamRescuer(low_winner).rescue("K-J", 30)
    high_consensus = BamRescuer(high_winner).rescue("K-J", 30)

    assert low_consensus is not None
    assert high_consensus is not None
    assert low_consensus == high_consensus
    assert low_consensus.supporting_record_minimum_kmer_depths == (5, 5, 5)
    assert high_consensus.supporting_record_minimum_kmer_depths == (8_704, None, 0)


def test_unequal_xd_cannot_break_an_equal_record_vote(tmp_path: Path) -> None:
    bam = _write_bam(
        tmp_path / "xd-tie.bam",
        _haplotypes("weak", 2, "20=1I20=", (5, 5)) + _haplotypes("strong", 2, "20=2I20=", (8_704, 8_704)),
    )
    assert BamRescuer(bam).rescue("K-J", 30) is None


def test_several_weak_xd_records_outvote_one_maximum_xd_record(tmp_path: Path) -> None:
    bam = _write_bam(
        tmp_path / "weak-majority.bam",
        _haplotypes("weak", 3, "20=1I20=", (5, 5, 5)) + _haplotypes("strong", 1, "20=2I20=", (2_147_483_647,)),
    )
    consensus = BamRescuer(bam).rescue("K-J", 30)
    assert consensus is not None
    assert consensus.inserted == 1
    assert consensus.supporting_haplotype_records == 3


def test_equal_xd_sums_do_not_override_the_record_count_winner(tmp_path: Path) -> None:
    bam = _write_bam(
        tmp_path / "equal-xd-sums.bam",
        _haplotypes("two-record-winner", 2, "20=1I20=", (5, 5))
        + _haplotypes("one-record-runner-up", 1, "20=2I20=", (10,)),
    )
    consensus = BamRescuer(bam).rescue("K-J", 30)
    assert consensus is not None
    assert consensus.inserted == 1
    assert consensus.supporting_haplotype_records == 2
    assert consensus.supporting_record_minimum_kmer_depths == (5, 5)


# ---------------------------------------------------------------------------
# Handle reuse
# ---------------------------------------------------------------------------


def test_one_handle_serves_many_loci(tmp_path: Path) -> None:
    """Re-opening per locus is the cost the design calls out; prove it does not."""
    bam = _write_bam(tmp_path / "reuse.bam", _haplotypes("variant", 6, "20=1I20="))
    rescuer = BamRescuer(bam)
    for _ in range(5):
        rescuer.rescue("K-J", 30)
    assert rescuer.opened is True
    assert rescuer.opens == 1
    assert rescuer.fetches == 5
    rescuer.close()
    assert rescuer.opened is False


def test_an_unknown_contig_is_not_an_error(tmp_path: Path) -> None:
    bam = _write_bam(tmp_path / "unknown.bam", _haplotypes("reference", 1, "41="))
    assert BamRescuer(bam).rescue("NOT-A-PAIR", 30) is None


def test_a_missing_bam_is_not_an_error(tmp_path: Path) -> None:
    assert BamRescuer(tmp_path / "gone.bam").rescue("K-J", 30) is None


# ---------------------------------------------------------------------------
# Consensus value object
# ---------------------------------------------------------------------------


def test_a_tie_yields_no_rescue_at_all(tmp_path: Path) -> None:
    """A tie has no winner. `Counter.most_common` breaks one by BAM order, so
    whichever haplotype record was written first would decide an allele that can
    override the VCF."""
    bam = _write_bam(
        tmp_path / "tied.bam",
        _haplotypes("one-base", 3, "20=1I20=", (5, 5, 5))
        + _haplotypes("two-base", 3, "20=2I20=", (8_704, 8_704, 8_704)),
    )
    assert BamRescuer(bam).rescue("K-J", 30) is None


def test_a_clear_winner_over_a_runner_up_is_still_rescued(tmp_path: Path) -> None:
    bam = _write_bam(
        tmp_path / "clear.bam",
        _haplotypes("one-base", 5, "20=1I20=") + _haplotypes("two-base", 2, "20=2I20="),
    )
    consensus = BamRescuer(bam).rescue("K-J", 30)
    assert consensus is not None
    assert consensus.supporting_haplotype_records == 5


def test_a_well_supported_consensus_is_not_thin() -> None:
    assert (
        BamConsensus(
            "delins",
            55,
            1,
            2,
            supporting_haplotype_records=9,
            fetched_haplotype_records=20,
            distinct_edit_count=2,
            supporting_record_minimum_kmer_depths=(5, 181, 7_416, 8_704, 5, 181, 7_416, 8_704, None),
        ).is_thin
        is False
    )


def test_xd_is_equality_and_hash_neutral_and_cannot_change_thinness() -> None:
    weak = BamConsensus(
        "delins",
        55,
        1,
        2,
        2,
        20,
        2,
        supporting_record_minimum_kmer_depths=(5, 5),
    )
    extreme = BamConsensus(
        "delins",
        55,
        1,
        2,
        2,
        20,
        2,
        supporting_record_minimum_kmer_depths=(2_147_483_647, None),
    )

    assert weak == extreme
    assert hash(weak) == hash(extreme)
    assert weak.is_thin is True
    assert extreme.is_thin is True
    assert weak.support == weak.supporting_haplotype_records
    assert weak.total == weak.fetched_haplotype_records
    assert weak.n_distinct == weak.distinct_edit_count


# ---------------------------------------------------------------------------
# VCF primary, BAM refines
# ---------------------------------------------------------------------------


def _named(name: str, event: str, source: str = "kestrel_vcf") -> Nomenclature:
    return Nomenclature(
        name=name,
        event=event,
        unit="X",
        tier="B",
        flags=(),
        ambiguity=None,
        repeat_form=None,
        net_length=1,
        source=source,
    )


def test_the_bam_supplies_an_allele_the_vcf_had_none_for() -> None:
    """The whole point of the rescue path."""
    undetermined = _tier("C")
    refined = refine(undetermined, _named("58_59insG", "insertion", "kestrel_bam"))
    assert refined.name == "58_59insG"
    assert refined.source == "kestrel_bam"


def test_from_bam_marks_thin_haplotype_record_support() -> None:
    consensus = BamConsensus(
        "insertion",
        60,
        0,
        1,
        supporting_haplotype_records=2,
        fetched_haplotype_records=10,
        distinct_edit_count=1,
        bases="A",
        supporting_record_minimum_kmer_depths=(5, 8_704),
    )
    named = from_bam("K-J", consensus)
    assert named is not None
    assert FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT in named.flags
    assert "low-read-support" not in named.flags


def test_a_thin_haplotype_consensus_cannot_settle_caller_disagreement() -> None:
    disagreement = Nomenclature(
        name=None,
        event="disagreement",
        unit="X",
        tier="C",
        flags=("caller-disagreement",),
        ambiguity=None,
        repeat_form=None,
        net_length=0,
        source="reconciled",
    )
    thin = _named("58_59insG", "insertion", "kestrel_bam")
    thin = Nomenclature(
        name=thin.name,
        event=thin.event,
        unit=thin.unit,
        tier=thin.tier,
        flags=(FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT,),
        ambiguity=thin.ambiguity,
        repeat_form=thin.repeat_form,
        net_length=thin.net_length,
        source=thin.source,
    )
    assert refine(disagreement, thin) is disagreement


def test_the_bam_may_not_veto_an_allele_the_vcf_already_has() -> None:
    """A thin haplotype-record consensus cannot overrule an existing VCF name."""
    vcf = _named("60dupA", "duplication")
    refined = refine(vcf, _named("59delC", "deletion", "kestrel_bam"))
    assert refined.name == "60dupA"
    assert "caller-disagreement" in refined.flags


def test_agreement_leaves_the_call_alone() -> None:
    vcf = _named("59dupC", "duplication")
    assert refine(vcf, _named("59dupC", "duplication", "kestrel_bam")) is vcf


def test_a_delins_from_haplotype_records_overrides_a_shape_the_vcf_cannot_hold() -> None:
    """Kestrel's VariantType has SNP, INSERTION and DELETION and nothing else, so
    whatever it wrote for a delins locus is the closest representable shape rather
    than the allele."""
    vcf = _named("55_56insA", "insertion")
    refined = refine(vcf, _named("55delinsAT", "delins", "kestrel_bam"))
    assert refined.name == "55delinsAT"
    assert "allele-unrepresentable-in-vcf" in refined.flags


def test_silence_from_haplotype_records_changes_nothing() -> None:
    vcf = _named("59dupC", "duplication")
    assert refine(vcf, None) is vcf
