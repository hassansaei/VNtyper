"""BAM rescue for alleles the Kestrel VCF cannot represent.

Kestrel's pinned jar emits only 1-vs-1, 1-vs-N and N-vs-1 records, so a delins is
unrepresentable in its VCF. ``output.bam`` -- already produced, it is the report's
IGV track -- keeps the reads aligned to the 120 bp pair reference with full
``=``/``X``/``I``/``D`` CIGARs, and a ``1X1I`` block *is* a delins.

These tests build small BAMs with pysam rather than reading benchmark files, so the
tier stays in ``make test-unit`` and needs no downloaded data.

Research use only.
"""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from vntyper.scripts.nomenclature import Nomenclature
from vntyper.scripts.nomenclature_bam import (
    BamConsensus,
    BamRescuer,
    is_candidate,
    merge_edits,
    refine,
    walk_cigar,
)

pytestmark = pytest.mark.unit

#: A 120 bp pair reference, long enough for the spans these tests use.
_REFERENCE = "ACGT" * 30


def _write_bam(path: Path, reads: list[tuple[str, int, str]]) -> Path:
    """Write a tiny indexed BAM.

    Args:
        path: Destination ``.bam``.
        reads: ``(name, 0-based start, cigar string)`` per read.

    Returns:
        Path: The written BAM.
    """
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "K-J", "LN": len(_REFERENCE)}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for name, start, cigar in reads:
            record = pysam.AlignedSegment(handle.header)
            record.query_name = name
            record.reference_id = 0
            record.reference_start = start
            record.mapping_quality = 255
            record.cigarstring = cigar
            tuples = record.cigartuples or []
            length = sum(count for op, count in tuples if op in (0, 1, 4, 7, 8))
            record.query_sequence = "A" * length
            handle.write(record)
    pysam.index(str(path))  # type: ignore[attr-defined]
    return path


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

    ``1X1I`` replaces one reference base with two read bases -- the shape of
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
# Voting
# ---------------------------------------------------------------------------


def test_only_length_changing_edits_are_voted_on(tmp_path: Path) -> None:
    """Substitutions in this reference are overwhelmingly motif polymorphism.

    A first attempt that voted on every edit was swamped by them, so the winner
    was whichever motif difference happened to be commonest.
    """
    bam = _write_bam(
        tmp_path / "subs.bam",
        [(f"r{n}", 10, "20=1X20=") for n in range(12)] + [(f"i{n}", 10, "20=1I20=") for n in range(3)],
    )
    consensus = BamRescuer(bam).rescue("K-J", 30)
    assert consensus is not None
    assert consensus.kind == "insertion"
    assert consensus.support == 3


def test_support_is_reported_not_hidden(tmp_path: Path) -> None:
    """A thin consensus must be visible as thin, never silently accepted."""
    bam = _write_bam(
        tmp_path / "thin.bam",
        [("only", 10, "20=1X1I20=")] + [(f"m{n}", 10, "41=") for n in range(9)],
    )
    consensus = BamRescuer(bam).rescue("K-J", 30)
    assert consensus is not None
    assert consensus.support == 1
    assert consensus.total == 10
    assert consensus.is_thin is True


def test_a_locus_with_no_length_changing_edit_yields_nothing(tmp_path: Path) -> None:
    bam = _write_bam(tmp_path / "flat.bam", [(f"r{n}", 10, "41=") for n in range(8)])
    assert BamRescuer(bam).rescue("K-J", 30) is None


def test_distinct_alleles_are_counted(tmp_path: Path) -> None:
    bam = _write_bam(
        tmp_path / "mixed.bam",
        [(f"a{n}", 10, "20=1I20=") for n in range(4)] + [(f"b{n}", 10, "20=2I20=") for n in range(2)],
    )
    consensus = BamRescuer(bam).rescue("K-J", 30)
    assert consensus is not None
    assert consensus.n_distinct == 2
    assert consensus.support == 4


# ---------------------------------------------------------------------------
# Handle reuse
# ---------------------------------------------------------------------------


def test_one_handle_serves_many_loci(tmp_path: Path) -> None:
    """Re-opening per locus is the cost the design calls out; prove it does not."""
    bam = _write_bam(tmp_path / "reuse.bam", [(f"r{n}", 10, "20=1I20=") for n in range(6)])
    rescuer = BamRescuer(bam)
    for _ in range(5):
        rescuer.rescue("K-J", 30)
    assert rescuer.opened is True
    assert rescuer.opens == 1
    assert rescuer.fetches == 5
    rescuer.close()
    assert rescuer.opened is False


def test_an_unknown_contig_is_not_an_error(tmp_path: Path) -> None:
    bam = _write_bam(tmp_path / "unknown.bam", [("r", 10, "41=")])
    assert BamRescuer(bam).rescue("NOT-A-PAIR", 30) is None


def test_a_missing_bam_is_not_an_error(tmp_path: Path) -> None:
    assert BamRescuer(tmp_path / "gone.bam").rescue("K-J", 30) is None


# ---------------------------------------------------------------------------
# Consensus value object
# ---------------------------------------------------------------------------


def test_a_tie_yields_no_rescue_at_all(tmp_path: Path) -> None:
    """A tie has no winner. `Counter.most_common` breaks one by BAM order, so
    whichever read was written first would decide an allele that can override the
    VCF."""
    bam = _write_bam(
        tmp_path / "tied.bam",
        [(f"a{n}", 10, "20=1I20=") for n in range(3)] + [(f"b{n}", 10, "20=2I20=") for n in range(3)],
    )
    assert BamRescuer(bam).rescue("K-J", 30) is None


def test_a_clear_winner_over_a_runner_up_is_still_rescued(tmp_path: Path) -> None:
    bam = _write_bam(
        tmp_path / "clear.bam",
        [(f"a{n}", 10, "20=1I20=") for n in range(5)] + [(f"b{n}", 10, "20=2I20=") for n in range(2)],
    )
    consensus = BamRescuer(bam).rescue("K-J", 30)
    assert consensus is not None
    assert consensus.support == 5


def test_a_well_supported_consensus_is_not_thin() -> None:
    assert BamConsensus("delins", 55, 1, 2, support=9, total=20, n_distinct=2).is_thin is False


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


def test_the_bam_may_not_veto_an_allele_the_vcf_already_has() -> None:
    """Measured: letting a thin read consensus overrule a VCF name cost more than
    it gained -- dupA fell from 6 correct to 1, insCCCC from 9 to 3."""
    vcf = _named("60dupA", "duplication")
    refined = refine(vcf, _named("59delC", "deletion", "kestrel_bam"))
    assert refined.name == "60dupA"
    assert "caller-disagreement" in refined.flags


def test_agreement_leaves_the_call_alone() -> None:
    vcf = _named("59dupC", "duplication")
    assert refine(vcf, _named("59dupC", "duplication", "kestrel_bam")) is vcf


def test_a_delins_from_the_reads_overrides_a_shape_the_vcf_cannot_hold() -> None:
    """Kestrel's VariantType has SNP, INSERTION and DELETION and nothing else, so
    whatever it wrote for a delins locus is the closest representable shape rather
    than the allele."""
    vcf = _named("55_56insA", "insertion")
    refined = refine(vcf, _named("55delinsAT", "delins", "kestrel_bam"))
    assert refined.name == "55delinsAT"
    assert "allele-unrepresentable-in-vcf" in refined.flags


def test_silence_from_the_reads_changes_nothing() -> None:
    vcf = _named("59dupC", "duplication")
    assert refine(vcf, None) is vcf
