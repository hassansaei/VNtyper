"""Recover MUC1 alleles the Kestrel VCF cannot represent, from resolved haplotypes.

Kestrel's pinned jar (1.0.1) emits only 1-vs-1, 1-vs-N and N-vs-1 records -- its
``VariantType`` has SNP, INSERTION and DELETION and nothing else -- so a delins is
unrepresentable in its VCF. ``output.bam`` is already produced (it is the report's
IGV track) and keeps resolved haplotype records aligned to the 120 bp pair
reference with full ``=``/``X``/``I``/``D`` CIGARs. A ``1X1I`` block **is** a
delins.

Kept separate from :mod:`vntyper.scripts.nomenclature` so the translator stays pure
and free of pysam.

Policy: **VCF primary, BAM refines.** The BAM is consulted only for rows the VCF
shape cannot represent, where the callers disagree, or where the tier logic would
otherwise emit tier C -- a minority of rows, so the common path pays nothing.

Research use only.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING

import pysam

from vntyper.scripts.molecular_identity import MolecularIdentity
from vntyper.scripts.nomenclature import (
    FLAG_ALLELE_UNREPRESENTABLE,
    FLAG_CALLER_DISAGREEMENT,
    FLAG_SEQUENCE_UNDETERMINED,
    FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT,
    Nomenclature,
    name_coding_pair_edit,
    nomenclature_config,
    pair_sequence,
    revcomp,
)
from vntyper.scripts.nomenclature_bam_evidence import (
    BamEditObservation,
    BamIdentityLocus,
    BamIdentityTranslator,
    BamLocusEvidence,
    bind_complete_winner_identity,
    collect_locus_evidence,
    project_record_observation,
)
from vntyper.scripts.nomenclature_evidence import resolve_bam_thin_haplotype_record_support

if TYPE_CHECKING:  # pragma: no cover - typing only
    from collections.abc import Iterable

logger = logging.getLogger(__name__)

__all__ = [
    "DEFAULT_FLANK",
    "BamConsensus",
    "BamRescuer",
    "Edit",
    "THIN_HAPLOTYPE_RECORD_SUPPORT",
    "from_bam",
    "is_candidate",
    "merge_edits",
    "minimum_kmer_depth",
    "refine",
    "walk_cigar",
]

#: Window either side of the called locus, in bases. Configured, not hard-coded.
DEFAULT_FLANK: int = nomenclature_config["thresholds"]["bam_flank"]

#: Resolved haplotype records below which a consensus is reported as thin.
THIN_HAPLOTYPE_RECORD_SUPPORT: int = resolve_bam_thin_haplotype_record_support(nomenclature_config["thresholds"])

#: Maximum XD emitted by the pinned Java implementation's signed integer field.
_MAXIMUM_MINIMUM_KMER_DEPTH = 2_147_483_647
_BAM_INTEGER_TAG_TYPES = frozenset({"c", "C", "s", "S", "i", "I"})

# CIGAR operation codes, as pysam returns them.
_OP_MATCH = 0
_OP_INS = 1
_OP_DEL = 2
_OP_SOFT_CLIP = 4
_OP_EQUAL = 7
_OP_DIFF = 8

_CONSUMES_REFERENCE = frozenset({_OP_MATCH, _OP_DEL, _OP_EQUAL, _OP_DIFF})
_CONSUMES_QUERY = frozenset({_OP_MATCH, _OP_INS, _OP_SOFT_CLIP, _OP_EQUAL, _OP_DIFF})
_NON_MATCH = frozenset({_OP_INS, _OP_DEL, _OP_DIFF})


@dataclass(frozen=True)
class Edit:
    """One non-matching CIGAR block, anchored on the reference.

    Attributes:
        kind: ``substitution`` | ``insertion`` | ``deletion`` | ``delins``.
        start: 0-based reference position the block starts at.
        ref_span: Reference bases consumed.
        inserted: Haplotype-record bases inserted.
        bases: The inserted haplotype-record bases. Empty when the sequence was not
            supplied -- the *count* is enough to vote, but naming the allele needs
            the sequence, which is the whole reason the BAM is consulted at all.
    """

    kind: str
    start: int
    ref_span: int
    inserted: int
    bases: str = ""

    @property
    def changes_length(self) -> bool:
        """Whether this edit alters the allele's length."""
        return self.ref_span != self.inserted


@dataclass(frozen=True)
class BamConsensus:
    """The winning length-changing edit at a locus, with its evidence.

    Attributes:
        kind: The edit class.
        start: 0-based reference position.
        ref_span: Reference bases consumed.
        inserted: Haplotype-record bases inserted.
        supporting_haplotype_records: Resolved haplotype records carrying this
            exact edit.
        fetched_haplotype_records: Resolved haplotype records examined at the locus.
        distinct_edit_count: Distinct length-changing edits seen.
        bases: Inserted bases on the genomic plus strand.
        supporting_record_minimum_kmer_depths: Optional XD minimum k-mer depth for
            each supporting haplotype record, in encounter order. Observational only.
        bound_identity: Existing complete Kestrel candidate identity proven for
            this exact-edit winner, or ``None``. Binding is equality-neutral and
            creates no representation or caller.
    """

    kind: str
    start: int
    ref_span: int
    inserted: int
    supporting_haplotype_records: int
    fetched_haplotype_records: int
    distinct_edit_count: int
    bases: str = ""
    supporting_record_minimum_kmer_depths: tuple[int | None, ...] = field(kw_only=True, compare=False)
    bound_identity: MolecularIdentity | None = field(default=None, kw_only=True, compare=False)

    @property
    def support(self) -> int:
        """Return the winning support through the pre-Phase-1 attribute name."""
        return self.supporting_haplotype_records

    @property
    def total(self) -> int:
        """Return fetched record count through the pre-Phase-1 attribute name."""
        return self.fetched_haplotype_records

    @property
    def n_distinct(self) -> int:
        """Return distinct edit count through the pre-Phase-1 attribute name."""
        return self.distinct_edit_count

    @property
    def is_thin(self) -> bool:
        """Whether the consensus rests on too few haplotype records to stand alone."""
        return self.supporting_haplotype_records < THIN_HAPLOTYPE_RECORD_SUPPORT


@dataclass(frozen=True)
class _LocusRecord:
    """One fetched record's merged edits and observational XD."""

    edits: tuple[Edit, ...]
    minimum_kmer_depth: int | None


def minimum_kmer_depth(record: pysam.AlignedSegment) -> int | None:
    """Parse optional XD minimum k-mer depth from a resolved haplotype record.

    XD is observational evidence only. Missing or invalid XD does not invalidate
    the record and cannot affect its vote.

    Args:
        record: One resolved haplotype record from Kestrel ``output.bam``.

    Returns:
        int | None: The validated signed Java integer, including zero, or ``None``
        when XD is absent or invalid.
    """
    try:
        value, value_type = record.get_tag("XD", with_value_type=True)
    except KeyError:
        return None
    except (TypeError, ValueError, OverflowError) as error:
        logger.debug("Invalid XD minimum k-mer depth on haplotype record %s: %s", record.query_name, error)
        return None

    if type(value) is not int or value_type not in _BAM_INTEGER_TAG_TYPES:
        logger.debug(
            "Invalid XD minimum k-mer depth on haplotype record %s: value=%r type=%r",
            record.query_name,
            value,
            value_type,
        )
        return None
    if value < 0 or value > _MAXIMUM_MINIMUM_KMER_DEPTH:
        logger.debug(
            "Invalid XD minimum k-mer depth on haplotype record %s: value=%r outside 0..%d",
            record.query_name,
            value,
            _MAXIMUM_MINIMUM_KMER_DEPTH,
        )
        return None
    return value


def from_bam(motifs: str, consensus: BamConsensus) -> Nomenclature | None:
    """Name the allele a resolved haplotype-record consensus describes.

    The consensus arrives in plus-strand pair coordinates, exactly like a Kestrel
    VCF record, and is handed to the same naming machinery -- a name recovered from
    the BAM must be produced the way a name from the VCF is, or the two would
    disagree for reasons that have nothing to do with the evidence.

    Args:
        motifs: The pair label from ``output.bed``, e.g. ``K-J``.
        consensus: The winning length-changing edit at the locus.

    Returns:
        Nomenclature | None: The named allele, or ``None`` when the pair is unknown.
        A thin consensus is named but flagged
        ``thin-haplotype-record-support``; it is never silently accepted or
        promoted on its own.
    """
    pair = pair_sequence(motifs)
    if pair is None:
        return None

    pair_length = len(pair)
    net = consensus.inserted - consensus.ref_span

    # `consensus.start` is 0-based on the plus strand; the replaced span is
    # [start+1, start+ref_span] in 1-based plus coordinates.
    plus_lo = consensus.start + 1
    plus_hi = consensus.start + consensus.ref_span

    if consensus.ref_span == 0:
        # A pure insertion sits in the gap after `consensus.start`. Reverse
        # complement swaps which side of the gap it lies on.
        coding_left = pair_length - consensus.start
        coding_start, coding_end = coding_left + 1, coding_left
    else:
        coding_start = pair_length + 1 - plus_hi
        coding_end = pair_length + 1 - plus_lo

    named = name_coding_pair_edit(
        motifs,
        pair,
        coding_start,
        coding_end,
        revcomp(consensus.bases),
        net,
        "kestrel_bam",
    )

    if consensus.is_thin:
        named = Nomenclature(
            name=named.name,
            event=named.event,
            unit=named.unit,
            tier=named.tier,
            flags=tuple(sorted({*named.flags, FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT})),
            ambiguity=named.ambiguity,
            repeat_form=named.repeat_form,
            net_length=named.net_length,
            source=named.source,
        )
    return named


def refine(call: Nomenclature, bam_call: Nomenclature | None) -> Nomenclature:
    """Apply the BAM evidence to a call. **VCF primary, BAM refines.**

    The BAM may supply an allele the VCF had none for, and may corroborate one the
    VCF already has. It may **not** veto one. Measured on the benchmark, letting a
    haplotype-record consensus overrule an existing name cost more than it gained:
    a locus often carries only 1-3 resolved records, and a thin consensus
    disagreeing with a well-supported VCF record destroyed names that were right
    -- ``dupA`` fell from 6 to 1 and ``insCCCC`` from 9 to 3.

    Args:
        call: The reconciled VCF/adVNTR call.
        bam_call: What the resolved haplotypes say, or ``None`` when they say nothing.

    Returns:
        Nomenclature: The refined call. Unchanged when the BAM adds nothing.
    """
    if bam_call is None or bam_call.name is None:
        return call

    if call.name is None:
        if FLAG_CALLER_DISAGREEMENT in call.flags and FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT in bam_call.flags:
            # There is no name here because two callers described different events,
            # which is a real conflict rather than a gap the BAM may fill. A thin
            # consensus -- one or two haplotype records -- is not
            # enough to settle it, and adopting it turned an honest "allele
            # undetermined" into a bare number backed by a single record.
            return call
        # The VCF had no allele to offer; the resolved haplotypes do. This is the
        # whole point of the rescue path -- it is where delins and the insG family
        # come from.
        return Nomenclature(
            name=bam_call.name,
            event=bam_call.event,
            unit=bam_call.unit,
            tier=bam_call.tier,
            flags=tuple(sorted({*call.flags, *bam_call.flags})),
            ambiguity=bam_call.ambiguity,
            repeat_form=bam_call.repeat_form,
            net_length=bam_call.net_length,
            source="kestrel_bam",
        )

    if bam_call.name == call.name:
        # Corroboration. Recorded as agreement, but tier promotion still belongs to
        # `reconcile`, which is the only place that sees the support quantity.
        return call

    if bam_call.event == "delins" and call.event != "delins":
        # The VCF shape cannot hold this. Kestrel's `VariantType` has SNP, INSERTION
        # and DELETION and nothing else, so whatever it wrote for this locus is the
        # closest representable shape rather than the allele -- the haplotype
        # records are the better evidence here even though they are the junior
        # source everywhere else.
        return Nomenclature(
            name=bam_call.name,
            event=bam_call.event,
            unit=bam_call.unit,
            tier=bam_call.tier,
            flags=tuple(sorted({*call.flags, *bam_call.flags, FLAG_ALLELE_UNREPRESENTABLE})),
            ambiguity=bam_call.ambiguity,
            repeat_form=bam_call.repeat_form,
            net_length=bam_call.net_length,
            source="kestrel_bam",
        )

    # Disagreement: keep the VCF allele and say so, rather than silently preferring
    # either. The flag is what stops this reaching tier A.
    return Nomenclature(
        name=call.name,
        event=call.event,
        unit=call.unit,
        tier=call.tier,
        flags=tuple(sorted({*call.flags, FLAG_CALLER_DISAGREEMENT})),
        ambiguity=call.ambiguity,
        repeat_form=call.repeat_form,
        net_length=call.net_length,
        source=call.source,
    )


def is_candidate(call: Nomenclature) -> bool:
    """Should the BAM be consulted for this call?

    The BAM is not scanned per row. It is consulted only where the VCF cannot hold
    what its resolved haplotype records show, where the callers disagree, or where
    the result would otherwise be tier C. An ordinary tier-B call -- the common
    outcome -- is not a candidate, which is what keeps the common path free.

    Args:
        call: The reconciled call.

    Returns:
        bool: True when the BAM may add something the VCF could not.
    """
    if call.tier == "A":
        return False
    if call.tier == "C":
        return True
    return bool({FLAG_CALLER_DISAGREEMENT, FLAG_ALLELE_UNREPRESENTABLE, FLAG_SEQUENCE_UNDETERMINED} & set(call.flags))


def walk_cigar(reference_start: int, cigar: Iterable[tuple[int, int]], sequence: str | None = None) -> list[Edit]:
    """Turn one haplotype record's CIGAR into reference-anchored non-matching blocks.

    A single pass, O(len(cigar)); no per-record reference slicing beyond the edit
    spans themselves.

    Args:
        reference_start: 0-based alignment start.
        cigar: ``(operation, length)`` pairs as pysam returns them.
        sequence: The haplotype sequence, when the inserted bases are wanted.
            Voting only needs lengths, but naming the allele needs the bases.

    Returns:
        list[Edit]: The non-matching blocks, in reference order.
    """
    edits: list[Edit] = []
    position = reference_start
    query = 0
    for operation, length in cigar:
        if operation in _NON_MATCH:
            bases = ""
            if sequence is not None and operation in (_OP_INS, _OP_DIFF):
                bases = sequence[query : query + length]
            if operation == _OP_INS:
                edits.append(Edit("insertion", position, 0, length, bases))
            elif operation == _OP_DEL:
                edits.append(Edit("deletion", position, length, 0, ""))
            else:
                edits.append(Edit("substitution", position, length, length, bases))
        if operation in _CONSUMES_REFERENCE:
            position += length
        if operation in _CONSUMES_QUERY:
            query += length
    return edits


def merge_edits(edits: Iterable[Edit]) -> list[Edit]:
    """Merge adjacent non-matching blocks into single events.

    This is what turns a ``1X1I`` pair into a delins -- one reference base replaced
    by two haplotype-record bases -- and it is the only reason the BAM can express
    an allele the VCF cannot.

    Args:
        edits: Reference-ordered blocks from :func:`walk_cigar`.

    Returns:
        list[Edit]: Merged events.
    """
    merged: list[Edit] = []
    for edit in edits:
        if merged:
            previous = merged[-1]
            if previous.start + previous.ref_span == edit.start:
                ref_span = previous.ref_span + edit.ref_span
                inserted = previous.inserted + edit.inserted
                kind = "delins" if ref_span and inserted and ref_span != inserted else previous.kind
                if ref_span and not inserted:
                    kind = "deletion"
                elif inserted and not ref_span:
                    kind = "insertion"
                merged[-1] = Edit(kind, previous.start, ref_span, inserted, previous.bases + edit.bases)
                continue
        merged.append(edit)
    return merged


class BamRescuer:
    """One open handle over a sample's ``output.bam``, reused across loci.

    Re-opening per locus is a real cost on a cohort, so the handle is opened lazily
    on the first fetch and kept. ``opens`` and ``fetches`` are exposed so a test can
    assert the BAM is never touched when there is no candidate -- the budget the
    design puts on this path is only meaningful if it is measured.
    """

    def __init__(self, bam_path: str | Path, flank: int = DEFAULT_FLANK) -> None:
        """
        Args:
            bam_path: Path to ``output.bam``.
            flank: Window either side of the locus, in bases.
        """
        self._path = Path(bam_path)
        self._flank = flank
        self._handle: pysam.AlignmentFile | None = None
        self.opens = 0
        self.fetches = 0

    @property
    def opened(self) -> bool:
        """Whether a handle is currently open."""
        return self._handle is not None

    def _open(self) -> pysam.AlignmentFile | None:
        """Open the BAM once, or report why it cannot be read.

        Returns:
            pysam.AlignmentFile | None: The handle, or ``None`` when unavailable.
        """
        if self._handle is not None:
            return self._handle
        if not self._path.is_file():
            logger.debug("BAM rescue skipped; no file at %s", self._path)
            return None
        try:
            self._handle = pysam.AlignmentFile(str(self._path), "rb")
        except (OSError, ValueError) as error:
            logger.debug("BAM rescue skipped; %s could not be opened: %s", self._path, error)
            return None
        self.opens += 1
        return self._handle

    def close(self) -> None:
        """Release the handle."""
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def __enter__(self) -> BamRescuer:
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def rescue(self, contig: str, position: int) -> BamConsensus | None:
        """Resolve the consensus length-changing edit at one locus.

        Only length-changing edits are voted on. Substitutions against this
        reference are overwhelmingly motif polymorphism -- a first attempt that
        voted on every edit was swamped by them, and the winner was whichever motif
        difference happened to be commonest rather than the variant.

        Args:
            contig: Pair reference name, e.g. ``K-J``, from ``output.bed``.
            position: 1-based locus position from ``output.bed``.

        Returns:
            BamConsensus | None: The winning edit with its evidence, or ``None`` when
            the locus is unreadable or carries no length-changing edit.
        """
        fetched = self._fetch_locus(contig, position)
        if fetched is None:
            return None
        records, start, end = fetched
        return self._compatibility_consensus(records, start, end, contig, position)

    def rescue_with_identity_evidence(
        self,
        contig: str,
        position: int,
        locus: BamIdentityLocus,
        component: BamIdentityTranslator,
    ) -> tuple[BamConsensus | None, BamLocusEvidence | None]:
        """Return the unchanged rescue verdict beside replayable identity evidence.

        This adapter fetches and parses each eligible record once. Its first result
        is produced by the established exact-edit record vote; its second binds the
        same records only to existing complete Kestrel candidate representations.
        A legacy tie can therefore return ``None`` while retaining every tied record
        for later replay.

        Args:
            contig: Pair reference name from ``output.bed``.
            position: One-based locus position.
            locus: Complete pair context and existing Kestrel candidates.
            component: Injected complete-context identity translator.

        Returns:
            ``(consensus, evidence)``. Missing or unreadable BAM evidence returns
            ``(None, None)`` rather than a fabricated zero-record locus.
        """
        fetched = self._fetch_locus(contig, position)
        if fetched is None:
            return None, None
        records, start, end = fetched
        consensus = self._compatibility_consensus(records, start, end, contig, position)
        observations = tuple(
            project_record_observation(
                tuple((edit.start, edit.ref_span, edit.inserted, edit.bases) for edit in record.edits),
                record.minimum_kmer_depth,
                start,
                end,
            )
            for record in records
        )
        evidence = collect_locus_evidence(observations, locus, component)
        if consensus is not None:
            winning_observation = BamEditObservation(
                consensus.start,
                consensus.ref_span,
                consensus.inserted,
                consensus.bases,
            )
            bound_identity = bind_complete_winner_identity(
                observations,
                evidence,
                winning_observation,
                locus,
                component,
                consensus.supporting_haplotype_records,
            )
            if bound_identity is not None:
                consensus = replace(consensus, bound_identity=bound_identity)
        return consensus, evidence

    def _fetch_locus(self, contig: str, position: int) -> tuple[tuple[_LocusRecord, ...], int, int] | None:
        """Fetch and parse one locus through the unchanged handle and window policy."""
        handle = self._open()
        if handle is None:
            return None

        start = max(0, position - 1 - self._flank)
        end = position + self._flank
        try:
            fetched = list(handle.fetch(contig, start, end))
        except (KeyError, ValueError) as error:
            logger.debug("BAM rescue skipped; %s:%s unavailable: %s", contig, position, error)
            return None
        self.fetches += 1

        records: list[_LocusRecord] = []
        for record in fetched:
            depth = minimum_kmer_depth(record)
            if record.cigartuples is None:
                records.append(_LocusRecord((), depth))
                continue
            edits = tuple(merge_edits(walk_cigar(record.reference_start, record.cigartuples, record.query_sequence)))
            records.append(_LocusRecord(edits, depth))
        return tuple(records), start, end

    @staticmethod
    def _compatibility_consensus(
        records: tuple[_LocusRecord, ...],
        start: int,
        end: int,
        contig: str,
        position: int,
    ) -> BamConsensus | None:
        """Apply the established unweighted exact-edit vote to parsed records."""
        votes: Counter[tuple[str, int, int, int, str]] = Counter()
        minimum_kmer_depths: dict[tuple[str, int, int, int, str], list[int | None]] = {}
        for record in records:
            for edit in record.edits:
                if edit.changes_length and start <= edit.start <= end:
                    edit_key = (edit.kind, edit.start, edit.ref_span, edit.inserted, edit.bases)
                    votes[edit_key] += 1
                    minimum_kmer_depths.setdefault(edit_key, []).append(record.minimum_kmer_depth)

        if not votes:
            return None

        ranked = votes.most_common()
        winning_edit, supporting_haplotype_records = ranked[0]
        kind, edit_start, ref_span, inserted, bases = winning_edit

        # A tie has no winner. `Counter.most_common` breaks one by insertion order,
        # which here is BAM order -- so whichever haplotype record happened to be
        # written first would decide the allele, and that allele can override the
        # VCF. Report no rescue instead of choosing arbitrarily.
        if len(ranked) > 1 and ranked[1][1] == supporting_haplotype_records:
            logger.debug(
                "BAM rescue declined at %s:%s; %d alleles tied on %d haplotype records",
                contig,
                position,
                sum(1 for _, count in ranked if count == supporting_haplotype_records),
                supporting_haplotype_records,
            )
            return None
        return BamConsensus(
            kind=kind,
            start=edit_start,
            ref_span=ref_span,
            inserted=inserted,
            bases=bases,
            supporting_haplotype_records=supporting_haplotype_records,
            fetched_haplotype_records=len(records),
            distinct_edit_count=len(votes),
            supporting_record_minimum_kmer_depths=tuple(minimum_kmer_depths[winning_edit]),
        )
