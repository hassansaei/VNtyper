"""
assembly_guard.py

Reconcile the reference assembly a caller *declares* with the one the alignment
header actually *describes*.

`--reference-assembly` is otherwise taken on faith. The MUC1 VNTR sits at
chr1:155,158,000-155,163,000 in GRCh37 and chr1:155,184,000-155,194,000 in
GRCh38 -- roughly 30 kb apart. Declaring the wrong build slices a region that
does not contain the VNTR, so Kestrel sees no supporting reads and the pipeline
reports a confident negative. That is a silent false negative on patient data,
which is the failure class this module exists to make visible.

This module is deliberately **pure and non-fatal**: it computes a verdict and
never raises, never logs an error, and never decides what should happen next.
Turning a mismatch into a hard failure is a behaviour change and belongs to the
caller, so that the verdict and the enforcement can be reviewed and reverted
independently.

Functions:
    reconcile_assembly: Compare a declared assembly against parsed contigs
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

from vntyper.scripts.chromosome_utils import (
    detect_assembly_from_chr1_length,
    detect_naming_convention,
)
from vntyper.scripts.reference_registry import get_coordinate_system

logger = logging.getLogger(__name__)

# Verdict statuses. Callers should import these rather than spell the strings.
STATUS_AGREE = "agree"
STATUS_MISMATCH = "mismatch"
STATUS_UNDETERMINED = "undetermined"

#: Value used for any build or convention that could not be determined.
UNKNOWN = "unknown"


@dataclass(frozen=True)
class AssemblyVerdict:
    """The outcome of reconciling a declared assembly against a header.

    Attributes:
        status: One of :data:`STATUS_AGREE`, :data:`STATUS_MISMATCH`,
            :data:`STATUS_UNDETERMINED`. `undetermined` means the question could
            not be answered -- it is neither a pass nor a failure.
        declared: The assembly name the caller declared, verbatim.
        declared_coordinate_system: The *build* `declared` resolves to --
            "GRCh37", "GRCh38", or :data:`UNKNOWN` if the name is not in the
            reference registry.
        coordinate_system: The *build* the contigs describe -- "GRCh37",
            "GRCh38", or :data:`UNKNOWN`. This is the detected value; compare it
            against `declared_coordinate_system`, not against a naming
            convention. The two are different questions.
        naming_convention: How the contigs are *named* -- "ucsc", "ensembl",
            "ncbi" or :data:`UNKNOWN`. Reported for diagnostics only; it never
            affects `status`, because a header can name its contigs
            inconsistently and still identify its build unambiguously.
        chr1_length: The chr1 length the verdict was reasoned from, or None if
            no usable chr1 was present *or* several contigs named chromosome 1
            with conflicting lengths. This is the evidence behind `status`; when
            it is None the reason is spelled out in `message`.
        message: A complete, human-readable sentence suitable both as a log
            line and as an exception message.
    """

    status: str
    declared: str
    declared_coordinate_system: str
    coordinate_system: str
    naming_convention: str
    chr1_length: int | None
    message: str


def _declared_coordinate_system(declared: str) -> str:
    """Resolve a declared assembly name to its build.

    Args:
        declared (str): Assembly name as the caller supplied it.

    Returns:
        str: "GRCh37", "GRCh38", or :data:`UNKNOWN`.
    """
    try:
        return get_coordinate_system(declared)
    except (ValueError, TypeError):
        logger.debug(f"Declared assembly '{declared}' is not in the reference registry")
        return UNKNOWN


def _usable_contigs(contigs: list[dict]) -> list[dict[str, Any]]:
    """Keep only entries that are dicts with a string name.

    The declared input shape is what `fastq_bam_processing.parse_contigs_from_header`
    returns, but this guard runs on whatever a real header happened to hold, so
    a malformed entry must be skipped rather than crash the pipeline.

    Args:
        contigs (list[dict]): Parsed contigs, possibly containing junk.

    Returns:
        list[dict[str, Any]]: The entries that can be reasoned about.
    """
    return [c for c in contigs if isinstance(c, dict) and isinstance(c.get("name"), str)]


def _chr1_lengths(contigs: list[dict[str, Any]]) -> list[int]:
    """Return every *distinct* integer chr1 length the contigs declare.

    A header should name chromosome 1 once, but nothing enforces that: a hybrid
    header can carry both `1` and `chr1`, and the two aliases can carry
    *different* lengths. Taking the first match would then let the order the
    contigs happen to appear in decide the build -- so
    `[{"1": 248956422}, {"chr1": 249250621}]` would read as GRCh38 and the same
    header written the other way round as GRCh37. Collect them all and let the
    caller decide deliberately.

    Non-integer and missing lengths are skipped rather than ending the search,
    so one malformed alias cannot mask a well-formed one.

    Args:
        contigs (list[dict[str, Any]]): Contigs already filtered by
            :func:`_usable_contigs`.

    Returns:
        list[int]: Distinct lengths, in the order the header declared them.
        Empty when no contig names chromosome 1 with an integer length.
    """
    from vntyper.scripts.chromosome_utils import is_chr1_name

    lengths: list[int] = []
    for contig in contigs:
        if not is_chr1_name(contig["name"]):
            continue
        length = contig.get("length")
        if isinstance(length, int) and not isinstance(length, bool) and length not in lengths:
            lengths.append(length)
    return lengths


def reconcile_assembly(declared: str, contigs: list[dict]) -> AssemblyVerdict:
    """Compare a declared reference assembly against a BAM/CRAM header's contigs.

    The chr1 length is the sole decider of the build: it is unique per build and
    unaffected by the ~170 alternate contigs a modern reference carries. The
    contig *naming* convention is reported alongside but never used to decide,
    because a header can name its contigs inconsistently and still identify its
    build unambiguously.

    This function never raises and never logs above INFO. It returns a verdict;
    the caller decides whether a mismatch is fatal.

    Args:
        declared (str): The assembly the caller declared, e.g. from
            `--reference-assembly`. Any name the reference registry knows is
            accepted ("hg19", "GRCh37", "hg38_ensembl", ...).
        contigs (list[dict]): Parsed contigs as
            `fastq_bam_processing.parse_contigs_from_header` returns them:
            `[{"name": "chr1", "length": 248956422}, ...]`. Entries that are not
            usable are skipped.

    A header that names chromosome 1 more than once with *conflicting* lengths
    -- a hybrid carrying both `1` and `chr1` from different builds -- is
    `undetermined`, never `mismatch`. Deciding from whichever alias came first
    would let contig order reject a usable input; refusing to decide keeps the
    guard from failing a run on the strength of an abnormal header, and the
    conflict is named in `message`.

    Returns:
        AssemblyVerdict: The verdict. `status` is `agree` when both builds are
        known and equal, `mismatch` when both are known and differ, and
        `undetermined` when either could not be determined -- including when
        several chr1 entries disagree.

    Examples:
        >>> v = reconcile_assembly("hg19", [{"name": "chr1", "length": 248956422}])
        >>> v.status
        'mismatch'
        >>> v.declared_coordinate_system, v.coordinate_system
        ('GRCh37', 'GRCh38')
    """
    usable = _usable_contigs(contigs)
    declared_build = _declared_coordinate_system(declared)
    naming_convention = detect_naming_convention([c["name"] for c in usable]) if usable else UNKNOWN
    chr1_lengths = _chr1_lengths(usable)

    # Exactly one chr1 length is the normal case and the only one worth deciding
    # from. Detection is fed a single synthesised contig rather than the raw list
    # so that contig order can never influence the verdict; the length -> build
    # table stays in chromosome_utils, which remains the one place that knows it.
    chr1_length = chr1_lengths[0] if len(chr1_lengths) == 1 else None
    if chr1_length is None:
        detected_build = UNKNOWN
    else:
        detected_build = detect_assembly_from_chr1_length([{"name": "chr1", "length": chr1_length}]) or UNKNOWN

    if declared_build != UNKNOWN and detected_build != UNKNOWN:
        if declared_build == detected_build:
            status = STATUS_AGREE
            message = (
                f"Declared reference assembly '{declared}' ({declared_build}) matches the alignment header "
                f"(chr1 length {chr1_length:,} bp)."
            )
        else:
            status = STATUS_MISMATCH
            message = (
                f"Declared reference assembly '{declared}' uses {declared_build} coordinates, but the alignment "
                f"header describes {detected_build} (chr1 length {chr1_length:,} bp). Extracting {declared_build} "
                f"coordinates from a {detected_build} alignment targets the wrong region and yields a false "
                f"negative. Re-run with --reference-assembly {detected_build}, or supply the matching input."
            )
    else:
        status = STATUS_UNDETERMINED
        reason = _undetermined_reason(declared_build, usable, chr1_lengths)
        message = (
            f"Could not verify the declared reference assembly '{declared}': {reason}. "
            f"The declared assembly is being used unchecked."
        )

    # Deliberately DEBUG for every outcome: this function computes a verdict, it
    # does not report one. The caller logs at the level its policy calls for.
    logger.debug(message)
    return AssemblyVerdict(
        status=status,
        declared=declared,
        declared_coordinate_system=declared_build,
        coordinate_system=detected_build,
        naming_convention=naming_convention,
        chr1_length=chr1_length,
        message=message,
    )


def _undetermined_reason(declared_build: str, usable: list[dict[str, Any]], chr1_lengths: list[int]) -> str:
    """Explain why no comparison was possible.

    Args:
        declared_build (str): Result of :func:`_declared_coordinate_system`.
        usable (list[dict[str, Any]]): Contigs that survived filtering.
        chr1_lengths (list[int]): Result of :func:`_chr1_lengths`.

    Returns:
        str: A clause naming the single reason, for embedding in a message.
    """
    if declared_build == UNKNOWN:
        return "it is not a name the reference registry knows"
    if not usable:
        return "the header contained no usable contigs"
    if not chr1_lengths:
        return "the header has no chr1 with an integer length"
    if len(chr1_lengths) > 1:
        # Sorted by contig name, not header order: the whole point of this branch
        # is that ordering is not evidence, so it must not reach the message either.
        from vntyper.scripts.chromosome_utils import is_chr1_name

        named = sorted(
            {
                (contig["name"], contig["length"])
                for contig in usable
                if is_chr1_name(contig["name"]) and contig.get("length") in chr1_lengths
            }
        )
        return (
            "the header names chromosome 1 more than once with conflicting lengths ("
            + ", ".join(f"{name} is {length:,} bp" for name, length in named)
            + "), so no single length identifies the build"
        )
    return f"chr1 is {chr1_lengths[0]:,} bp, which matches neither GRCh37 nor GRCh38"
