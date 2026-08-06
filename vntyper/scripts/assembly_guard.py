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
STATUS_CONFLICT = "conflict"
STATUS_UNDETERMINED = "undetermined"

#: Value used for any build or convention that could not be determined.
UNKNOWN = "unknown"


@dataclass(frozen=True)
class AssemblyVerdict:
    """The outcome of reconciling a declared assembly against a header.

    Attributes:
        status: One of :data:`STATUS_AGREE`, :data:`STATUS_MISMATCH`,
            :data:`STATUS_CONFLICT`, :data:`STATUS_UNDETERMINED`. `undetermined`
            means the question could not be answered -- it is neither a pass nor
            a failure. `conflict` means the header answered it *twice, and
            contradictorily*: two aliases of chromosome 1 carrying two different
            recognised build lengths. That is a defect in the input rather than
            an unanswered question, so callers treat it as they treat
            `mismatch`.
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

    A header that names chromosome 1 more than once with *conflicting* lengths
    -- a hybrid carrying both `1` and `chr1` at two lengths this module
    recognises as different builds -- is `conflict`. Contig order is still never
    consulted: no ordering of a self-contradictory header makes it consistent,
    and whichever alias downstream naming resolution picks, the declared build's
    MUC1 coordinates land on a contig of the other build's length. That is the
    wrong slice and a plausible false negative, so the evidence is reported as
    contradictory rather than as absent, and `message` names both contigs and
    both lengths.

    An **unrecognised** chr1 length is a different case and stays
    `undetermined`: a patched or non-human reference can declare a chromosome 1
    this module has no entry for, and that is an unanswered question, not a
    contradiction. Only mutually contradictory *recognised* evidence conflicts.

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

    Returns:
        AssemblyVerdict: The verdict. `status` is `agree` when both builds are
        known and equal, `mismatch` when both are known and differ, `conflict`
        when the header's own chr1 entries name two different recognised builds,
        and `undetermined` when the build could not be determined at all.

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

    recognised = _recognised_builds(chr1_lengths)
    if len({build for _, build in recognised}) > 1:
        # The header contradicts itself. Decided before the declared/detected
        # comparison because no declaration can reconcile it: there is no build
        # such a header agrees with, and `undetermined` would fail open onto
        # whichever alias downstream naming resolution happens to pick.
        conflict_message = _conflict_message(declared, usable, recognised)
        logger.debug(conflict_message)
        return AssemblyVerdict(
            status=STATUS_CONFLICT,
            declared=declared,
            declared_coordinate_system=declared_build,
            coordinate_system=UNKNOWN,
            naming_convention=naming_convention,
            chr1_length=None,
            message=conflict_message,
        )

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
        # Reaching here means at most one of these lengths is a recognised build:
        # two recognised ones are a `conflict` and never arrive. Ordering is not
        # evidence in either branch, so it must not reach the message either.
        return (
            "the header names chromosome 1 more than once with differing lengths ("
            + _describe_chr1_contigs(usable, chr1_lengths)
            + "), and no single one of them identifies a build"
        )
    return f"chr1 is {chr1_lengths[0]:,} bp, which matches neither GRCh37 nor GRCh38"


def _recognised_builds(chr1_lengths: list[int]) -> list[tuple[int, str]]:
    """Pair each chr1 length with the build it identifies, dropping the unrecognised.

    Each length is offered to detection on its own, synthesised into a single
    contig, so that neither header order nor the presence of a second chromosome 1
    can influence what any one length resolves to.

    Args:
        chr1_lengths (list[int]): Result of :func:`_chr1_lengths`.

    Returns:
        list[tuple[int, str]]: `(length, build)` for every length the length ->
        build table knows, in the order given. Lengths it does not know are
        omitted, which is what keeps an unrecognised chromosome 1 out of the
        conflict decision and on the `undetermined` path.
    """
    recognised: list[tuple[int, str]] = []
    for length in chr1_lengths:
        build = detect_assembly_from_chr1_length([{"name": "chr1", "length": length}])
        if build:
            recognised.append((length, build))
    return recognised


def _describe_chr1_contigs(usable: list[dict[str, Any]], lengths: list[int]) -> str:
    """Name the chromosome 1 contigs that carry the given lengths.

    Args:
        usable (list[dict[str, Any]]): Contigs that survived filtering.
        lengths (list[int]): The chr1 lengths to describe.

    Returns:
        str: A comma-separated clause such as `"1 is 248,956,422 bp, chr1 is
        249,250,621 bp"`, sorted by contig name so that header order cannot
        reach the message.
    """
    from vntyper.scripts.chromosome_utils import is_chr1_name

    named = sorted(
        {
            (contig["name"], contig["length"])
            for contig in usable
            if is_chr1_name(contig["name"]) and contig.get("length") in lengths
        }
    )
    return ", ".join(f"{name} is {length:,} bp" for name, length in named)


def _conflict_message(declared: str, usable: list[dict[str, Any]], recognised: list[tuple[int, str]]) -> str:
    """Compose the message for a self-contradictory header.

    This message becomes the exception the user reads, so it names both contigs,
    both lengths and both builds -- enough to see which alias is wrong without
    looking the length table up by hand.

    Args:
        declared (str): The assembly the caller declared, verbatim.
        usable (list[dict[str, Any]]): Contigs that survived filtering.
        recognised (list[tuple[int, str]]): Result of :func:`_recognised_builds`.

    Returns:
        str: A complete sentence suitable as both a log line and an exception
        message.
    """
    builds = sorted({build for _, build in recognised})
    described = _describe_chr1_contigs(usable, [length for length, _ in recognised])
    return (
        f"The alignment header contradicts itself: it names chromosome 1 more than once with conflicting "
        f"lengths ({described}), identifying {' and '.join(builds)} at the same time. Declared reference "
        f"assembly '{declared}' cannot be reconciled with a header that describes two builds: the MUC1 region "
        f"would be extracted from whichever alias contig resolution selects, and for one of these two that is "
        f"the wrong build's coordinates, which yields a false negative. Supply an alignment whose header names "
        f"chromosome 1 once."
    )
