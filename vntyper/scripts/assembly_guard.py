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
            no usable chr1 was present. This is the evidence behind `status`.
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


def _chr1_length(contigs: list[dict[str, Any]]) -> int | None:
    """Return the integer chr1 length present in the contigs, if any.

    Args:
        contigs (list[dict[str, Any]]): Contigs already filtered by
            :func:`_usable_contigs`.

    Returns:
        int | None: The length, or None if chr1 is absent or its length is not
        an integer.
    """
    from vntyper.scripts.chromosome_utils import is_chr1_name

    for contig in contigs:
        if is_chr1_name(contig["name"]):
            length = contig.get("length")
            return length if isinstance(length, int) and not isinstance(length, bool) else None
    return None


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

    Returns:
        AssemblyVerdict: The verdict. `status` is `agree` when both builds are
        known and equal, `mismatch` when both are known and differ, and
        `undetermined` when either could not be determined.

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
    chr1_length = _chr1_length(usable)

    detected_build = detect_assembly_from_chr1_length(usable) or UNKNOWN

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
        reason = _undetermined_reason(declared_build, usable, chr1_length)
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


def _undetermined_reason(declared_build: str, usable: list[dict[str, Any]], chr1_length: int | None) -> str:
    """Explain why no comparison was possible.

    Args:
        declared_build (str): Result of :func:`_declared_coordinate_system`.
        usable (list[dict[str, Any]]): Contigs that survived filtering.
        chr1_length (int | None): Result of :func:`_chr1_length`.

    Returns:
        str: A clause naming the single reason, for embedding in a message.
    """
    if declared_build == UNKNOWN:
        return "it is not a name the reference registry knows"
    if not usable:
        return "the header contained no usable contigs"
    if chr1_length is None:
        return "the header has no chr1 with an integer length"
    return f"chr1 is {chr1_length:,} bp, which matches neither GRCh37 nor GRCh38"
