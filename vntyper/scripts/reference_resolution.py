"""Pure CRAM reference-candidate policy and coverage decisions."""

from __future__ import annotations

import logging
from collections.abc import Iterable, Set

logger = logging.getLogger(__name__)

_DEFAULT_REFERENCE_ORDER = ("cli", "config_cram_reference", "config_bwa_reference", "htslib_resolved")
_EXPLICIT_REFERENCE_SOURCES = {"cli", "config_cram_reference", "config_bwa_reference"}
_TERMINAL_REFERENCE_SOURCE = "htslib_resolved"


def ordered_reference_candidates(
    config: dict,
    reference_assembly: str,
    reference_fasta: str | None,
) -> tuple[tuple[str, str | None], ...]:
    """Validate reference policy and map explicit sources to their paths.

    Args:
        config: Pipeline configuration. Missing candidate policy uses the shipped
            explicit-first order.
        reference_assembly: Assembly suffix used for configured reference keys.
        reference_fasta: Explicit CLI or web reference path, when supplied.

    Returns:
        Explicit ``(source, path)`` candidates in configured order. The required
        terminal htslib source is intentionally omitted because it has no path.

    Raises:
        ValueError: If the policy is not a list, contains unknown or duplicate
            sources, or lacks exactly one terminal htslib source at the end.
    """
    order = config.get("cram", {}).get("reference_candidate_order", list(_DEFAULT_REFERENCE_ORDER))
    if not isinstance(order, list):
        message = "cram.reference_candidate_order must be a list"
        logger.error(message)
        raise ValueError(message)
    if order.count(_TERMINAL_REFERENCE_SOURCE) != 1:
        message = "reference candidate order must contain exactly one terminal htslib_resolved"
        logger.error(message)
        raise ValueError(message)
    if order[-1] != _TERMINAL_REFERENCE_SOURCE:
        message = "reference candidate order must end with terminal htslib_resolved"
        logger.error(message)
        raise ValueError(message)
    for source in order:
        if not isinstance(source, str) or source not in _EXPLICIT_REFERENCE_SOURCES | {_TERMINAL_REFERENCE_SOURCE}:
            message = f"unknown reference candidate source: {source}"
            logger.error(message)
            raise ValueError(message)
    explicit_order = order[:-1]
    if len(explicit_order) != len(set(explicit_order)):
        message = "duplicate explicit reference candidate source"
        logger.error(message)
        raise ValueError(message)

    reference_data = config.get("reference_data", {})
    values = {
        "cli": reference_fasta,
        "config_cram_reference": reference_data.get(f"cram_reference_{reference_assembly}"),
        "config_bwa_reference": reference_data.get(f"bwa_reference_{reference_assembly}"),
    }
    return tuple((source, values[source]) for source in explicit_order)


def uncovered_reference_contigs(
    header_contigs: Iterable[str],
    reference_contigs: Set[str] | None,
) -> tuple[str, ...] | None:
    """Decide which alignment contigs a reference is known not to cover.

    Args:
        header_contigs: Contigs declared by the alignment header, in report order.
        reference_contigs: Contigs read from a FASTA index, or ``None`` when
            coverage evidence is unavailable.

    Returns:
        Missing contigs in header order, an empty tuple for known complete
        coverage, or ``None`` when coverage cannot be assessed.
    """
    if reference_contigs is None:
        return None
    return tuple(contig for contig in header_contigs if contig not in reference_contigs)
