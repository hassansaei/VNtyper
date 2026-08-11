"""Pure CRAM reference-candidate policy and coverage decisions."""

from __future__ import annotations

import logging
from collections.abc import Iterable, Set
from dataclasses import dataclass
from typing import Any

from vntyper.scripts.reference_registry import reference_keys

logger = logging.getLogger(__name__)

_DEFAULT_REFERENCE_ORDER = ("cli", "config_cram_reference", "config_bwa_reference", "htslib_resolved")
_EXPLICIT_REFERENCE_SOURCES = {"cli", "config_cram_reference", "config_bwa_reference"}
_TERMINAL_REFERENCE_SOURCE = "htslib_resolved"


@dataclass(frozen=True)
class ResolvedReference:
    """Which config key supplied a reference, and whether it was a fallback.

    Attributes:
        key: The `reference_data` key that was present.
        value: Its value. `None` means present-and-disabled, not missing.
        is_fallback: True when the UCSC-family key stood in for an absent physical
            key, which means the run uses UCSC sequence for a non-UCSC request.
    """

    key: str
    value: str | None
    is_fallback: bool


def resolve_from_mapping(kind: str, reference_assembly: str, mapping: dict) -> ResolvedReference | None:
    """Resolve a reference from a config mapping by key membership.

    Args:
        kind: One of `reference_registry.REFERENCE_KINDS`.
        reference_assembly: Supported assembly label.
        mapping: A dict of config keys to values, e.g. `config["reference_data"]`.

    Returns:
        The first key that is *present* in the mapping, or None when none are.
        Presence wins over truthiness so an explicit null stays authoritative.
        `is_fallback` is True only when the *last* (UCSC-family) key resolved and
        it was not the only key available - not merely when the resolved key is
        not the first one. The middle, physical-id tier that `hg19_ncbi`/`hg38_ncbi`
        share with `GRCh37`/`GRCh38` is that label's own physical file, not a UCSC
        stand-in, so matching it there must not be flagged as a fallback.

    Raises:
        ValueError: If the kind or the assembly is unknown.
    """
    keys = reference_keys(kind, reference_assembly)
    for index, key in enumerate(keys):
        if key in mapping:
            is_fallback = index == len(keys) - 1 and len(keys) > 1
            return ResolvedReference(key=key, value=mapping[key], is_fallback=is_fallback)
    return None


def configured_reference_candidates(
    config: dict,
    reference_assembly: str,
) -> tuple[tuple[str, Any], ...]:
    """Map configured CRAM reference sources without validating probe policy.

    Exact assembly-label keys take precedence even when their value is null;
    otherwise the matching UCSC-family key supplies the fallback.

    Args:
        config: Pipeline configuration containing reference paths.
        reference_assembly: Supported assembly label used for exact and family
            lookups.

    Returns:
        The raw configured CRAM and BWA reference-source values in stable order.
        Callers that interpret them as paths must require strings.
    """
    reference_data = config.get("reference_data", {})
    values: list[tuple[str, Any]] = []
    for source, kind in (("config_cram_reference", "cram"), ("config_bwa_reference", "bwa")):
        resolved = resolve_from_mapping(kind, reference_assembly, reference_data)
        values.append((source, resolved.value if resolved is not None else None))
    return tuple(values)


def ordered_reference_candidates(
    config: dict,
    reference_assembly: str,
    reference_fasta: str | None,
) -> tuple[tuple[str, str | None], ...]:
    """Validate reference policy and map explicit sources to their paths.

    Args:
        config: Pipeline configuration. Missing candidate policy uses the shipped
            explicit-first order.
        reference_assembly: Supported assembly label used first as an exact
            configured-key suffix, then by coordinate system for the UCSC
            family fallback.
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

    values = {"cli": reference_fasta, **dict(configured_reference_candidates(config, reference_assembly))}
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
