"""Validated fastp report cutoffs and their paired display labels."""

from __future__ import annotations

import logging
import numbers
from collections.abc import Mapping
from dataclasses import dataclass
from decimal import ROUND_HALF_UP, Decimal, InvalidOperation

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class FastpCutoff:
    """One numeric fastp decision cutoff and the label derived from it."""

    value: Decimal
    label: str


@dataclass(frozen=True)
class FastpCutoffs:
    """The four configured fastp report cutoffs with their display labels."""

    duplication_rate: FastpCutoff
    q20_rate: FastpCutoff
    q30_rate: FastpCutoff
    passed_filter_rate: FastpCutoff


def _validated_fraction(value: object, key: str) -> Decimal:
    """Validate one configured fastp fraction in the report's decimal domain."""
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        message = f"Config thresholds has invalid fastp cutoff {key!r}: expected a numeric fraction from 0 to 1."
        logger.error(message)
        raise ValueError(message)
    try:
        fraction = Decimal(str(value))
    except (InvalidOperation, ValueError) as error:
        message = f"Config thresholds has invalid fastp cutoff {key!r}: expected a finite fraction from 0 to 1."
        logger.error(message)
        raise ValueError(message) from error
    if not fraction.is_finite() or not Decimal(0) <= fraction <= Decimal(1):
        message = f"Config thresholds has invalid fastp cutoff {key!r}: expected a finite fraction from 0 to 1."
        logger.error(message)
        raise ValueError(message)
    return fraction


def _cutoff(value: object, key: str) -> FastpCutoff:
    """Build one validated fastp decision/display pair."""
    fraction = _validated_fraction(value, key)
    displayed_fraction = _displayed_fraction(fraction)
    percentage = displayed_fraction * Decimal(100)
    return FastpCutoff(value=displayed_fraction, label=f"{format(percentage.normalize(), 'f')}%")


def _displayed_fraction(fraction: Decimal) -> Decimal:
    """Round a fraction into the report's two-decimal percentage domain."""
    percentage = fraction * Decimal(100)
    return percentage.quantize(Decimal("0.01"), rounding=ROUND_HALF_UP) / Decimal(100)


def build_fastp_cutoffs(thresholds: Mapping[str, object]) -> FastpCutoffs:
    """Validate configured fastp cutoffs and pair every value with its display label.

    Args:
        thresholds: Report configuration's threshold mapping.

    Returns:
        The four fastp cutoff representations used by both icon decisions and labels.

    Raises:
        ValueError: If a required cutoff is missing or is not a finite 0..1 fraction.
    """
    keys = ("duplication_rate", "q20_rate", "q30_rate", "passed_filter_reads_rate")
    values: dict[str, FastpCutoff] = {}
    for key in keys:
        try:
            values[key] = _cutoff(thresholds[key], key)
        except KeyError as error:
            message = f"Config thresholds is missing required fastp cutoff {key!r}."
            logger.error(message)
            raise ValueError(message) from error
    return FastpCutoffs(
        duplication_rate=values["duplication_rate"],
        q20_rate=values["q20_rate"],
        q30_rate=values["q30_rate"],
        passed_filter_rate=values["passed_filter_reads_rate"],
    )


def fastp_threshold_rate(rate: float | None) -> Decimal | None:
    """Return a fastp fraction rounded to the two-decimal percentage readers see.

    Args:
        rate: Raw fastp fraction, or ``None`` when the metric was not measured.

    Returns:
        The displayed rate as a fraction, preserving a missing value.
    """
    if rate is None:
        return None
    return _displayed_fraction(Decimal(str(rate)))
