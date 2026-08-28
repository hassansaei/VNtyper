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


@dataclass(frozen=True)
class FastpMeasurement:
    """One validated fastp measurement shared by its display and icon decision."""

    value: Decimal | None
    display: str | None


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


def _validated_measured_fraction(value: object, key: str) -> Decimal | None:
    """Validate one output.json fastp fraction while preserving an absent measurement."""
    if value is None:
        return None
    message = f"Fastp output has invalid measured rate {key!r}: expected a finite numeric fraction from 0 to 1."
    if isinstance(value, bool) or not isinstance(value, (numbers.Real, Decimal)):
        logger.error(message)
        raise ValueError(message)
    try:
        fraction = Decimal(str(value))
    except (InvalidOperation, ValueError) as error:
        logger.error(message)
        raise ValueError(message) from error
    if not fraction.is_finite() or not Decimal(0) <= fraction <= Decimal(1):
        logger.error(message)
        raise ValueError(message)
    return fraction


def _validated_passed_filter_count(value: object, component: str) -> Decimal:
    """Validate one source count before deriving the passed-filter rate."""
    message = (
        f"Fastp output has invalid passed_filter_rate source count {component!r}: "
        "expected a finite non-negative numeric count."
    )
    if isinstance(value, bool) or not isinstance(value, (numbers.Real, Decimal)):
        logger.error(message)
        raise ValueError(message)
    try:
        count = Decimal(str(value))
    except (InvalidOperation, ValueError) as error:
        logger.error(message)
        raise ValueError(message) from error
    if not count.is_finite() or count < Decimal(0):
        logger.error(message)
        raise ValueError(message)
    return count


def calculate_passed_filter_rate(passed_filter_reads: object, total_reads: object) -> float | None:
    """Validate source counts and calculate the passed-filter rate.

    Args:
        passed_filter_reads: The count fastp reports as passing its filters.
        total_reads: The count fastp reports before filtering.

    Returns:
        The valid raw float fraction, or ``None`` when the valid total is zero.

    Raises:
        ValueError: If either count is malformed, or passed reads exceed total reads.
    """
    passed_count = _validated_passed_filter_count(passed_filter_reads, "passed_filter_reads")
    total_count = _validated_passed_filter_count(total_reads, "total_reads")
    if passed_count > total_count:
        message = (
            "Fastp output has invalid passed_filter_rate source count 'passed_filter_reads': "
            "cannot exceed 'total_reads'."
        )
        logger.error(message)
        raise ValueError(message)
    if total_count == Decimal(0):
        return None
    return float(passed_count / total_count)


def _cutoff(value: object, key: str) -> FastpCutoff:
    """Build one validated fastp decision/display pair."""
    fraction = _validated_fraction(value, key)
    displayed_fraction = _displayed_fraction(fraction)
    return FastpCutoff(value=displayed_fraction, label=_cutoff_label(displayed_fraction))


def _displayed_fraction(fraction: Decimal) -> Decimal:
    """Round a fraction into the report's two-decimal percentage domain."""
    percentage = fraction * Decimal(100)
    return percentage.quantize(Decimal("0.01"), rounding=ROUND_HALF_UP) / Decimal(100)


def _cutoff_label(fraction: Decimal) -> str:
    """Format a displayed fastp fraction as its concise cutoff label."""
    return f"{format((fraction * Decimal(100)).normalize(), 'f')}%"


def fastp_display_rate(rate: object) -> str | None:
    """Format one raw fastp metric with the icon's exact decision rounding.

    Args:
        rate: Raw fastp fraction, or ``None`` when the metric was not measured.

    Returns:
        The visible percentage string, preserving a missing value.
    """
    return build_fastp_measurement(rate, "fastp rate").display


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


def fastp_threshold_rate(rate: object) -> Decimal | None:
    """Return a fastp fraction rounded to the two-decimal percentage readers see.

    Args:
        rate: Raw fastp fraction, or ``None`` when the metric was not measured.

    Returns:
        The displayed rate as a fraction, preserving a missing value.
    """
    return build_fastp_measurement(rate, "fastp rate").value


def build_fastp_measurement(rate: object, key: str) -> FastpMeasurement:
    """Build one validated fastp value used for both its display and icon.

    Args:
        rate: Raw fraction from fastp ``output.json``, or ``None`` when unavailable.
        key: Report metric key used in an error message.

    Returns:
        One decimal decision value and matching display text, or a missing pair.

    Raises:
        ValueError: If a nonmissing rate is not a finite numeric 0..1 fraction.
    """
    fraction = _validated_measured_fraction(rate, key)
    if fraction is None:
        return FastpMeasurement(value=None, display=None)
    displayed_fraction = _displayed_fraction(fraction)
    percentage = format((displayed_fraction * Decimal(100)).normalize(), "f")
    if "." not in percentage:
        percentage = f"{percentage}.0"
    return FastpMeasurement(value=displayed_fraction, display=f"{percentage}%")
