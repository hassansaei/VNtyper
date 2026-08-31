"""Validated fastp report cutoffs and their paired display labels."""

from __future__ import annotations

import copy
import logging
import numbers
from collections.abc import Mapping
from dataclasses import dataclass
from decimal import ROUND_HALF_UP, Decimal, InvalidOperation, localcontext
from typing import Any, cast

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


class DecimalPreservingMapping(dict[str, Any]):
    """A legacy-valued JSON object with exact Decimal source values.

    The mapping exposes the values produced by ordinary ``json.loads`` so
    existing callers keep their dict/float contract. Exact values are sidecar
    provenance for the small number of report decisions that require it. If a
    caller mutates a value, the current value wins rather than allowing stale
    provenance to influence a later decision.
    """

    def __init__(self, raw: Mapping[str, Any], exact: Mapping[str, Any], *, recursive: bool = False) -> None:
        """Build a dict-compatible mapping from ordinary and exact JSON parses.

        Args:
            raw: Values parsed with the standard JSON numeric types.
            exact: Values parsed with ``parse_float=Decimal``.
            recursive: Whether nested object pairs should carry provenance too.
        """
        values: dict[str, Any] = {}
        for key, raw_value in raw.items():
            exact_value = exact.get(key, raw_value)
            if recursive and isinstance(raw_value, Mapping) and isinstance(exact_value, Mapping):
                values[key] = DecimalPreservingMapping(raw_value, exact_value, recursive=True)
            else:
                values[key] = raw_value
        super().__init__(values)
        self._source_values = dict(values)
        self._exact_values = dict(exact)

    def __setitem__(self, key: str, value: Any) -> None:
        """Store a public value and invalidate source provenance for that key."""
        super().__setitem__(key, value)
        if hasattr(self, "_source_values"):
            self._invalidate_provenance(key)

    def __delitem__(self, key: str) -> None:
        """Delete a public value and its exact source provenance."""
        super().__delitem__(key)
        self._invalidate_provenance(key)

    # Runtime ``dict.__ior__`` accepts the same mapping/iterable inputs as ``update``;
    # mypy incorrectly requires its signature to match narrower ``dict.__or__`` too.
    def __ior__(self, values: object) -> DecimalPreservingMapping:  # type: ignore[override,misc]
        """Merge public values through the provenance-invalidating update path."""
        self.update(values)
        return self

    def update(
        self,
        values: object = (),
        /,
        **kwargs: Any,
    ) -> None:
        """Update public values while invalidating their source provenance."""
        pending: dict[str, Any] = {}
        pending.update(cast(Any, values))
        pending.update(kwargs)
        for key, value in pending.items():
            self[key] = value

    def setdefault(self, key: str, default: Any = None) -> Any:
        """Insert a default through the invalidating path only when absent."""
        if key in self:
            return self[key]
        self[key] = default
        return default

    def pop(self, key: str, *default: Any) -> Any:
        """Remove a public value and its exact source provenance."""
        value = super().pop(key, *default)
        self._invalidate_provenance(key)
        return value

    def popitem(self) -> tuple[str, Any]:
        """Remove the last public item and its exact source provenance."""
        key, value = super().popitem()
        self._invalidate_provenance(key)
        return key, value

    def clear(self) -> None:
        """Remove every public value and all exact source provenance."""
        super().clear()
        self._source_values.clear()
        self._exact_values.clear()

    def __deepcopy__(self, memo: dict[int, Any]) -> DecimalPreservingMapping:
        """Copy public values and only the provenance that remains valid."""
        copied = type(self)(copy.deepcopy(dict(self), memo), copy.deepcopy(self.exact, memo))
        memo[id(self)] = copied
        return copied

    @property
    def exact(self) -> Mapping[str, Any]:
        """Return current values with still-valid exact provenance applied."""
        values: dict[str, Any] = {}
        for key, value in self.items():
            values[key] = value.exact if isinstance(value, DecimalPreservingMapping) else self.exact_value(key)
        return values

    def exact_value(self, key: str) -> object:
        """Return exact source provenance unless the public value was changed."""
        current = self[key]
        if key not in self._source_values:
            return current
        source = self._source_values[key]
        if type(current) is type(source) and current == source:
            return self._exact_values.get(key, current)
        return current

    def _invalidate_provenance(self, key: str) -> None:
        """Discard exact source metadata after one public-key mutation."""
        self._source_values.pop(key, None)
        self._exact_values.pop(key, None)


class FastpJsonPayload(DecimalPreservingMapping):
    """A dict-compatible fastp JSON document with exact-decimal provenance."""

    def __init__(self, raw: Mapping[str, Any], exact: Mapping[str, Any]) -> None:
        """Recursively retain exact rates while exposing legacy mapping values."""
        super().__init__(raw, exact, recursive=True)


def preserve_exact_fastp_thresholds(raw_config: object, exact_config: object) -> object:
    """Attach exact provenance only to the report's fastp threshold mapping.

    Args:
        raw_config: Complete config parsed with ordinary JSON numeric types.
        exact_config: The same JSON parsed with ``parse_float=Decimal``.

    Returns:
        The original top-level config object with a dict-compatible threshold
        mapping when both parses contain one.
    """
    if not isinstance(raw_config, dict) or not isinstance(exact_config, Mapping):
        return raw_config
    raw_thresholds = raw_config.get("thresholds")
    exact_thresholds = exact_config.get("thresholds")
    if isinstance(raw_thresholds, Mapping) and isinstance(exact_thresholds, Mapping):
        raw_config["thresholds"] = DecimalPreservingMapping(raw_thresholds, exact_thresholds)
    return raw_config


def exact_json_value(mapping: Mapping[str, object], key: str) -> object:
    """Return an exact JSON value when compatible provenance is available."""
    if isinstance(mapping, DecimalPreservingMapping):
        return mapping.exact_value(key)
    return mapping[key]


def _validated_fraction(value: object, key: str) -> Decimal:
    """Validate one configured fastp fraction in the report's decimal domain."""
    if isinstance(value, bool) or not isinstance(value, (numbers.Real, Decimal)):
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


def validated_fastp_fraction(value: object, key: str) -> Decimal | None:
    """Validate one output.json fastp fraction while preserving absence.

    Args:
        value: Parsed fastp rate, or ``None`` when that optional metric is absent.
        key: Stable metric key used in operator-facing failures.

    Returns:
        The exact finite Decimal fraction, or ``None``.

    Raises:
        ValueError: If a nonmissing value is not a finite numeric 0..1 fraction.
    """
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
        "expected a finite non-negative integer count."
    )
    if isinstance(value, bool) or not isinstance(value, (numbers.Real, Decimal)):
        logger.error(message)
        raise ValueError(message)
    try:
        count = Decimal(str(value))
    except (InvalidOperation, ValueError) as error:
        logger.error(message)
        raise ValueError(message) from error
    if not count.is_finite() or count < Decimal(0) or count != count.to_integral_value():
        logger.error(message)
        raise ValueError(message)
    return count


def calculate_passed_filter_rate(passed_filter_reads: object, total_reads: object) -> Decimal | None:
    """Validate source counts and calculate the passed-filter rate.

    Args:
        passed_filter_reads: The count fastp reports as passing its filters.
        total_reads: The count fastp reports before filtering.

    Returns:
        The exact Decimal fraction, or ``None`` when the valid total is zero.

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
    with localcontext() as context:
        context.prec = max(context.prec, total_count.adjusted() + 11)
        return passed_count / total_count


def calculate_passed_filter_rate_from_sources(
    before_filtering: Mapping[str, object], filtering_result: Mapping[str, object]
) -> Decimal | None:
    """Extract required fastp counts and calculate the passed-filter rate.

    Args:
        before_filtering: The fastp ``summary.before_filtering`` mapping.
        filtering_result: The fastp ``filtering_result`` mapping.

    Returns:
        The exact Decimal fraction, or ``None`` when both counts are zero.

    Raises:
        ValueError: If either required source key is missing or either count is invalid.
    """
    sources = (
        (before_filtering, "total_reads", "summary.before_filtering.total_reads"),
        (filtering_result, "passed_filter_reads", "filtering_result.passed_filter_reads"),
    )
    values: dict[str, object] = {}
    for source, key, path in sources:
        try:
            values[key] = source[key]
        except KeyError as error:
            message = f"Fastp output is missing required passed_filter_rate source key {path!r}."
            logger.error(message)
            raise ValueError(message) from error
    return calculate_passed_filter_rate(values["passed_filter_reads"], values["total_reads"])


def validated_fastp_mapping(value: object, path: str) -> Mapping[str, Any]:
    """Validate one object boundary in parsed fastp JSON.

    Args:
        value: The parsed value at the boundary.
        path: The schema path reported to the operator.

    Returns:
        The mapping after validation.

    Raises:
        ValueError: If the parsed value is not a JSON object.
    """
    if not isinstance(value, Mapping):
        message = f"Fastp output has invalid object at {path!r}: expected a dictionary."
        logger.error(message)
        raise ValueError(message)
    return value


def _cutoff(value: object, key: str) -> FastpCutoff:
    """Build one validated fastp decision/display pair."""
    fraction = _validated_fraction(value, key)
    displayed_fraction = _displayed_fraction(fraction)
    return FastpCutoff(value=displayed_fraction, label=_cutoff_label(displayed_fraction))


def _displayed_fraction(fraction: Decimal) -> Decimal:
    """Round a fraction into the report's two-decimal percentage domain."""
    displayed_fraction = fraction.quantize(Decimal("0.0001"), rounding=ROUND_HALF_UP)
    return Decimal(0) if displayed_fraction.is_zero() else displayed_fraction


def _cutoff_label(fraction: Decimal) -> str:
    """Format a displayed fastp fraction as its concise cutoff label."""
    return f"{format((fraction * Decimal(100)).normalize(), 'f')}%"


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
            values[key] = _cutoff(exact_json_value(thresholds, key), key)
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
    fraction = validated_fastp_fraction(rate, key)
    if fraction is None:
        return FastpMeasurement(value=None, display=None)
    displayed_fraction = _displayed_fraction(fraction)
    percentage = format((displayed_fraction * Decimal(100)).normalize(), "f")
    if "." not in percentage:
        percentage = f"{percentage}.0"
    return FastpMeasurement(value=displayed_fraction, display=f"{percentage}%")
