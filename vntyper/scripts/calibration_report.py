"""Strict data model and static renderer for calibration evidence reports."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType

from jinja2 import Environment, FileSystemLoader, StrictUndefined, select_autoescape

import vntyper

_PHASES = frozenset({"fitted", "validation", "held-out"})
_PROVENANCE_FIELDS = {
    "software_versions",
    "reference_versions",
    "sample_composition",
    "assays",
    "depths",
    "read_lengths",
    "independent_array_size",
    "mutation_classes",
    "manifest_hashes",
    "access_attempts",
    "boundary_coverage",
    "seeds",
}
_STATISTICS_FIELDS = {"intervals", "roc_rows", "pr_rows", "joint_surface_rows"}


@dataclass(frozen=True)
class CalibrationReport:
    """Complete static report context with no network or script dependency."""

    phase: str
    profile_sha256: str
    protocol_sha256: str
    evidence_sha256: str
    objective: str
    tier_metrics: tuple[Mapping[str, object], ...]
    abstentions: tuple[Mapping[str, object], ...]
    provenance: Mapping[str, tuple[str, ...]]
    statistics: Mapping[str, tuple[str, ...]]
    limitations: tuple[str, ...]


def decode_calibration_report(value: object) -> CalibrationReport:
    """Decode the closed version-1 report context.

    Args:
        value: Parsed JSON-compatible report data.

    Returns:
        A validated immutable report context.

    Raises:
        ValueError: If fields, hashes, phase, objective, rows, or limitations differ.
    """
    fields = {
        "schema_version",
        "phase",
        "profile_sha256",
        "protocol_sha256",
        "evidence_sha256",
        "objective",
        "tier_metrics",
        "abstentions",
        "provenance",
        "statistics",
        "limitations",
    }
    root = _exact(value, fields, "calibration report")
    if root["schema_version"] != "calibration-report-v1":
        raise ValueError("calibration report schema version must be calibration-report-v1")
    phase = root["phase"]
    if phase not in _PHASES:
        raise ValueError(f"unsupported calibration report phase: {phase!r}")
    objective = root["objective"]
    if objective != "lexicographic-safety-v1":
        raise ValueError("calibration report objective must be lexicographic-safety-v1")
    tiers = tuple(_tier_row(row) for row in _nonempty_sequence(root["tier_metrics"], "tier metrics"))
    if len({row["tier"] for row in tiers}) != len(tiers):
        raise ValueError("calibration report tier rows must be unique")
    abstentions = tuple(_abstention_row(row) for row in _sequence(root["abstentions"], "abstentions"))
    provenance = _string_sections(root["provenance"], _PROVENANCE_FIELDS, "provenance")
    statistics = _string_sections(root["statistics"], _STATISTICS_FIELDS, "statistics")
    limitations = _strings(root["limitations"], "limitations")
    if not limitations:
        raise ValueError("calibration report limitations must not be empty")
    return CalibrationReport(
        str(phase),
        _digest(root["profile_sha256"], "profile"),
        _digest(root["protocol_sha256"], "protocol"),
        _digest(root["evidence_sha256"], "evidence"),
        str(objective),
        tiers,
        abstentions,
        MappingProxyType(provenance),
        MappingProxyType(statistics),
        limitations,
    )


def render_calibration_report(report: CalibrationReport) -> str:
    """Render one deterministic, escaped, script-free HTML document."""
    if not isinstance(report, CalibrationReport):
        raise ValueError("calibration report renderer requires CalibrationReport")
    template_dir = Path(vntyper.__file__).resolve().parent / "templates"
    environment = Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=select_autoescape(["html"]),
        undefined=StrictUndefined,
        keep_trailing_newline=True,
    )
    return environment.get_template("calibration_report.html").render(report=report)


def write_calibration_report(destination: Path, report: CalibrationReport) -> None:
    """Write rendered report bytes into an atomically staged artifact directory."""
    if not isinstance(destination, Path):
        raise ValueError("calibration report destination must be a Path")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(render_calibration_report(report), encoding="utf-8")


def _tier_row(value: object) -> Mapping[str, object]:
    row = _exact(value, {"tier", "displayed", "exact", "wrong"}, "calibration tier metric")
    tier = row["tier"]
    if not isinstance(tier, str) or not tier:
        raise ValueError("calibration report tier must be a non-empty string")
    parsed: dict[str, object] = {"tier": tier}
    for field in ("displayed", "exact", "wrong"):
        item = row[field]
        if isinstance(item, bool) or not isinstance(item, int) or item < 0:
            raise ValueError(f"calibration report tier {field} must be a non-negative integer")
        parsed[field] = item
    return MappingProxyType(parsed)


def _abstention_row(value: object) -> Mapping[str, object]:
    row = _exact(value, {"split", "reason", "count", "rate"}, "calibration abstention row")
    parsed: dict[str, object] = {}
    for field in ("split", "reason", "rate"):
        item = row[field]
        if not isinstance(item, str) or not item:
            raise ValueError(f"calibration abstention {field} must be a non-empty string")
        parsed[field] = item
    count = row["count"]
    if isinstance(count, bool) or not isinstance(count, int) or count < 0:
        raise ValueError("calibration abstention count must be a non-negative integer")
    parsed["count"] = count
    return MappingProxyType(parsed)


def _string_sections(value: object, fields: set[str], label: str) -> dict[str, tuple[str, ...]]:
    raw = _exact(value, fields, f"calibration report {label}")
    return {field: _strings(raw[field], f"{label} {field}") for field in sorted(fields)}


def _strings(value: object, label: str) -> tuple[str, ...]:
    values = _sequence(value, label)
    if any(not isinstance(item, str) or not item for item in values):
        raise ValueError(f"calibration report {label} must contain non-empty strings")
    return tuple(values)  # type: ignore[arg-type]


def _nonempty_sequence(value: object, label: str) -> Sequence[object]:
    values = _sequence(value, label)
    if not values:
        raise ValueError(f"calibration report {label} must not be empty")
    return values


def _sequence(value: object, label: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"calibration report {label} must be a sequence")
    return value


def _digest(value: object, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"calibration report {label} SHA-256 must be lowercase hexadecimal")
    return value


def _exact(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        actual = sorted(value) if isinstance(value, Mapping) else type(value).__name__
        raise ValueError(f"{label} fields differ: expected {sorted(fields)}, got {actual}")
    return value
