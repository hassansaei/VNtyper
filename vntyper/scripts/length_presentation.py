"""Pure report projection for optional total VNTR length results."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import Literal, cast

from vntyper.scripts.canonical_json import canonical_sha256

LengthDisplayStatus = Literal["not-recorded", "disabled", "measured-only", "estimated", "unavailable"]
_STATUSES = ("not-recorded", "disabled", "measured-only", "estimated", "unavailable")
_CONFIG_FIELDS = {"label", "unit", "help", "status_labels", "feature_labels"}
_FEATURE_LABEL_FIELDS = {"A", "F"}
_SUMMARY_FIELDS = {
    "length_estimation_status",
    "estimated_total_repeat_count",
    "length_calibration_id",
    "length_estimation_reasons",
    "length_features",
    "length_features_sha256",
    "length_model_sha256",
}
_ROW_FIELDS = {
    "manifest_key",
    "assembly",
    "assay_class",
    "input_scope",
    "annotation_sha256",
    "counting_policy_sha256",
    "provenance_sha256",
    "regions",
    "A",
    "F",
    "status",
    "reasons",
}
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")


@dataclass(frozen=True)
class LengthPresentation:
    """Closed, immutable values consumed by the per-sample report template."""

    status: LengthDisplayStatus
    status_text: str
    label: str
    value: str
    unit: str
    help: str
    feature_values: tuple[tuple[str, str], ...]
    reasons: tuple[str, ...]
    calibration_id: str | None
    model_sha256: str | None
    sensitivity_tier: str | None = None
    warning_badge: str | None = None
    notice_text: str | None = None
    uncertainty_text: str | None = None


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value.strip() != value:
        raise ValueError(f"length presentation {label} must be non-empty trimmed text")
    return value


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise ValueError(f"length presentation {label} must be lowercase SHA-256")
    return value


def _number(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"length presentation {label} must be finite numeric")
    return float(value)


def _shown(value: float, precision: int) -> str:
    return f"{value:.{precision}f}".rstrip("0").rstrip(".")


def _configuration(
    report_config: Mapping[str, object],
) -> tuple[str, str, str, Mapping[str, object], Mapping[str, object]] | None:
    raw = report_config.get("length_estimation")
    if raw is None:
        return None
    if not isinstance(raw, Mapping) or set(raw) != _CONFIG_FIELDS:
        raise ValueError("length presentation configuration fields differ from the closed contract")
    statuses = raw["status_labels"]
    features = raw["feature_labels"]
    if not isinstance(statuses, Mapping) or set(statuses) != set(_STATUSES):
        raise ValueError("length presentation status label configuration differs")
    if not isinstance(features, Mapping) or set(features) != _FEATURE_LABEL_FIELDS:
        raise ValueError("length presentation feature label configuration differs")
    for key in _STATUSES:
        _text(statuses[key], f"status label {key}")
    for key in _FEATURE_LABEL_FIELDS:
        _text(features[key], f"feature label {key}")
    return (
        _text(raw["label"], "label"),
        _text(raw["unit"], "unit"),
        _text(raw["help"], "help"),
        statuses,
        features,
    )


def _feature_values(
    raw: object,
    expected_sha256: object,
    labels: Mapping[str, object],
) -> tuple[tuple[str, str], ...]:
    if not isinstance(raw, dict) or set(raw) != {"schema_version", "rows"}:
        raise ValueError("length presentation features differ from the closed artifact contract")
    if raw["schema_version"] != "length-features-v1":
        raise ValueError("length presentation features schema is unsupported")
    rows = raw["rows"]
    if not isinstance(rows, list) or len(rows) != 1 or not isinstance(rows[0], dict) or set(rows[0]) != _ROW_FIELDS:
        raise ValueError("length presentation requires exactly one closed feature row")
    if canonical_sha256(raw) != _digest(expected_sha256, "features digest"):
        raise ValueError("length presentation features digest differs from content")
    row = rows[0]
    values: list[tuple[str, str]] = []
    for name in ("A", "F"):
        value = row[name]
        if value is not None:
            values.append((_text(labels[name], f"feature label {name}"), _shown(_number(value, name), 3)))
    return tuple(values)


def _reasons(value: object) -> tuple[str, ...]:
    if not isinstance(value, list) or any(not isinstance(reason, str) or not reason for reason in value):
        raise ValueError("length presentation reasons must be a list of non-empty codes")
    if len(value) != len(set(value)):
        raise ValueError("length presentation reasons must be unique")
    return tuple(value)


def build_length_presentation(
    summary: Mapping[str, object],
    report_config: Mapping[str, object],
    *,
    is_positive: bool | None = None,
) -> LengthPresentation | None:
    """Validate and format optional length summary values for the report.

    Args:
        summary: Loaded pipeline summary mapping.
        report_config: Report wording configuration.
        is_positive: Whether the overall screening finding is positive.

    Returns:
        A frozen presentation, or ``None`` when an older config has no wording.

    Raises:
        ValueError: If configured wording or recorded length state is inconsistent.
    """
    if not isinstance(summary, Mapping) or not isinstance(report_config, Mapping):
        raise ValueError("length presentation inputs must be mappings")
    presentation = _build_base(summary, report_config)
    if presentation is None:
        return None
    from vntyper.scripts.length_sensitivity import build_sensitivity_view

    view = build_sensitivity_view(summary, report_config, is_positive=is_positive)
    if view is None:
        return presentation
    return replace(
        presentation,
        help=presentation.help if view.help is None else f"{presentation.help} {view.help}",
        sensitivity_tier=view.tier,
        warning_badge=view.badge,
        notice_text=view.notice,
        uncertainty_text=view.uncertainty,
    )


def _build_base(summary: Mapping[str, object], report_config: Mapping[str, object]) -> LengthPresentation | None:
    if "length_model_source" in summary:
        from vntyper.scripts.length_standard_presentation import build_standard_length_presentation

        return build_standard_length_presentation(summary, report_config)
    configured = _configuration(report_config)
    if configured is None:
        return None
    label, unit, help_text, status_labels, feature_labels = configured
    raw_status: object
    if "length_estimation_status" not in summary:
        partial = sorted((_SUMMARY_FIELDS - {"length_estimation_status"}) & set(summary))
        if partial:
            raise ValueError(f"legacy length presentation contains fields without status: {partial}")
        raw_status = "not-recorded"
    else:
        raw_status = summary["length_estimation_status"]
    if not isinstance(raw_status, str) or raw_status not in _STATUSES:
        raise ValueError("length presentation status is unsupported")
    status = cast(LengthDisplayStatus, raw_status)
    if status == "not-recorded":
        return LengthPresentation(
            status,
            _text(status_labels[status], "status text"),
            label,
            _text(status_labels[status], "legacy value"),
            unit,
            help_text,
            (),
            (),
            None,
            None,
        )

    estimate = summary.get("estimated_total_repeat_count")
    calibration_id = summary.get("length_calibration_id")
    model_sha256 = summary.get("length_model_sha256")
    reasons = _reasons(summary.get("length_estimation_reasons"))
    features = summary.get("length_features")
    features_sha256 = summary.get("length_features_sha256")
    if status == "disabled":
        if (
            any(value is not None for value in (estimate, calibration_id, model_sha256, features, features_sha256))
            or reasons
        ):
            raise ValueError("disabled length presentation contains measured or model values")
        feature_values: tuple[tuple[str, str], ...] = ()
    else:
        feature_values = _feature_values(features, features_sha256, feature_labels)
    if status == "measured-only":
        if estimate is not None or calibration_id is not None or model_sha256 is not None:
            raise ValueError("measurement-only length presentation contains model output")
    elif status in {"estimated", "unavailable"}:
        calibration_id = _text(calibration_id, "calibration ID")
        model_sha256 = _digest(model_sha256, "model digest")
        if status == "estimated" and reasons:
            raise ValueError("estimated length presentation cannot contain unavailability reasons")
        if status == "unavailable" and (estimate is not None or not reasons):
            raise ValueError("unavailable length presentation requires reasons and no estimate")
    if status == "estimated":
        prediction = _number(estimate, "estimate")
        if prediction <= 0:
            raise ValueError("length presentation estimate must be positive")
        value = _shown(prediction, 2)
    else:
        if estimate is not None:
            raise ValueError("non-estimated length presentation cannot contain an estimate")
        value = "Not calculated"
    return LengthPresentation(
        status,
        _text(status_labels[status], "status text"),
        label,
        value,
        unit,
        help_text,
        feature_values,
        reasons,
        cast(str | None, calibration_id),
        cast(str | None, model_sha256),
    )
