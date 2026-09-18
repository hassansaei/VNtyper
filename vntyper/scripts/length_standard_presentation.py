"""Config-driven wording for packaged and locally fitted research length models."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from vntyper.scripts.length_presentation import LengthPresentation


def build_standard_length_presentation(
    summary: Mapping[str, object],
    report_config: Mapping[str, object],
    *,
    is_positive: bool | None = None,
) -> LengthPresentation | None:
    """Project a research result with explicit source and repeat-count convention.

    Args:
        summary: Pipeline summary containing a standard research model result.
        report_config: Configured report wording, including standard model vocabulary.
        is_positive: Whether the overall screening finding is positive.

    Returns:
        Validated report presentation, or None for an older report configuration.

    Raises:
        ValueError: If the summary is inconsistent, forged as approved, or malformed.
    """
    from vntyper.scripts.length_presentation import (
        LengthPresentation,
        _configuration,
        _digest,
        _number,
        _reasons,
        _shown,
        _text,
    )

    words = report_config.get("standard_length_estimation")
    configured = _configuration(report_config)
    if words is None or configured is None:
        return None
    base_keys = {"help", "sources", "count_conventions", "extrapolation_help"}
    if (
        not isinstance(words, Mapping)
        or not base_keys.issubset(set(words))
        or not set(words).issubset(base_keys | {"sensitivity_warning"})
    ):
        raise ValueError("standard length report vocabulary differs")
    sources, conventions = words["sources"], words["count_conventions"]
    if not isinstance(sources, Mapping) or not isinstance(conventions, Mapping):
        raise ValueError("standard length report vocabularies must be mappings")
    sensitivity_cfg = words.get("sensitivity_warning")
    if sensitivity_cfg is not None:
        if not isinstance(sensitivity_cfg, Mapping) or set(sensitivity_cfg) != {
            "warning_code",
            "badge",
            "notice_negative",
            "notice_positive",
            "help",
        }:
            raise ValueError("standard length sensitivity warning vocabulary differs")
        for key in ("warning_code", "badge", "notice_negative", "notice_positive", "help"):
            _text(sensitivity_cfg[key], f"sensitivity warning {key}")
    source = summary.get("length_model_source")
    convention = summary.get("length_count_convention")
    if source not in ("packaged-research", "local-research") or source not in sources:
        raise ValueError("standard length model source differs")
    if not isinstance(convention, str) or convention not in conventions:
        raise ValueError("standard length count convention differs")
    if summary.get("length_model_evidence_status") != "research-development":
        raise ValueError("standard length evidence status differs")
    if any(summary.get(key) is not None for key in ("length_calibration_id", "length_portable_approval_sha256")):
        raise ValueError("standard research length result cannot claim approval")
    status = summary.get("length_estimation_status")
    if status not in ("estimated", "unavailable"):
        raise ValueError("standard length report status differs")
    reasons = _reasons(summary.get("length_estimation_reasons"))
    estimate = summary.get("estimated_total_repeat_count")
    if status == "unavailable" and (estimate is not None or not reasons):
        raise ValueError("unavailable standard length requires reasons and no estimate")
    model_sha = _digest(summary.get("length_model_sha256"), "standard model digest")
    label, unit, _, status_labels, feature_labels = configured
    features = _features(summary, feature_labels)
    if status == "estimated":
        value = _number(estimate, "standard estimate")
        if value <= 0 or summary.get("length_standard_features") is None:
            raise ValueError("standard length estimate requires positive count and measured features")
        shown = _shown(value, 2)
    else:
        shown = _text(status_labels[status], "standard unavailable value")
    warnings = _reasons(summary.get("length_estimation_warnings", []))
    has_extrapolation = any(w.startswith(("feature_", "feature-range")) for w in warnings)
    sensitivity_code = (
        sensitivity_cfg["warning_code"] if sensitivity_cfg is not None else "vntr_length_exceeds_sensitivity_cutoff"
    )
    has_sensitivity_warning = status == "estimated" and sensitivity_cfg is not None and sensitivity_code in warnings

    warning_badge: str | None = None
    notice_text: str | None = None

    if has_sensitivity_warning and sensitivity_cfg is not None:
        threshold_val = summary.get("length_warning_threshold", 110.0)
        threshold_num = _number(threshold_val, "length warning threshold")
        threshold_str = str(int(threshold_num)) if threshold_num.is_integer() else _shown(threshold_num, 1)
        badge_template = _text(sensitivity_cfg["badge"], "sensitivity warning badge")
        notice_template = _text(
            sensitivity_cfg["notice_positive"] if is_positive else sensitivity_cfg["notice_negative"],
            "sensitivity warning notice",
        )
        warning_badge = badge_template.format(threshold=threshold_str, estimate=shown)
        notice_text = notice_template.format(threshold=threshold_str, estimate=shown)

    help_parts = [
        _text(words["help"], "standard help"),
        _text(sources[source], "standard source"),
        _text(conventions[convention], "standard count convention"),
    ]
    if has_extrapolation:
        help_parts.append(_text(words["extrapolation_help"], "standard extrapolation help"))
    if has_sensitivity_warning and sensitivity_cfg is not None:
        threshold_val = summary.get("length_warning_threshold", 110.0)
        threshold_num = _number(threshold_val, "length warning threshold")
        threshold_str = str(int(threshold_num)) if threshold_num.is_integer() else _shown(threshold_num, 1)
        help_parts.append(
            _text(sensitivity_cfg["help"], "sensitivity help").format(threshold=threshold_str, estimate=shown)
        )

    help_text = " ".join(help_parts)
    return LengthPresentation(
        status,
        _text(status_labels[status], "standard status"),
        label,
        shown,
        unit,
        help_text,
        features,
        reasons,
        None,
        model_sha,
        warning_badge=warning_badge,
        notice_text=notice_text,
        has_sensitivity_warning=has_sensitivity_warning,
    )


def _features(summary: Mapping[str, object], labels: Mapping[str, object]) -> tuple[tuple[str, str], ...]:
    from vntyper.scripts.length_presentation import _digest, _shown, _text

    raw = summary.get("length_standard_features")
    expected = summary.get("length_standard_features_sha256")
    if raw is None:
        if expected is not None:
            raise ValueError("missing standard features cannot carry a digest")
        return ()
    from vntyper.scripts.length_standard_features import decode_standard_length_measurement

    measurement = decode_standard_length_measurement(raw)
    if measurement.sha256 != _digest(expected, "standard features digest"):
        raise ValueError("standard length feature digest differs")
    values: list[tuple[str, str]] = []
    for name in ("A", "F"):
        value = measurement.values[name]
        if value is not None:
            values.append((_text(labels[name], "standard feature label"), _shown(value, 3)))
    return tuple(values)
