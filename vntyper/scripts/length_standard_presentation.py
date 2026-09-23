"""Config-driven wording for packaged and locally fitted research length models."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from vntyper.scripts.length_presentation import LengthPresentation


def build_standard_length_presentation(
    summary: Mapping[str, object], report_config: Mapping[str, object]
) -> LengthPresentation | None:
    """Project a research result with explicit source and repeat-count convention.

    Args:
        summary: Pipeline summary containing a standard research model result.
        report_config: Configured report wording, including standard model vocabulary.

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
    if not isinstance(words, Mapping) or set(words) != {"help", "sources", "count_conventions", "extrapolation_help"}:
        raise ValueError("standard length report vocabulary differs")
    sources, conventions = words["sources"], words["count_conventions"]
    if not isinstance(sources, Mapping) or not isinstance(conventions, Mapping):
        raise ValueError("standard length report vocabularies must be mappings")
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
    help_text = " ".join(
        (
            _text(words["help"], "standard help"),
            _text(sources[source], "standard source"),
            _text(conventions[convention], "standard count convention"),
        )
    )
    from vntyper.scripts.length_sensitivity import CAUTION_CODE, HIGH_CODE

    # Sensitivity codes describe the array, not the features; only the rest mean extrapolation.
    if set(_reasons(summary.get("length_estimation_warnings", []))) - {CAUTION_CODE, HIGH_CODE}:
        help_text += " " + _text(words["extrapolation_help"], "standard extrapolation help")
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
