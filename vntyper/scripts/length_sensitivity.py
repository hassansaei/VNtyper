"""Two-tier sensitivity policy for the total MUC1 VNTR length estimate.

Short-read frameshift detection loses sensitivity as the total (diploid) array grows,
because Kestrel's depth ratio scales as 1/array length. This module turns a recorded
length estimate into a tier using a configured policy. It is pure: the pipeline applies
it once at the summary merge point, and the report re-checks what was recorded.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal, cast

CAUTION_CODE = "vntr_length_exceeds_sensitivity_cutoff"
HIGH_CODE = "vntr_length_high_sensitivity_risk"
SensitivityTier = Literal["below", "caution", "high", "not-assessed"]
TIERS: tuple[SensitivityTier, ...] = ("below", "caution", "high", "not-assessed")
_POLICY_FIELDS = ("caution_threshold", "high_threshold", "uncertainty_repeats")
#: Cutoffs are defined on complete counts (all units, including the nine invariant
#: terminal units per allele). A canonical-only estimate is moved into that frame by the
#: packaged ``canonical-only-plus-nine-terminals-v1`` conversion; a source-reported one
#: has no known frame and is not assessed. Approved-path models are complete by contract
#: and record no convention.
_COMPLETE_FRAME_OFFSET: dict[object, float] = {None: 0.0, "complete": 0.0, "canonical-only": 18.0}


@dataclass(frozen=True)
class LengthSensitivityPolicy:
    """Validated cutoffs (repeat units, strict ``>``) and the displayed estimate uncertainty."""

    caution_threshold: float
    high_threshold: float
    uncertainty_repeats: float

    def as_dict(self) -> dict[str, float]:
        """Return the policy as recorded in ``pipeline_summary.json``."""
        return {name: getattr(self, name) for name in _POLICY_FIELDS}


def _positive(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"length sensitivity {name} must be a finite positive number")
    return float(value)


def decode_length_sensitivity_policy(raw: object) -> LengthSensitivityPolicy:
    """Validate a policy mapping with exactly the three closed fields.

    Args:
        raw: The configured or recorded policy mapping.

    Returns:
        The validated policy.

    Raises:
        ValueError: On missing or extra keys, non-numeric, non-finite or non-positive
            values, or a caution cutoff not below the high cutoff.
    """
    if not isinstance(raw, Mapping) or set(raw) != set(_POLICY_FIELDS):
        raise ValueError("length sensitivity policy fields differ from the closed contract")
    policy = LengthSensitivityPolicy(*(_positive(raw[name], name) for name in _POLICY_FIELDS))
    if policy.caution_threshold >= policy.high_threshold:
        raise ValueError("length sensitivity caution threshold must be below the high threshold")
    return policy


def resolve_length_sensitivity_policy(config: Mapping[str, object]) -> LengthSensitivityPolicy | None:
    """Read ``length_estimation.sensitivity`` from the run configuration.

    Args:
        config: The loaded run configuration.

    Returns:
        The validated policy, or ``None`` for a configuration without the block.

    Raises:
        ValueError: If the section or the policy is malformed.
    """
    section = config.get("length_estimation", {})
    if not isinstance(section, Mapping):
        raise ValueError("length_estimation configuration must be an object")
    raw = section.get("sensitivity")
    return None if raw is None else decode_length_sensitivity_policy(raw)


def classify_length_sensitivity(
    status: object, estimate: object, policy: LengthSensitivityPolicy, convention: object = None
) -> SensitivityTier:
    """Tier a recorded estimate in the complete-count frame the cutoffs are defined on.

    Args:
        status: Recorded ``length_estimation_status``.
        estimate: Recorded ``estimated_total_repeat_count``.
        policy: The cutoffs to compare against.
        convention: Recorded ``length_count_convention``; ``None`` for the approved path.

    Returns:
        The sensitivity tier; ``not-assessed`` for anything but a finite ``estimated``
        value in a count frame that can be placed on the complete frame.
    """
    if (
        convention not in _COMPLETE_FRAME_OFFSET
        or status != "estimated"
        or isinstance(estimate, bool)
        or not isinstance(estimate, (int, float))
        or not math.isfinite(estimate)
    ):
        return "not-assessed"
    complete = _complete_tenth(estimate, convention)
    if complete > policy.high_threshold:
        return "high"
    if complete > policy.caution_threshold:
        return "caution"
    return "below"


def _complete_tenth(estimate: float, convention: object) -> float:
    """The estimate in the complete-count frame, at the tenth the report displays.

    Tiering on the displayed value means a report can never show "150.0 exceeds 150".
    """
    return round(estimate + _COMPLETE_FRAME_OFFSET[convention], 1)


def _expected_codes(tier: str) -> tuple[str, ...]:
    return {"caution": (CAUTION_CODE,), "high": (CAUTION_CODE, HIGH_CODE)}.get(tier, ())


def apply_length_sensitivity(fields: Mapping[str, object], policy: LengthSensitivityPolicy | None) -> dict[str, object]:
    """Return length summary fields with the tier, policy and warning codes added.

    Args:
        fields: Length summary fields from either length path.
        policy: The configured policy, or ``None`` when the configuration has none.

    Returns:
        A new mapping. Disabled runs and a missing policy leave the fields unchanged.

    Raises:
        ValueError: If recorded warnings are not a list of codes.
    """
    result = dict(fields)
    if policy is None or fields.get("length_estimation_status") == "disabled":
        return result
    tier = classify_length_sensitivity(
        fields.get("length_estimation_status"),
        fields.get("estimated_total_repeat_count"),
        policy,
        fields.get("length_count_convention"),
    )
    recorded = fields.get("length_estimation_warnings", [])
    if not isinstance(recorded, list) or any(not isinstance(code, str) for code in recorded):
        raise ValueError("length estimation warnings must be a list of codes")
    warnings = list(recorded)
    warnings.extend(code for code in _expected_codes(tier) if code not in warnings)
    result["length_estimation_warnings"] = warnings
    result["length_sensitivity_tier"] = tier
    result["length_sensitivity_policy"] = policy.as_dict()
    return result


_WORD_FIELDS: dict[str, frozenset[str]] = {
    "caution": frozenset({"badge", "help"}),
    "high": frozenset({"badge", "help", "notice_not_positive"}),
    "labels": frozenset({"notice_prefix", "cohort_kpi_label", "cohort_kpi_detail"}),
}


@dataclass(frozen=True)
class SensitivityView:
    """Report-ready sensitivity wording; ``None`` members render nothing."""

    tier: SensitivityTier
    badge: str | None
    notice: str | None
    help: str | None
    uncertainty: str | None
    value: str | None = None
    notice_prefix: str | None = None


def _words(report_config: Mapping[str, object]) -> dict[str, dict[str, str]] | None:
    raw = report_config.get("length_sensitivity")
    if raw is None:
        return None
    if not isinstance(raw, Mapping) or set(raw) != set(_WORD_FIELDS):
        raise ValueError("length sensitivity wording differs from the closed contract")
    words: dict[str, dict[str, str]] = {}
    for tier, fields in _WORD_FIELDS.items():
        block = raw[tier]
        if not isinstance(block, Mapping) or set(block) != fields:
            raise ValueError(f"length sensitivity {tier} wording differs from the closed contract")
        words[tier] = {}
        for name in fields:
            text = block[name]
            if not isinstance(text, str) or not text or text.strip() != text:
                raise ValueError(f"length sensitivity {tier} {name} must be non-empty trimmed text")
            words[tier][name] = text
    return words


def _count(value: float) -> str:
    return str(int(value)) if value.is_integer() else f"{value:.1f}"


def build_sensitivity_view(
    summary: Mapping[str, object], report_config: Mapping[str, object], *, is_positive: bool | None
) -> SensitivityView | None:
    """Re-check the recorded tier and project the configured wording for it.

    Args:
        summary: Loaded pipeline summary mapping.
        report_config: Report wording configuration.
        is_positive: Whether the overall screening finding is positive. The banner
            notice is produced only when it is not.

    Returns:
        The view, or ``None`` for a summary without a recorded tier or a report
        configuration without sensitivity wording.

    Raises:
        ValueError: If the recorded tier, policy or configured wording is malformed, or
            the tier differs from the recorded estimate and policy.
    """
    if "length_sensitivity_tier" not in summary:
        return None
    words = _words(report_config)
    if words is None:
        return None
    policy = decode_length_sensitivity_policy(summary.get("length_sensitivity_policy"))
    tier = summary["length_sensitivity_tier"]
    estimate = summary.get("estimated_total_repeat_count")
    expected = classify_length_sensitivity(
        summary.get("length_estimation_status"), estimate, policy, summary.get("length_count_convention")
    )
    if tier != expected:
        raise ValueError("recorded length sensitivity tier differs from the recorded estimate and policy")
    recorded = summary.get("length_estimation_warnings", [])
    recorded_codes = (
        {code for code in recorded if code in (CAUTION_CODE, HIGH_CODE)} if isinstance(recorded, list) else None
    )
    if recorded_codes != set(_expected_codes(expected)):
        raise ValueError("recorded length sensitivity warning codes differ from the recorded tier")
    if expected == "not-assessed":
        return SensitivityView(expected, None, None, None, None)
    # The recorded uncertainty is the packaged model's own held-out error; other models
    # never earned it, so their estimates are shown bare.
    uncertainty = (
        f"± {_count(float(round(policy.uncertainty_repeats)))}"
        if summary.get("length_model_source") == "packaged-research"
        else None
    )
    # Beside a ±14 error bar, a second decimal is precision the estimate does not have.
    value = None if uncertainty is None else f"{cast(float, estimate):.1f}"
    if expected == "below":
        return SensitivityView(expected, None, None, None, uncertainty, value)
    block = words[expected]
    threshold = policy.high_threshold if expected == "high" else policy.caution_threshold
    shown = _count(_complete_tenth(cast(float, estimate), summary.get("length_count_convention")))
    values = {
        "threshold": _count(threshold),
        "estimate": shown if uncertainty is None else f"{shown} {uncertainty}",
    }
    notice = block["notice_not_positive"].format(**values) if expected == "high" and not is_positive else None
    return SensitivityView(
        expected,
        block["badge"].format(**values),
        notice,
        block["help"].format(**values),
        uncertainty,
        value,
        words["labels"]["notice_prefix"] if notice is not None else None,
    )


def cohort_kpi_text(report_config: Mapping[str, object], counts: Mapping[str, int] | None) -> tuple[str, str] | None:
    """Configured label and detail line for the cohort's high-tier count.

    Args:
        report_config: Report wording configuration.
        counts: ``high``/``caution``/``assessed`` counts, or None when nothing was assessed.

    Returns:
        ``(label, detail)``, or None when there are no counts or no configured wording.
    """
    words = _words(report_config)
    if words is None or counts is None:
        return None
    labels = words["labels"]
    return labels["cohort_kpi_label"], labels["cohort_kpi_detail"].format(**counts)
