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
from typing import Literal

CAUTION_CODE = "vntr_length_exceeds_sensitivity_cutoff"
HIGH_CODE = "vntr_length_high_sensitivity_risk"
SensitivityTier = Literal["below", "caution", "high", "not-assessed"]
TIERS: tuple[SensitivityTier, ...] = ("below", "caution", "high", "not-assessed")
_POLICY_FIELDS = ("caution_threshold", "high_threshold", "uncertainty_repeats")


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


def classify_length_sensitivity(status: object, estimate: object, policy: LengthSensitivityPolicy) -> SensitivityTier:
    """Tier a recorded estimate; anything but a finite ``estimated`` value is not assessed.

    Args:
        status: Recorded ``length_estimation_status``.
        estimate: Recorded ``estimated_total_repeat_count``.
        policy: The cutoffs to compare against.

    Returns:
        The sensitivity tier.
    """
    if (
        status != "estimated"
        or isinstance(estimate, bool)
        or not isinstance(estimate, (int, float))
        or not math.isfinite(estimate)
    ):
        return "not-assessed"
    if estimate > policy.high_threshold:
        return "high"
    if estimate > policy.caution_threshold:
        return "caution"
    return "below"


def apply_length_sensitivity(fields: Mapping[str, object], policy: LengthSensitivityPolicy | None) -> dict[str, object]:
    """Return length summary fields with the tier, policy and warning codes added.

    Args:
        fields: Length summary fields from either length path.
        policy: The configured policy, or ``None`` when the configuration has none.

    Returns:
        A new mapping. Disabled runs and a missing policy leave the fields unchanged.
    """
    result = dict(fields)
    if policy is None or fields.get("length_estimation_status") == "disabled":
        return result
    tier = classify_length_sensitivity(
        fields.get("length_estimation_status"), fields.get("estimated_total_repeat_count"), policy
    )
    recorded = fields.get("length_estimation_warnings") or []
    warnings = [str(code) for code in recorded] if isinstance(recorded, list) else []
    codes = {"caution": (CAUTION_CODE,), "high": (CAUTION_CODE, HIGH_CODE)}.get(tier, ())
    warnings.extend(code for code in codes if code not in warnings)
    result["length_estimation_warnings"] = warnings
    result["length_sensitivity_tier"] = tier
    result["length_sensitivity_policy"] = policy.as_dict()
    return result
