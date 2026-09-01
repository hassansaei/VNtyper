"""Role-row dominance replay from allowlisted scalar calibration features."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction

from vntyper.scripts.calibration_objective import count_free_parameters


@dataclass(frozen=True)
class ScalarDecision:
    """Selected canonical identity or one closed whole-locus abstention reason."""

    selected_identity: str | None
    abstention_reason: str | None
    applicable: bool


def replay_scalar_dominance(features: Mapping[str, object], component: Mapping[str, object]) -> ScalarDecision:
    """Replay dominance fields that are losslessly represented by scalar features."""
    count_free_parameters(component)
    identity = features.get("canonical_identity")
    if not isinstance(identity, str) or not identity:
        return ScalarDecision(None, None, False)
    if component["enabled"] is not True:
        return ScalarDecision(None, None, True)
    required = (
        "haplotype_record_count_margin",
        "haplotype_record_share",
        "haplotype_record_share_margin",
        "haplotype_record_tie",
    )
    if any(name not in features for name in required):
        return ScalarDecision(None, None, False)
    if features["haplotype_record_tie"] is True:
        return ScalarDecision(None, "record-tie", True)
    count_margin = _nonnegative_fraction(features["haplotype_record_count_margin"], "record count margin")
    share = _unit_fraction(features["haplotype_record_share"], "record share")
    share_margin = _unit_fraction(features["haplotype_record_share_margin"], "record share margin")
    minimum_count = _nonnegative_fraction(component["minimum_record_count_margin"], "minimum record count margin")
    minimum_share = _unit_fraction(component["minimum_record_share"], "minimum record share")
    minimum_share_margin = _unit_fraction(component["minimum_record_share_margin"], "minimum record share margin")
    if count_margin < minimum_count or share < minimum_share or share_margin < minimum_share_margin:
        return ScalarDecision(None, "insufficient-dominance", True)
    veto = component["xd_veto"]
    if veto == "missingness" and _unit_fraction(features.get("xd_availability_fraction", 0), "XD availability") < 1:
        return ScalarDecision(None, "xd-missingness", True)
    if veto in {"concentration", "discordance"}:
        raise ValueError(f"calibration scalar evidence cannot evaluate the {veto} XD veto")
    if (
        component["abstain_on_inadmissible_advntr"] is True
        and features.get("advntr_evidence_disposition") != "admissible"
    ):
        return ScalarDecision(None, "inadmissible-advntr", True)
    return ScalarDecision(identity, None, True)


def _nonnegative_fraction(value: object, label: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"calibration {label} must be numeric")
    parsed = Fraction(str(value))
    if parsed < 0:
        raise ValueError(f"calibration {label} must be non-negative")
    return parsed


def _unit_fraction(value: object, label: str) -> Fraction:
    parsed = _nonnegative_fraction(value, label)
    if parsed > 1:
        raise ValueError(f"calibration {label} must not exceed one")
    return parsed
