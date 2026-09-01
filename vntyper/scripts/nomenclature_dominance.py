"""Pure record-vote dominance and closed whole-locus abstention decisions."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from vntyper.scripts.molecular_identity import EvidenceDisposition, MolecularIdentity
from vntyper.scripts.nomenclature_bam_evidence import BamLocusEvidence

DominanceOutcome = Literal["selected", "abstained", "not-applicable"]
AbstentionReason = Literal[
    "record-tie",
    "insufficient-dominance",
    "xd-missingness",
    "xd-concentration",
    "xd-discordance",
    "inadmissible-advntr",
]

_COMPONENT_FIELDS = {
    "enabled",
    "minimum_record_count_margin",
    "minimum_record_share",
    "minimum_record_share_margin",
    "xd_veto",
    "abstain_on_inadmissible_advntr",
}
_XD_VETOES = frozenset({"disabled", "missingness", "concentration", "discordance"})
_OUTCOMES = frozenset({"selected", "abstained", "not-applicable"})
_ABSTENTION_REASONS = frozenset(
    {
        "record-tie",
        "insufficient-dominance",
        "xd-missingness",
        "xd-concentration",
        "xd-discordance",
        "inadmissible-advntr",
    }
)


@dataclass(frozen=True)
class DominanceEvidence:
    """Replay evidence required by one dominance decision."""

    bam_evidence: BamLocusEvidence | None
    advntr_disposition: EvidenceDisposition

    def __post_init__(self) -> None:
        """Validate the optional BAM evidence and governed adVNTR state."""
        if self.bam_evidence is not None and not isinstance(self.bam_evidence, BamLocusEvidence):
            raise ValueError("dominance BAM evidence must be a BamLocusEvidence or None")
        if not isinstance(self.advntr_disposition, EvidenceDisposition):
            raise ValueError("dominance adVNTR disposition must be an EvidenceDisposition")


@dataclass(frozen=True)
class DominanceDecision:
    """Selected identity or one explicit closed abstention outcome."""

    outcome: DominanceOutcome
    identity: MolecularIdentity | None
    abstention_reason: AbstentionReason | None
    top_record_count: int = 0
    runner_up_record_count: int = 0
    record_count_margin: int = 0
    top_record_share: Fraction = Fraction(0)
    runner_up_record_share: Fraction = Fraction(0)
    record_share_margin: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        """Keep selected, abstained, and non-applicable states disjoint."""
        if self.outcome not in _OUTCOMES:
            raise ValueError(f"unsupported dominance outcome: {self.outcome}")
        if self.identity is not None and not isinstance(self.identity, MolecularIdentity):
            raise ValueError("dominance identity must be a MolecularIdentity or None")
        if self.abstention_reason is not None and self.abstention_reason not in _ABSTENTION_REASONS:
            raise ValueError(f"unsupported dominance abstention reason: {self.abstention_reason}")
        if self.outcome == "selected" and (self.identity is None or self.abstention_reason is not None):
            raise ValueError("selected dominance decisions require an identity and no abstention reason")
        if self.outcome == "abstained" and (self.identity is not None or self.abstention_reason is None):
            raise ValueError("abstained dominance decisions require one reason and no identity")
        if self.outcome == "not-applicable" and (self.identity is not None or self.abstention_reason is not None):
            raise ValueError("non-applicable dominance decisions cannot select or abstain")


@dataclass(frozen=True)
class _DominancePolicy:
    enabled: bool
    minimum_record_count_margin: int
    minimum_record_share: Fraction
    minimum_record_share_margin: Fraction
    xd_veto: str
    abstain_on_inadmissible_advntr: bool


def evaluate_dominance(evidence: DominanceEvidence, component: Mapping[str, object]) -> DominanceDecision:
    """Evaluate record dominance followed by governed whole-locus vetoes.

    Args:
        evidence: Replayable record evidence and governed adVNTR disposition.
        component: Complete resolved ``dominance`` profile component.

    Returns:
        A selected identity, explicit abstention, or non-applicable outcome.

    Raises:
        ValueError: If evidence or component values violate the closed contract.
    """
    if not isinstance(evidence, DominanceEvidence):
        raise ValueError("dominance evidence must be a DominanceEvidence")
    policy = _decode_policy(component)
    bam = evidence.bam_evidence
    if not policy.enabled or bam is None or bam.eligible_record_count == 0 or not bam.counts:
        return DominanceDecision("not-applicable", None, None)

    ranked = sorted(bam.counts.items(), key=lambda item: item[1], reverse=True)
    top_identity, top_count = ranked[0]
    runner_up_count = ranked[1][1] if len(ranked) > 1 else 0
    denominator = bam.eligible_record_count
    count_margin = top_count - runner_up_count
    top_share = Fraction(top_count, denominator)
    runner_up_share = Fraction(runner_up_count, denominator)
    share_margin = Fraction(count_margin, denominator)
    metrics = (top_count, runner_up_count, count_margin, top_share, runner_up_share, share_margin)
    if top_count == runner_up_count:
        return DominanceDecision("abstained", None, "record-tie", *metrics)
    if (
        count_margin < policy.minimum_record_count_margin
        or top_share < policy.minimum_record_share
        or share_margin < policy.minimum_record_share_margin
    ):
        return DominanceDecision("abstained", None, "insufficient-dominance", *metrics)
    if policy.abstain_on_inadmissible_advntr and evidence.advntr_disposition.value != "admissible":
        return DominanceDecision("abstained", None, "inadmissible-advntr", *metrics)

    xd_reason = _xd_abstention_reason(bam, top_identity, policy.xd_veto)
    if xd_reason is not None:
        return DominanceDecision("abstained", None, xd_reason, *metrics)
    return DominanceDecision("selected", top_identity, None, *metrics)


def _decode_policy(component: Mapping[str, object]) -> _DominancePolicy:
    if not isinstance(component, Mapping) or set(component) != _COMPONENT_FIELDS:
        actual = sorted(component) if isinstance(component, Mapping) else type(component).__name__
        raise ValueError(f"dominance component fields differ: {actual}")
    enabled = component["enabled"]
    count_margin = component["minimum_record_count_margin"]
    share = component["minimum_record_share"]
    share_margin = component["minimum_record_share_margin"]
    xd_veto = component["xd_veto"]
    advntr_veto = component["abstain_on_inadmissible_advntr"]
    if not isinstance(enabled, bool):
        raise ValueError("dominance enabled must be Boolean")
    if isinstance(count_margin, bool) or not isinstance(count_margin, int) or count_margin < 0:
        raise ValueError("dominance minimum record-count margin must be a non-negative integer")
    parsed_share = _unit_fraction(share, "minimum record share")
    parsed_share_margin = _unit_fraction(share_margin, "minimum record-share margin")
    if not isinstance(xd_veto, str) or xd_veto not in _XD_VETOES:
        raise ValueError(f"dominance XD veto is unsupported: {xd_veto!r}")
    if not isinstance(advntr_veto, bool):
        raise ValueError("dominance adVNTR veto must be Boolean")
    return _DominancePolicy(enabled, count_margin, parsed_share, parsed_share_margin, xd_veto, advntr_veto)


def _unit_fraction(value: object, label: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"dominance {label} must be numeric")
    parsed = Fraction(str(value))
    if not 0 <= parsed <= 1:
        raise ValueError(f"dominance {label} must be between zero and one")
    return parsed


def _xd_abstention_reason(
    evidence: BamLocusEvidence,
    winner: MolecularIdentity,
    mode: str,
) -> AbstentionReason | None:
    if mode == "disabled":
        return None
    winner_values = tuple(record.minimum_kmer_depth for record in evidence.records if winner in record.identities)
    if mode == "missingness":
        return "xd-missingness" if any(value is None for value in winner_values) else None
    available_winner = tuple(value for value in winner_values if value is not None)
    if mode == "concentration":
        total = sum(available_winner)
        return "xd-concentration" if total > 0 and max(available_winner, default=0) * 2 > total else None
    winner_total = sum(available_winner)
    runner_totals = (
        sum(record.minimum_kmer_depth or 0 for record in evidence.records if identity in record.identities)
        for identity in evidence.counts
        if identity != winner
    )
    return "xd-discordance" if any(total >= winner_total for total in runner_totals) else None
