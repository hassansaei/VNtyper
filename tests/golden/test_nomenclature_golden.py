"""Known-truth golden set: 200 simulated samples, both callers.

These are the numbers the design argues from, pinned so they cannot regress
silently. The benchmark tree is supplied out of band -- ``VNTYPER_SIM_ROOT`` for the
simulated cohort, ``VNTYPER_ADVNTR_ROOT`` for the adVNTR batch -- so no path to
either is committed and the tier skips wherever they are absent.

Not part of ``make test-unit``, which selects ``-m unit tests/unit``. Run with::

    VNTYPER_SIM_ROOT=... VNTYPER_ADVNTR_ROOT=... pytest -m golden tests/golden

Research use only.
"""

from __future__ import annotations

import csv
import os
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import pytest

from vntyper.scripts.nomenclature import Nomenclature, from_advntr, from_kestrel, reconcile, render
from vntyper.scripts.nomenclature_bam import BamRescuer, from_bam, is_candidate, refine

pytestmark = pytest.mark.golden

SIM_ROOT = os.environ.get("VNTYPER_SIM_ROOT")
ADVNTR_ROOT = os.environ.get("VNTYPER_ADVNTR_ROOT")

EXPERIMENTS = ("experiment1_dupC", "experiment2_atypical")

#: The name each truth class should carry, derived from the simulator's own mutation
#: definitions rather than guessed.
EXPECTED_NAME: dict[str, str] = {
    "dupC": "59dupC",
    "insCCCC": "56_59dupCCCC",
    "insG": "58_59insG",
    "dupA": "60dupA",
    "delinsAT": "54_56delinsAT",
    "delGCCCA": "1_5delGCCCA",
    "ins25bp": "30_31insCAGGCCGGCCCCGGGCTCCGGACAC",
    "insC_pos23": "23dupC",
    "insG_pos58": "57_58insG",
    "insG_pos54": "53_54insG",
    "insA_pos54": "53_54insA",
}

#: Spec §2.6, the VCF-only column: how many of each class the Kestrel path names
#: correctly. A floor, not a ceiling -- an improvement must not fail the build --
#: except ``dupC``, which is pinned exactly because it is the canonical allele.
#:
#: ``dupA`` earns its place here. Deliberately shifting the insertion anchor by one
#: base -- the off-by-one that reading an interbase anchor as a base would cause --
#: leaves ``dupC`` at 96/96, because the 7xC tract absorbs the shift when the name is
#: normalised. ``dupA`` sits at position 60 against a non-repeating base and cannot
#: roll, so it is what actually detects that class of bug. Verified by breaking the
#: rule and watching this floor fail while the canonical case stayed green.
VCF_ONLY_FLOOR: dict[str, int] = {
    "dupC": 96,
    "insCCCC": 9,
    "insC_pos23": 9,
    "dupA": 6,
    "delGCCCA": 3,
    "insG_pos54": 0,
    "delinsAT": 0,
    "insG": 0,
    "insA_pos54": 0,
    "insG_pos58": 0,
    "ins25bp": 0,
}


def _require_sim() -> Path:
    if not SIM_ROOT:
        pytest.skip("VNTYPER_SIM_ROOT unset")
    root = Path(SIM_ROOT)
    if not root.is_dir():
        pytest.skip(f"benchmark root not found: {root}")
    return root


def _require_advntr() -> Path:
    if not ADVNTR_ROOT:
        pytest.skip("VNTYPER_ADVNTR_ROOT unset")
    root = Path(ADVNTR_ROOT)
    if not root.is_dir():
        pytest.skip(f"adVNTR batch root not found: {root}")
    return root


def _truth(root: Path, experiment: str) -> dict[str, str]:
    """Map ``pair_id`` -> mutation class, for mutated samples only."""
    truth: dict[str, str] = {}
    with (root / experiment / "ground_truth.csv").open() as handle:
        for row in csv.DictReader(handle):
            if row["condition"] == "mutated":
                truth[row["pair_id"]] = row["mutation"]
    return truth


def _normal_pair_ids(root: Path, experiment: str) -> list[str]:
    """Return normal-control pair IDs from the simulator's truth table."""
    with (root / experiment / "ground_truth.csv").open() as handle:
        return [row["pair_id"] for row in csv.DictReader(handle) if row["condition"] == "normal"]


def _kestrel_records(path: Path) -> list[dict[str, str]]:
    """Read a Kestrel result TSV, dropping the negative placeholder row."""
    if not path.exists():
        return []
    with path.open() as handle:
        lines = [line for line in handle if not line.startswith("##")]
    return [
        row
        for row in csv.DictReader(lines, delimiter="\t")
        if not (
            row.get("Confidence") == "Negative"
            or row.get("VID") == "Negative"
            or (row.get("Motif") == "None" and "Motifs" not in row)
        )
    ]


def _advntr_records(path: Path) -> list[tuple[str, int | None]]:
    """Read kept adVNTR rows as ``(state, optional supporting reads)``."""
    if not path.exists():
        return []
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    out: list[tuple[str, int | None]] = []
    for row in rows:
        state = (row.get("Variant") or "").strip()
        if not state or state == "Not applicable":
            continue
        try:
            support = int(float(row.get("NumberOfSupportingReads") or ""))
        except (OverflowError, ValueError):
            support = None
        out.append((state, support))
    return out


@lru_cache(maxsize=1)
def _vcf_only() -> tuple[Counter, Counter, int]:
    """Name every mutated sample from the Kestrel record alone.

    Returns:
        tuple[Counter, Counter, int]: correct-per-class, called-per-class, and the
        number of samples seen.
    """
    root = _require_sim()
    correct: Counter = Counter()
    called: Counter = Counter()
    samples = 0
    for experiment in EXPERIMENTS:
        for pair_id, klass in sorted(_truth(root, experiment).items()):
            samples += 1
            records = _kestrel_records(
                root / experiment / "vntyper" / pair_id / "mutated" / "kestrel" / "kestrel_result.tsv"
            )
            for record in records:
                called[klass] += 1
                call = from_kestrel(record["Motifs"], int(record["POS"]), record["REF"], record["ALT"])
                if call.name == EXPECTED_NAME[klass]:
                    correct[klass] += 1
    return correct, called, samples


@dataclass(frozen=True)
class DisplayCounts:
    """Displayed, exact, and wrong names for one policy or tier."""

    displayed: int
    exact: int
    wrong: int


@dataclass(frozen=True)
class PolicyMetrics:
    """Display-aware outcomes for one reconciliation policy."""

    total: DisplayCounts
    by_tier: dict[str, DisplayCounts]
    exact_by_class: Counter[str]
    control_findings: int
    render_mismatches: tuple[tuple[str, str | None, str], ...]


@dataclass(frozen=True)
class GoldenReplay:
    """Both policies plus corpus-shape and BAM-consultation evidence."""

    without_bam: PolicyMetrics
    with_bam: PolicyMetrics
    mutated_samples: int
    normal_samples: int
    bam_eligible_samples: int
    bam_row_fetches: int
    kestrel_records: int
    usable_row_loci: int


@dataclass(frozen=True)
class _SourceEvidence:
    """Translated calls and source-bound support for one sample."""

    vcf_calls: list[Nomenclature]
    advntr_calls: list[Nomenclature]
    supports: dict[str, int | None]


def _lesser(left: int | None, right: int | None) -> int | None:
    """Return the smaller known value, propagating unknown support."""
    if left is None or right is None:
        return None
    return min(left, right)


def _optional_int(value: object) -> int | None:
    """Parse a production result-cell integer without inventing a zero."""
    try:
        return int(float(str(value)))
    except (OverflowError, TypeError, ValueError):
        return None


def _source_evidence(
    kestrel_records: list[dict[str, str]],
    advntr_records: list[tuple[str, int | None]],
) -> _SourceEvidence:
    """Translate caller rows in file order and bind support to its source."""
    vcf_calls = []
    advntr_calls = []
    supports: dict[str, int | None] = {}

    for record in kestrel_records:
        try:
            call = from_kestrel(record["Motifs"], int(record["POS"]), record["REF"], record["ALT"])
        except (KeyError, TypeError, ValueError):
            continue
        vcf_calls.append(call)
        depth = _optional_int(record.get("Estimated_Depth_AlternateVariant"))
        if "kestrel_vcf" not in supports:
            supports["kestrel_vcf"] = depth
        else:
            supports["kestrel_vcf"] = _lesser(supports["kestrel_vcf"], depth)

    for state, state_support in advntr_records:
        row_calls = list(from_advntr(state))
        advntr_calls.extend(row_calls)
        if not row_calls:
            continue
        if "advntr" not in supports:
            supports["advntr"] = state_support
        else:
            supports["advntr"] = _lesser(supports["advntr"], state_support)

    return _SourceEvidence(vcf_calls, advntr_calls, supports)


def _row_locus(record: dict[str, str]) -> tuple[str, int] | None:
    """Return the row's own Kestrel pair locus, if both cells are usable."""
    contig = record.get("Motif_fasta")
    position = record.get("POS_fasta")
    if not contig or position is None:
        return None
    try:
        return contig, int(position)
    except (TypeError, ValueError):
        return None


def _rescue_rows(
    records: list[dict[str, str]],
    kestrel_dir: Path,
    supports: dict[str, int | None],
) -> tuple[list[Nomenclature | None], int]:
    """Rescue each Kestrel row at its own locus with one BAM handle."""
    bam_path = kestrel_dir / "output.bam"
    if not records or not bam_path.is_file():
        return [], 0

    calls = []
    with BamRescuer(bam_path) as rescuer:
        for record in records:
            locus = _row_locus(record)
            if locus is None:
                calls.append(None)
                continue
            contig, position = locus
            consensus = rescuer.rescue(contig, position)
            if consensus is None:
                calls.append(None)
                continue
            call = from_bam(contig, consensus)
            calls.append(call)
            if call is None:
                continue
            support = consensus.supporting_haplotype_records
            existing = supports.get("kestrel_bam")
            supports["kestrel_bam"] = support if existing is None else min(existing, support)
        fetches = rescuer.fetches
    return calls, fetches


def _display_metrics(
    outcomes: list[tuple[str, str, Nomenclature]],
    control_findings: int,
) -> PolicyMetrics:
    """Measure only names whose emitted text equals the selected name."""
    displayed: Counter[str] = Counter()
    exact: Counter[str] = Counter()
    wrong: Counter[str] = Counter()
    exact_by_class: Counter[str] = Counter()
    mismatches: list[tuple[str, str | None, str]] = []

    for pair_id, klass, call in outcomes:
        shown = render(call)
        if call.name is not None and shown != call.name:
            mismatches.append((pair_id, call.name, shown))
        if call.name is None and shown and shown[0].isdigit():
            mismatches.append((pair_id, None, shown))
        is_displayed = call.name is not None and shown == call.name
        if not is_displayed:
            continue
        displayed[call.tier] += 1
        if shown == EXPECTED_NAME[klass]:
            exact[call.tier] += 1
            exact_by_class[klass] += 1
        else:
            wrong[call.tier] += 1

    by_tier = {tier: DisplayCounts(displayed[tier], exact[tier], wrong[tier]) for tier in ("A", "B", "C")}
    total = DisplayCounts(sum(displayed.values()), sum(exact.values()), sum(wrong.values()))
    return PolicyMetrics(total, by_tier, exact_by_class, control_findings, tuple(mismatches))


@lru_cache(maxsize=1)
def _replay() -> GoldenReplay:
    """Replay no-BAM and current unweighted row voting over the loaded cohort."""
    root = _require_sim()
    advntr_root = _require_advntr()
    without_bam_outcomes = []
    with_bam_outcomes = []
    mutated_samples = 0
    normal_samples = 0
    control_findings = 0
    eligible_samples = 0
    row_fetches = 0
    kestrel_record_count = 0
    usable_loci = 0

    for experiment in EXPERIMENTS:
        for pair_id, klass in sorted(_truth(root, experiment).items()):
            sample_dir = root / experiment / "vntyper" / pair_id
            kestrel_dir = sample_dir / "mutated" / "kestrel"
            kestrel_path = kestrel_dir / "kestrel_result.tsv"
            advntr_path = advntr_root / experiment / pair_id / "mutated" / "advntr" / "output_adVNTR_result.tsv"
            mutated_samples += int(kestrel_path.is_file() and advntr_path.is_file())
            records = _kestrel_records(kestrel_path)
            advntr_rows = _advntr_records(advntr_path)
            evidence = _source_evidence(records, advntr_rows)

            kestrel_record_count += len(records)
            usable_loci += sum(_row_locus(record) is not None for record in records)

            no_bam_call = reconcile(
                *evidence.vcf_calls,
                *evidence.advntr_calls,
                supports=dict(evidence.supports),
            )
            without_bam_outcomes.append((pair_id, klass, no_bam_call))

            with_bam_supports = dict(evidence.supports)
            preliminary = reconcile(
                *evidence.vcf_calls,
                *evidence.advntr_calls,
                supports=with_bam_supports,
            )
            eligible = is_candidate(preliminary)
            bam_calls = []
            if eligible:
                eligible_samples += 1
                bam_calls, sample_fetches = _rescue_rows(records, kestrel_dir, with_bam_supports)
                row_fetches += sample_fetches
            named_bam = [call for call in bam_calls if call is not None]
            with_bam_call = reconcile(
                *evidence.vcf_calls,
                *named_bam,
                *evidence.advntr_calls,
                supports=with_bam_supports,
            )
            for bam_call in named_bam:
                with_bam_call = refine(with_bam_call, bam_call)
            with_bam_outcomes.append((pair_id, klass, with_bam_call))

        for pair_id in sorted(_normal_pair_ids(root, experiment)):
            kestrel_path = root / experiment / "vntyper" / pair_id / "normal" / "kestrel" / "kestrel_result.tsv"
            advntr_path = advntr_root / experiment / pair_id / "normal" / "advntr" / "output_adVNTR_result.tsv"
            normal_samples += int(kestrel_path.is_file() and advntr_path.is_file())
            normal_records = _kestrel_records(kestrel_path)
            normal_advntr = _advntr_records(advntr_path)
            normal_evidence = _source_evidence(normal_records, normal_advntr)
            control_findings += int(bool(normal_evidence.vcf_calls or normal_evidence.advntr_calls))

    return GoldenReplay(
        without_bam=_display_metrics(without_bam_outcomes, control_findings),
        with_bam=_display_metrics(with_bam_outcomes, control_findings),
        mutated_samples=mutated_samples,
        normal_samples=normal_samples,
        bam_eligible_samples=eligible_samples,
        bam_row_fetches=row_fetches,
        kestrel_records=kestrel_record_count,
        usable_row_loci=usable_loci,
    )


# ---------------------------------------------------------------------------
# Cohort shape
# ---------------------------------------------------------------------------


def test_the_golden_set_is_two_hundred_samples() -> None:
    _, _, samples = _vcf_only()
    assert samples == 200


def test_the_kestrel_path_makes_one_hundred_and_seventy_eight_calls() -> None:
    """22 of the 200 are no-calls. The 178 is what spec §2.5 measured."""
    _, called, _ = _vcf_only()
    assert sum(called.values()) == 178


# ---------------------------------------------------------------------------
# The two hard gates
# ---------------------------------------------------------------------------


def test_display_metrics_are_pinned_for_both_policies() -> None:
    """Every displayed name is counted, including wrong lower-tier names."""
    replay = _replay()
    assert replay.without_bam.total == DisplayCounts(displayed=138, exact=119, wrong=19)
    assert replay.without_bam.by_tier == {
        "A": DisplayCounts(displayed=53, exact=53, wrong=0),
        "B": DisplayCounts(displayed=85, exact=66, wrong=19),
        "C": DisplayCounts(displayed=0, exact=0, wrong=0),
    }
    assert replay.with_bam.total == DisplayCounts(displayed=154, exact=136, wrong=18)
    assert replay.with_bam.by_tier == {
        "A": DisplayCounts(displayed=53, exact=53, wrong=0),
        "B": DisplayCounts(displayed=101, exact=83, wrong=18),
        "C": DisplayCounts(displayed=0, exact=0, wrong=0),
    }


def test_canonical_dupc_is_ninety_six_of_ninety_six() -> None:
    """Exact, not a floor. 96 of the 100 dupC samples produce a call; all 96 name."""
    correct, called, _ = _vcf_only()
    assert called["dupC"] == 96
    assert correct["dupC"] == 96


def test_no_tier_a_name_is_ever_wrong() -> None:
    """The one hard gate. A confident wrong name fails the build."""
    assert _replay().without_bam.by_tier["A"].wrong == 0


def test_no_tier_a_name_is_wrong_with_bam_rescue_either() -> None:
    """The rescue path must not buy coverage with a confident falsehood."""
    assert _replay().with_bam.by_tier["A"].wrong == 0


def test_the_bam_is_consulted_only_for_a_minority_of_samples() -> None:
    """Pin sample eligibility and row fetches, not an approximate BED locus."""
    replay = _replay()
    assert replay.bam_eligible_samples == 83
    assert replay.bam_row_fetches == 68
    assert replay.bam_row_fetches < 100


def test_no_negative_control_is_ever_named() -> None:
    """All 200 normal samples must stay silent in both callers."""
    replay = _replay()
    assert replay.without_bam.control_findings == 0
    assert replay.with_bam.control_findings == 0


def test_external_corpus_is_loaded_and_every_kestrel_row_has_a_locus() -> None:
    """Prove all samples and row loci entered the replay before asserting metrics."""
    replay = _replay()
    assert replay.mutated_samples == 200
    assert replay.normal_samples == 200
    assert replay.kestrel_records == 178
    assert replay.usable_row_loci == 178


def test_display_counts_partition_into_exact_and_wrong_names() -> None:
    """No lower-tier displayed name can fall outside the exact/wrong denominator."""
    for policy in (_replay().without_bam, _replay().with_bam):
        assert policy.total.displayed == policy.total.exact + policy.total.wrong
        assert sum(counts.displayed for counts in policy.by_tier.values()) == policy.total.displayed
        assert sum(counts.exact for counts in policy.by_tier.values()) == policy.total.exact
        assert sum(counts.wrong for counts in policy.by_tier.values()) == policy.total.wrong
        for counts in policy.by_tier.values():
            assert counts.displayed == counts.exact + counts.wrong


# ---------------------------------------------------------------------------
# Per-class regression floors
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("mutation,floor", sorted(VCF_ONLY_FLOOR.items()))
def test_the_vcf_only_per_class_baseline_does_not_regress(mutation: str, floor: int) -> None:
    correct, _, _ = _vcf_only()
    assert correct[mutation] >= floor, f"{mutation}: {correct[mutation]} < pinned {floor}"


def test_the_vcf_only_total_does_not_regress() -> None:
    """Spec §2.4 measured 123 of the 178 calls."""
    correct, _, _ = _vcf_only()
    assert sum(correct.values()) >= 123


def test_the_current_record_voting_total_is_exact() -> None:
    """Pin all displayed outcomes; current record voting names 154, 136 exactly."""
    assert _replay().with_bam.total == DisplayCounts(displayed=154, exact=136, wrong=18)


def test_the_bam_rescue_recovers_alleles_the_vcf_could_not() -> None:
    """Pin the exact classes gained by production-shaped record voting."""
    replay = _replay()
    expected_policy_two = {
        "delGCCCA": 3,
        "dupA": 6,
        "dupC": 98,
        "insCCCC": 6,
        "insC_pos23": 9,
        "insG": 6,
        "insG_pos54": 7,
        "insG_pos58": 1,
    }
    assert {klass: replay.with_bam.exact_by_class[klass] for klass in EXPECTED_NAME} == {
        klass: expected_policy_two.get(klass, 0) for klass in EXPECTED_NAME
    }
    expected_gains = {"insCCCC": 6, "insG": 4, "insG_pos54": 7}
    gains = {
        klass: replay.with_bam.exact_by_class[klass] - replay.without_bam.exact_by_class[klass]
        for klass in EXPECTED_NAME
    }
    assert gains == {klass: expected_gains.get(klass, 0) for klass in EXPECTED_NAME}


def test_the_emitted_cell_shows_the_name_whenever_one_was_computed() -> None:
    """The display-aware denominator must not hide an internal/render mismatch."""
    replay = _replay()
    assert replay.without_bam.render_mismatches == ()
    assert replay.with_bam.render_mismatches == ()


def test_reconciliation_produces_tier_a_names() -> None:
    """Tier A must be reachable in practice, not just in principle.

    Without this, a change that quietly made tier A unreachable would still pass
    the "no tier-A name is wrong" gate -- vacuously.
    """
    assert _replay().without_bam.by_tier["A"] == DisplayCounts(displayed=53, exact=53, wrong=0)
