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
from functools import lru_cache
from pathlib import Path

import pytest

from vntyper.scripts.nomenclature import from_advntr, from_kestrel, reconcile, render
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


def _kestrel_records(path: Path) -> list[dict[str, str]]:
    """Read a Kestrel result TSV, dropping the negative placeholder row."""
    if not path.exists():
        return []
    with path.open() as handle:
        lines = [line for line in handle if not line.startswith("##")]
    return [
        row
        for row in csv.DictReader(lines, delimiter="\t")
        if row.get("Confidence") != "Negative" and row.get("Motif") != "None"
    ]


def _advntr_records(path: Path) -> list[tuple[str, int]]:
    """Read an adVNTR result TSV as ``(state, supporting_reads)``."""
    if not path.exists():
        return []
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    out: list[tuple[str, int]] = []
    for row in rows:
        state = (row.get("Variant") or "").strip()
        if not state or state == "Not applicable":
            continue
        try:
            support = int(float(row.get("NumberOfSupportingReads") or 0))
        except ValueError:
            support = 0
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


@lru_cache(maxsize=1)
def _reconciled() -> tuple[Counter, Counter, list[tuple[str, str, str]], int]:
    """Name every mutated sample from both callers reconciled.

    Returns:
        tuple: correct-per-class, tier-A-correct-per-class, tier-A names that
        disagree with truth, and the number of negative controls that produced any
        call at all.
    """
    root = _require_sim()
    advntr_root = _require_advntr()

    correct: Counter = Counter()
    tier_a_correct: Counter = Counter()
    wrong: list[tuple[str, str, str]] = []
    controls_called = 0

    for experiment in EXPERIMENTS:
        for pair_id, klass in sorted(_truth(root, experiment).items()):
            sample_dir = root / experiment / "vntyper" / pair_id
            calls = [
                from_kestrel(record["Motifs"], int(record["POS"]), record["REF"], record["ALT"])
                for record in _kestrel_records(sample_dir / "mutated" / "kestrel" / "kestrel_result.tsv")
            ]
            # Support is bound to the source it came from, and an adVNTR agreement
            # is taken at its weakest state: a sample-wide maximum would let a
            # well-covered unrelated observation lend depth to a 1-read agreement.
            supports: dict[str, int | None] = {}
            # Every event adVNTR reported, not just the first: a locus showing three
            # simultaneous events is not one simple allele, and hiding the rest lets
            # a wrong name reach tier A.
            for state, state_support in _advntr_records(
                advntr_root / experiment / pair_id / "mutated" / "advntr" / "output_adVNTR_result.tsv"
            ):
                calls.extend(from_advntr(state))
                seen = supports.get("advntr")
                supports["advntr"] = state_support if seen is None else min(seen, state_support)
            for record in _kestrel_records(sample_dir / "mutated" / "kestrel" / "kestrel_result.tsv"):
                try:
                    supports["kestrel_vcf"] = int(float(record.get("Estimated_Depth_AlternateVariant") or 0))
                except ValueError:
                    supports["kestrel_vcf"] = None

            merged = reconcile(*calls, supports=supports)
            if merged.name == EXPECTED_NAME[klass]:
                correct[klass] += 1
                if merged.tier == "A":
                    tier_a_correct[klass] += 1
            elif merged.tier == "A":
                wrong.append((pair_id, str(merged.name), EXPECTED_NAME[klass]))

            if _kestrel_records(sample_dir / "normal" / "kestrel" / "kestrel_result.tsv") or _advntr_records(
                advntr_root / experiment / pair_id / "normal" / "advntr" / "output_adVNTR_result.tsv"
            ):
                controls_called += 1

    return correct, tier_a_correct, wrong, controls_called


@lru_cache(maxsize=1)
def _hybrid() -> tuple[Counter, Counter, list[tuple[str, str, str]], int]:
    """Name every sample from both callers plus BAM rescue.

    Returns:
        tuple: correct-per-class, tier-A-correct-per-class, tier-A names that
        disagree with truth, and the number of BAM region fetches performed.
    """
    root = _require_sim()
    advntr_root = _require_advntr()

    correct: Counter = Counter()
    tier_a_correct: Counter = Counter()
    wrong: list[tuple[str, str, str]] = []
    fetches = 0

    for experiment in EXPERIMENTS:
        for pair_id, klass in sorted(_truth(root, experiment).items()):
            kestrel_dir = root / experiment / "vntyper" / pair_id / "mutated" / "kestrel"
            calls = [
                from_kestrel(record["Motifs"], int(record["POS"]), record["REF"], record["ALT"])
                for record in _kestrel_records(kestrel_dir / "kestrel_result.tsv")
            ]
            supports: dict[str, int | None] = {}
            for state, state_support in _advntr_records(
                advntr_root / experiment / pair_id / "mutated" / "advntr" / "output_adVNTR_result.tsv"
            ):
                calls.extend(from_advntr(state))
                seen = supports.get("advntr")
                supports["advntr"] = state_support if seen is None else min(seen, state_support)
            for record in _kestrel_records(kestrel_dir / "kestrel_result.tsv"):
                try:
                    supports["kestrel_vcf"] = int(float(record.get("Estimated_Depth_AlternateVariant") or 0))
                except ValueError:
                    supports["kestrel_vcf"] = None

            merged = reconcile(*calls, supports=supports)

            if is_candidate(merged):
                bed = kestrel_dir / "output.bed"
                if bed.is_file():
                    fields = bed.read_text().split()
                    if len(fields) >= 2:
                        contig, position = fields[0], int(fields[1])
                        with BamRescuer(kestrel_dir / "output.bam") as rescuer:
                            consensus = rescuer.rescue(contig, position)
                            fetches += rescuer.fetches
                        if consensus is not None:
                            merged = refine(merged, from_bam(contig, consensus))

            if merged.name == EXPECTED_NAME[klass]:
                correct[klass] += 1
                if merged.tier == "A":
                    tier_a_correct[klass] += 1
            elif merged.tier == "A":
                wrong.append((pair_id, str(merged.name), EXPECTED_NAME[klass]))

    return correct, tier_a_correct, wrong, fetches


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


def test_canonical_dupc_is_ninety_six_of_ninety_six() -> None:
    """Exact, not a floor. 96 of the 100 dupC samples produce a call; all 96 name."""
    correct, called, _ = _vcf_only()
    assert called["dupC"] == 96
    assert correct["dupC"] == 96


def test_no_tier_a_name_is_ever_wrong() -> None:
    """The one hard gate. A confident wrong name fails the build."""
    _, _, wrong, _ = _reconciled()
    assert wrong == [], f"tier-A names disagreeing with truth: {wrong}"


def test_no_tier_a_name_is_wrong_with_bam_rescue_either() -> None:
    """The rescue path must not buy coverage with a confident falsehood."""
    _, _, wrong, _ = _hybrid()
    assert wrong == [], f"tier-A names disagreeing with truth after rescue: {wrong}"


def test_the_bam_is_consulted_only_for_a_minority_of_samples() -> None:
    """`VCF primary, BAM refines`: the common path must not pay for the rescue."""
    _, _, _, fetches = _hybrid()
    assert fetches < 100, f"{fetches} BAM fetches across 200 samples is not a minority"


def test_no_negative_control_is_ever_named() -> None:
    """All 200 normal samples must stay silent in both callers."""
    _, _, _, controls_called = _reconciled()
    assert controls_called == 0


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


def test_the_hybrid_total_does_not_regress() -> None:
    """Kestrel VCF + adVNTR + BAM rescue, one uniform policy, no truth consulted.

    Spec §4.3 asks for 136. That figure is the per-class **maximum** of the VCF-only
    and BAM-only columns in §2.6 -- 96+9+9+6+3+7+4+2 -- so reaching it requires
    knowing, per class, which of the two methods to believe. A policy that has to
    decide without truth cannot: preferring the reads costs more than it gains
    wherever the VCF was already right. 129 is what one uniform rule achieves.
    """
    correct, _, _, _ = _hybrid()
    assert sum(correct.values()) >= 129


def test_the_bam_rescue_recovers_alleles_the_vcf_could_not() -> None:
    """The rescue path must earn its cost: classes the VCF scores 0 on must gain."""
    correct, _, _, _ = _hybrid()
    vcf_correct, _, _ = _vcf_only()
    gained = {klass for klass in EXPECTED_NAME if correct[klass] > vcf_correct[klass]}
    assert gained, "BAM rescue changed nothing; the whole path is then dead weight"


def test_every_tier_a_cell_renders_as_the_bare_name_and_no_other_tier_does() -> None:
    """Assert the *emitted* cell, not the internal field.

    The golden numbers are computed from `Nomenclature.name`, which every tier
    carries. What reaches a user is `render`. Without this, replacing
    ``render(call)`` with ``call.name or ""`` in the production adapter would leak a
    bare number on every tier-B row while this whole tier stayed green.
    """
    root = _require_sim()
    advntr_root = _require_advntr()
    leaked: list[tuple[str, str]] = []

    for experiment in EXPERIMENTS:
        for pair_id, _ in sorted(_truth(root, experiment).items()):
            calls = [
                from_kestrel(record["Motifs"], int(record["POS"]), record["REF"], record["ALT"])
                for record in _kestrel_records(
                    root / experiment / "vntyper" / pair_id / "mutated" / "kestrel" / "kestrel_result.tsv"
                )
            ]
            supports: dict[str, int | None] = {}
            for state, state_support in _advntr_records(
                advntr_root / experiment / pair_id / "mutated" / "advntr" / "output_adVNTR_result.tsv"
            ):
                calls.extend(from_advntr(state))
                seen = supports.get("advntr")
                supports["advntr"] = state_support if seen is None else min(seen, state_support)
            if not calls:
                continue

            merged = reconcile(*calls, supports=supports)
            shown = render(merged)
            if merged.tier == "A":
                if shown != merged.name:
                    leaked.append((pair_id, shown))
            elif merged.name is not None and shown == merged.name:
                leaked.append((pair_id, shown))

    assert leaked == [], f"rendered cell disagreed with what the tier permits: {leaked[:5]}"


def test_reconciliation_produces_tier_a_names() -> None:
    """Tier A must be reachable in practice, not just in principle.

    Without this, a change that quietly made tier A unreachable would still pass
    the "no tier-A name is wrong" gate -- vacuously.
    """
    _, tier_a_correct, _, _ = _reconciled()
    assert sum(tier_a_correct.values()) >= 40
