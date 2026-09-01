"""Independent development-corpus snapshot for calibration golden checks.

This module may use the separately guarded truth oracle and Python's standard
library only. It must never import VNtyper decisions, canonicalizers, codecs,
grouping rules, calibration predicates, or profile resolution.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from tests.golden import identity_oracle
from tests.golden.identity_oracle import DisplayCounts, GoldenCorpus


@dataclass(frozen=True)
class DevelopmentCalibrationSnapshot:
    """Literal pre-fit projection and truthful evidence eligibility state."""

    sim_root: Path
    advntr_root: Path
    mutated_samples: int
    control_samples: int
    public_identity_rows: int
    selected_locus_rows: int
    total: DisplayCounts
    by_tier: dict[str, DisplayCounts]
    control_findings: int
    evidence_role: str
    eligible_for_locked_evaluate: bool
    ineligibility_reason: str


def load_development_snapshot(sim_root: Path, advntr_root: Path) -> DevelopmentCalibrationSnapshot:
    """Load both roots and freeze the exact shipped pre-fit projection."""
    corpus = identity_oracle.load_golden_corpus(sim_root, advntr_root)
    return snapshot_from_corpus(corpus)


def snapshot_from_corpus(corpus: GoldenCorpus) -> DevelopmentCalibrationSnapshot:
    """Project only independently observed corpus facts, never fitted decisions."""
    if not isinstance(corpus, GoldenCorpus):
        raise AssertionError("calibration oracle requires an independently loaded GoldenCorpus")
    public_rows = sum(len(rows) for rows in corpus.identity_projection_by_sample.values())
    selected_rows = sum(len(rows) for rows in corpus.selected_projection_by_sample.values())
    return DevelopmentCalibrationSnapshot(
        corpus.sim_root,
        corpus.advntr_root,
        corpus.mutated_samples,
        corpus.control_samples,
        public_rows,
        selected_rows,
        corpus.total,
        dict(corpus.by_tier),
        corpus.control_findings,
        "previously-examined-development-simulation",
        False,
        (
            "The simulations and paired caller outputs were previously examined development evidence; "
            "they are neither independent external validation nor a custodian-locked held-out cohort."
        ),
    )


def assert_independent_import_closure(repo_root: Path) -> frozenset[Path]:
    """Reject direct, recursive, or dynamic production imports in this oracle."""
    scanned = identity_oracle.assert_independent_import_closure(Path(__file__), repo_root)
    return frozenset(scanned)
