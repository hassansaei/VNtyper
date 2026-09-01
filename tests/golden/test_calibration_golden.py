"""Development-only calibration replay and external-evidence blocker."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pytest

from tests.golden import calibration_oracle
from tests.golden.calibration_oracle import load_development_snapshot
from tests.golden.identity_oracle import DisplayCounts

pytestmark = pytest.mark.golden

_SIM_ROOT = os.environ.get("VNTYPER_SIM_ROOT")
_ADVNTR_ROOT = os.environ.get("VNTYPER_ADVNTR_ROOT")
if not _SIM_ROOT or not _ADVNTR_ROOT:
    pytest.skip("VNTYPER_SIM_ROOT and VNTYPER_ADVNTR_ROOT benchmark roots are unset", allow_module_level=True)

SIM_ROOT = Path(_SIM_ROOT)
ADVNTR_ROOT = Path(_ADVNTR_ROOT)
REPO_ROOT = Path(__file__).parents[2]
PACKAGED_PROFILE_SHA256 = "be6329fb12107a1b6b65e425257be6233c7e2115e299e941c12a63a6a6d59718"
PACKAGED_PROJECTION_SHA256 = "338fe05d010f623e794dcf93393904fa13bd8713e2d074c8a5b6c72d6efd96fd"


@pytest.fixture(scope="module")
def snapshot():
    """Load both explicit roots with no per-test skip fallback."""
    return load_development_snapshot(SIM_ROOT, ADVNTR_ROOT)


def test_calibration_oracle_has_no_production_import_path() -> None:
    scanned = calibration_oracle.assert_independent_import_closure(REPO_ROOT)

    assert Path(calibration_oracle.__file__).resolve() in scanned


def test_both_roots_complete_population_and_row_locus_counts_are_loaded(snapshot) -> None:
    assert snapshot.sim_root == SIM_ROOT.resolve()
    assert snapshot.advntr_root == ADVNTR_ROOT.resolve()
    assert snapshot.mutated_samples == 200
    assert snapshot.control_samples == 200
    assert snapshot.public_identity_rows == 374
    assert snapshot.selected_locus_rows == 178


def test_shipped_projection_is_reproduced_before_any_candidate_claim(snapshot) -> None:
    assert snapshot.total == DisplayCounts(displayed=154, exact=136, wrong=18)
    assert snapshot.by_tier == {
        "A": DisplayCounts(displayed=53, exact=53, wrong=0),
        "B": DisplayCounts(displayed=101, exact=83, wrong=18),
        "C": DisplayCounts(displayed=0, exact=0, wrong=0),
    }
    assert snapshot.control_findings == 0

    profile = REPO_ROOT / "vntyper" / "profiles" / "decision_profile.json"
    projection = REPO_ROOT / "vntyper" / "profiles" / "decision_projection.json"
    assert hashlib.sha256(profile.read_bytes()).hexdigest() == PACKAGED_PROFILE_SHA256
    assert hashlib.sha256(projection.read_bytes()).hexdigest() == PACKAGED_PROJECTION_SHA256


def test_examined_simulations_are_not_misrepresented_as_locked_heldout(snapshot) -> None:
    assert snapshot.evidence_role == "previously-examined-development-simulation"
    assert snapshot.eligible_for_locked_evaluate is False
    assert "previously examined development evidence" in snapshot.ineligibility_reason
    assert "neither independent external validation nor a custodian-locked held-out cohort" in (
        snapshot.ineligibility_reason
    )
