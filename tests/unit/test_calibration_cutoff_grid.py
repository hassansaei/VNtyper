"""Automatic component grids exercise actual production cutoff decisions."""

from copy import deepcopy
from dataclasses import replace

import pytest

from tests.builders import kestrel_config
from tests.unit.test_calibration_caller_policy import policy_document, policy_values
from tests.unit.test_calibration_kestrel_replay import _candidate, _capture, _raw
from vntyper.scripts.calibration_caller_policy import decode_caller_policy_values

pytestmark = pytest.mark.unit

FLOOR = "/components/kestrel/confidence_assignment/reporting_floor"
LOW = "/components/kestrel/confidence_assignment/depth_score_thresholds/low"
HIGH = "/components/kestrel/confidence_assignment/depth_score_thresholds/high"
GG = "/components/kestrel/alt_filtering/gg_depth_score_threshold"
MODE = "/components/advntr/calibrated_calling/mode"


def custom(*, dual=True):
    doc = {
        "schema_version": "calibration-cutoff-grid-v1",
        "kestrel": {"reporting_floor": [0.002, 0.008], "gg_depth_score_threshold": [0.002]},
    }
    if dual:
        doc["advntr"] = {"cutoff": [0.0001, 0.005], "minimum_read_support": [2, 5]}
    return doc


def test_default_component_grids_preserve_baseline_and_joint_low_boundary():
    from vntyper.scripts.calibration_cutoff_grid import build_cutoff_grid

    raw = policy_document()
    policy_values(raw)[FLOOR] = 0.00469
    policy_values(raw)[MODE] = "legacy"
    base = decode_caller_policy_values(raw)
    grid = build_cutoff_grid(base)
    assert len(grid.kestrel) == 55 and len(grid.advntr) == 24
    assert grid.kestrel[0].candidate_id == grid.advntr[0].candidate_id == "baseline"
    assert grid.kestrel[0].policy is base and grid.advntr[0].policy is base
    for candidate in grid.kestrel:
        assert candidate.policy.values[LOW] == min(base.values[LOW], candidate.policy.values[FLOOR])
        assert candidate.policy.values[HIGH] == base.values[HIGH]
        assert set(candidate.parameters) <= {FLOOR, LOW, GG}
    assert len({c.policy.sha256 for c in grid.kestrel}) == len(grid.kestrel)
    with pytest.raises(TypeError):
        grid.kestrel[1].parameters[FLOOR] = 1


def test_joint_floor_relaxation_changes_production_disposition_and_floor_alone_does_not():
    from vntyper.scripts.calibration_cutoff_grid import build_cutoff_grid
    from vntyper.scripts.calibration_kestrel_replay import replay_kestrel_capture

    capture = _capture(_raw(depth_alt=10, depth_region=3333), kestrel_config())
    grid = build_cutoff_grid(capture.baseline_policy, custom(dual=False))
    joint = next(c.policy for c in grid.kestrel if c.policy.values[FLOOR] == 0.002)
    floor_only = _candidate(joint, **{LOW: capture.baseline_policy.values[LOW]})

    def replay(policy):
        return replay_kestrel_capture(capture, policy, capture_policy_sha256=capture.provenance.capture_policy_sha256)

    assert replay(capture.baseline_policy).disposition == "no-call"
    assert replay(floor_only).disposition == "no-call"
    assert replay(joint).disposition == "called"


def test_stricter_candidate_removes_weaker_finding_while_retaining_stronger_finding():
    from vntyper.scripts.calibration_cutoff_grid import build_cutoff_grid
    from vntyper.scripts.calibration_kestrel_replay import replay_kestrel_capture

    captures = [_capture(_raw(depth_alt=10, depth_region=depth), kestrel_config()) for depth in (2000, 1000)]
    grid = build_cutoff_grid(captures[0].baseline_policy, custom(dual=False))
    tight = next(c.policy for c in grid.kestrel if c.policy.values[FLOOR] == 0.008)
    for capture, expected in zip(captures, ("no-call", "called"), strict=True):
        commitment = capture.provenance.capture_policy_sha256
        assert (
            replay_kestrel_capture(capture, capture.baseline_policy, capture_policy_sha256=commitment).disposition
            == "called"
        )
        assert replay_kestrel_capture(capture, tight, capture_policy_sha256=commitment).disposition == expected


def test_custom_exact_baseline_is_retained_and_ad_alternatives_are_legacy():
    from vntyper.scripts.calibration_cutoff_grid import build_cutoff_grid

    base = decode_caller_policy_values(policy_document())
    grid = build_cutoff_grid(base, custom())
    assert len(grid.kestrel) == 3 and len(grid.advntr) == 5
    assert grid.advntr[0].policy.values[MODE] == "exact"
    assert all(c.policy.values[MODE] == "legacy" for c in grid.advntr[1:])
    assert (
        build_cutoff_grid(decode_caller_policy_values(policy_document(include_advntr=False)), custom(dual=False)).advntr
        == ()
    )


@pytest.mark.parametrize(
    "field,value",
    [
        ("reporting_floor", []),
        ("reporting_floor", [True]),
        ("reporting_floor", [float("nan")]),
        ("reporting_floor", [-0.1]),
        ("gg_depth_score_threshold", [1.1]),
        ("gg_depth_score_threshold", ["0.1"]),
    ],
)
def test_malformed_kestrel_grid_values_are_refused(field, value):
    from vntyper.scripts.calibration_cutoff_grid import build_cutoff_grid

    raw = custom()
    raw["kestrel"][field] = value
    with pytest.raises(ValueError):
        build_cutoff_grid(decode_caller_policy_values(policy_document()), raw)


@pytest.mark.parametrize(
    "field,value",
    [
        ("cutoff", [0]),
        ("cutoff", [1]),
        ("cutoff", [False]),
        ("minimum_read_support", [0]),
        ("minimum_read_support", [2.0]),
        ("minimum_read_support", [True]),
    ],
)
def test_malformed_advntr_grid_values_are_refused(field, value):
    from vntyper.scripts.calibration_cutoff_grid import build_cutoff_grid

    raw = custom()
    raw["advntr"][field] = value
    with pytest.raises(ValueError):
        build_cutoff_grid(decode_caller_policy_values(policy_document()), raw)


@pytest.mark.parametrize("cap", [True, 0, -1, 2.0, 14])
def test_candidate_product_cap_is_enforced(cap):
    from vntyper.scripts.calibration_cutoff_grid import build_cutoff_grid

    with pytest.raises(ValueError):
        build_cutoff_grid(decode_caller_policy_values(policy_document()), custom(), max_candidates=cap)


def test_custom_schema_is_closed_and_baseline_identity_is_revalidated():
    from vntyper.scripts.calibration_cutoff_grid import build_cutoff_grid

    base = decode_caller_policy_values(policy_document())
    for doc in ({}, {**custom(), "extra": 1}, {**custom(), "schema_version": "unknown"}, custom(dual=False)):
        with pytest.raises(ValueError):
            build_cutoff_grid(base, doc)
    with pytest.raises(ValueError):
        build_cutoff_grid(replace(base, sha256="f" * 64))
    doc = custom()
    doc["kestrel"]["extra"] = []
    with pytest.raises(ValueError):
        build_cutoff_grid(base, doc)


def test_duplicate_unsorted_custom_values_collapse_to_reproducible_policies():
    from vntyper.scripts.calibration_cutoff_grid import build_cutoff_grid

    base = decode_caller_policy_values(policy_document())
    raw = custom()
    duplicate = deepcopy(raw)
    duplicate["kestrel"]["reporting_floor"] = [0.008, 0.002, 0.002]
    a, b = build_cutoff_grid(base, raw), build_cutoff_grid(base, duplicate)
    assert [(c.candidate_id, c.policy.sha256) for c in a.kestrel] == [
        (c.candidate_id, c.policy.sha256) for c in b.kestrel
    ]


def test_huge_fraction_and_product_are_refused_before_materialization():
    from vntyper.scripts.calibration_cutoff_grid import build_cutoff_grid

    base = decode_caller_policy_values(policy_document())
    raw = custom()
    raw["kestrel"]["reporting_floor"] = [10**400]
    with pytest.raises(ValueError):
        build_cutoff_grid(base, raw)
    with pytest.raises(ValueError, match="product"):
        build_cutoff_grid(base, custom(), max_candidates=7)
