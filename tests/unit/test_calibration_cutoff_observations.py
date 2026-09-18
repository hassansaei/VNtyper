"""Replayed cutoff endpoints become truth-bound arms without inventing negatives."""

import json
from dataclasses import replace
from pathlib import Path
from types import MappingProxyType

import pytest

from vntyper.modules.advntr.advntr_calibration_policy import AdvntrCapabilities
from vntyper.scripts.calibration_caller_metrics import CallerObservation, calculate_caller_metrics
from vntyper.scripts.calibration_cohort_manifest import CohortSample
from vntyper.scripts.calibration_cutoff_advntr import (
    AdvntrCutoffGridResult,
    AdvntrCutoffPolicyResult,
    AdvntrCutoffSample,
)
from vntyper.scripts.calibration_cutoff_kestrel import KestrelGridObservation, KestrelGridReplay

pytestmark = pytest.mark.unit

DIGEST = "0" * 64
RIGHT = "MUC1-X-60-coding-v1|60|59|-|C"
WRONG = "MUC1-X-60-coding-v1|60|59|-|A"
CAPABILITIES = AdvntrCapabilities("2.3.0", "build", None, (), (2,), ("advntr-frameshift-policy-v1",), (), DIGEST)


def sample(key, genotype=True, truth_variant=None, group=None):
    return CohortSample(
        key,
        Path("/nonexistent") / f"{key}.bam",
        "GRCh38",
        genotype,
        None,
        group or f"sample:{key}",
        None,
        None,
        None,
        truth_variant,
    )


def endpoint(key, policy_id="baseline", *, call=True, identity=RIGHT, confidence="High_Precision", flag="Not flagged"):
    """Build one replayed Kestrel endpoint from the call it is meant to encode."""
    disposition = {True: "called", False: "no-call", None: "unassessable-no-candidates"}[call]
    called = disposition == "called"
    return KestrelGridObservation(
        key,
        policy_id,
        DIGEST,
        DIGEST,
        DIGEST,
        DIGEST,
        disposition,
        call,
        confidence if called else None,
        flag if called else None,
        identity if called else None,
        "translated" if called else None,
        None,
        False if called else None,
    )


def replay(arms):
    keys = tuple(sorted({row.key for rows in arms.values() for row in rows}))
    policy_ids = tuple(arms)
    return KestrelGridReplay(
        policy_ids,
        keys,
        MappingProxyType(
            {name: MappingProxyType({row.key: row for row in rows}) for name, rows in arms.items()},
        ),
        MappingProxyType(dict.fromkeys(policy_ids, DIGEST)),
        MappingProxyType(dict.fromkeys(policy_ids, DIGEST)),
        MappingProxyType(dict.fromkeys(keys, DIGEST)),
        MappingProxyType(dict.fromkeys(keys, None)),
        MappingProxyType(dict.fromkeys(keys, "capture-replay-authoritative")),
        DIGEST,
        1,
        DIGEST,
        DIGEST,
    )


def advntr(arms):
    policies = tuple(
        AdvntrCutoffPolicyResult(
            name,
            DIGEST,
            "advntr-" + DIGEST,
            tuple(AdvntrCutoffSample(key, assessable, call, ()) for key, assessable, call in rows),
        )
        for name, rows in arms.items()
    )
    return AdvntrCutoffGridResult(Path("/nonexistent"), "baseline", DIGEST, CAPABILITIES, policies, DIGEST)


def test_every_kestrel_disposition_maps_to_its_scientific_call():
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms

    rows = (endpoint("a", call=True), endpoint("b", call=False), endpoint("c", call=None))
    samples = tuple(sample(key) for key in ("a", "b", "c"))
    arm = kestrel_observation_arms(replay({"baseline": rows}), samples)["baseline"]
    assert [row.called_positive for row in arm] == [True, False, None]
    assert [row.called_variants for row in arm] == [(RIGHT,), (), ()]


def test_native_negative_placeholder_makes_an_empty_capture_a_real_negative():
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms

    proven = KestrelGridObservation(
        "a", "baseline", DIGEST, DIGEST, DIGEST, DIGEST, "unassessable-no-candidates", False, *(None,) * 6
    )
    arm = kestrel_observation_arms(replay({"baseline": (proven,)}), (sample("a"),))["baseline"]
    assert arm[0].called_positive is False


@pytest.mark.parametrize(
    "disposition,call",
    [
        ("called", False),
        ("called", None),
        ("no-call", None),
        ("no-call", True),
        ("unassessable-no-candidates", True),
        ("invented", None),
    ],
)
def test_kestrel_disposition_contradicting_its_call_is_refused(disposition, call):
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms

    row = KestrelGridObservation("a", "baseline", DIGEST, DIGEST, DIGEST, DIGEST, disposition, call, *(None,) * 6)
    with pytest.raises(ValueError):
        kestrel_observation_arms(replay({"baseline": (row,)}), (sample("a"),))


def test_called_endpoint_without_a_canonical_identity_keeps_the_call():
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms

    row = endpoint("a", identity=None)
    arm = kestrel_observation_arms(replay({"baseline": (row,)}), (sample("a"),))["baseline"]
    assert arm[0].called_positive is True and arm[0].called_variants == ()


def test_tier_a_stays_empty_because_replay_evidence_never_carries_the_tier():
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms

    arm = kestrel_observation_arms(replay({"baseline": (endpoint("a"),)}), (sample("a"),))["baseline"]
    assert arm[0].tier_a_variants == ()


def test_advntr_unassessable_evidence_never_becomes_a_negative():
    from vntyper.scripts.calibration_cutoff_observations import advntr_observation_arms

    result = advntr({"baseline": (("a", False, None), ("b", True, False), ("c", True, True))})
    arm = advntr_observation_arms(result, tuple(sample(key) for key in ("a", "b", "c")))["baseline"]
    assert [row.called_positive for row in arm] == [None, False, True]
    assert all(row.called_variants == () for row in arm)


@pytest.mark.parametrize("rows", [(("a", True, None),), (("a", False, True),), (("a", False, False),)])
def test_advntr_assessability_contradicting_its_call_is_refused(rows):
    from vntyper.scripts.calibration_cutoff_observations import advntr_observation_arms

    with pytest.raises(ValueError):
        advntr_observation_arms(advntr({"baseline": rows}), (sample("a"),))


def test_unknown_truth_survives_into_the_unknown_denominator():
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms

    samples = (sample("a", genotype=None), sample("b", genotype=False), sample("c", genotype=True))
    rows = tuple(endpoint(key) for key in ("a", "b", "c"))
    arm = kestrel_observation_arms(replay({"baseline": rows}), samples)["baseline"]
    assert [row.truth_positive for row in arm] == [None, False, True]
    metrics = calculate_caller_metrics(arm)
    assert metrics.unknown_truth_count == 1 and metrics.known_truth_count == 2


def test_truth_variants_distinguish_unavailable_identity_from_a_known_negative():
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms

    samples = (
        sample("a", genotype=True, truth_variant=RIGHT),
        sample("b", genotype=True),
        sample("c", genotype=False),
        sample("d", genotype=None),
    )
    rows = tuple(endpoint(key) for key in ("a", "b", "c", "d"))
    arm = kestrel_observation_arms(replay({"baseline": rows}), samples)["baseline"]
    assert [row.truth_variants for row in arm] == [(RIGHT,), None, (), None]


@pytest.mark.parametrize("genotype", [False, None])
def test_a_confirmed_variant_on_a_nonpositive_sample_is_refused(genotype):
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms

    samples = (sample("a", genotype=genotype, truth_variant=RIGHT),)
    with pytest.raises(ValueError, match="truth_variant"):
        kestrel_observation_arms(replay({"baseline": (endpoint("a"),)}), samples)


def test_exact_identity_recovery_separates_a_right_call_from_a_wrong_one():
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms

    samples = (sample("a", truth_variant=RIGHT), sample("b", truth_variant=RIGHT))
    rows = (endpoint("a", identity=RIGHT), endpoint("b", identity=WRONG))
    arm = kestrel_observation_arms(replay({"baseline": rows}), samples)["baseline"]
    metrics = calculate_caller_metrics(arm)
    assert metrics.true_positives == 2
    assert metrics.exact_variant_recovery.events == 1 and metrics.exact_variant_recovery.total == 2
    assert metrics.wrong_identity_groups == 1


def test_a_roster_sample_absent_from_the_replay_is_refused():
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms

    with pytest.raises(ValueError, match="roster mismatch"):
        kestrel_observation_arms(replay({"baseline": (endpoint("a"),)}), (sample("a"), sample("b")))


def test_a_replayed_sample_absent_from_the_roster_is_refused():
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms

    rows = (endpoint("a"), endpoint("b"))
    with pytest.raises(ValueError, match="roster mismatch"):
        kestrel_observation_arms(replay({"baseline": rows}), (sample("a"),))


def test_advntr_roster_mismatch_is_refused():
    from vntyper.scripts.calibration_cutoff_observations import advntr_observation_arms

    with pytest.raises(ValueError, match="roster mismatch"):
        advntr_observation_arms(advntr({"baseline": (("a", True, True),)}), (sample("a"), sample("b")))


@pytest.mark.parametrize("samples", [(), ("not a sample",), (sample("a"), sample("a"))])
def test_a_malformed_roster_is_refused(samples):
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms

    with pytest.raises(ValueError):
        kestrel_observation_arms(replay({"baseline": (endpoint("a"),)}), samples)


@pytest.mark.parametrize("name", ["kestrel_observation_arms", "advntr_observation_arms"])
def test_foreign_evidence_objects_are_refused(name):
    from vntyper.scripts import calibration_cutoff_observations as module

    with pytest.raises(ValueError):
        getattr(module, name)("not replayed evidence", (sample("a"),))


def test_endpoint_tallies_refuse_foreign_evidence():
    from vntyper.scripts.calibration_cutoff_observations import endpoint_tallies

    with pytest.raises(ValueError):
        endpoint_tallies("not a replay")


def test_a_foreign_kestrel_endpoint_is_refused():
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms

    forged = replace(replay({"baseline": (endpoint("a"),)}), observations=MappingProxyType({"baseline": {"a": "row"}}))
    with pytest.raises(ValueError, match="KestrelGridObservation"):
        kestrel_observation_arms(forged, (sample("a"),))


def test_a_foreign_advntr_outcome_is_refused():
    from vntyper.scripts.calibration_cutoff_observations import advntr_observation_arms

    policies = (AdvntrCutoffPolicyResult("baseline", DIGEST, "advntr-" + DIGEST, (endpoint("a"),)),)
    result = AdvntrCutoffGridResult(Path("/nonexistent"), "baseline", DIGEST, CAPABILITIES, policies, DIGEST)
    with pytest.raises(ValueError, match="AdvntrCutoffSample"):
        advntr_observation_arms(result, (sample("a"),))


@pytest.mark.parametrize("call", ["kestrel_observation_arms", "endpoint_tallies"])
def test_a_declared_policy_without_endpoints_is_refused(call):
    from vntyper.scripts import calibration_cutoff_observations as module

    forged = replace(replay({"baseline": (endpoint("a"),)}), policy_ids=("baseline", "ghost"))
    with pytest.raises(ValueError, match="no endpoints for policy ghost"):
        getattr(module, call)(forged) if call == "endpoint_tallies" else getattr(module, call)(forged, (sample("a"),))


def test_a_repeated_advntr_policy_identity_is_refused():
    from vntyper.scripts.calibration_cutoff_observations import advntr_observation_arms

    policy = AdvntrCutoffPolicyResult(
        "baseline", DIGEST, "advntr-" + DIGEST, (AdvntrCutoffSample("a", True, True, ()),)
    )
    result = AdvntrCutoffGridResult(Path("/nonexistent"), "baseline", DIGEST, CAPABILITIES, (policy, policy), DIGEST)
    with pytest.raises(ValueError, match="unique AdvntrCutoffPolicyResult"):
        advntr_observation_arms(result, (sample("a"),))


def test_a_repeated_advntr_sample_within_one_policy_is_refused():
    from vntyper.scripts.calibration_cutoff_observations import advntr_observation_arms

    with pytest.raises(ValueError, match="repeats a sample"):
        advntr_observation_arms(advntr({"baseline": (("a", True, True), ("a", True, False))}), (sample("a"),))


def test_every_arm_agrees_on_key_group_and_truth_and_feeds_the_evaluator():
    from vntyper.scripts.calibration_cutoff_evaluation import evaluate_cutoff_arms
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms
    from vntyper.scripts.calibration_cutoff_selection import SearchSpec

    samples = tuple(sample(key, genotype=key in {"a", "b"}) for key in ("a", "b", "c", "d"))
    arms = kestrel_observation_arms(
        replay(
            {
                "baseline": (endpoint("a"), endpoint("b", call=False), endpoint("c", call=False), endpoint("d")),
                "relaxed": tuple(endpoint(key, "relaxed", call=key in {"a", "b"}) for key in ("a", "b", "c", "d")),
            }
        ),
        samples,
    )
    shared = {
        name: tuple((r.key, r.group_key, r.truth_positive, r.truth_variants) for r in rows)
        for name, rows in arms.items()
    }
    assert len(set(shared.values())) == 1
    evaluation = evaluate_cutoff_arms(arms, spec=SearchSpec("balanced-accuracy"), folds=2, seed=7)
    assert evaluation["full_data_operating_points"]["relaxed"]["true_positives"] == 2


def test_arm_ordering_is_deterministic_regardless_of_roster_order():
    from vntyper.scripts.calibration_cutoff_observations import kestrel_observation_arms

    rows = tuple(endpoint(key) for key in ("a", "b", "c"))
    samples = tuple(sample(key) for key in ("a", "b", "c"))
    first = kestrel_observation_arms(replay({"baseline": rows}), samples)
    second = kestrel_observation_arms(replay({"baseline": tuple(reversed(rows))}), tuple(reversed(samples)))
    assert first == second
    assert [row.key for row in first["baseline"]] == ["a", "b", "c"]


def union_inputs():
    from vntyper.scripts.calibration_cutoff_observations import advntr_observation_arms, kestrel_observation_arms

    calls = [(left, right) for left in (None, False, True) for right in (None, False, True)]
    keys = tuple(f"s{index}" for index in range(len(calls)))
    samples = tuple(sample(key) for key in keys)
    kestrel = kestrel_observation_arms(
        replay({"baseline": tuple(endpoint(key, call=left) for key, (left, _) in zip(keys, calls, strict=True))}),
        samples,
    )
    native = advntr_observation_arms(
        advntr({"legacy": tuple((key, right is not None, right) for key, (_, right) in zip(keys, calls, strict=True))}),
        samples,
    )
    return calls, keys, kestrel, native


def test_union_propagates_unknown_across_all_nine_call_combinations():
    from vntyper.scripts.calibration_cutoff_observations import union_observation_arms

    calls, keys, kestrel, native = union_inputs()
    arm = union_observation_arms(kestrel, native)["baseline+legacy"]
    observed = {row.key: row.called_positive for row in arm}
    expected = {
        key: True if True in pair else (None if None in pair else False) for key, pair in zip(keys, calls, strict=True)
    }
    assert observed == expected
    identities = {row.key: row.called_variants for row in arm}
    assert identities["s6"] == identities["s7"] == identities["s8"] == (RIGHT,)
    assert identities["s2"] == identities["s5"] == ()


def test_union_arm_ids_name_both_components_and_are_deterministic():
    from vntyper.scripts.calibration_cutoff_observations import (
        advntr_observation_arms,
        kestrel_observation_arms,
        union_observation_arms,
    )

    samples = (sample("a"),)
    kestrel = kestrel_observation_arms(
        replay({"baseline": (endpoint("a"),), "relaxed": (endpoint("a", "relaxed"),)}), samples
    )
    native = advntr_observation_arms(advntr({"legacy": (("a", True, True),), "strict": (("a", True, False),)}), samples)
    arms = union_observation_arms(kestrel, native)
    assert list(arms) == ["baseline+legacy", "baseline+strict", "relaxed+legacy", "relaxed+strict"]
    assert arms == union_observation_arms(kestrel, native)


def test_union_refuses_arms_built_over_different_rosters():
    from vntyper.scripts.calibration_cutoff_observations import (
        advntr_observation_arms,
        kestrel_observation_arms,
        union_observation_arms,
    )

    kestrel = kestrel_observation_arms(replay({"baseline": (endpoint("a"),)}), (sample("a"),))
    native = advntr_observation_arms(advntr({"legacy": (("b", True, True),)}), (sample("b"),))
    with pytest.raises(ValueError, match="identical"):
        union_observation_arms(kestrel, native)


def test_union_refuses_component_ids_that_collide_after_joining():
    from vntyper.scripts.calibration_cutoff_observations import union_observation_arms

    rows = (CallerObservation("a", "sample:a", True, None, True, (RIGHT,), ()),)
    with pytest.raises(ValueError, match="unique"):
        union_observation_arms({"x": rows, "x+y": rows}, {"y": rows, "y+y": rows})


@pytest.mark.parametrize("arms", ["not a mapping", {}, {"": ()}, {" x": ()}])
def test_union_refuses_malformed_arm_inventories(arms):
    from vntyper.scripts.calibration_cutoff_observations import union_observation_arms

    rows = (CallerObservation("a", "sample:a", True, None, True, (RIGHT,), ()),)
    with pytest.raises(ValueError):
        union_observation_arms(arms, {"legacy": rows})


def test_endpoint_tallies_count_confidence_flags_and_dispositions_per_policy():
    from vntyper.scripts.calibration_cutoff_observations import endpoint_tallies

    rows = (
        endpoint("a", confidence="High_Precision", flag="Not flagged"),
        endpoint("b", confidence="Low_Precision", flag="Low_Precision_Depth"),
        endpoint("c", confidence="Low_Precision", flag="Not flagged"),
        endpoint("d", call=False),
        endpoint("e", call=None),
    )
    tallies = endpoint_tallies(
        replay({"baseline": rows, "relaxed": tuple(endpoint(row.key, "relaxed") for row in rows)})
    )
    assert dict(tallies.confidence["baseline"]) == {"High_Precision": 1, "Low_Precision": 2}
    assert dict(tallies.flags["baseline"]) == {"Low_Precision_Depth": 1, "Not flagged": 2}
    assert dict(tallies.dispositions["baseline"]) == {
        "called": 3,
        "no-call": 1,
        "unassessable-no-candidates": 1,
    }
    assert dict(tallies.dispositions["relaxed"]) == {"called": 5}
    with pytest.raises(TypeError):
        tallies.confidence["baseline"]["High_Precision"] = 2


def test_observation_arms_document_round_trips_through_json():
    from vntyper.scripts.calibration_cutoff_observations import (
        kestrel_observation_arms,
        observation_arms_document,
    )

    samples = (sample("a", truth_variant=RIGHT), sample("b", genotype=False))
    arms = kestrel_observation_arms(
        replay(
            {
                "baseline": (endpoint("a"), endpoint("b", call=False)),
                "relaxed": (endpoint("a", "relaxed"), endpoint("b", "relaxed", call=None)),
            }
        ),
        samples,
    )
    document = observation_arms_document(arms)
    assert document["schema_version"] == "calibration-cutoff-observations-v1"
    assert document["policy_ids"] == ["baseline", "relaxed"]
    assert document["keys"] == ["a", "b"]
    restored = json.loads(json.dumps(document))
    assert restored == document
    decoded = {
        arm["policy_id"]: tuple(
            CallerObservation(
                row["key"],
                row["group_key"],
                row["truth_positive"],
                None if row["truth_variants"] is None else tuple(row["truth_variants"]),
                row["called_positive"],
                tuple(row["called_variants"]),
                tuple(row["tier_a_variants"]),
            )
            for row in arm["observations"]
        )
        for arm in restored["arms"]
    }
    assert decoded == arms


def test_observation_arms_document_refuses_disagreeing_arms():
    from vntyper.scripts.calibration_cutoff_observations import observation_arms_document

    with pytest.raises(ValueError, match="identical"):
        observation_arms_document(
            {
                "baseline": (CallerObservation("a", "sample:a", True, None, True, (RIGHT,), ()),),
                "relaxed": (CallerObservation("a", "sample:a", False, (), False, (), ()),),
            }
        )
