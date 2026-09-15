"""Identity audit preserves biological linkage, technical arms and split safety."""

from copy import deepcopy
from dataclasses import replace
from hashlib import sha256
from importlib import import_module

import pytest

from tests.unit.test_calibration_intake_contract import synthetic_intake
from vntyper.scripts.calibration_intake_contract import decode_intake
from vntyper.scripts.calibration_read_fingerprints import LogicalReadFingerprint

pytestmark = pytest.mark.unit


def digest(text):
    return sha256(text.encode()).hexdigest()


def intake_for(count=2):
    template = synthetic_intake()
    raw = {
        "schema_version": template["schema_version"],
        **{key: [] for key in ("specimens", "artifacts", "aliases", "truth", "assignments")},
    }
    for i in range(count):
        key = f"sample-{i}"
        specimen = deepcopy(template["specimens"][0])
        specimen.update(key=key, individual_key=f"individual-{i}")
        artifact = deepcopy(template["artifacts"][0])
        artifact.update(
            key=f"artifact-{i}",
            specimen_key=key,
            path=f"/synthetic/input-{i}.bam",
            expected_sha256=None,
            replicate_group=f"replicate-{i}",
        )
        assignment = deepcopy(template["assignments"][0])
        assignment.update(specimen_key=key, groups={name: [f"{name}:{i}"] for name in assignment["groups"]})
        raw["specimens"].append(specimen)
        raw["artifacts"].append(artifact)
        raw["assignments"].append(assignment)
    return raw


def fingerprints_for(raw):
    m = import_module("vntyper.scripts.calibration_identity")
    return {
        artifact["key"]: m.ArtifactFingerprint(
            artifact["key"],
            digest(artifact["key"] + ":bytes"),
            None,
            LogicalReadFingerprint(
                digest(artifact["key"] + ":alignment"),
                digest(artifact["key"] + ":named"),
                digest(artifact["key"] + ":unnamed"),
                2,
                True,
                (),
            ),
        )
        for artifact in raw["artifacts"]
    }


def audit(raw, fingerprints=None, priority=("synthetic-preprocessing-v1",)):
    m = import_module("vntyper.scripts.calibration_identity")
    return m.resolve_identities(
        decode_intake(raw),
        fingerprints_for(raw) if fingerprints is None else fingerprints,
        preprocessing_priority=priority,
    )


def test_distinct_read_sets_do_not_override_explicit_biological_linkage():
    raw = intake_for()
    raw["specimens"][1]["individual_key"] = raw["specimens"][0]["individual_key"]
    result = audit(raw)
    assert result.specimen_groups["sample-0"] == result.specimen_groups["sample-1"]
    assert len(result.primary_artifact_by_group) == 1
    assert result.primary_artifact_by_specimen == {"sample-0": "artifact-0", "sample-1": "artifact-1"}
    assert result.quarantined_specimens == {}


def test_byte_mirrors_share_execution_but_reprocessing_arms_remain_declared():
    raw = intake_for(1)
    mirror = deepcopy(raw["artifacts"][0])
    mirror.update(key="artifact-1", path="/synthetic/mirror.bam")
    arm = deepcopy(mirror)
    arm.update(key="artifact-2", preprocessing_id="alternative")
    raw["artifacts"].extend([mirror, arm])
    fps = fingerprints_for(raw)
    for key in ("artifact-1", "artifact-2"):
        fps[key] = replace(fps["artifact-0"], artifact_key=key)
    result = audit(raw, fps, priority=("synthetic-preprocessing-v1", "alternative"))
    assert result.execution_representatives == {
        "artifact-0": "artifact-0",
        "artifact-1": "artifact-0",
        "artifact-2": "artifact-2",
    }
    assert result.primary_artifact_by_specimen == {"sample-0": "artifact-0"}


def test_logical_duplicates_do_not_merge_people_and_are_quarantined_without_linkage():
    raw = intake_for()
    fps = fingerprints_for(raw)
    fps["artifact-1"] = replace(fps["artifact-1"], logical=fps["artifact-0"].logical)
    result = audit(raw, fps)
    assert result.specimen_groups["sample-0"] == result.specimen_groups["sample-1"]
    assert result.primary_artifact_by_specimen.keys() == {"sample-0", "sample-1"}
    assert "duplicate_input_identity_requires_adjudication" in result.quarantined_specimens["sample-0"]
    assert result.execution_representatives["artifact-1"] == "artifact-1"


def test_sequence_only_collision_is_a_review_reason_not_a_silent_merge():
    raw = intake_for()
    fps = fingerprints_for(raw)
    fps["artifact-1"] = replace(
        fps["artifact-1"],
        logical=replace(
            fps["artifact-1"].logical, unnamed_sequence_sha256=fps["artifact-0"].logical.unnamed_sequence_sha256
        ),
    )
    result = audit(raw, fps)
    assert "sequence_only_collision_requires_adjudication" in result.quarantined_specimens["sample-1"]
    assert result.execution_representatives["artifact-1"] == "artifact-1"


@pytest.mark.parametrize("linkage", ["individual", "family", "replicate", "logical"])
def test_direct_and_transitive_evidence_linkage_cannot_cross_roles(linkage):
    raw = intake_for()
    raw["assignments"][1]["role"] = "policy-selection"
    fps = fingerprints_for(raw)
    if linkage == "individual":
        raw["specimens"][1]["individual_key"] = raw["specimens"][0]["individual_key"]
    elif linkage == "family":
        for specimen in raw["specimens"]:
            specimen["family_key"] = "shared-family"
    elif linkage == "replicate":
        raw["artifacts"][1]["replicate_group"] = raw["artifacts"][0]["replicate_group"]
    else:
        fps["artifact-1"] = replace(fps["artifact-1"], logical=fps["artifact-0"].logical)
    with pytest.raises(ValueError, match="cross roles"):
        audit(raw, fps)


def test_primary_selection_uses_predeclared_priority_then_artifact_key():
    raw = intake_for(1)
    second = deepcopy(raw["artifacts"][0])
    second.update(key="artifact-z", preprocessing_id="preferred")
    raw["artifacts"].append(second)
    result = audit(raw, priority=("preferred", "synthetic-preprocessing-v1"))
    assert result.primary_artifact_by_specimen["sample-0"] == "artifact-z"
    raw["artifacts"].reverse()
    assert audit(raw, priority=("preferred", "synthetic-preprocessing-v1")) == result


def test_conflicting_truth_is_preserved_as_quarantine_not_selected_by_method_order():
    raw = intake_for(1)
    truth = synthetic_intake()["truth"][0]
    truth.update(specimen_key="sample-0", source_row=1)
    conflict = deepcopy(truth)
    conflict.update(genotype="negative", variants=[], source_row=2)
    raw["truth"] = [truth, conflict]
    result = audit(raw)
    assert "genotype_truth_conflict" in result.quarantined_specimens["sample-0"]


def test_unreliable_sequence_identity_is_development_only():
    raw = intake_for(1)
    fps = fingerprints_for(raw)
    fps["artifact-0"] = replace(
        fps["artifact-0"],
        logical=replace(
            fps["artifact-0"].logical, sequence_identity_reliable=False, reasons=("hard_clipped_sequence",)
        ),
    )
    result = audit(raw, fps)
    assert "unreliable_sequence_identity" in result.quarantined_specimens["sample-0"]
    raw["assignments"][0]["role"] = "validation"
    with pytest.raises(ValueError, match="confirmatory"):
        audit(raw, fps)


def test_fingerprints_and_expected_bytes_must_match_the_exact_artifact_roster():
    raw = intake_for(1)
    fps = fingerprints_for(raw)
    with pytest.raises(ValueError, match="artifact"):
        audit(raw, {})
    raw["artifacts"][0]["expected_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="expected"):
        audit(raw, fps)


def test_partition_members_keep_artifact_keys_and_audit_is_revalidated():
    m = import_module("vntyper.scripts.calibration_identity")
    raw = intake_for()
    declaration = decode_intake(raw)
    result = audit(raw)
    members = m.build_partition_members(declaration, result)
    assert tuple(member.key for member in members) == ("artifact-0", "artifact-1")
    assert all(member.role == "training" for member in members)
    with pytest.raises(ValueError, match="audit"):
        m.build_partition_members(declaration, replace(result, sha256="0" * 64))


def test_conflicting_truth_across_two_labels_for_one_individual_is_quarantined():
    raw = intake_for()
    raw["specimens"][1]["individual_key"] = raw["specimens"][0]["individual_key"]
    for i, genotype in enumerate(("positive", "negative")):
        truth = synthetic_intake()["truth"][0]
        truth.update(specimen_key=f"sample-{i}", genotype=genotype, variants=[])
        raw["truth"].append(truth)
    result = audit(raw)
    assert all("genotype_truth_conflict" in result.quarantined_specimens[key] for key in ("sample-0", "sample-1"))


def test_different_truth_for_family_members_is_not_a_specimen_truth_conflict():
    raw = intake_for()
    for i, genotype in enumerate(("positive", "negative")):
        raw["specimens"][i]["family_key"] = "family"
        truth = synthetic_intake()["truth"][0]
        truth.update(specimen_key=f"sample-{i}", genotype=genotype, variants=[])
        raw["truth"].append(truth)
    assert audit(raw).quarantined_specimens == {}


@pytest.mark.parametrize(
    "field,value",
    [("byte_sha256", "invalid"), ("mate_byte_sha256", "a" * 64), ("artifact_key", "other"), ("logical", None)],
)
def test_malformed_physical_fingerprints_fail_before_grouping(field, value):
    raw = intake_for(1)
    fps = fingerprints_for(raw)
    fps["artifact-0"] = replace(fps["artifact-0"], **{field: value})
    with pytest.raises(ValueError):
        audit(raw, fps)


@pytest.mark.parametrize(
    "field,value",
    [
        ("primary_record_count", True),
        ("primary_record_count", 0),
        ("alignment_sha256", "invalid"),
        ("sequence_identity_reliable", False),
        ("reasons", ["missing_sequence"]),
        ("reasons", ("unknown_reason",)),
    ],
)
def test_malformed_logical_fingerprints_cannot_enter_audit(field, value):
    raw = intake_for(1)
    fps = fingerprints_for(raw)
    fps["artifact-0"] = replace(fps["artifact-0"], logical=replace(fps["artifact-0"].logical, **{field: value}))
    with pytest.raises(ValueError):
        audit(raw, fps)


def test_fastq_pairs_require_and_bind_both_file_hashes():
    raw = intake_for(1)
    raw["artifacts"][0].update(format="FASTQ_PAIR", mate_path="/synthetic/mate.fastq")
    fps = fingerprints_for(raw)
    with pytest.raises(ValueError, match="digest"):
        audit(raw, fps)
    fps["artifact-0"] = replace(fps["artifact-0"], mate_byte_sha256="a" * 64)
    result = audit(raw, fps)
    assert result.fingerprints["artifact-0"].mate_byte_sha256 == "a" * 64


def test_byte_identical_inputs_cannot_claim_different_logical_evidence():
    raw = intake_for()
    fps = fingerprints_for(raw)
    fps["artifact-1"] = replace(fps["artifact-1"], byte_sha256=fps["artifact-0"].byte_sha256)
    with pytest.raises(ValueError, match="byte-identical"):
        audit(raw, fps)


def test_equal_logical_hashes_cannot_have_different_occurrence_counts():
    raw = intake_for()
    fps = fingerprints_for(raw)
    fps["artifact-1"] = replace(fps["artifact-1"], logical=replace(fps["artifact-0"].logical, primary_record_count=3))
    with pytest.raises(ValueError, match="occurrence"):
        audit(raw, fps)


@pytest.mark.parametrize(
    "priority", [(), ("unknown",), ("synthetic-preprocessing-v1", "synthetic-preprocessing-v1"), (True,)]
)
def test_preprocessing_priority_is_complete_unique_and_frozen(priority):
    with pytest.raises(ValueError, match="priority"):
        audit(intake_for(1), priority=priority)


def test_locked_metadata_does_not_create_a_fictitious_input_member():
    m = import_module("vntyper.scripts.calibration_identity")
    raw = intake_for()
    raw["artifacts"].pop()
    raw["assignments"][1].update(role="locked-heldout", provenance="external-custodian")
    result = audit(raw)
    assert "sample-1" in result.specimen_groups
    assert "sample-1" not in result.primary_artifact_by_specimen
    assert tuple(member.key for member in m.build_partition_members(decode_intake(raw), result)) == ("artifact-0",)


@pytest.mark.parametrize(
    "kind,reason",
    [("disputed", "disputed_truth"), ("variant", "variant_truth_conflict"), ("length", "length_truth_conflict")],
)
def test_all_independent_truth_conflicts_require_adjudication(kind, reason):
    raw = intake_for(1)
    first = synthetic_intake()["truth"][0]
    first.update(specimen_key="sample-0")
    second = deepcopy(first)
    second["source_row"] = 2
    if kind == "disputed":
        second["status"] = "disputed"
    elif kind == "variant":
        second["variants"] = ["different-synthetic-variant"]
    else:
        for row, allele in ((first, 40), (second, 50)):
            row["length"] = {
                "allele_1": allele,
                "allele_2": 90,
                "unit": "repeat-count",
                "repeat_unit_bp": 60,
                "boundary_definition": "target-v1",
                "measurement": "exact",
                "lower_bound": None,
                "upper_bound": None,
                "conversion_id": None,
            }
    raw["truth"] = [first, second]
    assert reason in audit(raw).quarantined_specimens["sample-0"]


def test_unresolved_development_identity_is_visible_and_does_not_become_confirmed():
    raw = intake_for(1)
    raw["specimens"][0].update(individual_key=None, identity_status="unresolved")
    assert audit(raw).quarantined_specimens["sample-0"] == ("unresolved_biological_identity",)


def test_stale_intake_and_untyped_audit_are_refused():
    m = import_module("vntyper.scripts.calibration_identity")
    raw = intake_for(1)
    declaration = decode_intake(raw)
    with pytest.raises(ValueError, match="intake digest"):
        m.resolve_identities(
            replace(declaration, sha256="0" * 64),
            fingerprints_for(raw),
            preprocessing_priority=("synthetic-preprocessing-v1",),
        )
    with pytest.raises(ValueError, match="IdentityAudit"):
        m.build_partition_members(declaration, {})


def test_partition_projection_returns_canonical_order_from_replaced_intake():
    m = import_module("vntyper.scripts.calibration_identity")
    raw = intake_for()
    declaration = decode_intake(raw)
    result = audit(raw)
    reordered = replace(declaration, artifacts=tuple(reversed(declaration.artifacts)))
    assert tuple(member.key for member in m.build_partition_members(reordered, result)) == ("artifact-0", "artifact-1")
