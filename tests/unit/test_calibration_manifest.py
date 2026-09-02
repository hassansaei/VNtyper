"""Role-bound manifests and transitive leakage groups."""

import pytest

from tests.unit.test_calibration_contract import synthetic_protocol
from vntyper.scripts.calibration_manifest import (
    connected_leakage_groups,
    decode_partition_manifest,
    decode_study_declaration,
    require_operation_roles,
)

pytestmark = pytest.mark.unit

_GROUP_KEYS = (
    "individual-family",
    "simulated-pair",
    "backbone-seed-lineage",
    "replicate-rerun",
    "depth-series-source",
    "batch",
    "repeat-context",
)


def _member(
    key: str,
    role: str,
    *,
    groups: dict[str, list[str]] | None = None,
    assay_class: str = "capture-short-read",
) -> dict[str, object]:
    complete = {namespace: [f"{namespace}:{key}"] for namespace in _GROUP_KEYS}
    if groups:
        complete.update(groups)
    return {
        "key": key,
        "role": role,
        "provenance": "development",
        "assay_class": assay_class,
        "groups": complete,
    }


def _manifest(*members: dict[str, object]) -> dict[str, object]:
    return {"schema_version": "calibration-partitions-v1", "members": list(members)}


def test_connected_groups_are_transitive_not_only_adjacent() -> None:
    manifest = decode_partition_manifest(
        _manifest(
            _member("a", "training", groups={"individual-family": ["person-1"]}),
            _member(
                "b",
                "training",
                groups={"individual-family": ["person-1"], "replicate-rerun": ["rerun-9"]},
            ),
            _member("c", "training", groups={"replicate-rerun": ["rerun-9"]}),
            _member("d", "policy-selection"),
        )
    )

    groups = connected_leakage_groups(manifest)

    assert groups["a"] == groups["b"] == groups["c"]
    assert groups["d"] != groups["a"]
    assert len(set(groups.values())) == 2


@pytest.mark.parametrize("namespace", _GROUP_KEYS)
def test_every_group_namespace_blocks_cross_role_leakage(namespace: str) -> None:
    raw = _manifest(
        _member("a-train", "training", groups={namespace: ["shared"]}),
        _member("b-select", "policy-selection", groups={namespace: ["shared"]}),
    )

    with pytest.raises(ValueError, match="leakage|roles"):
        decode_partition_manifest(raw)


def test_transitive_chain_crossing_a_role_is_rejected() -> None:
    raw = _manifest(
        _member("a", "training", groups={"individual-family": ["person-1"]}),
        _member(
            "b",
            "training",
            groups={"individual-family": ["person-1"], "replicate-rerun": ["rerun-9"]},
        ),
        _member("c", "validation", groups={"replicate-rerun": ["rerun-9"]}),
    )

    with pytest.raises(ValueError, match="leakage|roles"):
        decode_partition_manifest(raw)


def test_group_namespace_cannot_be_dropped_or_renamed() -> None:
    missing = _member("a", "training")
    missing_groups = missing["groups"]
    assert isinstance(missing_groups, dict)
    del missing_groups["batch"]
    renamed = _member("b", "training")
    renamed_groups = renamed["groups"]
    assert isinstance(renamed_groups, dict)
    renamed_groups["lab"] = renamed_groups.pop("batch")

    for member in (missing, renamed):
        with pytest.raises(ValueError, match="group"):
            decode_partition_manifest(_manifest(member))


def test_fit_roles_are_distinct_and_validation_or_heldout_access_is_refused() -> None:
    manifest = decode_partition_manifest(
        _manifest(
            {**_member("held", "locked-heldout"), "provenance": "external-custodian"},
            _member("select", "policy-selection"),
            _member("train", "training"),
            _member("validate", "validation"),
        )
    )

    assert require_operation_roles(manifest, "fit") == ("training", "policy-selection")
    assert require_operation_roles(manifest, "validate") == ("validation",)
    assert require_operation_roles(manifest, "evaluate") == ("locked-heldout",)
    with pytest.raises(ValueError):
        require_operation_roles(manifest, "fit", requested_roles=("training", "validation"))
    with pytest.raises(ValueError):
        require_operation_roles(manifest, "fit", requested_roles=("policy-selection",))


def test_member_keys_roles_provenance_and_order_are_strict() -> None:
    duplicate = _manifest(_member("a", "training"), _member("a", "training"))
    reversed_members = _manifest(_member("b", "training"), _member("a", "training"))
    bad_role = _manifest(_member("a", "test"))
    relabeled = _member("a", "locked-heldout")

    for raw in (duplicate, reversed_members, bad_role, _manifest(relabeled)):
        with pytest.raises(ValueError):
            decode_partition_manifest(raw)


def test_study_declaration_hashes_protocol_and_all_four_roles() -> None:
    partitions = _manifest(
        {**_member("held", "locked-heldout"), "provenance": "external-custodian"},
        _member("select", "policy-selection"),
        _member("train", "training"),
        _member("validate", "validation"),
    )
    raw = {"schema_version": "calibration-study-v1", "protocol": synthetic_protocol(), "partitions": partitions}

    study = decode_study_declaration(raw)

    assert study.protocol.objective == "lexicographic-safety-v1"
    assert {member.role for member in study.partitions.members} == {
        "training",
        "policy-selection",
        "validation",
        "locked-heldout",
    }
    assert len(study.sha256) == 64

    incomplete = dict(raw)
    incomplete["partitions"] = _manifest(_member("train", "training"))
    with pytest.raises(ValueError, match="four"):
        decode_study_declaration(incomplete)


def test_study_declaration_binds_each_member_to_a_predeclared_assay_class() -> None:
    partitions = _manifest(
        {**_member("held", "locked-heldout"), "provenance": "external-custodian"},
        _member("select", "policy-selection"),
        _member("train", "training", assay_class="genome-short-read"),
        _member("validate", "validation"),
    )
    raw = {"schema_version": "calibration-study-v1", "protocol": synthetic_protocol(), "partitions": partitions}

    with pytest.raises(ValueError, match="assay class.*protocol|protocol.*assay class"):
        decode_study_declaration(raw)
