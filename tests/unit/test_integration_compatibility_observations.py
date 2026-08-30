"""Tests for versioned append-only integration observations."""

import copy
import importlib
from typing import Any

import pytest

from tests.unit.test_integration_compatibility import _contract, _manifest

pytestmark = pytest.mark.unit

VERSION = "2.0.24"
PROVENANCE = "f9e57f73e10d88d0c27cc4c4e8501c892594f0db"
CURRENT_REPORT = ["current report"]


def _module() -> Any:
    """Import the policy module while keeping the initial RED collectable."""
    try:
        return importlib.import_module("scripts.integration_compatibility_observations")
    except ModuleNotFoundError:
        return None


def _override(*, suite: str = "bam_tests", test_name: str = "case-a") -> dict[str, Any]:
    return {"suite": suite, "test_name": test_name, "report": [*CURRENT_REPORT]}


def _set(
    version: str = VERSION,
    *,
    extends: str | None = None,
    overrides: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    return {
        "version": version,
        "provenance_commit": PROVENANCE,
        "extends": extends,
        "report_overrides": [_override()] if overrides is None else overrides,
    }


def _v2(*contracts: dict, sets: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    return {
        "schema_version": 2,
        "contracts": list(contracts),
        "observation_sets": [_set()] if sets is None else sets,
    }


def test_selects_exact_current_version_without_mutating_historical_contracts() -> None:
    """A selected observation changes only the effective report projection."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"
    contract = _contract()
    original = copy.deepcopy(contract)
    indexed = {("bam_tests", "case-a"): contract}

    effective = module.effective_contracts(_v2(contract), indexed, VERSION)

    assert contract == original
    assert effective[("bam_tests", "case-a")]["outcomes"]["report"] == CURRENT_REPORT
    assert effective[("bam_tests", "case-a")]["outcomes"]["kestrel"] == original["outcomes"]["kestrel"]


def test_legacy_manifest_returns_an_independent_contract_copy() -> None:
    """Schema v1 remains supported without observations or shared mutable state."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"
    contract = _contract()
    indexed = {("bam_tests", "case-a"): contract}

    effective = module.effective_contracts(_manifest(contract), indexed, None)
    effective[("bam_tests", "case-a")]["outcomes"]["report"].append("changed")

    assert indexed[("bam_tests", "case-a")]["outcomes"]["report"] == ["High_Precision*"]


def test_public_observation_validator_dispatches_by_schema() -> None:
    """Schema v1 is a no-op while schema v2 validates nested observations."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"
    contract = _contract()

    assert module.validate_observation_sets(_manifest(contract), {("bam_tests", "case-a")}) is None
    malformed = _v2(contract, sets=[])
    with pytest.raises(ValueError) as error:
        module.validate_observation_sets(malformed, {("bam_tests", "case-a")})

    assert str(error.value) == "manifest observation_sets must be a non-empty list"


@pytest.mark.parametrize("version", [None, "2.0.23", "2.0.25", "2.0", "v2.0.24", "2.0.24-rc1"])
def test_current_version_must_resolve_exactly_to_the_final_observation_set(version: str | None) -> None:
    """Missing, unknown, prerelease, and stale selectors fail closed."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"
    contract = _contract()

    with pytest.raises(ValueError, match=r"version|final observation"):
        module.effective_contracts(_v2(contract), {("bam_tests", "case-a"): contract}, version)


def test_versions_are_ordered_numerically_and_extend_the_immediate_predecessor() -> None:
    """Semantic numeric ordering must not degrade to lexical string ordering."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"
    contract = _contract()
    manifest = _v2(
        contract,
        sets=[
            _set("2.0.9"),
            _set("2.0.10", extends="2.0.9", overrides=[]),
        ],
    )

    effective = module.effective_contracts(manifest, {("bam_tests", "case-a"): contract}, "2.0.10")

    assert effective[("bam_tests", "case-a")]["outcomes"]["report"] == CURRENT_REPORT
    malformed = copy.deepcopy(manifest)
    malformed["observation_sets"][1]["extends"] = "2.0.8"
    with pytest.raises(ValueError, match="immediately preceding"):
        module.effective_contracts(malformed, {("bam_tests", "case-a"): contract}, "2.0.10")


@pytest.mark.parametrize(
    ("sets", "message"),
    [
        ([_set(), _set()], "duplicate observation version"),
        ([_set("2.0.10"), _set("2.0.9", extends="2.0.10")], "strictly increasing"),
        ([_set("2.0.9", extends="2.0.8")], "first observation set"),
        ([_set("02.0.9")], "plain release version"),
    ],
)
def test_observation_set_chain_rejects_ambiguous_or_malformed_versions(
    sets: list[dict[str, Any]], message: str
) -> None:
    """The version chain is unique, strict, ordered, and rooted."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"
    contract = _contract()

    with pytest.raises(ValueError, match=message):
        module.effective_contracts(_v2(contract, sets=sets), {("bam_tests", "case-a"): contract}, VERSION)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda override: override.update(extra=True), "extra keys"),
        (lambda override: override.update(report=[]), "must not be empty"),
        (lambda override: override.update(report=["same", "same"]), "duplicate"),
        (lambda override: override.update(report=[1]), "non-empty string"),
        (lambda override: override.update(test_name="missing"), "unknown contract identity"),
    ],
)
def test_report_override_schema_is_strict(mutation: Any, message: str) -> None:
    """An observation overlay cannot modify protected non-report fields."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"
    contract = _contract()
    override = _override()
    mutation(override)

    with pytest.raises(ValueError, match=message):
        module.effective_contracts(
            _v2(contract, sets=[_set(overrides=[override])]),
            {("bam_tests", "case-a"): contract},
            VERSION,
        )


def test_report_overrides_reject_duplicate_identity() -> None:
    """Two competing report observations in one version are ambiguous."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"
    contract = _contract()

    with pytest.raises(ValueError, match="duplicate report override"):
        module.effective_contracts(
            _v2(contract, sets=[_set(overrides=[_override(), _override()])]),
            {("bam_tests", "case-a"): contract},
            VERSION,
        )


@pytest.mark.parametrize(
    ("manifest", "message"),
    [
        ([], "manifest must be an object"),
        ({"schema_version": 1}, "manifest has missing keys: ['contracts']"),
        ({"schema_version": 1, "contracts": {}}, "manifest contracts must be a list"),
        (
            {"schema_version": 2, "contracts": [], "observation_sets": {}},
            "manifest observation_sets must be a list",
        ),
    ],
)
def test_manifest_container_rejects_each_malformed_top_level_shape(manifest: object, message: str) -> None:
    """Top-level object, keys, and list containers fail with exact diagnostics."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"

    with pytest.raises(ValueError) as error:
        module.validate_manifest_container(manifest)

    assert str(error.value) == message


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda manifest: manifest.update(observation_sets=[]), "manifest observation_sets must be a non-empty list"),
        (
            lambda manifest: manifest["observation_sets"].__setitem__(0, "not-an-object"),
            "observation_sets[0] must be an object",
        ),
        (
            lambda manifest: manifest["observation_sets"][0].update(provenance_commit="ABC"),
            "observation_sets[0].provenance_commit must be a full lowercase Git SHA",
        ),
        (
            lambda manifest: manifest["observation_sets"][0].update(report_overrides={}),
            "observation_sets[0].report_overrides must be a list",
        ),
        (
            lambda manifest: manifest["observation_sets"][0]["report_overrides"][0].update(report="not-a-list"),
            "observation_sets[0].report_overrides[0].report must be a list",
        ),
    ],
)
def test_observation_sets_reject_each_malformed_nested_container(mutation: Any, message: str) -> None:
    """Observation containers, provenance, overrides, and reports fail closed."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"
    contract = _contract()
    manifest = _v2(contract)
    mutation(manifest)

    with pytest.raises(ValueError) as error:
        module.effective_contracts(manifest, {("bam_tests", "case-a"): contract}, VERSION)

    assert str(error.value) == message


@pytest.mark.parametrize(
    ("base_contracts", "current_contracts"),
    [({}, []), ([], {}), (None, []), ([], None)],
)
def test_append_only_history_rejects_non_list_contract_containers(
    base_contracts: object, current_contracts: object
) -> None:
    """History comparison requires lists on both sides before prefix checks."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"
    base = {"schema_version": 1, "contracts": base_contracts}
    current = {"schema_version": 1, "contracts": current_contracts}

    with pytest.raises(ValueError) as error:
        module.validate_append_only_history(base, current)

    assert str(error.value) == "compatibility contracts must be lists"


def test_schema_v1_to_v2_migration_preserves_contracts_as_an_ordered_prefix() -> None:
    """The only schema migration adds observations without rewriting history."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"
    first = _contract()
    second = _contract(name="case-b", path="tests/data/second.bam")
    base = _manifest(first, second)

    module.validate_append_only_history(base, _v2(first, second))

    for mutation in (
        lambda rows: rows.pop(),
        lambda rows: rows.reverse(),
        lambda rows: rows[0]["outcomes"]["report"].append("mutated"),
        lambda rows: rows.insert(0, _contract(name="inserted")),
    ):
        current = _v2(copy.deepcopy(first), copy.deepcopy(second))
        mutation(current["contracts"])
        with pytest.raises(ValueError, match=r"removed|mutated|ordered prefix"):
            module.validate_append_only_history(base, current)


@pytest.mark.parametrize(
    ("path", "replacement"),
    [
        (("outcomes", "kestrel", "Estimated_Depth_AlternateVariant", "value"), 5.0),
        (("artifacts", "archive"), 0),
    ],
)
def test_append_only_history_is_type_exact(path: tuple[str, ...], replacement: object) -> None:
    """Python-equivalent integer, float, and Boolean values are distinct JSON evidence."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"
    contract = _contract()
    current_contract = copy.deepcopy(contract)
    target = current_contract
    for part in path[:-1]:
        target = target[part]
    target[path[-1]] = replacement

    with pytest.raises(ValueError, match="compatibility contracts were mutated"):
        module.validate_append_only_history(_manifest(contract), _v2(current_contract))


def test_schema_v2_history_permits_only_appended_contracts_and_observation_sets() -> None:
    """Existing observation sets retain exact content and order."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"
    contract = _contract()
    base = _v2(contract)
    appended = copy.deepcopy(base)
    appended["observation_sets"].append(_set("2.0.25", extends=VERSION, overrides=[]))

    module.validate_append_only_history(base, appended)

    for mutation in (
        lambda sets: sets.clear(),
        lambda sets: sets[0].update(provenance_commit="a" * 40),
        lambda sets: sets.insert(0, _set("2.0.23")),
    ):
        current = copy.deepcopy(appended)
        mutation(current["observation_sets"])
        with pytest.raises(ValueError, match="observation sets were removed or mutated"):
            module.validate_append_only_history(base, current)


@pytest.mark.parametrize(
    ("base_schema", "current_schema"),
    [(2, 1), (1, 3), (2, 3), (3, 3)],
)
def test_only_schema_one_to_two_migration_is_supported(base_schema: int, current_schema: int) -> None:
    """Downgrades and unrecognized schema transitions fail closed."""
    module = _module()
    assert module is not None, "versioned observation policy is not implemented"
    contract = _contract()
    base = _manifest(contract) if base_schema == 1 else _v2(contract)
    current = _manifest(contract) if current_schema == 1 else _v2(contract)
    base["schema_version"] = base_schema
    current["schema_version"] = current_schema

    with pytest.raises(ValueError, match="schema transition"):
        module.validate_append_only_history(base, current)
