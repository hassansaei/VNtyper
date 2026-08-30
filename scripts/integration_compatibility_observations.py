"""Pure policy for append-only, versioned integration observations."""

from __future__ import annotations

import copy
import json
import re
from typing import Any

Manifest = dict[str, Any]
Identity = tuple[str, str]
ContractIndex = dict[Identity, dict[str, Any]]

LEGACY_SCHEMA_VERSION = 1
OBSERVATION_SCHEMA_VERSION = 2
SUPPORTED_SCHEMA_VERSIONS = {LEGACY_SCHEMA_VERSION, OBSERVATION_SCHEMA_VERSION}

_LEGACY_TOP_KEYS = {"schema_version", "contracts"}
_OBSERVATION_TOP_KEYS = {"schema_version", "contracts", "observation_sets"}
_OBSERVATION_SET_KEYS = {"version", "provenance_commit", "extends", "report_overrides"}
_REPORT_OVERRIDE_KEYS = {"suite", "test_name", "report"}
_PLAIN_RELEASE_RE = re.compile(r"(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)")
_SHA_RE = re.compile(r"[0-9a-f]{40}")


def _mapping(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    return value


def _exact_keys(value: dict[str, Any], expected: set[str], label: str) -> None:
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing:
        raise ValueError(f"{label} has missing keys: {missing}")
    if extra:
        raise ValueError(f"{label} has extra keys: {extra}")


def _nonempty_string(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} must be a non-empty string")
    return value


def _release_version(value: object, label: str) -> tuple[int, int, int]:
    if not isinstance(value, str) or _PLAIN_RELEASE_RE.fullmatch(value) is None:
        raise ValueError(f"{label} must be a plain release version MAJOR.MINOR.PATCH")
    return tuple(int(part) for part in value.split("."))  # type: ignore[return-value]


def _schema_version(manifest: Manifest, label: str) -> int:
    root = _mapping(manifest, label)
    version = root.get("schema_version")
    if isinstance(version, bool) or not isinstance(version, int) or version not in SUPPORTED_SCHEMA_VERSIONS:
        raise ValueError(f"{label} schema_version must be one of {sorted(SUPPORTED_SCHEMA_VERSIONS)}")
    return version


def validate_manifest_container(manifest: Manifest) -> int:
    """Validate the schema-specific top-level manifest container.

    Args:
        manifest: Parsed compatibility manifest.

    Returns:
        Validated schema version.

    Raises:
        ValueError: If the schema or top-level shape is unsupported.
    """
    version = _schema_version(manifest, "manifest")
    expected = _LEGACY_TOP_KEYS if version == LEGACY_SCHEMA_VERSION else _OBSERVATION_TOP_KEYS
    _exact_keys(manifest, expected, "manifest")
    if not isinstance(manifest["contracts"], list):
        raise ValueError("manifest contracts must be a list")
    if version == OBSERVATION_SCHEMA_VERSION and not isinstance(manifest["observation_sets"], list):
        raise ValueError("manifest observation_sets must be a list")
    return version


def _report(value: object, label: str) -> list[str]:
    if not isinstance(value, list):
        raise ValueError(f"{label} must be a list")
    result = [_nonempty_string(item, f"{label}[{index}]") for index, item in enumerate(value)]
    if not result:
        raise ValueError(f"{label} must not be empty")
    if len(set(result)) != len(result):
        raise ValueError(f"{label} contains a duplicate")
    return result


def _observation_sets(manifest: Manifest, identities: set[Identity]) -> list[dict[str, Any]]:
    values = manifest["observation_sets"]
    if not isinstance(values, list) or not values:
        raise ValueError("manifest observation_sets must be a non-empty list")
    result: list[dict[str, Any]] = []
    seen_versions: set[str] = set()
    previous_version: str | None = None
    previous_parsed: tuple[int, int, int] | None = None
    for set_index, raw_set in enumerate(values):
        observation = _mapping(raw_set, f"observation_sets[{set_index}]")
        _exact_keys(observation, _OBSERVATION_SET_KEYS, f"observation_sets[{set_index}]")
        version = observation["version"]
        parsed = _release_version(version, f"observation_sets[{set_index}].version")
        if version in seen_versions:
            raise ValueError(f"duplicate observation version: {version}")
        if previous_parsed is not None and parsed <= previous_parsed:
            raise ValueError("observation versions must be strictly increasing")
        extends = observation["extends"]
        if set_index == 0:
            if extends is not None:
                raise ValueError("first observation set must extend null")
        elif extends != previous_version:
            raise ValueError("observation set must extend the immediately preceding version")
        provenance = observation["provenance_commit"]
        if not isinstance(provenance, str) or _SHA_RE.fullmatch(provenance) is None:
            raise ValueError(f"observation_sets[{set_index}].provenance_commit must be a full lowercase Git SHA")
        overrides = observation["report_overrides"]
        if not isinstance(overrides, list):
            raise ValueError(f"observation_sets[{set_index}].report_overrides must be a list")
        seen_overrides: set[Identity] = set()
        for override_index, raw_override in enumerate(overrides):
            label = f"observation_sets[{set_index}].report_overrides[{override_index}]"
            override = _mapping(raw_override, label)
            _exact_keys(override, _REPORT_OVERRIDE_KEYS, label)
            identity = (
                _nonempty_string(override["suite"], f"{label}.suite"),
                _nonempty_string(override["test_name"], f"{label}.test_name"),
            )
            if identity not in identities:
                raise ValueError(f"{label} refers to unknown contract identity {identity}")
            if identity in seen_overrides:
                raise ValueError(f"duplicate report override for {identity}")
            _report(override["report"], f"{label}.report")
            seen_overrides.add(identity)
        result.append(observation)
        seen_versions.add(version)
        previous_version = version
        previous_parsed = parsed
    return result


def validate_observation_sets(manifest: Manifest, identities: set[Identity]) -> None:
    """Validate every schema-v2 observation against historical identities.

    Args:
        manifest: Parsed compatibility manifest.
        identities: Validated contract identities available for report overrides.

    Raises:
        ValueError: If an observation set is malformed or ambiguous.
    """
    schema = validate_manifest_container(manifest)
    if schema == OBSERVATION_SCHEMA_VERSION:
        _observation_sets(manifest, identities)


def effective_contracts(
    manifest: Manifest,
    indexed_contracts: ContractIndex,
    current_version: str | None,
) -> ContractIndex:
    """Resolve immutable contracts to the exact checkout-version observation.

    Args:
        manifest: Validated current compatibility manifest.
        indexed_contracts: Historical contracts indexed by identity.
        current_version: Authoritative package version from the current checkout.

    Returns:
        A deep-copied contract index with report observations applied.

    Raises:
        ValueError: If selection or any observation is malformed or ambiguous.
    """
    schema = validate_manifest_container(manifest)
    if schema == LEGACY_SCHEMA_VERSION:
        return copy.deepcopy(indexed_contracts)
    parsed_current = _release_version(current_version, "current package version")
    observations = _observation_sets(manifest, set(indexed_contracts))
    final_version = observations[-1]["version"]
    if parsed_current != _release_version(final_version, "final observation version"):
        raise ValueError(f"current package version must select the final observation set {final_version}")
    effective = copy.deepcopy(indexed_contracts)
    for observation in observations:
        for override in observation["report_overrides"]:
            identity = (override["suite"], override["test_name"])
            effective[identity]["outcomes"]["report"] = copy.deepcopy(override["report"])
    return effective


def _ordered_prefix(base: object, current: object, label: str) -> None:
    if not isinstance(base, list) or not isinstance(current, list):
        raise ValueError(f"{label} must be lists")
    if len(current) < len(base):
        if label == "compatibility contracts":
            raise ValueError("compatibility contracts were removed")
        raise ValueError(f"{label} were removed or mutated")
    canonical_base = json.dumps(base, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
    canonical_prefix = json.dumps(
        current[: len(base)], sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    )
    if canonical_prefix != canonical_base:
        if label == "compatibility contracts":
            raise ValueError("compatibility contracts were mutated; historical rows must remain an ordered prefix")
        raise ValueError(f"{label} were removed or mutated; historical entries must remain an ordered prefix")


def validate_append_only_history(base_manifest: Manifest, current_manifest: Manifest) -> None:
    """Prove that contracts and observation sets retain exact content and order.

    Args:
        base_manifest: Manifest read from the authoritative event base.
        current_manifest: Manifest from the current checkout.

    Raises:
        ValueError: If history was removed, changed, reordered, or downgraded.
    """
    try:
        base_schema = _schema_version(base_manifest, "base manifest")
        current_schema = _schema_version(current_manifest, "current manifest")
    except ValueError as error:
        raise ValueError(f"unsupported manifest schema transition: {error}") from error
    if (base_schema, current_schema) not in {(1, 1), (1, 2), (2, 2)}:
        raise ValueError(f"unsupported manifest schema transition {base_schema} -> {current_schema}")
    _ordered_prefix(base_manifest.get("contracts"), current_manifest.get("contracts"), "compatibility contracts")
    if base_schema == OBSERVATION_SCHEMA_VERSION:
        _ordered_prefix(
            base_manifest.get("observation_sets"),
            current_manifest.get("observation_sets"),
            "observation sets",
        )
