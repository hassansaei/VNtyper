"""Closed portable caller bundles verified before generated profiles can run."""

from __future__ import annotations

import hashlib
import logging
import os
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import NoReturn, cast

from vntyper.scripts.calibration_advntr_runtime_policy import (
    AdvntrRuntimePolicy,
    decode_advntr_runtime_policy,
    validate_advntr_runtime_policy,
)
from vntyper.scripts.calibration_caller_policy import CallerPolicyValues, decode_caller_policy_values
from vntyper.scripts.calibration_candidate import CandidateEnvelope, decode_candidate, validate_candidate_payload
from vntyper.scripts.calibration_payload import (
    CallerBundleDescriptor,
    PayloadManifest,
    decode_caller_bundle_descriptor,
    decode_payload_manifest,
    validate_payload_observations,
)
from vntyper.scripts.calibration_portable_background import validate_portable_background
from vntyper.scripts.calibration_portable_projection import (
    PortableApproval,
    decode_portable_approval,
    validate_portable_approval_candidate,
)
from vntyper.scripts.calibration_secure_io import SecureDirectoryReader
from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object
from vntyper.scripts.decision_profile import (
    ResolvedDecisionProfile,
    load_packaged_decision_profile,
    parse_decision_profile,
)

logger = logging.getLogger(__name__)
_HEADERS = {"candidate.json", "payload-manifest.json", "portable-approval.json", "checksums.json"}
_BASE = _HEADERS | {"caller-bundle.json", "decision-profile.json"}
_INVENTORIES = (_BASE, _BASE | {"advntr-policy.json"}, _BASE | {"advntr-policy.json", "background.json"})


@dataclass(frozen=True)
class CallerModelBundle:
    """Approved caller policy and immutable bytes of optional native background."""

    candidate: CandidateEnvelope
    payload: PayloadManifest
    approval: PortableApproval
    profile: ResolvedDecisionProfile
    caller_policy: CallerPolicyValues
    advntr_policy: AdvntrRuntimePolicy | None
    background_bytes: bytes | None
    sha256: str


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(raw: bytes) -> dict[str, object]:
    value = load_strict_json_object(raw)
    if canonical_json_bytes(value) != raw:
        _fail("runtime caller bundle requires canonical JSON bytes")
    return value


def _manifest(raw: bytes) -> PayloadManifest:
    value = load_strict_json_object(b'{"value":' + raw + b"}")["value"]
    if canonical_json_bytes(value) != raw:
        _fail("runtime caller bundle manifest requires canonical JSON bytes")
    return decode_payload_manifest(value)


def _checksums(raw: Mapping[str, bytes]) -> str:
    document = _object(raw["checksums.json"])
    if set(document) != {"schema_version", "files"} or document["schema_version"] != "calibration-checksums-v1":
        _fail("runtime caller bundle checksum fields or schema differ")
    values = document["files"]
    if not isinstance(values, Mapping) or set(values) != set(raw) - {"checksums.json"}:
        _fail("runtime caller bundle checksum inventory differs")
    for name, digest in values.items():
        if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
            _fail("runtime caller bundle checksum must be lowercase SHA-256")
        if hashlib.sha256(raw[name]).hexdigest() != digest:
            _fail("runtime caller bundle checksum differs from observed bytes")
    return hashlib.sha256(raw["checksums.json"]).hexdigest()


def validate_caller_payload(
    files: Mapping[str, bytes], candidate: CandidateEnvelope, descriptor: CallerBundleDescriptor
) -> tuple[ResolvedDecisionProfile, CallerPolicyValues, AdvntrRuntimePolicy | None]:
    """Validate profile composition and selected policy against the frozen candidate.

    Args:
        files: Verified exact payload bytes.
        candidate: Decoded candidate bound to that payload.
        descriptor: Decoded caller composition descriptor.

    Returns:
        Resolved decision profile, selected caller values, and optional native policy.

    Raises:
        ValueError: If profile, composition, partition, or native policy differs.
    """
    _object(files["decision-profile.json"])
    profile = parse_decision_profile(
        files["decision-profile.json"],
        packaged_document=load_packaged_decision_profile().document,
        allow_caller_generated=True,
    )
    document = _object(profile.canonical_bytes)
    metadata = document.get("generated_metadata")
    if document["schema_version"] != 2 or profile.profile_kind != "generated" or not isinstance(metadata, Mapping):
        _fail("calibration export caller payload requires a caller-generated schema-v2 decision profile")
    if (
        metadata["generation_target"] != "callers"
        or tuple(metadata["required_callers"]) != descriptor.required_callers
        or metadata["partition_manifest_hash"] != candidate.partition_sha256
    ):
        _fail("calibration export caller profile target, composition, or partition differs")
    inventory = cast(Mapping[str, Mapping[str, object]], document["inventory"])
    pointers = cast(list[str], metadata["generated_pointers"])
    policy = decode_caller_policy_values(
        {
            "schema_version": "calibration-caller-policy-values-v1",
            "required_callers": list(descriptor.required_callers),
            "values": {pointer: inventory[pointer]["value"] for pointer in pointers},
        }
    )
    runtime = None
    if "advntr" in descriptor.required_callers:
        runtime = decode_advntr_runtime_policy(_object(files["advntr-policy.json"]))
        background_sha256 = None
        if "background.json" in files:
            validate_portable_background(load_strict_json_object(files["background.json"]))
            background_sha256 = hashlib.sha256(files["background.json"]).hexdigest()
        validate_advntr_runtime_policy(runtime, policy, background_raw_sha256=background_sha256)

    return profile, policy, runtime


def load_caller_model_bundle(path: Path) -> CallerModelBundle:
    """Verify an exact approved portable caller bundle through pinned descriptors.

    Args:
        path: Portable directory with six, seven, or eight closed runtime files.

    Returns:
        Approved resolved caller settings and their transitive asset identities.

    Raises:
        ValueError: If inventory, canonical bytes, hashes, or approval bindings differ.
    """
    if not isinstance(path, Path):
        _fail("runtime caller bundle requires a Path")
    try:
        names = set(os.listdir(path))
        if names not in _INVENTORIES:
            _fail("runtime caller bundle file inventory differs")
        with SecureDirectoryReader.open(path, names) as reader:
            raw = reader.read_files(tuple(sorted(names)))
    except OSError as error:
        raise ValueError("runtime caller bundle is missing or unreadable") from error
    digest = _checksums(raw)
    candidate = decode_candidate(_object(raw["candidate.json"]))
    manifest = _manifest(raw["payload-manifest.json"])
    payload = {name: value for name, value in raw.items() if name not in _HEADERS}
    validate_payload_observations(
        manifest, {name: (len(value), hashlib.sha256(value).hexdigest()) for name, value in payload.items()}
    )
    descriptor = decode_caller_bundle_descriptor(_object(payload["caller-bundle.json"]))
    validate_candidate_payload(candidate, manifest, expected_target="callers", caller_descriptor=descriptor)
    approval = decode_portable_approval(_object(raw["portable-approval.json"]))
    validate_portable_approval_candidate(approval, candidate)
    profile, policy, native = validate_caller_payload(payload, candidate, descriptor)
    return CallerModelBundle(
        candidate, manifest, approval, profile, policy, native, payload.get("background.json"), digest
    )
