"""Verified portable projection of approved caller and length research profiles."""

from __future__ import annotations

import hashlib
import logging
import os
from argparse import Namespace
from pathlib import Path
from typing import NoReturn, cast

from vntyper.scripts.calibration_artifact_io import write_checksums, write_json
from vntyper.scripts.calibration_caller_bundle import load_caller_model_bundle, validate_caller_payload
from vntyper.scripts.calibration_candidate import (
    CalibrationTarget,
    CandidateEnvelope,
    decode_candidate,
    validate_candidate_payload,
)
from vntyper.scripts.calibration_payload import (
    CallerBundleDescriptor,
    PayloadManifest,
    decode_caller_bundle_descriptor,
    decode_payload_manifest,
    validate_payload_observations,
)
from vntyper.scripts.calibration_portable_projection import build_portable_approval, portable_approval_document
from vntyper.scripts.calibration_secure_io import SecureDirectoryReader, read_regular_path
from vntyper.scripts.calibration_target_attestation import (
    decode_target_locked_attestation,
    decode_target_validation_attestation,
)
from vntyper.scripts.calibration_target_authority import decode_target_custodian_authority
from vntyper.scripts.calibration_target_completion import decode_target_completion
from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object
from vntyper.scripts.length_model_bundle import load_length_model_bundle

logger = logging.getLogger(__name__)
_LENGTH_FILES = {"length-model.json", "length-annotation.json"}
_CALLER_FILES = {"caller-bundle.json", "decision-profile.json", "advntr-policy.json", "background.json"}
_HEADERS = {"candidate.json", "payload-manifest.json"}


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _path(value: object, label: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value):
        _fail(f"calibration export {label} requires an explicit path")
    return Path(value)


def _reader(path: Path, names: set[str] | None, *, parent_descriptor: int | None = None) -> SecureDirectoryReader:
    # Root research directories may include outcomes and subdirectories, but only
    # their two named headers are opened. Payload directories must be exactly closed.
    if not getattr(os, "O_NOFOLLOW", 0) or not getattr(os, "O_DIRECTORY", 0):
        _fail("calibration export requires no-follow directory descriptors")
    descriptor = None
    try:
        descriptor = os.open(
            path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC, dir_fd=parent_descriptor
        )
        actual = frozenset(os.listdir(descriptor))
        if names is None:
            names = set(actual)
        elif not names <= actual:
            _fail("calibration export research profile headers are missing")
        return SecureDirectoryReader(path, descriptor, frozenset(names))
    except BaseException:
        if descriptor is not None:
            os.close(descriptor)
        raise


def _object(raw: bytes) -> dict[str, object]:
    document = load_strict_json_object(raw)
    if canonical_json_bytes(document) != raw:
        _fail("calibration export payload and candidate require canonical JSON bytes")
    return document


def _manifest(raw: bytes) -> PayloadManifest:
    # Reuse the strict duplicate-key/nonfinite parser for the manifest's array root.
    value = load_strict_json_object(b'{"value":' + raw + b"}")["value"]
    if canonical_json_bytes(value) != raw:
        _fail("calibration export payload manifest requires canonical JSON bytes")
    return decode_payload_manifest(value)


def _profile(
    path: Path, target: CalibrationTarget
) -> tuple[CandidateEnvelope, PayloadManifest, dict[str, bytes], dict[str, bytes], CallerBundleDescriptor | None]:
    try:
        with _reader(path, _HEADERS) as root:
            headers = root.read_files(tuple(sorted(_HEADERS)))
            candidate = decode_candidate(_object(headers["candidate.json"]))
            manifest = _manifest(headers["payload-manifest.json"])
            names = {item.path for item in manifest.files}
            if target == "length" and names != _LENGTH_FILES:
                _fail("calibration export length payload must contain exactly model and annotation")
            if target == "callers" and not {"caller-bundle.json", "decision-profile.json"} <= names <= _CALLER_FILES:
                _fail("calibration export caller payload file inventory differs")
            with _reader(Path("payload"), None, parent_descriptor=root.descriptor) as payload:
                if payload.names != frozenset(names):
                    _fail("calibration export payload directory differs from its exact manifest")
                files = payload.read_files(tuple(sorted(names)))
    except OSError as error:
        raise ValueError("calibration export profile or payload is missing, unreadable, or a symlink") from error
    validate_payload_observations(
        manifest, {name: (len(raw), hashlib.sha256(raw).hexdigest()) for name, raw in files.items()}
    )
    descriptor = decode_caller_bundle_descriptor(_object(files["caller-bundle.json"])) if target == "callers" else None
    validate_candidate_payload(candidate, manifest, expected_target=target, caller_descriptor=descriptor)
    return candidate, manifest, headers, files, descriptor


def export_calibration_bundle(args: Namespace, output: Path) -> bool:
    """Export only the verified payload and aggregate-free portable approval.

    The CLI owns ``atomic_output`` and supplies a new empty staging directory.
    Attestations are checked under the trusted-operator contract; standalone
    completion documents do not establish current custody-ledger state.

    Args:
        args: target, profile, validation, evaluation, authority, completion paths.
        output: Empty atomic staging directory, not a published destination.

    Returns:
        True after the complete portable inventory has been written and checked.

    Raises:
        ValueError: If inputs, observed bytes, semantics, approvals, or staging differ.
    """
    target = getattr(args, "target", None)
    if target not in ("length", "callers"):
        _fail("calibration export target must be length or callers")
    if not isinstance(output, Path) or output.is_symlink() or not output.is_dir() or any(output.iterdir()):
        _fail("calibration export requires an empty nonsymlink atomic staging directory")
    candidate, _, headers, files, descriptor = _profile(
        _path(getattr(args, "profile", None), "profile"), cast(CalibrationTarget, target)
    )
    documents = {
        field: load_strict_json_object(read_regular_path(_path(getattr(args, field, None), field)))
        for field in ("validation", "evaluation", "authority", "completion")
    }
    approval = build_portable_approval(
        candidate,
        decode_target_validation_attestation(documents["validation"]),
        decode_target_locked_attestation(documents["evaluation"]),
        decode_target_custodian_authority(documents["authority"]),
        decode_target_completion(documents["completion"]),
    )
    if target == "callers":
        validate_caller_payload(files, candidate, cast(CallerBundleDescriptor, descriptor))
    for name, raw in {**headers, **files}.items():
        (output / name).write_bytes(raw)
    write_json(output / "portable-approval.json", portable_approval_document(approval))
    write_checksums(output)
    if target == "length":
        load_length_model_bundle(output)
    else:
        load_caller_model_bundle(output)
    return True
