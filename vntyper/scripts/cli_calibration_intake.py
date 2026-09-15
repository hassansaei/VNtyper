"""Strict CLI adapter for local calibration intake bundle production."""

from __future__ import annotations

import argparse
import logging
import re
from collections.abc import Mapping
from pathlib import Path
from types import MappingProxyType
from typing import NoReturn

from vntyper.scripts.calibration_intake_io import PinnedCramReference, prepare_intake_bundle
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.canonical_json import load_strict_json_object

logger = logging.getLogger(__name__)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _exact_object(value: object, fields: set[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != fields:
        _fail(f"{label} fields differ from the closed schema")
    return value


def load_cram_references(path: Path | None) -> Mapping[str, PinnedCramReference]:
    """Decode the optional closed assembly-to-reference JSON contract.

    Args:
        path: Strict reference JSON path, or None when no CRAM inputs are declared.

    Returns:
        Immutable assembly bindings. The intake producer checks this is the exact
        assembly set required by the decoded intake.

    Raises:
        ValueError: If the path, JSON, fields, assembly names, paths, or hashes are invalid.
    """
    if path is None:
        return MappingProxyType({})
    if not isinstance(path, Path):
        _fail("calibration CRAM references must be supplied as a Path")
    root = _exact_object(
        load_strict_json_object(read_regular_path(path)),
        {"schema_version", "references"},
        "calibration CRAM reference document",
    )
    if root["schema_version"] != "calibration-cram-references-v1":
        _fail("calibration CRAM reference schema version must be calibration-cram-references-v1")
    values = root["references"]
    if not isinstance(values, Mapping):
        _fail("calibration CRAM references must be an object")
    references: dict[str, PinnedCramReference] = {}
    for assembly in sorted(values):
        if not isinstance(assembly, str) or not assembly or assembly.strip() != assembly:
            _fail("calibration CRAM reference assembly must be a non-empty trimmed string")
        row = _exact_object(values[assembly], {"path", "sha256"}, "calibration CRAM reference")
        reference_path = row["path"]
        digest = row["sha256"]
        if not isinstance(reference_path, str) or not reference_path or not Path(reference_path).is_absolute():
            _fail("calibration CRAM reference path must be an absolute non-empty string")
        if not isinstance(digest, str) or _SHA256.fullmatch(digest) is None:
            _fail("calibration CRAM reference SHA256 must be lowercase hexadecimal")
        references[assembly] = PinnedCramReference(Path(reference_path), digest)
    return MappingProxyType(references)


def run_calibration_intake(args: argparse.Namespace) -> None:
    """Translate parsed CLI arguments into one atomic intake producer call.

    Args:
        args: Parsed ``calibrate intake`` arguments.

    Raises:
        ValueError: If arguments or the optional CRAM reference document are invalid.
        RuntimeError: If input evidence changes or atomic publication is unavailable.
    """
    if not isinstance(args, argparse.Namespace):
        _fail("calibration intake CLI requires parsed arguments")
    manifest = getattr(args, "manifest", None)
    output = getattr(args, "output", None)
    priority = getattr(args, "preprocessing_priority", None)
    cram_reference_path = getattr(args, "cram_references", None)
    temporary = getattr(args, "temporary_directory", None)
    if not isinstance(manifest, Path) or not isinstance(output, Path):
        _fail("calibration intake manifest and output must be Paths")
    if not isinstance(priority, list):
        _fail("calibration intake preprocessing priority must be an explicit ordered list")
    if cram_reference_path is not None and not isinstance(cram_reference_path, Path):
        _fail("calibration CRAM references must be supplied as a Path")
    if temporary is not None and not isinstance(temporary, Path):
        _fail("calibration intake temporary directory must be a Path")
    references = load_cram_references(cram_reference_path)
    prepare_intake_bundle(
        manifest,
        output,
        preprocessing_priority=tuple(priority),
        cram_references=references,
        temporary_parent=temporary,
    )
