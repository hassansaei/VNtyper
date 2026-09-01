"""Canonical JSON and direct-child checksum I/O for calibration bundles."""

from __future__ import annotations

import hashlib
from collections.abc import Mapping
from pathlib import Path
from types import MappingProxyType

from vntyper.scripts.canonical_json import canonical_json_bytes, load_strict_json_object


def write_json(path: Path, value: object) -> None:
    """Write canonical JSON into an already-staged calibration tree."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_json_bytes(value))


def load_object(path: Path, label: str) -> dict[str, object]:
    """Load one strict JSON object with a path-aware validation error."""
    try:
        return load_strict_json_object(path.read_bytes())
    except (OSError, ValueError) as error:
        raise ValueError(f"{label} is missing, unreadable, or invalid: {path}") from error


def write_checksums(root: Path) -> None:
    """Write SHA-256 bindings for every direct file in one artifact directory."""
    files = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.iterdir())
        if path.is_file() and path.name != "checksums.json"
    }
    write_json(root / "checksums.json", {"schema_version": "calibration-checksums-v1", "files": files})


def verify_checksums(root: Path) -> None:
    """Verify a closed direct-child checksum manifest without traversing roles."""
    document = load_object(root / "checksums.json", "calibration checksums")
    if set(document) != {"schema_version", "files"} or document["schema_version"] != "calibration-checksums-v1":
        raise ValueError("calibration checksum fields or schema version differ")
    files = document["files"]
    if not isinstance(files, Mapping) or not files:
        raise ValueError("calibration checksum manifest must be a non-empty object")
    for name, expected in files.items():
        relative = Path(name)
        if relative.is_absolute() or len(relative.parts) != 1 or relative.name == "checksums.json":
            raise ValueError("calibration checksum paths must be direct artifact names")
        if not isinstance(expected, str) or len(expected) != 64:
            raise ValueError("calibration checksum digest must be SHA-256")
        try:
            observed = hashlib.sha256((root / relative).read_bytes()).hexdigest()
        except OSError as error:
            raise ValueError(f"calibration checksummed artifact is unreadable: {relative}") from error
        if observed != expected:
            raise ValueError(f"calibration artifact checksum differs: {relative}")


def freeze_json(value: object) -> object:
    """Recursively make decoded JSON mappings and arrays immutable."""
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): freeze_json(child) for key, child in value.items()})
    if isinstance(value, list):
        return tuple(freeze_json(child) for child in value)
    return value


def thaw_json(value: object) -> object:
    """Return plain canonicalizable mappings and arrays from frozen JSON."""
    if isinstance(value, Mapping):
        return {str(key): thaw_json(child) for key, child in value.items()}
    if isinstance(value, tuple):
        return [thaw_json(child) for child in value]
    return value
