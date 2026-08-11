"""Shared fixtures for tests/unit/."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from vntyper.scripts.reference_provenance import build_record
from vntyper.scripts.reference_provenance import merge as merge_provenance

_REAL_INSTALL_CONFIG = Path(__file__).resolve().parents[2] / "vntyper" / "scripts" / "install_references_config.json"


def _touch(path: Path) -> None:
    """Create an empty file at `path`, making any missing parent directories."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()


@pytest.fixture
def install_config(tmp_path: Path) -> dict:
    """The real, shipped `install_references_config.json`, staged over `tmp_path`.

    Loads the shipped config rather than a hand-written stub, so a drift between the
    writer's schema and the one VNtyper actually ships cannot go unnoticed. `tmp_path`
    is the same directory pytest hands to the test itself (fixtures that request the
    same built-in fixture share the one instance for that test node), and every file
    `canonical_reference_keys` can name - each genome's `installed_path`, every
    `common_references` entry, every derivation's `output` - is created here (empty)
    and given a provenance record (see `reference_provenance`), so this fixture models
    a *verified* complete install, not merely a present one. `canonical_reference_keys`
    now requires a provenance record for a path, not just that it exists (the fix for
    the parked finding this fixture's own docstring used to warn about: a fixture that
    touches every path but records no provenance is exactly the shape that let an
    unverified retained file get blessed into `config.json` unnoticed); without staging
    both, every test built on this fixture would see an empty result on a fresh
    `tmp_path` regardless of what the config says.

    Args:
        tmp_path: Pytest's per-test temporary directory.

    Returns:
        dict: The parsed install config, unmodified.
    """
    config = json.loads(_REAL_INSTALL_CONFIG.read_text(encoding="utf-8"))
    records: dict[str, dict] = {}

    def _touch_and_record(relative: str) -> None:
        _touch(tmp_path / relative)
        records[relative] = build_record(sha256="0" * 64, source="from-source", source_url="https://example.com/test")

    for section in ("ucsc_references", "ncbi_references", "ensembl_references"):
        for entry in config.get(section, {}).values():
            _touch_and_record(entry["installed_path"])

    for entry in config.get("common_references", []):
        _touch_and_record(entry["installed_path"])

    for spec in config.get("derivations", []):
        _touch_and_record(spec["output"])

    merge_provenance(tmp_path, records)

    return config
