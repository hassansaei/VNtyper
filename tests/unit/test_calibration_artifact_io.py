"""Canonical calibration artifact I/O and checksum failures."""

from pathlib import Path

import pytest

from vntyper.scripts.calibration_artifact_io import (
    freeze_json,
    load_object,
    thaw_json,
    verify_checksums,
    write_checksums,
    write_json,
)

pytestmark = pytest.mark.unit


def test_canonical_json_and_direct_checksums_round_trip(tmp_path: Path) -> None:
    write_json(tmp_path / "nested" / "value.json", {"b": 2, "a": 1})
    write_checksums(tmp_path / "nested")

    assert load_object(tmp_path / "nested" / "value.json", "value") == {"a": 1, "b": 2}
    verify_checksums(tmp_path / "nested")


def test_checksum_tampering_fails_closed(tmp_path: Path) -> None:
    write_json(tmp_path / "value.json", {"value": 1})
    write_checksums(tmp_path)
    (tmp_path / "value.json").write_bytes(b"{}\n")

    with pytest.raises(ValueError, match="checksum differs"):
        verify_checksums(tmp_path)


def test_missing_or_invalid_json_has_a_path_aware_error(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="missing.*missing.json"):
        load_object(tmp_path / "missing.json", "missing value")


def test_json_freeze_and_thaw_are_recursive_and_canonicalizable() -> None:
    original = {"rows": [{"value": 1}], "label": "kept"}

    frozen = freeze_json(original)
    assert frozen["rows"] == ({"value": 1},)  # type: ignore[index]
    with pytest.raises(TypeError):
        frozen["label"] = "changed"  # type: ignore[index]
    assert thaw_json(frozen) == original
    assert freeze_json("scalar") == "scalar"
    assert thaw_json("scalar") == "scalar"
