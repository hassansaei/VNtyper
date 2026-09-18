"""Committed result bytes stay pinned during authorized extraction."""

from __future__ import annotations

import hashlib
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts.calibration_target_runs import TargetRunAsset

pytestmark = pytest.mark.unit


def commitment(path: Path, raw: bytes) -> TargetRunAsset:
    path.write_bytes(raw)
    return TargetRunAsset(path, hashlib.sha256(raw).hexdigest(), len(raw))


def test_reads_exact_committed_json_and_stream_verifies_large_inputs(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_target_asset_io import read_target_json, verify_target_asset

    asset = commitment(tmp_path / "result.json", b'{"value":17}\n')
    assert read_target_json(asset) == {"value": 17}
    assert verify_target_asset(asset) == asset.sha256


@pytest.mark.parametrize("raw", [b'{"value":17,"value":18}', b'{"value":NaN}', b"[]"])
def test_json_contract_does_not_coerce_or_overwrite_fields(tmp_path: Path, raw: bytes) -> None:
    from vntyper.scripts.calibration_target_asset_io import read_target_json

    with pytest.raises(ValueError):
        read_target_json(commitment(tmp_path / "result.json", raw))


def test_same_size_changed_bytes_and_size_changes_fail(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_target_asset_io import read_target_asset

    asset = commitment(tmp_path / "result.tsv", b"original")
    asset.path.write_bytes(b"modified")
    with pytest.raises(ValueError, match="digest"):
        read_target_asset(asset)
    asset.path.write_bytes(b"short")
    with pytest.raises(ValueError, match="size"):
        read_target_asset(asset)


def test_symlink_and_special_file_are_not_result_evidence(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_target_asset_io import read_target_asset

    asset = commitment(tmp_path / "result", b"valid")
    link = tmp_path / "link"
    link.symlink_to(asset.path)
    with pytest.raises((OSError, ValueError)):
        read_target_asset(TargetRunAsset(link, asset.sha256, asset.size_bytes))
    with pytest.raises(ValueError, match="regular"):
        read_target_asset(TargetRunAsset(tmp_path, asset.sha256, asset.size_bytes))


def test_predeclared_size_limit_precedes_open(tmp_path: Path) -> None:
    from vntyper.scripts.calibration_target_asset_io import read_target_asset

    asset = commitment(tmp_path / "result", b"123456")
    with (
        patch("vntyper.scripts.calibration_target_asset_io.os.open") as opened,
        pytest.raises(ValueError, match="limit"),
    ):
        read_target_asset(asset, maximum_bytes=5)
    opened.assert_not_called()


def test_path_replacement_during_read_is_detected(tmp_path: Path) -> None:
    from vntyper.scripts import calibration_target_asset_io as module

    asset = commitment(tmp_path / "result", b"original")
    real_read = module.os.read
    replaced = False

    def replacing_read(fd: int, size: int) -> bytes:
        nonlocal replaced
        raw = real_read(fd, size)
        if not replaced:
            replacement = tmp_path / "replacement"
            replacement.write_bytes(b"original")
            replacement.replace(asset.path)
            replaced = True
        return raw

    with patch.object(module.os, "read", side_effect=replacing_read), pytest.raises(ValueError, match="changed"):
        module.read_target_asset(asset)


def test_direct_typed_relative_path_is_refused_before_open() -> None:
    from vntyper.scripts.calibration_target_asset_io import read_target_asset

    asset = TargetRunAsset(Path("relative.json"), "a" * 64, 1)
    with (
        patch("vntyper.scripts.calibration_target_asset_io.os.open") as opened,
        pytest.raises(ValueError, match="absolute"),
    ):
        read_target_asset(asset)
    opened.assert_not_called()
