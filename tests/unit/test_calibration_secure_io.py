"""Descriptor-pinned calibration import reads."""

import os
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts import calibration_secure_io
from vntyper.scripts.calibration_secure_io import SecureDirectoryReader, read_regular_path

pytestmark = pytest.mark.unit


def test_header_and_checksum_replacement_after_open_cannot_change_pinned_bytes(tmp_path: Path) -> None:
    authority = tmp_path / "authority.json"
    checksums = tmp_path / "checksums.json"
    authority.write_bytes(b"original authority")
    checksums.write_bytes(b"original checksums")
    replacement = tmp_path.parent / "replacement"
    replacement.write_bytes(b"replacement checksums")
    real_read = calibration_secure_io._read_descriptor
    calls = 0

    def replace_after_descriptors_open(descriptor: int) -> bytes:
        nonlocal calls
        calls += 1
        if calls == 1:
            os.replace(replacement, checksums)
        return real_read(descriptor)

    with (
        SecureDirectoryReader.open(tmp_path, {"authority.json", "checksums.json"}) as reader,
        patch("vntyper.scripts.calibration_secure_io._read_descriptor", side_effect=replace_after_descriptors_open),
    ):
        observed = reader.read_files(("authority.json", "checksums.json"))

    assert observed == {
        "authority.json": b"original authority",
        "checksums.json": b"original checksums",
    }
    assert checksums.read_bytes() == b"replacement checksums"


def test_secure_reads_reject_symlink_directories_and_files(tmp_path: Path) -> None:
    real = tmp_path / "real"
    real.mkdir()
    (real / "payload.json").write_bytes(b"payload")
    linked_root = tmp_path / "linked-root"
    linked_root.symlink_to(real, target_is_directory=True)
    linked_file = tmp_path / "linked-file"
    linked_file.symlink_to(real / "payload.json")

    with pytest.raises(ValueError, match="symlink"):
        SecureDirectoryReader.open(linked_root, {"payload.json"})
    with pytest.raises(ValueError, match="symlink"):
        read_regular_path(linked_file)
