"""Pinned calibration read sources keep one verified file identity."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pytest

from vntyper.scripts.calibration_read_sources import PinnedReadSource

pytestmark = pytest.mark.unit


def test_digest_copy_and_parser_duplicates_rewind_the_same_pinned_descriptor(tmp_path: Path) -> None:
    source_path = tmp_path / "reads"
    source_path.write_bytes(b"abcdef")
    copied = tmp_path / "copied"

    with PinnedReadSource.open(source_path, "synthetic source") as source:
        os.lseek(source.descriptor, 4, os.SEEK_SET)
        assert source.digest() == hashlib.sha256(b"abcdef").hexdigest()
        duplicate = source.duplicate_descriptor()
        try:
            assert os.read(duplicate, 3) == b"abc"
            assert os.lseek(source.descriptor, 0, os.SEEK_CUR) == 3
        finally:
            os.close(duplicate)
        with copied.open("wb") as target:
            assert source.copy_and_digest(target) == hashlib.sha256(b"abcdef").hexdigest()
        source.verify()

    assert copied.read_bytes() == b"abcdef"
    assert source.descriptor == -1


def test_replaced_path_is_rejected_while_the_original_descriptor_remains_pinned(tmp_path: Path) -> None:
    source_path = tmp_path / "reads"
    replacement = tmp_path / "replacement"
    source_path.write_bytes(b"original")
    replacement.write_bytes(b"alternate")

    with PinnedReadSource.open(source_path, "synthetic source") as source:
        source.digest()
        replacement.replace(source_path)
        with pytest.raises(RuntimeError, match="changed during read fingerprinting"):
            source.verify()


@pytest.mark.parametrize("kind", ["directory", "symlink"])
def test_nonregular_and_symlink_sources_are_rejected_without_following(kind: str, tmp_path: Path) -> None:
    source_path = tmp_path / "source"
    if kind == "directory":
        source_path.mkdir()
        message = "regular file"
    else:
        target = tmp_path / "target"
        target.write_bytes(b"private")
        source_path.symlink_to(target)
        message = "symlink"

    with pytest.raises(ValueError, match=message):
        PinnedReadSource.open(source_path, "synthetic source")
