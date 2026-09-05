"""Unit tests for atomic artifact publishing and cleanup (#314)."""

from __future__ import annotations

from pathlib import Path

import pytest

from vntyper.scripts.artifact_publish import (
    PARTIAL_SUFFIX,
    discard_partial,
    partial_output,
    partial_path,
    publish_partial,
)

pytestmark = pytest.mark.unit


def test_partial_path_is_deterministic_and_sibling(tmp_path: Path) -> None:
    final = tmp_path / "output_sliced.bam"
    partial = partial_path(final)
    assert partial == tmp_path / f"output_sliced.bam{PARTIAL_SUFFIX}"
    assert partial.parent == final.parent
    assert partial_path(str(final)) == partial


def test_publish_partial_replaces_atomically(tmp_path: Path) -> None:
    final = tmp_path / "output_sorted.bam"
    partial = partial_path(final)
    partial.write_bytes(b"BAM_DATA")

    published = publish_partial(partial, final)
    assert published == final
    assert final.read_bytes() == b"BAM_DATA"
    assert not partial.exists()


def test_publish_partial_raises_on_missing_source(tmp_path: Path) -> None:
    final = tmp_path / "target.bam"
    partial = partial_path(final)
    with pytest.raises(FileNotFoundError, match="Partial artifact not found"):
        publish_partial(partial, final)


def test_discard_partial_removes_existing_file_and_symlink(tmp_path: Path) -> None:
    partial = tmp_path / "file.bam.partial"
    partial.write_bytes(b"temp")
    discard_partial(partial)
    assert not partial.exists()

    # Repeated discard is a no-op
    discard_partial(partial)
    discard_partial(None)

    # Discarding a broken symlink unlinks the symlink itself
    target = tmp_path / "nonexistent"
    symlink = tmp_path / "symlink.bam.partial"
    symlink.symlink_to(target)
    assert symlink.is_symlink()
    discard_partial(symlink)
    assert not symlink.is_symlink()


def test_partial_output_publishes_on_clean_exit(tmp_path: Path) -> None:
    final = tmp_path / "result.bam"
    with partial_output(final) as partial:
        assert partial == partial_path(final)
        partial.write_bytes(b"VALID_BAM")

    assert final.exists()
    assert final.read_bytes() == b"VALID_BAM"
    assert not partial.exists()


def test_partial_output_discards_on_exception(tmp_path: Path) -> None:
    final = tmp_path / "result.bam"
    with pytest.raises(RuntimeError, match="stage exploded"), partial_output(final) as partial:
        partial.write_bytes(b"CORRUPT_OR_PARTIAL")
        raise RuntimeError("stage exploded")

    assert not final.exists()
    assert not partial.exists()


def test_partial_output_raises_when_partial_was_not_created(tmp_path: Path) -> None:
    final = tmp_path / "missing.bam"
    with pytest.raises(FileNotFoundError, match="Expected partial artifact not created"), partial_output(final):
        pass

    assert not final.exists()


def test_partial_output_checks_non_empty_when_requested(tmp_path: Path) -> None:
    final = tmp_path / "empty.bai"
    with (
        pytest.raises(OSError, match="Generated artifact is empty"),
        partial_output(final, check_non_empty=True) as partial,
    ):
        partial.touch()  # 0 bytes

    assert not final.exists()
    assert not partial.exists()
