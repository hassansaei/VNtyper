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


def test_read_files_closes_new_descriptor_when_fstat_fails(tmp_path: Path) -> None:
    payload = tmp_path / "payload.json"
    payload.write_bytes(b"payload")

    with SecureDirectoryReader.open(tmp_path, {"payload.json"}) as reader:
        descriptors_before = frozenset(os.listdir("/proc/self/fd"))
        with (
            patch("vntyper.scripts.calibration_secure_io.os.fstat", side_effect=OSError("fstat failed")),
            pytest.raises(ValueError, match="changed|unreadable|symlink"),
        ):
            reader.read_file("payload.json")
        descriptors_after = frozenset(os.listdir("/proc/self/fd"))

    assert descriptors_after == descriptors_before


def test_read_files_rejects_duplicate_names_without_leaking_descriptors(tmp_path: Path) -> None:
    payload = tmp_path / "payload.json"
    payload.write_bytes(b"payload")

    with SecureDirectoryReader.open(tmp_path, {"payload.json"}) as reader:
        descriptors_before = frozenset(os.listdir("/proc/self/fd"))
        with pytest.raises(ValueError, match="duplicate"):
            reader.read_files(("payload.json", "payload.json"))
        descriptors_after = frozenset(os.listdir("/proc/self/fd"))

    assert descriptors_after == descriptors_before


def test_read_regular_path_refuses_to_degrade_without_no_follow_support(tmp_path: Path) -> None:
    """Mutation caught: a platform without O_NOFOLLOW silently follows symlinks on payload reads."""
    payload = tmp_path / "payload.json"
    payload.write_bytes(b"{}\n")

    with patch.object(calibration_secure_io, "_NOFOLLOW", 0), pytest.raises(ValueError, match="O_NOFOLLOW"):
        read_regular_path(payload)

    assert read_regular_path(payload) == b"{}\n"


def test_regular_path_and_pinned_child_request_nonblocking_admission(tmp_path: Path) -> None:
    payload = tmp_path / "payload.json"
    payload.write_bytes(b"{}\n")
    path_descriptor = os.open(payload, os.O_RDONLY)

    with patch.object(calibration_secure_io.os, "open", return_value=path_descriptor) as path_open:
        assert read_regular_path(payload) == b"{}\n"
    assert path_open.call_args.args[1] & os.O_NONBLOCK

    with SecureDirectoryReader.open(tmp_path, {"payload.json"}) as reader:
        real_open = os.open
        observed_flags = 0

        def record_open(path, flags, *args, **kwargs):
            nonlocal observed_flags
            observed_flags = flags
            return real_open(path, flags, *args, **kwargs)

        with patch.object(calibration_secure_io.os, "open", side_effect=record_open):
            assert reader.read_file("payload.json") == b"{}\n"
    assert observed_flags & os.O_NONBLOCK


def test_read_regular_path_rejects_a_fifo_without_a_writer(tmp_path: Path) -> None:
    fifo = tmp_path / "payload.fifo"
    os.mkfifo(fifo)

    with pytest.raises(ValueError, match="regular"):
        read_regular_path(fifo)


@pytest.mark.parametrize("directory_reader", [False, True])
def test_secure_reads_reject_in_place_mutation_during_chunked_read(tmp_path: Path, directory_reader: bool) -> None:
    """Pinned descriptors must not turn torn file contents into trusted evidence."""
    payload = tmp_path / "payload.json"
    chunk_size = 1024 * 1024
    payload.write_bytes(b"A" * (2 * chunk_size))
    before = payload.stat()
    real_read = os.read
    changed = False

    def mutate_after_first_chunk(descriptor: int, count: int) -> bytes:
        nonlocal changed
        data = real_read(descriptor, count)
        if data and not changed:
            changed = True
            with payload.open("r+b") as writer:
                writer.seek(chunk_size)
                writer.write(b"B" * chunk_size)
            os.utime(payload, ns=(before.st_atime_ns, before.st_mtime_ns + 1_000_000_000))
        return data

    with (
        patch.object(calibration_secure_io.os, "read", side_effect=mutate_after_first_chunk),
        pytest.raises(ValueError, match="changed during"),
    ):
        if directory_reader:
            with SecureDirectoryReader.open(tmp_path, {payload.name}) as reader:
                reader.read_file(payload.name)
        else:
            read_regular_path(payload)
    assert changed


def test_secure_read_rejects_mutation_even_when_mtime_is_restored(tmp_path: Path) -> None:
    """Restoring modification time must not hide a changed inode change time."""
    payload = tmp_path / "payload.json"
    payload.write_bytes(b"original")
    metadata = payload.stat()
    changed_metadata = type(
        "ChangedMetadata",
        (),
        {
            "st_dev": metadata.st_dev,
            "st_ino": metadata.st_ino,
            "st_size": metadata.st_size,
            "st_mtime_ns": metadata.st_mtime_ns,
            "st_ctime_ns": metadata.st_ctime_ns + 1,
            "st_mode": metadata.st_mode,
        },
    )()
    with (
        patch.object(calibration_secure_io.os, "fstat", side_effect=[metadata, metadata, changed_metadata]),
        pytest.raises(ValueError, match="changed during"),
    ):
        read_regular_path(payload)


@pytest.mark.parametrize("data", [b"", b"short"])
def test_secure_read_rejects_premature_end_of_file(tmp_path: Path, data: bytes) -> None:
    payload = tmp_path / "payload.json"
    payload.write_bytes(b"complete evidence")
    with (
        patch.object(calibration_secure_io.os, "read", side_effect=[data, b""]),
        pytest.raises(ValueError, match="changed during"),
    ):
        read_regular_path(payload)


def test_secure_read_accepts_a_stable_empty_file(tmp_path: Path) -> None:
    payload = tmp_path / "payload.json"
    payload.write_bytes(b"")
    assert read_regular_path(payload) == b""
