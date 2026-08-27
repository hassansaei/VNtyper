"""Unit tests for bounded, atomic reference downloads."""

from __future__ import annotations

import logging
import os
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier
from unittest.mock import patch

import pytest

from vntyper.scripts.install_references import download_file

pytestmark = pytest.mark.unit


class _StreamingResponse:
    """A bounded-read response double with optional mid-stream failure."""

    def __init__(
        self,
        chunks: list[bytes],
        *,
        fail_after: int | None = None,
        before_failure: Callable[[], None] | None = None,
    ) -> None:
        self._chunks = iter(chunks)
        self._fail_after = fail_after
        self._before_failure = before_failure
        self.read_sizes: list[int] = []

    def __enter__(self) -> _StreamingResponse:
        return self

    def __exit__(self, *exc_info: object) -> None:
        return None

    def read(self, size: int = -1) -> bytes:
        self.read_sizes.append(size)
        if self._fail_after is not None and len(self.read_sizes) > self._fail_after:
            if self._before_failure is not None:
                self._before_failure()
            raise OSError("connection reset after first chunk")
        return next(self._chunks, b"")


def test_download_streams_with_a_timeout_and_lands_via_os_replace(tmp_path: Path, monkeypatch) -> None:
    """Catch an unbounded response read or a direct write to the final destination."""
    from vntyper.scripts import install_references, reference_download

    response = _StreamingResponse([b"first ", b"second"])
    destination = tmp_path / "nested" / "reference.fa.gz"
    replacements: list[tuple[Path, Path]] = []

    def _replace(source: str | os.PathLike[str], target: str | os.PathLike[str]) -> None:
        source_path = Path(source)
        target_path = Path(target)
        assert source_path.read_bytes() == b"first second"
        assert not target_path.exists()
        replacements.append((source_path, target_path))
        source_path.rename(target_path)

    monkeypatch.setattr(reference_download.os, "replace", _replace)
    assert install_references.download_file is reference_download.download_file
    with patch("vntyper.scripts.reference_download.urlopen", return_value=response) as urlopen:
        assert download_file("https://example.invalid/reference.fa.gz", destination) is True

    assert destination.read_bytes() == b"first second"
    assert destination.stat().st_mode & 0o444 == 0o444, "the non-root image user must be able to read references"
    assert len(replacements) == 1
    partial_path, landed_path = replacements[0]
    assert landed_path == destination
    assert partial_path.parent == destination.parent
    assert partial_path.name.startswith(".reference.fa.gz.")
    assert partial_path.name.endswith(".partial")
    assert not partial_path.exists()
    assert len(response.read_sizes) == 3
    assert set(response.read_sizes) == {reference_download.DOWNLOAD_CHUNK_SIZE}
    urlopen.assert_called_once_with(
        "https://example.invalid/reference.fa.gz", timeout=reference_download.DOWNLOAD_TIMEOUT_SECONDS
    )


def test_an_existing_destination_is_reused_without_opening_the_url(tmp_path: Path) -> None:
    """Catch a regression that re-downloads multi-gigabyte references unconditionally."""
    destination = tmp_path / "reference.fa.gz"
    destination.write_bytes(b"existing")

    with patch("vntyper.scripts.reference_download.urlopen") as urlopen:
        assert download_file("https://example.invalid/reference.fa.gz", destination) is False

    assert destination.read_bytes() == b"existing"
    urlopen.assert_not_called()


@pytest.mark.parametrize("hostile_kind", ["file", "symlink"])
def test_a_predictable_partial_path_and_its_target_are_never_touched(tmp_path: Path, hostile_kind: str) -> None:
    """Catch reopening a predictable partial name that another process can pre-create."""
    destination = tmp_path / "reference.fa.gz"
    predictable_partial = destination.with_name("reference.fa.gz.partial")
    symlink_target = tmp_path / "owned-by-someone-else"
    symlink_target.write_bytes(b"target sentinel")
    if hostile_kind == "symlink":
        predictable_partial.symlink_to(symlink_target)
    else:
        predictable_partial.write_bytes(b"file sentinel")

    with patch("vntyper.scripts.reference_download.urlopen", return_value=_StreamingResponse([b"reference payload"])):
        assert download_file("https://example.invalid/reference.fa.gz", destination) is True

    assert destination.is_file()
    assert not destination.is_symlink()
    assert destination.read_bytes() == b"reference payload"
    assert symlink_target.read_bytes() == b"target sentinel"
    if hostile_kind == "symlink":
        assert predictable_partial.is_symlink()
        assert predictable_partial.resolve() == symlink_target
    else:
        assert predictable_partial.read_bytes() == b"file sentinel"


def test_concurrent_downloads_use_isolated_partial_files(tmp_path: Path, monkeypatch) -> None:
    """Catch two transfers sharing one partial pathname before either can publish."""
    from vntyper.scripts import reference_download

    destination = tmp_path / "reference.fa.gz"
    barrier = Barrier(2)
    payloads = (b"first payload", b"second payload")

    class _ConcurrentResponse(_StreamingResponse):
        def read(self, size: int = -1) -> bytes:
            if not self.read_sizes:
                barrier.wait(timeout=5)
            return super().read(size)

    responses = [_ConcurrentResponse([payload]) for payload in payloads]
    replacement_sources: list[Path] = []
    real_replace = os.replace

    def _replace(source: str | os.PathLike[str], target: str | os.PathLike[str]) -> None:
        replacement_sources.append(Path(source))
        real_replace(source, target)

    monkeypatch.setattr(reference_download.os, "replace", _replace)
    with (
        patch("vntyper.scripts.reference_download.urlopen", side_effect=responses),
        ThreadPoolExecutor(max_workers=2) as executor,
    ):
        results = list(
            executor.map(
                lambda _: download_file("https://example.invalid/reference.fa.gz", destination),
                range(2),
            )
        )

    assert results == [True, True]
    assert len(replacement_sources) == 2
    assert len(set(replacement_sources)) == 2
    assert all(not path.exists() for path in replacement_sources)
    assert destination.is_file()
    assert not destination.is_symlink()
    assert destination.read_bytes() in payloads


def test_a_mid_stream_failure_removes_the_partial_and_preserves_no_destination(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Catch cleanup that handles connection setup failures but leaves truncated bytes."""
    destination = tmp_path / "reference.fa.gz"
    observed_partials: list[Path] = []

    def _record_partial() -> None:
        observed_partials.extend(tmp_path.glob(".reference.fa.gz.*.partial"))

    response = _StreamingResponse([b"partial bytes"], fail_after=1, before_failure=_record_partial)
    expected = "Failed to download https://example.invalid/reference.fa.gz: connection reset after first chunk"

    with (
        patch("vntyper.scripts.reference_download.urlopen", return_value=response),
        caplog.at_level(logging.ERROR, logger="vntyper.scripts.reference_download"),
        pytest.raises(RuntimeError, match="connection reset after first chunk") as raised,
    ):
        download_file("https://example.invalid/reference.fa.gz", destination)

    assert str(raised.value) == expected
    assert expected in caplog.messages
    assert not destination.exists()
    assert len(observed_partials) == 1
    assert not observed_partials[0].exists()


def test_an_open_failure_raises_runtimeerror_with_the_logged_message(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Catch the legacy process exit and any mismatch between raised and logged text."""
    destination = tmp_path / "reference.fa.gz"
    observed_partials: list[Path] = []
    expected = "Failed to download https://example.invalid/reference.fa.gz: network unavailable"

    def _fail_open(*args: object, **kwargs: object) -> None:
        observed_partials.extend(tmp_path.glob(".reference.fa.gz.*.partial"))
        raise OSError("network unavailable")

    with (
        patch("vntyper.scripts.reference_download.urlopen", side_effect=_fail_open),
        caplog.at_level(logging.ERROR, logger="vntyper.scripts.reference_download"),
        pytest.raises(RuntimeError, match="network unavailable") as raised,
    ):
        download_file("https://example.invalid/reference.fa.gz", destination)

    assert str(raised.value) == expected
    assert expected in caplog.messages
    assert not destination.exists()
    assert len(observed_partials) == 1
    assert not observed_partials[0].exists()


def test_a_replace_failure_removes_its_exact_unique_partial(tmp_path: Path, monkeypatch) -> None:
    """Catch cleanup that misses the completed temp file when publication fails."""
    from vntyper.scripts import reference_download

    destination = tmp_path / "reference.fa.gz"
    observed_partials: list[Path] = []

    def _fail_replace(source: str | os.PathLike[str], target: str | os.PathLike[str]) -> None:
        observed_partials.append(Path(source))
        raise OSError("rename blocked")

    monkeypatch.setattr(reference_download.os, "replace", _fail_replace)
    with (
        patch(
            "vntyper.scripts.reference_download.urlopen",
            return_value=_StreamingResponse([b"complete payload"]),
        ),
        pytest.raises(RuntimeError, match="rename blocked"),
    ):
        download_file("https://example.invalid/reference.fa.gz", destination)

    assert len(observed_partials) == 1
    partial_path = observed_partials[0]
    assert partial_path.name.startswith(".reference.fa.gz.")
    assert partial_path.name.endswith(".partial")
    assert not partial_path.exists()
    assert not destination.exists()


def test_a_temp_allocation_failure_is_wrapped_without_attempting_cleanup(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Catch cleanup that dereferences an absent temp path and hides allocation failure."""
    destination = tmp_path / "reference.fa.gz"
    expected = "Failed to download https://example.invalid/reference.fa.gz: no space for temp file"

    with (
        patch(
            "vntyper.scripts.reference_download.tempfile.NamedTemporaryFile",
            side_effect=OSError("no space for temp file"),
        ),
        patch.object(Path, "unlink", autospec=True) as unlink,
        caplog.at_level(logging.ERROR, logger="vntyper.scripts.reference_download"),
        pytest.raises(RuntimeError, match="no space for temp file") as raised,
    ):
        download_file("https://example.invalid/reference.fa.gz", destination)

    assert str(raised.value) == expected
    assert expected in caplog.messages
    unlink.assert_not_called()
    assert not destination.exists()
