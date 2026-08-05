"""The upload-size boundary for the web service.

`/run-job/` writes two client-supplied files into a directory on the volume that
every job's input and output share. The rate limiter bounds how many requests a
client may make, not how many bytes each one carries, so the byte count needs a
limit of its own.

Two things are checked, and the second is the one that matters:

1. The declared `Content-Length` is refused when it already exceeds the limit.
   That header is chosen by the client, so it is a cheap early answer and
   nothing more.
2. The bytes actually read are counted while they are copied, and the copy stops
   the moment the running total passes the limit. A rejected upload must leave
   nothing behind, so the assertions below look at the filesystem rather than
   only at the status code.

Each endpoint test states which of the two paths it exercises, because a test
that cannot tell them apart would pass against a header-only check.

`docker` is put on `sys.path` by `tests/unit/web/conftest.py`, which pytest
imports before this module, so `app.uploads` is importable here.
"""

import io
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

from app.config import settings  # noqa: E402
from app.uploads import UPLOAD_CHUNK_SIZE, save_upload_bounded  # noqa: E402

# Small enough to keep test bodies tiny, comfortably larger than the few hundred
# bytes of multipart framing that ride along with them in the request body.
CAP = 4096

_BOUNDARY = "vntyperuploadlimitboundary"
_CONTENT_TYPE = f"multipart/form-data; boundary={_BOUNDARY}"


# ---------------------------------------------------------------------------
# The bounded copy on its own, with no app in the way.
# ---------------------------------------------------------------------------


def test_copy_under_the_cap_writes_every_byte(tmp_path: Path) -> None:
    """A source below the limit is copied through unchanged.

    Args:
        tmp_path: Scratch directory standing in for the job input directory.
    """
    destination = tmp_path / "sample.bam"
    payload = b"a" * 100

    written = save_upload_bounded(io.BytesIO(payload), str(destination), 1000)

    assert written == 100
    assert destination.read_bytes() == payload


def test_copy_of_exactly_the_cap_is_accepted(tmp_path: Path) -> None:
    """The limit itself is allowed; only passing it is refused.

    Args:
        tmp_path: Scratch directory standing in for the job input directory.
    """
    destination = tmp_path / "sample.bam"
    payload = b"b" * 64

    written = save_upload_bounded(io.BytesIO(payload), str(destination), 64)

    assert written == 64
    assert destination.read_bytes() == payload


def test_copy_one_byte_over_the_cap_is_refused(tmp_path: Path) -> None:
    """A single byte past the limit is enough to refuse the whole copy.

    Args:
        tmp_path: Scratch directory standing in for the job input directory.
    """
    destination = tmp_path / "sample.bam"

    with pytest.raises(ValueError):
        save_upload_bounded(io.BytesIO(b"c" * 65), str(destination), 64)


def test_refused_copy_leaves_no_partial_file(tmp_path: Path) -> None:
    """The partially written file is removed, not left on the volume.

    This is the assertion that separates a real cap from one that merely stops
    writing: refusing an upload has to reclaim what it already wrote.

    Args:
        tmp_path: Scratch directory standing in for the job input directory.
    """
    destination = tmp_path / "sample.bam"
    oversized = b"d" * (UPLOAD_CHUNK_SIZE * 3)

    with pytest.raises(ValueError):
        save_upload_bounded(io.BytesIO(oversized), str(destination), UPLOAD_CHUNK_SIZE)

    assert not destination.exists()
    assert list(tmp_path.iterdir()) == []


def test_copy_reads_in_bounded_chunks_never_the_whole_source(tmp_path: Path) -> None:
    """The source is drained a chunk at a time, so its size never sets the memory used.

    A `read()` with no argument would return the entire upload, which defeats
    the point of counting while copying.

    Args:
        tmp_path: Scratch directory standing in for the job input directory.
    """
    requested: list[int] = []

    class RecordingSource:
        """A byte source that remembers how much of it was asked for at a time."""

        def __init__(self, data: bytes) -> None:
            self._buffer = io.BytesIO(data)

        def read(self, size: int = -1) -> bytes:
            """Record the requested size and delegate to the buffer.

            Args:
                size: Number of bytes asked for.

            Returns:
                bytes: Up to `size` bytes from the buffer.
            """
            requested.append(size)
            return self._buffer.read(size)

    save_upload_bounded(RecordingSource(b"e" * (UPLOAD_CHUNK_SIZE * 2)), str(tmp_path / "s.bam"), CAP * 1000)

    assert requested, "the source was never read"
    assert all(0 < size <= UPLOAD_CHUNK_SIZE for size in requested)


def test_error_message_does_not_echo_the_upload(tmp_path: Path) -> None:
    """The refusal explains the limit without quoting what was sent.

    Args:
        tmp_path: Scratch directory standing in for the job input directory.
    """
    with pytest.raises(ValueError) as excinfo:
        save_upload_bounded(io.BytesIO(b"f" * 200), str(tmp_path / "s.bam"), 10)

    assert "maximum accepted size" in str(excinfo.value)


# ---------------------------------------------------------------------------
# The configured limit.
# ---------------------------------------------------------------------------


def test_configured_limit_is_a_documented_default() -> None:
    """Pin the shipped ceiling so changing it stays a deliberate act.

    The documented workflow uploads a MUC1-region subset measured in megabytes,
    so a gibibyte leaves three orders of magnitude of headroom while still
    bounding what one request can write.
    """
    assert settings.MAX_UPLOAD_BYTES == 1024 * 1024 * 1024


# ---------------------------------------------------------------------------
# Endpoint level: the limit has to be wired in, not merely available.
# ---------------------------------------------------------------------------


@pytest.fixture
def capped_app(web_app, monkeypatch: pytest.MonkeyPatch):
    """Shrink the service's upload ceiling so test bodies can stay small.

    Args:
        web_app: Patched `app.main` module from conftest.
        monkeypatch: Standard pytest fixture; restores the ceiling at teardown.

    Returns:
        object: The same `app.main` module, with a `CAP`-byte upload ceiling.
    """
    monkeypatch.setattr(web_app, "MAX_UPLOAD_BYTES", CAP, raising=False)
    return web_app


def _job_files(input_dir: Path) -> list[Path]:
    """List every file written under the job input tree.

    Args:
        input_dir: The service's configured input directory.

    Returns:
        list[Path]: Files (not directories) found anywhere beneath it.
    """
    return [path for path in input_dir.rglob("*") if path.is_file()]


def _multipart(files: list[tuple[str, str, bytes]]) -> bytes:
    """Build a `multipart/form-data` body by hand.

    The bodies here have to carry a `Content-Length` that disagrees with what is
    actually sent, or none at all, which the high-level `files=` helper will not
    do.

    Args:
        files: `(field_name, filename, payload)` triples to attach.

    Returns:
        bytes: A complete request body for `_CONTENT_TYPE`.
    """
    chunks: list[bytes] = []
    for name, value in (("thread", "1"), ("reference_assembly", "hg19")):
        chunks.append(f'--{_BOUNDARY}\r\nContent-Disposition: form-data; name="{name}"\r\n\r\n{value}\r\n'.encode())
    for field, filename, payload in files:
        header = (
            f"--{_BOUNDARY}\r\n"
            f'Content-Disposition: form-data; name="{field}"; filename="{filename}"\r\n'
            "Content-Type: application/octet-stream\r\n\r\n"
        )
        chunks.append(header.encode() + payload + b"\r\n")
    chunks.append(f"--{_BOUNDARY}--\r\n".encode())
    return b"".join(chunks)


def test_declared_length_over_the_cap_is_refused(capped_app, client, tmp_path: Path) -> None:
    """Exercises the header path only.

    The body carries a single byte, so the running byte count cannot possibly
    trip; a rejection here can only have come from the declared length.

    Args:
        capped_app: `app.main` with the shrunk ceiling.
        client: TestClient fixture from conftest.
        tmp_path: The directory the fixture configured as the input tree.
    """
    response = client.post(
        "/run-job/",
        content=_multipart([("bam_file", "sample.bam", b"x")]),
        headers={"Content-Type": _CONTENT_TYPE, "Content-Length": str(CAP * 1000)},
    )

    assert response.status_code == 413
    assert "maximum accepted size" in response.json()["detail"]
    assert _job_files(tmp_path / "input") == []
    capped_app.run_vntyper_job.delay.assert_not_called()
    capped_app.run_vntyper_job.apply_async.assert_not_called()


def test_body_over_the_cap_without_a_declared_length_is_stopped_mid_copy(capped_app, client, tmp_path: Path) -> None:
    """Exercises the streaming path only.

    The body is sent with chunked transfer encoding, so there is no
    `Content-Length` for a header check to look at. Only a cap enforced on the
    bytes actually read can refuse this.

    Args:
        capped_app: `app.main` with the shrunk ceiling.
        client: TestClient fixture from conftest.
        tmp_path: The directory the fixture configured as the input tree.
    """
    body = _multipart([("bam_file", "sample.bam", b"y" * (CAP * 2))])

    def stream():
        """Yield the body in pieces so httpx omits `Content-Length`."""
        yield body

    response = client.post("/run-job/", content=stream(), headers={"Content-Type": _CONTENT_TYPE})

    assert response.status_code == 413
    assert _job_files(tmp_path / "input") == []
    capped_app.run_vntyper_job.delay.assert_not_called()


def test_body_that_understates_its_length_is_stopped_mid_copy(capped_app, client, tmp_path: Path) -> None:
    """Exercises the streaming path only.

    The declared length is well under the ceiling while the body sent is well
    over it, so a check that trusts the header lets this through.

    Args:
        capped_app: `app.main` with the shrunk ceiling.
        client: TestClient fixture from conftest.
        tmp_path: The directory the fixture configured as the input tree.
    """
    response = client.post(
        "/run-job/",
        content=_multipart([("bam_file", "sample.bam", b"z" * (CAP * 2))]),
        headers={"Content-Type": _CONTENT_TYPE, "Content-Length": "5"},
    )

    assert response.status_code == 413
    assert _job_files(tmp_path / "input") == []


def test_rejected_upload_leaves_no_job_directories_behind(capped_app, client, tmp_path: Path) -> None:
    """The empty directories the request created are reclaimed too.

    Args:
        capped_app: `app.main` with the shrunk ceiling.
        client: TestClient fixture from conftest.
        tmp_path: The directory the fixture configured as the input tree.
    """
    response = client.post(
        "/run-job/",
        content=_multipart([("bam_file", "sample.bam", b"z" * (CAP * 2))]),
        headers={"Content-Type": _CONTENT_TYPE, "Content-Length": "5"},
    )

    assert response.status_code == 413
    assert list((tmp_path / "input").iterdir()) == []
    assert list((tmp_path / "output").iterdir()) == []


def test_oversized_index_upload_is_refused_as_well(capped_app, client, tmp_path: Path) -> None:
    """The index slot is bounded too, not only the alignment slot.

    The alignment part is tiny and the index part is over the ceiling, so the
    refusal can only have come from bounding the index copy. The alignment file
    written just before it must not survive.

    Args:
        capped_app: `app.main` with the shrunk ceiling.
        client: TestClient fixture from conftest.
        tmp_path: The directory the fixture configured as the input tree.
    """
    response = client.post(
        "/run-job/",
        content=_multipart(
            [
                ("bam_file", "sample.bam", b"a" * 16),
                ("bai_file", "sample.bam.bai", b"b" * (CAP * 2)),
            ]
        ),
        headers={"Content-Type": _CONTENT_TYPE, "Content-Length": "5"},
    )

    assert response.status_code == 413
    assert _job_files(tmp_path / "input") == []
    capped_app.run_vntyper_job.delay.assert_not_called()


def test_upload_under_the_cap_still_succeeds_end_to_end(capped_app, client, tmp_path: Path) -> None:
    """The legitimate path is unchanged: both files land and the job is enqueued.

    Args:
        capped_app: `app.main` with the shrunk ceiling.
        client: TestClient fixture from conftest.
        tmp_path: The directory the fixture configured as the input tree.
    """
    bam_payload = b"n" * 1024
    bai_payload = b"i" * 32

    response = client.post(
        "/run-job/",
        files={
            "bam_file": ("sample.bam", bam_payload, "application/octet-stream"),
            "bai_file": ("sample.bam.bai", bai_payload, "application/octet-stream"),
        },
        data={"thread": "1", "reference_assembly": "hg19"},
    )

    assert response.status_code == 200
    job_input_dir = tmp_path / "input" / response.json()["job_id"]
    assert (job_input_dir / "sample.bam").read_bytes() == bam_payload
    assert (job_input_dir / "sample.bam.bai").read_bytes() == bai_payload
    capped_app.run_vntyper_job.delay.assert_called_once()


def test_upload_just_under_the_cap_is_accepted(capped_app, client, tmp_path: Path) -> None:
    """A payload that approaches the ceiling without reaching it is still stored.

    Args:
        capped_app: `app.main` with the shrunk ceiling.
        client: TestClient fixture from conftest.
        tmp_path: The directory the fixture configured as the input tree.
    """
    payload = b"m" * (CAP - 512)

    response = client.post(
        "/run-job/",
        content=_multipart([("bam_file", "sample.bam", payload)]),
        headers={"Content-Type": _CONTENT_TYPE, "Content-Length": str(CAP - 1)},
    )

    assert response.status_code == 200
    job_input_dir = tmp_path / "input" / response.json()["job_id"]
    assert (job_input_dir / "sample.bam").read_bytes() == payload
