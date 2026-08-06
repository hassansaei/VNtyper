"""The whole-request size boundary for the web service.

`docker/app/uploads.py` bounds what a submission may write into the job
directory, and `/run-job/` applies that bound while it copies each part. That
check runs inside the endpoint function, which the framework only calls once it
has read the request body and buffered every uploaded part. So the ceiling that
governs the volume says nothing about the bytes the request costs on the way
there.

`app.request_limits` closes that distance by counting the body as it arrives,
one ASGI message at a time, and answering 413 the moment the running total
passes the ceiling. The tests below are written so that a header-only check
cannot pass them:

* the ASGI-level tests assert on how much of the body the wrapped application
  was handed and how much was pulled off the stream at all -- an implementation
  that reads the whole body and then rejects it fails on both counts;
* the endpoint-level test streams its body from a generator that records what it
  yielded, so "the request was cut off part-way" is an assertion about the
  producer rather than about a status code.

That last test does not use the `client` fixture. Starlette's `TestClient` calls
`request.read()` inside its `receive`, so the whole body is materialised before
the application is ever invoked; a "stopped early" assertion made through it
would be measuring the test transport, not the service, and would pass against
an implementation that read everything first. `httpx.ASGITransport` pulls one
chunk per `receive`, so it is used directly there. It is the same in-process,
socketless arrangement `TestClient` is built on.

Everything without a body has to stay untouched, which is why the GET routes and
the file-streaming download endpoint are pinned here as well: this sits in front
of every request the service answers.

`docker` is put on `sys.path` by `tests/unit/web/conftest.py`, which pytest
imports before this module, so `app.request_limits` is importable here.
"""

import asyncio
import json
from collections.abc import AsyncIterator, Awaitable, Callable, MutableMapping
from pathlib import Path
from typing import Any
from uuid import uuid4

import httpx
import pytest

pytestmark = pytest.mark.unit

from app.config import settings  # noqa: E402
from app.request_limits import (  # noqa: E402
    REQUEST_FRAMING_ALLOWANCE,
    RequestSizeLimitMiddleware,
    configured_max_request_bytes,
)

# Small enough that the bodies below stay cheap, large enough that a chunk of
# them is a meaningful fraction of the ceiling.
CEILING = 4096
CHUNK = 1024

_BOUNDARY = "vntyperrequestsizeboundary"
_CONTENT_TYPE = f"multipart/form-data; boundary={_BOUNDARY}"


# ---------------------------------------------------------------------------
# Helpers: a synthetic ASGI environment, so the middleware can be driven with no
# framework, no transport and no test client in the way.
# ---------------------------------------------------------------------------


def _scope(method: str = "POST", headers: list[tuple[bytes, bytes]] | None = None) -> dict[str, Any]:
    """Build an HTTP ASGI scope.

    Args:
        method: HTTP method for the request.
        headers: Raw ASGI header pairs. Defaults to a chunked body, which is the
            case a declared-length check cannot answer.

    Returns:
        dict[str, Any]: A scope accepted by the middleware.
    """
    if headers is None:
        headers = [(b"transfer-encoding", b"chunked")]
    return {"type": "http", "method": method, "path": "/run-job/", "headers": headers}


def _stream(chunks: list[bytes]) -> tuple[Callable[[], Awaitable[dict[str, Any]]], list[int]]:
    """Build a `receive` that hands over `chunks` one ASGI message at a time.

    Args:
        chunks: Body pieces, in order.

    Returns:
        tuple: The `receive` callable, and a list that grows by one entry for
            every chunk actually pulled from it. The second element is what makes
            "stopped early" observable.
    """
    pulled: list[int] = []
    remaining = list(chunks)

    async def receive() -> dict[str, Any]:
        """Yield the next body message, then report disconnection.

        Returns:
            dict[str, Any]: An `http.request` message, or `http.disconnect` once
                the body is exhausted.
        """
        if not remaining:
            return {"type": "http.disconnect"}
        chunk = remaining.pop(0)
        pulled.append(len(chunk))
        return {"type": "http.request", "body": chunk, "more_body": bool(remaining)}

    return receive, pulled


# The stand-in applications and `_Recorder` below are spelled with
# `MutableMapping[str, Any]`, not `dict[str, Any]`, because that is the ASGI
# contract Starlette publishes and `request_limits` declares. A parameter is
# contravariant: a double that accepts only a `dict` is *narrower* than the
# callable it stands in for, so it would not be a legal argument even though it
# behaves correctly at runtime (#194).
def _draining_app(received: list[int]) -> Callable[..., Awaitable[None]]:
    """Build an ASGI app that reads the whole body, then answers 200.

    Args:
        received: List the app appends each delivered body length to.

    Returns:
        Callable: An ASGI application. It raises on `http.disconnect`, which is
            what Starlette does (`ClientDisconnect`) when a body read is cut off,
            so the middleware is exercised against the same control flow the real
            application produces.
    """

    async def app(scope: MutableMapping[str, Any], receive: Callable, send: Callable) -> None:
        """Drain the request body and reply.

        Args:
            scope: ASGI scope.
            receive: ASGI receive callable.
            send: ASGI send callable.

        Raises:
            RuntimeError: If the body read is cut off by a disconnect.
        """
        while True:
            message = await receive()
            if message["type"] == "http.disconnect":
                msg = "client disconnected while the body was being read"
                raise RuntimeError(msg)
            received.append(len(message.get("body", b"")))
            if not message.get("more_body"):
                break
        await send({"type": "http.response.start", "status": 200, "headers": []})
        await send({"type": "http.response.body", "body": b"handled"})

    return app


class _Recorder:
    """Collect the ASGI messages an application sends."""

    def __init__(self) -> None:
        self.messages: list[MutableMapping[str, Any]] = []

    async def __call__(self, message: MutableMapping[str, Any]) -> None:
        """Record one outbound message.

        Args:
            message: The ASGI message being sent.
        """
        self.messages.append(message)

    @property
    def status(self) -> int | None:
        """The status of the response that was started, if any.

        Returns:
            int | None: The status code, or None if no response was started.
        """
        for message in self.messages:
            if message["type"] == "http.response.start":
                return int(message["status"])
        return None

    @property
    def headers(self) -> dict[bytes, bytes]:
        """The headers of the response that was started.

        Returns:
            dict[bytes, bytes]: Lowercased header names mapped to their values.
        """
        for message in self.messages:
            if message["type"] == "http.response.start":
                return {name.lower(): value for name, value in message["headers"]}
        return {}

    @property
    def body(self) -> bytes:
        """The concatenated response body.

        Returns:
            bytes: Everything sent in `http.response.body` messages.
        """
        return b"".join(m.get("body", b"") for m in self.messages if m["type"] == "http.response.body")


# ---------------------------------------------------------------------------
# The ceiling itself.
# ---------------------------------------------------------------------------


def test_the_request_ceiling_is_derived_from_the_upload_ceiling() -> None:
    """One setting governs both bounds, so they cannot drift apart.

    The request is larger than the file it carries by the multipart framing
    around it, so the request ceiling is the upload ceiling plus a fixed
    allowance -- never a second, independently configured number.
    """
    assert configured_max_request_bytes() == settings.MAX_UPLOAD_BYTES + REQUEST_FRAMING_ALLOWANCE


def test_the_framing_allowance_is_generous_but_not_a_second_ceiling() -> None:
    """The allowance covers framing without meaningfully widening the bound.

    A submission to this service carries two file parts and a handful of short
    form fields, whose boundaries and part headers come to a few hundred bytes.
    The allowance has to exceed that comfortably while staying negligible next to
    the upload ceiling itself.
    """
    assert 4096 <= REQUEST_FRAMING_ALLOWANCE <= 1024 * 1024
    assert REQUEST_FRAMING_ALLOWANCE < settings.MAX_UPLOAD_BYTES // 1000


def test_a_shrunk_upload_ceiling_shrinks_the_request_ceiling(monkeypatch: pytest.MonkeyPatch) -> None:
    """The derivation is live, not captured at import time.

    Args:
        monkeypatch: Standard pytest fixture; restores the setting at teardown.
    """
    monkeypatch.setattr(settings, "MAX_UPLOAD_BYTES", 1234)

    assert configured_max_request_bytes() == 1234 + REQUEST_FRAMING_ALLOWANCE


# ---------------------------------------------------------------------------
# ASGI level: what the middleware does to a body while it is still arriving.
# ---------------------------------------------------------------------------


def test_a_body_over_the_ceiling_is_refused_before_it_has_all_been_read() -> None:
    """The refusal lands mid-body, not after the last byte.

    Two independent measurements say so: the wrapped application is never handed
    more than the ceiling, and the stream itself is pulled from only a chunk
    beyond it. An implementation that reads the whole body before deciding fails
    the second assertion even if it passes the first.
    """
    received: list[int] = []
    receive, pulled = _stream([b"x" * CHUNK] * 64)
    recorder = _Recorder()
    middleware = RequestSizeLimitMiddleware(_draining_app(received), max_bytes=CEILING)

    asyncio.run(middleware(_scope(), receive, recorder))

    assert recorder.status == 413
    assert sum(received) <= CEILING
    assert sum(pulled) <= CEILING + CHUNK
    assert sum(pulled) < 64 * CHUNK


def test_the_refusal_is_a_json_detail_like_every_other_error() -> None:
    """The 413 has the same shape as the service's other error responses."""
    receive, _ = _stream([b"x" * CHUNK] * 64)
    recorder = _Recorder()
    middleware = RequestSizeLimitMiddleware(_draining_app([]), max_bytes=CEILING)

    asyncio.run(middleware(_scope(), receive, recorder))

    assert recorder.headers[b"content-type"] == b"application/json"
    assert "maximum accepted size" in json.loads(recorder.body)["detail"]


def test_the_refusal_asks_for_the_connection_to_be_closed() -> None:
    """A client still sending must not be left holding the connection open.

    Nothing will read the rest of that body, so the response says so explicitly
    rather than leaving the server to reuse a connection whose request was never
    finished.
    """
    receive, _ = _stream([b"x" * CHUNK] * 64)
    recorder = _Recorder()
    middleware = RequestSizeLimitMiddleware(_draining_app([]), max_bytes=CEILING)

    asyncio.run(middleware(_scope(), receive, recorder))

    assert recorder.headers[b"connection"] == b"close"


def test_only_one_response_is_sent_when_the_ceiling_is_passed() -> None:
    """The wrapped application cannot answer on top of the refusal.

    Once the refusal is out, whatever the application goes on to send is dropped:
    two `http.response.start` messages on one connection is a protocol error.
    """
    recorder = _Recorder()

    async def stubborn_app(scope: MutableMapping[str, Any], receive: Callable, send: Callable) -> None:
        """Read until cut off, then answer anyway.

        Args:
            scope: ASGI scope.
            receive: ASGI receive callable.
            send: ASGI send callable.
        """
        while True:
            message = await receive()
            if message["type"] == "http.disconnect":
                break
            if not message.get("more_body"):
                break
        await send({"type": "http.response.start", "status": 200, "headers": []})
        await send({"type": "http.response.body", "body": b"handled"})

    receive, _ = _stream([b"x" * CHUNK] * 64)
    middleware = RequestSizeLimitMiddleware(stubborn_app, max_bytes=CEILING)

    asyncio.run(middleware(_scope(), receive, recorder))

    starts = [m for m in recorder.messages if m["type"] == "http.response.start"]
    assert len(starts) == 1
    assert starts[0]["status"] == 413


def test_a_declared_length_over_the_ceiling_is_refused_without_reading_a_byte() -> None:
    """A request that announces too much is answered before the body is touched.

    `receive` fails outright here, so a passing run proves the answer came from
    the declared length alone.
    """
    recorder = _Recorder()

    async def never_called() -> dict[str, Any]:
        """Fail if the body is read at all.

        Returns:
            dict[str, Any]: Never returns.

        Raises:
            AssertionError: Always.
        """
        raise AssertionError("the body was read despite an over-large declared length")

    middleware = RequestSizeLimitMiddleware(_draining_app([]), max_bytes=CEILING)
    scope = _scope(headers=[(b"content-length", str(CEILING + 1).encode())])

    asyncio.run(middleware(scope, never_called, recorder))

    assert recorder.status == 413
    assert recorder.headers[b"connection"] == b"close"


def test_a_declared_length_at_the_ceiling_is_allowed_through() -> None:
    """The ceiling itself is acceptable; only passing it is refused."""
    received: list[int] = []
    receive, _ = _stream([b"x" * CEILING])
    recorder = _Recorder()
    middleware = RequestSizeLimitMiddleware(_draining_app(received), max_bytes=CEILING)
    scope = _scope(headers=[(b"content-length", str(CEILING).encode())])

    asyncio.run(middleware(scope, receive, recorder))

    assert recorder.status == 200
    assert sum(received) == CEILING


def test_a_body_under_the_ceiling_arrives_intact() -> None:
    """Counting must not perturb what the application receives."""
    received: list[int] = []
    chunks = [b"a" * CHUNK, b"b" * CHUNK, b"c" * 17]
    receive, _ = _stream(chunks)
    recorder = _Recorder()
    middleware = RequestSizeLimitMiddleware(_draining_app(received), max_bytes=CEILING)

    asyncio.run(middleware(_scope(), receive, recorder))

    assert recorder.status == 200
    assert received == [len(chunk) for chunk in chunks]


def test_a_request_with_no_body_is_handed_the_untouched_stream() -> None:
    """Bodiless requests are passed through with nothing wrapped around them.

    A GET carries neither `Content-Length` nor `Transfer-Encoding`, so there is
    nothing to count. Handing the original `receive` and `send` straight to the
    application is what keeps a streaming response out of this middleware's way
    entirely.
    """
    seen: dict[str, Any] = {}
    recorder = _Recorder()

    async def app(scope: MutableMapping[str, Any], receive: Callable, send: Callable) -> None:
        """Record the callables it was handed.

        Args:
            scope: ASGI scope.
            receive: ASGI receive callable.
            send: ASGI send callable.
        """
        seen["receive"] = receive
        seen["send"] = send

    receive, _ = _stream([])
    middleware = RequestSizeLimitMiddleware(app, max_bytes=CEILING)

    asyncio.run(middleware(_scope(method="GET", headers=[]), receive, recorder))

    assert seen["receive"] is receive
    assert seen["send"] is recorder


@pytest.mark.parametrize("method", ["POST", "PUT", "PATCH", "DELETE"])
def test_a_body_arriving_with_no_framing_headers_is_still_counted(method: str) -> None:
    """The count is enabled by the method, not by a header the client chooses.

    `Content-Length` and `Transfer-Encoding` are the client's to send. Deciding
    whether to count on their presence makes the ceiling something the sender can
    switch off, and an ASGI server that delivered body messages without either
    would bypass it without the sender doing anything at all. Any method that may
    carry a body is counted, framing or none.

    Args:
        method: The HTTP method under test.
    """
    received: list[int] = []
    receive, pulled = _stream([b"x" * CHUNK] * 64)
    recorder = _Recorder()
    middleware = RequestSizeLimitMiddleware(_draining_app(received), max_bytes=CEILING)

    asyncio.run(middleware(_scope(method=method, headers=[]), receive, recorder))

    assert recorder.status == 413
    assert sum(received) <= CEILING
    assert sum(pulled) <= CEILING + CHUNK


def test_a_body_arriving_with_no_framing_headers_is_delivered_when_it_fits() -> None:
    """Counting an unframed body must not truncate a legitimate one.

    Without this, the test above would also pass against a middleware that
    refused every unframed request outright.
    """
    received: list[int] = []
    chunks = [b"a" * CHUNK, b"b" * 17]
    receive, _ = _stream(chunks)
    recorder = _Recorder()
    middleware = RequestSizeLimitMiddleware(_draining_app(received), max_bytes=CEILING)

    asyncio.run(middleware(_scope(method="POST", headers=[]), receive, recorder))

    assert recorder.status == 200
    assert received == [len(chunk) for chunk in chunks]


def test_a_get_that_declares_a_body_is_counted_too() -> None:
    """The header stays useful as an early answer where the method allows none.

    A GET is passed through when it says it has no body, because that is what
    keeps this middleware out of the way of the streaming downloads. One that
    announces a body is counted like anything else.
    """
    recorder = _Recorder()

    async def never_called() -> dict[str, Any]:
        """Fail if the body is read at all.

        Returns:
            dict[str, Any]: Never returns.

        Raises:
            AssertionError: Always.
        """
        raise AssertionError("the body was read despite an over-large declared length")

    middleware = RequestSizeLimitMiddleware(_draining_app([]), max_bytes=CEILING)
    scope = _scope(method="GET", headers=[(b"content-length", str(CEILING + 1).encode())])

    asyncio.run(middleware(scope, never_called, recorder))

    assert recorder.status == 413


def test_a_non_http_scope_is_passed_straight_through() -> None:
    """Lifespan and websocket traffic has no request body to bound."""
    seen: list[str] = []
    recorder = _Recorder()

    async def app(scope: MutableMapping[str, Any], receive: Callable, send: Callable) -> None:
        """Record the scope type it was called with.

        Args:
            scope: ASGI scope.
            receive: ASGI receive callable.
            send: ASGI send callable.
        """
        seen.append(scope["type"])

    receive, _ = _stream([])
    middleware = RequestSizeLimitMiddleware(app, max_bytes=CEILING)

    asyncio.run(middleware({"type": "lifespan"}, receive, recorder))

    assert seen == ["lifespan"]


def test_an_unparsable_declared_length_falls_through_to_counting() -> None:
    """A malformed header is not trusted either way; the bytes still decide."""
    received: list[int] = []
    receive, pulled = _stream([b"x" * CHUNK] * 64)
    recorder = _Recorder()
    middleware = RequestSizeLimitMiddleware(_draining_app(received), max_bytes=CEILING)
    scope = _scope(headers=[(b"content-length", b"not-a-number")])

    asyncio.run(middleware(scope, receive, recorder))

    assert recorder.status == 413
    assert sum(pulled) <= CEILING + CHUNK


def test_the_application_error_is_not_swallowed_when_the_ceiling_held() -> None:
    """Only a refusal suppresses the application; real failures still surface."""
    recorder = _Recorder()

    async def failing_app(scope: MutableMapping[str, Any], receive: Callable, send: Callable) -> None:
        """Fail while reading a body that is comfortably under the ceiling.

        Args:
            scope: ASGI scope.
            receive: ASGI receive callable.
            send: ASGI send callable.

        Raises:
            RuntimeError: Always.
        """
        await receive()
        msg = "the application itself failed"
        raise RuntimeError(msg)

    receive, _ = _stream([b"x" * 8])
    middleware = RequestSizeLimitMiddleware(failing_app, max_bytes=CEILING)

    with pytest.raises(RuntimeError, match="the application itself failed"):
        asyncio.run(middleware(_scope(), receive, recorder))


# ---------------------------------------------------------------------------
# Endpoint level: the bound has to be installed on the service, not merely
# available to it, and everything else has to be unaffected.
# ---------------------------------------------------------------------------


@pytest.fixture
def capped_request_app(web_app, monkeypatch: pytest.MonkeyPatch):
    """Shrink the service's request ceiling so test bodies can stay small.

    The upload ceiling is patched on the settings object the middleware derives
    its bound from, so the bound under test is the shipped derivation rather than
    a value the test injected directly.

    Args:
        web_app: Patched `app.main` module from conftest.
        monkeypatch: Standard pytest fixture; restores the setting at teardown.

    Returns:
        object: The same `app.main` module, with a shrunk request ceiling.
    """
    monkeypatch.setattr(settings, "MAX_UPLOAD_BYTES", CEILING)
    return web_app


def _multipart_prologue() -> bytes:
    """Build everything a submission sends before the BAM part's bytes.

    Split out from the body so an oversized submission can be streamed as valid
    multipart: a body of arbitrary bytes is rejected by the form parser as
    malformed long before its size matters, which would answer 400 and tell us
    nothing about the ceiling.

    Returns:
        bytes: The form fields and the file part's headers.
    """
    chunks: list[bytes] = []
    for name, value in (("thread", "1"), ("reference_assembly", "hg19")):
        chunks.append(f'--{_BOUNDARY}\r\nContent-Disposition: form-data; name="{name}"\r\n\r\n{value}\r\n'.encode())
    chunks.append(
        (
            f"--{_BOUNDARY}\r\n"
            'Content-Disposition: form-data; name="bam_file"; filename="sample.bam"\r\n'
            "Content-Type: application/octet-stream\r\n\r\n"
        ).encode()
    )
    return b"".join(chunks)


def _multipart_epilogue() -> bytes:
    """Build the closing boundary of a submission.

    Returns:
        bytes: Everything that follows the BAM part's bytes.
    """
    return f"\r\n--{_BOUNDARY}--\r\n".encode()


def _multipart(payload: bytes) -> bytes:
    """Build a `multipart/form-data` body carrying one BAM part.

    Args:
        payload: The file bytes to attach.

    Returns:
        bytes: A complete request body for `_CONTENT_TYPE`.
    """
    return _multipart_prologue() + payload + _multipart_epilogue()


async def _post_streamed(asgi_app: Any, body: AsyncIterator[bytes]) -> httpx.Response:
    """Submit a job whose body is produced lazily, one chunk per read.

    Args:
        asgi_app: The FastAPI application to call.
        body: An async iterator of body chunks. `httpx.ASGITransport` advances it
            once per `receive`, so the application controls how much of it is
            ever produced.

    Returns:
        httpx.Response: The service's response.
    """
    transport = httpx.ASGITransport(app=asgi_app)
    async with httpx.AsyncClient(transport=transport, base_url="http://testserver") as http_client:
        return await http_client.post("/run-job/", content=body, headers={"Content-Type": _CONTENT_TYPE})


def test_a_streamed_oversized_submission_is_cut_off_part_way(capped_request_app, tmp_path: Path) -> None:
    """The producer is left holding most of the body it meant to send.

    The body has no declared length, so nothing but a running count of the bytes
    read can refuse it, and the producer records what it managed to yield -- so
    "cut off part way" is measured at the source rather than inferred from the
    status code. See this module's docstring for why `TestClient` cannot make
    that measurement.

    Args:
        capped_request_app: `app.main` with the shrunk request ceiling.
        tmp_path: The directory the fixture configured as the input tree.
    """
    ceiling = configured_max_request_bytes()
    piece = b"y" * 8192
    total_pieces = (ceiling // len(piece)) * 4
    yielded: list[int] = []

    async def body() -> AsyncIterator[bytes]:
        """Yield a valid but oversized submission, recording progress.

        The part headers go first so the form parser stays satisfied and keeps
        reading; the refusal then has to come from the size, not from a parse
        error. The closing boundary is never reached.

        Yields:
            bytes: One piece of the request body.
        """
        prologue = _multipart_prologue()
        yielded.append(len(prologue))
        yield prologue
        for _ in range(total_pieces):
            yielded.append(len(piece))
            yield piece

    response = asyncio.run(_post_streamed(capped_request_app.app, body()))

    assert response.status_code == 413
    assert "maximum accepted size" in response.json()["detail"]
    assert sum(yielded) < total_pieces * len(piece), "the whole body was consumed before it was refused"
    assert sum(yielded) <= ceiling + len(piece)
    assert list((tmp_path / "input").iterdir()) == []
    capped_request_app.run_vntyper_job.delay.assert_not_called()


def test_a_streamed_submission_under_the_ceiling_is_delivered_whole(capped_request_app, tmp_path: Path) -> None:
    """The same lazy transport still delivers a legitimate body in full.

    Without this, the test above would also pass against a middleware that
    refused every streamed request.

    Args:
        capped_request_app: `app.main` with the shrunk request ceiling.
        tmp_path: The directory the fixture configured as the input tree.
    """
    payload = b"n" * 1024
    parts = [_multipart(payload)[i : i + 256] for i in range(0, len(_multipart(payload)), 256)]

    async def body() -> AsyncIterator[bytes]:
        """Yield a legitimate multipart body in small pieces.

        Yields:
            bytes: One piece of the request body.
        """
        for part in parts:
            yield part

    response = asyncio.run(_post_streamed(capped_request_app.app, body()))

    assert response.status_code == 200
    job_input_dir = tmp_path / "input" / response.json()["job_id"]
    assert (job_input_dir / "sample.bam").read_bytes() == payload


def test_a_declared_length_over_the_request_ceiling_is_refused(capped_request_app, client, tmp_path: Path) -> None:
    """The cheap header answer is wired in on the real service too.

    Args:
        capped_request_app: `app.main` with the shrunk request ceiling.
        client: TestClient fixture from conftest.
        tmp_path: The directory the fixture configured as the input tree.
    """
    response = client.post(
        "/run-job/",
        content=_multipart(b"x"),
        headers={"Content-Type": _CONTENT_TYPE, "Content-Length": str(configured_max_request_bytes() + 1)},
    )

    assert response.status_code == 413
    assert list((tmp_path / "input").iterdir()) == []
    capped_request_app.run_vntyper_job.delay.assert_not_called()


def test_a_submission_under_the_ceiling_still_completes(capped_request_app, client, tmp_path: Path) -> None:
    """The legitimate path is unchanged: the file lands and the job is enqueued.

    Args:
        capped_request_app: `app.main` with the shrunk request ceiling.
        client: TestClient fixture from conftest.
        tmp_path: The directory the fixture configured as the input tree.
    """
    payload = b"n" * 1024

    response = client.post(
        "/run-job/",
        files={"bam_file": ("sample.bam", payload, "application/octet-stream")},
        data={"thread": "1", "reference_assembly": "hg19"},
    )

    assert response.status_code == 200
    job_input_dir = tmp_path / "input" / response.json()["job_id"]
    assert (job_input_dir / "sample.bam").read_bytes() == payload
    capped_request_app.run_vntyper_job.delay.assert_called_once()


def test_get_routes_are_unaffected(capped_request_app, client) -> None:
    """Requests with no body are answered exactly as before.

    Args:
        capped_request_app: `app.main` with the shrunk request ceiling.
        client: TestClient fixture from conftest.
    """
    assert client.get("/health/").json() == {"status": "ok"}
    assert client.get("/version/").status_code == 200
    assert client.get("/download/does-not-exist/").status_code == 404


def test_the_download_endpoint_still_streams_a_file_larger_than_the_ceiling(
    capped_request_app, client, tmp_path: Path
) -> None:
    """Response bytes are not request bytes, and must not be counted as such.

    The archive served here is many times the request ceiling. A middleware that
    counted in the wrong direction, or that wrapped `send` for bodiless requests,
    would truncate or refuse it.

    Args:
        capped_request_app: `app.main` with the shrunk request ceiling.
        client: TestClient fixture from conftest.
        tmp_path: The directory the fixture configured as the output tree.
    """
    archive = b"r" * (configured_max_request_bytes() * 3)
    job_id = str(uuid4())
    (tmp_path / "output" / f"{job_id}.zip").write_bytes(archive)

    response = client.get(f"/download/{job_id}/")

    assert response.status_code == 200
    assert response.content == archive
