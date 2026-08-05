"""A ceiling on the size of a request, applied while the request is arriving.

`uploads.save_upload_bounded` bounds what a submission may write into the job
directory, and `/run-job/` applies that bound as it copies each part. That check
runs inside the endpoint function, and by the time the framework calls an
endpoint it has already read the request body off the connection and buffered
every uploaded part. So the endpoint's ceiling governs the job volume and
nothing else: the bytes a request costs on the way to the endpoint are not
covered by it.

This module covers them. `RequestSizeLimitMiddleware` is plain ASGI -- it wraps
`receive` and adds up the body one message at a time, so the running total is
known while the body is still arriving rather than only once it has all been
read. The moment the total passes the ceiling the request is answered with 413
and the application is told the client is gone, so nothing further is read and
the endpoint never runs.

It is deliberately *not* a `BaseHTTPMiddleware` subclass. That base class reads
the request body itself before handing it on, which is precisely the cost this
module exists to avoid.

The ceiling is derived from `settings.MAX_UPLOAD_BYTES` rather than configured
separately, so the two bounds cannot drift apart. A request is larger than the
file it carries -- multipart boundaries, per-part headers and the submission's
other form fields all ride along -- so a fixed allowance is added on top. The
allowance is orders of magnitude more than this API's form actually needs and
orders of magnitude less than the upload ceiling, which is what lets one number
govern both bounds without either being wrong.

Requests with no body are handed the original `receive` and `send` untouched:
there is nothing to count, and passing the callables straight through keeps this
middleware out of the way of streamed responses such as the result downloads.
"""

import json
import logging
from collections.abc import Awaitable, Callable, MutableMapping
from typing import Any

from .config import settings

logger = logging.getLogger(__name__)

# ASGI's own vocabulary, spelled out so the signatures below read as protocol
# rather than as `dict`/`Callable` soup. Starlette exports equivalents, but they
# are typing-only aliases and this module does not otherwise need Starlette.
Message = MutableMapping[str, Any]
Scope = MutableMapping[str, Any]
Receive = Callable[[], Awaitable[Message]]
Send = Callable[[Message], Awaitable[None]]
ASGIApp = Callable[[Scope, Receive, Send], Awaitable[None]]

# How much larger than its payload a submission is allowed to be. A request to
# this service carries at most two file parts and a handful of short form
# fields, whose boundary lines and part headers come to a few hundred bytes; 64
# KiB is two orders of magnitude of headroom for that, while being a
# ten-thousandth of the default upload ceiling. It exists so that a file of
# exactly the permitted size is not refused for the framing wrapped around it --
# not as a second, independently meaningful limit.
REQUEST_FRAMING_ALLOWANCE = 64 * 1024

# A request carries a body if and only if it says so with one of these
# (RFC 9112 section 6). Anything else -- every GET this service answers, the
# health check, the download routes -- has no body to bound.
BODY_HEADERS: tuple[bytes, ...] = (b"content-length", b"transfer-encoding")


def configured_max_request_bytes() -> int:
    """Return the whole-request ceiling implied by the configured upload ceiling.

    Read from `settings` on every call rather than captured at import time, so
    the bound the middleware applies is always the one currently configured.

    Returns:
        int: The largest request body accepted, in bytes.
    """
    return settings.MAX_UPLOAD_BYTES + REQUEST_FRAMING_ALLOWANCE


def _declared_length(headers: list[tuple[bytes, bytes]]) -> int | None:
    """Read the body length a request declares for itself.

    Args:
        headers: Raw ASGI header pairs.

    Returns:
        int | None: The declared length, or None when it is absent or not a
            plain number. A value that cannot be read is not trusted in either
            direction -- the running byte count still decides.
    """
    for name, value in headers:
        if name.lower() == b"content-length":
            text = value.decode("latin-1").strip()
            return int(text) if text.isdigit() else None
    return None


def _has_body(headers: list[tuple[bytes, bytes]]) -> bool:
    """Report whether a request announces a body at all.

    Args:
        headers: Raw ASGI header pairs.

    Returns:
        bool: True if the request carries a body framing header.
    """
    return any(name.lower() in BODY_HEADERS for name, _ in headers)


async def _refuse(send: Send, max_bytes: int) -> None:
    """Answer a request with 413 and close the connection.

    The response is built here rather than raised as an `HTTPException`, because
    middleware sits outside the framework's exception handling and has to put the
    bytes on the wire itself. Its shape matches what every other error from this
    service looks like, so a client parses one thing.

    `Connection: close` is explicit: nothing is going to read the remainder of
    the body, so the connection cannot be reused, and a client still sending must
    not be left waiting on a server that has stopped listening.

    Args:
        send: The ASGI send callable for this request.
        max_bytes: The ceiling that was passed, quoted back to the caller.
    """
    body = json.dumps({"detail": f"Request exceeds the maximum accepted size of {max_bytes} bytes"}).encode()
    await send(
        {
            "type": "http.response.start",
            "status": 413,
            "headers": [
                (b"content-type", b"application/json"),
                (b"content-length", str(len(body)).encode()),
                (b"connection", b"close"),
            ],
        }
    )
    await send({"type": "http.response.body", "body": body, "more_body": False})


class RequestSizeLimitMiddleware:
    """Bound the number of body bytes a single request may deliver.

    Attributes:
        app: The ASGI application being wrapped.
    """

    def __init__(self, app: ASGIApp, max_bytes: int | None = None) -> None:
        """Wrap an application in a request-size ceiling.

        Args:
            app: The ASGI application to wrap.
            max_bytes: A fixed ceiling in bytes. Leave unset in the service so
                the ceiling is derived from the configured upload ceiling on
                every request; pass a value to bound one application
                independently of that setting.
        """
        self.app = app
        self._max_bytes = max_bytes

    def _limit(self) -> int:
        """Return the ceiling to apply to the request being handled.

        Returns:
            int: The fixed ceiling if one was given, otherwise the configured
                one.
        """
        return configured_max_request_bytes() if self._max_bytes is None else self._max_bytes

    async def __call__(self, scope: Scope, receive: Receive, send: Send) -> None:
        """Handle one ASGI event stream.

        Args:
            scope: The connection scope.
            receive: The inbound message stream.
            send: The outbound message stream.
        """
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return

        headers = list(scope.get("headers") or [])
        if not _has_body(headers):
            await self.app(scope, receive, send)
            return

        max_bytes = self._limit()

        declared = _declared_length(headers)
        if declared is not None and declared > max_bytes:
            logger.warning(f"Refused a request declaring {declared} bytes; the ceiling is {max_bytes}")
            await _refuse(send, max_bytes)
            return

        # `refused` records that the ceiling was passed; `answered` records that
        # this middleware, rather than the application, put the response on the
        # wire. They differ in the one case where the application has already
        # started responding while still reading: there the refusal cannot be
        # sent, and the application's own response is left to finish.
        state = {"refused": False, "answered": False, "started": False}
        counted = 0

        async def limited_receive() -> Message:
            """Pass the body through, stopping it once the ceiling is passed.

            Returns:
                Message: The next inbound message, or `http.disconnect` once the
                    request has been refused.
            """
            nonlocal counted

            if state["refused"]:
                return {"type": "http.disconnect"}

            message = await receive()
            if message["type"] != "http.request":
                return message

            counted += len(message.get("body", b""))
            if counted <= max_bytes:
                return message

            state["refused"] = True
            logger.warning(f"Refused a request whose body passed the ceiling of {max_bytes} bytes")
            if not state["started"]:
                await _refuse(send, max_bytes)
                state["answered"] = True
            # The over-large chunk is not handed on, and neither is anything
            # after it. Reporting a disconnect is how ASGI says "there is no
            # more of this request", which unwinds the application's body read
            # instead of leaving it waiting for bytes that will never come.
            return {"type": "http.disconnect"}

        async def guarded_send(message: Message) -> None:
            """Forward the application's response unless the request was refused.

            Args:
                message: The outbound ASGI message.
            """
            if state["answered"]:
                # The connection already has its response; a second one would be
                # a protocol error.
                return
            if message["type"] == "http.response.start":
                state["started"] = True
            await send(message)

        try:
            await self.app(scope, limited_receive, guarded_send)
        except Exception:
            # A refused request unwinds the application by design -- reading a
            # body that stops mid-way raises, and that exception is the expected
            # end of this request, not a failure to report. Anything raised
            # while the ceiling still held is the application's own and is left
            # to the error handling above this middleware.
            if not state["answered"]:
                raise
            logger.info("Application unwound after the request was refused")
