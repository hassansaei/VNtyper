"""Serving a file that exists only for the response carrying it.

`/cohort-download/` assembles every member's results into one temporary ZIP and
streams it back. That archive is scratch space owned by a single response: it is
not addressable by anything else, nothing will look for it again, and it holds
one cohort's results, so leaving it on a shared volume is both a leak and a
disclosure.

Starlette's usual answer is `background=BackgroundTask(...)`, but the background
task runs only after `FileResponse.__call__` has awaited every body send, and
there is no `finally` around them: a client that disconnects part-way through
the download raises out of `send` and the cleanup never happens. The one case
where the response fails is the one where nothing is reclaimed.

`ScratchFileResponse` owns the deletion instead of delegating it, in a `finally`
around the whole send, so the file goes whether the body was delivered, refused
or abandoned.
"""

import logging
from pathlib import Path

from starlette.responses import FileResponse

logger = logging.getLogger(__name__)


class ScratchFileResponse(FileResponse):
    """A `FileResponse` that removes its file once the response is over.

    Takes the same arguments as `FileResponse`. The file at `path` must belong
    to this response alone -- it is deleted unconditionally when the response
    ends, delivered or not.
    """

    async def __call__(self, scope, receive, send) -> None:
        """Send the response, then remove the file it was reading.

        Args:
            scope: ASGI connection scope.
            receive: ASGI receive callable.
            send: ASGI send callable.
        """
        try:
            await super().__call__(scope, receive, send)
        finally:
            # `missing_ok` because a repeated or racing cleanup must not raise
            # here: this runs while an exception may already be propagating, and
            # replacing a dropped connection with a FileNotFoundError would hide
            # the real failure.
            Path(self.path).unlink(missing_ok=True)
            logger.debug(f"Removed the scratch file served as {self.filename or self.path}")
