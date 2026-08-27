"""The per-job protected handoff spool and service-private result store.

`/run-job/` gives every submission an identifier and two directories named after
it: one in the protected API/worker handoff spool holding the accepted alignment
and index, and one in the service-private result store. They are separate
trust surfaces even though this scope owns their admission lifetime together.

A submission can still be refused after those directories exist -- the request
carries an email address and cohort credentials that are only meaningful once
parsed, and enqueueing the job can fail on its own. A refusal that has already
written to service-owned storage has to take what it wrote with it: the caller is
never given an identifier, so nothing afterwards can address the leftovers, and
the worker's own cleanup only ever runs for a job that was actually enqueued.

`job_workspace` states that as a scope. Directories created on the way in are
removed on the way out if, and only if, the body did not complete. Nothing else
in the service creates them, so the guarantee holds wherever a submission is
refused rather than only where a particular refusal was anticipated.

That scope ends where ownership changes hands. Once the broker has accepted the
task, a worker holds the input path and is reading it, and the request's later
failures are no longer evidence that nothing is running -- deleting the
directory then takes a live job's input away from it. `commit()` marks that
moment explicitly, and after it the reclaim does not run.

The roots are passed in rather than read from settings, so the caller decides
which volume these directories live on and the tests can point them at a scratch
directory.
"""

import logging
import os
import shutil
import stat
from collections.abc import Iterator
from contextlib import contextmanager, suppress

logger = logging.getLogger(__name__)


class JobWorkspace:
    """A job's two directories, and whether they have been handed to a worker.

    Iterating yields `(input_dir, output_dir)`, so `as (job_input, job_output)`
    keeps working; the attributes and `commit()` are what the commit point adds.

    Attributes:
        input_dir: The job's input directory.
        output_dir: The job's output directory.
    """

    def __init__(self, input_dir: str, output_dir: str, input_descriptor: int, output_descriptor: int) -> None:
        """Record the two directories, uncommitted.

        Args:
            input_dir: The job's input directory.
            output_dir: The job's output directory.
        """
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.input_descriptor = input_descriptor
        self.output_descriptor = output_descriptor
        input_metadata = os.fstat(input_descriptor)
        output_metadata = os.fstat(output_descriptor)
        self._input_identity = (input_metadata.st_dev, input_metadata.st_ino)
        self._output_identity = (output_metadata.st_dev, output_metadata.st_ino)
        self._committed = False

    @property
    def committed(self) -> bool:
        """bool: Whether ownership has passed to a worker."""
        return self._committed

    def commit(self) -> None:
        """Hand the directories to the worker; the reclaim will no longer run.

        Call this as soon as the task has been accepted, and never before: every
        line between the upload and this call is still covered by the reclaim.
        """
        self._committed = True

    def handoff(
        self,
        alignment_identity: tuple[int, int],
        index_identity: tuple[int, int] | None,
        alignment_sha256: str,
        index_sha256: str | None,
    ) -> dict[str, list[int] | str | None]:
        """Return the JSON-safe inode identities a worker must reopen."""
        _require_current_directory(self.input_dir, self._input_identity, "input")
        _require_current_directory(self.output_dir, self._output_identity, "output")
        return {
            "input_dir": list(self._input_identity),
            "output_dir": list(self._output_identity),
            "alignment": list(alignment_identity),
            "alignment_sha256": alignment_sha256,
            "index": None if index_identity is None else list(index_identity),
            "index_sha256": index_sha256,
        }

    def __iter__(self) -> Iterator[str]:
        """Yield the input and output directories, in that order.

        Returns:
            Iterator[str]: The two directories, so the object unpacks as a pair.
        """
        return iter((self.input_dir, self.output_dir))


def _reclaim(path: str) -> None:
    """Remove a per-job directory without deleting through a symlink.

    ``shutil.rmtree(..., ignore_errors=True)`` leaves a symlink to a directory
    in place. A link swapped in after creation belongs to the refused request,
    but its target does not: unlink the entry itself and never traverse it.

    Args:
        path: The per-job directory, or whatever now occupies its name.
    """
    with suppress(OSError):
        if os.path.islink(path):
            os.unlink(path)
        else:
            shutil.rmtree(path, ignore_errors=True)


def _same_identity(metadata: os.stat_result, expected: tuple[int, int]) -> bool:
    return (metadata.st_dev, metadata.st_ino) == expected


def _require_current_directory(path: str, expected: tuple[int, int], label: str) -> None:
    try:
        metadata = os.stat(path, follow_symlinks=False)
    except OSError as error:
        msg = f"Refusing job: {label} directory identity changed before handoff: {error}"
        logger.error(msg)
        raise RuntimeError(msg) from error
    if not stat.S_ISDIR(metadata.st_mode) or not _same_identity(metadata, expected):
        msg = f"Refusing job: {label} directory identity changed before handoff"
        logger.error(msg)
        raise RuntimeError(msg)


def _clear_directory(descriptor: int) -> None:
    for name in os.listdir(descriptor):
        metadata = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
        if not stat.S_ISDIR(metadata.st_mode):
            os.unlink(name, dir_fd=descriptor)
            continue
        child_descriptor = os.open(
            name,
            os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0),
            dir_fd=descriptor,
        )
        try:
            opened = os.fstat(child_descriptor)
            if not _same_identity(opened, (metadata.st_dev, metadata.st_ino)):
                raise OSError(f"job workspace child {name} changed while being opened")
            _clear_directory(child_descriptor)
            current = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
            if not _same_identity(current, (opened.st_dev, opened.st_ino)):
                raise OSError(f"job workspace child {name} changed during cleanup")
            os.rmdir(name, dir_fd=descriptor)
        except BaseException:
            try:
                os.close(child_descriptor)
            except OSError as close_error:
                logger.warning(f"Closing refused job workspace child also failed: {close_error}")
            raise
        else:
            os.close(child_descriptor)


def _reclaim_bound(path: str, descriptor: int, expected: tuple[int, int]) -> None:
    """Clear one captured directory and remove only its still-current name."""
    _clear_directory(descriptor)
    try:
        current = os.stat(path, follow_symlinks=False)
    except OSError:
        return
    if stat.S_ISDIR(current.st_mode) and _same_identity(current, expected):
        os.rmdir(path)


@contextmanager
def job_workspace(input_root: str, output_root: str, job_id: str) -> Iterator[JobWorkspace]:
    """Create a job's input and output directories for the duration of the block.

    Both roots are created if they do not exist yet; those are shared and are
    never removed. Only the two per-job directories belong to this scope, and
    they are removed together if the block does not complete -- unless
    :meth:`JobWorkspace.commit` has been called, which ends the scope early
    because a worker now owns them.

    Args:
        input_root: Directory holding every job's input directory.
        output_root: Directory holding every job's output directory.
        job_id: Identifier naming this job's two directories.

    Yields:
        JobWorkspace: The job's directories, unpackable as
        `(input_dir, output_dir)`.

    Raises:
        RuntimeError: If either per-job entry cannot be freshly created. A
            server-generated identifier finding its name already occupied is
            refused rather than adopting the existing directory or symlink.
    """
    os.makedirs(input_root, exist_ok=True)
    os.makedirs(output_root, exist_ok=True)

    job_input_dir = os.path.join(input_root, job_id)
    job_output_dir = os.path.join(output_root, job_id)
    try:
        os.mkdir(job_input_dir)
    except OSError as error:
        msg = f"Refusing job {job_id}: its input directory could not be freshly created: {error}"
        logger.error(msg)
        raise RuntimeError(msg) from error
    try:
        os.mkdir(job_output_dir)
    except OSError as error:
        _reclaim(job_input_dir)
        msg = f"Refusing job {job_id}: its output directory could not be freshly created: {error}"
        logger.error(msg)
        raise RuntimeError(msg) from error

    input_descriptor = -1
    output_descriptor = -1
    try:
        input_descriptor = os.open(job_input_dir, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC)
        output_descriptor = os.open(job_output_dir, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC)
        if not stat.S_ISDIR(os.fstat(input_descriptor).st_mode) or not stat.S_ISDIR(
            os.fstat(output_descriptor).st_mode
        ):
            raise OSError("created job workspace entry is not a directory")
    except OSError as error:
        close_errors: list[OSError] = []
        for descriptor in (input_descriptor, output_descriptor):
            if descriptor >= 0:
                try:
                    os.close(descriptor)
                except OSError as close_error:
                    close_errors.append(close_error)
        _reclaim(job_input_dir)
        _reclaim(job_output_dir)
        for secondary_error in close_errors:
            logger.warning(f"Closing a refused job workspace also failed: {secondary_error}")
        msg = f"Refusing job {job_id}: its directories could not be bound: {error}"
        logger.error(msg)
        raise RuntimeError(msg) from error

    try:
        workspace = JobWorkspace(job_input_dir, job_output_dir, input_descriptor, output_descriptor)
    except OSError as error:
        close_errors = []
        for descriptor in (input_descriptor, output_descriptor):
            try:
                os.close(descriptor)
            except OSError as close_error:
                close_errors.append(close_error)
        _reclaim(job_input_dir)
        _reclaim(job_output_dir)
        for secondary_error in close_errors:
            logger.warning(f"Closing a refused job workspace also failed: {secondary_error}")
        msg = f"Refusing job {job_id}: its directories could not be bound: {error}"
        logger.error(msg)
        raise RuntimeError(msg) from error
    body_error: BaseException | None = None
    try:
        yield workspace
    except BaseException as error:
        body_error = error
        # `BaseException` rather than `Exception`: a submission cancelled while
        # it is being served unwinds through the same path, and the directories
        # it created are no more wanted than a refused submission's.
        if workspace.committed:
            # A worker is reading `job_input_dir` right now. Whatever failed here
            # failed after the handover, so it says nothing about the task.
            logger.error(
                f"Job {job_id} was accepted but its submission then failed; "
                f"its directories are kept because a worker owns them"
            )
            raise
        for path, descriptor, identity in (
            (job_input_dir, input_descriptor, workspace._input_identity),
            (job_output_dir, output_descriptor, workspace._output_identity),
        ):
            try:
                _reclaim_bound(path, descriptor, identity)
            except OSError as cleanup_error:
                logger.error(f"Error reclaiming refused job directory {path}: {cleanup_error}")
        logger.info(f"Reclaimed the job directories for {job_id}; the submission was not accepted")
        raise
    finally:
        first_close_error: OSError | None = None
        for descriptor in (output_descriptor, input_descriptor):
            try:
                os.close(descriptor)
            except OSError as close_error:
                if first_close_error is None:
                    first_close_error = close_error
        if first_close_error is not None:
            if body_error is None:
                raise first_close_error
            logger.warning(f"Closing a failed job workspace also failed: {first_close_error}")
