"""The per-job directories a submission writes into, and how long they live.

`/run-job/` gives every submission an identifier and two directories named after
it: one holding the uploaded alignment and its index, one the worker writes its
results into. Both sit on a volume every job shares.

A submission can still be refused after those directories exist -- the request
carries an email address and cohort credentials that are only meaningful once
parsed, and enqueueing the job can fail on its own. A refusal that has already
written to the shared volume has to take what it wrote with it: the caller is
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
from collections.abc import Iterator
from contextlib import contextmanager

logger = logging.getLogger(__name__)


class JobWorkspace:
    """A job's two directories, and whether they have been handed to a worker.

    Iterating yields `(input_dir, output_dir)`, so `as (job_input, job_output)`
    keeps working; the attributes and `commit()` are what the commit point adds.

    Attributes:
        input_dir: The job's input directory.
        output_dir: The job's output directory.
    """

    def __init__(self, input_dir: str, output_dir: str) -> None:
        """Record the two directories, uncommitted.

        Args:
            input_dir: The job's input directory.
            output_dir: The job's output directory.
        """
        self.input_dir = input_dir
        self.output_dir = output_dir
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

    def __iter__(self) -> Iterator[str]:
        """Yield the input and output directories, in that order.

        Returns:
            Iterator[str]: The two directories, so the object unpacks as a pair.
        """
        return iter((self.input_dir, self.output_dir))


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
    """
    os.makedirs(input_root, exist_ok=True)
    os.makedirs(output_root, exist_ok=True)

    job_input_dir = os.path.join(input_root, job_id)
    job_output_dir = os.path.join(output_root, job_id)
    os.makedirs(job_input_dir, exist_ok=True)
    os.makedirs(job_output_dir, exist_ok=True)

    workspace = JobWorkspace(job_input_dir, job_output_dir)
    try:
        yield workspace
    except BaseException:
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
        shutil.rmtree(job_input_dir, ignore_errors=True)
        shutil.rmtree(job_output_dir, ignore_errors=True)
        logger.info(f"Reclaimed the job directories for {job_id}; the submission was not accepted")
        raise
