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


@contextmanager
def job_workspace(input_root: str, output_root: str, job_id: str) -> Iterator[tuple[str, str]]:
    """Create a job's input and output directories for the duration of the block.

    Both roots are created if they do not exist yet; those are shared and are
    never removed. Only the two per-job directories belong to this scope, and
    they are removed together if the block does not complete.

    Args:
        input_root: Directory holding every job's input directory.
        output_root: Directory holding every job's output directory.
        job_id: Identifier naming this job's two directories.

    Yields:
        tuple[str, str]: The job's input directory and its output directory.
    """
    os.makedirs(input_root, exist_ok=True)
    os.makedirs(output_root, exist_ok=True)

    job_input_dir = os.path.join(input_root, job_id)
    job_output_dir = os.path.join(output_root, job_id)
    os.makedirs(job_input_dir, exist_ok=True)
    os.makedirs(job_output_dir, exist_ok=True)

    try:
        yield job_input_dir, job_output_dir
    except BaseException:
        # `BaseException` rather than `Exception`: a submission cancelled while
        # it is being served unwinds through the same path, and the directories
        # it created are no more wanted than a refused submission's.
        shutil.rmtree(job_input_dir, ignore_errors=True)
        shutil.rmtree(job_output_dir, ignore_errors=True)
        logger.info(f"Reclaimed the job directories for {job_id}; the submission was not accepted")
        raise
