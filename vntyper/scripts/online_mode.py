import logging
import os
import subprocess
import tempfile
import time
from collections.abc import Callable
from pathlib import Path

import requests

from vntyper.scripts.region_utils import get_region_string_with_fallback

logger = logging.getLogger(__name__)

#: The only status that means the remote job produced results.
STATUS_COMPLETED = "completed"

#: Statuses the web service reports while a job is still on its way to an answer.
#: ``docker/app/main.py`` maps Celery's ``PENDING`` and ``STARTED`` explicitly and
#: passes anything else through ``status.lower()``, so a terminal Celery state such as
#: ``REVOKED`` arrives as an unrecognised string. Treating *only* these as
#: non-terminal is what stops the poller looping forever on one.
IN_PROGRESS_STATUSES = frozenset({"pending", "started", "retry"})

#: Seconds between polls, when the configuration does not say.
DEFAULT_POLL_INTERVAL_SECONDS = 10

#: Total seconds to keep polling before giving up. A VNtyper job including adVNTR can
#: take well over an hour on a small runner, so this is generous - but it is finite,
#: which the previous ``while True`` was not.
DEFAULT_POLL_TIMEOUT_SECONDS = 4 * 60 * 60

#: Short, bounded retries for idempotent API reads. The polling timeout owns
#: long-horizon patience; this policy only rides out a connection or gateway hiccup.
TRANSIENT_GET_ATTEMPTS = 3
TRANSIENT_GET_BACKOFF_SECONDS = 2.0

#: Downloaded archives are streamed in bounded pieces instead of being buffered in RAM.
DOWNLOAD_CHUNK_SIZE = 1 << 20


def _close_response(response: requests.Response, url: str) -> None:
    """Release an HTTP response without masking the request's primary outcome.

    Args:
        response: Response whose connection should be released.
        url: Request URL used in diagnostics.
    """
    try:
        response.close()
    except Exception as error:  # noqa: BLE001 - cleanup cannot mask the primary request or download outcome
        logger.warning(f"Failed to close response from {url}: {error}")


def _send_get_with_retries(
    url: str,
    send: Callable[[], requests.Response],
    *,
    attempts: int = TRANSIENT_GET_ATTEMPTS,
    backoff_seconds: float = TRANSIENT_GET_BACKOFF_SECONDS,
) -> requests.Response:
    """Send an idempotent GET, retrying transient failures with backoff.

    Connection errors, timeouts and 5xx responses are retried. A 4xx response is
    returned to the caller as a fact about the request. Discarded 5xx responses are
    closed here; the caller owns and must close the returned response.

    Args:
        url: URL used in diagnostics.
        send: Zero-argument callable that performs the GET.
        attempts: Maximum number of requests.
        backoff_seconds: Initial sleep before a retry; subsequent sleeps double.

    Returns:
        The first non-5xx response, or the final 5xx response.

    Raises:
        RuntimeError: If every attempt failed before returning a response.
    """
    for attempt in range(1, attempts + 1):
        final_attempt = attempt == attempts
        try:
            response = send()
        except (requests.ConnectionError, requests.Timeout) as error:
            if final_attempt:
                message = f"{url} failed after {attempts} attempts: {error}"
                logger.error(message)
                raise RuntimeError(message) from error
            logger.warning(f"{url} failed ({error}); attempt {attempt}/{attempts}")
        else:
            if response.status_code < 500 or final_attempt:
                return response
            logger.warning(f"{url} returned {response.status_code}; attempt {attempt}/{attempts}")
            _close_response(response, url)
        time.sleep(backoff_seconds * (2 ** (attempt - 1)))

    message = f"{url} failed after {attempts} attempts"
    logger.error(message)
    raise RuntimeError(message)


def subset_bam(bam_path, region, output_bam):
    """
    Subset the BAM file to a specified region using samtools and index it.

    Args:
        bam_path (str): Path to the original BAM file.
        region (str): Genomic region (e.g. chr1:1000-2000).
        output_bam (str): Path to the subsetted BAM file.

    Raises:
        RuntimeError: If subset or indexing fails.
    """
    cmd = ["samtools", "view", "-P", "-b", bam_path, region, "-o", output_bam]
    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Failed to subset BAM: {result.stderr}")

    # Index the subset bam
    cmd_index = ["samtools", "index", output_bam]
    logger.info(f"Running: {' '.join(cmd_index)}")
    result_idx = subprocess.run(cmd_index, capture_output=True, text=True, check=False)
    if result_idx.returncode != 0:
        raise RuntimeError(f"Failed to index subset BAM: {result_idx.stderr}")


def submit_job(
    api_url,
    subset_bam,
    subset_bai,
    reference_assembly,
    threads,
    email=None,
    cohort_id=None,
    passphrase=None,
):
    """
    Submit a job to the online API.

    Args:
        api_url (str): The base URL of the API.
        subset_bam (str): Path to subset BAM.
        subset_bai (str): Path to subset BAI.
        reference_assembly (str): Reference assembly used.
        threads (int): Number of threads.
        email (str or None): Email for notification.
        cohort_id (str or None): Cohort ID.
        passphrase (str or None): Passphrase for cohort.

    Returns:
        dict: JSON response with job_id and message.

    Raises:
        RuntimeError: If submission fails.
    """
    data = {
        "thread": threads,
        "reference_assembly": reference_assembly,
        "fast_mode": "false",
        "keep_intermediates": "false",
        "archive_results": "true",
    }
    if email:
        data["email"] = email
    if cohort_id:
        data["cohort_id"] = cohort_id
    if passphrase:
        data["passphrase"] = passphrase

    submit_url = f"{api_url}/run-job/"
    logger.info(f"Submitting job to {submit_url}")

    with open(subset_bam, "rb") as bam_file:
        files = {"bam_file": ("subset.bam", bam_file, "application/octet-stream")}
        if subset_bai:
            with open(subset_bai, "rb") as bai_file:
                files["bai_file"] = ("subset.bam.bai", bai_file, "application/octet-stream")
                resp = requests.post(submit_url, files=files, data=data, timeout=60)
        else:
            resp = requests.post(submit_url, files=files, data=data, timeout=60)

    try:
        if resp.status_code != 200:
            message = f"Failed to submit job. Status: {resp.status_code}, Detail: {resp.text}"
            logger.error(message)
            raise RuntimeError(message)
        return resp.json()
    finally:
        _close_response(resp, submit_url)


def poll_job_status(
    api_url,
    job_id,
    poll_interval=DEFAULT_POLL_INTERVAL_SECONDS,
    timeout=DEFAULT_POLL_TIMEOUT_SECONDS,
):
    """
    Poll the job status until it reaches a terminal state or the timeout expires.

    Only :data:`IN_PROGRESS_STATUSES` are treated as "keep asking". Anything else is
    returned to the caller, including a status this client does not recognise: the
    previous ``status in ["completed", "failed"]`` test meant a Celery ``REVOKED``
    job - which is terminal - polled until the operator killed the process.

    The client reads ``status`` for its decision and treats ``error`` as opaque text.
    ``/job-status/`` deliberately returns a **generic** failure message rather than the
    task's exception detail, because the endpoint is unauthenticated; nothing here may
    depend on parsing it.

    Args:
        api_url (str): The base URL of the API.
        job_id (str): Job ID to poll.
        poll_interval (float): Seconds to wait between polls.
        timeout (float): Total seconds to keep polling before giving up.

    Returns:
        str: The terminal status the server reported, e.g. ``'completed'`` or
        ``'failed'``.

    Raises:
        RuntimeError: If the status request fails, or the job never reaches a
            terminal status within ``timeout``.
    """
    status_url = f"{api_url}/job-status/{job_id}/"
    # Counted rather than measured against the wall clock: the bound then holds
    # however long the server takes to answer, and it is observable from a test that
    # does not have to sleep for it.
    max_polls = max(1, int(timeout // poll_interval))
    polls = 0
    while True:
        polls += 1
        logger.info(f"Checking job status for {job_id}")
        resp = _send_get_with_retries(status_url, lambda: requests.get(status_url, timeout=30))
        try:
            if resp.status_code != 200:
                message = f"Failed to get job status. Status: {resp.status_code}, Detail: {resp.text}"
                logger.error(message)
                raise RuntimeError(message)
            data = resp.json()
        finally:
            _close_response(resp, status_url)
        status = data.get("status", "")
        detail = data.get("error")

        if status not in IN_PROGRESS_STATUSES:
            if status == STATUS_COMPLETED:
                logger.info(f"Job {job_id} completed.")
            elif detail:
                logger.error(f"Job {job_id} ended with status {status!r}: {detail}")
            else:
                logger.warning(f"Job {job_id} ended with status {status!r}.")
            return status

        if polls >= max_polls:
            msg = (
                f"Job {job_id} did not reach a terminal status within {timeout} seconds "
                f"({max_polls} polls at {poll_interval}s; last status: {status!r}). "
                "Re-run with --resume to keep waiting."
            )
            logger.error(msg)
            raise RuntimeError(msg)

        logger.info(f"Job {job_id} status: {status}, waiting...")
        time.sleep(poll_interval)


def download_results(api_url, job_id, output_dir):
    """
    Download the job results once completed.

    Args:
        api_url (str): The base URL of the API.
        job_id (str): Job ID.
        output_dir (Path): Output directory to save results.

    Raises:
        RuntimeError: If download fails.
    """
    dl_url = f"{api_url}/download/{job_id}/"
    logger.info(f"Downloading results from {dl_url}")
    resp = _send_get_with_retries(dl_url, lambda: requests.get(dl_url, timeout=60, stream=True))
    output_path = Path(output_dir)
    zip_path = output_path / f"{job_id}.zip"
    temp_path = None
    try:
        if resp.status_code != 200:
            message = f"Failed to download results. Status: {resp.status_code}, Detail: {resp.text}"
            logger.error(message)
            raise RuntimeError(message)
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=output_path,
            prefix=f".{job_id}.",
            suffix=".tmp",
            delete=False,
        ) as temp_file:
            temp_path = Path(temp_file.name)
            for chunk in resp.iter_content(chunk_size=DOWNLOAD_CHUNK_SIZE):
                if chunk:
                    temp_file.write(chunk)
        os.replace(temp_path, zip_path)
    finally:
        try:
            _close_response(resp, dl_url)
        finally:
            if temp_path is not None:
                temp_path.unlink(missing_ok=True)
    logger.info(f"Results saved to {zip_path}")


def _finish_job(api_url, job_id, status, output_path):
    """
    Download the results of a completed job, or fail loudly.

    Both the fresh-submission and ``--resume`` branches ended in their own copy of
    this decision, and both used to log ``"Job failed or status unknown."`` and
    return - so a failed remote job produced exit code 0.

    Args:
        api_url (str): The base URL of the API.
        job_id (str): The job that was polled.
        status (str): The terminal status ``poll_job_status`` returned.
        output_path (Path): Directory to save results into.

    Raises:
        RuntimeError: If ``status`` is anything but ``completed``.
    """
    if status == STATUS_COMPLETED:
        download_results(api_url, job_id, output_path)
        logger.info("Job completed successfully.")
        return

    msg = (
        f"Online job {job_id} did not complete: status {status!r}. "
        "No results were downloaded. The server reports a generic failure message by design; "
        "ask the instance operator for the job's log if you need the detail."
    )
    logger.error(msg)
    raise RuntimeError(msg)


def run_online_mode(
    config,
    bam,
    output_dir,
    reference_assembly,
    threads,
    email=None,
    cohort_id=None,
    passphrase=None,
    resume=False,
):
    """
    Run the online mode functionality:
    1. Determine region from config based on reference assembly.
    2. Subset the BAM file.
    3. Submit job to online API.
    4. Poll for completion.
    5. Download results.

    Args:
        config (dict): Loaded configuration dictionary.
        bam (str): Input BAM file path.
        output_dir (str): Directory for outputs.
        reference_assembly (str): Reference assembly ("hg19" or "hg38").
        threads (int): Number of threads.
        email (str or None): Notification email.
        cohort_id (str or None): Cohort ID.
        passphrase (str or None): Cohort passphrase.
        resume (bool): Whether to resume a previous job.

    Raises:
        ValueError: If ``resume`` is requested without a stored non-empty job id.
        RuntimeError: If the API returns no job id, or the job reaches any terminal
            status other than ``completed``. This used to be logged and swallowed, so
            ``main()`` exited 0 and a wrapping ``subprocess.run(..., check=True)`` saw
            a failed genotyping run as a success.
    """

    # Get the API base URL from config
    api_settings = config.get("api", {})
    api_url = api_settings.get("base_url", "http://vntyper.org/api")
    poll_interval = api_settings.get("poll_interval_seconds", DEFAULT_POLL_INTERVAL_SECONDS)
    poll_timeout = api_settings.get("poll_timeout_seconds", DEFAULT_POLL_TIMEOUT_SECONDS)

    output_path = Path(output_dir)
    job_id_file = output_path / "job_id.txt"

    if resume:
        job_id = job_id_file.read_text(encoding="utf-8").strip() if job_id_file.exists() else ""
        if not job_id:
            message = (
                f"--resume was given but {job_id_file} contains no job to resume. "
                "Run without --resume to submit a new job."
            )
            logger.error(message)
            raise ValueError(message)
        status = poll_job_status(api_url, job_id, poll_interval=poll_interval, timeout=poll_timeout)
        _finish_job(api_url, job_id, status, output_path)
        return

    # Use dynamic region resolution with fallback to legacy format only when a new
    # job was requested. A rejected resume must not invoke local alignment tooling.
    region = get_region_string_with_fallback(
        bam_file=bam, reference_assembly=reference_assembly, region_type="bam_region", config=config
    )

    output_path.mkdir(parents=True, exist_ok=True)

    # Fresh submission
    subset_bam_path = output_path / "subset.bam"
    subset_bam(bam, region, str(subset_bam_path))
    subset_bai_path = output_path / "subset.bam.bai"

    resp = submit_job(
        api_url=api_url,
        subset_bam=str(subset_bam_path),
        subset_bai=str(subset_bai_path),
        reference_assembly=reference_assembly,
        threads=threads,
        email=email,
        cohort_id=cohort_id,
        passphrase=passphrase,
    )
    job_id = resp.get("job_id")
    if not job_id:
        msg = f"No job_id returned from API at {api_url}; the job was not submitted."
        logger.error(msg)
        raise RuntimeError(msg)

    with open(job_id_file, "w") as f:
        f.write(job_id)

    # Poll until completed
    status = poll_job_status(api_url, job_id, poll_interval=poll_interval, timeout=poll_timeout)
    _finish_job(api_url, job_id, status, output_path)
