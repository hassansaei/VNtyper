import logging
import subprocess
import time
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

    if resp.status_code != 200:
        raise RuntimeError(f"Failed to submit job. Status: {resp.status_code}, Detail: {resp.text}")
    return resp.json()


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
        resp = requests.get(status_url, timeout=30)
        if resp.status_code != 200:
            raise RuntimeError(f"Failed to get job status. Status: {resp.status_code}, Detail: {resp.text}")
        data = resp.json()
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
    resp = requests.get(dl_url, timeout=60)
    if resp.status_code != 200:
        raise RuntimeError(f"Failed to download results. Status: {resp.status_code}, Detail: {resp.text}")
    zip_path = output_dir / f"{job_id}.zip"
    with open(zip_path, "wb") as f:
        f.write(resp.content)
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

    # Use dynamic region resolution with fallback to legacy format
    region = get_region_string_with_fallback(
        bam_file=bam, reference_assembly=reference_assembly, region_type="bam_region", config=config
    )

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    job_id_file = output_path / "job_id.txt"

    if resume and job_id_file.exists():
        with open(job_id_file) as f:
            job_id = f.read().strip()
        status = poll_job_status(api_url, job_id, poll_interval=poll_interval, timeout=poll_timeout)
        _finish_job(api_url, job_id, status, output_path)
        return

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
