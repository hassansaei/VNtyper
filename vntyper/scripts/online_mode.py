import logging
import subprocess
import time
import requests
from pathlib import Path


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
    logging.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Failed to subset BAM: {result.stderr}")

    # Index the subset bam
    cmd_index = ["samtools", "index", output_bam]
    logging.info(f"Running: {' '.join(cmd_index)}")
    result_idx = subprocess.run(cmd_index, capture_output=True, text=True)
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
    files = {
        "bam_file": ("subset.bam", open(subset_bam, "rb"), "application/octet-stream")
    }
    if subset_bai:
        files["bai_file"] = (
            "subset.bam.bai",
            open(subset_bai, "rb"),
            "application/octet-stream",
        )

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
    logging.info(f"Submitting job to {submit_url}")
    try:
        resp = requests.post(submit_url, files=files, data=data, timeout=60)
    finally:
        for f in files.values():
            f[1].close()

    if resp.status_code != 200:
        raise RuntimeError(
            f"Failed to submit job. Status: {resp.status_code}, Detail: {resp.text}"
        )
    return resp.json()


def poll_job_status(api_url, job_id):
    """
    Poll the job status until completion or failure.

    Args:
        api_url (str): The base URL of the API.
        job_id (str): Job ID to poll.

    Returns:
        str: 'completed' or 'failed'.

    Raises:
        RuntimeError: If unable to get job status.
    """
    status_url = f"{api_url}/job-status/{job_id}/"
    while True:
        logging.info(f"Checking job status for {job_id}")
        resp = requests.get(status_url, timeout=30)
        if resp.status_code != 200:
            raise RuntimeError(
                f"Failed to get job status. Status: {resp.status_code}, Detail: {resp.text}"
            )
        data = resp.json()
        status = data.get("status", "")
        if status in ["completed", "failed"]:
            return status
        logging.info(f"Job {job_id} status: {status}, waiting...")
        time.sleep(10)


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
    logging.info(f"Downloading results from {dl_url}")
    resp = requests.get(dl_url, timeout=60)
    if resp.status_code != 200:
        raise RuntimeError(
            f"Failed to download results. Status: {resp.status_code}, Detail: {resp.text}"
        )
    zip_path = output_dir / f"{job_id}.zip"
    with open(zip_path, "wb") as f:
        f.write(resp.content)
    logging.info(f"Results saved to {zip_path}")


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
    """

    # Get the API base URL from config
    api_url = config.get("api", {}).get("base_url", "http://vntyper.org/api")

    # Determine region from config based on reference_assembly
    bam_processing = config.get("bam_processing", {})
    region_key = f"bam_region_{reference_assembly}"
    if region_key not in bam_processing:
        raise RuntimeError(f"No region configured for {reference_assembly}")
    region = bam_processing[region_key]

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    job_id_file = output_path / "job_id.txt"

    if resume and job_id_file.exists():
        with open(job_id_file, "r") as f:
            job_id = f.read().strip()
        status = poll_job_status(api_url, job_id)
        if status == "completed":
            download_results(api_url, job_id, output_path)
            logging.info("Job completed successfully.")
        else:
            logging.error("Job failed or status unknown.")
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
        logging.error("No job_id returned from API.")
        return

    with open(job_id_file, "w") as f:
        f.write(job_id)

    # Poll until completed
    status = poll_job_status(api_url, job_id)
    if status == "completed":
        download_results(api_url, job_id, output_path)
        logging.info("Job completed successfully.")
    else:
        logging.error("Job failed or status unknown.")
