# docker/app/tasks.py

from .celery_app import celery_app
import subprocess
import os
import shutil
import logging
from datetime import datetime, timedelta
import redis
import hashlib
from typing import Optional, List  # Added List for typing in new task

from celery.utils.log import get_task_logger

from .config import settings
from .utils import send_email

logger = get_task_logger(__name__)

# Environment variables for Redis configuration
REDIS_HOST = os.getenv("REDIS_HOST", "redis")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))

# Retrieve Redis password from environment variables
REDIS_PASSWORD = os.getenv("REDIS_PASSWORD", "qE3!#zjraRG*`X2g4%<x&J")

# Redis DBs
REDIS_DB = int(os.getenv("REDIS_DB", 1))  # Job mappings
COHORT_REDIS_DB = int(os.getenv("COHORT_REDIS_DB", 3))  # Cohort data
USAGE_REDIS_DB = settings.USAGE_REDIS_DB  # Usage statistics

# Initialize Redis clients with authentication
redis_client = redis.Redis(
    host=REDIS_HOST,
    port=REDIS_PORT,
    db=REDIS_DB,
    password=REDIS_PASSWORD,
    decode_responses=True
)
redis_cohort_client = redis.Redis(
    host=REDIS_HOST,
    port=REDIS_PORT,
    db=COHORT_REDIS_DB,
    password=REDIS_PASSWORD,
    decode_responses=True
)
redis_usage_client = redis.Redis(
    host=REDIS_HOST,
    port=REDIS_PORT,
    db=USAGE_REDIS_DB,
    password=REDIS_PASSWORD,
    decode_responses=True
)


@celery_app.task(bind=True, max_retries=3, default_retry_delay=60)
def send_email_task(self, to_email: str, subject: str, content: str):
    """
    Celery task to send an email via SMTP.
    Retries up to 3 times in case of failure.
    """
    try:
        send_email(to_email=to_email, subject=subject, content=content)
        logger.info(f"Email sent to {to_email} with subject '{subject}'")
    except Exception as e:
        logger.error(f"Failed to send email to {to_email}: {e}")
        raise self.retry(exc=e)


@celery_app.task(bind=True)
def run_vntyper_job(
    self,
    bam_path: str,
    output_dir: str,
    thread: int,
    reference_assembly: str,
    fast_mode: bool,
    keep_intermediates: bool,
    archive_results: bool,
    email: Optional[str] = None,
    cohort_key: Optional[str] = None,
    client_ip: Optional[str] = None,
    user_agent: Optional[str] = None,
    advntr_mode: bool = False,
):
    """
    Celery task to run VNtyper pipeline with parameters.
    Sends an email upon job completion or failure if email is provided.
    """
    task_id = self.request.id
    job_id = os.path.basename(output_dir)
    try:
        logger.info(f"Starting VNtyper job for BAM file: {bam_path}")

        # Generate a unique hash for the user
        user_data = f"{client_ip}-{user_agent}"
        user_hash = hashlib.sha256(user_data.encode("utf-8")).hexdigest()

        # Store initial usage data
        usage_data = {
            "user_hash": user_hash,
            "timestamp": datetime.utcnow().isoformat(),
            "job_id": job_id,
            "status": "started",
        }
        redis_usage_client.hset(f"usage:{job_id}", mapping=usage_data)
        redis_usage_client.expire(
            f"usage:{job_id}", settings.USAGE_DATA_RETENTION_SECONDS
        )

        # Ensure the BAM index (.bai) exists
        bai_path = f"{bam_path}.bai"
        if not os.path.exists(bai_path):
            logger.info(f"BAI index not found for {bam_path}. Generating index.")
            try:
                subprocess.run(["samtools", "index", bam_path], check=True)
                logger.info(f"Successfully generated BAI index at {bai_path}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Error generating BAI index: {e}")
                # Update usage data on failure
                redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
                raise

        # Build the base command for VNtyper
        command = [
            "conda",
            "run",
            "-n",
            "vntyper",
            "vntyper",
            "pipeline",
            "--bam",
            bam_path,
            "-o",
            output_dir,
            "--thread",
            str(thread),
            "--reference-assembly",
            reference_assembly,
        ]

        if fast_mode:
            command.append("--fast-mode")
        if keep_intermediates:
            command.append("--keep-intermediates")
        if archive_results:
            command.append("--archive-results")
        if advntr_mode:
            command.extend(["--extra-modules", "advntr", "--advntr-max-coverage", "300"])

        # Run the VNtyper pipeline
        try:
            subprocess.run(command, check=True)
            logger.info(f"VNtyper job completed for {bam_path}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running VNtyper: {e}")
            # Update usage data on failure
            redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
            # Send failure email if provided
            if email:
                subject = "VNtyper Job Failed"
                if cohort_key:
                    cohort_id = cohort_key.split(":", 1)[1]
                    content = f"""
                    <p>Your VNtyper job with Job ID <strong>{job_id}</strong> has failed.</p>
                    <p>Cohort ID: <strong>{cohort_id}</strong></p>
                    <p>Error Details:</p>
                    <pre>{str(e)}</pre>
                    """
                else:
                    content = f"""
                    <p>Your VNtyper job with Job ID <strong>{job_id}</strong> has failed.</p>
                    <p>Error Details:</p>
                    <pre>{str(e)}</pre>
                    """
                send_email_task.delay(to_email=email, subject=subject, content=content)
            raise

        # Optionally, archive results
        if archive_results:
            try:
                shutil.make_archive(output_dir, "zip", output_dir)
                shutil.rmtree(output_dir)
                logger.info(
                    f"Archived results to {output_dir}.zip and removed original directory"
                )
            except Exception as e:
                logger.error(f"Error archiving results: {e}")
                # Update usage data on failure
                redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
                raise

        # Update usage data on success
        redis_usage_client.hset(f"usage:{job_id}", "status", "completed")

        # **Removed**: Cohort assignment from the task
        # if cohort_key:
        #     redis_cohort_client.sadd(f"{cohort_key}:jobs", job_id)

        # **Cohort assignment is now handled only upon successful completion**
        if cohort_key:
            redis_cohort_client.sadd(f"{cohort_key}:jobs", job_id)
            logger.info(f"Assigned job {job_id} to cohort {cohort_key}")

        # Construct the download URL
        download_url = f"{settings.API_BASE_URL}/api/download/{job_id}/"

        # Send success email if provided
        if email:
            subject = "VNtyper Job Completed Successfully"
            if cohort_key:
                cohort_id = cohort_key.split(":", 1)[1]
                content = f"""
                <p>Your VNtyper job has been completed successfully.</p>
                <p>Job ID: <strong>{job_id}</strong></p>
                <p>Cohort ID: <strong>{cohort_id}</strong></p>
                <p>You can download your results <a href="{download_url}">here</a>.</p>
                """
            else:
                content = f"""
                <p>Your VNtyper job has been completed successfully.</p>
                <p>Job ID: <strong>{job_id}</strong></p>
                <p>You can download your results <a href="{download_url}">here</a>.</p>
                """
            send_email_task.delay(to_email=email, subject=subject, content=content)
            logger.info(f"Triggered email sending to {email} with download link")

    except Exception as e:
        logger.error(f"Error in VNtyper job: {e}")
        # Update usage data on failure
        redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
        raise
    finally:
        # Remove the task ID from the Redis list
        redis_client.lrem("vntyper_job_queue", 0, self.request.id)
        logger.info(f"Removed task ID {self.request.id} from vntyper_job_queue")

        # Extend cohort TTL if necessary
        if cohort_key:
            ttl_seconds = settings.COHORT_RETENTION_DAYS * 86400
            redis_cohort_client.expire(cohort_key, ttl_seconds)
            redis_cohort_client.expire(f"{cohort_key}:jobs", ttl_seconds)

        # Delete input BAM and BAI files
        try:
            if os.path.exists(bam_path):
                os.remove(bam_path)
                logger.info(f"Deleted BAM file: {bam_path}")
        except Exception as e:
            logger.error(f"Error deleting BAM file {bam_path}: {e}")

        try:
            if os.path.exists(bai_path):
                os.remove(bai_path)
                logger.info(f"Deleted BAI file: {bai_path}")
        except Exception as e:
            logger.error(f"Error deleting BAI file {bai_path}: {e}")


@celery_app.task
def delete_old_results():
    """
    Celery task to delete result ZIP files older than a specified age.
    Also deletes cohorts that have expired.
    """
    max_age_days = settings.MAX_RESULT_AGE_DAYS
    output_dir = settings.DEFAULT_OUTPUT_DIR
    cutoff_time = datetime.now() - timedelta(days=max_age_days)

    logger.info(
        f"Running delete_old_results task. Deleting files older than {max_age_days} days."
    )

    for filename in os.listdir(output_dir):
        if filename.endswith(".zip"):
            file_path = os.path.join(output_dir, filename)
            if os.path.isfile(file_path):
                file_creation_time = datetime.fromtimestamp(os.path.getctime(file_path))
                if file_creation_time < cutoff_time:
                    try:
                        os.remove(file_path)
                        logger.info(f"Deleted old result file: {file_path}")
                    except Exception as e:
                        logger.error(f"Error deleting file {file_path}: {e}")

    # Delete expired cohorts
    for key in redis_cohort_client.scan_iter("cohort:*"):
        if not redis_cohort_client.exists(key):
            # Cohort has expired
            cohort_jobs_key = f"{key}:jobs"
            job_ids = redis_cohort_client.smembers(cohort_jobs_key)
            for job_id in job_ids:
                # Delete associated job results
                zip_path = os.path.join(output_dir, f"{job_id}.zip")
                if os.path.exists(zip_path):
                    try:
                        os.remove(zip_path)
                        logger.info(f"Deleted old result file: {zip_path}")
                    except Exception as e:
                        logger.error(f"Error deleting file {zip_path}: {e}")
                # Remove job ID from Redis
                redis_client.delete(job_id)
                redis_client.delete(f"celery-task-meta-{job_id}")  # Remove Celery task meta
            # Remove cohort jobs set
            redis_cohort_client.delete(cohort_jobs_key)
            logger.info(f"Deleted expired cohort data: {key}")

    logger.info("Completed delete_old_results task.")


# ----------------------------------------------------------------------
# Feature #82: Joint Cohort Analysis Task
# ----------------------------------------------------------------------
@celery_app.task(bind=True)
def run_cohort_analysis_job(
    self,
    cohort_id: str,
    zip_paths: List[str],
    output_dir: str,
    user_ip: Optional[str] = None,
    user_agent: Optional[str] = None,
):
    """
    Celery task to run a joint cohort analysis using 'vntyper cohort'.
    zip_paths: list of full paths to individual .zip result files for cohort samples
    output_dir: directory to store the joint analysis output and final zip.

    This task:
      1) Creates a file listing all .zip paths.
      2) Runs 'vntyper cohort --input-file <list> -o <output_dir>'.
      3) Zips the combined results and updates usage data.
    """
    task_id = self.request.id
    job_id = os.path.basename(output_dir)
    logger.info(f"Starting joint cohort analysis for Cohort ID: {cohort_id}")

    # Generate a unique hash for the user
    user_data = f"{user_ip}-{user_agent}"
    user_hash = hashlib.sha256(user_data.encode("utf-8")).hexdigest()

    # Store initial usage data
    usage_data = {
        "user_hash": user_hash,
        "timestamp": datetime.utcnow().isoformat(),
        "job_id": job_id,
        "status": "started",
        "analysis_type": "cohort_analysis",
        "cohort_id": cohort_id,
    }
    redis_usage_client.hset(f"usage:{job_id}", mapping=usage_data)
    redis_usage_client.expire(f"usage:{job_id}", settings.USAGE_DATA_RETENTION_SECONDS)

    try:
        # 1) Create directory, input file listing all .zip files
        os.makedirs(output_dir, exist_ok=True)
        input_file = os.path.join(output_dir, "cohort_input.txt")
        with open(input_file, "w") as f:
            for zpath in zip_paths:
                f.write(f"{zpath}\n")

        # 2) Run the "vntyper cohort" command
        command = [
            "conda",
            "run",
            "-n",
            "vntyper",
            "vntyper",
            "cohort",
            "--input-file",
            input_file,
            "-o",
            output_dir,
        ]
        logger.info(f"Running command: {' '.join(command)}")
        subprocess.run(command, check=True)
        logger.info("Joint cohort analysis completed.")

        # 3) Zip the results
        try:
            shutil.make_archive(output_dir, "zip", output_dir)
            logger.info(f"Zipped results to {output_dir}.zip")
        except Exception as e:
            logger.error(f"Error zipping results for cohort analysis: {e}")
            redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
            raise

        # Update usage data on success
        redis_usage_client.hset(f"usage:{job_id}", "status", "completed")

    except Exception as e:
        logger.error(f"Error in joint cohort analysis for {cohort_id}: {e}")
        redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
        raise
    finally:
        # Remove the task ID from the Redis list
        redis_client.lrem("vntyper_job_queue", 0, task_id)
        logger.info(f"Removed cohort analysis task ID {task_id} from queue")

        # Extend cohort TTL if necessary
        if cohort_id:
            cohort_key = f"cohort:{cohort_id}"
            ttl_seconds = settings.COHORT_RETENTION_DAYS * 86400
            redis_cohort_client.expire(cohort_key, ttl_seconds)
            redis_cohort_client.expire(f"{cohort_key}:jobs", ttl_seconds)

        # Delete input ZIP listing file
        try:
            if os.path.exists(input_file):
                os.remove(input_file)
                logger.info(f"Deleted cohort input file: {input_file}")
        except Exception as e:
            logger.error(f"Error deleting cohort input file {input_file}: {e}")

        # Delete individual .zip files if archive_results was used
        for zpath in zip_paths:
            try:
                if os.path.exists(zpath):
                    os.remove(zpath)
                    logger.info(f"Deleted individual result file: {zpath}")
            except Exception as e:
                logger.error(f"Error deleting file {zpath}: {e}")

        # Optionally, delete the output directory if it's empty
        try:
            if os.path.exists(output_dir) and not os.listdir(output_dir):
                os.rmdir(output_dir)
                logger.info(f"Deleted empty output directory: {output_dir}")
        except Exception as e:
            logger.error(f"Error deleting directory {output_dir}: {e}")
