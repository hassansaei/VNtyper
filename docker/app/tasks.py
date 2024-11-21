from .celery_app import celery_app
import subprocess
import os
import shutil
import logging
from datetime import datetime, timedelta
import redis
import hashlib
from typing import Optional

from celery.utils.log import get_task_logger

from .config import settings
from .utils import send_email

logger = get_task_logger(__name__)

# Environment variables for Redis configuration
REDIS_HOST = os.getenv("REDIS_HOST", "redis")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))

# Redis DBs
REDIS_DB = int(os.getenv("REDIS_DB", 1))  # Job mappings
COHORT_REDIS_DB = int(os.getenv("COHORT_REDIS_DB", 3))  # Cohort data
USAGE_REDIS_DB = settings.USAGE_REDIS_DB  # Usage statistics

# Initialize Redis clients
redis_client = redis.Redis(
    host=REDIS_HOST, port=REDIS_PORT, db=REDIS_DB, decode_responses=True
)
redis_cohort_client = redis.Redis(
    host=REDIS_HOST, port=REDIS_PORT, db=COHORT_REDIS_DB, decode_responses=True
)
redis_usage_client = redis.Redis(
    host=REDIS_HOST, port=REDIS_PORT, db=USAGE_REDIS_DB, decode_responses=True
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
                subprocess.run(
                    ["samtools", "index", bam_path],
                    check=True,
                )
                logger.info(f"Successfully generated BAI index at {bai_path}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Error generating BAI index: {e}")
                # Update usage data on failure
                redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
                raise

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

        try:
            # Run the VNtyper pipeline
            subprocess.run(command, check=True)
            logger.info(f"VNtyper job completed for {bam_path}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running VNtyper: {e}")
            # Update usage data on failure
            redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
            # Send failure email if provided
            if email:
                subject = "VNtyper Job Failed"
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

        # Send success email if provided
        if email:
            subject = "VNtyper Job Completed Successfully"
            download_url = f"{settings.API_BASE_URL}/api/download/{job_id}/"
            content = f"""
            <p>Your VNtyper job has been completed successfully.</p>
            <p>Job ID: <strong>{job_id}</strong></p>
            <p>You can download your results <a href="{download_url}">here</a>.</p>
            """
            send_email_task.delay(to_email=email, subject=subject, content=content)
            logger.info(f"Triggered email sending to {email}")

    except Exception as e:
        logger.error(f"Error in VNtyper job: {e}")
        # Update usage data on failure
        redis_usage_client.hset(f"usage:{job_id}", "status", "failed")
        raise
    finally:
        # Remove the task ID from the Redis list
        redis_client.lrem("vntyper_job_queue", 0, task_id)
        logger.info(f"Removed task ID {task_id} from vntyper_job_queue")

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
