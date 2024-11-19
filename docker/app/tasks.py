# docker/app/tasks.py

from .celery_app import celery_app
import subprocess
import os
import shutil
import logging
from datetime import datetime, timedelta
import redis  # Import Redis

logger = logging.getLogger(__name__)

# Environment variables for Redis configuration
REDIS_HOST = os.getenv("REDIS_HOST", "redis")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))
REDIS_DB = int(os.getenv("REDIS_DB", 1))  # Use DB 1 for job mappings

# Initialize Redis client
redis_client = redis.Redis(
    host=REDIS_HOST, port=REDIS_PORT, db=REDIS_DB, decode_responses=True
)

@celery_app.task(bind=True)
def run_vntyper_job(
    self,
    bam_path: str,
    output_dir: str,
    thread: int,
    reference_assembly: str,
    fast_mode: bool,
    keep_intermediates: bool,
    archive_results: bool
):
    """
    Celery task to run VNtyper pipeline with parameters.
    Automatically generates BAI index if missing and cleans up input files.
    """
    task_id = self.request.id
    try:
        logger.info(f"Starting VNtyper job for BAM file: {bam_path}")

        # Ensure the BAM index (.bai) exists
        bai_path = f"{bam_path}.bai"
        if not os.path.exists(bai_path):
            logger.info(f"BAI index not found for {bam_path}. Generating index.")
            try:
                subprocess.run(
                    ["samtools", "index", bam_path],
                    check=True
                )
                logger.info(f"Successfully generated BAI index at {bai_path}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Error generating BAI index: {e}")
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
            reference_assembly
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
            raise

        # Optionally, archive results
        if archive_results:
            try:
                shutil.make_archive(output_dir, "zip", output_dir)
                shutil.rmtree(output_dir)
                logger.info(f"Archived results to {output_dir}.zip and removed original directory")
            except Exception as e:
                logger.error(f"Error archiving results: {e}")
                raise
    except Exception as e:
        logger.error(f"Error in VNtyper job: {e}")
        raise
    finally:
        # Remove the task ID from the Redis list
        redis_client.lrem('vntyper_job_queue', 0, task_id)
        logger.info(f"Removed task ID {task_id} from vntyper_job_queue")

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
    """
    max_age_days = int(os.getenv("MAX_RESULT_AGE_DAYS", "7"))  # Default to 7 days
    output_dir = os.getenv("DEFAULT_OUTPUT_DIR", "/opt/vntyper/output")
    cutoff_time = datetime.now() - timedelta(days=max_age_days)

    logger.info(f"Running delete_old_results task. Deleting files older than {max_age_days} days.")

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
