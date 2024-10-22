# docker/app/tasks.py

from .celery_app import celery_app
import subprocess
import os
import shutil
import logging

logger = logging.getLogger(__name__)

@celery_app.task
def run_vntyper_job(
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
    """
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
