# backend/docker/app/main.py

from fastapi import FastAPI, APIRouter, UploadFile, File, Form, HTTPException
from fastapi.responses import FileResponse
from uuid import uuid4
import os
from .tasks import run_vntyper_job
import shutil
import logging

from .config import settings

from celery.result import AsyncResult
import redis

app = FastAPI(
    title=settings.PROJECT_NAME,
    version=settings.PROJECT_VERSION
)

logger = logging.getLogger(__name__)

# Environment variables for default directories
DEFAULT_INPUT_DIR = os.getenv("DEFAULT_INPUT_DIR", "/opt/vntyper/input")
DEFAULT_OUTPUT_DIR = os.getenv("DEFAULT_OUTPUT_DIR", "/opt/vntyper/output")

# Redis configuration for job_id to task_id mapping
REDIS_HOST = os.getenv("REDIS_HOST", "redis")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))
REDIS_DB = int(os.getenv("REDIS_DB", 1))  # Use DB 1 for job mappings

# Initialize Redis client
redis_client = redis.Redis(host=REDIS_HOST, port=REDIS_PORT, db=REDIS_DB, decode_responses=True)

# Initialize APIRouter with /api prefix
router = APIRouter(prefix="/api")

@router.post("/run-job/")
async def run_vntyper(
    bam_file: UploadFile = File(..., description="BAM file to process"),
    bai_file: UploadFile = File(None, description="Optional BAI index file"),
    thread: int = Form(4),
    reference_assembly: str = Form("hg38"),
    fast_mode: bool = Form(False),
    keep_intermediates: bool = Form(False),
    archive_results: bool = Form(False)
):
    """
    Endpoint to run VNtyper job with additional parameters.
    Accepts a BAM file and an optional BAI index file.
    """
    logger.info("Received job submission")
    # Ensure input and output directories exist
    os.makedirs(DEFAULT_INPUT_DIR, exist_ok=True)
    os.makedirs(DEFAULT_OUTPUT_DIR, exist_ok=True)

    # Generate a unique job ID
    job_id = str(uuid4())
    job_input_dir = os.path.join(DEFAULT_INPUT_DIR, job_id)
    job_output_dir = os.path.join(DEFAULT_OUTPUT_DIR, job_id)

    os.makedirs(job_input_dir, exist_ok=True)
    os.makedirs(job_output_dir, exist_ok=True)

    # Save the uploaded BAM file
    bam_path = os.path.join(job_input_dir, bam_file.filename)
    with open(bam_path, "wb") as f:
        shutil.copyfileobj(bam_file.file, f)
    logger.info(f"Saved uploaded BAM file to {bam_path}")

    # Save the uploaded BAI file if provided
    if bai_file:
        bai_path = os.path.join(job_input_dir, bai_file.filename)
        with open(bai_path, "wb") as f:
            shutil.copyfileobj(bai_file.file, f)
        logger.info(f"Saved uploaded BAI file to {bai_path}")
    else:
        bai_path = None

    # Enqueue the Celery task
    task = run_vntyper_job.delay(
        bam_path=bam_path,
        output_dir=job_output_dir,
        thread=thread,
        reference_assembly=reference_assembly,
        fast_mode=fast_mode,
        keep_intermediates=keep_intermediates,
        archive_results=archive_results
    )
    logger.info(f"Enqueued job {job_id} with task ID {task.id}")

    # Store the mapping between job_id and task.id in Redis with a TTL (e.g., 7 days)
    redis_client.set(job_id, task.id, ex=604800)  # 7 days in seconds

    return {"message": "Job submitted", "job_id": job_id}

@router.get("/job-status/{job_id}/")
def get_job_status(job_id: str):
    """
    Endpoint to get the status of a job.
    """
    # Retrieve task ID from Redis using job_id
    task_id = redis_client.get(job_id)
    if not task_id:
        logger.warning(f"Job ID {job_id} not found in Redis")
        raise HTTPException(status_code=404, detail="Job ID not found")

    # Get the task result using Celery's AsyncResult
    task_result = AsyncResult(task_id)

    status = task_result.status
    logger.info(f"Job {job_id} (Task ID: {task_id}) status queried: {status}")

    if status == 'PENDING':
        return {"job_id": job_id, "status": "pending"}
    elif status == 'STARTED':
        return {"job_id": job_id, "status": "started"}
    elif status == 'SUCCESS':
        return {"job_id": job_id, "status": "completed"}
    elif status == 'FAILURE':
        return {"job_id": job_id, "status": "failed", "error": str(task_result.info)}
    else:
        return {"job_id": job_id, "status": status}

@router.get("/download/{job_id}/")
def download_result(job_id: str):
    """
    Endpoint to download the result (zipped file).
    """
    logger.info(f"Received download request for job {job_id}")
    zip_path = os.path.join(DEFAULT_OUTPUT_DIR, f"{job_id}.zip")
    if os.path.exists(zip_path):
        logger.info(f"Serving zip file {zip_path}")
        return FileResponse(
            zip_path,
            media_type="application/zip",
            filename=f"{job_id}.zip"
        )
    logger.warning(f"File not found: {zip_path}")
    raise HTTPException(status_code=404, detail="File not found")

@router.get("/health/")
def health_check():
    """
    Simple health check endpoint.
    """
    return {"status": "ok"}

# Include the router in the FastAPI app
app.include_router(router)
