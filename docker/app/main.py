# docker/app/main.py

import logging
import os
import shutil
import subprocess
from typing import Optional
from uuid import uuid4

import redis
import redis.asyncio as aioredis
from celery.result import AsyncResult
from fastapi import (
    APIRouter,
    Depends,
    FastAPI,
    File,
    Form,
    HTTPException,
    Request,
    UploadFile,
)
from fastapi.responses import FileResponse, JSONResponse
from fastapi_limiter import FastAPILimiter
from fastapi_limiter.depends import RateLimiter
from pydantic import EmailStr

from .config import settings
from .tasks import run_vntyper_job
from .version import API_VERSION

app = FastAPI(
    title=settings.PROJECT_NAME,
    version=API_VERSION,
    root_path="/api",
    docs_url="/docs",
    redoc_url="/redoc",
    openapi_url="/openapi.json",
)

logger = logging.getLogger(__name__)

# Environment variables for default directories
DEFAULT_INPUT_DIR = settings.DEFAULT_INPUT_DIR
DEFAULT_OUTPUT_DIR = settings.DEFAULT_OUTPUT_DIR

# Redis configuration for job_id to task_id mapping
REDIS_HOST = os.getenv("REDIS_HOST", "redis")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))
REDIS_DB = int(os.getenv("REDIS_DB", 1))  # Use DB 1 for job mappings

# Rate limiting Redis configuration
RATE_LIMITING_REDIS_DB = settings.RATE_LIMITING_REDIS_DB
RATE_LIMIT_REDIS_URL = f"redis://{REDIS_HOST}:{REDIS_PORT}/{RATE_LIMITING_REDIS_DB}"

# Initialize Redis client for job_id mapping
redis_client = redis.Redis(
    host=REDIS_HOST, port=REDIS_PORT, db=REDIS_DB, decode_responses=True
)

# Global variable to store tool version
TOOL_VERSION = "unknown"


@app.on_event("startup")
async def startup_event():
    """Initialize rate limiting and cache the VNtyper tool version."""
    # Initialize Redis client for rate limiting
    redis_rate_limit = aioredis.from_url(
        RATE_LIMIT_REDIS_URL, encoding="utf8", decode_responses=True
    )
    await FastAPILimiter.init(redis_rate_limit)

    # Cache the VNtyper tool version
    global TOOL_VERSION
    try:
        tool_version_output = subprocess.check_output(
            ['vntyper', '-v'],
            stderr=subprocess.STDOUT,
            text=True,
            timeout=5  # Timeout after 5 seconds to prevent hanging
        )
        TOOL_VERSION = tool_version_output.strip()
        logger.info(f"VNtyper tool version: {TOOL_VERSION}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error retrieving tool version: {e.output.strip()}")
        TOOL_VERSION = "error retrieving tool version"
    except FileNotFoundError:
        logger.error("VNtyper tool not found.")
        TOOL_VERSION = "VNtyper tool not installed"
    except subprocess.TimeoutExpired:
        logger.error("Timeout expired while retrieving tool version.")
        TOOL_VERSION = "timeout retrieving tool version"


# Initialize APIRouter without changing URL paths
router = APIRouter(
    prefix="/api",
)


@router.get(
    "/version/",
    dependencies=[
        Depends(
            RateLimiter(
                times=settings.RATE_LIMIT_TIMES, seconds=settings.RATE_LIMIT_SECONDS
            )
        )
    ],
    summary="Get API and Tool Version",
    description=(
        "Retrieve the current version of the API and the VNtyper tool.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_TIMES} requests per "
        f"{settings.RATE_LIMIT_SECONDS} seconds."
    ),
)
def get_versions():
    """
    Endpoint to get the API and VNtyper tool versions.
    """
    return {"api_version": API_VERSION, "tool_version": TOOL_VERSION}


@router.post(
    "/run-job/",
    dependencies=[
        Depends(
            RateLimiter(
                times=settings.RATE_LIMIT_TIMES, seconds=settings.RATE_LIMIT_SECONDS
            )
        )
    ],
    summary="Submit a VNtyper job",
    description=(
        "Submit a VNtyper job with additional parameters. "
        "This endpoint is rate-limited to prevent abuse.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_TIMES} requests per "
        f"{settings.RATE_LIMIT_SECONDS} seconds."
    ),
)
async def run_vntyper(
    bam_file: UploadFile = File(..., description="BAM file to process"),
    bai_file: UploadFile = File(None, description="Optional BAI index file"),
    thread: int = Form(4),
    reference_assembly: str = Form("hg38"),
    fast_mode: bool = Form(False),
    keep_intermediates: bool = Form(False),
    archive_results: bool = Form(False),
    email: Optional[EmailStr] = Form(None, description="Optional email to receive results"),
):
    """
    Endpoint to run VNtyper job with additional parameters.
    Accepts a BAM file, an optional BAI index file, and an optional email.

    **Rate Limit:** {settings.RATE_LIMIT_TIMES} requests per {settings.RATE_LIMIT_SECONDS} seconds.
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

    # Enqueue the Celery task with email parameter
    task = run_vntyper_job.delay(
        bam_path=bam_path,
        output_dir=job_output_dir,
        thread=thread,
        reference_assembly=reference_assembly,
        fast_mode=fast_mode,
        keep_intermediates=keep_intermediates,
        archive_results=archive_results,
        email=email,
    )
    logger.info(f"Enqueued job {job_id} with task ID {task.id}")

    # Store the mapping between job_id and task.id in Redis with a TTL (e.g., 7 days)
    redis_client.set(job_id, task.id, ex=604800)  # 7 days in seconds

    # Add the task ID to a Redis list to track the queue
    redis_client.rpush('vntyper_job_queue', task.id)

    return {"message": "Job submitted", "job_id": job_id}


@router.get(
    "/job-status/{job_id}/",
    dependencies=[
        Depends(
            RateLimiter(
                times=settings.RATE_LIMIT_TIMES, seconds=settings.RATE_LIMIT_SECONDS
            )
        )
    ],
    summary="Get the status of a VNtyper job",
    description=(
        "Retrieve the current status of a submitted VNtyper job using its job ID. "
        "This endpoint is rate-limited to prevent abuse.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_TIMES} requests per "
        f"{settings.RATE_LIMIT_SECONDS} seconds."
    ),
)
def get_job_status(job_id: str):
    """
    Endpoint to get the status of a job.

    **Rate Limit:** {settings.RATE_LIMIT_TIMES} requests per {settings.RATE_LIMIT_SECONDS} seconds.
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

    if status == "PENDING":
        return {"job_id": job_id, "status": "pending"}
    elif status == "STARTED":
        return {"job_id": job_id, "status": "started"}
    elif status == "SUCCESS":
        return {"job_id": job_id, "status": "completed"}
    elif status == "FAILURE":
        return {"job_id": job_id, "status": "failed", "error": str(task_result.info)}
    else:
        return {"job_id": job_id, "status": status}


@router.get(
    "/download/{job_id}/",
    dependencies=[
        Depends(
            RateLimiter(
                times=settings.RATE_LIMIT_TIMES, seconds=settings.RATE_LIMIT_SECONDS
            )
        )
    ],
    summary="Download the result of a VNtyper job",
    description=(
        "Download the zipped result files of a completed VNtyper job using its job ID. "
        "This endpoint is rate-limited to prevent abuse.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_TIMES} requests per "
        f"{settings.RATE_LIMIT_SECONDS} seconds."
    ),
)
def download_result(job_id: str):
    """
    Endpoint to download the result (zipped file).

    **Rate Limit:** {settings.RATE_LIMIT_TIMES} requests per {settings.RATE_LIMIT_SECONDS} seconds.
    """
    logger.info(f"Received download request for job {job_id}")
    zip_path = os.path.join(DEFAULT_OUTPUT_DIR, f"{job_id}.zip")
    if os.path.exists(zip_path):
        logger.info(f"Serving zip file {zip_path}")
        return FileResponse(
            zip_path,
            media_type="application/zip",
            filename=f"{job_id}.zip",
        )
    logger.warning(f"File not found: {zip_path}")
    raise HTTPException(status_code=404, detail="File not found")


@router.get(
    "/health/",
    summary="Health check endpoint",
    description="Endpoint to check the health status of the API.",
)
def health_check():
    """
    Simple health check endpoint.
    """
    return {"status": "ok"}


@router.get(
    "/job-queue/",
    dependencies=[
        Depends(
            RateLimiter(
                times=settings.RATE_LIMIT_TIMES, seconds=settings.RATE_LIMIT_SECONDS
            )
        )
    ],
    summary="Get job queue information",
    description=(
        "Retrieve the total number of jobs in the queue, or the position of a specific job.\n\n"
        "If no `job_id` is provided, returns the total number of jobs in the queue.\n"
        "If a `job_id` is provided, returns the position of that job in the queue.\n\n"
        "**Note:** This endpoint is rate-limited to prevent abuse.\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_TIMES} requests per "
        f"{settings.RATE_LIMIT_SECONDS} seconds."
    ),
)
def get_job_queue(
    job_id: Optional[str] = None,
):
    """
    Endpoint to get job queue information.

    - If `job_id` is provided, returns the position of the job in the queue.
    - If `job_id` is not provided, returns the total number of jobs in the queue.

    **Rate Limit:** {settings.RATE_LIMIT_TIMES} requests per {settings.RATE_LIMIT_SECONDS} seconds.
    """
    try:
        # Get the list of task IDs from the Redis list
        task_ids = redis_client.lrange('vntyper_job_queue', 0, -1)
        queue_length = len(task_ids)
    except Exception as e:
        logger.error(f"Error accessing the job queue: {e}")
        raise HTTPException(status_code=500, detail="Error accessing the job queue")

    if job_id:
        try:
            # Retrieve the task ID associated with the provided job_id
            task_id = redis_client.get(job_id)
            if not task_id:
                logger.warning(f"Job ID {job_id} not found")
                raise HTTPException(status_code=404, detail="Job ID not found")

            if task_id in task_ids:
                position = task_ids.index(task_id) + 1  # Positions start at 1
                return {
                    "job_id": job_id,
                    "position_in_queue": position,
                    "total_jobs_in_queue": queue_length,
                }
            else:
                # The job is not in the queue; it might be processing or completed
                return {
                    "job_id": job_id,
                    "position_in_queue": None,
                    "status": "Job not in queue (might be processing or completed)",
                }
        except Exception as e:
            logger.error(f"Error retrieving job position: {e}")
            raise HTTPException(
                status_code=500, detail="Error retrieving job position"
            )
    else:
        # Return the total number of jobs in the queue
        return {"total_jobs_in_queue": queue_length}


# Include the router in the FastAPI app
app.include_router(router)


@app.exception_handler(HTTPException)
async def custom_http_exception_handler(request: Request, exc: HTTPException):
    """
    Custom exception handler for rate limiting and other HTTP exceptions.
    """
    if exc.status_code == 429:
        # Customize the error message for rate limit exceeded
        logger.warning(f"Rate limit exceeded for IP: {request.client.host}")
        return JSONResponse(
            status_code=exc.status_code,
            content={"detail": "Rate limit exceeded. Please try again later."},
        )
    # Handle other HTTP exceptions
    return JSONResponse(
        status_code=exc.status_code,
        content={"detail": exc.detail},
    )
