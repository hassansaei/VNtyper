# docker/app/main.py

from fastapi import (
    FastAPI,
    APIRouter,
    UploadFile,
    File,
    Form,
    HTTPException,
    Depends,
    Request,
)
from fastapi.responses import FileResponse, JSONResponse
from fastapi.security import HTTPBasic, HTTPBasicCredentials
from typing import Optional
from uuid import uuid4
import os
import shutil
import logging
import json
import base64

from .tasks import run_vntyper_job
from .config import settings
from .version import API_VERSION  # Import the API version

from celery.result import AsyncResult
import redis
import redis.asyncio as aioredis
from fastapi_limiter import FastAPILimiter
from fastapi_limiter.depends import RateLimiter

app = FastAPI(
    title=settings.PROJECT_NAME,
    version=API_VERSION,  # Use API_VERSION here
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

security = HTTPBasic()


@app.on_event("startup")
async def startup():
    # Initialize Redis client for rate limiting
    redis_rate_limit = aioredis.from_url(
        RATE_LIMIT_REDIS_URL, encoding="utf8", decode_responses=True
    )
    await FastAPILimiter.init(redis_rate_limit)


# Initialize APIRouter without changing URL paths
router = APIRouter(
    prefix="/api",
)


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
        f"Submit a VNtyper job with additional parameters. "
        f"This endpoint is rate-limited to prevent abuse.\n\n"
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
):
    """
    Endpoint to run VNtyper job with additional parameters.
    Accepts a BAM file and an optional BAI index file.

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

    # Enqueue the Celery task
    task = run_vntyper_job.delay(
        bam_path=bam_path,
        output_dir=job_output_dir,
        thread=thread,
        reference_assembly=reference_assembly,
        fast_mode=fast_mode,
        keep_intermediates=keep_intermediates,
        archive_results=archive_results,
    )
    logger.info(f"Enqueued job {job_id} with task ID {task.id}")

    # Store the mapping between job_id and task.id in Redis with a TTL (e.g., 7 days)
    redis_client.set(job_id, task.id, ex=604800)  # 7 days in seconds

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
    credentials: HTTPBasicCredentials = Depends(security),
):
    """
    Endpoint to get job queue information.

    - If `job_id` is provided, returns the position of the job in the queue.
    - If `job_id` is not provided, returns the total number of jobs in the queue.

    **Rate Limit:** {settings.RATE_LIMIT_TIMES} requests per {settings.RATE_LIMIT_SECONDS} seconds.
    """
    # Security: Validate the credentials (placeholder)
    # TODO: Implement actual authentication logic
    # For now, assume all users are authorized

    # Connect to Redis broker (Celery uses db 0 by default)
    redis_broker = redis.Redis(host=REDIS_HOST, port=REDIS_PORT, db=0)

    queue_name = "celery"  # Default queue name

    try:
        # Retrieve the total number of tasks in the queue
        queue_length = redis_broker.llen(queue_name)
    except Exception as e:
        logger.error(f"Error accessing the Celery queue: {e}")
        raise HTTPException(status_code=500, detail="Error accessing the job queue")

    if job_id:
        try:
            # Retrieve all tasks in the queue
            tasks = redis_broker.lrange(queue_name, 0, -1)
            # Extract task IDs from the task messages
            task_ids = []
            for task_message in tasks:
                # Each task message is a Redis byte string
                # Decode and parse the message to extract the task ID
                task_message = task_message.decode("utf-8")
                task_data = json.loads(task_message)
                body = task_data.get("body")
                if body:
                    # The body is Base64 encoded JSON
                    decoded_body = base64.b64decode(body)
                    body_data = json.loads(decoded_body)
                    task_id_in_queue = body_data.get("id")
                    if task_id_in_queue:
                        task_ids.append(task_id_in_queue)

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

# Custom exception handler for rate limiting
@app.exception_handler(HTTPException)
async def custom_http_exception_handler(request: Request, exc: HTTPException):
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
