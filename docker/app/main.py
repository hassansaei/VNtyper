# docker/app/main.py

import logging
import os
import shutil
import subprocess
from typing import Optional
from uuid import uuid4
from datetime import datetime

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
    Query,
    Request,
    UploadFile,
)
from fastapi.responses import FileResponse, JSONResponse
from fastapi_limiter import FastAPILimiter
from fastapi_limiter.depends import RateLimiter
from email_validator import validate_email, EmailNotValidError
from passlib.context import CryptContext

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

# Redis configuration
REDIS_HOST = os.getenv("REDIS_HOST", "redis")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))

# Redis DBs
REDIS_DB = int(os.getenv("REDIS_DB", 1))  # Job mappings
RATE_LIMITING_REDIS_DB = settings.RATE_LIMITING_REDIS_DB
COHORT_REDIS_DB = int(os.getenv("COHORT_REDIS_DB", 3))  # Cohort data

# Redis clients
redis_client = redis.Redis(
    host=REDIS_HOST, port=REDIS_PORT, db=REDIS_DB, decode_responses=True
)
redis_cohort_client = redis.Redis(
    host=REDIS_HOST, port=REDIS_PORT, db=COHORT_REDIS_DB, decode_responses=True
)

# Rate limiting Redis URL
RATE_LIMIT_REDIS_URL = f"redis://{REDIS_HOST}:{REDIS_PORT}/{RATE_LIMITING_REDIS_DB}"

# Global variable to store tool version
TOOL_VERSION = "unknown"

# Initialize password hashing context
pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")


def hash_passphrase(passphrase: str) -> str:
    return pwd_context.hash(passphrase)


def verify_passphrase(passphrase: str, hashed_passphrase: str) -> bool:
    return pwd_context.verify(passphrase, hashed_passphrase)


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


# Initialize APIRouter without prefix
router = APIRouter()  # Removed prefix="/api"

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
    "/create-cohort/",
    dependencies=[
        Depends(
            RateLimiter(
                times=settings.RATE_LIMIT_TIMES, seconds=settings.RATE_LIMIT_SECONDS
            )
        )
    ],
    summary="Create a new cohort",
    description=(
        "Create a new cohort with an optional alias and passphrase.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_TIMES} requests per {settings.RATE_LIMIT_SECONDS} seconds."
    ),
)
def create_cohort(
    alias: Optional[str] = Form(None, description="Optional cohort alias"),
    passphrase: Optional[str] = Form(None, description="Optional passphrase to protect the cohort"),
):
    # Generate a unique cohort identifier
    cohort_id = str(uuid4())

    # Hash the passphrase if provided
    hashed_passphrase = hash_passphrase(passphrase) if passphrase else None

    # Store cohort metadata in Redis hash
    cohort_key = f"cohort:{cohort_id}"
    redis_cohort_client.hset(cohort_key, mapping={
        "alias": alias or "",
        "hashed_passphrase": hashed_passphrase or "",
        "created_at": datetime.utcnow().isoformat()
    })

    # Set a TTL for the cohort based on retention policy
    ttl_seconds = settings.COHORT_RETENTION_DAYS * 86400
    redis_cohort_client.expire(cohort_key, ttl_seconds)

    return {"cohort_id": cohort_id, "alias": alias}


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
    email: Optional[str] = Form(None, description="Optional email to receive results"),
    cohort_id: Optional[str] = Form(None, description="Optional cohort identifier to associate the job"),
    alias: Optional[str] = Form(None, description="Optional cohort alias"),
    passphrase: Optional[str] = Form(None, description="Passphrase if required by the cohort"),
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

    # Validate the email if provided
    if email:
        try:
            valid = validate_email(email)
            email = valid.email  # Get the normalized form
            logger.info(f"Validated email: {email}")
        except EmailNotValidError as e:
            logger.error(f"Invalid email address provided: {email} - {str(e)}")
            raise HTTPException(status_code=400, detail="Invalid email address provided.")

    # Cohort handling
    if cohort_id or alias:
        # Retrieve the cohort using cohort_id or alias
        cohort_key = None
        cohort_data = None
        if cohort_id:
            cohort_key = f"cohort:{cohort_id}"
            cohort_data = redis_cohort_client.hgetall(cohort_key)
            if not cohort_data:
                raise HTTPException(status_code=404, detail="Cohort ID not found")
        elif alias:
            # Search for cohort by alias
            for key in redis_cohort_client.scan_iter("cohort:*"):
                data = redis_cohort_client.hgetall(key)
                if data.get("alias") == alias:
                    cohort_key = key
                    cohort_data = data
                    cohort_id = key.split(":", 1)[1]  # Extract cohort_id from key
                    break
            if not cohort_key:
                raise HTTPException(status_code=404, detail="Cohort alias not found")
        else:
            raise HTTPException(status_code=400, detail="Cohort identifier or alias required")

        # Verify passphrase if required
        if cohort_data.get("hashed_passphrase"):
            if not passphrase:
                raise HTTPException(status_code=401, detail="Passphrase required for this cohort")
            if not verify_passphrase(passphrase, cohort_data["hashed_passphrase"]):
                raise HTTPException(status_code=401, detail="Incorrect passphrase")
    else:
        cohort_key = None  # Job is not associated with any cohort

    # Enqueue the Celery task with email parameter and cohort information
    task = run_vntyper_job.delay(
        bam_path=bam_path,
        output_dir=job_output_dir,
        thread=thread,
        reference_assembly=reference_assembly,
        fast_mode=fast_mode,
        keep_intermediates=keep_intermediates,
        archive_results=archive_results,
        email=email,
        cohort_key=cohort_key,
    )
    logger.info(f"Enqueued job {job_id} with task ID {task.id}")

    # Store the mapping between job_id and task.id in Redis with a TTL (e.g., 7 days)
    redis_client.set(job_id, task.id, ex=604800)  # 7 days in seconds

    # Add the task ID to a Redis list to track the queue
    redis_client.rpush('vntyper_job_queue', task.id)

    # If associated with a cohort, add the job to the cohort's job set
    if cohort_key:
        redis_cohort_client.sadd(f"{cohort_key}:jobs", job_id)
        # Ensure the TTL is updated
        ttl_seconds = settings.COHORT_RETENTION_DAYS * 86400
        redis_cohort_client.expire(f"{cohort_key}:jobs", ttl_seconds)

    return {"message": "Job submitted", "job_id": job_id, "cohort_id": cohort_id}


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


@router.get(
    "/cohort-jobs/",
    dependencies=[
        Depends(
            RateLimiter(
                times=settings.RATE_LIMIT_TIMES, seconds=settings.RATE_LIMIT_SECONDS
            )
        )
    ],
    summary="Get all jobs in a cohort",
    description=(
        "Retrieve all job IDs associated with a cohort.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_TIMES} requests per {settings.RATE_LIMIT_SECONDS} seconds."
    ),
)
def get_cohort_jobs(
    cohort_id: Optional[str] = Query(None, description="Cohort identifier"),
    alias: Optional[str] = Query(None, description="Cohort alias"),
    passphrase: Optional[str] = Query(None, description="Passphrase if required by the cohort"),
):
    # Retrieve the cohort using cohort_id or alias
    cohort_key = None
    cohort_data = None
    if cohort_id:
        cohort_key = f"cohort:{cohort_id}"
        cohort_data = redis_cohort_client.hgetall(cohort_key)
        if not cohort_data:
            raise HTTPException(status_code=404, detail="Cohort ID not found")
    elif alias:
        # Search for cohort by alias
        for key in redis_cohort_client.scan_iter("cohort:*"):
            data = redis_cohort_client.hgetall(key)
            if data.get("alias") == alias:
                cohort_key = key
                cohort_data = data
                cohort_id = key.split(":", 1)[1]  # Extract cohort_id from key
                break
        if not cohort_key:
            raise HTTPException(status_code=404, detail="Cohort alias not found")
    else:
        raise HTTPException(status_code=400, detail="Cohort identifier or alias required")

    # Verify passphrase if required
    if cohort_data.get("hashed_passphrase"):
        if not passphrase:
            raise HTTPException(status_code=401, detail="Passphrase required for this cohort")
        if not verify_passphrase(passphrase, cohort_data["hashed_passphrase"]):
            raise HTTPException(status_code=401, detail="Incorrect passphrase")

    # Get all job IDs associated with the cohort
    job_ids = redis_cohort_client.smembers(f"{cohort_key}:jobs")

    return {"cohort_id": cohort_id, "alias": cohort_data.get("alias"), "job_ids": list(job_ids)}


@router.get(
    "/cohort-status/",
    dependencies=[
        Depends(
            RateLimiter(
                times=settings.RATE_LIMIT_TIMES, seconds=settings.RATE_LIMIT_SECONDS
            )
        )
    ],
    summary="Get status of all jobs in a cohort",
    description=(
        "Retrieve the status of all jobs associated with a cohort.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_TIMES} requests per {settings.RATE_LIMIT_SECONDS} seconds."
    ),
)
def get_cohort_status(
    cohort_id: Optional[str] = Query(None, description="Cohort identifier"),
    alias: Optional[str] = Query(None, description="Cohort alias"),
    passphrase: Optional[str] = Query(None, description="Passphrase if required by the cohort"),
):
    # Reuse the get_cohort_jobs function to retrieve job_ids
    response = get_cohort_jobs(cohort_id=cohort_id, alias=alias, passphrase=passphrase)
    job_ids = response["job_ids"]

    # Get status for each job
    job_statuses = []
    for job_id in job_ids:
        task_id = redis_client.get(job_id)
        if not task_id:
            status = "unknown"
        else:
            task_result = AsyncResult(task_id)
            status = task_result.status.lower()
        job_statuses.append({"job_id": job_id, "status": status})

    return {"cohort_id": response["cohort_id"], "alias": response["alias"], "jobs": job_statuses}


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
