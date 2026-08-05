import logging
import os
import shutil
import subprocess
from collections import Counter
from pathlib import Path
from typing import Optional, Union, List
from uuid import uuid4
from enum import Enum

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
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi_limiter import FastAPILimiter
from fastapi_limiter.depends import RateLimiter
from email_validator import validate_email, EmailNotValidError
from pydantic import BaseModel, Field
from starlette.background import BackgroundTask

from .cohorts import cohort_job_ids, create_cohort_record, resolve_cohort
from .config import build_redis_url, get_redis_password, require_redis_password, settings
from .tasks import run_vntyper_job, run_cohort_analysis_job
from .uploads import INDEX_EXTENSIONS, safe_upload_path, save_upload_bounded
from .version import API_VERSION

logger = logging.getLogger(__name__)


# Define valid reference assembly options
class ReferenceAssembly(str, Enum):
    HG19 = "hg19"
    HG38 = "hg38"
    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"


# Environment variables for default directories
DEFAULT_INPUT_DIR = settings.DEFAULT_INPUT_DIR
DEFAULT_OUTPUT_DIR = settings.DEFAULT_OUTPUT_DIR

# Total bytes a single job submission may write into the input directory above.
MAX_UPLOAD_BYTES = settings.MAX_UPLOAD_BYTES

# Redis configuration
REDIS_HOST = os.getenv("REDIS_HOST", "redis")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))

# Redis credential, from the single accessor in config.py so the API, the worker
# and beat all resolve the same value. There is no fallback; startup_event()
# calls require_redis_password() so a deployment that never set it stops there.
REDIS_PASSWORD = get_redis_password()

# Redis DBs
REDIS_DB = int(os.getenv("REDIS_DB", 1))  # Job mappings
RATE_LIMITING_REDIS_DB = settings.RATE_LIMITING_REDIS_DB
COHORT_REDIS_DB = int(os.getenv("COHORT_REDIS_DB", 3))  # Cohort data
USAGE_REDIS_DB = settings.USAGE_REDIS_DB  # Usage statistics

# Redis clients
redis_client = redis.Redis(
    host=REDIS_HOST,
    port=REDIS_PORT,
    db=REDIS_DB,
    password=REDIS_PASSWORD,
    decode_responses=True,
)
redis_cohort_client = redis.Redis(
    host=REDIS_HOST,
    port=REDIS_PORT,
    db=COHORT_REDIS_DB,
    password=REDIS_PASSWORD,
    decode_responses=True,
)
redis_usage_client = redis.Redis(
    host=REDIS_HOST,
    port=REDIS_PORT,
    db=USAGE_REDIS_DB,
    password=REDIS_PASSWORD,
    decode_responses=True,
)

# Global variable to store tool version
TOOL_VERSION = "unknown"

# Client-facing text for a job that failed. The task's own exception carries the
# command line and the per-job container paths, which a job-status poller has no
# use for; it is logged instead of returned.
JOB_FAILURE_MESSAGE = (
    "The job failed during processing. Quote the job ID when reporting it so the "
    "failure can be looked up in the server logs."
)

app = FastAPI(
    title="VNtyper Online API",
    version=API_VERSION,
    description=(
        """
        VNtyper Online API is an Application Programming Interface designed to facilitate the genotyping of MUC1 Variable Number Tandem Repeats (VNTR) in Autosomal Dominant Tubulointerstitial Kidney Disease (ADTKD-MUC1) using Short-Read Sequencing (SRS) data.

        This API allows users to submit genomic data for VNTR analysis, check job statuses, download results, and access aggregated usage statistics.

        Features

        - Submit Jobs: Upload BAM files and initiate VNTR analysis.
        - Job Management: Check the status of submitted jobs and retrieve results.
        - Cohort Support: Group jobs into cohorts for collective analysis.
        - Usage Statistics: Access anonymized usage statistics of the API.
        - In-Browser Processing: Leverages powerful genomic data processing tools.
        """
    ),
    terms_of_service="https://vntyper.org/terms/",
    contact={
        "name": "Support Team",
        "url": "https://vntyper.org/support/",
        "email": "support@vntyper.org",
    },
    license_info={
        "name": "BSD 3-Clause License",
        "url": "https://github.com/hassansaei/vntyper/blob/main/LICENSE",
    },
    root_path="/api",
    docs_url="/docs",
    redoc_url="/redoc",
    openapi_url="/openapi.json",
)

# CORS configuration for development
# Only allow localhost origins when ENVIRONMENT is 'development' or 'local'
ENVIRONMENT = os.getenv("ENVIRONMENT", "production")
if ENVIRONMENT in ["development", "local"]:
    app.add_middleware(
        CORSMiddleware,
        allow_origins=[
            "http://localhost:3000",
            "http://127.0.0.1:3000",
            "http://localhost:8000",
            "http://127.0.0.1:8000",
        ],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )
    logger.info(f"CORS enabled for {ENVIRONMENT} environment with origins: localhost:3000, localhost:8000")


@app.on_event("startup")
async def startup_event():
    """Initialize rate limiting and cache the VNtyper tool version.

    Raises:
        RuntimeError: If REDIS_PASSWORD is unset. Checked first, before anything
            opens a connection, so a misconfigured deployment fails here rather
            than coming up half-authenticated.
    """
    redis_password = require_redis_password()

    # Initialize Redis client for rate limiting
    try:
        rate_limit_redis_url = build_redis_url(
            REDIS_HOST, REDIS_PORT, RATE_LIMITING_REDIS_DB, redis_password
        )
        redis_rate_limit = aioredis.from_url(
            rate_limit_redis_url, encoding="utf8", decode_responses=True
        )
        await FastAPILimiter.init(redis_rate_limit)
        logger.info("Rate limiting initialized successfully.")
    except Exception as e:
        logger.error(f"Failed to initialize rate limiting: {e}")
        raise

    # Cache the VNtyper tool version
    global TOOL_VERSION
    try:
        tool_version_output = subprocess.check_output(
            ["vntyper", "-v"],
            stderr=subprocess.STDOUT,
            text=True,
            timeout=5,  # Timeout after 5 seconds to prevent hanging
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


# Define separate RateLimiter dependencies
simple_rate_limiter = RateLimiter(
    times=settings.RATE_LIMIT_SIMPLE_TIMES, seconds=settings.RATE_LIMIT_SIMPLE_SECONDS
)
high_rate_limiter = RateLimiter(
    times=settings.RATE_LIMIT_HIGH_TIMES, seconds=settings.RATE_LIMIT_HIGH_SECONDS
)

# Initialize APIRouter without prefix
router = APIRouter()


@router.get(
    "/version/",
    tags=["General"],
    dependencies=[Depends(simple_rate_limiter)],
    summary="Get API and Tool Version",
    description=(
        "Retrieve the current version of the API and the VNtyper tool.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_SIMPLE_TIMES} requests per "
        f"{settings.RATE_LIMIT_SIMPLE_SECONDS} seconds."
    ),
)
def get_versions():
    """
    **Description:**

    This endpoint returns the version information of the API and the VNtyper tool.

    **Returns:**

    - **api_version**: The version of the API.
    - **tool_version**: The version of the VNtyper tool.
    """
    return {"api_version": API_VERSION, "tool_version": TOOL_VERSION}


@router.post(
    "/create-cohort/",
    tags=["Cohort Management"],
    dependencies=[Depends(simple_rate_limiter)],
    summary="Create a new cohort",
    description=(
        "Create a new cohort. A passphrase is required, and it is the only "
        "credential that will open the cohort afterwards. An alias may be given "
        "as a label; aliases are unique, and an alias never identifies a cohort "
        "on its own.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_SIMPLE_TIMES} requests per {settings.RATE_LIMIT_SIMPLE_SECONDS} seconds."
    ),
)
def create_cohort(
    passphrase: str = Form(..., description="Passphrase protecting the cohort"),
    alias: Optional[str] = Form(None, description="Optional cohort alias"),
):
    """
    **Description:**

    This endpoint allows users to create a new cohort for grouping jobs. The passphrase is required: it is the only credential that will open the cohort for reading, downloading or analysis. An optional alias labels the cohort; aliases are unique across cohorts.

    **Parameters:**

    - **passphrase**: The passphrase that will protect cohort access. Required.
    - **alias**: An optional user-defined label for the cohort.

    **Returns:**

    - **cohort_id**: A unique identifier for the created cohort. Keep it: it is
      required, together with the passphrase, for every later request.
    - **alias**: The alias provided for the cohort.
    """
    try:
        return create_cohort_record(
            redis_cohort_client,
            alias=alias,
            passphrase=passphrase,
            retention_seconds=settings.COHORT_RETENTION_DAYS * 86400,
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


class RunJobResponse(BaseModel):
    message: str = Field(
        ..., description="Confirmation message indicating job submission."
    )
    job_id: str = Field(..., description="Unique identifier for the submitted job.")
    cohort_id: Optional[str] = Field(
        None, description="Identifier of the associated cohort, if any."
    )


@router.post(
    "/run-job/",
    tags=["Job Management"],
    dependencies=[Depends(high_rate_limiter)],
    summary="Submit a VNtyper job",
    description=(
        "Submit a VNtyper job with additional parameters. "
        "You can upload BAM files and configure various settings for the analysis. "
        "An optional email can be provided to receive notifications upon job completion.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_HIGH_TIMES} requests per {settings.RATE_LIMIT_HIGH_SECONDS} seconds."
    ),
    response_model=RunJobResponse,
)
async def run_vntyper(
    request: Request,
    bam_file: UploadFile = File(..., description="BAM file to process"),
    bai_file: UploadFile = File(None, description="Optional BAI index file"),
    thread: int = Form(4),
    reference_assembly: ReferenceAssembly = Form(ReferenceAssembly.HG38),
    fast_mode: bool = Form(False),
    keep_intermediates: bool = Form(False),
    archive_results: bool = Form(False),
    email: Optional[str] = Form(None, description="Optional email to receive results"),
    cohort_id: Optional[str] = Form(
        None, description="Optional cohort identifier to associate the job"
    ),
    alias: Optional[str] = Form(None, description="Optional cohort alias"),
    passphrase: Optional[str] = Form(
        None, description="Passphrase if required by the cohort"
    ),
    # ----------------------------------------------------
    # ADDED: single option for advntr_mode (default False)
    # ----------------------------------------------------
    advntr_mode: bool = Form(False),
):
    """
    **Description:**

    This endpoint allows users to submit a VNtyper job for VNTR analysis. Users can upload BAM files and specify analysis parameters such as threading, reference assembly, and more.

    **Parameters:**    - **bam_file**: The BAM file to be processed.
    - **bai_file**: Optional BAI index file corresponding to the BAM file.
    - **thread**: Number of threads to use for processing.
    - **reference_assembly**: Reference genome assembly to use ('hg19', 'hg38', 'GRCh37', 'GRCh38').
    - **fast_mode**: Boolean flag to enable fast mode processing.
    - **keep_intermediates**: Boolean flag to keep intermediate files.
    - **archive_results**: Boolean flag to archive results.
    - **email**: Optional email address to receive job notifications.
    - **cohort_id**: Optional cohort identifier to associate the job.
    - **alias**: Optional alias for the cohort.
    - **passphrase**: Passphrase if required by the cohort.

    **Returns:**

    - **message**: Confirmation message indicating job submission.
    - **job_id**: Unique identifier for the submitted job.
    - **cohort_id**: Identifier of the associated cohort, if any.
    """
    logger.info("Received job submission")

    # A request that already declares more than the ceiling is answered before
    # anything is written. The header is chosen by the client, so this is only a
    # cheap early answer -- the limit that counts is applied to the bytes
    # actually read, further down. Content-Length covers the whole body, which
    # is never smaller than the files inside it, so answering on it here cannot
    # let an oversized submission through.
    declared_length = request.headers.get("content-length")
    if declared_length and declared_length.isdigit() and int(declared_length) > MAX_UPLOAD_BYTES:
        msg = f"Upload exceeds the maximum accepted size of {MAX_UPLOAD_BYTES} bytes"
        logger.error(msg)
        raise HTTPException(status_code=413, detail=msg)

    # Ensure input and output directories exist
    os.makedirs(DEFAULT_INPUT_DIR, exist_ok=True)
    os.makedirs(DEFAULT_OUTPUT_DIR, exist_ok=True)

    # Generate a unique job ID
    job_id = str(uuid4())
    job_input_dir = os.path.join(DEFAULT_INPUT_DIR, job_id)
    job_output_dir = os.path.join(DEFAULT_OUTPUT_DIR, job_id)

    # Constrain the client-supplied filenames before anything is created on
    # disk, so a rejected submission leaves nothing behind.
    try:
        bam_path = safe_upload_path(job_input_dir, bam_file.filename)
        bai_path = None
        if bai_file is not None and bai_file.filename:
            bai_path = safe_upload_path(job_input_dir, bai_file.filename, INDEX_EXTENSIONS)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    os.makedirs(job_input_dir, exist_ok=True)
    os.makedirs(job_output_dir, exist_ok=True)

    # Save the uploaded files, counting the bytes as they are written. One
    # budget covers the whole submission, so the two parts together stay within
    # the same ceiling the declared length was checked against. A refused copy
    # takes the directories this request created with it, so nothing of a
    # rejected submission stays on the shared volume.
    try:
        remaining = MAX_UPLOAD_BYTES - save_upload_bounded(
            bam_file.file, bam_path, MAX_UPLOAD_BYTES
        )
        logger.info(f"Saved uploaded BAM file to {bam_path}")

        if bai_path is not None:
            save_upload_bounded(bai_file.file, bai_path, remaining)
            logger.info(f"Saved uploaded BAI file to {bai_path}")
    except ValueError as exc:
        shutil.rmtree(job_input_dir, ignore_errors=True)
        shutil.rmtree(job_output_dir, ignore_errors=True)
        msg = f"Upload exceeds the maximum accepted size of {MAX_UPLOAD_BYTES} bytes"
        raise HTTPException(status_code=413, detail=msg) from exc

    # Validate the email if provided
    if email:
        try:
            valid = validate_email(email)
            email = valid.email  # Get the normalized form
            logger.info(f"Validated email: {email}")
        except EmailNotValidError as e:
            logger.error(f"Invalid email address provided: {email} - {str(e)}")
            raise HTTPException(
                status_code=400, detail="Invalid email address provided."
            )

    # Cohort handling logic. Joining a cohort writes into it, so it is
    # authorized on exactly the same terms as reading one.
    if cohort_id or alias:
        cohort = authorized_cohort(cohort_id, alias, passphrase)
        cohort_id = cohort["cohort_id"]
        cohort_key = cohort["cohort_key"]
        # Assign the job to the cohort immediately
        redis_cohort_client.sadd(f"{cohort_key}:jobs", job_id)
        logger.info(f"Job {job_id} is associated with cohort {cohort_id}")
    else:
        cohort_key = None  # Job is not associated with any cohort
        logger.info(
            f"Job {job_id} submitted as a single job without cohort assignment."
        )

    # Extract client IP and User-Agent
    client_ip = request.client.host
    user_agent = request.headers.get(
        "User-Agent", "unknown"
    )  # ---------------------------------------------------------------------
    # If advntr_mode is True, use vntyper_long_queue; else the default queue.
    # ---------------------------------------------------------------------
    if advntr_mode:
        task = run_vntyper_job.apply_async(
            kwargs={
                "bam_path": bam_path,
                "output_dir": job_output_dir,
                "thread": thread,
                "reference_assembly": reference_assembly.value,
                "fast_mode": fast_mode,
                "keep_intermediates": keep_intermediates,
                "archive_results": archive_results,
                "email": email,
                "cohort_key": cohort_key,
                "client_ip": client_ip,
                "user_agent": user_agent,
                "advntr_mode": True,
            },
            queue="vntyper_long_queue",
        )
        logger.info(
            f"Enqueued adVNTR job {job_id} in long queue with task ID {task.id}"
        )
    else:
        task = run_vntyper_job.delay(
            bam_path=bam_path,
            output_dir=job_output_dir,
            thread=thread,
            reference_assembly=reference_assembly.value,
            fast_mode=fast_mode,
            keep_intermediates=keep_intermediates,
            archive_results=archive_results,
            email=email,
            cohort_key=cohort_key,
            client_ip=client_ip,
            user_agent=user_agent,
            advntr_mode=False,
        )
        logger.info(f"Enqueued job {job_id} with task ID {task.id}")

    # Store the mapping between job_id and task.id in Redis with a TTL (e.g., 7 days)
    redis_client.set(job_id, task.id, ex=604800)  # 7 days in seconds

    # Add the task ID to a Redis list to track the queue
    redis_client.rpush("vntyper_job_queue", task.id)

    return RunJobResponse(message="Job submitted", job_id=job_id, cohort_id=cohort_id)


class JobStatusResponse(BaseModel):
    job_id: str = Field(..., description="Unique identifier for the job.")
    status: str = Field(..., description="Current status of the job.")
    error: Optional[str] = Field(None, description="Error message if the job failed.")


@router.get(
    "/job-status/{job_id}/",
    tags=["Job Management"],
    dependencies=[Depends(simple_rate_limiter)],
    summary="Get the status of a VNtyper job",
    description=(
        "Retrieve the current status of a submitted VNtyper job using its job ID. "
        "Possible statuses include 'pending', 'started', 'completed', and 'failed'.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_SIMPLE_TIMES} requests per {settings.RATE_LIMIT_SIMPLE_SECONDS} seconds."
    ),
    response_model=JobStatusResponse,
)
def get_job_status(job_id: str):
    """
    **Description:**

    This endpoint retrieves the status of a VNtyper job using its unique job ID.

    **Parameters:**

    - **job_id**: The unique identifier of the job.

    **Returns:**

    - **job_id**: The unique identifier of the job.
    - **status**: The current status of the job.
    - **error**: Error message if the job has failed.
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
        return JobStatusResponse(job_id=job_id, status="pending")
    elif status == "STARTED":
        return JobStatusResponse(job_id=job_id, status="started")
    elif status == "SUCCESS":
        return JobStatusResponse(job_id=job_id, status="completed")
    elif status == "FAILURE":
        # Log the task's own exception, return a generic message: the detail is
        # operator-facing, and this endpoint is unauthenticated.
        logger.error(
            f"Job {job_id} (Task ID: {task_id}) failed: "
            f"{type(task_result.info).__name__}: {task_result.info}"
        )
        return JobStatusResponse(
            job_id=job_id, status="failed", error=JOB_FAILURE_MESSAGE
        )
    else:
        return JobStatusResponse(job_id=job_id, status=status.lower())


@router.get(
    "/download/{job_id}/",
    tags=["Job Management"],
    dependencies=[Depends(high_rate_limiter)],
    summary="Download the result of a VNtyper job",
    description=(
        "Download the zipped result files of a completed VNtyper job using its job ID. "
        "This endpoint is rate-limited to prevent abuse.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_HIGH_TIMES} requests per {settings.RATE_LIMIT_HIGH_SECONDS} seconds."
    ),
    responses={
        200: {
            "content": {"application/zip": {}},
            "description": "A ZIP file containing the job results.",
        },
        404: {"description": "File not found."},
    },
)
def download_result(job_id: str):
    """
    **Description:**

    This endpoint allows users to download the results of a completed VNtyper job as a ZIP file.

    **Parameters:**

    - **job_id**: The unique identifier of the job.

    **Returns:**

    - **FileResponse**: A ZIP file containing the results of the job.
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
    tags=["General"],
    dependencies=[Depends(simple_rate_limiter)],
    summary="Health check endpoint",
    description=(
        "Endpoint to check the health status of the API.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_SIMPLE_TIMES} requests per {settings.RATE_LIMIT_SIMPLE_SECONDS} seconds."
    ),
)
def health_check():
    """
    **Description:**

    Simple health check endpoint to verify that the API is running.

    **Returns:**

    - **status**: A message indicating the health status of the API.
    """
    return {"status": "ok"}


class JobQueueResponse(BaseModel):
    total_jobs_in_queue: int = Field(
        ..., description="Total number of jobs in the queue."
    )


class JobQueuePositionResponse(BaseModel):
    job_id: str = Field(..., description="Unique identifier for the job.")
    position_in_queue: Optional[int] = Field(
        None, description="Position of the job in the queue."
    )
    total_jobs_in_queue: int = Field(
        ..., description="Total number of jobs in the queue."
    )
    status: Optional[str] = Field(
        None, description="Status message if job is not in the queue."
    )


@router.get(
    "/job-queue/",
    tags=["Job Management"],
    dependencies=[Depends(simple_rate_limiter)],
    summary="Get job queue information",
    description=(
        "Retrieve the total number of jobs in the queue, or the position of a specific job.\n\n"
        "If no `job_id` is provided, returns the total number of jobs in the queue.\n"
        "If a `job_id` is provided, returns the position of that job in the queue.\n\n"
        "**Note:** This endpoint is rate-limited to prevent abuse.\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_SIMPLE_TIMES} requests per {settings.RATE_LIMIT_SIMPLE_SECONDS} seconds."
    ),
    response_model=Union[JobQueueResponse, JobQueuePositionResponse],
)
def get_job_queue(
    job_id: Optional[str] = None,
):
    """
    **Description:**

    This endpoint provides information about the job queue.

    **Parameters:**

    - **job_id**: (Optional) The unique identifier of a specific job.

    **Returns:**

    - If `job_id` is not provided:
        - **total_jobs_in_queue**: Total number of jobs currently in the queue.
    - If `job_id` is provided:
        - **job_id**: The unique identifier of the job.
        - **position_in_queue**: The position of the job in the queue.
        - **total_jobs_in_queue**: Total number of jobs currently in the queue.
        - **status**: Status message if the job is not in the queue.
    """
    try:
        # Get the list of task IDs from the Redis list
        task_ids = redis_client.lrange("vntyper_job_queue", 0, -1)
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
                return JobQueuePositionResponse(
                    job_id=job_id,
                    position_in_queue=position,
                    total_jobs_in_queue=queue_length,
                )
            else:
                # The job is not in the queue; it might be processing or completed
                return JobQueuePositionResponse(
                    job_id=job_id,
                    position_in_queue=None,
                    total_jobs_in_queue=queue_length,
                    status="Job not in queue (might be processing or completed)",
                )
        except HTTPException:
            # The "unknown job id" answer above is an HTTPException, and
            # HTTPException is an Exception, so without this it would be caught
            # below and reported as a server error. Let the handler's own
            # deliberate status codes through untouched.
            raise
        except Exception as e:
            logger.error(f"Error retrieving job position: {e}")
            raise HTTPException(status_code=500, detail="Error retrieving job position")
    else:
        # Return the total number of jobs in the queue
        return JobQueueResponse(total_jobs_in_queue=queue_length)


@router.get(
    "/cohort-status/",
    tags=["Cohort Management"],
    dependencies=[Depends(simple_rate_limiter)],
    summary="Get status of all jobs in a cohort",
    description=(
        "Retrieve the status of all jobs associated with a cohort.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_SIMPLE_TIMES} requests per {settings.RATE_LIMIT_SIMPLE_SECONDS} seconds."
    ),
)
def get_cohort_status(
    cohort_id: Optional[str] = Query(None, description="Cohort identifier (required)"),
    alias: Optional[str] = Query(None, description="Cohort alias, checked against the cohort"),
    passphrase: Optional[str] = Query(
        None, description="Passphrase protecting the cohort (required)"
    ),
):
    """
    **Description:**

    This endpoint retrieves the status of all jobs associated with a specific cohort.

    **Parameters:**

    - **cohort_id**: The unique identifier of the cohort. Required.
    - **alias**: (Optional) The alias of the cohort, checked against it.
    - **passphrase**: The passphrase protecting the cohort. Required.

    **Returns:**

    - **cohort_id**: The unique identifier of the cohort.
    - **alias**: The alias of the cohort.
    - **jobs**: A list of job statuses within the cohort.
    """
    # Reuse the get_cohort_jobs function to retrieve job_ids
    response = get_cohort_jobs(cohort_id=cohort_id, alias=alias, passphrase=passphrase)
    job_ids = response["job_ids"]

    # Define mapping from Celery statuses to standardized statuses
    celery_status_mapping = {
        "PENDING": "pending",
        "STARTED": "started",
        "SUCCESS": "completed",
        "FAILURE": "failed",
        # Add more mappings if necessary
    }

    # Get status for each job
    job_statuses = []
    for job_id in job_ids:
        task_id = redis_client.get(job_id)
        if not task_id:
            status = "unknown"
        else:
            task_result = AsyncResult(task_id)
            # Map the Celery status to standardized status
            status = celery_status_mapping.get(
                task_result.status, task_result.status.lower()
            )
        job_statuses.append({"job_id": job_id, "status": status})

    return {
        "cohort_id": response["cohort_id"],
        "alias": response["alias"],
        "jobs": job_statuses,
    }


class UsageStatisticsResponse(BaseModel):
    total_jobs: int = Field(..., description="Total number of jobs submitted.")
    unique_users: int = Field(..., description="Number of unique users.")
    job_statuses: dict = Field(..., description="Counts of jobs by status.")


@router.get(
    "/usage-statistics/",
    tags=["Statistics"],
    dependencies=[Depends(simple_rate_limiter)],
    summary="Get Usage Statistics",
    description=(
        "Retrieve aggregated usage statistics.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_SIMPLE_TIMES} requests per {settings.RATE_LIMIT_SIMPLE_SECONDS} seconds."
    ),
    response_model=UsageStatisticsResponse,
)
def get_usage_statistics():
    """
    **Description:**

    This endpoint provides aggregated usage statistics of the API, including total jobs submitted, unique users, and job statuses.

    **Returns:**

    - **total_jobs**: Total number of jobs submitted.
    - **unique_users**: Number of unique users who have submitted jobs.
    - **job_statuses**: A dictionary with job statuses as keys and counts as values.
    """
    # Retrieve all usage data
    usage_keys = redis_usage_client.keys("usage:*")
    total_jobs = len(usage_keys)
    unique_users = set()
    job_statuses = Counter()

    for key in usage_keys:
        data = redis_usage_client.hgetall(key)
        unique_users.add(data.get("user_hash"))
        status = data.get("status", "unknown")
        job_statuses[status] += 1

    return UsageStatisticsResponse(
        total_jobs=total_jobs,
        unique_users=len(unique_users),
        job_statuses=dict(job_statuses),
    )


def _remove_temp_file(path: str) -> None:
    """Delete a temporary file the service built for a single response.

    Tolerates the file already being gone so that a repeated cleanup cannot turn
    into an error raised after the response has been sent.

    Args:
        path: Filesystem path of the temporary file to remove.
    """
    Path(path).unlink(missing_ok=True)


@router.get(
    "/cohort-download/",
    tags=["Cohort Management"],
    dependencies=[Depends(high_rate_limiter)],
    summary="Download all cohort results as a single zip",
    description=(
        "Generate and download a single ZIP file containing all job result `.zip` "
        "files for the specified cohort.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_HIGH_TIMES} requests per {settings.RATE_LIMIT_HIGH_SECONDS} seconds."
    ),
)
def cohort_download(
    cohort_id: Optional[str] = Query(None, description="Cohort identifier (required)"),
    alias: Optional[str] = Query(None, description="Cohort alias, checked against the cohort"),
    passphrase: Optional[str] = Query(
        None, description="Passphrase protecting the cohort (required)"
    ),
):
    """
    **Description:**

    This endpoint creates a single ZIP file containing all job `.zip` result files
    for a given cohort, then serves it as a downloadable file.

    **Parameters:**

    - **cohort_id**: The unique identifier of the cohort. Required.
    - **alias**: (Optional) The alias of the cohort, checked against it.
    - **passphrase**: The passphrase protecting the cohort. Required.

    **Returns:**

    - **FileResponse**: A ZIP file containing all `.zip` results for the cohort.
    """
    import tempfile
    import zipfile

    # Retrieve job IDs for the given cohort
    response = get_cohort_jobs(cohort_id=cohort_id, alias=alias, passphrase=passphrase)
    job_ids = response["job_ids"]

    # Create a temporary ZIP file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".zip") as temp_zip:
        with zipfile.ZipFile(temp_zip.name, "w", zipfile.ZIP_DEFLATED) as zf:
            # For each job, look for its existing .zip file in the output directory
            for job_id in job_ids:
                job_zip_path = os.path.join(DEFAULT_OUTPUT_DIR, f"{job_id}.zip")
                if os.path.exists(job_zip_path):
                    zf.write(job_zip_path, arcname=os.path.basename(job_zip_path))

        final_zip_path = temp_zip.name  # We'll pass this to FileResponse

    # Suggest a download filename
    download_name = f"cohort_{response['cohort_id']}.zip"

    # Return the zipped file as a download. The archive is scratch space owned by
    # this one response, so it is removed in a background task, which Starlette
    # runs after the body has been sent -- the caller still receives the complete
    # archive, and nothing of it is left on disk afterwards.
    return FileResponse(
        path=final_zip_path,
        media_type="application/zip",
        filename=download_name,
        background=BackgroundTask(_remove_temp_file, final_zip_path),
    )


# ----------------------------------------------------------------------
# Feature #82: Joint Cohort Analysis Endpoint
# ----------------------------------------------------------------------
@router.post(
    "/cohort-analysis/",
    tags=["Cohort Management"],
    dependencies=[Depends(high_rate_limiter)],
    summary="Run a joint cohort analysis",
    description=(
        "Run a VNtyper joint cohort analysis using all .zip results files "
        "from each job in the specified cohort. "
        "Returns a new 'analysis job ID' for checking status and downloading results.\n\n"
        f"**Rate Limit:** {settings.RATE_LIMIT_HIGH_TIMES} requests per "
        f"{settings.RATE_LIMIT_HIGH_SECONDS} seconds."
    ),
)
def run_cohort_analysis(
    request: Request,
    cohort_id: Optional[str] = Form(None, description="Cohort identifier (required)"),
    alias: Optional[str] = Form(None, description="Cohort alias, checked against the cohort"),
    passphrase: Optional[str] = Form(
        None, description="Passphrase protecting the cohort (required)"
    ),
):
    """
    **Description:**

    This endpoint gathers the `.zip` files for all jobs in a specified cohort
    and performs a joint cohort analysis via `vntyper cohort`.

    **Parameters:**

    - **cohort_id**: The unique identifier of the cohort. Required.
    - **alias**: (Optional) The alias of the cohort, checked against it.
    - **passphrase**: The passphrase protecting the cohort. Required.

    **Returns:**
    - **message**: Confirmation message.
    - **analysis_job_id**: A unique identifier for this cohort analysis job.
    """
    logger.info("Received cohort analysis request")

    # 1) Retrieve the set of job IDs for this cohort
    response = get_cohort_jobs(cohort_id=cohort_id, alias=alias, passphrase=passphrase)
    cid = response["cohort_id"]
    job_ids = response["job_ids"]
    if not job_ids:
        raise HTTPException(
            status_code=400, detail="No jobs found in the specified cohort."
        )

    # 2) Build the list of existing .zip result paths for these jobs
    zip_paths: List[str] = []
    for jid in job_ids:
        candidate_zip = os.path.join(DEFAULT_OUTPUT_DIR, f"{jid}.zip")
        if os.path.exists(candidate_zip):
            zip_paths.append(candidate_zip)
        else:
            logger.warning(
                f"Job {jid} missing .zip file. Skipping from cohort analysis."
            )

    if not zip_paths:
        raise HTTPException(
            status_code=400,
            detail="No .zip files found for the specified cohort's jobs.",
        )

    # 3) Create a unique "analysis job ID" (directory) for results
    analysis_job_id = str(uuid4())
    analysis_output_dir = os.path.join(DEFAULT_OUTPUT_DIR, analysis_job_id)

    # 4) Enqueue the new Celery task
    client_ip = request.client.host
    user_agent = request.headers.get("User-Agent", "unknown")
    task = run_cohort_analysis_job.delay(
        cohort_id=cid,
        zip_paths=zip_paths,
        output_dir=analysis_output_dir,
        user_ip=client_ip,
        user_agent=user_agent,
    )
    logger.info(
        f"Enqueued cohort analysis job: analysis_job_id={analysis_job_id}, task_id={task.id}"
    )

    # Store the mapping for job status checking
    redis_client.set(analysis_job_id, task.id, ex=604800)  # 7 days TTL
    redis_client.rpush("vntyper_job_queue", task.id)

    return {
        "message": "Cohort analysis started",
        "analysis_job_id": analysis_job_id,
    }


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


def authorized_cohort(
    cohort_id: Optional[str], alias: Optional[str], passphrase: Optional[str]
):
    """
    Resolve a cohort and authorize the caller, as HTTP.

    The rules live in `cohorts.resolve_cohort`, which signals with builtin
    exceptions so it stays free of FastAPI. This is the only place that
    vocabulary is translated into status codes, so every cohort route answers a
    given failure the same way.

    Args:
        cohort_id: The cohort's identifier. Required.
        alias: Optional label, cross-checked against the resolved cohort.
        passphrase: The cohort's passphrase. Required.

    Returns:
        dict: The resolved `cohort_id`, `cohort_key` and `alias`.

    Raises:
        HTTPException: 400 for a malformed request, 401 for a refused one, 404
            for a cohort that does not exist.
    """
    try:
        return resolve_cohort(
            redis_cohort_client,
            cohort_id=cohort_id,
            alias=alias,
            passphrase=passphrase,
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except PermissionError as exc:
        raise HTTPException(status_code=401, detail=str(exc)) from exc
    except LookupError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


def get_cohort_jobs(
    cohort_id: Optional[str], alias: Optional[str], passphrase: Optional[str]
):
    """
    Helper function to retrieve job IDs associated with a cohort.

    Args:
        cohort_id: The cohort's identifier. Required.
        alias: Optional label, cross-checked against the resolved cohort.
        passphrase: The cohort's passphrase. Required.

    Returns:
        dict: The cohort's `cohort_id`, `alias` and member `job_ids`.
    """
    cohort = authorized_cohort(cohort_id, alias, passphrase)
    return {
        "cohort_id": cohort["cohort_id"],
        "alias": cohort["alias"],
        "job_ids": cohort_job_ids(redis_cohort_client, cohort["cohort_key"]),
    }
