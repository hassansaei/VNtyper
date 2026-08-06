# VNtyper Docker Container

A production-ready Docker container for **VNtyper**, built as two images:

| Image | Contents | Rebuilt |
| --- | --- | --- |
| `docker/Dockerfile.base` | conda environments, adVNTR, reference genomes + BWA indexes | only when `conda/**`, `docker/requirements-web.txt`, `install_references*` or `vntyper/dependencies/advntr/**` change |
| `docker/Dockerfile` | the VNtyper application and the FastAPI app | every commit (~3 min) |

Splitting them keeps the per-commit build small: the expensive layers are published
once to `ghcr.io/hassansaei/vntyper-base` and reused.

## **Building the Docker Image**

### **Quick Start**

```bash
git clone https://github.com/hassansaei/VNtyper.git
cd VNtyper

# Builds the application on the published base image (pulls the base, ~3 min)
make docker-build
```

Equivalently, by hand:

```bash
DOCKER_BUILDKIT=1 docker build -f docker/Dockerfile \
  --build-arg BASE_IMAGE=ghcr.io/hassansaei/vntyper-base:latest \
  -t vntyper:latest .
```

### **Rebuilding the base image**

Only needed when you change conda environments, the reference configuration, the web
requirements, or `Dockerfile.base` itself. It downloads and BWA-indexes reference
genomes, so expect 20-30 minutes.

```bash
make docker-build-base                                   # -> vntyper-base:local
make docker-build DOCKER_BASE_IMAGE=vntyper-base:local   # build the app on it
```

In CI this happens automatically: `.github/workflows/docker-base.yml` publishes a base
tagged with a content hash of those inputs, and `docker-build.yml` resolves the same
hash. If you change a base input without publishing a new base, the app build fails
with an explicit "base image does not exist" error rather than silently using a stale one.

### **Pull Pre-built Image**

```bash
# From Docker Hub
docker pull saei/vntyper:latest

# From GitHub Container Registry
docker pull ghcr.io/hassansaei/vntyper:latest
```

### **Generate Apptainer Image**

```bash
apptainer pull docker://saei/vntyper:latest
```

## **Testing the Build**

### **Quick Test** (Recommended)

```bash
# Run all Docker tests with testcontainers (automatic container management)
make docker-test

# Run quick Docker tests (exclude slow adVNTR tests)
make docker-test-quick

# Or use pytest directly
pytest -m docker                    # All Docker tests
pytest -m "docker and not slow"     # Quick tests only
```

### **Test Specific Components**

```bash
# Test BAM pipeline only
pytest tests/docker/test_docker_pipeline.py::test_docker_bam_pipeline -v

# Test adVNTR module (slow test)
pytest tests/docker/test_docker_pipeline.py::test_docker_advntr_pipeline -v

# Test container health
pytest tests/docker/test_docker_pipeline.py::test_docker_container_health -v
```

### **Verify Installation**

```bash
# Check VNtyper version
docker run --rm vntyper:latest vntyper --version

# Check Java runtime
docker run --rm vntyper:latest java -version

# Check bioinformatics tools
docker run --rm vntyper:latest samtools --version
```

## **Configuration**

The CLI needs no configuration. The web service (API + Celery worker + beat) reads
its settings from the environment; `docker/.env.example` lists them, with defaults,
and is the file to copy:

```bash
cp docker/.env.example docker/.env
```

`REDIS_PASSWORD` is **required** and has no default. Redis, the API and the worker
must all be given the same value, so the application refuses to start without it and
`docker compose` refuses to bring the stack up without it. Generate a fresh secret:

```bash
python -c "import secrets; print(secrets.token_urlsafe(32))"
```

When bringing an existing deployment onto this release, generate a new value rather than
carrying one over.

Do not set `CELERY_BROKER_URL` or `CELERY_RESULT_BACKEND`. Celery prefers those over the
broker URL the application passes it, so setting them replaces the URL
`docker/app/celery_app.py` builds from `REDIS_PASSWORD`.

## **Running the Docker Container**

### **CLI Usage**

Run docker interactively:

```bash
   docker run -w /opt/vntyper --rm \
    -v /local/input/folder/:/opt/vntyper/input \
    -v /local/output/folder/:/opt/vntyper/output \
    vntyper:latest \
    vntyper pipeline --bam /local/input/folder/filename.bam \
    -o /local/output/folder/filename/
```
Run apptainer interactively:

```bash
    apptainer run --pwd /opt/vntyper \
    -B /local/input/folder/:/opt/vntyper/input \
    -B /local/output/folder/:/opt/vntyper/output \
    vntyper_latest.sif vntyper pipeline \
    --bam /opt/vntyper/input/filename.bam \
    -o /opt/vntyper/output/filename/ 
```

### **API Usage**

#### **1. Run the API Server**

Start the FastAPI server by running the container:

```bash
docker run -d -p 8000:8000 \
    -v /local/input/folder/:/opt/vntyper/input \
    -v /local/output/folder/:/opt/vntyper/output \
    vntyper:latest
```

#### **2. Submit a Job via API**

Use `curl` to upload the BAM file and pass the required parameters:

```bash
curl -X POST "http://localhost:8000/run-job/" \
    -F "file=@/local/input/folder/filename.bam" \
    -F "thread=4" \
    -F "reference_assembly=hg38" \
    -F "fast_mode=true" \
    -F "keep_intermediates=true" \
    -F "archive_results=true"
```

**Response:**

```json
{
  "message": "Job started",
  "output_dir": "/download/filename"
}
```

#### **3. Download Results**

After the job completes, download the zipped results:

```bash
curl -O "http://localhost:8000/download/filename.zip"
```

## **Notes**

- **Volume Mounts:**
  - Ensure that the local directories `/local/input/folder/` and `/local/output/folder//` exist and have appropriate read/write permissions.
  
- **API Parameters:**
  - `thread`: Number of threads to use (e.g., `4`).
  - `reference_assembly`: Reference genome assembly (e.g., `hg38`).
  - `fast_mode`: Enable fast mode (`true` or `false`).
  - `keep_intermediates`: Retain intermediate files (`true` or `false`).
  - `archive_results`: Archive results as a ZIP file (`true` or `false`).

- **Accessing the API:**
  - **Submit a Job:** Use the `/run-job/` endpoint to submit a BAM file for processing.
  - **Download Results:** Use the `/download/{output_dir}` endpoint to retrieve the processed results.

- **Health Check:**
  - The Docker container includes a health check to ensure it's running correctly. You can monitor the container's health status using:
    ```bash
    docker ps
    ```
  
- **Logs:**
  - To view the container logs for debugging:
    ```bash
    docker logs <container_id>
    ```

## **Testing Architecture**

VNtyper uses [testcontainers-python](https://testcontainers-python.readthedocs.io/) for Docker testing with pytest:

- Automatic container lifecycle management
- Consistent test behavior between local and Docker environments
- Standard pytest fixtures and markers
