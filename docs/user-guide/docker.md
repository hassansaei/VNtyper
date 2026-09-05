# Docker

VNtyper 2 provides pre-built Docker images with all dependencies included.

## Pull Pre-built Images

Pull the latest released image from GitHub Container Registry:

```bash
docker pull ghcr.io/hassansaei/vntyper:latest
```

Released images are published to GHCR as `latest` (the newest release), the immutable `vX.Y.Z` and `X.Y.Z` tags naming one exact release, and the moving `X` and `X.Y` series tags. Pin `vX.Y.Z` for a reproducible run. The `main` tag is rolling and tracks the default branch. Docker Hub artifacts are legacy, frozen, and unsupported.

## Build from Source

```bash
git clone https://github.com/hassansaei/VNtyper.git
cd VNtyper
DOCKER_BUILDKIT=1 docker build -f docker/Dockerfile -t vntyper:latest .
```

Or use Make:

```bash
make docker-build
```

## Run the Pipeline

```bash
docker run -w /opt/vntyper --rm \
    --user $(id -u):$(id -g) \
    -v /path/to/input:/opt/vntyper/input \
    -v /path/to/output:/opt/vntyper/output \
    ghcr.io/hassansaei/vntyper:latest \
    vntyper pipeline --bam /opt/vntyper/input/sample.bam \
    -o /opt/vntyper/output/sample/
```

!!! tip "File Permissions"
    Use `--user $(id -u):$(id -g)` so output files are owned by your host user, not root.

### Volume Mounts

| Host Path | Container Path | Purpose |
|-----------|----------------|---------|
| Input directory | `/opt/vntyper/input` | BAM/CRAM/FASTQ files |
| Output directory | `/opt/vntyper/output` | Pipeline results |

!!! warning "Input and output must be two different host directories"
    VNtyper never writes into the directory holding the patient alignment. Mounting one host directory at both container paths places the output root inside the input tree, causing pipeline termination:

    ```
    Alignment output root must stay outside the patient input tree:
    /opt/vntyper/output/sample lies under /opt/vntyper/output, which is the same
    directory as the patient input tree /opt/vntyper/input.
    Give the run separate input and output directories.
    ```

    Mount separate host directories:

    ```bash
    docker run -w /opt/vntyper --rm \
        --user $(id -u):$(id -g) \
        -v "$PWD/data":/opt/vntyper/input \
        -v "$PWD/results":/opt/vntyper/output \
        ghcr.io/hassansaei/vntyper:latest \
        vntyper pipeline --cram /opt/vntyper/input/sample.cram \
        -o /opt/vntyper/output/sample/
    ```

## Verify Installation

```bash
docker run --rm ghcr.io/hassansaei/vntyper:latest vntyper --version
docker run --rm ghcr.io/hassansaei/vntyper:latest samtools --version
docker run --rm ghcr.io/hassansaei/vntyper:latest java -version
```

## Health Checks

Inspect container status:

```bash
docker ps
```

View container logs:

```bash
docker logs <container_id>
```

## API Server

The API stores job state in Redis and dispatches tasks to Celery workers. Setting `REDIS_PASSWORD` is mandatory across API, worker, and Redis instances; the service refuses startup without it. Generate a secret using:

```bash
python -c "import secrets; print(secrets.token_urlsafe(32))"
```

Initialize the full stack (Redis, API, worker, beat) via Docker Compose:

```bash
cp docker/.env.example docker/.env
$EDITOR docker/.env
docker compose -f docker/docker-compose.yml up -d
```

Shipped Compose keeps ordinary-job outputs in the persistent, service-private `result_store`, mounted only into the API and worker. The protected handoff spool uses the same service boundary. Replacing `/opt/vntyper/output` with a host/shared bind-mount override may be useful for a custom deployment, but it gives that mount's actors access to the result namespace and weakens the shipped security boundary.

For an existing web deployment, use this migration sequence:

1. Pause new submissions. Drain both the regular and long queues and all active jobs to completion; verify both queues and the active job count are zero. Never purge queued messages, because doing so can strand their protected-spool uploads.
2. Stop the API, workers, and beat, then provision the named `result_store`.
3. While services remain stopped and detached from the legacy host bind, either (recommended) copy existing unexpired output into `result_store` and retain a backup, noting that legacy bytes cannot be retroactively integrity-attested; or explicitly accept retirement and unavailability of those results and archive or remove the legacy store.
4. Deploy the API, all workers, and beat; verify retained result access; then resume submissions.

Arbitrary same-UID code in either the API or worker service namespace is out of scope: such code can access the private volumes or worker descriptors directly.

To execute the API container independently against an existing Redis instance:

```bash
docker run -d -p 8000:8000 \
    -e REDIS_HOST=my-redis-host \
    -e REDIS_PASSWORD="$REDIS_PASSWORD" \
    -v /path/to/input:/opt/vntyper/input \
    -v vntyper_handoff_spool:/opt/vntyper/handoff \
    -v vntyper_result_store:/opt/vntyper/output \
    ghcr.io/hassansaei/vntyper:latest
```

Submit an analysis job:

```bash
curl -X POST "http://localhost:8000/run-job/" \
    -F "file=@sample.bam" \
    -F "thread=4" \
    -F "reference_assembly=hg38" \
    -F "fast_mode=true" \
    -F "archive_results=true"
```

Download completed results:

```bash
curl -O "http://localhost:8000/download/sample.zip"
```

## Apptainer / Singularity

Convert the container image to an Apptainer SIF:

```bash
apptainer pull docker://ghcr.io/hassansaei/vntyper:latest
```

Run with Apptainer:

```bash
apptainer run --pwd /opt/vntyper \
    -B /path/to/input:/opt/vntyper/input \
    -B /path/to/output:/opt/vntyper/output \
    vntyper_latest.sif vntyper pipeline \
    --bam /opt/vntyper/input/sample.bam \
    -o /opt/vntyper/output/sample/
```
