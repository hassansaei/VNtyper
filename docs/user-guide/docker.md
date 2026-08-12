# Docker

VNtyper 2 provides pre-built Docker images with all dependencies included.

## Pull Pre-built Images

Pull the latest released image from GitHub Container Registry:

```bash
docker pull ghcr.io/hassansaei/vntyper:latest
```

Released images are published to GHCR as `latest` (the newest release), the immutable
`vX.Y.Z` and `X.Y.Z` tags naming one exact release, and the moving `X` and `X.Y` series
tags. Pin `vX.Y.Z` for a reproducible run. The `main` tag is rolling and tracks the
default branch. Docker Hub artifacts are legacy, frozen, and unsupported.

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
    VNtyper never writes into the directory holding the patient alignment. Mounting one
    host directory at both container paths puts the output root inside the input tree,
    and the run is rejected:

    ```
    Alignment output root must stay outside the patient input tree:
    /opt/vntyper/output/sample lies under /opt/vntyper/output, which is the same
    directory as the patient input tree /opt/vntyper/input.
    Give the run separate input and output directories.
    ```

    This fires even though the two container paths look different, because the check
    compares the directories themselves rather than their names. Do **not** do this:

    ```bash
    # rejected: one host directory mounted at both paths
    docker run -w /opt/vntyper --rm \
        --user $(id -u):$(id -g) \
        -v "$PWD":/opt/vntyper/input \
        -v "$PWD":/opt/vntyper/output \
        ghcr.io/hassansaei/vntyper:latest \
        vntyper pipeline --cram /opt/vntyper/input/sample.cram \
        -o /opt/vntyper/output/sample/
    ```

    Mount separate directories instead. A `results/` subdirectory beside the data is
    fine, as long as it is not the directory the alignment itself sits in:

    ```bash
    docker run -w /opt/vntyper --rm \
        --user $(id -u):$(id -g) \
        -v "$PWD":/opt/vntyper/input \
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

The container includes a built-in health check. Monitor status with:

```bash
docker ps
```

View container logs:

```bash
docker logs <container_id>
```

## API Server

The API stores job state in Redis and hands work to a Celery worker, so it needs a
Redis instance and the password for it. `REDIS_PASSWORD` is required and has no
default: the API, the worker and Redis itself must all be given the same value, and
the service refuses to start without it. Generate a fresh secret with
`python -c "import secrets; print(secrets.token_urlsafe(32))"`.

The full stack (Redis, API, worker, beat) is easiest to bring up with Compose. Copy
the environment template, set the password in it, and start:

```bash
cp docker/.env.example docker/.env
$EDITOR docker/.env          # set REDIS_PASSWORD
docker compose -f docker/docker-compose.yml up -d
```

`docker/.env.example` documents the remaining settings. To run the API container on
its own against an existing Redis:

```bash
docker run -d -p 8000:8000 \
    -e REDIS_HOST=my-redis-host \
    -e REDIS_PASSWORD="$REDIS_PASSWORD" \
    -v /path/to/input:/opt/vntyper/input \
    -v /path/to/output:/opt/vntyper/output \
    ghcr.io/hassansaei/vntyper:latest
```

Submit a job:

```bash
curl -X POST "http://localhost:8000/run-job/" \
    -F "file=@sample.bam" \
    -F "thread=4" \
    -F "reference_assembly=hg38" \
    -F "fast_mode=true" \
    -F "archive_results=true"
```

Download results after completion:

```bash
curl -O "http://localhost:8000/download/sample.zip"
```

## Apptainer / Singularity

Convert the Docker image to an Apptainer SIF:

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
