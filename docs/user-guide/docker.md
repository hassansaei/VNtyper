# Docker

VNtyper 2 provides pre-built Docker images with all dependencies included.

## Pull Pre-built Images

=== "Docker Hub"

    ```bash
    docker pull saei/vntyper:latest
    ```

=== "GitHub Container Registry"

    ```bash
    docker pull ghcr.io/hassansaei/vntyper:latest
    ```

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
    vntyper:latest \
    vntyper pipeline --bam /opt/vntyper/input/sample.bam \
    -o /opt/vntyper/output/sample/
```

!!! tip "File Permissions"
    Use `--user $(id -u):$(id -g)` so output files are owned by your host user, not root.

### Volume Mounts

| Host Path | Container Path | Purpose |
|-----------|----------------|---------|
| Input directory | `/opt/vntyper/input` | BAM/FASTQ files |
| Output directory | `/opt/vntyper/output` | Pipeline results |

## Verify Installation

```bash
docker run --rm vntyper:latest vntyper --version
docker run --rm vntyper:latest samtools --version
docker run --rm vntyper:latest java -version
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

Start the FastAPI server for programmatic access:

```bash
docker run -d -p 8000:8000 \
    -v /path/to/input:/opt/vntyper/input \
    -v /path/to/output:/opt/vntyper/output \
    vntyper:latest
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
apptainer pull docker://saei/vntyper:latest
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
