# Docker Testing Guide

**VNtyper Docker Testing with Testcontainers**

## Overview

VNtyper uses [testcontainers-python](https://testcontainers-python.readthedocs.io/) for automated Docker integration testing. This ensures the Docker image works correctly in an isolated, reproducible environment.

**Key Benefits:**
- ✅ **100% Test Identity**: Docker tests use the same validation logic as local integration tests
- ✅ **Automatic Lifecycle**: Containers are automatically built, started, and cleaned up
- ✅ **Test Isolation**: Each test gets its own output directory to prevent contamination
- ✅ **CI/CD Ready**: Identical behavior in local development and GitHub Actions

## Prerequisites

### 1. Docker Daemon

Docker must be running on your system:

```bash
# Check if Docker is running
docker info

# Start Docker (if not running)
# macOS/Windows: Start Docker Desktop
# Linux: sudo systemctl start docker
```

### 2. Development Dependencies

Install VNtyper with development dependencies (includes testcontainers):

```bash
pip install -e .[dev]
```

**Key dependency**: `testcontainers>=4.0.0`

## Quick Start

### Run All Docker Tests (Full Suite)

```bash
make test-docker
```

**What it does:**
- Builds Docker image from `docker/Dockerfile`
- Runs 8 parameterized tests (6 BAM pipeline + 2 adVNTR + health checks)
- Takes ~10 minutes (first run with build) or ~5 minutes (cached)

### Run Quick Tests (Recommended for Development)

```bash
make test-docker-quick
```

**What it does:**
- Runs 4 tests: 1 BAM pipeline test + 3 health checks
- Takes ~2 minutes (faster feedback during development)
- Same validation as full suite, just fewer test cases

### Run Specific Test

```bash
# Run single BAM test
pytest tests/docker/test_docker_pipeline.py::test_docker_bam_pipeline[example_b178_hg19_subset_fast] -v

# Run health checks only
pytest tests/docker/test_docker_pipeline.py::test_docker_container_health -v

# Run adVNTR tests (slow)
pytest tests/docker/test_docker_pipeline.py::test_docker_advntr_pipeline -v
```

## Test Architecture

### Directory Structure

```
tests/
├── docker/                      # Docker-specific tests
│   ├── __init__.py
│   ├── conftest.py             # Docker fixtures (testcontainers)
│   └── test_docker_pipeline.py # Docker test functions
├── helpers.py                   # Shared validation functions
├── parametrization.py           # Shared test case definitions
├── test_orchestration.py        # Shared test orchestration
└── test_data_config.json        # Expected values and tolerances
```

### How It Works

1. **Session-Scoped Image Build** (`tests/docker/conftest.py:vntyper_image`)
   - Builds Docker image once per pytest session
   - Command: `docker build -f docker/Dockerfile -t vntyper:test .`
   - Image is reused across all tests (efficient)

2. **Module-Scoped Container** (`tests/docker/conftest.py:vntyper_container`)
   - Creates running container with volume mounts
   - Mounts `tests/data/` → `/opt/vntyper/input` (read-only)
   - Mounts temp directory → `/opt/vntyper/output` (read-write)
   - Container runs as non-root user (UID 1001)

3. **Test-Specific Subdirectories**
   - Each test creates isolated subdirectory: `output_dir/test_name/`
   - Prevents output contamination between tests
   - Permissions set to 777 for container user write access

4. **Shared Validation**
   - Docker tests call `run_bam_test_case()` from `tests/test_orchestration.py`
   - Uses same validation functions as local tests (`tests/helpers.py`)
   - Guarantees 100% test identity

## Available Tests

### BAM Pipeline Tests (6 tests)

```bash
pytest -m docker -k "bam_pipeline" -v
```

Tests various BAM files with different reference assemblies:
- `example_b178_hg19_subset_fast` - Fast test case (used in quick suite)
- `example_a5c1_hg19` - Standard hg19 test
- `example_66bf_hg38` - Standard hg38 test
- And 3 more variants

### adVNTR Module Tests (2 tests, slow)

```bash
pytest -m docker -m slow -v
```

Tests adVNTR integration (takes ~9 minutes per test):
- Requires `--extra-modules advntr` flag
- Validates adVNTR database extraction and genotyping

### Health Check Tests (3 tests)

```bash
pytest tests/docker/test_docker_pipeline.py -k "health or volume or dependencies" -v
```

Validates Docker environment:
- `test_docker_container_health` - vntyper command availability
- `test_docker_volume_mounts` - Volume mount accessibility
- `test_docker_dependencies` - Required tools (Java, samtools, BWA)

## Configuration

### Test Cases

Test cases are defined in `tests/test_data_config.json`:

```json
{
  "test_name": "example_b178_hg19_subset_fast",
  "bam_file": "tests/data/example_b178_hg19_subset.bam",
  "reference_assembly": "hg19",
  "expected_values": {
    "Kestrel_Inference": 416,
    "Estimated_Depth_AlternateVariant": {"value": 176.34, "tolerance": 0.05}
  }
}
```

**Tolerance Support:**
- Exact values: `"Kestrel_Inference": 416` (must match exactly)
- Tolerances: `{"value": 176.34, "tolerance": 0.05}` (±5% allowed)

### Pytest Markers

Defined in `pytest.ini`:

```ini
markers =
    docker: marks tests as Docker container tests (requires Docker daemon)
    slow: marks tests as slow-running tests (e.g., adVNTR module)
```

**Usage:**
```bash
# Run only Docker tests
pytest -m docker

# Exclude slow tests
pytest -m docker -m "not slow"

# Run slow tests only
pytest -m docker -m slow
```

## Troubleshooting

### Docker Daemon Not Running

```
Error: Cannot connect to Docker daemon
Solution: Start Docker Desktop (macOS/Windows) or docker service (Linux)
```

### Permission Denied on Volume Mounts

```
Error: Permission denied: '/opt/vntyper/output/...'
Solution: Tests automatically set chmod 777 on output directories.
          This is safe for isolated pytest temp directories.
```

### Testcontainers Not Installed

```
Error: ModuleNotFoundError: No module named 'testcontainers'
Solution: pip install -e .[dev]
```

### Image Build Fails

```bash
# Manually test Docker build
docker build -f docker/Dockerfile -t vntyper:test .

# Check build logs
make docker-build
```

### Port Conflicts

Testcontainers uses random ports - no manual port configuration needed.

### Cleanup Stuck Containers

```bash
# List VNtyper test containers
docker ps -a | grep vntyper

# Remove all VNtyper containers
docker rm -f $(docker ps -a -q --filter ancestor=vntyper:test)

# Remove VNtyper images
make docker-clean
```

## CI/CD Integration

Docker tests run automatically in GitHub Actions (`.github/workflows/docker-build.yml`):

**On Pull Requests:**
- Builds Docker image with layer caching (60-80% faster)
- Runs `make test-docker-quick` (4 tests, 15min timeout)

**On Main Branch:**
- Runs full `make test-docker` suite (8 tests, 60min timeout)
- Uploads test results as artifacts (30-day retention)

**Manual Trigger:**
- Use "Run workflow" button in GitHub Actions UI
- Runs full test suite on any branch

## Advanced Usage

### Running Tests in Parallel

```bash
# Not recommended for Docker tests (container conflicts)
# Use quick suite instead for faster feedback
make test-docker-quick
```

### Custom Output Directory

Tests use pytest's `tmp_path_factory` - output is automatically cleaned up.

### Debugging Container

```bash
# Start container interactively
docker run -it --entrypoint /bin/bash vntyper:test

# Inspect running test container (during test pause)
docker ps
docker exec -it <container_id> /bin/bash
```

### View Container Logs

```bash
# Enable verbose pytest output
pytest tests/docker/test_docker_pipeline.py -v -s

# Check pytest log file
cat tests/pytest.log
```

## Best Practices

1. **Use Quick Suite for Development**: `make test-docker-quick` provides fast feedback
2. **Run Full Suite Before PR**: Ensure all test cases pass
3. **Don't Modify test_data_config.json Without Review**: Affects local and Docker tests
4. **Keep Docker Image Lean**: Multi-stage builds reduce size and improve speed
5. **Let Testcontainers Manage Lifecycle**: Don't manually start/stop containers

## Related Documentation

- **Docker Image**: `docker/README.md` - Building and running VNtyper Docker image
- **Dockerfile**: `docker/Dockerfile` - Multi-stage build configuration
- **Local Testing**: Project README - Running local integration tests
- **Test Data**: `tests/test_data_config.json` - Test case definitions

## Support

If Docker tests fail but local integration tests pass:
1. Check Docker daemon is running
2. Verify `pip install -e .[dev]` completed successfully
3. Try rebuilding image: `docker build -f docker/Dockerfile -t vntyper:test .`
4. Check logs: `tests/pytest.log`

For issues, report at: https://github.com/hassansaei/VNtyper/issues
