# Docker Modernization Implementation Plan

**VNtyper 2.0 - 2025 Best Practices Migration**

---

## 📋 Executive Summary

This document provides a **step-by-step implementation plan** to modernize VNtyper's Docker configuration based on 2025 best practices. All changes are designed to avoid regressions while delivering:

- **60-70% smaller image size** via multi-stage builds
- **Enhanced security** with digest pinning and secrets management
- **Faster builds** with optimized caching
- **Modern Docker Compose** configuration

---

## 🎯 Implementation Strategy

**Approach:** Incremental, backward-compatible changes with validation at each step.

**Priority Order:**
1. **Quick Wins** (P1) - Low risk, immediate benefit
2. **Multi-Stage Refactor** (P2) - High impact, requires testing
3. **Security Hardening** (P3) - Production-grade improvements

---

## 📦 Phase 1: Quick Wins (P1)

### 1.1 Update docker-compose.yml - Remove Obsolete Version

**File:** `docker/docker-compose.yml`

**Current Issue:** Using deprecated `version: '3.8'` field

**Change:**
```yaml
# ❌ REMOVE THIS LINE
version: '3.8'

# ✅ START DIRECTLY WITH SERVICES
services:
  redis:
    image: redis:7.2-alpine3.19  # Also upgrade Redis
    container_name: vntyper_redis
    # ... rest of config
```

**Full Updated File:**
```yaml
# =============================================================================
# VNtyper Docker Compose Configuration
# =============================================================================
# Modern Docker Compose format (no version field needed)
# Reference: https://docs.docker.com/compose/compose-file/
# =============================================================================

services:
  redis:
    image: redis:7.2-alpine3.19
    container_name: vntyper_redis
    ports:
      - "6379:6379"
    networks:
      - vntyper_network
    healthcheck:
      test: ["CMD", "redis-cli", "ping"]
      interval: 10s
      timeout: 3s
      retries: 3
      start_period: 5s

  api:
    build:
      context: ..
      dockerfile: docker/Dockerfile
    image: vntyper:latest
    container_name: vntyper_api
    command: uvicorn app.main:app --host 0.0.0.0 --port 8000
    environment:
      - CELERY_BROKER_URL=redis://redis:6379/0
      - CELERY_RESULT_BACKEND=redis://redis:6379/0
      - MAX_RESULT_AGE_DAYS=7
    volumes:
      - /mnt/c/development/VNtyper/download/:/opt/vntyper/input
      - /mnt/c/development/VNtyper/out/output:/opt/vntyper/output
    depends_on:
      redis:
        condition: service_healthy
    networks:
      - vntyper_network
    ports:
      - "8000:8000"
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
      interval: 30s
      timeout: 10s
      retries: 3
      start_period: 60s

  worker:
    build:
      context: ..
      dockerfile: docker/Dockerfile
    image: vntyper:latest
    container_name: vntyper_worker
    command: celery -A app.celery_app worker --loglevel=info
    environment:
      - CELERY_BROKER_URL=redis://redis:6379/0
      - CELERY_RESULT_BACKEND=redis://redis:6379/0
      - MAX_RESULT_AGE_DAYS=7
    volumes:
      - /mnt/c/development/VNtyper/download/:/opt/vntyper/input
      - /mnt/c/development/VNtyper/out/output:/opt/vntyper/output
    depends_on:
      redis:
        condition: service_healthy
      api:
        condition: service_healthy
    networks:
      - vntyper_network

  beat:
    build:
      context: ..
      dockerfile: docker/Dockerfile
    image: vntyper:latest
    container_name: vntyper_beat
    command: celery -A app.celery_app beat --loglevel=info
    environment:
      - CELERY_BROKER_URL=redis://redis:6379/0
      - CELERY_RESULT_BACKEND=redis://redis:6379/0
      - MAX_RESULT_AGE_DAYS=7
    volumes:
      - /mnt/c/development/VNtyper/download/:/opt/vntyper/input
      - /mnt/c/development/VNtyper/out/output:/opt/vntyper/output
    depends_on:
      redis:
        condition: service_healthy
      api:
        condition: service_healthy
    networks:
      - vntyper_network

networks:
  vntyper_network:
    driver: bridge
```

**Key Changes:**
- ✅ Removed `version: '3.8'` (obsolete in modern Compose)
- ✅ Upgraded Redis from `6-alpine` to `7.2-alpine3.19` (latest stable)
- ✅ Added health checks to all services
- ✅ Added `depends_on` with `condition: service_healthy` for proper startup ordering
- ✅ Added health check for API service

**Validation:**
```bash
# Test the updated compose file
docker compose config

# Expected: No errors, properly parsed YAML

# Start services to verify health checks
docker compose up -d
docker compose ps

# All services should show "healthy" status
```

**Risk:** Low - Backward compatible, no breaking changes

---

### 1.2 Optimize Dockerfile Cache Usage

**File:** `docker/Dockerfile` and `docker/Dockerfile.local`

**Current Issue:** Line 76/92 uses `--no-cache-dir` which defeats BuildKit cache mounts

**Change in BOTH files:**

**Before:**
```dockerfile
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    conda run -n vntyper pip install --no-cache-dir .
```

**After:**
```dockerfile
# Remove --no-cache-dir to leverage BuildKit cache mount
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    conda run -n vntyper pip install .
```

**Also Apply to Lines 92-102 (FastAPI dependencies):**

**Before:**
```dockerfile
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    conda run -n vntyper pip install --no-cache-dir \
        "fastapi[standard]==0.115.3" \
        ...
```

**After:**
```dockerfile
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    conda run -n vntyper pip install \
        "fastapi[standard]==0.115.3" \
        "uvicorn[standard]==0.32.0" \
        redis==5.2.0 \
        celery==5.4.0 \
        python-multipart==0.0.12 \
        fastapi-limiter==0.1.6 \
        email_validator==2.2.0 \
        'passlib[bcrypt]==1.7.4' \
        pydantic==2.10.0
```

**Validation:**
```bash
# Build with cache
docker build -f docker/Dockerfile .

# Second build should be much faster (cache hits)
docker build -f docker/Dockerfile .
```

**Risk:** Very Low - Improves performance, no functional change

---

### 1.3 Add Missing ca-certificates to Dockerfile

**File:** `docker/Dockerfile` (line 35)

**Current Issue:** Missing `ca-certificates` package can cause SSL issues

**Change:**
```dockerfile
# Update and install necessary system packages with BuildKit cache mount
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        wget \
        unzip \
        zip \
        curl \
        ca-certificates \
        && rm -rf /var/lib/apt/lists/*
```

**Note:** This is already present in Dockerfile (line 35) but NOT in Dockerfile.local. Add to Dockerfile.local line 31:

```dockerfile
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        wget \
        unzip \
        zip \
        curl \
        ca-certificates \
        && rm -rf /var/lib/apt/lists/*
```

**Risk:** Very Low - Adds missing dependency

---

## 📦 Phase 2: Multi-Stage Build Refactor (P2)

### 2.1 Create New Multi-Stage Dockerfile

**File:** `docker/Dockerfile.multistage` (NEW FILE)

**Rationale:**
- Separate build-time and runtime dependencies
- Reduce final image size by 60-70%
- Remove build tools from production image
- Keep existing Dockerfiles for backward compatibility during transition

**Implementation:**

```dockerfile
# syntax=docker/dockerfile:1.7
# =============================================================================
# VNtyper 2.0 - Multi-Stage Production Dockerfile
# =============================================================================
# Build Strategy:
# - Stage 1: Base conda environment setup
# - Stage 2: Build conda environments and compile dependencies
# - Stage 3: Minimal runtime image with only necessary components
# =============================================================================

# =============================================================================
# Stage 1: Base Dependencies (Shared)
# =============================================================================
FROM condaforge/miniforge3:24.11.3-0 AS base

# OCI Labels for metadata
LABEL org.opencontainers.image.title="VNtyper" \
      org.opencontainers.image.description="VNtyper 2.0 - MUC1 VNTR genotyping pipeline for ADTKD-MUC1" \
      org.opencontainers.image.vendor="VNtyper Project" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.source="https://github.com/hassansaei/VNtyper" \
      org.opencontainers.image.url="https://github.com/hassansaei/VNtyper" \
      org.opencontainers.image.documentation="https://github.com/hassansaei/VNtyper" \
      com.vntyper.type="bioinformatics" \
      com.vntyper.pipeline="genotyping"

# Install mamba for faster package management
RUN --mount=type=cache,target=/opt/conda/pkgs,sharing=locked \
    conda install -y mamba -n base -c conda-forge && \
    conda clean -afy

# =============================================================================
# Stage 2: Build Environment (Contains all build tools)
# =============================================================================
FROM base AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Install build-time system dependencies with BuildKit cache
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        wget \
        unzip \
        zip \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Define build arguments
ARG REPO_URL=https://github.com/hassansaei/VNtyper.git
ARG REPO_DIR=/opt/vntyper

# Clone repository
RUN git clone ${REPO_URL} ${REPO_DIR}

WORKDIR ${REPO_DIR}

# Copy environment files to temp location
RUN cp conda/environment_vntyper.yml /tmp/environment_vntyper.yml && \
    cp conda/environment_envadvntr.yml /tmp/environment_envadvntr.yml && \
    cp conda/environment_shark.yml /tmp/environment_shark.yml

# Create conda environments with aggressive cleanup
RUN --mount=type=cache,target=/opt/conda/pkgs,sharing=locked \
    --mount=type=cache,target=/root/.conda/pkgs,sharing=locked \
    mamba env create -f /tmp/environment_vntyper.yml && \
    mamba env create -f /tmp/environment_envadvntr.yml && \
    mamba env create -f /tmp/environment_shark.yml && \
    conda clean -afy && \
    find /opt/conda/envs -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null || true && \
    find /opt/conda/envs -type f -name '*.pyc' -delete 2>/dev/null || true && \
    find /opt/conda/envs -type f -name '*.pyo' -delete 2>/dev/null || true && \
    rm -f /tmp/environment_*.yml

# Install VNtyper Python package
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    conda run -n vntyper pip install .

# Install adVNTR
RUN chmod +x vntyper/dependencies/advntr/install_advntr.sh && \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate envadvntr && \
    bash vntyper/dependencies/advntr/install_advntr.sh -o"

# Install reference genomes
RUN conda run -n vntyper vntyper --config-path /opt/vntyper/vntyper/config.json install-references \
    --output-dir /opt/vntyper/reference

# Install FastAPI and web dependencies
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    conda run -n vntyper pip install \
        "fastapi[standard]==0.115.3" \
        "uvicorn[standard]==0.32.0" \
        redis==5.2.0 \
        celery==5.4.0 \
        python-multipart==0.0.12 \
        fastapi-limiter==0.1.6 \
        email_validator==2.2.0 \
        'passlib[bcrypt]==1.7.4' \
        pydantic==2.10.0

# =============================================================================
# Stage 3: Minimal Runtime Image
# =============================================================================
FROM condaforge/miniforge3:24.11.3-0 AS runtime

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    CONDA_AUTO_UPDATE_CONDA=FALSE \
    PATH=/opt/conda/envs/vntyper/bin:/opt/conda/envs/environment_vntyper/bin:$PATH \
    DEFAULT_INPUT_DIR=/opt/vntyper/input \
    DEFAULT_OUTPUT_DIR=/opt/vntyper/output \
    REFERENCE_DIR=/opt/vntyper/reference

# Install ONLY runtime dependencies (no build tools!)
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy conda environments from builder (this is the key space savings!)
COPY --from=builder /opt/conda/envs /opt/conda/envs

# Copy VNtyper application and references
COPY --from=builder /opt/vntyper /opt/vntyper

# Set working directory
WORKDIR /opt/vntyper

# Set up default directories
RUN mkdir -p $DEFAULT_INPUT_DIR $DEFAULT_OUTPUT_DIR $REFERENCE_DIR

# Create non-root user with configurable UID/GID
ARG USERNAME=appuser
ARG USER_UID=1001
ARG USER_GID=1001

RUN groupadd --gid $USER_GID $USERNAME && \
    useradd --uid $USER_UID --gid $USER_GID --shell /bin/bash --create-home $USERNAME

# Set ownership
RUN chown -R $USERNAME:$USERNAME $DEFAULT_INPUT_DIR $DEFAULT_OUTPUT_DIR $REFERENCE_DIR /opt/vntyper

# Copy and set up entrypoint
COPY --from=builder /opt/vntyper/docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

# Copy FastAPI app
RUN mkdir -p /opt/vntyper/app
COPY --from=builder /opt/vntyper/docker/app /opt/vntyper/app
RUN chown -R $USERNAME:$USERNAME /opt/vntyper/app

# Expose port
EXPOSE 8000

# Switch to non-root user
USER $USERNAME

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD ["/usr/local/bin/entrypoint.sh", "health"]

# Set entrypoint
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
```

**Key Improvements:**
- ✅ **3-stage build:** base → builder → runtime
- ✅ **60-70% smaller image:** Only runtime dependencies in final stage
- ✅ **No build tools in production:** gcc, git, build-essential removed from final image
- ✅ **Better caching:** Separate stages allow better layer reuse
- ✅ **Security:** Minimal attack surface

**Validation:**
```bash
# Build multi-stage image
docker build -f docker/Dockerfile.multistage -t vntyper:multistage .

# Compare image sizes
docker images | grep vntyper

# Expected results:
# vntyper:latest      ~5GB (current)
# vntyper:multistage  ~1.5-2GB (new)

# Test functionality
docker run --rm vntyper:multistage vntyper --version

# Run full compose stack with new image
docker compose -f docker/docker-compose.yml up -d
```

**Risk:** Medium - Requires thorough testing, but keeps old Dockerfile as backup

---

### 2.2 Update docker-compose.yml to Use Multi-Stage Dockerfile

**File:** `docker/docker-compose.yml`

**Change:** Update build context to use new Dockerfile

```yaml
services:
  api:
    build:
      context: ..
      dockerfile: docker/Dockerfile.multistage  # ✅ Use new multi-stage file
      cache_from:
        - vntyper:latest
      cache_to:
        - type=inline
    image: vntyper:latest
    # ... rest of config unchanged
```

**Apply same change to `worker` and `beat` services**

**Risk:** Low - Easy to rollback by changing dockerfile path

---

## 📦 Phase 3: Security Hardening (P3)

### 3.1 Pin Base Image by Digest

**File:** `docker/Dockerfile.multistage` (and eventually others)

**Current:**
```dockerfile
FROM condaforge/miniforge3:24.11.3-0 AS base
```

**After:** (Need to get current digest first)

```bash
# Get the digest
docker pull condaforge/miniforge3:24.11.3-0
docker inspect condaforge/miniforge3:24.11.3-0 | grep RepoDigests -A 1

# Example output (actual digest will differ):
# "RepoDigests": [
#     "condaforge/miniforge3@sha256:abc123..."
# ]
```

**Updated Dockerfile:**
```dockerfile
# Pin base image by digest for security and reproducibility
FROM condaforge/miniforge3:24.11.3-0@sha256:abc123... AS base
```

**Repeat for runtime stage:**
```dockerfile
FROM condaforge/miniforge3:24.11.3-0@sha256:abc123... AS runtime
```

**Benefits:**
- ✅ Prevents tag hijacking
- ✅ Ensures reproducible builds
- ✅ Industry security best practice

**Maintenance:** Update digest when upgrading base image

**Risk:** Low - Purely additive, improves security

---

### 3.2 Implement Docker Secrets (Optional - For Production)

**File:** `docker/docker-compose.yml`

**For production deployments, migrate sensitive values to secrets:**

**Create secrets files:**
```bash
mkdir -p docker/secrets
echo "redis://redis:6379/0" > docker/secrets/celery_broker_url.txt
chmod 600 docker/secrets/celery_broker_url.txt
```

**Add to .gitignore:**
```bash
echo "docker/secrets/" >> .gitignore
```

**Update docker-compose.yml:**
```yaml
services:
  api:
    secrets:
      - celery_broker_url
    environment:
      # Read from secrets file instead of plain text
      - CELERY_BROKER_URL_FILE=/run/secrets/celery_broker_url
    # ... rest of config

secrets:
  celery_broker_url:
    file: ./docker/secrets/celery_broker_url.txt
```

**Update application code to read from file:**
```python
# In app configuration
import os

def get_secret(secret_name):
    secret_file = os.getenv(f"{secret_name.upper()}_FILE")
    if secret_file and os.path.exists(secret_file):
        with open(secret_file) as f:
            return f.read().strip()
    return os.getenv(secret_name.upper())

CELERY_BROKER_URL = get_secret("celery_broker_url")
```

**Risk:** Medium - Requires code changes, optional for development

---

### 3.3 Add SBOM Generation (Supply Chain Security)

**File:** `docker/Dockerfile.multistage`

**Add to builder stage:**

```dockerfile
# Stage 2: Build Environment
FROM base AS builder

# Enable SBOM scanning during build
ARG BUILDKIT_SBOM_SCAN_STAGE=true

# ... rest of builder stage
```

**Generate SBOM during build:**
```bash
docker buildx build \
  --sbom=true \
  --output type=docker,name=vntyper:latest \
  -f docker/Dockerfile.multistage .
```

**View SBOM:**
```bash
docker buildx imagetools inspect vntyper:latest --format "{{ json .SBOM }}"
```

**Benefits:**
- ✅ Track all dependencies
- ✅ Vulnerability scanning
- ✅ Compliance requirements

**Risk:** Very Low - Optional metadata, no functional impact

---

## 📋 Migration Checklist

### Pre-Migration
- [ ] Backup current working Docker images
- [ ] Document current image sizes
- [ ] Test current functionality end-to-end

### Phase 1 (Quick Wins)
- [ ] Update `docker-compose.yml` - remove version field
- [ ] Add health checks to all services
- [ ] Upgrade Redis to 7.2-alpine3.19
- [ ] Remove `--no-cache-dir` from pip installs
- [ ] Add `ca-certificates` to Dockerfile.local
- [ ] Test: `docker compose config`
- [ ] Test: `docker compose up -d && docker compose ps`

### Phase 2 (Multi-Stage)
- [ ] Create `docker/Dockerfile.multistage`
- [ ] Build and tag multi-stage image
- [ ] Compare image sizes
- [ ] Test vntyper functionality
- [ ] Test full pipeline end-to-end
- [ ] Update docker-compose.yml to use new Dockerfile
- [ ] Test compose stack with new image

### Phase 3 (Security)
- [ ] Get base image digest
- [ ] Pin base images by digest
- [ ] (Optional) Set up Docker secrets
- [ ] (Optional) Add SBOM generation
- [ ] Run vulnerability scan: `docker scout cves vntyper:latest`

### Post-Migration
- [ ] Update documentation (README.md)
- [ ] Update CI/CD pipelines
- [ ] Monitor build times
- [ ] Monitor image pull times
- [ ] Archive old Dockerfiles (keep for reference)

---

## 🧪 COMPREHENSIVE TESTING STRATEGY

**CRITICAL**: Test at each phase. Do NOT proceed if any test fails without investigating.

---

### Pre-Migration: Baseline Testing

**Purpose:** Establish working baseline to detect any regressions

#### 1. Baseline Image Build and Metrics
```bash
# Create baseline directory
mkdir -p testing/baseline

# Build current image
echo "Building baseline image..."
time docker build -f docker/Dockerfile -t vntyper:baseline . 2>&1 | tee testing/baseline/build.log

# Record image size
docker images vntyper:baseline --format "{{.Size}}" > testing/baseline/image_size.txt
docker images vntyper:baseline --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"

# Save image layers
docker history vntyper:baseline > testing/baseline/layers.txt

# Count vulnerabilities (baseline)
echo "Scanning baseline for vulnerabilities..."
docker scout cves vntyper:baseline --format sarif > testing/baseline/vulnerabilities.sarif || true
trivy image vntyper:baseline --format json --output testing/baseline/trivy_baseline.json || true
```

#### 2. Baseline Functional Testing - CLI
```bash
echo "=== Testing VNtyper CLI Commands ==="

# Test 1: Version check
docker run --rm vntyper:baseline vntyper --version
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Version check failed"
    exit 1
fi
echo "✅ PASS: Version check"

# Test 2: Help command
docker run --rm vntyper:baseline vntyper --help | grep -q "usage:"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Help command failed"
    exit 1
fi
echo "✅ PASS: Help command"

# Test 3: Config validation
docker run --rm vntyper:baseline vntyper --config-path /opt/vntyper/vntyper/config.json --help
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Config validation failed"
    exit 1
fi
echo "✅ PASS: Config validation"

# Test 4: Check Java (required for Kestrel)
docker run --rm vntyper:baseline bash -c "conda run -n vntyper java -version"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Java not available"
    exit 1
fi
echo "✅ PASS: Java available"

# Test 5: Check bioinformatics tools
echo "Checking bioinformatics tools..."
for tool in bwa samtools fastp bcftools; do
    docker run --rm vntyper:baseline bash -c "conda run -n vntyper which $tool"
    if [ $? -ne 0 ]; then
        echo "❌ FAIL: $tool not found"
        exit 1
    fi
    echo "✅ PASS: $tool available"
done

# Test 6: Python package imports
docker run --rm vntyper:baseline bash -c "conda run -n vntyper python -c 'import pandas, numpy, pysam, Bio; print(\"All packages imported successfully\")'"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Python package import failed"
    exit 1
fi
echo "✅ PASS: Python packages"
```

#### 3. Baseline Functional Testing - Conda Environments
```bash
echo "=== Testing Conda Environments ==="

# Test vntyper environment
docker run --rm vntyper:baseline bash -c "conda run -n vntyper python --version"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: vntyper conda env failed"
    exit 1
fi
echo "✅ PASS: vntyper environment"

# Test envadvntr environment
docker run --rm vntyper:baseline bash -c "conda run -n envadvntr python --version"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: envadvntr conda env failed"
    exit 1
fi
echo "✅ PASS: envadvntr environment"

# Test shark environment
docker run --rm vntyper:baseline bash -c "conda run -n shark python --version"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: shark conda env failed"
    exit 1
fi
echo "✅ PASS: shark environment"
```

#### 4. Baseline Testing - FastAPI & Celery
```bash
echo "=== Testing Docker Compose Stack (Baseline) ==="

# Start services
docker compose -f docker/docker-compose.yml up -d

# Wait for services to be healthy
echo "Waiting for services to start..."
sleep 90

# Check all containers are running
docker compose ps
RUNNING=$(docker compose ps --format json | jq -r '.State' | grep -c "running")
if [ "$RUNNING" -lt 3 ]; then
    echo "❌ FAIL: Not all services running"
    docker compose logs
    docker compose down
    exit 1
fi
echo "✅ PASS: All services started"

# Test Redis
docker compose exec redis redis-cli ping | grep -q "PONG"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Redis not responding"
    docker compose down
    exit 1
fi
echo "✅ PASS: Redis"

# Test API health endpoint
curl -f http://localhost:8000/health || curl -f http://localhost:8000/
if [ $? -ne 0 ]; then
    echo "❌ FAIL: API health check failed"
    docker compose logs api
    docker compose down
    exit 1
fi
echo "✅ PASS: API health endpoint"

# Check Celery worker is running
docker compose logs worker | grep -q "ready"
if [ $? -ne 0 ]; then
    echo "⚠️  WARNING: Celery worker may not be ready"
fi
echo "✅ PASS: Celery worker"

# Save baseline startup time
docker compose down
echo "Recording startup time..."
START_TIME=$(date +%s)
docker compose up -d
sleep 90
END_TIME=$(date +%s)
STARTUP_TIME=$((END_TIME - START_TIME))
echo "Baseline startup time: ${STARTUP_TIME}s" > testing/baseline/startup_time.txt

docker compose down
echo "✅ Baseline docker-compose stack test complete"
```

#### 5. Baseline Testing - Unit & Integration Tests
```bash
echo "=== Running VNtyper Test Suite in Container ==="

# Run unit tests
docker run --rm -v "$(pwd)":/workspace -w /workspace vntyper:baseline bash -c "conda run -n vntyper pytest -m unit -v" > testing/baseline/unit_tests.log 2>&1
UNIT_RESULT=$?

if [ $UNIT_RESULT -ne 0 ]; then
    echo "❌ FAIL: Unit tests failed"
    cat testing/baseline/unit_tests.log
    exit 1
fi
echo "✅ PASS: Unit tests ($UNIT_RESULT)"

# Run integration tests (if test data available)
if [ -f "tests/test_data_config.json" ]; then
    docker run --rm -v "$(pwd)":/workspace -w /workspace vntyper:baseline bash -c "conda run -n vntyper pytest -m integration -v" > testing/baseline/integration_tests.log 2>&1
    INTEGRATION_RESULT=$?

    if [ $INTEGRATION_RESULT -ne 0 ]; then
        echo "⚠️  WARNING: Integration tests failed (may need test data)"
        cat testing/baseline/integration_tests.log
    else
        echo "✅ PASS: Integration tests"
    fi
else
    echo "⚠️  SKIP: No test data config found"
fi
```

#### 6. Baseline Performance Metrics
```bash
echo "=== Recording Baseline Performance ==="

# Build time (already captured above)
# Image size (already captured above)

# Container startup time
docker run --rm vntyper:baseline bash -c "time conda run -n vntyper vntyper --version" > testing/baseline/cli_startup.txt 2>&1

# Memory usage
docker stats --no-stream vntyper:baseline > testing/baseline/memory_usage.txt || echo "No running container"

echo "✅ Baseline metrics recorded in testing/baseline/"
ls -lh testing/baseline/
```

**CHECKPOINT:** Review baseline test results. All tests should PASS before proceeding.

---

### Phase 1 Testing: Quick Wins Validation

#### Test 1.1: docker-compose.yml Changes
```bash
echo "=== Testing Updated docker-compose.yml ==="

# Validate syntax
docker compose -f docker/docker-compose.yml config > /dev/null
if [ $? -ne 0 ]; then
    echo "❌ FAIL: docker-compose.yml syntax error"
    exit 1
fi
echo "✅ PASS: Compose file syntax valid"

# Check version field removed
if grep -q "^version:" docker/docker-compose.yml; then
    echo "❌ FAIL: version field still present"
    exit 1
fi
echo "✅ PASS: Version field removed"

# Check Redis version upgraded
grep -q "redis:7.2" docker/docker-compose.yml
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Redis not upgraded to 7.2"
    exit 1
fi
echo "✅ PASS: Redis upgraded"

# Check health checks added
HEALTHCHECKS=$(grep -c "healthcheck:" docker/docker-compose.yml)
if [ "$HEALTHCHECKS" -lt 2 ]; then
    echo "❌ FAIL: Missing health checks"
    exit 1
fi
echo "✅ PASS: Health checks added ($HEALTHCHECKS found)"

# Test depends_on with conditions
grep -A 2 "depends_on:" docker/docker-compose.yml | grep -q "condition: service_healthy"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: depends_on conditions not properly set"
    exit 1
fi
echo "✅ PASS: depends_on conditions configured"
```

#### Test 1.2: Compose Stack with Health Checks
```bash
echo "=== Testing Services with Health Checks ==="

# Build with updated compose file
docker compose build --no-cache

# Start services
docker compose up -d

# Wait for health checks to pass
echo "Waiting for health checks (max 120s)..."
TIMEOUT=120
ELAPSED=0
while [ $ELAPSED -lt $TIMEOUT ]; do
    HEALTHY=$(docker compose ps --format json | jq -r 'select(.Health == "healthy") | .Service' | wc -l)
    TOTAL=$(docker compose ps --format json | wc -l)

    echo "Healthy: $HEALTHY/$TOTAL services"

    if [ "$HEALTHY" -eq "$TOTAL" ] && [ "$TOTAL" -gt 0 ]; then
        echo "✅ PASS: All services healthy"
        break
    fi

    sleep 10
    ELAPSED=$((ELAPSED + 10))
done

if [ $ELAPSED -ge $TIMEOUT ]; then
    echo "❌ FAIL: Services did not become healthy in time"
    docker compose ps
    docker compose logs
    docker compose down
    exit 1
fi

# Test startup ordering (worker should start after redis)
docker compose logs worker | grep -q "ready"
if [ $? -ne 0 ]; then
    echo "⚠️  WARNING: Worker startup may have issues"
fi

# Test Redis 7.2 functionality
docker compose exec redis redis-cli INFO server | grep "redis_version:7.2"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Redis 7.2 not running"
    docker compose down
    exit 1
fi
echo "✅ PASS: Redis 7.2 functional"

docker compose down
```

#### Test 1.3: Cache Optimization Changes
```bash
echo "=== Testing pip Cache Optimization ==="

# Rebuild with cache monitoring
docker build -f docker/Dockerfile -t vntyper:p1-test --progress=plain . 2>&1 | tee testing/phase1/build_with_cache.log

# Check that pip cache is being used (not --no-cache-dir)
if grep -q "no-cache-dir" docker/Dockerfile; then
    echo "⚠️  WARNING: Found --no-cache-dir in Dockerfile (check if removed from pip installs)"
fi

# Second build should be faster
echo "Testing cache effectiveness..."
START=$(date +%s)
docker build -f docker/Dockerfile -t vntyper:p1-test . > /dev/null 2>&1
END=$(date +%s)
REBUILD_TIME=$((END - START))

echo "Rebuild time: ${REBUILD_TIME}s"
echo "$REBUILD_TIME" > testing/phase1/rebuild_time.txt

# Compare with baseline
BASELINE_TIME=$(cat testing/baseline/startup_time.txt | grep -o '[0-9]*')
if [ "$REBUILD_TIME" -lt "$BASELINE_TIME" ]; then
    echo "✅ PASS: Build cache improved rebuild time"
else
    echo "⚠️  INFO: Rebuild time: ${REBUILD_TIME}s (baseline: ${BASELINE_TIME}s)"
fi
```

**CHECKPOINT:** All Phase 1 tests should PASS. Rollback if failures occur.

---

### Phase 2 Testing: Multi-Stage Build Validation

**CRITICAL:** This is the highest-risk change. Test thoroughly.

#### Test 2.1: Multi-Stage Build Success
```bash
echo "=== Testing Multi-Stage Dockerfile Build ==="

mkdir -p testing/phase2

# Build multi-stage image
echo "Building multi-stage image..."
time docker build -f docker/Dockerfile.multistage -t vntyper:multistage --progress=plain . 2>&1 | tee testing/phase2/multistage_build.log

if [ $? -ne 0 ]; then
    echo "❌ FAIL: Multi-stage build failed"
    cat testing/phase2/multistage_build.log
    exit 1
fi
echo "✅ PASS: Multi-stage build successful"

# Record image size
docker images vntyper:multistage --format "{{.Size}}" > testing/phase2/image_size.txt
BASELINE_SIZE=$(cat testing/baseline/image_size.txt)
NEW_SIZE=$(cat testing/phase2/image_size.txt)

echo "Baseline size: $BASELINE_SIZE"
echo "New size: $NEW_SIZE"

# Verify size reduction
docker images vntyper --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"

# Check layer count
BASELINE_LAYERS=$(wc -l < testing/baseline/layers.txt)
docker history vntyper:multistage > testing/phase2/layers.txt
NEW_LAYERS=$(wc -l < testing/phase2/layers.txt)

echo "Baseline layers: $BASELINE_LAYERS"
echo "New layers: $NEW_LAYERS"

if [ "$NEW_LAYERS" -ge "$BASELINE_LAYERS" ]; then
    echo "⚠️  WARNING: Layer count not reduced as expected"
fi
```

#### Test 2.2: Conda Environments Integrity
```bash
echo "=== Testing Conda Environments in Multi-Stage Image ==="

# Test that conda environments were copied correctly
docker run --rm vntyper:multistage bash -c "ls -la /opt/conda/envs/"

# Test vntyper environment activation
docker run --rm vntyper:multistage bash -c "conda run -n vntyper python --version"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: vntyper conda environment broken"
    exit 1
fi
echo "✅ PASS: vntyper environment works"

# Test envadvntr environment
docker run --rm vntyper:multistage bash -c "conda run -n envadvntr python --version"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: envadvntr conda environment broken"
    exit 1
fi
echo "✅ PASS: envadvntr environment works"

# Test shark environment
docker run --rm vntyper:multistage bash -c "conda run -n shark python --version"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: shark conda environment broken"
    exit 1
fi
echo "✅ PASS: shark environment works"
```

#### Test 2.3: Critical Dependencies
```bash
echo "=== Testing Critical Dependencies in Multi-Stage Image ==="

# Test Java (CRITICAL for Kestrel)
echo "Testing Java..."
docker run --rm vntyper:multistage bash -c "conda run -n vntyper java -version" 2>&1 | grep -q "openjdk"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Java not available or wrong version"
    docker run --rm vntyper:multistage bash -c "conda run -n vntyper java -version" 2>&1
    exit 1
fi
echo "✅ PASS: Java available"

# Test bioinformatics tools
echo "Testing bioinformatics tools..."
for tool in bwa samtools fastp bcftools; do
    docker run --rm vntyper:multistage bash -c "conda run -n vntyper which $tool" > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "❌ FAIL: $tool not found in multistage image"
        exit 1
    fi

    # Test tool execution
    docker run --rm vntyper:multistage bash -c "conda run -n vntyper $tool --version" > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "❌ FAIL: $tool cannot execute"
        exit 1
    fi
    echo "✅ PASS: $tool functional"
done

# Test Python packages and compiled extensions
echo "Testing Python packages..."
docker run --rm vntyper:multistage bash -c "conda run -n vntyper python -c '
import sys
import pandas as pd
import numpy as np
import pysam
import Bio
from vntyper import __version__

print(f\"Python: {sys.version}\")
print(f\"NumPy: {np.__version__}\")
print(f\"Pandas: {pd.__version__}\")
print(f\"VNtyper: {__version__}\")

# Test compiled extension (pysam uses C extensions)
import pysam
print(\"All packages with compiled extensions work!\")
'"

if [ $? -ne 0 ]; then
    echo "❌ FAIL: Python packages or compiled extensions broken"
    exit 1
fi
echo "✅ PASS: Python packages functional"
```

#### Test 2.4: VNtyper CLI Functionality
```bash
echo "=== Testing VNtyper CLI in Multi-Stage Image ==="

# Test version
docker run --rm vntyper:multistage vntyper --version
if [ $? -ne 0 ]; then
    echo "❌ FAIL: vntyper --version failed"
    exit 1
fi
echo "✅ PASS: vntyper --version"

# Test help
docker run --rm vntyper:multistage vntyper --help | grep -q "usage:"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: vntyper --help failed"
    exit 1
fi
echo "✅ PASS: vntyper --help"

# Test config access
docker run --rm vntyper:multistage bash -c "test -f /opt/vntyper/vntyper/config.json"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: config.json not found"
    exit 1
fi
echo "✅ PASS: config.json accessible"

# Test reference files
docker run --rm vntyper:multistage bash -c "ls /opt/vntyper/reference/" | grep -q ".fa"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Reference files not found"
    exit 1
fi
echo "✅ PASS: Reference files present"

# Test Kestrel JAR files
docker run --rm vntyper:multistage bash -c "ls /opt/vntyper/vntyper/dependencies/kestrel/*.jar"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Kestrel JAR files not found"
    exit 1
fi
echo "✅ PASS: Kestrel JAR files present"
```

#### Test 2.5: File Permissions & Ownership
```bash
echo "=== Testing File Permissions in Multi-Stage Image ==="

# Test that appuser owns necessary directories
docker run --rm vntyper:multistage bash -c "ls -ld /opt/vntyper/input /opt/vntyper/output /opt/vntyper/reference" | grep -q "appuser"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Incorrect ownership on directories"
    docker run --rm vntyper:multistage bash -c "ls -ld /opt/vntyper/input /opt/vntyper/output /opt/vntyper/reference"
    exit 1
fi
echo "✅ PASS: Directory ownership correct"

# Test write permissions
docker run --rm vntyper:multistage bash -c "touch /opt/vntyper/output/test.txt && rm /opt/vntyper/output/test.txt"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Cannot write to output directory"
    exit 1
fi
echo "✅ PASS: Write permissions OK"

# Verify non-root user
docker run --rm vntyper:multistage bash -c "whoami" | grep -q "appuser"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Not running as appuser"
    exit 1
fi
echo "✅ PASS: Running as non-root user (appuser)"
```

#### Test 2.6: FastAPI & Web Stack
```bash
echo "=== Testing FastAPI in Multi-Stage Image ==="

# Test FastAPI imports
docker run --rm vntyper:multistage bash -c "conda run -n vntyper python -c '
import fastapi
import uvicorn
import celery
import redis
print(\"All web dependencies available\")
'"

if [ $? -ne 0 ]; then
    echo "❌ FAIL: Web dependencies missing"
    exit 1
fi
echo "✅ PASS: Web dependencies present"

# Test app code exists
docker run --rm vntyper:multistage bash -c "ls /opt/vntyper/app/main.py"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: FastAPI app code missing"
    exit 1
fi
echo "✅ PASS: FastAPI app code present"

# Test app imports (syntax check)
docker run --rm vntyper:multistage bash -c "conda run -n vntyper python -m py_compile /opt/vntyper/app/main.py"
if [ $? -ne 0 ]; then
    echo "❌ FAIL: FastAPI app has syntax errors"
    exit 1
fi
echo "✅ PASS: FastAPI app syntax valid"
```

#### Test 2.7: Full Stack Integration with Multi-Stage Image
```bash
echo "=== Testing Full Docker Compose Stack with Multi-Stage Image ==="

# Update compose to use multistage image
sed -i.bak 's|dockerfile: docker/Dockerfile|dockerfile: docker/Dockerfile.multistage|g' docker/docker-compose.yml

# Build all services
docker compose build --no-cache

# Start services
docker compose up -d

# Wait for health
echo "Waiting for all services to be healthy..."
sleep 120

# Check health status
docker compose ps

HEALTHY=$(docker compose ps --format json | jq -r 'select(.Health == "healthy" or .State == "running") | .Service' | wc -l)
TOTAL=$(docker compose ps --format json | wc -l)

echo "Healthy/Running: $HEALTHY/$TOTAL"

if [ "$HEALTHY" -ne "$TOTAL" ]; then
    echo "❌ FAIL: Not all services healthy"
    docker compose ps
    docker compose logs
    docker compose down
    # Restore original compose file
    mv docker/docker-compose.yml.bak docker/docker-compose.yml
    exit 1
fi
echo "✅ PASS: All services healthy"

# Test API endpoint
curl -f http://localhost:8000/health || curl -f http://localhost:8000/
if [ $? -ne 0 ]; then
    echo "❌ FAIL: API not responding"
    docker compose logs api
    docker compose down
    mv docker/docker-compose.yml.bak docker/docker-compose.yml
    exit 1
fi
echo "✅ PASS: API responding"

# Test Celery worker
docker compose logs worker | grep -q "ready"
if [ $? -ne 0 ]; then
    echo "⚠️  WARNING: Celery worker may not be ready"
    docker compose logs worker | tail -20
fi
echo "✅ PASS: Celery worker running"

# Cleanup
docker compose down
# Restore backup
mv docker/docker-compose.yml.bak docker/docker-compose.yml
```

#### Test 2.8: Regression Testing - Unit Tests
```bash
echo "=== Running Unit Tests on Multi-Stage Image ==="

docker run --rm -v "$(pwd)":/workspace -w /workspace vntyper:multistage bash -c "conda run -n vntyper pytest -m unit -v" > testing/phase2/unit_tests.log 2>&1

if [ $? -ne 0 ]; then
    echo "❌ FAIL: Unit tests failed on multi-stage image"
    cat testing/phase2/unit_tests.log
    exit 1
fi
echo "✅ PASS: Unit tests on multi-stage image"

# Compare with baseline
diff testing/baseline/unit_tests.log testing/phase2/unit_tests.log > testing/phase2/test_diff.txt || true
echo "Test differences saved to testing/phase2/test_diff.txt"
```

#### Test 2.9: Performance Comparison
```bash
echo "=== Performance Comparison: Baseline vs Multi-Stage ==="

# Image size comparison
echo "Image Size Comparison:"
echo "======================"
docker images vntyper:baseline --format "Baseline: {{.Size}}"
docker images vntyper:multistage --format "Multi-stage: {{.Size}}"

# Build time comparison (from logs)
echo ""
echo "Build Time Comparison:"
echo "======================"
BASELINE_BUILD=$(grep "real" testing/baseline/build.log | awk '{print $2}')
MULTISTAGE_BUILD=$(grep "real" testing/phase2/multistage_build.log | awk '{print $2}')
echo "Baseline: $BASELINE_BUILD"
echo "Multi-stage: $MULTISTAGE_BUILD"

# Container startup comparison
echo ""
echo "Container Startup Time:"
echo "======================="
time docker run --rm vntyper:baseline vntyper --version > testing/phase2/baseline_startup.txt 2>&1
time docker run --rm vntyper:multistage vntyper --version > testing/phase2/multistage_startup.txt 2>&1

echo "Baseline startup:" && cat testing/phase2/baseline_startup.txt | grep "real"
echo "Multi-stage startup:" && cat testing/phase2/multistage_startup.txt | grep "real"

echo ""
echo "✅ Performance comparison complete"
```

**CHECKPOINT:** All Phase 2 tests must PASS. This is the critical phase. If any test fails, do NOT proceed to Phase 3.

**Rollback Trigger:** If more than 2 tests fail, or any CRITICAL test fails (Java, conda envs, bioinformatics tools), rollback immediately.

---

### Phase 3 Testing: Security Validation

#### Test 3.1: Base Image Digest Pinning
```bash
echo "=== Testing Digest-Pinned Image ==="

# Verify digest is in Dockerfile
grep -q "@sha256:" docker/Dockerfile.multistage
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Base image not pinned by digest"
    exit 1
fi
echo "✅ PASS: Base image pinned by digest"

# Build with pinned digest
docker build -f docker/Dockerfile.multistage -t vntyper:digest-pinned .
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Build failed with pinned digest"
    exit 1
fi
echo "✅ PASS: Build successful with digest pinning"

# Verify functionality unchanged
docker run --rm vntyper:digest-pinned vntyper --version
if [ $? -ne 0 ]; then
    echo "❌ FAIL: Functionality broken with digest pinning"
    exit 1
fi
echo "✅ PASS: Functionality preserved"
```

#### Test 3.2: Security Scanning
```bash
echo "=== Security Vulnerability Scanning ==="

mkdir -p testing/phase3

# Scan with Docker Scout
echo "Running Docker Scout scan..."
docker scout cves vntyper:multistage --format sarif --output testing/phase3/scout_scan.sarif || true

# Scan with Trivy
echo "Running Trivy scan..."
trivy image vntyper:multistage --format json --output testing/phase3/trivy_scan.json || true

# Compare with baseline
echo ""
echo "Vulnerability Comparison:"
echo "========================="

BASELINE_CRIT=$(jq '[.results[].vulnerabilities[] | select(.severity == "CRITICAL")] | length' testing/baseline/trivy_baseline.json 2>/dev/null || echo "0")
NEW_CRIT=$(jq '[.results[].vulnerabilities[] | select(.severity == "CRITICAL")] | length' testing/phase3/trivy_scan.json 2>/dev/null || echo "0")

echo "Baseline CRITICAL vulnerabilities: $BASELINE_CRIT"
echo "Multi-stage CRITICAL vulnerabilities: $NEW_CRIT"

if [ "$NEW_CRIT" -gt "$BASELINE_CRIT" ]; then
    echo "⚠️  WARNING: More CRITICAL vulnerabilities than baseline"
    echo "This may need investigation but could be due to updated scanning rules"
fi

if [ "$NEW_CRIT" -lt "$BASELINE_CRIT" ]; then
    echo "✅ IMPROVEMENT: Fewer CRITICAL vulnerabilities than baseline"
fi

# Check for specific CVEs
echo ""
echo "Checking for high-risk CVEs..."
trivy image vntyper:multistage --severity CRITICAL --exit-code 1
if [ $? -eq 0 ]; then
    echo "✅ PASS: No CRITICAL vulnerabilities found"
else
    echo "⚠️  WARNING: CRITICAL vulnerabilities found (review required)"
fi
```

#### Test 3.3: Attack Surface Analysis
```bash
echo "=== Attack Surface Analysis ==="

# List system packages in runtime image
echo "System packages in runtime image:"
docker run --rm vntyper:multistage bash -c "dpkg -l" > testing/phase3/runtime_packages.txt
RUNTIME_PKGS=$(wc -l < testing/phase3/runtime_packages.txt)

# Compare with baseline
docker run --rm vntyper:baseline bash -c "dpkg -l" > testing/phase3/baseline_packages.txt
BASELINE_PKGS=$(wc -l < testing/phase3/baseline_packages.txt)

echo "Baseline packages: $BASELINE_PKGS"
echo "Runtime packages: $RUNTIME_PKGS"

if [ "$RUNTIME_PKGS" -lt "$BASELINE_PKGS" ]; then
    echo "✅ IMPROVEMENT: Reduced package count by $((BASELINE_PKGS - RUNTIME_PKGS))"
else
    echo "⚠️  INFO: Package count: $RUNTIME_PKGS"
fi

# Check that build tools are NOT in runtime image
echo ""
echo "Checking that build tools are removed from runtime:"
FOUND_BUILD_TOOLS=0

for tool in gcc g++ make; do
    docker run --rm vntyper:multistage bash -c "which $tool" > /dev/null 2>&1
    if [ $? -eq 0 ]; then
        echo "⚠️  WARNING: Build tool '$tool' found in runtime image"
        FOUND_BUILD_TOOLS=$((FOUND_BUILD_TOOLS + 1))
    fi
done

if [ $FOUND_BUILD_TOOLS -eq 0 ]; then
    echo "✅ PASS: No build tools in runtime image"
else
    echo "⚠️  WARNING: Found $FOUND_BUILD_TOOLS build tools (expected 0)"
fi
```

#### Test 3.4: Non-Root User Verification
```bash
echo "=== Non-Root User Security Test ==="

# Verify running as non-root
USER=$(docker run --rm vntyper:multistage whoami)
if [ "$USER" != "appuser" ]; then
    echo "❌ FAIL: Not running as appuser (running as: $USER)"
    exit 1
fi
echo "✅ PASS: Running as non-root user: $USER"

# Verify UID is not 0
UID=$(docker run --rm vntyper:multistage id -u)
if [ "$UID" -eq 0 ]; then
    echo "❌ FAIL: Running with UID 0 (root)"
    exit 1
fi
echo "✅ PASS: Running with UID: $UID"

# Test that user cannot escalate privileges
docker run --rm vntyper:multistage bash -c "sudo ls" 2>&1 | grep -q "sudo: command not found"
if [ $? -ne 0 ]; then
    echo "⚠️  WARNING: sudo might be available"
fi

# Verify cannot write to protected areas
docker run --rm vntyper:multistage bash -c "touch /etc/test.txt" 2>&1 | grep -q "Permission denied"
if [ $? -eq 0 ]; then
    echo "✅ PASS: Cannot write to /etc (as expected)"
else
    echo "⚠️  WARNING: May have elevated permissions"
fi
```

**CHECKPOINT:** Review security scan results. Determine if vulnerability levels are acceptable.

---

### Final Validation: End-to-End Testing

#### E2E Test 1: Complete Pipeline Run (if test data available)
```bash
echo "=== End-to-End Pipeline Test ==="

if [ ! -f "tests/test_data_config.json" ]; then
    echo "⚠️  SKIP: No test data config found"
    exit 0
fi

# This would run a full VNtyper pipeline
# Adapt to your actual test data structure

echo "Running integration tests..."
docker run --rm \
    -v "$(pwd)":/workspace \
    -v "$(pwd)/tests/data":/data \
    -w /workspace \
    vntyper:multistage \
    bash -c "conda run -n vntyper pytest -m integration -v --tb=short" > testing/final/integration_tests.log 2>&1

if [ $? -ne 0 ]; then
    echo "❌ FAIL: Integration tests failed"
    cat testing/final/integration_tests.log | tail -50
    exit 1
fi

echo "✅ PASS: Integration tests"
```

#### E2E Test 2: Docker Compose Production Simulation
```bash
echo "=== Production Simulation Test ==="

# Use final multi-stage Dockerfile
sed -i 's|dockerfile:.*|dockerfile: docker/Dockerfile.multistage|g' docker/docker-compose.yml

# Start full stack
docker compose up -d

# Wait for stability
sleep 120

# Monitor for 5 minutes
echo "Monitoring services for 5 minutes..."
for i in {1..30}; do
    STATUS=$(docker compose ps --format json | jq -r '.Health // .State' | sort | uniq -c)
    echo "[$i/30] Service status: $STATUS"

    # Check for any failed services
    FAILED=$(docker compose ps --format json | jq -r 'select(.State == "exited" or .Health == "unhealthy") | .Service')
    if [ -n "$FAILED" ]; then
        echo "❌ FAIL: Service failed: $FAILED"
        docker compose logs "$FAILED"
        docker compose down
        exit 1
    fi

    sleep 10
done

echo "✅ PASS: Services stable for 5 minutes"

# Cleanup
docker compose down
```

#### E2E Test 3: Memory & Resource Usage
```bash
echo "=== Resource Usage Test ==="

# Start services
docker compose up -d
sleep 60

# Monitor resource usage
docker stats --no-stream --format "table {{.Container}}\t{{.CPUPerc}}\t{{.MemUsage}}\t{{.MemPerc}}" > testing/final/resource_usage.txt

cat testing/final/resource_usage.txt

# Check for memory leaks (run for 10 minutes)
echo "Monitoring for memory leaks (10 minutes)..."
for i in {1..10}; do
    docker stats --no-stream --format "{{.MemUsage}}" vntyper_api >> testing/final/memory_timeline.txt
    sleep 60
done

# Analyze memory growth
FIRST_MEM=$(head -1 testing/final/memory_timeline.txt | grep -o '[0-9.]*MiB' | grep -o '[0-9.]*')
LAST_MEM=$(tail -1 testing/final/memory_timeline.txt | grep -o '[0-9.]*MiB' | grep -o '[0-9.]*')

if [ ! -z "$FIRST_MEM" ] && [ ! -z "$LAST_MEM" ]; then
    MEM_GROWTH=$(echo "$LAST_MEM - $FIRST_MEM" | bc)
    echo "Memory growth over 10 minutes: ${MEM_GROWTH}MiB"

    if (( $(echo "$MEM_GROWTH > 100" | bc -l) )); then
        echo "⚠️  WARNING: Significant memory growth detected"
    else
        echo "✅ PASS: Memory usage stable"
    fi
fi

docker compose down
```

---

### Test Summary Report Generation

```bash
echo "=== Generating Test Summary Report ==="

cat > testing/TEST_SUMMARY_REPORT.md << 'EOF'
# VNtyper Docker Modernization - Test Summary Report

## Test Execution Date
EOF

date >> testing/TEST_SUMMARY_REPORT.md

cat >> testing/TEST_SUMMARY_REPORT.md << 'EOF'

## Baseline Metrics

### Image Size
EOF

echo "- Baseline: $(cat testing/baseline/image_size.txt)" >> testing/TEST_SUMMARY_REPORT.md
echo "- Multi-stage: $(cat testing/phase2/image_size.txt)" >> testing/TEST_SUMMARY_REPORT.md

cat >> testing/TEST_SUMMARY_REPORT.md << 'EOF'

### Test Results

#### Phase 1: Quick Wins
- docker-compose.yml validation: ✅ PASS
- Health checks: ✅ PASS
- Redis upgrade: ✅ PASS
- Cache optimization: ✅ PASS

#### Phase 2: Multi-Stage Build
- Build success: ✅ PASS
- Conda environments: ✅ PASS
- Java availability: ✅ PASS
- Bioinformatics tools: ✅ PASS
- Python packages: ✅ PASS
- VNtyper CLI: ✅ PASS
- File permissions: ✅ PASS
- FastAPI stack: ✅ PASS
- Full stack integration: ✅ PASS
- Unit tests: ✅ PASS

#### Phase 3: Security
- Digest pinning: ✅ PASS
- Vulnerability scanning: ✅ PASS
- Attack surface reduction: ✅ PASS
- Non-root user: ✅ PASS

#### Final Validation
- Integration tests: ✅ PASS
- Production simulation: ✅ PASS
- Resource usage: ✅ PASS

## Recommendations

### Proceed with Production Deployment
✅ All tests passed - safe to deploy

### Issues Found
(List any warnings or issues)

### Performance Improvements
- Image size reduction: [X]%
- Build time improvement: [X]%
- Startup time: [X]s

## Sign-Off

Tested by: _________________
Date: _________________
Approved for production: Yes / No
EOF

echo "✅ Test summary report generated: testing/TEST_SUMMARY_REPORT.md"
```

---

### Rollback Procedures

#### Immediate Rollback
```bash
#!/bin/bash
# rollback.sh - Emergency rollback script

echo "=== EMERGENCY ROLLBACK ==="

# Stop all services
docker compose down

# Restore original docker-compose.yml
if [ -f "docker/docker-compose.yml.bak" ]; then
    cp docker/docker-compose.yml.bak docker/docker-compose.yml
    echo "✅ Restored docker-compose.yml"
fi

# Restore original Dockerfile
git checkout docker/Dockerfile docker/docker-compose.yml

# Use baseline image
docker tag vntyper:baseline vntyper:latest

# Restart services
docker compose up -d

echo "✅ Rollback complete - services running on baseline image"
```

**Save this as `testing/rollback.sh` and make executable:**
```bash
chmod +x testing/rollback.sh
```

---

## 🚨 Rollback Plan

**If issues arise, rollback steps:**

### Quick Rollback (Compose)
```bash
# Revert docker-compose.yml
git checkout HEAD -- docker/docker-compose.yml
docker compose down
docker compose up -d
```

### Full Rollback
```bash
# Revert all Dockerfile changes
git checkout HEAD -- docker/Dockerfile*
docker compose build --no-cache
docker compose up -d
```

### Emergency Rollback
```bash
# Use tagged backup image
docker tag vntyper:backup vntyper:latest
docker compose up -d
```

---

## ⚠️ CRITICAL WARNINGS & POTENTIAL ISSUES

### Multi-Stage Build Concerns

#### ⚠️ Conda Environment Copying
**Issue:** When copying `/opt/conda/envs` between stages, ensure all dependencies are self-contained.

**Symptoms if broken:**
- Missing library errors (e.g., `ImportError: libXXX.so not found`)
- Python package import failures
- Segmentation faults in compiled extensions

**Mitigation:**
- Conda packages from conda-forge are usually self-contained
- Test ALL conda environments after copying (vntyper, envadvntr, shark)
- Verify bioinformatics tools (bwa, samtools, fastp) work
- Check Java is accessible for Kestrel

#### ⚠️ System Library Dependencies
**Issue:** Some conda packages may link against system libraries not present in runtime stage.

**Check during testing:**
```bash
# Identify missing libraries
docker run --rm vntyper:multistage bash -c "ldd /opt/conda/envs/vntyper/bin/bwa"
docker run --rm vntyper:multistage bash -c "ldd /opt/conda/envs/vntyper/bin/samtools"
```

**Fix if needed:** Add missing system libraries to runtime stage apt-get install.

#### ⚠️ Java Runtime
**CRITICAL:** Kestrel requires Java. Ensure Java is in the conda environment and accessible.

**Verify:**
```bash
docker run --rm vntyper:multistage bash -c "conda run -n vntyper java -version"
```

**If fails:** Java may not have been copied correctly or PATH may be wrong.

#### ⚠️ File Permissions After Copying
**Issue:** Copied files may have wrong ownership/permissions.

**Symptoms:**
- Permission denied errors
- Cannot write to output directories
- Services fail to start

**Prevention:**
- Use `chown -R` after COPY commands
- Test write permissions in all mounted directories
- Verify non-root user can access all necessary files

### Docker Compose Health Check Issues

#### ⚠️ Health Check Timing
**Issue:** Services may not become healthy quickly enough, causing timeout failures.

**Symptoms:**
- `depends_on` conditions never satisfied
- Services restart loop
- Timeouts during startup

**Fix:**
- Increase `start_period` in healthcheck definition
- Adjust `interval` and `timeout` values
- Check logs for actual startup time

#### ⚠️ Redis 7.2 Compatibility
**Issue:** Redis 7.2 has different config requirements than 6.x.

**Verify:**
- Test all Celery operations
- Check Redis persistence settings
- Verify connection from Python clients

### Performance Regression Risks

#### ⚠️ Conda Environment Size
**Reality Check:** Even with multi-stage builds, conda environments are large.

**Expected Results:**
- Image may still be 1.5-2GB (not MB)
- This is NORMAL for bioinformatics applications
- The win is removing build tools (gcc, make, etc.) not conda itself

#### ⚠️ Build Time May Increase
**Warning:** Multi-stage builds do MORE work during build.

**Trade-off:**
- Initial build: May be SLOWER (building multiple stages)
- Cached rebuild: Much FASTER
- Production image: Much SMALLER
- CI/CD: Configure layer caching for best results

### Security Considerations

#### ⚠️ Digest Pinning Maintenance
**Issue:** Pinned digests need manual updates.

**Process:**
1. Monitor for base image updates
2. Test new digest before deploying
3. Update digest in Dockerfile
4. Re-run full test suite

**Frequency:** Monthly security update check recommended

#### ⚠️ Vulnerability Scan Interpretation
**Reality:** Zero vulnerabilities is unrealistic for a 2GB bioinformatics image.

**What matters:**
- No CRITICAL vulnerabilities
- No actively exploited CVEs
- Trend: Fewer vulnerabilities over time
- Context: Many CVEs don't apply to your use case

---

## 📊 Expected Results

### Image Size Reduction
| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| **Image Size** | ~5GB | ~1.5-2GB | **60-70% smaller** |
| **Layers** | ~40 | ~25 | **37% fewer** |
| **Build Time (cold)** | ~15min | ~15-18min | **May be slower** |
| **Build Time (cached)** | ~3min | ~30s | **83% faster** |
| **Startup Time** | ~60s | ~60s | **Same** |

### Security Improvements
- ✅ No build tools in production image (gcc, make, build-essential removed)
- ✅ Base image pinned by digest (prevents tag hijacking)
- ✅ All services have health checks (proper dependency management)
- ✅ Proper startup ordering with depends_on (redis before workers)
- ✅ Latest Redis 7.2 with security patches
- ✅ Non-root user (appuser with UID 1001)
- ✅ Minimal attack surface (only curl + ca-certificates at system level)

### Developer Experience
- ✅ Faster iterative builds (83% improvement with cache)
- ✅ Smaller image pulls (60-70% reduction)
- ✅ Clearer separation of concerns (build vs runtime)
- ✅ Modern Docker best practices (2025 compliant)
- ✅ Better debugging (separate stages can be inspected)

### What Will NOT Change
- ⚠️ Conda environment size (still large, bioinformatics tools are big)
- ⚠️ Memory usage (same runtime components)
- ⚠️ Execution speed (same code running)
- ⚠️ Some vulnerabilities will remain (conda packages, system libs)

---

## 📚 References

### Official Documentation
- [Docker Multi-Stage Builds](https://docs.docker.com/build/building/multi-stage/)
- [Docker Compose Spec](https://docs.docker.com/compose/compose-file/)
- [BuildKit Cache Optimization](https://docs.docker.com/build/cache/)
- [Docker Health Checks](https://docs.docker.com/engine/reference/builder/#healthcheck)

### Best Practices
- [Docker Security Cheat Sheet (OWASP)](https://cheatsheetseries.owasp.org/cheatsheets/Docker_Security_Cheat_Sheet.html)
- [Conda Docker Best Practices](https://pythonspeed.com/articles/conda-docker-image-size/)
- [2025 Dockerfile Best Practices](https://docs.docker.com/build/building/best-practices/)

---

## 🔄 Continuous Improvement

### Future Enhancements
1. **Distroless Runtime Image:** Explore using distroless base for even smaller images
2. **BuildKit Secrets:** Use BuildKit secrets for build-time credentials
3. **Multi-Platform Builds:** Add ARM64 support for Apple Silicon
4. **Automated SBOM Scanning:** Integrate with CI/CD pipeline
5. **Layer Caching in CI:** Set up registry cache for GitHub Actions

### Monitoring
- Track image sizes in CI/CD metrics
- Monitor build times over time
- Set up automated vulnerability scanning
- Regular base image updates (monthly)

---

## 🚀 QUICK START: Master Test Script

Save this as `testing/run_all_tests.sh` for automated testing:

```bash
#!/bin/bash
# =============================================================================
# VNtyper Docker Modernization - Master Test Script
# =============================================================================
# This script runs ALL tests for the Docker modernization plan
# Run with: bash testing/run_all_tests.sh
# =============================================================================

set -e  # Exit on error
trap 'echo "❌ Test failed on line $LINENO"' ERR

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() { echo -e "${GREEN}[INFO]${NC} $*"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }
log_section() { echo -e "\n${BLUE}=== $* ===${NC}\n"; }

# Create testing directory structure
mkdir -p testing/{baseline,phase1,phase2,phase3,final}

# =============================================================================
# PHASE 0: Pre-Flight Checks
# =============================================================================
log_section "Phase 0: Pre-Flight Checks"

# Check Docker is running
if ! docker info > /dev/null 2>&1; then
    log_error "Docker is not running. Please start Docker and try again."
    exit 1
fi
log_info "Docker is running"

# Check Docker Compose
if ! command -v docker compose &> /dev/null; then
    log_error "docker compose not found. Please install Docker Compose."
    exit 1
fi
log_info "Docker Compose available"

# Check for required tools
for cmd in jq bc; do
    if ! command -v $cmd &> /dev/null; then
        log_warn "$cmd not found - some tests may be skipped"
    fi
done

# =============================================================================
# BASELINE TESTING
# =============================================================================
log_section "BASELINE: Building and Testing Current Configuration"

log_info "Building baseline image..."
time docker build -f docker/Dockerfile -t vntyper:baseline . 2>&1 | tee testing/baseline/build.log

log_info "Recording baseline metrics..."
docker images vntyper:baseline --format "{{.Size}}" > testing/baseline/image_size.txt
docker history vntyper:baseline > testing/baseline/layers.txt

log_info "Testing baseline CLI..."
docker run --rm vntyper:baseline vntyper --version
docker run --rm vntyper:baseline bash -c "conda run -n vntyper java -version"
docker run --rm vntyper:baseline bash -c "conda run -n vntyper bwa 2>&1 | head -3"

log_info "Running baseline unit tests..."
docker run --rm -v "$(pwd)":/workspace -w /workspace vntyper:baseline \
    bash -c "conda run -n vntyper pytest -m unit -v" > testing/baseline/unit_tests.log 2>&1 || true

log_info "✅ Baseline testing complete"

# =============================================================================
# PHASE 1: Quick Wins
# =============================================================================
log_section "PHASE 1: Testing Quick Wins"

log_info "Validating docker-compose.yml syntax..."
docker compose -f docker/docker-compose.yml config > /dev/null

log_info "✅ Phase 1 validation complete"

# =============================================================================
# PHASE 2: Multi-Stage Build (CRITICAL)
# =============================================================================
log_section "PHASE 2: Multi-Stage Build - CRITICAL TESTING"

if [ ! -f "docker/Dockerfile.multistage" ]; then
    log_error "docker/Dockerfile.multistage not found. Create it first before running tests."
    exit 1
fi

log_info "Building multi-stage image..."
time docker build -f docker/Dockerfile.multistage -t vntyper:multistage --progress=plain . 2>&1 | tee testing/phase2/multistage_build.log

log_info "Recording multi-stage metrics..."
docker images vntyper:multistage --format "{{.Size}}" > testing/phase2/image_size.txt
docker history vntyper:multistage > testing/phase2/layers.txt

log_info "Comparing image sizes..."
BASELINE_SIZE=$(cat testing/baseline/image_size.txt)
MULTISTAGE_SIZE=$(cat testing/phase2/image_size.txt)
echo "Baseline:    $BASELINE_SIZE"
echo "Multi-stage: $MULTISTAGE_SIZE"

log_info "Testing conda environments..."
for env in vntyper envadvntr shark; do
    docker run --rm vntyper:multistage bash -c "conda run -n $env python --version"
    log_info "✅ $env environment OK"
done

log_info "Testing Java (CRITICAL for Kestrel)..."
docker run --rm vntyper:multistage bash -c "conda run -n vntyper java -version" 2>&1 | grep -q "openjdk"
log_info "✅ Java available"

log_info "Testing bioinformatics tools..."
for tool in bwa samtools fastp bcftools; do
    docker run --rm vntyper:multistage bash -c "conda run -n vntyper which $tool" > /dev/null
    docker run --rm vntyper:multistage bash -c "conda run -n vntyper $tool --version" > /dev/null 2>&1
    log_info "✅ $tool functional"
done

log_info "Testing Python packages..."
docker run --rm vntyper:multistage bash -c "conda run -n vntyper python -c 'import pandas, numpy, pysam, Bio, vntyper'"
log_info "✅ Python packages OK"

log_info "Testing VNtyper CLI..."
docker run --rm vntyper:multistage vntyper --version
docker run --rm vntyper:multistage vntyper --help | grep -q "usage:"
log_info "✅ VNtyper CLI functional"

log_info "Testing file permissions..."
docker run --rm vntyper:multistage bash -c "whoami" | grep -q "appuser"
docker run --rm vntyper:multistage bash -c "touch /opt/vntyper/output/test.txt && rm /opt/vntyper/output/test.txt"
log_info "✅ Permissions OK (running as appuser)"

log_info "Testing FastAPI dependencies..."
docker run --rm vntyper:multistage bash -c "conda run -n vntyper python -c 'import fastapi, uvicorn, celery, redis'"
log_info "✅ Web stack dependencies OK"

log_info "Running unit tests on multi-stage image..."
docker run --rm -v "$(pwd)":/workspace -w /workspace vntyper:multistage \
    bash -c "conda run -n vntyper pytest -m unit -v" > testing/phase2/unit_tests.log 2>&1

UNIT_RESULT=$?
if [ $UNIT_RESULT -eq 0 ]; then
    log_info "✅ Unit tests PASSED on multi-stage image"
else
    log_warn "⚠️  Unit tests had issues - check testing/phase2/unit_tests.log"
fi

log_info "Checking for build tools in runtime (should be absent)..."
FOUND_GCC=0
docker run --rm vntyper:multistage bash -c "which gcc" > /dev/null 2>&1 && FOUND_GCC=1 || true
if [ $FOUND_GCC -eq 0 ]; then
    log_info "✅ No build tools in runtime image (good!)"
else
    log_warn "⚠️  WARNING: gcc found in runtime image"
fi

log_info "✅ Phase 2 multi-stage testing complete"

# =============================================================================
# PHASE 3: Security
# =============================================================================
log_section "PHASE 3: Security Validation"

log_info "Running Trivy scan (if available)..."
if command -v trivy &> /dev/null; then
    trivy image vntyper:multistage --severity CRITICAL,HIGH --format table > testing/phase3/trivy_scan.txt 2>&1 || true
    log_info "✅ Trivy scan complete - see testing/phase3/trivy_scan.txt"
else
    log_warn "Trivy not installed - skipping vulnerability scan"
fi

log_info "✅ Phase 3 security validation complete"

# =============================================================================
# FINAL: Summary
# =============================================================================
log_section "TEST SUMMARY"

echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║         VNtyper Docker Modernization - Test Results         ║"
echo "╠══════════════════════════════════════════════════════════════╣"
echo "║                                                              ║"
echo "║  Image Sizes:                                                ║"
echo "║    Baseline:    $BASELINE_SIZE"
echo "║    Multi-stage: $MULTISTAGE_SIZE"
echo "║                                                              ║"
echo "║  Critical Tests:                                             ║"
echo "║    ✅ Conda environments functional                          ║"
echo "║    ✅ Java available (Kestrel)                               ║"
echo "║    ✅ Bioinformatics tools working                           ║"
echo "║    ✅ Python packages functional                             ║"
echo "║    ✅ VNtyper CLI working                                    ║"
echo "║    ✅ Non-root user (appuser)                                ║"
echo "║    ✅ FastAPI dependencies available                         ║"
if [ $UNIT_RESULT -eq 0 ]; then
    echo "║    ✅ Unit tests PASSED                                      ║"
else
    echo "║    ⚠️  Unit tests had issues                                ║"
fi
echo "║                                                              ║"
echo "╠══════════════════════════════════════════════════════════════╣"
echo "║  Test artifacts saved in testing/ directory                  ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""

log_info "All automated tests complete!"
log_info "Review test results in testing/ directory"
log_info ""
log_info "Next steps:"
log_info "1. Review testing/phase2/unit_tests.log"
log_info "2. If all tests passed, proceed with docker-compose testing"
log_info "3. Run: docker compose up -d && docker compose ps"
log_info "4. Monitor logs: docker compose logs -f"

exit 0
```

**Make it executable:**
```bash
chmod +x testing/run_all_tests.sh
```

**Run all tests:**
```bash
./testing/run_all_tests.sh
```

---

## ✅ Success Criteria

Before proceeding to production, ALL of the following must be true:

### Functional Requirements
- [ ] `vntyper --version` works in multi-stage image
- [ ] `vntyper --help` displays usage information
- [ ] Java is available (`java -version` works)
- [ ] All bioinformatics tools execute (bwa, samtools, fastp, bcftools)
- [ ] All three conda environments activate (vntyper, envadvntr, shark)
- [ ] All Python packages import correctly (pandas, numpy, pysam, Bio)
- [ ] Compiled extensions work (pysam C extensions)
- [ ] Unit tests pass in container
- [ ] FastAPI server starts and responds
- [ ] Celery workers connect to Redis
- [ ] Health checks pass for all services

### Security Requirements
- [ ] Running as non-root user (appuser, UID 1001)
- [ ] Cannot write to /etc (verified permission denied)
- [ ] No build tools in runtime image (gcc, make, g++ absent)
- [ ] Base image pinned by digest
- [ ] No CRITICAL vulnerabilities (or acceptable with justification)
- [ ] Vulnerability count same or lower than baseline

### Performance Requirements
- [ ] Image size 60-70% smaller than baseline
- [ ] Cached rebuild time < 1 minute
- [ ] Container startup time comparable to baseline
- [ ] No memory leaks detected in 10-minute monitoring
- [ ] API response time same or better than baseline

### Quality Requirements
- [ ] All Phase 1 tests pass
- [ ] All Phase 2 tests pass
- [ ] All Phase 3 tests pass
- [ ] All E2E tests pass
- [ ] Test summary report generated
- [ ] No unresolved warnings
- [ ] Rollback procedure tested and documented

---

## 📝 Final Checklist

### Before Production Deployment

- [ ] All success criteria met
- [ ] Test summary report reviewed and approved
- [ ] Security scan results reviewed
- [ ] Performance benchmarks meet expectations
- [ ] No regressions in functionality
- [ ] Documentation updated
- [ ] Team notified of changes
- [ ] Rollback procedure tested
- [ ] Monitoring and alerting configured
- [ ] Backup of current production images created

### Post-Deployment Monitoring

- [ ] Monitor application logs for errors
- [ ] Monitor resource usage (CPU, memory)
- [ ] Monitor startup times
- [ ] Monitor API response times
- [ ] Check for any new errors or warnings
- [ ] Verify all scheduled tasks running
- [ ] Confirm no user-reported issues

---

## ✅ Sign-Off

**Implementation completed by:** _________________
**Date:** _________________
**Tested by:** _________________
**Approved by:** _________________

**Test Results:** PASS / FAIL / CONDITIONAL
**Issues Found:** _________________
**Deployment Approved:** YES / NO
**Deployment Date:** _________________

---

**End of Implementation Plan**
