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
cd docker && docker compose config > /dev/null && cd ..

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
echo "║    Baseline:    $BASELINE_SIZE                               "
echo "║    Multi-stage: $MULTISTAGE_SIZE                             "
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
log_info "3. Run: cd docker && docker compose up -d && docker compose ps"
log_info "4. Monitor logs: docker compose logs -f"

exit 0
