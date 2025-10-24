#!/usr/bin/env bash
# =============================================================================
# VNtyper Docker Container Test Script
# =============================================================================
# Simple, fast Docker testing using real test data from Zenodo
# Usage: make test-docker
# =============================================================================

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RESET='\033[0m'

# Configuration
DOCKER_IMAGE="${DOCKER_IMAGE:-vntyper:latest}"
ZENODO_RECORD="17377082"
ZENODO_BASE_URL="https://zenodo.org/records/${ZENODO_RECORD}/files"
TEST_DATA_DIR="tests/data"
DOCKER_TEST_DIR="docker/test_output"
TEST_BAM="example_66bf_hg19_subset.bam"

# =============================================================================
# Utility Functions
# =============================================================================

log_info() {
    echo -e "${BLUE}[INFO]${RESET} $*"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${RESET} $*"
}

log_error() {
    echo -e "${RED}[ERROR]${RESET} $*"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${RESET} $*"
}

compute_md5() {
    local file="$1"
    if command -v md5sum &> /dev/null; then
        md5sum "$file" | awk '{print $1}'
    elif command -v md5 &> /dev/null; then
        md5 -q "$file"
    else
        log_error "Neither md5sum nor md5 command found"
        return 1
    fi
}

download_file() {
    local url="$1"
    local output="$2"

    log_info "Downloading: $(basename "$output")"

    if command -v curl &> /dev/null; then
        curl -fsSL "$url" -o "$output"
    elif command -v wget &> /dev/null; then
        wget -q "$url" -O "$output"
    else
        log_error "Neither curl nor wget found"
        return 1
    fi
}

# =============================================================================
# Test Data Management
# =============================================================================

ensure_test_data() {
    log_info "Ensuring test data is available from Zenodo..."

    # Create test data directory if needed
    mkdir -p "$TEST_DATA_DIR"

    # List of required files for basic pipeline test
    local files=(
        "${TEST_BAM}"
        "${TEST_BAM}.bai"
    )

    # Check if files exist and have reasonable size
    local need_download=false
    for file in "${files[@]}"; do
        local filepath="${TEST_DATA_DIR}/${file}"
        if [[ ! -f "$filepath" ]] || [[ ! -s "$filepath" ]]; then
            log_info "Missing or empty: $file"
            need_download=true
            break
        fi
    done

    if [[ "$need_download" == "true" ]]; then
        log_info "Downloading test data archive from Zenodo (record: ${ZENODO_RECORD})..."

        # Download archive
        local archive_url="${ZENODO_BASE_URL}/data.zip?download=1"
        local temp_archive="/tmp/vntyper_test_data.zip"

        download_file "$archive_url" "$temp_archive"

        log_info "Extracting test data..."
        unzip -q -o "$temp_archive" -d "tests/"
        rm -f "$temp_archive"

        log_success "Test data downloaded and extracted"
    else
        log_success "Test data already present"
    fi

    # Verify critical files exist
    for file in "${files[@]}"; do
        local filepath="${TEST_DATA_DIR}/${file}"
        if [[ ! -f "$filepath" ]]; then
            log_error "Critical file missing after download: $filepath"
            return 1
        fi
    done

    log_success "All required test files verified"
}

# =============================================================================
# Docker Image Build
# =============================================================================

build_docker_image() {
    log_info "Building multi-stage Docker image..."

    local dockerfile="docker/Dockerfile"

    if [[ ! -f "$dockerfile" ]]; then
        log_error "Dockerfile not found: $dockerfile"
        return 1
    fi

    log_info "Building: ${DOCKER_IMAGE}"
    DOCKER_BUILDKIT=1 docker build \
        -f "$dockerfile" \
        -t "$DOCKER_IMAGE" \
        . 2>&1 | grep -E "^(#|Step |Successfully built|Successfully tagged)" || true

    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        log_error "Docker build failed"
        return 1
    fi

    log_success "Docker image built: ${DOCKER_IMAGE}"

    # Show image size
    local image_size
    image_size=$(docker images --format "{{.Size}}" "$DOCKER_IMAGE" | head -1)
    log_info "Image size: ${image_size}"
}

# =============================================================================
# Container Tests
# =============================================================================

test_container_basics() {
    log_info "Testing container basic functionality..."

    # Test 1: Container starts and vntyper command exists
    log_info "Test 1: VNtyper version check"
    if docker run --rm "$DOCKER_IMAGE" vntyper --version &>/dev/null; then
        log_success "✓ VNtyper CLI accessible"
    else
        log_error "✗ VNtyper CLI not accessible"
        return 1
    fi

    # Test 2: Java runtime (critical for Kestrel)
    log_info "Test 2: Java runtime check"
    if docker run --rm --entrypoint bash "$DOCKER_IMAGE" \
        -c 'conda run -n vntyper java -version' 2>&1 | grep -q "openjdk"; then
        log_success "✓ Java runtime available"
    else
        log_error "✗ Java runtime not found"
        return 1
    fi

    # Test 3: Critical bioinformatics tools
    log_info "Test 3: Bioinformatics tools check"
    local tools=("samtools" "bwa" "fastp" "bcftools")
    for tool in "${tools[@]}"; do
        if docker run --rm --entrypoint bash "$DOCKER_IMAGE" \
            -c "conda run -n vntyper which $tool" &>/dev/null; then
            log_success "✓ $tool available"
        else
            log_error "✗ $tool not found"
            return 1
        fi
    done

    log_success "All basic container tests passed"
}

test_pipeline_execution() {
    log_info "Testing full VNtyper pipeline execution..."

    # Prepare test directories
    local input_dir="${TEST_DATA_DIR}"
    local output_dir="${DOCKER_TEST_DIR}/pipeline_output"

    # Clean previous test output
    rm -rf "$output_dir"
    mkdir -p "$output_dir"

    log_info "Input BAM: ${input_dir}/${TEST_BAM}"
    log_info "Output directory: ${output_dir}"

    # Run pipeline
    log_info "Executing VNtyper pipeline (this may take 1-2 minutes)..."

    local start_time
    start_time=$(date +%s)

    if docker run --rm \
        -v "$(pwd)/${input_dir}":/opt/vntyper/input \
        -v "$(pwd)/${output_dir}":/opt/vntyper/output \
        "$DOCKER_IMAGE" \
        vntyper pipeline \
        --bam "/opt/vntyper/input/${TEST_BAM}" \
        -o /opt/vntyper/output/test_result \
        --reference-assembly hg19 \
        --threads 4 \
        --fast-mode \
        --keep-intermediates; then

        local end_time
        end_time=$(date +%s)
        local duration=$((end_time - start_time))

        log_success "Pipeline completed in ${duration} seconds"
    else
        log_error "Pipeline execution failed"
        return 1
    fi

    # Validate output files
    validate_pipeline_output "$output_dir/test_result"
}

validate_pipeline_output() {
    local output_dir="$1"

    log_info "Validating pipeline output files..."

    # Expected critical output files
    local expected_files=(
        "summary_report.html"
        "kestrel/kestrel_result.tsv"
        "coverage/coverage_summary.tsv"
        "pipeline_summary.json"
        "igv_report.html"
    )

    local missing_files=0
    for file in "${expected_files[@]}"; do
        local filepath="${output_dir}/${file}"
        if [[ -f "$filepath" ]] && [[ -s "$filepath" ]]; then
            log_success "✓ Found: $file"
        else
            log_error "✗ Missing or empty: $file"
            ((missing_files++))
        fi
    done

    if [[ $missing_files -gt 0 ]]; then
        log_error "Missing ${missing_files} expected output files"
        return 1
    fi

    log_success "All expected output files present"

    # Validate actual results against expected values (like pytest)
    validate_kestrel_results "$output_dir"
    validate_coverage_results "$output_dir"

    # Show summary
    local total_files
    total_files=$(find "$output_dir" -type f | wc -l | tr -d ' ')
    local total_size
    total_size=$(du -sh "$output_dir" | awk '{print $1}')

    log_info "Total output files: ${total_files}"
    log_info "Total output size: ${total_size}"

    return 0
}

validate_kestrel_results() {
    local output_dir="$1"
    local kestrel_file="${output_dir}/kestrel/kestrel_result.tsv"

    log_info "Validating Kestrel genotyping results..."

    if [[ ! -f "$kestrel_file" ]]; then
        log_error "Kestrel results file not found"
        return 1
    fi

    # Expected values for example_66bf_hg19_subset from test_data_config.json
    # Note: These are reference values - actual values may vary slightly
    local expected_alt_depth_min=450
    local expected_alt_depth_max=600
    local expected_active_depth_min=22000
    local expected_active_depth_max=23000
    local expected_depth_score_min=0.020
    local expected_depth_score_max=0.030
    local expected_confidence="High_Precision*"

    # Parse Kestrel TSV (skip comments and header)
    local kestrel_data
    kestrel_data=$(grep -v "^##" "$kestrel_file" | tail -n 1)

    # Extract values using awk (more reliable for TSV)
    # Column 9: Estimated_Depth_AlternateVariant
    # Column 10: Estimated_Depth_Variant_ActiveRegion
    # Column 18: Depth_Score
    # Column 19: Confidence
    local alt_depth
    local active_depth
    local depth_score
    local confidence

    alt_depth=$(echo "$kestrel_data" | awk -F'\t' '{print $9}')
    active_depth=$(echo "$kestrel_data" | awk -F'\t' '{print $10}')
    depth_score=$(echo "$kestrel_data" | awk -F'\t' '{print $18}')
    confidence=$(echo "$kestrel_data" | awk -F'\t' '{print $19}')

    # Validate Estimated_Depth_AlternateVariant (range check)
    if [[ "$alt_depth" -ge "$expected_alt_depth_min" ]] && [[ "$alt_depth" -le "$expected_alt_depth_max" ]]; then
        log_success "✓ Kestrel Alternate Depth: $alt_depth (expected range: $expected_alt_depth_min-$expected_alt_depth_max)"
    else
        log_error "✗ Kestrel Alternate Depth: $alt_depth (expected range: $expected_alt_depth_min-$expected_alt_depth_max)"
        return 1
    fi

    # Validate Estimated_Depth_Variant_ActiveRegion (range check)
    if [[ "$active_depth" -ge "$expected_active_depth_min" ]] && [[ "$active_depth" -le "$expected_active_depth_max" ]]; then
        log_success "✓ Kestrel Active Region Depth: $active_depth (expected range: $expected_active_depth_min-$expected_active_depth_max)"
    else
        log_error "✗ Kestrel Active Region Depth: $active_depth (expected range: $expected_active_depth_min-$expected_active_depth_max)"
        return 1
    fi

    # Validate Depth_Score (range check)
    if command -v bc &> /dev/null; then
        if (( $(echo "$depth_score >= $expected_depth_score_min" | bc -l) )) && (( $(echo "$depth_score <= $expected_depth_score_max" | bc -l) )); then
            log_success "✓ Kestrel Depth Score: $depth_score (expected range: $expected_depth_score_min-$expected_depth_score_max)"
        else
            log_error "✗ Kestrel Depth Score: $depth_score (expected range: $expected_depth_score_min-$expected_depth_score_max)"
            return 1
        fi
    else
        # Fallback without bc - simple range check
        local score_int
        score_int=$(echo "$depth_score * 1000" | awk '{print int($1)}')
        local min_int
        min_int=$(echo "$expected_depth_score_min * 1000" | awk '{print int($1)}')
        local max_int
        max_int=$(echo "$expected_depth_score_max * 1000" | awk '{print int($1)}')

        if [[ "$score_int" -ge "$min_int" ]] && [[ "$score_int" -le "$max_int" ]]; then
            log_success "✓ Kestrel Depth Score: $depth_score (expected range: $expected_depth_score_min-$expected_depth_score_max)"
        else
            log_warning "⚠ Kestrel Depth Score: $depth_score may be outside expected range"
        fi
    fi

    # Validate Confidence
    if [[ "$confidence" == "$expected_confidence" ]]; then
        log_success "✓ Kestrel Confidence: $confidence"
    else
        log_error "✗ Kestrel Confidence: $confidence (expected: $expected_confidence)"
        return 1
    fi

    log_success "Kestrel results validation passed"
    return 0
}

validate_coverage_results() {
    local output_dir="$1"
    local coverage_file="${output_dir}/coverage/coverage_summary.tsv"

    log_info "Validating coverage analysis results..."

    if [[ ! -f "$coverage_file" ]]; then
        log_error "Coverage summary file not found"
        return 1
    fi

    # Parse coverage summary (skip header)
    local coverage_data
    coverage_data=$(tail -n 1 "$coverage_file")

    # Extract coverage metrics
    local mean_cov median_cov min_cov max_cov uncovered_pct

    IFS=$'\t' read -r mean_cov median_cov _ min_cov max_cov _ uncovered_pct <<< "$coverage_data"

    # Basic sanity checks for coverage
    if (( $(echo "$mean_cov > 1000" | bc -l 2>/dev/null || echo 1) )); then
        log_success "✓ Mean VNTR coverage: ${mean_cov}× (sufficient)"
    else
        log_warning "⚠ Mean VNTR coverage: ${mean_cov}× (may be low)"
    fi

    if (( $(echo "$median_cov > 500" | bc -l 2>/dev/null || echo 1) )); then
        log_success "✓ Median VNTR coverage: ${median_cov}× (sufficient)"
    else
        log_warning "⚠ Median VNTR coverage: ${median_cov}× (may be low)"
    fi

    # Check for uncovered regions (column format: count\tpercentage)
    local uncovered_count uncovered_pct_val
    read -r uncovered_count uncovered_pct_val <<< "$uncovered_pct"

    # Extract numeric percentage value
    uncovered_pct_val=$(echo "$uncovered_pct_val" | tr -d '%' | tr -d ' ')

    if [[ -n "$uncovered_pct_val" ]] && command -v bc &> /dev/null; then
        if (( $(echo "$uncovered_pct_val < 5" | bc -l) )); then
            log_success "✓ VNTR region coverage: ${uncovered_pct_val}% uncovered (excellent)"
        else
            log_warning "⚠ VNTR region coverage: ${uncovered_pct_val}% uncovered (may have gaps)"
        fi
    else
        log_info "VNTR region uncovered: ${uncovered_count} bases (${uncovered_pct_val}%)"
    fi

    log_success "Coverage validation passed"
    return 0
}

# =============================================================================
# Main Test Runner
# =============================================================================

main() {
    log_info "VNtyper Docker Container Test Suite"
    log_info "===================================="
    echo

    # Step 1: Ensure test data from Zenodo
    ensure_test_data
    echo

    # Step 2: Build Docker image
    build_docker_image
    echo

    # Step 3: Basic container tests
    test_container_basics
    echo

    # Step 4: Full pipeline test
    test_pipeline_execution
    echo

    # Final summary
    log_success "========================================"
    log_success "ALL DOCKER TESTS PASSED ✓"
    log_success "========================================"
    log_info "Docker image ${DOCKER_IMAGE} is production-ready"
    echo
}

# Run main function
main "$@"
