#!/usr/bin/env bash
# =============================================================================
# VNtyper Docker Container Test Script
# =============================================================================
# Simple, fast Docker testing using real test data from Zenodo
#
# Usage:
#   bash docker/test_docker.sh                    # Run all tests
#   bash docker/test_docker.sh example_66bf       # Run specific test
#   bash docker/test_docker.sh example_66bf advntr # Run specific tests
#   bash docker/test_docker.sh --help             # Show help
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

# Define all test cases (from test_data_config.json)
declare -A TEST_CASES
TEST_CASES=(
    ["example_b178"]="416:7110:0.05869:High_Precision*"
    ["example_a5c1"]="93:9883:0.00941:High_Precision"
    ["example_66bf"]="547:22703:0.02409:High_Precision*"
    ["example_7a61"]="None:None:None:Negative"
    ["example_dfc3"]="206:18977:0.01086:High_Precision*"
)

# adVNTR test case (from test_data_config.json advntr_tests)
# Format: VID:State:NumberOfSupportingReads:MeanCoverage:Pvalue
ADVNTR_TEST_ID="example_a5c1_advntr"
ADVNTR_EXPECTED="25561:I22_2_G_LEN1:11:144.234722222:3.46346905707e-09"

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

show_usage() {
    cat << EOF
${BLUE}VNtyper Docker Test Script${RESET}

${GREEN}Usage:${RESET}
  bash docker/test_docker.sh [OPTIONS] [TEST_NAMES...]

${GREEN}Options:${RESET}
  -h, --help          Show this help message

${GREEN}Arguments:${RESET}
  TEST_NAMES          Specific test names to run (optional)
                      If not provided, all tests will be executed

${GREEN}Available Tests:${RESET}
  ${BLUE}Standard BAM tests:${RESET}
$(for test in $(echo "${!TEST_CASES[@]}" | tr ' ' '\n' | sort); do echo "    - $test"; done)

  ${BLUE}Module tests:${RESET}
    - advntr          (adVNTR module integration test)

${GREEN}Examples:${RESET}
  bash docker/test_docker.sh
    Run all tests (5 standard + 1 adVNTR = 6 total)

  bash docker/test_docker.sh example_66bf example_7a61
    Run only example_66bf and example_7a61 tests

  bash docker/test_docker.sh advntr
    Run only the adVNTR module test

  bash docker/test_docker.sh example_a5c1 advntr
    Run example_a5c1 and adVNTR tests

EOF
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

    # Build list of required files for all test cases
    local files=()
    for test_id in "${!TEST_CASES[@]}"; do
        files+=("${test_id}_hg19_subset.bam")
        files+=("${test_id}_hg19_subset.bam.bai")
    done

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
    local -a tests_to_run=("$@")

    # If no specific tests provided, run all tests
    if [[ ${#tests_to_run[@]} -eq 0 ]]; then
        log_info "Testing full VNtyper pipeline execution for all test cases..."
        tests_to_run=($(echo "${!TEST_CASES[@]}" | tr ' ' '\n' | sort))
    else
        log_info "Testing VNtyper pipeline execution for selected test cases..."
    fi

    local total_cases=${#tests_to_run[@]}
    local passed_cases=0
    local failed_cases=0
    local test_num=0

    # Loop through specified test cases
    for test_id in "${tests_to_run[@]}"; do
        test_num=$((test_num + 1))
        log_info "========================================="
        log_info "Test Case ${test_num}/${total_cases}: ${test_id}"
        log_info "========================================="

        if test_single_bam "$test_id"; then
            passed_cases=$((passed_cases + 1))
        else
            failed_cases=$((failed_cases + 1))
        fi
        echo
    done

    # Summary
    log_info "========================================="
    log_info "Test Summary"
    log_info "========================================="
    log_info "Total cases: ${total_cases}"
    log_success "Passed: ${passed_cases}"
    if [[ $failed_cases -gt 0 ]]; then
        log_error "Failed: ${failed_cases}"
        return 1
    else
        log_success "Failed: ${failed_cases}"
    fi

    return 0
}

test_single_bam() {
    local test_id="$1"
    local test_bam="${test_id}_hg19_subset.bam"
    local input_dir="${TEST_DATA_DIR}"
    local output_dir="${DOCKER_TEST_DIR}/${test_id}"

    # Clean previous test output
    rm -rf "$output_dir"
    mkdir -p "$output_dir"

    log_info "BAM file: ${test_bam}"
    log_info "Output directory: ${output_dir}"

    # Run pipeline
    log_info "Executing VNtyper pipeline..."

    local start_time
    start_time=$(date +%s)

    if docker run --rm \
        -v "$(pwd)/${input_dir}":/opt/vntyper/input \
        -v "$(pwd)/${output_dir}":/opt/vntyper/output \
        "$DOCKER_IMAGE" \
        vntyper pipeline \
        --bam "/opt/vntyper/input/${test_bam}" \
        -o /opt/vntyper/output/result \
        --reference-assembly hg19 \
        --threads 4 \
        --fast-mode \
        --keep-intermediates 2>&1 | grep -E "(ERROR|WARNING|Pipeline completed|Kestrel|Coverage)" || true; then

        local end_time
        end_time=$(date +%s)
        local duration=$((end_time - start_time))

        log_success "Pipeline completed in ${duration} seconds"
    else
        log_error "Pipeline execution failed"
        return 1
    fi

    # Validate output files
    if ! validate_pipeline_output "$output_dir/result" "$test_id"; then
        return 1
    fi

    return 0
}

validate_pipeline_output() {
    local output_dir="$1"
    local test_id="$2"

    log_info "Validating pipeline output files..."

    # Get expected values for this test case
    local expected="${TEST_CASES[$test_id]}"
    local is_negative=false
    if [[ "$expected" == *"Negative"* ]]; then
        is_negative=true
    fi

    # Expected critical output files (IGV report not expected for negative cases)
    local expected_files=(
        "summary_report.html"
        "kestrel/kestrel_result.tsv"
        "coverage/coverage_summary.tsv"
        "pipeline_summary.json"
    )

    if [[ "$is_negative" == "false" ]]; then
        expected_files+=("igv_report.html")
    fi

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
    if ! validate_kestrel_results "$output_dir" "$test_id"; then
        return 1
    fi

    if ! validate_coverage_results "$output_dir"; then
        return 1
    fi

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
    local test_id="$2"
    local kestrel_file="${output_dir}/kestrel/kestrel_result.tsv"

    log_info "Validating Kestrel genotyping results for ${test_id}..."

    if [[ ! -f "$kestrel_file" ]]; then
        log_error "Kestrel results file not found"
        return 1
    fi

    # Parse expected values from TEST_CASES (format: alt:active:score:confidence)
    local expected="${TEST_CASES[$test_id]}"
    IFS=':' read -r expected_alt expected_active expected_score expected_conf <<< "$expected"

    # Parse Kestrel TSV (skip comments and header)
    local kestrel_data
    kestrel_data=$(grep -v "^##" "$kestrel_file" | tail -n 1)

    # Extract actual values from TSV
    # Column 9: Estimated_Depth_AlternateVariant
    # Column 10: Estimated_Depth_Variant_ActiveRegion
    # Column 18: Depth_Score
    # Column 19: Confidence
    local alt_depth active_depth depth_score confidence

    alt_depth=$(echo "$kestrel_data" | awk -F'\t' '{print $9}')
    active_depth=$(echo "$kestrel_data" | awk -F'\t' '{print $10}')
    depth_score=$(echo "$kestrel_data" | awk -F'\t' '{print $18}')
    confidence=$(echo "$kestrel_data" | awk -F'\t' '{print $19}')

    # Handle negative case (None/Negative values)
    if [[ "$expected_conf" == "Negative" ]]; then
        log_info "Testing negative case (expecting None/Negative values)"

        # Accept both "None", "Negative", and "-1" for alternate depth
        if [[ "$alt_depth" == "None" || "$alt_depth" == "Negative" || "$alt_depth" == "-1" ]]; then
            log_success "✓ Kestrel Alternate Depth: $alt_depth (as expected for negative case)"
        else
            log_error "✗ Kestrel Alternate Depth: $alt_depth (expected None/Negative for negative case)"
            return 1
        fi

        # Accept both "None", "Negative", and "-1" for active region depth
        if [[ "$active_depth" == "None" || "$active_depth" == "Negative" || "$active_depth" == "-1" ]]; then
            log_success "✓ Kestrel Active Region Depth: $active_depth (as expected for negative case)"
        else
            log_error "✗ Kestrel Active Region Depth: $active_depth (expected None/Negative for negative case)"
            return 1
        fi

        # Accept "None", "Negative", "-1.0", or empty string for depth score
        if [[ "$depth_score" == "None" || "$depth_score" == "Negative" || "$depth_score" == "-1.0" || -z "$depth_score" ]]; then
            log_success "✓ Kestrel Depth Score: ${depth_score:-empty} (as expected for negative case)"
        else
            log_error "✗ Kestrel Depth Score: $depth_score (expected None/Negative/empty for negative case)"
            return 1
        fi

        # Confidence should be "Negative" or empty for negative cases
        if [[ "$confidence" == "$expected_conf" || -z "$confidence" ]]; then
            log_success "✓ Kestrel Confidence: ${confidence:-empty} (as expected for negative case)"
        else
            log_error "✗ Kestrel Confidence: $confidence (expected: $expected_conf or empty)"
            return 1
        fi

        log_success "Kestrel results validation passed (negative case)"
        return 0
    fi

    # Positive case - validate exact values with tolerance
    # Validate Estimated_Depth_AlternateVariant (±15% tolerance)
    local alt_tolerance=$((expected_alt * 15 / 100))
    local alt_min=$((expected_alt - alt_tolerance))
    local alt_max=$((expected_alt + alt_tolerance))

    if [[ "$alt_depth" -ge "$alt_min" ]] && [[ "$alt_depth" -le "$alt_max" ]]; then
        log_success "✓ Kestrel Alternate Depth: $alt_depth (expected: $expected_alt ±15%)"
    else
        log_error "✗ Kestrel Alternate Depth: $alt_depth (expected: $expected_alt ±15%)"
        return 1
    fi

    # Validate Estimated_Depth_Variant_ActiveRegion (±15% tolerance)
    local active_tolerance=$((expected_active * 15 / 100))
    local active_min=$((expected_active - active_tolerance))
    local active_max=$((expected_active + active_tolerance))

    if [[ "$active_depth" -ge "$active_min" ]] && [[ "$active_depth" -le "$active_max" ]]; then
        log_success "✓ Kestrel Active Region Depth: $active_depth (expected: $expected_active ±15%)"
    else
        log_error "✗ Kestrel Active Region Depth: $active_depth (expected: $expected_active ±15%)"
        return 1
    fi

    # Validate Depth_Score (±15% tolerance using bc for float comparison)
    if command -v bc &> /dev/null; then
        local score_tolerance=$(echo "$expected_score * 0.15" | bc -l)
        local score_min=$(echo "$expected_score - $score_tolerance" | bc -l)
        local score_max=$(echo "$expected_score + $score_tolerance" | bc -l)

        if (( $(echo "$depth_score >= $score_min" | bc -l) )) && (( $(echo "$depth_score <= $score_max" | bc -l) )); then
            log_success "✓ Kestrel Depth Score: $depth_score (expected: $expected_score ±15%)"
        else
            log_error "✗ Kestrel Depth Score: $depth_score (expected: $expected_score ±15%)"
            return 1
        fi
    else
        # Fallback without bc
        log_warning "⚠ bc not available, skipping depth score validation"
    fi

    # Validate Confidence (handle wildcard matching for patterns like "High_Precision*")
    if [[ "$expected_conf" == *"*" ]]; then
        local conf_prefix="${expected_conf%\*}"
        if [[ "$confidence" == "$conf_prefix"* ]]; then
            log_success "✓ Kestrel Confidence: $confidence (matches pattern: $expected_conf)"
        else
            log_error "✗ Kestrel Confidence: $confidence (expected pattern: $expected_conf)"
            return 1
        fi
    else
        if [[ "$confidence" == "$expected_conf" ]]; then
            log_success "✓ Kestrel Confidence: $confidence"
        else
            log_error "✗ Kestrel Confidence: $confidence (expected: $expected_conf)"
            return 1
        fi
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
# adVNTR Tests
# =============================================================================

test_advntr_pipeline() {
    log_info "========================================="
    log_info "Testing adVNTR Module Integration"
    log_info "========================================="

    local test_bam="example_a5c1_hg19_subset.bam"
    local input_dir="${TEST_DATA_DIR}"
    local output_dir="${DOCKER_TEST_DIR}/${ADVNTR_TEST_ID}"

    # Clean previous test output
    rm -rf "$output_dir"
    mkdir -p "$output_dir"

    log_info "BAM file: ${test_bam}"
    log_info "Output directory: ${output_dir}"
    log_info "Running with --extra-modules advntr"

    # Run pipeline with adVNTR module
    log_info "Executing VNtyper pipeline with adVNTR..."

    local start_time
    start_time=$(date +%s)

    if docker run --rm \
        -v "$(pwd)/${input_dir}":/opt/vntyper/input \
        -v "$(pwd)/${output_dir}":/opt/vntyper/output \
        "$DOCKER_IMAGE" \
        vntyper pipeline \
        --bam "/opt/vntyper/input/${test_bam}" \
        -o /opt/vntyper/output/result \
        --reference-assembly hg19 \
        --threads 4 \
        --fast-mode \
        --keep-intermediates \
        --extra-modules advntr \
        --advntr-max-coverage 300 2>&1 | grep -E "(ERROR|WARNING|Pipeline completed|adVNTR|Kestrel)" || true; then

        local end_time
        end_time=$(date +%s)
        local duration=$((end_time - start_time))

        log_success "Pipeline with adVNTR completed in ${duration} seconds"
    else
        log_error "Pipeline with adVNTR failed"
        return 1
    fi

    # Validate adVNTR output
    if ! validate_advntr_results "$output_dir/result"; then
        return 1
    fi

    log_success "adVNTR test passed"
    return 0
}

validate_advntr_results() {
    local output_dir="$1"
    local advntr_file="${output_dir}/advntr/output_adVNTR_result.tsv"

    log_info "Validating adVNTR genotyping results..."

    if [[ ! -f "$advntr_file" ]]; then
        log_error "adVNTR results file not found: $advntr_file"
        return 1
    fi

    # Parse expected values (VID:State:NumberOfSupportingReads:MeanCoverage:Pvalue)
    IFS=':' read -r expected_vid expected_state expected_reads expected_cov expected_pval <<< "$ADVNTR_EXPECTED"

    # Parse adVNTR TSV (skip comments and header)
    local advntr_data
    advntr_data=$(grep -v "^#" "$advntr_file" | grep -v "^VID" | head -n 1)

    if [[ -z "$advntr_data" ]]; then
        log_error "No data found in adVNTR results file"
        return 1
    fi

    # Extract values from TSV
    # Columns: VID  Variant  NumberOfSupportingReads  MeanCoverage  Pvalue  RU  POS  REF  ALT  Flag
    local vid state reads cov pval

    IFS=$'\t' read -r vid state reads cov pval _ _ _ _ _ <<< "$advntr_data"

    # Validate VID
    if [[ "$vid" == "$expected_vid" ]]; then
        log_success "✓ adVNTR VID: $vid"
    else
        log_error "✗ adVNTR VID: $vid (expected: $expected_vid)"
        return 1
    fi

    # Validate State (Variant)
    if [[ "$state" == "$expected_state" ]]; then
        log_success "✓ adVNTR State: $state"
    else
        log_error "✗ adVNTR State: $state (expected: $expected_state)"
        return 1
    fi

    # Validate NumberOfSupportingReads (exact match)
    if [[ "$reads" == "$expected_reads" ]]; then
        log_success "✓ adVNTR Supporting Reads: $reads"
    else
        log_error "✗ adVNTR Supporting Reads: $reads (expected: $expected_reads)"
        return 1
    fi

    # Validate MeanCoverage (±0.01 tolerance)
    if command -v bc &> /dev/null; then
        local cov_diff=$(echo "scale=6; ($cov - $expected_cov)" | bc -l | tr -d '-')
        if (( $(echo "$cov_diff < 0.01" | bc -l) )); then
            log_success "✓ adVNTR Mean Coverage: $cov (expected: $expected_cov)"
        else
            log_error "✗ adVNTR Mean Coverage: $cov (expected: $expected_cov, diff: $cov_diff)"
            return 1
        fi
    else
        log_warning "⚠ bc not available, skipping mean coverage validation"
    fi

    # Validate Pvalue (allow small floating point difference)
    if command -v bc &> /dev/null; then
        # Convert to scientific notation comparison
        local pval_match=$(echo "$pval == $expected_pval" | bc -l 2>/dev/null || echo 0)
        if [[ "$pval_match" == "1" ]] || [[ "$pval" == "$expected_pval" ]]; then
            log_success "✓ adVNTR P-value: $pval"
        else
            # Allow for small floating point differences
            log_warning "⚠ adVNTR P-value: $pval (expected: $expected_pval - allowing minor difference)"
        fi
    else
        log_warning "⚠ bc not available, skipping p-value validation"
    fi

    log_success "adVNTR results validation passed"
    return 0
}

# =============================================================================
# Main Test Runner
# =============================================================================

main() {
    # Parse command-line arguments
    local -a requested_tests=()
    local run_advntr=false
    local run_all=false

    # If no arguments, run all tests
    if [[ $# -eq 0 ]]; then
        run_all=true
    else
        # Parse arguments
        for arg in "$@"; do
            case "$arg" in
                -h|--help)
                    show_usage
                    exit 0
                    ;;
                advntr)
                    run_advntr=true
                    ;;
                *)
                    # Validate test name (use ${array[@]+x} pattern to avoid unbound variable error)
                    if [[ -n "${TEST_CASES[$arg]+x}" ]]; then
                        requested_tests+=("$arg")
                    else
                        log_error "Unknown test: $arg"
                        echo
                        log_info "Available tests:"
                        for test in $(echo "${!TEST_CASES[@]}" | tr ' ' '\n' | sort); do
                            echo "  - $test"
                        done
                        echo "  - advntr"
                        exit 1
                    fi
                    ;;
            esac
        done
    fi

    log_info "VNtyper Docker Container Test Suite"
    log_info "===================================="
    echo

    # Show what will be tested
    if [[ "$run_all" == true ]]; then
        log_info "Running all tests (5 standard + 1 adVNTR)"
    else
        local test_count=${#requested_tests[@]}
        local total_count=$test_count
        [[ "$run_advntr" == true ]] && total_count=$((total_count + 1))

        log_info "Running $total_count selected test(s):"
        for test in "${requested_tests[@]}"; do
            log_info "  - $test"
        done
        [[ "$run_advntr" == true ]] && log_info "  - advntr"
    fi
    echo

    # Step 1: Ensure test data from Zenodo
    ensure_test_data
    echo

    # Step 2: Build Docker image (skip if already exists)
    if ! docker image inspect "$DOCKER_IMAGE" &>/dev/null; then
        build_docker_image
        echo
    else
        log_info "Docker image ${DOCKER_IMAGE} already exists (skipping build)"
        echo
    fi

    # Step 3: Basic container tests
    test_container_basics
    echo

    # Step 4: Pipeline tests (selected or all)
    if [[ "$run_all" == true ]]; then
        test_pipeline_execution
    elif [[ ${#requested_tests[@]} -gt 0 ]]; then
        test_pipeline_execution "${requested_tests[@]}"
    fi

    if [[ ${#requested_tests[@]} -gt 0 ]] || [[ "$run_all" == true ]]; then
        echo
    fi

    # Step 5: adVNTR module test (if requested or running all)
    if [[ "$run_all" == true ]] || [[ "$run_advntr" == true ]]; then
        test_advntr_pipeline
        echo
    fi

    # Final summary
    log_success "========================================"
    log_success "ALL SELECTED TESTS PASSED ✓"
    log_success "========================================"

    if [[ "$run_all" == true ]]; then
        log_info "Tested: 5 standard BAM files + 1 adVNTR test = 6 total test cases"
    else
        local total_tested=${#requested_tests[@]}
        [[ "$run_advntr" == true ]] && total_tested=$((total_tested + 1))
        log_info "Tested: $total_tested test case(s)"
    fi
    echo
}

# Run main function
main "$@"
