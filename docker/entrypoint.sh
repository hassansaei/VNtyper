#!/bin/bash
# =============================================================================
# VNtyper Entrypoint Script
# =============================================================================
# Production-grade entrypoint with error handling, validation, and health checks
# =============================================================================

set -euo pipefail  # Exit on error, undefined variables, and pipe failures
IFS=$'\n\t'        # Set safer Internal Field Separator

# =============================================================================
# Configuration
# =============================================================================
readonly CONDA_ENV="vntyper"
readonly LOG_LEVEL="${LOG_LEVEL:-INFO}"

# Color codes for output
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly NC='\033[0m' # No Color

# =============================================================================
# Logging Functions
# =============================================================================
log_info() {
    echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*" >&2
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*" >&2
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*" >&2
}

log_fatal() {
    log_error "$*"
    exit 1
}

# =============================================================================
# Error Handler
# =============================================================================
error_handler() {
    local line_num=$1
    local bash_lineno=$2
    local last_command=$3
    log_fatal "Error on line ${line_num}: Command '${last_command}' failed with exit code $?"
}

trap 'error_handler ${LINENO} ${BASH_LINENO} "$BASH_COMMAND"' ERR

# =============================================================================
# Cleanup Handler
# =============================================================================
cleanup() {
    log_info "Cleaning up before exit..."
    # Add any cleanup tasks here
}

trap cleanup EXIT

# =============================================================================
# Validation Functions
# =============================================================================
check_conda_environment() {
    if [[ ! -d "/opt/conda/envs/${CONDA_ENV}" ]]; then
        log_fatal "Conda environment '${CONDA_ENV}' not found"
    fi
}

check_reference_files() {
    local ref_dir="${REFERENCE_DIR:-/opt/vntyper/reference}"

    if [[ ! -d "${ref_dir}" ]]; then
        log_warn "Reference directory not found: ${ref_dir}"
        return 0
    fi

    if [[ ! -f "${ref_dir}/All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa" ]]; then
        log_warn "Required reference file not found"
    fi
}

check_permissions() {
    local output_dir="${DEFAULT_OUTPUT_DIR:-/opt/vntyper/output}"

    if [[ -d "${output_dir}" ]] && [[ ! -w "${output_dir}" ]]; then
        log_warn "Output directory is not writable: ${output_dir}"
    fi
}

# =============================================================================
# Health Check Function
# =============================================================================
health_check() {
    log_info "Performing health check..."

    # Check conda
    if ! command -v conda &> /dev/null; then
        log_error "Conda command not found"
        return 1
    fi

    # Check vntyper command
    if ! conda run -n "${CONDA_ENV}" vntyper --version &> /dev/null; then
        log_error "VNtyper command failed"
        return 1
    fi

    log_info "Health check passed"
    return 0
}

# =============================================================================
# Startup Validation
# =============================================================================
startup_validation() {
    if [[ "${LOG_LEVEL}" == "INFO" ]] || [[ "${LOG_LEVEL}" == "DEBUG" ]]; then
        log_info "Starting VNtyper container..."
        log_info "User: $(whoami) (UID: $(id -u), GID: $(id -g))"
    fi

    check_conda_environment
    check_reference_files
    check_permissions
}

# =============================================================================
# Source Conda
# =============================================================================
source_conda() {
    if [[ -f "/opt/conda/etc/profile.d/conda.sh" ]]; then
        source /opt/conda/etc/profile.d/conda.sh
    else
        log_fatal "Conda initialization script not found"
    fi
}

# =============================================================================
# Main Execution
# =============================================================================
main() {
    # Source conda environment
    source_conda

    # Run startup validation
    startup_validation

    # Parse command
    local cmd="${1:-}"

    case "${cmd}" in
        vntyper)
            log_info "Running VNtyper command: $*"
            exec conda run -n "${CONDA_ENV}" "$@"
            ;;

        celery)
            log_info "Starting Celery worker..."
            exec conda run -n "${CONDA_ENV}" celery -A app.celery_app worker \
                --loglevel=info \
                --concurrency=1
            ;;

        beat)
            log_info "Starting Celery Beat..."
            exec conda run -n "${CONDA_ENV}" celery -A app.celery_app beat \
                --loglevel=info
            ;;

        health)
            # Health check endpoint
            health_check
            exit $?
            ;;

        bash|sh)
            log_info "Starting interactive shell..."
            exec /bin/bash
            ;;

        "")
            # Default: Start FastAPI web service
            log_info "Starting FastAPI web service on port 8000..."
            exec conda run -n "${CONDA_ENV}" uvicorn app.main:app \
                --host 0.0.0.0 \
                --port 8000 \
                --log-level info
            ;;

        *)
            log_error "Unknown command: ${cmd}"
            log_info "Available commands: vntyper, celery, beat, health, bash"
            exit 1
            ;;
    esac
}

# =============================================================================
# Entry Point
# =============================================================================
main "$@"
