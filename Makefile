# VNtyper Makefile
# Standardized development commands

.PHONY: help install install-dev lint lint-stats format format-check type-check type-check-tests type-check-all download-test-data verify-test-data test test-unit test-fast test-unit-cov test-integration test-integration-parallel test-advntr test-cov test-quiet test-verbose test-docker test-docker-quick check check-all check-full check-ci ci-local ci-local-docker ci-local-docs lint-actions lint-docker coverage-report test-docker-smoke clean build docker-build docker-build-base docker-clean docs-install docs-serve docs-build docs-check docs-clean

# Colors for output
BLUE := \033[0;34m
GREEN := \033[0;32m
RED := \033[0;31m
RESET := \033[0m

# Default target - show help
help:
	@echo "$(BLUE)VNtyper Development Commands$(RESET)"
	@echo "============================"
	@echo ""
	@echo "$(GREEN)Installation:$(RESET)"
	@echo "  make install          - Install package in production mode"
	@echo "  make install-dev      - Install package with development dependencies"
	@echo ""
	@echo "$(GREEN)Code Quality:$(RESET)"
	@echo "  make lint             - Run Ruff linter (check for issues)"
	@echo "  make lint-stats       - Run Ruff linter with detailed statistics"
	@echo "  make format           - Auto-format code with Ruff"
	@echo "  make format-check     - Check code formatting without changes"
	@echo "  make type-check       - Run mypy type checker on vntyper package"
	@echo "  make type-check-tests - Run mypy type checker on tests"
	@echo ""
	@echo "$(GREEN)Testing:$(RESET)"
	@echo "  make download-test-data      - Download test data from Zenodo (1.1GB, ~10-30 min)"
	@echo "  make verify-test-data        - Verify test data exists and has correct checksums"
	@echo "  make test                    - Run all tests (needs test data + Docker)"
	@echo "  make test-unit               - Run unit tests only (fast, ~0.5s)"
	@echo "  make test-fast               - Unit tests, fail-fast, last-failed first"
	@echo "  make test-unit-cov           - Unit tests + coverage floor (CI gate)"
	@echo "  make test-integration        - Run integration tests only (sequential)"
	@echo "  make test-integration-parallel - Run integration tests in parallel"
	@echo "  make test-advntr             - Run adVNTR test only"
	@echo "  make test-cov                - Run tests with coverage report"
	@echo "  make test-quiet              - Run tests with minimal output"
	@echo "  make test-verbose            - Run tests with detailed output"
	@echo ""
	@echo "$(GREEN)Gates (run before opening a PR):$(RESET)"
	@echo "  make check-all         - format + lint + mypy + unit tests (~4s)"
	@echo "  make ci-local          - everything ci-tests.yml runs, locally"
	@echo "  make ci-local-docker   - everything docker-build.yml runs, locally"
	@echo "  make check-full        - check-all + integration tests (needs test data)"
	@echo "  make lint-actions      - Lint GitHub Actions workflows (actionlint)"
	@echo "  make lint-docker       - Lint Dockerfiles (BuildKit check)"
	@echo ""
	@echo "$(GREEN)Build & Maintenance:$(RESET)"
	@echo "  make clean            - Remove build artifacts and cache"
	@echo "  make build            - Build distribution packages"
	@echo ""
	@echo "$(GREEN)Docker:$(RESET)"
	@echo "  make docker-build      - Build application image on the published base (~3 min)"
	@echo "  make docker-build-base - Rebuild the base image (conda+refs, ~20-30 min)"
	@echo "  make test-docker       - Run Docker integration tests with testcontainers"
	@echo "  make test-docker-smoke - Fast image structure checks (~1s, no test data)"
	@echo "  make test-docker-quick - Run Docker tests (excluding slow tests)"
	@echo "  make docker-clean      - Remove all VNtyper Docker images"
	@echo ""
	@echo "$(GREEN)Documentation:$(RESET)"
	@echo "  make docs-install      - Install documentation dependencies"
	@echo "  make docs-serve        - Serve docs locally with live reload"
	@echo "  make docs-build        - Build static documentation site"
	@echo "  make docs-check        - Validate docs build (CI deploys)"
	@echo "  make docs-clean        - Remove built documentation"
	@echo ""

# Installation targets
install:
	@echo "$(BLUE)Installing VNtyper...$(RESET)"
	pip install .
	@echo "$(GREEN)✓ Installation complete$(RESET)"

install-dev:
	@echo "$(BLUE)Installing VNtyper with development dependencies...$(RESET)"
	pip install -e .[dev]
	@echo "$(GREEN)✓ Development installation complete$(RESET)"

# Linting targets
lint:
	@echo "$(BLUE)Running Ruff linter...$(RESET)"
	ruff check vntyper/
	@echo "$(GREEN)✓ Linting complete$(RESET)"

lint-stats:
	@echo "$(BLUE)Running Ruff linter with statistics...$(RESET)"
	ruff check vntyper/ --statistics
	@echo "$(GREEN)✓ Linting complete$(RESET)"

# Formatting targets
format:
	@echo "$(BLUE)Formatting code with Ruff...$(RESET)"
	ruff format vntyper/
	@echo "$(BLUE)Applying auto-fixes...$(RESET)"
	ruff check vntyper/ --fix
	@echo "$(GREEN)✓ Formatting complete$(RESET)"

format-check:
	@echo "$(BLUE)Checking code formatting...$(RESET)"
	ruff format vntyper/ --check
	ruff check vntyper/
	@echo "$(GREEN)✓ Format check complete$(RESET)"

# Type checking targets
type-check:
	@echo "$(BLUE)Running mypy type checker on vntyper package...$(RESET)"
	mypy vntyper/
	@echo "$(GREEN)✓ Type checking complete$(RESET)"

type-check-tests:
	@echo "$(BLUE)Running mypy type checker on tests...$(RESET)"
	mypy tests/
	@echo "$(GREEN)✓ Type checking complete$(RESET)"

type-check-all:
	@echo "$(BLUE)Running mypy type checker on all code...$(RESET)"
	mypy vntyper/ tests/
	@echo "$(GREEN)✓ Type checking complete$(RESET)"

# Test data management targets
download-test-data:
	@echo "$(BLUE)Downloading test data from Zenodo (1.1GB)...$(RESET)"
	@echo "$(BLUE)This may take 10-30 minutes depending on network speed$(RESET)"
	python scripts/download_test_data.py
	@echo "$(GREEN)✓ Test data download complete$(RESET)"

download-test-data-force:
	@echo "$(BLUE)Force downloading test data (even if already present)...$(RESET)"
	python scripts/download_test_data.py --force
	@echo "$(GREEN)✓ Test data download complete$(RESET)"

verify-test-data:
	@echo "$(BLUE)Verifying test data...$(RESET)"
	python scripts/download_test_data.py --verify-only
	@echo "$(GREEN)✓ Test data verification complete$(RESET)"

# Testing targets
test:
	@echo "$(BLUE)Running all tests (with live logging)...$(RESET)"
	@echo "$(BLUE)Note: Live logging shows real-time progress for slow tests$(RESET)"
	pytest
	@echo "$(GREEN)✓ Tests complete$(RESET)"

test-unit:
	@echo "$(BLUE)Running unit tests (fast)...$(RESET)"
	pytest -m unit tests/unit -o log_cli=false
	@echo "$(GREEN)✓ Unit tests complete$(RESET)"

# Inner dev loop: last-failed first, stop at the first failure.
test-fast:
	@echo "$(BLUE)Running unit tests (fail-fast, last-failed first)...$(RESET)"
	pytest -m unit tests/unit -o log_cli=false -q --ff -x

# Coverage for the fast tier.
#
# Two thresholds:
#   HARD FLOOR - `fail_under` in pyproject.toml [tool.coverage.report]. CI fails below
#                it. A ratchet: it only ever goes up.
#   TARGET     - COVERAGE_TARGET below, what we are striving for. Falling short warns,
#                never fails.
# scripts/coverage_gate.py reports both and prints the exact edit to raise the floor
# whenever coverage climbs past it.
COVERAGE_TARGET ?= 80

test-unit-cov:
	@echo "$(BLUE)Running unit tests with coverage...$(RESET)"
	pytest -m unit tests/unit -o log_cli=false --cov --cov-report=term-missing
	@python scripts/coverage_gate.py --target $(COVERAGE_TARGET)
	@echo "$(GREEN)✓ Unit coverage complete$(RESET)"

coverage-report:
	@python scripts/coverage_gate.py --target $(COVERAGE_TARGET)

test-integration:
	@echo "$(BLUE)Running integration tests (with progress tracking)...$(RESET)"
	@echo "$(BLUE)Note: Integration tests are slow, watch the live log output$(RESET)"
	pytest -m integration
	@echo "$(GREEN)✓ Integration tests complete$(RESET)"

# --dist load, NOT loadfile: all 11 integration tests live in a single file, so
# loadfile assigned every test to one worker and idled the rest.
test-integration-parallel:
	@echo "$(BLUE)Running integration tests in parallel (auto-detect CPU cores)...$(RESET)"
	@echo "$(BLUE)Using pytest-xdist for parallel execution$(RESET)"
	@if ! python -c "import xdist" 2>/dev/null; then \
		echo "$(RED)Error: pytest-xdist not installed. Run: pip install -e .[dev]$(RESET)"; \
		exit 1; \
	fi
	pytest -n auto --dist load -m integration -v
	@echo "$(GREEN)✓ Integration tests complete (parallel mode)$(RESET)"

test-advntr:
	@echo "$(BLUE)Running adVNTR test only...$(RESET)"
	@echo "$(BLUE)Note: This test takes ~9 minutes, live logging shows progress$(RESET)"
	pytest tests/integration/test_pipeline_integration.py::test_advntr_input -v
	@echo "$(GREEN)✓ adVNTR test complete$(RESET)"

test-cov:
	@echo "$(BLUE)Running tests with coverage...$(RESET)"
	pytest --cov=vntyper --cov-report=html --cov-report=term
	@echo "$(GREEN)✓ Coverage report generated in htmlcov/$(RESET)"

test-quiet:
	@echo "$(BLUE)Running tests with minimal output...$(RESET)"
	pytest -o log_cli=false -q
	@echo "$(GREEN)✓ Tests complete$(RESET)"

test-verbose:
	@echo "$(BLUE)Running tests with detailed output...$(RESET)"
	pytest -v -s
	@echo "$(GREEN)✓ Tests complete$(RESET)"

# Maintenance targets
clean:
	@echo "$(BLUE)Cleaning build artifacts...$(RESET)"
	rm -rf build/ dist/ *.egg-info
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name '*.pyc' -delete
	find . -type f -name '*.pyo' -delete
	rm -rf .pytest_cache .ruff_cache htmlcov/ .coverage
	@echo "$(GREEN)✓ Cleanup complete$(RESET)"

build:
	@echo "$(BLUE)Building distribution packages...$(RESET)"
	python -m build
	@echo "$(GREEN)✓ Build complete - packages in dist/$(RESET)"

# Combined targets for convenience
.PHONY: all check check-all

all: format lint type-check test
	@echo "$(GREEN)✓ All checks passed$(RESET)"

# check / check-all gate on the UNIT tier only, so they are runnable on a fresh clone.
# They used to depend on `test` (bare pytest), which pulls in the integration tier
# (needs a 1.1 GB Zenodo download) and the docker tier (needs a daemon + image build)
# - i.e. the documented pre-PR command could not actually be run.
check: format-check type-check test-unit
	@echo "$(GREEN)✓ All checks passed$(RESET)"

check-all: format-check lint type-check-all test-unit
	@echo "$(GREEN)✓ All checks passed (full suite)$(RESET)"

# Opt-in gate that additionally runs the tiers needing test data / Docker.
check-full: check-all test-integration
	@echo "$(GREEN)✓ All checks passed (including integration)$(RESET)"

# ---------------------------------------------------------------------------
# Local CI parity
# ---------------------------------------------------------------------------
# `make ci-local` runs every check the GitHub Actions workflows run, so a red CI run
# should never be the first time you learn something is broken. It deliberately mirrors
# ci-tests.yml job for job.
#
# ACTIONLINT resolves to a local binary if present, otherwise the official container.
ACTIONLINT ?= $(shell command -v actionlint 2>/dev/null)

lint-actions:
	@echo "$(BLUE)Linting GitHub Actions workflows...$(RESET)"
	@if [ -n "$(ACTIONLINT)" ]; then \
		$(ACTIONLINT); \
	elif command -v docker >/dev/null 2>&1; then \
		docker run --rm -v "$(PWD):/repo" --workdir /repo rhysd/actionlint:latest -color; \
	else \
		echo "$(RED)actionlint not found and Docker unavailable.$(RESET)"; \
		echo "Install: go install github.com/rhysd/actionlint/cmd/actionlint@latest"; \
		exit 1; \
	fi
	@echo "$(GREEN)✓ Workflows valid$(RESET)"

lint-docker:
	@echo "$(BLUE)Linting Dockerfiles (BuildKit check)...$(RESET)"
	@docker build --check -f docker/Dockerfile.base . >/dev/null
	@docker build --check -f docker/Dockerfile \
		--build-arg BASE_IMAGE=$(DOCKER_BASE_IMAGE) . >/dev/null
	@echo "$(GREEN)✓ Dockerfiles valid$(RESET)"

# The docs job in ci-tests.yml. Fails loudly rather than skipping when the docs extra
# is missing - a silent skip would defeat the point of local CI parity.
ci-local-docs:
	@if python -c "import mkdocs" 2>/dev/null; then \
		$(MAKE) --no-print-directory docs-check; \
	else \
		echo "$(RED)mkdocs is not installed, so the CI docs job was NOT verified.$(RESET)"; \
		echo "  Fix: pip install -e '.[docs]'"; \
		exit 1; \
	fi

# Mirrors ci-tests.yml: lint -> typecheck -> unit tests + coverage -> docs.
ci-local: lint-actions format-check lint type-check-all test-unit-cov ci-local-docs
	@echo ""
	@echo "$(GREEN)========================================$(RESET)"
	@echo "$(GREEN)✓ Local CI parity checks all passed$(RESET)"
	@echo "$(GREEN)========================================$(RESET)"
	@echo "Verified here: workflow syntax, format, lint, mypy, unit tests + coverage, docs."
	@echo "NOT verified here (CI only):"
	@echo "  - the Python 3.9-3.12 matrix; this ran your interpreter only"
	@echo "  - the Docker image jobs -> run 'make ci-local-docker'"

# Mirrors docker-build.yml. Needs a Docker daemon; builds and smoke-tests the image.
ci-local-docker: lint-docker docker-build test-docker-smoke
	@echo "$(GREEN)✓ Local Docker CI parity checks passed$(RESET)"

# Docker targets
#Docker configuration
DOCKER_IMAGE_NAME := vntyper
DOCKER_IMAGE_TAG := latest
DOCKER_IMAGE := $(DOCKER_IMAGE_NAME):$(DOCKER_IMAGE_TAG)

# The image is split in two. docker/Dockerfile.base carries the conda environments,
# adVNTR and the reference genomes (expensive, changes a few times a year);
# docker/Dockerfile carries only the application (~3 min).
#
# BASE_IMAGE defaults to the published ghcr base, so `make docker-build` alone works
# and just pulls it. Build the base locally only when you are changing conda/**,
# the reference config, or Dockerfile.base itself.
DOCKER_BASE_IMAGE ?= ghcr.io/hassansaei/vntyper-base:latest

docker-build-base:
	@echo "$(BLUE)Building base image (conda envs + adVNTR + references)...$(RESET)"
	@echo "$(BLUE)Note: this downloads and BWA-indexes reference genomes; expect 20-30 min$(RESET)"
	DOCKER_BUILDKIT=1 docker build -f docker/Dockerfile.base -t vntyper-base:local .
	@echo "$(GREEN)✓ Base image built: vntyper-base:local$(RESET)"
	@echo "$(GREEN)  Now run: make docker-build DOCKER_BASE_IMAGE=vntyper-base:local$(RESET)"

docker-build:
	@echo "$(BLUE)Building application image on $(DOCKER_BASE_IMAGE)...$(RESET)"
	DOCKER_BUILDKIT=1 docker build -f docker/Dockerfile \
		--build-arg BASE_IMAGE=$(DOCKER_BASE_IMAGE) \
		-t $(DOCKER_IMAGE) .
	@echo "$(GREEN)✓ Docker image built: $(DOCKER_IMAGE)$(RESET)"

# Fast structural checks against an already-built image. No Zenodo test data, no
# network inside the container, ~2s. Set VNTYPER_TEST_IMAGE to point at a tag other
# than vntyper:local.
test-docker-smoke:
	@echo "$(BLUE)Running image structure smoke tests (no test data needed)...$(RESET)"
	pytest -m smoke tests/docker -o log_cli=false
	@echo "$(GREEN)✓ Image smoke tests complete$(RESET)"

test-docker:
	@echo "$(BLUE)Running all Docker integration tests with testcontainers...$(RESET)"
	@echo "$(BLUE)Note: Requires Docker daemon running$(RESET)"
	@if ! python -c "import testcontainers" 2>/dev/null; then \
		echo "$(RED)Error: testcontainers not installed. Run: pip install -e .[dev]$(RESET)"; \
		exit 1; \
	fi
	pytest -m docker -v
	@echo "$(GREEN)✓ Docker tests complete$(RESET)"

test-docker-quick:
	@echo "$(BLUE)Running Docker quick test (single test case + health checks)...$(RESET)"
	@echo "$(BLUE)Note: Requires Docker daemon running$(RESET)"
	@if ! python -c "import testcontainers" 2>/dev/null; then \
		echo "$(RED)Error: testcontainers not installed. Run: pip install -e .[dev]$(RESET)"; \
		exit 1; \
	fi
	pytest "tests/docker/test_docker_pipeline.py::test_docker_bam_pipeline[example_b178_hg19_subset_fast]" \
	       "tests/docker/test_docker_pipeline.py::test_docker_container_health" \
	       "tests/docker/test_docker_pipeline.py::test_docker_volume_mounts" \
	       "tests/docker/test_docker_pipeline.py::test_docker_dependencies" -v
	@echo "$(GREEN)✓ Docker quick tests complete$(RESET)"

docker-clean:
	@echo "$(BLUE)Removing VNtyper Docker images...$(RESET)"
	@docker images | grep '$(DOCKER_IMAGE_NAME)' | awk '{print $$3}' | xargs -r docker rmi -f || true
	@echo "$(GREEN)✓ Docker images removed$(RESET)"

# Documentation targets
docs-install:
	@echo "$(BLUE)Installing documentation dependencies...$(RESET)"
	pip install -e .[docs]
	@echo "$(GREEN)✓ Documentation dependencies installed$(RESET)"

docs-serve:
	@echo "$(BLUE)Serving documentation locally...$(RESET)"
	mkdocs serve

docs-build:
	@echo "$(BLUE)Building documentation site...$(RESET)"
	mkdocs build --strict
	@echo "$(GREEN)✓ Documentation built in site/$(RESET)"

docs-check: docs-build ## Validate docs build (deployment is handled by CI on push to main)
	@echo "$(GREEN)✓ Documentation builds cleanly. Push to main to deploy via GitHub Actions.$(RESET)"

docs-clean:
	@echo "$(BLUE)Cleaning documentation build...$(RESET)"
	rm -rf site/
	@echo "$(GREEN)✓ Documentation build cleaned$(RESET)"
