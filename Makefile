# VNtyper Makefile
# Standardized development commands

.PHONY: help install install-dev lint lint-stats format format-check test test-unit test-integration test-cov clean build docker-build docker-test docker-scan docker-scan-critical docker-clean

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
	@echo "  make typecheck        - Run mypy type checker on vntyper package"
	@echo "  make typecheck-tests  - Run mypy type checker on tests"
	@echo ""
	@echo "$(GREEN)Testing:$(RESET)"
	@echo "  make test             - Run all tests"
	@echo "  make test-unit        - Run unit tests only"
	@echo "  make test-integration - Run integration tests only"
	@echo "  make test-cov         - Run tests with coverage report"
	@echo ""
	@echo "$(GREEN)Build & Maintenance:$(RESET)"
	@echo "  make clean            - Remove build artifacts and cache"
	@echo "  make build            - Build distribution packages"
	@echo ""
	@echo "$(GREEN)Docker:$(RESET)"
	@echo "  make docker-build          - Build multi-stage Docker image (production-ready)"
	@echo "  make docker-test           - Test Docker container with Zenodo data (fast, local)"
	@echo "                               Run specific tests: bash docker/test_docker.sh example_66bf advntr"
	@echo "  make docker-scan           - Scan Docker image for vulnerabilities (all severities)"
	@echo "  make docker-scan-critical  - Scan Docker image for CRITICAL vulnerabilities only"
	@echo "  make docker-clean          - Remove all VNtyper Docker images"
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

# Linting targets (Ruff replaces flake8)
lint:
	@echo "$(BLUE)Running Ruff linter...$(RESET)"
	ruff check vntyper/
	@echo "$(GREEN)✓ Linting complete$(RESET)"

lint-stats:
	@echo "$(BLUE)Running Ruff linter with statistics...$(RESET)"
	ruff check vntyper/ --statistics
	@echo "$(GREEN)✓ Linting complete$(RESET)"

# Formatting targets (Ruff replaces black)
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
typecheck:
	@echo "$(BLUE)Running mypy type checker on vntyper package...$(RESET)"
	mypy vntyper/ --ignore-missing-imports
	@echo "$(GREEN)✓ Type checking complete$(RESET)"

typecheck-tests:
	@echo "$(BLUE)Running mypy type checker on tests...$(RESET)"
	mypy tests/ --ignore-missing-imports
	@echo "$(GREEN)✓ Type checking complete$(RESET)"

typecheck-all:
	@echo "$(BLUE)Running mypy type checker on all code...$(RESET)"
	mypy vntyper/ tests/ --ignore-missing-imports
	@echo "$(GREEN)✓ Type checking complete$(RESET)"

# Testing targets
test:
	@echo "$(BLUE)Running all tests...$(RESET)"
	pytest
	@echo "$(GREEN)✓ Tests complete$(RESET)"

test-unit:
	@echo "$(BLUE)Running unit tests...$(RESET)"
	pytest -m unit
	@echo "$(GREEN)✓ Unit tests complete$(RESET)"

test-integration:
	@echo "$(BLUE)Running integration tests...$(RESET)"
	pytest -m integration
	@echo "$(GREEN)✓ Integration tests complete$(RESET)"

test-cov:
	@echo "$(BLUE)Running tests with coverage...$(RESET)"
	pytest --cov=vntyper --cov-report=html --cov-report=term
	@echo "$(GREEN)✓ Coverage report generated in htmlcov/$(RESET)"

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

all: format lint typecheck test
	@echo "$(GREEN)✓ All checks passed$(RESET)"

check: format-check typecheck test
	@echo "$(GREEN)✓ All checks passed$(RESET)"

check-all: format-check lint typecheck-all test
	@echo "$(GREEN)✓ All checks passed (full suite)$(RESET)"

# Docker targets
#Docker configuration
DOCKER_IMAGE_NAME := vntyper
DOCKER_IMAGE_TAG := latest
DOCKER_IMAGE := $(DOCKER_IMAGE_NAME):$(DOCKER_IMAGE_TAG)

docker-build:
	@echo "$(BLUE)Building Docker image (multi-stage production build)...$(RESET)"
	DOCKER_BUILDKIT=1 docker build -f docker/Dockerfile -t $(DOCKER_IMAGE) .
	@echo "$(GREEN)✓ Docker image built: $(DOCKER_IMAGE)$(RESET)"
	@echo "$(GREEN)✓ Image uses multi-stage build (35% smaller, more secure)$(RESET)"

docker-test:
	@echo "$(BLUE)Testing Docker container with Zenodo test data...$(RESET)"
	@if [ ! -f docker/test_docker.sh ]; then \
		echo "$(RED)Error: docker/test_docker.sh not found$(RESET)"; \
		exit 1; \
	fi
	@bash docker/test_docker.sh
	@echo "$(GREEN)✓ Docker tests complete$(RESET)"

docker-scan:
	@echo "$(BLUE)Scanning Docker image for vulnerabilities...$(RESET)"
	@if ! command -v trivy >/dev/null 2>&1; then \
		echo "$(RED)Error: trivy is not installed$(RESET)"; \
		echo "Install with: curl -sfL https://raw.githubusercontent.com/aquasecurity/trivy/main/contrib/install.sh | sh -s -- -b /usr/local/bin"; \
		exit 1; \
	fi
	@trivy image --severity LOW,MEDIUM,HIGH,CRITICAL $(DOCKER_IMAGE)
	@echo "$(GREEN)✓ Vulnerability scan complete$(RESET)"

docker-scan-critical:
	@echo "$(BLUE)Scanning Docker image for CRITICAL vulnerabilities...$(RESET)"
	@if ! command -v trivy >/dev/null 2>&1; then \
		echo "$(RED)Error: trivy is not installed$(RESET)"; \
		echo "Install with: curl -sfL https://raw.githubusercontent.com/aquasecurity/trivy/main/contrib/install.sh | sh -s -- -b /usr/local/bin"; \
		exit 1; \
	fi
	@trivy image --severity CRITICAL --exit-code 1 $(DOCKER_IMAGE)
	@echo "$(GREEN)✓ No CRITICAL vulnerabilities found$(RESET)"

docker-clean:
	@echo "$(BLUE)Removing VNtyper Docker images...$(RESET)"
	@docker images | grep '$(DOCKER_IMAGE_NAME)' | awk '{print $$3}' | xargs -r docker rmi -f || true
	@echo "$(GREEN)✓ Docker images removed$(RESET)"
