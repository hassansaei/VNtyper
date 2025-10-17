# VNtyper Makefile
# Standardized development commands

.PHONY: help install install-dev lint lint-stats format format-check test test-unit test-integration test-cov clean build

# Colors for output
BLUE := \033[0;34m
GREEN := \033[0;32m
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
.PHONY: all check

all: format lint test
	@echo "$(GREEN)✓ All checks passed$(RESET)"

check: format-check test
	@echo "$(GREEN)✓ All checks passed$(RESET)"
