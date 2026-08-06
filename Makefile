# VNtyper Makefile
# Standardized development commands

.PHONY: help install install-dev lint lint-stats format format-check type-check type-check-tests type-check-all download-test-data verify-test-data test test-unit test-fast test-unit-cov test-integration test-integration-parallel test-advntr test-cov test-quiet test-verbose test-docker test-docker-quick test-docker-fast check check-all check-full check-ci ci-local ci-local-docker ci-local-docs ci-local-uv lint-actions lint-docker coverage-report patch-coverage mutation mutation-render test-docker-smoke clean build docker-build docker-build-base docker-clean docs-install docs-serve docs-build docs-check docs-clean

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
	@echo "  make lint             - Run Ruff linter over vntyper/ docker/app/ tests/ scripts/"
	@echo "  make lint-stats       - Run Ruff linter with detailed statistics"
	@echo "  make format           - Auto-format the same paths with Ruff"
	@echo "  make format-check     - Check formatting and lint without changes"
	@echo "  make type-check       - Run mypy type checker on vntyper/ and docker/app/"
	@echo "  make type-check-tests - Run mypy type checker on tests"
	@echo ""
	@echo "$(GREEN)Testing:$(RESET)"
	@echo "  make download-test-data      - Download test data from Zenodo (1.1GB, ~10-30 min)"
	@echo "  make verify-test-data        - Verify test data exists and has correct checksums"
	@echo "  make test                    - Run all tests (needs test data + Docker)"
	@echo "  make test-unit               - Run unit tests only (fast, ~0.5s)"
	@echo "  make test-fast               - Unit tests, fail-fast, last-failed first"
	@echo "  make test-unit-cov           - Unit tests + coverage floor (CI gate)"
	@echo "  make patch-coverage          - Coverage of changed lines only (CI gate)"
	@echo "  make mutation                - Advisory mutation score (not gated)"
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
	@echo "  make ci-local-uv       - replicate CI's uv install path in a temp venv"
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

# Linting and formatting targets
#
# RUFF_PATHS is the whole of the repository's first-party Python. Keep the four ruff
# targets below reading from this one variable: they used to scope to
# `vntyper/ docker/app/`, which left tests/ and scripts/ - thousands of lines, and the
# code that decides whether everything else is correct - with no linter and no
# formatter on them at all. The three checks that gate a PR (`check`, `check-all`,
# `ci-local`) all reach ruff through these targets, so widening the variable widens
# every gate at once.
RUFF_PATHS := vntyper/ docker/app/ tests/ scripts/

lint:
	@echo "$(BLUE)Running Ruff linter...$(RESET)"
	ruff check $(RUFF_PATHS)
	@echo "$(GREEN)✓ Linting complete$(RESET)"

lint-stats:
	@echo "$(BLUE)Running Ruff linter with statistics...$(RESET)"
	ruff check $(RUFF_PATHS) --statistics
	@echo "$(GREEN)✓ Linting complete$(RESET)"

format:
	@echo "$(BLUE)Formatting code with Ruff...$(RESET)"
	ruff format $(RUFF_PATHS)
	@echo "$(BLUE)Applying auto-fixes...$(RESET)"
	ruff check $(RUFF_PATHS) --fix
	@echo "$(GREEN)✓ Formatting complete$(RESET)"

format-check:
	@echo "$(BLUE)Checking code formatting...$(RESET)"
	ruff format $(RUFF_PATHS) --check
	ruff check $(RUFF_PATHS)
	@echo "$(GREEN)✓ Format check complete$(RESET)"

# Type checking targets
type-check:
	@echo "$(BLUE)Running mypy type checker on vntyper package and web service...$(RESET)"
	mypy vntyper/ docker/app/
	@echo "$(GREEN)✓ Type checking complete$(RESET)"

type-check-tests:
	@echo "$(BLUE)Running mypy type checker on tests...$(RESET)"
	mypy tests/
	@echo "$(GREEN)✓ Type checking complete$(RESET)"

# This is the target CI's `typecheck` job runs, because `type-check` alone leaves
# everything under tests/ outside the type gate.
#
# It runs mypy twice rather than over one combined argument list. Passing docker/app/
# and tests/ to the same invocation puts the `app` package on mypy's search path for
# the tests too, so the `from app import ...` lines in tests/unit/web/ stop resolving
# to Any and a batch of findings in those tests' hand-rolled ASGI doubles lands on this
# target. Keeping the two runs separate leaves each suite checked against the same
# scope CI checks it against; wiring the web tests up to the real types is worthwhile
# but is its own change, and `[tool.mypy]` in pyproject.toml records exactly what it
# would cost and what has to happen first.
type-check-all: type-check
	@echo "$(BLUE)Running mypy type checker on tests...$(RESET)"
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

# `--cov-report=xml` costs ~0.1s and writes coverage.xml, which `patch-coverage` below
# consumes. Emitting it here rather than from a second pytest run keeps the patch gate
# free in CI: the same measurement feeds both the whole-repo floor and the patch gate.
test-unit-cov:
	@echo "$(BLUE)Running unit tests with coverage...$(RESET)"
	pytest -m unit tests/unit -o log_cli=false --cov --cov-report=term-missing --cov-report=xml
	@python scripts/coverage_gate.py --target $(COVERAGE_TARGET)
	@echo "$(GREEN)✓ Unit coverage complete$(RESET)"

coverage-report:
	@python scripts/coverage_gate.py --target $(COVERAGE_TARGET)

# Patch coverage: the share of the lines THIS BRANCH CHANGED that the unit tier executes.
#
# The floor above is a whole-repo ratchet. It moves slowly and it is an average, so a PR
# can add a hundred untested lines and still not drop it a single point - which is why
# AGENTS.md rule 1 ("touch a file, add tests for it") was unenforceable in CI. This gate
# scores only the changed lines, so an untested new function fails its own PR no matter
# what the repo total is doing. It is also the only realistic route to the 80% target:
# every PR from here lands at >= 80%, so the average climbs on its own.
#
# Deliberately NOT a ratchet and deliberately NOT the same number as the floor: this is
# a fixed bar on new work, set at the project's COVERAGE_TARGET.
#
# diff-cover's default `...` range notation scores against the MERGE BASE, not the tip of
# the base branch, so commits landing on main while a PR is open are never charged to that
# PR. Finding that merge base needs real history - a shallow clone has none, so the CI job
# that runs this sets `fetch-depth: 0` (actions/checkout defaults to a depth of 1).
#
# Failure modes, all of which must pass rather than fail:
#   * docs-only PR      - ci-tests.yml's `python` path filter skips the whole job.
#   * deletion-only PR  - no ADDED lines, so there is nothing to score. Reports 100%.
#   * tests-only PR     - tests/ is outside `[tool.coverage.run] source`, so it never
#                         appears in coverage.xml and contributes no measurable lines.
#   * no Python at all  - same: nothing measurable in the diff, reports 100%.
# In every one of those cases diff-cover reports "No lines with coverage information"
# and exits 0.
PATCH_COVERAGE_TARGET ?= 80
PATCH_COVERAGE_BASE ?= origin/main

patch-coverage:
	@test -f coverage.xml || { \
		echo "$(RED)coverage.xml not found. Run 'make test-unit-cov' first.$(RESET)"; \
		exit 1; \
	}
	@echo "$(BLUE)Checking coverage of changed lines against $(PATCH_COVERAGE_BASE)...$(RESET)"
	diff-cover coverage.xml \
		--compare-branch=$(PATCH_COVERAGE_BASE) \
		--fail-under=$(PATCH_COVERAGE_TARGET)
	@echo "$(GREEN)✓ Patch coverage >= $(PATCH_COVERAGE_TARGET)%$(RESET)"

# Advisory mutation score over the five pure-decision modules.
#
# ADVISORY. Nothing gates on it and nothing should: equivalent mutants have not been
# hand-classified, so the printed score understates the truth by an unknown margin.
# It answers the question coverage cannot - "would a test have NOTICED if this line
# were wrong?" - which is the question that matters for a tool whose failure mode is a
# silently wrong genotype rather than a crash. confidence_assignment.py once had 100%
# line coverage and a 21% mutation score.
#
# PYTHONDONTWRITEBYTECODE=1 IS LOAD-BEARING. DO NOT REMOVE IT.
#
# Mutations here are mostly byte-length preserving (`<` -> `>`, `and` -> `or`) and
# CPython validates a cached .pyc against the source's (mtime, size) pair, with mtime
# at one-second granularity. A mutant written in the same second as the file it
# replaces is therefore indistinguishable from it to the cache validator: the
# interpreter loads the stale .pyc and runs the UNMUTATED code, every mutant "survives"
# and the score is fiction. Two sweeps on this branch were exactly that before it was
# found. scripts/mutation_test.py additionally deletes every __pycache__ before it
# starts - the env var stops new caches, the deletion stops old ones, and both are
# needed. The full explanation is in that file's module docstring.
#
# While this runs, vntyper/scripts/*.py is REWRITTEN IN PLACE, so do not build, package
# or install from the tree: a docker build, pip install or python -m build started
# mid-sweep bakes a live mutant into its artefact. One image built this way crashed in
# the container with a pandas KeyError that read exactly like a production bug. The
# harness's finally-restore protects the repo, not anything already built from it.
# Check with `git diff --quiet -- vntyper/` immediately before and after such a step.
#
# Takes ~15-30 min: every mutant is a separate pytest run. Use --module to scope it.
mutation:
	@echo "$(BLUE)Running advisory mutation testing (not a gate)...$(RESET)"
	PYTHONDONTWRITEBYTECODE=1 python scripts/mutation_test.py \
		--output docs/development/mutation-testing.md \
		--results-json docs/development/mutation-results.json
	@echo "$(GREEN)✓ Mutation run complete (advisory - nothing gates on this)$(RESET)"

# Re-render the page from the last sweep's saved results. Use this after adding an
# EQUIVALENT_MUTANTS entry to scripts/mutation_test.py: classifying a survivor changes
# how the measurement is PRESENTED, not the measurement, so re-running the 15-30 min
# sweep to pick it up would be waste. Seconds instead.
mutation-render:
	@echo "$(BLUE)Re-rendering the mutation page from saved results...$(RESET)"
	python scripts/mutation_test.py \
		--render-only docs/development/mutation-results.json \
		--output docs/development/mutation-testing.md
	@echo "$(GREEN)✓ Mutation page re-rendered$(RESET)"

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

# The quotes around $(ACTIONLINT) below are load-bearing. When the lookup finds
# nothing the variable expands to "", and an unquoted expansion would leave the
# `then` branch as a bare `;` - a shell *syntax* error, raised while parsing the
# whole `if ... fi` compound, so it fires before the `[ -n ... ]` guard can select
# the container fallback. That made this target fail on precisely the machines the
# fallback exists for. tests/unit/test_makefile_recipes.py pins it.
lint-actions:
	@echo "$(BLUE)Linting GitHub Actions workflows...$(RESET)"
	@if [ -n "$(ACTIONLINT)" ]; then \
		"$(ACTIONLINT)"; \
	elif command -v docker >/dev/null 2>&1; then \
		docker run --rm -v "$(PWD):/repo" --workdir /repo rhysd/actionlint:latest -color; \
	else \
		echo "$(RED)actionlint not found and Docker unavailable.$(RESET)"; \
		echo "Install: go install github.com/rhysd/actionlint/cmd/actionlint@latest"; \
		exit 1; \
	fi
	@echo "$(GREEN)✓ Workflows valid$(RESET)"

# `docker build --check` needs BuildKit and Docker Engine >= 25. Detect it and fail
# with the actual remedy, rather than surfacing an opaque "unknown flag" error.
lint-docker:
	@echo "$(BLUE)Linting Dockerfiles (BuildKit check)...$(RESET)"
	@if ! docker build --help 2>/dev/null | grep -q -- '--check'; then \
		echo "$(RED)This Docker Engine does not support 'docker build --check'.$(RESET)"; \
		echo "  It needs BuildKit and Docker Engine 25 or newer; CI runs 28.x."; \
		echo "  Upgrade Docker, or run 'make docker-build' alone to skip this lint."; \
		exit 1; \
	fi
	@DOCKER_BUILDKIT=1 docker build --check -f docker/Dockerfile.base . >/dev/null
	@DOCKER_BUILDKIT=1 docker build --check -f docker/Dockerfile \
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

# Replicates the INSTALL path CI uses (astral-sh/setup-uv -> uv venv -> uv pip install)
# in a throwaway venv, then runs the same test command.
#
# This exists because `ci-local` runs in whatever environment you already have, so it
# cannot catch breakage in how CI *builds* its environment. It missed exactly that once:
# `uv pip install --system` fails on Ubuntu's PEP 668 interpreter
# ("error: The interpreter at /usr is externally managed") and every CI job died at the
# install step while every local check was green.
CI_LOCAL_VENV := .ci-local-venv

ci-local-uv:
	@command -v uv >/dev/null 2>&1 || { \
		echo "$(RED)uv not installed - cannot verify the CI install path.$(RESET)"; \
		echo "  Fix: curl -LsSf https://astral.sh/uv/install.sh | sh"; \
		exit 1; \
	}
	@echo "$(BLUE)Replicating the CI install path with uv...$(RESET)"
	@rm -rf $(CI_LOCAL_VENV)
	uv venv $(CI_LOCAL_VENV)
	VIRTUAL_ENV=$(PWD)/$(CI_LOCAL_VENV) uv pip install -e ".[dev]"
	@echo "$(BLUE)Running the CI test command in that environment...$(RESET)"
	VIRTUAL_ENV=$(PWD)/$(CI_LOCAL_VENV) PATH="$(PWD)/$(CI_LOCAL_VENV)/bin:$$PATH" \
		$(MAKE) --no-print-directory test-unit-cov
	@rm -rf $(CI_LOCAL_VENV)
	@echo "$(GREEN)✓ CI install path verified$(RESET)"

# Mirrors ci-tests.yml: lint -> typecheck -> unit tests + coverage -> docs, plus a
# from-scratch install exactly as CI builds it.
ci-local: lint-actions format-check lint type-check-all test-unit-cov patch-coverage ci-local-docs ci-local-uv
	@echo ""
	@echo "$(GREEN)========================================$(RESET)"
	@echo "$(GREEN)✓ Local CI parity checks all passed$(RESET)"
	@echo "$(GREEN)========================================$(RESET)"
	@echo "Verified here: workflow syntax, format, lint, mypy, unit tests + coverage, docs."
	@echo "NOT verified here (CI only):"
	@echo "  - the Python 3.9-3.12 matrix; this ran your interpreter only"
	@echo "  - the Docker image jobs -> run 'make ci-local-docker'"

# Files whose contents define the base image. CI hashes exactly this set; keep in sync
# with the hashFiles() list in docker-base.yml (tests/unit/test_version_consistency.py
# and test_workflow_consistency.py guard the pieces that can drift).
BASE_INPUTS := conda docker/Dockerfile.base docker/requirements-web.txt \
	vntyper/scripts/install_references.py vntyper/scripts/install_references_config.json \
	vntyper/dependencies/advntr reference .dockerignore

# Mirrors docker-build.yml. Needs a Docker daemon.
#
# VNTYPER_TEST_IMAGE is forced to the image just built - otherwise the smoke tier would
# fall back to its default tag and silently test a different, older image.
#
# If you have edited a base input, building on the published :latest base tests
# something CI will not build, so refuse rather than give false assurance.
ci-local-docker: lint-docker
	@if git rev-parse --verify -q origin/main >/dev/null 2>&1 && \
	    ! git diff --quiet origin/main -- $(BASE_INPUTS) 2>/dev/null; then \
		if [ "$(DOCKER_BASE_IMAGE)" = "ghcr.io/hassansaei/vntyper-base:latest" ]; then \
			echo "$(RED)You changed a base image input, but DOCKER_BASE_IMAGE still points at the published :latest.$(RESET)"; \
			echo "  CI will build a new base from your changes; this would not. Run:"; \
			echo "    make docker-build-base"; \
			echo "    make ci-local-docker DOCKER_BASE_IMAGE=vntyper-base:local"; \
			exit 1; \
		fi; \
		echo "$(BLUE)Base inputs changed; using $(DOCKER_BASE_IMAGE)$(RESET)"; \
	fi
	@$(MAKE) --no-print-directory docker-build
	@$(MAKE) --no-print-directory test-docker-smoke VNTYPER_TEST_IMAGE=$(DOCKER_IMAGE)
	@echo "$(BLUE)Running the Docker tier CI runs on PRs...$(RESET)"
	@$(MAKE) --no-print-directory test-docker-quick VNTYPER_TEST_IMAGE=$(DOCKER_IMAGE)
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
	@echo "$(BLUE)Running image structure smoke tests on $(or $(VNTYPER_TEST_IMAGE),vntyper:local)...$(RESET)"
	VNTYPER_TEST_IMAGE=$(or $(VNTYPER_TEST_IMAGE),vntyper:local) \
		pytest -m smoke tests/docker -o log_cli=false
	@echo "$(GREEN)✓ Image smoke tests complete$(RESET)"

# Everything except the `slow` tier. This is what runs on main: adVNTR's HMM
# genotyping is 15-25 min on a 2-core runner and dominates the suite, so it is
# exercised on a schedule instead of blocking every merge.
# Note the marker composition: repeated -m flags OVERRIDE rather than combine.
test-docker-fast:
	@echo "$(BLUE)Running Docker tests excluding the slow tier...$(RESET)"
	$(if $(VNTYPER_TEST_IMAGE),VNTYPER_TEST_IMAGE=$(VNTYPER_TEST_IMAGE)) \
		pytest -m "docker and not slow" -v
	@echo "$(GREEN)✓ Docker tests (fast tier) complete$(RESET)"

test-docker:
	@echo "$(BLUE)Running all Docker integration tests with testcontainers...$(RESET)"
	@echo "$(BLUE)Note: Requires Docker daemon running$(RESET)"
	@if ! python -c "import testcontainers" 2>/dev/null; then \
		echo "$(RED)Error: testcontainers not installed. Run: pip install -e .[dev]$(RESET)"; \
		exit 1; \
	fi
	$(if $(VNTYPER_TEST_IMAGE),VNTYPER_TEST_IMAGE=$(VNTYPER_TEST_IMAGE)) pytest -m docker -v
	@echo "$(GREEN)✓ Docker tests complete$(RESET)"

test-docker-quick:
	@echo "$(BLUE)Running Docker quick test (single test case + health checks)...$(RESET)"
	@echo "$(BLUE)Note: Requires Docker daemon running$(RESET)"
	@if ! python -c "import testcontainers" 2>/dev/null; then \
		echo "$(RED)Error: testcontainers not installed. Run: pip install -e .[dev]$(RESET)"; \
		exit 1; \
	fi
	$(if $(VNTYPER_TEST_IMAGE),VNTYPER_TEST_IMAGE=$(VNTYPER_TEST_IMAGE)) \
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
