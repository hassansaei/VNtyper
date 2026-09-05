# Contributing

Thank you for contributing to VNtyper 2. Contributions help improve tools for the ADTKD-MUC1 research community. This guide explains how to report issues, propose enhancements, set up local development, and submit pull requests.

## Reporting Bugs

Before submitting a bug report, search existing [GitHub Issues](https://github.com/hassansaei/VNtyper/issues) to verify the defect has not been logged.

When filing an issue, provide:

- A clear, concise title.
- Exact steps to reproduce the behavior.
- Expected versus observed results.
- Environment details (OS, Python version, architecture, container runtime if applicable).
- Relevant command logs, tracebacks, or output excerpts.

## Suggesting Features

Feature proposals are welcomed. Open a [new issue](https://github.com/hassansaei/VNtyper/issues/new) describing:

- The proposed functionality and use case.
- Relevant scientific literature or pipeline context.
- Implementation ideas or API considerations.

## Development Setup

### Prerequisites

- Python 3.10 or higher (`requires-python = ">=3.10"`; CI evaluates 3.10 through 3.13)
- Git
- Make

Pipeline execution requires external bioinformatics binaries (bwa, samtools, fastp, bcftools, Java 11). Creating the provided conda environment provides all prerequisites:

```bash
mamba env create -f conda/environment_vntyper.yml
conda activate vntyper
```

### Installation

1. Fork and clone the repository:

    ```bash
    git clone https://github.com/your-username/VNtyper.git
    cd VNtyper
    ```

2. Install VNtyper in editable mode with development dependencies:

    ```bash
    make install-dev
    ```

    Alternatively:

    ```bash
    pip install -e ".[dev]"
    ```

## Code Quality and Static Checks

VNtyper 2 enforces formatting and linting with **Ruff** and static type checking with **mypy**.

| Command | Target Scope | Purpose |
| --- | --- | --- |
| `make format` | `vntyper/`, `docker/app/`, `tests/`, `scripts/`, `docs/` | Format code and apply safe autofixes |
| `make format-check` | `RUFF_PATHS` | Check formatting compliance |
| `make lint` | `RUFF_PATHS` | Run Ruff linter |
| `make type-check` | `vntyper/`, `docker/app/`, `scripts/` | Run mypy type checker |
| `make type-check-all` | `type-check`, then `vntyper/`, `tests/` | Complete typecheck run matching CI |

`RUFF_PATHS` in `Makefile` covers `vntyper/ docker/app/ tests/ scripts/ docs/`, matching `ruff format --check .` discovery.

!!! tip "Full pre-PR gate"
    Execute `make check-all` before opening a pull request to run formatting, linting, type checks, and unit tests in a single command.

## Running Tests

```bash
# Fast unit tests (no external data required)
make test-unit

# Unit tests with coverage floor enforcement
make test-unit-cov

# Integration tests (requires Zenodo test cohort)
make test-integration
```

!!! tip "Test data download"
    Integration tests require the 1.1 GB Zenodo test data archive. Retrieve it via `make download-test-data` and verify checksums with `make verify-test-data`.

### Test execution in CI

| Tier | Trigger | Workflow |
| --- | --- | --- |
| Unit (`make test-unit-cov`) | Pull requests and pushes to `main` | `ci-tests.yml` |
| Docker smoke and fast checks | Pull requests and pushes to `main` | `docker-build.yml` |
| Full Docker tier (incl. adVNTR) | Nightly schedule; manual dispatch | `docker-build.yml` |
| Integration (`make test-integration`) | Weekly schedule against `:main` | `scheduled-tests.yml` |

The integration tier runs on a weekly schedule rather than blocking pull requests. This schedule isolates external risks (such as Zenodo uptime, reference genome mirror changes, and conda channel shifts) from standard PR gates. When modifying end-to-end execution paths, run `make check-full` locally.

## Commit Conventions

VNtyper 2 adheres to [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/). Format commit messages as:

```text
type(scope): subject

body (optional)

footer (optional)
```

Standard types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation updates
- `refactor`: Refactoring without functional changes
- `test`: Adding or adjusting tests
- `build`: Packaging or dependency changes

Example:

```text
feat(pipeline): add GRCh38 assembly support

Register reference coordinates for GRCh38 and update assembly
normalization in reference_registry.py.

Closes #42
```

## Pull Request Process

1. Create a feature branch:

    ```bash
    git checkout -b feature/your-feature-name
    ```

2. Make targeted changes following commit conventions.
3. Verify all checks pass:

    ```bash
    make check-all
    ```

4. Push to your fork and submit a pull request against `hassansaei/VNtyper:main`.
5. Describe the problem solved and reference related issues (e.g. `Closes #123`).
6. Address review feedback constructively.

## Community Standards

- Maintain respectful and constructive communication.
- Focus code reviews on technical merit, correctness, and maintainability.
- Adhere to the project Code of Conduct.
