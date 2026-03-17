# Contributing

Thank you for considering contributing to VNtyper 2. Your contributions help improve the tool for the entire ADTKD-MUC1 research community. This page explains how to report issues, suggest features, set up a development environment, and submit changes.

## Reporting Bugs

Before opening a new bug report, search the existing [GitHub Issues](https://github.com/hassansaei/VNtyper/issues) to check whether the problem has already been reported.

When filing a new issue, please include:

- A short, descriptive title.
- Steps to reproduce the problem.
- Expected behavior versus actual behavior.
- Relevant environment details (OS, Python version, Docker version if applicable).
- Any error messages, logs, or screenshots that may help.

## Suggesting Features

Feature ideas are welcome. Open a [new issue](https://github.com/hassansaei/VNtyper/issues/new) and describe:

- What the feature would do and why it is useful.
- Any relevant references, examples, or prior art.
- A suggested approach, if you have one in mind.

## Development Setup

### Prerequisites

- Python 3.9 or higher
- Git
- Make (for running project workflows)

### Getting Started

1. Fork the repository on GitHub, then clone your fork:

    ```bash
    git clone https://github.com/your-username/VNtyper.git
    cd VNtyper
    ```

2. Install VNtyper 2 in editable mode with development dependencies:

    ```bash
    pip install -e ".[dev]"
    ```

    Or use the Makefile shortcut:

    ```bash
    make install-dev
    ```

!!! tip "Virtual environments"
    It is strongly recommended to work inside a virtual environment (`venv`, `conda`, or `mamba`) to avoid dependency conflicts with other projects.

## Code Quality

VNtyper 2 uses **Ruff** for linting and formatting, and **mypy** for static type checking. The following Makefile targets are available:

| Command              | Purpose                                    |
|----------------------|--------------------------------------------|
| `make format`        | Auto-format code and apply lint fixes       |
| `make format-check`  | Check formatting without modifying files    |
| `make lint`          | Run the Ruff linter                        |
| `make type-check`    | Run mypy on `vntyper/`                     |

!!! tip "Run all checks at once"
    Use `make all` to run formatting, linting, type checking, and tests in a single command.

## Running Tests

```bash
# Unit tests (fast, no external dependencies)
make test-unit

# Integration tests (requires downloaded test data)
make test-integration

# All tests with coverage report
make test-cov
```

!!! tip "Downloading test data"
    Integration tests require approximately 1.1 GB of test data from Zenodo. Run `make download-test-data` to fetch it and `make verify-test-data` to confirm checksums.

## Commit Conventions

VNtyper 2 follows the [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) specification. Each commit message should have the form:

```
type(scope): subject

body (optional)

footer (optional)
```

Common types:

`feat`
:   A new feature.

`fix`
:   A bug fix.

`docs`
:   Documentation changes only.

`refactor`
:   Code restructuring without changing behavior.

`test`
:   Adding or updating tests.

`build`
:   Changes to the build system or dependencies.

Example:

```
feat(pipeline): add GRCh38 assembly support

Added reference region coordinates for GRCh38 and updated the
assembly normalization logic in reference_registry.py.

Closes #42
```

## Pull Request Process

1. Create a feature branch from `main`:

    ```bash
    git checkout -b feature/your-feature-name
    ```

2. Make your changes and commit following the conventions above.

3. Ensure all checks pass locally:

    ```bash
    make check
    ```

4. Push to your fork and open a pull request against `hassansaei/VNtyper:main`.

5. In the pull request description, provide a clear summary of the changes and reference any related issues (e.g., `Closes #123`).

6. Respond to review feedback and update your branch as needed.

## Community Guidelines

- Be respectful and constructive in all interactions.
- Provide helpful, specific feedback in code reviews.
- Follow the project's Code of Conduct.
