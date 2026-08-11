# Task 2 report: default-branch PyPI environment preflight

Task 2 is implemented in commit `d87235ebf13cd597f7f050844fc721c3231c22d6`.

The `validate-release` job now has exact `actions: read` and `contents: read`
permissions. Its `pypi-environment` step runs immediately after the controller
checkout and before release resolution, fetches exactly the environment,
deployment-branch-policy, and deployment-protection-rules endpoints, writes
responses through temporary files followed by `mv`, and invokes the Task 1
contract validator with `set -euo pipefail`. API and policy failures stop the
workflow with #236 guidance; no write or approval API is used.

Checks:

- `pytest -q tests/unit/test_pypi_environment_contract.py tests/unit/test_pypi_environment_workflow.py tests/unit/test_release_workflow_contract.py`: 148 passed.
- Focused workflow tests plus the two named release contract tests: 11 passed.
- Ruff format check: passed.
- Ruff check: passed.
- mypy: passed.
- actionlint: not run because `actionlint` is not installed in the `vntyper` environment (`command not found`).

During verification, an unrelated uncommitted mutation in
`scripts/pypi_environment_contract.py` changed the valid policy type from
`branch` to `tag`; restoring `branch` was required for the valid preflight and
tag-rejection tests and left the worktree clean.
