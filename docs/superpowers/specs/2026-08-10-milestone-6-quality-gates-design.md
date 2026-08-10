# Milestone 6 quality-gates design

**Status:** Approved implementation specification

**Date:** 2026-08-10

**Baseline:** `main` at `ebb15b26631242a3295607e4eda4c68f688cd9a2`

**Issues:** [#204](https://github.com/hassansaei/VNtyper/issues/204), [#211](https://github.com/hassansaei/VNtyper/issues/211)

**Milestone:** 5. Gates, harness and release automation (parallel)

## Decision summary

Bring every Python file under `scripts/` into the same static and dynamic quality gates as
the application without weakening any existing threshold. First close #204 by fixing the one
known annotation and widening the runtime mypy invocation. Then add pure unit tests until
`scripts/` alone reaches at least 88.00% branch-inclusive coverage. Only after that local
target is met may `scripts` be appended to `[tool.coverage.run].source`. The final combined
unit-tier result must remain at or above the existing 86.00% hard floor on every supported
Python version (3.10, 3.11, 3.12, and 3.13). The 80% changed-line gate remains fixed.

There are no coverage omissions. In particular, the two currently unexecuted modules,
`scripts/download_test_data.py` and `scripts/advntr_len_differential.py`, are first-class
members of both the tests and the final coverage denominator.

This resolves two review concerns that supersede older issue comments: the repository floor
has since risen from 80 to 86 and must not move, and an 80% scripts-only target no longer gives
enough margin for a combined 86% result. The accepted target is therefore at least 88% for
`scripts/`, verified independently before the coverage source changes.

The reconstructed baseline is 69.01% scripts-only branch-inclusive coverage: 3,817
measurable statement/branch units with approximately 1,183 uncovered. Reaching 88.00%
requires covering roughly 725 additional units before allowing for newly added code.
This is the milestone's largest implementation tranche and is scheduled before the
final BLE001 inventory, not as cleanup after that policy is pinned.

## Linked issues and acceptance criteria

### Issue #204: `scripts/` is linted but not type-checked by any gate

The issue's requested acceptance criteria, preserved faithfully, are:

1. Annotate `dir_counts` in `scripts/download_test_data.py` as a dictionary with string keys
   and integer values.
2. Add `scripts/` to the `type-check` mypy invocation.
3. Place a note beside `RUFF_PATHS` explaining that the runtime Ruff and mypy scopes are kept
   in parity, so widening either scope requires reviewing the other.

The issue discussion further fixes these implementation details:

- Preserve the deliberate two-invocation structure: `type-check` checks production/runtime
  code, while `type-check-all` adds `tests/` in a separate invocation.
- Do not add a `scripts.*` mypy override or a `mypy_path`; neither is required.
- Retire the obsolete AGENTS.md trap that says scripts are not type-checked, and update its
  cross-reference only after the gate is real.
- Verify the exact widened invocation succeeds without incremental-cache assistance.

### Issue #211: test `scripts/`, then put it under coverage

The live issue direction says to add tests before widening coverage, include all of `scripts/`,
and never conceal the two zero-coverage entry points with an omit list. This specification
retains that direction and updates its numerical acceptance criteria for the current baseline:

1. `scripts/` alone measures **at least 88.00%**, including branch arcs.
2. `scripts/download_test_data.py` and `scripts/advntr_len_differential.py` both execute in
   unit tests and report non-zero coverage; neither may be omitted.
3. Only after criteria 1 and 2 pass, append `scripts` to coverage `source`, after the existing
   `vntyper` and `docker/app` entries.
4. The combined unit-tier coverage is at least **86.00% on each of Python 3.10, 3.11, 3.12,
   and 3.13**. Passing on average or on only one interpreter is insufficient.
5. `[tool.coverage.report].fail_under` remains **86** and the Makefile's advisory
   `COVERAGE_TARGET` remains **86**. They have different semantics and must not be collapsed.
6. `PATCH_COVERAGE_TARGET` remains **80**. Once `scripts` enters `source`, changed Python lines
   there are scored by the same merge-base patch gate as other measured Python.
7. No file, package, pattern, generated branch, CLI entry point, or difficult I/O path under
   `scripts/` is omitted from measurement.
8. A regression test pins branch measurement and complete source scope, because deleting
   `branch = true` or removing a poorly covered source can make the headline percentage rise.

Older comments proposing `fail_under = 74`, later comments retaining an 80 floor, and their
associated measured totals are historical evidence, not the current acceptance threshold.
The authoritative floor at this baseline is 86.

## Scope

### In scope

- `scripts/download_test_data.py`: the missing `dir_counts` annotation and pure tests for
  argument validation, manifest/path decisions, checksum decisions, archive/resource handling,
  and exit policy with all network and filesystem effects mocked or confined to `tmp_path`.
- `scripts/advntr_len_differential.py`: pure tests for parsing, command construction, case
  classification, result comparison, summary/exit policy, and failure reporting.
- The quality-gates track creates `tests/unit/test_golden_cohort_launcher.py`; the later
  exception-policy track modifies that shared file rather than recreating it and preserves the
  single `golden_cohort.launcher` module identity.
- Every other Python file under `scripts/`, prioritizing the currently weakest command-line and
  subprocess drivers: `golden_cohort/launcher.py`, `golden_cohort_gate.py`,
  `golden_cohort/artifacts.py`, `golden_cohort/runner.py`, `make_cram_fixtures.py`, and
  `coverage_gate.py`.
- Focused extraction of pure decisions from an I/O-heavy script when direct testing would
  otherwise require a network, Docker, conda, reference data, or a real subprocess.
- `Makefile`, `pyproject.toml`, `pytest.ini` only where needed for the new gates and marker
  behavior.
- `pytest.ini` remains the sole live pytest authority with strict markers and all five marker names;
  the dead duplicated `[tool.pytest.ini_options]` block is removed rather than maintained in parallel.
- Unit tests and gate-invariant tests under `tests/unit/`.
- CI workflow commentary or commands necessary to make the 3.10-3.13 matrix exercise the final
  configuration.
- Contributor-facing gate documentation, AGENTS.md trap retirement, and MkDocs exclusion of
  this contributor design area.

### Out of scope

- Lowering or rounding down the 86% hard floor.
- Raising or lowering the independently fixed 80% patch threshold.
- Excluding either zero-coverage module or adding any coverage `omit` entry for `scripts/`.
- Network, Docker, reference-data, or full pipeline execution in the unit tier.
- Rewriting whole scripts merely to make them testable; extract only the pure region under test.
- Changing genotype, scoring, calibration, clinical-sounding messages, archive contents, test
  data, or external command semantics except where a characterization test exposes an existing
  defect that is separately authorized.
- Folding the tests mypy invocation into the runtime invocation.
- Adding blanket mypy ignores, weakening strictness, or making `scripts/` a Python package solely
  for type checking.
- Mutation-score gating. `make mutation` remains advisory.

## Considered approaches

### Approach A: widen coverage now and lower the floor

This is rejected. It makes the untested code visible quickly, but it weakens the current 86%
ratchet and contradicts the maintainer's direction to earn coverage first. It also makes a
known transient state part of the repository's gate contract.

### Approach B: omit entry points and measure only library-like helpers

This is rejected. The two zero-coverage entry points and the subprocess drivers are precisely
where policy, exit codes, and evidence production concentrate. An omit list would preserve the
gap the issue is intended to close, and the combined metric would cease to describe the whole
gate harness.

### Approach C: test first, independently prove 88%, then append the source

This is selected. It keeps the 86% floor intact, gives the combined suite margin across four
interpreters, makes both zero modules visible, and lets the final source-list change act only as
enforcement. Small pure-decision extraction follows the repository's established architecture:
policy is importable and unit-testable; orchestration retains filesystem and process effects.

## Architecture and component boundaries

### Static-analysis boundary

`Makefile` remains the single command interface:

- `type-check` runs mypy once over `vntyper/ docker/app/ scripts/`.
- `type-check-all` depends on `type-check`, then runs the existing second mypy invocation over
  `vntyper/ tests/`.
- `RUFF_PATHS` continues to include `vntyper/ docker/app/ tests/ scripts/ docs/`. Its nearby
  comment states that the production portions of Ruff and mypy scope must be reviewed together;
  tests and docs intentionally differ because tests are checked separately and docs are not a
  mypy target.

The `scripts/` directory remains importable for tests using its current layout. No package marker
or import-path override is introduced unless an existing import proves impossible on Python 3.10;
the baseline combined mypy run already demonstrates that it is unnecessary.

### Coverage boundary

Coverage has two explicit views:

1. A temporary scripts-only measurement used while tests are being built. It runs the unit tier
   from the repository root with branch measurement enabled, `--cov=scripts`, and an 88% failure
   threshold. Its data file lives in a temporary directory and is removed on success or failure,
   so it cannot be confused with the combined report.
2. The canonical combined measurement used by `make test-unit-cov`, `make check-all`, CI, and
   `make patch-coverage`. Its final source order is `vntyper`, `docker/app`, then `scripts`; its
   report floor remains 86.

The scripts-only target is an implementation and review proof, not a replacement for the
combined gate. The combined matrix is authoritative for merge readiness.

### Testability boundary

Tests call pure functions with explicit values and capture structured returns. Code that owns
network calls, subprocesses, archive extraction, environment mutation, or writes remains a thin
adapter. When extraction is required, the new pure function stays in the script's existing
module unless that file is already beyond the repository's size guideline and the extracted
region has a coherent independent responsibility. Every added or modified function receives
value assertions and failure-path coverage.

The unit tier uses `tmp_path`, `unittest.mock`, synthetic manifests/VCFs/TSVs, and deterministic
subprocess result objects. It never downloads test data or invokes conda, Docker, Java, samtools,
adVNTR, or VNtyper pipeline stages.

### Configuration-invariant boundary

Dedicated unit tests guard the configuration rather than trusting the reported percentage:

- coverage branch mode is exactly enabled;
- coverage source contains `vntyper`, `docker/app`, and `scripts`, in that order, with no scripts
  omission;
- hard floor and advisory target are both 86 and remain separately declared;
- patch target is 80;
- the runtime mypy scope contains `scripts/` while the separate test invocation remains present;
- every new test module is marked `pytest.mark.unit`.

These checks extend the existing `tests/unit/test_coverage_gate.py`,
`tests/unit/test_makefile_recipes.py`, and marker/version consistency tests rather than creating
a second configuration authority.

## Interfaces, inputs, outputs, and invariants

| Interface | Inputs | Outputs | Invariants |
| --- | --- | --- | --- |
| `make type-check` | repository Python sources, `[tool.mypy]` | success or non-zero mypy status | Runtime scope is exactly application, web app, and scripts; no scripts override hides errors. |
| `make type-check-all` | `type-check` result plus tests | success or non-zero mypy status | Test checking remains a second invocation. |
| scripts-only coverage target | unit selection, `scripts/`, temporary coverage data path | precise branch-inclusive total and per-file report | Total is at least 88.00%; all script files appear; both former zero modules are non-zero. |
| `make test-unit-cov` | unit tests, final coverage source | combined XML/terminal report and gate status | Each supported interpreter reaches at least 86.00%; branch mode is on. |
| `make patch-coverage` | combined `coverage.xml`, merge-base diff | changed-line percentage and gate status | Target remains 80; merge-base semantics remain intact; no measurable change is a pass. |
| scope guard tests | `Makefile`, `pyproject.toml`, `pytest.ini` | actionable assertion failures | Removing `scripts`, branch mode, a threshold, or the second mypy pass fails deterministically. |

Global invariants:

- The version floor is Python 3.10; implementation uses no 3.11+ syntax or standard-library API.
- All commands run from the repository root because test collection reads a relative data config.
- Coverage figures use statements plus branch arcs, never the rounded `TOTAL` column alone.
- The exact percentage emitted by the coverage gate is used for evidence; thresholds are not
  derived from rounded output.
- Tests prove values and decisions, not merely that a function returned.

## Dependency and implementation ordering

### Phase 1: freeze the current contracts

1. Re-measure the baseline with a clean coverage data file and record the precise combined total,
   scripts-only total, per-file totals, and mypy error set.
2. Add or update invariant tests for the current 86 hard floor, 86 advisory target, 80 patch
   target, branch mode, and two mypy invocations. The complete-source assertion lands with the
   source append in Phase 4 so every intermediate commit remains green.
3. Add `superpowers/` to `mkdocs.yml`'s `exclude_docs` before the implementation PR's strict docs
   build, because these files are contributor working specifications rather than published user
   documentation.

### Phase 2: close the static-analysis gap (#204)

1. Add `dir_counts: dict[str, int] = {}` at the existing declaration.
2. Add `scripts/` to `type-check`, preserving `type-check-all`'s second invocation.
3. Add the Ruff/mypy scope-parity explanation and update command help.
4. Run no-cache mypy on the widened runtime scope and the complete two-step target.
5. Retire the obsolete AGENTS.md trap only after both commands pass.

### Phase 3: earn scripts coverage while scripts is not yet a configured source

1. Add unit tests for `download_test_data.py` and `advntr_len_differential.py` first. Each file
   must move above zero before work proceeds to aggregate optimization.
2. Add decision- and failure-path tests for the weakest harness modules in measured priority
   order. Extract only pure logic needed to isolate effects.
3. Run unchanged-production characterizations before mutation-sensitive negative checks; a passing
   characterization is not called RED.
4. Run the isolated scripts-only command after each tranche. Continue until the precise
   branch-inclusive figure is at least 88.00%, with no omitted module.
5. Run mutation testing selectively where practical to ensure assertions detect wrong values;
   mutation results remain advisory.

### Phase 4: close the dynamic-analysis gap (#211)

1. Append `scripts` last to `[tool.coverage.run].source`.
2. Activate the source-scope invariant test in the same change.
3. Generate a fresh combined coverage report and run patch coverage from a full-history checkout.
4. Exercise the unit coverage gate on Python 3.10, 3.11, 3.12, and 3.13. A failure on any one
   version returns the work to Phase 3; thresholds and source scope do not change.

### Phase 5: repository-wide verification and documentation

1. Run `make check-all`.
2. Because workflow files may be touched, run `make ci-local`, including its fresh uv virtual
   environment installation path.
3. Run `make docs-build` and verify this `docs/superpowers/` tree is excluded.
4. Update contributor documentation with the final precise measurements and command names.

## Success, partial failure, recovery, retry, and idempotency

### Success

The change succeeds only when widened mypy is clean, scripts-only branch coverage is at least
88.00%, both former zero modules have non-zero measured coverage, combined coverage is at least
86.00% on all four Python versions, patch coverage remains gated at 80%, and all repository
checks pass with no scripts omission.

### Partial failure

- If mypy finds additional errors after a rebase, fix those errors with local annotations or
  sound code changes and tests; do not add a broad ignore.
- If scripts-only coverage is below 88%, add behavior and branch tests to the lowest-value gaps;
  do not add exclusions or move the source change forward.
- If only one Python version falls below 86%, inspect version-dependent branches and mocks on
  that interpreter. The matrix is conjunctive, so results from other versions do not compensate.
- If combined coverage drops after `scripts` is appended, remove only that final source-list
  change while continuing Phase 3. Never lower the floor as recovery.
- If patch coverage fails, test the changed decisions or reduce unrelated refactoring. Do not
  modify the 80 target.
- If strict docs build includes these specifications, correct `exclude_docs`; do not rewrite
  contributor specifications as end-user documentation.

### Recovery and retry

All operations are repository-local and deterministic. Coverage commands delete or isolate old
data before execution, so retries cannot consume a stale `.coverage` file. Unit tests replace
external effects with deterministic fakes and use fresh temporary directories. A failed phase is
retried from its first verification command after the local defect is corrected; later phases do
not mask earlier failures.

### Idempotency

Repeated type checks and unit tests do not modify tracked files. Temporary coverage data is
created under a unique temporary directory and cleaned through a shell trap. Fixtures never
write to the repository's real test-data directory, user cache, conda environments, or external
services. Re-running the full verification sequence from a clean checkout yields the same scope
and thresholds.

## Compatibility and environment constraints

- Supported Python versions are 3.10 through 3.13. Type annotations use syntax available in
  Python 3.10.
- Ruff target version, mypy Python version, package `requires-python`, classifiers, and CI matrix
  remain aligned through the existing consistency test.
- CI creates an explicit uv virtual environment, exports `VIRTUAL_ENV`, and appends its `bin`
  directory to `GITHUB_PATH` before installing. No `uv pip install --system` is introduced.
- Pytest runs from the repository root and selects `-m unit`; every new unit test module declares
  `pytestmark = pytest.mark.unit`.
- Coverage is collected with branch mode and without parallel stale-data merging.
- No external pipeline binary is a prerequisite for the new unit tests.

## Security and permissions

- The implementation changes no GitHub token, workflow permission, secret, release, package, or
  registry access.
- Unit tests do not make network requests. Download behavior is tested through fakes and local
  temporary files; archive paths are synthetic and cannot overwrite repository or user files.
- Subprocess tests use argument vectors or exact mocked command strings and never interpolate
  untrusted sample data into a live shell.
- Tests do not log credentials, signed URLs, sample identifiers, genomic content, or environment
  secrets.
- CI requires only repository read access for these gates and normal check reporting. Patch
  coverage requires full git history but no write permission.

## Observability and diagnostic output

The implementation must make a failure self-locating:

- mypy reports file, line, error code, and final checked-file count;
- scripts-only coverage prints the exact aggregate plus every file's statements, branches, and
  misses, including explicit rows for both former zero modules;
- combined coverage prints the precise gate value rather than relying on the rounded table total;
- matrix jobs identify the Python version in the job name and artifact/report metadata;
- patch coverage reports the merge base, changed measured lines, achieved percentage, and 80%
  threshold;
- invariant tests name the missing source, changed threshold, disabled branch setting, or missing
  mypy scope directly.

No new telemetry or external reporting service is required.

## Test strategy and objective acceptance checks

The following checks are mandatory and conjunctive:

1. `mypy --no-incremental vntyper/ docker/app/ scripts/` exits 0.
2. `make type-check-all` exits 0 and demonstrably retains both mypy invocations.
3. The scripts-only branch-coverage target exits 0 at a precise result of at least 88.00%.
4. Its per-file report lists every Python module under `scripts/`; the two named zero-baseline
   modules each have executed statements and a non-zero percentage.
5. Unit tests that parse quality configuration fail when a mutation removes `scripts`, disables
   branch coverage, changes 86/86/80, or deletes the separate tests mypy pass.
6. `make test-unit-cov` exits 0 independently under Python 3.10, 3.11, 3.12, and 3.13, each at
   or above 86.00%.
7. `make patch-coverage` exits 0 at the unchanged 80 target using a real merge base.
8. `make check-all` exits 0.
9. `make ci-local` exits 0 if workflow files changed.
10. `make docs-build` exits 0 and does not publish or template any file below
    `docs/superpowers/`.

## Rollout

Land this as one reviewable PR with ordered commits that preserve a green branch:

1. invariant scaffolding, MkDocs exclusion, and #204's annotation/type-scope closure;
2. tests and focused pure-logic extraction for the two zero modules;
3. remaining scripts tests until the isolated result is at least 88%;
4. the final `coverage.source` append plus matrix evidence and contributor documentation.

The source append is deliberately last. There is no runtime feature flag or data migration: the
rollout changes development and CI enforcement only. Rollback, if required before merge, removes
the final source append while retaining all sound annotations and tests. After merge, any new
`scripts/` code is immediately subject to mypy, the combined 86% branch-inclusive repository
floor, and 80% changed-line coverage.

This specification contains no deferred design decisions. Numerical thresholds, scope, ordering,
failure behavior, and recovery policy are fixed above.

## Unresolved questions

None.
