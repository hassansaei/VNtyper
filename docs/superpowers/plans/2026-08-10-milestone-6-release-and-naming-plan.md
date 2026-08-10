# Milestone 6 Release Automation and Naming Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Release only an exact, fully tested `main` commit by promoting its existing GHCR digest, publish Python distributions through protected PyPI OIDC, make GHCR documentation truthful, and use “VNtyper 2” only for the nine approved generation-prose targets.

**Architecture:** `scripts/release_policy.py` contains pure, typed release decisions; GitHub workflow steps own Git, API, registry, artifact, wait, and summary I/O. Every `main` push produces substantive CI and Docker evidence plus a short-SHA image carrying full revision and package-version labels. Authenticated `repository_dispatch` type `vntyper_release` supplies an existing tag in `client_payload.tag`, loads `publish-pypi.yml` from the default branch, polls the tag's exact SHA checks, promotes that tested digest under a globally serialized job, then invokes an OIDC-only publisher. Manual dispatch exercises the same read/build path for an existing tag but has no route to production writes.

**Tech Stack:** Python 3.10–3.13, pytest, dataclasses and typed standard-library collections, GitHub Actions, `gh api`, Docker Buildx imagetools, GHCR, PyPI Trusted Publishing, PyPA `gh-action-pypi-publish`, actionlint, MkDocs.

## Global Constraints

- Implement the approved spec at [`../specs/2026-08-10-milestone-6-release-and-naming-design.md`](../specs/2026-08-10-milestone-6-release-and-naming-design.md); do not redesign it during execution.
- Production accepts only authenticated `repository_dispatch` type/action `vntyper_release` with a pre-existing strict `vMAJOR.MINOR.PATCH` tag in `client_payload.tag`; its peeled 40-character SHA is an ancestor of `origin/main` and its plain version equals all three version sources.
- Keep the existing `sha-<7 lowercase hexadecimal characters>` image tag; prove identity with `org.opencontainers.image.revision=<full SHA>` and `org.opencontainers.image.version=<package version>`.
- Every push to `main` runs `Lint (Ruff)`, `Type Check (mypy)`, all four `Unit Tests` matrix jobs, `Docs build (strict)`, `Build and test image`, `CI Success`, and `Docker Success`; path filters remain PR-only optimizations.
- Release polling requires all ten exact component/aggregate check names and treats a skipped component as terminal failure even when an aggregator is green.
- Poll at most 120 times with a 30-second interval. Missing/pending checks wait; terminal non-success fails immediately; exhaustion reports the last state and URLs.
- Promote only the immutable digest recorded by the successful exact-SHA Docker run's contract-v1 evidence artifact,
  after verifying the short `sha-<7>` reference (or explicitly detecting its prefix collision), full revision label,
  package-version label, and manifest. Never substitute `main` or rebuild tag-context source.
- Required aliases are exactly `vX.Y.Z`, `X.Y.Z`, `X.Y`, `X`, and `latest`. Exact aliases never move to a different digest; floating aliases never downgrade and skip unorderable existing labels with notice. `main` remains `main`.
- Serialize `promote-ghcr` across all versions with one fixed repository-wide concurrency group and
  `cancel-in-progress: false`. This is a mutual-exclusion lock, not an unbounded queue: GitHub may replace an older
  pending promotion with a newer pending run. The canceled version is explicitly rerun; it cannot have written before
  acquiring the lock, and all promotion/package operations are idempotent.
- `workflow_dispatch` accepts one existing strict version tag, is dry-run only, and can never create/move a tag, promote an alias, request OIDC, or publish a package.
- Production has no `push.tags` trigger. Every production job guards exact event name
  `repository_dispatch` and action `vntyper_release`; the workflow and controller therefore come
  from the default branch, never from a tagged historical commit.
- Keep the trusted-publisher identity exactly workflow `publish-pypi.yml`, environment `pypi`. Build without publication privilege; publish with OIDC only through the SHA-pinned action and preserve `skip-existing: true`.
- GHCR is authoritative. Do not publish to Docker Hub, restore Docker Hub credentials, delete external settings, push a tag, create a GitHub Release, send a production dispatch, publish a version, or delete `PYPI_API_TOKEN` in this implementation.
- During first-release migration, legacy token workflows stored in historical tagged commits are not
  inert merely because default-branch `push.tags` was removed. Do not create or push a tag pointing at
  a pre-milestone commit while `PYPI_API_TOKEN` exists. Delete that token only after the first successful
  production OIDC publication; only then are those legacy token workflows inert.
- Change exactly eight README generation-prose targets and one Dockerfile description from `VNtyper 2.0` to `VNtyper 2`; preserve every real version, historical reference, dependency, identifier, parsed banner, bare canonical name, and `snakemake/vntyper2*` filename.
- All new Python is fully annotated and Python 3.10 compatible. New unit files declare `pytestmark = pytest.mark.unit` and remain pure: no network, Docker, registry, or reference archive.
- Run pytest from the repository root. Never lower coverage floor 86, target 86, patch target 80, or disable branch coverage.
- Workflow jobs use explicit virtual environments, never `uv pip install --system`; shell steps use validated environment variables, quoted expansions, and `set -euo pipefail`.
- Workflow-level permissions remain `{}`; grant permissions per job only. Do not introduce `pull_request_target`.
- Because workflows change, `make ci-local` is mandatory. Because the Docker workflow changes, `make ci-local-docker` is mandatory when a Docker daemon is available and any unavailability must be reported explicitly.

---

## Requirements traceability

| Requirement | Implementation task | Objective evidence |
| --- | --- | --- |
| [#214](https://github.com/hassansaei/VNtyper/issues/214): publish working latest/semver aliases from tested image | Tasks 3, 7 | Alias-policy unit tests; digest-source and five-alias workflow tests |
| #214: Docker release guarantee equivalent to package version check | Tasks 4–7 | Every-main-push tests; exact-SHA component polling; image revision/version/digest checks |
| #214: decide Docker Hub and reconcile install docs | Task 9 | No active `saei/vntyper` command; GHCR-only install contract |
| #214: prevent future workflow/docs drift | Tasks 4, 6, 7, 9 | `tests/unit/test_release_workflow_contract.py` parses both workflows and install surfaces |
| [#218](https://github.com/hassansaei/VNtyper/issues/218): matching publisher/environment and OIDC | Task 8 | Exact filename/environment, least-privilege, no-token, pinned-action tests |
| #218: first successful release before secret deletion | Tasks 8, 10, 11 | Rollout documentation; implementation performs no settings deletion |
| [#220](https://github.com/hassansaei/VNtyper/issues/220): nine manual generation edits | Task 9 | Explicit target table assertions, not a repository-wide stale count |
| #220: grammar repair and version/identifier preservation | Task 9 | Exact grammar/protected-name/version assertions and existing consistency tests |
| Exact failure, rerun, security, and diagnostics contracts | Tasks 2–8, 10 | Unit decisions, workflow job graph/permissions, recovery messages, summary literals |

## File structure and ownership

| File | Action | Responsibility |
| --- | --- | --- |
| `scripts/release_policy.py` | Create | Pure typed version, check-run, and alias decisions; JSON-compatible values only; no GitHub/registry I/O |
| `scripts/release_evidence.py` | Create | Typed local-JSON adapter selecting exact-SHA Docker runs and validating contract-v1 evidence/image labels; no network or registry I/O |
| `tests/unit/test_release_policy.py` | Create | Exhaustive decision tests, including invalid inputs, terminal states, reruns, conflicts, and dry-run execution flags |
| `tests/unit/test_release_evidence.py` | Create | Exact-run selection, malformed/stale evidence, digest/label/collision, and CLI output tests |
| `tests/unit/test_release_workflow_contract.py` | Create | Parsed and executable workflow graph/permission/trigger/label checks plus GHCR documentation and naming invariants |
| `.github/workflows/ci-tests.yml` | Modify | Make path-derived `python` and `docs` outputs true for every main push while retaining PR filtering |
| `.github/workflows/docker-build.yml` | Modify | Make image output true for every main push; retain short SHA; add full revision/package labels; upload exact-SHA digest evidence; remove dead tag metadata |
| `.github/workflows/publish-pypi.yml` | Rewrite in place | Default-branch repository-dispatch boundary, exact-tag validation, bounded evidence polling, missing-image recovery, package artifact build, serialized digest promotion, dry-run, OIDC publish, summary |
| `README.md` | Modify | GHCR-only commands, stable/rolling tag explanation, eight generation-name edits, grammar repair |
| `docs/getting-started/installation.md` | Modify | Remove active Docker Hub command; document supported GHCR pull |
| `docs/user-guide/docker.md` | Modify | Remove active Docker Hub/Apptainer commands; document aliases and `main` semantics |
| `docker/Dockerfile` | Modify | One approved `VNtyper 2` OCI description target |
| `docs/development/ci-followups.md` | Modify | Mark B4 implemented, preserve first-live-release/token-deletion follow-up |
| `AGENTS.md` | Modify | Correct release, all-main-gates, OIDC, short-SHA/full-label, retry, and scripts type-check guidance |

## Frozen policy interfaces

Implement these names and signatures exactly so workflow and other milestone plans can consume one contract:

```python
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

REQUIRED_CHECK_NAMES: tuple[str, ...] = (
    "Lint (Ruff)",
    "Type Check (mypy)",
    "Unit Tests (Python 3.10)",
    "Unit Tests (Python 3.11)",
    "Unit Tests (Python 3.12)",
    "Unit Tests (Python 3.13)",
    "Docs build (strict)",
    "CI Success",
    "Build and test image",
    "Docker Success",
)


@dataclass(frozen=True)
class ReleaseVersion:
    tag: str
    plain: str
    major: int
    minor: int
    patch: int


@dataclass(frozen=True)
class CheckVerdict:
    name: str
    state: Literal["missing", "pending", "success", "failure"]
    conclusion: str | None
    check_run_id: int | None
    url: str | None
    reason: str


@dataclass(frozen=True)
class CheckPoll:
    action: Literal["wait", "success", "fail", "timeout"]
    verdicts: tuple[CheckVerdict, ...]
    attempt: int
    elapsed_seconds: int


@dataclass(frozen=True)
class AliasState:
    digest: str
    version: str | None


@dataclass(frozen=True)
class AliasUpdate:
    alias: str
    decision: Literal["create", "advance", "no-op", "skip-newer", "skip-unorderable", "fail-conflict"]
    execute: bool
    reason: str


```

- `parse_release_tag(tag: str) -> ReleaseVersion`
- `required_aliases(version: ReleaseVersion) -> tuple[str, ...]`
- `classify_check_runs(candidate_sha: str, check_runs: Sequence[Mapping[str, object]], *, attempt: int, max_attempts: int = 120) -> CheckPoll`
- `plan_alias_updates(version: ReleaseVersion, source_digest: str, current_aliases: Mapping[str, AliasState | None], *, dry_run: bool) -> tuple[AliasUpdate, ...]`

`classify_check_runs` is pure: the caller supplies one exact-SHA API snapshot and the one-based
attempt number. It filters `head_sha == candidate_sha` and `app.slug == "github-actions"`, selects
the highest check-run `id` for each exact name, and makes skipped/neutral terminal conclusions
failures. It records `elapsed_seconds = (attempt - 1) * 30`. `plan_alias_updates` expects all five aliases as mapping keys; absence is represented by
`None`. `dry_run=True` preserves each decision but forces every `execute` flag false.

---

### Task 1: Strict release versions and alias derivation

**Files:**
- Create: `scripts/release_policy.py`
- Create: `tests/unit/test_release_policy.py`

**Interfaces:**
- Produces: `ReleaseVersion`, `parse_release_tag(tag: str)`, and `required_aliases(version: ReleaseVersion)` exactly as frozen above.
- Consumes: Python standard library only.

- [ ] **Step 1: Write RED tests for valid tags and exact aliases (2–5 min)**

First create an importable `scripts/release_policy.py` scaffold containing the frozen dataclasses/signatures with
function bodies that raise `NotImplementedError("not implemented")`; this supplies no behavior and prevents collection
failure. Then create the marked test file:

```python
import pytest

from scripts.release_policy import parse_release_tag, required_aliases

pytestmark = pytest.mark.unit


def test_release_tag_produces_all_five_aliases_in_promotion_order() -> None:
    version = parse_release_tag("v2.10.3")
    assert (version.plain, version.major, version.minor, version.patch) == ("2.10.3", 2, 10, 3)
    assert required_aliases(version) == ("v2.10.3", "2.10.3", "2.10", "2", "latest")
```

- [ ] **Step 2: Add RED rejection cases with precise unsafe literals (2–5 min)**

```python
@pytest.mark.parametrize(
    "tag",
    ("2.0.10", "v2.0", "v2.0.10-rc1", "v2.0.10+local", "v02.0.10", "v2.00.10", "v2.0.010", "v2.0.10;echo pwned"),
)
def test_release_tag_rejects_every_non_strict_form(tag: str) -> None:
    with pytest.raises(ValueError, match="strict vMAJOR.MINOR.PATCH"):
        parse_release_tag(tag)
```

- [ ] **Step 3: Run RED and confirm behavioral failure (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_policy.py -q`
Expected: test execution fails with `NotImplementedError: not implemented`, not collection/import failure.

- [ ] **Step 4: Implement the minimal strict parser and alias tuple (2–5 min)**

Use `re.fullmatch(r"v(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)", tag)`, populate the frozen dataclass, and raise:

```python
msg = f"Release tag {tag!r} must be strict vMAJOR.MINOR.PATCH with no prerelease or build suffix."
raise ValueError(msg)
```

- [ ] **Step 5: Run GREEN, formatting, and type checking (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_policy.py -q && ruff check scripts/release_policy.py tests/unit/test_release_policy.py && mypy scripts/release_policy.py`
Expected: all tests pass; Ruff and mypy report no errors.

- [ ] **Step 6: Commit the independently usable parser (2–5 min)**

```bash
git add scripts/release_policy.py tests/unit/test_release_policy.py
git commit -m "feat(release): define strict version aliases"
```

### Task 2: Exact-SHA component and aggregator check policy

**Files:**
- Modify: `scripts/release_policy.py`
- Modify: `tests/unit/test_release_policy.py`

**Interfaces:**
- Consumes: raw Check Runs API mappings with `name`, `head_sha`, `id`, `status`, `conclusion`, `html_url`, and nested `app.slug`.
- Produces: `REQUIRED_CHECK_NAMES`, `CheckVerdict`, `CheckPoll`, and `classify_check_runs(...)` exactly as frozen above.

- [ ] **Step 1: Add a local check-run factory and RED exactness test (2–5 min)**

Before importing the new names in tests, add `CheckVerdict`, `CheckPoll`, `REQUIRED_CHECK_NAMES`, and a
`classify_check_runs` scaffold that raises `NotImplementedError("not implemented")` to the already importable module.

```python
def _check(name: str, *, sha: str = "a" * 40, run_id: int = 1, status: str = "completed", conclusion: str | None = "success", app: str = "github-actions") -> dict[str, object]:
    return {"name": name, "head_sha": sha, "id": run_id, "status": status, "conclusion": conclusion, "html_url": f"https://example.test/check/{run_id}", "app": {"slug": app}}


def test_all_ten_latest_exact_sha_github_actions_checks_are_required() -> None:
    runs = [_check(name) for name in REQUIRED_CHECK_NAMES]
    poll = classify_check_runs("a" * 40, runs, attempt=1)
    assert poll.action == "success"
    assert tuple(verdict.name for verdict in poll.verdicts) == REQUIRED_CHECK_NAMES
    assert (poll.attempt, poll.elapsed_seconds) == (1, 0)
```

- [ ] **Step 2: Add RED newest-attempt, skipped, missing, and timeout tests (2–5 min)**

Assert these precise outcomes:

```python
assert classify_check_runs(SHA, successful + [_check("Docker Success", run_id=99, conclusion="failure")], attempt=1).action == "fail"
assert classify_check_runs(SHA, successful_without_docs, attempt=119).action == "wait"
assert classify_check_runs(SHA, successful_without_docs, attempt=120).action == "timeout"
assert classify_check_runs(SHA, with_skipped_build_component, attempt=1).action == "fail"
```

Also inject successful same-name records with another `head_sha` and `app.slug="external-ci"`; assert neither can satisfy a missing required check.
Assert the timeout result records `(attempt, elapsed_seconds) == (120, 3570)`.

- [ ] **Step 3: Run RED and confirm the behavioral failure (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_policy.py -q`
Expected: test execution fails with `NotImplementedError: not implemented`.

- [ ] **Step 4: Implement latest-record selection and terminal-state classification (2–5 min)**

Implement the frozen function. Treat only `status == "completed" and conclusion == "success"` as success; `None`, queued/requested/waiting/pending/in-progress statuses are pending; every completed non-success conclusion, including `skipped` and `neutral`, is failure. Return `fail` before timeout, `success` only when all ten succeed, otherwise `timeout` at `attempt >= max_attempts`, else `wait`.

- [ ] **Step 5: Refactor messages into deterministic values and run GREEN (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_policy.py -q && ruff check scripts/release_policy.py tests/unit/test_release_policy.py && mypy scripts/release_policy.py`
Expected: exact-SHA, newest-ID, skipped, failure, wait, and timeout cases all pass.

- [ ] **Step 6: Commit check policy (2–5 min)**

```bash
git add scripts/release_policy.py tests/unit/test_release_policy.py
git commit -m "feat(release): require exact-sha workflow evidence"
```

### Task 3: Immutable and monotonic alias planning

**Files:**
- Modify: `scripts/release_policy.py`
- Modify: `tests/unit/test_release_policy.py`

**Interfaces:**
- Consumes: parsed `ReleaseVersion`, one verified `sha256:` source digest, and five `AliasState | None` observations.
- Produces: `AliasState`, `AliasUpdate`, and `plan_alias_updates(...)` exactly as frozen above.

- [ ] **Step 1: Write RED exact-alias immutability tests (2–5 min)**

Add importable `AliasState`, `AliasUpdate`, and `plan_alias_updates` scaffold definitions before adding the tests; the
planner raises `NotImplementedError("not implemented")` until Step 4.

```python
def test_exact_aliases_create_noop_or_fail_but_never_advance() -> None:
    version = parse_release_tag("v2.0.10")
    current = {alias: None for alias in required_aliases(version)}
    source = "sha256:" + "a" * 64
    current["v2.0.10"] = AliasState(source, "2.0.10")
    current["2.0.10"] = AliasState("sha256:" + "b" * 64, "2.0.10")
    updates = plan_alias_updates(version, source, current, dry_run=False)
    assert [(u.alias, u.decision, u.execute) for u in updates[:2]] == [
        ("v2.0.10", "no-op", False),
        ("2.0.10", "fail-conflict", False),
    ]
```

- [ ] **Step 2: Write RED floating anti-downgrade and dry-run tests (2–5 min)**

Cover `2.0` absent → `create`, `2` at `2.0.9` → `advance`, `latest` at `2.1.0` → `skip-newer`, and an existing alias with `version=None` → `skip-unorderable`. Add:

```python
def test_equal_floating_version_with_a_different_digest_is_a_hard_conflict() -> None:
    version = parse_release_tag("v2.0.10")
    source = "sha256:" + "a" * 64
    current = {alias: None for alias in required_aliases(version)}
    current["latest"] = AliasState("sha256:" + "b" * 64, "2.0.10")
    latest = plan_alias_updates(version, source, current, dry_run=False)[-1]
    assert (latest.decision, latest.execute) == ("fail-conflict", False)
```

Re-run the planner after applying each prefix of executable operations and assert all later plans contain only
`no-op`/skip decisions for completed aliases. With `dry_run=True`, assert decisions are identical but every
`execute is False`.

- [ ] **Step 3: Run RED and confirm behavioral planner failure (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_policy.py -q`
Expected: test execution fails with `NotImplementedError: not implemented`.

- [ ] **Step 4: Implement the minimal five-alias planner (2–5 min)**

Validate `source_digest` with `r"sha256:[0-9a-f]{64}"`, require exactly the `required_aliases()` key set, compare parsed numeric triples rather than strings, plan exact aliases first, and set `execute = decision in {"create", "advance"} and not dry_run`.

- [ ] **Step 5: Run mutation-sensitive GREEN assertions (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_policy.py -q --cov=scripts.release_policy --cov-branch --cov-report=term-missing`
Expected: unit tests pass with every decision branch in strict parsing, exact conflicts, downgrade comparison, and `execute` exercised; no listed missed branch belongs to the new helper.

- [ ] **Step 6: Commit alias policy (2–5 min)**

```bash
git add scripts/release_policy.py tests/unit/test_release_policy.py
git commit -m "feat(release): plan idempotent GHCR aliases"
```

### Task 4: Substantive CI and Docker evidence on every main push

**Files:**
- Create: `tests/unit/test_release_workflow_contract.py`
- Modify: `.github/workflows/ci-tests.yml`
- Modify: `.github/workflows/docker-build.yml`

**Interfaces:**
- Produces: exact check-run names in `REQUIRED_CHECK_NAMES`; short-SHA source tag; full revision/package-version OCI labels.
- Produces: run-attempt-scoped `docker-release-evidence-<full SHA>-<run attempt>` JSON with contract version, SHA, digest, run ID/attempt,
  revision, and version.
- Consumes: existing PR path filters, base-image hash, Docker test tiers, and aggregate job names unchanged.

- [ ] **Step 1: Write RED workflow parsing helpers and main-push assertions (2–5 min)**

```python
from pathlib import Path
from typing import Any

import pytest
import yaml

pytestmark = pytest.mark.unit
ROOT = Path(__file__).resolve().parents[2]


def _workflow(name: str) -> dict[str, Any]:
    return yaml.load((ROOT / ".github" / "workflows" / name).read_text(), Loader=yaml.BaseLoader)


def test_main_push_forces_every_substantive_ci_and_docker_scope_true() -> None:
    ci = _workflow("ci-tests.yml")
    docker = _workflow("docker-build.yml")
    assert ci["concurrency"]["group"] == "ci-${{ github.event_name == 'push' && github.sha || github.ref }}"
    assert docker["concurrency"]["group"] == "docker-${{ github.event_name == 'push' && github.sha || github.ref }}"
    assert ci["concurrency"]["cancel-in-progress"] == "${{ github.event_name == 'pull_request' }}"
    assert docker["concurrency"]["cancel-in-progress"] == "false"
    for output in ("python", "docs"):
        assert "github.event_name == 'push'" in ci["jobs"]["changes"]["outputs"][output]
    assert "github.event_name == 'push'" in docker["jobs"]["changes"]["outputs"]["image"]
    assert {ci["jobs"][job]["if"] for job in ("lint", "typecheck", "test-unit")} == {
        "needs.changes.outputs.python == 'true'"
    }
    assert ci["jobs"]["docs"]["if"] == "needs.changes.outputs.docs == 'true'"
    assert set(ci["jobs"]["ci-success"]["needs"]) == {"changes", "lint", "typecheck", "test-unit", "docs"}
    assert docker["jobs"]["base-status"]["if"] == "needs.changes.outputs.image == 'true'"
    assert "needs.changes.outputs.image == 'true'" in docker["jobs"]["build-and-test"]["if"]
    assert set(docker["jobs"]["docker-success"]["needs"]) == {
        "changes", "base-status", "build-base", "build-and-test"
    }
```

Add assertions that PR filters remain present and that static job names plus expansion of
`Unit Tests (Python ${{ matrix.python-version }})` across `3.10`, `3.11`, `3.12`, and `3.13`
equal `REQUIRED_CHECK_NAMES` without omissions. Assert neither workflow uses `pull_request_target`; fork PRs retain
read-only permissions and the existing explicit missing-base diagnostic. The per-SHA main concurrency groups prevent a
third main push from replacing a different SHA's pending evidence; PR runs remain ref-grouped.

- [ ] **Step 2: Write RED image identity assertions (2–5 min)**

Assert Docker metadata still includes `type=sha`, excludes `type=ref,event=tag`, and the build labels contain literal
keys `org.opencontainers.image.revision` and `org.opencontainers.image.version` but do not pass
`${{ steps.meta.outputs.labels }}` (which overwrites the Dockerfile description). Parse the `upload-artifact` step and
assert action `actions/upload-artifact@v5`, name
`docker-release-evidence-${{ github.sha }}-${{ github.run_attempt }}`, `if-no-files-found: error`, and
`retention-days: "90"`. Assert a
JSON-producing prior step whose payload uses the registry-reported `{{.Manifest.Digest}}` inspected after the short tag
is pushed, includes
`${{ github.run_id }}`, `${{ github.run_attempt }}`, and
`"contract_version": 1`.

- [ ] **Step 3: Run RED and observe path-filter/label failures (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_workflow_contract.py -q`
Expected: failures show main-push output overrides and explicit OCI identity labels are absent; dead tag metadata remains.

- [ ] **Step 4: Make CI path filters PR-only (2–5 min)**

Change `changes.outputs.python` and `.docs` to expressions equivalent to:

```yaml
python: ${{ github.event_name == 'push' || steps.filter.outputs.python == 'true' }}
docs: ${{ github.event_name == 'push' || steps.filter.outputs.docs == 'true' }}
```

Keep component `if:` expressions reading those normalized outputs. Update comments: skipped jobs are legitimate only on PRs; every main push executes lint, typecheck, four unit jobs, docs, and aggregate.
Replace the workflow concurrency block explicitly:

```yaml
concurrency:
  group: ci-${{ github.event_name == 'push' && github.sha || github.ref }}
  cancel-in-progress: ${{ github.event_name == 'pull_request' }}
```

- [ ] **Step 5: Make Docker path filters PR-only and attach image identity (2–5 min)**

Normalize `changes.outputs.image` the same way. Before the build, execute `vntyper/version.py` via `runpy` into step output `package`; produce RFC3339 `created`, full `${GITHUB_SHA}` `revision`, and package `version`; pass only those three label keys to `docker/build-push-action`. Retain metadata `type=ref,event=branch`, `type=ref,event=pr`, `type=sha`; remove `type=ref,event=tag`.

Replace the Docker workflow concurrency block explicitly:

```yaml
concurrency:
  group: docker-${{ github.event_name == 'push' && github.sha || github.ref }}
  # Preserve the existing fork/base-build starvation contract: even superseded PR
  # runs may be the only run able to finish and publish a newly required base.
  cancel-in-progress: false
```

Retain the live workflow's complete evaluation-timing/base-build-starvation rationale comment above this block; only
the group expression changes. Do not replace that load-bearing explanation with the abbreviated snippet comment.

After the build, write this exact schema to `docker-release-evidence.json` with `python -c`/`json.dump` using values
passed through environment variables, then upload it in the same run:

```json
{"contract_version": 1, "sha": "<40 hex>", "digest": "sha256:<64 hex>", "run_id": 1, "run_attempt": 1, "revision": "<40 hex>", "version": "X.Y.Z"}
```

Use these literal steps immediately after the existing explicit `docker push` loop; validate the pushed tag's
registry-reported manifest digest and OCI labels before upload. Do not interpolate
JSON with shell string concatenation:

```yaml
      - id: registry-digest
        if: ${{ github.event_name == 'push' && github.ref == 'refs/heads/main' }}
        env:
          IMAGE: ghcr.io/hassansaei/vntyper
          VERSION: ${{ steps.package.outputs.version }}
        run: |
          set -euo pipefail
          SHORT_SHA=${GITHUB_SHA:0:7}
          for ATTEMPT in 1 2 3 4 5 6; do
            if DIGEST=$(docker buildx imagetools inspect "$IMAGE:sha-$SHORT_SHA" \
                --format '{{.Manifest.Digest}}'); then break; fi
            if [ "$ATTEMPT" -eq 6 ]; then echo "pushed short-SHA manifest was not observable" >&2; exit 1; fi
            sleep 5
          done
          case "$DIGEST" in sha256:????????????????????????????????????????????????????????????????) ;; *) exit 1 ;; esac
          docker buildx imagetools inspect "$IMAGE:sha-$SHORT_SHA" --format '{{json .Image}}' > pushed-image.json
          python - "$GITHUB_SHA" "$VERSION" pushed-image.json <<'PY'
          import json, sys
          image = json.load(open(sys.argv[3], encoding="utf-8"))
          labels = image.get("config", {}).get("Labels", {}) or {}
          if labels.get("org.opencontainers.image.revision") != sys.argv[1]:
              raise SystemExit("pushed image revision label mismatch")
          if labels.get("org.opencontainers.image.version") != sys.argv[2]:
              raise SystemExit("pushed image version label mismatch")
          PY
          printf 'value=%s\n' "$DIGEST" >> "$GITHUB_OUTPUT"
      - id: release-evidence
        if: ${{ github.event_name == 'push' && github.ref == 'refs/heads/main' }}
        env:
          DIGEST: ${{ steps.registry-digest.outputs.value }}
          REVISION: ${{ github.sha }}
          RUN_ATTEMPT: ${{ github.run_attempt }}
          RUN_ID: ${{ github.run_id }}
          VERSION: ${{ steps.package.outputs.version }}
        run: |
          set -euo pipefail
          python - <<'PY'
          import json
          import os
          import re

          sha = os.environ["REVISION"]
          digest = os.environ["DIGEST"]
          if re.fullmatch(r"[0-9a-f]{40}", sha) is None:
              raise SystemExit("invalid full revision")
          if re.fullmatch(r"sha256:[0-9a-f]{64}", digest) is None:
              raise SystemExit("invalid pushed image digest")
          payload = {"contract_version": 1, "sha": sha, "digest": digest,
                     "run_id": int(os.environ["RUN_ID"]), "run_attempt": int(os.environ["RUN_ATTEMPT"]),
                     "revision": sha, "version": os.environ["VERSION"]}
          with open("docker-release-evidence.json", "w", encoding="utf-8") as handle:
              json.dump(payload, handle, separators=(",", ":"))
          PY
      - uses: actions/upload-artifact@v5
        if: ${{ github.event_name == 'push' && github.ref == 'refs/heads/main' }}
        with:
          name: docker-release-evidence-${{ github.sha }}-${{ github.run_attempt }}
          path: docker-release-evidence.json
          if-no-files-found: error
          retention-days: "90"
```

- [ ] **Step 6: Run GREEN and the repository workflow lint target (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_workflow_contract.py -q && make lint-actions`
Expected: main/PR scope, ten names, short SHA, full labels, description ownership, and YAML semantics pass.

- [ ] **Step 7: Commit main evidence contract (2–5 min)**

```bash
git add .github/workflows/ci-tests.yml .github/workflows/docker-build.yml tests/unit/test_release_workflow_contract.py
git commit -m "ci(release): build substantive evidence on main"
```

### Task 5: Default-branch existing-tag validation and dry-run boundary

**Files:**
- Modify: `.github/workflows/publish-pypi.yml`
- Modify: `tests/unit/test_release_workflow_contract.py`

**Interfaces:**
- Consumes: `parse_release_tag`, fetched `origin/main`, `repository_dispatch.client_payload.tag`, or manual dispatch's existing tag.
- Produces from `validate-release`: `mode`, `tag`, `version`, `sha`, and `short_sha` job outputs.

- [ ] **Step 1: Write RED trigger and production-guard tests (2–5 min)**

Assert the trigger set is exactly `repository_dispatch` and `workflow_dispatch`;
`repository_dispatch.types == ["vntyper_release"]`; no `push`/`push.tags` trigger exists;
`workflow_dispatch.inputs.tag.required == "true"`; no manual input is named `publish`, `dry_run`, `sha`, or
`version`; and no step contains `git tag`, `git push`, or `gh release create`. Production-job guards are added and
tested in the task that creates each job; Task 5 does not create placeholders.

- [ ] **Step 2: Write RED validation contract tests (2–5 min)**

Assert `validate-release` uses `actions/checkout@v7` with `fetch-depth: "0"` and `persist-credentials: "false"`; runs `git rev-parse "refs/tags/${TAG}^{commit}"`, `git merge-base --is-ancestor "$SHA" origin/main`, `parse_release_tag`, and the precise version consistency test `pytest -m unit tests/unit/test_version_consistency.py -q` before any downstream job.
Assert the first checkout path is `controller`, the second is `candidate` at `${{ steps.resolve.outputs.sha }}`;
`parse_release_tag` runs in controller while `runpy` and pytest run with `working-directory: candidate`. This prevents
manual old-tag validation from importing helpers out of the candidate or testing current-main metadata.
Assert the candidate step creates `.release-venv` and installs exactly the test's collection/runtime requirements
`pytest packaging PyYAML requests`; omitting `requests` is a collection failure because repository `tests/conftest.py`
imports it before the selected test runs. Execute a mismatch fixture and assert `summary_json` preserves the expected
tag version plus observed package, `CITATION.cff`, and changelog versions and a Boolean verdict for each. Execute the
entire candidate shell step with fake successful `pip`/`pytest` commands for four hostile fixtures: package-only false,
citation-only false, changelog-only false, and three internally consistent sources that all mismatch the tag. Every
fixture must exit nonzero while retaining `version_test_passed: true`. Execute `resolve` with a former `push` event, a
wrong repository-dispatch action, and shell-metacharacter tags; each must fail before the fake Git executable is called.

- [ ] **Step 3: Run RED against the token-based combined workflow (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_workflow_contract.py -q`
Expected: failures name the tag-push trigger, absent repository dispatch, permissive event fallback, missing
`validate-release`, full-history checkout, main-ancestry proof, and false metadata verdicts exiting zero.

- [ ] **Step 4: Replace the header and create `validate-release` minimally (2–5 min)**

Document authenticated production repository dispatch over an existing `vX.Y.Z` tag and `workflow_dispatch`
existing-tag dry run. Build the workflow header/job from this
literal skeleton, filling only the validated shell body described below:

```yaml
name: Publish PyPI and promote GHCR
on:
  repository_dispatch:
    types: [vntyper_release]
  workflow_dispatch:
    inputs:
      tag:
        description: Existing strict vMAJOR.MINOR.PATCH tag to inspect without writes
        required: true
        type: string
permissions: {}
concurrency:
  group: release-${{ inputs.tag || github.event.client_payload.tag }}
  cancel-in-progress: false
jobs:
  validate-release:
    runs-on: ubuntu-24.04
    permissions:
      contents: read
    outputs:
      mode: ${{ steps.resolve.outputs.mode }}
      tag: ${{ steps.resolve.outputs.tag }}
      version: ${{ steps.resolve.outputs.version }}
      sha: ${{ steps.resolve.outputs.sha }}
      short_sha: ${{ steps.resolve.outputs.short_sha }}
      summary_json: ${{ steps.validate-result.outputs.summary_json }}
    steps:
      - uses: actions/checkout@v7
        with:
          fetch-depth: "0"
          persist-credentials: "false"
          path: controller
      - id: resolve
        working-directory: controller
        env:
          DISPATCH_TAG: ${{ inputs.tag }}
          EVENT_ACTION: ${{ github.event.action }}
          EVENT_NAME: ${{ github.event_name }}
          PRODUCTION_TAG: ${{ github.event.client_payload.tag }}
        run: |
          set -euo pipefail
          if [ "$EVENT_NAME" = "workflow_dispatch" ]; then
            MODE=dry-run
            TAG=$DISPATCH_TAG
          elif [ "$EVENT_NAME" = "repository_dispatch" ] && [ "$EVENT_ACTION" = "vntyper_release" ]; then
            MODE=production
            TAG=$PRODUCTION_TAG
          else
            echo "::error::Unsupported release event: ${EVENT_NAME}/${EVENT_ACTION}"
            exit 1
          fi
          VERSION=$(python - "$TAG" <<'PY'
          import sys
          from scripts.release_policy import parse_release_tag
          print(parse_release_tag(sys.argv[1]).plain)
          PY
          )
          SHA=$(git rev-parse "refs/tags/${TAG}^{commit}")
          git fetch --no-tags origin main
          git merge-base --is-ancestor "$SHA" origin/main
          SHORT_SHA=${SHA:0:7}
          for pair in "mode=$MODE" "tag=$TAG" "version=$VERSION" "sha=$SHA" "short_sha=$SHORT_SHA"; do
            printf '%s\n' "$pair" >> "$GITHUB_OUTPUT"
          done
      - uses: actions/checkout@v7
        with:
          ref: ${{ steps.resolve.outputs.sha }}
          persist-credentials: "false"
          path: candidate
      - id: candidate
        working-directory: candidate
        env:
          MODE: ${{ steps.resolve.outputs.mode }}
          SHA: ${{ steps.resolve.outputs.sha }}
          TAG: ${{ steps.resolve.outputs.tag }}
          VERSION: ${{ steps.resolve.outputs.version }}
        run: |
          set -euo pipefail
          python - <<'PY'
          import json
          import os
          import re
          import runpy
          from pathlib import Path

          expected = os.environ["VERSION"]
          package = runpy.run_path("vntyper/version.py")["__version__"]
          citation_match = re.search(r'^version:\s*"([^"]+)"', Path("CITATION.cff").read_text(encoding="utf-8"), re.MULTILINE)
          changelog_match = re.search(r"^##\s+([0-9]+\.[0-9]+\.[0-9]+)(?:\s|$)", Path("docs/about/changelog.md").read_text(encoding="utf-8"), re.MULTILINE)
          if citation_match is None or changelog_match is None:
              raise SystemExit("candidate release metadata format is unreadable")
          observed = {
              "expected_version": expected,
              "package": {"observed": package, "matches": package == expected},
              "citation": {"observed": citation_match.group(1), "matches": citation_match.group(1) == expected},
              "changelog": {"observed": changelog_match.group(1), "matches": changelog_match.group(1) == expected},
              "version_test_exit_code": None,
              "version_test_passed": False,
          }
          Path("candidate-version-observations.json").write_text(
              json.dumps(observed, separators=(",", ":")), encoding="utf-8"
          )
          PY
          python3 -m venv .release-venv
          .release-venv/bin/pip install pytest packaging PyYAML requests
          set +e
          .release-venv/bin/pytest -m unit tests/unit/test_version_consistency.py -q
          TEST_STATUS=$?
          set -e
          python - "$TEST_STATUS" <<'PY'
          import json
          import sys
          from pathlib import Path

          path = Path("candidate-version-observations.json")
          observed = json.loads(path.read_text(encoding="utf-8"))
          test_status = int(sys.argv[1])
          observed["version_test_exit_code"] = test_status
          observed["version_test_passed"] = test_status == 0
          path.write_text(json.dumps(observed, separators=(",", ":")), encoding="utf-8")
          if not all(observed[source]["matches"] for source in ("package", "citation", "changelog")):
              raise SystemExit("candidate release metadata does not match the requested tag version")
          raise SystemExit(test_status)
          PY
      - id: validate-result
        if: ${{ always() }}
        env:
          CANDIDATE_OUTCOME: ${{ steps.candidate.outcome }}
          MODE: ${{ steps.resolve.outputs.mode }}
          RESOLVE_OUTCOME: ${{ steps.resolve.outcome }}
          SHA: ${{ steps.resolve.outputs.sha }}
          TAG: ${{ steps.resolve.outputs.tag }}
          VERSION: ${{ steps.resolve.outputs.version }}
        run: |
          python - <<'PY' >> "$GITHUB_OUTPUT"
          import json, os, pathlib
          observations_path = pathlib.Path("candidate/candidate-version-observations.json")
          observations = (
              json.loads(observations_path.read_text(encoding="utf-8"))
              if observations_path.is_file()
              else {"expected_version": os.environ["VERSION"] or None,
                    "package": {"observed": None, "matches": False},
                    "citation": {"observed": None, "matches": False},
                    "changelog": {"observed": None, "matches": False},
                    "version_test_exit_code": None, "version_test_passed": False}
          )
          payload = {"mode": os.environ["MODE"] or "unresolved", "tag": os.environ["TAG"] or None,
                     "version": os.environ["VERSION"] or None, "sha": os.environ["SHA"] or None,
                     "main_ancestor": os.environ["RESOLVE_OUTCOME"] == "success",
                     "resolve_outcome": os.environ["RESOLVE_OUTCOME"],
                     "candidate_outcome": os.environ["CANDIDATE_OUTCOME"],
                     "version_validation": observations}
          print("summary_json=" + json.dumps(payload, separators=(",", ":")))
          PY
```

The controller checkout always holds the current default-branch coordinator helpers/workflow. Resolve either dispatch
input only through `git rev-parse "refs/tags/${TAG}^{commit}"` there; require the tag already exists. Accept production
only when `github.event_name == 'repository_dispatch'` and `github.event.action == 'vntyper_release'`; `github.sha`
identifies the default-branch coordinator and is never treated as candidate identity. Validate the payload tag before
using it in shell commands, fetch `origin/main`, prove ancestry, then check out the resolved candidate separately. Read
the candidate's `vntyper/version.py` via `runpy`,
create the explicit Python environment inside `candidate/`, install only `pytest`, `packaging`, `PyYAML`, and
`requests`, run the candidate's version-consistency test, and export all six outputs. Before environment setup, persist
the observed package/citation/changelog values and per-source match verdicts; after pytest, add its exit code/verdict.
Exit nonzero if pytest fails or if any per-source match verdict is false; internal agreement among all three version
sources is not sufficient when they all mismatch the requested tag.
The always-run serializer retains those observations on a mismatch or uses explicit unavailable values if setup failed.
`summary_json` is single-line JSON written through `$GITHUB_OUTPUT`.

- [ ] **Step 5: Finish the one-job graph without placeholders (2–5 min)**

Keep Task 5's workflow valid with only `validate-release`; later tasks append fully specified jobs. Retain the explicit
workflow-level concurrency block from the skeleton exactly. Do not grant
`id-token: write` or `packages: write` yet. Assert the parsed job set is exactly `{"validate-release"}` at this commit.

- [ ] **Step 6: Run GREEN and test hostile manual tags (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_workflow_contract.py -q && make lint-actions`
Expected: the default-branch repository trigger, existing-tag-only inputs, ancestry/version proof, explicit false-verdict
failure, no tag creation, exact event/action boundary, and hostile event/tag probes pass.

- [ ] **Step 7: Commit validation boundary (2–5 min)**

```bash
git add .github/workflows/publish-pypi.yml tests/unit/test_release_workflow_contract.py
git commit -m "ci(release): validate existing tagged main commits"
```

### Task 6: Bounded exact-SHA gates and missing-image recovery

**Files:**
- Modify: `.github/workflows/publish-pypi.yml`
- Modify: `tests/unit/test_release_workflow_contract.py`
- Create: `scripts/release_evidence.py`
- Create: `tests/unit/test_release_evidence.py`

**Interfaces:**
- Consumes: validation outputs, `classify_check_runs`, GitHub Check Runs API, GHCR short-SHA manifest.
- Produces from `wait-for-release-gates`: `source_ref`, `source_digest`, `source_revision`, `source_version`, and a JSON check summary.
- Produces: `eligible_successful_runs(payload: Mapping[str, object] | Sequence[Mapping[str, object]], sha: str) -> tuple[Mapping[str, object], ...]`, ordered by descending numeric run ID after filtering exact SHA, completed/success, `push`, and `main`,
  `validate_evidence(payload: Mapping[str, object], *, sha: str, version: str, run_id: int, run_attempt: int) -> Mapping[str, object]`, and
  `validate_image(evidence: Mapping[str, object], image: Mapping[str, object], short_image: Mapping[str, object], *, immutable_digest: str, short_digest: str) -> Mapping[str, object]`.
- Produces: `summarize_evidence_state(...) -> Mapping[str, object]`; missing local phase files become explicit
  unavailable/failed fields while every present run URL/digest/revision/version is preserved.

- [ ] **Step 1: Write RED permissions, names, and finite-poll assertions (2–5 min)**

Assert `wait-for-release-gates.permissions` is exactly `contents: read`, `checks: read`, `actions: read`,
`packages: read`; `timeout-minutes == "70"`; workflow literals include all `REQUIRED_CHECK_NAMES`,
`RELEASE_POLL_ATTEMPTS: "120"`, `RELEASE_POLL_INTERVAL_SECONDS: "30"`, `RELEASE_API_ATTEMPTS: "3"`, and
`RELEASE_API_RETRY_SECONDS: "5"`; the polling step calls `classify_check_runs` with
`max_attempts=int(...)`, has `env.GH_TOKEN == "${{ secrets.GITHUB_TOKEN }}"`, retries `gh api` at most three times,
has no unbounded `while true`, and exits nonzero after the outer loop instead of falling through.
Assert the evidence-download step independently carries the same `GH_TOKEN`; job permissions alone are not CLI
authentication.
Assert a `docker/login-action@v4` step using `GITHUB_TOKEN` precedes every production GHCR inspection and is guarded by
exactly `github.event_name == 'repository_dispatch' && github.event.action == 'vntyper_release'`; manual dispatch must
still perform no registry login or write.

- [ ] **Step 2: Write RED source-image and recovery assertions (2–5 min)**

Assert source construction uses `sha-${SHORT_SHA}`, never `:main`; `eligible_successful_runs` returns every
completed-success `docker-build.yml` run with exact `head_sha`, `event == "push"`, and `head_branch == "main"` in
descending numeric run-ID order. The workflow checks the exact attempt-qualified artifact on each eligible run and
falls back in that order until it finds one; it downloads only the selected run's
attempt-qualified artifact named `docker-release-evidence-${SHA}-${RUN_ATTEMPT}`, and requires `contract_version == 1`, evidence `run_id/run_attempt`, exact SHA,
`sha256:` digest, revision, and version. Inspection is by `${IMAGE}@${EVIDENCE_DIGEST}`; the short tag is checked for
collision but is not the provenance authority. Missing/ineligible evidence text includes `Docker Build`, exact SHA,
run URL, and “rerun this existing Docker Build run”. Task 8, which creates `build-package`, owns the dependency assertion
that missing evidence/image stops before package build; do not reference a not-yet-created job in Task 6's GREEN suite.

Add an executable workflow-step harness rather than relying only on string presence:

```python
def _run_step(tmp_path: Path, workflow: dict[str, Any], job: str, step_id: str, env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    step = next(step for step in workflow["jobs"][job]["steps"] if step.get("id") == step_id)
    return subprocess.run(
        ["bash", "-euo", "pipefail", "-c", step["run"]],
        cwd=tmp_path,
        env={**os.environ, **env},
        text=True,
        capture_output=True,
        check=False,
    )
```

Prepend a temporary `bin/` containing fake `gh` and `docker` executables. Drive: successful exact-SHA evidence,
missing artifact, wrong contract version, digest mismatch, revision mismatch, version mismatch, short-prefix collision,
skipped component, two transient `gh api` poll failures followed by success, three exhausted API attempts, and policy
timeout. Include two eligible Docker runs where the newest lacks the exact artifact and the older supplies it; assert
the older run ID/attempt/URL is selected. Assert exit status, JSON outputs, exact run URL, and that the fake mutation log stays
empty. This executes the actual checked-in shell with hostile/missing external observations and catches quoting/dataflow
defects that YAML substring checks cannot.

Assert `preflight_summary_json` always has exactly `sha`, `state`, `reason`, `eligible_runs`, `selected_run`, and
`step_outcome`: no completed successful run produces `state="pending"`; candidates without an exact artifact produce
`state="ineligible"`; exhausted API retries produce `state="infrastructure-failure"`; and fallback success produces
`state="eligible"` with the older selected run while retaining both ordered candidates. When preflight deliberately
terminates as `ineligible`, the always-running poll serializer must emit `action="not-run"`, copy the exact preflight
state/reason, and omit the false diagnostic `poll step ended before a Check Runs snapshot`; reserve that diagnostic for
an unexplained missing poll after a nonterminal preflight.

- [ ] **Step 3: Run RED and record absent gate/image job (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_workflow_contract.py -q`
Expected: failures identify missing ten-check polling, finite constants, short-SHA digest inspection, OCI comparisons, and recovery output.

- [ ] **Step 4: Implement bounded Check Runs API snapshots (2–5 min)**

Append this fully guarded job (the Python snippet serializes each API snapshot and policy result rather than parsing JSON
with shell):

```yaml
  wait-for-release-gates:
    needs: validate-release
    runs-on: ubuntu-24.04
    timeout-minutes: 70
    permissions:
      actions: read
      checks: read
      contents: read
      packages: read
    outputs:
      source_ref: ${{ steps.evidence.outputs.source_ref }}
      source_digest: ${{ steps.evidence.outputs.source_digest }}
      source_revision: ${{ steps.evidence.outputs.source_revision }}
      source_version: ${{ steps.evidence.outputs.source_version }}
      preflight_summary_json: ${{ steps.evidence-preflight-result.outputs.preflight_summary_json }}
      check_summary_json: ${{ steps.poll-result.outputs.check_summary_json }}
      evidence_summary_json: ${{ steps.evidence-result.outputs.evidence_summary_json }}
      dry_run_alias_summary_json: ${{ steps.dry-run-aliases.outputs.dry_run_alias_summary_json }}
    env:
      COORDINATOR_ROOT: ${{ github.workspace }}/controller
      RELEASE_POLL_ATTEMPTS: "120"
      RELEASE_POLL_INTERVAL_SECONDS: "30"
      RELEASE_API_ATTEMPTS: "3"
      RELEASE_API_RETRY_SECONDS: "5"
    steps:
      - uses: actions/checkout@v7
        with:
          ref: ${{ github.workflow_sha }}
          persist-credentials: "false"
          path: controller
      - id: evidence-preflight
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
          SHA: ${{ needs.validate-release.outputs.sha }}
        run: |
          set -euo pipefail
          api_json() {
            local endpoint=$1 output=$2
            for API_ATTEMPT in $(seq 1 "$RELEASE_API_ATTEMPTS"); do
              if gh api --paginate --slurp "$endpoint" > "$output.tmp"; then
                mv -- "$output.tmp" "$output"
                return 0
              fi
              rm -f -- "$output.tmp"
              if [ "$API_ATTEMPT" -lt "$RELEASE_API_ATTEMPTS" ]; then sleep "$RELEASE_API_RETRY_SECONDS"; fi
            done
            return 1
          }
          write_preflight_state() {
            python - "$SHA" "$1" "$2" "${3:-}" "${4:-}" "${5:-}" <<'PY' > preflight-state.json
          import json
          import pathlib
          import sys

          sha, state, reason, run_id, run_attempt, run_url = sys.argv[1:]
          eligible_path = pathlib.Path("eligible-runs.json")
          eligible = json.loads(eligible_path.read_text(encoding="utf-8")) if eligible_path.is_file() else []
          selected = None if not run_id else {
              "id": int(run_id), "run_attempt": int(run_attempt), "html_url": run_url or None
          }
          print(json.dumps({"sha": sha, "state": state, "reason": reason,
                            "eligible_runs": eligible, "selected_run": selected}, separators=(",", ":")))
          PY
          }
          if ! api_json \
              "repos/${GITHUB_REPOSITORY}/actions/workflows/docker-build.yml/runs?head_sha=${SHA}&per_page=100" \
              preflight-runs.json; then
            write_preflight_state infrastructure-failure "GitHub workflow-runs API failed after ${RELEASE_API_ATTEMPTS} attempts"
            exit 1
          fi
          python controller/scripts/release_evidence.py eligible-runs preflight-runs.json --sha "$SHA" \
            > eligible-runs.json
          COUNT=$(python -c 'import json; print(len(json.load(open("eligible-runs.json"))))')
          if [ "$COUNT" -eq 0 ]; then
            write_preflight_state pending "no completed successful exact-SHA main-push Docker Build run exists yet"
            exit 0
          fi
          FOUND=false; SELECTED_ID=; SELECTED_ATTEMPT=; SELECTED_URL=
          python -c 'import json; [print("{}\t{}\t{}".format(r["id"], r["run_attempt"], r.get("html_url", ""))) for r in json.load(open("eligible-runs.json"))]' \
            > eligible-runs.tsv
          while IFS=$'\t' read -r RUN_ID RUN_ATTEMPT RUN_URL; do
            if ! api_json "repos/${GITHUB_REPOSITORY}/actions/runs/${RUN_ID}/artifacts?per_page=100" \
                "artifacts-${RUN_ID}.json"; then
              write_preflight_state infrastructure-failure \
                "GitHub artifacts API failed after ${RELEASE_API_ATTEMPTS} attempts for Docker Build run ${RUN_ID}" \
                "$RUN_ID" "$RUN_ATTEMPT" "$RUN_URL"
              exit 1
            fi
            if python controller/scripts/release_evidence.py artifact-present "artifacts-${RUN_ID}.json" \
                --name "docker-release-evidence-${SHA}-${RUN_ATTEMPT}"; then
              FOUND=true; SELECTED_ID=$RUN_ID; SELECTED_ATTEMPT=$RUN_ATTEMPT; SELECTED_URL=$RUN_URL
              break
            fi
          done < eligible-runs.tsv
          if [ "$FOUND" != true ]; then
            write_preflight_state ineligible \
              "successful exact-SHA main-push Docker Build runs exist but none has its exact attempt-qualified contract-v1 artifact"
            echo "Candidate $SHA is pre-contract/ineligible: successful main-push Docker Build runs exist but none " \
                 "contains contract-v1 evidence. Rerun the exact existing Docker Build run after the milestone " \
                 "workflow is present; no production write was attempted." >&2
            exit 1
          fi
          write_preflight_state eligible "exact attempt-qualified Docker evidence artifact found" \
            "$SELECTED_ID" "$SELECTED_ATTEMPT" "$SELECTED_URL"
      - id: evidence-preflight-result
        if: ${{ always() }}
        env:
          PREFLIGHT_OUTCOME: ${{ steps.evidence-preflight.outcome }}
          SHA: ${{ needs.validate-release.outputs.sha }}
        run: |
          set -euo pipefail
          python - "$PREFLIGHT_OUTCOME" "$SHA" <<'PY'
          import json
          import pathlib
          import sys

          path = pathlib.Path("preflight-state.json")
          payload = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {
              "sha": sys.argv[2], "state": "infrastructure-failure",
              "reason": "preflight ended before structured state was written",
              "eligible_runs": [], "selected_run": None,
          }
          payload["step_outcome"] = sys.argv[1]
          path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
          PY
          printf 'preflight_summary_json=%s\n' "$(cat preflight-state.json)" >> "$GITHUB_OUTPUT"
      - id: poll
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
          SHA: ${{ needs.validate-release.outputs.sha }}
        run: |
          set -euo pipefail
          for attempt in $(seq 1 "$RELEASE_POLL_ATTEMPTS"); do
            API_OK=false
            for API_ATTEMPT in $(seq 1 "$RELEASE_API_ATTEMPTS"); do
              if gh api --paginate --slurp \
                  "repos/${GITHUB_REPOSITORY}/commits/${SHA}/check-runs?per_page=100" \
                  > checks-pages.json.tmp; then
                mv -- checks-pages.json.tmp checks-pages.json
                API_OK=true
                break
              fi
              rm -f -- checks-pages.json.tmp
              if [ "$API_ATTEMPT" -lt "$RELEASE_API_ATTEMPTS" ]; then sleep "$RELEASE_API_RETRY_SECONDS"; fi
            done
            if [ "$API_OK" != true ]; then
              python - "$attempt" "$RELEASE_API_ATTEMPTS" <<'PY' > poll.json
          import json
          import sys
          print(json.dumps({"action": "fail", "attempt": int(sys.argv[1]), "elapsed_seconds": 0,
                            "verdicts": [], "infrastructure_error":
                            f"Check Runs API failed after {sys.argv[2]} attempts"}, separators=(",", ":")))
          PY
              cat poll.json >&2
              exit 1
            fi
            PYTHONPATH="$COORDINATOR_ROOT" python - "$SHA" "$attempt" "$RELEASE_POLL_ATTEMPTS" checks-pages.json poll.json <<'PY'
          import dataclasses, json, sys
          from scripts.release_policy import classify_check_runs
          pages = json.load(open(sys.argv[4], encoding="utf-8"))
          runs = [run for page in pages for run in page.get("check_runs", [])]
          poll = classify_check_runs(
              sys.argv[1], runs, attempt=int(sys.argv[2]), max_attempts=int(sys.argv[3])
          )
          with open(sys.argv[5], "w", encoding="utf-8") as handle:
              json.dump(dataclasses.asdict(poll), handle, separators=(",", ":"))
          PY
            ACTION=$(python -c 'import json; print(json.load(open("poll.json"))["action"])')
            if [ "$ACTION" = success ]; then
              exit 0
            fi
            if [ "$ACTION" = fail ] || [ "$ACTION" = timeout ]; then
              cat poll.json >&2
              exit 1
            fi
            sleep "$RELEASE_POLL_INTERVAL_SECONDS"
          done
          echo "release poll exhausted without a terminal classifier result" >&2
          test -s poll.json && cat poll.json >&2
          exit 1
      - id: poll-result
        if: ${{ always() }}
        env:
          POLL_OUTCOME: ${{ steps.poll.outcome }}
        run: |
          set -euo pipefail
          if [ ! -s poll.json ]; then
            python - "$POLL_OUTCOME" <<'PY' > poll.json
          import json
          import pathlib
          import sys

          preflight_path = pathlib.Path("preflight-state.json")
          preflight = json.loads(preflight_path.read_text(encoding="utf-8")) if preflight_path.is_file() else None
          preflight_state = None if preflight is None else preflight.get("state")
          if preflight_state in {"ineligible", "infrastructure-failure"}:
              payload = {
                  "action": "not-run",
                  "attempt": 0,
                  "elapsed_seconds": 0,
                  "verdicts": [],
                  "preflight_state": preflight_state,
                  "preflight_reason": preflight.get("reason"),
                  "step_outcome": sys.argv[1],
              }
          else:
              payload = {
                  "action": "fail",
                  "attempt": 0,
                  "elapsed_seconds": 0,
                  "verdicts": [],
                  "infrastructure_error": "poll step ended before a Check Runs snapshot",
                  "step_outcome": sys.argv[1],
              }
          print(json.dumps(payload, separators=(",", ":")))
          PY
          fi
          printf 'check_summary_json=%s\n' "$(cat poll.json)" >> "$GITHUB_OUTPUT"
      - name: Log in to GHCR for production inspection
        if: ${{ github.event_name == 'repository_dispatch' && github.event.action == 'vntyper_release' }}
        uses: docker/login-action@v4
        with:
          registry: ghcr.io
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}
```

For attempts `1..120`, retry each Check Runs API snapshot at most three times, flatten the page objects'
`check_runs` arrays, pass the combined JSON, attempt, and explicit `max_attempts=120` into `classify_check_runs`, export
one-line `check_summary_json`, exit on `success`/`fail`/`timeout`, and sleep exactly 30 seconds only for `wait`. The shell
also exits nonzero after the loop, so a classifier/configuration disagreement cannot fall through as success. Never
accept aggregate success as replacement for a missing/skipped component. The exact production repository dispatch
authenticates before GHCR reads; manual dispatch deliberately uses the public read path without login and still has no
registry-write authority.

The preflight deliberately precedes the long check poll. A candidate with no completed main-push Docker run is still
possibly building and proceeds to bounded polling. A candidate with completed successful main-push runs but no
attempt-qualified evidence artifact is deterministically pre-contract/ineligible and exits immediately; this is the
required pre-milestone manual dry-run result, not an import failure from checking out the old candidate.

- [ ] **Step 5: Inspect and validate the existing short-SHA image by digest (2–5 min)**

After polling, refresh the eligible-run list because a preflight `pending` candidate may have completed while checks
were polled. In descending run-ID order, select the first run that owns the exact attempt-qualified artifact and run
`gh run download "$RUN_ID" -n "docker-release-evidence-${SHA}-${RUN_ATTEMPT}"`. Retry that exact download at most three
times; never substitute a different artifact name or an unexamined run. Parse the single JSON file with Python;
validate contract version/run identity/SHA/digest/revision/version, then inspect `${IMAGE}@${DIGEST}` and compare OCI
labels. Inspect `sha-${SHORT_SHA}`: matching digest passes; a different full-revision label is a detected prefix
collision and the summary records it, while promotion remains bound to the evidence digest. Missing evidence,
wrong-contract evidence, or missing immutable manifest fails as explicitly ineligible and prints the exact run URL plus
“rerun this existing Docker Build run”. Never fall back to `main` or another SHA.

The evidence step uses this concrete data path (the Python validator also writes the five declared outputs):

```yaml
      - id: evidence
        env:
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
          SHA: ${{ needs.validate-release.outputs.sha }}
          SHORT_SHA: ${{ needs.validate-release.outputs.short_sha }}
          VERSION: ${{ needs.validate-release.outputs.version }}
          IMAGE: ghcr.io/hassansaei/vntyper
        run: |
          set -euo pipefail
          api_json() {
            local endpoint=$1 output=$2
            for API_ATTEMPT in $(seq 1 "$RELEASE_API_ATTEMPTS"); do
              if gh api --paginate --slurp "$endpoint" > "$output.tmp"; then
                mv -- "$output.tmp" "$output"
                return 0
              fi
              rm -f -- "$output.tmp"
              if [ "$API_ATTEMPT" -lt "$RELEASE_API_ATTEMPTS" ]; then sleep "$RELEASE_API_RETRY_SECONDS"; fi
            done
            return 1
          }
          if ! api_json \
              "repos/${GITHUB_REPOSITORY}/actions/workflows/docker-build.yml/runs?head_sha=${SHA}&per_page=100" \
              evidence-runs.json; then
            echo "GitHub workflow-runs API failed after ${RELEASE_API_ATTEMPTS} attempts" >&2
            exit 1
          fi
          python controller/scripts/release_evidence.py eligible-runs evidence-runs.json --sha "$SHA" \
            > evidence-eligible-runs.json
          python -c 'import json; [print("{}\t{}\t{}".format(r["id"], r["run_attempt"], r.get("html_url", ""))) for r in json.load(open("evidence-eligible-runs.json"))]' \
            > evidence-eligible-runs.tsv
          RUN_ID=; RUN_ATTEMPT=; RUN_URL=
          while IFS=$'\t' read -r CANDIDATE_RUN CANDIDATE_ATTEMPT CANDIDATE_URL; do
            if ! api_json "repos/${GITHUB_REPOSITORY}/actions/runs/${CANDIDATE_RUN}/artifacts?per_page=100" \
                "evidence-artifacts-${CANDIDATE_RUN}.json"; then
              echo "GitHub artifacts API failed after ${RELEASE_API_ATTEMPTS} attempts for Docker Build run ${CANDIDATE_RUN}" >&2
              exit 1
            fi
            if python controller/scripts/release_evidence.py artifact-present \
                "evidence-artifacts-${CANDIDATE_RUN}.json" \
                --name "docker-release-evidence-${SHA}-${CANDIDATE_ATTEMPT}"; then
              RUN_ID=$CANDIDATE_RUN; RUN_ATTEMPT=$CANDIDATE_ATTEMPT; RUN_URL=$CANDIDATE_URL
              break
            fi
          done < evidence-eligible-runs.tsv
          if [ -z "$RUN_ID" ]; then
            echo "No eligible exact-SHA Docker evidence artifact; rerun the existing main-push Docker Build run." >&2
            exit 1
          fi
          python - "$RUN_ID" "$RUN_ATTEMPT" "$RUN_URL" <<'PY' > selected-run.json
          import json
          import sys
          print(json.dumps({"id": int(sys.argv[1]), "run_attempt": int(sys.argv[2]),
                            "html_url": sys.argv[3] or None}, separators=(",", ":")))
          PY
          RUN_ID=$(python -c 'import json; print(json.load(open("selected-run.json"))["id"])')
          RUN_ATTEMPT=$(python -c 'import json; print(json.load(open("selected-run.json"))["run_attempt"])')
          DEST=; DOWNLOADED=false
          for API_ATTEMPT in $(seq 1 "$RELEASE_API_ATTEMPTS"); do
            CANDIDATE_DEST="evidence-${RUN_ID}-${RUN_ATTEMPT}-${API_ATTEMPT}"
            if gh run download "$RUN_ID" -n "docker-release-evidence-${SHA}-${RUN_ATTEMPT}" \
                -D "$CANDIDATE_DEST"; then
              DEST=$CANDIDATE_DEST
              DOWNLOADED=true
              break
            fi
            if [ "$API_ATTEMPT" -lt "$RELEASE_API_ATTEMPTS" ]; then sleep "$RELEASE_API_RETRY_SECONDS"; fi
          done
          if [ "$DOWNLOADED" != true ]; then
            echo "Exact Docker evidence download failed after ${RELEASE_API_ATTEMPTS} attempts for run ${RUN_ID}; rerun this existing Docker Build run." >&2
            exit 1
          fi
          EVIDENCE_PATH="$DEST/docker-release-evidence.json"
          python controller/scripts/release_evidence.py validate-evidence "$EVIDENCE_PATH" \
            --sha "$SHA" --version "$VERSION" --run-id "$RUN_ID" --run-attempt "$RUN_ATTEMPT" \
            > validated-evidence.json
          DIGEST=$(python -c 'import json; print(json.load(open("validated-evidence.json"))["digest"])')
          IMMUTABLE_DIGEST=$(docker buildx imagetools inspect "$IMAGE@$DIGEST" --format '{{.Manifest.Digest}}')
          test "$IMMUTABLE_DIGEST" = "$DIGEST"
          docker buildx imagetools inspect "$IMAGE@$DIGEST" --format '{{json .Image}}' > image-config.json
          SHORT_DIGEST=$(docker buildx imagetools inspect "$IMAGE:sha-$SHORT_SHA" --format '{{.Manifest.Digest}}')
          docker buildx imagetools inspect "$IMAGE:sha-$SHORT_SHA" --format '{{json .Image}}' > short-image-config.json
          python controller/scripts/release_evidence.py validate-image validated-evidence.json image-config.json \
            short-image-config.json --immutable-digest "$IMMUTABLE_DIGEST" --short-digest "$SHORT_DIGEST" \
            --image "$IMAGE" --output "$GITHUB_OUTPUT"
      - id: evidence-result
        if: ${{ always() }}
        env:
          EVIDENCE_OUTCOME: ${{ steps.evidence.outcome }}
          SHA: ${{ needs.validate-release.outputs.sha }}
        run: |
          set -euo pipefail
          python controller/scripts/release_evidence.py summarize \
            --sha "$SHA" --step-outcome "$EVIDENCE_OUTCOME" --selected-run selected-run.json \
            --validated-evidence validated-evidence.json --image image-config.json \
            --output evidence-summary.json
          printf 'evidence_summary_json=%s\n' "$(cat evidence-summary.json)" >> "$GITHUB_OUTPUT"
```

Create the named typed adapter in Task 6 with pure JSON-in/value-out functions and a marked unit file; its CLI reads only
the named local files. Import `Mapping` and `Sequence` from `collections.abc`, not `typing`.
`eligible_successful_runs` rejects non-success/wrong-SHA/non-main/non-push runs, requires numeric `id` and
`run_attempt`, preserves `html_url`, and returns every eligible mapping ordered by descending numeric run ID;
`validate_evidence` requires contract-v1/exact run/SHA/version/revision/digest; `validate_image` requires config labels,
an immutable-manifest digest equal to the evidence digest, and either a matching short-tag digest or a proven prefix
collision without changing the evidence digest. A mismatch is a prefix collision only when the short image's own
revision label is a different 40-hex SHA sharing exactly the candidate's first seven characters; missing, unrelated,
or same-revision/different-digest labels are unexplained drift and fail. The CLI passes the two registry-reported
manifest digests into the pure
validator. Include it in Task 6
mypy, Ruff, coverage, and commit.

Write these behavioral tests before the adapter. They use complete local payloads and assert values, not mere return:

```python
import pytest

from scripts.release_evidence import eligible_successful_runs, validate_evidence, validate_image

pytestmark = pytest.mark.unit
SHA = "a" * 40
DIGEST = "sha256:" + "b" * 64


def test_evidence_binds_successful_exact_sha_run_and_image_labels() -> None:
    runs = {
        "workflow_runs": [
            {"id": 41, "run_attempt": 2, "head_sha": SHA, "head_branch": "main", "event": "push",
             "status": "completed", "conclusion": "success"},
            {"id": 99, "run_attempt": 1, "head_sha": SHA, "head_branch": "main", "event": "schedule",
             "status": "completed", "conclusion": "success"},
        ]
    }
    eligible = eligible_successful_runs(runs, SHA)
    assert tuple(run["id"] for run in eligible) == (41,)
    run = eligible[0]
    assert (run["id"], run["run_attempt"]) == (41, 2)
    evidence = validate_evidence(
        {"contract_version": 1, "sha": SHA, "digest": DIGEST, "run_id": 41, "run_attempt": 2,
         "revision": SHA, "version": "2.0.10"},
        sha=SHA,
        version="2.0.10",
        run_id=41,
        run_attempt=2,
    )
    image = {"config": {"Labels": {"org.opencontainers.image.revision": SHA,
                                     "org.opencontainers.image.version": "2.0.10"}}}
    result = validate_image(evidence, image, image, immutable_digest=DIGEST, short_digest=DIGEST)
    assert result == {"source_digest": DIGEST, "source_revision": SHA, "source_version": "2.0.10",
                      "short_tag_collision": False}


def test_equal_short_prefix_with_other_digest_is_reported_not_substituted() -> None:
    evidence = {"digest": DIGEST, "revision": SHA, "version": "2.0.10"}
    image = {"config": {"Labels": {"org.opencontainers.image.revision": SHA,
                                     "org.opencontainers.image.version": "2.0.10"}}}
    other_sha = SHA[:7] + "d" * 33
    short_image = {"config": {"Labels": {"org.opencontainers.image.revision": other_sha,
                                           "org.opencontainers.image.version": "2.0.9"}}}
    result = validate_image(
        evidence,
        image,
        short_image,
        immutable_digest=DIGEST,
        short_digest="sha256:" + "d" * 64,
    )
    assert result["source_digest"] == DIGEST
    assert result["short_tag_collision"] is True
```

Also parameterize wrong `contract_version`, run ID, SHA, revision, version, malformed digest, missing labels, and
immutable-manifest digest. Add short-tag mismatch rows for absent revision, unrelated seven-character prefix, and
candidate revision with different digest; each raises `ValueError("unexplained short-tag drift")`. Each case raises
`ValueError` naming the mismatched field. The CLI test invokes `main([...])`
with `tmp_path` JSON files, captures `$GITHUB_OUTPUT`, and asserts the exact five declared one-line outputs.
Add eligible-run rows in shuffled order with IDs `40`, `42`, and `41`; assert the result is `(42, 41, 40)`, then execute
the workflow harness with run 42 missing its exact artifact and run 41 owning it and assert run 41 is selected. No
singular newest-only selector, import, or call remains.

- [ ] **Step 6: Run GREEN and focused policy regression (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_policy.py tests/unit/test_release_evidence.py tests/unit/test_release_workflow_contract.py -q && ruff check scripts/release_evidence.py tests/unit/test_release_evidence.py && mypy scripts/release_evidence.py && make lint-actions`
Expected: exact ten-check polling, bounded API retries, skipped rejection, finite timeout/post-loop refusal, ordered
artifact fallback, structured preflight state, authenticated production image inspection, image identity, and
missing-image recovery contracts pass.

- [ ] **Step 7: Commit release evidence gate (2–5 min)**

```bash
git add .github/workflows/publish-pypi.yml scripts/release_evidence.py tests/unit/test_release_evidence.py tests/unit/test_release_workflow_contract.py
git commit -m "ci(release): wait for exact-sha artifacts"
```

### Task 7: Globally serialized digest promotion and reruns

**Files:**
- Modify: `.github/workflows/publish-pypi.yml`
- Modify: `tests/unit/test_release_workflow_contract.py`

**Interfaces:**
- Consumes: verified source digest/version and `plan_alias_updates` observations.
- Produces: GHCR aliases `vX.Y.Z`, `X.Y.Z`, `X.Y`, `X`, `latest` pointing to the tested digest plus structured alias summary.

- [ ] **Step 1: Write RED concurrency and digest-source tests (2–5 min)**

Assert `promote-ghcr.concurrency.group == "vntyper-ghcr-release-promotion"`, `cancel-in-progress == "false"`, permissions are exactly `contents: read` and `packages: write`, and every production
`docker buildx imagetools create` contains `--prefer-index=false`, tags exactly one approved alias, and uses
`${IMAGE}@${SOURCE_DIGEST}` as its sole source, never a tag or local daemon image.

Also assert the workflow header documents that the fixed group is a lock, not an unbounded queue: a superseded pending
promotion performs no writes and must be rerun. Do not assert that three simultaneous releases all remain queued.

- [ ] **Step 2: Write RED exact/floating/dry-run workflow tests (2–5 min)**

Assert the job calls `plan_alias_updates`; aliases are exactly the required five; `main` is not a target;
`fail-conflict` (including floating equal-version/different-digest) exits nonzero before any write; `skip-newer` and
`skip-unorderable` emit `::notice::`; only `create`/`advance` execute. Assert `promote-ghcr.if` contains
exactly `github.event_name == 'repository_dispatch' && github.event.action == 'vntyper_release'`, so manual or unrelated
dispatches cannot log in or write.

Use the Task 6 `_run_step` harness with a fake `docker` that returns controlled alias digests/labels and records every
`imagetools create`. Make the fake reject a production create lacking `--prefer-index=false` and execute the checked-in
promotion step; do not satisfy this contract with a YAML substring assertion alone. Assert the
equal-version/different-digest case exits nonzero with an empty mutation log; an older
candidate records skip notices; a normal five-alias plan records only `${IMAGE}@${SOURCE_DIGEST}` sources; rerunning
after every prefix produces no conflicting write. Force reinspection to fail after a successful create and assert the
atomically written progress row says `attempted=true`, `write_succeeded=true`, `verified=false`; a fully successful row
says all three are true and records the final digest. These are behavior tests of the actual shell, not source-string tests.

Execute the dispatch `dry-run-aliases` step with the same fake and require exactly one read-only
`imagetools create --dry-run --prefer-index=false ${IMAGE}@${SOURCE_DIGEST}` probe with no `--tag`, plus zero
production create records. This proves Buildx accepts the carbon-copy contract while preserving the no-write boundary.

- [ ] **Step 3: Run RED and confirm promotion is absent (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_workflow_contract.py -q`
Expected: failures name missing fixed concurrency, five aliases, policy call, digest source, and notices.

- [ ] **Step 4: Implement read-plan-write under the fixed global lock (2–5 min)**

Append this job and use a local observation directory plus the typed policy. Raw manifest bytes determine digests;
`.Image` JSON supplies the package-version label. No registry response is evaluated as shell:

```yaml
  promote-ghcr:
    needs: [validate-release, wait-for-release-gates]
    if: ${{ github.event_name == 'repository_dispatch' && github.event.action == 'vntyper_release' }}
    runs-on: ubuntu-24.04
    concurrency:
      group: vntyper-ghcr-release-promotion
      cancel-in-progress: false
    permissions:
      contents: read
      packages: write
    outputs:
      alias_summary_json: ${{ steps.promote-result.outputs.alias_summary_json }}
    steps:
      - uses: actions/checkout@v7
        with:
          ref: ${{ github.workflow_sha }}
          persist-credentials: "false"
          path: controller
      - uses: docker/login-action@v4
        with:
          registry: ghcr.io
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}
      - id: promote
        env:
          IMAGE: ghcr.io/hassansaei/vntyper
          SOURCE_DIGEST: ${{ needs.wait-for-release-gates.outputs.source_digest }}
          VERSION: ${{ needs.validate-release.outputs.version }}
        run: |
          set -euo pipefail
          OBS_DIR=$(mktemp -d)
          trap 'rm -rf -- "$OBS_DIR"' EXIT
          for ALIAS in "v$VERSION" "$VERSION" "${VERSION%.*}" "${VERSION%%.*}" latest; do
            SAFE_ALIAS=${ALIAS//./_}
            if docker buildx imagetools inspect "$IMAGE:$ALIAS" --format '{{.Manifest.Digest}}' \
                > "$OBS_DIR/$SAFE_ALIAS.digest"; then
              docker buildx imagetools inspect "$IMAGE:$ALIAS" --format '{{json .Image}}' \
                > "$OBS_DIR/$SAFE_ALIAS.image.json"
            else
              rm -f -- "$OBS_DIR/$SAFE_ALIAS.digest"
            fi
          done
          python - "$VERSION" "$SOURCE_DIGEST" "$OBS_DIR" > plan.json <<'PY'
          import dataclasses
          import json
          import pathlib
          import sys
          sys.path.insert(0, str(pathlib.Path("controller").resolve()))
          from scripts.release_policy import AliasState, parse_release_tag, plan_alias_updates, required_aliases

          version = parse_release_tag("v" + sys.argv[1])
          root = pathlib.Path(sys.argv[3])
          observed = {}
          for alias in required_aliases(version):
              stem = alias.replace(".", "_")
              digest_path = root / f"{stem}.digest"
              if not digest_path.exists():
                  observed[alias] = None
                  continue
              image = json.loads((root / f"{stem}.image.json").read_text(encoding="utf-8"))
              labels = image.get("config", {}).get("Labels", {}) or {}
              digest = digest_path.read_text(encoding="utf-8").strip()
              if not digest.startswith("sha256:") or len(digest) != 71:
                  raise SystemExit(f"invalid registry digest for {alias}: {digest}")
              observed[alias] = AliasState(digest, labels.get("org.opencontainers.image.version"))
          plan = plan_alias_updates(version, sys.argv[2], observed, dry_run=False)
          summary = {
              "source_digest": sys.argv[2],
              "observed": {key: None if value is None else dataclasses.asdict(value)
                           for key, value in observed.items()},
              "plan": [dataclasses.asdict(item) for item in plan],
          }
          print(json.dumps(summary, separators=(",", ":")))
          PY
          python - plan.json > executable.tsv <<'PY'
          import json
          import sys

          plan = json.load(open(sys.argv[1], encoding="utf-8"))["plan"]
          conflicts = [item for item in plan if item["decision"] == "fail-conflict"]
          if conflicts:
              raise SystemExit("release alias conflict: " + json.dumps(conflicts, separators=(",", ":")))
          for item in plan:
              if item["decision"] in {"skip-newer", "skip-unorderable", "no-op"}:
                  print(f"::notice title=GHCR alias {item['alias']}::{item['reason']}", file=sys.stderr)
              if item["execute"]:
                  print(f"{item['alias']}\t{item['decision']}")
          PY
          printf '[]' > alias-progress.json
          record_progress() {
            python - alias-progress.json plan.json "$1" "$2" "$3" <<'PY'
          import json
          import os
          import sys
          import tempfile

          path, plan_path, alias, phase, final_digest = sys.argv[1:]
          progress = json.load(open(path, encoding="utf-8"))
          plan_payload = json.load(open(plan_path, encoding="utf-8"))
          planned = next(item for item in plan_payload["plan"] if item["alias"] == alias)
          observed = plan_payload["observed"].get(alias)
          row = next((item for item in progress if item["alias"] == alias), None)
          if row is None:
              row = {
                  "alias": alias,
                  "decision": planned["decision"],
                  "reason": planned["reason"],
                  "previous_digest": None if observed is None else observed["digest"],
                  "previous_version": None if observed is None else observed["version"],
                  "attempted": False,
                  "write_succeeded": False,
                  "verified": False,
                  "final_digest": None,
              }
              progress.append(row)
          if phase == "attempted":
              row["attempted"] = True
          elif phase == "write-succeeded":
              row["attempted"] = True
              row["write_succeeded"] = True
          elif phase == "verified":
              row["attempted"] = True
              row["write_succeeded"] = True
              row["verified"] = True
              row["final_digest"] = final_digest
          else:
              raise SystemExit(f"unknown promotion progress phase: {phase}")
          fd, temp = tempfile.mkstemp(prefix=".alias-progress.", dir=".")
          try:
              with os.fdopen(fd, "w", encoding="utf-8") as handle:
                  json.dump(progress, handle, separators=(",", ":"))
                  handle.flush()
                  os.fsync(handle.fileno())
              os.replace(temp, path)
          finally:
              if os.path.exists(temp):
                  os.unlink(temp)
          PY
          }
          while IFS=$'\t' read -r ALIAS DECISION; do
            case "$ALIAS" in "v$VERSION"|"$VERSION"|"${VERSION%.*}"|"${VERSION%%.*}"|latest) ;; *) exit 2 ;; esac
            case "$DECISION" in create|advance) ;; *) exit 2 ;; esac
            record_progress "$ALIAS" attempted ""
            docker buildx imagetools create --prefer-index=false --tag "$IMAGE:$ALIAS" "$IMAGE@$SOURCE_DIGEST"
            record_progress "$ALIAS" write-succeeded ""
            FINAL_DIGEST=$(docker buildx imagetools inspect "$IMAGE:$ALIAS" --format '{{.Manifest.Digest}}')
            test "$FINAL_DIGEST" = "$SOURCE_DIGEST"
            record_progress "$ALIAS" verified "$FINAL_DIGEST"
          done < executable.tsv
      - id: promote-result
        if: ${{ always() }}
        env:
          PROMOTE_OUTCOME: ${{ steps.promote.outcome }}
        run: |
          set -euo pipefail
          python - "$PROMOTE_OUTCOME" <<'PY' > alias-summary.json
          import json, pathlib, sys
          def load(path, fallback):
              file = pathlib.Path(path)
              return json.loads(file.read_text(encoding="utf-8")) if file.is_file() else fallback
          print(json.dumps({"step_outcome": sys.argv[1], "plan": load("plan.json", {}),
                            "alias_progress": load("alias-progress.json", [])}, separators=(",", ":")))
          PY
          printf 'alias_summary_json=%s\n' "$(cat alias-summary.json)" >> "$GITHUB_OUTPUT"
```

Inspect every target alias before mutation, record digest and package-version label or `None`, call
`plan_alias_updates(..., dry_run=False)`, and fail before writes if any `fail-conflict` exists. Execute exact aliases in
tuple order followed by eligible floating aliases using
`docker buildx imagetools create --prefer-index=false --tag "$IMAGE:$ALIAS" "$IMAGE@$SOURCE_DIGEST"`. Parse the policy's JSON output in
Python or `jq -er`; never `eval` it as shell.

- [ ] **Step 5: Make partial reruns observable and convergent (2–5 min)**

Export single-line `alias_summary_json` containing previous digest/version, decision, reason, and final expected digest
for all five aliases. Before each create, atomically persist `attempted=true`; immediately after Buildx returns zero,
atomically persist `write_succeeded=true`; only after digest reinspection atomically persist `verified=true` and the
final digest. Emit notices for anti-downgrade/unorderable skips and no-ops. Reinspect changed aliases and
hard-fail if any final digest differs from the verified source digest. The final summary job, not this job's local step
summary, owns the cross-job report.

For `workflow_dispatch`, perform the same five read-only alias inspections in
`wait-for-release-gates`, call `plan_alias_updates(..., dry_run=True)`, and summarize the plan.
Do not enter `promote-ghcr`, acquire write permission, log in for write, tag an image, or execute any plan item.
Append this concrete read-only step to `wait-for-release-gates`; it deliberately has no login or mutating create
command. Its untagged Buildx `--dry-run` probe exercises the single-manifest carbon-copy contract without a registry
write:

```yaml
      - id: dry-run-aliases
        if: ${{ github.event_name == 'workflow_dispatch' }}
        env:
          IMAGE: ghcr.io/hassansaei/vntyper
          SOURCE_DIGEST: ${{ steps.evidence.outputs.source_digest }}
          VERSION: ${{ needs.validate-release.outputs.version }}
        run: |
          set -euo pipefail
          OBS_DIR=$(mktemp -d)
          trap 'rm -rf -- "$OBS_DIR"' EXIT
          for ALIAS in "v$VERSION" "$VERSION" "${VERSION%.*}" "${VERSION%%.*}" latest; do
            SAFE_ALIAS=${ALIAS//./_}
            if docker buildx imagetools inspect "$IMAGE:$ALIAS" --format '{{.Manifest.Digest}}' \
                > "$OBS_DIR/$SAFE_ALIAS.digest"; then
              docker buildx imagetools inspect "$IMAGE:$ALIAS" --format '{{json .Image}}' \
                > "$OBS_DIR/$SAFE_ALIAS.image.json"
            else
              rm -f -- "$OBS_DIR/$SAFE_ALIAS.digest"
            fi
          done
          docker buildx imagetools create --dry-run --prefer-index=false \
            "$IMAGE@$SOURCE_DIGEST" > "$OBS_DIR/source-carbon-copy.json"
          test -s "$OBS_DIR/source-carbon-copy.json"
          SUMMARY=$(python - "$VERSION" "$SOURCE_DIGEST" "$OBS_DIR" <<'PY'
          import dataclasses
          import json
          import pathlib
          import sys
          sys.path.insert(0, str(pathlib.Path("controller").resolve()))
          from scripts.release_policy import AliasState, parse_release_tag, plan_alias_updates, required_aliases

          version = parse_release_tag("v" + sys.argv[1])
          root = pathlib.Path(sys.argv[3])
          observed = {}
          for alias in required_aliases(version):
              stem = alias.replace(".", "_")
              digest_path = root / f"{stem}.digest"
              if not digest_path.exists():
                  observed[alias] = None
                  continue
              image = json.loads((root / f"{stem}.image.json").read_text(encoding="utf-8"))
              labels = image.get("config", {}).get("Labels", {}) or {}
              digest = digest_path.read_text(encoding="utf-8").strip()
              if not digest.startswith("sha256:") or len(digest) != 71:
                  raise SystemExit(f"invalid registry digest for {alias}: {digest}")
              observed[alias] = AliasState(digest, labels.get("org.opencontainers.image.version"))
          plan = plan_alias_updates(version, sys.argv[2], observed, dry_run=True)
          if any(item.execute for item in plan):
              raise SystemExit("dry-run policy attempted a write")
          print(json.dumps({"source_digest": sys.argv[2],
                            "observed": {key: None if value is None else dataclasses.asdict(value)
                                         for key, value in observed.items()},
                            "plan": [dataclasses.asdict(item) for item in plan]}, separators=(",", ":")))
          PY
          )
          printf 'dry_run_alias_summary_json=%s\n' "$SUMMARY" >> "$GITHUB_OUTPUT"
```

Change the job output to `${{ steps.dry-run-aliases.outputs.dry_run_alias_summary_json }}`. The workflow contract test
executes this checked-in step with the fake Docker command and asserts exactly one untagged
`imagetools create --dry-run --prefer-index=false` record and zero mutating create records.

- [ ] **Step 6: Run GREEN plus policy prefix-rerun suite (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_policy.py tests/unit/test_release_workflow_contract.py -q && make lint-actions`
Expected: global serialization, exact immutability, `--prefer-index=false` production and read-only dry-run execution,
monotonic floating aliases, same-digest no-ops, atomic three-phase promotion progress, dry-run exclusion, and post-write
verification pass.

- [ ] **Step 7: Commit promotion (2–5 min)**

```bash
git add .github/workflows/publish-pypi.yml tests/unit/test_release_workflow_contract.py
git commit -m "ci(release): promote tested GHCR digests"
```

### Task 8: Unprivileged package build and OIDC-only publication

**Files:**
- Modify: `.github/workflows/publish-pypi.yml`
- Modify: `tests/unit/test_release_workflow_contract.py`

**Interfaces:**
- Consumes: validated SHA, successful gates, verified image, successful promotion, built `dist/` artifact.
- Produces: PyPI wheel/sdist and attestations through environment `pypi`; partial reruns retain skip-existing behavior.

- [ ] **Step 1: Write RED privilege-separation tests (2–5 min)**

Assert `build-package` has `contents: read`, no `id-token`, environment, package write, or secret; checkout has
`persist-credentials: false`; `actions/setup-python@v7` selects `3.12`; it runs `python -m build` and
`twine check dist/*`. Assert its upload is `actions/upload-artifact@v5`, path `dist/`, and `retention-days: "7"`.
Assert `publish-pypi` downloads that exact name with `actions/download-artifact@v5`, has
`environment.name == "pypi"`, only `id-token: write`, no checkout/build/install step, and needs `validate-release`,
`wait-for-release-gates`, `build-package`, and `promote-ghcr`.

- [ ] **Step 2: Write RED trusted-publisher and rerun tests (2–5 min)**

Assert literal action `pypa/gh-action-pypi-publish@dc37677b2e1c63e2034f94d8a5b11f265b73ba33`, adjacent comment
`# v1.14.2; verified 2026-08-10 against upstream release/v1`, `skip-existing: "true"`, and no occurrence of
`PYPI_API_TOKEN`, `TWINE_PASSWORD`, `TWINE_USERNAME`, `DOCKER_USERNAME`, or `DOCKER_PASSWORD`. Inspect that PyPA
step's own `with` mapping has no `password` or `user`; a separate GHCR login may use
`password: ${{ secrets.GITHUB_TOKEN }}`. Assert build and publish are downstream of image evidence; publish is
guarded by exact `repository_dispatch` event name and `vntyper_release` action only.

- [ ] **Step 3: Run RED against the remaining token path (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_workflow_contract.py -q`
Expected: failures identify combined build/publish privilege, live token literals, absent environment, unpinned action, and missing job dependencies.

- [ ] **Step 4: Implement the package artifact job (2–5 min)**

Implement the job with this exact identity/dataflow:

```yaml
  build-package:
    needs: [validate-release, wait-for-release-gates]
    runs-on: ubuntu-24.04
    permissions:
      contents: read
    outputs:
      artifact_name: ${{ steps.package.outputs.artifact_name }}
      package_summary_json: ${{ steps.package-result.outputs.package_summary_json }}
    steps:
      - uses: actions/checkout@v7
        with:
          ref: ${{ needs.validate-release.outputs.sha }}
          persist-credentials: "false"
          path: candidate
      - uses: actions/setup-python@v7
        with:
          python-version: "3.12"
      - id: package
        working-directory: candidate
        run: |
          set -euo pipefail
          python3 -m venv .release-venv
          .release-venv/bin/pip install build twine
          .release-venv/bin/python -m build
          .release-venv/bin/twine check dist/*
          ARTIFACT_NAME="python-dist-${{ needs.validate-release.outputs.version }}-${{ github.run_id }}-${{ github.run_attempt }}"
          PACKAGE_SUMMARY_JSON=$(python - "$ARTIFACT_NAME" dist <<'PY'
          import json
          import pathlib
          import sys
          files = sorted(path.name for path in pathlib.Path(sys.argv[2]).iterdir() if path.is_file())
          if len(files) != 2 or not any(name.endswith(".whl") for name in files) or not any(name.endswith(".tar.gz") for name in files):
              raise SystemExit(f"expected one wheel and one sdist, got {files!r}")
          print(json.dumps({"artifact_name": sys.argv[1], "files": files}, separators=(",", ":")))
          PY
          )
          printf 'artifact_name=%s\n' "$ARTIFACT_NAME" >> "$GITHUB_OUTPUT"
      - uses: actions/upload-artifact@v5
        with:
          name: ${{ steps.package.outputs.artifact_name }}
          path: candidate/dist/
          if-no-files-found: error
          retention-days: "7"
      - id: package-result
        if: ${{ always() }}
        env:
          ARTIFACT_NAME: ${{ steps.package.outputs.artifact_name }}
          PACKAGE_OUTCOME: ${{ steps.package.outcome }}
        run: |
          python - <<'PY' >> "$GITHUB_OUTPUT"
          import json, os, pathlib
          root = pathlib.Path("candidate/dist")
          files = sorted(path.name for path in root.iterdir() if path.is_file()) if root.is_dir() else []
          payload = {"artifact_name": os.environ["ARTIFACT_NAME"] or None, "files": files,
                     "step_outcome": os.environ["PACKAGE_OUTCOME"]}
          print("package_summary_json=" + json.dumps(payload, separators=(",", ":")))
          PY
```

List exact wheel/sdist names in `package_summary_json`. Add `build-package` to `promote-ghcr.needs`, preserving
`validate-release` and `wait-for-release-gates`, so promotion waits for a valid package artifact.

- [ ] **Step 5: Implement protected OIDC publish with skip-existing (2–5 min)**

Download the exact `needs.build-package.outputs.artifact_name` in this job; do not use a version-only wildcard:

```yaml
  publish-pypi:
    needs: [validate-release, wait-for-release-gates, build-package, promote-ghcr]
    if: ${{ github.event_name == 'repository_dispatch' && github.event.action == 'vntyper_release' }}
    runs-on: ubuntu-24.04
    environment:
      name: pypi
    permissions:
      id-token: write
    outputs:
      publish_summary_json: ${{ steps.result.outputs.publish_summary_json }}
    steps:
      - uses: actions/download-artifact@v5
        with:
          name: ${{ needs.build-package.outputs.artifact_name }}
          path: dist/
      - id: existing
        env:
          VERSION: ${{ needs.validate-release.outputs.version }}
        run: |
          set -euo pipefail
          STATUS=$(curl --silent --show-error --output /dev/null --write-out '%{http_code}' \
            "https://pypi.org/pypi/vntyper/${VERSION}/json")
          case "$STATUS" in 200) EXISTS=true ;; 404) EXISTS=false ;; *) echo "PyPI lookup returned HTTP $STATUS" >&2; exit 1 ;; esac
          echo "value=$EXISTS" >> "$GITHUB_OUTPUT"
      # v1.14.2; verified 2026-08-10 against upstream release/v1
      - id: publish
        uses: pypa/gh-action-pypi-publish@dc37677b2e1c63e2034f94d8a5b11f265b73ba33
        with:
          packages-dir: dist/
          skip-existing: "true"
      - id: result
        if: ${{ always() }}
        env:
          EXISTED_BEFORE: ${{ steps.existing.outputs.value }}
          PUBLISH_OUTCOME: ${{ steps.publish.outcome }}
        run: |
          python - "$EXISTED_BEFORE" "$PUBLISH_OUTCOME" <<'PY' >> "$GITHUB_OUTPUT"
          import json, sys
          existed = sys.argv[1] == "true"
          outcome = "already-existed-skip" if existed and sys.argv[2] == "success" else (
              "published" if not existed and sys.argv[2] == "success" else "failed"
          )
          print("publish_summary_json=" + json.dumps({"result": outcome, "step_outcome": sys.argv[2],
                                                       "existed_before": existed}, separators=(",", ":")))
          PY
```

Do not check out source or expose any repository/package secret.

- [ ] **Step 6: Run GREEN, actionlint, and secret-literal scan (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_workflow_contract.py -q && make lint-actions && ! rg -n 'PYPI_API_TOKEN|TWINE_PASSWORD|DOCKER_(USERNAME|PASSWORD)' .github/workflows/publish-pypi.yml`
Expected: permission graph, pinned action, OIDC environment, skip-existing, and absence scan all pass.

- [ ] **Step 7: Commit OIDC migration (2–5 min)**

```bash
git add .github/workflows/publish-pypi.yml tests/unit/test_release_workflow_contract.py
git commit -m "ci(release): publish packages through PyPI OIDC"
```

### Task 9: GHCR-only documentation and nine product-generation edits

**Files:**
- Modify: `README.md`
- Modify: `docs/getting-started/installation.md`
- Modify: `docs/user-guide/docker.md`
- Modify: `docker/Dockerfile`
- Modify: `tests/unit/test_release_workflow_contract.py`

**Interfaces:**
- Consumes: supported aliases from `required_aliases`; rolling `main` behavior from Docker workflow.
- Produces: GHCR-only current commands and exactly nine `VNtyper 2` generation targets without changing protected versions/names.

- [ ] **Step 1: Write RED GHCR documentation assertions (2–5 min)**

Read the three install surfaces and assert active `docker pull`, `docker run`, and `apptainer pull` lines contain
`ghcr.io/hassansaei/vntyper`; assert no active command matches `(?:docker://)?saei/vntyper`. Until a gated release has
created `latest`, runnable examples must use `:main` and call it rolling/development. Prose may describe `latest` and
exact `vX.Y.Z`/`X.Y.Z` as the post-release stable/immutable aliases, but must explicitly say they become available only
after the first gated release. This makes every command truthful at PR merge.

- [ ] **Step 2: Write RED nine-target and protected-identity assertions (2–5 min)**

Use this exact baseline-to-required mapping rather than line numbers:

```python
README_RENAMES = {
    "# VNtyper 2.0 - A Pipeline": "# VNtyper 2 - A Pipeline",
    "**VNtyper 2.0** is an advanced": "**VNtyper 2** is an advanced",
    "VNtyper 2.0 uses modern Python packaging": "VNtyper 2 uses modern Python packaging",
    "VNtyper 2.0 offers multiple subcommands": "VNtyper 2 offers multiple subcommands",
    "Docker image for VNtyper 2.0 is provided": "Docker image for VNtyper 2 is provided",
    "VNtyper 2.0 integrates multiple steps": "VNtyper 2 integrates multiple steps",
    "VNtyper 2.0 relies on several tools": "VNtyper 2 relies on several tools",
    "If you use VNtyper 2.0 in your research": "If you use VNtyper 2 in your research",
}

def test_generation_renames_and_grammar_are_exact() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    for before, after in README_RENAMES.items():
        assert before not in readme
        assert readme.count(after) == 1
    assert "This version is a refactored version of VNtyper v1 integrates" not in readme
    assert "This refactored version of VNtyper v1 integrates" in readme
```

Add the Dockerfile description assertion as the ninth target. Also assert exact protected literals
`title: "VNtyper"`, `site_name: VNtyper`, `## VNtyper Version:`, and both `snakemake/vntyper2.smk` filenames remain.

- [ ] **Step 3: Run RED and capture current stale registry/naming failures (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_workflow_contract.py -q`
Expected: failures list active `saei/vntyper` commands and all nine `VNtyper 2.0` presentation targets.

- [ ] **Step 4: Replace active Docker Hub commands with GHCR and explain aliases (2–5 min)**

Use `ghcr.io/hassansaei/vntyper:main` for all runnable container examples at PR merge and state it is
rolling/unreleased. Document the future stable `latest` and exact semantic aliases as unavailable until the first gated
post-milestone release; after that release a docs-only follow-up may switch stable examples to `latest`. Remove all
current Docker Hub pull/run/Apptainer commands; one prose transition note may state Docker Hub artifacts are legacy,
frozen, and unsupported, with no installation command.

- [ ] **Step 5: Hand-edit exactly the nine generation targets and grammar (2–5 min)**

Apply only the eight literal `README_RENAMES` pairs above and the Dockerfile description. Do not use line numbers or a
global substitution: the earlier registry edit changes line offsets. Make the opening read “This refactored version of
VNtyper v1 integrates …” and preserve every numeric version elsewhere.

- [ ] **Step 6: Run GREEN plus version and docs checks (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_workflow_contract.py tests/unit/test_version_consistency.py -q && make docs-build`
Expected: GHCR/naming/protected-identity tests pass and strict docs build succeeds.

- [ ] **Step 7: Inspect the zero-context diff for forbidden version edits (2–5 min)**

Run: `git diff -U0 -- README.md docs/getting-started/installation.md docs/user-guide/docker.md docker/Dockerfile | rg '^[+-]'` and `git diff -U0 | rg '^[+-].*\b2\.0\.[0-9]+\b'`
Expected: first output contains the intended registry, nine generation, and grammar lines; second command exits 1 with no matching version change.

- [ ] **Step 8: Commit docs and naming contract (2–5 min)**

```bash
git add README.md docs/getting-started/installation.md docs/user-guide/docker.md docker/Dockerfile tests/unit/test_release_workflow_contract.py
git commit -m "docs(release): make GHCR and VNtyper 2 canonical"
```

### Task 10: Contributor guidance, migration record, and release diagnostics

**Files:**
- Modify: `AGENTS.md`
- Modify: `docs/development/ci-followups.md`
- Modify: `.github/workflows/publish-pypi.yml` header and summary steps
- Modify: `tests/unit/test_release_workflow_contract.py`

**Interfaces:**
- Consumes: final job names, permissions, failure/recovery rules, rollout boundary.
- Produces: accurate agent/maintainer guidance and an always-run structured release summary.

- [ ] **Step 1: Write RED guidance and summary literal tests (2–5 min)**

Assert `AGENTS.md` names the ten checks, short-SHA/full-revision distinction, digest promotion, `main` versus `latest`, manual existing-tag dry run, fixed promotion serialization, OIDC environment, and “do not delete `PYPI_API_TOKEN` until the first successful OIDC release”. Assert B4 says `RESOLVED` for workflow migration and retains owner-only token deletion as follow-up.
Also require the exact `vntyper_release` repository-dispatch type/payload, default-branch controller semantics, absence of
production `push.tags`, exact event/action guards on every production job, and the first-release warning that historical
tag workflows become inert only after the post-OIDC token deletion.

- [ ] **Step 2: Assert precise diagnostic categories (2–5 min)**

Parse `release-summary` and require `if == "${{ always() }}"`, `needs` exactly
`[validate-release, wait-for-release-gates, build-package, promote-ghcr, publish-pypi]`, no OIDC/package permission,
and environment/dataflow for every upstream result plus its JSON output. Require rendered fields for mode, tag, full SHA,
main ancestry, observed package/citation/changelog versions and verdicts, structured preflight state/reason/candidates,
ten check verdicts and URLs, attempt/elapsed time, Docker run/evidence provenance, source
ref/digest/revision/version, alias attempted/write-succeeded/verified progress, all alias decisions
(`create`, `advance`, `no-op`, `skip-newer`, `skip-unorderable`,
`fail-conflict`), package filenames, PyPI result, and `dry run performed no production writes`.

- [ ] **Step 3: Run RED against stale token/tag guidance (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_workflow_contract.py -q`
Expected: failures cite AGENTS trap 12/“Never push” wording, absent repository-dispatch/default-branch guidance,
unresolved B4, stale workflow header, missing legacy-token migration warning, and incomplete summary categories.

- [ ] **Step 4: Update workflow header and always-run summary (2–5 min)**

Document the authenticated production `repository_dispatch` event type `vntyper_release`, `client_payload.tag`,
default-branch policy selection, dry-run syntax, immutable/floating aliases, missing-image rerun, canceled-pending
promotion retry, partial PyPI retry, and the prohibition on workflow-created tags. Add this job and implement its
Python renderer:

```yaml
  release-summary:
    needs: [validate-release, wait-for-release-gates, build-package, promote-ghcr, publish-pypi]
    if: ${{ always() }}
    runs-on: ubuntu-24.04
    permissions: {}
    steps:
      - env:
          VALIDATE_RESULT: ${{ needs.validate-release.result }}
          VALIDATE_JSON: ${{ needs.validate-release.outputs.summary_json }}
          GATES_RESULT: ${{ needs.wait-for-release-gates.result }}
          PREFLIGHT_JSON: ${{ needs.wait-for-release-gates.outputs.preflight_summary_json }}
          CHECKS_JSON: ${{ needs.wait-for-release-gates.outputs.check_summary_json }}
          EVIDENCE_JSON: ${{ needs.wait-for-release-gates.outputs.evidence_summary_json }}
          BUILD_RESULT: ${{ needs.build-package.result }}
          PACKAGE_JSON: ${{ needs.build-package.outputs.package_summary_json }}
          PROMOTE_RESULT: ${{ needs.promote-ghcr.result }}
          ALIASES_JSON: ${{ needs.promote-ghcr.outputs.alias_summary_json }}
          DRY_ALIASES_JSON: ${{ needs.wait-for-release-gates.outputs.dry_run_alias_summary_json }}
          PUBLISH_RESULT: ${{ needs.publish-pypi.result }}
          PUBLISH_JSON: ${{ needs.publish-pypi.outputs.publish_summary_json }}
          EVENT_NAME: ${{ github.event_name }}
        run: |
          set -euo pipefail
          python - <<'PY' >> "$GITHUB_STEP_SUMMARY"
          import json
          import os

          sections = (
              ("Validation", "VALIDATE_RESULT", "VALIDATE_JSON"),
              ("Docker evidence preflight", "GATES_RESULT", "PREFLIGHT_JSON"),
              ("Exact-SHA checks", "GATES_RESULT", "CHECKS_JSON"),
              ("Docker evidence", "GATES_RESULT", "EVIDENCE_JSON"),
              ("Package", "BUILD_RESULT", "PACKAGE_JSON"),
              ("GHCR aliases", "PROMOTE_RESULT", "ALIASES_JSON"),
              ("Dry-run aliases", "GATES_RESULT", "DRY_ALIASES_JSON"),
              ("PyPI", "PUBLISH_RESULT", "PUBLISH_JSON"),
          )
          print("# Release summary")
          if os.environ["EVENT_NAME"] == "workflow_dispatch":
              print("\n> Dry run performed no production writes.")
          for title, result_key, payload_key in sections:
              result = os.environ.get(result_key, "not reached") or "not reached"
              raw = os.environ.get(payload_key, "")
              print(f"\n## {title}\n\nJob result: `{result}`")
              if not raw:
                  print("\nStructured output: unavailable (job was skipped, failed before output, or was not reached).")
                  continue
              try:
                  payload = json.loads(raw)
              except json.JSONDecodeError as exc:
                  print(f"\nStructured output: invalid JSON ({exc.msg}); raw value intentionally not emitted.")
                  continue
              print("\n```json")
              print(json.dumps(payload, indent=2, sort_keys=True))
              print("```")
          PY
```

The renderer treats empty outputs as unavailable, includes `needs.*.result`, emits the dry-run no-write banner when
`EVENT_NAME == "workflow_dispatch"`, and never prints tokens, OIDC assertions, registry passwords, or signed URLs.

- [ ] **Step 5: Update AGENTS and ci-followups accurately (2–5 min)**

Replace “tag pushes publish immediately” with the existing-tag plus authenticated `vntyper_release` repository-dispatch
sequence; explain why the default-branch coordinator has no production `push.tags` trigger; expand version trap 12 to all three sources,
exact-SHA Docker evidence artifact, short-SHA/full label, aliases, anti-downgrade, prefix-collision, and rerun semantics;
state every main push is substantive while PRs retain filters. Reconcile the scripts type-check paragraph with the
already-landed milestone-6 quality-gate change rather than restoring stale “not type-checked” text. Mark B4 resolved
and leave secret deletion explicitly pending first live OIDC proof. State that historical tagged commits retain their
legacy tag-triggered token workflow, forbid new pre-milestone tags while the token exists, and call those workflows
inert only after the owner deletes the token following successful OIDC publication. Include a literal fenced workflow-state table in AGENTS
mapping validation/gates/build/promotion/publish/summary to permissions and retry behavior so the prose test has exact
content rather than substring-only coverage.

- [ ] **Step 6: Run GREEN and docs build (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_workflow_contract.py -q && make docs-build && make lint-actions`
Expected: guidance, diagnostics, strict docs, and workflow syntax pass.

- [ ] **Step 7: Commit contributor contract (2–5 min)**

```bash
git add AGENTS.md docs/development/ci-followups.md .github/workflows/publish-pypi.yml tests/unit/test_release_workflow_contract.py
git commit -m "docs(ci): record gated release recovery" -m "Refs #214" -m "Refs #218" -m "Refs #220"
```

Do not close #214, #218, or #220 from this implementation commit. Their specification requires links to the first
successful gated production release and, for #218, the separately authorized obsolete-secret deletion follow-up.

### Task 11: Local verification and handoff to master integration

**Files:**
- Verify only: all files above
- Do not modify: version sources, repository settings, tags, GitHub Releases, GHCR, PyPI, or Docker Hub

**Interfaces:**
- Consumes: complete branch implementation.
- Produces: local gate transcript and next-release owner checklist; no push, PR, workflow dispatch,
  repository dispatch, tag, or production release.

- [ ] **Step 1: Run focused tests freshly (2–5 min)**

Run: `pytest -m unit tests/unit/test_release_policy.py tests/unit/test_release_workflow_contract.py tests/unit/test_version_consistency.py tests/unit/test_workflow_consistency.py -q`
Expected: all focused tests pass with no skips caused by missing tracked workflows.

- [ ] **Step 2: Run formatting, lint, types, and unit gates (2–5 min command launch)**

Run: `make format-check && make lint && make type-check-all && make test-unit && make test-unit-cov && make patch-coverage`
Expected: every command exits 0; branch-inclusive floor remains at least 86 and patch coverage at least 80.

- [ ] **Step 3: Run workflow and docs gates (2–5 min command launch)**

Run: `make lint-actions && make docs-build && make ci-local`
Expected: actionlint, strict docs, every CI-equivalent job, and clean uv environment installation exit 0.

- [ ] **Step 4: Run the Docker-local equivalent (2–5 min command launch)**

Run: `make ci-local-docker`
Expected: image workflow-equivalent build/tests exit 0. If the daemon or required local base cannot run, capture the exact command/error in the PR; do not report the tier as passed and do not substitute a published mutable `latest` base after base inputs changed.

- [ ] **Step 5: Re-run invariant and secret scans (2–5 min)**

```bash
! rg -n 'PYPI_API_TOKEN|TWINE_PASSWORD|DOCKER_(USERNAME|PASSWORD)|@release/v1' .github/workflows/publish-pypi.yml
! rg -n '^\s*(docker (pull|run)|apptainer pull).*saei/vntyper' README.md docs/getting-started/installation.md docs/user-guide/docker.md
git diff origin/main...HEAD -- vntyper/version.py CITATION.cff docs/about/changelog.md
```

Expected: both negative scans are empty; version-source diff is empty.

- [ ] **Step 6: Review commit history and branch diff (2–5 min)**

Run: `git log --oneline origin/main..HEAD` and `git diff --check origin/main...HEAD`
Expected: focused Conventional Commits in task order and no whitespace errors, generated data, credential material, version bump, tag, or unrelated edit.

- [ ] **Step 7: Record the PR-check expectation without external mutation**

Write the ten exact expected check names into the SDD task report. Do not push, open/update a PR, wait on GitHub, push
a `v*.*.*` tag, manually dispatch a workflow, or send `repository_dispatch` here; master integration owns final reviews
before PR publication.

- [ ] **Step 8: Record the post-merge dry-run acceptance script without running it**

Record these owner commands/checks in the task report, explicitly marked not executed: manually dispatch a
pre-milestone strict tag to prove the ineligibility/no-write path; after the first post-milestone release, manually
dispatch that eligible existing tag to prove the full no-write summary. These are read-only `workflow_dispatch` runs,
not production repository dispatches. The checklist requires ten exact-SHA checks, Docker evidence
run/digest/revision/version, alias dry-run plan, package files, and “no production writes”.

- [ ] **Step 9: Hand off the first-live-release checklist without executing it**

Record: create the next version on `main`; wait for all exact-SHA checks; create/push its strict tag externally at the
post-milestone commit; send authenticated `repository_dispatch` with `event_type=vntyper_release` and that tag in
`client_payload.tag`; approve environment `pypi`; verify five GHCR aliases and PyPI attestations; only after that OIDC
publication succeeds, have an owner delete `PYPI_API_TOKEN`, record that legacy historical tag workflows are now inert,
and link evidence on #218. Do not create or push a tag at a pre-milestone commit while the token exists. No implementation
commit is created for verification-only evidence unless a failing gate requires a focused code/docs fix and a rerun.

## Completion definition

- Every traceability row has a passing named test and corresponding focused commit.
- Every `main` push produces non-skipped component and aggregator evidence plus a tested short-SHA image whose labels carry full SHA and package version.
- Production `vntyper_release` repository dispatch and manual existing-tag modes share validation/gates/build; only the
  exact production event/action can promote or request OIDC, and no production tag-push trigger exists.
- Missing/ineligible image evidence names the exact Docker run to rerun; partial alias/PyPI failures converge safely;
  floating aliases cannot downgrade; promotions cannot execute concurrently. A pending run canceled by GitHub's
  one-pending-slot concurrency semantics is rerun explicitly and has performed no writes.
- GHCR is the only active documented registry; all five supported aliases and rolling `main` semantics are truthful.
- PyPI build and publish privileges are separated; environment/action pin/skip-existing/token absence are executable contracts.
- Exactly nine generation-prose targets say `VNtyper 2`; versions, identifiers, historical prose, banners, and filenames remain protected.
- `make ci-local` and `make ci-local-docker` outcomes are recorded honestly, and no tag, release, registry, package, secret, or settings mutation occurs during implementation.
