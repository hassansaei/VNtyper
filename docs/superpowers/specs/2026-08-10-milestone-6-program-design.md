# Milestone 6 Gates, Harness, and Release Automation Program Design

**Date:** 2026-08-10
**Milestone:** [#6 — 5. Gates, harness and release automation (parallel)](https://github.com/hassansaei/VNtyper/milestone/6)
**Issues:** [#71](https://github.com/hassansaei/VNtyper/issues/71), [#204](https://github.com/hassansaei/VNtyper/issues/204), [#208](https://github.com/hassansaei/VNtyper/issues/208), [#211](https://github.com/hassansaei/VNtyper/issues/211), [#214](https://github.com/hassansaei/VNtyper/issues/214), [#218](https://github.com/hassansaei/VNtyper/issues/218), [#219](https://github.com/hassansaei/VNtyper/issues/219), [#220](https://github.com/hassansaei/VNtyper/issues/220), [#226](https://github.com/hassansaei/VNtyper/issues/226)

## Decision summary

Milestone 6 will land as one dependency-ordered program with five bounded tracks:

1. **Quality gates** — type-check and branch-cover `scripts/` without lowering any floor (#204, #211).
2. **Harness matrix** — make artifact-policy assertions meaningful, fix option precedence, add an independent non-SHARK FASTQ path, and close already-completed CRAM fixture work with evidence (#71, #226).
3. **Mutation safety** — execute mutants only in a disposable detached worktree and retain reports outside it (#208).
4. **Release and naming** — promote the tested image by digest, gate releases on the exact commit, publish PyPI with OIDC, align install surfaces, and distinguish the generation name from versions (#214, #218, #220).
5. **Exception policy** — establish an exact, classified BLE001 baseline and protect the stable fail-open decisions that can affect results (#219).

The selected implementation is the smallest coherent **risk-first, gate-completing** approach: type-check scripts first, remove in-place mutation risk second, implement harness and release contracts, classify the resulting final exception surface, then add sufficient decision-sensitive script tests before extending the permanent coverage source. This orders the data-loss-prevention fix before the substantial coverage expansion while still making the final tree enforce every gate.

## Source of truth and reconciled acceptance criteria

The governing evidence order is: issue acceptance criteria and live issue discussion; `AGENTS.md`; executable tests and workflows; current production behavior; historical plans. Where these disagree, this section makes the contradiction explicit.

| Issue | Exact acceptance criteria from the live issue | Reconciled obligation on current `main` |
| --- | --- | --- |
| #71 | Docker pipeline success; API success/error/JSON coverage; adVNTR CLI output; one passing scenario per `pipeline`, `report`, `cohort`, and `online`; flag combinations assert files present or absent; an alternate FASTQ runs without SHARK. | Criteria 1–4 already have executable coverage and are retained. Implement criteria 5–6. All nine normal BAM matrix cases now intentionally exit 1, so they cannot prove cleanup; use successful cases. The current second FASTQ is another view of sample 6449, not an alternate specimen. |
| #204 | Annotate `dir_counts`; include `scripts/` in `type-check`; document that Ruff and mypy production scopes match. | Implement verbatim, keep test mypy as a separate invocation, and update stale repository guidance. |
| #208 | Stop interruption from leaving a live mutant in production source; disposable copy/worktree is the stated fix direction. | A disposable detached Git worktree overlaid with the current tracked and non-ignored untracked working state is mandatory. Signal handlers alone are insufficient because SIGKILL, interpreter failure, and host loss cannot unwind `finally`; HEAD-only execution would silently exclude a newly written mutation-killing test. |
| #211 | Test `scripts/` and put it under coverage after a deliberate floor decision. | The issue title's 80% is stale: scripts at 80% would reduce combined branch coverage below the 86% hard floor. Require scripts at least 88.00% and combined coverage at least 86%, then add `scripts` to coverage source without changing floor 86, target 86, patch floor 80, or `branch = true`. |
| #214 | Publish working `latest` and semver GHCR aliases; decide Docker Hub; reconcile install docs; add a drift guard. | The metadata tag rule is currently unreachable because the workflow does not run on tags. Promote the exact already-tested `sha-<commit>` digest after exact-SHA gates, make GHCR authoritative, and remove active Docker Hub instructions. |
| #218 | Configure the PyPI publisher/environment; grant OIDC; eliminate token use; prove a successful OIDC release before secret deletion. | PyPI publisher and GitHub `pypi` environment are live, so the blocker and label are stale. Split unprivileged build from the OIDC-only publish job. The first successful real release is the only production proof; secret deletion remains a post-release owner action and is not performed by this PR. |
| #219 | Correct the stale count, make it self-checking, classify handlers, prioritize fail-open category, and do not globally narrow all catches. | Current Ruff BLE001 count is 103, not 67 or 46. Guard the exact aggregate and a targeted stable-symbol fail-open audit. Do not create a 100-file waiver churn or change clinical-adjacent behavior without an existing contract and behavior test. |
| #220 | Change exactly nine generation-name `VNtyper 2.0` lines to `VNtyper 2`; repair adjacent README grammar; preserve all real versions and machine-readable names. | Implement exactly those presentation changes. Current version is 2.0.10 and six legitimate historical `VNtyper 2.0.x` references remain, so the issue's old one-hit verification command is not authoritative. Bare `VNtyper` remains canonical for identifiers, citation, API, and machine-readable banners. |
| #226 | Implement the cited `build_reference_dependent_fixture` or replace the dangling citation. | Already implemented and tested on `main` by commits `80f42fbba63c5c045ee605e49949cdd970884320` and `edabd3eaf594b785906d8ba03c9cce60f1c6babd`. No duplicate code is allowed; close with traceable code and test evidence. |

## Scope

### In scope

- Production mypy coverage for `vntyper/`, `docker/app/`, and `scripts/` on Python 3.10–3.13.
- Branch-inclusive unit coverage for `scripts/`, with meaningful assertions for every new or modified decision and failure path.
- Disposable-worktree mutation execution, selected-target/output dirty guards, working-state overlays, isolation diagnostics, deterministic cleanup, and report round-trip stability.
- Data-driven positive and negative artifact expectations for successful pipeline cases, consistent forwarding of case options, delete-over-keep precedence, and a distinct non-SHARK FASTQ specimen.
- Release check-run coordination for an exact commit, digest promotion, semver alias and anti-downgrade rules, PyPI OIDC privilege separation, GHCR-only current install documentation, and the nine naming edits.
- A self-checking BLE001 policy and classification focused on fail-open behavior.
- Specifications, executable plans, traceability, verification evidence, adversarial reviews, and a pull request covering all nine issues.

### Explicitly out of scope

- Changing genotyping, scoring, confidence, flag definitions, or config-driven clinical wording.
- Upgrading Kestrel, external tools, package dependencies, the Python floor, or image environments.
- Raising or lowering coverage floors solely to make this milestone pass.
- Globally enabling BLE001, mechanically narrowing every broad handler, or changing G004 policy.
- Publishing to Docker Hub or reintroducing Docker Hub credentials.
- Pushing a production release tag, publishing a new package version, deleting `PYPI_API_TOKEN`, or merging the pull request.
- Rebuilding #226, broad CLI/API redesign, test archive replacement, or unrelated cleanup.

## Architecture and component boundaries

```text
current main commit
  |
  +--> static and unit gates --------> exact-SHA CI check-run
  |
  +--> Docker build/tests -----------> GHCR sha-<commit> digest
  |                                      |
  |                                      +--> release coordinator
  |                                             - validate tag/version/main ancestry
  |                                             - poll required exact-SHA check-runs
  |                                             - promote digest and aliases
  |                                             - invoke OIDC-only PyPI publish
  |
  +--> disposable mutation worktree
  |       - mutable source and pytest cwd
  |       - known-killed canary
  |       - no report output
  |     real checkout
  |       - preserved working source and guarded mutation target
  |       - generated reports only
  |
  +--> test-data case declarations --> shared artifact assertion policy
  |
  +--> exact BLE001 aggregate + classified stable-symbol audit
```

The five track specifications own detailed interfaces. Shared integration contracts are frozen here:

- `make type-check` checks production code in `vntyper/ docker/app/ scripts/`; `make type-check-all` adds tests separately.
- `[tool.coverage.run].source` ends as `vntyper`, `docker/app`, and `scripts`; branch measurement and all three numeric thresholds remain unchanged.
- Mutation code receives distinct immutable `REAL_REPO_ROOT`, disposable `sweep_root`, and pre-resolved persistent output paths. Subprocess `cwd`, not only `PYTHONPATH`, selects the disposable source.
- Test cases may declare literal artifact paths expected present or absent. Assertions run only after the case reaches its declared successful exit.
- `--delete-intermediates` wins when both cleanup flags are present, matching CLI help; this behavior change is recorded under Unreleased.
- A release never builds tag-context source. Every push to `main` runs the substantive CI and Docker jobs even for
  documentation-only diffs and uploads a run-attempt-scoped full-SHA/post-push-registry-digest evidence artifact. The release selects the
  successful exact-SHA Docker run, verifies that evidence plus the short tag/full revision/package labels, and promotes
  the immutable evidence digest only after every named component and aggregate check for that exact SHA succeeds.
- Production release aliases are `vX.Y.Z`, `X.Y.Z`, `X.Y`, `X`, and `latest`; ordinary `main` images remain `main` and SHA-tagged.
- The existing `publish-pypi.yml` filename and `pypi` environment name remain fixed because PyPI's trust configuration names them.

## Dependency ordering and concurrency

1. Commit the program and five track specifications after independent adversarial review.
2. Commit master and track plans plus requirements-to-task traceability after independent review.
3. Land #204 so every later script edit is type-checked.
4. Land #208 to remove source-corruption risk before mutation work is used for later tasks.
5. Land #71 harness changes; reconcile #226 using existing evidence.
6. Land #214/#218/#220 as a single release-surface track because they share workflows, docs, and image metadata.
7. Add #211 tests, remeasure, and only then include `scripts/` in the permanent source set.
8. Land #219 after the large script/test tranche stabilizes so its baseline is not instantly obsolete.
9. Integrate, run fresh final gates, conduct whole-branch reviews, push, open the PR, and monitor required checks.

Read-only research and review may run concurrently. The installed Superpowers 6.2.0
`subagent-driven-development` skill forbids concurrent implementation subagents, so all implementation tasks are
serialized even when their file sets are disjoint. `Makefile`, `pyproject.toml`, workflows, shared test helpers, and
integration contracts remain explicit serialization points if a future skill version permits implementation parallelism.

## Interfaces, inputs, outputs, and invariants

### Inputs

- A clean branch based on the verified current `origin/main` fork point.
- Live issue bodies, comments, milestone state, repository policies, and current workflow/check-run names.
- Existing 1.1 GB integration archive and Docker daemon when available.
- A semantic tag in strict `vMAJOR.MINOR.PATCH` form for release validation; workflow dispatch accepts an existing tag only in dry-run mode.

### Outputs

- Six reviewed specifications, six reviewed plans, and a requirements-to-task matrix under `docs/superpowers/`.
- Focused Conventional Commits and an auditable progress ledger.
- Updated gates, tests, harness, release workflows/docs, exception policy, and generated mutation evidence where the generator owns it.
- A pushed branch and pull request against `main`; no merge and no release tag.

### Program invariants

- The real checkout's production source is never overwritten by a mutant.
- A failing unmutated baseline can never be reported as a killed mutant.
- A test result is accepted only when imports and subprocess cwd resolve inside the disposable tree.
- A coverage increase never disables branch coverage or reduces floor 86, target 86, or patch target 80.
- New unit files carry `pytestmark = pytest.mark.unit` and run from repository root without network, Docker, or reference data.
- Code and workflow helpers remain Python 3.10 compatible and shell commands are quoted for paths with spaces.
- No release alias or PyPI artifact is published until every substantive component and aggregate
  check for the exact tag commit succeeds; a skipped component is not success evidence.
- `latest` and floating semver aliases never move to a lower release; an older rerun records a skip rather than failing or downgrading.
- Rerunning a partially completed release converges on the same immutable exact-version digest.
- The PyPI publishing job has `id-token: write` and no source checkout, build toolchain, repository token, password, or long-lived package credential.
- Product-generation edits never mutate `2.0.x` versions, package/import/CLI names, citation identity, or parsed banners.
- Fail-open behavior is explicit, classified, and protected by observable behavior tests.

## Success, partial failure, recovery, retry, and idempotency

| Operation | Success | Partial failure and recovery | Retry/idempotency |
| --- | --- | --- | --- |
| Type/coverage gates | Every configured path is checked and combined branch coverage is at least 86%. | A named file/error or precise coverage measurement fails the gate; add tests or annotations, never weaken configuration. | Deterministic on the same tree and supported Python version. |
| Mutation sweep | Baseline passes; canary is killed; requested mutants run in disposable source; reports are rendered to the real output root. | Failure before mutation leaves real source untouched. Interruption may leave only an expendable worktree/admin entry; log its path and provide cleanup. Report generation uses completed records only. | A rerun recreates a disposable snapshot of the then-current tracked and non-ignored working state over the same HEAD and can replace deterministic generated reports. |
| Harness case | Declared successful exit plus literal presence/absence assertions all pass. | Exit mismatch stops artifact assertions with captured stdout/stderr. Artifact mismatch names each missing or unexpected path; case output remains for diagnosis. | Case output isolation and cleanup rules make reruns independent. Archives and cleanup assertions cannot be satisfied by stale output. |
| Release coordination | Exact SHA is on main, version matches tag, required checks succeed, SHA image is present, aliases converge, then PyPI publishes. | Polls are bounded. Failure names missing/pending/failed checks or image. Exact alias can be reused only for the same revision. Floating-alias downgrade is skipped with notice. PyPI failure does not change the image digest. | Exact alias is immutable by policy; repeated promotion is a no-op/convergence. PyPI keeps existing skip-existing behavior for partial reruns. |
| Exception audit | Current aggregate and every targeted stable symbol match an allowed category and behavior. | New, moved, or reclassified handlers fail with precise path/symbol/category diagnostics. | Deterministic static audit plus behavior tests. |

## Compatibility and environment constraints

- Runtime and helper code supports Python 3.10 through 3.13. Do not use Python 3.11-only syntax or stdlib APIs.
- CI installation continues to use `uv venv`, `VIRTUAL_ENV`, and `.venv/bin` rather than `uv pip install --system`.
- Integration commands run from the repository root. Missing archive dependencies must be reported, with `VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1` used to prevent an unintended download.
- Existing conda environments, Java 11, Kestrel 1.0.1, Docker base-hash binding, shell-based process substitution, and relative configuration paths remain unchanged.
- Worktree creation requires Git and a resolvable committed HEAD. Modified/untracked working state is copied through a validated manifest; only selected mutation targets and persistent output paths are rejected when dirty.
- GitHub workflow shell is Bash on Ubuntu; commands quote repository, tag, digest, and artifact paths and use bounded retries.

## Security and GitHub permission model

- Default workflow permissions remain empty or read-only; permissions are granted per job.
- Pull requests never publish application images or packages through the milestone release paths. The preserved
  same-repository missing-base bootstrap may publish its content-hash-addressed base image; fork PR behavior remains
  the explicit read-only-token failure documented by the repository. No `pull_request_target` is introduced.
- Docker `build-and-test` retains its existing job-level `packages: write` because the same job publishes tested
  application tags on main pushes. PR execution contains no application-image registry-write step: cache export is
  disabled and every application push/evidence step is guarded by exact push-to-main authority. Release inspection
  uses read-only metadata access, and promotion receives package write only for its job.
- The PyPI build job has `contents: read`, uses `persist-credentials: false`, and uploads a workflow artifact. The publish job has only `id-token: write`, uses environment `pypi`, and invokes the SHA-pinned PyPA publisher without secrets.
- Workflow dispatch is a dry run over an already-existing tag. It must not create or move a production tag or publish aliases/packages.
- Release validation prevents an arbitrary tag/ref from publishing by verifying strict format, package version, main ancestry, exact-SHA check-runs, and source image revision.
- Mutation guards resolve and compare canonical paths before any write, rejecting the actual repository root and any source path outside the disposable worktree.

## Observability and diagnostics

- Gates print checked source scopes, Python version, test counts, exact branch-inclusive percentages, thresholds, and actionable failures.
- Mutation runs log immutable HEAD, real root, disposable root, output root, baseline/canary status, per-mutant outcome, cleanup result, and any abandoned path.
- Harness failures print case ID, command, declared exit, actual exit, stdout/stderr location, and missing/unexpected artifacts.
- Release jobs summarize tag, commit, main-ancestry result, each required check-run and conclusion, source image reference/digest/revision/version, aliases updated/skipped, anti-downgrade comparisons, and PyPI result. Poll timeouts include the last observed state.
- BLE001 audit failures name the exact path, stable symbol, observed category/count, and required remediation.
- Final PR evidence records commands, exit codes, test counts, coverage figures, skipped tiers, check URLs/statuses, and adversarial-review dispositions.

## Test strategy and objective acceptance checks

Every implementation task follows red → observed intended failure → minimal green → refactor while green. Assertions must detect wrong values, wrong paths, wrong branches, missing side effects, unsafe writes, and partial failures rather than merely successful calls.

### Focused checks

- #204: mypy on `vntyper/ docker/app/ scripts/` under the Python matrix; the previously failing annotation is accepted; tests remain a separate invocation.
- #208: unit tests for path guards, worktree lifecycle, subprocess cwd/env, cleanup failures, baseline failure, canary survival, interruption, and output-root separation; a scoped real sweep proves the real source is byte-identical before/after.
- #71: pure unit tests for artifact-policy parsing/enforcement and option forwarding; successful integration cases prove delete/keep/archive/no-archive outcomes; a b178 paired FASTQ case succeeds without SHARK. Existing Docker, API, adVNTR, and subcommand tests are cited and rerun.
- #214/#218/#220: executable workflow-contract helpers/tests cover tag parsing, exact-SHA gate decisions, bounded polling states, alias planning, anti-downgrade, digest/revision mismatch, permissions, artifact flow, GHCR documentation consistency, and naming preservation; actionlint validates workflow semantics.
- #219: exact aggregate check plus stable-symbol classification and fail-open/fail-closed observable behavior tests.
- #211: decision-sensitive unit tests cover both zero-coverage modules and the lowest coverage harness modules until scripts-only branch coverage is approximately 88% and final combined branch coverage is at least 86%; the source-membership guard fails if `scripts` is removed.
- #226: rerun existing unit and reference-contract integration tests; cite implementation commits and symbols.

### Final commands

Run freshly from repository root on the final tree:

```bash
make format-check
make lint
make type-check
make test-unit
make test-unit-cov
make patch-coverage
make check-all
make ci-local
make ci-local-docker
```

Run relevant mutation checks and the generator round-trip. Run `make check-full` only if the archive and external dependencies are available; otherwise run the available integration targets with fail-fast download suppression and report the unavailable tier exactly. Required GitHub PR checks must reach successful terminal states before completion is claimed.

## Rollout and release behavior

The PR changes automation but does not itself release. After merge, the next semantic release must be created from a commit already on `main`. The release workflow validates and waits for that exact commit's required CI and Docker results, promotes the tested image digest, then enters the protected `pypi` environment. Maintainer approval remains the irreversible publication gate. After the first successful OIDC release, an owner may remove the obsolete `PYPI_API_TOKEN`; that settings change is deliberately outside this PR.

Existing Docker Hub images remain historical and unsupported. Current install surfaces use GHCR. Existing exact release aliases are never repointed to different revisions; floating aliases advance monotonically.

## Objective definition of done

- All six specifications, the master plus five bounded implementation plans, and the separate execution prompt exist,
  have no placeholders, and have passed self-review plus independent `claude-opus-5` maximum-effort review with no
  unresolved Critical or High finding.
- The persistent task ledger maps every criterion above to implementation, tests, commits, and review verdicts.
- All nine issues are traced to current or new executable evidence; #226 is not redundantly reimplemented.
- Fresh local gates applicable to every changed file pass with captured results, and unavailable tiers are not claimed.
- Final whole-branch Superpowers and independent Opus reviews have no unresolved Critical or High finding.
- The branch is clean, pushed without force, a PR against `main` exists, and required GitHub checks are terminal and successful.
- The PR is not merged and no production version tag is pushed by this work.

## Unresolved questions

None. Live issue discussion resolves #218; repository behavior plus the user's authorized
milestone direction select the conservative #214/#219/#220 policies where those issues
have no maintainer comments. Bare canonical `VNtyper` names remain unchanged because
that is the explicit smallest-scope recommendation in #220, not because discussion
closed its `needs-decision` label.
