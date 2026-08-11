# Automated PyPI OIDC Publishing Design

**Date:** 2026-08-11
**Issue:** #236
**Status:** Implemented; live acceptance complete

## Problem

VNtyper's production controller already uses PyPI Trusted Publishing correctly: the
`publish-pypi` job keeps the `pypi` GitHub environment, grants only `id-token: write`,
and invokes the digest-pinned PyPA publish action without a username, password, or API
token. Release validation, exact-SHA checks, package construction, and GHCR promotion
all complete automatically.

At planning time, the live `pypi` environment had `required_reviewers = [hassansaei]` and no
deployment branch policy. Release v2.0.11 is therefore waiting after every automated
gate has passed. PyPI remains at 2.0.10. The active collaborator can push but has neither
`maintain` nor `admin` repository permission; GitHub requires repository Administration
write permission to change an environment or its deployment branch policies.

Removing only the reviewer would eliminate the delay but weaken the boundary: the
environment currently permits every ref, so a collaborator could run a modified branch
copy of the trusted workflow. The replacement must automate approval and narrow the
permitted ref at the same time.

## Goals

- Publish automatically after an authenticated production dispatch and the existing
  exact signed-tag, ancestry, metadata, ten-check, and immutable-image validations.
- Retain PyPI OIDC Trusted Publishing and the `pypi` environment claim.
- Remove every manual approval and wait timer from the publishing path.
- Permit the `pypi` environment only from `main`, which is the ref used by the
  default-branch `repository_dispatch` controller.
- Fail early with actionable diagnostics if the live environment drifts, instead of
  allowing a release to wait indefinitely after GHCR promotion.
- Complete v2.0.11 before publishing v2.0.12, then prove both releases end to end.

## Non-goals

- Do not reintroduce `PYPI_API_TOKEN` or any other long-lived package credential.
- Do not remove the `pypi` environment from the workflow or PyPI publisher identity.
- Do not weaken tag validation, main ancestry, required checks, Docker evidence,
  immutable alias behavior, concurrency, or idempotent retry.
- Do not automate creation of arbitrary release tags. The explicitly authorized
  v2.0.12 tag remains signed, annotated, immutable, and pinned to the already-green
  release commit.
- Do not add a GitHub App merely to auto-approve a gate that no longer expresses the
  desired policy.

## Evidence and platform constraints

PyPI documents that a GitHub environment is optional but strongly recommended for a
Trusted Publisher. When present, its name is part of the publisher configuration and the
OIDC identity. GitHub documents the corresponding OIDC subject as
`repo:OWNER/REPOSITORY:environment:ENVIRONMENT`. Therefore the environment stays named
`pypi`; only its GitHub-side protection rules change.

GitHub environment protection rules are evaluated before a job reaches a runner. A
required reviewer necessarily introduces human latency. GitHub's environment and branch
policy APIs require repository Administration write permission for mutation, while reads
require only Actions read access (and public repository state is readable without
authentication). The current account cannot perform the one-time mutation, so #236
contains exact instructions for `@hassansaei`.

Official references:

- <https://docs.pypi.org/trusted-publishers/using-a-publisher/>
- <https://docs.pypi.org/trusted-publishers/adding-a-publisher/>
- <https://docs.pypi.org/trusted-publishers/troubleshooting/>
- <https://docs.github.com/en/actions/reference/security/oidc>
- <https://docs.github.com/en/actions/reference/workflows-and-actions/deployments-and-environments>
- <https://docs.github.com/en/rest/deployments/environments>
- <https://docs.github.com/en/rest/deployments/branch-policies>

## Selected design

### 1. One-time live environment migration

A repository administrator updates `pypi` atomically to:

- no required reviewers;
- zero wait timer;
- custom deployment branch policies enabled;
- exactly one branch policy: branch `main`;
- no environment secrets.

The environment itself is not deleted or renamed. The PyPI Trusted Publisher continues
to identify owner `hassansaei`, repository `VNtyper`, workflow
`publish-pypi.yml`, and environment `pypi`.

### 2. Version-controlled fail-fast contract

Add a small pure Python contract module under `scripts/` that accepts three GitHub REST
responses: the environment object, its deployment branch-policy collection, and its
custom deployment-protection-rule collection. It returns
success only when the live shape is exactly the selected policy. It rejects reviewers,
wait timers, custom protection rules, an unrestricted environment, missing/extra policy
patterns, or a tag policy masquerading as `main`. The built-in `branch_policy` protection
rule is required; it is the API representation of the selected main-only branch boundary,
not a manual approval rule.

The default-branch release controller preflights exactly the environment,
deployment-branch-policy, and custom-deployment-protection-rule API responses during
`validate-release`, before package building or registry mutation, and invokes the pure validator.
API failures and policy mismatches fail with #236 and the owner command/UI path. This turns future
drift into an immediate diagnostic rather than a late waiting deployment. Zero environment secrets
is separately verified live administrator state, not controller-enforced.

The workflow adds only the read permission required for this inspection. The
`publish-pypi` job remains the sole job with `id-token: write`; it retains no `contents`,
`packages`, or repository write permission.

Maintainer guidance records the same contract: the environment is reviewer-free, has no wait
timer, uses custom branch policies for the exact branch `main`, and has no custom
deployment-protection rules. A mismatch cites #236 and fails before package or registry writes.
OIDC remains the only publisher; never reintroduce `PYPI_API_TOKEN`. Zero environment secrets is
separately verified live administrator state, not controller-enforced.

### 3. Completed v2.0.11 recovery

The production OIDC run `31465885545` completed green and PyPI exposed 2.0.11. Its
`skip-existing` package upload and idempotent GHCR promotion behaved as designed; no token
fallback was used.

### 4. Completed v2.0.12 release

The production OIDC run `31464328451` completed green. PyPI exposes 2.0.12, and the
repository owner deleted `PYPI_API_TOKEN` after these green OIDC releases. The configured
publisher remains owner `hassansaei`, repository `VNtyper`, workflow `publish-pypi.yml`, and
environment `pypi`.

## Error handling

- Live-policy API error: fail before any production write and identify the failed URL.
- Live-policy mismatch: fail before any production write and print the exact observed
  reviewers/timer/branch policies plus #236.
- Existing v2.0.11 run remains waiting: use the documented one-time start or idempotent
  redispatch; do not alter package or alias state manually.
- OIDC `invalid-publisher`: verify PyPI's owner/repository/workflow/environment tuple;
  never fall back to the legacy token.
- Existing or conflicting v2.0.12 tag: fail closed; never move or force-push it.
- Partial production release: rerun the same tag; existing skip/convergence rules remain
  authoritative.

## Tests and verification

TDD covers the pure live-policy validator with literal fixtures and mutations for every
forbidden rule, the required built-in branch-policy rule, every enabled custom deployment
protection rule, and every missing/extra/wrong-type `main` policy. Workflow contract tests
prove the validation runs before build/promotion, API failures propagate, permissions
remain read-only except for the isolated OIDC job, and neither token credentials nor
approval bypasses appear.

Repository gates are `make check-all`, `make ci-local`, workflow semantic tests,
actionlint, and the release workflow's existing mutation-sensitive safety suite. Live
acceptance verified the environment API, green OIDC runs `31465885545` (2.0.11) and
`31464328451` (2.0.12), PyPI 2.0.12, GHCR, GitHub Release, and deletion of
`PYPI_API_TOKEN`.

## Rejected alternatives

- **GitHub App auto-approver:** adds an installation, credential, and custom protection
  service solely to approve every release. It preserves ceremony without meaningful
  review and has a larger failure surface.
- **Repository token publishing:** reintroduces a long-lived credential and bypasses
  OIDC's short-lived, workflow-bound identity.
- **Delete the environment:** breaks the currently configured PyPI OIDC subject and
  removes the useful environment identity boundary.
- **Remove reviewers but allow every ref:** automates publication at the cost of allowing
  modified branch workflow copies to reach the trusted environment.
