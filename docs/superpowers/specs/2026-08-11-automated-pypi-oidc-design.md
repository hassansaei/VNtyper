# Automated PyPI OIDC Publishing Design

**Date:** 2026-08-11  
**Issue:** #236  
**Status:** Approved for implementation

## Problem

VNtyper's production controller already uses PyPI Trusted Publishing correctly: the
`publish-pypi` job keeps the `pypi` GitHub environment, grants only `id-token: write`,
and invokes the digest-pinned PyPA publish action without a username, password, or API
token. Release validation, exact-SHA checks, package construction, and GHCR promotion
all complete automatically.

The live `pypi` environment nevertheless has `required_reviewers = [hassansaei]` and no
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

The default-branch release controller reads those three public API resources during
`validate-release`, before package building or registry mutation, and invokes the pure
validator. API failures and policy mismatches fail with #236 and the owner command/UI
path. This turns future drift into an immediate diagnostic rather than a late waiting
deployment.

The workflow adds only the read permission required for this inspection. The
`publish-pypi` job remains the sole job with `id-token: write`; it retains no `contents`,
`packages`, or repository write permission.

### 3. Current v2.0.11 recovery

After the administrator saves the environment policy, inspect the existing waiting run.
If GitHub re-evaluates it, monitor it normally. If it remains pending, the administrator
uses **Start all waiting jobs** once, or the production `v2.0.11` dispatch is rerun. The
workflow already makes both paths safe: package upload uses `skip-existing`, and GHCR
promotion converges idempotently.

Success requires the workflow summary to report publication success and PyPI's JSON API
to expose version 2.0.11. If PyPI returns `invalid-publisher`, `@hassansaei` verifies the
four Trusted Publisher fields above; no token fallback is permitted.

### 4. v2.0.12 release completion

Once v2.0.11 is visible on PyPI and the policy PR is on `main`:

1. Create an annotated SSH-signed `v2.0.12` tag at
   `7f8583dd60565fec5d0297cbb26ccd2d7f439b22` and verify its signature and peeled commit.
2. Push that new tag without force.
3. Run the manual dry-run against the existing tag; require exact identity, all ten
   checks, Docker evidence, package build, and no production writes.
4. Send the authenticated `vntyper_release` repository dispatch.
5. Require terminal success for validation, gates, package build, GHCR promotion, PyPI
   publication, and release summary.
6. Verify PyPI 2.0.12 files, GHCR immutable and floating aliases, image revision/version
   labels, and the GitHub Release bound to the existing signed tag.
7. Remove only completed temporary release branches/worktrees.

After the first successful OIDC publication, the repository owner deletes the obsolete
`PYPI_API_TOKEN` secret separately, as the existing release policy requires.

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
acceptance then verifies the environment API, v2.0.11 recovery, signed-tag identity,
dry-run no-write evidence, production workflow, PyPI, GHCR, and GitHub Release.

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
