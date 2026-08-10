# Milestone 6 Release Automation and Naming Design

**Date:** 2026-08-10
**Status:** Approved design; implementation in progress through Task 5
**Milestone:** [#6 — 5. Gates, harness and release automation (parallel)](https://github.com/hassansaei/VNtyper/milestone/6)
**Issues:** [#214](https://github.com/hassansaei/VNtyper/issues/214),
[#218](https://github.com/hassansaei/VNtyper/issues/218),
[#220](https://github.com/hassansaei/VNtyper/issues/220)
**Evidence baseline:** `main` at `ebb15b26631242a3295607e4eda4c68f688cd9a2`

## 1. Outcome

A production release is a promotion of artifacts already tested for one exact commit. An
authenticated `repository_dispatch` with event type `vntyper_release` names a pre-existing strict
`vMAJOR.MINOR.PATCH` tag in `client_payload.tag`. GitHub loads the coordinator from the default
branch; the coordinator peels that tag to a commit already on `main`, waits for that commit's
`CI Success` and `Docker Success` check-runs, verifies the tag and package version, and promotes
the existing GHCR SHA image by digest. It does not rebuild source in tag context.

The promoted digest receives these aliases:

- `vMAJOR.MINOR.PATCH`
- `MAJOR.MINOR.PATCH`
- `MAJOR.MINOR`
- `MAJOR`
- `latest`

The first two are immutable release identities. The last three are floating aliases that move
only forward. `main` remains the rolling image built from `main`; a release never moves it.
GHCR is the sole current container registry. Existing Docker Hub images are legacy and all active
Docker Hub install instructions are removed.

The same workflow builds the Python distribution without publishing privilege and publishes it
only after the exact-SHA gates and GHCR promotion succeed. The publish job uses PyPI Trusted
Publishing, environment `pypi`, a SHA-pinned PyPA action, and no long-lived package credential.

Presentation prose calls the product generation **VNtyper 2**. Exactly the nine generation-name
targets enumerated by #220 change from `VNtyper 2.0` to `VNtyper 2`; package versions, historical
statements, dependencies, identifiers, filenames, parsed banners, and canonical bare `VNtyper`
names do not change.

## 2. Source of truth and exact acceptance criteria

The live issue body and discussion are the requirements source. Current code and registry state
are evidence used to make those requirements executable. Historical line numbers and counts are
not requirements when the underlying files have moved.

### 2.1 Issue #214 — GHCR aliases and registry truth

[The issue](https://github.com/hassansaei/VNtyper/issues/214) requires all of the following:

1. Working GHCR `latest` and semantic-version aliases, not documentation for tags CI never
   publishes.
2. An explicit Docker Hub decision: automate it, or remove it from current instructions and mark
   it legacy.
3. Install documentation reconciled with the registry that is actually published.
4. An executable drift guard so documentation and publication rules cannot silently diverge
   again.
5. A Docker release guarantee equivalent to the package tag/version check.

Acceptance for this design is stronger than adding metadata-action patterns: the workflow must
actually be invoked for a release and must alias the digest already tested for the tag commit.
This resolves the current contradiction in `.github/workflows/docker-build.yml`: `type=ref,event=tag`
exists under `steps.meta`, but the workflow trigger admits only pushes to `main`, pull requests,
schedules, and manual runs, so the tag rule is unreachable.

### 2.2 Issue #218 — PyPI Trusted Publishing

[The issue](https://github.com/hassansaei/VNtyper/issues/218) and the
[maintainer completion comment](https://github.com/hassansaei/VNtyper/issues/218#issuecomment-5215430366)
require this ordered migration:

1. PyPI trusts owner `hassansaei`, repository `VNtyper`, workflow `publish-pypi.yml`, and
   environment `pypi`.
2. GitHub has the matching `pypi` environment.
3. The publish job declares that environment and `id-token: write`.
4. PyPI upload uses `pypa/gh-action-pypi-publish` without `TWINE_PASSWORD` or a `password` input.
5. The long-lived `PYPI_API_TOKEN` is deleted only after a successful real OIDC release.

Items 1 and 2 are complete according to the maintainer's comment and the live GitHub environment.
The implementation owns items 3 and 4. Item 5 remains an owner-only rollout action after the first
production proof; neither this design nor its implementation deletes repository settings.

### 2.3 Issue #220 — generation prose, not versions

[The issue](https://github.com/hassansaei/VNtyper/issues/220) requires:

1. Exactly eight generation-name edits in `README.md` and one in `docker/Dockerfile`.
2. `VNtyper 2.0` becomes `VNtyper 2` only at those nine presentation targets.
3. The adjacent README grammar error is repaired while retaining the historical `VNtyper v1`
   reference.
4. Every real `2.0.x` version, version heading, dependency number, fixture, image/tag example,
   machine-readable name, and generation-2 filename is preserved.
5. Bare canonical names stay bare: citation title, distribution and command names, docs site
   identity, OCI title, CLI help, API names, and parsed `VNtyper Version:` banners remain
   `VNtyper` rather than gaining `2`.

The issue's old verification prediction of one remaining `VNtyper 2.0.x` line is stale. At this
baseline there are multiple legitimate post-2.0.7 historical references. Acceptance therefore
checks the nine explicit targets and protected version/name invariants rather than a repository-wide
one-line count.

## 3. Scope

### 3.1 In scope

- Change `.github/workflows/docker-build.yml` so every main push runs the substantive image job
  and its existing short-SHA tag carries an unambiguous full-revision and package-version label.
- Refactor `.github/workflows/publish-pypi.yml` into the release coordinator while preserving its
  filename because PyPI's trusted-publisher identity names it exactly.
- Trigger production only through authenticated `repository_dispatch` event type
  `vntyper_release`, with the candidate supplied as `client_payload.tag`, so the coordinator always
  runs from the default branch rather than from an arbitrary tagged historical commit.
- Add deterministic, pure release-policy logic and unit tests for tag parsing, check-run verdicts,
  aliases, immutability, and anti-downgrade decisions.
- Poll the exact commit's substantive CI/Docker component jobs plus `CI Success` and
  `Docker Success` with a finite deadline; a skipped component is not release evidence.
- Prove tag, commit, `main` ancestry, package version, source image revision, source image version,
  and source digest before any production mutation.
- Promote the already-tested SHA image by digest to all five production aliases.
- Split unprivileged Python build/check from OIDC-only PyPI publication and retain partial-rerun
  `skip-existing` behavior.
- Make `workflow_dispatch` an inspection/build dry run that cannot publish or create/move a tag.
- Make GHCR authoritative; remove Docker Hub from current install commands and describe its
  existing artifacts as unsupported legacy where a transition note is useful.
- Apply the nine explicit product-generation prose edits and the adjacent README grammar repair.
- Add workflow, registry-documentation, OCI-label, permission, alias, and naming invariant tests.
- Update the obsolete Trusted Publishing block in `docs/development/ci-followups.md`, the release
  workflow header, and the affected release/type-check guidance in `AGENTS.md` so contributor
  instructions describe the post-migration checks and permissions accurately.

### 3.2 Out of scope

- Pushing a production tag, creating a GitHub Release, publishing a new version, or running the
  first live release as part of the implementation PR.
- Publishing, mirroring, or restoring automation for Docker Hub.
- Deleting Docker Hub images, Docker Hub credentials, `PYPI_API_TOKEN`, or any repository setting.
- Rebuilding a release image from a tag or copying a locally built image into production.
- Changing the `main` alias into `latest`, or publishing `latest` from ordinary main builds.
- Supporting prerelease, build-metadata, branch, date, or non-semantic production tags.
- Changing `vntyper/version.py`, `CITATION.cff`, changelog release headings, dependency versions,
  fixtures, or the current package version.
- Renaming the `vntyper` distribution, import, executable, repository, API, citation, docs site,
  parsed summary banner, or `snakemake/vntyper2*` files.
- Redesigning Docker base-image publication, PR/fork image behavior, package contents, or runtime
  genotyping behavior.

## 4. Current evidence and contradictions

At the baseline commit:

- `.github/workflows/docker-build.yml` publishes only after `Build and test image` succeeds. Its
  terminal required job is `Docker Success`; `CI Tests` analogously exposes `CI Success`.
- Main publishes `main` and a metadata-action SHA tag. GHCR has matching `main` and
  `sha-ebb15b2` manifests, but no `latest`, `v2.0.10`, or `2.0.10` manifest.
- Docker Hub's `saei/vntyper:latest` is an older manually maintained artifact; Docker Hub
  automation was removed by commit `1ea4ddc`.
- `.github/workflows/publish-pypi.yml` is one combined job. It validates only tag versus
  `vntyper/version.py`, then uploads with `TWINE_PASSWORD=${{ secrets.PYPI_API_TOKEN }}` and
  `twine upload --skip-existing`.
- A tag-push event selects the workflow from the tagged commit. A coordinator fixed only on the
  current default branch therefore cannot guard a new tag that points at a pre-milestone commit;
  that historical workflow can still use `PYPI_API_TOKEN` until the secret is deleted after the
  first successful OIDC publication.
- Release `v2.0.10` was published to PyPI while the exact main CI and Docker workflows were still
  running. Tag/version agreement alone therefore does not prove the released commit passed.
- `tests/unit/test_version_consistency.py` protects metadata and environment version invariants,
  while `tests/unit/test_workflow_consistency.py` protects only base-image hash agreement. No
  executable test currently owns release aliases, OIDC permissions, Docker Hub drift, or product
  generation naming.
- `docker/metadata-action` labels override same-key Dockerfile labels. The main build must pass
  explicit `created`, `revision`, and package `version` labels without replacing the generation
  description defined by `docker/Dockerfile`.

These facts rule out “add four tag patterns to the existing main workflow” as a complete fix.
That approach either stays unreachable or requires tag builds, which creates a second, untested
artifact instead of promoting the image that earned `Docker Success`.

## 5. Architecture and boundaries

```text
push to main at exact SHA
  |-- CI Tests ------------------------------> check-run: CI Success
  `-- Docker Build -> build -> tests -> push -> check-run: Docker Success
                                    |
                                    `-> ghcr.io/hassansaei/vntyper:sha-<short SHA>
                                           labels: revision=<full SHA>
                                                   version=<package version>

pre-existing strict vX.Y.Z tag                 workflow_dispatch candidate
  |                                                        |
  +-- authenticated repository_dispatch                    |
  |      event_type=vntyper_release                         |
  |      client_payload.tag=<tag>                           |
  `---------------- default-branch publish-pypi.yml --------'
                         |
                         +-- validate tag/candidate, version, SHA on main
                         +-- bounded exact-SHA check-run polling
                         +-- inspect SHA image and resolve digest
                         +-- build/check Python distributions (unprivileged)
                         +-- compute alias plan
                         |
                         +-- repository dispatch: promote digest in GHCR
                         |              `-> exact aliases, then monotonic floating aliases
                         |              `-> protected OIDC PyPI publish
                         |
                         `-- dispatch: summary only; no production writes
```

### 5.1 Docker build boundary

`.github/workflows/docker-build.yml` continues to own source checkout, image construction, smoke
and Docker tests, and the first push. It still does not publish PR images. For a main SHA it emits:

- rolling `ghcr.io/hassansaei/vntyper:main`;
- existing-compatible `ghcr.io/hassansaei/vntyper:sha-<7 lowercase hexadecimal characters>`;
- `org.opencontainers.image.revision=<the same 40-character SHA>`;
- `org.opencontainers.image.version=<vntyper.version.__version__>`.
- a run-attempt-scoped `docker-release-evidence-<full SHA>-<run attempt>` artifact containing contract version, full
  SHA, the post-push registry-reported manifest digest, workflow run ID/attempt, revision label, and package-version label.
  It uses `actions/upload-artifact@v5`, fails when the file is absent, and has explicit 90-day retention so an
  operator-created tag can consume evidence after the originating main run.

The package version is read from `vntyper/version.py` without importing the application. The build
passes only the dynamic OCI keys it owns (`created`, `revision`, and `version`) so the Dockerfile's
title and renamed description remain authoritative. The dead tag metadata rule is removed. A main
build never emits a production semantic alias or `latest`.
Nightly and manual runs retain the full test tier, but scheduled and manual Docker validation never publish application tags.
Only an exact push to `refs/heads/main` publishes the rolling `main`, short-SHA tag, and matching evidence.

The `changes` job remains a pull-request cost optimization only. Every push to `main`
runs the substantive CI jobs and `build-and-test`, even when its diff is documentation
only, so a green aggregate can never stand in for skipped release evidence and every
post-milestone main SHA has a tested `sha-<short>` image. The main-push exemption applies
to `base-status`, `build-base` dependency evaluation, and `build-and-test` together; it
cannot bypass only the final job while leaving an upstream `needs` job skipped. CI likewise runs lint, type
checking, all four unit-test interpreters, and strict docs on every main push. PRs keep
the existing path-aware behavior.

### 5.2 Pure policy boundary

A focused Python policy helper under `scripts/` owns decisions but performs no GitHub, registry,
filesystem mutation, subprocess, or network I/O. Its tested interfaces cover:

- strict parsing of `vMAJOR.MINOR.PATCH` and plain `MAJOR.MINOR.PATCH`;
- selection and classification of the latest GitHub Actions check-run for each required name;
- derivation of the five aliases;
- exact-alias immutability decisions;
- floating-alias semver comparison and anti-downgrade decisions;
- a structured release plan suitable for `$GITHUB_OUTPUT` and the step summary.

Runtime collection protocols (`Mapping` and `Sequence`) come from `collections.abc`; only `Literal` comes from
`typing`, preserving the Python 3.10 floor without deprecated runtime typing aliases.

Workflow steps own `gh api`, Git fetches, `docker buildx imagetools inspect/create`, artifact
upload/download, sleep, and logging. They convert external observations into policy inputs and
execute only an approved plan. This keeps release judgment independently unit-testable without
shipping it inside the `vntyper` runtime package.

### 5.3 Release coordinator boundary

`.github/workflows/publish-pypi.yml` preserves its trusted-publisher filename and owns five jobs:

1. `validate-release` — strict event/action/payload/input validation; fetch full history; verify the
   existing tag, commit, `main` ancestry, version sources, and event mode.
2. `wait-for-release-gates` — bounded polling of the exact-SHA substantive component jobs and
   `CI Success`/`Docker Success`, then inspection of the SHA image's digest and OCI identity.
3. `build-package` — unprivileged checkout of the verified SHA, `python -m build`, `twine check`,
   and upload of `dist/` as a workflow artifact.
4. `promote-ghcr` — calculate and execute the alias plan from the resolved digest. It runs only
   for `repository_dispatch` action `vntyper_release` and only after gates and package build succeed.
5. `publish-pypi` — download the built distributions and invoke the protected OIDC publisher. It
   runs only for `repository_dispatch` action `vntyper_release` and needs validation, gates, package
   build, and GHCR promotion.

A final always-running summary job depends on every prior job with job-level `if: always()` and composes their explicit
JSON outputs plus `needs.<job>.result`; individual job summaries are not treated as cross-job state.
Each fallible phase has a separate `if: always()` result serializer. Polling persists the last snapshot before deciding;
the evidence preflight persists `pending`, `eligible`, `ineligible`, or `infrastructure-failure` plus ordered candidate
runs and the selected run; evidence persists selected-run/validation state; promotion atomically records `attempted`,
then `write_succeeded`, then `verified` for each alias so a successful write followed by failed inspection remains visible; publication
records pre-publish PyPI existence plus action outcome. Terminal failures therefore retain structured last state.
Whole-workflow concurrency is per version with `cancel-in-progress: false`. The `promote-ghcr` job additionally uses one fixed
repository-wide concurrency group with `cancel-in-progress: false`, so promotions cannot execute concurrently between
inspecting and updating floating aliases. GitHub retains at most one pending run per group: a newer pending promotion
may cancel an older pending promotion before it acquires write authority. That canceled version is retried explicitly;
idempotent alias/package rules make retry safe. This group is a mutual-exclusion lock, not an unbounded queue.

## 6. Interfaces, inputs, outputs, and invariants

### 6.1 Production repository-dispatch input

- Event: authenticated `repository_dispatch` with exact type/action `vntyper_release`. The workflow
  file and controller checkout come from the default branch; there is no production `push.tags`
  trigger.
- Candidate input: `github.event.client_payload.tag`, naming an already-existing tag. The payload
  must match `^v(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$` before it is used in Git or shell
  commands.
- API interface: an authorized release operator sends JSON
  `{"event_type":"vntyper_release","client_payload":{"tag":"vX.Y.Z"}}` to the repository-dispatch
  endpoint only after the tag exists. No boolean publish switch or candidate SHA is accepted.
- `github.sha` identifies the default-branch coordinator for this event and is not candidate
  identity. The candidate SHA is obtained only by peeling `refs/tags/<payload tag>^{commit}`.
- That full SHA must be an ancestor of fetched `origin/main`; a commit that exists only on another
  branch is rejected.
- The tag's plain version must equal the executed value of `vntyper/version.py::__version__`.
- The existing version-consistency test must also prove `CITATION.cff` and the current changelog
  version agree. A mismatch fails before polling or publication.
- Validation summary state retains the expected tag version, each observed package, `CITATION.cff`, and current
  changelog version, a Boolean match verdict for each, and the version-consistency pytest exit/verdict. An early setup
  failure uses explicit unavailable values rather than dropping the fields.
- Every false package, citation, or changelog match verdict fails validation even when
  `test_version_consistency.py` passes because the three repository sources agree with each other.

An operator creates or pushes the immutable candidate tag separately, then sends the authenticated
repository dispatch. The workflow never creates, force-updates, deletes, or pushes a Git tag and
never creates a GitHub Release.

### 6.2 Manual dry-run input

`workflow_dispatch` accepts one required existing strict `vMAJOR.MINOR.PATCH` tag. It resolves the
tag to a 40-character candidate SHA and plain version. Manual mode:

- keeps the current workflow/controller checkout separate from the candidate checkout, so coordinator helpers remain
  available even when inspecting a pre-milestone tag while version tests execute against candidate files;
- verifies the SHA is on `main` and the candidate version equals package metadata at that SHA;
- performs an evidence-contract preflight before the long check poll: completed successful main-push Docker runs with
  no attempt-qualified contract artifact are explicitly pre-contract/ineligible; no completed run remains pending.
  Eligible runs are ordered by descending numeric run ID and the exact attempt-qualified artifact search falls back in
  that order; the structured preflight state preserves every eligible run and the selected run;
- polls the same exact-SHA checks and inspects the same SHA image;
- builds and checks the distributions;
- computes the exact alias and anti-downgrade plan;
- writes the complete plan to the summary;
- skips GHCR login/write, alias creation, OIDC permission, PyPI publication, tag creation, and
  GitHub Release creation unconditionally.

There is no `publish` or `dry_run=false` input. Production behavior is reachable only from exact
`repository_dispatch` type `vntyper_release`; `workflow_dispatch` remains read-only for every tag.

### 6.3 Check-run interface

The required contract is the GitHub Actions check-runs named exactly:

- `Lint (Ruff)`
- `Type Check (mypy)`
- `Unit Tests (Python 3.10)`
- `Unit Tests (Python 3.11)`
- `Unit Tests (Python 3.12)`
- `Unit Tests (Python 3.13)`
- `Docs build (strict)`
- `CI Success`
- `Build and test image`
- `Docker Success`

For each name the coordinator reads check-runs for the full candidate SHA, filters to the GitHub
Actions app, and evaluates the newest run attempt. A success from a different SHA, older attempt,
or similarly named job is irrelevant.

Polling is bounded to 120 attempts at 30-second intervals: no more than 60 minutes of deliberate sleep. The job timeout
is at least 65 minutes so setup/API overhead cannot preempt the policy timeout.
Each Check Runs API snapshot is itself bounded to three attempts with five seconds between failures. The classifier
receives the configured maximum explicitly, and the shell exits nonzero after the outer loop even if a future
classifier/configuration mismatch returns `wait` on the last iteration.
Missing, queued, requested, waiting, or in-progress checks continue polling. `failure`,
`cancelled`, `timed_out`, `action_required`, `startup_failure`, `stale`, or another terminal
non-success conclusion fails immediately. Exhausting the bound fails with the last observed state
and check-run URLs. Every named check must be `completed/success` before image promotion or PyPI
publication. A `skipped` component is terminal failure for release purposes even if its aggregate
is green. The always-run-on-main workflow contract is guarded by unit tests so a path filter
cannot silently restore green-on-skip release behavior.

### 6.4 Image input and promotion output

The human-readable release image reference is:

`ghcr.io/hassansaei/vntyper:sha-<first 7 characters of candidate SHA>`

The authoritative input is the post-push registry digest from the successful exact-SHA Docker workflow's
run-attempt-scoped evidence artifact. After its explicit push, the Docker workflow reads the registry-reported
`{{.Manifest.Digest}}` for `sha-<7>` and verifies that registry object's revision/version labels before upload. The coordinator orders all
completed-success main-branch `push` runs of `docker-build.yml` whose `head_sha` equals the candidate by descending
numeric run ID, selects the first with its exact run-attempt-qualified artifact, downloads only that artifact, and
validates its contract version, run ID/attempt, SHA, digest, revision, and version. Before
promotion, inspection must prove:

- the tag exists;
- it resolves to one non-empty `sha256:` manifest digest;
- OCI `org.opencontainers.image.revision` equals the full candidate SHA;
- OCI `org.opencontainers.image.version` equals the plain release version;
- the short tag points to that same digest, or its own revision label is a different 40-character SHA sharing the
  candidate's seven-character prefix; only that case is a collision, while any other mismatch fails as unexplained
  drift, and the immutable evidence digest is used only after its full-SHA/label proof;
- the manifest is the artifact published by the successful Docker workflow for that SHA.

Every target is created from `ghcr.io/hassansaei/vntyper@sha256:...` with manifest-level tooling and
`docker buildx imagetools create --prefer-index=false`, preserving a single manifest as a carbon copy instead of
silently wrapping it in a new index.
No layer is rebuilt or pulled into a local daemon. All aliases therefore resolve to the identical
tested digest.

Manual dispatch executes one untagged `imagetools create --dry-run --prefer-index=false` probe against the immutable
source digest. The probe is read-only, has no `--tag`, and proves the Buildx contract without registry login or writes.

Alias invariants:

- `vX.Y.Z` and `X.Y.Z` are exact. Absent means create; same digest means no-op; a different digest
  is an invariant violation and hard failure. They are never silently repointed.
- `X.Y`, `X`, and `latest` are floating. Absent means create. An existing lower semantic version
  advances to the candidate. The same version/digest is a no-op.
- A floating alias with the candidate version but a different digest is `fail-conflict`; it is never treated as a
  no-op or silently overwritten.
- If a floating alias already reports a higher semantic package version, the candidate is an old
  rerun: skip that alias and emit an explicit anti-downgrade notice. Do not fail the otherwise
  valid rerun and do not move the alias backward.
- An existing floating alias labelled exactly `main` is the recognized legacy rolling `main` state
  and advances to the evidence-verified release digest during the first migration.
- Any other missing or unrecognized version label fails closed before every alias write. The error
  requires an operator to repair or remove that floating alias before retrying; a release cannot
  remain green while one of its five required aliases stays stale.
- Promotion never reads from or writes `main`.

Exact aliases execute before floating aliases. This ensures that a release always acquires its
immutable identity even when a newer release already owns all floating aliases.

If the evidence artifact or source manifest is absent, the release stops before package build or any write and names
the exact `Docker Build` run to rerun. An evidence payload with the wrong contract version is explicitly ineligible.
Every eligible post-milestone `main` push is required to build regardless of paths, so rerunning that existing run is a
real recovery path. Pre-milestone SHAs lack this artifact and are dry-run diagnostics only; production rejects them with
the explicit ineligibility reason. No alternate image or `main` tag may be substituted.

### 6.5 Python distribution output

The build job checks out the verified SHA with persisted Git credentials disabled, runs on Python
3.12 inside an explicit virtual environment, installs only build/check tooling, runs
`python -m build` and `twine check dist/*`, and uploads only the resulting wheel and sdist. It has
`contents: read` and no `id-token: write`, environment, package-write permission, or package secret.
It uses `actions/upload-artifact@v5` with seven-day retention; the OIDC job consumes the exact artifact name with
`actions/download-artifact@v5`.

The publish job has no source checkout or build toolchain. It downloads that workflow artifact,
uses GitHub environment `pypi`, and grants only `id-token: write` plus permissions required to read
the workflow artifact. It invokes:

`pypa/gh-action-pypi-publish@dc37677b2e1c63e2034f94d8a5b11f265b73ba33`

That SHA is the reviewed `v1.14.2` release at this design baseline. The action receives
`skip-existing: true`, preserving the current recovery contract. It receives no password,
username, API token, repository secret, or GitHub checkout credential.

### 6.6 Naming outputs and protected identities

The nine hand-edited targets are:

| Target at baseline | Required generation prose |
| --- | --- |
| `README.md:3` heading | `VNtyper 2` |
| `README.md:5` opening name | `VNtyper 2` |
| `README.md:61` packaging sentence | `VNtyper 2` |
| `README.md:89` subcommands sentence | `VNtyper 2` |
| `README.md:133` Docker sentence | `VNtyper 2` |
| `README.md:209` pipeline overview sentence | `VNtyper 2` |
| `README.md:222` dependencies sentence | `VNtyper 2` |
| `README.md:323` citation sentence | `VNtyper 2` |
| `docker/Dockerfile:24` OCI description | `VNtyper 2` |

The grammar at `README.md:5` becomes a valid sentence equivalent to “This refactored version of
VNtyper v1 integrates …”; it does not rename `VNtyper v1`.

Protected values include, but are not limited to:

- `vntyper/version.py`, `CITATION.cff`, and every `docs/about/changelog.md` version heading;
- dependency and environment numbers in `pyproject.toml`, `conda/`, and Docker requirements;
- historical statements such as “Before VNtyper 2.0.6/2.0.8/2.0.9” and their test fixtures;
- distribution/import/command `vntyper`, repository names, version tags, and image aliases;
- `CITATION.cff` title, `mkdocs.yml` site name, `pyproject.toml` description, OCI title, CLI/API
  titles, and report labels that use the canonical bare product name;
- `VNtyper Version:` machine banners and their exact golden-cohort normalization anchor;
- `snakemake/vntyper2.smk` and `snakemake/run_vntyper2.sh`.

## 7. Success, partial failure, recovery, retry, and idempotency

| Phase | Success | Partial failure and recovery | Retry behavior |
| --- | --- | --- | --- |
| Validation | Strict version/tag, exact SHA, main ancestry, and all metadata agree. | No external artifact changed. Correct the version/ref; do not retag an already published package version. | Deterministic for one ref and repository state. |
| Gate wait | Every named substantive and aggregate exact-SHA check finishes successfully within the bound. | Failed or skipped component stops immediately; missing/pending checks time out with links and last states. Re-run or repair the failing upstream workflow. | A rerun observes existing successes and proceeds without waiting. |
| Image identity | SHA tag exists and digest/revision/version all agree. | Nothing promoted. Re-run the main Docker workflow for that SHA; never substitute another tag. | Reinspection is read-only and deterministic. |
| Package build | Wheel and sdist build and pass `twine check`. | Nothing promoted or published because promotion depends on the build. Fix packaging and create a new version if immutable publication already occurred. | Rebuilds from the same verified SHA. |
| Exact promotion | Both exact aliases point to the tested digest. | One exact alias may exist before a later command fails. The summary names completed operations. | Same-digest aliases are no-ops; different-digest aliases hard-fail, so rerun converges without mutation. |
| Floating promotion | Eligible aliases and legacy `main` advance; newer aliases skip with notice; every other unorderable label hard-fails before writes. | Failure after some advances leaves every changed alias on the intended tested digest. | Re-evaluation makes same-version operations no-ops and never downgrades. |
| PyPI publish | OIDC action uploads missing distributions and attestations in environment `pypi`. | GHCR may already be complete. PyPI may contain one distribution from a partial upload. Do not roll back immutable aliases. | `skip-existing: true` ignores already-uploaded files and uploads only missing ones. |
| Manual dry run | Full validation, gate/image inspection, package check, and plan summary complete. | No production mutation exists to roll back. | Repeatable; it cannot cross into production mode. |

A PyPI failure after image promotion is an explicitly supported partial state. GHCR exact aliases
are immutable evidence for the release, so they remain. The next workflow rerun verifies the same
digest, treats promotion as a no-op, and retries OIDC publication with `skip-existing`.

## 8. Python and environment constraints

- All repository Python, including the pure release-policy helper, supports Python 3.10. Do not
  use 3.11-only syntax or standard-library APIs. Unit tests run across the existing 3.10–3.13
  matrix.
- GitHub build jobs may select Python 3.12, matching the image/test convention, but repository
  helper code remains compatible with the declared floor.
- CI creates an explicit virtual environment (`uv venv`/`VIRTUAL_ENV` and `.venv/bin`, or an
  equivalent `setup-python` venv). It never uses `uv pip install --system` on Ubuntu 24.04.
- Candidate validation's explicit venv installs `pytest`, `packaging`, `PyYAML`, and `requests`; `requests` is required
  at collection time by repository `tests/conftest.py` even when only the version-consistency node is selected.
- Tests run from the repository root. Every new unit test module declares
  `pytestmark = pytest.mark.unit`; pure policy tests use no network, Docker daemon, registry, or
  test-data archive.
- No conda environment, dependency pin, Java/Kestrel version, Docker base hash, or runtime package
  installation behavior changes in this track.
- Because `scripts/` becomes a production mypy and coverage source elsewhere in milestone 6, the
  release-policy helper must be fully annotated and decision branches must have direct unit tests.

## 9. Security and GitHub permission model

Workflow-level permissions remain `{}`. Each job receives only what it uses:

| Job | Permissions | Prohibited capability |
| --- | --- | --- |
| `validate-release` | `contents: read` | no checks, package write, OIDC, or secrets |
| `wait-for-release-gates` | `actions: read`, `contents: read`, `checks: read`, `packages: read` | no package write or OIDC |
| `build-package` | `contents: read`, workflow artifact access | no OIDC, environment, or package secret |
| `promote-ghcr` | `contents: read`, `packages: write` | no OIDC or PyPI credential |
| `publish-pypi` | `id-token: write`, workflow artifact access | no checkout credential, package write, or long-lived secret |

- No `pull_request_target` event is introduced. PRs and fork PRs cannot reach production jobs.
- The production event is authenticated `repository_dispatch` with exact type/action
  `vntyper_release`. A dispatch or tag alone is not authority: the candidate tag must already
  exist, and strict syntax, `main` ancestry, exact-SHA checks and evidence, package version, and
  image identity all fail closed.
- `GITHUB_TOKEN` authenticates GHCR. Docker Hub credentials are unused and must not be added.
- Production gate inspection performs an explicit `docker/login-action@v4` GHCR login before Buildx reads; promotion
  performs its own write-authorized login. Manual dispatch performs no GHCR login and remains read-only.
- Shell steps place validated values in environment variables, quote expansions, use
  `set -euo pipefail`, and never interpolate unvalidated tag/input strings into executable code.
- The PyPI job uses the protected `pypi` environment. Its OIDC token is short-lived and job-bound;
  the SHA-pinned action prevents a moving action tag from changing privileged code.
- `workflow_dispatch` has no publish switch. Every production job has the exact guard
  `github.event_name == 'repository_dispatch' && github.event.action == 'vntyper_release'`, so a
  maintainer cannot turn a dry run, an unrelated dispatch, or a push into publication.
- Historical tagged commits retain their historical workflow files. Removing `push.tags` from the
  default branch does not disable a tag-triggered token workflow stored in an older commit. During
  first-release migration, do not create or push a tag pointing at a pre-milestone commit; those
  legacy workflows become inert only after #218's `PYPI_API_TOKEN` deletion, which must happen only
  after a successful production OIDC publication.
- Existing fork behavior for missing GHCR base images remains unchanged: a read-only fork token
  cannot publish the base and receives the current explicit diagnostic.

## 10. Observability and diagnostic output

The final `release-summary` job writes a single structured GitHub step summary containing:

- event mode (`production repository dispatch` or `dry run`), tag/candidate version, full SHA, and main-ancestry
  result;
- package, citation, and changelog version verdicts;
- evidence-preflight state/reason, ordered eligible Docker runs, and selected run when one exists;
- required check name, selected check-run URL/attempt/status/conclusion, attempt count, and elapsed
  polling time;
- SHA image reference, immutable digest, OCI revision, and OCI package version;
- Docker workflow run ID/attempt/URL, evidence contract version, and provenance verdict;
- each alias, its previous digest/version if present, decision (`create`, `advance`, `no-op`,
  `skip-newer`, or `fail-conflict`), its attempted/write-succeeded/verified state, and final digest
  where known;
- package artifact filenames and `twine check` result;
- PyPI mode/result and whether the candidate version existed immediately before the OIDC action, allowing the summary
  to distinguish `already-existed-skip`, `published`, and `failed`;
- a prominent statement that dry run performed no registry, PyPI, tag, or release mutation.

Anti-downgrade decisions emit `::notice::` annotations. Unorderable-floating labels, identity
mismatches, exact-alias conflicts, failed checks, and poll exhaustion emit `::error::` annotations
with the exact offending name/ref and a recovery instruction. No token, OIDC assertion, secret,
registry password, or signed upload URL is logged.

## 11. Test strategy and objective acceptance checks

### 11.1 Pure policy tests

Unit tests must prove:

- strict tags accept ordinary `v2.0.10`-style versions and reject missing `v`, prereleases,
  build metadata, leading-zero components, extra components, and shell metacharacters;
- required-check evaluation is exact-SHA and exact-name, selects the newest GitHub Actions
  attempt, distinguishes pending/success/terminal failure, and reaches timeout after the stated
  bound;
- alias derivation returns exactly `vX.Y.Z`, `X.Y.Z`, `X.Y`, `X`, `latest` in deterministic order;
- exact aliases cover absent, same-digest, and conflicting-digest states;
- floating aliases cover absent, older, equal, newer, recognized legacy `main`, and missing/unparseable label states;
- recognized legacy `main` advances, while every other missing/unparseable label hard-fails before writes;
- equal-version floating aliases cover same-digest no-op and different-digest hard conflict;
- a rerun after every possible prefix of completed alias operations converges on the same digest;
- dry-run plans contain the same decisions but no executable mutation operations.

### 11.2 Workflow contract tests

Static executable tests must parse the workflows and prove:

- Docker main pushes publish `main` plus the existing short-SHA tag with a full revision label,
  never semantic aliases or `latest`; schedules and manual runs retain tests but publish no application tag;
- tag builds are not introduced into `docker-build.yml`;
- revision and package-version labels are explicit and the Dockerfile description is not
  overwritten by repository metadata;
- `publish-pypi.yml` has only `repository_dispatch.types: [vntyper_release]` and manual triggers,
  no `push.tags` trigger, and manual mode cannot satisfy any production job guard;
- production resolves only the existing strict `github.event.client_payload.tag`; the controller
  is the default-branch workflow checkout and every production job guard requires both exact
  event name and action;
- `publish-pypi.yml` retains its exact filename and environment `pypi`;
- Docker evidence and package uploads use `actions/upload-artifact@v5` with explicit 90-day and seven-day retention,
  respectively, and PyPI uses `actions/download-artifact@v5` with the exact package artifact name;
- all ten exact component/aggregate check-run names and the finite polling constants are present;
- Check Runs API retries are finite, the classifier receives the configured bound, and loop exhaustion cannot fall
  through as success;
- promotion consumes the immutable digest bound to the successful exact-SHA Docker evidence artifact, verifies the
  short-tag/full-revision/package-version relationship, and creates exactly five aliases with
  `--prefer-index=false`; an executable shell test and the read-only Buildx dry-run probe enforce this rather than a
  source-string assertion;
- `main` is absent from release-promotion destinations;
- package build and PyPI publish are separate jobs; publish needs validation, gates, build, and
  promotion;
- only the publish job has `id-token: write`; only promotion has `packages: write`;
- no live `PYPI_API_TOKEN`, `TWINE_PASSWORD`, Docker Hub secret, or unpinned PyPA action remains;
- the pinned PyPA action SHA and `skip-existing: true` are present;
- `actionlint` accepts both changed workflows.

### 11.3 Registry documentation and naming tests

Tests must prove:

- active install surfaces in `README.md`, `docs/getting-started/installation.md`, and
  `docs/user-guide/docker.md` use `ghcr.io/hassansaei/vntyper`, with no active
  `saei/vntyper` pull/run/Apptainer command;
- documented current aliases are members of the release policy and documented `main` semantics remain
  rolling/development rather than stable/latest. Until the first gated release creates `latest`, runnable install
  examples use `main` or an exact existing SHA and explicitly label the transition;
- the nine explicit targets contain `VNtyper 2` and no longer contain generation prose
  `VNtyper 2.0`;
- the README grammar repair retains `VNtyper v1`;
- authoritative version metadata still agrees through
  `tests/unit/test_version_consistency.py`;
- protected bare names, parsed `VNtyper Version:` anchor, historical version sentences, and
  `snakemake/vntyper2*` filenames remain unchanged.

The test must not assert a repository-wide count of `VNtyper 2.0.x` historical prose. It owns the
explicit rename targets and protected identities instead.

### 11.4 Required verification commands

From the repository root in the project environment:

```bash
make format-check
make lint
make type-check-all
make test-unit
make test-unit-cov
make patch-coverage
make docs-build
make ci-local
make lint-actions
```

`make ci-local` is mandatory because this track changes workflows and must exercise CI's clean
environment installer. Docker workflow changes additionally require `make ci-local-docker` when a
Docker daemon is available; unavailability is reported rather than silently skipped. The PR's
required `CI Success` and `Docker Success` checks must both complete successfully.

## 12. Approaches considered

### 12.1 Selected: default-branch repository dispatch in the existing PyPI workflow

This design preserves the configured trusted-publisher filename, guarantees that current default-branch
policy validates every production request, reuses the image already proven by Docker tests, separates
irreversible privileges, and makes partial reruns convergent. It adds coordination complexity, but every
wait is bounded and every decision is testable and summarized.

### 12.2 Rejected: add semver metadata rules and build on tag

Adding `type=semver` and a tag trigger to `docker-build.yml` is syntactically small, but it builds
a new tag-context artifact after the main image was tested. The tag build can differ in metadata,
base resolution, timing, or source, and its success is not the `Docker Success` already required
for the main commit. It also leaves PyPI able to publish before CI/Docker finish and provides no
cross-registry partial-failure contract.

### 12.3 Rejected: coordinate production directly from a tag push

GitHub selects a tag-push workflow from the tagged commit. A tag pointing at historical code can
therefore select historical publication policy and bypass every fail-closed check added later on the
default branch. Exact-SHA validation inside the current coordinator cannot repair a coordinator that
was never selected. The narrow repository-dispatch event keeps the candidate tag immutable while
loading release authority from the default branch.

### 12.4 Rejected: publish release aliases from every main build

Making default-branch metadata emit `latest` would cause unreleased commits to replace the stable
install target. Deriving major/minor aliases from the package version on every main push would
likewise publish release-looking tags before an operator creates a release. This collapses the
deliberate distinction between rolling `main` and production aliases and makes rollback/rerun
behavior unsafe.

### 12.5 Rejected: mirror GHCR to Docker Hub

Dual publication expands secrets, permission surfaces, rate limits, partial-failure combinations,
and documentation drift. Current automation has already removed Docker Hub, its live artifacts are
stale, and #214 explicitly permits removal/legacy status. GHCR-only publication is the smallest
coherent registry contract.

## 13. Dependencies

- The PyPI trusted publisher and GitHub `pypi` environment must remain named exactly as confirmed
  in #218's maintainer comment. Renaming `publish-pypi.yml` or the environment breaks OIDC trust.
- The main `Docker Build` workflow must publish its short-SHA image and package-version/full-revision
  labels before the release coordinator can promote it.
- Branch protection/check naming must retain `CI Success` and `Docker Success`, or the workflow
  contract and tests must change together before a release.
- The release-policy helper depends on milestone-6 work that brings `scripts/` under production
  mypy and branch coverage; it must satisfy those gates when introduced.
- A repository owner is needed after the first successful OIDC release to delete
  `PYPI_API_TOKEN`. The implementation PR performs no production dispatch; until deletion, operators
  must not create or push tags at commits containing the legacy token workflow.
- There are no formal GitHub dependency links among #214, #218, and #220. They are implemented
  together because they share release workflows, published image metadata, install prose, and
  tests—not because one issue blocks another.

## 14. Rollout

1. Land implementation and tests without changing the package version or pushing any tag.
2. Merge only after local workflow gates and GitHub `CI Success`/`Docker Success` pass.
3. Use `workflow_dispatch` with an existing pre-milestone strict tag only to prove the explicit ineligibility path and
   zero production mutations. Do not describe it as a full-path success. Do not create or push any new tag pointing at
   a pre-milestone commit while `PYPI_API_TOKEN` still exists: that commit's historical tag workflow, not the current
   default-branch coordinator, would be selected.
4. Prepare the next normal patch release on `main`, updating the repository's three authoritative
   version surfaces through the established process.
5. After that exact main SHA earns both required checks, an operator creates and pushes the strict
   candidate tag or creates a GitHub Release targeting that SHA. Because that commit contains the
   new coordinator without `push.tags`, tag creation alone does not publish. The workflow itself
   never creates or moves the tag.
6. Send an authenticated `repository_dispatch` with `event_type=vntyper_release` and
   `client_payload.tag` equal to that existing candidate tag. Observe exact-alias promotion,
   floating-alias decisions, protected `pypi` approval, OIDC upload,
   PyPI attestations, and the final summary. Verify all GHCR aliases resolve to the recorded digest
   and PyPI files correspond to the workflow artifact.
7. Manually dispatch that newly eligible existing tag and verify the complete dry-run path reports the same evidence and alias
   plan while performing no writes.
8. Only after the first successful real OIDC publication, a repository owner deletes the obsolete
   `PYPI_API_TOKEN` secret and records that action on #218. Legacy tag-triggered workflows stored in
   historical commits are considered inert only after this deletion.
9. Close #214, #218, and #220 only with links to passing contract tests, the production release
   run, published GHCR aliases, PyPI OIDC evidence, and the completed secret-removal follow-up.

Rollback of workflow or documentation code uses an ordinary PR. Published PyPI files and exact
GHCR aliases are immutable by policy and are not rolled back or overwritten; a defect requires a
new patch version. Floating aliases may advance to that new patch but never move backward.

## 15. Unresolved questions

None. The approved decisions are: authenticated default-branch `vntyper_release` repository dispatch
over a pre-existing tag; exact-SHA gated digest promotion; monotonic release aliases; rolling `main`;
GHCR as the authoritative registry; Docker Hub as unsupported legacy; separated OIDC publication
with the existing filename/environment; manual dry-run only; and exactly nine generation-prose
renames while versions and identifiers remain unchanged.
