# Milestone 6 End-to-End Superpowers Execution Prompt

Copy the prompt below into a fresh agent session rooted at
`/home/bernt-popp/development/VNtyper/.worktrees/ci-milestone-6-gates-harness-release` after the
reviewed planning commit is present. The driver is Codex: repository edits use `apply_patch`,
Superpowers agents use the Codex plugin skills, GitHub publication uses the named `github:*`
skills, and Claude Code is invoked only as the independent read-only Opus reviewer.

---

Act as the principal engineer responsible for executing VNtyper milestone #6 end to end and
opening a fully evidenced pull request. Continue without routine approval pauses. Do not merge
the PR, create/push a version tag, create a GitHub Release, publish PyPI/GHCR artifacts, delete
credentials, force-push, or discard unrelated user changes.

Repository:

```text
/home/bernt-popp/development/VNtyper/.worktrees/ci-milestone-6-gates-harness-release
```

Target branch and existing isolated worktree:

```text
branch:   ci/milestone-6-gates-harness-release
worktree: /home/bernt-popp/development/VNtyper/.worktrees/ci-milestone-6-gates-harness-release
base:     ebb15b26631242a3295607e4eda4c68f688cd9a2
```

Issues:

```text
#71 #204 #208 #211 #214 #218 #219 #220 #226
```

## Mandatory authority and skills

Before changing repository files:

1. Read the repository's complete `AGENTS.md`; it is the repository-level authority.
2. Invoke and completely read the current `superpowers:using-superpowers` skill.
3. Read and follow these current skills at their applicable gates:
   - `superpowers:using-git-worktrees`
   - `superpowers:dispatching-parallel-agents` for the explicitly independent read-only audit and review waves
   - `superpowers:subagent-driven-development`
   - `superpowers:test-driven-development`, including `writing-good-tests.md`
   - `superpowers:systematic-debugging` for every unexpected failure
   - `superpowers:receiving-code-review`
   - `superpowers:requesting-code-review`
   - `superpowers:verification-before-completion`
   - `superpowers:finishing-a-development-branch`
   - `github:yeet` for branch push and PR creation
   - `github:gh-fix-ci` for failing PR checks
   - `github:gh-address-comments` for actionable PR feedback
4. Do not rerun brainstorming or redesign the accepted architecture. The specifications and plans
   below have already passed brainstorming, reconstruction, and adversarial design review.

## Required reading

Read every file completely before execution:

Specifications:

```text
docs/superpowers/specs/2026-08-10-milestone-6-program-design.md
docs/superpowers/specs/2026-08-10-milestone-6-quality-gates-design.md
docs/superpowers/specs/2026-08-10-milestone-6-harness-matrix-design.md
docs/superpowers/specs/2026-08-10-milestone-6-mutation-safety-design.md
docs/superpowers/specs/2026-08-10-milestone-6-release-and-naming-design.md
docs/superpowers/specs/2026-08-10-milestone-6-exception-policy-design.md
```

Plans:

```text
docs/superpowers/plans/2026-08-10-milestone-6-integration-plan.md
docs/superpowers/plans/2026-08-10-milestone-6-quality-gates-plan.md
docs/superpowers/plans/2026-08-10-milestone-6-harness-matrix-plan.md
docs/superpowers/plans/2026-08-10-milestone-6-mutation-safety-plan.md
docs/superpowers/plans/2026-08-10-milestone-6-release-and-naming-plan.md
docs/superpowers/plans/2026-08-10-milestone-6-exception-policy-plan.md
```

The master integration plan governs ordering, shared interfaces, traceability, final verification,
and PR delivery. A bounded track plan governs the detailed task when it is more specific. If two
plan requirements genuinely conflict, batch all conflicts into one concise user question before
implementation, as required by the current SDD skill.

## Preflight

Perform read-only checks and record the evidence:

```bash
git -C /home/bernt-popp/development/VNtyper status --short --branch
git -C /home/bernt-popp/development/VNtyper/.worktrees/ci-milestone-6-gates-harness-release status --short --branch
git -C /home/bernt-popp/development/VNtyper/.worktrees/ci-milestone-6-gates-harness-release rev-parse HEAD
git -C /home/bernt-popp/development/VNtyper/.worktrees/ci-milestone-6-gates-harness-release merge-base HEAD origin/main
git -C /home/bernt-popp/development/VNtyper/.worktrees/ci-milestone-6-gates-harness-release remote -v
```

Preserve all unrelated changes. If the milestone worktree contains only the reviewed uncommitted
spec/plan package, continue by completing its prescribed review and documentation commits. Do not
reset or recreate it.

## SDD setup and persistent state

Use one current-SDD workspace per bounded plan, not one master ledger:

```bash
cd /home/bernt-popp/development/VNtyper/.worktrees/ci-milestone-6-gates-harness-release
M6_SDD_DIR=/home/bernt-popp/.codex/plugins/cache/openai-curated-remote/superpowers/6.2.0/skills/subagent-driven-development
for M6_PLAN in \
  docs/superpowers/plans/2026-08-10-milestone-6-quality-gates-plan.md \
  docs/superpowers/plans/2026-08-10-milestone-6-mutation-safety-plan.md \
  docs/superpowers/plans/2026-08-10-milestone-6-harness-matrix-plan.md \
  docs/superpowers/plans/2026-08-10-milestone-6-release-and-naming-plan.md \
  docs/superpowers/plans/2026-08-10-milestone-6-exception-policy-plan.md
do
  "$M6_SDD_DIR/scripts/sdd-workspace" "$M6_PLAN"
done
```

Using `apply_patch`, create `progress.md` in each of the five printed directories. Its first line
must name that directory's exact bounded plan, for example:

```text
# SDD ledger — plan: docs/superpowers/plans/2026-08-10-milestone-6-quality-gates-plan.md
```

On every task record:

- base SHA before dispatch;
- implementer agent/model and report path;
- commit(s) and focused test evidence;
- review-package path over exact base..head;
- specification verdict and quality verdict;
- each fix round and scoped re-review verdict;
- final task completion line.

After compaction or restart, trust this ledger plus `git log`; never redispatch a task with a
`Task N: complete` line.

## Execution order and maximum safe parallelism

The installed `superpowers:subagent-driven-development` v6.2.0 skill explicitly says **never
dispatch multiple implementation subagents in parallel**. Honor that current hard rule. Maximize
parallelism elsewhere: independent read-only audits, evidence collection, and external Opus review
may overlap, but each implementation task must complete its implementer → task review → fix loop →
ledger gate before the next implementation task begins.

At the start of each wave, use `superpowers:dispatching-parallel-agents` to run only disjoint read-only preflight
audits (requirements drift, owned-file status, and prior-interface verification) while the controller prepares the next
task brief. At whole-branch review, run the independent Superpowers reviewer and Opus process concurrently. Never let
either reviewer edit, and never overlap two implementers or an implementer with a reviewer of changing files.

Execute the plan family in this exact dependency order:

```text
Wave 0  reviewed specs and plans commits
Wave 1  quality-gates Tasks 1-2: measurement interface and #204 scripts mypy gate
Wave 2  mutation-safety Tasks 1-7: #208
Wave 3  harness-matrix Tasks 1-6: #71 and #226 evidence
Wave 4  release-and-naming Tasks 1-11: #214/#218/#220
Wave 5  resume quality-gates Tasks 3-10: #211 against final scripts
Wave 6  exception-policy Tasks 0-12: #219 against final Python/test tree
Wave 7  whole-branch verification, reviews, fixes, PR, CI monitoring
```

This order is load-bearing: #219 must follow #211 so its exact final inventory is not invalidated by
the largest test tranche. Do not move coverage source inclusion before scripts-only measurement.

For each bounded task, set `BOUNDED_PLAN` to that track's exact path and `N` to the task number:

1. Run `"$M6_SDD_DIR/scripts/task-brief" "$BOUNDED_PLAN" "$N"` and give the implementer only that brief plus frozen interfaces.
2. Require strict RED → observed intended failure → minimal GREEN → refactor while green for behavior changes.
   For a test-only characterization, record the live passing baseline, temporarily make the expectation wrong or apply
   the plan's controlled mutation, observe the intended failure, restore it, and rerun GREEN; never fabricate a RED by
   expecting collection/import failure.
3. Require pure unit tests with meaningful literal value/side-effect assertions and
   `pytestmark = pytest.mark.unit` in every new unit file.
4. Require a focused Conventional Commit and report file.
5. Run `"$M6_SDD_DIR/scripts/review-package" "$BOUNDED_PLAN" "$BASE" "$HEAD"` using the captured base, never `HEAD~1`.
6. Dispatch a fresh task reviewer with the brief, report, review package, and exact global
   constraints. Require both specification and quality verdicts.
7. Use the SDD five-round fix loop. Rounds 1–3 resume the implementer; rounds 4–5 use a fresh,
   stronger implementer. Every fix gets a scoped re-review. Never fix findings in the controller.
8. Advance only after the task review is clean or the skill's breaker is reached and every
   non-load-bearing residual is explicitly parked with a technical ruling. Stop on a load-bearing
   residual.

Always pass an explicit model to every agent dispatch: fast only for complete single-file
mechanical briefs, standard for multi-file work and all ordinary reviewers, most capable for
release/concurrency judgment and the final whole-branch reviewer. Record the model in the ledger.

## Frozen high-risk contracts

Do not regress these Opus-reviewed corrections:

- The actual cleanup artifact is exactly
  `fastq_bam_processing/output_unmapped.bam`. A positive keep test must first prove it exists;
  negative delete/precedence assertions reuse that proven literal. Never derive it from a fixture
  name.
- Every push to `main` runs substantive lint, mypy, unit 3.10–3.13, strict docs, and Docker
  build/test jobs even for docs-only diffs. Main concurrency is per full SHA so a third push cannot
  replace a different pending SHA's evidence. Release polling requires every named component and
  aggregate check to be `completed/success`; `skipped` is terminal failure.
- Release promotion selects a successful exact-SHA main-push Docker run, downloads its run-attempt-qualified full-SHA
  evidence artifact, verifies contract version/run/post-push registry digest/revision/package version, checks the existing
  `sha-<7>` tag for equality or proves a different full revision with the same seven-character prefix, then uses the immutable evidence digest. Do not
  introduce a breaking long-SHA tag scheme.
- Every `docker buildx imagetools create` promotion uses `--prefer-index=false`; the task must execute the plan's
  read-only Buildx contract probe, not merely assert that the string appears in YAML. A successful create is recorded
  before any fallible reinspection, then separately marked verified.
- The GHCR promotion job uses a fixed repository-wide concurrency group with
  `cancel-in-progress: false`; different versions may not execute a read/decide/write alias
  sequence concurrently. GitHub may replace an older pending run; that canceled-before-write
  version must be rerun explicitly. Never call this lock an unbounded queue.
- Mutation execution uses a detached worktree overlaid with the current tracked and non-ignored
  untracked working state plus tracked deletions. Overlay operations are distinct from the authoritative post-overlay
  baseline manifest, so a staged change reverted in the working tree can normalize cleanly relative to HEAD. This
  preserves the existing uncommitted test-killing-mutant loop while keeping mutation writes out of the real source tree.
- `tests/unit/test_golden_cohort_launcher.py` is created and owned by quality-gates Task 5. Exception-policy tasks may
  only modify it and must preserve all existing assertions; Wave 6 finishes by rerunning scripts-only coverage after
  its additions.
- BLE001 inventory failures name both the reviewed Ruff version and actual measuring version. A
  tool-version change alone passes when normalized diagnostic identities remain identical.
- Scripts baseline is 69.01%; scripts-only acceptance is at least 88.00%; combined acceptance is
  at least 86.00% on Python 3.10–3.13. Branch mode and 86/86/80 thresholds remain unchanged.

## Mandatory adversarial reviews

Use a real Claude Code process with exact model `claude-opus-5`, `--effort max`, and read-only plan
permission. Confirm `modelUsage.claude-opus-5.canonicalModel == "claude-opus-5"`; never silently
substitute another model.

Run Opus on actual artifacts, not favorable summaries:

1. After each high-risk implementation track (#208, release, #211, #219), review that track's
   base..head package, exact spec/plan, tests, and focused evidence.
2. After integrated implementation, review merge-base..HEAD plus all specs/plans/workflows/tests.
3. After fresh final gates and the PR diff, review exact results, permissions, shell quoting,
   reruns, partial failures, caches, fork behavior, Python 3.10, marker collection, coverage arcs,
   patch coverage, mutation sensitivity, release races, and undocumented behavior.
4. Resolve every Critical/High through `superpowers:receiving-code-review`; do not implement
   suggestions blindly. Use one coordinated final fix wave and one scoped re-review.

The reviewer must classify Critical/High/Medium/Low with exact paths/symbols. Do not advance with an
unresolved Critical/High.

## Mandatory final verification

Run freshly from repository root after all fixes:

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

Also run the relevant mutation canary/generator round trip. The primary checkout contains the verified git-ignored
archive, but the milestone worktree initially does not. Harness Task 3 must verify and copy that archive into the
worktree without network access, then run `make cram-fixtures` before any bare integration node. After that staging,
run the focused integration tiers and `make check-full` if all external tools remain available. Set
`VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1` so a missing tier fails fast instead of downloading. Never claim an unavailable
tier passed.

Capture command, exit code, test count, precise branch-inclusive coverage, scripts-only coverage,
patch coverage, skipped/unavailable tiers, and artifacts. Re-run all affected and final gates after
any review fix.

## Final review and PR delivery

After all five bounded-plan ledgers are complete:

1. Generate the whole-branch review package from `git merge-base origin/main HEAD` to HEAD.
2. Dispatch the current Superpowers final reviewer on the most capable available model.
3. Run exact Opus 5 maximum-effort whole-branch review.
4. Apply one coordinated final fix wave and one scoped re-review. There is no second final fix
   wave: after the scoped re-review, technically adjudicate residual findings under the SDD skill;
   stop on any unresolved load-bearing finding.
5. Verify every issue criterion against the traceability matrix and fresh evidence.
6. Ensure the working tree is clean and commits are coherent, reviewable Conventional Commits.
7. Use `github:yeet` to push `ci/milestone-6-gates-harness-release` without force and open a PR
   against `main` using the repository template.

The PR body must contain:

- milestone #6 and all nine issue links;
- architecture/dependency summary;
- criterion → implementation → test traceability;
- exact local commands/results, coverage, and patch coverage;
- integration/Docker tiers passed or accurately unavailable;
- compatibility, permission, migration, and release risk;
- Opus findings and resolutions;
- genuine residual limitations;
- `Closes #N` only for issues whose implementation is delivered by this PR. Cite #226's already
  landed implementation/test evidence and milestone disposition, but do not auto-close it from an
  unrelated diff; close that tracking issue separately with its existing commit evidence.

After opening the PR, inspect and monitor GitHub Actions to terminal state. Use
`github:gh-fix-ci` plus `superpowers:systematic-debugging` for failures and
`github:gh-address-comments` for actionable threads. Every repository fix must return to its
owning bounded-plan SDD ledger: append a `PR fix round`, capture the pre-fix SHA, dispatch an
implementer with the exact failure/comment and covering-test requirement, require report/commit,
generate a scoped review package, and obtain a task re-review. The controller never authors the
fix. Then rerun all affected gates plus the complete final gate list and regenerate the
whole-branch review package. Because the single final fix wave has already been consumed, run one
renewed Superpowers and Opus review of the post-PR-fix SHA and adjudicate without another
controller-authored fix wave; stop on any load-bearing finding. The final pushed SHA must equal the
freshly verified and reviewed SHA.

When required checks are terminal/successful, delete only the five git-ignored plan-scoped
scratch directories printed by `sdd-workspace`, invoke
`superpowers:finishing-a-development-branch`, choose the PR/leave-worktree path, and keep the Git
worktree in place for PR iteration.

Do not declare complete until the PR URL exists, required checks are terminal/successful, no
Critical/High finding remains, all applicable fresh gates pass, and the branch is clean/pushed.

Final response must lead with the PR URL and state branch/final commit, all nine issues, specs and
plans used, exact verification results, GitHub checks, Opus verdict, unavailable tiers, and residual
risks. Do not merge or release.

---
