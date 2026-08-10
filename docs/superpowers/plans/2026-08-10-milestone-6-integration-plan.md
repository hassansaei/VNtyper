# Milestone 6 Gates, Harness, and Release Automation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Execute five bounded implementation plans in dependency order, integrate their reviewed commits, verify all nine milestone issues, and open a fully evidenced pull request.

**Architecture:** This file is the program controller, traceability authority, and final-integration plan; it is not a substitute task brief for any bounded implementation plan. Each track owns its own SDD workspace, ledger, task briefs, implementers, review packages, and fix loops. The controller pauses the quality plan after its mypy foundation, executes the three intervening tracks, resumes quality against the final script tree, then executes the exception inventory last.

**Tech Stack:** Python 3.10–3.13, pytest, coverage.py/diff-cover, mypy, Ruff, Git worktrees, Bash, GitHub Actions, Docker Buildx/GHCR, PyPI OIDC, MkDocs, GitHub CLI.

## Global Constraints

- Work only in branch `ci/milestone-6-gates-harness-release`, based on `ebb15b26631242a3295607e4eda4c68f688cd9a2`; never implement on `main`.
- Read the complete `AGENTS.md` and every current required Superpowers/GitHub skill before its stage.
- Preserve Python 3.10 compatibility, Ruff line length 120, double quotes, public Google docstrings, and repository logging/error conventions.
- Run pytest from repository root; every new unit-test file declares `pytestmark = pytest.mark.unit` and remains pure.
- Every behavior-changing task follows RED → observed intended failure → minimal GREEN → refactor while green. A
  test-only characterization task first records the live passing baseline, then proves assertion sensitivity with a
  temporary wrong expectation or controlled mutation that fails, restores the test, and reruns GREEN. Every task then
  receives a focused commit/report/review/fix loop.
- Keep branch coverage enabled and thresholds 86 hard floor, 86 advisory target, and 80 patch floor; never weaken or collapse them.
- Scripts-only branch coverage must reach at least 88.00% before `scripts` enters canonical source; combined coverage remains at least 86.00% on Python 3.10–3.13.
- Never mutate real source during mutation testing; snapshot overlays preserve current tracked/non-ignored untracked test semantics.
- Release automation changes workflows but never publishes during this PR. Never create/push a version tag, create a GitHub Release, delete credentials, merge the PR, or force-push.
- Preserve unrelated user changes; never reset, restore-discard, broad-clean, or overwrite them.
- Current SDD v6.2.0 forbids concurrent implementation subagents. Parallelize read-only audits and independent reviews only; implementation tasks are serialized.

---

## Plan family and frozen interfaces

| Plan | Owned task range | Frozen integration interface |
| --- | --- | --- |
| `2026-08-10-milestone-6-quality-gates-plan.md` | Tasks 1–2, pause, then Tasks 3–10 | `make type-check` covers `vntyper/ docker/app/ scripts/`; final coverage source is `vntyper`, `docker/app`, `scripts`. |
| `2026-08-10-milestone-6-mutation-safety-plan.md` | Tasks 1–7 | `MutationWorkspace(real_root, sweep_root, head, overlay_operations, baseline_manifest, baseline_status, baseline_digests)`; post-overlay baseline-relative restore; explicit disposable subprocess root. |
| `2026-08-10-milestone-6-harness-matrix-plan.md` | Tasks 1–6 | Literal `fastq_bam_processing/output_unmapped.bam`; disjoint success-only artifact declarations; key-presence archive Boolean. |
| `2026-08-10-milestone-6-release-and-naming-plan.md` | Tasks 1–11 | `scripts/release_policy.py`; existing `sha-<7>` plus full revision label; ten exact check names; digest feeds five aliases. |
| `2026-08-10-milestone-6-exception-policy-plan.md` | Tasks 0–12 | `scripts/ble001_policy.py` plus JSON inventory; normalized path/symbol/category; reviewed/actual Ruff version attribution. |

The bounded plan is the only valid `PLAN_FILE` for its implementer and reviewer briefs. Never run
`task-brief` for a master summary and ask an implementer to infer a track plan.

## SDD workspaces and execution schedule

Start every execution session in the authorized feature worktree, verify the branch, then use the installed current
skill directory. All relative plan paths and SDD workspace metadata are rooted here, never in the `main` checkout:

```bash
cd /home/bernt-popp/development/VNtyper/.worktrees/ci-milestone-6-gates-harness-release
test "$(git branch --show-current)" = ci/milestone-6-gates-harness-release
M6_SDD_DIR=/home/bernt-popp/.codex/plugins/cache/openai-curated-remote/superpowers/6.2.0/skills/subagent-driven-development
```

For each bounded plan, run:

```bash
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

Use `apply_patch` to create `progress.md` in each printed workspace. Its first line is
`# SDD ledger — plan: ` followed by that command's exact plan path; for example:

```text
# SDD ledger — plan: docs/superpowers/plans/2026-08-10-milestone-6-quality-gates-plan.md
```

Execute without a user pause:

```text
1. quality-gates Tasks 1–2         scripts-only measurement target and #204 mypy foundation
2. mutation-safety Tasks 1–7       #208
3. harness-matrix Tasks 1–6        #71 and #226 evidence
4. release-and-naming Tasks 1–11   #214, #218, #220
5. quality-gates Tasks 3–10        #211 against the final scripts tree
6. exception-policy Tasks 0–12     #219 against the final Python/test tree
7. master Tasks 1–5 below          integration, reviews, PR, post-PR terminal checks
```

Pausing quality after Task 2 is deliberate: mutation and release add script code that must be in
the 88.00% measurement. Its ledger persists completed Tasks 1–2; resume at Task 3 and never
redispatch them.

For every bounded task, follow current SDD exactly:

1. Capture `BASE=$(git rev-parse HEAD)` before dispatch.
2. Run `task-brief BOUNDED_PLAN TASK_NUMBER`; pass only that brief, frozen prior interfaces, and a
   report path to a fresh implementer with an explicit model.
3. Require tests, focused commit(s), self-review, and report.
4. Run `review-package BOUNDED_PLAN BASE HEAD`; give brief/report/package to a fresh reviewer with
   an explicit model and require both specification and quality verdicts.
5. Execute up to five fix rounds with scoped re-review; rounds 1–3 resume the implementer, rounds
   4–5 use a fresh stronger model. The controller never edits implementation fixes.
6. Append exact commits/verdicts/fix rounds/completion to that plan's ledger before advancing.

Model selection is explicit on every dispatch: fast model only for complete single-file mechanical
briefs, standard model for multi-file integration, most capable model for architecture/release
judgment and final whole-branch review; reviewers use at least the standard tier.

### Task 1: Verify and commit the reviewed specification corrections

**Files:**
- Modify: the six existing `docs/superpowers/specs/2026-08-10-milestone-6-*-design.md` files when review corrections apply.
- Existing: `mkdocs.yml` excludes `superpowers/`; initial specification commit is `ce588b2`.

**Interfaces:**
- Consumes: live issues, current code/workflows, reconstruction evidence, three exact Opus review rounds.
- Produces: a focused follow-up commit containing only adversarial-review corrections to the frozen contracts.

- [ ] **Step 1: Verify six specs and exclusion exist**

Run: `find docs/superpowers/specs -maxdepth 1 -name '2026-08-10-milestone-6-*.md' -print | sort && rg -n '^  superpowers/$' mkdocs.yml`

Expected: exactly six spec paths plus the exclusion.

- [ ] **Step 2: Verify no actual placeholder prose remains**

Run: `rg -n -i '\b(TBD|TODO|FIXME|implement later)\b' docs/superpowers/specs`

Expected: exit 1 and no output.

- [ ] **Step 3: Verify final scoped Opus verdict**

Inspect the captured JSON result from exact `claude-opus-5 --effort max --permission-mode plan`.

Expected: canonical model `claude-opus-5`, no unresolved Critical/High, verdict `ADVANCE`.

- [ ] **Step 4: Run strict docs build**

Run: `make docs-build`

Expected: exit 0; no `docs/superpowers/` page exists below `site/`.

- [ ] **Step 5: Verify the initial commit and commit only reviewed corrections**

Run:

```bash
git show --stat --oneline ce588b2
git diff --check -- docs/superpowers/specs mkdocs.yml
git add docs/superpowers/specs
git diff --cached --quiet || git commit -m "docs(ci): harden milestone six contracts"
```

Expected: `ce588b2` is the intentional initial specification commit; if review changed a specification, one focused
`docs(ci): harden milestone six contracts` correction commit exists. `mkdocs.yml` is not recommitted when unchanged.

### Task 2: Verify, adversarially review, and commit all plans

**Files:**
- Create: this master plan, five bounded track plans, and
  `2026-08-10-milestone-6-EXECUTION-PROMPT.md` (seven files total).

**Interfaces:**
- Consumes: Task 1 specs.
- Produces: complete execution package and requirements-to-task matrix below.

- [ ] **Step 1: Verify all seven plan-package files**

Run: `find docs/superpowers/plans -maxdepth 1 -name '2026-08-10-milestone-6-*.md' -print | sort`

Expected: master, five bounded plans, and execution prompt.

- [ ] **Step 2: Verify mandatory plan headers and checkbox actions**

Run: `for f in docs/superpowers/plans/*-plan.md; do rg -q 'REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development' "$f" && rg -q '^- \[ \]' "$f" || exit 1; done`

Expected: exit 0.

- [ ] **Step 3: Verify frozen cross-plan literals**

Run: `rg -n 'output_unmapped\.bam|MutationWorkspace\(|sha-<7|88\.00%|ble001_policy' docs/superpowers/plans`

Expected: each literal appears in the relevant bounded plan, master, and execution prompt without a conflicting form.

- [ ] **Step 4: Run exact Opus 5 maximum-effort plan review**

Give Opus the seven plan files, six specs, live issues, current code/workflows, current skills, and
the prior spec review findings. Require exact evidence and Critical/High/Medium/Low classification.

Expected: canonical model `claude-opus-5`, no unresolved Critical/High, verdict `ADVANCE` after one scoped correction review if required.

- [ ] **Step 5: Commit all seven planning artifacts**

Run:

```bash
git add docs/superpowers/plans
git commit -m "docs(ci): plan milestone six execution"
```

Expected: one Conventional Commit; worktree clean before implementation begins.

### Task 3: Run fresh whole-tree verification after all bounded ledgers complete

**Files:** no production edit is authorized by this task; fixes return to the owning bounded plan's SDD loop.

**Interfaces:**
- Consumes: every bounded ledger completion and reviewed track commit.
- Produces: final command/test/coverage/tier evidence for the exact HEAD.

- [ ] **Step 1: Inspect the complete diff and status**

Run: `git diff --stat origin/main...HEAD && git diff --check origin/main...HEAD && git status --short`

Expected: milestone-only diff, no whitespace errors, clean tree.

- [ ] **Step 2: Run formatting, lint, and type checks**

Run: `make format-check && make lint && make type-check`

Expected: exit 0 for all three.

- [ ] **Step 3: Run unit, coverage, and patch gates**

Run: `make test-unit && make test-unit-cov && make patch-coverage`

Expected: exit 0; record test count, precise combined/scripts branch coverage, and patch coverage.

- [ ] **Step 4: Run aggregate and fresh-environment CI**

Run: `make check-all && make ci-local`

Expected: exit 0; `ci-local-uv` rebuilds the environment successfully.

- [ ] **Step 5: Run Docker CI mirror**

Run: `make ci-local-docker`

Expected: exit 0 with built image/test artifacts; if daemon is genuinely unavailable, capture the exact failure and do not claim pass.

- [ ] **Step 6: Run mutation and integration tiers**

Run the mutation plan's real canary/round-trip commands, then run `VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1 make check-full` when archive/external tools are available.

Expected: real source byte-identical; generator round trip stable; available integration tiers pass; unavailable dependency is reported exactly.

### Task 4: Obtain whole-branch Superpowers and Opus reviews

**Files:** merge-base..HEAD review package; any fix returns to one owning bounded SDD ledger.

**Interfaces:**
- Consumes: Task 3 evidence and every spec/plan.
- Produces: reviewed final SHA with no unresolved Critical/High.

- [ ] **Step 1: Generate whole-branch review package**

Run: `"$M6_SDD_DIR/scripts/review-package" docs/superpowers/plans/2026-08-10-milestone-6-integration-plan.md "$(git merge-base origin/main HEAD)" "$(git rev-parse HEAD)"`

Expected: printed package path contains full commits/stat/diff.

- [ ] **Step 2: Dispatch final Superpowers reviewer with explicit most-capable model**

Expected: specification and quality verdicts cover the whole branch plus parked/deferred ledger lines.

- [ ] **Step 3: Run exact Opus 5 integrated review**

Run Claude Code with exact `claude-opus-5`, `--effort max`, read-only plan mode over actual diff/specs/plans/tests/workflows/evidence.

Expected: canonical model confirmation and classified findings.

- [ ] **Step 4: Route one coordinated fix wave through SDD**

If findings exist, select the owning bounded plan, capture `FIX_BASE`, resume/fresh-dispatch its
implementer with the complete accepted findings, require report/tests/commit, then generate a scoped
review package and dispatch one re-review.

Expected: accepted findings addressed; technical rulings recorded for rejected feedback; no load-bearing residual.

- [ ] **Step 5: Rerun Task 3 and perform the single scoped re-review**

Run one scoped Superpowers re-review and one scoped Opus re-review over the accepted fix wave. There is no second final
fix wave. Technically adjudicate residual feedback under `receiving-code-review`; an unresolved load-bearing finding
stops the program. Expected: the reviewed SHA exactly equals the SHA whose fresh gates will be reported.

### Task 5: Push, open, monitor, and revalidate the pull request

**Files:** PR metadata; any post-PR repository edit returns to bounded SDD and Tasks 3–4.

**Interfaces:**
- Consumes: clean reviewed SHA and exact evidence.
- Produces: pushed branch, PR URL, terminal required checks, reviewed final SHA.

- [ ] **Step 1: Verify branch history and clean state**

Run: `git status --short && git log --oneline --decorate origin/main..HEAD`

Expected: clean tree and coherent Conventional Commits.

- [ ] **Step 2: Push without force and open PR using `github:yeet`**

Run: `git push -u origin ci/milestone-6-gates-harness-release`, then use the PR skill against `main`.

Expected: PR URL and body containing nine issues, traceability, exact evidence, risk, review dispositions, and accurate `Closes` references.

- [ ] **Step 3: Monitor checks to terminal state**

Use GitHub checks monitoring; diagnose failures with `github:gh-fix-ci` and `superpowers:systematic-debugging`.

Expected: every required check terminal; no failure is dismissed as flaky without evidence.

- [ ] **Step 4: Route every post-PR fix back through SDD**

For any CI/review fix, identify its owning bounded plan, append `PR fix round` plus the next integer to that ledger,
capture base SHA, dispatch implementer with exact failure/comment evidence, require regression test,
commit/report, scoped task review, and affected gate rerun.

Expected: no controller-authored or unreviewed post-PR commit.

- [ ] **Step 5: Revalidate every post-PR commit within the exhausted final-review budget**

Rerun Task 3, regenerate the whole-branch package, and obtain one renewed Superpowers/Opus review of the post-PR-fix
SHA. Do not open another final fix wave; stop on a load-bearing finding. Push without force and monitor again.

Expected: final reported SHA is exactly the final locally/externally reviewed SHA.

- [ ] **Step 6: Finish without merging or deleting the Git worktree**

After terminal-success checks, delete only the five scratch workspaces printed by `sdd-workspace`
as current SDD requires, invoke `superpowers:finishing-a-development-branch`, choose the
PR/leave-worktree option, and retain the Git worktree for iteration.

Expected: PR remains open, branch/worktree remain available, no tag/release/merge occurs.

## Requirements-to-task traceability matrix

| Requirement | Bounded plan/task(s) | Objective evidence |
| --- | --- | --- |
| #71 AC1 Docker | Harness plan evidence table + Task 6 | Existing Docker pipeline node IDs rerun/cited. |
| #71 AC2 API | Harness plan evidence table + Task 6 | Existing API status/error/JSON node IDs rerun/cited. |
| #71 AC3 adVNTR | Harness plan evidence table + Task 6 | Existing adVNTR output/max-coverage cases rerun when tier available. |
| #71 AC4 subcommands | Harness plan evidence table + Task 6 | Existing parser/handler cases cited with explicit narrowing disposition. |
| #71 AC5 parameter artifacts | Harness Tasks 1–5 | Positive literal existence then negative delete/precedence plus archive true/false. |
| #71 AC6 alternate FASTQ/no SHARK | Harness Task 6 | Paired b178 mates, no SHARK option, successful integration declaration. |
| #204 annotation/mypy/parity | Quality Task 2 | Widened no-cache mypy and executable Makefile contract test. |
| #208 interruption safety | Mutation Tasks 1–7 | Overlay/rename/path tests, provenance/canary, atomic output, byte-identity real run. |
| #211 scripts coverage | Quality Tasks 3–10 | Both zero modules non-zero; scripts >=88.00%; combined >=86.00% on 3.10–3.13; patch >=80. |
| #214 GHCR aliases/docs | Release Tasks 1–7, 9–11 | Ten check names, short-SHA/full revision, digest aliases, global concurrency, GHCR docs. |
| #218 PyPI OIDC | Release Tasks 5, 8, 10–11 | Existing-tag validation, unprivileged build, environment `pypi`, publish-only OIDC, pinned action. |
| #219 BLE001 policy | Exception Tasks 0–12 | Final dual inventory, A/B/C completeness, version-attributed drift, category-C behavior tests. |
| #220 generation naming | Release Task 9 | Exactly nine authorized edits; protected versions/identifiers/history unchanged. |
| #226 fixture citation | Harness Task 6 evidence | Existing commits/function/unit/reference-contract tests cited and rerun, no duplicate code. |
| Python 3.10–3.13 | All bounded plans; master Task 3 | Syntax/type constraints and GitHub unit matrix. |
| Branch/patch coverage | Quality Tasks 8–9; master Task 3 | Branch source invariant, precise combined result, merge-base diff-cover. |
| Fresh CI environment | Release Task 11; master Task 3 | `make ci-local` including `ci-local-uv`. |
| Docker workflow behavior | Release Tasks 4–11; master Task 3 | actionlint plus `make ci-local-docker`. |
| Partial failure/rerun/security | Mutation, release, exception plans | Negative tests, promotion mutual exclusion plus explicit canceled-pending retry, atomic outputs, job-scoped permissions. |
| Final reviews and PR | Master Tasks 4–5 | Reviewed final SHA, PR URL, terminal required checks, no unresolved Critical/High. |

## Program completion

The program is complete only when all five bounded ledgers name every task complete, Tasks 3–5
have been rerun for the final SHA, the branch is clean/pushed, the PR URL exists, required checks
are terminal and successful, and no unresolved Critical/High remains. Do not merge or release.
