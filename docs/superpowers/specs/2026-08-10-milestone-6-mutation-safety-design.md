# Milestone 6 mutation safety design

**Status:** approved design for
[#208](https://github.com/hassansaei/VNtyper/issues/208).

**Baseline:** `main` at `ebb15b26631242a3295607e4eda4c68f688cd9a2` (`v2.0.10`).

## Purpose

Make a mutation sweep incapable of leaving a live mutant in the operator's real VNtyper
source tree, including after SIGHUP, SIGQUIT, SIGKILL, interpreter crash or host loss.
The sweep will mutate only a disposable detached worktree created from the current HEAD
and overlaid with validated current working-state paths. Tests will execute with that worktree as their actual CWD and import root, while
the generated Markdown and JSON outputs remain rooted in the real repository.

This implements the worktree direction selected in the
[#208 evidence and recommendation](https://github.com/hassansaei/VNtyper/issues/208#issuecomment-5214049578).
Signal handlers alone are not the fix: they cannot cover SIGKILL or host loss.

## Source acceptance criteria

The [#208 issue body](https://github.com/hassansaei/VNtyper/issues/208) defines the
failure and the acceptable direction. This design turns it into the following exact,
linked criteria:

1. **AC208-1 — no real-source mutation.** No mutation path, bytecode cleanup or pytest
   child may write under the real repository's `vntyper/` tree.
2. **AC208-2 — disposable working snapshot.** A sweep uses a unique
   `git worktree add --detach <path> <head>` checkout, verifies its HEAD, overlays the
   validated current tracked and non-ignored untracked state, records both the overlay
   operations and the authoritative post-overlay baseline manifest,
   and mutates only there.
3. **AC208-3 — actual execution isolation.** Every baseline, canary, scoped and full-tier
   pytest subprocess uses the disposable worktree as `cwd`; import provenance must
   resolve `vntyper` and the selected target below the disposable root, never through
   the editable install of the real tree.
4. **AC208-4 — separate output root.** `--output`, `--results-json` and `--render-only`
   remain interpreted relative to the original real repository and are never redirected
   into the disposable worktree.
5. **AC208-5 — killed-mutant canary.** Before the full sweep, one committed, known-killed
   mutant is run in the disposable tree. If it is not loaded and killed, the run refuses
   before publishing any measurement.
6. **AC208-6 — crash-safe real tree.** Normal completion, Python exceptions, handled
   signals, failed cleanup and uncatchable termination all leave real production source
   byte-identical because it was never opened for mutation.
7. **AC208-7 — best-effort disposal.** Normal and handled-failure paths remove the
   temporary worktree and Git registration. Failure to remove it is visible and
   non-zero but never triggers destructive cleanup of the real tree.
8. **AC208-8 — measurement semantics preserved.** Green-baseline refusal, scoped kill,
   full-tier survivor escalation, bytecode defenses, equivalent-mutant reporting and
   advisory-only scoring retain their current meaning.

## In scope

- A focused, fully annotated `scripts/mutation_workspace.py` module owning Git worktree
  lifecycle, path confinement and import-provenance checks.
- Refactoring `scripts/mutation_test.py` to distinguish the immutable real repository
  root from the disposable sweep root.
- Passing the sweep root explicitly to pytest, baseline, bytecode and mutation
  operations rather than consulting one global `REPO_ROOT`.
- Preserving real-root cleanliness checks for selected source paths and generated
  outputs before creating a worktree.
- A known-killed canary based on the committed scoring target.
- Best-effort cleanup on success, failure, SIGINT, SIGTERM, SIGHUP and SIGQUIT.
- Pure unit tests with mocked Git/process boundaries plus one objective real-worktree
  canary test.
- Updating mutation harness documentation and Makefile warnings so they no longer claim
  production source is rewritten in place.

## Out of scope

- Changing mutation operators, target modules, equivalent-mutant classifications or
  the score threshold. Mutation testing remains advisory.
- Using a durable backup, recovery sentinel or supervisor around real-source writes.
- Installing a second environment inside the disposable worktree.
- Copying gitignored integration data into the worktree; mutation targets run unit
  tests only.
- Making mutation testing a PR gate or shortening its 15–30 minute runtime.
- Automatically removing unrelated worktrees or running broad `git worktree prune`.
- Recovering an orphaned disposable worktree after host loss as a condition of source
  safety. Orphans are bounded cleanup debt, not live-source corruption.
- Recommending `git checkout --`, `git restore`, `git reset --hard` or `git clean`.

## Architecture and boundaries

### 1. Real repository and disposable sweep repository are distinct capabilities

Rename the current conceptual root to `REAL_REPO_ROOT`. It is used only for:

- resolving the captured HEAD commit;
- checking selected real targets and requested output files for uncommitted changes;
- creating/removing the disposable worktree through Git; and
- writing requested report/result outputs after a valid measurement.

The disposable tree is a snapshot of the current invocation tree, not merely committed
HEAD. This preserves the existing developer loop in which a newly written, uncommitted
test is run against a mutant before the test is committed. The selected mutation target
and persistent outputs still must pass the existing dirty guard.

`REAL_REPO_ROOT` must never be passed to `generate_mutants`, `sweep_module`, mutation-target
`write_text`, target bytecode cleanup or pytest `cwd`.

The disposable `sweep_root` is used only for:

- reading and mutating target source;
- clearing `sweep_root / "vntyper"` bytecode caches;
- import provenance probes;
- the unmutated baseline;
- the killed-mutant canary;
- scoped and full-tier pytest runs; and
- the end-of-sweep post-overlay-baseline restore check inside the disposable worktree.

Both roots are resolved absolute paths. Startup refuses if they are equal, if either is
an ancestor of the other, or if a target resolved from `sweep_root` escapes that root.
The temporary parent is the system temporary directory, not the real checkout.

### 2. Extract worktree lifecycle from the oversized runner

`scripts/mutation_test.py` is already well beyond the repository's approximate 650-line
limit. Do not add Git lifecycle and path-policy code to it. Add
`scripts/mutation_workspace.py` with a small interface such as:

```python
@dataclass(frozen=True)
class MutationWorkspace:
    real_root: Path
    sweep_root: Path
    head: str

    def target_path(self, relative: str) -> Path: ...


@contextmanager
def detached_head_workspace(real_root: Path) -> Iterator[MutationWorkspace]: ...
```

Repository convention forbids custom exception classes. Git, path or provenance
failures log an exact message and raise `RuntimeError` with the same text; the mutation
CLI catches the stage-boundary error and returns 1.

Lifecycle:

1. Resolve `git -C <real_root> rev-parse --verify HEAD^{commit}`.
2. Create a unique system temporary directory.
3. Run `git -C <real_root> worktree add --detach <sweep_root> <head>`.
4. Verify the disposable HEAD exactly equals the captured commit.
5. Snapshot `git status --porcelain=v1 -z --untracked-files=all` in the real checkout.
   Overlay every modified tracked file and non-ignored untracked file into the
   disposable tree, and mirror tracked deletions there. Parse rename records as the
   two NUL-delimited paths emitted by porcelain v1: copy the new path and remove the
   old path. Reject absolute paths, `..`,
   paths below the disposable-worktree parent, and symlinks whose resolved target
   escapes the real checkout. Never copy `.git`, ignored caches, mutation outputs, or
   a registered worktree whose root is a strict descendant of `real_root`. A registered
   worktree that is an ancestor of `real_root` is the containing checkout and does not
   exclude ordinary paths inside `real_root`.
6. Record the deterministic `overlay_changes` applied from the real checkout. Then query
   status inside the disposable tree and record a separate deterministic
   `baseline_manifest` plus content/symlink/deletion digests. `baseline_manifest`, not
   `overlay_changes`, is authoritative for later restoration: a staged real change whose
   worktree bytes were reverted to HEAD still requires a copy operation, but that copy
   correctly produces no disposable-tree status entry. Verify selected target bytes in
   the disposable tree equal their guarded real-tree bytes.
7. Yield the workspace.
8. In `finally`, run `git -C <real_root> worktree remove --force <sweep_root>` and remove
   a remaining temporary directory only when it is the exact validated temporary path.

No unresolved glob, `$HOME`, repository root or broad directory may be a cleanup target.

### 3. Thread the sweep root through existing interfaces

The following current symbols must stop deriving mutation execution from global
`REPO_ROOT`:

- `run_pytest(test_paths, timeout, parallel=False, *, repo_root)`
- `run_tests(..., *, repo_root)`
- `baseline_refusal(targets, timeout, *, repo_root)`
- `sweep_module(path_str, scoped_tests, timeout, verbose, *, repo_root)`

Within `main`:

- `_refuse_if_dirty` continues to use `REAL_REPO_ROOT` and protects selected real target
  paths plus real output paths.
- `detached_head_workspace(REAL_REPO_ROOT)` provides `workspace.sweep_root`, the applied
  `overlay_changes`, and the authoritative post-overlay `baseline_manifest`.
- every baseline/canary/sweep call receives that root explicitly;
- `write_outputs` receives already-resolved real output paths and runs only after a
  complete valid sweep; and
- `--render-only` creates no worktree because it runs no tests and mutates no source.

Relative output arguments are resolved before the worktree is created and before any
CWD-sensitive subprocess. They retain today's meaning relative to the invocation's real
repository CWD.

### 4. CWD isolation and import provenance are separate invariants

Setting only `PYTHONPATH` is insufficient, and setting only `cwd` can silently lose to
an editable installation or a preloaded module. `run_pytest` must do both:

- `cwd=sweep_root` for every child;
- prepend `sweep_root` to a sanitized child `PYTHONPATH` while removing the resolved
  real repository root;
- retain `python -B` and `PYTHONDONTWRITEBYTECODE=1`; and
- never import target production modules into the long-lived parent process.

Before the baseline, run a child-process provenance probe from `sweep_root` that imports
`vntyper` and the selected target module, resolves their `__file__` paths, and refuses
unless both are beneath `sweep_root`. The probe prints the resolved paths on success and
failure so an editable-install leak is diagnosable.

The existing deletion of every `__pycache__` under `vntyper/` remains, but its root is
the disposable worktree. The real tree's caches are never touched.

The overlay means baseline, canary, scoped tests, and full-tier escalation see the same
tracked and non-ignored untracked tests and dependent source that the current in-place
runner sees. A unit test written immediately before `make mutation` is therefore part of
the measurement. Ignored/generated files remain environment inputs and are not copied;
tests that require them must address them through existing absolute configuration or be
outside the unit mutation target.

### 5. A known-killed mutant is the execution canary

Path provenance proves where Python says it imported from; the canary proves that the
actual mutated bytes affect the exact pytest decision used by the score.

Use the committed scoring mutant identified by:

```text
path:        vntyper/scripts/scoring.py
function:    split_depth_and_calculate_frame_score
line/token:  current line 74, '/' -> '*'
tests:       TARGETS["vntyper/scripts/scoring.py"]
```

The committed `docs/development/mutation-results.json` records scoring as 20 killed and
0 surviving mutants at the baseline. The canary lookup must fail closed if that exact
mutant no longer exists; line or source drift requires choosing and documenting a new
known-killed mutant rather than silently taking "the first" mutant.

Canary sequence:

1. The unmutated scoped and full-tier baseline passes in `sweep_root`.
2. Write only the selected canary mutant to its guarded worktree target.
3. Clear bytecode below `sweep_root / "vntyper"`.
4. Run its scoped tests with `cwd=sweep_root`.
5. Require a non-zero pytest result, meaning the canary was killed.
6. Restore the disposable target in `finally` and prove the disposable worktree status
   exactly equals the recorded post-overlay `baseline_manifest` before starting the full
   sweep. The baseline manifest may be empty even when `overlay_changes` is non-empty—for
   example, when staged content is reverted in the real worktree back to HEAD bytes.

If the canary passes, times out for an unrelated reason, cannot be located, or cannot be
restored, return non-zero and write no result files. A timeout is a kill during an
ordinary measured mutant, but it is not acceptable provenance proof for the canary; the
canary must fail by the expected scoped assertion signal rather than infrastructure
failure. Capture and classify the pytest result accordingly.

### 6. Output ownership remains with the real root

The default Makefile destinations remain:

- `docs/development/mutation-testing.md`
- `docs/development/mutation-results.json`

They are protected by the existing dirty-output preflight. They are not copied into the
worktree and are not part of its cleanup. The sweep accumulates results in memory and
does not call `write_outputs` until the baseline, provenance probe, canary and every
selected mutant have completed.

Each output should be staged in its destination directory and installed with
`os.replace`, so an interrupted single-file write leaves either the previous complete
file or the new complete file. The two files are not claimed to be a filesystem
transaction; a failure between their replacements returns non-zero and is repaired by a
retry or `make mutation-render` from the complete JSON result.

### 7. Cleanup is best effort, while real-source safety is unconditional

SIGINT, SIGTERM, SIGHUP and SIGQUIT should enter the same Python unwind path so the
worktree context manager attempts cleanup. This improves hygiene but is not the safety
boundary.

SIGKILL, interpreter crash and host loss may leave a registered temporary worktree that
contains a mutant. They still leave the real source byte-identical because no mutation
write targeted it. The next run creates a new unique worktree and does not reuse the
orphan.

If normal cleanup fails:

- name the exact disposable path and Git error;
- return non-zero;
- do not delete or reset anything in the real tree;
- do not run a broad prune automatically; and
- retain already completed result files, while saying that the measurement completed
  but cleanup did not.

## Inputs, outputs and invariants

### Inputs

- `REAL_REPO_ROOT`, resolved from the script location.
- The exact committed HEAD object ID.
- Existing CLI arguments: `--module`, `--timeout`, `--verbose`, `--output`,
  `--results-json`, and `--render-only`.
- `TARGETS` and their scoped unit test paths.
- The fixed canary identity.
- The active Python environment; no new dependency environment is built.

### Outputs

- Console progress and diagnostics.
- Optional real-root Markdown report.
- Optional real-root JSON measurement.
- Exit 0 only for a valid completed measurement and successful cleanup, or a successful
  render-only operation.
- Exit 1 for preflight, Git, provenance, canary, measurement, output or cleanup failure.

### Hard invariants

- A path opened for mutation resolves below `sweep_root` and outside
  `REAL_REPO_ROOT`.
- Every pytest child has `cwd == sweep_root`.
- Imported `vntyper` and target paths resolve below `sweep_root`.
- The canary is killed before any report can be published.
- The real source digest is unchanged across all outcomes, including cleanup failure.
- Bytecode from another revision cannot be loaded.
- A red unmutated baseline cannot score mutants as killed.
- A scoped survivor is still escalated to the full unit tier.
- Output paths never become worktree cleanup targets.
- Render-only never creates or mutates a worktree.

## Success, partial failure, recovery, retry and idempotency

### Success

A sweep is successful only when:

1. selected real targets and real outputs pass the dirty-path guard;
2. the detached worktree is at the captured HEAD, all validated overlay operations were
   applied, and its status equals the separately recorded post-overlay baseline manifest;
3. provenance resolves into the worktree;
4. the unmutated baseline is green;
5. the known-killed canary is killed and restored;
6. all selected mutants are measured with existing scoped/full semantics;
7. requested outputs are installed successfully;
8. the disposable tree is removed; and
9. the real production-source digest equals its startup digest.

### Partial failure

| Failure point | Result |
| --- | --- |
| Dirty or indeterminate real-root preflight | Refuse before worktree creation; write nothing. |
| Worktree add/verification failure | Return 1; write no measurement; remove only the exact temporary path if safe. |
| Import provenance leak | Return 1 before baseline/canary; write no measurement. |
| Red baseline | Return 1; no target is mutated and prior outputs remain. |
| Canary missing, surviving or infrastructurally failing | Return 1; restore disposable target; write no measurement. |
| Exception or handled signal during sweep | Restore is attempted in the disposable tree, cleanup is attempted, real source remains untouched. |
| SIGKILL/crash/host loss | Disposable mutant may remain orphaned; real source and prior outputs remain untouched unless output installation had already begun. |
| Output replacement failure | Return 1; each destination is individually old-complete or new-complete. |
| Worktree cleanup failure after output | Return 1 and name the orphan; completed measurement files remain available. |

### Recovery

No source restoration command is ever required. For an orphan, the diagnostic names the
exact temporary worktree and leaves removal to an explicit operator action. Existing
guards continue to prohibit advice that can discard real uncommitted work.

If only the Markdown output is stale after a partial output failure and the JSON result
is complete, `make mutation-render` is the supported recovery. Otherwise rerun the
selected sweep.

### Retry and idempotency

- Every retry uses a new unique detached worktree from the same current HEAD.
- A previous orphan is never reused.
- Mutation order and classification remain deterministic for the same commit and
  environment; elapsed time may differ.
- Successful output installation replaces the same destinations rather than appending.
- `--module` narrows the baseline, canary applicability and sweep exactly as today. If
  the selected module excludes the global scoring canary, run the canary as a separate
  provenance-only preflight; it does not enter the selected module's reported totals.
- Render-only remains idempotent for the same JSON input, apart from no measurement
  timestamp being introduced.

## Compatibility and environment constraints

- Python support remains 3.10 through 3.13. Use `tempfile`, `contextlib`, `Path` and
  `subprocess`; do not require Python 3.11 APIs.
- Git with worktree support is a precondition. Missing or failing Git is an explicit
  refusal, not permission to fall back to real-source mutation.
- The active venv/conda interpreter and installed dependencies are reused. The worktree
  contains captured-HEAD source/tests plus the validated working-state overlay; unit tests require
  no downloaded reference data.
- Pytest still runs from a repository root—now the disposable one—so relative
  `tests/test_data_config.json` collection remains valid.
- `python -B`, `PYTHONDONTWRITEBYTECODE=1`, per-mutant cache deletion,
  `-p no:cacheprovider`, unit marker selection and timeout behavior are preserved.
- New unit test files carry `pytestmark = pytest.mark.unit`.
- `scripts/` must be included in mypy before or in the same milestone dependency wave;
  new workspace code is fully annotated.

## Security and permissions

- The design requires no network, credentials, GitHub token or permission change.
- All Git operations are local and target the exact repository, captured commit and
  validated temporary path.
- Commands use argument arrays, not `shell=True`.
- The worktree path is generated by `tempfile`; it is not supplied through an
  environment variable, glob or unresolved command substitution.
- Target resolution rejects absolute paths, `..`, symlink escapes and paths outside the
  disposable root.
- Cleanup never targets `$HOME`, `~`, `/`, the real repository root or an unvalidated
  directory.
- The worktree contains source and tests only; real generated outputs are not exposed as
  mutation targets.
- A malicious or broken editable-install mapping is detected by the provenance probe
  and killed-mutant canary before publication.

## Observability and diagnostic output

At minimum, the harness prints:

- captured HEAD object ID;
- disposable worktree path;
- resolved `vntyper` and selected target import paths;
- baseline check labels and outcomes;
- canary identity and `killed (scoped)` outcome;
- existing per-mutant progress and score;
- real output destinations;
- cleanup success, or the exact orphan path and Git error; and
- final confirmation that the real source digest did not change.

Diagnostics distinguish "measurement failed" from "measurement completed, cleanup
failed". They never tell the operator to discard local work.

## Test strategy and objective acceptance checks

### Pure/unit tests

Add `tests/unit/test_mutation_workspace.py` for the extracted lifecycle and extend
`tests/unit/test_mutation_test.py` for orchestration:

| Test | Acceptance evidence |
| --- | --- |
| Worktree command construction uses `add --detach` and the captured full HEAD | AC208-2 |
| Added worktree with wrong HEAD is refused; post-overlay status must equal its separately captured baseline manifest rather than the overlay-operation list or an assumed empty status | AC208-2, AC208-8 |
| Porcelain rename records copy the new path and remove the old path in the snapshot | AC208-8 |
| Target resolver rejects absolute, parent, symlink-escape and real-root paths | AC208-1, AC208-6 |
| Every baseline/canary/scoped/full subprocess receives disposable `cwd` | AC208-3 |
| Child environment removes real-root import precedence and prepends sweep root | AC208-3 |
| Provenance probe accepts worktree paths and rejects editable real-root paths | AC208-3 |
| Output arguments resolve against the real invocation root before workspace entry | AC208-4 |
| Canary lookup requires the exact scoring mutant | AC208-5 |
| Canary pass, timeout, collection error or missing identity refuses publication | AC208-5 |
| Canary restoration returns the disposable worktree to its exact post-overlay baseline-manifest state | AC208-5, AC208-8 |
| Real-source digest is unchanged after success, raised exception and cleanup failure | AC208-1, AC208-6 |
| Cleanup is attempted on success, exception and handled signals | AC208-7 |
| Cleanup failure returns non-zero and never invokes broad/destructive Git commands | AC208-7 |
| Red baseline still aborts before canary, sweep or output | AC208-8 |
| Scoped survivors still escalate to the disposable full tier | AC208-8 |
| Bytecode caches are cleared only under the disposable root | AC208-1, AC208-8 |
| Render-only creates no workspace and writes only its declared real output | AC208-4, AC208-8 |

All Git and pytest process calls in unit tests are mocked. Assertions inspect complete
argument arrays, CWDs, environments, ordering and written paths rather than merely
checking that a call occurred.

### Objective real-worktree canary

Add a bounded test or documented verification mode that actually:

1. creates a detached worktree from HEAD;
2. proves the real and disposable target digests initially match;
3. applies the exact scoring canary only in the disposable tree;
4. runs its scoped tests from that worktree and observes the expected assertion failure;
5. proves the real target digest never changed;
6. restores/removes the worktree; and
7. proves `git worktree list --porcelain` no longer names it.

This is the acceptance canary for the silent failure mode; mocked unit tests alone
cannot prove editable-install isolation.

Required verification:

```bash
pytest -m unit tests/unit/test_mutation_workspace.py tests/unit/test_mutation_test.py
make test-unit-cov
make patch-coverage
python scripts/mutation_test.py --module scoring.py --verbose \
  --output /tmp/vntyper-mutation-canary.md \
  --results-json /tmp/vntyper-mutation-canary.json
make check-all
```

After the scoped run, compare the real scoring source digest captured before the command
and verify no new worktree registration remains. The temporary `/tmp` outputs avoid
rewriting committed evidence during acceptance.

## Rollout

1. Land the scripts mypy gate prerequisite or include it in the same dependency wave.
2. Add failing unit tests for lifecycle commands, root separation, CWD/import
   provenance, canary behavior and real-source immutability.
3. Extract `scripts/mutation_workspace.py`; do not grow `mutation_test.py` with the Git
   lifecycle.
4. Thread explicit `repo_root` parameters through baseline, pytest and sweep functions.
5. Add the provenance probe and exact known-killed canary.
6. Move mutation and bytecode writes exclusively to the disposable root.
7. Keep outputs rooted in the real repository and make individual replacements atomic.
8. Add best-effort signal cleanup and precise orphan diagnostics.
9. Update the mutation module narrative, generated documentation text and Makefile
   warning to describe worktree isolation rather than in-place real-source mutation.
10. Run focused unit tests, scripts-inclusive coverage, the real scoring sweep and the
    full pre-PR gate.
11. Re-render committed mutation documentation only from a complete accepted
    measurement; confirm the renderer round trip is byte-identical.

The rollout does not include an in-place fallback. If worktree creation or provenance
verification fails, the safe behavior is refusal.

## Unresolved questions

None. The isolation mechanism, root ownership, canary, signal posture, output location,
cleanup policy and failure semantics are fixed by this design.
