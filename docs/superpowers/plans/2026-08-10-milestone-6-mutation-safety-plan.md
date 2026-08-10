# Milestone 6 Mutation Safety Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close #208 by mutating only a disposable detached worktree that faithfully overlays the current non-ignored working state, proves pytest imports the mutant with a known-killed canary, writes reports atomically to the real root, and leaves real production source untouched under every outcome.

**Architecture:** A new fully annotated `scripts/mutation_workspace.py` owns porcelain parsing, confined overlay, detached-worktree lifecycle, restoration baselines and import provenance. `scripts/mutation_test.py` receives an explicit disposable root for every baseline/test/mutation action; a small `scripts/mutation_output.py` owns atomic output replacement so the already oversized runner does not absorb more I/O logic. The detached HEAD is overlaid with modified tracked files, non-ignored untracked files and tracked deletions so newly written uncommitted tests retain today's mutation semantics.

**Tech Stack:** Python 3.10–3.13, Git worktrees and porcelain v1 `-z`, `pathlib`, `tempfile`, `contextlib`, `subprocess`, pytest, Ruff, mypy, coverage.py and diff-cover.

## Global Constraints

- Follow `docs/superpowers/specs/2026-08-10-milestone-6-mutation-safety-design.md` exactly.
- Never mutate, clear bytecode below, or run pytest from `REAL_REPO_ROOT / "vntyper"`.
- Create a unique system-temp `git worktree add --detach <sweep_root> <captured-head>`; no in-place fallback exists.
- Overlay current tracked modifications, staged/unstaged tracked deletions and non-ignored untracked files before baseline tests; preserve an uncommitted test developer loop.
- Selected mutation targets and persistent outputs still pass the existing real-root dirty guard.
- Ordinary dirty non-target source and tests are overlaid, not refused; no whole-tree-clean precondition is permitted.
- Reject absolute, parent-traversal and symlink-escape overlay/target paths; never copy `.git`, ignored files, mutation outputs or another worktree.
- Every pytest/provenance child has `cwd=sweep_root`; sanitize `PYTHONPATH`, preserve `python -B` and `PYTHONDONTWRITEBYTECODE=1`, and clear caches only in the disposable tree.
- The exact canary is `vntyper/scripts/scoring.py`, current line 74, `/` to `*`, killed by its declared scoped tests.
- Reports remain rooted at the real invocation root and are individually installed with same-directory temporary files plus `os.replace`.
- SIGINT, SIGTERM, SIGHUP and SIGQUIT attempt cleanup; SIGKILL safety comes from root separation, not signal handling.
- Cleanup is best effort and exact-path only; never run broad prune or recommend destructive real-tree Git commands.
- Python support is 3.10–3.13; new scripts code is fully typed and must pass the scripts mypy gate.
- New unit files declare `pytestmark = pytest.mark.unit`; unit tests mock Git and pytest subprocesses.
- Preserve mutation operators, targets, equivalent classifications, baseline refusal, survivor escalation and advisory scoring.
- Treat each checkbox as one 2–5-minute edit, assertion addition, review, or command launch; record its result before starting the next checkbox.

---

## File map

| File | Responsibility |
| --- | --- |
| `scripts/mutation_workspace.py` | Porcelain parsing, safe overlay plan/application, detached worktree lifecycle, target confinement, baseline snapshot and provenance probe. |
| `scripts/mutation_output.py` | Atomic same-directory text installation and cleanup after failed replacement. |
| `scripts/mutation_test.py` | Mutation generation/scoring and top-level orchestration using explicit real/sweep roots. |
| `scripts/mutation_guard.py` | Existing real-root target/output dirty guard; wording updated only where it still claims real source is rewritten. |
| `tests/unit/test_mutation_workspace.py` | Pure/mock lifecycle, overlay, safety, provenance, cleanup and retry tests. |
| `tests/unit/test_mutation_output.py` | Atomic success/failure tests. |
| `tests/unit/test_mutation_test.py` | Root threading, canary, ordering, signals, render-only and real-source immutability tests. |
| `Makefile` | Remove in-place warning; document disposable worktree behavior. |
| `docs/development/mutation-testing.md` | Machine-generated explanation after renderer update. |
| `docs/development/mutation-results.json` | Measurement input/output; contents change only after an accepted real sweep. |

### Task 1: Parse and validate the current working-state overlay

**Files:**
- Create: `scripts/mutation_workspace.py`
- Create: `tests/unit/test_mutation_workspace.py`

**Interfaces:**
- Produces: `OverlayChange(action: Literal["copy", "delete"], path: str)`.
- Produces: `parse_porcelain_z(payload: bytes, real_root: Path) -> tuple[OverlayChange, ...]`; filesystem state resolves
  multiple porcelain records for the same path to one final copy/delete action.
- Produces: `_affected_paths_from_porcelain_z(payload: bytes) -> tuple[str, ...]`; pure token decoding/rename expansion.
- Produces: `confined_path(root: Path, relative: str, *, must_exist: bool) -> Path`.

- [ ] **Step 1: Create an importable scaffold and the marked RED porcelain test**

Create `scripts/mutation_workspace.py` with only the imports, `OverlayChange` dataclass, and public signatures whose
bodies raise `NotImplementedError("not implemented")`. This avoids a collection/import error while providing no
behavior. Then create the marked test module:

```python
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import mutation_workspace

pytestmark = pytest.mark.unit


def test_porcelain_z_maps_changes_renames_untracked_and_deletions(tmp_path: Path) -> None:
    payload = (
        b" M tests/unit/test_changed.py\0"
        b"?? tests/unit/test_new.py\0"
        b" D docs/removed.md\0"
        b"R  tests/unit/test_renamed.py\0tests/unit/test_old.py\0"
    )
    for relative in ("tests/unit/test_changed.py", "tests/unit/test_new.py", "tests/unit/test_renamed.py"):
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("current\n", encoding="utf-8")
    assert mutation_workspace.parse_porcelain_z(payload, tmp_path) == (
        mutation_workspace.OverlayChange("delete", "docs/removed.md"),
        mutation_workspace.OverlayChange("copy", "tests/unit/test_changed.py"),
        mutation_workspace.OverlayChange("copy", "tests/unit/test_new.py"),
        mutation_workspace.OverlayChange("delete", "tests/unit/test_old.py"),
        mutation_workspace.OverlayChange("copy", "tests/unit/test_renamed.py"),
    )
```

- [ ] **Step 2: Add RED confinement tests**

```python
@pytest.mark.parametrize("relative", ["", ".", "/absolute.py", "../escape.py", ".git/config", "a/.git/config"])
def test_confined_path_rejects_unsafe_names(tmp_path: Path, relative: str) -> None:
    with pytest.raises(ValueError, match="unsafe workspace path"):
        mutation_workspace.confined_path(tmp_path, relative, must_exist=False)


def test_confined_path_rejects_a_symlink_escape(tmp_path: Path) -> None:
    outside = tmp_path.parent / "outside.py"
    outside.write_text("outside\n", encoding="utf-8")
    (tmp_path / "link.py").symlink_to(outside)
    with pytest.raises(ValueError, match="escapes workspace root"):
        mutation_workspace.confined_path(tmp_path, "link.py", must_exist=True)


def test_porcelain_rejects_unsafe_path_before_filesystem_probe(tmp_path: Path, monkeypatch) -> None:
    probes: list[object] = []
    monkeypatch.setattr(mutation_workspace.os.path, "lexists", lambda path: probes.append(path) or False)
    with pytest.raises(ValueError, match="unsafe workspace path"):
        mutation_workspace.parse_porcelain_z(b"?? ../outside.py\0", tmp_path)
    assert probes == []


def test_delete_plus_untracked_replacement_has_one_final_copy(tmp_path: Path) -> None:
    replacement = tmp_path / "tests/unit/replaced.py"
    replacement.parent.mkdir(parents=True)
    replacement.write_text("replacement\n", encoding="utf-8")
    payload = b"D  tests/unit/replaced.py\0?? tests/unit/replaced.py\0"
    assert mutation_workspace.parse_porcelain_z(payload, tmp_path) == (
        mutation_workspace.OverlayChange("copy", "tests/unit/replaced.py"),
    )
```

- [ ] **Step 3: Run RED**

Run: `pytest -m unit tests/unit/test_mutation_workspace.py -v`

Expected: behavioral FAIL with `NotImplementedError: not implemented`, not a collection/import failure.

- [ ] **Step 4: Implement the minimal typed data and path boundary**

```python
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Literal


@dataclass(frozen=True)
class OverlayChange:
    action: Literal["copy", "delete"]
    path: str


def confined_path(root: Path, relative: str, *, must_exist: bool) -> Path:
    candidate = Path(relative)
    if relative in {"", "."} or candidate.is_absolute() or ".." in candidate.parts or ".git" in candidate.parts:
        raise ValueError(f"unsafe workspace path: {relative}")
    resolved_root = root.resolve()
    lexical = resolved_root / candidate
    resolved_parent = lexical.parent.resolve(strict=must_exist)
    if not resolved_parent.is_relative_to(resolved_root):
        raise ValueError(f"workspace path escapes workspace root: {relative}")
    if must_exist and not os.path.lexists(lexical):
        raise ValueError(f"workspace path does not exist: {relative}")
    if lexical.is_symlink() and not lexical.resolve(strict=True).is_relative_to(resolved_root):
        raise ValueError(f"workspace path escapes workspace root: {relative}")
    return lexical
```

- [ ] **Step 5: Implement exact porcelain v1 `-z` parsing**

Parse `XY path\0`; for `R`/`C`, consume the following original path. Collect affected
paths first, including both rename paths, through `_affected_paths_from_porcelain_z`. Validate every affected path
with `confined_path(real_root, path, must_exist=False)` before the first filesystem probe. Then resolve exactly one
final action per validated lexical path: `copy` when `os.path.lexists(validated_path)` and `delete` otherwise. Thus
`D  path\0?? path\0` copies the untracked replacement that exists in the current
filesystem instead of sorting a copy/delete collision. Reject `!!`, which should never
be returned without `--ignored`, reject conflicting rename encodings, and return
`tuple(sorted(changes, key=lambda change: change.path))`. The Step 2 tests cover staged-delete-plus-untracked
replacement; add the same pre-implementation RED assertions for a deleted path (one delete) and duplicate status
records (one final action) before editing this parser.

- [ ] **Step 6: Refactor the pure boundary**

Add Google-style docstrings. Keep `_affected_paths_from_porcelain_z` independent of filesystem access; the public
`parse_porcelain_z` is deliberately filesystem-dependent because it resolves the current final state. Keep
`confined_path` independent of Git subprocess behavior.

- [ ] **Step 7: Run GREEN and mypy**

```bash
pytest -m unit tests/unit/test_mutation_workspace.py -v
mypy scripts/mutation_workspace.py
```

Expected: PASS.

- [ ] **Step 8: Commit**

```bash
git add scripts/mutation_workspace.py tests/unit/test_mutation_workspace.py
git commit -m "test(mutation): define safe worktree snapshot policy"
```

### Task 2: Create, overlay, verify and dispose the detached worktree

**Files:**
- Modify: `scripts/mutation_workspace.py`
- Modify: `tests/unit/test_mutation_workspace.py`

**Interfaces:**
- Produces: `MutationWorkspace(real_root: Path, sweep_root: Path, head: str, overlay_changes: tuple[OverlayChange, ...], baseline_manifest: tuple[OverlayChange, ...], baseline_status: bytes, baseline_digests: dict[str, str])`.
- Produces: `detached_head_workspace(real_root: Path, targets: tuple[str, ...], excluded_outputs: tuple[Path, ...]) -> ContextManager[MutationWorkspace]`.
- Produces: `MutationWorkspace.target_path(relative: str) -> Path` and `verify_baseline() -> None`.
- Invariant: `overlay_changes` records source-side copy/delete operations; `baseline_manifest` is independently captured
  from disposable-tree status after those operations and is the only authoritative status restoration baseline.

- [ ] **Step 1: Add the lifecycle scaffold and write the exact RED Git command test**

Add `MutationWorkspace` and `detached_head_workspace` signatures with bodies that raise
`NotImplementedError("not implemented")`, plus the `contextlib`, `shutil`, `subprocess`, and `tempfile` imports those
signatures will own. Add `collections.Counter`, `shutil`, `subprocess`, and `unittest.mock` imports to the test file,
then add this command-aware double and test. It deliberately maps seven
distinct subprocess calls instead of consuming an order-only response iterator:

```python
GitCall = tuple[tuple[str, ...], Path]


def _install_workspace_git_double(
    monkeypatch: pytest.MonkeyPatch,
    real: Path,
    sweep: Path,
    *,
    disposable_head: str = "a" * 40,
    real_status: bytes = b"",
    sweep_status: bytes = b"",
) -> tuple[Counter[GitCall], dict[GitCall, tuple[int, bytes, bytes]]]:
    head = "a" * 40
    sweep.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(mutation_workspace.tempfile, "mkdtemp", lambda **_kwargs: str(sweep))
    worktrees = (
        f"worktree {real}\0HEAD {head}\0\0"
        f"worktree {sweep}\0HEAD {head}\0detached\0\0"
    ).encode()
    responses: dict[GitCall, tuple[int, bytes, bytes]] = {
        (("git", "rev-parse", "--verify", "HEAD^{commit}"), real): (0, f"{head}\n".encode(), b""),
        (("git", "worktree", "add", "--detach", str(sweep), head), real): (0, b"", b""),
        (("git", "-C", str(sweep), "rev-parse", "--verify", "HEAD^{commit}"), real): (
            0,
            f"{disposable_head}\n".encode(),
            b"",
        ),
        (("git", "status", "--porcelain=v1", "-z", "--untracked-files=all"), real): (
            0,
            real_status,
            b"",
        ),
        (("git", "worktree", "list", "--porcelain", "-z"), real): (0, worktrees, b""),
        (("git", "-C", str(sweep), "status", "--porcelain=v1", "-z", "--untracked-files=all"), real): (
            0,
            sweep_status,
            b"",
        ),
        (("git", "worktree", "remove", "--force", str(sweep)), real): (0, b"", b""),
    }
    seen: Counter[GitCall] = Counter()

    def fake_run(command: list[str], *, cwd: Path, **_kwargs: object) -> subprocess.CompletedProcess[bytes]:
        key = (tuple(command), Path(cwd))
        assert key in responses, f"unexpected subprocess call: {key}"
        seen[key] += 1
        returncode, stdout, stderr = responses[key]
        if command[1:3] == ["worktree", "remove"] and returncode == 0 and sweep.exists():
            shutil.rmtree(sweep)
        return subprocess.CompletedProcess(command, returncode, stdout, stderr)

    monkeypatch.setattr(mutation_workspace.subprocess, "run", fake_run)
    return seen, responses


def test_workspace_adds_the_captured_head_detached_and_removes_exact_path(tmp_path: Path, monkeypatch) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    seen, responses = _install_workspace_git_double(monkeypatch, real, sweep)
    with mutation_workspace.detached_head_workspace(real, (), ()) as workspace:
        assert workspace.sweep_root == sweep
        assert workspace.head == "a" * 40
        assert workspace.overlay_changes == ()
        assert workspace.baseline_manifest == ()
    assert seen == Counter({key: 1 for key in responses})
```

- [ ] **Step 2: Write wrong-HEAD and malformed-porcelain RED tests before lifecycle implementation**

Use the same command-aware helper; no implementation code may precede these tests:

```python
def test_workspace_refuses_a_disposable_head_mismatch(tmp_path: Path, monkeypatch) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    _install_workspace_git_double(monkeypatch, real, sweep, disposable_head="b" * 40)
    with pytest.raises(RuntimeError, match=r"disposable HEAD.*b{40}.*captured HEAD.*a{40}"):
        with mutation_workspace.detached_head_workspace(real, (), ()):
            pytest.fail("wrong-HEAD workspace must not yield")


def test_workspace_refuses_malformed_porcelain_before_overlay(tmp_path: Path, monkeypatch) -> None:
    real = tmp_path / "real"
    real.mkdir()
    sweep = tmp_path / "sweep"
    _install_workspace_git_double(monkeypatch, real, sweep, real_status=b"R  new.py\0")
    with pytest.raises(RuntimeError, match="malformed porcelain rename"):
        with mutation_workspace.detached_head_workspace(real, (), ()):
            pytest.fail("malformed status must not yield")
```

- [ ] **Step 3: Run the lifecycle RED**

Run: `pytest -m unit tests/unit/test_mutation_workspace.py -k 'captured_head or head_mismatch or malformed_porcelain' -v`

Expected: FAIL with `NotImplementedError: not implemented`; the wrong-HEAD and malformed rows do not pass against a
partially implemented lifecycle.

- [ ] **Step 4: Implement only capture/add/verify/status/remove lifecycle**

Use exactly the seven argument-array calls keyed by Step 1, all with `cwd=real_root`. Decode the two HEAD responses as
ASCII and require equality before requesting real status. Convert `ValueError` from porcelain parsing to
`RuntimeError` with the original diagnostic. On every exit after `worktree add`, call only
`git worktree remove --force <exact sweep_root>`; if that succeeds but the validated temporary directory remains,
remove that exact directory. Do not implement overlay copying or baseline hashing in this step.

- [ ] **Step 5: Run lifecycle GREEN**

Run: `pytest -m unit tests/unit/test_mutation_workspace.py -k 'captured_head or head_mismatch or malformed_porcelain' -v`

Expected: PASS with exactly seven subprocess calls in the successful row, four calls in the wrong-HEAD row
(capture/add/verify/remove), and five in the malformed row (capture/add/verify/real-status/remove).

- [ ] **Step 6: Write RED overlay semantics, including staged-plus-worktree-reverted state**

Create a fake real/worktree pair under `tmp_path` with a modified tracked test, an
untracked test, a tracked deletion and `tests/unit/test_reverted.py` whose porcelain status is `MM` while its current
bytes equal HEAD. Name the test
`test_overlay_operations_and_baseline_manifest_diverge_for_staged_revert`. Arrange the fake disposable status to
contain only the first three changes. Then assert:

```python
assert (workspace.sweep_root / "tests/unit/test_changed.py").read_text() == "changed\n"
assert (workspace.sweep_root / "tests/unit/test_new.py").read_text() == "new\n"
assert not (workspace.sweep_root / "docs/removed.md").exists()
assert mutation_workspace.OverlayChange("copy", "tests/unit/test_reverted.py") in workspace.overlay_changes
assert mutation_workspace.OverlayChange("copy", "tests/unit/test_reverted.py") not in workspace.baseline_manifest
assert workspace.baseline_manifest == mutation_workspace.parse_porcelain_z(
    workspace.baseline_status,
    workspace.sweep_root,
)
workspace.verify_baseline()
```

This explicitly rejects `baseline_manifest == overlay_changes`: the operation list and resulting status can differ.
Also assert `overlay_changes` is sorted and includes the copied internal symlink.

- [ ] **Step 7: Write RED exclusion and registered-worktree tests**

Before overlay implementation, add `test_excluded_overlay_paths_are_rejected`,
`test_registered_nested_worktree_is_rejected`, and `test_external_symlink_overlay_is_rejected` for `.git`, each
configured output, external symlink, and a registered nested worktree. Model `git worktree list --porcelain -z` as:

```python
worktrees = f"worktree {real}\0HEAD {'a' * 40}\0\0worktree {real / 'vendor/nested'}\0HEAD {'b' * 40}\0\0".encode()
```

An overlay entry below `vendor/nested/` raises `RuntimeError("nested registered worktree")`. A normal destination below
`workspace.sweep_root` is required and must not be rejected merely because the sweep itself lives under a temporary
parent. Add `test_registered_ancestor_worktree_does_not_exclude_real_root`: include
`worktree {real.parent}` before the exact `real` row, overlay `tests/unit/test_changed.py`, and assert it is copied.
This models the authorized milestone worktree nested below the primary checkout. Assert the Git double receives the
exact `worktree list` command.

- [ ] **Step 8: Run overlay/exclusion RED tests**

Run: `pytest -m unit tests/unit/test_mutation_workspace.py -k 'overlay or reverted or excluded or registered or external_symlink' -v`

Expected: FAIL on missing overlay application and independent baseline capture.

- [ ] **Step 9: Implement validated exclusions and overlay operations minimally**

Parse registered worktree roots from the exact Git output. Exclude only a registered worktree root that is a strict
descendant of `real_root`, and only when the relative source resolves at or below that descendant. Never exclude a
source because another registered worktree is a strict ancestor of `real_root`; that is the containing primary
checkout in the authorized layout. Do not classify `sweep_root` itself as an excluded overlay destination. For every
`overlay_changes` entry, validate the source under `real_root` and destination under
`sweep_root`, reject any `.git` component and exact configured output path/subtree, copy regular files with
`shutil.copy2`, reproduce only confined symlink text, and mirror a deletion only at the validated destination. A missing
deletion is an idempotent no-op. Never derive the eventual baseline manifest from this operation list.

- [ ] **Step 10: Run overlay GREEN**

Run: `pytest -m unit tests/unit/test_mutation_workspace.py -k 'overlay or reverted or excluded or registered or external_symlink' -v`

Expected: PASS; the staged-plus-reverted copy is present only in `overlay_changes`.

- [ ] **Step 11: Write RED baseline snapshot/verification tests**

Before implementing the snapshot, add parameterized
`test_baseline_detects_content_symlink_and_deletion_drift`: modify one copied destination byte, replace an internal
symlink, and recreate a declared deletion in separate rows. Each row calls `verify_baseline()` and expects
`RuntimeError` naming the relative path and mismatch kind. Add `test_baseline_accepts_unchanged_post_overlay_state`.

- [ ] **Step 12: Run baseline-verification RED**

Run: `pytest -m unit tests/unit/test_mutation_workspace.py -k 'baseline_detects or baseline_accepts' -v`

Expected: FAIL because `verify_baseline()` does not yet compare authoritative status and content state.

- [ ] **Step 13: Capture and verify the post-overlay baseline**

Run `git status --porcelain=v1 -z --untracked-files=all` inside the disposable tree, parse it with
`real_root=workspace.sweep_root` and store its parsed semantic final state as `baseline_manifest`, without comparing it
to `overlay_changes`. Store raw `baseline_status` for diagnostics plus SHA-256 hashes of every copied regular file and
selected target, symlink text for copied links, and explicit deletion sentinels. `verify_baseline()` repeats semantic
status and hash/symlink/deletion checks against those post-overlay values. Raw `XY` bytes remain diagnostic only.

- [ ] **Step 14: Write cleanup-failure and fresh-retry RED tests**

In `test_cleanup_failure_retains_only_the_exact_orphan`, set only the exact
`git worktree remove --force <sweep>` response to `(1, b"", b"busy")`; assert the context raises
`RuntimeError(f"orphaned worktree: {sweep}: busy")`, the exact sweep directory remains, and no `git worktree prune`,
`git reset`, `git checkout`, or `git clean` command was recorded. In `test_fresh_retry_never_reuses_an_orphan`, patch `mkdtemp` with
`side_effect=[str(first_sweep), str(second_sweep)]`, invoke the context twice, and assert the second add command names
`second_sweep`, never `first_sweep`.

- [ ] **Step 15: Run cleanup/retry RED**

Run: `pytest -m unit tests/unit/test_mutation_workspace.py -k 'cleanup_failure or fresh_retry' -v`

Expected: FAIL on missing orphan diagnostic or accidental path reuse.

- [ ] **Step 16: Implement exact cleanup failure and retry behavior**

On remove failure, retain the validated directory and raise `RuntimeError` containing its absolute path and Git stderr.
Never fall back to `shutil.rmtree` after Git reports failure. A successful remove may delete only the already validated
temporary path if it remains. Call `mkdtemp(prefix="vntyper-mutation-")` once per context invocation and retain no
module-level previous path.

- [ ] **Step 17: Run full workspace GREEN**

Run: `pytest -m unit tests/unit/test_mutation_workspace.py -v`

Expected: PASS.

- [ ] **Step 18: Commit**

```bash
git add scripts/mutation_workspace.py tests/unit/test_mutation_workspace.py
git commit -m "feat(mutation): build disposable overlaid worktree"
```

### Task 3: Force pytest CWD/import isolation and prove provenance

**Files:**
- Modify: `scripts/mutation_workspace.py`
- Modify: `scripts/mutation_test.py:430-590`
- Modify: `tests/unit/test_mutation_workspace.py`
- Modify: `tests/unit/test_mutation_test.py:155-335,843-1025`

**Interfaces:**
- Changes: `run_pytest(..., *, repo_root: Path) -> TestRun` and `run_tests(..., *, repo_root: Path) -> bool`.
- Changes: `TestRun(passed: bool, output: str, returncode: int | None, timed_out: bool)`.
- Produces: `verify_import_provenance(workspace: MutationWorkspace, modules: tuple[str, ...]) -> tuple[Path, ...]`.
- Renames the live module boundary `REPO_ROOT` to `REAL_REPO_ROOT`; after this task no `REPO_ROOT` identifier remains.

- [ ] **Step 1: Write RED subprocess-root test**

Extend the existing pytest subprocess test to capture kwargs and assert:

```python
run = mutation_test.run_pytest(("tests/unit/test_scoring.py",), 7, repo_root=sweep_root)
assert seen["cwd"] == sweep_root
assert seen["env"]["PYTHONDONTWRITEBYTECODE"] == "1"
assert seen["env"]["PYTHONPATH"].split(os.pathsep)[0] == str(sweep_root)
assert str(real_root) not in seen["env"]["PYTHONPATH"].split(os.pathsep)
assert run.returncode == 0
assert run.timed_out is False
```

- [ ] **Step 2: Write RED provenance success/leak tests**

Mock a child returning two newline-separated module paths below `sweep_root`; assert the
resolved tuple is returned. Then return one path below `real_root` and assert
`RuntimeError` contains `import escaped disposable worktree` and both roots.

- [ ] **Step 3: Run RED**

Run: `pytest -m unit tests/unit/test_mutation_workspace.py tests/unit/test_mutation_test.py -k 'provenance or subprocess_disables' -v`

Expected: FAIL on missing keyword/root interfaces and missing `TestRun` fields.

- [ ] **Step 4: Implement minimal child environment and result detail**

Rename module-level `REPO_ROOT` to `REAL_REPO_ROOT` first and update every real-root guard/output reference in the same
edit. Prepend `sweep_root` to `PYTHONPATH`, remove entries equal to or resolving below `REAL_REPO_ROOT`, retain `-B`,
and pass `cwd=repo_root`. A completed process records its
integer return code; `subprocess.TimeoutExpired` returns `returncode=None` and
`timed_out=True`.

- [ ] **Step 5: Implement provenance probe**

Run `sys.executable -B -c` from `sweep_root`; convert each target path by removing `.py` and joining components with
dots (`vntyper/scripts/scoring.py` -> `vntyper.scripts.scoring`), reject non-`.py` targets, import `vntyper` and each
exact module name, print resolved `__file__` paths one per line, and write no bytecode.
Require every path to be below `workspace.sweep_root` and outside `workspace.real_root`.

- [ ] **Step 6: Thread the root through baseline/refusal calls**

Update exact signatures and calls:

```python
baseline_refusal(targets, timeout, *, repo_root)
run_pytest(test_paths, timeout, parallel=False, *, repo_root)
run_tests(test_paths, timeout, parallel=False, *, repo_root)
```

The scoped and full-tier baseline must receive the identical disposable root.

- [ ] **Step 7: Migrate every constructor, fake, and explicit-root call site**

Search for every `run_pytest`, `run_tests` and `baseline_refusal` call and make the keyword `repo_root` explicit.
Within this task migrate `_prepare_main`, every fake `run_tests`/`run_pytest`, and every two-field `TestRun(...)` at
`tests/unit/test_mutation_test.py:196-335,557-1037`; fakes accept `repo_root`, and constructors supply actual
`returncode`/`timed_out`. Run `rg -n '\bREPO_ROOT\b|TestRun\(' scripts/mutation_test.py tests/unit/test_mutation_test.py`
and require no `REPO_ROOT` match and no two-argument constructor. Remove no bytecode defense.

- [ ] **Step 8: Run GREEN and mypy**

```bash
pytest -m unit tests/unit/test_mutation_workspace.py tests/unit/test_mutation_test.py -v
mypy scripts/mutation_workspace.py scripts/mutation_test.py
```

Expected: PASS; no test or implementation uses the real root as pytest CWD.

- [ ] **Step 9: Commit**

```bash
git add scripts/mutation_workspace.py scripts/mutation_test.py tests/unit/test_mutation_workspace.py tests/unit/test_mutation_test.py
git commit -m "fix(mutation): isolate pytest imports in sweep worktree"
```

### Task 4: Thread mutation writes to the worktree and add the killed canary

**Files:**
- Modify: `scripts/mutation_test.py:127-170,602-666,1136-1244`
- Modify: `tests/unit/test_mutation_test.py`

**Interfaces:**
- Changes: `sweep_module(..., *, repo_root: Path) -> ModuleResult`.
- Changes: `generate_mutants(path: Path, *, repo_root: Path) -> list[Mutant]`.
- Produces: `CANARY_KEY = ("vntyper/scripts/scoring.py", 74, "/", "*")`.
- Produces: `canary_refusal(timeout: int, *, repo_root: Path) -> str | None`.

- [ ] **Step 1: Write RED root-write regression**

Use separate `real_root` and `sweep_root`, write identical `sample.py`, monkeypatch all
mutant tests as killed, and seed `real_root/vntyper/__pycache__/sentinel.pyc`; then assert:

```python
before = (real_root / "sample.py").read_bytes()
result = mutation_test.sweep_module("sample.py", ("tests",), 5, False, repo_root=sweep_root)
assert result.killed == result.total > 0
assert (real_root / "sample.py").read_bytes() == before
assert (real_root / "vntyper/__pycache__/sentinel.pyc").read_bytes() == b"real-cache"
assert cleared_roots and set(cleared_roots) == {sweep_root / "vntyper"}
```

- [ ] **Step 2: Write RED canary cases**

Test exact identity lookup and these results:

- return code 1, `timed_out=False` -> no refusal;
- return code 0 -> refusal containing `canary survived`;
- return code 2 -> refusal containing `canary infrastructure failure`;
- `timed_out=True` -> refusal containing `canary timed out`;
- missing exact key -> refusal containing all four key fields.

After every row, assert the canary helper restored only the canary target to its exact pre-call bytes. The Task 6
caller performs `workspace.verify_baseline()` immediately afterward; this unit does not invent an unavailable
workspace parameter. Do not reset, checkout or reapply other overlay entries.

- [ ] **Step 3: Run RED**

Run: `pytest -m unit tests/unit/test_mutation_test.py -k 'worktree or canary' -v`

Expected: FAIL because sweep/root/canary interfaces are absent.

- [ ] **Step 4: Make mutation paths explicit**

Replace every mutation/cache use of the Task 3 boundary constant `REAL_REPO_ROOT` with the keyword `repo_root`.
Resolve targets through `workspace.target_path` before passing them to mutation logic;
never resolve a target against `REAL_REPO_ROOT` after preflight. Change `generate_mutants` to compute `Mutant.path`
relative to its explicit `repo_root`; migrate its five existing unit callers at
`tests/unit/test_mutation_test.py:79-145` to `generate_mutants(path, repo_root=tmp_path)`.

- [ ] **Step 5: Implement exact canary sequence**

Generate mutants from disposable scoring source, select only `mutant.key == CANARY_KEY`,
write it in `try`, clear disposable caches, run `TARGETS[CANARY_KEY[0]]`, accept only
pytest return code 1 without harness timeout, and restore the exact pre-call target bytes in `finally`. Restoration
writes back only the selected mutation target; it does not normalize or erase unrelated overlaid changes. The caller's
immediate `workspace.verify_baseline()` checks the authoritative post-overlay status and digests.
The canary never enters reported totals.

- [ ] **Step 6: Preserve ordinary measurement semantics**

Run existing tests proving a timeout kills an ordinary mutant, a scoped survivor
escalates to the full unit tier, a red baseline writes nothing, caches clear between
mutants, and equivalent classifications remain unchanged.

- [ ] **Step 7: Refactor away residual mutation globals**

Run `rg -n 'REPO_ROOT|clear_bytecode_caches' scripts/mutation_test.py` and ensure every
mutation/test/cache call consumes its explicit disposable root; `REAL_REPO_ROOT` remains
only at the guard/output/lifecycle boundary. Migrate the five existing `generate_mutants` callers and baseline/sweep
calls introduced by this task; Task 3 already migrated `_prepare_main`, pytest fakes, and every `TestRun` constructor.

- [ ] **Step 8: Run GREEN and mypy**

```bash
pytest -m unit tests/unit/test_mutation_test.py -v
mypy scripts/mutation_test.py
```

Expected: PASS.

- [ ] **Step 9: Commit**

```bash
git add scripts/mutation_test.py tests/unit/test_mutation_test.py
git commit -m "fix(mutation): require killed mutant canary"
```

### Task 5: Install report outputs atomically at the real root

**Files:**
- Create: `scripts/mutation_output.py`
- Create: `tests/unit/test_mutation_output.py`
- Modify: `scripts/mutation_test.py:1088-1135`
- Modify: `tests/unit/test_mutation_test.py:720-840`

**Interfaces:**
- Produces: `atomic_write_text(destination: Path, content: str) -> None`.
- `write_outputs` consumes already resolved real-root paths and calls the helper.

- [ ] **Step 1: Write RED atomic success/failure tests**

```python
def test_atomic_write_replaces_complete_text(tmp_path: Path) -> None:
    destination = tmp_path / "report.md"
    destination.write_text("old\n", encoding="utf-8")
    mutation_output.atomic_write_text(destination, "new\n")
    assert destination.read_text(encoding="utf-8") == "new\n"
    assert list(tmp_path.iterdir()) == [destination]


def test_failed_replace_preserves_old_text_and_removes_temp(tmp_path: Path, monkeypatch) -> None:
    destination = tmp_path / "report.md"
    destination.write_text("old\n", encoding="utf-8")
    monkeypatch.setattr(mutation_output.os, "replace", mock.Mock(side_effect=OSError("blocked")))
    with pytest.raises(OSError, match="blocked"):
        mutation_output.atomic_write_text(destination, "new\n")
    assert destination.read_text(encoding="utf-8") == "old\n"
    assert list(tmp_path.iterdir()) == [destination]
```

Create the file with imports, scripts-path insertion and `pytestmark = pytest.mark.unit`.

- [ ] **Step 2: Run RED**

Run: `pytest -m unit tests/unit/test_mutation_output.py -v`

Before running, create an importable `mutation_output.py` scaffold whose public function raises
`NotImplementedError("not implemented")`. Expected: behavioral FAIL with that exception, not collection failure.

- [ ] **Step 3: Implement minimal atomic write**

Use `tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)`, write
UTF-8 through `os.fdopen`, `flush`, `os.fsync`, then `os.replace(temp, destination)`.
In `finally`, unlink only the exact temporary path when it remains.

- [ ] **Step 4: Add RED partial-pair and render-only tests**

Call `write_outputs` directly and patch the second atomic call to fail. Assert it raises `OSError`, the first JSON file
is complete, and the Markdown file remains old-complete; Task 6 owns conversion of that exception to `main() == 1`.
Feed the complete JSON into the existing render-only entry path and assert it rebuilds Markdown without asking for a
workspace. Run the tests and observe failure before changing `write_outputs`.

- [ ] **Step 5: Use the helper from `write_outputs`**

Render Markdown/plain text and JSON completely in memory, then call
`atomic_write_text` for JSON first and Markdown second; assert this exact destination order. If the second replacement
fails, the complete JSON is authoritative for `make mutation-render` recovery. Resolve relative CLI output paths
against the invocation real root before entering the worktree context. Keep the existing
dirty-output preflight and `render-only` behavior.

- [ ] **Step 6: Run GREEN/refactor**

```bash
pytest -m unit tests/unit/test_mutation_output.py tests/unit/test_mutation_test.py -k 'output or render' -v
mypy scripts/mutation_output.py scripts/mutation_test.py
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add scripts/mutation_output.py scripts/mutation_test.py tests/unit/test_mutation_output.py tests/unit/test_mutation_test.py
git commit -m "fix(mutation): install reports atomically"
```

### Task 6: Orchestrate overlay, guards, signals, cleanup and retries

**Files:**
- Modify: `scripts/mutation_test.py:1136-1244`
- Modify: `scripts/mutation_guard.py:1-219`
- Modify: `tests/unit/test_mutation_test.py`
- Modify: `tests/unit/test_mutation_workspace.py`

**Interfaces:**
- Consumes: `detached_head_workspace`, provenance, explicit-root baseline/canary/sweep and atomic outputs.
- Produces: `main() -> int` with exit 0 only after measurement, outputs, cleanup and real-source digest verification.
- Produces: `_real_target_digests(real_root: Path, targets: Mapping[str, object]) -> dict[str, str]` and
  `_verify_real_target_digests(real_root: Path, expected: Mapping[str, str]) -> None`.
- Produces: `_terminate(signum: int, _frame: FrameType | None) -> NoReturn` and
  `_install_signal_handlers() -> tuple[int, ...]`.

- [ ] **Step 1: Add a complete orchestration double and write the RED success-order test**

Add `contextlib`, `signal`, `Callable`, `dataclass`, and `field` imports to `tests/unit/test_mutation_test.py`, then add this helper.
It defines every event list, verification label, root and fake result used by the test; it does not call `_prepare_main`
or a real process:

```python
@dataclass
class OrchestrationHarness:
    real_root: Path
    sweep_root: Path
    real_target: Path
    output: Path
    results_json: Path
    events: list[str] = field(default_factory=list)
    sweep_roots: list[Path] = field(default_factory=list)
    commands: list[list[str]] = field(default_factory=list)


def _install_orchestration_harness(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    failure: str | None = None,
) -> OrchestrationHarness:
    real_root = tmp_path / "real"
    sweep_root = tmp_path / "sweep"
    real_target = real_root / "sample.py"
    sweep_target = sweep_root / "sample.py"
    real_target.parent.mkdir(parents=True)
    sweep_target.parent.mkdir(parents=True)
    real_target.write_text("VALUE = 1\n", encoding="utf-8")
    sweep_target.write_text("VALUE = 1\n", encoding="utf-8")
    harness = OrchestrationHarness(
        real_root=real_root,
        sweep_root=sweep_root,
        real_target=real_target,
        output=real_root / "report.md",
        results_json=real_root / "results.json",
    )
    verify_labels = ["verify-after-canary", "verify-after-sweep"]

    @dataclass
    class FakeWorkspace:
        real_root: Path
        sweep_root: Path
        head: str = "a" * 40
        overlay_changes: tuple[object, ...] = ()
        baseline_manifest: tuple[object, ...] = ()

        def target_path(self, relative: str) -> Path:
            return self.sweep_root / relative

        def verify_baseline(self) -> None:
            label = verify_labels.pop(0)
            harness.events.append(label)
            if failure == label:
                message = "canary restore" if label == "verify-after-canary" else "ordinary restore"
                raise RuntimeError(message)

    workspace = FakeWorkspace(real_root, sweep_root)

    @contextlib.contextmanager
    def fake_detached(
        requested_root: Path,
        targets: tuple[str, ...],
        excluded_outputs: tuple[Path, ...],
    ):
        assert requested_root == real_root
        assert targets == ("sample.py",)
        assert excluded_outputs == (harness.output, harness.results_json)
        harness.events.append("workspace-enter")
        try:
            yield workspace
        finally:
            harness.events.append("workspace-exit")
            if failure == "cleanup":
                raise RuntimeError(f"orphaned worktree: {sweep_root}")

    def fake_dirty_guard(targets: object, outputs: object) -> None:
        assert tuple(targets) == ("sample.py",)
        assert tuple(outputs) == (harness.output, harness.results_json)
        harness.events.append("dirty-guard")
        return None

    def fake_capture(root: Path, targets: object) -> dict[str, str]:
        assert root == real_root
        assert tuple(targets) == ("sample.py",)
        harness.events.append("real-digest-captured")
        return {"sample.py": "startup-digest"}

    def fake_verify_real(root: Path, expected: object) -> None:
        assert root == real_root
        assert expected == {"sample.py": "startup-digest"}
        assert real_target.read_text(encoding="utf-8") == "VALUE = 1\n"
        harness.events.append("real-digest-verified")

    def fake_provenance(actual_workspace: object, modules: tuple[str, ...]) -> tuple[Path, ...]:
        assert actual_workspace is workspace
        assert modules == ("sample",)
        harness.events.append("provenance")
        if failure == "provenance":
            raise RuntimeError("escaped import")
        return (sweep_root / "vntyper/__init__.py", sweep_target)

    def fake_baseline(targets: object, timeout: int, *, repo_root: Path) -> str | None:
        assert tuple(targets) == ("sample.py",)
        assert timeout == 600
        harness.sweep_roots.append(repo_root)
        harness.events.append("baseline")
        if failure == "baseline":
            return "baseline red"
        return None

    def fake_canary(timeout: int, *, repo_root: Path) -> str | None:
        assert timeout == 600
        harness.sweep_roots.append(repo_root)
        harness.events.append("canary")
        if failure == "canary":
            return "canary survived"
        return None

    def fake_sweep(
        path: str,
        tests: tuple[str, ...],
        timeout: int,
        verbose: bool,
        *,
        repo_root: Path,
    ) -> mutation_test.ModuleResult:
        assert (path, tests, timeout, verbose) == ("sample.py", ("tests/unit/x.py",), 600, False)
        harness.sweep_roots.append(repo_root)
        harness.events.append("sweep")
        if failure == "sweep":
            raise RuntimeError("sweep failed")
        return mutation_test.ModuleResult(path=path, killed=1)

    def fake_outputs(results, elapsed, output, results_json) -> None:
        assert len(results) == 1 and results[0].path == "sample.py"
        assert elapsed == 2.0
        assert output == harness.output
        assert results_json == harness.results_json
        harness.events.append("outputs")
        if failure == "outputs":
            raise OSError("json replace")
        if failure == "cleanup":
            harness.results_json.write_bytes(b"new-json")
            harness.output.write_bytes(b"new-markdown")

    monotonic = iter((10.0, 12.0))
    monkeypatch.setattr(mutation_test, "REAL_REPO_ROOT", real_root)
    monkeypatch.setattr(mutation_test, "TARGETS", {"sample.py": ("tests/unit/x.py",)})
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "mutation_test.py",
            "--module",
            "sample.py",
            "--output",
            str(harness.output),
            "--results-json",
            str(harness.results_json),
        ],
    )
    monkeypatch.setattr(mutation_test.signal, "signal", lambda *_args: None)
    monkeypatch.setattr(mutation_test, "_refuse_if_dirty", fake_dirty_guard)
    monkeypatch.setattr(mutation_test, "_real_target_digests", fake_capture)
    monkeypatch.setattr(mutation_test, "_verify_real_target_digests", fake_verify_real)
    monkeypatch.setattr(mutation_test, "detached_head_workspace", fake_detached)
    monkeypatch.setattr(mutation_test, "verify_import_provenance", fake_provenance)
    monkeypatch.setattr(mutation_test, "baseline_refusal", fake_baseline)
    monkeypatch.setattr(mutation_test, "canary_refusal", fake_canary)
    monkeypatch.setattr(mutation_test, "sweep_module", fake_sweep)
    monkeypatch.setattr(mutation_test, "write_outputs", fake_outputs)
    monkeypatch.setattr(mutation_test.time, "monotonic", lambda: next(monotonic))
    return harness


def test_main_orders_roots_measurement_cleanup_and_digest_verification(tmp_path: Path, monkeypatch) -> None:
    harness = _install_orchestration_harness(tmp_path, monkeypatch)
    assert mutation_test.main() == 0
    assert harness.events == [
        "dirty-guard",
        "real-digest-captured",
        "workspace-enter",
        "provenance",
        "baseline",
        "canary",
        "verify-after-canary",
        "sweep",
        "verify-after-sweep",
        "outputs",
        "workspace-exit",
        "real-digest-verified",
    ]
    assert harness.sweep_roots == [harness.sweep_root] * 3
```

- [ ] **Step 2: Run the success-order RED**

Run: `pytest -m unit tests/unit/test_mutation_test.py::test_main_orders_roots_measurement_cleanup_and_digest_verification -v`

Expected: FAIL on the absent orchestration interfaces or wrong event order; it must not fail on an undefined helper
name, event list, root or result shape.

- [ ] **Step 3: Implement the minimal successful orchestration path**

Resolve relative output arguments against `REAL_REPO_ROOT`, call `_refuse_if_dirty`, capture
`_real_target_digests`, and enter `detached_head_workspace(REAL_REPO_ROOT, tuple(targets),
(output_path, results_path))`. Inside the context, call `verify_import_provenance(workspace,
tuple(Path(path).with_suffix("").as_posix().replace("/", ".") for path in targets))`, then baseline, canary,
`workspace.verify_baseline()`, each `sweep_module(..., repo_root=workspace.sweep_root)` followed by
`workspace.verify_baseline()`, and `write_outputs`. After context exit call `_verify_real_target_digests`. Implement
the two digest helpers with SHA-256 over exactly the selected real targets. Do not add signal registration or failure
catching in this step.

- [ ] **Step 4: Run the successful path GREEN**

Run: `pytest -m unit tests/unit/test_mutation_test.py::test_main_orders_roots_measurement_cleanup_and_digest_verification -v`

Expected: PASS with the exact event and root lists above.

- [ ] **Step 5: Verify the overlay/guard dependency before failure orchestration**

Run:

```bash
pytest -m unit \
  tests/unit/test_mutation_workspace.py::test_overlay_operations_and_baseline_manifest_diverge_for_staged_revert \
  tests/unit/test_mutation_test.py::test_a_dirty_target_file_stops_the_sweep_before_anything_is_rewritten \
  tests/unit/test_mutation_test.py::test_main_orders_roots_measurement_cleanup_and_digest_verification -v
```

Expected: PASS. The first node proves ordinary dirty non-target state is overlaid using distinct operation/baseline
records; the second proves a dirty selected target refuses before workspace creation; the third proves baseline,
canary and sweep receive the disposable root.

- [ ] **Step 6: Write RED stage-failure rows using the concrete harness**

```python
@pytest.mark.parametrize(
    ("failure", "required_text", "forbidden_event"),
    [
        ("provenance", "escaped import", "baseline"),
        ("baseline", "baseline red", "canary"),
        ("canary", "canary survived", "sweep"),
        ("verify-after-canary", "canary restore", "sweep"),
        ("sweep", "sweep failed", "outputs"),
        ("verify-after-sweep", "ordinary restore", "outputs"),
    ],
)
def test_orchestration_unwinds_without_publishing_on_stage_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    failure: str,
    required_text: str,
    forbidden_event: str,
) -> None:
    harness = _install_orchestration_harness(tmp_path, monkeypatch, failure=failure)
    harness.output.write_bytes(b"old-markdown")
    harness.results_json.write_bytes(b"old-json")
    assert mutation_test.main() == 1
    assert required_text in capsys.readouterr().err
    assert "workspace-exit" in harness.events
    assert "real-digest-verified" in harness.events
    assert forbidden_event not in harness.events
    assert harness.output.read_bytes() == b"old-markdown"
    assert harness.results_json.read_bytes() == b"old-json"
    assert harness.real_target.read_bytes() == b"VALUE = 1\n"
    assert not any(
        command[:2] in (["git", "reset"], ["git", "checkout"]) or "clean" in command
        for command in harness.commands
    )
```

The canary and every individual module sweep must call `workspace.verify_baseline()` before another sweep or any
output. A mismatch returns 1 and publishes nothing.

- [ ] **Step 7: Write RED output, cleanup-failure and retry rows**

```python
def test_main_reports_output_failure_and_still_unwinds(tmp_path: Path, monkeypatch, capsys) -> None:
    harness = _install_orchestration_harness(tmp_path, monkeypatch, failure="outputs")
    harness.output.write_bytes(b"old-markdown")
    harness.results_json.write_bytes(b"old-json")
    assert mutation_test.main() == 1
    assert "json replace" in capsys.readouterr().err
    assert harness.events[-2:] == ["workspace-exit", "real-digest-verified"]
    assert harness.real_target.read_bytes() == b"VALUE = 1\n"


def test_cleanup_failure_is_nonzero_and_a_fresh_retry_succeeds(tmp_path: Path, monkeypatch, capsys) -> None:
    failed = _install_orchestration_harness(tmp_path / "first", monkeypatch, failure="cleanup")
    assert mutation_test.main() == 1
    assert f"orphaned worktree: {failed.sweep_root}" in capsys.readouterr().err
    assert failed.output.read_bytes() == b"new-markdown"
    assert failed.results_json.read_bytes() == b"new-json"
    assert failed.real_target.read_bytes() == b"VALUE = 1\n"
    assert "real-digest-verified" in failed.events

    retry = _install_orchestration_harness(tmp_path / "retry", monkeypatch)
    assert retry.sweep_root != failed.sweep_root
    assert mutation_test.main() == 0
    assert "outputs" in retry.events
    assert retry.real_target.read_bytes() == b"VALUE = 1\n"
```

Retain Task 5's direct JSON-first/Markdown-second partial-pair test; Task 6 tests only conversion to exit 1, context
unwinding and digest verification, so it does not duplicate the atomic writer's filesystem-state matrix.

- [ ] **Step 8: Run every failure-path RED**

Run: `pytest -m unit tests/unit/test_mutation_test.py -k 'orchestration or provenance_failure or baseline_failure or canary_failure or restoration_failure or output_failure or cleanup_failure or retry' -v`

Expected: behavioral failures name the missing exit/unwind/digest paths; no collection error or undefined helper.

- [ ] **Step 9: Implement exact failure/unwind behavior**

Keep one `status = 0` and no early return after digest capture. Catch `(RuntimeError, OSError, KeyboardInterrupt)`
outside the worktree context, print the exact exception to stderr and set `status = 1`. In an outer `finally`, call
`_verify_real_target_digests` whenever startup digests were captured; if verification itself fails, print that exact
error and force status 1 without hiding the earlier failure. Output replacement remains inside the context so cleanup
failure can coexist with complete outputs. Return `status` only after this `finally`.

- [ ] **Step 10: Run failure-path GREEN**

Run the exact Step 8 command again.

Expected: PASS; every row exits the context and verifies the real target.

- [ ] **Step 11: Write signal-handler RED before widening registration**

Add this parameterized test before changing signal setup:

```python
@pytest.mark.parametrize("signal_name", ["SIGTERM", "SIGHUP", "SIGQUIT"])
def test_handled_signal_unwinds_workspace(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    signal_name: str,
) -> None:
    signum = getattr(signal, signal_name, None)
    if signum is None:
        pytest.skip(f"{signal_name} is unavailable on this platform")
    harness = _install_orchestration_harness(tmp_path, monkeypatch)
    registered: dict[int, Callable[[int, object], None]] = {}
    monkeypatch.setattr(mutation_test.signal, "signal", lambda number, handler: registered.__setitem__(number, handler))

    def interrupting_sweep(*_args: object, **_kwargs: object) -> mutation_test.ModuleResult:
        handler = registered[signum]
        handler(signum, None)
        raise AssertionError("signal handler must raise")

    monkeypatch.setattr(mutation_test, "sweep_module", interrupting_sweep)
    assert mutation_test.main() == 1
    assert "workspace-exit" in harness.events
    assert "terminated by signal" in capsys.readouterr().err
```

Also add `test_sigint_is_not_replaced`: capture registrations, call `_install_signal_handlers()`, and assert
`signal.SIGINT not in registered`.

- [ ] **Step 12: Run signal RED**

Run: `pytest -m unit tests/unit/test_mutation_test.py -k 'handled_signal or sigint_is_not_replaced' -v`

Expected: SIGHUP/SIGQUIT fail because only SIGTERM is currently registered; the test must run before the widened
registration implementation.

- [ ] **Step 13: Install best-effort signal unwinding minimally**

Move `_terminate` to module scope, annotate it, and have `_install_signal_handlers` register exactly this available set
and return the registered integer signal numbers:

```python
def _terminate(signum: int, _frame: FrameType | None) -> NoReturn:
    raise KeyboardInterrupt(f"terminated by signal {signum}")


def _install_signal_handlers() -> tuple[int, ...]:
    registered: list[int] = []
    for name in ("SIGTERM", "SIGHUP", "SIGQUIT"):
        signum = getattr(signal, name, None)
        if signum is not None:
            signal.signal(signum, _terminate)
            registered.append(signum)
    return tuple(registered)
```

Import `FrameType` from `types` and `NoReturn` from `typing`. Call `_install_signal_handlers()` once after the dirty
guard and before workspace entry. SIGINT keeps Python's `KeyboardInterrupt` behavior. Tests invoke the captured handlers
and assert workspace `finally` runs. Do not claim this covers SIGKILL.

- [ ] **Step 14: Run signal GREEN**

Run: `pytest -m unit tests/unit/test_mutation_test.py -k 'handled_signal or sigint_is_not_replaced' -v`

Expected: PASS for every signal available on the host.

- [ ] **Step 15: Refactor guard wording with fixed replacement claims**

Replace every claim that mutation targets are rewritten in the real tree with these two facts: selected targets must
be clean so the detached baseline corresponds to committed target bytes; requested real outputs must be clean because
they are replaced after measurement. `format_indeterminate_refusal` must contain `cannot verify selected targets and
requested outputs`; `format_dirty_tree_refusal` must contain `selected targets define the committed measurement
baseline` and `requested outputs are replaced in the real checkout`. Remove the obsolete end-of-run
`format_unrestored_warning` call from `mutation_test.py`; the real digest verifier and disposable baseline verifier now
own those failures. Keep all `DESTRUCTIVE_ADVICE` absence assertions.

- [ ] **Step 16: Run focused GREEN**

Run:

```bash
pytest -m unit tests/unit/test_mutation_workspace.py tests/unit/test_mutation_output.py tests/unit/test_mutation_test.py -v
mypy scripts/mutation_workspace.py scripts/mutation_output.py scripts/mutation_test.py scripts/mutation_guard.py
```

Expected: PASS.

- [ ] **Step 17: Commit and close #208**

```bash
git add scripts/mutation_test.py scripts/mutation_guard.py tests/unit/test_mutation_test.py tests/unit/test_mutation_workspace.py
git commit -m "fix(mutation): orchestrate crash-safe worktree sweeps" -m "Closes #208"
```

### Task 7: Update generated narrative and prove a real killed canary

**Files:**
- Modify: `scripts/mutation_test.py:1-115,780-980`
- Modify: `Makefile:285-325`
- Modify: `tests/unit/test_mutation_test.py`
- Modify: `docs/development/mutation-testing.md` through the generator only
- Modify: `docs/development/mutation-results.json` only after an accepted measurement

**Interfaces:**
- Consumes: completed safe harness.
- Produces: operator documentation matching worktree/overlay/canary/cleanup behavior.

- [ ] **Step 1: Write RED generated-copy assertions**

Replace old tests that require the in-place warning with assertions that generated
Markdown includes all of these exact phrases:

```python
for phrase in (
    "disposable detached worktree",
    "current non-ignored working state",
    "known-killed canary",
    "real production source is never mutated",
    "cleanup is best effort",
):
    assert phrase in page
```

Also assert the page no longer says `Writes each mutant over the real file`.

- [ ] **Step 2: Run RED narrative tests**

Run: `pytest -m unit tests/unit/test_mutation_test.py -k 'page or markdown' -v`

Expected: FAIL on missing worktree narrative and retained in-place wording.

- [ ] **Step 3: Update generator and Makefile minimally**

Update the module docstring and `format_markdown` source text—not the generated page by
hand. Replace Makefile's in-place/build warning with: the sweep uses a detached overlaid
worktree, the real tree is never mutated, and an uncatchable crash can leave only an
orphan disposable worktree.

- [ ] **Step 4: Run unit GREEN and renderer round trip**

```bash
pytest -m unit tests/unit/test_mutation_test.py -k 'page or markdown' -v
make mutation-render
cp docs/development/mutation-testing.md /tmp/mutation-testing.once.md
make mutation-render
cmp /tmp/mutation-testing.once.md docs/development/mutation-testing.md
```

Expected: PASS and byte-identical second render.

- [ ] **Step 5: Run the objective real-worktree scoring canary**

```bash
real_digest_before=$(sha256sum vntyper/scripts/scoring.py)
worktrees_before=$(git worktree list --porcelain)
set +e
python scripts/mutation_test.py --module scoring.py --verbose \
  --output /tmp/vntyper-mutation-canary.md \
  --results-json /tmp/vntyper-mutation-canary.json
mutation_status=$?
set -e
real_digest_after=$(sha256sum vntyper/scripts/scoring.py)
worktrees_after=$(git worktree list --porcelain)
test "$mutation_status" -eq 0
test "$real_digest_before" = "$real_digest_after"
test "$worktrees_before" = "$worktrees_after"
```

Expected: exit 0; output names captured HEAD, disposable path, import paths and canary
`killed (scoped)`; both equality checks pass. The overlaid uncommitted tests are used,
but no real source write occurs.

- [ ] **Step 6: Run full verification and inspect failure-free cleanup**

```bash
make type-check
make test-unit-cov
make patch-coverage
make check-all
git diff --check
git worktree list --porcelain
```

Expected: all gates PASS, branch-inclusive total stays at least 86%, patch coverage is
at least 80%, and no orphan from the acceptance run remains.

- [ ] **Step 7: Reconcile generated artifacts**

If the accepted scoped run intentionally replaces committed raw results, run one final
`make mutation-render` and verify the page matches that JSON. Otherwise leave the
committed full-sweep JSON unchanged and stage only generator-derived prose.

- [ ] **Step 8: Commit**

```bash
git add scripts/mutation_test.py Makefile tests/unit/test_mutation_test.py docs/development/mutation-testing.md
git add docs/development/mutation-results.json  # only when produced by the accepted measurement
git commit -m "docs(mutation): document isolated sweep lifecycle"
```
