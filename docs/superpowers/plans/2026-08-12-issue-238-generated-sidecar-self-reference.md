# #238 Generated-Sidecar Self-Reference Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop the CRAM preflight from replacing an htslib-generated reference sidecar with a symlink to this process's own descriptor for that same file, which makes the entry self-referential and unopenable on runtimes that resolve procfd links textually (#238).

**Architecture:** One behavioural change plus one defence-in-depth layer. `ReferenceBinding` records the generated sidecar in place instead of replacing it — the branch the code already took when a proc link could not be installed. Separately, a new `consumer_reachable_identity()` helper makes both binding classes prove that an installed run-local view can be opened *through its published pathname* and reaches the bound inode, falling back to the hardlink strategy when it cannot.

**Tech Stack:** Python 3.10+ (`py310` target), pytest (`unit` marker), Ruff (format + lint, line length 120), mypy, samtools/htslib 1.20 via the `vntyper` conda env.

## Global Constraints

- Line length 120, double quotes, `target-version = "py310"`. Ruff is the only formatter — never reformat to 88.
- Run pytest from the repo root; `tests/parametrization.py` reads a relative path at collection time.
- Put the env bin first on PATH or the gate lies: `export PATH="$(conda info --base)/envs/vntyper/bin:$PATH"`.
- The `vntyper` console script is not importable in this env — use `python -m vntyper.cli` from the repo root.
- `make check-all` must pass before the PR.
- Every public function keeps a Google-style docstring with `Args:`/`Returns:`/`Raises:`.
- No new dependencies.
- Behaviour for BAM and FASTQ inputs must not change: `ReferenceBinding` is reached only under `file_format == FORMAT_CRAM` (`vntyper/scripts/alignment_preflight.py:589-591`).

## File Structure

| File | Responsibility after the change |
| --- | --- |
| `vntyper/scripts/preflight_input_io.py` | Gains `consumer_reachable_identity()` — the only place that answers "what inode does a child process reach through this pathname?". Keeps every other bounded/nonblocking read helper. |
| `vntyper/scripts/reference_binding.py` | `_InodeView` stops self-replacing generated sidecars and proves its installed proc link. `ReferenceBinding` lifetime and cleanup are unchanged. |
| `vntyper/scripts/alignment_binding.py` | `_install_proc_view` proves the installed view before recording it; the hardlink install remains the fallback. |
| `tests/unit/test_reference_binding.py` | Pins that a retained generated sidecar keeps its own name and one link. |
| `tests/unit/test_alignment_binding_lifecycle.py` | Pins that an unreachable proc view is never published. |
| `docs/plans/2026-08-12-issue-238-generated-sidecar-self-reference.md` | The spec: root cause, measured evidence, and what was ruled out. |

No file is split: `reference_binding.py` is 391 lines and `alignment_binding.py` 232, both well inside the ~650-line guideline.

---

### Task 1: Pin the defect in the reference binding

**Files:**
- Modify: `vntyper/scripts/reference_binding.py:215-221` (`_replace_generated_entry`), `:171-192` (`_install_proc_link`)
- Test: `tests/unit/test_reference_binding.py`

**Interfaces:**
- Consumes: `ReferenceBinding(input_path, output_dir, output_name, position)`, `.consumer_path`, `.bind_generated_sidecars()`, `.close()`
- Produces: `_InodeView._retain_generated_entry()` replacing `_replace_generated_entry()`; `_InodeView._install_proc_link()` loses its `replace` keyword

- [x] **Step 1: Write the failing test**

```python
def test_generated_sidecar_keeps_its_name_instead_of_linking_to_its_own_descriptor(
    tmp_path: Path,
) -> None:
    """#238: a sidecar replaced by a link to its own FD only resolves where procfs jumps dentries."""
    reference = _reference(tmp_path / "reference.fa")
    output = tmp_path / "run"
    output.mkdir()
    binding = ReferenceBinding(str(reference), str(output), "sample", 1)
    try:
        generated = Path(f"{binding.consumer_path}.fai")
        generated.write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

        binding.bind_generated_sidecars()

        assert not generated.is_symlink()
        assert generated.stat(follow_symlinks=False).st_nlink == 1
        assert generated.read_text(encoding="utf-8") == "chr1\t4\t6\t4\t5\n"
    finally:
        binding.close()
    assert not os.path.lexists(output / ".sample_reference_1")
```

- [x] **Step 2: Run it to verify it fails for the right reason**

Run: `python -m pytest tests/unit/test_reference_binding.py::test_generated_sidecar_keeps_its_name_instead_of_linking_to_its_own_descriptor -q --tb=line`
Expected: `AssertionError: assert not True` — the entry is a symlink today.

- [x] **Step 3: Make the retained entry recorded, not replaced**

```python
    def _retain_generated_entry(self) -> None:
        # The entry already *is* the retained inode, so it is only recorded, never
        # replaced. Substituting a link to this descriptor would point the name at an
        # unlinked copy of itself, leaving an entry that resolves only where procfs
        # jumps dentries rather than re-walking the link text (#238). It also protects
        # nothing extra: whoever could replace the name could replace the link.
        metadata = os.lstat(self._destination)
        if not stat.S_ISREG(metadata.st_mode) or (metadata.st_dev, metadata.st_ino) != self._identity:
            raise OSError("generated reference sidecar changed before it could be retained")
        self._record_destination("regular")
```

Update the single call site in `_InodeView.__init__` from `self._replace_generated_entry()` to `self._retain_generated_entry()`. `_install_proc_link`'s `replace=True` branch is now unreachable — delete the keyword and the branch, leaving the `replace=False` body.

- [x] **Step 4: Run the test to verify it passes**

Run: `python -m pytest tests/unit/test_reference_binding.py -q`
Expected: PASS.

- [x] **Step 5: Repair the test that pinned the defect**

`test_generated_fai_is_bound_before_the_coverage_probe_reopens_it` asserted `index.is_symlink()` was `True` — it pinned the broken state as intended behaviour. Its real intent is that the coverage probe sees the retained inode under its own name, so assert the link count instead:

```python
            observed_bound_state.append((index.stat(follow_symlinks=False).st_nlink, index.read_bytes()))
...
        assert observed_bound_state == [(1, b"chr1\t4\t6\t4\t5\n")]
```

Change its annotation to `list[tuple[int, bytes]]` — mypy rejects the `bool` one.

---

### Task 2: Prove an installed view through its own pathname

**Files:**
- Modify: `vntyper/scripts/preflight_input_io.py` (new helper), `vntyper/scripts/alignment_binding.py:97-120`, `vntyper/scripts/reference_binding.py:171-192`
- Test: `tests/unit/test_alignment_binding_lifecycle.py`

**Interfaces:**
- Consumes: `_open_nonblocking(path)` (already in `preflight_input_io.py`); `AlignmentBinding._descriptor_identity`, `_InodeView._identity`
- Produces: `consumer_reachable_identity(path: str | Path) -> tuple[tuple[int, int] | None, str | None]` — the reached `(device, inode)` and `None`, or `None` and an actionable reason

- [x] **Step 1: Write the failing test**

```python
def test_unreachable_proc_view_falls_back_to_a_hardlink(tmp_path: Path) -> None:
    """#238: a view a consumer cannot open through its own pathname is never published."""
    alignment = tmp_path / "input.cram"
    alignment.write_bytes(b"CRAM input")
    view = tmp_path / "run" / "input.cram"
    view.parent.mkdir()
    binding = AlignmentBinding(str(alignment))
    try:
        with patch(
            "vntyper.scripts.alignment_binding.consumer_reachable_identity",
            return_value=(None, "Too many levels of symbolic links (errno 40)"),
        ):
            binding.install_view(view)

        assert not view.is_symlink()
        assert view.stat(follow_symlinks=False).st_ino == alignment.stat().st_ino
        assert view.read_bytes() == b"CRAM input"
    finally:
        binding.close()
    assert not os.path.lexists(view)
```

- [x] **Step 2: Run it to verify it fails for the right reason**

Run: `python -m pytest tests/unit/test_alignment_binding_lifecycle.py::test_unreachable_proc_view_falls_back_to_a_hardlink -q --tb=line`
Expected: `AttributeError: … does not have the attribute 'consumer_reachable_identity'`.

- [x] **Step 3: Add the helper**

```python
def consumer_reachable_identity(path: str | Path) -> tuple[tuple[int, int] | None, str | None]:
    """Return the inode a child process reaches by opening an installed run-local view.

    The installer's own checks inspect the descriptor and the directory entry. Neither
    proves that the published pathname resolves, which is the only thing an external
    tool can use (#238).

    Args:
        path: Installed run-local view pathname handed to external tools.

    Returns:
        The reached ``(device, inode)`` and ``None``, or ``None`` and an actionable
        reason when the pathname cannot be opened at all.
    """
    try:
        descriptor = _open_nonblocking(path)
    except OSError as error:
        return None, f"{error.strerror} (errno {error.errno})"
    try:
        metadata = os.fstat(descriptor)
    except OSError as error:
        return None, f"{error.strerror} (errno {error.errno})"
    finally:
        os.close(descriptor)
    return (metadata.st_dev, metadata.st_ino), None
```

- [x] **Step 4: Wire it into both binding classes**

In `AlignmentBinding._install_proc_view`, after the final `os.lstat(destination)` and before `_record_view`:

```python
        reachable, reason = consumer_reachable_identity(destination)
        if reachable != self._descriptor_identity:
            # An installed view an external tool cannot open through its own pathname is
            # not a view. Leave the entry for the hardlink install's atomic replace (#238).
            logger.warning(
                f"Run-local alignment view {destination} does not reach the bound alignment "
                f"through its own pathname ({reason or 'identity mismatch'}); using a hardlink view instead."
            )
            return False
```

In `_InodeView._install_proc_link`, after `os.symlink(...)` succeeds. This one *must* remove its own entry, because `_install_source_view` re-checks `os.path.lexists(self._destination)` before the hardlink path and would otherwise reject its own leftover:

```python
    def _remove_own_symlink(self) -> None:
        with suppress(OSError):
            if os.readlink(self._destination) == self._proc_target:
                os.unlink(self._destination)
```

```python
        reachable, reason = consumer_reachable_identity(self._destination)
        if reachable != self._identity:
            logger.warning(
                f"Run-local reference view {self._destination} does not reach the bound reference "
                f"through its own pathname ({reason or 'identity mismatch'}); using a hardlink view instead."
            )
            self._remove_own_symlink()
            return False
```

Import the helper in both modules. `preflight_input_io` imports only the standard library, so neither import creates a cycle.

- [x] **Step 5: Run the binding suites**

Run: `python -m pytest tests/unit/test_reference_binding.py tests/unit/test_alignment_binding_lifecycle.py -q`
Expected: 40 passed.

---

### Task 3: Gate, evidence, and commit

**Files:**
- Create: `docs/plans/2026-08-12-issue-238-generated-sidecar-self-reference.md`

- [x] **Step 1: Format and run the full gate**

```bash
export PATH="$(conda info --base)/envs/vntyper/bin:$PATH"
make format && make check-all
```

Expected: `✓ All checks passed (full suite)`, 5923 unit tests.

- [x] **Step 2: Prove it on real CRAM data, not only in unit tests**

The generated-sidecar path only runs when the reference has **no** `.fai` beside it. Build an input directory with `tests/data/cram/example_40cf_hg38_subset.cram` and `reference/alignment/chr1.hg38.fa` (copy the FASTA, *not* its `.fai`), then run the reporter's command shape against `vntyper:latest` with a watcher inside the container printing the sidecar entry, its target, and that descriptor's own path.

Expected on stock 2.0.12: `reference.fa.fai -> /proc/22/fd/7` whose own path is `.../reference.fa.fai (deleted)`.
Expected with the three changed files bind-mounted over `site-packages`: `TYPE: regular file, links=1` and `Pipeline finished successfully`.

- [x] **Step 3: Commit**

```bash
git checkout -b fix/issue-238-generated-sidecar-self-reference
git add vntyper/scripts/alignment_binding.py vntyper/scripts/preflight_input_io.py \
        vntyper/scripts/reference_binding.py tests/unit/test_alignment_binding_lifecycle.py \
        tests/unit/test_reference_binding.py docs/plans/2026-08-12-issue-238-generated-sidecar-self-reference.md
git commit   # conventional commit naming the mechanism and closing #238
```

---

### Task 4: Clear the second assertion that pinned the defect

The unit suite is not enough: `make check-all` gates on the unit tier only, and a second
assertion in the integration tier also required the sidecar to be a symlink. It is only
visible by running the CRAM integration tier and comparing the failure set against `main`.

**Files:**
- Modify: `tests/integration/test_read_only_alignment_preflight.py:326`

- [x] **Step 1: Run the CRAM integration tier on the branch and on `main`, and diff the failures**

```bash
export PATH="$(conda info --base)/envs/vntyper/bin:$PATH"
make cram-fixtures
T="tests/integration/test_cram_pipeline.py tests/integration/test_cram_reference_pipeline_contract.py tests/integration/test_read_only_alignment_preflight.py"
python -m pytest $T -m integration -q -p no:warnings | grep '^FAILED tests/' | sort -u
```

Compare against the same command on `main`. Only a difference is a regression: this
machine has 11 pre-existing failures because `reference/All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`
is absent after the milestone-5 reference externalisation, so Kestrel writes no VCF.

Expected before the fix: one extra failure,
`test_unindexed_read_only_reference_uses_a_run_local_index_for_every_consumer`, at
`assert Path(f"{plan.reference_path}.fai").is_symlink()`.

- [x] **Step 2: Assert the test's real subject instead**

The test's subject is that an unindexed read-only reference gets a run-local index every
consumer can use while the operator tree stays untouched — not that the entry is a link.
Keep line 325 (`Path(plan.reference_path).is_symlink()`): the reference FASTA view has a
different source and destination, so it stays a proc link.

```python
        # The generated index is retained under its own name, not replaced by a link to
        # this process's descriptor for it, which would point the entry at itself (#238).
        generated_index = Path(f"{plan.reference_path}.fai")
        assert not generated_index.is_symlink()
        assert generated_index.stat(follow_symlinks=False).st_nlink == 1
        assert generated_index.read_text(encoding="utf-8").startswith("chr1\t")
        assert not Path(f"{reference}.fai").exists()
```

- [x] **Step 3: Re-run both tiers and confirm the failure sets are identical**

Expected: 11 failures on each side, `comm` reporting no regressions, and
`make check-all` still `✓ All checks passed (full suite)`.

---

### Task 5: Close the patient-input-tree guard bypass

The same command exposes a second defect: `validate_alignment_output_root()` compares pathnames only, so one host directory mounted at two container paths defeats it. The identical layout with a **single** mount was already rejected, so this is a bypass, not a policy change.

**Files:**
- Modify: `vntyper/scripts/alignment_target_io.py` (new `_aliased_input_ancestor`, wired into `validate_alignment_output_root`)
- Test: `tests/unit/test_alignment_pipeline_input_safety.py`
- Docs: `README.md`, `docs/user-guide/docker.md`

- [x] **Step 1: Write the failing test**

A symlinked alias cannot express this — `resolve()` already unifies symlinks, so such a test passes on the old code. The distinguishing property is two real pathnames with one inode, which only a mount provides, so the unit test injects exactly that at the layer the code consults and Step 4 proves the real mount.

```python
    real_samefile = os.path.samefile
    aliased = {str(patient.resolve()), str(output_mount.resolve())}

    def one_host_directory_two_mounts(left: str | Path, right: str | Path) -> bool:
        pair = {str(Path(left).resolve()), str(Path(right).resolve())}
        if pair == aliased:
            return True
        return real_samefile(left, right)

    with mock.patch(
        "vntyper.scripts.alignment_target_io.os.path.samefile",
        side_effect=one_host_directory_two_mounts,
    ):
        with pytest.raises(ValueError, match="same directory as the patient input tree"):
            validate_alignment_output_root(output_root, alignment, "cram")
```

Add a companion asserting the documented separate-directory layout is still accepted, so the inode check cannot regress into rejecting everything.

- [x] **Step 2: Run it to verify it fails**

Expected: `Failed: DID NOT RAISE ValueError`.

- [x] **Step 3: Compare by inode, not by name**

```python
def _aliased_input_ancestor(root_absolute: Path, input_trees: tuple[Path, ...]) -> tuple[Path, Path] | None:
    for ancestor in (root_absolute, *root_absolute.parents):
        for input_tree in input_trees:
            if ancestor != input_tree and _same_file(ancestor, input_tree):
                return ancestor, input_tree
    return None
```

Wire it in after the existing lexical loop and reject with a message naming **both** pathnames — the whole difficulty is that they look unrelated.

- [x] **Step 4: Prove it against real bind mounts, not simulated ones**

```bash
docker run --rm -v .:/opt/vntyper/input -v .:/opt/vntyper/output …        # expect the rejection
docker run --rm -v .:/opt/vntyper/input -v ./results:/opt/vntyper/output … # expect success
```

- [x] **Step 5: Update the usage instructions**

`README.md` and the `docs/user-guide/docker.md` volume-mount section must state the requirement, show the rejected form and the accepted one. Note the contract test `test_current_container_commands_use_only_the_rolling_ghcr_main_image`: every ```bash block on those surfaces containing `docker run` must name `ghcr.io/hassansaei/vntyper:main`, so write complete commands, not elided ones.

---

## Self-review

**Spec coverage.** Spec §3.1 (stop self-replacing) is Task 1; §3.2 (prove every installed view) is Task 2; §2 (measured evidence) is Task 3 Step 2; §4 (the patient-input-tree guard bypass) is Task 5.

**Placeholders.** None — every step carries the real code and the real command.

**Type consistency.** `consumer_reachable_identity` returns `tuple[tuple[int, int] | None, str | None]` in Task 2 Step 3 and is compared against `self._descriptor_identity` / `self._identity`, both `tuple[int, int] | None`, in Step 4. `_retain_generated_entry` is the name used in both its definition and its call site. `_aliased_input_ancestor` returns `tuple[Path, Path] | None` and its caller unpacks exactly that pair.
