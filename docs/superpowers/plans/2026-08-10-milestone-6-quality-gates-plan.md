# Milestone 6 Quality Gates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Type-check every Python script and raise scripts-only branch-inclusive unit coverage from 69.01% to at least 88.00% before adding all of `scripts/` to the unchanged 86% repository coverage gate.

**Architecture:** Keep command/process/filesystem effects in the existing script entry points and test their decisions through pure helpers, `tmp_path`, and mocks. Add one isolated scripts-only coverage target for progress measurement, then append `scripts` to the canonical coverage source only after that target passes; the existing CI unit matrix supplies the final Python 3.10-3.13 proof.

**Tech Stack:** Python 3.10-3.13, pytest, pytest-cov/coverage.py with branch coverage, mypy, Ruff, GNU Make, GitHub Actions, pandas, stdlib `unittest.mock` and `zipfile`.

## Global Constraints

- Work from the repository root; `tests/parametrization.py` reads `tests/test_data_config.json` relative to the current directory.
- Keep `[tool.coverage.report].fail_under = 86`, `COVERAGE_TARGET ?= 86`, and `PATCH_COVERAGE_TARGET ?= 80` unchanged.
- Reach a precise scripts-only branch-inclusive result of at least 88.00%, up from the reconstructed 69.01% baseline (3,817 measurable units, approximately 1,183 uncovered, approximately 725 additional units required before new-code overhead).
- Include every Python file under `scripts/`; add no coverage omit entry and no package/file exclusion.
- `scripts/download_test_data.py` and `scripts/advntr_len_differential.py` must both move from 0% to non-zero measured coverage.
- Append `scripts` to `[tool.coverage.run].source` only after the isolated 88.00% measurement passes; its final position is after `vntyper` and `docker/app`.
- The combined branch-inclusive result must be at least 86.00% independently on Python 3.10, 3.11, 3.12, and 3.13.
- Preserve `type-check-all` as two mypy invocations: runtime sources first, then `mypy vntyper/ tests/`.
- Add no `scripts.*` mypy override, `mypy_path`, blanket ignore, network-dependent unit test, Docker dependency, or reference-data dependency.
- Every new unit test file declares `pytestmark = pytest.mark.unit`; every modified or extracted function has value and failure-path tests.
- Use Python 3.10-compatible syntax, Ruff's 120-column/double-quote style, and Google-style public docstrings.
- Do not change genotype, scoring, calibration, clinical wording, external command semantics, workflow permissions, secrets, or package/release behavior.
- `docs/superpowers/` remains contributor working material and must stay under `mkdocs.yml` `exclude_docs`.

---

## File map

- Modify `Makefile`: add `SCRIPTS_COVERAGE_TARGET`, `test-scripts-cov`, widened runtime mypy scope, accurate help and scope comments.
- Modify `scripts/download_test_data.py`: add the one required annotation and make `main(argv)` return the existing CLI status so it is unit-testable without a child process.
- Create `tests/unit/test_download_test_data.py`: pure download, archive, verification, and CLI-policy tests.
- Create `tests/unit/test_advntr_len_differential.py`: exhaustive decision/oracle tests and bounded sweep/CLI tests.
- Create `tests/unit/test_golden_cohort_launcher.py`: import-resolution, launch refusal, command-recording, and
  unhandled-error tests. This plan owns initial creation; the later exception-policy plan modifies the same file.
- Create `tests/unit/test_golden_cohort_artifacts.py`: missing/malformed artifact, HTML/table/donut, normalization, and sorting tests.
- Create `tests/unit/test_golden_cohort_gate.py`: argument splitting/parser and command exit-policy tests.
- Modify `tests/unit/test_golden_cohort_runner.py`: subprocess timeout/error, cohort, load, and parallel-order branches.
- Modify `tests/unit/test_make_cram_fixtures.py`: remaining CLI/build error and empty-discovery branches.
- Modify `tests/unit/test_coverage_gate.py`: scripts-only Make contract first; final exact coverage source/no-omit invariants last.
- Modify `tests/unit/test_makefile_recipes.py`: exact mypy two-pass and temporary coverage-recipe assertions.
- Modify `tests/unit/test_workflow_consistency.py`: prove the unit coverage command is inside the four-version matrix and patch coverage remains Python 3.12-only.
- Modify `tests/unit/test_marker_hygiene.py`: pin root `pytest.ini` as the sole live pytest authority, including strict markers and `smoke`.
- Modify `pyproject.toml`: append `scripts` to coverage source only in the final source-scope task; retain branch/floor/exclusions.
- Modify `.github/workflows/ci-tests.yml`: comments only if needed to state that `scripts/**` is now type- and coverage-gated; do not alter permissions or matrix behavior.
- Modify `AGENTS.md`: retire trap 16 and replace it with the final commands/measurements only after the gates pass.

### Task 1: Add a clean scripts-only measurement interface

**Files:**
- Modify: `Makefile:4,20-42,211-227`
- Modify: `tests/unit/test_makefile_recipes.py`

**Interfaces:**
- Consumes: existing unit selection `pytest -m unit tests/unit`, `[tool.coverage.run].branch = true`.
- Produces: `SCRIPTS_COVERAGE_TARGET ?= 88` and phony target `test-scripts-cov`; no canonical `coverage.xml` or persistent `.coverage` output.

- [ ] **Step 1: Add RED recipe contract tests**

Append these literal assertions to `tests/unit/test_makefile_recipes.py`:

```python
def test_scripts_coverage_target_is_isolated_and_fixed_at_88() -> None:
    text = MAKEFILE.read_text(encoding="utf-8")
    recipe = _recipes()["test-scripts-cov"]
    assert re.search(r"^SCRIPTS_COVERAGE_TARGET\s*\?=\s*88$", text, re.MULTILINE)
    assert "--cov=scripts" in recipe
    assert "--cov-fail-under=$(SCRIPTS_COVERAGE_TARGET)" in recipe
    assert "--cov-report=term-missing" in recipe
    assert "coverage.xml" not in recipe
    assert "mktemp -d" in recipe
    assert "COVERAGE_FILE=" in recipe
    assert "test -n \"$$tmp_dir\"" in recipe
    assert "rm -rf -- \"$$tmp_dir\"" in recipe


def test_repository_coverage_thresholds_remain_independent() -> None:
    text = MAKEFILE.read_text(encoding="utf-8")
    assert re.search(r"^COVERAGE_TARGET\s*\?=\s*86$", text, re.MULTILINE)
    assert re.search(r"^PATCH_COVERAGE_TARGET\s*\?=\s*80$", text, re.MULTILINE)
```

- [ ] **Step 2: Run the RED tests**

Run: `pytest -m unit tests/unit/test_makefile_recipes.py -q`

Expected: FAIL with `KeyError: 'test-scripts-cov'` or a missing `SCRIPTS_COVERAGE_TARGET` assertion.

- [ ] **Step 3: Add the minimal Make interface**

Add `test-scripts-cov` to `.PHONY`, add the help row, and add this recipe without changing `test-unit-cov`:

```make
SCRIPTS_COVERAGE_TARGET ?= 88

test-scripts-cov:
	@tmp_dir="$$(mktemp -d)" || exit $$?; test -n "$$tmp_dir" || exit 1; coverage_file="$$tmp_dir/.coverage"; \
	trap 'rm -rf -- "$$tmp_dir"' EXIT; \
	COVERAGE_FILE="$$coverage_file" pytest -m unit tests/unit -o log_cli=false \
		--cov=scripts --cov-config=pyproject.toml --cov-report=term-missing \
		--cov-fail-under=$(SCRIPTS_COVERAGE_TARGET)
```

- [ ] **Step 4: Verify GREEN recipe behavior and the intended baseline failure**

Run: `pytest -m unit tests/unit/test_makefile_recipes.py -q`

Expected: PASS.

Run: `make test-scripts-cov`

Expected: FAIL only because precise branch-inclusive scripts coverage is approximately 69.01%, below 88%; the terminal report must list both zero modules and must not create `coverage.xml`.

- [ ] **Step 5: Refactor comments without changing recipes**

Document that the scripts-only target is a pre-source proof and that `test-unit-cov` remains canonical. Add two
recipe-execution tests that prepend a temporary fake `pytest` to `PATH`: one exits 0 and one exits 17. In both cases
capture the `COVERAGE_FILE` parent written by the fake and assert that directory no longer exists after `make` returns;
the failure case must preserve exit 17. Add a fake `mktemp` returning failure and assert no `/.coverage` path is passed
to pytest. Run those tests and `git diff --check`.

- [ ] **Step 6: Commit**

```bash
git add Makefile tests/unit/test_makefile_recipes.py
git commit -m "test(ci): add isolated scripts coverage target"
```

### Task 2: Close #204's mypy scope gap

**Files:**
- Modify: `scripts/download_test_data.py:147`
- Modify: `Makefile:26,87-118,144-165`
- Modify: `tests/unit/test_makefile_recipes.py`

**Interfaces:**
- Consumes: current `dir_counts` keys (`str`) and values (`int`).
- Produces: runtime command `mypy vntyper/ docker/app/ scripts/`; retains `mypy vntyper/ tests/` as the second pass.

- [ ] **Step 1: Add RED exact-scope assertions**

```python
def test_mypy_runtime_and_test_scopes_are_deliberately_separate() -> None:
    recipes = _recipes()
    assert "mypy vntyper/ docker/app/ scripts/" in recipes["type-check"]
    assert "tests/" not in recipes["type-check"]
    assert "mypy vntyper/ tests/" in recipes["type-check-all"]
    assert "scripts/" not in recipes["type-check-all"]
```

- [ ] **Step 2: Run RED scope and mypy commands**

Run: `pytest -m unit tests/unit/test_makefile_recipes.py::test_mypy_runtime_and_test_scopes_are_deliberately_separate -q`

Expected: FAIL because `scripts/` is absent from `type-check`.

Run: `mypy --no-incremental vntyper/ docker/app/ scripts/`

Expected: FAIL only at `scripts/download_test_data.py:147` with `[var-annotated]` for `dir_counts`.

- [ ] **Step 3: Apply the minimal type fix and scope widening**

Change the declaration to:

```python
dir_counts: dict[str, int] = {}
```

Change the runtime recipe to exactly `mypy vntyper/ docker/app/ scripts/`, update its help/echo text, and add the scope
parity explanation in the comment block immediately above the `RUFF_PATHS := ...` assignment: Ruff's production
scope and mypy's runtime scope must be reviewed together; tests are a separate mypy pass and docs are Ruff-only.

- [ ] **Step 4: Run GREEN focused verification**

Run: `mypy --no-incremental vntyper/ docker/app/ scripts/`

Expected: `Success: no issues found` with no `scripts.*` override.

Run: `pytest -m unit tests/unit/test_makefile_recipes.py -q && make type-check-all`

Expected: PASS; output visibly contains two mypy runs.

- [ ] **Step 5: Refactor only comments/help and re-run no-cache mypy**

Run: `mypy --no-incremental vntyper/ docker/app/ scripts/ && mypy --no-incremental vntyper/ tests/`

Expected: both PASS.

- [ ] **Step 6: Commit**

```bash
git add Makefile scripts/download_test_data.py tests/unit/test_makefile_recipes.py
git commit -m "fix(ci): type-check scripts runtime scope" -m "Closes #204"
```

### Task 3: Move `download_test_data.py` above zero coverage

**Files:**
- Create: `tests/unit/test_download_test_data.py`
- Modify: `scripts/download_test_data.py:314-430`

**Interfaces:**
- Consumes: `compute_md5(Path) -> str`, `download_file_requests(str, Path, int) -> None`, `extract_archive(Path, Path) -> None`, `verify_test_data(Path, Path, bool) -> tuple[bool, list[str]]`.
- Produces: `main(argv: list[str] | None = None) -> int`; direct execution remains `sys.exit(main())` and CLI exit codes remain 0/1.

- [ ] **Step 1: Create the marked test module and RED pure-helper tests**

Create the file with imports and these literal cases:

```python
import hashlib
import json
import sys
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pytest

import scripts.download_test_data as dtd

pytestmark = pytest.mark.unit


def test_compute_md5_hashes_file_in_chunks(tmp_path: Path) -> None:
    path = tmp_path / "payload"
    path.write_bytes(b"abc" * 30000)
    assert dtd.compute_md5(path) == hashlib.md5(path.read_bytes()).hexdigest()


@pytest.mark.parametrize(
    ("members", "expected", "absent"),
    [
        ({**{f"data/f{i}.txt": str(i) for i in range(10)}, "README": "r"}, "f0.txt", "data"),
        ({"a/x": "x", "b/y": "y", "root": "r"}, "a/x", "not-stripped"),
    ],
)
def test_extract_archive_handles_dominant_and_mixed_layouts(tmp_path, members, expected, absent) -> None:
    archive = tmp_path / "data.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        for name, value in members.items():
            handle.writestr(name, value)
    output = tmp_path / "out"
    dtd.extract_archive(archive, output)
    if absent == "data":
        assert (output / expected).read_text() == "0"
        assert not (output / "data").exists()
    else:
        assert (output / expected).read_text() == "x"
        assert (output / "b/y").read_text() == "y"


@pytest.mark.parametrize("member", ["data/../../escaped", "data/link/payload"])
def test_extract_archive_rejects_members_that_escape_destination(tmp_path: Path, member: str) -> None:
    archive = tmp_path / "data.zip"
    output = tmp_path / "out"
    output.mkdir()
    if member == "data/link/payload":
        (output / "link").symlink_to(tmp_path / "outside", target_is_directory=True)
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr(member, "owned")
    with pytest.raises(ValueError, match="archive member escapes extraction directory"):
        dtd.extract_archive(archive, output)
    assert not (tmp_path / "escaped").exists()
    assert not (tmp_path / "outside/payload").exists()


def test_mixed_layout_rejects_preexisting_symlink_escape(tmp_path: Path) -> None:
    archive = tmp_path / "mixed.zip"
    output = tmp_path / "out"
    output.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (output / "link").symlink_to(outside, target_is_directory=True)
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("link/payload", "owned")
        handle.writestr("other/value", "ordinary")
    with pytest.raises(ValueError, match="archive member escapes extraction directory"):
        dtd.extract_archive(archive, output)
    assert not (outside / "payload").exists()
```

- [ ] **Step 2: Run the helper tests and verify the security RED**

Run: `pytest -m unit tests/unit/test_download_test_data.py -q`

Expected: ordinary layout tests PASS; both dominant-prefix malicious cases and the balanced mixed-layout symlink case
FAIL because neither extraction path yet enforces destination confinement.

- [ ] **Step 3: Add the minimal archive-member confinement helper**

This plan explicitly authorizes this narrow existing defect correction: the new characterization proves that
`extract_archive` can write outside its caller-provided destination through traversal or a pre-existing symlink. The
fix changes only extraction confinement and preserves archive layout, contents, permissions, and normal member names;
it does not broaden the quality-gates track into general downloader cleanup.

Add `_confined_archive_target(root: Path, relative: Path) -> Path`. Reject empty/absolute member paths, every `..`
component, and any target whose resolved parent is outside `root.resolve()` (including traversal through a pre-existing
symlink). Replace the mixed-layout `extractall` call with the same explicit per-member directory/file loop used by the
dominant path, calling the helper before directory creation or file opening in both branches. Preserve dominant-prefix
stripping, permissions already preserved by live code, and all normal member names. Re-run all three extraction tests;
expect PASS.

- [ ] **Step 4: Add verification and download negative cases**

Add tests whose exact assertions are:

```python
def test_verify_reports_missing_and_md5_mismatch(tmp_path: Path) -> None:
    data = tmp_path / "data"
    data.mkdir()
    (data / "bad.txt").write_text("bad")
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"file_resources": [
        {"local_path": "tests/data", "filename": "missing.txt", "md5sum": "0"},
        {"local_path": "tests/data", "filename": "bad.txt", "md5sum": "deadbeef"},
    ]}))
    ok, errors = dtd.verify_test_data(config, data)
    assert ok is False
    assert errors[0] == f"Missing: {data / 'missing.txt'}"
    assert errors[1].startswith(f"MD5 mismatch: {data / 'bad.txt'}")


def test_download_rejects_a_short_response(tmp_path, monkeypatch) -> None:
    response = SimpleNamespace(
        headers={"content-length": "4"},
        raise_for_status=lambda: None,
        iter_content=lambda chunk_size: iter([b"abc"]),
    )
    monkeypatch.setitem(sys.modules, "requests", SimpleNamespace(get=lambda *a, **k: response))
    with pytest.raises(RuntimeError, match="Incomplete download: 3/4 bytes"):
        dtd.download_file_requests("https://example.invalid/data.zip", tmp_path / "out")
```

- [ ] **Step 5: Add RED CLI-return tests**

Build a temporary repository layout, monkeypatch `dtd.__file__`, and assert:

```python
def test_verify_only_returns_existing_verification_status(tmp_path, monkeypatch) -> None:
    script = tmp_path / "scripts/download_test_data.py"
    script.parent.mkdir()
    script.touch()
    config = tmp_path / "tests/test_data_config.json"
    config.parent.mkdir()
    config.write_text(json.dumps({"archive_file": {"url": "https://example.invalid/x"}}))
    monkeypatch.setattr(dtd, "__file__", str(script))
    monkeypatch.setattr(dtd, "verify_test_data", lambda *a, **k: (False, ["missing"]))
    assert dtd.main(["--verify-only"]) == 1
```

Run: `pytest -m unit tests/unit/test_download_test_data.py::test_verify_only_returns_existing_verification_status -q`

Expected: FAIL with `TypeError: main() takes 0 positional arguments`.

- [ ] **Step 6: Make the minimal testability refactor**

Change `main` to accept `argv`, use `parser.parse_args(argv)`, replace internal `sys.exit(0/1)` with `return 0/1`, annotate `errors: list[str] = []`, and finish with:

```python
if __name__ == "__main__":
    sys.exit(main())
```

Do not change URL, extraction, checksum, or logging text.

- [ ] **Step 7: Cover success, skip-MD5, bad zip, missing config/archive, and cleanup**

Add exact assertions: skip-MD5 succeeds with a deliberately wrong checksum; `BadZipFile` propagates; absent config and absent `archive_file` return 1; successful mocked force-download calls download, extract, final verify, returns 0, and leaves no named temporary archive; post-extraction verification failure returns 1 and still removes the temporary archive.

- [ ] **Step 8: Verify and measure**

Run: `pytest -m unit tests/unit/test_download_test_data.py -q`

Expected: PASS with no network.

Run: `make test-scripts-cov SCRIPTS_COVERAGE_TARGET=0 | tee /tmp/scripts-cov-after-download.txt`

Expected: PASS; `scripts/download_test_data.py` is present and non-zero. Record its exact statement/branch result in the PR notes.

- [ ] **Step 9: Commit**

```bash
git add scripts/download_test_data.py tests/unit/test_download_test_data.py
git commit -m "test(scripts): cover test-data downloader policy"
```

### Task 4: Move `advntr_len_differential.py` above zero coverage

**Files:**
- Create: `tests/unit/test_advntr_len_differential.py`
- Modify only if a failing characterization exposes an import-time effect: `scripts/advntr_len_differential.py`

**Interfaces:**
- Consumes all existing public functions from `generate_states()` through `main(argv)`.
- Produces no new production interface; tests bound `sweep()` by monkeypatching `generate_states`.

- [ ] **Step 1: Add characterization decision-table tests**

Create a marked unit file and assert the exact grammar:

```python
@pytest.mark.parametrize(
    ("state", "category", "differs"),
    [
        ("D50_2", "no_len_token", False),
        ("I9_2_A_LEN3", "terminal_len_single_part", False),
        ("D50_2&I9_2_A_LEN3", "terminal_len_compound", False),
        ("I9_2_A_LEN3&D50_2", "material_after_first_len", True),
        ("I9_2_A_LENX&I50_2_A_LEN3", "stray_len_literal", True),
        ("I9_2_A_LEN0&D50_2", "material_after_first_len", False),
    ],
)
def test_classification_and_oracle(state: str, category: str, differs: bool) -> None:
    assert differential.classify(state) == category
    assert differential.oracle_predicts_difference(state) is differs
```

Add series assertions: historic length for `LEN3`, `LEN3&D`, and no-LEN is `[3, 0, 0]`; signed survival for insertion/deletion/zero nets is disjoint; absolute survival reproduces the retired mixed-state result.

- [ ] **Step 2: Run the unchanged-production characterization first**

Run: `pytest -m unit tests/unit/test_advntr_len_differential.py -q`

Expected: PASS for current pure functions. This is honest characterization, not a RED implementation failure. Any real
mismatch is resolved against the module docstring, not by weakening expected values.

- [ ] **Step 3: Add bounded sweep and production-cross-check negative tests**

Monkeypatch `generate_states` to the six-state table above and patch one production arm to return the wrong index set.
Assert `cross_check_against_production(...)` names `model-only` and `production-only`; this mutation-sensitive check is
added only after Step 2 has frozen the live values. Assert bounded `sweep(max_examples=1)` retains exhaustive counts
while each class has at most one difference record. Temporarily invert the patched wrong index set so it agrees with the
model, run this single negative node and observe it fail, then restore the fixture and confirm `git diff` contains no
temporary mutation.

- [ ] **Step 4: Add CLI exit/output characterizations**

Patch `sweep` with one clean result and one result containing each failure list. Assert clean `main(["--out", str(path)]) == 0`, JSON is sorted/parseable, and the failure case returns 1 with `REGRESSION` in captured output. Include a pathogenic residue key `1` and a pure moved state so both terminal checks execute.

- [ ] **Step 5: Minimal GREEN change only if required**

The current `main(argv: list[str] | None = None) -> int` is already testable. Do not refactor it unless a failing test proves a defect; prefer monkeypatching `sweep` and `print_summary`.

- [ ] **Step 6: Focused verification and measurement**

Run: `pytest -m unit tests/unit/test_advntr_len_differential.py -q`

Expected: PASS.

Run: `make test-scripts-cov SCRIPTS_COVERAGE_TARGET=0`

Expected: PASS; both formerly zero modules are listed with non-zero coverage.

- [ ] **Step 7: Commit**

```bash
git add tests/unit/test_advntr_len_differential.py scripts/advntr_len_differential.py
git commit -m "test(scripts): cover advntr differential oracle"
```

Omit `scripts/advntr_len_differential.py` from `git add` when it did not change.

### Task 5: Cover golden-cohort launch and artifact boundaries

**Files:**
- Create: `tests/unit/test_golden_cohort_launcher.py`
- Create: `tests/unit/test_golden_cohort_artifacts.py`
- Modify only for focused pure extraction if required: `scripts/golden_cohort/launcher.py`, `scripts/golden_cohort/artifacts.py`

**Interfaces:**
- Consumes: `launcher.resolve`, `_launch_line`, `launch`; `artifacts.read_tsv`, `read_json`, `read_delimited`, `_summary_boxes`, `read_report`, `_donut_totals`, `_sorted_records`, pipeline/cohort readers.
- Produces: no new interface unless a test requires a pure parser extraction; any extraction remains typed and in its owning module.
- Produces test imports exactly as `sys.path.insert(0, str(REPO_ROOT / "scripts")); from golden_cohort import artifacts, launcher`.
- Ownership: this task creates `tests/unit/test_golden_cohort_launcher.py`; the exception-policy plan may only modify it
  after this task lands and must preserve the single-module import identity.

- [ ] **Step 1: Write launcher behavior characterizations**

Insert the repository `scripts/` directory at the front of `sys.path` before importing `golden_cohort.launcher`, matching
the existing golden-cohort unit-test bootstrap. Use a synthetic `resolve` result and exact assertions: `_launch_line`
starts with `GATE-LAUNCH`, contains `marker_state=present`, and contains `expected_marker=absent`; `launch` returns 97
when tree/marker expectations disagree; it returns 98 when the CLI raises; it returns the CLI's 0/1 unchanged. Exercise
`_record_commands` in a child Python process, not in the pytest process: the child installs the recorder, runs a fake
`Popen`, and prints the JSONL. Assert exact `command` and `shell` fields and that an `OSError` writing the log does not
stop the fake process. This prevents its process-wide replacement of `subprocess.Popen.__init__` leaking to later tests.

Start both new files with this exact single-module import convention:

```python
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))
from golden_cohort import artifacts, launcher  # noqa: E402
```

Do not also import `scripts.golden_cohort.*`; the two names would create independently patchable module objects.

- [ ] **Step 2: Run launcher characterization tests**

Run: `pytest -m unit tests/unit/test_golden_cohort_launcher.py -q`

Expected: PASS for the current documented launch contract. If it fails on an import/global effect, treat that exact
failure as the RED extraction case and patch only the owning boundary; no real CLI or process runs.

Temporarily change the launch exception outcome from 98 to 99, run the exact exception node, and observe failure. Restore
the production line immediately and confirm `git diff` contains no temporary mutation before adding artifact tests.

- [ ] **Step 3: Write artifact parser characterizations**

Use `tmp_path` and `rules=[]`. Assert missing paths return `None`; malformed JSON returns `None` and logs a warning; short TSV/delimited rows are padded with `""`; long rows are truncated to the header; `sort_key` orders rows deterministically; report HTML returns the screening/cross-match emphasis booleans and tables; donut extraction returns numeric strings in document order; non-list records pass through `_sorted_records` unchanged.

```python
def test_delimited_rows_are_normalized_to_the_header_width(tmp_path: Path) -> None:
    path = tmp_path / "rows.tsv"
    path.write_text("a\tb\n1\n2\t3\t4\n", encoding="utf-8")
    assert artifacts.read_delimited(path, "\t", []) == {
        "columns": ["a", "b"],
        "rows": [{"a": "1", "b": ""}, {"a": "2", "b": "3"}],
    }
```

- [ ] **Step 4: Run artifact tests and make only minimal parser fixes**

Run: `pytest -m unit tests/unit/test_golden_cohort_artifacts.py -q`

Expected: PASS after any required focused fix. Do not change normalisation rules or artifact names.

- [ ] **Step 5: Measure recovered units**

Run: `make test-scripts-cov SCRIPTS_COVERAGE_TARGET=0`

Expected: PASS and a higher exact aggregate. If either target module remains below 80%, add branch cases only for displayed missing lines before committing.

- [ ] **Step 6: Commit**

```bash
git add tests/unit/test_golden_cohort_launcher.py tests/unit/test_golden_cohort_artifacts.py scripts/golden_cohort/launcher.py scripts/golden_cohort/artifacts.py
git commit -m "test(golden-cohort): cover launch and artifact failures"
```

### Task 6: Cover gate dispatch and runner partial failures

**Files:**
- Create: `tests/unit/test_golden_cohort_gate.py`
- Modify: `tests/unit/test_golden_cohort_runner.py`
- Modify only if tests expose a defect: `scripts/golden_cohort_gate.py`, `scripts/golden_cohort/runner.py`

**Interfaces:**
- Consumes: gate `_split_launch_argv`, `build_parser`, `cmd_matrix`, `cmd_probe`, `cmd_launch`, `cmd_run`, `cmd_compare`, `main`; runner `_run_one`, `run_side`, `_run_cohorts`, `load_side`.
- Produces: no new runtime interface.
- Produces gate import by `importlib.util.spec_from_file_location("golden_cohort_gate_under_test", REPO_ROOT / "scripts/golden_cohort_gate.py")` after inserting only `REPO_ROOT / "scripts"` into `sys.path`; runner stays `from golden_cohort import runner`.

- [ ] **Step 1: Add RED parser/dispatch table**

Assert `_split_launch_argv(["launch", "--", "--bam", "x"]) == (["launch"], ["--bam", "x"])`; no delimiter returns the entire input and `[]`; missing subcommand raises argparse exit 2; `cmd_matrix` returns 1 for `ValueError` and an empty matrix, writes exact JSON for a valid matrix; each main subcommand calls only its matching handler and returns its status.

```python
def test_main_dispatches_only_the_selected_handler(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(gate, "cmd_matrix", lambda _args: calls.append("matrix") or 7)
    assert gate.main(["matrix", "--data-dir", "tests/data"]) == 7
    assert calls == ["matrix"]
```

- [ ] **Step 2: Run gate tests**

Run: `pytest -m unit tests/unit/test_golden_cohort_gate.py -q`

Expected: PASS after mocks isolate matrix/probe/launch/run/compare; no worktree, process, or test data is touched.

- [ ] **Step 3: Add runner partial-failure tests**

In `test_golden_cohort_runner.py`, add exact assertions for: subprocess timeout produces an exit record with `expectation_met=False`; missing required artifact names each missing path; launch marker absence fails even when CLI exit is expected; parallel results are emitted in matrix case order; cohort execution is skipped only when `skip_cohort=True`; malformed/missing `side.json` makes `load_side` return the documented failed shape or raise the current documented error.

```python
def test_parallel_results_retain_matrix_case_order(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    cases = [_case("slow"), _case("fast")]

    def fake_run_one(*, case: dict, **_kwargs: object) -> dict:
        if case["case_id"] == "slow":
            time.sleep(0.02)
        return _ok_record(case["case_id"])

    with (
        mock.patch.object(runner, "_run_one", side_effect=fake_run_one),
        mock.patch.object(runner.admissibility, "describe_tree", return_value=CLEAN_REVISION),
    ):
        record = runner.run_side(
            matrix=_matrix(cases), tree=tmp_path, run_root=tmp_path / "run", side="after",
            marker="vntyper.scripts.cohort_rules", expect_marker=True, threads=4,
            advntr_threads=8, jobs=2, timeout=60, skip_cohort=True,
        )
    assert list(record["pipeline_results"]) == ["slow", "fast"]
```

Reuse the existing `_case`, `_ok_record`, `_matrix`, and `CLEAN_REVISION` helpers already defined in that file; add
`import time`. Do not introduce a second runner module import.

- [ ] **Step 4: Run RED then minimal GREEN**

Run: `pytest -m unit tests/unit/test_golden_cohort_gate.py tests/unit/test_golden_cohort_runner.py -q`

Expected: PASS. If a current branch has no defined error shape, preserve the observed exception and assert its exact type/message rather than inventing a recovery policy.

- [ ] **Step 5: Measure and close displayed branch misses**

Run: `make test-scripts-cov SCRIPTS_COVERAGE_TARGET=0`

Expected: PASS. Add only concrete cases for missing lines shown in `golden_cohort_gate.py` and `runner.py`, rerunning the two focused files after each 2-5 minute edit.

- [ ] **Step 6: Commit**

```bash
git add tests/unit/test_golden_cohort_gate.py tests/unit/test_golden_cohort_runner.py scripts/golden_cohort_gate.py scripts/golden_cohort/runner.py
git commit -m "test(golden-cohort): cover dispatch and runner recovery"
```

### Task 7: Close the fixed remaining scripts tranche

**Files:**
- Modify: `tests/unit/test_make_cram_fixtures.py`
- Modify: `tests/unit/test_coverage_gate.py`
- No discovery-time production or test file may be added to this task. If Tasks 3-6 plus the exact cases below do not
  reach 88.00%, stop and amend this reviewed plan with named files/cases before dispatching more implementation.

**Interfaces:**
- Consumes: per-file `term-missing` report from `make test-scripts-cov`.
- Produces: precise scripts-only branch coverage at or above 88.00%, every script listed, no omission, using only the
  two named test modules in this task after the mutation-track tests have already landed.

- [ ] **Step 1: Capture the fresh denominator and gap**

Run: `make test-scripts-cov SCRIPTS_COVERAGE_TARGET=0 | tee /tmp/scripts-coverage-before-final.txt`

Expected: PASS. Record exact total and calculate `ceil((88.00 * measurable_units / 100) - covered_units)` from statements plus branches; do not use the rounded TOTAL display.

- [ ] **Step 2: Add fixture-generator failure characterizations**

Add tests for `_run` non-zero status including stderr, empty discovery, `_select_source_bams` declared-only versus `--all`, `build_fixtures` collecting a skipped conversion without losing successful fixtures, `main` missing data/config returning 1, and manifest parent creation. Assert `LossyConversionError` remains fatal for the individual conversion and no failed fixture enters the manifest.

```python
def test_run_failure_includes_stderr(monkeypatch: pytest.MonkeyPatch) -> None:
    completed = subprocess.CompletedProcess(["samtools"], 2, stdout="", stderr="broken index")
    monkeypatch.setattr(cram_fixtures.subprocess, "run", mock.Mock(return_value=completed))
    with pytest.raises(RuntimeError, match="broken index"):
        cram_fixtures._run(["samtools"])
```

- [ ] **Step 3: Add coverage-gate failure characterizations**

Extend `test_coverage_gate.py`: corrupt/missing `.coverage` makes `read_total()` return `None`; `main()` returns 1 with `could not read`; target below floor is warning-only while total below floor is failure; `write_summary` is a no-op outside GitHub and appends exact Markdown under `GITHUB_STEP_SUMMARY` when set.

```python
def test_unreadable_coverage_fails_the_gate(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    monkeypatch.setattr(sys, "argv", ["coverage_gate.py"])
    monkeypatch.setattr(coverage_gate, "read_total", lambda: None)
    assert coverage_gate.main() == 1
    assert "could not read" in capsys.readouterr().out.lower()
```

- [ ] **Step 4: Run focused GREEN checks**

Run: `pytest -m unit tests/unit/test_make_cram_fixtures.py tests/unit/test_coverage_gate.py -q`

Expected: PASS against current documented behavior. A failure is a named conformance defect, not evidence that these
characterization assertions should be weakened.

- [ ] **Step 5: Run the fixed tranche and enforce the 88.00% gate**

Run after each small test addition: `make test-scripts-cov`

Expected: PASS at `>= 88.00%`. Confirm the report lists every `.py` discovered by
`find scripts -name '*.py'`, and specifically contains non-zero rows for `download_test_data.py` and
`advntr_len_differential.py`. A result below 88 is a plan-gate failure: record the exact missing lines and return to
plan review; do not let an implementer choose additional scope ad hoc and never alter 88, 86, 86, or 80.

Add this named, pure consistency test to `tests/unit/test_coverage_gate.py` before accepting the measurement:

```python
def test_scripts_only_coverage_scope_matches_the_repository_tree() -> None:
    scripts = sorted(
        path.relative_to(coverage_gate.REPO_ROOT).as_posix()
        for path in coverage_gate.REPO_ROOT.glob("scripts/**/*.py")
    )
    makefile = (coverage_gate.REPO_ROOT / "Makefile").read_text(encoding="utf-8")
    assert scripts
    assert "--cov=scripts" in makefile
    assert "--cov-omit" not in makefile
    assert all(path.startswith("scripts/") for path in scripts)
```

Run the named node and then compare the `term-missing` file rows to `find scripts -name '*.py'`; both checks must pass.

- [ ] **Step 6: Run marker, lint, and type checks**

Run: `pytest -m unit tests/unit/test_marker_hygiene.py -q && make format-check && make type-check-all`

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add tests/unit/test_make_cram_fixtures.py tests/unit/test_coverage_gate.py
git commit -m "test(scripts): reach 88 percent branch coverage"
```

Before committing, inspect `git diff --name-only`; remove unrelated paths from the index and confirm no coverage data or `/tmp` report is staged.

### Task 8: Append `scripts` to canonical coverage source last

**Files:**
- Modify: `pyproject.toml:350-365`
- Modify: `tests/unit/test_coverage_gate.py`
- Modify: `tests/unit/test_workflow_consistency.py`
- Modify: `tests/unit/test_marker_hygiene.py`
- Modify comments only if inaccurate: `.github/workflows/ci-tests.yml:38-52,150-222`

**Interfaces:**
- Consumes: proven scripts-only `>=88.00%` result from Task 7.
- Produces: canonical source `("vntyper", "docker/app", "scripts")`, combined matrix floor `>=86.00%`, patch measurement of scripts at 80%.

- [ ] **Step 1: Re-run the prerequisite immediately before source mutation**

Run: `make test-scripts-cov`

Expected: PASS at a precise value `>=88.00%`; stop here on any failure.

- [ ] **Step 2: Add RED source and matrix invariants**

Add to `test_coverage_gate.py`:

```python
def test_coverage_source_is_complete_and_has_no_script_omission() -> None:
    config = coverage.Coverage(config_file=str(coverage_gate.PYPROJECT)).config
    assert list(config.source) == ["vntyper", "docker/app", "scripts"]
    script_paths = [str(path) for path in sorted(coverage_gate.REPO_ROOT.glob("scripts/**/*.py"))]
    assert not config.run_omit
    assert not config.report_omit
    assert script_paths
```

Add this named source-consistency test to `test_workflow_consistency.py`:

```python
def test_unit_coverage_matrix_and_patch_coverage_version_are_fixed() -> None:
    workflow = (WORKFLOWS / "ci-tests.yml").read_text(encoding="utf-8")
    assert "python-version: ['3.10', '3.11', '3.12', '3.13']" in workflow
    assert "run: make test-unit-cov" in workflow
    assert "matrix.python-version == '3.12'" in workflow
    assert "PATCH_COVERAGE_BASE" in workflow
```

Add to `test_marker_hygiene.py`:

```python
REPO_ROOT = UNIT_DIR.parents[1]


def test_root_pytest_ini_is_the_single_live_marker_authority() -> None:
    pytest_ini = (REPO_ROOT / "pytest.ini").read_text(encoding="utf-8")
    pyproject = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert "--strict-markers" in pytest_ini
    for marker in ("unit", "integration", "docker", "smoke", "slow"):
        assert f"{marker}:" in pytest_ini
    assert "[tool.pytest.ini_options]" not in pyproject
```

- [ ] **Step 3: Run RED source test**

Run: `pytest -m unit tests/unit/test_coverage_gate.py::test_coverage_source_is_complete_and_has_no_script_omission -q`

Expected: FAIL with actual source `['vntyper', 'docker/app']`.

- [ ] **Step 4: Apply the one final source mutation**

Change only this line:

```toml
source = ["vntyper", "docker/app", "scripts"]
```

Update the adjacent comment to state that root `scripts` is a filesystem source containing the project's validation instruments. Leave `branch = true`, `fail_under = 86`, `exclude_also`, and every threshold untouched.
Remove the dead `[tool.pytest.ini_options]` block from `pyproject.toml`; `pytest.ini` remains byte-unchanged and is now
the only pytest configuration authority. Do not move or reword its live marker definitions.

- [ ] **Step 5: Run GREEN combined coverage locally**

Run: `pytest -m unit tests/unit/test_coverage_gate.py tests/unit/test_workflow_consistency.py tests/unit/test_marker_hygiene.py -q && make test-unit-cov`

Expected: PASS; precise combined branch coverage is `>=86.00%`; `coverage.xml` contains `scripts/` paths.

- [ ] **Step 6: Run patch coverage against a real merge base**

Run: `make patch-coverage PATCH_COVERAGE_BASE=origin/main`

Expected: PASS at the unchanged 80% patch threshold, or “No lines with coverage information” only if the eventual diff contains no measured production lines. A failure is recovered by adding tests, not changing the threshold.

- [ ] **Step 7: Commit**

```bash
git add pyproject.toml tests/unit/test_coverage_gate.py tests/unit/test_workflow_consistency.py tests/unit/test_marker_hygiene.py .github/workflows/ci-tests.yml
git commit -m "test(ci): measure scripts in coverage gate"
```

Omit the workflow path when no comment changed.

### Task 9: Verify the four-version matrix and recover partial failures

**Files:**
- Modify only when a version-specific test defect is proven: the owning test/production file.

**Interfaces:**
- Consumes: final source scope and unchanged thresholds.
- Produces: four independent combined-coverage successes.

- [ ] **Step 1: Run Python 3.10 in an isolated environment**

Run: `uv run --isolated --python 3.10 --with-editable '.[dev]' make test-unit-cov`

Expected: PASS at `>=86.00%`.

- [ ] **Step 2: Run Python 3.11**

Run: `uv run --isolated --python 3.11 --with-editable '.[dev]' make test-unit-cov`

Expected: PASS at `>=86.00%`.

- [ ] **Step 3: Run Python 3.12**

Run: `uv run --isolated --python 3.12 --with-editable '.[dev]' make test-unit-cov`

Expected: PASS at `>=86.00%`.

- [ ] **Step 4: Run Python 3.13**

Run: `uv run --isolated --python 3.13 --with-editable '.[dev]' make test-unit-cov`

Expected: PASS at `>=86.00%`.

- [ ] **Step 5: Recover any single-version failure**

If one command fails, rerun that interpreter with
`uv run --isolated --python 3.10 --with-editable '.[dev]' pytest -m unit tests/unit -x -vv`, substituting only the
failing version from the four literal commands above; replace version-sensitive ordering/time/mock assumptions with
deterministic assertions. Re-run all four commands after the fix. Do not average results or change source/thresholds.

- [ ] **Step 6: Commit version-specific fixes only if present**

```bash
git add <the exact owning test and production paths named by the diagnosed failure>
git commit -m "test(ci): stabilize coverage across python matrix"
```

Skip this commit when no file changed.

### Task 10: Documentation, full verification, and handoff

**Files:**
- Modify: `AGENTS.md:112,479-493`
- Modify contributor docs only where current commands or figures are stated.
- Verify: `mkdocs.yml:exclude_docs`

**Interfaces:**
- Consumes: exact final mypy file count, scripts-only percentage, and four combined percentages.
- Produces: accurate contributor instructions with no stale “scripts is not type-checked/measured” trap.

- [ ] **Step 1: Write RED documentation/source assertions if an existing test owns them**

Extend the appropriate source-consistency test to assert AGENTS no longer contains “scripts/ is linted and formatted but is not type-checked” and `mkdocs.yml` exclusion contains `superpowers/`.

- [ ] **Step 2: Update contributor text with measured facts**

Replace trap 16 with the final widened mypy scope, scripts-only exact percentage, full coverage source, and 86/86/80 semantics. Do not copy rounded totals and do not claim `scripts/` is in the tests mypy pass.

- [ ] **Step 3: Run focused documentation and quality checks**

Run: `make docs-build && make format-check && make lint && make type-check-all`

Expected: PASS; strict MkDocs output does not template or publish `docs/superpowers/`.

- [ ] **Step 4: Run the full local gate from clean coverage state**

Run: `make test-scripts-cov && make check-all && make test-unit-cov && make patch-coverage PATCH_COVERAGE_BASE=origin/main`

Expected: PASS with a freshly printed precise scripts-only result `>=88.00%`, combined result `>=86.00%`, and patch
result `>=80%`. Record the exact scripts-only percentage from this final run; do not reuse Task 7's earlier number.

- [ ] **Step 5: Run CI-equivalent installation when workflow content changed**

Run: `make ci-local`

Expected: PASS including actionlint, fresh uv venv install, format, lint, widened mypy, matrix-equivalent unit coverage, patch coverage, and strict docs.

- [ ] **Step 6: Inspect the final diff and retry all partial failures**

Run: `git diff --check && git status --short && git diff --stat`

Expected: no generated `.coverage`, `coverage.xml`, test data, or unrelated files staged. On a failed gate, fix the named test/branch and rerun that focused command followed by the full command; never weaken 88/86/86/80 or add an omission.

- [ ] **Step 7: Commit**

```bash
git add AGENTS.md docs mkdocs.yml
git commit -m "docs(ci): document scripts quality gates" -m "Closes #211"
```

Do not stage these implementation plan/specification files if they were already committed by the program-design setup.
