# Milestone 6 Exception Policy Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace BLE001's stale prose-only exemption with an exact dual-count, Ruff-version-attributed, stable-symbol policy inventory and behavior-tested audit of every fail-open handler.

**Architecture:** A typed read-only adapter in `scripts/ble001_policy.py` invokes Ruff in normal and `--ignore-noqa` modes, normalizes diagnostics, and maps each row to its enclosing AST symbol. A reviewed JSON artifact classifies every diagnostic and records every category-C contract; unit tests compare live discovery to that artifact and characterize the affected production behavior without globally enabling or mechanically narrowing BLE001.

**Tech Stack:** Python 3.10-3.13, Ruff 0.14.3 reviewed baseline plus actual `ruff --version`, stdlib `ast`/`argparse`/`dataclasses`/`json`/`re`/`subprocess`, pytest, mypy, GNU Make.

## Global Constraints

- This plan starts only after the #204/#211 quality-gates plan has landed: `scripts/` is already in mypy and `[tool.coverage.run].source`, scripts-only coverage is at least 88.00%, and combined coverage is at least 86.00% on Python 3.10-3.13.
- Baseline measurements at `ebb15b26631242a3295607e4eda4c68f688cd9a2` are exactly 103 normal BLE001 diagnostics and 108 with `--ignore-noqa`; five existing explicit suppressions explain the delta.
- Issue #219's acceptance prose says “update the comment to 67”, but live no-cache Ruff 0.14.3 measures 103/108. The
  executable inventory uses the current evidence; `67` is retained only as stale issue history and is not written as a
  current count.
- Record reviewed Ruff version `ruff 0.14.3` and actual `ruff --version` for every measurement; a tool-version change alone passes only when normalized diagnostic identities are unchanged.
- Remeasure and pin the exact final normal/all pair after all reviewed edits; every delta from 103/108 must map to named inventory changes.
- Keep `BLE001` absent from `[tool.ruff.lint].select`; add no blanket per-file ignore, mass `noqa`, or generated acceptance baseline.
- Classify every all-mode diagnostic exactly once as A controlled boundary, B enumerable candidate, or C fail-open/degraded continuation.
- Identify handlers by relative path and qualified AST symbol, never by a persistent line number; module-level handlers use `<module>`.
- Every category-C symbol has one disposition, exact outcome, domain rationale, behavior-test node ID, and observability channel.
- Freeze every category-C `fail_open` record before live validation. The live validator exits 1 and lists unresolved
  nodes until later behavior tasks create or strengthen every linked node; no intermediate policy may claim complete
  success or contain a C handler without its complete record.
- A behavior node first characterizes unchanged production, then proves mutation sensitivity, and only then satisfies
  the full policy gate. A characterization expected to pass is never described as RED.
- Preserve intentional fail-open behavior only with rationale and a behavior test. Change behavior only when existing documentation, configuration, tests, or caller contract supplies the replacement outcome.
- When no alternate outcome exists, use final disposition `preserved-no-authorized-alternative`; do not invent an exit code, exception, fallback, message, or clinical policy.
- Do not mechanically narrow category-B handlers, rework `cross_match.match_variants`, change G004/logging style, introduce custom exceptions, or change completed-outcome CLI codes 0/1 and argparse usage code 2.
- Policy discovery is source-read-only: it does not import, `eval`, or execute scanned repository code and requires no network, Docker, reference data, secrets, or expanded workflow permissions.
- All new Python uses 3.10-compatible syntax, complete annotations, Google-style public docstrings, 120-column/double-quote Ruff style, and unit tests marked `pytest.mark.unit`.
- Because #211 already measures `scripts/`, this new policy script must itself satisfy the 86 repository floor and 80 patch-coverage floor.

---

## Frozen artifacts and interfaces

Create exactly these policy artifacts:

- `scripts/ble001_policy.py`: typed measurement, AST mapping, JSON loading, validation, diagnostics, and CLI.
- `scripts/ble001_policy_validation.py`: pure Ruff-scope containment, schema-field, and behavior-node resolution.
- `scripts/ble001_policy.json`: hand-reviewed complete inventory and fail-open manifest.
- `tests/unit/test_ble001_policy.py`: adapter unit tests, live exact guards, inventory completeness, negative mutations, Ruff configuration checks, and CLI tests.

The Python interface is fixed as follows:

```python
@dataclass(frozen=True, order=True)
class Diagnostic:
    path: str
    row: int
    column: int
    code: str
    message: str


@dataclass(frozen=True)
class Measurement:
    ruff_version: str
    diagnostics: tuple[Diagnostic, ...]


@dataclass(frozen=True)
class HandlerPolicy:
    path: str
    symbol: str
    normal_count: int
    all_count: int
    category: str
    rationale: str


@dataclass(frozen=True)
class FailOpenPolicy:
    path: str
    symbol: str
    disposition: str
    outcome: str
    rationale: str
    behavior_test: str
    observable_via: str


@dataclass(frozen=True)
class Policy:
    reviewed_ruff_version: str
    expected_normal: int
    expected_all: int
    expected_identities: int
    expected_categories: tuple[int, int, int]
    handlers: tuple[HandlerPolicy, ...]
    fail_open: tuple[FailOpenPolicy, ...]
```

The exact callable signatures are:

- `read_ruff_paths(makefile: Path) -> tuple[str, ...]`
- `normalize_diagnostics(payload: str, repo_root: Path) -> tuple[Diagnostic, ...]`
- `measure_ble001(repo_root: Path, scope: tuple[str, ...], *, ignore_noqa: bool, ruff_executable: str = "ruff") -> Measurement`
- `enclosing_symbol(source: str, row: int) -> str`
- `load_policy(path: Path) -> Policy`
- `validate_policy(repo_root: Path, normal: Measurement, all_handlers: Measurement, policy: Policy) -> list[str]`
- `main(argv: list[str] | None = None) -> int`

The JSON shape is fixed:

```json
{
  "reviewed_ruff_version": "ruff 0.14.3",
  "expected_counts": {
    "normal": 103,
    "ignore_noqa": 108,
    "identities": 79,
    "categories": {"A": 30, "B": 16, "C": 33}
  },
  "handlers": [
    {
      "path": "relative/path.py",
      "symbol": "ClassName.method_or_function",
      "normal_count": 1,
      "all_count": 1,
      "category": "A",
      "rationale": "One sentence stating the controlled outcome."
    }
  ],
  "fail_open": [
    {
      "path": "relative/path.py",
      "symbol": "function_name",
      "disposition": "preserved-contract",
      "outcome": "Exact structured fallback or terminal outcome.",
      "rationale": "Why continuing with this result is intentional or why no alternate outcome is authorized.",
      "behavior_test": "tests/unit/test_module.py::test_exact_behavior",
      "observable_via": "warning log, error log, result field, exception, or exit status"
    }
  ]
}
```

The displayed one-element arrays illustrate the schema only. The committed artifact contains the complete live inventory; validation rejects illustrative, incomplete, duplicate, or empty records.

### Task 0: Required dependency preflight (read-only)

**Files:**
- Verify only: `pyproject.toml`, `Makefile`, `.github/workflows/ci-tests.yml`

**Interfaces:**
- Consumes: final quality-gates source scope and current Ruff executable.
- Produces: recorded command evidence `ruff 0.14.3`, 103 normal, 108 all-mode; no file change.

- [ ] **Step 1: Verify the quality-gate prerequisite**

Run: `make test-scripts-cov && make test-unit-cov && make type-check-all`

Expected: scripts-only `>=88.00%`, combined `>=86.00%`, and mypy success including `scripts/`. If any fails, stop and finish #211; do not begin #219 against an unmeasured policy script.

- [ ] **Step 2: Verify final source and thresholds literally**

Run: `python -c 'import coverage; c=coverage.Coverage(config_file="pyproject.toml").config; print(c.source, c.branch, c.fail_under)'`

Expected: source is `['vntyper', 'docker/app', 'scripts']`, branch is `True`, floor is `86.0`; Makefile retains 86 advisory and 80 patch thresholds.

- [ ] **Step 3: Capture reviewed and actual Ruff attribution**

Run: `ruff --version`

Expected at reconstruction: `ruff 0.14.3`. If actual version differs, retain `ruff 0.14.3` as reviewed metadata and record the actual version beside both measurements.

- [ ] **Step 4: Reproduce both baseline commands**

Run: `ruff check --no-cache --select BLE001 --statistics vntyper/ docker/app/ tests/ scripts/ docs/`

Expected: exactly `103 BLE001 blind-except` and Ruff status 1 because diagnostics exist.

Run: `ruff check --no-cache --select BLE001 --ignore-noqa --statistics vntyper/ docker/app/ tests/ scripts/ docs/`

Expected: exactly `108 BLE001 blind-except`; delta exactly 5.

### Task 1: Build the typed read-only Ruff/AST adapter

**Files:**
- Create: `scripts/ble001_policy.py`
- Create: `tests/unit/test_ble001_policy.py`

**Interfaces:**
- Consumes: the frozen signatures above and Ruff JSON fields `filename`, `code`, `message`, `location.row`, `location.column` (accept `cell.start` as a compatibility spelling if emitted by the installed Ruff).
- Produces: deterministic relative diagnostics, actual Ruff version, stable enclosing symbol, actionable validation errors, CLI status 0/1.

- [ ] **Step 1: Create the marked test module and add RED `RUFF_PATHS` parser tests**

First create an importable `scripts/ble001_policy.py` scaffold with the frozen dataclasses/signatures and public
functions raising `NotImplementedError("not implemented")`; it supplies no behavior. Then start the marked test file:

Start with:

```python
from pathlib import Path

import pytest

from scripts.ble001_policy import Diagnostic, enclosing_symbol, normalize_diagnostics, read_ruff_paths

pytestmark = pytest.mark.unit
REPO_ROOT = Path(__file__).resolve().parents[2]
```

```python
def test_read_ruff_paths_reads_the_single_make_authority(tmp_path: Path) -> None:
    makefile = tmp_path / "Makefile"
    makefile.write_text("RUFF_PATHS := vntyper/ docker/app/ tests/ scripts/ docs/\n")
    assert read_ruff_paths(makefile) == ("vntyper/", "docker/app/", "tests/", "scripts/", "docs/")


def test_read_ruff_paths_rejects_missing_or_duplicate_assignments(tmp_path: Path) -> None:
    makefile = tmp_path / "Makefile"
    makefile.write_text("all:\n\ttrue\n")
    with pytest.raises(ValueError, match="exactly one RUFF_PATHS"):
        read_ruff_paths(makefile)
    makefile.write_text("RUFF_PATHS := a/\nRUFF_PATHS := b/\n")
    with pytest.raises(ValueError, match="exactly one RUFF_PATHS"):
        read_ruff_paths(makefile)
```

- [ ] **Step 2: Add RED JSON normalization tests**

Use one absolute in-root filename and assert exact relative output; use `../outside.py` and malformed JSON and assert `ValueError`. Assert sorting is path/row/column/code/message independent of Ruff output order.

- [ ] **Step 3: Add RED AST symbol tests**

```python
SOURCE = """\
try:\n    pass\nexcept Exception:\n    pass\n\ndef outer():\n    def inner():\n        try:\n            pass\n        except Exception:\n            return None\n    return inner()\n\nclass Worker:\n    def run(self):\n        try:\n            pass\n        except Exception:\n            return 1\n"""

def test_enclosing_symbol_uses_qualified_names_not_lines() -> None:
    assert enclosing_symbol(SOURCE, 3) == "<module>"
    assert enclosing_symbol(SOURCE, 10) == "outer.inner"
    assert enclosing_symbol(SOURCE, 18) == "Worker.run"
    assert enclosing_symbol("\n\n" + SOURCE, 12) == "outer.inner"
```

The blank-line assertion deliberately moves the handler while preserving identity.

- [ ] **Step 4: Run RED adapter tests**

Run: `pytest -m unit tests/unit/test_ble001_policy.py -q`

Expected: behavioral FAIL with `NotImplementedError: not implemented`, not collection/import failure.

- [ ] **Step 5: Implement dataclasses, parsing, and AST mapping minimally**

Implement the frozen signatures. AST mapping walks nested `ClassDef`, `FunctionDef`, and `AsyncFunctionDef` nodes whose `lineno <= row <= end_lineno`, then returns the deepest qualified path; syntax errors raise `ValueError` naming the relative source path from the caller.

- [ ] **Step 6: Add RED measurement subprocess tests**

Mock `subprocess.run` for four calls: each independent `measure_ble001` call executes its own `ruff --version` and
one JSON check, so the order is version, normal, version, all-mode. Assert the check argv begins exactly:

```python
["ruff", "check", "--no-cache", "--select", "BLE001", "--output-format", "json"]
```

and appends `--ignore-noqa` only for all-mode, then `--`, followed by the five scope paths. Reject option-shaped,
non-path, absolute, missing, and repository-escaping `RUFF_PATHS` tokens before invoking Ruff. Ruff result codes 0 and
1 are accepted; result code 2 raises `RuntimeError` containing stderr and actual version.

- [ ] **Step 7: Implement measurement and diagnostics**

Use `subprocess.run(..., cwd=repo_root, capture_output=True, text=True, check=False)`. Never use `shell=True`. `Measurement.ruff_version` is stripped `ruff --version`; a malformed/empty version raises `RuntimeError`.

- [ ] **Step 8: Add RED CLI tests and minimal CLI**

Assert `main(["--repo-root", str(tmp_path), "--policy", str(policy)])` returns 1 and prints each validation error; a mocked clean validator returns 0 and prints reviewed/actual version plus normal/all/category counts. Implement only `--repo-root`, `--policy`, and `--ruff` options; no write/snapshot option may mutate the reviewed JSON.

- [ ] **Step 9: Run focused GREEN, type, and patch checks**

Run: `pytest -m unit tests/unit/test_ble001_policy.py -q && mypy scripts/ble001_policy.py tests/unit/test_ble001_policy.py && ruff check scripts/ble001_policy.py tests/unit/test_ble001_policy.py`

Expected: PASS.

- [ ] **Step 10: Commit**

```bash
git add scripts/ble001_policy.py tests/unit/test_ble001_policy.py
git commit -m "test(ci): add typed BLE001 measurement adapter"
```

### Task 2: Load and validate the frozen policy schema

**Files:**
- Modify: `scripts/ble001_policy.py`
- Modify: `tests/unit/test_ble001_policy.py`

**Interfaces:**
- Consumes: frozen JSON schema, `Policy`, both `Measurement` objects.
- Produces: `load_policy(Path) -> Policy`, `validate_policy(...) -> list[str]`; a partial illustrative policy always fails.

- [ ] **Step 1: Add RED schema tests with a complete one-handler fixture**

Write a temporary JSON fixture with counts 1/1, one A handler, no C entries. Assert all dataclass fields exactly. Parametrize invalid fixtures for: missing reviewed version, negative count, normal greater than all, category outside A/B/C, empty rationale, duplicate `(path, symbol)`, absolute path, `..` path, C without fail-open row, fail-open row for non-C, invalid disposition, empty test node ID, and category totals not matching expected all.

```python
def test_policy_rejects_a_category_c_handler_without_behavior_evidence(tmp_path: Path) -> None:
    payload = _policy_payload(
        handlers=[{"path": "module.py", "symbol": "run", "normal_count": 1, "all_count": 1,
                   "category": "C", "rationale": "Returns a fallback."}],
        fail_open=[],
    )
    path = tmp_path / "policy.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="category C.*fail_open"):
        load_policy(path)
```

- [ ] **Step 2: Run RED schema tests**

Run: `pytest -m unit tests/unit/test_ble001_policy.py -k 'policy or schema' -q`

Expected: FAIL because `load_policy`/validation are absent.

- [ ] **Step 3: Implement strict JSON loading**

Reject unknown top-level and record keys so misspellings cannot silently drop policy. Accept only dispositions
`preserved-contract`, `conformed-to-existing-contract`, and `preserved-no-authorized-alternative`. Validate behavior
node IDs with this Python 3.10-compatible expression, then separately reject `.` or `..` path components:

```python
TEST_NODE_ID = re.compile(
    r"^tests/unit/(?:[A-Za-z0-9_.-]+/)*test_[A-Za-z0-9_]+\.py"
    r"(?:::[A-Za-z_][A-Za-z0-9_]*)*::test_[A-Za-z_][A-Za-z0-9_]*"
    r"(?:\[[^\]\r\n]+\])?$"
)
```

This accepts top-level nodes, nested unit paths such as
`tests/unit/web/test_tasks.py::test_input_cleanup_logs_one_error_and_attempts_every_owned_path`, class-qualified nodes,
and parametrized suffixes. Add positive tests for all four shapes and negative tests for integration paths, missing
`test_` function names, `..` components, newlines, and a path without `::`.

- [ ] **Step 4: Add RED live-identity validation tests**

Construct synthetic source and measurements. Assert errors for unclassified diagnostics, multiply counted symbols, changed normal/all counts, stale policy symbols, and version-attributed mismatch text containing both `reviewed ruff 0.14.3` and `actual ruff 0.16.0`. Assert version-only difference returns `[]` when identities and counts match.

- [ ] **Step 5: Implement validation minimally**

For every all-mode diagnostic, read source as UTF-8, map to `(path, symbol)`, and aggregate actual normal/all counts.
Compare exact maps to policy. Validate sums, normal subset, `normal <= all`, one-to-one C manifest, and
`normal.ruff_version == all_handlers.ruff_version`; a mismatch fails before comparing suppression deltas. Error output
lists added/removed stable identities and rows only as observations.

- [ ] **Step 6: Verify unit GREEN for synthetic fixtures**

Run: `pytest -m unit tests/unit/test_ble001_policy.py -k 'not live' -q`

Expected: PASS.

- [ ] **Step 7: Commit strict schema validation without an incomplete live artifact**

```bash
git add scripts/ble001_policy.py tests/unit/test_ble001_policy.py
git commit -m "test(ci): validate BLE001 policy schema"
```

### Task 3: Verify and freeze the complete 108-diagnostic inventory

**Files:**
- Create: `scripts/ble001_policy.json`
- Modify: `tests/unit/test_ble001_policy.py`
- Modify: `pyproject.toml:249-254`

**Interfaces:**
- Consumes: live normalized normal/all diagnostics and AST identities.
- Produces: every all-mode diagnostic classified exactly once and every category-C symbol furnished with its complete
  frozen `fail_open` record; exact totals 103/108 unless a separately reviewed semantic edit already changed them.

- [ ] **Step 1: Capture machine-readable discovery without writing policy**

Run: `ruff check --no-cache --select BLE001 --output-format json vntyper/ docker/app/ tests/ scripts/ docs/ > /tmp/ble001-normal.json || test $? -eq 1`

Run: `ruff check --no-cache --select BLE001 --ignore-noqa --output-format json vntyper/ docker/app/ tests/ scripts/ docs/ > /tmp/ble001-all.json || test $? -eq 1`

Expected: JSON arrays of 103 and 108 diagnostics. These temporary files are discovery evidence only and are never committed.

- [ ] **Step 2: Create the policy artifact from the fixed Appendix A worklist**

Set `reviewed_ruff_version` to `ruff 0.14.3`, expected counts to 103/108, and enter all 79 Appendix A rows with their
exact normal/all counts and planned category. In the same edit, add one complete `fail_open` record for every C row,
using the fixed Tasks 4–5 and Appendix B node IDs, exact outcomes, dispositions, rationales, and observability channels.
Do not add an empty, illustrative, or deferred record. Some linked nodes are intentionally created in Tasks 4–9, so
run the full policy CLI and require exit 1 with every unresolved node listed; Task 3 proves diagnostic
identities/counts without claiming those future nodes resolve.

- [ ] **Step 3: Verify and write Docker plus top-level-script rationales (2–5 min)**

Process Appendix A rows from `docker/app/archive_delivery.py` through `scripts/golden_cohort/runner.py` in table order.

- [ ] **Step 4: Verify and write test/benchmark rationales (2–5 min)**

Process every `tests/benchmark`, `tests/docker`, `tests/support`, and `tests/unit` Appendix A row.

- [ ] **Step 5: Verify and write CLI/module/alignment rationales (2–5 min)**

Process `vntyper/cli.py`, `vntyper/modules/advntr/advntr_genotyping.py`, and alignment/archive binding rows through
`vntyper/scripts/archive_safety.py`.

- [ ] **Step 6: Verify and write cohort rationales (2–5 min)**

Process `cohort_exports.py`, every `cohort_inputs.py` row, and every `cohort_summary.py` row.

- [ ] **Step 7: Verify and write flag/report/install rationales (2–5 min)**

Process `flagging.py`, `generate_report.py`, and `install_references.py` rows.

- [ ] **Step 8: Verify and write pipeline/reference rationales (2–5 min)**

Process `pipeline_guards.py`, `preflight_error_io.py`, and all reference-binding/cache rows.

- [ ] **Step 9: Verify and write region/report/screening rationales (2–5 min)**

Process `region_utils.py`, `report_formatting.py`, and `screening_summary.py` rows.

- [ ] **Step 10: Verify and write summary/utils/variant rationales (2–5 min)**

Process every `summary.py` row, both `utils.py` rows, and `variant_parsing.py`.

For Steps 3–10, inspect the handler and direct caller and use exactly one category-specific rationale form, prefixed
by the stable `path::symbol`: A = `terminal, re-raise, cleanup, or explicit error-state boundary; no ordinary success
value continues`; B = `test/benchmark-only audit candidate; production behavior is unaffected`; C =
`fallback/continuation boundary; the caller-visible outcome is pinned by <behavior_test> and failure remains
observable via <observable_via>`. The angle-bracket values are the actual manifest node/severity, not placeholders in
the committed JSON. If source contradicts a locked category, stop on that exact row. When a symbol contains multiple
broad handlers, C is locked if any handler continues with ordinary-looking output; counts still include all handlers.

- [ ] **Step 11: Add all existing explicit suppressions to the all-mode inventory**

Confirm the five-count delta maps to the current suppressed symbols in `docker/app/cohorts.py`, `scripts/golden_cohort/launcher.py`, and `tests/unit/test_pipeline_guards.py`. The exact live mapping, not this prose list, is authoritative; a moved suppression must appear under its new stable symbol.

- [ ] **Step 12: Add the live diagnostic-identity/count test**

```python
def test_live_ble001_diagnostics_match_reviewed_handler_counts() -> None:
    scope = read_ruff_paths(REPO_ROOT / "Makefile")
    normal = measure_ble001(REPO_ROOT, scope, ignore_noqa=False)
    all_handlers = measure_ble001(REPO_ROOT, scope, ignore_noqa=True)
    policy = load_policy(REPO_ROOT / "scripts/ble001_policy.json")
    assert len(normal.diagnostics) == policy.expected_normal
    assert len(all_handlers.diagnostics) == policy.expected_all
    expected_normal = {(row.path, row.symbol): row.normal_count for row in policy.handlers}
    expected_all = {(row.path, row.symbol): row.all_count for row in policy.handlers}

    def counts(measurement: Measurement) -> dict[tuple[str, str], int]:
        found: Counter[tuple[str, str]] = Counter()
        for diagnostic in measurement.diagnostics:
            source = (REPO_ROOT / diagnostic.path).read_text(encoding="utf-8")
            found[(diagnostic.path, enclosing_symbol(source, diagnostic.row))] += 1
        return dict(found)

    assert counts(normal) == {key: value for key, value in expected_normal.items() if value}
    assert counts(all_handlers) == {key: value for key, value in expected_all.items() if value}
    assert {(row.path, row.symbol) for row in policy.handlers if row.category == "C"} == {
        (row.path, row.symbol) for row in policy.fail_open
    }
```

Import `Counter` from `collections`. Add exact assertions that the current policy begins at 103/108, sums
`normal_count`/`all_count` to those values, and `expected_all - expected_normal == 5` before later reviewed deltas. This
test deliberately does not require future behavior nodes to resolve; `validate_policy` already resolves and reports
them, while the first passing full live gate remains Task 9.

- [ ] **Step 13: Run live discovery until the fixed inventory matches**

Run: `pytest -m unit tests/unit/test_ble001_policy.py::test_live_ble001_diagnostics_match_reviewed_handler_counts -vv`

Expected: PASS with all 108 diagnostics mapped exactly once, all 79 identities pinned at A=30/B=16/C=33, and exact
one-to-one C records already present. Then run the CLI and require exit 1 enumerating unresolved behavior nodes. A
missing/stale symbol is a source-drift failure to reconcile against Appendix A, not permission to invent a new row
during implementation. Do not claim behavior-node resolution yet.

- [ ] **Step 14: Update Ruff rationale accurately**

Replace the stale “46 handlers” comment with text stating: BLE001 remains deliberately unselected; the reviewed baseline is 103 normal/108 including five explicit suppressions under Ruff 0.14.3; `scripts/ble001_policy.json` and `tests/unit/test_ble001_policy.py` are the executable inventory; not every handler is a process boundary.

- [ ] **Step 15: Add configuration RED/GREEN assertions**

Parse only the relevant `[tool.ruff.lint]` and `[tool.ruff.lint.per-file-ignores]` text with a Python 3.10-compatible targeted parser using `re` plus `ast.literal_eval`; do not import `tomllib` or `tomli`, because neither exists in the Python 3.10 unit environment. Assert `BLE001` is absent from `select`, no per-file ignore contains it, and `G004` remains unselected and unchanged. Run the focused test.

- [ ] **Step 16: Commit the complete inventory atomically**

```bash
git add scripts/ble001_policy.py scripts/ble001_policy.json tests/unit/test_ble001_policy.py pyproject.toml
git commit -m "test(ci): freeze reviewed BLE001 inventory"
```

### Task 4: Freeze fail-open behavior for flagging and variant parsing

**Files:**
- Modify: `scripts/ble001_policy.json`
- Modify: `tests/unit/test_flagging.py`
- Modify: `tests/unit/test_variant_parsing.py`
- Modify production only if an existing contract requires conformance: `vntyper/scripts/flagging.py`, `vntyper/scripts/variant_parsing.py`

**Interfaces:**
- Consumes: `flagging.regex_match`, `flagging.evaluate_condition`, `variant_parsing.read_vcf_without_comments`.
- Produces: exact C-manifest records and behavior-test node IDs for those three seed symbols.

- [ ] **Step 1: Add unchanged-production log/outcome characterizations**

Extend existing tests so invalid regex returns `False` and logs ERROR containing `Error in regex_match`; missing name returns `False` and logs WARNING containing the missing name; a non-NameError evaluation failure returns `False` and logs ERROR. Use a fixed synthetic Series and expression; do not expand eval builtins.

```python
def test_invalid_regex_is_false_and_observable(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.ERROR):
        assert regex_match("[", "value") is False
    assert "Error in regex_match" in caplog.text
```

- [ ] **Step 2: Run focused behavior tests**

Run: `pytest -m unit tests/unit/test_flagging.py -k 'invalid_regex or invalid_pattern or missing_column or evaluation_error' -q`

Expected: PASS after log assertions; current fallback is preserved because repository/config documentation explicitly describes it.

Temporarily change the `regex_match` exception fallback from `False` to `True`, run
`tests/unit/test_flagging.py::test_invalid_regex_is_false_and_observable`, and observe failure on the returned value.
Restore the exact production line immediately and use `git diff` to prove the mutation is gone before touching policy.

- [ ] **Step 3: Strengthen VCF three-way distinction**

Keep the existing valid-empty, missing-file, and unreadable-path cases. Assert all return empty frames, but missing logs `VCF file not found`, unexpected failure logs `Error reading VCF file`, and a valid empty VCF emits neither error. This characterizes current behavior without inventing an exception outcome.

- [ ] **Step 4: Verify the frozen fail-open JSON records**

For each of the three symbols, verify the Task 3 record has the deterministic disposition, structured outcome (`False`
or empty `DataFrame`), domain rationale, exact now-existing test node ID, and log observability. The inventory category
must be C. Change a record only when the characterization proves the frozen fact is wrong and the existing contract
authorizes the corrected outcome.

- [ ] **Step 5: Run policy and behavior GREEN**

Run: `pytest -m unit tests/unit/test_flagging.py tests/unit/test_variant_parsing.py tests/unit/test_ble001_policy.py -q`

Expected: PASS; live counts remain unchanged unless a contract-backed handler edit was made.

- [ ] **Step 6: Commit**

```bash
git add scripts/ble001_policy.json tests/unit/test_flagging.py tests/unit/test_variant_parsing.py vntyper/scripts/flagging.py vntyper/scripts/variant_parsing.py
git commit -m "test(errors): characterize flagging and VCF fallbacks"
```

Omit unchanged production files.

### Task 5: Freeze guard, report, summary, and coverage-gate degradation

**Files:**
- Modify: `scripts/ble001_policy.json`
- Modify: `tests/unit/test_pipeline_guards.py`
- Modify: `tests/unit/test_region_utils.py`
- Modify: `tests/unit/test_report_formatting.py`
- Modify: `tests/unit/test_screening_summary.py`
- Modify: `tests/unit/test_coverage_gate.py`

**Interfaces:**
- Consumes: `read_alignment_header`, `_bam_contig_names`, `parse_coverage_stats`, `load_report_config`, `build_screening_summary`, `coverage_gate.read_total/main`.
- Produces: explicit C policies and caller-visible behavior assertions for every named seed contract.

- [ ] **Step 1: Pin alignment-header and BAM-contig fallback outcomes**

Assert `read_alignment_header` catches `CalledProcessError`, `FileNotFoundError`, and `OSError`, returns `None`, logs WARNING, and causes `enforce_declared_assembly` to return the documented undetermined verdict rather than agreement. Assert `_bam_contig_names` returns `[]` and logs WARNING on the same dependency failure, then assert its caller rejects only when actual contrary evidence exists.

- [ ] **Step 2: Pin the existing partial coverage-stat fallback**

Pass a row whose `mean` converts but whose `stdev` raises conversion. Assert `parse_coverage_stats` preserves the
converted `mean`, leaves `stdev` and every later initialized field as `None`, returns the complete `COVERAGE_COLUMNS`
mapping, and logs ERROR. This is the live contract pinned by
`tests/unit/test_report_formatting.py::test_parse_coverage_stats_survives_an_unreadable_value`; do not change
`report_formatting.py` merely to make the fallback all-or-nothing.

```python
def test_parse_coverage_stats_preserves_values_before_a_conversion_failure(caplog) -> None:
    row = {column: "3.5" for column in COVERAGE_COLUMNS}
    row["stdev"] = object()
    with caplog.at_level(logging.ERROR):
        result = parse_coverage_stats([row])
    assert result["mean"] == 3.5
    assert result["stdev"] is None
    assert all(result[column] is None for column in COVERAGE_COLUMNS[2:])
    assert "coverage" in caplog.text.lower()
```

- [ ] **Step 3: Pin report-config and screening unavailable behavior**

Patch `open`/`json.load` failure for `load_report_config`: assert `{}` plus ERROR. Patch one internal dependency of `build_screening_summary` to raise: assert text equals `UNAVAILABLE_SUMMARY_MESSAGE`, `is_positive is False`, both algorithm results are empty, `quality_metrics_pass is False`, `matched_rule is False`, and ERROR is logged. Do not assert or rewrite config-driven clinical sentences.

- [ ] **Step 4: Pin coverage-gate fail-closed caller behavior**

Patch `coverage.Coverage.load`/`report` to raise and assert `read_total() is None`; run `main()` and assert status 1 plus actionable `could not read` output. Record `read_total` as C internally but state in its rationale that the caller converts the sentinel into a non-zero gate, so CI does not pass open.

- [ ] **Step 5: Verify exact manifest records and run behavior/schema GREEN**

Verify the complete Task 3 C rows now point to the exact behavior node IDs and observability assertions. Before policy
validation, temporarily invert one asserted fallback field in its production return, run that one node to observe
failure, and restore the mutation. Then run:

`pytest -m unit tests/unit/test_pipeline_guards.py tests/unit/test_region_utils.py tests/unit/test_report_formatting.py tests/unit/test_screening_summary.py tests/unit/test_coverage_gate.py tests/unit/test_ble001_policy.py -q`

Expected: PASS. If an existing caller outcome differs, follow the current documented/tested outcome; do not create a new one.

- [ ] **Step 6: Commit**

```bash
git add scripts/ble001_policy.json tests/unit/test_pipeline_guards.py tests/unit/test_region_utils.py tests/unit/test_report_formatting.py tests/unit/test_screening_summary.py tests/unit/test_coverage_gate.py
git commit -m "test(errors): pin fail-open boundary outcomes"
```

### Task 6: Complete task, launcher, and adVNTR category-C evidence

**Files:**
- Modify: `scripts/ble001_policy.json`
- Modify: `tests/unit/web/test_result_expiry.py`
- Modify: `tests/unit/web/test_tasks.py`
- Modify: `tests/unit/test_golden_cohort_launcher.py` (created by quality-gates Task 5; preserve its import identity and
  `pytestmark = pytest.mark.unit`)
- Modify: `tests/unit/test_advntr_output_parsing.py`
- Modify production only for an exact contract-backed conformance change.

**Interfaces:**
- Consumes: the first three Appendix B rows.
- Produces: evidence for task cleanup, golden-cohort resolution, and adVNTR output parsing; every C policy record is
  already complete, while linked behavior nodes for Tasks 7–9 remain deliberately unresolved.

- [ ] **Step 1: Verify only the Task 6–9 behavior nodes are unresolved**

Run: `python scripts/ble001_policy.py --repo-root . --policy scripts/ble001_policy.json`

Expected: exit 1 listing unresolved behavior-test node IDs for the Task 6–9 Appendix B records. It must not report a
missing category-C record, empty rationale, or one-to-one C-manifest error. Any symbol outside that fixed node set is a
plan-review failure.

- [ ] **Step 2: Verify the locked Task 6 policy outcomes**

Inspect each named symbol's docstring, caller, config, and existing test; verify the Appendix B disposition/outcome and
verify the already-frozen record. If evidence contradicts the locked row, stop for plan review. No “needs review”, empty
rationale, or future-work marker is permitted.

- [ ] **Step 3: Write the four exact Task 6 dependency-exception nodes**

Use the exception type and assertion listed in Appendix B rather than applying one generic `RuntimeError` to unrelated
interfaces. Assert exact return/exit/state, logger severity, and whether caller continues or stops.
The new launcher test inserts `Path(__file__).resolve().parents[2] / "scripts"` at `sys.path[0]`, imports
`from golden_cohort import launcher`, and declares `pytestmark = pytest.mark.unit`; it never spawns a process.

For example, the cleanup continuation row is implemented literally as:

```python
def test_input_cleanup_logs_one_error_and_attempts_every_owned_path(monkeypatch, caplog) -> None:
    attempted: list[str] = []
    monkeypatch.setattr(tasks.os.path, "exists", lambda _path: True)

    def remove(path: str) -> None:
        attempted.append(path)
        if path.endswith("sample.bam"):
            raise OSError("blocked")

    monkeypatch.setattr(tasks.os, "remove", remove)
    tasks.remove_job_input_files("/jobs/1/sample.bam", "/jobs/1/sample.bam.bai")
    assert attempted == list(dict.fromkeys(("/jobs/1/sample.bam", "/jobs/1/sample.bam.bai",
                                             *tasks.derived_index_paths("/jobs/1/sample.bam"))))
    assert "blocked" in caplog.text
```

- [ ] **Step 4: Run each Task 6 node before editing production**

Run: `pytest -m unit tests/unit/web/test_result_expiry.py tests/unit/web/test_tasks.py tests/unit/test_golden_cohort_launcher.py tests/unit/test_advntr_output_parsing.py -x -vv`

Expected: PASS when characterizing current contract, or FAIL with the documented alternate outcome when a conformance change is authorized.

After the live characterization passes, temporarily change `launcher.resolve`'s structured fallback so `in_tree=True`,
run `test_resolve_returns_structured_import_failure`, observe failure, and restore the exact production line. Confirm
`git diff` contains no temporary mutation before policy validation.

- [ ] **Step 5: Make only minimal authorized GREEN changes**

For a conformance failure, change only that handler/caller to the already named outcome. Do not narrow its caught
exception tuple unless API documentation plus tests prove the complete raisable set. Update the already-frozen record
only if the authorized outcome or node identity changes.

- [ ] **Step 6: Validate Task 6 nodes and remaining-set accounting**

Run: `pytest -m unit tests/unit/test_ble001_policy.py -q`

Expected: PASS for schema and the now-existing Task 6 nodes. The CLI still exits 1 and names only unresolved Task 7–9
behavior node IDs; it reports no missing C record.

- [ ] **Step 7: Run focused affected tests and coverage**

Run: `pytest -m unit tests/unit -q && make test-unit-cov && make patch-coverage PATCH_COVERAGE_BASE=origin/main`

Expected: PASS at 86 repository and 80 patch floors. If patch coverage fails, add assertions for the changed branch; do not suppress or reclassify it merely to clear the score.

- [ ] **Step 8: Commit**

```bash
git add scripts/ble001_policy.json tests/unit/web/test_result_expiry.py tests/unit/web/test_tasks.py tests/unit/test_golden_cohort_launcher.py tests/unit/test_advntr_output_parsing.py
git commit -m "test(errors): pin task and launcher fallbacks"
```

Inspect the staged diff and exclude `scripts/ble001_policy.py` unless it intentionally changed.

### Task 7: Complete cohort category-C evidence

**Files:**
- Modify: `scripts/ble001_policy.json`
- Modify: `tests/unit/test_cohort_exports.py`
- Modify: `tests/unit/test_cohort_inputs.py`
- Modify: `tests/unit/test_cohort_summary_oracle.py`

**Interfaces:**
- Consumes: the cohort export/input/summary Appendix B rows.
- Produces: exact OSError/BadZipFile/ValueError fallback and continuation evidence for all ten cohort symbols.

- [ ] **Step 1: Add cohort export/input characterization nodes (2–5 min each case)**

Inject `OSError("blocked")` for table write, summary read, identity read, and cleanup; inject `zipfile.BadZipFile` for
archive discovery. Assert the exact Appendix B return/continuation, ERROR/WARNING log, and that later directories/inputs
are still attempted where named.

- [ ] **Step 2: Add cohort-summary characterization nodes (2–5 min each case)**

Inject `OSError` into image open, `ValueError` into `Axes.pie`, and `json.JSONDecodeError` into config load. Assert
`encode_image_to_base64` and `generate_donut_chart` return literal `""`; `load_report_config` returns `{}`; each logs
ERROR and no output file is claimed.

- [ ] **Step 3: Run the cohort nodes before manifest edits**

Run: `pytest -m unit tests/unit/test_cohort_exports.py tests/unit/test_cohort_inputs.py tests/unit/test_cohort_summary_oracle.py -q`

Expected: PASS characterization or the exact Appendix B conformance failure; no network/files outside `tmp_path`.

- [ ] **Step 4: Prove one cohort fallback assertion is mutation-sensitive**

Temporarily change `encode_image_to_base64`'s exception fallback from `""` to `"unavailable"`, run
`tests/unit/test_cohort_summary_oracle.py::test_image_encoding_failure_returns_empty_string`, and observe failure on the
exact value. Restore the production line and confirm `git diff` contains no temporary mutation.

- [ ] **Step 5: Verify the frozen manifest records**

Verify each Task 3 record uses disposition `preserved-contract`, the fixed Appendix B node ID/outcome, and asserted
logger severity. No production edit or record change is authorized unless Step 3 proves a documented contradiction.

- [ ] **Step 6: Validate remaining-node accounting**

Run: `pytest -m unit tests/unit/test_ble001_policy.py -q && python scripts/ble001_policy.py --repo-root . --policy scripts/ble001_policy.json`

Expected: unit schema and completed-node checks PASS; CLI exit 1 names only unresolved Task 8–9 behavior node IDs and
no missing C record.

- [ ] **Step 7: Commit**

```bash
git add scripts/ble001_policy.json tests/unit/test_cohort_exports.py tests/unit/test_cohort_inputs.py tests/unit/test_cohort_summary_oracle.py
git commit -m "test(errors): pin cohort fallback contracts"
```

### Task 8: Complete report and installer category-C evidence

**Files:**
- Modify: `scripts/ble001_policy.json`
- Modify: `tests/unit/test_generate_report.py`
- Modify: `tests/unit/test_install_references.py`

**Interfaces:**
- Consumes: generate-report and executable-probe Appendix B rows.
- Produces: exact typed report fallbacks and installer probe `False` evidence.

- [ ] **Step 1: Add report loader characterization nodes (2–5 min each case)**

Inject `OSError` into IGV, fastp, pipeline-log, and pipeline-summary reads. Assert respectively `("", "", "")`, `{}`,
`"Failed to load pipeline log."`, and `{}` plus the existing ERROR message.

- [ ] **Step 2: Add Kestrel conversion-continuation characterization node**

Make `pd.to_numeric` raise `ValueError("bad depth")`; assert both returned frames retain the input row, matching frame
remains unformatted, display frame still applies confidence formatting/escaping, and WARNING contains `bad depth`.
This is continuation with existing frames, not an empty-frame fallback.

- [ ] **Step 3: Add executable-probe characterization node**

Make the probe subprocess raise `OSError("missing loader")`; assert `check_executable_available(...) is False` and the
documented diagnostic logger/severity. Do not change install behavior.

- [ ] **Step 4: Run focused nodes before manifest edits**

Run: `pytest -m unit tests/unit/test_generate_report.py tests/unit/test_install_references.py -q`

Expected: PASS characterization or an exact documented conformance failure.

- [ ] **Step 5: Prove the IGV tuple assertion is mutation-sensitive**

Temporarily change one element of `extract_igv_content`'s exception tuple from `""` to `"missing"`, run
`tests/unit/test_generate_report.py::test_igv_extraction_failure_returns_empty_fragment`, and observe failure against
the exact `("", "", "")` assertion. Restore the production line and confirm `git diff` contains no temporary mutation.

- [ ] **Step 6: Verify manifest records and validate remaining nodes**

Verify disposition `preserved-contract` plus exact Appendix B nodes/outcomes/severities, then run
`pytest -m unit tests/unit/test_ble001_policy.py -q`. The policy CLI may still exit 1 only for unresolved Task 9 behavior
node IDs and must report no missing C record.

- [ ] **Step 7: Commit**

```bash
git add scripts/ble001_policy.json tests/unit/test_generate_report.py tests/unit/test_install_references.py
git commit -m "test(errors): pin report fallback contracts"
```

### Task 9: Complete summary and tool-version category-C evidence

**Files:**
- Modify: `scripts/ble001_policy.json`
- Modify: `tests/unit/test_summary_parsers.py`
- Modify: `tests/unit/test_summary_record_step.py`
- Modify: `tests/unit/test_utils.py`

**Interfaces:**
- Consumes: final summary/utils Appendix B rows.
- Produces: complete one-to-one category-C evidence with no missing/orphan node.

- [ ] **Step 1: Add parser/hash characterization nodes (2–5 min each case)**

Inject `OSError` for `md5sum`, CSV and TSV open; `json.JSONDecodeError` for JSON. Assert `None`, exact error comment,
exact error mapping, and exact error comment respectively, including whether ERROR is returned versus logged.

- [ ] **Step 2: Add record-step continuation characterization node**

Make the selected parser raise `ValueError("parser failed")`; assert one step is appended with the exact error payload,
the supplied step name/result path remain present, and later summary construction can continue.

- [ ] **Step 3: Add tool-version characterization node**

Make the version subprocess raise `OSError("exec failed")`; assert literal `"unknown"` and ERROR containing the tool
name plus failure.

- [ ] **Step 4: Run focused nodes before manifest edits**

Run: `pytest -m unit tests/unit/test_summary_parsers.py tests/unit/test_summary_record_step.py tests/unit/test_utils.py -q`

Expected: PASS characterization or the exact documented conformance failure.

- [ ] **Step 5: Prove the record-step continuation assertion is mutation-sensitive**

Temporarily prevent `record_step` from appending its error payload, run
`tests/unit/test_summary_record_step.py::test_parser_failure_is_recorded_and_step_is_appended`, and observe failure on
the exact appended step. Restore the production line and confirm `git diff` contains no temporary mutation.

- [ ] **Step 6: Verify all manifest records and add the full live validator gate**

Verify disposition `preserved-contract` and fixed Appendix B node/outcome/severity. Now add the deferred full gate:

```python
def test_live_ble001_inventory_and_behavior_nodes_match_reviewed_policy() -> None:
    scope = read_ruff_paths(REPO_ROOT / "Makefile")
    normal = measure_ble001(REPO_ROOT, scope, ignore_noqa=False)
    all_handlers = measure_ble001(REPO_ROOT, scope, ignore_noqa=True)
    policy = load_policy(REPO_ROOT / "scripts/ble001_policy.json")
    assert validate_policy(REPO_ROOT, normal, all_handlers, policy) == []
```

Run: `pytest -m unit tests/unit/test_ble001_policy.py -q && python scripts/ble001_policy.py --repo-root . --policy scripts/ble001_policy.json`.
Expected: both PASS; `set(C symbols) == set(fail_open manifest symbols)` and every nested/top-level node resolves. This
is the first **passing** full live validator gate because it is the first point at which every frozen behavior node
exists; earlier Tasks 3–8 runs intentionally exit 1 and enumerate unresolved nodes.

- [ ] **Step 7: Commit**

```bash
git add scripts/ble001_policy.json tests/unit/test_summary_parsers.py tests/unit/test_summary_record_step.py tests/unit/test_utils.py
git commit -m "test(errors): complete fallback evidence"
```

### Task 10: Prove drift detection, version attribution, and negative cases

**Files:**
- Modify: `tests/unit/test_ble001_policy.py`
- Modify if diagnostics are inadequate: `scripts/ble001_policy.py`

**Interfaces:**
- Consumes: complete live policy.
- Produces: objective mutation-style tests for unsuppressed/suppressed growth, symbol drift, metadata loss, and Ruff-version change.

- [ ] **Step 1: Add unsuppressed-growth negative test**

In a temporary repository, write a function with `except Exception`, run/fixture a normal and all measurement containing it, and omit it from policy. Assert validation errors include its relative path/symbol and both expected/actual counts.

```python
def test_validation_rejects_unreviewed_handler_growth(tmp_path: Path) -> None:
    source = tmp_path / "module.py"
    source.write_text("try:\n    pass\nexcept Exception:\n    pass\n", encoding="utf-8")
    diagnostic = Diagnostic("module.py", 3, 1, "BLE001", "blind-except")
    errors = validate_policy(
        tmp_path, Measurement("ruff 0.14.3", (diagnostic,)),
        Measurement("ruff 0.14.3", (diagnostic,)), _empty_policy(),
    )
    assert any("module.py::<module>" in error and "expected 0" in error for error in errors)
```

- [ ] **Step 2: Add suppressed-growth negative test**

Add `# noqa: BLE001` to the temporary handler. Normal measurement excludes it and all-mode includes it. Assert only all-count/classification errors fire and suppression delta changes from the policy value.

- [ ] **Step 3: Add symbol rename and line-motion tests**

Assert adding blank lines leaves validation green; renaming `Worker.run` to `Worker.execute` produces removed `Worker.run` and added `Worker.execute`, never a permanent row-number identity.

- [ ] **Step 4: Add manifest-field negative tests**

Delete each C field in turn (`disposition`, `outcome`, `rationale`, `behavior_test`, `observable_via`) and assert `load_policy` names the exact field. Change a test node ID to a nonexistent path/function and assert validation fails.

- [ ] **Step 5: Add reviewed/actual version attribution tests**

Use identical diagnostics with reviewed `ruff 0.14.3` and actual `ruff 0.16.0`: expect success plus displayed attribution. Change one identity under actual 0.16.0: expect failure containing both versions and an instruction to reclassify, not an unqualified “source growth” claim.

- [ ] **Step 6: Run GREEN and mutation check**

Run: `pytest -m unit tests/unit/test_ble001_policy.py -q`

Expected: PASS.

Temporarily add one local blind handler and run the live node; expected FAIL. Revert that temporary edit immediately and rerun; expected PASS. Do not commit the mutation.

- [ ] **Step 7: Verify typing and patch coverage**

Run: `mypy scripts/ble001_policy.py tests/unit/test_ble001_policy.py && make test-unit-cov && make patch-coverage PATCH_COVERAGE_BASE=origin/main`

Expected: PASS.

- [ ] **Step 8: Commit**

```bash
git add scripts/ble001_policy.py tests/unit/test_ble001_policy.py
git commit -m "test(ci): guard BLE001 policy drift"
```

### Task 11: Final dual remeasurement and exact policy update

**Files:**
- Modify: `scripts/ble001_policy.json`
- Modify: `pyproject.toml:249-254`
- Modify only if final literal guards live there: `tests/unit/test_ble001_policy.py`

**Interfaces:**
- Consumes: final source after all reviewed behavior edits.
- Produces: final exact normal/all counts and actual-version attribution; offsetting additions/removals remain explicit inventory changes.

- [ ] **Step 1: Run final no-cache normal measurement**

Run: `ruff --version && ruff check --no-cache --select BLE001 --statistics vntyper/ docker/app/ tests/ scripts/ docs/`

Expected: Ruff version printed and a deterministic exact normal count. It may differ from 103 only by named reviewed edits.

- [ ] **Step 2: Run final no-cache all-mode measurement**

Run: `ruff check --no-cache --select BLE001 --ignore-noqa --statistics vntyper/ docker/app/ tests/ scripts/ docs/`

Expected: deterministic exact all count, never below normal. Its delta equals live explicit suppressions.

- [ ] **Step 3: Reconcile every delta rather than merely changing totals**

Run the policy CLI. For each added/removed stable identity, confirm the corresponding source commit and update handler/C manifest records. If counts changed with no intended source explanation, restore the accidental edit. An addition offset by a removal still requires both inventory changes.

- [ ] **Step 4: Pin final exact counts and rationale**

Set `expected_counts.normal` and `.ignore_noqa` to the two exact results. Update the pyproject comment to the same pair, reviewed Ruff `0.14.3`, actual version if different, and executable policy/test paths. Do not use raw grep counts.

```python
def test_policy_counts_equal_live_measurements() -> None:
    policy = load_policy(REPO_ROOT / "scripts/ble001_policy.json")
    normal, all_handlers = _live_measurements()
    assert (policy.expected_normal, policy.expected_all) == (
        len(normal.diagnostics), len(all_handlers.diagnostics)
    )
    assert normal.ruff_version == all_handlers.ruff_version
```

- [ ] **Step 5: Run final live equality guard twice**

Run: `python scripts/ble001_policy.py --repo-root . --policy scripts/ble001_policy.json && pytest -m unit tests/unit/test_ble001_policy.py -q`

Expected: PASS; output includes reviewed version, actual version, exact normal/all pair, suppression delta, and A/B/C totals.

- [ ] **Step 6: Commit**

```bash
git add scripts/ble001_policy.json pyproject.toml tests/unit/test_ble001_policy.py
git commit -m "docs(ci): record final BLE001 policy counts"
```

### Task 12: Full verification, partial-failure recovery, and rollout

**Files:**
- Modify contributor docs only if they still quote 46/67 or imply every broad handler is a boundary.
- Verify: `mkdocs.yml`, `AGENTS.md`, `.github/workflows/ci-tests.yml`.

**Interfaces:**
- Consumes: final policy, quality gates from #211, all behavior tests.
- Produces: merge-ready #219 with no runtime migration or new permissions.

- [ ] **Step 1: Scan for stale current-count claims and forbidden policy changes**

Run: `rg -n 'BLE001|blind-except|46 handlers|67' pyproject.toml AGENTS.md docs scripts tests`

Expected: historical counts are clearly labeled historical or removed; current claims match policy. Confirm no global BLE001 selection or broad per-file ignore.

Add this literal source guard before editing contributor prose:

```python
def test_contributor_exception_counts_are_generated_from_policy() -> None:
    agents = (REPO_ROOT / "AGENTS.md").read_text(encoding="utf-8")
    policy = load_policy(REPO_ROOT / "scripts/ble001_policy.json")
    assert f"{policy.expected_normal} normal/{policy.expected_all} including suppressions" in agents
    assert "current BLE001 count is 67" not in agents
```

- [ ] **Step 2: Run focused policy and behavior suite**

Run: `pytest -m unit tests/unit/test_ble001_policy.py tests/unit/test_flagging.py tests/unit/test_variant_parsing.py tests/unit/test_pipeline_guards.py tests/unit/test_region_utils.py tests/unit/test_report_formatting.py tests/unit/test_screening_summary.py tests/unit/test_coverage_gate.py -q`

Expected: PASS.

- [ ] **Step 3: Run repository gates**

Run: `make format-check && make lint && make type-check-all && make test-unit-cov && make patch-coverage PATCH_COVERAGE_BASE=origin/main`

Expected: PASS at unchanged 86 combined and 80 patch floors; `scripts/ble001_policy.py` is measured.

- [ ] **Step 4: Run full local and docs gates**

Run: `make check-all && make docs-build`

Expected: PASS; `docs/superpowers/` is excluded from publication/template processing.

- [ ] **Step 5: Run CI-equivalent environment verification when workflow content changed**

Run: `make ci-local`

Expected: PASS using an explicit uv venv; no workflow permission or secret change.

- [ ] **Step 6: Recover partial failures deterministically**

If normal changes but all-mode does not, inspect `noqa` changes. If both change, inspect stable-identity diff. If a symbol is unmapped, classify it as `<module>`. If a behavior test needs external state, mock/extract its effect boundary. If narrowing causes an unhandled crash or changed caller result, revert the narrowing but retain policy/characterization. After each fix rerun the focused failure, both Ruff measurements, then the full failed gate.

- [ ] **Step 7: Inspect final state**

Run: `git diff --check && git status --short && git log --oneline origin/main..HEAD`

Expected: only intentional policy/source/test/docs changes; no Ruff JSON captures, coverage data, test data, secrets, or generated artifacts.

- [ ] **Step 8: Final Conventional Commit if documentation changed after remeasurement**

```bash
git add AGENTS.md pyproject.toml scripts/ble001_policy.json
git commit -m "docs(errors): document reviewed exception policy" -m "Closes #219"
```

Skip this commit when Task 8 already captured all final documentation and the worktree is clean.

## Appendix A: fixed 108-diagnostic inventory worklist

This worklist is the base-SHA AST aggregation used by Task 3. `N/A` is normal/all-mode count. Categories are locked for
execution: A is an explicit terminal, re-raise, cleanup, or error-state boundary; B is a test/benchmark handler left as
an enumerable narrowing candidate; C returns an ordinary-looking fallback or continues after degradation. Source drift
must be reviewed as a plan change rather than silently absorbed.

| Path | Qualified symbol | N/A | Category |
| --- | --- | ---: | :---: |
| `docker/app/archive_delivery.py` | `OwnedFileResponse.__call__` | 2/2 | A |
| `docker/app/archive_delivery.py` | `_close_snapshot_directory` | 3/3 | A |
| `docker/app/archive_delivery.py` | `snapshot_owned_archives` | 1/1 | A |
| `docker/app/cohorts.py` | `_discard_partial_cohort` | 0/2 | A |
| `docker/app/tasks.py` | `delete_old_results` | 2/2 | C |
| `docker/app/tasks.py` | `remove_job_input_files` | 1/1 | C |
| `docker/app/tasks.py` | `run_cohort_analysis_job` | 7/7 | A |
| `docker/app/tasks.py` | `run_vntyper_job` | 7/7 | A |
| `scripts/coverage_gate.py` | `read_total` | 1/1 | C |
| `scripts/download_test_data.py` | `main` | 1/1 | A |
| `scripts/golden_cohort/launcher.py` | `launch` | 0/1 | A |
| `scripts/golden_cohort/launcher.py` | `resolve` | 0/1 | C |
| `scripts/golden_cohort/runner.py` | `_run_one` | 1/1 | A |
| `tests/benchmark/benchamrk_downsample.py` | `parse_pipeline_log` | 1/1 | B |
| `tests/benchmark/benchamrk_downsample.py` | `summarize_advntr_results` | 1/1 | B |
| `tests/benchmark/benchamrk_downsample.py` | `summarize_vntr_results` | 1/1 | B |
| `tests/benchmark/benchmark_vntyper.py` | `main` | 1/1 | B |
| `tests/benchmark/benchmark_vntyper.py` | `summarize_vntyper_results` | 1/1 | B |
| `tests/benchmark/plot_vntyper_summary.py` | `get_advntr_call` | 1/1 | B |
| `tests/benchmark/plot_vntyper_summary.py` | `main` | 1/1 | B |
| `tests/benchmark/plot_vntyper_summary.py` | `parse_advntr_result` | 1/1 | B |
| `tests/docker/image_probe.py` | `_run` | 1/1 | B |
| `tests/support/data_utils.py` | `ensure_test_data_downloaded` | 1/1 | B |
| `tests/unit/test_pipeline_guards.py` | `test_the_pure_verdict_function_still_has_no_policy` | 0/1 | B |
| `tests/unit/test_preflight_input_reads.py` | `_bed_read_worker` | 1/1 | B |
| `tests/unit/test_preflight_input_reads.py` | `_reference_reason_worker` | 1/1 | B |
| `tests/unit/test_ref_path_is_pinned.py` | `test_an_invalid_scope_does_not_block_the_next_thread.valid_scope` | 1/1 | B |
| `tests/unit/test_ref_path_is_pinned.py` | `test_overlapping_threads_cannot_reopen_ambient_reference_resolution.first_scope` | 1/1 | B |
| `tests/unit/test_ref_path_is_pinned.py` | `test_overlapping_threads_cannot_reopen_ambient_reference_resolution.second_scope` | 1/1 | B |
| `vntyper/cli.py` | `load_config` | 1/1 | A |
| `vntyper/cli.py` | `main` | 1/1 | A |
| `vntyper/modules/advntr/advntr_genotyping.py` | `process_advntr_output` | 2/2 | C |
| `vntyper/modules/advntr/advntr_genotyping.py` | `run_advntr` | 2/2 | A |
| `vntyper/scripts/alignment_binding.py` | `AlignmentBinding.__del__` | 1/1 | A |
| `vntyper/scripts/alignment_contract.py` | `AlignmentPlan.close` | 2/2 | A |
| `vntyper/scripts/alignment_index_binding.py` | `AlignmentIndexBinding.__del__` | 1/1 | A |
| `vntyper/scripts/archive_safety.py` | `create_safe_archive` | 3/3 | A |
| `vntyper/scripts/cohort_exports.py` | `write_pseudonymization_table` | 1/1 | C |
| `vntyper/scripts/cohort_inputs.py` | `_identity_from_summary` | 1/1 | C |
| `vntyper/scripts/cohort_inputs.py` | `cleanup_temp_dirs` | 1/1 | C |
| `vntyper/scripts/cohort_inputs.py` | `discover_sample_directories` | 1/1 | C |
| `vntyper/scripts/cohort_inputs.py` | `load_pipeline_summary_for_sample` | 1/1 | C |
| `vntyper/scripts/cohort_summary.py` | `encode_image_to_base64` | 1/1 | C |
| `vntyper/scripts/cohort_summary.py` | `generate_donut_chart` | 1/1 | C |
| `vntyper/scripts/cohort_summary.py` | `load_report_config` | 1/1 | C |
| `vntyper/scripts/flagging.py` | `evaluate_condition` | 1/1 | C |
| `vntyper/scripts/flagging.py` | `regex_match` | 1/1 | C |
| `vntyper/scripts/generate_report.py` | `build_kestrel_frames` | 1/1 | C |
| `vntyper/scripts/generate_report.py` | `extract_igv_content` | 1/1 | C |
| `vntyper/scripts/generate_report.py` | `load_fastp_output` | 1/1 | C |
| `vntyper/scripts/generate_report.py` | `load_pipeline_log` | 1/1 | C |
| `vntyper/scripts/generate_report.py` | `load_pipeline_summary` | 1/1 | C |
| `vntyper/scripts/install_references.py` | `calculate_md5` | 1/1 | A |
| `vntyper/scripts/install_references.py` | `check_executable_available` | 1/1 | C |
| `vntyper/scripts/install_references.py` | `download_file` | 1/1 | A |
| `vntyper/scripts/install_references.py` | `load_install_config` | 1/1 | A |
| `vntyper/scripts/install_references.py` | `main` | 1/1 | A |
| `vntyper/scripts/install_references.py` | `process_ucsc_references` | 3/3 | A |
| `vntyper/scripts/install_references.py` | `process_vntyper_references` | 2/2 | A |
| `vntyper/scripts/install_references.py` | `update_config` | 2/2 | A |
| `vntyper/scripts/install_references.py` | `write_md5_checksums` | 1/1 | A |
| `vntyper/scripts/pipeline_guards.py` | `read_alignment_header` | 1/1 | C |
| `vntyper/scripts/preflight_error_io.py` | `persist_preflight_failure` | 1/1 | A |
| `vntyper/scripts/reference_binding.py` | `ReferenceBinding.__del__` | 1/1 | A |
| `vntyper/scripts/reference_binding.py` | `ReferenceBinding.__init__` | 1/1 | A |
| `vntyper/scripts/reference_binding.py` | `ReferenceBinding.close` | 3/3 | A |
| `vntyper/scripts/reference_cache_binding.py` | `PrivateReferenceCache.close` | 2/2 | A |
| `vntyper/scripts/region_utils.py` | `_bam_contig_names` | 1/1 | C |
| `vntyper/scripts/report_formatting.py` | `parse_coverage_stats` | 1/1 | C |
| `vntyper/scripts/screening_summary.py` | `build_screening_summary` | 1/1 | C |
| `vntyper/scripts/screening_summary.py` | `load_report_config` | 1/1 | C |
| `vntyper/scripts/summary.py` | `md5sum` | 1/1 | C |
| `vntyper/scripts/summary.py` | `parse_csv` | 1/1 | C |
| `vntyper/scripts/summary.py` | `parse_json_file` | 1/1 | C |
| `vntyper/scripts/summary.py` | `parse_tsv` | 1/1 | C |
| `vntyper/scripts/summary.py` | `record_step` | 1/1 | C |
| `vntyper/scripts/utils.py` | `get_tool_version` | 1/1 | C |
| `vntyper/scripts/utils.py` | `load_config` | 1/1 | A |
| `vntyper/scripts/variant_parsing.py` | `read_vcf_without_comments` | 1/1 | C |

## Appendix B: fixed additional category-C behavior worklist

Tasks 4-5 own flagging, VCF, guard, region, report-formatting, screening, and coverage-gate rows. Task 6 owns the rows
below. Every node is added to the named existing unit file before its manifest record; test names are fixed so the JSON
can refer to stable node IDs.

| Symbols | Injected exception(s) | Disposition | Test node(s) and exact outcome |
| --- | --- | --- | --- |
| `docker/app/tasks.py::{delete_old_results,remove_job_input_files}` | `OSError("blocked")` from one owned-path deletion | `preserved-contract` | `tests/unit/web/test_result_expiry.py::test_cleanup_continues_after_one_delete_error` and `tests/unit/web/test_tasks.py::test_input_cleanup_logs_one_error_and_attempts_every_owned_path`: ERROR is logged and remaining owned paths are attempted; no exception replaces the task result. |
| `scripts/golden_cohort/launcher.py::resolve` | `ImportError("candidate import failed")` | `preserved-contract` | `tests/unit/test_golden_cohort_launcher.py::test_resolve_returns_structured_import_failure`: `error` names the exception, `in_tree=False`, and `vntyper_file=None`. |
| `vntyper/modules/advntr/advntr_genotyping.py::process_advntr_output` | `OSError("unreadable result")` | `preserved-contract` | `tests/unit/test_advntr_output_parsing.py::test_unreadable_advntr_output_logs_and_returns_none`: returns `None`, logs ERROR, writes no result. |
| `vntyper/scripts/cohort_exports.py::write_pseudonymization_table` | `OSError("blocked")` | `preserved-contract` | `tests/unit/test_cohort_exports.py::test_pseudonym_table_write_failure_is_logged`: returns `None`, logs ERROR, and leaves no claimed successful table. |
| `vntyper/scripts/cohort_inputs.py::{_identity_from_summary,cleanup_temp_dirs,discover_sample_directories,load_pipeline_summary_for_sample}` | `OSError` for reads/cleanup and `zipfile.BadZipFile` for archive discovery | `preserved-contract` | `tests/unit/test_cohort_inputs.py::{test_identity_read_failure_uses_directory_fallback,test_cleanup_attempts_all_directories_after_failure,test_bad_archive_is_skipped_and_other_inputs_continue,test_summary_read_failure_returns_three_empty_results}` with the literal fallback/continuation named by each test. |
| `vntyper/scripts/cohort_summary.py::{encode_image_to_base64,generate_donut_chart,load_report_config}` | `OSError`, `ValueError`, `json.JSONDecodeError` | `preserved-contract` | `tests/unit/test_cohort_summary_oracle.py::{test_image_encoding_failure_returns_empty_string,test_donut_failure_returns_empty_string,test_report_config_failure_returns_empty_mapping}`; each logs ERROR and returns respectively `""`, `""`, and `{}`. |
| `vntyper/scripts/generate_report.py::{build_kestrel_frames,extract_igv_content,load_fastp_output,load_pipeline_log,load_pipeline_summary}` | `ValueError` for numeric conversion; `OSError` for reads | `preserved-contract` | `tests/unit/test_generate_report.py::{test_kestrel_conversion_failure_preserves_both_frames,test_igv_extraction_failure_returns_empty_fragment,test_fastp_failure_returns_empty_mapping,test_pipeline_log_failure_returns_failure_message,test_pipeline_summary_failure_returns_empty_mapping}`; outcomes are retained display/matching frames, exact tuple `("", "", "")`, `{}`, literal `"Failed to load pipeline log."`, and `{}` with the existing log severity. |
| `vntyper/scripts/install_references.py::check_executable_available` | `OSError("missing loader")` | `preserved-contract` | `tests/unit/test_install_references.py::test_executable_probe_error_returns_false`: returns `False` and logs the live diagnostic severity. |
| `vntyper/scripts/summary.py::{md5sum,parse_csv,parse_json_file,parse_tsv,record_step}` | `OSError`, `json.JSONDecodeError`, `ValueError("parser failed")` | `preserved-contract` | `tests/unit/test_summary_parsers.py::{test_md5_failure_returns_none,test_csv_failure_returns_error_comment,test_json_failure_returns_error_mapping,test_tsv_failure_returns_error_comment}` and `tests/unit/test_summary_record_step.py::test_parser_failure_is_recorded_and_step_is_appended`; assert the exact structured fallback in each node. |
| `vntyper/scripts/utils.py::get_tool_version` | `OSError("exec failed")` | `preserved-contract` | `tests/unit/test_utils.py::test_tool_version_unexpected_failure_returns_unknown`: returns literal `"unknown"` and logs ERROR. |
