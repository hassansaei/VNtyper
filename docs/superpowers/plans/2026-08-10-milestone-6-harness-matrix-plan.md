# Milestone 6 Harness Matrix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete #71 by proving intermediate/archive presence and absence, correcting delete-overrides-keep behavior, and running an independent paired b178 FASTQ case without SHARK, while treating #71 criteria 1–4 and #226 as evidence-only.

**Architecture:** `tests/test_data_config.json` remains the declarative case source. Pure helpers in `tests/support/orchestration.py` validate confined output-relative artifacts and optional sibling ZIP state; successful single-end cases exercise the matrix, while `process_bam_to_fastq` owns the production precedence fix. The cleanup artifact is always `fastq_bam_processing/output_unmapped.bam` because `pipeline.py` passes `output_name="output"`.

**Tech Stack:** Python 3.10–3.13, pytest, `pathlib.Path`, JSON parametrization, VNtyper CLI integration fixtures, Ruff, mypy, coverage.py and diff-cover.

## Global Constraints

- Follow `docs/superpowers/specs/2026-08-10-milestone-6-harness-matrix-design.md` exactly.
- #71 criteria 1–4 are satisfied evidence: do not add Docker, API, adVNTR or subcommand scope.
- #226 is evidence-only, already implemented by `80f42fbba63c5c045ee605e49949cdd970884320` and `edabd3eaf594b785906d8ba03c9cce60f1c6babd`; do not modify its code or tests.
- Use the literal `fastq_bam_processing/output_unmapped.bam`; never interpolate a case name.
- Establish the keep-row presence proof before adding either absence row.
- `--delete-intermediates` overrides `--keep-intermediates`; record this above 2.0.10 under `Unreleased`.
- The independent FASTQ case is paired b178 R1/R2 and has no SHARK option.
- The issue comment claiming exactly one existing FASTQ case is stale: the current
  `example_6449_hg19_subset_single_fastq` already omits SHARK. AC71-6 still needs an independent *paired b178* case;
  Task 6 adds that missing specimen rather than claiming no non-SHARK case exists.
- Unit tests are pure, run from repository root and carry the `unit` marker.
- Preserve Python 3.10 compatibility, Ruff's 120-column limit and Google-style public docstrings.
- Keep the branch-inclusive floor and target at 86 and patch coverage at 80.
- Treat each checkbox as one 2–5-minute edit, assertion addition, review, or command launch; record its result before starting the next checkbox.

---

## File map and evidence-only criteria

| File | Responsibility |
| --- | --- |
| `tests/support/orchestration.py` | Validate/assert output-relative artifacts and optional sibling ZIP state. |
| `tests/unit/test_integration_outcome_contract.py` | Pure declaration, diagnostic, hygiene and exact-case tests. |
| `tests/unit/test_fastq_bam_command_wiring.py` | Four-combination cleanup precedence regression. |
| `vntyper/scripts/fastq_bam_processing.py` | Production delete-overrides-keep decision. |
| `tests/test_data_config.json` | Single-end matrix and independent paired b178 case. |
| `tests/integration/test_single_end_pipeline.py` | Successful matrix runner and archive assertion. |
| `tests/integration/test_pipeline_integration.py` | Shared archive assertion for ordinary BAM cases and FASTQ runner. |
| `docs/about/changelog.md` | User-visible precedence correction. |

AC71-1 remains covered by `tests/docker/test_docker_pipeline.py` and
`tests/docker/test_image_structure.py`; AC71-2 by `tests/unit/web/`; AC71-3 by the
successful `example_b178_hg19_bwa_advntr` case; and AC71-4 by the parser/handler/dispatch,
report and online unit suites. Close #226 separately with
`build_reference_dependent_fixture`, its two commits,
`test_reference_dependent_fixture_has_a_local_ur_target_that_can_be_removed`, and the
A209-1/A209-2/A209-3 integration tests. No task below changes those surfaces.

### Task 1: Confined artifact declarations

**Files:**
- Modify: `tests/unit/test_integration_outcome_contract.py`
- Modify: `tests/support/orchestration.py:35-130`

**Interfaces:**
- Consumes: optional `expected_present: list[str]` and `expected_absent: list[str]`.
- Produces: `assert_declared_artifacts(test_case: dict, output_dir: Path) -> None`.

- [ ] **Step 1: Write the RED presence/absence test**

```python
def test_declared_artifacts_assert_both_presence_and_absence(tmp_path: Path) -> None:
    kept = tmp_path / "kept.txt"
    kept.write_text("kept\n", encoding="utf-8")
    case = {
        "test_name": "artifact-case",
        "expected_present": ["kept.txt"],
        "expected_absent": ["removed.txt"],
    }
    orchestration.assert_declared_artifacts(case, tmp_path)

    kept.unlink()
    (tmp_path / "removed.txt").write_text("unexpected\n", encoding="utf-8")
    with pytest.raises(AssertionError) as exc_info:
        orchestration.assert_declared_artifacts(case, tmp_path)
    assert f"case=artifact-case field=expected_present missing: {kept.resolve()}" in str(exc_info.value)
    assert f"case=artifact-case field=expected_absent unexpectedly present: {(tmp_path / 'removed.txt').resolve()}" in str(
        exc_info.value
    )
```

- [ ] **Step 2: Write RED validation cases**

```python
@pytest.mark.parametrize(
    "case",
    [
        {"expected_present": [""]},
        {"expected_present": ["/absolute.txt"]},
        {"expected_present": ["../escape.txt"]},
        {"expected_present": ["nested/../../escape.txt"]},
        {"expected_present": ["same.txt", "same.txt"]},
        {"expected_present": ["same.txt"], "expected_absent": ["same.txt"]},
    ],
)
def test_declared_artifacts_reject_invalid_paths(tmp_path: Path, case: dict) -> None:
    with pytest.raises(ValueError, match="artifact|expected_"):
        orchestration.assert_declared_artifacts(case, tmp_path)
```

- [ ] **Step 3: Add RED symlink-confinement fixtures**

Add all filesystem-escape cases before the helper exists. The intermediate-directory case is distinct from a final
symlink: resolving only the leaf would miss it.

```python
def test_declared_artifacts_reject_leaf_and_intermediate_symlink_escapes(tmp_path: Path) -> None:
    output = tmp_path / "output"
    output.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "leaf.txt").write_text("outside\n", encoding="utf-8")
    (output / "leaf.txt").symlink_to(outside / "leaf.txt")
    (output / "linked-dir").symlink_to(outside, target_is_directory=True)

    for declared in ("leaf.txt", "linked-dir/leaf.txt"):
        with pytest.raises(ValueError, match="escapes output_dir"):
            orchestration.assert_declared_artifacts({"expected_present": [declared]}, output)


def test_declared_absence_rejects_a_dangling_symlink_entry(tmp_path: Path) -> None:
    output = tmp_path / "output"
    output.mkdir()
    (output / "broken.txt").symlink_to(tmp_path / "missing.txt")
    case = {"test_name": "broken-link", "expected_absent": ["broken.txt"]}
    with pytest.raises(AssertionError, match="case=broken-link field=expected_absent unexpectedly present"):
        orchestration.assert_declared_artifacts(case, output)
```

- [ ] **Step 4: Run RED**

Run: `pytest -m unit tests/unit/test_integration_outcome_contract.py -k declared_artifacts -v`

Expected: FAIL because `assert_declared_artifacts` does not exist; none of the traversal or symlink fixtures may pass
through a partially implemented helper.

- [ ] **Step 5: Implement the minimal helper**

Add `import os` beside the existing stdlib imports, then implement:

```python
def _declared_artifact_paths(test_case: dict, field: str, output_dir: Path) -> list[Path]:
    test_name = str(test_case.get("test_name", "<unnamed>"))
    raw_paths = test_case.get(field, [])
    if not isinstance(raw_paths, list):
        raise ValueError(f"case={test_name} field={field} must be a list of relative artifact paths")
    root = output_dir.resolve()
    result: list[Path] = []
    seen: set[str] = set()
    for raw in raw_paths:
        if not isinstance(raw, str) or not raw:
            raise ValueError(f"case={test_name} field={field} contains an empty or non-string artifact path")
        relative = Path(raw)
        if relative.is_absolute() or ".." in relative.parts or raw in seen:
            raise ValueError(f"case={test_name} field={field} invalid artifact: {raw}")
        seen.add(raw)
        candidate = root / relative
        resolved_parent = candidate.parent.resolve(strict=False)
        if not resolved_parent.is_relative_to(root):
            raise ValueError(f"case={test_name} field={field} artifact escapes output_dir: {raw}")
        if field == "expected_present" and not candidate.resolve(strict=False).is_relative_to(root):
            raise ValueError(f"case={test_name} field={field} artifact escapes output_dir: {raw}")
        result.append(candidate)
    return result


def assert_declared_artifacts(test_case: dict, output_dir: Path) -> None:
    """Assert every case-declared output-relative artifact state."""
    test_name = str(test_case.get("test_name", "<unnamed>"))
    present = _declared_artifact_paths(test_case, "expected_present", output_dir)
    absent = _declared_artifact_paths(test_case, "expected_absent", output_dir)
    overlap = set(present).intersection(absent)
    if overlap:
        raise ValueError(f"case={test_name} fields=expected_present,expected_absent overlap: {sorted(map(str, overlap))}")
    failures = [
        f"case={test_name} field=expected_present missing: {path}" for path in present if not path.exists()
    ]
    failures.extend(
        f"case={test_name} field=expected_absent unexpectedly present: {path}"
        for path in absent
        if os.path.lexists(path)
    )
    assert not failures, "Declared artifact mismatch:\n" + "\n".join(failures)
```

- [ ] **Step 6: Refactor helper names and documentation**

Add Google-style docstrings to both helpers and keep path parsing separate from
assertion aggregation. Do not change the tested messages.

- [ ] **Step 7: Run focused GREEN**

Run: `pytest -m unit tests/unit/test_integration_outcome_contract.py -v`

Expected: PASS, including direct and nested parent traversal, leaf and intermediate symlink escape, and dangling
`expected_absent` entry cases. Absence uses `os.path.lexists`, not `Path.exists`.

- [ ] **Step 8: Commit**

```bash
git add tests/support/orchestration.py tests/unit/test_integration_outcome_contract.py
git commit -m "test(harness): validate declared pipeline artifacts"
```

### Task 2: Successful-case and archive wiring

**Files:**
- Modify: `tests/unit/test_integration_outcome_contract.py`
- Modify: `tests/support/orchestration.py:69-130`
- Modify: `tests/integration/test_pipeline_integration.py:147-201`
- Modify: `tests/integration/test_single_end_pipeline.py:24-63`

**Interfaces:**
- Consumes: `assert_declared_artifacts` from Task 1.
- Produces: `assert_declared_archive(test_case: dict, output_dir: Path) -> None`.

- [ ] **Step 1: Write RED orchestration and archive tests**

```python
def test_successful_bam_enforces_declared_artifacts(tmp_path: Path) -> None:
    case = {
        "test_name": "bam-artifact",
        "bam": "clean.bam",
        "reference_assembly": "hg19",
        "kestrel_assertions": {},
        "expected_present": ["must-exist.txt"],
    }
    with (
        mock.patch.object(orchestration, "assert_required_files"),
        mock.patch.object(orchestration, "validate_kestrel_output"),
        mock.patch.object(orchestration, "validate_coverage_output", return_value={"mean_cov": 0}),
    ):
        with pytest.raises(AssertionError, match="case=bam-artifact field=expected_present"):
            orchestration.run_bam_test_case(case, mock.Mock(return_value=0), tmp_path)


def test_declared_archive_distinguishes_omitted_false_and_true(tmp_path: Path) -> None:
    output = tmp_path / "result"
    output.mkdir()
    archive = Path(f"{output}.zip")
    orchestration.assert_declared_archive({"test_name": "omitted"}, output)
    orchestration.assert_declared_archive({"test_name": "absent", "expected_archive": False}, output)
    with pytest.raises(AssertionError, match="case=present field=expected_archive"):
        orchestration.assert_declared_archive({"test_name": "present", "expected_archive": True}, output)
    archive.write_bytes(b"zip")
    orchestration.assert_declared_archive({"test_name": "present", "expected_archive": True}, output)
    with pytest.raises(AssertionError, match="case=absent field=expected_archive"):
        orchestration.assert_declared_archive({"test_name": "absent", "expected_archive": False}, output)


def test_declared_archive_rejects_invalid_boolean_and_broken_symlink(tmp_path: Path) -> None:
    output = tmp_path / "result"
    output.mkdir()
    with pytest.raises(ValueError, match="case=invalid field=expected_archive"):
        orchestration.assert_declared_archive({"test_name": "invalid", "expected_archive": "false"}, output)
    Path(f"{output}.zip").symlink_to(tmp_path / "missing.zip")
    with pytest.raises(AssertionError, match="case=absent field=expected_archive"):
        orchestration.assert_declared_archive({"test_name": "absent", "expected_archive": False}, output)
```

- [ ] **Step 2: Run RED**

Run: `pytest -m unit tests/unit/test_integration_outcome_contract.py -k 'successful_bam_enforces or declared_archive' -v`

Expected: the successful BAM test fails because the declared missing file is not enforced, and archive tests fail because the helper is missing.

- [ ] **Step 3: Implement minimal GREEN**

Call `assert_declared_artifacts(test_case, output_dir)` after successful coverage
validation in `run_bam_test_case`. Add:

```python
def assert_declared_archive(test_case: dict, output_dir: Path) -> None:
    """Assert the optional default sibling ZIP state for a successful case."""
    if "expected_archive" not in test_case:
        return
    case_id = str(test_case.get("test_name", "<unnamed>"))
    expected = test_case["expected_archive"]
    if not isinstance(expected, bool):
        raise ValueError(f"case={case_id} field=expected_archive must be a Boolean")
    archive = Path(f"{output_dir}.zip")
    if expected:
        assert archive.is_file(), f"case={case_id} field=expected_archive expected file is missing: {archive}"
    else:
        assert not os.path.lexists(archive), f"case={case_id} field=expected_archive unexpected entry is present: {archive}"
```

Add `import os` beside the existing standard-library imports.

- [ ] **Step 4: Wire both local integration callers**

Import and call `assert_declared_archive(case, output_dir)` after
`run_bam_test_case` in both integration modules. Replace the old truthiness-only
archive block in `test_bam_input_with_kestrel_checks`.

- [ ] **Step 5: Extend failure-case hygiene**

Add `"expected_present"` and `"expected_absent"` to the exact `success_only` set in
`test_declared_failures_do_not_retain_unreachable_success_expectations`.

- [ ] **Step 6: Run GREEN and collection**

```bash
pytest -m unit tests/unit/test_integration_outcome_contract.py -v
pytest --collect-only -q -m integration tests/integration/test_pipeline_integration.py tests/integration/test_single_end_pipeline.py
```

Expected: PASS with no marker or import error.

- [ ] **Step 7: Commit**

```bash
git add tests/support/orchestration.py tests/unit/test_integration_outcome_contract.py tests/integration/test_pipeline_integration.py tests/integration/test_single_end_pipeline.py
git commit -m "test(harness): enforce declared artifact states"
```

### Task 3: Prove the keep artifact exists before absence rows

**Files:**
- Modify: `tests/unit/test_integration_outcome_contract.py`
- Modify: `tests/test_data_config.json:867-920`

**Interfaces:**
- Produces: `example_b178_hg19_single_end_keep`.

- [ ] **Step 1: Write RED declaration test**

```python
def test_single_end_keep_case_names_the_real_unmapped_artifact() -> None:
    cases = {case["test_name"]: case for case in load_test_config()["integration_tests"]["single_end_bam_tests"]}
    keep = cases["example_b178_hg19_single_end_keep"]
    assert keep["cli_options"] == ["--keep-intermediates"]
    assert keep["expected_present"] == ["fastq_bam_processing/output_unmapped.bam"]
    assert keep["expected_archive"] is False
```

Import `load_test_config` from `tests.parametrization`.

- [ ] **Step 2: Run RED**

Run: `pytest -m unit tests/unit/test_integration_outcome_contract.py::test_single_end_keep_case_names_the_real_unmapped_artifact -v`

Expected: FAIL with the missing case key.

- [ ] **Step 3: Add minimal declaration**

Rename the existing non-fast `_default` case to
`example_b178_hg19_single_end_keep`; retain its genotype fields and add:

```json
"cli_options": ["--keep-intermediates"],
"expected_present": ["fastq_bam_processing/output_unmapped.bam"],
"expected_archive": false
```

- [ ] **Step 4: Run unit GREEN**

Run the Step 2 command. Expected: PASS.

- [ ] **Step 5: Stage verified worktree-local data and derive CRAM fixtures**

The primary checkout owns the already downloaded git-ignored archive, while the isolated milestone worktree starts
without it. Run these exact non-network commands from the milestone worktree; refuse to merge with or overwrite a
partial destination:

```bash
test ! -e tests/data
VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1 make -C /home/bernt-popp/development/VNtyper verify-test-data
cp -a --reflink=auto /home/bernt-popp/development/VNtyper/tests/data tests/data
VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1 make verify-test-data
VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1 make cram-fixtures
```

Expected: source and destination verification both pass, the copy is worktree-local, and all required derived CRAM
fixtures exist below this worktree's `tests/data/cram/`. Do not symlink the source tree and do not run fixture
derivation in the primary checkout.

- [ ] **Step 6: Run the keep integration node**

```bash
pytest -m integration tests/integration/test_single_end_pipeline.py::test_single_end_bam_produces_a_genotype[example_b178_hg19_single_end_keep] -v
```

Expected: PASS, objectively establishing the literal artifact before Task 5 adds
absence assertions. Stop if the measured path differs.

- [ ] **Step 7: Remove stale case IDs**

Use `rg` to find `_default`; update only references that mean the renamed keep case.

- [ ] **Step 8: Rerun collection**

Run: `rg -n 'single_end_default|single_end_keep' tests && pytest --collect-only -q -m integration tests/integration/test_single_end_pipeline.py`

Expected: no stale `_default` contract and a collected `_keep` node.

- [ ] **Step 9: Commit**

```bash
git add tests/test_data_config.json tests/unit/test_integration_outcome_contract.py
git commit -m "test(integration): pin retained pipeline intermediate"
```

### Task 4: Fix delete-overrides-keep

**Files:**
- Modify: `tests/unit/test_fastq_bam_command_wiring.py:45-210`
- Modify: `vntyper/scripts/fastq_bam_processing.py:81-115,237-249`
- Modify: `docs/about/changelog.md:1-8`

**Interfaces:**
- Produces: `(delete=False, keep=*) -> retain`; `(delete=True, keep=*) -> remove`.

- [ ] **Step 1: Write RED truth-table test**

```python
@pytest.mark.parametrize(
    ("delete_intermediates", "keep_intermediates", "expected_exists"),
    [(False, False, True), (False, True, True), (True, False, False), (True, True, False)],
)
def test_delete_intermediates_overrides_keep(tmp_path: Path, delete_intermediates: bool, keep_intermediates: bool, expected_exists: bool) -> None:
    unmapped = tmp_path / "output_unmapped.bam"
    unmapped.write_bytes(b"intermediate")
    _run_bam_to_fastq(tmp_path, fast_mode=False, delete_intermediates=delete_intermediates, keep_intermediates=keep_intermediates)
    assert unmapped.exists() is expected_exists
```

- [ ] **Step 2: Run RED**

Run: `pytest -m unit tests/unit/test_fastq_bam_command_wiring.py::test_delete_intermediates_overrides_keep -v`

Expected: only `(True, True, False)` FAILS.

- [ ] **Step 3: Implement minimal GREEN**

Replace `if delete_intermediates and not keep_intermediates:` with
`if delete_intermediates:` in `process_bam_to_fastq`; update both flag docstrings to
state that delete wins.

- [ ] **Step 4: Run GREEN**

Run the Step 2 command. Expected: four PASS.

- [ ] **Step 5: Add exact changelog text**

```markdown
## Unreleased

When both `--keep-intermediates` and `--delete-intermediates` are supplied,
`--delete-intermediates` now wins as the CLI help has always documented. Earlier
versions kept the intermediate BAM in this conflicting form.
```

- [ ] **Step 6: Verify touched surfaces**

Run: `pytest -m unit tests/unit/test_fastq_bam_command_wiring.py -v && make docs-build`

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add vntyper/scripts/fastq_bam_processing.py tests/unit/test_fastq_bam_command_wiring.py docs/about/changelog.md
git commit -m "fix(io): make delete intermediates override keep"
```

### Task 5: Add delete, precedence and archive rows

**Files:**
- Modify: `tests/unit/test_integration_outcome_contract.py`
- Modify: `tests/test_data_config.json:867-960`

**Interfaces:**
- Consumes: Task 3's proven literal and Task 4's behavior.
- Produces: `_delete`, `_delete_overrides_keep` and `_archive` cases.

- [ ] **Step 1: Write RED exact-case assertions**

```python
def test_single_end_cases_declare_the_negative_matrix_after_keep() -> None:
    cases = {case["test_name"]: case for case in load_test_config()["integration_tests"]["single_end_bam_tests"]}
    artifact = "fastq_bam_processing/output_unmapped.bam"
    assert cases["example_b178_hg19_single_end_delete"]["cli_options"] == ["--delete-intermediates"]
    assert cases["example_b178_hg19_single_end_delete"]["expected_absent"] == [artifact]
    assert cases["example_b178_hg19_single_end_delete_overrides_keep"]["cli_options"] == ["--keep-intermediates", "--delete-intermediates"]
    assert cases["example_b178_hg19_single_end_delete_overrides_keep"]["expected_absent"] == [artifact]
    assert cases["example_b178_hg19_single_end_archive"]["cli_options"] == ["--archive-results"]
    assert cases["example_b178_hg19_single_end_archive"]["expected_archive"] is True
```

- [ ] **Step 2: Run RED**

Run: `pytest -m unit tests/unit/test_integration_outcome_contract.py::test_single_end_cases_declare_the_negative_matrix_after_keep -v`

Expected: FAIL on the first missing case.

- [ ] **Step 3: Add minimal JSON rows**

Insert these complete literal objects; do not copy fields at execution time and do not carry `expected_present` from
the keep row:

```json
{
  "test_name": "example_b178_hg19_single_end_delete",
  "fixture_name": "example_b178_hg19_single_end",
  "reference_assembly": "hg19",
  "cli_options": ["--delete-intermediates"],
  "kestrel_assertions": {
    "Estimated_Depth_AlternateVariant": {"value": 416, "tolerance_percentage": 5},
    "Estimated_Depth_Variant_ActiveRegion": {"value": 6111, "tolerance_percentage": 5},
    "Depth_Score": {"value": 0.0677466863033873, "tolerance_percentage": 5},
    "Confidence": "High_Precision*"
  },
  "check_igv_report": false,
  "expected_absent": ["fastq_bam_processing/output_unmapped.bam"],
  "expected_archive": false
}
```

```json
{
  "test_name": "example_b178_hg19_single_end_delete_overrides_keep",
  "fixture_name": "example_b178_hg19_single_end",
  "reference_assembly": "hg19",
  "cli_options": ["--keep-intermediates", "--delete-intermediates"],
  "kestrel_assertions": {
    "Estimated_Depth_AlternateVariant": {"value": 416, "tolerance_percentage": 5},
    "Estimated_Depth_Variant_ActiveRegion": {"value": 6111, "tolerance_percentage": 5},
    "Depth_Score": {"value": 0.0677466863033873, "tolerance_percentage": 5},
    "Confidence": "High_Precision*"
  },
  "check_igv_report": false,
  "expected_absent": ["fastq_bam_processing/output_unmapped.bam"],
  "expected_archive": false
}
```

```json
{
  "test_name": "example_b178_hg19_single_end_archive",
  "fixture_name": "example_b178_hg19_single_end",
  "reference_assembly": "hg19",
  "cli_options": ["--archive-results"],
  "kestrel_assertions": {
    "Estimated_Depth_AlternateVariant": {"value": 416, "tolerance_percentage": 5},
    "Estimated_Depth_Variant_ActiveRegion": {"value": 6111, "tolerance_percentage": 5},
    "Depth_Score": {"value": 0.0677466863033873, "tolerance_percentage": 5},
    "Confidence": "High_Precision*"
  },
  "check_igv_report": false,
  "expected_archive": true
}
```

- [ ] **Step 4: Run unit GREEN**

Run the Step 2 command. Expected: PASS.

- [ ] **Step 5: Run the three real matrix nodes**

Run: `pytest -m integration tests/integration/test_single_end_pipeline.py -k 'single_end_delete or single_end_archive' -v`

Expected: three PASS.

- [ ] **Step 6: Validate JSON and collection**

Run: `python -m json.tool tests/test_data_config.json >/dev/null && pytest --collect-only -q -m integration tests/integration/test_single_end_pipeline.py`

Expected: valid JSON and explicit case IDs.

- [ ] **Step 7: Commit**

```bash
git add tests/test_data_config.json tests/unit/test_integration_outcome_contract.py
git commit -m "test(integration): cover cleanup and archive matrix"
```

### Task 6: Add paired b178 FASTQ without SHARK

**Files:**
- Modify: `tests/unit/test_integration_outcome_contract.py`
- Modify: `tests/test_data_config.json:967-1010`

**Interfaces:**
- Produces: `example_b178_hg19_subset_paired_fastq_no_shark`.

- [ ] **Step 1: Write RED exact declaration test**

```python
def test_alternate_paired_fastq_uses_b178_and_omits_shark() -> None:
    case = {case["test_name"]: case for case in get_fastq_test_cases()}["example_b178_hg19_subset_paired_fastq_no_shark"]
    assert case["fastq1"] == "tests/data/example_b178_hg19_subset_R1.fastq.gz"
    assert case["fastq2"] == "tests/data/example_b178_hg19_subset_R2.fastq.gz"
    assert case["reference_assembly"] == "hg19"
    assert case["expected_files"] == ["summary_report.html", "kestrel/kestrel_result.tsv"]
    runner = mock.Mock(return_value=0)
    with mock.patch.object(orchestration, "assert_required_files"):
        orchestration.run_fastq_test_case(case, runner, Path("output"))
    assert runner.call_args.args[4] == []
```

- [ ] **Step 2: Run RED**

Run: `pytest -m unit tests/unit/test_integration_outcome_contract.py::test_alternate_paired_fastq_uses_b178_and_omits_shark -v`

Expected: FAIL with the missing case key.

- [ ] **Step 3: Add minimal JSON declaration**

```json
{
  "test_name": "example_b178_hg19_subset_paired_fastq_no_shark",
  "fastq1": "tests/data/example_b178_hg19_subset_R1.fastq.gz",
  "fastq2": "tests/data/example_b178_hg19_subset_R2.fastq.gz",
  "reference_assembly": "hg19",
  "expected_files": ["summary_report.html", "kestrel/kestrel_result.tsv"]
}
```

Do not add `cli_options`.

- [ ] **Step 4: Run unit GREEN**

Run the Step 2 command. Expected: PASS.

- [ ] **Step 5: Run the real paired FASTQ node**

Run: `pytest -m integration tests/integration/test_pipeline_integration.py::test_fastq_input[example_b178_hg19_subset_paired_fastq_no_shark] -v`

Expected: exit 0 and both declared files exist.

- [ ] **Step 6: Run final verification gates**

```bash
pytest -m unit tests/unit/test_integration_outcome_contract.py tests/unit/test_fastq_bam_command_wiring.py -v
make test-unit-cov
make patch-coverage
make test-integration
make check-all
git diff --check
```

Expected: all PASS, branch-inclusive total stays at least 86%, patch coverage is at
least 80%, and no #226/Docker/API/subcommand implementation changed.

- [ ] **Step 7: Commit and close only #71**

```bash
git add tests/test_data_config.json tests/unit/test_integration_outcome_contract.py
git commit -m "test(integration): add paired b178 fastq path" -m "Closes #71"
```

Close #226 separately with existing evidence; never add `Closes #226` to these commits.
