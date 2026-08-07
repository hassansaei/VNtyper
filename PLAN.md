# Milestone 2 — "Correctness of reported numbers" Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development
> to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Read [SPEC.md](SPEC.md) first — it is the contract; this file is the order of operations.

**Goal:** Close #171, #172, #174, #203 and #212 so that every number VNtyper reports is
computed over the base set it claims, gates on the metrics it displays, excludes known
artifacts from calls, stops writing a coordinate nothing can resolve, and can never turn
an interrupted run into a silent negative.

**Architecture:** Four independent lanes in isolated git worktrees. Pure decision logic
goes into focused modules (`coverage_qc.py` new, `flagging.py` extended) and the I/O
stays where it is, per AGENTS.md rule 3. Every behaviour change is gated by a test that
fails before the change and passes after. Lanes are rebased onto the milestone branch in
a fixed order so history stays linear and one PR closes all five issues.

**Tech Stack:** Python 3.10 floor (CI matrix 3.10–3.13), pandas 2.2.3, pytest 9.1.1,
ruff (line length 120, double quotes), mypy, Jinja2, samtools, the vendored Kestrel 1.0.1 JAR.

## Global Constraints

Copied verbatim from AGENTS.md and SPEC.md. **Every task's requirements implicitly include
this section.**

- **Environment.** `conda activate vntyper`, then `export PATH="$CONDA_PREFIX/bin:$PATH"`.
  Without the PATH line, `~/.local/bin/pytest` (system Python 3.14) shadows the env and
  every test module fails collection with `ModuleNotFoundError: No module named 'pandas'`.
- **Run pytest from the repo root.** `tests/parametrization.py` opens
  `tests/test_data_config.json` by relative path at collection time.
- **Code must run on Python 3.10.** Local interpreter is 3.12; 3.11+ syntax passes locally
  and breaks CI. PEP 604 (`X | None`) and PEP 585 (`list[str]`) are available.
- **Ruff only. Line length 120, double quotes, `target-version = "py310"`.** Do not
  reformat to 88.
- **Every new test file declares `pytestmark = pytest.mark.unit`** after the imports. CI
  runs `pytest -m unit`; an unmarked file silently never runs, and
  `tests/unit/test_marker_hygiene.py` fails the build naming it.
- **Unit tests stay pure**: `tmp_path` + `unittest.mock`, no network, no Docker, no
  reference data.
- **Logging**: `logger = logging.getLogger(__name__)` after imports; f-strings in log calls
  are the established style. Never `logging.info(...)` on the root logger.
- **Errors**: no custom exception classes. `logger.error(msg)` then
  `raise ValueError(msg)` / `RuntimeError(msg)` with the same message.
- **Google-style docstrings** (`Args:` / `Returns:` / `Raises:`) on public functions.
- **Research use only. Clinical-sounding output text is config-driven — never invent or
  reword it.** No task in this plan adds or edits a message string in
  `report_config.json`.
- **Config-driven, never hardcoded.** Artifact flag names and thresholds live in
  `kestrel_config.json` / `config.json`. Python reads the list; it never names a flag
  inline.
- **Config is loaded at import time** (AGENTS.md trap 1). `kestrel_genotyping.py` reads
  `kestrel_config.json` into a module global; tests must patch the module global, not pass
  a `config` argument.
- **Conventional Commits**: `type(scope): lowercase subject` under ~72 chars, `Closes #N`
  in the footer. PRs are **not** squashed — every commit is permanent history.
- **Never push a `v*.*.*` tag.** It publishes to PyPI immediately and irreversibly.
- **Never claim tests pass without showing the command output.**

## Lanes and dependency graph

```
Lane A (sequential, the only chain):  T1 → T2 → T3 → T4 → T5 → T6
Lane B (independent):                 T7 → T8 → T9
Lane C (independent):                 T10 → T11
Lane D (independent):                 T12 → T13
                                       ↓
Integration (after all lanes):        T14 → T15 → T16
```

| Lane | Worktree | Branch | Issues |
| --- | --- | --- | --- |
| A | `.claude/worktrees/m2-lane-a` | `fix/issue-171-172-coverage-arithmetic-and-qc` | #171 then #172 |
| B | `.claude/worktrees/m2-lane-b` | `fix/issue-174-artifact-gate` | #174 |
| C | `.claude/worktrees/m2-lane-c` | `fix/issue-203-pos-fasta` | #203 + BED rider |
| D | `.claude/worktrees/m2-lane-d` | `fix/issue-212-kestrel-skip` | #212 |

**A, B, C and D run concurrently.** They touch `kestrel_genotyping.py` in three disjoint
regions (B: `filter_cols` ~line 826 and step 6.5 ~line 581; C: `generate_bed_file` ~line
607; D: `run_kestrel` ~line 294), so the rebase in T14 is expected to be conflict-free.

**All four worktrees already exist** at commit `989b160` and must be rebased onto the
current milestone branch tip before work starts:

```bash
cd /home/bernt-popp/development/VNtyper/.claude/worktrees/m2-lane-<X>
git rebase fix/milestone-2-correctness-of-reported-numbers
```

## File structure

| File | Status | Responsibility | Lane |
| --- | --- | --- | --- |
| `vntyper/scripts/coverage_stats.py` | modify | region-based arithmetic; the frozen schema | A |
| `vntyper/scripts/command_builders.py` | modify | the depth command gains `-a` | A |
| `vntyper/scripts/coverage_qc.py` | **create** | the pure QC verdict — one decision, no I/O | A |
| `vntyper/scripts/fastq_bam_processing.py` | modify | calls the verdict, merges it into the stats it writes | A |
| `vntyper/scripts/report_formatting.py` | modify | `COVERAGE_FIELD_TYPES` gains the ninth column | A |
| `vntyper/scripts/screening_summary.py` | modify | quality axis takes a `CoverageQC` | A |
| `vntyper/scripts/generate_report.py` | modify | evaluates the verdict once, feeds axis + context | A |
| `vntyper/templates/report_template.html` | modify | one new Coverage QC row | A |
| `vntyper/scripts/cohort_summary.py` | modify | third export: `cohort_stats` | A |
| `vntyper/scripts/flagging.py` | modify | `add_artifact_gate`, pure | B |
| `vntyper/scripts/kestrel_config.json` | modify | `artifact_flags` declaration | B |
| `vntyper/scripts/kestrel_genotyping.py` | modify | 6th gate (B), BED coords (C), no skip (D) | B, C, D |
| `vntyper/scripts/motif_processing.py` | modify | delete the dead rebase | C |
| `vntyper/scripts/summary.py` | modify | loud missing result file | D |
| `tests/builders.py` | modify | `STAGE_COLUMNS` gains the 6th gate | B |
| `tests/helpers.py` | modify | its own `COVERAGE_COLUMNS` copy gains the 9th | A |

---

## Lane A — #171 then #172

### Task 1: Coverage arithmetic over the region (#171)

**Lane A. Depends on: nothing. Blocks: T2, T3.**

**Files:**
- Modify: `vntyper/scripts/coverage_stats.py:20-26` (module docstring), `:130-179`
  (`summarise_coverage`)
- Test: `tests/unit/test_coverage_stats.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `summarise_coverage(coverage_values: list[int], total_region_length: int) -> dict`
  — signature **unchanged**; the returned values change. Keys stay exactly
  `COVERAGE_COLUMNS`.

- [ ] **Step 1: Write the failing tests**

Add to `tests/unit/test_coverage_stats.py`:

```python
def test_the_mean_is_over_the_region_not_over_covered_positions():
    """#171: ten positions covered at 20x inside a 100 bp region is a mean of 2.0.

    The old definition divided by the covered-base count and reported 20.0, which is
    the depth *where there were reads*, not the depth of the region.
    """
    stats = summarise_coverage([20] * 10, total_region_length=100)

    assert stats["mean"] == 2.0
    assert stats["uncovered_bases"] == 90
    assert stats["percent_uncovered"] == 90.0


def test_uncovered_bases_counts_zero_rows_when_samtools_emitted_every_position():
    """#171: with `-a` the file already carries the zeros, so subtraction reports 0.

    This is the regression that adding `-a` alone would have introduced: under the old
    formula `region_length - len(values)` is `4 - 4 = 0` here, and the one metric that
    was correct becomes the one that is always wrong.
    """
    stats = summarise_coverage([10, 0, 0, 10], total_region_length=4)

    assert stats["uncovered_bases"] == 2
    assert stats["percent_uncovered"] == 50.0
    assert stats["mean"] == 5.0
    assert stats["min"] == 0, "the true minimum of a partly uncovered region is zero"


def test_a_legacy_depth_file_without_a_is_padded_to_the_region():
    """A file with no zero rows still yields the region-wide statistics."""
    stats = summarise_coverage([10, 20, 30], total_region_length=6)

    assert stats["uncovered_bases"] == 3
    assert stats["mean"] == 10.0
    assert stats["min"] == 0
    assert stats["max"] == 30


@pytest.mark.parametrize(
    ("values", "length"),
    [([20] * 10, 100), ([10, 20, 30], 1501), ([5] * 999, 1000), ([7], 2)],
)
def test_the_closed_form_identity_reconciles_old_and_new_means(values, length):
    """#171's free regression test: `mean_new == mean_old * (1 - pct_old / 100)`.

    With `S = sum(v)`, `c = len(v)`, `T = length`: `(S/c) * (c/T) == S/T`. It lets a
    historical report be reconciled with a re-run without new fixtures. Guarded on
    `T > 0`, because an unparseable region forces `percent_uncovered = 0`.
    """
    assert length > 0
    old_mean = sum(values) / len(values)
    old_percent = (length - len(values)) / length * 100

    stats = summarise_coverage(values, total_region_length=length)

    assert stats["mean"] == pytest.approx(old_mean * (1 - old_percent / 100))


def test_a_region_shorter_than_the_depth_file_does_not_invent_negative_padding():
    """Degenerate but reachable if a region string and a BAM disagree."""
    stats = summarise_coverage([10, 20, 30, 40], total_region_length=2)

    assert stats["mean"] == 25.0, "no padding is added; the four observed values stand"
    assert stats["uncovered_bases"] == 0
```

Replace the body of the existing `test_partial_coverage_counts_the_missing_positions_as_uncovered`
(currently `test_coverage_stats.py:136`) — the new `test_the_mean_is_over_the_region...`
supersedes it — and change `test_a_zero_length_region_reports_zero_percent_rather_than_dividing_by_zero`
(`:173`) to:

```python
    assert stats["uncovered_bases"] == 0, (
        "#171: a negative uncovered count was nonsense. `absent` is floored at 0, so an "
        "unparseable region reports 0 uncovered rather than -2. The WARNING from "
        "parse_region_length is still the only signal that the region was unreadable."
    )
```

- [ ] **Step 2: Run the tests and verify they fail**

```bash
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vntyper \
  && export PATH="$CONDA_PREFIX/bin:$PATH" \
  && pytest tests/unit/test_coverage_stats.py -m unit -o log_cli=false -q
```
Expected: FAIL. `test_the_mean_is_over_the_region_not_over_covered_positions` reports
`assert 20.0 == 2.0`.

- [ ] **Step 3: Implement**

Replace `summarise_coverage`'s body (`coverage_stats.py:160-179`):

```python
    if not coverage_values:
        raise RuntimeError("No coverage data found.")

    # `samtools depth -a` emits one row per position in the region, zeros included, so
    # `coverage_values` *is* the region. A file written without `-a` - a legacy artefact,
    # or a region truncated at a contig end - is short by exactly the positions that had
    # no reads, so restoring them as zeros makes both cases the same base set. `absent` is
    # floored at 0: an unparseable region (length 0) must not subtract padding. See #171.
    absent = max(total_region_length - len(coverage_values), 0)
    base = coverage_values + [0] * absent

    uncovered_bases = sum(1 for depth in base if depth == 0)
    percent_uncovered = 0 if total_region_length <= 0 else uncovered_bases / total_region_length * 100

    return {
        "mean": sum(base) / len(base),
        "median": statistics.median(base),
        "stdev": statistics.stdev(base) if len(base) > 1 else 0,
        "min": min(base),
        "max": max(base),
        "region_length": total_region_length,
        "uncovered_bases": uncovered_bases,
        "percent_uncovered": percent_uncovered,
    }
```

Rewrite the module docstring at `:20-26`, which currently states the covered-position
definition as intentional:

```
One definition is worth stating explicitly because the opposite was true until #171:
every statistic here is over the **region**, not over the covered positions. ``mean``,
``median``, ``stdev``, ``min`` and ``max`` are computed after uncovered positions are
restored as zeros, so a sample covered at 30x across 10% of the VNTR reports
``mean = 3``, ``min = 0`` and ``percent_uncovered = 90``. The old ``mean = 30`` was the
depth where there happened to be reads, which is not a property of the region and was
systematically too high exactly where coverage was patchy.
```

Update `summarise_coverage`'s `Note:` (`:156-158`) to match, and its `Raises:` to say that
under `-a` an empty table now means the region matched no contig at all.

- [ ] **Step 4: Run the tests and verify they pass**

```bash
pytest tests/unit/test_coverage_stats.py -m unit -o log_cli=false -q
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/coverage_stats.py tests/unit/test_coverage_stats.py
git commit -m "fix(coverage): compute every statistic over the region, not covered bases

Closes #171"
```

---

### Task 2: `-a` on the depth command, and the wiring fixtures (#171)

**Lane A. Depends on: T1. Blocks: T3.**

**Files:**
- Modify: `vntyper/scripts/command_builders.py:383-386`
- Test: `tests/unit/test_command_builders.py:241`,
  `tests/unit/test_fastq_bam_command_wiring.py:379`

**Interfaces:**
- Consumes: T1's `summarise_coverage`.
- Produces: `build_samtools_depth_command(...) -> str` — signature unchanged, output gains
  `-a` immediately after `depth`.

- [ ] **Step 1: Write the failing test**

Replace `test_the_depth_command_is_pinned` in `tests/unit/test_command_builders.py:241`:

```python
def test_the_depth_command_requests_every_position_in_the_region():
    """#171: without `-a`, samtools emits only covered positions.

    Every statistic downstream is then over the covered bases rather than the region,
    and `uncovered_bases` - which is derived by subtraction - silently reads 0. The flag
    is asserted as a parsed argument, not as a substring, so `-all` could not satisfy it.
    """
    command = build_samtools_depth_command(
        samtools_path=SAMTOOLS,
        threads=4,
        region="chr1:155160500-155162000",
        bam_file="/data/sample.bam",
        coverage_output="/out/cov_vntr_coverage.txt",
    )

    assert command == (
        "samtools depth -a -@ 4 -r chr1:155160500-155162000 /data/sample.bam > /out/cov_vntr_coverage.txt"
    )
    assert "-a" in shlex.split(command.split(">")[0])
```

- [ ] **Step 2: Run it and verify it fails**

```bash
pytest tests/unit/test_command_builders.py -m unit -o log_cli=false -q
```
Expected: FAIL — the produced string has no `-a`.

- [ ] **Step 3: Implement**

`command_builders.py:383-386`:

```python
    return (
        f"{samtools_path} depth -a -@ {quote_path(threads)} -r {quote_path(region)} "
        f"{quote_path(bam_file)} > {quote_path(coverage_output)}"
    )
```

Add to the docstring's `Note:`:

```
``-a`` makes samtools emit a row for every position in the region, zero-depth ones
included. Without it the depth table is the covered positions only, and every statistic
computed from it is over a base set that varies per sample (#171).
```

- [ ] **Step 4: Run it and verify it passes; then fix the end-to-end fixture**

`tests/unit/test_fastq_bam_command_wiring.py:379`
(`test_calculate_vntr_coverage_writes_the_frozen_tsv_schema`) asserts exact TSV bytes.
Recompute them under T1's arithmetic for whatever depth rows the test's
`fake_run_command` writes, and add to the docstring: *"These bytes changed with #171; they
are the region-wide statistics, not the covered-position ones."* Do not weaken the
assertion to a substring check.

```bash
pytest tests/unit/test_command_builders.py tests/unit/test_fastq_bam_command_wiring.py -m unit -o log_cli=false -q
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/command_builders.py tests/unit/test_command_builders.py tests/unit/test_fastq_bam_command_wiring.py
git commit -m "fix(coverage): pass -a so samtools depth emits uncovered positions

Closes #171"
```

---

### Task 3: The `coverage_qc` verdict, as a pure module (#172)

**Lane A. Depends on: T1, T2. Blocks: T4, T5.**

**Files:**
- Create: `vntyper/scripts/coverage_qc.py`
- Test: `tests/unit/test_coverage_qc.py` (create)

**Interfaces:**
- Consumes: nothing at runtime.
- Produces:
  - `COVERAGE_QC_PASS: str = "PASS"`, `COVERAGE_QC_FAIL: str = "FAIL"`
  - `REASON_MEAN: str = "mean_vntr_coverage"`, `REASON_UNCOVERED: str = "percent_vntr_uncovered"`
  - `@dataclass(frozen=True) class CoverageQC` with `passed: bool`, `status: str`,
    `reasons: tuple[str, ...]`
  - `evaluate_coverage_qc(mean_vntr_coverage: float | None, percent_vntr_uncovered: float | None, mean_threshold: float, percent_threshold: float) -> CoverageQC`

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_coverage_qc.py`:

```python
"""
Unit tests for the coverage QC verdict (#172).

The verdict is the one thing that decides `quality_metrics_pass`, and it is also written
into the coverage TSV, so a disagreement between the two would put a FAIL column beside a
PASS sentence. These tests pin the rule, its boundaries and its rounding contract.
"""

import pytest

from vntyper.scripts.coverage_qc import (
    COVERAGE_QC_FAIL,
    COVERAGE_QC_PASS,
    REASON_MEAN,
    REASON_UNCOVERED,
    evaluate_coverage_qc,
)

pytestmark = pytest.mark.unit


def test_a_sample_above_both_thresholds_passes():
    qc = evaluate_coverage_qc(250.0, 5.0, 100, 50.0)

    assert qc.passed is True
    assert qc.status == COVERAGE_QC_PASS
    assert qc.reasons == ()


def test_a_low_mean_fails_and_names_the_mean():
    qc = evaluate_coverage_qc(99.0, 5.0, 100, 50.0)

    assert qc.passed is False
    assert qc.status == COVERAGE_QC_FAIL
    assert qc.reasons == (REASON_MEAN,)


def test_a_patchy_vntr_fails_even_with_acceptable_mean():
    """#172's headline case, and the reason the issue exists.

    Half the VNTR uncovered is precisely where a frameshift call can be missed, yet
    before this change the sample passed QC on its mean alone.
    """
    qc = evaluate_coverage_qc(250.0, 80.0, 100, 50.0)

    assert qc.passed is False
    assert qc.reasons == (REASON_UNCOVERED,)


def test_both_failures_are_reported_in_declaration_order():
    qc = evaluate_coverage_qc(10.0, 90.0, 100, 50.0)

    assert qc.reasons == (REASON_MEAN, REASON_UNCOVERED)


@pytest.mark.parametrize("mean", [None, 250.0])
@pytest.mark.parametrize("percent", [None, 5.0])
def test_a_metric_that_was_never_measured_does_not_fail_the_gate(mean, percent):
    """Preserved from the pre-#172 behaviour, deliberately.

    `test_screening_summary.py` pins that an unmeasured sample reports as passing. That
    is a displayed interpretation for every run with no Coverage Calculation step, so
    #172 does not change it - it only adds the second metric.
    """
    assert evaluate_coverage_qc(mean, percent, 100, 50.0).passed is True


def test_the_boundaries_pass_on_equality():
    """Mean fails strictly below; uncovered fails strictly above. Asymmetric on purpose:
    it matches `threshold_icon`'s `higher_better` argument in report_formatting."""
    assert evaluate_coverage_qc(100.0, 50.0, 100, 50.0).passed is True
    assert evaluate_coverage_qc(99.99, 50.0, 100, 50.0).passed is False
    assert evaluate_coverage_qc(100.0, 50.01, 100, 50.0).passed is False


def test_the_verdict_is_a_function_of_the_published_figures():
    """The rounding contract (#172, adversarial review A1).

    `format_coverage_summary` writes mean and percent with `:.2f`, and the report reads
    those strings back. If the writer evaluated the raw value and the reader the rounded
    one, a raw mean of 99.999 would emit FAIL into the TSV while the report recomputed
    PASS. Callers round before calling; this test pins that both sides then agree, and
    that the answer matches what the report prints.
    """
    raw, published = 99.999, round(99.999, 2)

    assert published == 100.0
    assert evaluate_coverage_qc(raw, 0.0, 100, 50.0).passed is False
    assert evaluate_coverage_qc(published, 0.0, 100, 50.0).passed is True
```

- [ ] **Step 2: Run and verify it fails**

```bash
pytest tests/unit/test_coverage_qc.py -m unit -o log_cli=false -q
```
Expected: FAIL — `ModuleNotFoundError: No module named 'vntyper.scripts.coverage_qc'`.

- [ ] **Step 3: Implement**

Create `vntyper/scripts/coverage_qc.py`:

```python
"""
coverage_qc.py

The sample-level coverage quality verdict (#172).

Before this module, ``quality_metrics_pass`` considered mean VNTR coverage only.
``percent_vntr_uncovered`` was configured with a threshold, computed on every run and
compared to nothing - it drove a red or green icon and no decision. A sample with
acceptable mean coverage and half the VNTR uncovered passed QC, which is the opposite of
the desirable failure mode: a patchy VNTR is exactly where a frameshift call can be
missed.

**The verdict is a function of the *published* figures, not the raw ones.**
``coverage_stats.format_coverage_summary`` writes ``mean`` and ``percent_uncovered`` with
two decimal places, and the report reads those strings back out of
``pipeline_summary.json``. Evaluating the raw value in one place and the rounded value in
the other lets the two disagree at a threshold boundary - a raw mean of 99.999 is below a
threshold of 100, but serialises as ``100.00``. Callers therefore round before calling,
and the report prints no ``FAIL`` beside a displayed ``100.00``.

Functions:
    evaluate_coverage_qc: Two metrics and two thresholds to a verdict
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

logger = logging.getLogger(__name__)

#: The verdict written into the coverage summary's ``coverage_qc`` column.
COVERAGE_QC_PASS = "PASS"
COVERAGE_QC_FAIL = "FAIL"

#: Reason identifiers, named after the ``config.json`` threshold keys they come from so a
#: consumer can look up the number that was applied.
REASON_MEAN = "mean_vntr_coverage"
REASON_UNCOVERED = "percent_vntr_uncovered"


@dataclass(frozen=True)
class CoverageQC:
    """The coverage quality verdict for one sample.

    Attributes:
        passed: Whether the sample met every configured coverage threshold.
        status: :data:`COVERAGE_QC_PASS` or :data:`COVERAGE_QC_FAIL`. This is the value
            written to the ``coverage_qc`` column.
        reasons: The threshold keys that failed, in declaration order. Empty when passed.
    """

    passed: bool
    status: str
    reasons: tuple[str, ...]


def evaluate_coverage_qc(
    mean_vntr_coverage: float | None,
    percent_vntr_uncovered: float | None,
    mean_threshold: float,
    percent_threshold: float,
) -> CoverageQC:
    """Decide whether a sample's VNTR coverage meets both configured thresholds.

    Args:
        mean_vntr_coverage (float | None): Mean depth over the VNTR region, as published
            (two decimal places). ``None`` when no coverage step ran.
        percent_vntr_uncovered (float | None): Percentage of the region at zero depth, as
            published. ``None`` when no coverage step ran.
        mean_threshold (float): ``config.json``'s ``thresholds.mean_vntr_coverage``.
        percent_threshold (float): ``config.json``'s ``thresholds.percent_vntr_uncovered``.

    Returns:
        CoverageQC: The verdict. A metric that is ``None`` never fails the gate - an
        unmeasured sample reported as failing would change a displayed interpretation for
        every run with no coverage step, which is out of scope for #172.

    Note:
        The comparisons are asymmetric on purpose, matching ``threshold_icon``'s
        ``higher_better`` argument: the mean fails strictly *below* its threshold, the
        uncovered fraction fails strictly *above* its own. A sample at exactly 100x and
        exactly 50.0% uncovered passes both.
    """
    reasons: list[str] = []

    if mean_vntr_coverage is not None and mean_vntr_coverage < mean_threshold:
        reasons.append(REASON_MEAN)
    if percent_vntr_uncovered is not None and percent_vntr_uncovered > percent_threshold:
        reasons.append(REASON_UNCOVERED)

    passed = not reasons
    if not passed:
        logger.info(f"Coverage QC failed on: {', '.join(reasons)}")

    return CoverageQC(
        passed=passed,
        status=COVERAGE_QC_PASS if passed else COVERAGE_QC_FAIL,
        reasons=tuple(reasons),
    )
```

- [ ] **Step 4: Run and verify it passes**

```bash
pytest tests/unit/test_coverage_qc.py -m unit -o log_cli=false -q
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/coverage_qc.py tests/unit/test_coverage_qc.py
git commit -m "feat(coverage): add the pure coverage QC verdict

Refs #172"
```

---

### Task 4: Emit `coverage_qc` and widen the quality axis (#172)

**Lane A. Depends on: T3. Blocks: T5.**

**Files:**
- Modify: `vntyper/scripts/coverage_stats.py:48-57` (`COVERAGE_COLUMNS`)
- Modify: `vntyper/scripts/fastq_bam_processing.py:325-348`
- Modify: `vntyper/scripts/report_formatting.py:117-126` (`COVERAGE_FIELD_TYPES`)
- Modify: `vntyper/scripts/screening_summary.py:23-32`, `:72-92`, `:241-319`
- Modify: `vntyper/scripts/generate_report.py:410-414`, `:449-452`, `:494-503`, context dict
- Modify: `vntyper/templates/report_template.html` (after the Percentage Uncovered row)
- Modify: `tests/helpers.py:405`
- Test: `tests/unit/test_screening_summary.py`, `tests/unit/test_report_formatting.py`,
  `tests/unit/test_coverage_stats.py`, `tests/unit/test_helpers.py`

**Interfaces:**
- Consumes: `evaluate_coverage_qc`, `CoverageQC` from T3.
- Produces: `build_screening_summary(kestrel_df, advntr_df, advntr_available, coverage_qc: CoverageQC, report_config) -> ScreeningSummary`
  — **five** positional arguments, down from six. `COVERAGE_COLUMNS` gains a ninth entry
  `"coverage_qc"` at the end.

- [ ] **Step 1: Write the failing tests**

```python
# tests/unit/test_coverage_stats.py — extend the frozen-schema literal
    assert COVERAGE_COLUMNS == (
        "mean", "median", "stdev", "min", "max",
        "region_length", "uncovered_bases", "percent_uncovered",
        "coverage_qc",
    )
```

```python
# tests/unit/test_screening_summary.py — new
def test_a_patchy_sample_fails_the_screening_quality_axis(report_config) -> None:
    """#172: acceptable mean, half the VNTR uncovered. Before this change it passed."""
    qc = evaluate_coverage_qc(250.0, 80.0, 100, 50.0)
    summary = ss.build_screening_summary(pd.DataFrame(), pd.DataFrame(), False, qc, report_config)

    assert summary.quality_metrics_pass is False


def test_widening_the_quality_axis_cannot_change_positivity(report_config) -> None:
    """`is_positive` derives from the two algorithm axes only, never from quality.

    Stated as a test because it is the one way #172 could have become
    genotype-affecting, and it is not.
    """
    failing = evaluate_coverage_qc(1.0, 99.0, 100, 50.0)
    passing = evaluate_coverage_qc(250.0, 1.0, 100, 50.0)

    assert (
        ss.build_screening_summary(pd.DataFrame(), pd.DataFrame(), False, failing, report_config).is_positive
        is ss.build_screening_summary(pd.DataFrame(), pd.DataFrame(), False, passing, report_config).is_positive
    )
```

Convert **all ten** existing `build_screening_summary` call sites in
`tests/unit/test_screening_summary.py` (lines 329, 351, 357, 364, 371, 378, 384, 392, 404,
417) from `(..., mean, threshold, report_config)` to `(..., CoverageQC, report_config)`.
Preserve each test's meaning: `:384`'s below-threshold mean becomes
`evaluate_coverage_qc(99.0, 0.0, 100, 50.0)`; `:392`'s unmeasured case becomes
`evaluate_coverage_qc(None, None, 100, 50.0)`; `:404` and `:417` keep driving the
exception path with `Exploding()`.

- [ ] **Step 2: Run and verify they fail**

```bash
pytest tests/unit/test_coverage_stats.py tests/unit/test_screening_summary.py -m unit -o log_cli=false -q
```
Expected: FAIL — the schema literal mismatches and `build_screening_summary` takes six args.

- [ ] **Step 3: Implement**

1. `coverage_stats.py`: append `"coverage_qc"` to `COVERAGE_COLUMNS`. It is **not** in
   `_TWO_DECIMAL_COLUMNS`, so it is written verbatim. Note in the docstring that the
   column is a verdict rather than a measurement and is filled by the caller.
2. `report_formatting.py:117-126`: add `"coverage_qc": str`.
3. `fastq_bam_processing.py`, inside `calculate_vntr_coverage` after `summarise_coverage`:

```python
        thresholds = config.get("thresholds", {})
        # `.get` with the shipped defaults rather than `[...]`: `--config-path` replaces
        # the whole config instead of merging (AGENTS.md trap 2), so a caller-supplied
        # config legitimately lacks these keys and a KeyError here would abort a run over
        # a display threshold.
        qc = evaluate_coverage_qc(
            round(stats["mean"], 2),
            round(stats["percent_uncovered"], 2),
            thresholds.get("mean_vntr_coverage", 100),
            thresholds.get("percent_vntr_uncovered", 50.0),
        )
        stats["coverage_qc"] = qc.status
```

   The `round` calls are the A1 contract: the verdict is evaluated on the same figures
   `format_coverage_summary` is about to write.
4. `screening_summary.py`: change `build_screening_summary`'s fourth and fifth parameters
   to a single `coverage_qc: CoverageQC`, set
   `quality_metrics_pass = coverage_qc.passed`, and update the module docstring's third
   axis (`:29`) and `ScreeningSummary.quality_metrics_pass`'s docstring (`:82`) — both
   still say "mean VNTR coverage against its threshold".
5. `generate_report.py`: build the verdict once from the parsed coverage and the
   thresholds it already loads, pass it to `build_screening_summary`, and add
   `"coverage_qc": qc.status` to the context. It recomputes rather than reading the
   stored column so `vntyper report` still works on a summary written before this change;
   both sides evaluate published figures, so they cannot disagree.
6. `report_template.html`: add after the Percentage Uncovered row:

```html
                    <tr>
                        <td>Coverage QC</td>
                        <td>{{ coverage_qc }}</td>
                    </tr>
```

7. `tests/helpers.py:405`: append `"coverage_qc"` to its own `COVERAGE_COLUMNS` copy.
   `test_helpers.py:83` pins it to production's tuple. `validate_coverage_output` floats
   only `mean`, `median` and `percent_uncovered`, so the string column needs no other change.

- [ ] **Step 4: Run the full unit tier**

```bash
pytest -m unit tests/unit -o log_cli=false -q
```
Expected: PASS, 2592 + the new tests.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(report): gate quality on both coverage metrics and emit coverage_qc

Closes #172"
```

---

### Task 5: The cohort statistics export (#172)

**Lane A. Depends on: T4. Blocks: T6.**

**Files:**
- Modify: `vntyper/scripts/cohort_summary.py:429-452`
- Test: `tests/unit/test_cohort_exports.py` or `tests/unit/test_cohort_tables.py`

**Interfaces:**
- Consumes: `additional_stats_frame` (`cohort_tables.py`), `write_cohort_frame`
  (`cohort_exports.py`) — both existing.
- Produces: `cohort_stats.{csv,tsv,json}` in the cohort output directory.

- [ ] **Step 1: Write the failing test**

```python
def test_the_cohort_statistics_export_carries_the_coverage_qc_verdict(tmp_path):
    """#172: a cohort consumer reads the CSV, not the HTML.

    `write_cohort_frame` was called for the Kestrel and adVNTR frames only; the
    statistics frame was rendered to HTML and never written, so no machine-readable
    cohort output carried any coverage at all.
    """
    frame = additional_stats_frame(
        [{"Sample": "s1", "coverage": {"mean": 250.0, "percent_uncovered": 80.0, "coverage_qc": "FAIL"}}]
    )

    write_cohort_frame(frame, tmp_path, "cohort_stats", "Statistics", ["csv"])

    written = (tmp_path / "cohort_stats.csv").read_text(encoding="utf-8")
    assert "cov_coverage_qc" in written
    assert "FAIL" in written
```

- [ ] **Step 2: Run and verify it fails**

Expected: FAIL — `cohort_stats.csv` does not exist.

- [ ] **Step 3: Implement**

In `cohort_summary.py`, hoist the statistics frame so both the HTML and the export use one
value, then add the third export beside the existing two:

```python
    additional_stats_df = additional_stats_frame(additional_stats_list) if additional_stats_list else pd.DataFrame()
    additional_stats_html = stats_table_html(additional_stats_df) if additional_stats_list else ""
    ...
    if additional_formats:
        formats = parse_output_formats(additional_formats)
        write_cohort_frame(kestrel_df, output_dir, "cohort_kestrel", "Kestrel", formats)
        write_cohort_frame(advntr_df, output_dir, "cohort_advntr", "adVNTR", formats)
        # #172: the statistics frame carries every cov_* column, coverage_qc among them.
        # Until now it reached the HTML table only, so no machine-readable cohort output
        # carried a coverage figure at all.
        write_cohort_frame(additional_stats_df, output_dir, "cohort_stats", "Statistics", formats)
```

- [ ] **Step 4: Run and verify it passes**

```bash
pytest -m unit tests/unit -o log_cli=false -q
```

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(cohort): export the per-sample statistics frame

Closes #172"
```

---

### Task 6: Lane A documentation

**Lane A. Depends on: T5.**

**Files:**
- Modify: `docs/pipeline/input-processing.md:47`, `docs/pipeline/reports.md:41-42`,
  `docs/user-guide/configuration.md:44-45`

- [ ] **Step 1: Update the three pages**

- `input-processing.md:47` — state that depth is requested with `-a` and every statistic
  is over the region, with `min = 0` wherever anything is uncovered.
- `reports.md:41-42` — the table already lists `Percent VNTR uncovered | samtools depth |
  <= 50%` as if it were enforced. It was not; #172 makes the row true. Add a sentence
  naming `coverage_qc` as the emitted verdict.
- `configuration.md:44-45` — document that **both** `thresholds` keys now decide
  `quality_metrics_pass`, and that a metric that was not measured does not fail it.

- [ ] **Step 2: Verify the docs build**

```bash
make docs-build
```
Expected: `mkdocs build --strict` succeeds. No `nav:` change is needed — all three pages
are already registered.

- [ ] **Step 3: Commit**

```bash
git add docs/
git commit -m "docs(coverage): describe region-based statistics and the two-metric QC gate

Refs #171 #172"
```

---

## Lane B — #174

### Task 7: The artifact gate, as pure logic

**Lane B. Depends on: nothing. Blocks: T8.**

**Files:**
- Modify: `vntyper/scripts/flagging.py` (append `add_artifact_gate`)
- Modify: `vntyper/scripts/kestrel_config.json` (add `artifact_flags`)
- Test: `tests/unit/test_flagging.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `add_artifact_gate(df: pd.DataFrame, artifact_flags: Sequence[str]) -> pd.DataFrame`
  — returns a copy carrying a boolean `flag_filter_pass` column on every row, always.

- [ ] **Step 1: Write the failing tests**

```python
def test_a_declared_artifact_flag_fails_the_gate():
    df = pd.DataFrame({"Flag": ["False_Positive_4bp_Insertion"]})

    assert add_artifact_gate(df, ["False_Positive_4bp_Insertion"])["flag_filter_pass"].tolist() == [False]


def test_an_advisory_flag_does_not_fail_the_gate():
    """#174 excludes artifacts, not flags.

    `Low_Depth_Conserved_Motifs` marks a call that warrants scrutiny, not one that is
    technically impossible. Excluding it would delete real low-depth calls, which is the
    failure mode this milestone exists to prevent.
    """
    df = pd.DataFrame({"Flag": ["Low_Depth_Conserved_Motifs", "Not flagged"]})

    assert add_artifact_gate(df, ["False_Positive_4bp_Insertion"])["flag_filter_pass"].tolist() == [True, True]


def test_a_row_carrying_both_an_artifact_and_an_advisory_flag_is_excluded():
    """`Flag` is comma-joined, so membership is per element."""
    df = pd.DataFrame({"Flag": ["Low_Depth_Conserved_Motifs, False_Positive_4bp_Insertion"]})

    assert add_artifact_gate(df, ["False_Positive_4bp_Insertion"])["flag_filter_pass"].tolist() == [False]


def test_a_flag_name_is_matched_whole_and_not_as_a_substring():
    """A substring test would exclude `False_Positive_4bp_Insertion_v2` given `..._Insertion`."""
    df = pd.DataFrame({"Flag": ["False_Positive_4bp_Insertion_v2"]})

    assert add_artifact_gate(df, ["False_Positive_4bp_Insertion"])["flag_filter_pass"].tolist() == [True]


def test_a_frame_without_a_flag_column_still_gains_the_gate():
    """A negative run's frame legitimately carries no `Flag` (report_formatting.py:68).

    The gate must still be present, or `filter_final_dataframe` raises on a missing
    required column - turning a safety gate into an abort.
    """
    out = add_artifact_gate(pd.DataFrame({"POS": [1, 2]}), ["False_Positive_4bp_Insertion"])

    assert out["flag_filter_pass"].tolist() == [True, True]


@pytest.mark.parametrize("value", [None, float("nan"), pd.NA])
def test_an_unknown_flag_value_is_not_a_declared_artifact(value):
    """`Flag=None` is an accepted input - test_haplo_count_and_selection.py:398 relies on it.

    Selection still deprioritises it. But absence of evidence is not evidence of an
    artifact, so it passes this gate rather than being excluded, and a bare `.split()`
    would have raised.
    """
    df = pd.DataFrame({"Flag": [value]})

    assert add_artifact_gate(df, ["False_Positive_4bp_Insertion"])["flag_filter_pass"].tolist() == [True]


def test_an_empty_frame_gains_the_gate_column():
    assert "flag_filter_pass" in add_artifact_gate(pd.DataFrame(), ["X"]).columns


def test_an_empty_artifact_list_excludes_nothing():
    """Emptying `artifact_flags` in config restores the pre-#174 behaviour with no code change."""
    df = pd.DataFrame({"Flag": ["False_Positive_4bp_Insertion"]})

    assert add_artifact_gate(df, [])["flag_filter_pass"].tolist() == [True]


def test_the_shipped_config_declares_the_4bp_insertion_as_an_artifact():
    """Ties the code to configuration: the flag name is never written inline in Python."""
    config = json.loads(Path("vntyper/scripts/kestrel_config.json").read_text(encoding="utf-8"))

    assert config["artifact_flags"] == ["False_Positive_4bp_Insertion"]
    assert "False_Positive_4bp_Insertion" in config["flagging_rules"]
```

- [ ] **Step 2: Run and verify they fail**

```bash
pytest tests/unit/test_flagging.py -m unit -o log_cli=false -q
```
Expected: FAIL — `ImportError: cannot import name 'add_artifact_gate'`.

- [ ] **Step 3: Implement**

Add `"artifact_flags": ["False_Positive_4bp_Insertion"]` to
`vntyper/scripts/kestrel_config.json`, beside `flagging_rules`.

Append to `flagging.py`:

```python
def add_artifact_gate(df: pd.DataFrame, artifact_flags: Sequence[str]) -> pd.DataFrame:
    """Mark rows carrying a declared artifact flag as failing the final filter.

    A flag is an annotation; most flags mean "this call warrants scrutiny". A few name a
    known technical artifact instead - a row that is not a candidate variant at all - and
    those must not reach ``kestrel_result.tsv``, where a consumer would read them as a
    call. Before #174, an artifact-only sample produced one flagged row that
    ``report_config.json`` mapped to ``High_Precision_flagged``, which ``is_finding``
    treats as a positive.

    Which flags are artifacts is **configuration** (``kestrel_config.json``'s
    ``artifact_flags``), never a name written into this function. Emptying that list
    restores the pre-#174 behaviour without a code change.

    This stage marks; it does not filter. ``filter_final_dataframe`` drops the rows, and
    ``kestrel_pre_result.tsv`` keeps every one of them with ``flag_filter_pass=False`` so
    the evidence for a sample is never destroyed.

    Args:
        df (pd.DataFrame): The variant frame. A missing ``Flag`` column is legitimate -
            a negative run carries none - and yields an all-True gate.
        artifact_flags (Sequence[str]): Flag names that disqualify a row.

    Returns:
        pd.DataFrame: A copy with a boolean ``flag_filter_pass`` column, present on every
        row including when the frame is empty. ``filter_final_dataframe`` raises on a
        missing gate column, so unconditional presence is the contract.

    Note:
        ``Flag`` is a comma-joined string, so membership is tested per element rather than
        by substring: an artifact named ``X`` must not exclude a flag named ``XY``. Values
        that are not strings - ``None``, ``NaN``, ``pd.NA``, all accepted inputs today -
        carry no declared artifact and therefore pass. ``select_single_best_variant``
        still deprioritises them.
    """
    result = df.copy()
    artifacts = set(artifact_flags)

    if "Flag" not in result.columns or not artifacts:
        result["flag_filter_pass"] = pd.Series(True, index=result.index, dtype=bool)
        return result

    def _is_clean(value: object) -> bool:
        if not isinstance(value, str):
            return True
        return not artifacts.intersection(part.strip() for part in value.split(","))

    result["flag_filter_pass"] = result["Flag"].map(_is_clean).astype(bool)
    return result
```

Add `from collections.abc import Sequence` to the imports.

- [ ] **Step 4: Run and verify they pass**

```bash
pytest tests/unit/test_flagging.py -m unit -o log_cli=false -q
```

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/flagging.py vntyper/scripts/kestrel_config.json tests/unit/test_flagging.py
git commit -m "feat(flagging): add the config-driven artifact gate

Refs #174"
```

---

### Task 8: Wire the sixth gate into the filter (#174)

**Lane B. Depends on: T7. Blocks: T9.**

**Files:**
- Modify: `vntyper/scripts/kestrel_genotyping.py:581-586` (step 6.5), `:826-833` (`filter_cols`)
- Modify: `tests/builders.py:33` (`STAGE_COLUMNS["flagged"]` and `["final"]`), `kestrel_row`
- Modify: `tests/unit/test_kestrel_filtering.py:44` (`GATE_COLUMNS`), `:72`
- Test: `tests/unit/test_kestrel_filtering.py`, `tests/unit/test_haplo_count_and_selection.py`

**Interfaces:**
- Consumes: `add_artifact_gate` from T7.
- Produces: `filter_cols` is six entries; `flag_filter_pass` present on every non-empty
  frame reaching `filter_final_dataframe`.

- [ ] **Step 1: Write the failing tests**

```python
def test_an_artifact_only_sample_is_not_reported_as_a_call(tmp_path):
    """#174's regression test, end to end through the real postprocessing.

    Before this change the sample yielded one `High_Precision_flagged` row, which
    `is_finding` reads as positive - a known technical artifact presented as a positive
    MUC1 call.
    """
    frame = kestrel_stage_frame("confidence", rows=1, ref="C", alt="CGGCA")

    out = process_kmer_results(frame, _merged_motifs(), str(tmp_path), _shipped_config())

    assert out.empty, "an artifact-only sample has no call"


def test_the_artifact_row_survives_in_the_pre_result_with_its_gate_false(tmp_path):
    """Evidence is never destroyed: stages mark, the final filter drops."""
    frame = kestrel_stage_frame("confidence", rows=1, ref="C", alt="CGGCA")

    process_kmer_results(frame, _merged_motifs(), str(tmp_path), _shipped_config())

    pre = pd.read_csv(tmp_path / "kestrel_pre_result.tsv", sep="\t")
    assert pre["flag_filter_pass"].tolist() == [False]
    assert "False_Positive_4bp_Insertion" in pre["Flag"].iloc[0]
```

Change `GATE_COLUMNS` at `test_kestrel_filtering.py:44` to include `"flag_filter_pass"`,
and `:72`'s assertion from `len(found) == 5` to `== 6`, renaming the test to
`test_the_gate_columns_are_exactly_the_six_this_file_pins` and recording in the docstring
that #174 added the sixth.

- [ ] **Step 2: Run and verify they fail**

```bash
pytest tests/unit/test_kestrel_filtering.py tests/unit/test_builders.py -m unit -o log_cli=false -q
```
Expected: FAIL. `test_builders.py:290` additionally fails with
`ValueError: Required filter column 'flag_filter_pass' is missing` — that is the second
tripwire and is expected here.

- [ ] **Step 3: Implement**

In `process_kmer_results` step 6.5, **outside** the existing conditional `add_flags` block
so the column is always written:

```python
    # (6.5b) #174: derive the artifact gate. Unconditional, unlike add_flags above: a
    # frame that reached the final filter without `flag_filter_pass` would abort the run
    # on a missing required gate column (#185).
    from vntyper.scripts.flagging import add_artifact_gate

    df = add_artifact_gate(df, kestrel_config.get("artifact_flags", []))
```

Add `"flag_filter_pass"` to `filter_cols` (`:827-833`), and extend the `filter_final_dataframe`
docstring's list of the boolean columns it ANDs.

Add `"flag_filter_pass"` to `STAGE_COLUMNS["flagged"]` and `STAGE_COLUMNS["final"]` in
`tests/builders.py`, and give `kestrel_row` a default of `True` for it.

- [ ] **Step 4: Run the full unit tier**

```bash
pytest -m unit tests/unit -o log_cli=false -q
```
Expected: PASS. Extend `test_haplo_count_and_selection.py:427` (`test_issue_145_scenario`)
with an artifact-flagged competitor, and add a docstring note to `:368`
(`test_all_flagged_selects_best`) that it exercises selection directly, which still
deprioritises rather than excludes — exclusion now happens upstream at the gate.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "fix(kestrel): exclude declared artifact flags at the final filter

Closes #174"
```

---

### Task 9: Lane B documentation

**Lane B. Depends on: T8.**

**Files:**
- Modify: `docs/pipeline/flagging.md:3`, `:16-20`, `:88-98`
- Modify: `docs/user-guide/configuration.md:87-88`

- [ ] **Step 1: Record the contract reversal explicitly**

`flagging.md:3` currently says flagging "identifies calls that may be technically valid
but warrant additional scrutiny", and `:88-98` describes the impact as **selection
priority only**. #174 makes that true of advisory flags and false of artifact flags.
State both classes, name `artifact_flags`, and say that an artifact-flagged row is
excluded from `kestrel_result.tsv` while remaining in `kestrel_pre_result.tsv` with
`flag_filter_pass=False`.

`configuration.md:87-88` gains the `artifact_flags` key with a note that emptying it
restores the previous behaviour.

- [ ] **Step 2: Verify and commit**

```bash
make docs-build
git add docs/
git commit -m "docs(flagging): distinguish advisory flags from artifact flags

Refs #174"
```

---

## Lane C — #203

### Task 10: Delete the dead `POS_fasta` rebase

**Lane C. Depends on: nothing. Blocks: T11.**

**Files:**
- Modify: `vntyper/scripts/motif_processing.py:490-498`
- Test: `tests/unit/test_motif_decisions.py:125`, `:140`

**Interfaces:**
- Consumes: nothing.
- Produces: no signature change. `POS_fasta` values are **unchanged** — this removes dead
  code only.

- [ ] **Step 1: Write the specification test**

Rename `test_the_original_pos_column_is_never_rebased_only_pos_fasta_is_derived`
(`:140`). Its current name asserts the opposite of the truth: `POS` is "never rebased"
only because the rebase landed on a discarded working copy, and `POS_fasta` is not
"derived" — it is a verbatim copy.

```python
    def test_pos_and_pos_fasta_are_both_the_vcf_position_in_the_pair_record(self):
        """#203: `Motif_fasta` names a 120 bp *pair* record, so `POS` is already right.

        Every record of All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa is two 60 bp
        motifs concatenated and named `<left>-<right>`; `1-2` is a record ID, `2` is not.
        `Motif_fasta` is a verbatim copy of the VCF `#CHROM`, i.e. the pair name, so the
        coordinate beside it in output.bed lives in the 120 bp record's space and the raw
        `POS` is the correct value. Subtracting `position_threshold` would make every row
        at or above 60 wrong by exactly 60 bp.

        `Motif` - the half-motif name - is what is repeat-unit-relative, and nothing
        writes a coordinate beside it.
        """
        out = motif_correction_and_annotation(_kestrel_frame(), _merged_motifs(), _shipped_config())

        assert out["POS"].tolist() == [54, 67, 67, 67, 68, 55]
        assert out.loc[1, "POS_fasta"] == 67, "not 7: no rebase, by design"
        assert out.loc[0, "POS_fasta"] == 54
```

- [ ] **Step 2: Run it and verify it passes already**

```bash
pytest tests/unit/test_motif_decisions.py -m unit -o log_cli=false -q
```
Expected: PASS. **This is the one task whose test passes before the change**, because the
change removes dead code rather than altering behaviour. The test's value is that it now
states *why* the value is what it is, so the next reader cannot re-derive the wrong
conclusion. Record that in the commit message.

- [ ] **Step 3: Delete the rebase**

`motif_processing.py:490-498` becomes:

```python
        # Adjust POS => create POS_fasta
        #
        # #203: there is nothing to adjust. `Motif_fasta` is a verbatim copy of `Motifs`,
        # which is the VCF `#CHROM` - a record of
        # All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa, every one of which is a
        # 120 bp *pair* of 60 bp motifs named `<left>-<right>`. So `POS` is already an
        # offset into the record `Motif_fasta` names, and that is exactly what
        # generate_bed_file writes beside it. A `- position_threshold` rebase here would
        # make every row at or above the threshold wrong by 60 bp.
        #
        # The code used to compute that rebase and assign it back to `POS` via
        # `DataFrame.update` - `Series.mask` keeps the name `POS`, so it never reached
        # `POS_fasta` - and `POS` is not read again. It read as working code.
        #
        # `position_threshold` keeps its live use in `split_left_right` above, which is
        # what decides *which half* of the pair name a row's `Motif` is.
        combined_df["POS"] = pd.to_numeric(combined_df["POS"], errors="coerce").fillna(-1).astype(int)
        combined_df["POS_fasta"] = combined_df["POS"]
```

- [ ] **Step 4: Run the full unit tier**

```bash
pytest -m unit tests/unit -o log_cli=false -q
```
Expected: PASS, unchanged count. **Any failure here means the rebase was not dead and the
change must be reverted.**

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/motif_processing.py tests/unit/test_motif_decisions.py
git commit -m "fix(motif): delete the POS_fasta rebase that wrote to a discarded column

Closes #203"
```

---

### Task 11: `output.bed` coordinates and the first test of `generate_bed_file`

**Lane C. Depends on: T10.**

**Files:**
- Modify: `vntyper/scripts/kestrel_genotyping.py:607-640`
- Test: `tests/unit/test_generate_bed_file.py` (create)

**Interfaces:**
- Consumes: `POS_fasta` semantics from T10.
- Produces: `generate_bed_file(df, output_dir) -> str | None` — signature unchanged;
  emitted interval shifts from `[p, p+1)` to `[p-1, p)`.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_generate_bed_file.py`:

```python
"""
Unit tests for `generate_bed_file` (#203).

This function had no test of its output content at all, which is how a coordinate
convention error survived: nothing asserted what it writes, only that it returned a path.
"""

import pandas as pd
import pytest

from vntyper.scripts.kestrel_genotyping import generate_bed_file

pytestmark = pytest.mark.unit


def test_it_writes_zero_based_half_open_intervals(tmp_path):
    """VCF POS is 1-based; BED is 0-based half-open, so POS p is the interval [p-1, p).

    The previous output was [p, p+1), naming the base *after* the variant - so IGV
    highlighted the wrong position for every row.
    """
    df = pd.DataFrame({"Motif_fasta": ["X-5"], "POS_fasta": [67]})

    generate_bed_file(df, str(tmp_path))

    assert (tmp_path / "output.bed").read_text(encoding="utf-8") == "X-5\t66\t67\n"


def test_it_names_the_pair_record_that_the_reference_fasta_actually_contains(tmp_path):
    """Column 1 must resolve against the FASTA the IGV report is given.

    `Motif_fasta` is the 120 bp pair name (`1-2`); `Motif` is the half name (`2`), which
    is not a record in that file. Writing the half name would produce a BED nothing can
    resolve.
    """
    df = pd.DataFrame({"Motif_fasta": ["1-2"], "POS_fasta": [54], "Motif": ["2"]})

    generate_bed_file(df, str(tmp_path))

    assert (tmp_path / "output.bed").read_text(encoding="utf-8").split("\t")[0] == "1-2"


def test_it_writes_one_interval_per_row(tmp_path):
    df = pd.DataFrame({"Motif_fasta": ["1-2", "X-5"], "POS_fasta": [54, 67]})

    generate_bed_file(df, str(tmp_path))

    assert (tmp_path / "output.bed").read_text(encoding="utf-8") == "1-2\t53\t54\nX-5\t66\t67\n"


@pytest.mark.parametrize("columns", [{"Motif_fasta": ["1-2"]}, {"POS_fasta": [54]}])
def test_it_returns_none_when_either_column_is_absent(columns, tmp_path):
    assert generate_bed_file(pd.DataFrame(columns), str(tmp_path)) is None
    assert not (tmp_path / "output.bed").exists()


def test_it_returns_none_for_an_empty_frame(tmp_path):
    df = pd.DataFrame({"Motif_fasta": [], "POS_fasta": []})

    assert generate_bed_file(df, str(tmp_path)) is None
```

- [ ] **Step 2: Run and verify they fail**

```bash
pytest tests/unit/test_generate_bed_file.py -m unit -o log_cli=false -q
```
Expected: FAIL — the first test reports `'X-5\t67\t68\n' != 'X-5\t66\t67\n'`.

- [ ] **Step 3: Implement**

`kestrel_genotyping.py:632-636`:

```python
    # Each row: motif_fasta, start=pos_fasta-1, end=pos_fasta.
    #
    # `POS_fasta` is the 1-based VCF position within the 120 bp pair record named by
    # `Motif_fasta` (#203). BED intervals are 0-based and half-open, so position p is the
    # interval [p-1, p). This used to write [p, p+1), which named the base after the
    # variant - IGV highlighted the wrong position on every row.
    with open(bed_file_path, "w") as bed_file:
        for _, row in df.iterrows():
            motif_fasta = row["Motif_fasta"]
            pos = int(row["POS_fasta"])
            bed_file.write(f"{motif_fasta}\t{pos - 1}\t{pos}\n")
```

- [ ] **Step 4: Run and verify they pass**

```bash
pytest -m unit tests/unit -o log_cli=false -q
```

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/kestrel_genotyping.py tests/unit/test_generate_bed_file.py
git commit -m "fix(kestrel): write output.bed as 0-based half-open intervals

Refs #203"
```

---

## Lane D — #212

### Task 12: Remove the unconditional Kestrel skip

**Lane D. Depends on: nothing. Blocks: T13.**

**Files:**
- Modify: `vntyper/scripts/kestrel_genotyping.py:294-315`
- Test: `tests/unit/test_run_kestrel_skip.py` (create)

**Interfaces:**
- Consumes: nothing.
- Produces: `run_kestrel(...) -> None` — signature unchanged; now raises `RuntimeError`
  when no k-mer size produced a result.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_run_kestrel_skip.py`:

```python
"""
Unit tests for run_kestrel's handling of a pre-existing output.vcf (#212).

The skip was unconditional, gated on no flag, undocumented and untested. Because it
`return`ed rather than falling through, it also skipped the two statements that turn a
VCF into a result - so re-running into a directory left by an interrupted run could
report a confident negative for a sample carrying a pathogenic variant.
"""
```

```python
def test_a_stale_vcf_does_not_skip_the_kestrel_run(tmp_path, monkeypatch):
    """The regression test #212 asks for: assert Kestrel ran, never that it silently didn't."""
    (tmp_path / "output.vcf").write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    executed = []
    monkeypatch.setattr(kg, "run_command", lambda cmd, log, **kw: executed.append(cmd) or True)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)
    monkeypatch.setattr(kg, "process_kestrel_output", lambda *a, **k: None)

    kg.run_kestrel(tmp_path / "output.vcf", str(tmp_path), "r1.fq", "r2.fq",
                   "ref.fa", "kestrel.jar", _config(), "sample")

    assert any("kestrel" in c.lower() for c in executed), "Kestrel was skipped"


def test_a_stale_vcf_is_removed_before_the_run(tmp_path, monkeypatch, caplog):
    """Unlinking makes the reason legible in the log rather than surfacing as a Java error."""
    vcf = tmp_path / "output.vcf"
    vcf.write_text("stale\n", encoding="utf-8")
    seen = {}
    monkeypatch.setattr(kg, "run_command", lambda cmd, log, **kw: seen.setdefault("existed", vcf.is_file()) or True)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: None)
    monkeypatch.setattr(kg, "process_kestrel_output", lambda *a, **k: None)

    with caplog.at_level(logging.WARNING):
        kg.run_kestrel(vcf, str(tmp_path), "r1.fq", "r2.fq", "ref.fa", "kestrel.jar", _config(), "sample")

    assert seen["existed"] is False, "the stale VCF was still present when Kestrel ran"
    assert "output.vcf" in caplog.text


def test_a_run_that_produces_no_vcf_raises_rather_than_returning_silently(tmp_path, monkeypatch):
    """Outcome 2 of #212, closed at the source.

    A Kestrel invocation that exits 0 without writing a VCF used to fall out of the loop
    and return None. `record_step` then wrote md5sum=None and an empty data list, which
    both the HTML report and cohort mode render as a negative.
    """
    monkeypatch.setattr(kg, "run_command", lambda *a, **kw: True)

    with pytest.raises(RuntimeError, match="produced no VCF"):
        kg.run_kestrel(tmp_path / "output.vcf", str(tmp_path), "r1.fq", "r2.fq",
                       "ref.fa", "kestrel.jar", _config(), "sample")


def test_post_processing_runs_after_a_successful_kestrel_run(tmp_path, monkeypatch):
    """The two statements the old `return` skipped are what produce kestrel_result.tsv."""
    vcf = tmp_path / "output.vcf"
    calls = []
    monkeypatch.setattr(kg, "run_command", lambda *a, **kw: vcf.write_text("ok\n", encoding="utf-8") or True)
    monkeypatch.setattr(kg, "convert_sam_to_bam_and_index", lambda *a, **k: calls.append("bam"))
    monkeypatch.setattr(kg, "process_kestrel_output", lambda *a, **k: calls.append("post"))

    kg.run_kestrel(vcf, str(tmp_path), "r1.fq", "r2.fq", "ref.fa", "kestrel.jar", _config(), "sample")

    assert calls == ["bam", "post"]
```

- [ ] **Step 2: Run and verify they fail**

```bash
pytest tests/unit/test_run_kestrel_skip.py -m unit -o log_cli=false -q
```
Expected: FAIL — the first test finds no executed Kestrel command.

- [ ] **Step 3: Implement**

Replace `kestrel_genotyping.py:294-297` and add the completion guard:

```python
    completed = False

    for kmer_size in kmer_sizes:
        ...
        # #212: a pre-existing output.vcf used to skip the whole stage and `return`,
        # which also skipped the two statements that turn a VCF into a result. Re-running
        # into a directory left by an interrupted run then either re-reported a stale
        # result or produced none at all - and a step with no result file is rendered as
        # a negative. Deliberate reuse belongs behind the --resume flag proposed in #20;
        # this ad-hoc, unflagged version is unsafe. The stale file is removed rather than
        # overwritten so the log says why.
        if vcf_path.is_file():
            logger.warning(f"Removing a VCF left by an earlier run before re-running Kestrel: {vcf_path}")
            vcf_path.unlink()

        logger.info(f"Launching Kestrel with k-mer size {kmer_size}...")
        ...
                process_kestrel_output(output_dir, vcf_path, reference_vntr, kestrel_config, config)
                completed = True
                break

    if not completed:
        msg = (
            "Kestrel produced no VCF for any configured k-mer size, so no result file was written. "
            "Reporting this as a negative would manufacture a confident negative genotype. See issue #212."
        )
        logger.error(msg)
        raise RuntimeError(msg)
```

Update the docstring's step 2 ("If the final VCF already exists, skip") and its `Raises:`.

- [ ] **Step 4: Run the full unit tier**

```bash
pytest -m unit tests/unit -o log_cli=false -q
```

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/kestrel_genotyping.py tests/unit/test_run_kestrel_skip.py
git commit -m "fix(kestrel): never skip the stage because output.vcf exists

Closes #212"
```

---

### Task 13: A missing result file is loud (#212)

**Lane D. Depends on: T12.**

**Files:**
- Modify: `vntyper/scripts/summary.py:209-236`
- Test: `tests/unit/test_summary_record_step.py` (create or extend the existing summary tests)

**Interfaces:**
- Consumes: nothing.
- Produces: `record_step(...)` adds `record["result_file_missing"] = True` **only** when the
  file is absent.

- [ ] **Step 1: Write the failing tests**

```python
def test_record_step_flags_a_missing_result_file(tmp_path, caplog):
    """#212's second half. parse_tsv swallows the FileNotFoundError into a comment and
    returns an empty `data` list, so a step that produced nothing is indistinguishable
    from a step that legitimately found nothing."""
    summary = {"steps": []}

    with caplog.at_level(logging.ERROR):
        record_step(summary, "Kestrel Genotyping", str(tmp_path / "absent.tsv"), "tsv", "cmd", _t0(), _t1())

    assert summary["steps"][0]["result_file_missing"] is True
    assert "absent.tsv" in caplog.text


def test_record_step_adds_no_key_when_the_file_exists(tmp_path):
    """Pins gate-neutrality: the golden harness compares whole step records, so an
    unconditional key would diff every case of every run."""
    path = tmp_path / "present.tsv"
    path.write_text("a\tb\n1\t2\n", encoding="utf-8")
    summary = {"steps": []}

    record_step(summary, "Kestrel Genotyping", str(path), "tsv", "cmd", _t0(), _t1())

    assert "result_file_missing" not in summary["steps"][0]
```

- [ ] **Step 2: Run and verify they fail**

- [ ] **Step 3: Implement**

In `record_step`, after the record dict is built and before the MD5:

```python
    # #212: a step whose result file is absent is a failure, not an empty result.
    # `parse_tsv` turns FileNotFoundError into a comment and returns `data: []`, and both
    # the HTML report and cohort mode render an empty `data` as a negative. The key is
    # added only when the file is missing, so a normal run's summary - and the
    # golden-cohort `pipeline_step_records` artefact - is byte-identical.
    if not os.path.exists(result_file):
        logger.error(f"Step '{step_name}' recorded a result file that does not exist: {result_file}")
        record["result_file_missing"] = True
```

- [ ] **Step 4: Run and verify they pass**

```bash
pytest -m unit tests/unit -o log_cli=false -q
```

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/summary.py tests/unit/test_summary_record_step.py
git commit -m "fix(summary): record a missing result file loudly instead of as no data

Closes #212"
```

---

## Integration

### Task 14: Rebase the lanes and update the benchmark note

**Depends on: T6, T9, T11, T13.**

- [ ] **Step 1: Rebase each lane onto the milestone branch, in order**

Order is A, B, C, D — A is the largest and B/C/D touch `kestrel_genotyping.py` in disjoint
regions, so conflicts should not arise. Resolve any by keeping both changes; they are
independent edits to different functions.

```bash
cd /home/bernt-popp/development/VNtyper
for lane in a b c d; do
  branch=$(git -C .claude/worktrees/m2-lane-$lane branch --show-current)
  git -C .claude/worktrees/m2-lane-$lane rebase fix/milestone-2-correctness-of-reported-numbers
  git merge --ff-only "$branch"
done
git log --oneline --graph -25
```

- [ ] **Step 2: Update `tests/benchmark/benchamrk_downsample.py:172-220`**

It carries a private copy of `calculate_vntr_coverage` with the covered-position mean and
a depth command without `-a`. Nothing gates it, but it becomes a misleading benchmark.
Either update it to match, or add a comment naming #171 and stating that it deliberately
reproduces the pre-#171 arithmetic for historical comparison. Do not leave it silently
divergent.

- [ ] **Step 3: Run the whole gate**

```bash
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vntyper \
  && export PATH="$CONDA_PREFIX/bin:$PATH" \
  && make check-all && make test-unit-cov && make patch-coverage
```
Expected: all three green. Paste the output. If `patch-coverage` is below 80%, add tests —
never lower the target.

- [ ] **Step 4: Commit**

```bash
git commit -am "test(benchmark): align the downsample benchmark with the #171 arithmetic"
```

---

### Task 15: Version, changelog and the golden-cohort gate

**Depends on: T14.**

- [ ] **Step 1: Bump the version in both places**

`vntyper/version.py` and `CITATION.cff` to `2.1.0` — #171, #172 and #174 change reported
values and the QC verdict, which is more than a patch.
`tests/unit/test_version_consistency.py` fails the build if they disagree.
**Do not create or push a `v*.*.*` tag.**

- [ ] **Step 2: Write the changelog**

`docs/about/changelog.md`, one entry per issue under a `2.1.0` heading. State plainly that
coverage numbers change, that #171's region harmonisation is *not* included, and that
`docs/pipeline/reports.md` previously documented a threshold that was not enforced.

- [ ] **Step 3: Run the golden-cohort gate**

```bash
python scripts/golden_cohort_gate.py run --advntr-max-coverage 300 ...
python scripts/golden_cohort_gate.py compare ... --expect-after-sha "$(git rev-parse HEAD)" --require-clean
```

Follow `docs/development/golden-cohort-gate.md`. Baseline is `cb593b6`. Every delta must
be attributed to a named change from SPEC.md, and every delta SPEC.md predicts must be
observed or explained. The predictions to check against:

| Artefact | Expected | From |
| --- | --- | --- |
| `coverage_summary` | values move + one new column, every case | #171, #172 |
| `executed_commands` | `-a`, every case | #171 |
| `kestrel_result`, `kestrel_pre_result` | `columns_added: flag_filter_pass`, every case reaching the filter | #174 |
| `kestrel_result` | rows removed on artifact-selected cases | #174 |
| `advntr_result` | **possible genotype delta** on 300× cases | #171 |
| `screening_summary` | where a corrected mean or the new uncovered rule crosses a threshold | #171, #172 |
| cohort output file set | `cohort_stats.{csv,tsv,json}` added | #172 |
| `pipeline_step_records` | **no delta** | — |
| anything attributable to #203 or #212 | **no delta** | — |

- [ ] **Step 4: Write up the run**

Add a row to the run table in `docs/development/golden-cohort-gate.md` and a result
section attributing every delta. Note that `output.bed` is not a gate artefact, so the
gate is **not** evidence for the #203 rider — `tests/unit/test_generate_bed_file.py` is.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "chore(release): 2.1.0"
```

---

### Task 16: One PR, and the issue bookkeeping

**Depends on: T15.**

- [ ] **Step 1: Verify before claiming**

Use superpowers:verification-before-completion. Re-run and paste:
`make check-all`, `make test-unit-cov`, `make patch-coverage`.

- [ ] **Step 2: Open one PR**

Not stacked. Body closes all five:
`Closes #171`, `Closes #172`, `Closes #174`, `Closes #203`, `Closes #212`.
State explicitly what is **not** included: #171's region harmonisation.

- [ ] **Step 3: Issue bookkeeping**

- Comment on #171 listing its four suggested items and which three are done.
- File the follow-up: "Harmonise the hg19/hg38 VNTR region spans (1501 vs 4501 bp)",
  naming the flank-convention question as the decision it needs.
- Use `gh api repos/hassansaei/VNtyper/...` — `gh issue view` 500s on this repo.

- [ ] **Step 4: Clean up the worktrees**

```bash
for lane in a b c d; do git worktree remove .claude/worktrees/m2-lane-$lane; done
git worktree prune
```

---

## Self-review

**Spec coverage.** Every SPEC.md section maps to a task: #171 → T1, T2 (+ T6 docs);
#172 → T3, T4, T5 (+ T6); #174 → T7, T8 (+ T9); #203 → T10, T11; #212 → T12, T13;
cross-cutting (version, changelog, docs, gate, PR) → T14–T16. The adversarial-review
resolutions A1 (rounding) → T3/T4, A3 (`kestrel_result` column) → T8 + T15 predictions,
A4 (cohort export) → T5, A7 (coverage gates in DoD) → T14, A13 (benchmark) → T14.

**Placeholders.** None. Every code step carries real code; every test step carries real
assertions. T15's gate command is elided to `...` only where
`docs/development/golden-cohort-gate.md` is the authority on the exact invocation, and it
is named as such.

**Type consistency.** `evaluate_coverage_qc` and `CoverageQC` are defined in T3 and used
with the same names in T4 and T5's tests. `add_artifact_gate(df, artifact_flags)` is
defined in T7 and called with that signature in T8. `flag_filter_pass` is spelled
identically in T7, T8, `tests/builders.py` and T15's prediction table. `COVERAGE_COLUMNS`
gains `"coverage_qc"` in T4 and is asserted with that spelling in T4's test and
`tests/helpers.py`.

**One known ordering hazard.** T10's test passes *before* its change — the only such task
in the plan, because #203 deletes dead code. It is called out in the task itself so an
executor following TDD does not conclude the plan is broken.
