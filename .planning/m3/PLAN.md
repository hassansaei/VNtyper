# Milestone 3 — Web service and cohort integrity: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:subagent-driven-development`.
> Steps use checkbox (`- [ ]`) syntax. Read `.planning/m3/SPEC.md` before starting your
> workstream, and `AGENTS.md` before touching any file.
>
> **`.planning/m3/REVIEW.md` supersedes this document wherever the two conflict.** This plan was
> adversarially reviewed after it was written; 15 findings were verified against the code and
> accepted or rejected there, and the corrections were not all folded back into the task bodies
> below. Read REVIEW.md for your workstream **before** starting, and treat its "Action:" lines as
> binding. The tasks below are still the right decomposition and the right order.
>
> **Subagents run no git commands at all.** Not `add`, not `commit`, not `checkout`. Four
> workstreams share one worktree, so `.git/index`, `index.lock` and the branch ref are shared,
> and concurrent staging can capture another stream's files into the wrong commit
> (REVIEW.md finding 15). Edit files, run pytest, report. The integrator stages and commits each
> workstream sequentially. The `git commit` blocks in the tasks below are the **integrator's**
> checklist — they record the intended message and file set, nothing more.

**Goal:** Close GitHub milestone 3 (#216 #207 #206 #205 #201 #167 #162) plus #210, so that a web
job leaves nothing behind, two cohort runs are byte-comparable, and no sample-derived string
reaches the DOM unescaped.

**Architecture:** Four workstreams over disjoint file sets, executed concurrently on one branch,
each strictly TDD. A — move every sample-derived value in the two templates onto DOM APIs and
re-serialise the two IGV fragments through `json.dumps` in Python. B — give every write that
currently derives its path from the input alignment an output-directory destination. C — widen
the cohort pseudonym, make the map injective-or-loud, and give ZIP samples an identity and an
order that do not contain a temp path. D — resolve the VCF `vntyper report` already has.

**Tech Stack:** Python 3.10 floor, ruff (line length 120, double quotes), mypy, pytest with
`-m unit`, Jinja2 templates, jQuery 3.5 + DataTables 1.10 in the templates.

## Global Constraints

Copied verbatim from `AGENTS.md`; every task's requirements include this section.

- **Branch:** `fix/milestone-3-web-cohort-integrity`, off `main`. One branch, one PR, never
  stacked. No force-push.
- **Commits:** Conventional Commits, `type(scope): lowercase subject` under ~72 chars, `Closes #N`
  in the footer. Every commit becomes permanent history — PRs are not squashed.
- **Python floor is 3.10.** PEP 604 (`X | None`) and PEP 585 (`list[str]`) are available; nothing
  newer is.
- **Ruff is the only formatter and linter.** Line length **120**, double quotes,
  `target-version = "py310"`. Never reformat to 88.
- **Every new unit test file must declare `pytestmark = pytest.mark.unit`** after its imports.
  CI runs `pytest -m unit`; an unmarked file silently never runs, and
  `tests/unit/test_marker_hygiene.py` fails the build naming it.
- **Unit tier only.** `tmp_path` + `unittest.mock`. No network, no Docker, no reference data, no
  1.1 GB Zenodo archive.
- **Run pytest from the repo root.** `tests/parametrization.py` opens
  `tests/test_data_config.json` by relative path at collection time; any other CWD breaks
  collection, including `-m unit`.
- **Google-style docstrings** (`Args:` / `Returns:` / `Raises:`) on public functions.
- **Logging:** `logger = logging.getLogger(__name__)` after imports; `logger.info(...)`, never
  `logging.info(...)`, never `basicConfig`. f-strings in log calls are the house style.
- **Errors:** no custom exception classes. `logger.error(msg)` then
  `raise ValueError(msg)` / `RuntimeError(msg)` with the same message.
- **Config-driven, never hardcoded.** Thresholds and rule tables live in the JSON configs.
  `--config-path` **replaces** the whole config rather than merging (AGENTS.md trap 2), so every
  new key is read through a `.get()` chain against a module-level default constant — never
  `config["x"]["y"]`.
- **Patch coverage ≥ 80%** on changed lines (`make patch-coverage`). Write tests that kill
  mutants: assert on values, not that a call returned.
- **Never** commit into `tests/data/`, `reference/`, `out/`, `download/`; never push a `v*.*.*`
  tag; never add a `docs/` page without a `mkdocs.yml` `nav:` entry.
- **Subagents never touch** `vntyper/version.py`, `CITATION.cff`, `docs/about/changelog.md`,
  `AGENTS.md`, `pyproject.toml`, or another workstream's files. The integrator owns those.

---

## File structure — the four sets, and their disjointness

| WS | Source files owned | Test files owned |
| --- | --- | --- |
| **A** | `vntyper/scripts/report_formatting.py`, `vntyper/scripts/generate_report.py`, `vntyper/templates/report_template.html`, `vntyper/templates/cohort_summary_template.html` | `tests/unit/test_report_formatting.py`, `tests/unit/test_generate_report.py`, `tests/unit/test_template_escaping.py` *(new)* |
| **B** | `vntyper/scripts/utils.py`, `vntyper/scripts/pipeline.py`, `vntyper/scripts/fastq_bam_processing.py`, `vntyper/scripts/command_builders.py`, `vntyper/scripts/alignment_processing.py` | `tests/unit/test_utils.py`, `tests/unit/test_command_builders.py`, `tests/unit/test_fastq_bam_command_wiring.py`, `tests/unit/test_pipeline_cwd.py`, `tests/unit/web/test_tasks.py`, `tests/unit/web/test_index_handoff.py`, `tests/docker/conftest.py`, `tests/unit/test_input_tree_is_never_written.py` *(new)* |
| **C** | `vntyper/scripts/cohort_inputs.py`, `vntyper/scripts/cohort_summary.py`, `vntyper/config.json`, `scripts/golden_cohort/compare.py` | `tests/unit/test_cohort_inputs.py`, `tests/unit/test_cohort_summary_oracle.py`, `tests/unit/test_golden_cohort_compare.py`, `tests/unit/test_cohort_identity.py` *(new)* |
| **D** | `vntyper/scripts/cli_report.py`, `vntyper/scripts/cli_parser.py` | `tests/unit/test_cli_report.py`, `tests/unit/test_cli_parser_contract.py` |

**Proof of disjointness.** Sort every path above and check for repeats — there are none. The four
source sets are pairwise disjoint (A∩B=A∩C=A∩D=B∩C=B∩D=C∩D=∅) and so are the four test sets. Two
adjacencies that are *not* overlaps:

- A and D both reach `generate_summary_report`. A edits its two template-context lines
  (`generate_report.py:587-588`); D edits its **caller** (`cli_report.py:128-139`). Different
  files.
- B and C both touch code `pipeline.py` calls, but only B edits `pipeline.py`.

**Two shared files are read-only for everyone:** `tests/support/pipeline_harness.py` and
`tests/builders.py`. If a workstream believes it must edit either, it stops and reports to the
integrator instead — those feed other workstreams' tests.

**`tests/unit/test_cohort_summary_oracle.py` is assigned to C, and A must not edit it.** Both
workstreams reach it, so the tie is broken explicitly rather than left to whoever gets there
first:

- C **must** edit it. Line 862-863 pins a literal MD5-derived pseudonym —
  `assert ">anon_65622</td>" in html` and `assert "anon_65622\tsample_one" in table` — which
  Task C1 changes by construction.
- A **must not**. A1 changes only the script body of `cohort_summary_template.html`, and the
  oracle's `_skeleton()` substitutes `<SCRIPT-BODY>` for every script body before hashing, so
  `EXPECTED_FINGERPRINT` (line 174) must not move. **If A finds the fingerprint has moved, A
  stops and reports to the integrator** — a moved fingerprint means the change reached non-script
  markup, which A1 is not supposed to do.

**Three more source-text tripwires fire on C's changes, in files C therefore owns:**

- `tests/unit/test_golden_cohort_compare.py:219` asserts the literal string
  `"return sorted(processed_dirs), temp_dirs"` is present in `cohort_inputs.py`. Task C3 changes
  that return statement.
- `tests/unit/test_golden_cohort_compare.py:230` asserts
  `"ZIP input extracts to a randomly-named tempdir"` is in `compare.COHORT_ORDER_WHY`. After C3
  that justification is **false** — the temp path is no longer in the sort key, which is the
  whole point of the fix.
- `tests/unit/test_golden_cohort_compare.py:200-210` asserts every `test_*` name cited in
  `COHORT_ORDER_WHY` still exists in `test_cohort_inputs.py`, so renaming a cited test there
  fails this file too.

Task **C4** closes all three.

**Files no workstream expects to edit but must re-run:** `tests/unit/web/test_tasks.py` already
has `test_a_successful_job_is_marked_completed_and_cleans_up_its_inputs` (line 454) and
`test_cleanup_leaves_a_nonempty_input_directory_and_warns_instead_of_deleting_it` (line 683).
B must confirm both still pass; if either needs an edit, B owns it and says so in the commit.

---

# Workstream A — XSS (#207, #216)

**Branch discipline:** commit after every task. Scope for commits: `report`, `templates`.

### Task A1: Kill the flag-tooltip concatenation in both templates

**Files:**
- Modify: `vntyper/templates/cohort_summary_template.html:152-157`
- Modify: `vntyper/templates/report_template.html:421-425`
- Test: `tests/unit/test_template_escaping.py` *(create)*

**Interfaces:**
- Consumes: nothing.
- Produces: nothing importable. The contract is the absence of a source-text pattern.

**Why a source-text test.** The unit tier cannot execute JavaScript (no network, no browser, no
node). A source-text tripwire is the honest gate and this repo already uses the pattern —
`tests/unit/test_kestrel_filtering.py` reads `kestrel_genotyping.py` as text and asserts on it.
The test's docstring must say plainly that it is a tripwire, not a behavioural test.

- [ ] **Step 1: Write the failing test**

```python
"""Neither report template may build markup by concatenating a sample-derived value.

SPECIFICATION (#207): the server escapes every ``Flag`` cell, and then
``updateFlagColumn()`` used to undo that escaping -- ``.text()`` returns the
*decoded* value, the value was concatenated into a ``title="..."`` attribute, and
``.html()`` reparsed the result. A ``Flag`` of
``"></span><img src=x onerror=alert(1)><span title="`` renders safely in the
initial document and executes on the first DataTables redraw.

This is a **tripwire over the template source text**, not a behavioural test. The
unit tier has no JavaScript engine, so nothing here can prove the page is safe by
running it; what it can prove is that the sink which made it unsafe is gone and
has not come back. The behavioural evidence lives in the Phase 5 rendered-HTML
grep recorded in ``.planning/m3/EVIDENCE.md``.
"""

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

TEMPLATES = (
    Path("vntyper/templates/cohort_summary_template.html"),
    Path("vntyper/templates/report_template.html"),
)

#: ``.html('<... ' + something + ' ...>')`` -- a variable spliced into a markup
#: string that is then parsed. This is the #207 sink in its exact shape.
_HTML_CONCAT = re.compile(r"\.html\(\s*'[^']*'\s*\+")


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_no_template_concatenates_a_variable_into_parsed_markup(template: Path) -> None:
    source = template.read_text(encoding="utf-8")
    offenders = _HTML_CONCAT.findall(source)
    assert offenders == [], f"{template} splices a variable into markup that .html() reparses: {offenders}"


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_the_flag_tooltip_title_is_set_through_attr(template: Path) -> None:
    source = template.read_text(encoding="utf-8")
    assert ".attr('title', originalText)" in source, (
        f"{template} must set the tooltip title with .attr(), which takes a value rather than a fragment"
    )


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_the_flag_cell_is_emptied_before_the_mark_is_appended(template: Path) -> None:
    """``.empty().append()`` replaces the cell without an HTML parse."""
    source = template.read_text(encoding="utf-8")
    assert "$flagCell.empty().append(" in source, f"{template} must rebuild the flag cell with DOM APIs"
```

- [ ] **Step 2: Run it and watch it fail**

Run: `python -m pytest tests/unit/test_template_escaping.py -v`
Expected: FAIL — `test_no_template_concatenates_a_variable_into_parsed_markup` finds two
offenders per template, and the two positive assertions find neither string.

- [ ] **Step 3: Rewrite both `updateFlagColumn` bodies**

Replace the `if/else` that calls `$flagCell.html(...)` in **both** files with:

```javascript
          // #207: the title is a sample-derived value. .attr() takes a value and
          // .text() takes a text node, so neither goes through an HTML parse --
          // which is what undid the server's escaping when this was built by
          // concatenating into a markup string and calling .html().
          var isClean = (originalText === "Not flagged" || originalText === "Not applicable" || originalText === "");
          var $mark = $('<span>')
              .attr('data-toggle', 'tooltip')
              .attr('title', originalText)
              .css({ 'color': isClean ? 'green' : 'red', 'font-size': '16px' })
              .text(isClean ? '\u2713' : '\u2716');
          $flagCell.empty().append($mark);
```

Indent to match the surrounding file: the cohort template uses 10 spaces at this depth, the
per-sample template uses 20. Do **not** touch the `data-original` block above it — it is already
`.attr()`-set, the DataTables search filter at `cohort_summary_template.html:103` reads it, and
it is what keeps the function idempotent across redraws.

- [ ] **Step 4: Run the test and the two suites it could disturb**

Run:
```bash
python -m pytest tests/unit/test_template_escaping.py tests/unit/test_generate_report.py \
  tests/unit/test_cohort_summary_oracle.py tests/unit/test_cohort_summary_escaping.py -v
```
Expected: all PASS. The oracle fingerprint substitutes `<SCRIPT-BODY>` before hashing, so a
script-body-only change must not move it. **If the oracle fails, stop and report** — that means
the change reached non-script markup, which this task must not do.

- [ ] **Step 5: Commit**

```bash
git add tests/unit/test_template_escaping.py vntyper/templates/cohort_summary_template.html \
        vntyper/templates/report_template.html
git commit -m "fix(templates): build the flag tooltip with DOM APIs, not markup

Closes #207"
```

---

### Task A2: Re-serialise the IGV fragments as script-safe JSON

**Files:**
- Modify: `vntyper/scripts/report_formatting.py:403-419` (replace `js_object_literal`)
- Modify: `vntyper/scripts/generate_report.py:42-49` (import), `:587-588` (call sites)
- Test: `tests/unit/test_report_formatting.py`

**Interfaces:**
- Consumes: `EMPTY_TABLE_JSON`, `EMPTY_SESSION_DICTIONARY` from `report_formatting.py:140-141`
  (unchanged).
- Produces: `js_json_literal(fragment: str, fallback: str) -> str`, replacing
  `js_object_literal(fragment: str, fallback: str) -> str`. `generate_report.py` is its only
  production caller.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_report_formatting.py` (the file already carries
`pytestmark = pytest.mark.unit`; do not add a second one):

```python
# ---------------------------------------------------------------------------
# js_json_literal (#216)
# ---------------------------------------------------------------------------
#
# SPECIFICATION: `report_template.html` interpolates these two fragments straight
# into a <script> block with `|safe`. They are lifted verbatim out of another
# tool's HTML by `extract_line_after`, and the only guard they used to get was
# `js_object_literal`, which checked that the fragment was non-empty. That is a
# *syntax* guard against `const tableJson = ;`, not a safety guard: a `</script>`
# anywhere in the fragment closes the block early and everything after it is
# parsed as HTML. The fragment is sample-derived -- it is the igv-reports variant
# table, built from VNtyper's own BED and VCF, whose REF, ALT and Motif_sequence
# values come from the sample's reads.


def test_a_well_formed_fragment_round_trips_as_json() -> None:
    result = js_json_literal('{"headers": ["a"], "rows": [[1]]}', EMPTY_TABLE_JSON)
    assert json.loads(result) == {"headers": ["a"], "rows": [[1]]}


def test_the_trailing_statement_terminator_is_stripped() -> None:
    """``extract_line_after`` returns the rest of the line, semicolon included."""
    result = js_json_literal('{"a": 1};', EMPTY_TABLE_JSON)
    assert json.loads(result) == {"a": 1}


def test_a_script_close_in_a_string_value_cannot_terminate_the_block() -> None:
    fragment = json.dumps({"rows": [["</script><img src=x onerror=alert(1)>"]]})
    result = js_json_literal(fragment, EMPTY_TABLE_JSON)
    assert "</" not in result
    assert "<\\/script>" in result
    # It is still the same data once the browser parses it.
    assert json.loads(result.replace("<\\/", "</"))["rows"][0][0].startswith("</script>")


@pytest.mark.parametrize("separator", ["\u2028", "\u2029"])
def test_javascript_line_terminators_are_escaped(separator: str) -> None:
    """U+2028 and U+2029 end a line in JavaScript but are legal inside a JSON string."""
    result = js_json_literal(json.dumps({"a": separator}), EMPTY_TABLE_JSON)
    assert separator not in result


@pytest.mark.parametrize(
    "fragment",
    ["", "   ", "not json at all", "{unquoted: 1}", "{'single': 'quotes'}", "{"],
)
def test_a_fragment_that_is_not_json_falls_back(fragment: str) -> None:
    assert js_json_literal(fragment, EMPTY_TABLE_JSON) == EMPTY_TABLE_JSON


def test_the_fallback_is_returned_verbatim_for_the_session_dictionary() -> None:
    assert js_json_literal("", EMPTY_SESSION_DICTIONARY) == EMPTY_SESSION_DICTIONARY


def test_the_output_is_deterministic_for_the_same_input() -> None:
    """Two runs over the same IGV page must emit byte-identical script (exit criterion E2)."""
    fragment = '{"b": 2, "a": 1}'
    assert js_json_literal(fragment, EMPTY_TABLE_JSON) == js_json_literal(fragment, EMPTY_TABLE_JSON)
    assert js_json_literal(fragment, EMPTY_TABLE_JSON) == '{"a":1,"b":2}'


def test_a_malformed_fragment_is_logged_at_warning(caplog) -> None:
    with caplog.at_level(logging.WARNING):
        js_json_literal("not json", EMPTY_TABLE_JSON)
    assert any("could not be parsed as JSON" in r.getMessage() for r in caplog.records)
```

Add `import json` and `import logging` to the file's imports if absent, and add
`js_json_literal` to the existing `from vntyper.scripts.report_formatting import (...)` block.
Delete any existing `js_object_literal` import and test — there are none today (grep confirms the
function has **no test at all**, which is how #216 survived).

- [ ] **Step 2: Run and watch it fail**

Run: `python -m pytest tests/unit/test_report_formatting.py -v`
Expected: FAIL at collection — `ImportError: cannot import name 'js_json_literal'`.

- [ ] **Step 3: Implement**

In `report_formatting.py`, replace `js_object_literal` entirely:

```python
def js_json_literal(fragment: str, fallback: str) -> str:
    """Re-serialise an extracted fragment as a literal that is safe in a ``<script>``.

    ``report_template.html`` interpolates the return value directly into a script
    block as ``const tableJson = {{ table_json|safe }};``. The fragment reaching
    here was lifted verbatim out of the igv-reports page by
    :func:`extract_line_after` and is sample-derived, so it is re-parsed and
    re-emitted rather than trusted: what the template receives is always the
    output of :func:`json.dumps`, never the extracted text.

    Three characters are escaped afterwards because they are legal inside a JSON
    string but not inside a script block: ``</`` would close the element early
    (``</script>`` in any string value used to terminate the block and turn
    everything after it into HTML), and U+2028/U+2029 are line terminators to a
    JavaScript parser.

    Keys are sorted and separators are minimised so that two runs over the same
    IGV page emit byte-identical script.

    Args:
        fragment: The extracted literal, possibly empty, possibly not JSON.
        fallback: A valid literal to use when the fragment cannot be parsed. An
            empty fragment would otherwise produce ``const tableJson = ;`` -- a
            syntax error that takes the whole script block down, and with it the
            variant table, the flag toggles and the coverage switch.

    Returns:
        str: A JSON literal that parses as JavaScript and cannot escape the
        script block.
    """
    candidate = fragment.strip().rstrip(";").strip()
    if not candidate:
        return fallback
    try:
        value = json.loads(candidate)
    except ValueError as e:
        logger.warning(f"IGV fragment could not be parsed as JSON and was discarded: {e}")
        return fallback

    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return encoded.replace("</", "<\\/").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
```

Add `import json` to `report_formatting.py` if it is not already imported.

In `generate_report.py`, change the import at line 49 from `js_object_literal` to
`js_json_literal` and both call sites at 587-588:

```python
        "table_json": js_json_literal(table_json, EMPTY_TABLE_JSON),
        "session_dictionary": js_json_literal(session_dictionary, EMPTY_SESSION_DICTIONARY),
```

- [ ] **Step 4: Run**

Run: `python -m pytest tests/unit/test_report_formatting.py tests/unit/test_generate_report.py -v`
Expected: PASS. Then `grep -rn "js_object_literal" vntyper/ tests/` must return nothing.

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/report_formatting.py vntyper/scripts/generate_report.py \
        tests/unit/test_report_formatting.py
git commit -m "fix(report): re-serialise IGV fragments as script-safe JSON

Closes #216"
```

---

### Task A3: `innerHTML` → `textContent` for the IGV variant table

**Files:**
- Modify: `vntyper/templates/report_template.html:515`, `:546`
- Test: `tests/unit/test_template_escaping.py`

**Interfaces:** none.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_template_escaping.py`:

```python
def test_the_igv_variant_table_writes_text_not_markup() -> None:
    """#216: both cells carry igv-reports values derived from the sample's reads.

    ``innerHTML`` parses; ``textContent`` does not. Neither cell needs markup.
    """
    source = Path("vntyper/templates/report_template.html").read_text(encoding="utf-8")
    assert "cell.innerHTML" not in source
    assert "cell.textContent = headers[j];" in source
    assert "cell.textContent = rowData[j];" in source
```

- [ ] **Step 2: Run and watch it fail**

Run: `python -m pytest tests/unit/test_template_escaping.py::test_the_igv_variant_table_writes_text_not_markup -v`
Expected: FAIL — `cell.innerHTML` is present twice.

- [ ] **Step 3: Change both lines**

`report_template.html:515`: `cell.innerHTML = headers[j];` → `cell.textContent = headers[j];`
`report_template.html:546`: `cell.innerHTML = rowData[j];` → `cell.textContent = rowData[j];`

Leave `report_template.html:482` (`document.getElementById("igvContainer").innerHTML = "<p>…"`)
alone — that value is a VNtyper-authored constant with no interpolation.

- [ ] **Step 4: Run**

Run: `python -m pytest tests/unit/test_template_escaping.py tests/unit/test_generate_report.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add vntyper/templates/report_template.html tests/unit/test_template_escaping.py
git commit -m "fix(templates): write igv table cells as text, not markup

Refs #216"
```

---

# Workstream B — nothing writes into the input directory (#201, #162, #210)

**Commit scopes:** `utils`, `pipeline`, `bam`.

### Task B1: Pin the mechanism — a log path in a read-only directory raises

**Files:**
- Test: `tests/unit/test_input_tree_is_never_written.py` *(create)*

**Interfaces:**
- Produces: nothing. This task exists so the *reason* the bug is a crash is pinned before it is
  fixed, and so a later reviewer can see that `run_command` — not `samtools` — is what fails.

- [ ] **Step 1: Write the test**

```python
"""No VNtyper process may write into the input directory (#162, #201, #210).

SPECIFICATION: the input directory holds patient-derived data and is routinely
mounted read-only to protect it. Two production writes derived an output path
from an input path:

1. ``validate_bam_file`` wrote ``<input>.quickcheck.log`` beside the alignment
   (#201). On the web service that is the per-job directory on the shared upload
   volume, so ``run_vntyper_job``'s cleanup found a non-empty directory and its
   ``os.rmdir`` never ran -- one leaked directory and inode per job, forever, and
   a "still holds files" warning that fired on 100% of jobs and therefore could
   no longer report a genuinely unexpected leftover.
2. ``process_bam_to_fastq`` ran ``samtools index <in_bam>`` for non-fast BAM runs
   whenever ``<in_bam>.bai`` was missing -- including when a usable ``<stem>.bai``
   was already there (#210).

Both are crashes, not annoyances: ``run_command`` opens its log file before it
runs anything, so on a read-only mount write (1) is an unhandled
``PermissionError`` raised before quickcheck even executes.
"""

import os
import stat
from pathlib import Path

import pytest

from vntyper.scripts.utils import run_command, validate_bam_file

pytestmark = pytest.mark.unit


@pytest.fixture
def readonly_dir(tmp_path: Path) -> Path:
    """A directory that cannot be written to, restored on teardown.

    Skipped when the test runs as root, for whom the mode bits are advisory.
    """
    if os.geteuid() == 0:
        pytest.skip("root ignores the write bit, so this fixture cannot deny anything")
    d = tmp_path / "input"
    d.mkdir()
    (d / "sample.bam").write_bytes(b"BAM\x01")
    d.chmod(stat.S_IRUSR | stat.S_IXUSR)
    yield d
    d.chmod(stat.S_IRWXU)


def test_run_command_raises_when_its_log_file_is_not_writable(readonly_dir: Path) -> None:
    """The mechanism: the log is opened before the child process starts."""
    with pytest.raises(PermissionError):
        run_command("true", str(readonly_dir / "sample.bam.quickcheck.log"))
```

- [ ] **Step 2: Run — this one passes immediately**

Run: `python -m pytest tests/unit/test_input_tree_is_never_written.py -v`
Expected: PASS. It is a characterisation test of `run_command`, which is not changing. It is
here so the next test's failure is unambiguous.

- [ ] **Step 3: Commit**

```bash
git add tests/unit/test_input_tree_is_never_written.py
git commit -m "test(utils): pin that a log in a read-only dir raises before the child runs

Refs #162"
```

---

### Task B2: Give `validate_bam_file` a log destination

**Files:**
- Modify: `vntyper/scripts/utils.py:307-368`
- Modify: `vntyper/scripts/pipeline.py:218-221`
- Test: `tests/unit/test_input_tree_is_never_written.py`, `tests/unit/test_utils.py`

**Interfaces:**
- Produces: `validate_bam_file(file_path, cwd=None, log_dir=None)`. `log_dir=None` keeps today's
  behaviour — the log lands beside the input. This default is the issue author's explicit
  recommendation ("the parameter should default to today's behaviour and be passed explicitly by
  the pipeline; that keeps this a contained change rather than a breaking one") and is followed,
  not re-litigated. `pipeline.py` passes `log_dir=output_dir` at both call sites.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_input_tree_is_never_written.py`:

```python
def test_validate_bam_file_writes_its_log_under_log_dir(readonly_dir: Path, tmp_path: Path) -> None:
    """The whole point: a read-only input mount must survive a run."""
    out = tmp_path / "out"
    out.mkdir()
    validate_bam_file(str(readonly_dir / "sample.bam"), log_dir=str(out))
    assert (out / "sample.bam.quickcheck.log").exists()
    assert list(readonly_dir.iterdir()) == [readonly_dir / "sample.bam"]


def test_validate_bam_file_still_logs_beside_the_input_without_log_dir(tmp_path: Path) -> None:
    """The default is unchanged, deliberately -- this stays a contained change."""
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM\x01")
    validate_bam_file(str(bam))
    assert (tmp_path / "sample.bam.quickcheck.log").exists()
```

Both call the real `samtools quickcheck` through `run_command` with `critical=False`. On a machine
without samtools the command fails and `validate_bam_file` raises `ValueError` — the log file is
still written, which is what is asserted, so wrap the call:

```python
    with contextlib.suppress(ValueError):
        validate_bam_file(str(readonly_dir / "sample.bam"), log_dir=str(out))
```

Use `contextlib.suppress(ValueError)` in **both** tests and add `import contextlib`. The
assertion is about where the log landed, not about whether the alignment is valid — a 4-byte
stub is not a real BAM and must not be expected to pass quickcheck.

Append to `tests/unit/test_utils.py` (it already has `pytestmark`):

```python
def test_validate_bam_file_passes_the_log_dir_through_to_run_command(tmp_path):
    """#201: the log is an artefact of the run and belongs with the run's output."""
    bam_file = tmp_path / "in" / "sample.bam"
    bam_file.parent.mkdir()
    bam_file.write_text("dummy")
    out = tmp_path / "out"
    out.mkdir()
    with mock.patch("vntyper.scripts.utils.run_command", return_value=True) as run:
        validate_bam_file(str(bam_file), log_dir=str(out))
    log_file = run.call_args[0][1]
    assert log_file == str(out / "sample.bam.quickcheck.log")
    assert str(bam_file.parent) not in log_file
```

- [ ] **Step 2: Run and watch it fail**

Run: `python -m pytest tests/unit/test_input_tree_is_never_written.py tests/unit/test_utils.py -v`
Expected: FAIL — `TypeError: validate_bam_file() got an unexpected keyword argument 'log_dir'`.

- [ ] **Step 3: Implement**

`utils.py` — change the signature to `def validate_bam_file(file_path, cwd=None, log_dir=None):`,
extend the docstring's `Args:` block with:

```
        log_dir (str, optional): Directory to write the ``samtools quickcheck`` log
            into. Defaults to None, which writes it beside the input alignment --
            today's behaviour, kept so this stays a contained change. The pipeline
            passes its output directory, because the log is an artefact of the run
            and the input directory is routinely mounted read-only (#162, #201).
            ``run_command`` opens the log before it runs anything, so deriving this
            path from the input made a read-only mount raise ``PermissionError``
            before quickcheck ever executed.
```

and replace line 362 with:

```python
    log_file = f"{file_path}.quickcheck.log" if log_dir is None else str(Path(log_dir) / f"{Path(file_path).name}.quickcheck.log")
```

That line exceeds 120 characters — split it:

```python
    if log_dir is None:
        log_file = f"{file_path}.quickcheck.log"
    else:
        log_file = str(Path(log_dir) / f"{Path(file_path).name}.quickcheck.log")
```

Confirm `from pathlib import Path` is imported in `utils.py`; add it if not.

`pipeline.py:218-221` — pass the run's output directory:

```python
        if input_type == "BAM":
            validate_bam_file(bam, cwd=project_root, log_dir=output_dir)
        elif input_type == "CRAM":
            validate_bam_file(cram, cwd=project_root, log_dir=output_dir)
```

Check that `output_dir` is in scope and already created at that point; if the directory tree is
built later, move the `validate_bam_file` calls after `dirs` is created rather than creating the
directory early — and say so in the commit message.

- [ ] **Step 4: Run**

Run:
```bash
python -m pytest tests/unit/test_input_tree_is_never_written.py tests/unit/test_utils.py \
  tests/unit/test_shell_quoting.py tests/unit/test_pipeline_guards.py tests/unit/test_pipeline_cwd.py \
  tests/unit/web/test_cram_alignment_handoff.py -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/utils.py vntyper/scripts/pipeline.py \
        tests/unit/test_input_tree_is_never_written.py tests/unit/test_utils.py
git commit -m "fix(utils): write the quickcheck log to the output dir, not beside the input

Closes #201"
```

---

### Task B3: Resolve the BAM index htslib's way, and build it into the output directory

**Files:**
- Modify: `vntyper/scripts/command_builders.py:215-226`
- Modify: `vntyper/scripts/fastq_bam_processing.py:160-177`
- Test: `tests/unit/test_command_builders.py`, `tests/unit/test_fastq_bam_command_wiring.py`,
  `tests/unit/test_input_tree_is_never_written.py`

**Interfaces:**
- Consumes: nothing from B1/B2.
- Produces: `build_samtools_index_command(*, samtools_path: str, bam_file: str | Path,
  output_bai: str | Path | None = None) -> str`. With `output_bai` set it emits
  `samtools index -o <out> <bam>`; without it, the existing `samtools index <bam>`.
  `-o` requires samtools ≥ 1.15; `conda/environment_vntyper.yml` pins 1.20.
- Also produces: `resolve_bam_index(in_bam: str | Path) -> str | None` in
  `fastq_bam_processing.py` — htslib's own order, `<file>.bai` then `<stem>.bai`, else None.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_command_builders.py`:

```python
def test_the_index_command_takes_an_output_path() -> None:
    """#210/#162: the index must be buildable somewhere other than beside the input."""
    command = build_samtools_index_command(
        samtools_path="samtools", bam_file="/in/sample.bam", output_bai="/out/sample.bam.bai"
    )
    assert command == "samtools index -o /out/sample.bam.bai /in/sample.bam"


def test_the_index_command_without_an_output_path_is_unchanged() -> None:
    command = build_samtools_index_command(samtools_path="samtools", bam_file="/in/sample.bam")
    assert command == "samtools index /in/sample.bam"


def test_the_index_output_path_is_quoted() -> None:
    command = build_samtools_index_command(
        samtools_path="samtools", bam_file="/in/sample.bam", output_bai="/out/has space.bai"
    )
    assert shlex.split(command) == ["samtools", "index", "-o", "/out/has space.bai", "/in/sample.bam"]
```

Append to `tests/unit/test_input_tree_is_never_written.py`:

```python
def test_an_alternate_index_name_is_found_and_no_second_index_is_built(tmp_path: Path) -> None:
    """#210: the upload path accepts ``sample.bai``; htslib resolves it, so must we."""
    from vntyper.scripts.fastq_bam_processing import resolve_bam_index

    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM\x01")
    alternate = tmp_path / "sample.bai"
    alternate.write_bytes(b"BAI\x01")
    assert resolve_bam_index(bam) == str(alternate)


def test_the_conventional_index_name_wins_over_the_alternate(tmp_path: Path) -> None:
    """htslib tries ``<file>.bai`` first."""
    from vntyper.scripts.fastq_bam_processing import resolve_bam_index

    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM\x01")
    (tmp_path / "sample.bam.bai").write_bytes(b"BAI\x01")
    (tmp_path / "sample.bai").write_bytes(b"BAI\x01")
    assert resolve_bam_index(bam) == str(tmp_path / "sample.bam.bai")


def test_no_index_at_all_resolves_to_none(tmp_path: Path) -> None:
    from vntyper.scripts.fastq_bam_processing import resolve_bam_index

    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM\x01")
    assert resolve_bam_index(bam) is None
```

Append to `tests/unit/test_fastq_bam_command_wiring.py` — a wiring test that the index the
non-fast BAM path builds lands under `output`, mocking `run_command` and
`extract_unmapped_reads_from_offset`. Follow whatever mocking shape that file already uses; the
assertion is:

```python
    index_commands = [c for c in issued_commands if " index " in c]
    assert index_commands, "the non-fast BAM path must index when no index exists"
    assert all(str(tmp_path / "input") not in c for c in index_commands), (
        "no samtools index command may name a path inside the input directory as its output"
    )
```

- [ ] **Step 2: Run and watch it fail**

Run:
```bash
python -m pytest tests/unit/test_command_builders.py tests/unit/test_input_tree_is_never_written.py \
  tests/unit/test_fastq_bam_command_wiring.py -v
```
Expected: FAIL — `build_samtools_index_command() got an unexpected keyword argument 'output_bai'`
and `ImportError: cannot import name 'resolve_bam_index'`.

- [ ] **Step 3: Implement**

`command_builders.py`:

```python
def build_samtools_index_command(
    *, samtools_path: str, bam_file: str | Path, output_bai: str | Path | None = None
) -> str:
    """
    Build the ``samtools index`` command for a BAM file.

    Args:
        samtools_path (str): samtools invocation from config.
        bam_file (str | Path): The BAM to index.
        output_bai (str | Path, optional): Where to write the index. Defaults to
            None, which lets samtools write ``<bam_file>.bai`` beside the input --
            correct for a BAM the pipeline itself produced in its output
            directory, wrong for the user's input alignment, whose directory is
            routinely read-only (#162, #210). ``-o`` requires samtools >= 1.15;
            ``conda/environment_vntyper.yml`` pins 1.20.

    Returns:
        str: The complete index command.
    """
    if output_bai is None:
        return f"{samtools_path} index {quote_path(bam_file)}"
    return f"{samtools_path} index -o {quote_path(output_bai)} {quote_path(bam_file)}"
```

`fastq_bam_processing.py` — add the resolver near the top of the module:

```python
def resolve_bam_index(in_bam: str | Path) -> str | None:
    """Find an existing BAM index the way htslib itself does.

    htslib tries ``<file>.bai`` and then ``<stem>.bai``. The pipeline used to
    reconstruct only the first, so given the ``sample.bai`` that the upload
    endpoint and the worker both deliberately accept it found nothing and built a
    second index beside the input that nothing else knew about (#210).

    Args:
        in_bam: The alignment whose index is wanted.

    Returns:
        str | None: The existing index, or None when there is none.
    """
    bam = Path(in_bam)
    for candidate in (Path(f"{bam}.bai"), bam.with_suffix(".bai")):
        if candidate.exists():
            return str(candidate)
    return None
```

and replace `fastq_bam_processing.py:160-177` with:

```python
        if file_format.lower() == "bam":
            # Use the offset-based extraction. The index is resolved htslib's way
            # and, when one has to be built, it is built into the *output*
            # directory: the input directory holds patient data and is routinely
            # mounted read-only (#162, #210).
            bam_bai = resolve_bam_index(in_bam)
            if bam_bai is None:
                bam_bai = str(Path(output) / f"{output_name}_input.bam.bai")
                index_cmd = build_samtools_index_command(
                    samtools_path=samtools_path, bam_file=in_bam, output_bai=bam_bai
                )
                log_file_index = Path(output) / f"{output_name}_unmapped_index.log"
                logger.info(f"Indexing BAM before extracting unmapped: {index_cmd}")
                success = run_command(str(index_cmd), str(log_file_index), critical=True)
                if not success:
                    raise RuntimeError("Indexing BAM file failed.")

            logger.info("Extracting unmapped reads using offset calculation...")
            extract_unmapped_reads_from_offset(
                bam_file=str(in_bam),
                bai_file=bam_bai,
                output_bam=str(unmapped_bam),
            )
```

Leave the *re*-index at `fastq_bam_processing.py:221-226` alone: `final_bam` is already inside
the output directory, so `samtools index <final_bam>` writes there.

- [ ] **Step 4: Run**

Run:
```bash
python -m pytest tests/unit/test_command_builders.py tests/unit/test_fastq_bam_command_wiring.py \
  tests/unit/test_input_tree_is_never_written.py tests/unit/test_extract_unmapped.py \
  tests/unit/test_file_processing.py tests/unit/web/test_tasks.py -v
```
Expected: PASS. `tests/unit/web/test_tasks.py` has two index tests (lines 603, 629) covering the
worker's own preflight index — they exercise `docker/app/tasks.py`, not this file, and must be
unaffected. If either fails, **stop and report**.

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/command_builders.py vntyper/scripts/fastq_bam_processing.py \
        tests/unit/test_command_builders.py tests/unit/test_fastq_bam_command_wiring.py \
        tests/unit/test_input_tree_is_never_written.py
git commit -m "fix(bam): resolve the index htslib's way and build it in the output dir

Closes #210
Closes #162"
```

---

### Task B4: Confirm the worker's `os.rmdir` is now reachable

**Files:**
- Test: `tests/unit/web/test_tasks.py` *(read; edit only if a test needs it)*

**Interfaces:** none. `docker/app/tasks.py` is deliberately **not** edited — its `os.rmdir` at
line 409 simply becomes reachable now that nothing else lands in the directory.

- [ ] **Step 1: Read the two existing tests**

`test_a_successful_job_is_marked_completed_and_cleans_up_its_inputs` (line 454) and
`test_cleanup_leaves_a_nonempty_input_directory_and_warns_instead_of_deleting_it` (line 683).

- [ ] **Step 2: Run them**

Run: `python -m pytest tests/unit/web/test_tasks.py -v`
Expected: PASS, unchanged. Both mock the pipeline subprocess, so neither ever produced a
quickcheck log — they were always testing the *intended* behaviour, which is now also the actual
behaviour.

- [ ] **Step 3: Add the missing half if it is absent**

If no test asserts that the "still holds files" warning does **not** fire on a clean job, add one
beside `test_a_successful_job_is_marked_completed_and_cleans_up_its_inputs`:

```python
def test_a_clean_job_does_not_fire_the_leftover_warning(...):
    """#201: the guard must keep working -- it must only stop firing on the normal path."""
    ...
    assert not any("still holds files" in r.getMessage() for r in caplog.records)
```

- [ ] **Step 4: Commit (only if step 3 added something)**

```bash
git add tests/unit/web/test_tasks.py
git commit -m "test(tasks): assert the leftover warning is silent on a clean job

Refs #201"
```

---

# Workstream C — cohort identity and order (#206, #205)

**Commit scopes:** `cohort`, `config`.

### Task C1: Widen the pseudonym digest, config-driven

**Files:**
- Modify: `vntyper/config.json` (add `cohort.pseudonym`)
- Modify: `vntyper/scripts/cohort_inputs.py:160-177`
- Test: `tests/unit/test_cohort_identity.py` *(create)*, `tests/unit/test_cohort_inputs.py`

**Interfaces:**
- Produces: `pseudonymized_sample_name(prefix: Any, original_sample: str, *, algorithm: str = DEFAULT_PSEUDONYM_ALGORITHM, length: int = DEFAULT_PSEUDONYM_LENGTH) -> str`
  in `cohort_inputs.py`, plus module constants `DEFAULT_PSEUDONYM_ALGORITHM = "sha256"` and
  `DEFAULT_PSEUDONYM_LENGTH = 12`. `cohort_summary.aggregate_cohort` reads the pair out of
  `config` and passes them.
- Config shape (`vntyper/config.json`, extending the existing `"cohort"` object):

```json
  "cohort": {
    "kestrel_result_file": "kestrel_result.tsv",
    "advntr_result_file": "output_adVNTR.tsv",
    "pseudonym": {
      "algorithm": "sha256",
      "digest_characters": 12
    }
  }
```

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_cohort_identity.py`:

```python
"""A cohort pseudonym must identify exactly one sample (#206).

SPECIFICATION: the pseudonym was ``<prefix>`` plus the first **five** hex
characters of the MD5 of the sample name -- 20 bits. ``sample_mapping`` is keyed
on the pseudonym, so a collision silently overwrote the earlier original name:
two patients' rows became indistinguishable across every export,
``sample_categories()`` counted them as one sample, and
``pseudonymization_table.tsv`` mapped the shared pseudonym to whichever original
was written last. Birthday probability of at least one collision is
``1 - exp(-n(n-1)/2**21)``: ~37.9% at 1,000 samples.

The verified first collision in ``sample_0..sample_19999`` was ``sample_42`` and
``sample_919``, both MD5-prefixing to ``168eb``. That exact pair is the probe
below.
"""

import json
from pathlib import Path

import pytest

from vntyper.scripts.cohort_inputs import (
    DEFAULT_PSEUDONYM_ALGORITHM,
    DEFAULT_PSEUDONYM_LENGTH,
    pseudonymized_sample_name,
)

pytestmark = pytest.mark.unit


def test_the_verified_md5_collision_no_longer_collides() -> None:
    assert pseudonymized_sample_name("anon_", "sample_42") != pseudonymized_sample_name("anon_", "sample_919")


def test_pseudonyms_are_injective_over_twenty_thousand_names() -> None:
    """The exact probe that found the bug, at the configured width."""
    names = [f"sample_{i}" for i in range(20000)]
    pseudonyms = {pseudonymized_sample_name("anon_", n) for n in names}
    assert len(pseudonyms) == len(names)


def test_the_pseudonym_is_stable_across_calls() -> None:
    """The pseudonymization table written beside the report has to stay meaningful."""
    assert pseudonymized_sample_name("anon_", "s1") == pseudonymized_sample_name("anon_", "s1")


def test_the_default_digest_is_twelve_characters_of_sha256() -> None:
    assert DEFAULT_PSEUDONYM_ALGORITHM == "sha256"
    assert DEFAULT_PSEUDONYM_LENGTH == 12
    assert pseudonymized_sample_name("anon_", "s1") == "anon_" + __import__("hashlib").sha256(b"s1").hexdigest()[:12]


def test_the_digest_width_is_honoured() -> None:
    assert len(pseudonymized_sample_name("", "s1", length=5)) == 5
    assert len(pseudonymized_sample_name("", "s1", length=32)) == 32


def test_a_non_string_prefix_is_interpolated_rather_than_concatenated() -> None:
    """Preserved behaviour: ``--pseudonymize-samples`` may carry a non-string."""
    assert pseudonymized_sample_name(True, "s1").startswith("True")


def test_an_unavailable_algorithm_is_refused_by_name() -> None:
    with pytest.raises(ValueError, match="not-a-hash"):
        pseudonymized_sample_name("anon_", "s1", algorithm="not-a-hash")


def test_the_shipped_config_declares_the_pseudonym_settings() -> None:
    """Config-driven, never hardcoded: the values live in config.json."""
    config = json.loads(Path("vntyper/config.json").read_text())
    assert config["cohort"]["pseudonym"] == {"algorithm": "sha256", "digest_characters": 12}
```

- [ ] **Step 2: Run and watch it fail**

Run: `python -m pytest tests/unit/test_cohort_identity.py -v`
Expected: FAIL at collection — `DEFAULT_PSEUDONYM_ALGORITHM` does not exist.

- [ ] **Step 3: Implement**

In `cohort_inputs.py`, after the `PIPELINE_SUMMARY_FILENAME` constant:

```python
#: Digest used to build a cohort pseudonym. Overridable through
#: ``config["cohort"]["pseudonym"]``; declared here so a config that omits the key
#: -- which ``--config-path`` produces, because it replaces rather than merges
#: (AGENTS.md trap 2) -- does not raise.
DEFAULT_PSEUDONYM_ALGORITHM = "sha256"

#: Hex characters of the digest a pseudonym carries. Twelve is 48 bits: at 10,000
#: samples the birthday probability of a collision is ~1.8e-10. The previous value
#: was five (20 bits), which collides with probability ~37.9% at 1,000 samples
#: (#206).
DEFAULT_PSEUDONYM_LENGTH = 12
```

and replace `pseudonymized_sample_name`:

```python
def pseudonymized_sample_name(
    prefix: Any,
    original_sample: str,
    *,
    algorithm: str = DEFAULT_PSEUDONYM_ALGORITHM,
    length: int = DEFAULT_PSEUDONYM_LENGTH,
) -> str:
    """Build the pseudonym a sample is reported under.

    The pseudonym is the caller's prefix followed by the first ``length`` hex
    digits of the digest of the original name, so it is stable across runs and the
    pseudonymization table written beside the report stays meaningful.

    MD5 at five characters was the original scheme and is gone for two reasons:
    20 bits collides at realistic cohort sizes (#206), and ``hashlib.md5()``
    raises on a FIPS-enabled build unless it is called with
    ``usedforsecurity=False``.

    Args:
        prefix: The value ``--pseudonymize-samples`` supplied. Interpolated rather
            than concatenated, so a non-string does not raise.
        original_sample: The sample's identity.
        algorithm: A ``hashlib`` algorithm name.
        length: How many hex characters of the digest to keep.

    Returns:
        str: The pseudonym.

    Raises:
        ValueError: If ``algorithm`` is not available in this interpreter's
            ``hashlib``. Refused by name rather than silently falling back,
            because a silent fallback changes every pseudonym in the report
            without saying so.
    """
    if algorithm not in hashlib.algorithms_available:
        msg = f"Unknown pseudonym digest algorithm: {algorithm}"
        logger.error(msg)
        raise ValueError(msg)
    digest = hashlib.new(algorithm, original_sample.encode()).hexdigest()[:length]
    return f"{prefix}{digest}"
```

In `cohort_summary.py:394-396`, read the config through `.get()` chains:

```python
    pseudonym_config = config.get("cohort", {}).get("pseudonym", {})
    pseudonym_algorithm = pseudonym_config.get("algorithm", DEFAULT_PSEUDONYM_ALGORITHM)
    pseudonym_length = pseudonym_config.get("digest_characters", DEFAULT_PSEUDONYM_LENGTH)
```

hoisted above the sample loop, and pass both into `pseudonymized_sample_name`.

Add the config block to `vntyper/config.json` exactly as shown in the Interfaces section.

- [ ] **Step 4: Run**

Run:
```bash
python -m pytest tests/unit/test_cohort_identity.py tests/unit/test_cohort_inputs.py \
  tests/unit/test_cohort_exports.py tests/unit/test_cohort_summary_oracle.py -v
```
Expected: PASS. `tests/unit/test_cohort_inputs.py` has a test pinning the *old* collision
(`test_different_sample_names_get_different_pseudonyms` and the collision test #199 added,
referencing #206) — rewrite it to pin the new behaviour and cite this milestone. Do not delete it.

- [ ] **Step 5: Commit**

```bash
git add vntyper/config.json vntyper/scripts/cohort_inputs.py vntyper/scripts/cohort_summary.py \
        tests/unit/test_cohort_identity.py tests/unit/test_cohort_inputs.py
git commit -m "fix(cohort): widen the pseudonym digest to 12 sha256 chars, config-driven

Closes #206"
```

---

### Task C2: Refuse to merge two samples onto one pseudonym

**Files:**
- Modify: `vntyper/scripts/cohort_summary.py:387-400`
- Test: `tests/unit/test_cohort_identity.py`

**Interfaces:**
- Consumes: `pseudonymized_sample_name` from Task C1.
- Produces: nothing importable. `aggregate_cohort` raises `ValueError` on a collision.

- [ ] **Step 1: Write the failing test**

```python
def test_two_originals_mapping_to_one_pseudonym_abort_the_cohort(tmp_path, monkeypatch) -> None:
    """A cohort that silently merges two patients is worse than one that refuses to run.

    At 48 bits this is unreachable in practice, so it is forced here by pinning the
    digest to a width narrow enough to collide. That is the point: no digest width
    makes a silent merge acceptable, so the guard does not depend on the width.
    """
    from vntyper.scripts import cohort_summary

    monkeypatch.setattr(cohort_summary, "pseudonymized_sample_name", lambda prefix, name, **kw: f"{prefix}same")
    ...  # two sample directories, each with a pipeline_summary.json
    with pytest.raises(ValueError) as excinfo:
        cohort_summary.aggregate_cohort(
            input_paths=[str(a), str(b)],
            output_dir=str(tmp_path / "out"),
            summary_file="cohort.html",
            config=config,
            pseudonymize_samples="anon_",
        )
    message = str(excinfo.value)
    assert "sample_a" in message and "sample_b" in message


def test_the_same_sample_reached_twice_is_not_a_collision(tmp_path) -> None:
    """De-duplication is not a collision: the same original mapping to itself is fine."""
```

Build the two sample directories with a minimal valid `pipeline_summary.json`
(`{"version": "x", "input_files": {}, "steps": []}`) — the cohort drops samples it cannot read
but must still have built the mapping first, so the guard fires either way.

- [ ] **Step 2: Run and watch it fail**

Run: `python -m pytest tests/unit/test_cohort_identity.py -k collision -v`
Expected: FAIL — `DID NOT RAISE ValueError`; the second sample silently overwrites the first.

- [ ] **Step 3: Implement**

In `cohort_summary.py`, inside the sample loop, replace the bare
`sample_mapping[pseudonym] = original_sample` with:

```python
        if pseudonymize_samples:
            pseudonym = pseudonymized_sample_name(
                pseudonymize_samples, original_sample, algorithm=pseudonym_algorithm, length=pseudonym_length
            )
            existing = sample_mapping.get(pseudonym)
            if existing is not None and existing != original_sample:
                # #206: sample_mapping is keyed on the pseudonym, so this used to
                # overwrite silently -- two patients' rows became one row, and
                # sample_categories() counted one result where there were two.
                msg = (
                    f"Pseudonym collision: {original_sample!r} and {existing!r} both map to {pseudonym!r}. "
                    "Widen cohort.pseudonym.digest_characters in the configuration and re-run."
                )
                logger.error(msg)
                raise ValueError(msg)
            sample_mapping[pseudonym] = original_sample
```

- [ ] **Step 4: Run**

Run: `python -m pytest tests/unit/test_cohort_identity.py tests/unit/test_cohort_summary_oracle.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/cohort_summary.py tests/unit/test_cohort_identity.py
git commit -m "fix(cohort): refuse to map two samples onto one pseudonym

Refs #206"
```

---

### Task C3: Give ZIP samples an identity and an order that contain no temp path

**Files:**
- Modify: `vntyper/scripts/cohort_inputs.py:56-140`
- Modify: `vntyper/scripts/cohort_summary.py:380-400`
- Test: `tests/unit/test_cohort_inputs.py`, `tests/unit/test_cohort_identity.py`

**Interfaces:**
- Produces: `discover_sample_directories(input_paths: list[str]) -> tuple[list[DiscoveredSample], list[str]]`,
  where

```python
@dataclass(frozen=True)
class DiscoveredSample:
    """One sample the cohort found, with an identity that does not move between runs."""

    directory: Path
    #: The name the sample is reported under, before any pseudonymisation.
    identity: str
    #: What the ordering is total on: (index of the input path on the command
    #: line, the sample's path relative to that input's root). Contains no
    #: temporary component, so two runs agree.
    sort_key: tuple[int, tuple[str, ...]]
```

  `cohort_summary.aggregate_cohort` consumes `.directory` and `.identity`, replacing
  `Path(sample_dir).name` at `cohort_summary.py:393`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_cohort_identity.py`:

```python
def _summary(input_files: dict) -> str:
    return json.dumps({"version": "2.0.0", "input_files": input_files, "steps": []})


def test_a_root_level_zip_sample_is_identified_by_the_run_s_own_input(tmp_path) -> None:
    """#205: this layout is what the web worker produces, so it is the normal path.

    ``processed_dirs.add(temp_path)`` added the ``tempfile.mkdtemp`` root itself,
    and ``cohort_summary`` took ``Path(sample_dir).name`` as the sample -- so the
    reported ``Sample`` was a ``cohort_zip_*`` string that differed on every run,
    and any pseudonym derived from it differed too.
    """
    src = tmp_path / "run"
    src.mkdir()
    (src / "pipeline_summary.json").write_text(_summary({"bam": "patient1.bam"}))
    archive = shutil.make_archive(str(tmp_path / "job7"), "zip", root_dir=src)

    samples, temp_dirs = discover_sample_directories([archive])
    try:
        assert [s.identity for s in samples] == ["patient1"]
        assert "cohort_zip_" not in samples[0].identity
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_root_level_zip_without_input_files_falls_back_to_the_archive_stem(tmp_path) -> None:
    src = tmp_path / "run"
    src.mkdir()
    (src / "pipeline_summary.json").write_text(_summary({}))
    archive = shutil.make_archive(str(tmp_path / "job7"), "zip", root_dir=src)

    samples, temp_dirs = discover_sample_directories([archive])
    try:
        assert [s.identity for s in samples] == ["job7"]
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_a_cram_run_is_identified_by_its_cram(tmp_path) -> None:
    ...  # input_files {"cram": "patient2.cram"} -> "patient2"


def test_a_fastq_run_is_identified_by_its_first_mate(tmp_path) -> None:
    ...  # input_files {"fastq1": "patient3_R1.fastq.gz", "fastq2": ...} -> "patient3_R1"


def test_a_sample_in_a_subdirectory_keeps_its_directory_name(tmp_path) -> None:
    """Only the root-level case changes: elsewhere the directory name is meaningful."""


def test_the_order_is_the_same_across_two_processes(tmp_path) -> None:
    """#205: the sort key used to contain the randomised mkdtemp component.

    Two processes, because ``Path.__hash__`` is the hash of the path string and
    Python randomises string hashing per process -- which is the mechanism the
    existing cross-process test in ``test_cohort_inputs.py`` was written for.
    """
    # Build three archives, run discover_sample_directories in two subprocesses
    # with different PYTHONHASHSEED, compare the identity lists.


def test_the_sort_key_contains_no_temporary_path_component(tmp_path) -> None:
    samples, temp_dirs = discover_sample_directories([archive_a, archive_b])
    try:
        assert all("cohort_zip_" not in part for s in samples for part in s.sort_key[1])
    finally:
        cleanup_temp_dirs(temp_dirs)


def test_two_zips_come_back_in_command_line_order(tmp_path) -> None:
    """The first element of the sort key is the input's position on the command line."""
```

Fill every `...` with real bodies following the shape of the two complete tests above. Every test
must call `cleanup_temp_dirs(temp_dirs)` in a `finally`, or `tmp_path` fills with extracted
archives.

- [ ] **Step 2: Run and watch it fail**

Run: `python -m pytest tests/unit/test_cohort_identity.py -v`
Expected: FAIL — `discover_sample_directories` returns `list[Path]`, so `.identity` raises
`AttributeError`.

- [ ] **Step 3: Implement**

In `cohort_inputs.py`:

1. Add the `DiscoveredSample` dataclass exactly as in the Interfaces block, with a module-level
   docstring paragraph explaining why the sort key exists.
2. Add a private helper:

```python
def _identity_from_summary(sample_dir: Path, fallback: str) -> str:
    """Name a sample from the input file its own run recorded.

    The directory name is the right identity almost everywhere, but not for a zip
    whose ``pipeline_summary.json`` sits at the root: that directory is the
    randomised ``tempfile.mkdtemp`` root, so the reported sample differed on every
    run (#205). That layout is what the web worker produces, so it is the normal
    path for web cohorts rather than an edge case.

    Args:
        sample_dir: The directory holding ``pipeline_summary.json``.
        fallback: The identity to use when the summary records no input files --
            the archive's filename stem.

    Returns:
        str: The sample's identity.
    """
```

   Read `input_files` and return the stem of `bam`, else `cram`, else `fastq1`. `Path(name).stem`
   strips one suffix, so `patient3_R1.fastq.gz` yields `patient3_R1.fastq` — call `.stem` twice
   for `.gz`, or strip known double suffixes. Pin whichever you choose with the FASTQ test above.
   Any failure to read the summary returns `fallback`; this helper never raises, matching the
   module's stated "one bad sample must not abort the cohort" contract.

3. Rework the loop to build `DiscoveredSample` values. Keep the de-duplicating set — key it on
   `directory` — and sort on `sort_key` on the way out. For a directory input, `sort_key` is
   `(index, sample_dir.relative_to(path).parts)`; for a ZIP, `(index, sample_dir.relative_to(temp_path).parts)`.
   A sample at the root gets `()` as its relative parts, which sorts first — that is fine and
   deterministic.

In `cohort_summary.py:392-393`:

```python
    for sample in processed_dirs:
        sample_dir = sample.directory
        original_sample = sample.identity
```

Update the `if not processed_dirs:` guard and the `logger.info(f"Processing sample directory:
{sample_dir} as {pseudonym}")` line to match.

- [ ] **Step 4: Run**

Run:
```bash
python -m pytest tests/unit/test_cohort_identity.py tests/unit/test_cohort_inputs.py \
  tests/unit/test_cohort_exports.py tests/unit/test_cohort_categories.py \
  tests/unit/test_cohort_tables.py tests/unit/test_cohort_summary_oracle.py \
  tests/unit/test_cohort_summary_escaping.py -v
```
Expected: PASS. The characterisation test #199 added to pin the old ZIP behaviour **will** fail —
rewrite it to pin the new behaviour, citing #205. Do not delete it.

Then run `make type-check` — the return type of `discover_sample_directories` changed and
`cohort_inputs.py` is fully annotated.

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/cohort_inputs.py vntyper/scripts/cohort_summary.py \
        tests/unit/test_cohort_identity.py tests/unit/test_cohort_inputs.py
git commit -m "fix(cohort): identify and order zip samples without the temp path

Closes #205"
```

---

### Task C4: Update the golden-cohort ordering justification and the oracle's pseudonym

**Files:**
- Modify: `scripts/golden_cohort/compare.py:70-108` (`COHORT_ORDER_WHY` and the note above it)
- Modify: `tests/unit/test_golden_cohort_compare.py:213-232`
- Modify: `tests/unit/test_cohort_summary_oracle.py:862-863`
- Modify: `tests/unit/test_cohort_inputs.py:182-186` (`_ORDER_PROBE`)

**Interfaces:** none. This task makes three existing tripwires tell the truth again.

**Why this is a task and not a cleanup.** `COHORT_ORDER_WHY` is the recorded justification for
the golden-cohort gate **not comparing** cohort sample order. Reason 2 of that justification is
"ZIP inputs, on any version … whose random suffix is part of the path and therefore part of the
sort key". Task C3 removes exactly that, so leaving the note in place would mean the gate keeps
normalising away a difference it no longer has a reason to normalise — a silent loss of coverage
dressed up as a documented decision.

- [ ] **Step 1: Run the three tripwires and watch them fail**

Run:
```bash
python -m pytest tests/unit/test_golden_cohort_compare.py tests/unit/test_cohort_summary_oracle.py \
  tests/unit/test_cohort_inputs.py -v
```
Expected: FAIL —
`test_the_cohort_order_justification_does_not_claim_discovery_is_unordered_today` (the return
statement it greps for is gone), `test_the_cohort_order_justification_keeps_the_zip_reason` (the
reason it greps for is now false), `test_pseudonymised_samples_reach_the_report_under_their_
pseudonyms` (`anon_65622` was the MD5 pseudonym), and the `_ORDER_PROBE` subprocess
(`AttributeError: 'DiscoveredSample' object has no attribute 'name'`).

- [ ] **Step 2: Fix `_ORDER_PROBE`**

`tests/unit/test_cohort_inputs.py:182-186`:

```python
    "print(json.dumps([s.identity for s in discover_sample_directories(sys.argv[1:])[0]]))"
```

`identity` rather than `directory.name`: identity is what the report actually shows, so the
cross-process test now pins the value that matters instead of a path component.

- [ ] **Step 3: Rewrite `COHORT_ORDER_WHY` and its note**

Reason 1 (cross-version) survives verbatim — a baseline predating the determinism fix still
iterates the set directly. Reason 2 must change from "ZIP order is irreproducible on any version"
to a version-boundary statement: ZIP order was irreproducible **before this change** because the
`mkdtemp` suffix was in the sort key, and is reproducible after it because the sort key is now
`(input index, path relative to that input's root)`. Keep the "what this costs" paragraph — the
gate still does not attest ordering — and update the cited test names to whichever tests exist in
`test_cohort_inputs.py` after C3.

Do **not** delete `tempfile.mkdtemp(prefix="cohort_zip_")` from `cohort_inputs.py`: extraction
still goes to a temp directory. Only its role in the sort key is gone.

- [ ] **Step 4: Update the two golden-cohort tripwires to assert the new truth**

`test_the_cohort_order_justification_does_not_claim_discovery_is_unordered_today` should grep for
whatever C3's return statement actually is, and its docstring should say why it changed.
`test_the_cohort_order_justification_keeps_the_zip_reason` becomes a test that the justification
states the ZIP reason is **version-bounded**, citing #205.

- [ ] **Step 5: Update the oracle's pseudonym assertions**

`tests/unit/test_cohort_summary_oracle.py:862-863` — replace the literal `anon_65622` with the
value computed from the new digest. Compute it, do not guess:

```bash
python -c "import hashlib; print('anon_' + hashlib.sha256(b'sample_one').hexdigest()[:12])"
```

Assert the computed literal, not a call to `pseudonymized_sample_name` — a test that recomputes
the value with the code under test asserts nothing.

Then check `EXPECTED_FINGERPRINT` (line 174). If the fingerprint test renders a pseudonymised
cohort, it moves and must be re-recorded with the reason in the commit message. If it does not,
it must **not** move — and if it does anyway, stop and report.

- [ ] **Step 6: Run**

Run:
```bash
python -m pytest tests/unit/test_golden_cohort_compare.py tests/unit/test_cohort_summary_oracle.py \
  tests/unit/test_cohort_inputs.py tests/unit/test_cohort_identity.py \
  tests/unit/test_golden_cohort_normalise.py tests/unit/test_golden_cohort_matrix.py -v
```
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add scripts/golden_cohort/compare.py tests/unit/test_golden_cohort_compare.py \
        tests/unit/test_cohort_summary_oracle.py tests/unit/test_cohort_inputs.py
git commit -m "docs(gate): the zip ordering reason is now version-bounded

Refs #205
Refs #206"
```

---

# Workstream D — `vntyper report` resolves its VCF (#167)

**Commit scope:** `cli`.

### Task D1: Add `--vcf-file` and resolve it from the run directory

**Files:**
- Modify: `vntyper/scripts/cli_parser.py:201-213`
- Modify: `vntyper/scripts/cli_report.py:26-32` (imports), `:100-139`
- Test: `tests/unit/test_cli_report.py`, `tests/unit/test_cli_parser_contract.py`

**Interfaces:**
- Consumes: `vntyper.scripts.artifact_names.select_best_vcf_file(kestrel_dir) -> str | None`,
  which already exists and is what `pipeline.py` uses. **Do not write a second resolver.**
- Produces: `resolve_vcf_file(output_dir: str, input_dir: Path | None, explicit: Path | None) -> str | None`
  in `cli_report.py`, mirroring the existing `resolve_log_file` shape.

**Already fixed, do not redo.** The signature mismatch the issue body opens with
(`input_files`, `pipeline_version`, `mean_vntr_coverage`) was closed by #179. `cli_report.py`
today passes exactly what `generate_summary_report` accepts. **The only live defect is the VCF.**

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_cli_report.py`, using the existing `bound_calls` fixture:

```python
# ---------------------------------------------------------------------------
# The IGV variant track (#167)
# ---------------------------------------------------------------------------
#
# SPECIFICATION: `handle_report` passed `vcf_file=None` unconditionally, so a
# regenerated report never had an IGV variant track -- even though the run being
# re-reported had written the VCF. `generate_report.run_igv_report` only adds the
# track when `vcf_file` is set and exists, so the omission was silent.


def _kestrel_vcf(run_dir: Path, name: str = "output_indel.vcf.gz") -> Path:
    kestrel = run_dir / "kestrel"
    kestrel.mkdir(parents=True, exist_ok=True)
    vcf = kestrel / name
    vcf.write_bytes(b"")
    return vcf


def test_the_vcf_is_discovered_under_the_output_dir(tmp_path, bound_calls) -> None:
    vcf = _kestrel_vcf(tmp_path)
    cli.main(["report", "-o", str(tmp_path)])
    assert bound_calls[0].arguments["vcf_file"] == str(vcf)


def test_the_vcf_is_discovered_under_the_input_dir_when_one_is_given(tmp_path, bound_calls) -> None:
    run = tmp_path / "run"
    vcf = _kestrel_vcf(run)
    cli.main(["report", "-o", str(tmp_path / "out"), "--input-dir", str(run)])
    assert bound_calls[0].arguments["vcf_file"] == str(vcf)


def test_the_compressed_vcf_wins_over_the_uncompressed_one(tmp_path, bound_calls) -> None:
    """`select_best_vcf_file`'s existing preference, reached from the report path."""
    _kestrel_vcf(tmp_path, "output_indel.vcf")
    gz = _kestrel_vcf(tmp_path, "output_indel.vcf.gz")
    cli.main(["report", "-o", str(tmp_path)])
    assert bound_calls[0].arguments["vcf_file"] == str(gz)


def test_an_explicit_vcf_file_beats_the_discovered_one(tmp_path, bound_calls) -> None:
    _kestrel_vcf(tmp_path)
    explicit = tmp_path / "mine.vcf.gz"
    explicit.write_bytes(b"")
    cli.main(["report", "-o", str(tmp_path), "--vcf-file", str(explicit)])
    assert bound_calls[0].arguments["vcf_file"] == str(explicit)


def test_no_vcf_anywhere_stays_none_and_does_not_raise(tmp_path, bound_calls) -> None:
    """bcftools is optional; a missing VCF is a warning, never an error."""
    cli.main(["report", "-o", str(tmp_path)])
    assert bound_calls[0].arguments["vcf_file"] is None


def test_an_explicit_vcf_file_that_does_not_exist_is_still_passed_through(tmp_path, bound_calls) -> None:
    """The user asked for it by name; `run_igv_report` warns and skips the track."""
    cli.main(["report", "-o", str(tmp_path), "--vcf-file", str(tmp_path / "absent.vcf.gz")])
    assert bound_calls[0].arguments["vcf_file"] == str(tmp_path / "absent.vcf.gz")
```

If `tests/unit/test_cli_parser_contract.py` enumerates the `report` subcommand's arguments, add
`--vcf-file` to that expectation in the same commit.

- [ ] **Step 2: Run and watch it fail**

Run: `python -m pytest tests/unit/test_cli_report.py -v`
Expected: FAIL — `unrecognized arguments: --vcf-file` on the explicit tests, and
`assert None == '<path>'` on the discovery tests.

- [ ] **Step 3: Implement**

`cli_parser.py`, after the `--bam-file` line at 202:

```python
    parser_report.add_argument(
        "--vcf-file",
        type=Path,
        default=None,
        help="Path to the VCF file for the IGV variant track. Discovered under the run directory when omitted.",
    )
```

`cli_report.py` — import the existing resolver and add the handler helper:

```python
from vntyper.scripts.artifact_names import select_best_vcf_file
```

```python
def resolve_vcf_file(output_dir: str, input_dir: Path | None, explicit: Path | None) -> str | None:
    """Choose the VCF the report's IGV panel should show.

    ``handle_report`` passed None unconditionally, so a regenerated report never
    had a variant track even though the run being re-reported had written the VCF
    (#167). :func:`~vntyper.scripts.artifact_names.select_best_vcf_file` is the
    resolver ``pipeline.py`` already uses -- it prefers the compressed VCF and
    returns None rather than raising when neither exists, because bcftools is
    optional.

    An explicit ``--vcf-file`` wins and is passed through even when it does not
    exist: the user named it, and ``run_igv_report`` warns and skips the track.

    The search directory is ``--input-dir`` when given, and ``--output-dir``
    otherwise. That is deliberately wider than the ``--bam-file``/``--bed-file``
    discovery below, which consults ``--input-dir`` only: ``--output-dir`` is
    documented as "Directory containing pipeline results", and
    :func:`resolve_log_file` already reads ``<output-dir>/pipeline.log`` on that
    basis. Widening bam/bed to match would change the behaviour of an argument
    #167 does not name, so it is left alone.

    Args:
        output_dir: The ``--output-dir`` value.
        input_dir: The ``--input-dir`` value, or None.
        explicit: The ``--vcf-file`` value, or None.

    Returns:
        str | None: Path to the VCF, or None when there is none.
    """
    if explicit is not None:
        logger.debug(f"vcf_file set to {explicit} (given on the command line)")
        return str(explicit)
    kestrel_dir = Path(input_dir if input_dir is not None else output_dir) / "kestrel"
    return select_best_vcf_file(kestrel_dir)
```

Then in `handle_report`, replace `vcf_file=None,` with:

```python
        vcf_file=resolve_vcf_file(args.output_dir, args.input_dir, getattr(args, "vcf_file", None)),
```

and update the comment block above the call, which currently explains why three keywords are
absent — add a sentence saying the VCF is now resolved rather than hardcoded.

- [ ] **Step 4: Run**

Run:
```bash
python -m pytest tests/unit/test_cli_report.py tests/unit/test_cli_parser.py \
  tests/unit/test_cli_parser_contract.py tests/unit/test_cli_dispatch.py \
  tests/unit/test_pipeline_artifact_paths.py -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add vntyper/scripts/cli_parser.py vntyper/scripts/cli_report.py \
        tests/unit/test_cli_report.py tests/unit/test_cli_parser_contract.py
git commit -m "fix(cli): resolve the vcf for vntyper report instead of passing none

Closes #167"
```

---

# Integration (the orchestrator, not a subagent)

- [ ] **I1.** `make format && make check-all`. Formatting is *not* gated in CI — run
  `make format-check` yourself.
- [ ] **I2.** `make patch-coverage` — must be ≥ 80% on changed lines.
- [ ] **I3.** Prove all three exit criteria with commands, recorded in `.planning/m3/EVIDENCE.md`:
  E1 a real run against a `chmod a-w` input mount; E2 a byte-compare of two cohort runs in two
  processes; E3 a grep of rendered HTML for an unescaped sample string.
- [ ] **I4.** Cohort output shape moved (C3 changes the `Sample` column for ZIP inputs), so run
  the golden-cohort gate: `python scripts/make_cram_fixtures.py` first, **never**
  `--allow-matrix-drift`. Assert that no genotype field (`Confidence`, `Flag`, `POS`, `REF`,
  `ALT`, `Motif*`, `Depth_Score`, `VID`) appears in `columns_added` or `cells_changed`.
- [ ] **I5.** `make ci-local` only if anything under `.github/workflows/` changed. Nothing in this
  plan does, so this should be a no-op — confirm with `git diff --name-only main...HEAD -- .github/`.
- [ ] **I6.** Open the PR. Address Copilot and Sourcery in **follow-up commits**, never a
  force-push.
- [ ] **I7.** Merge when `ci-success` is green. The 39-69 min Docker Build is not a useful
  pre-merge signal.
- [ ] **I8.** Patch bump to 2.0.x in `vntyper/version.py`, `CITATION.cff` and
  `docs/about/changelog.md`, then `gh release create v2.0.x --target <sha>`. `git tag` is denied
  by design — never route around it.

---

## Self-review

**Spec coverage.** #207 → A1. #216 → A2 (fragments) + A3 (`innerHTML`). #201 → B2. #162 → B2 + B3.
#210 → B3. Worker `os.rmdir` reachability → B4. #206 → C1 (width) + C2 (collision). #205 → C3.
#167 → D1. Exit criteria E1/E2/E3 → I3. No spec section is unclaimed.

**Placeholders.** Three tests in C3 and two in C2 are given as shapes with `...` bodies rather
than full code, with the two adjacent complete tests as the pattern to follow. That is a
deliberate exception to the no-placeholders rule for the repetitive cases (`cram`/`fastq1`
identity variants, sub-directory case) where writing all six out verbatim would obscure the two
that carry the design. Every other step is complete.

**Type consistency.** `js_json_literal(fragment, fallback) -> str` — one name, A2 only.
`validate_bam_file(file_path, cwd=None, log_dir=None)` — B2, consumed by nothing else.
`build_samtools_index_command(*, samtools_path, bam_file, output_bai=None)` and
`resolve_bam_index(in_bam) -> str | None` — B3. `pseudonymized_sample_name(prefix,
original_sample, *, algorithm, length)` — defined C1, consumed C2. `DiscoveredSample(directory,
identity, sort_key)` — defined C3, consumed by `cohort_summary.aggregate_cohort` in C3.
`resolve_vcf_file(output_dir, input_dir, explicit) -> str | None` — D1 only. No name is spelled
two ways.

**Ordering within workstreams.** A1→A2→A3 (independent, any order works; committed in this order
for a readable history). B1→B2→B3→B4 (B1's fixture is imported by B2 and B3, so B1 must land
first). C1→C2→C3 (C2 consumes C1's signature; C3 changes what C2's loop iterates). D1 alone.
**Across** workstreams there is no ordering dependency at all — that is what the disjointness
proof buys.
