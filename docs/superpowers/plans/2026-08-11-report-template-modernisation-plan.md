# Per-sample HTML Report Modernisation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Issue:** #242 · **Design:** `docs/superpowers/specs/2026-08-11-report-template-modernisation-design.md`

**Goal:** Make the per-sample report render the state the pipeline already computes, stop it hiding evidence it narrates, make it operable without a mouse, and make the archived artifact correct offline and on paper.

**Architecture:** `ScreeningSummary` gains a derived `verdict` — from `is_positive` and `matched_rule` only, never from QC — and its existing axes reach the Jinja2 context as a `screening_state` mapping rendered as a banner with three orthogonal state chips. Bootstrap, jQuery and DataTables are removed in favour of CSS custom properties, CSS grid, native `<details>`, and ~150 lines of vanilla sort JS; igv.js is inlined so the report stays one portable file. A shared Jinja2 partial holds the token layer for both templates. Every behavioural claim is pinned by a test that fails first — in the browser tier where the unit tier cannot observe it.

**Tech Stack:** Python 3.10, Jinja2, pandas, pytest, Ruff (line length 120), mypy, vanilla ES2020, CSS custom properties.

## Global Constraints

- Do not change how any clinical state is computed. `build_screening_summary`, `compute_algorithm_result`, `CoverageQC` and the `report_config.json` rule *conditions* keep their semantics and their tests.
- Do not rewrite the 40 configured clinical messages. Splitting one on `<br>` into `headline` / `detail` / `recommendation` is a mechanical migration; the words do not change.
- A `report_config.json` supplying only `message` must keep rendering. The new keys are optional everywhere, including in the loader's type hints.
- Autoescaping stays on. No new `|safe` interpolation beyond the fragments that already carry it. `tests/unit/test_template_escaping.py` governs the new JavaScript exactly as it governs the old.
- Templates obey the repository's 120-character line length; `test_no_template_line_exceeds_the_repository_line_length` enforces it and materially more template source now falls under it.
- The rendered report resolves zero external hosts by default.
- No new Python dependency, no bundler, no build step, no CSS or JS framework.
- Vendored assets are written through `contained_output_path` so `archive_safety.py` guarantees continue to hold.
- Every behavioural edit follows RED, causal failure, minimal GREEN, refactor, verification, commit.
- Verify with `make format && make test-unit` per task and `make check-all` before the PR. Put the `vntyper` env's bin first on `PATH` — the base env and `~/.local/bin` shadow it and `make check-all` will otherwise report on the wrong interpreter.

## File map

| File | Change |
| --- | --- |
| `vntyper/scripts/screening_summary.py` | Add `ScreeningSummary.verdict`; no rule changes |
| `vntyper/scripts/generate_report.py` | Pass `screening_state`, sample identity, run time; drop dead keys; route Kestrel through `escaped_table_html` |
| `vntyper/scripts/report_formatting.py` | Confidence pill class replaces inline colour; flag reason as text |
| `vntyper/scripts/report_config.json` | Split each `message` into an ordered `segments` list, keep `message` verbatim |
| `pyproject.toml`, `MANIFEST.in` | Ship `vntyper.templates` and the pinned igv.js asset |
| `vntyper/scripts/cli_parser.py`, `pipeline.py` | Wire `--report-assets` through both entry paths |
| `vntyper/templates/_report_base.html` | New: token layer, print rules, shared head |
| `vntyper/templates/report_template.html` | Rewrite of the presentation layer |
| `vntyper/templates/cohort_summary_template.html` | Consume the shared partial |
| `vntyper/scripts/report_assets.py` | New: vendor igv.js beside the report |
| `tests/unit/test_generate_report.py` | Update two pinned tests; add verdict, flag-row, identity, empty-state tests |
| `tests/unit/test_report_presentation.py` | New: contrast, no-external-host, print block, `<details>` structure |
| `docs/pipeline/reports.md`, `docs/user-guide/output-files.md` | Document the banner, the assets directory, `--report-assets` |

---

## Task 0: Preconditions (gate — no code)

An adversarial review of this plan returned **unsafe as written**, on the grounds that several tasks depend on decisions nobody has made. Do not start a gated task before its precondition is answered by the domain owners; record each answer in the design's "Unresolved questions" section.

- [ ] **P1 — intended use and labelling.** `README.md:315` says research use only. Is a top-level `finding`/`no-finding`/`indeterminate` label created at all, and with whose approval? *Gates Task 8.*
- [ ] **P2 — artifact shape.** Single self-contained HTML (+~1.3 MB), sibling asset directory, or CDN with SRI? *Gates Task 9.*
- [ ] **P3 — sample identity.** `--sample-name` with a conservative fallback, or basename derivation only? *Gates Task 3.*
- [ ] **P4 — cohort scope.** Does the cohort report keep filtering, search and paging? *Gates Tasks 7 and 9.*
- [ ] **P5 — provenance.** Persist declared assembly, detected assembly and resolved region (summary schema version + legacy fallback), or print only what is already recorded? *Gates Task 3.*
- [ ] **P6 — stage failure vs true negative.** Add per-algorithm execution states here, or defer to #173? *Gates Task 8.*

## Phase 1 — Clinical safety

### Task 1: The computed screening state reaches the template

**Files:**
- Modify: `vntyper/scripts/screening_summary.py`, `vntyper/scripts/generate_report.py`
- Modify: `tests/unit/test_generate_report.py`, `tests/unit/test_screening_summary.py`

**Interfaces:**
- Produces: `ScreeningSummary.verdict -> Literal["finding", "no-finding", "indeterminate"]`
- Produces: context key `screening_state: dict[str, object]` with `kestrel_result`, `advntr_result`, `quality_metrics_pass`, `matched_rule`, `verdict`

- [ ] **Step 1: RED — pin the three verdict states and the discarded flag**

`matched_rule` is the reason this task is first: `screening_summary.py:18-21` documents it as the guard against announcing a negative for a sample two algorithms called positive, and today it never leaves Python.

```python
def test_an_unmatched_rule_is_indeterminate_not_negative() -> None:
    summary = ScreeningSummary(text=FALLBACK_SUMMARY_MESSAGE, is_positive=False, matched_rule=False, ...)
    assert summary.verdict == "indeterminate"

def test_a_finding_with_failing_quality_metrics_is_still_a_finding() -> None:
    """QC is orthogonal. It never suppresses a call.

    The rule table describes exactly this combination as a pathogenic finding with
    low-quality metrics (report_config.json:102, :373), and is_positive is derived
    from the algorithm calls independently of QC by design
    (screening_summary.py:280, :305). An earlier draft of this plan let failed QC
    force "indeterminate", which would have silently reclassified a confirmed call.
    """
    assert ScreeningSummary(is_positive=True, quality_metrics_pass=False, matched_rule=True, ...).verdict == "finding"
```

`quality_metrics_pass` must not appear in the `verdict` expression. It renders as its own always-visible chip beside the call.

Then assert the context, not just the object:

```python
def test_the_screening_state_reaches_the_report(positive_summary) -> None:
    assert 'data-verdict="finding"' in positive_summary.read_text()
```

- [ ] **Step 2: GREEN — derive `verdict`, add `screening_state` to the context**

`verdict` composes existing fields only; it introduces no new axis and no new rule lookup. Add the mapping at `generate_report.py:579-627`. Render `data-verdict` on the summary element so the state is assertable and stylable without a second source of truth.

- [ ] **Step 3: Emit the state as a provenance line**

Under the summary, render `Kestrel: High_Precision · adVNTR: not performed · Coverage QC: PASS` from `screening_state`. This is what makes I2 visible to a reader rather than only to a test.

- [ ] **Step 4: Verify and commit** — `make format && make test-unit`; commit `feat(report): pass the computed screening state to the template`.

### Task 2: Flagged variant rows are never removed

**Files:**
- Modify: `vntyper/templates/report_template.html`, `vntyper/scripts/report_formatting.py`
- Modify: `tests/unit/test_generate_report.py`

- [ ] **Step 1: RED — a flagged variant is in the rendered HTML**

The defect is that a `High_Precision_flagged` report narrates a flagged pathogenic variant and renders an empty table. Measured on a three-variant sample, the same file shows **1 of 3 rows** with the CDNs reachable and **3 of 3** air-gapped, because jQuery fails and the filter never runs. Pin the contradiction directly:

**The obvious test does not work.** Asserting that the flag reason appears in the rendered HTML **passes against the unmodified code**: Kestrel rows are server-rendered by `to_html` (`generate_report.py:528`) and the hiding happens later in the client-side filter — which is exactly why the air-gapped render showed all 3 rows. A RED step that is already green pins nothing.

The defect is only observable as a **visible row count after initialisation**, so this test belongs in the browser tier:

```python
@pytest.mark.browser
@pytest.mark.parametrize("offline", [False, True])
def test_every_flagged_row_is_visible_after_initialisation(rendered_report, offline) -> None:
    page = open_report(rendered_report, offline=offline)
    visible = page.eval("[...document.querySelectorAll('#kestrel_table tbody tr')]"
                        ".filter(r => r.offsetParent !== null).length")
    assert visible == 3   # today: 1 online, 3 offline
```

Pair it with a static invariant stronger than banning one string — a renamed vanilla-JS filter must not be able to reintroduce the behaviour:

```python
def test_no_per_sample_result_row_enters_a_hiding_path(template) -> None:
    """No row-removal or row-hiding construct may touch a results table."""
```

- [ ] **Step 2: GREEN — delete the filter, keep the emphasis**

Remove the `$.fn.dataTable.ext.search.push` predicate (`report_template.html:370-384`). `#toggleFlagged` becomes `#highlightFlagged` and sets a class on the table; the label becomes "Highlight flagged values". The flag cell renders the reason as text plus a glyph, so the meaning survives print, screen readers and a failed script load.

- [ ] **Step 3: Round every displayed number server-side, once**

`applyRounding` (`report_template.html:396-403`) rewrites **every numeric cell in every initialised table** to 4 dp — not just Depth Score. Measured: Depth Score prints `0.01` when jQuery resolves and `0.010012` when it does not. A displayed number must not depend on network state.

Because the browser rule is table-agnostic, removing it changes adVNTR's displayed precision too (mean coverage, p-value — `report_formatting.py:94`). This is an output-format change that downstream HTML scrapers can observe, so it needs explicit per-column server-side formatters, not one Kestrel-specific rule:

```python
def test_no_displayed_number_is_computed_in_the_browser(template) -> None:
    assert "applyRounding" not in template.read_text()

@pytest.mark.parametrize("column", KESTREL_NUMERIC_COLUMNS + ADVNTR_NUMERIC_COLUMNS)
def test_each_numeric_column_matches_the_current_online_rendering(column) -> None:
    """Baseline captured from a report rendered with the CDNs reachable, before this change."""
```

Capture that baseline **before** deleting the browser function. Document any deliberate divergence in `docs/user-guide/output-files.md` as an output-format migration.

- [ ] **Step 5: Verify and commit** — commit `fix(report): stop hiding the flagged variant rows the summary narrates`.

### Task 3: The report identifies itself

**Files:**
- Modify: `vntyper/scripts/generate_report.py`, `vntyper/templates/report_template.html`
- Modify: `tests/unit/test_generate_report.py`

- [ ] **Step 1: RED — CRAM and single-end FASTQ render an input line**

`pipeline_inputs.py:163-169` produces four shapes; the template branches on two.

```python
@pytest.mark.parametrize("inputs", [{"cram": "S1.cram"}, {"fastq1": "S1.fq.gz"},
                                    {"bam": "S1.bam"}, {"fastq1": "a.fq", "fastq2": "b.fq"}])
def test_every_input_shape_names_its_files(tmp_path, inputs) -> None: ...

def test_the_title_carries_the_sample_name(tmp_path) -> None:
    assert "<title>MUC1 VNTR report — S1" in render(inputs={"cram": "S1.cram"})
```

- [ ] **Step 2: GREEN — iterate `input_files`, derive `sample_name`**

Derive from the input basename with suffixes stripped. Put it in `<title>`, `<h1>` and the print header, beside the assay name, reference assembly, VNTR region span, pipeline version, run time and render time.

- [ ] **Step 3: Separate run time from render time**

`report_date` is `datetime.now()` at `generate_report.py:589`, so re-running `vntyper report` over an archived run restamps it. Read the run timestamp from `pipeline_summary.json` and print both, labelled.

- [ ] **Step 4: Verify and commit** — commit `fix(report): name the sample, assay and assembly on every input shape`.

### Task 4: Every control is reachable by keyboard

**Files:**
- Modify: `vntyper/templates/report_template.html`, `vntyper/templates/cohort_summary_template.html`
- Create: `tests/unit/test_report_presentation.py`

**Root cause:** `input[type='checkbox'] { display: none; }` (`report_template.html:71-73`) is unscoped and removes all four inputs from the focus order and the accessibility tree. Patching the selector would fix the symptom; `<details>` removes the pattern.

- [ ] **Step 1: RED — no hidden checkbox, every disclosure is native**

```python
def test_no_interactive_input_is_removed_from_the_accessibility_tree(template) -> None:
    assert not re.search(r"input\[type=.checkbox.\][^{]*\{[^}]*display:\s*none", template.read_text())

def test_every_disclosure_is_a_details_element(template) -> None:
    source = template.read_text()
    assert source.count("<details") == source.count("<summary")
    assert 'class="toggle"' not in source
```

- [ ] **Step 2: GREEN — replace the checkbox accordions with `<details>`/`<summary>`**

Native disclosure is focusable, announced, expanded by find-in-page, and expanded by print. The remaining real checkbox (`#highlightFlagged`) uses a scoped visually-hidden pattern — `position:absolute; opacity:0; width:1px; height:1px` — never `display:none`, with a `:focus-visible` ring on the label.

- [ ] **Step 3: Text beside every glyph**

The flag mark carries its reason as text (Task 2). The coverage and fastp icons from `report_formatting.py:55-64` gain an accessible name. Confidence stops being colour-only.

- [ ] **Step 4: Make variant navigation and column sorting real controls**

"Every control" must mean every control. IGV variant rows are `<tr onclick>` with no keyboard path (`report_template.html:536,548`); changing `event.target` to `event.currentTarget` in Task 9 fixes a crash, not the accessibility. Render an actual button or link per row with an accessible name, focus styling and Enter/Space, and `aria-selected` on the selected row. The vanilla sort headers Task 9 introduces are new interactive controls and need the same: focusable, named, Enter/Space, and `aria-sort` reflecting the current state.

- [ ] **Step 5: Add the viewport meta, landmarks and reduced motion** — `<meta name="viewport">` is absent from **both** templates, so mobile lays out at 980px; add it to both in this commit. Wrap the body in `<main>`, give `header_warning` `role="alert"`, and add a `prefers-reduced-motion` rule — the disclosures animate unconditionally today (`:74`, `:105`) and this is the one identified defect most likely to survive a rewrite unnoticed.

- [ ] **Step 6: Browser-tier keyboard traversal**

```python
@pytest.mark.browser
def test_every_control_is_reachable_by_tab(rendered_report) -> None:
    """Includes the IGV-present and IGV-failure cases."""
```

Source assertions cannot see focusability; `tests/unit/test_template_escaping.py:17` says the unit tier has no browser at all.

- [ ] **Step 7: Verify and commit** — commit `fix(report): make every report control keyboard-operable`.

### Task 5: Authored empty states

**Files:**
- Modify: `vntyper/scripts/generate_report.py`, `vntyper/templates/report_template.html`
- Modify: `tests/unit/test_generate_report.py`

- [ ] **Step 1: RED — a negative run does not render `None`/`NaN`**

```python
def test_a_negative_kestrel_run_renders_a_sentence_not_a_none_row(tmp_path) -> None:
    html = render_negative()
    assert "None None" not in html and "NaN" not in html
    assert "No variant detected by Kestrel" in html

def test_an_absent_kestrel_step_does_not_emit_a_zero_column_table(tmp_path) -> None:
    assert '<thead><tr style="text-align: right;"></tr></thead>' not in render_without_kestrel()
```

- [ ] **Step 2: GREEN** — route the Kestrel table through `escaped_table_html` (`generate_report.py:528` bypasses the helper written for exactly this, `report_formatting.py:278-279`), suppress the `output_empty_result` placeholder row from display, and mirror the adVNTR branch's authored sentence.

- [ ] **Step 3: Remove the unreachable branch** — `generate_report.py:537` sets an adVNTR message the template's `{% if %}` at line 326 can never render, while line 331 prints a differently worded one. Delete one.

- [ ] **Step 4: Verify and commit** — commit `fix(report): author the empty states instead of printing None and NaN`.

### Task 6: Print is a build target

**Files:**
- Modify: `vntyper/templates/report_template.html`
- Modify: `tests/unit/test_report_presentation.py`

- [ ] **Step 1: RED — the print block exists and unsets the clamp**

```python
def test_the_report_has_a_print_stylesheet(template) -> None:
    assert "@media print" in template.read_text()

def test_print_unsets_the_truncation_clamp(template) -> None:
    block = extract_media_print(template.read_text())
    assert "max-width: none" in block and "white-space: normal" in block
```

- [ ] **Step 2: GREEN** — unset the 150px clamp and `nowrap`, hide interactive chrome, add a running header with sample, assay, assembly, version and page numbers, and print the log as a pointer rather than a wall of DEBUG output.

**CSS alone cannot open a closed `<details>`** — `open` is an attribute, not a style, and `content-visibility` tricks do not reliably restore it in print. Set `open` server-side for print-relevant sections, or toggle it in an `onbeforeprint` handler with a no-JS fallback of rendering those sections open by default. Decide which, and say so here.

- [ ] **Step 3: Browser-tier print assertion**

The substring test above is a tripwire, not proof: a stylesheet can satisfy every source assertion while the PDF still omits closed content, page numbers or long cells.

```python
@pytest.mark.browser
def test_the_printed_pdf_is_a_complete_record(rendered_report) -> None:
    text = print_to_pdf_and_extract_text(rendered_report)
    assert SAMPLE_NAME in text and FULL_MOTIF_SEQUENCE in text
    assert text.count(VARIANT_ROW_MARKER) == 3
    assert "Page 1 of" in text and DEBUG_LOG_LINE not in text
```

Name the supported print renderer and version in the docstring — page counters and running headers differ materially between engines.

- [ ] **Step 4: Verify and commit** — commit `feat(report): make the archived PDF a correct record`.

---

## Phase 2 — Presentation

### Task 7: One token layer, two reports

**Files:**
- Create: `vntyper/templates/_report_base.html`
- Modify: both templates
- Modify: `vntyper/scripts/report_formatting.py`, `tests/unit/test_report_presentation.py`

- [ ] **Step 1: RED — contrast is computed, not asserted from a table**

```python
def test_every_token_pair_meets_wcag_aa() -> None:
    for name, (fg, bg) in TOKEN_PAIRS.items():
        assert contrast_ratio(fg, bg) >= 4.5, f"{name}: {contrast_ratio(fg, bg):.2f}:1"
```

Seed it with the failures this replaces: white on `lightskyblue` (`report_template.html:34-38`), and the `Low_Precision` orange and `High_Precision` red from `CONFIDENCE_STYLES` (`report_formatting.py:109-113`).

- [ ] **Step 2: GREEN — define the scale in `_report_base.html`**

`--vn-finding`, `--vn-no-finding`, `--vn-indeterminate`, `--vn-risk-high`, `--vn-risk-low`, plus surface, text and border tokens. Confidence moves to a pill class, so red is reserved for "something is wrong" everywhere on the page.

- [ ] **Step 3: Delete the duplicated CSS** — the truncation block is byte-identical across the two templates (`report_template.html:44-63`, `cohort_summary_template.html:18-34`). `.table-container` is used in the per-sample report and defined only in the cohort one; both facts go away with the shared partial.

- [ ] **Step 4: Update the pinned test** — `test_the_confidence_column_is_colour_coded` asserts the inline `color:red` this task removes. It becomes an assertion on the pill class, in this commit.

- [ ] **Step 5: Verify and commit** — commit `refactor(report): share one token layer between both reports`.

### Task 8: The verdict banner

**Files:**
- Modify: `vntyper/scripts/report_config.json`, `vntyper/scripts/screening_summary.py`, `vntyper/templates/report_template.html`
- Modify: `tests/unit/test_generate_report.py`

- [ ] **Step 1: RED — each verdict renders its banner, and a message-only config still works**

```python
@pytest.mark.parametrize("verdict", ["finding", "no-finding", "indeterminate"])
def test_each_verdict_renders_its_banner(tmp_path, verdict) -> None: ...

def test_a_config_with_only_message_keys_still_renders(tmp_path) -> None:
    """The new headline/detail/recommendation keys are optional, forever."""
```

- [ ] **Step 2: Migrate the 40 messages into an ordered `segments` list**

Not into fixed `headline`/`detail`/`recommendation` fields. Audited: **31 rules carry one `<br>` (two segments) and 9 carry two (three segments)**. A three-field schema would force the migration to invent an empty field or to decide subjectively whether a combined sentence is "detail" or "recommendation" — reinterpreting clinical text, which the constraints forbid.

Split each `message` on `<br>` into `segments`, keep `message` verbatim, and assert the round trip on rendering rather than on shape:

```python
@pytest.mark.parametrize("rule", ALL_40_RULES, ids=rule_id)
def test_rendering_reproduces_the_exact_legacy_message(rule) -> None:
    assert render_segments(rule) == rule["message"]
```

Produce a rule-by-rule migration table for clinical review of the segment boundaries. Render each segment as a separate autoescaped element — do not mark assembled message HTML safe.

- [ ] **Step 3: GREEN — render the banner directly under `<h1>`**

Three states from `screening_state.verdict`, each with a distinct tint, a 4px left rule, an icon **with a text label**, `headline` at ≥1.5rem, `detail` as body copy, and `recommendation` as an explicit next step. Coverage and fastp move below it. Update `test_a_positive_screening_is_styled_as_a_finding` to assert the banner state in this commit.

- [ ] **Step 4: Verify and commit** — commit `feat(report): lead with the clinical verdict`.

### Task 9: Zero external hosts

**Files:**
- Modify: `vntyper/templates/report_template.html`, `vntyper/templates/cohort_summary_template.html`
- Create: `vntyper/scripts/report_assets.py`, `tests/unit/test_report_assets.py`
- Modify: `vntyper/scripts/cli_report.py`

- [ ] **Step 1: RED — the rendered report resolves nothing off-machine**

```python
def test_the_rendered_report_has_no_external_host(positive_summary) -> None:
    assert not re.search(r'(?:src|href)=["\']https?://', positive_summary.read_text())
```

Eleven such tags exist today across the two templates, over four origins, with no `integrity` and no fallback.

- [ ] **Step 2: GREEN — remove Bootstrap, jQuery and DataTables**

1.85 MB raw / 507 KB gzip today. DataTables 1.10.21 is the one genuinely vulnerable pin (CVE-2020-28458, CVE-2021-23445, SNYK-JS-DATATABLESNET-598806; all cleared only at ≥1.11.3) and is EOL. Bootstrap 4.5.2 is EOL with no live CVE. jQuery 3.5.1 is **not** affected by CVE-2020-11022/11023 — do not cite those; the argument for removing jQuery is D0 and weight, not a vulnerability.

Replace with CSS grid, custom properties and ~150 lines of vanilla JS for column sort and the highlight toggle. Pagination is not reintroduced; its absence is what makes print correct by default.

- [ ] **Step 3: Inline igv.js, and load it only when there is something to show** *(gated on P2)*

Measured: `igv.min.js` is 1.31 MB of the 1.85 MB total and is fetched unconditionally from `<head>` — verified loading on a sample with no IGV data, which renders the "not available" fallback. Load it from the `sessionDictionary` branch, so a report without alignments loads nothing.

**Do not copy it to a sibling `assets/` directory through `contained_output_path`.** That helper rejects any name where `basename(name) != name` (`output_paths.py:20`), so `assets/igv.min.js` raises `ValueError` — and a sibling directory would break `--report-file`'s promise that the report is one file (`cli_parser.py:197`, `docs/cli/report.md:53`). Inline it into a `<script>` element instead. `--report-assets {inline,cdn}`; `cdn` keeps today's delivery and gains `integrity` + `crossorigin`.

Add a failure state distinct from the empty-session state at `report_template.html:487-492`: a failed load currently throws `ReferenceError`, so `initIGV` is never reached and the reader sees blank space with no explanation.

- [ ] **Step 3b: Ship the new files, and wire the flag through both entry points**

Two release-only failures that editable-checkout tests cannot catch:

- `pyproject.toml:135` declares package data only for `vntyper.scripts`, and `MANIFEST.in:9-10` names exactly the two existing templates. Add `vntyper.templates` (base partial included) and the pinned `igv.min.js` with its version, source URL, SHA-256 and licence/notice. Verify by building a wheel **and** an sdist, installing each in a clean environment, and rendering both report types from the installed package.
- `--report-assets` cannot work by touching `cli_report.py` alone: the parser is `cli_parser.py:177` and normal pipeline runs reach the generator by a separate path (`pipeline.py:568`). Wire it through the parser, the handler, `run_pipeline` and `generate_summary_report`, with signature-binding tests on **both** entry paths (extend `test_cli_report.py:41`) and an archive-content test per mode. Include a CLI compatibility matrix proving `--report-file`, `--bed-file`, `--bam-file`, `--vcf-file`, `--reference-fasta` and `--flanking` all still reach the generator.

- [ ] **Step 4: Fix the row-click crash** — `event.target.parentElement.id` (`report_template.html:553`) yields `tbody`'s empty id when the click lands on the `<tr>` itself, so `getElementById("")` returns `null` and `classList` throws, killing IGV navigation for the session. Use `event.currentTarget`.

- [ ] **Step 5: Verify and commit** — commit `feat(report): make the archived report self-contained`.

### Task 10: Context and template agree

**Files:**
- Modify: `vntyper/scripts/generate_report.py`, `tests/unit/test_generate_report.py`

**Deprecate; do not delete.** `config.json:138` makes the template directory configurable and both renderers load from it (`generate_report.py:556`, `cohort_summary.py:295`), so a third-party template may reference exactly the keys the shipped one ignores. "Unused by the shipped template" is not "not part of the renderer API".

- [ ] **Step 1: RED — every *new* key is referenced, and the resolver follows includes**

```python
def test_every_new_context_key_is_referenced_by_a_shipped_template() -> None:
    unused = (context_keys() - DEPRECATED_KEYS) - jinja_referenced_names_recursive(REPORT_TEMPLATE)
    assert unused == set(), f"dead context keys: {sorted(unused)}"
```

`jinja_referenced_names_recursive` must resolve `include`/`extends`, or it reports false positives the moment `_report_base.html` lands.

- [ ] **Step 2: GREEN** — keep passing the six `*_color` keys and `igv_content`; consume the colour keys in the token layer where they fit. Record them in `DEPRECATED_KEYS`, document them with a removal release in `docs/pipeline/reports.md`, and add a test rendering a representative legacy custom template that uses them.

- [ ] **Step 3: Verify and commit** — commit `refactor(report): delete the context keys nothing renders`.

### Task 11: Documentation

**Files:**
- Modify: `docs/pipeline/reports.md`, `docs/user-guide/output-files.md`, `docs/cli/report.md`

- [ ] **Step 1** — document the verdict banner and its three states, the provenance line, the assets directory, `--report-assets`, and the print behaviour.
- [ ] **Step 2** — state explicitly that the per-sample report never hides a variant row, and that flag filtering remains a cohort-triage affordance.
- [ ] **Step 3: Verify and commit** — `make check-all`; commit `docs: describe the report verdict banner and self-contained assets`.

---

## Verification

- Per task: `make format && make test-unit`.
- Before the PR: `make check-all` with the `vntyper` env's bin first on `PATH`.
- **Determinism (the D0 acceptance check).** Render one `High_Precision_flagged` three-variant sample. Open it twice — once normally, once with the network blocked (`chrome --host-resolver-rules="MAP * 127.0.0.1:1,EXCLUDE 127.0.0.1"`, which is how the defect was measured). The rendered row count and every displayed number must be identical. Today they are not: 1 of 3 rows online versus 3 of 3 offline, and `0.01` versus `0.010012`.
- **Browser tier.** The unit tier has no browser and no JavaScript engine (`tests/unit/test_template_escaping.py:17`), so three central claims are invisible to it and must be asserted in a browser: visible flagged-row count (online and offline), printed-PDF completeness, and keyboard traversal. A browser-tier marker and its runner are a prerequisite of Tasks 2, 4 and 6 — decide in Task 0 whether it is Playwright, or the CDP-over-headless-Chrome approach already used to measure these defects.
- **Packaging.** Build a wheel and an sdist, install each in a clean environment, and render both report types from the installed package. Editable-checkout tests cannot catch a missing template or asset.
- **No unapproved network request.** Browser network interception over both templates and both asset modes, failing on any request not on the allowlist — a `src`/`href` grep misses CSS `url()`, `@import`, `srcset`, dynamic JS and IGV session URLs.
- Manual: confirm the narrated flagged variant is visible; print to PDF and confirm the motif sequence, every variant row and the sample name are present; tab through every control and reach all four disclosures; check the report at 390px.

## Rollout

One issue, one PR, commits as above. Phase 1 and Phase 2 do not split into separate PRs: they touch the same template and a split would leave `main` with a half-migrated stylesheet.

## Section preservation

Before Task 7 and Task 9, write the preservation matrix and give each row a test: screening summary, **cross-match summary** (`report_template.html:308`, `test_generate_report.py:433`), BAM-header warning, coverage basic and detailed, fastp metrics, Kestrel results, adVNTR results, IGV, pipeline log, report date, pipeline version, input files. A broad rewrite with no explicit preservation criterion drops sections silently.

Decide explicitly rather than dropping: DataTables supplies **search** and **"showing N of M"** on both reports today. Paging is deliberately not reintroduced — its absence is what makes print correct — but that does not require removing search from screen media. If search is kept, implement it accessibly and render the total/visible count server-side.

`tests/unit/test_cohort_summary_oracle.py:411` holds a whole-document fingerprint of the cohort report. Any cohort change must re-baseline it as a reviewed, deliberate step, not as an incidental test fix.

## Unresolved questions

The six blocking ones are now gated preconditions — see **Task 0**. One remains open and blocks nothing: land `confidence_tier` as `None` now so #173 is purely additive later?
