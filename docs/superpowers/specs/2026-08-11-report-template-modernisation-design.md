# Per-sample HTML report: clinical-safety fixes and presentation modernisation

**Date:** 2026-08-11
**Issue:** #242
**Status:** Proposed — blocked on preconditions P1-P6
**Review:** Adversarially reviewed by Codex gpt-5.6-sol (high) on 2026-08-11; verdict "unsafe as written". Findings verified and incorporated: the verdict definition no longer lets QC suppress a call, the message migration uses ordered segments, asset packaging replaces the sibling-directory plan, context keys are deprecated rather than deleted, and the browser test tier is now explicit.

## Problem

`vntyper/templates/report_template.html` is the artifact a reader consults to judge
whether a sample carries a pathogenic MUC1 VNTR frameshift, and it is the artifact that
gets zipped (`archive_safety.py:400`) and forwarded. The pipeline underneath it is
careful. The template is not: it is a 2015 Bootstrap-4 + DataTables table skin carrying
a burden it was never designed for.

**Intended use is a precondition, not a detail.** `README.md:315` states the tool is for
**research use only**. This design therefore does not create a clinical verdict surface.
It makes the report state what the pipeline computed, legibly and reproducibly, and it
carries the research-use statement on screen and in print. Every user-facing label
remains configuration owned by the domain authors, not wording invented in
implementation. Where this document says "verdict" it means *the state the existing rule
table resolved*, never a new clinical classification. See the preconditions in
"Unresolved questions" — P1 blocks the banner work.

One root cause produces most of the defects below. **The report computes a structured
clinical state and then discards the structure**, rendering it as one prose paragraph
in a grey box.

`screening_summary.py:248` resolves three axes —

| Axis | Values |
| --- | --- |
| `kestrel_result` | `High_Precision`, `High_Precision_flagged`, `Low_Precision`, `Low_Precision_flagged`, `negative` |
| `advntr_result` | `positive`, `positive flagged`, `negative`, `none` |
| `quality_metrics_pass` | `true`, `false` |

— against a 40-rule table in `report_config.json`. Every rule's payload is a ~240-character
sentence built from three `<br>`-joined clauses: the call, the QC state, the recommended
action. `generate_report.py:579-627` passes only the flattened sentence and a boolean.
The template (`report_template.html:304-306`) renders it as `<p class="summary-box">`
with `font-weight: bold` when the boolean is set. That bold is the entire visual
distinction between "this patient has a pathogenic MUC1 frameshift" and "this patient
does not."

### Defects

**D0 — The report is not a fixed record: its clinical content depends on network
reachability.** Measured on a rendered three-variant sample, opened twice from the same
file. With the CDNs reachable, the Kestrel table shows **1 of 3 variants** — the flag
filter runs and DataTables reports "(filtered from 3 total entries)". Air-gapped, jQuery
fails to resolve, the inline script throws `ReferenceError: $ is not defined` at its
first statement (`report_template.html:370`), the whole `<script>` block aborts, and the
same table shows **all 3 variants**. The Depth Score differs too: `0.01` online,
`0.010012` offline, because `applyRounding` (`:396-403`) never runs.

This is the defect that orders the rest. A clinical artifact whose visible evidence and
whose printed numbers depend on whether a hospital network reaches
`code.jquery.com` cannot be archived, cited or reproduced. D1 and D6 below are its two
halves; neither is fully fixable without the other.

**D1 — The default view removes evidence rows.** `#toggleFlagged` ships unchecked, and
the DataTables filter at `report_template.html:370-384` returns `false` for any row whose
Flag is not `Not flagged` / `Not applicable` / empty. The row leaves the DOM. A
`High_Precision_flagged` report therefore prints a sentence that narrates a flagged
pathogenic variant (`report_config.json`, rules 9-16) above an empty table, with no
explanation beyond DataTables' 12.8px "(filtered from 1 total entries)". `docs/pipeline/flagging.md`
documents advisory flags as *reported*; the default view contradicts that contract.

**D2 — The verdict is unranked.** The screening summary sits below the heading, date,
version, input files, the BAM-header block, two toggles, the coverage section and, when
fastp ran, a five-row QC table — roughly 40-60% down the page. Nothing in the first
viewport answers the question the reader came with. The accordion label for the pipeline
DEBUG log (`1.2rem`, uppercase, monospace, grey slab, `report_template.html:74-87`)
outweighs the diagnosis.

**D3 — The document cannot identify itself, and cannot print.** `<title>` and `<h1>` are
the constant string `Summary Report` (lines 5, 149). The only identifier is `input_files`,
and the template branches on two of the four shapes `pipeline_inputs.py:163-169` produces —
CRAM input and single-end FASTQ both render an empty `Input Files:` line. `<title>` is
the same string for every report the lab produces, which is precisely the affordance a
screen-reader user relies on to tell two open reports apart.

Neither template contains a single author media rule — confirmed in the browser, not
grepped: `document.styleSheets` reports `mediaRules: []`. So an archived PDF truncates
the ~60 bp motif sequence at `max-width: 150px` with a `:hover`-only reveal, prints 10 of
N variants at the DataTables default page length, and prints collapsed sections as
nothing. There is likewise no `prefers-reduced-motion` and no `prefers-color-scheme`.

The same absence breaks the phone: with no `<meta name="viewport">`, the layout viewport
falls back to 980px, so at a 390px device width the page is scaled by 0.398 and the
12.8px body text renders at an effective **~5.1px**. The 11-column table overflows by
237px because `.table-container` (`:322,328`) has no rule in this template at all — it is
defined only in the cohort one — so there is no `overflow-x` wrapper. `th { white-space:
pre-wrap }` then stacks headers one character per line: the `Motif` header measures 279px
tall in a 48px column.

**D4 — Every custom control is keyboard-dead and invisible to assistive technology.**
`input[type='checkbox'] { display: none; }` (line 71) is unscoped. It hits all four
interactive inputs: `#toggleFlagged`, `#toggleDetailedCoverage`, `#collapsible`,
`#logToggle`. `display: none` removes an element from the focus order *and* the
accessibility tree, and Bootstrap's switch still renders because its visuals are
`::before`/`::after` on the sibling label — so the defect is invisible in visual QA. The
consequence: "Show detailed coverage" is the only route to the Coverage QC verdict, the
min/max, and the uncovered-base count, and no keyboard or screen-reader user can reach it.

Measured in the browser: all four inputs report `width: 0, height: 0, focusable: false`.
**Neither template contains a single `tabindex`, `role`, `aria-*` or key handler.** The
IGV variant rows (`:548-559`) are `<tr>` elements with a click handler and no keyboard
path. The hover-to-reveal truncation (`:48-63`) has no `:focus` or `:focus-within`, so a
~60 bp motif sequence is unreachable without a pointer — and on touch. The flag tooltips
(`:428-436`) carry no `tabindex`, and Bootstrap moves `title` into `data-original-title`,
removing the native fallback as well. WCAG 2.1.1, unambiguously.

**D5 — Red means three contradictory things, and none of them passes contrast.**
`CONFIDENCE_STYLES` (`report_formatting.py:109-113`) paints `High_Precision` red — the
*most trustworthy* call. The flag glyph is red for `✖` — *untrustworthy*.
`header_warning` is red (`report_template.html:177`) — a *setup problem*. All three can
appear in one viewport. Measured ratios against AA's 4.5:1 for normal text, on a 12.8px
body (`body { font-size: 80% }`):

| Element | Pair | Ratio | AA |
| --- | --- | --- | --- |
| `th` label | `#FFFFFF` on `#87CEFA` | **1.72** | fail (fails 3:1 large too) |
| `Low_Precision` | `#FFA500` on `#FFFFFF` | **1.97** | fail |
| `Low_Precision`, striped row | `#FFA500` on `#F2F2F2` | **1.76** | fail |
| `header_warning` | `#FF0000` on `#F9F9F9` | **3.80** | fail |
| `High_Precision`, red `✖` | `#FF0000` on `#FFFFFF` | **4.00** | fail |
| same, striped row | `#FF0000` on `#F2F2F2` | **3.57** | fail |
| green `✓` | `#008000` on `#FFFFFF` | 5.14 | pass |

Every failing pair is one that carries clinical meaning. The one that passes is the
"nothing is wrong" glyph.

**D6 — The report depends on four uncontrolled hosts but ships as an archive.**
Eleven CDN tags across four origins, across both templates — **zero `integrity`, zero
`crossorigin`, zero fallback, no vendored copy** (`report_template.html:7,9,145,362,363,365`;
`cohort_summary_template.html:6,7,85,87,89`). Measured payload:

| Asset | Pinned | Raw | Gzip | Current |
| --- | --- | --- | --- | --- |
| igv.min.js | 3.0.2 | 1,310,337 | 358,265 | 3.8.5 |
| bootstrap.min.css | 4.5.2 | 160,302 | 24,900 | 5.3.8 |
| jquery-3.5.1.js *(unminified)* | 3.5.1 | 287,630 | 84,374 | 4.0.0 |
| jquery.dataTables.min.js | 1.10.21 | 84,647 | 27,751 | 3.0.1 |
| bootstrap.bundle.min.js | 4.5.2 | 80,927 | 22,350 | — |
| jquery.dataTables.min.css | 1.10.21 | 13,892 | 1,838 | — |
| **Total** | | **1.85 MB** | **507 KB** | |

`igv.min.js` is 1.31 MB of that and is fetched unconditionally from `<head>` — verified
loading on a sample with no IGV data at all, which renders the "not available" fallback.
jQuery is served unminified, costing 199 KB over the minified build for nothing.

Advisories, each verified against OSV/NVD/GitHub Advisory DB rather than assumed:
**DataTables 1.10.21 is the one genuinely vulnerable pin** (EOL Feb 2024) — CVE-2020-28458
(prototype pollution, CVSS 7.3), CVE-2021-23445 (XSS, fixed 1.11.3), and
SNYK-JS-DATATABLESNET-598806; clearing all three needs ≥1.11.3. Bootstrap 4.5.2 is EOL
(2023-01-01) with **no live CVE** — the 4.x CVEs were fixed at or before 4.3.1 and the
2024 pair was rescinded by NVD, though distro-keyed scanners still flag 4.5.2.
**jQuery 3.5.1 is not affected** by CVE-2020-11022/11023 — 3.5.0 is the fix. igv.js 3.0.2
has no advisories, and is eight minors stale.

Beyond D0, the air-gapped failure was measured end to end: both toggle switches render as
bare unclickable text (the inline `display: none` still applies but the Bootstrap CSS that
draws the control is gone), the flag glyphs revert to raw text, and the IGV fallback
message **never appears** — `initIGV` is never reached, so the reader sees blank space
with no explanation. The designed empty state at lines 487-492 covers *empty session
data*, not *failed library load*.

**D7 — Computed state is discarded at the last hop.** `ScreeningSummary.matched_rule`
exists, per its own docstring (`screening_summary.py:18-21`), so that a report cannot
announce a negative screening for a sample two algorithms called positive. It is logged
(`:297`) and never passed to the template, so the fallback sentence "The screening was
negative (no valid Kestrel or adVNTR data)." renders identically to a genuine negative.
Six `*_color` keys and `igv_content` are likewise computed in `generate_report.py:596-627`
and referenced nowhere in the template.

**D8 — Two authored empty states, and several unauthored ones.** The IGV fallback
(lines 482-493) is written by a human and is good. By contrast a negative sample renders
`None None None None None None None None NaN Negative` under the heading "Kestrel
Identified Variants", because `output_empty_result` writes the literal string `None`
into all nine columns and `build_kestrel_frames` coerces `Depth_Score` through
`pd.to_numeric(errors="coerce")` (`generate_report.py:264`). The commonest report in any
cohort looks like a crashed pipeline.

**D9 — The two templates are a copy-paste pair.** `report_template.html:44-63` and
`cohort_summary_template.html:18-34` are the same truncation block, comment included.
Both load the same four CDNs. `.table-container` is used twice in the per-sample report
(lines 322, 328) and defined only in the cohort template, so the variant tables sit in a
styling hook that does nothing.

## Goals

- The reader learns the call, its confidence, its concordance and its QC status in the
  first viewport, without scrolling and without reading a paragraph word by word.
- No configuration of the pipeline can produce a report that hides a variant row it
  simultaneously narrates in prose.
- Every control is operable by keyboard and announced by assistive technology; every
  meaning carried by colour or glyph is also carried by text.
- The archived artifact renders and prints correctly with no network access, and names
  the sample, assay, reference assembly, pipeline version and run time on paper.
- The template renders the state the pipeline computes, rather than re-deriving
  presentation from prose or discarding the state.
- One token layer and one base template serve both the per-sample and the cohort report.

## Non-goals

- Do not change how any clinical state is computed. `build_screening_summary`,
  `compute_algorithm_result`, `CoverageQC` and the `report_config.json` rule conditions
  keep their current semantics and their current tests. This work changes what the state
  is *allowed to express*, not what it is.
- Do not change the wording of the 40 configured clinical messages. They are the
  authors' clinical text; this work decomposes them into fields, it does not rewrite them.
- Do not add a graded confidence tier. That is #173 and it lands in the same context
  dictionary; this design leaves the seam for it and says where.
- Do not redesign the cohort summary's plots or its cohort-triage filtering semantics.
  Hiding flagged rows is defensible for cohort triage; it is not defensible for a
  single-patient diagnostic read.
- Do not introduce a build step, a bundler, a CSS framework or a JavaScript framework.

## Architecture

### 1. The template renders state, not prose

`ScreeningSummary` gains no new computation, but the context gains the fields it already
holds:

```python
"screening_state": {
    "kestrel_result": screening.kestrel_result,
    "advntr_result": screening.advntr_result,
    "quality_metrics_pass": screening.quality_metrics_pass,
    "matched_rule": screening.matched_rule,
    "verdict": screening.verdict,      # finding | no-finding | indeterminate
},
```

`verdict` is a derived property on `ScreeningSummary`, and its definition is
deliberately narrow:

```
indeterminate  when matched_rule is False        # no rule matched; the state is unknown
finding        when is_positive
no-finding     otherwise
```

**`quality_metrics_pass` must not enter this expression.** An earlier draft of this
design let failed QC force `indeterminate`, which would have silently reclassified a
confirmed pathogenic call with poor coverage as "state unknown" — contradicting the rule
table, which describes exactly that combination as a finding with low-quality metrics
(`report_config.json:102`, `:373`), and contradicting `is_positive`, which is derived
from the algorithm calls independently of QC by design (`screening_summary.py:280,305`).
Call status, QC status, concordance and rule-match status stay **orthogonal**: QC renders
as its own always-visible chip beside the call, never as a modifier that suppresses it.

Even this narrow `verdict` is a new top-level label on a research artifact, so it is
gated on precondition P1 below. If the answer there is "no new label", the banner renders
`headline` alone with the three-chip row and no verdict word; nothing else in this design
changes.

`report_config.json` rules gain one optional sibling key to `message`:

```json
{ "conditions": {...}, "segments": ["...", "..."], "message": "..." }
```

An ordered list, **not** fixed `headline`/`detail`/`recommendation` fields. The messages
are not uniformly three-clause: of the 40 rules, **31 carry one `<br>` (two segments) and
9 carry two (three segments)**. A fixed three-field schema would force the migration to
invent an empty field or to decide subjectively whether a combined sentence is "detail"
or "recommendation" — which is precisely the reinterpretation of clinical text this
design forbids. `segments[0]` is the headline; the remainder render as body copy in
order. `message` stays the authoritative fallback and the round-trip test asserts that
**rendering reproduces the exact legacy message**, not that any rule has a given number
of segments.

Compatibility note, corrected: `load_report_config()` (`screening_summary.py:102`) takes
no path and always loads the packaged file, so there is **no supported user-supplied
report-config API today**. The optional-key design is therefore forward compatibility for
a config that ships with the package, not preservation of an existing injection point. If
such a path is ever added it needs a schema and a trust boundary first — the current
`message|safe` interpolation would turn configuration into executable report content.

### 2. Zero external hosts

The presentation layer drops Bootstrap, jQuery and DataTables entirely. The report uses
a container, a table class, a custom switch and a grid from Bootstrap; sort, search and
pagination from DataTables; and a selector engine from jQuery. Against the ~1.2 MB those
three cost, the replacement is roughly 400 lines of CSS custom properties and grid, and
roughly 150 lines of vanilla JS for column sort and the filter toggle. Three specific
gains beyond weight:

- `<details>` / `<summary>` replaces the `input[type=checkbox]` + label accordion hack,
  which removes D4's root cause structurally rather than by patching the selector, and
  is expanded by the browser's own find-in-page and print paths.
- No pagination means print correctness is the default rather than an
  `onbeforeprint` workaround.
- The jQuery/DataTables advisory surface leaves the artifact permanently, which matters
  for a file that outlives the release that produced it.

igv.js stays — it is the scientific payload, not chrome — and it is loaded from the
`sessionDictionary` branch rather than from `<head>`, so a report with no alignments
fetches nothing at all. It is also the whole of the packaging problem, and the naive
answer does not work:

- **A sibling `assets/` directory cannot go through `contained_output_path`.** That
  helper rejects any name where `basename(name) != name` (`output_paths.py:20`), so
  `assets/igv.min.js` raises `ValueError`. A directory-aware confined-output helper with
  symlink checks and atomic writes would have to be written first.
- **A sibling directory also breaks the artifact.** `--report-file` names one HTML file
  today (`cli_parser.py:197`) and the docs present that file as *the* report
  (`docs/cli/report.md:53`). Making a sibling directory mandatory means renaming,
  moving, emailing or attaching the HTML alone silently breaks IGV — an output-format
  migration, not an implementation detail.

The design therefore takes **single-file inlining** as the default: igv.js is emitted
into a `<script>` element in the document, so the report stays exactly one portable file
and D0 is closed without a bundle contract. `--report-assets {inline,cdn}` selects; `cdn`
keeps today's delivery for size-sensitive batch runs and gains `integrity` and
`crossorigin`. Both modes get an authored failure state distinct from the empty-session
state. The cost is roughly 1.3 MB per report, which precondition P2 must accept.

Packaging is part of this, and is currently missing: `pyproject.toml:135` declares
package data only for `vntyper.scripts`, and `MANIFEST.in:9-10` names exactly the two
existing templates. A new base template and a vendored `igv.min.js` would be absent from
wheels and sdists while every editable-checkout test passed — a release-only failure. The
vendored asset needs a pinned version, source URL, SHA-256 and licence/notice, and both
packaging declarations need updating, verified by installing a built wheel *and* sdist in
a clean environment and rendering both report types from the installed package.

### 3. Evidence is never removed, only de-emphasised

In the per-sample report the flag toggle stops filtering. Flagged rows always render,
carrying a left rule and the flag reason as **text**, not as a tooltip on a glyph. The
toggle becomes "Highlight flagged values", which sets a class. Should filtering be
retained for any surface, the table must be preceded by an unmissable, non-dismissible
count line — "1 flagged variant is hidden" — rendered server-side from the row count,
not by DataTables' footer text.

### 4. Colour is a semantic scale, not an ad-hoc palette

One token layer in a shared partial, consumed by both templates:

- `--vn-finding`, `--vn-no-finding`, `--vn-indeterminate` for verdict states.
- `--vn-risk-high` / `--vn-risk-low` for QC and flag outcomes.
- Confidence stops using the alert palette. `High_Precision` becomes a weight and a
  pill label ("High precision"), so red is reserved for "something is wrong" everywhere
  on the page.

Every token pair is checked against WCAG 2.2 AA (4.5:1 for normal text) by a test, not
by eye.

### 5. Print is a build target

`@media print` unsets the truncation clamp and `nowrap`, forces every `<details>` open
except the pipeline log, hides interactive chrome, and adds a running header carrying
sample, assay, reference assembly, pipeline version, run time and page numbers. The
pipeline log prints as a one-line pointer to the HTML original rather than as a wall of
DEBUG output.

### 6. Identity is derived once

`input_files` is iterated rather than branched on, so CRAM and single-end FASTQ stop
rendering an empty line. A `sample_name` is derived from the input basename and reaches
`<title>`, `<h1>`, the print header and the archive. The report states the assay
("MUC1 VNTR genotyping"), the reference assembly and the VNTR region span — none of
which appear anywhere on the page today.

`report_date` is the render time (`generate_report.py:589`), so re-running
`vntyper report` over an archived run restamps it. The context gains a distinct run
timestamp read from `pipeline_summary.json`, and the report prints both.

## Data flow and invariants

- **I1.** The render context is a **versioned public interface, not an internal detail**.
  `config.json:138` makes the template directory configurable and both renderers load
  from it (`generate_report.py:556`, `cohort_summary.py:295`), so a third-party template
  may reference exactly the keys the shipped template ignores. The six `*_color` keys and
  `igv_content` are therefore **deprecated, not deleted**: they keep being passed, are
  documented as deprecated with a removal release, and a test asserts every *new* key is
  referenced by some shipped template. Any parity test must resolve `include`/`extends`
  recursively or it will report false positives the moment the base partial lands.
- **I2.** Emphasis is derived from `screening_state`, never from `summary_text`. This
  invariant already exists and is pinned by `test_the_template_no_longer_decides_emphasis_from_the_message_text`;
  it must survive the rewrite.
- **I3.** Autoescaping stays on, and no new `|safe` interpolation is introduced beyond
  the existing fragments. `test_no_inner_or_outer_html_is_assigned_from_a_non_literal`
  and `test_no_template_concatenates_a_variable_into_parsed_markup` govern the new JS.
- **I4.** The rendered report resolves zero external hosts by default. A test greps the
  output for `http://` and `https://` in `src`/`href` attributes.
- **I5.** A row present in the results frame is present in the rendered table.
- **I6.** The report is deterministic with respect to the network. The rows it shows and
  the numbers it prints are identical with and without connectivity. I5 supplies half of
  this by removing the filter; the other half is that no displayed value is computed in
  the browser — `applyRounding` moves server-side, so a printed Depth Score cannot depend
  on whether jQuery resolved.

## Test strategy

Existing tests that must keep passing unchanged: the escaping guards, the coverage
statistic assertions, the screening emphasis tests, and the IGV fragment tests.

Existing tests that change deliberately, and must be updated in the same commit as the
code they pin:

- `test_the_confidence_column_is_colour_coded` — pins the inline `color:red` style that
  D5 removes. It becomes an assertion on the confidence pill class.
- `test_a_positive_screening_is_styled_as_a_finding` — pins `.summary-positive`. It
  becomes an assertion on the verdict banner state.
- `test_no_template_line_exceeds_the_repository_line_length` — unchanged at 120 chars,
  and it now governs materially more template source.

**A browser tier is required, and the unit tier cannot substitute for it.**
`tests/unit/test_template_escaping.py:17` says so in its own words: the unit tier has no
browser and no JavaScript engine, so its assertions are source-text tripwires. Three of
this design's central claims are invisible to it, and an earlier draft of the plan
"pinned" the first with a test that passes against the unmodified code:

| Claim | Why source assertions cannot see it |
| --- | --- |
| Flagged rows are visible | The rows are already in the server-rendered HTML (`to_html`, `generate_report.py:528`); the hiding is a client-side DataTables filter. Grepping the HTML for the flag reason passes **today**. Only a rendered, post-initialisation **visible row count** observes the defect — online and offline. |
| Print produces a correct record | A stylesheet can satisfy every substring assertion while the PDF still omits closed `<details>`, page numbers or long cells. Needs print-to-PDF with extracted-text assertions, against a named renderer and version. |
| Every control is keyboard-operable | Focusability is computed, not declared. Needs real tab traversal, including the IGV-present and IGV-failure cases. |

Two further tests must be behavioural rather than syntactic: the no-external-host check
has to inventory *every* URL-bearing construct — CSS `url()`, `@import`, `srcset`,
dynamic JS, IGV session URLs — across both templates and both asset modes, backed by
browser network interception that fails on any unapproved request; and the contrast test
must parse tokens out of the **rendered** CSS, since a test-local palette passes happily
while the template drifts.

New tests:

- A flagged variant row is present in the rendered HTML for every `*_flagged` state.
- Each of the three verdict states renders its banner, and `matched_rule == False`
  renders `indeterminate` rather than a bare negative.
- The rendered report contains no external host in `src`/`href`.
- Every token pair in the palette meets 4.5:1; computed in the test, not asserted from
  a table of numbers.
- A `@media print` block exists and the truncation clamp is unset inside it.
- Every `<details>` has a `<summary>`; no `input[type=checkbox]` is hidden with
  `display: none`.
- A CRAM-input and a single-end-FASTQ context both render a non-empty input line.
- A negative Kestrel run renders an authored empty state and not a `None`/`NaN` row.

## Rollout

Phase 1 and Phase 2 ship as one PR against one issue, as separate commits, because they
touch the same template and splitting them would leave `main` with a half-migrated
stylesheet. Phase 3 is deliberately deferred.

- **Phase 1 — clinical safety, no visual redesign.** D1, D3 (identity and print), D4,
  D7, D8, and the server-side half of D0 (move `applyRounding` out of the browser). These
  are defects in a clinical artifact and none of them requires the dependency change.
- **Phase 2 — presentation.** D2, D5, D6, D9 and the remaining half of D0: verdict
  banner, token layer, shared base template, dependency removal, vendored igv.js.

D0 is what makes the two phases one PR rather than two. Phase 1 alone leaves a report
that no longer hides rows but still renders differently offline; Phase 2 alone leaves a
self-contained report that still hides them. Neither half is shippable as a clinical
artifact on its own.
- **Phase 3 — deferred.** Dark mode, CSV export, per-variant deep links, and applying
  the token layer to the cohort report's plots.

## Compatibility and environment

- A user-supplied `report_config.json` carrying only `message` keys renders unchanged.
- `--report-assets cdn` preserves today's delivery for anyone who needs the smaller
  output directory.
- No new Python dependency. No new build step. Jinja2, pandas and the standard library
  only.
- Output-directory contract: vendored assets land under a subdirectory of the existing
  output directory and go through `contained_output_path`, so the archive-safety
  guarantees in `archive_safety.py` continue to hold.

## Section preservation

A template rewrite with no explicit preservation criterion drops sections silently. Every
one of these renders today and must still render, each with a test naming it: screening
summary, **cross-match summary** (`report_template.html:308`, tested at
`test_generate_report.py:433`), BAM-header warning, coverage basic and detailed views,
fastp metrics, Kestrel results, adVNTR results, IGV, pipeline log, report date, pipeline
version, input files.

Two behaviours are being changed rather than preserved, and must be decided explicitly
rather than dropped:

- **Search and row counts.** DataTables supplies search, paging and "showing N of M"
  today on both reports. This design refuses to reintroduce paging, because its absence
  is what makes print correct. Search and a server-rendered visible/total count are a
  separate decision — removing paging does not require removing search from screen media.
- **Display rounding.** `applyRounding` rewrites **every numeric cell in every
  initialised table** to 4 dp, not just Depth Score. Moving it server-side therefore
  changes adVNTR's displayed precision too (mean coverage, p-value). This needs explicit
  per-column server-side formatters tested against the current online rendering, or a
  documented output-format migration. It is not a one-column change.

## Unresolved questions

P1-P3 are **preconditions**: the tasks that depend on them cannot start until they are
answered by the domain owners.

1. **P1 — Intended use and labelling.** `README.md:315` states research use only. Does a
   top-level `finding` / `no-finding` / `indeterminate` label get created at all, and if
   so with what wording and whose approval? If the answer is no, the banner renders the
   configured `segments[0]` plus the three orthogonal chips and no verdict word. *Blocks
   the banner task.*
2. **P2 — Artifact shape.** Single self-contained HTML at roughly +1.3 MB per report, or
   a report plus a sibling asset directory, or CDN with SRI? This design recommends
   single-file inlining because it is the only option that keeps `--report-file` honest.
   *Blocks the asset task.*
3. **P3 — Sample identity.** "Strip suffixes from the basename" does not specify
   `sample_R1.fastq.gz` versus `sample_R2.fastq.gz`, multiple dots, lane names, or
   conflicting mate stems, and can produce misleading or duplicate identifiers. Prefer an
   explicit operator-supplied `--sample-name` with a conservative documented fallback?
   *Blocks the identity task.*
4. **P4 — Cohort scope.** Does the cohort report keep flag filtering, search and paging?
   Hiding flagged rows is defensible for triage and indefensible for a single-patient
   read, so the two reports may legitimately diverge. Note
   `tests/unit/test_cohort_summary_oracle.py:411` holds a whole-document fingerprint that
   any cohort change must deliberately re-baseline.
5. **P5 — Provenance.** The report cannot honestly print the reference assembly and VNTR
   span today: the resolved region is computed at `pipeline.py:411` and never persisted,
   and `summary.py:26` stores only start time, version, inputs and steps. Reading
   `config["default_values"]["reference_assembly"]` would mislabel any `--reference-assembly`
   override and cannot reconstruct `--custom-regions` at all. Persisting declared
   assembly, detected assembly and resolved region needs a summary schema version and an
   explicit "not recorded" fallback for legacy runs. Do that here, or print only what is
   already known?
6. **P6 — Stage failure versus true negative.** `compute_algorithm_result` returns the
   configured negative for an empty frame (`screening_summary.py:166`), and report
   generation ignores the `result_file_missing` marker `summary.py:225` records. A missing
   Kestrel step and a completed negative one are therefore indistinguishable, and the
   banner would call both "no finding". Add per-algorithm execution states
   (performed / not-performed / failed / unavailable) here, or track separately with #173?
7. **#173 seam.** A graded `confidence_tier` belongs in the same banner. Land the field as
   `None` now so #173 is purely additive later?
