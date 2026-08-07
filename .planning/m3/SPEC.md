# Milestone 3 — Web service and cohort integrity

GitHub milestone 3, *"2. Web service and cohort integrity"* — issues #216 #207 #206 #205 #201
#167 #162. All seven open, none closed. Audited against `main` @ `b46da80`.

## Exit criteria (the three claims this milestone must be able to prove)

Restated after the adversarial review — the original three overclaimed, and `.planning/m3/REVIEW.md`
records what was rejected and why. What is *not* proven is stated here rather than left implied.

| # | Claim | Proof (Phase 5) |
| --- | --- | --- |
| E1 | On a web job whose cleanup succeeds, no file and no directory remains under the input tree — and cleanup no longer fails on the normal path, so `tasks.py:409`'s `os.rmdir`, dead code in production until now, runs | Run `vntyper pipeline` against a `chmod a-w` input mount; assert it completes and the input directory is byte-unchanged. **Not claimed:** that cleanup always succeeds — `tasks.py:351` marks the job completed before the `finally:` block, and removal failures are logged and swallowed by design (`test_tasks.py:715`, `:752`) |
| E2 | The machine-readable cohort exports (CSV/TSV/JSON) and `pseudonymization_table.tsv` are byte-identical across two runs; the HTML is identical after the oracle's existing normalisation | Two `aggregate_cohort` runs over the same ZIP inputs in two processes; `diff` every export. **Not claimed:** raw HTML byte-equality — it carries a report timestamp and Plotly UUIDs, which is why `_skeleton()` exists |
| E3 | No template builds markup from a value it did not author: every sample-derived value reaches the DOM through `.attr()`, `.text()` or `textContent`, and everything entering a `<script>` block is `json.dumps` output | A source-text tripwire over both templates, plus a Python-side test that a malicious `Flag` in a stored `pipeline_summary.json` is server-escaped in the rendered HTML. **Not claimed:** runtime proof — this repo has no browser test tier and `AGENTS.md` forbids one in the unit tier, so the residual risk is a future sink the tripwire's pattern does not match |

---

## Workstream A — XSS (#207, #216)

### A1. #207 — Flag tooltip undoes the server's escaping

**Root cause.** `updateFlagColumn()` exists **twice**, verbatim, in
`vntyper/templates/cohort_summary_template.html:135-162` and
`vntyper/templates/report_template.html:406-429`. Both do:

```js
originalText = $flagCell.text().trim();                       // decodes server escaping
$flagCell.html('<span data-toggle="tooltip" title="' + originalText + '" …>&#10003;</span>');
```

Jinja2 autoescapes the `Flag` cell server-side. `.text()` returns the **decoded** value; the
string concatenation puts it into an attribute with no attribute-context escaping; `.html()`
reparses it. A `Flag` of `"></span><img src=x onerror=alert(1)><span title="` executes on the
first DataTables `preDrawCallback`.

**Why the value is untrusted — the trust boundary is the stored artefact, not the rule engine.**
On a normal pipeline run `Flag` is safe by construction: `flagging.py:135-147` sets
`df_copy["Flag"] = ", ".join(flag_list)` where `flag_list` holds the **keys** of
`flagging_rules` — config-declared identifiers such as `False_Positive_4bp_Insertion`. No
DataFrame value is interpolated, so there is no read-derived payload path there. (An earlier
draft of this spec claimed otherwise; the review caught it.)

The untrusted path is the one #207 names: *"`Flag` values come from supplied
`pipeline_summary.json` files, so they are not trusted markup."* `vntyper cohort` consumes ZIPs
uploaded by the client and `vntyper report` regenerates over an output directory the user may
have received from elsewhere. The tests therefore plant a malicious value in a **stored
artefact**, which is the real boundary.

**Invariant.** No sample-derived string is ever concatenated into markup that is subsequently
parsed. Attribute values are set through DOM APIs that take a value, not a fragment.

**Fix.** In both templates, replace the two `.html('<span …>')` calls with a jQuery DOM build:

```js
var $mark = $('<span>')
    .attr('data-toggle', 'tooltip')
    .attr('title', originalText)
    .css({ color: isClean ? 'green' : 'red', 'font-size': '16px' })
    .text(isClean ? '✓' : '✖');
$flagCell.empty().append($mark);
```

`.attr()` and `.text()` set the property/text node directly — no HTML parse, so no escape to
undo. The `✓`/`✖` escapes replace the `&#10003;`/`&#10006;` entities, which only
existed because the value was going through an HTML parse.

**Do not change** the `data-original` mechanism: it is already `.attr()`-set, it is what the
DataTables search filter at `cohort_summary_template.html:103` reads, and it is what makes the
function idempotent across redraws (after draw 1 the cell's text is the tick, not the flag).

**Blast radius.** Both templates; the tooltip title and the tick/cross glyph are the only
user-visible surface. The `test_cohort_summary_oracle.py` fingerprint substitutes `<SCRIPT-BODY>`
for every script body before hashing, so this change alone must not move the fingerprint —
re-run it and confirm.

**The test that would have caught it.** A source-text tripwire over both templates
(`tests/unit/test_template_escaping.py`, new): assert neither template contains a `.html(` call
whose argument concatenates a variable into an attribute, and assert both contain
`.attr('title', originalText)`. This repo already uses source-text tripwires
(`tests/unit/test_kestrel_filtering.py` reads `kestrel_genotyping.py` as text). The unit tier
cannot execute JavaScript, so the tripwire is the honest gate — it is stated as a tripwire, not
as a behavioural test.

### A2. #216 — igv-reports fragments inlined into `<script>`

**Root cause.** `report_template.html:469-470`:

```html
const tableJson = {{ table_json | safe }};
const sessionDictionary = {{ session_dictionary | safe }};
```

Both come from `extract_igv_fragments` (`report_formatting.py:441-442`) via
`extract_line_after`, which lifts the raw remainder of the marker's line out of the
igv-reports page. The only guard is `js_object_literal` (`report_formatting.py:403-419`), which
returns the fragment unchanged if it is non-empty. That is a *syntax* guard against
`const tableJson = ;`, not a safety guard: it does not parse the fragment and does not
neutralise `</script>`.

**Why the data is untrusted.** `tableJson` is the igv-reports variant table, built from
VNtyper's own BED and VCF, whose `REF`, `ALT` and `Motif_sequence` values are derived from the
sample's reads. `vntyper report` also regenerates a report over an output directory the user may
have received from elsewhere.

**Invariant.** Nothing reaches a `<script>` block that VNtyper did not itself serialise. The
value emitted is always the output of `json.dumps`, escaped for script context.

**Fix.** Add a pure function to `report_formatting.py` and use it at the two `generate_report.py`
call sites (`generate_report.py:587-588`):

```python
def js_json_literal(fragment: str, fallback: str) -> str:
    """Re-serialise an extracted fragment as a script-safe JSON literal."""
```

Behaviour:
1. Strip surrounding whitespace and a single trailing `;`. **Verified against the installed
   igv-reports**: `templates/variant_template.html:155-156` is `const tableJson = "@TABLE_JSON@"`
   with *no* terminator, and `report.py:178-183` substitutes the placeholder including its
   quotes, so `extract_line_after` returns the bare literal. The strip is therefore defensive
   against an igv-reports version that does emit one, not a correction of today's output — say so
   in the docstring rather than claiming the terminator is present.
2. `json.loads`. **Verified the fragments are valid JSON**: `bedtable.BedTable.to_JSON`
   (`bedtable.py:34-37`) is a `json.dumps`, and the session dictionary is built the same way. On
   `ValueError`, log at `warning` and return `fallback` — a fragment that is not JSON is not
   something to emit.
3. `json.dumps(value, sort_keys=True, separators=(",", ":"))` — deterministic, so two runs over
   the same IGV page emit byte-identical script (feeds E2).
4. Escape for script context: `</` → `<\/`, `U+2028` → `U+2028`, `U+2029` → `U+2029`.
   The first closes the `</script>` hole; the latter two are line terminators in JS but not in
   JSON strings.

`js_object_literal` is **replaced**, not kept alongside — one function, one purpose. Its two
call sites and `tests/unit/test_report_formatting.py` move with it.

**Also in scope (named in #216).** `report_template.html:515` `cell.innerHTML = headers[j]` and
`:546` `cell.innerHTML = rowData[j]` become `textContent`.

Verified against the installed igv-reports that this is not a display regression. VNtyper calls
`create_report <bed> --flanking N --fasta X --tracks <vcf> <bam>` (`generate_report.py:124-159`),
so the sites table is a `BedTable`. `generate_bed_file` (`kestrel_genotyping.py:682-686`) writes
**three** columns — `Motif_fasta`, `POS_fasta - 1`, `POS_fasta` — so `has_name` is false, the
headers are `["unique_id", "Chrom", "Start", "End"]`, and the rendered cells are the motif record
name and two integers. None of those is HTML-escaped by igv-reports, and the motif record name is
sample-derived, so `innerHTML` is the sink and `textContent` is the fix.

The one column igv-reports *does* pre-escape is `Name` (`bedtable.py:32`,
`html.escape(feature.name)`), which VNtyper never emits. If a future BED gains a name column it
would render its entities literally under `textContent` — note that in the code comment so the
next reader is not surprised, but do not add an unescape step for a column that does not exist.

**Blast radius.** `report_formatting.py`, `generate_report.py` (two lines),
`report_template.html`. A sample with no IGV report already falls through to
`EMPTY_TABLE_JSON`/`EMPTY_SESSION_DICTIONARY`, and step 2 routes an unparseable fragment to the
same place — so the "no IGV" branch at `report_template.html:479` keeps working unchanged.

**The test that would have caught it.** `js_json_literal("{}</script>…", …)` — assert the
returned string contains no `</`. Plus a round-trip test that a well-formed fragment survives,
and a fallback test that a malformed one does not reach the output.

---

## Workstream B — the input directory is written to (#201, #162)

### Root cause

`vntyper/scripts/utils.py:361-363`:

```python
command = f"samtools quickcheck -v {quote_path(file_path)}"
log_file = f"{file_path}.quickcheck.log"
success = run_command(command, log_file, critical=False, cwd=cwd)
```

The log path is derived from the **input**. `run_command` opens it unconditionally
(`utils.py:51`) *before* running anything, so on a read-only input mount this is an unhandled
`PermissionError` raised before quickcheck executes. `validate_bam_file` is called at
`pipeline.py:219` and `:221`, so **a read-only input mount fails every BAM and CRAM run at
HEAD** (measured by the issue author).

On the web service the input directory is per-job on the shared upload volume.
`docker/app/tasks.py:403-410` removes the alignment and its index, then refuses to `os.rmdir` a
non-empty directory — and it is never empty, because the quickcheck log is still there. So the
`os.rmdir` is dead code in production, and the "still holds files and was left in place" warning
fires on 100% of jobs, which is exactly when a warning stops being able to report the
genuinely-unexpected leftover it was written for.

### Invariant

**No VNtyper process writes into the input directory.** An artefact of a run belongs with the
run's other artefacts, in the output directory.

### Fix

`validate_bam_file(file_path, cwd=None, log_dir=None)`. With `log_dir` set, the log is
`Path(log_dir) / f"{Path(file_path).name}.quickcheck.log"`. `log_dir=None` keeps today's
behaviour — this is the issue author's explicit recommendation ("the parameter should default to
today's behaviour and be passed explicitly by the pipeline; that keeps this a contained change
rather than a breaking one") and is followed rather than re-litigated. `pipeline.py` passes
`log_dir=output_dir` at both call sites.

The log is **kept**, not deleted: a job that dies on a corrupt upload must still leave the
diagnostic that explains why. Adding the log to the worker's removal tuple was considered and
rejected in the issue for three stated reasons; that decision stands.

### The second write — #210, pulled into scope (Fork B-1, resolved)

#162 is the umbrella and names **two** production writes:

1. the quickcheck log — this is #201;
2. a redundant `samtools index <in_bam>` writing `<in_bam>.bai` beside the input
   (`fastq_bam_processing.py:160-170`, non-fast **BAM only** — the CRAM branch writes no index
   beside its input) — this is #210, filed under milestone 4.

Fixing only (1) makes a read-only input mount work for every CRAM run and every `--fast-mode`
BAM run — which is what the web service submits — but a non-fast BAM run would still fail with
`RuntimeError("Indexing BAM file failed.")`. **#210 is pulled into this PR so that #162's own
acceptance criterion — no VNtyper process writes into the input directory — becomes true and
testable now.**

#210's fix has two halves, and both are needed:

*Resolution order.* `process_bam_to_fastq` reconstructs the index as `f"{in_bam}.bai"` only. The
upload endpoint and worker deliberately accept `sample.bai` as well as `sample.bam.bai`, so given
`sample.bai` the pipeline does not find it and builds a second index nothing else knows about.
Follow htslib's own order: `<file>.bai`, then `<stem>.bai`.

*Destination.* When neither exists, the index is built **into the output directory**, not beside
the input: `samtools index -o <output>/<name>_input.bam.bai <in_bam>` (`-o` requires samtools
≥1.15; `conda/environment_vntyper.yml` pins 1.20). `extract_unmapped_reads_from_offset` already
takes `bai_file` as an explicit argument, so the resolved path threads straight through.

### Blast radius

`utils.py` (one signature, one path), `pipeline.py` (two call sites),
`fastq_bam_processing.py:160-177`, `command_builders.py` (`build_samtools_index_command` gains an
optional output path). `docker/app/tasks.py` is **not edited** — its `os.rmdir` simply becomes
reachable. `tests/data/*.quickcheck.log` housekeeping is explicitly out of scope per the issue.

### The tests that would have caught it

- `run_command` with a log path inside a `chmod 0o555` directory raises `PermissionError` —
  pins the mechanism, no mock.
- `validate_bam_file(bam_in_readonly_dir, log_dir=tmp_path/"out")` completes and the log lands
  under `out/`, with a `run_command` spy asserting the exact path and that it is **not** under
  the input directory.
- `pipeline.py` passes `log_dir` — asserted through the existing `tests/support/pipeline_harness.py`
  spy, which already stubs `validate_bam_file`.
- `tests/unit/web/test_tasks.py`: an input directory holding only the alignment and its index is
  `os.rmdir`-ed and logs no warning; one holding an unexpected extra file still warns and is
  still left in place. The guard must keep working — it must only stop firing on the normal path.

---

## Workstream C — cohort identity and order (#206, #205)

Both issues carry `needs-decision`. The mechanisms below are settled; the two open questions are
Forks C-1 and C-2.

### C1. #206 — pseudonyms are 20 bits and collide

**Root cause.** `cohort_inputs.py:176`:

```python
hash_suffix = hashlib.md5(original_sample.encode()).hexdigest()[:5]
return f"{prefix}{hash_suffix}"
```

Five hex characters is 20 bits. `sample_mapping` is keyed on the pseudonym
(`cohort_summary.py:396`), so a collision silently overwrites the earlier original name. Verified
first collision in `sample_0..sample_19999`: `sample_42` and `sample_919` both → `168eb`.
Birthday probability of at least one collision is `1 - exp(-n(n-1)/2^21)`: **~37.9% at 1,000
samples**, ~1.9% at 200.

**Consequences.** Two patients' rows become indistinguishable across every export;
`sample_categories()` counts them as one sample, so a cohort with one positive and one negative
reports **one** result; `pseudonymization_table.tsv` maps the shared pseudonym to whichever
original was written last.

**Invariant.** The pseudonym map is **injective**. Two distinct originals never share a
pseudonym — and if the digest ever says otherwise, the run says so rather than merging two
patients' data.

**Fix, in two parts.**

*Part 1 — widen the digest, config-driven (Fork C-1, resolved: **SHA-256, 12 hex characters**).*
48 bits — collision probability ~1.8e-10 at 10,000 samples. MD5 is dropped as well as widened:
`hashlib.md5()` raises on a FIPS-enabled build unless called with `usedforsecurity=False`, and
since the pseudonym format changes either way this is the moment to stop depending on it.

Algorithm and length move into `vntyper/config.json` under a new `cohort.pseudonym` section, read
with `.get()` chains against module-level defaults (AGENTS.md trap 2: `--config-path` replaces the
whole config, so a missing key must not `KeyError`). No length or algorithm name is written
inline in Python. An algorithm not in `hashlib.algorithms_available` is refused by name rather
than falling back silently.

*Part 2 — detect, never merge (Fork C-3, resolved: **abort**).* Regardless of length,
`aggregate_cohort` refuses to map two different originals to one pseudonym: `logger.error` then
`raise ValueError` naming both colliding originals, per the repo's error convention. A cohort
that silently merges two patients' genotypes is worse than one that refuses to run. At 48 bits
this is a tripwire rather than an operational risk.

### C2. #205 — ZIP inputs have no stable identity or reproducible order

**Root cause.** `discover_sample_directories` (`cohort_inputs.py:110-140`) extracts each ZIP into
`tempfile.mkdtemp(prefix="cohort_zip_")` and returns `sorted(processed_dirs)`. The sort key is
the **full extracted path**, which contains the randomised temporary component. Two failures
follow:

1. **Order is not reproducible across runs.** Rows come out in a different order in every HTML,
   CSV, TSV and JSON artefact. (The #199 determinism fix is real but applies only to fixed
   directory inputs, where `sorted()` over distinct full paths is a total order.)
2. **A ZIP with `pipeline_summary.json` at its root gets the random directory name as its sample
   identity.** `processed_dirs.add(temp_path)` adds the temporary root itself
   (`cohort_inputs.py:119`), and `cohort_summary.py:393` takes `Path(sample_dir).name` as the
   sample. The reported `Sample` is therefore a `cohort_zip_*` string that differs every run, and
   the pseudonym derived from it differs too. **That root-level layout is what the web worker
   produces**, so this is the normal path for web cohorts.

**Invariant.** A sample's reported identity and its position in the output depend only on the
inputs, never on where they were unpacked.

**Fix — ordering (uncontroversial, lands independently).** `discover_sample_directories` returns
directories ordered by an **origin key** rather than by the extracted path: for each discovered
directory, `(index of the input path on the command line, path relative to that input's root)`.
For a directory input the relative path is the sample directory's path under it; for a ZIP it is
the member path inside the archive. The de-duplicating set stays (it is what makes a directory
reached twice appear once); only the sort key changes. `sorted()` on that key is a total order
that contains no temporary component.

**Fix — identity (Fork C-2, resolved).** For a ZIP whose `pipeline_summary.json` is at the root,
the identity is **the stem of the input file the run itself recorded**, read from that summary's
`input_files` (`bam`, else `cram`, else `fastq1`) — so `patient1.bam` yields `patient1`. It is
the run's own record of what it processed and is the string a clinician recognises. When
`input_files` is absent or empty (older summaries), fall back to the **archive filename stem**.

The identity is threaded as an explicit value out of `discover_sample_directories` alongside each
directory, **not** re-derived from `Path(sample_dir).name` at `cohort_summary.py:393`. Samples
found in subdirectories keep their directory name as today — only the root-level case, where the
directory name is the random temp root, changes.

### Blast radius

`cohort_inputs.py` (return type widens from `list[Path]` to carry the identity — every caller is
in `cohort_summary.py`), `cohort_summary.py:392-416`, `vntyper/config.json`. The characterisation
test #199 added to pin the current behaviour is expected to fail and is rewritten to pin the new
behaviour, citing this milestone.

### The tests that would have caught it

- `pseudonymized_sample_name` is injective over the first 20,000 `sample_N` names (the exact
  probe that found the bug), at the configured length.
- Two `aggregate_cohort` runs in **two separate processes** over the same ZIP inputs produce
  byte-identical exports. Cross-process matters: `PYTHONHASHSEED` randomisation is the mechanism
  the existing `test_cohort_inputs.py` cross-process test was written for.
- A ZIP whose `pipeline_summary.json` is at the root reports the same `Sample` on two runs, and
  that string contains no `cohort_zip_`.
- Two originals colliding on the digest raise, naming both.

---

## Workstream D — `vntyper report` passes no VCF (#167)

**Root cause.** `cli_report.py:137` hardcodes `vcf_file=None`. `generate_report.py:138-144` only
adds the IGV variant track when `vcf_file` is set and exists, so a regenerated report never has
one — even though the run being re-reported wrote the VCF.

**Already fixed, do not redo.** The signature mismatch that the issue body opens with
(`input_files`, `pipeline_version`, `mean_vntr_coverage` passed but not accepted) was closed by
#179; `cli_report.py` today passes exactly what `generate_summary_report` accepts, and
`log_file` is resolved properly by `resolve_log_file`. **The only live defect is the VCF.**

**Invariant.** `vntyper report -o <run>` produces the same IGV panel the run itself produced.

**Fix.** `artifact_names.select_best_vcf_file()` already exists and is exactly what
`pipeline.py` uses to pick `output_indel.vcf.gz` over `output_indel.vcf`. Reuse it — do not write
a second resolver.

1. Add `--vcf-file` to the `report` subparser (`cli_parser.py:201-207`, beside `--bed-file` and
   `--bam-file`), so the user can set it. The issue names "cannot be set by the user" as half the
   defect.
2. When unset, resolve `select_best_vcf_file(<dir>/"kestrel")`, where `<dir>` is `--input-dir` if
   given, else `--output-dir`. `--output-dir` is documented as *"Directory containing pipeline
   results"*, and `resolve_log_file` already reads `<output-dir>/pipeline.log` on that basis.
   This is deliberately **more** permissive than the `--bam-file`/`--bed-file` resolution, which
   consults `--input-dir` only; the asymmetry is documented in the handler rather than silently
   introduced. Widening bam/bed to match is **out of scope** — it would change the behaviour of
   an argument this issue does not name.
3. A missing VCF stays a warning, never an error: `select_best_vcf_file` returns `None` and
   `generate_report` skips the track. bcftools is optional.

**Blast radius.** `cli_report.py`, `cli_parser.py`. `tests/unit/test_cli_parser_contract.py` pins
the argument set and will need the new flag.

**The test that would have caught it.** The existing signature-binding spy on
`generate_summary_report` (AGENTS.md trap 11) asserts `vcf_file` is the resolved path, not
`None`, for: explicit `--vcf-file`; `--input-dir` with a `kestrel/output_indel.vcf.gz`;
`--output-dir` only; and `None` when no VCF exists anywhere.

---

## File ownership — the four sets are disjoint

| WS | Source | Tests |
| --- | --- | --- |
| **A** | `vntyper/templates/report_template.html`, `vntyper/templates/cohort_summary_template.html`, `vntyper/scripts/report_formatting.py`, `vntyper/scripts/generate_report.py` | `tests/unit/test_report_formatting.py`, `tests/unit/test_generate_report.py`, `tests/unit/test_template_escaping.py` *(new)* |
| **B** | `vntyper/scripts/utils.py`, `vntyper/scripts/pipeline.py`, `vntyper/scripts/fastq_bam_processing.py`, `vntyper/scripts/command_builders.py` | `tests/unit/test_utils.py`, `tests/unit/test_run_command_contract.py`, `tests/unit/test_command_builders.py`, `tests/unit/test_fastq_bam_command_wiring.py`, `tests/unit/web/test_tasks.py`, `tests/unit/test_input_tree_is_never_written.py` *(new)* |
| **C** | `vntyper/scripts/cohort_inputs.py`, `vntyper/scripts/cohort_summary.py`, `vntyper/config.json` | `tests/unit/test_cohort_inputs.py`, `tests/unit/test_cohort_summary_oracle.py`, `tests/unit/test_cohort_identity.py` *(new)* |
| **D** | `vntyper/scripts/cli_report.py`, `vntyper/scripts/cli_parser.py` | `tests/unit/test_cli_report.py`, `tests/unit/test_cli_parser_contract.py` |

Pairwise intersections are all empty. The integrator alone owns `vntyper/version.py`,
`CITATION.cff`, `docs/about/changelog.md`, `AGENTS.md` and `pyproject.toml` — no subagent edits
them.

Two adjacencies worth naming, neither an overlap:

- A and D both reach `generate_summary_report`, but from different sides: A edits its two
  template-context lines (`generate_report.py:587-588`), D edits its caller (`cli_report.py`).
- B and C both touch code `pipeline.py` calls, but B edits `pipeline.py` and C never does.

---

## Forks — all resolved 2026-08-07

| Fork | Question | Decision |
| --- | --- | --- |
| **B-1** | Does #162 close in this PR? | **Yes** — pull #210 in. #162, #201 and #210 all close. |
| **C-1** | Pseudonym digest width and algorithm | **SHA-256, 12 hex characters**, config-driven. |
| **C-2** | Identity of a root-level ZIP sample | **`pipeline_summary.json`'s own input-file stem**, falling back to the archive filename stem. |
| **C-3** | Behaviour on a pseudonym collision | **Abort**, naming both originals. |

## Issues closed by this PR

#216, #207, #206, #205, #201, #167, #162 (the whole milestone) and #210 (milestone 4, pulled in
by decision B-1).

## Deliberately out of scope

- The committed `tests/data/*.quickcheck.log` files (#201 names this as separate housekeeping
  against managed test data).
- Widening `--bam-file`/`--bed-file` resolution in `vntyper report` to consult `--output-dir`
  (#167 does not name it; changing it would alter an argument's behaviour unasked).
- The remaining five issues in milestone 4 beyond #210.
