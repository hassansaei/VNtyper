# Issues to file from the adversarial review of PR #199

An adversarial Codex review (gpt-5.6-sol, xhigh, five scoped lanes) of PR #199 returned 47
findings, 6 Critical. Most are being fixed on the branch. The ones below are **not**: each
either predates the branch, needs a decision that is not mine to make, or changes output in
a way that needs @hassansaei's sign-off.

Every claim here was verified against source before being written down. Where a number is
quoted, the command that produced it is given.

---

## R1 — Cohort ZIP inputs have no stable sample identity or ordering

**Labels**: bug, cohort
**Severity**: Important. Pre-existing; not introduced by #199.

`discover_sample_directories()` extracts each ZIP input into
`tempfile.mkdtemp(prefix="cohort_zip_")` and then returns `sorted(processed_dirs)`. The
sort key is the full extracted path, so it includes the randomised temporary component.
Two consequences:

1. **Ordering across multiple ZIPs is not reproducible.** Run the same cohort twice and the
   rows can come out in a different order in every HTML, CSV, TSV and JSON artefact.
2. **Worse — when `pipeline_summary.json` sits at the ZIP root, the sample identity becomes
   the random directory name.** `processed_dirs.add(temp_path)` adds the temporary root
   itself, so the reported `Sample` is a `cohort_zip_*` string that differs on every run,
   and any pseudonym derived from it differs too. **That root-level layout is what the web
   worker produces**, so this is the normal path for web cohorts, not an edge case.

`vntyper/scripts/cohort_inputs.py`, the ZIP branch of `discover_sample_directories`.

The determinism fix in #199 (`90f61fa`) is real but applies only to fixed directory
inputs, where `sorted()` on distinct full paths is a total order. The PR body has been
corrected to scope the claim accordingly.

**The decision needed**: what is the stable identity of a sample that arrives as a
root-level ZIP? The archive filename is the obvious candidate, but it is client-supplied.
A stable sort key of `(archive_path, member_path)` fixes ordering independently of the
identity question and could land first.

A characterisation test pinning the current behaviour was added in #199 so this stays
visible.

---

## R2 — Cohort pseudonyms are 20 bits and collide at realistic cohort sizes

**Labels**: bug, cohort, privacy
**Severity**: Important. Pre-existing; not introduced by #199.

A pseudonym is `<prefix>` plus the **first five hexadecimal characters** of the MD5 of the
sample name — 20 bits of entropy. `sample_mapping` is keyed on the pseudonym, so a
collision silently overwrites the earlier original name.

Verified collision:

```
$ python3 -c "import hashlib
for n in ('sample_42','sample_919'): print(n, hashlib.md5(n.encode()).hexdigest()[:5])"
sample_42 168eb
sample_919 168eb
```

Both become `anon_168eb`. Their rows are then indistinguishable by sample ID across every
export; `sample_categories()` groups them as one sample, so a cohort containing one
positive and one negative counts **one** result rather than two; and the surviving
`pseudonymization_table.tsv` entry maps the shared pseudonym to whichever original was
written last — which may be the wrong one.

Birthday probability of at least one collision: `1 - exp(-n(n-1)/2^21)`. At 1,000 samples
that is **~37.9%**. At 200 samples it is ~1.9%.

`test_different_sample_names_get_different_pseudonyms` asserted this injectivity while
checking only `"s1"` and `"s2"`. #199 renames it to what it actually checks and adds a test
pinning the collision, referencing this issue.

**The decision needed**: how many digest characters, and does changing the pseudonym format
break any downstream consumer or existing report? This changes output, so it needs
@hassansaei.

---

## R3 — The cohort report's Flag tooltip undoes the server's HTML escaping

**Labels**: bug, security, cohort
**Severity**: Important. Pre-existing; unchanged by #199.

The server escapes every `Flag` cell correctly. Then, in the browser,
`updateFlagColumn()` in `vntyper/templates/cohort_summary_template.html`:

```js
originalText = $flagCell.text().trim();          // decodes the escaping
...
$flagCell.html('<span data-toggle="tooltip" title="' + originalText + '" ...>&#10003;</span>');
```

`.text()` returns the **decoded** value, which is then concatenated into an HTML string and
reparsed by `.html()`. Server-side escaping does not survive that round trip, and the value
lands in an attribute with no attribute-context escaping.

A `Flag` value of `"></span><img src=x onerror=alert(1)><span title="` renders safely in
the initial document and executes after `updateFlagColumn()` runs. `Flag` values come from
supplied `pipeline_summary.json` files, so they are not trusted markup.

**Fix direction**: build the span with DOM APIs that keep attributes out of markup —
create the element, `.attr("title", originalText)`, `.text()` for the symbol — rather than
string concatenation into `.html()`.

The server-side escaping policy is sound and was independently confirmed; this is purely
the second-stage JavaScript.

---

## R4 — An interrupted mutation sweep can leave a live mutant in clinical source

**Labels**: bug, tooling
**Severity**: Critical for what it can produce, low probability.

`scripts/mutation_test.py` rewrites `vntyper/scripts/*.py` **in place**. Its `try/finally`
correctly restores on ordinary exceptions and on `KeyboardInterrupt`, and SIGTERM is
explicitly converted into `KeyboardInterrupt`. But no handler is installed for **SIGHUP**
or **SIGQUIT**, whose default actions terminate the interpreter without unwinding Python
`finally` blocks. SIGKILL, an interpreter crash, and host loss are uncatchable by any
handler. Restoration is itself another in-place `write_text()`, with no outer recovery if
it is interrupted.

Concretely: a long sweep run over SSH, the session drops, SIGHUP arrives during a pytest
run, and `vntyper/scripts/<target>.py` is left mutated. A later commit, build, editable
install or Docker build then consumes a deliberately introduced clinical decision defect as
ordinary source.

#199 adds a green-baseline preflight (a separate defect) but does **not** change the
in-place design.

**Fix direction**: run mutants in a disposable copy or git worktree rather than mutating
the working tree. If in-place is kept, write a durable backup and a recovery sentinel
before the first mutation, and have a supervising process own restoration. Converting
SIGHUP/SIGQUIT into an unwind helps but cannot make in-place mutation safe against SIGKILL.

---

## R5 — CRAM decoding has no reference contract

**Labels**: bug, cram
**Severity**: Important.

`process_bam_to_fastq()` sets `cram_ref_option = ""` unconditionally, even though both CRAM
command builders accept a `-T ref.fa` fragment. The CLI resolves `bwa_reference` from
configuration but uses it only for FASTQ alignment; it is never passed to CRAM decoding.
The web worker passes `--reference-assembly`, which selects coordinates, not a FASTA path.
There is no web form field for a CRAM reference and no `REF_PATH`/`REF_CACHE` configuration
anywhere in the repository.

A CRAM whose reference blocks are not embedded and whose header `UR:` names a path on the
producing institution's filesystem therefore cannot be decoded.

**Measured**, on a purpose-built reference-compressed CRAM with its reference then removed
(samtools 1.20 / htslib 1.23):

| Step | Result |
|---|---|
| `samtools quickcheck -v` — what `validate_bam_file()` runs | **exit 0** — the job is accepted |
| `samtools view -c` — counting records | **exit 0** — counting needs no reference |
| `samtools view -h` — the unmapped scan VNtyper runs | **exit 1**, slice decode failure |

So the failure is **loud but late and opaque**. Validation passes, the job is queued and
started, and it dies deep inside `process_bam_to_fastq` with
`[E::cram_decode_slice] Unable to fetch reference` in the stage log while the web surface
shows only a generic failed subprocess. It is *not* a silent wrong answer — that is the
good news and it narrows the severity — but a user gets no indication that the missing
piece is a reference they could have supplied.

Note also that htslib honours the header's `UR:` field. On the machine that produced the
CRAM the reference resolves and everything works, which is exactly why this defect
survives a sender's own testing and only appears on the receiving server.

Whether htslib's default `REF_PATH` reaches the network could not be settled here: the
synthetic CRAM's `M5` would 404 at EBI immediately, so a fast failure is consistent with
both "no attempt" and "attempt refused". A deployment behind a proxy, or one where the
`M5` does resolve at EBI, may behave differently. Settling it needs a CRAM whose `M5` is a
real assembly digest, on a host with the network blocked at the firewall.

This became user-visible when #199 routed web CRAM uploads through `--cram`.

**Cheapest useful fix**, independent of the full reference contract: run
`samtools view -h <cram> | head` (or decode one slice) during validation, so an
undecodable CRAM is rejected at submission with a message naming the reference, instead of
failing an hour into a cohort run.

**The decision needed**: resolve a full reference FASTA matching the CRAM and pass a quoted
`-T` to every decoding `samtools view`, or configure and verify an htslib reference
cache/refget before accepting CRAM at all. The shipped chromosome-1 BWA reference is not
sufficient for the whole-file unmapped scan, which can meet records from any contig.

Related: #188 (no CRAM fixture in the test cohort) blocks validating any of this.

---

## R6 — The pipeline builds a second BAM index it does not know about

**Labels**: bug, web
**Severity**: Important.

The upload endpoint and worker deliberately accept both `sample.bam.bai` and `sample.bai`,
and correctly skip the preflight `samtools index` when the alternate name was supplied. But
`process_bam_to_fastq()` reconstructs the index path as `f"{in_bam}.bai"` only. Given
`sample.bai`, it does not find it, and default non-fast processing runs `samtools index`
and creates `sample.bam.bai`.

#199 fixes the worker-side half — cleanup now removes every deterministically derived index
name. The remaining half is in the pipeline: it should resolve **both** htslib-supported
index names before creating another, or the resolved index should be passed through the
CLI boundary.

`vntyper/scripts/fastq_bam_processing.py`.

---

## R7 — `scripts/` is not measured by coverage

**Labels**: tooling, test-coverage

`[tool.coverage.run] source = ["vntyper", "docker/app"]`. The 2,637-line golden-cohort
harness and the mutation harness contribute nothing to the coverage total and no gate
notices if they are untested. Before #199 no test referenced `golden_cohort` at all:

```
$ grep -rln "golden_cohort" tests/ | wc -l
0
```

#199 adds unit tests for the harness fixes it makes, but the measurement gap remains, and
`scripts/` is also outside mypy (#204).

**The decision needed**: adding `scripts/` to `coverage.source` will move the total and
therefore the `fail_under` floor. That should be a deliberate step with its own ratchet,
not a side effect.
