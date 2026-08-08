# Milestone 4 — CRAM and input robustness

**Status:** spec, awaiting adversarial review
**Branch:** `fix/milestone-4-cram-input-robustness`
**Closes:** #213, #225, #209, #178, #165, #161
**Date:** 2026-08-08

---

## 1. The exit bar

A CRAM run either succeeds, or fails with a message naming exactly what is missing,
**before any stage does work**. No input is ever silently discarded. Every behaviour this
spec adds is config-driven; no threshold, contig name, reference path or scan strategy is
written inline in Python.

"Before any stage does work" is deliberately not "at submission", because the two entry
points differ and the difference is not something this milestone can paper over:

* **CLI** — the preflight is the first thing `run_pipeline` does, so a bad input fails in
  under a second with the message on stderr. This is "at submission" in every sense a CLI
  user means.
* **Web** — `docker/app/main.py` stores the upload, enqueues a Celery task and returns
  *Job submitted* before any worker looks at the file (`main.py:463`, `:494`, `:547`), so
  no check VNtyper adds inside `run_pipeline` can make the HTTP response name the problem.
  The honest scope is therefore: the preflight runs first in the worker, and the
  **message it raises is surfaced in the job status** instead of the generic text
  `main.py:605` substitutes today.

  "Surfaced" needs a mechanism, because there is currently no path for one.
  `run_pipeline` converts every exception into `SystemExit(1)` (`pipeline.py:721`), the
  worker runs the CLI with `subprocess.run` and does not capture stderr (`tasks.py:317`),
  and the status endpoint substitutes generic text **deliberately** — it is
  unauthenticated, and raw exception text would leak absolute paths and full samtools
  command lines. So the transport is neither "print the exception" nor "return it to the
  client":

  1. The preflight writes `preflight_error.json` into the run's output directory before
     raising — `{"code": ..., "message": ..., "candidates": [...]}` — built from the same
     `alignment_contract` helpers that compose the CLI message, and containing only values
     VNtyper itself composed.
  2. The worker reads that file when the pipeline exits non-zero and stores `code` and
     `message` on the job's Redis hash.
  3. The status endpoint returns that stored message when present, and its existing
     generic text otherwise. It is a curated string rather than an exception, so the
     endpoint's no-leak property is preserved rather than traded away.

  A-WEB-1 tests exactly this, including that no absolute path from the worker's filesystem
  appears in the response. Preflighting inside the HTTP request remains out of scope (§9).

## 2. The one root cause

The six issues are usually described as six bugs. Five of them are one bug.

**VNtyper has no input contract.** It discovers what an alignment needs — an index, a
reference, a read layout, a contig naming convention — *inside stage code, at the moment
it needs it*, and delegates the resolution of the first two to htslib's ambient
environment.

Ambient resolution has three properties that together produce every symptom in this
milestone:

1. **It is late.** The need is discovered after the run has been accepted, queued and
   started, so the user learns about it from a stage log rather than from the submission.
2. **It is opaque.** htslib reports what it could not fetch, not what the operator should
   have supplied.
3. **It is unbounded in time.** Reference resolution can block on a network endpoint with
   no timeout and no output. This is the difference between #209 (a loud, late error) and
   #178 (a hang), and it is why #178 has resisted diagnosis.

`#213` is genuinely separate: it is a missing flag that makes all three invisible. It
lands first for that reason.

## 3. Measured evidence

All measurements on **samtools 1.20 with htslib 1.23**, which is what
`conda/environment_vntyper.yml`'s `samtools=1.20` actually resolves to here:

```
$ conda list | grep -E "^(samtools|htslib)"
htslib    1.23   h566b1c6_0   bioconda
samtools  1.20   h50ea8bc_1   bioconda
```

`v2.0.3` pinned the same `samtools=1.20`, but conda resolves htslib independently and at
a different date, so **the image #178 was reported against may have carried a different
htslib than the one measured here.** That is a limitation of these measurements and it is
also the point of §3.4: the behaviour VNtyper inherits from "leave it to htslib" is not
pinned by anything in this repository.

Fixture: a purpose-built reference-compressed CRAM (10 kb single contig, 50 reads,
`M5` present, `UR:` present) with its reference then removed.

### 3.1 The #209 ladder — validation passes, decoding fails

| Step | Exit | Note |
| --- | --- | --- |
| `samtools quickcheck -v out.cram` | **0** | this is what `validate_bam_file` runs — the job is accepted |
| `samtools view -c out.cram` | **0** | counting reads needs no reference |
| `samtools view -h out.cram` | **1** | `[E::cram_decode_slice] Unable to fetch reference chr1:100-5049` |
| `samtools view -P -b out.cram chr1:100-500` | **1** | the region slice fails the same way |

So a CRAM that cannot be decoded passes validation. Confirms #209 as filed.

### 3.2 Ambient reference resolution can block forever — the #178 mechanism

With `REF_PATH` pointed at an endpoint that accepts TCP connections and never responds
(a firewalled or overloaded refget server, or a proxy that drops):

```
REF_PATH="http://127.0.0.1:18081/%s" timeout 25 samtools view -h out.cram > /dev/null
exit=124 elapsed=25s
```

**Exit 124 means the timeout fired.** samtools produced no output, no error and did not
exit.

**This is a proven hang mechanism on this stack. It is not proof that it is #178's
mechanism**, and the spec previously overstated that. `run_pipeline` calls
`create_output_directories` at `pipeline.py:161`, before validation and before any CRAM is
touched, and `cli.py:139` creates the log directory earlier still — so a hang anywhere
inside `process_bam_to_fastq` would leave an output directory and a `pipeline.log` behind.
The reporter states that *no* output directory was generated, which points at something
before that line rather than at the reference fetch.

What can be claimed: #213 is why the report contains no stage log to distinguish these,
this milestone removes the reference-fetch hang and the whole-file scan as *possible*
causes, and once #213 lands the next report of this shape arrives with the evidence needed
to close it. A-178-4 asks the reporter for `docker logs` from a post-#213 image.

The control confirms this is specific to a *network* REF_PATH, not to a missing reference:

```
REF_PATH="/nonexistent/%2s/%2s/%s" timeout 25 samtools view -h out.cram
exit=1   # fails immediately, naming the reference it could not open
```

A local-only `REF_PATH` fails fast. A network `REF_PATH` blocks indefinitely.

### 3.3 `-T` makes ambient resolution unreachable

With the same blackhole `REF_PATH` still set, but an explicit `-T`:

```
REF_PATH="http://127.0.0.1:18082/%s" timeout 25 samtools view -h -T ref.fa out.cram | wc -l
54
exit=0 elapsed=0s
```

**Passing `-T` short-circuits REF_PATH entirely.** This is the load-bearing fact of the
whole milestone: supplying the reference does not merely make the common case work, it
removes the code path on which the hang lives.

### 3.4 htslib's ambient default is version-dependent

```
strings libhts.so.1.23 | grep -iE "ena/cram|ebi.ac.uk"   # no matches
```

htslib 1.23 no longer carries the EBI reference server URL that older releases used as
the `REF_PATH` default. `strace -f -e trace=connect,socket` on a decode with `REF_PATH`
unset and no `UR:` recorded no network syscalls.

**This is a reason to pin, not a reason to relax.** The behaviour VNtyper gets from
"leave it to htslib" changed between releases without VNtyper changing, which is exactly
the property an input contract exists to remove. A user on an older image, or one with
`REF_PATH` set in their environment or by their site's module system, still reaches 3.2.

### 3.5 The unmapped scan is O(whole file) for CRAM and O(index) for BAM

`build_cram_unmapped_filter_command` streams the entire CRAM through `samtools view -h`
(as SAM **text**) and filters flag 12 downstream. The BAM path does not: it reads the
BAI's last chunk offset and seeks straight to the unplaced block.

The index-based equivalent for CRAM retrieves the identical read set:

| Method | Reads recovered |
| --- | --- |
| `samtools view -h big.cram \| samtools view -b -f 12` (today) | 4000 |
| `samtools view -b -f 12 big.cram '*'` (indexed) | 4000 |

On a 30x WGS CRAM the first decompresses and re-encodes the whole genome as text; the
second seeks. #178 was run under `--platform linux/amd64`, i.e. very likely qemu
emulation on an arm64 host, where the first is measured in hours. Both mechanisms — the
blocking reference fetch and the whole-file scan — produce "hangs indefinitely", and
#213 guarantees neither leaves a trace.

### 3.6 A chr1-only reference decodes the chr1 slice

VNtyper ships `reference/alignment/chr1.hg19.fa` and `chr1.hg38.fa` (config key
`bwa_reference_<assembly>`) and `install-references` never fetches a whole genome. On a
two-contig CRAM, slicing a chr1 region with `-T chr1only.fa` and `REF_PATH` pinned to a
nonexistent local path produced results identical to `-T full.fa`: exit 0, 135 reads, all
on chr1.

That result holds only for a slice **without** `-P`. See §3.7, which overturns the
optimistic reading of it.

### 3.7 `-P` fetches cross-contig mates, so a chr1-only reference is not sufficient

Re-run on a correctly coordinate-sorted two-contig CRAM containing one deliberate
cross-contig pair (`x1` on chr1:5000, mate on chr2:7000), with `REF_PATH` pinned to a
nonexistent local path and the header's `UR:` target removed so no fallback can
contaminate the result:

| Query | Reference | Exit | Reads | Contigs returned |
| --- | --- | --- | --- | --- |
| `view -P … chr1:4900-5100` | full | 0 | 12 | chr1, **chr2** |
| `view -P … chr1:4900-5100` | chr1-only | **1** | 11 (partial, then aborts) | chr1 |
| `view … chr1:4900-5100` (no `-P`) | chr1-only | 0 | 9 | chr1 |

`[E::cram_decode_slice] Unable to fetch reference chr2:1-12070`

**`samtools view -P` — which `build_samtools_slice_command` already passes — chases mates
onto other contigs and therefore needs those contigs' reference.** A chr1-only FASTA is
sufficient only for a CRAM whose reads in the MUC1 region have no cross-contig mate, which
is a property of the *sample*, not of the pipeline. §3.6's conclusion is withdrawn.

Two consequences, both load-bearing:

1. **The probe must have the same shape as the real slice.** A probe that omits `-P`
   passes with the chr1-only reference and the real slice then fails mid-run — the exact
   "fix that masks rather than removes a failure mode" this milestone exists to avoid. The
   probe therefore uses the same flags and the same region as the slice it authorises.
2. **The shipped chr1 FASTA is a last-resort candidate, not the answer.** It is kept in
   the order because it succeeds for many samples, but the run logs which candidate won,
   and warns when it is a reference that does not cover every contig in the CRAM header.

### 3.8 htslib will not silently return wrong bases

Three failure shapes were measured, and none of them is silent:

| Situation | Result |
| --- | --- |
| Reference **absent** entirely | exit 1, `Unable to fetch reference`, 0 reads |
| Reference present but **contig missing** from it | exit 1, `Unable to fetch reference <contig>` |
| Reference present, contig present, **wrong sequence** (same name, same length) | exit 1, `MD5 checksum reference mismatch`, prints CRAM/Ref/@SQ digests, 0 reads |

htslib verifies the per-slice `M5` before it emits anything, so a mismatched reference
cannot produce wrong bases. **This is what makes "prove it by decoding" a sound design**:
a probe that exits 0 has not merely produced output, it has produced output htslib
checksum-verified against the CRAM's own digest.

The partial-output case is worth naming: the failing `-P` run wrote 11 of 12 reads before
aborting. VNtyper is safe there only because `run_command(..., critical=True)` raises on
the non-zero exit and the partial file is never consumed. That is an existing invariant
this milestone must not weaken.


### 3.9 An index outside the alignment's directory is invisible without `-X`

The obvious reading of #225 — "build the missing index into the output directory" — does
not work on its own. Measured:

```
samtools index -o out/s.cram.crai in/s.cram          # builds fine
samtools view -c -T full.fa in/s.cram chr1:4900-5100
  [E::cram_index_load] Could not retrieve index file for 'in/s.cram'
  samtools view: Random alignment retrieval only works for indexed ... files.
samtools view -c -T full.fa -X in/s.cram out/s.cram.crai chr1:4900-5100
  9
```

Same result for BAM with a relocated `.bai`. **Putting the index path in a plan object
does not make htslib find it.** Either every random-retrieval command carries
`-X <alignment> <index>`, or the alignment and its index are made co-located by some
means that does not write to the input tree. §4.1a takes the second route, because
`samtools idxstats` — which §4.3(b) depends on — has no `-X` option at all:

```
samtools idxstats -X in/s.cram out/s.cram.crai
  idxstats: invalid option -- 'X'
```

### 3.10 `-P` and `-c` cannot be combined, so the probe writes instead of counting

```
samtools view -P -c ... :  samtools view: The options -P and -c cannot be combined
```

The probe therefore runs the real slice shape and discards the output:

```
samtools view -P -T <candidate> <aln> <region> -o /dev/null   # exit status is the proof
```

Measured against the §3.7 fixture with the `UR:` fallback suppressed: exit **1** for the
chr1-only candidate, exit **0** for the full reference. The probe rejects exactly the
reference the real slice would have failed on.

**A `UR:` target that still exists on disk silently rescues a candidate that would
otherwise fail.** Two of the measurements in this document were wrong until that was
noticed. It is correct behaviour — a resolvable `UR:` *is* a working reference — but it
makes candidate 4 succeed for reasons that have nothing to do with candidates 1–3, and
any test of this path must suppress it.

### 3.11 `idxstats` proves the indexed scan is lossless, from the index alone

`samtools view -f 12 <cram> '*'` returns *unplaced* reads. A flag-12 read (read unmapped
**and** mate unmapped) that is nonetheless *placed* would be missed. That is a
silent-data-loss risk, and a config escape hatch does not discharge it.

`samtools idxstats` reports, per contig, `name / length / mapped / placed-unmapped`, and
reads only the index:

```
samtools idxstats big.cram      wall=0.00s
chr1  20000  601  0
chr2  20000  301  0
*     0      0    4000
```

**If column 4 sums to zero across every contig, no placed unmapped read exists, so `'*'`
provably recovers every flag-12 read in the file.** That is a per-run proof costing one
index read, not an assumption. When it is non-zero the run falls back to the whole-file
stream scan automatically. §4.3(b) is rewritten around this.


### 3.12 The run-local alignment view works, measured end to end

§4.1a is a new mechanism, so it was measured before being specified rather than after.
Setup: an input directory made **read-only** (`chmod -w`) containing a CRAM with no index,
and an output directory holding a symlink to it.

| Operation, through the view | Result |
| --- | --- |
| `samtools index out/view.cram` (bare, no `-o`) | exit 0 — writes `out/view.cram.crai` **beside the symlink**, not beside the target |
| input directory afterwards | `s.cram` only — untouched, and it was read-only throughout |
| `samtools view -c -T ref out/view.cram chr1:4900-5100` | 9 — the index is found through the symlink |
| `samtools view -P -b -T ref out/view.cram <region> -o slice.bam` | exit 0 |
| `samtools idxstats out/view.cram` | full table — the operation that has no `-X` form |
| `samtools view -c -T ref out/view.cram '*'` | 0, matching the `*` row of `idxstats` |
| existing index beside the input, symlinked as `out/view.cram.crai` | resolves; `view` and `idxstats` both work |

Three things this settles:

1. **`samtools index` does not follow the symlink when choosing where to write.** So the
   view needs no `-o` and cannot write into the input tree even by accident — the
   milestone-3 invariant is preserved by the mechanism rather than by remembering a flag.
2. **`-X` is not needed anywhere.** Co-location is restored, so no builder grows an index
   argument and `idxstats` — which has no `-X` — works.
3. **An existing index is symlinked rather than rebuilt**, so the common case costs two
   symlinks and no samtools invocation at all.


### 3.13 The `idxstats` premise, tested against a CRAM that would lose reads

§3.11 assumed `idxstats` column 4 detects placed unmapped reads. Assumption tested on a
purpose-built fixture containing 25 flag-12 **pairs placed on chr1** (RNAME and POS set,
both mates unmapped) alongside 40 genuinely unplaced pairs:

```
samtools view -c -f 12 placed.bam        -> 130     # ground truth
samtools view -c -f 12 placed.bam '*'    ->  80     # what the indexed scan would recover
samtools idxstats placed.bam             -> chr1 20000 600 50   /   * 0 0 80
samtools idxstats placed.cram            -> chr1 20000 600 50   /   * 0 0 80
```

**The indexed scan would have silently dropped 50 of 130 reads on this file, and column 4
reports exactly those 50 — identically for BAM and CRAM.** So `sum(col4) > 0` is a correct
trigger, and `sum(col4) == 0` is a correct proof that `'*'` is complete.

Two further properties, both measured:

* **`idxstats` needs no reference.** With `REF_PATH` pointed at a nonexistent path and no
  `-T`, it returned the correct table for the CRAM and exited 0. The scan decision can
  therefore be taken **before** reference resolution and cannot itself hang.
* **It is cheap, but the spec does not claim it is O(index).** CRAI carries no per-contig
  counts, so samtools derives them from container headers rather than from the index
  alone. On the fixtures here it was not measurably distinguishable from zero; on a WGS
  CRAM it reads container headers, not slices. That is far cheaper than the whole-file
  decode it replaces, which is the only comparison that matters, and it is stated this way
  rather than as a complexity claim that was not measured.

## 4. Design

### 4.1 The seam

One preflight, run immediately after `validate_bam_file` and before any region is
resolved or any stage runs. It resolves and **proves** everything the run will need, and
returns a frozen `AlignmentPlan`. `process_bam_to_fastq` consumes the plan instead of
rediscovering its contents.

`fastq_bam_processing.py` is 649 LOC against a ~650 guideline, and AGENTS.md rule 3 says
to extract the region under change rather than grow the file. So the decisions go into
new **pure** modules and only the subprocess calls stay behind I/O:

| Module | Kind | Owns |
| --- | --- | --- |
| `alignment_contract.py` | pure | what each format requires; the `AlignmentPlan` dataclass; the text of every failure message |
| `reference_resolution.py` | pure | validated ordered candidate mapping and known or unavailable FASTA contig coverage |
| `read_layout.py` | pure | the `paired`/`single`/`mixed`/`empty` verdict from counts, and which FASTQ paths feed downstream |
| `alignment_preflight.py` | I/O | runs samtools: build the alignment view, resolve or build the index, choose the scan from `idxstats`, probe-decode each reference candidate. **It does not compute the read layout** — §4.4 explains why that cannot be known before conversion. |

The three pure modules are unit-testable with no filesystem and are held to the ~100%
branch coverage the existing pure modules (`scoring.py`, `region_utils.py`,
`cohort_rules.py`) already meet.

`fastq_bam_processing.py` measures **649 lines** today, not the 612 AGENTS.md's snapshot
table records. It is therefore already at the ~650 guideline before this milestone adds a
line, so every change to it must be net-negative: the index-resolution block moves into
the preflight and the CRAM branch moves into a helper.

### 4.1a The run-local alignment view

§3.9 measured that an index built into the output directory is **invisible** to
`samtools view`; the naive reading of #225 builds an index and then fails exactly as
before. Two mechanisms can fix that, and only one of them works everywhere:

* `-X <alignment> <index>` on every random-retrieval command. Correct for `view`, but
  `samtools idxstats` — which §4.3(b)'s losslessness proof depends on — **has no `-X`**
  (measured: `idxstats: invalid option -- 'X'`).
* Make the alignment and its index co-located somewhere writable.

So the preflight materialises a **run-local alignment view** inside the output directory:
a symlink to the input alignment, plus its index beside it — the existing index symlinked
when there is one, a freshly built one when there is not. Every subsequent samtools call
uses the view's path.

This costs two symlinks, writes nothing to the input tree (the milestone-3 invariant), and
makes co-location the normal case rather than a special case, so no builder needs an
index argument and `idxstats` works. When the input already has a co-located index the
view is still created, so there is exactly one code path rather than two.

### 4.2 The reference contract — proof, not inference

> **A reference is whatever decodes the target region. Prove it; do not infer it.**

Candidates are tried in this order, and the first one that **decodes a probe slice** wins:

1. `--reference-fasta` (new CLI flag; also a config key for the web worker)
2. config `reference_data.cram_reference_<assembly>`
3. config `reference_data.bwa_reference_<assembly>` — the shipped chr1 FASTA
4. the header's `UR:` field, if it names a readable file

**The probe is the real slice command with its output sent to `/dev/null`** — same `-P`,
same region *or the run's actual BED file*, same reference fragment, same alignment view.
Nothing else is an authorisation: §3.7 measured a `-P`-less probe passing on a reference
the real slice then failed on.

It cannot be a `-c` count: §3.10 measured that **`-P` and `-c` cannot be combined**. The
probe's exit status is the proof, and it discards its output:

```
samtools view -P -T <candidate.fa> <view.cram> <region|-L bed> -o /dev/null
```

(The `-T` is not optional in that line: without it the reference path is read as a second
input alignment.)

When the run has a BED file — `--bed`, `--custom-regions`, or the predefined regions
`pipeline.py` writes — the probe uses that BED, not a region string. A BED naming contigs
the candidate does not cover is exactly the case a region-only probe would miss.

Candidate 3 is retained but demoted in the logs: when the winner does not cover every
contig in the CRAM header, the run logs a warning naming the uncovered contigs, because
its success is a property of this sample's mate placement rather than of the pipeline.

The winner is passed as `-T` to **every** CRAM samtools invocation for the rest of the
run. If no candidate decodes, the run is rejected **at submission** with a message naming
the contig, its `M5`, and every candidate tried with the reason each failed.

Soundness rests on §3.8: htslib verifies each slice's `M5` against the reference before
emitting a record, so "the probe exited 0" means "htslib checksum-verified this reference
against this CRAM", not merely "some bytes came out".

This is deliberately empirical, and that is what makes it correct across every assembly
and naming convention VNtyper supports. A CRAM whose chr1 is called `1` or
`NC_000001.11`, or that was built against a soft-masked or hs38DH chr1 whose `M5`
disagrees with the shipped FASTA, is not reasoned about — it is probed, and it either
decodes or produces a named submission-time error. No ucsc/ensembl/ncbi inference and no
`M5` table are needed anywhere in this path.

**"No reference" and "header `UR:`" are the same probe, and are therefore one candidate,
tried last.** htslib consults the header's `UR:` automatically whenever it needs a
reference and no `-T` was given, so a probe with no `-T` cannot distinguish "this CRAM is
reference-free" from "htslib resolved the reference itself, via `UR:` or via the pinned
local `REF_PATH`". An earlier draft listed them as separate candidates and probed the
no-reference one *first*, which had two defects: `UR:` could win over an explicitly
supplied `--reference-fasta`, and a success was recorded as `reference_path=None` —
"reference-free" — when an external file had in fact done the decoding.

So the order is: **explicit candidates first** (`--reference-fasta`, config
`cram_reference_<assembly>`, config `bwa_reference_<assembly>`), then a final candidate
`no -T`, whose source is recorded as `htslib-resolved (header UR: or REF_PATH)` rather
than as `None`. A `no_ref=1` CRAM reaches that candidate and succeeds there, and because
the winning fragment is empty its commands are byte-identical to today's, so #199's
fixtures are unaffected. Passing `-T` to a reference-free CRAM is harmless, so nothing is
lost by trying the explicit candidates first.

**The `UR:` caveat is load-bearing for anyone testing this.** §3.10: a `UR:` target that
still exists on disk silently rescues a candidate that would otherwise fail, because
htslib falls back to it. That is correct behaviour in production — a resolvable `UR:` *is*
a working reference, and it is candidate 4 — but it means a test asserting that candidate
3 fails must delete or rename the `UR:` target first. Two measurements in §3 were wrong
until this was noticed.

### 4.3 Making the hang unreachable

Two independent changes, either of which alone would leave a path open:

**(a) Pin `REF_PATH`.** For the duration of a run, `REF_PATH` is set to a local-only
value, so no samtools invocation VNtyper makes — including any this milestone did not
think of — can block on a network fetch. Config-driven: `cram.allow_ambient_reference_resolution`,
default `false`. Set it `true` and an operator's deliberately configured refget server is
honoured, with a logged warning that the run can now block on it.

**(b) Bound the unmapped scan, and *prove* the bound is lossless.**
`samtools view -b -f 12 <view.cram> '*'` replaces the whole-file stream — but `'*'`
returns only *unplaced* reads, so a flag-12 read that is nonetheless placed would be
dropped. "It is unplaced in every aligner we tried" is not a guarantee, and a config
escape hatch does not discharge a silent-loss risk.

§3.11 supplies the guarantee. `samtools idxstats` reports placed-unmapped counts per
contig from the index alone, in 0.00 s. So the scan is **selected per run, from evidence**:

```
placed_unmapped = sum(column 4 of samtools idxstats <view.cram>)
if placed_unmapped == 0:  indexed scan   # provably recovers every flag-12 read
else:                     stream scan    # the O(file) form, because it is the correct one
```

Config-driven: `cram.unmapped_scan` is `auto` (default), `indexed` or `stream`.

**The explicit values cannot be used to discard reads.** `indexed` forced on a file whose
`idxstats` says reads would be lost **raises**, naming the count; it does not warn and
continue. An earlier draft allowed it to proceed "for reproducing a historical run", which
would have let a configuration flag defeat the exit bar this spec opens with — a run that
exits 0 having dropped 50 of 130 unmapped reads (§3.13) is precisely the outcome the
milestone exists to make impossible. Reproducing a historical run is what `stream` is for;
it is always safe, only slower.

This is the difference between removing a failure mode and hiding it: the fast path is
taken only where it is proven equivalent, and the proof costs one index read.

### 4.4 No read is silently discarded

`process_bam_to_fastq` produces four FASTQs and its two callers in `pipeline.py` bind two
of them to `_`. `read_layout.py` replaces that with an explicit verdict: every produced
FASTQ is either routed into the run or named in the failure. There is no branch in which
a non-empty FASTQ is dropped without the run saying so.

**The verdict is computed from the produced FASTQs, not from the input slice**, and that
correction matters twice over:

1. *Correctness.* In normal mode the FASTQ source is the merged slice **plus** the
   whole-input unmapped set (`fastq_bam_processing.py:157`, `:202`, `:258`). A perfectly
   paired slice can therefore still yield singletons once unmapped mates are merged in, so
   a layout inferred from the slice would plan for `paired` and then silently strand the
   `other`/`single` files — the very class this section closes.
2. *Constructibility.* An earlier draft had the preflight both return a frozen
   `AlignmentPlan` **and** compute the layout, while the layout was to be read from a slice
   that only exists after `process_bam_to_fastq` runs — a cycle. Splitting it removes the
   cycle: `AlignmentPlan` carries index and reference and is frozen at preflight; the
   layout is a separate verdict computed after conversion and before anything consumes a
   FASTQ.

**`mixed_tolerance` is removed.** A tolerance necessarily coerces some reads into a layout
they do not belong to without saying where they went, which contradicts the invariant this
section exists to state. The rule is exact: a non-empty FASTQ that nothing consumes fails
the run, naming the file and its read count.


### 4.5 Operational contracts the mechanisms need

Three of the new mechanisms read output, create files, or mutate process state. Each needs
a stated contract, because the default behaviour of the code around them is permissive.

**Parsing `idxstats` fails closed.** `run_command` returns a boolean and streams stdout
into a log file (`utils.py:18`), so the scan decision needs its own capture seam rather
than reusing it. The contract: exit status 0; every line exactly four tab-separated
fields; counts parse as non-negative integers; exactly one terminal `*` row. **Anything
else selects the stream scan**, logged with the offending line. An unparsable table must
never be read as "column 4 summed to zero", which is the reading that loses reads.

**The alignment view is created exclusively, and never reused blindly.**
`create_output_directories` reuses an existing output directory silently (`utils.py:128`),
and `--keep-intermediates` makes reruns normal, so a view from a previous run can already
be present. The contract: create the symlink with `os.symlink` and catch `FileExistsError`;
on collision, read the existing link's target with `os.readlink` and accept it **only if
it resolves to the same file as this run's input** (`os.path.samefile`); otherwise replace
it. A stale view pointing at a different alignment would silently genotype the wrong
sample, which is the worst failure this repository can have.

**`REF_PATH` is saved and restored.** Pinning mutates `os.environ`, which is
process-global and outlives a single `run_pipeline` call — and `run_pipeline` is imported
and called in-process by tests and by anything embedding VNtyper as a library. The
contract: capture the previous value (including "unset"), set the pinned one, and restore
it in a `finally`. Without that, one CRAM run silently reconfigures every later run in the
same process.

**Which index the run needs depends on the mode, and only one mode is fussy.** The
non-fast BAM path feeds `extract_unmapped_from_offset`, which parses the BAI container
itself and rejects anything whose first four bytes are not `BAI\x01` — so that mode needs
a BAI specifically, and `resolve_bam_index` stays BAI-only for it. Every other mode
(fast-mode BAM, and CRAM in both modes) hands the index to samtools, which accepts CSI and
CRAI too. The preflight therefore asks for "an index htslib accepts" in those modes and
symlinks the first compatible candidate under the same suffix beside the view: BAM tries
BAI then CSI spellings, while CRAM tries CRAI then CSI spellings. The CRAM CSI support was
measured directly with the pinned samtools 1.20/htslib stack. Only the non-fast BAM path
forces a BAI rebuild when only a CSI is present.

## 5. Per-issue changes

### #213 — `conda run` buffers all output

* **Root cause:** `docker/entrypoint.sh:158,163,168,186` invoke `conda run` without
  `--no-capture-output`, so conda buffers child stdout/stderr until exit. The Docker test
  suite bypasses `entrypoint.sh` via `bash -c` and *does* pass the flag, so no test
  exercises the path users take.
* **Change:** add the flag at all four sites — **and at the fifth, which the issue does
  not mention.** `docker/app/tasks.py:191` builds
  `["conda", "run", "-n", "vntyper", "vntyper", "pipeline", ...]` and runs it at `:312`.
  That is the invocation every *web* job goes through, and it has the same defect. Fixing
  only `entrypoint.sh` makes the Celery worker's own output stream while the pipeline it
  launches stays buffered, so the web path — the one #199 routed CRAM uploads onto — would
  still produce a job that runs for an hour and logs nothing.
* **Failure mode made impossible:** a running pipeline that produces no log output, by
  either entry point.
* **Test:** a source-level assertion that every `conda run` in the repository that launches
  a long-running child carries the flag (so a sixth site added later fails the build), plus
  a Docker test that runs the image **through its entrypoint** and asserts output appears.

### #225 — index resolved too late (BAM) or not at all (CRAM)

* **Root cause:** the region slice is random retrieval and requires an index;
  `resolve_bam_index` runs ~7 lines after it, and CRAM has no equivalent
  (`grep -rn crai vntyper/` returns nothing).
* **Change:** preflight resolves the index before the slice, for both formats, and builds
  a missing one into the **output** directory (`samtools index -o`, available at the
  pinned 1.20) — then exposes it through the run-local alignment view of §4.1a, because
  §3.9 measured that an index merely *placed* in the output directory is invisible to
  `samtools view`. `resolve_bam_index` is left BAI-only and unchanged — the offset
  extractor parses BAI directly and would raise on a CSI. The general
  `resolve_any_index` handles CRAI under both `<file>.crai` and `<stem>.crai`, followed
  by the two equivalent CSI spellings measured against the pinned toolchain.
* **An index that exists is not an index that works.** Resolution is by filename only
  (`Path.exists()`), so a truncated, stale or wrong-sample index resolves happily. The
  preflight does not add a separate validation step for this: the reference probe (§4.2)
  is a real indexed retrieval through the view, so for CRAM a broken index fails there,
  at submission. For BAM the equivalent is a probe slice with no reference fragment. That
  is cheaper and stronger than inspecting the index, because it tests the operation the
  run will actually perform.
* **Fast mode and CSI.** `resolve_bam_index` returning None causes a BAI rebuild even for
  a CSI-indexed BAM in fast mode, which never reaches the offset extractor and would have
  sliced fine through htslib. The preflight therefore asks for a BAI only when the run
  needs the offset extractor (non-fast BAM); otherwise any index htslib accepts is enough,
  which the probe decides.
* **Failure mode made impossible:** `Random alignment retrieval only works for indexed…`
  reaching a user mid-run.
* **Test:** unit tests for both resolvers over all four name spellings; a preflight test
  proving an unindexed CRAM is rejected or indexed before the slice command is built;
  `test_input_tree_is_never_written.py` extended to CRAM.

### #209 — no reference contract

* **Root cause:** `process_bam_to_fastq` sets `cram_ref_option = ""` unconditionally,
  even though both CRAM command builders accept a `-T` fragment.
* **Change:** §4.2 — **and the winning reference must reach every CRAM consumer, not just
  the slice.** `pipeline.py:468` sets `input_bam = Path(cram)` and runs
  `calculate_vntr_coverage` directly against the original CRAM, while
  `build_samtools_depth_command` (`command_builders.py:398`) emits no `-T` at all. So a
  reference-dependent CRAM would pass the probe, slice, convert — and then die at coverage;
  and pinning `REF_PATH` (§4.3a) removes the ambient fallback that hides this today. The
  depth command takes the reference, and coverage runs against the alignment view.
* **The reference is a path, not a preformatted flag.** `cram_ref_option` is interpolated
  **unquoted** by contract (`command_builders.py:187`, `:263`), so building `"-T " + path`
  from a user-supplied `--reference-fasta` would split on spaces and execute on a
  metacharacter under `run_command(shell=True)`. The builders take `reference_path` and
  quote it with `quote_path` themselves.
* **Failure mode made impossible:** an undecodable CRAM being accepted, queued, started
  and failing an hour into a cohort run; and a CRAM that decodes for the slice but not for
  coverage.
* **Test:** the §3.1 ladder as a unit test with samtools mocked; an integration test on a
  reference-dependent CRAM fixture with its reference removed, asserting the run is
  rejected at submission and the message names the reference.

### #178 — hangs indefinitely on `--cram`

* **Root cause:** ambient reference resolution can block without bound (§3.2), and the
  CRAM unmapped scan is O(whole file) (§3.5). #213 removes the evidence of both.
* **Change:** §4.3(a) and §4.3(b).
* **Failure mode made impossible:** a VNtyper run blocking on a network reference fetch.
* **Test:** a unit test asserting `REF_PATH` is pinned to a local-only value before any
  CRAM subprocess is spawned; a command-builder test pinning the `'*'` form; the
  golden-cohort gate proving indexed and stream scans agree.

### #165 — `detect_naming_convention` misresolves on decoy/alt/random contigs

* **Root cause:** `chromosome_utils.detect_naming_convention` divides each convention's
  match count by `len(contig_names)`, so contigs that match *no* convention
  (`chr11_gl000202_random`, `chrUn_gl0002`) inflate the denominator. The issue's example:
  25 UCSC matches out of 93 contigs is 0.27, below the 0.5 threshold, so a header whose
  every primary contig is UCSC returns `unknown`.
* **Change:** score against contigs that match *some* convention — non-matching contigs
  leave the denominator instead of poisoning it — with the primary-contig set and the
  threshold read from config.
* **The tie has to be fixed at the same time.** The existing code tests `>= 0.5` and
  returns the first convention that passes, so a header split exactly 50/50 between two
  conventions returns whichever is checked first (`ucsc`). Narrowing the denominator makes
  that reachable far more often, since the classified-only denominator is small. So a
  convention wins only if it is a **strict majority of classified contigs and no other
  convention ties it**; otherwise `unknown`. The test for this is a two-way 50/50 header —
  a three-way one-third split does not exercise it.
* **Failure mode made impossible:** a header whose primary contigs unanimously agree
  returning `unknown`.
* **Test:** the issue's exact 93-contig shape as a unit test; a genuinely ambiguous
  header still returning `unknown`; `assembly_guard`'s verdict pinned across the change.

### #161 — single-end BAM routes 100% of reads to `output_other.fastq.gz`

* **Root cause:** `samtools fastq -0` collects reads with READ1 and READ2 both unset —
  which is every single-end read — and both `pipeline.py` call sites discard that path
  with `_`. Downstream, `align_and_sort_fastq` and `construct_kestrel_command` require
  two FASTQs.
* **Change:** single-end becomes a supported layout. Preflight determines it from flag
  counts; `read_layout` decides which FASTQs feed downstream; `process_fastq`,
  `align_and_sort_fastq` and `construct_kestrel_command` accept a single FASTQ.
  `mixed` — R1/R2 **and** `other`/`single` all non-empty — is **rejected**, naming every
  file and its record count. It is not coerced and not silently half-processed: VNtyper
  has one k-mer genotyping path and no defined semantics for genotyping two read sets as
  one sample, so inventing one here would be a guess dressed as support. The `paired`
  verdict additionally requires R1 and R2 to hold **equal record counts**; unequal counts
  are `mixed`, because that is a truncated conversion rather than a layout.
  `cli_handlers.py:227` also rejects `--fastq1` without `--fastq2`, so single-end **FASTQ**
  input is refused at argument parsing. That is the same gap seen from the other side, and
  the issue asks about both, so it is fixed here too rather than left as a second report.
* **Failure mode made impossible:** a non-empty FASTQ produced by a stage and consumed by
  nothing.
* **Test:** single-end fixtures derived from the cohort; a unit test that fails if any
  produced FASTQ is neither routed nor named; an integration test genotyping a
  single-end BAM.

### Fixture deriver

`scripts/make_cram_fixtures.py` derives only the fixtures declared in the test-data
config by default, with `--all` for the golden-cohort gate. It gains a
reference-dependent fixture (the `build_reference_dependent_fixture` its own docstring
describes) because #209's path cannot be tested at all without one.

## 5b. Run speed

Robustness must not be bought with wall-clock, and the BAM path is the common one. Three
rules and four measured wins.

### Rules

1. **The preflight costs nothing on the BAM path.** Index resolution is two
   `Path.exists()` calls — no subprocess. The reference probe is **skipped entirely for
   BAM**, which needs no reference. That leaves at most one extra indexed region read.
2. **Read layout is decided from the FASTQs the conversion already wrote**, by reading
   their sizes and, only for a non-empty one, its record count. It is not a second pass
   over the input and not a pass over the slice — §4.4 explains why the slice is the wrong
   source. Missing *inputs* (index, reference) are still decided before any stage runs;
   those are the exit bar. Layout is a property of the output, and deciding it there still
   precedes bwa, fastp, Kestrel and adVNTR by a wide margin.
3. **The CRAM probe decodes one region and discards it.** It cannot be a `-c` count:
   §3.10 measured that `-P` and `-c` cannot be combined, and the probe must carry `-P`
   because the slice does. So it is
   `samtools view -P -T <candidate> <view> <region|-L bed> -o /dev/null`, whose exit
   status is the whole result. Milliseconds, once, CRAM only.

### Measured defects to fix in passing

These are in the files this milestone already opens, and each is a strict improvement.

| # | Defect | Evidence | Fix |
| --- | --- | --- | --- |
| P1 | The region slice runs **single-threaded**. `build_samtools_slice_command` takes no `threads` parameter and emits no `-@`, while merge, fastq, depth, bwa and the CRAM filter all thread. | `command_builders.py:208-212` — no `-@` in the returned string; contrast `:322`, `:361`, `:399`. | Thread the slice. |
| P2 | `samtools index` runs **single-threaded** everywhere. `build_samtools_index_command` has no `threads` parameter. | `command_builders.py:234-236`. | Thread the index; `-@` is supported at the pinned 1.20. |
| P3 | In the default (non-fast) mode the slice's index is **built, invalidated and rebuilt**. The slice command indexes `<name>_sliced.bam`; the merge then `os.replace`s the merged BAM onto that same path, leaving the index stale; the code re-indexes it. | `command_builders.py:211` builds it; `fastq_bam_processing.py:224` overwrites the BAM; `:231` rebuilds the index. | Make indexing an explicit argument of the slice builder and skip it when the merge will overwrite the file. Fast mode keeps it — `pipeline.py:580` consumes `output_sliced.bam` there. |
| P4 | The CRAM unmapped scan decodes and re-encodes the entire file as SAM text. | §3.5. | §4.3(b) — the indexed `'*'` fetch. Largest single win in the milestone, CRAM only. |

P1–P3 are BAM-path wins and apply to CRAM equally. P3 removes one whole index build from
every default run.

### Not changing

`samtools view -P` (`--fetch-pairs`) stays. It costs a second pass and mate seeks, but it
is what makes the slice contain complete pairs, and Kestrel's k-mer genotyping is what
consumes them. Removing it is a genotyping change, not a performance change.

### Gate

`A-PERF-1` below. The golden cohort is the measuring instrument. The first Wave 3 pass
also measured a property the earlier design reviews could not infer from BAM metadata:
32 of the 50 base cases produce a paired FASTQ set **and** at least one non-empty single
FASTQ. The 2.0.9 arm silently ignores those singles; the milestone arm correctly refuses
the mixed layout. Comparing their full wall-clocks would therefore compare completed
genotyping with an intentional early failure. The performance arm is the 18-case
intersection that completes on both revisions; the 32 fail-closed outcomes are recorded
separately and are not reclassified or suppressed to make the timing gate pass.

**One before/after run cannot decide this.** The slice covers 5–10 kb
(`config.json` `bam_region_*`), so `-@` adds thread setup to a job that may be too small
to benefit, and page cache makes a second run of anything faster than the first. The
measurement is therefore: three runs of each arm, alternating, on an otherwise idle host,
reporting median and range; a regression is called only when the slower arm's *best* run
is worse than the faster arm's *worst*. If P1/P2 fall inside the noise band they are kept
anyway — they are strictly less work — but the spec will say so rather than claim a win it
did not measure.

### Wave 3 measurement (2026-08-08)

The host was otherwise idle; the six runs used one harness job, four pipeline threads and
alternated `ddf49a1` (2.0.9) with `678b2d2` (milestone code; the later `4489e06` changes
only the CRAM gate harness). Every counted side verified its package marker, revision and
all 18 expectations.

| Arm | Runs (s) | Median | Range |
| --- | --- | --- | --- |
| 2.0.9 baseline | 88.43, 88.77, 88.68 | 88.68 | 88.43–88.77 |
| milestone | 87.41, 86.08, 85.92 | 86.08 | 85.92–87.41 |

The ranges do not overlap: the milestone arm's worst run is 1.02 s faster than the
baseline arm's best. On the work both revisions can complete, A-PERF-1 is a measured win,
not merely “inside noise.” The separate whole-50 run found 32 explicit mixed-layout
failures on the milestone arm; the baseline completed those cases only by discarding the
stranded FASTQ, so they are exit-bar evidence rather than performance samples.

The forced CRAM measurement also corrected an impossible premise in the original
A-178-2 wording. Both declared golden CRAM sources have non-zero placed-unmapped counts,
so the production guard refuses forced `indexed` before producing an unmapped BAM. Forced
`stream` recorded stable duplicate-case evidence: `7a61` produced 622,690 records with
digest `6677ba29…83cb8b91`; `b178` produced 4,478 with digest
`dad9a81a…4d90a8b`. A raw diagnostic `'*'` fetch on `7a61` produced only 2,690 records and
a different digest — 620,000 reads lost — while `b178` happened to match despite 329
placed-unmapped reads. Therefore equality is required only when preflight authorises both
strategies; otherwise the acceptance result is the indexed rejection plus the stream
count/digest and the measured would-be loss. Requiring the rejected strategy to emit an
equal read set would contradict A-SCAN-1 and reopen silent loss.

## 6. Config surface

Every behaviour above is reachable from configuration. New keys:

```jsonc
"cram": {
  // §4.3(a). false pins REF_PATH to `local_ref_path` for the run, so no samtools call
  // can block on a network refget endpoint. true honours the operator's own REF_PATH
  // and logs that the run may now block on it.
  "allow_ambient_reference_resolution": false,
  "local_ref_path": "%2s/%2s/%s",

  // §4.3(b). "auto" selects indexed when idxstats proves no placed-unmapped read
  // exists, and stream otherwise. The explicit values reproduce a historical run.
  "unmapped_scan": "auto",                      // "auto" | "indexed" | "stream"

  // Order in which reference candidates are probed. Editing this list is the supported
  // way to change resolution policy; the Python fallback is identical for replacement
  // configs that omit this new key.
  "reference_candidate_order": [
    "cli", "config_cram_reference", "config_bwa_reference", "htslib_resolved"
  ]
},
// No `read_layout` block. An earlier draft had `fail_on_unconsumed_reads`, defaulting
// true, "to reproduce pre-milestone-4 behaviour". There is no such key: a switch whose
// false value silently discards reads contradicts the exit bar, and "config-driven" means
// the *policy* is configurable, not that the guarantee is optional.
"reference_data": {
  "cram_reference_hg19": null,
  "cram_reference_hg38": null
},
"assembly_detection": {
  // §5 #165. The threshold was inline; the primary set is what the denominator is
  // computed over, so both belong here.
  "naming_convention_threshold": 0.5,
  "primary_contig_patterns": [
    "^chr[0-9XYM]+$", "^NC_\\d{6}\\.\\d+$", "^([0-9]+|X|Y|MT?)$"
  ]
}
```

`vntyper/config.json` is the shipped default. Because `--config-path` **replaces** the
whole config rather than merging it (trap 2), every read of these keys uses `.get` with
the value above as its default, and a config that omits them behaves exactly like one that
sets them to the shipped values.



## 7. Acceptance criteria

Each is a command whose output decides it. None is satisfied by reading code.

| ID | Criterion | Decided by |
| --- | --- | --- |
| A-213-1 | Every `conda run` that launches a long-running child carries `--no-capture-output`, in `docker/entrypoint.sh` **and** `docker/app/tasks.py`. | source assertion in the unit tier |
| A-213-2 | The image, run through its entrypoint, emits output before the process exits. | docker tier |
| A-225-1 | An unindexed BAM and an unindexed CRAM both complete a region slice, with the built index never written beside the input. | integration, read-only input mount |
| A-225-2 | The slice is performed through a path whose index htslib can actually find — proven by a run whose input directory contains **no** index at any point. | integration |
| A-225-3 | After a run over a read-only input mount, the input tree is byte-identical, for BAM and CRAM. | `test_input_tree_is_never_written.py` |
| A-209-1 | A reference-dependent CRAM with its reference removed **and its `UR:` target renamed** is rejected before any stage runs, with a message naming the contig, its `M5` and every candidate tried. | integration |
| A-209-2 | The same CRAM with `--reference-fasta` pointed at its reference runs to completion. | integration |
| A-209-3 | A `no_ref=1` CRAM runs to completion with no reference supplied, and its commands carry no `-T`. | unit + golden cohort |
| A-209-4 | The reference probe uses the same flags (`-P`) and the same region *or BED* as the slice it authorises. A `-P`-less probe must fail this test. | unit |
| A-209-5 | Coverage (`samtools depth`) on a reference-dependent CRAM carries the same reference the slice used. | unit + integration |
| A-209-6 | A reference path containing a space or a shell metacharacter is quoted, not executed. | unit |
| A-178-1 | With `REF_PATH` pointed at an endpoint that accepts TCP and never responds, a CRAM run **exits within 120 s**. The test asserts on elapsed time; without a deadline the criterion cannot fail. | integration |
| A-178-2 | For every golden-cohort CRAM sample that preflight authorises for both strategies, indexed and stream produce the **same read set** (`samtools view -c` plus sorted read-name digest). When placed-unmapped evidence makes indexed unsafe, forced indexed rejects before work and the gate records the stream read set plus the raw indexed loss; it never bypasses the guard to manufacture equality. | golden cohort |
| A-178-4 | The #178 reporter, on a post-#213 image, supplies `docker logs`; the stage the run reaches is recorded here. This is an evidence request, not a code change, and it does not gate the PR. | issue thread |
| A-SCAN-1 | On a CRAM containing placed flag-12 reads, `auto` selects `stream` and the run recovers all of them; forcing `indexed` **raises** naming the count rather than dropping them. | unit + fixture |
| A-SCAN-2 | A malformed, empty or non-zero-exit `idxstats` selects `stream`, never `indexed`. | unit |
| A-VIEW-1 | A view symlink left by a previous run pointing at a *different* alignment is replaced, not reused. | unit |
| A-VIEW-2 | `REF_PATH` holds its original value (including unset) after `run_pipeline` returns and after it raises. | unit |
| A-165-2 | A header split exactly 50/50 between two conventions returns `unknown`, not the first one checked. | unit |
| A-161-4 | R1 and R2 with unequal record counts are reported as `mixed` and the run fails, rather than being genotyped as a pair. | unit |
| A-178-3 | `auto` selects `stream` on a CRAM constructed to contain a placed flag-12 read, and `indexed` on one without. | unit + purpose-built fixture |
| A-165-1 | The issue's 93-contig header returns `ucsc`; a genuinely ambiguous header still returns `unknown`; a header with no classifiable contig returns `unknown` and does not divide by zero. | unit |
| A-161-1 | A single-end BAM produces a genotype rather than an empty R1/R2 pair. | integration, derived fixture |
| A-161-2 | A run that produces a non-empty FASTQ nothing consumes fails, naming the file and its read count. | unit |
| A-161-3 | `--fastq1` without `--fastq2` is accepted and genotyped rather than rejected at argument parsing. | unit + integration |
| A-PERF-1 | Golden-cohort wall-clock does not regress on BAM cases that complete on both revisions: three alternating runs per arm on an idle host, median and range reported; fail-closed mixed-layout cases are reported separately, never timed as if an early refusal were completed work. A regression is called only when the slower arm's best beats the faster arm's worst. | golden cohort |
| A-PERF-2 | A default (non-fast) BAM run builds the sliced BAM's index **once**, counted from the stage logs. | integration |
| A-WEB-1 | A preflight failure reaches the job status as its own message, not the generic text `main.py:605` substitutes. | web tier |
| A-ALL-1 | `make check-all` passes; `make patch-coverage` ≥ 80%; the coverage floor is not lowered; every function modified by this milestone gains a unit test for the behaviour that changed (AGENTS.md rule 1). | full gate |

## 8. Dependency order

```
wave 0   #213                       independent; lands first so every later failure is legible
wave 1   #225 + the preflight seam  keystone; every later change consumes AlignmentPlan
wave 2   {#209, #178}   #165   #161   fixtures + gate matrix   parallel worktrees
wave 3   golden-cohort gate, full check-all, one PR
```

**#209 and #178 are not independent**, contrary to the milestone brief: both edit the
CRAM branch of `process_bam_to_fastq` and both edit `command_builders.py`. They are one
task in one worktree. #165 (`chromosome_utils.py`) and #161 (`read_layout.py`,
`alignment_processing.py`, `kestrel_genotyping.py`) share no files with them or with each
other.

## 9. Out of scope

* Whole-genome reference download in `install-references`. The contract makes a reference
  supplyable and provable; shipping 3 GB is a separate decision.
* `REF_CACHE` population by `M5` digest. Considered — it would make the contract
  naming-independent without a probe — but the probe already achieves that outcome and
  costs one samtools call instead of a cache builder.
* A web form field for the CRAM reference. The worker reads the configured key; a UI for
  it is a separate change to `docker/app/`.
* #210's second BAM index and #188, both already filed.
* Preflighting inside the HTTP request. §1 explains why: it would make the API run
  samtools synchronously, changing the service's contract. What *is* in scope is surfacing
  the preflight's message in the job status (A-WEB-1).
* A timeout on `run_command`. `allow_ambient_reference_resolution: true` re-enables an
  unbounded wait by definition, and the default `false` removes it; a general subprocess
  deadline is a change to every stage's failure semantics and belongs in its own issue.
  Filed as a follow-up rather than smuggled in here.
* Extracting `construct_kestrel_command` out of the 936-line `kestrel_genotyping.py`.
  AGENTS.md rule 3 asks for it and this milestone touches that function; the extraction is
  included in #161's task, but a wider split of that module is not.

## 10. Adversarial review

Codex `gpt-5.6-sol` at `xhigh` reasoning, run against this spec, against the plan, and
against the final diff. Verdict recorded here. **No HIGH concern may survive into the
PR.**

| Round | Target | HIGH | MED | LOW | Outcome |
| --- | --- | --- | --- | --- | --- |
| 1 | spec | 14 | 12 | 2 | 12 accepted and fixed, 2 rebutted with evidence |
| 2 | spec + plan | 11 | 7 | 4 | 10 accepted and fixed, 1 rebutted by measurement |

### Wave 3 gate evidence before round 3

`make cram-fixtures` derived 50/50 lossless CRAMs with zero skips. The BAM timing result
is the non-overlapping 88.43–88.77 s baseline versus 85.92–87.41 s milestone range in
§5b. The full base matrix also exposed 32 inputs whose previously ignored single FASTQ is
now named and rejected. Forced indexed CRAM runs were rejected by the placed-unmapped
guard; forced stream recorded the 622,690- and 4,478-record read sets, and a raw indexed
diagnostic proved a 620,000-read loss on `7a61`. These are recorded as gate findings, not
hidden by changing scan policy or relaxing the no-discard rule.

### Round 2 — what changed as a result

Round 2's dominant finding was procedural and correct: **the spec had been revised and the
plan had not**, so the plan still contained the invalid `view -c -P` probe, `mixed_tolerance`,
`indexed` as an unconditional default, and layout counted in the preflight. Several
findings recorded as "only appears resolved" were that same staleness seen from different
angles. The plan is regenerated from this spec rather than patched.

Substantive changes: the web error transport now has a mechanism instead of a promise
(`preflight_error.json` -> Redis -> status endpoint, §1), because `run_pipeline` flattens
every exception to `SystemExit(1)` and the status endpoint's genericity is a deliberate
no-leak property rather than an oversight. Both "historical reproduction" escape hatches
are removed: a forced `indexed` scan now raises rather than dropping reads, and
`fail_on_unconsumed_reads` does not exist (§4.3(b), §6). "No reference" and "header `UR:`"
are one candidate tried last rather than two with the no-reference one first, which had let
`UR:` outrank an explicit `--reference-fasta` and had recorded an externally-decoded CRAM
as reference-free (§4.2). New §4.5 states the contracts the new mechanisms need: `idxstats`
parsing fails closed, the view symlink is created exclusively and validated against
`os.path.samefile` on collision, `REF_PATH` is restored in a `finally`, and which index a
run needs is decided per mode so a CSI-only BAM is not rebuilt in fast mode. #165 gains a
tie rule; #161 defines `mixed` as rejected and adds R1/R2 count parity. `AGENTS.md`'s docs
rule is actually corrected this time — round 2 was right that the previous entry claimed an
edit that had not been made.

### Round 2 — rebutted

* **"`auto` relies on a false index-only premise."** Tested rather than argued (§3.13). On
  a CRAM built to contain 25 placed flag-12 pairs, ground truth is 130 flag-12 reads and
  `'*'` recovers 80 — so the indexed scan would indeed have lost 50, exactly as the finding
  feared. `samtools idxstats` reports those 50 in column 4, identically for BAM and CRAM,
  and needs no reference to do it. The premise holds and the trigger fires. The finding's
  *secondary* point is accepted: CRAI stores no counts, so this is not O(index) for CRAM,
  and the spec no longer claims it is — it claims only that reading container headers is
  far cheaper than the whole-file decode it replaces, which is the comparison that matters.

### Round 1 — what changed as a result

Two findings would have shipped a fix that does not fix anything, and both were confirmed
by measurement rather than by argument:

* **An index in the output directory is invisible to `samtools view`** (§3.9). The naive
  #225 fix builds an index and then fails identically. → §4.1a, the run-local alignment
  view.
* **`-P` and `-c` cannot be combined** (§3.10), so the specified probe cannot be run at
  all. → the probe writes to `/dev/null` and its exit status is the proof.

Also accepted and specified: the second `conda run` in `docker/app/tasks.py` (#213's real
web path); the reference never reaching `samtools depth` (#209 × coverage); the web
service returning *Job submitted* before any check runs (§1's honest scope); layout
inferred from the pre-merge slice being both wrong and circular (§4.4); `'*'` losing
placed flag-12 reads (§4.3(b), now proven per run by `idxstats`); `mixed_tolerance`
contradicting the no-discard invariant (removed); an unquoted `-T` fragment being an
injection path (builders take a path); the missing reference-free candidate (candidate 0);
config keys promised but not specified (§6); acceptance criteria with no deadline and no
read-set comparison (§7); the stale 612-line figure (649).

### Round 1 — rebutted

* **"The measured htslib stack is not the repository's pinned stack."** The claim was that
  `samtools=1.20` pins htslib `<1.21`. In this environment `conda list` reports
  `samtools 1.20 h50ea8bc_1` alongside `htslib 1.23 h566b1c6_0`, and `samtools --version`
  reports *Using htslib 1.23* — the measurements were taken on exactly the pair the
  repository resolves. The *underlying* point is nonetheless right and is now stated in
  §3: conda resolves htslib independently of the pin, so the image #178 was reported
  against may differ, and that is an argument for pinning behaviour rather than inheriting
  it.
* **"The plan page violates the repository's docs rule."** AGENTS.md says never add a page
  under `docs/` without registering it in `mkdocs.yml` `nav:`. Measured: an unregistered
  page is an **INFO** in mkdocs 1.6.1, not a warning, so it does not fail `--strict`; what
  actually failed the build was the `macros` plugin reading the prose `{#209, #178}` as an
  unterminated Jinja comment. `exclude_docs: plans/` states the real intent — these are
  contributor working documents, not published pages — and keeps the strict build honest
  about the pages that *are* published. AGENTS.md's wording is updated to match rather
  than the plan being bent to a rule that measurement shows is misstated.
