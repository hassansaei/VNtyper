# Milestone 4 — CRAM and input robustness

**Status:** spec, awaiting adversarial review
**Branch:** `fix/milestone-4-cram-input-robustness`
**Closes:** #213, #225, #209, #178, #165, #161
**Date:** 2026-08-08

---

## 1. The exit bar

A CRAM run either succeeds, or fails **at submission** with a message naming exactly what
is missing. No input is ever silently discarded. Every behaviour this spec adds is
config-driven; no threshold, contig name, reference path or scan strategy is written
inline in Python.

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

All measurements on samtools 1.20 / htslib 1.23 — the pin in
`conda/environment_vntyper.yml`, and the same pin `v2.0.3` carried, so these results
apply to the environment #178 was reported against.

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
exit. That is #178's reported shape verbatim: *"The process does not exit, no error is
raised, and no output directory/files are generated."*

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

## 4. Design

### 4.1 The seam

One preflight, run immediately after `validate_bam_file` and before any region is
resolved or any stage runs. It resolves and **proves** everything the run will need, and
returns a frozen `AlignmentPlan`. `process_bam_to_fastq` consumes the plan instead of
rediscovering its contents.

`fastq_bam_processing.py` is 612 LOC against a ~650 guideline, and AGENTS.md rule 3 says
to extract the region under change rather than grow the file. So the decisions go into
new **pure** modules and only the subprocess calls stay behind I/O:

| Module | Kind | Owns |
| --- | --- | --- |
| `alignment_contract.py` | pure | what each format requires; the `AlignmentPlan` dataclass; the text of every failure message |
| `reference_resolution.py` | pure | the ordered candidate list and which candidate wins given probe outcomes |
| `read_layout.py` | pure | the `paired`/`single`/`mixed`/`empty` verdict from counts, and which FASTQ paths feed downstream |
| `alignment_preflight.py` | I/O | runs samtools: resolve or build the index, probe-decode each candidate, count the layout |

The three pure modules are unit-testable with no filesystem and are held to the ~100%
branch coverage the existing pure modules (`scoring.py`, `region_utils.py`,
`cohort_rules.py`) already meet.

### 4.2 The reference contract — proof, not inference

> **A reference is whatever decodes the target region. Prove it; do not infer it.**

Candidates are tried in this order, and the first one that **decodes a probe slice** wins:

1. `--reference-fasta` (new CLI flag; also a config key for the web worker)
2. config `reference_data.cram_reference_<assembly>`
3. config `reference_data.bwa_reference_<assembly>` — the shipped chr1 FASTA
4. the header's `UR:` field, if it names a readable file

**The probe is the real slice command with `-c` substituted for `-b -o`** — same `-P`,
same region, same reference fragment. Nothing else is an authorisation: §3.7 measured a
`-P`-less probe passing on a reference the real slice then failed on. The probe region
defaults to the run's own `bam_region` for the same reason.

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

Reference-free (`no_ref=1`) CRAMs, which #199's fixtures use, decode with no `-T` at all;
the probe records that and the run proceeds with no reference, unchanged.

### 4.3 Making the hang unreachable

Two independent changes, either of which alone would leave a path open:

**(a) Pin `REF_PATH`.** For the duration of a run, `REF_PATH` is set to a local-only
value, so no samtools invocation VNtyper makes — including any this milestone did not
think of — can block on a network fetch. Config-driven: `cram.allow_ambient_reference_resolution`,
default `false`. Set it `true` and an operator's deliberately configured refget server is
honoured, with a logged warning that the run can now block on it.

**(b) Bound the unmapped scan.** `samtools view -b -f 12 <cram> '*'` replaces the
whole-file stream. Config-driven: `cram.unmapped_scan`, `indexed` (default) or `stream`.
The `stream` path is kept because `'*'` retrieves *unplaced* reads, and although flag 12
(read unmapped **and** mate unmapped) is unplaced in every aligner output we have
measured, a CRAM that places such reads would lose them silently — so the escape hatch
stays, and the golden-cohort gate proves the default matches the stream on real data
before it ships.

### 4.4 No read is silently discarded

`process_bam_to_fastq` produces four FASTQs and its two callers in `pipeline.py` bind two
of them to `_`. `read_layout.py` replaces that with an explicit verdict: every produced
FASTQ is either routed into the run or named in the failure. There is no branch in which
a non-empty FASTQ is dropped without the run saying so.

## 5. Per-issue changes

### #213 — `conda run` buffers all output

* **Root cause:** `docker/entrypoint.sh:158,163,168,186` invoke `conda run` without
  `--no-capture-output`, so conda buffers child stdout/stderr until exit. The Docker test
  suite bypasses `entrypoint.sh` via `bash -c` and *does* pass the flag, so no test
  exercises the path users take.
* **Change:** add the flag at all four sites.
* **Failure mode made impossible:** a running pipeline that produces no log output.
* **Test:** a Docker test that runs the image **through its entrypoint** and asserts a
  stage log line appears before the process exits.

### #225 — index resolved too late (BAM) or not at all (CRAM)

* **Root cause:** the region slice is random retrieval and requires an index;
  `resolve_bam_index` runs ~7 lines after it, and CRAM has no equivalent
  (`grep -rn crai vntyper/` returns nothing).
* **Change:** preflight resolves the index before the slice, for both formats, and builds
  a missing one into the **output** directory (`samtools index -o`, available at the
  pinned 1.20). `resolve_bam_index` is left BAI-only and unchanged — the offset extractor
  parses BAI directly and would raise on a CSI. A new `resolve_cram_index` handles
  `.crai` under both `<file>.crai` and `<stem>.crai`.
* **Failure mode made impossible:** `Random alignment retrieval only works for indexed…`
  reaching a user mid-run.
* **Test:** unit tests for both resolvers over all four name spellings; a preflight test
  proving an unindexed CRAM is rejected or indexed before the slice command is built;
  `test_input_tree_is_never_written.py` extended to CRAM.

### #209 — no reference contract

* **Root cause:** `process_bam_to_fastq` sets `cram_ref_option = ""` unconditionally,
  even though both CRAM command builders accept a `-T` fragment.
* **Change:** §4.2.
* **Failure mode made impossible:** an undecodable CRAM being accepted, queued, started
  and failing an hour into a cohort run.
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
  `mixed` (both paired and unpaired reads present above a configured tolerance) is
  reported, not silently coerced.
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
2. **Read layout is decided from the slice VNtyper already produces**, not from a second
   pass over the input. The slice is a ~5 kb region that exists on disk by the time the
   layout matters, so counting flags in it is free. Missing *inputs* (index, reference)
   are still decided at submission — those are the exit bar. Layout is a property of the
   data, and deciding it from the slice still precedes bwa, fastp, Kestrel and adVNTR by
   a wide margin.
3. **The CRAM probe is a count, not a copy.** `samtools view -c -T <candidate> <cram>
   <probe_region>` decodes one region and writes nothing. Milliseconds, once, CRAM only.

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

`A-PERF-1` below. The golden cohort is the measuring instrument: it already runs the
whole pipeline over 50 samples, so its wall-clock before and after is the honest number.

## 6. Config surface

Every behaviour above is reachable from configuration. New keys:

```jsonc
"cram": {
  "allow_ambient_reference_resolution": false,  // §4.3(a)
  "unmapped_scan": "indexed",                   // "indexed" | "stream", §4.3(b)
  "reference_probe_region": null                // null = the run's own bam_region
},
"read_layout": {
  "mixed_tolerance": 0.0                        // §5 #161
},
"reference_data": {
  "cram_reference_hg19": null,
  "cram_reference_hg38": null
},
"assembly_detection": {
  "naming_convention_threshold": 0.5            // §5 #165, was inline
}
```

`--config-path` replaces the whole config rather than merging (trap 2), so every read of
these keys uses `.get` with the shipped default and no `KeyError` can abort a run over a
threshold.

## 7. Acceptance criteria

Each is a command whose output decides it. None is satisfied by reading code.

| ID | Criterion |
| --- | --- |
| A-213-1 | The image, run through its entrypoint, emits a stage log line before the process exits. |
| A-225-1 | An unindexed BAM and an unindexed CRAM both reach the slice with an index, or are rejected at submission naming the missing index. |
| A-225-2 | After a run over a read-only input mount, the input tree is byte-identical, for BAM and CRAM. |
| A-209-1 | A reference-dependent CRAM with its reference removed is rejected **before** any stage runs, with a message naming the contig, its `M5` and every candidate tried. |
| A-209-2 | The same CRAM with `--reference-fasta` pointed at its reference runs to completion. |
| A-209-3 | A `no_ref=1` CRAM runs to completion with no reference supplied. |
| A-209-4 | **Settled (§3.7): `-P` does fetch cross-contig mates.** The criterion is now that the reference probe uses the *same* flags and region as the slice it authorises — a `-P`-less probe must fail this test, since it passes on a chr1-only reference the real slice then rejects. |
| A-209-5 | A CRAM whose winning reference does not cover every header contig produces a logged warning naming the uncovered contigs. |
| A-178-1 | With `REF_PATH` pointed at an unresponsive endpoint, a CRAM run completes or fails; it does not block. |
| A-178-2 | Indexed and stream unmapped scans produce identical genotypes across the golden cohort. |
| A-165-1 | The issue's 93-contig header returns `ucsc`; a genuinely ambiguous header still returns `unknown`. |
| A-161-1 | A single-end BAM produces a genotype rather than an empty R1/R2 pair. |
| A-161-2 | No run produces a non-empty FASTQ that nothing consumes and nothing names. |
| A-PERF-1 | Golden-cohort wall-clock does not regress on the BAM path, measured on the same host before and after. P1–P3 are expected to improve it; the criterion is "no regression", and the measured delta is recorded in the plan. |
| A-PERF-2 | A default (non-fast) BAM run builds the sliced BAM's index **once**, proven by counting `samtools index` invocations in the stage logs. |
| A-ALL-1 | `make check-all` passes; `make patch-coverage` is ≥ 80%; the coverage floor is not lowered. |

## 8. Dependency order

```
wave 0   #213                       independent; lands first so every later failure is legible
wave 1   #225 + the preflight seam  keystone; every later change consumes AlignmentPlan
wave 2   {#209, #178}   #165   #161   fixtures        parallel worktrees
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

## 10. Adversarial review

Codex `gpt-5.6-sol` at `xhigh` reasoning, run against this spec, against the plan, and
against the final diff. Verdict recorded here. **No HIGH concern may survive into the
PR.**

| Round | Target | HIGH | MED | LOW | Verdict |
| --- | --- | --- | --- | --- | --- |
| 1 | this spec | _pending_ | | | |
