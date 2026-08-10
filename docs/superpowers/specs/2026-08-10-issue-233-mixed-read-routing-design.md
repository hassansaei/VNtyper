# Issue 233 Mixed-Read Routing Regression Design

**Date:** 2026-08-10<br>
**Issue:** [#233 — Mixed paired/singleton BAM and CRAM inputs fail after conversion](https://github.com/hassansaei/VNtyper/issues/233)<br>
**Severity:** High — a common, valid alignment topology aborts before genotyping<br>
**Status:** Design approved by real-data investigation; implementation pending

## Executive summary

VNtyper 2.0.11 rejects a BAM or CRAM whenever `samtools fastq` emits an equal R1/R2
pair plus one or more singleton or unpaired reads. The failure is not caused by invalid
FASTQ, corrupt pairing, CRAM reference resolution, Docker packaging, or Kestrel. It is
caused by VNtyper's routing policy: every layout containing more than one read category
is labelled `mixed`, and every `mixed` layout is rejected even when the R1/R2 counts are
equal and all produced files can be consumed losslessly.

The regression was introduced by commit
`2b4597d8b57f8986f008baad6c2505359cebea76` while adding single-end routing. Commit
`52c4146596fef2d1e2402991fbab062ba8021889` subsequently changed ten real integration
oracles from successful genotype/artifact assertions to expected exit code 1. The tests
therefore detected the regression but were changed to approve it.

The sustainable fix is to preserve paired-read preprocessing, keep strict R1/R2 parity
validation, and pass every non-empty conversion FASTQ to Kestrel 1.0.1 as positional
inputs belonging to one `-s` sample. Nothing is concatenated, recompressed, or silently
dropped. Unequal or one-sided mate outputs continue to fail closed.

## Scope

### In scope

- BAM and CRAM conversion results containing:
  - a valid equal R1/R2 pair plus singleton reads;
  - a valid equal R1/R2 pair plus `other` reads;
  - a valid equal R1/R2 pair plus both unpaired categories;
  - one or both unpaired categories without paired reads.
- An ordered, lossless read-set interface from conversion routing through Kestrel.
- Exact Kestrel command construction for one through four FASTQ files.
- Local and Docker end-to-end regression coverage using registered real files.
- Restoration of the successful BAM/adVNTR integration assertions that were replaced by
  mixed-layout refusal assertions.
- Correcting documentation that attributes the behavior to the wrong commit or calls it
  intentional after this design supersedes that policy.

### Out of scope

- Changing `samtools fastq` flags or classification semantics.
- Re-pairing orphaned reads, synthesizing mates, or changing read names.
- Concatenating or recompressing FASTQ files.
- Upgrading the pinned Kestrel 1.0.1 JAR or changing calibrated Kestrel parameters.
- Recalibrating confidence, flagging, or clinical-sounding report text.
- Changing adVNTR input behavior; adVNTR continues to consume the alignment.
- Widening CRAM reference or index resolution policy.
- Treating unequal R1/R2 counts as acceptable.

## User-visible failure

For the registered `example_b178_hg19_subset.bam`, VNtyper 2.0.11 emits:

```text
R1: 14690 records
R2: 14690 records
other: 0 records
single: 1 record
FASTQ layout 'mixed' cannot be consumed without dropping reads
```

The run exits 1 before Kestrel, `pipeline_summary.json`, `kestrel_result.tsv`, and the
HTML report are completed. The diagnostic is internally inconsistent: the reads could
be consumed without dropping them because Kestrel accepts multiple sequence files for
one sample.

The singleton is genuine. Its QNAME is
`NB501645:314:HTKKHBGXC:3:21610:22303:10988`; the source alignment contains a primary
READ1 record with flag 99 and only a secondary READ2 record with flag 403. `samtools
fastq` correctly places the unmatched primary in the singleton output.

## Reproduction evidence

All repository-source commands used the `vntyper` mamba environment from the repository
root. Docker commands used read-only input mounts and unique writable output directories.

### Version comparison

| Version | Input | Produced records `(R1,R2,other,single)` | Exit | Result |
| --- | --- | ---: | ---: | --- |
| v2.0.3, source extracted with `git archive` | registered b178 BAM | `(14690,14690,0,1)` | 0 | `High_Precision*`; alt depth 416; active-region depth 7088 |
| v2.0.11 / current `main` | same BAM | `(14690,14690,0,1)` | 1 | routing refusal before Kestrel |

Retained evidence:

- v2.0.3: `/tmp/vntyper-issue233-v203-run.uJ3psq`
- v2.0.11 source: `/tmp/vntyper-issue233-baseline.j8aqGU`

v2.0.3 did not solve the general case: it destructured the four conversion files and
passed only R1/R2 to Kestrel, silently dropping the singleton. That behavior proves the
regression boundary but is not the desired fix.

### Current source input matrix

| Input | Topology | Exit | Interpretation |
| --- | --- | ---: | --- |
| registered b178 BAM | `(14690,14690,0,1)` | 1 | reproduces #233 |
| equivalent original `no_ref=1` b178 CRAM | `(14690,14690,0,1)` | 1 | same routing defect after successful CRAM-container decode; no reference-decode claim |
| remapped b178 CRAM | `(14689,14689,0,0)` | 0 | proves pure-paired CRAM path works |
| b178 paired FASTQ | paired | 0 | paired FASTQ path works |

Retained source-run evidence includes
`/tmp/vntyper-issue233-fastq-local.WiLhmL`,
`/tmp/vntyper-issue233-cram-local.d2ZACH`, and
`/tmp/vntyper-issue233-cram-subset-local.mwDLQC`.

### Published Docker 2.0.11 matrix

The exact tested image was
`ghcr.io/hassansaei/vntyper@sha256:e46bb0e91bb2c526e9ea5fc6870afc5876b20c72f58ed16fe9dc9aeb570899e7`.
It reports VNtyper 2.0.11 and revision
`dc8ec2294d9c424a4df20f5fc8b507f380431c0c`.

| Docker input | Exit | Objective evidence |
| --- | ---: | --- |
| b178 paired FASTQ | 0 | Kestrel TSV and HTML report; `High_Precision*` insertion |
| 6449 single FASTQ | 0 | Kestrel TSV and HTML report; `High_Precision*` insertion |
| b178 BAM | 1 | exact `(14690,14690,0,1)` refusal; no final artifacts |
| equivalent `no_ref=1` b178 CRAM with an explicit image-reference argument | 1 | same routing refusal after CRAM container decoding; the file does not require the supplied reference |
| remapped pure-paired b178 CRAM | 0 | non-empty summary, Kestrel TSV, coverage, and HTML report |

Retained Docker outputs:

- `/tmp/vntyper-issue233-2.0.11-fastq-paired-sr9KQ7`
- `/tmp/vntyper-issue233-2.0.11-fastq-single-Z07bNq`
- `/tmp/vntyper-issue233-bam-bM3Fo5`
- `/tmp/vntyper-issue233-cram-original-s0Aea0`
- `/tmp/vntyper-issue233-cram-remapped-90e65E`

This isolates the defect from image construction and proves that the same mixed-read
routing failure is reachable from a CRAM container. It does **not** prove
reference-compressed CRAM decoding or image reference packaging: the observed b178 CRAM
was generated with `no_ref=1`, so the explicit reference argument was not load-bearing.
The fix therefore adds a separate reference-compressed b178 fixture whose decoded-record
digest must equal its registered source BAM when the exact image-shipped reference is
used, and must fail with a missing or wrong reference.

### Kestrel 1.0.1 capability proof

The vendored JAR's help declares `SEQ_FILE [SEQ_FILE...]`, and `-s` associates every
following sequence file with one sample. Two direct real-data probes used the exact
pinned JAR and calibrated arguments:

| Probe | Inputs under one sample | Exit | Comparison |
| --- | --- | ---: | --- |
| b178 | R1, R2, one-singleton file | 0 | same 6280-record VCF multiset as pair-only; VCF and SAM written |
| 40cf stress case | R1, R2, 93-singleton file | 0 | all 9069 pair-only records preserved; 9 additional records; final result remains `Negative` |

Evidence is retained at `/tmp/vntyper-issue233-kestrel-multifile.HnBW9I` and
`/tmp/vntyper-issue233-worst93.28nbgY`.

## Root cause

The failure chain is deterministic:

1. `process_bam_to_fastq()` emits four named files: R1, R2, unpaired `other`, and
   singleton.
2. `pipeline_read_routing.route_converted_fastqs()` counts every file.
3. `read_layout.classify_layout()` returns `mixed` for every combination not exactly
   pure paired, exactly one unpaired category, or empty. It does not distinguish
   lossless mixed layouts from unequal/truncated mate output.
4. `read_layout.route_fastqs()` routes no file for `mixed` and reports every non-empty
   file as stranded.
5. `route_converted_fastqs()` rejects `mixed` unconditionally and returns at most two
   paths.
6. `pipeline.run_pipeline()` and `kestrel_genotyping.run_kestrel()` expose only
   `fastq_1` and optional `fastq_2`, so the interface cannot represent all four
   conversion outputs.

The abstraction conflates two materially different states:

- **lossless mixed:** equal R1/R2 counts plus valid unpaired reads, or both unpaired
  categories without mates;
- **invalid paired output:** unequal R1/R2 counts or only one mate category.

Rejecting both was simpler than widening the two-file interface, but it converted valid
biological data into an application failure.

## How the regression escaped

It did not escape detection. The process failed at oracle adjudication.

- Commit `2b4597d8...` introduced the unconditional mixed refusal.
- Real integration cases immediately encountered it.
- Commit `52c41465...` changed nine BAM cases and one adVNTR case from exit-0 genotype,
  report, archive, and value assertions to exit 1 plus the exact refusal text.
- The milestone-4 specification explicitly recorded 32 of 50 golden-cohort cases as
  accepted mixed-layout failures and excluded early refusals from performance comparison.
- The Docker quick tier selects `example_b178_hg19_subset_fast`; because that case now
  expects exit 1, the merge-blocking Docker end-to-end check actively approves this
  regression.
- A new derived pure-single-end fixture proved the newly added topology but did not
  protect the existing, more common paired-plus-singleton topology.

The missing invariant was not “does the test run?” It was: **a registered real input
that previously completed may not be reclassified as an expected failure without an
explicit compatibility decision and value-level downstream comparison.**

The milestone-6 harness design also states that `52c41465...` introduced the behavior.
That is historically wrong: `2b4597d8...` introduced it; `52c41465...` normalized its
test expectations.

## Candidate approaches

### A. Restore v2.0.3 drop-with-warning behavior

Route R1/R2 and log that singleton/other reads are ignored.

Advantages: smallest code change; reproduces historical successful results.

Rejected: deliberately loses available evidence; retains the two-file abstraction; a
warning is not a safety property; contradicts the issue's requirement to handle rather
than discard valid reads.

### B. Route every non-empty file to Kestrel under one sample — selected

Validate mate parity, retain fixed R1/R2/other/single ordering, omit empty files, and
pass one through four files after a single Kestrel `-s` sample token.

Advantages: lossless; directly supported by Kestrel 1.0.1; no temporary merged file;
preserves paired FASTQ preprocessing; idempotent; smallest coherent interface change.

Risk: retained singleton evidence can legitimately change a call, so restored cohort
oracles must be derived from real runs rather than copied from v2.0.3.

### C. Concatenate or recompress all reads into one FASTQ

Rejected by experiment. A real b178 whole-pipeline probe retained 16,372 reads after
single-end fastp instead of 21,138 through paired preprocessing and changed selected
Kestrel motif/depth evidence. It also adds disk use, partial-write recovery, and quoting
surfaces without improving Kestrel compatibility.

## Selected architecture

```text
samtools fastq
  -> four files + exact record counts
  -> classify: paired | single | mixed | invalid | empty
  -> validate parity and choose every non-empty file in fixed order
  -> immutable tuple[str, ...]
  -> Kestrel command: one -s sample + N positional FASTQs
  -> existing VCF/SAM processing and report pipeline
```

### Read-layout contract

`classify_layout(r1, r2, other, single)` gains an `invalid` verdict.

| Counts | Verdict | Downstream files |
| --- | --- | --- |
| all zero | `empty` | none; fail |
| R1 = R2 > 0; unpaired zero | `paired` | R1, R2 |
| R1 = R2 > 0; either unpaired category nonzero | `mixed` | every non-empty file |
| R1 = R2 = 0; exactly one unpaired category nonzero | `single` | that file |
| R1 = R2 = 0; both unpaired categories nonzero | `mixed` | other, singleton |
| R1 != R2 or exactly one mate category populated | `invalid` | none; fail |
| any negative count | exception | none; fail |

`route_fastqs()` returns every selected path exactly once in the stable order R1, R2,
other, singleton. The routing layer must prove that no non-empty path is stranded.

`route_converted_fastqs()` returns `tuple[str, ...]` with length 1–4. It logs the layout,
each file, each record count, and the exact selected file list. `empty`, `invalid`, a
stranded non-empty file, or an inconsistent router result raises `ValueError` before
Kestrel.

### Kestrel interface

The Kestrel command builder and runner accept `Sequence[str | Path]`, normalize it once
to a non-empty immutable tuple, reject a scalar `str`/`Path` container and
empty-string/`None`/unsupported elements, quote each path independently, and preserve
order. A relative operand whose string starts with `-` is rejected because shell quoting
does not stop Kestrel from interpreting it as an option; absolute paths and relative
paths with a directory prefix are unaffected. The command has exactly one sample option
followed by all files:

```text
-sSAMPLE R1.fastq.gz R2.fastq.gz other.fastq.gz single.fastq.gz
```

There is no duplicated `-s`, no literal `None`, no shell glob, and no concatenated path
token. `additional_settings` may not override sample/input grouping with `-s`,
`--sample`, or `--filespersample`. Existing one- and two-file command bytes remain
unchanged apart from the internal Python calling convention.

### Component boundaries

- `read_layout.py`: pure verdict and selection logic.
- `pipeline_read_routing.py`: FASTQ I/O, record counting, diagnostics, and immutable
  selected tuple.
- `kestrel_command.py`: pure quoting and command construction for N inputs.
- a focused Kestrel execution-planning module: extracts the touched configuration and
  invocation planning from the over-limit `kestrel_genotyping.py`.
- a focused pipeline Kestrel-stage module: removes the touched orchestration region from
  the over-limit `pipeline.py` and forwards the tuple without changing it.
- shared integration orchestration: one typed semantic request and captured result for
  local and Docker runners; one argv builder owns input kind, threads, log level, CLI
  options, and normalized resource identity, while runners only map paths and execute.
  Shared validation owns routing evidence, archives, artifacts, and values.
- `scripts/integration_compatibility.py`: a tested executable core that validates live
  success contracts and append-only history together; its thin CLI receives the explicit
  event base from GitHub Actions.
- golden-cohort case expectations: previously measured mixed IDs become ordinary success
  cases; only genuinely invalid parity retains a causal failure oracle.

The implementation must not grow either over-limit module in place.

## Failure, recovery, retry, and idempotency

- Missing, unreadable, undecompressible, or incomplete four-line FASTQ output fails
  during counting with the path and causal error. Full FASTQ header/quality syntax
  validation is not introduced by this change.
- Unequal or one-sided mates fail before Kestrel with all four counts.
- Empty output fails before Kestrel.
- Any non-empty produced file not present in the selected tuple fails as an internal
  routing invariant violation.
- Kestrel failure behavior remains unchanged and names the k-mer/log file.
- A retry re-counts the actual conversion outputs and reconstructs the same ordered
  tuple. No merged input artifact exists to become stale or partial.
- Existing stale-VCF cleanup remains before invocation.
- Archive protection remains bound to original CLI inputs and explicit reference/BED/BWA
  paths. The routed tuple is never interpreted as one filesystem path.

## Compatibility and security

- Python 3.10-compatible type syntax only.
- Kestrel remains at 1.0.1 with unchanged k-mer/state/memory settings.
- All paths remain individually shell-quoted because `run_command` intentionally uses
  `shell=True`.
- Input mounts and source files remain read-only; routing performs no input writes.
- FASTQ and pure-paired BAM/CRAM behavior is preserved.
- SHARK, fastp, and BWA keep their existing one/two-file interfaces. The tuple is
  introduced only at the post-alignment Kestrel routing boundary.
- The internal two-argument Kestrel Python interface changes to a sequence; all in-repo
  callers and tests change atomically. It is not a documented public API.
- adVNTR still consumes the alignment. The change only prevents an earlier Kestrel-stage
  routing abort, after which existing adVNTR and cross-match behavior executes.
- Local and Docker real-data cases use the same frozen thread count, log level, normalized
  CLI options, and resource identities. Transport path prefixes may differ; no runner
  supplies a hidden thread or logging default.

## Observability

One INFO diagnostic must contain:

- the verdict;
- all four path/count pairs;
- all selected paths in execution order.

It is emitted as one canonical `READ_SET_ROUTING` JSON record with stable basenames. The
shared local/Docker validator requires exactly one record and compares counts and
selection to the case declaration. This proves no-drop behavior even when a genotype
value happens to remain within tolerance.

Invalid parity must use an ERROR diagnostic that says the mate outputs are inconsistent,
not that valid singleton reads “cannot be consumed.” Tests assert stable causal phrases
and values, not an entire path-dependent line where that would make harmless formatting
changes expensive.

## Test strategy

### Layer 1: pure and mutation-sensitive unit tests

- Complete count truth table, including equal pair + 1 singleton and +93 singletons.
- `invalid` for unequal/one-sided mates.
- Every non-empty file selected exactly once and in fixed order.
- Empty and negative counts fail.
- N-file Kestrel command exact string, hostile paths, empty sequence, empty elements, no
  duplicated sample option, and no `None` token.
- Scalar sequence-container rejection and leading-dash option confusion.
- Mutations that omit singleton, reorder files, accept unequal pairs, or emit a second
  sample option must be killed by named assertions.
- Mutations that enable `--filespersample`, pass a tuple into archive ownership, or
  restore the golden-cohort mixed-failure declaration must also be killed.

### Layer 2: pipeline wiring tests

- BAM and CRAM conversion tuples reach Kestrel unchanged.
- FASTQ preprocessing still occurs before post-alignment conversion.
- Single FASTQ and paired FASTQ remain one/two files after conversion.
- Kestrel-stage failure prevents false summary/report success.

### Layer 3: real local integration

At minimum:

- b178 BAM `(14690,14690,0,1)` completes and asserts Kestrel values, report, summary,
  coverage, archive state, and retained intermediates;
- 40cf BAM `(3474,3474,0,93)` completes and remains `Negative`;
- b178 and 40cf run in default mode as well as fast mode; the reporter's default-mode
  path has its own counts, values, and artifacts rather than inferred coverage;
- the existing `no_ref=1` original-derived b178 CRAM completes as a CRAM-container and
  mixed-routing control;
- a separately generated reference-compressed b178 CRAM completes only with the exact
  registered reference, has a decoded-record digest equal to the source BAM, and
  exercises the same `(14690,14690,0,1)` mixed topology; real missing-reference and
  wrong-reference invocations fail during decode before any success artifact;
- pure-paired remapped b178 CRAM still completes;
- b178 paired FASTQ and 6449 single FASTQ still complete;
- the mixed a5c1 adVNTR case completes and asserts both Kestrel and adVNTR values.

FASTQ controls assert Kestrel values, coverage, summary, and report content. The mixed
adVNTR case additionally asserts the cross-match step/result. Restored values must come
from fresh final-tree runs. A test may not assert only exit 0 or artifact existence.

### Layer 4: Docker end to end

- The quick Docker case remains b178, but it must now require success and non-empty
  Kestrel, coverage, summary, report, and archive artifacts.
- Docker must run the same real b178 BAM success oracle through shared orchestration.
- Add Docker coverage for b178 paired FASTQ, 6449 single FASTQ, the `no_ref=1` original
  mixed b178 CRAM, the reference-compressed mixed b178 CRAM, and the pure-paired
  remapped CRAM. The runner is the only Docker-specific layer.
- Docker fast/full preparation deterministically runs `cram-fixtures` before collecting
  CRAM parameters; a locally pre-existing ignored fixture is not accepted as CI proof.
- The reference-compressed CRAM uses the image-shipped reference explicitly and checks
  the decoded-record digest, proving packaging and decode; `no_ref=1` cases are not
  accepted as evidence for that claim.
- Main/full tier retains all registered BAM cases; scheduled full tier retains adVNTR.
- A YAML-semantic contract pins PR to quick, main push to fast, and schedule/manual to
  full. Testing a Make target that the workflow no longer calls is not sufficient.

### Layer 5: compatibility guard

A versioned success-baseline manifest, independent of the mutable test-case declarations,
records every established real-data success. Its identity is `(suite, test_name)`, so
the same input may legitimately appear in BAM and adVNTR suites. Each row fingerprints
the input/resource digest, reference, normalized CLI options/modules, routing evidence,
required artifacts/archive state, value fields, expected values, and bounded tolerances.
Two checks enforce it:

1. the executable validator requires each baseline identity to resolve to one matching
   live exit-0 case and every qualifying live real success to have a baseline row; and
2. the same executable rejects deletion or mutation of any baseline identity relative
   to the explicit event base.

The initial seed is independently reconstructed from
`52c4146596fef2d1e2402991fbab062ba8021889^` and pins the exact nine BAM plus a5c1
adVNTR successes; absence of a base manifest is never an unchecked bootstrap. The check
is part of `make check-all` and an explicit full-history GitHub Actions step. Pull
requests compare against `origin/${{ github.base_ref }}`; main pushes compare against
`${{ github.event.before }}`. Missing, all-zero, shallow, or unreachable bases fail
closed. This avoids the ineffective `origin/main == HEAD` comparison on direct pushes.

Changing `expected_exit_code`, invocation parameters, values/tolerances, renaming the
case, or deleting its baseline row cannot turn a regression green. A future deliberate
breaking change requires an explicit repository-policy change reviewed as such; it cannot
be encoded as free text beside the failing case. Issue #233 adds no exception mechanism.

## Acceptance criteria

| ID | Objective acceptance check |
| --- | --- |
| AC233-1 | Equal R1/R2 plus any non-empty unpaired category is routed losslessly, never rejected merely as `mixed`. |
| AC233-2 | Every non-empty conversion file reaches one Kestrel sample exactly once in R1/R2/other/single order. |
| AC233-3 | Unequal or one-sided mates, empty output, unreadable/undecompressible/incomplete FASTQ, and stranded paths still fail before Kestrel with causal diagnostics. |
| AC233-4 | Existing one- and two-file Kestrel commands retain their exact behavior; paths remain independently shell-quoted. |
| AC233-5 | Real b178 BAM in default/fast mode, the equivalent `no_ref=1` CRAM, and a separately generated reference-compressed b178 CRAM complete locally and in Docker with non-empty genotype, coverage, summary, and report artifacts; the reference-compressed case binds the reference SHA-256, proves exact image-reference decode by source-record digest, and fails real wrong/missing-reference controls. |
| AC233-6 | Real paired FASTQ, single FASTQ, pure-paired CRAM, 93-singleton BAM, and mixed-input adVNTR/cross-match cases complete with value-level assertions. |
| AC233-7 | Tests kill omission, reordering, parity-weakening, extra-sample, and exit-only false-positive mutations. |
| AC233-8 | The Docker quick merge gate and workflow mapping require a successful mixed b178 pipeline and cannot pass on the current exit-1 behavior. |
| AC233-9 | One executable live/bidirectional validator plus event-base append-only guard prevents changing invocation/outcomes or renaming/deleting an established success to normalize a regression. |
| AC233-10 | Current unit, coverage, patch, integration, Docker, formatting, lint, typing, and Python 3.10 gates pass without weakening thresholds. |
| AC233-11 | Documentation attributes introduction to `2b4597d8...`, oracle inversion to `52c41465...`, and marks the old unconditional refusal policy superseded. |

## Rollout and release behavior

The fix is backward-compatible for valid inputs but can change genotype evidence because
previously dropped singleton reads are now retained. The PR must report fresh golden
cohort diffs rather than promise byte-identical results. Any changed result requires a
read-level explanation and review; it must not be hidden by broadening tolerances.

After merge, publish through the repository's normal release automation. Do not create,
move, or push a tag from the implementation task. The release notes should call out that
valid paired alignments containing singleton/unpaired reads no longer abort and that
those reads are retained rather than discarded.

## Unresolved questions

None. The shipped Kestrel capability, desired no-drop policy, invalid-parity behavior,
and real-data acceptance matrix are all established by current evidence.
