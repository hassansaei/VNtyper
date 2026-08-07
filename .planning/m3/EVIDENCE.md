# Milestone 3 — exit-criterion evidence

Every claim below is a command and its output, run in the `vntyper` conda environment
(`conda run -n vntyper`). The default shell resolves to miniforge **base**, which has
plotly 6.5.2 against the environment's pinned 6.9.0 — the cohort fingerprint test fails
in base and passes in the environment, so base results are not evidence of anything.

## E1 — a read-only input mount survives a run, and nothing is left behind

Setup: a test BAM plus its index copied into a directory, `chmod a-w` applied, then a
full **non-fast** `vntyper pipeline` (non-fast is the path that indexes, #210).

### On `main` (a detached worktree at b46da80)

```
EXIT=1
PermissionError: [Errno 13] Permission denied:
  '.../e1c/roinput/example_6449_hg19_subset.bam.quickcheck.log'
```

The run dies in `validate_bam_file` before `samtools quickcheck` executes, because
`run_command` opens its log file first. Every BAM and CRAM run, as #162 reports.

### On this branch

```
EXIT=0
Pipeline completed in 0.08 minutes.

=== input dir unchanged? ===
.../roinput/example_6449_hg19_subset.bam: OK
.../roinput/example_6449_hg19_subset.bam.bai: OK
example_6449_hg19_subset.bam
example_6449_hg19_subset.bam.bai

=== output directory ===
advntr  alignment_processing  coverage
example_6449_hg19_subset.bam.quickcheck.log
fastq_bam_processing  igv_report.html  kestrel
pipeline.log  pipeline_summary.json  predefined_regions_hg19.bed
```

`md5sum -c` passes on both input files and the directory listing is unchanged — no new
file, no new directory. The quickcheck log is in the **output** directory, where a run's
artefacts belong, and is still available for diagnosis.

## Adjacent gap found while proving E1 — reported, not fixed

A read-only input mount holding a BAM with **no index at all** still fails, for a
different reason and not because VNtyper writes there:

```
[E::idx_find_and_load] Could not retrieve index file for '.../roinput/example_6449_hg19_subset.bam'
samtools view: Random alignment retrieval only works for indexed SAM.gz, BAM or CRAM files.
```

`fastq_bam_processing.py` slices the region (`samtools view -L`, ~line 171) **before** it
resolves or builds the index (~line 178), so the index this branch now builds into the
output directory arrives too late for the slice that needs it. The input directory stays
pristine throughout — this is samtools declining to do random retrieval without an index,
not a permission failure.

Not fixed here: #162's acceptance criterion is "no VNtyper process writes into the input
directory", which is met, and reordering the slice against the index build is a behaviour
change none of the seven issues asks for. It belongs with the rest of milestone 4
("CRAM and input robustness"), beside #210. Filing it is the honest move; widening the
milestone unilaterally is not.

## E2 — cohort exports are byte-identical across two runs

Two `aggregate_cohort` runs in two processes over the same two web-worker-shaped archives
(`pipeline_summary.json` at the ZIP root), then `diff -r -x '*.html'`:

```
IDENTICAL: every export, the pseudonym table, and every plot

=== sample identity ===
anon_7c32332fff2a   patient_a
anon_b24ff67f5a73   patient_b
```

On `main`, the same two archives, run twice:

```
--- run 1 ---            --- run 2 ---
anon_9b76d  cohort_zip_c_zi4ju1     anon_d5d5a  cohort_zip_0mwx7c3w
anon_77f1b  cohort_zip_er16szzt     anon_c29f2  cohort_zip_el4we5lm
cohort_kestrel.csv: DIFFER
```

The identity is now the input file the run itself recorded, not the extraction directory,
and it is the same on both runs. The HTML is excluded by design: it carries a report
timestamp and Plotly div UUIDs, which is why the oracle has a `_skeleton()` normaliser.

## E3 — no template builds markup from a value it did not author

A malicious `Flag` planted in a stored `pipeline_summary.json` and rendered through
`aggregate_cohort`:

```
PAYLOAD: "></span><img src=x onerror=alert(1)><span title="

raw payload present in rendered HTML? -> False
unescaped '<img src=x' present?       -> False
escaped form '&lt;img src=x' present? -> True
```

and in the shipped templates: zero `.html('<` string-concatenation sinks in either file,
zero `cell.innerHTML`, and `.attr('title', originalText)` in both.

Not proven by execution: this repo has no browser test tier and AGENTS.md forbids one in
the unit tier, so the runtime `.text()` -> `.html()` transition is not exercised. The
residual risk is a future sink the source-text tripwire's pattern does not match.

## Golden-cohort gate — no genotype field moved

Required because cohort output shape moved. Baseline `b46da80`, candidate `2f091fa`,
marker `vntyper.scripts.cohort_pseudonyms` (absent before, present after). 50 CRAM
fixtures derived first, so the matrix is the full **60 cases** and `--allow-matrix-drift`
was never passed.

```
verdict: DELTAS | attestation_grade: True
expectations_unmet: []   blocked_cases: []
launch_verified_both_sides: True   (67 runs per side)
columns_added: none   columns_removed: none   cells_changed: none
```

14 of 25 artefacts show no delta at all. All 11 that do are confined to **one** case,
`cohort_multi_pseudonymized` — and every one of them is identical once the pseudonym
string is stripped:

```
cohort_kestrel_csv    rows 60->60   identical once Sample is stripped: True
cohort_kestrel_tsv    rows 60->60   identical once Sample is stripped: True
cohort_kestrel_json   rows 60->60   identical once Sample is stripped: True
cohort_stats_{csv,tsv,json}         identical once Sample is stripped: True
cohort_advntr_{csv,tsv,json}        identical once Sample is stripped: True
cohort_tables         cells compared, identical once Sample column dropped: True
pseudonymization_table rows 60->60  identical once Pseudonym is stripped: True
```

`cohort_single` and `cohort_multi` — the two non-pseudonymised cohort cases — are byte
`same`. No `Confidence`, `Flag`, `POS`, `REF`/`ALT`, `Motif*`, `Depth_Score` or `VID`
value changed anywhere in the matrix.

### The first run of this gate was invalid, and why

It returned `EXPECTATIONS_UNMET` on six cases, every delta marked `added_after`. Cause:
`git worktree add` does not copy gitignored downloads, so the baseline had no
`reference/vntr_db_advntr/` and its three adVNTR cases died with *"adVNTR genotyping
results not found for cross-match"*, taking the three cohort cases that consume adVNTR
output with them. The candidate side, running in the main tree, had the databases. That
is a broken baseline, not a genotype change. The two downloaded directories were linked
in — the git-tracked reference files still come from the baseline commit, so the
baseline's own configuration stays in force — and the side was re-run.

The harness behaved correctly throughout: it refused to call the run clean, named all six
cases, and `launch_verified_both_sides` confirmed neither side reached the wrong package.
