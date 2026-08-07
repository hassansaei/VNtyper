# New issues to file — defects found on `fix/issue-181-197-followups` and left unfixed

Five bodies, one section each. **Do not file these until the branch merges**, so each can
cite the commit that characterised it. All five are defects the branch found and
deliberately did not fix; the eleven it *did* fix are in the PR body, not here.

Suggested titles are the section headings.

---

# 1. `vntyper cohort` over an empty directory exits 0 and writes no report

**Type:** bug — silent no-op on a successful exit code
**Component:** `vntyper/scripts/cohort_summary.py`, `vntyper/scripts/cli_handlers.py`
**Found by:** the golden-cohort harness's `cohort_empty` case (`cde3aca`)

## Summary

Point `vntyper cohort` at a directory tree containing no `pipeline_summary.json` and it
logs an error, writes nothing but `cohort.log`, and **exits 0**. A caller using
`subprocess.run(..., check=True)` cannot distinguish "the cohort was empty" from "the
cohort report was written".

## Mechanism

`vntyper/scripts/cohort_summary.py:383-385`:

```python
    processed_dirs, temp_dirs = discover_sample_directories(input_paths)

    if not processed_dirs:
        logger.error("No valid input directories or zip files found for cohort aggregation.")
        return
```

A bare `return`, not a raise. `vntyper/scripts/cli_handlers.py:362` then calls
`aggregate_cohort(...)` and **ignores the return value**, so `handle_cohort` falls off its
end and the CLI exits 0.

Note the mismatch inside the module itself: the message is logged at **ERROR**, i.e. the
code already judges this an error condition. Only the exit code disagrees.

## Why it matters

- Any automation wrapping `vntyper cohort` — a pipeline step, a Makefile, a web worker —
  treats exit 0 as "the report is there" and will fail later, further away, on a missing
  file.
- A mistyped `--input-dir`, a zip that failed to extract, or a directory of results from a
  version that did not write `pipeline_summary.json` all land here and all look like
  success.
- The repository's own convention is the opposite: AGENTS.md records `logger.error(msg)`
  followed by `raise ValueError(msg)` with the same message, and only exit codes 0 and 1
  are used. This path logs the error and skips the raise.

## Suggested fix

Raise `ValueError` with the same message, matching the house convention, and let
`cli.py`'s existing handler turn it into exit 1. `handle_cohort` needs no change if the
exception propagates.

Consider separately whether "zero samples" should ever be a *success* — e.g. an explicit
`--allow-empty` for callers that legitimately sweep a possibly-empty tree. If so it should
be opt-in, not the default.

## Suggested acceptance criteria

- `vntyper cohort` over a directory with no `pipeline_summary.json` exits non-zero and
  says why.
- The message is unchanged, so existing log-scraping keeps working.
- A cohort with at least one valid sample is unaffected.
- The golden-cohort harness's `cohort_empty` case is updated to expect the new exit code
  (it currently records the exit-0 behaviour).

---

# 2. Every web job leaves a `.quickcheck.log` in its input directory, so the directory is never reclaimed

**Type:** bug — resource leak on the shared upload volume
**Component:** `docker/app/tasks.py` (worker cleanup) / `vntyper/scripts/utils.py` (`validate_bam_file`)
**Found:** while gathering end-to-end CRAM evidence for #188. Affects **BAM and CRAM alike** and is
entirely independent of #188 — the #188 fixes neither cause nor cure it.

## Summary

`validate_bam_file` writes its `samtools quickcheck` log to a path derived from the **input
alignment**, which for a web job is inside the per-job upload directory on the volume every job
shares. `run_vntyper_job`'s cleanup removes the alignment and its index, then refuses to remove a
non-empty directory — and the directory is never empty, because the quickcheck log is still in it.

Every single web job therefore leaves behind one file and one directory, forever. Nothing else ever
removes them: `delete_old_results` ages out result `.zip` archives in the *output* tree only.

## The two halves

`vntyper/scripts/utils.py:347-349` — the log lands beside the input, not in the output directory:

```python
    command = f"samtools quickcheck -v {quote_path(file_path)}"
    log_file = f"{file_path}.quickcheck.log"
    success = run_command(command, log_file, critical=False, cwd=cwd)
```

`docker/app/tasks.py` cleanup — the removal list names exactly two paths, and the directory is only
removed when nothing is left:

```python
        # Delete the alignment and its index. Both are patient-derived and sit
        # on the volume every job shares, so this is the only thing that ever
        # removes them.
        for path in (bam_path, index_path):
            ...
        # The per-job input directory holds nothing but this submission's own
        # files, so it goes with them. Removing it only when it is empty means
        # anything unexpected left in there is reported rather than deleted.
        job_input_dir = os.path.dirname(bam_path)
        try:
            if os.path.isdir(job_input_dir):
                if os.listdir(job_input_dir):
                    logger.warning(f"Input directory {job_input_dir} still holds files and was left in place")
                else:
                    os.rmdir(job_input_dir)
```

The comment "holds nothing but this submission's own files" is true; the comment that the
`os.listdir` guard exists so "anything unexpected left in there is reported rather than deleted" is
also true and is good design. The bug is that the quickcheck log is **not** unexpected — it is
produced by every run — so the safety guard fires on every job and the intended `os.rmdir` is dead
code in production.

## Observed

Two runs through the real `run_vntyper_job` (only Redis faked), one BAM and one CRAM built from it,
`example_a5c1_hg19_subset`, hg19, `--fast-mode`. Survivors in each per-job input directory
afterwards:

```
{"cram": ["example_a5c1_hg19_subset.cram.quickcheck.log"],
 "bam":  ["example_a5c1_hg19_subset.bam.quickcheck.log"]}
```

with, in both cases:

```
WARNING  app.tasks: Input directory .../input/<job_id> still holds files and was left in place
```

Corroborating evidence that this is long-standing and not an artefact of the harness:
`tests/data/` contains a committed, zero-byte `example_*.quickcheck.log` beside **every one of the
eight** sample BAMs. The same mechanism has been dropping these next to input alignments for a long
time; on a developer's checkout it is litter, on the web service's shared volume it is a leak.

## Blast radius on a long-running shared volume

- **One directory + one inode per job, unbounded.** The count grows monotonically with total jobs
  ever submitted. There is no retention policy on the input tree at all — `delete_old_results`
  only walks `settings.DEFAULT_OUTPUT_DIR` for `.zip` files.
- **Inode exhaustion before disk exhaustion.** The logs are ~0 bytes each, so this will not show up
  as a disk-usage alert; it will show up as `No space left on device` with `df -h` reporting the
  volume nearly empty. That is a confusing outage to diagnose under time pressure.
- **The warning is not actionable and trains people to ignore it.** It fires on 100% of jobs, so
  the log line that was designed to mean "something unexpected happened, look at this" carries no
  information and will be filtered out — which is exactly when it stops being able to report a
  genuinely unexpected leftover, such as a partially written upload.
- **Not a PHI leak, but adjacent to one.** The log's *contents* are `samtools quickcheck -v`
  output, and on a successful check it is empty. Its *name*, however, is the client-supplied
  filename, which for clinical submissions is frequently a patient or specimen identifier. So the
  volume accumulates a permanent, unreclaimed record of every filename ever submitted, in a tree
  whose whole design intent — stated in the cleanup comment above — is that patient-derived input
  does not survive the job.

## Which half should change

**Recommendation: change the log's location, not the cleanup list.** The log is an artefact *of the
run*, so it belongs with the run's other artefacts.

Adding `f"{bam_path}.quickcheck.log"` to the worker's removal tuple would make the symptom go away,
but it is the wrong fix for three reasons:

1. It silently discards the only diagnostic that exists for a failed `samtools quickcheck`. Today,
   a job that dies on a corrupt upload at least leaves the log explaining why; deleting it removes
   the evidence at exactly the moment it is wanted.
2. It couples `docker/app/tasks.py` to an implementation detail of a `vntyper/scripts/utils.py`
   function three call-layers away. The next artefact written beside the input reintroduces the
   same bug, and the same guard swallows it the same way.
3. It leaves the CLI's behaviour untouched, so `tests/data/` keeps accumulating stray logs beside
   the managed sample data.

The right shape is for `validate_bam_file` to take the destination for its log — the run's output
directory, which `pipeline.py` already knows — rather than deriving it from the input path. Note
that `validate_bam_file`'s signature is public-ish and called from more than one place, so the
parameter should default to today's behaviour and be passed explicitly by the pipeline; that keeps
this a contained change rather than a breaking one.

Once the log no longer lands there, `run_vntyper_job`'s `os.rmdir` becomes reachable in production
for the first time, the per-job input directory is actually reclaimed, and the "still holds files"
warning goes back to meaning what it says.

## Suggested acceptance criteria

- A completed web job (BAM **and** CRAM) leaves no file and no directory under the input tree.
- The quickcheck log for a run is present in that run's output directory, and a run that fails
  quickcheck still has its log available afterwards.
- The "still holds files and was left in place" warning does **not** fire on a successful job.
  A test that plants an unexpected file in the input directory should still see it fire — the
  guard must keep working, it just must stop firing on the normal path.
- Running the CLI over `tests/data/*.bam` no longer writes into `tests/data/`.

## Out of scope

The existing committed `tests/data/*.quickcheck.log` files. Removing them is a separate
housekeeping change against managed test data (AGENTS.md), and `tests/test_data_config.json` should
be checked first for whether anything references them.

---

# 3. `Deletion_length` counts deletion *events* and is then used as a length in bp

**Type:** possible correctness defect in adVNTR frameshift filtering — **needs a decision
from @hassansaei**
**Component:** `vntyper/modules/advntr/advntr_genotyping.py`
**Related:** #192 (this was surfaced while implementing the `Insertion_len` summation and
is posted there; filing separately so it is not lost when #192 closes)

## Summary

`Deletion_length` is computed as the number of `D` characters in the `Variant` string:

```python
# advntr_genotyping.py:177 (deletion path) and :216 (insertion path)
df1["Deletion_length"] = df1["Variant"].str.count("D")
```

It is then used **arithmetically as a length in base pairs**, against a value that
genuinely is one:

```python
# :180-184 and :219-223
df1["Insertion_len"]   = df1["Variant"].map(sum_insertion_lengths).astype(int)   # bp, summed LEN tokens
df1["Deletion_length"] = df1["Deletion_length"].astype(int)                      # a COUNT of D events
df1["frame"]           = abs(df1["Insertion_len"] - df1["Deletion_length"]).astype(str)
```

`frame` then selects rows against `ins_frame` (3n+1) and `del_frame` (3n+2), so this
number decides which adVNTR calls are reported at all.

## The assumption this rests on

Counting events equals summing lengths **only if every adVNTR `D` event is exactly one
base.** That is consistent with the State grammar as observed (`D50_2`, `D17_2&D18_2&…`
carry no length suffix, unlike `I…_LEN<n>`), and consistent with @hassansaei's #192
example — but it has not been confirmed against adVNTR's own documentation or source.

If a single `D` token can ever denote more than one deleted base, then:

- `frame` is understated for that row by (deleted bp − 1);
- the row is filtered against the wrong residue class;
- and because 0 is in neither `ins_frame` nor `del_frame`, a row can be **dropped
  silently** rather than mis-scored — the same failure mode #192 fixed on the insertion
  side.

The asymmetry is the tell: insertions carry an explicit `LEN<n>` that is now summed, and
deletions carry no length at all and are counted. Two different quantities feed one
subtraction.

## What is needed

1. Confirmation from @hassansaei (or from adVNTR's `State` grammar) that a `D` event is
   always exactly one base.
2. If yes: rename `Deletion_length` to `deletion_event_count`, or add a comment at both
   sites stating the equivalence and why it holds, so the next reader does not have to
   re-derive it. A test should pin the assumption.
3. If no: `Deletion_length` must parse a length the way `sum_insertion_lengths` now does,
   and the change is genotype-affecting and must go through the golden-cohort gate.

## Notes

- The cohort's only compound adVNTR state, `D17_2&D18_2&D19_2&D20_2&D21_2`, is five
  single-base `D` events, so it is identical under either reading. The golden cohort
  cannot answer this question.
- `Insertion_length` (`str.count("I")`, `:178` / `:217`) has the same shape but is not used
  in the `frame` arithmetic, so it is informational only.

---

# 4. The `POS_fasta` threshold rebase updates a column nothing reads

**Type:** bug — a computed correction is written to the wrong column and discarded
**Component:** `vntyper/scripts/motif_processing.py`
**Found:** while extracting the motif decision layer for #195 (`11e2300`). Pinned by a
characterisation test, not fixed.

## Summary

`motif_correction_and_annotation` intends to derive a repeat-unit-relative coordinate
(`POS_fasta`) from the genome-ish `POS` by subtracting `position_threshold` from rows at or
above it. The comment says so. The code assigns the **un-rebased** `POS` to `POS_fasta` and
then applies the rebase to `POS` — a column no subsequent statement reads.

## Mechanism

`vntyper/scripts/motif_processing.py:490-498`:

```python
        # Adjust POS => create POS_fasta
        combined_df["POS"] = pd.to_numeric(combined_df["POS"], errors="coerce").fillna(-1).astype(int)
        combined_df["POS_fasta"] = combined_df["POS"]          # <- the raw value
        combined_df.update(
            combined_df["POS"].mask(                            # <- the Series is named "POS",
                combined_df["POS"] >= position_threshold,       #    so update() writes back to POS
                lambda x: x - position_threshold,
            )
        )
```

After line 498, `combined_df["POS"]` is never read again. The only columns copied out to
the caller's frame are `Motif_fasta`, `POS_fasta` and `Motif` (`:522-527`). So the rebase
is computed and thrown away, and `POS_fasta` carries the raw `POS`.

## Observed

Pinned today by
`tests/unit/test_motif_decisions.py::test_only_passing_rows_carry_motif_motif_fasta_and_pos_fasta`.
With the shipped `position_threshold = 60`, a row at `POS 67` yields `POS_fasta == 67`. If
the rebase had landed where the comment says, it would be `7`.

## Consequence

`POS_fasta` is consumed by `kestrel_genotyping.generate_bed_file`
(`kestrel_genotyping.py:629-636`), which writes one BED interval per row as
`{Motif_fasta}\t{POS_fasta}\t{POS_fasta + 1}` for IGV. The `Motif_fasta` name is a
repeat-unit identifier, so the coordinate beside it is being interpreted in that
coordinate space — and it is not in it. Any BED consumer that resolves `Motif_fasta` as a
FASTA sequence name gets an offset that is wrong by `position_threshold` for every row at
or above the threshold, and out of range for short repeat units.

`output.bed` is a visualisation aid only. **No genotype, confidence value or report field
depends on `POS_fasta`**, which is why this is filed rather than fixed inside a
genotype-gated branch.

## Which way to fix it

Two readings, and they are not equivalent — that is the decision this issue needs:

- **The comment is right and the code is wrong:** `POS_fasta` should be the rebased
  coordinate. Then `output.bed` changes for every row at or above the threshold, and the
  characterisation test above must be updated to a specification test.
- **The code is right and the comment is wrong:** `POS_fasta` is deliberately the genome
  coordinate, and the rebase is leftover. Then the rebase and its comment should be
  deleted, which is a genuine simplification and removes a mutation-testing blind spot.

Either way, the current state — computing a value and writing it somewhere nothing reads —
should not survive. It reads as working code.

---

# 5. `scripts/` is linted but not type-checked by any gate

**Type:** gap in the quality gates
**Component:** `Makefile`, `scripts/`

## Summary

`RUFF_PATHS` is `vntyper/ docker/app/ tests/ scripts/`, so everything under `scripts/` is
formatted and linted. But `make type-check-all` runs:

```make
type-check-all: type-check          # mypy vntyper/ docker/app/
	mypy vntyper/ tests/
```

`scripts/` appears in **neither invocation**, so nothing under it is type-checked by any
gate — not `make check-all`, not `make ci-local`, not CI's `typecheck` job.

## Why it matters now

`scripts/` is no longer a couple of helper files. It holds:

- `golden_cohort_gate.py` and `scripts/golden_cohort/*` — **the instrument backing every
  genotype claim on this project**
- `coverage_gate.py` — the gate that decides whether CI passes on coverage
- `mutation_test.py` and `mutation_guard.py` — the mutation sweep and the guard that stops
  it baking a mutant into a committed artefact
- `advntr_len_differential.py` — the differential sweep that attests #192
- `download_test_data.py`

The `#181–#197` follow-up branch alone added ~3200 lines there. The code that decides
whether everything else is correct is the code with the weakest gate on it.

## What blocks the one-line fix

Adding `scripts/` to the first mypy invocation:

```make
type-check:
	mypy vntyper/ docker/app/ scripts/
```

currently produces exactly one error, in a file the follow-up branch does not touch:

```
scripts/download_test_data.py:147: error: Need type annotation for "dir_counts"
  (hint: "dir_counts: dict[<type>, <type>] = ...")  [var-annotated]
Found 1 error in 1 file (checked 75 source files)
```

Everything the follow-up branch added under `scripts/` already type-checks clean.

## Suggested fix

1. Annotate `dir_counts` in `scripts/download_test_data.py:147`.
2. Add `scripts/` to `type-check`.
3. Add a note beside `RUFF_PATHS` in the `Makefile` that the ruff and mypy scopes are now
   the same, so the next person to widen one widens the other.

Recorded as AGENTS.md trap 16 in the meantime, so the gap is at least visible.
