# Mutation Testing

!!! warning "Advisory only - nothing gates on this number"

    This score is **not** a pass/fail threshold and CI does not bind against
    it. Read it as a map of which decisions are untested, not as a grade.

Line coverage answers *did a test execute this line?*. It does not answer
*would a test have noticed if this line were wrong?* - and for VNtyper the
second question is the one that matters, because the characteristic failure
is a silently wrong genotype call rather than a crash.

Mutation testing answers the second question directly: it introduces a
deliberate defect, runs the tests, and records whether anything failed. A
**surviving** mutant is a defect the suite cannot see.

## Result

**177 of 220 mutants killed - a raw mutation score of 80.5%.**

Of the 43 survivors, 0 are hand-classified as
*equivalent* (the mutation cannot change observable behaviour, so no test could
ever kill it) and 43 are genuine gaps. Excluding the equivalent
mutants the score is **80.5%** (177/220).

Both numbers are given because neither alone is honest: the raw score
understates the suite by counting unkillable mutants against it, and the
adjusted score depends on classifications that are a human judgement call.
Every classification is listed below with its reason so it can be checked.

| Module | Killed | Total | Raw score |
| --- | ---: | ---: | ---: |
| `vntyper/scripts/kestrel_genotyping.py` | 60 | 85 | 70.6% |
| `vntyper/scripts/motif_processing.py` | 43 | 59 | 72.9% |
| `vntyper/scripts/motif_decisions.py` | 7 | 8 | 87.5% |
| `vntyper/scripts/confidence_assignment.py` | 10 | 11 | 90.9% |
| `vntyper/scripts/scoring.py` | 20 | 20 | 100.0% |
| `vntyper/scripts/variant_parsing.py` | 7 | 7 | 100.0% |
| `vntyper/scripts/flagging.py` | 30 | 30 | 100.0% |

## Surviving mutants

### Genuine gaps

Each of these is a change to the source that the **entire** unit tier does
not notice. They are the actionable output of this exercise: a test that
kills one is a test that would have caught a real defect of that shape.

| Module | Line | Mutation |
| --- | ---: | --- |
| `vntyper/scripts/confidence_assignment.py` | 138 | `-` &rarr; `+` |
| `vntyper/scripts/kestrel_genotyping.py` | 182 | `and` &rarr; `or` |
| `vntyper/scripts/kestrel_genotyping.py` | 220 | `continue` &rarr; `break` |
| `vntyper/scripts/kestrel_genotyping.py` | 284 | `20` &rarr; `21` |
| `vntyper/scripts/kestrel_genotyping.py` | 285 | `30` &rarr; `31` |
| `vntyper/scripts/kestrel_genotyping.py` | 286 | `30` &rarr; `31` |
| `vntyper/scripts/kestrel_genotyping.py` | 368 | `break` &rarr; `continue` |
| `vntyper/scripts/kestrel_genotyping.py` | 492 | `False` &rarr; `True` |
| `vntyper/scripts/kestrel_genotyping.py` | 608 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 611 | `True` &rarr; `False` |
| `vntyper/scripts/kestrel_genotyping.py` | 635 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/kestrel_genotyping.py` | 639 | `or` &rarr; `and` |
| `vntyper/scripts/kestrel_genotyping.py` | 777 | `False` &rarr; `True` |
| `vntyper/scripts/kestrel_genotyping.py` | 876 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 933 | `3` &rarr; `4` |
| `vntyper/scripts/kestrel_genotyping.py` | 940 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 944 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 946 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 949 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 951 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 954 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 956 | `0` &rarr; `1` |
| `vntyper/scripts/kestrel_genotyping.py` | 963 | `1` &rarr; `2` |
| `vntyper/scripts/kestrel_genotyping.py` | 1031 | `False` &rarr; `True` |
| `vntyper/scripts/kestrel_genotyping.py` | 1083 | `==` &rarr; `!=` |
| `vntyper/scripts/kestrel_genotyping.py` | 1083 | `1` &rarr; `2` |
| `vntyper/scripts/motif_decisions.py` | 86 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 216 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 231 | `==` &rarr; `!=` |
| `vntyper/scripts/motif_processing.py` | 239 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 290 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 292 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 295 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 330 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 337 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 341 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 367 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 491 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 510 | `-` &rarr; `+` |
| `vntyper/scripts/motif_processing.py` | 510 | `1` &rarr; `2` |
| `vntyper/scripts/motif_processing.py` | 523 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/motif_processing.py` | 525 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/motif_processing.py` | 534 | `False` &rarr; `True` |

### Classified equivalent

None classified yet.

## How this compares to the 43.5% baseline

The experiment that motivated this work scored **43.5%** (27 of 62 mutants
killed) across the eight highest-coverage modules. That harness was never
committed, which is why this one exists.

!!! note "The two totals are not directly comparable"

    Different mutant population, different modules: 62 mutants over eight
    modules then, 220 over 7 modules now, generated by a different operator set.
    A higher or lower headline number would not by itself mean the suite has
    improved or regressed. Only per-module figures on the same module carry
    across, and even those only loosely.

The one comparison that is meaningful is `confidence_assignment.py`, the
module that motivated the whole effort: it had **100% line coverage and a 21%
mutation score**, i.e. four of five deliberate defects in it went undetected by
a fully green build.

| `confidence_assignment.py` | Then | Now |
| --- | ---: | ---: |
| Line coverage | 100% | 100% |
| Mutation score | 21% | 90.9% raw, 90.9% adjusted |

## Reproducing this

```bash
make mutation
```

The harness is `scripts/mutation_test.py`. It mutates one token at a time,
runs the module's own tests first and escalates anything that survives to
the full unit tier before recording it as a survivor, so the score is not
biased by the scoping.

Before it mutates anything it runs those same pytest invocations against the
**unmutated** tree and aborts unless they pass, printing the failure. A mutant
is counted as killed whenever pytest exits non-zero, and pytest exits non-zero
for an unrelated failure, a collection error or a missing dependency just as
readily - so without that preflight a broken checkout scores 100% and
overwrites this page with the result.

!!! danger "No child may load bytecode built from a different revision"

    That is the invariant. CPython validates a cached `.pyc` against the
    source's `(mtime, size)` pair with one-second mtime granularity, so a
    mutant written in the same second as the file it replaces **and of the
    same byte length** (`==` to `!=`, `1` to `2`) is indistinguishable from
    it to the cache validator: the interpreter loads the stale `.pyc`, runs
    the **unmutated** code, every such mutant *survives* and the score is
    fiction. Two sweeps produced exactly that before it was found.

    Two defences hold it, and both are required. `run_pytest()` passes
    `python -B` and sets `PYTHONDONTWRITEBYTECODE=1` in the child, so no
    `.pyc` is written during the run; and every `__pycache__` under
    `vntyper/` is deleted before the sweep starts and again after each
    mutant is written, so none left by an earlier run or an earlier mutant
    can be loaded. The flags stop new caches, the deletion stops old ones.

    The `PYTHONDONTWRITEBYTECODE=1` on the `make mutation` recipe is defence
    in depth for the parent process, which never imports a target module -
    it is not what holds the invariant, and the harness is safe without it.

!!! note "Each measurement runs in an isolated workspace"

    The harness captures HEAD in a disposable detached worktree and overlays
    the current non-ignored working state, except selected mutation targets
    and requested output paths. Selected targets therefore come from the
    captured commit, while ordinary edits and new tests participate in the
    measurement without being written back.

    Import provenance is proved against the pinned worktree before testing.
    A green baseline and a known-killed canary must then pass before ordinary
    mutants are measured, and the post-overlay baseline is verified after
    the canary and after every target.

    Every mutant and bytecode-cache write is confined to that workspace;
    real production source is never mutated. Requested report artifacts are
    built completely and installed atomically in the real checkout.
    The cleanup is best effort: SIGINT, SIGTERM, SIGHUP and SIGQUIT attempt
    the common unwind path, while SIGKILL or a host crash can leave only an
    orphan disposable worktree for later inspection and removal.

## Related: branch coverage, now enabled

Mutation testing and branch coverage were investigated together, because both
ask a sharper question than line coverage and they agreed on which modules are
weakest. The branch-coverage half of that work is recorded here so it is not
re-derived from scratch.

`[tool.coverage.run]` now sets `branch = true`, so an `if` that is entered but
never taken no longer counts as fully covered. It was enabled in **#196**,
measured on `fix/issue-181-197-followups` at `5bb2463`:

| Measure | Value |
| --- | ---: |
| Line (statement) coverage | 76.60% |
| **Branch-inclusive total** | **74.22%** |
| Branch-only coverage | 66.00% |
| Branch exits never taken | 512 of 1506 |

`fail_under` was raised **70 &rarr; 74** in the same commit, to the figure
`scripts/coverage_gate.py` printed for that run. **The floor was raised to meet
the measurement; the measurement was not weakened to fit the floor.** That
distinction is the whole point of the ratchet.

!!! warning "74 is a branch-inclusive floor, and nothing else notices if that changes"

    Deleting `branch = true` does not fail any gate on its own - it *raises*
    the reported total, because statement-only coverage of the same suite is
    76.60% against the branch-inclusive 74.22%. The build would go green while
    measuring strictly less. `tests/unit/test_coverage_gate.py::test_branch_coverage_is_enabled`
    exists solely to make that edit fail, and it is the only thing that does.

### Correction: the previously recorded prerequisite was wrong

This section formerly recorded the opposite conclusion, and it is kept here
rather than deleted, because a document that quietly rewrites its own history
stops being worth trusting.

The earlier measurement was **63.80%** branch-inclusive (76.60% is the current
line figure; it was 66.82% then) against a `fail_under` of **66**. Enabling
branch coverage at that point really would have failed CI on the enabling
commit, and the decision not to enable it - and specifically not to lower the
floor to admit it - was correct.

What was wrong was the stated route out. The old text identified
`cohort_summary.py` and `install_references.py` as holding 275 of the 685
untaken exits, both on the oversized-file list in `AGENTS.md`, and concluded:
*"Splitting them is the prerequisite, not writing more tests against them."*

**It was not a prerequisite.** Branch coverage cleared the floor with both
files still unsplit and still untested. The gap was closed instead by testing
five small, already-testable modules to 100% and a sixth to 98%:
`cross_match.py`, `utils.py`, `file_processing.py`,
`extract_unmapped_from_offset.py`, `variant_parsing.py`, and `docker/app/tasks.py`.

The generalisable mistake was reading "these two files hold the most untaken
exits" as "these two files are the ones that must be fixed". Concentration of
missing coverage is not the same as cheapness of covering it: the two oversized
modules are expensive precisely because they fuse I/O with logic, while the
same number of units was available across several modules that could be called
directly. Splitting `cohort_summary.py` and `install_references.py` remains
worth doing for the reasons `AGENTS.md` gives - it was simply never a blocker
for this.

## Raw output

```text
VNtyper mutation testing - advisory score
============================================================

Command:  make mutation
Total:    220 mutants, 177 killed, 43 survived
Score:    80.5%
Duration: 57.4 min

Per module
------------------------------------------------------------
  70.6%   60/ 85  vntyper/scripts/kestrel_genotyping.py
  72.9%   43/ 59  vntyper/scripts/motif_processing.py
  87.5%    7/  8  vntyper/scripts/motif_decisions.py
  90.9%   10/ 11  vntyper/scripts/confidence_assignment.py
 100.0%   20/ 20  vntyper/scripts/scoring.py
 100.0%    7/  7  vntyper/scripts/variant_parsing.py
 100.0%   30/ 30  vntyper/scripts/flagging.py

Surviving mutants  [E] = hand-classified equivalent, [ ] = genuine gap
------------------------------------------------------------
vntyper/scripts/confidence_assignment.py
  [ ] line  138  '-' -> '+'

vntyper/scripts/kestrel_genotyping.py
  [ ] line  182  'and' -> 'or'
  [ ] line  220  'continue' -> 'break'
  [ ] line  284  '20' -> '21'
  [ ] line  285  '30' -> '31'
  [ ] line  286  '30' -> '31'
  [ ] line  368  'break' -> 'continue'
  [ ] line  492  'False' -> 'True'
  [ ] line  608  '0' -> '1'
  [ ] line  611  'True' -> 'False'
  [ ] line  635  'not' -> ''
  [ ] line  639  'or' -> 'and'
  [ ] line  777  'False' -> 'True'
  [ ] line  876  '0' -> '1'
  [ ] line  933  '3' -> '4'
  [ ] line  940  '0' -> '1'
  [ ] line  944  '0' -> '1'
  [ ] line  946  '0' -> '1'
  [ ] line  949  '0' -> '1'
  [ ] line  951  '0' -> '1'
  [ ] line  954  '0' -> '1'
  [ ] line  956  '0' -> '1'
  [ ] line  963  '1' -> '2'
  [ ] line 1031  'False' -> 'True'
  [ ] line 1083  '==' -> '!='
  [ ] line 1083  '1' -> '2'

vntyper/scripts/motif_decisions.py
  [ ] line   86  'True' -> 'False'

vntyper/scripts/motif_processing.py
  [ ] line  216  'False' -> 'True'
  [ ] line  231  '==' -> '!='
  [ ] line  239  'True' -> 'False'
  [ ] line  290  'False' -> 'True'
  [ ] line  292  'False' -> 'True'
  [ ] line  295  'False' -> 'True'
  [ ] line  330  'False' -> 'True'
  [ ] line  337  'True' -> 'False'
  [ ] line  341  'False' -> 'True'
  [ ] line  367  'True' -> 'False'
  [ ] line  491  'True' -> 'False'
  [ ] line  510  '-' -> '+'
  [ ] line  510  '1' -> '2'
  [ ] line  523  'not' -> ''
  [ ] line  525  'not' -> ''
  [ ] line  534  'False' -> 'True'

```
