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

**70 of 131 mutants killed - a raw mutation score of 53.4%.**

Of the 61 survivors, 10 are hand-classified as
*equivalent* (the mutation cannot change observable behaviour, so no test could
ever kill it) and 51 are genuine gaps. Excluding the equivalent
mutants the score is **57.9%** (70/121).

Both numbers are given because neither alone is honest: the raw score
understates the suite by counting unkillable mutants against it, and the
adjusted score depends on classifications that are a human judgement call.
Every classification is listed below with its reason so it can be checked.

| Module | Killed | Total | Raw score |
| --- | ---: | ---: | ---: |
| `vntyper/scripts/motif_processing.py` | 21 | 68 | 30.9% |
| `vntyper/scripts/variant_parsing.py` | 3 | 7 | 42.9% |
| `vntyper/scripts/confidence_assignment.py` | 11 | 19 | 57.9% |
| `vntyper/scripts/flagging.py` | 15 | 17 | 88.2% |
| `vntyper/scripts/scoring.py` | 20 | 20 | 100.0% |

## Surviving mutants

### Genuine gaps

Each of these is a change to the source that the **entire** unit tier does
not notice. They are the actionable output of this exercise: a test that
kills one is a test that would have caught a real defect of that shape.

| Module | Line | Mutation |
| --- | ---: | --- |
| `vntyper/scripts/confidence_assignment.py` | 105 | `0` &rarr; `1` |
| `vntyper/scripts/confidence_assignment.py` | 110 | `-` &rarr; `+` |
| `vntyper/scripts/motif_processing.py` | 77 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 79 | `1` &rarr; `2` |
| `vntyper/scripts/motif_processing.py` | 79 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 80 | `-` &rarr; `+` |
| `vntyper/scripts/motif_processing.py` | 80 | `1` &rarr; `2` |
| `vntyper/scripts/motif_processing.py` | 81 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 104 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 106 | `1` &rarr; `2` |
| `vntyper/scripts/motif_processing.py` | 106 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 107 | `-` &rarr; `+` |
| `vntyper/scripts/motif_processing.py` | 107 | `1` &rarr; `2` |
| `vntyper/scripts/motif_processing.py` | 108 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 189 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 204 | `==` &rarr; `!=` |
| `vntyper/scripts/motif_processing.py` | 212 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 263 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 265 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 268 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 303 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 310 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 314 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 323 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 332 | `or` &rarr; `and` |
| `vntyper/scripts/motif_processing.py` | 336 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 337 | `-` &rarr; `+` |
| `vntyper/scripts/motif_processing.py` | 337 | `1` &rarr; `2` |
| `vntyper/scripts/motif_processing.py` | 340 | `<` &rarr; `>=` |
| `vntyper/scripts/motif_processing.py` | 344 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/motif_processing.py` | 346 | `in` &rarr; `not in` |
| `vntyper/scripts/motif_processing.py` | 347 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 374 | `in` &rarr; `not in` |
| `vntyper/scripts/motif_processing.py` | 375 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 415 | `True` &rarr; `False` |
| `vntyper/scripts/motif_processing.py` | 420 | `-` &rarr; `+` |
| `vntyper/scripts/motif_processing.py` | 420 | `1` &rarr; `2` |
| `vntyper/scripts/motif_processing.py` | 424 | `>=` &rarr; `<` |
| `vntyper/scripts/motif_processing.py` | 425 | `-` &rarr; `+` |
| `vntyper/scripts/motif_processing.py` | 434 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 439 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/motif_processing.py` | 441 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/motif_processing.py` | 443 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/motif_processing.py` | 447 | `False` &rarr; `True` |
| `vntyper/scripts/motif_processing.py` | 449 | `in` &rarr; `not in` |
| `vntyper/scripts/motif_processing.py` | 452 | `in` &rarr; `not in` |
| `vntyper/scripts/motif_processing.py` | 455 | `in` &rarr; `not in` |
| `vntyper/scripts/motif_processing.py` | 459 | `True` &rarr; `False` |
| `vntyper/scripts/variant_parsing.py` | 59 | `not` &rarr; `(deleted)` |
| `vntyper/scripts/variant_parsing.py` | 59 | `and` &rarr; `or` |
| `vntyper/scripts/variant_parsing.py` | 68 | `and` &rarr; `or` |

### Classified equivalent

These mutations cannot change behaviour that any test could legitimately
observe, so they are **not** gaps in the suite. Each reason below is a claim
you can check against the source; if one turns out to be wrong the entry
should be deleted, not the score explained away.

Most of them are `.get()` defaults on `kestrel_config.json` keys that the
shipped config always supplies, which makes the default value dead code.
Being precise about the scope of that claim: a `--config-path` omitting the
key *would* reach the default, so these are unreachable **with the shipped
configuration** rather than unreachable in principle. That is the right
standard here - `AGENTS.md` trap 2 records that `--config-path` replaces the
whole config rather than merging it, and that a partial config already fails
with `KeyError` elsewhere in the pipeline, so a config missing these keys is
not a supported input.

| Module | Line | Mutation | Why it cannot be killed |
| --- | ---: | --- | --- |
| `vntyper/scripts/confidence_assignment.py` | 82 | `0` &rarr; `1` | `.get()` default for `var_active_region_threshold`; the shipped config supplies 200 |
| `vntyper/scripts/confidence_assignment.py` | 91 | `0.2` &rarr; `1.2` | `.get()` default for `depth_score_thresholds.low`; the shipped config supplies 0.00469 |
| `vntyper/scripts/confidence_assignment.py` | 92 | `0.4` &rarr; `1.4` | `.get()` default for `depth_score_thresholds.high`; the shipped config supplies 0.00515 |
| `vntyper/scripts/confidence_assignment.py` | 95 | `5` &rarr; `6` | `.get()` default for `alt_depth_thresholds.low`; the shipped config supplies 20 |
| `vntyper/scripts/confidence_assignment.py` | 96 | `10` &rarr; `11` | `.get()` default for `alt_depth_thresholds.mid_low`; the shipped config supplies 21 |
| `vntyper/scripts/confidence_assignment.py` | 97 | `20` &rarr; `21` | `.get()` default for `alt_depth_thresholds.mid_high`; the shipped config supplies 100 |
| `vntyper/scripts/flagging.py` | 150 | `False` &rarr; `True` | `.get()` default for `duplicate_flagging.enabled`; the shipped config supplies it explicitly |
| `vntyper/scripts/flagging.py` | 238 | `False` &rarr; `True` | `itertuples(index=)` only adds an `Index` field; the loop reads `row_tuple.Flag` by name |
| `vntyper/scripts/motif_processing.py` | 315 | `60` &rarr; `61` | `.get()` default for `motif_filtering.position_threshold`; the shipped config supplies 60 |
| `vntyper/scripts/variant_parsing.py` | 114 | `0.0` &rarr; `1.0` | `.get()` default for `alt_filtering.gg_depth_score_threshold`; the shipped config supplies 0.00469 |

## How this compares to the 43.5% baseline

The experiment that motivated this work scored **43.5%** (27 of 62 mutants
killed) across the eight highest-coverage modules. That harness was never
committed, which is why this one exists.

!!! note "The two totals are not directly comparable"

    Different mutant population, different modules: 62 mutants over eight
    modules then, 131 over five now, generated by a different operator set.
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
| Mutation score | 21% | 57.9% raw, 84.6% adjusted |

## Reproducing this

```bash
make mutation
```

The harness is `scripts/mutation_test.py`. It mutates one token at a time,
runs the module's own tests first and escalates anything that survives to
the full unit tier before recording it as a survivor, so the score is not
biased by the scoping.

!!! danger "`PYTHONDONTWRITEBYTECODE=1` is load-bearing"

    Most mutations here are byte-length preserving (`<` to `>`, `and` to
    `or`), and CPython validates a cached `.pyc` against the source's
    `(mtime, size)` pair with one-second mtime granularity. A mutant written
    in the same second as the file it replaces is therefore
    indistinguishable from it to the cache validator: the interpreter loads
    the stale `.pyc`, runs the **unmutated** code, every mutant *survives*
    and the score is fiction. Two sweeps produced exactly that before it was
    found. The harness also deletes every `__pycache__` before it starts -
    the environment variable stops new caches, the deletion stops old ones,
    and both are required.

!!! danger "Nothing may build or install from the tree while a sweep runs"

    The harness rewrites `vntyper/scripts/*.py` **in place**, so for most of
    a run the working tree holds a deliberately broken module. Anything
    that snapshots source mid-sweep bakes that mutant into its artefact -
    a docker build, `pip install`, `python -m build`, a tarball.

    This has happened: an image built during a sweep crashed in the
    container at `motif_processing.py` with a pandas `KeyError`, which
    reads exactly like a production bug and cost a full diagnosis cycle
    before it was traced back to the sweep. Rebuilding from a clean tree
    passed.

    The `finally` restore protects the **repository**, not any artefact
    already produced from it. Run `git diff --quiet -- vntyper/`
    immediately before and after any build, package or install step; if it
    reports a difference you did not make, a sweep is running and the
    artefact is void.

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
Total:    131 mutants, 70 killed, 61 survived
Score:    53.4%
Duration: 9.2 min

Per module
------------------------------------------------------------
  30.9%   21/ 68  vntyper/scripts/motif_processing.py
  42.9%    3/  7  vntyper/scripts/variant_parsing.py
  57.9%   11/ 19  vntyper/scripts/confidence_assignment.py
  88.2%   15/ 17  vntyper/scripts/flagging.py
 100.0%   20/ 20  vntyper/scripts/scoring.py

Surviving mutants  [E] = hand-classified equivalent, [ ] = genuine gap
------------------------------------------------------------
vntyper/scripts/confidence_assignment.py
  [E] line   82  '0' -> '1'
          equivalent: `.get()` default for `var_active_region_threshold`; the shipped config supplies 200
  [E] line   91  '0.2' -> '1.2'
          equivalent: `.get()` default for `depth_score_thresholds.low`; the shipped config supplies 0.00469
  [E] line   92  '0.4' -> '1.4'
          equivalent: `.get()` default for `depth_score_thresholds.high`; the shipped config supplies 0.00515
  [E] line   95  '5' -> '6'
          equivalent: `.get()` default for `alt_depth_thresholds.low`; the shipped config supplies 20
  [E] line   96  '10' -> '11'
          equivalent: `.get()` default for `alt_depth_thresholds.mid_low`; the shipped config supplies 21
  [E] line   97  '20' -> '21'
          equivalent: `.get()` default for `alt_depth_thresholds.mid_high`; the shipped config supplies 100
  [ ] line  105  '0' -> '1'
  [ ] line  110  '-' -> '+'

vntyper/scripts/flagging.py
  [E] line  150  'False' -> 'True'
          equivalent: `.get()` default for `duplicate_flagging.enabled`; the shipped config supplies it explicitly
  [E] line  238  'False' -> 'True'
          equivalent: `itertuples(index=)` only adds an `Index` field; the loop reads `row_tuple.Flag` by name

vntyper/scripts/motif_processing.py
  [ ] line   77  'True' -> 'False'
  [ ] line   79  '1' -> '2'
  [ ] line   79  'True' -> 'False'
  [ ] line   80  '-' -> '+'
  [ ] line   80  '1' -> '2'
  [ ] line   81  'True' -> 'False'
  [ ] line  104  'True' -> 'False'
  [ ] line  106  '1' -> '2'
  [ ] line  106  'True' -> 'False'
  [ ] line  107  '-' -> '+'
  [ ] line  107  '1' -> '2'
  [ ] line  108  'True' -> 'False'
  [ ] line  189  'False' -> 'True'
  [ ] line  204  '==' -> '!='
  [ ] line  212  'True' -> 'False'
  [ ] line  263  'False' -> 'True'
  [ ] line  265  'False' -> 'True'
  [ ] line  268  'False' -> 'True'
  [ ] line  303  'False' -> 'True'
  [ ] line  310  'True' -> 'False'
  [ ] line  314  'False' -> 'True'
  [E] line  315  '60' -> '61'
          equivalent: `.get()` default for `motif_filtering.position_threshold`; the shipped config supplies 60
  [ ] line  323  'True' -> 'False'
  [ ] line  332  'or' -> 'and'
  [ ] line  336  'True' -> 'False'
  [ ] line  337  '-' -> '+'
  [ ] line  337  '1' -> '2'
  [ ] line  340  '<' -> '>='
  [ ] line  344  'not' -> ''
  [ ] line  346  'in' -> 'not in'
  [ ] line  347  'True' -> 'False'
  [ ] line  374  'in' -> 'not in'
  [ ] line  375  'True' -> 'False'
  [ ] line  415  'True' -> 'False'
  [ ] line  420  '-' -> '+'
  [ ] line  420  '1' -> '2'
  [ ] line  424  '>=' -> '<'
  [ ] line  425  '-' -> '+'
  [ ] line  434  'False' -> 'True'
  [ ] line  439  'not' -> ''
  [ ] line  441  'not' -> ''
  [ ] line  443  'not' -> ''
  [ ] line  447  'False' -> 'True'
  [ ] line  449  'in' -> 'not in'
  [ ] line  452  'in' -> 'not in'
  [ ] line  455  'in' -> 'not in'
  [ ] line  459  'True' -> 'False'

vntyper/scripts/variant_parsing.py
  [ ] line   59  'not' -> ''
  [ ] line   59  'and' -> 'or'
  [ ] line   68  'and' -> 'or'
  [E] line  114  '0.0' -> '1.0'
          equivalent: `.get()` default for `alt_filtering.gg_depth_score_threshold`; the shipped config supplies 0.00469

```
