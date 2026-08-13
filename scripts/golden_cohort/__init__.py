"""The golden-cohort gate harness (#179).

``docs/development/golden-cohort-gate.md`` describes a series of before-versus-after
comparisons over the whole local test cohort - one section per run, so this docstring
deliberately does not say how many there are; the sentence it replaces said "three" and was
wrong from the day run 4 was written up. The scripts that produced the first three lived in
``/tmp`` and are gone, so the instrument backing every genotype claim on this project had
to be rebuilt from prose each time, by a different person, with no guarantee it was the
same instrument. This package **is** that instrument, committed.

Entry point: ``scripts/golden_cohort_gate.py``. Nothing here is imported by
``vntyper``; the harness only ever drives it as a subprocess, through
:mod:`golden_cohort.launcher`, so that which package a run resolved is a recorded fact
rather than an assumption.

Modules:
    matrix: Derive the base cases from ``tests/data`` and apply the declared policies.
    launcher: The in-process wrapper that proves which tree a run executed.
    admissibility: The checks that decide whether a run is admissible as evidence at all -
        per-case exit and artefact expectations, side opposition, and revision recording.
    runner: Execute one side of the comparison.
    artifacts: Read one run's comparable artefacts off disk.
    normalise: The substitutions applied to everything that differs by construction.
    waiver: Which deltas are fatal, and the single one a caller may declare away.
    compare: Diff two sides and report.
"""

from __future__ import annotations

#: Bumped when the harness changes what it measures or how it normalises, so a written-up
#: run can name the instrument version as well as the two commits.
#:
#: 1.1.0 enforces each case's declared exit code and required artefacts, refuses a
#: comparison of two sides that are not opposed, records each side's git revision and
#: working-tree state, refuses a drifted or empty matrix, keeps the ``md5sum`` of step
#: result files that have no direct comparator, and folds a changed provenance banner into
#: a table's status. Runs recorded under 1.0.0 measured strictly less; the four runs on the
#: gate page were all produced by 1.0.0.
#:
#: 1.3.0 adds the parsed placed-unmapped guard count to the comparison surface, so a
#: forced-indexed exit remains causally attributable in the rendered gate result.
#:
#: 1.4.0 requires the observed extraction mode to match the declared forced mode and
#: makes an unsafe forced-indexed attestation prove its causal exit-before-work shape.
#:
#: 1.5.0 compares ``output.bed``, which the collector never read, and adds the declared
#: command-delta mode: a caller may state in advance that the executed command stream
#: changes on purpose, which keeps that one delta out of the verdict while leaving every
#: other delta - genotype, report, step, evidence and expectation - fatal. A run recorded
#: under 1.5.0 with the declaration in force attests strictly less about the command
#: stream than one without it, so the result document and the rendered report both say
#: which cases used it.
HARNESS_VERSION = "1.5.0"
