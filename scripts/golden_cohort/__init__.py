"""The golden-cohort gate harness (#179).

``docs/development/golden-cohort-gate.md`` describes three runs of a before-versus-after
comparison over the whole local test cohort. The scripts that produced those runs lived
in ``/tmp`` and are gone, so the instrument backing every genotype claim on this project
had to be rebuilt from prose each time, by a different person, with no guarantee it was
the same instrument. This package **is** that instrument, committed.

Entry point: ``scripts/golden_cohort_gate.py``. Nothing here is imported by
``vntyper``; the harness only ever drives it as a subprocess, through
:mod:`golden_cohort.launcher`, so that which package a run resolved is a recorded fact
rather than an assumption.

Modules:
    matrix: Derive the case matrix from what is actually in ``tests/data``.
    launcher: The in-process wrapper that proves which tree a run executed.
    runner: Execute one side of the comparison.
    artifacts: Read one run's comparable artefacts off disk.
    normalise: The substitutions applied to everything that differs by construction.
    compare: Diff two sides and report.
"""

from __future__ import annotations

#: Bumped when the harness changes what it measures or how it normalises, so a written-up
#: run can name the instrument version as well as the two commits.
HARNESS_VERSION = "1.0.0"
