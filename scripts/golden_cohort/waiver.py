"""Which deltas fail the gate, and which one a caller may declare away.

Every delta this harness finds is fatal. That is the whole design: a gate whose failure
set narrows quietly reports a pass it did not measure. The single exception lives here,
alone, because it is the only rule in the harness that can turn a check *off* - and a
policy that can turn checks off should be readable in one screen rather than found by
grep in the middle of the diffing code.

The rule has two halves and they must agree. :func:`fatal_deltas` decides the verdict;
:func:`waived_cases` decides what the rendered report says was waived. If the second ever
stopped tracking the first, a run could reach ``IDENTICAL`` by waiving a delta and never
name the case it waived - which is precisely the silent narrowing the rest of the harness
is built to prevent.
"""

from __future__ import annotations

from typing import Any

#: The one artefact whose delta a caller may declare in advance.
#:
#: Performance work changes the executed command stream by construction, so without a
#: declaration the gate can never report anything but ``DELTAS`` for it - and a gate that
#: always fails is one nobody reads. The exemption is exactly one artefact wide and it is
#: an allowlist rather than a "genotype-bearing only" rule, because ``report_tables``,
#: ``pipeline_steps``, ``pipeline_step_records``, the CRAM read-set evidence and the
#: ``EXPECTATION`` delta are all deliberate gate failures that such a rule would silently
#: turn off.
DECLARABLE_DELTA = "executed_commands"

#: The rendered report's heading for the cases whose command delta was waived.
#:
#: A distinct constant rather than prose, because a run that reaches ``IDENTICAL`` by
#: waiving a delta attests strictly less than one that reaches it outright - the same
#: distinction ``REDUCED`` exists to make - and the reader must be told which cases used
#: the waiver so they can go and read the diff.
DECLARED_DELTA_HEADING = "## Declared command-stream deltas, waived by --expect-command-delta"


def fatal_deltas(case: dict[str, Any], *, expect_command_delta: bool) -> list[str]:
    """Return the deltas that must fail the gate.

    Every delta is fatal except :data:`DECLARABLE_DELTA`, and that one only when the
    caller declares in advance that the command stream was changed on purpose.

    The delta itself is never removed from ``case["deltas"]``: the report must always
    name what changed, because the declaration is only worth anything if a human then
    reads the diff and confirms it is the change they intended. Only the verdict, and
    therefore the exit status, is affected.

    Args:
        case: One case's ``diff_case`` result, after the expectation fold.
        expect_command_delta: Whether a command-stream change was declared.

    Returns:
        list[str]: The fatal delta names, in the order ``diff_case`` reported them.
    """
    exempt = {DECLARABLE_DELTA} if expect_command_delta else set()
    return [name for name in case.get("deltas", []) if name not in exempt]


def waived_cases(result: dict[str, Any]) -> list[str]:
    """Return the cases that reported a command delta the declaration waived.

    Derived from the same two lists :func:`fatal_deltas` produced rather than from the
    caller's flag, so the report can only ever name a case whose delta was genuinely
    excused from the verdict.

    Args:
        result: The comparison document, which must carry ``cases``.

    Returns:
        list[str]: Case ids, sorted, so that two runs' reports diff cleanly.
    """
    return sorted(
        case_id
        for case_id, case in result["cases"].items()
        if DECLARABLE_DELTA in case.get("deltas", []) and DECLARABLE_DELTA not in case.get("fatal_deltas", [])
    )
