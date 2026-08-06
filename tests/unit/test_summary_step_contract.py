"""The pipeline-summary step-name contract (AGENTS.md trap 5).

Five step names are matched by exact string comparison across `pipeline.py`
(the producer) and `generate_report.py`, `cohort_summary.py` and
`cross_match.py` (the consumers). A typo in any one of them does not fail --
it silently drops a section from the report.

`summary_steps.py` is the single source of truth for those names. These tests
enforce three things:

1. The constants' values, pinned individually and byte for byte. Everything else
   here is a claim about which name is used *where*; this is the claim about
   *what* the names are, and adopting the constants must not have changed one.
2. No bare step-name literal survives in any of the four modules. The producer
   scan alone was not enough: a consumer comparing against a mistyped literal
   still drops a section, and nothing checked the consumers at all.
3. Every name the producer can record is a declared constant, and every declared
   constant is a name the producer records.

Both source scans are AST-based rather than regex-based, and both are guarded by
a positive control: a scan that matches nothing makes every assertion built on
it pass vacuously.
"""

import ast
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

from vntyper.scripts import summary_steps  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "vntyper" / "scripts"

#: The producer plus every consumer that matches a step name.
STEP_NAME_MODULES = ("pipeline.py", "generate_report.py", "cohort_summary.py", "cross_match.py")

#: How many distinct `STEP_*` constants each module must reference. These are the
#: counts of bare literals the modules carried before adoption; dropping below one
#: means a site was deleted rather than converted.
MINIMUM_CONSTANT_REFERENCES = {
    "pipeline.py": 5,
    "generate_report.py": 5,
    "cohort_summary.py": 4,
    "cross_match.py": 2,
}

# One producer call site builds its name from an f-string, guarded by
# `if input_type in ["BAM", "CRAM"]`. Expand it so the tests reason about the names
# that actually reach a pipeline_summary.json.
_DYNAMIC_EXPANSIONS = {
    "{input_type} to FASTQ Conversion": {
        "BAM to FASTQ Conversion",
        "CRAM to FASTQ Conversion",
    },
}


def _summary(*steps: str) -> dict:
    return {
        "steps": [
            {"step": name, "parsed_result": {"comments": [], "data": [{"n": str(i)}]}} for i, name in enumerate(steps)
        ]
    }


def _parse(module: str) -> ast.Module:
    """Parse one module under ``vntyper/scripts``.

    Args:
        module: The file name, e.g. ``pipeline.py``.

    Returns:
        ast.Module: The parsed tree.
    """
    return ast.parse((SCRIPTS / module).read_text(encoding="utf-8"))


def _bare_step_literals(tree: ast.Module) -> list[tuple[int, str]]:
    """Find string literals exactly equal to a declared step name.

    Comments are invisible to the AST and prose mentions never compare equal, so
    this finds load-bearing literals only.

    Args:
        tree: A parsed module.

    Returns:
        list[tuple[int, str]]: ``(line number, literal)`` for each hit.
    """
    return [
        (node.lineno, node.value)
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, str) and node.value in summary_steps.STEP_NAMES
    ]


def _step_constant_references(tree: ast.Module) -> set[str]:
    """Find every ``STEP_*`` constant a module refers to.

    Args:
        tree: A parsed module.

    Returns:
        set[str]: The referenced constant names.
    """
    declared = {name for name in vars(summary_steps) if name.startswith("STEP_")}
    referenced = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute) and node.attr in declared:
            referenced.add(node.attr)
        elif isinstance(node, ast.Name) and node.id in declared:
            referenced.add(node.id)
    return referenced - {"STEP_NAMES"}


def _resolve_step_name(node: ast.expr) -> set[str]:
    """Resolve a ``record_step`` name argument to the names it can produce.

    Args:
        node: The second positional argument of a ``record_step`` call.

    Returns:
        set[str]: The names, f-string templates expanded.
    """
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return _DYNAMIC_EXPANSIONS.get(node.value, {node.value})

    if isinstance(node, ast.Name) and hasattr(summary_steps, node.id):
        return {getattr(summary_steps, node.id)}

    if isinstance(node, ast.Attribute) and hasattr(summary_steps, node.attr):
        return {getattr(summary_steps, node.attr)}

    if isinstance(node, ast.JoinedStr):
        template = ""
        for part in node.values:
            if isinstance(part, ast.Constant):
                template += str(part.value)
            elif isinstance(part, ast.FormattedValue):
                template += "{" + ast.unparse(part.value) + "}"
        return _DYNAMIC_EXPANSIONS.get(template, {template})

    return set()


def _record_step_calls(tree: ast.Module) -> tuple[int, set[str]]:
    """Find every ``record_step(summary, <name>, ...)`` call and its name.

    Args:
        tree: A parsed module.

    Returns:
        tuple[int, set[str]]: How many call sites were found, and the set of
        names they can record.
    """
    call_count = 0
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "record_step":
            call_count += 1
            if len(node.args) >= 2:
                names |= _resolve_step_name(node.args[1])
    return call_count, names


def _recorded_step_names() -> set[str]:
    """Return every step name pipeline.py can record, f-strings expanded.

    Returns:
        set[str]: The recorded names.
    """
    return _record_step_calls(_parse("pipeline.py"))[1]


# ---------------------------------------------------------------------------
# 1. The constants themselves, pinned byte for byte
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "constant,value",
    [
        ("STEP_BAM_HEADER", "BAM Header Parsing"),
        ("STEP_COVERAGE", "Coverage Calculation"),
        ("STEP_KESTREL", "Kestrel Genotyping"),
        ("STEP_ADVNTR", "adVNTR Genotyping"),
        ("STEP_CROSS_MATCH", "Cross-Match Variant Comparison"),
    ],
)
def test_each_step_constant_has_its_exact_value(constant: str, value: str) -> None:
    """Adopting the constants must not have changed a single character.

    The set assertion below cannot catch two constants swapped with each other,
    and a swap would silently move a report section. This pins each one.
    """
    assert getattr(summary_steps, constant) == value


def test_step_names_are_exactly_the_five_matched_literals() -> None:
    expected = frozenset(
        {
            "BAM Header Parsing",
            "Coverage Calculation",
            "Kestrel Genotyping",
            "adVNTR Genotyping",
            "Cross-Match Variant Comparison",
        }
    )
    assert expected == summary_steps.STEP_NAMES


# ---------------------------------------------------------------------------
# 2. No bare literal survives, producer or consumer
# ---------------------------------------------------------------------------


def test_the_literal_scan_actually_finds_literals() -> None:
    """Guard the guard: a scan that matches nothing makes the next test vacuous.

    ``summary_steps.py`` is the one place the literals legitimately appear, so it
    is a positive control that cannot rot.
    """
    hits = _bare_step_literals(_parse("summary_steps.py"))
    found = {value for _, value in hits}
    assert found == set(summary_steps.STEP_NAMES), (
        f"the literal scanner found {sorted(found)} in summary_steps.py, expected all five; "
        "the scanner has drifted and every assertion built on it is vacuous"
    )


def test_the_literal_scan_ignores_prose_and_comments() -> None:
    """Docstrings and comments name the steps on purpose; only code counts."""
    source = (
        '"""The Kestrel Genotyping step feeds Cross-Match Variant Comparison."""\n'
        "# adVNTR Genotyping is optional.\n"
        'x = "Coverage Calculation"\n'
    )
    assert _bare_step_literals(ast.parse(source)) == [(3, "Coverage Calculation")]


@pytest.mark.parametrize("module", STEP_NAME_MODULES)
def test_no_bare_step_name_literal_survives(module: str) -> None:
    """A literal is a typo waiting to drop a report section silently."""
    hits = _bare_step_literals(_parse(module))
    assert not hits, (
        f"{module} still compares against bare step-name literals: "
        + ", ".join(f"line {line}: {value!r}" for line, value in hits)
        + ". Use the summary_steps.STEP_* constants instead."
    )


@pytest.mark.parametrize("module", STEP_NAME_MODULES)
def test_each_module_actually_references_the_constants(module: str) -> None:
    """Deleting a site also removes its literal; this catches that."""
    referenced = _step_constant_references(_parse(module))
    minimum = MINIMUM_CONSTANT_REFERENCES[module]
    assert len(referenced) >= minimum, (
        f"{module} refers to only {len(referenced)} STEP_* constants ({sorted(referenced)}), "
        f"expected at least {minimum}. A step site was removed rather than converted."
    )
    assert referenced <= {name for name in vars(summary_steps) if name.startswith("STEP_")}


# ---------------------------------------------------------------------------
# 3. Producer and consumers agree on the same set of names
# ---------------------------------------------------------------------------


def test_the_record_step_scan_actually_finds_call_sites() -> None:
    """Guard the guard: a scan that matches nothing makes the next test vacuous."""
    call_count, recorded = _record_step_calls(_parse("pipeline.py"))
    assert call_count >= 10, f"only found {call_count} record_step call sites in pipeline.py; the scan has drifted"
    assert len(recorded) >= call_count, (
        f"resolved only {len(recorded)} names from {call_count} call sites: {sorted(recorded)}. "
        "An argument form the resolver does not understand has appeared."
    )


def test_every_consumed_step_name_is_recorded_by_the_pipeline() -> None:
    """A constant no consumer can ever match is a section that never renders."""
    missing = sorted(summary_steps.STEP_NAMES - _recorded_step_names())
    assert not missing, f"summary_steps declares step names pipeline.py never records: {missing}."


def test_pipeline_records_only_known_step_names() -> None:
    """Every record_step name in pipeline.py must be a declared constant.

    Catches a step renamed at the producer without updating the consumers.
    """
    # Steps that exist but no consumer reads; allowed, but they must be declared here
    informational = {
        "SHARK Filtering",
        "FASTQ Quality Control",
        "BAM to FASTQ Conversion",
        "BAM to FASTQ Conversion (Post-alignment)",
        "CRAM to FASTQ Conversion",
        "FASTQ Alignment",
    }
    unknown = _recorded_step_names() - summary_steps.STEP_NAMES - informational
    assert not unknown, (
        f"pipeline.py records step names that no consumer knows about: {sorted(unknown)}. "
        "Add them to summary_steps.STEP_NAMES or to the informational set in this test."
    )


# ---------------------------------------------------------------------------
# The accessors
# ---------------------------------------------------------------------------


def test_get_step_data_returns_rows_for_a_present_step() -> None:
    summary = _summary(summary_steps.STEP_KESTREL, summary_steps.STEP_COVERAGE)
    assert summary_steps.get_step_data(summary, summary_steps.STEP_KESTREL) == [{"n": "0"}]


def test_get_step_data_returns_empty_for_an_absent_step() -> None:
    summary = _summary(summary_steps.STEP_KESTREL)
    assert summary_steps.get_step_data(summary, summary_steps.STEP_ADVNTR) == []


def test_get_step_data_tolerates_a_non_data_parsed_result() -> None:
    """BAM Header Parsing yields a flat dict, not {"data": [...]}"""
    summary = {"steps": [{"step": summary_steps.STEP_BAM_HEADER, "parsed_result": {"assembly_text": "hg38"}}]}
    assert summary_steps.get_step_data(summary, summary_steps.STEP_BAM_HEADER) == []
    assert summary_steps.get_step_result(summary, summary_steps.STEP_BAM_HEADER) == {"assembly_text": "hg38"}


def test_get_step_handles_a_summary_with_no_steps_key() -> None:
    assert summary_steps.get_step({}, summary_steps.STEP_KESTREL) is None
