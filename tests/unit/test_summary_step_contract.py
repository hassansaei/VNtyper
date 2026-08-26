"""The pipeline-summary step-name contract (AGENTS.md trap 5).

Five step names are matched by exact string comparison across `pipeline.py` and
`pipeline_kestrel.py` (the producers) and `generate_report.py`, `cohort_inputs.py` and
`cross_match.py` (the consumers). A typo in any one of them does not fail --
it silently drops a section from the report.

The cohort's four matches lived in `cohort_summary.py` until that file was split;
they are now in `cohort_inputs.py`, and `cohort_summary.py` remains under the
bare-literal scan so a reintroduced literal there is still caught.

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

#: The producer plus every module on the consuming side. `cohort_summary.py` stays on
#: this list although it no longer matches a name itself: it is the module a step-name
#: literal would most plausibly be reintroduced into, and the scan below is what says
#: it has none.
STEP_NAME_MODULES = (
    "pipeline.py",
    "pipeline_kestrel.py",
    "generate_report.py",
    "cohort_summary.py",
    "cohort_inputs.py",
    "cross_match.py",
)

#: How many distinct `STEP_*` constants each module must reference. These are the
#: counts of bare literals the modules carried before adoption; dropping below one
#: means a site was deleted rather than converted.
#:
#: `cohort_summary.py` is absent because its four matches moved wholesale into
#: `cohort_inputs.py` when the file was split (Task 22 of the #181-#197 follow-ups),
#: taking the count with them. Every module that matches a step name is still required
#: to reference the constants; a module that matches none is declared in
#: `MODULES_MATCHING_NO_STEP_NAME` below rather than being listed here with a zero,
#: because a zero is a passing assertion about nothing.
MINIMUM_CONSTANT_REFERENCES = {
    "pipeline.py": 4,
    "pipeline_kestrel.py": 1,
    "generate_report.py": 5,
    "cohort_inputs.py": 4,
    "cross_match.py": 2,
}

#: Modules that are scanned for bare literals but match no step name themselves, so they
#: have no minimum. Declaring them explicitly is what keeps the split above honest:
#: `test_every_scanned_module_is_classified` fails if a module is added to
#: `STEP_NAME_MODULES` without appearing in exactly one of the two, which is the
#: `KeyError` the old single-list parametrisation used to raise.
MODULES_MATCHING_NO_STEP_NAME = frozenset({"cohort_summary.py"})

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
    """Return every step name the pipeline producers can record, f-strings expanded.

    Returns:
        set[str]: The recorded names.
    """
    names: set[str] = set()
    for module in ("pipeline.py", "pipeline_kestrel.py"):
        names |= _record_step_calls(_parse(module))[1]
    return names


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


def test_every_scanned_module_is_classified() -> None:
    """Guard the split: a scanned module must be required to reference constants, or
    declared as matching no step name. Never neither.

    The two lists below used to be one, and the test that consumes the minimums used to
    raise `KeyError` for a module it had no entry for. Splitting them to avoid asserting
    a meaningless zero also removed that failure, so a module added to the scan list
    would silently get no minimum at all. This restores it.
    """
    classified = set(MINIMUM_CONSTANT_REFERENCES) | MODULES_MATCHING_NO_STEP_NAME

    assert classified == set(STEP_NAME_MODULES), (
        f"these scanned modules are unclassified: {sorted(set(STEP_NAME_MODULES) - classified)}; "
        f"these are classified but not scanned: {sorted(classified - set(STEP_NAME_MODULES))}. "
        "Add each to MINIMUM_CONSTANT_REFERENCES or to MODULES_MATCHING_NO_STEP_NAME."
    )
    assert not (set(MINIMUM_CONSTANT_REFERENCES) & MODULES_MATCHING_NO_STEP_NAME)


@pytest.mark.parametrize("module", sorted(MINIMUM_CONSTANT_REFERENCES))
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
    scanned = [_record_step_calls(_parse(module)) for module in ("pipeline.py", "pipeline_kestrel.py")]
    call_count = sum(count for count, _ in scanned)
    recorded = set().union(*(names for _, names in scanned))
    assert call_count >= 10, f"only found {call_count} pipeline record_step call sites; the scan has drifted"
    assert len(recorded) >= call_count, (
        f"resolved only {len(recorded)} names from {call_count} call sites: {sorted(recorded)}. "
        "An argument form the resolver does not understand has appeared."
    )


def test_every_consumed_step_name_is_recorded_by_the_pipeline() -> None:
    """A constant no consumer can ever match is a section that never renders."""
    missing = sorted(summary_steps.STEP_NAMES - _recorded_step_names())
    assert not missing, f"summary_steps declares step names the pipeline producers never record: {missing}."


def test_pipeline_records_only_known_step_names() -> None:
    """Every pipeline producer's record_step name must be a declared constant.

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
        f"pipeline producers record step names that no consumer knows about: {sorted(unknown)}. "
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


# ---------------------------------------------------------------------------
# Three states, not two - #212 read back
# ---------------------------------------------------------------------------


def test_a_step_that_is_not_recorded_is_absent() -> None:
    assert summary_steps.get_step_state({"steps": []}, summary_steps.STEP_KESTREL) == summary_steps.STEP_ABSENT


def test_a_step_that_produced_rows_was_read() -> None:
    summary = _summary(summary_steps.STEP_KESTREL)
    assert summary_steps.get_step_state(summary, summary_steps.STEP_KESTREL) == summary_steps.STEP_READ


def test_a_step_that_found_nothing_was_still_read() -> None:
    """The distinction this exists for: an empty result is a result.

    ``kestrel_result.tsv`` with a header and no rows is what a run that genotyped and
    called nothing writes, and it must stay distinguishable from one that wrote no file
    at all - otherwise the third state would swallow every negative report.
    """
    summary = {"steps": [{"step": summary_steps.STEP_KESTREL, "parsed_result": {"comments": [], "data": []}}]}

    assert summary_steps.get_step_state(summary, summary_steps.STEP_KESTREL) == summary_steps.STEP_READ


@pytest.mark.parametrize(
    ("step", "case"),
    [
        (
            {"step": summary_steps.STEP_KESTREL, "result_file_missing": True, "parsed_result": {"data": []}},
            "record_step's missing-file flag (#212)",
        ),
        ({"step": summary_steps.STEP_KESTREL, "parsed_result": None}, "record_step's initial value"),
        ({"step": summary_steps.STEP_KESTREL}, "no parsed_result at all"),
        (
            {"step": summary_steps.STEP_KESTREL, "parsed_result": {"error": "Error parsing file: boom"}},
            "record_step's parse-failure shape",
        ),
    ],
    ids=["missing-file", "null", "absent", "parse-error"],
)
def test_a_step_whose_result_could_not_be_read_says_so(step: dict, case: str) -> None:
    """Each signal is structural, and none of them is a message read for its wording.

    ``parse_tsv`` also records its failure as a *comment* on an otherwise ordinary
    result, and that is deliberately not matched here: recognising it would mean
    string-matching "Error parsing TSV file" against text a legitimate comment could
    carry. The missing-file flag covers the case that produces it.
    """
    assert summary_steps.get_step_state({"steps": [step]}, summary_steps.STEP_KESTREL) == (
        summary_steps.STEP_UNREADABLE
    ), case


def test_the_three_states_are_distinct_values() -> None:
    """A caller branches on these, so two of them being equal would be silent."""
    states = (summary_steps.STEP_ABSENT, summary_steps.STEP_READ, summary_steps.STEP_UNREADABLE)
    assert len(set(states)) == 3, states


# ---------------------------------------------------------------------------
# get_step_comments -- the channel #266's note travels on
# ---------------------------------------------------------------------------


class TestStepComments:
    """``summary.parse_tsv`` records a TSV's ``#`` lines here, with the hashes stripped.

    #266's below-reporting-floor note rides this channel precisely because it is *not*
    ``data``: a banner comment reaches the report without any consumer of the table
    being able to read it as a call.
    """

    def test_it_returns_the_comment_lines(self):
        summary = {
            "steps": [
                {
                    "step": summary_steps.STEP_KESTREL,
                    "parsed_result": {"comments": ["VNtyper Kestrel result", "Subthreshold candidate: x"], "data": []},
                }
            ]
        }

        assert summary_steps.get_step_comments(summary, summary_steps.STEP_KESTREL) == [
            "VNtyper Kestrel result",
            "Subthreshold candidate: x",
        ]

    def test_an_absent_step_has_no_comments(self):
        assert summary_steps.get_step_comments({"steps": []}, summary_steps.STEP_KESTREL) == []

    def test_a_step_that_recorded_none_has_no_comments(self):
        summary = {"steps": [{"step": summary_steps.STEP_KESTREL, "parsed_result": {"data": []}}]}

        assert summary_steps.get_step_comments(summary, summary_steps.STEP_KESTREL) == []

    def test_a_non_list_comments_value_is_ignored(self):
        """``parsed_result`` is read off disk, so its shape is not guaranteed."""
        summary = {"steps": [{"step": summary_steps.STEP_KESTREL, "parsed_result": {"comments": "oops"}}]}

        assert summary_steps.get_step_comments(summary, summary_steps.STEP_KESTREL) == []

    def test_an_unparsed_step_has_no_comments(self):
        summary = {"steps": [{"step": summary_steps.STEP_KESTREL, "parsed_result": None}]}

        assert summary_steps.get_step_comments(summary, summary_steps.STEP_KESTREL) == []

    def test_comments_and_data_are_separate_channels(self):
        """The property #266 depends on: nothing that reads rows can see the note."""
        summary = {
            "steps": [
                {
                    "step": summary_steps.STEP_KESTREL,
                    "parsed_result": {"comments": ["Subthreshold candidate: x"], "data": [{"Confidence": "Negative"}]},
                }
            ]
        }

        rows = summary_steps.get_step_data(summary, summary_steps.STEP_KESTREL)

        assert rows == [{"Confidence": "Negative"}]
        assert all("Subthreshold" not in str(value) for row in rows for value in row.values())
