"""The pipeline-summary step-name contract (AGENTS.md trap 5)."""

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

from vntyper.scripts import summary_steps  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]

# pipeline.py calls record_step(summary, "<name>", ...) with the name as the second
# argument, on its own line. Matching only `record_step("` finds nothing at all, which
# would make test_pipeline_records_only_known_step_names silently vacuous.
_RECORD_STEP = re.compile(r'record_step\(\s*summary,\s*f?"([^"]+)"')

# One call site builds its name from an f-string, guarded by
# `if input_type in ["BAM", "CRAM"]`. Expand it so the test reasons about the names
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


def _recorded_step_names() -> set[str]:
    """Return every step name pipeline.py can record, f-strings expanded.

    Returns:
        set[str]: The recorded names.
    """
    source = (REPO_ROOT / "vntyper" / "scripts" / "pipeline.py").read_text(encoding="utf-8")
    names: set[str] = set()
    for literal in _RECORD_STEP.findall(source):
        names |= _DYNAMIC_EXPANSIONS.get(literal, {literal})
    return names


def test_step_names_are_exactly_the_five_matched_literals() -> None:
    assert summary_steps.STEP_NAMES == frozenset(
        {
            "BAM Header Parsing",
            "Coverage Calculation",
            "Kestrel Genotyping",
            "adVNTR Genotyping",
            "Cross-Match Variant Comparison",
        }
    )


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


def test_the_record_step_scan_actually_finds_call_sites() -> None:
    """Guard the guard: a regex that matches nothing makes the next test vacuous."""
    recorded = _recorded_step_names()
    assert len(recorded) >= 10, f"Only found {sorted(recorded)}; the record_step regex has drifted from pipeline.py."


def test_every_consumed_step_name_is_recorded_by_the_pipeline() -> None:
    """A constant no consumer can ever match is a section that never renders."""
    missing = sorted(summary_steps.STEP_NAMES - _recorded_step_names())
    assert not missing, f"summary_steps declares step names pipeline.py never records: {missing}."


def test_pipeline_records_only_known_step_names() -> None:
    """Every record_step literal in pipeline.py must be a declared constant.

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
