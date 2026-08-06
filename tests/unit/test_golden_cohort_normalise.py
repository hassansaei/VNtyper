"""Unit tests for the golden-cohort gate's normalisation and its declared drops.

A normalisation is a claim that a difference does not matter, so every one has to be true
of what it is applied to. ``md5sum`` was dropped from *every* ``pipeline_summary.json``
step record on the justification that "the file's content is compared directly instead".
That was true of three of the six files this pipeline records a step for and false of the
other three - and two of those three carry decisions nothing else in the harness compares:

* ``fastq_bam_processing/pipeline_info.json`` holds the assembly guard's own verdict;
* ``advntr/cross_match_results.tsv`` holds the per-variant cross-match table.

So for those, dropping the checksum left the file compared by nothing at all.
"""

import inspect
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from golden_cohort import artifacts, normalise  # noqa: E402


def _step(result_file: str, **extra):
    """Build one ``pipeline_summary.json`` step record.

    Args:
        result_file: The absolute path the step recorded.
        **extra: Further keys to include.

    Returns:
        dict: The step record.
    """
    return {
        "step": "Some Step",
        "start": "2026-08-06T16:55:13.385444",
        "end": "2026-08-06T16:55:13.385747",
        "command": "do_something(...)",
        "result_file": result_file,
        "file_type": "tsv",
        "md5sum": "37e9ac6e1962d6c91ee2685e82e1c3a3",
        "parsed_result": {"data": [{"a": "1"}]},
        **extra,
    }


# --------------------------------------------------------------------------------------
# the conditional md5sum drop
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "result_file",
    [
        "/run/cases/x/kestrel/kestrel_result.tsv",
        "/run/cases/x/coverage/coverage_summary.tsv",
        "/run/cases/x/advntr/output_adVNTR_result.tsv",
    ],
)
def test_the_checksum_is_dropped_for_a_file_the_harness_parses_directly(result_file: str) -> None:
    """These three are read row by row by ``read_pipeline_case``, banner normalised.

    ``kestrel_result.tsv`` is the file the original justification was written for: its
    ``##`` banner carries the analysis timestamp and the VNtyper version, so its checksum
    differed on 59 of 59 cases in run 4 while its normalised content was identical.

    Args:
        result_file: A step result file with a direct comparator.
    """
    assert "md5sum" not in normalise.strip_step_record(_step(result_file))


@pytest.mark.parametrize(
    "result_file",
    [
        "/run/cases/x/fastq_bam_processing/pipeline_info.json",
        "/run/cases/x/fastq_bam_processing/output_R1.fastq.gz",
        "/run/cases/x/advntr/cross_match_results.tsv",
    ],
)
def test_the_checksum_is_kept_for_a_file_with_no_direct_comparator(result_file: str) -> None:
    """Dropping it here left the file compared by nothing at all.

    Measured rather than assumed: across run 4's 59 cases these three checksums are
    identical on both sides (59/59, 59/59 and 3/3), so restoring the comparison adds no
    delta to the recorded attestation.

    Args:
        result_file: A step result file the harness never opens.
    """
    assert normalise.strip_step_record(_step(result_file))["md5sum"]


def test_a_step_with_no_result_file_keeps_its_checksum() -> None:
    """No named file means no named comparator, so the checksum is all there is."""
    step = _step("/x/y.tsv")
    step.pop("result_file")
    assert "md5sum" in normalise.strip_step_record(step)


def test_the_directly_compared_set_is_matched_by_basename_not_by_path() -> None:
    """Every recorded path is per-side, and the rules erase the per-side prefix anyway."""
    left = normalise.strip_step_record(_step("/gate/before/cases/x/kestrel/kestrel_result.tsv"))
    right = normalise.strip_step_record(_step("/elsewhere/entirely/kestrel_result.tsv"))
    assert "md5sum" not in left
    assert "md5sum" not in right


def test_a_file_whose_name_merely_contains_a_compared_name_keeps_its_checksum() -> None:
    """The match is the whole basename; ``kestrel_result.tsv.bak`` is not that file."""
    assert "md5sum" in normalise.strip_step_record(_step("/run/kestrel/kestrel_result.tsv.bak"))


def test_the_wall_clock_keys_and_the_reparse_are_always_dropped() -> None:
    """``start``/``end`` are pure clock and ``parsed_result`` is the pipeline's own reparse."""
    stripped = normalise.strip_step_record(_step("/run/cases/x/fastq_bam_processing/pipeline_info.json"))
    assert "start" not in stripped
    assert "end" not in stripped
    assert "parsed_result" not in stripped
    assert stripped["step"] == "Some Step"
    assert stripped["command"] == "do_something(...)"


def test_stripping_does_not_mutate_the_step_it_was_given() -> None:
    """The caller iterates the parsed summary; a destructive strip would corrupt it."""
    step = _step("/run/cases/x/kestrel/kestrel_result.tsv")
    normalise.strip_step_record(step)
    assert "md5sum" in step
    assert "parsed_result" in step


def test_every_directly_compared_file_is_actually_opened_by_read_pipeline_case() -> None:
    """The exemption is only sound while the comparator it names still exists.

    If someone stops reading ``coverage_summary.tsv``, the set becomes a lie and its
    checksum goes on being discarded with nothing behind it. This is the guard against
    exactly the drift that produced the original over-broad claim.
    """
    source = inspect.getsource(artifacts.read_pipeline_case)
    for name in normalise.DIRECTLY_COMPARED_RESULT_FILES:
        assert name in source, f"{name} is exempted from the checksum comparison but is no longer read directly"


# --------------------------------------------------------------------------------------
# the manifest
# --------------------------------------------------------------------------------------


def test_the_manifest_declares_both_halves_of_the_checksum_decision() -> None:
    """A manifest that mentioned only the drop would still be the over-broad claim."""
    entries = normalise.manifest([])
    kinds = {entry["kind"] for entry in entries if entry["name"] == "md5sum"}
    assert kinds == {"dropped-key-conditional", "compared"}


def test_the_manifest_names_the_files_the_checksum_is_dropped_for() -> None:
    """A reader must be able to tell which files the exemption covers without reading code."""
    entry = next(e for e in normalise.manifest([]) if e["name"] == "md5sum" and e["kind"] == "dropped-key-conditional")
    for name in normalise.DIRECTLY_COMPARED_RESULT_FILES:
        assert name in entry["pattern"]


def test_the_manifest_no_longer_cites_a_test_that_does_not_exist() -> None:
    """The cohort-order note cited ``test_discovery_returns_an_unordered_set_today``.

    That test was renamed to ``test_the_discovered_directories_come_back_sorted`` when the
    determinism fix landed, so the citation named nothing. A justification pointing at a
    test that does not exist is worse than no justification: it reads as verified.
    """
    text = " ".join(entry["why"] for entry in normalise.manifest([]))
    assert "test_discovery_returns_an_unordered_set_today" not in text


def test_the_cohort_order_note_states_both_surviving_reasons() -> None:
    """Neither reason alone is the whole truth, and the old note stated a third that is false."""
    entry = next(e for e in normalise.manifest([]) if e["name"] == "cohort_row_order")
    assert "predating the determinism fix" in entry["why"]
    assert "ZIP" in entry["why"]


def test_the_manifest_renders_every_substitution_rule_it_was_given() -> None:
    """The point of the manifest is that no normalisation ships unlisted."""
    rules = normalise.build_rules(source_root=Path("/trees/cand"), run_root=Path("/gate/after"))
    entries = normalise.manifest(rules)
    rendered = {entry["name"] for entry in entries if entry["kind"] == "substitution"}
    assert {rule.name for rule in rules} == rendered
