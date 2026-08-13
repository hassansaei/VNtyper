"""How the report names the sample it is about, and what it says when it cannot.

Every report a lab produces used to be titled "Summary Report" -- the same
constant string in ``<title>`` and ``<h1>`` for every run -- so two reports open
in two browser tabs were indistinguishable, and a printed one carried no
identity at all (#242).

The rule these tests pin is deliberately conservative. An explicit
``--sample-name`` always wins. Otherwise a name is derived from the input
basename by stripping a *recognised* compound extension and then a single
trailing ``_R1``/``_R2``, and **anything the rule does not recognise is printed
verbatim** rather than guessed at: a wrong sample name on a report is worse than
an ugly one.

The four provenance helpers share the same discipline in the other direction: a
value the run did not record renders :data:`NOT_RECORDED`, never a configured
default. Falling back to ``config["default_values"]["reference_assembly"]``
would mislabel any ``--reference-assembly`` override and cannot reconstruct
``--custom-regions`` at all.
"""

import logging
from pathlib import Path

import pytest

from vntyper.scripts.report_identity import (
    NOT_RECORDED,
    UNNAMED_SAMPLE,
    derive_sample_name,
    format_input_files,
    format_region,
    format_run_timestamp,
    input_file_names,
    recorded_or_not,
    recorded_sample_name,
    resolve_sample_name,
)

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# The sample-name fallback rule
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "basename,expected",
    [
        ("example_b178_hg19_subset_R1.fastq.gz", "example_b178_hg19_subset"),
        ("example_b178_hg19_subset_R2.fastq.gz", "example_b178_hg19_subset"),
        ("example_b178_hg19_bwa.cram", "example_b178_hg19_bwa"),
        ("S1.lane3.L001_R1.fastq.gz", "S1.lane3.L001"),
    ],
)
def test_sample_name_fallback(basename, expected) -> None:
    """The four cases the plan states, verbatim.

    The third and fourth are the ones that make a naive implementation wrong:
    stripping two suffixes unconditionally turns ``example_b178_hg19_bwa.cram``
    into ``example_b178_hg19``, while refusing any name with more than one dot
    turns ``S1.lane3.L001_R1.fastq.gz`` into the whole file name.
    """
    assert derive_sample_name(basename) == expected


def test_an_explicit_sample_name_overrides_the_fallback() -> None:
    assert resolve_sample_name("--given--", Path("x_R1.fq.gz")) == "--given--"


@pytest.mark.parametrize(
    "basename",
    [
        "sample.txt",  # not an input format this pipeline reads
        "no_extension_at_all",
        "archive.fastq.gz.bak",
    ],
)
def test_an_unrecognised_extension_is_printed_verbatim(basename) -> None:
    """The rule strips only extensions it recognises, and otherwise stops."""
    assert derive_sample_name(basename) == basename


def test_a_mate_marker_that_is_not_trailing_is_left_alone() -> None:
    """The shape a looser rule gets wrong.

    ``S1_S1_L001_R1_001.fastq.gz`` has a recognised extension, so that much is
    stripped -- but an implementation that searched for ``_R1`` *anywhere* would
    go on to name the sample ``S1_S1_L001_001``, which is not what any other
    artefact of the run calls it. Only a trailing marker is a mate marker.
    """
    assert derive_sample_name("S1_S1_L001_R1_001.fastq.gz") == "S1_S1_L001_R1_001"


def test_a_directory_component_is_dropped() -> None:
    """``input_files`` holds basenames, but a caller may pass a path."""
    assert derive_sample_name("/runs/2026/S1_R1.fastq.gz") == "S1"


def test_the_extension_match_is_case_insensitive() -> None:
    assert derive_sample_name("S1_R1.FASTQ.GZ") == "S1"


def test_a_name_that_is_only_an_extension_stays_verbatim() -> None:
    """Stripping would leave nothing, and an empty title names no one."""
    assert derive_sample_name(".cram") == ".cram"


def test_a_stem_that_is_only_a_mate_marker_stays_verbatim() -> None:
    assert derive_sample_name("_R1.fastq.gz") == "_R1"


def test_an_empty_basename_is_named_as_unnamed() -> None:
    assert derive_sample_name("   ") == UNNAMED_SAMPLE


def test_agreeing_mates_resolve_to_the_shared_stem() -> None:
    assert resolve_sample_name(None, "S1_R1.fastq.gz", "S1_R2.fastq.gz") == "S1"


def test_disagreeing_mates_fall_back_to_the_first_basename_verbatim() -> None:
    """ "Mate stems that disagree" is the third ambiguous case, and the honest
    answer is the file name -- not one of the two names, presented as if the
    pipeline had decided between them."""
    assert resolve_sample_name(None, "a.fq", "b.fq") == "a.fq"


def test_no_inputs_at_all_resolve_to_unnamed() -> None:
    assert resolve_sample_name(None) == UNNAMED_SAMPLE


def test_a_blank_explicit_name_does_not_win() -> None:
    """An empty ``--sample-name`` is not a name; it must not blank the title."""
    assert resolve_sample_name("   ", "S1_R1.fastq.gz") == "S1"


def test_an_explicit_name_is_stripped_but_otherwise_untouched() -> None:
    assert resolve_sample_name("  PATIENT 042  ") == "PATIENT 042"


# ---------------------------------------------------------------------------
# The three levels of precedence
# ---------------------------------------------------------------------------


def test_the_run_s_own_name_beats_the_basename_derivation() -> None:
    """Level 2. `vntyper pipeline -s PATIENT_042 --bam foo.bam` embedded
    `PATIENT_042` in Kestrel's outputs and its VCF header while the report was
    titled `foo` -- one run, two identities, in the artefact people forward."""
    assert resolve_sample_name(None, "foo.bam", recorded="PATIENT_042") == "PATIENT_042"


def test_an_explicit_name_beats_the_run_s_own() -> None:
    """Level 1. Someone re-reporting an archived run gets the last word."""
    assert resolve_sample_name("RENAMED", "foo.bam", recorded="PATIENT_042") == "RENAMED"


@pytest.mark.parametrize("recorded", [None, "", "   "])
def test_a_summary_with_no_recorded_name_falls_through_to_the_derivation(recorded) -> None:
    """Level 3, and the regression risk that matters: every summary written
    before the field existed has no level 2 at all and must resolve exactly as
    it did before."""
    assert resolve_sample_name(None, "S1_R1.fastq.gz", recorded=recorded) == "S1"


def test_the_placeholder_name_does_not_win(caplog) -> None:
    """`cli_handlers` falls back to the literal `"sample"` when it can resolve
    nothing, and the run embeds that. It is the absence of a name, not a name,
    so it must not displace a basename that names the sample perfectly well."""
    with caplog.at_level(logging.INFO, logger="vntyper.scripts.report_identity"):
        assert resolve_sample_name(None, "S1_R1.fastq.gz", recorded="sample") == "S1"

    assert "placeholder" in caplog.text


def test_an_explicit_name_still_wins_over_the_placeholder_rule() -> None:
    """An explicit value at render time is unambiguous in a way a recorded one is
    not, so `vntyper report --sample-name sample` is honoured."""
    assert resolve_sample_name("sample", "S1_R1.fastq.gz", recorded="PATIENT_042") == "sample"


def test_the_run_s_own_name_wins_even_when_the_mates_disagree() -> None:
    """The disagreeing-mate fallback is level 3; it never overrides level 2."""
    assert resolve_sample_name(None, "a.fq", "b.fq", recorded="PATIENT_042") == "PATIENT_042"


def test_the_run_s_own_name_is_used_when_there_are_no_input_files_at_all() -> None:
    assert resolve_sample_name(None, recorded="PATIENT_042") == "PATIENT_042"


@pytest.mark.parametrize(
    "recorded,expected",
    [
        ("PATIENT_042", "PATIENT_042"),
        ("  PATIENT_042  ", "PATIENT_042"),
        ("sample", None),
        ("", None),
        ("   ", None),
        (None, None),
    ],
)
def test_the_recorded_name_predicate(recorded, expected) -> None:
    assert recorded_sample_name(recorded) == expected


# ---------------------------------------------------------------------------
# The input-file line, which must not branch on the shape
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "input_files,expected",
    [
        ({"fastq1": "S1.fq.gz"}, "S1.fq.gz"),
        ({"fastq1": "a.fq", "fastq2": "b.fq"}, "a.fq, b.fq"),
        ({"bam": "S1.bam"}, "S1.bam"),
        ({"cram": "S1.cram"}, "S1.cram"),
        # A shape nobody has written yet. The template used to branch on two of
        # the four `resolve_pipeline_input` produces, so CRAM and single-end
        # FASTQ rendered an empty line; iterating means a fifth cannot.
        ({"ubam": "S1.ubam"}, "S1.ubam"),
    ],
)
def test_every_input_shape_is_listed(input_files, expected) -> None:
    assert format_input_files(input_files) == expected


def test_no_input_files_says_so() -> None:
    assert format_input_files({}) == NOT_RECORDED


def test_a_blank_input_file_entry_is_dropped() -> None:
    assert format_input_files({"bam": "", "cram": "S1.cram"}) == "S1.cram"


def test_the_names_are_returned_in_the_order_the_summary_records_them() -> None:
    assert input_file_names({"fastq1": "a.fq", "fastq2": "b.fq"}) == ["a.fq", "b.fq"]


@pytest.mark.parametrize("supplied", [["a.fq"], "a.fq", None, 7])
def test_an_input_files_value_that_is_not_a_mapping_names_nothing(supplied, caplog) -> None:
    """``vntyper report`` consumes a *supplied* ``pipeline_summary.json`` (#207).

    A malformed one must cost the reader the input-file line, not the report: the
    module this feeds already holds that an empty section is a usable diagnostic
    and a traceback is not.
    """
    with caplog.at_level(logging.ERROR, logger="vntyper.scripts.report_identity"):
        assert input_file_names(supplied) == []
        assert format_input_files(supplied) == NOT_RECORDED

    assert "not a mapping" in caplog.text


# ---------------------------------------------------------------------------
# Provenance rendering
# ---------------------------------------------------------------------------


def test_a_region_is_printed_with_thousands_separators() -> None:
    assert format_region("chr1:155184000-155194000") == "chr1:155,184,000-155,194,000"


def test_an_absent_region_is_not_recorded() -> None:
    assert format_region(None) == NOT_RECORDED
    assert format_region("  ") == NOT_RECORDED


def test_a_region_that_is_not_a_span_is_printed_as_recorded() -> None:
    """Whatever the run wrote is what the report shows; it is not reformatted
    into a shape it does not have."""
    assert format_region("MUC1_custom_region") == "MUC1_custom_region"


def test_a_recorded_value_is_shown_and_an_absent_one_is_not_recorded() -> None:
    assert recorded_or_not("hg38_ensembl") == "hg38_ensembl"
    assert recorded_or_not(1) == "1"
    assert recorded_or_not(None) == NOT_RECORDED
    assert recorded_or_not("") == NOT_RECORDED
    assert recorded_or_not("   ") == NOT_RECORDED


def test_a_run_timestamp_is_labelled_utc() -> None:
    """``start_summary`` writes a naive ISO timestamp that *is* UTC; saying so
    is the difference between a comparable time and an ambiguous one."""
    assert format_run_timestamp("2026-08-13T09:15:33.123456") == "2026-08-13 09:15:33 UTC"


def test_an_offset_run_timestamp_is_converted_to_utc() -> None:
    assert format_run_timestamp("2026-08-13T11:15:33+02:00") == "2026-08-13 09:15:33 UTC"


def test_an_absent_run_timestamp_is_not_recorded() -> None:
    assert format_run_timestamp(None) == NOT_RECORDED
    assert format_run_timestamp("") == NOT_RECORDED


def test_an_unparseable_run_timestamp_is_printed_as_recorded() -> None:
    assert format_run_timestamp("last Tuesday") == "last Tuesday"
