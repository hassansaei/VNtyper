"""Characterisation for the Kestrel filtering seam.

Scope
-----
``filter_final_dataframe`` (the gate that decides which variants reach the
report), ``select_single_best_variant`` (which picks the one row that becomes
the genotype) and ``construct_kestrel_command`` (the exact command line the
vendored, version-pinned JAR is invoked with).

Everything here is **characterisation**. It records what the code does today so
a change becomes visible. Where the recorded behaviour looks dangerous the
docstring says so and stops there; changing which variants survive needs a
golden-cohort diff, not a unit test.

The builders are used only to supply the columns a stage requires. Every value
under test is written out here -- ``tests.builders.kestrel_stage_frame``
independently reimplements frameshift validity and hardcodes ``Confidence``,
the gate booleans and ``Flag``, so asserting production output against a
builder-computed value would be circular.
"""

import logging

import pandas as pd
import pytest

from tests.builders import kestrel_config, kestrel_stage_frame
from vntyper.scripts.kestrel_genotyping import (
    construct_kestrel_command,
    filter_final_dataframe,
    select_single_best_variant,
)

pytestmark = pytest.mark.unit

# The boolean columns ``filter_final_dataframe`` ANDs together, in source order.
# Read off the function, not off a builder.
GATE_COLUMNS = (
    "is_frameshift",
    "is_valid_frameshift",
    "depth_confidence_pass",
    "alt_filter_pass",
    "motif_filter_pass",
)


def _passing_frame(rows: int = 1) -> pd.DataFrame:
    """Build a frame in which every gate is explicitly True.

    The builder supplies the surrounding columns; the gate values are set here so
    nothing in this file asserts a production result against a builder-computed
    one.

    Args:
        rows: How many rows to build.

    Returns:
        pd.DataFrame: The frame, with all five gate columns True.
    """
    frame = kestrel_stage_frame("flagged", rows=rows).copy()
    for column in GATE_COLUMNS:
        frame[column] = True
    return frame


def test_the_gate_columns_are_exactly_the_five_this_file_pins():
    """Guard against a sixth gate being added without a test.

    ``filter_final_dataframe`` builds its list inline, so this reads the source
    rather than importing a constant that does not exist. The hit-count assertion
    matters: a scan that silently matches nothing would make every claim below
    vacuous.
    """
    from pathlib import Path

    source = Path("vntyper/scripts/kestrel_genotyping.py").read_text(encoding="utf-8")
    body = source.split("filter_cols = [", 1)
    assert len(body) == 2, "could not find the filter_cols list in filter_final_dataframe"
    declared = body[1].split("]", 1)[0]
    found = tuple(line.strip().strip(',"') for line in declared.strip().splitlines())
    assert len(found) == 5, f"expected 5 gate columns in the source, found {len(found)}: {found}"
    assert found == GATE_COLUMNS


@pytest.mark.parametrize("gate", GATE_COLUMNS)
def test_a_gate_that_is_present_and_true_lets_the_row_through(gate, tmp_path):
    """State 1 of 3 per gate: present and True.

    Args:
        gate: The boolean column under test.
        tmp_path: Destination for ``kestrel_pre_result.tsv``.
    """
    frame = _passing_frame()
    frame[gate] = True

    out = filter_final_dataframe(frame, str(tmp_path))

    assert len(out) == 1


@pytest.mark.parametrize("gate", GATE_COLUMNS)
def test_a_gate_that_is_present_and_false_drops_the_row(gate, tmp_path):
    """State 2 of 3 per gate: present and False.

    Args:
        gate: The boolean column under test.
        tmp_path: Destination for ``kestrel_pre_result.tsv``.
    """
    frame = _passing_frame()
    frame[gate] = False

    out = filter_final_dataframe(frame, str(tmp_path))

    assert out.empty, f"a False {gate!r} must drop the row"


@pytest.mark.parametrize("gate", GATE_COLUMNS)
def test_a_gate_that_is_absent_is_ignored_and_the_row_survives(gate, tmp_path, caplog):
    """State 3 of 3 per gate: absent -- and the gate **fails open**.

    ``filter_final_dataframe`` skips a gate column that is not in the frame. A
    stage that crashes, returns early, or is reordered out of the chain therefore
    does not block anything: its gate silently stops existing and every row sails
    past it.

    That is dangerous. The five gates are the entire safety net between a raw
    Kestrel call and a reported genotype, and four of the five encode a
    *pathogenicity* judgement (frameshift validity, depth confidence, ALT
    filtering, motif filtering). Losing one to an upstream failure turns a missing
    check into a permissive one, which is the wrong direction for a diagnostic.

    This test pins the fail-open behaviour rather than fixing it: making the gate
    fail closed changes which variants reach the report and needs a golden-cohort
    diff first.

    Args:
        gate: The boolean column removed from the frame.
        tmp_path: Destination for ``kestrel_pre_result.tsv``.
        caplog: Captures the skip message.
    """
    frame = _passing_frame().drop(columns=[gate])
    assert gate not in frame.columns

    with caplog.at_level(logging.INFO, logger="vntyper.scripts.kestrel_genotyping"):
        out = filter_final_dataframe(frame, str(tmp_path))

    assert len(out) == 1, f"a missing {gate!r} is ignored rather than treated as a failure"
    skipped = [
        record
        for record in caplog.records
        if record.levelno >= logging.INFO and record.getMessage() == f"Filter column '{gate}' not found; skipping."
    ]
    assert skipped, (
        f"the skip is announced at INFO, not merely at DEBUG; got {[r.getMessage() for r in caplog.records]}"
    )


def test_every_gate_missing_at_once_still_lets_everything_through(tmp_path):
    """The fail-open behaviour compounds: with no gates at all, nothing is filtered.

    Stated separately from the per-gate cases because this is the shape the failure
    actually takes in production -- an upstream stage that returns its input
    unchanged drops several gates at once, not one.
    """
    frame = _passing_frame(rows=3).drop(columns=list(GATE_COLUMNS))

    out = filter_final_dataframe(frame, str(tmp_path))

    assert len(out) == 1, "three ungated rows collapse to one only because of best-variant selection"


def test_the_pre_result_tsv_keeps_every_input_row_including_the_rejected_ones(tmp_path):
    """``kestrel_pre_result.tsv`` is the debuggability contract: it is pre-filter.

    Stages mark, they do not filter; ``filter_final_dataframe`` ANDs the marks at
    the end and writes the *unfiltered* frame to disk first. Without that file
    there is no way to see why a variant was dropped, which is the difference
    between "no variant found" and "a variant was found and discarded by gate X".
    """
    frame = _passing_frame(rows=4)
    frame.loc[1, "motif_filter_pass"] = False
    frame.loc[2, "is_valid_frameshift"] = False
    frame.loc[3, "depth_confidence_pass"] = False

    out = filter_final_dataframe(frame, str(tmp_path))
    assert len(out) == 1, "only the row with every gate True survives"

    pre_result = pd.read_csv(tmp_path / "kestrel_pre_result.tsv", sep="\t")
    assert len(pre_result) == len(frame) == 4, "the pre-filter file must carry every input row, not just the survivors"
    assert set(GATE_COLUMNS).issubset(pre_result.columns), "and every gate column, so the reason is recoverable"
    assert sorted(pre_result["POS"]) == sorted(frame["POS"]), "rows must be identifiable, not merely counted"
    assert int(pre_result["motif_filter_pass"].sum()) == 3
    assert int(pre_result["is_valid_frameshift"].sum()) == 3
    assert int(pre_result["depth_confidence_pass"].sum()) == 3


def test_a_non_gate_boolean_column_does_not_filter_anything(tmp_path):
    """Only the five listed columns gate. A sixth boolean is inert.

    Without this, "the gates are these five" would be an assertion about the five
    and silence about everything else.
    """
    frame = _passing_frame()
    frame["some_other_flag"] = False

    out = filter_final_dataframe(frame, str(tmp_path))

    assert len(out) == 1


def test_select_single_best_variant_does_not_mutate_the_callers_frame():
    """The caller's frame must not gain ``_priority`` or ``_is_unflagged``.

    ``select_single_best_variant`` adds those two sort keys, then drops them from
    the row it returns. If the copy at the top of the function were removed, the
    return value would still look right and only the caller's frame would carry
    the leak -- straight into ``kestrel_pre_result.tsv`` and any later consumer of
    the same object.
    """
    frame = pd.DataFrame(
        {
            "Confidence": ["High_Precision", "Low_Precision"],
            "haplo_count": [100, 50],
            "Depth_Score": [0.01, 0.005],
            "POS": [67, 54],
            "Flag": ["Not flagged", "Low_Depth_Conserved_Motifs"],
        }
    )
    columns_before = list(frame.columns)
    rows_before = len(frame)

    result = select_single_best_variant(frame)

    assert list(frame.columns) == columns_before, "the caller's frame gained a column"
    assert "_priority" not in frame.columns
    assert "_is_unflagged" not in frame.columns
    assert len(frame) == rows_before, "the caller's frame lost or gained rows"
    assert list(result.columns) == columns_before


def test_filter_final_dataframe_does_not_leak_sort_keys_into_the_result(tmp_path):
    """The same leak seen through the caller that actually invokes the selection."""
    frame = _passing_frame(rows=2)

    out = filter_final_dataframe(frame, str(tmp_path))

    assert len(out) == 1
    assert "_priority" not in out.columns
    assert "_is_unflagged" not in out.columns

    pre_result = pd.read_csv(tmp_path / "kestrel_pre_result.tsv", sep="\t")
    assert "_priority" not in pre_result.columns
    assert "_is_unflagged" not in pre_result.columns


class TestConstructKestrelCommand:
    """The exact command line the pinned Kestrel 1.0.1 JAR is invoked with.

    Kestrel is pinned because the scoring thresholds are calibrated to it, so the
    invocation is part of that calibration. The ``-s`` spacing in particular has
    flip-flopped twice in one day historically: Kestrel takes the sample name
    attached to the flag with **no space**, and inserting one silently changes what
    the sample is called in the VCF.
    """

    ARGS = {
        "kestrel_path": "/opt/kestrel.jar",
        "reference_vntr": "/ref/muc1.fa",
        "output_dir": "/out",
        "fastq_1": "/in/R1.fastq.gz",
        "fastq_2": "/in/R2.fastq.gz",
        "vcf_out": "/out/output.vcf",
        "java_path": "/usr/bin/java",
        "java_memory": "12g",
        "log_level": "info",
        "sample_name": "SAMPLE1",
    }

    def _command(self, **overrides) -> str:
        """Render the command with the shipped Kestrel settings.

        Args:
            **overrides: Replace any argument.

        Returns:
            str: The rendered command line.
        """
        settings = kestrel_config()["kestrel_settings"]
        args = {
            **self.ARGS,
            "kmer_size": settings["kmer_sizes"][0],
            "max_align_states": settings["max_align_states"],
            "max_hap_states": settings["max_hap_states"],
            "additional_settings": settings["additional_settings"],
        }
        args.update(overrides)
        return construct_kestrel_command(**args)

    def test_the_shipped_settings_are_the_calibrated_ones(self):
        """Pin k=20 and 40/40 states in the config itself, not only in the rendering."""
        settings = kestrel_config()["kestrel_settings"]
        assert settings["kmer_sizes"] == [20]
        assert settings["max_align_states"] == 40
        assert settings["max_hap_states"] == 40
        assert settings["java_memory"] == "12g"
        assert settings["additional_settings"] == ""

    def test_the_rendered_command_is_exactly_this(self):
        """One whole-string assertion, so any reordering or respacing fails here."""
        assert self._command() == (
            "/usr/bin/java -Xmx12g -jar /opt/kestrel.jar -k 20 "
            "--maxalignstates 40 --maxhapstates 40 "
            "-r /ref/muc1.fa -o /out/output.vcf "
            "-sSAMPLE1 /in/R1.fastq.gz /in/R2.fastq.gz "
            "--hapfmt sam -p /out/output.sam --logstderr --logstdout "
            "--loglevel INFO --temploc /out"
        )

    @pytest.mark.parametrize("fragment", ["-k 20", "--maxalignstates 40", "--maxhapstates 40"])
    def test_the_calibrated_parameters_reach_the_command_line(self, fragment):
        """The three tuned numeric parameters, each named individually.

        Args:
            fragment: The exact substring that must appear.
        """
        assert fragment in self._command()

    def test_the_sample_name_is_attached_to_the_s_flag_with_no_space(self):
        """``-sSAMPLE1``, never ``-s SAMPLE1``."""
        command = self._command()
        assert "-sSAMPLE1 " in command
        assert "-s SAMPLE1" not in command
        assert " -s " not in command

    def test_the_log_level_is_upper_cased(self):
        """``--loglevel`` takes the upper-cased name; ``run_kestrel`` passes a mixed-case one."""
        assert "--loglevel DEBUG" in self._command(log_level="debug")

    def test_additional_settings_are_appended_last(self):
        """An extra-settings string lands at the end, after ``--temploc``."""
        assert self._command(additional_settings="--flank 5").endswith("--temploc /out --flank 5")

    @pytest.mark.parametrize(
        ("fastq_1", "fastq_2"),
        [(None, "/in/R2.fastq.gz"), ("/in/R1.fastq.gz", None), (None, None), ("", "/in/R2.fastq.gz")],
    )
    def test_a_missing_fastq_raises_value_error(self, fastq_1, fastq_2):
        """Either FASTQ missing is a hard error, not a half-built command.

        Args:
            fastq_1: R1 path, or a falsy placeholder.
            fastq_2: R2 path, or a falsy placeholder.
        """
        with pytest.raises(ValueError, match="FASTQ input files are missing or invalid"):
            self._command(fastq_1=fastq_1, fastq_2=fastq_2)
