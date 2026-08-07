"""Characterisation for the Kestrel filtering seam.

Scope
-----
``filter_final_dataframe`` (the gate that decides which variants reach the
report), ``select_single_best_variant`` (which picks the one row that becomes
the genotype) and ``construct_kestrel_command`` (the exact command line the
vendored, version-pinned JAR is invoked with).

Everything here is **characterisation** -- it records what the code does today so
a change becomes visible -- *except* the gate-absence tests marked "#185
specification", each of which quotes the decision it encodes in its own
docstring. Where recorded behaviour merely looks dangerous the docstring says so
and stops there; changing which variants survive needs a golden-cohort diff, not
a unit test.

The builders are used only to supply the columns a stage requires. Every value
under test is written out here -- ``tests.builders.kestrel_stage_frame``
independently reimplements frameshift validity and hardcodes ``Confidence``,
the gate booleans and ``Flag``, so asserting production output against a
builder-computed value would be circular.
"""

import logging
from pathlib import Path
from unittest import mock

import pandas as pd
import pytest

from tests.builders import kestrel_config, kestrel_stage_frame
from tests.support.pipeline_harness import run_pipeline_under_harness
from vntyper.scripts import kestrel_genotyping
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
    "flag_filter_pass",
)


def _passing_frame(rows: int = 1) -> pd.DataFrame:
    """Build a frame in which every gate is explicitly True.

    The builder supplies the surrounding columns; the gate values are set here so
    nothing in this file asserts a production result against a builder-computed
    one.

    Args:
        rows: How many rows to build.

    Returns:
        pd.DataFrame: The frame, with all six gate columns True.
    """
    frame = kestrel_stage_frame("flagged", rows=rows).copy()
    for column in GATE_COLUMNS:
        frame[column] = True
    return frame


def test_the_gate_columns_are_exactly_the_six_this_file_pins():
    """Guard against a seventh gate being added without a test.

    ``filter_final_dataframe`` builds its list inline, so this reads the source
    rather than importing a constant that does not exist. The hit-count assertion
    matters: a scan that silently matches nothing would make every claim below
    vacuous.

    The count was five until #174 added ``flag_filter_pass`` -- the config-driven
    artifact gate. Changing this number is the deliberate act that adding a gate
    requires; it is not maintenance.
    """
    from pathlib import Path

    source = Path("vntyper/scripts/kestrel_genotyping.py").read_text(encoding="utf-8")
    body = source.split("filter_cols = [", 1)
    assert len(body) == 2, "could not find the filter_cols list in filter_final_dataframe"
    declared = body[1].split("]", 1)[0]
    found = tuple(line.strip().strip(',"') for line in declared.strip().splitlines())
    assert len(found) == 6, f"expected 6 gate columns in the source, found {len(found)}: {found}"
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
def test_a_gate_that_is_absent_from_a_non_empty_frame_raises(gate, tmp_path, caplog):
    """State 3 of 3 per gate: absent -- and the gate **fails closed** (#185 specification).

    @hassansaei on #185: "Prefer fail loud: a missing required gate column should
    raise (abort the run), not be skipped. [...] The current 2.x behaviour [...]
    can silently permit variants that should have been filtered. That is not
    acceptable for this pipeline. So: do not keep fail-open."

    The six gates are the entire safety net between a raw Kestrel call and a
    reported genotype. Four encode a *pathogenicity* judgement (frameshift
    validity, depth confidence, ALT filtering, motif filtering) and the sixth,
    ``flag_filter_pass`` (#174), encodes an *artifact* judgement -- the row is not
    a candidate variant at all. The stages mark rather than filter (AGENTS.md trap
    4), so a stage that stops emitting its column does not remove a check -- it
    converts that check into a permit. Every gate is parametrised here so none is
    left unguarded.

    Args:
        gate: The boolean column removed from the frame.
        tmp_path: Destination for ``kestrel_pre_result.tsv``.
        caplog: Captures the error message.
    """
    frame = _passing_frame().drop(columns=[gate])
    assert gate not in frame.columns

    with (
        caplog.at_level(logging.ERROR, logger="vntyper.scripts.kestrel_genotyping"),
        pytest.raises(ValueError, match=gate) as raised,
    ):
        filter_final_dataframe(frame, str(tmp_path))

    assert "issue #185" in str(raised.value)
    # Repo convention: logger.error(msg) immediately before raise ValueError(msg).
    errors = [
        record.getMessage()
        for record in caplog.records
        if record.levelno >= logging.ERROR and record.name == "vntyper.scripts.kestrel_genotyping"
    ]
    assert str(raised.value) in errors, f"the raise must be preceded by the same message at ERROR; got {errors}"


def test_every_gate_missing_at_once_raises_on_the_first_one(tmp_path):
    """The compounded case: an upstream stage that returns its input unchanged.

    #185 specification. Stated separately from the per-gate cases because this is
    the shape the failure actually takes in production -- a stage that returns its
    input untouched drops several gates at once, not one. The report names the
    first gate in source order, so the message points at the earliest stage that
    could be responsible.
    """
    frame = _passing_frame(rows=3).drop(columns=list(GATE_COLUMNS))

    with pytest.raises(ValueError, match=GATE_COLUMNS[0]):
        filter_final_dataframe(frame, str(tmp_path))


def test_an_empty_frame_is_the_documented_empty_result_path(tmp_path):
    """#185 specification: the carve-out for empty frames is explicit, not silent.

    @hassansaei: "If a softer path is ever needed for empty early-return frames,
    that should be an explicit, documented empty-result path -- not silent
    omission of safety gates on a non-empty candidate table."

    An empty frame carries no candidate variant, so there is nothing for a missing
    gate to permit. Structurally this is unreachable in production -- every stage
    in ``process_kmer_results`` ahead of the single call site ends in
    ``if df.empty: return df`` -- but the contract is stated here rather than left
    to that structure.
    """
    out = filter_final_dataframe(pd.DataFrame(), str(tmp_path))

    assert out.empty


def test_the_pre_result_tsv_is_still_written_before_the_raise(tmp_path):
    """``kestrel_pre_result.tsv`` is the debuggability artefact (AGENTS.md trap 4).

    #185 specification. It must be on disk before the raise, or the aborted run
    loses the evidence of what the frame held and which gates it still carried --
    which is exactly what a person diagnosing the abort needs.
    """
    frame = _passing_frame(rows=2).drop(columns=["alt_filter_pass"])

    with pytest.raises(ValueError, match="alt_filter_pass"):
        filter_final_dataframe(frame, str(tmp_path))

    pre_result = pd.read_csv(tmp_path / "kestrel_pre_result.tsv", sep="\t")
    assert len(pre_result) == 2, "the pre-filter file must carry every input row"
    assert "alt_filter_pass" not in pre_result.columns, "and record that the gate was genuinely absent"
    assert sorted(pre_result["POS"]) == sorted(frame["POS"])


def _merged_motifs() -> pd.DataFrame:
    """Return the motif annotation table ``motif_correction_and_annotation`` merges on.

    Returns:
        pd.DataFrame: Both halves of the builder's ``X-5`` motif pair.
    """
    return pd.DataFrame({"Motif": ["X", "5"], "Motif_sequence": ["SEQ1", "SEQ2"]})


def _run_kmer_results_with_a_stage_that_loses(gate: str, output_dir: str) -> pd.DataFrame:
    """Drive the real ``process_kmer_results`` with one stage no longer emitting ``gate``.

    Nothing is stubbed except that single omission, so the frame reaching
    ``filter_final_dataframe`` is one the real stages built.

    Args:
        gate: The column ``motif_correction_and_annotation`` is made to drop.
        output_dir: Where the Kestrel artefacts are written.

    Returns:
        pd.DataFrame: Whatever ``process_kmer_results`` returns, if it returns.
    """
    real = kestrel_genotyping.motif_correction_and_annotation

    def dropping(df, merged_motifs, config):
        return real(df, merged_motifs, config).drop(columns=[gate], errors="ignore")

    with mock.patch.object(kestrel_genotyping, "motif_correction_and_annotation", dropping):
        return kestrel_genotyping.process_kmer_results(
            kestrel_stage_frame("raw", rows=2), _merged_motifs(), output_dir, kestrel_config()
        )


def test_the_only_caller_does_not_swallow_the_raise(tmp_path):
    """#185 specification: ``process_kmer_results`` lets the abort through.

    ``filter_final_dataframe`` has exactly one caller. A raise that the caller
    turned into an empty result would be *worse* than the fail-open it replaces:
    the run would report nothing rather than reporting unfiltered variants, and
    would still look successful. This drives the real stage chain from a raw frame
    with one stage made to stop emitting its column -- the production failure the
    decision describes -- and requires the ValueError to reach the caller's caller.
    """
    with pytest.raises(ValueError, match="motif_filter_pass"):
        _run_kmer_results_with_a_stage_that_loses("motif_filter_pass", str(tmp_path))


def test_a_missing_gate_column_aborts_the_run_rather_than_reporting_nothing(tmp_path, caplog):
    """#185 specification: the abort survives the pipeline's stage boundary.

    ``run_pipeline`` wraps its whole body in ``except Exception`` (``pipeline.py``'s
    outer handler). If that turned the raise into a quiet no-result the pipeline
    would still exit 0 with an empty genotype, which is the failure mode the
    decision exists to prevent. The real ``run_pipeline`` runs here with the
    Kestrel stage replaced by the real Kestrel postprocessing over a frame that
    lost a gate; the run must end at exit code 1 with the reason in the log.
    """
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.kestrel_genotyping")

    def kestrel_stage(**kwargs):
        _run_kmer_results_with_a_stage_that_loses("motif_filter_pass", str(Path(kwargs["output_dir"])))

    harness = run_pipeline_under_harness(
        tmp_path / "out",
        expect_failure=True,
        stage_side_effects={"run_kestrel": kestrel_stage},
    )

    assert isinstance(harness.error, SystemExit), f"the run did not exit; it returned {harness.error!r}"
    assert harness.error.code == 1, "a missing gate must abort the run, not finish it with no result"
    reported = [
        record.getMessage()
        for record in caplog.records
        if record.levelno >= logging.ERROR and record.name == "vntyper.scripts.kestrel_genotyping"
    ]
    assert any("issue #185" in message for message in reported), f"the reason never reached the log; got {reported}"
    assert any("motif_filter_pass" in message for message in reported)


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


class TestTheArtifactGateEndToEnd:
    """#174, driven through the real ``process_kmer_results`` from a raw frame.

    Nothing here is stubbed. ``kestrel_config()`` is the shipped configuration, so
    these tests fail if ``artifact_flags`` is removed from it as surely as if the
    call site were removed from the code -- which is the point: the decision is
    configuration, and both halves have to hold for the gate to work.
    """

    def _artifact_frame(self) -> pd.DataFrame:
        """Build the raw frame for a sample whose only call is the 4 bp artifact.

        ``REF=C``/``ALT=CGGCA`` is the flagging rule's condition. Every other gate
        passes on this row -- it is a valid +4 frameshift at ample depth -- so the
        artifact gate is the only thing that can stop it.

        Returns:
            pd.DataFrame: A one-row ``raw``-stage frame.
        """
        return kestrel_stage_frame("raw", rows=1, ref="C", alt="CGGCA")

    def test_an_artifact_only_sample_is_not_reported_as_a_call(self, tmp_path):
        """#174's regression test, end to end through the real postprocessing.

        Before this change the sample yielded one ``High_Precision*`` row carrying
        ``False_Positive_4bp_Insertion``, which ``report_config.json`` maps to
        ``High_Precision_flagged`` and ``is_finding`` reads as positive -- a known
        technical artifact presented as a positive MUC1 call.
        """
        out = kestrel_genotyping.process_kmer_results(
            self._artifact_frame(), _merged_motifs(), str(tmp_path), kestrel_config()
        )

        assert out.empty, "an artifact-only sample has no call"

    def test_the_artifact_row_survives_in_the_pre_result_with_its_gate_false(self, tmp_path):
        """Evidence is never destroyed: stages mark, the final filter drops.

        This is what answers #131's symptom -- a flagged variant vanishing without
        trace. The row and the reason it was dropped are both on disk.
        """
        kestrel_genotyping.process_kmer_results(
            self._artifact_frame(), _merged_motifs(), str(tmp_path), kestrel_config()
        )

        pre = pd.read_csv(tmp_path / "kestrel_pre_result.tsv", sep="\t")
        assert pre["flag_filter_pass"].tolist() == [False]
        assert "False_Positive_4bp_Insertion" in pre["Flag"].iloc[0]
        assert pre["Confidence"].iloc[0] == "High_Precision*", "every other stage still judged it a strong call"

    def test_every_other_gate_still_passed_so_only_the_artifact_gate_dropped_it(self, tmp_path):
        """Locates the exclusion at the new gate rather than merely observing it.

        Without this, a change that broke, say, motif filtering for ``CGGCA`` would
        make the test above pass for entirely the wrong reason.
        """
        kestrel_genotyping.process_kmer_results(
            self._artifact_frame(), _merged_motifs(), str(tmp_path), kestrel_config()
        )

        pre = pd.read_csv(tmp_path / "kestrel_pre_result.tsv", sep="\t")
        for gate in GATE_COLUMNS:
            expected = gate != "flag_filter_pass"
            assert bool(pre[gate].iloc[0]) is expected, f"{gate} was {pre[gate].iloc[0]}, expected {expected}"

    def test_an_advisory_flag_still_reaches_the_result(self, tmp_path):
        """The other half of #174, and the one that would be silently lost.

        ``Low_Depth_Conserved_Motifs`` is advisory: it marks a call that warrants
        scrutiny, not one that is technically impossible. Suppressing those rows
        would delete real low-depth calls -- the exact failure mode this milestone
        exists to prevent. This passed before #174 too, and must keep passing.
        """
        frame = kestrel_stage_frame("raw", rows=1, motif="2", motifs="2-5", motif_sequence="SEQ3")
        merged = pd.DataFrame({"Motif": ["2", "5"], "Motif_sequence": ["SEQ3", "SEQ2"]})

        out = kestrel_genotyping.process_kmer_results(frame, merged, str(tmp_path), kestrel_config())

        assert len(out) == 1
        assert out["Flag"].iloc[0] == "Low_Depth_Conserved_Motifs"
        assert bool(out["flag_filter_pass"].iloc[0]) is True

    def test_emptying_artifact_flags_in_config_restores_the_previous_behaviour(self, tmp_path):
        """The reversibility property, exercised rather than argued.

        #174 rests on the decision living in ``kestrel_config.json``. If narrowing
        or withdrawing it needed a code change, the artifact rule could not be
        revised by whoever owns the domain judgement.
        """
        config = kestrel_config(artifact_flags=[])

        out = kestrel_genotyping.process_kmer_results(self._artifact_frame(), _merged_motifs(), str(tmp_path), config)

        assert len(out) == 1, "with no declared artifacts the row is reported exactly as it was before #174"
        assert out["Flag"].iloc[0] == "False_Positive_4bp_Insertion"

    def test_the_gate_is_present_even_when_no_flagging_rule_ran(self, tmp_path):
        """The reason the call sits *outside* the conditional ``add_flags`` block.

        With no flagging rules and duplicate flagging off, ``add_flags`` is skipped
        entirely and the frame carries no ``Flag`` column at all. A gate written
        inside that ``if`` would then be absent, and ``filter_final_dataframe``
        raises on a missing required gate (#185) -- turning the safety gate into an
        abort of the whole run.
        """
        config = kestrel_config(flagging_rules={}, artifact_flags=["False_Positive_4bp_Insertion"])
        config["duplicate_flagging"]["enabled"] = False

        out = kestrel_genotyping.process_kmer_results(self._artifact_frame(), _merged_motifs(), str(tmp_path), config)

        pre = pd.read_csv(tmp_path / "kestrel_pre_result.tsv", sep="\t")
        assert "Flag" not in pre.columns, "the premise of this test is that add_flags never ran"
        assert pre["flag_filter_pass"].tolist() == [True]
        assert len(out) == 1


def test_a_non_gate_boolean_column_does_not_filter_anything(tmp_path):
    """Only the six listed columns gate. A seventh boolean is inert.

    Without this, "the gates are these six" would be an assertion about the six
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

    # -- shell quoting -----------------------------------------------------------------
    #
    # ``run_command`` executes this string under ``shell=True`` (trap 9), which makes
    # quoting the caller's job. The web service accepts filenames that reach this line,
    # so a space alone splits one argument into two and a ``;`` is arbitrary command
    # execution. ``java_path`` stays unquoted: config.json holds a command *prefix*
    # there, which quoting would collapse into one token.

    HOSTILE = ["/data/my samples/R1.fastq.gz", "/data/o'brien/R1.fastq.gz", "/data/x; id/R1.fastq.gz"]

    @pytest.mark.parametrize("hostile", HOSTILE)
    @pytest.mark.parametrize("argument,flag", [("fastq_1", None), ("reference_vntr", "-r"), ("vcf_out", "-o")])
    def test_a_hostile_path_survives_as_one_shell_word(self, argument, flag, hostile):
        """Each operand, rendered and then re-split the way bash would split it."""
        import shlex

        words = shlex.split(self._command(**{argument: hostile}))

        if flag is None:
            assert hostile in words, f"{hostile!r} did not survive as one word"
        else:
            assert words[words.index(flag) + 1] == hostile

    def test_a_sample_name_with_a_space_stays_attached_to_the_s_flag(self):
        """``-s`` takes its value with no space, so quoting must not introduce one."""
        import shlex

        words = shlex.split(self._command(sample_name="patient 7"))

        assert "-spatient 7" in words
        assert "-s" not in words, "the sample name must not become a separate argument"

    def test_a_hostile_output_directory_survives_in_both_places(self):
        import shlex

        words = shlex.split(self._command(output_dir="/out/my run"))

        assert words[words.index("-p") + 1] == "/out/my run/output.sam"
        assert words[words.index("--temploc") + 1] == "/out/my run"

    def test_the_java_invocation_stays_several_words(self):
        """Quoting it would make bash look for one binary named `mamba run ... java`."""
        assert self._command(java_path="mamba run -n vntyper java").startswith("mamba run -n vntyper java -Xmx12g ")
