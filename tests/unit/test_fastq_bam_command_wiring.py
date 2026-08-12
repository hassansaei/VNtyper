"""
Pins the commands ``fastq_bam_processing`` and ``alignment_processing`` emit.

``tests/unit/test_command_builders.py`` pins the builders. This module pins the
**wiring**: that each stage calls the right builder with the right arguments, in
the right order, and that no stage has quietly kept an inline f-string. Those are
different failures - a builder can be perfectly quoted while the caller still
hand-rolls its own string next to it.

Every expected string here was captured from the code as it stood *before*
``command_builders`` existed, so this file doubles as the proof that the
extraction is behaviour-preserving. The only differences from that capture are the
five deliberate fixes:

* ``set -o pipefail; `` on the three multi-stage pipes,
* the CRAM unmapped-read extractor calling the configured samtools instead of a
  bare ``samtools``,
* ``shlex.quote`` around interpolated paths, which is a no-op for every path that
  needs no quoting - hence the byte-identical strings below,
* the CRAM unmapped-read extractor writing through a plain pipe rather than a
  ``tee >(...)`` process substitution the shell does not wait for - see the
  section on that path below for the measurement, and
* dedup-enabled fastp commands using one worker so the shared atomic duplicate
  table retains a reproducible representative.

``run_command`` is mocked throughout; nothing here starts a process.
"""

import shlex
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts import alignment_processing, fastq_bam_processing
from vntyper.scripts.alignment_contract import AlignmentPlan

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit

CONFIG = {
    "tools": {"samtools": "samtools", "bwa": "bwa", "fastp": "fastp"},
    "bam_processing": {
        "compression_level": 6,
        "disable_adapter_trimming": True,
        "deduplication": True,
        "dup_calc_accuracy": 3,
        "length_required": 50,
        "qualified_quality_phred": 20,
    },
}

REGION = "chr1:155158000-155163000"


class _Recorder:
    """A ``run_command`` stand-in that records the command strings it is given."""

    def __init__(self):
        self.commands: list[str] = []

    def __call__(self, command, log_file, critical=False, cwd=None):
        self.commands.append(command)
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        Path(log_file).write_text("")
        return True


def _run_bam_to_fastq(tmp_path, **overrides):
    """Drive ``process_bam_to_fastq`` with every external effect mocked out."""
    recorder = _Recorder()
    input_path = overrides.pop("input_path", "/data/sample.bam")
    kwargs = {
        "output": str(tmp_path),
        "output_name": "output",
        "threads": 4,
        "config": CONFIG,
        "fast_mode": True,
    }
    file_format = overrides.pop("file_format", "bam")
    kwargs.update(overrides)
    kwargs.setdefault(
        "plan",
        AlignmentPlan(
            input_path=str(input_path),
            view_path=str(input_path),
            file_format=file_format,
            index_path=f"{input_path}.{'bai' if file_format == 'bam' else 'crai'}",
            reference_path=None,
            reference_source="test",
            uncovered_contigs=(),
            unmapped_scan="stream" if file_format == "cram" else "indexed",
        ),
    )

    with (
        patch.object(fastq_bam_processing, "run_command", recorder),
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value=REGION),
        patch.object(fastq_bam_processing.os, "replace"),
    ):
        fastq_bam_processing.process_bam_to_fastq(**kwargs)

    return recorder.commands


# ---------------------------------------------------------------------------
# process_fastq
# ---------------------------------------------------------------------------


def test_process_fastq_emits_the_pinned_fastp_command(tmp_path):
    """The caller requests four threads, but dedup-enabled fastp is serialized."""
    recorder = _Recorder()

    with patch.object(fastq_bam_processing, "run_command", recorder):
        fastq_bam_processing.process_fastq("/data/in_R1.fq.gz", "/data/in_R2.fq.gz", 4, str(tmp_path), "output", CONFIG)

    assert recorder.commands == [
        f"fastp --thread 1 --in1 /data/in_R1.fq.gz --in2 /data/in_R2.fq.gz "
        f"--out1 {tmp_path}/output_R1.fastq.gz --out2 {tmp_path}/output_R2.fastq.gz "
        f"--compression 6 "
        f"--qualified_quality_phred 20 "
        f"--dup_calc_accuracy 3 "
        f"--length_required 50 "
        f"--html {tmp_path}/output.html "
        f"--json {tmp_path}/output.json "
        f" --disable_adapter_trimming --dedup"
    ]


@pytest.mark.parametrize(
    ("disable_adapter_trimming", "deduplication", "expected_threads"),
    [
        (True, True, 1),
        (True, False, 4),
        (False, True, 1),
        (False, False, 4),
    ],
)
def test_process_fastq_serializes_only_deduplication(
    tmp_path, disable_adapter_trimming, deduplication, expected_threads
):
    """Adapter trimming does not influence the fastp worker count."""
    recorder = _Recorder()
    config = {
        **CONFIG,
        "bam_processing": {
            **CONFIG["bam_processing"],
            "disable_adapter_trimming": disable_adapter_trimming,
            "deduplication": deduplication,
        },
    }

    with patch.object(fastq_bam_processing, "run_command", recorder):
        fastq_bam_processing.process_fastq("/data/in_R1.fq.gz", "/data/in_R2.fq.gz", 4, str(tmp_path), "output", config)

    assert recorder.commands[0].startswith(f"fastp --thread {expected_threads} ")


def test_process_fastq_raises_when_the_command_fails(tmp_path):
    """A failed fastp run must abort the stage rather than return silently."""
    with (
        patch.object(fastq_bam_processing, "run_command", return_value=False),
        pytest.raises(RuntimeError, match="FASTQ quality control failed"),
    ):
        fastq_bam_processing.process_fastq("/data/in_R1.fq.gz", "/data/in_R2.fq.gz", 4, str(tmp_path), "output", CONFIG)


def test_process_fastq_emits_the_single_input_fastp_form(tmp_path):
    recorder = _Recorder()

    with patch.object(fastq_bam_processing, "run_command", recorder):
        fastq_bam_processing.process_fastq("/data/single.fq.gz", None, 4, str(tmp_path), "output", CONFIG)

    assert "--in1 /data/single.fq.gz" in recorder.commands[0]
    assert "--out1" in recorder.commands[0]
    assert "--in2" not in recorder.commands[0]
    assert "--out2" not in recorder.commands[0]


# ---------------------------------------------------------------------------
# process_bam_to_fastq - one assertion per branch
# ---------------------------------------------------------------------------


def test_the_bam_fast_mode_path_slices_then_converts(tmp_path):
    """Fast mode skips unmapped-read extraction entirely: two commands, in order."""
    commands = _run_bam_to_fastq(tmp_path)

    assert commands == [
        f"samtools view -P -b -@ 4 -X /data/sample.bam /data/sample.bam.bai {REGION} "
        f"-o {tmp_path}/output_sliced.bam && "
        f"samtools index -@ 4 {tmp_path}/output_sliced.bam",
        f"set -o pipefail; samtools sort -n -@ 4 {tmp_path}/output_sliced.bam | "
        f"samtools fastq -@ 4 - -1 {tmp_path}/output_R1.fastq.gz "
        f"-2 {tmp_path}/output_R2.fastq.gz -0 {tmp_path}/output_other.fastq.gz "
        f"-s {tmp_path}/output_single.fastq.gz",
    ]


def test_the_bam_normal_path_indexes_extracts_merges_and_reindexes(tmp_path):
    """
    The full BAM path, in order.

    Preflight has already proved and supplied the input index. The temporary slice
    is not indexed because the merge immediately replaces it; the merged BAM is
    indexed once before conversion.
    """
    commands = _run_bam_to_fastq(tmp_path, fast_mode=False)

    assert commands == [
        f"samtools view -P -b -F 4 -@ 4 -X /data/sample.bam /data/sample.bam.bai {REGION} "
        f"-o {tmp_path}/output_sliced.bam",
        f"samtools view -b -f 4 -@ 4 -X /data/sample.bam /data/sample.bam.bai '*' -o {tmp_path}/output_unmapped.bam",
        f"samtools merge -f -@ 4 {tmp_path}/output_sliced_unmapped.bam "
        f"{tmp_path}/output_sliced.bam {tmp_path}/output_unmapped.bam",
        f"samtools index -@ 4 {tmp_path}/output_sliced.bam",
        f"set -o pipefail; samtools sort -n -@ 4 {tmp_path}/output_sliced.bam | "
        f"samtools fastq -@ 4 - -1 {tmp_path}/output_R1.fastq.gz "
        f"-2 {tmp_path}/output_R2.fastq.gz -0 {tmp_path}/output_other.fastq.gz "
        f"-s {tmp_path}/output_single.fastq.gz",
    ]


@pytest.mark.parametrize(
    ("delete_intermediates", "keep_intermediates", "expected_exists"),
    [(False, False, True), (False, True, True), (True, False, False), (True, True, False)],
)
def test_delete_intermediates_overrides_keep(
    tmp_path: Path,
    delete_intermediates: bool,
    keep_intermediates: bool,
    expected_exists: bool,
) -> None:
    unmapped = tmp_path / "output_unmapped.bam"
    unmapped.write_bytes(b"intermediate")
    _run_bam_to_fastq(
        tmp_path,
        fast_mode=False,
        delete_intermediates=delete_intermediates,
        keep_intermediates=keep_intermediates,
    )
    assert unmapped.exists() is expected_exists


def test_indexed_bam_recovery_uses_the_htslib_literal_star_command(tmp_path):
    """A valid all-unplaced BAM has no BAI chunk offset for the custom reader."""
    commands = _run_bam_to_fastq(tmp_path, fast_mode=False)

    assert [command for command in commands if "-f 4" in command] == [
        f"samtools view -b -f 4 -@ 4 -X /data/sample.bam /data/sample.bam.bai '*' -o {tmp_path}/output_unmapped.bam"
    ]


def _index_commands_for(commands: list[str], indexed_file: str) -> list[list[str]]:
    """Return the tokenised ``samtools index`` calls that index ``indexed_file``.

    Commands are split on ``&&`` first: the slice command ends in an ``index`` of
    the file it just wrote *and* names the input BAM earlier in the same string,
    so a substring test over the whole command matches it too.

    Args:
        commands: The shell command strings the stage emitted.
        indexed_file: The file whose index calls are wanted.

    Returns:
        list[list[str]]: One token list per matching ``samtools index`` call.
    """
    matches = []
    for command in commands:
        for segment in command.split("&&"):
            tokens = shlex.split(segment)
            if "index" in tokens and tokens[-1] == indexed_file:
                matches.append(tokens)
    return matches


def test_the_bam_stage_never_rebuilds_the_preflight_index(tmp_path):
    """
    Preflight owns input index resolution and construction. The conversion stage
    must not repeat that decision or write beside the patient alignment.
    """
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    in_bam = input_dir / "sample.bam"
    in_bam.write_bytes(b"BAM\x01")

    commands = _run_bam_to_fastq(tmp_path, fast_mode=False, input_path=str(in_bam), output=str(output_dir))

    index_commands = _index_commands_for(commands, str(in_bam))
    assert index_commands == []


@pytest.mark.parametrize("index_name", ["sample.bam.bai", "sample.bai"])
def test_an_existing_index_beside_the_input_is_reused_rather_than_rebuilt(tmp_path, index_name):
    """
    #210's other half: both names htslib resolves are honoured.

    ``sample.bai`` is a name the upload endpoint and the web worker both
    deliberately accept. Reconstructing only ``<file>.bai`` meant the pipeline did
    not see it, so it built a second index that nothing else knew about -- next to
    the patient's alignment.
    """
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    in_bam = input_dir / "sample.bam"
    in_bam.write_bytes(b"BAM\x01")
    (input_dir / index_name).write_bytes(b"BAI\x01")

    commands = _run_bam_to_fastq(tmp_path, fast_mode=False, input_path=str(in_bam), output=str(output_dir))

    assert _index_commands_for(commands, str(in_bam)) == [], (
        f"{index_name} was already there; indexing again writes a file nothing reads"
    )


def test_indexed_bam_recovery_uses_the_plan_view_and_exact_index_without_a_reference(tmp_path):
    """The htslib fetch consumes the proven view and its retained exact BAI."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    in_bam = input_dir / "sample.bam"
    in_bam.write_bytes(b"BAM\x01")

    plan = AlignmentPlan(
        input_path=str(in_bam),
        view_path=str(tmp_path / "view.bam"),
        file_format="bam",
        index_path=str(tmp_path / "view.bam.bai"),
        reference_path=None,
        reference_source="not-required",
        uncovered_contigs=(),
        unmapped_scan="indexed",
    )

    recorder = _Recorder()
    with (
        patch.object(fastq_bam_processing, "run_command", recorder),
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value=REGION),
        patch.object(fastq_bam_processing.os, "replace"),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=str(tmp_path),
            output_name="output",
            threads=4,
            config=CONFIG,
            plan=plan,
            fast_mode=False,
        )

    (unmapped_command,) = [command for command in recorder.commands if "-f 4" in command]
    assert unmapped_command == (
        f"samtools view -b -f 4 -@ 4 -X {plan.view_path} {plan.index_path} '*' -o {tmp_path}/output_unmapped.bam"
    )
    assert " -T " not in unmapped_command


@pytest.mark.parametrize("index_path", ["sample.bam.bai", "sample.bai", None])
def test_indexed_bam_recovery_uses_the_plan_index_operand(tmp_path, index_path):
    """Manual plans still pass their declared index explicitly to htslib."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    in_bam = input_dir / "sample.bam"
    in_bam.write_bytes(b"BAM\x01")
    if index_path is not None:
        (input_dir / index_path).write_bytes(b"BAI\x01")

    recorder = _Recorder()
    with (
        patch.object(fastq_bam_processing, "run_command", recorder),
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value=REGION),
        patch.object(fastq_bam_processing.os, "replace"),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=str(output_dir),
            output_name="output",
            threads=4,
            config=CONFIG,
            fast_mode=False,
            plan=AlignmentPlan(
                input_path=str(in_bam),
                view_path=str(in_bam),
                file_format="bam",
                index_path=str(input_dir / index_path) if index_path is not None else str(output_dir / "proven.bai"),
                reference_path=None,
                reference_source="test",
                uncovered_contigs=(),
                unmapped_scan="indexed",
            ),
        )

    (unmapped_command,) = [command for command in recorder.commands if "-f 4" in command]
    asserted_index = input_dir / index_path if index_path is not None else output_dir / "proven.bai"
    assert f"-X {in_bam} {asserted_index}" in unmapped_command
    assert in_bam.as_posix() in unmapped_command


# ---------------------------------------------------------------------------
# The CRAM unmapped-read path
#
# This branch used to be
#     ... | tee >(samtools view -b -f 4 - -o unmapped.bam) > /dev/null
# and bash does not wait for a ``>(...)`` process substitution: the shell returned
# as soon as ``tee`` exited, while the substituted samtools was still flushing
# ``unmapped.bam``. ``process_bam_to_fastq`` runs ``samtools merge`` against that
# file on the very next line.
#
# Measured against a synthetic 600k-read CRAM with samtools 1.20, ten trials each:
# with the substitution the file held 199,797 of 200,000 unmapped reads at the
# instant the shell returned and ``samtools merge`` accepted it with only a
# ``W::bam_hdr_read`` warning; with the plain pipe below it held 200,000 every time.
# The reads at stake are exactly the ones this stage exists to recover for Kestrel,
# so the old shape could under-call a sample with no pipeline error at all.
#
# The tests below pin the pipe. ``tee`` is not merely unnecessary, it is the defect:
# its own stdout went to ``/dev/null``, so the substitution was the only consumer.
# ---------------------------------------------------------------------------


def _cram_unmapped_command(tmp_path, **overrides):
    """
    Drive the CRAM branch and return its single unmapped-read extraction command.

    Args:
        tmp_path (Path): pytest temporary directory.
        **overrides: Passed through to ``_run_bam_to_fastq``.

    Returns:
        str: The one emitted command containing ``-f 4``.
    """
    config = {**CONFIG, "tools": {**CONFIG["tools"], "samtools": "/envs/vntyper/bin/samtools"}}
    kwargs = {
        "fast_mode": False,
        "file_format": "cram",
        "input_path": "/data/sample.cram",
        "config": config,
    }
    kwargs.update(overrides)

    commands = _run_bam_to_fastq(tmp_path, **kwargs)
    filters = [c for c in commands if "-f 4" in c]
    assert len(filters) == 1, f"expected exactly one unmapped-read extraction command, got {filters}"
    return filters[0]


def test_the_cram_unmapped_command_is_pinned_as_a_plain_pipe(tmp_path):
    """
    The whole command, byte for byte.

    With no resolved reference yet, the optional ``-T`` fragment is omitted.
    """
    assert _cram_unmapped_command(tmp_path) == (
        f"set -o pipefail; /envs/vntyper/bin/samtools view -@ 4 -h /data/sample.cram | "
        f"/envs/vntyper/bin/samtools view -b -f 4 -@ 4 - -o {tmp_path}/output_unmapped.bam"
    )


@pytest.mark.parametrize("scan", ["indexed", "stream"])
def test_a_no_ref_cram_plan_emits_no_reference_flag_in_any_conversion_command(tmp_path, scan):
    """A-209-3: the proven no-reference plan must remain byte-free of ``-T``."""
    input_path = "/data/no-ref.cram"
    plan = AlignmentPlan(
        input_path=input_path,
        view_path=input_path,
        file_format="cram",
        index_path=f"{input_path}.crai",
        reference_path=None,
        reference_source="htslib-resolved (header UR: or REF_PATH)",
        uncovered_contigs=(),
        unmapped_scan=scan,
    )

    commands = _run_bam_to_fastq(
        tmp_path,
        file_format="cram",
        input_path=input_path,
        fast_mode=False,
        plan=plan,
    )

    assert commands
    assert all(" -T " not in command for command in commands)


def test_the_cram_unmapped_command_has_no_process_substitution(tmp_path):
    """
    No ``>(...)``, and no ``tee`` for it to hang off.

    This is the assertion that fails if the race is reintroduced. ``tee``'s stdout
    went to ``/dev/null``, so the substitution was the pipeline's only consumer -
    there is nothing for a second output to feed.
    """
    command = _cram_unmapped_command(tmp_path)

    assert ">(" not in command, "a process substitution is not waited for; the merge would race the writer"
    assert "tee" not in command, "tee had exactly one consumer, so it was pure overhead"
    assert "/dev/null" not in command, "nothing is discarded any more; the writer consumes the whole stream"


def test_the_cram_unmapped_writer_is_a_pipeline_stage(tmp_path):
    """
    The writing samtools must be a **stage**, which is what makes bash wait for it.

    Two stages joined by one ``|``: the reader streams the CRAM, the writer filters
    flag 4 into the output BAM. Because the writer is in the pipeline, the shell
    does not return until it has exited, and ``pipefail`` covers its exit status.
    """
    command = _cram_unmapped_command(tmp_path)
    stages = command.removeprefix("set -o pipefail; ").split(" | ")

    assert len(stages) == 2, f"expected a two-stage pipeline, got {len(stages)}: {stages}"
    reader, writer = stages
    assert reader.endswith("-h /data/sample.cram"), f"the reader must end at the CRAM input: {reader}"
    assert writer.startswith("/envs/vntyper/bin/samtools view -b -f 4"), f"the writer is not stage two: {writer}"
    assert writer.endswith(f"- -o {tmp_path}/output_unmapped.bam"), (
        f"the writer must read stdin and write the unmapped BAM: {writer}"
    )


def test_the_cram_unmapped_command_sets_pipefail(tmp_path):
    """
    ``pipefail`` is only meaningful now that the writer is in the pipeline.

    Under the substitution it bought nothing: a writer that consumed its whole input
    and then failed never made ``tee`` see EPIPE, so the command exited 0 and the
    stage carried on with a BAM that was never finished.
    """
    command = _cram_unmapped_command(tmp_path)

    assert "|" in command, "this test is checking the wrong command if there is no pipe in it"
    assert command.startswith("set -o pipefail; "), "without pipefail a failed writer exits 0"


def test_the_cram_unmapped_command_uses_the_configured_samtools_in_both_stages(tmp_path):
    """
    D2 through the real call site.

    The stage passes ``config["tools"]["samtools"]`` to the builder for both stages.
    A bare ``samtools`` in either would resolve against whatever PATH ``mamba run``
    set up, so the unmapped reads would be extracted by a different build - or by
    nothing at all, while the pipeline reported success.
    """
    command = _cram_unmapped_command(tmp_path)

    assert command.count("/envs/vntyper/bin/samtools") == 2, f"both stages must use the configured samtools: {command}"
    assert [token for token in command.split() if token == "samtools"] == [], (
        f"a bare `samtools` token means the mismatch is back: {command}"
    )


def test_a_hostile_cram_path_survives_both_stages_of_the_pipe(tmp_path):
    """
    Every interpolated path stays one operand, on both sides of the ``|``.

    Now that the command is an ordinary pipeline, ``shlex.split`` parses it the way
    bash does with no special-casing - which the ``>(...)`` form needed.
    """
    hostile_cram = "/data/patient sample/it's a CRAM.cram"
    output_dir = tmp_path / "run one"
    output_dir.mkdir()

    command = _cram_unmapped_command(tmp_path, input_path=hostile_cram, output=str(output_dir))
    tokens = shlex.split(command.removeprefix("set -o pipefail; "))

    assert hostile_cram in tokens, f"the input CRAM did not survive as one operand: {command}"
    assert str(output_dir / "output_unmapped.bam") in tokens, (
        f"the output BAM did not survive as one operand: {command}"
    )


def test_the_merge_runs_immediately_after_the_cram_unmapped_extraction(tmp_path):
    """
    The adjacency that made the race matter, pinned.

    ``samtools merge`` is the next command the stage emits after the extraction, so
    any part of ``unmapped.bam`` not written by the time ``run_command`` returns is
    simply absent from the merged BAM. Nothing sleeps, polls or re-checks in between.
    """
    config = {**CONFIG, "tools": {**CONFIG["tools"], "samtools": "/envs/vntyper/bin/samtools"}}
    commands = _run_bam_to_fastq(
        tmp_path, fast_mode=False, file_format="cram", input_path="/data/sample.cram", config=config
    )

    extraction = next(i for i, c in enumerate(commands) if "-f 4" in c)
    assert commands[extraction + 1] == (
        f"/envs/vntyper/bin/samtools merge -f -@ 4 {tmp_path}/output_sliced_unmapped.bam "
        f"{tmp_path}/output_sliced.bam {tmp_path}/output_unmapped.bam"
    ), f"the merge must be the very next command: {commands[extraction + 1]}"


def test_the_bed_file_branch_passes_minus_l_instead_of_a_region(tmp_path):
    """A BED file replaces the region string in the slice command."""
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t155158000\t155163000\n")

    commands = _run_bam_to_fastq(tmp_path, bed_file=bed)

    assert commands[0] == (
        f"samtools view -P -b -@ 4 -X /data/sample.bam /data/sample.bam.bai -L {bed} "
        f"-o {tmp_path}/output_sliced.bam && "
        f"samtools index -@ 4 {tmp_path}/output_sliced.bam"
    )
    assert REGION not in commands[0]


def test_a_bam_path_with_a_space_survives_every_emitted_command(tmp_path):
    """
    The end-to-end quoting property, through the real stage rather than a builder.

    ``/data/patient sample.bam`` is an ordinary path on any real filesystem. Before
    the builders it split into two arguments and samtools read ``/data/patient``.
    """
    hostile_bam = "/data/patient sample/it's a BAM.bam"
    output_dir = tmp_path / "run one"
    output_dir.mkdir()

    commands = _run_bam_to_fastq(tmp_path, input_path=hostile_bam, fast_mode=False, output=str(output_dir))

    assert commands, "the stage emitted no commands; this test would pass vacuously"
    for command in commands:
        tokens = shlex.split(command.replace("set -o pipefail; ", ""))
        assert all(" " not in token or token.startswith("/") for token in tokens), (
            f"a space leaked out of quoting in: {command}"
        )
    assert any(hostile_bam in shlex.split(c.replace("set -o pipefail; ", "")) for c in commands), (
        "the hostile input path never survived as a single argument"
    )


# ---------------------------------------------------------------------------
# calculate_vntr_coverage - the C1 seam, end to end
# ---------------------------------------------------------------------------


def test_calculate_vntr_coverage_writes_the_frozen_tsv_schema(tmp_path):
    """
    Contract C1 through the real call site.

    ``generate_report.py`` reads ``region_length``, ``uncovered_bases`` and
    ``percent_uncovered`` from this file with ``0`` defaults, so a renamed column
    makes the report state zero uncovered bases for a sample with no coverage at
    all.

    These bytes changed with #171; they are the region-wide statistics, not the
    covered-position ones. Three covered positions ``[10, 20, 30]`` in a 1501 bp
    region means the base set is ``[10, 20, 30] + [0] * 1498``, so ``mean`` is
    ``60 / 1501 = 0.0399733...`` (``0.04``), ``median`` is one of the restored zeros,
    ``stdev`` is ``sqrt((1400 - 1501 * mean**2) / 1500) = 0.965263...`` (``0.97``) and
    ``min`` is ``0``. It read ``20.00 20.00 10.00 10 30 1501 1498 99.80`` before.

    The ninth field is #172's verdict. ``CONFIG`` above declares no ``thresholds``
    block at all - which is legitimate, because ``--config-path`` replaces the config
    rather than merging it - so the shipped defaults of 100 and 50.0 apply, and a mean
    of 0.04 with 99.80% uncovered fails both.
    """
    depth_file = tmp_path / "cov_vntr_coverage.txt"

    def fake_run_command(command, log_file, critical=False, cwd=None):
        assert command == (f"samtools depth -a -@ 4 -r chr1:155160500-155162000 /data/sample.bam > {depth_file}")
        Path(log_file).write_text("")
        depth_file.write_text("chr1\t155160500\t10\nchr1\t155160501\t20\nchr1\t155160502\t30\n")
        return True

    with patch.object(fastq_bam_processing, "run_command", fake_run_command):
        stats = fastq_bam_processing.calculate_vntr_coverage(
            bam_file="/data/sample.bam",
            region="chr1:155160500-155162000",
            threads=4,
            config=CONFIG,
            output_dir=str(tmp_path),
            output_name="cov",
        )

    assert (tmp_path / "cov_summary.tsv").read_text() == (
        "mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n"
        "0.04\t0.00\t0.97\t0\t30\t1501\t1498\t99.80\tFAIL\n"
    )
    assert stats["mean"] == pytest.approx(60 / 1501)
    assert stats["min"] == 0
    assert stats["region_length"] == 1501
    assert stats["uncovered_bases"] == 1498
    assert stats["percent_uncovered"] == pytest.approx(99.80013324450367)
    assert stats["coverage_qc"] == "FAIL"


def test_an_empty_depth_file_aborts_the_coverage_stage(tmp_path):
    """
    A region that matched nothing must not become ``mean = 0``.

    This is the silent-false-negative shape: a wrong contig name yields an empty
    depth table, and a zero mean would be reported as a real measurement.
    """
    depth_file = tmp_path / "cov_vntr_coverage.txt"

    def fake_run_command(command, log_file, critical=False, cwd=None):
        Path(log_file).write_text("")
        depth_file.write_text("")
        return True

    with (
        patch.object(fastq_bam_processing, "run_command", fake_run_command),
        pytest.raises(RuntimeError, match="No coverage data found"),
    ):
        fastq_bam_processing.calculate_vntr_coverage(
            bam_file="/data/sample.bam",
            region="chr1:155160500-155162000",
            threads=4,
            config=CONFIG,
            output_dir=str(tmp_path),
            output_name="cov",
        )


# ---------------------------------------------------------------------------
# align_and_sort_fastq
# ---------------------------------------------------------------------------


def test_align_and_sort_emits_the_pinned_pipeline_with_pipefail(tmp_path):
    """
    The three-stage alignment pipe, with the ``pipefail`` fix.

    Without it a ``bwa mem`` that aborted part-way still let ``samtools sort``
    write a valid but incomplete BAM and exit 0, and everything downstream
    genotyped a fraction of the reads.
    """
    reference = tmp_path / "ref.fa"
    reference.touch()
    for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"):
        (tmp_path / f"ref.fa{ext}").touch()
    sorted_bam = tmp_path / "out" / "output_sorted.bam"
    recorder = _Recorder()

    def run_and_create(command, log_file, critical=False, cwd=None):
        recorder(command, log_file, critical, cwd)
        sorted_bam.parent.mkdir(parents=True, exist_ok=True)
        sorted_bam.touch()
        sorted_bam.with_suffix(".bam.bai").touch()
        return True

    with patch.object(alignment_processing, "run_command", run_and_create):
        result = alignment_processing.align_and_sort_fastq(
            fastq1=Path("/data/r1.fq.gz"),
            fastq2=Path("/data/r2.fq.gz"),
            reference=reference,
            output_dir=tmp_path / "out",
            output_name="output",
            threads=4,
            config=CONFIG,
        )

    assert result == str(sorted_bam)
    assert recorder.commands == [
        f"set -o pipefail; bwa mem -t 4 {reference} /data/r1.fq.gz /data/r2.fq.gz | "
        f"samtools view -@ 4 -b | "
        f"samtools sort -@ 4 -o {sorted_bam}",
        f"samtools index -@ 4 {sorted_bam}",
    ]


def test_align_and_sort_emits_the_single_input_bwa_form(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.touch()
    for ext in (".amb", ".ann", ".bwt", ".pac", ".sa"):
        (tmp_path / f"ref.fa{ext}").touch()
    sorted_bam = tmp_path / "out" / "output_sorted.bam"
    recorder = _Recorder()

    def run_and_create(command, log_file, critical=False, cwd=None):
        recorder(command, log_file, critical, cwd)
        sorted_bam.parent.mkdir(parents=True, exist_ok=True)
        sorted_bam.touch()
        sorted_bam.with_suffix(".bam.bai").touch()
        return True

    with patch.object(alignment_processing, "run_command", run_and_create):
        result = alignment_processing.align_and_sort_fastq(
            fastq1=Path("/data/single.fq.gz"),
            fastq2=None,
            reference=reference,
            output_dir=tmp_path / "out",
            output_name="output",
            threads=4,
            config=CONFIG,
        )

    assert result == str(sorted_bam)
    assert recorder.commands[0] == (
        f"set -o pipefail; bwa mem -t 4 {reference} /data/single.fq.gz | "
        f"samtools view -@ 4 -b | samtools sort -@ 4 -o {sorted_bam}"
    )


def test_an_invalid_region_is_refused_before_samtools_runs(tmp_path, monkeypatch):
    """#224: the ordering is the fix, not the validation.

    ``samtools depth -a`` writes one row per *declared* base, so an oversized region exhausts
    disk before Python allocates anything. Validating after ``run_command`` would catch the
    bad value only once the damage was done. This asserts the command was never built or run --
    coverage alone cannot prove the ordering, only that the line executed.
    """
    launched = []
    monkeypatch.setattr(fastq_bam_processing, "run_command", lambda *a, **k: launched.append(a) or True)

    with pytest.raises(ValueError, match="above the"):
        fastq_bam_processing.calculate_vntr_coverage(
            bam_file=tmp_path / "in.bam",
            region="chr1:1-250000000",
            threads=1,
            config={"tools": {"samtools": "samtools"}},
            output_dir=tmp_path,
            output_name="coverage",
        )

    assert launched == [], "samtools ran for a region that should have been refused first"


def test_a_reversed_region_is_refused_before_samtools_runs(tmp_path, monkeypatch):
    """`samtools depth -r chr1:2000-1000` exits 1 anyway, but the refusal should not depend on
    that -- ``parse_region_length`` is public and reachable without samtools at all."""
    launched = []
    monkeypatch.setattr(fastq_bam_processing, "run_command", lambda *a, **k: launched.append(a) or True)

    with pytest.raises(ValueError, match="is before start"):
        fastq_bam_processing.calculate_vntr_coverage(
            bam_file=tmp_path / "in.bam",
            region="chr1:2000-1000",
            threads=1,
            config={"tools": {"samtools": "samtools"}},
            output_dir=tmp_path,
            output_name="coverage",
        )

    assert launched == []


def test_an_unparseable_region_never_reaches_samtools(tmp_path, monkeypatch):
    """#224: the span bound is worthless if a region without an end can walk past it.

    ``parse_region_length`` degrades an unparseable region to 0 and only warns, which is its
    documented contract. But ``samtools -r`` accepts a bare contig, so passing ``chr1`` through
    unchecked emits chromosome-wide depth -- recreating the exact disk exhaustion the bound
    exists to prevent, through the one input the bound cannot measure.
    """
    launched = []
    monkeypatch.setattr(fastq_bam_processing, "run_command", lambda *a, **k: launched.append(a) or True)

    with pytest.raises(ValueError, match="could not be parsed into coordinates"):
        fastq_bam_processing.calculate_vntr_coverage(
            bam_file=tmp_path / "in.bam",
            region="chr1",
            threads=1,
            config={"tools": {"samtools": "samtools"}},
            output_dir=tmp_path,
            output_name="coverage",
        )

    assert launched == [], "a region with no end reached samtools"
