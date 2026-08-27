"""
Contract tests for :mod:`vntyper.scripts.command_builders`.

Every external tool in the BAM/CRAM/FASTQ path is driven by a **shell string**
handed to :func:`vntyper.scripts.utils.run_command`, which runs it under
``/bin/bash`` with ``shell=True`` (trap 9 - the CRAM unmapped-read branch needs
process substitution, which no argv list can express). Because the command is a
string, the only place quoting can happen is at construction. That is what this
module pins.

Three properties are pinned, each because breaking it is silent:

* **Every interpolated value is quoted.** A sample name, output directory or BAM
  path is caller-supplied text that reaches a shell, so an unquoted one would
  split on whitespace or be read as syntax. The web service constrains names at
  its own boundary; this is the complement, not the alternative, because the CLI
  has no such boundary at all.
* **The CRAM path uses one samtools.** The outer command and the ``>(...)``
  substitution must both call ``{samtools_path}``: under ``mamba run`` a bare
  ``samtools`` resolves to a different binary, so unmapped-read extraction would
  run something other than the rest of the pipeline, or nothing.
* **Pipelines set ``pipefail``.** ``samtools sort -n | samtools fastq`` otherwise
  reports the exit status of ``samtools fastq`` alone, so a ``sort`` that died
  half-way exits 0 and the pipeline continues with a **truncated FASTQ**,
  genotyping a fraction of the reads. That is the silent-wrong-answer class this
  effort exists to remove.

The expected strings below are asserted byte for byte, so the extraction of these
builders is provably behaviour-preserving apart from the three properties above.

Nothing here starts a real process.
"""

import shlex

import pytest

from vntyper.scripts import command_builders
from vntyper.scripts.command_builders import (
    PIPEFAIL_PREFIX,
    build_bam_to_fastq_command,
    build_bwa_align_sort_command,
    build_cram_unmapped_filter_command,
    build_cram_unmapped_indexed_command,
    build_fastp_command,
    build_samtools_depth_command,
    build_samtools_index_command,
    build_samtools_merge_command,
    build_samtools_slice_command,
)

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit

SAMTOOLS = "samtools"
BWA = "bwa"
FASTP = "fastp"

#: A path that breaks every naive f-string: a space splits the argument in two,
#: the apostrophe breaks single-quoting and the double quote breaks double-quoting.
#: `shlex.quote` handles all three; nothing else does.
HOSTILE_DIR = '/data/patient sample/it\'s "weird"'
HOSTILE_NAME = "sample one"


def _tokens(command: str) -> list[str]:
    """
    Re-parse a built command the way bash would, for round-trip assertions.

    ``shlex.split`` does not understand bash process substitution: it turns
    ``>(samtools view ...)`` into the token ``>(samtools`` and leaves the closing
    paren stuck to the last argument, which would make a correctly quoted path
    look unquoted. The substitution is therefore split out and both halves are
    tokenised separately.

    Args:
        command (str): A built shell command.

    Returns:
        list[str]: The arguments bash would produce, with the process
        substitution flattened into the same list.
    """
    if ">(" not in command:
        return shlex.split(command)

    before, rest = command.split(">(", 1)
    inner, after = rest.rsplit(")", 1)
    return shlex.split(before) + shlex.split(inner) + shlex.split(after)


# ---------------------------------------------------------------------------
# Byte-exact pins - captured from the pre-extraction code
# ---------------------------------------------------------------------------


def test_the_fastp_command_is_pinned_with_both_flags_on():
    """
    Deduplication is serialized even when the caller requests more workers.

    The double space before ``--disable_adapter_trimming`` is not a typo: the base
    command ends with a trailing space and each flag is appended with a leading
    one. It is preserved deliberately so the extraction is byte-identical.
    """
    command = build_fastp_command(
        fastp_path=FASTP,
        threads=4,
        fastq_1="/data/in_R1.fq.gz",
        fastq_2="/data/in_R2.fq.gz",
        output="/out",
        output_name="output",
        compression_level=6,
        qualified_quality_phred=20,
        dup_calc_accuracy=3,
        length_required=50,
        disable_adapter_trimming=True,
        deduplication=True,
    )

    assert command == (
        "fastp --thread 1 --in1 /data/in_R1.fq.gz --in2 /data/in_R2.fq.gz "
        "--out1 /out/output_R1.fastq.gz --out2 /out/output_R2.fastq.gz "
        "--compression 6 "
        "--qualified_quality_phred 20 "
        "--dup_calc_accuracy 3 "
        "--length_required 50 "
        "--html /out/output.html "
        "--json /out/output.json "
        " --disable_adapter_trimming --dedup"
    )


def test_the_fastp_command_omits_both_optional_flags_when_disabled():
    """Without deduplication, fastp retains the caller's requested concurrency."""
    command = build_fastp_command(
        fastp_path=FASTP,
        threads=4,
        fastq_1="/data/in_R1.fq.gz",
        fastq_2="/data/in_R2.fq.gz",
        output="/out",
        output_name="output",
        compression_level=6,
        qualified_quality_phred=20,
        dup_calc_accuracy=3,
        length_required=50,
        disable_adapter_trimming=False,
        deduplication=False,
    )

    assert command.startswith("fastp --thread 4 ")
    assert command.endswith("--json /out/output.json ")
    assert "--disable_adapter_trimming" not in command
    assert "--dedup" not in command


def test_the_fastp_command_uses_its_single_input_form_when_mate_two_is_absent():
    command = build_fastp_command(
        fastp_path=FASTP,
        threads=4,
        fastq_1="/data/single.fq.gz",
        fastq_2=None,
        output="/out",
        output_name="output",
        compression_level=6,
        qualified_quality_phred=20,
        dup_calc_accuracy=3,
        length_required=50,
        disable_adapter_trimming=False,
        deduplication=False,
    )

    assert command == (
        "fastp --thread 4 --in1 /data/single.fq.gz --out1 /out/output_R1.fastq.gz "
        "--compression 6 --qualified_quality_phred 20 --dup_calc_accuracy 3 "
        "--length_required 50 --html /out/output.html --json /out/output.json "
    )
    assert "--in2" not in command
    assert "--out2" not in command


@pytest.mark.parametrize(
    ("disable_adapter_trimming", "deduplication", "expected_threads", "expected_suffix"),
    [
        (True, True, 1, " --disable_adapter_trimming --dedup"),
        (True, False, 4, " --disable_adapter_trimming"),
        (False, True, 1, " --dedup"),
        (False, False, 4, ""),
    ],
)
def test_every_fastp_flag_combination_is_pinned(
    disable_adapter_trimming, deduplication, expected_threads, expected_suffix
):
    """Only deduplication serializes fastp across all four flag combinations."""
    command = build_fastp_command(
        fastp_path=FASTP,
        threads=4,
        fastq_1="/a_1.fq",
        fastq_2="/a_2.fq",
        output="/o",
        output_name="s",
        compression_level=6,
        qualified_quality_phred=20,
        dup_calc_accuracy=3,
        length_required=50,
        disable_adapter_trimming=disable_adapter_trimming,
        deduplication=deduplication,
    )

    assert command.startswith(f"fastp --thread {expected_threads} ")
    assert command.endswith("--json /o/s.json " + expected_suffix)


def test_the_slice_command_with_a_region_is_pinned():
    """The region branch: ``view -P -b`` then ``index``, joined by ``&&``."""
    command = build_samtools_slice_command(
        samtools_path=SAMTOOLS,
        in_bam="/data/sample.bam",
        output_bam="/out/output_sliced.bam",
        region="chr1:155158000-155163000",
    )

    assert command == (
        "samtools view -P -b /data/sample.bam chr1:155158000-155163000 -o /out/output_sliced.bam && "
        "samtools index /out/output_sliced.bam"
    )


def test_the_slice_command_with_a_bed_file_is_pinned():
    """The BED branch passes ``-L`` instead of a region string."""
    command = build_samtools_slice_command(
        samtools_path=SAMTOOLS,
        in_bam="/data/sample.bam",
        output_bam="/out/output_sliced.bam",
        bed_file="/data/regions.bed",
    )

    assert command == (
        "samtools view -P -b /data/sample.bam -L /data/regions.bed -o /out/output_sliced.bam && "
        "samtools index /out/output_sliced.bam"
    )


def test_the_slice_command_needs_either_a_region_or_a_bed_file():
    """Neither given is a programming error, not a command that slices everything."""
    with pytest.raises(ValueError, match="region"):
        build_samtools_slice_command(
            samtools_path=SAMTOOLS,
            in_bam="/data/sample.bam",
            output_bam="/out/output_sliced.bam",
        )


def test_the_index_command_is_pinned():
    command = build_samtools_index_command(samtools_path=SAMTOOLS, bam_file="/out/output_sliced.bam")

    assert command == "samtools index /out/output_sliced.bam"


def test_the_index_command_takes_an_output_path():
    """#210/#162: the index must be buildable somewhere other than beside the input.

    ``-o`` requires samtools >= 1.15; ``conda/environment_vntyper.yml`` pins 1.20.
    """
    command = build_samtools_index_command(
        samtools_path=SAMTOOLS, bam_file="/in/sample.bam", output_bai="/out/sample.bam.bai"
    )

    assert command == "samtools index -o /out/sample.bam.bai /in/sample.bam"


def test_the_index_command_without_an_output_path_is_unchanged():
    """The default stays byte-identical, so the two callers indexing a BAM they
    produced in the output directory need no change at all."""
    command = build_samtools_index_command(samtools_path=SAMTOOLS, bam_file="/in/sample.bam")

    assert command == "samtools index /in/sample.bam"


def test_the_index_output_path_is_quoted():
    """The destination is as caller-supplied as the input, so it is quoted too."""
    command = build_samtools_index_command(
        samtools_path=SAMTOOLS, bam_file="/in/sample.bam", output_bai="/out/has space.bai"
    )

    assert shlex.split(command) == ["samtools", "index", "-o", "/out/has space.bai", "/in/sample.bam"]


def test_the_index_output_path_is_the_operand_after_minus_o():
    """The operand order is load-bearing and easy to invert silently.

    ``samtools index -o A B`` writes the index for **B** into **A**. Swapping them
    is not a syntax error -- samtools would try to index the destination -- so the
    positions are pinned rather than only the token set.
    """
    tokens = shlex.split(
        build_samtools_index_command(samtools_path=SAMTOOLS, bam_file="/in/sample.bam", output_bai="/out/sample.bai")
    )

    assert tokens[tokens.index("-o") + 1] == "/out/sample.bai"
    assert tokens[-1] == "/in/sample.bam"


# ---------------------------------------------------------------------------
# Milestone 4 - CRAM references and bounded unmapped extraction
# ---------------------------------------------------------------------------


def test_thread_flag_omits_samtools_threads_at_one_or_less():
    """A one-thread command retains the historical syntax without ``-@``."""
    assert command_builders._thread_flag(1) == ""
    assert command_builders._thread_flag(0) == ""


def test_thread_flag_adds_samtools_threads_above_one():
    """Multiple threads must reach samtools instead of only the callers' logs."""
    assert command_builders._thread_flag(8) == "-@ 8 "


def test_reference_flag_omits_none_and_quotes_the_path():
    """The optional reference is absent by default but remains one shell argument when set."""
    assert command_builders._reference_flag(None) == ""
    assert command_builders._reference_flag("/r/my genome.fa") == "-T '/r/my genome.fa' "


def test_slice_carries_threads_and_can_skip_indexing():
    """The preflight probe's slice leaves no unwanted index beside its throwaway output."""
    command = build_samtools_slice_command(
        samtools_path=SAMTOOLS,
        in_bam="/o/view.cram",
        output_bam="/o/s.bam",
        region="chr1:1-2",
        threads=8,
        index_output=False,
    )

    assert "-@ 8" in command
    assert "&&" not in command
    assert " index " not in command


def test_slice_uses_the_exact_custom_index_operand() -> None:
    """Dropping ``-X`` would make htslib rediscover a replaced public index path."""
    command = build_samtools_slice_command(
        samtools_path=SAMTOOLS,
        in_bam="/o/view.bam",
        index_path="/proc/123/fd/9",
        output_bam="/o/s.bam",
        region="chr1:1-2",
        index_output=False,
    )

    tokens = shlex.split(command)
    position = tokens.index("-X")
    assert tokens[position + 1 : position + 3] == ["/o/view.bam", "/proc/123/fd/9"]


def test_slice_indexes_by_default():
    """Ordinary conversion keeps creating the index expected by downstream stages."""
    command = build_samtools_slice_command(
        samtools_path=SAMTOOLS,
        in_bam="/o/view.cram",
        output_bam="/o/s.bam",
        region="chr1:1-2",
    )

    assert "&& samtools index /o/s.bam" in command


def test_index_carries_threads_but_one_thread_preserves_the_old_command():
    """Indexing is threaded only when doing so can increase parallelism."""
    threaded = build_samtools_index_command(samtools_path=SAMTOOLS, bam_file="/o/s.bam", threads=8)
    single_thread = build_samtools_index_command(samtools_path=SAMTOOLS, bam_file="/o/s.bam", threads=1)

    assert threaded == "samtools index -@ 8 /o/s.bam"
    assert single_thread == "samtools index /o/s.bam"


def test_depth_quotes_a_reference_path_with_spaces_and_shell_metacharacters():
    """Depth uses its supported long reference flag without exposing a shell injection path."""
    reference = "/r/my genome; touch /tmp/pwned.fa"
    command = build_samtools_depth_command(
        samtools_path=SAMTOOLS,
        threads=4,
        region="chr1:1-2",
        bam_file="/o/view.cram",
        coverage_output="/o/d.txt",
        reference_path=reference,
    )

    assert "--reference '/r/my genome; touch /tmp/pwned.fa'" in command
    assert "-T" not in command
    tokens = shlex.split(command)
    assert tokens[tokens.index("--reference") + 1] == reference
    assert "touch" not in tokens


def test_a_reference_path_with_a_space_is_quoted_not_split():
    """The builder owns shell quoting because its command runs with ``shell=True``."""
    command = build_samtools_slice_command(
        samtools_path=SAMTOOLS,
        in_bam="/o/v.cram",
        output_bam="/o/s.bam",
        region="chr1:1-2",
        reference_path="/r/my genome.fa",
    )

    assert "-T '/r/my genome.fa'" in command


def test_cram_filter_carries_the_reference_and_omits_threading_at_one():
    """The whole-file fallback must decode CRAM without adding needless ``-@ 1`` flags."""
    command = build_cram_unmapped_filter_command(
        samtools_path=SAMTOOLS,
        in_bam="/o/view.cram",
        unmapped_bam="/o/unmapped.bam",
        threads=1,
        reference_path="/r/g.fa",
    )

    assert "-T /r/g.fa" in command
    assert "-@" not in command


def test_the_indexed_unmapped_fetch_uses_star_without_a_pipe():
    """The bounded fast path must fetch literal unplaced reads, not stream the CRAM."""
    command = command_builders.build_cram_unmapped_indexed_command(
        samtools_path=SAMTOOLS,
        in_bam="/o/view.cram",
        unmapped_bam="/o/unmapped.bam",
        threads=4,
        reference_path="/r/g.fa",
    )

    assert "'*'" in command
    assert "-f 4" in command
    assert "|" not in command


def test_the_indexed_unmapped_fetch_uses_the_exact_custom_index() -> None:
    """Literal-star retrieval must not fall back to index discovery by basename."""
    command = command_builders.build_cram_unmapped_indexed_command(
        samtools_path=SAMTOOLS,
        in_bam="/o/view.cram",
        index_path="/proc/123/fd/10",
        unmapped_bam="/o/unmapped.bam",
        threads=4,
    )

    tokens = shlex.split(command)
    position = tokens.index("-X")
    assert tokens[position + 1 : position + 3] == ["/o/view.cram", "/proc/123/fd/10"]


@pytest.mark.parametrize(
    "builder",
    [build_cram_unmapped_filter_command, command_builders.build_cram_unmapped_indexed_command],
)
def test_every_unmapped_builder_selects_read_unmapped_flag_four_for_single_end(builder):
    """Unpaired unmapped reads have flag 4 but not paired-only flag 12."""
    command = builder(
        samtools_path=SAMTOOLS,
        in_bam="/o/view.cram",
        unmapped_bam="/o/unmapped.bam",
        threads=2,
        reference_path="/r/g.fa",
    )

    assert "-f 4" in command
    assert "-f 12" not in command


def test_the_probe_never_combines_P_with_c_because_samtools_rejects_that():
    """``samtools view`` rejects the apparently natural ``-P -c`` combination."""
    probe = command_builders.build_cram_reference_probe_command(
        samtools_path=SAMTOOLS,
        in_bam="/o/view.cram",
        region="chr1:1-2",
        reference_path="/r/g.fa",
    )

    assert "-P" in probe
    assert " -c" not in probe, "samtools: The options -P and -c cannot be combined"
    assert "-b" not in probe, "the probe discards SAM output and need not encode a BAM"
    assert "-o /dev/null" in probe


def test_the_stream_reference_probe_decodes_the_whole_cram_without_a_target():
    """A target-only probe cannot authorize the later whole-file stream consumer."""
    probe = command_builders.build_cram_stream_reference_probe_command(
        samtools_path=SAMTOOLS,
        in_bam="/o/view.cram",
        reference_path="/r/full genome.fa",
        threads=4,
    )

    assert probe == "samtools view -@ 4 -T '/r/full genome.fa' -h /o/view.cram -o /dev/null"
    assert " -P " not in f" {probe} "
    assert "chr" not in probe


def test_the_probe_has_the_same_shape_as_the_slice_it_authorises():
    """A passing probe must exercise the flags the subsequent slice will use."""
    kwargs = {"samtools_path": SAMTOOLS, "in_bam": "/o/view.cram", "reference_path": "/r/g.fa"}
    probe = command_builders.build_cram_reference_probe_command(bed_file="/o/r.bed", **kwargs)
    sliced = build_samtools_slice_command(output_bam="/o/s.bam", bed_file="/o/r.bed", **kwargs)

    for flag in ("-P", "-L /o/r.bed", "-T /r/g.fa"):
        assert flag in probe and flag in sliced


def test_the_probe_needs_a_region_or_bed_file():
    """A reference check without a target would decode the entire alignment."""
    with pytest.raises(ValueError, match="region"):
        command_builders.build_cram_reference_probe_command(
            samtools_path=SAMTOOLS,
            in_bam="/o/view.cram",
            reference_path="/r/g.fa",
        )


def test_idxstats_is_threaded_without_a_reference_flag():
    """Index statistics read index/container metadata and must not try to resolve FASTA."""
    command = command_builders.build_samtools_idxstats_command(
        samtools_path=SAMTOOLS,
        in_bam="/o/view.cram",
        threads=8,
    )

    assert command == "samtools idxstats -@ 8 /o/view.cram"
    assert "-T" not in command


def test_indexed_unmapped_count_probes_the_exact_literal_star_consumer():
    """Authorization must count the same flag and region used by indexed recovery."""
    command = command_builders.build_samtools_unmapped_indexed_count_command(
        samtools_path=SAMTOOLS,
        in_bam="/o/view.cram",
        threads=8,
    )

    assert command == "samtools view -c -f 4 -@ 8 /o/view.cram '*'"
    assert "-T" not in command


def test_the_merge_command_is_pinned():
    command = build_samtools_merge_command(
        samtools_path=SAMTOOLS,
        merged_bam="/out/output_sliced_unmapped.bam",
        sliced_bam="/out/output_sliced.bam",
        unmapped_bam="/out/output_unmapped.bam",
        threads=4,
    )

    assert command == (
        "samtools merge -f -@ 4 /out/output_sliced_unmapped.bam /out/output_sliced.bam /out/output_unmapped.bam"
    )


def test_the_depth_command_requests_every_position_in_the_region():
    """#171: without `-a`, samtools emits only covered positions.

    Every statistic downstream is then over the covered bases rather than the region,
    and `uncovered_bases` - which is derived by subtraction - silently reads 0. The flag
    is asserted as a parsed argument, not as a substring, so `-all` could not satisfy it.

    ``samtools depth`` redirects to a file; a redirect is not a pipe, so no pipefail.
    """
    command = build_samtools_depth_command(
        samtools_path=SAMTOOLS,
        threads=4,
        region="chr1:155160500-155162000",
        bam_file="/data/sample.bam",
        coverage_output="/out/cov_vntr_coverage.txt",
    )

    assert command == (
        "samtools depth -a -@ 4 -r chr1:155160500-155162000 /data/sample.bam > /out/cov_vntr_coverage.txt"
    )
    assert "-a" in shlex.split(command.split(">")[0])
    assert PIPEFAIL_PREFIX not in command, "a single command with a redirect already reports its own exit status"


def test_depth_uses_the_exact_custom_index_operand() -> None:
    """Coverage must retain the same index bytes as conversion through its final read."""
    command = build_samtools_depth_command(
        samtools_path=SAMTOOLS,
        threads=4,
        region="chr1:1-2",
        bam_file="/o/view.cram",
        index_path="/proc/123/fd/11",
        coverage_output="/o/depth.txt",
    )

    tokens = shlex.split(command.split(">", maxsplit=1)[0])
    position = tokens.index("-X")
    assert tokens[position + 1 : position + 3] == ["/o/view.cram", "/proc/123/fd/11"]


# ---------------------------------------------------------------------------
# D2 - the CRAM process-substitution binary mismatch
# ---------------------------------------------------------------------------


def test_the_cram_filter_runs_the_configured_samtools():
    """
    The extraction must use the samtools from ``config["tools"]["samtools"]``.

    It used to call a bare ``samtools`` in its writing stage. Under
    ``mamba run -n <env>`` a bare name resolves against a different PATH - so the
    unmapped reads were extracted by a different build of samtools, or by nothing
    at all while the pipeline reported success. With the pipe collapsed there is one
    invocation rather than two (#262), so this now asserts that the single one is the
    configured one and that no bare token survives anywhere.
    """
    configured = "/opt/conda/envs/vntyper/bin/samtools"

    command = build_cram_unmapped_filter_command(
        samtools_path=configured,
        in_bam="/data/sample.cram",
        unmapped_bam="/out/output_unmapped.bam",
        threads=4,
    )

    assert command.startswith(configured), f"the extraction must call {configured}, got: {command[:80]!r}"

    bare_calls = [
        part for part in command.split() if part == "samtools"
    ]  # a bare `samtools` token anywhere means the mismatch is back
    assert bare_calls == [], "no bare `samtools` may remain; every invocation comes from config"


def test_the_cram_filter_command_is_pinned():
    """The whole CRAM command, now one process rather than a pipe (#262).

    The expectation changed because the command genuinely changed: decoding to SAM
    text, piping and re-parsing is strictly more work than decoding once, and it does
    not scale with threads because each stage is single-threaded. The read set is
    identical -- the old first stage applied no filter and no region, so both forms
    select flag 4 over every record in the file.
    """
    command = build_cram_unmapped_filter_command(
        samtools_path=SAMTOOLS,
        in_bam="/data/sample.cram",
        unmapped_bam="/out/output_unmapped.bam",
        threads=4,
    )

    assert command == "samtools view -b -f 4 -@ 4 /data/sample.cram -o /out/output_unmapped.bam"


def test_the_cram_filter_uses_no_process_substitution():
    """
    The writer must be a pipeline stage, not a ``>(...)``, or the shell will not wait.

    bash does not wait for a process substitution to finish, so the old form let
    ``run_command`` return while samtools was still writing ``unmapped.bam`` - and
    ``process_bam_to_fastq`` then merged the partial file. Measured on a synthetic
    600k-read CRAM: 199,797 of 200,000 unmapped reads were still missing when the
    shell returned, and ``samtools merge`` accepted the truncated file with exit 0.
    A silently under-called sample, with no pipeline error.

    ``tee``'s stdout went to ``/dev/null``, so it had exactly one consumer and was
    never needed. A plain pipe puts the writer in the pipeline, so the shell waits
    for it and ``pipefail`` finally means something for this command.
    """
    command = build_cram_unmapped_filter_command(
        samtools_path=SAMTOOLS,
        in_bam="/data/sample.cram",
        unmapped_bam="/out/output_unmapped.bam",
        threads=4,
    )

    assert ">(" not in command, "a process substitution is not waited for; the merge would race it"
    assert " tee " not in command, "tee had a single consumer and is what forced the substitution"


def test_the_cram_filter_no_longer_needs_bash_for_itself():
    """This builder stopped being a reason for trap 9 (#262).

    Trap 9 pinned ``/bin/bash`` first for the process substitution, then for
    ``set -o pipefail``. With the pipe collapsed into one process this command needs
    neither: its exit status is samtools' own, and there is no non-POSIX syntax in it
    at all.

    **Trap 9 itself still holds** -- ``build_bam_to_fastq_command`` and
    ``build_bwa_align_sort_command`` are still pipes and still require pipefail, which
    ``test_every_multi_stage_pipe_sets_pipefail`` enforces over the whole inventory.
    What has changed is only which commands are the reason.
    """
    command = build_cram_unmapped_filter_command(
        samtools_path=SAMTOOLS,
        in_bam="/data/sample.cram",
        unmapped_bam="/out/output_unmapped.bam",
        threads=4,
    )

    assert not command.startswith(PIPEFAIL_PREFIX), "no pipeline means no stage whose status can be masked"
    assert "|" not in command
    assert ">(" not in command


# ---------------------------------------------------------------------------
# D3 - pipefail on every multi-stage pipe
# ---------------------------------------------------------------------------


def test_the_bam_to_fastq_command_is_pinned_with_pipefail():
    """
    Without ``pipefail`` a failed ``samtools sort`` produced a truncated FASTQ.

    The pipeline's exit status was ``samtools fastq``'s alone, and ``fastq`` is
    perfectly happy to consume a short stream and exit 0. The result was a FASTQ
    holding a fraction of the reads, genotyped without a single warning.
    """
    command = build_bam_to_fastq_command(
        samtools_path=SAMTOOLS,
        in_bam="/out/output_sliced.bam",
        threads=4,
        fastq_r1="/out/output_R1.fastq.gz",
        fastq_r2="/out/output_R2.fastq.gz",
        fastq_other="/out/output_other.fastq.gz",
        fastq_single="/out/output_single.fastq.gz",
        sort_tmp_prefix="/out/output_sort_tmp",
    )

    assert command == (
        "set -o pipefail; samtools sort -n -@ 4 -T /out/output_sort_tmp /out/output_sliced.bam | "
        "samtools fastq -@ 4 - -1 /out/output_R1.fastq.gz "
        "-2 /out/output_R2.fastq.gz -0 /out/output_other.fastq.gz "
        "-s /out/output_single.fastq.gz"
    )


def test_the_sort_stage_takes_an_explicit_tmp_prefix_operand():
    """D2: the name sort must spill beside its output, not in the launch CWD."""
    command = build_bam_to_fastq_command(
        samtools_path=SAMTOOLS,
        in_bam="/out/output_sliced.bam",
        threads=4,
        fastq_r1="/out/output_R1.fastq.gz",
        fastq_r2="/out/output_R2.fastq.gz",
        fastq_other="/out/output_other.fastq.gz",
        fastq_single="/out/output_single.fastq.gz",
        sort_tmp_prefix="/out/output_sort_tmp",
    )

    sort_stage = command.removeprefix(PIPEFAIL_PREFIX).split("|")[0]
    tokens = shlex.split(sort_stage)
    assert tokens[tokens.index("-T") + 1] == "/out/output_sort_tmp"


def test_the_bwa_align_sort_command_is_pinned_with_pipefail():
    """A three-stage pipe: a failing ``bwa mem`` must not be masked by ``samtools sort``."""
    command = build_bwa_align_sort_command(
        bwa_path=BWA,
        samtools_path=SAMTOOLS,
        threads=4,
        reference="/ref/ref.fa",
        fastq1="/data/r1.fq.gz",
        fastq2="/data/r2.fq.gz",
        sorted_bam="/out/output_sorted.bam",
    )

    assert command == (
        "set -o pipefail; bwa mem -t 4 /ref/ref.fa /data/r1.fq.gz /data/r2.fq.gz | "
        "samtools view -@ 4 -b | "
        "samtools sort -@ 4 -o /out/output_sorted.bam"
    )


def test_the_bwa_align_sort_command_accepts_one_fastq_without_a_none_operand():
    command = build_bwa_align_sort_command(
        bwa_path=BWA,
        samtools_path=SAMTOOLS,
        threads=4,
        reference="/ref/ref.fa",
        fastq1="/data/single.fq.gz",
        fastq2=None,
        sorted_bam="/out/output_sorted.bam",
    )

    assert command == (
        "set -o pipefail; bwa mem -t 4 /ref/ref.fa /data/single.fq.gz | "
        "samtools view -@ 4 -b | "
        "samtools sort -@ 4 -o /out/output_sorted.bam"
    )
    assert "None" not in command


@pytest.mark.parametrize(
    ("name", "command"),
    [
        (
            "bam_to_fastq",
            build_bam_to_fastq_command(
                samtools_path=SAMTOOLS,
                in_bam="/b.bam",
                threads=1,
                fastq_r1="/1.gz",
                fastq_r2="/2.gz",
                fastq_other="/0.gz",
                fastq_single="/s.gz",
                sort_tmp_prefix="/s_tmp",
            ),
        ),
        (
            "bwa_align_sort",
            build_bwa_align_sort_command(
                bwa_path=BWA,
                samtools_path=SAMTOOLS,
                threads=1,
                reference="/r.fa",
                fastq1="/1.fq",
                fastq2="/2.fq",
                sorted_bam="/o.bam",
            ),
        ),
        (
            "cram_unmapped_filter",
            build_cram_unmapped_filter_command(
                samtools_path=SAMTOOLS,
                in_bam="/s.cram",
                unmapped_bam="/u.bam",
                threads=1,
            ),
        ),
        (
            "cram_unmapped_indexed",
            build_cram_unmapped_indexed_command(
                samtools_path=SAMTOOLS,
                in_bam="/s.cram",
                unmapped_bam="/u.bam",
                threads=1,
            ),
        ),
        (
            "samtools_slice",
            build_samtools_slice_command(
                samtools_path=SAMTOOLS,
                in_bam="/s.cram",
                output_bam="/o.bam",
                region="chr1:1-2",
            ),
        ),
        (
            "samtools_merge",
            build_samtools_merge_command(
                samtools_path=SAMTOOLS,
                merged_bam="/m.bam",
                sliced_bam="/s.bam",
                unmapped_bam="/u.bam",
                threads=1,
            ),
        ),
        (
            "samtools_depth",
            build_samtools_depth_command(
                samtools_path=SAMTOOLS,
                threads=1,
                region="chr1:1-2",
                bam_file="/b.bam",
                coverage_output="/d.txt",
            ),
        ),
    ],
)
def test_pipefail_is_set_on_exactly_the_commands_that_pipe(name, command):
    """A command pipes if and only if it opts into pipefail. No exceptions either way.

    This used to enumerate the three builders that were pipes and assert each one
    carried the prefix, which meant the inventory had to be edited by hand whenever a
    builder stopped being a pipe -- and an entry left behind would fail while an entry
    *removed* would silently stop checking anything.

    Stating it as a biconditional over every shell-producing builder covers both
    directions at once: a new pipe without pipefail fails, and a prefix left on a
    command that no longer pipes fails too. `cram_unmapped_filter` is now the second
    kind (#262), and is deliberately still in this list rather than deleted from it.
    """
    assert ("|" in command) == command.startswith(PIPEFAIL_PREFIX), (
        f"{name}: pipefail and piping must agree, got pipe={'|' in command} "
        f"prefix={command.startswith(PIPEFAIL_PREFIX)}"
    )


# ---------------------------------------------------------------------------
# D1 - quoting, across every builder
# ---------------------------------------------------------------------------


def _hostile_commands() -> dict[str, tuple[str, list[str]]]:
    """
    Build every command with a hostile path, paired with the paths that must survive.

    Returns:
        dict: builder name -> (command string, list of paths that must each be a
        single token after re-parsing).
    """
    r1 = f"{HOSTILE_DIR}/{HOSTILE_NAME}_R1.fastq.gz"
    r2 = f"{HOSTILE_DIR}/{HOSTILE_NAME}_R2.fastq.gz"
    other = f"{HOSTILE_DIR}/{HOSTILE_NAME}_other.fastq.gz"
    single = f"{HOSTILE_DIR}/{HOSTILE_NAME}_single.fastq.gz"
    sort_tmp = f"{HOSTILE_DIR}/{HOSTILE_NAME}_sort_tmp"
    in_bam = f"{HOSTILE_DIR}/{HOSTILE_NAME}.bam"
    in_cram = f"{HOSTILE_DIR}/{HOSTILE_NAME}.cram"
    sliced = f"{HOSTILE_DIR}/{HOSTILE_NAME}_sliced.bam"
    unmapped = f"{HOSTILE_DIR}/{HOSTILE_NAME}_unmapped.bam"
    merged = f"{HOSTILE_DIR}/{HOSTILE_NAME}_sliced_unmapped.bam"
    bed = f"{HOSTILE_DIR}/{HOSTILE_NAME}.bed"
    reference = f"{HOSTILE_DIR}/{HOSTILE_NAME}.fa"
    depth = f"{HOSTILE_DIR}/{HOSTILE_NAME}_vntr_coverage.txt"

    return {
        "fastp": (
            build_fastp_command(
                fastp_path=FASTP,
                threads=4,
                fastq_1=r1,
                fastq_2=r2,
                output=HOSTILE_DIR,
                output_name=HOSTILE_NAME,
                compression_level=6,
                qualified_quality_phred=20,
                dup_calc_accuracy=3,
                length_required=50,
                disable_adapter_trimming=True,
                deduplication=True,
            ),
            [r1, r2, f"{HOSTILE_DIR}/{HOSTILE_NAME}.html", f"{HOSTILE_DIR}/{HOSTILE_NAME}.json"],
        ),
        "slice_region": (
            build_samtools_slice_command(
                samtools_path=SAMTOOLS,
                in_bam=in_bam,
                output_bam=sliced,
                region="chr1:155158000-155163000",
            ),
            [in_bam, sliced, "chr1:155158000-155163000"],
        ),
        "slice_bed": (
            build_samtools_slice_command(
                samtools_path=SAMTOOLS,
                in_bam=in_bam,
                output_bam=sliced,
                bed_file=bed,
            ),
            [in_bam, sliced, bed],
        ),
        "index": (
            build_samtools_index_command(samtools_path=SAMTOOLS, bam_file=sliced),
            [sliced],
        ),
        "merge": (
            build_samtools_merge_command(
                samtools_path=SAMTOOLS,
                merged_bam=merged,
                sliced_bam=sliced,
                unmapped_bam=unmapped,
                threads=4,
            ),
            [merged, sliced, unmapped],
        ),
        "bam_to_fastq": (
            build_bam_to_fastq_command(
                samtools_path=SAMTOOLS,
                in_bam=sliced,
                threads=4,
                fastq_r1=r1,
                fastq_r2=r2,
                fastq_other=other,
                fastq_single=single,
                sort_tmp_prefix=sort_tmp,
            ),
            [sliced, r1, r2, other, single, sort_tmp],
        ),
        "depth": (
            build_samtools_depth_command(
                samtools_path=SAMTOOLS,
                threads=4,
                region="chr1:155160500-155162000",
                bam_file=in_bam,
                coverage_output=depth,
            ),
            [in_bam, depth],
        ),
        "bwa_align_sort": (
            build_bwa_align_sort_command(
                bwa_path=BWA,
                samtools_path=SAMTOOLS,
                threads=4,
                reference=reference,
                fastq1=r1,
                fastq2=r2,
                sorted_bam=sliced,
            ),
            [reference, r1, r2, sliced],
        ),
        "cram_unmapped_filter": (
            build_cram_unmapped_filter_command(
                samtools_path=SAMTOOLS,
                in_bam=in_cram,
                unmapped_bam=unmapped,
                threads=4,
            ),
            [in_cram, unmapped],
        ),
    }


@pytest.mark.parametrize("builder", sorted(_hostile_commands()))
def test_a_path_with_a_space_and_a_quote_survives_as_one_argument(builder):
    """
    Every interpolated path must round-trip through ``shlex.split`` as one token.

    A path containing a space is the realistic case - ``/data/patient sample/`` -
    and it is enough to split one argument into two, so samtools reads the wrong
    file or writes to the wrong place. The apostrophe and double quote in the
    fixture make sure the quoting is `shlex.quote`, not a hand-rolled wrapper in
    one kind of quote.
    """
    command, must_survive = _hostile_commands()[builder]
    tokens = _tokens(command)

    for path in must_survive:
        assert path in tokens, (
            f"{builder}: {path!r} did not survive re-parsing as a single argument.\n"
            f"command: {command}\ntokens:  {tokens}"
        )


@pytest.mark.parametrize("builder", sorted(_hostile_commands()))
def test_a_hostile_path_cannot_inject_a_second_command(builder):
    """
    Metacharacters in a path must not become shell syntax.

    ``run_command`` hands the string to bash, so an unquoted ``;`` or ``$(...)``
    in a sample name is arbitrary command execution. The builders are the only
    place this can be stopped.
    """
    injected = "/data/x; touch /tmp/pwned #"
    r1 = f"{injected}/s_R1.fastq.gz"

    commands = {
        "fastp": lambda: build_fastp_command(
            fastp_path=FASTP,
            threads=1,
            fastq_1=r1,
            fastq_2=r1,
            output=injected,
            output_name="s",
            compression_level=6,
            qualified_quality_phred=20,
            dup_calc_accuracy=3,
            length_required=50,
            disable_adapter_trimming=False,
            deduplication=False,
        ),
        "slice_region": lambda: build_samtools_slice_command(
            samtools_path=SAMTOOLS, in_bam=injected, output_bam=injected, region="chr1:1-2"
        ),
        "slice_bed": lambda: build_samtools_slice_command(
            samtools_path=SAMTOOLS, in_bam=injected, output_bam=injected, bed_file=injected
        ),
        "index": lambda: build_samtools_index_command(samtools_path=SAMTOOLS, bam_file=injected),
        "merge": lambda: build_samtools_merge_command(
            samtools_path=SAMTOOLS, merged_bam=injected, sliced_bam=injected, unmapped_bam=injected, threads=1
        ),
        "bam_to_fastq": lambda: build_bam_to_fastq_command(
            samtools_path=SAMTOOLS,
            in_bam=injected,
            threads=1,
            fastq_r1=injected,
            fastq_r2=injected,
            fastq_other=injected,
            fastq_single=injected,
            sort_tmp_prefix=injected,
        ),
        "depth": lambda: build_samtools_depth_command(
            samtools_path=SAMTOOLS, threads=1, region="chr1:1-2", bam_file=injected, coverage_output=injected
        ),
        "bwa_align_sort": lambda: build_bwa_align_sort_command(
            bwa_path=BWA,
            samtools_path=SAMTOOLS,
            threads=1,
            reference=injected,
            fastq1=injected,
            fastq2=injected,
            sorted_bam=injected,
        ),
        "cram_unmapped_filter": lambda: build_cram_unmapped_filter_command(
            samtools_path=SAMTOOLS, in_bam=injected, unmapped_bam=injected, threads=1
        ),
    }

    command = commands[builder]()
    tokens = _tokens(command)

    assert "touch" not in tokens, f"{builder}: the injected command survived as a token: {command}"
    assert any(injected in token for token in tokens), f"{builder}: the hostile path vanished entirely: {command}"


def test_the_tool_invocation_is_not_quoted_because_it_may_be_a_command_prefix():
    """
    Tool entries in ``config.json`` are command prefixes, not paths.

    ``config["tools"]["advntr"]`` is literally ``"mamba run -n envadvntr advntr"``
    (trap 6), and ``--config-path`` replaces the whole config, so an operator can
    legitimately set ``samtools`` to a multi-word invocation. Running
    ``shlex.quote`` over it would turn the whole thing into one argument and bash
    would look for a binary with a space in its name.

    The split is deliberate: the **tool** side is operator-controlled
    configuration, the **path** side is user-controlled input. Only the latter is
    quoted.
    """
    command = build_samtools_index_command(
        samtools_path="mamba run -n vntyper samtools", bam_file="/out/output_sliced.bam"
    )

    assert command == "mamba run -n vntyper samtools index /out/output_sliced.bam"
    assert _tokens(command)[:5] == ["mamba", "run", "-n", "vntyper", "samtools"]


# ---------------------------------------------------------------------------
# Uncompressed intermediates (#262)
# ---------------------------------------------------------------------------
#
# `-u` is BGZF level 0: still a valid, indexable BAM, just not deflated. It is worth
# taking only where the file is re-read within milliseconds and then deleted or
# replaced, which is true of the region slice in non-fast mode and of the unmapped
# extraction, and false of everything else.


def test_the_slice_command_can_request_uncompressed_output():
    """The non-fast slice is consumed by the merge and then replaced by it."""
    command = build_samtools_slice_command(
        samtools_path=SAMTOOLS,
        in_bam="/data/sample.bam",
        output_bam="/out/output_sliced.bam",
        region="chr1:1-2",
        index_output=False,
        uncompressed=True,
    )

    assert command == "samtools view -P -b -u /data/sample.bam chr1:1-2 -o /out/output_sliced.bam"


def test_the_slice_command_is_compressed_by_default():
    """A default that silently wrote level 0 would reach files that survive the run."""
    command = build_samtools_slice_command(
        samtools_path=SAMTOOLS,
        in_bam="/data/sample.bam",
        output_bam="/out/output_sliced.bam",
        region="chr1:1-2",
        index_output=False,
    )

    assert "-u" not in _tokens(command)


@pytest.mark.parametrize(
    "builder",
    [build_cram_unmapped_filter_command, build_cram_unmapped_indexed_command],
)
def test_both_unmapped_builders_can_request_uncompressed_output(builder):
    """The unmapped BAM is merged milliseconds later and then deleted.

    Both scan modes write the same throwaway file, so a parameter on only one of them
    would make the saving depend on which scan preflight happened to prove.
    """
    command = builder(
        samtools_path=SAMTOOLS,
        in_bam="/data/sample.bam",
        unmapped_bam="/out/output_unmapped.bam",
        threads=4,
        uncompressed=True,
    )

    assert "-u" in _tokens(command)


@pytest.mark.parametrize(
    "builder",
    [build_cram_unmapped_filter_command, build_cram_unmapped_indexed_command],
)
def test_both_unmapped_builders_are_compressed_by_default(builder):
    command = builder(
        samtools_path=SAMTOOLS,
        in_bam="/data/sample.bam",
        unmapped_bam="/out/output_unmapped.bam",
        threads=4,
    )

    assert "-u" not in _tokens(command)


# ---------------------------------------------------------------------------
# One process, not a decode-to-SAM-text-and-reparse pipe (#262)
# ---------------------------------------------------------------------------


def test_unmapped_stream_extraction_is_a_single_process():
    """One samtools view replaces the decode-to-SAM-text-and-reparse pipe.

    The old shape decoded the whole alignment to SAM **text**, piped it, and re-parsed
    it. Measured on the 963,549-read 7a61 fixture with a warm page cache, three runs
    each: 0.12 s to 0.08 s at -@4 and 0.13 s to 0.08 s at -@8, peak RSS 10.6 MB to
    6.1 MB and 16.5 MB to 8.2 MB. The pipe gains nothing from the extra threads --
    SAM serialise and parse are each single-threaded per stage -- while the single
    view holds 0.08 s at both.
    """
    command = build_cram_unmapped_filter_command(
        samtools_path=SAMTOOLS,
        in_bam="/o/view.cram",
        unmapped_bam="/o/unmapped.bam",
        threads=4,
        reference_path="/r/g.fa",
        uncompressed=True,
    )

    assert command == "samtools view -b -f 4 -u -T /r/g.fa -@ 4 /o/view.cram -o /o/unmapped.bam"


def test_the_collapsed_extraction_needs_no_pipefail():
    """With no pipe there is no masked stage: the exit status is samtools' own.

    The `tee >(...)` flush hazard this builder used to document cannot occur either --
    it needed a writer the shell does not wait for, and there is no longer a second
    process at all.
    """
    command = build_cram_unmapped_filter_command(
        samtools_path=SAMTOOLS,
        in_bam="/o/view.cram",
        unmapped_bam="/o/unmapped.bam",
        threads=4,
    )

    assert "|" not in command
    assert PIPEFAIL_PREFIX not in command
    assert command.count("samtools") == 1


def test_the_collapsed_extraction_still_reads_the_whole_file():
    """It must not acquire a region or an index query: that would change the read set.

    A `'*'` query returns only *unplaced* unmapped reads and drops placed ones -- flag
    4 parked at a mapped mate's coordinate -- which is measurably wrong on this data:
    329, 3,732 and 129 reads on the b178, 6449 and 7a61 fixtures respectively --
    measured directly as `samtools view -c -f 4 <bam>` against
    `samtools view -c -f 4 <bam> '*'`, and as high as 5,806 on 6c28.
    """
    command = build_cram_unmapped_filter_command(
        samtools_path=SAMTOOLS,
        in_bam="/o/view.cram",
        unmapped_bam="/o/unmapped.bam",
        threads=4,
    )

    assert "'*'" not in command
    assert "-X" not in _tokens(command)
    assert "chr" not in command


def test_the_merged_bam_can_never_be_written_uncompressed():
    """It is renamed to <name>_sliced.bam, survives the run, and is shipped to users.

    An earlier draft of this work applied `-u` to the merge, which was the largest
    single saving available here. It is forfeited deliberately: the merge is
    renamed to `<name>_sliced.bam`; the intermediate cleanup removes only
    `<name>_unmapped.bam`, and only when `delete_intermediates` is true, whose CLI
    default is False; `<name>_sliced.bam` is an enumerated artifact in
    artifact_names.py; and both the CLI archive and the web service archive include the
    whole output directory recursively. Writing it uncompressed would ship a roughly 3x
    larger BAM -- gigabytes on a WGS input -- into user-downloadable archives.

    The builder therefore grows no `uncompressed` parameter at all, so no caller can
    reintroduce this by passing True.
    """
    import inspect

    assert "uncompressed" not in inspect.signature(build_samtools_merge_command).parameters

    command = build_samtools_merge_command(
        samtools_path=SAMTOOLS,
        merged_bam="/out/output_sliced_unmapped.bam",
        sliced_bam="/out/output_sliced.bam",
        unmapped_bam="/out/output_unmapped.bam",
        threads=4,
    )

    assert "-u" not in _tokens(command)
