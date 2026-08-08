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
    Both optional flags append after the base command.

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
        "fastp --thread 4 --in1 /data/in_R1.fq.gz --in2 /data/in_R2.fq.gz "
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
    """With both switches off the command ends at the trailing space, as before."""
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

    assert command.endswith("--json /out/output.json ")
    assert "--disable_adapter_trimming" not in command
    assert "--dedup" not in command


@pytest.mark.parametrize(
    ("disable_adapter_trimming", "deduplication", "expected_suffix"),
    [
        (True, True, " --disable_adapter_trimming --dedup"),
        (True, False, " --disable_adapter_trimming"),
        (False, True, " --dedup"),
        (False, False, ""),
    ],
)
def test_every_fastp_flag_combination_is_pinned(disable_adapter_trimming, deduplication, expected_suffix):
    """All four conditional-flag combinations, so neither switch can invert unnoticed."""
    command = build_fastp_command(
        fastp_path=FASTP,
        threads=1,
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


def test_depth_carries_the_reference_so_coverage_does_not_die_on_a_cram():
    """Coverage reads the original CRAM and therefore needs the resolved reference too."""
    command = build_samtools_depth_command(
        samtools_path=SAMTOOLS,
        threads=4,
        region="chr1:1-2",
        bam_file="/o/view.cram",
        coverage_output="/o/d.txt",
        reference_path="/r/g.fa",
    )

    assert "-T /r/g.fa" in command


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
    assert "-f 12" in command
    assert "|" not in command


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


# ---------------------------------------------------------------------------
# D2 - the CRAM process-substitution binary mismatch
# ---------------------------------------------------------------------------


def test_the_cram_filter_runs_the_configured_samtools_in_both_pipeline_stages():
    """
    The writing stage must use the **same** samtools as the reading stage.

    It used to call a bare ``samtools``. Everywhere else in the pipeline the
    binary comes from ``config["tools"]["samtools"]``, and under
    ``mamba run -n <env>`` a bare name resolves against a different PATH - so the
    unmapped reads were extracted by a different build of samtools, or by nothing
    at all while the pipeline reported success.
    """
    configured = "/opt/conda/envs/vntyper/bin/samtools"

    command = build_cram_unmapped_filter_command(
        samtools_path=configured,
        in_bam="/data/sample.cram",
        unmapped_bam="/out/output_unmapped.bam",
        threads=4,
    )

    writer = command.split("|", 1)[1].strip()
    assert writer.startswith(configured), f"the writing stage must call {configured}, got: {writer[:80]!r}"

    bare_calls = [
        part for part in command.split() if part == "samtools"
    ]  # a bare `samtools` token anywhere means the mismatch is back
    assert bare_calls == [], "no bare `samtools` may remain; every invocation comes from config"


def test_the_cram_filter_command_is_pinned():
    """The whole CRAM command, including the deliberate pipefail prefix."""
    command = build_cram_unmapped_filter_command(
        samtools_path=SAMTOOLS,
        in_bam="/data/sample.cram",
        unmapped_bam="/out/output_unmapped.bam",
        threads=4,
    )

    assert command == (
        "set -o pipefail; samtools view -@ 4 -h /data/sample.cram | "
        "samtools view -b -f 12 -@ 4 - -o /out/output_unmapped.bam"
    )


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


def test_the_cram_filter_still_requires_bash():
    """
    Trap 9 survives the pipe rewrite, for a different reason than before.

    ``run_command`` uses ``shell=True`` with ``executable="/bin/bash"``. That was
    justified by the process substitution; with the substitution gone, the reason
    is ``set -o pipefail``, which is not POSIX and is what makes a failure in the
    reading stage fail the whole command instead of being masked by the writer's
    exit status.
    """
    command = build_cram_unmapped_filter_command(
        samtools_path=SAMTOOLS,
        in_bam="/data/sample.cram",
        unmapped_bam="/out/output_unmapped.bam",
        threads=4,
    )

    assert command.startswith("set -o pipefail; "), "pipefail is why this still cannot run under plain sh"


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
    )

    assert command == (
        "set -o pipefail; samtools sort -n -@ 4 /out/output_sliced.bam | "
        "samtools fastq -@ 4 - -1 /out/output_R1.fastq.gz "
        "-2 /out/output_R2.fastq.gz -0 /out/output_other.fastq.gz "
        "-s /out/output_single.fastq.gz"
    )


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
    ],
)
def test_every_multi_stage_pipe_sets_pipefail(name, command):
    """Any command containing a ``|`` must opt into pipefail. No exceptions."""
    assert "|" in command, f"{name} is supposed to be a pipe; this test is checking the wrong thing"
    assert command.startswith(PIPEFAIL_PREFIX), f"{name} pipes without pipefail: an upstream failure exits 0"


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
            ),
            [sliced, r1, r2, other, single],
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
