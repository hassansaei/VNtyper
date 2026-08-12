"""
command_builders.py

Shell command construction for the external tools the BAM/CRAM/FASTQ path drives.

Every command in this module is handed to :func:`vntyper.scripts.utils.run_command`,
which runs it as a **single string** under ``/bin/bash`` with ``shell=True``. That
is deliberate and load-bearing (trap 9): every command here is shell syntax - pipes,
``&&``, output redirection - which no ``shell=False`` argv list can express, and the
``set -o pipefail`` prefix below is a non-POSIX ``set`` option, so pinning
``/bin/bash`` is what makes it available regardless of which shell the host installs
as ``/bin/sh``. The consequence is that quoting can only happen here, at construction
time.

Trap 9 used to be justified by bash **process substitution** specifically, in the CRAM
unmapped-read branch. That substitution is gone - see
:func:`build_cram_unmapped_filter_command` - but the requirement survives it.

Two rules follow, and they pull in opposite directions:

* **Paths, regions and sample names are quoted** with :func:`shlex.quote`. They
  come from the user - an input BAM path, an output directory, a sample name - and
  a space alone is enough to split one argument into two, while a ``;`` or a
  backtick is arbitrary command execution.
* **Tool invocations are not quoted.** ``config["tools"]`` holds *command
  prefixes*, not paths: ``advntr`` is literally ``"mamba run -n envadvntr advntr"``
  (trap 6), and ``--config-path`` replaces the whole config, so an operator may
  legitimately point ``samtools`` at a multi-word invocation. Quoting that would
  make bash look for a binary with spaces in its name. Tool entries are
  operator-controlled configuration, not user input.

The strings produced here are byte-identical to the ones the pipeline produced
before this module existed, apart from six deliberate fixes:

1. ``shlex.quote`` around interpolated paths (a no-op for paths that need no
   quoting, which is every path in the test suite and most real ones).
2. The CRAM unmapped-read extractor now calls the **configured** samtools instead
   of a bare ``samtools`` - under ``mamba run`` those are different binaries.
3. Multi-stage pipes are prefixed with ``set -o pipefail``. Without it a pipeline
   reports only its last stage's exit status, so a ``samtools sort`` that died
   half-way still exited 0 and the next stage genotyped a truncated FASTQ.
4. The CRAM unmapped-read extractor writes through a plain pipe rather than a
   ``tee >(...)`` process substitution, so the shell waits for the writer instead
   of returning while it is still flushing.
5. Non-fast target slicing excludes flag-4 reads (``-F 4``) before merging with
   complete recovery. On the registered b178 single-end fixture this prevents 329
   target reads from being duplicated: the merged count is 4,807 rather than 5,136.
6. Dedup-enabled fastp commands use one worker because fastp 0.23.4's shared
   atomic duplicate table retains whichever representative wins a scheduling race.
   Dedup-disabled commands retain the caller's requested concurrency.

Functions:
    build_fastp_command: FASTQ quality control
    build_samtools_slice_command: region or BED slicing, plus indexing
    build_samtools_index_command: BAM indexing
    build_cram_unmapped_filter_command: CRAM unmapped-read extraction
    build_samtools_merge_command: merge sliced and unmapped BAMs
    build_bam_to_fastq_command: name-sort then convert to paired FASTQ
    build_samtools_depth_command: per-base depth over a region
    build_bwa_align_sort_command: align, convert and coordinate-sort
"""

from __future__ import annotations

import logging
from pathlib import Path

from vntyper.scripts.samtools_command_fragments import customized_index_input, quote_path
from vntyper.scripts.samtools_command_fragments import reference_flag as _reference_flag
from vntyper.scripts.samtools_command_fragments import thread_flag as _thread_flag

logger = logging.getLogger(__name__)

#: Prefix for every command containing a ``|``.
#:
#: bash reports only the exit status of a pipeline's **last** stage. ``samtools
#: sort -n bad.bam | samtools fastq -`` therefore exits 0 when ``sort`` dies, and
#: the caller happily carries on with a FASTQ holding a fraction of the reads.
#: ``pipefail`` makes the pipeline exit non-zero if any stage does.
PIPEFAIL_PREFIX = "set -o pipefail; "


def build_fastp_command(
    *,
    fastp_path: str,
    threads: int,
    fastq_1: str | Path,
    fastq_2: str | Path | None,
    output: str | Path,
    output_name: str,
    compression_level: int,
    qualified_quality_phred: int,
    dup_calc_accuracy: int,
    length_required: int,
    disable_adapter_trimming: bool,
    deduplication: bool,
) -> str:
    """
    Build the fastp quality-control command for one or two FASTQ files.

    Args:
        fastp_path (str): fastp invocation from ``config["tools"]["fastp"]``.
        threads (int): Requested thread count. Deduplication is serialized to one
            fastp worker; without deduplication this value is used unchanged.
        fastq_1 (str | Path): Input R1.
        fastq_2 (str | Path | None): Optional input R2.
        output (str | Path): Output directory.
        output_name (str): Base name for every produced file.
        compression_level (int): gzip compression level for the outputs.
        qualified_quality_phred (int): Per-base quality threshold.
        dup_calc_accuracy (int): fastp duplicate-detection accuracy level.
        length_required (int): Minimum read length to keep.
        disable_adapter_trimming (bool): Append ``--disable_adapter_trimming``.
        deduplication (bool): Append ``--dedup``.

    Returns:
        str: The complete fastp command.

    Note:
        fastp 0.23.4 workers share an atomic Bloom-filter duplicate table. The
        first worker to set a key keeps that read, so parallel scheduling changes
        which representative survives and makes downstream coverage
        nondeterministic. Using one fastp worker when ``deduplication`` is enabled
        makes read membership reproducible; it does not change the thread count
        requested from BWA or samtools.

        The base command ends with a trailing space and each optional flag is
        appended with a leading one, so both flags on produces a double space
        before ``--disable_adapter_trimming``. That is preserved from the
        pre-extraction code so the refactor is byte-identical; bash ignores it.
    """
    out_1 = f"{output}/{output_name}_R1.fastq.gz"
    html = f"{output}/{output_name}.html"
    json_report = f"{output}/{output_name}.json"

    input_fragment = f"--in1 {quote_path(fastq_1)} "
    output_fragment = f"--out1 {quote_path(out_1)} "
    if fastq_2 is not None:
        out_2 = f"{output}/{output_name}_R2.fastq.gz"
        input_fragment += f"--in2 {quote_path(fastq_2)} "
        output_fragment += f"--out2 {quote_path(out_2)} "

    fastp_threads = 1 if deduplication else threads
    command = (
        f"{fastp_path} --thread {quote_path(fastp_threads)} "
        f"{input_fragment}"
        f"{output_fragment}"
        f"--compression {quote_path(compression_level)} "
        f"--qualified_quality_phred {quote_path(qualified_quality_phred)} "
        f"--dup_calc_accuracy {quote_path(dup_calc_accuracy)} "
        f"--length_required {quote_path(length_required)} "
        f"--html {quote_path(html)} "
        f"--json {quote_path(json_report)} "
    )

    if disable_adapter_trimming:
        command += " --disable_adapter_trimming"
    if deduplication:
        command += " --dedup"

    return command


def build_samtools_slice_command(
    *,
    samtools_path: str,
    in_bam: str | Path,
    output_bam: str | Path,
    region: str | None = None,
    bed_file: str | Path | None = None,
    reference_path: str | Path | None = None,
    index_path: str | Path | None = None,
    threads: int = 1,
    index_output: bool = True,
    exclude_unmapped: bool = False,
    uncompressed: bool = False,
) -> str:
    """
    Build the region-slicing command, followed by indexing the slice.

    Exactly one of ``region`` and ``bed_file`` is used; ``bed_file`` wins if both
    are given, matching the pre-extraction branch order.

    Args:
        samtools_path (str): samtools invocation from config.
        in_bam (str | Path): Input BAM or CRAM.
        output_bam (str | Path): Where the slice is written.
        region (str | None): Region string such as ``chr1:155158000-155163000``.
        bed_file (str | Path | None): BED file passed as ``-L``.
        reference_path (str | Path | None): Reference FASTA for CRAM decoding.
        index_path (str | Path | None): Exact custom index passed with ``-X``.
        threads (int): Thread count for view and index.
        index_output (bool): Whether to append indexing of the resulting slice.
        exclude_unmapped (bool): Exclude flag-4 reads recovered separately for a
            subsequent disjoint merge. Fast-mode slices leave this disabled.
        uncompressed (bool): Emit ``-u``, BGZF level 0. Still a valid, indexable
            BAM, just not deflated - worth taking only where the file is re-read
            within milliseconds and then replaced or deleted. Defaults to False,
            because a slice that **survives** the run is archived and shipped to
            users, and shipping it at level 0 would roughly triple its size.

    Returns:
        str: ``samtools view ...``, optionally followed by ``&& samtools index``.

    Raises:
        ValueError: If neither ``region`` nor ``bed_file`` is given. Slicing with
            no target would read the whole file, which is never what any caller
            wants and would silently cost hours.
    """
    if bed_file is not None:
        target = f"-L {quote_path(bed_file)}"
    elif region is not None:
        target = quote_path(region)
    else:
        message = "build_samtools_slice_command needs either a region or a bed_file"
        logger.error(message)
        raise ValueError(message)

    exclude_flag = "-F 4 " if exclude_unmapped else ""
    level_flag = "-u " if uncompressed else ""
    indexed_input = customized_index_input(in_bam, index_path)
    command = (
        f"{samtools_path} view -P -b {level_flag}{exclude_flag}{_thread_flag(threads)}"
        f"{_reference_flag(reference_path)}{indexed_input} {target} -o {quote_path(output_bam)}"
    )
    if not index_output:
        return command
    return f"{command} && {samtools_path} index {_thread_flag(threads)}{quote_path(output_bam)}"


def build_samtools_index_command(
    *, samtools_path: str, bam_file: str | Path, output_bai: str | Path | None = None, threads: int = 1
) -> str:
    """
    Build the ``samtools index`` command for a BAM file.

    Args:
        samtools_path (str): samtools invocation from config.
        bam_file (str | Path): The BAM to index.
        output_bai (str | Path, optional): Where to write the index. Defaults to
            None, which lets samtools write ``<bam_file>.bai`` beside the input --
            correct for a BAM the pipeline itself produced in its output
            directory, wrong for the user's input alignment, whose directory is
            routinely read-only (#162, #210). ``-o`` requires samtools >= 1.15;
            ``conda/environment_vntyper.yml`` pins 1.20.
        threads (int): Thread count for indexing.

    Returns:
        str: The complete index command.
    """
    if output_bai is None:
        return f"{samtools_path} index {_thread_flag(threads)}{quote_path(bam_file)}"
    return f"{samtools_path} index {_thread_flag(threads)}-o {quote_path(output_bai)} {quote_path(bam_file)}"


def build_cram_unmapped_filter_command(
    *,
    samtools_path: str,
    in_bam: str | Path,
    unmapped_bam: str | Path,
    threads: int,
    reference_path: str | Path | None = None,
    uncompressed: bool = False,
) -> str:
    """
    Build the whole-file unmapped-read extraction command.

    Stream plans for BAM and CRAM decode the whole file through samtools and pick
    out flag 4 (read unmapped). This also retains unpaired reads, which do not
    carry the mate-unmapped bit.

    Args:
        samtools_path (str): samtools invocation from config. This is used for
            **both** stages of the pipeline - the writing stage used to be a bare
            ``samtools``, which under ``mamba run`` resolves against a different
            PATH.
        in_bam (str | Path): The input BAM or CRAM.
        unmapped_bam (str | Path): Where the unmapped reads are written.
        threads (int): Thread count for both samtools invocations.
        reference_path (str | Path | None): Reference FASTA for CRAM decoding.
        uncompressed (bool): Emit ``-u``, BGZF level 0. The unmapped BAM is merged
            milliseconds later and then deleted, so deflating it is pure cost -
            measured 0.869 s to 0.081 s on a 963,549-read BAM.

    Returns:
        str: The complete command, prefixed with ``set -o pipefail``.

    Note:
        The writing samtools is a **pipeline stage**, not a ``tee >(...)`` process
        substitution, and that is the whole point of this builder's shape.

        bash does not wait for a process substitution and a trailing ``wait`` does
        not change that: the shell returns as soon as ``tee`` exits, while the
        substituted samtools is still flushing ``unmapped_bam``.
        :func:`~vntyper.scripts.fastq_bam_processing.process_bam_to_fastq` runs
        ``samtools merge`` against that file on the very next line, so it merged a
        BAM that was still being written - measured on a 600k-read CRAM, 199 797 of
        200 000 unmapped reads were present at the instant the shell returned, and
        ``samtools merge`` accepted the short file with only a ``W::bam_hdr_read``
        warning. Those are exactly the reads this stage exists to recover for
        Kestrel, so the pipeline reported success on an under-called sample.
        ``pipefail`` could not catch it either: a substitution that consumes its
        whole input and *then* fails is invisible, because ``tee`` never sees EPIPE
        and exits 0.

        A plain pipe fixes both halves at once. ``tee``'s own stdout went to
        ``/dev/null``, so the substitution was the pipeline's only consumer and
        ``tee`` was pure overhead; making the writer a real stage means bash waits
        for it before the shell returns, and ``pipefail`` becomes meaningful because
        the writer's exit status is now part of the pipeline's. The bytes reaching
        the writer are unchanged - SAM text on stdin either way.
    """
    level_flag = "-u " if uncompressed else ""
    return (
        f"{PIPEFAIL_PREFIX}"
        f"{samtools_path} view {_reference_flag(reference_path)}{_thread_flag(threads)}-h {quote_path(in_bam)} | "
        f"{samtools_path} view -b -f 4 {level_flag}{_thread_flag(threads)}- -o {quote_path(unmapped_bam)}"
    )


def build_cram_unmapped_indexed_command(
    *,
    samtools_path: str,
    in_bam: str | Path,
    unmapped_bam: str | Path,
    threads: int,
    reference_path: str | Path | None = None,
    index_path: str | Path | None = None,
    uncompressed: bool = False,
) -> str:
    """Build the indexed alignment command for unplaced unmapped reads.

    Args:
        samtools_path: samtools invocation from config.
        in_bam: Indexed BAM or CRAM alignment view.
        unmapped_bam: Destination BAM for read-unmapped records.
        threads: Thread count for samtools.
        reference_path: Reference FASTA for CRAM decoding.
        index_path: Exact custom index passed with ``-X``.
        uncompressed: Emit ``-u``, BGZF level 0. This writes the same throwaway file
            as the stream builder, so it takes the same decision - otherwise the
            saving would depend on which scan preflight happened to prove.

    Returns:
        A single ``samtools view`` command that requests literal ``'*'`` reads
        whose read-unmapped bit is set, including unpaired reads.
    """
    level_flag = "-u " if uncompressed else ""
    indexed_input = customized_index_input(in_bam, index_path)
    return (
        f"{samtools_path} view -b -f 4 {level_flag}{_thread_flag(threads)}{_reference_flag(reference_path)}"
        f"{indexed_input} {quote_path('*')} -o {quote_path(unmapped_bam)}"
    )


def build_cram_reference_probe_command(
    *,
    samtools_path: str,
    in_bam: str | Path,
    region: str | None = None,
    bed_file: str | Path | None = None,
    reference_path: str | Path | None = None,
    threads: int = 1,
) -> str:
    """Build a throwaway slice that proves a reference decodes the requested target.

    Args:
        samtools_path: samtools invocation from config.
        in_bam: Indexed BAM or CRAM alignment view.
        region: Region string to probe.
        bed_file: BED target to probe. Takes precedence over ``region``.
        reference_path: Candidate reference FASTA.
        threads: Thread count for samtools.

    Returns:
        A CRAM-decoding ``samtools view`` command which writes to ``/dev/null``.

    Raises:
        ValueError: If neither a region nor a BED target was supplied.
    """
    if bed_file is not None:
        target = f"-L {quote_path(bed_file)}"
    elif region is not None:
        target = quote_path(region)
    else:
        message = "build_cram_reference_probe_command needs either a region or a bed_file"
        logger.error(message)
        raise ValueError(message)

    return (
        f"{samtools_path} view -P {_thread_flag(threads)}{_reference_flag(reference_path)}"
        f"{quote_path(in_bam)} {target} -o /dev/null"
    )


def build_cram_stream_reference_probe_command(
    *,
    samtools_path: str,
    in_bam: str | Path,
    reference_path: str | Path | None = None,
    threads: int = 1,
) -> str:
    """Build the whole-file decode proof required by the stream consumer.

    Args:
        samtools_path: Samtools invocation from config.
        in_bam: CRAM alignment view to decode completely.
        reference_path: Candidate reference FASTA, or ``None`` for htslib resolution.
        threads: Thread count for samtools.

    Returns:
        A whole-file SAM decode whose output is discarded.
    """
    return (
        f"{samtools_path} view {_thread_flag(threads)}{_reference_flag(reference_path)}"
        f"-h {quote_path(in_bam)} -o /dev/null"
    )


def build_samtools_idxstats_command(*, samtools_path: str, in_bam: str | Path, threads: int = 1) -> str:
    """Build ``samtools idxstats`` without a reference argument.

    Args:
        samtools_path: samtools invocation from config.
        in_bam: Indexed BAM or CRAM alignment view.
        threads: Thread count for samtools.

    Returns:
        The complete ``idxstats`` command.
    """
    return f"{samtools_path} idxstats {_thread_flag(threads)}{quote_path(in_bam)}"


def build_samtools_unmapped_indexed_count_command(*, samtools_path: str, in_bam: str | Path, threads: int = 1) -> str:
    """Count records visible to the exact indexed unmapped-read consumer.

    Args:
        samtools_path: samtools invocation from config.
        in_bam: Indexed BAM or CRAM alignment view.
        threads: Thread count for samtools.

    Returns:
        A literal-``'*'`` flag-4 count command without a reference argument.
    """
    return f"{samtools_path} view -c -f 4 {_thread_flag(threads)}{quote_path(in_bam)} {quote_path('*')}"


def build_samtools_merge_command(
    *,
    samtools_path: str,
    merged_bam: str | Path,
    sliced_bam: str | Path,
    unmapped_bam: str | Path,
    threads: int,
) -> str:
    """
    Build the command merging the region slice with the extracted unmapped reads.

    Args:
        samtools_path (str): samtools invocation from config.
        merged_bam (str | Path): Output path for the merged BAM.
        sliced_bam (str | Path): The region slice.
        unmapped_bam (str | Path): The extracted unmapped reads.
        threads (int): Thread count.

    Returns:
        str: The complete merge command. ``-f`` overwrites an existing output.
    """
    return (
        f"{samtools_path} merge -f -@ {quote_path(threads)} "
        f"{quote_path(merged_bam)} {quote_path(sliced_bam)} {quote_path(unmapped_bam)}"
    )


def build_bam_to_fastq_command(
    *,
    samtools_path: str,
    in_bam: str | Path,
    threads: int,
    fastq_r1: str | Path,
    fastq_r2: str | Path,
    fastq_other: str | Path,
    fastq_single: str | Path,
) -> str:
    """
    Build the name-sort-then-convert pipeline that turns a BAM back into FASTQ.

    Args:
        samtools_path (str): samtools invocation from config.
        in_bam (str | Path): The BAM to convert.
        threads (int): Thread count for both stages.
        fastq_r1 (str | Path): Output R1.
        fastq_r2 (str | Path): Output R2.
        fastq_other (str | Path): Output for reads that are neither R1 nor R2.
        fastq_single (str | Path): Output for singleton reads.

    Returns:
        str: The complete pipeline, prefixed with ``set -o pipefail``.

    Note:
        The ``pipefail`` prefix is the point of this builder. ``samtools fastq``
        exits 0 on a short input stream, so before it a ``samtools sort`` that ran
        out of disk or hit a truncated BAM produced a **partial FASTQ that the
        pipeline then genotyped without a warning**. With ``pipefail`` the stage
        fails loudly instead.
    """
    return (
        f"{PIPEFAIL_PREFIX}"
        f"{samtools_path} sort -n -@ {quote_path(threads)} {quote_path(in_bam)} | "
        f"{samtools_path} fastq -@ {quote_path(threads)} - -1 {quote_path(fastq_r1)} "
        f"-2 {quote_path(fastq_r2)} -0 {quote_path(fastq_other)} "
        f"-s {quote_path(fastq_single)}"
    )


def build_samtools_depth_command(
    *,
    samtools_path: str,
    threads: int,
    region: str,
    bam_file: str | Path,
    coverage_output: str | Path,
    reference_path: str | Path | None = None,
    index_path: str | Path | None = None,
) -> str:
    """
    Build the per-base depth command for a region, redirected to a file.

    Args:
        samtools_path (str): samtools invocation from config.
        threads (int): Thread count.
        region (str): Region string such as ``chr1:155160500-155162000``.
        bam_file (str | Path): The BAM to measure.
        coverage_output (str | Path): Where the depth table is written.
        reference_path (str | Path | None): Reference FASTA for CRAM decoding.
        index_path (str | Path | None): Exact custom index passed with ``-X``.

    Returns:
        str: The complete command.

    Note:
        This is one command with a redirect, not a pipe, so its exit status is
        already samtools' own and ``pipefail`` would add nothing.

        ``-a`` makes samtools emit a row for every position in the region, zero-depth
        ones included. Without it the depth table is the covered positions only, and
        every statistic computed from it is over a base set that varies per sample
        (#171).
    """
    reference_flag = "" if reference_path is None else f"--reference {quote_path(reference_path)} "
    indexed_input = customized_index_input(bam_file, index_path)
    return (
        f"{samtools_path} depth -a {_thread_flag(threads)}{reference_flag}-r {quote_path(region)} "
        f"{indexed_input} > {quote_path(coverage_output)}"
    )


def build_bwa_align_sort_command(
    *,
    bwa_path: str,
    samtools_path: str,
    threads: int,
    reference: str | Path,
    fastq1: str | Path,
    fastq2: str | Path | None,
    sorted_bam: str | Path,
) -> str:
    """
    Build the three-stage align, convert and coordinate-sort pipeline.

    Args:
        bwa_path (str): bwa invocation from config.
        samtools_path (str): samtools invocation from config.
        threads (int): Thread count for all three stages.
        reference (str | Path): Reference FASTA, already BWA-indexed.
        fastq1 (str | Path): Input R1.
        fastq2 (str | Path | None): Optional input R2.
        sorted_bam (str | Path): Output coordinate-sorted BAM.

    Returns:
        str: The complete pipeline, prefixed with ``set -o pipefail``.

    Note:
        Three stages means two chances for a failure to be masked. A ``bwa mem``
        that aborts part-way still lets ``samtools sort`` write a valid,
        **incomplete** BAM and exit 0.
    """
    fastq_inputs = quote_path(fastq1)
    if fastq2 is not None:
        fastq_inputs += f" {quote_path(fastq2)}"

    return (
        f"{PIPEFAIL_PREFIX}"
        f"{bwa_path} mem -t {quote_path(threads)} {quote_path(reference)} "
        f"{fastq_inputs} | "
        f"{samtools_path} view -@ {quote_path(threads)} -b | "
        f"{samtools_path} sort -@ {quote_path(threads)} -o {quote_path(sorted_bam)}"
    )
