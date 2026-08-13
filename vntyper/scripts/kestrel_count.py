"""Build the command line for the vendored KAnalyze k-mer counter.

Kestrel 1.0.1 constructs a KAnalyze ``CountModule`` and configures almost nothing on
it. ``KestrelRunnerBase.getCountModule()`` calls ``configure(null)``, ``addLibraryURL``
once per ``--lib``/``--liburl``, ``setKSize``, ``setTempDirName``,
``addPostCountFilterDefinition("kmercount:" + minKmerCount)`` when that value is
positive, and ``setFreeSegment``; ``IkcCountMap``'s constructor then adds
``setOutputFormat("ikc")``, ``setMinimizerSize`` and ``setMinimizerMask``. It never
calls ``setKmerThreadCount``, ``setSplitThreadCount`` or ``setThreads``, all three of
which exist on ``CountModule``, so counting runs at the compile-time defaults of one
k-mer thread and one split thread -- while counting is roughly 88% of Kestrel's work on
an unmapped-dominated input.

Running the bundled ``kanalyze.jar`` as its own step and handing Kestrel the resulting
indexed k-mer count (IKC) file lets that step use the machine. Kestrel accepts a
supplied IKC: ``IkcCountMap.preModuleRun`` adopts the sample's single input file when
its format is ``ikc`` or when the format is ``auto`` and the file name matches
``.*\\.ikc``, sets ``rmLastTemp = false`` -- so Kestrel never deletes an IKC it was
handed, whatever ``--rmikc`` says -- and returns false so the count module does not run.

That last point is why every parameter here has to mirror the value Kestrel would have
used: once the IKC is supplied, nothing downstream re-derives any of them, and
``IkcReader`` validates only the k-size. A same-k IKC built with a different minimum
count or minimizer size is silently accepted and silently different. This module's job
is to make that impossible for the shipped configuration, and
``construct_kestrel_command``'s allowlist is what stops an operator reintroducing it
through ``additional_settings``.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from vntyper.scripts.command_builders import quote_path

logger = logging.getLogger(__name__)

#: Kestrel's ``--mincount`` default, applied as KAnalyze's ``kmercount`` post-count
#: filter. ``getCountModule`` adds it only when the value is positive, and
#: ``kestrel_config.json`` does not override it.
KESTREL_DEFAULT_MIN_COUNT = 5

#: Kestrel's ``--minsize`` default, the minimizer size ``IkcCountMap`` passes to
#: ``setMinimizerSize``. It must be non-zero: the constructor raises "K-mer utility was
#: not configured with a minimizer size" otherwise, and the IKC format needs the
#: minimizer grouping to be searchable.
KESTREL_DEFAULT_MIN_SIZE = 15


def _positive_int(value: object, name: str) -> int:
    """Validate a positive integer argument.

    Args:
        value: The supplied value.
        name: The argument name, for the message.

    Returns:
        int: The validated value.

    Raises:
        ValueError: If the value is not a positive, non-boolean integer.
    """
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        msg = f"{name} must be a positive integer."
        logger.error(msg)
        raise ValueError(msg)
    return value


def allocate_count_threads(threads: int) -> tuple[int, int, int]:
    """Split one thread budget across KAnalyze's three concurrent count stages.

    ``-d`` (k-mer), ``-l`` (split) and ``-t`` (sort) are independent concurrent stages,
    not three spellings of one worker count: ``kmerThreadCount`` sizes an array of
    ``KmerComponent``, ``splitThreadCount`` an array of ``CountSplitComponent``, and
    ``threads`` the ``SortSynchronizer``'s worker pool. Passing the caller's
    ``--threads`` to all three would start roughly 2.5x the requested concurrency -- on
    a four-core host ``--threads 4`` would run about ten count workers -- and this
    pipeline already reads ``-@ N`` as "N *additional* threads", so a budget that
    silently multiplied would be doubly surprising.

    The k-mer stage takes the largest share because it is the measured bottleneck:
    raising ``-d`` from 1 to 4 while leaving ``-t`` at KAnalyze's own default of 4 took
    the count step from 31.72 s to 16.62 s. Every measured configuration also used
    ``-d : -l`` of 2 : 1, which is preserved here.

    Args:
        threads: The run's total thread budget.

    Returns:
        tuple[int, int, int]: ``(-d, -l, -t)``. Every stage is floored at one worker,
        so the three sum to at most ``threads`` for ``threads >= 3`` and to exactly
        three below that -- a stage with no workers cannot run, and three is still far
        below the seven-odd workers stock Kestrel starts whatever ``--threads`` says.

    Raises:
        ValueError: If ``threads`` is not a positive integer. A caller passing zero has
            a bug, and flooring it silently would hide that.
    """
    budget = _positive_int(threads, "threads")
    # Split and sort take a quarter each and the k-mer stage takes the remainder, so
    # the three sum to exactly the budget above the floors rather than approximately.
    split_threads = max(1, budget // 4)
    sort_threads = max(1, budget // 4)
    kmer_threads = max(1, budget - split_threads - sort_threads)
    return kmer_threads, split_threads, sort_threads


def construct_kanalyze_count_command(
    *,
    kmer_size: int,
    kanalyze_path: str | Path,
    output_dir: str | Path,
    ikc_out: str | Path,
    fastq_files: Sequence[str | Path],
    java_path: str,
    java_memory: str,
    java_opts: str,
    threads: int,
    min_count: int = KESTREL_DEFAULT_MIN_COUNT,
    min_size: int = KESTREL_DEFAULT_MIN_SIZE,
) -> str:
    """Build the KAnalyze count command that produces Kestrel's IKC input.

    Args:
        kmer_size: K-mer size, which must match the Kestrel invocation that reads the
            IKC. ``IkcReader`` validates this one and rejects a mismatch.
        kanalyze_path: Path to ``kanalyze.jar``.
        output_dir: Directory KAnalyze offloads its sorted segments into.
        ikc_out: Destination IKC file.
        fastq_files: Non-empty ordered FASTQ paths, all belonging to one sample.
        java_path: Java invocation from config. Left unquoted because
            ``config["tools"]`` holds command *prefixes*, which quoting would collapse
            into a single token bash looks for as one binary.
        java_memory: Heap size such as ``"12g"``.
        java_opts: Extra JVM options, rendered between the heap size and ``-jar``.
            Empty keeps the G1 default, which the count step wants because it is
            genuinely parallel.
        threads: The run's total thread budget, allocated by
            :func:`allocate_count_threads`.
        min_count: Minimum k-mer count, mirroring Kestrel's ``--mincount``.
        min_size: Minimizer size, mirroring Kestrel's ``--minsize``.

    Returns:
        str: The complete command line.

    Raises:
        ValueError: If the FASTQ operands are missing, scalar, invalid or duplicated,
            or if a numeric argument is not a positive integer.
    """
    _positive_int(kmer_size, "k-mer size")
    _positive_int(min_count, "min_count")
    _positive_int(min_size, "min_size")
    kmer_threads, split_threads, sort_threads = allocate_count_threads(threads)

    if isinstance(fastq_files, (str, Path)) or not isinstance(fastq_files, Sequence):
        msg = "fastq_files must be a non-scalar sequence of FASTQ paths."
        logger.error(msg)
        raise ValueError(msg)
    normalized: list[str] = []
    for fastq in fastq_files:
        if not isinstance(fastq, (str, Path)) or not str(fastq):
            msg = "FASTQ input files are missing or invalid."
            logger.error(msg)
            raise ValueError(msg)
        normalized.append(str(fastq))
    if not normalized:
        msg = "FASTQ input files are missing or invalid."
        logger.error(msg)
        raise ValueError(msg)
    if len(set(normalized)) != len(normalized):
        # KAnalyze would count the same reads twice, doubling every k-mer count and
        # moving every variant across the calibrated Depth_Score thresholds.
        msg = "FASTQ input files contain duplicate paths."
        logger.error(msg)
        raise ValueError(msg)

    options = f" {java_opts}" if java_opts else ""
    inputs = " ".join(quote_path(path) for path in normalized)
    return (
        f"{java_path} -Xmx{java_memory}{options} -jar {quote_path(kanalyze_path)} count "
        f"-k {kmer_size} -c kmercount:{min_count} --minsize {min_size} "
        f"-m ikc -d {kmer_threads} -l {split_threads} -t {sort_threads} "
        f"--temploc {quote_path(output_dir)} -o {quote_path(ikc_out)} {inputs}"
    )
