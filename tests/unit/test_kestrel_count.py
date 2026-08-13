"""The standalone KAnalyze count command that produces Kestrel's IKC input (#262).

Kestrel 1.0.1 builds a KAnalyze ``CountModule`` and configures almost nothing on it.
``KestrelRunnerBase.getCountModule()`` calls exactly ``configure(null)``,
``addLibraryURL`` per ``--lib``/``--liburl``, ``setKSize``, ``setTempDirName``,
``addPostCountFilterDefinition("kmercount:" + minKmerCount)`` when that is positive, and
``setFreeSegment``; ``IkcCountMap``'s constructor adds ``setOutputFormat("ikc")``,
``setMinimizerSize`` and ``setMinimizerMask``. It never calls ``setKmerThreadCount``,
``setSplitThreadCount`` or ``setThreads``, all three of which exist on ``CountModule``,
so counting runs at the compile-time defaults of one k-mer and one split thread.

Handing Kestrel a pre-built IKC is a supported input: ``IkcCountMap.preModuleRun``
adopts the sample's single file when its format is ``ikc`` or when the format is
``auto`` and the name matches ``.*\\.ikc``, sets ``rmLastTemp = false`` -- so Kestrel
never deletes an IKC it was handed, whatever ``--rmikc`` says -- and returns false so
the count module never runs at all. That last part is why this module has to reproduce
Kestrel's counting parameters exactly: nothing downstream will re-derive them.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from vntyper.scripts.kestrel_count import allocate_count_threads, construct_kanalyze_count_command

pytestmark = pytest.mark.unit


def _command(**overrides) -> str:
    """Build a count command from a complete default argument set.

    Args:
        **overrides: Replacements for any default argument.

    Returns:
        str: The rendered command.
    """
    kwargs: dict[str, Any] = {
        "kmer_size": 20,
        "kanalyze_path": "k/kanalyze.jar",
        "output_dir": "/out",
        "ikc_out": "/out/kmer_20/kestrel_kmers.ikc",
        "fastq_files": ["/a_R1.fastq.gz", "/a_R2.fastq.gz"],
        "java_path": "java",
        "java_memory": "12g",
        "java_opts": "",
        "threads": 8,
    }
    kwargs.update(overrides)
    return construct_kanalyze_count_command(**kwargs)


# ---------------------------------------------------------------------------
# The rendered command
# ---------------------------------------------------------------------------


def test_the_count_command_mirrors_the_defaults_kestrel_would_have_used():
    """Byte for byte, because every token corresponds to a Kestrel default.

    ``-c kmercount:5`` mirrors ``--mincount`` 5, ``--minsize 15`` mirrors ``--minsize``
    15 (which ``IkcCountMap`` requires to be non-zero), and the charset, sequence
    filter and reverse-complement settings are left at KAnalyze's own defaults because
    Kestrel does not change them either.
    """
    assert _command() == (
        "java -Xmx12g -jar k/kanalyze.jar count -k 20 -c kmercount:5 --minsize 15 "
        "-m ikc -d 4 -l 2 -t 2 --temploc /out -o /out/kmer_20/kestrel_kmers.ikc "
        "/a_R1.fastq.gz /a_R2.fastq.gz"
    )


def test_the_count_command_carries_no_option_kestrel_leaves_at_its_default():
    """Every difference from Kestrel's own CountModule is a chance to change genotypes.

    Kestrel sets no charset, no sequence filter and no reverse-complement mode, so
    neither does this: KAnalyze's own defaults (UTF-8, ``--noseqfilter``,
    ``--noreverse``) are already what Kestrel would have got. ``--minmask`` is likewise
    omitted because Kestrel's default is 0 and so is KAnalyze's.
    """
    command = _command()

    for absent in ("--charset", "--seqfilter", "--quality", "-r ", "--minmask", "--memcount", "--lib"):
        assert absent not in command, f"{absent} diverges from Kestrel's own count configuration"


def test_the_kmer_size_is_the_one_the_call_step_will_use():
    """A k-size mismatch is the one divergence ``IkcReader`` itself rejects."""
    assert " -k 31 " in _command(kmer_size=31)


def test_java_options_are_placed_between_the_heap_size_and_the_jar():
    """Anything after ``-jar`` is an argument to the program, not to the JVM."""
    command = _command(java_opts="-XX:+UseSerialGC")

    assert command.startswith("java -Xmx12g -XX:+UseSerialGC -jar k/kanalyze.jar count ")


def test_no_empty_option_slot_is_left_when_there_are_no_java_options():
    """The shipped default is empty, and a stray double space is a rendering bug."""
    assert "  " not in _command(java_opts="")


def test_the_java_invocation_is_not_quoted_because_it_may_be_a_command_prefix():
    """``config["tools"]`` holds command prefixes; quoting collapses them to one token."""
    command = _command(java_path="mamba run -n vntyper java")

    assert command.startswith("mamba run -n vntyper java -Xmx12g -jar ")


def test_paths_with_spaces_and_metacharacters_are_quoted():
    """``run_command`` executes this as one string under bash, so quoting happens here."""
    command = _command(
        output_dir="/o dir",
        ikc_out="/o dir/x.ikc",
        kanalyze_path="/jars/kanalyze one.jar",
        fastq_files=["/a b.fastq.gz", "/c;touch pwned.fastq.gz"],
    )

    assert "'/o dir'" in command
    assert "'/o dir/x.ikc'" in command
    assert "'/jars/kanalyze one.jar'" in command
    assert "'/a b.fastq.gz'" in command
    assert "'/c;touch pwned.fastq.gz'" in command


def test_the_fastq_operands_keep_the_order_they_were_routed_in():
    """All files belong to one sample, and the sample is the unit KAnalyze counts."""
    command = _command(fastq_files=["/z.fastq.gz", "/a.fastq.gz", "/m.fastq.gz"])

    assert command.endswith("/z.fastq.gz /a.fastq.gz /m.fastq.gz")


def test_a_path_object_is_accepted_wherever_a_string_is():
    """The caller has Paths; requiring str() at every call site invites a missed one."""
    command = _command(ikc_out=Path("/out/x.ikc"), output_dir=Path("/out"), fastq_files=[Path("/a.fastq.gz")])

    assert "-o /out/x.ikc" in command
    assert command.endswith("/a.fastq.gz")


# ---------------------------------------------------------------------------
# Rejected inputs
# ---------------------------------------------------------------------------


def test_an_empty_fastq_sequence_is_rejected():
    """Counting nothing produces an empty IKC and a confident negative genotype."""
    with pytest.raises(ValueError, match="FASTQ input files are missing"):
        _command(fastq_files=[])


@pytest.mark.parametrize("bad", ["/a.fastq.gz", Path("/a.fastq.gz")])
def test_a_scalar_instead_of_a_sequence_is_rejected(bad):
    """A bare string is iterable, so this would otherwise count one file per character."""
    with pytest.raises(ValueError, match="non-scalar sequence"):
        _command(fastq_files=bad)


@pytest.mark.parametrize("bad", [[""], [None], [3]])
def test_an_invalid_fastq_entry_is_rejected(bad):
    with pytest.raises(ValueError, match="FASTQ input files are missing or invalid"):
        _command(fastq_files=bad)


def test_duplicate_fastq_paths_are_rejected():
    """KAnalyze would count the same reads twice, doubling every k-mer count."""
    with pytest.raises(ValueError, match="duplicate paths"):
        _command(fastq_files=["/a.fastq.gz", "/a.fastq.gz"])


@pytest.mark.parametrize("bad", [0, -1, True, "20", 20.0])
def test_an_invalid_kmer_size_is_rejected(bad):
    with pytest.raises(ValueError, match="k-mer size"):
        _command(kmer_size=bad)


@pytest.mark.parametrize("bad", [0, -1, True, "8", 8.0])
def test_an_invalid_thread_count_is_rejected(bad):
    with pytest.raises(ValueError, match="threads"):
        _command(threads=bad)


# ---------------------------------------------------------------------------
# The thread budget
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("threads", "expected"),
    [
        (1, (1, 1, 1)),
        (2, (1, 1, 1)),
        (3, (1, 1, 1)),
        (4, (2, 1, 1)),
        (6, (4, 1, 1)),
        (8, (4, 2, 2)),
        (12, (6, 3, 3)),
        (16, (8, 4, 4)),
        (32, (16, 8, 8)),
    ],
)
def test_the_budget_is_allocated_across_the_three_stages(threads: int, expected: tuple[int, int, int]):
    """``--threads`` is a total concurrency budget, not a per-stage multiplier.

    ``-d`` (k-mer), ``-l`` (split) and ``-t`` (sort) are independent concurrent
    KAnalyze stages -- ``kmerThreadCount`` sizes an array of ``KmerComponent``,
    ``splitThreadCount`` an array of ``CountSplitComponent``, and ``threads`` the
    ``SortSynchronizer``'s worker pool. Passing the caller's ``--threads`` to all three
    would start roughly 2.5x the requested concurrency, and this pipeline already reads
    ``-@ N`` as "N *additional* threads", so a budget that silently multiplied would be
    doubly surprising.
    """
    assert allocate_count_threads(threads) == expected


@pytest.mark.parametrize("threads", [1, 2, 3, 4, 5, 8, 13, 16, 31, 32, 64])
def test_no_stage_is_ever_allocated_zero_workers(threads: int):
    """KAnalyze rejects a zero thread count, and a stage with none would deadlock."""
    assert all(count >= 1 for count in allocate_count_threads(threads))


@pytest.mark.parametrize("threads", [3, 4, 5, 8, 13, 16, 31, 32, 64])
def test_the_allocation_never_exceeds_the_budget_above_two_threads(threads: int):
    """Below three, flooring every stage at one necessarily overshoots; say so.

    At one or two threads the three floors already cost three workers. That is
    deliberate and unavoidable -- a stage with zero workers cannot run -- and it is
    still far below the seven-ish workers stock Kestrel starts regardless of --threads.
    """
    assert sum(allocate_count_threads(threads)) <= threads


@pytest.mark.parametrize("threads", [4, 8, 16, 32])
def test_the_kmer_stage_gets_the_largest_share(threads: int):
    """It is the measured bottleneck, and the plan had this backwards.

    Raising ``-d`` from 1 to 4 while leaving ``-t`` at its default of 4 took the count
    step from 31.72 s to 16.62 s, so the k-mer stage -- not sort -- is what the budget
    should buy. An allocation that gave sort half the budget would hand ``--threads 4``
    fewer sort workers than stock KAnalyze uses by default, for no gain elsewhere.
    """
    kmer, split, sort = allocate_count_threads(threads)

    assert kmer >= split
    assert kmer >= sort


def test_a_thread_budget_below_one_is_rejected_rather_than_silently_floored():
    """A caller passing 0 has a bug; floooring it hides that from them."""
    with pytest.raises(ValueError, match="threads"):
        allocate_count_threads(0)
