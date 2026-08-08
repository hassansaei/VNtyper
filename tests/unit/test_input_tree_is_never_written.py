"""No code in the ``vntyper`` package writes into the input directory (#162, #201, #210).

SPECIFICATION: the input directory holds patient-derived data and is routinely
mounted read-only to protect it, which is what #162 reports. Two production
writes in the ``vntyper`` package derived an output path from an input path:

1. ``validate_bam_file`` wrote ``<input>.quickcheck.log`` beside the alignment
   (#201). On the web service that is the per-job directory on the shared upload
   volume, so ``run_vntyper_job``'s cleanup found a non-empty directory and its
   ``os.rmdir`` never ran -- one leaked directory and inode per job, forever, and
   a "still holds files" warning that fired on 100% of jobs and therefore could
   no longer report a genuinely unexpected leftover.
2. ``process_bam_to_fastq`` ran ``samtools index <in_bam>`` for non-fast BAM runs
   whenever ``<in_bam>.bai`` was missing -- including when a usable ``<stem>.bai``
   was already there, which is a name the upload endpoint and the worker both
   deliberately accept (#210).

Both are crashes, not annoyances: ``run_command`` opens its log file before it
runs anything, so on a read-only mount write (1) is an unhandled
``PermissionError`` raised before quickcheck even executes.

The rule deciding *which* index counts as already present lives in
``vntyper/scripts/alignment_index.py`` and is pinned by
``tests/unit/test_alignment_index.py``. What is tested here is the invariant it
serves, end to end through the real stage.

**Scope.** The invariant pinned here is about the ``vntyper`` package, not about
every VNtyper process. ``docker/app/tasks.py`` deliberately runs
``samtools index`` beside the upload before starting the pipeline and
deliberately removes that index in cleanup (#199); the upload volume is writable
by design, so that is not a violation and is not tested here. The read-only-input
guarantee is a guarantee about ``vntyper pipeline``.
"""

import os
import shlex
import stat
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

from vntyper.scripts import utils
from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.utils import run_command, validate_bam_file

pytestmark = pytest.mark.unit


@pytest.fixture
def readonly_dir(tmp_path: Path) -> Iterator[Path]:
    """A directory that cannot be written to, restored on teardown.

    Skipped when the test runs as root, for whom the mode bits are advisory.

    Args:
        tmp_path: Pytest temporary directory.

    Yields:
        Path: The read-only directory, holding one stub alignment.
    """
    if os.geteuid() == 0:
        pytest.skip("root ignores the write bit, so this fixture cannot deny anything")
    d = tmp_path / "input"
    d.mkdir()
    (d / "sample.bam").write_bytes(b"BAM\x01")
    d.chmod(stat.S_IRUSR | stat.S_IXUSR)
    yield d
    d.chmod(stat.S_IRWXU)


def test_run_command_raises_when_its_log_file_is_not_writable(readonly_dir: Path) -> None:
    """The mechanism: the log is opened before the child process starts.

    This is a **characterisation test, and it passes against ``main`` as well as
    against this branch** -- deliberately. ``run_command`` is not what #201
    changed; what changed is the path the pipeline hands it. So what this proves
    is only that opening the log is what fails on a read-only mount, and that it
    fails before the child runs. It proves *nothing* about where VNtyper chooses
    to write, which is the actual invariant and is the business of the two tests
    below. It is here so their failures are unambiguous when they do fail.

    It is also the one test in this module that calls ``run_command`` for real,
    and that too is deliberate: the command is ``true``, a shell builtin with no
    output and no dependency on any bioinformatics tool, so the test stays pure
    in the sense AGENTS.md requires -- no external binary, nothing downloaded,
    nothing environment-dependent about the result.

    Args:
        readonly_dir: The read-only input directory fixture.
    """
    with pytest.raises(PermissionError):
        run_command("true", str(readonly_dir / "sample.bam.quickcheck.log"))


@pytest.fixture
def recorded_quickcheck(monkeypatch: pytest.MonkeyPatch) -> list[dict[str, Any]]:
    """Capture the ``run_command`` call ``validate_bam_file`` makes, running nothing.

    These two tests are about *which path* ``validate_bam_file`` hands to
    ``run_command``, and that is decided before any process starts. Letting the
    real call through invoked the host's ``samtools`` -- verified: the log a real
    run left behind read ``sample.bam caused an error whilst reading its header``,
    which is samtools' own wording. On a host without samtools the shell's
    ``command not found`` landed in the same file and the assertion passed
    anyway, so the tier's purity was spent on a call whose outcome never mattered
    (AGENTS.md: the unit tier is ``tmp_path`` + ``unittest.mock``, no external
    binaries).

    The stub returns True, so ``validate_bam_file`` reaches its success path and
    no ``ValueError`` has to be suppressed -- another thing that used to hide
    which failure the test was actually tolerating.

    Args:
        monkeypatch: Pytest's monkeypatch fixture.

    Returns:
        list[dict[str, Any]]: One entry per ``run_command`` call, with the
        keyword names ``validate_bam_file`` uses.
    """
    calls: list[dict[str, Any]] = []

    def _record(command: str, log_file: str, critical: bool = False, cwd: str | None = None) -> bool:
        calls.append({"command": command, "log_file": log_file, "critical": critical, "cwd": cwd})
        return True

    monkeypatch.setattr(utils, "run_command", _record)
    return calls


def test_validate_bam_file_writes_its_log_under_log_dir(
    readonly_dir: Path, tmp_path: Path, recorded_quickcheck: list[dict[str, Any]]
) -> None:
    """The whole point: a read-only input mount must survive a run.

    ``run_command`` opens its log before it runs anything, so the log path is the
    whole of #201 -- asserting on the argument it received is asserting on the
    thing that broke, not on a side effect of it.

    Args:
        readonly_dir: The read-only input directory fixture.
        tmp_path: Pytest temporary directory.
        recorded_quickcheck: The ``run_command`` recorder.
    """
    out = tmp_path / "out"
    out.mkdir()

    validate_bam_file(str(readonly_dir / "sample.bam"), log_dir=str(out))

    (call,) = recorded_quickcheck
    assert call["log_file"] == str(out / "sample.bam.quickcheck.log")
    assert list(readonly_dir.iterdir()) == [readonly_dir / "sample.bam"]


def test_validate_bam_file_still_logs_beside_the_input_without_log_dir(
    tmp_path: Path, recorded_quickcheck: list[dict[str, Any]]
) -> None:
    """The default is unchanged, deliberately -- this stays a contained change.

    Args:
        tmp_path: Pytest temporary directory.
        recorded_quickcheck: The ``run_command`` recorder.
    """
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM\x01")

    validate_bam_file(str(bam))

    (call,) = recorded_quickcheck
    assert call["log_file"] == f"{bam}.quickcheck.log"


def _index_destinations(command: str) -> list[tuple[str, str]]:
    """Return ``(indexed_file, written_index)`` for every ``samtools index`` in a command.

    The naive assertion -- "no index command mentions the input directory" --
    fails after the fix as well as before it, because a correct
    ``samtools index -o <out> <in>`` must name the input as its operand. So the
    command is re-parsed with ``shlex.split`` and each ``index`` is resolved to
    the one path samtools actually writes: the argument immediately after ``-o``
    when there is one, and ``<operand>.bai`` when there is not.

    Args:
        command: A shell command, possibly several joined with ``&&``.

    Returns:
        list[tuple[str, str]]: The file being indexed and the index written for it.
    """
    destinations: list[tuple[str, str]] = []
    for segment in command.split("&&"):
        tokens = shlex.split(segment)
        if "index" not in tokens:
            continue
        operand = tokens[-1]
        if "-o" in tokens:
            written = tokens[tokens.index("-o") + 1]
            operand = tokens[tokens.index("-o") + 2]
        else:
            written = f"{operand}.bai"
        destinations.append((operand, written))
    return destinations


def test_no_index_the_bam_path_builds_is_written_outside_the_output_directory(tmp_path: Path) -> None:
    """#162/#210, end to end through the real stage.

    Every ``samtools index`` the non-fast BAM path emits is resolved to the path
    it writes, and every one of them must land in the run's output directory.
    The input index is already present in the proven plan and is never rebuilt.

    Args:
        tmp_path: Pytest temporary directory.
    """
    from unittest.mock import patch

    from vntyper.scripts import fastq_bam_processing

    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    in_bam = input_dir / "sample.bam"
    in_bam.write_bytes(b"BAM\x01")

    issued: list[str] = []

    def _record(command, log_file, critical=False, cwd=None):
        issued.append(command)
        Path(log_file).write_text("")
        return True

    with (
        patch.object(fastq_bam_processing, "run_command", _record),
        patch.object(fastq_bam_processing, "get_region_string_with_fallback", return_value="chr1:1-2"),
        patch.object(fastq_bam_processing, "extract_unmapped_reads_from_offset"),
        patch.object(fastq_bam_processing.os, "replace"),
    ):
        fastq_bam_processing.process_bam_to_fastq(
            output=str(output_dir),
            output_name="output",
            threads=4,
            config={"tools": {"samtools": "samtools"}},
            plan=AlignmentPlan(
                input_path=str(in_bam),
                view_path=str(output_dir / "input.bam"),
                file_format="bam",
                index_path=str(output_dir / "input.bam.bai"),
                reference_path=None,
                reference_source="not-required",
                uncovered_contigs=(),
                unmapped_scan="indexed",
            ),
            fast_mode=False,
        )

    indexed = [pair for command in issued for pair in _index_destinations(command)]
    assert indexed, "the non-fast BAM path emitted no index command at all; this test would pass vacuously"

    for operand, written in indexed:
        assert Path(written).parent == output_dir, (
            f"indexing {operand} wrote {written}, which is outside the run's output directory"
        )

    assert all(operand != str(in_bam) for operand, _ in indexed)
