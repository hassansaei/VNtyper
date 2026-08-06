"""Wiring ``reconcile_assembly`` into the run: the enforcement half of the guard.

``assembly_guard.reconcile_assembly`` shipped as a pure function with zero callers,
deliberately, so that the behaviour-changing half could land as its own revertible
commit. This is that half.

Why it matters: the MUC1 VNTR sits ~30 kb apart in GRCh37 and GRCh38. Declaring the
wrong ``--reference-assembly`` slices a region that does not contain the VNTR, Kestrel
sees no supporting reads, and the pipeline reports a confident negative. Nothing in the
run fails - the report simply says the sample is negative.

**This rejects inputs that succeed today.** The tests below pin exactly which:

* only a *decided* disagreement between two known builds is fatal;
* an alias pair (``hg19``/``GRCh37``) names one coordinate system and agrees;
* an unreadable header, an unknown declared name, a headerless input and a chr1 that
  matches neither build are all ``undetermined``, and ``undetermined`` warns and
  continues - it is neither a pass nor a failure;
* FASTQ input is not guarded at all, because the only header available after alignment
  describes the **BWA reference**, not the sample, and a mismatch there would send the
  user to fix an input file that is not wrong.
"""

import logging
import subprocess
from pathlib import Path
from unittest import mock

import pytest

from tests.support.pipeline_harness import run_pipeline_under_harness
from vntyper.scripts import pipeline_guards
from vntyper.scripts.assembly_guard import STATUS_AGREE, STATUS_MISMATCH, STATUS_UNDETERMINED

pytestmark = pytest.mark.unit

#: chr1 lengths, which are what decides the build.
GRCH37_CHR1 = 249250621
GRCH38_CHR1 = 248956422


def _header(name: str = "chr1", length: int = GRCH37_CHR1) -> str:
    """Build a minimal SAM header with one contig.

    Args:
        name: The contig name.
        length: The contig length.

    Returns:
        str: SAM header text.
    """
    return f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{name}\tLN:{length}\n"


# --------------------------------------------------------------------------------------
# read_alignment_header
# --------------------------------------------------------------------------------------


def test_the_header_is_read_once_and_returned() -> None:
    """The happy path: one samtools invocation, the header text back."""
    with mock.patch.object(pipeline_guards, "extract_bam_header", return_value=_header()) as stub:
        assert pipeline_guards.read_alignment_header("in.bam", {}) == _header()
    assert stub.call_count == 1


@pytest.mark.parametrize(
    "error",
    [
        pytest.param(subprocess.CalledProcessError(1, "samtools"), id="called_process_error"),
        pytest.param(FileNotFoundError("samtools"), id="missing_binary"),
        pytest.param(OSError("permission denied"), id="os_error"),
    ],
)
def test_an_unreadable_header_warns_and_yields_none(error: Exception, caplog) -> None:
    """``extract_bam_header`` uses ``check=True``; a CRAM without its reference raises.

    A guard that turned an unreadable header into a crash would be worse than the
    defect it closes, so the read is non-fatal and the verdict becomes undetermined.

    Args:
        error: The exception samtools raises.
        caplog: Pytest log capture.
    """
    with (
        mock.patch.object(pipeline_guards, "extract_bam_header", side_effect=error),
        caplog.at_level(logging.WARNING),
    ):
        assert pipeline_guards.read_alignment_header("in.cram", {}) is None

    assert any(record.levelno == logging.WARNING for record in caplog.records)


# --------------------------------------------------------------------------------------
# enforce_declared_assembly
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "declared",
    ["hg19", "GRCh37", "hg19_ncbi", "hg19_ensembl"],
)
def test_every_spelling_of_a_build_agrees_with_that_builds_header(declared: str, caplog) -> None:
    """``hg19`` and ``GRCh37`` name the same coordinate system; both must agree.

    Args:
        declared: The declared assembly name.
        caplog: Pytest log capture.
    """
    with caplog.at_level(logging.DEBUG):
        verdict = pipeline_guards.enforce_declared_assembly(declared, _header(length=GRCH37_CHR1))

    assert verdict.status == STATUS_AGREE
    assert not [record for record in caplog.records if record.levelno >= logging.WARNING]


@pytest.mark.parametrize(
    ("declared", "chr1_length"),
    [
        pytest.param("hg19", GRCH38_CHR1, id="declared_37_header_38"),
        pytest.param("GRCh37", GRCH38_CHR1, id="declared_37_alias_header_38"),
        pytest.param("hg38", GRCH37_CHR1, id="declared_38_header_37"),
        pytest.param("hg38_ensembl", GRCH37_CHR1, id="declared_38_ensembl_header_37"),
    ],
)
def test_a_decided_disagreement_is_fatal(declared: str, chr1_length: int, caplog) -> None:
    """The behaviour change: a wrong ``--reference-assembly`` now stops the run.

    Args:
        declared: The declared assembly name.
        chr1_length: The chr1 length in the header.
        caplog: Pytest log capture.
    """
    header = _header(length=chr1_length)
    with caplog.at_level(logging.ERROR), pytest.raises(ValueError) as excinfo:
        pipeline_guards.enforce_declared_assembly(declared, header)

    message = str(excinfo.value)
    assert declared in message
    errors = [record for record in caplog.records if record.levelno >= logging.ERROR]
    assert errors, "a mismatch must be logged at ERROR as well as raised"
    assert errors[0].message == message, "AGENTS.md: log and raise the same message"


@pytest.mark.parametrize(
    ("declared", "header"),
    [
        pytest.param("hg19", None, id="header_unreadable"),
        pytest.param("hg19", "@HD\tVN:1.6\n", id="no_contigs"),
        pytest.param("hg19", _header(name="chr2", length=GRCH37_CHR1), id="no_chr1"),
        pytest.param("hg19", _header(length=12345), id="chr1_matches_no_build"),
        pytest.param("not_an_assembly", _header(), id="declared_name_unknown"),
    ],
)
def test_an_undetermined_verdict_warns_and_continues(declared: str, header: str | None, caplog) -> None:
    """``undetermined`` is neither a pass nor a failure, so it must not stop the run.

    Args:
        declared: The declared assembly name.
        header: The header text, or None when it could not be read.
        caplog: Pytest log capture.
    """
    with caplog.at_level(logging.WARNING):
        verdict = pipeline_guards.enforce_declared_assembly(declared, header)

    assert verdict.status == STATUS_UNDETERMINED
    assert any(record.levelno == logging.WARNING for record in caplog.records)


def test_an_ncbi_named_header_is_still_decided() -> None:
    """Agent B closed the NCBI gap; the wiring must not reintroduce it."""
    header = _header(name="NC_000001.10", length=GRCH37_CHR1)
    assert pipeline_guards.enforce_declared_assembly("hg19_ncbi", header).status == STATUS_AGREE

    with pytest.raises(ValueError):
        pipeline_guards.enforce_declared_assembly("hg38_ncbi", header)


# --------------------------------------------------------------------------------------
# The wiring inside run_pipeline
# --------------------------------------------------------------------------------------


def _run(tmp_path: Path, **kwargs):
    """Run the harness with the guard's collaborators live.

    Args:
        tmp_path: Pytest temporary directory.
        **kwargs: Forwarded to ``run_pipeline_under_harness``.

    Returns:
        PipelineHarness: The recorded run.
    """
    return run_pipeline_under_harness(tmp_path / "out", **kwargs)


def test_a_matching_bam_still_runs(tmp_path: Path) -> None:
    """The guard must not change a run whose declared assembly is right.

    Args:
        tmp_path: Pytest temporary directory.
    """
    harness = _run(tmp_path, header=_header(length=GRCH37_CHR1), reference_assembly="hg19")
    assert harness.stages["run_kestrel"].called


def test_a_mismatched_bam_stops_before_any_region_is_resolved(tmp_path: Path) -> None:
    """The guard sits before region resolution, so no wrong region is ever sliced.

    Args:
        tmp_path: Pytest temporary directory.
    """
    harness = _run(
        tmp_path,
        header=_header(length=GRCH38_CHR1),
        reference_assembly="hg19",
        expect_failure=True,
    )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert not harness.stages["get_region_string_with_fallback"].called, "a region was resolved after a mismatch"
    assert not harness.stages["process_bam_to_fastq"].called
    assert not harness.stages["run_kestrel"].called


def test_a_mismatched_cram_stops_too(tmp_path: Path) -> None:
    """CRAM shares the BAM code path, at the same insertion point.

    Args:
        tmp_path: Pytest temporary directory.
    """
    cram = tmp_path / "in.cram"
    cram.parent.mkdir(parents=True, exist_ok=True)
    cram.touch()
    harness = _run(
        tmp_path,
        bam=None,
        cram=str(cram),
        header=_header(length=GRCH38_CHR1),
        reference_assembly="hg19",
        expect_failure=True,
    )

    assert isinstance(harness.error, SystemExit)
    assert not harness.stages["process_bam_to_fastq"].called


def test_a_cram_whose_header_cannot_be_read_still_runs(tmp_path: Path) -> None:
    """The case the try/except exists for: an unresolvable CRAM reference.

    ``extract_bam_header`` runs samtools with ``check=True``, so a CRAM whose
    reference cannot be found raises ``CalledProcessError`` - which is neither
    ``KeyError`` nor ``ValueError``, so the existing region fallback does not catch it.

    Args:
        tmp_path: Pytest temporary directory.
    """
    cram = tmp_path / "in.cram"
    cram.parent.mkdir(parents=True, exist_ok=True)
    cram.touch()
    harness = _run(
        tmp_path,
        bam=None,
        cram=str(cram),
        header_error=subprocess.CalledProcessError(1, "samtools view -H"),
        reference_assembly="hg19",
    )

    assert harness.stages["run_kestrel"].called


def test_the_header_is_read_once_for_a_bam(tmp_path: Path) -> None:
    """Hoisting the read must not add a second samtools invocation.

    ``parse_header_pipeline_info`` already needed the header; the guard reuses it.

    Args:
        tmp_path: Pytest temporary directory.
    """
    harness = _run(tmp_path, header=_header(length=GRCH37_CHR1), reference_assembly="hg19")

    guard_reads = harness.stages["read_alignment_header"].call_count
    fallback_reads = harness.stages["extract_bam_header"].call_count
    assert guard_reads == 1
    assert fallback_reads == 0, "the header was read a second time on the happy path"
    assert harness.stages["parse_header_pipeline_info"].called


def test_an_unreadable_bam_header_still_fails_the_header_parsing_step(tmp_path: Path) -> None:
    """The guard must not turn a samtools failure into a silently skipped step.

    Before the guard existed, an ``extract_bam_header`` failure propagated out of the
    BAM branch and ended the run. Swallowing it in the guard and then skipping
    ``parse_header_pipeline_info`` would be a new, quieter behaviour; instead the BAM
    branch re-reads and the original failure semantics are preserved.

    Args:
        tmp_path: Pytest temporary directory.
    """
    harness = _run(
        tmp_path,
        header_error=subprocess.CalledProcessError(1, "samtools view -H"),
        reference_assembly="hg19",
        expect_failure=True,
    )

    assert isinstance(harness.error, SystemExit)
    assert harness.stages["extract_bam_header"].called, "the BAM branch did not re-read the header"


def test_fastq_input_is_not_guarded(tmp_path: Path) -> None:
    """FASTQ has no header of its own, so the guard deliberately does not run.

    The only header available after alignment describes the reference that was
    indexed for BWA, not the sample. A mismatch there is a real misconfiguration
    with the same false-negative consequence, but the remediation is different, so
    enforcing it with this message would send users to fix an input that is not
    wrong. Left to a follow-up that can word it properly; pinned here so that
    adding it later is a visible change.

    Args:
        tmp_path: Pytest temporary directory.
    """
    harness = _run(
        tmp_path,
        bam=None,
        fastq1="/in/r1.fastq.gz",
        fastq2="/in/r2.fastq.gz",
        header=_header(length=GRCH38_CHR1),
        reference_assembly="hg19",
    )

    assert harness.stages["run_kestrel"].called, "a FASTQ run must not be stopped by the assembly guard"
    assert not harness.stages["read_alignment_header"].called


def test_the_guard_runs_after_input_validation(tmp_path: Path) -> None:
    """An invalid BAM must still be reported as an invalid BAM, not as a mismatch.

    Args:
        tmp_path: Pytest temporary directory.
    """
    harness = _run(
        tmp_path,
        header=_header(length=GRCH38_CHR1),
        reference_assembly="hg19",
        stage_side_effects={"validate_bam_file": mock.Mock(side_effect=ValueError("not a BAM"))},
        expect_failure=True,
    )

    assert isinstance(harness.error, SystemExit)
    assert not harness.stages["read_alignment_header"].called


def test_the_pure_verdict_function_still_has_no_policy() -> None:
    """``assembly_guard`` stays pure; the policy lives here and is revertible alone."""
    import ast
    import inspect

    from vntyper.scripts import assembly_guard

    tree = ast.parse(inspect.getsource(assembly_guard))
    raises = [node for node in ast.walk(tree) if isinstance(node, ast.Raise)]
    assert not raises, "assembly_guard must stay non-fatal; enforcement belongs in pipeline_guards"

    statuses = {STATUS_AGREE, STATUS_MISMATCH, STATUS_UNDETERMINED}
    assert len(statuses) == 3, "the three statuses collapsed; the tests above would be vacuous"
