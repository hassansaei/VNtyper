"""Which tools ``run_pipeline`` asks ``get_tool_versions`` to probe (#262).

``get_tool_versions`` used to probe every entry under ``config["tools"]``, so every
Kestrel-only run shelled out to ``mamba run -n envadvntr advntr`` (315 ms, adVNTR's own
import rather than mamba's overhead) and ``mamba run -n shark_env shark`` (36 ms) to
produce a value that is assigned once and logged once. ``tool_versions`` has exactly two
references in the codebase and reaches no summary, report or artefact, so narrowing the
probe set is structurally result-neutral.

The set is computed from the **input type and the enabled modules together**, not from
the modules alone: fastp and BWA belong to the FASTQ path and are never invoked for
BAM or CRAM input, which no module list can express.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from tests.support.pipeline_harness import MINIMAL_CONFIG, run_pipeline_under_harness

pytestmark = pytest.mark.unit


def _tools_in_use(harness) -> set[str]:
    """Return the tool-name set ``run_pipeline`` declared for this run.

    Args:
        harness: The harness returned by ``run_pipeline_under_harness``.

    Returns:
        set[str]: The declared names.
    """
    kwargs = harness.kwargs("get_tool_versions")
    assert "tools_in_use" in kwargs, "run_pipeline must declare which tools it will invoke"
    return set(kwargs["tools_in_use"])


def test_a_bam_run_declares_neither_fastp_nor_bwa(tmp_path: Path) -> None:
    """The BAM path never aligns, so neither binary is invoked."""
    harness = run_pipeline_under_harness(tmp_path / "out")

    assert _tools_in_use(harness) == {"samtools", "kestrel", "java_path"}


def test_a_cram_run_declares_the_same_set_as_a_bam_run(tmp_path: Path) -> None:
    """CRAM decodes through samtools; it does not add a tool."""
    # Outside the output tree: the pipeline refuses to write into the patient input tree.
    input_root = tmp_path / "input"
    input_root.mkdir()
    cram = input_root / "in.cram"
    cram.touch()

    harness = run_pipeline_under_harness(tmp_path / "out", cram=str(cram))

    assert _tools_in_use(harness) == {"samtools", "kestrel", "java_path"}


def test_a_fastq_run_declares_fastp_and_bwa(tmp_path: Path) -> None:
    """Modules alone cannot express this: it is the input type that adds them."""
    input_root = tmp_path / "input"
    input_root.mkdir()
    fastq1, fastq2 = input_root / "r1.fastq.gz", input_root / "r2.fastq.gz"
    fastq1.touch()
    fastq2.touch()

    harness = run_pipeline_under_harness(tmp_path / "out", fastq1=str(fastq1), fastq2=str(fastq2))

    assert _tools_in_use(harness) == {"samtools", "kestrel", "java_path", "fastp", "bwa"}


def test_requesting_advntr_declares_advntr(tmp_path: Path) -> None:
    """The probe is skipped because the module is off, not because it is never wanted.

    The shipped ``config.json`` declares ``tools.advntr``, which the harness's minimal
    config does not, so this test supplies it -- the pipeline names a module only when
    the config declares a tool of that name.
    """
    config = deepcopy(MINIMAL_CONFIG)
    config["tools"]["advntr"] = "mamba run -n envadvntr advntr"

    harness = run_pipeline_under_harness(tmp_path / "out", config=config, extra_modules=["advntr"])

    assert "advntr" in _tools_in_use(harness)


def test_a_module_the_config_does_not_declare_as_a_tool_is_not_named(tmp_path: Path) -> None:
    """``--extra-modules`` names modules; ``config["tools"]`` names tools.

    The harness config declares no ``shark`` entry, so requesting the module must not
    invent a tool name. ``get_tool_versions`` would skip it anyway, but the pipeline
    should not be asserting the existence of something it has not read.
    """
    harness = run_pipeline_under_harness(tmp_path / "out", extra_modules=["shark"])

    assert "shark" not in _tools_in_use(harness)


def test_the_declared_set_never_contains_kanalyze(tmp_path: Path) -> None:
    """kanalyze is a JAR with no version flag; it is versioned with kestrel."""
    harness = run_pipeline_under_harness(tmp_path / "out")

    assert "kanalyze" not in _tools_in_use(harness)


# ---------------------------------------------------------------------------
# The conversion stage learns whether adVNTR will read its output (#262)
# ---------------------------------------------------------------------------


def test_the_bam_conversion_is_told_adVNTR_is_off(tmp_path: Path) -> None:
    """The gate is decided in run_pipeline, which is the only place that knows."""
    harness = run_pipeline_under_harness(tmp_path / "out")

    assert harness.kwargs("process_bam_to_fastq")["needs_advntr"] is False


def test_the_bam_conversion_is_told_adVNTR_is_on(tmp_path: Path) -> None:
    harness = run_pipeline_under_harness(tmp_path / "out", extra_modules=["advntr"])

    assert harness.kwargs("process_bam_to_fastq")["needs_advntr"] is True


def test_the_post_alignment_conversion_is_told_too(tmp_path: Path) -> None:
    """The FASTQ path converts through the same stage and has the same consumer."""
    input_root = tmp_path / "input"
    input_root.mkdir()
    fastq1, fastq2 = input_root / "r1.fastq.gz", input_root / "r2.fastq.gz"
    fastq1.touch()
    fastq2.touch()

    harness = run_pipeline_under_harness(
        tmp_path / "out", fastq1=str(fastq1), fastq2=str(fastq2), extra_modules=["advntr"]
    )

    assert harness.kwargs("process_bam_to_fastq")["needs_advntr"] is True
