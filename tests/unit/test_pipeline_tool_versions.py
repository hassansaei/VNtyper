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

import json
import subprocess
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

import pytest

from tests.support.pipeline_harness import MINIMAL_CONFIG, advntr_stub, run_pipeline_under_harness
from vntyper.modules.advntr.model_provenance import AdvntrProbeStatus, AdvntrVersionOutcome

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
    # A stub reporting a span-aware adVNTR. The real `mamba run -n envadvntr advntr`
    # is refused before the run when the installed adVNTR predates the recorded genomic
    # end (#268) -- that check is the point, but this test is about tool declaration.
    config["tools"]["advntr"] = advntr_stub("2.0.4")

    harness = run_pipeline_under_harness(tmp_path / "out", config=config, extra_modules=["advntr"])

    assert "advntr" in _tools_in_use(harness)


def test_requesting_advntr_shares_the_verified_version_with_startup_logging(tmp_path: Path) -> None:
    """Dropping the override would make startup logging launch a redundant mamba probe."""
    config = deepcopy(MINIMAL_CONFIG)
    config["tools"]["advntr"] = advntr_stub("2.0.4")

    harness = run_pipeline_under_harness(tmp_path / "out", config=config, extra_modules=["advntr"])

    assert harness.kwargs("get_tool_versions")["version_overrides"] == {"advntr": "2.0.4"}


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


def test_the_model_a_run_used_is_recorded_in_the_summary(tmp_path: Path) -> None:
    """A result must be traceable to the model that produced it (#268).

    The fetch window comes from the model's own content, so which file was resolved
    decides which reads adVNTR could ever see. Recording it is what makes a past run
    interpretable.
    """
    config = deepcopy(MINIMAL_CONFIG)
    config["tools"]["advntr"] = advntr_stub("2.0.4")

    harness = run_pipeline_under_harness(tmp_path / "out", config=config, extra_modules=["advntr"])

    summary = json.loads((harness.output_dir / "pipeline_summary.json").read_text())
    model = summary["advntr_model"]
    assert model["schema_version"] == "v2"
    assert model["window_bp"] == 155192032 - 155188507
    assert model["genomic_interval"] == "chr1:155188507-155192032"
    assert len(model["sha256"]) == 64


def test_a_model_advntr_cannot_read_stops_the_run(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """An adVNTR predating the recorded end would fetch the old 840 bp window.

    It reports no error while doing so, so the run has to be refused up front rather
    than diagnosed from results that look confident.
    """
    config = deepcopy(MINIMAL_CONFIG)
    incompatible = subprocess.CompletedProcess(
        [config["tools"]["advntr"], "--version"],
        0,
        stdout="",
        stderr="2.0.3\n",
    )

    # run_pipeline swallows every exception into sys.exit(1), so what the harness records
    # is the exit; the refusal itself is in the log. Both matter: the run must stop, and
    # it must say why.
    caplog.set_level("ERROR")
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=incompatible,
    ) as runner:
        harness = run_pipeline_under_harness(
            tmp_path / "out",
            config=config,
            extra_modules=["advntr"],
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert runner.call_count == 1
    harness.stages["run_kestrel"].assert_not_called()
    assert "Install adVNTR >= 2.0.4" in caplog.text


def test_a_203_nonzero_tagged_banner_reaches_the_upgrade_refusal(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """The old CLI's banner outranks lock noise and stops before Kestrel."""
    config = deepcopy(MINIMAL_CONFIG)
    config["tools"]["advntr"] = "mamba run -n envadvntr advntr"
    incompatible = subprocess.CompletedProcess(
        ["mamba", "run", "-n", "envadvntr", "advntr", "--version"],
        2,
        stdout="",
        stderr=(
            "warning  libmamba Cannot lock '/home/test/.cache/mamba/proc'\n"
            "    Waiting for other mamba process to finish\n"
            "usage: \n"
            "=======================================================\n"
            "adVNTR 2.0.3: Genopyting tool for VNTRs\n"
            "=======================================================\n"
            "usage: advntr <command> [options]\n"
            "advntr: error: too few arguments\n"
        ),
    )

    caplog.set_level("ERROR")
    with (
        patch(
            "vntyper.modules.advntr.model_provenance.subprocess.run",
            return_value=incompatible,
        ) as runner,
        patch("vntyper.modules.advntr.model_provenance.time.sleep") as sleep,
    ):
        harness = run_pipeline_under_harness(
            tmp_path / "out",
            config=config,
            extra_modules=["advntr"],
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert runner.call_count == 1
    sleep.assert_not_called()
    harness.stages["run_kestrel"].assert_not_called()
    assert "Install adVNTR >= 2.0.4" in caplog.text
    assert "command exited with status 2" not in caplog.text


def test_a_malformed_startup_version_stops_before_kestrel_without_reprobing(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """The classified startup outcome reaches the early guard exactly once."""
    config = deepcopy(MINIMAL_CONFIG)
    argv = [config["tools"]["advntr"], "--version"]
    malformed = subprocess.CompletedProcess(
        argv,
        0,
        stdout="usage: advntr [options]\n",
        stderr="Python 3.12.1 required\n",
    )

    caplog.set_level("ERROR")
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=malformed,
    ) as runner:
        harness = run_pipeline_under_harness(
            tmp_path / "out",
            config=config,
            extra_modules=["advntr"],
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert runner.call_count == 1
    harness.stages["run_kestrel"].assert_not_called()
    expected = "adVNTR version command succeeded but its response was unparseable or ambiguous."
    assert any(record.getMessage() == expected for record in caplog.records)


def test_a_permanent_launch_exception_stops_before_kestrel_with_its_own_error(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """An environment launch failure remains distinct from a malformed tool answer."""
    caplog.set_level("ERROR")
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        side_effect=FileNotFoundError("missing"),
    ) as runner:
        harness = run_pipeline_under_harness(
            tmp_path / "out",
            extra_modules=["advntr"],
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert runner.call_count == 1
    harness.stages["run_kestrel"].assert_not_called()
    assert any(
        record.getMessage().startswith("adVNTR version launch failed: command not found:") for record in caplog.records
    )


def test_a_permanent_nonzero_exit_stops_before_kestrel_with_its_status(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A completed launcher failure carries its status without becoming malformed output."""
    config = deepcopy(MINIMAL_CONFIG)
    failed = subprocess.CompletedProcess([config["tools"]["advntr"], "--version"], 127, stdout="", stderr="missing\n")

    caplog.set_level("ERROR")
    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=failed,
    ) as runner:
        harness = run_pipeline_under_harness(
            tmp_path / "out",
            config=config,
            extra_modules=["advntr"],
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert runner.call_count == 1
    harness.stages["run_kestrel"].assert_not_called()
    expected = "adVNTR version launch failed: command exited with status 127."
    assert any(record.getMessage() == expected for record in caplog.records)


def test_transient_exhaustion_stops_before_kestrel_with_a_bounded_launch_count(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Exhaustion is its own outcome after exactly three attempts and two waits."""
    config = deepcopy(MINIMAL_CONFIG)
    config["tools"]["advntr"] = "mamba run -n envadvntr advntr"
    argv = ["mamba", "run", "-n", "envadvntr", "advntr", "--version"]
    transient = subprocess.CompletedProcess(
        argv,
        1,
        stdout="",
        stderr="error libmamba Could not set lock (Resource temporarily unavailable)\n",
    )

    caplog.set_level("ERROR")
    with (
        patch(
            "vntyper.modules.advntr.model_provenance.subprocess.run",
            return_value=transient,
        ) as runner,
        patch("vntyper.modules.advntr.model_provenance.time.sleep") as sleep,
    ):
        harness = run_pipeline_under_harness(
            tmp_path / "out",
            config=config,
            extra_modules=["advntr"],
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert runner.call_count == 3
    assert sleep.call_count == 2
    harness.stages["run_kestrel"].assert_not_called()
    expected = "adVNTR version detection exhausted 3 attempts after transient mamba launch failures."
    assert any(record.getMessage() == expected for record in caplog.records)


def test_a_verified_version_launches_once_and_reaches_kestrel(tmp_path: Path) -> None:
    """The successful 2.0.4 path remains unchanged by terminal-outcome classification."""
    config = deepcopy(MINIMAL_CONFIG)
    verified = subprocess.CompletedProcess(
        [config["tools"]["advntr"], "--version"],
        0,
        stdout="",
        stderr="2.0.4\n",
    )

    with patch(
        "vntyper.modules.advntr.model_provenance.subprocess.run",
        return_value=verified,
    ) as runner:
        harness = run_pipeline_under_harness(tmp_path / "out", config=config, extra_modules=["advntr"])

    assert harness.error is None
    assert runner.call_count == 1
    harness.stages["run_kestrel"].assert_called_once()


def test_a_forged_invalid_verified_outcome_stops_before_compatibility_and_kestrel(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Defense in depth rejects a producer invariant breach at the pipeline boundary."""
    invalid = object.__new__(AdvntrVersionOutcome)
    object.__setattr__(invalid, "status", AdvntrProbeStatus.VERIFIED)
    object.__setattr__(invalid, "version", (99,))
    object.__setattr__(invalid, "message", "")
    caplog.set_level("ERROR")

    with (
        patch(
            "vntyper.modules.advntr.model_provenance.detect_advntr_version",
            return_value=invalid,
        ),
        patch("vntyper.modules.advntr.model_provenance.require_compatible_advntr") as compatibility,
    ):
        harness = run_pipeline_under_harness(
            tmp_path / "out",
            extra_modules=["advntr"],
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    compatibility.assert_not_called()
    harness.stages["run_kestrel"].assert_not_called()
    assert "version must be exactly three non-negative integers" in caplog.text
