"""A crashed adVNTR command must stop the pipeline before output parsing."""

import logging
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

import pytest

from tests.support.pipeline_harness import MINIMAL_CONFIG, run_pipeline_under_harness
from vntyper.modules.advntr import advntr_genotyping as advntr
from vntyper.modules.advntr.model_provenance import AdvntrProbeStatus, AdvntrVersionOutcome

pytestmark = pytest.mark.unit

MAIN_CONFIG = {"tools": {"advntr": "mamba run -n envadvntr advntr"}}

STALE_ADVNTR_ARTIFACTS = (
    "output_adVNTR_result.tsv",
    "output_adVNTR.tsv",
    "output_adVNTR.vcf",
)


@pytest.fixture
def advntr_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Return existing database, alignment, and output paths."""
    database = tmp_path / "vntr.db"
    database.write_text("database", encoding="utf-8")
    alignment = tmp_path / "sample.bam"
    alignment.write_text("alignment", encoding="utf-8")
    output = tmp_path / "advntr"
    output.mkdir()
    return database, alignment, output


def test_run_advntr_propagates_a_critical_command_failure(
    advntr_inputs: tuple[Path, Path, Path], monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """The stage must not translate a failed critical command into an ignored integer."""
    database, alignment, output = advntr_inputs
    failure = RuntimeError("Critical command failed: injected adVNTR")
    monkeypatch.setattr(advntr, "run_command", lambda *_args, **_kwargs: (_ for _ in ()).throw(failure))

    with caplog.at_level(logging.ERROR, logger=advntr.logger.name), pytest.raises(RuntimeError) as raised:
        advntr.run_advntr(str(database), str(alignment), str(output), "output", MAIN_CONFIG)

    assert raised.value is failure
    assert "adVNTR genotyping failed" in caplog.text


def test_run_advntr_rejects_a_false_critical_command_result(
    advntr_inputs: tuple[Path, Path, Path], monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """A future non-raising critical runner must not restore the silent failure."""
    database, alignment, output = advntr_inputs
    monkeypatch.setattr(advntr, "run_command", lambda *_args, **_kwargs: False)

    with (
        caplog.at_level(logging.ERROR, logger=advntr.logger.name),
        pytest.raises(RuntimeError, match="adVNTR genotyping failed"),
    ):
        advntr.run_advntr(str(database), str(alignment), str(output), "output", MAIN_CONFIG)

    assert "adVNTR genotyping failed" in caplog.text


def test_a_nonzero_advntr_command_stops_the_pipeline_before_parsing(tmp_path: Path) -> None:
    """An injected header-writing crash must produce exit 1, never a Negative result."""
    fake_advntr = tmp_path / "fake_advntr"
    fake_advntr.write_text(
        """#!/usr/bin/env python3
import pathlib
import sys

if "--version" in sys.argv:
    print("2.0.4")
    raise SystemExit(0)

output = pathlib.Path(sys.argv[sys.argv.index("-o") + 1])
output.write_text("#VID\\tState\\tNumberOfSupportingReads\\tMeanCoverage\\tPvalue\\n", encoding="utf-8")
raise SystemExit(23)
""",
        encoding="utf-8",
    )
    fake_advntr.chmod(0o755)
    config = deepcopy(MINIMAL_CONFIG)
    config["tools"]["advntr"] = str(fake_advntr)
    real_run_advntr = advntr.run_advntr

    harness = run_pipeline_under_harness(
        tmp_path / "out",
        config=config,
        extra_modules=["advntr"],
        expect_failure=True,
        stage_side_effects={"run_advntr": real_run_advntr},
    )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert harness.stages["run_advntr"].called
    assert not harness.stages["process_advntr_output"].called
    assert (harness.output_dir / "advntr" / "output_adVNTR.vcf").read_text(encoding="utf-8").startswith("#VID")
    assert not (harness.output_dir / "advntr" / "output_adVNTR_result.tsv").exists()


def test_a_nonzero_advntr_return_status_stops_before_result_parsing(tmp_path: Path) -> None:
    """Pre-command validation failures return 1 and must not reach any result consumer."""
    harness = run_pipeline_under_harness(
        tmp_path / "out",
        extra_modules=["advntr"],
        expect_failure=True,
        stage_side_effects={"run_advntr": lambda *_args, **_kwargs: 1},
    )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert harness.stages["run_advntr"].called
    assert not harness.stages["process_advntr_output"].called
    assert not harness.stages["extract_results_from_pipeline_summary"].called
    assert not harness.stages["cross_match_variants"].called
    assert not harness.stages["generate_summary_report"].called


def test_a_cross_format_rerun_invalidates_all_preexisting_artifacts_before_execution(tmp_path: Path) -> None:
    """Neither raw format nor the derived call may survive into a new adVNTR attempt."""
    output_dir = tmp_path / "out"
    advntr_dir = output_dir / "advntr"
    artifacts = {
        "tsv": advntr_dir / "output_adVNTR.tsv",
        "vcf": advntr_dir / "output_adVNTR.vcf",
        "derived": advntr_dir / "output_adVNTR_result.tsv",
    }
    advntr_dir.mkdir(parents=True)
    patient_call = (
        "#VID\tState\tNumberOfSupportingReads\tMeanCoverage\tPvalue\n25561\tI22_2_G_LEN1\t11\t153.98\t0.0001\n"
    )
    artifacts["tsv"].write_text(patient_call, encoding="utf-8")
    artifacts["vcf"].write_text(patient_call, encoding="utf-8")
    artifacts["derived"].write_text("VID\tVariant\n25561\tI22_2_G_LEN1\n", encoding="utf-8")
    existence_at_execution: dict[str, bool] = {}

    def observe_then_fail(*_args, **_kwargs) -> int:
        existence_at_execution.update({name: path.exists() for name, path in artifacts.items()})
        return 1

    harness = run_pipeline_under_harness(
        output_dir,
        extra_modules=["advntr"],
        expect_failure=True,
        stage_side_effects={"run_advntr": observe_then_fail},
    )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert harness.stages["run_advntr"].called
    assert not harness.stages["process_advntr_output"].called
    assert existence_at_execution == {"tsv": False, "vcf": False, "derived": False}
    assert not any(path.exists() for path in artifacts.values())


@pytest.mark.parametrize(
    "outcome",
    [
        AdvntrVersionOutcome(
            AdvntrProbeStatus.VERSIONED_LAUNCH_FAILURE,
            version=(2, 0, 4),
            message="adVNTR version launch failed: command exited with status 2.",
        ),
        AdvntrVersionOutcome(
            AdvntrProbeStatus.LAUNCH_FAILURE,
            message="adVNTR version launch failed: command exited with status 127.",
        ),
        AdvntrVersionOutcome(
            AdvntrProbeStatus.UNPARSEABLE_SUCCESS,
            message="adVNTR version command succeeded but its response was unparseable or ambiguous.",
        ),
        AdvntrVersionOutcome(
            AdvntrProbeStatus.TRANSIENT_EXHAUSTED,
            message="adVNTR version detection exhausted 3 attempts after transient mamba launch failures.",
        ),
    ],
    ids=lambda outcome: outcome.status.value,
)
def test_every_nonverified_startup_outcome_invalidates_all_prior_advntr_artifacts(
    tmp_path: Path, outcome: AdvntrVersionOutcome
) -> None:
    """An early probe refusal must not leave a previous patient's call looking current."""
    output_dir = tmp_path / "out"
    advntr_dir = output_dir / "advntr"
    advntr_dir.mkdir(parents=True)
    artifacts = [advntr_dir / name for name in STALE_ADVNTR_ARTIFACTS]
    for artifact in artifacts:
        artifact.write_text("stale prior patient result\n", encoding="utf-8")

    with patch(
        "vntyper.modules.advntr.model_provenance.detect_advntr_version",
        return_value=outcome,
    ):
        harness = run_pipeline_under_harness(
            output_dir,
            extra_modules=["advntr"],
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert not any(artifact.exists() for artifact in artifacts)
    harness.stages["run_kestrel"].assert_not_called()


def test_verified_startup_also_invalidates_all_prior_advntr_artifacts(tmp_path: Path) -> None:
    """The cleanup move must preserve the established successful-run invalidation."""
    output_dir = tmp_path / "out"
    advntr_dir = output_dir / "advntr"
    advntr_dir.mkdir(parents=True)
    artifacts = [advntr_dir / name for name in STALE_ADVNTR_ARTIFACTS]
    for artifact in artifacts:
        artifact.write_text("stale prior patient result\n", encoding="utf-8")

    outcome = AdvntrVersionOutcome(AdvntrProbeStatus.VERIFIED, version=(2, 0, 4))
    with patch(
        "vntyper.modules.advntr.model_provenance.detect_advntr_version",
        return_value=outcome,
    ):
        harness = run_pipeline_under_harness(output_dir, extra_modules=["advntr"])

    assert harness.error is None
    assert not any(artifact.exists() for artifact in artifacts)


def test_input_ownership_refusal_precedes_advntr_artifact_invalidation(tmp_path: Path) -> None:
    """An input aliased into the output tree must be protected before any cleanup write."""
    output_dir = tmp_path / "out"
    advntr_dir = output_dir / "advntr"
    advntr_dir.mkdir(parents=True)
    input_bam = tmp_path / "patient.bam"
    input_bam.write_text("protected patient input\n", encoding="utf-8")
    protected_artifact = advntr_dir / STALE_ADVNTR_ARTIFACTS[0]
    protected_artifact.hardlink_to(input_bam)
    other_artifacts = [advntr_dir / name for name in STALE_ADVNTR_ARTIFACTS[1:]]
    for artifact in other_artifacts:
        artifact.write_text("stale prior patient result\n", encoding="utf-8")

    with patch("vntyper.modules.advntr.model_provenance.detect_advntr_version") as detector:
        harness = run_pipeline_under_harness(
            output_dir,
            bam=str(input_bam),
            extra_modules=["advntr"],
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    detector.assert_not_called()
    assert input_bam.read_text(encoding="utf-8") == "protected patient input\n"
    assert protected_artifact.exists()
    assert all(artifact.exists() for artifact in other_artifacts)


def test_a_success_status_without_current_raw_output_cannot_republish_a_preexisting_call(tmp_path: Path) -> None:
    """Status zero is insufficient when the producer leaves only a previous run's raw TSV."""
    output_dir = tmp_path / "out"
    raw_path = output_dir / "advntr" / "output_adVNTR.tsv"
    result_path = output_dir / "advntr" / "output_adVNTR_result.tsv"
    raw_path.parent.mkdir(parents=True)
    raw_path.write_text(
        "#VID\tState\tNumberOfSupportingReads\tMeanCoverage\tPvalue\n25561\tI22_2_G_LEN1\t11\t153.98\t0.0001\n",
        encoding="utf-8",
    )

    harness = run_pipeline_under_harness(
        output_dir,
        extra_modules=["advntr"],
        expect_failure=True,
        stage_side_effects={
            "run_advntr": lambda *_args, **_kwargs: 0,
            "process_advntr_output": advntr.process_advntr_output,
        },
    )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert harness.stages["run_advntr"].called
    assert harness.stages["process_advntr_output"].called
    assert not raw_path.exists(), "the new attempt must invalidate the previous producer artifact"
    assert not result_path.exists(), "the stale raw call must not be republished as this run's derived result"
