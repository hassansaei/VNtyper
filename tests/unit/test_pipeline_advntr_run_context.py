"""Run-owned adVNTR cleanup, model, and command bindings."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from copy import deepcopy
from pathlib import Path
from unittest.mock import DEFAULT, patch

import pytest

from tests.support.pipeline_harness import MINIMAL_CONFIG, advntr_stub, run_pipeline_under_harness
from vntyper.modules.advntr import advntr_genotyping
from vntyper.modules.advntr.model_provenance import AdvntrProbeStatus, AdvntrVersionOutcome

pytestmark = pytest.mark.unit


MODEL_SNAPSHOT_RELATIVE = Path("advntr/advntr_model.db")
PUBLIC_OUTPUTS = (
    Path("advntr/output_adVNTR_result.tsv"),
    Path("advntr/output_adVNTR.tsv"),
    Path("advntr/output_adVNTR.vcf"),
    Path("advntr/cross_match_results.tsv"),
    Path("advntr/output_advntr.log"),
    Path("pipeline_summary.json"),
    Path("pipeline_summary.csv"),
    Path("pipeline_summary.tsv"),
    Path("summary_report.html"),
)


def _write_paths(root: Path, relative_paths: tuple[Path, ...], content: bytes = b"stale patient output") -> list[Path]:
    paths = [root / relative_path for relative_path in relative_paths]
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    return paths


def _unparseable_outcome() -> AdvntrVersionOutcome:
    return AdvntrVersionOutcome(
        AdvntrProbeStatus.UNPARSEABLE_SUCCESS,
        message="adVNTR version command succeeded but its response was unparseable or ambiguous.",
    )


def test_an_early_refusal_revokes_the_exact_public_output_allowlist_only(tmp_path: Path) -> None:
    """Omitting any listed output can expose a prior patient's completed adVNTR run."""
    output = tmp_path / "out"
    stale = _write_paths(output, PUBLIC_OUTPUTS)
    preserved = _write_paths(
        output,
        (
            Path("advntr/unrelated.txt"),
            Path("kestrel/kestrel_result.tsv"),
            Path("pipeline.log"),
        ),
        b"must remain",
    )
    unselected_archives = _write_paths(tmp_path, (Path("out.zip"), Path("out.tar.gz")), b"unselected archive")

    with patch(
        "vntyper.modules.advntr.model_provenance.detect_advntr_version",
        return_value=_unparseable_outcome(),
    ):
        harness = run_pipeline_under_harness(
            output,
            extra_modules=["advntr"],
            expect_failure=True,
            log_file=str(output / "pipeline.log"),
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert [path.exists() for path in stale] == [False] * len(PUBLIC_OUTPUTS)
    assert {path.relative_to(tmp_path).as_posix(): path.read_bytes() for path in preserved} == {
        "out/advntr/unrelated.txt": b"must remain",
        "out/kestrel/kestrel_result.tsv": b"must remain",
        "out/pipeline.log": b"must remain",
    }
    assert [path.read_bytes() for path in unselected_archives] == [b"unselected archive", b"unselected archive"]
    harness.stages["run_kestrel"].assert_not_called()


def test_verified_preflight_revokes_every_public_output_before_the_probe(tmp_path: Path) -> None:
    """Cleanup must finish before a successful version/model path can proceed."""
    output = tmp_path / "out"
    stale = _write_paths(output, PUBLIC_OUTPUTS)
    existence_at_probe: list[bool] = []

    def observe_cleanup(*_args: object, **_kwargs: object) -> AdvntrVersionOutcome:
        existence_at_probe.extend(path.exists() for path in stale)
        return AdvntrVersionOutcome(AdvntrProbeStatus.VERIFIED, version=(2, 0, 4))

    with patch(
        "vntyper.modules.advntr.model_provenance.detect_advntr_version",
        side_effect=observe_cleanup,
    ):
        harness = run_pipeline_under_harness(output, extra_modules=["advntr"])

    assert harness.error is None
    assert existence_at_probe == [False] * len(PUBLIC_OUTPUTS)


def test_independently_reproduced_stale_cross_match_is_revoked_on_refusal(tmp_path: Path) -> None:
    """A cross-match from the previous run must not survive a refused new attempt."""
    output = tmp_path / "out"
    cross_match = output / "advntr" / "cross_match_results.tsv"
    cross_match.parent.mkdir(parents=True)
    cross_match.write_text("Match\nYes\n", encoding="utf-8")

    with patch(
        "vntyper.modules.advntr.model_provenance.detect_advntr_version",
        return_value=_unparseable_outcome(),
    ):
        harness = run_pipeline_under_harness(output, extra_modules=["advntr"], expect_failure=True)

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert not cross_match.exists()
    harness.stages["run_kestrel"].assert_not_called()


@pytest.mark.parametrize(
    ("archive_format", "selected_suffix", "other_suffix"),
    [("zip", ".zip", ".tar.gz"), ("tar.gz", ".tar.gz", ".zip")],
)
def test_only_the_requested_archive_destination_is_revoked(
    tmp_path: Path, archive_format: str, selected_suffix: str, other_suffix: str
) -> None:
    """A refused run withdraws its selected public archive without touching the other format."""
    output = tmp_path / "out"
    selected = Path(f"{output}{selected_suffix}")
    other = Path(f"{output}{other_suffix}")
    selected.write_bytes(b"selected stale archive")
    other.write_bytes(b"unselected archive")

    with patch(
        "vntyper.modules.advntr.model_provenance.detect_advntr_version",
        return_value=_unparseable_outcome(),
    ):
        harness = run_pipeline_under_harness(
            output,
            extra_modules=["advntr"],
            archive_results=True,
            archive_format=archive_format,
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert not selected.exists()
    assert other.read_bytes() == b"unselected archive"


def test_public_output_invalidation_failure_aborts_before_probe_or_stage_work(tmp_path: Path) -> None:
    """A public name that cannot be withdrawn must stop rather than advertise stale data."""
    output = tmp_path / "out"
    cross_match = output / "advntr" / "cross_match_results.tsv"
    cross_match.parent.mkdir(parents=True)
    cross_match.write_text("stale\n", encoding="utf-8")
    real_unlink = Path.unlink

    def fail_cross_match(path: Path, missing_ok: bool = False) -> None:
        if path == cross_match:
            raise PermissionError("injected cross-match refusal")
        real_unlink(path, missing_ok=missing_ok)

    with (
        patch.object(Path, "unlink", autospec=True, side_effect=fail_cross_match),
        patch("vntyper.modules.advntr.model_provenance.detect_advntr_version") as detector,
    ):
        harness = run_pipeline_under_harness(output, extra_modules=["advntr"], expect_failure=True)

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert cross_match.read_text(encoding="utf-8") == "stale\n"
    detector.assert_not_called()
    harness.stages["run_kestrel"].assert_not_called()


def test_selected_archive_invalidation_failure_aborts_before_probe_or_stage_work(tmp_path: Path) -> None:
    """A selected archive that cannot be withdrawn must stop before any new result work."""
    with (
        patch(
            "vntyper.scripts.pipeline_advntr_run_context.revoke_public_archive",
            side_effect=OSError("injected archive revocation failure"),
        ) as revoker,
        patch("vntyper.modules.advntr.model_provenance.detect_advntr_version") as detector,
    ):
        harness = run_pipeline_under_harness(
            tmp_path / "out",
            extra_modules=["advntr"],
            archive_results=True,
            archive_format="zip",
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    revoker.assert_called_once()
    detector.assert_not_called()
    harness.stages["run_kestrel"].assert_not_called()


def test_ownership_refusal_precedes_cross_match_invalidation(tmp_path: Path) -> None:
    """Cleanup must not unlink a public name that aliases the selected patient BAM."""
    output = tmp_path / "out"
    patient = tmp_path / "patient.bam"
    patient.write_bytes(b"operator patient input")
    cross_match = output / "advntr" / "cross_match_results.tsv"
    cross_match.parent.mkdir(parents=True)
    cross_match.hardlink_to(patient)

    with patch("vntyper.modules.advntr.model_provenance.detect_advntr_version") as detector:
        harness = run_pipeline_under_harness(
            output,
            bam=str(patient),
            extra_modules=["advntr"],
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert patient.read_bytes() == b"operator patient input"
    assert cross_match.read_bytes() == b"operator patient input"
    detector.assert_not_called()
    harness.stages["run_kestrel"].assert_not_called()


def test_model_snapshot_destination_is_owned_before_any_replacement(tmp_path: Path) -> None:
    """The fixed snapshot name cannot be a second link to the operator's model."""
    output = tmp_path / "out"
    source = tmp_path / "operator-model.db"
    source.write_bytes(Path(MINIMAL_CONFIG["reference_data"]["advntr_reference_vntr_hg19"]).read_bytes())
    snapshot = output / MODEL_SNAPSHOT_RELATIVE
    snapshot.parent.mkdir(parents=True)
    snapshot.hardlink_to(source)
    config = deepcopy(MINIMAL_CONFIG)
    config["reference_data"]["advntr_reference_vntr_hg19"] = str(source)

    with patch("vntyper.modules.advntr.model_provenance.detect_advntr_version") as detector:
        harness = run_pipeline_under_harness(
            output,
            config=config,
            extra_modules=["advntr"],
            expect_failure=True,
        )

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert source.stat().st_ino == snapshot.stat().st_ino
    detector.assert_not_called()
    harness.stages["run_kestrel"].assert_not_called()


def test_source_replacement_during_kestrel_cannot_change_snapshot_execution_or_provenance(tmp_path: Path) -> None:
    """The model bytes hashed at startup must be the same bytes handed to adVNTR later."""
    output = tmp_path / "out"
    source = tmp_path / "operator-model.db"
    original_bytes = Path(MINIMAL_CONFIG["reference_data"]["advntr_reference_vntr_hg19"]).read_bytes()
    source.write_bytes(original_bytes)
    replacement = tmp_path / "replacement.db"
    replacement.write_bytes(b"replacement bytes that are not the validated model")
    config = deepcopy(MINIMAL_CONFIG)
    config["reference_data"]["advntr_reference_vntr_hg19"] = str(source)

    def replace_source(*_args: object, **_kwargs: object) -> object:
        os.replace(replacement, source)
        return DEFAULT

    harness = run_pipeline_under_harness(
        output,
        config=config,
        extra_modules=["advntr"],
        stage_side_effects={"run_kestrel": replace_source},
    )

    snapshot = output / MODEL_SNAPSHOT_RELATIVE
    model_argument = Path(harness.positional("run_advntr")[0])
    summary = json.loads((output / "pipeline_summary.json").read_text(encoding="utf-8"))
    provenance = summary["advntr_model"]
    assert model_argument == snapshot
    assert snapshot.read_bytes() == original_bytes
    assert source.read_bytes() == b"replacement bytes that are not the validated model"
    assert provenance["path"] == str(snapshot)
    assert provenance["source_path"] == str(source)
    assert provenance["sha256"] == hashlib.sha256(original_bytes).hexdigest()
    assert hashlib.sha256(model_argument.read_bytes()).hexdigest() == provenance["sha256"]


def test_caller_config_mutation_during_kestrel_cannot_change_the_verified_command(tmp_path: Path) -> None:
    """Late execution must use the command whose 2.0.4 answer startup already logged."""
    output = tmp_path / "out"
    config = deepcopy(MINIMAL_CONFIG)
    original_command = advntr_stub("2.0.4", tmp_path)
    config["tools"]["advntr"] = original_command
    mutated_command = "/unexpected/replacement-advntr"
    launched: list[str] = []

    def mutate_config(*_args: object, **_kwargs: object) -> object:
        config["tools"]["advntr"] = mutated_command
        return DEFAULT

    def record_command(command: str, *_args: object, **_kwargs: object) -> bool:
        launched.append(command)
        return True

    with patch.object(advntr_genotyping, "run_command", side_effect=record_command):
        harness = run_pipeline_under_harness(
            output,
            config=config,
            extra_modules=["advntr"],
            stage_side_effects={
                "run_kestrel": mutate_config,
                "run_advntr": advntr_genotyping.run_advntr,
            },
        )

    assert harness.error is None
    assert config["tools"]["advntr"] == mutated_command
    assert len(launched) == 1
    assert launched[0].startswith(f"{original_command} genotype ")
    assert mutated_command not in launched[0]
    assert harness.kwargs("get_tool_versions")["version_overrides"] == {"advntr": "2.0.4"}
    execution_config = harness.kwargs("run_advntr")["config"]
    assert execution_config["tools"]["advntr"] == original_command


def test_model_copy_failure_aborts_before_probe_or_kestrel(tmp_path: Path) -> None:
    """A partial or failed copy can never fall back to executing the source model."""
    with (
        patch.object(shutil, "copyfileobj", side_effect=OSError("injected model copy failure")),
        patch("vntyper.modules.advntr.model_provenance.detect_advntr_version") as detector,
    ):
        harness = run_pipeline_under_harness(tmp_path / "out", extra_modules=["advntr"], expect_failure=True)

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    detector.assert_not_called()
    harness.stages["run_kestrel"].assert_not_called()


def test_model_snapshot_validation_failure_aborts_before_probe_or_kestrel(tmp_path: Path) -> None:
    """Only a validated snapshot can enter the run context or authorize stage work."""
    output = tmp_path / "out"
    source = tmp_path / "invalid-model.db"
    source.write_bytes(b"not an adVNTR SQLite model")
    config = deepcopy(MINIMAL_CONFIG)
    config["reference_data"]["advntr_reference_vntr_hg19"] = str(source)

    with patch("vntyper.modules.advntr.model_provenance.detect_advntr_version") as detector:
        harness = run_pipeline_under_harness(output, config=config, extra_modules=["advntr"], expect_failure=True)

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert not (output / MODEL_SNAPSHOT_RELATIVE).exists()
    detector.assert_not_called()
    harness.stages["run_kestrel"].assert_not_called()


def test_model_snapshot_install_failure_aborts_before_probe_or_kestrel(tmp_path: Path) -> None:
    """Failure to atomically land the validated bytes cannot reopen the source later."""
    output = tmp_path / "out"
    snapshot = output / MODEL_SNAPSHOT_RELATIVE
    real_replace = os.replace

    def fail_snapshot_install(
        source: str | bytes | os.PathLike[str] | os.PathLike[bytes],
        destination: str | bytes | os.PathLike[str] | os.PathLike[bytes],
    ) -> None:
        if Path(os.fsdecode(destination)) == snapshot:
            raise OSError("injected model install failure")
        real_replace(source, destination)

    with (
        patch.object(os, "replace", side_effect=fail_snapshot_install),
        patch("vntyper.modules.advntr.model_provenance.detect_advntr_version") as detector,
    ):
        harness = run_pipeline_under_harness(output, extra_modules=["advntr"], expect_failure=True)

    assert isinstance(harness.error, SystemExit)
    assert harness.error.code == 1
    assert not snapshot.exists()
    detector.assert_not_called()
    harness.stages["run_kestrel"].assert_not_called()
