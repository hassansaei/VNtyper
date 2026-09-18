"""Single-consumer routing preserves existing length contracts."""

from pathlib import Path
from unittest.mock import Mock

import pytest

from vntyper.scripts import pipeline_length_routing as subject
from vntyper.scripts.pipeline_standard_length import StandardLengthConfiguration, StandardLengthRunner

pytestmark = pytest.mark.unit


def test_legacy_default_configuration_remains_disabled() -> None:
    result = subject.validate_pipeline_length(None, ())
    assert result.measurement_enabled is False
    assert subject.validate_pipeline_length(result, (Path("model.json"),)) == result


@pytest.mark.parametrize("configuration,paths", [(object(), ()), (None, []), (None, (17,))])
def test_invalid_legacy_configuration_or_paths_fail(configuration, paths) -> None:
    with pytest.raises(ValueError):
        subject.validate_pipeline_length(configuration, paths)


@pytest.mark.parametrize("approved,standard", [(False, False), (True, False), (False, True)])
def test_router_selects_exactly_one_reader(approved: bool, standard: bool) -> None:
    factory = Mock()
    configuration = Mock(measurement_enabled=approved)
    standard_config = StandardLengthConfiguration(standard, None, None, "a" * 64)
    result = subject.build_length_consumer(
        configuration,
        standard_config,
        assembly="hg38",
        reference="reference.fa",
        project_root="/run",
        samtools="samtools",
        approved_factory=factory,
    )
    if approved:
        assert result is factory.return_value
        factory.assert_called_once_with(
            configuration=configuration, bwa_reference="reference.fa", project_root="/run", samtools_path="samtools"
        )
    elif standard:
        assert isinstance(result, StandardLengthRunner)
        assert result.project_root == Path("/run")
        factory.assert_not_called()
    else:
        assert result is None
        factory.assert_not_called()


def test_two_enabled_readers_are_rejected() -> None:
    with pytest.raises(ValueError, match="multiple"):
        subject.build_length_consumer(
            Mock(measurement_enabled=True),
            StandardLengthConfiguration(True, None, None, "a" * 64),
            assembly="hg38",
            reference=None,
            project_root="/run",
            samtools="samtools",
        )


def test_only_completed_standard_summary_can_be_published() -> None:
    config = subject.validate_pipeline_length(None, ())
    runner = StandardLengthRunner(StandardLengthConfiguration(True, None, None, "a" * 64), "hg38", None, Path("/run"))
    with pytest.raises(ValueError, match="did not complete"):
        subject.completed_length_summary(config, runner)
    runner.summary = {"length_estimation_status": "unavailable"}
    result = subject.completed_length_summary(config, runner)
    assert result == runner.summary
    assert result is not runner.summary


def test_completed_approved_summary_uses_existing_validation(monkeypatch: pytest.MonkeyPatch) -> None:
    config = subject.validate_pipeline_length(None, ())
    legacy_result = object()
    legacy = Mock(result=legacy_result)
    projector = Mock(return_value={"length_estimation_status": "measured-only"})
    assert subject.completed_length_summary(config, legacy, approved_projector=projector) == {
        "length_estimation_status": "measured-only"
    }
    projector.assert_called_once_with(config, legacy_result)
