"""Pure adVNTR preflight planning before any pipeline stage can run."""

from __future__ import annotations

import logging
from unittest.mock import patch

import pytest

from vntyper.scripts.pipeline_advntr_preflight import AdvntrPreflight, plan_advntr_preflight

pytestmark = pytest.mark.unit


CONFIG = {
    "reference_data": {
        "advntr_reference_vntr_hg19": "/refs/hg19.db",
        "advntr_reference_vntr_hg38": "/refs/hg38.db",
    }
}


def test_disabled_module_needs_no_reference_configuration() -> None:
    assert plan_advntr_preflight({}, [], {}, "hg38") == AdvntrPreflight(enabled=False, reference=None)


@pytest.mark.parametrize(
    ("reference_assembly", "expected"),
    [("hg19", "/refs/hg19.db"), ("hg38", "/refs/hg38.db")],
)
def test_unconfigured_override_uses_the_requested_assembly(reference_assembly: str, expected: str) -> None:
    assert plan_advntr_preflight(CONFIG, ["advntr"], {}, reference_assembly) == AdvntrPreflight(
        enabled=True,
        reference=expected,
    )


@pytest.mark.parametrize(("override", "expected"), [("hg19", "/refs/hg19.db"), ("hg38", "/refs/hg38.db")])
def test_named_override_resolves_to_its_explicit_database(override: str, expected: str) -> None:
    module_args = {"advntr": {"advntr_reference": override}}
    assert plan_advntr_preflight(CONFIG, ["advntr"], module_args, "hg19") == AdvntrPreflight(
        enabled=True,
        reference=expected,
    )


def test_literal_path_override_is_rejected_before_processing(caplog: pytest.LogCaptureFixture) -> None:
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.pipeline_advntr_preflight")

    with pytest.raises(ValueError, match="Invalid advntr_reference: /tmp/model.db"):
        plan_advntr_preflight(
            CONFIG,
            ["advntr"],
            {"advntr": {"advntr_reference": "/tmp/model.db"}},
            "hg19",
        )

    assert "Invalid advntr_reference: /tmp/model.db" in caplog.messages


@pytest.mark.parametrize("override", [["hg19"], {"build": "hg19"}, {"hg19"}, ("hg19",), True, 1])
def test_non_string_override_logs_a_keyed_value_error_before_reference_resolution(
    override: object, caplog: pytest.LogCaptureFixture
) -> None:
    """Unhashable config values must not escape the planner's explicit failure contract."""
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.pipeline_advntr_preflight")

    with (
        patch("vntyper.scripts.pipeline_advntr_preflight.select_advntr_reference") as resolver,
        pytest.raises(ValueError, match="Invalid advntr_reference"),
    ):
        plan_advntr_preflight(
            CONFIG,
            ["advntr"],
            {"advntr": {"advntr_reference": override}},
            "hg19",
        )

    resolver.assert_not_called()
    planner_messages = [
        record.getMessage()
        for record in caplog.records
        if record.name == "vntyper.scripts.pipeline_advntr_preflight" and record.levelno == logging.ERROR
    ]
    assert len(planner_messages) == 1
    assert planner_messages[0].startswith("Invalid advntr_reference:")


@pytest.mark.parametrize("override", [[], {}, set(), (), False, 0])
def test_falsy_non_string_override_is_not_mistaken_for_absence(override: object) -> None:
    """An explicitly configured non-string is invalid even when it is falsy."""
    with (
        patch("vntyper.scripts.pipeline_advntr_preflight.select_advntr_reference") as resolver,
        pytest.raises(ValueError, match="Invalid advntr_reference"),
    ):
        plan_advntr_preflight(
            CONFIG,
            ["advntr"],
            {"advntr": {"advntr_reference": override}},
            "hg19",
        )

    resolver.assert_not_called()


@pytest.mark.parametrize(
    ("config", "module_args", "assembly"),
    [
        ({}, {}, "hg19"),
        ({"reference_data": {"advntr_reference_vntr_hg19": None}}, {}, "hg19"),
        ({"reference_data": {}}, {"advntr": {"advntr_reference": "hg38"}}, "hg19"),
    ],
)
def test_missing_resolved_path_fails_loudly(
    config: dict[str, object],
    module_args: dict[str, dict[str, str]],
    assembly: str,
    caplog: pytest.LogCaptureFixture,
) -> None:
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.pipeline_advntr_preflight")

    with pytest.raises(ValueError, match="adVNTR reference path not found in configuration"):
        plan_advntr_preflight(config, ["advntr"], module_args, assembly)

    assert "adVNTR reference path not found in configuration." in caplog.messages


def test_non_string_resolved_path_fails_before_filesystem_io(caplog: pytest.LogCaptureFixture) -> None:
    """Stringifying malformed configuration must not create a surprising path."""
    config = {"reference_data": {"advntr_reference_vntr_hg19": ["/refs/hg19.db"]}}
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.pipeline_advntr_preflight")

    with pytest.raises(ValueError, match="adVNTR reference path must be a non-empty string"):
        plan_advntr_preflight(config, ["advntr"], {}, "hg19")

    assert "adVNTR reference path must be a non-empty string." in caplog.messages
