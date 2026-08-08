"""Unit tests for the run-scoped CRAM reference environment contract."""

from __future__ import annotations

import os

import pytest

from vntyper.scripts.alignment_preflight import pin_reference_resolution, restore_reference_resolution

pytestmark = pytest.mark.unit


def test_an_unset_ref_path_is_pinned_to_the_shipped_local_only_default(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("REF_PATH", raising=False)

    previous = pin_reference_resolution({})

    assert previous is None
    assert os.environ["REF_PATH"] == "%2s/%2s/%s"


def test_an_operators_network_ref_path_is_overridden_by_default(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("REF_PATH", "http://blackhole.invalid/%s")
    config = {"cram": {"local_ref_path": "/local/cache/%2s/%2s/%s"}}

    previous = pin_reference_resolution(config)

    assert previous == "http://blackhole.invalid/%s"
    assert os.environ["REF_PATH"] == "/local/cache/%2s/%2s/%s"


def test_the_override_is_skipped_when_ambient_resolution_is_explicitly_allowed(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    ambient = "https://refget.example/%s"
    monkeypatch.setenv("REF_PATH", ambient)
    config = {"cram": {"allow_ambient_reference_resolution": True}}

    with caplog.at_level("WARNING"):
        previous = pin_reference_resolution(config)

    assert previous == ambient
    assert os.environ["REF_PATH"] == ambient
    assert "block" in caplog.text.lower()


@pytest.mark.parametrize(
    "remote_path",
    [
        "http://refget.example/%s",
        "https://refget.example/%s",
        "/local/cache/%s:http://refget.example/%s",
        "s3://reference-bucket/%s",
    ],
)
def test_a_remote_local_ref_path_is_rejected_when_ambient_resolution_is_disabled(
    monkeypatch: pytest.MonkeyPatch, remote_path: str
) -> None:
    original = "/operator/original/%s"
    monkeypatch.setenv("REF_PATH", original)
    config = {"cram": {"allow_ambient_reference_resolution": False, "local_ref_path": remote_path}}

    with pytest.raises(ValueError, match="local_ref_path must not contain a remote URL"):
        pin_reference_resolution(config)

    assert os.environ["REF_PATH"] == original


@pytest.mark.parametrize("previous", [None, "/original/%s"])
def test_restore_puts_back_the_previous_value_including_unset(
    monkeypatch: pytest.MonkeyPatch, previous: str | None
) -> None:
    monkeypatch.setenv("REF_PATH", "/pinned/%s")

    restore_reference_resolution(previous)

    if previous is None:
        assert "REF_PATH" not in os.environ
    else:
        assert os.environ["REF_PATH"] == previous
