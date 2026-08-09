"""Unit tests for the run-scoped CRAM reference environment contract."""

from __future__ import annotations

import os
import threading

import pytest

from vntyper.scripts.alignment_preflight import pin_reference_resolution, restore_reference_resolution

pytestmark = pytest.mark.unit


def test_an_unset_ref_path_is_pinned_to_the_shipped_local_only_default(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("REF_PATH", raising=False)

    previous = pin_reference_resolution({})

    assert previous is None
    assert os.environ["REF_PATH"] == "%2s/%2s/%s"
    restore_reference_resolution(previous)


def test_an_operators_network_ref_path_is_overridden_by_default(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("REF_PATH", "http://blackhole.invalid/%s")
    config = {"cram": {"local_ref_path": "/local/cache/%2s/%2s/%s"}}

    previous = pin_reference_resolution(config)

    assert previous == "http://blackhole.invalid/%s"
    assert os.environ["REF_PATH"] == "/local/cache/%2s/%2s/%s"
    restore_reference_resolution(previous)


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
    restore_reference_resolution(previous)


@pytest.mark.parametrize("invalid", ["false", "true", 0, 1, None])
def test_ambient_reference_resolution_requires_an_actual_boolean(
    monkeypatch: pytest.MonkeyPatch, invalid: object
) -> None:
    """JSON lookalikes cannot become a truthiness-based network waiver."""
    monkeypatch.setenv("REF_PATH", "/operator/original/%s")

    with pytest.raises(ValueError, match="allow_ambient_reference_resolution must be true or false"):
        pin_reference_resolution({"cram": {"allow_ambient_reference_resolution": invalid}})

    assert os.environ["REF_PATH"] == "/operator/original/%s"


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


def test_overlapping_threads_cannot_reopen_ambient_reference_resolution(monkeypatch: pytest.MonkeyPatch) -> None:
    """A second CRAM scope waits until the first has restored process-global REF_PATH."""
    ambient = "http://ambient.example/%s"
    monkeypatch.setenv("REF_PATH", ambient)
    first_pinned = threading.Event()
    release_first = threading.Event()
    second_attempting = threading.Event()
    second_pinned = threading.Event()
    failures: list[BaseException] = []
    observations: list[tuple[str | None, str]] = []

    def first_scope() -> None:
        try:
            previous = pin_reference_resolution({"cram": {"local_ref_path": "/cache/first/%s"}})
            first_pinned.set()
            assert release_first.wait(timeout=2)
            restore_reference_resolution(previous)
        except BaseException as error:
            failures.append(error)

    def second_scope() -> None:
        try:
            assert first_pinned.wait(timeout=2)
            second_attempting.set()
            previous = pin_reference_resolution({"cram": {"local_ref_path": "/cache/second/%s"}})
            observations.append((previous, os.environ["REF_PATH"]))
            second_pinned.set()
            restore_reference_resolution(previous)
        except BaseException as error:
            failures.append(error)

    first = threading.Thread(target=first_scope)
    second = threading.Thread(target=second_scope)
    first.start()
    second.start()
    try:
        assert first_pinned.wait(timeout=2)
        assert second_attempting.wait(timeout=2)
        assert not second_pinned.wait(timeout=0.1), "the second CRAM scope mutated REF_PATH while the first was active"
        assert os.environ["REF_PATH"] == "/cache/first/%s"
        release_first.set()
        assert second_pinned.wait(timeout=2)
    finally:
        release_first.set()
        first.join(timeout=2)
        second.join(timeout=2)

    assert not first.is_alive() and not second.is_alive()
    assert failures == []
    assert observations == [(ambient, "/cache/second/%s")]
    assert os.environ["REF_PATH"] == ambient


def test_an_invalid_scope_does_not_block_the_next_thread(monkeypatch: pytest.MonkeyPatch) -> None:
    """Validation failures release the process-global REF_PATH scope lock."""
    ambient = "/operator/original/%s"
    monkeypatch.setenv("REF_PATH", ambient)

    with pytest.raises(ValueError, match="allow_ambient_reference_resolution must be true or false"):
        pin_reference_resolution({"cram": {"allow_ambient_reference_resolution": "false"}})

    finished = threading.Event()
    failures: list[BaseException] = []

    def valid_scope() -> None:
        try:
            previous = pin_reference_resolution({"cram": {"local_ref_path": "/cache/valid/%s"}})
            assert os.environ["REF_PATH"] == "/cache/valid/%s"
            restore_reference_resolution(previous)
            finished.set()
        except BaseException as error:
            failures.append(error)

    worker = threading.Thread(target=valid_scope)
    worker.start()
    worker.join(timeout=2)

    assert finished.is_set()
    assert not worker.is_alive()
    assert failures == []
    assert os.environ["REF_PATH"] == ambient
