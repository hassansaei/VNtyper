"""Atomic publication of derived adVNTR result tables."""

from pathlib import Path

import pandas as pd
import pytest

pytestmark = pytest.mark.unit


def _publisher():
    from vntyper.modules.advntr.advntr_result_io import publish_advntr_result

    return publish_advntr_result


def test_publish_advntr_result_replaces_the_destination_with_exact_tsv_bytes(tmp_path: Path) -> None:
    """A successful publication exposes the complete established TSV contract."""
    destination = tmp_path / "output_adVNTR_result.tsv"
    destination.write_text("stale prior result\n", encoding="utf-8")
    frame = pd.DataFrame([{"VID": "25561", "Variant": "I22_2_G_LEN1"}])

    _publisher()(frame, destination)

    assert destination.read_bytes() == b"VID\tVariant\n25561\tI22_2_G_LEN1\n"
    assert list(tmp_path.iterdir()) == [destination]


def test_publish_advntr_result_write_failure_preserves_destination_and_removes_temporary_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Partial candidate bytes never replace an already complete published table."""
    destination = tmp_path / "output_adVNTR_result.tsv"
    destination.write_bytes(b"complete prior result\n")
    frame = pd.DataFrame([{"VID": "25561"}])

    def fail_after_partial_write(self, path_or_buf, **_kwargs):
        path_or_buf.write("partial candidate")
        raise OSError("disk full")

    monkeypatch.setattr(pd.DataFrame, "to_csv", fail_after_partial_write)
    expected = f"Failed to publish adVNTR result {destination}: disk full"

    with caplog.at_level("ERROR"), pytest.raises(RuntimeError) as raised:
        _publisher()(frame, destination)

    assert str(raised.value) == expected
    assert expected in caplog.messages
    assert destination.read_bytes() == b"complete prior result\n"
    assert list(tmp_path.iterdir()) == [destination]


def test_publish_advntr_result_replace_failure_preserves_destination_and_removes_temporary_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A failed atomic landing leaves neither candidate bytes nor a damaged result."""
    from vntyper.modules.advntr import advntr_result_io

    destination = tmp_path / "output_adVNTR_result.tsv"
    destination.write_bytes(b"complete prior result\n")
    frame = pd.DataFrame([{"VID": "25561"}])
    monkeypatch.setattr(advntr_result_io.os, "replace", lambda *_args: (_ for _ in ()).throw(OSError("blocked")))
    expected = f"Failed to publish adVNTR result {destination}: blocked"

    with caplog.at_level("ERROR"), pytest.raises(RuntimeError) as raised:
        advntr_result_io.publish_advntr_result(frame, destination)

    assert str(raised.value) == expected
    assert expected in caplog.messages
    assert destination.read_bytes() == b"complete prior result\n"
    assert list(tmp_path.iterdir()) == [destination]


def test_publish_advntr_result_logs_cleanup_failure_without_replacing_primary_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A leaked candidate is observable while the publication failure stays authoritative."""
    from vntyper.modules.advntr import advntr_result_io

    destination = tmp_path / "output_adVNTR_result.tsv"
    destination.write_bytes(b"complete prior result\n")
    frame = pd.DataFrame([{"VID": "25561"}])
    original_unlink = Path.unlink

    def fail_after_partial_write(self, path_or_buf, **_kwargs):
        path_or_buf.write("partial candidate")
        raise OSError("disk full")

    def refuse_candidate_cleanup(path: Path, missing_ok: bool = False) -> None:
        if path.name.startswith(f".{destination.name}.") and path.suffix == ".tmp":
            raise OSError("cleanup denied")
        original_unlink(path, missing_ok=missing_ok)

    monkeypatch.setattr(pd.DataFrame, "to_csv", fail_after_partial_write)
    monkeypatch.setattr(Path, "unlink", refuse_candidate_cleanup)
    primary_message = f"Failed to publish adVNTR result {destination}: disk full"

    with caplog.at_level("ERROR", logger=advntr_result_io.logger.name), pytest.raises(RuntimeError) as raised:
        advntr_result_io.publish_advntr_result(frame, destination)

    assert str(raised.value) == primary_message
    assert primary_message in caplog.messages
    assert any(
        message.startswith("Failed to remove incomplete adVNTR result candidate ")
        and message.endswith(": cleanup denied")
        for message in caplog.messages
    )
    assert destination.read_bytes() == b"complete prior result\n"
    candidates = [path for path in tmp_path.iterdir() if path != destination]
    assert len(candidates) == 1
    assert candidates[0].read_text(encoding="utf-8") == "partial candidate"
