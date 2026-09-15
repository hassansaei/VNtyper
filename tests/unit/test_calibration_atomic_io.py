"""No-clobber atomic installation for calibration artifact directories."""

import stat
from importlib import import_module
from pathlib import Path
from unittest.mock import patch

import pytest

pytestmark = pytest.mark.unit


def test_atomic_output_installs_one_private_complete_tree(tmp_path: Path) -> None:
    atomic_output = import_module("vntyper.scripts.calibration_atomic_io").atomic_output
    output = tmp_path / "result"

    def produce(staging: Path) -> bool:
        assert staging.parent == output.parent
        assert stat.S_IMODE(staging.stat().st_mode) == 0o700
        (staging / "complete.json").write_bytes(b"{}\n")
        return False

    assert atomic_output(output, produce) is False
    assert (output / "complete.json").read_bytes() == b"{}\n"
    assert stat.S_IMODE(output.stat().st_mode) == 0o700
    assert not tuple(tmp_path.glob(".result.*"))


@pytest.mark.parametrize("kind", ["file", "directory", "symlink"])
def test_atomic_output_rejects_an_existing_destination_before_production(tmp_path: Path, kind: str) -> None:
    atomic_output = import_module("vntyper.scripts.calibration_atomic_io").atomic_output
    output = tmp_path / "result"
    if kind == "file":
        output.write_text("occupied", encoding="utf-8")
    elif kind == "directory":
        output.mkdir()
    else:
        target = tmp_path / "target"
        target.mkdir()
        output.symlink_to(target, target_is_directory=True)
    called = False

    def produce(_staging: Path) -> bool:
        nonlocal called
        called = True
        return True

    with pytest.raises(ValueError, match="already exists"):
        atomic_output(output, produce)
    assert called is False


def test_atomic_output_refuses_a_racing_empty_directory_without_deleting_it(tmp_path: Path) -> None:
    atomic_output = import_module("vntyper.scripts.calibration_atomic_io").atomic_output
    output = tmp_path / "result"

    def race(staging: Path) -> bool:
        (staging / "complete.json").write_bytes(b"{}\n")
        output.mkdir()
        return True

    with pytest.raises(ValueError, match="already exists"):
        atomic_output(output, race)

    assert output.is_dir()
    assert not tuple(output.iterdir())
    assert not tuple(tmp_path.glob(".result.*"))


@pytest.mark.parametrize("failure", [KeyboardInterrupt(), SystemExit(7)])
def test_atomic_output_cleans_its_staging_tree_after_base_exception(tmp_path: Path, failure: BaseException) -> None:
    atomic_output = import_module("vntyper.scripts.calibration_atomic_io").atomic_output
    output = tmp_path / "result"

    def interrupt(staging: Path) -> bool:
        (staging / "partial.json").write_text("partial", encoding="utf-8")
        raise failure

    with pytest.raises(type(failure)):
        atomic_output(output, interrupt)
    assert not output.exists()
    assert not tuple(tmp_path.glob(".result.*"))


@pytest.mark.parametrize("mode", ["empty", "non_boolean"])
def test_atomic_output_rejects_incomplete_producer_contract(tmp_path: Path, mode: str) -> None:
    atomic_output = import_module("vntyper.scripts.calibration_atomic_io").atomic_output
    output = tmp_path / "result"

    def produce(staging: Path):
        if mode != "empty":
            (staging / "result.json").write_bytes(b"{}\n")
        return True if mode == "empty" else 1

    with pytest.raises(ValueError, match="no artifacts|success value"):
        atomic_output(output, produce)
    assert not output.exists()
    assert not tuple(tmp_path.glob(".result.*"))


def test_atomic_output_fails_closed_when_renameat2_is_unavailable(tmp_path: Path) -> None:
    module = import_module("vntyper.scripts.calibration_atomic_io")
    output = tmp_path / "result"

    def produce(staging: Path) -> bool:
        (staging / "result.json").write_bytes(b"{}\n")
        return True

    with patch.object(module, "_renameat2", None), pytest.raises(RuntimeError, match="renameat2"):
        module.atomic_output(output, produce)
    assert not output.exists()
    assert not tuple(tmp_path.glob(".result.*"))
