"""Ownership and transition regressions for run-local alignment indexes."""

from __future__ import annotations

import os
import shlex
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts.alignment_index_provenance import (
    generated_index_is_owned,
    remove_generated_index_provenance,
    write_generated_index_provenance,
)
from vntyper.scripts.alignment_preflight import build_alignment_view

pytestmark = pytest.mark.unit


def _config() -> dict:
    """Return the minimal shipped-default-compatible tool configuration."""
    return {"tools": {"samtools": "samtools"}}


def _build_index(
    command: str,
    log_file: str,
    cwd: str | None = None,
    *,
    protected_paths: tuple[str | Path, ...] = (),
) -> tuple[bool, str]:
    """Materialize the temporary index named by a mocked samtools command."""
    del log_file, cwd, protected_paths
    arguments = shlex.split(command)
    Path(arguments[arguments.index("-o") + 1]).write_bytes(b"generated index")
    return True, ""


def _alignment(directory: Path, suffix: str, contents: bytes) -> Path:
    """Create one unindexed alignment fixture."""
    directory.mkdir(parents=True)
    path = directory / f"sample.{suffix}"
    path.write_bytes(contents)
    return path


def _provenance(index_path: str | Path) -> Path:
    """Return the reserved provenance sidecar beside a generated index."""
    return Path(f"{index_path}.vntyper.json")


def test_provenance_public_api_round_trip_and_fail_closed_contract(tmp_path: Path) -> None:
    """The ownership API records, validates, replaces, and removes one artifact."""
    index = tmp_path / "sample.cram.crai"
    index.write_bytes(b"generated index")

    write_generated_index_provenance(index, (), replace_owned=False)

    assert generated_index_is_owned(index, ()) is True
    with pytest.raises(ValueError, match="collision"):
        write_generated_index_provenance(index, (), replace_owned=False)
    write_generated_index_provenance(index, (), replace_owned=True)
    assert generated_index_is_owned(index, ()) is True
    remove_generated_index_provenance(index, ())
    assert not _provenance(index).exists()
    with pytest.raises(ValueError, match="disappeared"):
        remove_generated_index_provenance(index, ())
    with pytest.raises(ValueError, match="no generated-index provenance"):
        generated_index_is_owned(index, ())


@pytest.mark.parametrize("index_kind", ["missing", "symlink", "directory"])
def test_provenance_write_rejects_a_non_regular_index(tmp_path: Path, index_kind: str) -> None:
    """Only a present non-symlink regular index can receive ownership."""
    index = tmp_path / "sample.cram.crai"
    if index_kind == "symlink":
        target = tmp_path / "target.crai"
        target.write_bytes(b"target")
        index.symlink_to(target)
    elif index_kind == "directory":
        index.mkdir()

    with pytest.raises(ValueError, match="[Gg]enerated index"):
        write_generated_index_provenance(index, (), replace_owned=False)


def test_provenance_replacement_requires_an_existing_owned_record(tmp_path: Path) -> None:
    """Replacement mode cannot silently create a missing ownership record."""
    index = tmp_path / "sample.cram.crai"
    index.write_bytes(b"generated index")

    with pytest.raises(ValueError, match="disappeared"):
        write_generated_index_provenance(index, (), replace_owned=True)


@pytest.mark.parametrize("index_kind", ["missing", "symlink", "directory"])
def test_unpaired_provenance_fails_closed(tmp_path: Path, index_kind: str) -> None:
    """Metadata without its matching regular index never establishes ownership."""
    index = tmp_path / "sample.cram.crai"
    if index_kind == "symlink":
        target = tmp_path / "target.crai"
        target.write_bytes(b"target")
        index.symlink_to(target)
    elif index_kind == "directory":
        index.mkdir()
    _provenance(index).write_text("{}")

    with pytest.raises(ValueError, match="provenance"):
        generated_index_is_owned(index, ())


def test_a_provenance_symlink_is_rejected_even_when_its_target_is_not_protected(tmp_path: Path) -> None:
    """No-follow validation rejects a metadata symlink before parsing its target."""
    index = tmp_path / "sample.cram.crai"
    index.write_bytes(b"generated index")
    unrelated = tmp_path / "unrelated.json"
    unrelated.write_text("{}")
    _provenance(index).symlink_to(unrelated)

    with pytest.raises(ValueError, match="must not be a symlink"):
        generated_index_is_owned(index, ())


def test_a_one_shot_protected_path_iterable_still_guards_the_provenance_inode(tmp_path: Path) -> None:
    """Multiple safety checks must not exhaust a caller's protected-path iterable."""
    index = tmp_path / "sample.cram.crai"
    index.write_bytes(b"generated index")
    patient = tmp_path / "patient.cram"
    patient.write_text("{}")
    os.link(patient, _provenance(index))

    with pytest.raises(ValueError, match="aliases protected source"):
        generated_index_is_owned(index, (path for path in (patient,)))

    assert patient.read_text() == "{}"


def test_unreadable_provenance_encoding_fails_closed(tmp_path: Path) -> None:
    """A no-follow read rejects bytes that are not valid JSON text."""
    index = tmp_path / "sample.cram.crai"
    index.write_bytes(b"generated index")
    _provenance(index).write_bytes(b"\xff")

    with pytest.raises(ValueError, match="Unable to read"):
        generated_index_is_owned(index, ())


def test_a_failed_exclusive_provenance_install_cleans_its_temporary_file(tmp_path: Path) -> None:
    """An OS failure cannot leave a misleading final record or temporary file."""
    index = tmp_path / "sample.cram.crai"
    index.write_bytes(b"generated index")

    with (
        patch("vntyper.scripts.alignment_index_provenance.os.link", side_effect=OSError("link failed")),
        pytest.raises(ValueError, match="link failed"),
    ):
        write_generated_index_provenance(index, (), replace_owned=False)

    assert not _provenance(index).exists()
    assert not tuple(tmp_path.glob("*.tmp"))


def test_a_provenance_remove_oserror_is_reported_without_hiding_the_record(tmp_path: Path) -> None:
    """A failed record removal is diagnostic and leaves the record intact."""
    index = tmp_path / "sample.cram.crai"
    index.write_bytes(b"generated index")
    write_generated_index_provenance(index, (), replace_owned=False)

    with (
        patch("vntyper.scripts.alignment_index_provenance.os.unlink", side_effect=OSError("unlink failed")),
        pytest.raises(ValueError, match="unlink failed"),
    ):
        remove_generated_index_provenance(index, ())

    assert _provenance(index).is_file()


def test_a_local_build_records_regular_output_contained_provenance(tmp_path: Path) -> None:
    """A generated index receives a non-symlink ownership record beside it."""
    cram = _alignment(tmp_path / "input", "cram", b"patient alignment")
    output = tmp_path / "output"

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index):
        _, generated_index, binding = build_alignment_view(
            str(cram), str(output), "sample", "cram", _config(), threads=2
        )

    provenance = _provenance(generated_index)
    assert provenance.parent == output
    assert provenance.is_file()
    assert not provenance.is_symlink()


@pytest.mark.parametrize("provenance_state", ["absent", "malformed", "directory"])
def test_an_unowned_regular_index_fails_closed_before_view_mutation(tmp_path: Path, provenance_state: str) -> None:
    """Absent or invalid provenance never authorizes replacement of a regular index."""
    cram = _alignment(tmp_path / "input", "cram", b"patient alignment")
    output = tmp_path / "output"
    output.mkdir()
    index = output / "sample.cram.crai"
    index.write_bytes(b"unowned index")
    provenance = _provenance(index)
    if provenance_state == "malformed":
        provenance.write_text("{")
    elif provenance_state == "directory":
        provenance.mkdir()

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index),
        pytest.raises(ValueError, match="provenance"),
    ):
        build_alignment_view(str(cram), str(output), "sample", "cram", _config(), threads=2)

    assert index.read_bytes() == b"unowned index"
    assert not (output / "sample.cram").exists()
    assert cram.read_bytes() == b"patient alignment"


def test_stale_generated_index_provenance_fails_closed(tmp_path: Path) -> None:
    """Changing a generated index invalidates its ownership fingerprint."""
    cram = _alignment(tmp_path / "input", "cram", b"patient alignment")
    output = tmp_path / "output"

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index):
        view_path, generated_index, binding = build_alignment_view(
            str(cram), str(output), "sample", "cram", _config(), threads=2
        )

    index = Path(generated_index)
    index.write_bytes(b"changed index with a different fingerprint")
    assert _provenance(index).is_file()

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index),
        pytest.raises(ValueError, match="provenance"),
    ):
        build_alignment_view(str(cram), str(output), "sample", "cram", _config(), threads=2)

    assert Path(view_path).samefile(cram)
    assert index.read_bytes() == b"changed index with a different fingerprint"
    assert cram.read_bytes() == b"patient alignment"


@pytest.mark.parametrize("alias_kind", ["symlink", "hardlink"])
def test_a_provenance_alias_of_patient_data_is_rejected_without_following_it(tmp_path: Path, alias_kind: str) -> None:
    """Ownership metadata cannot read through or unlink a patient-data alias."""
    cram = _alignment(tmp_path / "input", "cram", b"patient alignment")
    output = tmp_path / "output"
    output.mkdir()
    index = output / "sample.cram.crai"
    index.write_bytes(b"unowned index")
    provenance = _provenance(index)
    if alias_kind == "symlink":
        provenance.symlink_to(cram)
    else:
        os.link(cram, provenance)

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index),
        pytest.raises(ValueError, match="provenance|protected"),
    ):
        build_alignment_view(str(cram), str(output), "sample", "cram", _config(), threads=2)

    assert cram.read_bytes() == b"patient alignment"
    assert provenance.read_bytes() == b"patient alignment"
    assert index.read_bytes() == b"unowned index"
    assert not (output / "sample.cram").exists()


def test_generated_cram_index_is_rebuilt_when_a_source_crai_appears(tmp_path: Path) -> None:
    """A newly appearing source CRAI is preserved but never trusted."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    cram = _alignment(input_dir, "cram", b"patient alignment")

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index):
        _, generated_index, binding = build_alignment_view(
            str(cram), str(output), "sample", "cram", _config(), threads=2
        )

    assert Path(generated_index).is_file()
    assert not Path(generated_index).is_symlink()
    source_crai = Path(f"{cram}.crai")
    source_crai.write_bytes(b"patient source index")

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index):
        view_path, selected_index, binding = build_alignment_view(
            str(cram), str(output), "sample", "cram", _config(), threads=2
        )

    assert Path(view_path).samefile(cram)
    assert Path(selected_index).is_file()
    assert not Path(selected_index).is_symlink()
    assert generated_index_is_owned(selected_index, (cram, source_crai)) is True
    assert cram.read_bytes() == b"patient alignment"
    assert source_crai.read_bytes() == b"patient source index"


def test_generated_cram_index_is_rebuilt_for_a_different_input_with_source_csi(tmp_path: Path) -> None:
    """Output reuse rebuilds from the new CRAM instead of trusting its CSI."""
    output = tmp_path / "output"
    first_cram = _alignment(tmp_path / "first", "cram", b"first patient alignment")
    second_cram = _alignment(tmp_path / "second", "cram", b"second patient alignment")

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index):
        _, stale_crai, binding = build_alignment_view(
            str(first_cram), str(output), "sample", "cram", _config(), threads=2
        )

    stale_provenance = _provenance(stale_crai)
    source_csi = Path(f"{second_cram}.csi")
    source_csi.write_bytes(b"second patient source index")

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index):
        view_path, selected_index, binding = build_alignment_view(
            str(second_cram), str(output), "sample", "cram", _config(), threads=2
        )

    assert Path(view_path).samefile(second_cram)
    assert Path(selected_index) == output / "sample.cram.crai"
    assert Path(selected_index).is_file()
    assert not Path(selected_index).is_symlink()
    assert generated_index_is_owned(selected_index, (second_cram, source_csi)) is True
    assert stale_provenance.is_file()
    assert first_cram.read_bytes() == b"first patient alignment"
    assert second_cram.read_bytes() == b"second patient alignment"
    assert source_csi.read_bytes() == b"second patient source index"


@pytest.mark.parametrize(
    ("bai_only", "source_suffix"),
    [(True, ".bai"), (False, ".csi")],
    ids=["non-fast-bai", "fast-csi"],
)
def test_generated_bam_index_is_rebuilt_even_when_an_accepted_source_index_appears(
    tmp_path: Path, bai_only: bool, source_suffix: str
) -> None:
    """Neither BAI-only nor fast mode trusts an unbound source index."""
    input_dir = tmp_path / "input"
    output = tmp_path / "output"
    bam = _alignment(input_dir, "bam", b"patient alignment")

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index):
        _, generated_index, binding = build_alignment_view(
            str(bam), str(output), "sample", "bam", _config(), threads=2, bai_only=bai_only
        )

    generated_provenance = _provenance(generated_index)
    source_index = Path(f"{bam}{source_suffix}")
    source_index.write_bytes(b"patient source index")

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index):
        view_path, selected_index, binding = build_alignment_view(
            str(bam), str(output), "sample", "bam", _config(), threads=2, bai_only=bai_only
        )

    assert Path(view_path).samefile(bam)
    assert Path(selected_index).is_file()
    assert not Path(selected_index).is_symlink()
    assert generated_index_is_owned(selected_index, (bam, source_index)) is True
    assert generated_provenance.is_file()
    assert bam.read_bytes() == b"patient alignment"
    assert source_index.read_bytes() == b"patient source index"


def test_nonfast_bam_rerun_ignores_new_source_csi_and_rebuilds_owned_bai(tmp_path: Path) -> None:
    """The BAI-only consumer keeps a generated BAI when only source CSI appears."""
    bam = _alignment(tmp_path / "input", "bam", b"patient alignment")
    output = tmp_path / "output"

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index):
        _, first_index, binding = build_alignment_view(
            str(bam), str(output), "sample", "bam", _config(), threads=2, bai_only=True
        )
        source_csi = Path(f"{bam}.csi")
        source_csi.write_bytes(b"patient source CSI")
        view_path, second_index, binding = build_alignment_view(
            str(bam), str(output), "sample", "bam", _config(), threads=2, bai_only=True
        )

    assert first_index == second_index
    assert Path(view_path).samefile(bam)
    assert Path(second_index).is_file()
    assert not Path(second_index).is_symlink()
    assert generated_index_is_owned(second_index, (bam,)) is True
    assert source_csi.read_bytes() == b"patient source CSI"
    assert not Path(f"{bam}.bai").exists()


def test_first_provenance_install_cleanup_failure_leaves_no_final_pair_and_rerun_recovers(
    tmp_path: Path,
) -> None:
    """A failed first ownership install cannot strand provenance without its index."""
    cram = _alignment(tmp_path / "input", "cram", b"patient alignment")
    output = tmp_path / "output"
    view_index = output / "sample.cram.crai"
    provenance = _provenance(view_index)
    real_unlink = os.unlink
    fault_injected = False

    def fail_first_temporary_provenance_unlink(path: str | os.PathLike[str]) -> None:
        nonlocal fault_injected
        candidate = Path(path)
        if (
            not fault_injected
            and candidate.parent == output
            and candidate.name.startswith(f".{provenance.name}.")
            and candidate.name.endswith(".tmp")
            and provenance.exists()
        ):
            fault_injected = True
            raise OSError("temporary provenance unlink failed")
        real_unlink(path)

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index),
        patch(
            "vntyper.scripts.alignment_index_provenance.os.unlink",
            side_effect=fail_first_temporary_provenance_unlink,
        ),
        pytest.raises(ValueError, match="temporary provenance unlink failed"),
    ):
        build_alignment_view(str(cram), str(output), "sample", "cram", _config(), threads=2)

    assert fault_injected
    assert not os.path.lexists(view_index)
    assert not os.path.lexists(provenance)

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index):
        _, recovered_index, binding = build_alignment_view(
            str(cram), str(output), "sample", "cram", _config(), threads=2
        )

    assert Path(recovered_index) == view_index
    assert generated_index_is_owned(recovered_index, (cram,)) is True


def test_provenance_update_failure_restores_the_old_pair_and_rerun_recovers(tmp_path: Path) -> None:
    """Rebuilding an owned index rolls its bytes back when metadata cannot update."""
    cram = _alignment(tmp_path / "input", "cram", b"patient alignment")
    output = tmp_path / "output"
    build_number = 0

    def build_numbered_index(
        command: str,
        log_file: str,
        cwd: str | None = None,
        *,
        protected_paths: tuple[str | Path, ...] = (),
    ) -> tuple[bool, str]:
        nonlocal build_number
        del log_file, cwd, protected_paths
        build_number += 1
        arguments = shlex.split(command)
        Path(arguments[arguments.index("-o") + 1]).write_bytes(f"generated index {build_number}".encode())
        return True, ""

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=build_numbered_index):
        _, generated_index, binding = build_alignment_view(
            str(cram), str(output), "sample", "cram", _config(), threads=2
        )

    index = Path(generated_index)
    provenance = _provenance(index)
    old_inode = index.stat().st_ino
    real_replace = os.replace
    fault_injected = False

    def fail_provenance_update(
        source: str | os.PathLike[str],
        destination: str | os.PathLike[str],
    ) -> None:
        nonlocal fault_injected
        if not fault_injected and Path(destination) == provenance:
            fault_injected = True
            raise OSError("provenance update failed")
        real_replace(source, destination)

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=build_numbered_index),
        patch("vntyper.scripts.alignment_index_provenance.os.replace", side_effect=fail_provenance_update),
        pytest.raises(ValueError, match="provenance update failed"),
    ):
        build_alignment_view(str(cram), str(output), "sample", "cram", _config(), threads=2)

    assert fault_injected
    assert index.read_bytes() == b"generated index 1"
    assert index.stat().st_ino == old_inode
    assert generated_index_is_owned(index, (cram,)) is True

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=build_numbered_index):
        _, recovered_index, binding = build_alignment_view(
            str(cram), str(output), "sample", "cram", _config(), threads=2
        )

    assert recovered_index == generated_index
    assert index.read_bytes() == b"generated index 3"
    assert generated_index_is_owned(index, (cram,)) is True


def test_rebuild_tombstone_cleanup_failure_keeps_the_new_owned_index(tmp_path: Path) -> None:
    """A post-commit tombstone cleanup failure cannot roll back a trusted rebuild."""
    first_cram = _alignment(tmp_path / "first", "cram", b"first patient alignment")
    second_cram = _alignment(tmp_path / "second", "cram", b"second patient alignment")
    output = tmp_path / "output"

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index):
        _, generated_index, binding = build_alignment_view(
            str(first_cram), str(output), "sample", "cram", _config(), threads=2
        )

    stale_index = Path(generated_index)
    stale_provenance = _provenance(stale_index)
    source_csi = Path(f"{second_cram}.csi")
    source_csi.write_bytes(b"patient source index")
    real_unlink = os.unlink
    fault_injected = False

    def fail_stale_index_unlink(path: str | os.PathLike[str]) -> None:
        nonlocal fault_injected
        candidate = Path(path)
        if (
            not fault_injected
            and candidate.is_file()
            and not candidate.is_symlink()
            and candidate.read_bytes() == b"generated index"
        ):
            fault_injected = True
            raise OSError("stale index unlink failed")
        real_unlink(path)

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=_build_index),
        patch("vntyper.scripts.alignment_index_provenance.os.unlink", side_effect=fail_stale_index_unlink),
    ):
        view_path, rebuilt_index, binding = build_alignment_view(
            str(second_cram), str(output), "sample", "cram", _config(), threads=2
        )

    assert fault_injected
    assert Path(view_path).samefile(second_cram)
    assert Path(rebuilt_index) == output / "sample.cram.crai"
    assert Path(rebuilt_index).is_file()
    assert not Path(rebuilt_index).is_symlink()
    assert stale_provenance.is_file()
    assert generated_index_is_owned(rebuilt_index, (second_cram, source_csi)) is True
    assert source_csi.read_bytes() == b"patient source index"
