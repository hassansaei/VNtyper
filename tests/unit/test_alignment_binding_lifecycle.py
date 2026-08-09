"""Ownership, cleanup, and portability tests for alignment descriptor bindings."""

from __future__ import annotations

import errno
import gc
import os
import shlex
from collections.abc import Iterable
from contextlib import suppress
from pathlib import Path
from unittest.mock import patch

import pytest

from vntyper.scripts import alignment_preflight, fastq_bam_processing
from vntyper.scripts.alignment_binding import AlignmentBinding
from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.alignment_index_binding import AlignmentIndexBinding
from vntyper.scripts.alignment_preflight import build_alignment_view
from vntyper.scripts.alignment_target_io import protect_alignment_inputs

pytestmark = pytest.mark.unit


def _build_view(
    tmp_path: Path,
    *,
    input_symlink: bool = False,
    file_format: str = "bam",
) -> tuple[Path, Path, AlignmentPlan]:
    """Build one mocked-index BAM view through the production binding seam."""
    input_dir = tmp_path / "patient-input"
    input_dir.mkdir(exist_ok=True)
    alignment = input_dir / f"sample.{file_format}"
    alignment.write_bytes(b"original alignment")
    input_path = alignment
    if input_symlink:
        input_path = tmp_path / "sample-link.bam"
        input_path.symlink_to(alignment)
    stage_dir = tmp_path / "run" / "fastq_bam_processing"

    def build_index(command: str, *_args: object, **_kwargs: object) -> tuple[bool, str]:
        arguments = shlex.split(command)
        Path(arguments[arguments.index("-o") + 1]).write_bytes(b"fresh index")
        return True, ""

    with patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=build_index):
        view_path, index_path, binding = build_alignment_view(
            str(input_path),
            str(stage_dir),
            "input",
            file_format,
            {"tools": {"samtools": "samtools"}},
            threads=1,
        )

    plan = AlignmentPlan(
        input_path=str(input_path),
        view_path=view_path,
        file_format=file_format,
        index_path=index_path,
        reference_path=None,
        reference_source="not-required",
        uncovered_contigs=(),
        unmapped_scan="indexed",
        binding=binding,
    )
    return alignment, stage_dir, plan


def _hide_proc_descriptor(
    path: str | os.PathLike[str],
    *,
    dir_fd: int | None = None,
    follow_symlinks: bool = True,
) -> os.stat_result:
    """Make only proc descriptor paths unavailable while preserving normal stat calls."""
    if os.fspath(path).startswith("/proc/"):
        raise FileNotFoundError("procfs unavailable")
    return _REAL_STAT(path, dir_fd=dir_fd, follow_symlinks=follow_symlinks)


_REAL_STAT = os.stat


def test_plan_close_removes_the_exact_owned_view_before_closing_its_descriptor(tmp_path: Path) -> None:
    """Index and alignment views disappear before their respective descriptors close."""
    alignment, _, plan = _build_view(tmp_path)
    view = Path(plan.view_path)
    index = Path(plan.index_path)
    close_observations: list[bool] = []
    real_close = os.close

    def observe_close(descriptor: int) -> None:
        close_observations.append(not os.path.lexists(view))
        real_close(descriptor)

    with patch("vntyper.scripts.alignment_binding.os.close", side_effect=observe_close):
        plan.close()

    assert close_observations == [False, True]
    assert not os.path.lexists(view)
    assert index.read_bytes() == b"fresh index"
    assert alignment.read_bytes() == b"original alignment"


def test_conversion_uses_the_bound_fresh_index_after_the_public_index_is_replaced(tmp_path: Path) -> None:
    """Removing `-X` must reproduce the accepted wrong-index empty-slice regression."""
    _, _, plan = _build_view(tmp_path)
    public_index = Path(plan.index_path)
    replacement = public_index.with_name("wrong-sample.bai")
    replacement.write_bytes(b"valid index from another sample")
    replacement.replace(public_index)
    bed = tmp_path / "target.bed"
    bed.write_text("chr1\t0\t10\n", encoding="utf-8")
    commands: list[str] = []

    def record(command: str, *_args: object, **_kwargs: object) -> bool:
        commands.append(command)
        return True

    try:
        with patch.object(fastq_bam_processing, "run_command", side_effect=record):
            fastq_bam_processing.process_bam_to_fastq(
                output=tmp_path / "conversion",
                output_name="output",
                threads=1,
                config={"tools": {"samtools": "samtools"}},
                plan=plan,
                fast_mode=True,
                bed_file=bed,
            )

        slice_tokens = shlex.split(commands[0].split("&&", maxsplit=1)[0])
        custom_index_position = slice_tokens.index("-X")
        assert slice_tokens[custom_index_position + 1] == plan.view_path
        bound_index = Path(slice_tokens[custom_index_position + 2])
        assert bound_index != public_index
        assert bound_index.read_bytes() == b"fresh index"
        assert public_index.read_bytes() == b"valid index from another sample"
    finally:
        plan.close()


@pytest.mark.parametrize("file_format", ["bam", "cram"])
def test_proc_index_binding_retains_fresh_bai_and_crai_bytes(tmp_path: Path, file_format: str) -> None:
    """Both generated index formats remain descriptor-bound after publication."""
    _, _, plan = _build_view(tmp_path, file_format=file_format)
    replacement = Path(plan.index_path).with_name(f"replacement.{file_format}.index")
    replacement.write_bytes(b"other alignment index")
    replacement.replace(plan.index_path)

    try:
        assert plan.stable_index_path.startswith("/proc/")
        assert Path(plan.stable_index_path).read_bytes() == b"fresh index"
        assert Path(plan.index_path).read_bytes() == b"other alignment index"
    finally:
        plan.close()


def test_index_descriptor_is_open_before_the_public_index_can_be_replaced(tmp_path: Path) -> None:
    """Binding after publication would retain replacement bytes instead of generated bytes."""
    real_install = alignment_preflight._install_generated_index

    def install_then_replace(
        temporary_index: str | Path,
        index_path: str | Path,
        protected_paths: Iterable[str | Path],
        *,
        replace_owned: bool,
    ) -> None:
        real_install(temporary_index, index_path, protected_paths, replace_owned=replace_owned)
        published = Path(index_path)
        replacement = published.with_name("replacement-after-publication.bai")
        replacement.write_bytes(b"replacement after publication")
        replacement.replace(published)

    with patch.object(alignment_preflight, "_install_generated_index", side_effect=install_then_replace):
        _, _, plan = _build_view(tmp_path)

    try:
        assert Path(plan.stable_index_path).read_bytes() == b"fresh index"
        assert Path(plan.index_path).read_bytes() == b"replacement after publication"
    finally:
        plan.close()


def test_alignment_binding_refuses_to_replace_an_owned_index_descriptor(tmp_path: Path) -> None:
    """A second bind cannot silently discard the exact generated index."""
    _, _, plan = _build_view(tmp_path)
    replacement = tmp_path / "second.bai"
    replacement.write_bytes(b"second index")
    binding = plan.binding
    assert binding is not None

    try:
        with pytest.raises(RuntimeError, match="already owns a generated-index descriptor"):
            binding.bind_index(replacement, tmp_path / ".second.bai.bound")
        assert Path(plan.stable_index_path).read_bytes() == b"fresh index"
    finally:
        plan.close()


def test_index_proc_unavailable_uses_and_cleans_a_verified_hardlink(tmp_path: Path) -> None:
    """The portable fallback keeps the original bytes and leaves no archiveable view."""
    with patch("vntyper.scripts.alignment_index_binding.os.stat", side_effect=_hide_proc_descriptor):
        _, _, plan = _build_view(tmp_path)
    bound_index = Path(plan.stable_index_path)
    replacement = Path(plan.index_path).with_name("replacement.bai")
    replacement.write_bytes(b"other alignment index")
    replacement.replace(plan.index_path)

    assert bound_index.read_bytes() == b"fresh index"
    assert bound_index != Path(plan.index_path)
    plan.close()
    assert not os.path.lexists(bound_index)
    assert Path(plan.index_path).read_bytes() == b"other alignment index"


def test_index_proc_identity_mismatch_uses_the_verified_hardlink(tmp_path: Path) -> None:
    """A descriptor filesystem with the wrong identity is treated as unavailable."""
    index = tmp_path / "fresh.bai"
    fallback = tmp_path / ".fresh.bai.bound"
    other = tmp_path / "other.bai"
    index.write_bytes(b"fresh index")
    other.write_bytes(b"other index")
    wrong_stat = other.stat()

    with patch("vntyper.scripts.alignment_index_binding.os.stat", return_value=wrong_stat):
        binding = AlignmentIndexBinding(index, fallback)

    assert binding.consumer_path == str(fallback)
    assert fallback.read_bytes() == b"fresh index"
    binding.close()
    assert not fallback.exists()


def test_index_fallback_link_collision_never_removes_the_colliding_entry(tmp_path: Path) -> None:
    """A link-creation race must not unlink an entry the binding did not create."""
    index = tmp_path / "fresh.bai"
    fallback = tmp_path / ".fresh.bai.bound"
    index.write_bytes(b"fresh index")

    def collide_then_fail(*_args: object, **_kwargs: object) -> None:
        fallback.write_bytes(b"other entry")
        raise FileExistsError("binding collision")

    with (
        patch("vntyper.scripts.alignment_index_binding.os.stat", side_effect=_hide_proc_descriptor),
        patch("vntyper.scripts.alignment_index_binding.os.link", side_effect=collide_then_fail),
        pytest.raises(RuntimeError, match="binding collision"),
    ):
        AlignmentIndexBinding(index, fallback)

    assert fallback.read_bytes() == b"other entry"


def test_replaced_index_fallback_is_preserved_and_keeps_both_bindings_open(tmp_path: Path) -> None:
    """Cleanup refusal preserves replacement bytes and prevents descriptor reuse."""
    with patch("vntyper.scripts.alignment_index_binding.os.stat", side_effect=_hide_proc_descriptor):
        _, _, plan = _build_view(tmp_path)
    bound_index = Path(plan.stable_index_path)
    bound_index.unlink()
    bound_index.write_bytes(b"replacement entry")
    binding = plan.binding
    assert binding is not None

    with pytest.raises(RuntimeError, match="generated-index binding because it was replaced"):
        plan.close()

    assert bound_index.read_bytes() == b"replacement entry"
    assert binding.is_open is True
    bound_index.unlink()
    plan.close()
    assert binding.is_open is False


def test_clean_hardlink_fallback_rerun_rebuilds_without_a_stale_binding(tmp_path: Path) -> None:
    """A normal close removes the deterministic fallback so the next run can bind it."""
    stat_target = "vntyper.scripts.alignment_index_binding.os.stat"
    with patch(stat_target, side_effect=_hide_proc_descriptor):
        _, _, first = _build_view(tmp_path)
    fallback = Path(first.stable_index_path)
    first.close()
    assert not os.path.lexists(fallback)

    with patch(stat_target, side_effect=_hide_proc_descriptor):
        _, _, second = _build_view(tmp_path)
    try:
        assert Path(second.stable_index_path).read_bytes() == b"fresh index"
    finally:
        second.close()


def test_index_binding_failure_removes_the_alignment_view_and_temporary_index(tmp_path: Path) -> None:
    """Failure before index publication closes the outer alignment binding exactly once."""
    alignment = tmp_path / "sample.bam"
    alignment.write_bytes(b"alignment")
    output = tmp_path / "run"

    def build_index(command: str, *_args: object, **_kwargs: object) -> tuple[bool, str]:
        arguments = shlex.split(command)
        Path(arguments[arguments.index("-o") + 1]).write_bytes(b"fresh index")
        return True, ""

    with (
        patch("vntyper.scripts.alignment_preflight.capture_command", side_effect=build_index),
        patch.object(AlignmentBinding, "bind_index", side_effect=RuntimeError("index bind failed")),
        pytest.raises(RuntimeError, match="index bind failed"),
    ):
        build_alignment_view(
            str(alignment),
            str(output),
            "input",
            "bam",
            {"tools": {"samtools": "samtools"}},
            threads=1,
        )

    assert not os.path.lexists(output / "input.bam")
    assert not os.path.lexists(output / "input.bam.bai")
    assert not tuple(output.glob("*.tmp"))


def test_nonregular_index_binding_closes_its_rejected_descriptor(tmp_path: Path) -> None:
    """A directory cannot be retained as an index, and its descriptor must not leak."""
    index_directory = tmp_path / "index-directory"
    index_directory.mkdir()
    closed: list[int] = []
    real_close = os.close

    def record_close(descriptor: int) -> None:
        closed.append(descriptor)
        real_close(descriptor)

    with (
        patch("vntyper.scripts.alignment_index_binding.os.close", side_effect=record_close),
        pytest.raises(RuntimeError, match="not a regular file"),
    ):
        AlignmentIndexBinding(index_directory, tmp_path / ".bound")

    assert len(closed) == 1


def test_missing_hardlink_fallback_still_releases_the_index_descriptor(tmp_path: Path) -> None:
    """Already-absent owned fallback cleanup is idempotent and closes the FD."""
    index = tmp_path / "fresh.bai"
    fallback = tmp_path / ".fresh.bai.bound"
    index.write_bytes(b"fresh index")
    with patch("vntyper.scripts.alignment_index_binding.os.stat", side_effect=_hide_proc_descriptor):
        binding = AlignmentIndexBinding(index, fallback)
    fallback.unlink()

    binding.close()
    binding.close()

    assert binding.is_open is False


def test_index_fallback_unlink_failure_keeps_the_descriptor_open(tmp_path: Path) -> None:
    """A fallback that cannot be removed must remain backed by its retained FD."""
    index = tmp_path / "fresh.bai"
    fallback = tmp_path / ".fresh.bai.bound"
    index.write_bytes(b"fresh index")
    with patch("vntyper.scripts.alignment_index_binding.os.stat", side_effect=_hide_proc_descriptor):
        binding = AlignmentIndexBinding(index, fallback)
    real_unlink = os.unlink

    def fail_fallback_unlink(path: str | os.PathLike[str]) -> None:
        if Path(path) == fallback:
            raise PermissionError("fallback locked")
        real_unlink(path)

    with (
        patch("vntyper.scripts.alignment_index_binding.os.unlink", side_effect=fail_fallback_unlink),
        pytest.raises(RuntimeError, match="before release"),
    ):
        binding.close()

    assert binding.is_open is True
    binding.close()
    assert binding.is_open is False


@pytest.mark.parametrize("replacement_kind", ["file", "different-symlink", "same-target-symlink"])
def test_plan_close_never_unlinks_a_replaced_view(tmp_path: Path, replacement_kind: str) -> None:
    """Cleanup ownership is the created directory entry, not merely its pathname or target."""
    _, _, plan = _build_view(tmp_path)
    view = Path(plan.view_path)
    owned_target = os.readlink(view)
    view.unlink()
    if replacement_kind == "file":
        view.write_bytes(b"operator replacement")
    elif replacement_kind == "different-symlink":
        replacement_target = tmp_path / "operator-target.bam"
        replacement_target.write_bytes(b"operator replacement")
        view.symlink_to(replacement_target)
    else:
        view.symlink_to(owned_target)

    with pytest.raises(RuntimeError, match="owned alignment view was replaced"):
        plan.close()

    assert os.path.lexists(view)
    if replacement_kind == "file":
        assert view.read_bytes() == b"operator replacement"
    elif replacement_kind == "different-symlink":
        assert view.readlink() == replacement_target
    else:
        assert os.readlink(view) == owned_target
    view.unlink()
    plan.close()


def test_closed_plan_allows_the_same_output_directory_to_pass_the_next_prewrite_guard(tmp_path: Path) -> None:
    """A completed run must not strand a dangling view that rejects a safe rerun."""
    alignment, _, plan = _build_view(tmp_path)
    plan.close()

    protect_alignment_inputs(
        tmp_path / "run",
        alignment,
        "bam",
        None,
        None,
        {"tools": {"samtools": "samtools"}},
        "hg19",
    )


def test_owned_view_unlink_failure_keeps_the_descriptor_open_and_stops_cleanup(tmp_path: Path) -> None:
    """A view that cannot be removed must never become a rebindable closed-FD symlink."""
    _, _, plan = _build_view(tmp_path)
    view = Path(plan.view_path)
    binding = plan.binding
    assert binding is not None
    real_unlink = os.unlink

    def reject_owned_view(path: str | os.PathLike[str], *, dir_fd: int | None = None) -> None:
        if Path(path) == view:
            raise PermissionError("view is locked")
        real_unlink(path, dir_fd=dir_fd)

    with (
        patch("vntyper.scripts.alignment_binding.os.unlink", side_effect=reject_owned_view),
        pytest.raises(RuntimeError, match="before descriptor release"),
    ):
        plan.close()

    assert binding.is_open is True
    assert os.path.lexists(view)
    plan.close()
    assert binding.is_open is False


def test_a_binding_rejects_a_second_view_without_leaking_or_replacing_the_first(tmp_path: Path) -> None:
    """One binding owns exactly one directory entry until explicit close."""
    alignment = tmp_path / "patient.bam"
    alignment.write_bytes(b"patient alignment")
    first_view = tmp_path / "first.bam"
    second_view = tmp_path / "second.bam"
    binding = AlignmentBinding(str(alignment))
    binding.install_view(first_view)

    try:
        with pytest.raises(RuntimeError, match="already owns an alignment view"):
            binding.install_view(second_view)

        assert first_view.read_bytes() == b"patient alignment"
        assert not os.path.lexists(second_view)
    finally:
        binding.close()
        with suppress(OSError):
            first_view.unlink()
        with suppress(OSError):
            second_view.unlink()


def test_gc_preserves_the_descriptor_when_a_same_target_replacement_blocks_cleanup(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Finalization must leak safely instead of letting the proc target rebind to a reused FD."""
    alignment, _, plan = _build_view(tmp_path)
    view = Path(plan.view_path)
    binding = plan.binding
    assert binding is not None
    target = os.readlink(view)
    descriptor = int(target.rsplit("/", maxsplit=1)[1])
    view.unlink()
    view.symlink_to(target)
    reuse_candidate = tmp_path / "reuse-candidate.bin"
    reuse_candidate.write_bytes(b"must never appear through the alignment view")
    reused_descriptor: int | None = None

    try:
        del plan
        del binding
        gc.collect()
        reused_descriptor = os.open(reuse_candidate, os.O_RDONLY)

        assert reused_descriptor != descriptor
        assert view.read_bytes() == alignment.read_bytes()
        assert "Preserving alignment descriptor" in caplog.text
    finally:
        if reused_descriptor is not None:
            with suppress(OSError):
                os.close(reused_descriptor)
        with suppress(OSError):
            os.close(descriptor)
        with suppress(OSError):
            view.unlink()


@pytest.mark.parametrize("input_symlink", [False, True])
def test_proc_unavailable_falls_back_to_a_verified_inode_hardlink(tmp_path: Path, input_symlink: bool) -> None:
    """The fallback must bind the opened inode, including through a supplied symlink."""
    with patch("vntyper.scripts.alignment_binding.os.stat", side_effect=_hide_proc_descriptor):
        alignment, _, plan = _build_view(tmp_path, input_symlink=input_symlink)

    view = Path(plan.view_path)
    assert view.is_symlink() is False
    assert view.read_bytes() == b"original alignment"
    replacement = tmp_path / "replacement.bam"
    replacement.write_bytes(b"replacement alignment")
    replacement.replace(alignment)
    assert view.read_bytes() == b"original alignment"

    plan.close()
    assert not os.path.lexists(view)
    assert alignment.read_bytes() == b"replacement alignment"


def test_hardlink_fallback_refuses_a_path_replaced_between_open_and_link(tmp_path: Path) -> None:
    """Post-link inode verification closes the fallback pathname TOCTOU."""
    alignment = tmp_path / "sample.bam"
    alignment.write_bytes(b"original alignment")
    replacement = tmp_path / "replacement.bam"
    replacement.write_bytes(b"replacement alignment")
    real_link = os.link

    def replace_then_link(
        source: str | os.PathLike[str],
        destination: str | os.PathLike[str],
        *,
        src_dir_fd: int | None = None,
        dst_dir_fd: int | None = None,
        follow_symlinks: bool = True,
    ) -> None:
        replacement.replace(alignment)
        real_link(
            source,
            destination,
            src_dir_fd=src_dir_fd,
            dst_dir_fd=dst_dir_fd,
            follow_symlinks=follow_symlinks,
        )

    with (
        patch("vntyper.scripts.alignment_binding.os.stat", side_effect=_hide_proc_descriptor),
        patch("vntyper.scripts.alignment_binding.os.link", side_effect=replace_then_link),
        patch("vntyper.scripts.alignment_preflight.capture_command") as capture,
        pytest.raises(RuntimeError, match="hardlink does not identify the already-open alignment"),
    ):
        build_alignment_view(
            str(alignment),
            str(tmp_path / "run"),
            "input",
            "bam",
            {"tools": {"samtools": "samtools"}},
            threads=1,
        )

    capture.assert_not_called()
    assert alignment.read_bytes() == b"replacement alignment"
    assert not os.path.lexists(tmp_path / "run" / "input.bam")
    assert not tuple((tmp_path / "run").glob("*.tmp"))


def test_proc_and_cross_filesystem_hardlink_unavailable_refuses_before_index_work(tmp_path: Path) -> None:
    """A platform without either immutable binding must fail precisely before samtools."""
    alignment = tmp_path / "sample.bam"
    alignment.write_bytes(b"original alignment")
    capture_target = "vntyper.scripts.alignment_preflight.capture_command"

    with (
        patch("vntyper.scripts.alignment_binding.os.stat", side_effect=_hide_proc_descriptor),
        patch(
            "vntyper.scripts.alignment_binding.os.link",
            side_effect=OSError(errno.EXDEV, "cross-device link"),
        ),
        patch(capture_target) as capture,
        pytest.raises(
            RuntimeError,
            match="procfs descriptor binding is unavailable and a same-filesystem hardlink could not be created",
        ),
    ):
        build_alignment_view(
            str(alignment),
            str(tmp_path / "run"),
            "input",
            "bam",
            {"tools": {"samtools": "samtools"}},
            threads=1,
        )

    capture.assert_not_called()
    assert alignment.read_bytes() == b"original alignment"
