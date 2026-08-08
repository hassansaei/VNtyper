"""Safe installation of generated alignment target BED files."""

from __future__ import annotations

import logging
import os
import stat
import tempfile
from contextlib import suppress
from pathlib import Path

from vntyper.scripts.alignment_contract import index_candidate_names

logger = logging.getLogger(__name__)


def _reject(message: str) -> None:
    logger.error(message)
    raise ValueError(message)


def _absolute(path: str | Path) -> Path:
    return Path(os.path.abspath(path))


def _same_file(left: Path, right: Path) -> bool:
    try:
        return os.path.samefile(left, right)
    except OSError:
        return False


def _protected_alignment_paths(input_path: str | Path, file_format: str) -> tuple[Path, ...]:
    paths = [_absolute(input_path)]
    paths.extend(
        _absolute(candidate)
        for candidate in index_candidate_names(str(input_path), file_format)
        if os.path.lexists(candidate)
    )
    return tuple(paths)


def validate_alignment_output_root(
    output_root: str | Path,
    input_path: str | Path,
    file_format: str,
) -> None:
    """Reject an output root that could modify an alignment's input tree.

    Args:
        output_root: Intended pipeline or generated-target output directory.
        input_path: Protected BAM/CRAM input.
        file_format: Alignment format used to enumerate protected indexes.

    Raises:
        ValueError: If the root is unsafe or aliases protected input state.
    """
    root = Path(output_root)
    root_absolute = _absolute(root)
    root_variants = (root_absolute, root_absolute.resolve(strict=False))
    logical_input_tree = _absolute(input_path).parent
    resolved_input_tree = _absolute(input_path).resolve(strict=False).parent
    for root_variant in root_variants:
        if root_variant in (logical_input_tree, resolved_input_tree):
            _reject(f"Alignment output root must stay outside the patient input tree: {root}")

    if os.path.lexists(root) and not root.is_dir():
        _reject(f"Alignment output root must be a directory: {root}")
    for protected in _protected_alignment_paths(input_path, file_format):
        if root_absolute == protected or _same_file(root, protected):
            _reject(f"Alignment output root aliases protected alignment input or index: {root}")


def _validate_generated_destination(
    destination: Path,
    input_path: str | Path | None,
    file_format: str | None,
) -> None:
    destination_absolute = _absolute(destination)
    if os.path.lexists(destination):
        if destination.is_symlink():
            _reject(f"Generated BED destination must not be a symlink: {destination}")
        if not stat.S_ISREG(os.lstat(destination).st_mode):
            _reject(f"Generated BED destination must be a regular file: {destination}")

    if input_path is None or file_format is None:
        return

    validate_alignment_output_root(destination.parent, input_path, file_format)
    protected_paths = _protected_alignment_paths(input_path, file_format)
    for protected in protected_paths:
        if destination_absolute == protected or _same_file(destination, protected):
            _reject(f"Generated BED destination aliases protected alignment input or index: {destination}")


def install_generated_bed(
    destination: str | Path,
    bed_text: str,
    *,
    input_path: str | Path | None,
    file_format: str | None,
) -> None:
    """Atomically install generated BED text without following a destination link.

    Args:
        destination: Final generated BED path.
        bed_text: Fully formatted BED contents.
        input_path: Protected BAM/CRAM path, or ``None`` for FASTQ input.
        file_format: Protected alignment format, or ``None`` for FASTQ input.

    Raises:
        OSError: If the temporary write or atomic replacement fails.
        ValueError: If the output could mutate the input tree or an unsafe entry.
    """
    final_path = Path(destination)
    _validate_generated_destination(final_path, input_path, file_format)
    previous_mode = stat.S_IMODE(os.lstat(final_path).st_mode) if os.path.lexists(final_path) else None
    descriptor, temporary_name = tempfile.mkstemp(
        dir=final_path.parent,
        prefix=f".{final_path.name}.",
        suffix=".tmp",
    )
    temporary_path = Path(temporary_name)
    try:
        if previous_mode is not None:
            os.fchmod(descriptor, previous_mode)
        with os.fdopen(descriptor, "w", encoding="utf-8") as temporary_file:
            descriptor = -1
            temporary_file.write(bed_text)
            temporary_file.flush()
            os.fsync(temporary_file.fileno())
        os.replace(temporary_path, final_path)
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        with suppress(OSError):
            os.unlink(temporary_path)
