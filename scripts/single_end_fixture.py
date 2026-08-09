"""Derive a single-end BAM from a declared paired-end cohort fixture."""

from __future__ import annotations

import logging
import os
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import pysam

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SingleEndFixture:
    """A validated, portable single-end fixture declaration."""

    name: str
    source_bam: Path
    output_bam: Path
    root: Path


def _index_path(output_bam: Path) -> Path:
    return output_bam.with_suffix(".bam.bai")


def _lexists(path: Path) -> bool:
    return os.path.lexists(path)


def _require_contained(root: Path, path: Path, label: str) -> None:
    if not path.resolve(strict=False).is_relative_to(root.resolve()):
        raise ValueError(f"Invalid single-end fixture declaration: {label} must be contained by the repository root.")


def _validate_output_pair(output_bam: Path) -> bool:
    output_index = _index_path(output_bam)
    for artifact in (output_bam, output_index):
        if _lexists(artifact) and (artifact.is_symlink() or not artifact.is_file()):
            raise ValueError(f"Single-end fixture output artifact must be a regular file, not a link: {artifact}")
    bam_exists = _lexists(output_bam)
    index_exists = _lexists(output_index)
    if bam_exists != index_exists:
        raise ValueError(f"Single-end fixture output must be absent or a complete BAM/BAI pair: {output_bam}")
    return bam_exists


def _reject_source_alias(source_bam: Path, output_bam: Path) -> None:
    if _lexists(source_bam) and _lexists(output_bam):
        try:
            aliases_source = os.path.samefile(source_bam, output_bam)
        except OSError as exc:
            raise ValueError(f"Could not validate single-end fixture source/output identity: {exc}") from exc
        if aliases_source:
            raise ValueError("Single-end fixture source_bam and output_bam resolve to the same file.")


def _unique_sibling(artifact: Path, suffix: str) -> Path:
    descriptor, name = tempfile.mkstemp(prefix=f".{artifact.name}.", suffix=suffix, dir=artifact.parent)
    os.close(descriptor)
    path = Path(name)
    path.unlink()
    return path


def _remove_if_present(path: Path) -> None:
    if _lexists(path):
        path.unlink()


def _install_pair(
    temporary_bam: Path,
    temporary_index: Path,
    output_bam: Path,
    output_index: Path,
    *,
    replacing_pair: bool,
) -> None:
    backup_bam = _unique_sibling(output_bam, ".backup.bam") if replacing_pair else None
    backup_index = _unique_sibling(output_index, ".backup.bai") if replacing_pair else None
    backups_ready = False
    try:
        if backup_bam is not None and backup_index is not None:
            os.link(output_bam, backup_bam)
            os.link(output_index, backup_index)
            backups_ready = True
        os.replace(temporary_bam, output_bam)
        os.replace(temporary_index, output_index)
    except BaseException:
        if backups_ready and backup_bam is not None and backup_index is not None:
            os.replace(backup_bam, output_bam)
            os.replace(backup_index, output_index)
        elif not replacing_pair:
            _remove_if_present(output_bam)
            _remove_if_present(output_index)
        raise
    finally:
        for backup in (backup_bam, backup_index):
            if backup is not None:
                _remove_if_present(backup)


def parse_single_end_fixture(
    entry: Mapping[str, object],
    *,
    root: Path = Path("."),
) -> SingleEndFixture:
    """Validate one ``single_end_bam`` fixture declaration.

    Args:
        entry: Mapping from ``tests/test_data_config.json``.
        root: Repository root used to resolve the declaration's relative paths.

    Returns:
        The validated fixture name and resolved paths.

    Raises:
        ValueError: If the declaration is incomplete, has the wrong kind, uses
            absolute paths, or would overwrite its own source.
    """
    name = entry.get("name")
    kind = entry.get("kind")
    source_value = entry.get("source_bam")
    output_value = entry.get("output_bam")
    values = (name, source_value, output_value)
    if kind != "single_end_bam" or any(not isinstance(value, str) or not value.strip() for value in values):
        raise ValueError(
            "Invalid single-end fixture declaration: kind must be single_end_bam and name/source_bam/output_bam "
            "must be non-empty strings."
        )

    assert isinstance(name, str)
    assert isinstance(source_value, str)
    assert isinstance(output_value, str)
    source_relative = Path(source_value)
    output_relative = Path(output_value)
    paths = (source_relative, output_relative)
    if Path(name).name != name:
        raise ValueError("Invalid single-end fixture declaration: name must be a portable basename.")
    if any(path.is_absolute() or ".." in path.parts or path.suffix != ".bam" for path in paths):
        raise ValueError(
            "Invalid single-end fixture declaration: source_bam and output_bam must be relative BAM paths "
            "contained by the repository."
        )

    root_resolved = root.resolve()
    source_bam = root / source_relative
    output_bam = root / output_relative
    _require_contained(root_resolved, source_bam, "source_bam")
    _require_contained(root_resolved, output_bam, "output_bam")
    if source_bam == output_bam:
        raise ValueError("Invalid single-end fixture declaration: output_bam must differ from source_bam.")
    _reject_source_alias(source_bam, output_bam)
    _validate_output_pair(output_bam)
    return SingleEndFixture(name=name, source_bam=source_bam, output_bam=output_bam, root=root_resolved)


def derive_single_end_bam(spec: SingleEndFixture) -> int:
    """Copy every alignment while clearing all paired-end relationship flags.

    The read name, sequence, qualities, alignment and every unrelated flag stay
    unchanged. The output is indexed after the complete BAM has been written.

    Args:
        spec: Validated fixture declaration.

    Returns:
        The number of records copied to the single-end BAM.

    Raises:
        FileNotFoundError: If the declared source BAM is absent.
        ValueError: If source and destination resolve to the same file.
    """
    _require_contained(spec.root, spec.source_bam, "source_bam")
    _require_contained(spec.root, spec.output_bam, "output_bam")
    if not spec.source_bam.is_file():
        raise FileNotFoundError(f"Single-end fixture source BAM does not exist: {spec.source_bam}")
    _reject_source_alias(spec.source_bam, spec.output_bam)
    replacing_pair = _validate_output_pair(spec.output_bam)

    spec.output_bam.parent.mkdir(parents=True, exist_ok=True)
    _require_contained(spec.root, spec.output_bam, "output_bam")
    temporary_bam = _unique_sibling(spec.output_bam, ".bam")
    temporary_index = _unique_sibling(_index_path(spec.output_bam), ".bai")
    records = 0
    try:
        with (
            pysam.AlignmentFile(str(spec.source_bam), "rb") as source,
            pysam.AlignmentFile(str(temporary_bam), "wb", template=source) as output,
        ):
            for read in source.fetch(until_eof=True):
                read.is_paired = False
                read.is_read1 = False
                read.is_read2 = False
                read.is_proper_pair = False
                read.mate_is_unmapped = False
                read.mate_is_reverse = False
                output.write(read)
                records += 1
        pysam.samtools.index("-o", str(temporary_index), str(temporary_bam))
        _reject_source_alias(spec.source_bam, spec.output_bam)
        replacing_pair = _validate_output_pair(spec.output_bam)
        _install_pair(
            temporary_bam,
            temporary_index,
            spec.output_bam,
            _index_path(spec.output_bam),
            replacing_pair=replacing_pair,
        )
    finally:
        _remove_if_present(temporary_bam)
        _remove_if_present(temporary_index)
    logger.info("Derived %s with %d single-end records at %s", spec.name, records, spec.output_bam)
    return records
