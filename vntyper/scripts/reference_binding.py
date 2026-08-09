"""Lifetime-owned run-local views of explicit CRAM reference files."""

from __future__ import annotations

import logging
import math
import os
import secrets
import stat
import tempfile
from contextlib import suppress
from pathlib import Path

from vntyper.scripts.preflight_input_io import (
    read_bounded_regular_text,
    regular_file_unavailable_reason,
    try_read_bounded_regular_text,
)

logger = logging.getLogger(__name__)

REFERENCE_SIDECAR_SUFFIXES = (".fai", ".gzi")
MAX_REFERENCE_PROBE_TIMEOUT_SECONDS = 120.0


def reference_probe_timeout_seconds(config: dict) -> float:
    """Return a validated reference-probe deadline no greater than 120 seconds.

    Args:
        config: Pipeline configuration containing optional CRAM probe policy.

    Returns:
        The validated timeout in seconds.

    Raises:
        ValueError: If the configured timeout is not finite and within bounds.
    """
    configured = config.get("cram", {}).get("reference_probe_timeout_seconds", MAX_REFERENCE_PROBE_TIMEOUT_SECONDS)
    if (
        isinstance(configured, bool)
        or not isinstance(configured, (int, float))
        or not math.isfinite(configured)
        or configured <= 0
        or configured > MAX_REFERENCE_PROBE_TIMEOUT_SECONDS
    ):
        message = "cram.reference_probe_timeout_seconds must be a finite number greater than 0 and at most 120"
        logger.error(message)
        raise ValueError(message)
    return float(configured)


def reference_contigs(reference_path: str, max_bytes: int) -> set[str] | None:
    """Read bounded contig names from a run-local reference FAI.

    Args:
        reference_path: Stable run-local FASTA path.
        max_bytes: Maximum accepted FAI byte count.

    Returns:
        The indexed contig names, or ``None`` when the optional FAI is unavailable.
    """
    index_text, reason = try_read_bounded_regular_text(
        Path(f"{reference_path}.fai"),
        max_bytes=max_bytes,
        description="reference FASTA index",
    )
    if reason is not None or index_text is None:
        return None
    return {line.split("\t", 1)[0] for line in index_text.splitlines() if line.strip()}


def reference_unavailable_reason(reference_path: str) -> str | None:
    """Return why an operator reference is unusable.

    Args:
        reference_path: Operator-owned FASTA candidate.

    Returns:
        An actionable reason, or ``None`` when the file is readable and regular.
    """
    return regular_file_unavailable_reason(reference_path, description="reference FASTA")


def target_contig(
    region: str | None,
    bed_file: str | Path | None,
    header_contigs: tuple[str, ...],
    max_bytes: int,
) -> str:
    """Choose the target contig used in terminal reference diagnostics.

    Args:
        region: Region target when no BED is supplied.
        bed_file: Optional BED target, which takes precedence.
        header_contigs: Contigs declared by the CRAM header.
        max_bytes: Maximum accepted BED byte count.

    Returns:
        The first active BED contig, region contig, first header contig, or ``unknown``.
    """
    if bed_file is not None:
        bed_text = read_bounded_regular_text(
            bed_file,
            max_bytes=max_bytes,
            description="alignment target BED",
        )
        for line in bed_text.splitlines():
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                return stripped.split(maxsplit=1)[0]
    elif region:
        return region.split(":", 1)[0]
    return header_contigs[0] if header_contigs else "unknown"


class _InodeView:
    """Retain one regular-file inode behind an owned directory entry."""

    def __init__(self, source: str | Path, destination: str | Path, *, generated: bool = False) -> None:
        self._descriptor: int | None = None
        self._identity: tuple[int, int] | None = None
        self._destination = Path(destination)
        self._destination_identity: tuple[int, int] | None = None
        self._destination_kind: str | None = None
        self._proc_target: str | None = None
        flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NONBLOCK", 0)
        if generated:
            flags |= getattr(os, "O_NOFOLLOW", 0)
        try:
            descriptor = os.open(source, flags)
            self._descriptor = descriptor
            metadata = os.fstat(descriptor)
            if not stat.S_ISREG(metadata.st_mode):
                raise OSError("bound reference entry is not a regular file")
            self._identity = (metadata.st_dev, metadata.st_ino)
            self._proc_target = f"/proc/{os.getpid()}/fd/{descriptor}"
            if generated:
                self._replace_generated_entry()
            else:
                self._install_source_view(source)
        except OSError as error:
            self._close_descriptor()
            message = f"Unable to retain CRAM reference input {source}: {error}"
            logger.error(message)
            raise RuntimeError(message) from error

    @property
    def is_open(self) -> bool:
        """Return whether this exact input inode remains retained."""
        return self._descriptor is not None

    def _temporary_path(self) -> Path:
        for _ in range(100):
            temporary = self._destination.with_name(f".{self._destination.name}.{secrets.token_hex(8)}.tmp")
            if not os.path.lexists(temporary):
                return temporary
        raise OSError(f"unable to allocate a temporary reference view beside {self._destination}")

    def _proc_target_is_exact(self) -> bool:
        assert self._proc_target is not None
        try:
            metadata = os.stat(self._proc_target)
        except OSError:
            return False
        return (metadata.st_dev, metadata.st_ino) == self._identity

    def _record_destination(self, kind: str) -> None:
        metadata = os.lstat(self._destination)
        self._destination_identity = (metadata.st_dev, metadata.st_ino)
        self._destination_kind = kind

    def _install_proc_link(self, *, replace: bool) -> bool:
        if not self._proc_target_is_exact():
            return False
        assert self._proc_target is not None
        if not replace:
            try:
                os.symlink(self._proc_target, self._destination)
            except OSError:
                return False
            self._record_destination("symlink")
            return True

        temporary = self._temporary_path()
        try:
            os.symlink(self._proc_target, temporary)
            os.replace(temporary, self._destination)
        except OSError:
            with suppress(OSError):
                os.unlink(temporary)
            return False
        self._record_destination("symlink")
        return True

    def _install_source_view(self, source: str | Path) -> None:
        if os.path.lexists(self._destination):
            raise OSError(f"run-local reference view already exists: {self._destination}")
        if self._install_proc_link(replace=False):
            return
        if os.path.lexists(self._destination):
            raise OSError(f"run-local reference view already exists: {self._destination}")
        temporary = self._temporary_path()
        try:
            os.link(os.path.realpath(source), temporary, follow_symlinks=True)
            metadata = os.stat(temporary, follow_symlinks=False)
            if not stat.S_ISREG(metadata.st_mode) or (metadata.st_dev, metadata.st_ino) != self._identity:
                raise OSError("hardlink does not identify the already-open reference input")
            os.link(temporary, self._destination, follow_symlinks=False)
            os.unlink(temporary)
        except OSError:
            with suppress(OSError):
                os.unlink(temporary)
            raise
        self._record_destination("regular")

    def _replace_generated_entry(self) -> None:
        metadata = os.lstat(self._destination)
        if not stat.S_ISREG(metadata.st_mode) or (metadata.st_dev, metadata.st_ino) != self._identity:
            raise OSError("generated reference sidecar changed before it could be retained")
        if self._install_proc_link(replace=True):
            return
        self._record_destination("regular")

    def _remove_destination(self) -> None:
        if self._destination_identity is None or self._destination_kind is None:
            return
        try:
            metadata = os.lstat(self._destination)
        except FileNotFoundError:
            self._destination_identity = None
            self._destination_kind = None
            return
        except OSError as error:
            raise RuntimeError(f"Unable to inspect owned CRAM reference view {self._destination}: {error}") from error
        correct_kind = (
            stat.S_ISLNK(metadata.st_mode) if self._destination_kind == "symlink" else stat.S_ISREG(metadata.st_mode)
        )
        correct_target = (
            self._destination_kind != "symlink"
            or correct_kind
            and self._proc_target is not None
            and os.readlink(self._destination) == self._proc_target
        )
        if (metadata.st_dev, metadata.st_ino) != self._destination_identity or not correct_kind or not correct_target:
            raise RuntimeError(f"Refusing to remove replaced CRAM reference view: {self._destination}")
        try:
            os.unlink(self._destination)
        except OSError as error:
            raise RuntimeError(f"Unable to remove owned CRAM reference view {self._destination}: {error}") from error
        self._destination_identity = None
        self._destination_kind = None

    def _close_descriptor(self) -> None:
        descriptor = self._descriptor
        if descriptor is None:
            return
        self._descriptor = None
        with suppress(OSError):
            os.close(descriptor)

    def close(self) -> None:
        """Remove the exact owned entry before releasing its descriptor."""
        self._remove_destination()
        self._close_descriptor()


class ReferenceBinding:
    """Hold one explicit reference and its htslib sidecars for a plan lifetime."""

    def __init__(self, input_path: str, output_dir: str | Path, output_name: str, position: int) -> None:
        """Install a stable run-local reference namespace without copying FASTA bytes.

        Args:
            input_path: Operator-owned explicit or configured FASTA path.
            output_dir: Pipeline-owned preflight directory.
            output_name: Validated alignment output basename.
            position: One-based reference-candidate position.

        Raises:
            RuntimeError: If the reference cannot be bound to pipeline-owned storage.
        """
        self.input_path = input_path
        self._directory: Path | None = None
        self._directory_identity: tuple[int, int] | None = None
        self._reference: _InodeView | None = None
        self._sidecars: dict[str, _InodeView] = {}
        suffixes = "".join(Path(input_path).suffixes) or ".fa"
        try:
            directory = Path(tempfile.mkdtemp(prefix=f".{output_name}_reference_{position}.", dir=output_dir))
            self._directory = directory
            directory_stat = os.stat(directory, follow_symlinks=False)
            self._directory_identity = (directory_stat.st_dev, directory_stat.st_ino)
            self.consumer_path = str(directory / f"reference{suffixes}")
            self._reference = _InodeView(input_path, self.consumer_path)
            for suffix in REFERENCE_SIDECAR_SUFFIXES:
                source = Path(f"{input_path}{suffix}")
                if regular_file_unavailable_reason(source, description=f"reference FASTA {suffix} index") is None:
                    destination = Path(f"{self.consumer_path}{suffix}")
                    self._sidecars[suffix] = _InodeView(source, destination)
        except BaseException:
            self._close_after_failed_construction()
            raise

    @property
    def is_open(self) -> bool:
        """Return whether the bound reference descriptor remains available."""
        return self._reference is not None and self._reference.is_open

    def bind_generated_sidecars(self) -> None:
        """Retain `.fai` and `.gzi` files htslib generated in the private namespace."""
        for suffix in REFERENCE_SIDECAR_SUFFIXES:
            if suffix in self._sidecars:
                continue
            destination = Path(f"{self.consumer_path}{suffix}")
            if os.path.lexists(destination):
                self._sidecars[suffix] = _InodeView(destination, destination, generated=True)

    def _remove_directory(self) -> None:
        if self._directory is None or self._directory_identity is None:
            return
        try:
            metadata = os.stat(self._directory, follow_symlinks=False)
        except FileNotFoundError:
            self._directory = None
            self._directory_identity = None
            return
        except OSError as error:
            raise RuntimeError(
                f"Unable to inspect owned CRAM reference directory {self._directory}: {error}"
            ) from error
        if not stat.S_ISDIR(metadata.st_mode) or (metadata.st_dev, metadata.st_ino) != self._directory_identity:
            raise RuntimeError(f"Refusing to remove replaced CRAM reference directory: {self._directory}")
        try:
            os.rmdir(self._directory)
        except OSError as error:
            raise RuntimeError(f"Unable to remove owned CRAM reference directory {self._directory}: {error}") from error
        self._directory = None
        self._directory_identity = None

    def _close_after_failed_construction(self) -> None:
        for binding in reversed(tuple(self._sidecars.values())):
            with suppress(Exception):
                binding.close()
        self._sidecars.clear()
        if self._reference is not None:
            with suppress(Exception):
                self._reference.close()
            self._reference = None
        with suppress(Exception):
            self._remove_directory()

    def close(self) -> None:
        """Remove retained sidecars and reference view, then their private directory."""
        failure: Exception | None = None

        def record_failure(error: Exception) -> None:
            nonlocal failure
            if failure is None:
                failure = error
            else:
                logger.error(f"Additional CRAM reference cleanup failure: {error}")

        for suffix in reversed(REFERENCE_SIDECAR_SUFFIXES):
            binding = self._sidecars.get(suffix)
            if binding is not None:
                try:
                    binding.close()
                except Exception as error:
                    record_failure(error)
                else:
                    del self._sidecars[suffix]
        if self._reference is not None:
            try:
                self._reference.close()
            except Exception as error:
                record_failure(error)
            else:
                self._reference = None
        try:
            self._remove_directory()
        except Exception as error:
            record_failure(error)
        if failure is not None:
            raise failure

    def __del__(self) -> None:
        """Preserve descriptors when replaced entries make safe cleanup impossible."""
        try:
            self.close()
        except Exception as error:
            with suppress(Exception):
                logger.error(f"Preserving CRAM reference descriptors because cleanup was refused: {error}")
