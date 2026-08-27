"""Private M5 cache bindings for safe terminal CRAM reference resolution."""

from __future__ import annotations

import gzip
import hashlib
import logging
import os
import re
import stat
from collections.abc import Callable, Iterable, Iterator
from contextlib import contextmanager, suppress
from pathlib import Path
from typing import BinaryIO
from urllib.parse import unquote, urlsplit

from vntyper.scripts.reference_binding import _InodeView
from vntyper.scripts.reference_resolution_environment import activate_private_reference_cache
from vntyper.scripts.reference_uri_policy import (
    LocalHeaderReference,
    allow_ambient_reference_resolution,
    remote_ref_path_suffix,
)

logger = logging.getLogger(__name__)

_CACHE_TOKEN = re.compile(r"%(?P<width>\d*)s")


def cache_entry_path(pattern: str, digest: str) -> Path:
    """Expand htslib's progressive ``%Ns`` digest-cache template.

    Args:
        pattern: Local htslib cache template.
        digest: Complete reference sequence digest.

    Returns:
        Filesystem path for that digest.

    Raises:
        ValueError: If the template cannot consume the digest exactly.
    """
    cursor = 0
    pieces: list[str] = []
    previous = 0
    for match in _CACHE_TOKEN.finditer(pattern):
        pieces.append(pattern[previous : match.start()])
        width_text = match.group("width")
        width = int(width_text) if width_text else len(digest) - cursor
        if width <= 0 or cursor + width > len(digest):
            raise ValueError("Invalid local CRAM reference cache template")
        pieces.append(digest[cursor : cursor + width])
        cursor += width
        previous = match.end()
    pieces.append(pattern[previous:])
    if cursor != len(digest):
        raise ValueError("Invalid local CRAM reference cache template")
    return Path("".join(pieces))


def _digest_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _local_cache_patterns(value: str) -> tuple[str, ...]:
    """Split a local REF_PATH list while preserving file-scheme entries."""
    tokens = value.split(os.pathsep)
    result: list[str] = []
    position = 0
    while position < len(tokens):
        token = tokens[position]
        if token.lower() == "file" and position + 1 < len(tokens) and tokens[position + 1].startswith("//"):
            parsed = urlsplit(f"file:{tokens[position + 1]}")
            if parsed.netloc.lower() not in {"", "localhost"} or parsed.query or parsed.fragment:
                result.append("")
            else:
                result.append(unquote(parsed.path))
            position += 2
            continue
        if (
            re.fullmatch(r"[A-Za-z][A-Za-z0-9+.-]*", token)
            and position + 1 < len(tokens)
            and tokens[position + 1].startswith("//")
        ):
            break
        result.append(token)
        position += 1
    return tuple(result)


class PrivateReferenceCache:
    """Retain verified M5 entries in a pipeline-owned cache namespace."""

    def __init__(
        self,
        required: Iterable[tuple[str, str]],
        configured_pattern: str,
        header_references: Iterable[LocalHeaderReference],
        output_dir: str | Path,
        output_name: str,
    ) -> None:
        """Build a complete private cache without modifying source trees.

        Args:
            required: Contig and M5 pairs required by downstream consumers.
            configured_pattern: Pinned local htslib lookup template.
            header_references: Local SQ UR records already parsed at the owned header boundary.
            output_dir: Pipeline-owned preflight directory.
            output_name: Validated alignment output basename.

        Raises:
            ValueError: If a required reference sequence cannot be resolved and verified.
            RuntimeError: If the private namespace cannot be owned or cleaned safely.
        """
        self._root: Path | None = None
        self._directory_identities: dict[Path, tuple[int, int]] = {}
        self._entries: dict[str, _InodeView] = {}
        try:
            root = Path(output_dir) / f".{output_name}_reference_cache"
            try:
                os.mkdir(root, mode=0o700)
            except FileExistsError as error:
                raise RuntimeError(f"reserved CRAM reference cache directory already exists: {root}") from error
            self._root = root
            self._record_directory(root)
            self.pattern = str(root / "%2s" / "%2s" / "%s")
            required_map = dict(required)
            self._bind_configured_entries(required_map, configured_pattern)
            missing = {contig: digest for contig, digest in required_map.items() if digest not in self._entries}
            if missing:
                self._materialize_header_entries(missing, tuple(header_references))
            unresolved = [(contig, digest) for contig, digest in required_map.items() if digest not in self._entries]
            if unresolved:
                contig, digest = unresolved[0]
                raise ValueError(f"Private CRAM reference cache could not resolve contig={contig}, M5={digest}")
        except BaseException:
            try:
                self.close()
            except Exception:
                logger.exception("Private CRAM reference cache cleanup failed while preserving construction failure.")
            raise

    def _record_directory(self, directory: Path) -> None:
        metadata = os.stat(directory, follow_symlinks=False)
        if not stat.S_ISDIR(metadata.st_mode):
            raise RuntimeError(f"Private CRAM reference cache entry is not a directory: {directory}")
        self._directory_identities[directory] = (metadata.st_dev, metadata.st_ino)

    def _ensure_parent(self, destination: Path) -> None:
        if self._root is None:
            message = "private CRAM reference cache used before its root directory was created"
            logger.error(message)
            raise ValueError(message)
        relative = destination.parent.relative_to(self._root)
        current = self._root
        for part in relative.parts:
            current /= part
            if current in self._directory_identities:
                continue
            try:
                os.mkdir(current, mode=0o700)
            except FileExistsError as error:
                raise RuntimeError(f"Private CRAM reference cache path already exists: {current}") from error
            self._record_directory(current)

    def _destination(self, digest: str) -> Path:
        return cache_entry_path(self.pattern, digest)

    def _bind_source(self, source: Path, digest: str) -> bool:
        if not source.is_file():
            return False
        destination = self._destination(digest)
        self._ensure_parent(destination)
        # Own the view before it installs anything: a failed install then still has a
        # reachable owner that can retry its exact-entry removal.
        binding = _InodeView(source, destination)
        self._entries[digest] = binding
        binding.install()
        try:
            matches = _digest_file(destination).lower() == digest.lower()
            if not matches:
                binding.close()
                del self._entries[digest]
            return matches
        except BaseException:
            try:
                binding.close()
            except Exception:
                logger.exception("Verified cache-entry cleanup failed while preserving the digest-read outcome.")
            else:
                del self._entries[digest]
            raise

    def _bind_configured_entries(self, required: dict[str, str], configured_pattern: str) -> None:
        for digest in dict.fromkeys(required.values()):
            for pattern in _local_cache_patterns(configured_pattern):
                try:
                    source = cache_entry_path(pattern, digest)
                except ValueError:
                    continue
                if self._bind_source(source, digest):
                    break

    @contextmanager
    def _open_fasta(self, path: str) -> Iterator[BinaryIO | gzip.GzipFile]:
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            os.close(descriptor)
            raise OSError("local CRAM header reference is not a regular file")
        raw = os.fdopen(descriptor, "rb")
        try:
            if path.lower().endswith((".gz", ".bgz", ".bgzf")):
                with gzip.GzipFile(fileobj=raw) as fasta:
                    yield fasta
            else:
                yield raw
        finally:
            raw.close()

    def _materialize_one_fasta(self, path: str, expected: dict[str, str]) -> None:
        active_contig: str | None = None
        active_digest: str | None = None
        active_path: Path | None = None
        active_handle = None
        hasher = hashlib.md5()

        def finish() -> None:
            nonlocal active_contig, active_digest, active_path, active_handle, hasher
            if active_handle is None or active_digest is None or active_path is None:
                return
            active_handle.close()
            active_handle = None
            if hasher.hexdigest().lower() != active_digest.lower():
                active_path.unlink()
            else:
                self._entries[active_digest] = _InodeView(active_path, active_path, generated=True)
                self._entries[active_digest].install()
            active_contig = None
            active_digest = None
            active_path = None
            hasher = hashlib.md5()

        try:
            with self._open_fasta(path) as fasta:
                for line in fasta:
                    if line.startswith(b">"):
                        finish()
                        contig_bytes = line[1:].split(None, 1)[0]
                        contig = contig_bytes.decode("utf-8", errors="strict")
                        digest = expected.get(contig)
                        if digest is None or digest in self._entries:
                            continue
                        destination = self._destination(digest)
                        self._ensure_parent(destination)
                        active_handle = destination.open("xb")
                        active_contig = contig
                        active_digest = digest
                        active_path = destination
                        continue
                    if active_contig is not None and active_handle is not None:
                        sequence = b"".join(line.split()).upper()
                        active_handle.write(sequence)
                        hasher.update(sequence)
                finish()
        finally:
            if active_handle is not None:
                active_handle.close()
            if active_path is not None and active_digest not in self._entries:
                with suppress(OSError):
                    active_path.unlink()

    def _materialize_header_entries(
        self,
        missing: dict[str, str],
        references: tuple[LocalHeaderReference, ...],
    ) -> None:
        by_path: dict[str, dict[str, str]] = {}
        for reference in references:
            digest = missing.get(reference.contig)
            if digest is None or reference.m5 not in {None, digest}:
                continue
            by_path.setdefault(reference.path, {})[reference.contig] = digest
        for path, expected in by_path.items():
            if all(digest in self._entries for digest in expected.values()):
                continue
            try:
                self._materialize_one_fasta(path, expected)
            except (OSError, UnicodeError):
                logger.warning("Unable to read a local CRAM header reference while building the private cache.")

    def _remove_directories(self) -> None:
        for directory in sorted(self._directory_identities, key=lambda path: len(path.parts), reverse=True):
            identity = self._directory_identities[directory]
            metadata = os.stat(directory, follow_symlinks=False)
            if not stat.S_ISDIR(metadata.st_mode) or (metadata.st_dev, metadata.st_ino) != identity:
                raise RuntimeError(f"Refusing to remove replaced CRAM reference cache directory: {directory}")
            os.rmdir(directory)
            del self._directory_identities[directory]
        self._root = None

    def close(self) -> None:
        """Remove exact private cache entries and their owned namespace."""
        failure: Exception | None = None
        for digest, binding in tuple(self._entries.items())[::-1]:
            try:
                binding.close()
            except Exception as error:
                if failure is None:
                    failure = error
                else:
                    logger.error(f"Additional CRAM reference cache cleanup failure: {error}")
            else:
                del self._entries[digest]
        try:
            self._remove_directories()
        except Exception as error:
            if failure is None:
                failure = error
            else:
                logger.error(f"Additional CRAM reference cache cleanup failure: {error}")
        if failure is not None:
            raise failure


def probe_private_reference_cache(
    *,
    header_contigs: tuple[str, ...],
    header_m5s: Iterable[tuple[str, str]],
    header_references: tuple[LocalHeaderReference, ...],
    config: dict,
    has_remote_header_reference: bool,
    output_dir: str,
    output_name: str,
    position: int,
    probe: Callable[[str | None, int], tuple[bool, str]],
) -> tuple[PrivateReferenceCache | None, str | None]:
    """Build, activate, and prove one complete private terminal cache.

    Returns:
        A retained cache on success, otherwise a path-safe failure reason.
    """
    configured_pattern = os.environ.get("REF_PATH", config.get("cram", {}).get("local_ref_path", "%2s/%2s/%s"))
    ambient_resolution = allow_ambient_reference_resolution(config)
    remote_reference_path = remote_ref_path_suffix(configured_pattern) if ambient_resolution else None
    if ambient_resolution and (remote_reference_path is not None or has_remote_header_reference):
        unique_names = tuple(dict.fromkeys(reference.contig for reference in header_references))
    else:
        unique_names = header_contigs
    checksum_by_contig = dict(header_m5s)
    missing_m5 = next((contig for contig in unique_names if contig not in checksum_by_contig), None)
    if missing_m5 is not None:
        return None, f"header M5 is missing for contig {missing_m5}"
    required = tuple((contig, checksum_by_contig[contig]) for contig in unique_names)
    cache: PrivateReferenceCache | None = None
    try:
        cache = PrivateReferenceCache(
            required,
            configured_pattern,
            header_references,
            output_dir,
            output_name,
        )
        activate_private_reference_cache(cache.pattern, remote_reference_path=remote_reference_path)
        exit_ok, output = probe(None, position)
        if exit_ok:
            return cache, None
        cache.close()
        return None, output.strip() or "probe exited non-zero"
    except ValueError as error:
        return None, str(error)
    except BaseException:
        if cache is not None:
            try:
                cache.close()
            except Exception:
                logger.exception("Private reference cache cleanup failed while preserving the primary outcome.")
        raise
