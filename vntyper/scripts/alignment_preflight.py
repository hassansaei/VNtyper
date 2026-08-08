"""I/O preflight for alignment inputs before any genotyping stage runs."""

from __future__ import annotations

import logging
import os
import re
import subprocess
import tempfile
from collections.abc import Iterable
from contextlib import suppress
from pathlib import Path
from typing import NoReturn

from vntyper.scripts import preflight_error_io as error_io
from vntyper.scripts.alignment_contract import (
    FORMAT_BAM,
    FORMAT_CRAM,
    AlignmentPlan,
    ReferenceAttempt,
    index_candidate_names,
    missing_index_message,
    unresolvable_reference_message,
)
from vntyper.scripts.alignment_index import resolve_any_index, resolve_bam_index
from vntyper.scripts.alignment_index_provenance import (
    _atomic_symlink,
    _install_generated_index,
    _remove_stale_view_indexes,
    _replace_owned_index_with_symlink,
    generated_index_is_owned,
)
from vntyper.scripts.command_builders import (
    build_cram_reference_probe_command,
    build_samtools_idxstats_command,
    build_samtools_index_command,
)
from vntyper.scripts.idxstats_parsing import choose_scan
from vntyper.scripts.reference_resolution import ordered_reference_candidates, uncovered_reference_contigs

logger = logging.getLogger(__name__)

HTSLIB_REFERENCE_SOURCE = "htslib-resolved (header UR: or REF_PATH)"
_REMOTE_REFERENCE_URL = re.compile(r"(?:^|:)(?!file://)[A-Za-z][A-Za-z0-9+.-]*://", re.IGNORECASE)


def _reject(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _same_file(left: str | Path, right: str | Path) -> bool:
    try:
        return os.path.samefile(left, right)
    except OSError:
        return False


def _protected_collision(path: Path, protected_paths: Iterable[str | Path]) -> str | None:
    path_absolute = os.path.abspath(path)
    for protected in protected_paths:
        protected_absolute = os.path.abspath(protected)
        if path_absolute == protected_absolute:
            return protected_absolute
        if os.path.lexists(path) and _same_file(path, protected):
            return protected_absolute
    return None


def _validate_symlink_entry(
    path: Path,
    protected_paths: Iterable[str | Path],
    *,
    allow_regular: bool = False,
) -> None:
    path_absolute = os.path.abspath(path)
    for protected_path in protected_paths:
        protected_absolute = os.path.abspath(protected_path)
        if path_absolute == protected_absolute:
            message = (
                f"Derived output {path} resolves to the same path or same file as protected source {protected_absolute}"
            )
            _reject(message)
    protected = _protected_collision(path, protected_paths)
    if protected is not None and not path.is_symlink():
        message = f"Derived output {path} resolves to the same path or same file as protected source {protected}"
        _reject(message)
    if os.path.lexists(path) and not path.is_symlink():
        if allow_regular and path.is_file():
            return
        message = f"Derived symlink path has a wrong type and will not be replaced: {path}"
        _reject(message)


def _validate_log_entry(path: Path, protected_paths: Iterable[str | Path]) -> None:
    protected = _protected_collision(path, protected_paths)
    if protected is not None:
        message = f"Log path {path} aliases protected source {protected}"
        _reject(message)
    if os.path.lexists(path) and path.is_dir() and not path.is_symlink():
        message = f"Log path has a wrong type and will not be replaced: {path}"
        _reject(message)


def _safe_output_paths(output_dir: str, output_name: str, file_format: str) -> tuple[Path, Path]:
    if not output_name or output_name in {".", ".."} or Path(output_name).name != output_name:
        message = f"Alignment output_name must be a single basename: {output_name}"
        _reject(message)
    output = Path(output_dir)
    view_path = output / f"{output_name}.{file_format}"
    if view_path.parent.resolve(strict=False) != output.resolve(strict=False):
        message = f"Alignment view must stay inside output directory: {view_path}"
        _reject(message)
    if os.path.lexists(output) and not output.is_dir():
        message = f"Alignment output directory has a wrong type: {output}"
        _reject(message)
    return output, view_path


def capture_command(
    command: str,
    log_file: str,
    cwd: str | None = None,
    *,
    protected_paths: Iterable[str | Path] = (),
) -> tuple[bool, str]:
    """Run a shell command while retaining its combined output for parsing.

    Args:
        command: Complete shell command to run under Bash.
        log_file: File receiving the command's combined stdout and stderr.
        cwd: Working directory for the child, or ``None`` to inherit it.
        protected_paths: Paths whose inodes the final log must not alias.

    Returns:
        A success flag and complete output or safe diagnostic. Non-zero exits and
        OS failures are returned so preflight can try another candidate.
    """
    logger.debug(f"Running captured command: {command}")
    final_log = Path(log_file)
    try:
        if not final_log.parent.is_dir():
            raise OSError(f"log directory does not exist: {final_log.parent}")
        _validate_log_entry(final_log, protected_paths)
    except (OSError, ValueError) as error:
        diagnostic = f"Unable to write command log safely: {error}"
        logger.error(diagnostic)
        return False, diagnostic
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=final_log.parent,
            prefix=f".{final_log.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary_log:
            temporary_path = Path(temporary_log.name)
            try:
                completed = subprocess.run(
                    command,
                    shell=True,
                    executable="/bin/bash",
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    cwd=cwd,
                    check=False,
                )
                output = completed.stdout or ""
                return_code = completed.returncode
            except OSError as error:
                output = f"Unable to run command: {error}"
                return_code = 1
            temporary_log.write(output)
        os.replace(temporary_path, final_log)
        temporary_path = None
    except OSError as error:
        diagnostic = f"Unable to write command log safely: {error}"
        logger.error(diagnostic)
        return False, diagnostic
    finally:
        if temporary_path is not None:
            with suppress(OSError):
                os.unlink(temporary_path)
    if return_code != 0:
        logger.error(f"Command failed: {command}")
        return False, output
    return True, output


def build_alignment_view(
    in_path: str,
    output_dir: str,
    output_name: str,
    file_format: str,
    config: dict,
    threads: int,
    *,
    bai_only: bool = False,
) -> tuple[str, str]:
    """Create a run-local alignment symlink and a co-located usable index.

    Args:
        in_path: Input BAM or CRAM path.
        output_dir: Writable run output directory.
        output_name: Base name for the alignment view and its logs.
        file_format: Alignment format, either ``"bam"`` or ``"cram"``.
        config: Pipeline configuration; missing tool configuration uses ``samtools``.
        threads: Thread count for an index build.
        bai_only: Require BAI for a BAM consumer that parses BAI bytes directly.

    Returns:
        The alignment view path and its co-located index path.

    Raises:
        RuntimeError: If samtools cannot build a missing index or does not create it.
        ValueError: If the format, output name, or a derived entry is unsafe.
    """
    input_index_candidates = index_candidate_names(in_path, file_format, bai_only=bai_only)
    output, view_path = _safe_output_paths(output_dir, output_name, file_format)
    existing_index = (
        resolve_bam_index(in_path)
        if bai_only and file_format == FORMAT_BAM
        else resolve_any_index(in_path, file_format)
    )
    built_suffix = {FORMAT_BAM: ".bai", FORMAT_CRAM: ".crai"}[file_format]
    selected_suffix = Path(existing_index).suffix if existing_index is not None else built_suffix
    view_index = Path(f"{view_path}{selected_suffix}")
    view_index_candidates = tuple(Path(candidate) for candidate in index_candidate_names(str(view_path), file_format))
    protected_paths: tuple[str | Path, ...] = (in_path, existing_index) if existing_index is not None else (in_path,)
    _validate_symlink_entry(view_path, (in_path,))
    owned_indexes: set[Path] = set()
    for candidate in view_index_candidates:
        if generated_index_is_owned(candidate, protected_paths):
            owned_indexes.add(candidate)
        _validate_symlink_entry(
            candidate,
            protected_paths,
            allow_regular=candidate in owned_indexes,
        )
    log_file = output / f"{output_name}_index.log"
    if existing_index is None:
        _validate_log_entry(log_file, protected_paths)
    output.mkdir(parents=True, exist_ok=True)
    _atomic_symlink(in_path, view_path)
    if existing_index is not None:
        _remove_stale_view_indexes(str(view_path), file_format, str(view_index), owned_indexes, protected_paths)
        if view_index in owned_indexes:
            _replace_owned_index_with_symlink(view_index, existing_index, protected_paths)
        else:
            _atomic_symlink(existing_index, view_index)
        return str(view_path), str(view_index)
    _remove_stale_view_indexes(str(view_path), file_format, str(view_index), owned_indexes, protected_paths)
    samtools_path = config.get("tools", {}).get("samtools", "samtools")
    temporary_index: Path | None = None
    try:
        descriptor, temporary_name = tempfile.mkstemp(
            dir=output,
            prefix=f".{view_index.name}.",
            suffix=".tmp",
        )
        os.close(descriptor)
        temporary_index = Path(temporary_name)
        command = build_samtools_index_command(
            samtools_path=samtools_path,
            bam_file=view_path,
            output_bai=temporary_index,
            threads=threads,
        )
        exit_ok, _ = capture_command(command, str(log_file), protected_paths=protected_paths)
        if not exit_ok or temporary_index.stat().st_size == 0:
            raise OSError("samtools did not create a non-empty index")
        _install_generated_index(
            temporary_index,
            view_index,
            protected_paths,
            replace_owned=view_index in owned_indexes,
        )
        temporary_index = None
    except OSError as error:
        message = missing_index_message(in_path, file_format, input_index_candidates)
        logger.error(message)
        raise RuntimeError(message) from error
    finally:
        if temporary_index is not None:
            with suppress(OSError):
                os.unlink(temporary_index)
    return str(view_path), str(view_index)


def choose_unmapped_scan(
    view_path: str,
    config: dict,
    threads: int,
    output_dir: str,
    output_name: str,
) -> str:
    """Choose a lossless CRAM unmapped-read scan from captured idxstats evidence.

    Args:
        view_path: Run-local indexed alignment path.
        config: Pipeline configuration. Missing keys use shipped defaults.
        threads: Thread count for idxstats.
        output_dir: Run output directory for the idxstats log.
        output_name: Base name for the idxstats log.

    Returns:
        ``"indexed"`` only when idxstats proves it lossless; otherwise ``"stream"``.

    Raises:
        ValueError: If scan mode is invalid or forced indexed mode would lose reads.
    """
    samtools_path = config.get("tools", {}).get("samtools", "samtools")
    configured = config.get("cram", {}).get("unmapped_scan", "auto")
    command = build_samtools_idxstats_command(
        samtools_path=samtools_path,
        in_bam=view_path,
        threads=threads,
    )
    log_file = str(Path(output_dir) / f"{output_name}_idxstats.log")
    exit_ok, output = capture_command(command, log_file)
    scan, reason = choose_scan(configured, output, exit_ok)
    logger.info(f"Selected {scan} unmapped-read scan: {reason}")
    return scan


def pin_reference_resolution(config: dict) -> str | None:
    """Pin htslib reference lookup to local paths unless ambient lookup is allowed.

    Args:
        config: Pipeline configuration. Missing CRAM keys use shipped defaults.

    Returns:
        The prior ``REF_PATH``, including ``None``; restore it in a ``finally``.

    Raises:
        ValueError: If the supposedly local ``REF_PATH`` contains a remote URL.
    """
    previous = os.environ.get("REF_PATH")
    cram_config = config.get("cram", {})
    if cram_config.get("allow_ambient_reference_resolution", False):
        logger.warning("Ambient CRAM reference resolution is enabled and may block on a network endpoint.")
        return previous
    local_ref_path = cram_config.get("local_ref_path", "%2s/%2s/%s")
    if _REMOTE_REFERENCE_URL.search(local_ref_path):
        message = (
            "cram.local_ref_path must not contain a remote URL when ambient reference resolution is disabled: "
            f"{local_ref_path}"
        )
        _reject(message)
    os.environ["REF_PATH"] = local_ref_path
    return previous


def restore_reference_resolution(previous: str | None) -> None:
    """Restore the ``REF_PATH`` value captured before a run.

    Args:
        previous: Original ``REF_PATH`` value, or ``None`` when it was unset.
    """
    if previous is None:
        os.environ.pop("REF_PATH", None)
    else:
        os.environ["REF_PATH"] = previous


def _reference_contigs(reference_path: str) -> set[str] | None:
    fai_path = Path(f"{reference_path}.fai")
    try:
        if fai_path.is_file():
            with fai_path.open() as index_handle:
                return {line.split("\t", 1)[0] for line in index_handle if line.strip()}
    except OSError:
        return None
    return None


def _reference_unavailable_reason(reference_path: str) -> str | None:
    try:
        with open(reference_path, "rb") as reference_handle:
            reference_handle.read(1)
    except FileNotFoundError:
        return "reference FASTA not found"
    except OSError as error:
        return f"reference FASTA unreadable: {error}"
    return None


def _target_contig(
    region: str | None,
    bed_file: str | Path | None,
    header_contigs: tuple[str, ...],
) -> str:
    if bed_file is not None:
        try:
            with open(bed_file) as bed_handle:
                for line in bed_handle:
                    stripped = line.strip()
                    if stripped and not stripped.startswith("#"):
                        return stripped.split(maxsplit=1)[0]
        except OSError:
            pass
    elif region:
        return region.split(":", 1)[0]
    return header_contigs[0] if header_contigs else "unknown"


def resolve_reference(
    view_path: str,
    candidates: Iterable[tuple[str, str | None]],
    region: str | None,
    bed_file: str | Path | None,
    config: dict,
    threads: int,
    output_dir: str,
    output_name: str,
    header_contigs: Iterable[str],
    m5: str | None,
    *,
    failure_context: error_io.PreflightErrorContext | None = None,
) -> tuple[str | None, str, tuple[str, ...]]:
    """Probe explicit CRAM references before one final htslib-resolved candidate.

    Args:
        view_path: Run-local indexed CRAM path.
        candidates: Explicit ``(source, path)`` candidates in policy order.
        region: This run's region string when no BED file is used.
        bed_file: This run's BED file. It takes precedence over ``region``.
        config: Pipeline configuration; missing tool configuration uses ``samtools``.
        threads: Thread count for each probe.
        output_dir: Run output directory for probe logs.
        output_name: Base name for probe logs.
        header_contigs: Contigs declared by the alignment header.
        m5: Header M5 checksum for the target contig, when available.
        failure_context: Boundary diagnostics receiving structured candidate details.

    Returns:
        The winning path, source, and uncovered header contigs. An htslib-resolved
        winner has path ``None`` because the no-``-T`` probe cannot identify it.

    Raises:
        ValueError: If no candidate decodes the selected target.
    """
    samtools_path = config.get("tools", {}).get("samtools", "samtools")
    attempts: list[ReferenceAttempt] = []
    explicit_candidates = list(candidates)
    for position, (source, reference_path) in enumerate(explicit_candidates, start=1):
        if reference_path is None:
            attempts.append((source, None, "not supplied"))
            continue
        unavailable_reason = _reference_unavailable_reason(reference_path)
        if unavailable_reason is not None:
            attempts.append((source, reference_path, unavailable_reason))
            continue
        command = build_cram_reference_probe_command(
            samtools_path=samtools_path,
            in_bam=view_path,
            region=region,
            bed_file=bed_file,
            reference_path=reference_path,
            threads=threads,
        )
        log_file = str(Path(output_dir) / f"{output_name}_reference_probe_{position}.log")
        exit_ok, output = capture_command(command, log_file)
        if not exit_ok:
            attempts.append((source, reference_path, output.strip() or "probe exited non-zero"))
            continue

        uncovered_decision = uncovered_reference_contigs(header_contigs, _reference_contigs(reference_path))
        uncovered: tuple[str, ...]
        if uncovered_decision is None:
            uncovered = ()
            logger.warning(
                f"Resolved CRAM reference from {source}, but reference-contig coverage unavailable because "
                f"{reference_path}.fai is missing or unreadable"
            )
        else:
            uncovered = uncovered_decision
        if uncovered:
            logger.warning(
                f"Resolved CRAM reference from {source}, but it does not cover header contigs: {', '.join(uncovered)}"
            )
        logger.info(f"Resolved CRAM reference from {source}: {reference_path}")
        return reference_path, source, uncovered
    command = build_cram_reference_probe_command(
        samtools_path=samtools_path,
        in_bam=view_path,
        region=region,
        bed_file=bed_file,
        reference_path=None,
        threads=threads,
    )
    ambient_position = len(explicit_candidates) + 1
    log_file = str(Path(output_dir) / f"{output_name}_reference_probe_{ambient_position}.log")
    exit_ok, output = capture_command(command, log_file)
    if exit_ok:
        logger.info(f"Resolved CRAM reference through {HTSLIB_REFERENCE_SOURCE}")
        return None, HTSLIB_REFERENCE_SOURCE, ()
    attempts.append((HTSLIB_REFERENCE_SOURCE, None, output.strip() or "probe exited non-zero"))
    contigs = tuple(header_contigs)
    target_contig = _target_contig(region, bed_file, contigs)
    message = unresolvable_reference_message(view_path, target_contig, m5, attempts)
    if failure_context is not None:
        failure_context.payload = error_io.public_reference_error_payload(target_contig, m5, attempts)
    _reject(message)


def _validate_preflight_logs(
    in_path: str,
    output_dir: str,
    output_name: str,
    file_format: str,
    candidates: tuple[tuple[str, str | None], ...],
    *,
    bai_only: bool,
) -> None:
    output, _ = _safe_output_paths(output_dir, output_name, file_format)
    existing_index = (
        resolve_bam_index(in_path)
        if bai_only and file_format == FORMAT_BAM
        else resolve_any_index(in_path, file_format)
    )
    protected_paths: tuple[str | Path, ...] = (in_path, existing_index) if existing_index is not None else (in_path,)
    if file_format == FORMAT_CRAM:
        log_paths = [output / f"{output_name}_idxstats.log"]
        log_paths.extend(
            output / f"{output_name}_reference_probe_{position}.log" for position in range(1, len(candidates) + 2)
        )
    else:
        log_paths = [output / f"{output_name}_alignment_probe.log"]
    for log_path in log_paths:
        _validate_log_entry(log_path, protected_paths)


def run_preflight(
    in_path: str,
    output_dir: str,
    output_name: str,
    file_format: str,
    config: dict,
    threads: int,
    *,
    region: str | None = None,
    bed_file: str | Path | None = None,
    reference_assembly: str = "hg19",
    reference_fasta: str | None = None,
    header_contigs: Iterable[str] = (),
    m5: str | None = None,
    fast_mode: bool = False,
    error_output_dir: str | Path | None = None,
) -> AlignmentPlan:
    """Resolve and prove every alignment prerequisite used by later stages.

    Args:
        in_path: Input BAM or CRAM path.
        output_dir: Writable run output directory.
        output_name: Base name for the alignment view and preflight logs.
        file_format: Alignment format, either ``"bam"`` or ``"cram"``.
        config: Pipeline configuration; absent new keys use shipped defaults.
        threads: Thread count for samtools commands.
        region: This run's target region when no BED file is used.
        bed_file: This run's BED target. It takes precedence over ``region``.
        reference_assembly: Assembly suffix for configured reference paths.
        reference_fasta: Explicit CLI or web reference candidate.
        header_contigs: Contigs declared by the alignment header.
        m5: Header M5 checksum for the target contig, when available.
        fast_mode: Whether downstream BAM processing can use any htslib index.
        error_output_dir: Run output root for the curated failure artifact.

    Returns:
        A frozen plan whose index, scan, and CRAM reference have been exercised.

    Raises:
        RuntimeError: If an indexed BAM cannot retrieve the requested target.
        ValueError: If format, scan, candidate policy, target, or reference is invalid.
    """
    failure_context = error_io.PreflightErrorContext(error_output_dir if error_output_dir is not None else output_dir)
    with error_io.persist_preflight_failure(failure_context):
        bai_only = file_format == FORMAT_BAM and not fast_mode
        failure_context.phase = error_io.REFERENCE_POLICY_FAILURE
        candidates = (
            ordered_reference_candidates(config, reference_assembly, reference_fasta)
            if file_format == FORMAT_CRAM
            else ()
        )
        failure_context.phase = error_io.OUTPUT_SAFETY_FAILURE
        _validate_preflight_logs(in_path, output_dir, output_name, file_format, candidates, bai_only=bai_only)
        failure_context.phase = error_io.VIEW_INDEX_FAILURE
        view_path, index_path = build_alignment_view(
            in_path,
            output_dir,
            output_name,
            file_format,
            config,
            threads,
            bai_only=bai_only,
        )
        if file_format == FORMAT_CRAM:
            failure_context.phase = error_io.SCAN_SELECTION_FAILURE
            unmapped_scan = choose_unmapped_scan(view_path, config, threads, output_dir, output_name)
            failure_context.phase = error_io.REFERENCE_PROBE_FAILURE
            reference_path, reference_source, uncovered_contigs = resolve_reference(
                view_path,
                candidates,
                region,
                bed_file,
                config,
                threads,
                output_dir,
                output_name,
                tuple(header_contigs),
                m5,
                failure_context=failure_context,
            )
        else:
            failure_context.phase = error_io.BAM_PROBE_FAILURE
            command = build_cram_reference_probe_command(
                samtools_path=config.get("tools", {}).get("samtools", "samtools"),
                in_bam=view_path,
                region=region,
                bed_file=bed_file,
                reference_path=None,
                threads=threads,
            )
            log_file = str(Path(output_dir) / f"{output_name}_alignment_probe.log")
            exit_ok, output = capture_command(command, log_file)
            if not exit_ok:
                reason = output.strip() or "probe exited non-zero"
                message = f"BAM preflight probe failed: {reason}"
                logger.error(message)
                raise RuntimeError(message)
            unmapped_scan, reference_path, reference_source, uncovered_contigs = "indexed", None, "not-required", ()

        return AlignmentPlan(
            input_path=in_path,
            view_path=view_path,
            file_format=file_format,
            index_path=index_path,
            reference_path=reference_path,
            reference_source=reference_source,
            uncovered_contigs=uncovered_contigs,
            unmapped_scan=unmapped_scan,
        )
