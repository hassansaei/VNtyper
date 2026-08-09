"""I/O preflight for alignment inputs before any genotyping stage runs."""

from __future__ import annotations

import logging
import os
import tempfile
from collections.abc import Iterable
from contextlib import nullcontext, suppress
from pathlib import Path
from typing import NoReturn

from vntyper.scripts import preflight_error_io as error_io
from vntyper.scripts.alignment_binding import AlignmentBinding
from vntyper.scripts.alignment_contract import (
    FORMAT_BAM,
    FORMAT_CRAM,
    AlignmentPlan,
    ReferenceAttempt,
    index_candidate_names,
    missing_index_message,
    missing_reference_contig,
    unresolvable_reference_message,
)
from vntyper.scripts.alignment_index import resolve_any_index, resolve_bam_index
from vntyper.scripts.alignment_index_provenance import (
    _install_generated_index,
    _remove_stale_view_indexes,
    generated_index_is_owned,
)
from vntyper.scripts.alignment_preflight_logs import preflight_log_paths
from vntyper.scripts.alignment_scan import select_unmapped_scan
from vntyper.scripts.command_builders import (
    build_cram_reference_probe_command,
    build_cram_stream_reference_probe_command,
    build_samtools_depth_command,
    build_samtools_index_command,
)
from vntyper.scripts.preflight_command_io import (
    capture_command,
)
from vntyper.scripts.preflight_command_io import (
    protected_collision as _protected_collision,
)
from vntyper.scripts.preflight_command_io import (
    validate_log_entry as _validate_log_entry,
)
from vntyper.scripts.preflight_input_io import (
    configured_preflight_text_limit,
    read_bounded_regular_text,
)
from vntyper.scripts.reference_binding import ReferenceBinding
from vntyper.scripts.reference_binding import reference_contigs as _reference_contigs
from vntyper.scripts.reference_binding import reference_probe_timeout_seconds as _reference_probe_timeout_seconds
from vntyper.scripts.reference_binding import reference_unavailable_reason as _reference_unavailable_reason
from vntyper.scripts.reference_binding import target_contig as _target_contig
from vntyper.scripts.reference_resolution import ordered_reference_candidates, uncovered_reference_contigs
from vntyper.scripts.reference_resolution_environment import (
    pin_reference_resolution as pin_reference_resolution,
)
from vntyper.scripts.reference_resolution_environment import (
    restore_reference_resolution as restore_reference_resolution,
)

logger = logging.getLogger(__name__)

HTSLIB_REFERENCE_SOURCE = "htslib-resolved (header UR: or REF_PATH)"
HEADER_UR_REFERENCE_SOURCE = "header_ur"


def _reject(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


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


def build_alignment_view(
    in_path: str,
    output_dir: str,
    output_name: str,
    file_format: str,
    config: dict,
    threads: int,
    *,
    bai_only: bool = False,
    binding: AlignmentBinding | None = None,
    bound_view_path: str | Path | None = None,
) -> tuple[str, str, AlignmentBinding]:
    """Create a bound run-local alignment and freshly built co-located index.

    Args:
        in_path: Input BAM or CRAM path.
        output_dir: Writable run output directory.
        output_name: Base name for the alignment view and its logs.
        file_format: Alignment format, either ``"bam"`` or ``"cram"``.
        config: Pipeline configuration; missing tool configuration uses ``samtools``.
        threads: Thread count for an index build.
        bai_only: Use the legacy BAI-only BAM preflight/protection contract.
        binding: Binding installed before validation, when present.
        bound_view_path: Early-installed view owned by ``binding``.

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
    view_index = Path(f"{view_path}{built_suffix}")
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
    _validate_log_entry(log_file, protected_paths)
    owns_new_binding = binding is None
    if binding is None:
        binding = AlignmentBinding(in_path)
    elif str(bound_view_path) != str(view_path) or binding.view_path != str(view_path):
        message = "Pre-bound alignment view does not match the preflight destination."
        logger.error(message)
        raise ValueError(message)
    temporary_index: Path | None = None
    try:
        output.mkdir(parents=True, exist_ok=True)
        if owns_new_binding:
            binding.install_view(view_path)
        _remove_stale_view_indexes(str(view_path), file_format, str(view_index), owned_indexes, protected_paths)
        samtools_path = config.get("tools", {}).get("samtools", "samtools")
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
        binding.bind_index(temporary_index, output / f".{view_index.name}.bound")
        _install_generated_index(
            temporary_index,
            view_index,
            protected_paths,
            replace_owned=view_index in owned_indexes,
        )
        temporary_index = None
    except OSError as error:
        binding.close()
        message = missing_index_message(in_path, file_format, input_index_candidates)
        logger.error(message)
        raise RuntimeError(message) from error
    except BaseException:
        binding.close()
        raise
    finally:
        if temporary_index is not None:
            with suppress(OSError):
                os.unlink(temporary_index)
    return str(view_path), str(view_index), binding


def choose_unmapped_scan(
    view_path: str,
    config: dict,
    threads: int,
    output_dir: str,
    output_name: str,
    *,
    file_format: str = FORMAT_CRAM,
) -> str:
    """Choose a lossless alignment unmapped-read scan from captured idxstats evidence.

    Args:
        view_path: Run-local indexed alignment path.
        config: Pipeline configuration. Missing keys use shipped defaults.
        threads: Thread count for idxstats.
        output_dir: Run output directory for the idxstats log.
        output_name: Base name for the idxstats log.
        file_format: Alignment format whose scan policy is being selected.

    Returns:
        ``"indexed"`` only when idxstats proves it lossless; otherwise ``"stream"``.

    Raises:
        ValueError: If scan mode is invalid or forced indexed mode would lose reads.
    """
    return select_unmapped_scan(
        view_path,
        config,
        threads,
        output_dir,
        output_name,
        file_format=file_format,
        capture=capture_command,
    )


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
    coverage_region: str | None = None,
    header_m5s: Iterable[tuple[str, str]] = (),
    unmapped_scan: str = "indexed",
    failure_context: error_io.PreflightErrorContext | None = None,
) -> tuple[str | None, str, tuple[str, ...], ReferenceBinding | None]:
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
        coverage_region: Independently resolved region used by the later depth consumer.
        header_m5s: Header contig and M5 pairs used to diagnose a later
            whole-file stream failure.
        unmapped_scan: Downstream CRAM unmapped-read scan. Stream mode requires
            an additional whole-file decode proof from the same candidate.
        failure_context: Boundary diagnostics receiving structured candidate details.

    Returns:
        The winning run-local path, source, uncovered header contigs, and retained
        reference binding. An htslib-resolved winner has path and binding ``None``
        because the no-``-T`` probe cannot identify its backing file.

    Raises:
        ValueError: If no candidate decodes the selected target.
    """
    samtools_path = config.get("tools", {}).get("samtools", "samtools")
    timeout_seconds = _reference_probe_timeout_seconds(config)
    max_text_bytes = configured_preflight_text_limit(config)
    attempts: list[ReferenceAttempt] = []
    explicit_candidates = list(candidates)

    def probe_candidate(
        reference_path: str | None,
        position: int,
        reference_binding: ReferenceBinding | None = None,
    ) -> tuple[bool, str]:
        target_command = build_cram_reference_probe_command(
            samtools_path=samtools_path,
            in_bam=view_path,
            region=region,
            bed_file=bed_file,
            reference_path=reference_path,
            threads=threads,
        )
        target_log = str(Path(output_dir) / f"{output_name}_reference_probe_{position}.log")
        exit_ok, probe_output = capture_command(
            target_command,
            target_log,
            timeout_seconds=timeout_seconds,
        )
        if reference_binding is not None:
            reference_binding.bind_generated_sidecars()
        if not exit_ok:
            return exit_ok, probe_output
        if coverage_region is not None:
            coverage_command = build_samtools_depth_command(
                samtools_path=samtools_path,
                threads=threads,
                region=coverage_region,
                bam_file=view_path,
                coverage_output="/dev/null",
                reference_path=reference_path,
            )
            coverage_log = str(Path(output_dir) / f"{output_name}_reference_coverage_probe_{position}.log")
            exit_ok, probe_output = capture_command(
                coverage_command,
                coverage_log,
                timeout_seconds=timeout_seconds,
            )
        if not exit_ok or unmapped_scan != "stream":
            return exit_ok, probe_output
        stream_command = build_cram_stream_reference_probe_command(
            samtools_path=samtools_path,
            in_bam=view_path,
            reference_path=reference_path,
            threads=threads,
        )
        stream_log = str(Path(output_dir) / f"{output_name}_reference_stream_probe_{position}.log")
        return capture_command(
            stream_command,
            stream_log,
            timeout_seconds=timeout_seconds,
        )

    for position, (source, reference_path) in enumerate(explicit_candidates, start=1):
        if reference_path is None:
            attempts.append((source, None, "not supplied"))
            continue
        unavailable_reason = _reference_unavailable_reason(reference_path)
        if unavailable_reason is not None:
            attempts.append((source, reference_path, unavailable_reason))
            continue
        reference_binding = ReferenceBinding(reference_path, output_dir, output_name, position)
        try:
            exit_ok, output = probe_candidate(reference_binding.consumer_path, position, reference_binding)
        except BaseException:
            try:
                reference_binding.close()
            except Exception:
                logger.exception("Reference binding cleanup failed while preserving the primary probe outcome.")
            raise
        if not exit_ok:
            attempts.append((source, reference_path, output.strip() or "probe exited non-zero"))
            reference_binding.close()
            continue

        uncovered_decision = uncovered_reference_contigs(
            header_contigs,
            _reference_contigs(reference_binding.consumer_path, max_text_bytes),
        )
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
        return reference_binding.consumer_path, source, uncovered, reference_binding
    if any(source == HEADER_UR_REFERENCE_SOURCE for source, _reference_path in explicit_candidates):
        contigs = tuple(header_contigs)
        target_contig = _target_contig(region, bed_file, contigs, max_text_bytes)
        failure_contig = next(
            (
                parsed
                for _, _, reason in reversed(attempts)
                if (parsed := missing_reference_contig(reason, contigs)) is not None
            ),
            target_contig,
        )
        failure_m5 = dict(header_m5s).get(failure_contig, m5 if failure_contig == target_contig else None)
        message = unresolvable_reference_message(view_path, failure_contig, failure_m5, attempts)
        if failure_context is not None:
            failure_context.payload = error_io.public_reference_error_payload(failure_contig, failure_m5, attempts)
        _reject(message)

    ambient_position = len(explicit_candidates) + 1
    exit_ok, output = probe_candidate(None, ambient_position)
    if exit_ok:
        logger.info(f"Resolved CRAM reference through {HTSLIB_REFERENCE_SOURCE}")
        return None, HTSLIB_REFERENCE_SOURCE, (), None
    attempts.append((HTSLIB_REFERENCE_SOURCE, None, output.strip() or "probe exited non-zero"))
    contigs = tuple(header_contigs)
    target_contig = _target_contig(region, bed_file, contigs, max_text_bytes)
    failure_contig = next(
        (
            parsed
            for _, _, reason in reversed(attempts)
            if (parsed := missing_reference_contig(reason, contigs)) is not None
        ),
        target_contig,
    )
    failure_m5 = dict(header_m5s).get(failure_contig, m5 if failure_contig == target_contig else None)
    message = unresolvable_reference_message(view_path, failure_contig, failure_m5, attempts)
    if failure_context is not None:
        failure_context.payload = error_io.public_reference_error_payload(failure_contig, failure_m5, attempts)
    _reject(message)


def _validate_preflight_logs(
    in_path: str,
    output_dir: str,
    output_name: str,
    file_format: str,
    candidates: tuple[tuple[str, str | None], ...],
    *,
    bai_only: bool,
    fast_mode: bool,
    coverage_region: str | None,
) -> None:
    output, _ = _safe_output_paths(output_dir, output_name, file_format)
    existing_index = (
        resolve_bam_index(in_path)
        if bai_only and file_format == FORMAT_BAM
        else resolve_any_index(in_path, file_format)
    )
    protected_paths: tuple[str | Path, ...] = (in_path, existing_index) if existing_index is not None else (in_path,)
    for log_path in preflight_log_paths(
        output,
        output_name,
        file_format,
        candidate_count=len(candidates),
        fast_mode=fast_mode,
        coverage_region=coverage_region,
    ):
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
    coverage_region: str | None = None,
    reference_assembly: str = "hg19",
    reference_fasta: str | None = None,
    header_contigs: Iterable[str] = (),
    header_reference_paths: Iterable[str] = (),
    m5: str | None = None,
    header_m5s: Iterable[tuple[str, str]] = (),
    fast_mode: bool = False,
    error_output_dir: str | Path | None = None,
    failure_context: error_io.PreflightErrorContext | None = None,
    binding: AlignmentBinding | None = None,
    bound_view_path: str | Path | None = None,
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
        coverage_region: Independently resolved region used by the later depth consumer.
        reference_assembly: Assembly suffix for configured reference paths.
        reference_fasta: Explicit CLI or web reference candidate.
        header_contigs: Contigs declared by the alignment header.
        header_reference_paths: Local FASTA paths parsed from CRAM SQ UR tags.
        m5: Header M5 checksum for the target contig, when available.
        header_m5s: Header contig and M5 pairs for terminal stream diagnostics.
        fast_mode: Whether downstream BAM processing can use any htslib index.
        error_output_dir: Run output root for the curated failure artifact.
        failure_context: Optional outer-owned failure context.
        binding: Binding installed before pipeline validation, when present.
        bound_view_path: Run-local view already owned by ``binding``.

    Raises:
        RuntimeError: If an indexed BAM cannot retrieve the requested target.
        ValueError: If format, scan, candidate policy, target, or reference is invalid.
    """
    owns_failure_context = failure_context is None
    if failure_context is None:
        failure_context = error_io.PreflightErrorContext(
            error_output_dir if error_output_dir is not None else output_dir
        )
    persistence = error_io.persist_preflight_failure(failure_context) if owns_failure_context else nullcontext()
    with persistence:
        max_text_bytes = configured_preflight_text_limit(config)
        if bed_file is not None:
            read_bounded_regular_text(
                bed_file,
                max_bytes=max_text_bytes,
                description="alignment target BED",
            )
        bai_only = file_format == FORMAT_BAM and not fast_mode
        failure_context.phase = error_io.REFERENCE_POLICY_FAILURE
        candidates: tuple[tuple[str, str | None], ...] = ()
        if file_format == FORMAT_CRAM:
            candidates = (
                *ordered_reference_candidates(config, reference_assembly, reference_fasta),
                *((HEADER_UR_REFERENCE_SOURCE, path) for path in header_reference_paths),
            )
        if file_format == FORMAT_CRAM:
            _reference_probe_timeout_seconds(config)
        failure_context.phase = error_io.OUTPUT_SAFETY_FAILURE
        _validate_preflight_logs(
            in_path,
            output_dir,
            output_name,
            file_format,
            candidates,
            bai_only=bai_only,
            fast_mode=fast_mode,
            coverage_region=coverage_region,
        )
        failure_context.phase = error_io.VIEW_INDEX_FAILURE
        view_path, index_path, binding = build_alignment_view(
            in_path,
            output_dir,
            output_name,
            file_format,
            config,
            threads,
            bai_only=bai_only,
            binding=binding,
            bound_view_path=bound_view_path,
        )
        reference_binding: ReferenceBinding | None = None
        try:
            unmapped_scan = "not-required"
            if not fast_mode:
                failure_context.phase = error_io.SCAN_SELECTION_FAILURE
                unmapped_scan = choose_unmapped_scan(
                    view_path,
                    config,
                    threads,
                    output_dir,
                    output_name,
                    file_format=file_format,
                )
            if file_format == FORMAT_CRAM:
                failure_context.phase = error_io.REFERENCE_PROBE_FAILURE
                reference_path, reference_source, uncovered_contigs, reference_binding = resolve_reference(
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
                    coverage_region=coverage_region,
                    header_m5s=header_m5s,
                    unmapped_scan=unmapped_scan,
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
                reference_path, reference_source, uncovered_contigs = None, "not-required", ()

            return AlignmentPlan(
                input_path=in_path,
                view_path=view_path,
                file_format=file_format,
                index_path=index_path,
                reference_path=reference_path,
                reference_source=reference_source,
                uncovered_contigs=uncovered_contigs,
                unmapped_scan=unmapped_scan,
                binding=binding,
                reference_binding=reference_binding,
            )
        except BaseException:
            if reference_binding is not None:
                try:
                    reference_binding.close()
                except Exception:
                    logger.exception("Reference binding cleanup failed while preserving the primary preflight outcome.")
            try:
                binding.close()
            except Exception:
                logger.exception("Alignment binding cleanup failed while preserving the primary preflight outcome.")
            raise
