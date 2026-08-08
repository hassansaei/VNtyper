"""I/O preflight for alignment inputs before any genotyping stage runs."""

from __future__ import annotations

import logging
import os
import re
import subprocess
from collections.abc import Iterable
from pathlib import Path

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
from vntyper.scripts.command_builders import (
    build_cram_reference_probe_command,
    build_samtools_idxstats_command,
    build_samtools_index_command,
)
from vntyper.scripts.idxstats_parsing import choose_scan
from vntyper.scripts.utils import run_command

logger = logging.getLogger(__name__)

HTSLIB_REFERENCE_SOURCE = "htslib-resolved (header UR: or REF_PATH)"
_AMBIENT_SOURCE_ALIASES = {"none", "header_ur", "htslib", "htslib_resolved", HTSLIB_REFERENCE_SOURCE}
_DEFAULT_REFERENCE_ORDER = ("cli", "config_cram_reference", "config_bwa_reference", "htslib_resolved")
_EXPLICIT_REFERENCE_SOURCES = {"cli", "config_cram_reference", "config_bwa_reference"}
_REMOTE_REFERENCE_URL = re.compile(r"(?:^|:)(?!file://)[A-Za-z][A-Za-z0-9+.-]*://", re.IGNORECASE)


def capture_command(command: str, log_file: str, cwd: str | None = None) -> tuple[bool, str]:
    """Run a shell command while retaining its combined output for parsing.

    Args:
        command: Complete shell command to run under Bash.
        log_file: File receiving the command's combined stdout and stderr.
        cwd: Working directory for the child, or ``None`` to inherit it.

    Returns:
        A success flag and the complete captured output. A non-zero exit is
        returned rather than raised so preflight can try another candidate.
    """
    logger.debug(f"Running captured command: {command}")
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
    with open(log_file, "w") as log_handle:
        log_handle.write(output)
    if completed.returncode != 0:
        logger.error(f"Command failed: {command}")
        return False, output
    return True, output


def _ensure_symlink(target: str, link_path: str) -> None:
    """Create a symlink, retaining only an existing link to the same file."""
    absolute_target = os.path.abspath(target)
    try:
        os.symlink(absolute_target, link_path)
        return
    except FileExistsError:
        pass

    same_target = False
    try:
        existing_target = os.readlink(link_path)
        if not os.path.isabs(existing_target):
            existing_target = os.path.join(os.path.dirname(link_path), existing_target)
        same_target = os.path.samefile(existing_target, absolute_target)
    except OSError:
        same_target = False

    if same_target:
        return
    os.unlink(link_path)
    os.symlink(absolute_target, link_path)


def _remove_stale_view_indexes(view_path: str, file_format: str, keep: str | None) -> None:
    """Remove co-located indexes other than the one selected for this run."""
    for candidate in index_candidate_names(view_path, file_format):
        if candidate != keep and os.path.lexists(candidate):
            os.unlink(candidate)


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
        config: Pipeline configuration. Missing tool configuration uses the shipped
            ``samtools`` default.
        threads: Thread count for an index build.
        bai_only: Require BAI for a BAM consumer that parses BAI bytes directly.

    Returns:
        The alignment view path and its co-located index path.

    Raises:
        RuntimeError: If samtools cannot build a missing index or does not create it.
        ValueError: If ``file_format`` is unsupported, the view escapes the output
            directory, or the input and view have the same path.
    """
    input_index_candidates = index_candidate_names(in_path, file_format, bai_only=bai_only)
    output = Path(output_dir)
    view_path = output / f"{output_name}.{file_format}"
    output_absolute = os.path.abspath(output)
    view_absolute = os.path.abspath(view_path)
    input_absolute = os.path.abspath(in_path)
    try:
        view_is_contained = os.path.commonpath((output_absolute, view_absolute)) == output_absolute
    except ValueError:
        view_is_contained = False
    if not view_is_contained:
        message = f"Alignment view must stay inside output directory: {view_absolute}"
        logger.error(message)
        raise ValueError(message)
    if input_absolute == view_absolute:
        message = f"Input alignment and alignment view resolve to the same path: {input_absolute}"
        logger.error(message)
        raise ValueError(message)

    output.mkdir(parents=True, exist_ok=True)
    _ensure_symlink(in_path, str(view_path))

    existing_index = (
        resolve_bam_index(in_path)
        if bai_only and file_format == FORMAT_BAM
        else resolve_any_index(in_path, file_format)
    )
    if existing_index is not None:
        index_suffix = Path(existing_index).suffix
        view_index = Path(f"{view_path}{index_suffix}")
        _remove_stale_view_indexes(str(view_path), file_format, str(view_index))
        _ensure_symlink(existing_index, str(view_index))
        return str(view_path), str(view_index)

    built_suffix = {FORMAT_BAM: ".bai", FORMAT_CRAM: ".crai"}[file_format]
    view_index = Path(f"{view_path}{built_suffix}")
    _remove_stale_view_indexes(str(view_path), file_format, None)
    samtools_path = config.get("tools", {}).get("samtools", "samtools")
    command = build_samtools_index_command(
        samtools_path=samtools_path,
        bam_file=view_path,
        threads=threads,
    )
    log_file = output / f"{output_name}_index.log"
    if not run_command(command, str(log_file), critical=False) or not view_index.exists():
        message = missing_index_message(in_path, file_format, input_index_candidates)
        logger.error(message)
        raise RuntimeError(message)
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
        ``"indexed"`` only when idxstats proves it lossless, otherwise
        ``"stream"``.

    Raises:
        ValueError: If the configured scan mode is invalid or explicitly forcing
            indexed mode would discard placed-unmapped reads.
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
        The previous ``REF_PATH`` value, including ``None`` when it was unset.
        Pass this value to :func:`restore_reference_resolution` in a ``finally``.

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
        logger.error(message)
        raise ValueError(message)
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
    """Read contig names from a FASTA index without scanning a whole genome."""
    fai_path = Path(f"{reference_path}.fai")
    try:
        if fai_path.is_file():
            with fai_path.open() as index_handle:
                return {line.split("\t", 1)[0] for line in index_handle if line.strip()}
    except OSError:
        return None
    return None


def _target_contig(
    region: str | None,
    bed_file: str | Path | None,
    header_contigs: tuple[str, ...],
) -> str:
    """Identify the target contig used by the real probe for diagnostics."""
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
) -> tuple[str | None, str, tuple[str, ...]]:
    """Probe explicit CRAM references before one final htslib-resolved candidate.

    Args:
        view_path: Run-local indexed CRAM path.
        candidates: Explicit ``(source, path)`` candidates in policy order. A
            missing path is recorded as not supplied and is not probed.
        region: This run's region string when no BED file is used.
        bed_file: This run's BED file. It takes precedence over ``region``.
        config: Pipeline configuration. Missing tool configuration uses the
            shipped ``samtools`` default.
        threads: Thread count for each probe.
        output_dir: Run output directory for probe logs.
        output_name: Base name for probe logs.
        header_contigs: Contigs declared by the alignment header.
        m5: Header M5 checksum for the target contig, when available.

    Returns:
        The winning reference path, its source, and header contigs not covered by
        an explicit winning reference. An htslib-resolved winner has path
        ``None`` because a no-``-T`` probe cannot identify which fallback won.

    Raises:
        ValueError: If no candidate decodes the target, or neither a region nor
            BED file was supplied.
    """
    samtools_path = config.get("tools", {}).get("samtools", "samtools")
    attempts: list[ReferenceAttempt] = []
    explicit_candidates = [(source, path) for source, path in candidates if source not in _AMBIENT_SOURCE_ALIASES]

    for position, (source, reference_path) in enumerate(explicit_candidates, start=1):
        if reference_path is None:
            attempts.append((source, None, "not supplied"))
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

        reference_contigs = _reference_contigs(reference_path)
        uncovered = (
            tuple(contig for contig in header_contigs if contig not in reference_contigs)
            if reference_contigs is not None
            else ()
        )
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
    logger.error(message)
    raise ValueError(message)


def _ordered_reference_candidates(
    config: dict,
    reference_assembly: str,
    reference_fasta: str | None,
) -> tuple[tuple[str, str | None], ...]:
    """Build explicit candidates in configured order, leaving ambient lookup implicit."""
    order = config.get("cram", {}).get("reference_candidate_order", list(_DEFAULT_REFERENCE_ORDER))
    reference_data = config.get("reference_data", {})
    values = {
        "cli": reference_fasta,
        "config_cram_reference": reference_data.get(f"cram_reference_{reference_assembly}"),
        "config_bwa_reference": reference_data.get(f"bwa_reference_{reference_assembly}"),
    }
    candidates: list[tuple[str, str | None]] = []
    for source in order:
        if source in _AMBIENT_SOURCE_ALIASES:
            continue
        if source not in _EXPLICIT_REFERENCE_SOURCES:
            message = f"unknown reference candidate source: {source}"
            logger.error(message)
            raise ValueError(message)
        candidates.append((source, values[source]))
    return tuple(candidates)


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
) -> AlignmentPlan:
    """Resolve and prove every alignment prerequisite used by later stages.

    Args:
        in_path: Input BAM or CRAM path.
        output_dir: Writable run output directory.
        output_name: Base name for the alignment view and preflight logs.
        file_format: Alignment format, either ``"bam"`` or ``"cram"``.
        config: Pipeline configuration. Every new key uses its shipped default
            when absent.
        threads: Thread count for samtools commands.
        region: This run's target region when no BED file is used.
        bed_file: This run's BED target. It takes precedence over ``region``.
        reference_assembly: Assembly suffix for configured reference paths.
        reference_fasta: Explicit CLI or web reference candidate.
        header_contigs: Contigs declared by the alignment header.
        m5: Header M5 checksum for the target contig, when available.
        fast_mode: Whether downstream BAM processing can use any htslib index.

    Returns:
        A frozen alignment plan whose index, scan, and CRAM reference have been
        exercised before stage work starts.

    Raises:
        RuntimeError: If an indexed BAM cannot retrieve the requested target.
        ValueError: If the format, scan mode, candidate policy, target, or CRAM
            reference is invalid.
    """
    bai_only = file_format == FORMAT_BAM and not fast_mode
    view_path, index_path = build_alignment_view(
        in_path,
        output_dir,
        output_name,
        file_format,
        config,
        threads,
        bai_only=bai_only,
    )
    contigs = tuple(header_contigs)

    if file_format == FORMAT_CRAM:
        unmapped_scan = choose_unmapped_scan(view_path, config, threads, output_dir, output_name)
        candidates = _ordered_reference_candidates(config, reference_assembly, reference_fasta)
        reference_path, reference_source, uncovered_contigs = resolve_reference(
            view_path,
            candidates,
            region,
            bed_file,
            config,
            threads,
            output_dir,
            output_name,
            contigs,
            m5,
        )
    else:
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
        unmapped_scan = "indexed"
        reference_path = None
        reference_source = "not-required"
        uncovered_contigs = ()

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
