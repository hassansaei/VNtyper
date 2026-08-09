"""Pure destination planning for alignment-preflight command logs."""

from pathlib import Path

from vntyper.scripts.alignment_contract import FORMAT_CRAM


def preflight_log_paths(
    output_dir: str | Path,
    output_name: str,
    file_format: str,
    *,
    candidate_count: int,
    fast_mode: bool,
    coverage_region: str | None,
) -> tuple[Path, ...]:
    """Return every command log destination reachable by a preflight plan.

    Args:
        output_dir: Alignment-preflight output directory.
        output_name: Base name shared by generated logs.
        file_format: Alignment format, either BAM or CRAM.
        candidate_count: Number of explicit CRAM reference candidates.
        fast_mode: Whether the unused unmapped-read scan is skipped.
        coverage_region: Independently resolved depth-consumer region, when any.

    Returns:
        Ordered log paths that the selected preflight mode can write.
    """
    output = Path(output_dir)
    if file_format != FORMAT_CRAM:
        paths = [output / f"{output_name}_alignment_probe.log"]
        if not fast_mode:
            paths.insert(0, output / f"{output_name}_idxstats.log")
        return tuple(paths)

    paths = [] if fast_mode else [output / f"{output_name}_idxstats.log"]
    positions = range(1, candidate_count + 2)
    paths.extend(output / f"{output_name}_reference_probe_{position}.log" for position in positions)
    if coverage_region is not None:
        positions = range(1, candidate_count + 2)
        paths.extend(output / f"{output_name}_reference_coverage_probe_{position}.log" for position in positions)
    if not fast_mode:
        positions = range(1, candidate_count + 2)
        paths.extend(output / f"{output_name}_reference_stream_probe_{position}.log" for position in positions)
    return tuple(paths)
