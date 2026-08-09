"""Consume a proven alignment plan for coverage and release its binding."""

from __future__ import annotations

from collections.abc import Callable

from vntyper.scripts.alignment_contract import AlignmentPlan


def calculate_alignment_coverage(
    *,
    plan: AlignmentPlan,
    region: str | None,
    reference_assembly: str,
    threads: int,
    config: dict,
    output_dir: str,
    coverage_calculator: Callable[..., object],
    region_resolver: Callable[..., str],
) -> str:
    """Run the final alignment consumer and release its descriptor binding.

    Args:
        plan: Proven alignment paths, reference, and descriptor binding.
        region: Independently preflighted coverage region, when already resolved.
        reference_assembly: Declared assembly for legacy fallback resolution.
        threads: Samtools thread count.
        config: Pipeline configuration.
        output_dir: Coverage artifact directory.
        coverage_calculator: Coverage stage callable.
        region_resolver: Legacy region fallback callable.

    Returns:
        The exact region consumed by the coverage stage.
    """
    if region is None:
        region = region_resolver(
            bam_file=plan.view_path,
            reference_assembly=reference_assembly,
            region_type="vntr_region",
            config=config,
        )
    try:
        coverage_calculator(
            bam_file=plan.view_path,
            region=region,
            threads=threads,
            config=config,
            output_dir=output_dir,
            output_name="coverage",
            reference_path=plan.reference_path,
        )
    finally:
        plan.close()
    return region
