"""Consume a proven alignment plan for coverage and release its binding."""

from __future__ import annotations

from collections.abc import Callable

from vntyper.scripts.alignment_contract import AlignmentPlan
from vntyper.scripts.pipeline_cleanup import close_alignment_plan


def _assembly_geometry(config: dict, reference_assembly: str) -> dict | None:
    """The `bam_processing.assemblies` entry for this run, or None.

    The alias has to be resolved first. `assemblies` is keyed by coordinate system -
    `GRCh37`, `GRCh38` - while the caller carries the user-facing name, and the default
    is `hg19`. Looking the alias up directly missed on the default and on six of the
    eight supported names, so #222's columns wrote "not measured" on almost every real
    run; the tests did not catch it because they pass the assembly dict in directly.

    `.get` chains rather than subscripts throughout: `--config-path` replaces the whole
    config instead of merging it (AGENTS.md trap 2), so a caller-supplied config may
    carry no assemblies at all. Absent, the columns record as not measured - a reporting
    figure must never abort a run.
    """
    from vntyper.scripts.region_utils import resolve_assembly_alias

    assemblies = config.get("bam_processing", {}).get("assemblies", {})
    if not assemblies:
        return None
    try:
        coordinate_system = resolve_assembly_alias(reference_assembly)
    except (KeyError, ValueError):
        return assemblies.get(reference_assembly)
    return assemblies.get(coordinate_system) or assemblies.get(reference_assembly)


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
    primary_failure = False
    try:
        if region is None:
            region = region_resolver(
                bam_file=plan.view_path,
                reference_assembly=reference_assembly,
                region_type="vntr_region",
                config=config,
            )
        coverage_calculator(
            bam_file=plan.view_path,
            region=region,
            threads=threads,
            config=config,
            output_dir=output_dir,
            output_name="coverage",
            reference_path=plan.reference_path,
            index_path=plan.stable_index_path,
            assembly_config=_assembly_geometry(config, reference_assembly),
        )
    except BaseException:
        primary_failure = True
        raise
    finally:
        close_alignment_plan(plan, preserve_primary=primary_failure)
    return region
