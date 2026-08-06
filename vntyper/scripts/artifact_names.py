"""The single declaration of the basename every pipeline artefact is named from.

``vntyper pipeline`` accepts ``-n/--output-name`` and ``config.json`` carries a
``default_values.output_name``. Neither has ever reached the pipeline: ``cli.py``
resolved the value and then dropped it, because ``run_pipeline`` has no such
parameter. Every artefact is named from a literal instead.

That is not a one-line fix, and this module exists to say exactly why. Enumerating the
surface (``tests/unit/test_pipeline_artifact_paths.py``) splits it three ways:

**1. Movable.** ``process_bam_to_fastq``, ``process_fastq``, ``align_and_sort_fastq``,
``run_advntr`` and ``process_advntr_output`` each take an ``output_name`` argument.
``run_pipeline`` passes the literal ``"output"`` to all five.

**2. Rebuilt from a literal.** ``pipeline.py`` reconstructs paths rather than consuming
what a stage returned: the FASTQ pair, the sliced BAM handed to adVNTR, the sorted BAM
handed to coverage, and the adVNTR result filenames (frozen by contract C4). Those all
move together only because they all spell the same literal.

**3. Out of reach.** Three readers hardcode the basename in modules ``run_pipeline``
never passes it to:

* ``generate_report.py`` opens ``fastq_bam_processing/output.json`` for the fastp
  quality metrics;
* ``cli_report.py`` rediscovers ``kestrel/output.bam`` and ``kestrel/output.bed`` for
  the IGV panel, and ``vntyper report`` has no ``--output-name`` flag at all, so a
  custom-named run could never be re-reported;
* ``kestrel_genotyping.py`` derives ``output.sam``, ``output.bam``, ``output.bed``,
  ``output_indel.vcf``, ``output_insertion.vcf`` and ``output_deletion.vcf`` from its
  ``output_dir`` alone. It takes a ``vcf_path`` but does not name the rest from it, so
  moving ``vcf_path`` desynchronises Kestrel's own outputs from each other.

A basename that moves for group 1 and not group 3 is not a broken run. It is a
**silently negative** one: ``summary.record_step`` records a missing file with
``md5sum=None`` rather than raising, ``build_kestrel_frames`` turns that into an empty
frame, and an empty Kestrel frame is the report's configured negative default. Group 3
therefore cannot be threaded from ``run_pipeline``, and a partial threading is worse
than the defect it closes.

Until every reader takes the basename as an argument, the pipeline's basename is fixed.
:data:`PIPELINE_BASENAME` is that fixed value, :func:`pipeline_artifact_paths` is the
enumeration, and ``cli_handlers`` refuses an ``--output-name`` it cannot honour instead
of dropping it silently.

Attributes:
    PIPELINE_BASENAME: The basename every pipeline artefact is named from.
    COVERAGE_BASENAME: The separate basename the coverage stage uses.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

logger = logging.getLogger(__name__)

#: The basename every pipeline artefact is named from. Changing this changes every
#: path in :func:`pipeline_artifact_paths` *and* the literals in
#: ``generate_report.py``, ``cli_report.py`` and ``kestrel_genotyping.py``, which are
#: not derived from it. ``tests/unit/test_pipeline_artifact_paths.py`` fails if the two
#: drift apart.
PIPELINE_BASENAME = "output"

#: The coverage stage uses its own basename, so ``coverage_summary.tsv`` is stable
#: regardless of :data:`PIPELINE_BASENAME`. That file's schema is frozen by
#: contract C1 and read by ``generate_report.py``.
COVERAGE_BASENAME = "coverage"


def pipeline_artifact_paths(
    output_dir: str | os.PathLike[str],
    *,
    basename: str = PIPELINE_BASENAME,
    coverage_basename: str = COVERAGE_BASENAME,
) -> dict[str, str]:
    """Enumerate every artefact path a pipeline run derives from a basename.

    The keys are stable logical names; the values are the paths the run produces.
    A caller that wants to know "what did this run write, and where" reads this
    rather than re-deriving the layout, which is how the three-way split in the
    module docstring came about in the first place.

    Args:
        output_dir: The run's output directory.
        basename: The artefact basename. Defaults to :data:`PIPELINE_BASENAME`;
            overriding it is for tests that assert every path moves together.
        coverage_basename: The coverage stage's separate basename.

    Returns:
        dict[str, str]: Logical artefact name -> absolute path.
    """
    root = Path(output_dir)
    fastq_dir = root / "fastq_bam_processing"
    align_dir = root / "alignment_processing"
    kestrel_dir = root / "kestrel"
    advntr_dir = root / "advntr"
    coverage_dir = root / "coverage"

    return {
        # Produced by process_bam_to_fastq / process_fastq
        "fastq_r1": str(fastq_dir / f"{basename}_R1.fastq.gz"),
        "fastq_r2": str(fastq_dir / f"{basename}_R2.fastq.gz"),
        "fastq_other": str(fastq_dir / f"{basename}_other.fastq.gz"),
        "fastq_single": str(fastq_dir / f"{basename}_single.fastq.gz"),
        "sliced_bam": str(fastq_dir / f"{basename}_sliced.bam"),
        # Read back by generate_report.py, which spells this literally.
        "fastp_json": str(fastq_dir / f"{basename}.json"),
        "fastp_html": str(fastq_dir / f"{basename}.html"),
        # Produced by align_and_sort_fastq on the FASTQ path.
        "sorted_bam": str(align_dir / f"{basename}_sorted.bam"),
        # Kestrel. Only ``kestrel_vcf`` is passed in; the rest are names
        # kestrel_genotyping.py derives independently from the same literal.
        "kestrel_vcf": str(kestrel_dir / f"{basename}.vcf"),
        "kestrel_indel_vcf": str(kestrel_dir / f"{basename}_indel.vcf"),
        "kestrel_indel_vcf_gz": str(kestrel_dir / f"{basename}_indel.vcf.gz"),
        "kestrel_bam": str(kestrel_dir / f"{basename}.bam"),
        "kestrel_bed": str(kestrel_dir / f"{basename}.bed"),
        # adVNTR. Contract C4: pipeline.py rebuilds these rather than consuming
        # what the module returns, so the two must stay spelled the same.
        "advntr_tsv": str(advntr_dir / f"{basename}_adVNTR.tsv"),
        "advntr_vcf": str(advntr_dir / f"{basename}_adVNTR.vcf"),
        "advntr_result": str(advntr_dir / f"{basename}_adVNTR_result.tsv"),
        # Basename-independent, listed so the enumeration is complete.
        "coverage_summary": str(coverage_dir / f"{coverage_basename}_summary.tsv"),
        "kestrel_result": str(kestrel_dir / "kestrel_result.tsv"),
        "cross_match": str(advntr_dir / "cross_match_results.tsv"),
        "pipeline_info": str(fastq_dir / "pipeline_info.json"),
        "pipeline_summary": str(root / "pipeline_summary.json"),
    }


def validate_output_name(output_name: str | None) -> str:
    """Accept only the basename the pipeline can actually honour.

    ``--output-name`` used to be resolved and then dropped, so asking for one had no
    effect and no error. Silently ignoring it is the worse half of the defect: a user
    who asks for ``-n mysample`` and then looks for ``mysample_*`` finds nothing, and
    a caller that greps the output directory for its own basename gets an empty
    result that looks like a negative genotype.

    Args:
        output_name: The resolved value, or None to accept the default.

    Returns:
        str: The basename the pipeline will use.

    Raises:
        ValueError: If ``output_name`` is anything the pipeline cannot honour.
    """
    if output_name is None or output_name == PIPELINE_BASENAME:
        return PIPELINE_BASENAME

    msg = (
        f"--output-name {output_name!r} cannot be honoured: the pipeline's artefact basename is fixed at "
        f"{PIPELINE_BASENAME!r}. generate_report.py, cli_report.py and kestrel_genotyping.py each name their "
        "files from that literal and take no basename argument, so moving only the stages that do accept one "
        "would leave the report reading files nothing wrote - which this pipeline reports as a negative "
        "genotype rather than an error. Drop the flag, or rename the output directory instead."
    )
    logger.error(msg)
    raise ValueError(msg)
