"""Drive :func:`vntyper.scripts.pipeline.run_pipeline` with every collaborator stubbed.

``run_pipeline`` is 600 lines of orchestration wrapped around a dozen external tools.
Nothing about the *orchestration* needs a real samtools, a real Kestrel JAR or a real
FASTA - but until this harness existed there was no way to observe the orchestration
without them, which is why ``pipeline.py`` sat at 10% unit coverage (AGENTS.md rule 2).

The harness replaces every stage call with a recorder, runs the real ``run_pipeline``,
and hands back what each stage was asked to do. Two things it deliberately does **not**
stub:

* ``record_step`` / ``write_summary`` - the summary is the pipeline's own record of which
  artefact each step produced, and it is exactly what the artefact-path enumeration reads.
* ``create_output_directories`` - the real directory layout is part of every path.

``run_pipeline`` swallows every exception into ``sys.exit(1)`` (``pipeline.py``'s outer
``except Exception``), so a failure surfaces here as ``SystemExit(1)`` with the real cause
in the captured log rather than as a traceback. :func:`run_pipeline_under_harness`
re-raises the recorded exception instead, so a broken stub is visible.
"""

from __future__ import annotations

import logging
from contextlib import ExitStack
from pathlib import Path
from typing import Any
from unittest import mock

from vntyper.scripts import pipeline as pipeline_module
from vntyper.scripts import pipeline_alignment as pipeline_alignment_module
from vntyper.scripts.alignment_contract import AlignmentPlan

logger = logging.getLogger(__name__)

#: Minimal configuration ``run_pipeline`` reads. Every key here is one the pipeline
#: dereferences without ``.get`` (AGENTS.md trap 2), so dropping one is a ``KeyError``
#: deep inside the run rather than a clear failure.
MINIMAL_CONFIG: dict[str, Any] = {
    "tools": {
        "samtools": "samtools",
        "fastp": "fastp",
        "bwa": "bwa",
        "kestrel": "kestrel.jar",
        "java_path": "java",
    },
    "reference_data": {
        "muc1_reference_vntr": "/refs/muc1.fa",
        "bwa_reference_hg19": "/refs/hg19.fa",
        "bwa_reference_hg38": "/refs/hg38.fa",
        "advntr_reference_vntr_hg19": "/refs/advntr_hg19.db",
        "advntr_reference_vntr_hg38": "/refs/advntr_hg38.db",
    },
    "bam_processing": {"bam_region_hg19": "chr1:155158000-155163000"},
    "default_values": {"flanking": 50},
    "paths": {"template_dir": "vntyper/templates"},
}

#: The stage functions the harness replaces. Header reading is owned by the input
#: alignment boundary; the remaining attributes live on ``pipeline``.
PIPELINE_STAGE_ATTRS = (
    "get_tool_versions",
    "validate_bam_file",
    "validate_fastq_file",
    "get_region_string_with_fallback",
    "process_bam_to_fastq",
    "route_converted_fastqs",
    "process_fastq",
    "align_and_sort_fastq",
    "extract_bam_header",
    "read_alignment_header",
    "run_preflight",
    "parse_header_pipeline_info",
    "calculate_vntr_coverage",
    "run_kestrel",
    "generate_summary_report",
    "cross_match_variants",
    "extract_results_from_pipeline_summary",
    "write_results_tsv",
)

#: adVNTR is imported *inside* ``run_pipeline`` (an optional Python-2.7-adjacent module,
#: AGENTS.md trap 6), so it has to be patched on its own module rather than on
#: ``pipeline``'s namespace.
ADVNTR_MODULE = "vntyper.modules.advntr.advntr_genotyping"
ADVNTR_STAGE_ATTRS = ("run_advntr", "process_advntr_output", "load_advntr_config")


class PipelineHarness:
    """Recorded stage calls from one ``run_pipeline`` invocation.

    Attributes:
        stages: Stage name -> :class:`unittest.mock.MagicMock` that replaced it.
        output_dir: The directory the pipeline was pointed at.
        error: The exception ``run_pipeline`` swallowed, or None.
    """

    def __init__(self, stages: dict[str, mock.MagicMock], output_dir: Path) -> None:
        self.stages = stages
        self.output_dir = output_dir
        self.error: BaseException | None = None

    def call(self, stage: str) -> mock.MagicMock:
        """Return the mock that replaced ``stage``, asserting it was reached.

        Args:
            stage: The stage name, as listed in :data:`PIPELINE_STAGE_ATTRS`.

        Returns:
            unittest.mock.MagicMock: The replacement, already asserted as called.
        """
        stub = self.stages[stage]
        assert stub.called, f"the pipeline never reached {stage}()"
        return stub

    def kwargs(self, stage: str) -> dict[str, Any]:
        """Return the keyword arguments ``stage`` was last called with.

        Args:
            stage: The stage name.

        Returns:
            dict[str, Any]: The recorded keyword arguments.
        """
        return dict(self.call(stage).call_args.kwargs)

    def positional(self, stage: str) -> tuple:
        """Return the positional arguments ``stage`` was last called with.

        Args:
            stage: The stage name.

        Returns:
            tuple: The recorded positional arguments.
        """
        return tuple(self.call(stage).call_args.args)


def _fastq_pair(output_dir: Path, basename: str = "output") -> tuple[str, str, str, str]:
    """Return the four FASTQ paths ``process_bam_to_fastq`` reports back.

    Args:
        output_dir: The pipeline's output directory.
        basename: The artefact basename the pipeline asked for.

    Returns:
        tuple[str, str, str, str]: R1, R2, other and single FASTQ paths.
    """
    stage_dir = Path(output_dir) / "fastq_bam_processing"
    return (
        str(stage_dir / f"{basename}_R1.fastq.gz"),
        str(stage_dir / f"{basename}_R2.fastq.gz"),
        str(stage_dir / f"{basename}_other.fastq.gz"),
        str(stage_dir / f"{basename}_single.fastq.gz"),
    )


def _alignment_plan(*args: Any, **kwargs: Any) -> AlignmentPlan:
    """Return the proven plan a stubbed preflight call describes.

    Args:
        *args: Unused because ``run_pipeline`` calls preflight by keyword.
        **kwargs: The real ``run_preflight`` keyword arguments.

    Returns:
        A plan rooted at the requested output directory.
    """
    file_format = kwargs["file_format"]
    output_dir = Path(kwargs["output_dir"])
    output_name = kwargs["output_name"]
    view_path = output_dir / f"{output_name}.{file_format}"
    index_suffix = "bai" if file_format == "bam" else "crai"
    return AlignmentPlan(
        input_path=str(kwargs["in_path"]),
        view_path=str(view_path),
        file_format=file_format,
        index_path=f"{view_path}.{index_suffix}",
        reference_path="/refs/hg19.fa" if file_format == "cram" else None,
        reference_source="harness",
        uncovered_contigs=(),
        unmapped_scan="indexed",
    )


def run_pipeline_under_harness(
    output_dir: Path,
    *,
    pipeline_output_dir: str | Path | None = None,
    create_output_dir: bool = True,
    config: dict[str, Any] | None = None,
    extra_modules: list[str] | None = None,
    header: str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:249250621\n",
    header_error: BaseException | None = None,
    expect_failure: bool = False,
    stage_side_effects: dict[str, Any] | None = None,
    **run_pipeline_kwargs: Any,
) -> PipelineHarness:
    """Run the real ``run_pipeline`` with every stage replaced by a recorder.

    Args:
        output_dir: Directory the pipeline writes into; use ``tmp_path``.
        pipeline_output_dir: Optional raw output argument passed to the pipeline.
            This preserves caller spelling such as a trailing path separator while
            ``output_dir`` remains the canonical path used by the harness.
        create_output_dir: Create the output root before invoking the pipeline.
        config: Configuration mapping. Defaults to :data:`MINIMAL_CONFIG`.
        extra_modules: ``--extra-modules`` list as ``run_pipeline`` receives it.
        header: The SAM header text the header reads return. Parsed for real by the
            assembly guard, so it decides the guard's verdict.
        header_error: If given, ``extract_bam_header`` raises this and the guard's
            own read returns None, as it would in production.
        expect_failure: When True, a swallowed exception is recorded on the
            returned harness instead of being re-raised.
        stage_side_effects: Stage name -> ``side_effect`` to attach to that stage's
            recorder. Return ``unittest.mock.DEFAULT`` to keep the stubbed return
            value; anything else replaces it.
        **run_pipeline_kwargs: Forwarded to ``run_pipeline``; ``bam`` defaults to a
            BAM in a separate sibling input tree so the common case needs no arguments.

    Returns:
        PipelineHarness: The recorded stage calls.

    Raises:
        AssertionError: If ``run_pipeline`` failed and ``expect_failure`` is False.
    """
    output_dir = Path(output_dir)
    if create_output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
    resolved_config = MINIMAL_CONFIG if config is None else config
    basename = "output"

    kwargs: dict[str, Any] = {
        "bwa_reference": "/refs/hg19.fa",
        "output_dir": output_dir if pipeline_output_dir is None else pipeline_output_dir,
        "extra_modules": extra_modules or [],
        "module_args": {"advntr": {}},
        "config": resolved_config,
        "reference_assembly": "hg19",
        "sample_name": "sample",
    }
    kwargs.update(run_pipeline_kwargs)
    if not any(kwargs.get(key) for key in ("bam", "cram", "fastq1", "fastq2")):
        input_root = output_dir.parent / f"{output_dir.name}_input"
        input_root.mkdir()
        bam_path = input_root / "in.bam"
        bam_path.touch()
        kwargs["bam"] = str(bam_path)

    # The adVNTR branch reads back the sliced BAM that the (stubbed) BAM-to-FASTQ stage
    # would have written, and aborts if it is missing. Stand it in so the branch is
    # reachable at all.
    if "advntr" in (extra_modules or []):
        sliced = output_dir / "fastq_bam_processing" / f"{basename}_sliced.bam"
        sliced.parent.mkdir(parents=True, exist_ok=True)
        sliced.touch()

    stages: dict[str, mock.MagicMock] = {}
    with ExitStack() as stack:
        for attr in PIPELINE_STAGE_ATTRS:
            owner = pipeline_alignment_module if attr == "read_alignment_header" else pipeline_module
            stages[attr] = stack.enter_context(mock.patch.object(owner, attr, autospec=True))
        stack.enter_context(
            mock.patch.object(
                pipeline_alignment_module,
                "get_region_string_with_fallback",
                stages["get_region_string_with_fallback"],
            )
        )
        stack.enter_context(mock.patch.object(pipeline_alignment_module, "run_preflight", stages["run_preflight"]))
        stack.enter_context(
            mock.patch.object(
                pipeline_alignment_module,
                "pin_reference_resolution",
                pipeline_module.pin_reference_resolution,
            )
        )
        stack.enter_context(
            mock.patch.object(
                pipeline_alignment_module,
                "restore_reference_resolution",
                pipeline_module.restore_reference_resolution,
            )
        )
        for attr in ADVNTR_STAGE_ATTRS:
            stages[attr] = stack.enter_context(mock.patch(f"{ADVNTR_MODULE}.{attr}", autospec=True))

        stages["get_tool_versions"].return_value = {"samtools": "1.19"}
        stages["get_region_string_with_fallback"].return_value = "chr1:155158000-155163000"
        stages["run_preflight"].side_effect = _alignment_plan
        stages["process_bam_to_fastq"].return_value = _fastq_pair(output_dir, basename)
        stages["route_converted_fastqs"].return_value = _fastq_pair(output_dir, basename)[:2]
        stages["align_and_sort_fastq"].return_value = str(
            output_dir / "alignment_processing" / f"{basename}_sorted.bam"
        )
        # The assembly guard's own samtools call. It is stubbed rather than the
        # verdict, so ``reconcile_assembly`` and the contig parser run for real
        # against the header text a test supplies.
        if header_error is not None:
            stages["extract_bam_header"].side_effect = header_error
            stages["read_alignment_header"].return_value = None
        else:
            stages["extract_bam_header"].return_value = header
            stages["read_alignment_header"].return_value = header
        stages["calculate_vntr_coverage"].return_value = {"mean": 100.0}
        stages["load_advntr_config"].return_value = {"advntr_settings": {"output_format": "tsv"}}
        stages["extract_results_from_pipeline_summary"].return_value = ([{"a": 1}], [{"b": 2}])
        stages["cross_match_variants"].return_value = {"matches": [], "overall_match": True}

        for stage_name, side_effect in (stage_side_effects or {}).items():
            stages[stage_name].side_effect = side_effect

        harness = PipelineHarness(stages, output_dir)
        try:
            pipeline_module.run_pipeline(**kwargs)
        except SystemExit as exc:  # run_pipeline's outer handler
            harness.error = exc
            if not expect_failure:
                raise AssertionError(
                    "run_pipeline exited instead of completing; the real cause is in the captured log"
                ) from exc
        except BaseException as exc:  # a stage raised before the outer handler
            harness.error = exc
            if not expect_failure:
                raise

    return harness


def artifact_paths_from_summary(output_dir: Path) -> set[str]:
    """Return every ``result_file`` the run recorded, relative to ``output_dir``.

    The summary is the pipeline's own claim about what it produced, and it is what
    ``generate_report.py`` reads back. A path recorded here that no stage actually
    wrote is a silent empty section in the report, never an error (contract C5).

    Args:
        output_dir: The directory the pipeline ran in.

    Returns:
        set[str]: POSIX-style relative paths.
    """
    import json

    summary_path = Path(output_dir) / "pipeline_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    return {
        Path(step["result_file"]).resolve().relative_to(Path(output_dir).resolve()).as_posix()
        for step in summary["steps"]
    }
