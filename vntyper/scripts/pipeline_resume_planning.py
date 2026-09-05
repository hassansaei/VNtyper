"""Pure resume planning and compatibility evaluation.

Extracted from pipeline.py (#20).
"""

from __future__ import annotations

import logging
import os
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

from vntyper.scripts import summary_steps
from vntyper.scripts.resume import (
    caller_advntr_matches,
    caller_kestrel_matches,
    caller_shark_matches,
    fingerprint_file,
    fingerprint_runtime,
    make_reused_step_record,
    reference_content_matches,
    step_is_reusable,
)
from vntyper.scripts.run_configuration import RunConfiguration

logger = logging.getLogger(__name__)


def _record_reused_step(summary: dict[str, Any], record: dict[str, Any]) -> None:
    st_name = record.get("step")
    idx = next((i for i, s in enumerate(summary["steps"]) if s.get("step") == st_name), None)
    if idx is not None:
        summary["steps"][idx] = record
    else:
        summary["steps"].append(record)


def resolve_effective_tool_identity(
    raw_command: str | None,
) -> dict[str, str | None]:
    """Resolve command, executable path, and content fingerprint for a tool command."""
    if not raw_command:
        return {"command": None, "executable": None, "fingerprint": None}
    cmd = str(raw_command).strip()
    exe_path: str | None = None
    p = Path(cmd)
    if p.is_file():
        exe_path = str(p.resolve())
    else:
        which = shutil.which(cmd)
        if which:
            exe_path = str(Path(which).resolve())
    fp = fingerprint_file(Path(exe_path)) if exe_path and Path(exe_path).is_file() else None
    return {
        "command": cmd,
        "executable": exe_path,
        "fingerprint": fp,
    }


def resolve_effective_preprocessing_tools(
    config: Mapping[str, Any] | None,
    input_type: str | None,
    extra_modules: tuple[str, ...] | Sequence[str] = (),
) -> dict[str, Any] | None:
    """Resolve effective tool identities for FASTQ preprocessing (fastp, bwa, and optionally shark)."""
    if input_type != "FASTQ" or not isinstance(config, Mapping):
        return None
    tools = config.get("tools", {})
    if not isinstance(tools, Mapping):
        return None
    res: dict[str, Any] = {
        "fastp": resolve_effective_tool_identity(tools.get("fastp")),
        "bwa": resolve_effective_tool_identity(tools.get("bwa")),
    }
    if "shark" in extra_modules:
        res["shark"] = resolve_effective_tool_identity(tools.get("shark", "shark"))
    return res


def build_analysis_settings(
    *,
    reference_assembly: str | None,
    fast_mode: bool,
    custom_regions: str | None,
    bed_file: str | None,
    advntr_reference: str | None,
    module_args: Mapping[str, Any] | None,
    config: Mapping[str, Any] | None,
    extra_modules: tuple[str, ...],
    input_type: str | None = None,
) -> dict[str, Any]:
    """Build canonical analysis settings dictionary for run identity tracking."""
    advntr_max_coverage = None
    if isinstance(module_args, dict) and "advntr" in module_args and isinstance(module_args["advntr"], dict):
        advntr_max_coverage = module_args["advntr"].get("max_coverage")

    bam_processing_settings = None
    if isinstance(config, dict) and "bam_processing" in config and isinstance(config["bam_processing"], dict):
        bam_processing_settings = dict(config["bam_processing"])

    preprocessing_tools = resolve_effective_preprocessing_tools(config, input_type, extra_modules)

    return {
        "reference_assembly": reference_assembly,
        "fast_mode": bool(fast_mode),
        "custom_regions": str(custom_regions) if custom_regions else None,
        "bed_file": str(Path(bed_file).resolve()) if bed_file else None,
        "advntr_model": str(Path(advntr_reference).resolve()) if advntr_reference else None,
        "advntr_max_coverage": advntr_max_coverage,
        "bam_processing": bam_processing_settings,
        "extra_modules": sorted(extra_modules),
        "preprocessing_tools": preprocessing_tools,
    }


def build_canonical_inputs_and_fingerprints(
    input_type: str,
    fastq1: str | None,
    fastq2: str | None,
    bam: str | None,
    cram: str | None,
    bed_file: str | None,
) -> tuple[dict[str, str], dict[str, str]]:
    """Resolve absolute canonical input paths and calculate their content fingerprints."""
    canonical_input_files: dict[str, str] = {}
    input_fingerprints: dict[str, str] = {}
    if input_type == "BAM" and bam is not None:
        p = Path(bam).resolve()
        canonical_input_files["bam"] = str(p)
        if p.is_file():
            input_fingerprints["bam"] = fingerprint_file(p)
    elif input_type == "CRAM" and cram is not None:
        p = Path(cram).resolve()
        canonical_input_files["cram"] = str(p)
        if p.is_file():
            input_fingerprints["cram"] = fingerprint_file(p)
    elif input_type == "FASTQ" and fastq1 is not None:
        f1_path = Path(fastq1).resolve()
        canonical_input_files["fastq1"] = str(f1_path)
        if f1_path.is_file():
            input_fingerprints["fastq1"] = fingerprint_file(f1_path)
        if fastq2 is not None:
            f2_path = Path(fastq2).resolve()
            canonical_input_files["fastq2"] = str(f2_path)
            if f2_path.is_file():
                input_fingerprints["fastq2"] = fingerprint_file(f2_path)

    if bed_file is not None:
        bed_path = Path(bed_file).resolve()
        if bed_path.is_file():
            input_fingerprints["bed_file"] = fingerprint_file(bed_path)

    return canonical_input_files, input_fingerprints


def resolve_effective_kestrel_runtime(
    run_configuration: RunConfiguration,
    config: Mapping[str, Any],
    project_root: str,
) -> tuple[str, dict[str, Any], str]:
    """Resolve effective Kestrel runtime mapping and fingerprint including executable identity."""
    from vntyper.scripts.kestrel_counting import DEFAULT_KANALYZE_PATH
    from vntyper.scripts.pipeline_kestrel import resolve_kestrel_counting_mode

    kestrel_counting_mode = resolve_kestrel_counting_mode(run_configuration.kestrel_runtime, config)
    raw_kestrel = config.get("tools", {}).get("kestrel")
    kestrel_path = (
        str(Path(os.path.join(project_root, raw_kestrel)).resolve())
        if raw_kestrel and not os.path.isabs(raw_kestrel)
        else (str(Path(raw_kestrel).resolve()) if raw_kestrel else None)
    )
    kestrel_fp = fingerprint_file(Path(kestrel_path)) if kestrel_path and Path(kestrel_path).is_file() else None
    raw_kanalyze = config.get("tools", {}).get("kanalyze", DEFAULT_KANALYZE_PATH)
    kanalyze_path = (
        str(Path(os.path.join(project_root, raw_kanalyze)).resolve())
        if raw_kanalyze and not os.path.isabs(raw_kanalyze)
        else (str(Path(raw_kanalyze).resolve()) if raw_kanalyze else None)
    )
    kanalyze_fp = (
        fingerprint_file(Path(kanalyze_path))
        if kanalyze_path and Path(kanalyze_path).is_file() and kestrel_counting_mode == "split"
        else None
    )
    effective_kestrel_runtime = {
        **dict(run_configuration.kestrel_runtime),
        "kestrel_counting_mode": kestrel_counting_mode,
        "kanalyze": kanalyze_path,
        "kanalyze_fingerprint": kanalyze_fp,
        "kestrel_executable": kestrel_path,
        "kestrel_executable_fingerprint": kestrel_fp,
    }
    kestrel_runtime_fingerprint = fingerprint_runtime(effective_kestrel_runtime)
    return kestrel_counting_mode, effective_kestrel_runtime, cast(str, kestrel_runtime_fingerprint)


def resolve_effective_advntr_runtime(
    run_configuration: RunConfiguration,
    config: Mapping[str, Any],
    advntr_version: str | tuple[int, int, int] | None = None,
) -> tuple[dict[str, Any], str]:
    """Resolve effective adVNTR runtime mapping and fingerprint including command and tool identity."""
    raw_advntr = config.get("tools", {}).get("advntr")
    advntr_command = str(raw_advntr).strip() if raw_advntr is not None else None
    advntr_fp = None
    if advntr_command:
        p = Path(advntr_command)
        if p.is_file():
            advntr_fp = fingerprint_file(p)

    ver_str = None
    if advntr_version is not None:
        ver_str = (
            ".".join(str(part) for part in advntr_version) if isinstance(advntr_version, tuple) else str(advntr_version)
        )

    effective_advntr_runtime = {
        **dict(run_configuration.advntr_runtime),
        "advntr_command": advntr_command,
        "advntr_command_fingerprint": advntr_fp,
        "advntr_version": ver_str,
    }
    advntr_runtime_fingerprint = fingerprint_runtime(effective_advntr_runtime)
    return effective_advntr_runtime, cast(str, advntr_runtime_fingerprint)


def resolve_effective_shark_runtime(
    run_configuration: RunConfiguration,
    config: Mapping[str, Any],
) -> tuple[dict[str, Any], str]:
    """Resolve effective SHARK runtime mapping and fingerprint including command and tool identity."""
    raw_shark = config.get("tools", {}).get("shark", "shark") if isinstance(config.get("tools"), Mapping) else "shark"
    shark_tool = resolve_effective_tool_identity(raw_shark)
    effective_shark_runtime = {
        **dict(run_configuration.shark_runtime),
        "shark_command": shark_tool["command"],
        "shark_executable": shark_tool["executable"],
        "shark_executable_fingerprint": shark_tool["fingerprint"],
    }
    shark_runtime_fingerprint = fingerprint_runtime(effective_shark_runtime)
    return effective_shark_runtime, cast(str, shark_runtime_fingerprint)


@dataclass(frozen=True)
class ResumeCompatibility:
    """Structured compatibility assessment across pipeline stages."""

    kestrel_ref_matches: bool
    advntr_model_matches: bool
    shark_ref_matches: bool
    bwa_ref_matches: bool
    cram_ref_matches: bool
    inval_align: bool
    inval_cram: bool
    inval_qc: bool = False


def evaluate_resume_compatibility(
    prior_summary: Mapping[str, Any] | None,
    *,
    input_type: str,
    kestrel_reference_path: str | None,
    kestrel_reference_fingerprint: str | None,
    kestrel_motifs_path: str | None,
    kestrel_motifs_fingerprint: str | None,
    kestrel_runtime_fingerprint: str | None,
    kestrel_counting_mode: str | None,
    advntr_model_sha: str | None,
    advntr_rus_path: str | None,
    advntr_rus_fingerprint: str | None,
    advntr_runtime_fingerprint: str | None,
    shark_reference_path: str | None,
    shark_reference_fingerprint: str | None,
    shark_runtime_fingerprint: str | None,
    effective_reference_path: str | None,
    effective_reference_fingerprint: str | None,
    advntr_version: str | None = None,
    current_preprocessing_tools: Mapping[str, Any] | None = None,
) -> ResumeCompatibility:
    """Evaluate whether callers, references, and alignments match prior checkpoint."""
    kestrel_ref_matches = caller_kestrel_matches(
        prior_summary,
        kestrel_reference_path=kestrel_reference_path,
        kestrel_reference_fingerprint=kestrel_reference_fingerprint,
        kestrel_motifs_path=kestrel_motifs_path,
        kestrel_motifs_fingerprint=kestrel_motifs_fingerprint,
        kestrel_runtime_fingerprint=kestrel_runtime_fingerprint,
        kestrel_counting_mode=kestrel_counting_mode,
    )

    advntr_model_matches = caller_advntr_matches(
        prior_summary,
        curr_model_sha=advntr_model_sha,
        advntr_rus_path=advntr_rus_path,
        advntr_rus_fingerprint=advntr_rus_fingerprint,
        advntr_runtime_fingerprint=advntr_runtime_fingerprint,
        advntr_version=advntr_version,
    )

    shark_ref_matches = caller_shark_matches(
        prior_summary,
        shark_reference_path=shark_reference_path,
        shark_reference_fingerprint=shark_reference_fingerprint,
        shark_runtime_fingerprint=shark_runtime_fingerprint,
    )
    bwa_ref_matches = reference_content_matches(
        prior_summary,
        reference_fingerprint=effective_reference_fingerprint,
    )
    cram_ref_matches = True
    if input_type == "CRAM" and prior_summary is not None:
        prior_ref_fp = prior_summary.get("reference_fingerprint")
        prior_ref_path = prior_summary.get("persistent_reference_path") or prior_summary.get("reference_path")
        if (
            (prior_ref_fp is not None or effective_reference_fingerprint is not None)
            and prior_ref_fp != effective_reference_fingerprint
            or (
                prior_ref_path is not None
                and effective_reference_path is not None
                and str(Path(prior_ref_path).resolve()) != str(Path(effective_reference_path).resolve())
            )
        ):
            cram_ref_matches = False

    fastp_matches = True
    bwa_tool_matches = True
    shark_tool_matches = True
    if input_type == "FASTQ" and prior_summary is not None:
        prior_tools = prior_summary.get("analysis_settings", {}).get("preprocessing_tools")
        if prior_tools is not None and current_preprocessing_tools is not None:
            fastp_matches = prior_tools.get("fastp") == current_preprocessing_tools.get("fastp")
            bwa_tool_matches = prior_tools.get("bwa") == current_preprocessing_tools.get("bwa")
            if "shark" in prior_tools or "shark" in current_preprocessing_tools:
                shark_tool_matches = prior_tools.get("shark") == current_preprocessing_tools.get("shark")

    inval_qc = not fastp_matches or not shark_tool_matches
    inval_align = (
        not shark_ref_matches
        or not shark_tool_matches
        or not bwa_ref_matches
        or not fastp_matches
        or not bwa_tool_matches
    )
    inval_cram = not cram_ref_matches
    if input_type == "FASTQ" and inval_align:
        kestrel_ref_matches = False
        advntr_model_matches = False
    if input_type == "CRAM" and inval_cram:
        kestrel_ref_matches = False
        advntr_model_matches = False

    return ResumeCompatibility(
        kestrel_ref_matches=kestrel_ref_matches,
        advntr_model_matches=advntr_model_matches,
        shark_ref_matches=shark_ref_matches,
        bwa_ref_matches=bwa_ref_matches,
        cram_ref_matches=cram_ref_matches,
        inval_align=inval_align,
        inval_cram=inval_cram,
        inval_qc=inval_qc,
    )


def record_reused_stage(
    summary: dict[str, Any],
    prior_summary: dict[str, Any],
    step_name: str,
) -> None:
    """Record one reused stage record and carry forward its artifact MD5s."""
    prior_step = next((s for s in prior_summary.get("steps", []) if s.get("step") == step_name), None)
    if prior_step is not None:
        _record_reused_step(summary, make_reused_step_record(prior_step, prior_summary.get("pipeline_start")))
    if prior_summary.get("stage_artifact_md5s", {}).get(step_name):
        summary.setdefault("stage_artifact_md5s", {})[step_name] = dict(prior_summary["stage_artifact_md5s"][step_name])


def initial_stage_carry_forward(
    summary: dict[str, Any],
    prior_summary: dict[str, Any],
    output_dir: str | Path,
    *,
    compatibility: ResumeCompatibility,
    needs_advntr: bool,
    extra_modules: tuple[str, ...],
) -> None:
    """Populate reused steps and stage artifact hashes from prior summary before stage execution."""
    for s in prior_summary.get("steps", []):
        st = s.get("step")
        if not st or not step_is_reusable(prior_summary, st, output_dir, needs_advntr=needs_advntr):
            continue
        if st == summary_steps.STEP_FASTQ_QC and compatibility.inval_qc:
            continue
        if (
            st
            in (
                summary_steps.STEP_FASTQ_ALIGNMENT,
                summary_steps.STEP_BAM_TO_FASTQ_POST_ALIGNMENT,
                summary_steps.STEP_SHARK,
            )
            and compatibility.inval_align
        ):
            continue
        if st == summary_steps.STEP_CRAM_TO_FASTQ and compatibility.inval_cram:
            continue
        if st == summary_steps.STEP_KESTREL and not compatibility.kestrel_ref_matches:
            continue
        if st == summary_steps.STEP_ADVNTR and (
            "advntr" not in extra_modules or not compatibility.advntr_model_matches
        ):
            continue
        if st == summary_steps.STEP_CROSS_MATCH and (
            not compatibility.kestrel_ref_matches
            or "advntr" not in extra_modules
            or not compatibility.advntr_model_matches
        ):
            continue
        _record_reused_step(summary, make_reused_step_record(s, prior_summary.get("pipeline_start")))

    if (
        compatibility.advntr_model_matches
        and any(s.get("step") == summary_steps.STEP_ADVNTR for s in summary.get("steps", []))
        and prior_summary.get("advntr_model")
    ):
        summary["advntr_model"] = dict(prior_summary["advntr_model"])

    if prior_summary.get("stage_artifact_md5s"):
        for st, md5s in prior_summary["stage_artifact_md5s"].items():
            if st == summary_steps.STEP_FASTQ_QC and compatibility.inval_qc:
                continue
            if (
                st
                in (
                    summary_steps.STEP_FASTQ_ALIGNMENT,
                    summary_steps.STEP_BAM_TO_FASTQ_POST_ALIGNMENT,
                    summary_steps.STEP_SHARK,
                )
                and compatibility.inval_align
            ):
                continue
            if st == summary_steps.STEP_CRAM_TO_FASTQ and compatibility.inval_cram:
                continue
            if st == summary_steps.STEP_KESTREL and not compatibility.kestrel_ref_matches:
                continue
            if st == summary_steps.STEP_ADVNTR and (
                "advntr" not in extra_modules or not compatibility.advntr_model_matches
            ):
                continue
            if st == summary_steps.STEP_CROSS_MATCH and (
                not compatibility.kestrel_ref_matches
                or "advntr" not in extra_modules
                or not compatibility.advntr_model_matches
            ):
                continue
            summary.setdefault("stage_artifact_md5s", {})[st] = dict(md5s)
