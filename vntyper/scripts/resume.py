"""
vntyper/scripts/resume.py

Module Purpose:
---------------
Pure logic and checkpoint inspection for pipeline resumption (--resume, #20).

Provides:
- load_prior_summary: Load and validate an existing pipeline_summary.json.
- resume_refusals: Enforce fatal run-identity invariant checks.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import logging
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Final

from vntyper.scripts import summary_steps

logger = logging.getLogger(__name__)

#: Table of required extra sibling files per reusable stage
STEP_OUTPUT_SIBLINGS: Final[dict[str, tuple[str, ...]]] = {
    summary_steps.STEP_KESTREL: (
        "output.vcf",
        "output_indel.vcf",
        "output.bam",
        "kestrel_pre_result.tsv",
    ),
    summary_steps.STEP_ADVNTR: (),
    summary_steps.STEP_BAM_TO_FASTQ: (
        "output_sliced.bam",
        "output_R1.fastq.gz",
        "output_R2.fastq.gz",
        "output_single.fastq.gz",
        "output_other.fastq.gz",
    ),
    summary_steps.STEP_CRAM_TO_FASTQ: (
        "output_sliced.bam",
        "output_R1.fastq.gz",
        "output_R2.fastq.gz",
        "output_single.fastq.gz",
        "output_other.fastq.gz",
    ),
    summary_steps.STEP_BAM_TO_FASTQ_POST_ALIGNMENT: (
        "output_sliced.bam",
        "output_R1.fastq.gz",
        "output_R2.fastq.gz",
        "output_single.fastq.gz",
        "output_other.fastq.gz",
    ),
    summary_steps.STEP_FASTQ_ALIGNMENT: (),
    summary_steps.STEP_FASTQ_QC: (),
}

_REUSABLE_STEPS: Final[frozenset[str]] = frozenset(STEP_OUTPUT_SIBLINGS.keys())


def fingerprint_file(path: str | Path) -> str:
    """Compute a lightweight, deterministic content fingerprint of an input file.

    Combines file size with SHA-256 of the initial 64 KB and trailing 64 KB of the
    file. Fast on multi-gigabyte BAM/CRAM files while detecting in-place modifications
    or content replacements.
    """
    p = Path(path)
    size = p.stat().st_size
    hasher = hashlib.sha256()
    hasher.update(str(size).encode("utf-8"))
    with p.open("rb") as handle:
        head = handle.read(65536)
        hasher.update(head)
        if size > 65536:
            handle.seek(max(0, size - 65536))
            tail = handle.read(65536)
            hasher.update(tail)
    return hasher.hexdigest()


def fingerprint_runtime(settings: Mapping[str, Any] | None) -> str | None:
    """Compute canonical sha256 digest of thawed runtime configuration settings."""
    if settings is None:
        return None
    from vntyper.scripts.canonical_json import canonical_sha256
    from vntyper.scripts.decision_profile import _thaw

    return canonical_sha256(_thaw(settings))


def caller_kestrel_matches(
    prior_summary: Mapping[str, Any] | None,
    *,
    kestrel_reference_path: str | None = None,
    kestrel_reference_fingerprint: str | None = None,
    kestrel_motifs_path: str | None = None,
    kestrel_motifs_fingerprint: str | None = None,
    kestrel_runtime_fingerprint: str | None = None,
) -> bool:
    """Return whether prior summary Kestrel reference, motifs, and runtime match current run."""
    if not prior_summary:
        return True

    prior_k_ref = prior_summary.get("kestrel_reference_path")
    if (prior_k_ref is None) != (kestrel_reference_path is None) or (
        prior_k_ref is not None
        and kestrel_reference_path is not None
        and str(Path(prior_k_ref).resolve()) != str(Path(kestrel_reference_path).resolve())
    ):
        return False

    prior_k_fp = prior_summary.get("kestrel_reference_fingerprint")
    if (
        prior_k_fp is not None or kestrel_reference_fingerprint is not None
    ) and prior_k_fp != kestrel_reference_fingerprint:
        return False

    prior_k_motifs = prior_summary.get("kestrel_motifs_path")
    if (prior_k_motifs is None) != (kestrel_motifs_path is None) or (
        prior_k_motifs is not None
        and kestrel_motifs_path is not None
        and str(Path(prior_k_motifs).resolve()) != str(Path(kestrel_motifs_path).resolve())
    ):
        return False

    prior_k_motifs_fp = prior_summary.get("kestrel_motifs_fingerprint")
    if (
        prior_k_motifs_fp is not None or kestrel_motifs_fingerprint is not None
    ) and prior_k_motifs_fp != kestrel_motifs_fingerprint:
        return False

    prior_k_rt_fp = prior_summary.get("kestrel_runtime_fingerprint")
    return prior_k_rt_fp == kestrel_runtime_fingerprint


def caller_advntr_matches(
    prior_summary: Mapping[str, Any] | None,
    *,
    curr_model_sha: str | None = None,
    advntr_rus_path: str | None = None,
    advntr_rus_fingerprint: str | None = None,
    advntr_runtime_fingerprint: str | None = None,
) -> bool:
    """Return whether prior summary adVNTR model, repeat units, and runtime match current run."""
    if not prior_summary:
        return True

    prior_has_advntr = any(s.get("step") == summary_steps.STEP_ADVNTR for s in prior_summary.get("steps", []))
    if prior_has_advntr:
        prior_model_sha = prior_summary.get("advntr_model", {}).get("sha256")
        if not prior_model_sha or not curr_model_sha or prior_model_sha != curr_model_sha:
            return False

    prior_rus_path = prior_summary.get("advntr_rus_path")
    if (prior_rus_path is None) != (advntr_rus_path is None) or (
        prior_rus_path is not None
        and advntr_rus_path is not None
        and str(Path(prior_rus_path).resolve()) != str(Path(advntr_rus_path).resolve())
    ):
        return False

    prior_rus_fp = prior_summary.get("advntr_rus_fingerprint")
    if (prior_rus_fp is not None or advntr_rus_fingerprint is not None) and prior_rus_fp != advntr_rus_fingerprint:
        return False

    prior_adv_rt_fp = prior_summary.get("advntr_runtime_fingerprint")
    return prior_adv_rt_fp == advntr_runtime_fingerprint


def _compute_md5(path: Path) -> str | None:
    """Compute the MD5 checksum of a file in 64 KiB chunks."""
    if not path.is_file():
        return None
    md5 = hashlib.md5()
    try:
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                md5.update(chunk)
        return md5.hexdigest()
    except OSError as exc:
        logger.warning("Error calculating MD5 for %s: %s", path, exc)
        return None


def load_prior_summary(path: str | Path) -> dict[str, Any] | None:
    """Load a prior pipeline summary from a file or directory path.

    Args:
        path: Path to pipeline_summary.json or the run directory containing it.

    Returns:
        dict[str, Any] | None: The parsed JSON dictionary, or None if missing/invalid.
    """
    p = Path(path)
    if p.is_dir():
        p = p / "pipeline_summary.json"

    data = None
    if p.is_file():
        try:
            parsed = json.loads(p.read_text(encoding="utf-8"))
            if isinstance(parsed, dict):
                data = parsed
            else:
                logger.warning("Prior summary at %s is not a dictionary", p)
        except (json.JSONDecodeError, OSError) as exc:
            logger.warning("Failed to parse prior summary at %s: %s", p, exc)

    donor_path = p.parent / "pipeline_summary.donor.json"
    if donor_path.is_file():
        with contextlib.suppress(json.JSONDecodeError, OSError):
            donor_data = json.loads(donor_path.read_text(encoding="utf-8"))
            if isinstance(donor_data, dict):
                if data is None:
                    data = donor_data
                else:
                    donor_profile = (
                        donor_data.get("decision_profile", {}).get("digest")
                        if isinstance(donor_data.get("decision_profile"), dict)
                        else None
                    )
                    data_profile = (
                        data.get("decision_profile", {}).get("digest")
                        if isinstance(data.get("decision_profile"), dict)
                        else None
                    )
                    donor_matches = (
                        donor_data.get("sample_name") == data.get("sample_name")
                        and donor_data.get("input_files") == data.get("input_files")
                        and donor_data.get("canonical_input_files") == data.get("canonical_input_files")
                        and donor_data.get("reference_assembly_requested") == data.get("reference_assembly_requested")
                        and donor_profile == data_profile
                    )
                    if not donor_matches:
                        logger.warning("Ignoring incompatible donor checkpoint at %s", donor_path)
                    else:
                        k_ref_matches = caller_kestrel_matches(
                            donor_data,
                            kestrel_reference_path=data.get("kestrel_reference_path"),
                            kestrel_reference_fingerprint=data.get("kestrel_reference_fingerprint"),
                            kestrel_motifs_path=data.get("kestrel_motifs_path"),
                            kestrel_motifs_fingerprint=data.get("kestrel_motifs_fingerprint"),
                            kestrel_runtime_fingerprint=data.get("kestrel_runtime_fingerprint"),
                        )
                        adv_model_matches = caller_advntr_matches(
                            donor_data,
                            curr_model_sha=data.get("advntr_model", {}).get("sha256"),
                            advntr_rus_path=data.get("advntr_rus_path"),
                            advntr_rus_fingerprint=data.get("advntr_rus_fingerprint"),
                            advntr_runtime_fingerprint=data.get("advntr_runtime_fingerprint"),
                        )
                        existing_steps = {s.get("step") for s in data.get("steps", [])}
                        for s in donor_data.get("steps", []):
                            st = s.get("step")
                            if not st or st in existing_steps:
                                continue
                            if st == summary_steps.STEP_KESTREL and not k_ref_matches:
                                continue
                            if st == summary_steps.STEP_ADVNTR and not adv_model_matches:
                                continue
                            if st == summary_steps.STEP_CROSS_MATCH and (not k_ref_matches or not adv_model_matches):
                                continue
                            data.setdefault("steps", []).append(s)
                        if donor_data.get("stage_artifact_md5s"):
                            for st, md5s in donor_data["stage_artifact_md5s"].items():
                                if st == summary_steps.STEP_KESTREL and not k_ref_matches:
                                    continue
                                if st == summary_steps.STEP_ADVNTR and not adv_model_matches:
                                    continue
                                if st == summary_steps.STEP_CROSS_MATCH and (
                                    not k_ref_matches or not adv_model_matches
                                ):
                                    continue
                                data.setdefault("stage_artifact_md5s", {}).setdefault(st, md5s)

    if data is None:
        logger.debug("Prior summary not found at %s", p)
    return data


def resume_refusals(
    prior: dict[str, Any],
    *,
    version: str,
    input_files: dict[str, Any],
    sample_name: str,
    reference_key_used: str | None,
    decision_profile_sha256: str,
    canonical_input_files: dict[str, str] | None = None,
    reference_assembly: str | None = None,
    analysis_settings: dict[str, Any] | None = None,
    reference_path: str | None = None,
    input_fingerprints: dict[str, str] | None = None,
) -> list[str]:
    """Validate that run identity invariants match the prior summary.

    Refusals concern run identity only: version, input files, sample name,
    reference key, decision-profile digest, reference assembly, analysis settings,
    and effective reference path.

    Args:
        prior: Prior summary dictionary.
        version: Pipeline version for the current run.
        input_files: Input file paths mapping for the current run (display basenames).
        sample_name: Resolved sample name for the current run.
        reference_key_used: Reference key used for the current run.
        decision_profile_sha256: SHA256 digest of the resolved decision profile.
        canonical_input_files: Optional canonical resolved input file paths mapping.
        reference_assembly: Optional requested reference assembly.
        analysis_settings: Optional result-affecting analysis settings dictionary.
        reference_path: Optional effective canonical reference path for alignment or decoding.

    Returns:
        list[str]: Descriptions of all detected mismatches; empty if allowed.
    """
    refusals: list[str] = []

    prior_version = prior.get("version")
    if prior_version != version:
        refusals.append(f"version differs (prior: {prior_version!r}, current: {version!r})")

    prior_inputs = prior.get("input_files")
    if prior_inputs != input_files:
        refusals.append(f"input files differ (prior: {prior_inputs!r}, current: {input_files!r})")

    if canonical_input_files is not None:
        prior_canonical = prior.get("canonical_input_files")
        if prior_canonical is None:
            refusals.append("checkpoint lacks canonical input identities (legacy or incomplete checkpoint)")
        elif prior_canonical != canonical_input_files:
            refusals.append(
                f"canonical input files differ (prior: {prior_canonical!r}, current: {canonical_input_files!r})"
            )

    if input_fingerprints is not None:
        prior_fingerprints = prior.get("input_fingerprints")
        if prior_fingerprints is None:
            refusals.append("checkpoint lacks input content fingerprints (legacy checkpoint)")
        elif prior_fingerprints != input_fingerprints:
            refusals.append(
                f"input file contents differ from checkpoint (prior: {prior_fingerprints!r}, current: {input_fingerprints!r})"
            )

    prior_sample = prior.get("sample_name")
    if prior_sample != sample_name:
        refusals.append(f"sample name differs (prior: {prior_sample!r}, current: {sample_name!r})")

    if reference_assembly is not None:
        prior_assembly = prior.get("reference_assembly_requested")
        if prior_assembly is not None and prior_assembly != reference_assembly:
            refusals.append(f"reference assembly differs (prior: {prior_assembly!r}, current: {reference_assembly!r})")

    prior_ref_key = prior.get("reference_key_used")
    if prior_ref_key != reference_key_used:
        refusals.append(f"reference key differs (prior: {prior_ref_key!r}, current: {reference_key_used!r})")

    if reference_path is not None:
        prior_ref_path = prior.get("persistent_reference_path") or prior.get("reference_path")
        if prior_ref_path is not None:
            prior_resolved = str(Path(prior_ref_path).resolve()) if prior_ref_path else None
            curr_resolved = str(Path(reference_path).resolve()) if reference_path else None
            if prior_resolved != curr_resolved:
                refusals.append(f"reference path differs (prior: {prior_ref_path!r}, current: {reference_path!r})")

    prior_profile_sha = prior.get("decision_profile_sha256")
    if prior_profile_sha != decision_profile_sha256:
        refusals.append(
            f"decision profile digest differs (prior: {prior_profile_sha!r}, current: {decision_profile_sha256!r})"
        )

    if analysis_settings is not None:
        prior_settings = prior.get("analysis_settings")
        if prior_settings is None:
            refusals.append("checkpoint lacks recorded analysis settings")
        elif prior_settings != analysis_settings:
            for k, v in analysis_settings.items():
                if prior_settings.get(k) != v:
                    refusals.append(
                        f"analysis setting {k!r} differs (prior: {prior_settings.get(k)!r}, current: {v!r})"
                    )

    return refusals


def step_is_reusable(
    prior: dict[str, Any],
    step_name: str,
    output_root: str | Path,
) -> bool:
    """Determine whether a step from a prior run is valid for reuse.

    Checks:
    1. The stage is in the reusable set (conversion/alignment, Kestrel, adVNTR).
    2. The step was recorded in the prior summary without result_file_missing.
    3. The recorded result_file exists on disk.
    4. The file's current MD5 matches the recorded checksum.
    5. Every required sibling output file exists on disk.

    Args:
        prior: Prior summary dictionary.
        step_name: Step name to evaluate.
        output_root: Root directory of the output.

    Returns:
        bool: True if all artifacts exist and are uncorrupted; False otherwise.
    """
    if step_name not in _REUSABLE_STEPS:
        return False

    prior_steps = prior.get("steps", [])
    step_record = next((s for s in prior_steps if isinstance(s, dict) and s.get("step") == step_name), None)
    if step_record is None:
        logger.debug("Step %r was not recorded in prior summary", step_name)
        return False

    if step_record.get("result_file_missing"):
        logger.debug("Step %r recorded a missing result file in prior summary", step_name)
        return False

    raw_result_file = step_record.get("result_file")
    if not raw_result_file:
        logger.debug("Step %r has no result_file recorded", step_name)
        return False

    recorded_md5 = step_record.get("md5sum")
    if not recorded_md5:
        logger.debug("Step %r has no md5sum recorded", step_name)
        return False

    result_path = Path(raw_result_file)
    if not result_path.is_file():
        # Fallback relative to output_root if directory was relocated
        candidate = Path(output_root) / result_path.name
        if candidate.is_file():
            result_path = candidate
        else:
            logger.debug("Result file %s for step %r does not exist", result_path, step_name)
            return False

    actual_md5 = _compute_md5(result_path)
    if actual_md5 != recorded_md5:
        logger.warning(
            "MD5 mismatch for step %r (%s): recorded %s, current %s",
            step_name,
            result_path,
            recorded_md5,
            actual_md5,
        )
        return False

    # Check required sibling files and verify checksums if recorded
    stage_dir = result_path.parent
    recorded_siblings = prior.get("stage_artifact_md5s", {}).get(step_name, {})
    for sibling in STEP_OUTPUT_SIBLINGS.get(step_name, ()):
        sibling_path = stage_dir / sibling
        if not sibling_path.is_file():
            logger.debug("Required sibling %s for step %r does not exist", sibling_path, step_name)
            return False
        expected_md5 = recorded_siblings.get(sibling)
        if expected_md5 is not None:
            actual_sib_md5 = _compute_md5(sibling_path)
            if actual_sib_md5 != expected_md5:
                logger.warning(
                    "MD5 mismatch for sibling %s of step %r: recorded %s, current %s",
                    sibling,
                    step_name,
                    expected_md5,
                    actual_sib_md5,
                )
                return False

    # Check all recorded sibling checksums and existence
    for sibling, expected_md5 in recorded_siblings.items():
        sibling_path = stage_dir / sibling
        if not sibling_path.is_file():
            logger.debug("Recorded sibling %s for step %r does not exist", sibling_path, step_name)
            return False
        actual_sib_md5 = _compute_md5(sibling_path)
        if actual_sib_md5 != expected_md5:
            logger.warning(
                "MD5 mismatch for recorded sibling %s of step %r: recorded %s, current %s",
                sibling,
                step_name,
                expected_md5,
                actual_sib_md5,
            )
            return False

    # For Kestrel, validate retained BAM replay artifact if present or if result is identity-aware
    if step_name == summary_steps.STEP_KESTREL:
        replay_path = stage_dir / "bam_identity_replay.v1.json"
        if replay_path.is_file():
            try:
                from vntyper.scripts.nomenclature_bam_replay import read_bam_replay_artifact

                read_bam_replay_artifact(stage_dir)
            except (OSError, ValueError) as err:
                logger.warning("BAM replay artifact corrupted or invalid in %s: %s", stage_dir, err)
                return False
        is_identity = False
        parsed_data = step_record.get("parsed_result", {}).get("data", [])
        if parsed_data:
            from vntyper.scripts.nomenclature_dominance_runtime import rows_carry_identity_metadata

            is_identity = rows_carry_identity_metadata(parsed_data)
        if not is_identity and result_path.is_file():
            try:
                with result_path.open("r", encoding="utf-8") as handle:
                    header = handle.readline()
                    if any(col in header for col in ("__Identity_", "ReconciledHaplotype")):
                        is_identity = True
            except (OSError, UnicodeDecodeError):
                pass
        if is_identity and not replay_path.is_file():
            logger.warning("Identity-aware Kestrel result missing required replay artifact %s", replay_path)
            return False

    if step_name == summary_steps.STEP_FASTQ_QC:
        try:
            from vntyper.scripts.generate_report import load_fastp_output

            load_fastp_output(result_path)
        except (ValueError, OSError, json.JSONDecodeError):
            logger.warning("QC report %s is corrupted or invalid", result_path)
            return False

    # Check mate FASTQ for BAM/CRAM conversion if R1 is present
    if "_R1" in result_path.name:
        mate_name = result_path.name.replace("_R1", "_R2")
        mate_path = stage_dir / mate_name
        # If mate was created in the stage, it must still exist
        if not mate_path.is_file():
            logger.debug("Expected mate FASTQ %s for step %r does not exist", mate_path, step_name)
            return False
        expected_mate_md5 = recorded_siblings.get(mate_name)
        if expected_mate_md5 is not None:
            actual_mate_md5 = _compute_md5(mate_path)
            if actual_mate_md5 != expected_mate_md5:
                logger.warning(
                    "MD5 mismatch for mate FASTQ %s of step %r: recorded %s, current %s",
                    mate_name,
                    step_name,
                    expected_mate_md5,
                    actual_mate_md5,
                )
                return False

    return True


def make_reused_step_record(
    prior_step: dict[str, Any],
    prior_pipeline_start: str | None,
) -> dict[str, Any]:
    """Construct a carry-forward step record with reused_from provenance.

    Args:
        prior_step: The exact step record from the prior summary.
        prior_pipeline_start: The pipeline_start timestamp of the donor run.

    Returns:
        dict[str, Any]: Verbatim copy with reused_from set.
    """
    record = dict(prior_step)
    record["reused_from"] = prior_pipeline_start
    return record
