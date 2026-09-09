from __future__ import annotations

import hashlib
import hmac
import json
import logging
import os
import uuid
import zipfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from vntyper.scripts.canonical_json import canonical_sha256

logger = logging.getLogger(__name__)

REPORT_INTEGRITY_VERSION = "2.0"

CANDIDATE_DECISION_FILES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("kestrel", ("kestrel/kestrel_result.tsv", "kestrel_result.tsv")),
    ("coverage", ("coverage/coverage_summary.tsv", "coverage_summary.tsv")),
    ("advntr", ("advntr/output_adVNTR_result.tsv", "output_adVNTR_result.tsv")),
    ("cross_match", ("advntr/cross_match_results.tsv", "cross_match_results.tsv")),
    ("provenance", ("provenance/decision_profile.json", "decision_profile.json")),
)


def _canonical_file_bytes(content: bytes, filename: str) -> bytes:
    """Normalize line endings for text files to ensure deterministic hashing across OSes."""
    if filename.endswith((".tsv", ".csv", ".json", ".txt")):
        text = content.decode("utf-8", errors="replace")
        text = text.replace("\r\n", "\n").replace("\r", "\n").rstrip("\n") + "\n"
        return text.encode("utf-8")
    return content


def _hash_file_content(content: bytes, filename: str) -> str:
    canonical = _canonical_file_bytes(content, filename)
    return hashlib.sha256(canonical).hexdigest()


def _hash_file_content_v1(content: bytes, filename: str) -> str:
    """Legacy v1.0 hashing without trailing newline normalization."""
    hasher = hashlib.sha256()
    if filename.endswith((".tsv", ".csv", ".json", ".txt")):
        content = content.replace(b"\r\n", b"\n")
        lines = content.split(b"\n")
        lines = [line.rstrip() for line in lines]
        content = b"\n".join(lines)
    hasher.update(content)
    return hasher.hexdigest()


def _hash_file(filepath: str | Path) -> str:
    path = Path(filepath)
    content = path.read_bytes()
    return _hash_file_content(content, path.name)


def compute_decision_digest(output_dir_or_files: Path | Mapping[str, Path | str] | str) -> str:
    """Computes deterministic SHA-256 digest over the key decision-bearing files."""
    hasher = hashlib.sha256()

    files_to_hash: list[tuple[str, Path]] = []
    if isinstance(output_dir_or_files, (str, Path)):
        base_dir = Path(output_dir_or_files)
        for key, relative_paths in CANDIDATE_DECISION_FILES:
            for rel_path in relative_paths:
                candidate = base_dir / rel_path
                if candidate.is_file():
                    files_to_hash.append((key, candidate))
                    break
    else:
        for key in sorted(output_dir_or_files.keys()):
            p = Path(output_dir_or_files[key])
            if p.is_file():
                files_to_hash.append((key, p))

    for key, path in files_to_hash:
        file_hash = _hash_file(path)
        hasher.update(f"{key}:{file_hash}\n".encode())

    return hasher.hexdigest()


def _compute_dir_decision_digest_v1(
    output_dir_or_files: Path | Mapping[str, Path | str] | str,
) -> str:
    """Computes legacy v1.0 SHA-256 digest over kestrel and coverage files."""
    hasher = hashlib.sha256()
    if isinstance(output_dir_or_files, (str, Path)):
        base_dir = Path(output_dir_or_files)
        kestrel_file = base_dir / "kestrel_result.tsv"
        if not kestrel_file.exists():
            kestrel_file = base_dir / "kestrel" / "output_indel.vcf.gz"

        coverage_file = base_dir / "coverage" / "coverage_summary.tsv"
        if not coverage_file.exists():
            coverage_file = base_dir / "coverage_summary.tsv"

        files_to_hash: list[tuple[str, Path]] = [("kestrel", kestrel_file), ("coverage", coverage_file)]
    else:
        files_to_hash = []
        for key in sorted(output_dir_or_files.keys()):
            files_to_hash.append((key, Path(output_dir_or_files[key])))

    for key, path in files_to_hash:
        if path.exists():
            file_hash = _hash_file_content_v1(path.read_bytes(), path.name)
            hasher.update(f"{key}:{file_hash}\n".encode())

    return hasher.hexdigest()


def _compute_zip_decision_digest_v1(zip_ref: zipfile.ZipFile, zip_dir: str) -> str:
    """Recomputes legacy v1.0 digest from zip archive."""
    hasher = hashlib.sha256()
    kestrel_path = os.path.join(zip_dir, "kestrel_result.tsv") if zip_dir else "kestrel_result.tsv"
    if kestrel_path not in zip_ref.namelist():
        kestrel_path = (
            os.path.join(zip_dir, "kestrel/output_indel.vcf.gz") if zip_dir else "kestrel/output_indel.vcf.gz"
        )

    coverage_path = (
        os.path.join(zip_dir, "coverage/coverage_summary.tsv") if zip_dir else "coverage/coverage_summary.tsv"
    )
    if coverage_path not in zip_ref.namelist():
        coverage_path = os.path.join(zip_dir, "coverage_summary.tsv") if zip_dir else "coverage_summary.tsv"

    for key, path in [("kestrel", kestrel_path), ("coverage", coverage_path)]:
        if path in zip_ref.namelist():
            with zip_ref.open(path) as f:
                content = f.read()
            file_hash = _hash_file_content_v1(content, os.path.basename(path))
            hasher.update(f"{key}:{file_hash}\n".encode())

    return hasher.hexdigest()


def _compute_zip_decision_digest_v2(zip_ref: zipfile.ZipFile, zip_dir: str) -> str:
    """Recomputes canonical v2.0 digest from zip archive."""
    hasher = hashlib.sha256()
    for key, relative_paths in CANDIDATE_DECISION_FILES:
        for rel_path in relative_paths:
            candidate = os.path.join(zip_dir, rel_path) if zip_dir else rel_path
            if candidate in zip_ref.namelist():
                with zip_ref.open(candidate) as f:
                    content = f.read()
                file_hash = _hash_file_content(content, os.path.basename(candidate))
                hasher.update(f"{key}:{file_hash}\n".encode())
                break
    return hasher.hexdigest()


def _build_integrity_payload_v1(
    run_id: str,
    version: str,
    tool_version: str,
    sample_name: str,
    decision_files_digest: str,
    decision_profile_id: str,
    decision_profile_digest: str,
) -> str:
    return (
        f"{run_id}:{version}:{tool_version}:{sample_name}:"
        f"{decision_files_digest}:{decision_profile_id}:{decision_profile_digest}"
    )


def _build_integrity_payload_v2(
    run_id: str,
    version: str,
    tool_version: str,
    sample_name: str,
    decision_files_digest: str,
    pre_anchor_summary_digest: str,
    decision_profile_id: str,
    decision_profile_sha256: str,
) -> str:
    return (
        f"{run_id}:{version}:{tool_version}:{sample_name}:"
        f"{decision_files_digest}:{pre_anchor_summary_digest}:"
        f"{decision_profile_id}:{decision_profile_sha256}"
    )


_build_integrity_payload = _build_integrity_payload_v2


def anchor_pipeline_summary(
    summary: dict[str, Any],
    output_dir: str | Path,
    secret_key: str | None = None,
) -> dict[str, str]:
    """Calculates and writes integrity anchor fields into the summary dictionary."""
    version = REPORT_INTEGRITY_VERSION
    existing_integrity = summary.get("report_integrity")
    run_id: str = (
        str(existing_integrity.get("run_id"))
        if isinstance(existing_integrity, dict) and existing_integrity.get("run_id")
        else str(uuid.uuid4())
    )

    decision_files_digest = compute_decision_digest(output_dir)

    clean_summary = {k: v for k, v in summary.items() if k != "report_integrity"}
    pre_anchor_summary_digest = canonical_sha256(clean_summary)

    tool_version = str(summary.get("version") or "")
    sample_name = str(summary.get("sample_name") or "")
    decision_profile_id = str(summary.get("decision_profile_id") or "")
    decision_profile_sha256 = str(
        summary.get("decision_profile_sha256") or summary.get("decision_profile_digest") or ""
    )

    payload = _build_integrity_payload(
        run_id,
        version,
        tool_version,
        sample_name,
        decision_files_digest,
        pre_anchor_summary_digest,
        decision_profile_id,
        decision_profile_sha256,
    )

    key = secret_key or os.environ.get("VNTYPER_INTEGRITY_KEY")
    if key:
        mac = hmac.new(key.encode("utf-8"), payload.encode("utf-8"), hashlib.sha256)
        report_integrity_digest = mac.hexdigest()
    else:
        report_integrity_digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()

    summary["report_integrity"] = {
        "report_integrity_version": version,
        "run_id": run_id,
        "decision_files_digest": decision_files_digest,
        "pre_anchor_summary_digest": pre_anchor_summary_digest,
        "report_integrity_digest": report_integrity_digest,
    }

    return summary["report_integrity"]


def verify_report_integrity(archive_path_or_dir: str | Path, secret_key: str | None = None) -> dict[str, Any]:
    """Takes a result ZIP archive or extracted output directory and verifies its cryptographic integrity."""
    is_zip = str(archive_path_or_dir).endswith(".zip")

    summary: dict[str, Any] | None = None
    recomputed_decision_files_digest: str | None = None

    if is_zip:
        with zipfile.ZipFile(archive_path_or_dir, "r") as zip_ref:
            summary_path = next((name for name in zip_ref.namelist() if name.endswith("pipeline_summary.json")), None)
            if not summary_path:
                return {"valid": False, "error": "pipeline_summary.json not found in archive", "signed": False}

            with zip_ref.open(summary_path) as f:
                summary = json.load(f)

            if not summary or not isinstance(summary, dict):
                return {"valid": False, "error": "invalid summary format", "signed": False}

            integrity = summary.get("report_integrity")
            if not integrity or not isinstance(integrity, dict):
                return {"valid": False, "error": "report_integrity missing from summary", "signed": False}

            version = str(integrity.get("report_integrity_version") or "1.0")
            zip_dir = os.path.dirname(summary_path)

            if version == "1.0":
                recomputed_decision_files_digest = _compute_zip_decision_digest_v1(zip_ref, zip_dir)
            elif version == "2.0":
                recomputed_decision_files_digest = _compute_zip_decision_digest_v2(zip_ref, zip_dir)
            else:
                return {
                    "valid": False,
                    "error": f"unsupported report integrity version: {version}",
                    "signed": False,
                }
    else:
        dir_path = Path(archive_path_or_dir)
        summary_path_obj = dir_path / "pipeline_summary.json"
        if not summary_path_obj.exists():
            return {"valid": False, "error": "pipeline_summary.json not found", "signed": False}

        with open(summary_path_obj, encoding="utf-8") as f:
            summary = json.load(f)

        if not summary or not isinstance(summary, dict):
            return {"valid": False, "error": "invalid summary format", "signed": False}

        integrity = summary.get("report_integrity")
        if not integrity or not isinstance(integrity, dict):
            return {"valid": False, "error": "report_integrity missing from summary", "signed": False}

        version = str(integrity.get("report_integrity_version") or "1.0")
        if version == "1.0":
            recomputed_decision_files_digest = _compute_dir_decision_digest_v1(dir_path)
        elif version == "2.0":
            recomputed_decision_files_digest = compute_decision_digest(dir_path)
        else:
            return {
                "valid": False,
                "error": f"unsupported report integrity version: {version}",
                "signed": False,
            }

    run_id = str(integrity.get("run_id") or "")
    decision_files_digest = str(integrity.get("decision_files_digest") or "")
    report_integrity_digest = str(integrity.get("report_integrity_digest") or "")

    tool_version = str(summary.get("version") or "")
    sample_name = str(summary.get("sample_name") or "")
    decision_profile_id = str(summary.get("decision_profile_id") or "")
    decision_profile_sha256 = str(
        summary.get("decision_profile_sha256") or summary.get("decision_profile_digest") or ""
    )

    if recomputed_decision_files_digest != decision_files_digest:
        return {
            "valid": False,
            "run_id": run_id,
            "tool_version": tool_version,
            "decision_files_digest": recomputed_decision_files_digest,
            "error": "decision_files_digest mismatch",
            "signed": False,
        }

    if version == "1.0":
        payload = _build_integrity_payload_v1(
            run_id,
            version,
            tool_version,
            sample_name,
            decision_files_digest,
            decision_profile_id,
            decision_profile_sha256,
        )
    else:
        recorded_pre_anchor_summary_digest = integrity.get("pre_anchor_summary_digest")
        if not recorded_pre_anchor_summary_digest:
            return {
                "valid": False,
                "run_id": run_id,
                "tool_version": tool_version,
                "decision_files_digest": decision_files_digest,
                "error": "pre_anchor_summary_digest missing from report_integrity",
                "signed": False,
            }

        clean_summary = {k: v for k, v in summary.items() if k != "report_integrity"}
        recomputed_pre_anchor_summary_digest = canonical_sha256(clean_summary)

        if recorded_pre_anchor_summary_digest != recomputed_pre_anchor_summary_digest:
            return {
                "valid": False,
                "run_id": run_id,
                "tool_version": tool_version,
                "decision_files_digest": decision_files_digest,
                "error": "summary content tampered (pre_anchor_summary_digest mismatch)",
                "signed": False,
            }

        payload = _build_integrity_payload_v2(
            run_id,
            version,
            tool_version,
            sample_name,
            decision_files_digest,
            str(recorded_pre_anchor_summary_digest),
            decision_profile_id,
            decision_profile_sha256,
        )

    integrity_key = secret_key or os.environ.get("VNTYPER_INTEGRITY_KEY")
    signed = False

    if integrity_key:
        mac = hmac.new(integrity_key.encode("utf-8"), payload.encode("utf-8"), hashlib.sha256)
        expected_digest = mac.hexdigest()
        signed = True
    else:
        hasher = hashlib.sha256(payload.encode("utf-8"))
        expected_digest = hasher.hexdigest()
        if expected_digest != report_integrity_digest:
            return {
                "valid": False,
                "run_id": run_id,
                "tool_version": tool_version,
                "decision_files_digest": decision_files_digest,
                "error": "report_integrity_digest mismatch (possibly signed with a key)",
                "signed": False,
            }

    if expected_digest != report_integrity_digest:
        return {
            "valid": False,
            "run_id": run_id,
            "tool_version": tool_version,
            "decision_files_digest": decision_files_digest,
            "error": "report_integrity_digest mismatch",
            "signed": signed,
        }

    return {
        "valid": True,
        "run_id": run_id,
        "tool_version": tool_version,
        "decision_files_digest": decision_files_digest,
        "error": None,
        "signed": signed,
    }
