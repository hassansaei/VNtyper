import hashlib
import hmac
import json
import logging
import os
import uuid
import zipfile
from collections.abc import Mapping
from pathlib import Path

logger = logging.getLogger(__name__)


def _hash_file_content(content: bytes, filename: str) -> str:
    hasher = hashlib.sha256()

    if filename.endswith((".tsv", ".csv", ".json", ".txt")):
        content = content.replace(b"\r\n", b"\n")
        lines = content.split(b"\n")
        lines = [line.rstrip() for line in lines]
        content = b"\n".join(lines)

    hasher.update(content)
    return hasher.hexdigest()


def _hash_file(filepath: str | Path) -> str:
    with open(filepath, "rb") as f:
        content = f.read()

    filepath_str = filepath.name if isinstance(filepath, Path) else os.path.basename(filepath)
    return _hash_file_content(content, filepath_str)


def compute_decision_digest(output_dir_or_files: Path | Mapping[str, Path] | str) -> str:
    """Computes SHA-256 digest over the key decision-bearing files."""
    hasher = hashlib.sha256()

    if isinstance(output_dir_or_files, (str, Path)):
        base_dir = Path(output_dir_or_files)

        # We need to correctly find the key decision-bearing files based on the problem description
        kestrel_file = base_dir / "kestrel_result.tsv"
        if not kestrel_file.exists():
            kestrel_file = base_dir / "kestrel" / "output_indel.vcf.gz"

        coverage_file = base_dir / "coverage" / "coverage_summary.tsv"
        if not coverage_file.exists():
            coverage_file = base_dir / "coverage_summary.tsv"

        files_to_hash = [("kestrel", kestrel_file), ("coverage", coverage_file)]
    else:
        files_to_hash = []
        for key in sorted(output_dir_or_files.keys()):
            files_to_hash.append((key, output_dir_or_files[key]))

    for key, path in files_to_hash:
        if path and Path(path).exists():
            file_hash = _hash_file(path)
            hasher.update(f"{key}:{file_hash}\n".encode())

    return hasher.hexdigest()


def anchor_pipeline_summary(summary: dict, output_dir: str | Path, secret_key: str | None = None) -> dict:
    """Calculates and writes integrity anchor fields into the summary dictionary."""
    version = "1.0"
    run_id = str(uuid.uuid4())

    decision_files_digest = compute_decision_digest(output_dir)

    tool_version = summary.get("version", "")
    sample_name = summary.get("sample_name", "")
    decision_profile_id = summary.get("decision_profile_id", "")
    decision_profile_digest = summary.get("decision_profile_digest", "")

    payload = f"{run_id}:{version}:{tool_version}:{sample_name}:{decision_files_digest}:{decision_profile_id}:{decision_profile_digest}"

    key = secret_key or os.environ.get("VNTYPER_INTEGRITY_KEY")
    if key:
        mac = hmac.new(key.encode("utf-8"), payload.encode("utf-8"), hashlib.sha256)
        report_integrity_digest = mac.hexdigest()
    else:
        hasher = hashlib.sha256(payload.encode("utf-8"))
        report_integrity_digest = hasher.hexdigest()

    summary["report_integrity"] = {
        "report_integrity_version": version,
        "run_id": run_id,
        "decision_files_digest": decision_files_digest,
        "report_integrity_digest": report_integrity_digest,
    }

    return summary["report_integrity"]


def verify_report_integrity(archive_path_or_dir: str | Path, secret_key: str | None = None) -> dict:
    """
    Takes a result ZIP archive or extracted output directory, reads pipeline_summary.json,
    extracts the decision-bearing files, recomputes the digest, and verifies it.
    """
    is_zip = str(archive_path_or_dir).endswith(".zip")

    summary = None
    recomputed_decision_files_digest = None

    if is_zip:
        with zipfile.ZipFile(archive_path_or_dir, "r") as zip_ref:
            # find pipeline_summary.json
            summary_path = next((name for name in zip_ref.namelist() if name.endswith("pipeline_summary.json")), None)
            if not summary_path:
                return {"valid": False, "error": "pipeline_summary.json not found in archive", "signed": False}

            with zip_ref.open(summary_path) as f:
                summary = json.load(f)

            # Compute decision digest
            hasher = hashlib.sha256()
            base_dir = os.path.dirname(summary_path)

            kestrel_path = os.path.join(base_dir, "kestrel_result.tsv")
            if kestrel_path not in zip_ref.namelist():
                kestrel_path = os.path.join(base_dir, "kestrel/output_indel.vcf.gz")

            coverage_path = os.path.join(base_dir, "coverage/coverage_summary.tsv")
            if coverage_path not in zip_ref.namelist():
                coverage_path = os.path.join(base_dir, "coverage_summary.tsv")

            files_to_hash = [("kestrel", kestrel_path), ("coverage", coverage_path)]

            for key, path in files_to_hash:
                if path in zip_ref.namelist():
                    with zip_ref.open(path) as f:
                        content = f.read()
                    file_hash = _hash_file_content(content, os.path.basename(path))
                    hasher.update(f"{key}:{file_hash}\n".encode())

            recomputed_decision_files_digest = hasher.hexdigest()

    else:
        # It's a directory
        base_dir = Path(archive_path_or_dir)
        summary_path = base_dir / "pipeline_summary.json"
        if not summary_path.exists():
            return {"valid": False, "error": "pipeline_summary.json not found", "signed": False}

        with open(summary_path) as f:
            summary = json.load(f)

        recomputed_decision_files_digest = compute_decision_digest(base_dir)

    integrity = summary.get("report_integrity")
    if not integrity:
        return {"valid": False, "error": "report_integrity missing from summary", "signed": False}

    run_id = integrity.get("run_id")
    version = integrity.get("report_integrity_version")
    decision_files_digest = integrity.get("decision_files_digest")
    report_integrity_digest = integrity.get("report_integrity_digest")

    tool_version = summary.get("version", "")
    sample_name = summary.get("sample_name", "")
    decision_profile_id = summary.get("decision_profile_id", "")
    decision_profile_digest = summary.get("decision_profile_digest", "")

    if recomputed_decision_files_digest != decision_files_digest:
        return {
            "valid": False,
            "run_id": run_id,
            "tool_version": tool_version,
            "decision_files_digest": recomputed_decision_files_digest,
            "error": "decision_files_digest mismatch",
            "signed": False,
        }

    payload = f"{run_id}:{version}:{tool_version}:{sample_name}:{decision_files_digest}:{decision_profile_id}:{decision_profile_digest}"
    key = secret_key or os.environ.get("VNTYPER_INTEGRITY_KEY")
    signed = False

    if key:
        mac = hmac.new(key.encode("utf-8"), payload.encode("utf-8"), hashlib.sha256)
        expected_digest = mac.hexdigest()
        signed = True
    else:
        # Check if it was signed with HMAC
        hasher = hashlib.sha256(payload.encode("utf-8"))
        expected_digest = hasher.hexdigest()
        if expected_digest != report_integrity_digest:
            # We don't have the key, but it might be HMAC signed.
            # But the requirement says "Test HMAC verification with secret key vs without secret key"
            # It also says: `signed: True/False`
            # If no key is provided, we can only verify SHA256. If it fails, maybe it was signed?
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
