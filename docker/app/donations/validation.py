"""Validation and integrity verification for research data donations."""

from __future__ import annotations

import hashlib
import json
import logging
import tempfile
import zipfile
from pathlib import Path
from typing import Any

from vntyper.scripts.report_integrity import verify_report_integrity

logger = logging.getLogger(__name__)

FORBIDDEN_EXTENSIONS = {
    ".bam",
    ".sam",
    ".cram",
    ".bai",
    ".crai",
    ".fastq",
    ".fq",
    ".fastq.gz",
    ".fq.gz",
    ".fasta",
    ".fa",
    ".fna",
}

ALLOWED_EXTENSIONS = {".json", ".tsv", ".csv", ".txt", ".html", ".log", ".bed", ".vcf", ".vcf.gz"}


def pseudonymize_sample_name(sample_name: str) -> str:
    """Derive an irreversible 16-character hexadecimal pseudonym from a sample name."""
    salt = "vntyper_donation_pseudonym_salt"
    h = hashlib.sha256(f"{salt}:{sample_name.strip()}".encode())
    return f"sample_{h.hexdigest()[:16]}"


def validate_donation_archive(archive_path_or_bytes: Path | bytes) -> dict[str, Any]:
    """Inspect and verify a stripped VNtyper result archive.

    Enforces:
    1. No raw sequences (BAM, CRAM, FASTQ) remain in the archive.
    2. Decision-bearing files are present.
    3. Report integrity cryptographic anchor validates against decision files.

    Returns:
        dict with extracted and validated data:
        - run_id
        - tool_version
        - sample_hash
        - decision_files_digest
        - report_integrity_digest
        - mean_coverage
        - flank_mean_depth
        - kestrel_result
        - is_positive_pipeline
    """
    if isinstance(archive_path_or_bytes, bytes):
        with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as tmp:
            tmp.write(archive_path_or_bytes)
            tmp_path = Path(tmp.name)
    else:
        tmp_path = Path(archive_path_or_bytes)

    try:
        if not zipfile.is_zipfile(tmp_path):
            raise ValueError("Uploaded file is not a valid ZIP archive.")

        with zipfile.ZipFile(tmp_path, "r") as zf:
            namelist = zf.namelist()

            # 1. Check for forbidden sequence files
            for name in namelist:
                lower = name.lower()
                for ext in FORBIDDEN_EXTENSIONS:
                    if lower.endswith(ext):
                        raise ValueError(
                            f"Archive contains raw sequence file '{name}' with forbidden extension '{ext}'. "
                            "You must strip all BAM/CRAM/FASTQ files prior to donating data."
                        )

            # 2. Check for required summary file
            summary_member = None
            for name in namelist:
                if name.endswith("pipeline_summary.json"):
                    summary_member = name
                    break

            if not summary_member:
                raise ValueError("Missing 'pipeline_summary.json' in result archive.")

            summary_bytes = zf.read(summary_member)
            try:
                summary = json.loads(summary_bytes.decode("utf-8"))
            except (json.JSONDecodeError, UnicodeDecodeError, ValueError) as e:
                raise ValueError(f"Corrupted pipeline_summary.json: {e}") from e

        # 3. Verify report integrity anchor
        verification = verify_report_integrity(tmp_path)
        if not verification.get("valid"):
            err = verification.get("error") or "Cryptographic checksum mismatch"
            raise ValueError(f"Report integrity verification failed: {err}. Archive may have been altered.")

        run_id = verification["run_id"]
        tool_version = verification.get("tool_version") or summary.get("version", "unknown")
        decision_digest = verification["decision_files_digest"]
        integrity_digest = summary.get("report_integrity", {}).get("report_integrity_digest", "")

        sample_name = summary.get("sample_name", "unknown")
        sample_hash = pseudonymize_sample_name(sample_name)

        # Extract coverage metrics
        coverage_dict = summary.get("coverage", {})
        mean_cov = coverage_dict.get("mean")
        flank_mean = coverage_dict.get("vntr_flank_mean_depth")

        # Extract algorithm calls
        screening = summary.get("screening", {})
        is_positive = screening.get("is_positive", False)
        kestrel_res = screening.get("kestrel_result", "unknown")

        return {
            "run_id": run_id,
            "tool_version": tool_version,
            "sample_hash": sample_hash,
            "decision_files_digest": decision_digest,
            "report_integrity_digest": integrity_digest,
            "mean_coverage": float(mean_cov) if mean_cov is not None else None,
            "flank_mean_depth": float(flank_mean) if flank_mean is not None else None,
            "kestrel_result": str(kestrel_res),
            "is_positive_pipeline": bool(is_positive),
        }

    finally:
        if isinstance(archive_path_or_bytes, bytes) and tmp_path.exists():
            tmp_path.unlink()
