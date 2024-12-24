"""
An example of Pytest-based integration tests for VNtyper that load test
scenarios from a JSON config file, making tests more adaptable and extensible.

This version includes a fixture to auto-download files if missing or MD5 mismatch,
and references the new structure: "unit_tests" -> "fastq_tests",
and "integration_tests" -> "bam_tests".
"""

import os
import csv
import json
import math
import shutil
import pytest
import logging
import hashlib
import requests
import subprocess
from pathlib import Path

logging.basicConfig(level=logging.INFO)


@pytest.fixture(scope="session")
def test_config():
    """
    Fixture to load the test scenarios from JSON config.
    This file (test_data_config.json) must exist in tests/ or wherever you prefer.
    """
    config_path = Path("tests/test_data_config.json")
    if not config_path.exists():
        pytest.exit(f"Config file {config_path} not found!", returncode=1)

    with config_path.open("r") as f:
        return json.load(f)


@pytest.fixture(scope="session")
def ensure_test_data(test_config):
    """
    Session-scoped fixture that ensures all test data is present and valid (by MD5).
    Downloads them if missing or if the MD5 does not match.
    """
    logger = logging.getLogger(__name__)
    file_resources = test_config.get("file_resources", [])

    for resource in file_resources:
        local_path = Path(resource["filename"])
        url = resource["url"]
        expected_md5 = resource["md5sum"]

        # Check if file exists and if MD5 matches
        if local_path.exists():
            current_md5 = compute_md5(local_path)
            if current_md5.lower() == expected_md5.lower():
                logger.info("File %s exists and MD5 verified.", local_path)
                continue
            else:
                logger.warning(
                    "File %s exists but MD5 mismatch. Expected=%s, Got=%s.\n"
                    "Re-downloading...",
                    local_path, expected_md5, current_md5
                )
                local_path.unlink()

        # If missing or removed, download now
        logger.info("Downloading file from %s to %s", url, local_path)
        download_file(url, local_path)

        # Verify MD5 after download
        final_md5 = compute_md5(local_path)
        if final_md5.lower() != expected_md5.lower():
            logger.error(
                "Downloaded file %s has incorrect MD5.\n"
                "Expected=%s, Got=%s",
                local_path, expected_md5, final_md5
            )
            pytest.exit(f"MD5 check failed for {local_path}", returncode=1)
        else:
            logger.info("Successfully downloaded and verified %s", local_path)

    return  # This fixture only ensures data is valid.


def compute_md5(file_path: Path) -> str:
    """
    Compute the MD5 hash of a file, returning the hex digest string.
    """
    hasher = hashlib.md5()
    with file_path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def download_file(url: str, dest_path: Path):
    """
    Download a file from 'url' and save it to 'dest_path'.
    Raises an exception on download failure.
    """
    resp = requests.get(url, stream=True, timeout=60)
    resp.raise_for_status()
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(dest_path, "wb") as f:
        for chunk in resp.iter_content(chunk_size=65536):
            f.write(chunk)


#
# 1) FASTQ Tests
#
@pytest.mark.integration
def test_fastq_input(tmp_path, test_config, ensure_test_data, fastq_case):
    """
    Parametrized test for "fastq_tests" items from the JSON config,
    now located under "unit_tests" in the config.

    This test depends on:
      - ensure_test_data => ensures all needed files are downloaded/verified
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting test_fastq_input for case: %s", fastq_case["test_name"])

    # 1) Prepare input and output
    fastq1 = fastq_case["fastq1"]
    fastq2 = fastq_case["fastq2"]
    cli_options = fastq_case["cli_options"]  # e.g. ["--enable-shark"]
    expected_files = fastq_case["expected_files"]
    output_dir = tmp_path / fastq_case["test_name"]

    # Clean up old output
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 2) Build the command
    command = [
        "vntyper",
        "--config-path", "vntyper/config.json",
        "pipeline",
        "--fastq1", fastq1,
        "--fastq2", fastq2,
        "--threads", "4",
        "--output-dir", str(output_dir),
    ] + cli_options

    logger.info("Command to execute: %s", " ".join(command))
    result = subprocess.run(command, capture_output=True, text=True)

    logger.info("Return code: %d", result.returncode)
    logger.info("STDOUT:\n%s", result.stdout)
    logger.info("STDERR:\n%s", result.stderr)

    # 3) Check return code
    assert result.returncode == 0, (
        f"Pipeline returned non-zero exit code.\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    # 4) Check the expected output files
    for rel_path in expected_files:
        check_path = output_dir / rel_path
        logger.info("Checking file: %s", check_path)
        assert check_path.exists(), (
            f"Expected file {rel_path} not found in output directory.\n"
            f"Contents: {list(output_dir.rglob('*'))}"
        )


#
# 2) BAM Tests
#
@pytest.mark.integration
def test_bam_input_with_kestrel_checks(tmp_path, test_config, ensure_test_data, bam_case):
    """
    Parametrized test for all "bam_tests" items from the JSON config,
    which is now under "integration_tests".

    This test depends on:
      - ensure_test_data => ensures all needed files are downloaded
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting test_bam_input_with_kestrel_checks for case: %s", bam_case["test_name"])

    bam_path = bam_case["bam"]
    cli_options = bam_case["cli_options"]
    output_dir = tmp_path / bam_case["test_name"]

    # Clean up old output
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    command = [
        "vntyper",
        "-l", "DEBUG",
        "pipeline",
        "--bam", bam_path,
        "--threads", "4",
        "--reference-assembly", "hg19",
        "-o", str(output_dir)
    ] + cli_options

    logger.info("Command to execute: %s", " ".join(command))
    result = subprocess.run(command, capture_output=True, text=True)

    logger.info("Return code: %d", result.returncode)
    logger.info("STDOUT:\n%s", result.stdout)
    logger.info("STDERR:\n%s", result.stderr)

    # 1) Confirm success
    assert result.returncode == 0, (
        f"Pipeline returned non-zero exit code.\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    # 2) Check the archive if needed
    if bam_case.get("expected_archive"):
        archive_zip = f"{output_dir}.zip"
        logger.info("Looking for %s", archive_zip)
        assert os.path.exists(archive_zip), (
            "No archive created despite --archive-results.\n"
            f"Files in tmp_path: {list(tmp_path.iterdir())}"
        )

    # 3) If we have Kestrel checks, parse kestrel_result.tsv
    if "kestrel_assertions" in bam_case:
        kestrel_tsv = output_dir / "kestrel" / "kestrel_result.tsv"
        logger.info("Looking for %s", kestrel_tsv)
        assert kestrel_tsv.exists(), (
            "kestrel_result.tsv not found in the kestrel subfolder.\n"
            f"Folder contents: {list((output_dir / 'kestrel').iterdir())}"
        )

        with open(kestrel_tsv, "r") as f:
            reader = csv.DictReader((row for row in f if not row.startswith("#")), delimiter="\t")
            rows = list(reader)

        assert len(rows) > 0, "kestrel_result.tsv is empty after skipping comments."

        # Check the first row by default
        row = rows[0]

        if "Estimated_Depth_AlternateVariant" in bam_case["kestrel_assertions"]:
            expected_alt = bam_case["kestrel_assertions"]["Estimated_Depth_AlternateVariant"]
            alt_depth = int(row["Estimated_Depth_AlternateVariant"])
            assert alt_depth == expected_alt, f"Expected alt depth={expected_alt}, got {alt_depth}"

        if "Estimated_Depth_Variant_ActiveRegion" in bam_case["kestrel_assertions"]:
            expected_var = bam_case["kestrel_assertions"]["Estimated_Depth_Variant_ActiveRegion"]
            var_depth = int(row["Estimated_Depth_Variant_ActiveRegion"])
            assert var_depth == expected_var, f"Expected region depth={expected_var}, got {var_depth}"

        if "Depth_Score" in bam_case["kestrel_assertions"]:
            ds_info = bam_case["kestrel_assertions"]["Depth_Score"]
            ds_value = ds_info["value"]
            tolerance_pct = ds_info.get("tolerance_percentage", 5)
            depth_score = float(row["Depth_Score"])

            allowed_variation = ds_value * (tolerance_pct / 100.0)
            diff = abs(depth_score - ds_value)
            assert diff <= allowed_variation, (
                f"Depth_Score mismatch. Got={depth_score}, Expected ~{ds_value} ±{allowed_variation}"
            )

        if "Confidence" in bam_case["kestrel_assertions"]:
            expected_conf = bam_case["kestrel_assertions"]["Confidence"]
            actual_conf = row["Confidence"]
            assert actual_conf == expected_conf, (
                f"Expected Confidence='{expected_conf}', got '{actual_conf}'"
            )

    # 4) If we need igv_report.html
    if bam_case.get("check_igv_report"):
        igv_report = output_dir / "igv_report.html"
        logger.info("Looking for %s", igv_report)
        assert igv_report.exists(), (
            "Expected igv_report.html not found.\n"
            f"Files are: {list(output_dir.iterdir())}"
        )


#
# Parametrization
#

@pytest.fixture(scope="function")
def fastq_case(request, test_config):
    """
    Yields one item from test_config["unit_tests"]["fastq_tests"] at a time.
    """
    return request.param


@pytest.fixture(scope="function")
def bam_case(request, test_config):
    """
    Yields one item from test_config["integration_tests"]["bam_tests"] at a time.
    """
    return request.param


def pytest_generate_tests(metafunc):
    """
    Hook to generate dynamic parameters from 'test_data_config.json'.
    We have new structure: "unit_tests" -> "fastq_tests",
                          "integration_tests" -> "bam_tests".
    """
    # Load config data once
    config_data = metafunc.config._store.get("test_config_data", None)
    if not config_data:
        config_path = Path("tests/test_data_config.json")
        with config_path.open("r") as f:
            config_data = json.load(f)
        metafunc.config._store["test_config_data"] = config_data

    # For FASTQ tests
    if "fastq_case" in metafunc.fixturenames:
        # We get them from "unit_tests" -> "fastq_tests"
        fastq_cases = config_data.get("unit_tests", {}).get("fastq_tests", [])
        metafunc.parametrize(
            "fastq_case", fastq_cases,
            ids=[c["test_name"] for c in fastq_cases]
        )

    # For BAM tests
    if "bam_case" in metafunc.fixturenames:
        # We get them from "integration_tests" -> "bam_tests"
        bam_cases = config_data.get("integration_tests", {}).get("bam_tests", [])
        metafunc.parametrize(
            "bam_case", bam_cases,
            ids=[c["test_name"] for c in bam_cases]
        )