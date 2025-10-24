"""
An example of Pytest-based integration tests for VNtyper that load test
scenarios from a JSON config file, making tests more adaptable and extensible.

This version includes:
 - A fixture to auto-download files if missing or MD5 mismatch,
 - Top-level "server_base_url" plus per-file "local_path", "filename", "url_suffix" in the config,
 - Handling "None"/negative values in the Kestrel checks,
 - Acceptance of multi-argument "--extra-modules",
 - A new 'adVNTR' test for advanced VNTR detection.
"""

import csv
import hashlib
import json
import logging
import os
import shutil
import subprocess
from pathlib import Path

import pytest
import requests

# Configure logging for the entire module.
logging.basicConfig(level=logging.INFO)


@pytest.fixture(scope="session")
def test_config():
    """
    Fixture to load the test scenarios from JSON config.
    This file (test_data_config.json) must exist in the `tests/` folder or
    another location if so configured.

    Returns:
        dict: The parsed JSON configuration dictionary. This config contains
              references to files, test scenarios, and more.
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

    This version supports both:
    1. Archive-based downloads (data.zip from Zenodo)
    2. Individual file downloads (legacy support)

    If 'archive_file' is specified in config, it will download and extract the archive
    when any file is missing or has MD5 mismatch. Otherwise, it falls back to individual
    file downloads.

    Args:
        test_config (dict): The loaded JSON config with:
            - archive_file (optional): dict with url, extract_to
            - file_resources: list of file dicts with local_path, filename, md5sum
    """
    logger = logging.getLogger(__name__)
    import tempfile
    import zipfile

    file_resources = test_config.get("file_resources", [])
    archive_config = test_config.get("archive_file")

    # Check if we need to download anything
    need_download = False
    for resource in file_resources:
        local_dir = Path(resource["local_path"])
        filename = resource["filename"]
        local_path = local_dir / filename
        expected_md5 = resource["md5sum"]

        if not local_path.exists():
            logger.info("File %s is missing.", local_path)
            need_download = True
            break

        current_md5 = compute_md5(local_path)
        if current_md5.lower() != expected_md5.lower():
            logger.warning("File %s has MD5 mismatch. Expected=%s, Got=%s", local_path, expected_md5, current_md5)
            need_download = True
            break

    if not need_download:
        logger.info("All test data files verified. No download needed.")
        return

    # Need to download - use archive if configured
    if archive_config:
        archive_url = archive_config["url"]
        extract_to = Path(archive_config["extract_to"])

        logger.info("Downloading test data archive from %s", archive_url)

        with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as tmp_file:
            tmp_path = Path(tmp_file.name)

        try:
            # Download archive
            download_file(archive_url, tmp_path)
            logger.info("Archive downloaded to %s", tmp_path)

            # Extract archive
            logger.info("Extracting archive to %s", extract_to)
            extract_to.mkdir(parents=True, exist_ok=True)

            with zipfile.ZipFile(tmp_path, "r") as zip_ref:
                zip_ref.extractall(extract_to)

            logger.info("Archive extracted successfully")

            # Verify all files after extraction
            for resource in file_resources:
                local_dir = Path(resource["local_path"])
                filename = resource["filename"]
                local_path = local_dir / filename
                expected_md5 = resource["md5sum"]

                if not local_path.exists():
                    pytest.exit(f"File {local_path} not found after archive extraction!", returncode=1)

                current_md5 = compute_md5(local_path)
                if current_md5.lower() != expected_md5.lower():
                    pytest.exit(
                        f"File {local_path} MD5 mismatch after extraction.\nExpected={expected_md5}, Got={current_md5}",
                        returncode=1,
                    )
                logger.info("Verified %s", local_path)
        finally:
            # Clean up temp file
            if tmp_path.exists():
                tmp_path.unlink()
    else:
        # Legacy: download individual files
        logger.info("Using legacy individual file download mode")
        server_base_url = test_config.get("server_base_url", "")

        for resource in file_resources:
            local_dir = Path(resource["local_path"])
            filename = resource["filename"]
            local_path = local_dir / filename
            url_suffix = resource.get("url_suffix", "")
            expected_md5 = resource["md5sum"]

            if local_path.exists():
                current_md5 = compute_md5(local_path)
                if current_md5.lower() == expected_md5.lower():
                    logger.info("File %s exists and MD5 verified.", local_path)
                    continue
                else:
                    logger.warning("File %s MD5 mismatch, re-downloading.", local_path)
                    local_path.unlink()

            if server_base_url and url_suffix:
                full_url = server_base_url.rstrip("/") + "/" + url_suffix.lstrip("/")
                logger.info("Downloading %s from %s", filename, full_url)
                download_file(full_url, local_path)

                final_md5 = compute_md5(local_path)
                if final_md5.lower() != expected_md5.lower():
                    pytest.exit(
                        f"Downloaded file {local_path} MD5 mismatch.\nExpected={expected_md5}, Got={final_md5}",
                        returncode=1,
                    )
            else:
                pytest.exit(f"File {local_path} not found and no download URL configured!", returncode=1)

    return


def compute_md5(file_path: Path) -> str:
    """
    Compute the MD5 hash of a file, returning the hex digest string.

    Args:
        file_path (Path): Path to the file whose MD5 should be computed.

    Returns:
        str: Hex digest string of the file's MD5 sum.
    """
    hasher = hashlib.md5()
    with file_path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def download_file(url: str, dest_path: Path):
    """
    Download a file from 'url' and save it to 'dest_path'.
    Raises an exception if the download fails (non-200 code).

    Args:
        url (str): The remote URL to download from.
        dest_path (Path): The local path where the downloaded file is saved.
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
    Parametrized test for "fastq_tests" items from the JSON config, which
    reside under "unit_tests" -> "fastq_tests".

    This test depends on:
      - ensure_test_data => ensures all needed files are downloaded/verified

    The 'fastq_case' fixture yields a single test scenario with:
      - fastq1 / fastq2
      - cli_options
      - expected_files
      - test_name
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting test_fastq_input for case: %s", fastq_case["test_name"])

    # 1) Prepare input (FASTQ files) and output directory
    fastq1 = fastq_case["fastq1"]
    fastq2 = fastq_case["fastq2"]
    cli_options = fastq_case["cli_options"]
    expected_files = fastq_case["expected_files"]
    output_dir = tmp_path / fastq_case["test_name"]

    # Clean up old output to ensure a fresh start
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 2) Build the CLI command for VNtyper
    command = [
        "vntyper",
        "--config-path",
        "vntyper/config.json",
        "pipeline",
        "--fastq1",
        fastq1,
        "--fastq2",
        fastq2,
        "--threads",
        "4",
        "--output-dir",
        str(output_dir),
    ] + cli_options

    logger.info("Command to execute: %s", " ".join(command))
    # 3) Execute the CLI command
    result = subprocess.run(command, capture_output=True, text=True)

    # Log output for debugging
    logger.info("Return code: %d", result.returncode)
    logger.info("STDOUT:\n%s", result.stdout)
    logger.info("STDERR:\n%s", result.stderr)

    # Check the return code (0 => success)
    assert result.returncode == 0, (
        f"Pipeline returned non-zero exit code.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )

    # 4) Check for the expected output files in the output directory
    for rel_path in expected_files:
        check_path = output_dir / rel_path
        logger.info("Checking file: %s", check_path)
        assert check_path.exists(), (
            f"Expected file {rel_path} not found in output directory.\nContents: {list(output_dir.rglob('*'))}"
        )


#
# 2) BAM Tests
#
@pytest.mark.integration
def test_bam_input_with_kestrel_checks(tmp_path, test_config, ensure_test_data, bam_case):
    """
    Parametrized test for all "bam_tests" items from the JSON config,
    which is now under "integration_tests" -> "bam_tests".

    This test ensures:
      - The pipeline runs successfully with the specified BAM,
      - The (optional) archive is created if requested,
      - Kestrel checks (like depth, confidence, etc.) match expected values,
      - IGV report is present if required.

    The 'bam_case' fixture yields a single test scenario with:
      - bam
      - cli_options
      - expected_archive (optional)
      - kestrel_assertions (optional)
      - check_igv_report (optional)
      - test_name
    """
    logger = logging.getLogger(__name__)
    logger.info(
        "Starting test_bam_input_with_kestrel_checks for case: %s",
        bam_case["test_name"],
    )

    bam_path = bam_case["bam"]
    cli_options = bam_case["cli_options"]
    reference_assembly = bam_case.get("reference_assembly", "hg19")
    output_dir = tmp_path / bam_case["test_name"]

    # Clean up old output
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build the command, including any extra CLI options from the config
    command = [
        "vntyper",
        "-l",
        "DEBUG",
        "pipeline",
        "--bam",
        bam_path,
        "--threads",
        "4",
        "--reference-assembly",
        reference_assembly,
        "-o",
        str(output_dir),
    ] + cli_options

    logger.info("Command to execute: %s", " ".join(command))
    result = subprocess.run(command, capture_output=True, text=True)

    # Log outputs
    logger.info("Return code: %d", result.returncode)
    logger.info("STDOUT:\n%s", result.stdout)
    logger.info("STDERR:\n%s", result.stderr)

    # 1) Confirm success
    assert result.returncode == 0, (
        f"Pipeline returned non-zero exit code.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )

    # 2) Check the archive if needed
    if bam_case.get("expected_archive"):
        archive_zip = f"{output_dir}.zip"
        logger.info("Looking for %s", archive_zip)
        assert os.path.exists(archive_zip), (
            f"No archive created despite --archive-results.\nFiles in tmp_path: {list(tmp_path.iterdir())}"
        )

    # 3) If we have Kestrel checks, parse kestrel_result.tsv and compare
    if "kestrel_assertions" in bam_case:
        kestrel_tsv = output_dir / "kestrel" / "kestrel_result.tsv"
        logger.info("Looking for %s", kestrel_tsv)
        assert kestrel_tsv.exists(), (
            "kestrel_result.tsv not found in the kestrel subfolder.\n"
            f"Folder contents: {list((output_dir / 'kestrel').iterdir())}"
        )

        with open(kestrel_tsv) as f:
            # Skip comment lines (those starting with '#')
            reader = csv.DictReader((row for row in f if not row.startswith("#")), delimiter="\t")
            rows = list(reader)

        assert len(rows) > 0, "kestrel_result.tsv is empty after skipping comments."
        row = rows[0]  # examine the first row by default

        def parse_int_allow_none(val):
            """
            Returns None if val == "None", otherwise returns int(val).
            This handles negative or standard integer strings.
            """
            if val is None or val == "None":
                return None
            return int(val)

        def parse_float_allow_none(val):
            """
            Returns None if val == "None", otherwise returns float(val).
            This handles negative or standard float strings.
            """
            if val is None or val == "None":
                return None
            return float(val)

        # Check Estimated_Depth_AlternateVariant
        if "Estimated_Depth_AlternateVariant" in bam_case["kestrel_assertions"]:
            expected_alt = bam_case["kestrel_assertions"]["Estimated_Depth_AlternateVariant"]
            alt_str = row["Estimated_Depth_AlternateVariant"]
            alt_val = parse_int_allow_none(alt_str)

            # If expected_alt == "None", we expect alt_val to be None
            if isinstance(expected_alt, str) and expected_alt == "None":
                assert alt_val is None, f"Expected None, got {alt_val}"
            elif isinstance(expected_alt, dict):
                # Dict format with tolerance support: {"value": 416, "tolerance_percentage": 5}
                expected_alt_int = int(expected_alt["value"])
                tolerance_pct = expected_alt.get("tolerance_percentage", 5)

                if expected_alt_int is None:
                    assert alt_val is None, f"Expected None, got {alt_val}"
                else:
                    assert alt_val is not None, "Row Estimated_Depth_AlternateVariant is None, but config says numeric"
                    allowed_variation = abs(expected_alt_int) * (tolerance_pct / 100.0)
                    diff = abs(alt_val - expected_alt_int)
                    assert diff <= allowed_variation, (
                        f"Estimated_Depth_AlternateVariant mismatch. Got={alt_val}, "
                        f"Expected ~{expected_alt_int} ±{allowed_variation:.2f}"
                    )
            else:
                # Simple integer format - exact match (backward compatible)
                expected_alt_int = int(expected_alt)
                assert alt_val == expected_alt_int, f"Expected alt depth={expected_alt_int}, got {alt_val}"

        # Check Estimated_Depth_Variant_ActiveRegion
        if "Estimated_Depth_Variant_ActiveRegion" in bam_case["kestrel_assertions"]:
            expected_var = bam_case["kestrel_assertions"]["Estimated_Depth_Variant_ActiveRegion"]
            var_str = row["Estimated_Depth_Variant_ActiveRegion"]
            var_val = parse_int_allow_none(var_str)

            if isinstance(expected_var, str) and expected_var == "None":
                assert var_val is None, f"Expected None, got {var_val}"
            elif isinstance(expected_var, dict):
                # Dict format with tolerance support: {"value": 416, "tolerance_percentage": 5}
                expected_var_int = int(expected_var["value"])
                tolerance_pct = expected_var.get("tolerance_percentage", 5)

                if expected_var_int is None:
                    assert var_val is None, f"Expected None, got {var_val}"
                else:
                    assert var_val is not None, (
                        "Row Estimated_Depth_Variant_ActiveRegion is None, but config says numeric"
                    )
                    allowed_variation = abs(expected_var_int) * (tolerance_pct / 100.0)
                    diff = abs(var_val - expected_var_int)
                    assert diff <= allowed_variation, (
                        f"Estimated_Depth_Variant_ActiveRegion mismatch. Got={var_val}, "
                        f"Expected ~{expected_var_int} ±{allowed_variation:.2f}"
                    )
            else:
                # Simple integer format - exact match (backward compatible)
                expected_var_int = int(expected_var)
                assert var_val == expected_var_int, f"Expected region depth={expected_var_int}, got {var_val}"

        # Check Depth_Score
        if "Depth_Score" in bam_case["kestrel_assertions"]:
            ds_info = bam_case["kestrel_assertions"]["Depth_Score"]
            ds_val_str = ds_info["value"]
            tolerance_pct = ds_info.get("tolerance_percentage", 5)

            row_ds_str = row["Depth_Score"]
            row_ds = parse_float_allow_none(row_ds_str)
            config_ds = parse_float_allow_none(ds_val_str)

            if config_ds is None:
                # Expect "None"
                assert row_ds is None, f"Expected Depth_Score=None, got {row_ds}"
            else:
                # numeric check with tolerance
                assert row_ds is not None, "Row Depth_Score is None, but config says numeric"
                allowed_variation = abs(config_ds) * (tolerance_pct / 100.0)
                diff = abs(row_ds - config_ds)
                assert diff <= allowed_variation, (
                    f"Depth_Score mismatch. Got={row_ds}, Expected ~{config_ds} ±{allowed_variation}"
                )

        # Check Confidence
        if "Confidence" in bam_case["kestrel_assertions"]:
            expected_conf = bam_case["kestrel_assertions"]["Confidence"]
            actual_conf = row["Confidence"]
            if expected_conf == "Negative":
                logger.info("Test expects a 'Negative' confidence => skipping strict check.")
            else:
                # Some are "High_Precision*" => partial match check
                if expected_conf.endswith("*"):
                    prefix = expected_conf[:-1]  # remove trailing '*'
                    assert actual_conf.startswith(prefix), (
                        f"Expected Confidence to start with '{prefix}', got '{actual_conf}'"
                    )
                else:
                    # Exact match check
                    assert actual_conf == expected_conf, f"Expected Confidence='{expected_conf}', got '{actual_conf}'"

    # 4) Check for IGV report if requested
    if bam_case.get("check_igv_report"):
        igv_report = output_dir / "igv_report.html"
        logger.info("Looking for %s", igv_report)
        assert igv_report.exists(), f"Expected igv_report.html not found.\nFiles are: {list(output_dir.iterdir())}"


#
# 3) adVNTR Tests
#
@pytest.mark.integration
def test_advntr_input(tmp_path, test_config, ensure_test_data, advntr_case):
    """
    Integration test for the adVNTR module.

    This test validates the filtered adVNTR results from output_adVNTR_result.tsv,
    which contains the best/most significant variant (lowest p-value) after adVNTR
    filtering. This ensures we test the final processed output rather than raw VCF.

    The 'advntr_case' fixture yields a single test scenario with:
      - bam: path to input BAM file
      - cli_options: CLI flags including extra modules for adVNTR
      - expected_vcf: output filename (e.g., output_adVNTR_result.tsv)
      - advntr_assertions: expected values for VID, State, reads, coverage, p-value
      - test_name: unique identifier for this test case
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting test_advntr_input for case: %s", advntr_case["test_name"])

    bam_path = advntr_case["bam"]
    cli_options = advntr_case["cli_options"]
    expected_vcf = advntr_case["expected_vcf"]
    reference_assembly = advntr_case.get("reference_assembly", "hg19")
    output_dir = tmp_path / advntr_case["test_name"]

    # Clean up old output
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build the command for adVNTR usage
    command = [
        "vntyper",
        "-l",
        "DEBUG",
        "pipeline",
        "--bam",
        bam_path,
        "--threads",
        "4",
        "--reference-assembly",
        reference_assembly,
        "-o",
        str(output_dir),
    ] + cli_options

    logger.info("Command to execute: %s", " ".join(command))
    result = subprocess.run(command, capture_output=True, text=True)

    logger.info("Return code: %d", result.returncode)
    logger.info("STDOUT:\n%s", result.stdout)
    logger.info("STDERR:\n%s", result.stderr)

    # Check that the process completed successfully
    assert result.returncode == 0, (
        f"Pipeline returned non-zero exit code.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )

    # Look for the filtered adVNTR results TSV inside the 'advntr' subfolder
    advntr_dir = output_dir / "advntr"
    result_path = advntr_dir / expected_vcf

    logger.info("Looking for adVNTR result file: %s", result_path)
    assert result_path.exists(), (
        f"Expected adVNTR result file {expected_vcf} was not generated in the 'advntr/' folder.\n"
        f"Output folder contents: {list(output_dir.iterdir())}\n"
        f"adVNTR folder contents: {list(advntr_dir.iterdir()) if advntr_dir.exists() else 'N/A'}"
    )

    # Parse the adVNTR output, skipping comment lines starting with '#'
    data_lines = []
    with open(result_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue  # skip comment lines
            data_lines.append(line)

    # Skip the header line (VID Variant NumberOfSupportingReads MeanCoverage Pvalue ...)
    if data_lines and data_lines[0].startswith("VID"):
        logger.info("Skipping adVNTR header line: %s", data_lines[0])
        data_lines.pop(0)

    assert len(data_lines) > 0, f"No data lines found in {result_path} after skipping header."

    # Parse the first data line (should be the best/filtered result)
    columns = data_lines[0].split("\t")
    logger.info("Parsed adVNTR result columns: %s", columns)
    # Typical line example from output_adVNTR_result.tsv:
    # 25561  I22_2_G_LEN1  13  144.234722222  3.46346905707e-09  2  22  T  TG  Not flagged
    # We need at least the first 5 columns: VID, Variant, NumberOfSupportingReads, MeanCoverage, Pvalue
    assert len(columns) >= 5, "Expected at least 5 columns in adVNTR output."

    advntr_expected = advntr_case["advntr_assertions"]

    # Compare each field as needed.
    # Column indices: 0 => VID, 1 => Variant (called "State" in assertions),
    #                 2 => NumberOfSupportingReads, 3 => MeanCoverage, 4 => Pvalue
    # Note: TSV has additional columns after index 4 (RU, POS, REF, ALT, Flag)
    actual_vid = columns[0]
    assert actual_vid == advntr_expected["VID"], f"Expected VID={advntr_expected['VID']}, got {actual_vid}"

    actual_state = columns[1]
    assert actual_state == advntr_expected["State"], f"Expected State={advntr_expected['State']}, got {actual_state}"

    # NumberOfSupportingReads is an integer in the example; parse as float to handle decimals
    actual_num_reads = float(columns[2])
    assert abs(actual_num_reads - advntr_expected["NumberOfSupportingReads"]) < 1e-7, (
        f"Expected NumberOfSupportingReads={advntr_expected['NumberOfSupportingReads']}, got {actual_num_reads}"
    )

    actual_mean_cov = float(columns[3])
    # Handle MeanCoverage with optional tolerance (similar to Kestrel depth assertions)
    if isinstance(advntr_expected["MeanCoverage"], dict):
        expected_val = advntr_expected["MeanCoverage"]["value"]
        tolerance_pct = advntr_expected["MeanCoverage"].get("tolerance_percentage", 0)
        tolerance = expected_val * (tolerance_pct / 100.0)
        assert abs(actual_mean_cov - expected_val) <= tolerance, (
            f"MeanCoverage mismatch. Got={actual_mean_cov}, Expected ~{expected_val} ±{tolerance} ({tolerance_pct}%)"
        )
    else:
        # Backward compatibility: exact match with tiny float tolerance
        assert abs(actual_mean_cov - advntr_expected["MeanCoverage"]) < 1e-7, (
            f"Expected MeanCoverage={advntr_expected['MeanCoverage']}, got {actual_mean_cov}"
        )

    # Compare p-value with order-of-magnitude tolerance (p-values are stochastic)
    actual_pval = float(columns[4])
    if isinstance(advntr_expected["Pvalue"], dict):
        # New format: structured p-value with log10 tolerance
        expected_val = advntr_expected["Pvalue"]["value"]
        log10_tol = advntr_expected["Pvalue"].get("log10_tolerance", 2)
        import math

        # Compare on log10 scale to allow order-of-magnitude variation
        assert abs(math.log10(actual_pval) - math.log10(expected_val)) <= log10_tol, (
            f"P-value mismatch. Got={actual_pval}, Expected ~{expected_val} "
            f"(must be within {log10_tol} orders of magnitude)"
        )
    else:
        # Backward compatibility: exact match with tiny float tolerance
        assert abs(actual_pval - advntr_expected["Pvalue"]) < 1e-12, (
            f"Expected Pvalue={advntr_expected['Pvalue']}, got {actual_pval}"
        )


#
# Parametrization
#


@pytest.fixture(scope="function")
def fastq_case(request, test_config):
    """
    Yields exactly one item from test_config["integration_tests"]["fastq_tests"].
    Each item is a dict describing a particular FASTQ-based scenario.

    Returns:
        dict: Keys include "test_name", "fastq1", "fastq2", "cli_options",
              "expected_files", etc.
    """
    return request.param


@pytest.fixture(scope="function")
def bam_case(request, test_config):
    """
    Yields exactly one item from test_config["integration_tests"]["bam_tests"].
    Each item is a dict describing a particular BAM-based scenario.

    Returns:
        dict: Keys include "test_name", "bam", "cli_options",
              "expected_archive" (bool), "kestrel_assertions" (dict),
              and "check_igv_report" (bool).
    """
    return request.param


@pytest.fixture(scope="function")
def advntr_case(request, test_config):
    """
    Yields exactly one item from test_config["integration_tests"]["advntr_tests"].
    Each item is a dict describing a particular adVNTR-based scenario.

    Returns:
        dict: Keys include "test_name", "bam", "cli_options",
              "expected_vcf", "advntr_assertions" (dict).
    """
    return request.param


def pytest_generate_tests(metafunc):
    """
    A custom hook that generates dynamic parameters from 'test_data_config.json'.
    This allows each test function (e.g., test_fastq_input, test_bam_input_with_kestrel_checks,
    test_advntr_input) to be parametrized with an array of cases from the JSON file.

    We handle three categories (all under integration_tests):
      - "integration_tests" -> "fastq_tests"
      - "integration_tests" -> "bam_tests"
      - "integration_tests" -> "advntr_tests"
    """
    config_data = metafunc.config._store.get("test_config_data", None)
    if not config_data:
        # Load the JSON config once and store it
        config_path = Path("tests/test_data_config.json")
        with config_path.open("r") as f:
            config_data = json.load(f)
        metafunc.config._store["test_config_data"] = config_data

    # For FASTQ tests
    if "fastq_case" in metafunc.fixturenames:
        fastq_cases = config_data.get("integration_tests", {}).get("fastq_tests", [])
        metafunc.parametrize("fastq_case", fastq_cases, ids=[c["test_name"] for c in fastq_cases])

    # For BAM tests
    if "bam_case" in metafunc.fixturenames:
        bam_cases = config_data.get("integration_tests", {}).get("bam_tests", [])
        metafunc.parametrize("bam_case", bam_cases, ids=[c["test_name"] for c in bam_cases])

    # For adVNTR tests
    if "advntr_case" in metafunc.fixturenames:
        advntr_cases = config_data.get("integration_tests", {}).get("advntr_tests", [])
        metafunc.parametrize("advntr_case", advntr_cases, ids=[c["test_name"] for c in advntr_cases])
