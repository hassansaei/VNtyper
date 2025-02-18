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

import os
import csv
import json
import shutil
import pytest
import logging
import hashlib
import requests
import subprocess
from pathlib import Path

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
    If a file is missing, or if the existing file's MD5 does not match the
    expected value from the config, the file is re-downloaded.

    This fixture only ensures data integrity before tests run.
    No return value is needed.

    Args:
        test_config (dict): The loaded JSON config, which includes a
            "file_resources" list of dicts, each containing:
            - local_path
            - filename
            - url_suffix
            - md5sum
    """
    logger = logging.getLogger(__name__)

    # Obtain the server base URL and the file resources from the config
    server_base_url = test_config.get("server_base_url", "")
    file_resources = test_config.get("file_resources", [])

    # Loop over each resource and verify or download
    for resource in file_resources:
        local_dir = Path(resource["local_path"])
        filename = resource["filename"]
        local_path = local_dir / filename

        url_suffix = resource.get("url_suffix", "")
        expected_md5 = resource["md5sum"]

        logger.info("Handling resource: %s", local_path)

        # If the file exists, check its MD5. If it matches, skip download.
        # Otherwise, remove it to trigger a fresh download.
        if local_path.exists():
            current_md5 = compute_md5(local_path)
            if current_md5.lower() == expected_md5.lower():
                logger.info("File %s exists and MD5 verified.", local_path)
                continue
            else:
                logger.warning(
                    "File %s exists but MD5 mismatch. Expected=%s, Got=%s.\n"
                    "Re-downloading...",
                    local_path,
                    expected_md5,
                    current_md5,
                )
                local_path.unlink()

        # If file is missing (or removed) and we have a valid URL, download it.
        if server_base_url and url_suffix:
            full_url = server_base_url.rstrip("/") + "/" + url_suffix.lstrip("/")
            logger.info("Downloading file from %s to %s", full_url, local_path)
            download_file(full_url, local_path)

            # Verify the MD5 of the newly downloaded file.
            final_md5 = compute_md5(local_path)
            if final_md5.lower() != expected_md5.lower():
                logger.error(
                    "Downloaded file %s has incorrect MD5.\n" "Expected=%s, Got=%s",
                    local_path,
                    expected_md5,
                    final_md5,
                )
                pytest.exit(f"MD5 check failed for {local_path}", returncode=1)
            else:
                logger.info("Successfully downloaded and verified %s", local_path)
        else:
            # If there's no valid download URL and the file doesn't exist,
            # we cannot proceed.
            if not local_path.exists():
                pytest.exit(
                    f"File {local_path} not found and no base_url/url_suffix to download!\n"
                    f"server_base_url={server_base_url}, url_suffix={url_suffix}",
                    returncode=1,
                )

    # This fixture only ensures data integrity; it doesn't return anything.
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
        f"Pipeline returned non-zero exit code.\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    # 4) Check for the expected output files in the output directory
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
def test_bam_input_with_kestrel_checks(
    tmp_path, test_config, ensure_test_data, bam_case
):
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
        "hg19",
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

    # 3) If we have Kestrel checks, parse kestrel_result.tsv and compare
    if "kestrel_assertions" in bam_case:
        kestrel_tsv = output_dir / "kestrel" / "kestrel_result.tsv"
        logger.info("Looking for %s", kestrel_tsv)
        assert kestrel_tsv.exists(), (
            "kestrel_result.tsv not found in the kestrel subfolder.\n"
            f"Folder contents: {list((output_dir / 'kestrel').iterdir())}"
        )

        with open(kestrel_tsv, "r") as f:
            # Skip comment lines (those starting with '#')
            reader = csv.DictReader(
                (row for row in f if not row.startswith("#")), delimiter="\t"
            )
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
            expected_alt = bam_case["kestrel_assertions"][
                "Estimated_Depth_AlternateVariant"
            ]
            alt_str = row["Estimated_Depth_AlternateVariant"]
            alt_val = parse_int_allow_none(alt_str)

            # If expected_alt == "None", we expect alt_val to be None
            if isinstance(expected_alt, str) and expected_alt == "None":
                assert alt_val is None, f"Expected None, got {alt_val}"
            else:
                expected_alt_int = int(expected_alt)
                assert (
                    alt_val == expected_alt_int
                ), f"Expected alt depth={expected_alt_int}, got {alt_val}"

        # Check Estimated_Depth_Variant_ActiveRegion
        if "Estimated_Depth_Variant_ActiveRegion" in bam_case["kestrel_assertions"]:
            expected_var = bam_case["kestrel_assertions"][
                "Estimated_Depth_Variant_ActiveRegion"
            ]
            var_str = row["Estimated_Depth_Variant_ActiveRegion"]
            var_val = parse_int_allow_none(var_str)

            if isinstance(expected_var, str) and expected_var == "None":
                assert var_val is None, f"Expected None, got {var_val}"
            else:
                expected_var_int = int(expected_var)
                assert (
                    var_val == expected_var_int
                ), f"Expected region depth={expected_var_int}, got {var_val}"

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
                assert (
                    row_ds is not None
                ), "Row Depth_Score is None, but config says numeric"
                allowed_variation = abs(config_ds) * (tolerance_pct / 100.0)
                diff = abs(row_ds - config_ds)
                assert diff <= allowed_variation, (
                    f"Depth_Score mismatch. Got={row_ds}, "
                    f"Expected ~{config_ds} ±{allowed_variation}"
                )

        # Check Confidence
        if "Confidence" in bam_case["kestrel_assertions"]:
            expected_conf = bam_case["kestrel_assertions"]["Confidence"]
            actual_conf = row["Confidence"]
            if expected_conf == "Negative":
                logger.info(
                    "Test expects a 'Negative' confidence => skipping strict check."
                )
            else:
                # Some are "High_Precision*" => partial match check
                if expected_conf.endswith("*"):
                    prefix = expected_conf[:-1]  # remove trailing '*'
                    assert actual_conf.startswith(prefix), (
                        f"Expected Confidence to start with '{prefix}', "
                        f"got '{actual_conf}'"
                    )
                else:
                    # Exact match check
                    assert (
                        actual_conf == expected_conf
                    ), f"Expected Confidence='{expected_conf}', got '{actual_conf}'"

    # 4) Check for IGV report if requested
    if bam_case.get("check_igv_report"):
        igv_report = output_dir / "igv_report.html"
        logger.info("Looking for %s", igv_report)
        assert igv_report.exists(), (
            "Expected igv_report.html not found.\n"
            f"Files are: {list(output_dir.iterdir())}"
        )


#
# 3) adVNTR Tests
#
@pytest.mark.integration
def test_advntr_input(tmp_path, test_config, ensure_test_data, advntr_case):
    """
    A new integration test specifically for the adVNTR module.
    We expect a specialized VCF output (e.g. output_adVNTR.vcf)
    and then validate the columns.

    The 'advntr_case' fixture yields a single test scenario with:
      - bam
      - cli_options (including extra modules for adVNTR)
      - expected_vcf
      - advntr_assertions (the expected columns and values)
      - test_name
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting test_advntr_input for case: %s", advntr_case["test_name"])

    bam_path = advntr_case["bam"]
    cli_options = advntr_case["cli_options"]
    expected_vcf = advntr_case["expected_vcf"]
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
        "hg19",
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
        f"Pipeline returned non-zero exit code.\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    # MINIMAL CHANGE: now look for the VCF inside the 'advntr' subfolder
    advntr_dir = output_dir / "advntr"
    vcf_path = advntr_dir / expected_vcf

    logger.info("Looking for %s", vcf_path)
    assert vcf_path.exists(), (
        f"Expected adVNTR VCF {expected_vcf} was not generated in the 'advntr/' folder.\n"
        f"Output folder contents: {list(output_dir.iterdir())}"
    )

    # Parse the adVNTR output, skipping lines starting with '#'
    data_lines = []
    with open(vcf_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue  # skip comment lines
            data_lines.append(line)

    # If the next line is the header (VID State NumberOfSupportingReads ...), skip it too.
    if data_lines and data_lines[0].startswith("VID"):
        logger.info("Skipping adVNTR header line: %s", data_lines[0])
        data_lines.pop(0)

    assert (
        len(data_lines) > 0
    ), f"No data lines found in {vcf_path} after skipping header."

    # Now parse the first data line
    columns = data_lines[0].split("\t")
    # Typical line example:
    # 25561  I22_2_G_LEN1  11  153.986111111  6.78296229901e-07
    assert len(columns) >= 5, "Expected at least 5 columns in adVNTR output."

    advntr_expected = advntr_case["advntr_assertions"]

    # Compare each field as needed.
    # 0 => VID, 1 => State, 2 => NumberOfSupportingReads, 3 => MeanCoverage, 4 => Pvalue
    actual_vid = columns[0]
    assert (
        actual_vid == advntr_expected["VID"]
    ), f"Expected VID={advntr_expected['VID']}, got {actual_vid}"

    actual_state = columns[1]
    assert (
        actual_state == advntr_expected["State"]
    ), f"Expected State={advntr_expected['State']}, got {actual_state}"

    # NumberOfSupportingReads is an integer in the example; parse as float to handle decimals
    actual_num_reads = float(columns[2])
    assert abs(actual_num_reads - advntr_expected["NumberOfSupportingReads"]) < 1e-7, (
        f"Expected NumberOfSupportingReads={advntr_expected['NumberOfSupportingReads']}, "
        f"got {actual_num_reads}"
    )

    actual_mean_cov = float(columns[3])
    assert (
        abs(actual_mean_cov - advntr_expected["MeanCoverage"]) < 1e-7
    ), f"Expected MeanCoverage={advntr_expected['MeanCoverage']}, got {actual_mean_cov}"

    # Compare p-value with a suitable floating tolerance
    actual_pval = float(columns[4])
    assert (
        abs(actual_pval - advntr_expected["Pvalue"]) < 1e-12
    ), f"Expected Pvalue={advntr_expected['Pvalue']}, got {actual_pval}"


#
# Parametrization
#


@pytest.fixture(scope="function")
def fastq_case(request, test_config):
    """
    Yields exactly one item from test_config["unit_tests"]["fastq_tests"].
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

    We handle three categories:
      - "unit_tests" -> "fastq_tests"
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
        fastq_cases = config_data.get("unit_tests", {}).get("fastq_tests", [])
        metafunc.parametrize(
            "fastq_case", fastq_cases, ids=[c["test_name"] for c in fastq_cases]
        )

    # For BAM tests
    if "bam_case" in metafunc.fixturenames:
        bam_cases = config_data.get("integration_tests", {}).get("bam_tests", [])
        metafunc.parametrize(
            "bam_case", bam_cases, ids=[c["test_name"] for c in bam_cases]
        )

    # For adVNTR tests
    if "advntr_case" in metafunc.fixturenames:
        advntr_cases = config_data.get("integration_tests", {}).get("advntr_tests", [])
        metafunc.parametrize(
            "advntr_case", advntr_cases, ids=[c["test_name"] for c in advntr_cases]
        )
