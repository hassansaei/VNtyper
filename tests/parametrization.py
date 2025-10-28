"""
Shared test parametrization for VNtyper.

This module provides utilities to ensure local and Docker tests
execute EXACTLY the same test cases with EXACTLY the same parameters.

Critical for test identity guarantee.
"""

import json
from pathlib import Path
from typing import Any


def load_test_config() -> dict[str, Any]:
    """
    Load test configuration from test_data_config.json.

    Single source of truth for all tests.

    Returns:
        Dict: Test configuration

    Raises:
        FileNotFoundError: If config file not found
    """
    config_path = Path("tests/test_data_config.json")
    if not config_path.exists():
        raise FileNotFoundError(f"Test config not found: {config_path}")

    with config_path.open("r") as f:
        return json.load(f)


def get_bam_test_cases() -> list[dict[str, Any]]:
    """
    Get BAM test cases for parametrization.

    Returns:
        List of test case dicts from integration_tests.bam_tests
    """
    config = load_test_config()
    return config.get("integration_tests", {}).get("bam_tests", [])


def get_bam_test_ids() -> list[str]:
    """
    Get BAM test IDs for parametrization.

    Returns:
        List of test IDs (test_name values)
    """
    return [case["test_name"] for case in get_bam_test_cases()]


def get_advntr_test_cases() -> list[dict[str, Any]]:
    """
    Get adVNTR test cases for parametrization.

    Returns:
        List of test case dicts from integration_tests.advntr_tests
    """
    config = load_test_config()
    return config.get("integration_tests", {}).get("advntr_tests", [])


def get_advntr_test_ids() -> list[str]:
    """
    Get adVNTR test IDs for parametrization.

    Returns:
        List of test IDs
    """
    return [case["test_name"] for case in get_advntr_test_cases()]


def get_fastq_test_cases() -> list[dict[str, Any]]:
    """
    Get FASTQ test cases for parametrization.

    Returns:
        List of test case dicts from integration_tests.fastq_tests
    """
    config = load_test_config()
    return config.get("integration_tests", {}).get("fastq_tests", [])


def get_fastq_test_ids() -> list[str]:
    """
    Get FASTQ test IDs for parametrization.

    Returns:
        List of test IDs
    """
    return [case["test_name"] for case in get_fastq_test_cases()]
