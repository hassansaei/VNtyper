"""
Configuration loader for Docker tests.

Loads test configuration from test_data_config.json to ensure single source
of truth for expected values. Fixes data drift issue between shell script and pytest.
"""

import json
from pathlib import Path
from typing import Any, Optional

from .utils.logging import log


class TestConfig:
    """Configuration loader for Docker tests."""

    def __init__(self, config_path: Optional[Path] = None):
        """
        Initialize configuration loader.

        Args:
            config_path: Path to test_data_config.json. Defaults to tests/test_data_config.json
        """
        if config_path is None:
            config_path = Path("tests/test_data_config.json")

        self.config_path = config_path
        self._config: Optional[dict[str, Any]] = None
        self._load_config()

    def _load_config(self) -> None:
        """Load configuration from JSON file."""
        if not self.config_path.exists():
            raise FileNotFoundError(f"Config file not found: {self.config_path}")

        with self.config_path.open("r") as f:
            self._config = json.load(f)

        log.success(f"Loaded test configuration from {self.config_path}")

    def get_bam_tests(self) -> list[dict[str, Any]]:
        """
        Get BAM test cases configuration.

        Returns:
            List of test configurations
        """
        if self._config is None:
            raise RuntimeError("Configuration not loaded")

        integration_tests = self._config.get("integration_tests", {})
        return integration_tests.get("bam_tests", [])

    def get_advntr_tests(self) -> list[dict[str, Any]]:
        """
        Get adVNTR test cases configuration.

        Returns:
            List of adVNTR test configurations
        """
        if self._config is None:
            raise RuntimeError("Configuration not loaded")

        integration_tests = self._config.get("integration_tests", {})
        return integration_tests.get("advntr_tests", [])

    def get_bam_test(self, test_name: str) -> dict[str, Any]:
        """
        Get specific BAM test configuration.

        Args:
            test_name: Test case name (e.g., "example_b178_hg19_subset_fast")

        Returns:
            Test configuration dictionary

        Raises:
            KeyError: If test_name not found
        """
        bam_tests = self.get_bam_tests()
        for test in bam_tests:
            if test.get("test_name") == test_name:
                return test

        available = [t.get("test_name") for t in bam_tests]
        raise KeyError(f"Test name not found: {test_name}. Available: {available}")

    def get_advntr_test(self, test_name: str) -> dict[str, Any]:
        """
        Get specific adVNTR test configuration.

        Args:
            test_name: Test case name (e.g., "example_a5c1_hg19_subset_advntr")

        Returns:
            adVNTR test configuration dictionary

        Raises:
            KeyError: If test_name not found
        """
        advntr_tests = self.get_advntr_tests()
        for test in advntr_tests:
            if test.get("test_name") == test_name:
                return test

        available = [t.get("test_name") for t in advntr_tests]
        raise KeyError(f"adVNTR test name not found: {test_name}. Available: {available}")

    def list_bam_test_ids(self) -> list[str]:
        """Get list of all BAM test names."""
        return sorted([t["test_name"] for t in self.get_bam_tests() if "test_name" in t])

    def list_advntr_test_ids(self) -> list[str]:
        """Get list of all adVNTR test names."""
        return sorted([t["test_name"] for t in self.get_advntr_tests() if "test_name" in t])

    @property
    def zenodo_record(self) -> str:
        """Get Zenodo record ID for test data."""
        if self._config is None:
            raise RuntimeError("Configuration not loaded")

        return self._config.get("zenodo_record", "17377082")

    @property
    def test_data_dir(self) -> Path:
        """Get test data directory path."""
        return Path("tests/data")

    @property
    def docker_test_dir(self) -> Path:
        """Get Docker test output directory path."""
        return Path("docker/test_output")
