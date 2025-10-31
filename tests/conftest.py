"""
Shared fixtures (like `test_config`) for all tests.
Place this `conftest.py` in your `tests/` directory.

Includes pytest hooks for better progress visibility and timing information.
"""

import json
import logging
import time
from pathlib import Path

import pytest

from tests.test_data_utils import ensure_test_data_downloaded

# ============================================================================
# Logging Configuration
# ============================================================================
logger = logging.getLogger(__name__)


# ============================================================================
# Session-level Fixtures
# ============================================================================
@pytest.fixture(scope="session")
def test_config():
    """
    Single fixture for test configuration used in all tests.

    Loads 'tests/test_data_config.json' by default or
    raises an error if not found.
    """
    config_path = Path("tests/test_data_config.json")
    if not config_path.exists():
        pytest.exit(f"Config file {config_path} not found!", returncode=1)

    with config_path.open("r") as f:
        return json.load(f)


@pytest.fixture(scope="session")
def test_data_dir():
    """
    Session-scoped fixture for test data directory path.

    Computed once per test session and reused across all tests.
    """
    return Path("tests/data")


@pytest.fixture(scope="session")
def reference_dir():
    """
    Session-scoped fixture for reference directory path.

    Computed once per test session and reused across all tests.
    """
    return Path("reference")


@pytest.fixture(scope="session")
def kestrel_jar_path():
    """
    Session-scoped fixture for Kestrel JAR path validation.

    Validates once per session instead of once per test.
    """
    jar_path = Path("vntyper/dependencies/kestrel/kestrel.jar")
    if not jar_path.exists():
        pytest.exit(f"Kestrel JAR not found: {jar_path}", returncode=1)
    return jar_path


@pytest.fixture(scope="session")
def ensure_test_data(test_config):
    """
    Session-scoped fixture that ensures all test data is present and valid.

    Downloads test data from Zenodo if files are missing or MD5 checksums don't match.
    This fixture is used by integration and Docker tests.

    Args:
        test_config: Test configuration from test_data_config.json

    Returns:
        None (side effect: downloads and validates test data)
    """
    ensure_test_data_downloaded(test_config)


# ============================================================================
# Pytest Hooks for Progress Visibility
# ============================================================================


def pytest_collection_modifyitems(config, items):
    """
    Hook called after test collection. Shows how many tests were collected.
    """
    logger.info(f"Collected {len(items)} test(s)")

    # Count integration vs unit tests
    integration_count = sum(1 for item in items if "integration" in [m.name for m in item.iter_markers()])
    unit_count = sum(1 for item in items if "unit" in [m.name for m in item.iter_markers()])

    if integration_count > 0:
        logger.info(f"  - Integration tests: {integration_count}")
    if unit_count > 0:
        logger.info(f"  - Unit tests: {unit_count}")


def pytest_runtest_setup(item):
    """
    Hook called before each test setup phase. Logs test start with timing.
    """
    # Store start time for the test
    item.test_start_time = time.time()

    # Get test markers
    markers = [m.name for m in item.iter_markers()]
    marker_str = f" [{', '.join(markers)}]" if markers else ""

    logger.info(f"{'=' * 80}")
    logger.info(f"STARTING: {item.nodeid}{marker_str}")
    logger.info(f"{'=' * 80}")


def pytest_runtest_call(item):
    """
    Hook called during test execution phase. Logs that test is running.
    """
    logger.info(f"EXECUTING: {item.nodeid}")


def pytest_runtest_teardown(item, nextitem):
    """
    Hook called during test teardown phase. Logs test completion and duration.
    """
    if hasattr(item, "test_start_time"):
        duration = time.time() - item.test_start_time
        logger.info(f"DURATION: {duration:.2f}s for {item.nodeid}")


@pytest.hookimpl(hookwrapper=True, tryfirst=True)
def pytest_runtest_makereport(item, call):
    """
    Hook wrapper to log test outcomes (PASSED/FAILED/SKIPPED) with timing info.
    """
    outcome = yield
    rep = outcome.get_result()

    # Only log at the end of the call phase (actual test execution)
    if rep.when == "call":
        status_emoji = {"passed": "✓", "failed": "✗", "skipped": "⊘"}.get(rep.outcome, "?")

        duration_str = f" ({rep.duration:.2f}s)" if hasattr(rep, "duration") else ""

        if rep.outcome == "passed":
            logger.info(f"{status_emoji} PASSED{duration_str}: {item.nodeid}")
        elif rep.outcome == "failed":
            logger.error(f"{status_emoji} FAILED{duration_str}: {item.nodeid}")
            if hasattr(rep, "longrepr"):
                logger.error(f"  Error: {rep.longrepr}")
        elif rep.outcome == "skipped":
            logger.warning(f"{status_emoji} SKIPPED{duration_str}: {item.nodeid}")


def pytest_sessionstart(session):
    """
    Hook called at the very start of the test session.
    """
    logger.info("=" * 80)
    logger.info("PYTEST SESSION STARTED")
    logger.info("=" * 80)


def pytest_sessionfinish(session, exitstatus):
    """
    Hook called at the end of the test session. Shows final summary.
    """
    logger.info("=" * 80)
    logger.info("PYTEST SESSION FINISHED")
    logger.info(f"Exit status: {exitstatus}")
    logger.info("=" * 80)
