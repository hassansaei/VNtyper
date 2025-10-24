#!/usr/bin/env python3
"""
VNtyper Docker Test Script (Python version).

Replaces test_docker.sh with Python implementation using test_data_config.json.
FIXES data drift by loading expected values from config instead of hardcoding.

Usage:
    python docker/test_docker.py                    # Run all tests
    python docker/test_docker.py example_66bf       # Run specific test
    python docker/test_docker.py example_66bf advntr  # Run specific tests
    python docker/test_docker.py --help             # Show help
"""

import argparse
import sys
from pathlib import Path

# Add docker_test to path
sys.path.insert(0, str(Path(__file__).parent))

from docker_test.config import TestConfig
from docker_test.container import DockerError, DockerOperations
from docker_test.runners.pipeline import PipelineRunner
from docker_test.utils.logging import log


def show_help():
    """Show usage information."""
    help_text = """
VNtyper Docker Test Script

Usage:
  python docker/test_docker.py [OPTIONS] [TEST_NAMES...]

Options:
  -h, --help          Show this help message

Arguments:
  TEST_NAMES          Specific test names to run (optional)
                      If not provided, all tests will be executed

Available Tests:
  Standard BAM tests:
    - example_b178
    - example_a5c1
    - example_66bf
    - example_7a61 (negative case)
    - example_dfc3

  Module tests:
    - advntr          (adVNTR module integration test)

Examples:
  python docker/test_docker.py
    Run all tests (5 standard + 1 adVNTR = 6 total)

  python docker/test_docker.py example_66bf example_7a61
    Run only example_66bf and example_7a61 tests

  python docker/test_docker.py advntr
    Run only the adVNTR module test

  python docker/test_docker.py example_a5c1 advntr
    Run example_a5c1 and adVNTR tests
"""
    print(help_text)


def ensure_test_data(config: TestConfig) -> bool:
    """
    Ensure test data exists.

    Args:
        config: Test configuration

    Returns:
        True if data exists
    """
    test_data_dir = config.test_data_dir

    # Check if at least one test file exists
    test_files = list(test_data_dir.glob("*.bam"))

    if not test_files:
        log.error(f"No test data found in {test_data_dir}")
        log.error("Please run: make download-test-data")
        return False

    log.success(f"Test data available in {test_data_dir}")
    return True


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="VNtyper Docker Test Script", add_help=False)
    parser.add_argument("-h", "--help", action="store_true", help="Show help message")
    parser.add_argument("tests", nargs="*", help="Test names to run")

    args = parser.parse_args()

    if args.help:
        show_help()
        return 0

    log.info("VNtyper Docker Container Test Suite")
    log.info("=" * 60)
    print()

    # Load configuration
    try:
        config = TestConfig()
    except FileNotFoundError as e:
        log.error(f"Configuration error: {e}")
        return 1

    # Ensure test data exists
    if not ensure_test_data(config):
        return 1

    # Initialize Docker operations
    try:
        docker = DockerOperations(image_name="vntyper:latest")
    except DockerError as e:
        log.error(f"Docker error: {e}")
        return 1

    # Build image if needed
    if not docker.image_exists():
        log.info("Docker image not found, building...")
        try:
            docker.build_image(dockerfile=Path("docker/Dockerfile"), context=Path("."), timeout=900)
        except DockerError as e:
            log.error(f"Build failed: {e}")
            return 1
        print()
    else:
        info = docker.get_image_info()
        log.info(f"Using existing image (size: {info['size']})")
        print()

    # Run health checks
    try:
        docker.check_health()
    except DockerError as e:
        log.error(f"Health check failed: {e}")
        return 1
    print()

    # Parse test arguments
    requested_tests = args.tests if args.tests else []
    run_advntr = "advntr" in requested_tests
    if run_advntr:
        requested_tests.remove("advntr")

    run_all = len(requested_tests) == 0 and not run_advntr

    # Show what will be tested
    if run_all:
        log.info("Running all tests (5 standard + 1 adVNTR)")
    else:
        total = len(requested_tests) + (1 if run_advntr else 0)
        log.info(f"Running {total} selected test(s):")
        for test in requested_tests:
            log.info(f"  - {test}")
        if run_advntr:
            log.info("  - advntr")
    print()

    # Initialize runner
    runner = PipelineRunner(config, docker)

    # Run BAM tests
    if run_all:
        results = runner.run_multiple_tests()
    elif requested_tests:
        results = runner.run_multiple_tests(requested_tests)
    else:
        results = {}

    # Run adVNTR test if requested
    advntr_pass = None
    if run_all or run_advntr:
        try:
            advntr_pass = runner.run_advntr_test("example_a5c1_hg19_subset_advntr")
        except Exception as e:
            log.error(f"adVNTR test crashed: {e!r}")
            advntr_pass = False
        print()

    # Final summary
    log.info("=" * 60)
    total_tests = len(results) + (1 if advntr_pass is not None else 0)
    passed_tests = sum(1 for v in results.values() if v) + (1 if advntr_pass else 0)
    failed_tests = total_tests - passed_tests

    if failed_tests == 0:
        log.success("=" * 60)
        log.success("ALL TESTS PASSED ✓")
        log.success("=" * 60)
        log.info(f"Tested: {total_tests} test case(s)")
        return 0
    else:
        log.error("=" * 60)
        log.error("SOME TESTS FAILED ✗")
        log.error("=" * 60)
        log.info(f"Passed: {passed_tests}/{total_tests}")
        log.error(f"Failed: {failed_tests}/{total_tests}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
