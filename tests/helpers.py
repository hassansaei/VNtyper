"""
Shared test helpers for VNtyper.

These helpers are used by BOTH:
- Pytest integration tests (tests/integration/test_pipeline_integration.py)
- Docker integration tests (tests/docker/test_docker_pipeline.py)

Following principles:
- DRY: Single source of truth for validation logic
- KISS: Simple functions, no classes
- SOLID: Small, focused functions with single responsibility
"""

import csv
import math
from pathlib import Path
from typing import Any, Optional, Union

# ============================================================================
# Parsing Utilities
# ============================================================================


def parse_int_allow_none(value: str) -> Optional[int]:
    """
    Parse string to int, return None if value is None/null/empty.

    Args:
        value: String value from TSV

    Returns:
        int or None

    Examples:
        >>> parse_int_allow_none("123")
        123
        >>> parse_int_allow_none("none")
        None
    """
    if not value or value.lower() in ("none", "null", ""):
        return None
    try:
        return int(value)
    except ValueError:
        return None


def parse_float_allow_none(value: str) -> Optional[float]:
    """
    Parse string to float, return None if value is None/null/empty.

    Args:
        value: String value from TSV

    Returns:
        float or None

    Examples:
        >>> parse_float_allow_none("123.45")
        123.45
        >>> parse_float_allow_none("none")
        None
    """
    if not value or value.lower() in ("none", "null", ""):
        return None
    try:
        return float(value)
    except ValueError:
        return None


# ============================================================================
# File Validation
# ============================================================================


def assert_file_exists(filepath: Path, description: str = "File") -> None:
    """
    Assert file exists and is non-empty.

    Args:
        filepath: Path to check
        description: Description for error message

    Raises:
        AssertionError: If file missing or empty
    """
    assert filepath.exists(), f"{description} not found: {filepath}"
    assert filepath.stat().st_size > 0, f"{description} is empty: {filepath}"


def assert_required_files(output_dir: Path, required_files: list) -> None:
    """
    Assert all required files exist in output directory.

    Args:
        output_dir: Output directory path
        required_files: List of relative file paths to check

    Raises:
        AssertionError: If any file missing or empty
    """
    for file_path in required_files:
        full_path = output_dir / file_path
        assert_file_exists(full_path, f"Required file: {file_path}")


# ============================================================================
# Tolerance-Based Assertions
# ============================================================================


def assert_value_with_tolerance(
    actual: Union[int, float, None],
    expected: Union[int, float, str, dict[str, Any]],
    field_name: str,
    default_tolerance_pct: float = 5.0,
) -> None:
    """
    Assert value matches expected with optional tolerance.

    Supports two formats:
    1. Simple: expected = 123 (exact match)
    2. Tolerance: expected = {"value": 123, "tolerance_percentage": 10}
    3. None: expected = "None" or None

    Args:
        actual: Actual value from output
        expected: Expected value or dict with tolerance
        field_name: Field name for error messages
        default_tolerance_pct: Default tolerance if not specified

    Raises:
        AssertionError: If value outside tolerance

    Examples:
        >>> assert_value_with_tolerance(100, 100, "depth")  # Exact
        >>> assert_value_with_tolerance(95, {"value": 100, "tolerance_percentage": 10}, "depth")
    """
    # Handle expected None
    if expected is None or (isinstance(expected, str) and expected.lower() == "none"):
        assert actual is None, f"{field_name}: Expected None, got {actual}"
        return

    # Extract value and tolerance
    if isinstance(expected, dict):
        expected_val = expected["value"]
        if isinstance(expected_val, str) and expected_val.lower() == "none":
            assert actual is None, f"{field_name}: Expected None, got {actual}"
            return
        tolerance_pct = expected.get("tolerance_percentage", default_tolerance_pct)
    else:
        expected_val = expected
        tolerance_pct = 0  # Exact match for simple values

    # Check actual is not None when expecting a value
    assert actual is not None, f"{field_name}: Expected {expected_val}, got None"

    # Calculate tolerance
    if tolerance_pct == 0:
        # Exact match
        assert actual == expected_val, f"{field_name}: Expected {expected_val}, got {actual}"
    else:
        # Tolerance match
        tolerance = abs(expected_val) * (tolerance_pct / 100.0)
        diff = abs(actual - expected_val)
        assert diff <= tolerance, (
            f"{field_name}: Expected {expected_val} ±{tolerance_pct}% "
            f"(tolerance={tolerance:.2f}), got {actual} (diff={diff:.2f})"
        )


def assert_log10_tolerance(
    actual: float,
    expected: Union[float, dict[str, Any]],
    field_name: str,
    default_log10_tol: float = 2.0,
) -> None:
    """
    Assert P-value using log10 difference tolerance.

    Used for P-values where order of magnitude matters more than exact value.

    Args:
        actual: Actual P-value
        expected: Expected P-value or dict with log10_tolerance
        field_name: Field name for error messages
        default_log10_tol: Default log10 tolerance (default: 2 orders of magnitude)

    Raises:
        AssertionError: If log10 difference exceeds tolerance

    Examples:
        >>> assert_log10_tolerance(1e-7, 6.78e-7, "Pvalue")  # Within 2 orders
    """
    if isinstance(expected, dict):
        expected_val = expected["value"]
        log10_tol = expected.get("log10_tolerance", default_log10_tol)
    else:
        expected_val = expected
        log10_tol = default_log10_tol

    assert actual > 0, f"{field_name}: P-value must be positive, got {actual}"
    assert expected_val > 0, f"{field_name}: Expected P-value must be positive, got {expected_val}"

    log10_diff = abs(math.log10(actual) - math.log10(expected_val))
    assert log10_diff <= log10_tol, (
        f"{field_name}: P-value log10 difference too large. "
        f"Expected {expected_val:.2e}, got {actual:.2e} "
        f"(log10_diff={log10_diff:.2f}, tolerance={log10_tol})"
    )


def assert_pattern_match(actual: str, expected: str, field_name: str) -> None:
    """
    Assert string matches pattern (supports wildcard * suffix).

    Args:
        actual: Actual string value
        expected: Expected pattern (use "prefix*" for wildcard)
        field_name: Field name for error messages

    Raises:
        AssertionError: If pattern doesn't match

    Examples:
        >>> assert_pattern_match("High_Precision*", "High_Precision*", "Confidence")
        >>> assert_pattern_match("High_Precision_v2", "High_Precision*", "Confidence")
    """
    if expected.endswith("*"):
        # Wildcard match
        prefix = expected[:-1]
        assert actual.startswith(prefix), f"{field_name}: Expected to start with '{prefix}', got '{actual}'"
    else:
        # Exact match
        assert actual == expected, f"{field_name}: Expected '{expected}', got '{actual}'"


# ============================================================================
# Kestrel Output Validation
# ============================================================================


def validate_kestrel_output(output_dir: Path, expected: dict[str, Any]) -> None:
    """
    Validate Kestrel genotyping output.

    Args:
        output_dir: Output directory containing kestrel/ subdirectory
        expected: Expected values from test_data_config.json

    Raises:
        AssertionError: If validation fails

    Example expected dict:
        {
            "Estimated_Depth_AlternateVariant": 416,
            "Estimated_Depth_Variant_ActiveRegion": {"value": 7110, "tolerance_percentage": 5},
            "Depth_Score": {"value": 0.0586, "tolerance_percentage": 5},
            "Confidence": "High_Precision*"
        }
    """
    kestrel_file = output_dir / "kestrel" / "kestrel_result.tsv"
    assert_file_exists(kestrel_file, "Kestrel result file")

    # Read TSV
    with kestrel_file.open("r") as f:
        # Skip comment lines starting with #
        lines = [line for line in f if not line.startswith("#")]
        reader = csv.DictReader(lines, delimiter="\t")
        rows = list(reader)

    assert len(rows) > 0, "Kestrel TSV has no data rows"
    row = rows[0]  # First data row

    # Validate each field
    for field, expected_val in expected.items():
        actual_str = row.get(field)
        assert actual_str is not None, f"Kestrel output missing field: {field}"

        if field == "Confidence":
            # String pattern match
            assert_pattern_match(actual_str, expected_val, field)

        elif field in ("Estimated_Depth_AlternateVariant", "Estimated_Depth_Variant_ActiveRegion"):
            # Integer with optional tolerance
            actual_int = parse_int_allow_none(actual_str)
            assert_value_with_tolerance(actual_int, expected_val, field)

        elif field == "Depth_Score":
            # Float with optional tolerance
            actual_float = parse_float_allow_none(actual_str)
            assert_value_with_tolerance(actual_float, expected_val, field)

        else:
            # Generic field
            assert actual_str == str(expected_val), f"{field}: Expected {expected_val}, got {actual_str}"


# ============================================================================
# adVNTR Output Validation
# ============================================================================


def validate_advntr_output(output_dir: Path, expected: dict[str, Any]) -> None:
    """
    Validate adVNTR genotyping output.

    Args:
        output_dir: Output directory containing output_adVNTR_result.tsv
        expected: Expected values from test_data_config.json

    Raises:
        AssertionError: If validation fails

    Example expected dict:
        {
            "VID": "25561",
            "State": "I22_2_G_LEN1",
            "NumberOfSupportingReads": 11,
            "MeanCoverage": {"value": 153.986, "tolerance_percentage": 10},
            "Pvalue": {"value": 6.78e-7, "log10_tolerance": 2}
        }
    """
    advntr_file = output_dir / "advntr" / "output_adVNTR_result.tsv"
    assert_file_exists(advntr_file, "adVNTR result file")

    # Read TSV
    with advntr_file.open("r") as f:
        # Skip comment lines starting with #
        lines = [line for line in f if not line.startswith("#")]
        reader = csv.DictReader(lines, delimiter="\t")
        rows = list(reader)

    assert len(rows) > 0, "adVNTR TSV has no data rows"
    row = rows[0]  # First data row

    # Validate each field
    for field, expected_val in expected.items():
        # Handle "State" → "Variant" column name mapping
        # The TSV header uses "Variant" but test config uses "State" (legacy name)
        actual_field = "Variant" if field == "State" else field
        actual_str = row.get(actual_field)
        assert actual_str is not None, f"adVNTR output missing field: {actual_field} (looking for test field: {field})"

        if field == "VID":
            # String exact match
            assert actual_str == expected_val, f"VID: Expected {expected_val}, got {actual_str}"

        elif field == "State":
            # String exact match (reads from "Variant" column)
            assert actual_str == expected_val, f"State (Variant): Expected {expected_val}, got {actual_str}"

        elif field == "NumberOfSupportingReads":
            # Integer exact match
            actual_int = int(actual_str)
            assert actual_int == expected_val, f"{field}: Expected {expected_val}, got {actual_int}"

        elif field == "MeanCoverage":
            # Float with tolerance
            actual_float = float(actual_str)
            assert_value_with_tolerance(actual_float, expected_val, field, default_tolerance_pct=10)

        elif field == "Pvalue":
            # P-value with log10 tolerance
            actual_pvalue = float(actual_str)
            assert_log10_tolerance(actual_pvalue, expected_val, field)

        else:
            # Generic field
            assert actual_str == str(expected_val), f"{field}: Expected {expected_val}, got {actual_str}"


# ============================================================================
# Coverage Output Validation
# ============================================================================


def validate_coverage_output(output_dir: Path) -> dict[str, float]:
    """
    Validate coverage summary exists and parse basic metrics.

    Args:
        output_dir: Output directory containing coverage/ subdirectory

    Returns:
        Dict with mean_cov, median_cov, uncovered_pct

    Raises:
        AssertionError: If file missing or unparseable
    """
    coverage_file = output_dir / "coverage" / "coverage_summary.tsv"
    assert_file_exists(coverage_file, "Coverage summary file")

    with coverage_file.open("r") as f:
        lines = [line for line in f if line.strip() and not line.startswith("#")]
        reader = csv.DictReader(lines, delimiter="\t")
        rows = list(reader)

    assert len(rows) > 0, "Coverage summary has no data rows"

    # Return first row metrics (basic validation)
    row = rows[0]
    return {
        "mean_cov": float(row.get("Mean", 0)),
        "median_cov": float(row.get("Median", 0)),
        "uncovered_pct": float(row.get("Uncovered%", "0").rstrip("%")),
    }
