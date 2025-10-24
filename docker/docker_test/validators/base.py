"""
Base validator class with tolerance helpers.

Implements DRY principle for tolerance-based assertions used across all validators.
"""

import math
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Union

from ..utils.logging import log


class ValidationError(Exception):
    """Raised when validation fails."""

    pass


class BaseValidator(ABC):
    """
    Abstract base validator with shared tolerance methods.

    Follows SOLID principles:
    - Single Responsibility: Validation logic only
    - Open/Closed: Extend for specific validators, closed for modification
    - Liskov Substitution: All validators can replace this base
    - Interface Segregation: Small, focused interface
    - Dependency Inversion: Depends on abstractions (pathlib.Path, not concrete paths)
    """

    def __init__(self, output_dir: Path):
        """
        Initialize validator.

        Args:
            output_dir: Pipeline output directory to validate
        """
        self.output_dir = output_dir

    @abstractmethod
    def validate(self) -> bool:
        """
        Validate results.

        Returns:
            True if validation passes

        Raises:
            ValidationError: If validation fails
        """
        pass

    @staticmethod
    def assert_percentage_tolerance(
        actual: float, expected: Union[float, dict[str, Any]], field_name: str, default_tolerance: float = 15.0
    ) -> None:
        """
        Assert value is within percentage tolerance.

        Handles both simple float expected values and dict with {value, tolerance_percentage}.

        Args:
            actual: Actual value from results
            expected: Expected value (float) or dict with {"value": float, "tolerance_percentage": float}
            field_name: Name of field being validated (for error messages)
            default_tolerance: Default tolerance percentage if not specified (default: 15%)

        Raises:
            ValidationError: If value outside tolerance or expected value invalid
        """
        if isinstance(expected, dict):
            expected_val = expected["value"]
            tolerance_pct = expected.get("tolerance_percentage", default_tolerance)
        else:
            expected_val = expected
            tolerance_pct = default_tolerance

        # Edge case: zero expected value
        if expected_val == 0:
            if actual != 0:
                raise ValidationError(f"{field_name}: Got {actual}, Expected 0 (exact match required for zero)")
            log.success(f"✓ {field_name}: {actual} (exact match for zero)")
            return

        # Edge case: negative expected value (invalid)
        if expected_val < 0:
            raise ValidationError(f"{field_name}: Negative expected value {expected_val} is invalid")

        # Calculate tolerance bounds
        tolerance = abs(expected_val) * (tolerance_pct / 100.0)
        min_val = expected_val - tolerance
        max_val = expected_val + tolerance

        # Validate
        if min_val <= actual <= max_val:
            log.success(f"✓ {field_name}: {actual} (expected: {expected_val} ±{tolerance_pct}%)")
        else:
            raise ValidationError(
                f"{field_name}: Got {actual}, Expected {expected_val} ±{tolerance_pct}% "
                f"(valid range: {min_val:.2f} - {max_val:.2f})"
            )

    @staticmethod
    def assert_log10_tolerance(
        actual: float, expected: Union[float, dict[str, Any]], field_name: str, default_log10_tol: float = 2.0
    ) -> None:
        """
        Assert log10 difference is within tolerance.

        Used for P-values where relative differences matter more than absolute.

        Args:
            actual: Actual P-value
            expected: Expected P-value (float) or dict with {"value": float, "log10_tolerance": float}
            field_name: Name of field (for error messages)
            default_log10_tol: Default log10 tolerance if not specified (default: 2.0)

        Raises:
            ValidationError: If log10 difference exceeds tolerance or p-value invalid
        """
        if isinstance(expected, dict):
            expected_val = expected["value"]
            log10_tol = expected.get("log10_tolerance", default_log10_tol)
        else:
            expected_val = expected
            log10_tol = default_log10_tol

        # Edge case: zero or negative p-values (invalid)
        if actual <= 0:
            raise ValidationError(f"{field_name}: P-value must be positive. Got {actual}")

        if expected_val <= 0:
            raise ValidationError(f"{field_name}: Expected P-value must be positive. Got {expected_val}")

        # Calculate log10 difference
        try:
            log10_diff = abs(math.log10(actual) - math.log10(expected_val))
        except ValueError as e:
            raise ValidationError(f"{field_name}: Cannot compute log10. actual={actual}, expected={expected_val}: {e}") from e

        # Validate
        if log10_diff <= log10_tol:
            log.success(f"✓ {field_name}: {actual:.2e} (expected: {expected_val:.2e}, log10_diff: {log10_diff:.2f})")
        else:
            raise ValidationError(
                f"{field_name}: P-value log10 difference too large. "
                f"Got {actual:.2e}, Expected {expected_val:.2e} "
                f"(log10_diff: {log10_diff:.2f} > {log10_tol})"
            )

    @staticmethod
    def assert_exact_or_tolerance(
        actual: Any,
        expected: Union[Any, dict[str, Any]],
        field_name: str,
        tolerance_mode: str = "exact",
        default_tolerance: float = 1e-7,
    ) -> None:
        """
        Assert exact match or tolerance-based match.

        Handles backward compatibility with simple values vs dict with tolerance specs.

        Args:
            actual: Actual value
            expected: Expected value or dict with tolerance spec
            field_name: Field name for error messages
            tolerance_mode: "exact", "percentage", or "log10"
            default_tolerance: Default tolerance if not specified

        Raises:
            ValidationError: If validation fails
        """
        if isinstance(expected, dict):
            # Use appropriate tolerance method based on mode
            if tolerance_mode == "percentage":
                BaseValidator.assert_percentage_tolerance(actual, expected, field_name)
            elif tolerance_mode == "log10":
                BaseValidator.assert_log10_tolerance(actual, expected, field_name)
            else:
                # Exact match from dict
                expected_val = expected.get("value", expected)
                if abs(actual - expected_val) <= default_tolerance:
                    log.success(f"✓ {field_name}: {actual}")
                else:
                    raise ValidationError(f"{field_name}: Got {actual}, Expected {expected_val}")
        else:
            # Simple value comparison
            if tolerance_mode == "exact":
                if isinstance(actual, float) and isinstance(expected, float):
                    if abs(actual - expected) <= default_tolerance:
                        log.success(f"✓ {field_name}: {actual}")
                    else:
                        raise ValidationError(f"{field_name}: Got {actual}, Expected {expected}")
                else:
                    if actual == expected:
                        log.success(f"✓ {field_name}: {actual}")
                    else:
                        raise ValidationError(f"{field_name}: Got {actual}, Expected {expected}")
            else:
                raise ValueError(f"Tolerance mode '{tolerance_mode}' requires expected to be a dict")

    @staticmethod
    def assert_file_exists(filepath: Path, description: str = "File") -> None:
        """
        Assert file exists and is non-empty.

        Args:
            filepath: Path to file
            description: Description for error message

        Raises:
            ValidationError: If file missing or empty
        """
        if not filepath.exists():
            raise ValidationError(f"{description} not found: {filepath}")

        if filepath.stat().st_size == 0:
            raise ValidationError(f"{description} is empty: {filepath}")

        log.success(f"✓ Found: {filepath.name}")

    @staticmethod
    def assert_in_range(
        actual: Union[int, float], min_val: Union[int, float], max_val: Union[int, float], field_name: str
    ) -> None:
        """
        Assert value is within range [min_val, max_val].

        Args:
            actual: Actual value
            min_val: Minimum allowed value
            max_val: Maximum allowed value
            field_name: Field name for error message

        Raises:
            ValidationError: If value outside range
        """
        if min_val <= actual <= max_val:
            log.success(f"✓ {field_name}: {actual} (range: {min_val} - {max_val})")
        else:
            raise ValidationError(f"{field_name}: Got {actual}, Expected range [{min_val}, {max_val}]")

    @staticmethod
    def assert_pattern_match(actual: str, expected_pattern: str, field_name: str) -> None:
        """
        Assert string matches pattern (supports wildcard * suffix).

        Args:
            actual: Actual string value
            expected_pattern: Expected pattern (supports "prefix*" wildcard)
            field_name: Field name for error message

        Raises:
            ValidationError: If pattern doesn't match
        """
        if expected_pattern.endswith("*"):
            # Wildcard match
            prefix = expected_pattern[:-1]
            if actual.startswith(prefix):
                log.success(f"✓ {field_name}: {actual} (matches pattern: {expected_pattern})")
            else:
                raise ValidationError(f"{field_name}: Got '{actual}', Expected pattern '{expected_pattern}'")
        else:
            # Exact match
            if actual == expected_pattern:
                log.success(f"✓ {field_name}: {actual}")
            else:
                raise ValidationError(f"{field_name}: Got '{actual}', Expected '{expected_pattern}'")
