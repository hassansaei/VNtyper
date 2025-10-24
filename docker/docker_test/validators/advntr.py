"""
adVNTR result validator.

Validates adVNTR module TSV output against expected values from config.
CRITICAL: Fixes data drift issue by loading P-value from test_data_config.json
instead of hardcoded shell script values.
"""

from pathlib import Path
from typing import Any

from ..utils.logging import log
from .base import BaseValidator, ValidationError


class adVNTRValidator(BaseValidator):
    """
    Validate adVNTR genotyping results.

    CRITICAL FIX: Loads expected P-value from config instead of hardcoded value.
    Shell script had: 3.46346905707e-09
    Correct config value: 6.78296229901e-07 (~200x difference!)

    Follows SOLID principles:
    - Single Responsibility: Validate adVNTR output only
    - Open/Closed: Extends BaseValidator without modifying it
    - Liskov Substitution: Can replace BaseValidator anywhere
    """

    def __init__(self, output_dir: Path, expected: dict[str, Any]):
        """
        Initialize adVNTR validator.

        Args:
            output_dir: Pipeline output directory
            expected: Expected adVNTR values from test config
        """
        super().__init__(output_dir)
        self.expected = expected
        self.advntr_file = output_dir / "advntr" / "output_adVNTR_result.tsv"

    def validate(self) -> bool:
        """
        Validate adVNTR results.

        Returns:
            True if validation passes

        Raises:
            ValidationError: If validation fails
        """
        log.info("Validating adVNTR genotyping results...")

        # Check file exists
        self.assert_file_exists(self.advntr_file, "adVNTR results")

        # Parse TSV
        data_line = self._parse_advntr_tsv()

        # Extract values
        fields = data_line.strip().split("\t")

        # Validate we have enough columns
        # Expected: VID, Variant, NumberOfSupportingReads, MeanCoverage, Pvalue, RU, POS, REF, ALT, Flag
        if len(fields) < 5:
            raise ValidationError(
                f"adVNTR TSV has {len(fields)} columns, expected at least 5. " "File may be corrupted or incomplete."
            )

        # Extract key metrics (tab-separated)
        vid = fields[0]
        state = fields[1]  # Variant/State column
        reads = fields[2]  # NumberOfSupportingReads
        mean_cov = fields[3]  # MeanCoverage
        pvalue = fields[4]  # Pvalue

        # Validate each field
        self._validate_vid(vid)
        self._validate_state(state)
        self._validate_reads(reads)
        self._validate_mean_coverage(mean_cov)
        self._validate_pvalue(pvalue)

        log.success("adVNTR results validation passed")
        return True

    def _parse_advntr_tsv(self) -> str:
        """
        Parse adVNTR TSV file.

        Returns:
            Data line (non-comment, non-header)

        Raises:
            ValidationError: If parsing fails
        """
        try:
            with self.advntr_file.open() as f:
                lines = [line for line in f if not line.startswith("#")]

            if not lines:
                raise ValidationError("adVNTR TSV has no data rows (only comments)")

            # First non-comment, non-header line is the data
            # Skip "VID" header if present
            data_lines = [line for line in lines if not line.startswith("VID")]

            if not data_lines:
                raise ValidationError("adVNTR TSV has header but no data rows")

            data_line = data_lines[0]  # First data line

            if not data_line.strip():
                raise ValidationError("adVNTR TSV data line is empty")

            return data_line

        except OSError as e:
            raise ValidationError(f"Failed to read adVNTR TSV: {e}") from e

    def _validate_vid(self, actual: str) -> None:
        """
        Validate VID (exact match).

        Args:
            actual: Actual VID from results

        Raises:
            ValidationError: If VID doesn't match
        """
        expected = self.expected.get("VID")
        if expected is None:
            log.warning("No expected VID in config, skipping validation")
            return

        if actual == str(expected):
            log.success(f"✓ adVNTR VID: {actual}")
        else:
            raise ValidationError(f"adVNTR VID mismatch. Got '{actual}', Expected '{expected}'")

    def _validate_state(self, actual: str) -> None:
        """
        Validate State/Variant (exact match).

        Args:
            actual: Actual state from results

        Raises:
            ValidationError: If state doesn't match
        """
        expected = self.expected.get("State")
        if expected is None:
            log.warning("No expected State in config, skipping validation")
            return

        if actual == expected:
            log.success(f"✓ adVNTR State: {actual}")
        else:
            raise ValidationError(f"adVNTR State mismatch. Got '{actual}', Expected '{expected}'")

    def _validate_reads(self, actual: str) -> None:
        """
        Validate NumberOfSupportingReads (exact match).

        Args:
            actual: Actual reads from results

        Raises:
            ValidationError: If reads don't match
        """
        expected = self.expected.get("NumberOfSupportingReads")
        if expected is None:
            log.warning("No expected NumberOfSupportingReads in config, skipping validation")
            return

        try:
            actual_int = int(actual)
        except ValueError as e:
            raise ValidationError(f"Failed to parse NumberOfSupportingReads '{actual}': {e}")from e

        if actual_int == expected:
            log.success(f"✓ adVNTR Supporting Reads: {actual_int}")
        else:
            raise ValidationError(
                f"adVNTR Supporting Reads mismatch. Got {actual_int}, Expected {expected}"
            )

    def _validate_mean_coverage(self, actual: str) -> None:
        """
        Validate MeanCoverage with percentage tolerance.

        Uses tolerance from config (typically ±10% for adVNTR stochastic variation).

        Args:
            actual: Actual mean coverage from results

        Raises:
            ValidationError: If coverage outside tolerance
        """
        expected = self.expected.get("MeanCoverage")
        if expected is None:
            log.warning("No expected MeanCoverage in config, skipping validation")
            return

        try:
            actual_float = float(actual)
        except ValueError as e:
            raise ValidationError(f"Failed to parse MeanCoverage '{actual}': {e}") from e

        # Use DRY helper from BaseValidator
        self.assert_percentage_tolerance(
            actual=actual_float, expected=expected, field_name="adVNTR Mean Coverage", default_tolerance=10.0
        )

    def _validate_pvalue(self, actual: str) -> None:
        """
        Validate P-value with log10 tolerance.

        CRITICAL FIX: Uses P-value from config, not hardcoded shell script value.

        Shell script hardcoded: 3.46346905707e-09 (WRONG!)
        Config correct value: 6.78296229901e-07
        Difference: ~200x (log10_diff ≈ 2.3)

        Uses log10 tolerance because P-values span many orders of magnitude.
        Tolerance of 2.0 means we accept up to 100x difference (reasonable for
        stochastic adVNTR algorithm).

        Args:
            actual: Actual P-value from results

        Raises:
            ValidationError: If P-value log10 difference exceeds tolerance
        """
        expected = self.expected.get("Pvalue")
        if expected is None:
            log.warning("No expected Pvalue in config, skipping validation")
            return

        try:
            actual_float = float(actual)
        except ValueError as e:
            raise ValidationError(f"Failed to parse Pvalue '{actual}': {e}") from e

        # Use DRY log10 helper from BaseValidator
        # This FIXES the data drift issue!
        self.assert_log10_tolerance(
            actual=actual_float, expected=expected, field_name="adVNTR P-value", default_log10_tol=2.0
        )
