"""
Kestrel result validator.

Validates Kestrel genotyping TSV output against expected values from config.
Handles both positive (variant detected) and negative (no variant) cases.
"""

from pathlib import Path
from typing import Any

from ..utils.logging import log
from .base import BaseValidator, ValidationError


class KestrelValidator(BaseValidator):
    """
    Validate Kestrel genotyping results.

    Follows SOLID principles:
    - Single Responsibility: Validate Kestrel output only
    - Open/Closed: Extends BaseValidator without modifying it
    - Liskov Substitution: Can replace BaseValidator anywhere
    """

    def __init__(self, output_dir: Path, expected: dict[str, Any]):
        """
        Initialize Kestrel validator.

        Args:
            output_dir: Pipeline output directory
            expected: Expected Kestrel values from test config
        """
        super().__init__(output_dir)
        self.expected = expected
        self.kestrel_file = output_dir / "kestrel" / "kestrel_result.tsv"

    def validate(self) -> bool:
        """
        Validate Kestrel results.

        Returns:
            True if validation passes

        Raises:
            ValidationError: If validation fails
        """
        log.info("Validating Kestrel genotyping results...")

        # Check file exists
        self.assert_file_exists(self.kestrel_file, "Kestrel results")

        # Parse TSV
        data_line = self._parse_kestrel_tsv()

        # Extract values
        fields = data_line.strip().split("\t")

        # Determine if this is a negative case (10 columns) or positive case (26 columns)
        is_negative_format = len(fields) == 10
        is_positive_format = len(fields) >= 19

        if not (is_negative_format or is_positive_format):
            raise ValidationError(
                f"Kestrel TSV has {len(fields)} columns. "
                f"Expected 10 (negative case) or ≥19 (positive case). "
                "File may be corrupted or incomplete."
            )

        # Extract key metrics (column indices differ between negative and positive formats)
        if is_negative_format:
            # Negative case format: 10 columns
            alt_depth = fields[6]  # Estimated_Depth_AlternateVariant
            active_depth = fields[7]  # Estimated_Depth_Variant_ActiveRegion
            depth_score = fields[8]  # Depth_Score
            confidence = fields[9]  # Confidence
        else:
            # Positive case format: 26 columns
            alt_depth = fields[8]  # Estimated_Depth_AlternateVariant
            active_depth = fields[9]  # Estimated_Depth_Variant_ActiveRegion
            depth_score = fields[17]  # Depth_Score
            confidence = fields[18]  # Confidence

        # Validate based on test type (detect negative case by checking if Confidence == "Negative")
        expected_confidence = self.expected.get("Confidence", "")
        is_negative = expected_confidence == "Negative" or is_negative_format

        if is_negative:
            self._validate_negative_case(alt_depth, active_depth, depth_score, confidence)
        else:
            self._validate_positive_case(alt_depth, active_depth, depth_score, confidence)

        log.success("Kestrel results validation passed")
        return True

    def _parse_kestrel_tsv(self) -> str:
        """
        Parse Kestrel TSV file.

        Returns:
            Data line (non-comment, non-header)

        Raises:
            ValidationError: If parsing fails
        """
        try:
            with self.kestrel_file.open() as f:
                lines = [line for line in f if not line.startswith("##")]

            if len(lines) < 2:
                raise ValidationError("Kestrel TSV has no data rows (only comments/header)")

            # Last non-comment line is the data
            data_line = lines[-1]

            if not data_line.strip():
                raise ValidationError("Kestrel TSV data line is empty")

            return data_line

        except OSError as e:
            raise ValidationError(f"Failed to read Kestrel TSV: {e}") from e

    def _validate_positive_case(
        self, alt_depth: str, active_depth: str, depth_score: str, confidence: str
    ) -> None:
        """
        Validate positive (variant detected) case.

        Args:
            alt_depth: Estimated_Depth_AlternateVariant
            active_depth: Estimated_Depth_Variant_ActiveRegion
            depth_score: Depth_Score
            confidence: Confidence level

        Raises:
            ValidationError: If any validation fails
        """
        # Convert to appropriate types
        try:
            alt_depth_int = int(alt_depth)
            active_depth_int = int(active_depth)
            depth_score_float = float(depth_score)
        except ValueError as e:
            raise ValidationError(f"Failed to parse Kestrel values: {e}") from e

        # Validate Alternate Depth (±15% tolerance)
        expected_alt = self.expected.get("Estimated_Depth_AlternateVariant")
        if expected_alt:
            self.assert_percentage_tolerance(
                actual=alt_depth_int, expected=expected_alt, field_name="Kestrel Alternate Depth", default_tolerance=15.0
            )

        # Validate Active Region Depth (±15% tolerance)
        expected_active = self.expected.get("Estimated_Depth_Variant_ActiveRegion")
        if expected_active:
            self.assert_percentage_tolerance(
                actual=active_depth_int,
                expected=expected_active,
                field_name="Kestrel Active Region Depth",
                default_tolerance=15.0,
            )

        # Validate Depth Score (±15% tolerance)
        expected_score = self.expected.get("Depth_Score")
        if expected_score:
            self.assert_percentage_tolerance(
                actual=depth_score_float, expected=expected_score, field_name="Kestrel Depth Score", default_tolerance=15.0
            )

        # Validate Confidence (pattern matching with wildcard support)
        expected_conf = self.expected.get("Confidence")
        if expected_conf:
            self.assert_pattern_match(actual=confidence, expected_pattern=expected_conf, field_name="Kestrel Confidence")

    def _validate_negative_case(
        self, alt_depth: str, active_depth: str, depth_score: str, confidence: str
    ) -> None:
        """
        Validate negative (no variant detected) case.

        For negative cases, Kestrel outputs various representations:
        - "None", "Negative", "-1", or empty string

        Args:
            alt_depth: Alternate depth value
            active_depth: Active region depth value
            depth_score: Depth score value
            confidence: Confidence value

        Raises:
            ValidationError: If values don't indicate negative case
        """
        log.info("Testing negative case (expecting None/Negative values)")

        # Valid negative indicators
        negative_indicators = ["None", "Negative", "-1", ""]

        # Validate Alternate Depth
        if alt_depth in negative_indicators:
            log.success(f"✓ Kestrel Alternate Depth: {alt_depth or 'empty'} (as expected for negative case)")
        else:
            raise ValidationError(
                f"Negative case: Expected negative indicator for Alternate Depth, " f"got '{alt_depth}'"
            )

        # Validate Active Region Depth
        if active_depth in negative_indicators:
            log.success(f"✓ Kestrel Active Region Depth: {active_depth or 'empty'} (as expected for negative case)")
        else:
            raise ValidationError(
                f"Negative case: Expected negative indicator for Active Region Depth, " f"got '{active_depth}'"
            )

        # Validate Depth Score (may also be "-1.0" for float representation)
        depth_score_indicators = negative_indicators + ["-1.0"]
        if depth_score in depth_score_indicators:
            log.success(f"✓ Kestrel Depth Score: {depth_score or 'empty'} (as expected for negative case)")
        else:
            raise ValidationError(f"Negative case: Expected negative indicator for Depth Score, " f"got '{depth_score}'")

        # Validate Confidence
        if confidence in negative_indicators or confidence == "Negative":
            log.success(f"✓ Kestrel Confidence: {confidence or 'empty'} (as expected for negative case)")
        else:
            raise ValidationError(f"Negative case: Expected 'Negative' or empty for Confidence, " f"got '{confidence}'")
