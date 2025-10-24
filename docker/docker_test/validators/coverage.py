"""
Coverage analysis validator.

Validates coverage summary TSV output with sanity checks.
Uses warnings instead of hard failures since coverage varies by sample quality.
"""

import contextlib
from pathlib import Path

from ..utils.logging import log
from .base import BaseValidator, ValidationError


class CoverageValidator(BaseValidator):
    """
    Validate coverage analysis results.

    Note: Uses warnings instead of hard failures because coverage metrics
    vary based on sequencing depth, read quality, and VNTR complexity.
    Kestrel and adVNTR provide the definitive pass/fail validation.

    Follows SOLID principles:
    - Single Responsibility: Validate coverage output only
    - Open/Closed: Extends BaseValidator without modifying it
    """

    # Coverage thresholds for warnings
    MIN_MEAN_COVERAGE = 1000  # Minimum mean coverage for "sufficient"
    MIN_MEDIAN_COVERAGE = 500  # Minimum median coverage for "sufficient"
    MAX_UNCOVERED_PERCENT = 5.0  # Maximum uncovered percentage for "excellent"

    def __init__(self, output_dir: Path):
        """
        Initialize coverage validator.

        Args:
            output_dir: Pipeline output directory
        """
        super().__init__(output_dir)
        self.coverage_file = output_dir / "coverage" / "coverage_summary.tsv"

    def validate(self) -> bool:
        """
        Validate coverage results.

        Returns:
            True if file exists and parseable (always passes, uses warnings)

        Raises:
            ValidationError: Only if file missing or unparseable
        """
        log.info("Validating coverage analysis results...")

        # Check file exists
        self.assert_file_exists(self.coverage_file, "Coverage summary")

        # Parse TSV
        coverage_data = self._parse_coverage_tsv()

        # Extract metrics
        mean_cov = coverage_data.get("mean_cov")
        median_cov = coverage_data.get("median_cov")
        uncovered_pct = coverage_data.get("uncovered_pct")

        # Validate each metric (with warnings, not failures)
        if mean_cov is not None:
            self._check_mean_coverage(mean_cov)

        if median_cov is not None:
            self._check_median_coverage(median_cov)

        if uncovered_pct is not None:
            self._check_uncovered_percentage(uncovered_pct)

        log.success("Coverage validation completed")
        return True

    def _parse_coverage_tsv(self) -> dict:
        """
        Parse coverage summary TSV file.

        Expected format (tab-separated):
        MeanCoverage  MedianCoverage  Q1  MinCoverage  MaxCoverage  Q3  UncoveredBases

        Returns:
            Dictionary with coverage metrics

        Raises:
            ValidationError: If parsing fails
        """
        try:
            with self.coverage_file.open() as f:
                lines = [line for line in f if line.strip() and not line.startswith("#")]

            if len(lines) < 2:
                raise ValidationError("Coverage summary has no data rows (header only)")

            # Last line is the data (skip header)
            data_line = lines[-1].strip()
            fields = data_line.split("\t")

            if len(fields) < 7:
                raise ValidationError(
                    f"Coverage TSV has {len(fields)} columns, expected at least 7. " "File may be corrupted."
                )

            # Parse values
            try:
                mean_cov = float(fields[0])
                median_cov = float(fields[1])
                # fields[2] = Q1
                # fields[3] = MinCoverage
                # fields[4] = MaxCoverage
                # fields[5] = Q3

                # UncoveredBases can be "count\tpercentage" or just count
                uncovered_field = fields[6] if len(fields) > 6 else ""

                # Try to extract percentage
                uncovered_pct = None
                if len(fields) > 7:
                    # Percentage is in field 7
                    uncovered_pct_str = fields[7].strip().rstrip("%")
                    try:
                        uncovered_pct = float(uncovered_pct_str)
                    except ValueError:
                        log.warning(f"Could not parse uncovered percentage: '{fields[7]}'")
                else:
                    # Try to extract from field 6 if it contains "%"
                    if "%" in uncovered_field:
                        uncovered_pct_str = uncovered_field.split()[-1].strip().rstrip("%")
                        with contextlib.suppress(ValueError):
                            uncovered_pct = float(uncovered_pct_str)

                return {"mean_cov": mean_cov, "median_cov": median_cov, "uncovered_pct": uncovered_pct}

            except (ValueError, IndexError) as e:
                raise ValidationError(f"Failed to parse coverage values: {e}") from e

        except OSError as e:
            raise ValidationError(f"Failed to read coverage summary: {e}") from e

    def _check_mean_coverage(self, mean_cov: float) -> None:
        """
        Check mean coverage against threshold.

        Args:
            mean_cov: Mean VNTR coverage
        """
        if mean_cov >= self.MIN_MEAN_COVERAGE:
            log.success(f"✓ Mean VNTR coverage: {mean_cov:.1f}× (sufficient)")
        else:
            log.warning(f"⚠ Mean VNTR coverage: {mean_cov:.1f}× (below {self.MIN_MEAN_COVERAGE}×, may be low)")

    def _check_median_coverage(self, median_cov: float) -> None:
        """
        Check median coverage against threshold.

        Args:
            median_cov: Median VNTR coverage
        """
        if median_cov >= self.MIN_MEDIAN_COVERAGE:
            log.success(f"✓ Median VNTR coverage: {median_cov:.1f}× (sufficient)")
        else:
            log.warning(
                f"⚠ Median VNTR coverage: {median_cov:.1f}× (below {self.MIN_MEDIAN_COVERAGE}×, may be low)"
            )

    def _check_uncovered_percentage(self, uncovered_pct: float) -> None:
        """
        Check uncovered percentage against threshold.

        Args:
            uncovered_pct: Percentage of VNTR region uncovered
        """
        if uncovered_pct <= self.MAX_UNCOVERED_PERCENT:
            log.success(f"✓ VNTR region coverage: {uncovered_pct:.2f}% uncovered (excellent)")
        else:
            log.warning(
                f"⚠ VNTR region coverage: {uncovered_pct:.2f}% uncovered "
                f"(above {self.MAX_UNCOVERED_PERCENT}%, may have gaps)"
            )
