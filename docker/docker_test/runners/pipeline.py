"""Pipeline test runner - orchestrates BAM and adVNTR tests."""

import shutil
import time
from pathlib import Path
from typing import Optional

from ..config import TestConfig
from ..container import DockerExecutionError, DockerOperations, DockerTimeoutError
from ..utils.logging import log
from ..validators.advntr import adVNTRValidator
from ..validators.coverage import CoverageValidator
from ..validators.kestrel import KestrelValidator


class PipelineRunner:
    """Run and validate VNtyper pipeline tests."""

    def __init__(self, config: TestConfig, docker_ops: DockerOperations):
        self.config = config
        self.docker = docker_ops
        # Resolve ALL paths at initialization to completely avoid getcwd() calls later
        import os
        self.original_cwd = Path(os.getcwd()).resolve()
        self.test_data_dir = config.test_data_dir.resolve()
        self.output_base_dir = config.docker_test_dir.resolve()

    def run_bam_test(self, test_name: str) -> bool:
        """
        Run single BAM test.

        Args:
            test_name: Test case name (e.g., "example_b178_hg19_subset_fast")

        Returns:
            True if test passes
        """
        log.section(f"Test Case: {test_name}")

        # Get test configuration
        test_config = self.config.get_bam_test(test_name)

        # Setup paths from config - use pre-resolved base paths to avoid getcwd() issues
        bam_path_str = test_config["bam"]
        bam_file = (
            Path(bam_path_str)
            if Path(bam_path_str).is_absolute()
            else Path(self.original_cwd) / bam_path_str
        )
        output_dir = self.output_base_dir / test_name
        reference_assembly = test_config.get("reference_assembly", "hg19")

        if not bam_file.exists():
            log.error(f"BAM file not found: {bam_file}")
            return False

        # Clean previous output
        if output_dir.exists():
            try:
                shutil.rmtree(output_dir)
            except Exception as e:
                log.warning(f"Failed to remove old output dir: {e}")

        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            log.error(f"Failed to create output directory {output_dir}: {e}")
            raise

        log.info(f"BAM: {bam_file.name}")
        log.info(f"Reference: {reference_assembly}")
        log.info(f"Output: {output_dir}")

        # Run pipeline
        start_time = time.time()
        try:
            self.docker.run_pipeline(
                bam_file=bam_file,
                output_dir=output_dir,
                reference_assembly=reference_assembly,
                threads=4
            )
            duration = time.time() - start_time
            log.success(f"Pipeline completed in {duration:.1f}s")
        except (DockerExecutionError, DockerTimeoutError) as e:
            log.error(f"Pipeline failed: {e}")
            return False

        # Validate output files exist
        result_dir = output_dir / "result"
        required_files = [
            "summary_report.html",
            "kestrel/kestrel_result.tsv",
            "coverage/coverage_summary.tsv",
            "pipeline_summary.json",
        ]

        # Check if IGV report is expected
        check_igv_report = test_config.get("check_igv_report", False)
        if check_igv_report:
            required_files.append("igv_report.html")

        for file_path in required_files:
            full_path = result_dir / file_path
            if not full_path.exists() or full_path.stat().st_size == 0:
                log.error(f"Missing or empty: {file_path}")
                return False
            log.success(f"✓ Found: {file_path}")

        # Validate results
        try:
            # Kestrel validation
            kestrel_assertions = test_config.get("kestrel_assertions")
            if kestrel_assertions:
                validator = KestrelValidator(result_dir, kestrel_assertions)
                validator.validate()

            # Coverage validation
            coverage_validator = CoverageValidator(result_dir)
            coverage_validator.validate()

        except Exception as e:
            log.error(f"Validation failed: {e}")
            return False

        log.success(f"Test {test_name} PASSED")
        return True

    def run_advntr_test(self, test_name: str) -> bool:
        """
        Run adVNTR module test.

        Args:
            test_name: adVNTR test name (e.g., "example_a5c1_hg19_subset_advntr")

        Returns:
            True if test passes
        """
        log.section(f"adVNTR Test: {test_name}")

        # Get test configuration
        test_config = self.config.get_advntr_test(test_name)

        # Setup paths from config - use pre-resolved base paths to avoid getcwd() issues
        bam_path_str = test_config["bam"]
        bam_file = (
            Path(bam_path_str)
            if Path(bam_path_str).is_absolute()
            else Path(self.original_cwd) / bam_path_str
        )
        output_dir = self.output_base_dir / test_name
        reference_assembly = test_config.get("reference_assembly", "hg19")

        if not bam_file.exists():
            log.error(f"BAM file not found: {bam_file}")
            return False

        # Clean previous output
        if output_dir.exists():
            try:
                shutil.rmtree(output_dir)
            except Exception as e:
                log.warning(f"Failed to remove old output dir: {e}")

        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            log.error(f"Failed to create output directory {output_dir}: {e}")
            raise

        log.info(f"BAM: {bam_file.name}")
        log.info(f"Reference: {reference_assembly}")
        log.info(f"Output: {output_dir}")
        log.info("Running with --extra-modules advntr")

        # Run pipeline with adVNTR
        start_time = time.time()
        try:
            self.docker.run_pipeline(
                bam_file=bam_file,
                output_dir=output_dir,
                reference_assembly=reference_assembly,
                threads=4,
                extra_modules=["advntr"],
            )
            duration = time.time() - start_time
            log.success(f"Pipeline with adVNTR completed in {duration:.1f}s")
        except (DockerExecutionError, DockerTimeoutError) as e:
            log.error(f"Pipeline failed: {e}")
            return False

        # Validate adVNTR output
        result_dir = output_dir / "result"
        advntr_assertions = test_config.get("advntr_assertions", {})
        try:
            validator = adVNTRValidator(result_dir, advntr_assertions)
            validator.validate()
        except Exception as e:
            log.error(f"adVNTR validation failed: {e}")
            return False

        log.success(f"adVNTR test {test_name} PASSED")
        return True

    def run_multiple_tests(self, test_ids: Optional[list[str]] = None) -> dict[str, bool]:
        """
        Run multiple BAM tests.

        Args:
            test_ids: List of test IDs, or None for all

        Returns:
            Dict mapping test_id -> pass/fail
        """
        if test_ids is None:
            test_ids = self.config.list_bam_test_ids()

        results = {}
        passed = 0
        failed = 0

        for i, test_id in enumerate(test_ids, 1):
            log.info("=" * 60)
            log.info(f"Running test {i}/{len(test_ids)}: {test_id}")
            log.info("=" * 60)

            try:
                success = self.run_bam_test(test_id)
                results[test_id] = success
                if success:
                    passed += 1
                else:
                    failed += 1
            except Exception as e:
                import traceback
                log.error(f"Test {test_id} crashed: {e}")
                log.error(f"Traceback: {traceback.format_exc()}")
                results[test_id] = False
                failed += 1

            print()

        # Summary
        log.info("=" * 60)
        log.info("Test Summary")
        log.info("=" * 60)
        log.info(f"Total: {len(test_ids)}")
        log.success(f"Passed: {passed}")
        if failed > 0:
            log.error(f"Failed: {failed}")
        else:
            log.success(f"Failed: {failed}")

        return results
