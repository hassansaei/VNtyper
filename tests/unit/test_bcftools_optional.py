"""
Unit tests for bcftools optional functionality (Issue: bcftools dependency fix).

Tests the helper functions that make bcftools optional for IGV report generation:
- _try_compress_vcf_with_bcftools() in kestrel_genotyping.py
- _select_best_vcf_file() in pipeline.py

Test Coverage:
- bcftools availability detection
- Graceful fallback when bcftools unavailable
- VCF file selection logic (compressed vs uncompressed)
- Edge cases: missing files, command failures

Related:
- Commit: f002814 - fix: make bcftools optional for IGV report generation
- BCFTOOLS_FIX_COMPLETED.md
"""

import os
import tempfile
from unittest.mock import MagicMock, patch

import pytest

from vntyper.scripts.kestrel_genotyping import _try_compress_vcf_with_bcftools
from vntyper.scripts.pipeline import _select_best_vcf_file


@pytest.mark.unit
class TestTryCompressVcfWithBcftools:
    """Test suite for _try_compress_vcf_with_bcftools function."""

    def test_bcftools_not_available_returns_false(self):
        """
        When bcftools is not in PATH, should return False without attempting compression.

        This is the primary use case for the fix - tests without bcftools installed.
        """
        with patch('shutil.which', return_value=None):
            with tempfile.TemporaryDirectory() as tmpdir:
                input_vcf = os.path.join(tmpdir, "input.vcf")
                output_vcf_gz = os.path.join(tmpdir, "output.vcf.gz")

                # Create dummy input file
                with open(input_vcf, 'w') as f:
                    f.write("##fileformat=VCFv4.2\n")

                result = _try_compress_vcf_with_bcftools(input_vcf, output_vcf_gz, tmpdir)

                assert result is False, "Should return False when bcftools not available"
                assert not os.path.exists(output_vcf_gz), "Should not create output file"

    def test_bcftools_available_but_command_fails_returns_false(self):
        """
        When bcftools exists but command fails, should return False.

        Handles cases like corrupted VCF files or bcftools errors.
        """
        with patch('shutil.which', return_value='/usr/bin/bcftools'):
            with patch('vntyper.scripts.kestrel_genotyping.run_command', return_value=False):
                with tempfile.TemporaryDirectory() as tmpdir:
                    input_vcf = os.path.join(tmpdir, "input.vcf")
                    output_vcf_gz = os.path.join(tmpdir, "output.vcf.gz")

                    with open(input_vcf, 'w') as f:
                        f.write("##fileformat=VCFv4.2\n")

                    result = _try_compress_vcf_with_bcftools(input_vcf, output_vcf_gz, tmpdir)

                    assert result is False, "Should return False when bcftools command fails"

    def test_bcftools_available_and_succeeds_returns_true(self):
        """
        When bcftools is available and compression succeeds, should return True.

        This is the optimal case - bcftools installed and working.
        """
        with patch('shutil.which', return_value='/usr/bin/bcftools'):
            with patch('vntyper.scripts.kestrel_genotyping.run_command', return_value=True):
                with tempfile.TemporaryDirectory() as tmpdir:
                    input_vcf = os.path.join(tmpdir, "input.vcf")
                    output_vcf_gz = os.path.join(tmpdir, "output.vcf.gz")

                    with open(input_vcf, 'w') as f:
                        f.write("##fileformat=VCFv4.2\n")

                    result = _try_compress_vcf_with_bcftools(input_vcf, output_vcf_gz, tmpdir)

                    assert result is True, "Should return True when compression succeeds"

    def test_correct_bcftools_command_called(self):
        """Verify the correct bcftools command is constructed and executed."""
        with patch('shutil.which', return_value='/usr/bin/bcftools'):
            with patch('vntyper.scripts.kestrel_genotyping.run_command', return_value=True) as mock_run:
                with tempfile.TemporaryDirectory() as tmpdir:
                    input_vcf = os.path.join(tmpdir, "input.vcf")
                    output_vcf_gz = os.path.join(tmpdir, "output.vcf.gz")

                    with open(input_vcf, 'w') as f:
                        f.write("##fileformat=VCFv4.2\n")

                    _try_compress_vcf_with_bcftools(input_vcf, output_vcf_gz, tmpdir)

                    # Verify run_command was called with correct command
                    mock_run.assert_called_once()
                    call_args = mock_run.call_args[0][0]
                    assert "bcftools sort" in call_args
                    assert input_vcf in call_args
                    assert output_vcf_gz in call_args
                    assert "-O z" in call_args  # gzip compression


@pytest.mark.unit
class TestSelectBestVcfFile:
    """Test suite for _select_best_vcf_file function."""

    def test_compressed_vcf_exists_returns_compressed(self):
        """
        When .vcf.gz exists, should prefer compressed file (optimal case).
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            vcf_gz = os.path.join(tmpdir, "output_indel.vcf.gz")
            vcf = os.path.join(tmpdir, "output_indel.vcf")

            # Create both files
            with open(vcf_gz, 'w') as f:
                f.write("compressed")
            with open(vcf, 'w') as f:
                f.write("uncompressed")

            result = _select_best_vcf_file(tmpdir)

            assert result == vcf_gz, "Should prefer compressed .vcf.gz"

    def test_only_uncompressed_vcf_exists_returns_uncompressed(self):
        """
        When only .vcf exists (bcftools unavailable), should return uncompressed file.

        This is the fallback case when bcftools is not installed.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            vcf = os.path.join(tmpdir, "output_indel.vcf")

            # Create only uncompressed file
            with open(vcf, 'w') as f:
                f.write("uncompressed")

            result = _select_best_vcf_file(tmpdir)

            assert result == vcf, "Should return uncompressed .vcf when .vcf.gz missing"

    def test_neither_file_exists_returns_none(self):
        """
        When neither VCF file exists, should return None.

        Handles edge case where VCF generation failed completely.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # No VCF files created
            result = _select_best_vcf_file(tmpdir)

            assert result is None, "Should return None when no VCF files exist"

    def test_empty_directory_returns_none(self):
        """Empty kestrel directory should return None."""
        with tempfile.TemporaryDirectory() as tmpdir:
            kestrel_dir = os.path.join(tmpdir, "kestrel")
            os.makedirs(kestrel_dir)

            result = _select_best_vcf_file(kestrel_dir)

            assert result is None, "Should return None for empty directory"

    def test_preference_order_compressed_over_uncompressed(self):
        """
        Verify explicit preference: .vcf.gz > .vcf > None

        Even if uncompressed is newer/larger, compressed should win.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            vcf_gz = os.path.join(tmpdir, "output_indel.vcf.gz")
            vcf = os.path.join(tmpdir, "output_indel.vcf")

            # Create uncompressed first (older timestamp)
            with open(vcf, 'w') as f:
                f.write("uncompressed - created first")

            # Create compressed second (newer timestamp)
            with open(vcf_gz, 'w') as f:
                f.write("gz")

            result = _select_best_vcf_file(tmpdir)

            assert result == vcf_gz, "Should always prefer .vcf.gz regardless of timestamp"


@pytest.mark.unit
class TestBcftoolsIntegration:
    """Integration tests for bcftools optional workflow."""

    def test_full_workflow_bcftools_unavailable(self):
        """
        Test complete workflow when bcftools is not available.

        Workflow: compress fails → select uncompressed → IGV uses uncompressed
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            input_vcf = os.path.join(tmpdir, "output_indel.vcf")
            output_vcf_gz = os.path.join(tmpdir, "output_indel.vcf.gz")

            # Create input VCF
            with open(input_vcf, 'w') as f:
                f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\n")

            # Simulate bcftools unavailable
            with patch('shutil.which', return_value=None):
                compress_result = _try_compress_vcf_with_bcftools(
                    input_vcf, output_vcf_gz, tmpdir
                )

            assert compress_result is False, "Compression should fail"
            assert not os.path.exists(output_vcf_gz), "Compressed file should not exist"

            # Select best file should return uncompressed
            selected = _select_best_vcf_file(tmpdir)
            assert selected == input_vcf, "Should select uncompressed VCF"

    def test_full_workflow_bcftools_available(self):
        """
        Test complete workflow when bcftools is available and works.

        Workflow: compress succeeds → select compressed → IGV uses compressed
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            input_vcf = os.path.join(tmpdir, "output_indel.vcf")
            output_vcf_gz = os.path.join(tmpdir, "output_indel.vcf.gz")

            # Create input VCF
            with open(input_vcf, 'w') as f:
                f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\n")

            # Simulate bcftools available and succeeds
            with patch('shutil.which', return_value='/usr/bin/bcftools'):
                with patch('vntyper.scripts.kestrel_genotyping.run_command', return_value=True):
                    # Also create the output file to simulate successful compression
                    with open(output_vcf_gz, 'w') as f:
                        f.write("compressed")

                    compress_result = _try_compress_vcf_with_bcftools(
                        input_vcf, output_vcf_gz, tmpdir
                    )

            assert compress_result is True, "Compression should succeed"

            # Select best file should return compressed
            selected = _select_best_vcf_file(tmpdir)
            assert selected == output_vcf_gz, "Should select compressed VCF"


# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit
