#!/usr/bin/env python3
"""
Unit tests for region_utils.py

Tests region string construction and chromosome name caching functionality.
"""

from unittest.mock import patch

import pytest

from vntyper.scripts.region_utils import (
    build_region_string,
    clear_chromosome_cache,
    get_cache_info,
    get_region_string,
    get_region_string_with_fallback,
    resolve_assembly_alias,
)

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit


class TestBuildRegionString:
    """Test basic region string construction."""

    def test_build_ucsc_region(self):
        """Test building region string with UCSC chromosome name."""
        result = build_region_string("chr1", "155158000-155163000")
        assert result == "chr1:155158000-155163000"

    def test_build_ncbi_region(self):
        """Test building region string with NCBI accession."""
        result = build_region_string("NC_000001.10", "155158000-155163000")
        assert result == "NC_000001.10:155158000-155163000"

    def test_build_simple_region(self):
        """Test building region string with simple numeric chromosome."""
        result = build_region_string("1", "155158000-155163000")
        assert result == "1:155158000-155163000"

    def test_empty_chromosome_name(self):
        """Test error handling for empty chromosome name."""
        with pytest.raises(ValueError, match="Invalid inputs"):
            build_region_string("", "155158000-155163000")

    def test_empty_coordinates(self):
        """Test error handling for empty coordinates."""
        with pytest.raises(ValueError, match="Invalid inputs"):
            build_region_string("chr1", "")

    def test_invalid_coordinate_format(self):
        """Test error handling for invalid coordinate format."""
        with pytest.raises(ValueError, match="Invalid coordinate format"):
            build_region_string("chr1", "155158000")  # Missing hyphen

    def test_invalid_coordinate_values(self):
        """Test error handling for non-numeric coordinates."""
        with pytest.raises(ValueError, match="Invalid coordinate values"):
            build_region_string("chr1", "abc-def")


class TestResolveAssemblyAlias:
    """Test assembly alias resolution to coordinate systems."""

    def test_resolve_hg19_to_grch37(self):
        """Test resolving hg19 to GRCh37 coordinate system."""
        assert resolve_assembly_alias("hg19") == "GRCh37"

    def test_resolve_grch37(self):
        """Test resolving GRCh37 (no change)."""
        assert resolve_assembly_alias("GRCh37") == "GRCh37"

    def test_resolve_hg38_to_grch38(self):
        """Test resolving hg38 to GRCh38 coordinate system."""
        assert resolve_assembly_alias("hg38") == "GRCh38"

    def test_resolve_grch38(self):
        """Test resolving GRCh38 (no change)."""
        assert resolve_assembly_alias("GRCh38") == "GRCh38"

    def test_unknown_assembly(self):
        """Test handling of unknown assembly (defaults to GRCh37)."""
        assert resolve_assembly_alias("unknown") == "GRCh37"


class TestGetRegionString:
    """Test dynamic region string generation."""

    @patch('vntyper.scripts.chromosome_utils.get_chromosome_name_from_bam')
    def test_get_region_hg19_ucsc(self, mock_get_chr):
        """Test getting region string for hg19 with UCSC naming."""
        mock_get_chr.return_value = "chr1"

        config = {
            "bam_processing": {
                "assemblies": {
                    "GRCh37": {
                        "bam_region_coords": "155158000-155163000",
                        "vntr_region_coords": "155160500-155162000",
                        "chromosome": 1
                    }
                }
            }
        }

        result = get_region_string(
            "test.bam", "hg19", "bam_region_coords", config
        )
        assert result == "chr1:155158000-155163000"

    @patch('vntyper.scripts.chromosome_utils.get_chromosome_name_from_bam')
    def test_get_region_grch37_ncbi(self, mock_get_chr):
        """Test getting region string for GRCh37 with NCBI naming."""
        mock_get_chr.return_value = "NC_000001.10"

        config = {
            "bam_processing": {
                "assemblies": {
                    "GRCh37": {
                        "bam_region_coords": "155158000-155163000",
                        "chromosome": 1
                    }
                },
                "known_chromosome_naming": {
                    "GRCh37": {"ncbi": "NC_000001.10"}
                }
            }
        }

        result = get_region_string(
            "test.bam", "GRCh37", "bam_region_coords", config
        )
        assert result == "NC_000001.10:155158000-155163000"

    @patch('vntyper.scripts.chromosome_utils.get_chromosome_name_from_bam')
    def test_get_region_hg38(self, mock_get_chr):
        """Test getting region string for hg38."""
        mock_get_chr.return_value = "chr1"

        config = {
            "bam_processing": {
                "assemblies": {
                    "GRCh38": {
                        "bam_region_coords": "155184000-155194000",
                        "chromosome": 1
                    }
                }
            }
        }

        result = get_region_string(
            "test.bam", "hg38", "bam_region_coords", config
        )
        assert result == "chr1:155184000-155194000"

    @patch('vntyper.scripts.chromosome_utils.get_chromosome_name_from_bam')
    def test_caching(self, mock_get_chr):
        """Test that chromosome names are cached."""
        mock_get_chr.return_value = "chr1"

        config = {
            "bam_processing": {
                "assemblies": {
                    "GRCh37": {
                        "bam_region_coords": "155158000-155163000",
                        "vntr_region_coords": "155160500-155162000",
                        "chromosome": 1
                    }
                }
            }
        }

        # Clear cache first
        clear_chromosome_cache()

        # First call should invoke get_chromosome_name_from_bam
        result1 = get_region_string(
            "test.bam", "hg19", "bam_region_coords", config
        )
        assert result1 == "chr1:155158000-155163000"
        assert mock_get_chr.call_count == 1

        # Second call should use cached value
        result2 = get_region_string(
            "test.bam", "hg19", "vntr_region_coords", config
        )
        assert result2 == "chr1:155160500-155162000"
        assert mock_get_chr.call_count == 1  # Still 1, used cache

    def test_missing_assembly_config(self):
        """Test error when assembly not found in config."""
        config = {"bam_processing": {}}

        with pytest.raises(KeyError, match="Configuration missing"):
            get_region_string(
                "test.bam", "hg19", "bam_region_coords", config
            )

    def test_missing_region_type(self):
        """Test error when region type not found."""
        config = {
            "bam_processing": {
                "assemblies": {
                    "GRCh37": {
                        "chromosome": 1
                    }
                }
            }
        }

        with pytest.raises(KeyError, match="Region type"):
            get_region_string(
                "test.bam", "hg19", "nonexistent_region", config
            )


class TestGetRegionStringWithFallback:
    """Test region string resolution with fallback to legacy format."""

    @patch('vntyper.scripts.region_utils.get_region_string')
    def test_new_format_success(self, mock_get_region):
        """Test successful resolution with new format."""
        mock_get_region.return_value = "chr1:155158000-155163000"

        config = {
            "bam_processing": {
                "assemblies": {
                    "GRCh37": {
                        "bam_region_coords": "155158000-155163000",
                        "chromosome": 1
                    }
                }
            }
        }

        result = get_region_string_with_fallback(
            "test.bam", "hg19", "bam_region", config
        )
        assert result == "chr1:155158000-155163000"

    @patch('vntyper.scripts.region_utils.get_region_string')
    def test_fallback_to_legacy(self, mock_get_region):
        """Test fallback to legacy config format."""
        mock_get_region.side_effect = KeyError("Missing config")

        config = {
            "bam_processing": {
                "bam_region_hg19": "chr1:155158000-155163000"
            }
        }

        result = get_region_string_with_fallback(
            "test.bam", "hg19", "bam_region", config
        )
        assert result == "chr1:155158000-155163000"

    @patch('vntyper.scripts.region_utils.get_region_string')
    def test_both_formats_fail(self, mock_get_region):
        """Test error when both new and legacy formats fail."""
        mock_get_region.side_effect = KeyError("Missing config")

        config = {"bam_processing": {}}

        with pytest.raises(ValueError, match="Neither new nor legacy format"):
            get_region_string_with_fallback(
                "test.bam", "hg19", "bam_region", config
            )


class TestCacheManagement:
    """Test chromosome name cache management."""

    def test_clear_cache(self):
        """Test clearing the cache."""
        # First, add something to cache by calling get_cache_info
        clear_chromosome_cache()
        info = get_cache_info()
        assert info["size"] == 0

    def test_get_cache_info(self):
        """Test getting cache information."""
        clear_chromosome_cache()
        info = get_cache_info()
        assert "size" in info
        assert "entries" in info
        assert isinstance(info["entries"], list)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
