#!/usr/bin/env python3
"""
test_reference_registry.py

Unit tests for the reference_registry module.

Tests cover:
- Assembly name normalization
- Coordinate system resolution
- Reference source resolution
- Metadata retrieval
- Validation functions
"""

import pytest

from vntyper.scripts.reference_registry import (
    # Data structures (for validation)
    ASSEMBLY_ALIASES,
    ASSEMBLY_METADATA,
    COORDINATE_SYSTEMS,
    REFERENCE_SOURCES,
    get_all_coordinate_systems,
    get_all_reference_sources,
    get_assembly_metadata,
    get_coordinate_system,
    get_coordinates,
    get_reference_source,
    get_reference_source_info,
    is_deprecated,
    list_assemblies,
    # Main functions
    normalize_assembly_name,
    resolve_chromosome_name,
    validate_registry,
)

# Mark all tests in this module as unit tests
pytestmark = pytest.mark.unit


# =============================================================================
# Test normalize_assembly_name()
# =============================================================================


class TestNormalizeAssemblyName:
    """Tests for normalize_assembly_name() function."""

    @pytest.mark.parametrize(
        "input_name,expected",
        [
            # All canonical names (should pass through)
            ("hg19", "hg19"),
            ("hg38", "hg38"),
            ("GRCh37", "GRCh37"),
            ("GRCh38", "GRCh38"),
            ("hg19_ncbi", "hg19_ncbi"),
            ("hg38_ncbi", "hg38_ncbi"),
            ("hg19_ensembl", "hg19_ensembl"),
            ("hg38_ensembl", "hg38_ensembl"),
        ],
    )
    def test_canonical_names(self, input_name, expected):
        """Test that canonical names pass through unchanged."""
        result = normalize_assembly_name(input_name, warn_deprecated=False)
        assert result == expected

    def test_invalid_assembly_raises_error(self):
        """Test that invalid assembly names raise ValueError."""
        with pytest.raises(ValueError, match="Unknown assembly"):
            normalize_assembly_name("invalid_assembly")


# =============================================================================
# Test get_coordinate_system()
# =============================================================================


class TestGetCoordinateSystem:
    """Tests for get_coordinate_system() function."""

    @pytest.mark.parametrize(
        "assembly_name,expected_coord_system",
        [
            ("hg19", "GRCh37"),
            ("hg38", "GRCh38"),
            ("GRCh37", "GRCh37"),
            ("GRCh38", "GRCh38"),
            ("hg19_ncbi", "GRCh37"),
            ("hg38_ncbi", "GRCh38"),
            ("hg19_ensembl", "GRCh37"),
            ("hg38_ensembl", "GRCh38"),
        ],
    )
    def test_coordinate_system_resolution(self, assembly_name, expected_coord_system):
        """Test that assemblies resolve to correct coordinate systems."""
        result = get_coordinate_system(assembly_name)
        assert result == expected_coord_system

    def test_invalid_assembly_raises_error(self):
        """Test that invalid assembly raises ValueError."""
        with pytest.raises(ValueError):
            get_coordinate_system("invalid_assembly")


# =============================================================================
# Test get_reference_source()
# =============================================================================


class TestGetReferenceSource:
    """Tests for get_reference_source() function."""

    @pytest.mark.parametrize(
        "assembly_name,expected_source",
        [
            ("hg19", "ucsc"),
            ("hg38", "ucsc"),
            ("GRCh37", "ncbi"),
            ("GRCh38", "ncbi"),
            ("hg19_ncbi", "ncbi"),
            ("hg38_ncbi", "ncbi"),
            ("hg19_ensembl", "ensembl"),
            ("hg38_ensembl", "ensembl"),
        ],
    )
    def test_reference_source_resolution(self, assembly_name, expected_source):
        """Test that assemblies resolve to correct reference sources."""
        result = get_reference_source(assembly_name)
        assert result == expected_source

    def test_invalid_assembly_raises_error(self):
        """Test that invalid assembly raises ValueError."""
        with pytest.raises(ValueError):
            get_reference_source("invalid_assembly")


# =============================================================================
# Test get_coordinates()
# =============================================================================


class TestGetCoordinates:
    """Tests for get_coordinates() function."""

    @pytest.mark.parametrize(
        "assembly_name,region_type,expected_coords",
        [
            # GRCh37 assemblies
            ("hg19", "bam_region_coords", "155158000-155163000"),
            ("hg19", "vntr_region_coords", "155160500-155162000"),
            ("GRCh37", "bam_region_coords", "155158000-155163000"),
            ("GRCh37", "vntr_region_coords", "155160500-155162000"),
            ("hg19_ncbi", "bam_region_coords", "155158000-155163000"),
            ("hg19_ensembl", "vntr_region_coords", "155160500-155162000"),
            # GRCh38 assemblies
            ("hg38", "bam_region_coords", "155184000-155194000"),
            ("hg38", "vntr_region_coords", "155188000-155192500"),
            ("GRCh38", "bam_region_coords", "155184000-155194000"),
            ("GRCh38", "vntr_region_coords", "155188000-155192500"),
            ("hg38_ncbi", "bam_region_coords", "155184000-155194000"),
            ("hg38_ensembl", "vntr_region_coords", "155188000-155192500"),
        ],
    )
    def test_coordinate_resolution(self, assembly_name, region_type, expected_coords):
        """Test that coordinates are correctly retrieved."""
        result = get_coordinates(assembly_name, region_type)
        assert result == expected_coords

    def test_invalid_region_type_raises_error(self):
        """Test that invalid region type raises ValueError."""
        with pytest.raises(ValueError, match="Region type"):
            get_coordinates("hg19", "invalid_region_type")

    def test_invalid_assembly_raises_error(self):
        """Test that invalid assembly raises ValueError."""
        with pytest.raises(ValueError):
            get_coordinates("invalid_assembly", "bam_region_coords")


# =============================================================================
# Test get_assembly_metadata()
# =============================================================================


class TestGetAssemblyMetadata:
    """Tests for get_assembly_metadata() function."""

    def test_metadata_structure(self):
        """Test that metadata has expected structure."""
        metadata = get_assembly_metadata("hg19")

        assert "coordinate_system" in metadata
        assert "reference_source" in metadata
        assert "description" in metadata
        assert "deprecated" in metadata

    @pytest.mark.parametrize(
        "assembly_name",
        ["hg19", "hg38", "GRCh37", "GRCh38", "hg19_ncbi", "hg38_ncbi", "hg19_ensembl", "hg38_ensembl"],
    )
    def test_metadata_retrieval(self, assembly_name):
        """Test that metadata can be retrieved for all assemblies."""
        metadata = get_assembly_metadata(assembly_name)
        assert metadata is not None
        assert metadata["deprecated"] is False

    def test_metadata_is_copy(self):
        """Test that returned metadata is a copy, not reference."""
        metadata1 = get_assembly_metadata("hg19")
        metadata2 = get_assembly_metadata("hg19")

        # Modify one copy
        metadata1["test_key"] = "test_value"

        # Other copy should be unchanged
        assert "test_key" not in metadata2

    def test_invalid_assembly_raises_error(self):
        """Test that invalid assembly raises ValueError."""
        with pytest.raises(ValueError):
            get_assembly_metadata("invalid_assembly")


# =============================================================================
# Test get_reference_source_info()
# =============================================================================


class TestGetReferenceSourceInfo:
    """Tests for get_reference_source_info() function."""

    @pytest.mark.parametrize(
        "source_name,expected_format",
        [
            ("ucsc", "chr{n}"),
            ("ncbi", "NC_{accession}"),
            ("ensembl", "{n}"),
        ],
    )
    def test_source_info_retrieval(self, source_name, expected_format):
        """Test that source info is correctly retrieved."""
        info = get_reference_source_info(source_name)
        assert info["chr_format"] == expected_format

    def test_source_info_is_copy(self):
        """Test that returned info is a copy, not reference."""
        info1 = get_reference_source_info("ucsc")
        info2 = get_reference_source_info("ucsc")

        # Modify one copy
        info1["test_key"] = "test_value"

        # Other copy should be unchanged
        assert "test_key" not in info2

    def test_invalid_source_raises_error(self):
        """Test that invalid source raises ValueError."""
        with pytest.raises(ValueError, match="Unknown reference source"):
            get_reference_source_info("invalid_source")


# =============================================================================
# Test list_assemblies()
# =============================================================================


class TestListAssemblies:
    """Tests for list_assemblies() function."""

    def test_list_all_canonical_names(self):
        """Test that all canonical names are returned."""
        assemblies = list_assemblies(include_deprecated=False)

        assert "hg19" in assemblies
        assert "hg38" in assemblies
        assert "GRCh37" in assemblies
        assert "GRCh38" in assemblies
        assert "hg19_ncbi" in assemblies
        assert "hg38_ncbi" in assemblies
        assert "hg19_ensembl" in assemblies
        assert "hg38_ensembl" in assemblies
        # Should have exactly 8 assemblies
        assert len(assemblies) == 8

    def test_list_is_sorted(self):
        """Test that returned list is sorted."""
        assemblies = list_assemblies()
        assert assemblies == sorted(assemblies)


# =============================================================================
# Test is_deprecated()
# =============================================================================


class TestIsDeprecated:
    """Tests for is_deprecated() function."""

    @pytest.mark.parametrize(
        "assembly_name",
        ["hg19", "hg38", "GRCh37", "GRCh38", "hg19_ncbi", "hg38_ncbi", "hg19_ensembl", "hg38_ensembl"],
    )
    def test_no_deprecated_assemblies(self, assembly_name):
        """Test that no assemblies are deprecated (all are canonical)."""
        assert is_deprecated(assembly_name) is False


# =============================================================================
# Test get_all_coordinate_systems()
# =============================================================================


class TestGetAllCoordinateSystems:
    """Tests for get_all_coordinate_systems() function."""

    def test_coordinate_systems_list(self):
        """Test that all coordinate systems are returned."""
        systems = get_all_coordinate_systems()

        assert "GRCh37" in systems
        assert "GRCh38" in systems
        assert len(systems) == 2


# =============================================================================
# Test get_all_reference_sources()
# =============================================================================


class TestGetAllReferenceSources:
    """Tests for get_all_reference_sources() function."""

    def test_reference_sources_list(self):
        """Test that all reference sources are returned."""
        sources = get_all_reference_sources()

        assert "ucsc" in sources
        assert "ncbi" in sources
        assert "ensembl" in sources
        assert len(sources) == 3


# =============================================================================
# Test resolve_chromosome_name()
# =============================================================================


class TestResolveChromosomeName:
    """Tests for resolve_chromosome_name() function."""

    @pytest.mark.parametrize(
        "assembly_name,chromosome_number,expected",
        [
            # UCSC naming
            ("hg19", 1, "chr1"),
            ("hg19", 2, "chr2"),
            ("hg19", 23, "chrX"),
            ("hg19", 24, "chrY"),
            ("hg19", 25, "chrM"),
            # ENSEMBL naming
            ("hg19_ensembl", 1, "1"),
            ("hg19_ensembl", 2, "2"),
            ("hg19_ensembl", 23, "X"),
            ("hg19_ensembl", 24, "Y"),
            ("hg19_ensembl", 25, "MT"),
            # NCBI naming (chr1 only implemented)
            ("hg19_ncbi", 1, "NC_000001.10"),
            ("hg38_ncbi", 1, "NC_000001.11"),
        ],
    )
    def test_chromosome_name_resolution(self, assembly_name, chromosome_number, expected):
        """Test that chromosome names are correctly resolved."""
        result = resolve_chromosome_name(assembly_name, chromosome_number)
        assert result == expected

    def test_detected_convention_override(self):
        """Test that detected_convention overrides assembly default."""
        # hg19 normally uses UCSC, but override with ensembl
        result = resolve_chromosome_name("hg19", 1, detected_convention="ensembl")
        assert result == "1"

        # hg19_ensembl normally uses ENSEMBL, but override with ucsc
        result = resolve_chromosome_name("hg19_ensembl", 1, detected_convention="ucsc")
        assert result == "chr1"

    def test_invalid_chromosome_number_raises_error(self):
        """Test that invalid chromosome numbers raise ValueError."""
        with pytest.raises(ValueError, match="Invalid chromosome number"):
            resolve_chromosome_name("hg19", 0)

        with pytest.raises(ValueError, match="Invalid chromosome number"):
            resolve_chromosome_name("hg19", 26)

    def test_ncbi_chr2_not_implemented(self):
        """Test that NCBI naming for chr2+ raises NotImplementedError."""
        with pytest.raises(NotImplementedError):
            resolve_chromosome_name("hg19_ncbi", 2)


# =============================================================================
# Test validate_registry()
# =============================================================================


class TestValidateRegistry:
    """Tests for validate_registry() function."""

    def test_registry_is_valid(self):
        """Test that the registry is internally consistent."""
        is_valid, errors = validate_registry()

        if not is_valid:
            # Print errors for debugging
            for error in errors:
                print(f"  - {error}")

        assert is_valid, f"Registry validation failed with {len(errors)} errors"
        assert len(errors) == 0


# =============================================================================
# Test Data Structure Integrity
# =============================================================================


class TestDataStructures:
    """Tests for data structure integrity and consistency."""

    def test_all_canonical_names_in_metadata(self):
        """Test that all canonical names have metadata."""
        canonical_names = set(ASSEMBLY_METADATA.keys())

        for alias, canonical in ASSEMBLY_ALIASES.items():
            if canonical not in canonical_names:
                # This is only allowed if it's a self-reference
                assert canonical == alias, f"Canonical name '{canonical}' not in ASSEMBLY_METADATA"

    def test_all_coordinate_systems_referenced(self):
        """Test that all coordinate systems are used by at least one assembly."""
        used_coord_systems = set()

        for metadata in ASSEMBLY_METADATA.values():
            used_coord_systems.add(metadata["coordinate_system"])

        for coord_system in COORDINATE_SYSTEMS.keys():
            assert coord_system in used_coord_systems, f"Coordinate system '{coord_system}' is unused"

    def test_all_reference_sources_referenced(self):
        """Test that all reference sources are used by at least one assembly."""
        used_sources = set()

        for metadata in ASSEMBLY_METADATA.values():
            used_sources.add(metadata["reference_source"])

        for source in REFERENCE_SOURCES.keys():
            assert source in used_sources, f"Reference source '{source}' is unused"


    def test_coordinate_systems_have_required_fields(self):
        """Test that coordinate systems have all required fields."""
        required_fields = {"chromosome", "bam_region_coords", "vntr_region_coords"}

        for coord_system, data in COORDINATE_SYSTEMS.items():
            for field in required_fields:
                assert field in data, f"Coordinate system '{coord_system}' missing field '{field}'"

    def test_assembly_metadata_has_required_fields(self):
        """Test that assembly metadata has all required fields."""
        required_fields = {"coordinate_system", "reference_source", "description", "deprecated"}

        for assembly, metadata in ASSEMBLY_METADATA.items():
            for field in required_fields:
                assert field in metadata, f"Assembly '{assembly}' missing field '{field}'"

    def test_reference_sources_have_required_fields(self):
        """Test that reference sources have all required fields."""
        required_fields = {"description", "chr_format", "example", "download_base_url"}

        for source, data in REFERENCE_SOURCES.items():
            for field in required_fields:
                assert field in data, f"Reference source '{source}' missing field '{field}'"


# =============================================================================
# Integration Tests
# =============================================================================


class TestIntegration:
    """Integration tests for common workflows."""

    def test_full_workflow_canonical_name(self):
        """Test full workflow with canonical assembly name."""
        assembly = "hg19_ensembl"

        # Normalize (should pass through)
        canonical = normalize_assembly_name(assembly, warn_deprecated=False)
        assert canonical == "hg19_ensembl"

        # Get coordinate system
        coord_system = get_coordinate_system(canonical)
        assert coord_system == "GRCh37"

        # Get reference source
        ref_source = get_reference_source(canonical)
        assert ref_source == "ensembl"

        # Get coordinates
        coords = get_coordinates(canonical, "bam_region_coords")
        assert coords == "155158000-155163000"

        # Resolve chromosome name
        chr_name = resolve_chromosome_name(canonical, 1)
        assert chr_name == "1"

    def test_coordinate_consistency_across_sources(self):
        """Test that all GRCh37 assemblies share the same coordinates."""
        grch37_assemblies = ["hg19", "hg19_ncbi", "hg19_ensembl", "GRCh37"]

        expected_coords = "155158000-155163000"

        for assembly in grch37_assemblies:
            coords = get_coordinates(assembly, "bam_region_coords")
            assert coords == expected_coords, f"Assembly '{assembly}' has incorrect coordinates"

    def test_coordinate_consistency_grch38(self):
        """Test that all GRCh38 assemblies share the same coordinates."""
        grch38_assemblies = ["hg38", "hg38_ncbi", "hg38_ensembl", "GRCh38"]

        expected_coords = "155184000-155194000"

        for assembly in grch38_assemblies:
            coords = get_coordinates(assembly, "bam_region_coords")
            assert coords == expected_coords, f"Assembly '{assembly}' has incorrect coordinates"
