"""The single mapping from an assembly label to the config keys that name its files.

Eight labels denote six physical files: `GRCh37`/`hg19_ncbi` are one NCBI file and
`GRCh38`/`hg38_ncbi` are another. Keying on the label cannot work, because
`install_references_config.json` only ever produces a `GRCh37` and a `GRCh38` entry -
a label-keyed writer would never emit `bwa_reference_hg38_ncbi`, and that run would
silently fall back to the UCSC FASTA.
"""

import pytest

from vntyper.scripts.reference_registry import (
    list_assemblies,
    physical_reference_id,
    reference_keys,
    ucsc_family,
)

pytestmark = pytest.mark.unit

PHYSICAL_ID = {
    "hg19": "hg19",
    "hg38": "hg38",
    "GRCh37": "GRCh37",
    "hg19_ncbi": "GRCh37",
    "GRCh38": "GRCh38",
    "hg38_ncbi": "GRCh38",
    "hg19_ensembl": "hg19_ensembl",
    "hg38_ensembl": "hg38_ensembl",
}

UCSC_FAMILY = {
    "hg19": "hg19",
    "GRCh37": "hg19",
    "hg19_ncbi": "hg19",
    "hg19_ensembl": "hg19",
    "hg38": "hg38",
    "GRCh38": "hg38",
    "hg38_ncbi": "hg38",
    "hg38_ensembl": "hg38",
}


class TestPhysicalIdentity:
    @pytest.mark.parametrize("label,expected", sorted(PHYSICAL_ID.items()))
    def test_each_label_maps_to_its_physical_file(self, label, expected):
        assert physical_reference_id(label) == expected

    def test_the_two_ncbi_aliases_collapse_to_one_file(self):
        assert physical_reference_id("GRCh38") == physical_reference_id("hg38_ncbi")
        assert physical_reference_id("GRCh37") == physical_reference_id("hg19_ncbi")

    def test_every_registry_label_is_covered(self):
        """A label added to the registry without a physical id is a silent wrong file."""
        assert set(list_assemblies()) == set(PHYSICAL_ID)

    def test_there_are_exactly_six_physical_files(self):
        assert len(set(PHYSICAL_ID.values())) == 6


class TestUcscFamily:
    @pytest.mark.parametrize("label,expected", sorted(UCSC_FAMILY.items()))
    def test_family_is_the_ucsc_name_of_the_coordinate_system(self, label, expected):
        assert ucsc_family(label) == expected


class TestPhysicalKeyedKinds:
    """bwa and cram: contig naming differs by source, so the source is load-bearing."""

    @pytest.mark.parametrize("kind,prefix", [("bwa", "bwa_reference"), ("cram", "cram_reference")])
    def test_exact_key_first_then_ucsc_family(self, kind, prefix):
        assert reference_keys(kind, "hg38_ensembl") == (
            f"{prefix}_hg38_ensembl",
            f"{prefix}_hg38",
        )

    @pytest.mark.parametrize("kind,prefix", [("bwa", "bwa_reference"), ("cram", "cram_reference")])
    def test_the_existing_label_override_contract_survives(self, kind, prefix):
        """test_reference_resolution.py pins that a config may specialise hg19_ncbi."""
        assert reference_keys(kind, "hg19_ncbi")[0] == f"{prefix}_hg19_ncbi"

    @pytest.mark.parametrize("label", ["GRCh38", "hg38_ncbi"])
    def test_both_ncbi_labels_reach_the_same_physical_key(self, label):
        assert "bwa_reference_GRCh38" in reference_keys("bwa", label)

    def test_a_label_specific_key_is_offered_before_the_physical_one(self):
        assert reference_keys("bwa", "hg38_ncbi") == (
            "bwa_reference_hg38_ncbi",
            "bwa_reference_GRCh38",
            "bwa_reference_hg38",
        )

    def test_a_ucsc_label_yields_a_single_key_not_a_duplicate_pair(self):
        assert reference_keys("bwa", "hg38") == ("bwa_reference_hg38",)


class TestCoordinateKeyedKinds:
    """advntr and shark: only two databases and two regions exist, by coordinate system."""

    @pytest.mark.parametrize(
        "label,expected",
        [(label, UCSC_FAMILY[label]) for label in sorted(UCSC_FAMILY)],
    )
    def test_advntr_key_follows_the_coordinate_system(self, label, expected):
        assert reference_keys("advntr", label) == (f"advntr_reference_vntr_{expected}",)

    @pytest.mark.parametrize(
        "label,expected",
        [(label, UCSC_FAMILY[label]) for label in sorted(UCSC_FAMILY)],
    )
    def test_shark_key_follows_the_coordinate_system(self, label, expected):
        assert reference_keys("shark", label) == (f"muc1_region_fasta_{expected}",)

    @pytest.mark.parametrize("label", ["hg38_ncbi", "hg38_ensembl"])
    def test_the_four_labels_that_used_to_fall_through_now_reach_grch38(self, label):
        """pipeline.py's raw dict sent these to the hg19 database. Regression for F2."""
        assert reference_keys("advntr", label) == ("advntr_reference_vntr_hg38",)


class TestRejection:
    def test_an_unknown_kind_is_rejected_by_name(self):
        with pytest.raises(ValueError, match="minimap2"):
            reference_keys("minimap2", "hg19")

    def test_an_unknown_assembly_is_rejected_with_the_supported_list(self):
        with pytest.raises(ValueError, match="hg38_ensembl"):
            reference_keys("bwa", "b37")
