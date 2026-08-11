"""Which adVNTR database a run uses, for every assembly label.

pipeline.py used a four-entry dict with `.get(label, "hg19")`. The four
source-qualified labels were absent, so hg38_ncbi and hg38_ensembl silently loaded the
GRCh37 database - the wrong coordinate system, with no warning. Genotype-affecting.
"""

import pytest

from vntyper.scripts.reference_registry import list_assemblies
from vntyper.scripts.reference_resolution import resolve_from_mapping

pytestmark = pytest.mark.unit

EXPECTED = {
    "hg19": "hg19",
    "GRCh37": "hg19",
    "hg19_ncbi": "hg19",
    "hg19_ensembl": "hg19",
    "hg38": "hg38",
    "GRCh38": "hg38",
    "hg38_ncbi": "hg38",
    "hg38_ensembl": "hg38",
}


@pytest.mark.parametrize("label", sorted(EXPECTED))
def test_every_label_reaches_its_own_coordinate_systems_database(label):
    mapping = {
        "advntr_reference_vntr_hg19": "/refs/vntr_db_advntr/hg19_muc1.db",
        "advntr_reference_vntr_hg38": "/refs/vntr_db_advntr/hg38_muc1.db",
    }
    resolved = resolve_from_mapping("advntr", label, mapping)
    assert resolved.value.endswith(f"{EXPECTED[label]}_muc1.db")


@pytest.mark.parametrize("label", ["hg38_ncbi", "hg38_ensembl"])
def test_the_two_labels_that_used_to_load_the_wrong_database(label):
    mapping = {"advntr_reference_vntr_hg19": "/hg19.db", "advntr_reference_vntr_hg38": "/hg38.db"}
    assert resolve_from_mapping("advntr", label, mapping).value == "/hg38.db"


def test_no_label_is_left_unmapped():
    assert set(list_assemblies()) == set(EXPECTED)


@pytest.mark.parametrize("label", sorted(EXPECTED))
def test_the_pipeline_itself_selects_the_right_database(label):
    from vntyper.scripts.pipeline import select_advntr_reference

    config = {
        "reference_data": {
            "advntr_reference_vntr_hg19": "/refs/hg19_muc1.db",
            "advntr_reference_vntr_hg38": "/refs/hg38_muc1.db",
        }
    }
    assert select_advntr_reference(config, label).endswith(f"{EXPECTED[label]}_muc1.db")


class TestMissingOrDisabledReferenceData:
    """`resolve_from_mapping`'s null/absent behaviour was, until now, exercised only for
    the `bwa` kind (`test_reference_resolution.py::TestResolveFromMapping`).
    `select_advntr_reference` is the reader whose original bug - the four-entry dict with
    `.get(label, "hg19")` this module's docstring describes - motivated the milestone, so
    it needs its own direct coverage of the same two membership-vs-truthiness cases
    rather than relying on the generic `resolve_from_mapping` tests to stand in for it.
    """

    @pytest.mark.parametrize("label", sorted(EXPECTED))
    def test_a_missing_reference_data_section_resolves_to_none(self, label):
        from vntyper.scripts.pipeline import select_advntr_reference

        assert select_advntr_reference({}, label) is None

    @pytest.mark.parametrize("label", sorted(EXPECTED))
    def test_a_present_but_null_database_key_resolves_to_none_not_the_string_null(self, label):
        from vntyper.scripts.pipeline import select_advntr_reference
        from vntyper.scripts.reference_registry import reference_keys

        (key,) = reference_keys("advntr", label)
        config = {"reference_data": {key: None}}

        assert select_advntr_reference(config, label) is None
