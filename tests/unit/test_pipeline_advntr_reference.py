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
