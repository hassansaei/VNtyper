"""The VCF column contract both preprocessors depend on.

12 of ``motif_processing.py``'s 46 surviving mutants lived in
``preprocessing_insertion`` and ``preprocessing_deletion``, and every one was
column plumbing that no test observed: the ``#CHROM`` -> ``Motifs`` rename, the
five-element drop list, and the positional "the last column is the sample"
assumption. The existing tests all assert on the *merged* output, so none of
them could see a column being renamed wrongly.

Each test below kills at least one of those mutants. They are
CHARACTERISATION tests: they record the column contract the pipeline already
relies on (``kestrel_genotyping.py`` reads ``Sample``, ``Motifs`` and
``Variant`` straight out of these frames), they do not claim it is the
contract anyone would design today.
"""

import pandas as pd
import pytest

from vntyper.scripts.motif_processing import (
    VCF_METADATA_COLUMNS,
    _preprocess_vcf_frame,
    preprocessing_deletion,
    preprocessing_insertion,
)

pytestmark = pytest.mark.unit


def _vcf_frame():
    """A frame shaped like a filtered Kestrel VCF, with the sample column last."""
    return pd.DataFrame(
        {
            "#CHROM": ["X-5"],
            "POS": [67],
            "ID": ["."],
            "REF": ["C"],
            "ALT": ["CC"],
            "QUAL": ["."],
            "FILTER": ["PASS"],
            "INFO": ["."],
            "FORMAT": ["GT:AD"],
            "SAMPLE_1": ["Ins:50:5000"],
        }
    )


def _muc1_ref():
    return pd.DataFrame({"Motifs": ["X-5"], "Motif_sequence": ["ACGTACGT"]})


@pytest.mark.parametrize(
    ("func", "expected_label"),
    [(preprocessing_insertion, "Insertion"), (preprocessing_deletion, "Deletion")],
)
def test_chrom_is_renamed_to_motifs(func, expected_label):
    """The '#CHROM' -> 'Motifs' rename is the merge key; without it the merge cannot run."""
    out = func(_vcf_frame(), _muc1_ref())
    assert "Motifs" in out.columns
    assert "#CHROM" not in out.columns
    assert out["Variant"].iloc[0] == expected_label


@pytest.mark.parametrize("func", [preprocessing_insertion, preprocessing_deletion])
def test_exactly_the_five_vcf_metadata_columns_are_dropped(func):
    """The drop list is exactly the five VCF metadata columns - no more, no fewer."""
    out = func(_vcf_frame(), _muc1_ref())
    for dropped in ("ID", "QUAL", "FILTER", "INFO", "FORMAT"):
        assert dropped not in out.columns
    for kept in ("Motifs", "POS", "REF", "ALT"):
        assert kept in out.columns


def test_the_drop_list_is_shared_by_both_preprocessors():
    """One list, one place to change it: the twins drifting apart was the original risk."""
    assert VCF_METADATA_COLUMNS == ["ID", "QUAL", "FILTER", "INFO", "FORMAT"]


@pytest.mark.parametrize("func", [preprocessing_insertion, preprocessing_deletion])
def test_the_last_column_becomes_sample_and_it_is_the_last_one(func):
    """``df.columns[-1]`` is a positional assumption nothing downstream re-checks.

    The sample column is whatever ends up last after the drop. A ``-1`` -> ``+1``
    or ``-1`` -> ``-2`` mutation renames a *different* column to 'Sample', and the
    genotype/depth string the whole caller chain parses is silently replaced by a
    position or an allele.
    """
    out = func(_vcf_frame(), _muc1_ref())
    assert "Sample" in out.columns
    assert out["Sample"].iloc[0] == "Ins:50:5000"
    assert "SAMPLE_1" not in out.columns


@pytest.mark.parametrize("func", [preprocessing_insertion, preprocessing_deletion])
def test_a_second_sample_column_means_only_the_last_is_renamed(func):
    """A multi-sample VCF: the positional rule must take the final column, not any sample."""
    frame = _vcf_frame()
    frame.insert(len(frame.columns) - 1, "SAMPLE_0", ["Ins:1:2"])
    out = func(frame, _muc1_ref())
    assert out["Sample"].iloc[0] == "Ins:50:5000"
    assert "SAMPLE_0" in out.columns


@pytest.mark.parametrize("func", [preprocessing_insertion, preprocessing_deletion])
def test_the_merge_key_is_motifs_and_an_unmatched_motif_yields_na(func):
    """A left merge: an unknown motif ID produces NA rather than dropping the variant."""
    out = func(_vcf_frame(), pd.DataFrame({"Motifs": ["NOPE"], "Motif_sequence": ["TTTT"]}))
    assert out["Motif_sequence"].isna().all()
    assert len(out) == 1, "a left merge keeps the row even when the motif is unknown"


@pytest.mark.parametrize("variant_label", ["Insertion", "Deletion"])
def test_the_shared_body_differs_only_in_the_variant_label(variant_label):
    """The two public functions are ``_preprocess_vcf_frame`` with one literal changed."""
    out = _preprocess_vcf_frame(_vcf_frame(), _muc1_ref(), variant_label)
    assert out["Variant"].tolist() == [variant_label]


def test_both_public_functions_produce_identical_frames_apart_from_variant():
    """The merge cannot fix on one path and stay broken on the other."""
    insertion = preprocessing_insertion(_vcf_frame(), _muc1_ref())
    deletion = preprocessing_deletion(_vcf_frame(), _muc1_ref())

    assert insertion.columns.tolist() == deletion.columns.tolist()
    pd.testing.assert_frame_equal(insertion.drop(columns=["Variant"]), deletion.drop(columns=["Variant"]))
