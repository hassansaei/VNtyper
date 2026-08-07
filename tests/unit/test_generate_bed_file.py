"""Unit tests for ``generate_bed_file`` (#203).

This function had no test of its output content at all, which is how a coordinate
convention error survived: nothing asserted what it writes, only that it returned a path.
"""

import pandas as pd
import pytest

from vntyper.scripts.kestrel_genotyping import generate_bed_file

pytestmark = pytest.mark.unit


def test_it_writes_zero_based_half_open_intervals(tmp_path):
    """VCF POS is 1-based; BED is 0-based half-open, so POS p is the interval [p-1, p).

    The previous output was [p, p+1), naming the base *after* the variant - so IGV
    highlighted the wrong position for every row.
    """
    df = pd.DataFrame({"Motif_fasta": ["X-5"], "POS_fasta": [67]})

    generate_bed_file(df, str(tmp_path))

    assert (tmp_path / "output.bed").read_text(encoding="utf-8") == "X-5\t66\t67\n"


def test_it_names_the_pair_record_that_the_reference_fasta_actually_contains(tmp_path):
    """Column 1 must resolve against the FASTA the IGV report is given.

    ``Motif_fasta`` is the 120 bp pair name (``1-2``); ``Motif`` is the half name
    (``2``), which is not a record in that file. Writing the half name would produce a
    BED nothing can resolve.
    """
    df = pd.DataFrame({"Motif_fasta": ["1-2"], "POS_fasta": [54], "Motif": ["2"]})

    generate_bed_file(df, str(tmp_path))

    assert (tmp_path / "output.bed").read_text(encoding="utf-8").split("\t")[0] == "1-2"


def test_it_writes_one_interval_per_row(tmp_path):
    df = pd.DataFrame({"Motif_fasta": ["1-2", "X-5"], "POS_fasta": [54, 67]})

    generate_bed_file(df, str(tmp_path))

    assert (tmp_path / "output.bed").read_text(encoding="utf-8") == "1-2\t53\t54\nX-5\t66\t67\n"


@pytest.mark.parametrize("columns", [{"Motif_fasta": ["1-2"]}, {"POS_fasta": [54]}])
def test_it_returns_none_when_either_column_is_absent(columns, tmp_path):
    assert generate_bed_file(pd.DataFrame(columns), str(tmp_path)) is None
    assert not (tmp_path / "output.bed").exists()


def test_it_returns_none_for_an_empty_frame(tmp_path):
    df = pd.DataFrame({"Motif_fasta": [], "POS_fasta": []})

    assert generate_bed_file(df, str(tmp_path)) is None
