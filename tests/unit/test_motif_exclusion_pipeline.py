"""Production-path coverage for the conserved-motif exclusion."""

import pandas as pd
import pytest

from tests.builders import kestrel_config, kestrel_stage_frame
from vntyper.scripts import kestrel_genotyping

pytestmark = pytest.mark.unit


def test_the_measured_low_depth_artifact_remains_final_empty_by_independent_gates(tmp_path):
    """C1 adds an earlier exclusion without weakening the two existing defences.

    The three measured ``8-7 / POS 70 / C>CGGCA`` rows are already Negative:
    both the depth gate and the config-driven artifact gate independently reject
    them. The motif exclusion is defence in depth, not the sole reason no final
    call is emitted.
    """
    frame = kestrel_stage_frame(
        "raw",
        rows=1,
        motifs="8-7",
        motif="8",
        pos=70,
        ref="C",
        alt="CGGCA",
        depth_alt=1,
        depth_region=1800,
    )
    merged = pd.DataFrame({"Motif": ["8", "7"], "Motif_sequence": ["SEQ8", "SEQ7"]})

    out = kestrel_genotyping.process_kmer_results(frame, merged, str(tmp_path), kestrel_config())

    pre = pd.read_csv(tmp_path / "kestrel_pre_result.tsv", sep="\t")
    assert out.empty
    assert pre["Confidence"].tolist() == ["Negative"]
    assert pre["depth_confidence_pass"].tolist() == [False]
    assert pre["flag_filter_pass"].tolist() == [False]
    assert pre["motif_filter_pass"].tolist() == [False]


def test_an_otherwise_passing_conserved_right_motif_is_excluded_without_a_gg_sibling(tmp_path):
    """The C1 semantic change is the motif gate itself, not the measured artifact's other gates.

    No such row became the selected call in the measured 400-sample cohort, but a future
    otherwise-passing right-half motif 8 call must not regain the legacy dependency on
    an unrelated ``ALT=GG`` sibling.
    """
    frame = kestrel_stage_frame(
        "raw",
        rows=1,
        motifs="8-7",
        motif="8",
        pos=70,
        ref="G",
        alt="GC",
        depth_alt=51,
        depth_region=5000,
    )
    merged = pd.DataFrame({"Motif": ["8", "7"], "Motif_sequence": ["SEQ8", "SEQ7"]})

    out = kestrel_genotyping.process_kmer_results(frame, merged, str(tmp_path), kestrel_config())

    pre = pd.read_csv(tmp_path / "kestrel_pre_result.tsv", sep="\t")
    assert out.empty
    assert pre["Confidence"].tolist() == ["High_Precision"]
    assert pre["depth_confidence_pass"].tolist() == [True]
    assert pre["flag_filter_pass"].tolist() == [True]
    assert pre["motif_filter_pass"].tolist() == [False]
