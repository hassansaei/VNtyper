"""The positive-result publication tail of ``process_kestrel_output``.

Ten unit tests monkeypatch ``process_kestrel_output`` whole, ``run_kestrel`` discards
its return value, and the one test that reached it patched
``process_kmer_results`` to return an empty frame -- so the lines that write a
*positive* ``kestrel_result.tsv`` never executed under test, and a mutant that wrote
``processed_df.iloc[0:0]`` to the file survived the entire unit tier. This file
drives the real ``process_kestrel_output`` through the non-empty path and asserts on
values parsed back from the written file with the repo's own reader.
"""

from pathlib import Path
from unittest import mock

import pandas as pd
import pytest

from tests.builders import kestrel_stage_frame
from vntyper.scripts import kestrel_genotyping as kg
from vntyper.scripts.identity_candidates import IdentityTranslationComponent
from vntyper.scripts.summary import parse_tsv

pytestmark = pytest.mark.unit

META = "##fileformat=VCFv4.2\n"
HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
RECORD = "chr1\t100\t.\tC\tCC\t.\t.\tDP=10\tGT:GDP:DP\t1:120:12000\n"

EXPECTED_SCHEMA = [
    "Motifs",
    "Variant",
    "POS",
    "REF",
    "ALT",
    "Sample",
    "Motif_sequence",
    "Del",
    "Estimated_Depth_AlternateVariant",
    "Estimated_Depth_Variant_ActiveRegion",
    "ref_len",
    "alt_len",
    "Frame_Score",
    "is_frameshift",
    "direction",
    "frameshift_amount",
    "is_valid_frameshift",
    "Depth_Score",
    "Confidence",
    "depth_confidence_pass",
    "haplo_count",
    "alt_filter_pass",
    "motif_filter_pass",
    "Motif",
    "Flag",
    "Motif_fasta",
    "POS_fasta",
    "flag_filter_pass",
    "Nomenclature",
    "Nomenclature_Tier",
    "Nomenclature_Flags",
    "Ambiguity_Interval",
    "Repeat_Form",
    "Nomenclature_Note",
    "Nomenclature_Kestrel",
    "Nomenclature_adVNTR",
]

EXPECTED_ROW = {
    "Motifs": "X-5",
    "Variant": "Insertion",
    "POS": "67",
    "REF": "G",
    "ALT": "GG",
    "Sample": "Del:120:12000",
    "Motif_sequence": "SEQ1",
    "Del": "Del",
    "Estimated_Depth_AlternateVariant": "120",
    "Estimated_Depth_Variant_ActiveRegion": "12000",
    "ref_len": "1",
    "alt_len": "2",
    "Frame_Score": "0.3333333333333333",
    "is_frameshift": "True",
    "direction": "1",
    "frameshift_amount": "1",
    "is_valid_frameshift": "True",
    "Depth_Score": "0.01",
    "Confidence": "High_Precision*",
    "depth_confidence_pass": "True",
    "haplo_count": "1",
    "alt_filter_pass": "True",
    "motif_filter_pass": "True",
    "Motif": "X",
    "Flag": "Not flagged",
    "Motif_fasta": "X",
    "POS_fasta": "67",
    "flag_filter_pass": "True",
    "Nomenclature": "59dupC",
    "Nomenclature_Tier": "B",
    "Nomenclature_Flags": "known-variant;position-ambiguous",
    "Ambiguity_Interval": "53_59",
    "Repeat_Form": "53C[7]>53C[8]",
    "Nomenclature_Note": "matches a described MUC1 variant (Kirby et al. 2013, PMID:23396133); requires validation",
    "Nomenclature_Kestrel": "59dupC",
    "Nomenclature_adVNTR": "",
}


def test_a_positive_result_publishes_its_data_row_to_the_tsv(tmp_path: Path) -> None:
    """The written ``kestrel_result.tsv`` must carry the complete processed row.

    The raw VCF, the indel filters, the derived-VCF guard, the preprocessing merge,
    the real header, the real nomenclature annotation and the final write all run
    unmocked. Only the reference-file loaders, the bcftools compression convenience
    and the scoring chain are replaced; the scored frame handed back is the
    builders' canonical positive call. The file is then parsed back with
    ``summary.parse_tsv`` so the assertions are on what a consumer reads, not on
    what the writer intended to write.
    """
    vcf = tmp_path / "output.vcf"
    vcf.write_text(META + HEADER + RECORD, encoding="utf-8")
    seen: dict[str, object] = {}
    runtime_identity_config = {
        "motifs": {"runtime-X": "A"},
        "advntr": {"mappable_repeat_units": {"runtime-2": "runtime-X"}, "rotation_offset": 17},
    }

    def fake_process_kmer_results(
        combined_df,
        merged_motifs,
        output_dir,
        config,
        compiled_flag_rules=None,
        identity_component=None,
    ):
        seen["combined_rows"] = len(combined_df)
        seen["alts"] = list(combined_df["ALT"])
        seen["identity_component"] = identity_component
        return kestrel_stage_frame("final")

    with (
        mock.patch.object(kg, "_try_compress_vcf_with_bcftools", lambda *args: None),
        mock.patch.object(
            kg,
            "load_muc1_reference",
            return_value=pd.DataFrame({"Motifs": ["chr1"], "Motif_sequence": ["ACGT"]}),
        ),
        mock.patch.object(kg, "load_additional_motifs", return_value=pd.DataFrame()),
        mock.patch.object(kg, "process_kmer_results", side_effect=fake_process_kmer_results),
        mock.patch.object(
            kg,
            "load_nomenclature_config",
            return_value=runtime_identity_config,
            create=True,
        ) as load_identity_config,
    ):
        returned = kg.process_kestrel_output(str(tmp_path), vcf, "ref.fa", {}, {})

    assert seen["combined_rows"] == 1 and seen["alts"] == ["CC"], "the real VCF record must reach the scoring seam"
    assert isinstance(seen["identity_component"], IdentityTranslationComponent)
    load_identity_config.assert_called_once_with()
    assert seen["identity_component"].kestrel_motifs == {"runtime-X": "A"}
    assert seen["identity_component"].advntr_rotation_offset == 17
    assert returned is not None and len(returned) == 1

    parsed = parse_tsv(str(tmp_path / "kestrel_result.tsv"))
    assert parsed["comments"][0] == "VNtyper Kestrel result"
    rows = parsed["data"]
    assert len(rows) == 1, "the published TSV must carry the data row, not just headers"
    row = rows[0]
    assert list(row) == EXPECTED_SCHEMA
    assert row == EXPECTED_ROW
