"""#223, second half: the derived VCFs must not fail open into a Negative.

``process_kestrel_output`` splits Kestrel's VCF into ``output_insertion.vcf`` and
``output_deletion.vcf`` and parses both. ``read_vcf_without_comments`` converts a missing or
unreadable path into an empty frame, and two empty frames are rendered as the ``Negative``
placeholder -- so a derived file that fails independently of the raw VCF is a second, distinct
route to a manufactured negative. ``run_kestrel``'s check on the raw file cannot reach it.

The issue states this as its second requirement: "Kestrel ran and found nothing" and "Kestrel
produced nothing readable" must not render identically.

``read_vcf_without_comments`` is deliberately not changed. Its empty-frame fallback is a
reviewed disposition in ``scripts/ble001_policy.json`` with tests asserting it as correct, so
the check lives at the call site.
"""

import pytest

from vntyper.scripts import kestrel_genotyping as kg

pytestmark = pytest.mark.unit

META = "##fileformat=VCFv4.2\n"
HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _raw_vcf(tmp_path, body=""):
    """Write a well-formed raw Kestrel VCF, so only the derived files are under test."""
    path = tmp_path / "output.vcf"
    path.write_text(META + HEADER + body, encoding="utf-8")
    return path


def _call(tmp_path, vcf):
    return kg.process_kestrel_output(str(tmp_path), vcf, "ref.fa", {}, {})


def test_a_header_only_vcf_still_produces_the_legitimate_negative(tmp_path, monkeypatch):
    """The guard must not turn a true negative into an aborted run.

    A Kestrel VCF with a valid header and no indel records is a real empty result. It goes
    through the *real* filters and the *real* derived guard here -- nothing is mocked -- and
    must still reach ``output_empty_result``.
    """
    monkeypatch.setattr(kg, "generate_header", lambda reference_vntr: ["##fileformat=VCFv4.2"])

    assert _call(tmp_path, _raw_vcf(tmp_path)) is None

    result = tmp_path / "kestrel_result.tsv"
    assert result.is_file(), "the legitimate empty result must still be written"
    assert "Negative" in result.read_text(encoding="utf-8")


def test_a_missing_derived_vcf_raises_rather_than_reporting_a_negative(tmp_path, monkeypatch):
    """If the split writes nothing, both frames come back empty and look like "no variants"."""
    monkeypatch.setattr(kg, "filter_indel_vcf", lambda indel_vcf, output_ins, output_del: None)

    with pytest.raises(RuntimeError) as excinfo:
        _call(tmp_path, _raw_vcf(tmp_path))

    assert "could not be read" in str(excinfo.value)
    assert not (tmp_path / "kestrel_result.tsv").exists(), "no Negative may be written"


def test_a_headerless_derived_vcf_raises_rather_than_reporting_a_negative(tmp_path, monkeypatch):
    """The insertion file exists and carries records, but lost its header on the way."""

    def bad_split(indel_vcf, output_ins, output_del):
        with open(output_ins, "w") as handle:
            handle.write(META + "chr1\t1\t.\tC\tCG\t.\t.\tDP=1\n")
        with open(output_del, "w") as handle:
            handle.write(META + HEADER)

    monkeypatch.setattr(kg, "filter_indel_vcf", bad_split)

    with pytest.raises(RuntimeError) as excinfo:
        _call(tmp_path, _raw_vcf(tmp_path))

    assert "output_insertion.vcf" in str(excinfo.value)
    assert "no #CHROM header line" in str(excinfo.value)
    assert not (tmp_path / "kestrel_result.tsv").exists()


def test_the_deletion_file_is_checked_too(tmp_path, monkeypatch):
    """Both derived files are parsed, so both are guarded. Checking only the first would leave
    half the path open."""

    def bad_split(indel_vcf, output_ins, output_del):
        with open(output_ins, "w") as handle:
            handle.write(META + HEADER)
        with open(output_del, "w") as handle:
            handle.write(META)

    monkeypatch.setattr(kg, "filter_indel_vcf", bad_split)

    with pytest.raises(RuntimeError) as excinfo:
        _call(tmp_path, _raw_vcf(tmp_path))

    assert "output_deletion.vcf" in str(excinfo.value)
    assert not (tmp_path / "kestrel_result.tsv").exists()
