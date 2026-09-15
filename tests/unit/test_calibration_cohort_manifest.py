"""Development cohort input validation keeps optional truth and identities separate."""

from pathlib import Path

import pytest

pytestmark = pytest.mark.unit


def manifest(tmp_path, rows, header="sample_id\tbam\tassembly\tgenotype\tallele_1\tallele_2\tgroup_id"):
    path = tmp_path / "samples.tsv"
    path.write_text(header + "\n" + "\n".join(rows) + "\n")
    return path


def test_optional_truth_and_explicit_groups_are_preserved(tmp_path):
    from vntyper.scripts.calibration_cohort_manifest import read_cohort_manifest

    path = manifest(
        tmp_path,
        [
            "a\ta.bam\tGRCh38\tpositive\t10\t20\tfamily",
            "b\tb.bam\tGRCh38\tnegative\t\t\tfamily",
            "c\tc.bam\thg19\t\t\t\t",
        ],
    )
    rows = read_cohort_manifest(path)
    assert rows[0].length_total == 30
    assert rows[1].length_total is None
    assert rows[2].genotype is None
    assert rows[0].group_id == rows[1].group_id
    assert rows[2].group_id != rows[0].group_id
    assert rows[0].bam == tmp_path / "a.bam"


@pytest.mark.parametrize(
    "rows",
    [
        ["a\ta.bam\tGRCh38\tpositive\tnan\t20\t"],
        ["a\ta.bam\tGRCh38\tPOS\t10\t20\t"],
        ["a\ta.bam\tGRCh38\tnegative\t10\t20\t", "a\tb.bam\tGRCh38\tpositive\t10\t20\t"],
        ["a\ta.bam\tGRCh38\tnegative\t10\t20\t", "b\ta.bam\tGRCh38\tnegative\t10\t20\t"],
    ],
)
def test_invalid_or_duplicate_truth_is_refused(tmp_path, rows):
    from vntyper.scripts.calibration_cohort_manifest import read_cohort_manifest

    with pytest.raises(ValueError):
        read_cohort_manifest(manifest(tmp_path, rows))


def test_absent_truth_variant_column_leaves_every_identity_unavailable(tmp_path):
    from vntyper.scripts.calibration_cohort_manifest import read_cohort_manifest

    rows = read_cohort_manifest(manifest(tmp_path, ["a\ta.bam\tGRCh38\tpositive\t\t\t"]))
    assert rows[0].truth_variant is None


def test_declared_truth_variant_is_preserved_and_an_empty_cell_stays_unavailable(tmp_path):
    from vntyper.scripts.calibration_cohort_manifest import read_cohort_manifest

    header = "sample_id\tbam\tassembly\tgenotype\ttruth_variant"
    rows = read_cohort_manifest(
        manifest(
            tmp_path,
            ["a\ta.bam\tGRCh38\tpositive\tMUC1-X-60-coding-v1|60|59|-|C", "b\tb.bam\tGRCh38\tpositive\t"],
            header,
        )
    )
    assert rows[0].truth_variant == "MUC1-X-60-coding-v1|60|59|-|C"
    assert rows[1].truth_variant is None


@pytest.mark.parametrize("genotype", ["negative", "unknown", ""])
def test_a_confirmed_variant_without_a_positive_genotype_is_refused(tmp_path, genotype):
    from vntyper.scripts.calibration_cohort_manifest import read_cohort_manifest

    header = "sample_id\tbam\tassembly\tgenotype\ttruth_variant"
    with pytest.raises(ValueError, match="truth_variant"):
        read_cohort_manifest(
            manifest(tmp_path, [f"a\ta.bam\tGRCh38\t{genotype}\tMUC1-X-60-coding-v1|60|59|-|C"], header)
        )


def test_minimal_manifest_and_unknown_column_refusal(tmp_path):
    from vntyper.scripts.calibration_cohort_manifest import read_cohort_manifest

    rows = read_cohort_manifest(manifest(tmp_path, ["x\tx.bam\tGRCh38"], "sample_id\tbam\tassembly"))
    assert rows[0].length_total is None
    with pytest.raises(ValueError, match="columns"):
        read_cohort_manifest(manifest(tmp_path, ["x\tx.bam\tGRCh38\tlabel"], "sample_id\tbam\tassembly\tpredictor"))


def test_partial_length_keeps_mutation_truth_and_explicit_reason(tmp_path):
    from vntyper.scripts.calibration_cohort_manifest import read_cohort_manifest

    row = read_cohort_manifest(manifest(tmp_path, ["a\ta.bam\tGRCh38\tpositive\t10\t\t"]))[0]
    assert row.genotype is True
    assert row.length_total is None
    assert row.length_reason == "incomplete_allele_pair"


def test_byte_copied_alignments_cannot_contribute_twice(tmp_path):
    from vntyper.scripts.calibration_cohort_manifest import audit_alignment_duplicates, read_cohort_manifest

    (tmp_path / "a.bam").write_bytes(b"invented same read artifact")
    (tmp_path / "b.bam").write_bytes(b"invented same read artifact")
    rows = read_cohort_manifest(
        manifest(tmp_path, ["a\ta.bam\tGRCh38\tpositive\t\t\t", "b\tb.bam\tGRCh38\tnegative\t\t\t"])
    )
    with pytest.raises(ValueError, match="duplicate alignment content"):
        audit_alignment_duplicates(rows)


def test_content_audit_ignores_access_time_changes(tmp_path, monkeypatch):
    import os

    from vntyper.scripts.calibration_cohort_manifest import audit_alignment_duplicates, read_cohort_manifest

    path = tmp_path / "a.bam"
    path.write_bytes(b"invented alignment")
    rows = read_cohort_manifest(manifest(tmp_path, ["a\ta.bam\tGRCh38\tpositive\t\t\t"]))
    real_stat = Path.stat
    count = 0

    def changed_atime(self, *args, **kwargs):
        nonlocal count
        result = real_stat(self, *args, **kwargs)
        if self == path:
            count += 1
            values = list(result)
            values[7] += count
            return os.stat_result(values)
        return result

    monkeypatch.setattr(Path, "stat", changed_atime)
    assert audit_alignment_duplicates(rows)["a"] is not None


@pytest.mark.parametrize(
    "text",
    [
        "sample_id\tbam\tassembly\n",
        "sample_id\tbam\tassembly\na\ta.bam\tGRCh38\textra\n",
        "sample_id\tbam\tassembly\na\ta.bam\n",
        "sample_id\tbam\tassembly\na \ta.bam\tGRCh38\n",
        "sample_id\tbam\tassembly\n\ta.bam\tGRCh38\n",
        "sample_id\tbam\tassembly\tbam\na\ta.bam\tGRCh38\tb.bam\n",
    ],
)
def test_malformed_tsv_rows_and_headers_are_refused(tmp_path, text):
    from vntyper.scripts.calibration_cohort_manifest import read_cohort_manifest

    path = tmp_path / "bad.tsv"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_cohort_manifest(path)


def test_hardlinked_alignment_duplicate_is_refused_before_hashing(tmp_path):
    import os

    from vntyper.scripts.calibration_cohort_manifest import read_cohort_manifest

    (tmp_path / "a.bam").write_bytes(b"invented")
    os.link(tmp_path / "a.bam", tmp_path / "b.bam")
    with pytest.raises(ValueError, match="inode"):
        read_cohort_manifest(
            manifest(tmp_path, ["a\ta.bam\tGRCh38\tpositive\t\t\t", "b\tb.bam\tGRCh38\tpositive\t\t\t"])
        )
