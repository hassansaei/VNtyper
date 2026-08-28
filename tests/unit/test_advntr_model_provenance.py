"""The adVNTR model a run used must be recorded, and mismatches must fail loudly.

adVNTR derives its read-fetch window from the model's own content, so which model file
a run resolved decides which reads it could ever see. Before #268 that was untraceable:
`pipeline_summary.json` recorded neither the file nor its content.

Two failure modes are covered here.

Silent truncation: a v2 model records the array's genomic end, and an adVNTR predating
that field ignores it and reads the old 840 bp window. That produces wrong answers with
no error, so it must be refused before the run rather than diagnosed afterwards.

Silent staleness: an installation still carrying only the v1 database must be told to
reinstall references, never quietly genotyped against the defective model.

These tests use temporary SQLite fixtures rather than the installed bundle, so the unit
tier keeps running on a fresh clone.
"""

import sqlite3

import pytest

from vntyper.modules.advntr.model_provenance import (
    AdvntrModelError,
    describe_model,
    require_compatible_advntr,
)

pytestmark = pytest.mark.unit


LEGACY_SQL = (
    "CREATE TABLE vntrs(id INTEGER PRIMARY KEY, nonoverlapping TEXT, chromosome TEXT, "
    "ref_start INTEGER, gene_name TEXT, annotation TEXT, pattern TEXT, "
    "left_flanking TEXT, right_flanking TEXT, repeats TEXT, scaled_score REAL DEFAULT 0)"
)
V2_SQL = (
    "CREATE TABLE vntrs_v2(id INTEGER PRIMARY KEY, nonoverlapping TEXT, chromosome TEXT, "
    "ref_start INTEGER, gene_name TEXT, annotation TEXT, pattern TEXT, "
    "left_flanking TEXT, right_flanking TEXT, repeats TEXT, scaled_score REAL DEFAULT 0, "
    "ref_end INTEGER)"
)

SEGMENTS = ["A" * 60, "C" * 60, "A" * 60, "G" * 48]


def _write(path, sql, ref_end=None, segments=SEGMENTS, ref_start=155188507):
    table = sql.split("CREATE TABLE ")[1].split("(")[0]
    values = [
        25561,
        "True",
        "chr1",
        ref_start,
        "MUC1",
        "Coding",
        "A" * 60,
        "L" * 500,
        "R" * 500,
        ",".join(segments),
        0.0,
    ]
    if ref_end is not None:
        values.append(ref_end)
    db = sqlite3.connect(str(path))
    db.execute(sql)
    db.execute(f"INSERT INTO {table} VALUES ({','.join(['?'] * len(values))})", values)
    db.commit()
    db.close()
    return path


class TestDescribeModel:
    def test_v2_model_reports_the_recorded_genomic_interval(self, tmp_path):
        path = _write(tmp_path / "v2.db", V2_SQL, ref_end=155192032)
        info = describe_model(path)
        assert info["schema_version"] == "v2"
        assert info["genomic_interval"] == "chr1:155188507-155192032"
        assert info["window_bp"] == 3525

    def test_legacy_model_reports_the_derived_interval(self, tmp_path):
        path = _write(tmp_path / "v1.db", LEGACY_SQL)
        info = describe_model(path)
        assert info["schema_version"] == "v1"
        # 228 bp of stored units, used as a coordinate -- the defect, made visible.
        assert info["window_bp"] == 228

    def test_segment_shape_is_recorded_not_summarised_as_a_span(self, tmp_path):
        path = _write(tmp_path / "v2.db", V2_SQL, ref_end=155192032)
        info = describe_model(path)
        assert info["n_segments"] == 4
        assert info["n_distinct_segments"] == 3
        assert info["max_segment_len"] == 60
        # A sum of segment lengths must never be presented as a genomic span: that is
        # the conflation this whole change exists to remove.
        assert "span_bp" not in info

    def test_digest_identifies_the_file(self, tmp_path):
        a = _write(tmp_path / "a.db", V2_SQL, ref_end=155192032)
        b = _write(tmp_path / "b.db", V2_SQL, ref_end=155192032, segments=SEGMENTS[:2])
        assert describe_model(a)["sha256"] != describe_model(b)["sha256"]
        assert len(describe_model(a)["sha256"]) == 64

    def test_a_file_with_no_vntr_table_is_rejected(self, tmp_path):
        path = tmp_path / "empty.db"
        sqlite3.connect(str(path)).execute("CREATE TABLE other(x INTEGER)")
        with pytest.raises(AdvntrModelError):
            describe_model(path)

    def test_a_missing_file_is_rejected(self, tmp_path):
        with pytest.raises(AdvntrModelError):
            describe_model(tmp_path / "absent.db")

    def test_an_empty_table_is_rejected(self, tmp_path):
        path = tmp_path / "norows.db"
        sqlite3.connect(str(path)).execute(V2_SQL)
        with pytest.raises(AdvntrModelError) as exc:
            describe_model(path)
        assert "install-references" in str(exc.value)

    def test_a_multi_row_database_selects_the_muc1_model(self, tmp_path):
        # A general adVNTR database carries thousands of VNTRs. Selecting MUC1 beats
        # assuming the file is the single-row bundle.
        path = _write(tmp_path / "many.db", V2_SQL, ref_end=155192032)
        db = sqlite3.connect(str(path))
        db.execute(
            "INSERT INTO vntrs_v2 VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            (999, "True", "chr7", 1000, "OTHER", "Coding", "T" * 60, "L" * 500, "R" * 500, "T" * 60, 0.0, 1600),
        )
        db.commit()
        db.close()
        assert describe_model(path)["vid"] == 25561

    def test_a_multi_row_database_without_muc1_is_rejected(self, tmp_path):
        path = _write(tmp_path / "nomuc1.db", V2_SQL, ref_end=155192032)
        db = sqlite3.connect(str(path))
        db.execute("UPDATE vntrs_v2 SET id = 999")
        db.execute(
            "INSERT INTO vntrs_v2 VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            (998, "True", "chr7", 1000, "OTHER", "Coding", "T" * 60, "L" * 500, "R" * 500, "T" * 60, 0.0, 1600),
        )
        db.commit()
        db.close()
        with pytest.raises(AdvntrModelError) as exc:
            describe_model(path)
        assert "25561" in str(exc.value)

    @pytest.mark.parametrize("ref_end", [155188507, 155188000])
    def test_an_end_at_or_before_the_start_is_rejected(self, tmp_path, ref_end):
        path = _write(tmp_path / "backwards.db", V2_SQL, ref_end=ref_end)
        with pytest.raises(AdvntrModelError) as exc:
            describe_model(path)
        assert "ref_start" in str(exc.value)


class TestCompatibility:
    def test_v2_model_with_span_aware_advntr_is_accepted(self, tmp_path):
        path = _write(tmp_path / "v2.db", V2_SQL, ref_end=155192032)
        require_compatible_advntr(describe_model(path), (2, 0, 4))

    def test_v2_model_with_old_advntr_is_refused(self, tmp_path):
        # The dangerous case: 2.0.3 selects the legacy columns by name, so it would
        # ignore ref_end and silently genotype against the 840 bp window.
        path = _write(tmp_path / "v2.db", V2_SQL, ref_end=155192032)
        with pytest.raises(AdvntrModelError) as exc:
            require_compatible_advntr(describe_model(path), (2, 0, 3))
        assert "2.0.4" in str(exc.value)

    def test_unknown_advntr_version_with_a_v2_model_is_refused(self, tmp_path):
        path = _write(tmp_path / "v2.db", V2_SQL, ref_end=155192032)
        with pytest.raises(AdvntrModelError):
            require_compatible_advntr(describe_model(path), None)

    def test_legacy_model_names_the_defect_and_says_how_to_fix_it(self, tmp_path):
        path = _write(tmp_path / "v1.db", LEGACY_SQL)
        with pytest.raises(AdvntrModelError) as exc:
            require_compatible_advntr(describe_model(path), (2, 0, 4))
        message = str(exc.value)
        assert "install-references" in message
        assert "268" in message

    def test_a_legacy_model_is_never_silently_accepted(self, tmp_path):
        # Falling back to the defective model would reintroduce the 24% window with no
        # signal to the user, which is exactly what #268 is about.
        path = _write(tmp_path / "v1.db", LEGACY_SQL)
        for version in [(2, 0, 3), (2, 0, 4), None]:
            with pytest.raises(AdvntrModelError):
                require_compatible_advntr(describe_model(path), version)
