# tests/unit/test_advntr_output_parsing.py

"""
Unit tests for adVNTR output parsing (``vntyper/modules/advntr/advntr_genotyping.py``).

These pin two things the rest of the pipeline silently depends on:

1. **The 10-column result schema and its filename.** ``pipeline.py`` does not consume the
   path :func:`process_advntr_output` writes -- it *reconstructs* it from ``output_name``
   (``os.path.join(dirs["advntr"], "output_adVNTR_result.tsv")``). ``generate_report.py``
   and ``cohort_summary.py`` then read fixed column names out of that file. A renamed
   column or a renamed file does not fail; it silently empties the adVNTR half of the
   report.

2. **The "Negative" placeholder contract.** ``report_config.json`` classifies an adVNTR
   result by testing ``VID != "Negative"`` and ``Flag == "Not flagged"``. The placeholder
   row is therefore the only thing standing between "adVNTR found nothing" and a crash
   in the report layer -- it must always be written, with exactly these values.

They also cover the three paths on which :func:`process_advntr_output` returns *without*
writing a file. A caller that gets no file today does not see the real cause; it sees a
misleading downstream failure much later.
"""

import logging
from pathlib import Path

import pandas as pd
import pytest

from vntyper.modules.advntr import advntr_genotyping as advntr

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# The frozen output contract (see wave-2 contract C4 -- do not edit casually).
# ---------------------------------------------------------------------------

FINAL_COLUMNS = [
    "VID",
    "Variant",
    "NumberOfSupportingReads",
    "MeanCoverage",
    "Pvalue",
    "RU",
    "POS",
    "REF",
    "ALT",
    "Flag",
]

RESULT_SUFFIX = "_adVNTR_result.tsv"

#: adVNTR writes a ``#``-prefixed header; the parser normalises ``#VID`` to ``VID``.
ADVNTR_HEADER = "#VID\tState\tNumberOfSupportingReads\tMeanCoverage\tPvalue\n"

#: The variant the shipped integration fixture actually reports (a 1 bp insertion).
CANONICAL_VARIANT = "I22_2_G_LEN1"
CANONICAL_ROW = f"25561\t{CANONICAL_VARIANT}\t11\t153.98\t0.0001\n"


def write_advntr_output(tmp_path: Path, blob: str, name: str = "output_adVNTR.vcf") -> Path:
    """Write a literal adVNTR output blob and return its path."""
    path = tmp_path / name
    path.write_text(blob)
    return path


def read_result(tmp_path: Path, output_name: str = "output") -> pd.DataFrame:
    """Read the produced result TSV back as raw strings (this is a file-format contract)."""
    path = tmp_path / f"{output_name}{RESULT_SUFFIX}"
    assert path.exists(), f"expected {path} to be written, directory holds {sorted(p.name for p in tmp_path.iterdir())}"
    return pd.read_csv(path, sep="\t", dtype=str)


def assert_is_negative_placeholder(df: pd.DataFrame) -> None:
    """Assert the frame is exactly the one-row negative placeholder the report expects."""
    assert list(df.columns) == FINAL_COLUMNS
    assert len(df) == 1
    row = df.iloc[0]
    assert row["VID"] == "Negative"
    for column in FINAL_COLUMNS[1:]:
        assert row[column] == "Not applicable", f"{column} should be 'Not applicable', got {row[column]!r}"


@pytest.fixture
def ru_fasta(tmp_path: Path) -> Path:
    """A minimal repeat-unit FASTA covering the RUs used in these tests."""
    path = tmp_path / "code-adVNTR_RUs.fa"
    path.write_text("".join(f">RU{n}\n{'ACGTACGTAC' * 7}\n" for n in (2, 6, 7)))
    return path


@pytest.fixture
def ru_config(ru_fasta: Path) -> dict:
    """A main-config shape carrying the RU FASTA, as ``pipeline.py`` passes it."""
    return {"reference_data": {"code_adVNTR_RUs": str(ru_fasta)}}


# ---------------------------------------------------------------------------
# The happy path and the schema
# ---------------------------------------------------------------------------


class TestCanonicalParsing:
    """A well-formed adVNTR output produces the frozen 10-column result."""

    def test_a_canonical_variant_row_produces_the_ten_column_result(self, tmp_path):
        source = write_advntr_output(tmp_path, ADVNTR_HEADER + CANONICAL_ROW)

        advntr.process_advntr_output(str(source), str(tmp_path), "output")

        df = read_result(tmp_path)
        assert list(df.columns) == FINAL_COLUMNS
        assert len(df) == 1
        row = df.iloc[0]
        assert row["VID"] == "25561"
        assert row["Variant"] == CANONICAL_VARIANT
        assert row["NumberOfSupportingReads"] == "11"
        assert row["Flag"] == "Not flagged"

    def test_the_result_filename_is_built_from_output_name(self, tmp_path):
        """Contract C4: ``pipeline.py`` reconstructs this name rather than consuming it."""
        source = write_advntr_output(tmp_path, ADVNTR_HEADER + CANONICAL_ROW)

        advntr.process_advntr_output(str(source), str(tmp_path), "sample_42")

        assert (tmp_path / f"sample_42{RESULT_SUFFIX}").exists()
        assert not (tmp_path / f"output{RESULT_SUFFIX}").exists()

    def test_the_hash_prefixed_header_is_normalised_in_the_source_file(self, tmp_path):
        """``#VID`` is rewritten to ``VID`` in place, otherwise ``comment='#'`` eats the header."""
        source = write_advntr_output(tmp_path, ADVNTR_HEADER + CANONICAL_ROW)

        advntr.process_advntr_output(str(source), str(tmp_path), "output")

        assert source.read_text().startswith("VID\tState\t")

    def test_rows_identical_in_vid_variant_and_read_count_are_deduplicated(self, tmp_path):
        source = write_advntr_output(
            tmp_path,
            ADVNTR_HEADER + CANONICAL_ROW + f"25561\t{CANONICAL_VARIANT}\t11\t153.98\t0.0002\n",
        )

        advntr.process_advntr_output(str(source), str(tmp_path), "output")

        assert len(read_result(tmp_path)) == 1

    def test_ru_annotation_fills_the_ru_pos_ref_and_alt_columns(self, tmp_path, ru_config):
        source = write_advntr_output(tmp_path, ADVNTR_HEADER + CANONICAL_ROW)

        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=ru_config)

        row = read_result(tmp_path).iloc[0]
        assert row["RU"] == "2"
        assert row["POS"] == "22"
        assert row["REF"] != "Not applicable"
        assert row["ALT"].startswith(row["REF"])

    def test_without_an_ru_fasta_the_annotation_columns_fall_back_to_not_applicable(self, tmp_path):
        source = write_advntr_output(tmp_path, ADVNTR_HEADER + CANONICAL_ROW)

        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=None)

        row = read_result(tmp_path).iloc[0]
        for column in ("RU", "POS", "REF", "ALT"):
            assert row[column] == "Not applicable"

    def test_a_configured_but_absent_ru_fasta_is_tolerated(self, tmp_path):
        source = write_advntr_output(tmp_path, ADVNTR_HEADER + CANONICAL_ROW)
        config = {"reference_data": {"code_adVNTR_RUs": str(tmp_path / "does_not_exist.fa")}}

        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=config)

        assert read_result(tmp_path).iloc[0]["RU"] == "Not applicable"


# ---------------------------------------------------------------------------
# The negative placeholder -- the normal outcome for most samples
# ---------------------------------------------------------------------------


class TestNegativePlaceholder:
    """Every "adVNTR found nothing" route must still write the placeholder row."""

    def test_a_header_only_file_produces_the_negative_placeholder(self, tmp_path):
        """The normal negative case: adVNTR ran, called nothing."""
        source = write_advntr_output(tmp_path, ADVNTR_HEADER)

        advntr.process_advntr_output(str(source), str(tmp_path), "output")

        assert_is_negative_placeholder(read_result(tmp_path))

    def test_a_file_whose_rows_are_all_filtered_out_produces_the_negative_placeholder(self, tmp_path):
        """``I10_2_A_LEN3`` is an in-frame 3 bp insertion, so the frameshift filter drops it."""
        source = write_advntr_output(tmp_path, ADVNTR_HEADER + "25561\tI10_2_A_LEN3\t11\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output")

        assert_is_negative_placeholder(read_result(tmp_path))

    def test_the_placeholder_warns_at_warning_level(self, tmp_path, caplog):
        source = write_advntr_output(tmp_path, ADVNTR_HEADER)

        with caplog.at_level(logging.WARNING, logger=advntr.logger.name):
            advntr.process_advntr_output(str(source), str(tmp_path), "output")

        warnings = [r for r in caplog.records if r.levelno >= logging.WARNING]
        assert warnings, "an empty adVNTR result must be announced above INFO"


# ---------------------------------------------------------------------------
# The three returns that write nothing at all
# ---------------------------------------------------------------------------


class TestPathsThatWriteNoFile:
    """
    Characterisation of the three ``return``-without-writing branches.

    Every one of these leaves ``pipeline.py`` pointing at a file that does not exist.
    The stage does not fail; the report simply loses its adVNTR half much later.
    """

    def test_a_missing_input_file_writes_nothing_and_logs_at_error(self, tmp_path, caplog):
        missing = tmp_path / "never_written.vcf"

        with caplog.at_level(logging.ERROR, logger=advntr.logger.name):
            advntr.process_advntr_output(str(missing), str(tmp_path), "output")

        assert not (tmp_path / f"output{RESULT_SUFFIX}").exists()
        errors = [r for r in caplog.records if r.levelno >= logging.ERROR]
        assert errors, "a missing adVNTR output must be reported at ERROR"
        assert "not found" in errors[0].getMessage()

    def test_a_zero_byte_input_file_writes_nothing_and_logs_at_error(self, tmp_path, caplog):
        source = write_advntr_output(tmp_path, "")

        with caplog.at_level(logging.ERROR, logger=advntr.logger.name):
            advntr.process_advntr_output(str(source), str(tmp_path), "output")

        assert not (tmp_path / f"output{RESULT_SUFFIX}").exists()
        assert [r for r in caplog.records if r.levelno >= logging.ERROR]

    @pytest.mark.parametrize(
        ("label", "blob"),
        [
            ("ragged_rows", ADVNTR_HEADER + f"25561\t{CANONICAL_VARIANT}\t11\t153.98\t0.0001\ta\tb\tc\n"),
            ("not_a_table", "this is not a table at all\njust prose\n"),
        ],
    )
    def test_a_malformed_input_file_writes_nothing_and_logs_at_error(self, tmp_path, caplog, label, blob):
        source = write_advntr_output(tmp_path, blob)

        with caplog.at_level(logging.ERROR, logger=advntr.logger.name):
            advntr.process_advntr_output(str(source), str(tmp_path), "output")

        assert not (tmp_path / f"output{RESULT_SUFFIX}").exists(), label
        assert [r for r in caplog.records if r.levelno >= logging.ERROR], label


# ---------------------------------------------------------------------------
# E5 -- compound variants must not destroy the whole sample
# ---------------------------------------------------------------------------


class TestCompoundVariants:
    """
    A ``State`` naming two insertions (``I9_2_A_LEN9&I50_2_A_LEN3``) used to raise inside
    ``advntr_processing_del``: the greedy ``extract("(LEN.*)")`` matched
    ``LEN9&I50_2_A_LEN3``, and splitting that on ``LEN`` yields *three* parts for a
    two-column assignment. The ``ValueError`` was swallowed by the broad handler, so the
    stage returned having written nothing -- discarding every other variant in the sample.
    """

    COMPOUND = "I9_2_A_LEN9&I50_2_A_LEN3"

    def test_a_compound_variant_does_not_discard_the_rest_of_the_sample(self, tmp_path):
        source = write_advntr_output(
            tmp_path,
            ADVNTR_HEADER + CANONICAL_ROW + f"25561\t{self.COMPOUND}\t11\t153.98\t0.0001\n",
        )

        advntr.process_advntr_output(str(source), str(tmp_path), "output")

        df = read_result(tmp_path)
        assert list(df.columns) == FINAL_COLUMNS
        assert CANONICAL_VARIANT in set(df["Variant"]), "the unrelated canonical variant must survive"

    def test_a_compound_variant_that_passes_the_filter_is_reported(self, tmp_path):
        source = write_advntr_output(tmp_path, ADVNTR_HEADER + "25561\tI9_2_A_LEN1&I50_2_A_LEN3\t11\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output")

        df = read_result(tmp_path)
        assert list(df["Variant"]) == ["I9_2_A_LEN1&I50_2_A_LEN3"]

    def test_a_compound_variant_is_annotated_part_by_part(self, tmp_path, ru_config):
        source = write_advntr_output(tmp_path, ADVNTR_HEADER + "25561\tI9_2_A_LEN1&I50_2_A_LEN3\t11\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=ru_config)

        row = read_result(tmp_path).iloc[0]
        assert row["RU"] == "2,2"
        assert row["POS"] == "9,50"

    def test_only_the_first_len_feeds_the_frameshift_length(self, tmp_path):
        """
        Characterisation, *not* an endorsement. The repaired extraction keeps the historic
        semantics -- ``Insertion_len`` is the first ``LEN`` in the state string, not the sum
        over all of its parts. For ``I9_2_A_LEN9&I50_2_A_LEN3`` the biologically meaningful
        inserted length is 12; this pins the 9 the code actually uses so that a deliberate
        change to summing shows up as a failing test rather than a silent reclassification.
        """
        df = pd.DataFrame(
            {
                "VID": [25561],
                "State": [self.COMPOUND],
                "NumberOfSupportingReads": [11],
                "MeanCoverage": [153.98],
                "Pvalue": [0.0001],
            }
        )

        processed = advntr.advntr_processing_ins(df.copy())
        # frame 9 is a multiple of 3, so the row is filtered out; recompute unfiltered.
        unfiltered = advntr.advntr_processing_ins(
            pd.DataFrame({**df.to_dict("list"), "State": ["I9_2_A_LEN1&I50_2_A_LEN3"]})
        )

        assert processed.empty, "first-LEN 9 is in frame, so the compound row is dropped"
        assert list(unfiltered["Insertion_len"]) == [1], "first LEN wins; the sum (4) is not used"


# ---------------------------------------------------------------------------
# Repeat-unit annotation helpers
# ---------------------------------------------------------------------------


class TestRepeatUnitAnnotation:
    def test_load_ru_sequences_strips_the_ru_prefix_and_joins_wrapped_lines(self, tmp_path):
        path = tmp_path / "ru.fa"
        path.write_text(">RU2\nACGT\nTTTT\n>RU7\nGGGG\n")

        assert advntr.load_ru_sequences(str(path)) == {"2": "ACGTTTTT", "7": "GGGG"}

    def test_an_insertion_is_annotated_with_the_reference_base_and_the_inserted_run(self, tmp_path):
        path = tmp_path / "ru.fa"
        path.write_text(">RU2\nACGTACGTAC\n")

        ru, pos, ref, alt = advntr.annotate_advntr_variants(pd.Series(["I3_2_T_LEN2"]), str(path))

        assert (ru, pos, ref, alt) == (["2"], ["3"], ["G"], ["GTT"])

    def test_a_deletion_is_annotated_with_the_preceding_base(self, tmp_path):
        path = tmp_path / "ru.fa"
        path.write_text(">RU2\nACGTACGTAC\n")

        ru, pos, ref, alt = advntr.annotate_advntr_variants(pd.Series(["D3_2"]), str(path))

        assert (ru, pos, ref, alt) == (["2"], ["3"], ["CG"], ["C"])

    def test_an_unparseable_state_string_annotates_as_dots(self, tmp_path):
        path = tmp_path / "ru.fa"
        path.write_text(">RU2\nACGTACGTAC\n")

        ru, pos, ref, alt = advntr.annotate_advntr_variants(pd.Series(["NOT_A_VARIANT"]), str(path))

        assert (ru, pos, ref, alt) == (["."], ["."], ["."], ["."])

    def test_an_unknown_repeat_unit_annotates_as_dots_but_keeps_the_ru_id(self, tmp_path):
        path = tmp_path / "ru.fa"
        path.write_text(">RU2\nACGTACGTAC\n")

        ru, pos, ref, alt = advntr.annotate_advntr_variants(pd.Series(["D3_99"]), str(path))

        assert (ru, pos, ref, alt) == (["99"], ["3"], ["."], ["."])
