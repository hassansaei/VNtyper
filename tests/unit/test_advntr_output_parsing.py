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
writing a file. A caller that gets no file does not see the real cause at that point; it
sees a misleading downstream failure much later, so which paths reach it is pinned.
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
    # MUC1 nomenclature, added for every caller in one pass. Nullable: written as
    # the empty string, never as "Not applicable", on rows that carry no name.
    "Nomenclature",
    "Nomenclature_Tier",
    "Nomenclature_Flags",
    "Ambiguity_Interval",
    "Repeat_Form",
    "Nomenclature_Note",
    "Nomenclature_Kestrel",
    "Nomenclature_adVNTR",
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


#: The nomenclature columns are nullable by contract: on a row that carries no name
#: they are written empty, never padded with the "Not applicable" the older columns
#: use. Read back with ``dtype=str`` an empty cell arrives as NaN.
NULLABLE_COLUMNS = (
    "Nomenclature",
    "Nomenclature_Tier",
    "Nomenclature_Flags",
    "Ambiguity_Interval",
    "Repeat_Form",
    "Nomenclature_Note",
    "Nomenclature_Kestrel",
    "Nomenclature_adVNTR",
)


def assert_is_negative_placeholder(df: pd.DataFrame) -> None:
    """Assert the frame is exactly the one-row negative placeholder the report expects."""
    assert list(df.columns) == FINAL_COLUMNS
    assert len(df) == 1
    row = df.iloc[0]
    assert row["VID"] == "Negative"
    for column in FINAL_COLUMNS[1:]:
        if column in NULLABLE_COLUMNS:
            assert pd.isna(row[column]) or row[column] == "", (
                f"{column} is nullable and must be empty, got {row[column]!r}"
            )
            continue
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


def test_header_rewrite_oserror_logs_and_returns_none(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A failed header rewrite remains observable without claiming a result."""
    source = write_advntr_output(tmp_path, ADVNTR_HEADER + CANONICAL_ROW)
    result_path = tmp_path / f"output{RESULT_SUFFIX}"
    real_open = open

    def fail_header_rewrite(path, mode="r", *args, **kwargs):
        if Path(path) == source and mode == "w":
            raise OSError("header rewrite blocked")
        return real_open(path, mode, *args, **kwargs)

    monkeypatch.setattr("builtins.open", fail_header_rewrite)

    with caplog.at_level(logging.ERROR, logger=advntr.logger.name):
        result = advntr.process_advntr_output(str(source), str(tmp_path), "output")

    assert result is None
    assert not result_path.exists()
    errors = [
        (record.levelno, record.getMessage())
        for record in caplog.records
        if record.name == advntr.logger.name and record.levelno >= logging.ERROR
    ]
    assert errors == [(logging.ERROR, "Error reading adVNTR output: header rewrite blocked")]


@pytest.mark.parametrize(
    ("failure_stage", "expected_message"),
    [
        ("dataframe_load", "Error loading data into DataFrame: unreadable result"),
        ("variant_processing", "Error during processing of deletions and insertions: unreadable result"),
    ],
)
def test_unreadable_advntr_output_logs_and_returns_none(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
    failure_stage: str,
    expected_message: str,
) -> None:
    """Each broad parsing catch returns None without claiming a result artifact."""
    source = write_advntr_output(tmp_path, ADVNTR_HEADER + CANONICAL_ROW)

    def unreadable(*_args, **_kwargs):
        raise OSError("unreadable result")

    if failure_stage == "dataframe_load":
        monkeypatch.setattr(advntr.pd, "read_csv", unreadable)
    else:
        monkeypatch.setattr(advntr, "advntr_processing_del", unreadable)

    cleanup_calls: list[tuple[str, str]] = []
    monkeypatch.setattr(
        advntr,
        "cleanup_files",
        lambda output, output_name: cleanup_calls.append((output, output_name)),
    )

    with caplog.at_level(logging.ERROR, logger=advntr.logger.name):
        assert advntr.process_advntr_output(str(source), str(tmp_path), "output") is None

    assert not (tmp_path / f"output{RESULT_SUFFIX}").exists()
    assert cleanup_calls == []
    errors = [
        record for record in caplog.records if record.name == advntr.logger.name and record.levelno >= logging.ERROR
    ]
    assert [record.levelno for record in errors] == [logging.ERROR]
    assert errors[0].getMessage() == expected_message


# ---------------------------------------------------------------------------
# E5 -- compound variants must not destroy the whole sample
# ---------------------------------------------------------------------------


class TestCompoundVariants:
    """
    A ``State`` naming two insertions (``I9_2_A_LEN9&I50_2_A_LEN3``) once raised inside
    ``advntr_processing_del``: the greedy ``extract("(LEN.*)")`` matched
    ``LEN9&I50_2_A_LEN3``, and splitting *that* on ``LEN`` yields *three* parts for a
    two-column assignment. The ``ValueError`` was swallowed by the broad handler, so the
    stage returned having written nothing -- discarding every other variant in the sample.

    The first repair was crash-only: the extraction stayed greedy and the split was bounded
    to one, which stopped the raise by giving the shape the same zero every other multi-part
    remainder already got. #192 replaced that zero with the sum over every ``LEN`` token,
    so the shape now carries a real length. See :class:`TestSummedInsertionLength` for the
    decision and :class:`TestInsertionLenAndFrameAfterFiltering` for the table it produces.
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

    @pytest.mark.parametrize(
        ("state", "reported"),
        [
            # 9 + 3 = 12, a multiple of 3: in frame, so both halves drop it.
            ("I9_2_A_LEN9&I50_2_A_LEN3", False),
            # 1 + 3 = 4 == 1 mod 3: a frameshift, so the insertion half reports it. Under
            # the pre-#192 zero this fell to the negative placeholder like the row above.
            ("I9_2_A_LEN1&I50_2_A_LEN3", True),
        ],
    )
    def test_the_two_len_shape_neither_raises_nor_is_forced_to_the_placeholder(self, tmp_path, caplog, state, reported):
        """
        The crash shape: two ``LEN`` tokens in one ``State``. Before the repair both filters
        raised ``ValueError`` and :func:`process_advntr_output` returned having written no
        file at all, leaving ``pipeline.py`` pointing at a path that does not exist. It now
        parses, and since #192 its two lengths are summed rather than read as zero -- so
        whether it reaches the report is decided by the #182 frameshift rule on the summed
        length, not by the parse collapsing.
        """
        source = write_advntr_output(tmp_path, ADVNTR_HEADER + f"25561\t{state}\t11\t153.98\t0.0001\n")

        with caplog.at_level(logging.ERROR, logger=advntr.logger.name):
            advntr.process_advntr_output(str(source), str(tmp_path), "output")

        df = read_result(tmp_path)
        if reported:
            assert list(df["Variant"]) == [state]
        else:
            assert_is_negative_placeholder(df)
        assert not [r for r in caplog.records if r.levelno >= logging.ERROR], "the parse must not fail any more"

    def test_a_compound_variant_is_annotated_part_by_part(self, tmp_path, ru_config):
        """``D2_2&I2_2_C_LEN5`` survives the insertion filter identically before and after."""
        source = write_advntr_output(tmp_path, ADVNTR_HEADER + "25561\tD2_2&I2_2_C_LEN5\t11\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=ru_config)

        row = read_result(tmp_path).iloc[0]
        assert row["Variant"] == "D2_2&I2_2_C_LEN5"
        assert row["RU"] == "2,2"
        assert row["POS"] == "2,2"


def state_frame(state: str) -> pd.DataFrame:
    """One adVNTR row carrying ``state``, in the shape ``pd.read_csv`` hands the filters."""
    return pd.DataFrame(
        {
            "VID": [25561],
            "State": [state],
            "NumberOfSupportingReads": [11],
            "MeanCoverage": [153.98],
            "Pvalue": [0.0001],
        }
    )


def reported_rows(state: str) -> int:
    """Rows ``state`` contributes to ``output_adVNTR_result.tsv``.

    ``process_advntr_output`` concatenates the deletion half and the insertion half, so a
    state is reported when *either* filter keeps it. Which half keeps it is invisible in
    the result file; whether one of them does is the whole of the reported outcome.
    """
    return len(advntr.advntr_processing_del(state_frame(state))) + len(advntr.advntr_processing_ins(state_frame(state)))


class TestSummedInsertionLength:
    """
    Specification (#192, decided 2026-08-06): every ``LEN`` token in a ``State`` counts.

    @hassansaei: "For compound adVNTR states [...] use the **sum** of inserted lengths and
    the **sum** of deleted lengths when computing the net length that feeds the pathogenic
    frameshift filter (``frame = abs(Insertion_len - Deletion_length)``, with the same
    3n+1 / 3n+2 rule as #182). Example: ``I9_2_A_LEN9&I50_2_A_LEN3`` to Insertion_len =
    9 + 3 = 12 (not first-LEN-only). [...] Do not keep first-LEN-wins as the defined
    semantics."

    **The diagnosis in #192 and in the #198 PR body is wrong, and is corrected here.**
    Both call the historic behaviour "first-LEN-wins". Measured against the real parser,
    *no input behaved that way*: the historic expression extracted greedily from the first
    ``LEN`` to the end of the string, so any material following that first ``LEN`` left a
    non-numeric remainder that ``pd.to_numeric(errors="coerce")`` turned into **zero**.
    ``I9_2_A_LEN9&I50_2_A_LEN3`` gave 0, not 9. The precise historic rule was: the value
    collapses to zero exactly when material follows the first ``LEN``; a compound whose
    only ``LEN`` is terminal parsed correctly. Where it collapsed, a pure-insertion
    compound got ``Net_indel_length == 0``, which is in frame and therefore not the
    pathogenic ``Delta % 3 == 1`` -- so those states were dropped in silence. That is a
    worse defect than the issue describes.
    """

    @pytest.mark.parametrize(
        ("state", "expected"),
        [
            # --- unchanged by this commit: a single terminal LEN already parsed ---
            ("I22_2_G_LEN1", 1),
            ("I50_2&I9_2_A_LEN3", 3),
            ("I50_2&D9_2&I80_2_A_LEN7", 7),
            ("D8_2&D9_2&I9_2_A_LEN9", 9),
            # --- unchanged by this commit: no LEN token at all ---
            ("D17_2&D18_2&D19_2&D20_2&D21_2", 0),  # the golden cohort's only compound call
            ("D50_2", 0),
            ("", 0),
            ("NOT_A_VARIANT", 0),
            ("I9_2_A_LENX", 0),  # a malformed LEN is not a LEN token
            # --- changed: material follows the first LEN, so this used to be 0 ---
            ("I9_2_A_LEN9&I50_2_A_LEN3", 12),  # his worked example: 9 + 3
            ("I9_2_A_LEN9&I50_2_A_LEN3&I80_2_A_LEN1", 13),
            ("I9_2_A_LEN2&D50_2", 2),
            ("I9_2_A_LEN9&", 9),  # a trailing & no longer collapses the value to 0
        ],
    )
    def test_every_len_token_is_summed(self, state, expected):
        assert advntr.sum_insertion_lengths(state) == expected

    @pytest.mark.parametrize("variant", [None, float("nan"), 12])
    def test_a_non_string_state_sums_to_zero_rather_than_raising(self, variant):
        """The historic pipeline tolerated a NaN ``Variant`` (``str.extract`` -> NaN -> 0).

        Characterisation of that tolerance, carried across deliberately so the helper is a
        total function. It is not a claim that a NaN ``State`` is supported input: the
        surrounding ``Deletion_length`` cast raises on one long before this is reached.
        """
        assert advntr.sum_insertion_lengths(variant) == 0

    def test_the_maintainers_worked_example_does_not_by_itself_change_survival(self):
        """Recorded so the example is not mistaken for the regression test.

        9 + 3 = 12 and 12 % 3 == 0, so this state is filtered out under BOTH semantics --
        it was dropped at ``Insertion_len == 0`` before and is dropped as an in-frame
        12 bp insertion now. It pins the arithmetic @hassansaei specified; it does not
        demonstrate a changed outcome. The three tests below do that.
        """
        assert advntr.sum_insertion_lengths("I9_2_A_LEN9&I50_2_A_LEN3") == 12
        assert 12 % 3 == 0
        assert reported_rows("I9_2_A_LEN9&I50_2_A_LEN3") == 0

    def test_a_compound_insertion_summation_admits_was_dropped_before(self):
        """The insertion path: summation can only ADD rows there, never remove them.

        ``Net_indel_length = Insertion_len - Deletion_length``; ``Deletion_length`` is
        ``count("D") == 0`` here; the #182 filter keeps ``Net_indel_length % 3 == 1``,
        which for a pure insertion means a length of 3n+1::

            I9_2_A_LEN9&I50_2_A_LEN1
              before:    Insertion_len = 0,  Delta = 0,  in frame     -> DROPPED
              summation: Insertion_len = 10, Delta = +10, 10 % 3 == 1 -> KEPT

        Every state whose value collapsed had ``Insertion_len == 0``, hence
        ``Delta == -Deletion_length``, which for a pure insertion is 0 -- in frame, so no
        collapsed pure-insertion state survived at all. On the insertion side summation is
        therefore monotone-additive: it makes real pathogenic compound calls visible and
        cannot lose one.
        """
        kept = advntr.advntr_processing_ins(state_frame("I9_2_A_LEN9&I50_2_A_LEN1"))

        assert len(kept) == 1
        assert kept["Insertion_len"].iloc[0] == 10
        assert list(kept["frame"]) == ["10"]

    def test_the_deletion_half_can_lose_a_row_the_insertion_half_then_keeps(self):
        """The deletion path is NOT monotone -- but this state is still reported.

        ``I9_2_A_LEN3&D50_2&D51_2`` (``Deletion_length = count("D") = 2``)::

            before:    Insertion_len = 0, Delta = 0-2 = -2, a 3n+2 net deletion -> kept by del
            summation: Insertion_len = 3, Delta = 3-2 = +1, a 3n+1 net insertion -> kept by ins

        Both are ``Delta % 3 == 1``, so both are the pathogenic frame; only the sign of
        the net change moves. The row moves from the deletion half to the insertion half.
        Because
        ``process_advntr_output`` concatenates the two halves, the reported output is
        unchanged; only which filter accounts for the row changes, and that is not a
        column of ``output_adVNTR_result.tsv``.

        The task brief predicted this state as the "a reported call disappears" case. It
        is not: measured, the insertion half picks it up. The state that genuinely stops
        being reported is the net-in-frame one below.
        """
        assert len(advntr.advntr_processing_del(state_frame("I9_2_A_LEN3&D50_2&D51_2"))) == 0
        assert len(advntr.advntr_processing_ins(state_frame("I9_2_A_LEN3&D50_2&D51_2"))) == 1
        assert reported_rows("I9_2_A_LEN3&D50_2&D51_2") == 1

    def test_a_net_in_frame_compound_stops_being_reported_at_all(self):
        """The one direction that removes a reported call. Asserted, not discovered later.

        ``I9_2_A_LEN2&D50_2&D51_2`` -- 2 inserted bases against 2 deleted bases::

            before:    Insertion_len = 0, Delta = 0-2 = -2, -2 % 3 == 1 -> REPORTED
            summation: Insertion_len = 2, Delta = 2-2 =  0,  0 % 3 == 0 -> DROPPED

        This is a reported adVNTR call disappearing, so it is pinned here and called out
        in the commit message and the PR body. @hassansaei's decision on #192 authorises
        it: the quantity the frameshift filter tests is the *net* indel length, and a net
        change of 0 is in frame, so the row was only ever reported because the inserted
        bases were being read as zero.
        """
        assert reported_rows("I9_2_A_LEN2&D50_2&D51_2") == 0


class TestInsertionLenAndFrameAfterFiltering:
    """
    The whole table each filter produces for a state: derivation, frame, and survival.

    Two rules meet here, both decided, so the table is a specification rather than a
    characterisation:

    * the ``Insertion_len`` **derivation** is specified by #192 (see
      :class:`TestSummedInsertionLength`);
    * the **pathogenic-frame rule** is specified by #182: a row survives exactly when its
      *signed* net change ``Net_indel_length = Insertion_len - Deletion_length`` satisfies
      ``Net_indel_length % 3 == 1``. The deletion arm serves the negative half of that
      (a net loss of 3n+2 bases) and the insertion arm the positive half (a net gain of
      3n+1 bases).

    ``frame`` is ``|Net_indel_length|`` -- the magnitude only. It is what the configurable
    ``3n+1``/``3n+2`` series are matched against, and each arm tests the sign separately.
    An earlier revision of this class recorded the ``>= 1`` presence guards as a "live
    defect" left characterised; those guards are gone, replaced by the sign test, and the
    row this table used to pin as reported through the wrong arm
    (``D49_2&I49_2_A_LEN12``) is corrected below.

    Every expectation was re-measured against the real functions at this commit.
    """

    #: ``(state, function, surviving Insertion_len values, surviving frame values)``.
    MEASURED_BEHAVIOUR = [
        # A single-part insertion. Unchanged by #192: one terminal LEN always parsed.
        ("I22_2_G_LEN1", "ins", [1], ["1"]),  # Delta = +1
        ("I22_2_G_LEN1", "del", [], []),
        # A compound whose only LEN sits in the last part. Also unchanged by #192.
        ("D8_2&D9_2&I9_2_A_LEN9", "ins", [9], ["7"]),  # 9 - 2 = +7, 7 % 3 == 1
        ("D2_2&I2_2_C_LEN5", "ins", [5], ["4"]),  # 5 - 1 = +4, 4 % 3 == 1
        # 12 inserted against 1 deleted -> Delta = +11, and 11 % 3 == 2, so this is a
        # frameshift in the frame VNtyper does not report. It used to be reported by the
        # DELETION arm, on the magnitude |+11| = 11 being in the 3n+2 series while the arm
        # guarded only on Deletion_length >= 1 -- a net insertion reported as a deletion
        # call. Both arms now drop it. (This state is a real one: it is named in
        # advntr_config.json's Polymorphic_Call list.)
        ("D49_2&I49_2_A_LEN12", "del", [], []),
        ("D49_2&I49_2_A_LEN12", "ins", [], []),
        # No LEN token at all -- 0 before and after.
        ("D58_2&D59_2", "del", [0], ["2"]),  # Delta = -2, -2 % 3 == 1
        ("D8_2", "ins", [], []),
        ("NOT_A_VARIANT", "ins", [], []),
        ("NOT_A_VARIANT", "del", [], []),
        # #192 changed these. Before, material after the first LEN collapsed the value to
        # zero; now the LEN tokens are summed, so the net length decides.
        ("I9_2_A_LEN2&D50_2", "ins", [2], ["1"]),  # was [] -- 2 inserted, 1 deleted, Delta = +1
        ("I9_2_A_LEN2&D50_2", "del", [], []),  # unchanged: a net insertion is never a deletion call
        ("I9_2_A_LEN2&D50_2&D51_2", "ins", [], []),  # unchanged: Delta = 0 is in frame
        ("I9_2_A_LEN2&D50_2&D51_2", "del", [], []),  # was [0]/["2"] -- Delta = 0 is in frame
        ("I9_2_A_LEN9&I50_2_A_LEN3", "ins", [], []),  # Delta = +12, in frame
        ("I9_2_A_LEN1&I50_2_A_LEN3", "ins", [4], ["4"]),  # Delta = +4, the pathogenic frame
        # Mixed states in the non-pathogenic frame, both signs. Neither arm may take them.
        ("I9_2_A_LEN3&D50_2", "del", [], []),  # Delta = +2; the branch-introduced case
        ("I9_2_A_LEN3&D50_2", "ins", [], []),
        ("D8_2&D9_2&I9_2_A_LEN1", "ins", [], []),  # Delta = -1
        ("D8_2&D9_2&I9_2_A_LEN1", "del", [], []),
    ]

    @pytest.mark.parametrize(("state", "which", "lengths", "frames"), MEASURED_BEHAVIOUR)
    def test_the_filters_produce_the_measured_insertion_len_and_frame(self, state, which, lengths, frames):
        run = advntr.advntr_processing_ins if which == "ins" else advntr.advntr_processing_del

        out = run(state_frame(state))

        assert list(out["Insertion_len"]) == lengths, f"{state} ({which}) Insertion_len"
        assert list(out["frame"]) == frames, f"{state} ({which}) frame"

    @pytest.mark.parametrize(
        ("state", "delta"),
        [
            # Net insertions reported by the deletion arm on |Delta| alone.
            ("I9_2_A_LEN3&D50_2", 2),
            ("D49_2&I49_2_A_LEN12", 11),
            # A net deletion reported by the insertion arm on |Delta| alone.
            ("D8_2&D9_2&I9_2_A_LEN1", -1),
        ],
    )
    def test_a_mixed_state_in_the_other_frame_reaches_the_result_file_from_neither_arm(self, tmp_path, state, delta):
        """End to end: the corrected verdict is what ``output_adVNTR_result.tsv`` shows.

        Every state here has ``Delta % 3 == 2`` -- a genuine frameshift, but not the
        pathogenic +1 frame -- and each used to be admitted by the arm whose *magnitude*
        series ``|Delta|`` happened to land in, because that arm's only other guard was a
        presence check a mixed state satisfies on both sides. The filter unit tests assert
        the arms; this asserts the file a clinician's report is built from.
        """
        assert delta % 3 == 2, "these states are the non-pathogenic frame"
        source = write_advntr_output(tmp_path, ADVNTR_HEADER + f"25561\t{state}\t11\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output")

        assert_is_negative_placeholder(read_result(tmp_path))

    def test_a_mixed_state_in_the_pathogenic_frame_still_reaches_the_result_file(self, tmp_path):
        """The control for the test above: ``Delta = +1`` from a compound state is reported.

        ``I9_2_A_LEN3&D50_2&D51_2`` -- 3 inserted against 2 deleted. Without this, the
        test above would also pass against a filter that had simply stopped reporting
        mixed states altogether.
        """
        source = write_advntr_output(tmp_path, ADVNTR_HEADER + "25561\tI9_2_A_LEN3&D50_2&D51_2\t11\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output")

        assert list(read_result(tmp_path)["Variant"]) == ["I9_2_A_LEN3&D50_2&D51_2"]


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
