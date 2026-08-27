# tests/unit/test_advntr_flagging_rules.py

"""
Contract tests for the adVNTR flagging rules in ``vntyper/modules/advntr/advntr_config.json``.

Those rule strings are handed to :func:`vntyper.scripts.flagging.evaluate_condition`, which
``eval``s them with the row's column names as locals. Two failure modes follow, and neither
of them is loud:

* a name that is not a column raises ``NameError``, which the evaluator downgrades to a
  ``logger.warning`` and a ``False`` result -- the rule is disabled and nothing fails;
* a comparison against the wrong *type* is simply ``False`` for every row -- the rule is
  dead and nothing fails.

This repository has already been bitten by the first: ``Poylmorhic_Call`` sat misspelled in
this file until commit ``59a85de``. The tests below retire both classes for this module.

Why it matters downstream: ``report_config.json`` reads the resulting ``Flag`` column and
classifies an adVNTR call as ``positive`` when ``Flag == "Not flagged"`` and
``positive flagged`` otherwise. A dead rule therefore reports a suspect call as a clean
positive.
"""

import ast
import json
from pathlib import Path

import pandas as pd
import pytest

from vntyper.modules.advntr import advntr_genotyping as advntr
from vntyper.scripts.flagging import evaluate_condition

pytestmark = pytest.mark.unit

#: Names the evaluator injects itself, so they are legal without being columns.
EVALUATOR_BUILTINS = {"regex_match"}

#: The flag names this module is expected to emit. Pinning them catches a typo in a rule
#: *key* -- the exact defect that shipped as ``Poylmorhic_Call``.
EXPECTED_FLAG_NAMES = {"Low_Coverage", "Repeat_Unit_7", "Polymorphic_Call"}

ADVNTR_HEADER = "#VID\tState\tNumberOfSupportingReads\tMeanCoverage\tPvalue\n"


@pytest.fixture(scope="module")
def advntr_config_on_disk() -> dict:
    """The shipped ``advntr_config.json``, read from disk rather than from the import-time global."""
    path = Path(advntr.__file__).parent / "advntr_config.json"
    assert path.exists(), f"advntr_config.json not found at {path}"
    with path.open() as handle:
        return json.load(handle)


@pytest.fixture(scope="module")
def flagging_rules(advntr_config_on_disk) -> dict:
    return advntr_config_on_disk["flagging_rules"]


@pytest.fixture
def ru_config(tmp_path: Path) -> dict:
    """A main-config carrying a minimal RU FASTA, so the RU/POS/REF/ALT columns get built."""
    fasta = tmp_path / "code-adVNTR_RUs.fa"
    fasta.write_text("".join(f">RU{n}\n{'ACGTACGTAC' * 7}\n" for n in (2, 6, 7)))
    return {"reference_data": {"code_adVNTR_RUs": str(fasta)}}


def free_names(expression: str) -> set[str]:
    """Every bare name an expression reads, i.e. every name ``eval`` must resolve."""
    tree = ast.parse(expression, mode="eval")
    return {node.id for node in ast.walk(tree) if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)}


# ---------------------------------------------------------------------------
# E3 -- every name in every rule must resolve to a column the parser produces
# ---------------------------------------------------------------------------


class TestEveryRuleNameResolves:
    def test_the_scan_found_a_plausible_number_of_rules(self, flagging_rules):
        """Guard against a vacuous pass: everything below iterates this dict."""
        assert len(flagging_rules) >= 3, f"expected at least 3 adVNTR flagging rules, found {len(flagging_rules)}"
        assert all(isinstance(expression, str) and expression.strip() for expression in flagging_rules.values())

    def test_the_flag_names_are_the_expected_set(self, flagging_rules):
        """A typo in a rule *key* renames the flag silently; ``Poylmorhic_Call`` did exactly that."""
        assert set(flagging_rules) == EXPECTED_FLAG_NAMES

    def test_the_polymorphic_call_list_is_the_size_the_cleanup_left(self, flagging_rules):
        """#267 removed 7 unreachable entries and 1 duplicate from the 32 that shipped.

        ``tests/unit/test_advntr_polymorphic_calls.py`` owns the reachability argument and
        the provenance record; this is the cheap tripwire in the file that already reads
        these rules, so an edit to the list is noticed here even if that file is not run.
        """
        node = ast.parse(flagging_rules["Polymorphic_Call"], mode="eval").body
        states = ast.literal_eval(node.comparators[0])

        assert len(states) == 24
        assert len(set(states)) == 24

    def test_every_rule_is_syntactically_valid_python(self, flagging_rules):
        for name, expression in flagging_rules.items():
            try:
                ast.parse(expression, mode="eval")
            except SyntaxError as exc:  # pragma: no cover - only on a broken config
                pytest.fail(f"flagging rule {name!r} is not parseable: {exc}")

    def test_the_scan_finds_names_to_check(self, flagging_rules):
        """Second vacuity guard: an expression scan that matches nothing proves nothing."""
        all_names = set().union(*(free_names(expression) for expression in flagging_rules.values()))
        columns = all_names - EVALUATOR_BUILTINS

        assert len(columns) >= 3, f"expected the rules to reference at least 3 columns, found {sorted(columns)}"

    def test_every_name_a_rule_reads_is_a_column_the_parser_produces(self, flagging_rules, produced_columns):
        """
        The real guard. ``produced_columns`` is captured from the DataFrame the parser
        actually hands to ``add_flags``, so renaming a column in ``advntr_genotyping.py``
        without updating this JSON fails here instead of silently disabling a flag.
        """
        for name, expression in flagging_rules.items():
            unresolved = free_names(expression) - EVALUATOR_BUILTINS - produced_columns
            assert not unresolved, (
                f"flagging rule {name!r} reads {sorted(unresolved)}, which the adVNTR parser "
                f"never produces. Available columns: {sorted(produced_columns)}. "
                "A missing name evaluates to False, so this rule is silently disabled."
            )

    def test_every_name_a_rule_reads_resolves_without_an_ru_fasta(
        self, flagging_rules, produced_columns_without_ru_fasta
    ):
        """No missing column may silently disable a rule on the no-FASTA branch."""
        for name, expression in flagging_rules.items():
            unresolved = free_names(expression) - EVALUATOR_BUILTINS - produced_columns_without_ru_fasta
            assert not unresolved, (
                f"flagging rule {name!r} reads {sorted(unresolved)}, which the adVNTR parser "
                "does not produce when no RU FASTA resolves"
            )


@pytest.fixture
def produced_columns(tmp_path, ru_config) -> set[str]:
    """
    The exact column set ``add_flags`` is called with, captured from a real parser run.

    Derived rather than hardcoded: a hardcoded list would drift out of step with the
    production code and the test would keep passing while the flags went dead.
    """
    captured: dict[str, set[str]] = {}
    real_add_flags = None

    import vntyper.scripts.flagging as flagging_module

    real_add_flags = flagging_module.add_flags

    def spy(df, flag_rules, duplicates_config=None):
        captured["columns"] = set(df.columns)
        return real_add_flags(df, flag_rules, duplicates_config)

    flagging_module.add_flags = spy
    try:
        source = tmp_path / "output_adVNTR.vcf"
        source.write_text(ADVNTR_HEADER + "25561\tI22_2_G_LEN1\t11\t153.98\t0.0001\n")
        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=ru_config)
    finally:
        flagging_module.add_flags = real_add_flags

    assert "columns" in captured, "the parser never reached the flagging step"
    return captured["columns"]


@pytest.fixture
def produced_columns_without_ru_fasta(tmp_path, monkeypatch) -> set[str]:
    """Capture the columns flagging receives when no RU FASTA resolves."""
    captured: dict[str, set[str]] = {}

    import vntyper.scripts.flagging as flagging_module

    real_add_flags = flagging_module.add_flags

    def spy(df, flag_rules, duplicates_config=None):
        captured["columns"] = set(df.columns)
        return real_add_flags(df, flag_rules, duplicates_config)

    monkeypatch.setattr(flagging_module, "add_flags", spy)
    source = tmp_path / "output_adVNTR.vcf"
    source.write_text(ADVNTR_HEADER + "25561\tI22_2_G_LEN1\t11\t153.98\t0.0001\n")
    advntr.process_advntr_output(str(source), str(tmp_path), "output", config=None)

    assert "columns" in captured, "the parser never reached the flagging step"
    return captured["columns"]


# ---------------------------------------------------------------------------
# E2 -- the Repeat_Unit_7 rule compares a string against an int
# ---------------------------------------------------------------------------


class TestRepeatUnitSeven:
    """
    ``RU`` is assembled by ``",".join(ru_parts)`` in :func:`annotate_advntr_variants`, so it
    is always a ``str``. The shipped rule read ``RU == 7``, and ``"7" == 7`` is ``False`` in
    Python for every possible row -- the rule could never fire.
    """

    def test_the_annotator_produces_repeat_unit_ids_as_strings(self, tmp_path):
        """This is *why* the comparison must be against a string literal."""
        fasta = tmp_path / "ru.fa"
        fasta.write_text(">RU7\nACGTACGTAC\n")

        ru_values, _, _, _ = advntr.annotate_advntr_variants(pd.Series(["I3_7_A_LEN1"]), str(fasta))

        assert ru_values == ["7"]
        assert all(isinstance(value, str) for value in ru_values)

    def test_the_rule_fires_for_a_repeat_unit_seven_row(self, flagging_rules):
        row = pd.Series({"RU": "7", "Variant": "I3_7_A_LEN1", "NumberOfSupportingReads": 42})

        assert evaluate_condition(row, flagging_rules["Repeat_Unit_7"]) is True

    def test_the_rule_does_not_fire_for_another_repeat_unit(self, flagging_rules):
        row = pd.Series({"RU": "2", "Variant": "I3_2_A_LEN1", "NumberOfSupportingReads": 42})

        assert evaluate_condition(row, flagging_rules["Repeat_Unit_7"]) is False

    def test_a_compound_call_spanning_two_units_does_not_fire(self, flagging_rules):
        row = pd.Series({"RU": "7,2", "Variant": "I3_7_A_LEN1&I9_2_A_LEN1", "NumberOfSupportingReads": 42})

        assert evaluate_condition(row, flagging_rules["Repeat_Unit_7"]) is False

    #: The expression as it shipped before the repair. Kept as a literal so the "it could
    #: never fire" claim is *asserted* rather than asserted-about-the-past in prose.
    HISTORIC_EXPRESSION = "RU == 7"

    #: Row shapes spanning every ``RU`` the annotator can emit: the target unit, another
    #: unit, compound joins in both orders, a look-alike, the no-config placeholder and a
    #: missing value.
    RU_VALUES = ["7", "2", "6", "7,2", "2,7", "77", "7.0", "", "Not applicable", None]

    def probe_rows(self) -> pd.DataFrame:
        """The full result schema, so a rule reading any column resolves it."""
        import itertools

        return pd.DataFrame(
            [
                {
                    "VID": "25561",
                    "Variant": variant,
                    "NumberOfSupportingReads": reads,
                    "MeanCoverage": 153.98,
                    "Pvalue": 0.0001,
                    "RU": ru,
                    "POS": "1",
                    "REF": "A",
                    "ALT": "AA",
                }
                for ru, variant, reads in itertools.product(
                    self.RU_VALUES,
                    ["I22_2_G_LEN1", "I10_2_A_LEN1", "I26_7_A_LEN25", "I3_7_A_LEN1"],
                    [3, 9, 10, 11, 42],
                )
            ]
        )

    def test_the_historic_expression_could_not_fire_for_any_row(self, flagging_rules):
        """``RU`` is always a ``str``, so ``RU == 7`` was dead code, not a live rule."""
        rows = self.probe_rows()
        assert len(rows) >= 100, f"vacuity guard: only {len(rows)} probe rows"

        fired = [
            index for index, row in rows.iterrows() if evaluate_condition(row, self.HISTORIC_EXPRESSION) is not False
        ]

        assert not fired, f"the historic integer comparison fired for rows {fired}"

    def test_the_repair_only_ever_adds_a_flag_and_moves_no_other_column(self, flagging_rules):
        """
        The safety argument for shipping this change: reviving the rule is monotone. No row
        loses a flag, the only flag any row gains is ``Repeat_Unit_7``, and no genotype
        column moves -- ``add_flags`` appends a column and filters nothing.
        """
        from vntyper.scripts.flagging import add_flags

        historic_rules = {**flagging_rules, "Repeat_Unit_7": self.HISTORIC_EXPRESSION}
        rows = self.probe_rows()

        before = add_flags(rows.copy(), historic_rules)
        after = add_flags(rows.copy(), flagging_rules)

        def flags(frame, index):
            return set(str(frame.iloc[index]["Flag"]).split(", ")) - {"Not flagged"}

        gained = set()
        for index in range(len(rows)):
            was, now = flags(before, index), flags(after, index)
            assert not was - now, f"row {index} ({rows.iloc[index]['RU']!r}) lost flags {sorted(was - now)}"
            gained |= now - was

        assert gained == {"Repeat_Unit_7"}, f"unexpected new flags: {sorted(gained)}"
        pd.testing.assert_frame_equal(before.drop(columns=["Flag"]), after.drop(columns=["Flag"]))

    def test_end_to_end_a_repeat_unit_seven_call_is_flagged_in_the_result_file(self, tmp_path, ru_config):
        source = tmp_path / "output_adVNTR.vcf"
        source.write_text(ADVNTR_HEADER + "25561\tI26_7_A_LEN1\t11\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=ru_config)

        result = pd.read_csv(tmp_path / "output_adVNTR_result.tsv", sep="\t", dtype=str)
        assert result.iloc[0]["RU"] == "7"
        assert "Repeat_Unit_7" in result.iloc[0]["Flag"]

    def test_a_repeat_unit_seven_call_is_flagged_without_an_ru_fasta(self, tmp_path):
        """Repeat-unit flagging depends on the state string, not a reference file."""
        source = tmp_path / "output_adVNTR.vcf"
        source.write_text(ADVNTR_HEADER + "25561\tI26_7_A_LEN1\t42\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=None)

        result = pd.read_csv(tmp_path / "output_adVNTR_result.tsv", sep="\t", dtype=str)
        row = result.iloc[0]
        assert row["RU"] == "7"
        assert row["POS"] == "26"
        assert row["REF"] == "Not applicable"
        assert row["ALT"] == "Not applicable"
        assert "Repeat_Unit_7" in row["Flag"]

    def test_a_compound_repeat_unit_seven_call_stays_unflagged_without_an_ru_fasta(self, tmp_path):
        """A compound RU annotation remains ``7,2`` and does not become an RU-7 flag."""
        source = tmp_path / "output_adVNTR.vcf"
        source.write_text(ADVNTR_HEADER + "25561\tI3_7_A_LEN1&I9_2_A_LEN3\t42\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=None)

        row = pd.read_csv(tmp_path / "output_adVNTR_result.tsv", sep="\t", dtype=str).iloc[0]
        assert row["RU"] == "7,2"
        assert row["POS"] == "3,9"
        assert row["Flag"] == "Not flagged"


# ---------------------------------------------------------------------------
# No rule may be dead
# ---------------------------------------------------------------------------


class TestNoRuleIsDead:
    """
    The generalisation of E2: for every configured rule there must exist a row that makes
    it fire. A rule that cannot fire for any input is indistinguishable from a deleted one.
    """

    #: One row per rule, each designed to trigger exactly that rule.
    TRIGGER_ROWS = [
        {"Variant": "I3_2_A_LEN1", "NumberOfSupportingReads": 3, "RU": "2"},  # Low_Coverage
        {"Variant": "I3_7_A_LEN1", "NumberOfSupportingReads": 42, "RU": "7"},  # Repeat_Unit_7
        {"Variant": "I10_2_A_LEN1", "NumberOfSupportingReads": 42, "RU": "2"},  # Polymorphic_Call
    ]

    def test_every_configured_rule_fires_for_at_least_one_row(self, flagging_rules):
        assert len(self.TRIGGER_ROWS) >= len(flagging_rules), "add a trigger row for each new rule"

        frame = pd.DataFrame(self.TRIGGER_ROWS)
        fired = {
            name
            for name, expression in flagging_rules.items()
            if frame.apply(lambda row, e=expression: evaluate_condition(row, e), axis=1).any()
        }

        assert fired == set(flagging_rules), f"these rules can never fire: {sorted(set(flagging_rules) - fired)}"


class TestTheOtherRules:
    def test_low_coverage_fires_below_ten_supporting_reads(self, flagging_rules):
        rule = flagging_rules["Low_Coverage"]

        assert evaluate_condition(pd.Series({"NumberOfSupportingReads": 9}), rule) is True
        assert evaluate_condition(pd.Series({"NumberOfSupportingReads": 10}), rule) is False

    def test_polymorphic_call_fires_for_a_listed_variant(self, flagging_rules):
        rule = flagging_rules["Polymorphic_Call"]

        assert evaluate_condition(pd.Series({"Variant": "I10_2_A_LEN1"}), rule) is True
        assert evaluate_condition(pd.Series({"Variant": "I22_2_G_LEN1"}), rule) is False

    def test_several_rules_combine_into_one_comma_separated_flag(self, tmp_path, ru_config):
        source = tmp_path / "output_adVNTR.vcf"
        # I26_7_A_LEN25 is both an RU-7 call and a listed polymorphic call, at low coverage.
        source.write_text(ADVNTR_HEADER + "25561\tI26_7_A_LEN25\t4\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=ru_config)

        flag = pd.read_csv(tmp_path / "output_adVNTR_result.tsv", sep="\t", dtype=str).iloc[0]["Flag"]
        assert set(flag.split(", ")) == EXPECTED_FLAG_NAMES

    def test_an_unremarkable_call_is_not_flagged(self, tmp_path, ru_config):
        source = tmp_path / "output_adVNTR.vcf"
        source.write_text(ADVNTR_HEADER + "25561\tI22_2_G_LEN1\t42\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=ru_config)

        result = pd.read_csv(tmp_path / "output_adVNTR_result.tsv", sep="\t", dtype=str)
        assert result.iloc[0]["Flag"] == "Not flagged"


class TestFlaggingReadsTheOriginalConfigGlobal:
    """
    H1: ``process_advntr_output`` reads ``advntr_config["flagging_rules"]`` -- the *raw*
    import-time global -- while the frameshift filter reads the derived ``advntr_settings``.
    A test that patches the wrong one is a silent no-op.
    """

    def test_patching_the_raw_config_global_changes_the_flags(self, tmp_path, monkeypatch, ru_config):
        monkeypatch.setattr(
            advntr,
            "advntr_config",
            {"flagging_rules": {"Custom_Flag": "NumberOfSupportingReads > 1"}},
        )
        source = tmp_path / "output_adVNTR.vcf"
        source.write_text(ADVNTR_HEADER + "25561\tI22_2_G_LEN1\t42\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=ru_config)

        result = pd.read_csv(tmp_path / "output_adVNTR_result.tsv", sep="\t", dtype=str)
        assert result.iloc[0]["Flag"] == "Custom_Flag"

    def test_patching_the_derived_settings_global_does_not_change_the_flags(self, tmp_path, monkeypatch, ru_config):
        monkeypatch.setattr(advntr, "advntr_settings", {"flagging_rules": {"Custom_Flag": "True"}})
        source = tmp_path / "output_adVNTR.vcf"
        source.write_text(ADVNTR_HEADER + "25561\tI22_2_G_LEN1\t42\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=ru_config)

        result = pd.read_csv(tmp_path / "output_adVNTR_result.tsv", sep="\t", dtype=str)
        assert result.iloc[0]["Flag"] == "Not flagged"

    def test_an_empty_rule_set_leaves_the_flag_column_as_not_applicable(self, tmp_path, monkeypatch, ru_config):
        monkeypatch.setattr(advntr, "advntr_config", {"flagging_rules": {}})
        source = tmp_path / "output_adVNTR.vcf"
        source.write_text(ADVNTR_HEADER + "25561\tI22_2_G_LEN1\t42\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=ru_config)

        result = pd.read_csv(tmp_path / "output_adVNTR_result.tsv", sep="\t", dtype=str)
        assert result.iloc[0]["Flag"] == "Not applicable"
