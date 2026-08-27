# tests/unit/test_advntr_flagging_rules.py

"""
Contract tests for the adVNTR flagging rules in ``vntyper/modules/advntr/advntr_config.json``.

Those rules are validated as structured comparator data before any row is processed. Two
historical failure modes motivate this live contract:

* a name that is not a column used to disable a string expression without failing the run;
* a comparison against the wrong *type* is simply ``False`` for every row -- the rule is
  dead and nothing fails.

This repository has already been bitten by the first: ``Poylmorhic_Call`` sat misspelled in
this file until commit ``59a85de``. The tests below retire both classes for this module.

Why it matters downstream: ``report_config.json`` reads the resulting ``Flag`` column and
classifies an adVNTR call as ``positive`` when ``Flag == "Not flagged"`` and
``positive flagged`` otherwise. A dead rule therefore reports a suspect call as a clean
positive.
"""

import json
from collections.abc import Mapping
from pathlib import Path

import pandas as pd
import pytest

from vntyper.modules.advntr import advntr_genotyping as advntr
from vntyper.scripts.flagging import add_flags

pytestmark = pytest.mark.unit

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


def configured_columns(configured: object) -> set[str]:
    """Extract column operands from one structured rule without parsing source text."""
    assert isinstance(configured, Mapping)
    assert set(configured) == {"all"}
    predicates = configured["all"]
    assert isinstance(predicates, list) and predicates
    columns: set[str] = set()
    for predicate in predicates:
        assert isinstance(predicate, Mapping)
        assert set(predicate) == {"left", "operator", "right"}
        for side in ("left", "right"):
            operand = predicate[side]
            assert isinstance(operand, Mapping)
            assert len(operand) == 1
            if "column" in operand:
                column = operand["column"]
                assert isinstance(column, str) and column
                columns.add(column)
            else:
                assert "literal" in operand
    return columns


def flag_fires(name: str, configured: object, row: Mapping[str, object]) -> bool:
    """Evaluate one rule through ``add_flags`` and report whether its flag was emitted."""
    result = add_flags(pd.DataFrame([dict(row)]), {name: configured})
    return result.iloc[0]["Flag"] == name


# ---------------------------------------------------------------------------
# E3 -- every name in every rule must resolve to a column the parser produces
# ---------------------------------------------------------------------------


class TestEveryRuleNameResolves:
    def test_the_scan_found_a_plausible_number_of_rules(self, flagging_rules):
        """Guard against a vacuous pass: everything below iterates this dict."""
        assert len(flagging_rules) >= 3, f"expected at least 3 adVNTR flagging rules, found {len(flagging_rules)}"
        assert all(isinstance(configured, dict) and configured for configured in flagging_rules.values())

    def test_the_flag_names_are_the_expected_set(self, flagging_rules):
        """A typo in a rule *key* renames the flag silently; ``Poylmorhic_Call`` did exactly that."""
        assert set(flagging_rules) == EXPECTED_FLAG_NAMES

    def test_the_polymorphic_call_list_is_the_size_the_cleanup_left(self, flagging_rules):
        """#267 removed 7 unreachable entries and 1 duplicate from the 32 that shipped.

        ``tests/unit/test_advntr_polymorphic_calls.py`` owns the reachability argument and
        the provenance record; this is the cheap tripwire in the file that already reads
        these rules, so an edit to the list is noticed here even if that file is not run.
        """
        predicate = flagging_rules["Polymorphic_Call"]["all"][0]
        assert predicate["left"] == {"column": "Variant"}
        assert predicate["operator"] == "in"
        states = predicate["right"]["literal"]

        assert len(states) == 24
        assert len(set(states)) == 24

    def test_the_scan_finds_columns_to_check(self, flagging_rules):
        """Second vacuity guard: a structured scan that matches nothing proves nothing."""
        columns = set().union(*(configured_columns(configured) for configured in flagging_rules.values()))

        assert len(columns) >= 3, f"expected the rules to reference at least 3 columns, found {sorted(columns)}"

    def test_every_name_a_rule_reads_is_a_column_the_parser_produces(self, flagging_rules, produced_columns):
        """
        The real guard. ``produced_columns`` is captured from the DataFrame the parser
        actually hands to ``add_flags``, so renaming a column in ``advntr_genotyping.py``
        without updating this JSON fails here instead of silently disabling a flag.
        """
        for name, configured in flagging_rules.items():
            unresolved = configured_columns(configured) - produced_columns
            assert not unresolved, (
                f"flagging rule {name!r} reads {sorted(unresolved)}, which the adVNTR parser "
                f"never produces. Available columns: {sorted(produced_columns)}. "
                "A missing column must be rejected before flagging begins."
            )

    def test_every_name_a_rule_reads_resolves_without_an_ru_fasta(
        self, flagging_rules, produced_columns_without_ru_fasta
    ):
        """No missing column may silently disable a rule on the no-FASTA branch."""
        for name, configured in flagging_rules.items():
            unresolved = configured_columns(configured) - produced_columns_without_ru_fasta
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
        row = {"RU": "7", "Variant": "I3_7_A_LEN1", "NumberOfSupportingReads": 42}

        assert flag_fires("Repeat_Unit_7", flagging_rules["Repeat_Unit_7"], row) is True

    def test_the_rule_does_not_fire_for_another_repeat_unit(self, flagging_rules):
        row = {"RU": "2", "Variant": "I3_2_A_LEN1", "NumberOfSupportingReads": 42}

        assert flag_fires("Repeat_Unit_7", flagging_rules["Repeat_Unit_7"], row) is False

    def test_a_compound_call_spanning_two_units_does_not_fire(self, flagging_rules):
        row = {"RU": "7,2", "Variant": "I3_7_A_LEN1&I9_2_A_LEN1", "NumberOfSupportingReads": 42}

        assert flag_fires("Repeat_Unit_7", flagging_rules["Repeat_Unit_7"], row) is False

    def test_integer_repeat_unit_seven_is_rejected(self, flagging_rules):
        with pytest.raises(ValueError, match="compatible families"):
            flag_fires("Repeat_Unit_7", flagging_rules["Repeat_Unit_7"], {"RU": 7})

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

    def test_the_historic_integer_expression_is_not_a_supported_migration(self, flagging_rules):
        """Only the exact last-release string is eligible for migration."""
        rows = self.probe_rows()
        assert len(rows) >= 100, f"vacuity guard: only {len(rows)} probe rows"

        with pytest.raises(ValueError, match="unsupported legacy expression"):
            add_flags(rows, {"Repeat_Unit_7": self.HISTORIC_EXPRESSION})

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
            for name, configured in flagging_rules.items()
            if add_flags(frame, {name: configured})["Flag"].eq(name).any()
        }

        assert fired == set(flagging_rules), f"these rules can never fire: {sorted(set(flagging_rules) - fired)}"


class TestTheOtherRules:
    def test_low_coverage_fires_below_ten_supporting_reads(self, flagging_rules):
        rule = flagging_rules["Low_Coverage"]

        assert flag_fires("Low_Coverage", rule, {"NumberOfSupportingReads": 9}) is True
        assert flag_fires("Low_Coverage", rule, {"NumberOfSupportingReads": 10}) is False

    def test_polymorphic_call_fires_for_a_listed_variant(self, flagging_rules):
        rule = flagging_rules["Polymorphic_Call"]

        assert flag_fires("Polymorphic_Call", rule, {"Variant": "I10_2_A_LEN1"}) is True
        assert flag_fires("Polymorphic_Call", rule, {"Variant": "I22_2_G_LEN1"}) is False

    def test_one_character_polymorphism_mutation_does_not_fire(self, flagging_rules):
        rule = flagging_rules["Polymorphic_Call"]

        assert flag_fires("Polymorphic_Call", rule, {"Variant": "I10_2_A_LEN2"}) is False

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
            {
                "flagging_rules": {
                    "Custom_Flag": {
                        "all": [
                            {
                                "left": {"literal": 1},
                                "operator": "lt",
                                "right": {"column": "NumberOfSupportingReads"},
                            }
                        ]
                    }
                }
            },
        )
        source = tmp_path / "output_adVNTR.vcf"
        source.write_text(ADVNTR_HEADER + "25561\tI22_2_G_LEN1\t42\t153.98\t0.0001\n")

        advntr.process_advntr_output(str(source), str(tmp_path), "output", config=ru_config)

        result = pd.read_csv(tmp_path / "output_adVNTR_result.tsv", sep="\t", dtype=str)
        assert result.iloc[0]["Flag"] == "Custom_Flag"

    def test_patching_the_derived_settings_global_does_not_change_the_flags(self, tmp_path, monkeypatch, ru_config):
        monkeypatch.setattr(
            advntr,
            "advntr_settings",
            {"flagging_rules": {"Custom_Flag": {"all": [{"invalid": True}]}}},
        )
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
