# tests/unit/test_advntr_polymorphic_calls.py

"""The adVNTR ``Polymorphic_Call`` denylist: reachable, satisfiable, provenanced (#267).

Flagging runs *after* the pathogenic-frame filter, so an entry whose signed net indel delta
is not ``1 (mod 3)`` is removed before ``add_flags`` is ever reached and can never fire.
Seven of the thirty-two shipped entries were in that state, and one was listed twice.

Neither failure was loud. ``evaluate_condition`` downgrades a ``NameError`` to a warning and
returns False, so a broken rule disables itself in silence -- which is how
``Poylmorhic_Call`` shipped misspelled until ``742b872`` and how ``Repeat_Unit_7`` shipped
comparing a string column against an int until ``52f822e``. A dead denylist entry is the
same class of defect one level down, and this module retires it.

Reachability here is decided by running the **production** filter functions, not by
restating ``delta % 3 == 1``. A change to the filter that stranded the list would then fail
this file rather than agree with it.
"""

import ast
import json
from collections import Counter
from pathlib import Path

import pandas as pd
import pytest

from vntyper.modules.advntr import advntr_genotyping as advntr
from vntyper.scripts.flagging import evaluate_condition

pytestmark = pytest.mark.unit

#: Entries the owner confirmed as established artifacts on #267 (2026-08-24):
#: "I would recommend keeping D58_2&D59_2 and anything reported in RU7 flagged".
CONFIRMED = {"D58_2&D59_2"}

#: How many distinct states the list carries after the #267 cleanup: 32 shipped, minus 7
#: that the pathogenic-frame filter removes before flagging, minus 1 duplicate. Changing
#: this is the deliberate act that editing the denylist requires; it is not maintenance.
EXPECTED_ENTRY_COUNT = 24

#: The three flags this module's configuration is expected to emit.
EXPECTED_FLAG_NAMES = {"Low_Coverage", "Repeat_Unit_7", "Polymorphic_Call"}

MODULE_DIR = Path(advntr.__file__).parent


@pytest.fixture(scope="module")
def config() -> dict:
    return json.loads((MODULE_DIR / "advntr_config.json").read_text())


@pytest.fixture(scope="module")
def calibration() -> dict:
    return json.loads((MODULE_DIR / "advntr_calibration.json").read_text())


@pytest.fixture(scope="module")
def states(config) -> list[str]:
    """The state strings the shipped ``Polymorphic_Call`` rule tests membership against."""
    expression = config["flagging_rules"]["Polymorphic_Call"]
    node = ast.parse(expression, mode="eval").body
    assert isinstance(node, ast.Compare), "Polymorphic_Call is expected to be a membership test"
    return ast.literal_eval(node.comparators[0])


def raw_frame(states: list[str]) -> pd.DataFrame:
    """The shape both filter arms consume, straight out of adVNTR's VCF."""
    return pd.DataFrame(
        {
            "VID": ["25561"] * len(states),
            "State": list(states),
            "NumberOfSupportingReads": [50] * len(states),
            "MeanCoverage": [100.0] * len(states),
            "Pvalue\n": [0.0001] * len(states),
        }
    )


def survivors(states: list[str]) -> set[str]:
    """The states the production pathogenic-frame filter keeps."""
    frame = raw_frame(states)
    return set(advntr.advntr_processing_ins(frame.copy())["Variant"]) | set(
        advntr.advntr_processing_del(frame.copy())["Variant"]
    )


class TestTheListItself:
    def test_the_scan_found_the_list(self, states):
        """Guard the guard: an empty list makes every assertion below vacuous."""
        assert len(states) == EXPECTED_ENTRY_COUNT

    def test_no_entry_is_listed_twice(self, states):
        repeated = sorted(state for state, count in Counter(states).items() if count > 1)

        assert not repeated, f"duplicated denylist entries: {repeated}"

    def test_every_entry_is_a_non_empty_string(self, states):
        assert all(isinstance(state, str) and state.strip() for state in states)

    def test_the_flag_names_are_the_expected_set(self, config):
        assert set(config["flagging_rules"]) == EXPECTED_FLAG_NAMES


class TestReachability:
    def test_every_entry_survives_the_pathogenic_frame_filter(self, states):
        """The #267 defect. Run the production arms; anything they drop can never be
        flagged, so listing it is a claim the code cannot honour."""
        unreachable = sorted(set(states) - survivors(states))

        assert not unreachable, (
            f"{len(unreachable)} denylist entries are removed by the pathogenic-frame filter "
            f"before add_flags runs, so they can never fire: {unreachable}"
        )

    def test_the_filter_really_does_drop_something(self):
        """Guard the guard: if both arms passed everything, the test above proves nothing.
        Delta = +2 is a genuine frameshift, in the non-pathogenic frame."""
        assert not survivors(["I17_2_G_LEN2"])

    @pytest.mark.parametrize(
        "state",
        [
            "D17_2&I17_2_G_LEN3",
            "D46_2&D47_2&I47_2_A_LEN4",
            "D49_2&I49_2_A_LEN12",
            "I22_2_C_LEN38",
            "I26_7_A_LEN24",
            "I60_2_A_LEN2",
            "I9_2_C_LEN5",
        ],
    )
    def test_the_seven_removed_entries_really_were_unreachable(self, state, states):
        """Pins *why* the cleanup is behaviour-preserving rather than asserting that it is:
        each removed entry is one the filter drops, so no Flag value can have depended on it."""
        assert state not in states, f"{state} was removed by the #267 cleanup"
        assert not survivors([state]), f"{state} would have been reachable after all"


class TestSatisfiability:
    """Every configured rule must be able to be True for some row the parser can produce.

    A rule that parses cleanly and resolves every name can still be dead: ``Repeat_Unit_7``
    was ``RU == 7`` against a column that is always a string.
    """

    def probes(self, states: list[str]) -> list[pd.Series]:
        return [
            pd.Series({"Variant": states[0], "RU": "2", "NumberOfSupportingReads": 50}),
            pd.Series({"Variant": "I3_7_A_LEN1", "RU": "7", "NumberOfSupportingReads": 50}),
            pd.Series({"Variant": "I1_2_A_LEN1", "RU": "2", "NumberOfSupportingReads": 3}),
        ]

    def test_polymorphic_call_fires_for_every_listed_state(self, config, states):
        rule = config["flagging_rules"]["Polymorphic_Call"]

        for state in states:
            row = pd.Series({"Variant": state, "RU": "2", "NumberOfSupportingReads": 50})
            assert evaluate_condition(row, rule) is True, state

    def test_polymorphic_call_does_not_fire_for_an_unlisted_state(self, config):
        row = pd.Series({"Variant": "I1_2_A_LEN1", "RU": "2", "NumberOfSupportingReads": 50})

        assert evaluate_condition(row, config["flagging_rules"]["Polymorphic_Call"]) is False

    def test_polymorphic_call_does_not_fire_for_a_removed_state(self, config):
        """The cleanup's visible consequence, stated as behaviour rather than as a count."""
        row = pd.Series({"Variant": "I22_2_C_LEN38", "RU": "2", "NumberOfSupportingReads": 50})

        assert evaluate_condition(row, config["flagging_rules"]["Polymorphic_Call"]) is False

    def test_every_rule_fires_for_at_least_one_row(self, config, states):
        probes = self.probes(states)

        for name, expression in config["flagging_rules"].items():
            assert any(evaluate_condition(row, expression) for row in probes), (
                f"flagging rule {name!r} is never True for any probe row, so it is dead. "
                "See #267 and the Repeat_Unit_7 history."
            )

    def test_every_rule_is_false_for_at_least_one_row(self, config, states):
        """The other half: a rule that is always True is not a flag, it is a constant."""
        probes = self.probes(states)

        for name, expression in config["flagging_rules"].items():
            assert not all(evaluate_condition(row, expression) for row in probes), (
                f"flagging rule {name!r} is True for every probe row, so it distinguishes nothing."
            )


class TestProvenance:
    def test_the_calibration_file_covers_exactly_the_live_list(self, calibration, states):
        documented = {entry["state"] for entry in calibration["polymorphic_calls"]["entries"]}

        assert documented == set(states)

    def test_it_documents_every_entry_once(self, calibration):
        documented = [entry["state"] for entry in calibration["polymorphic_calls"]["entries"]]

        assert len(documented) == len(set(documented)) == EXPECTED_ENTRY_COUNT

    def test_every_entry_declares_a_known_status(self, calibration):
        allowed = set(calibration["polymorphic_calls"]["statuses"])

        assert allowed, "the guard is vacuous with no statuses declared"
        for entry in calibration["polymorphic_calls"]["entries"]:
            assert entry["status"] in allowed, entry

    def test_the_owner_confirmed_entries_are_recorded_as_confirmed(self, calibration):
        by_state = {entry["state"]: entry for entry in calibration["polymorphic_calls"]["entries"]}

        for state in CONFIRMED:
            assert by_state[state]["status"] == "confirmed_artifact"

    def test_the_rest_are_recorded_as_awaiting_revalidation(self, calibration):
        """The owner asked for these to be re-measured against the renome cohort and
        decided case by case; recording that is the point of the file."""
        pending = [
            entry["state"]
            for entry in calibration["polymorphic_calls"]["entries"]
            if entry["status"] == "pending_renome_revalidation"
        ]

        assert len(pending) == EXPECTED_ENTRY_COUNT - len(CONFIRMED)

    def test_the_recorded_delta_is_the_one_the_code_computes(self, calibration):
        for entry in calibration["polymorphic_calls"]["entries"]:
            state = entry["state"]
            computed = advntr.sum_insertion_lengths(state) - state.count("D")

            assert entry["net_indel_delta"] == computed, state

    def test_every_recorded_delta_is_in_the_pathogenic_frame(self, calibration):
        """The same claim as the reachability test, over the documented data, so a
        hand-edited calibration file cannot disagree with the list it documents."""
        for entry in calibration["polymorphic_calls"]["entries"]:
            assert entry["net_indel_delta"] % 3 == 1, entry

    def test_the_colliding_entries_carry_a_note(self, calibration):
        """#267's substance: three live entries have the shape of a pathogenic variant the
        simulation produces, and a State string cannot separate the two. Losing that from
        the record is how the list became unexplainable in the first place."""
        by_state = {entry["state"]: entry for entry in calibration["polymorphic_calls"]["entries"]}

        for state in ("I23_6_G_LEN1", "I21_2_T_LEN1", "D17_2&D18_2&D19_2&D20_2&D21_2"):
            assert by_state[state].get("note"), state

    def test_the_criterion_the_cohort_and_the_limitation_are_recorded(self, calibration):
        block = calibration["polymorphic_calls"]

        assert block["cohort"] == "renome"
        assert block["criterion"]
        assert block["source"]
        assert block["limitation"]
        assert block["observation_counts_retained"] is False

    def test_the_other_two_rules_are_documented_too(self, calibration, config):
        documented = set(calibration["other_rules"])
        rules = set(config["flagging_rules"]) - {"Polymorphic_Call"}

        assert documented == rules

    def test_production_code_never_reads_this_file(self):
        """It is documentation with a test behind it, exactly like `calibration.json`."""
        hits = sorted(
            str(path)
            for path in Path("vntyper").rglob("*.py")
            if "advntr_calibration" in path.read_text(encoding="utf-8")
        )

        assert not hits, f"advntr_calibration.json is read by production code: {hits}"
