# tests/unit/test_advntr_options.py

"""
Unit tests for adVNTR command options and capability validation
(``vntyper/modules/advntr/advntr_options.py``).
"""

import pytest

from vntyper.modules.advntr import advntr_options as opt

pytestmark = pytest.mark.unit


def test_managed_and_extra_options_are_disjoint():
    """An option cannot be both managed by VNtyper and passed via additional_commands."""
    overlap = opt.MANAGED_ADVNTR_OPTIONS.intersection(opt.ADVNTR_EXTRA_OPTIONS)
    assert not overlap, f"Managed and extra options overlap: {overlap}"


def test_every_managed_option_names_an_owner():
    """Every managed option must define an owner description for informative refusals."""
    assert set(opt.MANAGED_ADVNTR_OPTIONS) == set(opt.MANAGED_ADVNTR_OPTION_OWNERS)


def test_v22_options_are_part_of_extra_options():
    """All v2.2.0 options must be allowlisted in ADVNTR_EXTRA_OPTIONS."""
    expected = {
        "--prune-reverse",
        "--exact-frameshift-caller",
        "--frameshift-background",
        "--frameshift-calibration-out",
        "--rare-unit-coverage-guard",
    }
    assert expected == opt.ADVNTR_V22_OPTIONS
    assert opt.ADVNTR_V22_OPTIONS.issubset(opt.ADVNTR_EXTRA_OPTIONS)


@pytest.mark.parametrize(
    "token, expected",
    [
        ("--prune-reverse", "--prune-reverse"),
        ("--rare-unit-coverage-guard=0.15", "--rare-unit-coverage-guard"),
        ("--frameshift-background=/path/to/bg", "--frameshift-background"),
        ("-aln", "-aln"),
        ("-c", "-c"),
        ("0.15", None),
        ("/path/to/bg", None),
        ("-1", None),
        ("-1.5", None),
        ("plain_string", None),
    ],
)
def test_advntr_option_token(token, expected):
    assert opt.advntr_option_token(token) == expected


@pytest.mark.parametrize(
    "cmd",
    [
        "--prune-reverse",
        "--rare-unit-coverage-guard 0.15",
        "--rare-unit-coverage-guard=0.15",
        "--rare-unit-coverage-guard",
        "--frameshift-calibration-out /tmp/calib.jsonl",
        "--exact-frameshift-caller --frameshift-background /tmp/bg.json",
        "--prune-reverse --rare-unit-coverage-guard 0.15",
    ],
)
def test_v22_options_pass_validation(cmd):
    result = opt.resolve_additional_commands({"additional_commands": cmd}, advntr_version=(2, 2, 0))
    assert result == cmd


def test_exact_frameshift_caller_without_background_raises():
    with pytest.raises(ValueError, match="without '--frameshift-background'"):
        opt.resolve_additional_commands({"additional_commands": "--exact-frameshift-caller"})


def test_v22_option_refused_on_older_advntr():
    with pytest.raises(ValueError, match="requires adVNTR >= 2.2.0, but installed adVNTR version is 2.0.4"):
        opt.resolve_additional_commands(
            {"additional_commands": "--prune-reverse"},
            advntr_version=(2, 0, 4),
        )


def test_legacy_options_pass_on_older_advntr():
    assert (
        opt.resolve_additional_commands(
            {"additional_commands": "-aln --haploid"},
            advntr_version=(2, 0, 4),
        )
        == "-aln --haploid"
    )
