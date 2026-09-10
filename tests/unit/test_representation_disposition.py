"""Unit tests for pure representation disposition decisions."""

import pytest

from vntyper.scripts.representation_disposition import (
    DISPOSITION_REPRESENTATION_LIMITED,
    FLAG_ALLELE_UNREPRESENTABLE,
    NOTE_REPRESENTATION_LIMITED,
    format_masthead_representation_limited,
    format_representation_limited_name,
    is_unrepresentable_molecular_class,
)

pytestmark = pytest.mark.unit


def test_is_unrepresentable_molecular_class():
    assert is_unrepresentable_molecular_class("delins") is True
    assert is_unrepresentable_molecular_class("insertion", ref_span=1, inserted=2) is True
    assert is_unrepresentable_molecular_class("deletion", flags=[FLAG_ALLELE_UNREPRESENTABLE]) is True
    assert is_unrepresentable_molecular_class("insertion", ref_span=0, inserted=1) is False
    assert is_unrepresentable_molecular_class("deletion", ref_span=1, inserted=0) is False
    assert is_unrepresentable_molecular_class("duplication", ref_span=0, inserted=1) is False


def test_format_representation_limited_name():
    assert format_representation_limited_name(1) == "frameshift +1, representation-limited"
    assert format_representation_limited_name(-2) == "frameshift -2, representation-limited"
    assert format_representation_limited_name(3) == "in-frame +3, representation-limited"
    assert format_representation_limited_name(-6) == "in-frame -6, representation-limited"
    assert format_representation_limited_name(0) == DISPOSITION_REPRESENTATION_LIMITED


def test_format_masthead_representation_limited():
    assert (
        format_masthead_representation_limited(1)
        == "Frameshift detected (+1 bp, positional naming withheld: representation-limited)"
    )
    assert (
        format_masthead_representation_limited(-1)
        == "Frameshift detected (-1 bp, positional naming withheld: representation-limited)"
    )
    assert (
        format_masthead_representation_limited(3)
        == "In-frame variant detected (+3 bp, positional naming withheld: representation-limited)"
    )
    assert (
        format_masthead_representation_limited(-3)
        == "In-frame variant detected (-3 bp, positional naming withheld: representation-limited)"
    )
    assert (
        format_masthead_representation_limited(0)
        == "Variant detected (positional naming withheld: representation-limited)"
    )
    assert NOTE_REPRESENTATION_LIMITED != ""
