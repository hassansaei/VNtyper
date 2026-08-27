"""Tests for filesystem-free adVNTR state annotation."""

import pandas as pd
import pytest

from vntyper.modules.advntr.advntr_variant_annotations import derive_ru_and_pos

pytestmark = pytest.mark.unit


@pytest.mark.parametrize(
    ("variant", "expected_ru", "expected_pos"),
    [
        ("I3_7_A_LEN1", "7", "3"),
        ("I22_2_G_LEN1", "2", "22"),
        ("D3_2", "2", "3"),
        ("D58_2&D59_2", "2,2", "58,59"),
        ("I3_7_A_LEN1&I9_2_A_LEN1", "7,2", "3,9"),
        ("NOT_A_VARIANT", ".", "."),
        ("I3_7_A_LEN1&NOT_A_VARIANT", "7,.", "3,."),
        ("", ".", "."),
    ],
)
def test_derive_ru_and_pos_reads_each_state_part(
    variant: str,
    expected_ru: str,
    expected_pos: str,
) -> None:
    """Each state part contributes one ordered RU and position annotation."""
    ru, pos = derive_ru_and_pos([variant])

    assert (ru, pos) == ([expected_ru], [expected_pos])


def test_derive_ru_and_pos_accepts_a_one_pass_iterable() -> None:
    """The public iterable contract does not require a pandas Series or a sequence."""
    variants = (variant for variant in ["I3_7_A_LEN1", "D58_2&D59_2"])

    assert derive_ru_and_pos(variants) == (["7", "2,2"], ["3", "58,59"])


@pytest.mark.parametrize("variant", [25561, None, pd.NA, False])
def test_derive_ru_and_pos_maps_non_strings_to_dots(variant: object) -> None:
    """Malformed DataFrame cells cannot abort unconditional annotation."""
    assert derive_ru_and_pos([variant]) == (["."], ["."])
