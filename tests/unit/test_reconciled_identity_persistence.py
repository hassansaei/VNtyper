"""Closed codec for typed whole-locus reconciliation selections."""

from __future__ import annotations

import pytest

from vntyper.scripts.molecular_identity import make_coding_edit, make_molecular_identity
from vntyper.scripts.reconciled_identity_persistence import (
    decode_reconciled_identity,
    encode_reconciled_identity,
)

pytestmark = pytest.mark.unit


def test_reconciled_identity_codec_round_trips_selection_and_abstention() -> None:
    identity = make_molecular_identity((make_coding_edit(60, 59, "", "C"),))

    assert decode_reconciled_identity(encode_reconciled_identity(identity)) == identity
    assert decode_reconciled_identity(encode_reconciled_identity(None)) is None


@pytest.mark.parametrize("value", [None, 1, "not-an-identity"])
def test_reconciled_identity_decoder_rejects_noncanonical_cells(value: object) -> None:
    with pytest.raises(ValueError, match="identity|string|canonical"):
        decode_reconciled_identity(value)


def test_reconciled_identity_encoder_rejects_untyped_values() -> None:
    with pytest.raises(ValueError, match="MolecularIdentity"):
        encode_reconciled_identity("MUC1-X")  # type: ignore[arg-type]
