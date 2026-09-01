"""Closed persistence codec for the typed whole-locus identity selection."""

from __future__ import annotations

from vntyper.scripts.molecular_identity import MolecularIdentity, parse_molecular_identity, serialize_molecular_identity

RECONCILED_IDENTITY_COLUMN = "__Reconciled_Molecular_Identity"


def encode_reconciled_identity(identity: MolecularIdentity | None) -> str:
    """Encode a typed reconciler selection for final caller artifacts.

    Args:
        identity: Selected canonical identity, or ``None`` after abstention.

    Returns:
        Canonical serialized identity or an empty abstention cell.

    Raises:
        ValueError: If ``identity`` is not a typed molecular identity or ``None``.
    """
    if identity is None:
        return ""
    if not isinstance(identity, MolecularIdentity):
        raise ValueError("reconciled identity persistence requires a MolecularIdentity or None")
    return serialize_molecular_identity(identity)


def decode_reconciled_identity(value: object) -> MolecularIdentity | None:
    """Decode one exact typed reconciler selection cell.

    Args:
        value: Final caller artifact cell.

    Returns:
        Parsed canonical identity, or ``None`` for the empty abstention cell.

    Raises:
        ValueError: If the cell is non-string or not canonically serialized.
    """
    if not isinstance(value, str):
        raise ValueError("reconciled identity persistence cell must be a string")
    if not value:
        return None
    identity = parse_molecular_identity(value)
    if serialize_molecular_identity(identity) != value:
        raise ValueError("reconciled identity persistence cell must use canonical serialization")
    return identity
