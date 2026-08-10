from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

_PORCELAIN_V1_STATUS_PAIRS = frozenset(
    {
        b" M",
        b" T",
        b" A",
        b" D",
        b" R",
        b" C",
        b"M ",
        b"MM",
        b"MT",
        b"MD",
        b"T ",
        b"TM",
        b"TT",
        b"TD",
        b"A ",
        b"AM",
        b"AT",
        b"AD",
        b"D ",
        b"R ",
        b"RM",
        b"RT",
        b"RD",
        b"C ",
        b"CM",
        b"CT",
        b"CD",
        b"DD",
        b"AU",
        b"UD",
        b"UA",
        b"DU",
        b"AA",
        b"UU",
        b"??",
        b"!!",
    }
)


@dataclass(frozen=True)
class OverlayChange:
    """One file operation needed to reproduce a working-tree change."""

    action: Literal["copy", "delete"]
    path: str


def confined_path(root: Path, relative: str, *, must_exist: bool) -> Path:
    """Resolve a repository-relative path without allowing a worktree escape.

    Args:
        root: Root of the workspace containing the path.
        relative: Workspace-relative path to validate.
        must_exist: Whether the path must currently exist, including symlinks.

    Returns:
        The validated lexical path beneath the resolved workspace root.

    Raises:
        ValueError: If the name is unsafe, escapes the root, or must exist but does not.
    """
    candidate = Path(relative)
    if (
        relative in {"", "."}
        or candidate.is_absolute()
        or candidate.as_posix() != relative
        or ".." in candidate.parts
        or ".git" in candidate.parts
    ):
        raise ValueError(f"unsafe workspace path: {relative}")
    resolved_root = root.resolve()
    lexical = resolved_root / candidate
    try:
        resolved_parent = lexical.parent.resolve(strict=must_exist)
    except FileNotFoundError:
        raise ValueError(f"workspace path does not exist: {relative}") from None
    if not resolved_parent.is_relative_to(resolved_root):
        raise ValueError(f"workspace path escapes workspace root: {relative}")
    if ".git" in resolved_parent.relative_to(resolved_root).parts:
        raise ValueError(f"unsafe workspace path: {relative}")
    if must_exist and not os.path.lexists(lexical):
        raise ValueError(f"workspace path does not exist: {relative}")
    if lexical.is_symlink():
        try:
            resolved = lexical.resolve(strict=True)
        except FileNotFoundError:
            raise ValueError(f"workspace path does not exist: {relative}") from None
        if not resolved.is_relative_to(resolved_root):
            raise ValueError(f"workspace path escapes workspace root: {relative}")
        if ".git" in resolved.relative_to(resolved_root).parts:
            raise ValueError(f"unsafe workspace path: {relative}")
    return lexical


def _affected_paths_from_porcelain_z(payload: bytes) -> tuple[str, ...]:
    """Decode paths affected by porcelain v1 ``-z`` records.

    Rename and copy records contribute both their current path and their original path.
    This function does not inspect the filesystem.

    Args:
        payload: Raw output from ``git status --porcelain=v1 -z``.

    Returns:
        Affected paths in record order, including rename and copy origins.

    Raises:
        ValueError: If the payload is malformed, ignored, or has ambiguous rename metadata.
    """
    if not payload:
        return ()
    if not payload.endswith(b"\0"):
        raise ValueError("malformed porcelain payload: missing NUL terminator")

    records = payload[:-1].split(b"\0")
    affected: list[str] = []
    index = 0
    while index < len(records):
        record = records[index]
        if len(record) < 4 or record[2:3] != b" " or not record[3:]:
            raise ValueError("malformed porcelain record")
        status = record[:2]
        if status == b"!!":
            raise ValueError("ignored porcelain record is not allowed")

        rename_fields = sum(field in b"RC" for field in status)
        if rename_fields > 1:
            raise ValueError("conflicting rename encoding")
        if status not in _PORCELAIN_V1_STATUS_PAIRS:
            raise ValueError(f"invalid porcelain status: {os.fsdecode(status)}")

        affected.append(os.fsdecode(record[3:]))
        if rename_fields == 1:
            index += 1
            if index >= len(records) or not records[index]:
                raise ValueError("rename record is missing original path")
            affected.append(os.fsdecode(records[index]))
        index += 1
    return tuple(affected)


def parse_porcelain_z(payload: bytes, real_root: Path) -> tuple[OverlayChange, ...]:
    """Resolve porcelain records into a safe snapshot of current filesystem state.

    Args:
        payload: Raw output from ``git status --porcelain=v1 -z``.
        real_root: Root of the working tree whose current files decide each action.

    Returns:
        One path-sorted copy or delete action per affected path.

    Raises:
        ValueError: If any porcelain record or affected path is unsafe.
    """
    affected_paths = _affected_paths_from_porcelain_z(payload)
    validated_paths = {path: confined_path(real_root, path, must_exist=False) for path in dict.fromkeys(affected_paths)}
    changes = tuple(
        OverlayChange("copy" if os.path.lexists(validated) else "delete", path)
        for path, validated in validated_paths.items()
    )
    return tuple(sorted(changes, key=lambda change: change.path))
