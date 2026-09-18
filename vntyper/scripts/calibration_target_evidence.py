"""Outcome-free extraction of target-v2 calibration evidence metadata."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import NoReturn

from vntyper.scripts.calibration_artifact_io import write_checksums, write_json
from vntyper.scripts.calibration_role_source import RoleSource, decode_role_source, role_source_document
from vntyper.scripts.calibration_secure_io import read_regular_path
from vntyper.scripts.calibration_target_contract import (
    LengthBaselinePlan,
    TargetStudy,
    decode_target_study,
    target_study_document,
)
from vntyper.scripts.calibration_target_runs import TargetRuns, decode_target_runs, target_runs_document
from vntyper.scripts.canonical_json import load_strict_json_object
from vntyper.scripts.length_annotation import LengthAnnotation, decode_length_annotation, encode_length_annotation

logger = logging.getLogger(__name__)

_ROLES = ("training", "policy-selection", "validation", "locked-heldout")


def _fail(message: str) -> NoReturn:
    logger.error(message)
    raise ValueError(message)


def _object(path: Path, label: str) -> dict[str, object]:
    try:
        return load_strict_json_object(read_regular_path(path))
    except ValueError as error:
        raise ValueError(f"{label} is missing, unreadable, or invalid") from error


def _sources(root: Path, study: TargetStudy, runs: TargetRuns) -> tuple[RoleSource, ...]:
    roles_root = root / "roles"
    if (
        not isinstance(root, Path)
        or not roles_root.is_dir()
        or roles_root.is_symlink()
        or {path.name for path in roles_root.iterdir()} != set(_ROLES)
    ):
        _fail("target evidence sources require the exact four-role directory inventory")
    sources = []
    for role in _ROLES:
        role_root = roles_root / role
        if not role_root.is_dir() or role_root.is_symlink() or not (role_root / "source.json").is_file():
            _fail("target evidence role source inventory differs")
        source = decode_role_source(
            _object(role_root / "source.json", f"target evidence {role} source"),
            study=study,
            runs=runs,
            expected_role=role,
        )
        local_truth = source.truth_asset.path.parent == role_root.absolute()
        expected_names = {"source.json", source.truth_asset.path.name} if local_truth else {"source.json"}
        if {path.name for path in role_root.iterdir()} != expected_names:
            _fail("target evidence role source inventory differs")
        sources.append(source)
    return tuple(sources)


def _annotation(study: TargetStudy, path: Path | None) -> LengthAnnotation | None:
    if study.target == "callers":
        if path is not None:
            _fail("caller target evidence forbids a length annotation")
        return None
    if path is None or not isinstance(study.baseline, LengthBaselinePlan):
        _fail("length target evidence requires an explicit annotation")
    annotation = decode_length_annotation(_object(path, "length target evidence annotation"))
    if annotation.sha256 != study.baseline.annotation_sha256:
        _fail("length target evidence annotation differs from the study baseline")
    return annotation


def extract_target_evidence(
    study_path: Path,
    runs_path: Path,
    sources_root: Path,
    output: Path,
    *,
    expected_target: str,
    length_annotation_path: Path | None = None,
) -> bool:
    """Validate and copy target metadata without opening sealed outcomes.

    Args:
        study_path: Frozen target study document.
        runs_path: Complete target run commitments.
        sources_root: Directory containing ``roles/<role>/source.json``.
        output: Empty staging directory supplied by atomic CLI publication.
        expected_target: Explicit CLI target, exactly callers or length.
        length_annotation_path: Required for length and forbidden for callers.

    Returns:
        True after writing a fit-ready canonical evidence tree.

    Raises:
        ValueError: If target, inventory, source, run, or annotation bindings differ.
    """
    if expected_target not in {"callers", "length"}:
        _fail("target evidence expected target must be callers or length")
    if any(not isinstance(path, Path) for path in (study_path, runs_path, sources_root, output)):
        _fail("target evidence paths must be Path values")
    study = decode_target_study(_object(study_path, "target evidence study"))
    runs = decode_target_runs(_object(runs_path, "target evidence runs"))
    if study.target != expected_target or runs.target != expected_target:
        _fail("target evidence study or runs differ from the expected target")
    annotation = _annotation(study, length_annotation_path)
    sources = _sources(sources_root, study, runs)

    write_json(output / "study.json", target_study_document(study))
    write_json(output / "runs.json", target_runs_document(runs))
    if annotation is not None:
        write_json(output / "annotation.json", encode_length_annotation(annotation))
    for source in sources:
        write_json(
            output / "roles" / source.role / "source.json",
            role_source_document(source, study=study, runs=runs),
        )
    write_checksums(output)
    return True
