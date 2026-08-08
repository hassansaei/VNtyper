"""Select the BAMs that may become derived CRAM fixtures.

The test-data manifest is the normal fixture contract. Discovery remains useful for the
golden cohort, but an unrelated BAM under ``tests/data`` must not silently grow the
fixture set a developer creates with the ordinary command.
"""

from __future__ import annotations

import json
from pathlib import Path

KNOWN_DERIVED_FIXTURE_KINDS = frozenset({"single_end_bam"})


def _validate_derived_fixtures(payload: dict[object, object]) -> None:
    """Validate declarations owned by neighboring fixture derivation tasks.

    Task 9 owns construction of ``single_end_bam`` fixtures. CRAM selection still parses
    that declaration so a new kind cannot silently be left out of a future dispatcher.
    """
    declarations = payload.get("derived_fixtures", [])
    if not isinstance(declarations, list):
        raise ValueError("derived_fixtures must be a list")
    for declaration in declarations:
        if not isinstance(declaration, dict):
            raise ValueError("each derived fixture declaration must be an object")
        kind = declaration.get("kind")
        if kind not in KNOWN_DERIVED_FIXTURE_KINDS:
            raise ValueError(f"unsupported derived fixture kind {kind!r}")
        for field in ("name", "source_bam", "output_bam"):
            if not isinstance(declaration.get(field), str) or not declaration[field]:
                raise ValueError(f"derived fixture kind {kind!r} must declare non-empty {field!r}")


def declared_bam_paths(data_config: Path, data_root: Path) -> set[Path]:
    """Return the BAM paths declared in a test-data manifest.

    Args:
        data_config: JSON test-data manifest containing ``file_resources``.
        data_root: Root that corresponds to ``tests/data`` for this invocation.

    Returns:
        The absolute paths of declared BAMs beneath ``data_root``.
    """
    payload = json.loads(data_config.read_text())
    if not isinstance(payload, dict):
        raise ValueError("test-data config must be an object")
    _validate_derived_fixtures(payload)
    declared: set[Path] = set()
    for resource in payload.get("file_resources", []):
        filename = resource.get("filename")
        local_path = resource.get("local_path")
        if not isinstance(filename, str) or not filename.endswith(".bam") or not isinstance(local_path, str):
            continue
        location = Path(local_path)
        try:
            relative_directory = location.relative_to("tests/data")
        except ValueError:
            continue
        declared.add((data_root / relative_directory / filename).resolve())
    return declared


def select_source_bams(discovered: list[Path], *, data_config: Path, data_root: Path, include_all: bool) -> list[Path]:
    """Choose discovered BAMs according to the requested fixture scope.

    Args:
        discovered: All source BAMs found below ``data_root``.
        data_config: Test-data manifest declaring the normal fixture set.
        data_root: Root passed to discovery.
        include_all: Whether an intentional all-cohort derivation was requested.

    Returns:
        A deterministic list of discovered BAMs to derive.
    """
    declared = declared_bam_paths(data_config, data_root)
    if include_all:
        return discovered
    return [bam for bam in discovered if bam.resolve() in declared]
