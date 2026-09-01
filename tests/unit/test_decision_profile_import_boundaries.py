"""Structural boundaries for per-run decision-profile ownership."""

from __future__ import annotations

import ast
import json
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

REPOSITORY_ROOT = Path(__file__).parents[2]
DECISION_MODULES = {
    REPOSITORY_ROOT / "vntyper/scripts/kestrel_genotyping.py": {"kestrel_config"},
    REPOSITORY_ROOT / "vntyper/modules/advntr/advntr_genotyping.py": {
        "advntr_config",
        "advntr_settings",
    },
    REPOSITORY_ROOT / "vntyper/modules/shark/shark_filtering.py": {
        "shark_config",
        "shark_settings",
    },
}
NOMENCLATURE_CONSUMERS = (
    REPOSITORY_ROOT / "vntyper/scripts/kestrel_genotyping.py",
    REPOSITORY_ROOT / "vntyper/modules/advntr/advntr_genotyping.py",
    REPOSITORY_ROOT / "vntyper/scripts/nomenclature_annotate.py",
    REPOSITORY_ROOT / "vntyper/scripts/nomenclature_bam_adapter.py",
)


def _module_assignments(tree: ast.Module) -> set[str]:
    """Return names assigned directly in a module body."""
    assigned: set[str] = set()
    for statement in tree.body:
        targets: list[ast.expr] = []
        if isinstance(statement, ast.Assign):
            targets = statement.targets
        elif isinstance(statement, ast.AnnAssign):
            targets = [statement.target]
        for target in targets:
            if isinstance(target, ast.Name):
                assigned.add(target.id)
    return assigned


def test_stage_modules_do_not_own_packaged_decision_globals() -> None:
    """A custom run must not compete with state loaded when a module was imported."""
    violations: list[str] = []
    for path, forbidden_names in DECISION_MODULES.items():
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        present = sorted(_module_assignments(tree) & forbidden_names)
        if present:
            violations.append(f"{path.relative_to(REPOSITORY_ROOT)}: {', '.join(present)}")

    assert violations == []


def test_nomenclature_consumers_do_not_reload_the_legacy_sidecar() -> None:
    """Resolved run policy must not compete with a second filesystem read."""
    violations: list[str] = []
    for path in NOMENCLATURE_CONSUMERS:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        if any(
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "load_nomenclature_config"
            for node in ast.walk(tree)
        ):
            violations.append(str(path.relative_to(REPOSITORY_ROOT)))

    assert violations == []


def test_shark_sidecar_contains_runtime_references_only() -> None:
    """An empty SHARK profile remains complete while its sidecar has paths only."""
    sidecar = json.loads((REPOSITORY_ROOT / "vntyper/modules/shark/shark_config.json").read_text(encoding="utf-8"))

    assert set(sidecar) == {"shark_settings"}
    assert sidecar["shark_settings"]
    assert all(key.startswith("muc1_region_fasta_") for key in sidecar["shark_settings"])
    assert all(isinstance(value, str) and value for value in sidecar["shark_settings"].values())
