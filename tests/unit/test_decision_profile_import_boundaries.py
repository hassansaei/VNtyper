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
    REPOSITORY_ROOT / "vntyper/scripts/pipeline.py": {
        "advntr_config",
        "advntr_settings",
        "kestrel_config",
        "shark_config",
        "shark_settings",
    },
    REPOSITORY_ROOT / "vntyper/scripts/pipeline_kestrel.py": {"kestrel_config"},
}
NOMENCLATURE_CONSUMERS = (
    REPOSITORY_ROOT / "vntyper/scripts/kestrel_genotyping.py",
    REPOSITORY_ROOT / "vntyper/modules/advntr/advntr_genotyping.py",
    REPOSITORY_ROOT / "vntyper/scripts/nomenclature_annotate.py",
    REPOSITORY_ROOT / "vntyper/scripts/nomenclature_bam_adapter.py",
)


def _module_assignments(tree: ast.Module) -> set[str]:
    """Return names bound at module scope, including imports and nested blocks."""

    class ModuleBindingVisitor(ast.NodeVisitor):
        def __init__(self) -> None:
            self.assigned: set[str] = set()

        def visit_Name(self, node: ast.Name) -> None:  # noqa: N802 - ast visitor API
            if isinstance(node.ctx, ast.Store):
                self.assigned.add(node.id)

        def visit_Import(self, node: ast.Import) -> None:  # noqa: N802 - ast visitor API
            self.assigned.update(alias.asname or alias.name.split(".")[0] for alias in node.names)

        def visit_ImportFrom(self, node: ast.ImportFrom) -> None:  # noqa: N802 - ast visitor API
            self.assigned.update(alias.asname or alias.name for alias in node.names)

        def visit_FunctionDef(self, node: ast.FunctionDef) -> None:  # noqa: N802 - ast visitor API
            self.assigned.add(node.name)

        def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:  # noqa: N802 - ast visitor API
            self.assigned.add(node.name)

        def visit_ClassDef(self, node: ast.ClassDef) -> None:  # noqa: N802 - ast visitor API
            self.assigned.add(node.name)

    visitor = ModuleBindingVisitor()
    visitor.visit(tree)
    return visitor.assigned


def _called_name(node: ast.Call) -> str | None:
    """Return the terminal name of a direct or module-qualified call."""
    if isinstance(node.func, ast.Name):
        return node.func.id
    if isinstance(node.func, ast.Attribute):
        return node.func.attr
    return None


def test_module_binding_scan_covers_imports_and_nested_module_blocks_but_not_locals() -> None:
    tree = ast.parse("import package as imported\nif enabled:\n    nested = 1\ndef function():\n    local = 2\n")

    assert _module_assignments(tree) == {"function", "imported", "nested"}


def test_call_scan_covers_direct_and_module_qualified_calls() -> None:
    direct, qualified = (
        node for node in ast.walk(ast.parse("load_config()\nmodule.load_config()\n")) if isinstance(node, ast.Call)
    )

    assert _called_name(direct) == _called_name(qualified) == "load_config"


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
            isinstance(node, ast.Call) and _called_name(node) == "load_nomenclature_config" for node in ast.walk(tree)
        ):
            violations.append(str(path.relative_to(REPOSITORY_ROOT)))

    assert violations == []


def test_shark_consumers_do_not_reload_the_runtime_sidecar_as_decision_policy() -> None:
    """SHARK's empty profile component must be resolved, not inferred from its path sidecar."""
    consumers = (
        REPOSITORY_ROOT / "vntyper/scripts/pipeline.py",
        REPOSITORY_ROOT / "vntyper/modules/shark/shark_filtering.py",
    )
    violations: list[str] = []
    for path in consumers:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        if any(isinstance(node, ast.Call) and _called_name(node) == "load_shark_config" for node in ast.walk(tree)):
            violations.append(str(path.relative_to(REPOSITORY_ROOT)))

    assert violations == []


def test_shark_has_no_independent_decision_constant_or_literal_comparison() -> None:
    """Any new SHARK decision must be inventoried instead of hard-coded beside the empty component."""
    path = REPOSITORY_ROOT / "vntyper/modules/shark/shark_filtering.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    uppercase_bindings = sorted(name for name in _module_assignments(tree) if name.isupper())
    comparisons = {ast.unparse(node) for node in ast.walk(tree) if isinstance(node, ast.Compare)}

    assert uppercase_bindings == []
    assert comparisons == {
        "config_path is None",
        "key in settings",
        "kept_reads_r1 != kept_reads_r2",
        "legacy_key in settings",
        "resolved is not None",
    }


def test_shark_sidecar_contains_runtime_references_only() -> None:
    """An empty SHARK profile remains complete while its sidecar has paths and runtime parameters only."""
    sidecar = json.loads((REPOSITORY_ROOT / "vntyper/modules/shark/shark_config.json").read_text(encoding="utf-8"))

    assert set(sidecar) == {"shark_settings"}
    settings = sidecar["shark_settings"]
    assert settings
    fasta_keys = {key for key in settings if key.startswith("muc1_region_fasta_")}
    assert fasta_keys == {"muc1_region_fasta_hg19", "muc1_region_fasta_hg38"}
    assert all(isinstance(settings[key], str) and settings[key] for key in fasta_keys)
    assert settings["kmer_size"] == 17
    assert settings["confidence"] == 0.6
    assert isinstance(settings.get("_comment_search_parameters"), str)
