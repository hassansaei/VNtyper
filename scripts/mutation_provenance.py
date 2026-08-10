from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Protocol

from mutation_workspace_fs import RootCapability, git_capability_path, using_root_capability


class WorkspaceRoots(Protocol):
    """Root capabilities needed by the import-provenance probe."""

    @property
    def real_root(self) -> Path:
        """Return the immutable real checkout root."""
        ...

    @property
    def sweep_root(self) -> Path:
        """Return the immutable disposable checkout root."""
        ...

    @property
    def root_capability(self) -> RootCapability | None:
        """Return the lifecycle-pinned disposable root when available."""
        ...


def verify_import_provenance(workspace: WorkspaceRoots, modules: tuple[str, ...]) -> tuple[Path, ...]:
    """Prove a child imports VNtyper and selected modules from the sweep root.

    Args:
        workspace: Real and disposable roots that define the isolation boundary.
        modules: Repository-relative Python module paths selected for mutation.

    Returns:
        Resolved package and module paths printed by the child probe.

    Raises:
        RuntimeError: If a target is not Python, the child fails, its output is
            malformed, or any import resolves outside the disposable root.
    """
    module_names: list[str] = []
    for module in modules:
        candidate = Path(module)
        if candidate.is_absolute() or candidate.suffix != ".py" or ".." in candidate.parts:
            raise RuntimeError(f"provenance target must be a Python module: {module}")
        module_names.append(".".join(candidate.with_suffix("").parts))

    with using_root_capability(workspace.root_capability or workspace.sweep_root) as capability:
        sweep_root = workspace.sweep_root.resolve(strict=True)
        real_root = workspace.real_root.resolve(strict=True)
        child_root = git_capability_path(capability)
        environment = dict(os.environ)
        retained_python_path: list[str] = []
        for entry in environment.get("PYTHONPATH", "").split(os.pathsep):
            if not entry:
                continue
            resolved_entry = Path(entry).resolve(strict=False)
            if resolved_entry == real_root or resolved_entry.is_relative_to(real_root):
                continue
            retained_python_path.append(entry)
        environment["PYTHONPATH"] = os.pathsep.join((str(child_root), *retained_python_path))
        environment["PYTHONDONTWRITEBYTECODE"] = "1"

        probe = """
import importlib
from pathlib import Path
import sys

package = importlib.import_module("vntyper")
print(Path(package.__file__).resolve())
for module_name in sys.argv[1:]:
    module = importlib.import_module(module_name)
    print(Path(module.__file__).resolve())
""".strip()
        command = [sys.executable, "-B", "-c", probe, *module_names]
        try:
            completed = subprocess.run(
                command,
                cwd=child_root,
                env=environment,
                capture_output=True,
                text=True,
                check=False,
                pass_fds=(capability.descriptor,),
            )
        except OSError as exc:
            raise RuntimeError(f"import provenance probe failed: {exc}") from exc
        if completed.returncode != 0:
            diagnostic = ((completed.stdout or "") + (completed.stderr or "")).strip() or "no child output"
            raise RuntimeError(f"import provenance probe failed: {diagnostic}")

        raw_paths = tuple(line for line in (completed.stdout or "").splitlines() if line)
        expected_count = len(module_names) + 1
        if len(raw_paths) != expected_count:
            raise RuntimeError(
                f"import provenance probe returned unexpected path count: expected {expected_count}, got {len(raw_paths)}"
            )
        try:
            paths = tuple(Path(raw_path).resolve(strict=True) for raw_path in raw_paths)
        except OSError as exc:
            raise RuntimeError(f"import provenance probe returned an invalid path: {exc}") from exc
        for path in paths:
            if not path.is_relative_to(sweep_root) or path.is_relative_to(real_root):
                raise RuntimeError(
                    f"import escaped disposable worktree: {path}; sweep root: {sweep_root}; real root: {real_root}"
                )
            print(f"Import provenance: {path}")
        return paths
