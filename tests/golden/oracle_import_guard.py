"""Standard-library guard for independent golden-oracle import closures."""

from __future__ import annotations

import ast
import sys
from pathlib import Path


def assert_independent_import_closure(entrypoint: Path, repository_root: Path) -> tuple[Path, ...]:
    """Recursively AST-scan an oracle's stdlib-and-local-only import closure.

    Args:
        entrypoint: Python file at the root of the closure.
        repository_root: Directory within which imports count as local.

    Returns:
        Every recursively scanned local source file.

    Raises:
        AssertionError: If a file is absent or imports production or third-party code.
    """
    root = repository_root.resolve()
    resolved_entrypoint = entrypoint.resolve()
    oracle_root = resolved_entrypoint.parent
    pending = [resolved_entrypoint]
    scanned: set[Path] = set()
    forbidden: list[str] = []

    while pending:
        source = pending.pop()
        if source in scanned:
            continue
        if not source.is_file():
            raise AssertionError(f"oracle import-closure source does not exist: {source}")
        try:
            source.relative_to(root)
        except ValueError as error:
            raise AssertionError(f"oracle import-closure source escapes repository root: {source}") from error
        scanned.add(source)
        tree = ast.parse(source.read_text(encoding="utf-8"), filename=str(source))
        importlib_aliases = _import_bindings(tree, "importlib")
        builtins_aliases = _import_bindings(tree, "builtins")
        import_module_aliases = {
            alias.asname or alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module == "importlib"
            for alias in node.names
            if alias.name == "import_module"
        }

        for node in ast.walk(tree):
            modules: list[tuple[str, int, bool]] = []
            if isinstance(node, ast.Import):
                modules.extend((alias.name, 0, False) for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                modules.append((node.module or "", node.level, False))
                modules.extend(
                    (f"{node.module}.{alias.name}" if node.module else alias.name, node.level, True)
                    for alias in node.names
                    if alias.name != "*"
                )
            elif isinstance(node, ast.Call):
                dynamic_module = _literal_dynamic_import(
                    node,
                    importlib_aliases,
                    import_module_aliases,
                    builtins_aliases,
                )
                if dynamic_module is not None:
                    modules.append((dynamic_module, 0, False))
            else:
                continue

            for module, level, imported_symbol_allowed in modules:
                if module == "vntyper" or module.startswith("vntyper."):
                    forbidden.append(f"{source}: imports {module}")
                    continue
                local_paths = _resolve_local_imports(source, root, module, level)
                if local_paths:
                    if _is_oracle_local(local_paths, oracle_root):
                        pending.extend(path for path in local_paths if path not in scanned)
                    else:
                        forbidden.append(f"{source}: imports non-oracle local module {module}")
                    continue
                if imported_symbol_allowed and "." in module:
                    parent_module = module.rsplit(".", 1)[0]
                    parent_paths = _resolve_local_imports(source, root, parent_module, level)
                    if parent_paths and _is_oracle_local(parent_paths, oracle_root):
                        pending.extend(path for path in parent_paths if path not in scanned)
                        continue
                if _is_oracle_namespace(source, root, module, level, oracle_root):
                    continue
                if level == 0 and module.split(".", 1)[0] in sys.stdlib_module_names:
                    continue
                forbidden.append(f"{source}: imports non-stdlib module {module or '<unresolved relative import>'}")

    if forbidden:
        raise AssertionError("forbidden import in oracle closure: " + "; ".join(sorted(forbidden)))
    return tuple(sorted(scanned))


def _literal_dynamic_import(
    node: ast.Call,
    importlib_aliases: set[str],
    import_module_aliases: set[str],
    builtins_aliases: set[str],
) -> str | None:
    """Return the literal target of a supported dynamic-import spelling."""
    argument = node.args[0] if node.args else next((item.value for item in node.keywords if item.arg == "name"), None)
    if not isinstance(argument, ast.Constant) or not isinstance(argument.value, str):
        return None
    if isinstance(node.func, ast.Name) and node.func.id == "__import__":
        return argument.value
    if (
        isinstance(node.func, ast.Attribute)
        and node.func.attr == "__import__"
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id in builtins_aliases
    ):
        return argument.value
    if isinstance(node.func, ast.Name) and node.func.id in import_module_aliases:
        return argument.value
    if (
        isinstance(node.func, ast.Attribute)
        and node.func.attr == "import_module"
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id in importlib_aliases
    ):
        return argument.value
    return None


def _import_bindings(tree: ast.AST, module: str) -> set[str]:
    bindings: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Import):
            continue
        for alias in node.names:
            if alias.name == module:
                bindings.add(alias.asname or module)
            elif alias.name.startswith(f"{module}.") and alias.asname is None:
                bindings.add(module)
    return bindings


def _resolve_local_imports(source: Path, root: Path, module: str, level: int) -> tuple[Path, ...]:
    candidates = _local_candidates(source, root, module, level)

    for base, candidate_parts in candidates:
        candidate = base.joinpath(*candidate_parts)
        file_candidate = candidate.with_suffix(".py")
        package_candidate = candidate / "__init__.py"
        if file_candidate.is_file():
            return _package_initializers(base, candidate_parts[:-1]) + (file_candidate.resolve(),)
        if package_candidate.is_file():
            return _package_initializers(base, candidate_parts)
    return ()


def _local_candidates(source: Path, root: Path, module: str, level: int) -> list[tuple[Path, list[str]]]:
    parts = module.split(".") if module else []
    if not level:
        return [(root, parts), (source.parent, parts)]
    base = source.parent
    for _ in range(level - 1):
        base = base.parent
    return [(base, parts)]


def _package_initializers(base: Path, parts: list[str]) -> tuple[Path, ...]:
    initializers: list[Path] = []
    for index in range(1, len(parts) + 1):
        initializer = base.joinpath(*parts[:index], "__init__.py")
        if initializer.is_file():
            initializers.append(initializer.resolve())
    return tuple(initializers)


def _is_oracle_local(paths: tuple[Path, ...], oracle_root: Path) -> bool:
    return any(path.is_relative_to(oracle_root) for path in paths)


def _is_oracle_namespace(source: Path, root: Path, module: str, level: int, oracle_root: Path) -> bool:
    return any(
        base.joinpath(*parts).is_dir() and base.joinpath(*parts).resolve().is_relative_to(oracle_root)
        for base, parts in _local_candidates(source, root, module, level)
    )
