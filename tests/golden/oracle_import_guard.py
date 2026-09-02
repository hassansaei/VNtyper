"""Standard-library guard for independent golden-oracle import closures."""

from __future__ import annotations

import ast
from pathlib import Path


def assert_independent_import_closure(entrypoint: Path, repository_root: Path) -> tuple[Path, ...]:
    """Recursively AST-scan local imports and reject every production import.

    Args:
        entrypoint: Python file at the root of the closure.
        repository_root: Directory within which imports count as local.

    Returns:
        Every recursively scanned local source file.

    Raises:
        AssertionError: If a file is absent or any closure member imports ``vntyper``.
    """
    root = repository_root.resolve()
    pending = [entrypoint.resolve()]
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
        importlib_aliases = {
            alias.asname or alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
            if alias.name == "importlib"
        }
        import_module_aliases = {
            alias.asname or alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module == "importlib"
            for alias in node.names
            if alias.name == "import_module"
        }

        for node in ast.walk(tree):
            modules: list[tuple[str, int]] = []
            if isinstance(node, ast.Import):
                modules.extend((alias.name, 0) for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                modules.append((node.module or "", node.level))
                modules.extend(
                    (f"{node.module}.{alias.name}" if node.module else alias.name, node.level)
                    for alias in node.names
                    if alias.name != "*"
                )
            elif isinstance(node, ast.Call):
                dynamic_module = _literal_dynamic_import(node, importlib_aliases, import_module_aliases)
                if dynamic_module is not None:
                    modules.append((dynamic_module, 0))
            else:
                continue

            for module, level in modules:
                if module == "vntyper" or module.startswith("vntyper."):
                    forbidden.append(f"{source}: imports {module}")
                    continue
                pending.extend(
                    path for path in _resolve_local_imports(source, root, module, level) if path not in scanned
                )

    if forbidden:
        raise AssertionError("forbidden production import in oracle closure: " + "; ".join(sorted(forbidden)))
    return tuple(sorted(scanned))


def _literal_dynamic_import(
    node: ast.Call,
    importlib_aliases: set[str],
    import_module_aliases: set[str],
) -> str | None:
    """Return the literal target of a supported dynamic-import spelling."""
    argument = node.args[0] if node.args else next((item.value for item in node.keywords if item.arg == "name"), None)
    if not isinstance(argument, ast.Constant) or not isinstance(argument.value, str):
        return None
    if isinstance(node.func, ast.Name) and node.func.id == "__import__":
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


def _resolve_local_imports(source: Path, root: Path, module: str, level: int) -> tuple[Path, ...]:
    if level:
        base = source.parent
        for _ in range(level - 1):
            base = base.parent
        parts = module.split(".") if module else []
        candidates = [(base, parts)]
    else:
        parts = module.split(".") if module else []
        candidates = [(root, parts), (source.parent, parts)]

    for base, candidate_parts in candidates:
        candidate = base.joinpath(*candidate_parts)
        file_candidate = candidate.with_suffix(".py")
        package_candidate = candidate / "__init__.py"
        if file_candidate.is_file():
            return _package_initializers(base, candidate_parts[:-1]) + (file_candidate.resolve(),)
        if package_candidate.is_file():
            return _package_initializers(base, candidate_parts)
    return ()


def _package_initializers(base: Path, parts: list[str]) -> tuple[Path, ...]:
    initializers: list[Path] = []
    for index in range(1, len(parts) + 1):
        initializer = base.joinpath(*parts[:index], "__init__.py")
        if initializer.is_file():
            initializers.append(initializer.resolve())
    return tuple(initializers)
