"""Public compatibility contract for per-sample report template context.

The report template directory is configurable, so values passed to Jinja are an API
for third-party templates even when VNtyper's shipped template no longer reads them.
This module records the deliberately retained deprecated paths and provides the pure
template-graph audit used to stop new context values becoming silent debt.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable, Mapping

from jinja2 import Environment, meta, nodes

logger = logging.getLogger(__name__)

#: Context paths retained for custom-template compatibility until the next major
#: release. Dotted paths name fields within a top-level context mapping.
DEPRECATED_KEYS = frozenset(
    {
        "percent_vntr_uncovered_color",
        "mean_vntr_coverage_color",
        "duplication_rate_color",
        "q20_color",
        "q30_color",
        "passed_filter_color",
        # Shipped HTML now reads Decimal-derived display strings, but configured
        # custom templates may still read these original raw fastp fractions.
        "duplication_rate",
        "q20_rate",
        "q30_rate",
        "passed_filter_rate",
        "igv_content",
        "screening_state.kestrel_result",
        "screening_state.advntr_result",
        # The cross-match sentence and its computed state. The shipped per-sample report
        # no longer renders a `Cross-Match Summary` section: the masthead already
        # carries a `Concordance` chip built from the same computed state, and for a run
        # that named one allele the allele panel prints both callers' own names beside
        # each other - so the section was a heading and a bordered box restating a fact
        # the reader had already been given twice, 1,300px further up. Both values are
        # still passed, because a configured template may render them.
        "cross_match_message",
        "cross_match_is_positive",
    }
)

#: First release allowed to remove :data:`DEPRECATED_KEYS` from the render context.
DEPRECATED_KEYS_REMOVAL_RELEASE = "3.0.0"

#: Context paths declared for custom templates and report consumers that are rendered
#: via composite containers (e.g. state_chips) in the shipped report template.
DECLARED_CONTEXT_KEYS = frozenset(
    {
        "screening_state.confidence_grade",
    }
)


def _static_access_path(node: nodes.Expr, undeclared_names: set[str]) -> str | None:
    """Return one dotted external path, collapsing non-string indexing to its root."""
    current = node
    segments: list[str] = []
    opaque = False

    while isinstance(current, (nodes.Getattr, nodes.Getitem)):
        if isinstance(current, nodes.Getattr):
            if not opaque:
                segments.append(current.attr)
        elif isinstance(current.arg, nodes.Const) and isinstance(current.arg.value, str):
            if not opaque:
                segments.append(current.arg.value)
        else:
            # Integer, slice and dynamic access address sequence or otherwise opaque
            # internals. The root value is consumed, but its runtime contents are not
            # part of the render-context contract.
            opaque = True
            segments.clear()
        current = current.node

    if not isinstance(current, nodes.Name) or current.name not in undeclared_names:
        return None
    if opaque or not segments:
        return current.name
    return ".".join((current.name, *reversed(segments)))


def _referenced_paths(parsed: nodes.Template) -> set[str]:
    """Collect maximal external access paths from one parsed Jinja template."""
    undeclared_names = set(meta.find_undeclared_variables(parsed))
    referenced: set[str] = set()

    def visit(node: nodes.Node) -> None:
        if isinstance(node, (nodes.Getattr, nodes.Getitem)):
            path = _static_access_path(node, undeclared_names)
            if path is not None:
                referenced.add(path)

            # Do not visit the access base: doing so would add every parent of a
            # maximal chain. Dynamic index expressions are independent references,
            # however, and must still participate in the audit.
            current: nodes.Expr = node
            while isinstance(current, (nodes.Getattr, nodes.Getitem)):
                if isinstance(current, nodes.Getitem) and not (
                    isinstance(current.arg, nodes.Const) and isinstance(current.arg.value, str)
                ):
                    visit(current.arg)
                current = current.node
            if not isinstance(current, nodes.Name):
                visit(current)
            return

        if isinstance(node, nodes.Name):
            if node.ctx == "load" and node.name in undeclared_names:
                referenced.add(node.name)
            return

        for child in node.iter_child_nodes():
            visit(child)

    visit(parsed)
    return referenced


def jinja_referenced_paths_recursive(environment: Environment, template_name: str) -> set[str]:
    """Return external context paths used by a template and static dependencies.

    Includes and inheritance are traversed with one visited set, so nested partials,
    base templates and cycles are all handled consistently. A dynamic include or
    inheritance target cannot be audited statically and therefore fails closed.
    Attribute access and constant string-key access become dotted mapping paths;
    sequence, slice and dynamic indexing consume only their top-level root.

    Args:
        environment: Jinja environment whose loader owns the template graph.
        template_name: Loader-relative entry template name.

    Returns:
        set[str]: Maximal context paths referenced anywhere in the graph.

    Raises:
        ValueError: If the environment has no loader or a template dependency is
            selected dynamically.
        TemplateNotFound: If the loader cannot resolve a static template reference.
    """
    loader = environment.loader
    if loader is None:
        msg = "A Jinja environment loader is required to audit referenced template paths."
        logger.error(msg)
        raise ValueError(msg)

    referenced_paths: set[str] = set()
    pending = [template_name]
    visited: set[str] = set()

    while pending:
        current = pending.pop()
        if current in visited:
            continue
        visited.add(current)

        source, _, _ = loader.get_source(environment, current)
        parsed = environment.parse(source)
        referenced_paths.update(_referenced_paths(parsed))

        for dependency in meta.find_referenced_templates(parsed):
            if dependency is None:
                msg = f"Template {current!r} contains a dynamic Jinja template reference that cannot be audited."
                logger.error(msg)
                raise ValueError(msg)
            pending.append(dependency)

    return referenced_paths


def _flatten_runtime_value(path: str, value: object) -> set[str]:
    """Flatten string-keyed mappings; treat every other runtime value as opaque."""
    if isinstance(value, Mapping) and all(isinstance(key, str) for key in value):
        if not value:
            return {path}
        flattened: set[str] = set()
        for key, child in value.items():
            flattened.update(_flatten_runtime_value(f"{path}.{key}", child))
        return flattened
    return {path}


def flatten_runtime_context_paths(context: Mapping[str, object]) -> set[str]:
    """Return auditable leaf paths from a report render context.

    String-keyed mappings are recursively flattened because Jinja can reference their
    members statically. Empty mappings remain leaf values. Sequences and mappings with
    non-string keys stay opaque at their containing path: requiring their element
    internals would confuse iteration data with the renderer's named context contract.

    Args:
        context: Top-level mapping passed to ``Template.render``.

    Returns:
        set[str]: Dotted mapping leaves and opaque container paths.
    """
    flattened: set[str] = set()
    for key, value in context.items():
        flattened.update(_flatten_runtime_value(key, value))
    return flattened


def unreferenced_runtime_context_paths(
    context: Mapping[str, object],
    *,
    referenced_paths: Iterable[str],
    deprecated_paths: Iterable[str],
    declared_paths: Iterable[str] = DECLARED_CONTEXT_KEYS,
) -> set[str]:
    """Return nondeprecated runtime paths not consumed by a template graph.

    A direct reference to a mapping parent consumes all of that mapping's leaves.
    Deprecation exemptions remain exact: deprecating one dotted member never exempts
    its siblings or descendants. Declared paths identify public custom-template
    API extensions rendered via composite structures in the shipped template.

    Args:
        context: Top-level mapping passed to ``Template.render``.
        referenced_paths: Static paths collected from the template graph.
        deprecated_paths: Exact compatibility paths exempt from the parity gate.
        declared_paths: Declared public API paths exempt from the dead-context gate.

    Returns:
        set[str]: Runtime paths that represent silent render-context debt.
    """
    referenced = set(referenced_paths)
    deprecated = set(deprecated_paths)
    declared = set(declared_paths)
    return {
        path
        for path in flatten_runtime_context_paths(context)
        if path not in deprecated
        and path not in declared
        and not any(path == reference or path.startswith(f"{reference}.") for reference in referenced)
    }
