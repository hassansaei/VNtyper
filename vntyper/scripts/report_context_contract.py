"""Public compatibility contract for per-sample report template context.

The report template directory is configurable, so values passed to Jinja are an API
for third-party templates even when VNtyper's shipped template no longer reads them.
This module records the deliberately retained deprecated paths and provides the pure
template-graph audit used to stop new context values becoming silent debt.
"""

from __future__ import annotations

import logging

from jinja2 import Environment, meta

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
        "igv_content",
        "screening_state.kestrel_result",
        "screening_state.advntr_result",
    }
)

#: First release allowed to remove :data:`DEPRECATED_KEYS` from the render context.
DEPRECATED_KEYS_REMOVAL_RELEASE = "3.0.0"


def jinja_referenced_names_recursive(environment: Environment, template_name: str) -> set[str]:
    """Return undeclared names used by a template and every static dependency.

    Includes and inheritance are traversed with one visited set, so nested partials,
    base templates and cycles are all handled consistently. A dynamic include or
    inheritance target cannot be audited statically and therefore fails closed.

    Args:
        environment: Jinja environment whose loader owns the template graph.
        template_name: Loader-relative entry template name.

    Returns:
        set[str]: Top-level context names referenced anywhere in the graph.

    Raises:
        ValueError: If the environment has no loader or a template dependency is
            selected dynamically.
        TemplateNotFound: If the loader cannot resolve a static template reference.
    """
    loader = environment.loader
    if loader is None:
        msg = "A Jinja environment loader is required to audit referenced template names."
        logger.error(msg)
        raise ValueError(msg)

    referenced_names: set[str] = set()
    pending = [template_name]
    visited: set[str] = set()

    while pending:
        current = pending.pop()
        if current in visited:
            continue
        visited.add(current)

        source, _, _ = loader.get_source(environment, current)
        parsed = environment.parse(source)
        referenced_names.update(meta.find_undeclared_variables(parsed))

        for dependency in meta.find_referenced_templates(parsed):
            if dependency is None:
                msg = f"Template {current!r} contains a dynamic Jinja template reference that cannot be audited."
                logger.error(msg)
                raise ValueError(msg)
            pending.append(dependency)

    return referenced_names
