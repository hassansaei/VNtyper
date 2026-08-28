"""The public report-template context and its recursive shipped-template audit."""

from __future__ import annotations

import pytest
from jinja2 import DictLoader, Environment, TemplateNotFound

from vntyper.scripts import report_context_contract as contract

pytestmark = pytest.mark.unit


def test_recursive_path_resolution_follows_nested_includes_and_inheritance() -> None:
    """Stopping at the entry template misses paths supplied to shared template layers."""
    environment = Environment(
        loader=DictLoader(
            {
                "report.html": ('{% extends "layout.html" %}{% include "section.html" %}{{ report_value.title }}'),
                "layout.html": '{{ layout_value["heading"] }}{% block body %}{% endblock %}',
                "section.html": '{% include "leaf.html" %}{{ section_value.label }}',
                "leaf.html": "{{ leaf_value }}",
            }
        )
    )

    paths = contract.jinja_referenced_paths_recursive(environment, "report.html")

    assert paths == {
        "report_value.title",
        "layout_value.heading",
        "section_value.label",
        "leaf_value",
    }


def test_recursive_path_resolution_terminates_on_an_include_cycle() -> None:
    """A custom template cycle is visited once rather than recursed into forever."""
    environment = Environment(
        loader=DictLoader(
            {
                "first.html": '{% include "second.html" %}{{ first_value }}',
                "second.html": '{% include "first.html" %}{{ second_value }}',
            }
        )
    )

    paths = contract.jinja_referenced_paths_recursive(environment, "first.html")

    assert paths == {"first_value", "second_value"}


def test_recursive_path_resolution_refuses_a_dynamic_template_reference() -> None:
    """An unresolved include cannot silently make the context audit incomplete."""
    environment = Environment(loader=DictLoader({"report.html": "{% include chosen_template %}"}))

    with pytest.raises(ValueError, match="dynamic Jinja template reference"):
        contract.jinja_referenced_paths_recursive(environment, "report.html")


def test_recursive_path_resolution_requires_an_environment_loader() -> None:
    """Without a loader there is no template graph to audit, so the resolver fails closed."""
    with pytest.raises(ValueError, match="loader"):
        contract.jinja_referenced_paths_recursive(Environment(), "report.html")


@pytest.mark.parametrize("directive", ['{% include "gone.html" %}', '{% extends "gone.html" %}'])
def test_recursive_path_resolution_propagates_a_missing_static_dependency(directive) -> None:
    """A missing static include or base must fail the audit rather than be skipped."""
    environment = Environment(loader=DictLoader({"report.html": directive}))

    with pytest.raises(TemplateNotFound, match="gone.html"):
        contract.jinja_referenced_paths_recursive(environment, "report.html")


def test_referenced_paths_track_mapping_members_without_parents_or_sequence_internals() -> None:
    """Nested mappings need leaf precision while sequence contents remain opaque."""
    environment = Environment(
        loader=DictLoader(
            {
                "report.html": (
                    '{{ screening_state.matched_rule }}{{ screening_state["emphasis"] }}'
                    "{% for row in rows %}{{ row.value }}{% endfor %}"
                    '{{ rows[0].label }}{{ rows.values[0].label }}{{ rows["values"][0].label }}'
                    "{{ factory().value }}"
                )
            }
        )
    )

    paths = contract.jinja_referenced_paths_recursive(environment, "report.html")

    assert paths == {"screening_state.matched_rule", "screening_state.emphasis", "rows", "factory"}


def test_runtime_context_paths_flatten_string_mappings_but_not_container_internals() -> None:
    """Only statically addressable mapping leaves participate in nested parity."""
    context = {
        "screening_state": {"matched_rule": "rule", "empty": {}},
        "rows": [{"value": 1}],
        "integer_keyed": {0: "opaque"},
        "plain": "value",
    }

    paths = contract.flatten_runtime_context_paths(context)

    assert paths == {
        "screening_state.matched_rule",
        "screening_state.empty",
        "rows",
        "integer_keyed",
        "plain",
    }


def test_unreferenced_paths_exempt_exact_deprecations_and_accept_a_parent_reference() -> None:
    """A whole-mapping consumer covers leaves, but a dotted exemption covers only itself."""
    context = {
        "screening_state": {
            "matched_rule": "rule",
            "legacy": {"child": "old"},
            "unused": "debt",
        },
        "metadata": {"sample": "sample-1", "run": "run-1"},
    }

    unused = contract.unreferenced_runtime_context_paths(
        context,
        referenced_paths={"screening_state.matched_rule", "metadata"},
        deprecated_paths={"screening_state.legacy"},
    )

    assert unused == {"screening_state.legacy.child", "screening_state.unused"}


def test_every_deprecated_context_path_has_an_explicit_major_release() -> None:
    """Legacy API keys stay available through 2.x and have one unambiguous removal release."""
    assert (
        frozenset(
            {
                "percent_vntr_uncovered_color",
                "mean_vntr_coverage_color",
                "duplication_rate_color",
                "q20_color",
                "q30_color",
                "passed_filter_color",
                "duplication_rate",
                "q20_rate",
                "q30_rate",
                "passed_filter_rate",
                "igv_content",
                "screening_state.kestrel_result",
                "screening_state.advntr_result",
                "cross_match_message",
                "cross_match_is_positive",
            }
        )
        == contract.DEPRECATED_KEYS
    )
    assert contract.DEPRECATED_KEYS_REMOVAL_RELEASE == "3.0.0"
