"""The public report-template context and its recursive shipped-template audit."""

from __future__ import annotations

import pytest
from jinja2 import DictLoader, Environment

from vntyper.scripts import report_context_contract as contract

pytestmark = pytest.mark.unit


def test_recursive_name_resolution_follows_nested_includes_and_inheritance() -> None:
    """Stopping at the entry template misses names supplied to shared template layers."""
    environment = Environment(
        loader=DictLoader(
            {
                "report.html": '{% extends "layout.html" %}{% include "section.html" %}{{ report_value }}',
                "layout.html": "{{ layout_value }}{% block body %}{% endblock %}",
                "section.html": '{% include "leaf.html" %}{{ section_value }}',
                "leaf.html": "{{ leaf_value }}",
            }
        )
    )

    names = contract.jinja_referenced_names_recursive(environment, "report.html")

    assert names == {"report_value", "layout_value", "section_value", "leaf_value"}


def test_recursive_name_resolution_terminates_on_an_include_cycle() -> None:
    """A custom template cycle is visited once rather than recursed into forever."""
    environment = Environment(
        loader=DictLoader(
            {
                "first.html": '{% include "second.html" %}{{ first_value }}',
                "second.html": '{% include "first.html" %}{{ second_value }}',
            }
        )
    )

    names = contract.jinja_referenced_names_recursive(environment, "first.html")

    assert names == {"first_value", "second_value"}


def test_recursive_name_resolution_refuses_a_dynamic_template_reference() -> None:
    """An unresolved include cannot silently make the context audit incomplete."""
    environment = Environment(loader=DictLoader({"report.html": "{% include chosen_template %}"}))

    with pytest.raises(ValueError, match="dynamic Jinja template reference"):
        contract.jinja_referenced_names_recursive(environment, "report.html")


def test_recursive_name_resolution_requires_an_environment_loader() -> None:
    """Without a loader there is no template graph to audit, so the resolver fails closed."""
    with pytest.raises(ValueError, match="loader"):
        contract.jinja_referenced_names_recursive(Environment(), "report.html")


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
                "igv_content",
                "screening_state.kestrel_result",
                "screening_state.advntr_result",
            }
        )
        == contract.DEPRECATED_KEYS
    )
    assert contract.DEPRECATED_KEYS_REMOVAL_RELEASE == "3.0.0"
