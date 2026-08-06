"""The substitutions applied to everything that differs between the two sides by construction.

A gate that reports "58 of 58 cases differ" because each side wrote its own output
directory into its own logs has measured nothing. This module holds the rules that remove
exactly those differences, and only those - and, because a normalisation is a claim that a
difference does not matter, :func:`manifest` renders the rules that were applied so the
claim ships with the result rather than living in someone's memory.

Every rule here is either a path that is per-side by definition, or a clock. Nothing here
touches a genotype field, a motif, a confidence, a flag or a coverage number.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

#: Keys dropped from ``pipeline_summary.json`` step records before comparison, with why.
#: ``md5sum`` is the checksum of a file whose ``##`` provenance banner carries the analysis
#: timestamp, so it differs on every run of the same code; the file's *content* is compared
#: directly instead, banner removed, which is what the gate page describes.
DROPPED_KEYS: dict[str, str] = {
    "md5sum": "checksum of a file whose provenance banner carries the run timestamp; content is compared instead",
    "start": "per-run wall-clock",
    "end": "per-run wall-clock",
}


@dataclass(frozen=True)
class Rule:
    """One normalisation.

    Attributes:
        name: Short identifier, used in the manifest and in reports.
        description: Why this difference is by construction.
        pattern: What to replace.
        replacement: What to replace it with.
    """

    name: str
    description: str
    pattern: re.Pattern[str]
    replacement: str


def _path_rules(name: str, description: str, path: Path, placeholder: str) -> list[Rule]:
    """Build the rules that erase one directory, following symlinks as well.

    ``tests/data`` and ``reference/`` are symlinked into the gate's trees, so a run can
    record either the symlinked path or its target. Both are erased to the same placeholder.

    Args:
        name: Rule name stem.
        description: Why this path is per-side.
        path: The directory to erase.
        placeholder: What to put in its place.

    Returns:
        list[Rule]: One rule per distinct spelling of the path, longest first.
    """
    spellings = {str(path), str(path.resolve())}
    return [
        Rule(
            name=f"{name}[{index}]" if index else name,
            description=description,
            pattern=re.compile(re.escape(spelling)),
            replacement=placeholder,
        )
        for index, spelling in enumerate(sorted(spellings, key=len, reverse=True))
    ]


def build_rules(
    *,
    source_root: Path,
    run_root: Path,
    data_root: Path | None = None,
    extra: dict[str, Path] | None = None,
) -> list[Rule]:
    """Build the rule list for one side.

    Order matters twice over. Paths go before clocks, because a path can contain digits a
    timestamp rule would otherwise chew into. And ``data_root`` goes before ``source_root``,
    because both sides read the same ``tests/data`` by absolute path while only one side's
    source tree contains it - normalising the source tree first would erase the shared data
    path on that side and leave it literal on the other, manufacturing a delta in every
    recorded command. The harness's own smoke test found exactly that.

    Args:
        source_root: The side's source tree.
        run_root: The side's gate output directory.
        data_root: The shared test-data directory, which is **not** per-side.
        extra: Further per-side paths to erase, name -> path.

    Returns:
        list[Rule]: The rules, in application order.
    """
    rules: list[Rule] = []
    if data_root is not None:
        rules.extend(
            _path_rules("data_root", "the shared test-data directory, read by both sides", data_root, "<DATA_ROOT>")
        )
    rules.extend(_path_rules("run_root", "the per-side gate output directory", run_root, "<RUN_ROOT>"))
    rules.extend(_path_rules("source_root", "the per-side source tree", source_root, "<SOURCE_ROOT>"))
    for name, path in sorted((extra or {}).items()):
        rules.extend(_path_rules(name, f"per-side path {name}", path, f"<{name.upper()}>"))

    rules.extend(
        [
            Rule(
                name="tempdir",
                description="per-run temporary directories",
                pattern=re.compile(r"/tmp/[A-Za-z0-9_./-]+"),
                replacement="<TMP>",
            ),
            Rule(
                name="iso_timestamp",
                description="wall-clock timestamps, ISO-8601 with or without a T separator",
                pattern=re.compile(r"\b\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}(?:\.\d+)?\b"),
                replacement="<TIMESTAMP>",
            ),
            Rule(
                name="date",
                description="wall-clock dates left over once timestamps are gone",
                pattern=re.compile(r"\b\d{4}-\d{2}-\d{2}\b"),
                replacement="<DATE>",
            ),
            Rule(
                name="runtime",
                description="the cohort's computed per-sample runtime, which is a clock difference",
                pattern=re.compile(r"\b\d+\.\d{2} seconds\b"),
                replacement="<RUNTIME>",
            ),
            Rule(
                name="vntyper_version",
                description="the version string, which necessarily differs between two commits",
                pattern=re.compile(r"(?<=VNtyper Version: )\S+"),
                replacement="<VERSION>",
            ),
            Rule(
                name="png_data_uri",
                description=(
                    "base64 PNG payloads: matplotlib output is not byte-stable across runs and the "
                    "chart's numbers are compared separately, from the plotly figure JSON"
                ),
                pattern=re.compile(r"(?<=data:image/png;base64,)[A-Za-z0-9+/=]+"),
                replacement="<PNG>",
            ),
            Rule(
                name="uuid",
                description="plotly's per-figure div identifiers",
                pattern=re.compile(r"\b[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}\b"),
                replacement="<UUID>",
            ),
        ]
    )
    return rules


def apply(text: str, rules: list[Rule]) -> str:
    """Apply every rule to one string, in order.

    Args:
        text: The string to normalise.
        rules: The rules from :func:`build_rules`.

    Returns:
        str: The normalised string.
    """
    for rule in rules:
        text = rule.pattern.sub(rule.replacement, text)
    return text


def apply_deep(value: Any, rules: list[Rule]) -> Any:
    """Apply the rules to every string inside a JSON-shaped value.

    Args:
        value: A string, list, dict or scalar.
        rules: The rules from :func:`build_rules`.

    Returns:
        Any: The same shape, with every string normalised.
    """
    if isinstance(value, str):
        return apply(value, rules)
    if isinstance(value, list):
        return [apply_deep(item, rules) for item in value]
    if isinstance(value, dict):
        return {key: apply_deep(item, rules) for key, item in value.items()}
    return value


def manifest(rules: list[Rule]) -> list[dict[str, str]]:
    """Render the rules that were applied, for the result document.

    Args:
        rules: The rules from :func:`build_rules`.

    Returns:
        list[dict[str, str]]: One entry per rule, plus the structurally dropped keys.
    """
    entries = [
        {"name": rule.name, "kind": "substitution", "pattern": rule.pattern.pattern, "why": rule.description}
        for rule in rules
    ]
    entries.extend(
        {"name": key, "kind": "dropped-key", "pattern": key, "why": why} for key, why in sorted(DROPPED_KEYS.items())
    )
    entries.append(
        {
            "name": "cohort_row_order",
            "kind": "sort",
            "pattern": "cohort tables and exports",
            "why": (
                "vntyper cohort discovers samples into a set, so their order is not reproducible even "
                "between two runs of the same code (see tests/unit/test_cohort_inputs.py::"
                "test_discovery_returns_an_unordered_set_today). Rows are sorted before comparison and "
                "each side's raw order is recorded, uncompared, as cohort_sample_order_raw."
            ),
        }
    )
    return entries
