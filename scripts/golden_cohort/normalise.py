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

#: Keys dropped from every ``pipeline_summary.json`` step record, unconditionally, with
#: why. Both are pure wall-clock.
DROPPED_KEYS: dict[str, str] = {
    "start": "per-run wall-clock",
    "end": "per-run wall-clock",
}

#: Basenames of the step result files :func:`golden_cohort.artifacts.read_pipeline_case`
#: parses and compares directly, field by field. Only for these is a step's ``md5sum``
#: genuinely redundant.
DIRECTLY_COMPARED_RESULT_FILES: frozenset[str] = frozenset(
    {
        "kestrel_result.tsv",
        "coverage_summary.tsv",
        "output_adVNTR_result.tsv",
    }
)

#: Why ``md5sum`` is dropped where it is dropped.
#:
#: It used to be dropped from *every* step record, justified by "the file's content is
#: compared directly instead". That justification was true of three of the six step result
#: files this pipeline writes and false of the other three, so a third of the checksums
#: were being discarded with no comparator behind them - and two of those three carry
#: decisions:
#:
#: * ``fastq_bam_processing/pipeline_info.json`` holds the assembly guard's own verdict
#:   (``assembly_text``, ``assembly_contig``, ``alignment_pipeline``, ``warning``). Nothing
#:   in :func:`golden_cohort.artifacts.read_pipeline_case` opens that file, so with the
#:   checksum gone the guard's verdict was compared by nothing at all.
#: * ``advntr/cross_match_results.tsv`` holds the per-variant cross-match table. The
#:   harness reads only the rendered cross-match *sentence* out of the HTML report, so the
#:   table behind it was likewise uncompared.
#: * ``fastq_bam_processing/output_R1.fastq.gz`` is the BAM-to-FASTQ product, read by
#:   nothing downstream of the run itself.
#:
#: The checksum is therefore kept for anything without a named direct comparator, and this
#: is measured rather than hoped: across run 4's 59 cases the ``md5sum`` of
#: ``pipeline_info.json`` (59/59), ``output_R1.fastq.gz`` (59/59) and
#: ``cross_match_results.tsv`` (3/3) is identical between the two sides, so restoring the
#: comparison adds no delta to the recorded attestation. ``kestrel_result.tsv`` is the one
#: file whose checksum differs on 59 of 59 - its ``##`` banner carries the analysis
#: timestamp and the VNtyper version - and it is also directly compared, banner normalised,
#: so it is exactly the file the original justification was written for.
MD5SUM_WHY_DROPPED = (
    "checksum of a file the harness parses and compares field by field, whose ## provenance banner carries "
    "the run timestamp and so differs on every run of the same code"
)

#: Why it is kept everywhere else.
MD5SUM_WHY_KEPT = (
    "checksum of a step result file with no direct comparator in read_pipeline_case; it is the only thing "
    "comparing that file's content at all"
)


def strip_step_record(step: dict[str, Any]) -> dict[str, Any]:
    """Remove the keys of one step record that must not take part in the comparison.

    Args:
        step: One entry of ``pipeline_summary.json``'s ``steps``.

    Returns:
        dict[str, Any]: A new mapping without the wall-clock keys, without
        ``parsed_result``, and without ``md5sum`` only when the step's ``result_file`` is
        one the harness compares directly.

    Note:
        ``parsed_result`` is dropped because it is the pipeline's own re-parse of the same
        result file, so for the three directly-compared files it duplicates the comparison.
        For the other three it is *not* duplicated - which is why their ``md5sum`` is now
        kept, and the checksum rather than the parse is what stands in for them. Comparing
        ``parsed_result`` itself would be the wider fix; it is not taken here because its
        effect on the recorded runs has not been measured, and an unmeasured widening of a
        gate is how this harness acquired the claims it is being corrected for.
    """
    result_file = step.get("result_file")
    basename = Path(str(result_file)).name if result_file else ""
    drop = set(DROPPED_KEYS) | {"parsed_result"}
    if basename in DIRECTLY_COMPARED_RESULT_FILES:
        drop.add("md5sum")
    return {key: value for key, value in step.items() if key not in drop}


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
            "name": "md5sum",
            "kind": "dropped-key-conditional",
            "pattern": f"md5sum, only on steps writing {', '.join(sorted(DIRECTLY_COMPARED_RESULT_FILES))}",
            "why": MD5SUM_WHY_DROPPED,
        }
    )
    entries.append(
        {
            "name": "md5sum",
            "kind": "compared",
            "pattern": "md5sum, on every other step",
            "why": MD5SUM_WHY_KEPT,
        }
    )
    entries.append(
        {
            "name": "parsed_result",
            "kind": "dropped-key",
            "pattern": "parsed_result",
            "why": (
                "the pipeline's own re-parse of the step's result file. For the directly-compared files it "
                "duplicates the row-set comparison; for the rest the retained md5sum stands in for it. See "
                "normalise.strip_step_record."
            ),
        }
    )
    entries.append(
        {
            "name": "cohort_row_order",
            "kind": "sort",
            "pattern": "cohort tables and exports",
            "why": (
                "cohort sample order is not comparable across a version boundary: a baseline predating the "
                "determinism fix iterates the discovery set directly, and a ZIP input extracts to a "
                "randomly-named tempdir on any version. Rows are sorted before comparison and each side's raw "
                "order is recorded, uncompared, as cohort_sample_order_raw. See "
                "golden_cohort.compare.COHORT_ORDER_WHY for the full statement and the tests that do attest "
                "the ordering."
            ),
        }
    )
    return entries
