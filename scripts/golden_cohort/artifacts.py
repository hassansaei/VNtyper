"""Read one run's comparable artefacts off disk, normalised.

Everything the gate page names as compared is read here, plus the cohort artefacts the
page's "What this gate does not cover" section says it never saw. An artefact that does not
exist reads as ``None`` and is reported as missing on that side, which is a delta - a file
that stops being written is exactly the kind of change a gate is for.

The parsers are deliberately stdlib-only. ``pandas.read_html`` needs lxml or bs4 and would
make the instrument depend on which of those happened to be installed on the day.
"""

from __future__ import annotations

import json
import logging
import re
from html.parser import HTMLParser
from pathlib import Path
from typing import Any

from golden_cohort import normalise
from golden_cohort.normalise import Rule

logger = logging.getLogger(__name__)

#: The columns ``kestrel_result.tsv`` rows are keyed on, as the gate page specifies.
KESTREL_KEY = ("Motifs", "POS", "REF", "ALT", "Variant")

#: adVNTR rows are keyed on the VNTR id and the called state.
ADVNTR_KEY = ("VID", "State")

#: ``<p class="summary-box ...">`` - the screening summary box and, when adVNTR ran, the
#: cross-match box. The class list is what carries the computed emphasis.
SUMMARY_BOX_RE = re.compile(r'<p class="summary-box([^"]*)">\s*(.*?)\s*</p>', re.DOTALL)

#: Plotly writes the donut's category counts into the figure JSON. Extracting them is how
#: the cohort's category counts are compared without depending on a base64 PNG.
DONUT_VALUES_RE = re.compile(r'"values":\s*\[([0-9,.\s]*)\]')

#: Every ``"text"`` value in the embedded figure JSON. Plotly escapes markup so it cannot
#: break out of the ``<script>`` block (``<b>2</b>`` is written ``<b>2</b>``),
#: so the value is decoded as a JSON string before the bold total is read out of it.
DONUT_TEXT_RE = re.compile(r'"text":\s*("(?:[^"\\]|\\.)*")')

#: The donut's centre annotation, i.e. the cohort total for that algorithm.
DONUT_TOTAL_RE = re.compile(r"<b>([^<]*)</b>")


class _TableExtractor(HTMLParser):
    """Pull every ``<table>`` out of a rendered report as header plus rows."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.tables: list[dict[str, Any]] = []
        self._table: dict[str, Any] | None = None
        self._row: list[str] | None = None
        self._cell: list[str] | None = None
        self._is_header = False

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "table":
            self._table = {"header": [], "rows": []}
        elif tag == "tr" and self._table is not None:
            self._row = []
            self._is_header = False
        elif tag in ("td", "th") and self._row is not None:
            self._cell = []
            self._is_header = self._is_header or tag == "th"

    def handle_endtag(self, tag: str) -> None:
        if tag in ("td", "th") and self._cell is not None and self._row is not None:
            self._row.append("".join(self._cell).strip())
            self._cell = None
        elif tag == "tr" and self._row is not None and self._table is not None:
            if self._is_header and not self._table["header"]:
                self._table["header"] = self._row
            else:
                self._table["rows"].append(self._row)
            self._row = None
        elif tag == "table" and self._table is not None:
            self.tables.append(self._table)
            self._table = None

    def handle_data(self, data: str) -> None:
        if self._cell is not None:
            self._cell.append(data)


def read_tsv(path: Path, rules: list[Rule]) -> dict[str, Any] | None:
    """Read a pipeline TSV into columns, rows and its ``##`` provenance banner.

    The banner is split out rather than dropped: it carries the analysis timestamp and the
    VNtyper version, both of which differ by construction, but a banner that changes shape
    is a real change and should be visible as its own artefact.

    Args:
        path: The TSV.
        rules: The normalisation rules for this side.

    Returns:
        dict[str, Any] | None: ``columns``, ``rows`` and ``provenance``, or None if the file
        does not exist.
    """
    if not path.is_file():
        return None
    provenance: list[str] = []
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.rstrip("\r")
        if not line.strip():
            continue
        if line.startswith("#"):
            provenance.append(line)
            continue
        if header is None:
            header = line.split("\t")
            continue
        values = line.split("\t")
        if len(values) < len(header):
            values = values + [""] * (len(header) - len(values))
        rows.append(dict(zip(header, values[: len(header)], strict=True)))
    return {
        "columns": header or [],
        "rows": normalise.apply_deep(rows, rules),
        "provenance": normalise.apply_deep(provenance, rules),
    }


def read_json(path: Path) -> Any | None:
    """Read a JSON file, returning None when it is absent or unparseable.

    Args:
        path: The JSON file.

    Returns:
        Any | None: The parsed value, or None.
    """
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        logger.warning(f"Could not parse {path}: {exc}")
        return None


def read_delimited(path: Path, sep: str, rules: list[Rule], sort_key: str | None = None) -> dict[str, Any] | None:
    """Read a cohort CSV/TSV export into columns and rows.

    Args:
        path: The export.
        sep: The field separator.
        rules: The normalisation rules for this side.
        sort_key: A column to sort rows by, because cohort sample order is not reproducible.

    Returns:
        dict[str, Any] | None: ``columns`` and ``rows``, or None if absent.
    """
    if not path.is_file():
        return None
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        return {"columns": [], "rows": []}
    header = lines[0].split(sep)
    rows = []
    for line in lines[1:]:
        values = line.split(sep)
        if len(values) < len(header):
            values = values + [""] * (len(header) - len(values))
        rows.append(dict(zip(header, values[: len(header)], strict=True)))
    rows = normalise.apply_deep(rows, rules)
    if sort_key and sort_key in header:
        rows.sort(key=lambda row: (row.get(sort_key, ""), json.dumps(row, sort_keys=True)))
    return {"columns": header, "rows": rows}


def _summary_boxes(html: str, rules: list[Rule]) -> list[dict[str, Any]]:
    """Extract each ``summary-box`` paragraph's text and computed emphasis.

    Args:
        html: The rendered per-sample report.
        rules: The normalisation rules for this side.

    Returns:
        list[dict[str, Any]]: ``text`` and ``is_positive``, in document order.
    """
    boxes = []
    for classes, text in SUMMARY_BOX_RE.findall(html):
        boxes.append(
            {
                "text": normalise.apply(" ".join(text.split()), rules),
                "is_positive": "summary-positive" in classes,
            }
        )
    return boxes


def read_report(path: Path, rules: list[Rule]) -> dict[str, Any] | None:
    """Read the per-sample HTML report's summary boxes and rendered tables.

    The tables are read because three of run 3's nine adjudicated presentation deltas live
    in them - the ``Motif`` cell showing the annotated motif rather than the raw pair, the
    ``Motif`` column header, and the ``None`` cell added to the negative-case table - and a
    gate that reads only the TSVs cannot see any of them. Row order is left alone: unlike
    the cohort's, the per-sample report's ordering is deterministic.

    Args:
        path: ``summary_report.html``.
        rules: The normalisation rules for this side.

    Returns:
        dict[str, Any] | None: ``screening``, ``cross_match`` and ``tables``, or None if the
        report was not produced.
    """
    if not path.is_file():
        return None
    html = path.read_text(encoding="utf-8", errors="replace")
    boxes = _summary_boxes(html, rules)
    extractor = _TableExtractor()
    extractor.feed(html)
    return {
        "screening": boxes[0] if boxes else None,
        "cross_match": boxes[1] if len(boxes) > 1 else None,
        "box_count": len(boxes),
        "tables": normalise.apply_deep(extractor.tables, rules),
    }


def read_pipeline_case(output_dir: Path, log_dir: Path, rules: list[Rule]) -> dict[str, Any]:
    """Read every artefact the gate compares for one per-sample run.

    Args:
        output_dir: The run's ``--output-dir``.
        log_dir: Where the harness put the run's launch line, logs and command record.
        rules: The normalisation rules for this side.

    Returns:
        dict[str, Any]: Artefact name -> value. A ``None`` value means the artefact was not
        produced, which is itself compared.
    """
    result = read_json(log_dir / "result.json") or {}
    summary = read_json(output_dir / "pipeline_summary.json")
    steps = summary.get("steps", []) if isinstance(summary, dict) else []

    kestrel = read_tsv(output_dir / "kestrel" / "kestrel_result.tsv", rules)
    kestrel_pre = read_tsv(output_dir / "kestrel" / "kestrel_pre_result.tsv", rules)
    advntr = read_tsv(output_dir / "advntr" / "output_adVNTR_result.tsv", rules)
    coverage = read_tsv(output_dir / "coverage" / "coverage_summary.tsv", rules)
    report = read_report(output_dir / "summary_report.html", rules)

    commands: list[str] = []
    commands_log = log_dir / "commands.jsonl"
    if commands_log.is_file():
        commands = [
            normalise.apply(json.loads(line)["command"], rules)
            for line in commands_log.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]

    step_records = []
    for step in steps:
        record = {key: value for key, value in step.items() if key not in normalise.DROPPED_KEYS}
        record.pop("parsed_result", None)
        step_records.append(normalise.apply_deep(record, rules))

    return {
        "exit_code": result.get("exit_code"),
        "launch_line": normalise.apply(result.get("launch_line") or "", rules),
        "kestrel_result": kestrel,
        "kestrel_pre_result": kestrel_pre,
        "advntr_result": advntr,
        "coverage_summary": coverage,
        "screening_summary": (report or {}).get("screening"),
        "cross_match_summary": (report or {}).get("cross_match"),
        "report_tables": (report or {}).get("tables"),
        "pipeline_steps": [step.get("step") for step in steps],
        "pipeline_step_records": step_records,
        "executed_commands": commands,
    }


def read_cohort_case(output_dir: Path, log_dir: Path, rules: list[Rule]) -> dict[str, Any]:
    """Read every artefact the gate compares for one ``vntyper cohort`` run.

    Args:
        output_dir: The cohort run's ``--output-dir``.
        log_dir: Where the harness put the run's launch line, logs and command record.
        rules: The normalisation rules for this side.

    Returns:
        dict[str, Any]: Artefact name -> value.
    """
    result = read_json(log_dir / "result.json") or {}
    html_path = output_dir / "cohort_summary.html"
    tables: list[dict[str, Any]] | None = None
    donut_values: list[list[str]] | None = None
    donut_totals: list[str] | None = None
    sample_order: list[str] | None = None

    if html_path.is_file():
        html = html_path.read_text(encoding="utf-8", errors="replace")
        extractor = _TableExtractor()
        extractor.feed(html)
        tables = []
        for table in extractor.tables:
            rows = normalise.apply_deep(table["rows"], rules)
            if sample_order is None and table["header"] and table["header"][0] == "Sample":
                sample_order = [row[0] for row in rows if row]
            rows.sort()
            tables.append({"header": normalise.apply_deep(table["header"], rules), "rows": rows})
        donut_values = [value.strip() for value in DONUT_VALUES_RE.findall(html)]
        donut_totals = _donut_totals(html)

    listing = None
    if output_dir.is_dir():
        listing = sorted(str(path.relative_to(output_dir)) for path in output_dir.rglob("*") if path.is_file())

    return {
        "exit_code": result.get("exit_code"),
        "launch_line": normalise.apply(result.get("launch_line") or "", rules),
        "cohort_tables": tables,
        "cohort_category_counts": donut_values,
        "cohort_category_totals": donut_totals,
        "cohort_kestrel_csv": read_delimited(output_dir / "cohort_kestrel.csv", ",", rules, "Sample"),
        "cohort_kestrel_tsv": read_delimited(output_dir / "cohort_kestrel.tsv", "\t", rules, "Sample"),
        "cohort_kestrel_json": _sorted_records(read_json(output_dir / "cohort_kestrel.json"), rules),
        "cohort_advntr_csv": read_delimited(output_dir / "cohort_advntr.csv", ",", rules, "Sample"),
        "cohort_advntr_tsv": read_delimited(output_dir / "cohort_advntr.tsv", "\t", rules, "Sample"),
        "cohort_advntr_json": _sorted_records(read_json(output_dir / "cohort_advntr.json"), rules),
        "pseudonymization_table": read_delimited(output_dir / "pseudonymization_table.tsv", "\t", rules, "Pseudonym"),
        "cohort_output_files": listing,
        "cohort_sample_order_raw": sample_order,
    }


def _donut_totals(html: str) -> list[str]:
    """Read the donut charts' centre totals out of the embedded plotly figure JSON.

    Args:
        html: The rendered cohort report.

    Returns:
        list[str]: One total per chart, in document order.
    """
    totals: list[str] = []
    for literal in DONUT_TEXT_RE.findall(html):
        try:
            decoded = json.loads(literal)
        except ValueError:
            continue
        totals.extend(DONUT_TOTAL_RE.findall(decoded))
    return totals


def _sorted_records(records: Any, rules: list[Rule]) -> Any:
    """Normalise a list of JSON records and sort it, cohort order not being reproducible.

    Args:
        records: The parsed JSON, expected to be a list of mappings.
        rules: The normalisation rules for this side.

    Returns:
        Any: The normalised records, sorted, or the value unchanged if it is not a list.
    """
    if not isinstance(records, list):
        return records
    normalised = normalise.apply_deep(records, rules)
    return sorted(normalised, key=lambda record: json.dumps(record, sort_keys=True))
