"""The behaviour oracle for the cohort report.

`cohort_summary.py` was split into focused modules (`cohort_rules`, `cohort_categories`,
`cohort_tables`, `cohort_inputs`, `cohort_exports`) as a **pure refactor**: every extracted
function had to produce byte-identical output to the code it replaced. This file is the
instrument that says so.

Why it is load-bearing rather than ceremonial: the golden-cohort gate
(`docs/development/golden-cohort-gate.md`) runs **no cohort-mode case at all**, so it
produces no delta for any change to this file. This oracle and
`tests/unit/test_cohort_summary_escaping.py` are the only two instruments on the cohort
report. If a refactor makes either of them need editing, the refactor changed behaviour.

What the fingerprint pins
-------------------------
* The three rendered tables, byte for byte - column selection, column order, cell values,
  the pandas scaffolding, and the per-column escaping exemptions.
* Both donut charts' data: the three segment values, their labels, and the total in the
  centre. That is where the sample-level categorisation (Positive / Positive_Flagged /
  Negative) becomes observable.
* The page skeleton the Jinja2 template produces around them, so a fragment routed into
  the wrong context slot shows up.

What it deliberately excludes, and why
--------------------------------------
* The embedded plotly.js bundle (4.8 MB of it, twice) and the base64 PNG payloads. Both
  are byte-stable for a given matplotlib/plotly build but change on a dependency bump,
  and a hash that turns red on `pip install -U plotly` stops being read. The chart *data*
  is extracted separately and is pinned exactly, so a miscounted cohort still fails.
* The report timestamp and plotly's per-figure UUIDs, which are new on every render.

What this oracle does **not** cover
-----------------------------------
* Zip-file inputs. `aggregate_cohort` extracts them to a temporary directory; the
  discovery half of that is pinned in `tests/unit/test_cohort_inputs.py` instead.
* The order samples appear in. `aggregate_cohort` collects sample directories into a
  `set`, so the row order of both the report and the CSV/TSV/JSON exports varies between
  processes. That is pinned as a defect, not as intent - see
  `test_cohort_inputs.py::test_discovery_returns_an_unordered_set_today`.
* Anything the pipeline writes *into* `pipeline_summary.json`. The oracle starts from a
  summary file, so a change in what `pipeline.py` records is invisible here.
"""

from __future__ import annotations

import hashlib
import json
import logging
import re
from pathlib import Path

import pandas as pd
import pytest

from vntyper.cli import load_config
from vntyper.scripts import cohort_summary
from vntyper.scripts.summary_steps import STEP_ADVNTR, STEP_BAM_HEADER, STEP_KESTREL

pytestmark = pytest.mark.unit

#: The recorded fingerprint of the two-sample cohort report below. A refactor that
#: changes this changed the report; that is the whole point of the number.
EXPECTED_FINGERPRINT = "d82eb3745e5a8f1f118659d8e5853492ad1b5b7b2b424be11a54e1c79c1c28ee"

_UUID = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")
_TIMESTAMP = re.compile(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}")
_SCRIPT_BODY = re.compile(r"(<script\b[^>]*>).*?(</script>)", re.DOTALL)
_BASE64_PNG = re.compile(r"data:image/png;base64,[A-Za-z0-9+/=]*")
_TABLE = re.compile(r"<table\b.*?</table>", re.DOTALL)
_CHART_VALUES = re.compile(r'"values":\s*\[[^\]]*\]')
_CHART_LABELS = re.compile(r'"labels":\s*\[[^\]]*\]')
#: Plotly serialises the centre annotation with its angle brackets JSON-escaped.
_CHART_TOTAL = re.compile(r'"text":"\\u003cb\\u003e(\d+)\\u003c\\u002fb\\u003e"')


def _skeleton(html: str) -> str:
    """Strip every version-volatile fragment out of the page.

    Args:
        html: The rendered cohort report.

    Returns:
        str: The page with script bodies, base64 payloads, UUIDs and the report
        timestamp replaced by fixed markers.
    """
    html = _SCRIPT_BODY.sub(r"\1<SCRIPT-BODY>\2", html)
    html = _BASE64_PNG.sub("data:image/png;base64,<PNG>", html)
    html = _UUID.sub("<UUID>", html)
    return _TIMESTAMP.sub("<TIMESTAMP>", html)


def _fingerprint(html: str) -> tuple[str, str]:
    """Reduce a rendered report to a canonical document and its digest.

    Args:
        html: The rendered cohort report.

    Returns:
        tuple[str, str]: The canonical document (for diffing on failure) and its
        SHA-256 hex digest.
    """
    parts = [
        "[TABLES]",
        *_TABLE.findall(html),
        "[CHART-VALUES]",
        *_CHART_VALUES.findall(html),
        "[CHART-LABELS]",
        *_CHART_LABELS.findall(html),
        "[CHART-TOTALS]",
        *_CHART_TOTAL.findall(html),
        "[IMAGES]",
        f"count={len(_BASE64_PNG.findall(html))}",
        f"nonempty={sum(1 for src in _BASE64_PNG.findall(html) if len(src) > 64)}",
        "[SKELETON]",
        _skeleton(html),
    ]
    document = "\n".join(parts)
    return document, hashlib.sha256(document.encode("utf-8")).hexdigest()


def _kestrel_frame() -> pd.DataFrame:
    """Build the Kestrel half of the fixture cohort.

    Three rows over two samples: one clean high-precision call, one flagged
    low-precision call, and one row whose ``Confidence`` matches no rule, so all three
    sample-level categories and both colour branches are exercised. ``Del`` is present
    and is not a display column, which pins the column selection.

    Several cells carry a literal ``&``. It is not decoration: an ``&`` renders as
    ``&amp;`` when the column is escaped and as itself when it is not, so widening the
    ``html_columns`` exemption to any of these columns moves the fingerprint. The
    ``&`` in ``Flag`` and ``Confidence`` leaves both rule outcomes unchanged - those
    two values match no rule either way.

    Returns:
        pd.DataFrame: A fresh frame - the renderer mutates the one it is handed.
    """
    return pd.DataFrame(
        [
            {
                "Sample": "sample_one",
                "Motif": "5",
                "Variant": "insC",
                "POS": "155188205",
                "REF": "G",
                "ALT": "GC",
                "Motif_sequence": "ACGT&ACGT",
                "Estimated_Depth_AlternateVariant": "12",
                "Estimated_Depth_Variant_ActiveRegion": "40",
                "Depth_Score": "0.3",
                "Confidence": "High_Precision",
                "Flag": "Not flagged",
                "Del": "dropped",
            },
            {
                "Sample": "sample_two",
                "Motif": "8",
                "Variant": "delG",
                "POS": "155188300",
                "REF": "GG",
                "ALT": "G",
                "Motif_sequence": "TTTT",
                "Estimated_Depth_AlternateVariant": "3",
                "Estimated_Depth_Variant_ActiveRegion": "9",
                "Depth_Score": "0.33",
                "Confidence": "Low_Precision",
                "Flag": "Depth_Score_flagged&",
                "Del": "dropped",
            },
            {
                "Sample": "sample_two",
                "Motif": "9&",
                "Variant": "insT&",
                "POS": "155188400",
                "REF": "A&",
                "ALT": "AT&",
                "Motif_sequence": "GGGG&",
                "Estimated_Depth_AlternateVariant": "1",
                "Estimated_Depth_Variant_ActiveRegion": "2",
                "Depth_Score": "0.5",
                "Confidence": "Unrecognised&",
                "Flag": "Not flagged",
                "Del": "dropped",
            },
        ]
    )


def _advntr_frame() -> pd.DataFrame:
    """Build the adVNTR half of the fixture cohort.

    The adVNTR table has no escaping exemption at all, so the first row marks every
    display column with an ``&`` for the same reason the Kestrel fixture does.
    ``25561&`` is still ``!= "Negative"``, so the rule outcome is unchanged.

    Returns:
        pd.DataFrame: A fresh frame - the renderer mutates the one it is handed.
    """
    return pd.DataFrame(
        [
            {
                "Sample": "sample_one",
                "VID": "25561&",
                "Variant": "I22_G_GG&",
                "NumberOfSupportingReads": "9&",
                "MeanCoverage": "30&",
                "Pvalue": "0.001&",
                "RU": "GG&",
                "POS": "155188205&",
                "REF": "G&",
                "ALT": "GG&",
                "Flag": "Not flagged",
                "Unlisted": "dropped",
            },
            {
                "Sample": "sample_two",
                "VID": "Negative",
                "Variant": "Negative",
                "NumberOfSupportingReads": "0",
                "MeanCoverage": "12",
                "Pvalue": "1.0",
                "RU": "None",
                "POS": "None",
                "REF": "None",
                "ALT": "None",
                "Flag": "Coverage_flagged",
                "Unlisted": "dropped",
            },
        ]
    )


def _stats_frame() -> pd.DataFrame:
    """Build the per-sample statistics half of the fixture cohort.

    Every cell here is read straight out of a sample's ``pipeline_summary.json``, so
    this table has no exemption either; the ``&`` marks say so.

    Returns:
        pd.DataFrame: One row per sample, ``Sample`` first.
    """
    return pd.DataFrame(
        [
            {
                "Sample": "sample_one",
                "runtime": "90.00 seconds",
                "version": "2.0.6&",
                "assembly": "hg38&",
                "pipeline": "bwa-mem&",
                "cov_mean": "31.2",
            },
            {
                "Sample": "sample_two",
                "runtime": "N/A",
                "version": "2.0.6",
                "assembly": "hg19",
                "pipeline": "N/A",
                "cov_mean": "12.0",
            },
        ]
    )


def _render(tmp_path: Path) -> str:
    """Render the fixture cohort and return the page.

    Args:
        tmp_path: The pytest per-test directory.

    Returns:
        str: The rendered HTML.
    """
    cohort_summary.generate_cohort_summary_report(
        output_dir=str(tmp_path),
        kestrel_df=_kestrel_frame(),
        advntr_df=_advntr_frame(),
        summary_file="cohort_summary.html",
        config=load_config(None),
        additional_stats_html=cohort_summary.stats_table_html(_stats_frame()),
    )
    return (tmp_path / "cohort_summary.html").read_text(encoding="utf-8")


def test_the_cohort_report_matches_its_recorded_fingerprint(tmp_path) -> None:
    """The whole report, reduced to a digest. Characterisation, not specification.

    This records what the cohort report produced before the split and makes no claim
    that any of it is right. It exists so that a refactor which changes a single cell,
    drops a column, reorders one, or miscounts a donut segment cannot land silently.

    On failure the canonical document is written next to the report so the change can be
    diffed rather than guessed at.
    """
    document, digest = _fingerprint(_render(tmp_path))

    if digest != EXPECTED_FINGERPRINT:
        dump = tmp_path / "fingerprint.txt"
        dump.write_text(document, encoding="utf-8")
        pytest.fail(
            f"the cohort report changed: {digest} != {EXPECTED_FINGERPRINT}. "
            f"The canonical document is at {dump}. If this was a refactor, the refactor "
            "was not behaviour-preserving."
        )


# ---------------------------------------------------------------------------
# The same facts stated readably, so a fingerprint failure can be localised
# ---------------------------------------------------------------------------


def test_the_kestrel_table_keeps_its_column_order(tmp_path) -> None:
    """Twelve display columns in a fixed order; `Del` is not one of them."""
    html = _render(tmp_path)
    table = _TABLE.findall(html)[0]

    headings = re.findall(r"<th>(.*?)</th>", table)

    assert headings == [
        "Sample",
        "Motif",
        "Variant",
        "POS",
        "REF",
        "ALT",
        "Motif_sequence",
        "Estimated_Depth_AlternateVariant",
        "Estimated_Depth_Variant_ActiveRegion",
        "Depth_Score",
        "Confidence",
        "Flag",
    ]


def test_the_advntr_table_keeps_its_column_order(tmp_path) -> None:
    """Eleven display columns in a fixed order; `Unlisted` is not one of them."""
    html = _render(tmp_path)
    table = _TABLE.findall(html)[1]

    headings = re.findall(r"<th>(.*?)</th>", table)

    assert headings == [
        "Sample",
        "VID",
        "Variant",
        "NumberOfSupportingReads",
        "MeanCoverage",
        "Pvalue",
        "RU",
        "POS",
        "REF",
        "ALT",
        "Flag",
    ]


def test_the_donut_charts_carry_the_sample_level_counts(tmp_path) -> None:
    """Where the categorisation becomes visible.

    Kestrel: `sample_one` is a clean high-precision call (Positive) and `sample_two`
    has a flagged low-precision call plus an unrecognised one, so it aggregates to
    Positive_Flagged - one Positive, one Positive_Flagged, no Negative.

    adVNTR: `sample_one` is positive and `sample_two` is flagged.
    """
    html = _render(tmp_path)

    assert _CHART_VALUES.findall(html) == ['"values":[1,1,0]', '"values":[1,1,0]']
    assert _CHART_TOTAL.findall(html) == ["2", "2"]


def test_the_report_is_written_under_the_output_directory(tmp_path) -> None:
    """`--summary-file` is a name, and the plots land beside it."""
    _render(tmp_path)

    assert (tmp_path / "cohort_summary.html").is_file()
    assert (tmp_path / "plots" / "kestrel_summary_plot.png").is_file()
    assert (tmp_path / "plots" / "advntr_summary_plot.png").is_file()


def test_the_renderer_annotates_the_frame_it_was_handed_today(tmp_path) -> None:
    """Characterisation of a live defect: the input frame is mutated in place.

    `generate_cohort_summary_report` writes `__row_result` and `__unified` onto the
    caller's DataFrame rather than onto its own copy. `aggregate_cohort` renders the
    report *before* it exports CSV/TSV/JSON, so both internal columns reach every
    machine-readable output - which the code comment two lines above the copy says is
    exactly what the copy exists to prevent.

    Pinned rather than fixed: the split that this oracle guards is not authorised to
    change what any output contains. See
    `.superpowers/sdd/2026-08-06-issue-181-197-followups-plan/issue-cohort-internal-columns-leak-into-exports.md`.
    """
    frame = _kestrel_frame()

    cohort_summary.generate_cohort_summary_report(
        output_dir=str(tmp_path),
        kestrel_df=frame,
        advntr_df=_advntr_frame(),
        summary_file="cohort_summary.html",
        config=load_config(None),
    )

    assert "__row_result" in frame.columns
    assert "__unified" in frame.columns
    assert list(frame["__unified"]) == ["Positive", "Positive_Flagged", "Negative"]


# ---------------------------------------------------------------------------
# `aggregate_cohort`, the seam the five extracted modules are wired back into
# ---------------------------------------------------------------------------


def _cohort_on_disk(root: Path, *, advntr: bool = True) -> Path:
    """Write a two-sample cohort under ``root``.

    Args:
        root: The directory to create the samples in.
        advntr: Whether the first sample recorded an adVNTR step.

    Returns:
        Path: ``root``, ready to be passed to `aggregate_cohort`.
    """
    steps_one = [
        {
            "step": STEP_KESTREL,
            "parsed_result": {"data": [{"Motif": "5", "Confidence": "High_Precision", "Flag": "Not flagged"}]},
        },
        {"step": STEP_BAM_HEADER, "parsed_result": {"assembly_text": "hg38", "alignment_pipeline": "bwa-mem"}},
    ]
    if advntr:
        steps_one.append(
            {"step": STEP_ADVNTR, "parsed_result": {"data": [{"VID": "25561", "Flag": "Not flagged"}]}},
        )
    (root / "sample_one").mkdir(parents=True)
    (root / "sample_one" / "pipeline_summary.json").write_text(
        json.dumps({"version": "2.0.6", "steps": steps_one}), encoding="utf-8"
    )
    (root / "sample_two").mkdir(parents=True)
    (root / "sample_two" / "pipeline_summary.json").write_text(
        json.dumps(
            {
                "version": "2.0.6",
                "steps": [
                    {
                        "step": STEP_KESTREL,
                        "parsed_result": {"data": [{"Motif": "8", "Confidence": "Low_Precision", "Flag": "flagged"}]},
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    return root


def test_a_cohort_run_writes_the_report_and_every_requested_export(tmp_path) -> None:
    """The end-to-end path through all five extracted modules at once."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(_cohort_on_disk(tmp_path / "cohort"))],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=load_config(None),
        additional_formats="csv,tsv,json",
    )

    written = {p.name for p in output_dir.iterdir()}

    assert written == {
        "cohort_summary.html",
        "plots",
        "cohort_kestrel.csv",
        "cohort_kestrel.tsv",
        "cohort_kestrel.json",
        "cohort_advntr.csv",
        "cohort_advntr.tsv",
        "cohort_advntr.json",
    }


def test_both_samples_reach_the_report(tmp_path) -> None:
    """The frames are built from the discovered directories, in whatever order the
    discovery yields them, so this asserts on membership rather than on order."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(_cohort_on_disk(tmp_path / "cohort"))],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=load_config(None),
    )

    html = (output_dir / "cohort_summary.html").read_text(encoding="utf-8")

    assert ">sample_one</td>" in html
    assert ">sample_two</td>" in html
    assert '"values":[1,1,0]' in html


def test_an_algorithm_no_sample_ran_produces_no_export_for_it(tmp_path, caplog) -> None:
    caplog.set_level(logging.WARNING, logger="vntyper.scripts.cohort_summary")
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(_cohort_on_disk(tmp_path / "cohort", advntr=False))],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=load_config(None),
        additional_formats="csv",
    )

    assert not (output_dir / "cohort_advntr.csv").exists()
    assert (output_dir / "cohort_kestrel.csv").exists()
    assert "No adVNTR data found in any sample." in caplog.text


def test_a_cohort_with_no_valid_input_writes_nothing_at_all(tmp_path, caplog) -> None:
    """`aggregate_cohort` returns before the render, so there is no empty report and no
    traceback either."""
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cohort_summary")
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(tmp_path / "absent")],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=load_config(None),
    )

    assert list(output_dir.iterdir()) == []
    assert "No valid input directories or zip files found" in caplog.text


def test_pseudonymised_samples_reach_the_report_under_their_pseudonyms(tmp_path) -> None:
    """The sample name is the one value in the report a caller fully controls, and
    through the web service it is the pseudonym recorded against an uploaded job."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(_cohort_on_disk(tmp_path / "cohort"))],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=load_config(None),
        pseudonymize_samples="anon_",
    )

    html = (output_dir / "cohort_summary.html").read_text(encoding="utf-8")
    table = (output_dir / "pseudonymization_table.tsv").read_text(encoding="utf-8")

    assert ">sample_one</td>" not in html
    assert ">anon_65622</td>" in html
    assert "anon_65622\tsample_one" in table


def test_the_pseudonymization_table_is_not_written_when_nothing_was_pseudonymised(tmp_path) -> None:
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(_cohort_on_disk(tmp_path / "cohort"))],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=load_config(None),
    )

    assert not (output_dir / "pseudonymization_table.tsv").exists()


def test_an_export_written_after_the_report_carries_the_internal_columns_today(tmp_path) -> None:
    """The consequence of the in-place annotation above, at the boundary a user sees.

    A cohort run with `--output-format csv` produces a `cohort_kestrel.csv` whose
    header carries `__row_result` and `__unified`. Characterisation of the same defect;
    see the issue file named in the test above.
    """
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(_cohort_on_disk(tmp_path / "cohort"))],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=load_config(None),
        additional_formats="csv",
    )

    header = (output_dir / "cohort_kestrel.csv").read_text(encoding="utf-8").splitlines()[0]

    assert "__row_result" in header
    assert "__unified" in header
