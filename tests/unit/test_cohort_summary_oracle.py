"""The behaviour oracle for the cohort report.

`cohort_summary.py` was split into focused modules (`cohort_rules`, `cohort_categories`,
`cohort_tables`, `cohort_inputs`, `cohort_exports`) as a **pure refactor**: every extracted
function had to produce byte-identical output to the code it replaced. This file is the
instrument that says so.

Why it is load-bearing rather than ceremonial: it is the only instrument that pins the
whole assembled page at once, and the only one that does so against a recorded constant.

The golden-cohort gate (`docs/development/golden-cohort-gate.md`) **does** now run cohort
mode - `scripts/golden_cohort/matrix.py::build_cohort_cases` defines four cases,
`cohort_multi`, `cohort_multi_pseudonymized`, `cohort_single` and `cohort_empty`, and run 4
of the gate is the first in this project's history to cover cohort mode at all. Two limits
on what that buys, both recorded on the gate page: it is a before-versus-after comparison
of two commits, run by hand and wired into no CI job or `make` target, so it produces
nothing until someone runs it and it attests one candidate commit and nothing after it;
and it **sorts** cohort rows before comparing them, so sample ordering is normalised away
there and is attested only by unit tests.

The per-module suites - `test_cohort_tables.py`, `test_cohort_categories.py`,
`test_cohort_rules.py`, `test_cohort_inputs.py`, `test_cohort_exports.py` - and
`test_cohort_summary_escaping.py` each cover their own seam. What none of them covers is
the assembled page: which fragment reached which slot of the template, and what the whole
thing hashes to. That is this file. If a refactor makes it need editing, the refactor
changed behaviour.

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
**Read this list as the limit of what a green run here proves.** It is not exhaustive by
construction - it is exhaustive as far as anyone has looked - so treat it as the floor of
the oracle's blind spots rather than the ceiling.

* **Decision-logic paths the fixture never exercises. This is the important one.**
  The fixture below produces 3 of the 5 Kestrel verdicts (`High_Precision`,
  `Low_Precision_flagged`, and the `negative` default) and 2 of the 3 adVNTR verdicts
  (`positive`, `positive flagged`). It never produces a bare `Low_Precision`, a
  `High_Precision_flagged`, or an adVNTR `negative`. A change to how an unexercised
  verdict is categorised therefore **passes this oracle and passes the escaping suite** -
  demonstrated: mutating `unify_kestrel_result` so `Low_Precision` maps to `Negative`
  instead of `Positive`, which would move a sample out of the positive segment of every
  cohort donut chart in the field, is caught only by
  `tests/unit/test_cohort_categories.py`.

  So state the guarantee precisely: **this oracle proves equivalence over the paths the
  fixture covers, plus - for the whole-function moves it was built to attest - byte
  identity of the moved bodies.** A branch the fixture never reaches cannot have been
  mis-transcribed if its body is unchanged, which is why the split it guards is sound
  despite the gap. It is *not* blanket behavioural equivalence, and it must not be leant
  on for a change that is anything other than a pure move. For those, the per-module
  suites are the instrument.

  Do **not** close this gap by adding fixture rows. It would churn the fingerprint for no
  gain in what the oracle is actually for, and the categorisation paths already have
  exhaustive direct coverage in `test_cohort_categories.py` and `test_cohort_rules.py`.

  **Worked example, from this file's own history.** Which verdicts the fixture happens to
  reach is not a detail - it decides whether a rule change is visible here at all. When
  `report_config.json`'s adVNTR rule 2 gained the `VID != "Negative"` guard that rule 1
  always had, the change was genotype-neutral across all five `(VID, Flag)` pairs the
  pipeline can produce, and it still broke this fingerprint - because the fixture's
  `sample_two` row was `VID="Negative"` with `Flag="Coverage_flagged"`, a state the
  pipeline cannot produce and the one input the guard did move. The reverse holds just as
  well: had the fixture used a different adVNTR row, a rule change that *was* observable
  could have gone unseen. The row is now a reachable one, and the comment on it says so.

* Zip-file inputs. `aggregate_cohort` extracts them to a temporary directory; the
  discovery half of that is pinned in `tests/unit/test_cohort_inputs.py` instead.
* The order samples appear in - as far as *this fingerprint* goes. It is taken from
  `generate_cohort_summary_report`, which is handed frames directly, so the discovery
  that decides the row order never runs. For **directory** inputs the order is now
  deterministic and is pinned at its source in
  `test_cohort_inputs.py::test_the_discovered_directories_come_back_sorted` and, end to
  end through `aggregate_cohort`, by `test_both_samples_reach_the_report_in_sorted_order`
  below. For **zip** inputs it is not: each archive extracts to a fresh
  `tempfile.mkdtemp(prefix="cohort_zip_")` whose random component is part of the sort key,
  so the order of two zips is not reproducible. That is characterised, not fixed, in
  `test_cohort_inputs.py::test_two_zip_inputs_are_ordered_by_their_random_temporary_directories_today`.
* Anything the pipeline writes *into* `pipeline_summary.json`. The oracle starts from a
  summary file, so a change in what `pipeline.py` records is invisible here.
* Export column order. The fingerprint is taken from the HTML only; the CSV/TSV/JSON
  headers are a different order entirely (`Sample` is twelfth, not first) and are pinned
  in `tests/unit/test_cohort_exports.py`.
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
#:
#: It has moved once, and the reason is recorded here because a changed fingerprint with
#: no explanation should be read as the worst case:
#:
#: * **Old**: ``d82eb3745e5a8f1f118659d8e5853492ad1b5b7b2b424be11a54e1c79c1c28ee``
#: * **New**: ``9889773ac381a6d0f33c2394c1f3d4f6a795cbc5bb38c5cbc9773f2e3a615645``
#:
#: **Cause: a fixture correction, not a behaviour change.** `sample_two`'s adVNTR row was
#: `VID="Negative"` with `Flag="Coverage_flagged"` - a state the pipeline cannot produce
#: (the only producer of `VID == "Negative"` sets `Flag = "Not applicable"` and never
#: reaches `add_flags`; `Coverage_flagged` is not a flag `add_flags` can emit at all). It
#: is now a reachable row: a real `VID="25561"` call flagged `Low_Coverage`. The trigger
#: was `report_config.json`'s adVNTR rule 2 gaining the `VID != "Negative"` guard that
#: rule 1 always had - genotype-neutral on all five `(VID, Flag)` pairs the pipeline can
#: produce, and it moved only that unreachable row. See the comment on the row itself.
#:
#: **The refactor attestation does not rest on this constant and is not weakened by the
#: re-derivation.** The old hash was independently verified against the *pre-split*
#: 911-line `cohort_summary.py`: the reviewer reverted the file to its original state and
#: ran this oracle unmodified, 27 passed. That proof is banked and stands on its own -
#: `d82eb374...` is what the original implementation produced, and the five extracted
#: modules reproduced it byte for byte.
#:
#: **Re-derived, and it did not move, when the two cohort defects were fixed** (the
#: in-place annotation and the non-deterministic sample order). It was re-derived twice
#: under `PYTHONHASHSEED=0` and `PYTHONHASHSEED=12345` in separate interpreters and both
#: produced `9889773a...`, which is the value already recorded above. Each fix was also
#: applied on its own and neither moved it. That is the expected result rather than a
#: surprise, and it says something worth keeping:
#:
#: * The internal columns never reached the *page*. `cohort_tables` builds its tables
#:   from the fixed `KESTREL_DISPLAY_COLUMNS` / `ADVNTR_DISPLAY_COLUMNS` lists, so
#:   `__row_result` and `__unified` were dropped there. The leak was into the
#:   CSV/TSV/JSON exports only, and this fingerprint is taken from the HTML - so the fix
#:   is pinned by `test_an_export_written_after_the_report_carries_no_internal_columns`
#:   and by `test_the_renderer_leaves_the_frames_it_was_handed_untouched`, not here.
#: * The sample order never reached this fingerprint either: it is measured from
#:   `generate_cohort_summary_report`, which is handed frames the test built, so
#:   `discover_sample_directories` does not run. See the exclusion list above.
#:
#: **What the two-seed re-derivation establishes, and what it does not.** It establishes
#: that this constant does not itself depend on the process hash seed - worth knowing,
#: because a fingerprint that moved with the seed would be useless as an instrument. It is
#: **not** evidence for the ordering fix, and an earlier version of this comment claimed
#: it was. The fingerprint calls `generate_cohort_summary_report` directly and never runs
#: `discover_sample_directories`, which is the only code path the ordering fix touched, so
#: it would have been stable across seeds before the fix as well - the two-seed result is
#: consistent with the fix and with its absence alike, which is what makes it evidence for
#: neither. The evidence for that fix is
#: `test_cohort_inputs.py::test_processes_with_different_hash_seeds_discover_the_same_order`,
#: which spawns five interpreters under five seeds, and it covers directory inputs only.
EXPECTED_FINGERPRINT = "9889773ac381a6d0f33c2394c1f3d4f6a795cbc5bb38c5cbc9773f2e3a615645"

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
        pd.DataFrame: A fresh frame each call, so the equality check in
        `test_the_renderer_leaves_the_frames_it_was_handed_untouched` has something
        independent to compare against.
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
        pd.DataFrame: A fresh frame each call, for the same reason as
        :func:`_kestrel_frame`.
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
            # A real call that the shipped flagging rules flagged. Every value here is
            # one the pipeline can actually produce: `VID` is the hardcoded MUC1 VNTR
            # locus id, and `Low_Coverage` is a flag name read straight off
            # `advntr_config.json`'s `flagging_rules`, earned by the four supporting
            # reads (its rule is `NumberOfSupportingReads < 10`).
            #
            # This row used to be `VID="Negative"` with `Flag="Coverage_flagged"` - a
            # state the pipeline cannot produce twice over: the only producer of a
            # `VID == "Negative"` row sets `Flag = "Not applicable"` and never reaches
            # `add_flags`, and `Coverage_flagged` is not a string `add_flags` can emit
            # at all. Modelling an impossible state held this fingerprint hostage to a
            # rule about impossible inputs, and it duly broke when that rule was
            # corrected. Keep this row reachable.
            {
                "Sample": "sample_two",
                "VID": "25561",
                "Variant": "I22_G_GG",
                "NumberOfSupportingReads": "4",
                "MeanCoverage": "12",
                "Pvalue": "0.04",
                "RU": "GG",
                "POS": "155188205",
                "REF": "G",
                "ALT": "GG",
                "Flag": "Low_Coverage",
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

    adVNTR: `sample_one` is a clean call (Positive) and `sample_two` is a real call
    carrying the `Low_Coverage` flag (Positive_Flagged) - the same one/one/zero split.
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


def test_the_renderer_leaves_the_frames_it_was_handed_untouched(tmp_path) -> None:
    """**Specification**: rendering reads its inputs and writes only the report.

    `generate_cohort_summary_report` used to write `__row_result` and `__unified` onto
    the caller's DataFrames rather than onto its own copies, and `aggregate_cohort`
    renders *before* it exports CSV/TSV/JSON, so both internal columns reached every
    machine-readable output - which the code comment above the one copy this module does
    take says is exactly what a copy exists to prevent. Fixed in
    `cohort_categories.sample_categories`, at the source rather than at the export.

    Both frames are asserted, not just the Kestrel one: the reduction runs twice and a
    fix applied to one call site only would leave adVNTR's exports still carrying them.
    See
    `.superpowers/sdd/2026-08-06-issue-181-197-followups-plan/issue-cohort-internal-columns-leak-into-exports.md`.
    """
    kestrel = _kestrel_frame()
    advntr = _advntr_frame()

    cohort_summary.generate_cohort_summary_report(
        output_dir=str(tmp_path),
        kestrel_df=kestrel,
        advntr_df=advntr,
        summary_file="cohort_summary.html",
        config=load_config(None),
    )

    pd.testing.assert_frame_equal(kestrel, _kestrel_frame())
    pd.testing.assert_frame_equal(advntr, _advntr_frame())


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


def test_both_samples_reach_the_report_in_sorted_order(tmp_path) -> None:
    """**Specification**: the row order of a cohort report is the sorted sample order.

    The frames are built by iterating the discovered directories, so this is the end of
    the chain `test_cohort_inputs.py` starts: discovery sorts, and the report inherits
    it. It used to assert membership only, because the directories arrived in a `set`
    and the order was whatever the process's hash seed produced. Asserting the
    *positions* is what makes two runs of `vntyper cohort` comparable byte for byte.

    **Scope: directory inputs.** The two samples here are directories under `tmp_path`,
    where the sort is total and reproducible. It does not extend to zip inputs, whose
    sort key carries the random `cohort_zip_*` extraction directory - see
    `test_cohort_inputs.py::test_two_zip_inputs_are_ordered_by_their_random_temporary_directories_today`.
    """
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(_cohort_on_disk(tmp_path / "cohort"))],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=load_config(None),
    )

    html = (output_dir / "cohort_summary.html").read_text(encoding="utf-8")

    assert html.index(">sample_one</td>") < html.index(">sample_two</td>")
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


def test_an_export_written_after_the_report_carries_no_internal_columns(tmp_path) -> None:
    """**Specification**, at the boundary a user sees: the same fix, end to end.

    A cohort run with `--output-format csv` used to produce a `cohort_kestrel.csv`
    whose header carried `__row_result` and `__unified`, because the render mutated the
    frame the export was about to write. Both exports are checked, and the whole header
    is asserted rather than the two names alone - the columns and their order are the
    contract a downstream parser reads. `Sample` is last because the export order is the
    insertion order of the dicts read out of `pipeline_summary.json` with `Sample`
    appended, not `cohort_tables`' display order; that is unchanged and is not this
    fix's business. See the issue file named in the test above.
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

    kestrel_header = (output_dir / "cohort_kestrel.csv").read_text(encoding="utf-8").splitlines()[0]
    advntr_header = (output_dir / "cohort_advntr.csv").read_text(encoding="utf-8").splitlines()[0]

    assert kestrel_header == "Motif,Confidence,Flag,Sample"
    assert advntr_header == "VID,Flag,Sample"
