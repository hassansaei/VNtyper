"""
subthreshold.py

The below-reporting-floor signal (#266).

`filter_final_dataframe` ANDs six boolean gates. Five encode a *structural or artifact*
judgement -- "this row is not a candidate variant". The sixth, :data:`DEPTH_GATE`, encodes
a *magnitude* one: it is ``Confidence != "Negative"``, which reduces to
``Depth_Score >= reporting_floor``. A row that fails **only** that gate is a
well-formed pathogenic-frame candidate that was too faint to call, and until #266 it left
no trace at all: the sample rendered identically to one where nothing was ever found.

This module decides what that class is and what may be said about it. It says very little
on purpose. @hassansaei on #266: "add a single line noting that a below-threshold variant
with its depth_score was identified and filtered out, without printing its full statistics
as if it were a candidate."

Three properties are load-bearing:

* **Eligibility is strict.** A row killed by ``flag_filter_pass`` carries a flag that
  ``kestrel_config.json`` declares an artifact (#174) and its ``Depth_Score`` may be
  excellent. Calling that "subthreshold" would say *weak signal* where the truth is
  *strong signal, deliberately discarded*.
* **A depth-driven ALT verdict is still a depth verdict.** ``alt_filter_pass`` also fails a
  dupC row (``ALT == alt_filtering.gg_alt_value``) whose ``Depth_Score`` is below
  ``gg_depth_score_threshold``; shipped, that threshold *is* the reporting floor, so every
  dupC candidate below the floor fails two gates. Given ``alt_filtering``, such a row is
  eligible when its ALT failure is explainable by depth alone -- an excluded ALT never is.
  A dupC row at or above the floor that fails only a stricter custom GG threshold is not
  subthreshold: it is not below the reporting floor.
* **Every gate verdict must be explicit.** See :func:`_verdict`: "not True" and "False"
  are different claims, and only an explicit one counts in either direction.
* **Events, not rows.** The pre-result carries the same event against many motif contexts
  -- measured at 1,250 rows for one carrier, 2 of which survive the structural gates -- so
  a row count would misstate how much was found.

Nothing here filters, calls or ranks a variant. The reporting floor is unchanged, and a
subthreshold candidate is not a call.
"""

from __future__ import annotations

import json
import logging
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Final

import pandas as pd

logger = logging.getLogger(__name__)

#: The one gate whose failure -- alone -- defines a subthreshold candidate.
DEPTH_GATE: Final[str] = "depth_confidence_pass"

#: The gate a dupC row also fails on depth alone; see :func:`_alt_failed_on_depth_only`.
ALT_GATE: Final[str] = "alt_filter_pass"

#: The prefix every emitted note carries, so the report can find it among the ``##``
#: banner lines without matching on prose.
#:
#: This is a machine contract and therefore **code, not configuration**: the writer
#: (`kestrel_genotyping`) and the reader (`generate_report`) must agree on it, and they do
#: not load the same configuration file. Only the human sentence after it is configurable,
#: in ``kestrel_config.json``'s ``subthreshold_note.template``.
NOTE_MARKER: Final[str] = "Subthreshold candidate:"

#: The sidecar beside ``kestrel_result.tsv`` that keeps the signal in numbers, so
#: ``pipeline_summary.json`` need not parse the configurable sentence back.
SIGNAL_FILENAME: Final[str] = "subthreshold_signal.json"

#: The pre-result column holding a row's alternate-allele read support.
SUPPORT_COLUMN: Final[str] = "Estimated_Depth_AlternateVariant"

#: Columns the signal is computed from, beyond the gates themselves.
_REQUIRED_COLUMNS: Final[tuple[str, ...]] = ("Depth_Score", "POS", "REF", "ALT")

#: String spellings a boolean column can carry after a TSV round trip. Anything outside
#: these two sets is *unknown*, not False -- see :func:`_verdict`.
_TRUE_TOKENS: Final[frozenset[str]] = frozenset({"true", "1", "t", "yes"})
_FALSE_TOKENS: Final[frozenset[str]] = frozenset({"false", "0", "f", "no"})


@dataclass(frozen=True)
class SubthresholdSignal:
    """What one sample's eligible rows amount to.

    Attributes:
        events: Distinct ``(POS, REF, ALT)`` tuples among the eligible rows. This is what
            gets printed; it is the number of *events* seen, not of rows.
        rows: Eligible rows. Kept for the log line only -- printing it would overstate what
            was found, because one event appears once per motif context.
        best_depth_score: The highest ``Depth_Score`` among the eligible rows, i.e. how
            close the strongest suppressed candidate came to the floor.
        floor: The reporting floor it fell below.
        support_count: Alternate-allele read support of the row with the best
            ``Depth_Score``, or None when the pre-result does not carry it.
    """

    events: int
    rows: int
    best_depth_score: float
    floor: float
    support_count: int | None = None


def _verdict(value: object) -> bool | None:
    """The explicit boolean one gate cell carries, or None when it carries none.

    A gate column is ``bool`` in memory and after a round trip through
    ``kestrel_pre_result.tsv`` whenever it holds nothing else -- but a single missing value
    makes pandas read the whole column back as ``object``, where the cells are the strings
    ``"True"`` and ``"False"``. ``bool("False")`` is ``True``, so a cast would turn a
    failing gate into a passing one and widen eligibility in silence.

    The return is deliberately three-valued rather than two. "Not True" and "False" are
    different claims: a row whose depth verdict was never recorded is not a row known to be
    subthreshold, and reading its absence as a failure would let a malformed evidence file
    manufacture exactly the misleading note this module exists to prevent.

    Args:
        value: One cell of a gate column. Python bool, numpy scalar, string, number,
            ``None``, ``pd.NA`` and ``NaN`` are all reachable.

    Returns:
        bool | None: The verdict, or None when the cell carries no recognisable one.
    """
    if isinstance(value, str):
        token = value.strip().lower()
        if token in _TRUE_TOKENS:
            return True
        if token in _FALSE_TOKENS:
            return False
        return None
    # numpy scalars do not subclass their Python counterparts -- `isinstance(np.True_,
    # bool)` is False -- so unwrap before every type test below.
    if hasattr(value, "item") and not isinstance(value, bytes):
        try:
            value = value.item()
        except (AttributeError, ValueError):
            return None
    if isinstance(value, bool):
        return value
    if value is None or value is pd.NA:
        return None
    if isinstance(value, (int, float)):
        # `pd.isna` cannot raise for a plain int or float, and by here the value is one:
        # every numpy scalar was unwrapped above. An `except` around this would be dead.
        if pd.isna(value):
            return None
        if value == 1:
            return True
        if value == 0:
            return False
    return None


def _explicitly(frame: pd.DataFrame, column: str, wanted: bool) -> pd.Series:
    """Rows whose gate cell explicitly carries ``wanted``.

    Args:
        frame: The pre-filter frame.
        column: The gate column's name.
        wanted: The verdict a row must carry to be kept.

    Returns:
        pd.Series: A boolean mask. A cell with no recognisable verdict is False whatever
        ``wanted`` is, so an unknown disqualifies the row in both directions.
    """
    verdicts = [_verdict(value) for value in frame[column]]
    unknown = sum(1 for verdict in verdicts if verdict is None)
    if unknown:
        logger.warning(
            "%d of %d rows carry no recognisable %r verdict; they cannot be subthreshold.",
            unknown,
            len(verdicts),
            column,
        )
    return pd.Series([verdict is wanted for verdict in verdicts], index=frame.index, dtype=bool)


def _alt_failed_on_depth_only(frame: pd.DataFrame, alt_filtering: Mapping[str, Any]) -> pd.Series:
    """Rows whose ALT gate failure is the GG depth rule and nothing else.

    Recomputed from the row's own ALT and Depth_Score with the same rule
    `variant_parsing.filter_by_alt_values_and_finalize` applies, so an excluded ALT (which
    fails whatever its depth) is never mistaken for a weak dupC.

    Args:
        frame: The pre-filter frame.
        alt_filtering: ``kestrel_config.json``'s ``alt_filtering`` section.

    Returns:
        pd.Series: A boolean mask.
    """
    gg = alt_filtering.get("gg_alt_value", "GG")
    excluded = list(alt_filtering.get("exclude_alts", []) or [])
    try:
        threshold = float(alt_filtering["gg_depth_score_threshold"])
    except (KeyError, TypeError, ValueError):
        return pd.Series(False, index=frame.index)
    scores = pd.to_numeric(frame["Depth_Score"], errors="coerce")
    alts = frame["ALT"].astype(str)
    return (alts == gg) & ~alts.isin(excluded) & (scores < threshold)


def detect(
    frame: pd.DataFrame,
    gate_columns: Sequence[str],
    floor: float,
    alt_filtering: Mapping[str, Any] | None = None,
) -> SubthresholdSignal | None:
    """Reduce a pre-filter frame to its subthreshold signal, or to nothing.

    Args:
        frame: The frame `filter_final_dataframe` writes to ``kestrel_pre_result.tsv`` --
            every candidate row, every gate column, nothing dropped.
        gate_columns: The gate columns to require, including :data:`DEPTH_GATE`. Passed in
            rather than restated so that a seventh gate added to
            `kestrel_genotyping.FILTER_COLUMNS` cannot silently widen eligibility here.
        floor: The reporting floor, for the record. It decides nothing -- the decision is
            the gate's, already made upstream by `calculate_depth_score_and_assign_confidence`.
        alt_filtering: The Kestrel ``alt_filtering`` config. When given, a row failing
            :data:`ALT_GATE` stays eligible if that failure is the GG depth rule alone.
            When None, :data:`ALT_GATE` is required to pass like any structural gate.

    Returns:
        SubthresholdSignal | None: The signal, or None when no row is eligible or the frame
        cannot support the question.
    """
    if frame is None or frame.empty:
        return None
    if DEPTH_GATE not in gate_columns:
        logger.warning(
            "Subthreshold detection skipped: %r is not among the gate columns %s, so there "
            "is no magnitude judgement to separate from the structural ones.",
            DEPTH_GATE,
            list(gate_columns),
        )
        return None

    missing = sorted({*gate_columns, *_REQUIRED_COLUMNS} - set(frame.columns))
    if missing:
        logger.warning("Subthreshold detection skipped: the pre-result frame is missing %s.", missing)
        return None

    mask = pd.Series(True, index=frame.index)
    for column in gate_columns:
        if column == DEPTH_GATE:
            continue
        passed = _explicitly(frame, column, True)
        if column == ALT_GATE and alt_filtering is not None:
            passed |= _explicitly(frame, column, False) & _alt_failed_on_depth_only(frame, alt_filtering)
        mask &= passed
    # Explicitly False, not "not True": see :func:`_verdict`.
    mask &= _explicitly(frame, DEPTH_GATE, False)

    eligible = frame[mask]
    if eligible.empty:
        return None

    scores = pd.to_numeric(eligible["Depth_Score"], errors="coerce")
    usable = eligible[scores.notna()]
    if usable.empty:
        logger.warning(
            "%d row(s) fail only the depth gate but carry no numeric Depth_Score; no note.",
            len(eligible),
        )
        return None

    events = usable[["POS", "REF", "ALT"]].astype(str).drop_duplicates()
    support = float("nan")
    if SUPPORT_COLUMN in usable.columns:
        support = float(pd.to_numeric(usable[SUPPORT_COLUMN], errors="coerce")[scores.idxmax()])
    signal = SubthresholdSignal(
        events=int(len(events)),
        rows=int(len(usable)),
        best_depth_score=float(scores.max()),
        floor=float(floor),
        support_count=None if math.isnan(support) else int(support),
    )
    logger.info(
        "Subthreshold signal: %d event(s) over %d row(s), best Depth_Score %g, floor %g.",
        signal.events,
        signal.rows,
        signal.best_depth_score,
        signal.floor,
    )
    return signal


def detect_from_file(
    path: str | Path,
    gate_columns: Sequence[str],
    floor: float,
    alt_filtering: Mapping[str, Any] | None = None,
) -> SubthresholdSignal | None:
    """:func:`detect`, reading the frame from a written ``kestrel_pre_result.tsv``.

    The pre-result is the documented evidence artefact and is written by
    `filter_final_dataframe` before it can raise, so reading it back means the printed line
    describes exactly what is on disk. Any read failure yields no note: the note is an
    annotation, and losing it must never cost a result.

    Args:
        path: Path to ``kestrel_pre_result.tsv``.
        gate_columns: As :func:`detect`.
        floor: As :func:`detect`.
        alt_filtering: As :func:`detect`.

    Returns:
        SubthresholdSignal | None: The signal, or None.
    """
    location = Path(path)
    if not location.is_file():
        logger.warning("No pre-result at %s; no subthreshold note.", location)
        return None
    try:
        frame = pd.read_csv(location, sep="\t")
    except (OSError, ValueError) as error:
        # Narrow on purpose rather than a bare `except Exception`: this is a local read,
        # not a process boundary, so it does not belong in the reviewed BLE001 inventory
        # (AGENTS.md). The two cover what this call can actually raise -- `OSError` for the
        # filesystem, and `ValueError` for both `pandas.errors.ParserError` (a ragged file)
        # and `UnicodeDecodeError` (a non-UTF-8 one), each of which subclasses it.
        logger.warning("Could not read %s for the subthreshold note: %s", location, error)
        return None
    return detect(frame, gate_columns, floor, alt_filtering)


def format_note(signal: SubthresholdSignal, template: str) -> str | None:
    """Render one signal as the single line that goes into the result banner.

    Args:
        signal: The signal to describe.
        template: The configured sentence, using ``{marker}``, ``{events}``, ``{noun}``,
            ``{best_depth_score}``, ``{floor}`` and ``{rows}``.

    Returns:
        str | None: One line, with every newline and tab collapsed to a space so it cannot
        split the TSV it is written into. None when the template names a field that does
        not exist -- a misconfigured sentence must not abort a run.
    """
    fields: dict[str, object] = {
        "marker": NOTE_MARKER,
        "events": signal.events,
        "rows": signal.rows,
        "noun": "candidate variant" if signal.events == 1 else "candidate variants",
        "best_depth_score": f"{signal.best_depth_score:.5g}",
        "floor": f"{signal.floor:.5g}",
    }
    try:
        rendered = template.format(**fields)
    except (AttributeError, IndexError, KeyError, ValueError) as error:
        # The four `str.format` raises for an operator-supplied template, measured rather
        # than guessed and pinned one per class by `TestRendering`:
        #   AttributeError  {marker.upper.foo}   an attribute lookup that does not resolve
        #   IndexError      {marker[9]}, {0}     an index, or a positional field
        #   KeyError        {nope}               a field name that does not exist
        #   ValueError      {marker, {events:s}  unterminated, or a bad format spec
        # `AttributeError` was missing until Sourcery caught it on #276, and it was the one
        # that mattered: this function's contract is that a misconfigured sentence costs the
        # note and never the run, and an escaping AttributeError would have aborted report
        # generation instead.
        logger.error("subthreshold_note.template is not renderable (%s); no note will be written.", error)
        return None
    return " ".join(rendered.split())


def write_signal(signal: SubthresholdSignal, output_dir: str | Path) -> None:
    """Record the signal as :data:`SIGNAL_FILENAME` beside the note it belongs to.

    Args:
        signal: The signal the note was rendered from.
        output_dir: The Kestrel output directory.
    """
    try:
        (Path(output_dir) / SIGNAL_FILENAME).write_text(json.dumps(asdict(signal), indent=4) + "\n", encoding="utf-8")
    except OSError as error:
        # An annotation: losing it must never cost a result.
        logger.warning("Could not write %s: %s", SIGNAL_FILENAME, error)


def summary_record(state: str, kestrel_dir: str | Path) -> dict[str, object]:
    """The ``subthreshold_signal`` mapping ``pipeline_summary.json`` carries.

    Args:
        state: The Kestrel screening result when it is the subthreshold state, else ``"none"``.
        kestrel_dir: The Kestrel output directory holding :data:`SIGNAL_FILENAME`.

    Returns:
        dict: ``state``, ``best_depth_score``, ``event_count`` and ``support_count``. The
        numbers are read only for a subthreshold state, so a sidecar an earlier run left
        behind can never describe a sample that has no note.
    """
    record: dict[str, object] = {"state": state, "best_depth_score": None, "event_count": 0, "support_count": None}
    if state == "none":
        return record
    try:
        saved = json.loads((Path(kestrel_dir) / SIGNAL_FILENAME).read_text(encoding="utf-8"))
        record.update(
            best_depth_score=float(saved["best_depth_score"]),
            event_count=int(saved["events"]),
            support_count=None if saved.get("support_count") is None else int(saved["support_count"]),
        )
    except (OSError, ValueError, KeyError, TypeError) as error:
        logger.warning("No readable %s for the subthreshold record: %s", SIGNAL_FILENAME, error)
    return record


def find_note(comments: Iterable[str] | None) -> str | None:
    """Return the marked line among a result file's comment lines.

    Args:
        comments: The ``comments`` list `summary.parse_tsv` records, whose entries have
            already had their leading ``#`` stripped -- though this tolerates them either
            way, because a caller reading the file directly has not.

    Returns:
        str | None: The first marked line, stripped, or None.
    """
    for comment in comments or ():
        text = str(comment).lstrip("#").strip()
        if text.startswith(NOTE_MARKER):
            return text
    return None
