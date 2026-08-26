"""
subthreshold.py

The below-reporting-floor signal (#266).

`filter_final_dataframe` ANDs six boolean gates. Five encode a *structural or artifact*
judgement -- "this row is not a candidate variant". The sixth, :data:`DEPTH_GATE`, encodes
a *magnitude* one: it is ``Confidence != "Negative"``, which reduces to
``Depth_Score >= depth_score_thresholds.low``. A row that fails **only** that gate is a
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
* **Every gate verdict must be explicit.** See :func:`_verdict`: "not True" and "False"
  are different claims, and only an explicit one counts in either direction.
* **Events, not rows.** The pre-result carries the same event against many motif contexts
  -- measured at 1,250 rows for one carrier, 2 of which survive the structural gates -- so
  a row count would misstate how much was found.

Nothing here filters, calls or ranks a variant. The reporting floor is unchanged, and a
subthreshold candidate is not a call.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import pandas as pd

logger = logging.getLogger(__name__)

#: The one gate whose failure -- alone -- defines a subthreshold candidate.
DEPTH_GATE: Final[str] = "depth_confidence_pass"

#: The prefix every emitted note carries, so the report can find it among the ``##``
#: banner lines without matching on prose.
#:
#: This is a machine contract and therefore **code, not configuration**: the writer
#: (`kestrel_genotyping`) and the reader (`generate_report`) must agree on it, and they do
#: not load the same configuration file. Only the human sentence after it is configurable,
#: in ``kestrel_config.json``'s ``subthreshold_note.template``.
NOTE_MARKER: Final[str] = "Subthreshold candidate:"

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
    """

    events: int
    rows: int
    best_depth_score: float
    floor: float


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
        try:
            if pd.isna(value):
                return None
        except (TypeError, ValueError):
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


def detect(
    frame: pd.DataFrame,
    gate_columns: Sequence[str],
    floor: float,
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
        mask &= _explicitly(frame, column, True)
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
    signal = SubthresholdSignal(
        events=int(len(events)),
        rows=int(len(usable)),
        best_depth_score=float(scores.max()),
        floor=float(floor),
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

    Returns:
        SubthresholdSignal | None: The signal, or None.
    """
    location = Path(path)
    if not location.is_file():
        logger.warning("No pre-result at %s; no subthreshold note.", location)
        return None
    try:
        frame = pd.read_csv(location, sep="\t")
    except Exception as error:  # noqa: BLE001 - an unreadable annotation source must not abort a run
        logger.warning("Could not read %s for the subthreshold note: %s", location, error)
        return None
    return detect(frame, gate_columns, floor)


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
    except (KeyError, IndexError, ValueError) as error:
        logger.error("subthreshold_note.template is not renderable (%s); no note will be written.", error)
        return None
    return " ".join(rendered.split())


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
