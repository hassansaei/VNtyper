"""
screening_summary.py

Module Purpose:
---------------
Turns the Kestrel and adVNTR result frames plus one quality gate into the
screening state, and looks that state up in ``report_config.json`` to get the
sentence the report prints.

Extracted from ``generate_report.py`` (861 LOC, 4% covered) under AGENTS.md rule
3. It is pure: two DataFrames, two numbers and a config mapping in, a
:class:`ScreeningSummary` out. No filesystem beyond reading the packaged config,
no rendering.

**The interpretive text is config-driven and must stay that way** (AGENTS.md). Nothing
here composes, rewords or conditions a sentence; it computes a state and reads the
message for that state out of the configuration. A state with no configured
message is marked by ``ScreeningSummary.matched_rule`` and presented with
``UNAVAILABLE_SUMMARY_MESSAGE``: otherwise the fallback default could announce a
negative screening for a sample two algorithms called positive.

The state has three axes:

* ``kestrel_result`` -- one of the ``result`` values under
  ``algorithm_logic.kestrel``, or that block's ``default``.
* ``advntr_result`` -- likewise for ``algorithm_logic.advntr``, plus
  :data:`NOT_PERFORMED` when the adVNTR stage did not run.
* ``quality_metrics_pass`` -- the coverage QC verdict: **both** mean VNTR
  coverage and the uncovered fraction against their configured thresholds
  (#172). The caller evaluates it; this module only reads ``.passed``, so the
  axis and the ``coverage_qc`` column written into the coverage summary are the
  same value and cannot disagree.

``is_positive`` is derived from the first two by comparing against each block's
declared ``default``, not by looking for a word in the rendered sentence. The
template used to test ``'negative' not in summary_text`` to decide whether to
style the box as a positive finding, and **not one** of the configured messages
contains that word -- only the fallback default does -- so every configured
message, including "No variant detected by either genotyping method", rendered as
a positive finding.

A fourth thing the state model has to carry is *whether each stage ran at all*, which
is not one of the three axes and cannot be: ``compute_algorithm_result`` reduces an empty
frame to the block's ``default``, and a stage that ran and lost its result hands it exactly
that empty frame. So the two execution axes are separate, they come from
``summary_steps.get_step_state``, and they decide whether the message the rule table
selected may be presented at all -- see :attr:`ScreeningSummary.state_is_established`.

Functions:
    load_report_config: Read the packaged ``report_config.json``
    compute_algorithm_result: One results frame + one rule block to a state value
    execution_state: One pipeline-summary step state to one execution state
    algorithm_state_text: One algorithm's state as the provenance line prints it
    message_segments: One configured message as the parts the report renders
    render_segments: Those parts reassembled, which must equal the message
    supports_subthreshold: Whether this configuration can describe a subthreshold negative
    state_chips: The computed state as the masthead's chip row
    build_screening_summary: The three state axes to the configured message
"""

from __future__ import annotations

import json
import logging
import os
from collections.abc import Sequence
from dataclasses import dataclass, replace
from typing import Any, Final, Literal

import pandas as pd

from vntyper.scripts.coverage_qc import CoverageQC
from vntyper.scripts.summary_steps import STEP_ABSENT, STEP_READ, STEP_UNREADABLE

logger = logging.getLogger(__name__)

#: The ``advntr_result`` value recorded when the adVNTR stage did not run at all.
#: Distinct from that block's ``default`` ("negative"), which means adVNTR ran and
#: found nothing.
NOT_PERFORMED = "none"

#: The Kestrel state for a sample that called nothing but held a candidate failing only the
#: depth gate (#266). It is a *described negative*, never a call: it is declared in
#: ``report_config.json``'s ``algorithm_logic.kestrel.non_finding_results``, so
#: :func:`is_finding` returns False for it and nothing styles it as a finding.
#:
#: The vocabulary is M6's (`.planning`'s milestone-7 reshape design), and only the
#: vocabulary: ``state_reason``, the top-level ``algorithms`` block and execution-state
#: pre-emption are that design's other half and remain deferred to #173/#223/#224.
SUBTHRESHOLD_RESULT: Final[str] = "negative_subthreshold"

#: How a stage's *execution* is reported, as distinct from what it called. The three
#: values are :func:`~vntyper.scripts.summary_steps.get_step_state`'s three states in the
#: vocabulary the report presents, and they exist because "recorded" is not "produced a
#: readable result": ``record_step`` writes the step either way and flags the missing file
#: (#212), leaving exactly the shape a stage that genotyped and called nothing produces.
EXECUTION_PERFORMED = "performed"
EXECUTION_NOT_PERFORMED = "not-performed"
EXECUTION_FAILED = "failed"

#: The one translation from step state to execution state. An unrecognised step state maps
#: to :data:`EXECUTION_FAILED`, which is the conservative direction: a state this module
#: does not know has established nothing about the sample, and presenting it as a result
#: is the defect the whole distinction exists to remove.
_EXECUTION_BY_STEP_STATE = {
    STEP_READ: EXECUTION_PERFORMED,
    STEP_ABSENT: EXECUTION_NOT_PERFORMED,
    STEP_UNREADABLE: EXECUTION_FAILED,
}

#: Used when ``report_config.json`` declares no default for an algorithm block.
FALLBACK_ALGORITHM_RESULT = "none"

#: Used when ``report_config.json`` declares no ``screening_summary_default``.
FALLBACK_SUMMARY_MESSAGE = "The screening was negative (no valid Kestrel or adVNTR data)."

#: What the report says when the summary could not be computed at all, or when it was
#: computed from a state one of the stages never established.
UNAVAILABLE_SUMMARY_MESSAGE = "No summary available."

#: The separator every configured message uses between its parts. It is markup, and the
#: report no longer interpolates it as markup: each part is rendered as an element of its
#: own, autoescaped like everything else. ``segments`` in ``report_config.json`` is that
#: split, stored beside the ``message`` it came from; the message stays verbatim and
#: authoritative, and :func:`render_segments` is what proves the two still agree.
SEGMENT_SEPARATOR = "<br>"


def message_segments(rule: dict[str, Any]) -> tuple[str, ...]:
    """The ordered parts one configured rule's message is rendered as.

    A rule that declares ``segments`` is taken at its word; one that declares only
    ``message`` is split on :data:`SEGMENT_SEPARATOR`. The second path is not a migration
    artefact to be removed later - ``report_config.json`` is configuration, a deployment
    may ship its own, and a config carrying only ``message`` keys has to keep rendering.

    Args:
        rule: One entry of ``screening_summary_rules``.

    Returns:
        tuple[str, ...]: The parts, in the order they are rendered.
    """
    segments = rule.get("segments")
    if isinstance(segments, list) and segments:
        return tuple(str(segment) for segment in segments)
    return tuple(str(rule.get("message", "")).split(SEGMENT_SEPARATOR))


def render_segments(rule: dict[str, Any]) -> str:
    """Reassemble one rule's message from the parts the report actually renders.

    The round trip is the contract: ``render_segments(rule) == rule["message"]`` for every
    configured rule. Asserting it on *rendering* rather than on any rule having a given
    number of segments is what keeps the migration honest - the configured wording is
    clinical text nothing here may reword, only split.

    Args:
        rule: One entry of ``screening_summary_rules``.

    Returns:
        str: The message as the reader receives it, separators restored.
    """
    return SEGMENT_SEPARATOR.join(message_segments(rule))


def execution_state(step_state: str) -> str:
    """Translate one pipeline-summary step state into an execution state.

    Args:
        step_state: One of :mod:`vntyper.scripts.summary_steps`' ``STEP_*`` state values.

    Returns:
        str: :data:`EXECUTION_PERFORMED`, :data:`EXECUTION_NOT_PERFORMED` or
        :data:`EXECUTION_FAILED`. An unrecognised state reads as failed.
    """
    return _EXECUTION_BY_STEP_STATE.get(step_state, EXECUTION_FAILED)


def algorithm_state_text(execution: str, result: str) -> str:
    """One algorithm's state as the report's raw-state provenance line prints it.

    The execution state comes first because the computed ``result`` is only a result when
    the stage produced one: an absent or unreadable stage hands
    :func:`compute_algorithm_result` an empty frame and gets back the block's ``default``,
    so a line that printed the result would read "Kestrel: negative" for a stage that
    never called anything.

    Args:
        execution: One of the ``EXECUTION_*`` values.
        result: The computed ``kestrel_result``/``advntr_result``.

    Returns:
        str: The computed result token when the stage produced one, else the execution
        state in words.
    """
    if execution == EXECUTION_PERFORMED:
        return result
    return execution.replace("-", " ")


@dataclass(frozen=True)
class ScreeningSummary:
    """The screening state and the sentence configured for it.

    Attributes:
        text: The message the report may present. It is the matched rule's message, or
            :data:`UNAVAILABLE_SUMMARY_MESSAGE` when :attr:`state_is_established` or
            ``matched_rule`` is False. A message chosen for an unestablished state, and
            a fallback chosen because no rule exists, are not evidence about the sample.
        is_positive: Whether either algorithm reported a finding. Drives the
            template's emphasis; computed from the state, never from ``text``.
        kestrel_result: The computed Kestrel state value.
        advntr_result: The computed adVNTR state value.
        quality_metrics_pass: Whether the sample met every configured coverage
            threshold - mean depth and uncovered fraction alike (#172).
        matched_rule: Whether a configured rule matched the computed state. False means
            the state has no message of its own. It says nothing about whether the state
            was established: that is what the two execution axes are for.
        kestrel_execution: Whether the Kestrel stage produced a readable result, was
            never asked to run, or ran and lost it.
        advntr_execution: The same for the adVNTR stage.
    """

    text: str
    is_positive: bool
    kestrel_result: str
    advntr_result: str
    quality_metrics_pass: bool
    matched_rule: bool
    kestrel_execution: str = EXECUTION_PERFORMED
    advntr_execution: str = EXECUTION_PERFORMED
    segments: tuple[str, ...] = ()

    @property
    def rendered_segments(self) -> tuple[str, ...]:
        """The parts the report renders, one element each.

        Falls back to splitting :attr:`text` so a summary built without ``segments`` -
        by a test, or by a caller that has only the message - still renders as elements
        rather than putting a literal ``<br>`` in front of the reader.

        Returns:
            tuple[str, ...]: One or more parts, never empty.
        """
        return self.segments or tuple(self.text.split(SEGMENT_SEPARATOR))

    @property
    def state_is_established(self) -> bool:
        """Whether every axis the rule table matched on describes something the run did.

        The two algorithms are treated differently, and the asymmetry is read off the
        configuration rather than chosen. ``advntr_result`` has a value for a stage that
        never ran -- :data:`NOT_PERFORMED` -- and ten of the forty configured rules are
        keyed on it, so an absent adVNTR stage is a state the table covers *in words*.
        ``kestrel_result`` has no such value: a Kestrel stage that produced nothing hands
        :func:`compute_algorithm_result` an empty frame, which returns the block's
        ``default`` -- the same ``negative`` token a stage that genotyped and called
        nothing produces. So a Kestrel stage that did not produce a result leaves the
        state unestablished whichever way it failed to, while an adVNTR stage does so only
        when it ran and lost its result.

        Returns:
            bool: True when the configured message describes states the run reached.
        """
        return self.kestrel_execution == EXECUTION_PERFORMED and self.advntr_execution != EXECUTION_FAILED

    @property
    def emphasis(self) -> Literal["finding", "no-finding", "indeterminate"]:
        """How the report should style this state.

        :attr:`state_is_established` is checked first, then ``matched_rule``, and both win
        over ``is_positive``. A state one stage never established is unknown, and so is a
        state with no configured rule -- neither is a negative, and neither is a finding
        the report may assert.

        ``quality_metrics_pass`` never appears in this expression -- call status, QC
        status and rule-match status are orthogonal. The rule table already describes
        a positive call with failing coverage QC as a finding, and ``is_positive`` is
        derived from the algorithm calls independently of QC by design; letting failed
        QC downgrade a matched positive rule to "indeterminate" would silently
        reclassify a confirmed pathogenic call with poor coverage as "state unknown",
        contradicting both the rule table and ``is_positive``.

        In practice ``matched_rule`` is False only through the ``except Exception``
        path in :func:`build_screening_summary`: every one of the 40 reachable
        ``(kestrel_result, advntr_result, quality_metrics_pass)`` combinations
        resolves to a configured rule (``tests/unit/test_screening_summary.py::
        test_every_reachable_state_has_its_own_message``). So a ``"indeterminate"``
        reached that way is genuinely exceptional and cannot mislabel an ordinary
        all-negative report.

        Returns:
            Literal["finding", "no-finding", "indeterminate"]: The emphasis to apply.
        """
        if not self.state_is_established:
            return "indeterminate"
        if self.matched_rule is False:
            return "indeterminate"
        if self.is_positive:
            return "finding"
        return "no-finding"


#: The four tones a state chip can be drawn in. Each names a token in the shared token
#: layer; the chip's *value* is always a word, so the tone is never the only carrier.
TONE_FINDING = "finding"
TONE_CAUTION = "caution"
TONE_OK = "ok"
TONE_NONE = "none"

#: What a chip says for a stage that produced no result. Neither is a statement about the
#: sample: one says the stage was never asked to run, the other that it ran and nothing
#: could be read from it. Neither may ever be replaced by a result word.
NOT_PERFORMED_CHIP = "Not performed"
NOT_AVAILABLE_CHIP = "Not available"

#: What the concordance chip says when there is nothing to compare. A stage that produced
#: no result cannot agree or disagree with one that did, and a chip claiming either would
#: be the same defect as reporting an absence as a negative.
NOT_ASSESSABLE_CHIP = "Not assessable"

#: The two words the cross-match state has. They compress the sentence
#: ``build_cross_match_summary`` writes; the sentence itself is still printed in full.
MATCH_CHIP = "Match"
NO_MATCH_CHIP = "No match"

#: The chip labels, named once so the chip row and its tests cannot drift apart.
KESTREL_LABEL = "Kestrel"
ADVNTR_LABEL = "adVNTR"
CONCORDANCE_LABEL = "Concordance"


@dataclass(frozen=True)
class StateChip:
    """One computed fact, as the masthead's chip row shows it.

    Attributes:
        label: What the fact is about.
        value: The fact, always in words.
        tone: Which state token colours it. Emphasis only - every chip is readable in
            greyscale, in print and with any colour vision.
    """

    label: str
    value: str
    tone: str


def result_word(result: str) -> str:
    """One computed result token as a chip prints it.

    A pure transformation of the token rather than a lookup table, and deliberately: a
    table would have to be extended by hand whenever ``report_config.json`` declares a new
    ``result``, and until someone did the chip would either fall back to a wrong word or
    crash. ``High_Precision`` reads "High precision"; ``positive flagged`` reads "Positive
    flagged"; nothing is composed, translated or interpreted.

    Args:
        result: A computed ``kestrel_result``/``advntr_result`` value.

    Returns:
        str: The same word, spaced and sentence-cased.
    """
    return result.replace("_", " ").capitalize()


def algorithm_chip(
    label: str,
    execution: str,
    result: str,
    default: str,
    non_finding: Sequence[str] = (),
) -> StateChip:
    """One algorithm's chip: what it called, or that it called nothing.

    This is a **second, independent** :func:`is_finding` call site --
    :func:`build_screening_summary` makes the first, for ``is_positive``. They must be
    given the same configuration or the masthead chip and the summary sentence can
    disagree about the same sample, which is how #266's non-calling state would have been
    toned as a finding while the sentence said otherwise.

    Args:
        label: :data:`KESTREL_LABEL` or :data:`ADVNTR_LABEL`.
        execution: One of the ``EXECUTION_*`` values.
        result: The computed result token.
        default: The algorithm block's configured ``default``.
        non_finding: The block's configured ``non_finding_results``, if any.

    Returns:
        StateChip: The chip to render.
    """
    if execution == EXECUTION_NOT_PERFORMED:
        return StateChip(label=label, value=NOT_PERFORMED_CHIP, tone=TONE_NONE)
    if execution != EXECUTION_PERFORMED:
        return StateChip(label=label, value=NOT_AVAILABLE_CHIP, tone=TONE_CAUTION)
    tone = TONE_FINDING if is_finding(result, default, non_finding) else TONE_NONE
    return StateChip(label=label, value=result_word(result), tone=tone)


def concordance_chip(
    summary: ScreeningSummary, cross_match_available: bool, cross_match_is_positive: bool
) -> StateChip:
    """Whether the two algorithms agreed, or whether that is a question at all.

    Concordance is the cross-match stage's own comparison, not a second opinion computed
    here. It needs two results to compare, so a stage that produced none makes the chip
    unassessable rather than negative - "the two did not agree" and "there was nothing to
    agree with" are different facts, and only one of them is about the sample.

    Args:
        summary: The computed screening state.
        cross_match_available: Whether the cross-match stage produced a comparison.
        cross_match_is_positive: Whether that comparison found a match.

    Returns:
        StateChip: The chip to render.
    """
    both_performed = summary.kestrel_execution == EXECUTION_PERFORMED and (
        summary.advntr_execution == EXECUTION_PERFORMED
    )
    if not both_performed or not cross_match_available:
        return StateChip(label=CONCORDANCE_LABEL, value=NOT_ASSESSABLE_CHIP, tone=TONE_NONE)
    if cross_match_is_positive:
        return StateChip(label=CONCORDANCE_LABEL, value=MATCH_CHIP, tone=TONE_FINDING)
    return StateChip(label=CONCORDANCE_LABEL, value=NO_MATCH_CHIP, tone=TONE_NONE)


def state_chips(
    summary: ScreeningSummary,
    report_config: dict[str, Any],
    *,
    cross_match_available: bool,
    cross_match_is_positive: bool,
) -> list[StateChip]:
    """The computed screening state as the masthead's chips.

    Args:
        summary: The computed screening state.
        report_config: The parsed ``report_config.json``, for each block's ``default``
            and ``non_finding_results``.
        cross_match_available: Whether the cross-match stage produced a comparison.
        cross_match_is_positive: Whether that comparison found a match.

    Returns:
        list[StateChip]: Kestrel, adVNTR and concordance, in scan order.
    """
    algorithm_logic = report_config.get("algorithm_logic", {})
    kestrel_logic = algorithm_logic.get("kestrel", {})
    advntr_logic = algorithm_logic.get("advntr", {})
    return [
        algorithm_chip(
            KESTREL_LABEL,
            summary.kestrel_execution,
            summary.kestrel_result,
            kestrel_logic.get("default", FALLBACK_ALGORITHM_RESULT),
            kestrel_logic.get("non_finding_results", ()),
        ),
        algorithm_chip(
            ADVNTR_LABEL,
            summary.advntr_execution,
            summary.advntr_result,
            advntr_logic.get("default", FALLBACK_ALGORITHM_RESULT),
            advntr_logic.get("non_finding_results", ()),
        ),
        concordance_chip(summary, cross_match_available, cross_match_is_positive),
    ]


def load_report_config() -> dict[str, Any]:
    """Load ``report_config.json`` from beside this module.

    Returns:
        dict[str, Any]: The report configuration, or ``{}`` if it cannot be read.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(script_dir, "report_config.json")
    try:
        with open(config_path) as f:
            report_config = json.load(f)
        logger.info("Loaded report config from %s", config_path)
        return report_config
    except Exception as e:
        logger.error("Failed to load report config: %s", e)
        return {}


def _condition_holds(actual: str, expected: Any) -> bool:
    """Evaluate one configured condition against one cell value.

    Two spellings are supported. A mapping is an explicit
    ``{"operator": ..., "value": ...}`` with operators ``==``, ``!=``, ``in`` and
    ``not in``; anything else is an implicit equality, or membership when it is a
    list.

    Args:
        actual: The cell value, already stringified and stripped.
        expected: The configured expectation.

    Returns:
        bool: Whether the condition holds. An unsupported operator is False.
    """
    if isinstance(expected, dict):
        op = expected.get("operator")
        exp_val = expected.get("value")
        if op == "==":
            return actual == str(exp_val).strip()
        if op == "!=":
            return actual != str(exp_val).strip()
        if op in ("in", "not in"):
            options = exp_val if isinstance(exp_val, list) else [exp_val]
            return (actual in options) if op == "in" else (actual not in options)
        logger.debug("Unsupported operator %r; condition fails.", op)
        return False
    if isinstance(expected, list):
        return actual in expected
    return actual == str(expected)


def compute_algorithm_result(df: pd.DataFrame, logic_config: dict[str, Any]) -> str:
    """Reduce one algorithm's results frame to a single state value.

    Only the first row is examined: the frame is sorted so the most confident
    call leads, and the screening summary describes the sample, not each variant.

    Args:
        df: The results frame. An empty frame yields the configured default.
        logic_config: One block of ``report_config.json``'s ``algorithm_logic``.

    Returns:
        str: The ``result`` of the first matching rule, else the block's
        ``default``.
    """
    default = logic_config.get("default", FALLBACK_ALGORITHM_RESULT)
    if df.empty:
        logger.debug("DataFrame is empty; returning default result %r.", default)
        return default

    row = df.iloc[0]
    logger.debug("Data row for evaluation: %s", row.to_dict())

    for idx, rule in enumerate(logic_config.get("rules", [])):
        conditions = rule.get("conditions", {})
        for col, expected in conditions.items():
            if col not in row:
                logger.debug("Rule %s: column %r not found; rule fails.", idx, col)
                break
            actual = str(row.get(col, "")).strip()
            if not _condition_holds(actual, expected):
                logger.debug("Rule %s: condition on %r not met (actual=%r).", idx, col, actual)
                break
        else:
            result = rule.get("result")
            logger.debug("Rule %s PASSED; returning result: %s", idx, result)
            return result

    logger.debug("No rule matched; returning default result %r.", default)
    return default


def supports_subthreshold(report_config: dict[str, Any]) -> bool:
    """Whether this configuration can describe a below-reporting-floor negative (#266).

    True only when the Kestrel block declares :data:`SUBTHRESHOLD_RESULT` under
    ``non_finding_results``. A ``report_config.json`` written before #266 does not, and
    two things then have to be suppressed **together**:

    * the promotion, because :func:`is_finding` would classify the promoted token as a
      finding under that configuration and the report would style a suppressed candidate
      as a call -- the inversion #266 exists to prevent;
    * the note's rendering, because the eight configured sentences that explain it are
      absent too, and the legacy ones say "No variant detected by either genotyping
      method" -- a report printing that beside a line announcing a candidate contradicts
      itself.

    One predicate for both is what keeps them from diverging. An older configuration
    therefore produces exactly the report it produced before #266; the ``##`` line still
    reaches ``kestrel_result.tsv``, where it cannot contradict anything, since the row it
    sits above says ``Negative`` and the line itself says the candidate is not a call.

    Args:
        report_config: The parsed ``report_config.json``.

    Returns:
        bool: Whether the subthreshold state may be used.
    """
    kestrel_logic = report_config.get("algorithm_logic", {}).get("kestrel", {})
    return SUBTHRESHOLD_RESULT in tuple(kestrel_logic.get("non_finding_results", ()))


def is_finding(result: str, default: str, non_finding: Sequence[str] = ()) -> bool:
    """Whether one algorithm's state value represents a finding.

    Derived from configuration rather than from a list written into this function, so a
    renamed ``result`` in ``report_config.json`` cannot silently invert this. That was
    ``default`` alone until #266 added a state which is neither the default nor a finding
    -- :data:`SUBTHRESHOLD_RESULT`, a negative the run can say something about. It is named
    in the block's own ``non_finding_results``, so the derivation is unchanged in kind: the
    configuration says which values are findings, and this reads it.

    ``tests/unit/test_screening_summary.py`` asserts every shipped ``result`` classifies as
    a finding unless it is declared here, and that every ``default`` does not.

    Args:
        result: The computed state value.
        default: The block's configured ``default``.
        non_finding: The block's configured ``non_finding_results``, if any. Empty
            restores the pre-#266 behaviour exactly.

    Returns:
        bool: True when the algorithm reported something.
    """
    return result not in (default, NOT_PERFORMED, "", *non_finding)


def rule_matches(current: dict[str, Any], conditions: dict[str, Any]) -> bool:
    """Whether every condition of one screening rule holds.

    Args:
        current: The three state axes.
        conditions: The rule's ``conditions`` mapping.

    Returns:
        bool: True when every named axis matches.
    """
    for key, rule_value in conditions.items():
        actual = current.get(key)
        if isinstance(rule_value, list):
            if actual not in rule_value:
                return False
        elif actual != rule_value:
            return False
    return True


def find_screening_rule(report_config: dict[str, Any], state: dict[str, Any]) -> dict[str, Any] | None:
    """Return the first configured screening rule covering ``state``.

    Args:
        report_config: The parsed ``report_config.json``.
        state: ``kestrel_result``, ``advntr_result`` and ``quality_metrics_pass``.

    Returns:
        dict[str, Any] | None: The matching rule, or None when the state has no
        message of its own and would fall back to ``screening_summary_default``.
    """
    for rule in report_config.get("screening_summary_rules", []):
        if rule_matches(state, rule.get("conditions", {})):
            return rule
    return None


def build_screening_summary(
    kestrel_df: pd.DataFrame,
    advntr_df: pd.DataFrame,
    advntr_available: bool,
    coverage_qc: CoverageQC,
    report_config: dict[str, Any],
    *,
    kestrel_execution: str | None = None,
    advntr_execution: str | None = None,
    kestrel_subthreshold: bool = False,
) -> ScreeningSummary:
    """Compute the screening state and look up its configured message.

    Args:
        kestrel_df: Kestrel results, unformatted (no HTML in the cells).
        advntr_df: adVNTR results.
        advntr_available: Whether the adVNTR stage produced a result to match on.
        coverage_qc: The coverage QC verdict, already evaluated by the caller from
            the published figures. Taking the verdict rather than the raw numbers is
            what keeps this axis and the ``coverage_qc`` column in the coverage
            summary from disagreeing - they are the same value (#172).
        report_config: The parsed ``report_config.json``.
        kestrel_execution: Whether the Kestrel stage produced a readable result. ``None``
            means "the two-state model this argument replaces", i.e. performed - which is
            all the older signature could say, so the default asserts nothing new.
        advntr_execution: The same for adVNTR. ``None`` derives it from
            ``advntr_available``, again restating the model that existed before.
        kestrel_subthreshold: Whether the Kestrel stage recorded a candidate that failed
            only the depth gate (#266). It promotes a *negative* to
            :data:`SUBTHRESHOLD_RESULT` and can do nothing else: the promotion is guarded
            on the computed result already equalling the block's ``default``, so a called
            sample is unreachable however this flag arrives.

    Returns:
        ScreeningSummary: The state and its message.
    """
    if kestrel_execution is None:
        kestrel_execution = EXECUTION_PERFORMED
    if advntr_execution is None:
        advntr_execution = EXECUTION_PERFORMED if advntr_available else EXECUTION_NOT_PERFORMED
    try:
        default_message = report_config.get("screening_summary_default", FALLBACK_SUMMARY_MESSAGE)
        algorithm_logic = report_config.get("algorithm_logic", {})
        kestrel_logic = algorithm_logic.get("kestrel", {})
        advntr_logic = algorithm_logic.get("advntr", {})

        kestrel_result = compute_algorithm_result(kestrel_df, kestrel_logic)
        advntr_result = compute_algorithm_result(advntr_df, advntr_logic) if advntr_available else NOT_PERFORMED
        logger.debug("Computed Kestrel result: %s; adVNTR result: %s", kestrel_result, advntr_result)

        kestrel_default = kestrel_logic.get("default", FALLBACK_ALGORITHM_RESULT)
        kestrel_non_finding = tuple(kestrel_logic.get("non_finding_results", ()))
        if kestrel_subthreshold and kestrel_result == kestrel_default:
            # Fail-safe against an older `report_config.json`. A deployment may ship its
            # own, and one written before #266 declares no `non_finding_results` -- under
            # it `is_finding(SUBTHRESHOLD_RESULT, "negative")` is True and the promoted
            # state would be styled as a *finding*, the exact inversion #266 forbids. So
            # promote only when the loaded configuration says the value is not a finding.
            # The `##` note in `kestrel_result.tsv` is unaffected either way: the report's
            # Kestrel section renders it independently of the screening state.
            if supports_subthreshold(report_config):
                kestrel_result = SUBTHRESHOLD_RESULT
                logger.info(
                    "Kestrel called nothing but recorded a candidate below the reporting "
                    "floor; the screening state is %s, which is not a finding.",
                    SUBTHRESHOLD_RESULT,
                )
            else:
                logger.warning(
                    "A below-reporting-floor candidate was recorded, but report_config.json "
                    "does not declare %r under algorithm_logic.kestrel.non_finding_results, "
                    "so the screening state stays %r rather than risk rendering it as a "
                    "finding.",
                    SUBTHRESHOLD_RESULT,
                    kestrel_result,
                )

        quality_metrics_pass = coverage_qc.passed

        current = {
            "kestrel_result": kestrel_result,
            "advntr_result": advntr_result,
            "quality_metrics_pass": quality_metrics_pass,
        }
        logger.debug("Unified screening conditions: %s", current)

        rule = find_screening_rule(report_config, current)
        text = rule.get("message", "") if rule is not None else ""
        segments = message_segments(rule) if rule is not None else ()
        if rule is not None:
            logger.debug("Unified rule matched: %s", rule.get("conditions"))

        matched_rule = bool(text)
        if not matched_rule:
            text = default_message
            segments = tuple(default_message.split(SEGMENT_SEPARATOR))
            logger.warning(
                "No screening rule covers kestrel_result=%r, advntr_result=%r, quality_metrics_pass=%r; "
                "the report will withhold the default message.",
                kestrel_result,
                advntr_result,
                quality_metrics_pass,
            )

        is_positive = is_finding(kestrel_result, kestrel_default, kestrel_non_finding) or is_finding(
            advntr_result,
            advntr_logic.get("default", FALLBACK_ALGORITHM_RESULT),
            advntr_logic.get("non_finding_results", ()),
        )
    except Exception as ex:
        logger.error("Exception in build_screening_summary: %s", ex)
        return ScreeningSummary(
            text=UNAVAILABLE_SUMMARY_MESSAGE,
            is_positive=False,
            kestrel_result="",
            advntr_result="",
            quality_metrics_pass=False,
            matched_rule=False,
            kestrel_execution=kestrel_execution,
            advntr_execution=advntr_execution,
            segments=(UNAVAILABLE_SUMMARY_MESSAGE,),
        )

    summary = ScreeningSummary(
        text=text,
        is_positive=is_positive,
        kestrel_result=kestrel_result,
        advntr_result=advntr_result,
        quality_metrics_pass=quality_metrics_pass,
        matched_rule=matched_rule,
        kestrel_execution=kestrel_execution,
        advntr_execution=advntr_execution,
        segments=segments,
    )
    # A message is evidence only when a rule matched a state every stage established.
    # Withholding anything else is the only honest option, and the vocabulary for saying
    # so already exists. The state axes above are left exactly as computed: they are what
    # the provenance line and tests are about, and rewriting them would hide the reason.
    if not summary.state_is_established or not summary.matched_rule:
        logger.warning(
            "Withholding the screening message: kestrel_execution=%r, advntr_execution=%r, matched_rule=%r. "
            "The computed state was not established or no configured rule described it.",
            kestrel_execution,
            advntr_execution,
            matched_rule,
        )
        summary = replace(summary, text=UNAVAILABLE_SUMMARY_MESSAGE, segments=(UNAVAILABLE_SUMMARY_MESSAGE,))

    logger.debug("Final screening summary: %s", summary.text)
    return summary
