"""
report_formatting.py

Module Purpose:
---------------
The pure presentation half of the HTML report: turning values into the strings
the Jinja2 template interpolates. No filesystem, no subprocess, no pandas I/O --
every function here takes values and returns values, which is what makes the
report's formatting testable at all.

It was extracted from ``generate_report.py`` (861 LOC, 4% covered) under
AGENTS.md rule 3. What stayed behind is the part that cannot be called without a
filesystem: loading ``pipeline_summary.json``, shelling out to ``create_report``
and writing the rendered HTML.

Four groups live here:

* **Status icons.** A value, a cutoff and a direction produce a coloured glyph
  and a colour name. The report shows eight of them and they all shared one
  ``higher_better`` flag plus two hand-rolled copies with a different treatment
  of a missing value; :func:`threshold_icon` takes that treatment as an argument
  instead.
* **Table formatting.** Column selection and renaming for the Kestrel table, the
  ``Confidence`` column's class, the ``Flag`` cell, and the per-column number
  formatting the browser used to do (#242).
* **Row counts.** The visible/total statement printed beside each results table,
  computed from the frame rather than read back out of DataTables' footer.
* **IGV fragment extraction.** ``create_report`` writes a standalone HTML page;
  the report splices three pieces out of it. The splicing is string work and is
  pure, so it lives here; opening the file does not.

Functions:
    threshold_icon: Value + cutoff to (icon, colour)
    select_display_columns: Project and rename a results frame for display
    confidence_html: Classify one ``Confidence`` value
    flag_html: One ``Flag`` value to a glyph plus its reason in words
    format_number_columns: Render a frame's numbers with the formatter each column declares
    flagged_row_count: How many rows of a results frame carry a flag
    row_count_statement: The visible/total sentence shown beside a results table
    numeric_headings: Which display headings hold numbers
    annotate_table_columns: Per-column alignment, face and column help on a rendered table
    variant_identity: What the run named, for the top of the report
    nomenclature_legend: The coded values this report prints, in words
    extract_line_after: Text following a marker, up to the next newline
    extract_igv_fragments: IGV page text to (container, tableJson, sessionDictionary)
"""

from __future__ import annotations

import html
import json
import logging
import numbers
import re
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from decimal import Decimal
from typing import Any

import pandas as pd

# Contract C1: the coverage TSV field names are declared once, by the module that
# produces them. Re-typing the strings here is how the report silently loses
# coverage when a column is renamed - `.get(name, 0)` raises nothing.
from vntyper.scripts import nomenclature

# Contract C1: the coverage TSV field names are declared once, by the module that
# produces them. Re-typing the strings here is how the report silently loses
# coverage when a column is renamed - `.get(name, 0)` raises nothing.
from vntyper.scripts.coverage_stats import _BUILD_COMPARABLE_COLUMNS, COVERAGE_COLUMNS, COVERAGE_NULL_TOKEN
from vntyper.scripts.fastp_cutoffs import calculate_passed_filter_rate_from_sources

logger = logging.getLogger(__name__)

#: Red warning triangle, shown when a metric fails its threshold.
#:
#: Both glyphs carry an accessible name (#242). A bare ``&#9888;`` in a table cell
#: is announced as its code point or skipped entirely, so the ``Status`` column of
#: the fastp table read as empty to a screen reader while looking complete on
#: screen. ``role="img"`` is what makes the browser use ``aria-label`` as the
#: element's name instead of its text. The names are deliberately about the mark
#: rather than about the direction of the comparison.
#:
#: Both colours are tokens from ``templates/_report_base.html`` rather than the CSS
#: keywords they used to be. ``red`` measured 4.00:1 against the page, under the 4.5:1
#: a foreground needs; ``green`` measured 5.14:1 and was the only pair in either report
#: that met it. A custom property resolves in an inline ``style`` exactly as it does in
#: a stylesheet, so these two fragments still carry their own presentation and now
#: follow the light, dark and print palettes with everything else.
WARNING_ICON = (
    '<span style="color:var(--state-caution);font-weight:bold;" role="img" aria-label="Warning">&#9888;</span>'
)

#: Green tick, shown when a metric passes its threshold.
OK_ICON = '<span style="color:var(--state-ok);font-weight:bold;" role="img" aria-label="No warning">&#10004;</span>'

#: What a metric with no value at all renders as. Absence is not a passing result.
MISSING_AS_BLANK: tuple[str, str] = ("", "")

#: Kestrel result column -> report heading, in display order. Columns absent from
#: the results frame are dropped, so this is a superset of any one run's output:
#: a negative run's ``kestrel_result.tsv`` carries neither ``Flag`` nor the depth
#: columns.
#:
#: The motif key is ``Motif``, not ``Motifs`` (contract C3). Both columns exist and
#: both are load-bearing upstream -- the shipped Kestrel flagging expressions
#: ``eval`` against singular ``Motif`` while duplicate ordering uses plural
#: ``Motifs``, and a missing name there evaluates to False rather than raising
#: (AGENTS.md trap 3), so neither may be renamed. ``Motifs`` is the raw
#: ``left-right`` pair Kestrel emits; ``Motif`` is the motif the variant was
#: annotated onto, which is what ``cohort_summary.py`` shows and what the heading
#: has always claimed to be.
#: ``Motif_sequence`` is **last**, and that is an observable output-format change
#: recorded in ``docs/user-guide/output-files.md``. It used to be sixth. Before its
#: value was corrected to the selected 60 bp half, the column held the 120 bp pair
#: record, whose long unbroken text pushed ``Confidence`` and ``Flag`` off the right
#: edge of the table at 1280px. Last, it is the column that scrolls.
KESTREL_DISPLAY_COLUMNS: dict[str, str] = {
    "Motif": "Motif",
    "Variant": "Variant",
    "POS": "Position",
    "REF": "REF",
    "ALT": "ALT",
    "Estimated_Depth_AlternateVariant": "Depth (Variant)",
    "Estimated_Depth_Variant_ActiveRegion": "Depth (Region)",
    "Depth_Score": "Depth Score",
    "Confidence": "Confidence",
    "Flag": "Flag",
    # MUC1 nomenclature. The heading says "name" rather than a coordinate system
    # because the tier decides how precise the available name may be.
    "Nomenclature": "MUC1 Name",
    "Nomenclature_Tier": "Tier",
    "Nomenclature_Flags": "Flags",
    # Both callers are shown, not just the verdict: where they disagree, which one
    # said what is the evidence a reader needs in order to weigh that verdict.
    "Nomenclature_Kestrel": "Kestrel Name",
    "Nomenclature_adVNTR": "adVNTR Name",
    "Ambiguity_Interval": "Ambiguity",
    "Repeat_Form": "Repeat Form",
    "Nomenclature_Note": "Naming Note",
    # Keep the long unbroken motif sequence last even as later output fields are added: the two
    # confidence/flag columns must not be displaced by the widest value (#242).
    "Motif_sequence": "Motif Sequence",
}

#: adVNTR result columns, in display order. Unlike Kestrel these are not renamed.
ADVNTR_DISPLAY_COLUMNS: tuple[str, ...] = (
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
    "Nomenclature",
    "Nomenclature_Tier",
    "Nomenclature_Flags",
    "Nomenclature_Kestrel",
    "Nomenclature_adVNTR",
    "Ambiguity_Interval",
    "Repeat_Form",
    "Nomenclature_Note",
)

#: adVNTR result column -> report heading. Applied to the *rendered* frame only.
#:
#: The two results tables sat one above the other in the same document with two
#: different conventions for naming a column: Kestrel's were English (``Depth
#: (Variant)``, ``MUC1 Name``, ``Naming Note``) and adVNTR's were the dataframe
#: identifiers as written (``NumberOfSupportingReads``, ``Nomenclature_Tier``,
#: ``Ambiguity_Interval``). It read as two maturity levels of the same report stapled
#: together, and it made the nomenclature columns - which are the same eight fields in
#: both tables - look like different fields.
#:
#: :data:`ADVNTR_DISPLAY_COLUMNS` is unchanged and stays the projection order, and
#: ``advntr_df`` itself keeps the source names because the screening rules match on
#: them (``VID``, ``Flag``). This renames the copy that is rendered, and nothing else:
#: it is an HTML display change, and the TSV remains the parsing source.
ADVNTR_DISPLAY_HEADINGS: dict[str, str] = {
    "VID": "VID",
    "Variant": "Variant",
    "NumberOfSupportingReads": "Supporting Reads",
    "MeanCoverage": "Mean Coverage",
    "Pvalue": "P-value",
    "RU": "Repeat Unit",
    "POS": "Position",
    "REF": "REF",
    "ALT": "ALT",
    "Flag": "Flag",
    "Nomenclature": "MUC1 Name",
    "Nomenclature_Tier": "Tier",
    "Nomenclature_Flags": "Flags",
    "Nomenclature_Kestrel": "Kestrel Name",
    "Nomenclature_adVNTR": "adVNTR Name",
    "Ambiguity_Interval": "Ambiguity",
    "Repeat_Form": "Repeat Form",
    "Nomenclature_Note": "Naming Note",
}

#: How a displayed value is turned into text. Every column of both results tables
#: names one of these, so adding a display column without deciding how its value is
#: rendered fails in :func:`format_number_columns` rather than silently inheriting
#: whatever pandas would have printed.
#:
#: Until #242 there was no such decision to make: the template shipped a jQuery
#: ``applyRounding()`` that rewrote every numeric cell of every initialised table to
#: four decimal places and stripped trailing zeroes. It ran in the reader's browser,
#: so the figure on screen depended on whether three CDNs resolved, and it was
#: table-agnostic, so one rule governed a genomic position, a depth score and a
#: p-value alike.
FORMAT_TEXT = "text"
FORMAT_INTEGER = "integer"
FORMAT_TWO_DECIMALS = "two-decimals"
FORMAT_SIX_DECIMALS = "six-decimals"
FORMAT_THREE_SIGNIFICANT = "three-significant"

#: Kestrel result column -> how its value is rendered. Keys are exactly
#: :data:`KESTREL_DISPLAY_COLUMNS`'s.
#:
#: ``Depth_Score`` is six decimal places rather than the browser's four: the
#: confidence calibration it is judged against is stated to five
#: (``kestrel_config.json``: 0.00469 and 0.00515), and four decimal places printed a
#: score of 1.234e-05 as the value ``0``.
KESTREL_CELL_FORMATS: dict[str, str] = {
    "Motif": FORMAT_TEXT,
    "Variant": FORMAT_TEXT,
    "POS": FORMAT_INTEGER,
    "REF": FORMAT_TEXT,
    "ALT": FORMAT_TEXT,
    "Estimated_Depth_AlternateVariant": FORMAT_INTEGER,
    "Estimated_Depth_Variant_ActiveRegion": FORMAT_INTEGER,
    "Depth_Score": FORMAT_SIX_DECIMALS,
    "Confidence": FORMAT_TEXT,
    "Flag": FORMAT_TEXT,
    "Nomenclature": FORMAT_TEXT,
    "Nomenclature_Tier": FORMAT_TEXT,
    "Nomenclature_Flags": FORMAT_TEXT,
    "Nomenclature_Kestrel": FORMAT_TEXT,
    "Nomenclature_adVNTR": FORMAT_TEXT,
    "Ambiguity_Interval": FORMAT_TEXT,
    "Repeat_Form": FORMAT_TEXT,
    "Nomenclature_Note": FORMAT_TEXT,
    "Motif_sequence": FORMAT_TEXT,
}

#: The same decisions keyed by the heading each Kestrel column is renamed to, which
#: is what the display frame carries by the time it is formatted.
KESTREL_DISPLAY_CELL_FORMATS: dict[str, str] = {
    KESTREL_DISPLAY_COLUMNS[column]: cell_format for column, cell_format in KESTREL_CELL_FORMATS.items()
}


#: adVNTR result column -> how its value is rendered. Keys are exactly
#: :data:`ADVNTR_DISPLAY_COLUMNS`'s.
#:
#: ``Pvalue`` is three significant figures rather than four decimal places, because
#: ``parseFloat((1e-9).toFixed(4)).toString()`` is ``"0"`` -- the online report
#: displayed a highly significant p-value as zero. Significant figures keep the
#: magnitude of a value whose whole meaning is its magnitude.
ADVNTR_CELL_FORMATS: dict[str, str] = {
    "VID": FORMAT_INTEGER,
    "Variant": FORMAT_TEXT,
    "NumberOfSupportingReads": FORMAT_INTEGER,
    "MeanCoverage": FORMAT_TWO_DECIMALS,
    "Pvalue": FORMAT_THREE_SIGNIFICANT,
    "RU": FORMAT_TEXT,
    "POS": FORMAT_INTEGER,
    "REF": FORMAT_TEXT,
    "ALT": FORMAT_TEXT,
    "Flag": FORMAT_TEXT,
    "Nomenclature": FORMAT_TEXT,
    "Nomenclature_Tier": FORMAT_TEXT,
    "Nomenclature_Flags": FORMAT_TEXT,
    "Nomenclature_Kestrel": FORMAT_TEXT,
    "Nomenclature_adVNTR": FORMAT_TEXT,
    "Ambiguity_Interval": FORMAT_TEXT,
    "Repeat_Form": FORMAT_TEXT,
    "Nomenclature_Note": FORMAT_TEXT,
}

#: What a nomenclature tier means, keyed by the letter the results frame carries.
#:
#: The tier is the single most compressed thing either results table says, it is a bare
#: letter, and until now the report printed it with nothing anywhere saying what it
#: meant - a reader who did not already know had no way to find out from the artefact in
#: front of them. Each entry is ``(short label, what it took to earn it)``, and both are
#: written from :func:`vntyper.scripts.nomenclature.reconcile`'s own conditions rather
#: than paraphrased: tier A requires two independent sources agreeing after
#: normalisation, a name present in ``KNOWN_VARIANTS``, support at or above
#: ``MIN_SUPPORT_FOR_TIER_A``, and none of the context or disagreement flags.
#:
#: The wording states what was established, never how much to trust the sample. A tier
#: is a statement about how far a *name* has been checked.
NOMENCLATURE_TIERS: dict[str, tuple[str, str]] = {
    "A": (
        "corroborated name",
        "Two independent callers name the same allele, that name matches a MUC1 variant "
        "described in the literature, read support meets the threshold, and no context or "
        "disagreement flag is set. Tier A is the only tier that states a bare position.",
    ),
    "B": (
        "qualified name",
        "A name was computed, but at least one tier-A condition is unmet - a single caller, "
        "an allele nobody has described before, thin read support, a diverging motif context, "
        "or the two callers disagreeing. The name is shown so it can be weighed; the flags "
        "beside it say which condition was missed.",
    ),
    "C": (
        "no position",
        "No coordinate could be computed. What is stated is the event and its net length "
        "change, without a position - not a negative, and not a name withheld.",
    ),
}

#: Which flag evidences which unmet tier-A condition, in the words
#: :func:`vntyper.scripts.nomenclature.reconcile` decides them by.
#:
#: **This is the question the tier definition alone cannot answer.** A tier-B call is
#: one where "at least one tier-A condition is unmet", and the report printed the
#: letter, printed the flags, and left the reader to work out which of five conditions
#: applied to *this* call. The person who wrote the caller asked the question, which is
#: the clearest possible evidence that a reader could not answer it from the artefact.
#:
#: The keys are the flag constants that :func:`reconcile` tests by name when it decides
#: tier A, plus :data:`~vntyper.scripts.nomenclature.FLAG_REPRESENTATION_ONLY`, which is
#: set on exactly the names absent from ``KNOWN_VARIANTS`` - the condition
#: ``chosen_name in KNOWN_VARIANTS`` fails on. They are imported rather than retyped, so
#: a renamed flag fails at import instead of silently explaining nothing.
#:
#: Two of tier A's conditions leave no flag behind and are therefore deliberately not
#: inferred here: a single uncorroborated caller, and read support that is *unknown*
#: rather than low. When the flags evidence no blocker, the report says the tier is B
#: and does not guess why - which is the honest reading of a call whose reason is not in
#: the record.
TIER_A_BLOCKERS: dict[str, str] = {
    nomenclature.FLAG_MOTIF_CONTEXT_DIVERGES: (
        "the motif context around the call diverges from the canonical X repeat unit"
    ),
    nomenclature.FLAG_SEQUENCE_UNDETERMINED: "the inserted or deleted sequence could not be determined",
    nomenclature.FLAG_CALLER_DISAGREEMENT: "the two callers did not name the same allele",
    nomenclature.FLAG_LOW_READ_SUPPORT: (
        f"read support is below the {nomenclature.MIN_SUPPORT_FOR_TIER_A} reads the corroborated tier requires"
    ),
    nomenclature.FLAG_REPRESENTATION_ONLY: "no MUC1 variant described in the literature matches this name",
}

#: The closed nomenclature flag vocabulary, in words.
#:
#: The keys are exactly the kebab-case constants declared in
#: :mod:`vntyper.scripts.nomenclature`; that module states the vocabulary is closed and
#: that a consumer may match on it, which is what makes writing the meanings out here
#: safe. ``length-truncated`` is declared there and not yet set by any path: it is
#: documented anyway, because a reserved token a reader may one day meet is worth a
#: sentence, and an unexplained one is the defect this table exists to close.
NOMENCLATURE_FLAG_MEANINGS: dict[str, str] = {
    "position-ambiguous": (
        "The edit can be written at more than one position in the repeat unit. The ambiguity "
        "interval gives the range over which the placements are indistinguishable."
    ),
    "spans-unit-junction": (
        "The event crosses the boundary between two repeat units, so it cannot be expressed "
        "as one span on a single unit."
    ),
    "motif-context-diverges": (
        "The sequence around the call in the motif the caller assigned differs from the "
        "canonical X unit, so the coordinate projected onto X is less certain."
    ),
    "allele-unrepresentable-in-vcf": (
        "The allele cannot be written in the caller's VCF shape. The name comes from the "
        "reads, which are the better evidence for this locus."
    ),
    "low-read-support": "Fewer reads support the call than the corroborated tier requires.",
    "caller-disagreement": "Kestrel and adVNTR did not name the same allele.",
    "length-truncated": "The reported length is a lower bound rather than the full extent.",
    "sequence-undetermined": (
        "The inserted or deleted sequence itself could not be determined, so no position is given."
    ),
    "known-variant": (
        "The name matches a MUC1 variant described in the literature. The table is used to "
        "check a name, never to produce one - a match says somebody has described this allele "
        "before, which is weaker than saying the call is correct."
    ),
    "representation-of-caller-call": (
        "The name represents what the caller reported and matches no described variant. It requires validation."
    ),
}

#: One line per results-table heading, shown as the column's own hover text and printed
#: in the reading key beneath each table.
#:
#: Keyed by the *heading*, which is what both frames carry by the time they are rendered
#: - the Kestrel columns are renamed by :data:`KESTREL_DISPLAY_COLUMNS` and the adVNTR
#: ones are not, so a shared key would have to be one or the other and could only ever
#: explain half the report.
COLUMN_HELP: dict[str, str] = {
    "Supporting Reads": "Reads adVNTR counted in support of this call.",
    "Mean Coverage": "Mean read depth adVNTR measured over the locus. Not the region mean above.",
    "P-value": "adVNTR's significance for the call, to three significant figures.",
    "Repeat Unit": "The adVNTR repeat unit the call was made in.",
    "Motif": "The MUC1 repeat motif this variant was annotated onto.",
    "Motifs": "The raw left-right motif pair Kestrel emitted for the record.",
    "Variant": "The event class the caller reported.",
    "Position": "Position within the 120 bp Kestrel motif pair, not a genomic coordinate.",
    "POS": "Position within the caller's own coordinate frame, not a genomic coordinate.",
    "REF": "The reference allele at the reported position.",
    "ALT": "The alternate allele the reads support.",
    "Depth (Variant)": "Reads Kestrel estimates support the alternate allele.",
    "Depth (Region)": "Reads Kestrel estimates over the active region around the call.",
    "Depth Score": (
        "Variant depth over region depth. It scales with the inverse of the array length, so "
        "it is comparable within a sample and not between assemblies."
    ),
    "Confidence": "Kestrel's own calibrated label for the call.",
    "Flag": "Whether a configured flagging rule fired on this row, and which one.",
    "MUC1 Name": "The reconciled MUC1 name for the allele, on the canonical X repeat unit.",
    "Nomenclature": "The reconciled MUC1 name for the allele, on the canonical X repeat unit.",
    "Tier": "How far the name beside it has been checked. See the reading key below the table.",
    "Nomenclature_Tier": "How far the name beside it has been checked. See the reading key below the table.",
    "Flags": "Closed-vocabulary qualifiers on the name. Each one present is spelled out below the table.",
    "Nomenclature_Flags": (
        "Closed-vocabulary qualifiers on the name. Each one present is spelled out below the table."
    ),
    "Kestrel Name": "The name derived from Kestrel's record alone, before reconciliation.",
    "Nomenclature_Kestrel": "The name derived from Kestrel's record alone, before reconciliation.",
    "adVNTR Name": "The name derived from adVNTR's record alone, before reconciliation.",
    "Nomenclature_adVNTR": "The name derived from adVNTR's record alone, before reconciliation.",
    "Ambiguity": "The range of positions over which this edit is indistinguishable.",
    "Ambiguity_Interval": "The range of positions over which this edit is indistinguishable.",
    "Repeat Form": "The affected tract before and after the edit, written as a repeat count.",
    "Repeat_Form": "The affected tract before and after the edit, written as a repeat count.",
    "Naming Note": "What the name is, and how far it has been checked.",
    "Nomenclature_Note": "What the name is, and how far it has been checked.",
    "Motif Sequence": "The 60 bp repeat-unit half named by the Motif column.",
    "VID": "adVNTR's own variant identifier.",
    "NumberOfSupportingReads": "Reads adVNTR counted in support of this call.",
    "MeanCoverage": "Mean read depth adVNTR measured over the locus.",
    "Pvalue": "adVNTR's significance for the call, to three significant figures.",
    "RU": "The adVNTR repeat unit the call was made in.",
}

#: Headings whose values are identifiers, sequences or coordinates rather than prose,
#: set in the mono face so a reader can compare them character by character.
MONO_COLUMNS: frozenset[str] = frozenset(
    {
        "REF",
        "ALT",
        "MUC1 Name",
        "Nomenclature",
        "Kestrel Name",
        "Nomenclature_Kestrel",
        "adVNTR Name",
        "Nomenclature_adVNTR",
        "Ambiguity",
        "Ambiguity_Interval",
        "Repeat Form",
        "Repeat_Form",
        "Motif Sequence",
        "RU",
    }
)

#: The columns a results table shows before the reader asks for the rest.
#:
#: **This hides columns. It never hides a row, and the distinction is the whole design.**
#: A hidden row is evidence the reader cannot know exists - the defect #242 is named
#: after. A column set to ``display: none`` by a control that says so, defaults to the
#: shorter set, states the full count on its own label, prints in full whichever way it
#: is left, and shows *everything* when no script runs, is the opposite: it is the
#: reader choosing how much of one row to look at.
#:
#: What is left out is chosen so that nothing leaves the page. Every optional column
#: here is either already stated in the masthead's allele panel for a run that named one
#: allele (the ambiguity interval, the repeat form, both callers' own names, the naming
#: note) or is a figure the column beside it derives (``Depth (Region)``, which
#: ``Depth Score`` is the ratio over) or the 60 bp sequence, which is the widest value
#: in the document and the least scanned. Nineteen columns measured 1,946px inside a
#: 1,130px frame; twelve fit.
KESTREL_ESSENTIAL_COLUMNS: frozenset[str] = frozenset(
    {
        "Motif",
        "Variant",
        "Position",
        "REF",
        "ALT",
        "Depth (Variant)",
        "Depth Score",
        "Confidence",
        "Flag",
        "MUC1 Name",
        "Tier",
        "Flags",
    }
)

#: The same decision for the adVNTR table, in its own headings.
ADVNTR_ESSENTIAL_COLUMNS: frozenset[str] = frozenset(
    {
        "VID",
        "Variant",
        "Supporting Reads",
        "Mean Coverage",
        "P-value",
        "Position",
        "REF",
        "ALT",
        "Flag",
        "MUC1 Name",
        "Tier",
        "Flags",
    }
)

#: The columns holding an unbroken biological sequence, which is the only kind of value
#: allowed to break inside itself when printed.
SEQUENCE_COLUMNS: frozenset[str] = frozenset({"Motif Sequence", "Motif_sequence"})

#: Headings whose values are several tokens or a sentence, given a floor so the
#: auto-layout cannot squeeze them into a column narrower than one token. The flag list
#: is the reason: at its natural width it wrapped to six lines and made a one-variant
#: row 244px tall, with every other cell in that row sitting over 180px of white.
WIDE_COLUMNS: frozenset[str] = frozenset({"Flags", "Naming Note", "Nomenclature_Note", "Confidence"})

#: The heading that absorbs the table's slack width, per results table. ``.col-grow``
#: has been declared in the shared stylesheet since #242 and emitted by nothing, so
#: every column sized to its content and the table never filled its container.
GROW_COLUMNS: frozenset[str] = frozenset({"Motif Sequence", "Naming Note", "Nomenclature_Note"})

#: ``Flag`` values that mean "nothing to report". These are the three the shipped
#: client-side code treated as clean, kept verbatim so the glyph a reader sees does
#: not change meaning with this release.
FLAG_CLEAN_VALUES: tuple[str, ...] = ("Not flagged", "Not applicable", "")

#: Heavy check mark, shown beside a clean ``Flag``.
FLAG_OK_GLYPH = "&#10004;"

#: Heavy multiplication X, shown beside a flagged one.
FLAG_WARNING_GLYPH = "&#10006;"

#: Classes on the ``Flag`` cell's wrapper. ``FLAG_FLAGGED_CLASS`` is what the
#: "Highlight flagged values" switch styles; it is emphasis only, and no rule
#: anywhere may use it to remove a row (#242).
FLAG_FLAGGED_CLASS = "flag-flagged"
FLAG_CLEAN_CLASS = "flag-clean"

#: CSS classes per ``Confidence`` value. A value not listed renders as bare text.
#:
#: **There is no hue here, and that is the change.** This used to be an inline style
#: per value: ``High_Precision`` - the *most* trustworthy call the pipeline makes - was
#: ``color:red``, one column away from a red ``Flag`` glyph that means the row is *not*
#: to be trusted, and it measured 4.00:1 against white and 3.57:1 against a striped
#: cohort row. ``Low_Precision`` was ``color:orange`` at 1.97:1 and 1.76:1. A
#: transitional underline was added so the two were not separated by hue alone; both
#: the hue and the underline go together here, because the honest channel was always
#: the label, and the label is text in the cell in every branch.
#:
#: What is left is a class. The report's stylesheet renders it as plain text in the
#: page's own ink (``templates/_report_base.html``), and red now means "something is
#: wrong" and nothing else, anywhere on the page. A deployment that wants a treatment
#: has a hook to hang one on; VNtyper ships none.
CONFIDENCE_CLASSES: dict[str, str] = {
    "Low_Precision": "confidence confidence-low-precision",
    "High_Precision": "confidence confidence-high-precision",
    "High_Precision*": "confidence confidence-high-precision",
}

#: And the same for adVNTR, keyed by :data:`ADVNTR_DISPLAY_HEADINGS`. The rendered frame
#: is renamed after it is formatted, so the alignment decision has to be looked up under
#: the heading rather than under the source name.
ADVNTR_DISPLAY_CELL_FORMATS: dict[str, str] = {
    ADVNTR_DISPLAY_HEADINGS[column]: cell_format for column, cell_format in ADVNTR_CELL_FORMATS.items()
}


#: How each coverage field is coerced for display. The keys must be exactly
#: :data:`~vntyper.scripts.coverage_stats.COVERAGE_COLUMNS`; a contract test
#: enforces that, so adding a coverage column fails here rather than in the HTML.
COVERAGE_FIELD_TYPES: dict[str, type] = {
    "mean": float,
    "median": float,
    "stdev": float,
    "min": int,
    "max": int,
    "region_length": int,
    "uncovered_bases": int,
    "percent_uncovered": float,
    # #222's build-comparable columns. `depth_counting_policy` is a token, not a number:
    # the array sum is only comparable across builds under the policy it names.
    "vntr_array_length": int,
    "vntr_array_depth_sum": int,
    "vntr_array_depth_sum_per_unit_length": float,
    "depth_sum_reference_length": int,
    "vntr_flank_bases": int,
    "vntr_flank_mean_depth": float,
    "depth_counting_policy": str,
    # The QC verdict is a label, not a measurement: `str` keeps "PASS"/"FAIL" intact
    # and keeps it out of the two-decimal rendering the numeric fields get (#172).
    "coverage_qc": str,
}

#: Markers delimiting the three fragments spliced out of the IGV report page.
IGV_CONTAINER_MARKER = '<div id="container"'
IGV_BODY_END_MARKER = "</body>"
IGV_TABLE_JSON_MARKER = "const tableJson = "
IGV_SESSION_DICTIONARY_MARKER = "const sessionDictionary = "

#: Valid JavaScript literals for a report with no IGV panel. The template
#: interpolates the extracted fragments directly into a ``<script>`` block, so an
#: empty one is a syntax error rather than an empty table.
EMPTY_TABLE_JSON = '{"headers": [], "rows": []}'
EMPTY_SESSION_DICTIONARY = "{}"


def threshold_icon(
    value: float | Decimal | None,
    cutoff: float | Decimal,
    *,
    higher_better: bool = True,
    on_missing: tuple[str, str] = MISSING_AS_BLANK,
) -> tuple[str, str]:
    """Return the status glyph and colour for one metric.

    Args:
        value: The measured value, or None when the metric was not computed.
        cutoff: The threshold to compare against.
        higher_better: True when a value below ``cutoff`` is the problem; False
            when a value above it is.
        on_missing: What to return when ``value`` is None.

    Returns:
        tuple[str, str]: An HTML glyph and a CSS colour name.
    """
    if value is None:
        logger.debug("threshold_icon called with value=None; returning %s.", on_missing)
        return on_missing
    failed = value < cutoff if higher_better else value > cutoff
    logger.debug(
        "Value %s against cutoff %s (higher_better=%s): %s.",
        value,
        cutoff,
        higher_better,
        "failed" if failed else "passed",
    )
    return (WARNING_ICON, "red") if failed else (OK_ICON, "green")


def select_display_columns(df: pd.DataFrame, columns: dict[str, str]) -> pd.DataFrame:
    """Project a results frame onto its display columns and rename them.

    Columns absent from ``df`` are skipped rather than raising, because the
    Kestrel result schema differs between a positive and a negative run.

    Args:
        df: The results frame.
        columns: Source column name -> display heading, in display order.

    Returns:
        pd.DataFrame: A new frame carrying only the columns that were present,
        renamed and in the order ``columns`` declares.
    """
    present = [column for column in columns if column in df.columns]
    missing = [column for column in columns if column not in df.columns]
    if missing:
        logger.debug("Display columns absent from the results frame: %s", missing)
    return df[present].rename(columns={column: columns[column] for column in present})


def confidence_html(value: Any) -> str:
    """Classify one ``Confidence`` value for display.

    The label is escaped either way. This is the one cell of the Kestrel table
    that legitimately carries markup, which is why the table is rendered with
    ``escape=False``; a value with no configured class therefore used to reach
    the HTML untouched.

    Args:
        value: The confidence label.

    Returns:
        str: The escaped label, wrapped in a classed span when it has a
        configured class. The span carries no colour: see
        :data:`CONFIDENCE_CLASSES`.
    """
    key = value if isinstance(value, str) else str(value)
    text = escape_html(key)
    css_class = CONFIDENCE_CLASSES.get(key)
    if css_class is None:
        return text
    return f'<span class="{css_class}">{text}</span>'


def flag_html(value: Any) -> str:
    """Render one ``Flag`` value as a glyph followed by its reason in words.

    The reason used to be visible only as a Bootstrap tooltip built in the browser:
    ``updateFlagColumn`` replaced the cell with a bare tick or cross and put the
    reason in a ``title`` attribute. That made it invisible in print, invisible to a
    screen reader once Bootstrap moves ``title`` into ``data-original-title``, and
    absent entirely when the script did not run (#242). Rendering it here means the
    cell says why with no script at all.

    The value is sample-derived -- ``vntyper report`` and ``vntyper cohort`` both
    consume a supplied ``pipeline_summary.json`` (#207) -- so it is escaped before it
    becomes markup, and the column it goes into must be named in ``html_columns``
    when the table is rendered so that :func:`escape_frame_cells` leaves it alone
    and escapes everything else.

    The glyph used to carry ``style="color:red"`` / ``color:green`` of its own. Both
    are now taken from the wrapper's class in ``templates/_report_base.html``, so the
    two marks follow the light, dark and print palettes instead of measuring 4.00:1
    against the page in the one that mattered.

    Args:
        value: The ``Flag`` cell value.

    Returns:
        str: A span carrying the glyph and the escaped reason.
    """
    key = value if isinstance(value, str) else str(value)
    clean = key in FLAG_CLEAN_VALUES
    glyph = FLAG_OK_GLYPH if clean else FLAG_WARNING_GLYPH
    css_class = FLAG_CLEAN_CLASS if clean else FLAG_FLAGGED_CLASS
    return (
        f'<span class="{css_class}">'
        f'<span class="flag-glyph" aria-hidden="true">{glyph}</span> '
        f'<span class="flag-reason">{escape_html(key)}</span>'
        f"</span>"
    )


def _numeric(value: Any) -> float | None:
    """Coerce one cell to a float, or None when it is not a number.

    ``numbers.Real`` rather than ``(int, float)``: pandas hands out numpy scalars,
    and ``numpy.int64`` is not a Python ``int`` on this platform. ``bool`` is
    excluded explicitly because it *is* a ``numbers.Integral`` and rendering True as
    ``1`` would be a silent type change.

    Args:
        value: The cell value.

    Returns:
        float | None: The number, or None for a missing or non-numeric value.
    """
    if isinstance(value, bool):
        return None
    if isinstance(value, str):
        try:
            return float(value.strip())
        except ValueError:
            return None
    if isinstance(value, numbers.Real):
        return None if pd.isna(value) else float(value)
    return None


def _as_text(value: Any) -> Any:
    """Leave a value exactly as it is."""
    return value


def _as_integer(value: Any) -> Any:
    """Render a number with no decimal part; pass anything else through."""
    number = _numeric(value)
    return value if number is None else str(int(round(number)))


def _as_two_decimals(value: Any) -> Any:
    """Render a number to two decimal places, trailing zeroes kept."""
    number = _numeric(value)
    return value if number is None else f"{number:.2f}"


def _as_six_decimals(value: Any) -> Any:
    """Render a number to six decimal places, trailing zeroes kept."""
    number = _numeric(value)
    return value if number is None else f"{number:.6f}"


def _as_three_significant(value: Any) -> Any:
    """Render a number to three significant figures, in scientific form when small."""
    number = _numeric(value)
    return value if number is None else f"{number:.3g}"


#: Format name -> the function that applies it.
CELL_FORMATTERS: dict[str, Callable[[Any], Any]] = {
    FORMAT_TEXT: _as_text,
    FORMAT_INTEGER: _as_integer,
    FORMAT_TWO_DECIMALS: _as_two_decimals,
    FORMAT_SIX_DECIMALS: _as_six_decimals,
    FORMAT_THREE_SIGNIFICANT: _as_three_significant,
}


def format_number_columns(df: pd.DataFrame, formats: Mapping[str, str]) -> pd.DataFrame:
    """Render every cell of a results frame with the formatter its column declares.

    Args:
        df: The projected display frame. Not modified.
        formats: Column -> format name. Must cover every column of ``df``.

    Returns:
        pd.DataFrame: A new frame whose formatted columns hold display strings. A
        value that is not a number - the string ``"None"`` a negative Kestrel run
        writes into every numeric column, or a missing value - passes through
        untouched.

    Raises:
        ValueError: If ``df`` carries a column ``formats`` says nothing about. This
            is the point of declaring the table per column: a new column that
            inherited a default would be rendered at whatever precision pandas chose
            for the rest of its column, which is how ``40.0`` came out as ``40.00``
            and ``1e-9`` as ``1.000000e-09`` in the same table.
    """
    undeclared = [str(column) for column in df.columns if column not in formats]
    if undeclared:
        msg = f"No display format declared for report column(s): {', '.join(sorted(undeclared))}"
        logger.error(msg)
        raise ValueError(msg)

    formatted = df.copy()
    for column, name in formats.items():
        if column in formatted.columns:
            formatted[column] = formatted[column].map(CELL_FORMATTERS[name])
    return formatted


#: What ``kestrel_genotyping.output_empty_result`` writes into every column but
#: ``Confidence`` when Kestrel ran and called nothing: the literal four characters
#: ``None``, not a JSON null. ``summary.parse_tsv`` splits the file on tabs and coerces
#: nothing, so this is what reaches the report (AGENTS.md trap 5).
EMPTY_RESULT_TOKEN = "None"

#: The ``Confidence`` that placeholder row carries. It is not a confidence value - no
#: rule in ``report_config.json`` matches it - it is the row saying there is no call.
EMPTY_RESULT_CONFIDENCE = "Negative"


def _names_nothing(value: Any) -> bool:
    """Whether one cell of a candidate placeholder row states no fact.

    Args:
        value: The cell value.

    Returns:
        bool: True for the placeholder token, an empty string, ``None`` and NaN.
    """
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip() in ("", EMPTY_RESULT_TOKEN)
    return isinstance(value, float) and pd.isna(value)


def is_empty_result_row(row: Mapping[str, Any]) -> bool:
    """Whether one Kestrel row is the placeholder a negative run writes, not a variant.

    ``output_empty_result`` writes one row of ``"None"`` with ``Confidence`` set to
    ``"Negative"`` so that ``kestrel_result.tsv`` has a body rather than a header alone.
    Rendered as a table row it read ``None None None None None None None None NaN
    Negative`` under the heading "Kestrel Identified Variants", which is what the
    commonest report in any cohort said (#242).

    **Both halves of the test are required.** The ``Confidence`` value alone is not
    enough: a report is a record, and a rule that removed every row carrying that label
    would remove a real call the moment a future calibration used the word. So the row
    must also name nothing - no motif, no position, no REF, no ALT.

    Args:
        row: One row of the ``Kestrel Genotyping`` step's data.

    Returns:
        bool: True when the row is the empty-result placeholder.
    """
    confidence = row.get("Confidence")
    if not isinstance(confidence, str) or confidence.strip() != EMPTY_RESULT_CONFIDENCE:
        return False
    return all(_names_nothing(value) for key, value in row.items() if key != "Confidence")


def drop_empty_result_rows(rows: Sequence[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    """Return the rows that state a result, in the order they arrived.

    Args:
        rows: The ``Kestrel Genotyping`` step's data. Not modified.

    Returns:
        list[Mapping[str, Any]]: Every row that is not the empty-result placeholder.
    """
    kept = [row for row in rows if not is_empty_result_row(row)]
    if len(kept) != len(rows):
        logger.info("Suppressed %d Kestrel empty-result placeholder row(s) from display.", len(rows) - len(kept))
    return kept


def flagged_row_count(df: pd.DataFrame) -> int:
    """Count the rows of a results frame whose ``Flag`` states a reason.

    Args:
        df: A results frame, formatted or not, possibly without a ``Flag`` column -
            a Kestrel run with no flagging rules configured produces none
            (AGENTS.md trap 4).

    Returns:
        int: The number of flagged rows; 0 when the column is absent.
    """
    if df.empty or "Flag" not in df.columns:
        return 0
    return sum(
        1
        for value in df["Flag"]
        if not (isinstance(value, float) and pd.isna(value))
        and (value if isinstance(value, str) else str(value)) not in FLAG_CLEAN_VALUES
    )


def row_count_statement(total: int, flagged: int, *, noun: str) -> str:
    """State how many rows of a results table the reader is being shown.

    Every row the pipeline wrote is rendered, so the visible count and the total are
    the same number - and printing both is the point: a reader who has been handed a
    report with rows silently removed has no way to know, and this sentence is the
    one that says nothing was withheld. It is computed here rather than read out of
    DataTables' "Showing 1 to 3 of 3 entries" footer, which only exists when the
    reader's browser could reach a CDN (#242).

    Args:
        total: Rows in the frame.
        flagged: How many of them carry a flag.
        noun: The algorithm's name, e.g. ``Kestrel``.

    Returns:
        str: The sentence to print beside the table.
    """
    rows = "row" if total == 1 else "rows"
    count = str(flagged) if flagged else "none"
    return f"Showing {total} of {total} {noun} {rows}; {count} flagged."


def escape_html(value: Any) -> str:
    """HTML-escape one value, including quotes.

    Args:
        value: Any value; non-strings are stringified first.

    Returns:
        str: The escaped text.
    """
    return html.escape(value if isinstance(value, str) else str(value), quote=True)


def escape_frame_cells(df: pd.DataFrame, html_columns: tuple[str, ...] = ()) -> pd.DataFrame:
    """Escape every string cell of a frame that is about to be rendered raw.

    Only ``str`` cells are touched. Numbers and NaN cannot carry markup, and
    stringifying them here would take pandas' own float and NA formatting out of
    the rendered table.

    Args:
        df: The frame to escape.
        html_columns: Columns already holding markup we constructed ourselves,
            left alone.

    Returns:
        pd.DataFrame: A new frame; the input is not modified.
    """
    escaped = df.copy()
    for column in escaped.columns:
        if column in html_columns:
            continue
        escaped[column] = escaped[column].map(
            lambda cell: html.escape(cell, quote=True) if isinstance(cell, str) else cell
        )
    return escaped


def escaped_table_html(
    df: pd.DataFrame,
    classes: str,
    html_columns: tuple[str, ...] = (),
    table_id: str | None = None,
) -> str:
    """Render a frame as an HTML table with every sample-derived cell escaped.

    ``DataFrame.to_html(escape=False)`` is needed whenever *any* column holds markup
    VNtyper built, and it disables escaping for the whole table - so the columns holding
    a sample's own strings go out as HTML too. This pairs it with
    :func:`escape_frame_cells`, which escapes everything except the columns named as
    ours, so the exemption is stated per column instead of applying to all of them.

    Args:
        df: The frame to render. Not modified.
        classes: The CSS classes for the ``<table>`` element.
        html_columns: Columns already holding markup this codebase constructed.
            Anything not named here is escaped.
        table_id: The ``id`` for the ``<table>`` element, when something addresses it.
            The per-sample Kestrel table is ``#kestrel_table`` to the template's
            DataTables initialisation and to the browser tier; the cohort tables are
            addressed by class and pass nothing.

    Returns:
        str: The table markup, or "" for an empty frame - ``to_html`` on one produces a
            headerless table that renders as a stray empty box. Callers treat that empty
            string as the hook for an authored empty state (#242).
    """
    if df.empty:
        return ""
    return escape_frame_cells(df, html_columns=html_columns).to_html(
        table_id=table_id,
        classes=classes,
        index=False,
        escape=False,
        # `border=0`, because pandas' default writes `border="1"` into the markup and a
        # presentation attribute is a colour no stylesheet can see: it rendered as a
        # 1px `inset` grey (3.95:1) on every cell, square-cornered, inside the report's
        # own hairline frame with its 6px radius. The report's rules paint the table.
        border=0,
    )


#: Matches the three table tags whose column position this module tracks. ``<table>``,
#: ``<thead>`` and ``<tbody>`` do not match: ``\b`` after ``th`` fails against ``thead``,
#: and neither ``tr`` nor ``td`` is a prefix of the other two.
_TABLE_TAG = re.compile(r"<(tr|td|th)\b([^>]*)>", re.IGNORECASE)

#: The inline alignment pandas puts on its header row. It is dropped rather than
#: fought: alignment is decided per column below, and an inline style on the row is a
#: declaration this stylesheet cannot see.
_PANDAS_HEADER_ALIGN = re.compile(r'\s*style="text-align:\s*right;?"')

#: Cell formats whose values are numbers, so their column is right-aligned on the
#: figures' own axis rather than centred. ``FORMAT_TEXT`` is everything else.
_NUMERIC_FORMATS: frozenset[str] = frozenset(
    {FORMAT_INTEGER, FORMAT_TWO_DECIMALS, FORMAT_SIX_DECIMALS, FORMAT_THREE_SIGNIFICANT}
)

#: The nomenclature fields a results row carries, and the key each becomes.
#:
#: Two column names per field, because the two frames name them differently and both
#: are read here. ``build_kestrel_frames`` projects through
#: :data:`KESTREL_DISPLAY_COLUMNS`, so its matching frame already carries the display
#: headings; ``build_advntr_frame`` projects through :data:`ADVNTR_DISPLAY_COLUMNS`,
#: which renames nothing. Reading one spelling would have summarised whichever caller
#: happened to use it and silently ignored the other.
_IDENTITY_FIELDS: dict[str, tuple[str, ...]] = {
    "name": ("Nomenclature", "MUC1 Name"),
    "tier": ("Nomenclature_Tier", "Tier"),
    "ambiguity": ("Ambiguity_Interval", "Ambiguity"),
    "repeat_form": ("Repeat_Form", "Repeat Form"),
    "kestrel_name": ("Nomenclature_Kestrel", "Kestrel Name"),
    "advntr_name": ("Nomenclature_adVNTR", "adVNTR Name"),
    "note": ("Nomenclature_Note", "Naming Note"),
}

#: The flag field, under both spellings. ``"Flags"`` is the Kestrel table's heading for
#: ``Nomenclature_Flags``; ``"Flag"`` - singular - is a different column entirely, the
#: configured flagging rule's verdict, and is deliberately not listed here.
_FLAG_FIELDS: tuple[str, ...] = ("Nomenclature_Flags", "Flags")


def numeric_headings(cell_formats: Mapping[str, str]) -> frozenset[str]:
    """Return the headings whose values are numbers.

    Args:
        cell_formats: Heading -> the format constant its values are rendered with.

    Returns:
        frozenset[str]: The headings to align on the figures' axis.
    """
    return frozenset(heading for heading, fmt in cell_formats.items() if fmt in _NUMERIC_FORMATS)


def annotate_table_columns(
    markup: str,
    headings: Sequence[str],
    *,
    numeric: frozenset[str] = frozenset(),
    essential: frozenset[str] = frozenset(),
    caption: str = "",
) -> str:
    """Give every cell of a rendered results table its column's presentation.

    ``DataFrame.to_html`` writes one ``<td>`` per column in order and, with
    ``index=False``, no ``<th>`` outside the header row - so counting cells from each
    ``<tr>`` addresses a column exactly. Four things are added, none of which changes a
    value:

    * ``class="num"`` on a numeric column, so the figures share one right-hand axis.
      The stylesheet has carried the rule since #242 and nothing ever emitted the class,
      so every number in both results tables was centred;
    * ``class="mono-cell"`` on identifiers, sequences and coordinates, which are read
      character by character rather than as words;
    * ``class="col-grow"`` on the one heading that absorbs the table's slack, which is
      the other half of the same omission;
    * ``scope="col"`` and the column's one-line explanation as ``title`` on each
      heading. The explanation is *also* printed under the table by
      :func:`nomenclature_legend` and :data:`COLUMN_HELP`, because a report that is read
      on paper may put nothing behind a pointer.

    A ``caption`` is inserted as the table's first child when one is given. Neither
    results table had one, nor an ``aria-label``: the two tables carrying the findings
    were the only tables in the document with no accessible name at all, while the small
    IGV selector beside them had both. A caption is the name a screen reader's table
    list reads, and it is visible text, so it is also the heading a reader has on paper
    once the ``<h2>`` above it is on the previous sheet.

    Args:
        markup: The table markup, possibly ``""`` for a frame with no rows.
        headings: The display headings, in column order.
        numeric: Which of them hold numbers.
        caption: The table's name. Escaped here; not inserted when empty.

    Returns:
        str: The same table with per-column presentation, or ``""`` unchanged.
    """
    if not markup:
        return markup

    if caption:
        markup = re.sub(
            r"(<table\b[^>]*>)",
            lambda match: f"{match.group(1)}<caption>{escape_html(caption)}</caption>",
            markup,
            count=1,
        )

    column = 0

    def annotate(match: re.Match[str]) -> str:
        nonlocal column
        tag = match.group(1).lower()
        attributes = match.group(2)
        if tag == "tr":
            column = 0
            return f"<tr{_PANDAS_HEADER_ALIGN.sub('', attributes)}>"

        index = column
        column += 1
        if index >= len(headings):
            return match.group(0)

        heading = headings[index]
        classes = []
        if heading in numeric:
            classes.append("num")
        if heading in MONO_COLUMNS:
            classes.append("mono-cell")
        # The one value with genuinely nowhere to break: the unbroken 60 bp motif half. It is the
        # only cell the print block lets break mid-token, because every other value that
        # did so was a number split mid-digit.
        if heading in SEQUENCE_COLUMNS:
            classes.append("sequence-cell")
        if heading in WIDE_COLUMNS:
            classes.append("wide-cell")
        # Marked, not removed. The cell is in the document either way; the class is the
        # hook the reader's own column control folds it behind, and the shared print
        # block brings every one of them back.
        if essential and heading not in essential:
            classes.append("col-optional")
        if tag == "th" and heading in GROW_COLUMNS:
            classes.append("col-grow")

        added = f' class="{" ".join(classes)}"' if classes else ""
        if tag == "th":
            added += ' scope="col"'
            explanation = COLUMN_HELP.get(heading)
            if explanation:
                added += f' title="{escape_html(explanation)}"'
        return f"<{tag}{attributes}{added}>"

    return _TABLE_TAG.sub(annotate, markup)


def folded_record_html(df: pd.DataFrame, essential: frozenset[str], *, noun: str) -> str:
    """Render the columns the table folds away, one labelled block per row, for paper.

    **The printed table cannot hold nineteen columns and never could.** Measured in
    Chromium at A4 portrait, the Kestrel table lays out at 1,442px on a 718px sheet:
    757px - ten of nineteen columns, the motif sequence among them - fell off the
    paper, and nothing on the sheet said so. Restoring every column for print was the
    wrong fix twice over: it is what produced the overflow, and it makes the archived
    PDF, the artefact that outlives the HTML, the *least* complete rendering of the run.

    So the printed table is the essential columns, which fit, and every folded value is
    printed here beneath it - labelled, per row, in the reading order the row has. The
    record is complete on paper for the first time.

    On screen this is not shown, and that is not a hidden value: the same cells are one
    click away in the table itself, under a control that names the full column count.
    Nothing in this block is reachable only here.

    Args:
        df: The display frame, formatted, as rendered into the table.
        essential: The headings the table shows by default; the rest are printed here.
        noun: The algorithm's name, for the block's own heading.

    Returns:
        str: The markup, or "" when the frame is empty or folds nothing.
    """
    folded = [column for column in df.columns if column not in essential]
    if df.empty or not folded:
        return ""

    blocks = []
    for position, (_, row) in enumerate(df.iterrows(), start=1):
        pairs = "".join(
            f"<dt>{escape_html(column)}</dt><dd>{escape_html(row[column])}</dd>"
            for column in folded
            if not (row[column] is None or (isinstance(row[column], float) and pd.isna(row[column])))
        )
        if pairs:
            blocks.append(f'<div class="record-row"><h3>{escape_html(noun)} row {position}</h3><dl>{pairs}</dl></div>')

    if not blocks:
        return ""
    return '<div class="record-appendix">' + "".join(blocks) + "</div>"


def _row_values(df: pd.DataFrame, columns: Sequence[str]) -> list[str]:
    """Return the first present column's non-empty values as text, in row order.

    Args:
        df: A results frame.
        columns: Accepted spellings of one field, most specific first.

    Returns:
        list[str]: The values, empty when no spelling is present or all are blank.
    """
    column = next((name for name in columns if name in df.columns), None)
    if column is None:
        return []
    values = []
    for value in df[column]:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            continue
        text = (value if isinstance(value, str) else str(value)).strip()
        if text and text.lower() != "nan":
            values.append(text)
    return values


def _split_flags(values: Sequence[str]) -> list[str]:
    """Split the semicolon-joined flag fields of several rows into unique tokens."""
    tokens: list[str] = []
    for value in values:
        for token in value.replace(",", ";").split(";"):
            token = token.strip()
            if token and token not in tokens:
                tokens.append(token)
    return tokens


def space_flag_tokens(df: pd.DataFrame) -> pd.DataFrame:
    """Give a semicolon-joined flag cell somewhere to wrap between its tokens.

    ``known-variant;motif-context-diverges;position-ambiguous`` is one unbroken string
    to a line breaker, so a column narrow enough to need a break got one wherever the
    hyphens fell: measured, it wrapped as ``known-variant;motif-context-`` /
    ``diverges;position-`` / ``ambiguous`` - three lines that read as five flags, two of
    which are not words. A space after each separator moves every break to a token
    boundary, so a wrapped cell lists the flags it holds and nothing else.

    This is a display-format change and is recorded in ``docs/user-guide/output-files.md``
    with the others. The separator itself is unchanged, the TSV is untouched, and a
    consumer splitting on ``;`` and stripping is unaffected.

    Args:
        df: A display frame. Not modified.

    Returns:
        pd.DataFrame: A new frame with the flag cells spaced, or the input when it
        carries no flag column.
    """
    column = next((name for name in _FLAG_FIELDS if name in df.columns), None)
    if column is None:
        return df
    spaced = df.copy()
    spaced[column] = spaced[column].map(
        lambda cell: "; ".join(part.strip() for part in cell.split(";")) if isinstance(cell, str) else cell
    )
    return spaced


def tier_reason(tier: str, flags: Sequence[str]) -> str:
    """Say which tier-A condition this call actually missed.

    "At least one tier-A condition is unmet" is a definition, not an answer. A reader
    holding a tier-B call wants to know *which* one, and the flags carry it - but only
    for somebody who already knows that ``motif-context-diverges`` is one of the five
    names :func:`~vntyper.scripts.nomenclature.reconcile` tests before it promotes a
    call. This states the connection the report was leaving to be inferred.

    Only flag-evidenced blockers are named. Two of tier A's conditions leave no flag - a
    single uncorroborated caller, and read support that is unknown rather than low - so
    a call held below tier A by one of those returns "" and the report says the tier
    without inventing a reason for it.

    Args:
        tier: The tier letter the call carries.
        flags: Its flags, already split into single tokens.

    Returns:
        str: The reason clause, or "" for tier A, for tier C - whose whole meaning is
        that no position could be computed - and for a tier-B call whose blocker left no
        flag behind.
    """
    if tier != "B":
        return ""
    reasons = [TIER_A_BLOCKERS[flag] for flag in flags if flag in TIER_A_BLOCKERS]
    if not reasons:
        return ""
    if len(reasons) == 1:
        return f"Held below the corroborated tier because {reasons[0]}."
    return "Held below the corroborated tier because " + "; and ".join(reasons) + "."


def variant_identity(*frames: pd.DataFrame) -> dict[str, Any] | None:
    """Summarise what the run named, for the top of the report.

    The nomenclature fields were readable only by scrolling to column eleven of a
    nineteen-column table and reading across it. They are the answer to the second
    question anybody opens this report with - *what is it called* - so they belong
    where the first answer already is.

    **This never chooses between names.** Where the rows carry more than one distinct
    name, every one of them is returned and the caller renders them all: picking one
    would be the report deciding which variant the reader is shown, which is the whole
    defect #242 exists to remove. The qualifiers are attached only when the rows agree
    on a single name, because a tier or an ambiguity interval belongs to one name and
    attributing it to a list of them would be a claim no row made.

    Args:
        *frames: The results frames, unformatted - the displayed Kestrel frame carries
            markup in some cells and would be summarised from the markup.

    Returns:
        dict[str, Any] | None: The named allele and its qualifiers, or None when no row
        carries a name at all.
    """
    collected: dict[str, list[str]] = {key: [] for key in _IDENTITY_FIELDS}
    flags: list[str] = []
    for frame in frames:
        if frame is None or frame.empty:
            continue
        for key, columns in _IDENTITY_FIELDS.items():
            collected[key].extend(_row_values(frame, columns))
        flags.extend(_row_values(frame, _FLAG_FIELDS))

    names = list(dict.fromkeys(collected["name"]))
    if not names:
        return None

    identity: dict[str, Any] = {
        "name": names[0],
        # Every further name the rows carry. The template lists them; nothing here
        # decides that one of them is the answer.
        "other_names": names[1:],
        "flags": [],
        "tier": "",
        "tier_label": "",
        "tier_meaning": "",
        "ambiguity": "",
        "repeat_form": "",
        "kestrel_name": "",
        "advntr_name": "",
        "note": "",
        "tier_reason": "",
        "callers_agree": False,
    }
    if len(names) > 1:
        return identity

    def single(key: str) -> str:
        """The one value the rows agree on, or "" when they do not carry exactly one."""
        values = list(dict.fromkeys(collected[key]))
        return values[0] if len(values) == 1 else ""

    tier = single("tier")
    label, meaning = NOMENCLATURE_TIERS.get(tier, ("", ""))
    identity["tier"] = tier
    identity["tier_label"] = label
    identity["tier_meaning"] = meaning
    identity["tier_reason"] = tier_reason(tier, _split_flags(flags))
    identity["ambiguity"] = single("ambiguity")
    identity["repeat_form"] = single("repeat_form")
    identity["kestrel_name"] = single("kestrel_name")
    identity["advntr_name"] = single("advntr_name")
    identity["note"] = single("note")
    identity["callers_agree"] = bool(identity["kestrel_name"] and identity["kestrel_name"] == identity["advntr_name"])
    identity["flags"] = [
        {"token": token, "meaning": NOMENCLATURE_FLAG_MEANINGS.get(token, "")} for token in _split_flags(flags)
    ]
    return identity


def nomenclature_legend(*frames: pd.DataFrame) -> list[dict[str, str]]:
    """Explain every coded nomenclature value this report actually prints.

    Only the terms present in these rows are returned. A key listing all three tiers and
    all ten flags under a table that uses two of them is a wall the reader has to filter
    themselves, and the terms they do meet are the ones worth the space.

    Args:
        *frames: The results frames, unformatted.

    Returns:
        list[dict[str, str]]: ``term``/``meaning`` pairs, tiers first, then flags in the
        order the rows carry them.
    """
    tiers: list[str] = []
    flags: list[str] = []
    for frame in frames:
        if frame is None or frame.empty:
            continue
        tiers.extend(_row_values(frame, _IDENTITY_FIELDS["tier"]))
        flags.extend(_row_values(frame, _FLAG_FIELDS))

    entries: list[dict[str, str]] = []
    for tier in dict.fromkeys(tiers):
        label, meaning = NOMENCLATURE_TIERS.get(tier, ("", ""))
        if meaning:
            entries.append({"term": f"Tier {tier}", "label": label, "meaning": meaning})
    for token in _split_flags(flags):
        # `""` rather than `None`, so the name is not rebound to a wider type than the
        # tier branch above gives it - and so an unrecognised token is skipped by the
        # same falsy test either way.
        flag_meaning = NOMENCLATURE_FLAG_MEANINGS.get(token, "")
        if flag_meaning:
            entries.append({"term": token, "label": "", "meaning": flag_meaning})
    return entries


def parse_coverage_stats(data: list[dict[str, Any]]) -> dict[str, Any]:
    """Coerce the first coverage row into the values the report renders.

    Args:
        data: The ``Coverage Calculation`` step's rows. Only the first is used.

    Returns:
        dict[str, Any]: One entry per
        :data:`~vntyper.scripts.coverage_stats.COVERAGE_COLUMNS`. A field is None when
        there is no coverage row, when the column is absent, when it carries
        :data:`~vntyper.scripts.coverage_stats.COVERAGE_NULL_TOKEN`, or when its value
        cannot be coerced.

    Note:
        Coercion is per field. It used to abort the whole row on the first failure, which
        took every later column with it - and ``coverage_qc`` is the last one, so a single
        malformed number discarded the QC verdict (#222).

        A missing column reads as ``None``, never ``0``. Zero is a measurement: it says the
        region was looked at and no reads were found. A summary written before a column
        existed, or by a run that could not compute it, has said nothing at all.
    """
    stats: dict[str, Any] = dict.fromkeys(COVERAGE_COLUMNS)
    if not data or not isinstance(data, list):
        return stats

    row = data[0]
    for name in COVERAGE_COLUMNS:
        raw = row.get(name)
        if name in _BUILD_COMPARABLE_COLUMNS:
            # #222's columns are isolated from the originals in both directions. Absent or
            # not-measured reads as None, never 0 - zero would say the region was looked at
            # and found empty - and an unreadable one is logged and skipped rather than
            # aborting the row, so a column added in 2026 cannot cost a reader the eight
            # statistics and the QC verdict that summaries have always carried.
            if raw is None or raw in (COVERAGE_NULL_TOKEN, ""):
                continue
            try:
                stats[name] = COVERAGE_FIELD_TYPES[name](raw)
            except (TypeError, ValueError) as e:
                logger.error("Error parsing VNTR coverage statistic %s=%r: %s", name, raw, e)
            continue

        # The original columns keep their pre-#222 behaviour exactly: absent coerces from 0,
        # and the first unreadable value stops the row so the fields after it stay None
        # rather than fabricating zeroes. Changing that would change what every existing
        # consumer reads out of a summary it already has.
        try:
            stats[name] = COVERAGE_FIELD_TYPES[name](raw if raw is not None else 0)
        except Exception as e:
            logger.error("Error parsing VNTR coverage statistics: %s", e)
            return stats

    logger.debug("All VNTR coverage statistics extracted successfully: %s", stats)
    return stats


@dataclass(frozen=True)
class FastpMetrics:
    """The four fastp rates the report shows, plus the sequencing setup line.

    Attributes:
        available: Whether fastp output was found at all. The report hides the
            whole section when it was not.
        duplication_rate: Read duplication rate, or None.
        q20_rate: Post-filtering Q20 rate, or None.
        q30_rate: Post-filtering Q30 rate, or None.
        passed_filter_rate: Reads passing the filter over reads before it, or
            None when no reads were seen.
        sequencing: The free-text sequencing setup, e.g. "paired end (150 cycles)".
    """

    available: bool
    duplication_rate: float | None
    q20_rate: float | None
    q30_rate: float | None
    passed_filter_rate: float | None
    sequencing: str


def summarise_fastp(fastp_data: dict[str, Any]) -> FastpMetrics:
    """Reduce parsed fastp JSON to the metrics the report shows.

    Args:
        fastp_data: The parsed ``output.json`` fastp wrote, or ``{}``.

    Returns:
        FastpMetrics: The metrics, with ``available`` False for empty input.
    """
    if not fastp_data:
        return FastpMetrics(False, None, None, None, None, "")

    summary = fastp_data.get("summary", {})
    after_filtering = summary.get("after_filtering", {})
    before_filtering = summary.get("before_filtering", {})
    filtering_result = fastp_data.get("filtering_result", {})

    passed_filter_rate = calculate_passed_filter_rate_from_sources(before_filtering, filtering_result)
    if passed_filter_rate is not None:
        logger.debug("Passed filter rate calculated: %.2f", passed_filter_rate)
    else:
        logger.debug("Total reads before filtering is zero; passed filter rate set to None.")

    return FastpMetrics(
        available=True,
        duplication_rate=fastp_data.get("duplication", {}).get("rate", None),
        q20_rate=after_filtering.get("q20_rate", None),
        q30_rate=after_filtering.get("q30_rate", None),
        passed_filter_rate=passed_filter_rate,
        sequencing=summary.get("sequencing", ""),
    )


def extract_line_after(content: str, marker: str) -> str:
    """Return the text between ``marker`` and the end of its line.

    Both edge cases are load-bearing and both used to be wrong:

    * an absent marker made ``find`` return -1, and ``-1 + len(marker)`` is a
      valid index, so the old code sliced from character 17 and returned
      arbitrary page text where a JavaScript object literal belonged;
    * a marker on the last line with no trailing newline made ``find("\\n", ...)``
      return -1, and slicing to -1 dropped the final character.

    Args:
        content: The document to search.
        marker: The literal that precedes the wanted text.

    Returns:
        str: The stripped remainder of the marker's line, or "" when the marker
        is absent.
    """
    start = content.find(marker)
    if start == -1:
        logger.debug("Marker %r not found.", marker)
        return ""
    start += len(marker)
    end = content.find("\n", start)
    if end == -1:
        end = len(content)
    return content[start:end].strip()


def js_json_literal(fragment: str, fallback: str) -> str:
    """Re-serialise an extracted fragment as a literal that is safe in a ``<script>``.

    ``report_template.html`` interpolates the return value directly into a script
    block as ``const tableJson = {{ table_json|safe }};``. The fragment reaching
    here was lifted verbatim out of the igv-reports page by
    :func:`extract_line_after` and is sample-derived, so it is re-parsed and
    re-emitted rather than trusted: what the template receives is always the
    output of :func:`json.dumps`, never the extracted text.

    A single trailing ``;`` is stripped before parsing. This is defensive against
    a *future* igv-reports version, not a correction of today's: verified against
    the installed igv-reports 1.16.0, ``templates/variant_template.html:155-156``
    is ``const tableJson = "@TABLE_JSON@"`` with no terminator, and
    ``report.py:178-183`` substitutes the placeholder including its quotes, so the
    fragment reaching here never carries a trailing ``;`` today.

    ``json.dumps`` is called with ``ensure_ascii=True`` (the stdlib default, kept
    explicit here because it is load-bearing): it escapes every non-ASCII
    codepoint, which includes U+2028 and U+2029 -- line terminators to a
    JavaScript parser that are legal inside a JSON string. That leaves ``<`` as
    the one remaining script-context hazard, and **every** literal ``<`` is
    escaped rather than only the ``</`` of a closing tag. Escaping ``</`` alone
    stops a direct ``</script>`` but not the HTML5 tokenizer's *double-escaped*
    script state: ``<!--`` followed by ``<script`` inside a script element -- a
    sequence containing no ``</`` at all -- puts the parser into a state where the
    real ``</script>`` no longer terminates the element, and the remainder of the
    document is consumed as script text. Escaping ``<`` subsumes both routes.

    ``\\u003c`` is a JSON string escape, so what the browser parses back is the
    original ``<`` and the data reaching the page is unchanged. Only string values
    can contain a ``<`` -- JSON's structural characters do not -- so the
    replacement cannot touch the literal's syntax.

    Keys are sorted and separators are minimised so that two runs over the same
    IGV page emit byte-identical script.

    Args:
        fragment: The extracted literal, possibly empty, possibly not JSON.
        fallback: A valid literal to use when the fragment cannot be parsed. An
            empty fragment would otherwise produce ``const tableJson = ;`` -- a
            syntax error that takes the whole script block down, and with it the
            variant table, the flag toggles and the coverage switch.

    Returns:
        str: A JSON literal that parses as JavaScript and cannot escape the
        script block.
    """
    candidate = fragment.strip().removesuffix(";").strip()
    if not candidate:
        return fallback
    try:
        value = json.loads(candidate)
    except ValueError as e:
        logger.warning(f"IGV fragment could not be parsed as JSON and was discarded: {e}")
        return fallback

    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return encoded.replace("<", "\\u003c")


def extract_igv_fragments(content: str) -> tuple[str, str, str]:
    """Split the three fragments the report needs out of an IGV report page.

    Args:
        content: The full text of the page ``create_report`` wrote.

    Returns:
        tuple[str, str, str]: The ``#container`` markup, the ``tableJson``
        literal and the ``sessionDictionary`` literal. All three are empty when
        the container cannot be located.
    """
    igv_start = content.find(IGV_CONTAINER_MARKER)
    igv_end = content.find(IGV_BODY_END_MARKER)

    if igv_start == -1 or igv_end == -1:
        logger.error("Failed to extract IGV content from report.")
        return "", "", ""

    igv_content = content[igv_start:igv_end].strip()
    table_json = extract_line_after(content, IGV_TABLE_JSON_MARKER)
    session_dictionary = extract_line_after(content, IGV_SESSION_DICTIONARY_MARKER)

    logger.info("Successfully extracted IGV content, tableJson, and sessionDictionary.")
    return igv_content, table_json, session_dictionary
