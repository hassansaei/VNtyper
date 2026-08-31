"""Contain Kestrel BAM evidence terminology to its truthful ontology."""

import ast
import inspect
import re
import textwrap
from collections.abc import Callable
from pathlib import Path

import pytest

from vntyper.scripts import (
    nomenclature,
    nomenclature_annotate,
    nomenclature_bam,
    nomenclature_evidence,
    nomenclature_presentation,
    report_formatting,
)

pytestmark = pytest.mark.unit

ROOT = Path(__file__).resolve().parents[2]


KESTREL_BAM_FUNCTIONS: tuple[Callable[..., object], ...] = (
    nomenclature_bam.BamRescuer.rescue,
    nomenclature_bam.from_bam,
    nomenclature_bam.refine,
    nomenclature_annotate.annotate_kestrel_frame,
    nomenclature_annotate._open_rescuer,
    nomenclature_annotate._row_haplotype_call,
    nomenclature_annotate._row_verdicts,
    nomenclature_annotate._haplotype_calls,
)

READ_NAMED_IDENTIFIER = re.compile(r"(?:^|_)(?:read|reads)(?:_|$)")
BANNED_BAM_PROSE = (
    "row's reads",
    "reads per source",
    "what the reads say",
    "the reads, as a third source",
)


def _identifiers(source: str) -> set[str]:
    tree = ast.parse(textwrap.dedent(source))
    return {node.id for node in ast.walk(tree) if isinstance(node, ast.Name)} | {
        node.arg for node in ast.walk(tree) if isinstance(node, ast.arg)
    }


def _read_named_identifiers(source: str) -> set[str]:
    return {identifier for identifier in _identifiers(source) if READ_NAMED_IDENTIFIER.search(identifier)}


def _forbidden_bam_prose(source: str) -> list[str]:
    lowered = source.lower()
    return [phrase for phrase in BANNED_BAM_PROSE if phrase in lowered]


def _forbidden_published_kestrel_bam_prose(source: str) -> list[str]:
    lowered = " ".join(source.lower().split())
    return [
        phrase
        for phrase in (
            "where the reads are consulted",
            "the reads may supply an allele",
            "thin read consensus",
            "one to three reads",
            "alignment splits reads",
            "sample` field as read counts",
            "read depth supporting the alternate allele",
            "total read depth in the variant active region",
        )
        if phrase in lowered
    ]


@pytest.mark.parametrize("function", KESTREL_BAM_FUNCTIONS, ids=lambda function: function.__name__)
def test_kestrel_bam_functions_do_not_name_haplotype_records_as_reads(function: Callable[..., object]) -> None:
    """A read-named BAM local would fail even if output stayed unchanged."""
    assert not _read_named_identifiers(inspect.getsource(function))


def test_identifier_guard_catches_underscore_delimited_read_names_but_not_file_reading() -> None:
    """Catch broad BAM-local regressions without rejecting pandas/path reading APIs."""
    mutation = """
def rescue(path):
    fetched_reads = path.read_text()
    bam_reads = fetched_reads
    read_count = len(bam_reads)
    frame = pd.read_csv(path)
    return read_count, frame
"""

    assert _read_named_identifiers(mutation) == {"bam_reads", "fetched_reads", "read_count"}


@pytest.mark.parametrize("function", KESTREL_BAM_FUNCTIONS, ids=lambda function: function.__name__)
def test_kestrel_bam_function_prose_does_not_call_haplotype_records_reads(
    function: Callable[..., object],
) -> None:
    """Inspect each full function so old BAM-specific prose cannot return."""
    assert not _forbidden_bam_prose(inspect.getsource(function))


@pytest.mark.parametrize(
    "phrase",
    ("row's reads", "reads per source", "What the reads say", "The reads, as a third source"),
)
def test_bam_prose_guard_catches_each_deliberate_regression(phrase: str) -> None:
    assert _forbidden_bam_prose(f'def rescue():\n    """{phrase}."""\n') == [phrase.lower()]


def test_old_read_named_bam_helpers_cannot_return() -> None:
    """The former private names encoded the wrong evidence ontology in interfaces."""
    source = inspect.getsource(nomenclature_annotate)

    assert "_row_read_call" not in source
    assert "_read_calls" not in source


def test_reconciliation_prose_does_not_turn_haplotype_support_into_reads() -> None:
    """Comments and docstrings must not reintroduce the old mixed-unit explanation."""
    source = "\n".join(
        (
            inspect.getsource(nomenclature._is_corroborated),
            inspect.getsource(nomenclature.reconcile),
        )
    ).lower()
    forbidden_phrases = (
        "caller read twice",
        "row's reads",
        "reads per source",
        "adVNTR and the reads",
        "1-read agreement",
        "one source's reads",
    )

    assert not [phrase for phrase in forbidden_phrases if phrase.lower() in source]


def test_bam_tie_log_names_haplotype_records() -> None:
    """Restoring read-count wording in the abstention log must fail this source contract."""
    source = inspect.getsource(nomenclature_bam.BamRescuer.rescue)

    assert re.search(r"alleles tied on %d haplotype records", source)
    assert not re.search(r"alleles tied on %d reads", source)


def test_evidence_flags_are_contained_in_the_public_nomenclature_vocabulary() -> None:
    """Every evidence token has one authoritative public spelling and report meaning."""
    evidence_flags = {name: value for name, value in vars(nomenclature_evidence).items() if name.startswith("FLAG_")}
    public_flags = {name: value for name, value in vars(nomenclature).items() if name.startswith("FLAG_")}

    assert evidence_flags == {
        "FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT": "thin-haplotype-record-support",
        "FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT": "low-haplotype-record-support",
        "FLAG_LOW_KMER_PATH_SUPPORT": "low-kmer-path-support",
        "FLAG_LOW_READ_SUPPORT": "low-read-support",
        "FLAG_LOW_EVIDENCE_SUPPORT": "low-evidence-support",
    }
    assert evidence_flags.items() <= public_flags.items()
    assert set(public_flags.values()) == set(nomenclature_presentation.NOMENCLATURE_FLAG_MEANINGS)


def test_kestrel_bam_presentation_never_labels_record_support_as_reads() -> None:
    """BAM-specific meanings may contrast records with reads, but cannot equate them."""
    bam_specific_text = (
        nomenclature_presentation.NOMENCLATURE_FLAG_MEANINGS[nomenclature.FLAG_ALLELE_UNREPRESENTABLE],
        nomenclature_presentation.NOMENCLATURE_FLAG_MEANINGS[nomenclature.FLAG_THIN_HAPLOTYPE_RECORD_SUPPORT],
        nomenclature_presentation.NOMENCLATURE_FLAG_MEANINGS[nomenclature.FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT],
        nomenclature_presentation.TIER_A_BLOCKERS[nomenclature.FLAG_LOW_HAPLOTYPE_RECORD_SUPPORT],
    )

    assert all(not re.search(r"\breads?\b", text, flags=re.IGNORECASE) for text in bam_specific_text)
    assert "not sequencing reads" in nomenclature_presentation.KESTREL_BAM_SEMANTICS
    assert "resolved haplotype records" in nomenclature_presentation.KESTREL_BAM_SEMANTICS


def test_presentation_reexports_stay_compatible_and_read_wording_stays_source_specific() -> None:
    """The report keeps its public tables and genuine adVNTR read-count vocabulary."""
    assert report_formatting.NOMENCLATURE_FLAG_MEANINGS is nomenclature_presentation.NOMENCLATURE_FLAG_MEANINGS
    assert report_formatting.TIER_A_BLOCKERS is nomenclature_presentation.TIER_A_BLOCKERS
    assert report_formatting.COLUMN_HELP is nomenclature_presentation.COLUMN_HELP
    assert "Reads adVNTR counted" in report_formatting.COLUMN_HELP["Supporting Reads"]
    assert "Reads adVNTR counted" in report_formatting.COLUMN_HELP["NumberOfSupportingReads"]


def test_published_kestrel_bam_pages_do_not_call_resolved_records_reads() -> None:
    """Scan only Kestrel evidence pages so genuine input/adVNTR reads remain allowed."""
    paths = (
        "docs/pipeline/nomenclature.md",
        "docs/pipeline/reports.md",
        "docs/user-guide/output-files.md",
        "docs/cli/report.md",
        "docs/getting-started/quickstart.md",
        "docs/pipeline/kestrel.md",
        "docs/pipeline/scoring-and-confidence.md",
        "vntyper/scripts/README.md",
    )

    failures = {
        path: found
        for path in paths
        if (found := _forbidden_published_kestrel_bam_prose((ROOT / path).read_text(encoding="utf-8")))
    }
    assert not failures


@pytest.mark.parametrize(
    "phrase",
    (
        "Where the reads are consulted",
        "the reads may supply an allele",
        "thin read consensus",
        "one to three reads",
        "the alignment splits reads",
        "Sample` field as read counts",
        "Read depth supporting the alternate allele",
        "Total read depth in the variant active region",
    ),
)
def test_published_prose_guard_catches_each_deliberate_regression(phrase: str) -> None:
    assert _forbidden_published_kestrel_bam_prose(phrase)
