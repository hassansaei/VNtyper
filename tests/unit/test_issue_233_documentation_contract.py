"""Documentation contracts for #233 and the bounded published #295 exceptions.

Three tests were removed when planning documents left the repository. They asserted
that historical design documents under `docs/superpowers/specs/` and `docs/plans/`
led with a dated supersession note, so that a reader meeting an old design was told
its routing policy had been replaced. Those historical documents now live in the
untracked `.planning/` workspace. The #295 program separately publishes its reviewed
design and implementation plans under `docs/superpowers/`; they are built by MkDocs and
do not reopen a general planning-docs exception.

What survives here is the half that was always about published output: the changelog
and the golden-cohort page state the current routing policy, and AGENTS.md states the
contracts an agent has to honour. Those are read by people who never see `.planning/`.

The changelog assertions below were previously entangled with the two planning
documents in one loop. They are kept, applied to the changelog alone.
"""

import hashlib
import re
from pathlib import Path

import pytest

from vntyper.scripts import nomenclature_bam
from vntyper.scripts.nomenclature_presentation import NOMENCLATURE_FLAG_MEANINGS

pytestmark = pytest.mark.unit

ROOT = Path(__file__).resolve().parents[2]


def _read(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def _section(page: str, start: str, end: str) -> str:
    return page[page.index(start) : page.index(end)]


def test_agents_records_the_exact_issue_295_published_planning_exceptions() -> None:
    """The exceptional public pages stay exact while ordinary planning stays untracked."""
    page = _read("AGENTS.md")
    normalized = " ".join(page.split())
    assert "docs/superpowers/specs/2026-08-31-kestrel-bam-evidence-semantics-design.md" in page
    assert "docs/superpowers/plans/2026-08-31-kestrel-bam-evidence-semantics.md" in page
    assert "docs/superpowers/specs/2026-08-31-issue-295-completion-program-design.md" in page
    assert "docs/superpowers/plans/2026-08-31-issue-295-completion-program.md" in page
    assert "No other planning artifact under `docs/` is allowed" in normalized
    assert "untracked `.planning/` workspace" in normalized


def test_agents_inventory_includes_both_phase_1_focused_modules() -> None:
    """Catch stale sixteen-module wording after the evidence split has landed."""
    page = _read("AGENTS.md")
    layout = _section(page, "## Layout", "## Code style")
    normalized = " ".join(layout.split())

    assert "These eighteen focused modules" in normalized
    assert "`nomenclature_evidence.py` — source-specific evidence-unit flags and the" in normalized
    assert "BAM thin-support configuration-key resolver" in normalized
    assert "`nomenclature_bam.py` separately owns XD parsing, resolved haplotype-record voting, and" in normalized
    assert "`BamConsensus`" in normalized
    evidence_line = normalized[
        normalized.index("`nomenclature_evidence.py`") : normalized.index("`nomenclature_presentation.py`")
    ]
    assert "XD parsing" not in evidence_line
    assert "`nomenclature_presentation.py` — source-specific report flag meanings" in normalized
    assert "reserves two further focused destinations" not in normalized


def test_changelog_states_the_current_lossless_selection_and_invalid_parity_policy() -> None:
    """The published record of #233 must name the policy, not merely that it changed.

    Previously this also scanned two planning documents in the same loop; they are no
    longer tracked. The changelog is the shipped statement and keeps every assertion.
    """
    changelog = _read("docs/about/changelog.md")
    section = " ".join(changelog[changelog.index("## 2.0.12") : changelog.index("## 2.0.11")].split())

    assert "singleton" in section
    assert "`other`" in section
    assert "every non-empty" in section
    assert "exactly once" in section
    assert "one Kestrel" in section
    assert "Unequal or one-sided" in section
    assert "invalid" in section
    assert "not concatenated or recompressed" in section


def test_changelog_preserves_released_history_and_records_the_fix_in_2_0_12() -> None:
    page = _read("docs/about/changelog.md")
    assert page.index("## Unreleased") < page.index("## 2.0.12") < page.index("## 2.0.11")
    unreleased = page[page.index("## Unreleased") : page.index("## 2.0.12")]
    release = page[page.index("## 2.0.12") : page.index("## 2.0.11")]
    # The section may legitimately be empty or carry entries; what #233 pins is that
    # released history below it is preserved, not that nothing is ever in flight.
    assert "No unreleased changes." in unreleased or "###" in unreleased
    assert "Valid mixed alignment conversions are now routed losslessly" in release
    assert "Superseded on 2026-08-11 by #233" in page
    assert "## 2.0.11 (Current)" not in page


def test_phase_1_changelog_entry_does_not_rewrite_released_history() -> None:
    """Catch edits to released history while allowing the Unreleased section to evolve.

    The Phase 1 entry shipped in 2.0.27, so everything from that heading down is released
    history and is pinned byte-for-byte. A release bump legitimately moves this pin once:
    it demotes the previous ``(Current)`` marker and adds the new release above it.
    """
    page = _read("docs/about/changelog.md")
    released = page[page.index("## 2.0.27") :]
    assert hashlib.sha256(released.encode()).hexdigest() == (
        "b67335664e7041f2cba903e3c5e615f8daef3cfff504461d14da4baf8108cd72"
    )

    phase_1_release = _section(page, "## 2.0.27", "## 2.0.26")
    normalized = " ".join(phase_1_release.split())
    assert "Phase 1" in normalized
    assert "Refs #295" in normalized
    assert "resolved haplotype records" in normalized
    assert "minimum k-mer depth" in normalized
    assert "does not weight votes or alter names or tiers" in normalized
    assert "bam_thin_haplotype_record_support" in phase_1_release
    assert "bam_thin_support" in phase_1_release
    assert "supporting_haplotype_records" in phase_1_release
    assert "fetched_haplotype_records" in phase_1_release
    assert "distinct_edit_count" in phase_1_release
    assert "read-only `support`, `total`, and `n_distinct` compatibility properties" in normalized
    assert "replace the pre-Phase-1 `low-read-support` token" in normalized
    assert "BAM-specific `thin-haplotype-record-support` and `low-haplotype-record-support`" in normalized
    assert all(f"#{issue}" in phase_1_release for issue in (270, 267, 269))


def test_nomenclature_page_documents_every_authoritative_flag_and_source_unit() -> None:
    """Catch a deleted flag row or a source-specific evidence unit relabelled as reads."""
    page = _read("docs/pipeline/nomenclature.md")
    flags = _section(page, "## Flags", "## The two companion fields")
    documented = set(re.findall(r"^\| `([^`]+)` \|", flags, flags=re.MULTILINE))

    assert documented == set(NOMENCLATURE_FLAG_MEANINGS)
    assert "all 14" in flags
    assert "Kestrel VCF" in flags and "k-mer-path depth" in flags
    assert "Kestrel `output.bam`" in flags and "haplotype-record support" in flags
    assert "adVNTR" in flags and "sequencing reads" in flags
    assert "`min_support_for_high_confidence`" in flags
    assert "`bam_thin_haplotype_record_support`" in flags
    assert "`bam_thin_support`" in flags


def test_nomenclature_page_pins_the_measured_bam_consultation_counts() -> None:
    """Catch restoration of the stale approximate benchmark claim."""
    page = _read("docs/pipeline/nomenclature.md")
    section = _section(page, "## Where the resolved haplotype records are consulted", "## When two sources")
    normalized = " ".join(section.split())

    assert "83/200 samples are eligible" in normalized
    assert "15 have no Kestrel result row" in normalized
    assert "68/200 (34%) produce one BAM row fetch" in normalized
    assert "about a fifth" not in section
    assert "minimum k-mer depth" in normalized
    assert "does not weight votes or alter names or tiers" in normalized


def test_nomenclature_page_states_the_typed_xd_contract() -> None:
    """Catch collapsing retained zero/extreme integers into unavailable XD evidence."""
    page = " ".join(_read("docs/pipeline/nomenclature.md").split()).lower()

    maximum = f"{nomenclature_bam._MAXIMUM_MINIMUM_KMER_DEPTH:,}"
    assert f"integer values from 1 through {maximum} are retained exactly" in page
    assert "zero is retained as zero" in page
    assert "missing or malformed values are unavailable" in page
    assert f"negative values and unsigned integers above {maximum} are unavailable" in page
    assert "every resolved haplotype record still contributes one unweighted vote" in page


def test_report_and_output_pages_explain_the_kestrel_bam_artifact() -> None:
    """Catch omission of the concise HTML and IGV semantics from a public surface."""
    pages = {
        path: " ".join(_read(path).split())
        for path in (
            "docs/pipeline/reports.md",
            "docs/user-guide/output-files.md",
            "docs/cli/report.md",
            "docs/getting-started/quickstart.md",
            "vntyper/scripts/README.md",
        )
    }
    for path, page in pages.items():
        assert "resolved haplotype records" in page, path
        assert "output.bam" in page, path

    for path in ("docs/pipeline/reports.md", "docs/user-guide/output-files.md", "docs/cli/report.md"):
        page = pages[path]
        assert "not sequencing reads" in page, path
        assert "minimum k-mer depth" in page, path
        assert "does not weight votes or alter names or tiers" in page, path


def test_kestrel_depth_pages_use_kmer_units_without_changing_the_rules() -> None:
    """Catch either Kestrel Sample depth being called reads or a threshold rewrite."""
    output_page = _read("docs/user-guide/output-files.md")
    scoring = " ".join(_read("docs/pipeline/scoring-and-confidence.md").split())
    kestrel = _read("docs/pipeline/kestrel.md")

    assert "k-mer-path depth supporting the alternate allele" in output_page
    assert "total k-mer depth in the variant active region" in output_page
    assert "Alternate depths come from Kestrel's `Sample` field as k-mer-path depths" in scoring
    assert "total active-region k-mer depth of 150" in scoring
    assert "total active-region k-mer depth of 5000" in scoring
    assert "same alternate k-mer-path depth" in kestrel
    assert "0.00469" in scoring and "0.00515" in scoring


def test_genuine_input_and_advntr_read_language_remains_public() -> None:
    """Prevent an ontology correction from globally deleting legitimate read counts."""
    reports = _read("docs/pipeline/reports.md")
    quickstart = _read("docs/getting-started/quickstart.md")
    scripts_readme = _read("vntyper/scripts/README.md")

    assert "supporting read count" in reports
    assert "Extracted FASTQ reads" in quickstart
    assert "BAM/CRAM" in scripts_readme and "FASTQ" in scripts_readme


def test_golden_page_records_the_issue_233_run_and_retained_evidence() -> None:
    page = _read("docs/development/golden-cohort-gate.md")
    normalized = " ".join(page.split())
    assert "Current routing policy (2026-08-11, #233)" in page
    assert "32 base cases plus 10 repeat/derived cases, 42 total" in normalized
    assert "comparison.json" in page
    assert "| 7 | `19c8acd` | `4678851` |" in page
    assert "5b8dc9199cd19fc1142e0a6ba7bd2740d4c0a97b0cdd9e5f8f4b08e51330e88e" in page
    assert "6808936b98be8b8d79decd17c76f89f5f4519a6e1fa9acc3f96c0c9eb6d14cbd" in page
    assert "d3b17029f55c4a610d708764bf4b9c5298f2caad3f0f3114ce532b79b43b41a3" in page
    assert "8a0c1a0460934cecf9db19b659c7f219f964bf685bf5f818bf12b2b3b69bac10" in page
    assert "6f09f9350d152ab1b69aa07cf2096aad895a01cca08f23c297608cf772029dd0" in page
    assert "67/67 pipeline/probe outcomes" in normalized
    assert "4/4 cohort outcomes" in normalized
    assert "This page also records the final gated candidate SHA" in normalized
    assert "all three adVNTR successes" in normalized
    assert "25 cases ran successfully on both sides" in normalized
    assert "preliminary issue #233 run at `49b0cc6`" in normalized
    assert "No issue #233 run result is recorded here" not in page


def test_agents_scopes_the_host_reference_prerequisite_to_integration() -> None:
    page = _read("AGENTS.md")
    testing = page[page.index("## Testing") : page.index("## Git and PRs")]
    normalized = " ".join(testing.split())
    assert "Integration tests require" in normalized
    assert "pinned host `reference/alignment/chr1.hg19.fa`" in normalized
    assert "Docker tests use the pinned reference bundled in the image" in normalized


def test_agents_records_the_fail_closed_append_only_success_contract() -> None:
    page = _read("AGENTS.md")
    trap = page[page.index("17. **Real-data successes are append-only compatibility contracts.") :]
    assert "tests/compatibility/real_success_baseline.json" in trap
    assert "(suite, test_name)" in trap
    assert "baseline-to-live" in trap
    assert "live-to-baseline" in trap
    assert "origin/${{ github.base_ref }}" in trap
    assert "${{ github.event.before }}" in trap
    assert "missing, empty, all-zero, shallow, or unreachable bases fail closed" in trap
    assert "there is no exception or free-text waiver" in trap
