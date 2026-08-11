from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

ROOT = Path(__file__).resolve().parents[2]
INTRODUCER = "2b4597d8b57f8986f008baad6c2505359cebea76"
ORACLE_INVERSION = "52c4146596fef2d1e2402991fbab062ba8021889"
DATE = "2026-08-11"


def _read(relative: str) -> str:
    return (ROOT / relative).read_text(encoding="utf-8")


def test_historical_designs_lead_with_the_dated_issue_233_supersession() -> None:
    for relative, historical_phrase in (
        (
            "docs/superpowers/specs/2026-08-10-milestone-6-harness-matrix-design.md",
            "intentionally failing mixed",
        ),
        (
            "docs/plans/2026-08-08-milestone-4-cram-input-robustness-spec.md",
            "is **rejected**",
        ),
    ):
        page = _read(relative)
        note = page.index("Superseded routing policy")
        assert note < page.index(historical_phrase)
        supersession = page[note : note + 1_500]
        assert DATE in supersession
        assert f"`{INTRODUCER}` introduced" in supersession
        assert f"`{ORACLE_INVERSION}` later" in supersession
        assert "nine BAM and one adVNTR" in supersession


def test_current_policy_names_lossless_selection_and_invalid_parity() -> None:
    milestone_6 = _read("docs/superpowers/specs/2026-08-10-milestone-6-harness-matrix-design.md")
    milestone_4 = _read("docs/plans/2026-08-08-milestone-4-cram-input-robustness-spec.md")
    changelog = _read("docs/about/changelog.md")
    current_sections = (
        milestone_6[milestone_6.index("Superseded routing policy") : milestone_6.index("## Purpose")],
        milestone_4[milestone_4.index("Superseded routing policy") : milestone_4.index("---")],
        changelog[changelog.index("## Unreleased") : changelog.index("## 2.0.11 (Current)")],
    )
    normalized_sections = tuple(" ".join(section.split()) for section in current_sections)
    for section in normalized_sections:
        assert "singleton" in section
        assert "`other`" in section
        assert "every non-empty" in section
        assert "exactly once" in section
        assert "one Kestrel" in section
        assert "Unequal or one-sided" in section
        assert "invalid" in section
    assert "R1/R2/other/single order" in normalized_sections[0]
    assert "R1/R2/other/single order" in normalized_sections[1]
    assert "not concatenated or recompressed" in normalized_sections[1]
    assert "not concatenated or recompressed" in normalized_sections[2]


def test_docker_quick_is_named_as_a_real_b178_success_sentinel() -> None:
    page = _read("docs/superpowers/specs/2026-08-10-milestone-6-harness-matrix-design.md")
    assert "example_b178_hg19_subset_default" in page
    assert "Docker quick" in page
    assert "required success sentinel" in page


def test_changelog_preserves_released_history_and_records_the_fix_unreleased() -> None:
    page = _read("docs/about/changelog.md")
    assert page.index("## Unreleased") < page.index("## 2.0.11 (Current)")
    unreleased = page[page.index("## Unreleased") : page.index("## 2.0.11 (Current)")]
    assert "Valid mixed alignment conversions are now routed losslessly" in unreleased
    assert "Superseded on 2026-08-11 by #233" in page
    assert "## 2.0.11 (Current)" in page


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
