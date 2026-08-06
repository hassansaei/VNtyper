# #181–#197 Follow-ups Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close all seventeen #179 follow-up issues, enable branch coverage, and kill the measured mutation and coverage gaps — without moving a reported genotype except where @hassansaei has decided one should move.

**Architecture:** Seven phases (0–6) ordered so that every genotype-affecting change lands **last**, after the instruments that would detect a regression are in place and after the refactors that would otherwise invalidate the attestation. Phases 1 and 2 fan out to six parallel agents each over disjoint file sets. Phase 5 is strictly serial: one genotype-affecting commit at a time, each with its own golden-cohort run naming its own candidate SHA, because a joint run can neither attribute a failure nor detect two changes cancelling each other out.

**Tech Stack:** Python 3.10-floor (`pyproject.toml:10`; AGENTS.md's "3.9" is stale), pandas, pytest, ruff, mypy, coverage.py, Jinja2. Conda env `vntyper`.

**Spec:** `specs/2026-08-06-issue-181-197-followups-design.md`

---

## Global Constraints

Every task's requirements implicitly include this section.

- **Run every command as** `conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" <cmd>`. A stale `~/.local/bin/pytest` shadows the env's and silently breaks `--cov`.
- **Run pytest from the repo root.** `tests/parametrization.py` opens `tests/test_data_config.json` by relative path at collection time.
- **Failing test first, observed failing.** A test never seen to fail is not written. Paste the failure output into the commit or the task record.
- **`pytestmark = pytest.mark.unit` immediately after the imports of every new test file.** CI runs `pytest -m unit`; an unmarked file silently never runs.
- **Never lower `fail_under`.** It is a ratchet, currently `70` at `pyproject.toml:341`. Raise it only to the number `scripts/coverage_gate.py` prints — never the rounded `TOTAL` column.
- **No value in `vntyper/scripts/kestrel_config.json` changes** except `duplicate_flagging.sort_by` in Task 7, which quotes @hassansaei's decision on #197.
- **No confidence, motif or frameshift decision logic changes** except Tasks 9, 10, 17, 18, 19 and 6a — each quoting a recorded decision, each its own revertible commit, **each with its own golden-cohort run naming its own candidate SHA**.
- **Every characterisation-to-specification upgrade names the individual test it upgrades.** Never a section banner, never a module docstring. `test_scoring.py`'s banner also covers sample-field parsing, `test_advntr_frameshift_filter.py`'s also covers parsing columns and import-time config, and `test_confidence_boundaries.py`'s covers the whole 54-cell matrix **including the cell Task 17 changes**. Upgrading a banner ratifies behaviour nobody decided.
- **Line length 120, double quotes, ruff `target-version = "py310"`. Code must run on Python 3.10** — `pyproject.toml:10` is `requires-python = ">=3.10"` and `pyproject.toml:157` is `py310`. **AGENTS.md is wrong on this**: it says "Code must run on Python 3.9 (CI matrix: 3.9–3.12)" and `target-version = "py39"`. The floor moved and the CI matrix is now 3.10–3.13 (verified on #198's checks). Task 24 corrects AGENTS.md; until then, trust `pyproject.toml`. Do not reformat to 88 columns.
- **Logging:** `logger = logging.getLogger(__name__)` after imports; f-strings in log calls are the house style. No `basicConfig`, no per-module handlers.
- **Errors:** no custom exception classes. `logger.error(msg)` then `raise ValueError(msg)` with the same message.
- **Do not touch base-image content-hash inputs:** `conda/**`, `docker/Dockerfile.base`, `docker/requirements-web.txt`, `vntyper/scripts/install_references.py`, `install_references_config.json`, `vntyper/dependencies/advntr/**`, `reference/**`, `.dockerignore`. Editing any of them — even a comment — turns the Docker build red until a new base publishes. **Task 23 is the exception and is explicitly deferred.**
- **Never add a page under `docs/` without registering it in `mkdocs.yml` `nav:`.** `mkdocs build --strict` fails on a dangling entry or broken internal link.
- **`tests/unit/test_marker_hygiene.py` fails under any filtered run, `-k` included.** Scope with a path (`pytest tests/unit/test_x.py`), never with `-k`, and run the unfiltered tier before committing.

## Parallel execution model

The git index is shared across the worktree and subagents cannot leave it. Two agents staging concurrently produce commits containing each other's files.

**Rule: agents author, the orchestrator commits.**

- Each agent owns a disjoint `{source file, test file, doc file}` set, listed in its task.
- **No agent runs `git add`, `git commit`, `git stash`, or `git checkout`.** The agent's deliverable is files on disk plus the pasted failing-then-passing test output.
- The orchestrator commits each lane serially, in task order, after the wave completes.
- **Orchestrator-only, never inside an agent:** `make format` (ruff `--fix` rewrites the whole tree), `make mutation` (rewrites `vntyper/scripts/*.py` in place — nothing may build, package or install while it runs), and any full `pytest -m unit tests/unit`.
- Parallel coverage runs must each set a distinct `COVERAGE_FILE`, or they corrupt each other's data.

| Phase | Lanes | Tasks | Gate |
| --- | ---: | --- | --- |
| 0 | 1 | 0 | CI on #198 |
| 1 | **6** | 1, 2, 3, 4, 5, 7, 6b (#194 only) | `make check-all` |
| 2 | **6** | 8, 11, 12, 13, 14, 15 | `make check-all` |
| 3 | 1 | 16 | `make ci-local` |
| 4 | 1 | 22 (`install_references.py` stays parked as Task 23) | `make check-all` |
| 5 | **1, strictly serial** | 9, 10, 17, 18, 19, 6a (#188) — **one golden-cohort run each** | golden cohort per commit |
| 6 | 1 | 21, 24 | golden cohort at the PR tip; `make ci-local` |

**Four sequencing corrections applied after the adversarial review. Read these before executing — they changed which phase several tasks live in.**

1. **The gate attests one commit and nothing after it.** The first draft ran the gate in phase 4 and then added file splits in phase 5, which would have invalidated it — the same class of error a prior review caught on this project. The splits move *before* the gated work, and the gated work is now last.
2. **One run per genotype-affecting commit, not one cumulative run.** A joint run cannot attribute a failure, and two changes can produce offsetting deltas that read as clean. Five runs at ~25 min is ~2 hours; that is the price of attribution.
3. **Tasks 9 and 10 (the `motif_processing` refactor) move from phase 2 to the gated phase.** They extract the code that computes `motif_filter_pass` (`motif_processing.py:430`), which `filter_final_dataframe` requires (`kestrel_genotyping.py:806`). "Intended to be behaviour-preserving" is not evidence; 58 real samples are.
4. **Task 6 splits.** #194 (mypy) stays in phase 1 as **Task 6b**. #188 (CRAM) moves to the gated phase as **Task 6a**, because it is neither one line nor zero-risk — see the task.
5. **Task 20 (`.get()` calibration defaults) is removed from this plan.** The review is right that #185 decided gate columns, not calibration keys, and this plan's own constraint forbids changing confidence code without an explicit decision. It becomes an issue to file in Task 24, not a commit.

---

## File Structure

**Created**

| Path | Responsibility |
| --- | --- |
| `tests/unit/test_cross_match.py` | The five `cross_match.py` functions, including the `eval`'d rule string |
| `tests/unit/test_frameshift_convention_parity.py` | Asserts Kestrel and adVNTR encode one frameshift convention |
| `tests/unit/test_motif_config_guard.py` | The #186 armed-combination guard |
| `tests/unit/test_motif_preprocessing.py` | The VCF column contract shared by both preprocessors |
| `tests/unit/test_motif_decisions.py` | The extracted pure motif decision layer |
| `tests/unit/test_file_processing.py` | `file_processing.py` |
| `tests/unit/test_extract_unmapped.py` | `extract_unmapped_from_offset.py` |
| `vntyper/scripts/motif_decisions.py` | Pure motif filtering decisions, extracted from `motif_correction_and_annotation` |

**Modified**

| Path | Change |
| --- | --- |
| `vntyper/scripts/scoring.py` | Docstring: the +1-frame biology (#181) |
| `vntyper/scripts/confidence_assignment.py` | Mid-band ordering (#184) **only**. The `.get()` default removal was **cut** — Task 20 is filed as an issue, not implemented |
| `vntyper/scripts/kestrel_genotyping.py` | `filter_final_dataframe` fails closed (#185) |
| `vntyper/scripts/motif_processing.py` | Config guard (#186); twin preprocessors merged; decisions extracted |
| `vntyper/scripts/flagging.py` | `kind="stable"` on the duplicate sort (#197) |
| `vntyper/scripts/kestrel_config.json` | `duplicate_flagging.sort_by` only (#197) |
| `vntyper/modules/advntr/advntr_genotyping.py` | Docstring (#182); `LEN` summation (#192) |
| `vntyper/modules/shark/shark_filtering.py` | Docstring + warning (#187) |
| `docker/app/tasks.py` | `--cram` for CRAM uploads (#188) |
| `pyproject.toml` | `mypy_path` (#194); `branch = true` and `fail_under` (#196) |
| `AGENTS.md` | Refresh the stale LOC/coverage table |

---

# Phase 0 — Land #198

### Task 0: Merge PR #198 and open the successor branch

**Files:** none in this repo tree.

**Interfaces:**
- Produces: branch `fix/issue-181-197-followups` off the post-merge `main`; the golden-cohort baseline for Task 21.

- [ ] **Step 1: Confirm #198 is still green**

```bash
gh pr checks 198 --repo hassansaei/VNtyper
```

Expected: `CI Success` pass, all four `Unit Tests` pass, `Docs build (strict)` pass, `Docker Success` pass.

- [ ] **Step 2: Merge without squashing**

PRs on this repo are **not** squashed — every commit becomes permanent history.

```bash
gh pr merge 198 --repo hassansaei/VNtyper --merge
```

- [ ] **Step 3: Close the issues #198 resolves**

`#179` closes automatically via the `Closes #179` footer. `#190` does not — close it manually with the evidence:

```bash
gh issue close 190 --repo hassansaei/VNtyper --comment \
  "Fixed on #198 at 50d7968: \`cohort_summary.py:657\` now builds its environment with \`autoescape=select_autoescape([\"html\", \"xml\"])\`, and \`report_formatting.escaped_table_html\` replaces the three \`to_html(escape=False)\` calls with per-column exemptions. Pinned by \`tests/unit/test_cohort_summary_escaping.py\`."
```

- [ ] **Step 4: Create the successor branch**

```bash
git fetch origin
git switch -c fix/issue-181-197-followups origin/main
git log --oneline -1
```

Expected: the merge commit of #198.

- [ ] **Step 5: Record the baseline commit for the gate**

```bash
git rev-parse HEAD > /tmp/golden-baseline.sha
cat /tmp/golden-baseline.sha
```

This SHA is the "before" side in Task 21. Nothing else in this plan may use a different one.

---

# Phase 1 — Zero-risk issue closeouts (6 parallel lanes)

**Tasks 1, 2, 3, 4, 5, 6b, 7.** Task 6a (#188) is written below beside 6b for readability but **executes in phase 5**, not here.

Every task in this phase is docs, tests, or a guard that cannot fire on the shipped config. None can move a genotype. None needs the golden-cohort gate.

### Task 1: #181 — record the +1-frame biology as specification

**Lane:** 1 (owns `vntyper/scripts/scoring.py`, `tests/unit/test_scoring.py`, `docs/pipeline/scoring-and-confidence.md`)

**Files:**
- Modify: `vntyper/scripts/scoring.py:121-141` (docstring of `extract_frameshifts`)
- Modify: `tests/unit/test_scoring.py:225-232` (the characterisation banner)
- Modify: `docs/pipeline/scoring-and-confidence.md`
- Test: `tests/unit/test_scoring.py`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing importable. Task 2 writes an independent parity test; it must not import from this file.

**Decision being encoded**, quoted verbatim from @hassansaei on #181:

> "Yes — this was a considered MUC1-specific choice, not accidental asymmetry. Insertions of length 3n+1 and deletions of length 3n+2 are frame-equivalent (both Δ ≡ +1 mod 3). That is the pathogenic ADTKD-MUC1 reading frame that produces the toxic MUC1-fs neo-protein. […] only the +1 frame yields the disease-associated MUC1-fs product. The other frame has not been established as pathogenic in patients and is treated as unknown / not clinically identified for ADTKD-MUC1."

**The rule itself does not change.** Only what the tests claim about it changes.

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_scoring.py`:

```python
def test_a_one_bp_deletion_is_not_a_valid_frameshift_by_design():
    """A deletion of 1 bp is ``frameshift_amount == 1`` and is deliberately rejected.

    Specification, not characterisation. @hassansaei on #181: insertions of
    3n+1 and deletions of 3n+2 are both Delta = +1 mod 3, the pathogenic
    ADTKD-MUC1 reading frame that produces the toxic MUC1-fs neo-protein.
    The opposite pair shifts into a frame not established as pathogenic in
    patients. Rejecting a (3n+1)-bp deletion is therefore correct, not a
    lost call. Changing this rule requires a new decision on #181.
    """
    df = pd.DataFrame(
        {
            "direction": [-1, -1, 1, 1],
            "frameshift_amount": [1, 2, 1, 2],
        }
    )

    out = extract_frameshifts(df)

    assert out["is_valid_frameshift"].tolist() == [False, True, True, False]
```

- [ ] **Step 2: Run it and confirm it passes for the right reason**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/test_scoring.py::test_a_one_bp_deletion_is_not_a_valid_frameshift_by_design -v
```

Expected: PASS. This one pins existing behaviour, so it cannot fail first. **Prove it can fail** by temporarily flipping `scoring.py:151` to `== 1` and re-running:

Expected: FAIL, `assert [True, False, True, False] == [False, True, True, False]`. Revert the flip before continuing.

- [ ] **Step 3: Leave the section banner ALONE; upgrade only the frame-convention tests**

**Do not touch `tests/unit/test_scoring.py:225-232`.** That banner covers everything below it, including the sample-field parsing tests at `:235` and the direction-zero behaviour — none of which @hassansaei decided. Replacing it would ratify behaviour nobody signed off on, which is the precise failure this plan exists to catch in others.

Instead, add a one-line pointer under the existing banner, and put the decision in the **docstring of each frame-convention test individually**:

```python
# Note: the frameshift-validity tests below carry their own SPECIFICATION
# docstrings (issue #181, decided 2026-08-06). Everything else in this
# section remains characterisation, as this banner says.
```

The tests to upgrade are the ones asserting the 3n+1 / 3n+2 rule — the new test from Step 1 and any existing test whose assertion is *only* about `is_valid_frameshift` given `direction` and `frameshift_amount`. **List them by name in the commit message.** If a test asserts anything else as well, leave it characterisation and do not split it.

- [ ] **Step 4: Add the biology to the source docstring**

In `vntyper/scripts/scoring.py`, `extract_frameshifts`'s docstring, after the existing "Steps:" block:

```
    Why the asymmetry (issue #181):
      An insertion of 3n+1 bases and a deletion of 3n+2 bases are frame-
      equivalent -- both shift the reading frame by Delta = +1 (mod 3). That
      is the pathogenic ADTKD-MUC1 frame, which produces the toxic MUC1-fs
      neo-protein (classically exemplified by dupC). The opposite pair
      (insertion 3n+2 / deletion 3n+1) shifts into the other frame, which has
      not been established as pathogenic in patients. Rejecting a (3n+1)-bp
      deletion is therefore intentional, not a lost call.
```

- [ ] **Step 5: State it in the docs**

Add to `docs/pipeline/scoring-and-confidence.md`, in the frameshift section, the same explanation in prose plus a link to issue #181. Do **not** add a new page — this file is already in `mkdocs.yml:124`.

- [ ] **Step 6: Run the full unfiltered tier**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest -m unit tests/unit -q
```

Expected: all pass. Do not use `-k` — `test_marker_hygiene.py` fails under any filtered run.

- [ ] **Step 7: Hand off to the orchestrator**

Report the three modified paths and paste the Step 2 induced-failure output. **Do not commit.**

Orchestrator commit:

```bash
git add vntyper/scripts/scoring.py tests/unit/test_scoring.py docs/pipeline/scoring-and-confidence.md
git commit -m "docs(scoring): record the +1-frame rule as specification, per #181

Closes #181"
```

---

### Task 2: #182 — adVNTR frame filter specification and cross-module parity

**Lane:** 2 (owns `vntyper/modules/advntr/advntr_genotyping.py` docstrings, `tests/unit/test_advntr_frameshift_filter.py`, and creates `tests/unit/test_frameshift_convention_parity.py`)

**Files:**
- Create: `tests/unit/test_frameshift_convention_parity.py`
- Modify: `vntyper/modules/advntr/advntr_genotyping.py:174-178` and `:213-217` (docstrings/comments around `del_frame`/`ins_frame`)
- Modify: `tests/unit/test_advntr_frameshift_filter.py` (banner)
- Test: both files above

**Interfaces:**
- Consumes: `vntyper.scripts.scoring.extract_frameshifts`, `vntyper.modules.advntr.advntr_genotyping` module constants. Read-only — Task 1 owns `scoring.py`, this task must not edit it.
- Produces: nothing importable.

**Decision**, quoted from @hassansaei on #182:

> "Yes, keep the same 3n+1 / 3n+2 rule for adVNTR as for Kestrel (#181). This is intentional shared convention, not something to relax independently. […] We should keep the filtering logic harmonized between Kestrel and adVNTR, as already implemented since v1.3."

The risk his answer creates is **drift**: two modules, one convention, nothing asserting they agree. This task closes that.

- [ ] **Step 1: Write the failing parity test**

Create `tests/unit/test_frameshift_convention_parity.py`:

```python
"""One frameshift convention, two modules — asserted to stay in agreement.

@hassansaei on #182: "We should keep the filtering logic harmonized between
Kestrel and adVNTR." Kestrel encodes the rule as a pair of pandas conditions
in ``scoring.extract_frameshifts``; adVNTR encodes it as two numpy arrays in
``advntr_genotyping``. Nothing previously asserted the two agree, so either
could have been changed alone. This file is that assertion.
"""

import numpy as np
import pandas as pd
import pytest

from vntyper.modules.advntr import advntr_genotyping
from vntyper.scripts.scoring import extract_frameshifts

pytestmark = pytest.mark.unit

MAX_FRAMESHIFT = 100
FRAMESHIFT_MULTIPLIER = 3


def _kestrel_accepts(direction, amount):
    """True when Kestrel would mark this (direction, frameshift_amount) valid."""
    out = extract_frameshifts(pd.DataFrame({"direction": [direction], "frameshift_amount": [amount]}))
    return bool(out["is_valid_frameshift"].iloc[0])


def test_kestrel_accepts_insertions_at_3n_plus_1_and_deletions_at_3n_plus_2():
    assert _kestrel_accepts(1, 1) is True
    assert _kestrel_accepts(1, 2) is False
    assert _kestrel_accepts(-1, 2) is True
    assert _kestrel_accepts(-1, 1) is False


def test_advntr_frame_sets_encode_the_same_convention_as_kestrel():
    """adVNTR's ins_frame/del_frame residues must match Kestrel's amounts."""
    ins_frame = (np.arange(MAX_FRAMESHIFT) * FRAMESHIFT_MULTIPLIER + 1).astype(int)
    del_frame = (np.arange(MAX_FRAMESHIFT) * FRAMESHIFT_MULTIPLIER + 2).astype(int)

    assert {int(v) % 3 for v in ins_frame} == {1}
    assert {int(v) % 3 for v in del_frame} == {2}

    # And those residues are exactly the ones Kestrel accepts.
    assert _kestrel_accepts(1, int(ins_frame[0]) % 3) is True
    assert _kestrel_accepts(-1, int(del_frame[0]) % 3) is True


def test_the_shipped_advntr_config_uses_the_multiplier_the_convention_requires():
    """A multiplier other than 3 would silently decouple the two modules."""
    config = advntr_genotyping.advntr_config["advntr_settings"]
    assert config["frameshift_multiplier"] == FRAMESHIFT_MULTIPLIER
```

- [ ] **Step 2: Run it and confirm it fails**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/test_frameshift_convention_parity.py -v
```

Expected: the third test FAILs with `KeyError` if the config key is named differently, or PASSes. If it passes, **induce a failure** by temporarily setting `frameshift_multiplier` to `4` in a local copy and confirm the test catches it. If the key name differs from `frameshift_multiplier`, read `vntyper/modules/advntr/advntr_config.json` and correct the test to the real name — do **not** change the config.

- [ ] **Step 3: Make it pass**

Correct any key-name mismatch found in Step 2. No production code changes.

- [ ] **Step 4: Run it and confirm it passes**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/test_frameshift_convention_parity.py -v
```

Expected: 3 passed.

- [ ] **Step 5: Upgrade individual test docstrings — not the file's framing**

`tests/unit/test_advntr_frameshift_filter.py` also characterises parsing columns (`:80`) and import-time configuration globals (`:212`), which #182 does not decide. **Do not rewrite the module docstring or any section banner.** Upgrade only the docstrings of tests whose assertion is solely about `ins_frame` / `del_frame` membership, and name them in the commit message.

In `advntr_genotyping.py` near lines 176 and 215, replace the bare comment with:

```python
    # Pathogenic-frame filter (#182, decided 2026-08-06). Deletions must be
    # 3n+2 and insertions 3n+1: both are Delta = +1 mod 3, the ADTKD-MUC1
    # frame that yields the toxic MUC1-fs product. This is the SAME rule
    # Kestrel applies in scoring.extract_frameshifts (#181); the two are
    # asserted to agree by tests/unit/test_frameshift_convention_parity.py.
```

- [ ] **Step 6: Run the full unfiltered tier**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest -m unit tests/unit -q
```

Expected: all pass, count up by 3.

- [ ] **Step 7: Hand off**

Orchestrator commit:

```bash
git add tests/unit/test_frameshift_convention_parity.py tests/unit/test_advntr_frameshift_filter.py \
        vntyper/modules/advntr/advntr_genotyping.py
git commit -m "test(advntr): assert Kestrel and adVNTR share one frameshift convention

Closes #182"
```

---

### Task 3: #183 — confidence-tier precedence is specified as last-wins

**Lane:** 3 (owns `tests/unit/test_confidence_boundaries.py`; **shares `docs/pipeline/scoring-and-confidence.md` with Task 1, so the two must run serially — Task 1 first**)

**Files:**
- Modify: `tests/unit/test_confidence_boundaries.py:1-20` (module docstring) and `:266` (the per-test caveat)
- Modify: `docs/pipeline/scoring-and-confidence.md:75` — it still calls tier precedence **"unspecified"**. #183 decided it. The first draft of the spec claimed this page already agreed; it does not. Replace with the decided last-wins behaviour, quoting @hassansaei and linking #183.
- Test: `tests/unit/test_confidence_boundaries.py`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing. **Must not edit `vntyper/scripts/confidence_assignment.py`** — Task 17 owns it in phase 4.

**Decision**, quoted from @hassansaei on #183:

> "Keep the current 2.x last-wins logic — do not restore the absolute region-depth ≤200 cap from 1.3. In practice it is very unlikely that a variant with region depth <200 is later promoted to High_Precision by a subsequent rule. Where that pattern can appear, it is mostly for early (beginning) and late conserved motifs; we already have a flagging rule when Depth_Score is far from ~0.5 (50%) […] the intentional behaviour going forward is the 2.x sequential assignment as implemented today."

- [ ] **Step 1: Write the failing test**

Append to `tests/unit/test_confidence_boundaries.py` the exact case from the issue body:

```python
def test_a_low_region_depth_row_is_deliberately_not_capped_at_low_precision():
    """alt=50 on a 150-read region reports High_Precision, and that is intended.

    Specification, not characterisation. v1.3 made region depth an absolute
    first-wins cap; 2.x applies the tiers as sequential ``df.loc`` assignments,
    so ``cond1``'s demotion is overwritten by a later match. @hassansaei on
    #183 decided this is the intended behaviour going forward and that the
    1.3 cap must not be restored: the pattern is rare, and the cases where it
    appears are already caught by the Depth_Score flagging rule.
    """
    df = pd.DataFrame(
        {
            "Estimated_Depth_AlternateVariant": [50],
            "Estimated_Depth_Variant_ActiveRegion": [150],
        }
    )

    out = calculate_depth_score_and_assign_confidence(df, _shipped_config())

    assert out["Confidence"].iloc[0] == "High_Precision"
    assert out["Estimated_Depth_Variant_ActiveRegion"].iloc[0] <= 200
```

Use whatever helper this file already provides for the shipped config; if it has none, read it with `json.load(open("vntyper/scripts/kestrel_config.json"))`.

- [ ] **Step 2: Run it**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/test_confidence_boundaries.py -v
```

Expected: PASS. **Prove it can fail** by temporarily moving the `cond1` assignment in `confidence_assignment.py` to last; expected FAIL with `assert 'Low_Precision' == 'High_Precision'`. Revert immediately — you do not own that file.

- [ ] **Step 3: Amend the module docstring — carefully, because #184 changes one of its cells**

The existing text at `:12-16` says everything in the file is characterisation. That is now *too broad in one direction and exactly right in another*: the ordering is decided (#183), the 54-cell threshold matrix is not, and Task 17 will change one of those cells. Amend it to distinguish the two rather than flipping the whole file:

```
The 54-cell boundary matrix below is **characterisation**: it records the
threshold arithmetic so a changed comparison operator becomes visible, and
makes no claim that any individual cutoff is right.

The *ordering* is different. @hassansaei decided on #183 (2026-08-06) that
2.x's last-wins sequential assignment is the intended behaviour and that
1.3's absolute region-depth cap must not be restored. Tests that pin the
order of the ``df.loc`` assignments are therefore **specification**, and
changing that order requires a new decision on #183.
```

Update the per-test caveat at `:266` the same way.

- [ ] **Step 4: Run the full unfiltered tier**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest -m unit tests/unit -q
```

Expected: all pass, count up by 1.

- [ ] **Step 5: Hand off**

Orchestrator commit:

```bash
git add tests/unit/test_confidence_boundaries.py
git commit -m "test(confidence): specify last-wins tier precedence, per #183

Closes #183"
```

---

### Task 4: #186 — refuse the armed uniform-filtering combination

**Lane:** 4 (owns `vntyper/scripts/motif_processing.py` in phase 1, creates `tests/unit/test_motif_config_guard.py`)

**Files:**
- Create: `tests/unit/test_motif_config_guard.py`
- Modify: `vntyper/scripts/motif_processing.py:313-320`
- Test: `tests/unit/test_motif_config_guard.py`

**Interfaces:**
- Consumes: nothing.
- Produces: the guard block inside `motif_correction_and_annotation`. **Tasks 9 and 10 refactor this same file in phase 2 and must preserve this guard.**

**Decision**, quoted from @hassansaei on #186:

> "**This is a high risk issue.** It can remove the true positive from results. Especially dangerous if someone enables uniform filtering while motifs_for_alt_gg is still []. […] Do not repopulate `motifs_for_alt_gg` with `["X"]` for the active path. Do not change the already-implemented v2.x motif logic on this point. Leave `use_uniform_filtering` off unless/until that branch is redesigned and validated; with an empty allowlist it would delete every GG, including canonical dupC."

**This guard is an inference from his hazard sentence, not a quoted instruction.** He described the danger and prescribed no mitigation. It is genotype-neutral by construction: the shipped config sets `use_uniform_filtering: false`, so the guard cannot fire. The commit message must say this in those words.

**Do not** repopulate `motifs_for_alt_gg`. **Do not** delete the uniform branch or `tests/unit/test_motif_filtering_issue_136.py` — he wrote "unless/until that branch is redesigned and validated", which means it stays.

- [ ] **Step 1: Write the failing test**

Create `tests/unit/test_motif_config_guard.py`:

```python
"""The uniform-filtering branch may not run with an empty GG allowlist.

@hassansaei on #186: enabling ``use_uniform_filtering`` while
``motifs_for_alt_gg`` is ``[]`` "would delete every GG, including canonical
dupC" -- the pathogenic call this tool exists to find. The shipped config
sets the flag to false, so this guard cannot fire in production; it exists
so that flipping the flag fails loudly instead of silently deleting calls.
"""

import copy
import json
from pathlib import Path

import pandas as pd
import pytest

from vntyper.scripts.motif_processing import motif_correction_and_annotation

pytestmark = pytest.mark.unit

CONFIG_PATH = Path("vntyper/scripts/kestrel_config.json")


def _shipped_config():
    return json.loads(CONFIG_PATH.read_text())


def _one_row():
    return pd.DataFrame(
        {
            "Motifs": ["X-5"],
            "POS": [67],
            "REF": ["C"],
            "ALT": ["CC"],
            "Estimated_Depth_AlternateVariant": [50],
            "Estimated_Depth_Variant_ActiveRegion": [5000],
            "Depth_Score": [0.01],
            "Confidence": ["High_Precision"],
        }
    )


def test_the_shipped_config_leaves_uniform_filtering_off():
    """If this fails, the guard below is no longer inert and #186 is reopened."""
    config = _shipped_config()
    assert config["motif_filtering"]["use_uniform_filtering"] is False
    assert config["motif_filtering"]["motifs_for_alt_gg"] == []


def test_uniform_filtering_with_an_empty_gg_allowlist_is_refused():
    config = _shipped_config()
    config["motif_filtering"]["use_uniform_filtering"] = True

    with pytest.raises(ValueError, match="motifs_for_alt_gg"):
        motif_correction_and_annotation(_one_row(), pd.DataFrame({"Motif": [], "Motif_sequence": []}), config)


def test_uniform_filtering_with_a_populated_allowlist_is_permitted():
    """The guard blocks the armed combination only, not the branch itself."""
    config = _shipped_config()
    config["motif_filtering"]["use_uniform_filtering"] = True
    config["motif_filtering"]["motifs_for_alt_gg"] = ["X"]

    out = motif_correction_and_annotation(_one_row(), pd.DataFrame({"Motif": ["5"], "Motif_sequence": ["ACGT"]}), config)

    assert "motif_filter_pass" in out.columns
```

- [ ] **Step 2: Run it and confirm it fails**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/test_motif_config_guard.py -v
```

Expected: `test_uniform_filtering_with_an_empty_gg_allowlist_is_refused` FAILs with `DID NOT RAISE <class 'ValueError'>`. The other two pass.

- [ ] **Step 3: Add the guard**

In `vntyper/scripts/motif_processing.py`, immediately after line 320 (`exclude_motifs_combined = mf.get(...)`):

```python
    # #186: with an empty allowlist the uniform branch treats *every* GG alt as
    # disallowed and deletes it, including the canonical dupC representation at
    # POS 67. @hassansaei: "with an empty allowlist it would delete every GG,
    # including canonical dupC, which is exactly the high-risk failure mode we
    # are avoiding." The shipped config sets the flag false, so this cannot
    # fire in production; it exists so that flipping it fails loudly.
    if use_uniform_filtering and not motifs_for_alt_gg:
        msg = (
            "use_uniform_filtering is enabled while motifs_for_alt_gg is empty. "
            "The uniform branch applies motifs_for_alt_gg as an unconditional allowlist, "
            "so an empty list deletes every GG alternate -- including the canonical dupC "
            "call at POS 67. Populate motifs_for_alt_gg or leave use_uniform_filtering off. "
            "See issue #186."
        )
        logger.error(msg)
        raise ValueError(msg)
```

- [ ] **Step 4: Run it and confirm it passes**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/test_motif_config_guard.py -v
```

Expected: 3 passed.

- [ ] **Step 5: Confirm the 588-LOC #136 suite is unaffected**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/test_motif_filtering_issue_136.py -v
```

Expected: all pass. That file imports `_apply_uniform_filtering_right_motif` directly rather than going through `motif_correction_and_annotation`, so the guard does not reach it. If any test there fails, **stop and report** — the guard is in the wrong place.

- [ ] **Step 6: Run the full unfiltered tier**

Expected: all pass, count up by 3.

- [ ] **Step 7: Hand off**

Orchestrator commit:

```bash
git add vntyper/scripts/motif_processing.py tests/unit/test_motif_config_guard.py
git commit -m "fix(motif): refuse uniform filtering with an empty GG allowlist

Inert on the shipped config, which sets use_uniform_filtering false. The guard
exists so that flipping the flag fails loudly rather than deleting every GG alt
including canonical dupC. Inferred from the hazard @hassansaei named on #186,
not from an instruction he gave; the config and the motif logic are unchanged.

Closes #186"
```

---

### Task 5: #187 — SHARK is sequence-based; say so and warn

**Lane:** 5 (owns `vntyper/modules/shark/shark_filtering.py`, `tests/unit/test_shark_filtering.py`, `docs/pipeline/optional-modules.md`)

**Files:**
- Modify: `vntyper/modules/shark/shark_filtering.py:37,50` (signature docstring) and the function body
- Modify: `tests/unit/test_shark_filtering.py` (`TestReferenceAssemblyIsAccceptedAndIgnored`)
- Modify: `docs/pipeline/optional-modules.md`
- Test: `tests/unit/test_shark_filtering.py`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing. The `reference_assembly` parameter keeps its name and default for API compatibility.

**Decision**, quoted from @hassansaei on #187:

> "SHARK is sequence-based, not coordinate-based; keep one MUC1 region FASTA; document that assembly does not change SHARK; optionally warn or ignore reference_assembly there instead of building a second FASTA unless GRCh37/GRCh38 sequences at MUC1 VNTR differ enough to matter, like the number of motif could impact the number of reads being retained by the SHARK model."

**No hg38 region FASTA is built.** The parameter stays, is documented as inert, and warns.

- [ ] **Step 1: Write the failing test**

Add to `tests/unit/test_shark_filtering.py`:

The entry point is `run_shark_filter(fastq_1, fastq_2, output_dir, config, main_config, sample_name, reference_assembly="hg19", threads=4)` at `shark_filtering.py:31`.

```python
@pytest.mark.parametrize("assembly", ["hg38", "GRCh38", "hg38_ensembl"])
def test_a_non_hg19_assembly_is_accepted_and_warned_about(assembly, caplog, tmp_path, monkeypatch):
    """SHARK filters against one sequence FASTA regardless of assembly.

    Specification (#187). @hassansaei: "SHARK is sequence-based, not
    coordinate-based; keep one MUC1 region FASTA; document that assembly does
    not change SHARK; optionally warn or ignore reference_assembly there
    instead of building a second FASTA." The parameter is kept for API
    compatibility and now says so at runtime.
    """
    caplog.set_level(logging.WARNING, logger="vntyper.modules.shark.shark_filtering")
    monkeypatch.setattr(shark_filtering, "run_command", lambda *a, **k: (True, "", ""))

    run_shark_filter(
        fastq_1=str(tmp_path / "r1.fastq"),
        fastq_2=str(tmp_path / "r2.fastq"),
        output_dir=str(tmp_path),
        config={"shark_settings": {"muc1_region_fasta": str(tmp_path / "muc1.fa")}},
        main_config={"tools": {"shark": "shark"}},
        sample_name="sample",
        reference_assembly=assembly,
    )

    assert any("reference_assembly" in r.message and assembly in r.message for r in caplog.records)


def test_hg19_does_not_warn():
    """The warning marks a mismatch, not every call. Same invocation, hg19."""
    # ... same body with reference_assembly="hg19"
    assert not [r for r in caplog.records if "reference_assembly" in r.message]


def test_the_region_fasta_is_the_same_one_at_every_assembly(tmp_path, monkeypatch):
    """The substantive claim: assembly does not select a region.

    Captures the command SHARK is actually given at hg19 and at hg38 and
    asserts the `-r` operand is byte-identical. A warning alone would not
    prove the parameter is inert; this does.
    """
    captured = []
    monkeypatch.setattr(shark_filtering, "run_command", lambda cmd, *a, **k: captured.append(cmd) or (True, "", ""))
    # ... invoke twice, once per assembly, then:
    assert captured[0] == captured[1]
```

Match `run_command`'s real return contract — read `vntyper/scripts/utils.py` rather than assuming the `(ok, stdout, stderr)` shape above. `pytest.ini` sets `log_level = DEBUG`, so a naive `caplog` assertion passes against invisible logs: pin the level **and** the logger name, as above.

- [ ] **Step 2: Run it and confirm it fails**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/test_shark_filtering.py -v
```

Expected: FAIL, `assert False` — no warning is emitted today.

- [ ] **Step 3: Emit the warning**

In `shark_filtering.py`, at the top of the function body:

```python
    if reference_assembly and reference_assembly.lower() not in ("hg19", "grch37"):
        logger.warning(
            f"reference_assembly={reference_assembly!r} does not select a region for SHARK. "
            "SHARK is sequence-based, not coordinate-based: it filters against the single "
            "MUC1 region FASTA regardless of assembly. The parameter is retained for API "
            "compatibility only. See issue #187."
        )
```

- [ ] **Step 4: Run it and confirm it passes**

Expected: PASS.

- [ ] **Step 5: Update the docstring and the class name's framing**

`shark_filtering.py:50` currently reads `reference_assembly (str): Reference assembly (hg19 or hg38).` Replace with:

```
        reference_assembly (str): Accepted for API compatibility and **ignored**.
            SHARK matches k-mers against a single MUC1 region FASTA and does not
            select a region by coordinate, so the assembly does not change what
            it retains. Passing anything other than hg19/GRCh37 logs a warning.
            See issue #187.
```

Rename `TestReferenceAssemblyIsAccceptedAndIgnored`'s docstring from "characterisation of a live bug" to specification citing #187 and quoting the sentence above. (Leave the class name's spelling alone — renaming it is churn.)

- [ ] **Step 6: Document it**

Add a short subsection to `docs/pipeline/optional-modules.md` under SHARK stating that `--reference-assembly` does not affect SHARK, why, and linking #187. The file is already in `mkdocs.yml`; do not add a page.

- [ ] **Step 7: Strict docs build**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" make docs-build
```

Expected: exit 0.

- [ ] **Step 8: Run the full unfiltered tier, then hand off**

Orchestrator commit:

```bash
git add vntyper/modules/shark/shark_filtering.py tests/unit/test_shark_filtering.py docs/pipeline/optional-modules.md
git commit -m "docs(shark): state that reference_assembly does not select a region

Closes #187"
```

---

### Task 6a (**gated phase**) / Task 6b (**phase 1**) — the docker/app lane

**Split after review.** **Task 6b (#194, mypy) stays in phase 1** — it is annotations and configuration, zero genotype risk. **Task 6a (#188, CRAM) moves to the gated phase** and grows two more fixes, because the "one-line change" introduces a regression on its own:

1. **The sample name would regress to the literal `"sample"`.** `cli_handlers.py:276-281` derives it from `args.bam` or `args.fastq1` and **never** from `args.cram`. Today a CRAM reaches the pipeline through `--bam`, so it gets its file stem; after switching to `--cram` every CRAM run without an explicit `--sample-name` would be named `sample`. That name reaches the report and the output filenames.
2. **The generated index name is wrong for CRAM.** `tasks.py:104` falls back to `f"{bam_path}.bai"`, but `samtools index` on a CRAM writes `.crai`. The cleanup path then removes a file that does not exist and leaves the one that does.

**The golden cohort contains no CRAM input** — its own page says so under "What this gate does not cover". So the gate cannot attest #188 either. **If no CRAM fixture can be produced, #188 parks on the same trigger as #191** rather than shipping on an argument. Decide that before starting, not after.

**Task 6a and 6b share `docker/app/` and must never run concurrently.**

**Files:**
- Modify: `docker/app/tasks.py:143`
- Modify: `docker/app/*.py` (22 annotation fixes)
- Modify: `pyproject.toml` `[tool.mypy]`
- Test: `tests/unit/web/test_job_submission_lifecycle.py`

**Interfaces:**
- Consumes: `docker/app/uploads.py`'s validated extension allowlist (`.bam`/`.cram` for the alignment).
- Produces: nothing importable. **Must not touch `pyproject.toml`'s `[tool.coverage.*]`** — Task 16 owns that.

**#188 is a one-line fix, not the three-file change the issue claims.** Verified: `--cram` already exists at `vntyper/scripts/cli_parser.py:92`, `cli_handlers.py:302` threads it into `run_pipeline`, and `pipeline.py:180/221/313/328/465` branches on it. Only `docker/app/tasks.py:143` hardcodes `--bam`. **Correct the issue text when closing it.**

- [ ] **Step 1: Write the failing test**

In `tests/unit/web/test_job_submission_lifecycle.py`:

The command is built **inline** inside `run_vntyper_job` at `docker/app/tasks.py:135-159`, so it cannot be tested without Celery. Extract it first, verbatim, to a module-level function:

```python
def build_vntyper_command(
    alignment_path: str,
    output_dir: str,
    thread: int,
    reference_assembly: str,
    fast_mode: bool = False,
    keep_intermediates: bool = False,
    archive_results: bool = False,
    advntr_mode: bool = False,
) -> list[str]:
    """Assemble the vntyper CLI invocation for one job.

    Extracted from run_vntyper_job so the flag selection is testable without a
    Celery worker. Behaviour is unchanged apart from the alignment flag (#188).
    """
```

Then the tests:

```python
def test_a_cram_upload_is_passed_to_the_cli_as_cram_not_bam():
    """The web service accepts .cram; it must not hand it to the BAM flag.

    Issue #188. The CLI has had a --cram flag since cli_parser.py:92,
    cli_handlers.py:302 threads it into run_pipeline, and pipeline.py:180
    branches on it. Only the task layer hardcoded --bam, so every accepted
    CRAM took the BAM code path.
    """
    command = build_vntyper_command(
        alignment_path="/data/sample.cram", output_dir="/out", thread=4, reference_assembly="hg38"
    )

    assert "--cram" in command
    assert "--bam" not in command
    assert command[command.index("--cram") + 1] == "/data/sample.cram"


def test_a_bam_upload_still_uses_the_bam_flag():
    command = build_vntyper_command(
        alignment_path="/data/sample.bam", output_dir="/out", thread=4, reference_assembly="hg38"
    )

    assert "--bam" in command
    assert "--cram" not in command


def test_the_flag_is_chosen_case_insensitively():
    """uploads.py accepts SAMPLE.CRAM, so the flag choice must match on case."""
    command = build_vntyper_command(
        alignment_path="/data/SAMPLE.CRAM", output_dir="/out", thread=4, reference_assembly="hg38"
    )

    assert "--cram" in command


def test_every_other_flag_is_unchanged_by_the_extraction():
    """The extraction must be behaviour-preserving apart from the flag choice."""
    command = build_vntyper_command(
        alignment_path="/data/s.bam",
        output_dir="/out",
        thread=8,
        reference_assembly="hg19",
        fast_mode=True,
        keep_intermediates=True,
        archive_results=True,
        advntr_mode=True,
    )

    assert command[:6] == ["conda", "run", "-n", "vntyper", "vntyper", "pipeline"]
    assert "--fast-mode" in command
    assert "--keep-intermediates" in command
    assert "--archive-results" in command
    assert command[command.index("--extra-modules") + 1] == "advntr"
    assert command[command.index("--advntr-max-coverage") + 1] == "300"
    assert command[command.index("--threads") + 1] == "8"
```

- [ ] **Step 2: Run it and confirm it fails**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/web/test_job_submission_lifecycle.py -v
```

Expected: the CRAM tests FAIL with `assert '--cram' in [...]`.

- [ ] **Step 3: Write the two failing tests the review found**

```python
def test_a_cram_run_derives_its_sample_name_from_the_file_stem():
    """cli_handlers.py:276-281 derives the name from --bam or --fastq1 only.

    Switching the worker to --cram without this would name every CRAM run
    "sample" -- the fallback at cli_handlers.py:283. Today CRAM arrives via
    --bam and gets its stem, so the one-line #188 fix would INTRODUCE this
    regression. The two must land together.
    """
    args = _parsed_args(cram="/data/NA12878.cram", sample_name=None)
    assert _derive_sample_name(args) == "NA12878"


def test_a_cram_upload_with_no_index_falls_back_to_the_crai_name(tmp_path):
    """tasks.py:104 falls back to f"{bam_path}.bai" for every format.

    samtools index writes .crai for a CRAM, so the fallback names a file that
    is never created: cleanup then removes nothing and leaves the real index
    on the shared volume.
    """
    assert _resolve_index_path("/data/s.cram", None) == "/data/s.cram.crai"
    assert _resolve_index_path("/data/s.bam", None) == "/data/s.bam.bai"
    assert _resolve_index_path("/data/s.cram", "/data/given.crai") == "/data/given.crai"
```

- [ ] **Step 4: Run them and confirm they fail**

Expected: `assert 'sample' == 'NA12878'` and `assert '/data/s.cram.bai' == '/data/s.cram.crai'`.

- [ ] **Step 5: Make all three fixes together**

At `docker/app/tasks.py:143`:

```python
    alignment_flag = "--cram" if str(alignment_path).lower().endswith(".cram") else "--bam"
```

At `docker/app/tasks.py:104`, choose the index suffix by format:

```python
    # samtools index writes .crai for a CRAM and .bai for a BAM. The fallback
    # must name the file that will actually exist, or cleanup misses it (#188).
    index_path = index_path or f"{bam_path}{'.crai' if str(bam_path).lower().endswith('.cram') else '.bai'}"
```

At `vntyper/scripts/cli_handlers.py:277`, add the CRAM arm **before** the FASTQ one so it mirrors the BAM case:

```python
        if args.bam:
            sample_name_val = Path(args.bam).stem
        elif args.cram:
            sample_name_val = Path(args.cram).stem
            logger.debug(f"sample_name set from CRAM file: {sample_name_val}")
        elif args.fastq1:
```

- [ ] **Step 6: Run them and confirm they pass**

Expected: 6 passed.

- [ ] **Step 7: Run a CRAM end-to-end through the worker path**

The golden cohort has no CRAM sample, so it cannot attest this. Convert one cohort BAM to CRAM against the matching reference, submit it through the worker path, and assert: the command carries `--cram`, the derived sample name is the stem, the `.crai` is produced and cleaned up, and the reported genotype matches the BAM run of the same sample.

**If a CRAM fixture cannot be produced, stop.** Park #188 on the #191 trigger, revert Steps 5–6, and report. Do not ship it on the argument that the code path "should" work.

- [ ] **Step 5: Turn on mypy over tests/ (#194)**

Read the three preconditions already recorded as comments in `pyproject.toml`'s `[tool.mypy]` block. Add `mypy_path = "docker"`, then:

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" make type-check-all
```

Expected: ~22 errors in `docker/app`. Fix each — ASGI callables typed as protocols rather than concretely, `create_cohort_record` returning a `TypedDict` rather than `dict[str, str | None]`. **Fix the code; do not add `# type: ignore` and do not relax a mypy setting.**

- [ ] **Step 6: Re-run type checking until clean**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" make type-check-all
```

Expected: `Success: no issues found`.

- [ ] **Step 7: Run the full unfiltered tier, then hand off**

Orchestrator commits, **two separate commits**:

```bash
git add docker/app/tasks.py tests/unit/web/test_job_submission_lifecycle.py
git commit -m "fix(web): pass --cram for CRAM uploads instead of --bam

Closes #188"

git add pyproject.toml docker/app/
git commit -m "chore(types): gate tests/ under mypy and fix docker/app annotations

Closes #194"
```

---

### Task 7: #197 — duplicate ordering falls back to Depth_Score

**Lane:** 7 (owns `vntyper/scripts/kestrel_config.json`, `vntyper/scripts/flagging.py`, `tests/unit/test_flagging.py`)

**Files:**
- Modify: `vntyper/scripts/kestrel_config.json` (`duplicate_flagging.sort_by` only)
- Modify: `vntyper/scripts/flagging.py:217` (`sort_values` call)
- Modify: `tests/unit/test_flagging.py` — `TestDuplicateSortColumnMustExist`, **and `:400`, which pins the old three-key sort and will fail otherwise**
- Modify: `docs/pipeline/flagging.md:38` — still documents the three-key sort; update it to the one-key rule, quoting @hassansaei and linking #197
- Test: `tests/unit/test_flagging.py`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing importable.

**Decision**, quoted from @hassansaei on #197:

> "Motif is not important for duplicate ordering. Fall back to the 1.3 Depth_Score-only rule: sort by `Depth_Score` descending when marking duplicates. Keep the 2.x behaviour of **flagging** duplicates (`Potential_Duplicate`) rather than dropping rows. Do not use `Motifs` or `Motif` as a sort key. Leave `duplicate_flagging.enabled` as `false` in the shipped config (as now). We have already tested with this setup. I do not know what will happen if we turn it on!"

**This is the one `kestrel_config.json` edit in the plan.** It is provably inert: `duplicate_flagging.enabled` stays `false`, and `kestrel_genotyping.py:588` only calls `add_flags` with a duplicates config when that flag is true. No genotype can move.

`kind="stable"` is an implementation detail beyond his decision — with a single sort key, pandas' default quicksort makes *which* row is flagged non-deterministic. Say so in the commit message.

- [ ] **Step 1: Write the failing test**

In `tests/unit/test_flagging.py`:

```python
def test_duplicate_ordering_uses_depth_score_only_and_no_motif_column():
    """#197: Motif is not a duplicate sort key; Depth_Score descending is.

    Specification. @hassansaei: "Fall back to the 1.3 Depth_Score-only rule
    [...] Do not use Motifs or Motif as a sort key." The previous config named
    the plural `Motifs`, which does not exist at step 6.5 and raised KeyError.
    """
    config = json.loads(Path("vntyper/scripts/kestrel_config.json").read_text())
    sort_by = config["duplicate_flagging"]["sort_by"]

    assert sort_by == [{"column": "Depth_Score", "ascending": False}]
    assert all(entry["column"] not in ("Motif", "Motifs") for entry in sort_by)


def test_duplicate_flagging_stays_disabled_in_the_shipped_config():
    """@hassansaei: "Leave duplicate_flagging.enabled as false (as now)."."""
    config = json.loads(Path("vntyper/scripts/kestrel_config.json").read_text())
    assert config["duplicate_flagging"]["enabled"] is False


def test_the_flagged_row_is_deterministic_when_depth_scores_tie():
    """A single-key sort_values is not stable; ties must not reorder.

    Implementation detail, not part of the #197 decision: with Depth_Score as
    the only key, pandas' default quicksort makes which tied row keeps the
    first position arbitrary, and therefore which row gets flagged.
    """
    df = pd.DataFrame(
        {
            "REF": ["C", "C", "C"],
            "ALT": ["CC", "CC", "CC"],
            "Depth_Score": [0.5, 0.5, 0.5],
            "POS": [10, 20, 30],
            "Flag": ["Not flagged", "Not flagged", "Not flagged"],
        }
    )

    first = mark_potential_duplicates(df.copy(), ["REF", "ALT"], ["Depth_Score"], [False])
    second = mark_potential_duplicates(df.copy(), ["REF", "ALT"], ["Depth_Score"], [False])

    assert first["Flag"].tolist() == second["Flag"].tolist()
    assert first.loc[first["POS"] == 10, "Flag"].iloc[0] == "Not flagged"
```

Match `mark_potential_duplicates`'s real signature at `flagging.py:179`.

- [ ] **Step 2: Run it and confirm it fails**

Expected: the first test FAILs — the shipped config still has the three-key `sort_by` including `Motifs`.

- [ ] **Step 3: Edit the config**

In `vntyper/scripts/kestrel_config.json`, replace `duplicate_flagging.sort_by` with:

```json
      "sort_by": [
        { "column": "Depth_Score", "ascending": false }
      ],
```

Leave `enabled`, `flag_name` and `group_by` untouched.

- [ ] **Step 4: Make the sort stable**

At `flagging.py:217`:

```python
    # kind="stable" because sort_by is now a single key (#197): pandas' default
    # quicksort would make which tied row keeps the first position -- and
    # therefore which row is flagged -- arbitrary between runs.
    df_copy.sort_values(by=sort_cols, ascending=sort_ascending, inplace=True, kind="stable")
```

- [ ] **Step 5: Run it and confirm it passes**

Expected: 3 passed.

- [ ] **Step 6: Update `TestDuplicateSortColumnMustExist`**

That class currently asserts both the `KeyError` and the old config value. The `KeyError` assertion must now be rewritten to construct the stale config explicitly rather than reading the shipped one, so it still documents the trap without depending on a value that has changed. Cite #197 in the docstring.

- [ ] **Step 7: Prove the config hash changed in exactly one key**

```bash
git diff --stat vntyper/scripts/kestrel_config.json
git diff vntyper/scripts/kestrel_config.json
```

Expected: only lines inside `duplicate_flagging.sort_by`.

- [ ] **Step 8: Run the full unfiltered tier, then hand off**

Orchestrator commit:

```bash
git add vntyper/scripts/kestrel_config.json vntyper/scripts/flagging.py tests/unit/test_flagging.py
git commit -m "fix(config): order duplicate flagging by Depth_Score alone, per #197

@hassansaei on #197: \"Fall back to the 1.3 Depth_Score-only rule [...] Do not
use Motifs or Motif as a sort key. Leave duplicate_flagging.enabled as false.\"

Provably inert: enabled stays false, and kestrel_genotyping.py:588 only reaches
mark_potential_duplicates when it is true. No genotype can move.

kind=\"stable\" is beyond his decision and is an implementation detail: with one
sort key, the default quicksort makes which tied row is flagged arbitrary.

Closes #197"
```

---

# Phase 2 — Coverage gaps (6 parallel lanes)

**Tasks 8, 11, 12, 13, 14, 15.** Tasks 9 and 10 are written below in this section because they belong to the same body of coverage work, but **they execute in phase 5, gated** — they touch the code that computes `motif_filter_pass`. Do not run them here.

No task in this phase changes production behaviour. Every one adds tests, and Tasks 9, 10 refactor without changing outcomes.

### Task 8: `cross_match.py` — the largest untested wrong-answer surface

**Lane:** 8 (owns `vntyper/scripts/cross_match.py`, creates `tests/unit/test_cross_match.py`)

**Files:**
- Create: `tests/unit/test_cross_match.py`
- Test: `tests/unit/test_cross_match.py`

**Interfaces:**
- Consumes: `determine_variant_type(ref, alt)`, `compute_allele_change(ref, alt, variant_type)`, `cross_match_variants(kestrel_records, advntr_records, config=None)`, `write_results_tsv(results, output_path)`, `extract_results_from_pipeline_summary(summary)` — all in `vntyper/scripts/cross_match.py`.
- Produces: nothing importable.

**Why this task exists.** 190 LOC at **13.9% line / 9.4% branch-inclusive coverage with no test file at all** — 96 uncovered units, the worst ratio in the repo. It `eval()`s a `match_logic` rule string from config at `cross_match.py:137`, which is AGENTS.md trap 3 — the class of defect that silently disabled a flag for months. It also produces the `"Cross-Match Variant Comparison"` step matched by exact string in `generate_report.py` and `cohort_summary.py` (trap 5).

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_cross_match.py` with, at minimum:

```python
"""Unit tests for cross_match.py.

This module had no test file. It evaluates a rule string from config against
DataFrame-derived names (AGENTS.md trap 3) and emits a summary step matched by
exact string literal downstream (trap 5). Both are silent-failure surfaces.
"""

import json

import pytest

from vntyper.scripts.cross_match import (
    compute_allele_change,
    cross_match_variants,
    determine_variant_type,
    extract_results_from_pipeline_summary,
    write_results_tsv,
)

pytestmark = pytest.mark.unit


@pytest.mark.parametrize(
    ("ref", "alt", "expected"),
    [
        ("C", "CC", "Insertion"),
        ("C", "CGGCA", "Insertion"),
        ("CC", "C", "Deletion"),
        # Equal length returns "Other" -- NOT "Substitution". cross_match.py:48-49.
        ("C", "G", "Other"),
        ("C", "C", "Other"),
        # Non-string inputs are coerced at :42-43, so a numeric REF must not crash.
        (1, 22, "Insertion"),
    ],
)
def test_determine_variant_type(ref, alt, expected):
    assert determine_variant_type(ref, alt) == expected


def test_compute_allele_change_strips_the_common_prefix_for_an_insertion():
    """cross_match.py:70-72 -- the return is the suffix, with no sign character."""
    assert compute_allele_change("C", "CGG", "Insertion") == "GG"


def test_compute_allele_change_strips_the_common_prefix_for_a_deletion():
    """cross_match.py:74-76."""
    assert compute_allele_change("CGG", "C", "Deletion") == "GG"


def test_compute_allele_change_returns_the_whole_allele_when_there_is_no_common_prefix():
    """cross_match.py:72 and :76 -- the fallback branch, otherwise uncovered."""
    assert compute_allele_change("C", "GGA", "Insertion") == "GGA"
    assert compute_allele_change("GGA", "C", "Deletion") == "GGA"


def test_duplication_is_treated_as_an_insertion():
    """cross_match.py:69 lists "duplication" beside "insertion"."""
    assert compute_allele_change("C", "CC", "Duplication") == "C"


def test_compute_allele_change_returns_empty_for_an_unknown_type():
    """cross_match.py:77 -- the bare `return ""`."""
    assert compute_allele_change("C", "G", "Nonsense") == ""


def test_a_matching_pair_is_reported_as_a_match():
    """overall_match is the STRING "Yes", not a boolean. cross_match.py:146."""
    result = cross_match_variants(
        kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        advntr_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
    )
    assert result["overall_match"] == "Yes"
    assert result["matches"][0]["Match"] == "Yes"


def test_a_non_matching_pair_is_not_reported_as_a_match():
    result = cross_match_variants(
        kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        advntr_records=[{"REF": "C", "ALT": "CGGCA", "POS": 67}],
    )
    assert result["overall_match"] == "No"
    assert result["matches"][0]["Match"] == "No"


def test_an_empty_advntr_record_set_produces_no_comparisons():
    """The nested loop at :121-122 never runs, so `matches` is empty."""
    result = cross_match_variants(kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}], advntr_records=[])
    assert result["overall_match"] == "No"
    assert result["matches"] == []


def test_an_explicit_kestrel_variant_field_wins_over_the_inferred_type():
    """cross_match.py:110 -- `Variant` is used when present and non-blank."""
    result = cross_match_variants(
        kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67, "Variant": "Duplication"}],
        advntr_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
    )
    assert result["matches"][0]["Kestrel_Variant_Type"] == "Duplication"


def test_a_rule_naming_a_column_that_does_not_exist_is_silently_no_match_today(caplog):
    """CHARACTERISATION of a live trap-3 fail-open. Do not "fix" this here.

    The eval at cross_match.py:137 is wrapped in `except Exception` at :138-140,
    which sets `match = False` and logs at ERROR. A rule naming a column the
    records do not carry therefore reports "no match" rather than failing --
    the same shape as the `Poylmorhic_Call` typo that disabled a flag for
    months, and the `RU == 7` rule that never fired in its life.

    Turning it fail-closed changes what the cross-match step reports, so it
    needs a domain decision. Filed; pinned here so it cannot drift silently.
    """
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.cross_match")

    result = cross_match_variants(
        kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        advntr_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        config={"cross_match": {"match_logic": "Nonexistent_Column == 1"}},
    )

    assert result["overall_match"] == "No"
    assert any("Error evaluating match logic" in r.message for r in caplog.records)


def test_the_default_match_logic_is_used_when_no_config_is_given():
    """cross_match.py:100-103 -- both arms of the config branch."""
    with_none = cross_match_variants(
        kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        advntr_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        config=None,
    )
    with_empty = cross_match_variants(
        kestrel_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        advntr_records=[{"REF": "C", "ALT": "CC", "POS": 67}],
        config={},
    )
    assert with_none["overall_match"] == with_empty["overall_match"] == "Yes"
```

Then read `write_results_tsv` (`:150`) and `extract_results_from_pipeline_summary`
(`:169`) and write the same shape of test for each, including the empty-input path and
the `"Kestrel Genotyping"` / `"adVNTR Genotyping"` step literals (trap 5 — import them
from `vntyper/scripts/summary_steps.py`, never as bare strings).

**These are characterisation tests of an untested module. Where a test disagrees with the code, the code wins** — unless you can show it is a defect, in which case pin the current behaviour, name the test `..._today`, say so in the docstring, and file an issue.

- [ ] **Step 2: Run them and record which fail**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/test_cross_match.py -v
```

Expected: the tests above were written against the source and should pass. Any that fail reveal a behaviour the source has and this plan got wrong — **fix the test, not the module.** Record which ones, because that list is itself a finding about how little was known about this module.

Then prove each test can fail. Induce the defect it guards — flip `>` to `<` at `cross_match.py:44`, change `"Other"` to `"Substitution"` at `:49`, drop the `alt.startswith(ref)` guard at `:70` — and confirm the matching test fails. Paste each induced failure. Revert every change; `git diff --quiet -- vntyper/` must be clean.

- [ ] **Step 3: File the trap-3 fail-open, do not fix it**

`cross_match.py:138-140` swallows every `eval` failure into `match = False`. It is pinned above as characterisation. File it:

```bash
gh issue create --repo hassansaei/VNtyper \
  --title "cross_match match_logic fails open: an eval error reports 'no match' rather than failing" \
  --label bug
```

Body: quote `:135-141`, name the two precedents (`Poylmorhic_Call`, `RU == 7`), state that the fix is behaviour-changing because it alters what the cross-match step reports, and cite the pinning test by name. Report the issue number to the orchestrator so Task 24 can reference it.

- [ ] **Step 4: Measure the gain**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  COVERAGE_FILE=.coverage.cross_match \
  pytest -m unit tests/unit -q --cov=vntyper/scripts/cross_match.py --cov-branch --cov-report=term
```

Expected: `cross_match.py` well above 9.4%. Record the exact figure for Task 16.

- [ ] **Step 5: Run the full unfiltered tier, then hand off**

Orchestrator commit:

```bash
git add tests/unit/test_cross_match.py
git commit -m "test(cross-match): cover the untested cross-match module"
```

---

### Task 9: Merge the twin VCF preprocessors and pin their column contract

**Lane:** 9 (owns `vntyper/scripts/motif_processing.py` for phase 2 — Tasks 9 and 10 run **serially** in this lane)

**Files:**
- Modify: `vntyper/scripts/motif_processing.py:60-111`
- Create: `tests/unit/test_motif_preprocessing.py`
- Test: `tests/unit/test_motif_preprocessing.py`

**Interfaces:**
- Consumes: the guard added by Task 4 — **preserve it**.
- Produces: `_preprocess_vcf_frame(df, muc1_ref, variant_label) -> pd.DataFrame`. `preprocessing_insertion(df, muc1_ref)` and `preprocessing_deletion(df, muc1_ref)` keep their public signatures and delegate.

**Why.** `preprocessing_insertion` (60–84) and `preprocessing_deletion` (87–111) are identical except the `Variant` literal, and they hold **12 of the module's 46 surviving mutants** — lines 77, 79×2, 80×2, 81 and 104, 106×2, 107×2, 108. Every one is column plumbing: the `#CHROM`→`Motifs` rename, the five-element drop list, the last-column→`Sample` rename. They survive because the existing tests assert on merged output, never on the column contract.

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_motif_preprocessing.py`:

```python
"""The VCF column contract both preprocessors depend on.

12 of motif_processing.py's 46 surviving mutants live in these two functions
and every one is column plumbing that no test observed: the #CHROM rename,
the drop list, and the positional last-column-is-Sample assumption. Each test
below kills at least one of them.
"""

import pandas as pd
import pytest

from vntyper.scripts.motif_processing import preprocessing_deletion, preprocessing_insertion

pytestmark = pytest.mark.unit


def _vcf_frame():
    """A frame shaped like a filtered Kestrel VCF, with the sample column last."""
    return pd.DataFrame(
        {
            "#CHROM": ["X-5"],
            "POS": [67],
            "ID": ["."],
            "REF": ["C"],
            "ALT": ["CC"],
            "QUAL": ["."],
            "FILTER": ["PASS"],
            "INFO": ["."],
            "FORMAT": ["GT:AD"],
            "SAMPLE_1": ["Ins:50:5000"],
        }
    )


def _muc1_ref():
    return pd.DataFrame({"Motifs": ["X-5"], "Motif_sequence": ["ACGTACGT"]})


@pytest.mark.parametrize(
    ("func", "expected_label"),
    [(preprocessing_insertion, "Insertion"), (preprocessing_deletion, "Deletion")],
)
def test_chrom_is_renamed_to_motifs(func, expected_label):
    """Kills the line 77 / 104 rename mutants."""
    out = func(_vcf_frame(), _muc1_ref())
    assert "Motifs" in out.columns
    assert "#CHROM" not in out.columns
    assert out["Variant"].iloc[0] == expected_label


@pytest.mark.parametrize("func", [preprocessing_insertion, preprocessing_deletion])
def test_exactly_the_five_vcf_metadata_columns_are_dropped(func):
    """Kills the line 78-79 / 105-106 drop-list mutants."""
    out = func(_vcf_frame(), _muc1_ref())
    for dropped in ("ID", "QUAL", "FILTER", "INFO", "FORMAT"):
        assert dropped not in out.columns
    for kept in ("Motifs", "POS", "REF", "ALT"):
        assert kept in out.columns


@pytest.mark.parametrize("func", [preprocessing_insertion, preprocessing_deletion])
def test_the_last_column_becomes_sample_and_it_is_the_last_one(func):
    """Kills the line 80-81 / 107-108 positional mutants.

    ``df.columns[-1]`` is a positional assumption: the sample column is
    whatever ends up last after the drop. A `-1` -> `-2` or `+1` mutation
    renames the wrong column, and nothing downstream re-checks it.
    """
    frame = _vcf_frame()
    out = func(frame, _muc1_ref())
    assert "Sample" in out.columns
    assert out["Sample"].iloc[0] == "Ins:50:5000"
    assert "SAMPLE_1" not in out.columns


@pytest.mark.parametrize("func", [preprocessing_insertion, preprocessing_deletion])
def test_a_second_sample_column_means_only_the_last_is_renamed(func):
    """A multi-sample VCF: the positional rule must take the final column."""
    frame = _vcf_frame()
    frame.insert(len(frame.columns) - 1, "SAMPLE_0", ["Ins:1:2"])
    out = func(frame, _muc1_ref())
    assert out["Sample"].iloc[0] == "Ins:50:5000"
    assert "SAMPLE_0" in out.columns


@pytest.mark.parametrize("func", [preprocessing_insertion, preprocessing_deletion])
def test_the_merge_key_is_motifs_and_an_unmatched_motif_yields_na(func):
    out = func(_vcf_frame(), pd.DataFrame({"Motifs": ["NOPE"], "Motif_sequence": ["TTTT"]}))
    assert out["Motif_sequence"].isna().all()
```

- [ ] **Step 2: Run them and confirm they pass, then prove each can fail**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/test_motif_preprocessing.py -v
```

Expected: all pass. Then, one at a time, apply the exact surviving mutation and confirm a failure:

| Line | Mutation | Test that must fail |
| --- | --- | --- |
| 77 | `True` → `False` (the `inplace=`) | `test_chrom_is_renamed_to_motifs` |
| 79 | `1` → `2` (the `axis=`) | `test_exactly_the_five_vcf_metadata_columns_are_dropped` |
| 80 | `-1` → `+1` | `test_the_last_column_becomes_sample_and_it_is_the_last_one` |
| 81 | `True` → `False` | `test_the_last_column_becomes_sample_and_it_is_the_last_one` |

Paste each induced failure. **Revert every mutation before continuing** — `git diff --quiet -- vntyper/` must be clean of anything you did not intend.

- [ ] **Step 3: Extract the shared function**

Replace lines 60–111 with:

```python
#: VCF metadata columns dropped before the motif merge. Both preprocessors
#: drop exactly this set; a change here changes both.
VCF_METADATA_COLUMNS = ["ID", "QUAL", "FILTER", "INFO", "FORMAT"]


def _preprocess_vcf_frame(df, muc1_ref, variant_label):
    """Shared body of preprocessing_insertion and preprocessing_deletion.

    The two differed only in the ``Variant`` literal. Keeping one body means a
    fix to the column plumbing cannot land on one path and miss the other.

    Args:
        df (pd.DataFrame): Variants from a filtered VCF, sample column last.
        muc1_ref (pd.DataFrame): MUC1 reference with 'Motifs' & 'Motif_sequence'.
        variant_label (str): Value written to the 'Variant' column.

    Returns:
        pd.DataFrame: Merged frame with 'Motifs', 'Sample' and 'Variant'.
    """
    df.rename(columns={"#CHROM": "Motifs"}, inplace=True)
    df.drop(columns=VCF_METADATA_COLUMNS, axis=1, inplace=True)
    last_column_name = df.columns[-1]
    df.rename(columns={last_column_name: "Sample"}, inplace=True)
    df = pd.merge(df, muc1_ref, on="Motifs", how="left")
    df["Variant"] = variant_label
    return df


def preprocessing_insertion(df, muc1_ref):
    """Preprocess insertion variants by merging them with the MUC1 reference motifs.

    Args:
        df (pd.DataFrame): Insertion variants from a filtered VCF.
        muc1_ref (pd.DataFrame): MUC1 reference DataFrame.

    Returns:
        pd.DataFrame: Updated with 'Variant' = 'Insertion'.
    """
    return _preprocess_vcf_frame(df, muc1_ref, "Insertion")


def preprocessing_deletion(df, muc1_ref):
    """Preprocess deletion variants by merging them with the MUC1 reference motifs.

    Args:
        df (pd.DataFrame): Deletion variants from a filtered VCF.
        muc1_ref (pd.DataFrame): MUC1 reference DataFrame.

    Returns:
        pd.DataFrame: Updated with 'Variant' = 'Deletion'.
    """
    return _preprocess_vcf_frame(df, muc1_ref, "Deletion")
```

- [ ] **Step 4: Run the full unfiltered tier**

Expected: all pass. The two public functions kept their signatures, so `kestrel_genotyping.py:434-435` is unaffected.

- [ ] **Step 5: Hand off**

Orchestrator commit:

```bash
git add vntyper/scripts/motif_processing.py tests/unit/test_motif_preprocessing.py
git commit -m "test(motif): pin the VCF column contract and merge the twin preprocessors"
```

---

### Task 10: Extract the pure motif decision layer

**Lane:** 9 (runs **after** Task 9 — same file)

**Files:**
- Create: `vntyper/scripts/motif_decisions.py`
- Create: `tests/unit/test_motif_decisions.py`
- Modify: `vntyper/scripts/motif_processing.py:282-463`
- Test: `tests/unit/test_motif_decisions.py`

**Interfaces:**
- Consumes: `_preprocess_vcf_frame` from Task 9 (unused here, but same file — do not revert it). The Task 4 guard must survive.
- Produces:
  - `split_left_right(df, position_threshold) -> tuple[pd.DataFrame, pd.DataFrame]`
  - `apply_right_motif_exclusions(motif_right, exclude_motifs_right) -> pd.DataFrame`
  - `apply_gg_alt_rule(motif_right, alt_for_motif_right_gg, motifs_for_alt_gg) -> pd.DataFrame`
  - `apply_combined_exclusions(df, exclude_alts_combined, exclude_motifs_combined) -> pd.DataFrame`

**Why.** `motif_correction_and_annotation` is 182 lines holding **28 of the 46 surviving mutants**, because it fuses config read, column plumbing, the positional split and the motif filtering decisions into one function that cannot be called without a full frame. AGENTS.md rule 3: extract the pure part, leave the I/O.

**Binding constraint.** This is a **pure refactor**. Every extracted function must produce byte-identical output to the code it replaces. It is not authorised to change a motif decision — #186 reserves that. If you find a defect while extracting, pin it with a characterisation test and file an issue; do not fix it.

**This task is GATED.** It moved out of phase 2 after the adversarial review: it extracts the code that computes `motif_filter_pass` (`motif_processing.py:430`), which `filter_final_dataframe` requires as a gate column (`kestrel_genotyping.py:806`). "Intended to be behaviour-preserving" is not evidence — a unit oracle over one synthetic frame is far weaker than 58 real samples. **Run Task 21 against this commit before moving on**, and do not batch it with any other gated commit.

- [ ] **Step 1: Pin current behaviour end to end before touching anything**

Write a characterisation test that drives `motif_correction_and_annotation` with the shipped config and a frame covering: a left-motif row (`POS < 60`), a right-motif row (`POS >= 60`), an excluded right motif (one of `8 9 7 6p 6`), a GG alt, and a combined-exclusion ALT (one of `CCGCC CGGCG CGGCC`). Assert the exact resulting row set and the exact `motif_filter_pass` values. Run it, record the output as the oracle.

- [ ] **Step 2: Run it and confirm it passes**

Expected: PASS. This is the invariant Steps 3–6 must preserve.

- [ ] **Step 3: Write the failing tests for the extracted functions**

Create `tests/unit/test_motif_decisions.py` with table tests for each of the four functions above, covering at minimum the mutants at lines 340 (`<` → `>=`), 344/439/441/443 (`not` deleted), 346/449/452/455 (`in` → `not in`), 374 (`in` → `not in`), 424 (`>=` → `<`). One test per mutant, named for what it protects.

- [ ] **Step 4: Run them and confirm they fail**

Expected: `ModuleNotFoundError: No module named 'vntyper.scripts.motif_decisions'`.

- [ ] **Step 5: Create the module by moving code, not rewriting it**

Cut each decision block out of `motif_correction_and_annotation` into `motif_decisions.py` verbatim, parameterising only what it read from `mf`. Fully annotate the new module — `scoring.py`, `region_utils.py` and `flagging.py` are the shape to copy.

- [ ] **Step 6: Run the Step 1 oracle and the new tests together**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/test_motif_decisions.py tests/unit/test_motif_preprocessing.py \
         tests/unit/test_motif_filtering_issue_136.py tests/unit/test_kestrel_filtering.py -v
```

Expected: all pass, the Step 1 oracle included and unmodified. **If the oracle needs editing, the refactor changed behaviour — stop and report.**

- [ ] **Step 6b: Register the new module with the mutation instrument — do this BEFORE the sweep**

`scripts/mutation_test.py`'s `TARGETS` dict (line ~96) names `vntyper/scripts/motif_processing.py` but would not name `motif_decisions.py`. **A module absent from `TARGETS` is not mutated at all.** Extracting the decision logic without registering it would make `motif_processing.py`'s score rise for the worst possible reason: the hard decisions left the measurement.

```python
    "vntyper/scripts/motif_decisions.py": (
        "tests/unit/test_motif_decisions.py",
        "tests/unit/test_motif_filtering_issue_136.py",
        "tests/unit/test_kestrel_filtering.py",
    ),
```

Keep the test list wide, per that file's own comment — every file that imports the module belongs there. Update `tests/unit/test_mutation_test.py` if it pins the `TARGETS` keys.

**Report the two modules' scores together.** A `motif_processing.py` figure quoted alone after this extraction is not comparable to the 30.9% baseline and must not be presented as if it were.

- [ ] **Step 6c: Add the line-332 characterisation test**

Task 24 files `motif_processing.py:332` as a defect and says it is pinned by a test added here. It is not, unless you add it — the Step 1 oracle has no malformed motif ID and none of the mutant tests reaches that line.

```python
def test_one_malformed_motif_id_nulls_every_row_in_the_sample_today():
    """CHARACTERISATION of a live defect. Do not fix it here.

    motif_processing.py:332 aggregates with `.max()` over the WHOLE column:

        working_df["Motifs"].str.count("-").max() != 1

    so a single row whose Motifs value lacks exactly one dash makes the
    condition true for every row, and the entire sample gets
    motif_filter_pass=False. One bad motif ID suppresses the whole call, and
    the pipeline exits 0. This is W8 from the #179 audit; it was never filed.
    The `or` -> `and` mutant at this line also survives the whole unit tier.

    Fixing it changes which rows survive motif filtering, which is a
    reported-genotype change with no decision behind it. Filed, not fixed.
    """
    df = pd.DataFrame(
        {
            "Motifs": ["X-5", "BROKEN"],  # the second has zero dashes
            "POS": [67, 68],
            "REF": ["C", "C"],
            "ALT": ["CC", "CC"],
            "Estimated_Depth_AlternateVariant": [50, 50],
            "Estimated_Depth_Variant_ActiveRegion": [5000, 5000],
            "Depth_Score": [0.01, 0.01],
            "Confidence": ["High_Precision", "High_Precision"],
        }
    )

    out = motif_correction_and_annotation(df, _merged_motifs(), _shipped_config())

    assert not out["motif_filter_pass"].any(), "the valid row is suppressed by the invalid one"
```

Run it, confirm it passes, then prove it can fail by making the check per-row — and **revert that immediately**.

- [ ] **Step 7: Run the full unfiltered tier and the scoped mutation sweep**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest -m unit tests/unit -q
```

Then, **orchestrator-only and with nothing else running**:

```bash
git diff --quiet -- vntyper/ && echo "tree clean, safe to sweep"
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  make mutation
git diff --quiet -- vntyper/ && echo "tree restored"
```

Expected: `motif_processing.py` well above 30.9%. `PYTHONDONTWRITEBYTECODE=1` is in the recipe and is load-bearing — do not work around it.

- [ ] **Step 8: Hand off**

Orchestrator commits:

```bash
git add vntyper/scripts/motif_decisions.py vntyper/scripts/motif_processing.py tests/unit/test_motif_decisions.py
git commit -m "refactor(motif): extract the pure decision layer from motif_correction_and_annotation

Behaviour-preserving: the pre-existing end-to-end oracle passes unmodified."

git add docs/development/mutation-testing.md docs/development/mutation-results.json
git commit -m "docs(mutation): re-sweep after the motif_processing tests"
```

---

### Task 11: `utils.py` — 147 uncovered units

**Lane:** 11 (owns `vntyper/scripts/utils.py`, `tests/unit/test_utils.py`)

**Files:**
- Modify: `tests/unit/test_utils.py`
- Test: `tests/unit/test_utils.py`

**Interfaces:**
- Consumes: whatever `utils.py` exports. **Do not change `run_command`'s `shell=True`** — it is deliberate, for process substitution in the CRAM unmapped-read path.
- Produces: nothing importable.

`utils.py` is 377 LOC at 41.0% line / 38.8% branch-inclusive — the largest gap in a file under the 650-LOC limit. Uncovered regions per the branch report: lines 134-136, 153-199, 213-235, 250-260, 280-296, 317-318, 326-327, 336-337, 352-377, plus the partial exits at 89→92 and 95→102.

- [ ] **Step 1: List what is uncovered**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  COVERAGE_FILE=.coverage.utils \
  pytest -m unit tests/unit -q --cov=vntyper/scripts/utils.py --cov-branch --cov-report=term-missing
```

Record the missing-line list. Read each region and note which are pure (testable with `tmp_path` + `unittest.mock`) and which need a subprocess boundary mocked.

- [ ] **Step 2: Write failing tests for the pure regions first**

One test per behaviour, each with a docstring naming what it protects. Mock at the `subprocess` boundary with `unittest.mock.patch`; never invoke a real binary. Include the **failure** path of every function you touch — a `run_command` that exits non-zero, a missing file, a malformed header.

Note: `run_command` logs failures at `logger.debug`, so a failure produces nothing at default INFO. If you assert on a log, pin the level and the logger name — `pytest.ini` sets `log_level = DEBUG`, so a naive `caplog` assertion passes against invisible logs.

- [ ] **Step 3: Run them and confirm they fail**

Expected: `AttributeError`/`AssertionError` per test, before implementation. Where you are covering existing code, prove each test can fail by inducing the defect it guards, and paste the output.

- [ ] **Step 4: Run them and confirm they pass**

- [ ] **Step 5: Re-measure**

Same command as Step 1. Target: `utils.py` above 75% branch-inclusive, i.e. at least 90 of its 147 missing units recovered. Record the figure.

- [ ] **Step 6: Run the full unfiltered tier, then hand off**

```bash
git add tests/unit/test_utils.py
git commit -m "test(utils): cover the command, logging and validation paths"
```

---

### Task 12: `file_processing.py` — 10.5% covered

**Lane:** 12 (owns `vntyper/scripts/file_processing.py`, creates `tests/unit/test_file_processing.py`)

**Files:**
- Create: `tests/unit/test_file_processing.py`
- Test: `tests/unit/test_file_processing.py`

**Interfaces:**
- Consumes: whatever `file_processing.py` exports.
- Produces: nothing importable.

38 statements, 34 missing units, 10.5% branch-inclusive. Two functions, `filter_vcf` (`:7`) and `filter_indel_vcf` (`:34`), both plain file I/O over VCF text and both trivially testable with `tmp_path`. They decide which variants ever reach Kestrel post-processing, so a defect here is a silently missing call.

- [ ] **Step 1: Write the failing tests**

```python
"""Unit tests for file_processing.py -- the indel split that feeds Kestrel.

`filter_vcf` decides what counts as an indel and `filter_indel_vcf` decides
which side of the insertion/deletion split it lands on. Both were untested.
A row dropped here never reaches motif processing, scoring or the report.
"""

import pytest

from vntyper.scripts.file_processing import filter_indel_vcf, filter_vcf

pytestmark = pytest.mark.unit

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _vcf(tmp_path, name, *rows):
    path = tmp_path / name
    path.write_text(HEADER + "".join(f"X-5\t{i}\t.\t{ref}\t{alt}\t.\tPASS\t.\n" for i, (ref, alt) in enumerate(rows)))
    return path


def _rows(path):
    return [line for line in path.read_text().splitlines() if not line.startswith(("##", "#CHROM"))]


def test_the_header_is_carried_through_verbatim(tmp_path):
    out = tmp_path / "out.vcf"
    filter_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "CC"))), str(out))
    assert out.read_text().startswith("##fileformat=VCFv4.2\n#CHROM")


def test_a_one_base_to_many_insertion_is_kept(tmp_path):
    out = tmp_path / "out.vcf"
    filter_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "CGGCA"))), str(out))
    assert len(_rows(out)) == 1


def test_a_many_to_one_base_deletion_is_kept(tmp_path):
    out = tmp_path / "out.vcf"
    filter_vcf(str(_vcf(tmp_path, "in.vcf", ("CGGCA", "C"))), str(out))
    assert len(_rows(out)) == 1


def test_a_single_base_substitution_is_dropped(tmp_path):
    out = tmp_path / "out.vcf"
    filter_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "G"))), str(out))
    assert _rows(out) == []


def test_an_indel_where_neither_allele_is_one_base_is_dropped_today(tmp_path):
    """CHARACTERISATION of a real gap, not a specification.

    file_processing.py:28-30 keeps a row only when EXACTLY ONE of REF/ALT has
    length snv_length (1). REF="AC", ALT="ACGGG" is a genuine 3-base insertion
    and is silently discarded before Kestrel post-processing ever sees it.
    Whether Kestrel can emit such a record is a domain question; pinned here so
    the behaviour cannot drift, and filed rather than fixed.
    """
    out = tmp_path / "out.vcf"
    filter_vcf(str(_vcf(tmp_path, "in.vcf", ("AC", "ACGGG"))), str(out))
    assert _rows(out) == []


def test_a_line_with_too_few_columns_raises(tmp_path):
    """file_processing.py:27 unpacks six fields; a short line is not tolerated."""
    path = tmp_path / "short.vcf"
    path.write_text(HEADER + "X-5\t1\t.\tC\n")
    with pytest.raises(ValueError):
        filter_vcf(str(path), str(tmp_path / "out.vcf"))


def test_insertions_and_deletions_are_written_to_separate_files(tmp_path):
    ins, dele = tmp_path / "ins.vcf", tmp_path / "del.vcf"
    filter_indel_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "CGG"), ("CGG", "C"))), str(ins), str(dele))
    assert len(_rows(ins)) == 1
    assert len(_rows(dele)) == 1


def test_both_output_files_get_the_header(tmp_path):
    ins, dele = tmp_path / "ins.vcf", tmp_path / "del.vcf"
    filter_indel_vcf(str(_vcf(tmp_path, "in.vcf", ("C", "CGG"))), str(ins), str(dele))
    assert dele.read_text().startswith("##fileformat")
    assert _rows(dele) == []


def test_a_multi_base_ref_insertion_is_classified_as_a_deletion_today(tmp_path):
    """CHARACTERISATION of a real mis-classification, not a specification.

    file_processing.py:60-63 tests only `len(ref) == 1 and len(alt) > 1` for
    the insertion side; the `else` is an unconditional catch-all. REF="AC",
    ALT="ACGGG" is an insertion and lands in the DELETION file, where the
    deletion frame rule (3n+2) is applied to it instead of the insertion rule
    (3n+1) -- the two are not interchangeable (#181). It cannot arrive today
    because `filter_vcf` drops it first (see above), so the two gaps mask each
    other. Pinned and filed; not fixed.
    """
    ins, dele = tmp_path / "ins.vcf", tmp_path / "del.vcf"
    path = tmp_path / "indel.vcf"
    path.write_text(HEADER + "X-5\t1\t.\tAC\tACGGG\t.\tPASS\t.\n")
    filter_indel_vcf(str(path), str(ins), str(dele))
    assert _rows(ins) == []
    assert len(_rows(dele)) == 1
```

- [ ] **Step 2: Run them and confirm they fail**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/test_file_processing.py -v
```

Expected: all pass, since they characterise existing code. **Prove each can fail**: flip `!=` to `==` at `:28`, change `>` to `>=` at `:60`, and drop the header branch at `:24`. Each must break a named test. Paste the failures and revert.

- [ ] **Step 3: File the two characterised gaps**

Both `..._today` tests above describe real defects that mask each other. File one issue covering both, quoting `:28-30` and `:60-63`, noting that the insertion/deletion misrouting matters because #181 makes the two frame rules non-interchangeable. Label `needs-domain-decision`. **Do not fix either** — changing what reaches Kestrel changes reported genotypes.

- [ ] **Step 4: Confirm the full unfiltered tier is green**

- [ ] **Step 5: Re-measure**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  COVERAGE_FILE=.coverage.file_processing \
  pytest -m unit tests/unit -q --cov=vntyper/scripts/file_processing.py --cov-branch --cov-report=term-missing
```

Target: above 85% branch-inclusive. Record the figure.

- [ ] **Step 6: Run the full unfiltered tier, then hand off**

```bash
git add tests/unit/test_file_processing.py
git commit -m "test(file-processing): cover the module end to end"
```

---

### Task 13: `extract_unmapped_from_offset.py` — 17.0% covered

**Lane:** 13 (owns `vntyper/scripts/extract_unmapped_from_offset.py`, creates `tests/unit/test_extract_unmapped.py`)

**Files:**
- Create: `tests/unit/test_extract_unmapped.py`
- Test: `tests/unit/test_extract_unmapped.py`

**Interfaces:**
- Consumes: whatever the module exports.
- Produces: nothing importable.

47 statements, 39 missing units. Four functions: `read_uint32` (`:16`), `read_uint64` (`:23`), `get_last_chunk_end` (`:30`) and `extract_unmapped_reads_from_offset` (`:56`). The first three are pure binary parsing over a file object and need no BAM — synthesise the BAI bytes with `struct.pack`. Only the fourth needs `pysam`, and it can be mocked.

This module decides where the unmapped-read scan starts. A parse error here means reads are silently missed, and the genotype is computed on a fraction of the input.

- [ ] **Step 1: Write the failing tests**

```python
"""Unit tests for extract_unmapped_from_offset.py.

get_last_chunk_end walks a BAI index by hand to find where the mapped chunks
end, and the unmapped scan starts there. If it returns too high an offset,
reads are silently skipped and the genotype is computed on a fraction of the
input -- with no error anywhere.
"""

import io
import struct

import pytest

from vntyper.scripts.extract_unmapped_from_offset import get_last_chunk_end, read_uint32, read_uint64

pytestmark = pytest.mark.unit


def _bai(references):
    """Build a minimal BAI. `references` is a list of lists of (beg, end) chunks."""
    out = [b"BAI\x01", struct.pack("<I", len(references))]
    for chunks in references:
        out.append(struct.pack("<I", 1 if chunks else 0))  # n_bins
        if chunks:
            out.append(struct.pack("<I", 4681))            # bin number, unused
            out.append(struct.pack("<I", len(chunks)))     # n_chunks
            for beg, end in chunks:
                out.append(struct.pack("<QQ", beg, end))
        out.append(struct.pack("<I", 2))                   # n_intv
        out.append(struct.pack("<QQ", 0, 0))               # the linear index, skipped
    return b"".join(out)


def test_read_uint32_is_little_endian_and_unsigned():
    assert read_uint32(io.BytesIO(b"\x01\x00\x00\x00")) == 1
    assert read_uint32(io.BytesIO(b"\xff\xff\xff\xff")) == 4294967295


def test_read_uint64_is_little_endian_and_unsigned():
    assert read_uint64(io.BytesIO(b"\x01" + b"\x00" * 7)) == 1
    assert read_uint64(io.BytesIO(b"\xff" * 8)) == 18446744073709551615


def test_the_maximum_chunk_end_is_returned_not_the_last_one(tmp_path):
    """Chunks are not sorted by end offset; taking the last would truncate."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 900), (100, 500)]]))
    assert get_last_chunk_end(str(path)) == 900


def test_the_maximum_is_taken_across_every_reference(tmp_path):
    """The linear index of reference 1 must be skipped correctly to reach 2."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 100)], [(0, 7000)], [(0, 300)]]))
    assert get_last_chunk_end(str(path)) == 7000


def test_a_reference_with_no_bins_contributes_nothing(tmp_path):
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[], [(0, 42)]]))
    assert get_last_chunk_end(str(path)) == 42


def test_an_index_with_no_chunks_at_all_returns_zero(tmp_path):
    """max_vo starts at 0 (line 35); an empty index must not raise."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[]]))
    assert get_last_chunk_end(str(path)) == 0


def test_a_truncated_index_does_not_silently_return_a_low_offset(tmp_path):
    """CHARACTERISATION. `f.read(4)` on EOF returns b"" and int.from_bytes
    gives 0, so a truncated BAI parses as "no more references" rather than
    raising -- and the caller then scans from offset 0. Pinned so the
    behaviour cannot drift; whether it should raise is a separate question.
    """
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 900)]])[:12])
    assert get_last_chunk_end(str(path)) == 0
```

Then cover `extract_unmapped_reads_from_offset` (`:56`), which holds most of the remaining branches. Prose is not enough here — this is the I/O boundary:

```python
def test_the_scan_seeks_to_the_offset_the_index_reported(tmp_path, monkeypatch):
    """The whole point of the module: start where the mapped chunks ended."""
    path = tmp_path / "x.bai"
    path.write_bytes(_bai([[(0, 7000)]]))
    inbam, outbam = MagicMock(), MagicMock()
    inbam.__enter__.return_value = inbam
    outbam.__enter__.return_value = outbam
    inbam.fetch.return_value = iter([])
    monkeypatch.setattr(pysam, "AlignmentFile", MagicMock(side_effect=[inbam, outbam]))

    extract_unmapped_reads_from_offset(str(tmp_path / "in.bam"), str(path), str(tmp_path / "out.bam"))

    inbam.seek.assert_called_once_with(7000)


def test_only_unmapped_reads_are_written(tmp_path, monkeypatch):
    """A mapped read reaching the output would double-count it downstream."""
    mapped, unmapped = MagicMock(is_unmapped=False), MagicMock(is_unmapped=True)
    # ... same mock scaffolding, with the input yielding [mapped, unmapped]
    assert outbam.write.call_count == 1
    outbam.write.assert_called_once_with(unmapped)


def test_both_files_are_closed_when_the_scan_raises(tmp_path, monkeypatch):
    """`with` must cover the write side too, or a failed run leaks a handle
    and leaves a truncated BAM that looks complete."""
    # ... input iterator raises OSError mid-scan
    with pytest.raises(OSError):
        extract_unmapped_reads_from_offset(...)
    inbam.__exit__.assert_called_once()
    outbam.__exit__.assert_called_once()
```

Read `:69-80` first and match the real call shape — whether it uses `fetch(until_eof=True)`, an iterator, or a plain loop, and whether the output file is opened with a `template=` or a `header=`.

- [ ] **Step 2: Run them and confirm they pass, then prove each can fail**

Induce: change `max(max_vo, chunk_end)` at `:49` to `chunk_end`, change `read_uint64` to `read_uint32` at `:48`, and drop the `bai.seek` at `:52`. Each must break a named test. Paste the failures and revert.

- [ ] **Step 3: File the truncated-index finding if it is judged a defect**

The `..._does_not_silently_return_a_low_offset` test names a real fail-open. Report it to the orchestrator; do not fix it here.

- [ ] **Step 4: Run the full unfiltered tier**
- [ ] **Step 5: Re-measure with `COVERAGE_FILE=.coverage.extract_unmapped`. Target above 85% branch-inclusive**
- [ ] **Step 6: Hand off**

```bash
git add tests/unit/test_extract_unmapped.py
git commit -m "test(unmapped): cover the offset extractor"
```

---

### Task 14: `variant_parsing.py` — 3 surviving mutants and 29 uncovered units

**Lane:** 14 (owns `vntyper/scripts/variant_parsing.py`, `tests/unit/test_variant_parsing.py`)

**Files:**
- Modify: `tests/unit/test_variant_parsing.py`
- Test: `tests/unit/test_variant_parsing.py`

**Interfaces:**
- Consumes: whatever `variant_parsing.py` exports.
- Produces: nothing importable.

Three genuine survivors, all boolean-operator mutations that change which variants are filtered:

| Line | Mutation |
| --- | --- |
| 59 | `not` deleted |
| 59 | `and` → `or` |
| 68 | `and` → `or` |

Lines 50–72 are entirely uncovered.

- [ ] **Step 1: Read lines 50-72 and write down the truth table**

For each compound condition, enumerate the input combinations that distinguish `and` from `or`, and the one that distinguishes `not X` from `X`. There is at least one per mutant.

- [ ] **Step 2: Write one test per mutant, named for it**

```python
def test_both_conditions_are_required_not_either_of_them():
    """Kills the line 59 `and` -> `or` mutant.

    An input satisfying exactly one of the two conditions must NOT pass. With
    `or` it would, and no existing test distinguishes the two.
    """
```

- [ ] **Step 3: Run them and confirm they pass, then prove each can fail**

Apply each of the three mutations in turn and confirm the matching test fails. Paste each induced failure. Revert every mutation; `git diff --quiet -- vntyper/` must be clean.

- [ ] **Step 4: Cover lines 50-72**

- [ ] **Step 5: Re-measure with `COVERAGE_FILE=.coverage.variant_parsing`**

Target: above 90% branch-inclusive.

- [ ] **Step 6: Run the full unfiltered tier, then hand off**

```bash
git add tests/unit/test_variant_parsing.py
git commit -m "test(variant-parsing): kill the three surviving boolean mutants"
```

---

### Task 15: `docker/app/tasks.py` — 98 uncovered units

**Lane:** 15 (owns `docker/app/tasks.py` in phase 2 and `tests/unit/web/test_tasks.py`)

**Files:**
- Create or modify: `tests/unit/web/test_tasks.py`
- Test: `tests/unit/web/`

**Interfaces:**
- Consumes: `build_vntyper_command(...)` if Task 6 extracted it. **Task 6 must land first** — same file.
- Produces: nothing importable.

63.0% line / 61.4% branch-inclusive, 98 missing units. Mock Celery and the subprocess boundary; do not start a worker.

- [ ] **Step 1: List what is uncovered**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  COVERAGE_FILE=.coverage.tasks \
  pytest -m unit tests/unit -q --cov=docker/app/tasks.py --cov-branch --cov-report=term-missing
```

- [ ] **Step 2: Write failing tests for the uncovered paths, failure paths included**

Cover at minimum: a job that raises, a job whose subprocess exits non-zero, the `index_path=None` compatibility default, and cleanup on both success and failure.

- [ ] **Step 3: Run them and confirm they fail**
- [ ] **Step 4: Run them and confirm they pass**
- [ ] **Step 5: Re-measure. Target: above 85% branch-inclusive**
- [ ] **Step 6: Run the full unfiltered tier, then hand off**

```bash
git add tests/unit/web/test_tasks.py
git commit -m "test(web): cover the Celery task failure and cleanup paths"
```

---

# Phase 3 — Enable branch coverage

### Task 16: Turn on `branch = true` and raise the floor (#196)

**Files:**
- Modify: `pyproject.toml` `[tool.coverage.run]` and `[tool.coverage.report]`
- Modify: `docs/development/mutation-testing.md` (the "branch coverage is measured but not enabled" section)
- Test: `tests/unit/test_coverage_gate.py`

**Interfaces:**
- Consumes: the coverage gained in phase 2.
- Produces: a branch-inclusive `fail_under`.

**#196 says splitting `cohort_summary.py` and `install_references.py` is "likely a prerequisite". It is not**, and phase 2 was chosen on that basis. Measured at `9104f64`: 4530 of 6700 units covered = 67.6119%, floor 70, so **160 units** are needed.

**Budget check, because Tasks 9 and 10 moved out of phase 2 into the gated phase.** The six phase-2 lanes still hold far more than 160 uncovered units between them:

| Task | Module | Missing units |
| --- | --- | ---: |
| 11 | `utils.py` | 147 |
| 8 | `cross_match.py` | 96 |
| 15 | `docker/app/tasks.py` | 98 |
| 13 | `extract_unmapped_from_offset.py` | 39 |
| 12 | `file_processing.py` | 34 |
| 14 | `variant_parsing.py` | 29 |
| | **available** | **443** |

Capturing 36% of that clears the gap. `motif_processing.py`'s 68 units arrive later, in phase 5, and are not counted here.

**Two numbers, not one — do not conflate them.** 160 units reaches exactly `4690/6700 = 70.00%`, which lets branch coverage be **enabled** against the existing floor of 70. It does **not** raise the floor: `scripts/coverage_gate.py:166` advises a higher integer only when `int(total) > floor`, so **227** units are needed before the floor can go to 71. Enabling is the goal of this task; raising is a bonus if the phase-2 lanes overshoot.

- [ ] **Step 1: Measure with branch coverage on, before changing the config**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  COVERAGE_FILE=.coverage.branchcheck \
  pytest -m unit tests/unit -o log_cli=false -q --cov --cov-branch --cov-report=term
```

Record the TOTAL. **If it is below 70, stop** — return to phase 2 and report the shortfall in units, not in percent.

- [ ] **Step 2: Enable branch coverage**

In `pyproject.toml` `[tool.coverage.run]`, beside `source` and `relative_files`:

```toml
# Branch coverage. An `if` that is entered but never taken is not "covered" for
# a codebase whose failure mode is silently taking the wrong branch. Enabled in
# #196 once the branch-inclusive total cleared the floor; the floor below is the
# branch-inclusive figure and was never lowered to admit this measurement.
branch = true
```

- [ ] **Step 3: Run the gate and read the figure it prints**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  make test-unit-cov
```

Expected: FAIL if the printed figure is below the current floor, otherwise the gate prints the exact line to paste. **Use the number `scripts/coverage_gate.py` prints — never the rounded `TOTAL` column**, which rounds 70.57 and 69.51 to the same integer and would make CI fail on the very run that produced it.

- [ ] **Step 4: Raise the floor to that exact figure**

Edit `fail_under` at `pyproject.toml:341` to the printed value. **Never lower it.**

- [ ] **Step 5: Re-run and confirm green**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  make test-unit-cov
```

Expected: `✓ Unit coverage complete`, no FAIL line.

- [ ] **Step 6: Update the artefact**

Rewrite the "Related: branch coverage is measured but not enabled" section of `docs/development/mutation-testing.md` to record that it is now enabled, at what figure, on what commit, and that the floor was raised rather than the measurement lowered. Update the test that pins those numbers.

- [ ] **Step 7: Run the whole CI gate locally**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  make ci-local
```

Expected: exit 0 — actionlint, format, lint, mypy ×2, tests + coverage, patch coverage, strict docs, and the from-scratch `uv` install path.

- [ ] **Step 8: Commit**

```bash
git add pyproject.toml docs/development/mutation-testing.md tests/unit/test_coverage_gate.py
git commit -m "chore(coverage): enable branch coverage and raise the floor to the printed figure

Closes #196"
```

---

# Phase 5 — Genotype-affecting changes (strictly serial, one gate run each)

**Execution order: Task 9, Task 10, Task 17, Task 18, Task 19, Task 6a.** Tasks 9, 10 and 6a are written earlier in this document; they execute here.

Nothing in this phase runs in parallel. Each task is one commit, and **each commit gets its own Task 21 run naming its own candidate SHA before the next task starts.** Do not batch them: a joint run cannot attribute a failure, and two changes can produce offsetting deltas that read as clean.

This phase is last because the gate attests one commit and nothing after it. **Nothing may land on the branch after the final run except Task 24's documentation commit, which is re-attested by the tip run in phase 6.**

### Task 17: #184 — mid-band Depth_Score demotion takes precedence

**Files:**
- Modify: `vntyper/scripts/confidence_assignment.py:128-154`
- Modify: `tests/unit/test_confidence_boundaries.py`
- Test: `tests/unit/test_confidence_boundaries.py`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing importable. **Task 20 was cut** — nothing else in this plan touches `confidence_assignment.py`, so this commit is the file's only change and its gate run attests it alone.

**Decision**, quoted from @hassansaei on #184:

> "No — do not remove this intent. Any variant with Depth_Score between 0.00469 and 0.00515 (inclusive) must be Low_Precision, even when alternate depth is ≥21 (or higher). That mid-range Depth_Score demotion from 1.3 is still the desired behaviour. […] Please keep / restore the rule so mid-band Depth_Score cannot be promoted by later conditions — either by making this rule take precedence, or by excluding Depth_Score ≤ 0.00515 from the High tiers."

**What the code actually does, and where it disagrees with him.** The assignments at lines 149–154 are last-wins and `cond6` (`low < DS < high` → `Low_Precision`) is the **last** one. So the *interior* of the band is already demoted at every alt depth, and `cond1` covers `DS == low` with nothing overwriting it. The single divergence is:

```
Depth_Score == 0.00515 exactly, alt_depth >= 21
  cond5 (alt in [21,100), DS >= high) -> High_Precision
  cond2 (alt >= 100,      DS >= high) -> High_Precision*
```

**His second suggested implementation is unsafe as written.** Changing `cond2`/`cond5` to `> high` *alone* leaves `alt_depth >= 100, DS == 0.00515` matching no condition, falling through to `NEGATIVE_LABEL` — a reported call becomes a non-call. Do not do that.

- [ ] **Step 1: Write the failing tests**

`0.00515` is exactly `103/20000`, so **every integer depth pair** reaching the boundary has `alt` a multiple of 103 and `region` the matching multiple of 20000. The minimum is `(103, 20000)`, no such pair has `alt` in `[21, 100)`, so `cond5` is unreachable at the boundary and the only label that moves is `High_Precision*`.

**That argument holds for integer depths only, and the plan must say so.** `confidence_assignment.py:105` uses `pd.to_numeric(...).fillna(0)` with **no integer cast**, and the existing boundary matrix in this file deliberately supplies *fractional* region depths — where `0.00515` is reachable at alt depths 21–99 and plain `High_Precision` moves too. Production depths come from splitting Kestrel's `Sample` field, which carries read counts, so they are integer-valued in practice; but that is an input property, not an enforced invariant.

Add a test that pins the invariant rather than assuming it, and scope every claim in the commit message to integer depths:

```python
def test_production_depths_arrive_as_whole_numbers():
    """The 103/20000 blast-radius argument depends on this and nothing enforces it.

    confidence_assignment.py:105 coerces with pd.to_numeric and never casts to
    int. Kestrel emits read counts, so the values are integral in practice --
    but if a fractional depth ever reaches this function, Depth_Score can hit
    0.00515 at alt depths 21-99 and plain High_Precision moves as well, which
    the #184 analysis does not cover.
    """
    frame = split_depth_and_calculate_frame_score(kestrel_stage_frame("raw"))
    for col in ("Estimated_Depth_AlternateVariant", "Estimated_Depth_Variant_ActiveRegion"):
        values = pd.to_numeric(frame[col], errors="coerce").dropna()
        assert (values == values.astype(int)).all(), f"{col} carries a fractional depth"
```

```python
#: Every integer (alt, region) pair whose quotient is exactly the shipped
#: high_threshold. 0.00515 == 103/20000, so alt is always a multiple of 103
#: and therefore always >= 103 > alt_mid_high (100). cond5 (alt in [21,100))
#: is arithmetically unreachable here; only cond2 fires.
EXACT_HIGH_THRESHOLD_PAIRS = [(103, 20000), (206, 40000), (515, 100000), (5150, 1000000)]


@pytest.mark.parametrize(("alt", "region"), EXACT_HIGH_THRESHOLD_PAIRS)
def test_the_top_of_the_mid_band_is_low_precision(alt, region):
    """Depth_Score == 0.00515 exactly must be Low_Precision (#184).

    Specification. @hassansaei: "Any variant with Depth_Score between 0.00469
    and 0.00515 (inclusive) must be Low_Precision, even when alternate depth
    is >= 21 (or higher)."

    cond6 already demoted the OPEN interval at every alt depth because it was
    applied last, so this boundary was the only divergence from his intent.
    Before this change cond2 promoted these rows to High_Precision*.
    """
    assert alt / region == 0.00515, "this pair does not reach the boundary in IEEE 754"

    df = pd.DataFrame(
        {
            "Estimated_Depth_AlternateVariant": [alt],
            "Estimated_Depth_Variant_ActiveRegion": [region],
        }
    )

    out = calculate_depth_score_and_assign_confidence(df, _shipped_config())

    assert out["Confidence"].iloc[0] == "Low_Precision"


def test_the_bottom_of_the_mid_band_is_low_precision():
    """DS == 0.00469 exactly. cond1 already handled this; pin it so it stays."""
    df = pd.DataFrame(
        {
            "Estimated_Depth_AlternateVariant": [469],
            "Estimated_Depth_Variant_ActiveRegion": [100000],
        }
    )
    assert 469 / 100000 == 0.00469

    out = calculate_depth_score_and_assign_confidence(df, _shipped_config())

    assert out["Confidence"].iloc[0] == "Low_Precision"


def test_one_ulp_above_the_mid_band_is_still_promoted():
    """The change must not swallow the High tiers -- only the boundary moves.

    np.nextafter gives the smallest representable step above the threshold, so
    this is the tightest possible proof that the new `inclusive="both"` band
    does not extend past its top edge.
    """
    region = 1000000
    alt = 5151  # 5151/1000000 = 0.005151 > 0.00515
    assert alt / region > np.nextafter(0.00515, 1)

    df = pd.DataFrame(
        {
            "Estimated_Depth_AlternateVariant": [alt],
            "Estimated_Depth_Variant_ActiveRegion": [region],
        }
    )

    out = calculate_depth_score_and_assign_confidence(df, _shipped_config())

    assert out["Confidence"].iloc[0] == "High_Precision*"


@pytest.mark.parametrize(("alt", "region"), EXACT_HIGH_THRESHOLD_PAIRS)
def test_no_row_at_the_boundary_falls_through_to_negative(alt, region):
    """Guards the trap in his second suggested implementation.

    Changing cond2/cond5 to `> high` ALONE would leave these rows matching no
    condition at all, so they would keep NEGATIVE_LABEL -- turning a reported
    call into a non-call, which is strictly worse than the bug. The chosen
    implementation (mid-band applied last, inclusive at both ends) cannot do
    that, and this test is what holds it to that.
    """
    df = pd.DataFrame(
        {
            "Estimated_Depth_AlternateVariant": [alt],
            "Estimated_Depth_Variant_ActiveRegion": [region],
        }
    )

    out = calculate_depth_score_and_assign_confidence(df, _shipped_config())

    assert out["Confidence"].iloc[0] != "Negative"


def test_a_boundary_demotion_changes_which_variant_is_reported():
    """The label is not the only thing that moves -- selection ranks on it.

    select_single_best_variant sorts by Confidence first
    (kestrel_genotyping.py:756, priority documented at :680), so demoting a
    boundary row from High_Precision* to Low_Precision can hand the reported
    genotype to a different candidate. This is the part of #184's blast radius
    the first draft of the spec missed, and it is why this change is gated
    rather than treated as cosmetic.
    """
    df = kestrel_stage_frame("final", rows=2)
    # Row 0: exactly at the boundary -- High_Precision* before, Low_Precision after.
    df.loc[0, ["Estimated_Depth_AlternateVariant", "Estimated_Depth_Variant_ActiveRegion"]] = (103, 20000)
    # Row 1: comfortably above it, unaffected by the change.
    df.loc[1, ["Estimated_Depth_AlternateVariant", "Estimated_Depth_Variant_ActiveRegion"]] = (200, 20000)

    scored = calculate_depth_score_and_assign_confidence(df, _shipped_config())
    best = select_single_best_variant(scored)

    assert len(best) == 1
    assert best["Confidence"].iloc[0] == "High_Precision*"
    assert best.index[0] == 1, "the boundary row must no longer win selection"


def test_cond5_is_unreachable_at_the_boundary():
    """Records why only High_Precision* can move, and proves it arithmetically.

    0.00515 == 103/20000 in lowest terms, so any integer pair reaching it has
    alt as a multiple of 103. The smallest is 103, which is already above
    alt_mid_high (100), so cond5's [21, 100) window can never contain a
    boundary row. If this ever fails, the config thresholds changed and the
    #184 blast-radius analysis must be redone.
    """
    from fractions import Fraction

    ratio = Fraction(0.00515).limit_denominator(10**6)
    assert ratio == Fraction(103, 20000)
    assert ratio.numerator > _shipped_config()["confidence_assignment"]["alt_depth_thresholds"]["mid_high"]
```

- [ ] **Step 2: Run them and confirm they fail**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/test_confidence_boundaries.py -v
```

Expected: FAIL, `assert 'High_Precision' == 'Low_Precision'`.

- [ ] **Step 3: Reorder the assignments**

Replace lines 128–134 and 147–154 with:

```python
    # Condition 3: the mid-band demotion. Retained as the expression that names
    # the intent -- @hassansaei on #184: "do not remove this intent".
    cond3 = df["Estimated_Depth_AlternateVariant"].between(
        alt_mid_low,
        alt_mid_high,
        inclusive="left",
    ) & df["Depth_Score"].between(low_threshold, high_threshold)

    # ...

    # Apply conditions in order (later conditions can overwrite earlier ones).
    df.loc[cond1 & above_min_threshold, "Confidence"] = low_prec_label
    df.loc[cond2 & above_min_threshold, "Confidence"] = high_prec_star_label
    df.loc[cond3 & above_min_threshold, "Confidence"] = low_prec_label
    df.loc[cond4 & above_min_threshold, "Confidence"] = low_prec_label
    df.loc[cond5 & above_min_threshold, "Confidence"] = high_prec_label

    # #184: the mid-band demotion must survive the High tiers, at every alt
    # depth. @hassansaei: "Any variant with Depth_Score between 0.00469 and
    # 0.00515 (inclusive) must be Low_Precision, even when alternate depth is
    # >= 21 (or higher)." cond6 already did this for the OPEN interval because
    # it is applied last; the only divergence was Depth_Score == high_threshold
    # exactly, which cond2/cond5 promoted. Making the band inclusive at the top
    # closes that and changes nothing else -- in particular it does not leave
    # any row unassigned, which the alternative (cond2/cond5 -> `> high`) would.
    cond_midband = df["Depth_Score"].between(low_threshold, high_threshold, inclusive="both")
    df.loc[cond_midband & above_min_threshold, "Confidence"] = low_prec_label
```

Delete the old `cond6` definition and its assignment; `cond_midband` subsumes it.

- [ ] **Step 4: Run them and confirm they pass**

Expected: all pass, and every pre-existing test in `test_confidence_boundaries.py` and `test_confidence_assignment.py` still passes. **If any pre-existing boundary test fails, stop and report** — the blast radius is larger than one point and this spec is wrong.

- [ ] **Step 5: Run the full unfiltered tier**

Expected: all 1744+ pass.

- [ ] **Step 6: Commit alone**

```bash
git add vntyper/scripts/confidence_assignment.py tests/unit/test_confidence_boundaries.py
git commit -m "fix(confidence): keep mid-band Depth_Score at Low_Precision, per #184

@hassansaei on #184: \"Any variant with Depth_Score between 0.00469 and 0.00515
(inclusive) must be Low_Precision, even when alternate depth is >= 21.\"

Measured blast radius: Depth_Score == 0.00515 exactly. cond6 was already the
last assignment and demoted the open interval at every alt depth, so only the
top boundary diverged from the decision. Gated by the golden-cohort run at the
tip of this phase, which is expected to show no delta -- exact equality with
0.00515 is essentially unreachable on real data, so the boundary table test is
the load-bearing evidence, not the gate.

Closes #184"
```

- [ ] **Step 7: Post the finding to #184**

```bash
gh issue comment 184 --repo hassansaei/VNtyper --body "..."
```

State that the band interior was already correct because `cond6` is applied last, that the divergence was the `0.00515` boundary only, and that the `> high` variant he suggested would have dropped `alt >= 100, DS == 0.00515` to `Negative`. His answer goes on the record beside what the code did.

---

### Task 18: #185 — a missing gate column raises

**Files:**
- Modify: `vntyper/scripts/kestrel_genotyping.py:806-830`
- Modify: `tests/unit/test_kestrel_filtering.py`
- Test: `tests/unit/test_kestrel_filtering.py`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing importable.

**Decision**, quoted from @hassansaei on #185:

> "Prefer fail loud: a missing required gate column should raise (abort the run), not be skipped. […] 1.3 had no boolean gates, but its final filters accessed columns directly and crashed with KeyError if anything was missing — it never failed open. The current 2.x behaviour […] can silently permit variants that should have been filtered. That is not acceptable for this pipeline. So: do not keep fail-open. […] If a softer path is ever needed for empty early-return frames, that should be an explicit, documented empty-result path — not silent omission of safety gates on a non-empty candidate table."

**Blast radius, established by reading the call graph.** `filter_final_dataframe` has exactly one caller, `kestrel_genotyping.py:594`, and the six stages before it each end in `if df.empty: return df` (lines 550, 555, 560, 565, 573, 578). A frame reaching line 594 is therefore non-empty and has traversed every stage. His empty-frame carve-out is already structurally satisfied.

- [ ] **Step 1: Write the failing tests**

```python
def test_a_missing_gate_column_raises_on_a_non_empty_frame():
    """#185: absence of a gate is an error, not a permit.

    Specification. @hassansaei: "a missing required gate column should raise
    (abort the run), not be skipped [...] That is not acceptable for this
    pipeline." Previously the column was silently skipped, so an upstream
    stage that stopped emitting its column turned a check into a permit.
    """
    df = kestrel_stage_frame("final").drop(columns=["motif_filter_pass"])

    with pytest.raises(ValueError, match="motif_filter_pass"):
        filter_final_dataframe(df, str(tmp_path))


@pytest.mark.parametrize(
    "column",
    ["is_frameshift", "is_valid_frameshift", "depth_confidence_pass", "alt_filter_pass", "motif_filter_pass"],
)
def test_every_gate_column_is_required(column, tmp_path):
    df = kestrel_stage_frame("final").drop(columns=[column])
    with pytest.raises(ValueError, match=column):
        filter_final_dataframe(df, str(tmp_path))


def test_an_empty_frame_is_the_documented_empty_result_path(tmp_path):
    """His carve-out: empty early-return frames are not an error.

    Structurally this cannot be reached in production -- the six stages before
    the only caller (kestrel_genotyping.py:594) each return early on an empty
    frame -- but the contract is stated here rather than left implicit.
    """
    out = filter_final_dataframe(pd.DataFrame(), str(tmp_path))
    assert out.empty


def test_the_pre_result_tsv_is_still_written_before_the_raise(tmp_path):
    """kestrel_pre_result.tsv is the debuggability artefact (AGENTS.md trap 4).

    It must be on disk before the raise, or a failing run loses the evidence
    of what was in the frame.
    """
    df = kestrel_stage_frame("final").drop(columns=["alt_filter_pass"])
    with pytest.raises(ValueError):
        filter_final_dataframe(df, str(tmp_path))
    assert (tmp_path / "kestrel_pre_result.tsv").exists()
```

Use `tests/builders.py`'s `kestrel_stage_frame("final")` — it returns exactly the columns present at that point.

- [ ] **Step 2: Run them and confirm they fail**

Expected: `DID NOT RAISE <class 'ValueError'>` for the first two, PASS for the third and fourth.

- [ ] **Step 3: Fail closed**

Replace the `else` at `kestrel_genotyping.py:828-829`:

```python
        else:
            # #185: a missing gate is an error, not a permit. @hassansaei:
            # "a missing required gate column should raise (abort the run), not
            # be skipped [...] That is not acceptable for this pipeline."
            # Reachability: this function's only caller is at line 594, behind
            # six `if df.empty: return df` guards, so any frame arriving here is
            # non-empty and has traversed every stage that adds a gate column.
            # An empty frame short-circuits before this point -- that is the
            # explicit empty-result path his decision carves out.
            msg = (
                f"Required filter column '{col}' is missing from a non-empty Kestrel result frame. "
                "An upstream stage stopped emitting it, so its safety gate would silently become a "
                "permit. Aborting rather than reporting unfiltered variants. See issue #185."
            )
            logger.error(msg)
            raise ValueError(msg)
```

Guard the whole loop with `if df.empty: return df` at the top of the function, **after** the `kestrel_pre_result.tsv` write, so the empty path stays documented and the artefact is still produced.

- [ ] **Step 4: Run them and confirm they pass**

- [ ] **Step 5: Prove the raise reaches exit 1 rather than being swallowed**

This is the named risk. `run_pipeline` uses `except Exception` at stage boundaries by convention; a raise swallowed there and converted into a silent no-result is **worse** than the fail-open it replaces.

Write an end-to-end test that drives the pipeline entry point with a mocked Kestrel stage returning a frame missing `motif_filter_pass`, and assert the process exits 1 and the message reaches the log at ERROR:

```python
def test_a_missing_gate_column_aborts_the_run_rather_than_reporting_nothing(caplog, tmp_path):
    caplog.set_level(logging.ERROR, logger="vntyper.scripts.kestrel_genotyping")
    with pytest.raises(SystemExit) as exc:
        ...
    assert exc.value.code == 1
    assert any("issue #185" in r.message for r in caplog.records)
```

**If the exception is swallowed, stop and report before committing.** Do not weaken the raise to work around it.

- [ ] **Step 6: Run the full unfiltered tier**

- [ ] **Step 7: Commit alone**

```bash
git add vntyper/scripts/kestrel_genotyping.py tests/unit/test_kestrel_filtering.py
git commit -m "fix(kestrel): raise when a required gate column is missing, per #185

@hassansaei on #185: \"a missing required gate column should raise (abort the
run), not be skipped [...] That is not acceptable for this pipeline.\"

Reachability: the only caller is kestrel_genotyping.py:594, behind six empty
guards, so every frame arriving here is non-empty and has traversed all six
stages. The empty-result path his answer carves out short-circuits earlier and
is now stated explicitly. kestrel_pre_result.tsv is still written before the
raise so a failing run keeps its evidence.

Closes #185"
```

---

### Task 19: #192 — sum every LEN token in a compound adVNTR state

**Files:**
- Modify: `vntyper/modules/advntr/advntr_genotyping.py:14-42` (the pattern constants), `:163-169`, `:202-208`
- Modify: `tests/unit/test_advntr_output_parsing.py`
- Create: a differential-sweep script under `scripts/` (not shipped in the package)
- Test: `tests/unit/test_advntr_output_parsing.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `sum_insertion_lengths(variant: str) -> int` in `advntr_genotyping.py`.

**Decision**, quoted from @hassansaei on #192:

> "For compound adVNTR states […] use the **sum** of inserted lengths and the **sum** of deleted lengths when computing the net length that feeds the pathogenic frameshift filter (`frame = abs(Insertion_len - Deletion_length)`, with the same 3n+1 / 3n+2 rule as #182). Example: `I9_2_A_LEN9&I50_2_A_LEN3` to Insertion_len = 9 + 3 = 12 (not first-LEN-only). […] If a full multi-event parse is too brittle for some malformed strings, a acceptable fallback is: still require that the **summed** insertion/deletion lengths satisfy the pathogenic frameshift rule, and if they do, report the variant/state as adVNTR emitted it (do not invent a simplified single event). Do not keep first-LEN-wins as the defined semantics."

**This changes reported output**, because `Insertion_len` feeds `frame = abs(Insertion_len - Deletion_length)`, which is the value the #182 filter tests.

**The golden cohort cannot exercise it.** The only compound state in the cohort is `example_dfc3`'s `D17_2&D18_2&D19_2&D20_2&D21_2`, which contains no `LEN` token; `Insertion_len` is `0` under both semantics. A PASS in Task 21 is **not** evidence for this change. Step 6 is.

**Open sub-question to record, not to resolve.** `Deletion_length` is `Variant.str.count("D")` — a count of deletion *events* used arithmetically as a length in bp. That equals "the sum of deleted lengths" only if every adVNTR `D` event is a single base. That reading is consistent with the State grammar and with his example, but it is an inference. Put it in the commit message and post it to #192.

- [ ] **Step 1: Write the failing tests**

```python
**Before writing these, know what the current behaviour actually is.** #192's text, the #198 PR body and this plan's first draft all say "first-LEN-wins". **That is wrong.** Measured by running `advntr_genotyping.py:163-169` directly:

| State | `Insertion_len` today |
| --- | ---: |
| `I22_2_G_LEN1` | 1 |
| `I9_2_A_LEN9&I50_2_A_LEN3` | **0** — not 9 |
| `I9_2_A_LEN2&D50_2` | **0** — not 2 |
| `I9_2_A_LEN9&` | **0** — not 9 |
| `I50_2&I9_2_A_LEN3` | **3** — the `LEN` is terminal, so it parses |
| `I50_2&D9_2&I80_2_A_LEN7` | **7** — likewise |

**The precise rule: the value collapses to zero exactly when text follows the first `LEN` token.** Not "every compound state" — a compound whose only `LEN` is terminal parses correctly. Where it does collapse, `frame == 0` for a pure-insertion compound, and 0 is in neither `ins_frame` (3n+1) nor `del_frame` (3n+2), so those states are dropped in silence.

**"first-LEN-wins" describes no input at all.** State that in the commit message and in the #192 close comment; it is a worse defect than the issue describes, and leaving the wrong diagnosis in place is how the next reader inherits it.

**Expected-difference oracle for Step 6, in these terms:** a state differs before/after **iff material follows its first `LEN` token** — a second `LEN`, a further `&` part, or a trailing `&`. States with no `LEN`, a single terminal `LEN`, or a `LEN` only in the last part must be identical. Anything else differing is a regression; stop and report.

```python
@pytest.mark.parametrize(
    ("state", "expected"),
    [
        ("I22_2_G_LEN1", 1),
        ("I9_2_A_LEN9&I50_2_A_LEN3", 12),          # his worked example: 9 + 3
        ("I9_2_A_LEN2&D50_2", 2),
        ("D17_2&D18_2&D19_2&D20_2&D21_2", 0),      # the cohort's compound call, no LEN
        ("I9_2_A_LEN9&I50_2_A_LEN3&I80_2_A_LEN1", 13),
        ("I9_2_A_LEN9&", 9),                        # trailing & no longer collapses to 0
        ("I50_2&I9_2_A_LEN3", 3),                   # terminal LEN: already 3 today, unchanged
        ("I50_2&D9_2&I80_2_A_LEN7", 7),             # terminal LEN in the last part, unchanged
        ("", 0),
        ("D50_2", 0),
    ],
)
def test_every_len_token_is_summed(state, expected):
    """#192: sum, not first-LEN-wins -- and not the zero it actually produces.

    Specification. @hassansaei: "use the sum of inserted lengths [...]
    Example: I9_2_A_LEN9&I50_2_A_LEN3 to Insertion_len = 9 + 3 = 12 (not
    first-LEN-only). [...] Do not keep first-LEN-wins as the defined
    semantics."
    """
    assert sum_insertion_lengths(state) == expected


def test_the_maintainers_worked_example_does_not_by_itself_change_survival():
    """Recorded so the example is not mistaken for the regression test.

    9 + 3 = 12 and 12 % 3 == 0, so this state is filtered out under BOTH
    semantics. It pins the arithmetic he specified; it does not demonstrate a
    changed outcome. The two tests below do that.
    """
    assert sum_insertion_lengths("I9_2_A_LEN9&I50_2_A_LEN3") == 12
    assert 12 % 3 == 0


def test_a_compound_insertion_that_summation_admits_is_dropped_today():
    """The insertion path: summation can only ADD rows, never remove them.

    frame = abs(Insertion_len - Deletion_length); Deletion_length is
    count("D") = 0 here; the #182 insertion filter keeps 3n+1.

      I9_2_A_LEN9&I50_2_A_LEN1
        today:     Insertion_len = 0,  frame = 0,  0 in neither set -> DROPPED
        summation: Insertion_len = 10, frame = 10, 10 % 3 == 1      -> KEPT

    Because every compound state is 0 today and 0 is in neither frame set, no
    compound insertion currently survives at all. Summation is therefore
    monotone-additive on this path: it makes real pathogenic compound calls
    visible and cannot lose one.
    """
    kept = advntr_processing_ins(_frame_with_state("I9_2_A_LEN9&I50_2_A_LEN1"), _advntr_config())
    assert len(kept) == 1
    assert kept["Insertion_len"].iloc[0] == 10


def test_a_compound_deletion_state_can_LOSE_a_row_under_summation():
    """The deletion path is NOT monotone. Asserted so it is not discovered later.

      I9_2_A_LEN3&D50_2&D51_2   (Deletion_length = count("D") = 2)
        today:     Insertion_len = 0, frame = |0-2| = 2, 2 in del_frame -> KEPT
        summation: Insertion_len = 3, frame = |3-2| = 1, 1 not in del   -> DROPPED

    @hassansaei's decision on #192 authorises this -- "the biologically
    meaningful quantity for keep/ignore is the net indel length on that state"
    -- and under #181/#182 a net Delta of +1 in the deletion direction is the
    non-pathogenic frame. But it is a reported call disappearing, so it is
    pinned here and called out in the commit message and the PR body.
    """
    kept = advntr_processing_del(_frame_with_state("I9_2_A_LEN3&D50_2&D51_2"), _advntr_config())
    assert len(kept) == 0
```

Both directions must be asserted, and they are **different code paths**: the insertion path can only add rows, the deletion path can add or remove. A suite that only exercised `advntr_processing_ins` would leave the removing direction unmeasured.

Verify each expected value by running the current parser before asserting it — the "first-LEN-wins" framing in #192 and in the #198 PR body is wrong, so any expectation derived from it is wrong too.

- [ ] **Step 2: Run them and confirm they fail**

Expected: `NameError: name 'sum_insertion_lengths' is not defined`, then value mismatches (`12 != 9`).

- [ ] **Step 3: Implement the summation**

```python
#: Every ``LEN<n>`` token in an adVNTR ``State`` string. #192 defines the
#: inserted length of a compound state as the SUM of these, not the first.
#: @hassansaei: "I9_2_A_LEN9&I50_2_A_LEN3 to Insertion_len = 9 + 3 = 12".
LEN_TOKEN_PATTERN = re.compile(r"LEN(\d+)")


def sum_insertion_lengths(variant: str) -> int:
    """Total inserted bases named by every LEN token in an adVNTR State string.

    Args:
        variant (str): An adVNTR ``State``/``Variant`` string, possibly compound
            (parts joined by ``&``).

    Returns:
        int: The sum of every ``LEN<n>``; 0 when there is none.
    """
    if not isinstance(variant, str):
        return 0
    return sum(int(match) for match in LEN_TOKEN_PATTERN.findall(variant))
```

Replace lines 163–169 and 202–208 with `df1["Insertion_len"] = df1["Variant"].map(sum_insertion_lengths).astype(int)`. Delete `GREEDY_INSERTION_LEN_PATTERN` and `INSERTION_LEN_SPLIT_LIMIT` and the `I` scratch column they produced, along with their now-obsolete comments.

- [ ] **Step 4: Run them and confirm they pass**

- [ ] **Step 5: Run the full unfiltered tier**

Expected: all pass. `test_advntr_output_parsing.py` and `test_advntr_frameshift_filter.py` in particular.

- [ ] **Step 6: Run the differential sweep — this is the evidence, not the gate**

Write `scripts/advntr_len_differential.py` that generates State strings across the grammar (single insertion, single deletion, 2–5-part compounds, mixed I/D, trailing `&`, malformed `LEN`, no `LEN`), runs both the old and new parser over each, and reports every input whose `Insertion_len`, `frame`, or keep/drop outcome differs.

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  python scripts/advntr_len_differential.py --out /tmp/advntr-len-diff.json
```

Expected, per the oracle above: differences **iff material follows the first `LEN` token**. **Any difference on a no-`LEN` state, a single terminal-`LEN` state, or a compound whose only `LEN` is in its last part is a regression — stop and report.** The generator must include all three of those unchanged classes, or the sweep proves nothing about them. Record the counts in the commit message, in the form the `d144505` sweep used (N probes, M previously-parsing inputs identical).

- [ ] **Step 7: Commit alone**

```bash
git add vntyper/modules/advntr/advntr_genotyping.py tests/unit/test_advntr_output_parsing.py scripts/advntr_len_differential.py
git commit -m "fix(advntr): sum every LEN token in a compound state, per #192

@hassansaei on #192: \"use the sum of inserted lengths [...] I9_2_A_LEN9&
I50_2_A_LEN3 to Insertion_len = 9 + 3 = 12 (not first-LEN-only). [...] Do not
keep first-LEN-wins as the defined semantics.\"

Insertion_len feeds frame = abs(Insertion_len - Deletion_length), so this
changes which compound rows survive the #182 filter -- reported output.

Evidence is a differential sweep, NOT the golden cohort: the cohort's only
compound state (dfc3's D17_2&...&D21_2) carries no LEN token, so Insertion_len
is 0 under both semantics and the gate cannot reach this code. <N> probes,
<M>/<M> single-LEN and no-LEN inputs byte-identical, differences confined to
multi-LEN states and the trailing-& case.

Recorded inference: Deletion_length remains Variant.str.count(\"D\"), i.e. a
count of deletion events used as a length in bp. That equals \"the sum of
deleted lengths\" only if every adVNTR D event is a single base. Consistent with
the State grammar and with his example, but posted to #192 for confirmation.

Closes #192"
```

- [ ] **Step 8: Post the `Deletion_length` inference to #192**

---

### Task 20: **CUT — file it as an issue instead**

**Do not implement this task.** The first draft removed the divergent `.get()` calibration defaults from `confidence_assignment.py` and justified it as an extension of #185's "prefer fail loud". The adversarial review rejected that, correctly: **#185 decided gate columns, not calibration keys**, and this plan's own constraint forbids changing confidence code without an explicit decision. Stretching one decision to cover a different site is exactly the overreach this work exists to catch in others.

**What to do instead:** file it as a question in Task 24, with the table below as the body, and implement it only once answered. The six mutants stay hand-classified in `scripts/mutation_test.py` in the meantime. Do not touch `confidence_assignment.py:79-97`.

The case, for the issue body — today a missing calibration key silently substitutes a value up to **78× wrong**:

| Key | Shipped | `.get()` default | Ratio |
| --- | ---: | ---: | ---: |
| `depth_score_thresholds.low` | 0.00469 | 0.2 | 43× |
| `depth_score_thresholds.high` | 0.00515 | 0.4 | 78× |
| `alt_depth_thresholds.low` | 20 | 5 | — |
| `alt_depth_thresholds.mid_low` | 21 | 10 | — |
| `alt_depth_thresholds.mid_high` | 100 | 20 | — |
| `var_active_region_threshold` | 200 | 0 | — |

Behaviour is identical under the shipped config, which supplies all six keys; it differs only for a partial config, which AGENTS.md trap 2 already documents as unsupported input that fails with `KeyError` elsewhere. It would **delete** six mutants the artefact currently hand-classifies as equivalent, rather than explaining them away.

The steps below are retained **only** as the implementation sketch to attach to the issue, so that whoever answers it can see exactly what is proposed. **None of them is executed in this plan.**

  * ~~**Step 1: Write the failing test**~~

```python
def test_a_config_missing_a_calibration_constant_raises_rather_than_guessing():
    """A partial config must not silently substitute a different calibration.

    Extension of the principle @hassansaei stated on #185 ("prefer fail loud"),
    applied to the calibrated constants rather than the gate columns. Not a
    decision he gave. AGENTS.md trap 2 records that --config-path replaces the
    whole config rather than merging it and that a partial config already fails
    with KeyError elsewhere, so a config missing these keys is not supported
    input -- it should fail here too, and loudly.
    """
    config = _shipped_config()
    del config["confidence_assignment"]["depth_score_thresholds"]["low"]

    with pytest.raises(KeyError, match="low"):
        calculate_depth_score_and_assign_confidence(_one_row(), config)
```

  * ~~Step 2: Run it and confirm it fails~~  *(not executed)*

Expected: `DID NOT RAISE` — today it silently uses `0.2`.

  * ~~Step 3: Require the keys~~  *(not executed)*

Replace lines 79–97's `.get(key, DEFAULT)` calls with direct subscripts, leaving the `confidence_levels` label fallbacks alone (those are cosmetic strings, not calibration):

```python
    conf_assign = kestrel_config["confidence_assignment"]
    thresholds = conf_assign["depth_score_thresholds"]
    alt_thresholds = conf_assign["alt_depth_thresholds"]
    var_region_threshold = conf_assign["var_active_region_threshold"]
    confidence_levels = conf_assign.get("confidence_levels", {})

    # Labels are cosmetic and may default; the CALIBRATION may not. A missing
    # calibration key previously substituted a value up to 78x wrong (0.4 for
    # a shipped 0.00515) and silently relabelled every variant in the run.
    low_prec_label = confidence_levels.get("low_precision", "Low_Precision")
    high_prec_label = confidence_levels.get("high_precision", "High_Precision")
    high_prec_star_label = confidence_levels.get("high_precision_star", "High_Precision*")

    low_threshold = thresholds["low"]
    high_threshold = thresholds["high"]
    alt_low = alt_thresholds["low"]
    alt_mid_low = alt_thresholds["mid_low"]
    alt_mid_high = alt_thresholds["mid_high"]
```

  * ~~Step 4: Run it and confirm it passes~~  *(not executed)*

  * ~~Step 5: Run the full unfiltered tier~~  *(not executed)*

Expected: all pass. Any test that relied on a partial config now raises — fix the **test** to supply the shipped config via `tests/builders.py`'s `kestrel_config(**dotted_overrides)`.

  * ~~Step 6: Delete the six equivalent-mutant classifications~~  *(not executed)*

Remove the `confidence_assignment.py` lines 82/91/92/95/96/97 entries from `EQUIVALENT_MUTANTS` in `scripts/mutation_test.py`. They are no longer equivalent — they no longer exist.

  * ~~Step 7: Re-render the mutation page~~  *(not executed)*

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  make mutation-render
```

Seconds, not a re-sweep — classification changes how the measurement is presented, not the measurement.

  * ~~Step 8: Commit alone~~  *(not executed)*

```bash
git add vntyper/scripts/confidence_assignment.py tests/unit/test_calibration_consistency.py \
        scripts/mutation_test.py docs/development/mutation-testing.md
# NOT RUN. Reproduced only so the issue can show what was proposed.
# git commit -m "fix(confidence): require the calibration constants instead of defaulting them

The .get() fallbacks encoded a second, wrong calibration -- 0.4 where the
shipped config says 0.00515, a factor of 78 -- reachable by dropping one key.
Behaviour is unchanged under the shipped config.

This deletes six mutants the artefact had to hand-classify as equivalent rather
than explaining them away.

Extension of the principle @hassansaei stated on #185 (\"prefer fail loud\"),
not a decision he gave on this code."
```

---

### Task 21: Golden-cohort attestation — **one run per gated commit, plus one at the PR tip**

**Files:**
- Modify: `docs/development/golden-cohort-gate.md`

**Interfaces:**
- Consumes: each gated commit in turn (Tasks 9, 10, 17, 18, 19, 6a) and the baseline SHA recorded in Task 0 Step 5.
- Produces: the attestations of record for this branch.

**Read `docs/development/golden-cohort-gate.md` in full before starting.** The harness exists; do not rebuild it.

**This task runs six-to-seven times, not once.** The first draft ran it once over four commits jointly. Two reasons that was wrong, both from the adversarial review and both accepted:

- **A joint run cannot attribute a failure.** "Individually revertible" tells you how to undo a commit, not which one broke the cohort.
- **Two changes can produce offsetting deltas** that a joint run reports as clean.

Each run's steps are identical; only the candidate SHA changes. At ~25 minutes each, six runs is about 2.5 hours — that is the price of attribution, and it is affordable.

**Also fix the existing page while you are in it.** `golden-cohort-gate.md:121`, `:140` and `:219` state that runs 2 and 3 compared `Insertion_len` across the adVNTR rows — `:140` goes further and says the row is identical "in every column, `Insertion_len` included". All three are unsupported. That column is **not** in the adVNTR output schema — `final_columns` at `advntr_genotyping.py:365` is `VID, Variant, NumberOfSupportingReads, MeanCoverage, Pvalue, RU, POS, REF, ALT, Flag`, and `base_columns` at `:432` is shorter still. The row-set comparison remains valid; the `Insertion_len` claim specifically is unsupported. Strike it rather than repeat it — an attestation that over-claims one column is the reason the whole page has to be read sceptically.

- [ ] **Step 1: Confirm the tree is clean and name the candidate**

```bash
git status --porcelain
git rev-parse HEAD
```

The candidate is this SHA, and **only** this SHA. **A verdict attests one candidate commit and nothing after it.** Record it in the run's own section before doing anything else — a run whose candidate is written down afterwards is a run whose candidate can drift.

**Nothing may be committed to this branch between Step 1 and Step 7.** If something is, the run is void and restarts.

- [ ] **Step 2: Build the baseline worktree at the Task 0 SHA**

```bash
git worktree add /tmp/gate4/before "$(cat /tmp/golden-baseline.sha)"
```

Symlink only the untracked data and reference directories in; the baseline must carry its **own** tracked `vntyper/config.json` and report/adVNTR configuration.

- [ ] **Step 3: Prove which code each side runs, before running anything**

This is the step that decides whether the exercise is worth anything. The `vntyper` console script resolves through setuptools' editable finder, which is *appended* to `sys.meta_path` and points at whichever worktree the editable install was made from — irrespective of CWD. Demonstrate the failure mode, then defeat it:

```bash
cd /tmp/gate4/before && python -c "import vntyper, importlib.util; print(vntyper.__file__, importlib.util.find_spec('vntyper.scripts.motif_decisions') is not None)"
```

Expected: it resolves to the **candidate** worktree — that is the trap. Then pin `PYTHONPATH` per side and confirm each resolves to its own tree, using `vntyper.scripts.motif_decisions` (created in Task 10) as the marker module: absent on the before side, present on the after side.

Every run must launch through a wrapper that prints its resolved `vntyper.__file__` and marker state as its first line and `sys.exit`s before dispatch unless both agree with its side.

- [ ] **Step 4: Run the 58-case matrix on both sides**

Reproduce the matrix from the existing gate runs verbatim: 7 multi-reference samples × 6 assemblies, their hg19 subsets, `example_40cf` at hg38, 5 non-fast-mode repeats, 3 `--extra-modules advntr` cases at `--advntr-max-coverage 300`, plus the 3 probes. 122 runs total. ~25 minutes.

- [ ] **Step 5: Compare**

Per case: the full `kestrel_result.tsv` row set keyed on `Motifs`/`POS`/`REF`/`ALT`/`Variant`, `kestrel_pre_result.tsv`, `output_adVNTR_result.tsv`, `coverage_summary.tsv`, the screening-summary sentence and its computed emphasis, the recorded pipeline steps, the executed shell command strings, and the exit code.

Expected: **0/58 deltas** on the Kestrel variant set, `Confidence`, `Flag`, adVNTR genotype fields and `Flag`.

- [ ] **Step 6: Write the run up as run 4**

Add a `## Result — run 4, candidate <sha>` section following the existing structure, and add the row to the runs table at the top. It **must** state, in the "What this gate does not cover" section:

- Task 19 (#192) is **not exercised**: no cohort state carries a `LEN` token, so `Insertion_len` is 0 on both sides. Its evidence is the differential sweep, cited by name.
- Task 17 (#184) is **almost certainly not exercised**: `Depth_Score == 0.00515` exactly is essentially unreachable. Report the observed `Depth_Score` range across the cohort so the claim is measured rather than asserted.
- Task 18 (#185) is exercised only in the negative — no cohort case is missing a gate column, so the run shows the raise does not fire on healthy input, which is what it is for.

**Do not write "no regression" without naming what the run could not see.**

- [ ] **Step 7: Remove the baseline worktree and commit**

```bash
git worktree remove /tmp/gate4/before
git add docs/development/golden-cohort-gate.md
git commit -m "docs(gate): attest the golden cohort at <sha>, and name what it could not reach"
```

---

# Phase 4 — Split the oversized files

**Executes BEFORE phase 5, not after.** The first draft put the splits after the gated phase, which would have invalidated its attestation — the gate attests one commit and nothing after it. Only Task 22 runs; Task 23 stays parked on the base-image trigger, so this phase is one lane, not two.

### Task 22: Split `cohort_summary.py` (911 LOC)

**Lane:** 22

**Files:**
- Modify: `vntyper/scripts/cohort_summary.py`
- Create: focused modules for the pure parts it extracts
- Test: new test files per extracted module

**Interfaces:**
- Consumes: `vntyper/scripts/summary_steps.py`'s step-name constants (trap 5 — never a bare literal).
- Produces: whatever it extracts. Name each in the commit message.

346 missing units at 37.8% branch-inclusive; 911 LOC against a ~650 limit. Pull the **pure** logic out first — aggregation, formatting, the per-column escaping exemptions — and leave the filesystem and Jinja2 rendering behind. `scoring.py`, `confidence_assignment.py` and `region_utils.py` are the shape to copy.

**Do not attempt a whole-file rewrite.** Split out one region at a time, test it, move on.

- [ ] **Step 1: Pin the current output before touching anything**

Build a two-sample cohort under `tmp_path` from `tests/builders.py`'s `kestrel_output_dir`, render the cohort report, and hash the result. That hash is the oracle for every step below. If it has to change, the refactor was not behaviour-preserving.

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  pytest tests/unit/test_cohort_summary_escaping.py -v
```

Expected: green, unmodified, at every step from here on.

- [ ] **Step 2: List the pure regions and pick the first**

```bash
grep -n "^def \|^class " vntyper/scripts/cohort_summary.py
```

Classify each as pure (aggregation, formatting, column selection, the escaping exemptions) or impure (filesystem walk, Jinja2 render, argument handling). Extract pure ones only. Take the one with the most uncovered units first — get that list from `--cov-report=term-missing`.

- [ ] **Step 3: Write the failing tests for the extracted module**

Against the new module path, so they fail with `ModuleNotFoundError` before the move.

- [ ] **Step 4: Run them and confirm they fail**

Expected: `ModuleNotFoundError`.

- [ ] **Step 5: Move the code verbatim; do not rewrite it**

Cut, do not retype. Parameterise only what it read from an enclosing scope. Fully annotate the new module.

- [ ] **Step 6: Run the new tests, the Step 1 oracle, and the escaping suite together**

All three must pass, and **the oracle must pass unmodified**. If it needs editing, stop and report.

- [ ] **Step 7: Repeat Steps 2–6 per region until `cohort_summary.py` is under 650 LOC**

Commit after each region, so a bad extraction is revertible on its own.

- [ ] **Step 8: Re-measure, raise the floor to the printed figure, run `make ci-local`, hand off**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" make test-unit-cov
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" make ci-local
wc -l vntyper/scripts/cohort_summary.py
```

Expected: `ci-local` exit 0, `cohort_summary.py` under 650 lines, floor raised to the printed figure.

---

### Task 23: Split `install_references.py` (901 LOC) — **deferred, base-image gated**

**Lane:** 23

**`vntyper/scripts/install_references.py` is a base-image content-hash input.** Editing it — even a comment — turns `Build base image` and `Docker Success` from skipped to failed until a new base publishes. That is why `9104f64` reverted the aligner-index quoting.

**This task does not run until a base-image rebuild is scheduled.** When it is, it lands in the same push as:

- **#193** — declare `bcrypt` directly in `docker/requirements-web.txt` and drop `passlib[bcrypt]` if nothing else needs it. `docker/app/utils.py` imports `bcrypt` directly; the dependency is currently undeclared and satisfied transitively.
- **The aligner-index quoting** — `execute_aligner_index` interpolates `ref_path`, `index_dir` and `index_base` into a template run with `shell=True`. Apply `command_builders.quote_path` per path. `tests/unit/test_shell_quoting.py`'s three `execute_aligner_index` cases characterise the unquoted behaviour and **will fail loudly** once the quoting lands — that is by design; update them to specification then.

- [ ] **Step 1: Confirm a base rebuild is scheduled and merge to `main` (or run `docker-base.yml` via `workflow_dispatch`) so the base publishes first**
- [ ] **Step 2: Locally, `make docker-build-base && make docker-build DOCKER_BASE_IMAGE=vntyper-base:local`**
- [ ] **Step 3: Land #193, the quoting, and the split in one push; close #193**

---

# Phase 6 — Close the loop

### Task 24: Refresh AGENTS.md, file what was found, close what is closed

**Files:**
- Modify: `AGENTS.md`

- [ ] **Step 1: Re-measure every row of the LOC/coverage table**

```bash
wc -l vntyper/scripts/*.py docker/app/*.py | sort -rn | head -20
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" \
  make test-unit-cov
```

The table is stale and its central claim — "every module over 650 lines is under 30% covered" — no longer holds: `generate_report.py` is 574 LOC (was 861) at 64.8%, `cohort_summary.py` is 911 not 856, `kestrel_genotyping.py` 843 at 41.9%. A table that is wrong stops being trusted, and this one is the argument for rule 2.

- [ ] **Step 2: Fix the four other stale statements in AGENTS.md**

| Line | Says | Reality |
| --- | --- | --- |
| Code style | "Code must run on Python 3.9 (CI matrix: 3.9–3.12)"; `target-version = "py39"` | `pyproject.toml:10` is `>=3.10` and `:157` is `py310`; #198's CI ran **3.10, 3.11, 3.12, 3.13** |
| `:151` | equivalent mutants "have not been hand-classified" | `docs/development/mutation-testing.md:21` hand-classifies ten, each with its reason |
| Trap 11 | `vntyper report` is broken | fixed on #198 — remove it and renumber, or mark resolved |
| Traps 14, 14 | two traps share the number 14 | renumber |

The Python floor matters most: it is the constraint every contributor writes code against, and it has been wrong since the floor moved.

- [ ] **Step 3: File the line-332 defect**

`motif_processing.py:332` — `working_df["Motifs"].str.count("-").max() != 1` is a column-wide gate: **one malformed motif ID sets `motif_filter_pass = False` for every row of the sample.** This is W8 from #179 §3.1 and was never filed. The `or` → `and` mutant at that line survives.

```bash
gh issue create --repo hassansaei/VNtyper \
  --label needs-domain-decision \
  --title "One malformed motif ID nulls the whole sample: the Motifs dash check is column-wide" \
  --body '`vntyper/scripts/motif_processing.py:332`:

```python
if "Motifs" not in working_df.columns or working_df["Motifs"].str.count("-").max() != 1:
    logger.error("Cannot split '"'"'Motifs'"'"' into left-right. Old code returns empty.")
    combined_df = pd.DataFrame(columns=working_df.columns)
```

`.max()` aggregates over the **whole column**. A single row whose `Motifs`
value does not contain exactly one dash therefore makes the condition true for
every row, so `motif_filter_pass` becomes `False` for the entire sample and the
run reports `Negative` — not for the malformed row, but for all of them.

Blast radius: one bad motif ID anywhere in a sample suppresses that sample'"'"'s
whole call. The pipeline exits 0.

This is W1/W8 from the #179 audit (§3.1) and was never filed. Mutation testing
also reaches it: the `or` → `and` mutant at line 332 survives the entire unit
tier.

Characterised by a test added alongside the motif decision-layer extraction, and
**deliberately not fixed**: making the check per-row changes which rows survive
motif filtering, which is a reported-genotype change and needs a domain decision.

What to decide: should a malformed `Motifs` value fail only its own row, fail the
sample loudly (raise), or continue to null the sample silently?'
```

- [ ] **Step 3b: File the six other defects this work characterised but did not fix**

Every one changes what is reported, and none has a decision. All are pinned by a test named `..._today` with a docstring saying so. File each with the source quoted, the blast radius, the pinning test by name, and label `needs-domain-decision`.

| Site | Defect |
| --- | --- |
| `cross_match.py:138-140` | The `eval` of `match_logic` is wrapped in `except Exception`, so a rule naming a missing column reports **"no match"** instead of failing. Trap 3, live. Precedents: `Poylmorhic_Call`, `RU == 7`. |
| `file_processing.py:28-30` | A row is kept only if **exactly one** of REF/ALT has length 1, so `REF="AC", ALT="ACGGG"` — a real 3-base insertion — is discarded before Kestrel post-processing sees it. |
| `file_processing.py:60-63` | The insertion test is `len(ref) == 1 and len(alt) > 1` with an unconditional `else`, so a multi-base-REF insertion is routed to the **deletion** file, where 3n+2 is applied instead of 3n+1. Under #181 those are not interchangeable. Currently masked by the row above — file them as one issue. |
| `extract_unmapped_from_offset.py:38-53` | `f.read(4)` at EOF returns `b""` and `int.from_bytes` gives `0`, so a truncated BAI parses as "no more references" and the unmapped scan starts at offset 0 rather than raising. |
| `motif_processing.py:332` | Step 3 above. |
| `confidence_assignment.py:79-97` | The `.get()` calibration defaults encode a second calibration up to **78× wrong** (`0.4` where the shipped config says `0.00515`), reachable by dropping one key. **Ask, do not assume:** #185 decided *gate columns*, not *calibration keys*, and this plan does not have authority to extend that. Attach Task 20's implementation sketch. |

- [ ] **Step 4: Close every issue this branch resolves, each with its evidence**

`#181 #182 #183 #184 #185 #186 #187 #192 #194 #195 #196 #197`, plus `#188` **only if a CRAM fixture substantiated it** — otherwise #188 joins the parked list below. Each close comment names the commit, the test that pins it, and for the gated ones the golden-cohort run and its candidate SHA.

**Correct #188's and #192's issue text when closing them.** #188 diagnoses a three-production-file change; `--cram` already existed. #192 describes the current behaviour as "first-LEN-wins"; it actually coerces every compound state to zero, which is a worse defect. Leaving a wrong diagnosis in a closed issue is how the next reader inherits it.

- [ ] **Step 5: Comment on the three that stay open with their trigger**

- **#189** — needs a decision on whether per-job access requires a credential distinct from the job id; issuing a per-job token breaks the `vntyper online` CLI contract.
- **#191** — needs a CRAM sample in `tests/data/` and in the golden-cohort matrix before the process-substitution change can be verified.
- **#193** — unparks at the next base-image rebuild, alongside Task 23.

- [ ] **Step 6: Final gate**

```bash
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" make check-all
conda run -n vntyper env PATH="/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH" make ci-local
```

Both must exit 0, with the output shown. **Never claim tests pass without showing the command output.**

- [ ] **Step 7: Write the PR body**

`specs/pr-body.md` does not exist and nothing else creates it. Write it here, and it must contain:

- **A table of every gated commit**, enumerated by task and SHA — not a hardcoded count. The gated set is Tasks 9, 10, 17, 18, 19 and 6a; #188 drops out if it parked. Each row: task, commit SHA, the golden-cohort run that attests it, and **what that run did not exercise**.
- The final tip-run SHA and its verdict.
- For #192: that the golden cohort could not reach it and the differential sweep is the evidence, with its counts.
- For #188: the CRAM end-to-end evidence, or the fact that it parked.
- The six newly filed issues, each with its number.
- The three issues still open with their triggers (#189, #191, #193).
- The corrected diagnoses posted to #188 and #192.

- [ ] **Step 8: Commit and open the PR**

```bash
git add AGENTS.md specs/pr-body.md
git commit -m "docs(agents): refresh the LOC/coverage table, the Python floor and trap 11 from measurement"
gh pr create --repo hassansaei/VNtyper --base main \
  --title "fix: close the #179 follow-ups and enable branch coverage" \
  --body-file specs/pr-body.md
```
