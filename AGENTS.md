# AGENTS.md

Instructions for coding agents working in this repository. Human contributors should
read [CONTRIBUTING.md](CONTRIBUTING.md); everything below is the agent-facing summary
plus the repo-specific traps that are easy to get wrong.

## Project

VNtyper genotypes the MUC1 coding VNTR to diagnose ADTKD-MUC1 from short-read data
(BAM/CRAM/FASTQ). It runs mapping-free k-mer genotyping (Kestrel, vendored JAR),
optionally validates with adVNTR, and emits confidence-scored calls plus an HTML report.

Python package `vntyper`, one console entry point: `vntyper`. Version lives only in
`vntyper/version.py`. Research use only — clinical-sounding output text is
config-driven, never invent or reword it.

## Setup

The pipeline needs external binaries (bwa, samtools, fastp, bcftools, Java 11), so
`pip install` alone does **not** give you a working pipeline:

```bash
mamba env create -f conda/environment_vntyper.yml   # env name: vntyper
conda activate vntyper
make install-dev                                    # pip install -e .[dev]
```

`CONTRIBUTING.md` still references a root `environment.yml` and a `vntyper-dev` env.
Neither exists — use the commands above.

## Commands

| Task | Command |
| --- | --- |
| Full pre-PR gate | `make check-all` |
| Everything CI runs, locally | `make ci-local` |
| Everything the image CI runs | `make ci-local-docker` |
| Format + autofix | `make format` |
| Lint | `make lint` |
| Type check | `make type-check` |
| Fast tests (what CI runs) | `make test-unit` |
| Browser tests (needs a browser engine) | `make test-browser` |
| Inner loop (fail-fast) | `make test-fast` |
| Unit coverage + floor | `make test-unit-cov` |
| Scripts-only coverage proof | `make test-scripts-cov` |
| Coverage of changed lines | `make patch-coverage` |
| Advisory mutation score | `make mutation` |
| Integration tests (needs 1.1 GB data) | `make test-integration` |
| Full gate incl. integration | `make check-full` |
| Docs preview | `make docs-serve` |
| Docs build (CI-equivalent) | `make docs-build` |

Always run `make check-all` before opening a PR. It gates on the **unit** tier only
so it is runnable on a fresh clone; use `make check-full` when you also want the tiers
that need the 1.1 GB Zenodo archive.

**If you changed anything under `.github/workflows/`, `make check-all` is not enough —
run `make ci-local`.** It mirrors `ci-tests.yml` job for job (actionlint, format, lint,
mypy, unit tests + coverage, strict docs build) *and* rebuilds the environment from
scratch through CI's own installer via `make ci-local-uv`. That last part matters:
`ci-local`'s other steps run in whatever environment you already have, so they cannot
catch breakage in how CI *builds* its environment — which is exactly how
`uv pip install --system` shipped green locally and then failed every CI job (see trap
14). `make ci-local-docker` covers the image workflow. Neither target skips silently:
if a required tool is missing they fail and say so.

Run pytest **from the repo root** —
`tests/parametrization.py` opens `tests/test_data_config.json` by relative path at
collection time, so any other CWD breaks collection, including `-m unit`.

## Layout

- `vntyper/cli.py` — argparse entry point; subcommands `pipeline`, `report`, `cohort`,
  `install-references`, `online`. Sole place logging is configured. The parser lives in
  `cli_parser.py`, the subcommand bodies in `cli_handlers.py`, and `report` in
  `cli_report.py`.
- `vntyper/scripts/pipeline.py` — `run_pipeline()`, orchestrates every stage.
- `vntyper/scripts/` — stage modules: `fastq_bam_processing`, `alignment_processing`,
  `kestrel_genotyping`, `variant_parsing`, `scoring`, `confidence_assignment`,
  `motif_processing`, `flagging`, `summary`, `cross_match`, `generate_report`,
  `cohort_summary`, `reference_registry`, `region_utils`, `chromosome_utils`.
- Pure-logic modules extracted from the stage files that grew past the limit. They hold
  the decisions; the file they came from keeps the I/O:
  - `motif_decisions.py` — the motif half/exclusion/GG-allowlist rules, split out of
    `motif_processing.py` (#195).
  - `algorithm_rules.py`, `cohort_categories.py`, `cohort_tables.py`, `cohort_inputs.py`,
    `cohort_exports.py` — the rule table, per-row categorisation, display tables, sample
    discovery and export writers, split out of `cohort_summary.py`.
  - `coverage_presentation.py` — config-driven presentation decisions for coverage QC,
    split out of `screening_summary.py` and `generate_report.py`.
  - `fastp_cutoffs.py` — validation of configured fastp fractions, their paired numeric
    decision values and exact percentage labels, and report-display rounding before icon
    comparison, split out of `report_formatting.py` and `generate_report.py` (#290).
  - `cross_match_presentation.py` — structural cross-match assessability and fixed
    verdict presentation, split out of `generate_report.py`.
  - `advntr_variant_annotations.py` — filesystem-free parsing of repeat-unit identity and
    position from adVNTR state strings, split out of `advntr_genotyping.py`.
  - `cohort_pseudonyms.py` — the pseudonym digest and the validation of the algorithm and
    width it is configured with, split out of `cohort_inputs.py`. Discovery keeps
    `DiscoveredSample`, `discover_sample_directories` and `duplicate_identity`: an
    identity collision is upstream of any digest, so no width fixes it.
  - `alignment_index.py` — resolving an *existing* **BAI** under either name it can carry
    (`<file>.bai`, then `<stem>.bai`), split out of `fastq_bam_processing.py`, which keeps
    every `run_command`/samtools call including the one that builds a missing index. It is
    **not** htslib's resolution order — htslib tries CSI first, and a CSI is ignored here
    on purpose. Non-fast preflight still uses this BAI-only result while protecting
    supplied indexes and builds a fresh co-located BAI. Production indexed recovery now
    lets htslib fetch literal `'*'` reads; it does not call the optional legacy offset
    extractor in `extract_unmapped_from_offset.py`. Widening the resolver remains a
    separate preflight-contract change.
  - `reference_uri_policy.py` — typed, path-free remote-scheme detection with distinct
    parsers for colon-separated `REF_PATH` and complete CRAM header `UR` values, plus
    strict validation of the ambient-resolution boolean.
  - `alignment_preflight_logs.py` — pure planning of every command-log destination
    reachable by a format, candidate count and fast/normal preflight mode.
  - `pipeline_advntr_preflight.py` — pure, typed planning of optional adVNTR enablement
    and model-reference resolution before the pipeline performs model or alignment I/O.
  - `comparator_rules.py` — validation and evaluation of the deliberately small,
    non-executable comparator DSL shared by flagging and cross-match rules.
  - `nomenclature_evidence.py` — source-specific evidence-unit flags and the BAM
    thin-support configuration-key resolver.
  - `nomenclature_presentation.py` — source-specific report flag meanings, tier blockers,
    column help and the Kestrel BAM semantics clarification.

`nomenclature_bam.py` separately owns XD parsing, resolved haplotype-record voting, and
the `BamConsensus` interface; those BAM-facing responsibilities do not belong in the
source-vocabulary helper. `reference_resolution_environment.py` separately owns CRAM-only
process-environment pin/restore I/O.

These eighteen focused modules keep pure decisions independently testable; measure their
current branch coverage rather than assuming a fixed percentage. Put new pure logic there
rather than back in the file it came from.
- `vntyper/modules/{advntr,shark}/` — optional `--extra-modules` stages.
- `docker/app/` — the FastAPI + Celery web service. It is *not* part of the `vntyper`
  package, but it **is** gated: `RUFF_PATHS` covers it and `make type-check` runs
  `mypy vntyper/ docker/app/ scripts/` (#194, #204). The root `scripts/` validation
  instruments are also part of the canonical coverage source — see trap 16.
- `vntyper/dependencies/kestrel/` — vendored JARs, never hand-edit.

## Code style

- Ruff is the only formatter and linter. **Line length 120**, double quotes,
  `target-version = "py310"`. Do not reformat to black's 88.
- `[tool.ruff.lint]` uses an explicit `select`, deliberately **not** `extend-select`.
  Ruff's defaults shift between releases (0.16 turned on ~440 rules at once and took
  CI from green to 740 errors with no code change). Add rules to `select` explicitly;
  never rely on defaults. `BLE001` and `G004` are omitted on purpose — see the
  rationale comment in `pyproject.toml`.
  The reviewed BLE001 policy is 95 normal/102 including suppressions; its executable
  inventory is `scripts/ble001_policy.json` and the policy tests. Not every broad
  handler is a process boundary, so do not globally select or mechanically narrow it.
- mypy is configured in `[tool.mypy]` in `pyproject.toml`, not via Makefile flags.
  Its `python_version`, ruff's `target-version` and `requires-python` are pinned to the
  same floor, and `tests/unit/test_version_consistency.py` fails the build if they drift.
- **Code must run on Python 3.10** (`requires-python = ">=3.10"`, CI matrix 3.10–3.13).
  Your local interpreter is newer, so 3.11+ syntax can pass locally and break CI.
  PEP 604 (`X | None`) and PEP 585 (`list[str]`) are both available at this floor.
- Google-style docstrings (`Args:` / `Returns:` / `Raises:`) on public functions.
- Logging: every module declares `logger = logging.getLogger(__name__)` after its
  imports and calls `logger.info(...)`. Never call `logging.info(...)` on the root
  logger, and never add `basicConfig` or per-module handlers — `setup_logging()` in
  `utils.py` configures the root logger once, and module loggers propagate to it.
  f-strings in log calls are the established style here.
- Errors: no custom exception classes. The convention is `logger.error(msg)` followed
  by `raise ValueError(msg)` / `RuntimeError(msg)` with the same message.
  `except Exception` at stage and process boundaries is intentional, not an oversight.
- Exit codes: **the code that ran to a conclusion uses 0 and 1** — every `sys.exit()` in
  `cli.py` and `cli_handlers.py` is one of those two. **Argument and usage errors exit
  with argparse's 2**, and there are two ways to get there: an unknown subcommand or a
  bad option, rejected by `parse_args()` (`vntyper fastq` → `invalid choice`, exit 2);
  and `parser.error(...)` in `handle_pipeline()` for conflicting or missing inputs
  (`--bam` together with `--cram`, or FASTQ without both mates), which also exits 2.
  Both were verified by running them. So: do not claim the CLI only ever exits 0 or 1 —
  claim that *completed application outcomes* are 0/1 and that argparse owns 2.
  (`parser.error()` never returns, so the `sys.exit(1)` lines immediately after the two
  calls in `cli_handlers.py` are unreachable. Leave them or remove them deliberately;
  do not "fix" the exit code by routing usage errors through them without deciding that
  the contract should change.)
- Type hints are partial. Newer modules (`reference_registry`, `region_utils`,
  `scoring`, `flagging`) are fully annotated — match that when editing them.

## Changing existing code

These three rules apply whenever you touch a file, not just when you add a feature.

**1. Touch a file, add tests for it.** This is not conditional on the module already
having tests, and not something to defer to a follow-up PR.

- Every function you **add** gets unit tests, including its failure paths.
- Every function you **modify** gets a test that would have caught the bug you fixed, or
  that pins the behaviour you changed. Write it *before* the fix where you can — a test
  that fails first and passes after is worth more than one written to match the code.
- Every function you **touch incidentally** (rename, refactor, re-indent) should leave
  with at least the coverage it had.

Tests go in the **unit** tier, so they must be pure: `tmp_path` + `unittest.mock`, no
network, no Docker, no reference data, and `pytestmark = pytest.mark.unit` after the
imports. If a change is only reachable through the integration tier, that is a signal the
logic needs extracting — see rule 3.

Three thresholds enforce this, and they are deliberately different:

| | Where | Behaviour |
| --- | --- | --- |
| **Hard floor: 86** | `fail_under` in `pyproject.toml` | CI **fails** below it. A ratchet — raise it only in a dedicated ratchet change after the rounded integer is sustained by the Python 3.10–3.13 matrix; never lower it to make a build pass. |
| **Patch gate: 80%** | `PATCH_COVERAGE_TARGET` in the `Makefile` | CI **fails a PR** whose *changed lines* fall below it. Not a ratchet, and not an average — it scores your diff and nothing else. |
| **Target: 86%** | `COVERAGE_TARGET` in the `Makefile` | **Warns** only. This is what the project is working towards. |

**The floor is branch-inclusive.** `branch = true` was enabled in `[tool.coverage.run]`
by #196, so `fail_under` is measured against statements *and* branch arcs. That makes the
number strictly harder to move than the statement-only figure the older floors were set
against. The final milestone-4 gate measured **86.14% branch-inclusive across 4,303 unit
tests**. Historically, the #171-#212 suite measured 81.70% after the floor moved from the
earlier 80.24%; statement-only coverage at that earlier branch-coverage rollout was
80.77%.
So deleting `branch = true` *raises* the reported total while covering strictly less —
the ratchet cannot catch that regression, because the number moves the wrong way, and
`tests/unit/test_coverage_gate.py::test_branch_coverage_is_enabled` is the only thing
that can.

The floor and the target now read the same number, and that is deliberate rather than a
merge accident: they measure the same figure for different purposes, so the target stays
a warning and keeps its headroom above the gate. Do not collapse them.

**The patch gate is what makes this rule enforceable.** The hard floor is an average over
thousands of statement and branch units, so a PR can add untested lines without moving
the rounded repository total enough to fail. `make patch-coverage` runs `diff-cover` over
the lines your branch changed and fails below 80%, so an untested new function fails its
own PR regardless of what the repo total is doing. The patch bar remains independently
fixed at 80%; making every PR meet it helps the whole-repository average climb toward the
current 86% target.

It scores against the **merge base**, so commits landing on `main` while your PR is open
are never charged to you. A PR that only deletes code, only touches docs, or changes no
measured Python has nothing to score and passes. The CI job checks out with
`fetch-depth: 0` because finding a merge base needs real history.

**A covered line is not a tested line.** Coverage records that a line *ran*, not that
anything would have failed had it been wrong — and this codebase's characteristic bug is
a silently wrong call, not a crash. `make mutation` measures the difference by breaking
the code on purpose and checking whether the suite notices; `confidence_assignment.py`
once scored 100% line coverage and 21% on that measure. It is **advisory and never
gated**, it runs on the weekly schedule, and the current score with every surviving
mutant named is in `docs/development/mutation-testing.md`. That page *does* hand-classify
the equivalent mutants — `[E]` marks one, with the reason beside it, and the headline
score is quoted both including and excluding them. Read the reason before you write a
test to kill an `[E]`: it cannot be killed, and the honest response is usually to delete
the code that makes it equivalent (the six `.get()` calibration defaults in
`confidence_assignment.py` were removed for exactly that reason in #184).

That page is **machine-generated** by `scripts/mutation_test.py`. Editing the Markdown
directly is silently reverted by the next `make mutation-render`; change the generator and
prove the round trip leaves the committed file byte-identical.

When you add tests to clear the patch gate, write them to kill mutants — assert on the
values, not just that the call returned.

`make test-unit-cov` reports both and prints the candidate edit when local coverage crosses
an integer. Raise the floor only in a dedicated ratchet change after that integer is
sustained by the Python 3.10–3.13 matrix. Never lower the floor to make a build pass — add
the test instead. The gap between the two is real work, and rule 2 explains why it exists.

**When raising the floor, use the number `make test-unit-cov` prints — never the `TOTAL`
column of the coverage table.** Both that column and `coverage report --format=total`
round to an integer, so a true 25.68% displays as `26%`; setting the floor from it makes
CI fail on the very run that produced the number. The gate prints the precise figure and
the exact line to paste.

**2. Keep files under ~650 LOC.** This is a real constraint here, but the reason is not
the line count. **Snapshot, 2026-08-06**, `wc -l` and the branch-inclusive unit-tier
figure at the tip of the `#181–#197` follow-up branch. Re-measure before quoting it — the
table that used to sit here has been wrong twice, once on the numbers and once on the
conclusion drawn from them:

| File | LOC | Unit coverage (branch-inclusive) |
| --- | --- | --- |
| `docker/app/main.py` | 1151 | 88.8% |
| `install_references.py` | 901 | 26.0% |
| `kestrel_genotyping.py` | 877 | 51.7% |
| `pipeline.py` | 721 | 68.5% |
| `fastq_bam_processing.py` | 612 | 60.7% |
| `generate_report.py` | 574 | 64.8% |
| `motif_processing.py` | 534 | 86.3% |
| `docker/app/tasks.py` | 531 | 98.5% |
| `cohort_summary.py` | 456 | 84.2% |
| `region_utils.py` | 334 | 98.0% |
| `confidence_assignment.py` | 194 | 100% |
| `cli.py` | 173 | 85.9% |
| `scoring.py` | 166 | 100% |

**The old claim that "every module over 650 lines is under 30% covered" is false, and so
is its converse.** `docker/app/main.py` is the largest file in the repository at 1151
lines and is 88.8% covered; `pipeline.py` is 721 lines at 68.5%; `kestrel_genotyping.py`
is 877 lines at 51.7%. Only `install_references.py` still fits the old rule.

What the numbers actually say is that **coupling to I/O predicts coverage and size does
not.** `main.py` is large but is a stack of thin FastAPI route handlers, every one
callable through `TestClient` with nothing else running — so it is large *and* testable.
`install_references.py` is the same size and downloads, unpacks and checksums files, so
its logic cannot be reached without a network and a filesystem, and it sits at 26% for
that reason rather than for its length. `docker/app/tasks.py` is 531 lines of Celery task
at 98.5%; `cohort_summary.py` was 911 lines at 38% and is 456 at 84.2% — not because 456
is under some threshold, but because the split moved the decisions into
`algorithm_rules`/`cohort_categories`/`cohort_tables`/`cohort_inputs`/`cohort_exports`, all
at 100%, and left the matplotlib and Jinja2 calls behind.

So keep the ~650 line guideline, but keep it for the right reason: **a file grows past
650 lines here by fusing I/O, orchestration and pure logic into functions that cannot be
called without a filesystem.** Size is the symptom you can measure cheaply; fusion is the
cause. A large file that is already decoupled does not need splitting, and a 300-line
file that shells out in the middle of its only public function is still untestable.
Splitting is how the coverage gets written — see rule 3 for which part to split.

**3. Refactor the part you touch.** When you edit inside a file over the limit, extract
the region you are working on into a focused module rather than growing the file further.
Pull the *pure* logic out first (parsing, scoring, filtering, formatting) and leave the
I/O behind — the pure part is the part that gets a test. `scoring.py`,
`confidence_assignment.py` and `region_utils.py` are the shape to copy. Do not attempt a
whole-file rewrite as a side quest: split out the region under change, test it, move on.

Still over the limit, worst first (2026-08-06): `docker/app/main.py` (1151),
`install_references.py` (901), `kestrel_genotyping.py` (877), `pipeline.py` (721). The
first two are the ones worth splitting next, and only one of them is a coverage problem.
`generate_report.py` (574), `cohort_summary.py` (456), `cli.py` (173) and
`tests/integration/test_pipeline_integration.py` (243) have all come back under the
limit.

## Testing

- `pytest.ini` at the repo root is the live config. The `[tool.pytest.ini_options]`
  block in `pyproject.toml` is **dead** — pytest.ini takes precedence. Edit `pytest.ini`.
- Markers: `unit`, `integration`, `docker`, `browser`, `golden`, `smoke`, `slow`. `smoke` is the fast image
  tier (`make test-docker-smoke`, ~1 s): it asserts the built image's *structure* — that
  every reference path `config.json` declares actually exists in it, that the three conda
  envs and their interpreters are present, that adVNTR imports under Python 2.7, and that
  the image stays under its size budget. It needs a Docker daemon but **no test data**,
  and it derives its assertions from the config inside the container rather than from a
  hardcoded list, so it cannot drift. Smoke tests carry **only** the `smoke` marker —
  adding `docker` too would make the slow tier re-run them.
- **`browser` is one of five directory tiers (`tests/browser`, `make test-browser`) and nothing gates
  on it.** The report is an HTML file people open in a browser, and its DataTables flag
  filter, its client-side rounding and its three CDN `<script>` tags do not exist until a
  JavaScript engine has run — so `tests/unit/test_generate_report.py` can assert on
  everything the report *is* and nothing about what it *does*. The tier renders one real
  report through `generate_summary_report` and opens it twice: once with the network
  reachable, once with every non-`file://` request route-blocked in the browser rather
  than trusted to be unreachable. It needs `pytest-playwright` (declared in the `dev`
  extra) plus a browser binary (`playwright install chromium`) and the network, so it is
  deliberately outside `make check-all`, outside CI and outside the coverage floor —
  `fail_under` stays a unit-tier figure. `test_the_report_reads_identically_online_and_offline`
  is the standing acceptance check for issue #242 and is **green**: the report now reads
  the same with and without a network. It was red when the tier landed (1 of 3 Kestrel
  rows online against 3 of 3 offline, at different precision), and it is kept so that the
  next change cannot reintroduce that. `make test-browser` is expected to pass.
  **Two marker rules that look inconsistent and are not.** `make test-unit` needs no
  `-m "not browser"`: it selects `-m unit tests/unit`, by path *and* by marker, so it
  cannot reach the tier, and the predicate would wrongly imply browser tests carry the
  `unit` marker. `make test`, `test-cov`, `test-quiet` and `test-verbose` **do** carry it:
  they invoke a bare `pytest`, `testpaths = tests` then collects `tests/browser`, and
  `make all` depends on `test` — so without the predicate the convenience gate demands a
  browser binary from everyone. State the exclusion where it does something.
  The online pass must *prove* it was online (a working jQuery/DataTables in the page,
  and no error statuses), because a machine with no network, or a CDN answering with
  something unusable, otherwise turns it into a second offline pass and the comparison
  compares a report with itself.
- **`golden` is the fifth directory tier** (`tests/golden`,
  `pytest -m golden tests/golden`). It compares against a known-truth simulated cohort
  supplied through `VNTYPER_SIM_ROOT` and `VNTYPER_ADVNTR_ROOT`, skips without those
  inputs, and deliberately remains outside `make check-all`.
- **Every new unit test file must declare `pytestmark = pytest.mark.unit`.** CI runs
  `pytest -m unit`, so an unmarked file silently never runs. This is enforced by
  `tests/unit/test_marker_hygiene.py`, which fails the build naming the offending file;
  `--strict-markers` additionally turns a *misspelled* marker into an error.
- Unit tests must stay pure: `tmp_path` + `unittest.mock`, no network, no reference data.
- Integration tests require both the ~1.1 GB Zenodo archive and the pinned host
  `reference/alignment/chr1.hg19.fa`. Install that reference with
  `vntyper install-references -d reference --references hg19`; the test-data download
  does not install it. Docker tests use the pinned reference bundled in the image and
  require only the archive on the host. The test bootstrap MD5s all 114 `file_resources`
  entries, and one missing file triggers a full re-download.
  **The test-data bootstrap has no skip fallback**: `tests/conftest.py` and
  `tests/support/data_utils.py` call `pytest.exit(..., returncode=1)` on missing config,
  a missing Kestrel JAR or a failed download, which ends the whole session — unit tests
  included — rather than skipping the tier that needs the data. Set
  `VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1` to fail fast instead of downloading.
  (Skips do exist elsewhere and are fine: `tests/docker/test_image_structure.py` carries
  a module-level `pytest.mark.skipif` on the Docker daemon plus a `pytest.skip` when the
  image is not present locally, and `test_version_consistency.py`,
  `test_workflow_consistency.py` and `test_integration_tier_hygiene.py` skip on files
  that are absent from a partial checkout. The old blanket claim that there are no
  `pytest.skip` guards *anywhere* was simply false — do not restore it.)
- Correct marker composition is `-m "docker and not slow"`; repeated `-m` flags override
  rather than combine.
- The Docker tier runs at three depths, chosen by how much signal each buys for its
  runtime: PRs get `test-docker-quick` (4 tests, ~10 s), pushes to `main` get
  `test-docker-fast` (everything except `slow`, ~1.5 min locally), and a nightly
  schedule plus every `workflow_dispatch` runs `test-docker` including adVNTR.
  adVNTR is off the merge path on purpose: one test costing 15-25 min on a
  2-core runner - more than the rest of the pipeline combined - for an optional
  module. Both adVNTR tests carry `@pytest.mark.timeout(ADVNTR_TIMEOUT_SECONDS)`,
  overriding the global 600 s; that global timeout is right for everything else and CI
  hardware is several times slower than a workstation. The value and the reasoning
  behind it live in `tests/support/orchestration.py` — the module both adVNTR tests already
  share — so recalibrating for new runner hardware is a one-line change there.

## Git and PRs

- Branch from `main` (there is no `develop` despite the CI trigger list). Observed
  naming: `fix/issue-145-flagging-before-selection`, `test/issue-156-gdp-inflation-guard`,
  `docs/mkdocs-site`. Use `<type>/issue-<N>-<slug>` or `<type>/<slug>`.
- Conventional Commits: `type(scope): lowercase subject` under ~72 chars, with
  `Closes #N` in the footer. Scopes are module or area names (`kestrel`, `flagging`,
  `config`, `ci`, `docker`).
- PRs are **not** squashed — every commit becomes permanent history. Keep them clean.
- Copilot and Sourcery AI review PRs. The established pattern is a follow-up commit
  addressing their comments, not a force-push.
- CI gates on PRs: `make lint`, `make type-check-all`, `make test-unit` across 3.10–3.13,
  plus a Docker image build and quick Docker tests. PR path filters remain, while every
  push to `main` is treated as substantive and runs both required-check workflows.
  Formatting is *not* gated in CI — run `make format-check` yourself.

## Release workflow

Production release policy lives on the default branch in
`.github/workflows/publish-pypi.yml`. An operator first creates the existing tag outside
the workflow, then sends an authenticated `repository_dispatch` of type `vntyper_release` whose
`client_payload.tag` is that strict `vX.Y.Z` tag. There is deliberately no production
`push.tags` trigger: the default-branch controller validates the tagged commit, but the
workflow never creates or moves a tag and historical workflow code cannot select the
current production policy. A manual invocation such as
`gh workflow run publish-pypi.yml -f tag=vX.Y.Z` validates an existing tag as a dry run
and performs no production writes.

The release accepts only a tag on `main` ancestry after these exact ten checks succeed
for the full tagged SHA: `Lint (Ruff)`, `Type Check (mypy)`, `Unit Tests (Python 3.10)`,
`Unit Tests (Python 3.11)`, `Unit Tests (Python 3.12)`, `Unit Tests (Python 3.13)`,
`Docs build (strict)`, `CI Success`, `Build and test image`, and `Docker Success`.
For each exact name and SHA, the newest GitHub Actions check is authoritative even when it
came from a scheduled or manual rerun. A newer failure intentionally supersedes an older
success; recover by rerunning that exact SHA to a newer green check, never by selecting the
stale success.
The Docker `sha-<7>` tag is only a convenience reference: contract-v1 evidence binds
the full SHA and immutable digest, and the image labels
`org.opencontainers.image.revision` and `org.opencontainers.image.version` must match
that full revision and release version. A proven short-prefix collision continues from
the evidence digest, while ambiguous short-tag drift fails closed.
The push of `main` is the only event allowed to publish the application `main` and
short-SHA tags and their evidence; scheduled and manual Docker validation never publish application tags.

Promotion copies the evidence-verified digest, never rebuilds it. Exact `vX.Y.Z` and
`X.Y.Z` aliases are immutable; floating `X.Y`, `X`, and `latest` advance monotonically
and never downgrade. `main` remains rolling/unreleased and is never a release alias.
For floating aliases, legacy rolling `main` is the one recognized pre-semver version
state and advances to the verified release digest. A missing or unrecognized version label
fails closed before any alias write and requires the operator to repair or remove that
floating alias before retrying.
GHCR mutation is serialized by the fixed `vntyper-ghcr-release-promotion` group. If a
pending promotion is canceled or superseded before it acquires that lock, rerun it; if
Docker evidence is missing, rerun the existing exact-SHA Docker Build run. A partial
promotion is idempotent, and a PyPI retry uses `skip-existing` while reporting whether
the version was newly published or already existed.

PyPI publication uses Trusted Publishing in the protected `pypi` environment with only
`id-token: write`; there is no package token in the current controller.
The environment must be reviewer-free, have no wait timer, use custom branch policies for the
exact branch `main`, and have no custom deployment-protection rules. The controller preflights
exactly the environment, deployment-branch-policy, and custom-deployment-protection-rule API
responses and fails before package or registry writes when they drift; follow #236 to repair the
environment rather than bypassing it. The absence of environment secrets is separately verified
live administrator state, not controller-enforced. OIDC is the only publisher: never reintroduce
`PYPI_API_TOKEN`. `PYPI_API_TOKEN` has been deleted after the green OIDC releases, so historical
tagged commits cannot retrieve the obsolete credential.

```text
phase | job | permissions | retry/recovery
validation | validate-release | actions: read, contents: read | fix identity or ancestry, then rerun the existing tag
gates | wait-for-release-gates | actions: read, checks: read, contents: read, packages: read | rerun missing exact-SHA checks or the existing Docker Build run
build | build-package | contents: read | rebuild the exact candidate wheel and sdist artifact
promotion | promote-ghcr | contents: read, packages: write | rerun; verified aliases no-op and partial progress converges
publish | publish-pypi | id-token: write | rerun safely; skip-existing distinguishes an existing release
summary | release-summary | none | always records success, failure, skipped jobs, and unavailable partial outputs
```

## Traps

1. **Config is loaded at import time.** `kestrel_genotyping.py`, `advntr_genotyping.py`
   and `shark_filtering.py` read their JSON into module globals on import, and
   `run_kestrel()` reads the module-level `kestrel_config` rather than its `config`
   argument. `--config-path` cannot override them; tests must patch the module global.
2. **`--config-path` replaces the whole config, it does not merge.** Missing keys raise
   `KeyError` deep in the pipeline (`config["tools"]["java_path"]`, no `.get`).
3. **Both configuration-driven rule consumers share one deliberately small comparator
   security boundary.** `flagging.py` and `cross_match.py` validate structured rules
   through `comparator_rules.py` before processing any row. Grep the JSON configs before
   renaming any column.
   - The only operators are `eq`, `lt`, `in`, and `casefold_eq`; predicates compose through
     explicit non-empty `all`/`any` nodes and single-child `not` nodes, with nesting capped
     at 32 boolean nodes. Null is false, types are not coerced, booleans are accepted only
     by same-family `eq`, and configured or non-null runtime numbers must be finite.
     Missing columns, malformed rules, unsafe flag names, and incompatible row values
     log and raise `ValueError` instead of silently producing "not flagged" or "no match".
   - Migration is intentionally finite. Flagging accepts only the five flag-name-scoped,
     byte-exact expressions shipped immediately before #286; cross-match accepts only its
     one byte-exact historical expression. Modified whitespace, custom expressions, and
     an expression assigned to the wrong flag are rejected. There is no general-purpose
     expression parser and no `eval`/`exec` path.
   - Packaged configs still load at their existing import-time locations. Validation
     occurs when the consuming flagging or cross-match stage begins, before copying or
     iterating a DataFrame, rather than at module import.
4. **Stages mark, they do not filter.** Kestrel stages append **six** boolean columns
   (`is_frameshift`, `is_valid_frameshift`, `depth_confidence_pass`, `alt_filter_pass`,
   `motif_filter_pass`, `flag_filter_pass`); `filter_final_dataframe()` ANDs them at the
   end. Preserve that contract — dropping rows early breaks `kestrel_pre_result.tsv`
   debuggability.

   `flag_filter_pass` arrived with #174 and behaves differently from the other five in one
   way worth knowing: it is written **unconditionally** in step 6.5, outside the `if
   flagging_rules or duplicates_config` block that guards `add_flags`. That is deliberate.
   The gate contract raises on a missing required column, and a run with no flagging rules
   configured legitimately produces no `Flag` column at all — so a conditional gate would
   turn a config choice into an aborted run. It reads its artifact list from
   `kestrel_config.json`'s `artifact_flags`; **no flag name is written inline in Python**,
   and emptying that list restores the pre-#174 behaviour with no code change.

   **Two tripwires fire on any change to the gate list, and they are not in the same
   file.** `tests/unit/test_kestrel_filtering.py` reads `kestrel_genotyping.py` as *source
   text* and asserts the exact count; `tests/builders.py`'s `STAGE_COLUMNS["flagged"]` and
   `["final"]` feed `kestrel_stage_frame()` to the real `filter_final_dataframe`, so a gate
   missing from those tuples raises `ValueError` rather than failing an assertion. Change
   both, deliberately, in the same commit.
5. **Summary step names are string literals.** `"Kestrel Genotyping"`,
   `"adVNTR Genotyping"`, `"Coverage Calculation"`, `"BAM Header Parsing"`,
   `"Cross-Match Variant Comparison"` are matched exactly by `generate_report.py`,
   `cohort_summary.py` and `cross_match.py`. All `parsed_result` values are strings.
6. **Conda env names are load-bearing.** `config.json` shells out to
   `mamba run -n envadvntr advntr` and `mamba run -n shark_env shark`. `envadvntr` is
   Python 2.7 and can never be merged into the main env.
7. **All tool and reference paths in `config.json` are relative to the process CWD.**
   `pipeline.py` pins `project_root = os.getcwd()` and threads it as `cwd=` into Java
   and samtools calls — do not remove that plumbing.
8. **Kestrel is pinned to 1.0.1.** Scoring thresholds are calibrated to it; 1.0.2 output
   differs. Do not upgrade the JAR.
9. **`run_command` uses `shell=True`** deliberately, for process substitution in the CRAM
   unmapped-read path. Converting it to `shell=False` breaks that branch.
10. **The image is split in two, and the halves are bound by a content hash.**
    `docker/Dockerfile.base` holds the conda envs, adVNTR and the reference genomes and
    is rebuilt only by `docker-base.yml` when `conda/**`, `requirements-web.txt`,
    `install_references*`, `reference_bundle.py` or `dependencies/advntr/**` change.
    `docker/Dockerfile` holds only the application and builds on top in ~3 min. Both
    workflows compute the same `hashFiles()` tag, so a base input can never be changed
    and then silently built against a stale base. What happens when the tag is missing
    depends on **who is running the workflow**, and the two paths are very different:
    - **Same-repository run** (push, or a PR from a branch in this repo): the
      `base-status` job records `missing=true`, the `build-base` job calls the reusable
      `./.github/workflows/docker-base.yml` and publishes the base, and `build-and-test`
      waits on it. No manual step, no failure — just a slow run.
    - **Fork PR**: a fork's `GITHUB_TOKEN` is read-only and cannot push to GHCR, so
      `base-status` fails the run with an explicit error. That is the only case where
      "a maintainer must run `docker-base.yml` via `workflow_dispatch`, or merge the base
      change first" applies.

    Locally, neither path exists — build the base yourself:
    `make docker-build-base && make docker-build DOCKER_BASE_IMAGE=vntyper-base:local`.
    `make ci-local-docker` refuses to run against the published `:latest` base when you
    have edited a base input, for the same reason.

    (The split replaced a `RUN git clone` of GitHub `main`, which meant PR builds never
    tested PR code and a cached layer could serve an arbitrarily stale checkout.)

    `reference/**` is **not** a base input, and has not been since milestone 5: reference
    data comes from a published, checksummed bundle in `berntpopp/vntyper-data` rather
    than from tracked files, and `reference/` holds only `README.md`, `pseudonymize.py`
    and `pseudonymize_config.json`. Changing a reference byte now means publishing a new
    bundle release and updating the `asset`/`asset_sha256` pins in
    `install_references_config.json` — which *is* still a base input, so that pin change
    is what triggers the rebuild instead.

    `reference_download.py` owns bounded, atomic reference transport and is imported by
    `install_references.py`; the refs-stage COPY list, all base-image hash expressions,
    the Docker push-path filter and `BASE_INPUTS` must include it together.
    `reference_integrity.py` separately owns the digest-gated reuse decision: only a
    checksum mismatch on a pre-existing file is re-downloaded, exactly once. Missing or
    contradictory digest configuration and a mismatch on fresh bytes fail immediately.
    It is another refs-stage dependency, so the same COPY/hash/input mirrors include it.

    The `refs` stage installs **all six** physical BWA references (`hg19`, `hg38`,
    `GRCh37`, `GRCh38`, `hg19_ensembl`, `hg38_ensembl`), not just the two UCSC ones a bare
    `install-references` run defaults to, and its `config.json` is the repo's file
    unmodified — no `--config-path` is passed at build time. That matters together: every
    one of the six `bwa_reference_*` keys `config.json` ships as a real relative path
    therefore resolves to a file that actually exists in the image, so
    `reference_registry`'s UCSC-family fallback (`bwa_reference_hg38_ensembl` missing,
    falling back to `bwa_reference_hg38`) should never fire inside the container for any
    accepted `--reference-assembly` label — every label gets its own physical file, not a
    UCSC stand-in. `tests/docker/test_image_structure.py::test_every_declared_reference_exists`
    is what enforces this: it reads `config.json`'s declared paths from *inside* the image
    and fails if any of them is missing, so it cannot drift from what the image actually
    ships. It is also the reason `MAX_IMAGE_BYTES` was raised from 6 GiB to 9 GiB — the
    four newly-shipped genomes add roughly 2.57 GiB uncompressed.
11. **The report's presentation logic lives outside `generate_report.py`.**
    `screening_summary.py` owns the screening state and the `report_config.json` rule
    table; `fastp_cutoffs.py` owns validated fastp decision values, their paired labels
    and report-display rounding; `report_formatting.py` owns the icons, the column
    projections and the IGV fragment splicing. Put new pure logic in the focused owner,
    not back in `generate_report.py`.
    Two rules that are easy to break: emphasis in the report comes from the computed
    state (`screening_state.emphasis`, and `cross_match_is_positive` for the
    cross-match box), never from searching the message text, and the
    Jinja2 environment autoescapes — anything you mark `|safe` must be a fragment
    VNtyper built, not a value read from a sample. `tests/unit/test_generate_report.py`
    enforces both. (`vntyper report` itself was broken until #179 — the handler passed
    arguments `generate_summary_report()` does not accept — and now lives in
    `cli_report.py` with its call contract pinned by a signature-binding spy.)
12. **Version lives in three places**: `vntyper/version.py` (authoritative),
    `CITATION.cff`, and `docs/about/changelog.md`. The default-branch release controller
    requires all three to equal the existing strict tag, requires the tag commit to be
    on `main`, waits for the exact ten full-SHA checks, and validates contract-v1 Docker
    evidence before any write. The `sha-<7>` image is accepted only when its digest and
    full `org.opencontainers.image.revision`/`org.opencontainers.image.version` labels
    match that evidence. A short-prefix collision proven by those identities continues
    from the evidence digest; ambiguous drift fails closed. Promotion copies the immutable
    digest to exact and monotonic floating aliases. Reruns converge after
    partial GHCR progress and PyPI uses `skip-existing`, but an identity conflict never
    overwrites an existing exact alias.
13. **The package and the image must declare the same versions, and a test enforces it.**
    `conda/environment_vntyper.yml` is what the Docker image *runs*; `pyproject.toml` is
    what the PyPI package *declares*. When they disagree, `pip install .` inside the
    image resolves something different from the recipe — that is not hypothetical:
    pyproject pinned `numpy>=1.26.0,<2.0.0` while the env installed `numpy=2.0.2`, so the
    published image silently ran 1.26.4 on a different numerical stack.
    `tests/unit/test_version_consistency.py` fails the build if
    a conda pin violates a pyproject specifier, if the image's Python is not in the CI
    matrix, if the classifiers and matrix disagree, if ruff's `target-version` or mypy's
    `python_version` drifts from `requires-python`, or if a binary the image smoke tier
    requires is missing from the environment. Change versions in **both** files.
14. **CI installs with `uv` into an explicit venv, never `--system`.** GitHub's
    `ubuntu-24.04` image ships a PEP 668 "externally managed" interpreter, so
    `uv pip install --system` fails outright with
    `error: The interpreter at /usr is externally managed`. Every job therefore runs
    `uv venv`, exports `VIRTUAL_ENV`, and appends `.venv/bin` to `$GITHUB_PATH` before
    installing. Keep that three-line block if you add a job; `make ci-local-uv`
    reproduces it locally, and it is the only local check that exercises the install
    layer at all.
15. **The Dockerfile installs the package with `pip install --no-deps`, deliberately.**
    Every runtime dependency is conda-managed by `conda/environment_vntyper.yml`; without
    `--no-deps`, pip resolves the same requirements again from PyPI and lays a second
    numpy over the conda-managed one, mixing provenance on a numerical stack. The two
    files no longer *conflict* (`pyproject.toml` pins `numpy>=1.26.0`, the env installs
    `numpy=2.0.2`, which satisfies it) — but they are still resolved by two different
    solvers, so `--no-deps` stays.
16. **`scripts/` is in both runtime quality gates.** `RUFF_PATHS` covers
    `vntyper/ docker/app/ tests/ scripts/ docs/`, and `make type-check` runs
    `mypy vntyper/ docker/app/ scripts/`; the final no-cache run checked 134 source files.
    `make type-check-all` then adds the deliberately separate `mypy vntyper/ tests/`
    pass, which checked 295 test sources. Do not claim `scripts/` is part of that tests
    invocation or collapse the two.

    Coverage uses ``source = [`vntyper`, `docker/app`, `scripts`]``, with root `scripts`
    appended only after the scripts aggregate exceeded the separate 88% threshold. The
    2026-08-11 verification measured 39 Python files at 7,249 of 7,781 measured
    units, or 93.16% aggregate scripts-only branch-inclusive coverage, and 17,310 of
    19,392 measured units, or 89.26% combined branch-inclusive coverage. All 5,432 unit
    tests passed with no skips, in the maintained Python 3.12.13 environment. **Do not
    quote a fixed warning count here** - it is not pinned by any gate and drifts with
    dependency versions, so a number recorded on one date reads as a false regression on
    another. Measure your own baseline (`pytest -m unit`) before a change and treat a
    *wider* count than that freshly measured baseline as the signal, never a count quoted
    in this file. (For calibration only, not as a target: on 2026-08-11, `main` measured
    5,669 passed / 464 warnings on this machine, and this branch measured the same 464
    warnings over more tests - both already superseding whatever fixed count was quoted
    here before, which is exactly the drift this note now warns against instead of
    restating a new one.)
    Adding `scripts/bundle_release.py` (milestone 5, the reference bundle
    builder) took the directory to 40 Python files and `make test-scripts-cov` to
    7,806 of 8,340 measured units, or 93.60%, over 5,665 passing unit tests.
    Adding `scripts/golden_cohort/waiver.py` (#262, the gate's declared-delta waiver
    policy, extracted out of `compare.py`) took it to all 43 Python files and
    `make test-scripts-cov` to 7,869 of 8,392 measured units, or 93.77%, over 6,299
    passing unit tests. Adding `scripts/integration_compatibility_observations.py`
    (#293, the version-selected append-only report-observation policy) took it to all
    42 Python files and `make test-scripts-cov` to 94.00% over 8,138 collected unit
    tests. **A new file under `scripts/` must update this sentence** -
    `tests/unit/test_coverage_gate.py::test_contributor_docs_match_the_scripts_quality_scope`
    counts root `scripts/**/*.py` and fails until it does, which is the tripwire working,
    not a flaky test. Package modules under `vntyper/scripts/` do not change that root-only
    count; adding `vntyper/scripts/reference_download.py` therefore leaves it at 42.
    `ci-local`'s clean Python 3.13.6 rebuild and the Python 3.10–3.13
    GitHub matrix remain the authoritative cross-version gates. These figures do not
    change the independent gate semantics:
    `[tool.coverage.report].fail_under = 86` is the hard floor,
    `COVERAGE_TARGET ?= 86` is advisory, `PATCH_COVERAGE_TARGET ?= 80` scores changed
    lines, and `SCRIPTS_COVERAGE_TARGET ?= 88` is the isolated pre-source proof. The
    existing `exclude_also` entry ignores only mechanical direct-execution bootstrap
    guards; callable `main(...)` functions, exit policy and substantive CLI branches
    remain measured. Add no scripts omit entry to improve any of these numbers.
17. **Real-data successes are append-only compatibility contracts.**
    `tests/compatibility/real_success_baseline.json` is independent of the mutable
    declarations in `tests/test_data_config.json`. Identity is `(suite, test_name)`, so
    one input may legitimately occur in more than one suite. `make
    check-integration-compatibility` validates baseline-to-live and qualifying
    live-to-baseline coverage in the same process, then compares the complete historical
    row and every versioned observation set with the explicit Git event base. Schema 2
    preserves the schema-1 `contracts` list as an exact ordered prefix and adds an
    append-only `observation_sets` chain selected from the exact package version in
    `vntyper/version.py`, never from mutable test JSON. Observation sets may override
    only a non-empty report assertion for an existing `(suite, test_name)`; their full
    provenance commit must resolve, be ancestral to the checkout, and carry the same
    package version. Missing, stale, ambiguous, reordered, or mutated history fails
    closed. Never normalize a regression by changing exit
    0 to exit 1, renaming or deleting a case or row, or weakening its input/reference
    digest, threads, log level, CLI options, modules, routing counts or selection,
    required artifacts, archive state, value fields, centers, or tolerances. A new
    qualifying real success must add a complete row. Pull requests use
    `origin/${{ github.base_ref }}` and pushes use `${{ github.event.before }}` with full
    history; missing, empty, all-zero, shallow, or unreachable bases fail closed. The
    only no-base bootstrap is the checker's pinned ten historical successes
    reconstructed from the commit before `52c4146596fef2d1e2402991fbab062ba8021889`;
    there is no exception or free-text waiver.

18. **A characterisation fingerprint must exclude a dependency's scaffolding, not just
    its payload.** `tests/unit/test_cohort_summary_oracle.py` hashes the assembled cohort
    report. It excluded the embedded base64 PNG *payloads* on the stated grounds that "a
    hash that turns red on `pip install -U plotly` stops being read" — and then pinned
    their **count**, which is the same quantity one level up. Measured 2026-08-26: plotly
    6.9.0 emits four 114-character scaffolding images and 7.0.0 emits two, so the oracle
    went red on `main` itself the first time CI resolved dependencies afresh, with no
    change to this project. The two `count=`/`nonempty=` lines were the *entire*
    difference between the two documents.

    Two rules follow. **Before adding a section to that document, ask whether a
    dependency bump alone can move it** — a count, a length, an ordering and an
    attribute set all can, not only a payload. And **a green run of that oracle proves
    nothing about a fresh dependency resolution**: reproduce a CI-only fingerprint
    failure with `uv venv --python 3.13 && uv pip install -e ".[dev]"` in a worktree at
    the base commit before assuming the branch caused it. The conda `vntyper` env pins
    older packages than CI resolves.

    The replacement states the property rather than the quantity — no embedded image is
    left without a payload — and is deliberately vacuous if plotly stops emitting any,
    because their existence is plotly's business and not this project's.

## Never

- Never commit into `tests/data/`, `reference/`, `out/`, `download/`, `vntyper.egg-info/`
  or the local `adVNTR/` clone — all are generated or downloaded.
- Never hand-edit `vntyper/dependencies/kestrel/*.jar` or anything in `vntyper.egg-info/`.
- Never add a **published** page under `docs/` without registering it in `mkdocs.yml`
  `nav:`. Two corrections to the older wording of this rule, both measured on
  mkdocs 1.6.1:
  - An unregistered page is logged at **INFO** ("The following pages exist in the docs
    directory, but are not included in the nav configuration"), not WARNING, so it does
    **not** fail `--strict` on its own. A dangling nav entry or a broken internal link
    does.
  - What actually breaks the build on an unregistered page is usually the `macros`
    plugin, which Jinja-templates every page: prose containing `{#...}` (a list of GitHub
    issues such as `{#209, #178}` is enough) fails as an unterminated comment.
  **`docs/` is strictly the published site — put nothing else there.** Contributor
  working documents (design specs, implementation plans, execution prompts, milestone
  gate artifacts) belong in an untracked `.planning/` workspace, conventionally
  `specs/`, `plans/`, `prompts/` and `milestones/` beneath it. That directory is
  gitignored, so **a fresh clone does not have one** and nothing may depend on its
  contents. The #295 program has one bounded, published exception for these exact paths:
  `docs/superpowers/specs/2026-08-31-kestrel-bam-evidence-semantics-design.md` and
  `docs/superpowers/plans/2026-08-31-kestrel-bam-evidence-semantics.md`, plus the umbrella
  `docs/superpowers/specs/2026-08-31-issue-295-completion-program-design.md` and
  `docs/superpowers/plans/2026-08-31-issue-295-completion-program.md`.
  No other planning artifact under `docs/` is allowed. All listed pages are registered in
  `mkdocs.yml` and tested against raw Jinja opening delimiters; this exception does not
  restore `exclude_docs:` or permit an unpublished page beneath `docs/`.

  Because `.planning/` is untracked, **the durable record of why a change was made is the
  commit message, the pull request and the issue** — not a file only the author has.
- Never claim tests pass without showing the command output.
