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
| Inner loop (fail-fast) | `make test-fast` |
| Unit coverage + floor | `make test-unit-cov` |
| Integration tests (needs 1.1 GB data) | `make test-integration` |
| Full gate incl. integration | `make check-full` |
| Docs preview | `make docs-serve` |
| Docs build (CI-equivalent) | `make docs-build` |

Always run `make check-all` before opening a PR. It gates on the **unit** tier only
(~0.5 s) so it is runnable on a fresh clone; use `make check-full` when you also want the
tiers that need the 1.1 GB Zenodo archive.

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

- `vntyper/cli.py` — argparse CLI; subcommands `pipeline`, `report`, `cohort`,
  `install-references`, `online`. Sole place logging is configured.
- `vntyper/scripts/pipeline.py` — `run_pipeline()`, orchestrates every stage.
- `vntyper/scripts/` — stage modules: `fastq_bam_processing`, `alignment_processing`,
  `kestrel_genotyping`, `variant_parsing`, `scoring`, `confidence_assignment`,
  `motif_processing`, `flagging`, `summary`, `cross_match`, `generate_report`,
  `cohort_summary`, `reference_registry`, `region_utils`, `chromosome_utils`.
- `vntyper/modules/{advntr,shark}/` — optional `--extra-modules` stages.
- `docker/app/` — the FastAPI + Celery web service. It is *not* part of the `vntyper`
  package and is not covered by `make lint` / `make type-check`.
- `vntyper/dependencies/kestrel/` — vendored JARs, never hand-edit.

## Code style

- Ruff is the only formatter and linter. **Line length 120**, double quotes,
  `target-version = "py39"`. Do not reformat to black's 88.
- `[tool.ruff.lint]` uses an explicit `select`, deliberately **not** `extend-select`.
  Ruff's defaults shift between releases (0.16 turned on ~440 rules at once and took
  CI from green to 740 errors with no code change). Add rules to `select` explicitly;
  never rely on defaults. `BLE001` and `G004` are omitted on purpose — see the
  rationale comment in `pyproject.toml`.
- mypy is configured in `[tool.mypy]` in `pyproject.toml`, not via Makefile flags.
  It targets 3.10 because mypy 2.x dropped 3.9; ruff's `target-version` still guards
  the 3.9 runtime floor.
- Code must run on Python 3.9 (CI matrix: 3.9–3.12). Your local interpreter is newer,
  so 3.10+ syntax can pass locally and break CI.
- Google-style docstrings (`Args:` / `Returns:` / `Raises:`) on public functions.
- Logging: every module declares `logger = logging.getLogger(__name__)` after its
  imports and calls `logger.info(...)`. Never call `logging.info(...)` on the root
  logger, and never add `basicConfig` or per-module handlers — `setup_logging()` in
  `utils.py` configures the root logger once, and module loggers propagate to it.
  f-strings in log calls are the established style here.
- Errors: no custom exception classes. The convention is `logger.error(msg)` followed
  by `raise ValueError(msg)` / `RuntimeError(msg)` with the same message. Only exit
  codes 0 and 1 are used. `except Exception` at stage and process boundaries is
  intentional, not an oversight.
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

Two thresholds enforce this, and they are deliberately different:

| | Where | Behaviour |
| --- | --- | --- |
| **Hard floor** | `fail_under` in `pyproject.toml` | CI **fails** below it. A ratchet — raise it when coverage climbs, never lower it to make a build pass. |
| **Target: 80%** | `COVERAGE_TARGET` in the `Makefile` | **Warns** only. This is what the project is working towards. |

`make test-unit-cov` reports both and prints the exact edit to raise the floor whenever
coverage exceeds it. Never lower the floor to make a build pass — add the test instead.
The gap between the two is real work, and rule 2 explains why it exists.

**When raising the floor, use the number `make test-unit-cov` prints — never the `TOTAL`
column of the coverage table.** Both that column and `coverage report --format=total`
round to an integer, so a true 25.68% displays as `26%`; setting the floor from it makes
CI fail on the very run that produced the number. The gate prints the precise figure and
the exact line to paste.

**2. Keep files under ~650 LOC.** This is a real constraint here, not style preference.
Measured on this repo, the correlation is total:

| File | LOC | Unit coverage |
| --- | --- | --- |
| `cohort_summary.py` | 856 | 0% |
| `cli.py` | 700 | 0% |
| `generate_report.py` | 861 | 4% |
| `pipeline.py` | 735 | 10% |
| `install_references.py` | 901 | 24% |
| `kestrel_genotyping.py` | 835 | 28% |
| `region_utils.py` | 246 | 98% |
| `scoring.py` | 176 | 100% |
| `confidence_assignment.py` | 190 | 100% |

Every module over 650 lines is under 30% covered; every well-covered module is under 650.
Oversized files here are oversized because they fuse I/O, orchestration and pure logic
into functions that cannot be called without a filesystem — which is exactly what makes
them untestable. Splitting them is how the coverage gets written.

**3. Refactor the part you touch.** When you edit inside a file over the limit, extract
the region you are working on into a focused module rather than growing the file further.
Pull the *pure* logic out first (parsing, scoring, filtering, formatting) and leave the
I/O behind — the pure part is the part that gets a test. `scoring.py`,
`confidence_assignment.py` and `region_utils.py` are the shape to copy. Do not attempt a
whole-file rewrite as a side quest: split out the region under change, test it, move on.

Known offenders, worst first: `docker/app/main.py` (1081), `install_references.py` (901),
`generate_report.py` (861), `cohort_summary.py` (856), `kestrel_genotyping.py` (835),
`pipeline.py` (735), `cli.py` (700), `tests/integration/test_pipeline_integration.py` (667).

## Testing

- `pytest.ini` at the repo root is the live config. The `[tool.pytest.ini_options]`
  block in `pyproject.toml` is **dead** — pytest.ini takes precedence. Edit `pytest.ini`.
- Markers: `unit`, `integration`, `docker`, `smoke`, `slow`. `smoke` is the fast image
  tier (`make test-docker-smoke`, ~1 s): it asserts the built image's *structure* — that
  every reference path `config.json` declares actually exists in it, that the three conda
  envs and their interpreters are present, that adVNTR imports under Python 2.7, and that
  the image stays under its size budget. It needs a Docker daemon but **no test data**,
  and it derives its assertions from the config inside the container rather than from a
  hardcoded list, so it cannot drift. Smoke tests carry **only** the `smoke` marker —
  adding `docker` too would make the slow tier re-run them.
- **Every new unit test file must declare `pytestmark = pytest.mark.unit`.** CI runs
  `pytest -m unit`, so an unmarked file silently never runs. This is enforced by
  `tests/unit/test_marker_hygiene.py`, which fails the build naming the offending file;
  `--strict-markers` additionally turns a *misspelled* marker into an error.
- Unit tests must stay pure: `tmp_path` + `unittest.mock`, no network, no reference data.
- Integration and docker tests pull a ~1.1 GB Zenodo archive and MD5 all 114 files; one
  missing file triggers a full re-download. There are **no** `pytest.skip` guards
  anywhere — missing data calls `pytest.exit()` and kills the whole session, unit tests
  included. Set `VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1` to fail fast instead of downloading.
- Correct marker composition is `-m "docker and not slow"`; repeated `-m` flags override
  rather than combine.
- The Docker tier runs at three depths, chosen by how much signal each buys for its
  runtime: PRs get `test-docker-quick` (4 tests, ~10 s), pushes to `main` get
  `test-docker-fast` (everything except `slow`, ~1.5 min locally), and a nightly
  schedule plus `workflow_dispatch` with `full: true` runs `test-docker` including
  adVNTR. adVNTR is off the merge path on purpose: one test costing 15-25 min on a
  2-core runner - more than the rest of the pipeline combined - for an optional
  module. Both adVNTR tests carry `@pytest.mark.timeout(2700)`, overriding the
  global 600 s; that global timeout is right for everything else and CI hardware is
  several times slower than a workstation.

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
- CI gates on PRs: `make lint`, `make type-check`, `make test-unit` across 3.9–3.12,
  plus a Docker image build and quick Docker tests. Formatting is *not* gated in CI —
  run `make format-check` yourself.

## Traps

1. **Config is loaded at import time.** `kestrel_genotyping.py`, `advntr_genotyping.py`
   and `shark_filtering.py` read their JSON into module globals on import, and
   `run_kestrel()` reads the module-level `kestrel_config` rather than its `config`
   argument. `--config-path` cannot override them; tests must patch the module global.
2. **`--config-path` replaces the whole config, it does not merge.** Missing keys raise
   `KeyError` deep in the pipeline (`config["tools"]["java_path"]`, no `.get`).
3. **Rule strings are `eval()`d against DataFrame column names.** `flagging.py` and
   `cross_match.py` evaluate expressions from `kestrel_config.json`,
   `advntr_config.json` and `config.json`. Renaming a column silently turns flags off —
   a `NameError` only logs a warning. Grep the JSON configs before renaming any column.
4. **Stages mark, they do not filter.** Kestrel stages append boolean columns
   (`is_frameshift`, `is_valid_frameshift`, `depth_confidence_pass`, `alt_filter_pass`,
   `motif_filter_pass`); `filter_final_dataframe()` ANDs them at the end. Preserve that
   contract — dropping rows early breaks `kestrel_pre_result.tsv` debuggability.
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
    `install_references*` or `dependencies/advntr/**` change. `docker/Dockerfile` holds
    only the application and builds on top in ~3 min. Both workflows compute the same
    `hashFiles()` tag, so changing a base input without publishing a new base makes the
    app build **fail loudly** with a missing tag rather than silently using a stale base.
    To change a base input: merge to `main` (or run `docker-base.yml` via
    `workflow_dispatch`) so the base publishes first. Locally:
    `make docker-build-base && make docker-build DOCKER_BASE_IMAGE=vntyper-base:local`.
    (This replaced a `RUN git clone` of GitHub `main`, which meant PR builds never
    tested PR code and a cached layer could serve an arbitrarily stale checkout.)
11. **`vntyper report` is currently broken** — `cli.py` passes arguments
    `generate_summary_report()` does not accept. Do not copy that call site as an example.
12. **Version lives in three places**: `vntyper/version.py` (authoritative),
    `CITATION.cff`, and `docs/about/changelog.md`. `publish-pypi.yml` refuses to publish
    if the pushed tag disagrees with `vntyper/version.py`.
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
14. **`pyproject.toml` and the conda env disagree on numpy.** `pyproject.toml` pins
    `numpy>=1.26.0,<2.0.0`; `conda/environment_vntyper.yml` installs `numpy=2.0.2`. The
    Dockerfile therefore installs the package with `pip install --no-deps`, or pip
    downgrades conda's numpy from PyPI and mixes provenance. Reconcile the two before
    removing `--no-deps`.

## Never

- Never push a `v*.*.*` tag. It publishes to PyPI immediately and irreversibly.
- Never commit into `tests/data/`, `reference/`, `out/`, `download/`, `vntyper.egg-info/`
  or the local `adVNTR/` clone — all are generated or downloaded.
- Never hand-edit `vntyper/dependencies/kestrel/*.jar` or anything in `vntyper.egg-info/`.
- Never add a page under `docs/` without registering it in `mkdocs.yml` `nav:` —
  the docs workflow runs `mkdocs build --strict`, so a dangling nav entry or broken
  internal link fails CI.
- Never claim tests pass without showing the command output.
