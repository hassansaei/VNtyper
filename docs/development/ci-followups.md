# CI/CD follow-ups

Open work left by the pipeline rebuild (PR #176, #180). Two categories require distinct handling: **behaviour this work changed**, which users may notice, and **pre-existing debt** the work surfaced but deliberately left open.

All entries are evidence-backed with measurements and file references rather than speculation.

!!! warning "Every measurement on this page is dated, and several are now historical"

    This page records state at PR #176/#180. Line counts and coverage figures shift with each merge. None of the numbers below represents a live claim: each carries the date it was taken, and resolved entries name their resolving pull requests. Several original figures (`generate_report.py` at 861 lines / 4%, `cohort_summary.py` at 856, `cli.py` at 700) were hundreds of lines out of date before this note was added.

    For the current picture, read **"Changing existing code" in `AGENTS.md`**: that section contains the maintained inventory and requires fresh measurements before quotation.

## A. Behaviour this work changed

These changes are intentional, verified, and catalogued to prevent downstream surprises.

| # | Change | Evidence of safety | Residual risk |
| --- | --- | --- | --- |
| A1 | Image interpreter 3.9.19 -> 3.12.13 | Numerical path packages resolve identically; genotype for `example_b178_hg19_subset` reproduces exactly (alt depth 416, active region 7088, `High_Precision*`) | Reporting stack moved (matplotlib, plotly, igv-reports, jinja2); HTML report rendering was initially unasserted |
| A2 | Image numpy 1.26.4 -> 2.0.2 | Legacy image ran 1.26.4 because `pip install .` downgraded conda's numpy; 2.0.2 is what the recipe always declared | numpy 2 semantics differ in edge cases; unit tier and pipeline smoke tests verify numerical stability |
| A3 | `requires-python` >=3.9 -> >=3.10 | Python 3.9 reached end of life (2025-10-31) | Users pinned to 3.9 cannot `pip install` newer releases |
| A4 | numpy upper bound `<2.0.0` removed | Prevents pip from overwriting conda's numpy; suite passes on numpy 2.x | PyPI installations resolve modern numpy distributions tested in CI |
| A5 | adVNTR moved off merge path | Nightly cron and manual workflow dispatches run the full Docker tier | adVNTR regressions are caught within 24 hours rather than blocking every pull request |
| A6 | Local conda environment rebuild | `conda/environment_vntyper.yml` updated interpreter | Run `mamba env create -f conda/environment_vntyper.yml` followed by `pip install -e ".[dev,docs]"` |

**A1 residual risk: RESOLVED (#179).** When written, `generate_report.py` was 861 lines and 4% covered while the interpreter bump shifted the reporting stack. The fix was implemented in `tests/unit/test_generate_report.py` across 43 test functions: each writes a `pipeline_summary.json` fixture into `tmp_path`, invokes `generate_summary_report()`, and verifies rendered HTML elements (coverage rows, Kestrel table headings, annotated `Motif` column, confidence colour coding, negative placeholder row, and screening summary styling). `generate_report.py` is now 574 lines (`wc -l vntyper/scripts/generate_report.py`); see `AGENTS.md` for current coverage figures.

## B. Pre-existing debt

Ordered by expected value rather than size.

### B1. Test-data fetching is all-or-nothing (HIGH)

`scripts/download_test_data.py` redownloads the full 1120 MB Zenodo archive when a single file is missing. In addition, `tests/support/data_utils.py` calls `pytest.exit()` on failure, terminating the entire test session including unit tests that require no external data.

Observed fault: two of 114 files were absent locally because config advanced ahead of the local copy, triggering a 1.1 GB refetch to recover 3.5 MB.

The metadata for fine-grained retrieval exists: `tests/test_data_config.json` stores per-file `md5sum` entries for all 114 items.

Plan:
1. Replace `pytest.exit()` with `pytest.skip(allow_module_level=True)` so fresh clones report unit tests passing and data tests skipped rather than aborting. Maintain hard failures when `GITHUB_ACTIONS=true` where missing files indicate genuine cache corruption.
2. Cache the downloaded archive rather than streaming into a temporary file discarded in a `finally` block, ensuring missing files require only extraction.
3. Store archive versions alongside data, enabling clear mismatch warnings when configuration updates.

### B2. Six modules over 650 LOC, all under 30% covered: SUPERSEDED (conclusion was incorrect)

**Historical snapshot, PR #176/#180.** Preserved for historical audit because the conclusion drawn from them proved false.

| Module | LOC (then) | Coverage (then) |
| --- | --- | --- |
| `cohort_summary.py` | 856 | 0% |
| `cli.py` | 700 | 0% |
| `generate_report.py` | 861 | 4% |
| `pipeline.py` | 735 | 10% |
| `install_references.py` | 901 | 24% |
| `kestrel_genotyping.py` | 835 | 28% |
| `docker/app/main.py` | 1081 | not measured |

The original claim that every module over 650 LOC has under 30% coverage is false, as is its converse. `docker/app/main.py` is the largest file in the repository and among the best covered because it comprises modular FastAPI handlers testable via `TestClient`. `install_references.py` is similar in size but exhibits low coverage because it downloads, unpacks, and checksums filesystem assets. Coupling to I/O predicts coverage; line count does not. `AGENTS.md` ("Changing existing code") maintains current metrics and design guidelines.

Three of the seven rows have since been resolved through refactoring:

| Module | LOC (then) | LOC (now) | Modifications |
| --- | ---: | ---: | --- |
| `cohort_summary.py` | 856 | 456 | Split into `cohort_rules`, `cohort_categories`, `cohort_tables`, `cohort_inputs`, `cohort_exports` |
| `cli.py` | 700 | 173 | Parser moved to `cli_parser.py`, handlers to `cli_handlers.py`, report logic to `cli_report.py` (#179) |
| `generate_report.py` | 861 | 574 | Presentation logic moved to `screening_summary.py` and `report_formatting.py` |
| `pipeline.py` | 735 | 721 | Retains legacy orchestration scope |
| `install_references.py` | 901 | 901 | High I/O coupling persists |
| `kestrel_genotyping.py` | 835 | 877 | Complex variant calling logic |
| `docker/app/main.py` | 1081 | 1151 | Modular route handlers with high unit test coverage |

Four production files remain above 650 LOC: `docker/app/main.py`, `install_references.py`, `kestrel_genotyping.py`, and `pipeline.py`. Follow-up work applies the separation rule documented in `AGENTS.md`: isolate pure logic into dedicated modules, verify thoroughly, and leave I/O wrappers thin.

### B3. Test suite linting and formatting: RESOLVED

`RUFF_PATHS` in `Makefile` covers `vntyper/ docker/app/ tests/ scripts/ docs/`. All Ruff targets (`lint`, `lint-stats`, `format`, `format-check`) consume this variable, ensuring `tests/` is linted and formatted under identical rules to production code. The path set matches what bare `ruff format --check .` discovers.

### B4. PyPI Trusted Publishing: RESOLVED, cleanup complete

The default-branch `publish-pypi.yml` controller implements PyPI Trusted Publishing. Its production publisher executes in the protected `pypi` environment with minimal `id-token: write` permissions, consumes the exact candidate artifact, and invokes a commit-pinned PyPA publisher with `skip-existing`. It does not read `PYPI_API_TOKEN` or any long-lived secret.

The `pypi` environment requires no reviewers, enforces no wait timer, and applies custom branch policies restricted to `main` without custom deployment protection rules. The release controller preflights environment definitions, branch policies, and deployment rules via GitHub APIs, aborting before package or registry mutations upon any mismatch. OIDC remains the sole authorized publisher; long-lived package tokens are strictly prohibited.

Cleanup complete: production OIDC runs `31465885545` (2.0.11) and `31464328451` (2.0.12) succeeded, PyPI 2.0.12 is published, and `PYPI_API_TOKEN has been deleted`. Obsolete credentials cannot be retrieved by historical workflows.

### B5. Subcommand `vntyper report` invocation: RESOLVED (#179)

`cli.py` previously supplied arguments rejected by `generate_summary_report()`, causing `TypeError` before execution. Issue #179 corrected this: handler dispatch moved to `vntyper/scripts/cli_report.py`, rejected keyword arguments were removed (the generator extracts them directly from `pipeline_summary.json`), and the required `log_file` parameter is passed explicitly. Regression testing executes `cli.main(["report", ...])` against an `inspect.signature` wrapper to ensure interface stability.

### B6. Action workflow validation is local-only (MEDIUM)

`actionlint` runs in `make ci-local` but is not enforced on GitHub Actions workflows. Workflow syntax errors or shell script regressions in `run:` blocks can evade pre-merge checks.

A past bug in `lint-actions` expanded an unquoted `$(ACTIONLINT)` variable into command position, producing a shell syntax error on systems without local `actionlint` installations. Quoting the expansion resolved the issue, backed by `tests/unit/test_makefile_recipes.py`.

Plan: integrate `actionlint` into `ci-tests.yml` triggered on changes within `.github/workflows/**`.

### B7. Base container image rebuild display anomaly (LOW, verification needed)

During dependency updates for `python-multipart`, the `Check base image` step resolved a new content digest (`base-c5dec6b6f780ab4f`) while `Build base image` reported `skipped`. Inspection of GHCR confirmed the new digest was pushed and the application built upon the updated base.

The runtime outcome was correct, but the skipped job status indicator warrants verification on future base recipe changes to prevent unobserved build skips.

### B8. adVNTR `Insertion_len` extraction rule requires domain validation (HIGH)

`advntr_processing_del` and `advntr_processing_ins` in `vntyper/modules/advntr/advntr_genotyping.py` extract `Insertion_len` from adVNTR `State` strings using:

```python
df1["Insertion_len"] = df1["Variant"].str.extract(r"(LEN.*)")[0]      # greedy match
df1["Insertion_len"] = df1["Insertion_len"].fillna("LEN")
df1[["I", "Insertion_len"]] = df1["Insertion_len"].str.split("LEN", n=1, expand=True)
df1["Insertion_len"] = pd.to_numeric(df1["Insertion_len"], errors="coerce").fillna(0).astype(int)
```

`Insertion_len` feeds `frame = abs(Insertion_len - Deletion_length)`, retaining rows only when `frame` falls into configured frameshift intervals. In compound calls with multiple segments separated by `&`, extracting trailing characters often yields non-numeric strings:

| `State` | Extracted string | Resulting `Insertion_len` | Mechanism |
| --- | --- | --- | --- |
| `I22_2_G_LEN1` | `1` | 1 | Single segment containing integer |
| `D8_2&D9_2&I9_2_A_LEN9` | `9` | 9 | `LEN` occurs in terminal segment |
| `I9_2_A_LEN2&D50_2` | `2&D50_2` | **0** | Trailing segments follow `LEN` |
| `I9_2_A_LEN2&D50_2&D51_2` | `2&D50_2&D51_2` | **0** | Trailing segments follow `LEN` |
| `I9_2_A_LEN9&I50_2_A_LEN3` | `9&I50_2_A_LEN3` | **0** | Multiple `LEN` tokens present |

Under the current implementation, any `&` segment following `LEN` forces the insertion length to zero. Candidate policies behave differently on compound calls:

| Policy | `I9_2_A_LEN2&D50_2` | `I9_2_A_LEN9&I50_2_A_LEN3` |
| --- | --- | --- |
| **NaN->0** (current) | 0 (insertion filter drops row) | 0 (dropped) |
| **First LEN wins** | 2 (frame 1, row passes filter) | 9 (frame 9, dropped) |
| **Sum of LENs** | 2 (passes) | 12 (frame 12, dropped) |

A sample whose sole adVNTR call is `I9_2_A_LEN2&D50_2` yields a negative placeholder under current rules, but produces a positive adVNTR call under the "first LEN wins" rule.

Current behavior is locked byte-for-byte by `tests/unit/test_advntr_output_parsing.py::TestInsertionLenIsCharacterised`. Modifying this behavior requires a formal clinical domain decision regarding compound MUC1 VNTR calls.

Related nuances:
* `Deletion_length` evaluates `State.str.count("D")`, counting the character `D` across the string rather than discrete deletion segments.
* `Insertion_length` (`str.count("I")`) is calculated but unconsumed by filtering logic.

## C. Deliberately excluded optimizations

- **Unit test runner sharding**: When evaluated at PR #176/#180, `-n auto` was 4.8x slower (2.08 s vs 0.43 s) due to interpreter boot overhead and re-importing pandas, numpy, and pysam across worker processes.
- **`pytest-testmon`**: Optimizes execution by skipping tests based on file modifications. Because unit tests complete in under 10 ms, database cache invalidation risks silently skipping tests without offering performance benefits.
- **arm64 container images**: Doubles base image build overhead. Bioconda `linux-aarch64` package parity requires testing before enabling multi-architecture builds.
