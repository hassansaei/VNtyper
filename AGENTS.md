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
| Format + autofix | `make format` |
| Lint | `make lint` |
| Type check | `make type-check` |
| Fast tests (what CI runs) | `make test-unit` |
| Integration tests | `make test-integration` |
| Coverage | `make test-cov` |
| Docs preview | `make docs-serve` |
| Docs build (CI-equivalent) | `make docs-build` |

Always run `make check-all` before opening a PR. Run pytest **from the repo root** —
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
- Code must run on Python 3.9 (CI matrix: 3.9–3.12). Your local interpreter is newer,
  so 3.10+ syntax can pass locally and break CI.
- Google-style docstrings (`Args:` / `Returns:` / `Raises:`) on public functions.
- Logging: modules call `logging.info(...)` on the root logger; configuration happens
  once in `cli.py`. Do not add `basicConfig` or per-module handlers.
- Errors: no custom exception classes. The convention is `logging.error(msg)` followed
  by `raise ValueError(msg)` / `RuntimeError(msg)` with the same message. Only exit
  codes 0 and 1 are used.
- Type hints are partial. Newer modules (`reference_registry`, `region_utils`,
  `scoring`, `flagging`) are fully annotated — match that when editing them.

## Testing

- `pytest.ini` at the repo root is the live config. The `[tool.pytest.ini_options]`
  block in `pyproject.toml` is **dead** — pytest.ini takes precedence. Edit `pytest.ini`.
- Markers: `unit`, `integration`, `docker`, `slow`.
- **Every new unit test file must declare `pytestmark = pytest.mark.unit`.** CI runs
  `pytest -m unit`, so an unmarked file silently never runs.
  (`tests/unit/test_haplo_count_and_selection.py` is currently unmarked — 30 tests
  invisible to CI.)
- Unit tests must stay pure: `tmp_path` + `unittest.mock`, no network, no reference data.
- Integration and docker tests pull a ~1.1 GB Zenodo archive and MD5 all 114 files; one
  missing file triggers a full re-download. There are **no** `pytest.skip` guards
  anywhere — missing data calls `pytest.exit()` and kills the whole session, unit tests
  included. Set `VNTYPER_TEST_DATA_SKIP_DOWNLOAD=1` to fail fast instead of downloading.
- Correct marker composition is `-m "docker and not slow"`; repeated `-m` flags override
  rather than combine.

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
   `run_kestrel()` uses `global kestrel_config` rather than its `config` argument.
   `--config-path` cannot override them; tests must patch the module global.
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
10. **The Docker build clones from GitHub**, it does not build your local checkout. Only
    `docker/entrypoint.sh` and `docker/app/` come from the local context.
11. **`vntyper report` is currently broken** — `cli.py` passes arguments
    `generate_summary_report()` does not accept. Do not copy that call site as an example.
12. **Version lives in three places**: `vntyper/version.py` (authoritative),
    `CITATION.cff`, and `docs/about/changelog.md` (currently stale at 2.0.1).

## Never

- Never push a `v*.*.*` tag. It publishes to PyPI immediately and irreversibly.
- Never commit into `tests/data/`, `reference/`, `out/`, `download/`, `vntyper.egg-info/`
  or the local `adVNTR/` clone — all are generated or downloaded.
- Never hand-edit `vntyper/dependencies/kestrel/*.jar` or anything in `vntyper.egg-info/`.
- Never add a page under `docs/` without registering it in `mkdocs.yml` `nav:` —
  the docs workflow runs `mkdocs build --strict`, so a dangling nav entry or broken
  internal link fails CI.
- Never claim tests pass without showing the command output.
