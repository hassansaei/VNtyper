# CI/CD follow-ups

Open work left by the pipeline rebuild (PR #176, #180). Two categories, because they
need different handling: **behaviour this work changed**, which users may notice, and
**pre-existing debt** the work surfaced but deliberately did not fix.

Everything here is evidence-backed — measurements and file references, not guesses.

## A. Behaviour this work changed

These are intentional, verified, and listed so nothing arrives as a surprise.

| # | Change | Evidence it is safe | Residual risk |
| --- | --- | --- | --- |
| A1 | Image interpreter 3.9.19 → 3.12.13 | Every package on the numerical path resolves identically; genotype for `example_b178_hg19_subset` reproduces exactly (alt depth 416, active region 7088, `High_Precision*`) | Reporting stack moved (matplotlib, plotly, igv-reports, jinja2) — HTML report rendering is not asserted by any test |
| A2 | Image numpy 1.26.4 → 2.0.2 | The old image only ran 1.26.4 because `pip install .` downgraded conda's numpy; 2.0.2 is what the recipe always declared | numpy 2 semantics differ in edge cases; only the unit tier and one pipeline case are covered |
| A3 | `requires-python` >=3.9 → >=3.10 | 3.9 is EOL (2025-10-31) | Users pinned to 3.9 cannot `pip install` the new version |
| A4 | numpy upper bound `<2.0.0` removed | Needed to stop pip clobbering conda's numpy; 330 tests pass on numpy 2.5.1 | PyPI users may now resolve numpy versions CI does not test |
| A5 | adVNTR moved off the merge path | Nightly cron + `workflow_dispatch full: true` still run it | An adVNTR regression is caught within a day, not at merge |
| A6 | Local conda env must be recreated | `conda/environment_vntyper.yml` changed interpreter | `mamba env create -f conda/environment_vntyper.yml` then `pip install -e ".[dev,docs]"` |

**A1 residual risk is the one worth acting on.** `generate_report.py` is 4% covered and
the reporting stack is exactly what moved. A single test that renders a report from a
fixture `pipeline_summary.json` and asserts the HTML contains the expected variant row
would close it cheaply.

## B. Pre-existing debt

Ordered by expected value, not size.

### B1. Test-data fetching is all-or-nothing — HIGH

`scripts/download_test_data.py` re-downloads the full **1120 MB** Zenodo archive when a
single file is missing, and `tests/support/data_utils.py` calls `pytest.exit()` on failure,
killing the whole session including unit tests that need no data.

Observed: two of 114 files were absent locally (the config had moved ahead of the local
copy), costing a full 1.1 GB refetch to recover 3.5 MB.

The metadata for something better already exists — `tests/test_data_config.json` carries
a per-file `md5sum` for all 114 entries.

Plan:
1. Replace `pytest.exit()` with `pytest.skip(allow_module_level=True)` so a fresh clone
   yields `330 passed, N skipped` instead of a session abort. Keep the hard failure when
   `GITHUB_ACTIONS=true`, where a missing file means a genuine cache fault.
2. Cache the archive rather than streaming it to a `NamedTemporaryFile` that is deleted
   in `finally`, so a single missing file costs one `zip_ref.open(member)`.
3. Record the archive version alongside the data, so a config bump reports "your data is
   from v2.0, config wants v2.1" rather than "2 files missing".

### B2. Six modules over 650 LOC, all under 30% covered — HIGH

Every module above the size limit in `AGENTS.md` is barely tested, and none of the
well-covered ones exceed it:

| Module | LOC | Coverage |
| --- | --- | --- |
| `cohort_summary.py` | 856 | 0% |
| `cli.py` | 700 | 0% |
| `generate_report.py` | 861 | 4% |
| `pipeline.py` | 735 | 10% |
| `install_references.py` | 901 | 24% |
| `kestrel_genotyping.py` | 835 | 28% |
| `docker/app/main.py` | 1081 | not measured |

They are untested *because* they fuse I/O, orchestration and pure logic into functions
that cannot be called without a filesystem. The route from 25.68% to the 80% target runs
through splitting them, one region at a time, under the rule already in `AGENTS.md`:
extract the pure part, test that, leave the I/O behind.

### B3. `tests/` is not linted or formatted — MEDIUM

`ruff check tests/` reports **30 fixable issues**; `ruff format --check tests/` wants to
reformat **12 of 35 files**. CI lints `vntyper/` only, so this is invisible.

Deliberately out of scope for the pipeline work — expanding lint scope inside a CI PR
would have mixed a cleanup with an infrastructure change. It is a self-contained PR:
run the formatter, fix or explicitly ignore the 30, then add `tests/` to `make lint`.

### B4. PyPI Trusted Publishing — MEDIUM

`publish-pypi.yml` still authenticates with a long-lived `PYPI_API_TOKEN`. The migration
is documented in the workflow header and is a two-part change that **must** be done in
order, or the next release fails:

1. On pypi.org: Manage project → Publishing → add a GitHub publisher
   (owner `hassansaei`, repo `VNtyper`, workflow `publish-pypi.yml`, environment `pypi`).
2. In the workflow: add `id-token: write` and the `environment:` block, replace the twine
   step with `pypa/gh-action-pypi-publish@release/v1`, delete the secret.

### B5. `vntyper report` is broken — MEDIUM

`cli.py` passes arguments `generate_summary_report()` does not accept, so the subcommand
raises. Recorded as trap 11 in `AGENTS.md` and untouched here. Needs its own fix plus the
regression test that would have caught it.

### B6. Workflow linting is local-only — MEDIUM

`actionlint` runs in `make ci-local` and in no GitHub Actions workflow, so nothing in CI
checks workflow syntax or shellchecks the `run:` blocks. A finding can sit on `main`
indefinitely, and one did: SC2129 in `docker-base.yml` predates this work.

That gap was widened by a bug in the gate itself. `lint-actions` expanded an empty
`$(ACTIONLINT)` straight into command position, leaving a bare `;` — a shell *syntax*
error, raised while parsing the whole `if ... fi` compound, so it fired before the
`[ -n ... ]` guard could select the container fallback. `make ci-local` therefore died at
its first step on every machine without a local `actionlint`, which is precisely the
population the fallback was written to serve. Fixed by quoting the expansion and pinned by
`tests/unit/test_makefile_recipes.py`, which renders each recipe with the tool variable
forced empty and syntax-checks it.

Plan: add an `actionlint` step to `ci-tests.yml`, gated on `.github/workflows/**` changing,
so the check is not contingent on a contributor happening to run the local target.

### B7. Unexplained: base rebuild reported as skipped — LOW, but verify

On the run that bumped `python-multipart`, `Check base image` resolved a new content hash
(`base-c5dec6b6f780ab4f`) and `Build base image` displayed **skipped**, yet that tag
exists in ghcr with a digest distinct from the previous base — so a rebuild demonstrably
happened and the app built on the correct base.

The outcome is right; the job-status display is not explained. Most likely a
reusable-workflow reporting artifact. Worth confirming once more on the next base-input
change before trusting the `skipped` signal, since "base silently not rebuilt" would be a
serious failure mode.

## C. Deliberately not doing

- **Sharding the unit suite.** It runs in 0.44 s. Measured: `-n auto` is **4.8× slower**
  (2.08 s vs 0.43 s) because worker bootstrap re-imports pandas/numpy/pysam per worker.
- **`pytest-testmon`.** Its value is skipping expensive tests; every test here is <10 ms,
  and its stale-database failure mode silently skips tests — the exact bug class that hid
  30 tests for months.
- **arm64 images.** Roughly doubles base build cost; bioconda's `linux-aarch64` coverage
  would need a `CONDA_SUBDIR=linux-aarch64 --dry-run` probe first, and demand is unconfirmed.
