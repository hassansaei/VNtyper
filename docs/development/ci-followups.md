# CI/CD follow-ups

Open work left by the pipeline rebuild (PR #176, #180). Two categories, because they
need different handling: **behaviour this work changed**, which users may notice, and
**pre-existing debt** the work surfaced but deliberately did not fix.

Everything here is evidence-backed — measurements and file references, not guesses.

!!! warning "Every measurement on this page is dated, and several are now historical"

    This page records the state at PR #176/#180. Line counts and coverage figures move
    with every merge, so **none of the numbers below are a live claim** — each is
    labelled with when it was taken, and an entry that has since been resolved says so
    and names what resolved it. Several of the original figures (`generate_report.py` at
    861 lines / 4%, `cohort_summary.py` at 856, `cli.py` at 700) were hundreds of lines
    out of date before this warning was added.

    For the current picture, read **"Changing existing code" in `AGENTS.md`** — that is
    the maintained table, and it tells you to re-measure before quoting it too.

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

**A1 residual risk — RESOLVED (#179).** At the time this was written `generate_report.py`
was 861 lines and 4% covered, and the reporting stack was exactly what the interpreter
bump moved. The proposed fix — "a single test that renders a report from a fixture
`pipeline_summary.json` and asserts the HTML contains the expected variant row" — is what
`tests/unit/test_generate_report.py` now does, across 43 test functions
(`grep -c '^def test' tests/unit/test_generate_report.py`): each writes a
`pipeline_summary.json` into `tmp_path`, calls the real `generate_summary_report()`, and
asserts against the rendered HTML (coverage rows, the Kestrel table headings, the
annotated `Motif` column, the confidence colour coding, the negative placeholder row and
the screening-summary styling). `generate_report.py` is 574 lines today
(`wc -l vntyper/scripts/generate_report.py`); for its coverage see the table in
`AGENTS.md`.

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
   yields `<all unit tests> passed, N skipped` instead of a session abort (the suite was
   330 tests when this was written and is several times that now). Keep the hard failure
   when `GITHUB_ACTIONS=true`, where a missing file means a genuine cache fault.
2. Cache the archive rather than streaming it to a `NamedTemporaryFile` that is deleted
   in `finally`, so a single missing file costs one `zip_ref.open(member)`.
3. Record the archive version alongside the data, so a config bump reports "your data is
   from v2.0, config wants v2.1" rather than "2 files missing".

### B2. Six modules over 650 LOC, all under 30% covered — SUPERSEDED; the conclusion was wrong

**Historical snapshot, PR #176/#180.** These are the figures as recorded then. They are
kept rather than deleted because the *conclusion* drawn from them turned out to be false,
and a page that quietly rewrites its own history stops being worth trusting.

| Module | LOC (then) | Coverage (then) |
| --- | --- | --- |
| `cohort_summary.py` | 856 | 0% |
| `cli.py` | 700 | 0% |
| `generate_report.py` | 861 | 4% |
| `pipeline.py` | 735 | 10% |
| `install_references.py` | 901 | 24% |
| `kestrel_genotyping.py` | 835 | 28% |
| `docker/app/main.py` | 1081 | not measured |

The heading's claim — *every module over 650 LOC is under 30% covered* — **is false, and
so is its converse.** `docker/app/main.py` is the largest file in the repository and one
of the best covered, because it is a stack of thin FastAPI handlers reachable through
`TestClient`. `install_references.py` is the same size and sits low because it downloads,
unpacks and checksums files. **Coupling to I/O predicts coverage; size does not.**
`AGENTS.md`, "Changing existing code", carries the current table and the full argument.

Three of the seven rows have since been closed, and the line counts are now very
different. Measured at the tip of the `#181–#197` follow-up branch with
`wc -l`:

| Module | LOC (then) | LOC (now) | What changed |
| --- | ---: | ---: | --- |
| `cohort_summary.py` | 856 | 456 | split into `cohort_rules`, `cohort_categories`, `cohort_tables`, `cohort_inputs`, `cohort_exports` |
| `cli.py` | 700 | 173 | parser → `cli_parser.py`, handlers → `cli_handlers.py`, `report` → `cli_report.py` (#179) |
| `generate_report.py` | 861 | 574 | presentation logic → `screening_summary.py`, `report_formatting.py`; see A1 |
| `pipeline.py` | 735 | 721 | still over the limit |
| `install_references.py` | 901 | 901 | still over the limit, and still the one module the old rule describes |
| `kestrel_genotyping.py` | 835 | 877 | still over the limit, and larger than it was |
| `docker/app/main.py` | 1081 | 1151 | still over the limit, and well covered |

So **four** production files remain over 650 LOC, not six: `docker/app/main.py`,
`install_references.py`, `kestrel_genotyping.py` and `pipeline.py`. The still-open work
is the splitting rule already in `AGENTS.md` — extract the pure part, test that, leave
the I/O behind — applied to those four, worst first.

### B3. `tests/` is not linted or formatted — RESOLVED

`RUFF_PATHS` in the `Makefile` is now
`vntyper/ docker/app/ tests/ scripts/ docs/`, and every ruff target (`lint`,
`lint-stats`, `format`, `format-check`) reads from it, so `tests/` is linted **and**
formatted by the same gates that cover `vntyper/`. The set is deliberately equal to what
a bare `ruff format --check .` discovers, so the Makefile target and the obvious command
cannot report different things; the comment above `RUFF_PATHS` gives the command that
checks that.

### B4. PyPI Trusted Publishing — RESOLVED; owner cleanup remains

The current default-branch `publish-pypi.yml` controller uses PyPI Trusted Publishing:
its production publisher runs in the protected `pypi` environment with only
`id-token: write`, consumes the exact candidate artifact, and invokes a commit-pinned
PyPA publisher with `skip-existing`. It does not read `PYPI_API_TOKEN` or another
long-lived package credential.

One owner-only rollout action remains intentionally pending. Do not delete
`PYPI_API_TOKEN` until the first successful OIDC release proves that the configured
publisher (owner `hassansaei`, repository `VNtyper`, workflow `publish-pypi.yml`,
environment `pypi`) works live. After that proof, the owner must delete the obsolete
secret separately and record the release run. Until deletion, do not create a release
tag pointing at a pre-milestone commit: historical tagged commits still contain their
legacy token workflow. Those old workflows become inert only after the token is
removed, so the current workflow migration is resolved while secret retirement is not.

### B5. `vntyper report` is broken — RESOLVED (#179)

`cli.py` passed arguments `generate_summary_report()` does not accept, so the subcommand
raised `TypeError` before doing any work. Fixed in #179: the handler moved to
`vntyper/scripts/cli_report.py`, the three unaccepted keywords were dropped (the report
generator reads all three out of `pipeline_summary.json` itself) and the required
`log_file` is now passed. The regression test drives the real
`cli.main(["report", ...])` into a spy that binds
`inspect.signature(generate_summary_report)`, so a call the real function would reject
fails in the same place.

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

### B8. adVNTR `Insertion_len` is derived by a rule nobody has ratified — HIGH, needs a domain decision

**This one is not a code question. It decides which adVNTR genotypes get reported, so it
is reserved for the domain owner and is deliberately left unchanged in #179.**

`advntr_processing_del` and `advntr_processing_ins` in
`vntyper/modules/advntr/advntr_genotyping.py` derive `Insertion_len` from the adVNTR
`State` string like this:

```python
df1["Insertion_len"] = df1["Variant"].str.extract(r"(LEN.*)")[0]      # greedy, to end of string
df1["Insertion_len"] = df1["Insertion_len"].fillna("LEN")
df1[["I", "Insertion_len"]] = df1["Insertion_len"].str.split("LEN", n=1, expand=True)
df1["Insertion_len"] = pd.to_numeric(df1["Insertion_len"], errors="coerce").fillna(0).astype(int)
```

`Insertion_len` then feeds `frame = abs(Insertion_len - Deletion_length)`, and the row
survives only if `frame` is in the configured frameshift series. A compound call names
several parts joined by `&`, so the extracted remainder is often not a number at all:

| `State` | extracted remainder | `Insertion_len` | why |
| --- | --- | --- | --- |
| `I22_2_G_LEN1` | `1` | 1 | single part — a bare number |
| `D8_2&D9_2&I9_2_A_LEN9` | `9` | 9 | the `LEN` is in the *last* part |
| `I9_2_A_LEN2&D50_2` | `2&D50_2` | **0** | a further part follows the `LEN` |
| `I9_2_A_LEN2&D50_2&D51_2` | `2&D50_2&D51_2` | **0** | as above |
| `I9_2_A_LEN9&I50_2_A_LEN3` | `9&I50_2_A_LEN3` | **0** | two `LEN` tokens |

So the effective rule is: *if any `&` part follows the `LEN` token, the inserted length is
treated as zero.* That is almost certainly not what anyone intended. There are three
candidate policies and they do not agree:

| Policy | `I9_2_A_LEN2&D50_2` | `I9_2_A_LEN9&I50_2_A_LEN3` |
| --- | --- | --- |
| **NaN→0** (today) | 0 → ins filter drops the row | 0 → dropped |
| **first LEN wins** | 2 → `frame` 1 → **row survives the insertion filter** | 9 → `frame` 9 → dropped |
| **sum of LENs** | 2 → survives | 12 → `frame` 12 → dropped |

The middle column is the reason this is filed rather than fixed. Through the
empty-result branch, a sample whose only adVNTR call is `I9_2_A_LEN2&D50_2` reports the
`Negative` placeholder under today's rule and a **positive adVNTR call** under
"first LEN wins"; a sample whose only call is `I9_2_A_LEN2&D50_2&D51_2` moves the other
way. That is a reported-genotype change either way round.

`#179` briefly shipped "first LEN wins" as a side effect of a crash fix
(`a7c3d9e`, reverted to crash-only in the follow-up commit) — which is how the question
was found. Today's behaviour is pinned byte for byte, against measurements taken from the
pre-fix commit, by
`tests/unit/test_advntr_output_parsing.py::TestInsertionLenIsCharacterised`; changing the
policy means changing that test deliberately, which is the point.

**What is needed:** a decision from the domain owner on which of the three is correct for
MUC1 VNTR compound calls, ideally checked against a cohort that actually contains one.
No sample in `tests/data/` does — see the golden-cohort gate's "what this gate does not
cover".

Two smaller questions ride along with it, both currently invisible for the same reason:

* `Deletion_length` is `State.str.count("D")` — a count of the letter `D` anywhere in the
  string, not of deletion *parts*. It happens to agree for the shapes seen so far.
* `Insertion_length` (plural, `str.count("I")`) is computed and never used by either
  filter.

## C. Deliberately not doing

- **Sharding the unit suite.** Measured at PR #176/#180, when the tier ran in 0.44 s:
  `-n auto` was **4.8× slower** (2.08 s vs 0.43 s) because worker bootstrap re-imports
  pandas/numpy/pysam per worker. The tier is much larger now, so the *ratio* is the part
  that carries; re-measure before acting on it either way.
- **`pytest-testmon`.** Its value is skipping expensive tests; every test here is <10 ms,
  and its stale-database failure mode silently skips tests — the exact bug class that hid
  30 tests for months.
- **arm64 images.** Roughly doubles base build cost; bioconda's `linux-aarch64` coverage
  would need a `CONDA_SUBDIR=linux-aarch64 --dry-run` probe first, and demand is unconfirmed.
