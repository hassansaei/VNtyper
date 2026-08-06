# Changelog

All notable changes to VNtyper 2 are documented on this page.

## 2.0.6 (Current)

**Unit coverage went from 25.68% to 70.57%, across 1,744 tests.** The number was not the
point: a mutation experiment had shown it overstated protection badly —
`confidence_assignment.py` had 100% line coverage and a 21% mutation score, meaning its
central condition could be inverted with the whole suite still green. That module now
scores 57.9% raw / 84.6% adjusted at the same 100% line coverage.

**No reported genotype or confidence label changes.** That is measured, not asserted: the
golden-cohort gate ran 122 pipeline runs before and after, and the Kestrel variant set,
`Confidence`, `Flag`, adVNTR genotype and `Flag` are identical on every sample and every
assembly. `vntyper/scripts/kestrel_config.json` is byte-identical to its pre-branch state.
See `docs/development/golden-cohort-gate.md`.

### ⚠️ Required before deploying the image

- **Rotate the two Redis passwords.** They are in this repository's public history;
  removing them from source does not invalidate them.
- **Set `REDIS_PASSWORD`.** The web service now refuses to start without it.
- **Deploy the worker before the API.** `run_vntyper_job` gained an `index_path`
  parameter defaulting to `None`; old-API-to-new-worker is safe, new-API-to-old-worker is
  not.

### Breaking changes (web service only; the CLI is unaffected)

- Alias-only cohort access no longer works — `cohort_id` **and** passphrase are both
  required. Cohorts with no stored passphrase hash are unopenable; treating a missing
  hash as "no passphrase" would be the bypass, and no cohort ever had a working hash
  (see below).
- `/job-status/` returns a generic message for a failed job; the detail is logged
  server-side rather than returned.
- `--output-name` is now refused rather than silently ignored.

### Web service hardening

`docker/app` — 1,679 LOC, internet-facing, processing patient-derived genomic data — had
no tests, no lint, no type checking and no coverage. No automated check had ever read it.
It is now formatted, linted, type-checked and 87% covered, and its path is in CI's filter
so those gates fire.

- **The passphrase check never worked.** passlib 1.7.4's bcrypt backend probe is
  incompatible with bcrypt ≥ 4.1 and raised `ValueError` for every input, so no cohort
  could ever have been protected. Now on `bcrypt` directly.
- **Compose set `CELERY_BROKER_URL`**, which Celery prefers over the constructed URL, so
  the worker's real broker was the unauthenticated one. Removed, with a test blocking its
  return. Two *different* hardcoded Redis defaults were also committed, so with the
  variable unset the API and the worker authenticated differently against the same
  instance.
- Upload filenames are constrained by an allowlist before they build a path; uploads and
  request bodies are bounded; job ids are validated as UUIDs; cohort analysis no longer
  deletes its members' archives; `/cohort-download/` cleans up its temporary archive.

### Pipeline fixes

Each with a test written first and observed failing.

- **`vntyper report` was dead** — an uncaught `TypeError`.
- **The `Motif` column was missing from every per-sample report**; the pipeline emits
  `Motif` and the report keyed on `Motifs`, and unmatched columns are silently dropped.
- **Every negative report was styled as positive**, and screening states with no
  configured message fell through to "the screening was negative" — including
  Kestrel-positive-and-flagged.
- **`--extra-modules advntr,shark` silently disabled both modules** (exact-string
  membership, no comma split), so a typo yielded a silent Kestrel-only run.
- **`_construct_ncbi_accession` compared against `"hg19"`** while its only caller passes
  `"GRCh37"`/`"GRCh38"`, so every non-chr1 GRCh37 NCBI name was built from the GRCh38
  table.
- **`online_mode` exited 0 on a failed job**, so a wrapping `subprocess.run(check=True)`
  saw success.
- **`set -o pipefail` added to three multi-stage pipes** — without it a failing
  `samtools sort` upstream of `samtools fastq` still exits 0, genotyping a truncated
  FASTQ as if complete.
- **The `RU == 7` adVNTR flag now fires.** It was added three hours and 45 minutes
  *before* the `RU` column existed, so it had never once evaluated. RU-7 calls move from
  `positive` to `positive flagged`; no row loses a flag and no genotype field changes.
- Jinja2 autoescaping on both the per-sample and cohort reports; `shlex.quote` on
  interpolated paths, regions and sample names.

### Test infrastructure

- **A patch-coverage gate at 80% of changed lines**, which is what makes "touch a file,
  add tests for it" enforceable — demonstrated failing on an untested one-line change.
- **A mutation-testing harness** with the results and the equivalent-mutant
  classifications committed (`docs/development/mutation-testing.md`).
- **A golden-cohort gate** (`docs/development/golden-cohort-gate.md`).
- `tests/builders.py` domain-object builders, `summary_steps.py` for the five step-name
  literals four modules matched by exact string, and `calibration.json` recording the
  provenance of the ten calibrated constants.

### Characterised, not fixed

Seventeen findings needing a domain decision or exceeding a hardening PR's scope were
each pinned by a test and filed as #181–#197. Three worth reading: **#186**
(`motifs_for_alt_gg` ships empty, inert in the active branch but destructive in the
uniform one), **#181** (every (3n+1)-bp deletion discarded), and **#185**
(`filter_final_dataframe` fails open).

Closes #179.

## 2.0.5

**Docker image now runs Python 3.12.** Python 3.9 reached end of life on 2025-10-31.
Every package on the numerical path is unchanged - bwa, samtools, fastp, bcftools,
openjdk (Kestrel), pysam, numpy, pandas and biopython all resolve to the same versions -
and the pipeline reproduces the expected genotype exactly. Only the reporting stack
moved (matplotlib, plotly, igv-reports, jinja2). The package now requires Python 3.10+
and is tested on 3.10-3.13.

**CI/CD pipeline rebuilt.** The Docker image is split into a rarely-rebuilt base
(conda environments, adVNTR, reference genomes) and a per-commit application image, so a
commit now rebuilds in about 3 minutes instead of 40-70. Images are published only after
their tests pass, and the build uses the repository's own source rather than cloning
GitHub, so a pull request now tests the code in that pull request.

- Fixed `pyproject.toml` declaring `numpy<2.0.0` while the conda environment installed
  numpy 2.0.2: `pip install .` inside the image downgraded conda's numpy, so the
  published image ran a different numerical stack than its own recipe declared.
- Fixed the TSV/CSV summary parsers silently truncating rows whose field count did not
  match the header; a malformed row is now logged and skipped without discarding the
  rest of the file.
- Recovered 30 unit tests that carried no `unit` marker and had therefore never run in
  CI, covering the Issue #136/#145 genotype tie-breaking logic.
- Added a fast image smoke tier that verifies every reference path `config.json`
  declares actually exists in the image.
- Added coverage reporting with a ratcheting floor and an 80% target.
- `make check-all` now runs in seconds and works on a fresh clone; it previously
  required a 1.1 GB download and a Docker daemon.

## 2.0.4

- Migrated all logging to per-module loggers (`logging.getLogger(__name__)`); log
  records now carry the emitting module name instead of `root`.
- Replaced deprecated `datetime.utcnow()` calls; emitted timestamps are unchanged.
- Ruff lint rules are now selected explicitly instead of extending ruff's defaults, so
  lint results no longer change when ruff widens its default rule set.
- mypy configuration moved from Makefile flags into `pyproject.toml`.
- Fixed `.gitignore` ignoring the tracked `docs/` directory, which silently hid newly
  added documentation pages from Git.
- Fixed `make test-quiet`, which used a `--log-cli=false` flag pytest does not accept.
- Added `AGENTS.md` and `CLAUDE.md` with repository instructions for coding agents.

## 2.0.3

- Raised Kestrel max align and haplotype states to 40, fixing GDP inflation
  ([#156](https://github.com/hassansaei/VNtyper/issues/156)).
- Handled `pd.NA` and `None` in flagging condition evaluation to prevent `TypeError`s.
- Added a GDP inflation guard test with an anonymized hg38 sample.
- Added `CITATION.cff` with the Zenodo concept DOI.

## 2.0.2

- Removed specific motifs from `exclude_motifs_right` to prevent unwanted motif filtering.
- Clarified depth score conditions in confidence assignment
  ([#154](https://github.com/hassansaei/VNtyper/issues/154)).
- Added a dynamic version variable to the docs via `mkdocs-macros-plugin`.
- Modernized docs deployment to use GitHub Actions Pages.
- Added a CFLAGS workaround for adVNTR compilation errors on GCC 14+.

## 2.0.1

- Disabled duplicate flagging in Kestrel configuration.
- Cleared `motifs_for_alt_gg` array to prevent unintended variant filtering.
- Fixed flagging to occur before variant selection, preventing selection of flagged variants ([#145](https://github.com/hassansaei/VNtyper/issues/145)).
- Fixed confidence assignment to filter sub-threshold variants at the root rather than via override.
- Raised `Low_Depth_Conserved_Motifs` threshold from 0.2 to 0.4.
- Updated variant selection to prioritize depth score before haplotype count.
- Extracted frameshift deduplication into a helper function for clarity.
- Added comprehensive unit tests for the flagging module.

## 2.0.0

A complete rewrite and modernization of the VNtyper 2 pipeline.

- **Modern Python packaging** using `pyproject.toml` (PEP 517/518/621).
- **Refactored Kestrel postprocessing** with configurable thresholds via `kestrel_config.json`.
- **HTML reports** with embedded IGV integration for interactive variant inspection.
- **Cohort analysis** command with built-in pseudonymization for multi-sample studies.
- **Docker multi-stage build** for reproducible, lightweight container images.
- **Comprehensive test suite** with unit and integration tests, including downloadable test data.
- **Multiple reference assemblies** supported: hg19, hg38, GRCh37, and GRCh38.
- **Modular architecture** separating variant parsing, motif processing, scoring, flagging, and confidence assignment into individual modules.
- **Snakemake workflow** for scalable batch processing.

## 1.x

The original VNtyper release as described in Saei et al., *iScience* 2023. This version provided the initial alignment-free genotyping approach for MUC1 VNTR using Kestrel and adVNTR.
