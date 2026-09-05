# Repository Quality, Validation & CI Harnesses (`scripts/`)

This directory contains repository-level engineering infrastructure, validation instruments,
testing harnesses, and CI/CD policy controllers.

> [!NOTE]
> **Pipeline source code lives under `vntyper/scripts/`, not here.**
> The Python package's operational stages (such as `pipeline.py`, `fastq_bam_processing.py`,
> `kestrel_genotyping.py`, and `scoring.py`) live in [`vntyper/scripts/`](../vntyper/scripts/README.md).
> The files in this root `scripts/` directory are **authoritative quality instruments** covered by
> unit tests, type checking (`mypy`), and the 86% branch-inclusive coverage gate.

---

## Directory Overview

### 1. Quality & Compatibility Gates
- **`check_integration_compatibility.py` / `integration_compatibility.py`**: Automated gate enforcing fail-closed backward compatibility against historical real-data success baselines.
- **`integration_compatibility_observations.py`**: Append-only versioned report observation ledger.
- **`coverage_gate.py`**: Enforces the hard 86% branch-inclusive test coverage floor and ratchet.
- **`ble001_policy.py` / `ble001_policy_validation.py` / `ble001_policy.json`**: AST-based verification of repository broad-exception handlers against the reviewed policy inventory.

### 2. Golden Cohort Regression Gate
- **`golden_cohort_gate.py`**: Entry point for the golden cohort validation gate (`make golden-gate`).
- **`golden_cohort/`**: Subpackage implementing the multi-case test matrix, CRAM/BAM evidence extraction, outcome normalization, artifact hashing, and delta waivers (`admissibility.py`, `runner.py`, `compare.py`, `matrix.py`, `waiver.py`, etc.).

### 3. Advisory Mutation Testing
- **`mutation_test.py`**: Advisory mutation testing engine (`make mutation`, `make mutation-render`).
- **`mutation_workspace.py` / `mutation_workspace_fs.py` / `mutation_workspace_git.py`**: Safe worktree isolation for AST mutation injection.
- **`mutation_guard.py` / `mutation_output.py` / `mutation_provenance.py` / `mutation_source.py`**: Mutant tracking, test run evaluation, and Markdown ledger generation.

### 4. Release Automation & Registry Controllers
- **`release_policy.py`**: Enforces strict semantic versioning, git tag consistency, and required CI check names before release execution.
- **`release_evidence.py` / `release_manifest.py` / `release_registry.py`**: GHCR container digest verification, image alias promotion (`vX.Y.Z`, `latest`, `X.Y`, `X`), and PyPI publishing automation.
- **`pypi_environment_contract.py`**: Validates PyPI release environment assumptions and metadata.
- **`bundle_release.py`**: Builds, verifies, and packages the authoritative MUC1 reference bundle (`make bundle-release`).

### 5. Test Data & Fixture Infrastructure
- **`download_test_data.py`**: Downloads and verifies the 1.1 GB Zenodo test dataset (`make download-data`).
- **`cram_reference_contract.py` / `make_cram_fixtures.py` / `cram_fixture_selection.py` / `single_end_fixture.py`**: Generates and validates synthetic CRAM and single-end test fixtures.
- **`advntr_len_differential.py`**: Helper for validating adVNTR repeat-unit length differential calculations.
- **`render_decision_profile.py`**: Validates and renders opt-in calibration profiles.
