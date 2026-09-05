# Repository Quality, Validation & CI Harnesses (`scripts/`)

This directory contains repository-level engineering infrastructure, validation instruments,
testing harnesses, release controllers, and fixture generators.

> [!IMPORTANT]
> **Pipeline source code lives under [`vntyper/scripts/`](../vntyper/scripts/README.md), not here.**
> The Python package's operational stages (such as `pipeline.py`, `fastq_bam_processing.py`,
> `kestrel_genotyping.py`, and `scoring.py`) execute the MUC1 genotyping workflow.
> In contrast, the files in this root `scripts/` directory are **authoritative quality instruments**
> covered by unit tests, type checking (`mypy`), and the 86% branch-inclusive test coverage floor.

---

## Complete Script Inventory

### 1. Quality, Coverage & Compatibility Gates

- **`coverage_gate.py`**: Enforces the hard 86% branch-inclusive test coverage floor and ratchet across the repository (`make test-unit-cov`).
- **`check_integration_compatibility.py`**: Command-line entry point for the backward compatibility gate run in CI on every PR.
- **`integration_compatibility.py`**: Validation engine enforcing fail-closed backward compatibility against immutable historical baselines.
- **`integration_compatibility_observations.py`**: Append-only, versioned ledger recording verified pipeline output observations.
- **`ble001_policy.py`**: AST analyzer verifying that all repository broad exception handlers (`except Exception:`) adhere to the reviewed policy inventory.
- **`ble001_policy_validation.py`**: Pure path and test-node validation logic for the BLE001 policy adapter.
- **`ble001_policy.json`**: Machine-readable inventory of approved exception handler locations, boundaries, and rationales.

### 2. Golden Cohort Regression Gate (`golden_cohort_gate.py` & `golden_cohort/`)

- **`golden_cohort_gate.py`**: Master entry point for reproducing the golden-cohort regression gate (`make golden-gate`).
- **`golden_cohort/__init__.py`**: Subpackage initialization and exported symbols.
- **`golden_cohort/admissibility.py`**: Verifies precondition checks (clean worktree, branch validity) before a gate run is admissible as evidence.
- **`golden_cohort/artifacts.py`**: Disk I/O reader that extracts and parses pipeline output artifacts for comparison.
- **`golden_cohort/case_expectations.py`**: Pure outcome declarations and expected outputs for each golden cohort case.
- **`golden_cohort/cohort_cases.py`**: Cohort case declarations and metadata mapping for the test matrix.
- **`golden_cohort/compare.py`**: Compares test outputs against baselines, diffing artifacts side-by-side with structured reporting.
- **`golden_cohort/cram_cases.py`**: CRAM-specific case policy and fixture-path derivation.
- **`golden_cohort/cram_evidence.py`**: Fail-closed decision logic over recorded CRAM read-set evidence.
- **`golden_cohort/launcher.py`**: Subprocess launcher proving exact worktree execution without path cross-contamination.
- **`golden_cohort/matrix.py`**: Case matrix builder dynamically discovering test datasets available in `tests/data`.
- **`golden_cohort/normalise.py`**: Normalizes dynamic values (timestamps, temporary paths, run durations) for bitwise artifact comparison.
- **`golden_cohort/read_set_commands.py`**: Memory-bounded command builder for extracting CRAM read-set evidence via samtools.
- **`golden_cohort/read_sets.py`**: Read-set data structures and equivalence checks for CRAM inputs.
- **`golden_cohort/runner.py`**: Test orchestrator executing comparison sides sequentially or concurrently.
- **`golden_cohort/waiver.py`**: Evaluates, validates, and records waivers for intentional, reviewed output deltas.

### 3. Advisory Mutation Testing Framework (`mutation_*.py`)

- **`mutation_test.py`**: Mutation scoring engine and CLI for VNtyper's decision logic (`make mutation`, `make mutation-render`).
- **`mutation_guard.py`**: Pre-flight safety checks preventing mutation runs on dirty or unrestored git worktrees.
- **`mutation_output.py`**: Atomic installation of mutation results and automated generation of `docs/development/mutation-testing.md`.
- **`mutation_provenance.py`**: Tracks mutant generation parameters, commit anchors, and test results.
- **`mutation_source.py`**: AST-based mutant generation, code injection, and syntax verification.
- **`mutation_workspace.py`**: Orchestrates isolated worktree creation for running mutated test suites safely.
- **`mutation_workspace_fs.py`**: Filesystem copy and cleanup operations for temporary mutation sandboxes.
- **`mutation_workspace_git.py`**: Git worktree branching and restoration logic for zero-footprint mutation testing.

### 4. Release Automation & Registry Controllers

- **`release_policy.py`**: Enforces strict semantic versioning, git tag consistency, and required CI check completion prior to release promotion.
- **`release_manifest.py`**: Inspects and parses OCI image manifest descriptors from GHCR.
- **`release_evidence.py`**: Verifies exact SHA-256 image digests and generates cryptographic release attestations.
- **`release_registry.py`**: Fail-closed error classifier for GHCR manifest queries and token authentications.
- **`pypi_environment_contract.py`**: Validates the live GitHub PyPI release environment configuration and security contracts.
- **`bundle_release.py`**: Builds, verifies with SHA-256, and packages the authoritative MUC1 reference bundle (`make bundle-release`).

### 5. Test Fixtures & Reference Infrastructure

- **`download_test_data.py`**: CLI utility that downloads, unpacks, and validates the 1.1 GB Zenodo test data archive (`make download-data`).
- **`make_cram_fixtures.py`**: Derives lossless, synthetic CRAM fixtures from BAM test cohorts.
- **`cram_fixture_selection.py`**: Selects candidate BAM files eligible for CRAM fixture generation.
- **`cram_reference_contract.py`**: Validates reference-compressed CRAM fixture invariants and URI schemes.
- **`single_end_fixture.py`**: Derives single-end BAM test fixtures from paired-end BAM inputs for mixed-layout testing.
- **`advntr_len_differential.py`**: Diagnostic utility verifying adVNTR repeat-unit `LEN` token calculations.
- **`render_decision_profile.py`**: Compiles and verifies the packaged calibration decision profile from legacy parameters.
- **`ucsc_hg19_primary_contigs.tsv`**: Canonical contig ordering and length table for UCSC hg19 primary assemblies.
