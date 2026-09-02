"""Independent development-corpus snapshot for calibration golden checks.

This module may use the separately guarded truth oracle and Python's standard
library only. It must never import VNtyper decisions, canonicalizers, codecs,
grouping rules, calibration predicates, or profile resolution.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from tests.golden import identity_oracle
from tests.golden.identity_oracle import DisplayCounts, GoldenCorpus
from tests.golden.oracle_import_guard import assert_independent_import_closure as assert_oracle_import_closure


@dataclass(frozen=True)
class DevelopmentCalibrationSnapshot:
    """Literal pre-fit projection and truthful evidence eligibility state."""

    sim_root: Path
    advntr_root: Path
    mutated_samples: int
    control_samples: int
    public_identity_rows: int
    selected_locus_rows: int
    total: DisplayCounts
    by_tier: dict[str, DisplayCounts]
    control_findings: int
    evidence_role: str
    eligible_for_independent_validation: bool
    eligible_for_locked_evaluate: bool
    ineligibility_reason: str


@dataclass(frozen=True)
class SourceArtifactBinding:
    """Exact external artifact consumed by the development-only bridge."""

    root_kind: str
    relative_path: str
    sha256: str


@dataclass(frozen=True)
class DevelopmentFixtureMember:
    """Literal source facts materialized as one current schema-3 run."""

    key: str
    role: str
    pair_id: str
    source_summary_schema: int
    motif_context: str
    pair_context: str
    support: int
    active_region_depth: int
    depth_score: float
    confidence: str
    flag: str
    tier: str


@dataclass(frozen=True)
class DevelopmentCalibrationFixture:
    """Paths and independent facts for one minimal development-only study."""

    truth_path: Path
    study_path: Path
    runs_path: Path
    classification_path: Path
    members: tuple[DevelopmentFixtureMember, ...]
    source_bindings: tuple[SourceArtifactBinding, ...]


_PROFILE_SHA256 = "be6329fb12107a1b6b65e425257be6233c7e2115e299e941c12a63a6a6d59718"
_IDENTITY = "MUC1-X-60-coding-v1|60|59|-|C"
_DISPLAY_NAME = "59dupC"
_SOURCE_MEMBERS = (
    ("held-pair-3004", "locked-heldout", "pair_3004"),
    ("select-pair-3001", "policy-selection", "pair_3001"),
    ("train-pair-3000", "training", "pair_3000"),
    ("validate-pair-3003", "validation", "pair_3003"),
)
_GROUP_NAMESPACES = (
    "individual-family",
    "simulated-pair",
    "backbone-seed-lineage",
    "replicate-rerun",
    "depth-series-source",
    "batch",
    "repeat-context",
)
_CORE_ARTIFACTS = (
    "pipeline_summary.json",
    "provenance/decision_profile.json",
    "kestrel/kestrel_pre_result.tsv",
    "kestrel/bam_identity_replay.v1.json",
    "kestrel/kestrel_result.tsv",
)
_ROOT_ANCHORS = {
    "VNTYPER_SIM_ROOT": ("experiment1_dupC/ground_truth.csv",),
    "VNTYPER_ADVNTR_ROOT": (
        "experiment1_dupC/pair_3000/mutated/pipeline_summary.json",
        "experiment1_dupC/pair_3000/mutated/kestrel/kestrel_result.tsv",
    ),
}


def require_explicit_roots(environment: Mapping[str, str]) -> tuple[Path, Path]:
    """Return both explicit roots or fail the golden tier without skipping.

    Args:
        environment: Environment mapping supplying the two required variables.

    Returns:
        The simulation and paired adVNTR roots in their declared order.

    Raises:
        AssertionError: If either variable is absent or either root is missing or incomplete.
    """
    values: list[Path] = []
    for name in ("VNTYPER_SIM_ROOT", "VNTYPER_ADVNTR_ROOT"):
        value = environment.get(name)
        if not value:
            raise AssertionError(f"{name} must name an explicit golden corpus root; skips are forbidden")
        root = Path(value).resolve()
        if not root.is_dir():
            raise AssertionError(f"{name} golden corpus root not found: {root}")
        values.append(root)
    for name, root in zip(("VNTYPER_SIM_ROOT", "VNTYPER_ADVNTR_ROOT"), values, strict=True):
        missing = [relative for relative in _ROOT_ANCHORS[name] if not (root / relative).is_file()]
        if missing:
            raise AssertionError(f"{name} golden corpus root is incomplete; missing: {', '.join(missing)}")
    return values[0], values[1]


def load_development_snapshot(sim_root: Path, advntr_root: Path) -> DevelopmentCalibrationSnapshot:
    """Load both roots and freeze the exact shipped pre-fit projection."""
    corpus = identity_oracle.load_golden_corpus(sim_root, advntr_root)
    return snapshot_from_corpus(corpus)


def snapshot_from_corpus(corpus: GoldenCorpus) -> DevelopmentCalibrationSnapshot:
    """Project only independently observed corpus facts, never fitted decisions."""
    if not isinstance(corpus, GoldenCorpus):
        raise AssertionError("calibration oracle requires an independently loaded GoldenCorpus")
    public_rows = sum(len(rows) for rows in corpus.identity_projection_by_sample.values())
    selected_rows = sum(len(rows) for rows in corpus.selected_projection_by_sample.values())
    return DevelopmentCalibrationSnapshot(
        sim_root=corpus.sim_root,
        advntr_root=corpus.advntr_root,
        mutated_samples=corpus.mutated_samples,
        control_samples=corpus.control_samples,
        public_identity_rows=public_rows,
        selected_locus_rows=selected_rows,
        total=corpus.total,
        by_tier=dict(corpus.by_tier),
        control_findings=corpus.control_findings,
        evidence_role="previously-examined-development-simulation",
        eligible_for_independent_validation=False,
        eligible_for_locked_evaluate=False,
        ineligibility_reason=(
            "The simulations and paired caller outputs were previously examined development evidence; "
            "they are neither independent external validation nor a custodian-locked held-out cohort."
        ),
    )


def assert_independent_import_closure(repo_root: Path) -> frozenset[Path]:
    """Reject direct, recursive, or dynamic production imports in this oracle."""
    scanned = assert_oracle_import_closure(Path(__file__), repo_root)
    return frozenset(scanned)


def materialize_development_fixture(
    sim_root: Path,
    advntr_root: Path,
    repo_root: Path,
    output_root: Path,
) -> DevelopmentCalibrationFixture:
    """Bridge literal schema-2 corpus facts into a minimal schema-3 study.

    The historical roots are read-only and are never described as current runs.
    Four synthetic current-run envelopes are written under ``output_root`` from
    independently checked truth and paired public-result facts.

    Args:
        sim_root: Explicit simulation corpus root.
        advntr_root: Explicit paired caller corpus root.
        repo_root: Checkout containing the byte-pinned packaged profile.
        output_root: Empty test-owned directory for generated inputs.

    Returns:
        Paths, source hashes, and literal member facts for CLI assertions.

    Raises:
        AssertionError: If a source is incomplete, changed, or not historical schema 2.
    """
    corpus = identity_oracle.load_golden_corpus(sim_root, advntr_root)
    output_root.mkdir(parents=True, exist_ok=False)
    profile_path = repo_root / "vntyper" / "profiles" / "decision_profile.json"
    profile_bytes = _require_bytes(profile_path)
    if hashlib.sha256(profile_bytes).hexdigest() != _PROFILE_SHA256:
        raise AssertionError("packaged profile differs from the independently pinned development baseline")
    profile = _json_object(profile_bytes, "packaged profile")

    truth_relative = "experiment1_dupC/ground_truth.csv"
    truth_path = sim_root / truth_relative
    bindings = [SourceArtifactBinding("simulation", truth_relative, _sha256_file(truth_path))]
    truth_rows = _read_csv(truth_path)
    members = []
    run_declarations: dict[str, object] = {}
    labels = []
    partition_members = []
    for key, role, pair_id in _SOURCE_MEMBERS:
        source_key = f"experiment1_dupC/{pair_id}"
        expectation = corpus.expected_by_sample.get(source_key)
        if expectation is None or expectation.mutation != "dupC":
            raise AssertionError(f"selected development source is not the expected dupC truth member: {source_key}")
        if expectation.expected.identity != _IDENTITY or expectation.expected.name != _DISPLAY_NAME:
            raise AssertionError(f"selected development truth changed for {source_key}")
        matching_truth = [
            row
            for row in truth_rows
            if row.get("pair_id") == pair_id and row.get("condition") == "mutated" and row.get("mutation") == "dupC"
        ]
        if len(matching_truth) != 1:
            raise AssertionError(f"selected simulation truth row is missing or ambiguous: {source_key}")

        source_base = Path("experiment1_dupC") / pair_id / "mutated"
        summary_relative = (source_base / "pipeline_summary.json").as_posix()
        result_relative = (source_base / "kestrel" / "kestrel_result.tsv").as_posix()
        summary_path = advntr_root / summary_relative
        result_path = advntr_root / result_relative
        summary = _json_object(_require_bytes(summary_path), f"historical summary {source_key}")
        schema = summary.get("schema_version")
        if schema != 2:
            raise AssertionError(f"development bridge source must remain historical schema 2: {source_key}")
        source_rows = _read_tsv(result_path)
        positive = [row for row in source_rows if row.get("Confidence") != "Negative"]
        if len(positive) != 1:
            raise AssertionError(f"development bridge requires one literal selected Kestrel row: {source_key}")
        row = positive[0]
        if row.get("Nomenclature") != _DISPLAY_NAME or row.get("POS") != "67":
            raise AssertionError(f"development source verdict changed for {source_key}")
        member = DevelopmentFixtureMember(
            key=key,
            role=role,
            pair_id=pair_id,
            source_summary_schema=schema,
            motif_context=_required_text(row.get("Motif"), f"{source_key} motif"),
            pair_context=_required_text(row.get("Motifs"), f"{source_key} motif pair"),
            support=_integer_text(row.get("Estimated_Depth_AlternateVariant"), f"{source_key} support"),
            active_region_depth=_integer_text(
                row.get("Estimated_Depth_Variant_ActiveRegion"), f"{source_key} active-region depth"
            ),
            depth_score=_float_text(row.get("Depth_Score"), f"{source_key} depth score"),
            confidence=_required_text(row.get("Confidence"), f"{source_key} confidence"),
            flag=_required_text(row.get("Flag"), f"{source_key} flag"),
            tier=_required_text(row.get("Nomenclature_Tier"), f"{source_key} tier"),
        )
        members.append(member)
        bindings.extend(
            (
                SourceArtifactBinding("advntr", summary_relative, _sha256_file(summary_path)),
                SourceArtifactBinding("advntr", result_relative, _sha256_file(result_path)),
            )
        )
        run_declarations[key] = _write_schema_three_run(output_root / "runs" / key, member, profile, profile_bytes)
        if role != "locked-heldout":
            labels.append(
                {
                    "label_key": f"label-{key}",
                    "manifest_key": key,
                    "truth_status": "mutated",
                    "expected_identity": _IDENTITY,
                    "expected_display_name": _DISPLAY_NAME,
                    "mutation_class": "duplication",
                }
            )
        partition_members.append(
            {
                "key": key,
                "role": role,
                "provenance": "external-custodian" if role == "locked-heldout" else "development",
                "assay_class": "development-schema2-bridge",
                "groups": {namespace: [f"{namespace}:{key}"] for namespace in _GROUP_NAMESPACES},
            }
        )

    truth_document = {
        "schema_version": "calibration-truth-v1",
        "labels": {"schema_version": "calibration-labels-v1", "rows": labels},
    }
    study_document = {
        "schema_version": "calibration-study-v1",
        "protocol": {
            "objective": "lexicographic-safety-v1",
            "bootstrap_iterations": 10000,
            "bootstrap_interval": "percentile",
            "multiplicity_method": "holm",
            "seed": 295,
            "maximum_free_parameters": 4,
            "minimum_stratum_count": 1,
            "maximum_abstention_fraction": 0.25,
            "assay_classes": ["development-schema2-bridge"],
            "mutation_classes": ["duplication"],
            "candidate_grid": {
                "minimum_record_count_margin": [1],
                "minimum_record_share": [0.5],
                "minimum_record_share_margin": [0.0],
                "xd_veto": ["disabled"],
            },
        },
        "partitions": {"schema_version": "calibration-partitions-v1", "members": partition_members},
    }
    runs_document = {"schema_version": "calibration-runs-v1", "runs": run_declarations}
    classification = {
        "schema_version": "calibration-development-evidence-classification-v1",
        "classification": "previously-examined-development-simulation",
        "eligible_for_independent_validation": False,
        "eligible_for_locked_evaluate": False,
        "statement": (
            "This development fixture is neither an independent external cohort nor custodian-locked heldout evidence."
        ),
    }
    generated_truth = output_root / "truth.json"
    generated_study = output_root / "study.json"
    generated_runs = output_root / "runs.json"
    classification_path = output_root / "development_evidence_classification.json"
    _write_json(generated_truth, truth_document)
    _write_json(generated_study, study_document)
    _write_json(generated_runs, runs_document)
    _write_json(classification_path, classification)
    return DevelopmentCalibrationFixture(
        generated_truth,
        generated_study,
        generated_runs,
        classification_path,
        tuple(members),
        tuple(bindings),
    )


def verify_checksum_tree(root: Path) -> None:
    """Independently verify one direct-child calibration checksum manifest."""
    document = _json_object(_require_bytes(root / "checksums.json"), "checksum manifest")
    if set(document) != {"schema_version", "files"} or document["schema_version"] != "calibration-checksums-v1":
        raise AssertionError(f"invalid calibration checksum manifest: {root}")
    files = document["files"]
    if not isinstance(files, dict) or not files:
        raise AssertionError(f"empty calibration checksum manifest: {root}")
    actual_names = {path.name for path in root.iterdir() if path.is_file() and path.name != "checksums.json"}
    if set(files) != actual_names:
        raise AssertionError(f"calibration checksum inventory differs: {root}")
    for name, expected in files.items():
        if not isinstance(name, str) or not isinstance(expected, str):
            raise AssertionError(f"invalid calibration checksum entry: {root}")
        if _sha256_file(root / name) != expected:
            raise AssertionError(f"calibration checksum differs: {root / name}")


def verify_exact_child_directories(root: Path, expected: tuple[str, ...]) -> None:
    """Require an exact direct-child directory inventory with no files or extras."""
    actual = {path.name for path in root.iterdir() if path.is_dir()}
    non_directories = {path.name for path in root.iterdir() if not path.is_dir()}
    if actual != set(expected) or non_directories:
        raise AssertionError(f"calibration child-directory inventory differs: {root}")


def _write_schema_three_run(
    root: Path,
    member: DevelopmentFixtureMember,
    profile: dict[str, object],
    profile_bytes: bytes,
) -> dict[str, object]:
    kestrel = root / "kestrel"
    provenance = root / "provenance"
    kestrel.mkdir(parents=True)
    provenance.mkdir(parents=True)
    (provenance / "decision_profile.json").write_bytes(profile_bytes)
    raw_key = json.dumps(
        {"source": "kestrel", "values": [member.pair_context, 67, "G", "GG"]},
        separators=(",", ":"),
        sort_keys=True,
    )
    capture = {
        "__Identity_Raw_Representation_Key": raw_key,
        "__Identity_Molecular_Identity": _IDENTITY,
        "__Identity_Translation_Status": "resolved",
        "__Identity_Translation_Failure": "absent",
        "__Identity_Context_Diverges": "false",
        "__Identity_Observation_Ordinal": "0",
    }
    pre_row = {
        "Motifs": member.pair_context,
        "POS": "67",
        "REF": "G",
        "ALT": "GG",
        "Motif": member.motif_context,
        "Motif_sequence": "A" * 120,
        "Estimated_Depth_AlternateVariant": str(member.support),
        "Estimated_Depth_Variant_ActiveRegion": str(member.active_region_depth),
        "Depth_Score": str(member.depth_score),
        "Confidence": member.confidence,
        "is_frameshift": "True",
        "is_valid_frameshift": "True",
        "depth_confidence_pass": "True",
        "alt_filter_pass": "True",
        "motif_filter_pass": "True",
        "flag_filter_pass": "True",
        "Molecular_Identity": _IDENTITY,
        "Molecular_Identity_Translation_Status": "resolved",
        "Molecular_Identity_Translation_Failure": "",
        "Molecular_Identity_Context_Diverges": "False",
        **capture,
    }
    _write_tsv(kestrel / "kestrel_pre_result.tsv", [pre_row])
    selection = {
        "__Identity_Selected_Raw_Representation_Key": raw_key,
        "__Identity_Equivalent_Representation_Count": "1",
        "__Identity_Hypothesis_Count": "1",
        "__Identity_Group_Blocking_Gates": "[]",
        "__Identity_Group_Flags": "[]",
        "__Identity_Selected_Observation_Ordinal": "0",
        "__Identity_Group_Context_Diverges": "false",
    }
    final_row = {
        "Motifs": member.pair_context,
        "POS": "67",
        "REF": "G",
        "ALT": "GG",
        "Motif": member.motif_context,
        "Estimated_Depth_AlternateVariant": str(member.support),
        "Estimated_Depth_Variant_ActiveRegion": str(member.active_region_depth),
        "Depth_Score": str(member.depth_score),
        "Confidence": member.confidence,
        "Flag": member.flag,
        "Nomenclature": _DISPLAY_NAME,
        "Nomenclature_Tier": member.tier,
        "Nomenclature_Flags": "",
        "Nomenclature_Kestrel": _DISPLAY_NAME,
        "Nomenclature_adVNTR": "",
        "Molecular_Identity": _IDENTITY,
        "Molecular_Identity_Status": "unique",
        "Equivalent_Representation_Count": "1",
        "Identity_Hypothesis_Count": "1",
        **capture,
        **selection,
    }
    result_path = kestrel / "kestrel_result.tsv"
    _write_tsv(result_path, [final_row])
    replay = {
        "schema_version": "bam-identity-replay-v1",
        "loci": [
            {
                "candidate_observation_ordinals": [0],
                "state": "observed",
                "evidence": {
                    "eligible_record_count": 2,
                    "records": [
                        {
                            "identities": [_IDENTITY],
                            "candidate_observation_ordinals": [[0]],
                            "minimum_kmer_depth": depth,
                        }
                        for depth in (8, 12)
                    ],
                    "counts": [{"identity": _IDENTITY, "record_count": 2}],
                },
            }
        ],
    }
    _write_json(kestrel / "bam_identity_replay.v1.json", replay)
    result_bytes = result_path.read_bytes()
    summary = {
        "schema_version": 3,
        "decision_policy": "legacy-selection-v1",
        "advntr_evidence_digest": None,
        "decision_profile_id": profile["profile_id"],
        "decision_profile_revision": profile["profile_revision"],
        "decision_profile_kind": profile["profile_kind"],
        "decision_profile_source": "package",
        "decision_profile_sha256": _PROFILE_SHA256,
        "decision_profile_snapshot": "provenance/decision_profile.json",
        "version": "2.0.26",
        "reference_assembly_requested": "hg19",
        "steps": [
            {
                "step": "Kestrel Genotyping",
                "result_file": f"/historical-schema2-bridge/{member.key}/kestrel_result.tsv",
                "file_type": "tsv",
                "md5sum": hashlib.md5(result_bytes, usedforsecurity=False).hexdigest(),
                "parsed_result": {"comments": [], "data": [final_row]},
            }
        ],
    }
    _write_json(root / "pipeline_summary.json", summary)
    artifacts = {relative: _sha256_file(root / relative) for relative in _CORE_ARTIFACTS}
    return {"member_key": member.key, "root": f"runs/{member.key}", "artifacts": artifacts}


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fields = list(rows[0])
    lines = ["\t".join(fields)]
    lines.extend("\t".join(row[field] for field in fields) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(_canonical_json_bytes(value))


def _canonical_json_bytes(value: object) -> bytes:
    return (json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True) + "\n").encode("utf-8")


def _json_object(raw: bytes, label: str) -> dict[str, object]:
    try:
        value = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise AssertionError(f"invalid {label}") from error
    if not isinstance(value, dict):
        raise AssertionError(f"{label} must be an object")
    return value


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def _required_text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise AssertionError(f"{label} must be a non-empty string")
    return value


def _integer_text(value: object, label: str) -> int:
    text = _required_text(value, label)
    if not text.isdigit():
        raise AssertionError(f"{label} must be a non-negative integer")
    return int(text)


def _float_text(value: object, label: str) -> float:
    text = _required_text(value, label)
    try:
        number = float(text)
    except ValueError as error:
        raise AssertionError(f"{label} must be numeric") from error
    if not 0 <= number < float("inf"):
        raise AssertionError(f"{label} must be finite and non-negative")
    return number


def _require_bytes(path: Path) -> bytes:
    try:
        return path.read_bytes()
    except OSError as error:
        raise AssertionError(f"required golden artifact not found: {path}") from error


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(_require_bytes(path)).hexdigest()
