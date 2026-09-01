"""Complete immutable schema-3 run fixtures for calibration unit tests."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from vntyper.modules.advntr.artifact_evidence import load_packaged_artifact_evidence
from vntyper.scripts.canonical_json import canonical_json_bytes
from vntyper.scripts.decision_profile import load_packaged_decision_profile
from vntyper.scripts.profile_provenance import profile_summary_fields

IDENTITY = "MUC1-X-60-coding-v1|60|59|-|C"
OTHER_IDENTITY = "MUC1-X-60-coding-v1|59|58|-|G"

CORE_ARTIFACTS = (
    "pipeline_summary.json",
    "provenance/decision_profile.json",
    "kestrel/kestrel_pre_result.tsv",
    "kestrel/bam_identity_replay.v1.json",
    "kestrel/kestrel_result.tsv",
)
ADVNTR_ARTIFACTS = (
    "provenance/advntr_artifact_evidence.json",
    "advntr/output_adVNTR_result.tsv",
)


def write_schema_three_run(
    root: Path,
    key: str,
    *,
    identity: str = IDENTITY,
    name: str = "59dupC",
    support: int = 10,
    with_advntr: bool = False,
    no_kestrel_finding: bool = False,
    retain_pre_result: bool = False,
) -> dict[str, object]:
    """Write one complete current-run artifact set and its immutable declaration."""
    kestrel = root / "kestrel"
    provenance = root / "provenance"
    kestrel.mkdir(parents=True)
    provenance.mkdir(parents=True)
    profile = load_packaged_decision_profile()
    (provenance / "decision_profile.json").write_bytes(profile.canonical_bytes)

    raw_key = json.dumps(
        {"source": "kestrel", "values": ["X-X", 67, "G", "GG"]},
        separators=(",", ":"),
        sort_keys=True,
    )
    capture = {
        "__Identity_Raw_Representation_Key": raw_key,
        "__Identity_Molecular_Identity": identity,
        "__Identity_Translation_Status": "resolved",
        "__Identity_Translation_Failure": "absent",
        "__Identity_Context_Diverges": "false",
        "__Identity_Observation_Ordinal": "0",
    }
    pre_row = {
        "Motifs": "X-X",
        "POS": "67",
        "REF": "G",
        "ALT": "GG",
        "Motif": "X",
        "Motif_sequence": "A" * 120,
        "Estimated_Depth_AlternateVariant": str(support),
        "Estimated_Depth_Variant_ActiveRegion": "1000",
        "Depth_Score": "0.01",
        "Confidence": "High_Precision",
        "is_frameshift": "True",
        "is_valid_frameshift": "True",
        "depth_confidence_pass": "True",
        "alt_filter_pass": "True",
        "motif_filter_pass": "True",
        "flag_filter_pass": "True",
        "Molecular_Identity": identity,
        "Molecular_Identity_Translation_Status": "resolved",
        "Molecular_Identity_Translation_Failure": "",
        "Molecular_Identity_Context_Diverges": "False",
        **capture,
    }
    pre_rows = [pre_row] if not no_kestrel_finding or retain_pre_result else []
    _write_tsv(kestrel / "kestrel_pre_result.tsv", pre_rows, list(pre_row))

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
        "Motifs": "X-X",
        "POS": "67",
        "REF": "G",
        "ALT": "GG",
        "Motif": "X",
        "Estimated_Depth_AlternateVariant": str(support),
        "Estimated_Depth_Variant_ActiveRegion": "1000",
        "Depth_Score": "0.01",
        "Confidence": "High_Precision",
        "Flag": "Not flagged",
        "Nomenclature": name,
        "Nomenclature_Tier": "A",
        "Nomenclature_Flags": "",
        "Nomenclature_Kestrel": name,
        "Nomenclature_adVNTR": "",
        "Molecular_Identity": identity,
        "Molecular_Identity_Status": "unique",
        "Equivalent_Representation_Count": "1",
        "Identity_Hypothesis_Count": "1",
        **capture,
        **selection,
    }
    result = kestrel / "kestrel_result.tsv"
    final_rows = [] if no_kestrel_finding else [final_row]
    _write_tsv(result, final_rows, list(final_row))

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
                            "identities": [identity],
                            "candidate_observation_ordinals": [[0]],
                            "minimum_kmer_depth": depth,
                        }
                        for depth in (8, 12)
                    ],
                    "counts": [{"identity": identity, "record_count": 2}],
                },
            }
        ]
        if pre_rows
        else [],
    }
    (kestrel / "bam_identity_replay.v1.json").write_bytes(canonical_json_bytes(replay))

    result_bytes = result.read_bytes()
    advntr_step: dict[str, object] | None = None
    evidence_digest: str | None = None
    artifact_names = list(CORE_ARTIFACTS)
    if with_advntr:
        evidence = load_packaged_artifact_evidence()
        evidence_path = provenance / "advntr_artifact_evidence.json"
        evidence_path.write_bytes(evidence.canonical_bytes)
        advntr_row = {
            "VID": "25561",
            "Variant": "NOT_GOVERNED",
            "NumberOfSupportingReads": "7",
            "Pvalue": "0.001",
            "Coverage": "42",
            "Evidence_Disposition": "admissible",
            "Nomenclature": name,
            "Nomenclature_Tier": "A",
            "Molecular_Identity": identity,
            "Molecular_Identity_Status": "unique",
            "Equivalent_Representation_Count": "1",
            "Identity_Hypothesis_Count": "1",
        }
        advntr_result = root / "advntr" / "output_adVNTR_result.tsv"
        advntr_result.parent.mkdir()
        _write_tsv(advntr_result, [advntr_row])
        advntr_bytes = advntr_result.read_bytes()
        advntr_step = {
            "step": "adVNTR Genotyping",
            "result_file": f"/untrusted/operator/path/{key}/advntr.tsv",
            "file_type": "tsv",
            "md5sum": hashlib.md5(advntr_bytes).hexdigest(),
            "parsed_result": {"comments": [], "data": [advntr_row]},
        }
        evidence_digest = evidence.digest
        artifact_names.extend(ADVNTR_ARTIFACTS)

    steps = [
        {
            "step": "Kestrel Genotyping",
            "result_file": f"/untrusted/operator/path/{key}/result.tsv",
            "file_type": "tsv",
            "md5sum": hashlib.md5(result_bytes).hexdigest(),
            "parsed_result": {"comments": [], "data": final_rows},
        }
    ]
    if advntr_step is not None:
        steps.append(advntr_step)
    summary = {
        "schema_version": 3,
        "decision_policy": "legacy-selection-v1",
        "advntr_evidence_digest": evidence_digest,
        **profile_summary_fields(profile),
        "version": "2.0.0",
        "reference_assembly_requested": "hg19",
        "steps": steps,
    }
    (root / "pipeline_summary.json").write_bytes(canonical_json_bytes(summary))
    hashes = {relative: hashlib.sha256((root / relative).read_bytes()).hexdigest() for relative in artifact_names}
    return {"root": str(root), "artifacts": hashes}


def refresh_run_hashes(entry: dict[str, object]) -> None:
    """Refresh a fixture declaration after an intentional test mutation."""
    root = Path(str(entry["root"]))
    artifacts = entry["artifacts"]
    assert isinstance(artifacts, dict)
    for relative in tuple(artifacts):
        artifacts[relative] = hashlib.sha256((root / relative).read_bytes()).hexdigest()


def _write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str] | None = None) -> None:
    fields = list(rows[0]) if fields is None else fields
    lines = ["\t".join(fields)]
    lines.extend("\t".join(row[field] for field in fields) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
