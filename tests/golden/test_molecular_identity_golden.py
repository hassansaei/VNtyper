"""Independent molecular-identity and historical Phase-1 golden contracts."""

from __future__ import annotations

import copy
import hashlib
import json
import os
import shutil
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from tests.golden import identity_oracle
from tests.golden.identity_oracle import (
    IDENTITY_COLUMNS,
    SELECTED_PROJECTION_COLUMNS,
    DisplayCounts,
    ExpectedIdentity,
    GoldenCorpus,
    IdentityCounts,
    assert_independent_import_closure,
    load_golden_corpus,
)
from vntyper.modules.advntr.advntr_genotyping import process_advntr_output
from vntyper.scripts import nomenclature_annotate, nomenclature_bam_adapter
from vntyper.scripts.identity_candidate_persistence import (
    IDENTITY_CAPTURE_COLUMNS,
    candidate_capture_cells,
    selected_candidate_cells,
)
from vntyper.scripts.identity_candidates import (
    capture_kestrel_observations,
    overlay_legacy_projection,
    translation_component_from_config,
    with_candidate_evidence,
)
from vntyper.scripts.kestrel_genotyping import (
    FILTER_COLUMNS,
    add_artifact_gate,
    filter_final_dataframe,
)
from vntyper.scripts.molecular_identity_presentation import (
    IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS,
    identity_translation_diagnostic_cells,
)
from vntyper.scripts.nomenclature import Nomenclature
from vntyper.scripts.nomenclature_bam import BamConsensus, BamRescuer
from vntyper.scripts.run_configuration import resolve_run_configuration

pytestmark = pytest.mark.golden

_SIM_ROOT = os.environ.get("VNTYPER_SIM_ROOT")
_ADVNTR_ROOT = os.environ.get("VNTYPER_ADVNTR_ROOT")
if not _SIM_ROOT or not _ADVNTR_ROOT:
    pytest.skip("VNTYPER_SIM_ROOT and VNTYPER_ADVNTR_ROOT benchmark roots are unset", allow_module_level=True)

SIM_ROOT = Path(_SIM_ROOT)
ADVNTR_ROOT = Path(_ADVNTR_ROOT)

PACKAGE_PROFILES = Path(__file__).parents[2] / "vntyper" / "profiles"
PR_B_DECISION_PROJECTION_SHA256 = "338fe05d010f623e794dcf93393904fa13bd8713e2d074c8a5b6c72d6efd96fd"
PACKAGED_DECISION_PROFILE_SHA256 = "be6329fb12107a1b6b65e425257be6233c7e2115e299e941c12a63a6a6d59718"


@dataclass(frozen=True)
class UutReplay:
    """Real reconciliation outputs and directly observed BAM decision counts."""

    corpus: GoldenCorpus
    eligible: int
    bam_fetches: int
    mutated_replays: int
    control_replays: int
    eligible_keys: frozenset[str]
    fetches_by_key: dict[str, int]
    bam_input_symlinks: int


def _plain_json(value: object) -> object:
    """Round-trip a recursively immutable component into plain JSON containers."""
    if isinstance(value, Mapping):
        return {str(key): _plain_json(child) for key, child in value.items()}
    if isinstance(value, tuple):
        return [_plain_json(child) for child in value]
    return value


EXPECTED_CLASS_IDENTITIES = {
    "delGCCCA": ExpectedIdentity("MUC1-X-60-coding-v1|1|5|GCCCA|-", "1_5delGCCCA"),
    "delinsAT": ExpectedIdentity("MUC1-X-60-coding-v1|54|56|CCC|AT", "54_56delinsAT"),
    "dupA": ExpectedIdentity("MUC1-X-60-coding-v1|60|59|-|A", "60dupA"),
    "dupC": ExpectedIdentity("MUC1-X-60-coding-v1|60|59|-|C", "59dupC"),
    "ins25bp": ExpectedIdentity(
        "MUC1-X-60-coding-v1|31|30|-|CAGGCCGGCCCCGGGCTCCGGACAC",
        "30_31insCAGGCCGGCCCCGGGCTCCGGACAC",
    ),
    "insA_pos54": ExpectedIdentity("MUC1-X-60-coding-v1|54|53|-|A", "53_54insA"),
    "insCCCC": ExpectedIdentity("MUC1-X-60-coding-v1|60|59|-|CCCC", "56_59dupCCCC"),
    "insC_pos23": ExpectedIdentity("MUC1-X-60-coding-v1|24|23|-|C", "23dupC"),
    "insG": ExpectedIdentity("MUC1-X-60-coding-v1|59|58|-|G", "58_59insG"),
    "insG_pos54": ExpectedIdentity("MUC1-X-60-coding-v1|54|53|-|G", "53_54insG"),
    "insG_pos58": ExpectedIdentity("MUC1-X-60-coding-v1|58|57|-|G", "57_58insG"),
}

PHASE1_TIER_A_KEYS = frozenset(
    {
        "experiment1_dupC/pair_3000",
        "experiment1_dupC/pair_3001",
        "experiment1_dupC/pair_3003",
        "experiment1_dupC/pair_3005",
        "experiment1_dupC/pair_3007",
        "experiment1_dupC/pair_3008",
        "experiment1_dupC/pair_3011",
        "experiment1_dupC/pair_3012",
        "experiment1_dupC/pair_3014",
        "experiment1_dupC/pair_3017",
        "experiment1_dupC/pair_3020",
        "experiment1_dupC/pair_3023",
        "experiment1_dupC/pair_3025",
        "experiment1_dupC/pair_3026",
        "experiment1_dupC/pair_3028",
        "experiment1_dupC/pair_3030",
        "experiment1_dupC/pair_3031",
        "experiment1_dupC/pair_3033",
        "experiment1_dupC/pair_3034",
        "experiment1_dupC/pair_3036",
        "experiment1_dupC/pair_3037",
        "experiment1_dupC/pair_3038",
        "experiment1_dupC/pair_3039",
        "experiment1_dupC/pair_3040",
        "experiment1_dupC/pair_3043",
        "experiment1_dupC/pair_3044",
        "experiment1_dupC/pair_3047",
        "experiment1_dupC/pair_3051",
        "experiment1_dupC/pair_3052",
        "experiment1_dupC/pair_3054",
        "experiment1_dupC/pair_3055",
        "experiment1_dupC/pair_3057",
        "experiment1_dupC/pair_3059",
        "experiment1_dupC/pair_3065",
        "experiment1_dupC/pair_3068",
        "experiment1_dupC/pair_3069",
        "experiment1_dupC/pair_3073",
        "experiment1_dupC/pair_3074",
        "experiment1_dupC/pair_3075",
        "experiment1_dupC/pair_3076",
        "experiment1_dupC/pair_3077",
        "experiment1_dupC/pair_3078",
        "experiment1_dupC/pair_3086",
        "experiment1_dupC/pair_3087",
        "experiment1_dupC/pair_3088",
        "experiment1_dupC/pair_3089",
        "experiment1_dupC/pair_3091",
        "experiment1_dupC/pair_3092",
        "experiment1_dupC/pair_3093",
        "experiment1_dupC/pair_3094",
        "experiment1_dupC/pair_3095",
        "experiment2_atypical/pair_4045",
        "experiment2_atypical/pair_4083",
    }
)

PHASE1_TIER_B_KEYS = frozenset(
    {
        "experiment1_dupC/pair_3002",
        "experiment1_dupC/pair_3004",
        "experiment1_dupC/pair_3006",
        "experiment1_dupC/pair_3009",
        "experiment1_dupC/pair_3010",
        "experiment1_dupC/pair_3013",
        "experiment1_dupC/pair_3015",
        "experiment1_dupC/pair_3016",
        "experiment1_dupC/pair_3018",
        "experiment1_dupC/pair_3021",
        "experiment1_dupC/pair_3022",
        "experiment1_dupC/pair_3024",
        "experiment1_dupC/pair_3027",
        "experiment1_dupC/pair_3029",
        "experiment1_dupC/pair_3032",
        "experiment1_dupC/pair_3035",
        "experiment1_dupC/pair_3041",
        "experiment1_dupC/pair_3042",
        "experiment1_dupC/pair_3045",
        "experiment1_dupC/pair_3046",
        "experiment1_dupC/pair_3048",
        "experiment1_dupC/pair_3049",
        "experiment1_dupC/pair_3050",
        "experiment1_dupC/pair_3053",
        "experiment1_dupC/pair_3056",
        "experiment1_dupC/pair_3058",
        "experiment1_dupC/pair_3061",
        "experiment1_dupC/pair_3062",
        "experiment1_dupC/pair_3063",
        "experiment1_dupC/pair_3064",
        "experiment1_dupC/pair_3066",
        "experiment1_dupC/pair_3067",
        "experiment1_dupC/pair_3070",
        "experiment1_dupC/pair_3071",
        "experiment1_dupC/pair_3072",
        "experiment1_dupC/pair_3079",
        "experiment1_dupC/pair_3080",
        "experiment1_dupC/pair_3081",
        "experiment1_dupC/pair_3082",
        "experiment1_dupC/pair_3083",
        "experiment1_dupC/pair_3084",
        "experiment1_dupC/pair_3085",
        "experiment1_dupC/pair_3090",
        "experiment1_dupC/pair_3096",
        "experiment1_dupC/pair_3097",
        "experiment1_dupC/pair_3098",
        "experiment1_dupC/pair_3099",
        "experiment2_atypical/pair_4001",
        "experiment2_atypical/pair_4002",
        "experiment2_atypical/pair_4003",
        "experiment2_atypical/pair_4004",
        "experiment2_atypical/pair_4005",
        "experiment2_atypical/pair_4007",
        "experiment2_atypical/pair_4008",
        "experiment2_atypical/pair_4009",
        "experiment2_atypical/pair_4010",
        "experiment2_atypical/pair_4011",
        "experiment2_atypical/pair_4012",
        "experiment2_atypical/pair_4014",
        "experiment2_atypical/pair_4015",
        "experiment2_atypical/pair_4017",
        "experiment2_atypical/pair_4018",
        "experiment2_atypical/pair_4032",
        "experiment2_atypical/pair_4033",
        "experiment2_atypical/pair_4035",
        "experiment2_atypical/pair_4037",
        "experiment2_atypical/pair_4038",
        "experiment2_atypical/pair_4039",
        "experiment2_atypical/pair_4041",
        "experiment2_atypical/pair_4042",
        "experiment2_atypical/pair_4043",
        "experiment2_atypical/pair_4044",
        "experiment2_atypical/pair_4046",
        "experiment2_atypical/pair_4047",
        "experiment2_atypical/pair_4048",
        "experiment2_atypical/pair_4049",
        "experiment2_atypical/pair_4050",
        "experiment2_atypical/pair_4053",
        "experiment2_atypical/pair_4057",
        "experiment2_atypical/pair_4058",
        "experiment2_atypical/pair_4060",
        "experiment2_atypical/pair_4061",
        "experiment2_atypical/pair_4062",
        "experiment2_atypical/pair_4063",
        "experiment2_atypical/pair_4065",
        "experiment2_atypical/pair_4066",
        "experiment2_atypical/pair_4067",
        "experiment2_atypical/pair_4068",
        "experiment2_atypical/pair_4069",
        "experiment2_atypical/pair_4076",
        "experiment2_atypical/pair_4080",
        "experiment2_atypical/pair_4082",
        "experiment2_atypical/pair_4087",
        "experiment2_atypical/pair_4088",
        "experiment2_atypical/pair_4090",
        "experiment2_atypical/pair_4092",
        "experiment2_atypical/pair_4093",
        "experiment2_atypical/pair_4094",
        "experiment2_atypical/pair_4097",
        "experiment2_atypical/pair_4098",
        "experiment2_atypical/pair_4099",
    }
)

PHASE1_SELECTED_PROJECTION_FINGERPRINT = "1289888132922e697c6083d164721da87086db0bfcdf85bea68e597820b4b05c"
DELINS_TRUTH_KEYS = frozenset(
    {
        "experiment2_atypical/pair_4020",
        "experiment2_atypical/pair_4021",
        "experiment2_atypical/pair_4022",
        "experiment2_atypical/pair_4023",
        "experiment2_atypical/pair_4024",
        "experiment2_atypical/pair_4025",
        "experiment2_atypical/pair_4026",
        "experiment2_atypical/pair_4027",
        "experiment2_atypical/pair_4028",
        "experiment2_atypical/pair_4029",
    }
)
DELINS_SELECTED_REPRESENTATION_KEYS = frozenset(
    {
        "experiment2_atypical/pair_4020",
        "experiment2_atypical/pair_4022",
        "experiment2_atypical/pair_4023",
        "experiment2_atypical/pair_4024",
        "experiment2_atypical/pair_4026",
        "experiment2_atypical/pair_4028",
        "experiment2_atypical/pair_4029",
    }
)
DELINS_NO_KESTREL_RESULT_KEYS = frozenset(
    {
        "experiment2_atypical/pair_4021",
        "experiment2_atypical/pair_4025",
        "experiment2_atypical/pair_4027",
    }
)
PAIR_4092_KEY = "experiment2_atypical/pair_4092"
PRA_CURRENT_TIER_FINGERPRINT = "a5271e8eed5a69cc391e4f95dc276c6827b8f25d47601890b54d99c2035939ef"
PRA_IDENTITY_OUTCOME_FINGERPRINT = "0b9b54eb0cf48a03b302a56c5170344228a471deaee628e1c48eac16a51cc9a3"
PRA_PUBLIC_TRUTH_IDENTITY_FINGERPRINT = "bd8f8b56cf8d5516affcff6116b6509213ecfdd4e32d5bdaa7efed0c35369ce9"
PRA_BAM_ELIGIBLE_FINGERPRINT = "0baf09e74569ed1e51710c8a2f844baab5ed015a96e32751a3f0a14d5f53b4b5"
PRA_BAM_FETCH_FINGERPRINT = "ec9f859ab3c2ad78399ee373a708a15d8aa8b4c09b8a0f134a6f2420947fa72f"
PRA_COMPLETE_BAM_FINGERPRINT = "a90dcdf5fa13873f29f272d30b7c7de3fe9b4ac39f9f0dad8b29e1dbf5d12839"
PRA_PUBLIC_IDENTITY_QUARTET_FINGERPRINT = "960fe7316a7f5f1d4de684f6f5c4ef5510a965e63face10d504b4a2ac9d33e32"
PRB_POLICY_STABLE_PUBLIC_FINGERPRINT = "0dd25343376d388aafd772bd0ec4a3a0653f0f2dbecbda2fdc1589ca95c3d048"
PRB_UNAFFECTED_PUBLIC_FINGERPRINT = "abd3404291407cbd0d95cf534f75443c79ae6721be417190e233d2cfa1d2da5f"
PRB_EXPECTED_EXPLANATION_FLAGS = {
    "experiment2_atypical/pair_4010/mutated": "known-variant",
    "experiment2_atypical/pair_4012/mutated": (
        "known-variant;low-haplotype-record-support;thin-haplotype-record-support"
    ),
    "experiment2_atypical/pair_4014/mutated": (
        "known-variant;low-haplotype-record-support;motif-context-diverges;"
        "position-ambiguous;thin-haplotype-record-support"
    ),
    "experiment2_atypical/pair_4015/mutated": (
        "known-variant;low-haplotype-record-support;motif-context-diverges;thin-haplotype-record-support"
    ),
    "experiment2_atypical/pair_4017/mutated": (
        "known-variant;low-haplotype-record-support;thin-haplotype-record-support"
    ),
    "experiment2_atypical/pair_4083/mutated": "known-variant;position-ambiguous",
    "experiment2_atypical/pair_4087/mutated": (
        "low-haplotype-record-support;motif-context-diverges;position-ambiguous;"
        "representation-of-caller-call;spans-unit-junction;thin-haplotype-record-support"
    ),
    "experiment2_atypical/pair_4088/mutated": (
        "low-haplotype-record-support;motif-context-diverges;position-ambiguous;"
        "representation-of-caller-call;spans-unit-junction;thin-haplotype-record-support"
    ),
}
PRA_BAM_TRUTH_MATCH_KEYS = frozenset(
    {
        "experiment2_atypical/pair_4033",
        "experiment2_atypical/pair_4035",
        "experiment2_atypical/pair_4038",
        "experiment2_atypical/pair_4039",
    }
)
PRA_BAM_POSITIVE_CONTROL_KEYS = frozenset(
    {
        "experiment2_atypical/pair_4033",
        "experiment2_atypical/pair_4035",
    }
)


@pytest.fixture(scope="module")
def corpus() -> GoldenCorpus:
    """Load both owner-supplied roots without a skip fallback."""
    return load_golden_corpus(SIM_ROOT, ADVNTR_ROOT)


@pytest.fixture(scope="module")
def uut_replay(corpus: GoldenCorpus, tmp_path_factory: pytest.TempPathFactory) -> UutReplay:
    """Run real reconciliation on isolated TSV copies with original BAM inputs."""
    output_root = tmp_path_factory.mktemp("molecular_identity_replay")
    eligible = 0
    bam_fetches = 0
    eligible_keys: set[str] = set()
    fetches_by_key: dict[str, int] = {}
    source_bam_stats: dict[Path, tuple[int, int]] = {}
    bam_input_symlinks: list[Path] = []
    active_key = ""
    replay_phase = ""
    run_configuration = resolve_run_configuration()
    kestrel_component = run_configuration.kestrel
    artifact_flags = kestrel_component["artifact_flags"]
    selection = kestrel_component["selection"]
    assert isinstance(artifact_flags, tuple)
    assert isinstance(selection, Mapping)
    original_is_candidate = nomenclature_bam_adapter.is_candidate
    original_rescue = BamRescuer.rescue_with_identity_evidence
    identity_component = translation_component_from_config(run_configuration.nomenclature)

    def counted_is_candidate(call: Nomenclature) -> bool:
        nonlocal eligible
        result = original_is_candidate(call)
        if replay_phase == "eligibility-probe" and result:
            eligible += 1
            eligible_keys.add(active_key)
        return result

    def counted_rescue(
        self: BamRescuer,
        contig: str,
        position: int,
        locus: Any,
        component: Any,
    ) -> tuple[BamConsensus | None, Any]:
        nonlocal bam_fetches
        if replay_phase == "reconcile":
            bam_fetches += 1
            fetches_by_key[active_key] = fetches_by_key.get(active_key, 0) + 1
        return original_rescue(self, contig, position, locus, component)

    def rebuild_kestrel_result(source_dir: Path, output_dir: Path) -> None:
        """Replay current capture/selection/publication from the frozen pre-result."""
        source_result = source_dir / "kestrel_result.tsv"
        raw_pre_result = pd.read_csv(
            source_dir / "kestrel_pre_result.tsv",
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )
        pre_result = pd.read_csv(
            source_dir / "kestrel_pre_result.tsv",
            sep="\t",
            dtype={"Motif_fasta": str, "POS_fasta": str, "Motif": str},
            keep_default_na=False,
        )
        pre_result = add_artifact_gate(pre_result, artifact_flags)
        records = pre_result.to_dict("records")
        for record in records:
            record["Motif_sequence"] = str(record["Motif_sequence"])
        candidates = capture_kestrel_observations(records, identity_component)
        capture_rows = [candidate_capture_cells(candidate) for candidate in candidates.candidates]
        diagnostic_rows = [
            identity_translation_diagnostic_cells(candidate.observation.translation)
            for candidate in candidates.candidates
        ]
        for column in IDENTITY_CAPTURE_COLUMNS:
            pre_result[column] = [cells[column] for cells in capture_rows]
        for column in IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS:
            pre_result[column] = [cells[column] for cells in diagnostic_rows]
        evidenced = with_candidate_evidence(candidates, pre_result.to_dict("records"))
        passing_mask = pre_result[list(FILTER_COLUMNS)].all(axis=1)
        passing_ordinals = tuple(int(value) for value in pre_result.loc[passing_mask, IDENTITY_CAPTURE_COLUMNS[5]])
        selected = filter_final_dataframe(pre_result, str(output_dir), selection=selection)

        if selected.empty:
            shutil.copyfile(source_result, output_dir / "kestrel_result.tsv")
            return

        selected = selected.drop(columns=list(IDENTITY_TRANSLATION_DIAGNOSTIC_COLUMNS))
        selected_ordinal = int(selected.iloc[0][IDENTITY_CAPTURE_COLUMNS[5]])
        selected_candidates = overlay_legacy_projection(evidenced, passing_ordinals, selected_ordinal)
        for column, value in selected_candidate_cells(selected_candidates).items():
            selected[column] = value
        for column in SELECTED_PROJECTION_COLUMNS:
            selected[column] = raw_pre_result.iloc[selected_ordinal][column]
        annotated = nomenclature_annotate.annotate_kestrel_frame(
            selected,
            output_dir,
            identity_component=identity_component,
        )
        header = [line for line in source_result.read_text(encoding="utf-8").splitlines() if line.startswith("##")]
        with (output_dir / "kestrel_result.tsv").open("w", encoding="utf-8") as handle:
            if header:
                handle.write("\n".join(header) + "\n")
            annotated.to_csv(handle, sep="\t", index=False)

    def rebuild_advntr_result(source_dir: Path, output_dir: Path) -> Path:
        """Replay the current adVNTR processor from its frozen raw caller VCF."""
        raw = output_dir / "input_adVNTR.vcf"
        shutil.copyfile(source_dir / "output_adVNTR.vcf", raw)
        process_advntr_output(str(raw), str(output_dir), "output", None)
        return output_dir / "output_adVNTR_result.tsv"

    def probe_candidate(kestrel_path: Path, advntr_path: Path) -> None:
        """Observe unchanged BAM eligibility before later production early returns."""
        nonlocal replay_phase
        kestrel = pd.read_csv(kestrel_path, sep="\t", comment="#", dtype=str)
        advntr = pd.read_csv(advntr_path, sep="\t", dtype=str)
        supports: dict[str, int | None] = {}
        kestrel_rows = [row for _, row in kestrel.iterrows() if not nomenclature_annotate._is_negative(row)]
        vcf_calls = []
        for row in kestrel_rows:
            call = nomenclature_annotate._kestrel_row_call(row)
            if call is None:
                continue
            vcf_calls.append(call)
            depth = nomenclature_annotate._as_int(row.get("Estimated_Depth_AlternateVariant"))
            existing = supports.get("kestrel_vcf")
            supports["kestrel_vcf"] = (
                depth if "kestrel_vcf" not in supports else nomenclature_annotate._lesser(existing, depth)
            )
        advntr_keep = [not nomenclature_annotate._is_negative(row) for _, row in advntr.iterrows()]
        advntr_calls_by_row = nomenclature_annotate._advntr_calls_by_row(advntr, advntr_keep, supports)
        advntr_calls = [call for calls in advntr_calls_by_row for call in calls]
        replay_phase = "eligibility-probe"
        try:
            nomenclature_bam_adapter.is_candidate(
                nomenclature_annotate.reconcile(*vcf_calls, *advntr_calls, supports=supports)
            )
        finally:
            replay_phase = ""

    def copy_and_reconcile(key: str, condition: str) -> None:
        nonlocal active_key, replay_phase
        experiment, pair_id = key.split("/", maxsplit=1)
        source_kestrel = SIM_ROOT / experiment / "vntyper" / pair_id / condition / "kestrel"
        source_advntr = ADVNTR_ROOT / experiment / pair_id / condition / "advntr"
        output_sample = output_root / experiment / pair_id / condition
        output_kestrel = output_sample / "kestrel"
        output_advntr = output_sample / "advntr"
        output_kestrel.mkdir(parents=True)
        output_advntr.mkdir(parents=True)
        for source_name in ("output.bam", "output.bam.bai"):
            source = source_kestrel / source_name
            if source.is_file():
                resolved_source = source.resolve()
                stat = resolved_source.stat()
                source_bam_stats[resolved_source] = (stat.st_size, stat.st_mtime_ns)
                link = output_kestrel / source_name
                link.symlink_to(resolved_source)
                bam_input_symlinks.append(link)

        active_key = key
        replay_phase = "kestrel"
        rebuild_kestrel_result(source_kestrel, output_kestrel)
        advntr_output = rebuild_advntr_result(source_advntr, output_advntr)
        kestrel_output = output_kestrel / "kestrel_result.tsv"
        if condition == "mutated":
            probe_candidate(kestrel_output, advntr_output)
        replay_phase = "reconcile"
        nomenclature_annotate.reconcile_caller_outputs(
            kestrel_output,
            advntr_output,
            output_kestrel,
        )
        replay_phase = ""

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(nomenclature_bam_adapter, "is_candidate", counted_is_candidate)
        patch.setattr(BamRescuer, "rescue_with_identity_evidence", counted_rescue)
        for key in sorted(corpus.mutated_keys):
            copy_and_reconcile(key, "mutated")
        for key in sorted(corpus.control_keys):
            copy_and_reconcile(key, "normal")

    observed = load_golden_corpus(SIM_ROOT, output_root)
    for link in bam_input_symlinks:
        resolved_source = link.resolve(strict=True)
        stat = resolved_source.stat()
        assert link.is_symlink()
        assert not resolved_source.is_relative_to(output_root)
        assert (stat.st_size, stat.st_mtime_ns) == source_bam_stats[resolved_source]
    return UutReplay(
        corpus=observed,
        eligible=eligible,
        bam_fetches=bam_fetches,
        mutated_replays=len(corpus.mutated_keys),
        control_replays=len(corpus.control_keys),
        eligible_keys=frozenset(eligible_keys),
        fetches_by_key=fetches_by_key,
        bam_input_symlinks=len(bam_input_symlinks),
    )


def test_oracle_import_closure_is_independent() -> None:
    """A production helper entering any recursive local import must fail the guard."""
    oracle_path = Path(identity_oracle.__file__)
    scanned = assert_independent_import_closure(oracle_path, Path(__file__).parents[2])
    assert oracle_path.resolve() in scanned


def test_policy_projection_excludes_only_the_internal_reconciled_identity_handoff() -> None:
    """Private identity persistence is not public, while an arbitrary new field remains locked."""

    def projections(
        row: dict[str, str],
    ) -> tuple[
        dict[str, tuple[tuple[tuple[str, str], ...], ...]],
        dict[str, tuple[tuple[tuple[str, str], ...], ...]],
    ]:
        stable: dict[str, tuple[tuple[tuple[str, str], ...], ...]] = {}
        unaffected: dict[str, tuple[tuple[tuple[str, str], ...], ...]] = {}
        explanations: dict[str, tuple[tuple[str, str, str, str], ...]] = {}
        identity_oracle._record_policy_projections(  # noqa: SLF001 - directly tests the oracle boundary
            stable,
            unaffected,
            explanations,
            "experiment/sample",
            "mutated",
            {"kestrel": (row,), "advntr": ()},
        )
        return stable, unaffected

    baseline = projections({"Nomenclature": "59dupC"})
    with_internal_handoff = projections(
        {
            "Nomenclature": "59dupC",
            "__Reconciled_Molecular_Identity": "MUC1-X-60-coding-v1|60|59|-|C",
        }
    )
    with_new_public_field = projections({"Nomenclature": "59dupC", "New_Public_Field": "changed"})

    assert with_internal_handoff == baseline
    assert identity_oracle.public_projection_fingerprint(with_new_public_field[0]) != (
        identity_oracle.public_projection_fingerprint(baseline[0])
    )
    assert identity_oracle.public_projection_fingerprint(with_new_public_field[1]) != (
        identity_oracle.public_projection_fingerprint(baseline[1])
    )


def test_recursive_guard_rejects_nested_production_import(tmp_path: Path) -> None:
    """The guard must recurse; scanning only the entry file would miss this import."""
    entry = tmp_path / "entry.py"
    helper = tmp_path / "helper.py"
    entry.write_text("from helper import value\n", encoding="utf-8")
    helper.write_text("from vntyper.scripts.nomenclature import reconcile\nvalue = reconcile\n", encoding="utf-8")

    with pytest.raises(AssertionError, match=r"helper\.py.*vntyper\.scripts\.nomenclature"):
        assert_independent_import_closure(entry, tmp_path)


@pytest.mark.parametrize(
    "source",
    [
        "import importlib\nvalue = importlib.import_module('vntyper.scripts.nomenclature')\n",
        "from importlib import import_module\nvalue = import_module('vntyper.scripts.nomenclature')\n",
        "import importlib as il\nvalue = il.import_module('vntyper.scripts.nomenclature')\n",
        "value = __import__('vntyper.scripts.nomenclature')\n",
    ],
)
def test_recursive_guard_rejects_literal_dynamic_production_import(tmp_path: Path, source: str) -> None:
    """Literal dynamic imports must not bypass the independent-oracle boundary."""
    entry = tmp_path / "entry.py"
    entry.write_text(source, encoding="utf-8")

    with pytest.raises(AssertionError, match=r"entry\.py.*vntyper\.scripts\.nomenclature"):
        assert_independent_import_closure(entry, tmp_path)


def test_truth_oracle_derives_every_literal_identity_and_name(corpus) -> None:
    """A coordinate, allele, normalization, or display projection change must fail."""
    assert corpus.expected_by_class == EXPECTED_CLASS_IDENTITIES
    assert len(corpus.expected_by_sample) == 200
    assert {sample.expected for sample in corpus.expected_by_sample.values() if sample.mutation == "dupC"} == {
        ExpectedIdentity("MUC1-X-60-coding-v1|60|59|-|C", "59dupC")
    }


def test_oracle_rejects_wrong_identity_under_a_truth_exact_caller_name() -> None:
    """One repeated nonempty identity cannot satisfy independently exact caller rows."""
    expected = EXPECTED_CLASS_IDENTITIES["delGCCCA"]
    public_rows = {
        "kestrel": (
            {
                "Nomenclature_Kestrel": expected.name,
                "Molecular_Identity": EXPECTED_CLASS_IDENTITIES["dupC"].identity,
                "Molecular_Identity_Status": "unique",
                "Equivalent_Representation_Count": "1",
                "Identity_Hypothesis_Count": "1",
            },
        ),
        "advntr": (),
    }

    violations, counts = identity_oracle._identity_observations("sample", "mutated", public_rows, expected)

    assert counts == identity_oracle.IdentityCounts(resolved=1, exact=0, wrong=1)
    assert len(violations) == 1
    assert f"expected {expected.identity!r}" in violations[0]


@pytest.mark.parametrize(
    ("status", "hypotheses"),
    [
        ("unique", "2"),
        ("legacy-selected-among-multiple", "1"),
    ],
)
def test_oracle_rejects_statuses_inconsistent_with_hypothesis_count(status: str, hypotheses: str) -> None:
    """Resolved status must encode whether the caller had one or multiple identities."""
    identity = EXPECTED_CLASS_IDENTITIES["dupC"].identity
    public_rows = {
        "kestrel": (
            {
                "Molecular_Identity": identity,
                "Molecular_Identity_Status": status,
                "Equivalent_Representation_Count": "1",
                "Identity_Hypothesis_Count": hypotheses,
            },
        ),
        "advntr": (),
    }

    violations, _ = identity_oracle._identity_observations("sample", "mutated", public_rows, None)

    assert len(violations) == 1
    assert "hypothesis count" in violations[0]


def test_oracle_rejects_inconsistent_caller_wide_and_advntr_quartet_counts() -> None:
    """All caller rows share hypotheses and adVNTR equivalence is exactly observable."""
    identity = EXPECTED_CLASS_IDENTITIES["dupC"].identity
    rows = (
        {
            "Molecular_Identity": identity,
            "Molecular_Identity_Status": "unique",
            "Equivalent_Representation_Count": "1",
            "Identity_Hypothesis_Count": "1",
        },
        {
            "Molecular_Identity": identity,
            "Molecular_Identity_Status": "legacy-selected-among-multiple",
            "Equivalent_Representation_Count": "1",
            "Identity_Hypothesis_Count": "2",
        },
    )

    violations, _ = identity_oracle._identity_observations(
        "sample",
        "mutated",
        {"kestrel": (), "advntr": rows},
        None,
    )

    assert any("caller-wide hypothesis count" in violation for violation in violations)
    assert any("equivalent representation count" in violation for violation in violations)


def _valid_bam_replay_payload() -> dict[str, Any]:
    dupc = EXPECTED_CLASS_IDENTITIES["dupC"].identity
    dupa = EXPECTED_CLASS_IDENTITIES["dupA"].identity
    return {
        "schema_version": "bam-identity-replay-v1",
        "loci": [
            {
                "candidate_observation_ordinals": [2, 7],
                "state": "observed",
                "evidence": {
                    "eligible_record_count": 3,
                    "records": [
                        {
                            "identities": [dupc],
                            "candidate_observation_ordinals": [[2]],
                            "minimum_kmer_depth": 5,
                        },
                        {
                            "identities": [dupc, dupa],
                            "candidate_observation_ordinals": [[2], [7]],
                            "minimum_kmer_depth": None,
                        },
                        {
                            "identities": [dupc],
                            "candidate_observation_ordinals": [[2]],
                            "minimum_kmer_depth": 0,
                        },
                    ],
                    "counts": [
                        {"identity": dupc, "record_count": 3},
                        {"identity": dupa, "record_count": 1},
                    ],
                },
            },
            {"candidate_observation_ordinals": [11], "state": "unavailable", "evidence": None},
        ],
    }


def _write_bam_replay(tmp_path: Path, payload: object) -> Path:
    sample = tmp_path / "kestrel"
    sample.mkdir()
    (sample / "bam_identity_replay.v1.json").write_text(json.dumps(payload), encoding="utf-8")
    return sample


def test_bam_oracle_derives_counts_and_rejects_a_fabricated_serialized_winner(tmp_path: Path) -> None:
    """Serialized count claims cannot override independently counted record identities."""
    payload = _valid_bam_replay_payload()
    evidence = payload["loci"][0]["evidence"]
    evidence["counts"][0]["record_count"] = 99
    sample = _write_bam_replay(tmp_path, payload)

    with pytest.raises(AssertionError, match="derived record counts"):
        identity_oracle._bam_replay_observation(sample, EXPECTED_CLASS_IDENTITIES["dupC"].identity)


@pytest.mark.parametrize(
    "defect",
    [
        "unknown-root-key",
        "unknown-locus-key",
        "unknown-evidence-key",
        "unknown-record-key",
        "unknown-count-key",
        "overlapping-group",
        "invalid-state-pairing",
        "eligible-record-mismatch",
        "identity-binding-cardinality",
        "out-of-group-binding",
        "noncanonical-binding-order",
        "invalid-identity",
        "duplicate-record-identity",
        "noncanonical-count-order",
        "zero-count",
    ],
)
def test_bam_oracle_rejects_closed_schema_binding_and_identity_defects(tmp_path: Path, defect: str) -> None:
    """Unknown keys, malformed bindings, and noncanonical identities fail closed."""
    payload = copy.deepcopy(_valid_bam_replay_payload())
    if defect == "unknown-root-key":
        payload["extra"] = True
    elif defect == "unknown-locus-key":
        payload["loci"][0]["extra"] = True
    elif defect == "unknown-evidence-key":
        payload["loci"][0]["evidence"]["extra"] = True
    elif defect == "unknown-record-key":
        payload["loci"][0]["evidence"]["records"][0]["extra"] = True
    elif defect == "unknown-count-key":
        payload["loci"][0]["evidence"]["counts"][0]["extra"] = True
    elif defect == "overlapping-group":
        payload["loci"][1]["candidate_observation_ordinals"] = [7]
    elif defect == "invalid-state-pairing":
        payload["loci"][0]["state"] = "unavailable"
    elif defect == "eligible-record-mismatch":
        payload["loci"][0]["evidence"]["eligible_record_count"] = 2
    elif defect == "identity-binding-cardinality":
        payload["loci"][0]["evidence"]["records"][0]["candidate_observation_ordinals"] = []
    elif defect == "out-of-group-binding":
        payload["loci"][0]["evidence"]["records"][0]["candidate_observation_ordinals"] = [[99]]
    elif defect == "noncanonical-binding-order":
        payload["loci"][0]["evidence"]["records"][1]["candidate_observation_ordinals"] = [[7], [2]]
    elif defect == "invalid-identity":
        payload["loci"][0]["evidence"]["records"][0]["identities"] = ["59dupC"]
    elif defect == "duplicate-record-identity":
        record = payload["loci"][0]["evidence"]["records"][1]
        record["identities"] = [record["identities"][0], record["identities"][0]]
    elif defect == "noncanonical-count-order":
        payload["loci"][0]["evidence"]["counts"].reverse()
    else:
        payload["loci"][0]["evidence"]["counts"][0]["record_count"] = 0
    sample = _write_bam_replay(tmp_path, payload)

    with pytest.raises(AssertionError):
        identity_oracle._bam_replay_observation(sample, EXPECTED_CLASS_IDENTITIES["dupC"].identity)


def test_bam_oracle_rejects_duplicate_json_keys(tmp_path: Path) -> None:
    """Duplicate object keys cannot be hidden by the standard decoder's last-key rule."""
    payload = json.dumps(_valid_bam_replay_payload())
    payload = payload.replace(
        '"schema_version": "bam-identity-replay-v1",',
        '"schema_version": "bam-identity-replay-v1", "schema_version": "bam-identity-replay-v1",',
        1,
    )
    sample = tmp_path / "kestrel"
    sample.mkdir()
    (sample / "bam_identity_replay.v1.json").write_text(payload, encoding="utf-8")

    with pytest.raises(AssertionError, match="duplicate"):
        identity_oracle._bam_replay_observation(sample, EXPECTED_CLASS_IDENTITIES["dupC"].identity)


def test_identity_quartet_fingerprint_changes_for_one_cell() -> None:
    """One public quartet-cell drift must change its literal compatibility digest."""
    projection = {"sample/kestrel": (("identity", "unique", "1", "1"),)}
    changed = {"sample/kestrel": (("identity", "unique", "1", "2"),)}

    assert identity_oracle.selected_projection_fingerprint(projection) == (
        "55bec88299f6d83eb5b3e36e9bf970d6d8a077201a8750f8a44c3f69c4458238"
    )
    assert identity_oracle.selected_projection_fingerprint(changed) != identity_oracle.selected_projection_fingerprint(
        projection
    )


def test_complete_public_identity_quartet_projection_is_literal(uut_replay: UutReplay) -> None:
    """Every public identity/status/representation/hypothesis cell is compatibility-locked."""
    projection = getattr(uut_replay.corpus, "identity_projection_by_sample", {})
    assert len(projection) == 800
    assert sum(len(rows) for rows in projection.values()) == 374
    assert all(len(cells) == len(IDENTITY_COLUMNS) for rows in projection.values() for cells in rows)
    assert identity_oracle.selected_projection_fingerprint(projection) == PRA_PUBLIC_IDENTITY_QUARTET_FINGERPRINT


def test_both_roots_and_every_sample_are_loaded(corpus: GoldenCorpus, uut_replay: UutReplay) -> None:
    """A missing root, sample, condition, or caller artifact must fail rather than skip."""
    assert corpus.sim_root == SIM_ROOT.resolve()
    assert corpus.advntr_root == ADVNTR_ROOT.resolve()
    assert corpus.mutated_samples == 200
    assert corpus.control_samples == 200
    assert uut_replay.corpus.sim_root == SIM_ROOT.resolve()
    assert uut_replay.mutated_replays == 200
    assert uut_replay.control_replays == 200


def test_pr_a_recurrent_state_collision_baseline_is_literal(corpus: GoldenCorpus) -> None:
    """Freeze every B1 collision before governed disposition changes any decision."""
    assert corpus.recurrent_state_collisions == (
        ("experiment2_atypical/pair_4010/mutated", 0, "I23_6_G_LEN1"),
        ("experiment2_atypical/pair_4012/mutated", 0, "I21_2_T_LEN1"),
        ("experiment2_atypical/pair_4012/mutated", 1, "I23_6_G_LEN1"),
        ("experiment2_atypical/pair_4014/mutated", 0, "I21_2_T_LEN1"),
        ("experiment2_atypical/pair_4015/mutated", 0, "I21_2_T_LEN1"),
        ("experiment2_atypical/pair_4015/mutated", 1, "I23_6_G_LEN1"),
        ("experiment2_atypical/pair_4017/mutated", 0, "I21_2_T_LEN1"),
        ("experiment2_atypical/pair_4017/mutated", 1, "I23_6_G_LEN1"),
        ("experiment2_atypical/pair_4083/mutated", 0, "D17_2&D18_2&D19_2&D20_2&D21_2"),
        ("experiment2_atypical/pair_4087/mutated", 0, "D17_2&D18_2&D19_2&D20_2&D21_2"),
        ("experiment2_atypical/pair_4088/mutated", 0, "D17_2&D18_2&D19_2&D20_2&D21_2"),
    )
    assert len(corpus.recurrent_state_collisions) == 11
    assert len({key for key, _, _ in corpus.recurrent_state_collisions}) == 8


def test_pr_b_changes_only_governed_explanations_while_every_public_finding_stays_fixed(
    corpus: GoldenCorpus,
    uut_replay: UutReplay,
) -> None:
    """Compare against the measured B1 replay, never expectations derived from B2."""
    observed = uut_replay.corpus

    assert identity_oracle.public_projection_fingerprint(observed.policy_stable_projection_by_sample) == (
        PRB_POLICY_STABLE_PUBLIC_FINGERPRINT
    )
    assert identity_oracle.public_projection_fingerprint(observed.unaffected_public_projection_by_sample) == (
        PRB_UNAFFECTED_PUBLIC_FINGERPRINT
    )

    collision_rows = 0
    admissible_rows = 0
    affected_decisions: set[str] = set()
    for projection_key, after_rows in observed.policy_explanation_by_sample.items():
        decision_key, caller = projection_key.rsplit("/", maxsplit=1)
        expected_flags = PRB_EXPECTED_EXPLANATION_FLAGS.get(decision_key)

        for variant, flags, note, disposition in after_rows:
            governed = caller == "advntr" and variant in identity_oracle.RECURRENT_STATE_EVIDENCE
            affected = governed or (caller == "kestrel" and expected_flags is not None)
            if not affected:
                if caller == "advntr":
                    assert disposition == "admissible"
                    admissible_rows += 1
                else:
                    assert disposition == "<missing>"
                continue

            affected_decisions.add(decision_key)
            if governed:
                collision_rows += 1
                assert disposition == "identity-insufficient"
            else:
                assert disposition == "<missing>"
            assert flags == expected_flags
            assert note.endswith(f"; {identity_oracle.GOVERNED_ASSERTION}")

    assert collision_rows == 11
    assert admissible_rows == 185
    assert affected_decisions == frozenset(PRB_EXPECTED_EXPLANATION_FLAGS)


def test_historical_phase1_projection_is_literal_per_sample(corpus: GoldenCorpus) -> None:
    """Aggregate swaps and disappearance of a lower tier must fail on exact sample keys."""
    observed = corpus
    assert observed.tier_keys == {
        "A": PHASE1_TIER_A_KEYS,
        "B": PHASE1_TIER_B_KEYS,
        "C": frozenset(),
    }
    assert observed.total.displayed == 154
    assert observed.total.exact == 136
    assert observed.total.wrong == 18
    assert observed.by_tier["A"].displayed == 53
    assert observed.by_tier["A"].exact == 53
    assert observed.by_tier["A"].wrong == 0
    assert observed.by_tier["B"].displayed == 101
    assert observed.by_tier["B"].exact == 83
    assert observed.by_tier["B"].wrong == 18
    assert observed.by_tier["C"].displayed == 0
    assert observed.by_tier["C"].exact == 0
    assert observed.by_tier["C"].wrong == 0
    assert observed.control_findings == 0


def test_historical_bam_consultation_and_canonical_dupc_are_literal(
    corpus: GoldenCorpus,
    uut_replay: UutReplay,
) -> None:
    """Eligibility, fetch cardinality, or canonical dupC naming drift must fail."""
    assert uut_replay.eligible == 83
    assert uut_replay.bam_fetches == 68
    assert corpus.dupc_vcf_names == ("59dupC",) * 96


def test_selected_kestrel_projection_is_exact_and_ordered(corpus: GoldenCorpus, uut_replay: UutReplay) -> None:
    """Every selected raw row, empty result, and deterministic tie winner is frozen."""
    historical = corpus.selected_projection_by_sample
    current = uut_replay.corpus.selected_projection_by_sample
    assert identity_oracle.selected_projection_fingerprint(historical) == PHASE1_SELECTED_PROJECTION_FINGERPRINT
    assert len(historical) == 400
    assert sum(bool(rows) for rows in historical.values()) == 178
    mismatches = {
        key: {"expected": historical[key], "observed": current.get(key)}
        for key in historical
        if current.get(key) != historical[key]
    }
    assert mismatches == {}


def test_structural_delins_and_pair_4092_are_never_fabricated(
    corpus: GoldenCorpus,
    uut_replay: UutReplay,
) -> None:
    """Absent complete candidate representations cannot be synthesized from fragments."""
    truth_delins = frozenset(
        key for key, expectation in corpus.expected_by_sample.items() if expectation.mutation == "delinsAT"
    )
    assert truth_delins == DELINS_TRUTH_KEYS
    assert DELINS_SELECTED_REPRESENTATION_KEYS | DELINS_NO_KESTREL_RESULT_KEYS == DELINS_TRUTH_KEYS
    historical_projection = corpus.selected_projection_by_sample
    assert all(historical_projection[f"{key}/mutated"] for key in DELINS_SELECTED_REPRESENTATION_KEYS)
    assert all(not historical_projection[f"{key}/mutated"] for key in DELINS_NO_KESTREL_RESULT_KEYS)

    observed = uut_replay.corpus
    protected = DELINS_SELECTED_REPRESENTATION_KEYS | {PAIR_4092_KEY}
    fabricated = (observed.public_truth_identity_keys & protected) - observed.bam_truth_match_keys
    assert fabricated == frozenset()
    assert not (observed.public_truth_identity_keys & DELINS_TRUTH_KEYS)
    assert PAIR_4092_KEY not in observed.public_truth_identity_keys


def test_positive_caller_rows_publish_the_complete_identity_quartet(uut_replay) -> None:
    """A positive caller row missing any identity field must fail semantically."""
    observed = uut_replay.corpus
    assert observed.identity_contract_violations == (), (
        f"{len(observed.identity_contract_violations)} positive caller rows violate the identity quartet "
        f"{IDENTITY_COLUMNS}: {observed.identity_contract_violations[:5]}"
    )
    identity_counts = observed.identity_on_truth_exact_names
    assert identity_counts.wrong == 0
    assert identity_counts.exact == identity_counts.resolved


def test_pr_a_measured_identity_policy_baseline_is_literal(
    corpus: GoldenCorpus,
    uut_replay: UutReplay,
) -> None:
    """Pin policy deltas separately from the frozen Phase-1 caller projection."""
    observed = uut_replay.corpus
    expected_total = DisplayCounts(displayed=154, exact=136, wrong=18)
    assert observed.total == corpus.total == expected_total
    assert observed.by_tier == {
        "A": DisplayCounts(0, 0, 0),
        "B": expected_total,
        "C": DisplayCounts(0, 0, 0),
    }
    assert {tier: len(keys) for tier, keys in observed.tier_keys.items()} == {"A": 0, "B": 154, "C": 0}
    assert identity_oracle.sample_sets_fingerprint(observed.tier_keys) == PRA_CURRENT_TIER_FINGERPRINT
    assert frozenset().union(*observed.tier_keys.values()) == frozenset().union(*corpus.tier_keys.values())
    assert observed.control_findings == corpus.control_findings == 0

    assert observed.identity_on_truth_exact_names == IdentityCounts(resolved=241, exact=241, wrong=0)
    assert observed.identity_on_truth_exact_names_by_tier == {
        "A": IdentityCounts(0, 0, 0),
        "B": IdentityCounts(230, 230, 0),
        "C": IdentityCounts(11, 11, 0),
    }
    assert {outcome: len(keys) for outcome, keys in observed.identity_outcome_keys.items()} == {
        "agreement": 98,
        "disagreement": 26,
        "unresolved": 76,
    }
    assert identity_oracle.sample_sets_fingerprint(observed.identity_outcome_keys) == PRA_IDENTITY_OUTCOME_FINGERPRINT
    assert corpus.identity_outcome_keys == {
        "agreement": frozenset(),
        "disagreement": frozenset(),
        "unresolved": corpus.mutated_keys,
    }
    assert len(observed.public_truth_identity_keys) == 148
    assert (
        identity_oracle.sample_sets_fingerprint({"truth": observed.public_truth_identity_keys})
        == PRA_PUBLIC_TRUTH_IDENTITY_FINGERPRINT
    )


def test_bam_positive_controls_pin_eligibility_fetch_evidence_resolution_and_truth(
    uut_replay: UutReplay,
) -> None:
    """Two literal controls prove a fetch is credited only through complete truth evidence."""
    observed = uut_replay.corpus
    fetch_keys = frozenset(uut_replay.fetches_by_key)
    assert len(uut_replay.eligible_keys) == 83
    assert (
        identity_oracle.sample_sets_fingerprint({"eligible": uut_replay.eligible_keys}) == PRA_BAM_ELIGIBLE_FINGERPRINT
    )
    assert len(fetch_keys) == uut_replay.bam_fetches == 68
    assert identity_oracle.sample_sets_fingerprint({"fetch": fetch_keys}) == PRA_BAM_FETCH_FINGERPRINT
    assert len(observed.complete_bam_evidence_keys) == 68
    assert (
        identity_oracle.sample_sets_fingerprint({"complete": observed.complete_bam_evidence_keys})
        == PRA_COMPLETE_BAM_FINGERPRINT
    )
    assert observed.bam_truth_match_keys == PRA_BAM_TRUTH_MATCH_KEYS
    assert uut_replay.bam_input_symlinks == 800

    for key in PRA_BAM_POSITIVE_CONTROL_KEYS:
        assert key in uut_replay.eligible_keys
        assert uut_replay.fetches_by_key[key] == 1
        assert key in observed.complete_bam_evidence_keys
        assert key in observed.bam_truth_match_keys
        assert key in observed.public_truth_identity_keys
        assert key in observed.identity_outcome_keys["unresolved"]
        assert observed.expected_by_sample[key].mutation == "insCCCC"


def test_packaged_profile_provenance_and_pr_b_projection_are_frozen(corpus: GoldenCorpus) -> None:
    """Pin the complete packaged policy independently of its renderer and schema code."""
    assert len(corpus.mutated_keys) == 200
    assert len(corpus.control_keys) == 200

    profile_bytes = (PACKAGE_PROFILES / "decision_profile.json").read_bytes()
    projection_bytes = (PACKAGE_PROFILES / "decision_projection.json").read_bytes()
    sidecar_bytes = (PACKAGE_PROFILES / "decision_profile.sha256").read_bytes()
    assert hashlib.sha256(profile_bytes).hexdigest() == PACKAGED_DECISION_PROFILE_SHA256
    assert hashlib.sha256(projection_bytes).hexdigest() == PR_B_DECISION_PROJECTION_SHA256
    assert sidecar_bytes == f"{PACKAGED_DECISION_PROFILE_SHA256}\n".encode("ascii")

    projection = json.loads(projection_bytes)
    run = resolve_run_configuration()
    resolved = run.decision_profile
    assert (
        resolved.profile_id,
        resolved.profile_revision,
        resolved.profile_kind,
        resolved.source,
        resolved.digest,
    ) == (
        "vntyper-packaged-default",
        "1",
        "packaged",
        "package",
        PACKAGED_DECISION_PROFILE_SHA256,
    )
    assert {
        "advntr": _plain_json(run.advntr),
        "cross_match": _plain_json(run.cross_match),
        "dominance": _plain_json(run.dominance),
        "kestrel": _plain_json(run.kestrel),
        "nomenclature": _plain_json(run.nomenclature),
        "shark": _plain_json(run.shark),
    } == projection
