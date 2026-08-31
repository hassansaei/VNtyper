"""Independent molecular-identity and historical Phase-1 golden contracts."""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pytest

from tests.golden import identity_oracle
from tests.golden.identity_oracle import (
    IDENTITY_COLUMNS,
    ExpectedIdentity,
    GoldenCorpus,
    assert_independent_import_closure,
    load_golden_corpus,
)
from vntyper.scripts import nomenclature_annotate
from vntyper.scripts.nomenclature import Nomenclature
from vntyper.scripts.nomenclature_bam import BamConsensus, BamRescuer

pytestmark = pytest.mark.golden

SIM_ROOT = Path(os.environ["VNTYPER_SIM_ROOT"])
ADVNTR_ROOT = Path(os.environ["VNTYPER_ADVNTR_ROOT"])


@dataclass(frozen=True)
class UutReplay:
    """Real reconciliation outputs and directly observed BAM decision counts."""

    corpus: GoldenCorpus
    eligible: int
    bam_fetches: int
    mutated_replays: int
    control_replays: int


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
    count_eligibility_outcome = False
    original_is_candidate = nomenclature_annotate.is_candidate
    original_rescue = BamRescuer.rescue

    def counted_is_candidate(call: Nomenclature) -> bool:
        nonlocal eligible
        result = original_is_candidate(call)
        if count_eligibility_outcome:
            eligible += int(result)
        return result

    def counted_rescue(self: BamRescuer, contig: str, position: int) -> BamConsensus | None:
        nonlocal bam_fetches
        bam_fetches += 1
        return original_rescue(self, contig, position)

    def probe_candidate(kestrel_path: Path, advntr_path: Path) -> None:
        """Exercise the production candidate seam even across later I/O early returns."""
        nonlocal count_eligibility_outcome
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
        count_eligibility_outcome = True
        try:
            nomenclature_annotate._haplotype_calls(
                kestrel_rows,
                output_root / "eligibility-probe-without-bam",
                vcf_calls,
                advntr_calls,
                supports,
            )
        finally:
            count_eligibility_outcome = False

    def copy_and_reconcile(key: str, condition: str) -> None:
        experiment, pair_id = key.split("/", maxsplit=1)
        source_kestrel = SIM_ROOT / experiment / "vntyper" / pair_id / condition / "kestrel"
        source_advntr = ADVNTR_ROOT / experiment / pair_id / condition / "advntr"
        output_sample = output_root / experiment / pair_id / condition
        source_and_relative_paths = (
            (source_kestrel / "kestrel_result.tsv", Path("kestrel/kestrel_result.tsv")),
            (source_advntr / "output_adVNTR_result.tsv", Path("advntr/output_adVNTR_result.tsv")),
        )
        for source, relative_path in source_and_relative_paths:
            destination = output_sample / relative_path
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
        kestrel_output = output_sample / source_and_relative_paths[0][1]
        advntr_output = output_sample / source_and_relative_paths[1][1]
        if condition == "mutated":
            probe_candidate(kestrel_output, advntr_output)
        nomenclature_annotate.reconcile_caller_outputs(
            kestrel_output,
            advntr_output,
            source_kestrel,
        )

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(nomenclature_annotate, "is_candidate", counted_is_candidate)
        patch.setattr(BamRescuer, "rescue", counted_rescue)
        for key in sorted(corpus.mutated_keys):
            copy_and_reconcile(key, "mutated")

    for key in sorted(corpus.control_keys):
        copy_and_reconcile(key, "normal")

    observed = load_golden_corpus(SIM_ROOT, output_root)
    return UutReplay(
        corpus=observed,
        eligible=eligible,
        bam_fetches=bam_fetches,
        mutated_replays=len(corpus.mutated_keys),
        control_replays=len(corpus.control_keys),
    )


def test_oracle_import_closure_is_independent() -> None:
    """A production helper entering any recursive local import must fail the guard."""
    oracle_path = Path(identity_oracle.__file__)
    scanned = assert_independent_import_closure(oracle_path, Path(__file__).parents[2])
    assert oracle_path.resolve() in scanned


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


def test_both_roots_and_every_sample_are_loaded(corpus: GoldenCorpus, uut_replay: UutReplay) -> None:
    """A missing root, sample, condition, or caller artifact must fail rather than skip."""
    assert corpus.sim_root == SIM_ROOT.resolve()
    assert corpus.advntr_root == ADVNTR_ROOT.resolve()
    assert corpus.mutated_samples == 200
    assert corpus.control_samples == 200
    assert uut_replay.corpus.sim_root == SIM_ROOT.resolve()
    assert uut_replay.mutated_replays == 200
    assert uut_replay.control_replays == 200


def test_historical_phase1_projection_is_literal_per_sample(uut_replay) -> None:
    """Aggregate swaps and disappearance of a lower tier must fail on exact sample keys."""
    observed = uut_replay.corpus
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


def test_historical_bam_consultation_and_canonical_dupc_are_literal(uut_replay) -> None:
    """Eligibility, fetch cardinality, or canonical dupC naming drift must fail."""
    assert uut_replay.eligible == 83
    assert uut_replay.bam_fetches == 68
    assert uut_replay.corpus.dupc_vcf_names == ("59dupC",) * 96


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
