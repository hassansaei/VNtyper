"""Independent molecular-identity and historical Phase-1 golden contracts."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from tests.golden.identity_oracle import (
    IDENTITY_COLUMNS,
    ExpectedIdentity,
    assert_independent_import_closure,
    load_golden_corpus,
)

pytestmark = pytest.mark.golden

SIM_ROOT = Path(os.environ["VNTYPER_SIM_ROOT"])
ADVNTR_ROOT = Path(os.environ["VNTYPER_ADVNTR_ROOT"])

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

PHASE1_BAM_ELIGIBLE_KEYS = frozenset(
    {
        "experiment1_dupC/pair_3019",
        "experiment1_dupC/pair_3060",
        "experiment2_atypical/pair_4000",
        "experiment2_atypical/pair_4001",
        "experiment2_atypical/pair_4002",
        "experiment2_atypical/pair_4003",
        "experiment2_atypical/pair_4004",
        "experiment2_atypical/pair_4006",
        "experiment2_atypical/pair_4008",
        "experiment2_atypical/pair_4009",
        "experiment2_atypical/pair_4010",
        "experiment2_atypical/pair_4012",
        "experiment2_atypical/pair_4013",
        "experiment2_atypical/pair_4014",
        "experiment2_atypical/pair_4015",
        "experiment2_atypical/pair_4016",
        "experiment2_atypical/pair_4017",
        "experiment2_atypical/pair_4019",
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
        "experiment2_atypical/pair_4030",
        "experiment2_atypical/pair_4031",
        "experiment2_atypical/pair_4032",
        "experiment2_atypical/pair_4033",
        "experiment2_atypical/pair_4034",
        "experiment2_atypical/pair_4035",
        "experiment2_atypical/pair_4036",
        "experiment2_atypical/pair_4037",
        "experiment2_atypical/pair_4038",
        "experiment2_atypical/pair_4039",
        "experiment2_atypical/pair_4040",
        "experiment2_atypical/pair_4051",
        "experiment2_atypical/pair_4052",
        "experiment2_atypical/pair_4053",
        "experiment2_atypical/pair_4054",
        "experiment2_atypical/pair_4055",
        "experiment2_atypical/pair_4056",
        "experiment2_atypical/pair_4057",
        "experiment2_atypical/pair_4058",
        "experiment2_atypical/pair_4059",
        "experiment2_atypical/pair_4060",
        "experiment2_atypical/pair_4061",
        "experiment2_atypical/pair_4062",
        "experiment2_atypical/pair_4063",
        "experiment2_atypical/pair_4064",
        "experiment2_atypical/pair_4065",
        "experiment2_atypical/pair_4066",
        "experiment2_atypical/pair_4067",
        "experiment2_atypical/pair_4068",
        "experiment2_atypical/pair_4069",
        "experiment2_atypical/pair_4070",
        "experiment2_atypical/pair_4071",
        "experiment2_atypical/pair_4072",
        "experiment2_atypical/pair_4073",
        "experiment2_atypical/pair_4074",
        "experiment2_atypical/pair_4075",
        "experiment2_atypical/pair_4077",
        "experiment2_atypical/pair_4078",
        "experiment2_atypical/pair_4079",
        "experiment2_atypical/pair_4081",
        "experiment2_atypical/pair_4084",
        "experiment2_atypical/pair_4085",
        "experiment2_atypical/pair_4086",
        "experiment2_atypical/pair_4087",
        "experiment2_atypical/pair_4088",
        "experiment2_atypical/pair_4089",
        "experiment2_atypical/pair_4090",
        "experiment2_atypical/pair_4091",
        "experiment2_atypical/pair_4093",
        "experiment2_atypical/pair_4094",
        "experiment2_atypical/pair_4095",
        "experiment2_atypical/pair_4096",
        "experiment2_atypical/pair_4097",
        "experiment2_atypical/pair_4098",
        "experiment2_atypical/pair_4099",
    }
)


@pytest.fixture(scope="module")
def corpus():
    """Load both owner-supplied roots without a skip fallback."""
    return load_golden_corpus(SIM_ROOT, ADVNTR_ROOT)


def test_oracle_import_closure_is_independent() -> None:
    """A production helper entering any recursive local import must fail the guard."""
    scanned = assert_independent_import_closure(Path(__file__), Path(__file__).parents[2])
    assert Path(__file__).resolve() in scanned
    assert (Path(__file__).parents[1] / "__init__.py").resolve() in scanned
    assert (Path(__file__).parent / "identity_oracle.py").resolve() in scanned


def test_recursive_guard_rejects_nested_production_import(tmp_path: Path) -> None:
    """The guard must recurse; scanning only the entry file would miss this import."""
    entry = tmp_path / "entry.py"
    helper = tmp_path / "helper.py"
    entry.write_text("from helper import value\n", encoding="utf-8")
    helper.write_text("from vntyper.scripts.nomenclature import reconcile\nvalue = reconcile\n", encoding="utf-8")

    with pytest.raises(AssertionError, match=r"helper\.py.*vntyper\.scripts\.nomenclature"):
        assert_independent_import_closure(entry, tmp_path)


def test_truth_oracle_derives_every_literal_identity_and_name(corpus) -> None:
    """A coordinate, allele, normalization, or display projection change must fail."""
    assert corpus.expected_by_class == EXPECTED_CLASS_IDENTITIES
    assert len(corpus.expected_by_sample) == 200
    assert {sample.expected for sample in corpus.expected_by_sample.values() if sample.mutation == "dupC"} == {
        ExpectedIdentity("MUC1-X-60-coding-v1|60|59|-|C", "59dupC")
    }


def test_both_roots_and_every_sample_are_loaded(corpus) -> None:
    """A missing root, sample, condition, or caller artifact must fail rather than skip."""
    assert corpus.sim_root == SIM_ROOT.resolve()
    assert corpus.advntr_root == ADVNTR_ROOT.resolve()
    assert corpus.mutated_samples == 200
    assert corpus.control_samples == 200


def test_historical_phase1_projection_is_literal_per_sample(corpus) -> None:
    """Aggregate swaps and disappearance of a lower tier must fail on exact sample keys."""
    assert corpus.tier_keys == {
        "A": PHASE1_TIER_A_KEYS,
        "B": PHASE1_TIER_B_KEYS,
        "C": frozenset(),
    }
    assert corpus.total.displayed == 154
    assert corpus.total.exact == 136
    assert corpus.total.wrong == 18
    assert corpus.by_tier["A"].displayed == 53
    assert corpus.by_tier["A"].exact == 53
    assert corpus.by_tier["A"].wrong == 0
    assert corpus.by_tier["B"].displayed == 101
    assert corpus.by_tier["B"].exact == 83
    assert corpus.by_tier["B"].wrong == 18
    assert corpus.by_tier["C"].displayed == 0
    assert corpus.by_tier["C"].exact == 0
    assert corpus.by_tier["C"].wrong == 0
    assert corpus.control_findings == 0


def test_historical_bam_consultation_and_canonical_dupc_are_literal(corpus) -> None:
    """Eligibility, fetch cardinality, or canonical dupC naming drift must fail."""
    assert len(PHASE1_BAM_ELIGIBLE_KEYS) == 83
    assert corpus.mutated_keys >= PHASE1_BAM_ELIGIBLE_KEYS
    bam_fetch_keys = PHASE1_BAM_ELIGIBLE_KEYS & corpus.kestrel_positive_keys
    assert len(bam_fetch_keys) == 68
    assert corpus.dupc_vcf_names == ("59dupC",) * 96


def test_positive_caller_rows_publish_the_complete_identity_quartet(corpus) -> None:
    """A positive caller row missing any identity field must fail semantically."""
    assert corpus.identity_contract_violations == (), (
        f"{len(corpus.identity_contract_violations)} positive caller rows violate the identity quartet "
        f"{IDENTITY_COLUMNS}: {corpus.identity_contract_violations[:5]}"
    )
