"""The anonymous research donation record (schema ``vntyper-donation/2``).

Every field is an enum, a bin or a boolean built in the browser from the donor's own
result files. Nothing here can name a run, a file, a sample or a precise time, and any
field outside this list is refused.
"""

from __future__ import annotations

from typing import Literal, get_args

from pydantic import BaseModel, ConfigDict, Field, model_validator

SCHEMA = "vntyper-donation/2"

DepthBin = Literal[
    "<20", "20-40", "40-60", "60-80", "80-100", "100-150", "150-200", "200-300", "300-500", ">=500", "NA"
]
DepthScoreBin = Literal[
    "<0.002",
    "0.002-0.004",
    "0.004-0.0045",
    "0.0045-0.00469",
    "0.00469-0.005",
    "0.005-0.00515",
    "0.00515-0.006",
    "0.006-0.008",
    "0.008-0.01",
    "0.01-0.02",
    "0.02-0.05",
    "0.05-0.1",
    ">=0.1",
]
# Copied verbatim from the frontend's donate.html kit list.
Kit = Literal[
    "Twist Exome 2.0",
    "Twist Comprehensive Exome",
    "Twist Human Core Exome",
    "Agilent SureSelect Human All Exon V8",
    "Agilent SureSelect Human All Exon V7",
    "Agilent SureSelect Human All Exon V6",
    "Illumina TruSeq Exome",
    "Illumina DNA Prep with Enrichment",
    "IDT xGen Exome Research Panel v2",
    "WGS (PCR-free)",
    "WGS (standard)",
    "Other",
]
Truth = Literal["confirmed_positive", "confirmed_negative", "unconfirmed"]

KITS: tuple[str, ...] = get_args(Kit)
TRUTHS: tuple[str, ...] = get_args(Truth)


class Donation(BaseModel):
    """One donated result. Strict: no coercion, no extra fields."""

    model_config = ConfigDict(extra="forbid", strict=True, frozen=True, populate_by_name=True)

    schema_: Literal["vntyper-donation/2"] = Field(alias="schema")
    vntyper_version: str = Field(pattern=r"^\d+\.\d+\.\d+$", max_length=20)
    coverage_qc: Literal["PASS", "REDUCED", "FAIL", "NA"]
    flank_depth_bin: DepthBin
    array_depth_bin: DepthBin
    kestrel_call: bool
    kestrel_flagged: bool
    kestrel_confidence: Literal["High_Precision", "Low_Precision", "none"]
    depth_score_bin: DepthScoreBin | None
    variant_class: Literal["none", "dupC", "other"]
    kit: Kit
    platform_class: Literal["illumina_short", "other_short", "long_read"]
    truth: Truth
    confirmation_class: Literal["orthogonal_short", "long_read", "none"]
    referral: Literal["suspected_adtkd", "unexplained_ckd", "other"] | None = None

    @model_validator(mode="after")
    def _consistent(self) -> Donation:
        no_row = not self.kestrel_call and not self.kestrel_flagged
        if no_row and self.depth_score_bin is not None:
            raise ValueError("depth_score_bin must be null when there is neither a call nor a flagged call")
        if self.kestrel_call and self.depth_score_bin is None:
            raise ValueError("depth_score_bin is required when kestrel_call is true")
        if (self.variant_class == "none") == self.kestrel_call:
            raise ValueError("variant_class must be 'none' exactly when kestrel_call is false")
        if (self.confirmation_class == "none") != (self.truth == "unconfirmed"):
            raise ValueError("confirmation_class must be 'none' exactly when truth is 'unconfirmed'")
        return self
