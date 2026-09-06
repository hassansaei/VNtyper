"""Pydantic models for anonymous research data donation."""

from __future__ import annotations

import re
from typing import Literal

from pydantic import BaseModel, Field, field_validator

CONTROLLED_KITS = [
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

CONTROLLED_PLATFORMS = [
    "Illumina NovaSeq 6000",
    "Illumina NovaSeq X / X Plus",
    "Illumina NextSeq 500 / 550",
    "Illumina NextSeq 1000 / 2000",
    "Illumina HiSeq 2500 / 4000",
    "Illumina MiSeq",
    "Element AVITI",
    "Ultima UG100",
    "PacBio Revio / Sequel IIe",
    "Oxford Nanopore PromethION / MinION",
    "Other",
]

CONTROLLED_CONFIRMATIONS = [
    "Sanger sequencing",
    "ddPCR / digital PCR",
    "Long-read sequencing (PacBio / ONT)",
    "Mass spectrometry",
    "Targeted amplicon NGS",
    "None / Research only",
]

HPO_REGEX = re.compile(r"^HP:\d{7}$")
COARSE_DATE_REGEX = re.compile(r"^\d{4}-(0[1-9]|1[0-2])$")


class DonationMetadata(BaseModel):
    """Metadata supplied alongside a stripped VNtyper result archive."""

    consent: bool = Field(
        ...,
        description="Explicit informed consent under GDPR Art. 9 for anonymous scientific data donation.",
    )
    kit: str = Field(..., description="Target capture kit or library preparation kit.")
    sequencing_platform: str = Field(..., description="Sequencing instrument/platform used.")
    phenotype_hpo: list[str] = Field(
        default_factory=list,
        description="List of Human Phenotype Ontology (HPO) terms (e.g., HP:0000112).",
    )
    positive_call: bool = Field(
        ...,
        description="Whether the sample has a positive finding for MUC1 VNTR mutation.",
    )
    confirmation_method: str | None = Field(
        None,
        description="Method used to confirm the mutation. Required when positive_call is True.",
    )
    depth_counting_policy: str = Field(
        default="vntr_flank_mean_depth",
        description="Depth counting policy used (e.g., vntr_flank_mean_depth).",
    )
    sex: Literal["XX", "XY", "other", "unknown"] | None = Field(
        None,
        description="Chromosomal/genetic sex (coarse category only).",
    )
    collection_month: str | None = Field(
        None,
        description="Coarse date of sample collection in YYYY-MM format. No day or time permitted.",
    )

    @field_validator("consent")
    @classmethod
    def validate_consent(cls, v: bool) -> bool:
        if not v:
            raise ValueError("Explicit consent under GDPR Art. 9 is required to donate data.")
        return v

    @field_validator("collection_month")
    @classmethod
    def validate_coarse_date(cls, v: str | None) -> str | None:
        if v is None or v == "":
            return None
        if not COARSE_DATE_REGEX.match(v):
            raise ValueError("Collection date must be in YYYY-MM format (e.g. 2024-05). Exact days are not permitted.")
        return v

    @field_validator("phenotype_hpo")
    @classmethod
    def validate_hpo_terms(cls, terms: list[str]) -> list[str]:
        cleaned = []
        for t in terms:
            t_str = str(t).strip().upper()
            if not t_str:
                continue
            if not HPO_REGEX.match(t_str):
                raise ValueError(f"Invalid HPO term '{t}'. Must match format HP:XXXXXXX.")
            cleaned.append(t_str)
        return cleaned

    @field_validator("confirmation_method")
    @classmethod
    def validate_confirmation(cls, v: str | None, info) -> str | None:
        positive = info.data.get("positive_call")
        if positive and (not v or not v.strip()):
            raise ValueError("Confirmation method is required for positive donations.")
        return v


class DonationResponse(BaseModel):
    status: str = "accepted"
    donation_id: str
    run_id: str
    tool_version: str
    message: str = "Donation accepted and verified successfully. Thank you for contributing to ADTKD research!"


class KitAggregate(BaseModel):
    total_samples: int
    positive_count: int
    negative_count: int
    negative_power_sufficient: bool = Field(
        ...,
        description="Whether negative arm power is sufficient (n >= 5) to establish kit-specific thresholds.",
    )
    mean_coverage: float | None = None


class DonationAggregatesResponse(BaseModel):
    """Aggregate statistics for data donations with strict minimum cell-size suppression."""

    total_donations: int
    positive_count: int
    negative_count: int
    by_kit: dict[str, KitAggregate]
    by_platform: dict[str, int]
    top_hpo_terms: dict[str, int]
    minimum_cell_size_threshold: int = 5
