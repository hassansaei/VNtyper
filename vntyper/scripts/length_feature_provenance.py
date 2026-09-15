"""Strict immutable provenance for VNTR length feature measurements."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from types import MappingProxyType
from typing import Literal, cast

from vntyper.scripts.canonical_json import canonical_sha256

LENGTH_FEATURE_CONTEXT_SCHEMA_VERSION = "length-feature-measurement-context-v1"
LENGTH_FEATURE_PROVENANCE_SCHEMA_VERSION = "length-feature-provenance-v1"
COUNTING_POLICY_ID = "primary-mapq0-baseq0-overlap-count-v1"
DENOMINATOR_NAMES = ("CORE", "INVARIANT", "ARRAY", "BOTH_FLANKS")
DenominatorEvidenceKind = Literal["unavailable", "read-pair-identity-qc-proxy"]
_SHA256_LENGTH = 64


@dataclass(frozen=True)
class QueryInterval:
    """One zero-based, half-open interval queried for depth."""

    start: int
    end: int


@dataclass(frozen=True)
class AlignerIdentity:
    """Aligner identity and primary/secondary-record policy."""

    name: str
    version: str
    arguments_sha256: str
    primary_secondary_marking: str


@dataclass(frozen=True)
class CountingPolicy:
    """Complete base-depth policy and exact queried interval union."""

    policy_id: str
    samtools_revision: str
    minimum_mapping_quality: int
    minimum_base_quality: int
    excluded_alignment_flags: tuple[str, ...]
    supplementary_alignment_policy: str
    overlap_policy: str
    base_counting_policy: str
    zero_coverage_policy: str
    queried_intervals: tuple[QueryInterval, ...]


@dataclass(frozen=True)
class LengthFeatureContext:
    """Validated measurement inputs before denominator support is computed."""

    schema_version: str
    manifest_key: str
    input_sha256: str
    assembly: str
    assay_class: str
    input_scope: str
    original_contig: str
    reference_fasta_sha256: str
    annotation_sha256: str
    aligner: AlignerIdentity
    preprocessing_id: str
    counting_policy: CountingPolicy
    sha256: str

    @property
    def counting_policy_sha256(self) -> str:
        """Return the canonical digest of the complete counting policy."""
        return canonical_sha256(_encode_counting_policy(self.counting_policy))


@dataclass(frozen=True)
class DenominatorQc:
    """Distinct read-pair support counts for every A/F numerator or denominator."""

    evidence_kind: DenominatorEvidenceKind
    supporting_fragment_counts: Mapping[str, int | None]

    def __post_init__(self) -> None:
        if isinstance(self.supporting_fragment_counts, Mapping):
            object.__setattr__(
                self, "supporting_fragment_counts", MappingProxyType(dict(self.supporting_fragment_counts))
            )


@dataclass(frozen=True)
class LengthFeatureProvenance:
    """Completed measurement provenance bound to computed denominator evidence."""

    schema_version: str
    measurement_context: LengthFeatureContext
    denominator_qc: DenominatorQc
    sha256: str

    @property
    def manifest_key(self) -> str:
        """Return the intake manifest key."""
        return self.measurement_context.manifest_key

    @property
    def assembly(self) -> str:
        """Return the declared reference assembly."""
        return self.measurement_context.assembly

    @property
    def assay_class(self) -> str:
        """Return the declared assay class."""
        return self.measurement_context.assay_class

    @property
    def input_scope(self) -> str:
        """Return whether the input was full or regional."""
        return self.measurement_context.input_scope

    @property
    def original_contig(self) -> str:
        """Return the exact contig spelling observed in the input."""
        return self.measurement_context.original_contig

    @property
    def reference_fasta_sha256(self) -> str:
        """Return the reference FASTA digest."""
        return self.measurement_context.reference_fasta_sha256

    @property
    def annotation_sha256(self) -> str:
        """Return the annotation digest."""
        return self.measurement_context.annotation_sha256

    @property
    def counting_policy_sha256(self) -> str:
        """Return the canonical digest of the complete counting policy."""
        return self.measurement_context.counting_policy_sha256


def decode_length_feature_context(value: object) -> LengthFeatureContext:
    """Decode a closed pre-measurement context and compute its canonical digest.

    Args:
        value: Parsed JSON-compatible context object.

    Returns:
        The validated immutable context.

    Raises:
        ValueError: If fields, types, hashes, or policy values are invalid.
    """
    root = _exact_object(
        value,
        {
            "schema_version",
            "manifest_key",
            "input_sha256",
            "assembly",
            "assay_class",
            "input_scope",
            "original_contig",
            "reference_fasta_sha256",
            "annotation_sha256",
            "aligner",
            "preprocessing_id",
            "counting_policy",
        },
        "length feature measurement context",
    )
    if root["schema_version"] != LENGTH_FEATURE_CONTEXT_SCHEMA_VERSION:
        raise ValueError(f"length feature context schema version must be {LENGTH_FEATURE_CONTEXT_SCHEMA_VERSION}")
    input_scope = root["input_scope"]
    if input_scope not in {"full", "regional"}:
        raise ValueError("length feature input scope must be full or regional")
    aligner = _decode_aligner(root["aligner"])
    policy = _decode_counting_policy(root["counting_policy"])
    context = LengthFeatureContext(
        schema_version=LENGTH_FEATURE_CONTEXT_SCHEMA_VERSION,
        manifest_key=_text(root["manifest_key"], "length feature manifest key"),
        input_sha256=_digest(root["input_sha256"], "length feature input digest"),
        assembly=_text(root["assembly"], "length feature assembly"),
        assay_class=_text(root["assay_class"], "length feature assay class"),
        input_scope=cast(Literal["full", "regional"], input_scope),
        original_contig=_text(root["original_contig"], "length feature original contig"),
        reference_fasta_sha256=_digest(root["reference_fasta_sha256"], "length feature reference FASTA digest"),
        annotation_sha256=_digest(root["annotation_sha256"], "length feature annotation digest"),
        aligner=aligner,
        preprocessing_id=_text(root["preprocessing_id"], "length feature preprocessing ID"),
        counting_policy=policy,
        sha256="",
    )
    document = _encode_context_unchecked(context)
    return replace(context, sha256=canonical_sha256(document))


def encode_length_feature_context(context: LengthFeatureContext) -> dict[str, object]:
    """Encode and revalidate a measurement context.

    Args:
        context: Context returned by :func:`decode_length_feature_context`.

    Returns:
        The closed JSON-compatible context document.

    Raises:
        ValueError: If a direct or replaced object is invalid or has a stale digest.
    """
    if not isinstance(context, LengthFeatureContext):
        raise ValueError("length feature context must be a LengthFeatureContext")
    if not isinstance(context.aligner, AlignerIdentity) or not isinstance(context.counting_policy, CountingPolicy):
        raise ValueError("length feature context contains invalid typed policy values")
    document = _encode_context_unchecked(context)
    decoded = decode_length_feature_context(document)
    if context != decoded:
        raise ValueError("length feature context digest does not match its canonical content")
    return document


def bind_length_feature_provenance(
    context: LengthFeatureContext, denominator_qc: DenominatorQc
) -> LengthFeatureProvenance:
    """Bind a validated measurement context to extraction-computed QC evidence.

    Args:
        context: Validated pre-measurement context.
        denominator_qc: Denominator counts computed by feature extraction.

    Returns:
        Completed provenance carrying its canonical digest.

    Raises:
        ValueError: If either typed input is invalid or stale.
    """
    context_document = encode_length_feature_context(context)
    qc_document = _encode_denominator_qc(denominator_qc)
    document: dict[str, object] = {
        "schema_version": LENGTH_FEATURE_PROVENANCE_SCHEMA_VERSION,
        "measurement_context": context_document,
        "denominator_qc": qc_document,
    }
    return LengthFeatureProvenance(
        schema_version=LENGTH_FEATURE_PROVENANCE_SCHEMA_VERSION,
        measurement_context=context,
        denominator_qc=denominator_qc,
        sha256=canonical_sha256(document),
    )


def decode_length_feature_provenance(value: object) -> LengthFeatureProvenance:
    """Decode a closed completed provenance document and compute its digest.

    Args:
        value: Parsed JSON-compatible provenance object.

    Returns:
        Validated completed provenance.

    Raises:
        ValueError: If fields, types, hashes, or denominator evidence are invalid.
    """
    root = _exact_object(
        value,
        {"schema_version", "measurement_context", "denominator_qc"},
        "length feature provenance",
    )
    if root["schema_version"] != LENGTH_FEATURE_PROVENANCE_SCHEMA_VERSION:
        raise ValueError(f"length feature provenance schema version must be {LENGTH_FEATURE_PROVENANCE_SCHEMA_VERSION}")
    context = decode_length_feature_context(root["measurement_context"])
    denominator_qc = _decode_denominator_qc(root["denominator_qc"])
    return bind_length_feature_provenance(context, denominator_qc)


def encode_length_feature_provenance(provenance: LengthFeatureProvenance) -> dict[str, object]:
    """Encode and revalidate completed length feature provenance.

    Args:
        provenance: Completed provenance returned by extraction or its decoder.

    Returns:
        The closed JSON-compatible provenance document.

    Raises:
        ValueError: If a direct or replaced object is invalid or has a stale digest.
    """
    if not isinstance(provenance, LengthFeatureProvenance):
        raise ValueError("length feature provenance must be a LengthFeatureProvenance")
    if provenance.schema_version != LENGTH_FEATURE_PROVENANCE_SCHEMA_VERSION:
        raise ValueError(f"length feature provenance schema version must be {LENGTH_FEATURE_PROVENANCE_SCHEMA_VERSION}")
    document: dict[str, object] = {
        "schema_version": provenance.schema_version,
        "measurement_context": encode_length_feature_context(provenance.measurement_context),
        "denominator_qc": _encode_denominator_qc(provenance.denominator_qc),
    }
    decoded = decode_length_feature_provenance(document)
    if provenance != decoded:
        raise ValueError("length feature provenance digest does not match its canonical content")
    return document


def _decode_aligner(value: object) -> AlignerIdentity:
    raw = _exact_object(
        value,
        {"name", "version", "arguments_sha256", "primary_secondary_marking"},
        "length feature aligner",
    )
    return AlignerIdentity(
        name=_text(raw["name"], "length feature aligner name"),
        version=_text(raw["version"], "length feature aligner version"),
        arguments_sha256=_digest(raw["arguments_sha256"], "length feature aligner arguments digest"),
        primary_secondary_marking=_text(raw["primary_secondary_marking"], "length feature primary/secondary marking"),
    )


def _decode_counting_policy(value: object) -> CountingPolicy:
    raw = _exact_object(
        value,
        {
            "policy_id",
            "samtools_revision",
            "minimum_mapping_quality",
            "minimum_base_quality",
            "excluded_alignment_flags",
            "supplementary_alignment_policy",
            "overlap_policy",
            "base_counting_policy",
            "zero_coverage_policy",
            "queried_intervals",
        },
        "length feature counting policy",
    )
    if raw["policy_id"] != COUNTING_POLICY_ID:
        raise ValueError(f"length feature counting policy ID must be {COUNTING_POLICY_ID}")
    for field in ("minimum_mapping_quality", "minimum_base_quality"):
        if not _is_int(raw[field]) or raw[field] != 0:
            raise ValueError(f"length feature {field.replace('_', ' ')} must be integer zero")
    flags = raw["excluded_alignment_flags"]
    if not isinstance(flags, list) or flags != ["UNMAP", "SECONDARY", "QCFAIL", "DUP"]:
        raise ValueError("length feature excluded alignment flags must be UNMAP, SECONDARY, QCFAIL, DUP")
    if raw["supplementary_alignment_policy"] != "included-unless-excluded-by-another-flag":
        raise ValueError("length feature supplementary alignment policy is invalid")
    if raw["overlap_policy"] != "count-overlapping-mates-independently":
        raise ValueError("length feature overlap policy is invalid")
    if raw["base_counting_policy"] != "one-per-aligned-covered-base":
        raise ValueError("length feature base counting policy is invalid")
    if raw["zero_coverage_policy"] != "emit-zero-for-every-queried-position":
        raise ValueError("length feature zero-coverage policy is invalid")
    intervals_raw = raw["queried_intervals"]
    if not isinstance(intervals_raw, list):
        raise ValueError("length feature queried intervals must be a list")
    intervals = tuple(_decode_query_interval(item) for item in intervals_raw)
    for previous, current in zip(intervals, intervals[1:], strict=False):
        if current.start <= previous.end:
            raise ValueError("length feature queried intervals must be sorted, disjoint, and consolidated")
    return CountingPolicy(
        policy_id=COUNTING_POLICY_ID,
        samtools_revision=_text(raw["samtools_revision"], "length feature samtools revision"),
        minimum_mapping_quality=0,
        minimum_base_quality=0,
        excluded_alignment_flags=("UNMAP", "SECONDARY", "QCFAIL", "DUP"),
        supplementary_alignment_policy="included-unless-excluded-by-another-flag",
        overlap_policy="count-overlapping-mates-independently",
        base_counting_policy="one-per-aligned-covered-base",
        zero_coverage_policy="emit-zero-for-every-queried-position",
        queried_intervals=intervals,
    )


def _decode_query_interval(value: object) -> QueryInterval:
    raw = _exact_object(value, {"start", "end"}, "length feature queried interval")
    start = raw["start"]
    end = raw["end"]
    if not _is_int(start) or not _is_int(end):
        raise ValueError("length feature queried interval must be non-empty zero-based half-open integers")
    start_int = cast(int, start)
    end_int = cast(int, end)
    if start_int < 0 or end_int <= start_int:
        raise ValueError("length feature queried interval must be non-empty zero-based half-open integers")
    return QueryInterval(start=start_int, end=end_int)


def _decode_denominator_qc(value: object) -> DenominatorQc:
    raw = _exact_object(value, {"evidence_kind", "supporting_fragment_counts"}, "length feature denominator QC")
    evidence_kind = raw["evidence_kind"]
    if evidence_kind not in {"unavailable", "read-pair-identity-qc-proxy"}:
        raise ValueError("length feature denominator QC evidence kind is invalid")
    counts_raw = _exact_object(
        raw["supporting_fragment_counts"], set(DENOMINATOR_NAMES), "length feature denominator fragment counts"
    )
    counts: dict[str, int | None] = {}
    for name in DENOMINATOR_NAMES:
        count = counts_raw[name]
        if count is not None:
            if not _is_int(count) or cast(int, count) < 0:
                raise ValueError("length feature denominator fragment counts must be null or non-negative integers")
            count = cast(int, count)
        counts[name] = count
    if evidence_kind == "unavailable" and any(count is not None for count in counts.values()):
        raise ValueError("unavailable length feature denominator evidence requires all counts to be null")
    return DenominatorQc(
        evidence_kind=cast(DenominatorEvidenceKind, evidence_kind),
        supporting_fragment_counts=counts,
    )


def _encode_context_unchecked(context: LengthFeatureContext) -> dict[str, object]:
    return {
        "schema_version": context.schema_version,
        "manifest_key": context.manifest_key,
        "input_sha256": context.input_sha256,
        "assembly": context.assembly,
        "assay_class": context.assay_class,
        "input_scope": context.input_scope,
        "original_contig": context.original_contig,
        "reference_fasta_sha256": context.reference_fasta_sha256,
        "annotation_sha256": context.annotation_sha256,
        "aligner": {
            "name": context.aligner.name,
            "version": context.aligner.version,
            "arguments_sha256": context.aligner.arguments_sha256,
            "primary_secondary_marking": context.aligner.primary_secondary_marking,
        },
        "preprocessing_id": context.preprocessing_id,
        "counting_policy": _encode_counting_policy(context.counting_policy),
    }


def _encode_counting_policy(policy: CountingPolicy) -> dict[str, object]:
    if not isinstance(policy.excluded_alignment_flags, tuple) or not isinstance(policy.queried_intervals, tuple):
        raise ValueError("length feature counting policy collections must be immutable tuples")
    if any(not isinstance(interval, QueryInterval) for interval in policy.queried_intervals):
        raise ValueError("length feature queried intervals must contain only QueryInterval values")
    return {
        "policy_id": policy.policy_id,
        "samtools_revision": policy.samtools_revision,
        "minimum_mapping_quality": policy.minimum_mapping_quality,
        "minimum_base_quality": policy.minimum_base_quality,
        "excluded_alignment_flags": list(policy.excluded_alignment_flags),
        "supplementary_alignment_policy": policy.supplementary_alignment_policy,
        "overlap_policy": policy.overlap_policy,
        "base_counting_policy": policy.base_counting_policy,
        "zero_coverage_policy": policy.zero_coverage_policy,
        "queried_intervals": [{"start": interval.start, "end": interval.end} for interval in policy.queried_intervals],
    }


def _encode_denominator_qc(qc: DenominatorQc) -> dict[str, object]:
    if not isinstance(qc, DenominatorQc):
        raise ValueError("length feature denominator QC must be a DenominatorQc")
    if not isinstance(qc.supporting_fragment_counts, Mapping):
        raise ValueError("length feature denominator fragment counts must be a mapping")
    document = {
        "evidence_kind": qc.evidence_kind,
        "supporting_fragment_counts": dict(qc.supporting_fragment_counts),
    }
    return {
        "evidence_kind": _decode_denominator_qc(document).evidence_kind,
        "supporting_fragment_counts": {name: qc.supporting_fragment_counts[name] for name in DENOMINATOR_NAMES},
    }


def _exact_object(value: object, keys: set[str], label: str) -> dict[str, object]:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    actual = set(value)
    if actual != keys:
        raise ValueError(f"{label} fields must be exactly {sorted(keys)}")
    return value


def _digest(value: object, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != _SHA256_LENGTH
        or any(char not in "0123456789abcdef" for char in value)
    ):
        raise ValueError(f"{label} must be a lowercase SHA-256 digest")
    return value


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value.strip() != value:
        raise ValueError(f"{label} must be a non-empty trimmed string")
    return value


def _is_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)
