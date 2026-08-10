"""Pure evidence contract for the registered reference-compressed CRAM fixture."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

HG19_CHR1_REFERENCE_SHA256 = "0c19925c13b1312f0cbdc2b804f62da260345589b8f9e8ad655abfb5d6e99338"
REFERENCE_VALIDATION_REGION = "chr1:155160500-155162000"
REGISTERED_B178_INDEXED_REGION_RECORDS = 13_868
REGISTERED_B178_INDEXED_REGION_DIGEST = "748b2aea22748d0277857b5ecbe641f85fed0e3dfd280aa62def0b71728e07a0"

# These identities are bound to the committed UCSC provenance artifact by a whole-file
# digest test. Keeping the runtime table path-free lets direct script execution work from
# any current directory without weakening the source attestation.
HG19_PRIMARY_CONTIGS: dict[str, tuple[int, str]] = {
    "chr1": (249250621, "1b22b98cdeb4a9304cb5d48026a85128"),
    "chr2": (243199373, "a0d9851da00400dec1098a9255ac712e"),
    "chr3": (198022430, "641e4338fa8d52a5b781bd2a2c08d3c3"),
    "chr4": (191154276, "23dccd106897542ad87d2765d28a19a1"),
    "chr5": (180915260, "0740173db9ffd264d728f32784845cd7"),
    "chr6": (171115067, "1d3a93a248d92a729ee764823acbbc6b"),
    "chr7": (159138663, "618366e953d6aaad97dbe4777c29375e"),
    "chr8": (146364022, "96f514a9929e410c6651697bded59aec"),
    "chr9": (141213431, "3e273117f15e0a400f01055d9f393768"),
    "chr10": (135534747, "988c28e000e84c26d552359af1ea2e1d"),
    "chr11": (135006516, "98c59049a2df285c76ffb1c6db8f8b96"),
    "chr12": (133851895, "51851ac0e1a115847ad36449b0015864"),
    "chr13": (115169878, "283f8d7892baa81b510a015719ca7b0b"),
    "chr14": (107349540, "98f3cae32b2a2e9524bc19813927542e"),
    "chr15": (102531392, "e5645a794a8238215b2cd77acb95a078"),
    "chr16": (90354753, "fc9b1a7b42b97a864f56b348b06095e6"),
    "chr17": (81195210, "351f64d4f4f9ddd45b35336ad97aa6de"),
    "chr18": (78077248, "b15d4b2d29dde9d3e4f93d1d0f2cbc9c"),
    "chr19": (59128983, "1aacd71f30db8e561810913e0b72636d"),
    "chr20": (63025520, "0dec9660ec1efaaf33281c0d5ea2560f"),
    "chr21": (48129895, "2979a6085bfe28e3ad6f552f361ed74d"),
    "chr22": (51304566, "a718acaa6135fdca8357d5bfe94211dd"),
    "chrX": (155270560, "7e0e2e580297b7764e31dbc80c2540dd"),
    "chrY": (59373566, "1e86411d73e6f00a10590f976be01623"),
    "chrM": (16571, "d2ed829b8a1628d16cbeee88e88e39eb"),
}


class LossyConversionError(RuntimeError):
    """A derived CRAM did not satisfy the record-equivalence contract."""


@dataclass(frozen=True)
class ReferenceCompressedFixture:
    """One explicitly referenced CRAM and its reproducibility evidence."""

    source_bam: Path
    cram: Path
    reference: Path
    records: int
    source_record_digest: str
    decoded_record_digest: str
    indexed_region: str
    indexed_region_records: int
    source_indexed_region_digest: str
    decoded_indexed_region_digest: str
    reference_sha256: str
    source_bytes: int
    cram_bytes: int

    def as_manifest_entry(self) -> dict[str, object]:
        """Return stable JSON-ready provenance for the fixture."""
        return {
            "encoding": "reference-compressed",
            "source_bam": str(self.source_bam),
            "cram": str(self.cram),
            "reference_fasta": str(self.reference),
            "records": self.records,
            "source_record_digest": self.source_record_digest,
            "decoded_record_digest": self.decoded_record_digest,
            "indexed_region": self.indexed_region,
            "indexed_region_records": self.indexed_region_records,
            "source_indexed_region_digest": self.source_indexed_region_digest,
            "decoded_indexed_region_digest": self.decoded_indexed_region_digest,
            "reference_sha256": self.reference_sha256,
            "source_bytes": self.source_bytes,
            "cram_bytes": self.cram_bytes,
        }


def normalize_sam_record(line: str) -> bytes:
    """Return stable SAM bytes with semantically unordered optional fields sorted."""
    fields = line.rstrip("\n").split("\t")
    return ("\t".join(fields[:11] + sorted(fields[11:])) + "\n").encode()


def header_with_hg19_m5(header_text: str) -> str:
    """Add verified UCSC hg19 M5 tags without changing the sequence dictionary.

    Raises:
        ValueError: If a sequence is unknown, has the wrong length, or conflicts
            with its pinned M5 identity.
    """
    output: list[str] = []
    sequences = 0
    for line in header_text.splitlines():
        if not line.startswith("@SQ\t"):
            output.append(line)
            continue
        fields = line.split("\t")
        tags = {field[:2]: field[3:] for field in fields[1:] if len(field) > 3 and field[2] == ":"}
        name = tags.get("SN", "")
        if name not in HG19_PRIMARY_CONTIGS:
            raise ValueError(f"Reference-compressed source header has unsupported hg19 sequence: {name or line}")
        expected_length, expected_m5 = HG19_PRIMARY_CONTIGS[name]
        if tags.get("LN") != str(expected_length):
            raise ValueError(
                f"Reference-compressed source header length mismatch for {name}: "
                f"expected {expected_length}, observed {tags.get('LN', '<missing>')}"
            )
        if "M5" in tags and tags["M5"] != expected_m5:
            raise ValueError(
                f"Reference-compressed source header M5 mismatch for {name}: "
                f"expected {expected_m5}, observed {tags['M5']}"
            )
        if "M5" not in tags:
            fields.append(f"M5:{expected_m5}")
        output.append("\t".join(fields))
        sequences += 1
    if sequences == 0:
        raise ValueError("Reference-compressed source header contains no @SQ sequences")
    return "\n".join(output) + "\n"


def validate_registered_b178_index_evidence(
    source_digest: str,
    source_records: int,
    decoded_digest: str,
    decoded_records: int,
) -> None:
    """Require the pinned nonempty source locus and an identical CRAI query result.

    Raises:
        LossyConversionError: If the source is not registered b178 evidence, or the
            decoded indexed query differs from that source evidence.
    """
    if source_records != REGISTERED_B178_INDEXED_REGION_RECORDS:
        raise LossyConversionError(
            f"Indexed proof expected {REGISTERED_B178_INDEXED_REGION_RECORDS} registered source records "
            f"at {REFERENCE_VALIDATION_REGION}, observed {source_records}; source is empty or is not the registered source locus"
        )
    if source_digest != REGISTERED_B178_INDEXED_REGION_DIGEST:
        raise LossyConversionError(
            f"Indexed proof source digest does not identify the registered source locus at {REFERENCE_VALIDATION_REGION}: "
            f"expected {REGISTERED_B178_INDEXED_REGION_DIGEST}, observed {source_digest}"
        )
    if decoded_records != source_records or decoded_digest != source_digest:
        raise LossyConversionError(
            f"indexed region is not lossless at {REFERENCE_VALIDATION_REGION}: "
            f"{source_records} source records digest {source_digest[:16]}, "
            f"{decoded_records} decoded records digest {decoded_digest[:16]}"
        )
