"""Pure policy tests for local and remote CRAM reference URIs."""

from __future__ import annotations

from pathlib import Path

import pytest

from vntyper.scripts import reference_uri_policy
from vntyper.scripts.reference_uri_policy import (
    RemoteHeaderReference,
    enforce_header_reference_policy,
    first_remote_header_reference,
    ref_path_remote_scheme,
    remote_ref_path_suffix,
)

pytestmark = pytest.mark.unit


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("https://refget.example/reference.fa", "https"),
        ("S3://reference-bucket/reference.fa", "s3"),
        ("/local/%s:http://127.0.0.1/reference.fa", "http"),
        ("reference.fa", None),
        ("../reference.fa", None),
        ("/local/reference.fa", None),
        ("file:///local/reference.fa", None),
        ("file://reference-host/local/reference.fa", None),
    ],
)
def test_ref_path_remote_scheme_distinguishes_network_and_local_search_entries(
    value: str, expected: str | None
) -> None:
    """A non-file URI scheme is the only shape governed by the network waiver."""
    assert ref_path_remote_scheme(value) == expected


def test_remote_ref_path_suffix_removes_local_cache_entries_before_ambient_lookup() -> None:
    assert remote_ref_path_suffix("/operator/%s:https://refget.example/%s") == "https://refget.example/%s"
    assert remote_ref_path_suffix("/operator/%s") is None


def test_first_remote_header_reference_reports_the_sq_contig_and_normalised_scheme() -> None:
    """Tag order and unrelated header lines cannot hide a remote SQ reference."""
    header = (
        "@HD\tVN:1.6\n"
        "@SQ\tUR:file:///local/chr1.fa\tBROKEN\tSN:chr1\tLN:100\n"
        "@SQ\tM5:digest\tUR:HTTPS://refget.example/secret/path.fa\tLN:200\tSN:chr2\n"
    )

    assert first_remote_header_reference(header) == RemoteHeaderReference(contig="chr2", scheme="https")


def test_remote_header_reference_percent_encodes_path_like_contig_characters() -> None:
    header = "@SQ\tSN:chr1/alt\\patch\tLN:10\tUR:https://reference.invalid/ref.fa\n"

    assert first_remote_header_reference(header) == RemoteHeaderReference(
        contig="chr1%2Falt%5Cpatch",
        scheme="https",
    )


def test_every_duplicate_ur_field_is_inspected_instead_of_the_last_one_winning() -> None:
    """A trailing local UR cannot hide an earlier remote reference on the same SQ line."""
    header = "@SQ\tSN:chr7\tLN:100\tUR:http://refget.example/private.fa\tUR:file:///local/reference.fa\n"

    assert first_remote_header_reference(header) == RemoteHeaderReference(contig="chr7", scheme="http")


def test_remote_ur_without_an_sq_name_uses_the_path_free_unknown_contig() -> None:
    """A malformed SQ line still rejects its remote target without exposing that target."""
    header = "@SQ\tLN:100\tUR:https://refget.example/private.fa\n"

    assert first_remote_header_reference(header) == RemoteHeaderReference(contig="unknown", scheme="https")


@pytest.mark.parametrize("uri", ["http:refget.example/private.fa", "S3:reference-bucket/reference.fa"])
def test_header_uri_schemes_are_recognised_without_slashes(uri: str) -> None:
    """RFC-style htslib schemes are anchored at the URI start and do not require ``//``."""
    header = f"@SQ\tSN:chr2\tLN:100\tUR:{uri}\n"

    assert first_remote_header_reference(header) == RemoteHeaderReference(
        contig="chr2", scheme=uri.partition(":")[0].lower()
    )


def test_header_and_ref_path_detection_have_context_specific_public_interfaces() -> None:
    """Header URI syntax and colon-separated REF_PATH syntax cannot share ambiguous parsing."""
    assert reference_uri_policy.header_reference_scheme("http:opaque-reference") == "http"
    assert reference_uri_policy.header_reference_scheme("/local/cache:http://refget.example/%s") is None
    assert reference_uri_policy.ref_path_remote_scheme("/local/cache:http://refget.example/%s") == "http"
    assert reference_uri_policy.ref_path_remote_scheme("file:///local/reference.fa:/local/cache/%s") is None


@pytest.mark.parametrize(
    "uri",
    ["reference.fa", "../reference.fa", "/reference/reference.fa", "file:///reference/reference.fa"],
)
def test_default_header_policy_allows_local_relative_and_file_uris(uri: str) -> None:
    """Every explicitly local header form remains available in default mode."""
    enforce_header_reference_policy(f"@SQ\tSN:chr1\tLN:100\tUR:{uri}\n", allow_ambient=False)


def test_explicit_ambient_waiver_allows_a_remote_header_uri() -> None:
    """The existing opt-in preserves operator-controlled network resolution."""
    enforce_header_reference_policy(
        "@SQ\tSN:chr1\tLN:100\tUR:https://refget.example/reference.fa\n",
        allow_ambient=True,
    )


def test_default_header_policy_rejects_remote_uri_without_disclosing_its_path() -> None:
    """The diagnostic identifies the decision while omitting the remote target path."""
    uri = "http://127.0.0.1:8765/private/reference.fa"

    with pytest.raises(ValueError) as raised:
        enforce_header_reference_policy(f"@SQ\tSN:chr7\tLN:100\tUR:{uri}\n", allow_ambient=False)

    message = str(raised.value)
    assert "contig=chr7" in message
    assert "scheme=http" in message
    assert "allow_ambient_reference_resolution=true" in message
    assert "/" not in message
    assert "file-scheme reference" in message
    assert uri not in message
    assert "/private/reference.fa" not in message


def test_local_header_reference_paths_resolve_against_the_input_and_deduplicate(tmp_path: Path) -> None:
    """Local UR spellings become one ordered filesystem candidate per lexical path."""
    input_cram = tmp_path / "patient" / "input.cram"
    relative_reference = input_cram.parent / "reference one.fa"
    absolute_reference = tmp_path / "reference two.fa"
    header = (
        "@SQ\tSN:chr1\tUR:reference one.fa\tLN:100\n"
        f"@SQ\tSN:chr2\tUR:{relative_reference.as_uri()}\tLN:100\n"
        f"@SQ\tSN:chr3\tUR:{absolute_reference.as_uri()}\tLN:100\n"
        "@SQ\tSN:chr4\tUR:https://reference.invalid/remote.fa\tLN:100\n"
    )

    assert reference_uri_policy.local_header_reference_paths(header, input_cram) == (
        str(relative_reference),
        str(absolute_reference),
    )


def test_local_header_references_preserve_each_sq_contig_digest_association(tmp_path: Path) -> None:
    """Per-contig local FASTAs remain associated with the M5 they can populate."""
    header = "@SQ\tSN:chr1\tM5:first\tUR:chr1.fa\n@SQ\tSN:chr2\tM5:second\tUR:chr2.fa\n"

    assert reference_uri_policy.local_header_references(header, tmp_path / "input.cram") == (
        reference_uri_policy.LocalHeaderReference("chr1", "first", str(tmp_path / "chr1.fa")),
        reference_uri_policy.LocalHeaderReference("chr2", "second", str(tmp_path / "chr2.fa")),
    )


def test_bare_header_reference_percent_sequences_remain_literal(tmp_path: Path) -> None:
    """Percent decoding applies only to file URIs, never ordinary filesystem names."""
    input_cram = tmp_path / "input.cram"

    assert reference_uri_policy.local_header_reference_paths("@SQ\tSN:chr1\tUR:ref%20one.fa\n", input_cram) == (
        str(tmp_path / "ref%20one.fa"),
    )


@pytest.mark.parametrize(
    "uri",
    [
        "file://reference-host/reference.fa",
        "file:///reference.fa?version=1",
        "file:///reference.fa#fragment",
        "file:///reference%00.fa",
        "file:",
    ],
)
def test_ambiguous_or_nonlocal_file_header_uris_are_rejected(tmp_path: Path, uri: str) -> None:
    """A file URI must identify one unambiguous local filesystem path."""
    with pytest.raises(ValueError, match="Invalid local CRAM header reference URI") as raised:
        reference_uri_policy.local_header_reference_paths(f"@SQ\tSN:chr1\tUR:{uri}\n", tmp_path / "input.cram")

    assert uri not in str(raised.value)
