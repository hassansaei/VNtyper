"""The vendored igv.js: its size, its digests, and what the verification actually covers.

Three separable claims, and the third is the one an audit asks about.

**Size.** The embedded library must cost about 497 KB in the document, not 1.31 MB.
That is the whole argument for storing it gzipped and expanding it in the browser,
and it is a number rather than an opinion.

**Reproducibility.** The gzip is pinned by digest, so it has to be byte-identical on
every rebuild. ``gzip.compress(..., mtime=0)`` is what makes that true: the gzip
header carries a modification time by default, so without it the file's bytes -- and
its digest -- differ every time anyone regenerates it.

**What the digest covers.** Checking the file on disk proves it was not corrupted. It
proves nothing about *which* igv.js is inside it, because any gzip of anything at all
would satisfy a digest taken over the gzip. The check that matters runs over the
**decompressed** bytes, which are exactly the source the report's ``eval`` sees, and
:func:`~vntyper.scripts.report_assets.igv_payload` encodes the same compressed bytes
it decompressed to get there. The last two tests here are what say so: one substitutes
a *valid, well-formed* gzip of different content and requires the payload to be
refused, the other rebuilds the payload independently and checks it round-trips to the
pinned source.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
from pathlib import Path

import pytest

from vntyper.scripts import report_assets

pytestmark = pytest.mark.unit

#: Measured on the fetched artifact, 2026-08-13. Every one of these is a fact about
#: files on disk rather than a target: raw 1,310,337 -> gzip 372,690 -> base64 496,920.
RAW_BYTES = 1_310_337
GZIP_BYTES = 372_690
BASE64_CHARS = 496_920


def test_the_vendored_asset_matches_its_pinned_digest() -> None:
    """The file this package ships is the file whose digest is recorded."""
    report_assets.verify_asset(report_assets.IGV_ASSET_PATH, report_assets.IGV_GZIP_SHA256)
    assert report_assets.IGV_ASSET_PATH.stat().st_size == GZIP_BYTES


def test_verify_asset_refuses_a_file_whose_bytes_moved(tmp_path: Path) -> None:
    """One flipped byte is a different asset, and the report must not embed it."""
    corrupted = tmp_path / "igv.min.js.gz"
    corrupted.write_bytes(report_assets.IGV_ASSET_PATH.read_bytes() + b"\x00")

    with pytest.raises(ValueError, match="does not match its pinned SHA-256"):
        report_assets.verify_asset(corrupted, report_assets.IGV_GZIP_SHA256)


def test_verify_asset_says_so_when_the_asset_is_not_installed(tmp_path: Path) -> None:
    """A distribution built without the asset must fail loudly, naming the packaging.

    This is the release-only failure the packaging guard exists for, stated from the
    other side: if package-data ever drops the file again, this is the message.
    """
    with pytest.raises(ValueError, match="is missing"):
        report_assets.verify_asset(tmp_path / "absent.gz", report_assets.IGV_GZIP_SHA256)


def test_payload_says_so_when_the_vendored_asset_is_not_installed(tmp_path: Path, monkeypatch) -> None:
    """The single-read payload path preserves the package-data failure message."""
    monkeypatch.setattr(report_assets, "IGV_ASSET_PATH", tmp_path / "absent.gz")

    with pytest.raises(ValueError, match="is missing"):
        report_assets.igv_payload(report_assets.REPORT_IGV_EMBEDDED)


def test_the_embedded_library_is_compressed() -> None:
    """Measured: 1,310,337 raw -> 372,690 gzip -> 496,920 base64.

    Inlining the raw library would roughly triple what the alignment browser costs
    the archived file. The bounds are wide enough to survive a re-pin of the version
    and narrow enough to fail if the payload ever became the raw source.
    """
    payload = report_assets.igv_payload(report_assets.REPORT_IGV_EMBEDDED)

    assert payload is not None
    assert 450_000 < len(payload) < 550_000, f"the embedded payload is {len(payload):,} characters"
    assert len(payload) == BASE64_CHARS


@pytest.mark.parametrize("mode", [report_assets.REPORT_IGV_SIDECAR, report_assets.REPORT_IGV_OFF])
def test_only_the_embedded_mode_produces_a_payload(mode: str) -> None:
    """``sidecar`` and ``off`` carry no library in the report itself."""
    assert report_assets.igv_payload(mode) is None


def test_an_unknown_mode_is_refused_rather_than_treated_as_off() -> None:
    """A typo must not silently produce a report with no alignment browser in it."""
    with pytest.raises(ValueError, match="Unknown --report-igv mode 'embeded'"):
        report_assets.igv_payload("embeded")

    with pytest.raises(ValueError, match="Unknown --report-igv mode"):
        report_assets.igv_provenance("embeded")


def test_the_gzip_is_reproducible_from_the_source_it_decompresses_to() -> None:
    """``mtime=0`` is required, not tidy.

    Recompressing the decompressed source at level 9 with ``mtime=0`` must reproduce
    the shipped file byte for byte. Without ``mtime=0`` the gzip header carries the
    time of compression, so the file's digest changes on every rebuild and pinning it
    pins nothing.
    """
    source = report_assets.igv_library_source()
    rebuilt = gzip.compress(source, compresslevel=9, mtime=0)

    assert rebuilt == report_assets.IGV_ASSET_PATH.read_bytes()
    assert hashlib.sha256(rebuilt).hexdigest() == report_assets.IGV_GZIP_SHA256

    stamped = gzip.compress(source, compresslevel=9, mtime=1_700_000_000)
    assert stamped != rebuilt, "the gzip header does not carry mtime here, so this test asserts nothing"


def test_the_decompressed_source_is_the_library_the_provenance_names() -> None:
    """The provenance digest is over the source, and the source is what it says it is."""
    source = report_assets.igv_library_source()

    assert len(source) == RAW_BYTES
    assert hashlib.sha256(source).hexdigest() == report_assets.IGV_SHA256
    assert b"igv" in source[:200_000].lower()


def test_the_payload_decodes_back_to_the_pinned_source() -> None:
    """The chain, end to end: what is embedded expands to what the digest names.

    This is the claim the report's Provenance section makes to a reader, checked the
    way the reader's browser would check it -- base64 decode, gunzip, digest.
    """
    payload = report_assets.igv_payload(report_assets.REPORT_IGV_EMBEDDED)
    assert payload is not None

    expanded = gzip.decompress(base64.b64decode(payload))

    assert hashlib.sha256(expanded).hexdigest() == report_assets.IGV_SHA256
    assert len(expanded) == RAW_BYTES


def test_the_payload_hashes_and_embeds_one_read_of_the_vendored_gzip(monkeypatch) -> None:
    """The gzip digest must cover the exact compressed buffer that is embedded.

    Reading the path once to check the gzip digest and again to build the payload
    leaves a replacement window between those operations. Counting only reads of the
    vendored asset keeps this assertion about the trust chain rather than incidental
    reads elsewhere in the renderer.
    """
    original_read_bytes = Path.read_bytes
    asset_reads = 0

    def _counted_read_bytes(path: Path) -> bytes:
        nonlocal asset_reads
        if path == report_assets.IGV_ASSET_PATH:
            asset_reads += 1
        return original_read_bytes(path)

    monkeypatch.setattr(Path, "read_bytes", _counted_read_bytes)

    payload = report_assets.igv_payload(report_assets.REPORT_IGV_EMBEDDED)

    assert payload is not None
    assert asset_reads == 1, "the gzip pin and embedded payload came from different reads"


def test_a_well_formed_gzip_of_the_wrong_library_is_refused(tmp_path: Path, monkeypatch) -> None:
    """The check that makes the digest mean something rather than merely pass.

    A digest taken over the gzip alone is satisfied by *any* gzip whose bytes happen
    to match -- and the way that goes wrong in practice is not corruption but
    substitution: someone regenerates the asset from a different igv.js and re-pins
    the gzip digest without touching the source digest. So this substitutes a valid,
    well-formed, correctly-digested gzip of different JavaScript and requires
    ``igv_payload`` to refuse it on the decompressed bytes.

    Without the second check this test fails: the file verifies, the payload is
    produced, and the report embeds a library nobody vouched for while printing a
    Provenance line that says otherwise.
    """
    impostor_source = b"window.igv = {};\n" * 1000
    impostor = tmp_path / "igv-3.0.2.min.js.gz"
    impostor.write_bytes(gzip.compress(impostor_source, compresslevel=9, mtime=0))

    monkeypatch.setattr(report_assets, "IGV_ASSET_PATH", impostor)
    monkeypatch.setattr(
        report_assets,
        "IGV_GZIP_SHA256",
        hashlib.sha256(impostor.read_bytes()).hexdigest(),
    )

    # The file itself now passes: this is a genuine substitution, not a corrupt file.
    report_assets.verify_asset(report_assets.IGV_ASSET_PATH, report_assets.IGV_GZIP_SHA256)

    with pytest.raises(ValueError, match="Decompressed igv.js 3.0.2 does not match its pinned SHA-256"):
        report_assets.igv_payload(report_assets.REPORT_IGV_EMBEDDED)


def test_the_provenance_line_states_the_version_and_the_source_digest() -> None:
    """What the report prints, and the reader checks upstream.

    The digest is printed in full. A truncated one is something a reader can compare
    at a glance and cannot actually verify, which is the wrong half of the trade for
    a line whose entire purpose is verification.
    """
    embedded = report_assets.igv_provenance(report_assets.REPORT_IGV_EMBEDDED)
    assert report_assets.IGV_VERSION in embedded
    assert report_assets.IGV_SHA256 in embedded
    assert "embedded in this file" in embedded

    sidecar = report_assets.igv_provenance(report_assets.REPORT_IGV_SIDECAR)
    assert report_assets.IGV_SHA256 in sidecar
    assert "igv_report.html" in sidecar

    off = report_assets.igv_provenance(report_assets.REPORT_IGV_OFF)
    assert report_assets.IGV_SHA256 not in off, "a report with no library must not claim one by digest"
    assert "off" in off


def test_the_modes_are_the_three_the_cli_offers() -> None:
    """The tuple the parser's ``choices`` is built from, and the default it uses."""
    assert report_assets.REPORT_IGV_MODES == ("embedded", "sidecar", "off")
    assert report_assets.DEFAULT_REPORT_IGV == "embedded"
    assert report_assets.DEFAULT_REPORT_IGV in report_assets.REPORT_IGV_MODES


def test_the_notice_records_the_source_url_and_both_digests() -> None:
    """Redistributing a third-party library means saying which one, and under what terms."""
    notice = (report_assets.ASSET_DIR / "NOTICE.md").read_text(encoding="utf-8")

    assert report_assets.IGV_SOURCE_URL in notice
    assert report_assets.IGV_SHA256 in notice
    assert report_assets.IGV_GZIP_SHA256 in notice
    assert "The MIT License (MIT)" in notice
    assert "Permission is hereby granted, free of charge" in notice
