"""
report_assets.py

Module Purpose:
---------------
The vendored assets the per-sample report can carry, and the verification that
decides whether they are allowed into it.

There is exactly one asset: igv.js, pinned at 3.0.2, stored gzipped under
``vntyper/assets/`` and base64-encoded into the document at render time. The
report used to fetch it from a CDN, along with Bootstrap, jQuery and DataTables,
which meant an archived report read differently -- and in three places said
something different -- depending on whether the reader's machine could reach four
third-party hosts years later (#242).

Why gzip and base64 rather than the raw file
--------------------------------------------
Measured, not assumed. The raw library is 1,310,337 bytes; gzipped it is 372,690;
base64 of that gzip is 496,920 characters, which is 37.9% of the raw size. Inlining
the raw library would roughly triple what the embedded browser costs the archive.
The expansion cost in the reader's browser is one ``DecompressionStream``, which
every engine from Chrome 80 / Safari 16.4 / Firefox 113 onwards has.

Why the version is pinned *down*
--------------------------------
3.8.5 measures 575,932 base64 characters -- 79 KB **worse**. The newest release is
not the right one for an artifact whose defining property is that it is small
enough to archive and to forward, so 3.0.2 stays until something other than
recency argues for a move.

The two digests, and which one matters
--------------------------------------
Both are pinned and both are checked, because they answer different questions:

* :data:`IGV_GZIP_SHA256` is what the file on disk carries. Checking it proves the
  file was not corrupted or truncated, and nothing else.
* :data:`IGV_SHA256` is the digest of the **decompressed** library -- the exact
  bytes a reader's browser later evaluates, and the number a reader can check
  against upstream. It is the one printed in the report's Provenance section.

Verifying only the file on disk would prove the gzip is intact while proving
nothing about *which* igv.js is inside it. :func:`igv_payload` therefore decompresses
and verifies the source before encoding the compressed bytes it came from, so the
check covers exactly what reaches the document -- and, from there, ``eval``.

The gzip is reproducible: start with
``gzip.compress(raw, compresslevel=9, mtime=0)`` and set byte 9 of the result to
``3`` (RFC 1952's Unix OS value). ``mtime=0`` fixes the timestamp; setting the OS
byte fixes a Python-version difference where supported interpreters emit either
Unix (3) or unknown (255). Without both steps, a pinned digest can vary with the
build environment rather than only with the library source.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

#: The pinned igv.js release. Printed in the report and used to name the asset file.
IGV_VERSION = "3.0.2"

#: Where the vendored bytes were fetched from, recorded so the digests below can be
#: checked against upstream without guessing which artifact was taken.
IGV_SOURCE_URL = "https://cdn.jsdelivr.net/npm/igv@3.0.2/dist/igv.min.js"

#: SHA-256 of the **decompressed** library: 1,310,337 bytes of ``igv.min.js`` exactly
#: as fetched. This is the provenance claim a reader can verify against upstream, and
#: it is what :func:`igv_payload` checks the payload against before encoding it.
IGV_SHA256 = "ab1aa79c514ee3a0d66a0ffc788b6d37803910e62cf6d114d9b2909d96b5e790"

#: SHA-256 of the gzipped file this package ships: 372,690 bytes.
IGV_GZIP_SHA256 = "0d8b512654b2ef588009453c403c8f1329dce88eedca90ba9e60888af6b2f79f"

#: The vendored assets directory, inside the installed package. It is declared in
#: pyproject's package-data and in MANIFEST.in; ``tests/unit/test_packaging_consistency.py``
#: fails the build for a file here that is in neither.
ASSET_DIR = Path(__file__).resolve().parent.parent / "assets"

#: The controlled template passed to igv-reports. It carries no external resource
#: URL; :func:`~vntyper.scripts.igv_report.run_igv_report` inserts the verified
#: library into the generated sidecar after igv-reports has populated its data.
IGV_REPORT_TEMPLATE_PATH = Path(__file__).resolve().parent.parent / "templates" / "igv_report_template.html"

#: The packaged report templates, resolved from the installed module instead of
#: the process CWD. Package-data and sdist coverage are enforced by
#: ``test_packaging_consistency.py``.
PACKAGE_TEMPLATE_DIR = Path(__file__).resolve().parent.parent / "templates"

#: The historical ``config.json`` spelling for the packaged directory. It is not
#: an operator override: following it relative to the CWD made installed report
#: generation fail and allowed an unrelated writable directory to replace the
#: report entry template.
LEGACY_TEMPLATE_DIR_LITERAL = "vntyper/templates"

#: The one literal replaced in a generated sidecar. Keeping it distinctive makes a
#: changed or partial template fail closed instead of emitting a browser with no library.
IGV_REPORT_LIBRARY_MARKER = "@VNTYPER_IGV_LIBRARY@"

#: The vendored library, gzipped.
IGV_ASSET_PATH = ASSET_DIR / f"igv-{IGV_VERSION}.min.js.gz"

#: ``--report-igv embedded``: the library travels inside the report, so the archived
#: file is a complete alignment browser with no second file and no network.
REPORT_IGV_EMBEDDED = "embedded"

#: ``--report-igv sidecar``: the report carries no library, and the self-contained
#: ``igv_report.html`` built from VNtyper's local template and vendored library is the
#: alignment browser. Chosen when the ~500 KB matters more than the single file.
REPORT_IGV_SIDECAR = "sidecar"

#: ``--report-igv off``: no alignment browser is produced at all.
REPORT_IGV_OFF = "off"

#: Every accepted ``--report-igv`` value, in the order the CLI lists them.
REPORT_IGV_MODES = (REPORT_IGV_EMBEDDED, REPORT_IGV_SIDECAR, REPORT_IGV_OFF)

#: What a run does when nothing asked for anything else. Embedding is the default
#: because a self-contained archive is the point of #242; the other two modes exist
#: for operators who have measured that they want something else.
DEFAULT_REPORT_IGV = REPORT_IGV_EMBEDDED


def template_search_paths(
    template_dir: str | Path | None,
    *,
    entry_template: str | None = None,
) -> list[str]:
    """Resolve report-template search paths without depending on the process CWD.

    The absent setting, an empty setting, and the shipped historical literal all
    select the installed package directory. An explicit operator directory is
    searched first, with packaged partials available as fallback. When a renderer
    names its entry template, that file must exist in the operator directory so a
    typo cannot silently substitute the packaged top-level report.

    Args:
        template_dir: Configured ``paths.template_dir`` value, or None.
        entry_template: Required top-level template name for an explicit override.

    Returns:
        list[str]: Jinja2 search paths in priority order.

    Raises:
        ValueError: If an explicit override lacks ``entry_template``.
    """
    configured = "" if template_dir is None else str(template_dir)
    packaged = str(PACKAGE_TEMPLATE_DIR)
    if configured in ("", LEGACY_TEMPLATE_DIR_LITERAL, packaged):
        return [packaged]

    if entry_template is not None and not (Path(configured) / entry_template).is_file():
        msg = (
            f"operator template directory '{template_dir}' does not contain required entry template '{entry_template}'."
        )
        logger.error(msg)
        raise ValueError(msg)

    logger.info(f"Using operator template directory '{template_dir}' with packaged templates as fallback.")
    return [configured, packaged]


def _verify_bytes(payload: bytes, expected_sha256: str, description: str) -> None:
    """Fail unless ``payload`` hashes to ``expected_sha256``.

    Args:
        payload: The bytes to digest.
        expected_sha256: The pinned lowercase hex digest.
        description: What these bytes are, for the message on failure.

    Raises:
        ValueError: If the digest does not match.
    """
    actual = hashlib.sha256(payload).hexdigest()
    if actual != expected_sha256:
        msg = (
            f"{description} does not match its pinned SHA-256: expected {expected_sha256}, "
            f"got {actual} over {len(payload)} bytes. The vendored asset has been modified "
            "or replaced; nothing is embedded into a report until it matches."
        )
        logger.error(msg)
        raise ValueError(msg)
    logger.debug("%s matches its pinned SHA-256 (%d bytes).", description, len(payload))


def verify_asset(path: Path, expected_sha256: str) -> None:
    """Fail unless the file at ``path`` hashes to ``expected_sha256``.

    Args:
        path: The vendored file to check.
        expected_sha256: The pinned lowercase hex digest.

    Raises:
        ValueError: If the file is missing or its digest does not match.
    """
    if not path.is_file():
        msg = (
            f"Vendored asset {path} is missing. A distribution built without it cannot embed "
            "the alignment browser; check pyproject's package-data and MANIFEST.in."
        )
        logger.error(msg)
        raise ValueError(msg)
    _verify_bytes(path.read_bytes(), expected_sha256, f"Vendored asset {path.name}")


def _read_verified_asset() -> tuple[bytes, bytes]:
    """Read the vendored asset and verify both forms of it.

    Both digests are checked here, in the order that makes the second one mean
    something: the file on disk against :data:`IGV_GZIP_SHA256`, and then what comes
    out of it against :data:`IGV_SHA256`. The second check is the one that says
    *which* igv.js this is, and it runs on the bytes the report's ``eval`` will see.

    The compressed bytes are returned alongside the source, and they are the ones
    that get encoded into the document. That is what closes the chain: these exact
    compressed bytes were shown to decompress to the verified source, so encoding
    them embeds the thing that was checked rather than something re-derived from it.

    Returns:
        tuple[bytes, bytes]: The compressed asset and its decompressed source.

    Raises:
        ValueError: If either digest does not match, the asset is missing, or the
            verified source could close the sidecar's raw-text script element.
    """
    if not IGV_ASSET_PATH.is_file():
        msg = (
            f"Vendored asset {IGV_ASSET_PATH} is missing. A distribution built without it cannot embed "
            "the alignment browser; check pyproject's package-data and MANIFEST.in."
        )
        logger.error(msg)
        raise ValueError(msg)
    compressed = IGV_ASSET_PATH.read_bytes()
    _verify_bytes(compressed, IGV_GZIP_SHA256, f"Vendored asset {IGV_ASSET_PATH.name}")
    source = gzip.decompress(compressed)
    _verify_bytes(source, IGV_SHA256, f"Decompressed igv.js {IGV_VERSION}")
    if b"</script" in source.lower():
        msg = (
            f"Decompressed igv.js {IGV_VERSION} contains a case-insensitive script-closing sequence and cannot "
            "be inserted as raw text into the controlled IGV sidecar template."
        )
        logger.error(msg)
        raise ValueError(msg)
    return compressed, source


def igv_library_source() -> bytes:
    """Return the verified, decompressed igv.js source.

    Returns:
        bytes: The decompressed library.

    Raises:
        ValueError: If either digest does not match, the asset is missing, or the
            verified source cannot be inserted safely into the sidecar.
    """
    return _read_verified_asset()[1]


def igv_payload(mode: str) -> str | None:
    """Return the base64 of the gzipped igv.js, or None when the mode wants no library.

    The compressed bytes are encoded only after their decompressed form has been
    verified against :data:`IGV_SHA256`, so the digest printed in the report's
    Provenance section covers exactly the source the document later evaluates rather
    than merely the file that happened to be on disk.

    Args:
        mode: One of :data:`REPORT_IGV_MODES`.

    Returns:
        str | None: The base64 payload for ``embedded``; None for ``sidecar`` and
        ``off``, which carry no library in the report itself.

    Raises:
        ValueError: If ``mode`` is not a recognised mode, or the asset fails
            verification.
    """
    if mode not in REPORT_IGV_MODES:
        msg = f"Unknown --report-igv mode {mode!r}; expected one of {', '.join(REPORT_IGV_MODES)}."
        logger.error(msg)
        raise ValueError(msg)
    if mode != REPORT_IGV_EMBEDDED:
        logger.debug("--report-igv %s: the report carries no alignment library.", mode)
        return None

    compressed, source = _read_verified_asset()
    payload = base64.b64encode(compressed).decode("ascii")
    logger.info(
        "Embedding igv.js %s: %d source bytes -> %d base64 characters (%.1f%%).",
        IGV_VERSION,
        len(source),
        len(payload),
        100.0 * len(payload) / len(source),
    )
    return payload


def igv_provenance(mode: str) -> str:
    """Describe, for the report's Provenance section, what carries the alignment browser.

    Args:
        mode: One of :data:`REPORT_IGV_MODES`.

    Returns:
        str: A line naming the version and the verified source digest, and saying
        where the library is. The digest is printed in full: a truncated one is a
        thing a reader can compare at a glance and cannot actually check.

    Raises:
        ValueError: If ``mode`` is not a recognised mode.
    """
    if mode not in REPORT_IGV_MODES:
        msg = f"Unknown --report-igv mode {mode!r}; expected one of {', '.join(REPORT_IGV_MODES)}."
        logger.error(msg)
        raise ValueError(msg)
    if mode == REPORT_IGV_OFF:
        return "not included (--report-igv off)"
    where = (
        "embedded in this file, gzipped"
        if mode == REPORT_IGV_EMBEDDED
        else "in igv_report.html beside this report (--report-igv sidecar)"
    )
    return f"{IGV_VERSION} · sha256 {IGV_SHA256} · {where}"
