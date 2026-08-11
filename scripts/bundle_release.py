#!/usr/bin/env python3
"""Turn a built reference tree into the assets of an immutable release.

`.github/workflows/build-reference-bundles.yml` runs `install-references --from-source`
to produce a complete reference tree and then calls this script to group that tree into
seven `.tar.gz` assets, describe them, and prove that every digest it records was taken
from the bytes it names.

Why this is written to fail rather than to cope
-----------------------------------------------
A published GitHub release is immutable. A bundle that omits a file, ships one twice, or
records a digest that does not describe what it names cannot be patched: it forces a new
release version, a full rebuild and reverify, and new digests in every consumer pin. So
the governing rule here is a **closed inventory** - every file under the reference tree
is assigned to exactly one asset or matches an explicit, reasoned exclusion, and anything
else stops the build by name. Two omissions (`filter_config.json` never being staged and
`vntr_db_advntr.zip` being left unassigned) were found by hand before that rule existed.

Where the metadata goes, and why it is not a loose file
-------------------------------------------------------
`release-manifest.json` and `BUILD_INFO.json` are written **inside every archive**, not
only beside them. Each asset's SHA-256 is committed in
`vntyper/scripts/install_references_config.json`, which is the consumer's trust anchor; a
metadata file downloaded separately has no committed digest to check it against, so an
actor who can replace a release asset can replace the loose JSON alongside it and the
installer would believe both. Inside the archive, the metadata is covered by the
archive's committed digest.

The release *also* publishes top-level `SHA256SUMS`, `release-manifest.json`,
`BUILD_INFO.json` and `verification-report.json`. Those are for the maintainer reviewing
the draft; the installer never reads them.

One consequence for the installer, which is easier to design around than to debug: every
archive carries its metadata at the *same* two paths, so extracting several assets into
one directory leaves only the last archive's `release-manifest.json` behind.
`BUILD_INFO.json` is identical in all of them and does not care, but a manifest read out
of a merged tree describes one asset and not the others. Extract each asset into its own
staging directory, or read the manifest out of the archive; the manifest's `asset` field
names the archive it belongs to, so a consumer that gets this wrong can detect it.
"""

from __future__ import annotations

import argparse
import fnmatch
import gzip
import hashlib
import io
import json
import logging
import platform
import re
import subprocess
import sys
import tarfile
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

_CHUNK = 1024 * 1024


def sha256_of(path: Path) -> str:
    """Return the hex SHA-256 of a file, read in chunks.

    A deliberate three-line copy of :func:`vntyper.scripts.reference_bundle.sha256_of`,
    not an import of it. The workflow runs this file by path out of a *second* checkout
    (`python vntyper/scripts/bundle_release.py`), where `sys.path[0]` is `scripts/` and
    not the repository root, so `import vntyper` resolves only if the package happens to
    be installed - which makes a release build depend on step ordering rather than on its
    own arguments. `TestDigestAgreement` in `tests/unit/test_bundle_release.py` pins that
    the two helpers return the same digest for the same bytes.

    Args:
        path: File to digest.

    Returns:
        str: Lowercase hex digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


#: Every asset is named `vntyper-references-<tag>-<suffix>.tar.gz`.
ASSET_PREFIX = "vntyper-references"

#: The trust anchor, read **by path** rather than by importing `vntyper`. The workflow
#: runs this file as `python vntyper/scripts/bundle_release.py` out of a second checkout,
#: where `sys.path[0]` is `scripts/` and `import vntyper` resolves only if the package
#: happens to be installed - which would make the check below depend on step ordering.
#: Assigned to a module attribute so a test can point it at a fixture.
COMMITTED_CONFIG = Path(__file__).resolve().parents[1] / "vntyper" / "scripts" / "install_references_config.json"

#: Where `install_references_config.json` keeps its downloadable genomes.
CONFIG_GENOME_SECTIONS = ("ucsc_references", "ncbi_references", "ensembl_references")

#: The fields a spec source and a config entry must agree on, if both declare them.
PINNED_SOURCE_FIELDS = ("url", "source_sha256")

#: Physical reference id -> asset suffix. The suffix is **source-prefixed**, which is not
#: the same string as the physical id (`hg19_ensembl` -> `ensembl-hg19`), so this mapping
#: is load-bearing: Task 10 commits an `asset_sha256` against each of these names, and a
#: rename after publication means a new release version.
GENOME_ASSETS: dict[str, str] = {
    "hg19": "ucsc-hg19",
    "hg38": "ucsc-hg38",
    "GRCh37": "ncbi-GRCh37",
    "GRCh38": "ncbi-GRCh38",
    "hg19_ensembl": "ensembl-hg19",
    "hg38_ensembl": "ensembl-hg38",
}

#: The MUC1 asset is not selectable - every installation needs it, whatever assembly it
#: asked for - so it has no physical id of its own.
MUC1_ASSET_SUFFIX = "muc1"

#: UCSC and Ensembl ship `.fa`, NCBI ships `.fna`. Exactly one must be present per genome.
FASTA_SUFFIXES = ("fa", "fna")

#: What `bwa index` writes beside a FASTA, and what `samtools faidx` writes.
BWA_SIDECARS = (".amb", ".ann", ".bwt", ".pac", ".sa")
GENOME_SIDECARS = (".fai", *BWA_SIDECARS)

#: The three MUC1 motif FASTAs plus both SHARK region cuts. Every one carries a `.fai`.
MUC1_FASTAS = (
    "All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa",
    "MUC1_motifs_Rev_com.fa",
    "code-adVNTR_RUs.fa",
    "muc1_region_hg19.fa",
    "muc1_region_hg38.fa",
)

#: Extracted out of `vntr_db_advntr.zip`; the zip itself is a build input, not an asset.
ADVNTR_DATABASES = ("vntr_db_advntr/hg19_muc1.db", "vntr_db_advntr/hg38_muc1.db")

MUC1_MEMBERS: tuple[str, ...] = (
    *(name for fasta in MUC1_FASTAS for name in (fasta, f"{fasta}.fai")),
    *ADVNTR_DATABASES,
)

#: Committed in the data repository rather than derived or downloaded.
SEED_FASTAS = ("MUC1_motifs_Rev_com.fa", "code-adVNTR_RUs.fa")

#: Every non-derivable artefact the data repository commits, and therefore everything the
#: release spec has to pin by digest: a build can reproduce anything else from an upstream
#: source, but these exist only because someone committed them.
REQUIRED_SEEDS = (*SEED_FASTAS, "vntr_db_advntr.zip", "filter_config.json")

#: The MUC1 FASTAs that are not seeds are, by construction, derivation outputs - so the
#: spec must declare a derivation for each. Derived from the frozen member set rather than
#: written out again, so the two cannot disagree.
DERIVED_FASTAS: tuple[str, ...] = tuple(name for name in MUC1_FASTAS if name not in SEED_FASTAS)

#: Written into every archive alongside the reference files, at the archive root.
METADATA_MEMBERS = ("release-manifest.json", "BUILD_INFO.json")

#: Exact relative paths that must never reach an asset, each with the reason it does not.
#: A file matching none of these and belonging to no asset fails the build.
EXCLUSIONS: tuple[tuple[str, str], ...] = (
    (
        "install_references.log",
        "installer log; setup_logging writes it into the output directory on every run",
    ),
    (
        "md5_checksums.txt",
        "build artefact of the legacy download path, superseded by release-manifest.json",
    ),
    (
        "filter_config.json",
        "build input for the merged-motif derivation; the pipeline never reads it",
    ),
    (
        "vntr_db_advntr.zip",
        "source archive for vntr_db_advntr/*.db, which ship extracted; shipping both duplicates them",
    ),
)

#: Glob-matched exclusions, same contract as EXCLUSIONS.
EXCLUSION_GLOBS: tuple[tuple[str, str], ...] = (
    (
        "alignment/*.gz",
        "downloaded source archive; the extracted FASTA ships instead, and keeping both "
        "roughly doubles every genome asset",
    ),
)

_BWA_VERSION = re.compile(r"^Version:\s*(\S+)", re.MULTILINE)
_SAMTOOLS_VERSION = re.compile(r"^samtools\s+(\S+)", re.MULTILINE)

#: gzip level 6 rather than 9. The bulk of a genome asset is the BWA `.bwt`/`.sa`/`.pac`,
#: which barely compress; level 9 buys a percent or two of size for a large multiple of
#: the CPU time, against a 240-minute job budget.
_COMPRESS_LEVEL = 6


@dataclass(frozen=True)
class FileEntry:
    """One file in the reference tree, already digested.

    Attributes:
        path: Path relative to the reference-tree root, POSIX-separated.
        size: Size in bytes.
        sha256: Lowercase hex SHA-256 of the file's contents.
    """

    path: str
    size: int
    sha256: str


@dataclass(frozen=True)
class Asset:
    """One release asset and the reference-tree files it carries.

    Attributes:
        name: Asset file name, e.g. `vntyper-references-refs-v1-ucsc-hg38.tar.gz`.
        reference_id: Physical reference id for a genome asset, None for the MUC1 asset.
        fasta: The genome's extracted chromosome FASTA, None for the MUC1 asset.
        members: Tree-relative paths carried by this asset.
    """

    name: str
    reference_id: str | None
    fasta: str | None
    members: tuple[str, ...]


def asset_name(tag: str, suffix: str) -> str:
    """Compose an asset file name.

    Args:
        tag: Release tag, e.g. `refs-v1`. Taken from the caller's spec, never hardcoded.
        suffix: Source-prefixed asset suffix, e.g. `ucsc-hg38`.

    Returns:
        str: The asset's file name.
    """
    return f"{ASSET_PREFIX}-{tag}-{suffix}.tar.gz"


def exclusion_reason(relative: str) -> str | None:
    """Say why a file is deliberately kept out of every asset.

    Args:
        relative: Tree-relative POSIX path.

    Returns:
        str | None: The documented reason, or None if the file is not excluded - in
        which case it must belong to an asset or the build fails.
    """
    for name, reason in EXCLUSIONS:
        if relative == name:
            return reason
    for pattern, reason in EXCLUSION_GLOBS:
        if fnmatch.fnmatch(relative, pattern):
            return reason
    return None


def genome_fasta(refs: Path, reference_id: str) -> str:
    """Locate the extracted chromosome FASTA of one assembly.

    Args:
        refs: Reference-tree root.
        reference_id: Physical reference id, e.g. `GRCh37`.

    Returns:
        str: Tree-relative path of the FASTA.

    Raises:
        ValueError: If neither spelling is present, or if both are - either way the
            build must stop rather than guess which bytes belong in the release.
    """
    candidates = [f"alignment/chr1.{reference_id}.{suffix}" for suffix in FASTA_SUFFIXES]
    present = [candidate for candidate in candidates if (refs / candidate).is_file()]
    if len(present) != 1:
        message = (
            f"reference '{reference_id}': expected exactly one extracted chromosome FASTA among "
            f"{', '.join(candidates)}; found {len(present)} ({', '.join(present) or 'none'})"
        )
        logger.error(message)
        raise ValueError(message)
    return present[0]


def genome_members(refs: Path, reference_id: str, fasta: str) -> tuple[str, ...]:
    """List everything one genome asset carries, failing on any missing index file.

    Args:
        refs: Reference-tree root.
        reference_id: Physical reference id, for the error message.
        fasta: Tree-relative path of the extracted FASTA.

    Returns:
        tuple[str, ...]: The FASTA, its `.fai` and the five BWA sidecars.

    Raises:
        ValueError: If any index file is absent. A genome asset without its BWA index is
            unusable, and the omission cannot be repaired after publication.
    """
    members = [fasta, *(f"{fasta}{suffix}" for suffix in GENOME_SIDECARS)]
    missing = [member for member in members if not (refs / member).is_file()]
    if missing:
        message = (
            f"reference '{reference_id}' is missing {len(missing)} index file(s): {', '.join(missing)}. "
            "Was the build run with --skip-indexing?"
        )
        logger.error(message)
        raise ValueError(message)
    return tuple(members)


def muc1_members(refs: Path) -> tuple[str, ...]:
    """List everything the MUC1 asset carries, failing on any missing file.

    Args:
        refs: Reference-tree root.

    Returns:
        tuple[str, ...]: Five MUC1 FASTAs with their `.fai`, and both adVNTR databases.

    Raises:
        ValueError: If any of them is absent.
    """
    missing = [member for member in MUC1_MEMBERS if not (refs / member).is_file()]
    if missing:
        message = (
            f"the MUC1 asset is missing {len(missing)} file(s): {', '.join(missing)}. "
            "Every installation depends on this asset regardless of assembly."
        )
        logger.error(message)
        raise ValueError(message)
    return MUC1_MEMBERS


def plan_assets(refs: Path, tag: str, references: Sequence[str]) -> list[Asset]:
    """Decide which file goes in which asset, before anything is written.

    Args:
        refs: Reference-tree root.
        tag: Release tag.
        references: Physical reference ids this build packages.

    Returns:
        list[Asset]: One asset per requested genome, then the implicit MUC1 asset.

    Raises:
        ValueError: On an unknown reference id, or on any missing expected file.
    """
    assets: list[Asset] = []
    for reference_id in references:
        suffix = GENOME_ASSETS.get(reference_id)
        if suffix is None:
            message = f"unknown reference id '{reference_id}'; the asset map covers {', '.join(sorted(GENOME_ASSETS))}"
            logger.error(message)
            raise ValueError(message)
        fasta = genome_fasta(refs, reference_id)
        assets.append(
            Asset(
                name=asset_name(tag, suffix),
                reference_id=reference_id,
                fasta=fasta,
                members=genome_members(refs, reference_id, fasta),
            )
        )
    assets.append(
        Asset(
            name=asset_name(tag, MUC1_ASSET_SUFFIX),
            reference_id=None,
            fasta=None,
            members=muc1_members(refs),
        )
    )
    return assets


def assign_every_file(present: Iterable[str], assets: Sequence[Asset]) -> dict[str, str]:
    """Prove the inventory is closed: every file is in one asset or explicitly excluded.

    This is the check that makes "we thought of the files we thought of" into a
    guarantee. It is deliberately bidirectional - an asset naming a file the tree does
    not have fails just as loudly as a file no asset claims.

    Args:
        present: Every tree-relative path in the reference tree.
        assets: The planned assets.

    Returns:
        dict[str, str]: Excluded path -> the reason it is excluded.

    Raises:
        ValueError: If a file is claimed by two assets, if an asset claims a file that is
            not in the tree, or if a file belongs to no asset and matches no exclusion.
            Every offending path is named.
    """
    owner: dict[str, str] = {}
    for asset in assets:
        for member in asset.members:
            if member in owner:
                message = f"'{member}' is claimed by both {owner[member]} and {asset.name}"
                logger.error(message)
                raise ValueError(message)
            owner[member] = asset.name

    paths = set(present)
    absent = sorted(set(owner) - paths)
    if absent:
        message = f"{len(absent)} planned asset member(s) are not in the reference tree: {', '.join(absent)}"
        logger.error(message)
        raise ValueError(message)

    excluded: dict[str, str] = {}
    unassigned: list[str] = []
    for path in sorted(paths - set(owner)):
        reason = exclusion_reason(path)
        if reason is None:
            unassigned.append(path)
        else:
            excluded[path] = reason

    if unassigned:
        message = (
            f"{len(unassigned)} file(s) under the reference tree belong to no asset and match no "
            f"exclusion: {', '.join(unassigned)}. Assign each one or add it to EXCLUSIONS with a "
            "reason; a release is not published with files nobody decided about."
        )
        logger.error(message)
        raise ValueError(message)
    return excluded


def inventory(refs: Path) -> dict[str, FileEntry]:
    """Digest every file in the reference tree once, up front.

    Args:
        refs: Reference-tree root.

    Returns:
        dict[str, FileEntry]: Tree-relative path -> its size and SHA-256.

    Raises:
        ValueError: If the tree contains a symbolic link. Nothing in a `--from-source`
            build creates one, and an archived link is a file whose recorded digest
            describes something other than what extraction produces.
    """
    entries: dict[str, FileEntry] = {}
    for path in sorted(refs.rglob("*")):
        relative = path.relative_to(refs).as_posix()
        if path.is_symlink():
            message = f"'{relative}' is a symbolic link; a reference bundle ships regular files only"
            logger.error(message)
            raise ValueError(message)
        if path.is_file():
            entries[relative] = FileEntry(relative, path.stat().st_size, sha256_of(path))
    return entries


def load_spec(spec_path: Path) -> dict[str, Any]:
    """Read the caller repository's committed release spec.

    Args:
        spec_path: Path to `releases/<tag>.json` in the data repository checkout.

    Returns:
        dict[str, Any]: The parsed document.

    Raises:
        ValueError: If the file is missing, unreadable, not JSON, or not a JSON object.
    """
    try:
        document = json.loads(spec_path.read_text(encoding="utf-8"))
    except OSError as error:
        message = f"cannot read release spec {spec_path}: {error}"
        logger.error(message)
        raise ValueError(message) from error
    except json.JSONDecodeError as error:
        message = f"release spec {spec_path} is not valid JSON: {error}"
        logger.error(message)
        raise ValueError(message) from error
    if not isinstance(document, dict):
        message = f"release spec {spec_path} must be a JSON object, got {type(document).__name__}"
        logger.error(message)
        raise ValueError(message)
    return document


def spec_seed_digests(spec: dict[str, Any]) -> dict[str, str]:
    """Read the spec's `seeds` block, accepting either spelling of an entry.

    Args:
        spec: A parsed release spec.

    Returns:
        dict[str, str]: Seed name -> expected SHA-256.

    Raises:
        ValueError: If the block is absent or empty, or an entry declares no digest.
    """
    seeds = spec.get("seeds")
    if not isinstance(seeds, dict) or not seeds:
        message = "release spec declares no non-empty 'seeds' block"
        logger.error(message)
        raise ValueError(message)

    digests: dict[str, str] = {}
    for name, entry in seeds.items():
        digest = entry.get("sha256") if isinstance(entry, dict) else entry
        if not isinstance(digest, str) or not digest:
            message = f"release spec seed '{name}' declares no sha256"
            logger.error(message)
            raise ValueError(message)
        digests[name] = digest
    return digests


def spec_source(spec: dict[str, Any], reference_id: str) -> dict[str, Any]:
    """Read one genome's pinned source, keyed exactly as the installer reads it.

    `install_references.resolve_source_location` looks the block up by **physical id**
    (`hg19_ensembl`, not `ensembl_hg19`) and reads `url` and `source_sha256`. The spec
    must use those names or the installer silently falls back to the digests committed in
    `install_references_config.json` and this build records provenance that does not
    describe what was downloaded.

    Args:
        spec: A parsed release spec.
        reference_id: Physical reference id.

    Returns:
        dict[str, Any]: That reference's source block, guaranteed to carry a non-empty
        `url` and `source_sha256`.

    Raises:
        ValueError: If the block is absent, is keyed by anything other than the physical
            id, declares no URL, or spells the digest field `sha256`.
    """
    sources = spec.get("sources")
    if not isinstance(sources, dict):
        message = "release spec declares no 'sources' block"
        logger.error(message)
        raise ValueError(message)
    entry = sources.get(reference_id)
    if not isinstance(entry, dict):
        message = (
            f"release spec names no source for reference '{reference_id}'. The 'sources' block is "
            "keyed by physical reference id, which is what install_references reads."
        )
        logger.error(message)
        raise ValueError(message)
    if not entry.get("url"):
        message = f"release spec source for '{reference_id}' declares no url"
        logger.error(message)
        raise ValueError(message)
    if not entry.get("source_sha256"):
        # The failure this rejects is silent, which is why it is worth a dedicated
        # message. `resolve_source_location` falls back field by field to
        # install_references_config.json, so a spec spelling the digest `sha256` still
        # builds - against the config's pin, not the spec's, while this manifest would
        # record the spec as the authority. Wrong bytes cannot ship (verify_tree is fatal
        # on a digest disagreement), but "the spec governs the build" would not hold.
        found = "'sha256'" if entry.get("sha256") else "no digest field"
        message = (
            f"release spec source for '{reference_id}' declares {found}; "
            "install_references.resolve_source_location reads 'source_sha256'. Rename the field: "
            "as written the spec's digest is ignored and the build silently falls back to the pin "
            "in install_references_config.json."
        )
        logger.error(message)
        raise ValueError(message)
    return entry


def committed_source_pins(config_path: Path | None = None) -> dict[str, dict[str, Any]]:
    """Read the URL and digest VNtyper commits for every downloadable genome.

    Args:
        config_path: Override for `COMMITTED_CONFIG`, for tests.

    Returns:
        dict[str, dict[str, Any]]: Physical reference id -> that reference's section of
        `install_references_config.json`.

    Raises:
        ValueError: If the file is missing, unreadable or not a JSON object. An
            unreadable trust anchor is a broken checkout, not a reason to skip the
            cross-check below and publish an unreviewed pin.
    """
    path = COMMITTED_CONFIG if config_path is None else config_path
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except OSError as error:
        message = f"cannot read VNtyper's committed reference config {path}: {error}"
        logger.error(message)
        raise ValueError(message) from error
    except json.JSONDecodeError as error:
        message = f"VNtyper's committed reference config {path} is not valid JSON: {error}"
        logger.error(message)
        raise ValueError(message) from error
    if not isinstance(document, dict):
        message = f"VNtyper's committed reference config {path} must be a JSON object"
        logger.error(message)
        raise ValueError(message)

    return {
        reference_id: entry
        for section in CONFIG_GENOME_SECTIONS
        for reference_id, entry in (document.get(section) or {}).items()
        if isinstance(entry, dict)
    }


def cross_check_sources(spec: dict[str, Any], references: Sequence[str], config_path: Path | None = None) -> None:
    """Refuse a spec whose sources contradict the pins committed in VNtyper.

    `install_references.resolve_source_location` merges the two field by field and
    `verify_tree` then checks the downloaded bytes against the **spec's** digest, so
    without this both layers agree with each other and neither ever consults
    `install_references_config.json`. A spec naming Ensembl release-115 while VNtyper
    documents release-116 would build, verify green, and freeze a provenance statement
    that does not describe the published bytes into an immutable release.

    §4.1 of the design says why the direction of authority is this way round: the trust
    anchor lives in VNtyper, not next to the assets. So moving a future release to a
    newer upstream requires a reviewed VNtyper commit - which is a base-image content
    hash change anyway.

    A reference the config does not pin keeps whatever the spec says; there is nothing
    to contradict.

    Args:
        spec: A parsed release spec.
        references: Physical reference ids this build packages.
        config_path: Override for `COMMITTED_CONFIG`, for tests.

    Raises:
        ValueError: On any disagreement, naming the reference, both values and the file
            each came from.
    """
    pins = committed_source_pins(config_path)
    anchor = COMMITTED_CONFIG if config_path is None else config_path

    disagreements: list[str] = []
    for reference_id in references:
        entry = spec_source(spec, reference_id)
        committed = pins.get(reference_id)
        if committed is None:
            continue
        for field in PINNED_SOURCE_FIELDS:
            declared = entry.get(field)
            pinned = committed.get(field)
            if not declared or not pinned or declared == pinned:
                continue
            disagreements.append(
                f"reference '{reference_id}' field '{field}': the release spec says {declared!r} "
                f"but {anchor.name} pins {pinned!r}"
            )

    if disagreements:
        message = (
            f"{len(disagreements)} source pin(s) in the release spec contradict the values committed in "
            f"{anchor}:\n  " + "\n  ".join(disagreements) + "\nThe trust anchor is VNtyper's committed "
            "config, not the spec beside the assets: a spec that silently wins would publish an immutable "
            "release whose provenance does not describe its bytes. Move the upstream in "
            "install_references_config.json first, in a reviewed commit, then match the spec to it."
        )
        logger.error(message)
        raise ValueError(message)


def validate_spec(spec: dict[str, Any], tag: str, references: Sequence[str]) -> None:
    """Reject a spec that cannot describe this build, before anything is downloaded.

    The workflow runs this as a preflight (`--check-spec-only`) in its staging step, so a
    malformed spec costs a minute instead of the three hours it takes to download and
    index six genomes.

    Args:
        spec: A parsed release spec.
        tag: The tag this build publishes under.
        references: Physical reference ids this build packages.

    Raises:
        ValueError: If the spec omits or contradicts the tag, names no usable source for
            a requested reference, contradicts a source pin committed in VNtyper, does
            not pin all four seeds, or does not fully declare every derivation. The spec
            is hand-written once and then frozen into an immutable release, so this is
            strict to the point of pedantry: a permissive preflight only moves the
            discovery of a typo from minute one to hour three, or past publication.
    """
    declared = spec.get("release_tag")
    if not declared:
        message = f"release spec declares no release_tag; this build publishes '{tag}'"
        logger.error(message)
        raise ValueError(message)
    if declared != tag:
        message = f"release spec declares tag '{declared}' but this build publishes '{tag}'"
        logger.error(message)
        raise ValueError(message)

    for reference_id in references:
        spec_source(spec, reference_id)
    cross_check_sources(spec, references)

    missing_seeds = [seed for seed in REQUIRED_SEEDS if seed not in spec_seed_digests(spec)]
    if missing_seeds:
        message = (
            f"release spec pins no digest for seed(s) {', '.join(missing_seeds)}. All "
            f"{len(REQUIRED_SEEDS)} committed seeds must be pinned - they are the only reference "
            "bytes a build cannot reproduce from an upstream source."
        )
        logger.error(message)
        raise ValueError(message)

    derivations = spec.get("derivations")
    if not isinstance(derivations, list) or not derivations:
        message = "release spec declares no non-empty 'derivations' list"
        logger.error(message)
        raise ValueError(message)
    for item in derivations:
        missing = (
            ["not a JSON object"]
            if not isinstance(item, dict)
            else [field for field in ("output", "command", "expected_sha256") if not item.get(field)]
        )
        if missing:
            # `command` is required, not optional: it is the only place the manifest's
            # `produced_by` can come from. Note the trap it guards - the derivations in
            # install_references_config.json carry `tool`, not `command`, so a spec copied
            # from the config's shape would leave every derived file's provenance blank.
            message = f"release spec derivation {item!r} is missing {', '.join(missing)}"
            logger.error(message)
            raise ValueError(message)

    undeclared = sorted(set(DERIVED_FASTAS) - {item["output"] for item in derivations})
    if undeclared:
        message = (
            f"release spec declares no derivation for {', '.join(undeclared)}. Every MUC1 FASTA "
            "that is not a committed seed is derived, and an undeclared one ships with neither a "
            "digest check nor a provenance record."
        )
        logger.error(message)
        raise ValueError(message)


def _capture(argv: list[str]) -> str:
    """Run a version probe and return its combined output, or "" if it cannot run.

    Args:
        argv: Command and arguments. `bwa` prints its banner to stderr and exits 1, so
            both streams are captured and the exit status is ignored.

    Returns:
        str: stdout and stderr concatenated, or "" when the binary is absent or hangs.
    """
    try:
        completed = subprocess.run(argv, capture_output=True, text=True, check=False, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return ""
    return f"{completed.stdout}\n{completed.stderr}"


def observed_toolchain() -> dict[str, str | None]:
    """Read the versions of the tools that produced the index and derivation files.

    Returns:
        dict[str, str | None]: `bwa_version` and `samtools_version`, None when the
        banner could not be read. Task 11's installer compares `bwa_version` against the
        local one and re-indexes with a warning when they differ.
    """
    bwa = _BWA_VERSION.search(_capture(["bwa"]))
    samtools = _SAMTOOLS_VERSION.search(_capture(["samtools", "--version"]))
    return {
        "bwa_version": bwa.group(1) if bwa else None,
        "samtools_version": samtools.group(1) if samtools else None,
    }


def _check(name: str, kind: str, expected: Any, actual: Any, ok: bool, verdict: str, *, fatal: bool = True) -> dict:
    return {
        "check": kind,
        "name": name,
        "expected": expected,
        "actual": actual,
        "ok": ok,
        "fatal": fatal,
        "verdict": verdict,
    }


def verify_tree(
    entries: dict[str, FileEntry],
    spec: dict[str, Any],
    assets: Sequence[Asset],
    observed: dict[str, str | None],
) -> list[dict[str, Any]]:
    """Re-check every digest the spec pins against the bytes actually in the tree.

    `install-references --from-source` already verified each download and each derivation
    against `install_references_config.json`. This is the independent second reading, and
    it is against the *spec* rather than the config, so a disagreement between the two
    surfaces here instead of being frozen into a release.

    Args:
        entries: The digested reference tree.
        spec: A parsed release spec.
        assets: The planned assets.
        observed: What `observed_toolchain` read from the environment.

    Returns:
        list[dict[str, Any]]: One record per check, each with `ok` and `fatal`. Nothing
        is raised - `main` writes the report before deciding, so a failed build still
        leaves evidence behind.
    """
    checks: list[dict[str, Any]] = []

    for name, expected in spec_seed_digests(spec).items():
        entry = entries.get(name)
        if entry is None:
            checks.append(
                _check(
                    name,
                    "seed",
                    expected,
                    None,
                    ok=True,
                    fatal=False,
                    verdict="not staged into the reference tree; verified in the seeds directory instead",
                )
            )
            continue
        matched = entry.sha256 == expected
        checks.append(
            _check(
                name,
                "seed",
                expected,
                entry.sha256,
                ok=matched,
                verdict="matches the release spec" if matched else "DOES NOT match the release spec",
            )
        )

    for item in spec.get("derivations", []):
        name = item["output"]
        expected = item["expected_sha256"]
        entry = entries.get(name)
        if entry is None:
            checks.append(
                _check(
                    name,
                    "derivation",
                    expected,
                    None,
                    ok=False,
                    verdict="declared by the release spec but absent from the reference tree",
                )
            )
            continue
        matched = entry.sha256 == expected
        checks.append(
            _check(
                name,
                "derivation",
                expected,
                entry.sha256,
                ok=matched,
                verdict="reproduces the digest the release spec pins"
                if matched
                else "DOES NOT reproduce the digest the release spec pins",
            )
        )

    for asset in assets:
        if asset.reference_id is None or asset.fasta is None:
            continue
        # `spec_source` guarantees a non-empty `source_sha256`; reading any other spelling
        # here is what let a spec using the wrong field name through the preflight.
        pinned: str = spec_source(spec, asset.reference_id)["source_sha256"]
        archive = f"{asset.fasta}.gz"
        entry = entries.get(archive)
        if entry is None:
            checks.append(
                _check(
                    archive,
                    "source-archive",
                    pinned,
                    None,
                    ok=True,
                    fatal=False,
                    verdict="download not retained in the reference tree; verified before decompression",
                )
            )
        else:
            matched = entry.sha256 == pinned
            checks.append(
                _check(
                    archive,
                    "source-archive",
                    pinned,
                    entry.sha256,
                    ok=matched,
                    verdict="matches the upstream source the release spec pins"
                    if matched
                    else "DOES NOT match the upstream source the release spec pins",
                )
            )

    for key in ("bwa_version", "samtools_version"):
        declared: str | None = spec.get(key)
        actual = observed[key]
        if not declared:
            continue
        matched = declared == actual
        checks.append(
            _check(
                key,
                "toolchain",
                declared,
                actual,
                ok=matched,
                fatal=False,
                verdict="matches the release spec"
                if matched
                else "differs from the release spec; every digest assertion still passed, "
                "but review the draft before publishing",
            )
        )

    assigned = sum(len(asset.members) for asset in assets)
    checks.append(
        _check(
            "reference tree",
            "inventory",
            len(entries),
            assigned,
            ok=True,
            verdict=f"{assigned} of {len(entries)} files assigned to {len(assets)} assets, "
            f"{len(entries) - assigned} explicitly excluded",
        )
    )
    return checks


def file_provenance(relative: str, asset: Asset, spec: dict[str, Any]) -> dict[str, Any]:
    """Describe how one file in an asset came to exist.

    Args:
        relative: Tree-relative path of the file.
        asset: The asset carrying it.
        spec: A parsed release spec, for source URLs and derivation commands.

    Returns:
        dict[str, Any]: Provenance fields to merge into the file's manifest entry.

    Raises:
        ValueError: If the file matches no provenance rule. `validate_spec` requires a
            declared derivation for every non-seed MUC1 FASTA, so with the frozen member
            set this cannot happen - reaching it means the member set and these rules
            have drifted apart, and a manifest entry with no provenance would ship.
    """
    if asset.reference_id is not None and asset.fasta is not None:
        if relative == asset.fasta:
            source = spec_source(spec, asset.reference_id)
            return {
                "source_url": source["url"],
                "source_sha256": source["source_sha256"],
                "produced_by": f"gunzip of {Path(asset.fasta).name}.gz",
            }
        if relative.endswith(".fai"):
            return {"produced_by": f"samtools faidx {Path(asset.fasta).name}"}
        return {"produced_by": f"bwa index {Path(asset.fasta).name}"}

    for item in spec.get("derivations", []):
        if item.get("output") == relative:
            return {"produced_by": item["command"], "expected_sha256": item["expected_sha256"]}
    if relative in SEED_FASTAS:
        return {"produced_by": "seed committed in the data repository"}
    if relative in ADVNTR_DATABASES:
        return {"produced_by": "extracted from vntr_db_advntr.zip"}
    if relative.endswith(".fai"):
        return {"produced_by": f"samtools faidx {relative[: -len('.fai')]}"}

    message = (
        f"no provenance rule covers '{relative}' in {asset.name}; refusing to record a manifest "
        "entry that does not say where the file came from"
    )
    logger.error(message)
    raise ValueError(message)


def build_info(tag: str, data_sha: str, builder_sha: str, observed: dict[str, str | None], now: str) -> dict[str, Any]:
    """Describe the machine and toolchain that produced this bundle.

    Args:
        tag: Release tag.
        data_sha: Resolved commit of the data repository checkout.
        builder_sha: Resolved commit of the VNtyper checkout this workflow was pinned at.
        observed: What `observed_toolchain` read.
        now: ISO-8601 UTC timestamp for this build.

    Returns:
        dict[str, Any]: The `BUILD_INFO.json` document. `bwa_version` and
        `samtools_version` are flat keys because the installer reads them by name.
    """
    return {
        "release_tag": tag,
        "generated_at": now,
        "builder_repository": "hassansaei/VNtyper",
        "builder_commit": builder_sha,
        "data_commit": data_sha,
        "bwa_version": observed["bwa_version"],
        "samtools_version": observed["samtools_version"],
        "python_version": platform.python_version(),
        # The runner image label itself is not read here: GitHub publishes it only as
        # `ImageOS`, and the repository's ruff configuration selects SIM112, which rejects
        # a mixed-case environment variable name. `runs-on` in the workflow and the
        # builder commit recorded above pin the image just as firmly.
        "runner": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
    }


def _normalise(info: tarfile.TarInfo) -> tarfile.TarInfo:
    """Strip build-host detail from an archive member.

    The uid, gid, user name and mtime of a file on a throwaway runner say nothing about
    the reference it holds, and leaving them in means two builds of identical bytes
    produce different archives. Removing them lets a rebuild be diffed against the
    published one down to the two metadata files.

    Args:
        info: The member header tarfile constructed.

    Returns:
        tarfile.TarInfo: The same header, normalised.
    """
    info.uid = 0
    info.gid = 0
    info.uname = ""
    info.gname = ""
    info.mtime = 0
    info.mode = 0o755 if info.isdir() else 0o644
    return info


def write_archive(refs: Path, archive: Path, members: Sequence[str], extras: dict[str, bytes]) -> None:
    """Write one asset, with tree-relative member paths and its metadata at the root.

    Member paths are relative to the reference-tree root, so extracting an asset into a
    target directory reproduces the layout directly - `alignment/chr1.hg38.fa`, never
    `refs/alignment/chr1.hg38.fa` and never an absolute path. `reference_bundle.safe_extract`
    rejects both of those, so a wrong root fails loudly at install time.

    Args:
        refs: Reference-tree root.
        archive: `.tar.gz` to create.
        members: Tree-relative paths to add, in a stable order.
        extras: Archive-root file name -> contents, added after the members.
    """
    archive.parent.mkdir(parents=True, exist_ok=True)
    with (
        archive.open("wb") as raw,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=_COMPRESS_LEVEL, mtime=0) as compressed,
        tarfile.open(fileobj=compressed, mode="w") as tar,
    ):
        for member in sorted(members):
            tar.add(refs / member, arcname=member, filter=_normalise)
        for name in sorted(extras):
            payload = extras[name]
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            tar.addfile(_normalise(info), io.BytesIO(payload))


def _json_bytes(document: dict[str, Any]) -> bytes:
    return (json.dumps(document, indent=2) + "\n").encode("utf-8")


def _write_json(path: Path, document: dict[str, Any]) -> None:
    path.write_bytes(_json_bytes(document))


def _prune(refs: Path, members: Iterable[str]) -> None:
    """Delete files whose archive has already been written and digested.

    Six indexed genomes are ~4.2 GB and their tarballs another ~2.7 GB, against roughly
    14 GB of free space on an `ubuntu-24.04` runner. Releasing each genome as soon as its
    asset exists keeps peak usage at the tree's own size instead of the sum of both.

    Args:
        refs: Reference-tree root.
        members: Tree-relative paths to remove; absent ones are ignored.
    """
    for member in members:
        (refs / member).unlink(missing_ok=True)


def build(args: argparse.Namespace) -> int:
    """Assemble, verify and describe the release.

    Args:
        args: Parsed command line.

    Returns:
        int: 0 on success, 1 if any fatal check failed.
    """
    now = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    refs: Path = args.refs
    out: Path = args.out

    spec = load_spec(args.spec)
    validate_spec(spec, args.tag, args.references)
    if args.check_spec_only:
        logger.info(f"Release spec {args.spec} describes tag '{args.tag}' and all {len(args.references)} references")
        return 0

    assets = plan_assets(refs, args.tag, args.references)
    entries = inventory(refs)
    excluded = assign_every_file(entries, assets)
    logger.info(
        f"{sum(len(asset.members) for asset in assets)} file(s) assigned to {len(assets)} asset(s), "
        f"{len(excluded)} explicitly excluded"
    )

    observed = observed_toolchain()
    checks = verify_tree(entries, spec, assets, observed)
    failures = [check for check in checks if not check["ok"] and check["fatal"]]

    out.mkdir(parents=True, exist_ok=True)
    info = build_info(args.tag, args.data_sha, args.builder_sha, observed, now)
    _write_json(
        out / "verification-report.json",
        {
            "release_tag": args.tag,
            "generated_at": now,
            "builder_commit": args.builder_sha,
            "data_commit": args.data_sha,
            "ok": not failures,
            "checks": checks,
        },
    )
    if failures:
        for failure in failures:
            logger.error(f"{failure['check']} {failure['name']}: {failure['verdict']}")
        logger.error(f"{len(failures)} fatal check(s) failed; no assets were written")
        return 1

    info_bytes = _json_bytes(info)
    # Every manifest is built before the first archive is written, so a provenance rule
    # that cannot classify a member aborts with nothing on disk - which matters under
    # `--prune`, where a mid-loop failure would have already deleted an earlier asset's
    # source files.
    planned: list[tuple[Asset, list[dict[str, Any]], dict[str, Any]]] = []
    for asset in assets:
        files = [
            {
                "path": member,
                "size": entries[member].size,
                "sha256": entries[member].sha256,
                **file_provenance(member, asset, spec),
            }
            for member in asset.members
        ]
        planned.append(
            (
                asset,
                files,
                {
                    "release_tag": args.tag,
                    "generated_at": now,
                    "asset": asset.name,
                    "reference_id": asset.reference_id,
                    "builder_commit": args.builder_sha,
                    "data_commit": args.data_sha,
                    "metadata": list(METADATA_MEMBERS),
                    "files": files,
                },
            )
        )

    records: list[dict[str, Any]] = []
    for asset, files, manifest in planned:
        archive = out / asset.name
        write_archive(
            refs,
            archive,
            asset.members,
            {"release-manifest.json": _json_bytes(manifest), "BUILD_INFO.json": info_bytes},
        )
        record = {
            "name": asset.name,
            "reference_id": asset.reference_id,
            "size": archive.stat().st_size,
            "sha256": sha256_of(archive),
            "files": files,
        }
        records.append(record)
        logger.info(f"{asset.name}: {len(files)} file(s), {record['size']} bytes, sha256 {record['sha256']}")
        if args.prune:
            _prune(refs, [*asset.members, *([f"{asset.fasta}.gz"] if asset.fasta else [])])

    _write_json(
        out / "release-manifest.json",
        {
            "release_tag": args.tag,
            "generated_at": now,
            "builder_repository": "hassansaei/VNtyper",
            "builder_commit": args.builder_sha,
            "data_commit": args.data_sha,
            "assets": records,
            "excluded": [
                {
                    "path": path,
                    "size": entries[path].size,
                    "sha256": entries[path].sha256,
                    "reason": reason,
                }
                for path, reason in sorted(excluded.items())
            ],
        },
    )
    _write_json(out / "BUILD_INFO.json", info)
    # Written here rather than by a shell `sha256sum *.tar.gz`: these are the same digests
    # the manifest records, so the two cannot disagree.
    (out / "SHA256SUMS").write_text(
        "".join(f"{record['sha256']}  {record['name']}\n" for record in records), encoding="utf-8"
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Parse arguments and run the build.

    Args:
        argv: Command-line arguments; `sys.argv[1:]` when None.

    Returns:
        int: 0 on success, 1 on any failure. Usage errors exit with argparse's 2.
    """
    parser = argparse.ArgumentParser(description="Assemble reference release assets from a built reference tree.")
    parser.add_argument("--refs", type=Path, required=True, help="Reference tree produced by --from-source.")
    parser.add_argument("--spec", type=Path, required=True, help="Committed release spec (releases/<tag>.json).")
    parser.add_argument("--tag", required=True, help="Release tag, e.g. refs-v1.")
    parser.add_argument("--out", type=Path, required=True, help="Directory to write the release assets into.")
    parser.add_argument("--data-sha", required=True, help="Resolved commit of the data repository checkout.")
    parser.add_argument("--builder-sha", required=True, help="Resolved commit of the VNtyper checkout.")
    parser.add_argument(
        "--references",
        nargs="+",
        default=list(GENOME_ASSETS),
        metavar="ID",
        help="Physical reference ids to package (default: all six).",
    )
    parser.add_argument(
        "--prune",
        action="store_true",
        help="Delete each asset's files once its archive is written, to cap peak disk usage.",
    )
    parser.add_argument(
        "--check-spec-only",
        action="store_true",
        help="Validate the release spec and exit, without reading the reference tree.",
    )
    args = parser.parse_args(argv)

    if not logging.getLogger().handlers:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    try:
        return build(args)
    except (ValueError, RuntimeError, OSError) as error:
        logger.error(f"bundle build failed: {error}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
