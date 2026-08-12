"""Where reference data may come from, as a config-shape invariant (#253).

**Every reference file's source of truth is `berntpopp/vntyper-data`, never the code
repository.** That is the milestone-5 design: reference data comes from a published,
checksummed bundle rather than from tracked files, and `reference/` keeps only
`README.md`, `pseudonymize.py` and `pseudonymize_config.json`.

One entry outlived that change. `own_repository_references.raw_files` carried
`All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa` pointed at
`raw.githubusercontent.com/hassansaei/VNtyper/main/reference/...`, which 404s because
milestone 5 emptied that directory.

It was harmless only by accident. `_install_source_seeds` drops any `raw_files` entry whose
`target_path` is also a `derivations` output, and that file is one -- so the URL was never
requested. It was also the only entry with no `source_sha256`, and `--from-source` refuses an
undigested seed, so even had it been reached the run would have aborted rather than 404'd.

**A 404 was the mildest way that could have failed.** A URL that resolved to a *stale but
live* file would have installed silently wrong reference data, with nothing to signal it.
These tests pin the rule rather than the one entry.
"""

import json
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

CONFIG_PATH = Path(__file__).resolve().parents[2] / "vntyper" / "scripts" / "install_references_config.json"

#: Hosts a reference file may legitimately be fetched from: this project's data repository,
#: and the three upstream genome providers.
ALLOWED_SOURCES = (
    "raw.githubusercontent.com/berntpopp/vntyper-data/",
    "https://github.com/berntpopp/vntyper-data/",
    "hgdownload.soe.ucsc.edu/",
    "ftp.ncbi.nlm.nih.gov/",
    "ftp.ensembl.org/",
)

#: The code repository. Nothing under it may be a reference-data source.
CODE_REPOSITORY = "hassansaei/VNtyper"


def _config() -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def _every_url(node, found=None) -> list[str]:
    """Collect every ``url`` value anywhere in the config, at any nesting depth."""
    found = [] if found is None else found
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "url" and isinstance(value, str):
                found.append(value)
            else:
                _every_url(value, found)
    elif isinstance(node, list):
        for item in node:
            _every_url(item, found)
    return found


def test_no_reference_file_is_sourced_from_the_code_repository():
    """The rule, checked over every URL at any depth rather than the one entry that broke it.

    A reference file fetched from the repository that consumes it is what milestone 5
    removed. Re-introducing one would not necessarily 404 -- it might resolve to a stale copy
    and install quietly wrong data.
    """
    offenders = [url for url in _every_url(_config()) if CODE_REPOSITORY in url]

    assert offenders == [], f"reference data must come from berntpopp/vntyper-data, not the code repo: {offenders}"


def test_every_source_url_names_an_allowed_host():
    """Positive form of the same rule, so a *new* unexpected host is also caught."""
    unexpected = [url for url in _every_url(_config()) if not any(allowed in url for allowed in ALLOWED_SOURCES)]

    assert unexpected == [], f"unrecognised reference source host: {unexpected}"


def test_no_raw_file_entry_is_also_a_derivation_output():
    """The dead-entry shape: config that looks like a fetch target but is never fetched.

    ``_install_source_seeds`` filters these out, so such an entry is unreachable by
    construction -- its URL and digest are never exercised and can rot unnoticed. That is
    exactly how a 404 survived in the tree.
    """
    config = _config()
    derived = {spec.get("output") for spec in config.get("derivations", [])}
    raw_files = config.get("own_repository_references", {}).get("raw_files", [])

    both = sorted({entry["target_path"] for entry in raw_files} & derived)

    assert both == [], f"these are declared as both a downloadable seed and a derivation output: {both}"


def test_every_raw_file_entry_carries_a_source_digest():
    """The property whose absence marked the dead entry as different from its neighbours.

    ``install_from_source`` already refuses a seed with no ``source_sha256`` at runtime; this
    moves that from a refusal that fires during a multi-gigabyte install to a guarantee about
    the config's shape.
    """
    raw_files = _config().get("own_repository_references", {}).get("raw_files", [])
    undigested = [entry["target_path"] for entry in raw_files if not entry.get("source_sha256")]

    assert undigested == [], f"a seed with no source_sha256 cannot be verified: {undigested}"


def test_the_motif_fasta_is_derived_rather_than_fetched():
    """The specific regression, pinned by name.

    ``All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`` is built from
    ``MUC1_motifs_Rev_com.fa`` + ``filter_config.json`` and shipped in the published bundle.
    It must not reacquire a download entry -- two sources for the same bytes means only one
    of them is digest-checked.
    """
    config = _config()
    name = "All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa"
    raw_files = config.get("own_repository_references", {}).get("raw_files", [])

    assert name not in {entry["target_path"] for entry in raw_files}
    assert name in {spec.get("output") for spec in config.get("derivations", [])}
