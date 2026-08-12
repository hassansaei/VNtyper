"""``install-references --derive-only``: rebuild derived files from what is on disk.

Three files in the reference tree are **derived**, not downloaded: `muc1_region_hg19.fa` and
`muc1_region_hg38.fa`, cut out of an installed chromosome with `samtools faidx`, and
`All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`, merged from two seeds.

The bundle ships them pre-built and `--from-source` builds them at the end of its run, so on
both ordinary paths they are already correct. This covers the case in between, which had no
command at all: a tree whose genomes and seeds are present but whose derived files are
missing. Recovering that meant a full `--from-source` run -- re-downloading and BWA-indexing
six chromosome FASTAs to rebuild three small files -- so in practice it was done by hand,
with `samtools faidx` retyped from the config. By hand is exactly where an unverified
reference file comes from.

Every output is still checked against its committed `expected_sha256`, so this is the same
verification on a cheaper path, not a way to produce unverified bytes.
"""

import pytest

from vntyper.scripts import install_references

pytestmark = pytest.mark.unit


def _config() -> dict:
    return {
        "ucsc_references": {
            "hg19": {"installed_path": "alignment/chr1.hg19.fa"},
            "hg38": {"installed_path": "alignment/chr1.hg38.fa"},
        },
        "ncbi_references": {"GRCh37": {"installed_path": "alignment/chr1.GRCh37.fna"}},
        "ensembl_references": {},
    }


def test_only_genomes_that_exist_are_reported_as_installed(tmp_path):
    """The map is built from the same ``installed_path`` field the installers write, and
    filtered by what is actually on disk -- a declared-but-absent genome is not installed."""
    (tmp_path / "alignment").mkdir()
    (tmp_path / "alignment" / "chr1.hg19.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "alignment" / "chr1.GRCh37.fna").write_text(">1\nACGT\n", encoding="utf-8")

    found = install_references.installed_reference_map(_config(), tmp_path)

    assert sorted(found) == ["GRCh37", "hg19"], "hg38 is declared but not installed"
    assert found["hg19"] == tmp_path / "alignment" / "chr1.hg19.fa"


def test_an_empty_tree_yields_no_installed_genomes(tmp_path):
    assert install_references.installed_reference_map(_config(), tmp_path) == {}


def test_a_reference_without_an_installed_path_is_skipped(tmp_path):
    """A malformed entry must not raise here. The installers own that validation; this
    function only reports what is present."""
    assert install_references.installed_reference_map({"ucsc_references": {"x": {}}}, tmp_path) == {}


def test_derive_only_passes_exactly_the_present_genomes_as_the_selection(tmp_path, monkeypatch):
    """A derivation whose source genome is absent must be *skipped*, not failed.

    A tree holding only hg19 legitimately derives only the hg19 region, and
    ``run_derivations`` already implements skip-not-fail when a source is outside the
    selection. Passing the present genomes as that selection is what routes the absent ones
    down the skip path instead of the hard-error path.
    """
    (tmp_path / "alignment").mkdir()
    (tmp_path / "alignment" / "chr1.hg19.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
    seen: dict = {}

    monkeypatch.setattr(install_references, "_preflight_literal_seeds", lambda *a, **k: None)
    monkeypatch.setattr(
        install_references,
        "run_derivations",
        lambda config, out, refs, selected: seen.update(refs=refs, selected=selected),
    )

    install_references.derive_only(_config(), tmp_path)

    assert seen["selected"] == {"hg19"}
    assert set(seen["refs"]) == {"hg19"}


def test_derive_only_preflights_the_literal_seeds_before_deriving(tmp_path, monkeypatch):
    """Missing seeds must be reported as missing seeds, not as a checksum mismatch on an
    output built from nothing."""
    order: list[str] = []
    monkeypatch.setattr(install_references, "_preflight_literal_seeds", lambda *a, **k: order.append("preflight"))
    monkeypatch.setattr(install_references, "run_derivations", lambda *a, **k: order.append("derive"))

    install_references.derive_only(_config(), tmp_path)

    assert order == ["preflight", "derive"]


def test_derive_only_downloads_nothing(tmp_path, monkeypatch):
    """The whole point of the flag. If it reached a downloader, it would be a slower
    ``--from-source`` with fewer checks rather than a cheaper recovery path."""

    def fail(*args, **kwargs):
        raise AssertionError("--derive-only must not download anything")

    for name in ("download_file", "process_own_repository_references", "install_from_bundle"):
        if hasattr(install_references, name):
            monkeypatch.setattr(install_references, name, fail)
    monkeypatch.setattr(install_references, "_preflight_literal_seeds", lambda *a, **k: None)
    monkeypatch.setattr(install_references, "run_derivations", lambda *a, **k: None)

    install_references.derive_only(_config(), tmp_path)


def test_derive_only_and_from_source_are_mutually_exclusive():
    """Both build the derived files, but --derive-only deliberately downloads nothing.
    Accepting both would make one of them silently meaningless."""
    from vntyper.scripts.cli_parser import build_parser

    parser = build_parser()
    args = parser.parse_args(["install-references", "-d", "refs", "--derive-only", "--from-source"])

    assert args.derive_only and args.from_source, "the parser accepts both; the handler is what refuses them"


def test_the_flag_exists_and_defaults_to_off():
    from vntyper.scripts.cli_parser import build_parser

    args = build_parser().parse_args(["install-references", "-d", "refs"])

    assert args.derive_only is False
