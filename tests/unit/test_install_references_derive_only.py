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

import logging

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

    def record(config, out, refs, selected):
        seen.update(refs=refs, selected=selected)
        return []

    monkeypatch.setattr(install_references, "_preflight_literal_seeds", lambda *a, **k: None)
    monkeypatch.setattr(install_references, "run_derivations", record)

    install_references.derive_only(_config(), tmp_path)

    assert seen["selected"] == {"hg19"}
    assert set(seen["refs"]) == {"hg19"}


def test_derive_only_preflights_the_literal_seeds_before_deriving(tmp_path, monkeypatch):
    """Missing seeds must be reported as missing seeds, not as a checksum mismatch on an
    output built from nothing."""
    order: list[str] = []

    def derive(*args, **kwargs):
        order.append("derive")
        return []

    monkeypatch.setattr(install_references, "_preflight_literal_seeds", lambda *a, **k: order.append("preflight"))
    monkeypatch.setattr(install_references, "run_derivations", derive)

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
    monkeypatch.setattr(install_references, "run_derivations", lambda *a, **k: [])

    install_references.derive_only(_config(), tmp_path)


def _config_with_two_shark_derivations() -> dict:
    config = _config()
    config["derivations"] = [
        {"output": "muc1_region_hg19.fa", "kind": "shark", "from": "hg19"},
        {"output": "muc1_region_hg38.fa", "kind": "shark", "from": "hg38"},
    ]
    return config


def test_a_skipped_derivation_is_not_reported_as_verified(tmp_path, monkeypatch, caplog):
    """The summary must not claim a digest match for a derivation that was skipped.

    A skip writes nothing, so no unverified file is *produced*. But a file may still sit at
    that path from an earlier install, and this run never read it. Reporting a blanket
    success would assert that those bytes match their committed digest on the strength of
    nothing -- the same defect as inferring samtools' exit status from a file's existence
    (#255). What is skipped has to be named.
    """
    (tmp_path / "alignment").mkdir()
    (tmp_path / "alignment" / "chr1.hg19.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
    # A stale leftover from an earlier install, whose source genome is no longer present.
    (tmp_path / "muc1_region_hg38.fa").write_text(">stale\nAAAA\n", encoding="utf-8")

    monkeypatch.setattr(install_references, "_preflight_literal_seeds", lambda *a, **k: None)
    monkeypatch.setattr(install_references, "run_derivations", lambda *a, **k: ["muc1_region_hg19.fa"])

    with caplog.at_level(logging.INFO):
        install_references.derive_only(_config_with_two_shark_derivations(), tmp_path)

    assert "not verified: muc1_region_hg38.fa" in caplog.text
    assert "verified all" not in caplog.text, "a run that skipped a derivation must not claim all of them"


def test_deriving_every_configured_file_does_report_all_of_them(tmp_path, monkeypatch, caplog):
    """The counterpart, so the warning above is a real branch rather than the only outcome.

    Without this, deleting the success message entirely would leave the skip test green.
    """
    monkeypatch.setattr(install_references, "_preflight_literal_seeds", lambda *a, **k: None)
    monkeypatch.setattr(
        install_references,
        "run_derivations",
        lambda *a, **k: ["muc1_region_hg19.fa", "muc1_region_hg38.fa"],
    )

    with caplog.at_level(logging.INFO):
        install_references.derive_only(_config_with_two_shark_derivations(), tmp_path)

    assert "verified all 2 reference file(s)" in caplog.text
    assert "not verified" not in caplog.text


def test_run_derivations_reports_only_what_it_verified(tmp_path):
    """The return value that the summary rests on, checked against the real function.

    A skipped ``shark`` derivation must be absent from it. If ``run_derivations`` appended
    before ``verify_sha256`` instead of after, or appended on the skip path, the honest
    summary above would go back to being a blanket success with extra steps.
    """
    config = {
        "derivations": [
            {
                "output": "muc1_region_hg38.fa",
                "kind": "shark",
                "from": "hg38",
                "region": "chr1:1-10",
                "expected_sha256": "a" * 64,
            }
        ]
    }

    assert install_references.run_derivations(config, tmp_path, {}, selected={"hg19"}) == []
    assert not (tmp_path / "muc1_region_hg38.fa").exists(), "a skip must not write anything"


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


def test_combining_derive_only_with_from_source_is_a_usage_error(tmp_path):
    """Refused rather than silently ignored. Both build the derived files, but
    --derive-only deliberately downloads nothing, so honouring both would make one of them
    meaningless. Usage errors exit 2, as argparse does."""
    import argparse

    from vntyper.scripts import cli_handlers

    parser = argparse.ArgumentParser()
    args = argparse.Namespace(
        output_dir=tmp_path,
        config_path=None,
        skip_indexing=False,
        threads=1,
        aligners=None,
        references=None,
        derive_only=True,
        from_source=True,
        release_spec=None,
    )

    with pytest.raises(SystemExit) as excinfo:
        cli_handlers.handle_install_references(args, {}, parser, 0, None)

    assert excinfo.value.code == 2


def test_main_takes_the_derive_only_path_and_installs_nothing(tmp_path, monkeypatch):
    """`main` must return after deriving: reaching an installer would download."""
    called: list[str] = []
    monkeypatch.setattr(install_references, "derive_only", lambda cfg, out: called.append("derived"))
    monkeypatch.setattr(
        install_references,
        "install_from_bundle",
        lambda *a, **k: called.append("downloaded"),
    )

    install_references.main(output_dir=tmp_path, derive_only_mode=True)

    assert called == ["derived"], "derive-only must not reach an installer"
