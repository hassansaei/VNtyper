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

import hashlib
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


def _bundle_tree_config() -> dict:
    """The shape a bundle-installed tree actually has, which is what the image ships.

    The published bundle ships the merged motif FASTA pre-built and does **not** stage
    ``filter_config.json`` beside it, so the literal derivation has no seeds and never will.
    """
    return {
        "ucsc_references": {"hg19": {"installed_path": "alignment/chr1.hg19.fa"}},
        "own_repository_references": {"raw_files": [{"target_path": "filter_config.json", "source_sha256": "a" * 64}]},
        "derivations": [
            {
                "output": "All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa",
                "kind": "literal",
                "from_seeds": ["MUC1_motifs_Rev_com.fa", "filter_config.json"],
                "expected_sha256": GOOD_DIGEST,
            }
        ],
    }


def test_a_literal_derivation_without_its_seeds_is_skipped_not_failed(tmp_path, caplog):
    """The regression a run inside the Docker image found.

    ``run_derivations`` raises for a literal derivation whose seeds are absent, which is
    right for ``--from-source``: it stages the seeds itself, so a missing one is a real
    fault. ``--derive-only`` stages nothing, so on a bundle-installed tree -- the tree the
    image carries, and the commonest in existence -- an absent ``filter_config.json`` is
    simply that tree's shape. Failing there made the command unusable exactly where it would
    most often be reached for.

    It is skipped like a shark derivation whose genome is absent, and the file already at
    that path is verified, so nothing is taken on trust.
    """
    (tmp_path / "MUC1_motifs_Rev_com.fa").write_text(">seed\nAC\n", encoding="utf-8")
    (tmp_path / "All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa").write_text(GOOD_BYTES, encoding="utf-8")

    with caplog.at_level(logging.INFO):
        install_references.derive_only(_bundle_tree_config(), tmp_path)

    assert "seed(s) filter_config.json are not in this tree" in caplog.text
    assert "already present and matching their committed digests" in caplog.text


def test_a_skipped_literal_derivation_still_has_its_output_verified(tmp_path):
    """Skipping must not become a way to carry a wrong reference forward.

    The seeds being absent is why it cannot be rebuilt; it is not a reason to stop looking
    at what is there. A bundle-installed tree whose motif FASTA was corrupted would
    otherwise pass this command silently.
    """
    stale = tmp_path / "All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa"
    stale.write_text(">stale\nTTTT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa"):
        install_references.derive_only(_bundle_tree_config(), tmp_path)

    assert not stale.exists()


def test_a_literal_derivation_with_its_seeds_present_is_still_built(tmp_path, monkeypatch):
    """The counterpart, so the skip above is a branch rather than the only behaviour.

    Without this, filtering out *every* literal derivation would leave both tests green.
    """
    config = _bundle_tree_config()
    (tmp_path / "MUC1_motifs_Rev_com.fa").write_text(">seed\nAC\n", encoding="utf-8")
    (tmp_path / "filter_config.json").write_text("{}", encoding="utf-8")
    handed: list[list] = []

    def record(cfg, out, refs, selected):
        handed.append(cfg.get("derivations", []))
        return [spec["output"] for spec in cfg.get("derivations", [])]

    monkeypatch.setattr(install_references, "run_derivations", record)

    install_references.derive_only(config, tmp_path)

    assert handed == [config["derivations"]], "a literal derivation whose seeds are present must be built"


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


#: sha256 of ">good\nACGT\n", so a fixture can be made to match or mismatch on purpose.
GOOD_BYTES = ">good\nACGT\n"
GOOD_DIGEST = hashlib.sha256(GOOD_BYTES.encode("utf-8")).hexdigest()


def _config_with_two_shark_derivations(hg38_digest: str = GOOD_DIGEST) -> dict:
    config = _config()
    config["derivations"] = [
        {"output": "muc1_region_hg19.fa", "kind": "shark", "from": "hg19", "expected_sha256": GOOD_DIGEST},
        {"output": "muc1_region_hg38.fa", "kind": "shark", "from": "hg38", "expected_sha256": hg38_digest},
    ]
    return config


def test_a_skipped_derivation_is_not_reported_as_verified(tmp_path, monkeypatch, caplog):
    """The summary must not claim a digest match for a derivation that was skipped.

    A skip writes nothing, so no unverified file is *produced*. But reporting a blanket
    success would assert a digest match on the strength of nothing -- the same defect as
    inferring samtools' exit status from a file's existence (#255). What was not rebuilt
    has to be named as not rebuilt.
    """
    (tmp_path / "alignment").mkdir()
    (tmp_path / "alignment" / "chr1.hg19.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "muc1_region_hg38.fa").write_text(GOOD_BYTES, encoding="utf-8")

    monkeypatch.setattr(install_references, "_preflight_literal_seeds", lambda *a, **k: None)
    monkeypatch.setattr(install_references, "run_derivations", lambda *a, **k: ["muc1_region_hg19.fa"])

    with caplog.at_level(logging.INFO):
        install_references.derive_only(_config_with_two_shark_derivations(), tmp_path)

    assert "Not rebuilt" in caplog.text and "muc1_region_hg38.fa" in caplog.text
    assert "verified all" not in caplog.text, "a run that skipped a derivation must not claim all of them"


def test_a_file_that_cannot_be_rebuilt_is_still_checked_against_its_digest(tmp_path, monkeypatch, caplog):
    """The command's whole purpose is answering "are my derived files right?".

    Answering with silence for exactly the files it could not rebuild, while exiting 0, is
    the weakest useful thing it could do. The digest is committed and the file is small.
    """
    (tmp_path / "muc1_region_hg38.fa").write_text(GOOD_BYTES, encoding="utf-8")
    monkeypatch.setattr(install_references, "_preflight_literal_seeds", lambda *a, **k: None)
    monkeypatch.setattr(install_references, "run_derivations", lambda *a, **k: ["muc1_region_hg19.fa"])

    with caplog.at_level(logging.INFO):
        install_references.derive_only(_config_with_two_shark_derivations(), tmp_path)

    assert "already present and matching their committed digests: muc1_region_hg38.fa" in caplog.text


def test_a_stale_file_that_cannot_be_rebuilt_is_discarded_rather_than_left(tmp_path, monkeypatch):
    """The case the check exists for, and the one that reaches a genotyping run if missed.

    A tree that lost its hg38 genome but kept a `muc1_region_hg38.fa` from an older install
    would otherwise carry that file forward silently. A wrong reference produces a plausible
    result rather than an obvious failure, so it is discarded, exactly as `run_derivations`
    discards a derived output that fails its digest.
    """
    stale = tmp_path / "muc1_region_hg38.fa"
    stale.write_text(">stale\nAAAA\n", encoding="utf-8")
    monkeypatch.setattr(install_references, "_preflight_literal_seeds", lambda *a, **k: None)
    monkeypatch.setattr(install_references, "run_derivations", lambda *a, **k: ["muc1_region_hg19.fa"])

    with pytest.raises(ValueError, match="muc1_region_hg38.fa"):
        install_references.derive_only(_config_with_two_shark_derivations(), tmp_path)

    assert not stale.exists(), "a reference that fails its digest must not be left in the tree"


def test_a_file_that_cannot_be_rebuilt_and_is_absent_is_reported_as_missing(tmp_path, monkeypatch, caplog):
    """Absent is not a digest failure -- there is nothing wrong in the tree, only nothing
    there. The message has to name the genome that would let it be rebuilt."""
    monkeypatch.setattr(install_references, "_preflight_literal_seeds", lambda *a, **k: None)
    monkeypatch.setattr(install_references, "run_derivations", lambda *a, **k: ["muc1_region_hg19.fa"])

    with caplog.at_level(logging.INFO):
        install_references.derive_only(_config_with_two_shark_derivations(), tmp_path)

    assert "missing from the tree: muc1_region_hg38.fa" in caplog.text


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


def test_the_parser_accepts_both_flags_and_leaves_the_refusal_to_the_handler():
    """Named for what it checks. The refusal itself is
    ``test_combining_derive_only_with_from_source_is_a_usage_error`` below; this only pins
    where that refusal lives, so the two are not both assumed to be somebody else's job."""
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
