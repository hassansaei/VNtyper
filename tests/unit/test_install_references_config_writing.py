"""What `install-references --config-path` writes into config.json.

It used to write only ucsc_*/ncbi_*/ensembl_*/vntyper_*/own_repo_* keys, which nothing
in vntyper/ reads, while the seven keys the pipeline does read kept their shipped
relative paths. After an install into a custom --output-dir the run then died at
pipeline.py:154-156. Two of the old keys also named files that cannot be used: the
.gz rather than the extracted FASTA, and the .zip rather than the two .db files.
"""

import json

import pytest

from vntyper.scripts.install_references import canonical_reference_keys

pytestmark = pytest.mark.unit


class TestCanonicalKeys:
    def test_a_genome_entry_writes_the_key_the_pipeline_reads(self, tmp_path, install_config):
        keys = canonical_reference_keys(install_config, tmp_path)
        assert "bwa_reference_hg38" in keys
        assert "ucsc_hg38" not in keys

    def test_the_written_path_is_the_extracted_fasta_not_the_gz(self, tmp_path, install_config):
        keys = canonical_reference_keys(install_config, tmp_path)
        assert str(keys["bwa_reference_hg38"]).endswith("chr1.hg38.fa")

    def test_the_advntr_keys_name_the_databases_not_the_zip(self, tmp_path, install_config):
        keys = canonical_reference_keys(install_config, tmp_path)
        assert str(keys["advntr_reference_vntr_hg19"]).endswith("vntr_db_advntr/hg19_muc1.db")
        assert str(keys["advntr_reference_vntr_hg38"]).endswith("vntr_db_advntr/hg38_muc1.db")

    def test_both_shark_regions_are_written(self, tmp_path, install_config):
        keys = canonical_reference_keys(install_config, tmp_path)
        assert set(keys) >= {"muc1_region_fasta_hg19", "muc1_region_fasta_hg38"}

    def test_every_written_key_is_one_the_registry_knows(self, tmp_path, install_config):
        from vntyper.scripts.reference_registry import REFERENCE_KINDS, list_assemblies, reference_keys

        known = {k for kind in REFERENCE_KINDS for a in list_assemblies() for k in reference_keys(kind, a)}
        known |= {"muc1_reference_vntr", "code_adVNTR_RUs", "muc1_motifs_rev_com"}
        assert set(canonical_reference_keys(install_config, tmp_path)) <= known


class TestSchemaValidationRaisesOnAMissingField:
    """A malformed entry must fail loudly rather than silently drop its reference.

    `canonical_reference_keys` reads its required fields by `[]`, not `.get()`, on
    purpose: a config entry missing `installed_path`, `config_key`, `kind`, `assembly`
    or `output` is a malformed schema, and skipping it quietly would recreate the same
    class of bug #163 was - a reference the operator believes was installed silently
    has no key at all.
    """

    def test_a_genome_entry_without_an_installed_path_raises(self, tmp_path):
        config = {"ucsc_references": {"hg19": {"kind": "bwa"}}}
        with pytest.raises(KeyError):
            canonical_reference_keys(config, tmp_path)

    def test_a_common_reference_without_an_installed_path_raises(self, tmp_path):
        config = {"common_references": [{"config_key": "muc1_reference_vntr"}]}
        with pytest.raises(KeyError):
            canonical_reference_keys(config, tmp_path)

    def test_a_common_reference_without_a_config_key_raises(self, tmp_path):
        config = {"common_references": [{"installed_path": "x.fa"}]}
        with pytest.raises(KeyError):
            canonical_reference_keys(config, tmp_path)

    def test_a_derivation_without_a_kind_raises(self, tmp_path):
        config = {"derivations": [{"output": "x.fa", "config_key": "y"}]}
        with pytest.raises(KeyError):
            canonical_reference_keys(config, tmp_path)

    def test_a_shark_derivation_without_an_assembly_raises(self, tmp_path):
        config = {"derivations": [{"kind": "shark", "output": "x.fa"}]}
        with pytest.raises(KeyError):
            canonical_reference_keys(config, tmp_path)

    def test_a_literal_derivation_without_a_config_key_raises(self, tmp_path):
        config = {"derivations": [{"kind": "literal", "output": "x.fa"}]}
        with pytest.raises(KeyError):
            canonical_reference_keys(config, tmp_path)


class TestAnUnverifiedRetainedFileIsNotBlessed:
    """Parked finding (was issue #244): `canonical_reference_keys` used to write
    `config.json[key]` for *any* path that existed under `output_dir`, with no check
    that this or any earlier `install-references` run had ever verified it.

    `staged_install` seeds a new run's staging directory from whatever tree already
    exists there, precisely so a partial install (`--references hg19`) does not erase a
    different assembly a previous run installed. That seeding step also used to carry
    forward a file nobody had ever verified: an old, hand-copied or tampered
    `alignment/chr1.GRCh38.fna` sitting in the output directory. A
    `--reference-assembly GRCh38` run would then read that file as authoritative,
    entirely unverified.
    """

    def test_a_present_but_unverified_genome_file_is_not_written_as_a_key(self, tmp_path, install_config):
        """RED before the fix: the fixture's own hg38 entry is real and verified, but a
        second, GRCh38 entry is added here pointing at a file that exists on disk with
        no install_provenance.json record at all - exactly the "old file sitting in the
        output directory" scenario. Before item 1's fix this was blessed into the
        returned mapping identically to every verified entry; after it, it must be
        omitted.
        """
        install_config["ucsc_references"]["GRCh38_unverified"] = {"installed_path": "alignment/unverified.fa"}
        unverified = tmp_path / "alignment" / "unverified.fa"
        unverified.parent.mkdir(parents=True, exist_ok=True)
        unverified.write_text(">chr1\nACGT\n", encoding="utf-8")

        keys = canonical_reference_keys(install_config, tmp_path)

        assert "bwa_reference_GRCh38_unverified" not in keys, (
            "a file with no provenance record must never be written into config.json, even though it exists on disk"
        )
        # The fixture's real, provenance-backed entries are unaffected.
        assert "bwa_reference_hg38" in keys

    def test_the_same_file_is_written_once_a_provenance_record_exists(self, tmp_path, install_config):
        """Complement: once the file is actually recorded as verified (what a real
        install does), it resolves exactly like every other entry."""
        from vntyper.scripts.reference_provenance import build_record, merge

        install_config["ucsc_references"]["GRCh38_now_verified"] = {"installed_path": "alignment/verified.fa"}
        verified = tmp_path / "alignment" / "verified.fa"
        verified.parent.mkdir(parents=True, exist_ok=True)
        verified.write_text(">chr1\nACGT\n", encoding="utf-8")
        merge(
            tmp_path,
            {"alignment/verified.fa": build_record(sha256="c" * 64, source="from-source", source_url="https://x")},
        )

        keys = canonical_reference_keys(install_config, tmp_path)

        assert "bwa_reference_GRCh38_now_verified" in keys

    def test_the_warning_names_the_file_and_the_remedy(self, tmp_path, install_config, caplog):
        install_config["ucsc_references"]["hg19_unverified"] = {"installed_path": "alignment/stale.fa"}
        stale = tmp_path / "alignment" / "stale.fa"
        stale.parent.mkdir(parents=True, exist_ok=True)
        stale.write_text(">chr1\nACGT\n", encoding="utf-8")

        with caplog.at_level("WARNING"):
            canonical_reference_keys(install_config, tmp_path)

        message = " ".join(record.getMessage() for record in caplog.records)
        assert "stale.fa" in message
        assert "--references hg19_unverified" in message

    def test_a_common_reference_with_no_provenance_is_not_blessed(self, tmp_path, install_config):
        """The same gate applies to a common/adVNTR-style entry, not only a genome."""
        install_config["common_references"].append(
            {"config_key": "some_unverified_common_asset", "installed_path": "stray.db"}
        )
        (tmp_path / "stray.db").write_text("not actually installed by anything", encoding="utf-8")

        keys = canonical_reference_keys(install_config, tmp_path)

        assert "some_unverified_common_asset" not in keys


class TestUpdateConfigIsAtomic:
    def test_a_write_failure_leaves_the_previous_config_intact(self, tmp_path, monkeypatch):
        from vntyper.scripts import install_references

        config_path = tmp_path / "config.json"
        config_path.write_text(json.dumps({"reference_data": {"bwa_reference_hg19": "old"}}))
        monkeypatch.setattr(
            install_references.json, "dump", lambda *a, **k: (_ for _ in ()).throw(OSError("disk full"))
        )
        with pytest.raises((OSError, SystemExit)):
            install_references.update_config(config_path, {"bwa_reference_hg19": "new"})
        assert json.loads(config_path.read_text())["reference_data"]["bwa_reference_hg19"] == "old"
