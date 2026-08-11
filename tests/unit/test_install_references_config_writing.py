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
