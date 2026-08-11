"""The MUC1 region FASTAs are derived, not shipped - and the derivation is verified.

`samtools faidx chr1.hg19.fa chr1:155158000-155163000` reproduces the FASTA this
repository used to track, byte for byte (md5 c9737129069d4855b433b178ebb21e1c). That
is why the file could be removed from git: the bundle build asserts the derivation
against a digest committed in the release spec, so a silent change in a UCSC
chromosome file turns the build red instead of publishing different sequence under an
unchanged name.

The same argument covers `All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`, which
is a pure function of two seeds; the merge is ported into `install_references.py` so
that `--from-source` - the path the bundle build workflow itself runs - produces it
without a helper script beside it.
"""

import gzip
import hashlib
import json
import logging
import shlex
import subprocess
from pathlib import Path

import pytest

from vntyper.scripts import install_references
from vntyper.scripts.cli_parser import build_parser
from vntyper.scripts.samtools_command_fragments import quote_path

pytestmark = pytest.mark.unit


#: The digests the controller measured against the tracked reference tree. A test that
#: disagrees with one of these is reporting a real regression, not a stale constant.
HG19_REGION_SHA256 = "d6007957193f9686e2a8b02becb8a7d02a8aa5dae4f087d233e24e4145b88a1d"
HG38_REGION_SHA256 = "60026035da804b4a7234ae051e4f8440c4a7f05f517affc47684ca85b5c25a6a"
MERGED_MOTIFS_SHA256 = "bd6a895ed9fee6aafdd2742f7d295f0b5b26d1a7bb613b9de7793eeba02df441"

SOURCE_SHA256 = {
    "hg19": "8d99e436ac1b8ee27999917cc4384a0bfec5a75467b59bbbea8ad314adb553d7",
    "hg38": "baf03990d81ed87d66a976779c2eb4d8a0e5243829ca913556d3819fe0fba696",
    "GRCh37": "cb58b0a7eb8d478d6f3b803182d6de7654e621a42223ba90c33d02c592dbdfb9",
    "GRCh38": "191c0c9f25748e67eb251a5cbeb2289e10dadc98746b4ddf0170c559db219a2e",
    "hg19_ensembl": "2f13339ef278b05fd7d80d60052ba12dbd661b6396024b2043a2c5ada07aab96",
    "hg38_ensembl": "05d4d42ed292962055afb9774a1a29a691d7b285fd85dbf9069c050a273e0d3c",
}


@pytest.fixture(autouse=True)
def _restore_root_logging():
    """`setup_logging` re-points the root logger at a file; put it back afterwards.

    `install_references.main` calls it, so without this every later test in the session
    would keep writing into a deleted `tmp_path`.
    """
    root = logging.getLogger()
    handlers = root.handlers[:]
    level = root.level
    yield
    for handler in root.handlers[:]:
        root.removeHandler(handler)
        if handler not in handlers:
            handler.close()
    for handler in handlers:
        root.addHandler(handler)
    root.setLevel(level)


def _shipped_config() -> dict:
    """Return the parsed `install_references_config.json` that actually ships."""
    path = Path(install_references.__file__).parent / "install_references_config.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _fake_samtools(store: dict[str, bytes], calls: list[str]):
    """Return a `subprocess.run` stand-in that implements just enough of samtools.

    Two forms are honoured, which are the only two the module builds:
    `samtools faidx SRC REGION > DEST` writes `store[REGION]` to DEST, and
    `samtools faidx FILE` writes a `.fai` beside FILE.

    Args:
        store: Region string to the bytes a faidx of that region should produce.
        calls: List the raw command strings are appended to, for assertions.

    Returns:
        Callable: A drop-in for `subprocess.run` with `shell=True`.
    """

    def _run(command, **kwargs):
        calls.append(command)
        tokens = shlex.split(command)
        if ">" in tokens:
            region, destination = tokens[-3], tokens[-1]
            Path(destination).write_bytes(store[region])
        else:
            Path(tokens[-1] + ".fai").write_text("stub\n", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    return _run


# =============================================================================
# _quote
# =============================================================================


class TestQuote:
    """The inlined quoting helper must not drift from the one it copies."""

    @pytest.mark.parametrize(
        "value",
        ["plain.fa", "od d's dir/ref.fa", "chr1:1-2", "a;rm -rf /", "", 8],
    )
    def test_it_agrees_with_the_helper_it_copies(self, value):
        """`_quote` is inlined only to keep the base-image hash small, not to differ."""
        assert install_references._quote(value) == quote_path(value)

    def test_a_path_with_a_space_survives_as_one_operand(self, tmp_path):
        """The property that matters: bash sees one word, not two."""
        path = tmp_path / "my refs" / "chr1.fa"
        assert shlex.split(f"samtools faidx {install_references._quote(path)}")[-1] == str(path)


# =============================================================================
# derive_region_fasta
# =============================================================================


class TestDeriveRegionFasta:
    def test_the_command_names_the_region_and_redirects_to_the_destination(self, tmp_path, monkeypatch):
        calls: list[str] = []
        monkeypatch.setattr(
            install_references.subprocess,
            "run",
            lambda cmd, **kw: calls.append(cmd) or subprocess.CompletedProcess(cmd, 0, "", ""),
        )

        destination = tmp_path / "muc1_region_hg38.fa"
        install_references.derive_region_fasta(
            source_fasta=tmp_path / "chr1.hg38.fa",
            region="chr1:155184000-155194000",
            destination=destination,
            samtools="samtools",
        )

        assert "faidx" in calls[0]
        assert "chr1:155184000-155194000" in calls[0]
        assert str(destination) in calls[0]

    def test_paths_with_spaces_are_quoted_into_single_operands(self, tmp_path, monkeypatch):
        """run_command runs this under bash, so an unquoted path with a space splits."""
        calls: list[str] = []
        monkeypatch.setattr(
            install_references.subprocess,
            "run",
            lambda cmd, **kw: calls.append(cmd) or subprocess.CompletedProcess(cmd, 0, "", ""),
        )

        source = tmp_path / "my refs" / "chr1.hg38.fa"
        install_references.derive_region_fasta(source, "chr1:1-2", tmp_path / "o.fa", "samtools")

        assert str(source) in shlex.split(calls[0])

    def test_a_failed_derivation_raises_rather_than_producing_an_empty_file(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            install_references.subprocess,
            "run",
            lambda cmd, **kw: subprocess.CompletedProcess(cmd, 1, "", "no such region"),
        )
        with pytest.raises(RuntimeError, match="faidx|deriv"):
            install_references.derive_region_fasta(tmp_path / "c.fa", "chr1:1-2", tmp_path / "o.fa", "samtools")

    def test_the_failure_message_names_the_file_and_repeats_samtools_stderr(self, tmp_path, monkeypatch):
        """An operator has to be able to act on this without re-running by hand."""
        monkeypatch.setattr(
            install_references.subprocess,
            "run",
            lambda cmd, **kw: subprocess.CompletedProcess(cmd, 1, "", "[faidx] Truncated file\n"),
        )
        with pytest.raises(RuntimeError) as excinfo:
            install_references.derive_region_fasta(tmp_path / "c.fa", "chr1:1-2", tmp_path / "o.fa", "samtools")

        assert "o.fa" in str(excinfo.value)
        assert "Truncated file" in str(excinfo.value)

    def test_the_samtools_prefix_is_left_unquoted_so_it_can_be_several_words(self, tmp_path, monkeypatch):
        """`samtools_path` is a command prefix like `mamba run -n env samtools`."""
        calls: list[str] = []
        monkeypatch.setattr(
            install_references.subprocess,
            "run",
            lambda cmd, **kw: calls.append(cmd) or subprocess.CompletedProcess(cmd, 0, "", ""),
        )

        install_references.derive_region_fasta(
            tmp_path / "c.fa", "chr1:1-2", tmp_path / "o.fa", "mamba run -n vntyper samtools"
        )

        assert calls[0].startswith("mamba run -n vntyper samtools faidx ")


# =============================================================================
# index_fasta_with_samtools
# =============================================================================


class TestIndexFastaWithSamtools:
    def test_it_faidxes_the_file(self, tmp_path, monkeypatch):
        calls: list[str] = []
        monkeypatch.setattr(
            install_references.subprocess,
            "run",
            lambda cmd, **kw: calls.append(cmd) or subprocess.CompletedProcess(cmd, 0, "", ""),
        )
        fasta = tmp_path / "my refs" / "out.fa"

        install_references.index_fasta_with_samtools(fasta, "samtools")

        assert shlex.split(calls[0]) == ["samtools", "faidx", str(fasta)]

    def test_a_failure_raises_and_names_the_file(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            install_references.subprocess,
            "run",
            lambda cmd, **kw: subprocess.CompletedProcess(cmd, 1, "", "not a fasta"),
        )
        with pytest.raises(RuntimeError) as excinfo:
            install_references.index_fasta_with_samtools(tmp_path / "out.fa", "samtools")

        assert "out.fa" in str(excinfo.value)
        assert "not a fasta" in str(excinfo.value)


# =============================================================================
# merge_pairwise_motifs
# =============================================================================


class TestMergePairwiseMotifs:
    """The literal derivation, ported out of `reference/generate_vntr_reference.py`."""

    @staticmethod
    def _seeds(tmp_path: Path, filter_config: dict) -> tuple[Path, Path]:
        seed = tmp_path / "MUC1_motifs_Rev_com.fa"
        seed.write_text(">A\nAAA\n>B\nBBB\n>C\nCCC\n", encoding="utf-8")
        config = tmp_path / "filter_config.json"
        config.write_text(json.dumps(filter_config), encoding="utf-8")
        return seed, config

    def test_every_ordered_pair_including_self_pairs_is_written_in_file_order(self, tmp_path):
        """`itertools.product` over an insertion-ordered dict is what makes this stable."""
        seed, config = self._seeds(tmp_path, {})
        destination = tmp_path / "merged.fa"

        install_references.merge_pairwise_motifs(seed, config, destination)

        assert destination.read_text(encoding="utf-8") == (
            ">A-A\nAAAAAA\n>A-B\nBBBAAA\n>A-C\nCCCAAA\n"
            ">B-A\nAAABBB\n>B-B\nBBBBBB\n>B-C\nCCCBBB\n"
            ">C-A\nAAACCC\n>C-B\nBBBCCC\n>C-C\nCCCCCC\n"
        )

    def test_the_sequence_order_is_the_reverse_of_the_header_order(self, tmp_path):
        """`>h1-h2` carries `seq(h2) + seq(h1)`; the inversion is load-bearing for the digest."""
        seed, config = self._seeds(tmp_path, {})
        destination = tmp_path / "merged.fa"

        install_references.merge_pairwise_motifs(seed, config, destination)

        records = destination.read_text(encoding="utf-8").splitlines()
        assert records[records.index(">A-B") + 1] == "BBBAAA"

    def test_a_disallowed_partner_drops_exactly_that_ordered_pair(self, tmp_path):
        """The filter is directional: banning B from A says nothing about A from B."""
        seed, config = self._seeds(tmp_path, {"A": ["B"]})
        destination = tmp_path / "merged.fa"

        install_references.merge_pairwise_motifs(seed, config, destination)

        headers = [line for line in destination.read_text(encoding="utf-8").splitlines() if line.startswith(">")]
        assert ">A-B" not in headers
        assert ">B-A" in headers
        assert len(headers) == 8

    def test_a_contig_absent_from_the_filter_config_keeps_every_partner(self, tmp_path):
        seed, config = self._seeds(tmp_path, {"A": ["B"]})
        destination = tmp_path / "merged.fa"

        install_references.merge_pairwise_motifs(seed, config, destination)

        headers = [line for line in destination.read_text(encoding="utf-8").splitlines() if line.startswith(">")]
        assert [h for h in headers if h.startswith(">C-")] == [">C-A", ">C-B", ">C-C"]

    def test_the_output_is_not_line_wrapped(self, tmp_path):
        """The tracked file has one line of sequence per record; wrapping changes the digest."""
        seed = tmp_path / "seed.fa"
        seed.write_text(">A\n" + "A" * 120 + "\n", encoding="utf-8")
        config = tmp_path / "filter_config.json"
        config.write_text("{}", encoding="utf-8")
        destination = tmp_path / "merged.fa"

        install_references.merge_pairwise_motifs(seed, config, destination)

        assert destination.read_text(encoding="utf-8") == ">A-A\n" + "A" * 240 + "\n"

    def test_a_seed_whose_headers_and_sequences_do_not_pair_up_is_rejected(self, tmp_path):
        """Silently zipping a ragged seed would emit a shorter file under the same name."""
        seed = tmp_path / "seed.fa"
        seed.write_text(">A\nAAA\n>B\n", encoding="utf-8")
        config = tmp_path / "filter_config.json"
        config.write_text("{}", encoding="utf-8")

        with pytest.raises(ValueError, match="2 headers"):
            install_references.merge_pairwise_motifs(seed, config, tmp_path / "merged.fa")

    def test_a_seed_with_no_records_is_rejected(self, tmp_path):
        seed = tmp_path / "seed.fa"
        seed.write_text("\n", encoding="utf-8")
        config = tmp_path / "filter_config.json"
        config.write_text("{}", encoding="utf-8")

        with pytest.raises(ValueError, match="no records"):
            install_references.merge_pairwise_motifs(seed, config, tmp_path / "merged.fa")


# =============================================================================
# resolve_source_location
# =============================================================================


class TestResolveSourceLocation:
    def test_without_a_spec_the_committed_url_and_digest_are_used(self):
        url, digest = install_references.resolve_source_location(
            "hg19", {"url": "https://ucsc/chr1.fa.gz", "source_sha256": "a" * 64}, None
        )
        assert (url, digest) == ("https://ucsc/chr1.fa.gz", "a" * 64)

    def test_a_release_spec_overrides_both_fields(self):
        spec = {"sources": {"hg19": {"url": "https://pinned/chr1.fa.gz", "source_sha256": "b" * 64}}}
        url, digest = install_references.resolve_source_location(
            "hg19", {"url": "https://ucsc/chr1.fa.gz", "source_sha256": "a" * 64}, spec
        )
        assert (url, digest) == ("https://pinned/chr1.fa.gz", "b" * 64)

    def test_a_reference_the_spec_does_not_name_falls_back_to_the_config(self):
        spec = {"sources": {"hg38": {"url": "https://pinned/other.gz"}}}
        url, digest = install_references.resolve_source_location(
            "hg19", {"url": "https://ucsc/chr1.fa.gz", "source_sha256": "a" * 64}, spec
        )
        assert (url, digest) == ("https://ucsc/chr1.fa.gz", "a" * 64)

    def test_a_missing_digest_is_refused_rather_than_installed_unverified(self):
        with pytest.raises(ValueError) as excinfo:
            install_references.resolve_source_location("hg19", {"url": "https://ucsc/chr1.fa.gz"}, None)

        assert "hg19" in str(excinfo.value)
        assert "source_sha256" in str(excinfo.value)

    def test_a_missing_url_is_refused(self):
        with pytest.raises(ValueError) as excinfo:
            install_references.resolve_source_location("hg19", {"source_sha256": "a" * 64}, None)

        assert "hg19" in str(excinfo.value)


# =============================================================================
# run_derivations
# =============================================================================


class TestRunDerivations:
    @staticmethod
    def _shark_config(expected: str) -> dict:
        return {
            "samtools_path": "samtools",
            "derivations": [
                {
                    "kind": "shark",
                    "assembly": "hg19",
                    "output": "muc1_region_hg19.fa",
                    "from": "hg19",
                    "region": "chr1:155158000-155163000",
                    "tool": "samtools faidx",
                    "expected_sha256": expected,
                }
            ],
        }

    def test_a_shark_derivation_is_cut_verified_and_indexed(self, tmp_path, monkeypatch):
        payload = b">chr1:155158000-155163000\nACGT\n"
        calls: list[str] = []
        monkeypatch.setattr(
            install_references.subprocess,
            "run",
            _fake_samtools({"chr1:155158000-155163000": payload}, calls),
        )
        source = tmp_path / "chr1.hg19.fa"
        source.write_text(">chr1\nACGT\n", encoding="utf-8")

        install_references.run_derivations(
            self._shark_config(hashlib.sha256(payload).hexdigest()), tmp_path, {"hg19": source}
        )

        assert (tmp_path / "muc1_region_hg19.fa").read_bytes() == payload
        assert (tmp_path / "muc1_region_hg19.fa.fai").exists()
        assert len(calls) == 2

    def test_a_derived_file_that_misses_its_digest_raises(self, tmp_path, monkeypatch):
        """This assertion is the whole reason the tracked FASTAs can be deleted."""
        calls: list[str] = []
        monkeypatch.setattr(
            install_references.subprocess,
            "run",
            _fake_samtools({"chr1:155158000-155163000": b">wrong\nTTTT\n"}, calls),
        )
        source = tmp_path / "chr1.hg19.fa"
        source.write_text(">chr1\nACGT\n", encoding="utf-8")

        with pytest.raises(ValueError, match="checksum mismatch"):
            install_references.run_derivations(self._shark_config("f" * 64), tmp_path, {"hg19": source})

        assert not (tmp_path / "muc1_region_hg19.fa.fai").exists(), "indexing must not follow a failed digest"

    def test_a_source_that_was_never_installed_fails_loudly(self, tmp_path):
        with pytest.raises(RuntimeError) as excinfo:
            install_references.run_derivations(self._shark_config("a" * 64), tmp_path, {})

        assert "muc1_region_hg19.fa" in str(excinfo.value)
        assert "hg19" in str(excinfo.value)

    def test_a_source_recorded_but_missing_from_disk_fails_loudly(self, tmp_path):
        with pytest.raises(RuntimeError, match="not installed"):
            install_references.run_derivations(self._shark_config("a" * 64), tmp_path, {"hg19": tmp_path / "gone.fa"})

    def test_a_literal_derivation_merges_its_seeds(self, tmp_path, monkeypatch):
        calls: list[str] = []
        monkeypatch.setattr(install_references.subprocess, "run", _fake_samtools({}, calls))
        (tmp_path / "MUC1_motifs_Rev_com.fa").write_text(">A\nAAA\n", encoding="utf-8")
        (tmp_path / "filter_config.json").write_text("{}", encoding="utf-8")
        payload = b">A-A\nAAAAAA\n"
        config = {
            "derivations": [
                {
                    "kind": "literal",
                    "config_key": "muc1_reference_vntr",
                    "output": "merged.fa",
                    "from_seeds": ["MUC1_motifs_Rev_com.fa", "filter_config.json"],
                    "tool": "generate_vntr_reference",
                    "expected_sha256": hashlib.sha256(payload).hexdigest(),
                }
            ]
        }

        install_references.run_derivations(config, tmp_path, {})

        assert (tmp_path / "merged.fa").read_bytes() == payload
        assert (tmp_path / "merged.fa.fai").exists()

    def test_a_literal_derivation_with_a_missing_seed_names_the_seed(self, tmp_path):
        (tmp_path / "MUC1_motifs_Rev_com.fa").write_text(">A\nAAA\n", encoding="utf-8")
        config = {
            "derivations": [
                {
                    "kind": "literal",
                    "output": "merged.fa",
                    "from_seeds": ["MUC1_motifs_Rev_com.fa", "filter_config.json"],
                    "expected_sha256": "a" * 64,
                }
            ]
        }

        with pytest.raises(RuntimeError) as excinfo:
            install_references.run_derivations(config, tmp_path, {})

        assert "filter_config.json" in str(excinfo.value)
        assert "merged.fa" in str(excinfo.value)

    def test_a_literal_derivation_needs_exactly_two_seeds(self, tmp_path):
        config = {
            "derivations": [{"kind": "literal", "output": "m.fa", "from_seeds": ["a.fa"], "expected_sha256": "a" * 64}]
        }
        with pytest.raises(ValueError, match="two seeds"):
            install_references.run_derivations(config, tmp_path, {})

    def test_an_unknown_kind_is_refused_instead_of_silently_skipped(self, tmp_path):
        config = {"derivations": [{"kind": "sorcery", "output": "m.fa", "expected_sha256": "a" * 64}]}
        with pytest.raises(ValueError, match="sorcery"):
            install_references.run_derivations(config, tmp_path, {})

    def test_a_derivation_without_a_committed_digest_is_refused(self, tmp_path):
        config = {"derivations": [{"kind": "shark", "output": "m.fa", "from": "hg19", "region": "chr1:1-2"}]}
        with pytest.raises(ValueError, match="expected_sha256"):
            install_references.run_derivations(config, tmp_path, {"hg19": tmp_path})

    def test_no_derivations_configured_is_a_no_op(self, tmp_path):
        install_references.run_derivations({}, tmp_path, {})
        assert list(tmp_path.iterdir()) == []


# =============================================================================
# install_from_source
# =============================================================================


def _gz(path: Path, payload: bytes) -> None:
    """Write `payload` gzip-compressed to `path`, creating parents.

    `mtime=0` and an empty `filename` keep the bytes deterministic: gzip otherwise
    stamps the current time and the member name into its header, so the same payload
    would hash differently on every call and no digest could be asserted.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
        handle.write(payload)


class TestInstallFromSource:
    """The definition of what a published bundle contains."""

    GENOME = b">chr1\nACGTACGTAC\n"

    def _config(self, digest: str, *, derivations: list | None = None) -> dict:
        return {
            "samtools_path": "samtools",
            "ucsc_references": {
                "hg19": {
                    "url": "https://hgdownload.example/hg19/chr1.fa.gz",
                    "target_path": "alignment/chr1.hg19.fa.gz",
                    "source_sha256": digest,
                }
            },
            "derivations": derivations or [],
        }

    def _download(self, payload: bytes, urls: list[str]):
        def _fake(url, dest_path):
            urls.append(url)
            _gz(dest_path, payload)

        return _fake

    def _gz_digest(self, tmp_path: Path, payload: bytes) -> str:
        archive = tmp_path / "probe.gz"
        _gz(archive, payload)
        digest = hashlib.sha256(archive.read_bytes()).hexdigest()
        archive.unlink()
        return digest

    def test_it_downloads_verifies_decompresses_indexes_and_returns_the_mapping(self, tmp_path, monkeypatch):
        urls: list[str] = []
        calls: list[str] = []
        digest = self._gz_digest(tmp_path, self.GENOME)
        monkeypatch.setattr(install_references, "download_file", self._download(self.GENOME, urls))
        monkeypatch.setattr(install_references.subprocess, "run", _fake_samtools({}, calls))
        indexed: list[Path] = []
        monkeypatch.setattr(
            install_references,
            "index_reference_with_aligners",
            lambda ref_path, aligners, threads=4, force_reindex=False: indexed.append(ref_path) or {},
        )

        installed = install_references.install_from_source(
            self._config(digest), tmp_path, ["hg19"], {"bwa": {}}, index_threads=2
        )

        fasta = tmp_path / "alignment" / "chr1.hg19.fa"
        assert installed == {"hg19": fasta}
        assert fasta.read_bytes() == self.GENOME
        assert indexed == [fasta]
        assert urls == ["https://hgdownload.example/hg19/chr1.fa.gz"]
        assert (tmp_path / "alignment" / "chr1.hg19.fa.fai").exists()

    def test_a_release_spec_overrides_the_config_url(self, tmp_path, monkeypatch):
        urls: list[str] = []
        digest = self._gz_digest(tmp_path, self.GENOME)
        monkeypatch.setattr(install_references, "download_file", self._download(self.GENOME, urls))
        monkeypatch.setattr(install_references.subprocess, "run", _fake_samtools({}, []))

        install_references.install_from_source(
            self._config("a" * 64),
            tmp_path,
            ["hg19"],
            {},
            index_threads=1,
            release_spec={"sources": {"hg19": {"url": "https://pinned.example/chr1.fa.gz", "source_sha256": digest}}},
        )

        assert urls == ["https://pinned.example/chr1.fa.gz"]

    def test_a_digest_mismatch_fails_before_any_decompression_or_indexing(self, tmp_path, monkeypatch):
        urls: list[str] = []
        calls: list[str] = []
        monkeypatch.setattr(install_references, "download_file", self._download(self.GENOME, urls))
        monkeypatch.setattr(install_references.subprocess, "run", _fake_samtools({}, calls))
        aligner_calls: list[Path] = []
        monkeypatch.setattr(
            install_references,
            "index_reference_with_aligners",
            lambda ref_path, aligners, threads=4, force_reindex=False: aligner_calls.append(ref_path) or {},
        )

        with pytest.raises(ValueError, match="checksum mismatch"):
            install_references.install_from_source(
                self._config("f" * 64), tmp_path, ["hg19"], {"bwa": {}}, index_threads=1
            )

        assert not (tmp_path / "alignment" / "chr1.hg19.fa").exists(), "decompression ran despite a bad digest"
        assert calls == [], "samtools ran despite a bad digest"
        assert aligner_calls == [], "an aligner ran despite a bad digest"

    def test_the_mismatch_message_names_the_file_and_both_digests(self, tmp_path, monkeypatch):
        digest = self._gz_digest(tmp_path, self.GENOME)
        monkeypatch.setattr(install_references, "download_file", self._download(self.GENOME, []))

        with pytest.raises(ValueError) as excinfo:
            install_references.install_from_source(self._config("f" * 64), tmp_path, ["hg19"], {}, index_threads=1)

        message = str(excinfo.value)
        assert "chr1.hg19.fa.gz" in message
        assert "f" * 64 in message
        assert digest in message

    def test_skipping_indexing_is_expressed_as_no_aligners(self, tmp_path, monkeypatch):
        digest = self._gz_digest(tmp_path, self.GENOME)
        monkeypatch.setattr(install_references, "download_file", self._download(self.GENOME, []))
        monkeypatch.setattr(install_references.subprocess, "run", _fake_samtools({}, []))
        called: list[Path] = []
        monkeypatch.setattr(
            install_references,
            "index_reference_with_aligners",
            lambda ref_path, aligners, threads=4, force_reindex=False: called.append(ref_path) or {},
        )

        install_references.install_from_source(self._config(digest), tmp_path, ["hg19"], {}, index_threads=1)

        assert called == []

    def test_a_reference_that_is_not_selected_is_not_downloaded(self, tmp_path, monkeypatch):
        urls: list[str] = []
        monkeypatch.setattr(install_references, "download_file", self._download(self.GENOME, urls))
        monkeypatch.setattr(install_references.subprocess, "run", _fake_samtools({}, []))

        installed = install_references.install_from_source(
            self._config("a" * 64), tmp_path, ["hg38"], {}, index_threads=1
        )

        assert installed == {}
        assert urls == []

    def test_an_uncompressed_source_is_used_where_it_lands(self, tmp_path, monkeypatch):
        """Not every upstream ships `.gz`; only `.gz` gets a decompression step."""
        payload = b">chr1\nACGT\n"
        digest = hashlib.sha256(payload).hexdigest()

        def _fake(url, dest_path):
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            dest_path.write_bytes(payload)

        monkeypatch.setattr(install_references, "download_file", _fake)
        monkeypatch.setattr(install_references.subprocess, "run", _fake_samtools({}, []))
        config = {
            "ucsc_references": {
                "hg19": {"url": "https://x/chr1.fa", "target_path": "alignment/chr1.fa", "source_sha256": digest}
            }
        }

        installed = install_references.install_from_source(config, tmp_path, ["hg19"], {}, index_threads=1)

        assert installed == {"hg19": tmp_path / "alignment" / "chr1.fa"}

    def test_a_failed_download_raises_instead_of_ending_the_process(self, tmp_path, monkeypatch):
        """`download_file` calls `sys.exit`; a bundle build needs an exception it can unwind."""

        def _fail(url, dest_path):
            raise SystemExit(1)

        monkeypatch.setattr(install_references, "download_file", _fail)

        with pytest.raises(RuntimeError) as excinfo:
            install_references.install_from_source(self._config("a" * 64), tmp_path, ["hg19"], {}, index_threads=1)

        assert "hg19" in str(excinfo.value)
        assert "https://hgdownload.example/hg19/chr1.fa.gz" in str(excinfo.value)

    def test_an_entry_without_a_target_path_is_refused(self, tmp_path, monkeypatch):
        monkeypatch.setattr(install_references, "download_file", self._download(self.GENOME, []))
        config = {"ucsc_references": {"hg19": {"url": "https://x/chr1.fa.gz", "source_sha256": "a" * 64}}}

        with pytest.raises(ValueError, match="target_path"):
            install_references.install_from_source(config, tmp_path, ["hg19"], {}, index_threads=1)

    def test_the_seeds_the_derivations_produce_are_not_also_downloaded(self, tmp_path, monkeypatch):
        """`All_Pairwise_...` is derived now, so fetching it first would be a wasted round trip."""
        requested: list[str] = []
        monkeypatch.setattr(
            install_references,
            "process_own_repository_references",
            lambda refs, out, skip, md5: requested.extend(f["target_path"] for f in refs["raw_files"]),
        )
        monkeypatch.setattr(install_references.subprocess, "run", _fake_samtools({}, []))
        config = {
            "own_repository_references": {
                "raw_files": [
                    {"url": "https://x/merged.fa", "target_path": "merged.fa"},
                    {"url": "https://x/seed.fa", "target_path": "MUC1_motifs_Rev_com.fa"},
                ]
            },
            # The literal derivation owns merged.fa; its digest is deliberately wrong so
            # the run stops right after the merge, once the seed list has been settled.
            "derivations": [
                {"kind": "literal", "output": "merged.fa", "from_seeds": ["a", "b"], "expected_sha256": "a" * 64}
            ],
        }
        (tmp_path / "a").write_text(">A\nAAA\n", encoding="utf-8")
        (tmp_path / "b").write_text("{}", encoding="utf-8")

        with pytest.raises(ValueError, match="checksum mismatch"):
            install_references.install_from_source(config, tmp_path, [], {}, index_threads=1)

        assert requested == ["MUC1_motifs_Rev_com.fa"]

    def test_a_raw_file_list_that_is_entirely_derived_fetches_nothing(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            install_references, "process_own_repository_references", _forbidden("process_own_repository_references")
        )
        config = {
            "own_repository_references": {"raw_files": [{"url": "https://x/merged.fa", "target_path": "merged.fa"}]},
            "derivations": [{"kind": "literal", "output": "merged.fa", "from_seeds": ["a", "b"]}],
        }

        with pytest.raises(ValueError, match="expected_sha256"):
            install_references.install_from_source(config, tmp_path, [], {}, index_threads=1)

    def test_a_selected_vntyper_reference_is_installed_alongside_the_genomes(self, tmp_path, monkeypatch):
        """The adVNTR database zip is selected by the same `--references` list."""
        seen: list[str] = []
        monkeypatch.setattr(
            install_references,
            "process_vntyper_references",
            lambda refs, out, bwa, skip, md5: seen.extend(refs),
        )
        config = {
            "vntyper_references": {
                "vntr_db_advntr": {"url": "https://x/db.zip", "target_path": "db.zip", "extract_to": "."}
            }
        }

        install_references.install_from_source(config, tmp_path, ["vntr_db_advntr"], {}, index_threads=1)

        assert seen == ["vntr_db_advntr"]


# =============================================================================
# main routing and the CLI surface
# =============================================================================


def _forbidden(name: str, exception: type[Exception] = AssertionError):
    """Return a callable that fails the test if the wrong code path calls it."""

    def _raise(*args, **kwargs):
        raise exception(f"{name} must not be called on this path")

    return _raise


class TestMainRouting:
    def test_from_source_takes_the_source_path_and_not_the_legacy_one(self, tmp_path, monkeypatch):
        seen: dict = {}
        monkeypatch.setattr(
            install_references,
            "install_from_source",
            lambda cfg, out, refs, aligners, index_threads, release_spec=None: (
                seen.update(references=refs, output_dir=out, release_spec=release_spec) or {}
            ),
        )
        for legacy in ("process_ucsc_references", "process_vntyper_references", "process_own_repository_references"):
            monkeypatch.setattr(install_references, legacy, _forbidden(legacy))

        install_references.main(tmp_path, skip_indexing=True, references_to_process=["hg19"], from_source=True)

        assert seen["references"] == ["hg19"]
        assert seen["output_dir"] == tmp_path
        assert seen["release_spec"] is None

    def test_without_the_flag_the_legacy_path_is_unchanged(self, tmp_path, monkeypatch):
        seen: list[str] = []
        monkeypatch.setattr(
            install_references, "install_from_source", _forbidden("install_from_source", exception=AssertionError)
        )
        monkeypatch.setattr(
            install_references,
            "process_ucsc_references",
            lambda refs, out, bwa, skip, md5, aligners=None, index_threads=4: seen.extend(refs),
        )
        monkeypatch.setattr(install_references, "process_vntyper_references", lambda refs, out, bwa, skip, md5: None)
        monkeypatch.setattr(install_references, "process_own_repository_references", lambda refs, out, skip, md5: None)

        install_references.main(tmp_path, skip_indexing=True, references_to_process=["hg19"])

        assert seen == ["hg19"]

    def test_a_release_spec_path_is_loaded_and_handed_over(self, tmp_path, monkeypatch):
        spec_path = tmp_path / "spec.json"
        spec_path.write_text(json.dumps({"sources": {"hg19": {"url": "https://pinned/x.gz"}}}), encoding="utf-8")
        seen: dict = {}
        monkeypatch.setattr(
            install_references,
            "install_from_source",
            lambda cfg, out, refs, aligners, index_threads, release_spec=None: seen.update(spec=release_spec) or {},
        )

        install_references.main(
            tmp_path, skip_indexing=True, references_to_process=["hg19"], from_source=True, release_spec_path=spec_path
        )

        assert seen["spec"] == {"sources": {"hg19": {"url": "https://pinned/x.gz"}}}


class TestCliSurface:
    def test_from_source_defaults_off_and_parses(self):
        parser = build_parser()
        assert parser.parse_args(["install-references", "-d", "refs"]).from_source is False
        assert parser.parse_args(["install-references", "-d", "refs", "--from-source"]).from_source is True

    def test_release_spec_parses_as_a_path(self):
        parser = build_parser()
        assert parser.parse_args(["install-references", "-d", "refs"]).release_spec is None
        args = parser.parse_args(["install-references", "-d", "refs", "--release-spec", "spec.json"])
        assert args.release_spec == Path("spec.json")

    def test_the_handler_threads_both_through_to_main(self, monkeypatch):
        import argparse

        from vntyper.scripts import cli_handlers

        seen: dict = {}
        monkeypatch.setattr(cli_handlers, "install_references_main", lambda **kwargs: seen.update(kwargs))
        args = argparse.Namespace(
            output_dir=Path("refs"),
            config_path=None,
            skip_indexing=False,
            threads=4,
            aligners=None,
            references=None,
            from_source=True,
            release_spec=Path("spec.json"),
        )

        with pytest.raises(SystemExit):
            cli_handlers.handle_install_references(args, {}, None, 20, None)

        assert seen["from_source"] is True
        assert seen["release_spec_path"] == Path("spec.json")


# =============================================================================
# The shipped configuration
# =============================================================================


class TestShippedConfig:
    """The committed digests are the contract; a drifting one must fail here first."""

    def test_the_three_derivations_carry_the_measured_digests(self):
        derivations = {spec["output"]: spec for spec in _shipped_config()["derivations"]}

        assert derivations["muc1_region_hg19.fa"]["expected_sha256"] == HG19_REGION_SHA256
        assert derivations["muc1_region_hg38.fa"]["expected_sha256"] == HG38_REGION_SHA256
        assert (
            derivations["All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa"]["expected_sha256"]
            == MERGED_MOTIFS_SHA256
        )

    def test_every_derivation_declares_the_fields_the_dispatch_reads(self):
        for spec in _shipped_config()["derivations"]:
            assert spec["kind"] in {"shark", "literal"}
            assert spec["output"] and spec["expected_sha256"] and spec["tool"]
            if spec["kind"] == "shark":
                assert spec["assembly"] and spec["from"] and spec["region"]
            else:
                assert spec["config_key"] and len(spec["from_seeds"]) == 2

    def test_the_shark_regions_are_the_ones_the_digests_were_measured_from(self):
        regions = {spec["output"]: spec["region"] for spec in _shipped_config()["derivations"] if "region" in spec}
        assert regions == {
            "muc1_region_hg19.fa": "chr1:155158000-155163000",
            "muc1_region_hg38.fa": "chr1:155184000-155194000",
        }

    def test_every_genome_carries_a_source_digest(self):
        config = _shipped_config()
        digests = {}
        for section in ("ucsc_references", "ncbi_references", "ensembl_references"):
            for ref_id, entry in config[section].items():
                digests[ref_id] = entry.get("source_sha256")

        assert digests == SOURCE_SHA256

    def test_ensembl_is_pinned_to_an_explicit_release_not_current(self):
        """A digest beside a mutable `current` URL rots at the next Ensembl release."""
        ensembl = _shipped_config()["ensembl_references"]

        for ref_id, entry in ensembl.items():
            assert "/current" not in entry["url"], f"{ref_id} still points at a mutable path"
            assert "release-116" in entry["url"], f"{ref_id} is not pinned to the measured release"

    def test_the_samtools_prefix_is_configurable(self):
        assert _shipped_config()["samtools_path"] == "samtools"

    def test_every_common_reference_names_a_config_key_and_an_installed_path(self):
        for entry in _shipped_config()["common_references"]:
            assert entry["config_key"]
            assert entry["installed_path"]
