"""BWA reference selection, and the safety guard that shares it.

Two readers used to collapse every assembly onto bwa_reference_hg19/hg38 independently:
cli_handlers for the run, and cli_logging_safety._selected_bwa_reference for the check
that stops --log-file naming an operator input. Because the guard runs before
setup_logging opens the log in append mode, a guard looking at the wrong file means
logging can append into a reference FASTA.
"""

import json
import logging
from pathlib import Path
from unittest import mock

import pytest

from vntyper.scripts import cli_handlers
from vntyper.scripts.artifact_names import PIPELINE_BASENAME
from vntyper.scripts.cli_handlers import select_bwa_reference
from vntyper.scripts.cli_parser import build_parser

pytestmark = pytest.mark.unit


def config_with(**reference_data):
    return {"reference_data": dict(reference_data)}


class TestSelection:
    def test_the_physical_key_is_preferred(self, tmp_path):
        ensembl = tmp_path / "ens.fa"
        ucsc = tmp_path / "ucsc.fa"
        ensembl.touch()
        ucsc.touch()
        cfg = config_with(bwa_reference_hg38_ensembl=str(ensembl), bwa_reference_hg38=str(ucsc))
        assert select_bwa_reference(cfg, "hg38_ensembl") == str(ensembl)

    @pytest.mark.parametrize("label", ["GRCh38", "hg38_ncbi"])
    def test_both_ncbi_labels_select_the_ncbi_reference(self, tmp_path, label):
        ncbi = tmp_path / "ncbi.fna"
        ucsc = tmp_path / "ucsc.fa"
        ncbi.touch()
        ucsc.touch()
        cfg = config_with(bwa_reference_GRCh38=str(ncbi), bwa_reference_hg38=str(ucsc))
        assert select_bwa_reference(cfg, label) == str(ncbi)

    def test_the_ucsc_family_fallback_warns_and_names_both_keys(self, tmp_path, caplog):
        ucsc = tmp_path / "ucsc.fa"
        ucsc.touch()
        cfg = config_with(bwa_reference_hg38=str(ucsc))
        with caplog.at_level("WARNING"):
            assert select_bwa_reference(cfg, "hg38_ensembl") == str(ucsc)
        message = " ".join(record.getMessage() for record in caplog.records)
        assert "hg38_ensembl" in message and "bwa_reference_hg38" in message and "ucsc" in message.lower()

    def test_an_exact_null_is_authoritative_and_fails_closed(self, tmp_path):
        ucsc = tmp_path / "ucsc.fa"
        ucsc.touch()
        cfg = config_with(bwa_reference_hg38_ensembl=None, bwa_reference_hg38=str(ucsc))
        with pytest.raises(ValueError, match="hg38_ensembl"):
            select_bwa_reference(cfg, "hg38_ensembl")

    def test_an_exact_null_is_distinguished_from_absent_in_the_message(self):
        """Item 11: a key present with an explicit ``null`` was deliberately disabled,
        which is a different operator story from a key that was never configured at
        all - the message must say so, the way `select_muc1_region_fasta` already does
        for SHARK."""
        cfg = config_with(bwa_reference_hg38_ensembl=None, bwa_reference_hg38=None)
        with pytest.raises(ValueError, match="disabled") as excinfo:
            select_bwa_reference(cfg, "hg38_ensembl")
        assert "bwa_reference_hg38_ensembl" in str(excinfo.value)

    def test_nothing_configured_fails_closed_naming_every_key_tried(self):
        with pytest.raises(ValueError) as excinfo:
            select_bwa_reference(config_with(), "hg38_ensembl")
        assert "bwa_reference_hg38_ensembl" in str(excinfo.value)
        assert "bwa_reference_hg38" in str(excinfo.value)
        assert "disabled" not in str(excinfo.value), "an absent key is not the same story as a disabled one"

    def test_a_malformed_present_value_fails_closed_with_value_error_not_type_error(self):
        """A present-but-non-string config value (an int here) must not reach `Path()`
        further down the call chain and crash with a bare `TypeError`; the required run
        path must fail closed the same way a missing key does."""
        cfg = config_with(bwa_reference_hg19=7)
        with pytest.raises(ValueError):
            select_bwa_reference(cfg, "hg19")


class TestThePresentButMissingFileFailsClosed:
    """Important 1: presence beats truthiness, but a configured path that names no file
    on disk is not a usable reference. A default `install-references` run only installs
    hg19/hg38, yet the shipped `config.json` also declares `bwa_reference_GRCh38` etc. as
    real relative paths - so a `--reference-assembly GRCh38` run used to resolve that
    unindexed path with `is_fallback=False` and no warning, and die stages later with a
    message naming neither the assembly nor the remedy.
    """

    def test_a_present_but_nonexistent_physical_path_fails_closed_with_an_actionable_message(self, tmp_path):
        missing = tmp_path / "reference" / "alignment" / "chr1.GRCh38.fna"
        cfg = config_with(bwa_reference_GRCh38=str(missing), bwa_reference_hg38=str(tmp_path / "ucsc.fa"))

        with pytest.raises(ValueError) as excinfo:
            select_bwa_reference(cfg, "GRCh38")

        message = str(excinfo.value)
        assert "GRCh38" in message, "must name the requested assembly"
        assert "bwa_reference_GRCh38" in message, "must name the key that resolved"
        assert str(missing) in message, "must name the missing path"
        assert "install-references" in message and "--references GRCh38" in message, (
            "must give the exact remedy, not a generic pointer"
        )

    def test_it_does_not_fall_through_to_the_ucsc_family_key(self, tmp_path):
        """The defect this milestone exists to kill: a GRCh38-labelled run must never
        quietly align against UCSC sequence just because the UCSC file happens to exist."""
        missing = tmp_path / "chr1.GRCh38.fna"
        ucsc = tmp_path / "chr1.hg38.fa"
        ucsc.write_text(">chr1\nACGT\n", encoding="utf-8")
        cfg = config_with(bwa_reference_GRCh38=str(missing), bwa_reference_hg38=str(ucsc))

        with pytest.raises(ValueError, match="GRCh38"):
            select_bwa_reference(cfg, "GRCh38")

    def test_the_same_config_resolves_normally_once_the_file_exists(self, tmp_path):
        """The complement of the two tests above: once `install-references` actually
        wrote the file, resolution must succeed exactly as before this fix."""
        present = tmp_path / "chr1.GRCh38.fna"
        present.write_text(">chr1\nACGT\n", encoding="utf-8")
        cfg = config_with(bwa_reference_GRCh38=str(present), bwa_reference_hg38=str(tmp_path / "ucsc.fa"))

        assert select_bwa_reference(cfg, "GRCh38") == str(present)

    def test_a_present_but_nonexistent_path_with_required_false_returns_none_not_raise(self, tmp_path):
        """The `required=False` logging-safety guard path (`cli_logging_safety`) must
        degrade to "no reference to protect" rather than raising - a missing file is not
        a log-destination collision to guard against."""
        missing = tmp_path / "chr1.GRCh38.fna"
        cfg = config_with(bwa_reference_GRCh38=str(missing))
        assert select_bwa_reference(cfg, "GRCh38", required=False) is None

    def test_the_ucsc_family_fallback_still_works_when_the_physical_key_is_genuinely_absent(self, tmp_path):
        """Unchanged path, must not regress: when the physical key is not present in the
        config at all (not present-with-a-broken-path), the family fallback still
        resolves to a real, installed UCSC file."""
        ucsc = tmp_path / "chr1.hg38.fa"
        ucsc.write_text(">chr1\nACGT\n", encoding="utf-8")
        cfg = config_with(bwa_reference_hg38=str(ucsc))

        assert select_bwa_reference(cfg, "GRCh38") == str(ucsc)


class TestThePartialInstallScenarioAgainstTheShippedConfig:
    """`tests/unit/conftest.py`'s `install_config` fixture touches every path
    `install_references_config.json` can write, so it only ever models a *complete*
    install - which is exactly why no existing test caught Important 1's defect. This
    class models the real, incomplete one instead: the shipped `vntyper/config.json`
    (not `install_config`, which never reaches the pipeline's own config) ships all six
    `bwa_reference_*` keys as real relative paths unconditionally, but a default
    `vntyper install-references -d reference` run only ever installs hg19 and hg38.
    """

    @staticmethod
    def _shipped_config() -> dict:
        config_path = Path(__file__).resolve().parents[2] / "vntyper" / "config.json"
        return json.loads(config_path.read_text(encoding="utf-8"))

    def test_the_shipped_config_ships_all_six_physical_keys_as_real_paths(self) -> None:
        """Pins the precondition the two tests below depend on: if `config.json` ever
        stops shipping these as real, non-null paths, this is the test that should fail,
        not the ones exercising the consequence."""
        reference_data = self._shipped_config()["reference_data"]
        for key in (
            "bwa_reference_hg19",
            "bwa_reference_hg38",
            "bwa_reference_GRCh37",
            "bwa_reference_GRCh38",
            "bwa_reference_hg19_ensembl",
            "bwa_reference_hg38_ensembl",
        ):
            assert isinstance(reference_data[key], str) and reference_data[key], key

    def test_a_default_hg19_hg38_only_install_fails_closed_for_grch38(self, tmp_path, monkeypatch) -> None:
        """The exact repro from the bug report:
        ``vntyper install-references -d reference`` (default: hg19 hg38), then
        ``vntyper pipeline --reference-assembly GRCh38 ...`` - must refuse, naming
        GRCh38 and the missing file, rather than silently aligning against UCSC hg38
        sequence with ``is_fallback=False`` and no warning."""
        config = self._shipped_config()
        monkeypatch.chdir(tmp_path)
        # Model the default install: only the two UCSC-tier physical files exist, at the
        # exact relative paths config.json names - resolved from the process CWD, the
        # same base the pipeline itself uses (AGENTS.md trap 7).
        for relative in (
            config["reference_data"]["bwa_reference_hg19"],
            config["reference_data"]["bwa_reference_hg38"],
        ):
            path = tmp_path / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()

        with pytest.raises(ValueError) as excinfo:
            select_bwa_reference(config, "GRCh38")

        message = str(excinfo.value)
        assert "GRCh38" in message, "must name the requested assembly"
        assert "install-references" in message and "--references GRCh38" in message, (
            "must give the exact remedy, not a generic pointer"
        )

    def test_the_same_shipped_config_resolves_hg38_normally(self, tmp_path, monkeypatch) -> None:
        """The complement: a label whose physical file the default install did write
        still resolves - this fix must not turn a correct default run into a failure."""
        config = self._shipped_config()
        monkeypatch.chdir(tmp_path)
        for relative in (
            config["reference_data"]["bwa_reference_hg19"],
            config["reference_data"]["bwa_reference_hg38"],
        ):
            path = tmp_path / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()

        assert select_bwa_reference(config, "hg38") == config["reference_data"]["bwa_reference_hg38"]


NO_BWA_CONFIG: dict = {
    "default_values": {
        "output_dir": "out",
        "threads": 4,
        "reference_assembly": "hg19",
        "output_name": PIPELINE_BASENAME,
        "archive_format": "zip",
    },
    "reference_data": {},
}


def _run_pipeline_handler(argv: list[str]) -> mock.MagicMock:
    """Parse ``argv`` and run ``handle_pipeline`` against :data:`NO_BWA_CONFIG`.

    Args:
        argv: Arguments after ``vntyper``, starting with ``pipeline``.

    Returns:
        unittest.mock.MagicMock: The ``run_pipeline`` stub.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    with mock.patch.object(cli_handlers, "run_pipeline", autospec=True) as stub:
        cli_handlers.handle_pipeline(
            args,
            config=NO_BWA_CONFIG,
            parser=parser,
            log_level_value=logging.INFO,
            log_file_str=None,
        )
    return stub


class TestTheRunPathOnlyRequiresBwaWhenItWillAlign:
    """`pipeline_inputs.py` raises for a missing BWA reference only when
    ``input_type == "FASTQ"``; the BAM and CRAM branches never read it
    (`if input_type == "FASTQ" and not bwa_reference: raise ValueError(...)`,
    `pipeline_inputs.py:153`). Resolving with `required=True` unconditionally would make
    `handle_pipeline` stricter than the consumer that actually needs the value, aborting
    a BAM/CRAM run a config with no BWA keys used to complete.
    """

    def test_a_bam_run_with_no_bwa_reference_configured_does_not_raise(self, tmp_path) -> None:
        stub = _run_pipeline_handler(["pipeline", "-o", str(tmp_path), "--bam", "in.bam"])
        assert stub.call_count == 1
        assert stub.call_args.kwargs["bwa_reference"] is None

    def test_a_cram_run_with_no_bwa_reference_configured_does_not_raise(self, tmp_path) -> None:
        stub = _run_pipeline_handler(["pipeline", "-o", str(tmp_path), "--cram", "in.cram"])
        assert stub.call_count == 1
        assert stub.call_args.kwargs["bwa_reference"] is None

    def test_a_fastq_run_with_no_bwa_reference_configured_still_fails_closed(self, tmp_path) -> None:
        with pytest.raises(ValueError, match="bwa_reference_hg19"):
            _run_pipeline_handler(
                ["pipeline", "-o", str(tmp_path), "--fastq1", "r1.fq", "--fastq2", "r2.fq"],
            )
