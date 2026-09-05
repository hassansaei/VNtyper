"""Unit tests for pure SHARK version discovery and parsing (#312).

SHARK 1.2.0 prints no version information to stdout or stderr when invoked with
``--version`` or ``--help``. For environments driven through ``mamba run -n <env>``
or ``conda run -n <env>``, version discovery inspects package metadata via
``conda list --json``.
"""

from __future__ import annotations

import pytest

from vntyper.scripts.shark_version import (
    build_shark_version_command,
    extract_conda_env_from_command,
    parse_shark_conda_list_json,
    parse_shark_help_output,
    resolve_shark_version_from_output,
)

pytestmark = pytest.mark.unit

# Captured verbatim from `mamba list -n shark_env shark --json` on 2026-09-05
PROBE_MAMBA_LIST_JSON = """[
    {
        "base_url": "https://conda.anaconda.org/bioconda",
        "build_number": 5,
        "build_string": "h077b44d_5",
        "channel": "bioconda",
        "dist_name": "shark-1.2.0-h077b44d_5",
        "md5": "b40d380390c2749c31c552c19b766165",
        "name": "shark",
        "platform": "linux-64",
        "sha256": "919e8f06a0020c167ead925016fa76bd3b5600fe0d44db0816ff234a1f85f890",
        "url": "https://conda.anaconda.org/bioconda/linux-64/shark-1.2.0-h077b44d_5.tar.bz2",
        "version": "1.2.0"
    }
]"""

PROBE_MAMBA_LIST_NO_BUILD_JSON = """[
    {
        "name": "shark",
        "version": "1.2.0"
    }
]"""

PROBE_EMPTY_LIST_JSON = "[]"

PROBE_NON_SHARK_LIST_JSON = """[
    {
        "name": "other_tool",
        "version": "0.9.1",
        "build_string": "py310_0"
    }
]"""

PROBE_HELP_OUTPUT = """Usage: shark -r <references> -1 <sample1> [OPTIONAL ARGUMENTS]

Arguments:
      -r, --reference                   reference sequences in FASTA format (can be gzipped)
      -1, --sample1                     sample in FASTQ (can be gzipped)

Optional arguments:
      -h, --help                        display this help and exit
      -2, --sample2                     second sample in FASTQ (optional, can be gzipped)
      -o, --out1                        first output sample in FASTQ (default: sharked_sample.1)
      -p, --out2                        second output sample in FASTQ (default: sharked_sample.2)
      -k, --kmer-size                   size of the kmers to index (default:17, max:31)
      -c, --confidence                  confidence for associating a read to a gene (default:0.6)
      -b, --bf-size                     bloom filter size in GB (default:1)
      -q, --min-base-quality            minimum base quality (assume FASTQ Illumina 1.8+ Phred scale, default:0, i.e., no filtering)
      -s, --single                      report an association only if a single gene is found
      -t, --threads                     number of threads (default:1)
      -v, --verbose                     verbose mode
"""


class TestExtractCondaEnvFromCommand:
    def test_extract_from_mamba_short_flag(self) -> None:
        result = extract_conda_env_from_command("mamba run -n shark_env shark")
        assert result == ("mamba", "shark_env")

    def test_extract_from_conda_long_flag(self) -> None:
        result = extract_conda_env_from_command("conda run --name custom_shark /usr/local/bin/shark")
        assert result == ("conda", "custom_shark")

    def test_extract_from_micromamba(self) -> None:
        result = extract_conda_env_from_command("/opt/conda/bin/micromamba run -n shark_env shark")
        assert result == ("/opt/conda/bin/micromamba", "shark_env")

    def test_returns_none_for_bare_command(self) -> None:
        assert extract_conda_env_from_command("shark") is None
        assert extract_conda_env_from_command("/usr/bin/shark") is None

    def test_returns_none_for_non_run_subcommand(self) -> None:
        assert extract_conda_env_from_command("mamba install -n shark_env shark") is None

    def test_returns_none_for_missing_env_flag(self) -> None:
        assert extract_conda_env_from_command("mamba run shark") is None

    def test_returns_none_for_too_few_tokens(self) -> None:
        assert extract_conda_env_from_command("mamba run") is None
        assert extract_conda_env_from_command("") is None

    def test_handles_shlex_parse_error(self) -> None:
        assert extract_conda_env_from_command('mamba run -n "unclosed quote') is None


class TestBuildSharkVersionCommand:
    def test_builds_argv_for_mamba(self) -> None:
        argv = build_shark_version_command("mamba run -n shark_env shark")
        assert argv == ["mamba", "list", "-n", "shark_env", "shark", "--json"]

    def test_builds_argv_for_conda(self) -> None:
        argv = build_shark_version_command("conda run --name my_env shark")
        assert argv == ["conda", "list", "-n", "my_env", "shark", "--json"]

    def test_returns_none_for_unsupported_prefix(self) -> None:
        assert build_shark_version_command("shark") is None
        assert build_shark_version_command("/usr/bin/shark") is None


class TestParseSharkCondaListJson:
    def test_parses_version_and_build_string(self) -> None:
        assert parse_shark_conda_list_json(PROBE_MAMBA_LIST_JSON) == "1.2.0+h077b44d_5"

    def test_parses_version_without_build_string(self) -> None:
        assert parse_shark_conda_list_json(PROBE_MAMBA_LIST_NO_BUILD_JSON) == "1.2.0"

    def test_returns_none_for_empty_list(self) -> None:
        assert parse_shark_conda_list_json(PROBE_EMPTY_LIST_JSON) is None

    def test_returns_none_for_invalid_json(self) -> None:
        assert parse_shark_conda_list_json("not valid json") is None

    def test_returns_none_for_empty_string(self) -> None:
        assert parse_shark_conda_list_json("") is None
        assert parse_shark_conda_list_json("   \n") is None

    def test_falls_back_to_first_entry_if_shark_name_not_matched(self) -> None:
        assert parse_shark_conda_list_json(PROBE_NON_SHARK_LIST_JSON) == "0.9.1+py310_0"


class TestParseSharkHelpOutput:
    def test_returns_none_for_standard_1_2_0_help(self) -> None:
        assert parse_shark_help_output(PROBE_HELP_OUTPUT) is None

    def test_extracts_version_when_present(self) -> None:
        synthetic = "shark version 1.2.0\nUsage: shark [options]"
        assert parse_shark_help_output(synthetic) == "1.2.0"

    def test_returns_none_for_empty_string(self) -> None:
        assert parse_shark_help_output("") is None


class TestResolveSharkVersionFromOutput:
    def test_resolves_valid_json_to_formatted_version(self) -> None:
        assert resolve_shark_version_from_output(PROBE_MAMBA_LIST_JSON) == "1.2.0+h077b44d_5"

    def test_falls_back_to_unknown_for_empty_or_failed_json(self) -> None:
        assert resolve_shark_version_from_output(PROBE_EMPTY_LIST_JSON) == "unknown"
        assert resolve_shark_version_from_output("") == "unknown"
        assert resolve_shark_version_from_output("invalid") == "unknown"
