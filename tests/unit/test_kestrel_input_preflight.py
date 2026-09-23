"""The Kestrel file preflight: fail before any work, and say how to fix it (#338).

Config paths resolve against the CWD (AGENTS.md trap 7), so a source checkout without
``vntyper install-references``, or a job script that ``cd``s elsewhere, has no motif
reference. In fast mode with a CRAM nothing else under ``reference/`` is read, so before
this check the first sign was Kestrel's silent exit 0 at the end of every sample's run.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from vntyper.scripts.failure_help import TROUBLESHOOTING_URL
from vntyper.scripts.kestrel_counting import DEFAULT_KANALYZE_PATH
from vntyper.scripts.kestrel_input_preflight import (
    REFERENCE_KEYS,
    KestrelInput,
    check_kestrel_inputs,
    describe_missing_inputs,
    missing_kestrel_inputs,
    required_kestrel_inputs,
)

pytestmark = pytest.mark.unit

SPLIT = {"kestrel_settings": {"split_counting": True}}
INTERNAL = {"kestrel_settings": {"split_counting": False}}


def _config() -> dict:
    return {
        "tools": {
            "kestrel": "vntyper/dependencies/kestrel/kestrel.jar",
            "kanalyze": "vntyper/dependencies/kestrel/kanalyze.jar",
        },
        "reference_data": {
            "muc1_reference_vntr": "reference/All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa",
            "muc1_motifs_rev_com": "reference/MUC1_motifs_Rev_com.fa",
        },
    }


def _install(root: Path, config: dict) -> None:
    """Create every file ``config`` declares under ``root``."""
    for section in ("tools", "reference_data"):
        for value in config[section].values():
            path = root / value
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("x", encoding="utf-8")


# ---------------------------------------------------------------------------
# required_kestrel_inputs
# ---------------------------------------------------------------------------


def test_split_counting_requires_both_jars_and_both_references(tmp_path):
    keys = [item.key for item in required_kestrel_inputs(_config(), tmp_path, split=True)]
    assert keys == [
        "tools.kestrel",
        "tools.kanalyze",
        "reference_data.muc1_reference_vntr",
        "reference_data.muc1_motifs_rev_com",
    ]


def test_internal_counting_does_not_require_kanalyze(tmp_path):
    keys = [item.key for item in required_kestrel_inputs(_config(), tmp_path, split=False)]
    assert "tools.kanalyze" not in keys
    assert len(keys) == 1 + len(REFERENCE_KEYS)


def test_relative_paths_resolve_against_the_project_root(tmp_path):
    item = required_kestrel_inputs(_config(), tmp_path, split=False)[1]
    assert item.resolved == (tmp_path / "reference/All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa").resolve()
    assert item.configured == "reference/All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa"


def test_an_absolute_path_is_used_as_given(tmp_path):
    config = _config()
    absolute = tmp_path / "elsewhere" / "motifs.fa"
    config["reference_data"]["muc1_reference_vntr"] = str(absolute)
    item = required_kestrel_inputs(config, tmp_path / "cwd", split=False)[1]
    assert item.resolved == absolute.resolve()


def test_an_undeclared_kanalyze_falls_back_to_the_shipped_default(tmp_path):
    """The stage uses the same default, so the preflight must check the same file."""
    config = _config()
    del config["tools"]["kanalyze"]
    kanalyze = required_kestrel_inputs(config, tmp_path, split=True)[1]
    assert kanalyze.configured == DEFAULT_KANALYZE_PATH


@pytest.mark.parametrize("value", [None, "", 7])
def test_an_unset_or_non_string_value_is_reported_not_skipped(tmp_path, value):
    config = _config()
    config["reference_data"]["muc1_reference_vntr"] = value
    item = required_kestrel_inputs(config, tmp_path, split=False)[1]
    assert item.configured == ""
    assert missing_kestrel_inputs([item]) == [item]


# ---------------------------------------------------------------------------
# missing_kestrel_inputs
# ---------------------------------------------------------------------------


def test_nothing_is_missing_when_everything_is_installed(tmp_path):
    config = _config()
    _install(tmp_path, config)
    assert missing_kestrel_inputs(required_kestrel_inputs(config, tmp_path, split=True)) == []


def test_the_missing_reference_is_the_one_reported(tmp_path):
    """The reported bug: the JARs ship in the checkout, the reference bundle does not."""
    config = _config()
    _install(tmp_path, config)
    (tmp_path / config["reference_data"]["muc1_reference_vntr"]).unlink()
    missing = missing_kestrel_inputs(required_kestrel_inputs(config, tmp_path, split=True))
    assert [item.key for item in missing] == ["reference_data.muc1_reference_vntr"]


def test_a_directory_is_not_a_reference(tmp_path):
    config = _config()
    _install(tmp_path, config)
    target = tmp_path / config["reference_data"]["muc1_motifs_rev_com"]
    target.unlink()
    target.mkdir()
    missing = missing_kestrel_inputs(required_kestrel_inputs(config, tmp_path, split=False))
    assert [item.key for item in missing] == ["reference_data.muc1_motifs_rev_com"]


# ---------------------------------------------------------------------------
# describe_missing_inputs
# ---------------------------------------------------------------------------


def test_a_missing_reference_message_names_path_cwd_and_install_command(tmp_path):
    item = KestrelInput("reference_data.muc1_reference_vntr", "reference/m.fa", tmp_path / "reference/m.fa")
    msg = describe_missing_inputs([item], tmp_path)
    assert f"current working directory ({tmp_path})" in msg
    assert f"reference_data.muc1_reference_vntr = 'reference/m.fa' -> {tmp_path / 'reference/m.fa'} (not found)" in msg
    assert "vntyper install-references --output-dir reference" in msg
    assert "vntyper/dependencies/kestrel" not in msg
    assert msg.endswith(f"Help: {TROUBLESHOOTING_URL}")


def test_a_missing_jar_message_points_at_the_source_tree(tmp_path):
    item = KestrelInput("tools.kestrel", "k.jar", tmp_path / "k.jar")
    msg = describe_missing_inputs([item], tmp_path)
    assert "vntyper/dependencies/kestrel/" in msg
    assert "install-references" not in msg


def test_both_kinds_missing_gives_both_fixes(tmp_path):
    items = [
        KestrelInput("tools.kanalyze", "a.jar", tmp_path / "a.jar"),
        KestrelInput("reference_data.muc1_motifs_rev_com", "b.fa", tmp_path / "b.fa"),
    ]
    msg = describe_missing_inputs(items, tmp_path)
    assert "install-references" in msg and "vntyper/dependencies/kestrel/" in msg


def test_an_unset_key_says_so(tmp_path):
    msg = describe_missing_inputs([KestrelInput("reference_data.muc1_reference_vntr", "", tmp_path)], tmp_path)
    assert "reference_data.muc1_reference_vntr is not set in the config" in msg


def test_describing_nothing_is_a_caller_bug(tmp_path):
    with pytest.raises(ValueError, match="at least one missing input"):
        describe_missing_inputs([], tmp_path)


# ---------------------------------------------------------------------------
# check_kestrel_inputs
# ---------------------------------------------------------------------------


def test_the_check_passes_on_a_complete_install(tmp_path):
    config = _config()
    _install(tmp_path, config)
    check_kestrel_inputs(config, tmp_path, runtime_component=SPLIT)


def test_the_check_raises_and_logs_on_a_missing_reference(tmp_path, caplog):
    config = _config()
    _install(tmp_path, config)
    (tmp_path / config["reference_data"]["muc1_reference_vntr"]).unlink()
    with caplog.at_level(logging.ERROR), pytest.raises(ValueError, match="Kestrel cannot run") as excinfo:
        check_kestrel_inputs(config, tmp_path, runtime_component=SPLIT)
    assert "muc1_reference_vntr" in str(excinfo.value)
    assert str(excinfo.value) in caplog.text


def test_split_mode_from_the_runtime_decides_whether_kanalyze_is_required(tmp_path):
    config = _config()
    _install(tmp_path, config)
    (tmp_path / config["tools"]["kanalyze"]).unlink()
    check_kestrel_inputs(config, tmp_path, runtime_component=INTERNAL)
    with pytest.raises(ValueError, match="tools.kanalyze"):
        check_kestrel_inputs(config, tmp_path, runtime_component=SPLIT)


def test_the_packaged_runtime_splits_counting(tmp_path):
    """With no runtime given, the shipped sidecar applies, and it splits (#262)."""
    config = _config()
    _install(tmp_path, config)
    (tmp_path / config["tools"]["kanalyze"]).unlink()
    with pytest.raises(ValueError, match="tools.kanalyze"):
        check_kestrel_inputs(config, tmp_path)
