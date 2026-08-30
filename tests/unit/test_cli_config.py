"""CLI configuration loading contracts needed by exact fastp cutoffs."""

from __future__ import annotations

import copy
import json
from decimal import Decimal
from pathlib import Path

import pytest

import vntyper
from vntyper.cli import load_config
from vntyper.scripts.fastp_cutoffs import build_fastp_cutoffs

pytestmark = pytest.mark.unit


def _config_with_exact_q20(tmp_path: Path) -> Path:
    """Write the complete shipped config with one adversarial exact JSON decimal."""
    shipped_path = Path(vntyper.__file__).resolve().parent / "config.json"
    shipped = shipped_path.read_text(encoding="utf-8")
    original = '"q20_rate": 0.8'
    assert shipped.count(original) == 1
    custom = shipped.replace(original, '"q20_rate": 0.60044999999999999')
    config_path = tmp_path / "config.json"
    config_path.write_text(custom, encoding="utf-8")
    return config_path


def test_load_config_preserves_exact_fastp_threshold_provenance_without_changing_public_values(
    tmp_path: Path,
) -> None:
    """Binary float parsing cannot move a JSON cutoff across its visible half-tie."""
    config = load_config(_config_with_exact_q20(tmp_path))

    assert type(config) is dict
    thresholds = config["thresholds"]
    assert isinstance(thresholds, dict)
    assert type(thresholds["q20_rate"]) is float
    assert thresholds["q20_rate"] == 0.60045
    assert json.loads(json.dumps(config))["thresholds"]["q20_rate"] == thresholds["q20_rate"]

    cutoff = build_fastp_cutoffs(thresholds).q20_rate
    assert cutoff.value == Decimal("0.6004")
    assert cutoff.label == "60.04%"


def test_exact_fastp_threshold_provenance_survives_deepcopy(tmp_path: Path) -> None:
    """Pipeline copies must not silently fall back to the rounded binary float."""
    copied = copy.deepcopy(load_config(_config_with_exact_q20(tmp_path)))

    cutoff = build_fastp_cutoffs(copied["thresholds"]).q20_rate
    assert cutoff.value == Decimal("0.6004")
    assert cutoff.label == "60.04%"


def test_programmatic_threshold_mutation_supersedes_loaded_json_provenance(tmp_path: Path) -> None:
    """A caller's current mapping value wins over stale source-side metadata."""
    config = load_config(_config_with_exact_q20(tmp_path))
    config["thresholds"]["q20_rate"] = Decimal("0.81234999999999999")

    cutoff = build_fastp_cutoffs(config["thresholds"]).q20_rate
    assert cutoff.value == Decimal("0.8123")
    assert cutoff.label == "81.23%"


@pytest.mark.parametrize("mutation_api", ["assignment", "update"])
def test_same_float_mutation_invalidates_loaded_json_provenance(tmp_path: Path, mutation_api: str) -> None:
    """Even an equal-looking float assignment is a new programmatic source value."""
    config = load_config(_config_with_exact_q20(tmp_path))
    thresholds = config["thresholds"]
    if mutation_api == "assignment":
        thresholds["q20_rate"] = 0.60045
    else:
        thresholds.update({"q20_rate": 0.60045})

    cutoff = build_fastp_cutoffs(thresholds).q20_rate
    assert cutoff.value == Decimal("0.6005")
    assert cutoff.label == "60.05%"


def test_replacing_the_threshold_mapping_keeps_plain_dict_compatibility(tmp_path: Path) -> None:
    """Programmatic configs need no provenance wrapper to retain existing behavior."""
    config = load_config(_config_with_exact_q20(tmp_path))
    config["thresholds"] = {
        "duplication_rate": 0.1,
        "q20_rate": 0.81234,
        "q30_rate": 0.7,
        "passed_filter_reads_rate": 0.8,
    }

    cutoff = build_fastp_cutoffs(config["thresholds"]).q20_rate
    assert cutoff.value == Decimal("0.8123")
    assert cutoff.label == "81.23%"
