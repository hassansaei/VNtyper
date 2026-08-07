"""The uniform-filtering branch may not run with an empty GG allowlist.

@hassansaei on #186: enabling ``use_uniform_filtering`` while
``motifs_for_alt_gg`` is ``[]`` "would delete every GG, including canonical
dupC" -- the pathogenic call this tool exists to find. The shipped config
sets the flag to false, so this guard cannot fire in production; it exists
so that flipping the flag fails loudly instead of silently deleting calls.
"""

import json
from pathlib import Path

import pandas as pd
import pytest

from vntyper.scripts.motif_processing import motif_correction_and_annotation

pytestmark = pytest.mark.unit

CONFIG_PATH = Path("vntyper/scripts/kestrel_config.json")


def _shipped_config():
    return json.loads(CONFIG_PATH.read_text())


def _one_row():
    # "Variant" is required downstream by motif_correction_and_annotation's keep_cols
    # (populated by mark_insertion/mark_deletion in the real pipeline); the drafted
    # fixture omitted it, which raised an unrelated pandas KeyError before the guard
    # under test ever ran. See task-4-report.md item 7.
    return pd.DataFrame(
        {
            "Motifs": ["X-5"],
            "Variant": ["67:C>CC"],
            "POS": [67],
            "REF": ["C"],
            "ALT": ["CC"],
            "Estimated_Depth_AlternateVariant": [50],
            "Estimated_Depth_Variant_ActiveRegion": [5000],
            "Depth_Score": [0.01],
            "Confidence": ["High_Precision"],
        }
    )


def test_the_shipped_config_leaves_uniform_filtering_off():
    """If this fails, the guard below is no longer inert and #186 is reopened."""
    config = _shipped_config()
    assert config["motif_filtering"]["use_uniform_filtering"] is False
    assert config["motif_filtering"]["motifs_for_alt_gg"] == []


def test_uniform_filtering_with_an_empty_gg_allowlist_is_refused():
    config = _shipped_config()
    config["motif_filtering"]["use_uniform_filtering"] = True

    with pytest.raises(ValueError, match="motifs_for_alt_gg"):
        motif_correction_and_annotation(_one_row(), pd.DataFrame({"Motif": [], "Motif_sequence": []}), config)


def test_uniform_filtering_with_a_populated_allowlist_is_permitted():
    """The guard blocks the armed combination only, not the branch itself."""
    config = _shipped_config()
    config["motif_filtering"]["use_uniform_filtering"] = True
    config["motif_filtering"]["motifs_for_alt_gg"] = ["X"]

    out = motif_correction_and_annotation(
        _one_row(), pd.DataFrame({"Motif": ["5"], "Motif_sequence": ["ACGT"]}), config
    )

    assert "motif_filter_pass" in out.columns
