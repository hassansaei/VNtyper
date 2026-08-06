# tests/unit/test_frameshift_convention_parity.py

"""One frameshift convention, two modules -- asserted to stay in agreement.

@hassansaei on #182: "Yes, keep the same 3n+1 / 3n+2 rule for adVNTR as for Kestrel
(#181). This is intentional shared convention, not something to relax independently.
[...] We should keep the filtering logic harmonized between Kestrel and adVNTR, as
already implemented since v1.3."

Kestrel encodes the rule as a pair of pandas conditions in
``scoring.extract_frameshifts``; adVNTR encodes it as two numpy arrays built inside
``advntr_processing_del``/``advntr_processing_ins`` from the shipped
``advntr_config.json``. Nothing previously asserted the two agree, so either could have
drifted alone and nothing would have noticed. This file is that assertion.

The tests below drive each module's real code path rather than reconstructing its
arithmetic locally: a local reconstruction of ``ins_frame``/``del_frame`` (e.g.
``np.arange(100) * 3 + 1``) would only prove properties of its own construction, not of
``advntr_genotyping`` -- it could not fail no matter what the module actually does.
"""

import pandas as pd
import pytest

from vntyper.modules.advntr import advntr_genotyping
from vntyper.scripts.scoring import extract_frameshifts

pytestmark = pytest.mark.unit

#: The multiplier the shared convention requires: both modules gate on ``3n+1``/``3n+2``.
FRAMESHIFT_MULTIPLIER = 3


def _kestrel_accepts(direction, amount):
    """True when Kestrel's real ``extract_frameshifts`` marks this row a valid frameshift.

    ``amount`` is passed as the already-mod-3-reduced residue: production
    ``scoring.split_frame_score`` computes ``frameshift_amount = abs(alt_len - ref_len)
    % 3`` (``vntyper/scripts/scoring.py:114``) before ``extract_frameshifts`` ever sees
    it, so ``extract_frameshifts`` itself only ever compares a residue in ``{0, 1, 2}``
    against the literals 1 (insertion) / 2 (deletion).
    """
    out = extract_frameshifts(pd.DataFrame({"direction": [direction], "frameshift_amount": [amount]}))
    return bool(out["is_valid_frameshift"].iloc[0])


def _advntr_frame(state: str) -> pd.DataFrame:
    """A single-row adVNTR frame in the shape ``process_advntr_output`` hands over."""
    return pd.DataFrame(
        {
            "VID": [25561],
            "State": [state],
            "NumberOfSupportingReads": [11],
            "MeanCoverage": [153.98],
            "Pvalue": [0.0001],
        }
    )


def _advntr_accepts_deletion(length: int) -> bool:
    """True when adVNTR's real ``advntr_processing_del`` keeps a pure ``length``-bp deletion.

    Drives the actual function -- not a local copy of ``del_frame`` -- on a pure
    deletion, where ``Deletion_length`` is always >= 1 (so that guard never decides the
    outcome) and the answer is therefore governed entirely by the module's real
    ``del_frame`` membership check.
    """
    state = "&".join(f"D{8 + offset}_2" for offset in range(length))
    return len(advntr_genotyping.advntr_processing_del(_advntr_frame(state))) == 1


def _advntr_accepts_insertion(length: int) -> bool:
    """True when adVNTR's real ``advntr_processing_ins`` keeps a pure ``length``-bp insertion.

    Mirror of :func:`_advntr_accepts_deletion`: ``Insertion_len`` is always >= 1 for a
    pure insertion, so the outcome is governed entirely by the module's real
    ``ins_frame`` membership check.
    """
    state = f"I10_2_A_LEN{length}"
    return len(advntr_genotyping.advntr_processing_ins(_advntr_frame(state))) == 1


def test_kestrel_accepts_insertions_at_3n_plus_1_and_deletions_at_3n_plus_2():
    """Anchors the comparison in Kestrel's own, real behaviour (#181)."""
    assert _kestrel_accepts(1, 1) is True
    assert _kestrel_accepts(1, 2) is False
    assert _kestrel_accepts(-1, 2) is True
    assert _kestrel_accepts(-1, 1) is False


@pytest.mark.parametrize("length", range(1, 10))
def test_advntr_real_filtering_behaviour_matches_kestrel_for_every_residue(length):
    """adVNTR's actual accept/reject decision must match Kestrel's for the same residue.

    Specification (#182, decided 2026-08-06): the two modules' filters are required to
    stay harmonized. This calls ``advntr_processing_del``/``advntr_processing_ins`` --
    the real code path -- for lengths 1-9 (covering every residue mod 3 more than once)
    and compares each outcome to ``extract_frameshifts``'s decision for the matching
    ``(direction, residue)``. Because both sides call real production code, this fails
    if either module's multiplier or offset drifts from the other's -- see the report
    for an induced-failure demonstration.
    """
    residue = length % 3

    assert _advntr_accepts_insertion(length) == _kestrel_accepts(1, residue)
    assert _advntr_accepts_deletion(length) == _kestrel_accepts(-1, residue)


def test_the_shipped_advntr_config_uses_the_multiplier_the_convention_requires():
    """A multiplier other than 3 would silently decouple the two modules.

    Reads ``advntr_genotyping.advntr_config`` -- the module-level global built from the
    shipped ``advntr_config.json`` at import time -- rather than restating the number,
    so a config edit is caught here even though ``advntr_processing_del``/``_ins``
    actually read the *derived* ``advntr_settings`` global at call time (see
    ``tests/unit/test_advntr_frameshift_filter.py::TestSettingsComeFromTheDerivedGlobal``
    for that distinction).
    """
    config = advntr_genotyping.advntr_config["advntr_settings"]
    assert config["frameshift_multiplier"] == FRAMESHIFT_MULTIPLIER
