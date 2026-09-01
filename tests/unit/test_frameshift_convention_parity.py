# tests/unit/test_frameshift_convention_parity.py

"""One frameshift convention, two modules -- asserted to stay in agreement.

@hassansaei on #182: "Yes, keep the same 3n+1 / 3n+2 rule for adVNTR as for Kestrel
(#181). This is intentional shared convention, not something to relax independently.
[...] We should keep the filtering logic harmonized between Kestrel and adVNTR, as
already implemented since v1.3."

The convention, stated once
---------------------------
Let ``Delta`` be the **signed** net change in bases. A row enters the pathogenic
ADTKD-MUC1 frame exactly when ``Delta % 3 == 1`` -- Python's ``%`` returns a non-negative
residue, so this covers a net insertion of ``3n+1`` bases (``+1 % 3 == 1``) and a net
deletion of ``3n+2`` bases (``-2 % 3 == 1``) and nothing else.

The two modules reach it from different inputs:

* Kestrel derives ``Delta`` from the REF/ALT lengths in ``scoring``:
  ``split_frame_score`` splits it into ``direction = sign(Delta)`` and
  ``frameshift_amount = |Delta| % 3``, and ``extract_frameshifts`` accepts
  ``direction > 0 & amount == 1`` or ``direction < 0 & amount == 2``.
* adVNTR derives ``Delta = Insertion_len - Deletion_length`` from the ``State`` string in
  ``advntr_genotyping.derive_indel_columns``, and its two arms test the sign of that
  against magnitude series built from the shipped ``advntr_config.json``.

Nothing previously asserted the two agree, so either could have drifted alone and nothing
would have noticed. This file is that assertion.

Mixed states are covered, not only pure ones
--------------------------------------------
An earlier revision compared only *pure* insertions and *pure* deletions, where the sign
of ``Delta`` is never in doubt. That is precisely the region in which adVNTR's filter was
correct: it discarded the sign (it tested ``|Delta|``) and guarded only on "names at least
one deletion" / "names at least one insertion", which a **mixed** state satisfies on both
arms. So the two implementations agreed on every case the parity test looked at while
disagreeing on the ones it did not. The mixed grid below closes that gap: it asserts the
two modules agree on the *signed* rule, which is the rule they are actually required to
share.

The tests drive each module's real code path rather than reconstructing its arithmetic
locally: a local reconstruction of the magnitude series (e.g. ``np.arange(100) * 3 + 1``)
would only prove properties of its own construction, not of ``advntr_genotyping`` -- it
could not fail no matter what the module actually does.
"""

import pandas as pd
import pytest

from vntyper.modules.advntr import advntr_genotyping
from vntyper.scripts.scoring import (
    extract_frameshifts,
    split_depth_and_calculate_frame_score,
    split_frame_score,
)

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


def _kestrel_accepts_net_change(delta: int) -> bool:
    """True when Kestrel's real pipeline calls a net change of ``delta`` bases pathogenic.

    Drives the whole production chain -- ``split_depth_and_calculate_frame_score`` ->
    ``split_frame_score`` -> ``extract_frameshifts`` -- from a REF/ALT pair whose length
    difference is ``delta``, so ``direction`` and ``frameshift_amount`` are computed by
    ``scoring`` rather than restated here. Kestrel's own input is a REF/ALT pair, and a
    mixed indel reaches it as exactly this: a single net length change.

    Args:
        delta (int): ``alt_len - ref_len``, the signed net change in bases.

    Returns:
        bool: ``is_valid_frameshift`` for that row.
    """
    ref = "A" * (1 + max(0, -delta))
    alt = "A" * (1 + max(0, delta))
    assert len(alt) - len(ref) == delta, "the constructed REF/ALT pair does not carry the intended delta"

    frame = pd.DataFrame({"Sample": ["Del:10:100"], "REF": [ref], "ALT": [alt]})
    scored = extract_frameshifts(split_frame_score(split_depth_and_calculate_frame_score(frame)))
    return bool(scored["is_valid_frameshift"].iloc[0])


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


def _advntr_reports(state: str) -> bool:
    """True when ``state`` reaches ``output_adVNTR_result.tsv``.

    ``process_advntr_output`` concatenates the deletion half and the insertion half, so a
    state is reported when *either* real filter keeps it. This drives both, which is what
    makes the comparison below one about the module's decision rather than about one arm
    of it -- and it is the only formulation that can catch an arm claiming a row of the
    wrong sign, since which arm kept a row is not a column of the result file.
    """
    frame = _advntr_frame(state)
    kept_del = len(advntr_genotyping.advntr_processing_del(frame.copy()))
    kept_ins = len(advntr_genotyping.advntr_processing_ins(frame.copy()))
    assert kept_del + kept_ins <= 1, f"{state} was claimed by both arms"
    return bool(kept_del + kept_ins)


def _advntr_net_change(state: str) -> int:
    """``Insertion_len - Deletion_length`` for ``state``, from the production derivation."""
    return int(advntr_genotyping.derive_indel_columns(_advntr_frame(state))["Net_indel_length"].iloc[0])


def _pure_deletion(length: int) -> str:
    """An adVNTR state naming ``length`` consecutive single-base deletions."""
    return "&".join(f"D{8 + offset}_2" for offset in range(length))


def _pure_insertion(length: int) -> str:
    """An adVNTR state naming a single insertion of ``length`` bases."""
    return f"I10_2_A_LEN{length}"


def _mixed(inserted: int, deleted: int) -> str:
    """A state naming ``deleted`` single-base deletions and one ``inserted``-bp insertion."""
    return "&".join([f"I9_2_A_LEN{inserted}"] + [f"D{50 + offset}_2" for offset in range(deleted)])


def _advntr_accepts_deletion(length: int) -> bool:
    """True when adVNTR's real ``advntr_processing_del`` keeps a pure ``length``-bp deletion.

    Drives the actual function -- not a local copy of the magnitude series -- on a pure
    deletion, where ``Net_indel_length`` is always negative (so the arm's sign test never
    decides the outcome) and the answer is therefore governed entirely by the module's
    real ``3n+2`` membership check.
    """
    return len(advntr_genotyping.advntr_processing_del(_advntr_frame(_pure_deletion(length)))) == 1


def _advntr_accepts_insertion(length: int) -> bool:
    """True when adVNTR's real ``advntr_processing_ins`` keeps a pure ``length``-bp insertion.

    Mirror of :func:`_advntr_accepts_deletion`: ``Net_indel_length`` is always positive
    for a pure insertion, so the outcome is governed entirely by the module's real
    ``3n+1`` membership check.
    """
    return len(advntr_genotyping.advntr_processing_ins(_advntr_frame(_pure_insertion(length)))) == 1


def test_kestrel_accepts_insertions_at_3n_plus_1_and_deletions_at_3n_plus_2():
    """Anchors the comparison in Kestrel's own, real behaviour (#181)."""
    assert _kestrel_accepts(1, 1) is True
    assert _kestrel_accepts(1, 2) is False
    assert _kestrel_accepts(-1, 2) is True
    assert _kestrel_accepts(-1, 1) is False


@pytest.mark.parametrize("delta", range(-9, 10))
def test_kestrels_own_chain_implements_the_signed_rule(delta):
    """Kestrel's real REF/ALT chain accepts exactly ``Delta % 3 == 1``.

    Anchors :func:`_kestrel_accepts_net_change` -- the oracle the mixed comparison below
    uses -- against the closed form, so a change in ``scoring`` that broke the signed rule
    would fail here as well as in the parity comparison, and it would be obvious which
    module moved.
    """
    assert _kestrel_accepts_net_change(delta) is (delta % 3 == 1)


@pytest.mark.parametrize("length", range(1, 10))
def test_advntr_real_filtering_behaviour_matches_kestrel_for_every_residue(length):
    """adVNTR's actual accept/reject decision must match Kestrel's for the same residue.

    Specification (#182, decided 2026-08-06): the two modules' filters are required to
    stay harmonized. This calls ``advntr_processing_del``/``advntr_processing_ins`` --
    the real code path -- for pure states of length 1-9 (covering every residue mod 3 more
    than once) and compares each outcome to ``extract_frameshifts``'s decision for the
    matching ``(direction, residue)``. Because both sides call real production code, this
    fails if either module's multiplier or offset drifts from the other's -- see the
    report for an induced-failure demonstration.
    """
    residue = length % 3

    assert _advntr_accepts_insertion(length) == _kestrel_accepts(1, residue)
    assert _advntr_accepts_deletion(length) == _kestrel_accepts(-1, residue)


@pytest.mark.parametrize("inserted", range(1, 8))
@pytest.mark.parametrize("deleted", range(1, 8))
def test_advntr_and_kestrel_agree_on_mixed_states_too(inserted, deleted):
    """The signed rule, asserted over a grid where only the *net* change decides.

    Every cell names both an insertion and a deletion, so adVNTR's two presence guards --
    ``Deletion_length >= 1`` and ``Insertion_len >= 1`` -- are both satisfied and neither
    can decide the outcome. What is left is the sign and the magnitude of
    ``Net_indel_length``, and Kestrel's answer for the same net change is the oracle.

    This is the region the pure-state comparison could not see. Before the sign repair,
    ``I9_2_A_LEN3&D50_2`` (``Delta = +2``) was reported by adVNTR's *deletion* arm while
    Kestrel rejects a net ``+2``; the grid below contains that cell at
    ``(inserted, deleted) == (3, 1)``.
    """
    state = _mixed(inserted, deleted)
    delta = inserted - deleted

    assert _advntr_net_change(state) == delta, "the adVNTR derivation does not carry the intended delta"
    assert _advntr_reports(state) == _kestrel_accepts_net_change(delta)


def test_the_two_modules_agree_on_the_named_regression_states():
    """The three states the review measured, by name, on both implementations.

    Kept separate from the grid so they cannot fall out of a parametrisation range, and so
    the arithmetic that makes each one non-pathogenic is written down beside it.
    """
    named = {
        # 3 inserted, 1 deleted -> Delta = +2; +2 % 3 == 2. Was reported via the deletion
        # arm because |+2| = 2 is in the 3n+2 series.
        "I9_2_A_LEN3&D50_2": 2,
        # 12 inserted, 1 deleted -> Delta = +11; 11 % 3 == 2. Was reported via the
        # deletion arm because |+11| = 11 is in the 3n+2 series.
        "D49_2&I49_2_A_LEN12": 11,
        # 1 inserted, 2 deleted -> Delta = -1; -1 % 3 == 2. Was reported via the insertion
        # arm because |-1| = 1 is in the 3n+1 series.
        "D8_2&D9_2&I9_2_A_LEN1": -1,
        # The control: 3 inserted, 2 deleted -> Delta = +1, the pathogenic frame. Reported
        # before and after, so the repair narrows the filter without closing it.
        "I9_2_A_LEN3&D50_2&D51_2": 1,
    }

    for state, delta in named.items():
        assert _advntr_net_change(state) == delta, state
        assert _advntr_reports(state) == _kestrel_accepts_net_change(delta), state
        assert _advntr_reports(state) is (delta % 3 == 1), state


def test_the_shipped_advntr_config_uses_the_multiplier_the_convention_requires():
    """A multiplier other than 3 would silently decouple the two modules.

    Reads ``advntr_genotyping.advntr_config`` -- the module-level global built from the
    shipped ``advntr_config.json`` at import time -- rather than restating the number,
    so a config edit is caught here even though ``advntr_processing_del``/``_ins``
    actually read the *derived* ``advntr_settings`` global at call time (see
    ``tests/unit/test_advntr_frameshift_filter.py::TestSettingsComeFromTheDerivedGlobal``
    for that distinction).
    """
    config = advntr_genotyping.load_advntr_config()["advntr_settings"]
    assert config["frameshift_multiplier"] == FRAMESHIFT_MULTIPLIER
