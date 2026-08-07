"""The name a cohort sample is reported under when pseudonyms are asked for.

`vntyper/scripts/cohort_pseudonyms.py` turns a sample's identity into `<prefix><digest>`,
and `pseudonymization_table.tsv` beside the report is the only way back from it. These
tests are **characterisation** - they record what a pseudonymised cohort run produces
today, including the exact literals, because `sha256` has no per-process salt and a
recorded value is therefore the whole cross-process stability guarantee.

They moved here with the code when the pseudonym left `cohort_inputs.py`, which locates
samples on disk and has nothing to do with naming them. The *specification* the digest
width answers - that the map from patient to reported sample must be injective or the run
must stop - is in `test_cohort_identity.py`, which also exercises the validation of the
algorithm and width read out of `config["cohort"]["pseudonym"]`.
"""

from __future__ import annotations

import logging

import pytest

from vntyper.scripts.cohort_pseudonyms import (
    DEFAULT_PSEUDONYM_ALGORITHM,
    DEFAULT_PSEUDONYM_LENGTH,
    pseudonym_settings,
    pseudonymized_sample_name,
)

pytestmark = pytest.mark.unit


def test_a_pseudonym_is_the_prefix_and_twelve_hex_digits() -> None:
    """#206: five hex characters of MD5 became twelve of SHA-256, and the value moved."""
    assert pseudonymized_sample_name("anon_", "sample_one") == "anon_c788e939395d"


def test_the_same_sample_name_always_gets_the_same_pseudonym() -> None:
    """The mapping has to be stable so a cohort re-run stays comparable to its
    predecessor and to the pseudonymization table written beside it.

    Two calls in one interpreter, so on its own this shows only that the function is not
    stateful; it cannot see a mapping that varied between processes. What establishes
    cross-process stability is `test_a_pseudonym_is_the_prefix_and_twelve_hex_digits`
    above, which pins an exact literal - `sha256` is a fixed digest with no per-process
    salt, unlike `hash()`, so a recorded value is the whole guarantee.
    """
    assert pseudonymized_sample_name("x", "s1") == pseudonymized_sample_name("x", "s1")


def test_two_particular_sample_names_get_different_pseudonyms() -> None:
    """`s1` and `s2` do not collide. That is all this shows; injectivity over a realistic
    cohort is measured in `test_cohort_identity.py`."""
    assert pseudonymized_sample_name("x", "s1") != pseudonymized_sample_name("x", "s2")


def test_the_two_sample_names_that_used_to_share_a_pseudonym_no_longer_do() -> None:
    """#206; this replaces the characterisation of the same pair as a defect.

    The pseudonym was the first **five** hex characters of an MD5 - 20 bits, about a
    million values - so a cohort of a few thousand samples was already odds-on to contain
    a collision by the birthday bound. `sample_42` and `sample_919` are the first
    colliding pair among `sample_0`..`sample_19999`, and both landed on `168eb`.

    Two consequences, both real rather than theoretical, and both closed by the widening:

    * `aggregate_cohort` builds `sample_mapping[pseudonym] = original_sample`, so the
      second sample's entry **overwrote** the first in `pseudonymization_table.tsv` and
      one of the two originals could not be recovered from it;
    * both samples' rows were reported under the same `Sample` value, so the cohort's
      Kestrel table showed them as one sample with two calls.

    Widening changes every pseudonym in every existing report; that was the recorded
    reason to leave it, and this milestone made the decision the other way. The map is now
    injective-or-loud: `aggregate_cohort` refuses a collision rather than merging, which
    is what makes no width a silent risk. See `test_cohort_identity.py`.
    """
    assert pseudonymized_sample_name("anon_", "sample_42") == "anon_55b67a5ddb12"
    assert pseudonymized_sample_name("anon_", "sample_919") == "anon_177ad4b75617"


def test_the_prefix_may_be_any_value_the_cli_accepted() -> None:
    """`--pseudonymize` takes a string, and the CLI has passed `True` through in the
    past; the prefix is interpolated rather than concatenated so neither raises."""
    assert pseudonymized_sample_name(True, "s1").startswith("True")


# ---------------------------------------------------------------------------
# pseudonym_settings
# ---------------------------------------------------------------------------
#
# SPECIFICATION: `aggregate_cohort` used to read these two settings with
# `config.get("cohort", {}).get("pseudonym", {})`, which defends against a key
# being absent and not against it being present and null. `--config-path`
# replaces the whole configuration rather than merging it (AGENTS.md trap 2), so
# a hand-written document can carry `"cohort": null` - and `None.get` is an
# `AttributeError` that names neither the key nor the file, raised even by a run
# that never asked for pseudonyms. The end-to-end consequences (including the
# extraction directories it leaked) are in `test_cohort_identity.py`; what is
# pinned here is the read itself, level by level.


def test_the_settings_default_when_the_configuration_says_nothing() -> None:
    assert pseudonym_settings({}) == (DEFAULT_PSEUDONYM_ALGORITHM, DEFAULT_PSEUDONYM_LENGTH)


def test_the_settings_are_read_from_a_well_formed_configuration() -> None:
    config = {"cohort": {"pseudonym": {"algorithm": "sha1", "digest_characters": 4}}}

    assert pseudonym_settings(config) == ("sha1", 4)


def test_one_setting_may_be_given_without_the_other() -> None:
    """The block is a partial override, not a replacement: each key defaults on its own."""
    assert pseudonym_settings({"cohort": {"pseudonym": {"algorithm": "sha512"}}}) == (
        "sha512",
        DEFAULT_PSEUDONYM_LENGTH,
    )
    assert pseudonym_settings({"cohort": {"pseudonym": {"digest_characters": 20}}}) == (
        DEFAULT_PSEUDONYM_ALGORITHM,
        20,
    )


def test_an_absent_block_is_not_worth_a_warning(caplog: pytest.LogCaptureFixture) -> None:
    """Every configuration written before #206 omits it; that is not a misconfiguration."""
    with caplog.at_level(logging.WARNING, logger="vntyper.scripts.cohort_pseudonyms"):
        assert pseudonym_settings({"cohort": {}}) == (DEFAULT_PSEUDONYM_ALGORITHM, DEFAULT_PSEUDONYM_LENGTH)
        assert pseudonym_settings({}) == (DEFAULT_PSEUDONYM_ALGORITHM, DEFAULT_PSEUDONYM_LENGTH)

    # `caplog` collects every logger, this suite's own progress lines included, so the
    # module under test is named rather than asserting the whole record list is empty.
    assert [r for r in caplog.records if r.name == "vntyper.scripts.cohort_pseudonyms"] == []


@pytest.mark.parametrize(
    ("config", "named"),
    [
        ({"cohort": None}, "cohort"),
        ({"cohort": "sha256"}, "cohort"),
        ({"cohort": []}, "cohort"),
        ({"cohort": {"pseudonym": None}}, "cohort.pseudonym"),
        ({"cohort": {"pseudonym": 12}}, "cohort.pseudonym"),
        ({"cohort": {"pseudonym": ["sha256", 12]}}, "cohort.pseudonym"),
    ],
    ids=["null", "string", "list", "null-inner", "number-inner", "list-inner"],
)
def test_a_level_that_is_not_a_mapping_defaults_loudly(
    config: dict, named: str, caplog: pytest.LogCaptureFixture
) -> None:
    """Loudly, because the digest settings decide every pseudonym in the report.

    Refusing outright is the wrong answer at this level: the settings are only *used* when
    `--pseudonymize-samples` was given, and a malformed key nobody asked to use must not
    abort a forty-sample cohort. Refusing by name is what `pseudonymized_sample_name` does
    with a value that is well-formed JSON and still unusable.
    """
    with caplog.at_level(logging.WARNING, logger="vntyper.scripts.cohort_pseudonyms"):
        assert pseudonym_settings(config) == (DEFAULT_PSEUDONYM_ALGORITHM, DEFAULT_PSEUDONYM_LENGTH)

    assert any(named in record.getMessage() for record in caplog.records), (
        f"the fallback must name {named!r} so the operator can find the key at fault"
    )


@pytest.mark.parametrize("config", [None, "config.json", ["cohort"], 7])
def test_a_configuration_that_is_not_a_mapping_at_all_defaults_loudly(
    config: object, caplog: pytest.LogCaptureFixture
) -> None:
    """`--config-path` points at an arbitrary JSON file, and JSON's top level is not
    required to be an object."""
    with caplog.at_level(logging.WARNING, logger="vntyper.scripts.cohort_pseudonyms"):
        assert pseudonym_settings(config) == (DEFAULT_PSEUDONYM_ALGORITHM, DEFAULT_PSEUDONYM_LENGTH)

    assert any("rather than a mapping" in record.getMessage() for record in caplog.records)


def test_the_configured_values_are_passed_through_without_validation() -> None:
    """Validation belongs to `pseudonymized_sample_name`, which is also reached directly.

    Duplicating it here would report the same failure from two places and let them drift;
    what this read owes the caller is that a malformed *shape* cannot crash it, not that
    the values are usable.
    """
    config = {"cohort": {"pseudonym": {"algorithm": "not-a-hash", "digest_characters": "twelve"}}}

    algorithm, length = pseudonym_settings(config)

    assert (algorithm, length) == ("not-a-hash", "twelve")
    with pytest.raises(ValueError, match="not-a-hash"):
        pseudonymized_sample_name("anon_", "s1", algorithm=algorithm, length=length)
