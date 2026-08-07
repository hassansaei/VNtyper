"""The digest a cohort pseudonym is built from, and the configuration that chooses it (#206).

**SPECIFICATION.** The pseudonym was ``<prefix>`` plus the first **five** hex characters of
the MD5 of the sample name - 20 bits. ``sample_mapping`` is keyed on the pseudonym, so a
collision silently overwrote the earlier original name: two patients' rows became
indistinguishable across every export, ``sample_categories()`` counted them as one sample,
and ``pseudonymization_table.tsv`` mapped the shared pseudonym to whichever original was
written last. Birthday probability of at least one collision is ``1 - exp(-n(n-1)/2**21)``:
~37.9% at 1,000 samples. The verified first collision in ``sample_0..sample_19999`` was
``sample_42`` and ``sample_919``, both MD5-prefixing to ``168eb``; that exact pair is the
probe below.

Both settings - the algorithm and the width - are read out of
``config["cohort"]["pseudonym"]``, and ``--config-path`` replaces the whole configuration
rather than merging it (AGENTS.md trap 2). So every value here can be any JSON value at
all, a malformed *shape* must not abort a cohort that never asked for pseudonyms, and a
value that is well-formed JSON and still unusable must be refused **by name** rather than
silently substituted: a silent substitution changes every pseudonym in a report without
saying so.

Split out of ``test_cohort_identity.py``, which reached 1,210 lines. What is left there is
the identity the digest is taken *of*; the ZIP half is in ``test_cohort_zip_identity.py``
and the de-duplication half in ``test_cohort_deduplication.py``.
"""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import pytest

from tests.support.cohort_samples import sample_on_disk
from vntyper.cli import load_config
from vntyper.scripts import cohort_summary
from vntyper.scripts.cohort_pseudonyms import (
    DEFAULT_PSEUDONYM_ALGORITHM,
    DEFAULT_PSEUDONYM_LENGTH,
    pseudonymized_sample_name,
)

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_the_verified_md5_collision_no_longer_collides() -> None:
    """`sample_42` and `sample_919` were the first colliding pair under the old scheme."""
    assert pseudonymized_sample_name("anon_", "sample_42") != pseudonymized_sample_name("anon_", "sample_919")


def test_pseudonyms_are_injective_over_twenty_thousand_names() -> None:
    """The exact probe that found the bug, at the configured width."""
    names = [f"sample_{index}" for index in range(20000)]
    pseudonyms = {pseudonymized_sample_name("anon_", name) for name in names}

    assert len(pseudonyms) == len(names)


def test_the_default_digest_is_twelve_characters_of_sha256() -> None:
    """The literal is recorded rather than recomputed with the code under test.

    ``sha256`` has no per-process salt, so a recorded value is the whole cross-process
    stability guarantee - the same reason the MD5 literal was recorded before it.

    A ``test_the_pseudonym_is_stable_across_calls`` stood beside this and compared two
    identical calls to each other. It duplicated the same idea in
    ``test_cohort_pseudonyms.py``, it passed on ``main``, and any deterministic but wrong
    digest satisfied it - so the recorded literals here subsume it entirely and it is gone.
    """
    assert DEFAULT_PSEUDONYM_ALGORITHM == "sha256"
    assert DEFAULT_PSEUDONYM_LENGTH == 12
    assert pseudonymized_sample_name("anon_", "s1") == "anon_e8bc163c82ee"
    assert pseudonymized_sample_name("anon_", "sample_one") == "anon_c788e939395d"


def test_the_digest_width_is_honoured() -> None:
    assert pseudonymized_sample_name("", "s1", length=5) == "e8bc1"
    assert len(pseudonymized_sample_name("", "s1", length=5)) == 5
    assert len(pseudonymized_sample_name("", "s1", length=64)) == 64


def test_a_non_string_prefix_is_interpolated_rather_than_concatenated() -> None:
    """Preserved behaviour: ``--pseudonymize-samples`` may carry a non-string."""
    assert pseudonymized_sample_name(True, "s1") == "Truee8bc163c82ee"


def test_an_unavailable_algorithm_is_refused_by_name() -> None:
    """A silent fallback would change every pseudonym in the report without saying so."""
    with pytest.raises(ValueError, match="not-a-hash"):
        pseudonymized_sample_name("anon_", "s1", algorithm="not-a-hash")


@pytest.mark.parametrize(
    "algorithm",
    [["sha256"], {"algorithm": "sha256"}, 256, None, "not-a-hash"],
    ids=["list", "dict", "number", "null", "unknown-name"],
)
def test_an_algorithm_that_is_not_a_usable_name_is_refused_the_same_way(algorithm: object) -> None:
    """`algorithm` comes out of a JSON file, so it can be any JSON value at all.

    `algorithm not in hashlib.algorithms_available` is a `set` membership test, and a
    `list` or a `dict` is unhashable: it raised `TypeError: unhashable type: 'list'` from
    inside the guard that exists to produce a `ValueError`. The documented `Raises:` was
    wrong for those two shapes, the repository's `logger.error` + `raise` convention was
    bypassed, and the message named neither the setting nor the configuration key. The
    width beside it was already checked by type first; the algorithm now is too.
    """
    with pytest.raises(ValueError, match="Unknown pseudonym digest algorithm"):
        pseudonymized_sample_name("anon_", "s1", algorithm=algorithm)  # type: ignore[arg-type]


@pytest.mark.parametrize("length", [0, -1, 2.5, "12", None])
def test_a_digest_length_that_is_not_a_positive_integer_is_refused(length: object) -> None:
    """``digest_characters`` comes out of a JSON file, so it can be anything at all."""
    with pytest.raises(ValueError, match="positive integer"):
        pseudonymized_sample_name("anon_", "s1", length=length)  # type: ignore[arg-type]


def test_a_digest_length_wider_than_the_digest_is_refused() -> None:
    """`sha256` produces 64 hex characters; asking for 65 would silently give 64."""
    with pytest.raises(ValueError, match="65"):
        pseudonymized_sample_name("anon_", "s1", length=65)


def test_a_boolean_length_is_refused_rather_than_read_as_one() -> None:
    """``True`` is an ``int`` in Python, and a one-character pseudonym is not a pseudonym."""
    with pytest.raises(ValueError, match="positive integer"):
        pseudonymized_sample_name("anon_", "s1", length=True)  # type: ignore[arg-type]


def test_an_algorithm_that_needs_a_digest_length_is_refused_by_name() -> None:
    """`shake_128` is in ``algorithms_available`` but its ``hexdigest()`` takes an argument.

    Left unguarded that arrives as a ``TypeError`` from inside ``hashlib``, which is not
    the repository's error convention and does not name the configuration key at fault.
    """
    assert "shake_128" in hashlib.algorithms_available
    with pytest.raises(ValueError, match="shake_128"):
        pseudonymized_sample_name("anon_", "s1", algorithm="shake_128")


def test_an_algorithm_the_backend_refuses_is_reported_by_name(monkeypatch, caplog) -> None:
    """``hashlib.algorithms_available`` lists what is *known*, not what will be computed.

    Under a FIPS-enforcing OpenSSL the provider refuses a non-approved digest at
    construction time and ``hashlib.new`` raises ``ValueError`` -- ``md5`` is listed there
    and still unusable. Only ``TypeError`` was translated, so that ``ValueError`` escaped
    with the backend's own wording and named neither the configured algorithm nor the
    configuration key it came from. The backend is stubbed here because a FIPS provider
    cannot be assumed on a developer machine, and the assertion is about the translation
    rather than about OpenSSL.
    """
    from vntyper.scripts import cohort_pseudonyms

    def _refused_by_the_provider(name: str, data: bytes = b"", **kwargs: object) -> object:
        raise ValueError("[digital envelope routines] unsupported")

    monkeypatch.setattr(cohort_pseudonyms.hashlib, "new", _refused_by_the_provider)

    with (
        caplog.at_level(logging.ERROR, logger="vntyper.scripts.cohort_pseudonyms"),
        pytest.raises(ValueError, match="sha256"),
    ):
        pseudonymized_sample_name("anon_", "s1", algorithm="sha256")

    assert any("sha256" in record.getMessage() for record in caplog.records), (
        "the refusal must be logged at error naming the algorithm, per AGENTS.md"
    )


def test_the_shipped_config_declares_the_pseudonym_settings() -> None:
    """Config-driven, never hardcoded: the values live in config.json."""
    config = json.loads((REPO_ROOT / "vntyper" / "config.json").read_text(encoding="utf-8"))

    assert config["cohort"]["pseudonym"] == {"algorithm": "sha256", "digest_characters": 12}


def test_the_cohort_reads_the_digest_settings_out_of_the_configuration(tmp_path) -> None:
    """`--config-path` is how an operator changes the width, so the read has to work.

    The narrowed digest is asserted through the pseudonymization table rather than
    through the function, because the point is that ``aggregate_cohort`` threads the
    configured pair down to it.
    """
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    config = load_config(None)
    config["cohort"]["pseudonym"] = {"algorithm": "sha1", "digest_characters": 4}

    cohort_summary.aggregate_cohort(
        input_paths=[str(sample_on_disk(tmp_path / "cohort" / "sample_one"))],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config=config,
        pseudonymize_samples="anon_",
    )

    table = (output_dir / "pseudonymization_table.tsv").read_text(encoding="utf-8")
    expected = "anon_" + hashlib.sha1(b"sample_one").hexdigest()[:4]  # noqa: S324 - not a security use
    assert f"{expected}\tsample_one" in table


def test_a_configuration_without_the_pseudonym_block_falls_back_to_the_defaults(tmp_path) -> None:
    """AGENTS.md trap 2: ``--config-path`` replaces the whole config rather than merging it.

    A config that never heard of ``cohort.pseudonym`` - which is every config written
    before this milestone - must produce the default pseudonym rather than a ``KeyError``.
    """
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    cohort_summary.aggregate_cohort(
        input_paths=[str(sample_on_disk(tmp_path / "cohort" / "sample_one"))],
        output_dir=str(output_dir),
        summary_file="cohort_summary.html",
        config={"paths": {"template_dir": "vntyper/templates"}},
        pseudonymize_samples="anon_",
    )

    table = (output_dir / "pseudonymization_table.tsv").read_text(encoding="utf-8")
    assert "anon_c788e939395d\tsample_one" in table


@pytest.mark.parametrize(
    ("cohort_block", "warned_key"),
    [
        (None, "cohort"),
        ({"pseudonym": None}, "cohort.pseudonym"),
        ("sha256", "cohort"),
        ({"pseudonym": ["sha256", 12]}, "cohort.pseudonym"),
        ({"pseudonym": 12}, "cohort.pseudonym"),
    ],
    ids=["cohort-is-null", "pseudonym-is-null", "cohort-is-a-string", "pseudonym-is-a-list", "pseudonym-is-a-number"],
)
def test_a_pseudonym_block_that_is_not_a_mapping_falls_back_to_the_defaults(
    tmp_path, caplog, cohort_block: object, warned_key: str
) -> None:
    """A hand-written config may carry ``"cohort": null``, and JSON has no schema.

    ``config.get("cohort", {}).get("pseudonym", {})`` reads ``.get`` off whatever the two
    keys hold. ``.get("cohort", {})`` only defends against the key being *absent*: present
    and null it returns ``None``, and ``None.get`` is an ``AttributeError`` with no log
    line, no mention of the key at fault, and - because it was raised before the cleanup
    ``try`` - a leaked extraction directory per zip input. It fired whether or not
    pseudonymisation had been asked for, so a config that never intended to use the
    feature still aborted the cohort.

    Every non-mapping shape a JSON document can put in either position is treated the same
    way: fall back to the module defaults and say so at warning, naming the key. A silent
    fallback is not acceptable here - the digest settings decide every pseudonym in the
    report - but neither is a crash on a setting the run was not going to use.
    """
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    config = {"paths": {"template_dir": "vntyper/templates"}, "cohort": cohort_block}

    with caplog.at_level(logging.WARNING, logger="vntyper.scripts.cohort_pseudonyms"):
        cohort_summary.aggregate_cohort(
            input_paths=[str(sample_on_disk(tmp_path / "cohort" / "sample_one"))],
            output_dir=str(output_dir),
            summary_file="cohort_summary.html",
            config=config,
            pseudonymize_samples="anon_",
        )

    table = (output_dir / "pseudonymization_table.tsv").read_text(encoding="utf-8")
    assert "anon_c788e939395d\tsample_one" in table
    assert any(warned_key in record.getMessage() for record in caplog.records), (
        f"the fallback must be logged at warning naming {warned_key!r}, not applied silently"
    )
