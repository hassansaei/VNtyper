"""Pure contracts for fail-closed GHCR manifest-absence classification."""

import importlib
from collections.abc import Callable
from pathlib import Path
from typing import Literal

import pytest

pytestmark = pytest.mark.unit

ROOT = Path(__file__).resolve().parents[2]
REGISTRY_SCRIPT = ROOT / "scripts" / "release_registry.py"


def _classifier() -> Callable[[str, str], Literal["absent", "ambiguous"]]:
    assert REGISTRY_SCRIPT.is_file(), "the focused release registry classifier has not been implemented"
    module = importlib.import_module("scripts.release_registry")
    return module.classify_manifest_absence


@pytest.mark.parametrize(
    ("error_text", "expected"),
    (
        (
            "unexpected status from HEAD request to "
            "https://ghcr.io/v2/hassansaei/vntyper/manifests/latest: 404 Not Found",
            "absent",
        ),
        ("ghcr.io/hassansaei/vntyper:latest: manifest unknown", "absent"),
        ("ghcr.io/hassansaei/vntyper:latest: NOT-FOUND", "absent"),
        ("manifest unknown: not found", "ambiguous"),
        (
            "unexpected status from HEAD request to "
            "https://ghcr.io/v2/hassansaei/vntyper/manifests/other: 404 Not Found",
            "ambiguous",
        ),
        (
            "unexpected status from HEAD request to "
            "https://ghcr.io/v2/hassansaei/vntyper/manifests/latest-evil: 404 Not Found",
            "ambiguous",
        ),
        (
            "HEAD https://ghcr.io/v2/hassansaei/vntyper/manifests/latest failed\ntoken endpoint: 404 Not Found",
            "ambiguous",
        ),
        (
            "ghcr.io/hassansaei/vntyper:latest:\nmanifest unknown",
            "ambiguous",
        ),
        (
            'ghcr.io/hassansaei/vntyper:latest: error getting credentials: executable file not found in $PATH, out: ""',
            "ambiguous",
        ),
        (
            "ghcr.io/hassansaei/vntyper:latest: 401 Unauthorized: manifest unknown",
            "ambiguous",
        ),
    ),
)
def test_manifest_absence_requires_one_exact_reference_bound_to_its_status(
    error_text: str,
    expected: Literal["absent", "ambiguous"],
) -> None:
    """Removing exact-reference or same-line status binding changes a literal verdict."""
    classify_manifest_absence = _classifier()

    assert classify_manifest_absence("ghcr.io/hassansaei/vntyper:latest", error_text) == expected


@pytest.mark.parametrize(
    "reference",
    (
        "docker.io/hassansaei/vntyper:latest",
        "https://ghcr.io/hassansaei/vntyper:latest",
        "ghcr.io/HassanSaei/vntyper:latest",
        "ghcr.io/hassansaei/vntyper",
        "ghcr.io/hassansaei/vntyper:latest/extra",
        "ghcr.io/hassansaei/vntyper:latest\n404 Not Found",
        "ghcr.io/hassansaei/vntyper:latest\rmanifest unknown",
    ),
)
def test_manifest_absence_rejects_malformed_or_multiline_requested_references(reference: str) -> None:
    classify_manifest_absence = _classifier()

    with pytest.raises(ValueError, match="exact tagged GHCR reference"):
        classify_manifest_absence(reference, f"{reference}: 404 Not Found")


def test_manifest_absence_treats_malformed_error_records_as_ambiguous() -> None:
    classify_manifest_absence = _classifier()

    assert (
        classify_manifest_absence(
            "ghcr.io/hassansaei/vntyper:latest",
            "ghcr.io/hassansaei/vntyper:latest:\x00 404 Not Found",
        )
        == "ambiguous"
    )


@pytest.mark.parametrize(
    ("reference", "error_bytes"),
    (
        ("ghcr.io/hassansaei/vntyper:latest\nspoof", b"404 Not Found"),
        ("ghcr.io/hassansaei/vntyper:latest", b"\xff"),
    ),
)
def test_classifier_cli_rejects_malformed_inputs(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    reference: str,
    error_bytes: bytes,
) -> None:
    module = importlib.import_module("scripts.release_registry")
    error_path = tmp_path / "inspect.stderr"
    error_path.write_bytes(error_bytes)

    assert module.main(["classify-absence", reference, str(error_path)]) == 2
    assert "release registry classification failed" in capsys.readouterr().err


def test_classifier_cli_rejects_an_unreadable_error_path(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    module = importlib.import_module("scripts.release_registry")

    assert module.main(["classify-absence", "ghcr.io/hassansaei/vntyper:latest", str(tmp_path / "missing")]) == 2
    assert "release registry classification failed" in capsys.readouterr().err
