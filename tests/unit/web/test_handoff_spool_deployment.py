"""Protected ordinary-job spool and coordinated deployment contract."""

import posixpath
from copy import deepcopy
from pathlib import Path

import pytest
import yaml
from app.config import Settings

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[3]


def _compose() -> dict:
    return yaml.safe_load((REPO_ROOT / "docker" / "docker-compose.yml").read_text(encoding="utf-8"))


def _assert_private_job_mounts(compose: dict) -> None:
    """Normalize Compose syntax and enforce exact protected target ownership."""
    expected_by_target = {
        "/opt/vntyper/handoff": "handoff_spool",
        "/opt/vntyper/output": "result_store",
    }
    assert {"handoff_spool", "result_store"} <= compose["volumes"].keys()
    for service_name, service in compose["services"].items():
        normalized: list[tuple[object, object]] = []
        for mount in service.get("volumes", []):
            if isinstance(mount, str):
                parts = mount.split(":")
                source, target = (None, parts[0]) if len(parts) == 1 else (parts[0], parts[1])
            else:
                source, target = mount.get("source"), mount.get("target")
            if isinstance(target, str):
                target = "/" + posixpath.normpath("/" + target.lstrip("/")).lstrip("/")
            if target in expected_by_target:
                normalized.append((source, target))
        if service_name in {"api", "worker"}:
            assert len(normalized) == len(expected_by_target)
            assert set(normalized) == {(source, target) for target, source in expected_by_target.items()}
        else:
            assert normalized == []


def _assert_result_store_migration_contract(text: str) -> None:
    normalized = " ".join(text.split()).lower()
    for required in (
        "Pause new submissions",
        "drain both the regular and long queues",
        "all active jobs",
        "verify both queues and the active job count are zero",
        "Never purge queued messages",
        "stop the API, workers, and beat",
        "provision the named `result_store`",
        "services remain stopped and detached from the legacy host bind",
        "copy existing unexpired output into `result_store`",
        "retain a backup",
        "cannot be retroactively integrity-attested",
        "accept retirement and unavailability",
        "archive or remove the legacy store",
        "deploy the API, all workers, and beat",
        "verify retained result access",
        "resume submissions",
        "arbitrary same-UID code in either the API or worker service namespace is out of scope",
    ):
        assert required.lower() in normalized


def test_handoff_spool_has_a_dedicated_setting() -> None:
    """API and worker agree on one non-public cross-Celery storage root."""
    assert Settings.DEFAULT_HANDOFF_SPOOL_DIR == "/opt/vntyper/handoff"


def test_only_api_and_worker_mount_the_protected_handoff_spool() -> None:
    """Beat and shared-mount actors receive no spool capability."""
    _assert_private_job_mounts(_compose())


def test_only_api_and_worker_mount_the_service_private_result_store() -> None:
    """Shipped Compose has no host/shared output mount to race before binding."""
    _assert_private_job_mounts(_compose())


@pytest.mark.parametrize(
    ("service", "replacement"),
    [
        ("beat", "/attacker:/opt/vntyper/handoff"),
        ("worker", {"type": "bind", "source": "/attacker", "target": "/opt/vntyper/output"}),
        ("beat", "/opt/vntyper/handoff"),
        ("api", "/opt/vntyper/output"),
        ("beat", "handoff_spool:/opt/vntyper/handoff/"),
        ("beat", "result_store:/opt/vntyper/./output"),
        (
            "beat",
            {"type": "volume", "source": "handoff_spool", "target": "//opt//vntyper//handoff/"},
        ),
    ],
)
def test_private_mount_contract_rejects_any_alternative_source_for_a_protected_target(
    service: str,
    replacement: str | dict[str, str],
) -> None:
    """Exact-string checks cannot hide a second short- or long-form protected mount."""
    compose = deepcopy(_compose())
    compose["services"][service].setdefault("volumes", []).append(replacement)

    with pytest.raises(AssertionError):
        _assert_private_job_mounts(compose)


def test_private_mount_contract_accepts_equivalent_long_named_volume_form() -> None:
    """Security checks normalize Compose short and long mount syntax."""
    compose = deepcopy(_compose())
    replacements = {
        "handoff_spool:/opt/vntyper/handoff": {
            "type": "volume",
            "source": "handoff_spool",
            "target": "/opt/vntyper/handoff",
        },
        "result_store:/opt/vntyper/output": {
            "type": "volume",
            "source": "result_store",
            "target": "/opt/vntyper/output",
        },
    }
    for service_name in ("api", "worker"):
        compose["services"][service_name]["volumes"] = [
            replacements.get(mount, mount) for mount in compose["services"][service_name]["volumes"]
        ]

    _assert_private_job_mounts(compose)


def test_image_creates_the_default_spool_with_service_ownership() -> None:
    """The named volume inherits a path writable by the non-root API and worker."""
    dockerfile = (REPO_ROOT / "docker" / "Dockerfile.base").read_text(encoding="utf-8")

    assert "DEFAULT_HANDOFF_SPOOL_DIR=/opt/vntyper/handoff" in dockerfile
    assert "$DEFAULT_HANDOFF_SPOOL_DIR" in dockerfile


def test_changelog_requires_a_drained_coordinated_deployment() -> None:
    """The incompatible handoff protocol is never described as rolling-safe."""
    changelog = (REPO_ROOT / "docs" / "about" / "changelog.md").read_text(encoding="utf-8")
    unreleased = changelog.split("## Unreleased", 1)[1].split("\n## ", 1)[0]
    release_206 = changelog.split("## 2.0.6", 1)[1].split("\n## ", 1)[0]
    normalized = " ".join(unreleased.split())
    normalized_206 = " ".join(release_206.split())

    assert "Pause new submissions" in unreleased
    assert "resume submissions" in unreleased
    assert "shared legacy input mount access" in normalized
    assert "service-private result store" in normalized
    assert "operator override" in normalized
    assert "weakens this security boundary" in normalized
    assert "Deploy the worker before the API" in normalized_206
    assert "old-API-to-new-worker is safe" in normalized_206
    assert "new-API-to-old-worker is not" in normalized_206
    assert "Pause new submissions" not in release_206
    _assert_result_store_migration_contract(unreleased)


def test_docker_operator_guide_documents_the_private_result_store_boundary() -> None:
    """Operators are warned that replacing the named store changes the threat model."""
    guide = (REPO_ROOT / "docs" / "user-guide" / "docker.md").read_text(encoding="utf-8")
    api_section = guide.split("## API Server", 1)[1]
    normalized = " ".join(api_section.split())

    assert "service-private `result_store`" in normalized
    assert "mounted only into the API and worker" in normalized
    assert "bind-mount override" in normalized
    assert "weakens the shipped security boundary" in normalized
    _assert_result_store_migration_contract(api_section)


def test_docker_readme_api_example_preserves_the_private_result_store_boundary() -> None:
    """The second operator-facing API example does not silently recommend a host result bind."""
    readme = (REPO_ROOT / "docker" / "README.md").read_text(encoding="utf-8")
    api_section = readme.split("### **API Usage**", 1)[1]
    normalized = " ".join(api_section.split())

    assert "vntyper_result_store:/opt/vntyper/output" in api_section
    assert "vntyper_handoff_spool:/opt/vntyper/handoff" in api_section
    assert "bind-mount override" in normalized
    assert "weakens the shipped security boundary" in normalized
    _assert_result_store_migration_contract(api_section)
