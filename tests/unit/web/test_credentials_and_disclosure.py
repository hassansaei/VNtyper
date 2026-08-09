"""Credential handling and error-disclosure guards for the web service.

Two independent properties are pinned here.

**No literal credential defaults.** `docker/app` is public source. A module that
falls back to a literal password when its environment variable is missing turns
that literal into a shared secret, and two modules that fall back to *different*
literals authenticate the API and the worker differently against the same Redis
instance. `test_no_module_ships_a_literal_credential_default` walks every module
under `docker/app` with `ast` and fails on either shape.

**No internal detail on the wire.** `/job-status/` is unauthenticated. A failed
Celery task carries the exception it died on, which for the pipeline is a
`CalledProcessError` holding the full argv and container paths. The endpoint must
return a generic message and log the detail instead.

Pure unit tier: `ast` over the source tree, plus the in-process `client` fixture.
No sockets, no Docker, no reference data.
"""

from __future__ import annotations

import ast
import asyncio
import logging
import re
import subprocess
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[3]
APP_DIR = REPO_ROOT / "docker" / "app"
COMPOSE_FILE = REPO_ROOT / "docker" / "docker-compose.yml"

# Celery reads these from the process environment and they take precedence over
# the broker/backend passed to `Celery(...)`. Setting either one outside the
# application silently replaces the URL `celery_app.py` builds - credential,
# percent-encoding and all - with whatever string is in the variable.
CELERY_URL_OVERRIDES = ("CELERY_BROKER_URL", "CELERY_RESULT_BACKEND")

# A name is treated as naming a credential when one of these words is its final
# underscore-separated component. Anchoring at the end keeps `PASSWORD_HASH_SCHEME`
# (a bcrypt scheme name, not a secret) out while keeping `SMTP_PASSWORD` in.
CREDENTIAL_NAME = re.compile(r"(?:^|_)(PASSWORD|PASSWD|SECRET|TOKEN|KEY|CREDENTIAL|CREDENTIALS|PASSPHRASE)$")

# Literals that are allowed to remain as defaults.
#
# The choice made here: an *unmistakable* placeholder passes, a plausible-looking
# literal fails. A value like "your_smtp_password" cannot be mistaken for a working
# credential - nothing authenticates with it, so it fails loudly at first use and
# documents the variable's existence. A value that looks like a real password does
# the opposite: it authenticates, so a deployment that forgot to set the variable
# comes up silently on a secret that is public in the source tree.
#
# The allowance is an explicit enumeration rather than a pattern on purpose. Adding
# an entry is a deliberate, reviewable act, and a pattern such as "starts with
# your_" would be trivially satisfied by naming a real secret `your_...`.
ALLOWED_PLACEHOLDER_DEFAULTS = frozenset({"your_smtp_password"})

# Calls whose second positional argument is an environment-variable fallback.
ENV_LOOKUPS = frozenset({"os.getenv", "os.environ.get", "os.environ.setdefault", "getenv"})

# A job identifier in the form the service issues, which is the only form the
# job routes accept -- see `test_job_identifiers.py`.
FAILED_JOB_ID = "3f2a7c18-5b6d-4e9a-8c31-0d4f7a2b9e65"

# What the pipeline's own failure looks like when it reaches /job-status/: the
# Celery task shells out through conda, so the exception carries the argv and the
# per-job paths inside the container.
FAILED_TASK_INFO = subprocess.CalledProcessError(
    returncode=1,
    cmd=[
        "conda",
        "run",
        "-n",
        "vntyper",
        "vntyper",
        "pipeline",
        "--bam",
        f"/opt/vntyper/input/{FAILED_JOB_ID}/sample.bam",
        "-o",
        f"/opt/vntyper/output/{FAILED_JOB_ID}/",
    ],
)


def _dotted_name(node: ast.AST) -> str | None:
    """Render an attribute/name expression as a dotted string.

    Args:
        node: The expression forming a call's `func`.

    Returns:
        str | None: e.g. `"os.environ.get"`, or None if the expression is not a
            plain chain of names and attributes.
    """
    parts: list[str] = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if not isinstance(node, ast.Name):
        return None
    parts.append(node.id)
    return ".".join(reversed(parts))


def _nonempty_str_constant(node: ast.AST) -> str | None:
    """Return the value of a non-empty string literal, else None.

    Args:
        node: Any expression node.

    Returns:
        str | None: The literal's value, or None if the node is not a non-empty
            string constant.
    """
    if isinstance(node, ast.Constant) and isinstance(node.value, str) and node.value:
        return node.value
    return None


def _module_level_string_constants(tree: ast.Module) -> dict[str, str]:
    """Map module-level names to the string literals they are bound to.

    Used to resolve `os.getenv("X", SOME_NAME)` back to the literal `SOME_NAME`
    holds, which is the shape the defect actually took.

    Args:
        tree: A parsed module.

    Returns:
        dict: Name to literal value, for module-level bindings only.
    """
    constants: dict[str, str] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, ast.AnnAssign):
            targets = [node.target]
        else:
            continue
        value = _nonempty_str_constant(node.value) if node.value is not None else None
        if value is None:
            continue
        for target in targets:
            if isinstance(target, ast.Name):
                constants[target.id] = value
    return constants


def _credential_defaults(path: Path) -> list[tuple[str, str]]:
    """Find literal credential defaults in one module.

    Two shapes are covered:

    1. An assignment (`ast.Assign` / `ast.AnnAssign`) whose target name ends in a
       credential word and whose value is a non-empty string literal.
    2. An environment lookup (`os.getenv` and friends) whose *variable name* ends
       in a credential word and whose fallback is a non-empty string literal, or a
       module-level name bound to one.

    Knowingly not covered, because none of them is the shape in play and each
    would cost more in false positives than it buys: credentials assembled at
    runtime (f-strings, concatenation), defaults in function signatures, values
    inside dict or list literals, secrets loaded from a committed data file, and
    names that do not read as credentials at all.

    Args:
        path: Module to scan.

    Returns:
        list: `(name, literal)` pairs, one per finding.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    module_constants = _module_level_string_constants(tree)
    findings: list[tuple[str, str]] = []
    pattern = CREDENTIAL_NAME

    def _record(name: str, literal: str | None) -> None:
        if literal and literal not in ALLOWED_PLACEHOLDER_DEFAULTS:
            findings.append((name, literal))

    for node in ast.walk(tree):
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for target in targets:
                if isinstance(target, ast.Name) and pattern.search(target.id.upper()):
                    _record(target.id, _nonempty_str_constant(node.value) if node.value else None)
        elif isinstance(node, ast.Call) and len(node.args) >= 2:
            if _dotted_name(node.func) not in ENV_LOOKUPS:
                continue
            var = node.args[0]
            if not (isinstance(var, ast.Constant) and isinstance(var.value, str)):
                continue
            if not pattern.search(var.value.upper()):
                continue
            fallback = node.args[1]
            literal = _nonempty_str_constant(fallback)
            if literal is None and isinstance(fallback, ast.Name):
                literal = module_constants.get(fallback.id)
            _record(var.value, literal)

    return findings


def test_no_module_ships_a_literal_credential_default() -> None:
    """No module under docker/app may fall back to a literal credential.

    Fails naming the file, the symbol and - deliberately - not the value, so the
    failure message never becomes another copy of the secret.
    """
    modules = sorted(APP_DIR.glob("*.py"))
    assert modules, f"no modules found under {APP_DIR}; has the layout changed?"

    offenders = [f"{path.relative_to(REPO_ROOT)}: {name}" for path in modules for name, _ in _credential_defaults(path)]
    assert not offenders, (
        "literal credential default(s) in the source tree: "
        + "; ".join(offenders)
        + ". Read the value from the environment with no fallback and fail fast when it is unset."
    )


def test_the_scan_catches_a_planted_credential_default(tmp_path: Path) -> None:
    """The scan must actually detect both shapes, not vacuously pass."""
    planted = tmp_path / "planted.py"
    planted.write_text(
        "import os\n"
        "FALLBACK = 'not-a-real-secret'\n"
        "DEFAULT_REDIS_PASSWORD = 'not-a-real-secret'\n"
        "REDIS_PASSWORD = os.getenv('REDIS_PASSWORD', FALLBACK)\n"
        "SMTP_PASSWORD = os.environ.get('SMTP_PASSWORD', 'not-a-real-secret')\n",
        encoding="utf-8",
    )
    names = {name for name, _ in _credential_defaults(planted)}
    assert names == {"DEFAULT_REDIS_PASSWORD", "REDIS_PASSWORD", "SMTP_PASSWORD"}


def test_the_scan_ignores_names_that_are_not_credentials(tmp_path: Path) -> None:
    """A non-credential name holding a literal must not trip the scan."""
    planted = tmp_path / "innocent.py"
    planted.write_text(
        "PASSWORD_HASH_SCHEME = 'bcrypt'\nSMTP_USERNAME = 'your_smtp_username'\nREDIS_HOST = 'redis'\n",
        encoding="utf-8",
    )
    assert _credential_defaults(planted) == []


def test_celery_prefers_the_environment_over_the_url_the_application_builds(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Pin the Celery behaviour the compose assertion below depends on.

    `Celery(broker=...)` looks like the last word on the broker URL. It is not:
    `CELERY_BROKER_URL` in the environment wins. This test states that, so the
    next reader does not have to rediscover it.
    """
    from celery import Celery

    monkeypatch.setenv("CELERY_BROKER_URL", "redis://from-environment:6379/9")
    probe = Celery("probe", broker="redis://:secret@from-constructor:6379/0")
    assert probe.conf.broker_url == "redis://from-environment:6379/9"


def test_compose_declares_the_credential_and_does_not_override_the_broker_url() -> None:
    """The shipped compose file must not undo the credential handling.

    Text-based rather than YAML-parsed, matching `test_version_consistency.py`:
    the assertion is about which variable names appear at all, so a parser buys
    nothing and adds a dependency the unit tier does not otherwise need. Whole-line
    comments are stripped first, so the file can explain the rule it obeys.
    """
    lines = COMPOSE_FILE.read_text(encoding="utf-8").splitlines()
    compose = "\n".join(line for line in lines if not line.lstrip().startswith("#"))

    assert "REDIS_PASSWORD" in compose, (
        f"{COMPOSE_FILE.name} must declare REDIS_PASSWORD for Redis and for every service that "
        "talks to it; there is no default any more."
    )

    present = [name for name in CELERY_URL_OVERRIDES if name in compose]
    assert not present, (
        f"{COMPOSE_FILE.name} sets {', '.join(present)}. Celery reads those from the environment "
        "in preference to the broker/backend passed to Celery(...), so they replace the "
        "credential-carrying URL app/celery_app.py builds. Remove them and let the application "
        "build the URL from REDIS_PASSWORD."
    )


def test_redis_password_has_no_default_and_is_read_lazily(monkeypatch: pytest.MonkeyPatch) -> None:
    """`get_redis_password` reads the environment on each call, with no fallback."""
    from app.config import get_redis_password

    monkeypatch.delenv("REDIS_PASSWORD", raising=False)
    assert get_redis_password() is None

    monkeypatch.setenv("REDIS_PASSWORD", "  ")
    assert get_redis_password() == "  "

    monkeypatch.setenv("REDIS_PASSWORD", "set-after-import")
    assert get_redis_password() == "set-after-import"


def test_require_redis_password_fails_fast_when_unset(monkeypatch: pytest.MonkeyPatch) -> None:
    """`require_redis_password` raises RuntimeError naming the variable."""
    from app.config import require_redis_password

    monkeypatch.delenv("REDIS_PASSWORD", raising=False)
    with pytest.raises(RuntimeError, match="REDIS_PASSWORD"):
        require_redis_password()

    monkeypatch.setenv("REDIS_PASSWORD", "")
    with pytest.raises(RuntimeError, match="REDIS_PASSWORD"):
        require_redis_password()

    monkeypatch.setenv("REDIS_PASSWORD", "configured")
    assert require_redis_password() == "configured"


def test_redis_url_percent_encodes_the_credential() -> None:
    """A password containing URL delimiters must survive URL construction.

    `@`, `:`, `/` and `#` all mean something inside a URL's authority. Left raw
    they redirect the client to a different host, port or path.
    """
    from urllib.parse import unquote, urlsplit

    import redis
    from app.config import build_redis_url

    password = "p@ss:w/rd#1"
    url = build_redis_url("redis", 6379, 2, password)

    parts = urlsplit(url)
    assert parts.scheme == "redis"
    assert parts.hostname == "redis"
    assert parts.port == 6379
    assert parts.path == "/2"
    # `SplitResult.password` is `str | None`. Asserting it is present first is the
    # difference between "the URL carries no credential at all" failing as a clear
    # assertion here and failing as a TypeError inside unquote() several frames down.
    assert parts.password is not None, "the URL carries no credential section"
    assert unquote(parts.password) == password
    assert password not in url, "the credential must be encoded, not interpolated raw"

    # The real consumer: redis-py (and, through kombu, Celery) parses the URL back.
    kwargs = redis.Redis.from_url(url).connection_pool.connection_kwargs
    assert kwargs["host"] == "redis"
    assert kwargs["port"] == 6379
    assert kwargs["db"] == 2
    assert kwargs["password"] == password


def test_redis_url_omits_credentials_when_no_password_is_configured() -> None:
    """With no password the URL carries no empty credential section."""
    from app.config import build_redis_url

    assert build_redis_url("redis", 6379, 0, None) == "redis://redis:6379/0"


def test_api_and_worker_read_the_credential_from_the_same_accessor(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`main` and `celery_app` must resolve to one value, not two copies."""
    from app import celery_app as app_celery
    from app import config as app_config
    from app import main as app_main
    from app import tasks as app_tasks

    monkeypatch.setenv("REDIS_PASSWORD", "one-shared-value")
    expected = app_config.get_redis_password()

    for module in (app_main, app_celery, app_tasks):
        assert module.get_redis_password is app_config.get_redis_password
        assert module.get_redis_password() == expected


def test_api_startup_fails_fast_without_a_configured_credential(web_app, monkeypatch: pytest.MonkeyPatch) -> None:
    """The API startup event refuses to come up with REDIS_PASSWORD unset.

    The check must be the first thing startup does: everything after it opens a
    connection or spawns a subprocess, neither of which the unit tier allows.
    """
    monkeypatch.delenv("REDIS_PASSWORD", raising=False)
    with pytest.raises(RuntimeError, match="REDIS_PASSWORD"):
        asyncio.run(web_app.startup_event())


@pytest.mark.parametrize("signal_name", ["celeryd_init", "beat_init"])
def test_worker_startup_fails_fast_without_a_configured_credential(
    signal_name: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Celery worker and beat refuse to come up with REDIS_PASSWORD unset.

    Driven through `Signal.send`, not by calling the receiver directly, because
    that dispatch is the part that can go wrong: Celery catches `Exception` from
    a signal receiver and merely logs it, so a receiver that raises
    `RuntimeError` does not stop the process. Only `BaseException` gets through.
    """
    import app.celery_app  # noqa: F401  (registers the receivers)
    from celery import signals

    signal = getattr(signals, signal_name)
    assert signal.receivers, f"{signal_name} has no receiver; the startup check is not wired up"

    monkeypatch.delenv("REDIS_PASSWORD", raising=False)
    with pytest.raises(SystemExit, match="REDIS_PASSWORD"):
        signal.send(sender=f"test-{signal_name}")

    monkeypatch.setenv("REDIS_PASSWORD", "configured")
    signal.send(sender=f"test-{signal_name}")


def _register_failed_job(fake_redis, web_app, monkeypatch: pytest.MonkeyPatch) -> None:
    """Point a job identifier at a task that Celery reports as FAILURE.

    Args:
        fake_redis: The in-process Redis stand-in the app is patched onto.
        web_app: The imported `app.main` module.
        monkeypatch: Used to replace `AsyncResult` in `app.main`.
    """
    fake_redis.set(FAILED_JOB_ID, "task-1")
    monkeypatch.setattr(
        web_app,
        "AsyncResult",
        lambda task_id: SimpleNamespace(status="FAILURE", info=FAILED_TASK_INFO),
    )


def test_job_status_of_a_failed_job_does_not_disclose_internal_detail(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed job reports that it failed, without the argv or container paths."""
    _register_failed_job(fake_redis, web_app, monkeypatch)

    response = client.get(f"/job-status/{FAILED_JOB_ID}/")

    assert response.status_code == 200
    payload = response.json()
    assert payload["job_id"] == FAILED_JOB_ID
    assert payload["status"] == "failed"

    error = payload["error"]
    assert "/opt/" not in error
    assert "conda" not in error.lower()
    assert "vntyper pipeline" not in error
    assert "returned non-zero exit status" not in error
    assert "sample.bam" not in error


def test_job_status_of_a_failed_job_still_says_something_useful(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The generic message must be a message, not an empty string.

    Asserting only on absence would pass trivially if the endpoint stopped
    returning anything at all.
    """
    _register_failed_job(fake_redis, web_app, monkeypatch)

    error = client.get(f"/job-status/{FAILED_JOB_ID}/").json()["error"]

    assert error, "a failed job must still explain that it failed"
    assert "fail" in error.lower()
    assert len(error.split()) >= 4, f"not an informative message: {error!r}"


def test_job_status_returns_the_stored_curated_preflight_message(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A known preflight failure is more useful than the generic failure text."""
    _register_failed_job(fake_redis, web_app, monkeypatch)
    message = "Unable to resolve CRAM reference: contig=chr1, M5=digest."
    fake_redis.hset(
        f"usage:{FAILED_JOB_ID}",
        mapping={"code": "reference_unresolved", "message": message},
    )

    response = client.get(f"/job-status/{FAILED_JOB_ID}/")

    assert response.json()["error"] == message


def test_remote_header_policy_failure_survives_artifact_redis_status_transport(
    client,
    web_app,
    fake_redis,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The exact curated policy message reaches the API without exposing its remote URI."""
    from app.job_failures import read_preflight_failure

    from vntyper.scripts.pipeline_alignment import prepare_input_alignment_preflight

    _register_failed_job(fake_redis, web_app, monkeypatch)
    output = tmp_path / "pipeline-output"
    output.mkdir()
    remote_uri = "http://127.0.0.1:8765/private/reference.fa"
    header = f"@SQ\tSN:chr7\tLN:100\tUR:{remote_uri}\n"
    expected_message = (
        "Remote CRAM header reference is disabled by policy: contig=chr7, scheme=http. Replace the @SQ UR with "
        "a local path, relative path, or file-scheme reference, or set "
        "cram.allow_ambient_reference_resolution=true to accept "
        "network access."
    )

    with (
        mock.patch("vntyper.scripts.pipeline_alignment.read_alignment_header", return_value=header),
        pytest.raises(ValueError, match="Remote CRAM header reference is disabled") as raised,
    ):
        prepare_input_alignment_preflight(
            in_path=tmp_path / "patient-input" / "sample.cram",
            input_type="CRAM",
            output_dir=output,
            config={"cram": {"local_ref_path": "/local/cache/%s"}},
            threads=1,
            reference_assembly="hg19",
            bed_file=None,
            custom_regions=None,
            reference_fasta=None,
            fast_mode=False,
        )

    assert str(raised.value) == expected_message
    transported = read_preflight_failure(output)
    assert transported == {"code": "reference_policy_invalid", "message": expected_message}
    assert remote_uri not in expected_message
    fake_redis.hset(f"usage:{FAILED_JOB_ID}", mapping=transported)

    response = client.get(f"/job-status/{FAILED_JOB_ID}/")

    assert response.json()["error"] == expected_message


def test_job_status_uses_the_generic_message_when_no_curated_message_is_stored(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Non-preflight failures keep the existing no-disclosure response."""
    _register_failed_job(fake_redis, web_app, monkeypatch)

    response = client.get(f"/job-status/{FAILED_JOB_ID}/")

    assert response.json()["error"] == web_app.JOB_FAILURE_MESSAGE


def test_job_status_refuses_a_stored_message_with_an_absolute_worker_path(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Even corrupt Redis state cannot turn the unauthenticated route into a path oracle."""
    _register_failed_job(fake_redis, web_app, monkeypatch)
    worker_path = f"/opt/vntyper/output/{FAILED_JOB_ID}/private.cram"
    fake_redis.hset(
        f"usage:{FAILED_JOB_ID}",
        mapping={"code": "reference_unresolved", "message": f"Cannot decode {worker_path}"},
    )

    response = client.get(f"/job-status/{FAILED_JOB_ID}/")

    assert response.json()["error"] == web_app.JOB_FAILURE_MESSAGE
    assert worker_path not in response.text


def test_job_status_failure_detail_is_still_logged_server_side(
    client, web_app, fake_redis, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """Withholding the detail from the client must not lose it for the operator."""
    _register_failed_job(fake_redis, web_app, monkeypatch)

    with caplog.at_level(logging.ERROR, logger="app.main"):
        client.get(f"/job-status/{FAILED_JOB_ID}/")

    logged = [record.getMessage() for record in caplog.records if record.levelno >= logging.ERROR]
    assert logged, "the failure detail was neither returned nor logged"
    assert any(FAILED_JOB_ID in message for message in logged)
    assert any(f"/opt/vntyper/input/{FAILED_JOB_ID}/sample.bam" in message for message in logged)
