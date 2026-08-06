"""Reading the peer address off a request that may not have one.

`/run-job/`, `/run-cohort-analysis/` and the rate-limit branch of the HTTP
exception handler all want the caller's address, for a usage-statistics hash and
for a log line. `Request.client` comes from the optional `client` key of the ASGI
scope, so it is `None` whenever the server did not supply one, and reading
`.host` off it unguarded raises `AttributeError` - which, inside an endpoint,
becomes a 500 on a request that was otherwise fine.

`docker` is put on `sys.path` by `tests/unit/web/conftest.py`, which pytest
imports before this module, so `app.utils` is importable here.
"""

import pytest
from starlette.requests import Request

pytestmark = pytest.mark.unit


def _request(scope_extras: dict) -> Request:
    """Build a request from a minimal HTTP scope.

    Args:
        scope_extras: Keys merged into the scope, so a test can supply a
            `client` entry or deliberately leave it out.

    Returns:
        Request: A request backed by that scope.
    """
    scope = {"type": "http", "method": "POST", "path": "/run-job/", "headers": []}
    scope.update(scope_extras)
    return Request(scope)


def test_the_peer_address_is_returned_when_the_scope_carries_one() -> None:
    """The ordinary case: a served request names its peer."""
    from app.utils import client_host

    assert client_host(_request({"client": ("203.0.113.7", 54321)})) == "203.0.113.7"


def test_a_scope_without_a_client_answers_none_rather_than_raising() -> None:
    """`client` is optional in the ASGI scope, so its absence must be an answer.

    The callers pass the result on as an optional value - into the usage hash
    and into a log line - so None is usable where an `AttributeError` would
    have turned the whole request into a 500.
    """
    from app.utils import client_host

    assert client_host(_request({})) is None
