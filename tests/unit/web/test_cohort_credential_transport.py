"""How a cohort passphrase reaches the two routes that read a cohort.

`/cohort-status/` and `/cohort-download/` are GETs, and they took the cohort's
passphrase as a query parameter. A query string is part of the request line: it
is written to the access log of every server and proxy the request passes
through, it is kept in browser history, and it is sent in `Referer` headers. The
passphrase is the only credential this service has, so putting it there means
the credential is recorded in places that are not credential stores and are not
managed as such.

The `X-Cohort-Passphrase` request header carries it instead. Headers are not
logged by default by any of those, and the change is additive: the query
parameter keeps working, because the existing web UI and the `vntyper online`
CLI subcommand both use it and breaking them to fix where a value is written
would be a worse trade. It is marked deprecated in the API description so a
client author is told which one to reach for.

Precedence is stated rather than left to fall out of the implementation: the
header wins. "Whichever verifies" would turn two credentials into two attempts
at one, and a caller sending both should be told plainly which one was used.

`docker` is put on `sys.path` by `tests/unit/web/conftest.py`, which pytest
imports before this module.
"""

from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

from app.cohorts import COHORT_KEY_PREFIX, preferred_passphrase  # noqa: E402

PASSPHRASE = "correct-horse-battery-staple"
WRONG_PASSPHRASE = "not-the-passphrase"
HEADER = "X-Cohort-Passphrase"

COHORT_READ_ROUTES = ["/cohort-status/", "/cohort-download/"]


def _cohort_with_one_result(client, fake_redis, tmp_path: Path) -> str:
    """Create a protected cohort holding one job that has a result archive.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
        tmp_path: The directory the `web_app` fixture configured as the job tree.

    Returns:
        str: The identifier of the created cohort.
    """
    created = client.post("/create-cohort/", data={"alias": "study", "passphrase": PASSPHRASE})
    assert created.status_code == 200, created.text
    cohort_id = created.json()["cohort_id"]

    member = "5d0b8e21-4c37-4f96-b3a5-9e1c7f804d62"
    fake_redis.sadd(f"{COHORT_KEY_PREFIX}{cohort_id}:jobs", member)
    (tmp_path / "output" / f"{member}.zip").write_bytes(b"member-result")
    return cohort_id


# ---------------------------------------------------------------------------
# The choice itself, with no app in the way
# ---------------------------------------------------------------------------


def test_the_header_is_used_when_it_carries_a_value() -> None:
    """A header with content is the credential, whatever the query said."""
    assert preferred_passphrase("from-header", "from-query") == "from-header"


def test_the_query_parameter_is_used_when_the_header_is_absent() -> None:
    """The older way in still works when it is the only one used."""
    assert preferred_passphrase(None, "from-query") == "from-query"


@pytest.mark.parametrize("header_value", [None, "", "   "])
def test_an_empty_header_does_not_displace_the_query_parameter(header_value: str | None) -> None:
    """A header sent empty means "not supplied", not "supplied as nothing".

    A proxy that adds the header unconditionally would otherwise blank out a
    credential the caller did send.

    Args:
        header_value: The empty header form under test.
    """
    assert preferred_passphrase(header_value, "from-query") == "from-query"


def test_neither_way_in_yields_nothing() -> None:
    """With no credential at all the answer is None, for the caller to refuse."""
    assert preferred_passphrase(None, None) is None


def test_the_chosen_value_is_never_trimmed() -> None:
    """Surrounding whitespace is part of the passphrase, not framing around it.

    The stored hash was made from the passphrase as it was typed, so trimming
    here would lock out anyone whose passphrase begins or ends with a space.
    """
    assert preferred_passphrase(" spaced ", None) == " spaced "


# ---------------------------------------------------------------------------
# Endpoint level: both ways in, on both routes
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("route", COHORT_READ_ROUTES)
def test_the_header_opens_a_cohort(client, fake_redis, tmp_path: Path, route: str) -> None:
    """The passphrase can be supplied without appearing in the request line.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
        route: The cohort read route under test.
    """
    cohort_id = _cohort_with_one_result(client, fake_redis, tmp_path)

    response = client.get(route, params={"cohort_id": cohort_id}, headers={HEADER: PASSPHRASE})

    assert response.status_code == 200, response.text


@pytest.mark.parametrize("route", COHORT_READ_ROUTES)
def test_the_query_parameter_still_opens_a_cohort(client, fake_redis, tmp_path: Path, route: str) -> None:
    """The older way in is unchanged; existing clients keep working.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
        route: The cohort read route under test.
    """
    cohort_id = _cohort_with_one_result(client, fake_redis, tmp_path)

    response = client.get(route, params={"cohort_id": cohort_id, "passphrase": PASSPHRASE})

    assert response.status_code == 200, response.text


@pytest.mark.parametrize("route", COHORT_READ_ROUTES)
def test_a_wrong_header_is_refused_even_beside_a_correct_query_parameter(
    client, fake_redis, tmp_path: Path, route: str
) -> None:
    """The header wins, so a wrong one is a refusal rather than a first attempt.

    This is the assertion that distinguishes "the header takes precedence" from
    "either will do": with a correct query parameter beside it, an implementation
    that tried both would answer 200.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
        route: The cohort read route under test.
    """
    cohort_id = _cohort_with_one_result(client, fake_redis, tmp_path)

    response = client.get(
        route,
        params={"cohort_id": cohort_id, "passphrase": PASSPHRASE},
        headers={HEADER: WRONG_PASSPHRASE},
    )

    assert response.status_code == 401
    assert response.json()["detail"] == "Incorrect passphrase"


@pytest.mark.parametrize("route", COHORT_READ_ROUTES)
def test_a_correct_header_wins_over_a_wrong_query_parameter(client, fake_redis, tmp_path: Path, route: str) -> None:
    """The precedence holds in the other direction too.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
        route: The cohort read route under test.
    """
    cohort_id = _cohort_with_one_result(client, fake_redis, tmp_path)

    response = client.get(
        route,
        params={"cohort_id": cohort_id, "passphrase": WRONG_PASSPHRASE},
        headers={HEADER: PASSPHRASE},
    )

    assert response.status_code == 200, response.text


@pytest.mark.parametrize("route", COHORT_READ_ROUTES)
def test_neither_way_in_is_still_refused(client, fake_redis, tmp_path: Path, route: str) -> None:
    """Adding a second way to supply the credential does not make it optional.

    Args:
        client: TestClient fixture from conftest.
        fake_redis: The store backing the app's cohort client.
        tmp_path: The directory the `web_app` fixture configured as the job tree.
        route: The cohort read route under test.
    """
    cohort_id = _cohort_with_one_result(client, fake_redis, tmp_path)

    response = client.get(route, params={"cohort_id": cohort_id})

    assert response.status_code == 401
    assert response.json()["detail"] == "Passphrase required for this cohort"


# ---------------------------------------------------------------------------
# The API description has to say which one to use
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("route", COHORT_READ_ROUTES)
def test_the_schema_offers_the_header_and_deprecates_the_query_parameter(client, route: str) -> None:
    """A client author reading the schema is told which way in to use.

    Args:
        client: TestClient fixture from conftest.
        route: The cohort read route under test.
    """
    schema = client.get("/openapi.json").json()
    parameters = schema["paths"][route]["get"]["parameters"]
    by_name = {(parameter["in"], parameter["name"]): parameter for parameter in parameters}

    header = by_name[("header", HEADER)]
    assert not header.get("deprecated", False)

    query = by_name[("query", "passphrase")]
    assert query["deprecated"] is True
    assert HEADER in query["description"]
