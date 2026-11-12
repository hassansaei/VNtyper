"""Anonymous research donations (schema vntyper-donation/2): model, endpoint, storage, logs."""

import json
import logging
import sqlite3

import pytest
from app.config import settings
from app.donations import DonationRepository
from app.donations.repository import K_DEFAULT, K_SENSITIVE, suppress
from fastapi import HTTPException
from starlette.testclient import TestClient

pytestmark = pytest.mark.unit

VALID = {
    "schema": "vntyper-donation/2",
    "vntyper_version": "2.0.49",
    "coverage_qc": "PASS",
    "flank_depth_bin": "100-150",
    "array_depth_bin": "40-60",
    "kestrel_call": True,
    "kestrel_flagged": False,
    "kestrel_confidence": "High_Precision",
    "depth_score_bin": "0.006-0.008",
    "variant_class": "dupC",
    "kit": "Twist Exome 2.0",
    "platform_class": "illumina_short",
    "truth": "confirmed_positive",
    "confirmation_class": "orthogonal_short",
    "referral": "suspected_adtkd",
}
NEGATIVE = {
    **VALID,
    "kestrel_call": False,
    "kestrel_confidence": "none",
    "depth_score_bin": None,
    "variant_class": "none",
    "truth": "unconfirmed",
    "confirmation_class": "none",
}
FORBIDDEN_COLUMNS = {"created_at", "run_id", "sample_hash", "ip", "user_agent", "job_id", "archive"}


@pytest.fixture
def repo(monkeypatch: pytest.MonkeyPatch, tmp_path) -> DonationRepository:
    monkeypatch.setattr(settings, "ENABLE_DONATIONS", True)
    repo = DonationRepository(db_url=str(tmp_path / "donations.db"))
    monkeypatch.setattr("app.main.get_donation_repo", lambda: repo)
    return repo


def _post(client: TestClient, payload: object, content_type: str = "application/json"):
    return client.post("/donations/", content=json.dumps(payload), headers={"content-type": content_type})


def test_status_reports_flag_and_schema(client: TestClient) -> None:
    assert client.get("/donations/status/").json() == {
        "enabled": settings.ENABLE_DONATIONS,
        "schema": "vntyper-donation/2",
    }


def test_valid_donations_are_recorded_without_any_identifier(client: TestClient, repo) -> None:
    for payload in (VALID, NEGATIVE, {k: v for k, v in NEGATIVE.items() if k != "referral"}):
        res = _post(client, payload)
        assert res.status_code == 201, res.text
        assert res.json() == {"status": "recorded"}


@pytest.mark.parametrize(
    "change",
    [
        {"schema": "vntyper-donation/1"},
        {"vntyper_version": "2.0"},
        {"coverage_qc": "NOT_EVALUATED"},
        {"flank_depth_bin": "123"},
        {"array_depth_bin": 80},
        {"kestrel_call": "true"},  # strict: no string coercion
        {"kestrel_confidence": "High_Precision*"},
        {"depth_score_bin": "0.0047"},
        {"variant_class": "c.428dupC"},
        {"kit": "My lab kit"},
        {"platform_class": "Illumina NovaSeq 6000"},
        {"truth": "positive"},
        {"confirmation_class": "Sanger sequencing"},
        {"referral": "ADTKD family"},
        {"sample_name": "DNA23-01234"},  # extra field
        {"run_id": "0b6c"},  # extra field
        {"depth_score_bin": None},  # a call needs its score bin
        {"variant_class": "none"},  # a call needs a class
        {"confirmation_class": "none"},  # confirmed needs a method
        {"kestrel_call": False, "kestrel_confidence": "none"},  # variant_class must then be none
    ],
)
def test_invalid_payloads_are_422_without_echoing_values(client: TestClient, repo, change: dict) -> None:
    res = _post(client, {**VALID, **change})
    assert res.status_code == 422, res.text
    assert "DNA23-01234" not in res.text


@pytest.mark.parametrize(
    "change",
    [
        {"depth_score_bin": "<0.002"},  # neither call nor flag: no score
        {"truth": "confirmed_negative"},  # confirmed needs a method
    ],
)
def test_negative_cross_field_rules(client: TestClient, repo, change: dict) -> None:
    assert _post(client, {**NEGATIVE, **change}).status_code == 422


def test_missing_field_is_422(client: TestClient, repo) -> None:
    assert _post(client, {k: v for k, v in VALID.items() if k != "kit"}).status_code == 422


def test_multipart_and_non_json_are_415(client: TestClient, repo) -> None:
    res = client.post("/donations/", files={"archive": ("r.zip", b"PK", "application/zip")}, data={"kit": "Other"})
    assert res.status_code == 415
    assert _post(client, VALID, content_type="text/plain").status_code == 415


def test_oversize_body_is_413(client: TestClient, repo) -> None:
    assert _post(client, {**VALID, "padding": "x" * 3000}).status_code == 413


def test_disabled_is_403(client: TestClient, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(settings, "ENABLE_DONATIONS", False)
    assert _post(client, VALID).status_code == 403
    assert client.get("/donations/aggregates/").status_code == 403


def test_repository_without_database_url_fails_loudly(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(settings, "POSTGRES_URL", None)
    with pytest.raises(RuntimeError, match="POSTGRES_URL"):
        DonationRepository()


def test_table_holds_only_whitelist_quarter_and_id(client: TestClient, repo, tmp_path) -> None:
    _post(client, VALID)
    conn = sqlite3.connect(tmp_path / "donations.db")
    columns = {row[1] for row in conn.execute("PRAGMA table_info(donations_v2)")}
    tables = {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    (quarter,) = conn.execute("SELECT received_quarter FROM donations_v2").fetchone()
    conn.close()

    assert columns == (set(VALID) - {"schema"}) | {"id", "received_quarter"}
    assert not columns & FORBIDDEN_COLUMNS
    assert tables == {"donations_v2"}
    assert len(quarter) == 7 and quarter[4:6] == "-Q"


def test_suppress_hides_small_cells_and_their_complement() -> None:
    assert suppress({"a": 20, "b": 7, "c": 3}, {}) == {"a": 20, "b": None, "c": None}
    # The sensitive cell needs k=10, and the hidden pair must sum to it.
    assert suppress({"pos": 9, "neg": 30, "unc": 6}, {"pos": K_SENSITIVE}) == {"pos": None, "neg": 30, "unc": None}
    assert suppress({"a": 1, "b": 6, "c": 50}, {}) == {"a": None, "b": None, "c": 50}
    assert suppress({"a": 1, "b": 2, "c": 50, "d": 60}, {}) == {"a": None, "b": None, "c": None, "d": 60}
    assert suppress({"a": 10, "b": 12}, {}) == {"a": 10, "b": 12}
    assert suppress({"a": 0, "b": 15}, {}) == {"a": 0, "b": 15}


@pytest.mark.parametrize("cells", [[0, 4, 30], [1, 1, 1, 40], [9, 30, 31], [5, 5, 0], [2, 100]])
def test_no_hidden_cell_can_be_recovered_from_the_total(cells: list[int]) -> None:
    counts = {str(i): n for i, n in enumerate(cells)}
    out = suppress(counts, {"0": K_SENSITIVE})
    hidden = [c for c, v in out.items() if v is None]
    assert len(hidden) != 1
    if hidden:
        assert sum(counts[c] for c in hidden) >= max(K_SENSITIVE if c == "0" else K_DEFAULT for c in hidden)
    for c, v in out.items():
        assert v is None or v == 0 or v >= (K_SENSITIVE if c == "0" else K_DEFAULT)


def _in_closed_quarter(monkeypatch: pytest.MonkeyPatch, client: TestClient, payload: dict, n: int) -> None:
    from app.donations import repository

    current = repository.received_quarter()
    monkeypatch.setattr(repository, "received_quarter", lambda: "2000-Q1")
    for _ in range(n):
        _post(client, payload)
    monkeypatch.setattr(repository, "received_quarter", lambda: current)


def test_aggregates_are_suppressed_counts_only(client: TestClient, repo, monkeypatch) -> None:
    _in_closed_quarter(monkeypatch, client, VALID, 3)
    small = client.get("/donations/aggregates/").json()
    assert small["total"] is None
    assert all(v is None for name in ("by_truth", "by_kestrel_call", "by_kit") for v in small[name].values())

    _in_closed_quarter(monkeypatch, client, NEGATIVE, 12)
    agg = client.get("/donations/aggregates/").json()
    assert agg["total"] == 15
    # 3 positives sit under k=10, so the 12 unconfirmed hide with them; zeros stay visible.
    assert agg["by_truth"] == {"confirmed_positive": None, "confirmed_negative": 0, "unconfirmed": None}
    assert agg["by_kestrel_call"] == {"true": None, "false": None}
    assert agg["by_kit"]["Twist Exome 2.0"] == 15
    assert not any(isinstance(v, list) for v in agg.values())


def test_aggregates_ignore_the_open_quarter(client: TestClient, repo, monkeypatch) -> None:
    """Donations of the running quarter cannot be isolated by polling before and after."""
    _in_closed_quarter(monkeypatch, client, NEGATIVE, 20)
    before = client.get("/donations/aggregates/").json()
    for _ in range(20):
        _post(client, VALID)
    assert client.get("/donations/aggregates/").json() == before
    assert before["total"] == 20


def _v1_db(tmp_path, rows: int) -> str:
    path = str(tmp_path / "v1.db")
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE donations (id TEXT, run_id TEXT)")
    conn.executemany("INSERT INTO donations VALUES (?, ?)", [(str(i), "r") for i in range(rows)])
    conn.commit()
    conn.close()
    return path


def _tables(path: str) -> set[str]:
    conn = sqlite3.connect(path)
    names = {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    conn.close()
    return names


def test_empty_v1_table_is_dropped(tmp_path) -> None:
    path = _v1_db(tmp_path, 0)
    DonationRepository(db_url=path)
    assert _tables(path) == {"donations_v2"}


def test_v1_table_with_rows_is_kept_with_an_error(tmp_path, caplog) -> None:
    path = _v1_db(tmp_path, 2)
    DonationRepository(db_url=path)
    assert _tables(path) == {"donations", "donations_v2"}
    assert [r.levelname for r in caplog.records if "Legacy table" in r.getMessage()] == ["ERROR"]


def test_donation_requests_log_no_address_agent_or_payload(client: TestClient, repo, web_app, caplog) -> None:
    caplog.set_level(logging.DEBUG)
    sentinel = {**VALID, "kit": "WGS (PCR-free)"}
    _post(client, sentinel)
    _post(client, {**sentinel, "truth": "NOT-A-VALUE"})

    async def _limited() -> None:
        raise HTTPException(status_code=429)

    web_app.app.dependency_overrides[web_app.high_rate_limiter] = _limited
    assert _post(client, sentinel).status_code == 429
    logging.getLogger("uvicorn.access").info(
        '%s - "%s %s HTTP/%s" %d', "203.0.113.7:5000", "POST", "/api/donations/", "1.1", 201
    )

    for needle in ("testclient", "203.0.113.7", "WGS (PCR-free)", "NOT-A-VALUE", "python-httpx"):
        assert needle not in caplog.text
    assert "the donation endpoint" in caplog.text


def test_other_paths_keep_their_access_log(caplog) -> None:
    caplog.set_level(logging.INFO)
    import app.main  # noqa: F401  (installs the filter)

    logging.getLogger("uvicorn.access").info(
        '%s - "%s %s HTTP/%s" %d', "203.0.113.7:5000", "GET", "/api/health/", "1.1", 200
    )
    assert "/api/health/" in caplog.text
