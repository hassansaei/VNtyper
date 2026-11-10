"""`/job-summary/{job_id}/` serves the archive's pipeline_summary.json, or 404 like `/download/`."""

import json
import uuid
import zipfile
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

JOB_ID = str(uuid.uuid4())
SUMMARY = {"coverage_qc": {"mean_coverage": 312.0}, "confidence_grade": "high", "version": "2.0.47"}


def _archive(web_app, members: dict[str, bytes]) -> Path:
    path = Path(web_app.DEFAULT_OUTPUT_DIR) / f"{JOB_ID}.zip"
    with zipfile.ZipFile(path, "w") as archive:
        for name, data in members.items():
            archive.writestr(name, data)
    return path


def test_returns_the_summary_from_the_result_archive(client, web_app) -> None:
    _archive(web_app, {f"{JOB_ID}/pipeline_summary.json": json.dumps(SUMMARY).encode(), "x.log": b"log"})

    response = client.get(f"/job-summary/{JOB_ID}/")

    assert response.status_code == 200
    assert response.json() == SUMMARY


def test_prefers_the_run_summary_over_nested_or_lookalike_members(client, web_app) -> None:
    _archive(
        web_app,
        {
            f"{JOB_ID}/advntr/pipeline_summary.json": b'{"nested": true}',
            f"{JOB_ID}/old_pipeline_summary.json": b'{"lookalike": true}',
            f"{JOB_ID}/pipeline_summary.json": json.dumps(SUMMARY).encode(),
        },
    )

    assert client.get(f"/job-summary/{JOB_ID}/").json() == SUMMARY


@pytest.mark.parametrize(
    "members",
    [None, {"x.log": b"no summary"}, {"pipeline_summary.json": b"{not json"}],
    ids=["expired-or-missing", "no-summary-member", "corrupt-summary"],
)
def test_missing_or_unusable_summary_is_404(client, web_app, members) -> None:
    if members is not None:
        _archive(web_app, members)

    response = client.get(f"/job-summary/{JOB_ID}/")

    assert response.status_code == 404
    assert response.json()["detail"] == "File not found"


def test_refuses_an_identifier_the_service_did_not_issue(client) -> None:
    assert client.get("/job-summary/..%2Fetc/").status_code == 404


def test_refuses_a_symlinked_archive(client, web_app, tmp_path: Path) -> None:
    target = tmp_path / "elsewhere.zip"
    with zipfile.ZipFile(target, "w") as archive:
        archive.writestr("pipeline_summary.json", json.dumps(SUMMARY))
    (Path(web_app.DEFAULT_OUTPUT_DIR) / f"{JOB_ID}.zip").symlink_to(target)

    assert client.get(f"/job-summary/{JOB_ID}/").status_code == 404
