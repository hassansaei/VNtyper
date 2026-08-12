"""The retention window, and the invariant that a cohort never outlives its members.

`/download/{job_id}/` requires no credential, so `MAX_RESULT_AGE_DAYS` *is* the exposure
window for a completed job: the id is the capability. Shortening it is a mitigation for
#189, not a fix -- it reduces exposure duration, not exposure -- so what is pinned here is
the invariant, not the number.

The clamp matters because the two lifetimes are driven differently.
`tasks.delete_old_results` reclaims each archive on a fixed clock from *its own* creation,
while `extend_cohort_retention` fires on every job completion and every cohort analysis, so
an actively-used cohort rolls its TTL forward from its *most recent* job. Comparing the two
constants understates the gap: an active cohort can stay openable with almost none of its
members' archives intact, listing downloads that 404.

``docker`` is put on ``sys.path`` by ``tests/unit/web/conftest.py``, which pytest imports
before this module, so ``app.config`` is importable here.
"""

import logging

import pytest
from app.config import Settings

pytestmark = pytest.mark.unit


def _settings(cohort_days: int, archive_days: int) -> Settings:
    """Build a Settings carrying one pair of retention values.

    Args:
        cohort_days: The configured ``COHORT_RETENTION_DAYS``.
        archive_days: The configured ``MAX_RESULT_AGE_DAYS``.

    Returns:
        Settings: An instance with those two values overridden on the instance.
    """
    settings = Settings()
    settings.COHORT_RETENTION_DAYS = cohort_days
    settings.MAX_RESULT_AGE_DAYS = archive_days
    return settings


def test_a_cohort_never_outlives_the_archives_it_lists():
    """The invariant, stated independently of the shipped numbers.

    Days 3 to 14 were a window in which ``/cohort-status/`` listed member job ids whose
    ``/download/{job_id}/`` returned 404 -- a cohort page advertising downloads that do not
    work.
    """
    assert _settings(cohort_days=14, archive_days=3).cohort_retention_days() == 3


def test_a_cohort_shorter_than_the_archive_window_is_left_alone():
    """The clamp is a ceiling, not an assignment. An operator who deliberately wants a
    shorter-lived cohort than its archives keeps it."""
    assert _settings(cohort_days=2, archive_days=3).cohort_retention_days() == 2


def test_equal_lifetimes_are_not_clamped_and_do_not_warn(caplog):
    """The boundary. Equal is already consistent, so there is nothing to say about it."""
    with caplog.at_level(logging.WARNING):
        assert _settings(cohort_days=3, archive_days=3).cohort_retention_days() == 3

    assert [record for record in caplog.records if record.levelno >= logging.WARNING] == []


def test_the_clamp_names_both_configured_values_when_it_bites(caplog):
    """An operator who set 14 has to learn that 14 is not what they got.

    Silently returning a different number than the one configured is the same class of defect
    as the code fallbacks removed in #247: the configuration says one thing and the runtime
    does another, with nothing in the log to connect them.
    """
    with caplog.at_level(logging.WARNING):
        _settings(cohort_days=14, archive_days=3).cohort_retention_days()

    warnings = [record.getMessage() for record in caplog.records if record.levelno == logging.WARNING]
    assert warnings, "clamping a configured value without saying so is a silent behaviour change"
    assert "14" in warnings[0]
    assert "3" in warnings[0]


def test_the_shipped_defaults_satisfy_the_invariant():
    """A regression guard on the shipped pair rather than on either number alone.

    This is deliberately not `== 3`: the window is an operational decision that may move
    again, while "a cohort outlives its members" is always wrong.
    """
    settings = Settings()

    assert settings.cohort_retention_days() <= settings.MAX_RESULT_AGE_DAYS


@pytest.mark.parametrize(("cohort_days", "archive_days"), [(30, 1), (365, 7), (8, 7)])
def test_the_invariant_holds_for_values_this_repository_never_sees(cohort_days, archive_days):
    """The deployed numbers do not come from this repository.

    ``vntyper-online-backend`` loads ``env_file: ${ENV_FILE:-.env.local}``, which is untracked
    and server-side, so whatever an operator sets there wins over the defaults here. The clamp
    therefore has to hold for arbitrary pairs, not just the shipped one.

    (Env vars are read when the class body is evaluated at import, so setting them from a test
    would not reach an already-imported ``Settings``. These override the instance instead,
    which is what the method actually reads.)
    """
    assert _settings(cohort_days, archive_days).cohort_retention_days() == archive_days
