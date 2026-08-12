"""The retention window, and the bound on how far a cohort can outlive its members.

`/download/{job_id}/` requires no credential, so `MAX_RESULT_AGE_DAYS` *is* the exposure
window for a completed job: the id is the capability. Shortening it is a mitigation for
#189, not a fix -- it reduces exposure duration, not exposure -- so what is pinned here is
the invariant, not the number.

The clamp matters because the two lifetimes are driven differently.
`tasks.delete_old_results` reclaims each archive on a fixed clock from *its own* creation,
while `extend_cohort_retention` fires on every job completion and every cohort analysis, so
an actively-used cohort rolls its TTL forward from its *most recent* job. So equalising the
two durations **bounds** the gap -- from up to 11 days to at most one archive window -- and
does not close it. A cohort whose second member finishes on day 2 is open until day 5 while
its first member became cleanup-eligible on day 3.

Closing it means deriving the cohort deadline from its *oldest* member, or filtering
`/cohort-status/` to members whose archives still exist. Both are design changes beyond a
retention bound, and neither is attempted here.

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


def test_a_configured_cohort_lifetime_is_bounded_by_the_archive_window():
    """Named for what it checks, which is arithmetic on two durations.

    An earlier version of this test was called
    ``test_a_cohort_never_outlives_the_archives_it_lists`` and asserted exactly this line.
    It did not test that: equal *durations* are not equal *expiry timestamps*, because the
    cohort TTL is reset to this duration from *now* on every job completion while each
    archive ages from its own timestamp. Naming a test for a property it does not verify is
    the defect this whole change set exists to remove, so the name was corrected rather than
    the assertion strengthened -- strengthening it would require a Redis TTL clock and
    staggered archives, which is a different test.
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


def test_the_shipped_default_window_is_three_days():
    """Pin the number, not only the relation.

    An earlier version asserted only ``effective <= maximum``, which stays green if the
    default reverts to 7 -- so it could not catch the regression it existed to prevent. The
    window is an operational decision and may move again; when it does, this line moves with
    it deliberately rather than silently.
    """
    assert Settings().MAX_RESULT_AGE_DAYS == 3


def test_the_shipped_pair_is_internally_consistent():
    """The relation, kept separately from the number so each failure names its own cause."""
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
