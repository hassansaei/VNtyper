import pytest

from scripts.release_policy import parse_release_tag, required_aliases

pytestmark = pytest.mark.unit


def test_release_tag_produces_all_five_aliases_in_promotion_order() -> None:
    version = parse_release_tag("v2.10.3")
    assert (version.plain, version.major, version.minor, version.patch) == ("2.10.3", 2, 10, 3)
    assert required_aliases(version) == ("v2.10.3", "2.10.3", "2.10", "2", "latest")


@pytest.mark.parametrize(
    "tag",
    ("2.0.10", "v2.0", "v2.0.10-rc1", "v2.0.10+local", "v02.0.10", "v2.00.10", "v2.0.010", "v2.0.10;echo pwned"),
)
def test_release_tag_rejects_every_non_strict_form(tag: str) -> None:
    with pytest.raises(ValueError, match="strict vMAJOR.MINOR.PATCH"):
        parse_release_tag(tag)
