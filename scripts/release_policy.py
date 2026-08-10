"""Pure release-version policy for VNtyper release automation."""

import re
from dataclasses import dataclass


@dataclass(frozen=True)
class ReleaseVersion:
    """A strict release tag and its numeric components."""

    tag: str
    plain: str
    major: int
    minor: int
    patch: int


def parse_release_tag(tag: str) -> ReleaseVersion:
    """Parse a strict VNtyper release tag.

    Args:
        tag: Candidate tag to parse.

    Returns:
        The parsed release version.

    Raises:
        ValueError: If the tag is not strict ``vMAJOR.MINOR.PATCH`` syntax.
    """
    match = re.fullmatch(r"v(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)", tag)
    if match is None:
        msg = f"Release tag {tag!r} must be strict vMAJOR.MINOR.PATCH with no prerelease or build suffix."
        raise ValueError(msg)

    major, minor, patch = (int(component) for component in match.groups())
    return ReleaseVersion(tag=tag, plain=tag[1:], major=major, minor=minor, patch=patch)


def required_aliases(version: ReleaseVersion) -> tuple[str, ...]:
    """Return the container aliases required for a release version.

    Args:
        version: Parsed release version.

    Returns:
        Aliases in promotion order, from immutable release tag to ``latest``.
    """
    return (
        version.tag,
        version.plain,
        f"{version.major}.{version.minor}",
        str(version.major),
        "latest",
    )
