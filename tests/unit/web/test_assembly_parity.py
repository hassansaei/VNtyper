"""The CLI and the server must accept the same assemblies.

`vntyper online --reference-assembly hg38_ensembl` subsets and uploads the BAM, then
gets a 422 because the server enum has only four values - reported to the user as a
generic submission failure, after the upload.
"""

import pytest

from vntyper.scripts.reference_registry import list_assemblies

pytestmark = pytest.mark.unit


def test_the_server_enum_matches_the_registry():
    # Not `from docker.app.main import ...`: `docker` resolves to the installed
    # docker-py package (a regular package, always wins the import search over the
    # local docker/ namespace package) rather than this repo's docker/ directory.
    # `tests/unit/web/conftest.py` puts docker/ itself on sys.path for exactly this
    # reason, so every other test in this directory imports `app.main`, not
    # `docker.app.main`; follow the same pattern here.
    from app.main import ReferenceAssembly

    assert {member.value for member in ReferenceAssembly} == set(list_assemblies())
