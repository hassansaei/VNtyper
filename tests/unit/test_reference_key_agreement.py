"""The writer's output must be resolvable by every reader, for every label.

A test that only compares expected key *tuples* would have passed a design that keys
on the assembly label - which can never resolve hg19_ncbi/hg38_ncbi, because those are
aliases of GRCh37/GRCh38 and the writer only ever emits the alias-free name.
"""

import pytest

from vntyper.scripts.install_references import canonical_reference_keys
from vntyper.scripts.reference_registry import list_assemblies
from vntyper.scripts.reference_resolution import resolve_from_mapping

pytestmark = pytest.mark.unit

KINDS = ["bwa", "advntr", "shark"]


@pytest.mark.parametrize("label", sorted(list_assemblies()))
@pytest.mark.parametrize("kind", KINDS)
def test_a_full_installation_resolves_without_a_fallback(tmp_path, install_config, label, kind):
    written = {key: str(path) for key, path in canonical_reference_keys(install_config, tmp_path).items()}
    resolved = resolve_from_mapping(kind, label, written)
    assert resolved is not None, f"{kind}/{label} resolves to nothing after a full install"
    assert resolved.value, f"{kind}/{label} resolved to an empty value"
    assert not resolved.is_fallback, (
        f"{kind}/{label} fell back to {resolved.key} after a full install - "
        "a complete installation must never degrade the requested reference source"
    )
