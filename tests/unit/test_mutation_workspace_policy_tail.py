from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import mutation_workspace

pytestmark = pytest.mark.unit


def test_confined_path_reports_a_missing_required_parent_as_a_value_error(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="workspace path does not exist"):
        mutation_workspace.confined_path(tmp_path, "missing/child.py", must_exist=True)


def test_porcelain_rejects_duplicate_lexical_spellings_before_lexists(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    changed = tmp_path / "dir/changed.py"
    changed.parent.mkdir()
    changed.write_text("changed\n", encoding="utf-8")
    probes: list[object] = []

    def record_probe(path: object) -> bool:
        probes.append(path)
        return True

    monkeypatch.setattr(mutation_workspace.os.path, "lexists", record_probe)

    with pytest.raises(ValueError, match="unsafe workspace path"):
        mutation_workspace.parse_porcelain_z(b"?? dir/changed.py\0?? dir/./changed.py\0", tmp_path)

    assert probes == []
