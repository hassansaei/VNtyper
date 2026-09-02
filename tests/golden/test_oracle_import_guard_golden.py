"""Golden contracts for the standard-library oracle import-closure guard."""

from pathlib import Path

import pytest

from tests.golden import identity_oracle
from tests.golden.oracle_import_guard import assert_independent_import_closure

pytestmark = pytest.mark.golden


def test_identity_oracle_import_closure_is_independent() -> None:
    """A production helper entering any recursive local import must fail the guard."""
    oracle_path = Path(identity_oracle.__file__)
    scanned = assert_independent_import_closure(oracle_path, Path(__file__).parents[2])
    assert oracle_path.resolve() in scanned


def test_guard_rejects_nested_production_import(tmp_path: Path) -> None:
    """The guard must recurse; scanning only the entry file would miss this import."""
    entry = tmp_path / "entry.py"
    helper = tmp_path / "helper.py"
    entry.write_text("from helper import value\n", encoding="utf-8")
    helper.write_text("from vntyper.scripts.nomenclature import reconcile\nvalue = reconcile\n", encoding="utf-8")

    with pytest.raises(AssertionError, match=r"helper\.py.*vntyper\.scripts\.nomenclature"):
        assert_independent_import_closure(entry, tmp_path)


def test_guard_rejects_direct_production_import(tmp_path: Path) -> None:
    """A direct production import in the entrypoint must fail the closure guard."""
    entry = tmp_path / "entry.py"
    entry.write_text("from vntyper.scripts.nomenclature import reconcile\n", encoding="utf-8")

    with pytest.raises(AssertionError, match=r"entry\.py.*vntyper\.scripts\.nomenclature"):
        assert_independent_import_closure(entry, tmp_path)


def test_guard_rejects_direct_third_party_import(tmp_path: Path) -> None:
    """A direct dependency outside the standard library must fail the oracle closure."""
    entry = tmp_path / "entry.py"
    entry.write_text("import pandas\n", encoding="utf-8")

    with pytest.raises(AssertionError, match=r"entry\.py.*pandas"):
        assert_independent_import_closure(entry, tmp_path)


def test_guard_rejects_nested_third_party_import(tmp_path: Path) -> None:
    """A recursively reached helper must not hide a third-party dependency."""
    entry = tmp_path / "entry.py"
    helper = tmp_path / "helper.py"
    entry.write_text("from helper import value\n", encoding="utf-8")
    helper.write_text("import pandas\nvalue = pandas.DataFrame\n", encoding="utf-8")

    with pytest.raises(AssertionError, match=r"helper\.py.*pandas"):
        assert_independent_import_closure(entry, tmp_path)


@pytest.mark.parametrize(
    "source",
    [
        "import importlib\nvalue = importlib.import_module('pandas')\n",
        "import importlib.util\nvalue = importlib.import_module('pandas')\n",
        "from importlib import import_module\nvalue = import_module('pandas')\n",
        "value = __import__('pandas')\n",
        "import builtins\nvalue = builtins.__import__('pandas')\n",
        "import importlib\nvalue = importlib.import_module(name='pandas')\n",
        "import builtins\nvalue = builtins.__import__(name='pandas')\n",
        "value = __import__(name='pandas')\n",
    ],
)
def test_guard_rejects_literal_dynamic_third_party_import(tmp_path: Path, source: str) -> None:
    """Literal dynamic imports cannot bypass the standard-library-only boundary."""
    entry = tmp_path / "entry.py"
    entry.write_text(source, encoding="utf-8")

    with pytest.raises(AssertionError, match=r"entry\.py.*pandas"):
        assert_independent_import_closure(entry, tmp_path)


def test_guard_allows_standard_library_and_recursive_local_imports(tmp_path: Path) -> None:
    """The strict boundary still permits stdlib modules and recursively scanned local helpers."""
    entry = tmp_path / "entry.py"
    helper = tmp_path / "helper.py"
    entry.write_text("import json\nfrom helper import value\nresult = json.dumps(value)\n", encoding="utf-8")
    helper.write_text("from collections import Counter\nvalue = Counter('aba')\n", encoding="utf-8")

    assert assert_independent_import_closure(entry, tmp_path) == (entry.resolve(), helper.resolve())


@pytest.mark.parametrize(
    "source",
    [
        "import importlib\nvalue = importlib.import_module('vntyper.scripts.nomenclature')\n",
        "import importlib.util\nvalue = importlib.import_module('vntyper.scripts.nomenclature')\n",
        "from importlib import import_module\nvalue = import_module('vntyper.scripts.nomenclature')\n",
        "import importlib as il\nvalue = il.import_module('vntyper.scripts.nomenclature')\n",
        "import builtins\nvalue = builtins.__import__('vntyper.scripts.nomenclature')\n",
        "value = __import__('vntyper.scripts.nomenclature')\n",
    ],
)
def test_guard_rejects_positional_dynamic_production_import(tmp_path: Path, source: str) -> None:
    """Literal positional dynamic imports must not bypass production-import rejection."""
    entry = tmp_path / "entry.py"
    entry.write_text(source, encoding="utf-8")

    with pytest.raises(AssertionError, match=r"entry\.py.*vntyper\.scripts\.nomenclature"):
        assert_independent_import_closure(entry, tmp_path)


@pytest.mark.parametrize(
    "source",
    [
        "import importlib\nvalue = importlib.import_module(name='vntyper.scripts.nomenclature')\n",
        "value = __import__(name='vntyper.scripts.nomenclature')\n",
    ],
)
def test_guard_rejects_keyword_name_dynamic_production_import(tmp_path: Path, source: str) -> None:
    """Using the valid ``name=`` spelling must not bypass production-import rejection."""
    entry = tmp_path / "entry.py"
    entry.write_text(source, encoding="utf-8")

    with pytest.raises(AssertionError, match=r"entry\.py.*vntyper\.scripts\.nomenclature"):
        assert_independent_import_closure(entry, tmp_path)
