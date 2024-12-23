"""
Shared fixtures (like `test_config`) for all tests.
Place this `conftest.py` in your `tests/` directory.
"""

import pytest
import json
from pathlib import Path


@pytest.fixture(scope="session")
def test_config():
    """
    Single fixture for test configuration used in all tests.

    Loads 'tests/test_data_config.json' by default or
    raises an error if not found.
    """
    config_path = Path("tests/test_data_config.json")
    if not config_path.exists():
        pytest.exit(f"Config file {config_path} not found!", returncode=1)

    with config_path.open("r") as f:
        return json.load(f)
