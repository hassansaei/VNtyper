"""Execution-time explicit-root resolution for selected golden tests."""

from __future__ import annotations

import os
from pathlib import Path

from tests.golden.calibration_oracle import require_explicit_roots


def resolve_explicit_roots() -> tuple[Path, Path]:
    """Resolve the two evidence roots only when a selected golden test executes."""
    return require_explicit_roots(os.environ)
