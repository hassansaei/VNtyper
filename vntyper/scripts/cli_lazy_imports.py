"""Import the subcommand-only stacks on first call, not at CLI import time.

``cohort_summary`` pulls matplotlib, plotly and seaborn; ``online_mode`` pulls requests;
``install_references`` pulls its own downloader stack. ``vntyper pipeline`` paid for all
three and used none of them -- measured startup 0.35-0.39 s and 115.5-116.2 MB peak RSS
against 0.19-0.20 s and 82.2-82.7 MB without.

These are **wrappers**, not function-local imports inside the handlers, and not a module
``__getattr__``. The suite patches ``cli_handlers.run_online_mode`` and
``cli_handlers.install_references_main`` as module attributes, so the name has to stay a
real attribute for ``monkeypatch.setattr`` to replace; and the handlers reach these names
through in-module ``LOAD_GLOBAL``, which never consults a module ``__getattr__`` and would
raise ``NameError`` in production. A wrapper satisfies both, and importing these three
names into ``cli_handlers`` keeps both properties there too: an imported name is a module
attribute and a module global alike.

They live here rather than in ``cli_handlers`` because that file was already at the
~650-line guideline and this block took it over (AGENTS.md rule 3).
"""

from __future__ import annotations

from typing import Any

# The loader is its own function per name so a test can replace it without reaching into
# another module's namespace.


def _load_aggregate_cohort():
    """Import the cohort aggregator, which pulls matplotlib, plotly and seaborn.

    Returns:
        Callable: :func:`vntyper.scripts.cohort_summary.aggregate_cohort`.
    """
    from vntyper.scripts.cohort_summary import aggregate_cohort as implementation

    return implementation


def _load_run_online_mode():
    """Import the online-mode client, which pulls requests.

    Returns:
        Callable: :func:`vntyper.scripts.online_mode.run_online_mode`.
    """
    from vntyper.scripts.online_mode import run_online_mode as implementation

    return implementation


def _load_install_references_main():
    """Import the reference installer.

    Returns:
        Callable: :func:`vntyper.scripts.install_references.main`.
    """
    from vntyper.scripts.install_references import main as implementation

    return implementation


def aggregate_cohort(*args: Any, **kwargs: Any) -> Any:
    """Run the cohort aggregation, importing its plotting stack on first call.

    Args:
        *args: Forwarded unchanged.
        **kwargs: Forwarded unchanged.

    Returns:
        Any: Whatever the implementation returns.
    """
    return _load_aggregate_cohort()(*args, **kwargs)


def run_online_mode(*args: Any, **kwargs: Any) -> Any:
    """Run the online mode, importing requests on first call.

    Args:
        *args: Forwarded unchanged.
        **kwargs: Forwarded unchanged.

    Returns:
        Any: Whatever the implementation returns.
    """
    return _load_run_online_mode()(*args, **kwargs)


def install_references_main(*args: Any, **kwargs: Any) -> Any:
    """Install reference data, importing the installer on first call.

    Args:
        *args: Forwarded unchanged.
        **kwargs: Forwarded unchanged.

    Returns:
        Any: Whatever the implementation returns.
    """
    return _load_install_references_main()(*args, **kwargs)
