"""What ``vntyper pipeline`` pays to import before it does any work (#262).

``cli_handlers`` imported ``cohort_summary`` (matplotlib, plotly, seaborn),
``online_mode`` (requests) and ``install_references`` at module scope. A ``pipeline``
run pays for all three and uses none of them; measured startup was 317-348 ms against
205-239 ms without them.

The deferral cannot be a plain function-local import. The suite patches
``cli_handlers.run_online_mode`` and ``cli_handlers.install_references_main`` as
**module attributes**, and moving the import inside the handler removes the attribute
those patches replace -- so they would either raise ``AttributeError`` or be silently
bypassed while the real function ran. A module ``__getattr__`` (PEP 562) does not work
either: it is consulted only for *external* attribute access, while the handlers reach
these names through in-module ``LOAD_GLOBAL``, which would raise ``NameError`` in
production. Thin module-level wrappers keep a real attribute for the patches, satisfy
the global lookup, and defer the import to call time.

Every subprocess here runs a fresh interpreter, because ``sys.modules`` is process-wide
and any earlier test that imported matplotlib would make an in-process assertion pass
for the wrong reason.
"""

from __future__ import annotations

import subprocess
import sys

import pytest

from vntyper.scripts import cli_handlers

pytestmark = pytest.mark.unit

#: Imported by the cohort, online and install-references subcommands, by nobody else.
HEAVY_MODULES = ("matplotlib", "plotly", "seaborn", "requests")


def _in_fresh_interpreter(script: str) -> str:
    """Run one snippet in a new interpreter and return its stdout.

    Args:
        script: The Python source to run with ``-c``.

    Returns:
        str: The stripped standard output.
    """
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, check=True)
    return result.stdout.strip()


def test_importing_the_cli_does_not_pull_the_plotting_or_http_stack() -> None:
    """``vntyper pipeline`` imports the CLI and then never touches any of these."""
    eager = _in_fresh_interpreter(
        f"import sys, vntyper.cli; print(','.join(name for name in {HEAVY_MODULES!r} if name in sys.modules))"
    )

    assert eager == "", f"eagerly imported: {eager}"


def test_importing_cli_handlers_does_not_pull_them_either() -> None:
    """The handlers module is where the imports were; importing it must stay cheap."""
    eager = _in_fresh_interpreter(
        "import sys, vntyper.scripts.cli_handlers; "
        f"print(','.join(name for name in {HEAVY_MODULES!r} if name in sys.modules))"
    )

    assert eager == "", f"eagerly imported: {eager}"


@pytest.mark.parametrize(
    ("wrapper", "module"),
    [
        ("aggregate_cohort", "matplotlib"),
        ("run_online_mode", "requests"),
        ("install_references_main", "vntyper.scripts.install_references"),
    ],
)
def test_calling_the_wrapper_is_what_performs_the_import(wrapper: str, module: str) -> None:
    """Deferral must actually defer, and the wrapper must actually import.

    The call is made with no arguments and its ``TypeError`` swallowed: the import
    happens on the first line of the wrapper, before the arguments are forwarded, so
    reaching the real function's signature check already proves the module loaded.
    Merely *reading* the attribute would prove nothing -- attribute access does not
    call anything.
    """
    loaded = _in_fresh_interpreter(
        "import sys\n"
        "import vntyper.scripts.cli_handlers as handlers\n"
        f"assert {module!r} not in sys.modules, 'imported before the wrapper ran'\n"
        "try:\n"
        f"    handlers.{wrapper}()\n"
        "except TypeError:\n"
        "    pass\n"
        f"print({module!r} in sys.modules)\n"
    )

    assert loaded == "True"


@pytest.mark.parametrize("wrapper", ["aggregate_cohort", "run_online_mode", "install_references_main"])
def test_the_module_attribute_patch_seam_still_works(wrapper: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """The suite patches these names as module attributes; deferral must not break it."""
    calls: list[tuple] = []
    monkeypatch.setattr(cli_handlers, wrapper, lambda *args, **kwargs: calls.append((args, kwargs)))

    getattr(cli_handlers, wrapper)("x", key="y")

    assert calls == [(("x",), {"key": "y"})]


@pytest.mark.parametrize("wrapper", ["aggregate_cohort", "run_online_mode", "install_references_main"])
def test_the_wrapper_forwards_every_argument_and_returns_the_result(
    wrapper: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A wrapper that dropped an argument or a return value would be a silent defect."""
    seen: dict[str, object] = {}

    def implementation(*args, **kwargs):
        seen["args"] = args
        seen["kwargs"] = kwargs
        return "returned"

    monkeypatch.setattr(cli_handlers, f"_load_{wrapper}", lambda: implementation)

    assert getattr(cli_handlers, wrapper)(1, 2, three=3) == "returned"
    assert seen == {"args": (1, 2), "kwargs": {"three": 3}}
