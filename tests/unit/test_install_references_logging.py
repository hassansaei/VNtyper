"""Installer logging must extend, and never replace, the CLI configuration."""

from __future__ import annotations

import logging
import shutil
import stat
from pathlib import Path

import pytest

from vntyper import cli
from vntyper.scripts import install_references, install_references_logging

pytestmark = pytest.mark.unit


@pytest.fixture()
def clean_root_logger():
    """Give each test a real, isolated root-handler lifecycle."""
    root = logging.getLogger()
    saved_handlers = root.handlers[:]
    saved_level = root.level
    for handler in root.handlers[:]:
        root.removeHandler(handler)
    yield root
    for handler in root.handlers[:]:
        root.removeHandler(handler)
        if handler not in saved_handlers:
            handler.close()
    for handler in saved_handlers:
        root.addHandler(handler)
    root.setLevel(saved_level)


def _config() -> dict[str, object]:
    """Return the complete shape `main` reads before installer dispatch."""
    return {
        "ucsc_references": {"hg19": {}},
        "ncbi_references": {},
        "ensembl_references": {},
        "vntyper_references": {},
        "aligners": {},
    }


def _prepare_main(monkeypatch: pytest.MonkeyPatch, install) -> None:
    """Replace only the external install body; logging and finalisation stay real."""
    monkeypatch.setattr(install_references, "load_install_config", lambda _path: _config())
    monkeypatch.setattr(install_references, "_install_references", install)


def test_setup_logging_preserves_cli_handlers_and_debug_level(tmp_path: Path, clean_root_logger) -> None:
    """Removing the CLI handler or forcing INFO would discard `-f` or `-l DEBUG`."""
    root = clean_root_logger
    cli_handler = logging.StreamHandler()
    root.addHandler(cli_handler)
    root.setLevel(logging.DEBUG)

    install_handler = install_references.setup_logging(tmp_path)
    try:
        logging.getLogger("vntyper.scripts.install_references").debug("digest probe detail")
        install_handler.flush()

        assert cli_handler in root.handlers
        assert root.level == logging.DEBUG
        assert install_handler in root.handlers
        assert install_handler.level == logging.NOTSET
        assert "digest probe detail" in (tmp_path / "install_references.log").read_text(encoding="utf-8")
    finally:
        root.removeHandler(install_handler)
        install_handler.close()


def test_setup_logging_provides_console_info_for_standalone_use(tmp_path: Path, clean_root_logger) -> None:
    """The Docker build invokes the installer module without CLI configuration."""
    root = clean_root_logger
    # Pytest installs its capture handler after fixture setup; direct module execution
    # genuinely has none, so reproduce that boundary immediately before the call.
    for handler in root.handlers[:]:
        root.removeHandler(handler)

    install_handler = install_references.setup_logging(tmp_path)

    assert root.level == logging.INFO
    assert install_handler in root.handlers
    assert any(type(handler) is logging.StreamHandler for handler in root.handlers)


def test_failed_standalone_file_open_leaves_root_logging_untouched(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, clean_root_logger
) -> None:
    """A file-open failure must not leave half a standalone configuration behind."""
    root = clean_root_logger
    for handler in root.handlers[:]:
        root.removeHandler(handler)
    root.setLevel(logging.ERROR)

    def refuse_file(*_args, **_kwargs):
        raise OSError("read-only log destination")

    monkeypatch.setattr(install_references_logging, "InstallLogHandler", refuse_file)

    with pytest.raises(OSError, match="read-only log destination"):
        install_references.setup_logging(tmp_path)

    assert root.handlers == []
    assert root.level == logging.ERROR


def test_main_preserves_a_preexisting_file_handler_on_success(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, clean_root_logger
) -> None:
    """Installer cleanup must remove its file handler, not an unrelated CLI file."""
    root = clean_root_logger
    cli_log = tmp_path / "cli.log"
    cli_handler = logging.FileHandler(cli_log)
    root.addHandler(cli_handler)
    root.setLevel(logging.DEBUG)
    _prepare_main(monkeypatch, lambda *_args: logging.getLogger(__name__).debug("successful install detail"))

    install_references.main(tmp_path / "refs", references_to_process=["hg19"])

    assert cli_handler in root.handlers
    assert root.level == logging.DEBUG
    logging.getLogger(__name__).debug("after successful install")
    cli_handler.flush()
    contents = cli_log.read_text(encoding="utf-8")
    assert "successful install detail" in contents
    assert "after successful install" in contents


def test_main_preserves_a_preexisting_file_handler_and_primary_exception(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, clean_root_logger
) -> None:
    """Cleanup may not replace the install failure or discard the CLI handler."""
    root = clean_root_logger
    cli_log = tmp_path / "cli.log"
    cli_handler = logging.FileHandler(cli_log)
    root.addHandler(cli_handler)
    root.setLevel(logging.DEBUG)

    def fail(*_args) -> None:
        logging.getLogger(__name__).debug("detail before failure")
        raise ValueError("primary install failure")

    _prepare_main(monkeypatch, fail)

    with pytest.raises(ValueError, match="primary install failure"):
        install_references.main(tmp_path / "refs", references_to_process=["hg19"])

    assert cli_handler in root.handlers
    assert root.level == logging.DEBUG
    logging.getLogger(__name__).debug("after failed install")
    cli_handler.flush()
    contents = cli_log.read_text(encoding="utf-8")
    assert "detail before failure" in contents
    assert "after failed install" in contents


def test_cli_same_target_log_survives_activation_and_keeps_writing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, clean_root_logger
) -> None:
    """`-f refs/install_references.log` must survive replacement of `refs`."""
    output_dir = tmp_path / "refs"
    final_log = output_dir / "install_references.log"
    captured: dict[str, logging.FileHandler] = {}

    monkeypatch.setattr(cli, "load_config", lambda _path=None: {"cli_defaults": {}})
    monkeypatch.setattr(install_references, "load_install_config", lambda _path: _config())
    real_cli_setup = cli.setup_logging

    def capture_cli_handler(**kwargs) -> None:
        real_cli_setup(**kwargs)
        matches = [
            handler
            for handler in logging.getLogger().handlers
            if isinstance(handler, logging.FileHandler) and Path(handler.baseFilename) == final_log
        ]
        assert len(matches) == 1
        captured["handler"] = matches[0]
        final_log.chmod(0o600)
        assert stat.S_IMODE(final_log.stat().st_mode) == 0o600

    monkeypatch.setattr(cli, "setup_logging", capture_cli_handler)

    def activate(*_args) -> None:
        previous = tmp_path / ".refs.previous"
        output_dir.rename(previous)
        output_dir.mkdir()
        logging.getLogger(__name__).debug("detail after directory activation")
        shutil.rmtree(previous)

    monkeypatch.setattr(install_references, "_install_references", activate)

    with pytest.raises(SystemExit) as exc_info:
        cli.main(
            [
                "-l",
                "DEBUG",
                "-f",
                str(final_log),
                "install-references",
                "-d",
                str(output_dir),
                "--references",
                "hg19",
            ]
        )

    assert exc_info.value.code == 0
    cli_file_handlers = [
        handler
        for handler in clean_root_logger.handlers
        if isinstance(handler, logging.FileHandler) and Path(handler.baseFilename) == final_log
    ]
    assert len(cli_file_handlers) == 1, "the CLI-owned handler must survive by identity"
    cli_handler = cli_file_handlers[0]
    assert cli_handler is captured["handler"]
    logging.getLogger(__name__).debug("detail after CLI dispatch")
    cli_handler.flush()

    contents = final_log.read_text(encoding="utf-8")
    assert stat.S_IMODE(final_log.stat().st_mode) == 0o600, "installer finalisation changed the operator's mode"
    assert "Logging has been set up with level 10" in contents, "early CLI records were overwritten"
    assert contents.count("detail after directory activation") == 1
    assert "detail after CLI dispatch" in contents, "the CLI handler still targets the removed inode"
    assert not (tmp_path / ".refs.install-references.log").exists()


def test_main_detaches_and_closes_its_handler_before_finalising(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, clean_root_logger
) -> None:
    """A still-attached closed handler can reopen the transient sibling log."""
    root = clean_root_logger
    root.addHandler(logging.NullHandler())
    seen: dict[str, logging.FileHandler] = {}
    real_setup = install_references.setup_logging
    real_finalize = install_references._finalize_install_log

    def capture_setup(*args, **kwargs):
        handler = real_setup(*args, **kwargs)
        seen["handler"] = handler
        return handler

    def inspect_finalize(log_file: Path, output_dir: Path) -> None:
        handler = seen["handler"]
        assert handler not in root.handlers, "installer handler was still attached during rename"
        assert handler.stream is None, "installer handler was not closed before rename"
        real_finalize(log_file, output_dir)

    monkeypatch.setattr(install_references, "setup_logging", capture_setup)
    monkeypatch.setattr(install_references, "_finalize_install_log", inspect_finalize)
    _prepare_main(monkeypatch, lambda *_args: None)

    install_references.main(tmp_path / "refs", references_to_process=["hg19"])


def test_cleanup_failure_does_not_mask_primary_install_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, clean_root_logger
) -> None:
    """A finalisation error during unwinding must leave the install error authoritative."""
    clean_root_logger.addHandler(logging.NullHandler())

    def fail(*_args) -> None:
        raise ValueError("primary install failure")

    _prepare_main(monkeypatch, fail)
    monkeypatch.setattr(
        install_references,
        "_finalize_install_log",
        lambda *_args: (_ for _ in ()).throw(RuntimeError("cleanup failure")),
    )

    with pytest.raises(ValueError, match="primary install failure"):
        install_references.main(tmp_path / "refs", references_to_process=["hg19"])


def test_cleanup_failure_is_reported_when_install_succeeded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, clean_root_logger
) -> None:
    """Suppressing cleanup errors unconditionally would report a run with no durable log."""
    clean_root_logger.addHandler(logging.NullHandler())
    _prepare_main(monkeypatch, lambda *_args: None)
    monkeypatch.setattr(
        install_references,
        "_finalize_install_log",
        lambda *_args: (_ for _ in ()).throw(RuntimeError("cleanup failure")),
    )

    with pytest.raises(RuntimeError, match="cleanup failure"):
        install_references.main(tmp_path / "refs", references_to_process=["hg19"])
