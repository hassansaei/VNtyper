"""
Colored logging utilities for Docker tests.

Provides consistent, colored output matching the shell script style.
"""

from enum import Enum


class LogLevel(Enum):
    """Log level types with associated colors."""

    INFO = ("\033[0;34m", "[INFO]")  # Blue
    SUCCESS = ("\033[0;32m", "[SUCCESS]")  # Green
    ERROR = ("\033[0;31m", "[ERROR]")  # Red
    WARNING = ("\033[1;33m", "[WARNING]")  # Yellow


class Logger:
    """Simple colored logger for test output."""

    RESET = "\033[0m"

    @staticmethod
    def _log(level: LogLevel, message: str) -> None:
        """Internal logging method."""
        color, prefix = level.value
        print(f"{color}{prefix}{Logger.RESET} {message}")

    @staticmethod
    def info(message: str) -> None:
        """Log informational message in blue."""
        Logger._log(LogLevel.INFO, message)

    @staticmethod
    def success(message: str) -> None:
        """Log success message in green."""
        Logger._log(LogLevel.SUCCESS, message)

    @staticmethod
    def error(message: str) -> None:
        """Log error message in red."""
        Logger._log(LogLevel.ERROR, message)

    @staticmethod
    def warning(message: str) -> None:
        """Log warning message in yellow."""
        Logger._log(LogLevel.WARNING, message)

    @staticmethod
    def separator(char: str = "=", length: int = 40) -> None:
        """Print a separator line."""
        print(char * length)

    @staticmethod
    def section(title: str, char: str = "=", length: int = 40) -> None:
        """Print a section header."""
        Logger.separator(char, length)
        Logger.info(title)
        Logger.separator(char, length)


# Convenience instance for module-level use
log = Logger()
