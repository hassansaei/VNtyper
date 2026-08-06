# docker/app/utils.py

import logging
import smtplib
from email.message import EmailMessage

import bcrypt
from starlette.requests import Request

from .config import settings

logger = logging.getLogger(__name__)

# bcrypt reads at most 72 bytes of input and ignores everything after that.
# The limit is stated here so a longer passphrase can be refused outright
# rather than silently reduced to a prefix of itself.
MAX_PASSPHRASE_BYTES = 72


def client_host(request: Request) -> str | None:
    """Return the peer address of a request, or None when there is not one.

    `request.client` is populated from the `client` key of the ASGI scope, and
    that key is optional: a server that cannot name the peer, and any caller
    that drives the application directly rather than over a socket, both leave
    it out. Reading `.host` off it unguarded turns those requests into a 500
    from an AttributeError.

    Both callers use the result only to build a usage-statistics hash and a log
    line, so an unknown peer is answered with None rather than an error.

    Args:
        request: The incoming request.

    Returns:
        str | None: The peer's address, or None if the scope did not carry one.
    """
    return request.client.host if request.client else None


def hash_passphrase(passphrase: str) -> str:
    """Hash a passphrase for storage.

    Args:
        passphrase: The passphrase to hash.

    Returns:
        str: A bcrypt hash, salt included, safe to store.

    Raises:
        ValueError: If the passphrase is longer than `MAX_PASSPHRASE_BYTES`
            once encoded.
    """
    encoded = passphrase.encode("utf-8")
    if len(encoded) > MAX_PASSPHRASE_BYTES:
        msg = f"Passphrase must be at most {MAX_PASSPHRASE_BYTES} bytes"
        logger.error(msg)
        raise ValueError(msg)
    return bcrypt.hashpw(encoded, bcrypt.gensalt()).decode("ascii")


def verify_passphrase(passphrase: str, hashed_passphrase: str) -> bool:
    """Check a passphrase against a stored hash.

    `bcrypt.checkpw` compares in constant time; no comparison is written here.

    An input that could never have produced the stored hash, and a stored value
    that is not a usable hash, both answer False rather than raising: the caller
    has to be able to tell a refusal from a fault, and an unopenable record is a
    refusal.

    Args:
        passphrase: The candidate passphrase.
        hashed_passphrase: The stored hash to check against.

    Returns:
        bool: True only if the passphrase matches the stored hash.
    """
    encoded = passphrase.encode("utf-8")
    if len(encoded) > MAX_PASSPHRASE_BYTES:
        return False
    try:
        return bcrypt.checkpw(encoded, hashed_passphrase.encode("utf-8"))
    except ValueError:
        logger.error("Stored passphrase hash is not in a usable form")
        return False


def send_email(to_email: str, subject: str, content: str):
    """
    Utility function to send an email via SMTP.
    """
    msg = EmailMessage()
    msg["Subject"] = subject
    msg["From"] = settings.EMAIL_FROM
    msg["To"] = to_email
    msg.set_content(content, subtype="html")

    # Connect to the SMTP server
    with smtplib.SMTP(settings.SMTP_HOST, settings.SMTP_PORT) as server:
        server.starttls()  # Upgrade the connection to secure TLS
        server.login(settings.SMTP_USERNAME, settings.SMTP_PASSWORD)
        server.send_message(msg)
