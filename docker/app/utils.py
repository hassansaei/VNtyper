# docker/app/utils.py

import os
import smtplib
from email.message import EmailMessage
from passlib.context import CryptContext
from typing import Optional

from .config import settings

pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")


def hash_passphrase(passphrase: str) -> str:
    return pwd_context.hash(passphrase)


def verify_passphrase(passphrase: str, hashed_passphrase: str) -> bool:
    return pwd_context.verify(passphrase, hashed_passphrase)


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
