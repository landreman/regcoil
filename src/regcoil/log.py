from __future__ import annotations

from datetime import datetime
import logging
from typing import TextIO

_LOG_HANDLER_NAME = "regcoil-info-handler"


class _RegcoilLogFormatter(logging.Formatter):
    def formatTime(self, record, datefmt=None):
        del datefmt
        return datetime.fromtimestamp(record.created).strftime("%H:%M:%S.%f")[:-4]


def log(level=logging.INFO, stream: TextIO | None = None):
    """Enable REGCOIL package logging with timestamps."""
    package_logger = logging.getLogger("regcoil")
    package_logger.setLevel(level)
    package_logger.propagate = False

    handler: logging.StreamHandler[TextIO] | None = None
    for existing in package_logger.handlers:
        if existing.get_name() == _LOG_HANDLER_NAME and isinstance(existing, logging.StreamHandler):
            handler = existing
            break

    if handler is None:
        handler = logging.StreamHandler(stream)
        handler.set_name(_LOG_HANDLER_NAME)
        package_logger.addHandler(handler)
    elif stream is not None:
        handler.setStream(stream)

    handler.setLevel(level)
    handler.setFormatter(_RegcoilLogFormatter("%(asctime)s %(message)s"))
