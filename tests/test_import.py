"""Smoke test: package is importable after install / editable layout."""

from io import StringIO
import re

import regcoil


def test_log_enables_timestamped_package_logging():
    stream = StringIO()
    regcoil.log(stream=stream)

    try:
        import logging

        logging.getLogger("regcoil.regcoil").info("hello")
        output = stream.getvalue().strip()
        assert re.fullmatch(r"\d{2}:\d{2}:\d{2}\.\d{2} hello", output)
    finally:
        package_logger = logging.getLogger("regcoil")
        package_logger.handlers.clear()
        package_logger.propagate = True


def test_version():
    assert regcoil.__version__
