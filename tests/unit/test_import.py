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
        lines = output.splitlines()
        assert re.fullmatch(r"\d{2}:\d{2}:\d{2}\.\d{2} REGCOIL: starting package logging", lines[0])
        assert re.fullmatch(r"\d{2}:\d{2}:\d{2}\.\d{2} OpenMP max threads: \d+", lines[1])
        assert re.fullmatch(r"\d{2}:\d{2}:\d{2}\.\d{2} hello", lines[2])
    finally:
        package_logger = logging.getLogger("regcoil")
        package_logger.handlers.clear()
        package_logger.propagate = True


def test_omp_max_threads():
    assert regcoil.omp_max_threads() >= 1


def test_version():
    assert regcoil.__version__
